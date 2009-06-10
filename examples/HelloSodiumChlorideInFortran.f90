! -----------------------------------------------------------------------------
!     OpenMM(tm) HelloSodiumChloride example in Fortran 95 (June 2009)
! ------------------------------------------------------------------------------
! This is a complete, self-contained "hello world" example demonstrating 
! GPU-accelerated constant temperature simulation of a very simple system with
! just nonbonded forces, consisting of several sodium (Na+) and chloride (Cl-) 
! ions in implicit solvent. A multi-frame PDB file is written to stdout which 
! can be read by VMD or other visualization tool to produce an animation of the 
! resulting trajectory.
!
! Pay particular attention to the handling of units in this example. Incorrect
! handling of units is a very common error; this example shows how you can
! continue to work with Amber-style units like Angstroms, kCals, and van der
! Waals radii while correctly communicating with OpenMM in nm, kJ, and sigma.
!
! This example is written entirely in Fortran 95, using a Fortran interface 
! module which is NOT official parts of the OpenMM distribution.
! ------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!                ATOM, FORCE FIELD, AND SIMULATION PARAMETERS
!-------------------------------------------------------------------------------
! We'll define this module as a simplified example of the kinds of data
! structures that may already be in an MD program that is to be converted
! to use OpenMM.
MODULE MyAtomInfo
    ! Simulation parameters:
    real*8 Temperature, FrictionInPerPs, SolventDielectric, SoluteDielectric
    parameter(Temperature        = 300) !Kelvins
    parameter(FrictionInPerPs    = 91)  !collisions per picosecond
    parameter(SolventDielectric  = 80)  !typical for water
    parameter(SoluteDielectric   = 2)   !typical for protein
    
    real*8 StepSizeInFs, ReportIntervalInFs, SimulationTimeInPs
    parameter(StepSizeInFs       = 2)   !integration step size (fs)
    parameter(ReportIntervalInFs = 50)  !how often for PDB frame (fs)
    parameter(SimulationTimeInPs = 100) !total simulation time (ps)
    
    ! Currently energy calculation is not available in the GPU kernels so
    ! asking for it requires slow Reference Platform computation at 
    ! reporting intervals.
    logical, parameter :: WantEnergy = .true.

    ! Atom and force field information:
    type Atom
        character*4       pdb
        real*8            mass, charge, vdwRadiusInAng, vdwEnergyInKcal
        real*8            gbsaRadiusInAng, gbsaScaleFactor
        real*8            initPosInAng(3)
        real*8            posInAng(3) ! leave room for runtime state info
    end type
    integer, parameter :: NumAtoms = 6
    type (Atom) :: atoms(NumAtoms) = (/ &
    ! pdb   mass charge vdwRad vdwEnergy  gbRad gbScale   initPos      pos
Atom(' NA ',22.99,  1,  1.8680, 0.00277,  1.992,  0.8, (/ 8, 0,  0/), (/0,0,0/)),&
Atom(' CL ',35.45, -1,  2.4700, 0.1000,   1.735,  0.8, (/-8, 0,  0/), (/0,0,0/)),&
Atom(' NA ',22.99,  1,  1.8680, 0.00277,  1.992,  0.8, (/ 0, 9,  0/), (/0,0,0/)),&
Atom(' CL ',35.45, -1,  2.4700, 0.1000,   1.735,  0.8, (/ 0,-9,  0/), (/0,0,0/)),&
Atom(' NA ',22.99,  1,  1.8680, 0.00277,  1.992,  0.8, (/ 0, 0,-10/), (/0,0,0/)),&
Atom(' CL ',35.45, -1,  2.4700, 0.1000,   1.735,  0.8, (/ 0, 0, 10/), (/0,0,0/)) &
    /)
END MODULE

!-------------------------------------------------------------------------------
!                                MAIN PROGRAM
!-------------------------------------------------------------------------------
PROGRAM HelloSodiumChloride
use OpenMM
use MyAtomInfo
implicit none

!-------------------------------------------------------------------------------
!                   MODELING AND SIMULATION PARAMETERS
!-------------------------------------------------------------------------------

integer NumReports, NumSilentSteps
parameter(NumReports = (SimulationTimeInPs*1000 / ReportIntervalInFs + 0.5))
parameter(NumSilentSteps = (ReportIntervalInFs / StepSizeInFs + 0.5))

type (OpenMM_Objects) omm
character*10 platformName
real*8 timeInPs, energyInKcal
integer frame

call myInitializeOpenMM(omm, platformName)

print "('REMARK  Using OpenMM platform ', A)", platformName
call myGetOpenMMState(omm, timeInPs, energyInKcal)
call myWritePDBFrame(0, timeInPs, energyInKcal)

do frame = 1, NumReports
    call myStepWithOpenMM(omm, NumSilentSteps)
    call myGetOpenMMState(omm, timeInPs, energyInKcal)
    call myWritePDBFrame(frame, timeInPs, energyInKcal)
end do

! Clean up top-level heap allocated objects that we're done with now.
call myTerminateOpenMM(omm)

END PROGRAM

!-------------------------------------------------------------------------------
!                                PDB FILE WRITER
!-------------------------------------------------------------------------------
! Given state data, output a single frame (pdb "model") of the trajectory
SUBROUTINE myWritePDBFrame(frameNum, timeInPs, energyInKcal)
    use MyAtomInfo
    implicit none
    integer frameNum
    real*8  timeInPs, energyInKcal
    integer n
    
    print "('MODEL',5X,I0)", frameNum
    print "('REMARK 250 time=', F0.3, ' picoseconds; Energy=', F0.3, ' kcal/mole')", &
        timeInPs, energyInKcal
    do n = 1,NumAtoms
        print "('ATOM  ', I5, ' ', A4, ' SLT     1    ', 3F8.3, '  1.00  0.00')",    &
            n, atoms(n)%pdb, atoms(n)%posInAng
    end do
    
    print "('ENDMDL')"
END SUBROUTINE 

!-------------------------------------------------------------------------------
!                            OpenMM-USING CODE
!-------------------------------------------------------------------------------

SUBROUTINE myInitializeOpenMM(omm, platformName)
    use OpenMM; use MyAtomInfo
    implicit none
    type (OpenMM_Objects) omm
    character*10 platformName

    ! These are temporary OpenMM objects used and discarded here.
    type(OpenMM_Vec3Array)          initialPosInNm
    type(OpenMM_NonbondedForce)     nonbond
    type(OpenMM_Force)              nonbondAsForce
    type(OpenMM_GBSAOBCForce)       gbsa
    type(OpenMM_Force)              gbsaAsForce
    type(OpenMM_LangevinIntegrator) langevin
    real*8 posInNm(3)
    integer n

    character*50 name
    type(OpenMM_String) dir
    
    call OpenMM_String_create(dir, '')
    call OpenMM_Platform_getDefaultPluginsDirectory(dir)
    call OpenMM_String_get(dir, name)
    print *,'dir="',name,'"'
    call OpenMM_Platform_loadPluginsFromDirectory(dir)
    
    call OpenMM_System_create(omm%system)
    call OpenMM_NonbondedForce_create(nonbond)
    call OpenMM_GBSAOBCForce_create(gbsa)
    
    ! Convert specific force types to generic OpenMM_Force so that we can
    ! add them to the OpenMM_System.
    call OpenMM_NonbondedForce_asForce(nonbond, nonbondAsForce)
    call OpenMM_GBSAOBCForce_asForce(gbsa, gbsaAsForce)
    call OpenMM_System_addForce(omm%system, nonbondAsForce)
    call OpenMM_System_addForce(omm%system, gbsaAsForce)

    ! Specify dielectrics for GBSA implicit solvation.
    call OpenMM_GBSAOBCForce_setSolventDielectric(gbsa, SolventDielectric)
    call OpenMM_GBSAOBCForce_setSoluteDielectric(gbsa, SoluteDielectric)
    
    call OpenMM_Vec3Array_create(initialPosInNm, 0)
    do n=1,NumAtoms
        print *,'atom ',n,atoms(n)
        call OpenMM_System_addParticle(omm%system, atoms(n)%mass)
    
        call OpenMM_NonbondedForce_addParticle(nonbond,                 &
                atoms(n)%charge,                                        &
                atoms(n)%vdwRadiusInAng  * OpenMM_NmPerAngstrom         &
                                         * OpenMM_SigmaPerVdwRadius,    &
                atoms(n)%vdwEnergyInKcal * OpenMM_KJPerKcal)
    
        call OpenMM_GBSAOBCForce_addParticle(gbsa,                      &
                atoms(n)%charge,                                        &
                atoms(n)%gbsaRadiusInAng * OpenMM_NmPerAngstrom,        &
                atoms(n)%gbsaScaleFactor)
    
        ! Convert this atom's initial position from Angstroms to nm
        call OpenMM_Vec3_scale(atoms(n)%initPosInAng, OpenMM_NmPerAngstrom, posInNm)
        call OpenMM_Vec3Array_append(initialPosInNm, posInNm)
    end do
    
    call OpenMM_LangevinIntegrator_create(langevin,                     &
                                          Temperature, FrictionInPerPs, &
                                          StepSizeInFs * OpenMM_PsPerFs)
    call OpenMM_LangevinIntegrator_asIntegrator(langevin, omm%integrator)
    call OpenMM_Context_create(omm%context, omm%system, omm%integrator)
    call OpenMM_Context_setPositions(omm%context, initialPosInNm)
    
    call OpenMM_Context_getPlatformName(omm%context, platformName)
END SUBROUTINE


!-------------------------------------------------------------------------------
!                     COPY STATE BACK TO CPU FROM OpenMM
!-------------------------------------------------------------------------------
SUBROUTINE myGetOpenMMState(omm, timeInPs, energyInKcal)
use OpenMM; use MyAtomInfo
implicit none
type (OpenMM_Objects) omm
real*8  timeInPs, energyInKcal

type (OpenMM_State)     state
type (OpenMM_Vec3Array) posArray
integer                 infoMask
integer                 n
integer npos, i, j
real*8 energy
real*8 posInAng(3), posInNm(3)

infoMask = OpenMM_State_Positions
if (WantEnergy) then
    infoMask = infoMask + OpenMM_State_Velocities ! for KE (cheap)
    infoMask = infoMask + OpenMM_State_Energy     ! for PE (very expensive)
end if
! Forces are also available (and cheap).

! Don't forget to destroy this State when you're done with it.
call OpenMM_Context_createState(omm%context, infoMask, state) 
timeInPs = OpenMM_State_getTime(state) ! OpenMM time is in ps already.


! Positions are maintained as a Vec3Array inside the State. This will give
! us access, but don't destroy it yourself -- it will go away with the State.
call OpenMM_State_getPositions(state, posArray)
do n = 1, NumAtoms
    call OpenMM_Vec3Array_get(posArray, n, posInNm)
    call OpenMM_Vec3_scale(posInNm, OpenMM_AngstromsPerNm, atoms(n)%posInAng)
end do

energyInKcal = 0
if (WantEnergy) then
    energyInKcal = (  OpenMM_State_getPotentialEnergy(state)   &
                    + OpenMM_State_getKineticEnergy(state))    &
                   * OpenMM_KcalPerKJ
end if

! Clean up the State memory
call OpenMM_State_destroy(state)
END SUBROUTINE

SUBROUTINE myStepWithOpenMM(omm, numSteps)
    use OpenMM
    implicit none
    type (OpenMM_Objects) omm
    integer numSteps
    
    call OpenMM_Integrator_step(omm%integrator, numSteps)
END SUBROUTINE



SUBROUTINE myTerminateOpenMM(omm)
    use OpenMM
    implicit none
    type (OpenMM_Objects) omm
    call OpenMM_Objects_destroy(omm)
END SUBROUTINE
