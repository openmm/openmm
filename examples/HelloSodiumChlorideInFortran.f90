! OpenMM HelloSodiumChloride example in Fortran 95
program HelloSodiumChloride
use OpenMM
implicit none

!-------------------------------------------------------------------------------
!                   MODELING AND SIMULATION PARAMETERS
!-------------------------------------------------------------------------------
real*8 StepSizeInFs, ReportIntervalInFs, SimulationTimeInPs
parameter(StepSizeInFs = 2)
parameter(ReportIntervalInFs = 10)
parameter(SimulationTimeInPs = 100)
integer NumSilentSteps
parameter(NumSilentSteps = (ReportIntervalInFs / StepSizeInFs + 0.5))

type Atom
    character*4       pdb
    real*8            mass
    real*8            charge
    real*8            vdwRadiusInAng
    real*8            vdwEnergyInKcal
    real*8            initPosInAng(3)
end type

integer NAtom
parameter(NAtom=6)
type(Atom) atoms(NAtom)
character*50 name
type(OpenMM_String) dir

atoms(1) = Atom(' NA ', 22.99,  1, 1.8680, 0.00277, (/  8, 0, 0  /))
atoms(2) = Atom(' CL ', 35.45, -1, 2.4700, 0.1000,  (/ -8, 0, 0  /))
atoms(3) = Atom(' NA ', 22.99,  1, 1.8680, 0.00277, (/  0, 9, 0  /))
atoms(4) = Atom(' CL ', 35.45, -1, 2.4700, 0.1000,  (/  0,-9, 0  /))
atoms(5) = Atom(' NA ', 22.99,  1, 1.8680, 0.00277, (/  0, 0,-10 /))
atoms(6) = Atom(' CL ', 35.45, -1, 2.4700, 0.1000,  (/  0, 0, 10 /))

call OpenMM_String_create(dir, '')
call OpenMM_Platform_getDefaultPluginsDirectory(dir)
call OpenMM_String_get(dir, name)
print *,'dir="',name,'"'
call OpenMM_Platform_loadPluginsFromDirectory(dir)

call simulateNaCl(atoms)

! Allow subroutines to inherit type definitions from program level
CONTAINS

!-------------------------------------------------------------------------------
!                               NaCl SIMULATION
!-------------------------------------------------------------------------------
subroutine simulateNaCl(atoms)
implicit none
type(Atom)                      atoms(NAtom)
type(OpenMM_System)             system
type(OpenMM_NonbondedForce)     nonbond
type(OpenMM_Force)              nonbondAsForce
type(OpenMM_Context)            context
type(OpenMM_VerletIntegrator)   verlet
type(OpenMM_Integrator)         integrator
type(OpenMM_Vec3Array)          initialPositionsInNm
character*10 platformName
real*8 posInNm(3),a(3),b(3),c(3), cutoff
! Periodic box size and cutoff in nm
parameter(a=(/5,0,0/), b=(/0,5,0/), c=(/0,0,5/), cutoff=2)
integer i

call OpenMM_System_create(system)
call OpenMM_NonbondedForce_create(nonbond)
call OpenMM_NonbondedForce_asForce(nonbond, nonbondAsForce)
call OpenMM_System_addForce(system, nonbondAsForce)

call OpenMM_NonbondedForce_setNonbondedMethod(nonbond, &
        OpenMM_NonbondedForce_CutoffPeriodic)
! cutoff distance here is 
call OpenMM_NonbondedForce_setCutoffDistance(nonbond, cutoff)
call OpenMM_NonbondedForce_setPeriodicBoxVectors(nonbond, a, b, c)

call OpenMM_Vec3Array_create(initialPositionsInNm, 0)
do i=1,NAtom
    call OpenMM_System_addParticle(system, atoms(i)%mass)
    call OpenMM_NonbondedForce_addParticle(nonbond,                 &
            atoms(i)%charge,                                        &
            atoms(i)%vdwRadiusInAng  * OpenMM_NmPerAngstrom         &
                                     * OpenMM_SigmaPerVdwRadius,    &
            atoms(i)%vdwEnergyInKcal * OpenMM_KJPerKcal)
    ! Convert this atom's initial position from Angstroms to nm
    call OpenMM_Vec3_scale(atoms(i)%initPosInAng, OpenMM_NmPerAngstrom, posInNm)
    call OpenMM_Vec3Array_append(initialPositionsInNm, posInNm)
end do

call OpenMM_VerletIntegrator_create(verlet, StepSizeInFs * OpenMM_PsPerFs)
call OpenMM_VerletIntegrator_asIntegrator(verlet, integrator)
call OpenMM_Context_create(context, system, integrator)
call OpenMM_Context_setPositions(context, initialPositionsInNm)

call OpenMM_Context_getPlatformName(context, platformName)

print "('REMARK  Using OpenMM platform ', A)", platformName
call writePDB(atoms, context)

do
    call OpenMM_Integrator_step(integrator, NumSilentSteps)
    call writePDB(atoms, context)
    if (OpenMM_Context_getTime(context) >= SimulationTimeInPs) exit
end do

! Clean up top-level heap allocated objects that we're done with now.

call OpenMM_Vec3Array_destroy(initialPositionsInNm)
call OpenMM_Context_destroy(context)
call OpenMM_Integrator_destroy(integrator)

end subroutine

!-------------------------------------------------------------------------------
!                                PDB FILE WRITER
!-------------------------------------------------------------------------------
subroutine writePDB(atoms, context)
implicit none
type(Atom) atoms(NAtom)
type(OpenMM_Context) context
integer, save :: modelFrameNumber = 0

type (OpenMM_State) state
type (OpenMM_Vec3Array) posArray
integer npos, i, j
real*8 energy
real*8 posInAng(3), posInNm(3)

! Caution: at the moment asking for energy requires use of slow Reference 
! platform calculation.

! Don't forget to destroy this State when you're done with it.
call OpenMM_Context_createState(context,                                           &
        OpenMM_State_Positions + OpenMM_State_Velocities + OpenMM_State_Energy,    &
        state)

energy = OpenMM_State_getPotentialEnergy(state)        &
         + OpenMM_State_getKineticEnergy(state)

! Positions are maintained as a Vec3Array inside the State. This will give
! us access, but don't destroy it yourself -- it will go away with the State.
call OpenMM_State_getPositions(state, posArray)
npos = OpenMM_Vec3Array_size(posArray)
modelFrameNumber = modelFrameNumber + 1
print "('MODEL',5X,I0)", modelFrameNumber
print "('REMARK 250 time=', F0.3, ' picoseconds; Energy=', F0.3, ' kilojoules/mole')", &
    OpenMM_State_getTime(state), energy
do i = 1,npos
    call OpenMM_Vec3Array_get(posArray, i, posInNm)
    call OpenMM_Vec3_scale(posInNm, OpenMM_AngstromsPerNm, posInAng)
    print "('ATOM  ', I5, ' ', A4, ' SLT     1    ', 3F8.3, '  1.00  0.00            ')", &
        i, atoms(i)%pdb, posInAng
end do

print "('ENDMDL')"

! Clean up the State memory
call OpenMM_State_destroy(state)
end subroutine

END PROGRAM

