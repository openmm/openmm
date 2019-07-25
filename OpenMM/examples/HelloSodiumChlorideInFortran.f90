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
! This example is written entirely in Fortran 95, using the OpenMM Fortran 
! interface module.
! ------------------------------------------------------------------------------

INCLUDE 'OpenMMFortranModule.f90'

!-------------------------------------------------------------------------------
!                ATOM, FORCE FIELD, AND SIMULATION PARAMETERS
!-------------------------------------------------------------------------------
! We'll define this module as a simplified example of the kinds of data
! structures that may already be in an MD program that is to be converted
! to use OpenMM. Note that we're using data in Angstrom and kcal units; we'll
! show how to safely convert to and from OpenMM's internal units as we go.
MODULE MyAtomInfo
    ! Simulation parameters
    ! ---------------------
    real*8 Temperature, FrictionInPerPs, SolventDielectric, SoluteDielectric
    parameter(Temperature        = 300) !Kelvins
    parameter(FrictionInPerPs    = 91)  !collisions per picosecond
    parameter(SolventDielectric  = 80)  !typical for water
    parameter(SoluteDielectric   = 2)   !typical for protein
    
    real*8 StepSizeInFs, ReportIntervalInFs, SimulationTimeInPs
    parameter(StepSizeInFs       = 2)   !integration step size (fs)
    parameter(ReportIntervalInFs = 50)  !how often for PDB frame (fs)
    parameter(SimulationTimeInPs = 100) !total simulation time (ps)
    
    logical, parameter :: WantEnergy = .true.

    ! Atom and force field information
    ! --------------------------------
    type Atom
        character*4       pdb
        real*8            mass, charge, vdwRadiusInAng, vdwEnergyInKcal
        real*8            gbsaRadiusInAng, gbsaScaleFactor
        real*8            initPosInAng(3)
        real*8            posInAng(3) ! leave room for runtime state info
    end type
    integer, parameter :: NumAtoms = 6
    type (Atom) :: atoms(NumAtoms) = (/ &
    ! pdb   mass charge vdwRad vdwEnergy  gbRad gbScale   initPos      runtime
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
! This makes use of four subroutines that encapsulate all the OpenMM calls:
!     myInitializeOpenMM
!     myStepWithOpenMM
!     myGetOpenMMState
!     myTerminateOpenMM
! and one minimalist PDB file writer that has nothing to do with OpenMM:
!     myWritePDBFrame
! All of these subroutines can be found later in this file. For use in a real
! MD code you would need to write your own interface routines along these lines.
! Note that the main program does NOT include the OpenMM module.
PROGRAM HelloSodiumChloride
    use MyAtomInfo

    ! Calculate the number of PDB frames we want to write out and how
    ! many steps to take on the GPU in between.
    integer NumReports, NumSilentSteps
    parameter(NumReports = (SimulationTimeInPs*1000 / ReportIntervalInFs + 0.5))
    parameter(NumSilentSteps = (ReportIntervalInFs / StepSizeInFs + 0.5))
    
    character*10 platformName; real*8 timeInPs, energyInKcal; integer frame

    ! This is an opaque handle to a container that holds the OpenMM runtime
    ! objects. You can use any type for this purpose as long as it is 
    ! big enough to hold a pointer. (If you use a pointer type you'll have
    ! to declare the subroutine interfaces before calling them.)
    integer*8 ommHandle

    ! Set up OpenMM data structures; returns platform name and handle.
    call myInitializeOpenMM(ommHandle, platformName)
 
    ! Run the simulation:
    !  (1) Write the first line of the PDB file and the initial configuration.
    !  (2) Run silently entirely within OpenMM between reporting intervals.
    !  (3) Write a PDB frame when the time comes.
    print "('REMARK  Using OpenMM platform ', A)", platformName
    call myGetOpenMMState(ommHandle, timeInPs, energyInKcal)
    call myWritePDBFrame(0, timeInPs, energyInKcal)
    
    do frame = 1, NumReports
        call myStepWithOpenMM(ommHandle, NumSilentSteps)
        call myGetOpenMMState(ommHandle, timeInPs, energyInKcal)
        call myWritePDBFrame(frame, timeInPs, energyInKcal)
    end do
    
    ! Clean up OpenMM data structures.
    call myTerminateOpenMM(ommHandle)
END PROGRAM

!-------------------------------------------------------------------------------
!                                PDB FILE WRITER
!-------------------------------------------------------------------------------
! Given state data which was written into the atoms array of the MyAtomInfo
! module, output a single frame (pdb "model") of the trajectory. This has
! nothing to do with OpenMM.
SUBROUTINE myWritePDBFrame(frameNum, timeInPs, energyInKcal)
    use MyAtomInfo; implicit none
    integer frameNum; real*8 timeInPs, energyInKcal

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
! The OpenMM Fortran interface module is included only at this point and below. 
! Normally these subroutines would be in a separate compilation module; we're 
! including them here for simplicity. We suggest that you write them in C++ 
! if possible, using extern "C" functions to make then callable from your
! Fortran main program. (See the C++ version of this example program for 
! an implementation of very similar routines.) However, these routines are 
! reimplemented entirely in Fortran 95 below in case you prefer.


! ------------------------------------------------------------------------------
!                      INITIALIZE OpenMM DATA STRUCTURES
! ------------------------------------------------------------------------------
! We take these actions here:
! (1) Load any available OpenMM plugins, e.g. Cuda and Brook.
! (2) Fill the OpenMM::System with the force field parameters we want to
!     use and the particular set of atoms to be simulated.
! (3) Create an Integrator and a Context associating the Integrator with
!     the System.
! (4) Select the OpenMM platform to be used.
! (5) Return an opaque handle to the Context object and the name of the 
!     Platform in use.
!
! Note that this routine must understand the calling MD code's molecule and
! force field data structures so will need to be customized for each MD code.

SUBROUTINE myInitializeOpenMM(ommHandle, platformName)
    use OpenMM; use MyAtomInfo; implicit none
    integer*8,    intent(out) :: ommHandle
    character*10, intent(out) :: platformName

    ! These are the objects we'll create here thare are stored in the
    ! Context for later access. Don't forget to delete them at the end.
    type (OpenMM_System)             system
    type (OpenMM_LangevinIntegrator) langevin
    type (OpenMM_Context)            context

    ! These are temporary OpenMM objects used and discarded here.
    type (OpenMM_StringArray)        pluginList
    type (OpenMM_Vec3Array)          initialPosInNm
    type (OpenMM_NonbondedForce)     nonbond
    type (OpenMM_GBSAOBCForce)       gbsa
    type (OpenMM_Platform)           platform
    character*100                    dirName
    integer*4                        n, ix
    real*8                           posInNm(3)

    ! Get the name of the default plugins directory, 
    ! and then load all the plugins found there.
    call OpenMM_Platform_getDefaultPluginsDirectory(dirName)
    call OpenMM_Platform_loadPluginsFromDirectory(dirName, pluginList)
    call OpenMM_StringArray_destroy(pluginList)    
    
    ! Create a System and Force objects for it. The System will take
    ! over ownership of the Forces; don't destroy them yourself.
    call OpenMM_System_create(system)
    call OpenMM_NonbondedForce_create(nonbond)
    call OpenMM_GBSAOBCForce_create(gbsa)
    
    ! Convert specific force types to generic OpenMM_Force so that we can
    ! add them to the OpenMM_System.
    ix = OpenMM_System_addForce(system, transfer(nonbond, OpenMM_Force(0)))
    ix = OpenMM_System_addForce(system, transfer(gbsa,    OpenMM_Force(0)))

    ! Specify dielectrics for GBSA implicit solvation.
    call OpenMM_GBSAOBCForce_setSolventDielectric(gbsa, SolventDielectric)
    call OpenMM_GBSAOBCForce_setSoluteDielectric(gbsa, SoluteDielectric)
    
    ! Specify the atoms and their properties:
    !  (1) System needs to know the masses.
    !  (2) NonbondedForce needs charges,van der Waals properties (in MD units!).
    !  (3) GBSA needs charge, radius, and scale factor.
    !  (4) Collect default positions for initializing the simulation later.
    call OpenMM_Vec3Array_create(initialPosInNm, NumAtoms)
    do n=1,NumAtoms
        ix = OpenMM_System_addParticle(system, atoms(n)%mass)
    
        ix = OpenMM_NonbondedForce_addParticle(nonbond,                 &
                atoms(n)%charge,                                        &
                atoms(n)%vdwRadiusInAng  * OpenMM_NmPerAngstrom         &
                                         * OpenMM_SigmaPerVdwRadius,    &
                atoms(n)%vdwEnergyInKcal * OpenMM_KJPerKcal)
    
        ix = OpenMM_GBSAOBCForce_addParticle(gbsa,                      &
                atoms(n)%charge,                                        &
                atoms(n)%gbsaRadiusInAng * OpenMM_NmPerAngstrom,        &
                atoms(n)%gbsaScaleFactor)
    
        ! Sets initPos(n) = atoms(n)%initPos * nm/Angstrom.
        call OpenMM_Vec3_scale(atoms(n)%initPosInAng,                   &
                               OpenMM_NmPerAngstrom, posInNm)
        call OpenMM_Vec3Array_set(initialPosInNm, n, posInNm)
    end do
    
    ! Choose an Integrator for advancing time, and a Context connecting the
    ! System with the Integrator for simulation. Let the Context choose the
    ! best available Platform. Initialize the configuration from the default
    ! positions we collected above. Initial velocities will be zero but could
    ! have been set here.
    call OpenMM_LangevinIntegrator_create(langevin,                     &
                                          Temperature, FrictionInPerPs, &
                                          StepSizeInFs * OpenMM_PsPerFs)

    ! Convert LangevinIntegrator to generic Integrator type for this call.
    call OpenMM_Context_create(context, system,                         &
                               transfer(langevin, OpenMM_Integrator(0)))
    call OpenMM_Context_setPositions(context, initialPosInNm)

    ! Get the platform name to return.
    call OpenMM_Context_getPlatform(context, platform)
    call OpenMM_Platform_getName(platform, platformName)
    
    ! References to the System and Integrator are in the Context, so
    ! we can extract them later for stepping and cleanup. Return an opaque 
    ! reference to the Context for use by the main program.
    ommHandle = transfer(context, ommHandle)
END SUBROUTINE


!-------------------------------------------------------------------------------
!                     COPY STATE BACK TO CPU FROM OpenMM
!-------------------------------------------------------------------------------
SUBROUTINE myGetOpenMMState(ommHandle, timeInPs, energyInKcal)
    use OpenMM; use MyAtomInfo; implicit none
    integer*8, intent(in) :: ommHandle
    real*8, intent(out)   :: timeInPs, energyInKcal
    
    type (OpenMM_State)     state
    type (OpenMM_Vec3Array) posArrayInNm
    integer                 infoMask, n
    real*8                  posInNm(3)

    type (OpenMM_Context)         context

    context = transfer(ommHandle, context)
    
    infoMask = OpenMM_State_Positions
    if (WantEnergy) then
        infoMask = infoMask + OpenMM_State_Velocities ! for KE (cheap)
        infoMask = infoMask + OpenMM_State_Energy     ! for PE (very expensive)
    end if
    ! Forces are also available (and cheap).
    
    ! Don't forget to destroy this State when you're done with it.
    call OpenMM_Context_getState(context, infoMask, OpenMM_False, state) 
    timeInPs = OpenMM_State_getTime(state) ! OpenMM time is in ps already.
    
    ! Positions are maintained as a Vec3Array inside the State. This will give
    ! us access, but don't destroy it yourself -- it will go away with the State.
    call OpenMM_State_getPositions(state, posArrayInNm)
    do n = 1, NumAtoms
        ! Sets atoms(n)%pos = posArray(n) * Angstroms/nm.
        call OpenMM_Vec3Array_get(posArrayInNm, n, posInNm)
        call OpenMM_Vec3_scale(posInNm, OpenMM_AngstromsPerNm, &
                               atoms(n)%posInAng)
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


!-------------------------------------------------------------------------------
!                     TAKE MULTIPLE STEPS USING OpenMM
!-------------------------------------------------------------------------------
SUBROUTINE myStepWithOpenMM(ommHandle, numSteps)
    use OpenMM; implicit none
    integer*8, intent(in) :: ommHandle
    integer,   intent(in) :: numSteps

    type (OpenMM_Context)    context
    type (OpenMM_Integrator) integrator
    
    context = transfer(ommHandle, context)
    call OpenMM_Context_getIntegrator(context, integrator)
    
    call OpenMM_Integrator_step(integrator, numSteps)
END SUBROUTINE


!-------------------------------------------------------------------------------
!                      DEALLOCATE ALL OpenMM OBJECTS
!-------------------------------------------------------------------------------
SUBROUTINE myTerminateOpenMM(ommHandle)
    use OpenMM; implicit none
    integer*8, intent(inout) :: ommHandle

    type (OpenMM_Context)    context
    type (OpenMM_Integrator) integrator
    type (OpenMM_System)     system
    
    context = transfer(ommHandle, context)
    call OpenMM_Context_getIntegrator(context, integrator)
    call OpenMM_Context_getSystem(context, system)
    
    call OpenMM_Context_destroy(context)
    call OpenMM_Integrator_destroy(integrator)
    call OpenMM_System_destroy(system)
END SUBROUTINE
