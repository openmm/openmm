! -----------------------------------------------------------------------------
!         OpenMM(tm) HelloArgon example in Fortran 95 (June 2009)
! -----------------------------------------------------------------------------
! This program demonstrates a simple molecular simulation using the OpenMM
! API for GPU-accelerated molecular dynamics simulation. The primary goal is
! to make sure you can compile, link, and run with OpenMM and view the output.
! The example is available in C++, C, and Fortran 95.
!
! The system modeled here is a small number of argon atoms in a vacuum.
! A multi-frame PDB file is written to stdout which  can be read by VMD or 
! other visualization tool to produce an animation of the resulting trajectory.
! -----------------------------------------------------------------------------

INCLUDE 'OpenMMFortranModule.f90'

PROGRAM HelloArgon
    use OpenMM; implicit none
    type(OpenMM_System)           system
    type(OpenMM_VerletIntegrator) verlet
    type(OpenMM_Context)          context
    type(OpenMM_Platform)         platform
    type(OpenMM_NonbondedForce)   nonbond 
    type(OpenMM_Vec3Array)        initPosInNm
    type(OpenMM_State)            state
    type(OpenMM_StringArray)      pluginList
    real*8                        timeInPs
    integer*4                     a, ix, frameNum
    character*10                  platformName
    character*100                 dirName

    ! Load any shared libraries containing GPU implementations.
    call OpenMM_Platform_getDefaultPluginsDirectory(dirName)    
    call OpenMM_Platform_loadPluginsFromDirectory(dirName, pluginList)
    call OpenMM_StringArray_destroy(pluginList)    

    ! Create a system with nonbonded forces. System takes ownership
    ! of Force; don't destroy it yourself. (We're using transfer here
    ! to recast the specific NonbondedForce to a general Force.)
    call OpenMM_System_create(system)
    call OpenMM_NonbondedForce_create(nonbond)
    ix = OpenMM_System_addForce(system, transfer(nonbond, OpenMM_Force(0)))
    
    ! Create three atoms.
    call OpenMM_Vec3Array_create(initPosInNm, 3)
    do a=1,3
        ! Space the atoms out evenly by atom index.
        call OpenMM_Vec3Array_set(initPosInNm, a, (/ 0.5d0*(a-1), 0d0, 0d0 /))
        
        ix = OpenMM_System_addParticle(system, 39.95d0) !mass of Ar, grams/mole

        ! charge, L-J sigma (nm), well depth (kJ) (vdWRad(Ar)=.188 nm)
        ix = OpenMM_NonbondedForce_addParticle(nonbond, 0d0, 0.3350d0, 0.996d0)
    end do

    ! Create particular integrator, and recast to generic one.
    call OpenMM_VerletIntegrator_create(verlet, 4d-3) !step size in ps

    ! Let OpenMM Context choose best platform.
    call OpenMM_Context_create(context, system, &
              transfer(verlet, OpenMM_Integrator(0)))
    call OpenMM_Context_getPlatform(context, platform)
    call OpenMM_Platform_getName(platform, platformName)
    print "('REMARK  Using OpenMM platform ', A)", platformName

    ! Set starting positions of the atoms. Leave time and velocity zero.
    call OpenMM_Context_setPositions(context, initPosInNm)

    ! Simulate.
    frameNum = 1
    do
        ! Output current state information.
        call OpenMM_Context_getState(context, OpenMM_State_Positions, OpenMM_False, state)
        timeInPs = OpenMM_State_getTime(state)
        call writePdbFrame(frameNum, state) !output coordinates
        call OpenMM_State_destroy(state)

        if (timeInPs .ge. 10.) then
            exit
        end if

        ! Advance state many steps at a time, for efficient use of OpenMM.
        ! (use a lot more than 10 normally) 
        call OpenMM_VerletIntegrator_step(verlet, 10)
        frameNum = frameNum + 1
    end do

    ! Free heap space for all the objects created above.
    call OpenMM_Vec3Array_destroy(initPosInNm)
    call OpenMM_Context_destroy(context)
    call OpenMM_VerletIntegrator_destroy(verlet)
    call OpenMM_System_destroy(system)
END PROGRAM

! Handy homebrew PDB writer for quick-and-dirty trajectory output.
SUBROUTINE writePDBFrame(frameNum, state)
    use OpenMM; implicit none
    integer            frameNum
    type(OpenMM_State) state

    type(OpenMM_Vec3Array) allPosInNm
    real*8                 posInNm(3), posInAng(3)
    integer                n
    
    ! Reference atomic positions in the OpenMM State.
    call OpenMM_State_getPositions(state, allPosInNm)

    print "('MODEL',5X,I0)", frameNum   ! start of frame
    do n = 1,OpenMM_Vec3Array_getSize(allPosInNm)
        call OpenMM_Vec3Array_get(allPosInNm, n, posInNm)
        call OpenMM_Vec3_scale(posInNm, 10d0, posInAng)
        print "('ATOM  ', I5, '  AR   AR     1    ', 3F8.3, '  1.00  0.00')",    &
            n, posInAng
    end do
    print "('ENDMDL')"
END SUBROUTINE 
