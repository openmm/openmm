! -----------------------------------------------------------------------------
!             OpenMM(tm) PROTOTYPE Fortran 95 Interface (June 2009)
! -----------------------------------------------------------------------------
! This is a Fortran 95 interface module providing access to the OpenMM API
! which is written in C++. At link time this module requires that the OpenMM
! C wrapper library (or object file) is available since that provides a
! simplified Fortran-style set of access methods that can be described
! adequately here without using any Fortran 2003 features.
!
! This is experimental and is not part of the OpenMM release. Improvements in
! substance and style would be greatly appreciated. If you have ideas (or 
! better code) please post to the OpenMM forum on simtk.org/home/openmm or
! if you're shy you can email Michael Sherman at msherman@stanford.edu.
! 
! Below we define two modules
!    OpenMM_Types
!    OpenMM
! Only "use OpenMM" need be included in Fortran program units since 
! that modules includes the other one.
! -----------------------------------------------------------------------------

! We use defined types containing opaque pointers as a way of getting
! a modicum of type safety without having to expose any of the OpenMM
! data structures here. You never have to do anything with those 
! pointers to deal with these objects; they get created by the API
! for you and you just pass them back to the API when you want to
! do something with them.
module OpenMM_Types
    implicit none

    ! The System, Integrator, and Context must persist between calls.
    ! They can be conveniently grouped in a RuntimeObjects structure.
    type OpenMM_System
        character, pointer :: handle => NULL()
    end type
    ! This is the generic Integrator type; it represents one of 
    ! the concrete integrators like Verlet or Langevin.
    type OpenMM_Integrator
        character, pointer :: handle => NULL()
    end type
    type OpenMM_Context
        character, pointer :: handle => NULL()
    end type

    ! This data structure can be used to hold the set of OpenMM objects
    ! that must persist from call to call while running a simulation.
    ! It contains an OpenMM_System, _Integrator, and _Context.
    type OpenMM_RuntimeObjects
        character, pointer :: handle => NULL()
    end type

    type OpenMM_State
        character, pointer :: handle => NULL()
    end type
    type OpenMM_Vec3Array
        character, pointer :: handle => NULL()
    end type
    type OpenMM_String
        character, pointer :: handle => NULL()
    end type
    ! This is the generic Force type.
    type OpenMM_Force
        character, pointer :: handle => NULL()
    end type
    type OpenMM_NonbondedForce
        character, pointer :: handle => NULL()
    end type
    type OpenMM_GBSAOBCForce
        character, pointer :: handle => NULL()
    end type
    type OpenMM_VerletIntegrator
        character, pointer :: handle => NULL()
    end type
    type OpenMM_LangevinIntegrator
        character, pointer :: handle => NULL()
    end type

    ! OpenMM::State enumerations
    integer*4 OpenMM_State_Positions, OpenMM_State_Velocities
    integer*4 OpenMM_State_Forces, OpenMM_State_Energy
    integer*4 OpenMM_State_Parameters
    parameter(OpenMM_State_Positions=1, OpenMM_State_Velocities=2)
    parameter(OpenMM_State_Forces=4, OpenMM_State_Energy=8)
    parameter(OpenMM_State_Parameters=16)

    !OpenMM::NonbondedForce enumerations
    integer*4 OpenMM_NonbondedForce_NoCutoff, OpenMM_NonbondedForce_CutoffNonPeriodic
    integer*4 OpenMM_NonbondedForce_CutoffPeriodic, OpenMM_NonbondedForce_Ewald
    parameter(OpenMM_NonbondedForce_NoCutoff=0, OpenMM_NonbondedForce_CutoffNonPeriodic=1)
    parameter(OpenMM_NonbondedForce_CutoffPeriodic=2, OpenMM_NonbondedForce_Ewald=3)

    !OpenMM units conversion constants
    real*8 OpenMM_NmPerAngstrom, OpenMM_AngstromsPerNm, OpenMM_PsPerFs, OpenMM_FsPerPs
    real*8 OpenMM_KJPerKcal, OpenMM_KcalPerKJ, OpenMM_RadiansPerDegree, OpenMM_DegreesPerRadian
    real*8 OpenMM_SigmaPerVdwRadius
    parameter(OpenMM_NmPerAngstrom=0.1, OpenMM_AngstromsPerNm=10.0)
    parameter(OpenMM_PsPerFs=0.001, OpenMM_FsPerPs=1000.0)
    parameter(OpenMM_KJPerKcal=4.184, OpenMM_KcalPerKJ=1.0/4.184)
    parameter(OpenMM_RadiansPerDegree=3.1415926535897932385/180.0)
    parameter(OpenMM_DegreesPerRadian=180.0/3.1415926535897932385)
    parameter(OpenMM_SigmaPerVdwRadius=1.78179743628068)

end module OpenMM_Types

module OpenMM
    use OpenMM_Types
    interface
        ! OpenMM_Vec3Array is an interface to the std::vector<Vec3>
        ! arrays used in various contexts by OpenMM. It is not the
        ! same as a Fortran array of Vec3s would be.
        subroutine OpenMM_Vec3Array_create(array, n)
            use OpenMM_Types
            type (OpenMM_Vec3Array) array
            integer*4 n
        end
        function OpenMM_Vec3Array_size(array)
            use OpenMM_Types
            integer*4 OpenMM_Vec3Array_size
            type (OpenMM_Vec3Array) array
        end
        subroutine OpenMM_Vec3Array_resize(array, n)
            use OpenMM_Types
            type (OpenMM_Vec3Array) array
            integer*4 n
        end
        subroutine OpenMM_Vec3Array_destroy(array)
            use OpenMM_Types
            type (OpenMM_Vec3Array) array
        end
        subroutine OpenMM_Vec3Array_append(array, v3)
            use OpenMM_Types
            type (OpenMM_Vec3Array) array
            real*8 v3(3)
        end
        subroutine OpenMM_Vec3Array_get(array, i, v3)
            use OpenMM_Types
            type (OpenMM_Vec3Array) array
            integer*4 i
            real*8, intent(out) :: v3(3)
        end
        subroutine OpenMM_Vec3Array_getScaled(array, i, s, v3)
            use OpenMM_Types
            type (OpenMM_Vec3Array) array
            integer*4 i
            real*8    s
            real*8, intent(out) :: v3(3)
        end
        subroutine OpenMM_Vec3Array_set(array, i, v3)
            use OpenMM_Types
            type (OpenMM_Vec3Array) array
            integer*4 i
            real*8, intent(in) :: v3(3)
        end
        subroutine OpenMM_Vec3Array_setScaled(array, i, v3, s)
            use OpenMM_Types
            type (OpenMM_Vec3Array) array
            integer*4 i
            real*8, intent(in) :: v3(3)
            real*8    s
        end
        subroutine OpenMM_Vec3_scale(v3in, s, v3out)
            real*8, intent(in)  :: v3in(3)
            real*8  s
            real*8, intent(out) :: v3out(3)
        end

        ! OpenMM_String is an interface to std::string, with some
        ! crude ability to copy from and out to fixed-size Fortran
        ! character arrays (with blank padding).
        subroutine OpenMM_String_create(string, initVal)
            use OpenMM_Types
            type (OpenMM_String) string
            character(*) initVal
        end
        subroutine OpenMM_String_destroy(string)
            use OpenMM_Types
            type (OpenMM_String) string
        end
        function OpenMM_String_length(string)
            use OpenMM_Types
            integer*4 OpenMM_String_length
            type (OpenMM_String) string
        end
        subroutine OpenMM_String_get(string, fstring)
            use OpenMM_Types
            type (OpenMM_String) string
            character(*) fstring
        end
        subroutine OpenMM_String_set(string, fstring)
            use OpenMM_Types
            type (OpenMM_String) string
            character(*) fstring
        end
        

        ! OpenMM::Platform
        subroutine OpenMM_Platform_loadPluginsFromDirectory(dirName)
            use OpenMM_Types
            type (OpenMM_String) dirName
        end
        subroutine OpenMM_Platform_getDefaultPluginsDirectory(dirName)
            use OpenMM_Types
            type (OpenMM_String) dirName
        end

        ! OpenMM::System
        subroutine OpenMM_System_create(system)
            use OpenMM_Types
            type (OpenMM_System) system
        end
        subroutine OpenMM_System_destroy(system)
            use OpenMM_Types
            type (OpenMM_System) system
        end
        subroutine OpenMM_System_addForce(system, force)
            use OpenMM_Types
            type (OpenMM_System) system
            type (OpenMM_Force) force
        end
        subroutine OpenMM_System_addParticle(system, mass)
            use OpenMM_Types
            type (OpenMM_System) system
            real*8 mass
        end
  
        ! OpenMM::NonbondedForce
        subroutine OpenMM_NonbondedForce_create(nonbond)
            use OpenMM_Types
            type (OpenMM_NonbondedForce) nonbond
        end
        subroutine OpenMM_NonbondedForce_destroy(nonbond)
            use OpenMM_Types
            type (OpenMM_NonbondedForce) nonbond
        end
        subroutine OpenMM_NonbondedForce_asForce(nonbond, force)
            use OpenMM_Types
            type (OpenMM_NonbondedForce) nonbond
            type (OpenMM_Force)          force
        end
        subroutine OpenMM_NonbondedForce_setNonbondedMethod(nonbond, method)
            use OpenMM_Types
            type (OpenMM_NonbondedForce) nonbond
            integer*4 method
        end
        subroutine OpenMM_NonbondedForce_setCutoffDistance(nonbond, distanceInNm)
            use OpenMM_Types
            type (OpenMM_NonbondedForce) nonbond
            real*8 distanceInNm
        end
        subroutine OpenMM_NonbondedForce_setPeriodicBoxVectors(nonbond, a, b, c)
            use OpenMM_Types
            type (OpenMM_NonbondedForce) nonbond
            real*8 a(3), b(3), c(3)
        end
        subroutine OpenMM_NonbondedForce_addParticle &
                            (nonbond, charge, sigmaInNm, vdwEnergyInKJ)
            use OpenMM_Types
            type (OpenMM_NonbondedForce) nonbond
            real*8 charge, sigmaInNm, vdwEnergyInKJ
        end

        ! OpenMM::GBSAOBCForce
        subroutine OpenMM_GBSAOBCForce_create(gbsa)
            use OpenMM_Types
            type (OpenMM_GBSAOBCForce) gbsa
        end
        subroutine OpenMM_GBSAOBCForce_destroy(gbsa)
            use OpenMM_Types
            type (OpenMM_GBSAOBCForce) gbsa
        end
        subroutine OpenMM_GBSAOBCForce_asForce(gbsa, force)
            use OpenMM_Types
            type (OpenMM_GBSAOBCForce) gbsa
            type (OpenMM_Force)          force
        end
        subroutine OpenMM_GBSAOBCForce_setSolventDielectric(gbsa, d)
            use OpenMM_Types
            type (OpenMM_GBSAOBCForce) gbsa
            real*8 d
        end
        subroutine OpenMM_GBSAOBCForce_setSoluteDielectric(gbsa, d)
            use OpenMM_Types
            type (OpenMM_GBSAOBCForce) gbsa
            real*8 d
        end
        subroutine OpenMM_GBSAOBCForce_addParticle &
                            (gbsa, charge, radiusInNm, scalingFactor)
            use OpenMM_Types
            type (OpenMM_GBSAOBCForce) gbsa
            real*8 charge, radiusInNm, scalingFactor
        end

        ! OpenMM::Integrator
        subroutine OpenMM_Integrator_step(integrator, numSteps)
            use OpenMM_Types
            type (OpenMM_Integrator) integrator
            integer*4 numSteps
        end
        subroutine OpenMM_Integrator_destroy(integrator)
            use OpenMM_Types
            type (OpenMM_Integrator) integrator
        end

        ! OpenMM::VerletIntegrator
        subroutine OpenMM_VerletIntegrator_create(verlet, stepSzInPs)
            use OpenMM_Types
            type (OpenMM_VerletIntegrator) verlet
            real*8 stepSzInPs
        end
        subroutine OpenMM_VerletIntegrator_destroy(verlet)
            use OpenMM_Types
            type (OpenMM_VerletIntegrator) verlet
        end
        subroutine OpenMM_VerletIntegrator_asIntegrator(verlet, integ)
            use OpenMM_Types
            type (OpenMM_VerletIntegrator) verlet
            type (OpenMM_Integrator)       integ
        end
        subroutine OpenMM_VerletIntegrator_step(verlet, numSteps)
            use OpenMM_Types
            type (OpenMM_VerletIntegrator) verlet
            integer*4 numSteps
        end

        ! OpenMM::LangevinIntegrator
        subroutine OpenMM_LangevinIntegrator_create &
                            (langevin, temperature, frictionInPerPs, stepSzInPs)
            use OpenMM_Types
            type (OpenMM_LangevinIntegrator) langevin
            real*8 temperature, frictionInPerPs, stepSzInPs
        end
        subroutine OpenMM_LangevinIntegrator_destroy(langevin)
            use OpenMM_Types
            type (OpenMM_LangevinIntegrator) langevin
        end
        subroutine OpenMM_LangevinIntegrator_asIntegrator(langevin, integ)
            use OpenMM_Types
            type (OpenMM_LangevinIntegrator) langevin
            type (OpenMM_Integrator)         integ
        end
        subroutine OpenMM_LangevinIntegrator_step(langevin, numSteps)
            use OpenMM_Types
            type (OpenMM_LangevinIntegrator) langevin
            integer*4 numSteps
        end

        ! OpenMM::Context
        subroutine OpenMM_Context_create(context, system, integrator)
            use OpenMM_Types
            type (OpenMM_Context) context
            type (OpenMM_System) system
            type (OpenMM_Integrator) integrator
        end
        subroutine OpenMM_Context_destroy(context)
            use OpenMM_Types
            type (OpenMM_Context) context
        end
        subroutine OpenMM_Context_setPositions(context, positions)
            use OpenMM_Types
            type (OpenMM_Context) context
            type (OpenMM_Vec3Array) positions
        end
        subroutine OpenMM_Context_setVelocities(context, velocities)
            use OpenMM_Types
            type (OpenMM_Context) context
            type (OpenMM_Vec3Array) velocities
        end
        subroutine OpenMM_Context_createState(context, types, state)
            use OpenMM_Types
            type (OpenMM_Context) context
            integer*4 types
            type (OpenMM_State) state
        end
        subroutine OpenMM_Context_getPlatformName(context, platformName)
            use OpenMM_Types
            type (OpenMM_Context) context
            character(*) platformName
        end

        ! OpenMM::State
        subroutine OpenMM_State_destroy(state)
            use OpenMM_Types
            type (OpenMM_State) state
        end
        function OpenMM_State_getTime(state)
            use OpenMM_Types
            type (OpenMM_State) state
            real*8 OpenMM_State_getTime
        end
        function OpenMM_State_getPotentialEnergy(state)
            use OpenMM_Types
            type (OpenMM_State) state
            real*8 OpenMM_State_getPotentialEnergy
        end
        function OpenMM_State_getKineticEnergy(state)
            use OpenMM_Types
            type (OpenMM_State) state
            real*8 OpenMM_State_getKineticEnergy
        end
        subroutine OpenMM_State_getPositions(state, positions)
            use OpenMM_Types
            type (OpenMM_State) state
            type (OpenMM_Vec3Array) positions
        end
        subroutine OpenMM_State_getVelocities(state, velocities)
            use OpenMM_Types
            type (OpenMM_State) state
            type (OpenMM_Vec3Array) velocities
        end
        subroutine OpenMM_RuntimeObjects_create(omm)
            use OpenMM_Types
            type (OpenMM_RuntimeObjects) omm
        end
        subroutine OpenMM_RuntimeObjects_clear(omm)
            use OpenMM_Types
            type (OpenMM_RuntimeObjects) omm
        end
        subroutine OpenMM_RuntimeObjects_destroy(omm)
            use OpenMM_Types
            type (OpenMM_RuntimeObjects) omm
        end
        subroutine OpenMM_RuntimeObjects_setSystem(omm,sys)
            use OpenMM_Types
            type (OpenMM_RuntimeObjects) omm
            type (OpenMM_System) sys
        end
        subroutine OpenMM_RuntimeObjects_setIntegrator(omm,integ)
            use OpenMM_Types
            type (OpenMM_RuntimeObjects) omm
            type (OpenMM_Integrator) integ
        end
        subroutine OpenMM_RuntimeObjects_setContext(omm,context)
            use OpenMM_Types
            type (OpenMM_RuntimeObjects) omm
            type (OpenMM_Context) context
        end
        subroutine OpenMM_RuntimeObjects_getSystem(omm,sys)
            use OpenMM_Types
            type (OpenMM_RuntimeObjects) omm
            type (OpenMM_System) sys
        end
        subroutine OpenMM_RuntimeObjects_getIntegrator(omm,integ)
            use OpenMM_Types
            type (OpenMM_RuntimeObjects) omm
            type (OpenMM_Integrator) integ
        end
        subroutine OpenMM_RuntimeObjects_getContext(omm,context)
            use OpenMM_Types
            type (OpenMM_RuntimeObjects) omm
            type (OpenMM_Context) context
        end

    end interface
end module OpenMM
