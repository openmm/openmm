
MODULE OpenMM_Types
    implicit none

    ! Global Constants
 
    real*8 OpenMM_KJPerKcal
    real*8 OpenMM_KcalPerKJ
    real*8 OpenMM_PsPerFs
    real*8 OpenMM_AngstromsPerNm
    real*8 OpenMM_FsPerPs
    real*8 OpenMM_RadiansPerDegree
    real*8 OpenMM_NmPerAngstrom
    real*8 OpenMM_SigmaPerVdwRadius
    real*8 OpenMM_VdwRadiusPerSigma
    real*8 OpenMM_DegreesPerRadian
    parameter(OpenMM_KJPerKcal=4.184)
    parameter(OpenMM_KcalPerKJ=0.2390057361376673)
    parameter(OpenMM_PsPerFs=0.001)
    parameter(OpenMM_AngstromsPerNm=10)
    parameter(OpenMM_FsPerPs=1000)
    parameter(OpenMM_RadiansPerDegree=0.017453292519943295)
    parameter(OpenMM_NmPerAngstrom=0.1)
    parameter(OpenMM_SigmaPerVdwRadius=1.7817974362806785)
    parameter(OpenMM_VdwRadiusPerSigma=0.5612310241546865)
    parameter(OpenMM_DegreesPerRadian=57.29577951308232)

    ! Type Declarations
 
    type OpenMM_HarmonicBondForce
        integer*8 :: handle = 0
    end type
 
    type OpenMM_BrownianIntegrator
        integer*8 :: handle = 0
    end type
 
    type OpenMM_OpenMMException
        integer*8 :: handle = 0
    end type
 
    type OpenMM_NonbondedForce
        integer*8 :: handle = 0
    end type
 
    type OpenMM_VariableLangevinIntegrator
        integer*8 :: handle = 0
    end type
 
    type OpenMM_GBVIForce
        integer*8 :: handle = 0
    end type
 
    type OpenMM_Context
        integer*8 :: handle = 0
    end type
 
    type OpenMM_GBSAOBCForce
        integer*8 :: handle = 0
    end type
 
    type OpenMM_VariableVerletIntegrator
        integer*8 :: handle = 0
    end type
 
    type OpenMM_CMMotionRemover
        integer*8 :: handle = 0
    end type
 
    type OpenMM_VerletIntegrator
        integer*8 :: handle = 0
    end type
 
    type OpenMM_ContextImpl
        integer*8 :: handle = 0
    end type
 
    type OpenMM_RBTorsionForce
        integer*8 :: handle = 0
    end type
 
    type OpenMM_LangevinIntegrator
        integer*8 :: handle = 0
    end type
 
    type OpenMM_Force
        integer*8 :: handle = 0
    end type
 
    type OpenMM_HarmonicAngleForce
        integer*8 :: handle = 0
    end type
 
    type OpenMM_AndersenThermostat
        integer*8 :: handle = 0
    end type
 
    type OpenMM_ForceImpl
        integer*8 :: handle = 0
    end type
 
    type OpenMM_Platform
        integer*8 :: handle = 0
    end type
 
    type OpenMM_State
        integer*8 :: handle = 0
    end type
 
    type OpenMM_PeriodicTorsionForce
        integer*8 :: handle = 0
    end type
 
    type OpenMM_Integrator
        integer*8 :: handle = 0
    end type
 
    type OpenMM_System
        integer*8 :: handle = 0
    end type
 
    type OpenMM_Vec3Array
        integer*8 :: handle = 0
    end type

    type OpenMM_StringArray
        integer*8 :: handle = 0
    end type

    type OpenMM_BondArray
        integer*8 :: handle = 0
    end type

    type OpenMM_ParameterArray
        integer*8 :: handle = 0
    end type

    ! Enumerations

    integer*4 OpenMM_False
    integer*4 OpenMM_True
    parameter(OpenMM_False=0)
    parameter(OpenMM_True=1)
 
    integer*4 OpenMM_NonbondedForce_NoCutoff
    integer*4 OpenMM_NonbondedForce_CutoffNonPeriodic
    integer*4 OpenMM_NonbondedForce_CutoffPeriodic
    integer*4 OpenMM_NonbondedForce_Ewald
    integer*4 OpenMM_NonbondedForce_PME
    parameter(OpenMM_NonbondedForce_NoCutoff=0)
    parameter(OpenMM_NonbondedForce_CutoffNonPeriodic=1)
    parameter(OpenMM_NonbondedForce_CutoffPeriodic=2)
    parameter(OpenMM_NonbondedForce_Ewald=3)
    parameter(OpenMM_NonbondedForce_PME=4)

    integer*4 OpenMM_State_Positions
    integer*4 OpenMM_State_Velocities
    integer*4 OpenMM_State_Forces
    integer*4 OpenMM_State_Energy
    integer*4 OpenMM_State_Parameters
    parameter(OpenMM_State_Positions=1)
    parameter(OpenMM_State_Velocities=2)
    parameter(OpenMM_State_Forces=4)
    parameter(OpenMM_State_Energy=8)
    parameter(OpenMM_State_Parameters=16)

END MODULE OpenMM_Types

MODULE OpenMM
    use OpenMM_Types; implicit none
    interface

        ! OpenMM_Vec3
        subroutine OpenMM_Vec3_scale(vec, scale, result)
            use OpenMM_Types; implicit none
            real*8 vec(3)
            real*8 scale
            real*8 result(3)
        end

        ! OpenMM_Vec3Array
        subroutine OpenMM_Vec3Array_create(result, size)
            use OpenMM_Types; implicit none
            integer*4 size
            type (OpenMM_Vec3Array) result
        end
        subroutine OpenMM_Vec3Array_destroy(destroy)
            use OpenMM_Types; implicit none
            type (OpenMM_Vec3Array) destroy
        end
        function OpenMM_Vec3Array_getSize(target)
            use OpenMM_Types; implicit none
            type (OpenMM_Vec3Array) target
            integer*4 OpenMM_Vec3Array_getSize
        end
        subroutine OpenMM_Vec3Array_resize(target, size)
            use OpenMM_Types; implicit none
            type (OpenMM_Vec3Array) target
            integer*4 size
        end
        subroutine OpenMM_Vec3Array_append(target, vec)
            use OpenMM_Types; implicit none
            type (OpenMM_Vec3Array) target
            real*8 vec(3)
        end
        subroutine OpenMM_Vec3Array_set(target, index, vec)
            use OpenMM_Types; implicit none
            type (OpenMM_Vec3Array) target
            integer*4 index
            real*8 vec(3)
        end
        subroutine OpenMM_Vec3Array_get(target, index, result)
            use OpenMM_Types; implicit none
            type (OpenMM_Vec3Array) target
            integer*4 index
            real*8 result(3)
        end

        ! OpenMM_StringArray
        subroutine OpenMM_StringArray_create(result, size)
            use OpenMM_Types; implicit none
            integer*4 size
            type (OpenMM_StringArray) result
        end
        subroutine OpenMM_StringArray_destroy(destroy)
            use OpenMM_Types; implicit none
            type (OpenMM_StringArray) destroy
        end
        function OpenMM_StringArray_getSize(target)
            use OpenMM_Types; implicit none
            type (OpenMM_StringArray) target
            integer*4 OpenMM_StringArray_getSize
        end
        subroutine OpenMM_StringArray_resize(target, size)
            use OpenMM_Types; implicit none
            type (OpenMM_StringArray) target
            integer*4 size
        end
        subroutine OpenMM_StringArray_append(target, str)
            use OpenMM_Types; implicit none
            type (OpenMM_StringArray) target
            character(*) str
        end
        subroutine OpenMM_StringArray_set(target, index, str)
            use OpenMM_Types; implicit none
            type (OpenMM_StringArray) target
            integer*4 index
            character(*) str
        end
        subroutine OpenMM_StringArray_get(target, index, result)
            use OpenMM_Types; implicit none
            type (OpenMM_StringArray) target
            integer*4 index
            character(*) result
        end

        ! OpenMM_BondArray
        subroutine OpenMM_BondArray_create(result, size)
            use OpenMM_Types; implicit none
            integer*4 size
            type (OpenMM_BondArray) result
        end
        subroutine OpenMM_BondArray_destroy(destroy)
            use OpenMM_Types; implicit none
            type (OpenMM_BondArray) destroy
        end
        function OpenMM_BondArray_getSize(target)
            use OpenMM_Types; implicit none
            type (OpenMM_BondArray) target
            integer*4 OpenMM_BondArray_getSize
        end
        subroutine OpenMM_BondArray_resize(target, size)
            use OpenMM_Types; implicit none
            type (OpenMM_BondArray) target
            integer*4 size
        end
        subroutine OpenMM_BondArray_append(target, particle1, particle2)
            use OpenMM_Types; implicit none
            type (OpenMM_BondArray) target
            integer*4 particle1
            integer*4 particle2
        end
        subroutine OpenMM_BondArray_set(target, index, particle1, particle2)
            use OpenMM_Types; implicit none
            type (OpenMM_BondArray) target
            integer*4 index
            integer*4 particle1
            integer*4 particle2
        end
        subroutine OpenMM_BondArray_get(target, index, particle1, particle2)
            use OpenMM_Types; implicit none
            type (OpenMM_BondArray) target
            integer*4 index
            integer*4 particle1
            integer*4 particle2
        end

        ! OpenMM_ParameterArray
        function OpenMM_ParameterArray_getSize(target)
            use OpenMM_Types; implicit none
            type (OpenMM_ParameterArray) target
            integer*4 OpenMM_ParameterArray_getSize
        end
        subroutine OpenMM_ParameterArray_get(target, name, result)
            use OpenMM_Types; implicit none
            type (OpenMM_ParameterArray) target
            character(*) name
            character(*) result
        end

 
 

        ! OpenMM::HarmonicBondForce
        subroutine OpenMM_HarmonicBondForce_create(result)
            use OpenMM_Types; implicit none
            type (OpenMM_HarmonicBondForce) result
        end
        subroutine OpenMM_HarmonicBondForce_destroy(destroy)
            use OpenMM_Types; implicit none
            type (OpenMM_HarmonicBondForce) destroy
        end
        function OpenMM_HarmonicBondForce_getNumBonds(target)
            use OpenMM_Types; implicit none
            type (OpenMM_HarmonicBondForce) target
            integer*4  OpenMM_HarmonicBondForce_getNumBonds
        end
        function OpenMM_HarmonicBondForce_addBond(target, particle1, particle2, length, k)
            use OpenMM_Types; implicit none
            type (OpenMM_HarmonicBondForce) target
            integer*4 particle1
            integer*4 particle2
            real*8 length
            real*8 k
            integer*4  OpenMM_HarmonicBondForce_addBond
        end
        subroutine OpenMM_HarmonicBondForce_getBondParameters(target, index, particle1, particle2, length, k)
            use OpenMM_Types; implicit none
            type (OpenMM_HarmonicBondForce) target
            integer*4 index
            integer*4 particle1
            integer*4 particle2
            real*8 length
            real*8 k
        end
        subroutine OpenMM_HarmonicBondForce_setBondParameters(target, index, particle1, particle2, length, k)
            use OpenMM_Types; implicit none
            type (OpenMM_HarmonicBondForce) target
            integer*4 index
            integer*4 particle1
            integer*4 particle2
            real*8 length
            real*8 k
        end

        ! OpenMM::BrownianIntegrator
        subroutine OpenMM_BrownianIntegrator_create(result, temperature, frictionCoeff, stepSize)
            use OpenMM_Types; implicit none
            type (OpenMM_BrownianIntegrator) result
            real*8 temperature
            real*8 frictionCoeff
            real*8 stepSize
        end
        subroutine OpenMM_BrownianIntegrator_destroy(destroy)
            use OpenMM_Types; implicit none
            type (OpenMM_BrownianIntegrator) destroy
        end
        function OpenMM_BrownianIntegrator_getTemperature(target)
            use OpenMM_Types; implicit none
            type (OpenMM_BrownianIntegrator) target
            real*8  OpenMM_BrownianIntegrator_getTemperature
        end
        subroutine OpenMM_BrownianIntegrator_setTemperature(target, temp)
            use OpenMM_Types; implicit none
            type (OpenMM_BrownianIntegrator) target
            real*8 temp
        end
        function OpenMM_BrownianIntegrator_getFriction(target)
            use OpenMM_Types; implicit none
            type (OpenMM_BrownianIntegrator) target
            real*8  OpenMM_BrownianIntegrator_getFriction
        end
        subroutine OpenMM_BrownianIntegrator_setFriction(target, coeff)
            use OpenMM_Types; implicit none
            type (OpenMM_BrownianIntegrator) target
            real*8 coeff
        end
        function OpenMM_BrownianIntegrator_getRandomNumberSeed(target)
            use OpenMM_Types; implicit none
            type (OpenMM_BrownianIntegrator) target
            integer*4  OpenMM_BrownianIntegrator_getRandomNumberSeed
        end
        subroutine OpenMM_BrownianIntegrator_setRandomNumberSeed(target, seed)
            use OpenMM_Types; implicit none
            type (OpenMM_BrownianIntegrator) target
            integer*4 seed
        end
        subroutine OpenMM_BrownianIntegrator_step(target, steps)
            use OpenMM_Types; implicit none
            type (OpenMM_BrownianIntegrator) target
            integer*4 steps
        end

        ! OpenMM::OpenMMException
        subroutine OpenMM_OpenMMException_create(result, message)
            use OpenMM_Types; implicit none
            type (OpenMM_OpenMMException) result
            character(*) message
        end
        subroutine OpenMM_OpenMMException_destroy(destroy)
            use OpenMM_Types; implicit none
            type (OpenMM_OpenMMException) destroy
        end
        subroutine OpenMM_OpenMMException_what(target, result)
            use OpenMM_Types; implicit none
            type (OpenMM_OpenMMException) target
            character(*) result
        end

        ! OpenMM::NonbondedForce
        subroutine OpenMM_NonbondedForce_create(result)
            use OpenMM_Types; implicit none
            type (OpenMM_NonbondedForce) result
        end
        subroutine OpenMM_NonbondedForce_destroy(destroy)
            use OpenMM_Types; implicit none
            type (OpenMM_NonbondedForce) destroy
        end
        function OpenMM_NonbondedForce_getNumParticles(target)
            use OpenMM_Types; implicit none
            type (OpenMM_NonbondedForce) target
            integer*4  OpenMM_NonbondedForce_getNumParticles
        end
        function OpenMM_NonbondedForce_getNumExceptions(target)
            use OpenMM_Types; implicit none
            type (OpenMM_NonbondedForce) target
            integer*4  OpenMM_NonbondedForce_getNumExceptions
        end
        subroutine OpenMM_NonbondedForce_getNonbondedMethod(target, result)
            use OpenMM_Types; implicit none
            type (OpenMM_NonbondedForce) target
            integer*4 result
        end
        subroutine OpenMM_NonbondedForce_setNonbondedMethod(target, method)
            use OpenMM_Types; implicit none
            type (OpenMM_NonbondedForce) target
            integer*4 method
        end
        function OpenMM_NonbondedForce_getCutoffDistance(target)
            use OpenMM_Types; implicit none
            type (OpenMM_NonbondedForce) target
            real*8  OpenMM_NonbondedForce_getCutoffDistance
        end
        subroutine OpenMM_NonbondedForce_setCutoffDistance(target, distance)
            use OpenMM_Types; implicit none
            type (OpenMM_NonbondedForce) target
            real*8 distance
        end
        function OpenMM_NonbondedForce_getReactionFieldDielectric(target)
            use OpenMM_Types; implicit none
            type (OpenMM_NonbondedForce) target
            real*8  OpenMM_NonbondedForce_getReactionFieldDielectric
        end
        subroutine OpenMM_NonbondedForce_setReactionFieldDielectric(target, dielectric)
            use OpenMM_Types; implicit none
            type (OpenMM_NonbondedForce) target
            real*8 dielectric
        end
        function OpenMM_NonbondedForce_getEwaldErrorTolerance(target)
            use OpenMM_Types; implicit none
            type (OpenMM_NonbondedForce) target
            real*8  OpenMM_NonbondedForce_getEwaldErrorTolerance
        end
        subroutine OpenMM_NonbondedForce_setEwaldErrorTolerance(target, tol)
            use OpenMM_Types; implicit none
            type (OpenMM_NonbondedForce) target
            real*8 tol
        end
        subroutine OpenMM_NonbondedForce_getPeriodicBoxVectors(target, a, b, c)
            use OpenMM_Types; implicit none
            type (OpenMM_NonbondedForce) target
            real*8 a(3)
            real*8 b(3)
            real*8 c(3)
        end
        subroutine OpenMM_NonbondedForce_setPeriodicBoxVectors(target, a, b, c)
            use OpenMM_Types; implicit none
            type (OpenMM_NonbondedForce) target
            real*8 a(3)
            real*8 b(3)
            real*8 c(3)
        end
        function OpenMM_NonbondedForce_addParticle(target, charge, sigma, epsilon)
            use OpenMM_Types; implicit none
            type (OpenMM_NonbondedForce) target
            real*8 charge
            real*8 sigma
            real*8 epsilon
            integer*4  OpenMM_NonbondedForce_addParticle
        end
        subroutine OpenMM_NonbondedForce_getParticleParameters(target, index, charge, sigma, epsilon)
            use OpenMM_Types; implicit none
            type (OpenMM_NonbondedForce) target
            integer*4 index
            real*8 charge
            real*8 sigma
            real*8 epsilon
        end
        subroutine OpenMM_NonbondedForce_setParticleParameters(target, index, charge, sigma, epsilon)
            use OpenMM_Types; implicit none
            type (OpenMM_NonbondedForce) target
            integer*4 index
            real*8 charge
            real*8 sigma
            real*8 epsilon
        end
        function OpenMM_NonbondedForce_addException(target, particle1, particle2, chargeProd, sigma, epsilon, replace)
            use OpenMM_Types; implicit none
            type (OpenMM_NonbondedForce) target
            integer*4 particle1
            integer*4 particle2
            real*8 chargeProd
            real*8 sigma
            real*8 epsilon
            integer*4 replace
            integer*4  OpenMM_NonbondedForce_addException
        end
        subroutine OpenMM_NonbondedForce_getExceptionParameters(target, index, particle1, particle2, chargeProd, sigma, epsilon)
            use OpenMM_Types; implicit none
            type (OpenMM_NonbondedForce) target
            integer*4 index
            integer*4 particle1
            integer*4 particle2
            real*8 chargeProd
            real*8 sigma
            real*8 epsilon
        end
        subroutine OpenMM_NonbondedForce_setExceptionParameters(target, index, particle1, particle2, chargeProd, sigma, epsilon)
            use OpenMM_Types; implicit none
            type (OpenMM_NonbondedForce) target
            integer*4 index
            integer*4 particle1
            integer*4 particle2
            real*8 chargeProd
            real*8 sigma
            real*8 epsilon
        end
        subroutine OpenMM_NonbondedForce_createExceptionsFromBonds(target, bonds, coulomb14Scale, lj14Scale)
            use OpenMM_Types; implicit none
            type (OpenMM_NonbondedForce) target
            type (OpenMM_BondArray) bonds
            real*8 coulomb14Scale
            real*8 lj14Scale
        end

        ! OpenMM::VariableLangevinIntegrator
        subroutine OpenMM_VariableLangevinIntegrator_create(result, temperature, frictionCoeff, errorTol)
            use OpenMM_Types; implicit none
            type (OpenMM_VariableLangevinIntegrator) result
            real*8 temperature
            real*8 frictionCoeff
            real*8 errorTol
        end
        subroutine OpenMM_VariableLangevinIntegrator_destroy(destroy)
            use OpenMM_Types; implicit none
            type (OpenMM_VariableLangevinIntegrator) destroy
        end
        function OpenMM_VariableLangevinIntegrator_getTemperature(target)
            use OpenMM_Types; implicit none
            type (OpenMM_VariableLangevinIntegrator) target
            real*8  OpenMM_VariableLangevinIntegrator_getTemperature
        end
        subroutine OpenMM_VariableLangevinIntegrator_setTemperature(target, temp)
            use OpenMM_Types; implicit none
            type (OpenMM_VariableLangevinIntegrator) target
            real*8 temp
        end
        function OpenMM_VariableLangevinIntegrator_getFriction(target)
            use OpenMM_Types; implicit none
            type (OpenMM_VariableLangevinIntegrator) target
            real*8  OpenMM_VariableLangevinIntegrator_getFriction
        end
        subroutine OpenMM_VariableLangevinIntegrator_setFriction(target, coeff)
            use OpenMM_Types; implicit none
            type (OpenMM_VariableLangevinIntegrator) target
            real*8 coeff
        end
        function OpenMM_VariableLangevinIntegrator_getErrorTolerance(target)
            use OpenMM_Types; implicit none
            type (OpenMM_VariableLangevinIntegrator) target
            real*8  OpenMM_VariableLangevinIntegrator_getErrorTolerance
        end
        subroutine OpenMM_VariableLangevinIntegrator_setErrorTolerance(target, tol)
            use OpenMM_Types; implicit none
            type (OpenMM_VariableLangevinIntegrator) target
            real*8 tol
        end
        function OpenMM_VariableLangevinIntegrator_getRandomNumberSeed(target)
            use OpenMM_Types; implicit none
            type (OpenMM_VariableLangevinIntegrator) target
            integer*4  OpenMM_VariableLangevinIntegrator_getRandomNumberSeed
        end
        subroutine OpenMM_VariableLangevinIntegrator_setRandomNumberSeed(target, seed)
            use OpenMM_Types; implicit none
            type (OpenMM_VariableLangevinIntegrator) target
            integer*4 seed
        end
        subroutine OpenMM_VariableLangevinIntegrator_step(target, steps)
            use OpenMM_Types; implicit none
            type (OpenMM_VariableLangevinIntegrator) target
            integer*4 steps
        end
        subroutine OpenMM_VariableLangevinIntegrator_stepTo(target, time)
            use OpenMM_Types; implicit none
            type (OpenMM_VariableLangevinIntegrator) target
            real*8 time
        end

        ! OpenMM::GBVIForce
        subroutine OpenMM_GBVIForce_create(result)
            use OpenMM_Types; implicit none
            type (OpenMM_GBVIForce) result
        end
        subroutine OpenMM_GBVIForce_destroy(destroy)
            use OpenMM_Types; implicit none
            type (OpenMM_GBVIForce) destroy
        end
        function OpenMM_GBVIForce_getNumParticles(target)
            use OpenMM_Types; implicit none
            type (OpenMM_GBVIForce) target
            integer*4  OpenMM_GBVIForce_getNumParticles
        end
        function OpenMM_GBVIForce_addParticle(target, charge, radius, gamma)
            use OpenMM_Types; implicit none
            type (OpenMM_GBVIForce) target
            real*8 charge
            real*8 radius
            real*8 gamma
            integer*4  OpenMM_GBVIForce_addParticle
        end
        subroutine OpenMM_GBVIForce_getParticleParameters(target, index, charge, radius, gamma)
            use OpenMM_Types; implicit none
            type (OpenMM_GBVIForce) target
            integer*4 index
            real*8 charge
            real*8 radius
            real*8 gamma
        end
        subroutine OpenMM_GBVIForce_setParticleParameters(target, index, charge, radius, gamma)
            use OpenMM_Types; implicit none
            type (OpenMM_GBVIForce) target
            integer*4 index
            real*8 charge
            real*8 radius
            real*8 gamma
        end
        function OpenMM_GBVIForce_getSolventDielectric(target)
            use OpenMM_Types; implicit none
            type (OpenMM_GBVIForce) target
            real*8  OpenMM_GBVIForce_getSolventDielectric
        end
        subroutine OpenMM_GBVIForce_setSolventDielectric(target, dielectric)
            use OpenMM_Types; implicit none
            type (OpenMM_GBVIForce) target
            real*8 dielectric
        end
        function OpenMM_GBVIForce_getSoluteDielectric(target)
            use OpenMM_Types; implicit none
            type (OpenMM_GBVIForce) target
            real*8  OpenMM_GBVIForce_getSoluteDielectric
        end
        subroutine OpenMM_GBVIForce_setSoluteDielectric(target, dielectric)
            use OpenMM_Types; implicit none
            type (OpenMM_GBVIForce) target
            real*8 dielectric
        end

        ! OpenMM::Context
        subroutine OpenMM_Context_create(result, system, integrator)
            use OpenMM_Types; implicit none
            type (OpenMM_Context) result
            type (OpenMM_System) system
            type (OpenMM_Integrator) integrator
        end
        subroutine OpenMM_Context_create_2(result, system, integrator, platform)
            use OpenMM_Types; implicit none
            type (OpenMM_Context) result
            type (OpenMM_System) system
            type (OpenMM_Integrator) integrator
            type (OpenMM_Platform) platform
        end
        subroutine OpenMM_Context_destroy(destroy)
            use OpenMM_Types; implicit none
            type (OpenMM_Context) destroy
        end
        subroutine OpenMM_Context_getSystem(target, result)
            use OpenMM_Types; implicit none
            type (OpenMM_Context) target
            type (OpenMM_System) result
        end
        subroutine OpenMM_Context_getIntegrator(target, result)
            use OpenMM_Types; implicit none
            type (OpenMM_Context) target
            type (OpenMM_Integrator) result
        end
        subroutine OpenMM_Context_getPlatform(target, result)
            use OpenMM_Types; implicit none
            type (OpenMM_Context) target
            type (OpenMM_Platform) result
        end
        subroutine OpenMM_Context_getState(target, types, result)
            use OpenMM_Types; implicit none
            type (OpenMM_Context) target
            integer*4 types
            type (OpenMM_State) result
        end
        subroutine OpenMM_Context_setTime(target, time)
            use OpenMM_Types; implicit none
            type (OpenMM_Context) target
            real*8 time
        end
        subroutine OpenMM_Context_setPositions(target, positions)
            use OpenMM_Types; implicit none
            type (OpenMM_Context) target
            type (OpenMM_Vec3Array) positions
        end
        subroutine OpenMM_Context_setVelocities(target, velocities)
            use OpenMM_Types; implicit none
            type (OpenMM_Context) target
            type (OpenMM_Vec3Array) velocities
        end
        function OpenMM_Context_getParameter(target, name)
            use OpenMM_Types; implicit none
            type (OpenMM_Context) target
            character(*) name
            real*8  OpenMM_Context_getParameter
        end
        subroutine OpenMM_Context_setParameter(target, name, value)
            use OpenMM_Types; implicit none
            type (OpenMM_Context) target
            character(*) name
            real*8 value
        end
        subroutine OpenMM_Context_reinitialize(target)
            use OpenMM_Types; implicit none
            type (OpenMM_Context) target
        end

        ! OpenMM::GBSAOBCForce
        subroutine OpenMM_GBSAOBCForce_create(result)
            use OpenMM_Types; implicit none
            type (OpenMM_GBSAOBCForce) result
        end
        subroutine OpenMM_GBSAOBCForce_destroy(destroy)
            use OpenMM_Types; implicit none
            type (OpenMM_GBSAOBCForce) destroy
        end
        function OpenMM_GBSAOBCForce_getNumParticles(target)
            use OpenMM_Types; implicit none
            type (OpenMM_GBSAOBCForce) target
            integer*4  OpenMM_GBSAOBCForce_getNumParticles
        end
        function OpenMM_GBSAOBCForce_addParticle(target, charge, radius, scalingFactor)
            use OpenMM_Types; implicit none
            type (OpenMM_GBSAOBCForce) target
            real*8 charge
            real*8 radius
            real*8 scalingFactor
            integer*4  OpenMM_GBSAOBCForce_addParticle
        end
        subroutine OpenMM_GBSAOBCForce_getParticleParameters(target, index, charge, radius, scalingFactor)
            use OpenMM_Types; implicit none
            type (OpenMM_GBSAOBCForce) target
            integer*4 index
            real*8 charge
            real*8 radius
            real*8 scalingFactor
        end
        subroutine OpenMM_GBSAOBCForce_setParticleParameters(target, index, charge, radius, scalingFactor)
            use OpenMM_Types; implicit none
            type (OpenMM_GBSAOBCForce) target
            integer*4 index
            real*8 charge
            real*8 radius
            real*8 scalingFactor
        end
        function OpenMM_GBSAOBCForce_getSolventDielectric(target)
            use OpenMM_Types; implicit none
            type (OpenMM_GBSAOBCForce) target
            real*8  OpenMM_GBSAOBCForce_getSolventDielectric
        end
        subroutine OpenMM_GBSAOBCForce_setSolventDielectric(target, dielectric)
            use OpenMM_Types; implicit none
            type (OpenMM_GBSAOBCForce) target
            real*8 dielectric
        end
        function OpenMM_GBSAOBCForce_getSoluteDielectric(target)
            use OpenMM_Types; implicit none
            type (OpenMM_GBSAOBCForce) target
            real*8  OpenMM_GBSAOBCForce_getSoluteDielectric
        end
        subroutine OpenMM_GBSAOBCForce_setSoluteDielectric(target, dielectric)
            use OpenMM_Types; implicit none
            type (OpenMM_GBSAOBCForce) target
            real*8 dielectric
        end

        ! OpenMM::VariableVerletIntegrator
        subroutine OpenMM_VariableVerletIntegrator_create(result, errorTol)
            use OpenMM_Types; implicit none
            type (OpenMM_VariableVerletIntegrator) result
            real*8 errorTol
        end
        subroutine OpenMM_VariableVerletIntegrator_destroy(destroy)
            use OpenMM_Types; implicit none
            type (OpenMM_VariableVerletIntegrator) destroy
        end
        function OpenMM_VariableVerletIntegrator_getErrorTolerance(target)
            use OpenMM_Types; implicit none
            type (OpenMM_VariableVerletIntegrator) target
            real*8  OpenMM_VariableVerletIntegrator_getErrorTolerance
        end
        subroutine OpenMM_VariableVerletIntegrator_setErrorTolerance(target, tol)
            use OpenMM_Types; implicit none
            type (OpenMM_VariableVerletIntegrator) target
            real*8 tol
        end
        subroutine OpenMM_VariableVerletIntegrator_step(target, steps)
            use OpenMM_Types; implicit none
            type (OpenMM_VariableVerletIntegrator) target
            integer*4 steps
        end
        subroutine OpenMM_VariableVerletIntegrator_stepTo(target, time)
            use OpenMM_Types; implicit none
            type (OpenMM_VariableVerletIntegrator) target
            real*8 time
        end

        ! OpenMM::CMMotionRemover
        subroutine OpenMM_CMMotionRemover_create(result, frequency)
            use OpenMM_Types; implicit none
            type (OpenMM_CMMotionRemover) result
            integer*4 frequency
        end
        subroutine OpenMM_CMMotionRemover_destroy(destroy)
            use OpenMM_Types; implicit none
            type (OpenMM_CMMotionRemover) destroy
        end
        function OpenMM_CMMotionRemover_getFrequency(target)
            use OpenMM_Types; implicit none
            type (OpenMM_CMMotionRemover) target
            integer*4  OpenMM_CMMotionRemover_getFrequency
        end
        subroutine OpenMM_CMMotionRemover_setFrequency(target, freq)
            use OpenMM_Types; implicit none
            type (OpenMM_CMMotionRemover) target
            integer*4 freq
        end

        ! OpenMM::VerletIntegrator
        subroutine OpenMM_VerletIntegrator_create(result, stepSize)
            use OpenMM_Types; implicit none
            type (OpenMM_VerletIntegrator) result
            real*8 stepSize
        end
        subroutine OpenMM_VerletIntegrator_destroy(destroy)
            use OpenMM_Types; implicit none
            type (OpenMM_VerletIntegrator) destroy
        end
        subroutine OpenMM_VerletIntegrator_step(target, steps)
            use OpenMM_Types; implicit none
            type (OpenMM_VerletIntegrator) target
            integer*4 steps
        end

        ! OpenMM::RBTorsionForce
        subroutine OpenMM_RBTorsionForce_create(result)
            use OpenMM_Types; implicit none
            type (OpenMM_RBTorsionForce) result
        end
        subroutine OpenMM_RBTorsionForce_destroy(destroy)
            use OpenMM_Types; implicit none
            type (OpenMM_RBTorsionForce) destroy
        end
        function OpenMM_RBTorsionForce_getNumTorsions(target)
            use OpenMM_Types; implicit none
            type (OpenMM_RBTorsionForce) target
            integer*4  OpenMM_RBTorsionForce_getNumTorsions
        end
        function OpenMM_RBTorsionForce_addTorsion(target, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5)
            use OpenMM_Types; implicit none
            type (OpenMM_RBTorsionForce) target
            integer*4 particle1
            integer*4 particle2
            integer*4 particle3
            integer*4 particle4
            real*8 c0
            real*8 c1
            real*8 c2
            real*8 c3
            real*8 c4
            real*8 c5
            integer*4  OpenMM_RBTorsionForce_addTorsion
        end
        subroutine OpenMM_RBTorsionForce_getTorsionParameters(target, index, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5)
            use OpenMM_Types; implicit none
            type (OpenMM_RBTorsionForce) target
            integer*4 index
            integer*4 particle1
            integer*4 particle2
            integer*4 particle3
            integer*4 particle4
            real*8 c0
            real*8 c1
            real*8 c2
            real*8 c3
            real*8 c4
            real*8 c5
        end
        subroutine OpenMM_RBTorsionForce_setTorsionParameters(target, index, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5)
            use OpenMM_Types; implicit none
            type (OpenMM_RBTorsionForce) target
            integer*4 index
            integer*4 particle1
            integer*4 particle2
            integer*4 particle3
            integer*4 particle4
            real*8 c0
            real*8 c1
            real*8 c2
            real*8 c3
            real*8 c4
            real*8 c5
        end

        ! OpenMM::LangevinIntegrator
        subroutine OpenMM_LangevinIntegrator_create(result, temperature, frictionCoeff, stepSize)
            use OpenMM_Types; implicit none
            type (OpenMM_LangevinIntegrator) result
            real*8 temperature
            real*8 frictionCoeff
            real*8 stepSize
        end
        subroutine OpenMM_LangevinIntegrator_destroy(destroy)
            use OpenMM_Types; implicit none
            type (OpenMM_LangevinIntegrator) destroy
        end
        function OpenMM_LangevinIntegrator_getTemperature(target)
            use OpenMM_Types; implicit none
            type (OpenMM_LangevinIntegrator) target
            real*8  OpenMM_LangevinIntegrator_getTemperature
        end
        subroutine OpenMM_LangevinIntegrator_setTemperature(target, temp)
            use OpenMM_Types; implicit none
            type (OpenMM_LangevinIntegrator) target
            real*8 temp
        end
        function OpenMM_LangevinIntegrator_getFriction(target)
            use OpenMM_Types; implicit none
            type (OpenMM_LangevinIntegrator) target
            real*8  OpenMM_LangevinIntegrator_getFriction
        end
        subroutine OpenMM_LangevinIntegrator_setFriction(target, coeff)
            use OpenMM_Types; implicit none
            type (OpenMM_LangevinIntegrator) target
            real*8 coeff
        end
        function OpenMM_LangevinIntegrator_getRandomNumberSeed(target)
            use OpenMM_Types; implicit none
            type (OpenMM_LangevinIntegrator) target
            integer*4  OpenMM_LangevinIntegrator_getRandomNumberSeed
        end
        subroutine OpenMM_LangevinIntegrator_setRandomNumberSeed(target, seed)
            use OpenMM_Types; implicit none
            type (OpenMM_LangevinIntegrator) target
            integer*4 seed
        end
        subroutine OpenMM_LangevinIntegrator_step(target, steps)
            use OpenMM_Types; implicit none
            type (OpenMM_LangevinIntegrator) target
            integer*4 steps
        end

        ! OpenMM::Force

        ! OpenMM::HarmonicAngleForce
        subroutine OpenMM_HarmonicAngleForce_create(result)
            use OpenMM_Types; implicit none
            type (OpenMM_HarmonicAngleForce) result
        end
        subroutine OpenMM_HarmonicAngleForce_destroy(destroy)
            use OpenMM_Types; implicit none
            type (OpenMM_HarmonicAngleForce) destroy
        end
        function OpenMM_HarmonicAngleForce_getNumAngles(target)
            use OpenMM_Types; implicit none
            type (OpenMM_HarmonicAngleForce) target
            integer*4  OpenMM_HarmonicAngleForce_getNumAngles
        end
        function OpenMM_HarmonicAngleForce_addAngle(target, particle1, particle2, particle3, angle, k)
            use OpenMM_Types; implicit none
            type (OpenMM_HarmonicAngleForce) target
            integer*4 particle1
            integer*4 particle2
            integer*4 particle3
            real*8 angle
            real*8 k
            integer*4  OpenMM_HarmonicAngleForce_addAngle
        end
        subroutine OpenMM_HarmonicAngleForce_getAngleParameters(target, index, particle1, particle2, particle3, angle, k)
            use OpenMM_Types; implicit none
            type (OpenMM_HarmonicAngleForce) target
            integer*4 index
            integer*4 particle1
            integer*4 particle2
            integer*4 particle3
            real*8 angle
            real*8 k
        end
        subroutine OpenMM_HarmonicAngleForce_setAngleParameters(target, index, particle1, particle2, particle3, angle, k)
            use OpenMM_Types; implicit none
            type (OpenMM_HarmonicAngleForce) target
            integer*4 index
            integer*4 particle1
            integer*4 particle2
            integer*4 particle3
            real*8 angle
            real*8 k
        end

        ! OpenMM::AndersenThermostat
        subroutine OpenMM_AndersenThermostat_create(result, defaultTemperature, defaultCollisionFrequency)
            use OpenMM_Types; implicit none
            type (OpenMM_AndersenThermostat) result
            real*8 defaultTemperature
            real*8 defaultCollisionFrequency
        end
        subroutine OpenMM_AndersenThermostat_destroy(destroy)
            use OpenMM_Types; implicit none
            type (OpenMM_AndersenThermostat) destroy
        end
        subroutine OpenMM_AndersenThermostat_Temperature(result)
            use OpenMM_Types; implicit none
            character(*) result
        end
        subroutine OpenMM_AndersenThermostat_CollisionFrequency(result)
            use OpenMM_Types; implicit none
            character(*) result
        end
        function OpenMM_AndersenThermostat_getDefaultTemperature(target)
            use OpenMM_Types; implicit none
            type (OpenMM_AndersenThermostat) target
            real*8  OpenMM_AndersenThermostat_getDefaultTemperature
        end
        function OpenMM_AndersenThermostat_getDefaultCollisionFrequency(target)
            use OpenMM_Types; implicit none
            type (OpenMM_AndersenThermostat) target
            real*8  OpenMM_AndersenThermostat_getDefaultCollisionFrequency
        end
        function OpenMM_AndersenThermostat_getRandomNumberSeed(target)
            use OpenMM_Types; implicit none
            type (OpenMM_AndersenThermostat) target
            integer*4  OpenMM_AndersenThermostat_getRandomNumberSeed
        end
        subroutine OpenMM_AndersenThermostat_setRandomNumberSeed(target, seed)
            use OpenMM_Types; implicit none
            type (OpenMM_AndersenThermostat) target
            integer*4 seed
        end

        ! OpenMM::Platform
        subroutine OpenMM_Platform_getName(target, result)
            use OpenMM_Types; implicit none
            type (OpenMM_Platform) target
            character(*) result
        end
        function OpenMM_Platform_getSpeed(target)
            use OpenMM_Types; implicit none
            type (OpenMM_Platform) target
            real*8  OpenMM_Platform_getSpeed
        end
        subroutine OpenMM_Platform_supportsDoublePrecision(target, result)
            use OpenMM_Types; implicit none
            type (OpenMM_Platform) target
            integer*4 result
        end
        subroutine OpenMM_Platform_getPropertyNames(target, result)
            use OpenMM_Types; implicit none
            type (OpenMM_Platform) target
            type (OpenMM_StringArray) result
        end
        subroutine OpenMM_Platform_getPropertyValue(target, context, property, result)
            use OpenMM_Types; implicit none
            type (OpenMM_Platform) target
            type (OpenMM_Context) context
            character(*) property
            character(*) result
        end
        subroutine OpenMM_Platform_setPropertyValue(target, context, property, value)
            use OpenMM_Types; implicit none
            type (OpenMM_Platform) target
            type (OpenMM_Context) context
            character(*) property
            character(*) value
        end
        subroutine OpenMM_Platform_getPropertyDefaultValue(target, property, result)
            use OpenMM_Types; implicit none
            type (OpenMM_Platform) target
            character(*) property
            character(*) result
        end
        subroutine OpenMM_Platform_setPropertyDefaultValue(target, property, value)
            use OpenMM_Types; implicit none
            type (OpenMM_Platform) target
            character(*) property
            character(*) value
        end
        subroutine OpenMM_Platform_contextCreated(target, context)
            use OpenMM_Types; implicit none
            type (OpenMM_Platform) target
            type (OpenMM_ContextImpl) context
        end
        subroutine OpenMM_Platform_contextDestroyed(target, context)
            use OpenMM_Types; implicit none
            type (OpenMM_Platform) target
            type (OpenMM_ContextImpl) context
        end
        subroutine OpenMM_Platform_supportsKernels(target, kernelNames, result)
            use OpenMM_Types; implicit none
            type (OpenMM_Platform) target
            type (OpenMM_StringArray) kernelNames
            integer*4 result
        end
        subroutine OpenMM_Platform_registerPlatform(platform)
            use OpenMM_Types; implicit none
            type (OpenMM_Platform) platform
        end
        function OpenMM_Platform_getNumPlatforms()
            use OpenMM_Types; implicit none
            integer*4  OpenMM_Platform_getNumPlatforms
        end
        subroutine OpenMM_Platform_getPlatform(index, result)
            use OpenMM_Types; implicit none
            integer*4 index
            type (OpenMM_Platform) result
        end
        subroutine OpenMM_Platform_findPlatform(kernelNames, result)
            use OpenMM_Types; implicit none
            type (OpenMM_StringArray) kernelNames
            type (OpenMM_Platform) result
        end
        subroutine OpenMM_Platform_loadPluginLibrary(file)
            use OpenMM_Types; implicit none
            character(*) file
        end
        subroutine OpenMM_Platform_loadPluginsFromDirectory(directory, result)
            use OpenMM_Types; implicit none
            character(*) directory
            type (OpenMM_StringArray) result
        end
        subroutine OpenMM_Platform_getDefaultPluginsDirectory(result)
            use OpenMM_Types; implicit none
            character(*) result
        end

        ! OpenMM::State
        subroutine OpenMM_State_destroy(destroy)
            use OpenMM_Types; implicit none
            type (OpenMM_State) destroy
        end
        function OpenMM_State_getTime(target)
            use OpenMM_Types; implicit none
            type (OpenMM_State) target
            real*8  OpenMM_State_getTime
        end
        subroutine OpenMM_State_getPositions(target, result)
            use OpenMM_Types; implicit none
            type (OpenMM_State) target
            type (OpenMM_Vec3Array) result
        end
        subroutine OpenMM_State_getVelocities(target, result)
            use OpenMM_Types; implicit none
            type (OpenMM_State) target
            type (OpenMM_Vec3Array) result
        end
        subroutine OpenMM_State_getForces(target, result)
            use OpenMM_Types; implicit none
            type (OpenMM_State) target
            type (OpenMM_Vec3Array) result
        end
        function OpenMM_State_getKineticEnergy(target)
            use OpenMM_Types; implicit none
            type (OpenMM_State) target
            real*8  OpenMM_State_getKineticEnergy
        end
        function OpenMM_State_getPotentialEnergy(target)
            use OpenMM_Types; implicit none
            type (OpenMM_State) target
            real*8  OpenMM_State_getPotentialEnergy
        end
        subroutine OpenMM_State_getParameters(target, result)
            use OpenMM_Types; implicit none
            type (OpenMM_State) target
            type (OpenMM_ParameterArray) result
        end

        ! OpenMM::PeriodicTorsionForce
        subroutine OpenMM_PeriodicTorsionForce_create(result)
            use OpenMM_Types; implicit none
            type (OpenMM_PeriodicTorsionForce) result
        end
        subroutine OpenMM_PeriodicTorsionForce_destroy(destroy)
            use OpenMM_Types; implicit none
            type (OpenMM_PeriodicTorsionForce) destroy
        end
        function OpenMM_PeriodicTorsionForce_getNumTorsions(target)
            use OpenMM_Types; implicit none
            type (OpenMM_PeriodicTorsionForce) target
            integer*4  OpenMM_PeriodicTorsionForce_getNumTorsions
        end
        function OpenMM_PeriodicTorsionForce_addTorsion(target, particle1, particle2, particle3, particle4, periodicity, phase, k)
            use OpenMM_Types; implicit none
            type (OpenMM_PeriodicTorsionForce) target
            integer*4 particle1
            integer*4 particle2
            integer*4 particle3
            integer*4 particle4
            integer*4 periodicity
            real*8 phase
            real*8 k
            integer*4  OpenMM_PeriodicTorsionForce_addTorsion
        end
        subroutine OpenMM_PeriodicTorsionForce_getTorsionParameters(target, index, particle1, particle2, particle3, particle4, periodicity, phase, k)
            use OpenMM_Types; implicit none
            type (OpenMM_PeriodicTorsionForce) target
            integer*4 index
            integer*4 particle1
            integer*4 particle2
            integer*4 particle3
            integer*4 particle4
            integer*4 periodicity
            real*8 phase
            real*8 k
        end
        subroutine OpenMM_PeriodicTorsionForce_setTorsionParameters(target, index, particle1, particle2, particle3, particle4, periodicity, phase, k)
            use OpenMM_Types; implicit none
            type (OpenMM_PeriodicTorsionForce) target
            integer*4 index
            integer*4 particle1
            integer*4 particle2
            integer*4 particle3
            integer*4 particle4
            integer*4 periodicity
            real*8 phase
            real*8 k
        end

        ! OpenMM::Integrator
        function OpenMM_Integrator_getStepSize(target)
            use OpenMM_Types; implicit none
            type (OpenMM_Integrator) target
            real*8  OpenMM_Integrator_getStepSize
        end
        subroutine OpenMM_Integrator_setStepSize(target, size)
            use OpenMM_Types; implicit none
            type (OpenMM_Integrator) target
            real*8 size
        end
        function OpenMM_Integrator_getConstraintTolerance(target)
            use OpenMM_Types; implicit none
            type (OpenMM_Integrator) target
            real*8  OpenMM_Integrator_getConstraintTolerance
        end
        subroutine OpenMM_Integrator_setConstraintTolerance(target, tol)
            use OpenMM_Types; implicit none
            type (OpenMM_Integrator) target
            real*8 tol
        end
        subroutine OpenMM_Integrator_step(target, steps)
            use OpenMM_Types; implicit none
            type (OpenMM_Integrator) target
            integer*4 steps
        end

        ! OpenMM::System
        subroutine OpenMM_System_create(result)
            use OpenMM_Types; implicit none
            type (OpenMM_System) result
        end
        subroutine OpenMM_System_destroy(destroy)
            use OpenMM_Types; implicit none
            type (OpenMM_System) destroy
        end
        function OpenMM_System_getNumParticles(target)
            use OpenMM_Types; implicit none
            type (OpenMM_System) target
            integer*4  OpenMM_System_getNumParticles
        end
        function OpenMM_System_addParticle(target, mass)
            use OpenMM_Types; implicit none
            type (OpenMM_System) target
            real*8 mass
            integer*4  OpenMM_System_addParticle
        end
        function OpenMM_System_getParticleMass(target, index)
            use OpenMM_Types; implicit none
            type (OpenMM_System) target
            integer*4 index
            real*8  OpenMM_System_getParticleMass
        end
        subroutine OpenMM_System_setParticleMass(target, index, mass)
            use OpenMM_Types; implicit none
            type (OpenMM_System) target
            integer*4 index
            real*8 mass
        end
        function OpenMM_System_getNumConstraints(target)
            use OpenMM_Types; implicit none
            type (OpenMM_System) target
            integer*4  OpenMM_System_getNumConstraints
        end
        function OpenMM_System_addConstraint(target, particle1, particle2, distance)
            use OpenMM_Types; implicit none
            type (OpenMM_System) target
            integer*4 particle1
            integer*4 particle2
            real*8 distance
            integer*4  OpenMM_System_addConstraint
        end
        subroutine OpenMM_System_getConstraintParameters(target, index, particle1, particle2, distance)
            use OpenMM_Types; implicit none
            type (OpenMM_System) target
            integer*4 index
            integer*4 particle1
            integer*4 particle2
            real*8 distance
        end
        subroutine OpenMM_System_setConstraintParameters(target, index, particle1, particle2, distance)
            use OpenMM_Types; implicit none
            type (OpenMM_System) target
            integer*4 index
            integer*4 particle1
            integer*4 particle2
            real*8 distance
        end
        function OpenMM_System_addForce(target, force)
            use OpenMM_Types; implicit none
            type (OpenMM_System) target
            type (OpenMM_Force) force
            integer*4  OpenMM_System_addForce
        end
        function OpenMM_System_getNumForces(target)
            use OpenMM_Types; implicit none
            type (OpenMM_System) target
            integer*4  OpenMM_System_getNumForces
        end
        subroutine OpenMM_System_getForce(target, index, result)
            use OpenMM_Types; implicit none
            type (OpenMM_System) target
            integer*4 index
            type (OpenMM_Force) result
        end
    end interface
END MODULE OpenMM
