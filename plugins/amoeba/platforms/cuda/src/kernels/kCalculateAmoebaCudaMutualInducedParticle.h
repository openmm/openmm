
struct MutualInducedParticle {

    float x;
    float y;
    float z;

    float inducedDipole[3];
    float inducedDipolePolar[3];

    float thole;
    float damp;

    float field[3];
    float fieldPolar[3];

#ifdef GK
    float bornRadius;

    float inducedDipoleS[3];
    float inducedDipolePolarS[3];

    float fieldS[3];
    float fieldPolarS[3];
#else
//    float padding;
#endif

#ifdef INCLUDE_MI_FIELD_BUFFERS
    float tempBuffer[3];
    float tempBufferP[3];
#endif
};

__device__ static void loadMutualInducedShared( MutualInducedParticle* sA, unsigned int atomI )
{
    // coordinates & charge

    float4 posq                  = cSim.pPosq[atomI];
    sA->x                        = posq.x;
    sA->y                        = posq.y;
    sA->z                        = posq.z;

    // dipole

    sA->inducedDipole[0]         = cAmoebaSim.pInducedDipole[atomI*3];
    sA->inducedDipole[1]         = cAmoebaSim.pInducedDipole[atomI*3+1];
    sA->inducedDipole[2]         = cAmoebaSim.pInducedDipole[atomI*3+2];

    // dipole polar

    sA->inducedDipolePolar[0]    = cAmoebaSim.pInducedDipolePolar[atomI*3];
    sA->inducedDipolePolar[1]    = cAmoebaSim.pInducedDipolePolar[atomI*3+1];
    sA->inducedDipolePolar[2]    = cAmoebaSim.pInducedDipolePolar[atomI*3+2];

    float2 dampingFactorAndThole = cAmoebaSim.pDampingFactorAndThole[atomI];
    sA->damp                     = dampingFactorAndThole.x;
    sA->thole                    = dampingFactorAndThole.y;

#ifdef GK

    sA->bornRadius               =  cSim.pBornRadii[atomI];

    // dipoleS

    sA->inducedDipoleS[0]        = cAmoebaSim.pInducedDipoleS[atomI*3];
    sA->inducedDipoleS[1]        = cAmoebaSim.pInducedDipoleS[atomI*3+1];
    sA->inducedDipoleS[2]        = cAmoebaSim.pInducedDipoleS[atomI*3+2];

    // dipole polar S

    sA->inducedDipolePolarS[0]   = cAmoebaSim.pInducedDipolePolarS[atomI*3];
    sA->inducedDipolePolarS[1]   = cAmoebaSim.pInducedDipolePolarS[atomI*3+1];
    sA->inducedDipolePolarS[2]   = cAmoebaSim.pInducedDipolePolarS[atomI*3+2];

#endif
}

__device__ static void zeroMutualInducedParticleSharedField( MutualInducedParticle* sA )

{
    // zero shared fields

    sA->field[0]              = 0.0f;
    sA->field[1]              = 0.0f;
    sA->field[2]              = 0.0f;

    sA->fieldPolar[0]         = 0.0f;
    sA->fieldPolar[1]         = 0.0f;
    sA->fieldPolar[2]         = 0.0f;

#ifdef GK
    sA->fieldS[0]             = 0.0f;
    sA->fieldS[1]             = 0.0f;
    sA->fieldS[2]             = 0.0f;

    sA->fieldPolarS[0]        = 0.0f;
    sA->fieldPolarS[1]        = 0.0f;
    sA->fieldPolarS[2]        = 0.0f;
#endif

}
