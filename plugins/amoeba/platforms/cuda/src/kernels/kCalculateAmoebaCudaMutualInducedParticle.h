
struct MutualInducedParticle {

    float x;
    float y;
    float z;
    float q;

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
#endif
};

__device__ void loadMutualInducedShared( struct MutualInducedParticle* sA, unsigned int atomI,
                                         float4* atomCoord,
                                         float* inducedDipole, float* inducedDipolePolar,
                                         float2* dampingFactorAndThole 
#ifdef GK
                                         , float* bornRadii, float* inducedDipoleS, float* inducedDipolePolarS
#endif
)

{
    // coordinates & charge

    sA->x                        = atomCoord[atomI].x;
    sA->y                        = atomCoord[atomI].y;
    sA->z                        = atomCoord[atomI].z;
    sA->q                        = atomCoord[atomI].w;

    // dipole

    sA->inducedDipole[0]         = inducedDipole[atomI*3];
    sA->inducedDipole[1]         = inducedDipole[atomI*3+1];
    sA->inducedDipole[2]         = inducedDipole[atomI*3+2];

    // dipole polar

    sA->inducedDipolePolar[0]    = inducedDipolePolar[atomI*3];
    sA->inducedDipolePolar[1]    = inducedDipolePolar[atomI*3+1];
    sA->inducedDipolePolar[2]    = inducedDipolePolar[atomI*3+2];

    sA->damp                     = dampingFactorAndThole[atomI].x;
    sA->thole                    = dampingFactorAndThole[atomI].y;

#ifdef GK

    sA->bornRadius               = bornRadii[atomI];

    // dipoleS

    sA->inducedDipoleS[0]        = inducedDipoleS[atomI*3];
    sA->inducedDipoleS[1]        = inducedDipoleS[atomI*3+1];
    sA->inducedDipoleS[2]        = inducedDipoleS[atomI*3+2];

    // dipole polar S

    sA->inducedDipolePolarS[0]   = inducedDipolePolarS[atomI*3];
    sA->inducedDipolePolarS[1]   = inducedDipolePolarS[atomI*3+1];
    sA->inducedDipolePolarS[2]   = inducedDipolePolarS[atomI*3+2];

#endif
}

__device__ void zeroMutualInducedParticleSharedField( struct MutualInducedParticle* sA )

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

// load struct and arrays w/ shared data in sA

__device__ void loadMutualInducedData( struct MutualInducedParticle* sA,  float4* jCoord,
                                       float* jInducedDipole,  float* jInducedDipolePolar
#ifdef GK
                                      , float* jBornRadius, float* jInducedDipoleS, float* jInducedDipolePolarS
#endif
)
{

    // load coords, charge, ...

    jCoord->x                = sA->x;
    jCoord->y                = sA->y;
    jCoord->z                = sA->z;
    jCoord->w                = sA->q;
 
    jInducedDipole[0]        = sA->inducedDipole[0];
    jInducedDipole[1]        = sA->inducedDipole[1];
    jInducedDipole[2]        = sA->inducedDipole[2];
 
    jInducedDipolePolar[0]   = sA->inducedDipolePolar[0];
    jInducedDipolePolar[1]   = sA->inducedDipolePolar[1];
    jInducedDipolePolar[2]   = sA->inducedDipolePolar[2];

#ifdef GK
    *jBornRadius             = sA->bornRadius;

    jInducedDipoleS[0]       = sA->inducedDipoleS[0];
    jInducedDipoleS[1]       = sA->inducedDipoleS[1];
    jInducedDipoleS[2]       = sA->inducedDipoleS[2];
 
    jInducedDipolePolarS[0]  = sA->inducedDipolePolarS[0];
    jInducedDipolePolarS[1]  = sA->inducedDipolePolarS[1];
    jInducedDipolePolarS[2]  = sA->inducedDipolePolarS[2];

#endif
 
}

