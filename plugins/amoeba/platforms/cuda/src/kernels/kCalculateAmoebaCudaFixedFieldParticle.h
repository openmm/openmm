
struct FixedFieldParticle {

    // coordinates charge

    float x;
    float y;
    float z;
    float q;

    // lab frame dipole

    float labFrameDipole_X;
    float labFrameDipole_Y;
    float labFrameDipole_Z;

    // lab frame quadrupole

    float labFrameQuadrupole_XX;
    float labFrameQuadrupole_XY;
    float labFrameQuadrupole_XZ;
    float labFrameQuadrupole_YY;
    float labFrameQuadrupole_YZ;
    float labFrameQuadrupole_ZZ;

    // scaling factor

    float thole;
    float damp;

    // field accumulators

    float eField[3];
    float eFieldP[3];

#ifdef GK

    // Born radius

    float bornR;

    // GK field

    float gkField[3];

#endif
};

__device__ void loadFixedFieldShared( struct FixedFieldParticle* sA, unsigned int atomI,
                                      float4* atomCoord, float* labDipole, float* labQuadrupole,
                                      float2* dampingFactorAndThole
#ifdef GK
    , float* bornR
#endif
)
{
    // coordinates & charge

    sA->x                        = atomCoord[atomI].x;
    sA->y                        = atomCoord[atomI].y;
    sA->z                        = atomCoord[atomI].z;
    sA->q                        = atomCoord[atomI].w;

    // lab dipole

    sA->labFrameDipole_X         = labDipole[atomI*3];
    sA->labFrameDipole_Y         = labDipole[atomI*3+1];
    sA->labFrameDipole_Z         = labDipole[atomI*3+2];

    // lab quadrupole

    sA->labFrameQuadrupole_XX    = labQuadrupole[atomI*9];
    sA->labFrameQuadrupole_XY    = labQuadrupole[atomI*9+1];
    sA->labFrameQuadrupole_XZ    = labQuadrupole[atomI*9+2];
    sA->labFrameQuadrupole_YY    = labQuadrupole[atomI*9+4];
    sA->labFrameQuadrupole_YZ    = labQuadrupole[atomI*9+5];
    sA->labFrameQuadrupole_ZZ    = labQuadrupole[atomI*9+8];

    sA->damp                     = dampingFactorAndThole[atomI].x;
    sA->thole                    = dampingFactorAndThole[atomI].y;
#ifdef GK
    sA->bornR                    = bornR[atomI];
#endif

}

// load struct and arrays w/ shared data in sA

__device__ void loadFixedFieldParticleData( struct FixedFieldParticle* sA, 
                                            float4* jCoord, float* jDipole, float* jQuadrupole
#ifdef GK
, float* bornR
#endif
)
{

    // load coords, charge, ...

    jCoord->x               = sA->x;
    jCoord->y               = sA->y;
    jCoord->z               = sA->z;
    jCoord->w               = sA->q;

    jDipole[0]              = sA->labFrameDipole_X;
    jDipole[1]              = sA->labFrameDipole_Y;
    jDipole[2]              = sA->labFrameDipole_Z;

    jQuadrupole[0]          = sA->labFrameQuadrupole_XX;
    jQuadrupole[1]          = sA->labFrameQuadrupole_XY;
    jQuadrupole[2]          = sA->labFrameQuadrupole_XZ;

    jQuadrupole[3]          = sA->labFrameQuadrupole_XY;
    jQuadrupole[4]          = sA->labFrameQuadrupole_YY;
    jQuadrupole[5]          = sA->labFrameQuadrupole_YZ;

    jQuadrupole[6]          = sA->labFrameQuadrupole_XZ;
    jQuadrupole[7]          = sA->labFrameQuadrupole_YZ;
    jQuadrupole[8]          = sA->labFrameQuadrupole_ZZ;
 
#ifdef GK
    *bornR                  = sA->bornR;
#endif
}

// zero fields

__device__ void zeroFixedFieldParticleSharedField( struct FixedFieldParticle* sA ) 
{

    sA->eField[0]    = 0.0f;
    sA->eField[1]    = 0.0f;
    sA->eField[2]    = 0.0f;

    sA->eFieldP[0]   = 0.0f;
    sA->eFieldP[1]   = 0.0f;
    sA->eFieldP[2]   = 0.0f;

#ifdef GK
    sA->gkField[0]   = 0.0f;
    sA->gkField[1]   = 0.0f;
    sA->gkField[2]   = 0.0f;
#endif
}

// body of fixed E-field calculation

__device__ void calculateFixedEFieldPairIxn_kernel( float4 atomCoordinatesI, float4 atomCoordinatesJ,
                                                    float dampingFactorI,    float dampingFactorJ,
                                                    float tholeI,            float tholeJ,
                                                    float* labDipoleI, float* labDipoleJ,
                                                    float* labQuadrupoleI, float* labQuadrupoleJ,
                                                    float scalingDistanceCutoff,
                                                    float field[2][3]
#ifdef AMOEBA_DEBUG
                                                    , float4 debugArray[12]
#endif
)
{

 
    // ---------------------------------------------------------------------------------------
 
    // get deltaR and r between 2 atoms
 
    float deltaR[3];
    deltaR[0]           = atomCoordinatesJ.x - atomCoordinatesI.x;
    deltaR[1]           = atomCoordinatesJ.y - atomCoordinatesI.y;
    deltaR[2]           = atomCoordinatesJ.z - atomCoordinatesI.z;

    float r             =  SQRT( deltaR[0]*deltaR[0] + deltaR[1]*deltaR[1] + deltaR[2]*deltaR[2] );
    float rI            =  1.0f/r;
    float r2I           =  rI*rI;

    float rr3           =  rI*r2I;
    float rr5           =  3.0f*rr3*r2I;
    float rr7           =  5.0f*rr5*r2I;
 
    // get scaling factors, if needed
    
    float damp          = dampingFactorI*dampingFactorJ;
    float dampExp;
    if( damp != 0.0f && r < scalingDistanceCutoff ){

        // get scaling factors
      
        float ratio     = r/damp;
        float pGamma    = tholeJ > tholeI ? tholeI : tholeJ; 
        damp            = ratio*ratio*ratio*pGamma;
        dampExp         = EXP( -damp );
    } else {
        dampExp         = 0.0f;
    }
      
    rr3                *= 1.0f - dampExp;
    rr5                *= 1.0f - ( 1.0f + damp )*dampExp;
    rr7                *= 1.0f - ( 1.0f + damp + (0.6f*damp*damp))*dampExp;
      
    float rr5_2         = rr5*2.0f;
 
#ifdef AMOEBA_DEBUG
    int index           = 0;

    // 0-2
    debugArray[index].x   = r;
    debugArray[index].y   = rr3;
    debugArray[index].z   = rr5;
    index++;
#endif

    float* dipole       = labDipoleJ;
    float* quadrupole   = labQuadrupoleJ;
    float  qDotDelta[3];
    qDotDelta[0]        = deltaR[0]*quadrupole[0] + deltaR[1]*quadrupole[1] + deltaR[2]*quadrupole[2];
    qDotDelta[1]        = deltaR[0]*quadrupole[3] + deltaR[1]*quadrupole[4] + deltaR[2]*quadrupole[5];
    qDotDelta[2]        = deltaR[0]*quadrupole[6] + deltaR[1]*quadrupole[7] + deltaR[2]*quadrupole[8];
 
    float dotdd         = deltaR[0]*dipole[0]    + deltaR[1]*dipole[1]    + deltaR[2]*dipole[2];
    float dotqd         = deltaR[0]*qDotDelta[0] + deltaR[1]*qDotDelta[1] + deltaR[2]*qDotDelta[2];
    float factor        = -rr3*atomCoordinatesJ.w + rr5*dotdd - rr7*dotqd;
 
#ifdef AMOEBA_DEBUG
    // 3-5
    debugArray[index].x   = dotdd;
    debugArray[index].y   = dotqd;
    debugArray[index].z   = factor;
    index++;
#endif

    field[0][0]         = deltaR[0]*factor - rr3*dipole[0] + rr5_2*qDotDelta[0];
    field[0][1]         = deltaR[1]*factor - rr3*dipole[1] + rr5_2*qDotDelta[1];
    field[0][2]         = deltaR[2]*factor - rr3*dipole[2] + rr5_2*qDotDelta[2];
 
    dipole              = labDipoleI;
    quadrupole          = labQuadrupoleI;
    qDotDelta[0]        = deltaR[0]*quadrupole[0] + deltaR[1]*quadrupole[1] + deltaR[2]*quadrupole[2];
    qDotDelta[1]        = deltaR[0]*quadrupole[3] + deltaR[1]*quadrupole[4] + deltaR[2]*quadrupole[5];
    qDotDelta[2]        = deltaR[0]*quadrupole[6] + deltaR[1]*quadrupole[7] + deltaR[2]*quadrupole[8];
 
    dotdd               = deltaR[0]*dipole[0]    + deltaR[1]*dipole[1]    + deltaR[2]*dipole[2];
    dotqd               = deltaR[0]*qDotDelta[0] + deltaR[1]*qDotDelta[1] + deltaR[2]*qDotDelta[2];
    factor              = rr3*atomCoordinatesI.w + rr5*dotdd + rr7*dotqd;
 
#ifdef AMOEBA_DEBUG
    // 6-8
    debugArray[index].x = dotdd;
    debugArray[index].y = dotqd;
    debugArray[index].z = factor;
    index++;
#endif

    field[1][0]         = deltaR[0]*factor - rr3*dipole[0] - rr5_2*qDotDelta[0];
    field[1][1]         = deltaR[1]*factor - rr3*dipole[1] - rr5_2*qDotDelta[1];
    field[1][2]         = deltaR[2]*factor - rr3*dipole[2] - rr5_2*qDotDelta[2];

#if 0
    float testValue  = 1.0f;
    field[0][0]      = testValue;
    field[0][1]      = testValue;
    field[0][2]      = testValue;
    field[1][0]      = testValue;
    field[1][1]      = testValue;
    field[1][2]      = testValue;
#endif
 
}
