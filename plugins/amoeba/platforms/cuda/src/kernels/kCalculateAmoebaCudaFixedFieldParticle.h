
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

#ifdef INCLUDE_FIXED_FIELD_BUFFERS
    float tempBuffer[3];
    float tempBufferP[3];
#endif
};

__device__ static void loadFixedFieldShared( struct FixedFieldParticle* sA, unsigned int atomI 
#ifdef GK
    , float* bornR
#endif
)
{
    // coordinates & charge

    float4 posq                  = cSim.pPosq[atomI];
    sA->x                        = posq.x;
    sA->y                        = cSim.pPosq[atomI].y;
    sA->z                        = cSim.pPosq[atomI].z;
    sA->q                        = cSim.pPosq[atomI].w;

    // lab dipole

    sA->labFrameDipole_X         = cAmoebaSim.pLabFrameDipole[atomI*3];
    sA->labFrameDipole_Y         = cAmoebaSim.pLabFrameDipole[atomI*3+1];
    sA->labFrameDipole_Z         = cAmoebaSim.pLabFrameDipole[atomI*3+2];

    // lab quadrupole

    sA->labFrameQuadrupole_XX    = cAmoebaSim.pLabFrameQuadrupole[atomI*9];
    sA->labFrameQuadrupole_XY    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+1];
    sA->labFrameQuadrupole_XZ    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+2];
    sA->labFrameQuadrupole_YY    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+4];
    sA->labFrameQuadrupole_YZ    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+5];
    sA->labFrameQuadrupole_ZZ    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+8];

    float2 dampingFactorAndThole = cAmoebaSim.pDampingFactorAndThole[atomI];
    sA->damp                     = dampingFactorAndThole.x;
    sA->thole                    = dampingFactorAndThole.y;
#ifdef GK
    sA->bornR                    = bornR[atomI];
#endif

}

// load struct and arrays w/ shared data in sA

__device__ static void loadFixedFieldParticleData( struct FixedFieldParticle* sA, 
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

__device__ static void zeroFixedFieldParticleSharedField( struct FixedFieldParticle* sA )
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

__device__ static void calculateFixedEFieldPairIxn_kernel( FixedFieldParticle& atomI, FixedFieldParticle& atomJ,
                                                           float field[2][3]
#ifdef AMOEBA_DEBUG
                                                    , float4 debugArray[12]
#endif
)
{

 
    // ---------------------------------------------------------------------------------------
 
    // get deltaR and r between 2 atoms
 
    float deltaR[3];
    deltaR[0]           = atomJ.x - atomI.x;
    deltaR[1]           = atomJ.y - atomI.y;
    deltaR[2]           = atomJ.z - atomI.z;

    float r             =  SQRT( deltaR[0]*deltaR[0] + deltaR[1]*deltaR[1] + deltaR[2]*deltaR[2] );
    float rI            =  1.0f/r;
    float r2I           =  rI*rI;

    float rr3           =  rI*r2I;
    float rr5           =  3.0f*rr3*r2I;
    float rr7           =  5.0f*rr5*r2I;
 
    // get scaling factors, if needed
    
    float damp          = atomI.damp*atomJ.damp;
    float dampExp;
    if( damp != 0.0f && r < cAmoebaSim.scalingDistanceCutoff ){

        // get scaling factors
      
        float ratio     = r/damp;
        float pGamma    = atomJ.thole > atomI.thole ? atomI.thole : atomJ.thole; 
        damp            = ratio*ratio*ratio*pGamma;
        dampExp         = EXP( -damp );
    } else {
        dampExp         = 0.0f;
    }
      
    rr3                *= 1.0f - dampExp;
    rr5                *= 1.0f - ( 1.0f + damp )*dampExp;
    rr7                *= 1.0f - ( 1.0f + damp + (0.6f*damp*damp))*dampExp;
      
    float rr5_2         = rr5*2.0f;
 
    float  qDotDelta[3];
    qDotDelta[0]        = deltaR[0]*atomJ.labFrameQuadrupole_XX + deltaR[1]*atomJ.labFrameQuadrupole_XY + deltaR[2]*atomJ.labFrameQuadrupole_XZ;
    qDotDelta[1]        = deltaR[0]*atomJ.labFrameQuadrupole_XY + deltaR[1]*atomJ.labFrameQuadrupole_YY + deltaR[2]*atomJ.labFrameQuadrupole_YZ;
    qDotDelta[2]        = deltaR[0]*atomJ.labFrameQuadrupole_XZ + deltaR[1]*atomJ.labFrameQuadrupole_YZ + deltaR[2]*atomJ.labFrameQuadrupole_ZZ;
 
    float dotdd         = deltaR[0]*atomJ.labFrameDipole_X      + deltaR[1]*atomJ.labFrameDipole_Y      + deltaR[2]*atomJ.labFrameDipole_Z;
    float dotqd         = deltaR[0]*qDotDelta[0]                + deltaR[1]*qDotDelta[1]                + deltaR[2]*qDotDelta[2];

    float factor        = -rr3*atomJ.q + rr5*dotdd - rr7*dotqd;
 
    field[0][0]         = deltaR[0]*factor - rr3*atomJ.labFrameDipole_X + rr5_2*qDotDelta[0];
    field[0][1]         = deltaR[1]*factor - rr3*atomJ.labFrameDipole_Y + rr5_2*qDotDelta[1];
    field[0][2]         = deltaR[2]*factor - rr3*atomJ.labFrameDipole_Z + rr5_2*qDotDelta[2];
 
    qDotDelta[0]        = deltaR[0]*atomI.labFrameQuadrupole_XX + deltaR[1]*atomI.labFrameQuadrupole_XY + deltaR[2]*atomI.labFrameQuadrupole_XZ;
    qDotDelta[1]        = deltaR[0]*atomI.labFrameQuadrupole_XY + deltaR[1]*atomI.labFrameQuadrupole_YY + deltaR[2]*atomI.labFrameQuadrupole_YZ;
    qDotDelta[2]        = deltaR[0]*atomI.labFrameQuadrupole_XZ + deltaR[1]*atomI.labFrameQuadrupole_YZ + deltaR[2]*atomI.labFrameQuadrupole_ZZ;
 
    dotdd               = deltaR[0]*atomI.labFrameDipole_X    + deltaR[1]*atomI.labFrameDipole_Y    + deltaR[2]*atomI.labFrameDipole_Z;
    dotqd               = deltaR[0]*qDotDelta[0] + deltaR[1]*qDotDelta[1] + deltaR[2]*qDotDelta[2];
    factor              = rr3*atomI.q + rr5*dotdd + rr7*dotqd;
 
    field[1][0]         = deltaR[0]*factor - rr3*atomI.labFrameDipole_X - rr5_2*qDotDelta[0];
    field[1][1]         = deltaR[1]*factor - rr3*atomI.labFrameDipole_Y - rr5_2*qDotDelta[1];
    field[1][2]         = deltaR[2]*factor - rr3*atomI.labFrameDipole_Z - rr5_2*qDotDelta[2];

}
