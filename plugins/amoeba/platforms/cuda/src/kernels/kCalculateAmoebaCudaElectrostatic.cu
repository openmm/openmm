//-----------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------

#include "amoebaGpuTypes.h"
#include "amoebaCudaKernels.h"
#include "kCalculateAmoebaCudaUtilities.h"

//#define AMOEBA_DEBUG

static __constant__ cudaGmxSimulation cSim;
static __constant__ cudaAmoebaGmxSimulation cAmoebaSim;

void SetCalculateAmoebaElectrostaticSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status         = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaElectrostaticSim: cudaMemcpyToSymbol: SetSim copy to cSim failed");
    status         = cudaMemcpyToSymbol(cAmoebaSim, &amoebaGpu->amoebaSim, sizeof(cudaAmoebaGmxSimulation));    
    RTERROR(status, "SetCalculateAmoebaElectrostaticSim: cudaMemcpyToSymbol: SetSim copy to cAmoebaSim failed");
}

void GetCalculateAmoebaElectrostaticSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));    
    RTERROR(status, "GetCalculateAmoebaElectrostaticSim: cudaMemcpyFromSymbol: SetSim copy from cSim failed");
    status = cudaMemcpyFromSymbol(&amoebaGpu->amoebaSim, cAmoebaSim, sizeof(cudaAmoebaGmxSimulation));    
    RTERROR(status, "GetCalculateAmoebaElectrostaticSim: cudaMemcpyFromSymbol: SetSim copy from cAmoebaSim failed");
}

static int const PScaleIndex            =  0; 
static int const DScaleIndex            =  1; 
static int const UScaleIndex            =  2; 
static int const MScaleIndex            =  3;
static int const Scale3Index            =  4;
static int const Scale5Index            =  5;
static int const Scale7Index            =  6;
static int const Scale9Index            =  7;
static int const Ddsc30Index            =  8;
//static int const Ddsc31Index            =  9;
//static int const Ddsc32Index            = 10; 
static int const Ddsc50Index            = 11;
//static int const Ddsc51Index            = 12;
//static int const Ddsc52Index            = 13; 
static int const Ddsc70Index            = 14;
//static int const Ddsc71Index            = 15;
//static int const Ddsc72Index            = 16;
static int const LastScalingIndex       = 17;

#define DOT3_4(u,v) ((u[0])*(v[0]) + (u[1])*(v[1]) + (u[2])*(v[2]))

#define MATRIXDOT31(u,v) u[0]*v[0] + u[1]*v[1] + u[2]*v[2] + \
  u[3]*v[3] + u[4]*v[4] + u[5]*v[5] + \
  u[6]*v[6] + u[7]*v[7] + u[8]*v[8]

#define DOT31(u,v) ((u[0])*(v[0]) + (u[1])*(v[1]) + (u[2])*(v[2]))

#define i35 0.257142857f
#define one 1.0f

__device__ void acrossProductVector3(   float* vectorX, float* vectorY, float* vectorZ ){
    vectorZ[0]  = vectorX[1]*vectorY[2] - vectorX[2]*vectorY[1];
    vectorZ[1]  = vectorX[2]*vectorY[0] - vectorX[0]*vectorY[2];
    vectorZ[2]  = vectorX[0]*vectorY[1] - vectorX[1]*vectorY[0];
}

__device__ void amatrixProductVector3(   float* matrixX, float* vectorY, float* vectorZ ){
    vectorZ[0]  = matrixX[0]*vectorY[0] + matrixX[3]*vectorY[1] + matrixX[6]*vectorY[2];
    vectorZ[1]  = matrixX[1]*vectorY[0] + matrixX[4]*vectorY[1] + matrixX[7]*vectorY[2];
    vectorZ[2]  = matrixX[2]*vectorY[0] + matrixX[5]*vectorY[1] + matrixX[8]*vectorY[2];
}

__device__ void amatrixCrossProductMatrix3( float* matrixX, float* matrixY, float* vectorZ ){
  
    float* xPtr[3];
    float* yPtr[3];
        
    xPtr[0]    = matrixX;
    xPtr[1]    = matrixX + 3;
    xPtr[2]    = matrixX + 6;
    
    yPtr[0]    = matrixY;
    yPtr[1]    = matrixY + 3;
    yPtr[2]    = matrixY + 6;
          
    vectorZ[0] = DOT31( xPtr[1], yPtr[2] ) - DOT31( xPtr[2], yPtr[1] );
    vectorZ[1] = DOT31( xPtr[2], yPtr[0] ) - DOT31( xPtr[0], yPtr[2] );
    vectorZ[2] = DOT31( xPtr[0], yPtr[1] ) - DOT31( xPtr[1], yPtr[0] );
  
}

struct ElectrostaticParticle {

    // coordinates charge

    float x;
    float y;
    float z;
    float q;

    // lab frame dipole

    float labFrameDipole[3];

    // lab frame quadrupole

    float labFrameQuadrupole[9];

    // induced dipole

    float inducedDipole[3];

    // polar induced dipole

    float inducedDipoleP[3];

    // scaling factors

    float thole;
    float damp;

    float force[3];

    float torque[3];
    float padding;

};

__device__ void calculateElectrostaticPairIxn_kernel( ElectrostaticParticle& atomI,   ElectrostaticParticle& atomJ,
                                                      float* scalingFactors, float4*  outputForce, float4  outputTorque[2]
#ifdef AMOEBA_DEBUG
                                                      ,float4* debugArray 
#endif
 ){
  
    float deltaR[3];
    
    // ---------------------------------------------------------------------------------------
    
    // ---------------------------------------------------------------------------------------

    float* ddsc3                    =  scalingFactors + Ddsc30Index;
    float* ddsc5                    =  scalingFactors + Ddsc50Index;
    float* ddsc7                    =  scalingFactors + Ddsc70Index;

    deltaR[0]                       = atomJ.x - atomI.x;
    deltaR[1]                       = atomJ.y - atomI.y;
    deltaR[2]                       = atomJ.z - atomI.z;

    float r2                        = DOT31( deltaR, deltaR );
    float r                         = sqrtf( r2 );
    float rr1                       = 1.0f/r;
    float rr2                       = rr1*rr1;
    float rr3                       = rr1*rr2;
    float rr5                       = 3.0f*rr3*rr2;
    float rr7                       = 5.0f*rr5*rr2;
    float rr9                       = 7.0f*rr7*rr2;
    float rr11                      = 9.0f*rr9*rr2;

    //-------------------------------------------

    if( atomI.damp != 0.0f && atomJ.damp != 0.0 && r < cAmoebaSim.scalingDistanceCutoff ){
   
        float distanceIJ, r2I;
        distanceIJ                    = r;
        r2I                           = rr2;
        
        float ratio                   = distanceIJ/(atomI.damp*atomJ.damp);
        float pGamma                  = atomJ.thole > atomI.thole ? atomI.thole : atomJ.thole;

        float damp                          = ratio*ratio*ratio*pGamma;
        float dampExp                 = expf( -damp );
        float damp1                   = damp + one;
        float damp2                   = damp*damp;
        float damp3                   = damp2*damp;

        scalingFactors[Scale3Index]   = one - dampExp;
        scalingFactors[Scale5Index]   = one - damp1*dampExp;
        scalingFactors[Scale7Index]   = one - ( damp1 + 0.6f*damp2)*dampExp;
        scalingFactors[Scale9Index]   = one - ( damp1 + ( 2.0f*damp2 + damp3 )*i35)*dampExp;

        float factor                  = 3.0f*damp*dampExp*r2I;
        float factor7                 = -0.2f + 0.6f*damp;
        
        for( int ii = 0; ii < 3; ii++ ){
            scalingFactors[Ddsc30Index + ii] = factor*deltaR[ii];
            scalingFactors[Ddsc50Index + ii] = scalingFactors[Ddsc30Index + ii]*damp;
            scalingFactors[Ddsc70Index + ii] = scalingFactors[Ddsc50Index + ii]*factor7;
        }

    }
      
    float scaleI0 = scalingFactors[Scale3Index]*scalingFactors[UScaleIndex];
    float dsc0    = scalingFactors[Scale3Index]*scalingFactors[DScaleIndex];
    float psc0    = scalingFactors[Scale3Index]*scalingFactors[PScaleIndex];
    float scaleI1 = scalingFactors[Scale3Index+1]*scalingFactors[UScaleIndex];
    float dsc1    = scalingFactors[Scale3Index+1]*scalingFactors[DScaleIndex];
    float psc1    = scalingFactors[Scale3Index+1]*scalingFactors[PScaleIndex];
    float dsc2    = scalingFactors[Scale3Index+2]*scalingFactors[DScaleIndex];
    float psc2    = scalingFactors[Scale3Index+2]*scalingFactors[PScaleIndex];
                       
    float qIr[3], qJr[3];

    amatrixProductVector3( atomJ.labFrameQuadrupole,      deltaR,      qJr);
    amatrixProductVector3( atomI.labFrameQuadrupole,      deltaR,      qIr);

    float sc2     = DOT3_4(        atomI.labFrameDipole,  atomJ.labFrameDipole );
    float sc3     = DOT3_4(        atomI.labFrameDipole,  deltaR  );
    float sc4     = DOT3_4(        atomJ.labFrameDipole,  deltaR  );
    
    float sc5     = DOT3_4(        qIr, deltaR  );
    float sc6     = DOT3_4(        qJr, deltaR  );
    
    float sc7     = DOT3_4(        qIr, atomJ.labFrameDipole );
    float sc8     = DOT3_4(        qJr, atomI.labFrameDipole );
    
    float sc9     = DOT3_4(        qIr, qJr );
    
    float sc10    = MATRIXDOT31( atomI.labFrameQuadrupole, atomJ.labFrameQuadrupole );
    
    float sci1    = DOT3_4(        atomI.inducedDipole,  atomJ.labFrameDipole ) +
                    DOT3_4(        atomJ.inducedDipole,  atomI.labFrameDipole );
        
    float sci3    = DOT3_4(        atomI.inducedDipole,  deltaR  );
    float sci4    = DOT3_4(        atomJ.inducedDipole,  deltaR  );
    
    float sci7    = DOT3_4(        qIr, atomJ.inducedDipole );
    float sci8    = DOT3_4(        qJr, atomI.inducedDipole );
    
    float scip1   = DOT3_4(        atomI.inducedDipoleP, atomJ.labFrameDipole ) +
                    DOT3_4(        atomJ.inducedDipoleP, atomI.labFrameDipole );
    
    float scip2   = DOT3_4(        atomI.inducedDipole,  atomJ.inducedDipoleP) +
                    DOT3_4(        atomJ.inducedDipole,  atomI.inducedDipoleP);
    
    float scip3   = DOT3_4(        atomI.inducedDipoleP, deltaR );
    float scip4   = DOT3_4(        atomJ.inducedDipoleP, deltaR );
    
    float scip7   = DOT3_4(        qIr, atomJ.inducedDipoleP );
    float scip8   = DOT3_4(        qJr, atomI.inducedDipoleP );

    float scaleF             = 0.5f*scalingFactors[UScaleIndex];
    float inducedFactor3     = scip2*rr3*scaleF;
    float inducedFactor5     = (sci3*scip4+scip3*sci4)*rr5*scaleF;
    float findmp_0           = inducedFactor3*ddsc3[0] - inducedFactor5*ddsc5[0];
    float findmp_1           = inducedFactor3*ddsc3[1] - inducedFactor5*ddsc5[1];
    float findmp_2           = inducedFactor3*ddsc3[2] - inducedFactor5*ddsc5[2];

    float gli1               = atomJ.q*sci3 - atomI.q*sci4;
    float gli2               = -sc3*sci4 - sci3*sc4;
    float gli3               = sci3*sc6 - sci4*sc5;
    float gli6               = sci1;
    float gli7               = 2.0f*(sci7-sci8);
    
    float glip1              = atomJ.q*scip3 - atomI.q*scip4;
    float glip2              = -sc3*scip4 - scip3*sc4;
    float glip3              = scip3*sc6 - scip4*sc5;
    float glip6              = scip1;
    float glip7              = 2.0f*(scip7-scip8);
    
    float factor3            = rr3*(( gli1  +  gli6)*scalingFactors[PScaleIndex] + (glip1  + glip6)*scalingFactors[DScaleIndex]);
    float factor5            = rr5*(( gli2  +  gli7)*scalingFactors[PScaleIndex] + (glip2  + glip7)*scalingFactors[DScaleIndex]);
    float factor7            = rr7*( gli3*scalingFactors[PScaleIndex] + glip3*scalingFactors[DScaleIndex]);
      
    float fridmp_0           = 0.5f*(factor3*ddsc3[0] + factor5*ddsc5[0] + factor7*ddsc7[0]);
    float fridmp_1           = 0.5f*(factor3*ddsc3[1] + factor5*ddsc5[1] + factor7*ddsc7[1]);
    float fridmp_2           = 0.5f*(factor3*ddsc3[2] + factor5*ddsc5[2] + factor7*ddsc7[2]);
      
    float gl0 = atomI.q*atomJ.q;
    float gl1 = atomJ.q*sc3 - atomI.q*sc4;
    float gl2 = atomI.q*sc6 + atomJ.q*sc5 - sc3*sc4;
    float gl3 = sc3*sc6 - sc4*sc5;
    float gl4 = sc5*sc6;
    float gl6 = sc2;
    float gl7 = 2.0f*(sc7-sc8);
    float gl8 = 2.0f*sc10;
    float gl5 = -4.0f*sc9;
    
    float gf1 = rr3*gl0 + rr5*(gl1+gl6) + rr7*(gl2+gl7+gl8) + rr9*(gl3+gl5) + rr11*gl4;
    float gf2 = -atomJ.q*rr3 + sc4*rr5 - sc6*rr7;
    float gf3 =  atomI.q*rr3 + sc3*rr5 + sc5*rr7;
    float gf4 = 2.0f*rr5;
    float gf5 = 2.0f*(-atomJ.q*rr5+sc4*rr7-sc6*rr9);
    float gf6 = 2.0f*(-atomI.q*rr5-sc3*rr7-sc5*rr9);
    float gf7 = 4.0f*rr7;

    // energy

    float em                 = scalingFactors[MScaleIndex]*(rr1*gl0 + rr3*(gl1+gl6) + rr5*(gl2+gl7+gl8) + rr7*(gl3+gl5) + rr9*gl4);
    float ei                 = 0.5f*(rr3*(gli1+gli6)*psc0 + rr5*(gli2+gli7)*psc1 + rr7*gli3*psc2);
    outputForce->w           = em+ei;
    
#ifdef AMOEBA_DEBUG
#if 0
if( 1 ){
    int debugIndex           = 0;
    debugArray[debugIndex].x = em;
    debugArray[debugIndex].y = ei;
    debugArray[debugIndex].z = rr1;
    debugArray[debugIndex].w = rr3;

    debugIndex++;
    debugArray[debugIndex].x = gl0;
    debugArray[debugIndex].y = gl1;
    debugArray[debugIndex].z = gl6;
    debugArray[debugIndex].w = gl2;

    debugIndex++;
    debugArray[debugIndex].x = gli1;
    debugArray[debugIndex].y = gli3;
    debugArray[debugIndex].z = gli2;
    debugArray[debugIndex].w = gli7;

    debugIndex++;
    debugArray[debugIndex].x = psc0;
    debugArray[debugIndex].y = psc1;
    debugArray[debugIndex].z = psc2;
    debugArray[debugIndex].w = scalingFactors[MScaleIndex];

}
#endif
#endif

    float temp1[3],temp2[3],temp3[3];
    float qIqJr[3], qJqIr[3], qIdJ[3], qJdI[3];
    amatrixProductVector3( atomI.labFrameQuadrupole,      atomJ.labFrameDipole,     qIdJ );//MK
    amatrixProductVector3( atomJ.labFrameQuadrupole,      atomI.labFrameDipole,     qJdI );//MK

    amatrixProductVector3( atomI.labFrameQuadrupole,      qJr,    qIqJr );//MK
    amatrixProductVector3( atomJ.labFrameQuadrupole,      qIr,    qJqIr );//MK
    amatrixProductVector3( atomJ.labFrameQuadrupole,      qIr,    temp1 );
    amatrixProductVector3( atomJ.labFrameQuadrupole,      atomI.labFrameDipole,     temp2 );

    float ftm2_0 = gf1*deltaR[0] +
                     gf2*atomI.labFrameDipole[0] + gf3*atomJ.labFrameDipole[0]  +
                     gf4*(temp2[0]  - qIdJ[0])   +
                     gf5*qIr[0]    + gf6*qJr[0]  +
                     gf7*(qIqJr[0] + temp1[0]);
    
    float ftm2_1 = gf1*deltaR[1]                 +
                     gf2*atomI.labFrameDipole[1] + gf3*atomJ.labFrameDipole[1]  +
                     gf4*(temp2[1]  - qIdJ[1])   +
                     gf5*qIr[1]    + gf6*qJr[1]  +
                     gf7*(qIqJr[1] + temp1[1]);
    
    float ftm2_2 = gf1*deltaR[2]                 +
                     gf2*atomI.labFrameDipole[2] + gf3*atomJ.labFrameDipole[2]  +
                     gf4*(temp2[2]  - qIdJ[2])   +
                     gf5*qIr[2]    + gf6*qJr[2]  +
                     gf7*(qIqJr[2] + temp1[2]);
    

    // get the induced force;

    // intermediate variables for the induced-permanent terms;
    
    float gfi1 = rr5*0.5f*((gli1+gli6)*psc0 + (glip1+glip6)*dsc0 + scip2*scaleI0) + rr7*((gli7+gli2)*psc1 + (glip7+glip2)*dsc1 -
                                                       (sci3*scip4+scip3*sci4)*scaleI1)*0.5f + rr9*(gli3*psc2+glip3*dsc2)*0.5f;
    float gfi4 = 2.0f*rr5;
    float gfi5 = rr7* (sci4*psc2 + scip4*dsc2);
    float gfi6 = -rr7*(sci3*psc2 + scip3*dsc2);


    float temp4[3];
    float temp5[3];
    float temp6[3];
    float temp7[3];
    float temp8[3];
    float temp9[3];
    float temp10[3];
    float temp11[3];
    float temp12[3];
    float temp13[3];
    float temp14[3];
    float temp15[3];
    float qIuJp[3], qJuIp[3];
    float qIuJ[3], qJuI[3];

    amatrixProductVector3(atomJ.labFrameQuadrupole,      atomI.inducedDipoleP,    temp4);

    amatrixProductVector3(atomI.labFrameQuadrupole,      atomJ.inducedDipoleP,    qIuJp);//MK
    amatrixProductVector3(atomJ.labFrameQuadrupole,      atomI.inducedDipoleP,    qJuIp);//MK
    amatrixProductVector3(atomJ.labFrameQuadrupole,      atomI.inducedDipole ,    qJuI);//MK

    amatrixProductVector3(atomJ.labFrameQuadrupole,      atomI.inducedDipole,    temp5);
    amatrixProductVector3(atomI.labFrameQuadrupole,      atomJ.inducedDipole ,     qIuJ);//MK

/*
               ftm2i(1) = gfi(1)*xr + 0.5d0*
     &           (- rr3*ck*(uind(1,i)*psc3+uinp(1,i)*dsc3)
     &            + rr5*sc(4)*(uind(1,i)*psc5+uinp(1,i)*dsc5)
     &            - rr7*sc(6)*(uind(1,i)*psc7+uinp(1,i)*dsc7))
     &            + (rr3*ci*(uind(1,k)*psc3+uinp(1,k)*dsc3)
     &            + rr5*sc(3)*(uind(1,k)*psc5+uinp(1,k)*dsc5)
     &            + rr7*sc(5)*(uind(1,k)*psc7+uinp(1,k)*dsc7))*0.5d0
     &            + rr5*scale5i*(sci(4)*uinp(1,i)+scip(4)*uind(1,i)
     &            + sci(3)*uinp(1,k)+scip(3)*uind(1,k))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(1)
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(1)
     &            + 0.5d0*gfi(4)*((qkui(1)-qiuk(1))*psc5
     &            + (qkuip(1)-qiukp(1))*dsc5)
     &            + gfi(5)*qir(1) + gfi(6)*qkr(1)
*/

    float ftm2i_0 = gfi1*deltaR[0] +
                    0.5f*(-rr3*atomJ.q*(atomI.inducedDipole[0]*psc0 + atomI.inducedDipoleP[0]*dsc0) +
                    rr5*sc4*(atomI.inducedDipole[0]*psc1 + atomI.inducedDipoleP[0]*dsc1) -
                    rr7*sc6*(atomI.inducedDipole[0]*psc2 + atomI.inducedDipoleP[0]*dsc2)) +
      
                   (rr3*atomI.q*(atomJ.inducedDipole[0]*psc0+atomJ.inducedDipoleP[0]*dsc0) +
                     rr5*sc3*(atomJ.inducedDipole[0]*psc1 +atomJ.inducedDipoleP[0]*dsc1) +
                     rr7*sc5*(atomJ.inducedDipole[0]*psc2 +atomJ.inducedDipoleP[0]*dsc2))*0.5f +
                     rr5*scaleI1*(sci4*atomI.inducedDipoleP[0]+scip4*atomI.inducedDipole[0] +
                     sci3*atomJ.inducedDipoleP[0]+scip3*atomJ.inducedDipole[0])*0.5f +
      
                    0.5f*(sci4*psc1+scip4*dsc1)*rr5*atomI.labFrameDipole[0] +
                    0.5f*(sci3*psc1+scip3*dsc1)*rr5*atomJ.labFrameDipole[0] +
                    0.5f*gfi4*((temp5[0]-qIuJ[0])*psc1 +
                    (temp4[0]-qIuJp[0])*dsc1) + gfi5*qIr[0] + gfi6*qJr[0];
      
    float ftm2i_1  = gfi1*deltaR[1] +
                    0.5f*(-rr3*atomJ.q*(atomI.inducedDipole[1]*psc0 + atomI.inducedDipoleP[1]*dsc0) +
                    rr5*sc4*(atomI.inducedDipole[1]*psc1 + atomI.inducedDipoleP[1]*dsc1) -
                    rr7*sc6*(atomI.inducedDipole[1]*psc2 + atomI.inducedDipoleP[1]*dsc2)) +
      
                    (rr3*atomI.q*(atomJ.inducedDipole[1]*psc0+atomJ.inducedDipoleP[1]*dsc0) +
                     rr5*sc3*(atomJ.inducedDipole[1]*psc1 +atomJ.inducedDipoleP[1]*dsc1) +
                     rr7*sc5*(atomJ.inducedDipole[1]*psc2 +atomJ.inducedDipoleP[1]*dsc2))*0.5f +
                     rr5*scaleI1*(sci4*atomI.inducedDipoleP[1]+scip4*atomI.inducedDipole[1] +
                     sci3*atomJ.inducedDipoleP[1]+scip3*atomJ.inducedDipole[1])*0.5f +
      
                    0.5f*(sci4*psc1+scip4*dsc1)*rr5*atomI.labFrameDipole[1] +
                    0.5f*(sci3*psc1+scip3*dsc1)*rr5*atomJ.labFrameDipole[1] +
                    0.5f*gfi4*((temp5[1]-qIuJ[1])*psc1 +
                    (temp4[1]-qIuJp[1])*dsc1) + gfi5*qIr[1] + gfi6*qJr[1];
      
    float ftm2i_2  = gfi1*deltaR[2] +
                    0.5f*(-rr3*atomJ.q*(atomI.inducedDipole[2]*psc0 + atomI.inducedDipoleP[2]*dsc0) +
                    rr5*sc4*(atomI.inducedDipole[2]*psc1 + atomI.inducedDipoleP[2]*dsc1) -
                    rr7*sc6*(atomI.inducedDipole[2]*psc2 + atomI.inducedDipoleP[2]*dsc2)) +
      
                    (rr3*atomI.q*(atomJ.inducedDipole[2]*psc0+atomJ.inducedDipoleP[2]*dsc0) +
                     rr5*sc3*(atomJ.inducedDipole[2]*psc1 +atomJ.inducedDipoleP[2]*dsc1) +
                     rr7*sc5*(atomJ.inducedDipole[2]*psc2 +atomJ.inducedDipoleP[2]*dsc2))*0.5f +
                     rr5*scaleI1*(sci4*atomI.inducedDipoleP[2]+scip4*atomI.inducedDipole[2] +
                     sci3*atomJ.inducedDipoleP[2]+scip3*atomJ.inducedDipole[2])*0.5f +
      
                    0.5f*(sci4*psc1+scip4*dsc1)*rr5*atomI.labFrameDipole[2] +
                    0.5f*(sci3*psc1+scip3*dsc1)*rr5*atomJ.labFrameDipole[2] +
                    0.5f*gfi4*((temp5[2]-qIuJ[2])*psc1 +
                    (temp4[2]-qIuJp[2])*dsc1) + gfi5*qIr[2] + gfi6*qJr[2];

    // handle of scaling for partially excluded interactions;
    // correction to convert mutual to direct polarization force;
    
    ftm2i_0 -= (fridmp_0 + findmp_0);
    ftm2i_1 -= (fridmp_1 + findmp_1);
    ftm2i_2 -= (fridmp_2 + findmp_2);
    
    // now perform the torque calculation;
    // intermediate terms for torque between multipoles i and j;
    
    float gti2 = 0.5f*(sci4*psc1+scip4*dsc1)*rr5;
    float gti3 = 0.5f*(sci3*psc1+scip3*dsc1)*rr5;
    float gti4 = gfi4;
    float gti5 = gfi5;
    float gti6 = gfi6;

    // get the permanent (ttm2, ttm3) and induced interaction torques (ttm2i, ttm3i)
    
    acrossProductVector3(atomI.labFrameDipole,      atomJ.labFrameDipole,      temp1);
    acrossProductVector3(atomI.labFrameDipole,      atomJ.inducedDipole ,      temp2);
    acrossProductVector3(atomI.labFrameDipole,      atomJ.inducedDipoleP,     temp3);
    acrossProductVector3(atomI.labFrameDipole,      deltaR,       temp4);
    acrossProductVector3(deltaR,       qIuJp,   temp5);
    acrossProductVector3(deltaR,       qIr,     temp6);
    acrossProductVector3(deltaR,       qIuJ,    temp7);
    acrossProductVector3(atomJ.inducedDipole ,     qIr,     temp8);
    acrossProductVector3(atomJ.inducedDipoleP,     qIr,     temp9);
    acrossProductVector3(atomI.labFrameDipole,     qJr,     temp10);
    acrossProductVector3(atomJ.labFrameDipole,     qIr,     temp11);
    acrossProductVector3(deltaR,       qIqJr,   temp12);
    acrossProductVector3(deltaR,       qIdJ,    temp13);

    amatrixCrossProductMatrix3(atomI.labFrameQuadrupole,      atomJ.labFrameQuadrupole,      temp14);
    acrossProductVector3(qJr, qIr,     temp15);

    float ttm2_0  = -rr3*temp1[0] + gf2*temp4[0]-gf5*temp6[0] + gf4*(temp10[0] + temp11[0] + temp13[0]-2.0f*temp14[0]) - gf7*(temp12[0] + temp15[0]);
    float ttm2i_0 = -rr3*(temp2[0]*psc0+temp3[0]*dsc0)*0.5f + gti2*temp4[0] + gti4*((temp8[0]+ temp7[0])*psc1 + (temp9[0] + temp5[0])*dsc1)*0.5f - gti5*temp6[0];
    float ttm2_1  = -rr3*temp1[1] + gf2*temp4[1]-gf5*temp6[1] + gf4*(temp10[1] + temp11[1] + temp13[1]-2.0f*temp14[1]) - gf7*(temp12[1] + temp15[1]);
    float ttm2i_1 = -rr3*(temp2[1]*psc0+temp3[1]*dsc0)*0.5f + gti2*temp4[1] + gti4*((temp8[1]+ temp7[1])*psc1 + (temp9[1] + temp5[1])*dsc1)*0.5f - gti5*temp6[1];
    float ttm2_2  = -rr3*temp1[2] + gf2*temp4[2]-gf5*temp6[2] + gf4*(temp10[2] + temp11[2] + temp13[2]-2.0f*temp14[2]) - gf7*(temp12[2] + temp15[2]);
    float ttm2i_2 = -rr3*(temp2[2]*psc0+temp3[2]*dsc0)*0.5f + gti2*temp4[2] + gti4*((temp8[2]+ temp7[2])*psc1 + (temp9[2] + temp5[2])*dsc1)*0.5f - gti5*temp6[2];

    acrossProductVector3(atomJ.labFrameDipole,      deltaR,       temp2  );
    acrossProductVector3(deltaR,       qJr,     temp3  );
    acrossProductVector3(atomI.labFrameDipole,      qJr,     temp4  );
    acrossProductVector3(atomJ.labFrameDipole,      qIr,     temp5  );
    acrossProductVector3(deltaR,       qJdI,    temp6  );
    acrossProductVector3(deltaR,       qJqIr,   temp7  );
    acrossProductVector3(qJr,     qIr,     temp8  ); // _qJrxqIr
    acrossProductVector3(atomJ.labFrameDipole,      atomI.inducedDipole ,      temp9  ); // _dJxuI
    acrossProductVector3(atomJ.labFrameDipole,      atomI.inducedDipoleP,     temp10 ); // _dJxuIp

    acrossProductVector3(atomI.inducedDipoleP,     qJr,     temp11 ); // _uIxqJrp
    acrossProductVector3(atomI.inducedDipole ,     qJr,     temp12 ); // _uIxqJr
    acrossProductVector3(deltaR,       qJuIp,   temp13 ); // _rxqJuIp
    acrossProductVector3(deltaR,       qJuI,    temp15 ); // _rxqJuI

    float ttm3_0 = rr3*temp1[0] + gf3*temp2[0] - gf6*temp3[0] - gf4*(temp4[0] + temp5[0] + temp6[0] - 2.0f*temp14[0]) - gf7*(temp7[0] - temp8[0]);
    float ttm3i_0 = -rr3*(temp9[0]*psc0+ temp10[0]*dsc0)*0.5f + gti3*temp2[0] - gti4*((temp12[0] + temp15[0])*psc1 + (temp11[0] + temp13[0])*dsc1)*0.5f - gti6*temp3[0];
    float ttm3_1 = rr3*temp1[1] + gf3*temp2[1] - gf6*temp3[1] - gf4*(temp4[1] + temp5[1] + temp6[1] - 2.0f*temp14[1]) - gf7*(temp7[1] - temp8[1]);
    float ttm3i_1 = -rr3*(temp9[1]*psc0+ temp10[1]*dsc0)*0.5f + gti3*temp2[1] - gti4*((temp12[1] + temp15[1])*psc1 + (temp11[1] + temp13[1])*dsc1)*0.5f - gti6*temp3[1];
    float ttm3_2 = rr3*temp1[2] + gf3*temp2[2] - gf6*temp3[2] - gf4*(temp4[2] + temp5[2] + temp6[2] - 2.0f*temp14[2]) - gf7*(temp7[2] - temp8[2]);
    float ttm3i_2 = -rr3*(temp9[2]*psc0+ temp10[2]*dsc0)*0.5f + gti3*temp2[2] - gti4*((temp12[2] + temp15[2])*psc1 + (temp11[2] + temp13[2])*dsc1)*0.5f - gti6*temp3[2];

    if( scalingFactors[MScaleIndex] < 1.0f ){
    
        ftm2_0 *= scalingFactors[MScaleIndex];
        ftm2_1 *= scalingFactors[MScaleIndex];
        ftm2_2 *= scalingFactors[MScaleIndex];
        
        ttm2_0 *= scalingFactors[MScaleIndex];
        ttm2_1 *= scalingFactors[MScaleIndex];
        ttm2_2 *= scalingFactors[MScaleIndex];
        
        ttm3_0 *= scalingFactors[MScaleIndex];
        ttm3_1 *= scalingFactors[MScaleIndex];
        ttm3_2 *= scalingFactors[MScaleIndex];
    
    }


#ifdef AMOEBA_DEBUG
if( 0 ){
int debugIndex               = 0;
    debugArray[debugIndex].x = scalingFactors[DScaleIndex];
    debugArray[debugIndex].y = scalingFactors[PScaleIndex];
    debugArray[debugIndex].z = scalingFactors[MScaleIndex];
    debugArray[debugIndex].w = scalingFactors[UScaleIndex];

    debugIndex++;
    debugArray[debugIndex].x = ftm2i_0 + (fridmp_0 + findmp_0);
    debugArray[debugIndex].y = ftm2i_1 + (fridmp_1 + findmp_1);
    debugArray[debugIndex].z = ftm2i_2 + (fridmp_2 + findmp_2);
    debugArray[debugIndex].w = 1.5;

/*
    debugIndex++;
    debugArray[debugIndex].x = temp2[0];
    debugArray[debugIndex].y = temp2[1];
    debugArray[debugIndex].z = temp2[2];
    debugArray[debugIndex].w = 2.0f;

    debugIndex++;
    debugArray[debugIndex].x = temp3[0];
    debugArray[debugIndex].y = temp3[1];
    debugArray[debugIndex].z = temp3[2];
    debugArray[debugIndex].w = 3.0f;

    debugIndex++;
    debugArray[debugIndex].x = temp4[0];
    debugArray[debugIndex].y = temp4[1];
    debugArray[debugIndex].z = temp4[2];
    debugArray[debugIndex].w = 4.0f;

    debugIndex++;
    debugArray[debugIndex].x = temp5[0];
    debugArray[debugIndex].y = temp5[1];
    debugArray[debugIndex].z = temp5[2];
    debugArray[debugIndex].w = 5.0f;

    debugIndex++;
    debugArray[debugIndex].x = temp6[0];
    debugArray[debugIndex].y = temp6[1];
    debugArray[debugIndex].z = temp6[2];
    debugArray[debugIndex].w = 6.0f;

    debugIndex++;
    debugArray[debugIndex].x = temp14[0];
    debugArray[debugIndex].y = temp14[1];
    debugArray[debugIndex].z = temp14[2];
    debugArray[debugIndex].w = 14.0f;

    debugIndex++;
    debugArray[debugIndex].x = temp7[0];
    debugArray[debugIndex].y = temp7[1];
    debugArray[debugIndex].z = temp7[2];
    debugArray[debugIndex].w = 7.0f;


    debugIndex++;
    debugArray[debugIndex].x = temp8[0];
    debugArray[debugIndex].y = temp8[1];
    debugArray[debugIndex].z = temp8[2];
    debugArray[debugIndex].w = 8.0f;

    debugIndex++;
    debugArray[debugIndex].x = rr3;
    debugArray[debugIndex].y = gf3;
    debugArray[debugIndex].z = gf6;
    debugArray[debugIndex].w = 20.0f;

    debugIndex++;
    debugArray[debugIndex].x = gf4;
    debugArray[debugIndex].y = gf7;
    debugArray[debugIndex].z = 0.0f;
    debugArray[debugIndex].w = 21.0f;

    debugIndex++;
    debugArray[debugIndex].x = atomJ.labFrameDipole[0];
    debugArray[debugIndex].y = atomJ.labFrameDipole[1];
    debugArray[debugIndex].z = atomJ.labFrameDipole[2];
    debugArray[debugIndex].w = 22.0f;

    debugIndex++;
    debugArray[debugIndex].x = deltaR[0];
    debugArray[debugIndex].y = deltaR[1];
    debugArray[debugIndex].z = deltaR[2];
    debugArray[debugIndex].w = 23.0f;
*/

}
#endif

    outputForce->x       = -(ftm2_0 + ftm2i_0);
    outputForce->y       = -(ftm2_1 + ftm2i_1);
    outputForce->z       = -(ftm2_2 + ftm2i_2);
    
    outputTorque[0].x    =  (ttm2_0 + ttm2i_0);
    outputTorque[0].y    =  (ttm2_1 + ttm2i_1);
    outputTorque[0].z    =  (ttm2_2 + ttm2i_2);

    outputTorque[1].x    =  (ttm3_0 + ttm3i_0);
    outputTorque[1].y    =  (ttm3_1 + ttm3i_1);
    outputTorque[1].z    =  (ttm3_2 + ttm3i_2);

    return;

}

__device__ void loadElectrostaticShared( struct ElectrostaticParticle* sA, unsigned int atomI,
                                         float4* atomCoord, float* labFrameDipoleJ, float* labQuadrupole,
                                         float* inducedDipole, float* inducedDipolePolar, float2* dampingFactorAndThole )
{
    // coordinates & charge

    sA->x                        = atomCoord[atomI].x;
    sA->y                        = atomCoord[atomI].y;
    sA->z                        = atomCoord[atomI].z;
    sA->q                        = atomCoord[atomI].w;

    // lab dipole

    sA->labFrameDipole[0]         = labFrameDipoleJ[atomI*3];
    sA->labFrameDipole[1]         = labFrameDipoleJ[atomI*3+1];
    sA->labFrameDipole[2]         = labFrameDipoleJ[atomI*3+2];

    // lab quadrupole

    sA->labFrameQuadrupole[0]    = labQuadrupole[atomI*9];
    sA->labFrameQuadrupole[1]    = labQuadrupole[atomI*9+1];
    sA->labFrameQuadrupole[2]    = labQuadrupole[atomI*9+2];
    sA->labFrameQuadrupole[3]    = labQuadrupole[atomI*9+3];
    sA->labFrameQuadrupole[4]    = labQuadrupole[atomI*9+4];
    sA->labFrameQuadrupole[5]    = labQuadrupole[atomI*9+5];
    sA->labFrameQuadrupole[6]    = labQuadrupole[atomI*9+6];
    sA->labFrameQuadrupole[7]    = labQuadrupole[atomI*9+7];
    sA->labFrameQuadrupole[8]    = labQuadrupole[atomI*9+8];

    // induced dipole

    sA->inducedDipole[0]          = inducedDipole[atomI*3];
    sA->inducedDipole[1]          = inducedDipole[atomI*3+1];
    sA->inducedDipole[2]          = inducedDipole[atomI*3+2];

    // induced dipole polar

    sA->inducedDipoleP[0]         = inducedDipolePolar[atomI*3];
    sA->inducedDipoleP[1]         = inducedDipolePolar[atomI*3+1];
    sA->inducedDipoleP[2]         = inducedDipolePolar[atomI*3+2];

    sA->damp                     = dampingFactorAndThole[atomI].x;
    sA->thole                    = dampingFactorAndThole[atomI].y;

}

// Include versions of the kernels for N^2 calculations.

#undef USE_OUTPUT_BUFFER_PER_WARP
#define METHOD_NAME(a, b) a##N2##b
#include "kCalculateAmoebaCudaElectrostatic.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2ByWarp##b
#include "kCalculateAmoebaCudaElectrostatic.h"

// reduce psWorkArray_3_1 -> force
// reduce psWorkArray_3_2 -> torque

static void kReduceForceTorque(amoebaGpuContext amoebaGpu )
{
    kReduceFields_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->fieldReduceThreadsPerBlock>>>(
                               amoebaGpu->paddedNumberOfAtoms*3, amoebaGpu->outputBuffers,
                               amoebaGpu->psWorkArray_3_1->_pDevData, amoebaGpu->psForce->_pDevData );
    LAUNCHERROR("kReduceElectrostaticForce");
    kReduceFields_kernel<<<amoebaGpu->nonbondBlocks, amoebaGpu->fieldReduceThreadsPerBlock>>>(
                               amoebaGpu->paddedNumberOfAtoms*3, amoebaGpu->outputBuffers,
                               amoebaGpu->psWorkArray_3_2->_pDevData, amoebaGpu->psTorque->_pDevData );
    LAUNCHERROR("kReduceElectrostaticTorque");
}


/**---------------------------------------------------------------------------------------

   Compute Amoeba electrostatic force & torque

   @param amoebaGpu        amoebaGpu context
   @param gpu              OpenMM gpu Cuda context

   --------------------------------------------------------------------------------------- */

void cudaComputeAmoebaElectrostatic( amoebaGpuContext amoebaGpu )
{
  
   // ---------------------------------------------------------------------------------------

    static unsigned int threadsPerBlock = 0;

#ifdef AMOEBA_DEBUG
    static const char* methodName = "cudaComputeAmoebaElectrostatic";
    static int timestep = 0;
    std::vector<int> fileId;
    timestep++;
    fileId.resize( 2 );
    fileId[0] = timestep;
    fileId[1] = 1;
#endif

    // ---------------------------------------------------------------------------------------

    gpuContext gpu = amoebaGpu->gpuContext;

    // apparently debug array can take up nontrivial no. registers

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){
        (void) fprintf( amoebaGpu->log, "%s %d maxCovalentDegreeSz=%d ZZZ\n",
                        methodName, gpu->natoms, amoebaGpu->maxCovalentDegreeSz );
    }   
    static const int maxSlots                 =20;
    int paddedNumberOfAtoms                   = amoebaGpu->gpuContext->sim.paddedNumberOfAtoms;
    CUDAStream<float4>* debugArray            = new CUDAStream<float4>(maxSlots*paddedNumberOfAtoms, 1, "DebugArray");
    memset( debugArray->_pSysData,      0, sizeof( float )*4*maxSlots*paddedNumberOfAtoms);
    debugArray->Upload();
    unsigned int targetAtom                   = 237;
#endif

    // on first pass, set threads/block

    if( threadsPerBlock == 0 ){
        unsigned int maxThreads;
        if (gpu->sm_version >= SM_20)
            maxThreads = 384;
        else if (gpu->sm_version >= SM_12)
            maxThreads = 128;
        else
            maxThreads = 64;
        threadsPerBlock = std::min(getThreadsPerBlock(amoebaGpu, sizeof(ElectrostaticParticle)), maxThreads);
    }

    kClearFields_3( amoebaGpu, 2 );
    LAUNCHERROR("kClearFields_3 kCalculateAmoebaCudaElectrostatic");

    if (gpu->bOutputBufferPerWarp){

#ifdef AMOEBA_DEBUG
        (void) fprintf( amoebaGpu->log, "kCalculateAmoebaCudaElectrostaticN2Forces warp:  numBlocks=%u numThreads=%u bufferPerWarp=%u atm=%lu shrd=%lu ixnCt=%lu workUnits=%u\n",
                        amoebaGpu->nonbondBlocks, threadsPerBlock, amoebaGpu->bOutputBufferPerWarp,
                        sizeof(ElectrostaticParticle), sizeof(ElectrostaticParticle)*threadsPerBlock, (*gpu->psInteractionCount)[0], gpu->sim.workUnits ); (void) fflush( amoebaGpu->log );
#endif


        kCalculateAmoebaCudaElectrostaticN2ByWarpForces_kernel<<<amoebaGpu->nonbondBlocks, threadsPerBlock, sizeof(ElectrostaticParticle)*threadsPerBlock>>>(
                                                                           amoebaGpu->psWorkUnit->_pDevData,
                                                                           gpu->psPosq4->_pDevData,
                                                                           amoebaGpu->psLabFrameDipole->_pDevData,
                                                                           amoebaGpu->psLabFrameQuadrupole->_pDevData,
                                                                           amoebaGpu->psInducedDipole->_pDevData,
                                                                           amoebaGpu->psInducedDipolePolar->_pDevData,
                                                                           amoebaGpu->psWorkArray_3_1->_pDevData,
#ifdef AMOEBA_DEBUG
                                                                           amoebaGpu->psWorkArray_3_2->_pDevData,
                                                                           debugArray->_pDevData, targetAtom );
#else
                                                                           amoebaGpu->psWorkArray_3_2->_pDevData );
#endif

    } else {

#ifdef AMOEBA_DEBUG
        (void) fprintf( amoebaGpu->log, "kCalculateAmoebaCudaElectrostaticN2Forces no warp:  numBlocks=%u numThreads=%u bufferPerWarp=%u atm=%u shrd=%u xnCt=%u workUnits=%u\n",
                        amoebaGpu->nonbondBlocks, threadsPerBlock, amoebaGpu->bOutputBufferPerWarp,
                        sizeof(ElectrostaticParticle), sizeof(ElectrostaticParticle)*threadsPerBlock, (*gpu->psInteractionCount)[0], gpu->sim.workUnits );
        (void) fflush( amoebaGpu->log );
#endif

        kCalculateAmoebaCudaElectrostaticN2Forces_kernel<<<amoebaGpu->nonbondBlocks, threadsPerBlock, sizeof(ElectrostaticParticle)*threadsPerBlock>>>(
                                                                           amoebaGpu->psWorkUnit->_pDevData,
                                                                           gpu->psPosq4->_pDevData,
                                                                           amoebaGpu->psLabFrameDipole->_pDevData,
                                                                           amoebaGpu->psLabFrameQuadrupole->_pDevData,
                                                                           amoebaGpu->psInducedDipole->_pDevData,
                                                                           amoebaGpu->psInducedDipolePolar->_pDevData,
                                                                           amoebaGpu->psWorkArray_3_1->_pDevData,
#ifdef AMOEBA_DEBUG
                                                                           amoebaGpu->psWorkArray_3_2->_pDevData,
                                                                           debugArray->_pDevData, targetAtom );
#else
                                                                           amoebaGpu->psWorkArray_3_2->_pDevData );
#endif
    }
    LAUNCHERROR("kCalculateAmoebaCudaElectrostaticN2Forces");

    kReduceForceTorque( amoebaGpu );
    LAUNCHERROR("kReduceForceTorque");

#ifdef AMOEBA_DEBUG
    if( amoebaGpu->log ){

        amoebaGpu->psForce->Download();
        amoebaGpu->psTorque->Download();
        debugArray->Download();

        (void) fprintf( amoebaGpu->log, "Finished Electrostatic kernel execution\n" ); (void) fflush( amoebaGpu->log );

        int maxPrint        = 1400;
        for( int ii = 0; ii < gpu->natoms; ii++ ){
           (void) fprintf( amoebaGpu->log, "%5d ", ii); 

            int indexOffset     = ii*3;
    
           // force

           (void) fprintf( amoebaGpu->log,"ElectrostaticF [%16.9e %16.9e %16.9e] ",
                           amoebaGpu->psForce->_pSysData[indexOffset],
                           amoebaGpu->psForce->_pSysData[indexOffset+1],
                           amoebaGpu->psForce->_pSysData[indexOffset+2] );
    
           // torque

           (void) fprintf( amoebaGpu->log,"ElectrostaticT [%16.9e %16.9e %16.9e] ",
                           amoebaGpu->psTorque->_pSysData[indexOffset],
                           amoebaGpu->psTorque->_pSysData[indexOffset+1],
                           amoebaGpu->psTorque->_pSysData[indexOffset+2] );

           (void) fprintf( amoebaGpu->log,"\n" );
           if( ii == maxPrint && (gpu->natoms - maxPrint) > ii ){
                ii = gpu->natoms - maxPrint;
           }
        }
        if( 1 ){
            (void) fprintf( amoebaGpu->log,"DebugElec\n" );
            int paddedNumberOfAtoms = amoebaGpu->gpuContext->sim.paddedNumberOfAtoms;
            for( int jj = 0; jj < gpu->natoms; jj++ ){
                int debugIndex = jj;
                for( int kk = 0; kk < 8; kk++ ){
                    float conversion = kk >= 1 && kk <= 8 ? 1.0f/4.184f : 1.0;
                    (void) fprintf( amoebaGpu->log,"%5d %5d [%16.9e %16.9e %16.9e %16.9e] E11\n", targetAtom, jj,
                                    conversion*debugArray->_pSysData[debugIndex].x, conversion*debugArray->_pSysData[debugIndex].y,
                                    conversion*debugArray->_pSysData[debugIndex].z, debugArray->_pSysData[debugIndex].w );
                    debugIndex += paddedNumberOfAtoms;
                }
                (void) fprintf( amoebaGpu->log,"\n" );
            }
        }
        if( 1 ){
            (void) fprintf( amoebaGpu->log,"DebugElec\n" );
            int paddedNumberOfAtoms = amoebaGpu->gpuContext->sim.paddedNumberOfAtoms;
            for( int jj = 0; jj < gpu->natoms; jj++ ){
                int debugIndex1 = jj + paddedNumberOfAtoms;
                int debugIndex2 = jj + 5*paddedNumberOfAtoms;
                int debugIndex3 = jj + 6*paddedNumberOfAtoms;
                int debugIndex4 = jj + 4*paddedNumberOfAtoms;
                int debugIndex5 = jj + 7*paddedNumberOfAtoms;
                float conversion = 1.0f/4.184f;
                int i1,i2;
                if( jj < targetAtom ){
                    i1 = jj;
                    i2 = targetAtom;
                } else {
                    i1 = targetAtom;
                    i2 = jj;
                }
                (void) fprintf( amoebaGpu->log,"%5d %5d %16.9e %16.9e %16.9e    %16.9e %16.9e %16.9e   %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e F11\n", i1,i2,
                                conversion*debugArray->_pSysData[debugIndex1].x,
                                conversion*debugArray->_pSysData[debugIndex1].y,
                                conversion*debugArray->_pSysData[debugIndex1].z, 
                                conversion*debugArray->_pSysData[debugIndex2].x,
                                conversion*debugArray->_pSysData[debugIndex2].y,
                                conversion*debugArray->_pSysData[debugIndex2].z, 
                                conversion*debugArray->_pSysData[debugIndex3].x,
                                conversion*debugArray->_pSysData[debugIndex3].y,
                                conversion*debugArray->_pSysData[debugIndex3].z,
                                conversion*debugArray->_pSysData[debugIndex5].x,
                                conversion*debugArray->_pSysData[debugIndex5].y,
                                conversion*debugArray->_pSysData[debugIndex5].z );
            }
        }
        (void) fflush( amoebaGpu->log );

        if( 0 ){
            (void) fprintf( amoebaGpu->log, "%s Tiled F & T\n", methodName ); fflush( amoebaGpu->log );
            int maxPrint = 12;
            for( int ii = 0; ii < gpu->natoms; ii++ ){
    
                // print cpu & gpu reductions
    
                int offset  = 3*ii;
    
                (void) fprintf( amoebaGpu->log,"%6d F[%16.7e %16.7e %16.7e] T[%16.7e %16.7e %16.7e]\n", ii,
                                amoebaGpu->psForce->_pSysData[offset],
                                amoebaGpu->psForce->_pSysData[offset+1],
                                amoebaGpu->psForce->_pSysData[offset+2],
                                amoebaGpu->psTorque->_pSysData[offset],
                                amoebaGpu->psTorque->_pSysData[offset+1],
                                amoebaGpu->psTorque->_pSysData[offset+2] );
                if( (ii == maxPrint) && (ii < (gpu->natoms - maxPrint)) )ii = gpu->natoms - maxPrint; 
            }   
        }   

        if( 1 ){
            std::vector<int> fileId;
            //fileId.push_back( 0 );
            VectorOfDoubleVectors outputVector;
            cudaLoadCudaFloat4Array( gpu->natoms, 3, gpu->psPosq4,            outputVector, NULL, 1.0f );
            cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psForce,      outputVector, NULL, 1.0f );
            cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psTorque,     outputVector, NULL, 1.0f );
            cudaWriteVectorOfDoubleVectorsToFile( "CudaForceTorque", fileId, outputVector );
         }

    }   
    delete debugArray;

#endif

    if( 0 ){
        std::vector<int> fileId;
        //fileId.push_back( 0 );
        VectorOfDoubleVectors outputVector;
        //cudaLoadCudaFloat4Array( gpu->natoms, 3, gpu->psPosq4,            outputVector, NULL, 1.0f );
        cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psForce,      outputVector, NULL, 1.0f/4.184 );
        cudaLoadCudaFloatArray( gpu->natoms,  3, amoebaGpu->psTorque,     outputVector, NULL, 1.0f/4.184 );
        cudaWriteVectorOfDoubleVectorsToFile( "CudaForceTorque", fileId, outputVector );
     }

   // ---------------------------------------------------------------------------------------
}
