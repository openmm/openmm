/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Scott Le Grand, Peter Eastman                                     *
 * Contributors:                                                              *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 * -------------------------------------------------------------------------- */

#include "amoebaGpuTypes.h"
#include "amoebaCudaKernels.h"
#include "kCalculateAmoebaCudaUtilities.h"
#include <stdio.h>

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
static int const LastScalingIndex       =  4;

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

    //float torque[3];
    //float padding;

};

#ifdef Original

#define i35 0.257142857f
#define DOT3_4(u,v) ((u[0])*(v[0]) + (u[1])*(v[1]) + (u[2])*(v[2]))

#define MATRIXDOT31(u,v) u[0]*v[0] + u[1]*v[1] + u[2]*v[2] + \
  u[3]*v[3] + u[4]*v[4] + u[5]*v[5] + \
  u[6]*v[6] + u[7]*v[7] + u[8]*v[8]

#define DOT31(u,v) ((u[0])*(v[0]) + (u[1])*(v[1]) + (u[2])*(v[2]))

#define one 1.0f

__device__ void calculateElectrostaticPairIxnOrig_kernel( ElectrostaticParticle& atomI,   ElectrostaticParticle& atomJ,
                                                      float* scalingFactors, float4*  outputForce, float4  outputTorque[2]){
  
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

        float damp                    = ratio*ratio*ratio*pGamma;
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
    
    if( cAmoebaSim.polarizationType )
    {
        float gfd     = 0.5*(rr5*scip2*scaleI0 - rr7*(scip3*sci4+sci3*scip4)*scaleI1);
        float temp5   = 0.5*rr5*scaleI1;
        float fdir_0  = gfd*deltaR[0] + temp5*(sci4*atomI.inducedDipoleP[0] + scip4*atomI.inducedDipole[0] + sci3*atomJ.inducedDipoleP[0] + scip3*atomJ.inducedDipole[0]);
        float fdir_1  = gfd*deltaR[1] + temp5*(sci4*atomI.inducedDipoleP[1] + scip4*atomI.inducedDipole[1] + sci3*atomJ.inducedDipoleP[1] + scip3*atomJ.inducedDipole[1]);
        float fdir_2  = gfd*deltaR[2] + temp5*(sci4*atomI.inducedDipoleP[2] + scip4*atomI.inducedDipole[2] + sci3*atomJ.inducedDipoleP[2] + scip3*atomJ.inducedDipole[2]);
        ftm2i_0      -= fdir_0 - findmp_0;
        ftm2i_1      -= fdir_1 - findmp_1;
        ftm2i_2      -= fdir_2 - findmp_2;

    }
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
#endif

static __device__ void loadElectrostaticParticle( volatile struct ElectrostaticParticle* sA, unsigned int atomI ){

    // coordinates & charge

    sA->x                        = cSim.pPosq[atomI].x;
    sA->y                        = cSim.pPosq[atomI].y;
    sA->z                        = cSim.pPosq[atomI].z;
    sA->q                        = cSim.pPosq[atomI].w;

    // lab dipole

    sA->labFrameDipole[0]        = cAmoebaSim.pLabFrameDipole[atomI*3];
    sA->labFrameDipole[1]        = cAmoebaSim.pLabFrameDipole[atomI*3+1];
    sA->labFrameDipole[2]        = cAmoebaSim.pLabFrameDipole[atomI*3+2];

    // lab quadrupole

    sA->labFrameQuadrupole[0]    = cAmoebaSim.pLabFrameQuadrupole[atomI*9];
    sA->labFrameQuadrupole[1]    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+1];
    sA->labFrameQuadrupole[2]    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+2];
    sA->labFrameQuadrupole[3]    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+3];
    sA->labFrameQuadrupole[4]    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+4];
    sA->labFrameQuadrupole[5]    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+5];
    sA->labFrameQuadrupole[6]    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+6];
    sA->labFrameQuadrupole[7]    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+7];
    sA->labFrameQuadrupole[8]    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+8];

    // induced dipole

    sA->inducedDipole[0]         = cAmoebaSim.pInducedDipole[atomI*3];
    sA->inducedDipole[1]         = cAmoebaSim.pInducedDipole[atomI*3+1];
    sA->inducedDipole[2]         = cAmoebaSim.pInducedDipole[atomI*3+2];

    // induced dipole polar

    sA->inducedDipoleP[0]        = cAmoebaSim.pInducedDipolePolar[atomI*3];
    sA->inducedDipoleP[1]        = cAmoebaSim.pInducedDipolePolar[atomI*3+1];
    sA->inducedDipoleP[2]        = cAmoebaSim.pInducedDipolePolar[atomI*3+2];

    sA->damp                     = cAmoebaSim.pDampingFactorAndThole[atomI].x;
    sA->thole                    = cAmoebaSim.pDampingFactorAndThole[atomI].y;

}

static __device__ void zeroElectrostaticParticle( volatile struct ElectrostaticParticle* sA ){
    sA->force[0]                 = 0.0f;
    sA->force[1]                 = 0.0f;
    sA->force[2]                 = 0.0f;
}

#undef SUB_METHOD_NAME
#undef F1
#define SUB_METHOD_NAME(a, b) a##F1##b
#define F1
#include "kCalculateAmoebaCudaElectrostatic_b.h"
#undef F1
#undef SUB_METHOD_NAME

#undef SUB_METHOD_NAME
#undef F2
#define SUB_METHOD_NAME(a, b) a##F2##b
#define F2
//#include "kCalculateAmoebaCudaElectrostatic_b.h"
#undef F2
#undef SUB_METHOD_NAME

#undef SUB_METHOD_NAME
#undef T1
#define SUB_METHOD_NAME(a, b) a##T1##b
#define T1
#include "kCalculateAmoebaCudaElectrostatic_b.h"
#undef T1
#undef SUB_METHOD_NAME

#undef SUB_METHOD_NAME
#undef T3
#define SUB_METHOD_NAME(a, b) a##T3##b
#define T3
#include "kCalculateAmoebaCudaElectrostatic_b.h"
#undef T3
#undef SUB_METHOD_NAME

__device__ void calculateElectrostaticPairIxn_kernel( ElectrostaticParticle& atomI,   ElectrostaticParticle& atomJ,
                                                      float* scalingFactors, float4*  outputForce, float4 outputTorque[2], float forceFactor){
#ifdef Orig
    return calculateElectrostaticPairIxn_kernel( atomI, atomJ, scalingFactors, outputForce, outputTorque);
#else

    float force[3];
    float energy;
    calculateElectrostaticPairIxnF1_kernel( atomI,  atomJ, scalingFactors, &energy, force);
    outputForce->x = force[0];
    outputForce->y = force[1];
    outputForce->z = force[2];
    outputForce->w = energy;

    calculateElectrostaticPairIxnT1_kernel( atomI,  atomJ, scalingFactors, force);
    outputTorque[0].x = force[0];
    outputTorque[0].y = force[1];
    outputTorque[0].z = force[2];

    calculateElectrostaticPairIxnT3_kernel( atomI,  atomJ, scalingFactors, force);
    outputTorque[1].x = force[0];
    outputTorque[1].y = force[1];
    outputTorque[1].z = force[2];

    return;
#endif

}

// Include versions of the kernels for N^2 calculations.

#undef USE_OUTPUT_BUFFER_PER_WARP
#define METHOD_NAME(a, b) a##N2##b
#include "kCalculateAmoebaCudaElectrostatic.h"
#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##N2ByWarp##b
#include "kCalculateAmoebaCudaElectrostatic.h"

// reduce psWorkArray_3_1 -> torque

static void kReduceTorque(amoebaGpuContext amoebaGpu ){
    gpuContext gpu = amoebaGpu->gpuContext;
    kReduceFields_kernel<<<gpu->sim.nonbond_blocks, gpu->sim.bsf_reduce_threads_per_block>>>(
                               gpu->sim.paddedNumberOfAtoms*3, gpu->sim.outputBuffers,
                               amoebaGpu->psWorkArray_3_1->_pDevData, amoebaGpu->psTorque->_pDevData, 0 );
    LAUNCHERROR("kReduceElectrostaticTorque");
}

/**---------------------------------------------------------------------------------------

   Compute Amoeba electrostatic force & torque

   @param amoebaGpu        amoebaGpu context
   @param addTorqueToForce if set, then add force resulting from torque to force array

   --------------------------------------------------------------------------------------- */

void cudaComputeAmoebaElectrostatic( amoebaGpuContext amoebaGpu, int addTorqueToForce ){
  
   // ---------------------------------------------------------------------------------------

    gpuContext gpu = amoebaGpu->gpuContext;

    // on first pass, set threads/block

    static unsigned int threadsPerBlock = 0;
    if( threadsPerBlock == 0 ){
        unsigned int maxThreads;
        if (gpu->sm_version >= SM_20)
            //maxThreads = 384;
            maxThreads = 512;
        else if (gpu->sm_version >= SM_12)
            maxThreads = 128;
        else
            maxThreads = 64;
        threadsPerBlock = std::min(getThreadsPerBlock(amoebaGpu, sizeof(ElectrostaticParticle), gpu->sharedMemoryPerBlock), maxThreads);
    }

    kClearFields_3( amoebaGpu, 1 );

    if (gpu->bOutputBufferPerWarp){
        kCalculateAmoebaCudaElectrostaticN2ByWarpForces_kernel<<<gpu->sim.nonbond_blocks, threadsPerBlock, sizeof(ElectrostaticParticle)*threadsPerBlock>>>(
                                                                           gpu->psWorkUnit->_pDevData, amoebaGpu->psWorkArray_3_1->_pDevData );
    } else {
        kCalculateAmoebaCudaElectrostaticN2Forces_kernel<<<gpu->sim.nonbond_blocks, threadsPerBlock, sizeof(ElectrostaticParticle)*threadsPerBlock>>>(
                                                                           gpu->psWorkUnit->_pDevData, amoebaGpu->psWorkArray_3_1->_pDevData );
    }
    LAUNCHERROR("kCalculateAmoebaCudaElectrostaticN2Forces");

    if( addTorqueToForce ){
        kReduceTorque( amoebaGpu );
        cudaComputeAmoebaMapTorqueAndAddToForce( amoebaGpu, amoebaGpu->psTorque );
    }

   // ---------------------------------------------------------------------------------------
}

struct ElectrostaticPotentialParticle {

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

};

/**---------------------------------------------------------------------------------------

   Load data for particle w/ index=atomI

   @param sa        address to store atomI's coordinates and multipole moments
   @param atomI     index of atom whose data is to be stored

   --------------------------------------------------------------------------------------- */

static __device__ void loadElectrostaticPotentialParticle( volatile struct ElectrostaticPotentialParticle* sA, unsigned int atomI ){

    // coordinates & charge

    sA->x                        = cSim.pPosq[atomI].x;
    sA->y                        = cSim.pPosq[atomI].y;
    sA->z                        = cSim.pPosq[atomI].z;
    sA->q                        = cSim.pPosq[atomI].w;

    // lab dipole

    sA->labFrameDipole[0]        = cAmoebaSim.pLabFrameDipole[atomI*3];
    sA->labFrameDipole[1]        = cAmoebaSim.pLabFrameDipole[atomI*3+1];
    sA->labFrameDipole[2]        = cAmoebaSim.pLabFrameDipole[atomI*3+2];

    // lab quadrupole

    sA->labFrameQuadrupole[0]    = cAmoebaSim.pLabFrameQuadrupole[atomI*9];
    sA->labFrameQuadrupole[1]    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+1];
    sA->labFrameQuadrupole[2]    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+2];
    sA->labFrameQuadrupole[3]    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+3];
    sA->labFrameQuadrupole[4]    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+4];
    sA->labFrameQuadrupole[5]    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+5];
    sA->labFrameQuadrupole[6]    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+6];
    sA->labFrameQuadrupole[7]    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+7];
    sA->labFrameQuadrupole[8]    = cAmoebaSim.pLabFrameQuadrupole[atomI*9+8];

    // induced dipole

    sA->inducedDipole[0]         = cAmoebaSim.pInducedDipole[atomI*3];
    sA->inducedDipole[1]         = cAmoebaSim.pInducedDipole[atomI*3+1];
    sA->inducedDipole[2]         = cAmoebaSim.pInducedDipole[atomI*3+2];

}

/**---------------------------------------------------------------------------------------

   Calculate potential at grid point due atomI
   Code adapted from TINKER routine potpoint in potpoint.f

   @param atomI     atomI's coordinates and multipole moments
   @param gridPoint grid coordinates
   @param potential output potential

   --------------------------------------------------------------------------------------- */

__device__ void calculateElectrostaticPotentialForAtomGridPoint_kernel( volatile ElectrostaticPotentialParticle& atomI, volatile float4& gridPoint, float* potential ){
  
    float xr                 = atomI.x - gridPoint.x;
    float yr                 = atomI.y - gridPoint.y;
    float zr                 = atomI.z - gridPoint.z;
   
    float r2                 = xr*xr + yr*yr + zr*zr;
    float r                  = sqrtf( r2 );

    float rr1                = 1.0f/r;
    *potential               = atomI.q*rr1;
    float rr2                = rr1*rr1;
    float rr3                = rr1*rr2;

    float scd                = atomI.labFrameDipole[0]*xr     +  atomI.labFrameDipole[1]*yr    + atomI.labFrameDipole[2]*zr;
    float scu                =  atomI.inducedDipole[0]*xr     +   atomI.inducedDipole[1]*yr    +  atomI.inducedDipole[2]*zr;
    *potential              -= (scd + scu)*rr3;

    float rr5                = 3.0f*rr3*rr2;
    float scq                = xr*(atomI.labFrameQuadrupole[0]*xr + atomI.labFrameQuadrupole[1]*yr + atomI.labFrameQuadrupole[2]*zr);
          scq               += yr*(atomI.labFrameQuadrupole[1]*xr + atomI.labFrameQuadrupole[4]*yr + atomI.labFrameQuadrupole[5]*zr);
          scq               += zr*(atomI.labFrameQuadrupole[2]*xr + atomI.labFrameQuadrupole[5]*yr + atomI.labFrameQuadrupole[8]*zr);
    *potential              += scq*rr5;

    return;

}

// Include versions of the kernels for N x PotentialGridSize calculations.

#undef USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME 
#define METHOD_NAME(a, b) a##NxG##b
#include "kCalculateAmoebaCudaElectrostaticPotential.h"

#define USE_OUTPUT_BUFFER_PER_WARP
#undef METHOD_NAME
#define METHOD_NAME(a, b) a##NxGByWarp##b
#include "kCalculateAmoebaCudaElectrostaticPotential.h"

// Kernel to reduce potential

__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 120)
__launch_bounds__(GT2XX_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_THREADS_PER_BLOCK, 1)
#endif
void kReducePotential_kernel()
{
    unsigned int pos             = (blockIdx.x * blockDim.x + threadIdx.x);
    float conversionFactor       = (cAmoebaSim.electric/cAmoebaSim.dielec);
   
    // Reduce potential
    while (pos < cAmoebaSim.paddedPotentialGridSize)
    {
        float totalPotential         = 0.0f;
        float* pFt                   = cAmoebaSim.pPotential + pos;
        int i                        = cSim.outputBuffers;
        while (i >= 4)
        {
            float f1             = *pFt;
            pFt                 += cAmoebaSim.paddedPotentialGridSize;
            float f2             = *pFt;
            pFt                 += cAmoebaSim.paddedPotentialGridSize;
            float f3             = *pFt;
            pFt                 += cAmoebaSim.paddedPotentialGridSize;
            float f4             = *pFt;
            pFt                 += cAmoebaSim.paddedPotentialGridSize;
            totalPotential      += f1 + f2 + f3 + f4;
            i                   -= 4;
        }
        if (i >= 2)
        {
            float f1             = *pFt;
            pFt                 += cAmoebaSim.paddedPotentialGridSize;
            float f2             = *pFt;
            pFt                 += cAmoebaSim.paddedPotentialGridSize;
            totalPotential      += f1 + f2;
            i                   -= 2;
        }
        if (i > 0)
        {
            totalPotential += *pFt;
        }
        totalPotential *= conversionFactor;
        pFt             = cAmoebaSim.pPotential + pos;
        *pFt            = totalPotential;
        pos            += gridDim.x*blockDim.x;
    }   
}

/**---------------------------------------------------------------------------------------

   Reduce Amoeba electrostatic potential

   @param gpu        gpu context

   --------------------------------------------------------------------------------------- */

void kReducePotential(gpuContext gpu)
{
    kReducePotential_kernel<<<gpu->sim.blocks, gpu->sim.bsf_reduce_threads_per_block>>>();
    LAUNCHERROR("kReducePotential");
}

/**---------------------------------------------------------------------------------------

   Compute Amoeba electrostatic potential

   @param amoebaGpu        amoebaGpu context

   --------------------------------------------------------------------------------------- */

void cudaComputeAmoebaElectrostaticPotential( amoebaGpuContext amoebaGpu ){
  
   // ---------------------------------------------------------------------------------------

    gpuContext gpu = amoebaGpu->gpuContext;

    // on first pass, set threads/block

    static unsigned int threadsPerBlock = 0;
    if( threadsPerBlock == 0 ){
        unsigned int maxThreads;
        if (gpu->sm_version >= SM_20)
            //maxThreads = 384;
            maxThreads = 512;
        else if (gpu->sm_version >= SM_12)
            maxThreads = 128;
        else
            maxThreads = 64;
        threadsPerBlock = std::min(getThreadsPerBlock(amoebaGpu, sizeof(ElectrostaticPotentialParticle), gpu->sharedMemoryPerBlock), maxThreads);
    }

    if (gpu->bOutputBufferPerWarp){
        kCalculateAmoebaCudaElectrostaticPotentialNxGByWarp_kernel<<<gpu->sim.nonbond_blocks, threadsPerBlock, sizeof(ElectrostaticPotentialParticle)*threadsPerBlock>>>( );
    } else {
        kCalculateAmoebaCudaElectrostaticPotentialNxG_kernel<<<gpu->sim.nonbond_blocks, threadsPerBlock, sizeof(ElectrostaticPotentialParticle)*threadsPerBlock>>>( );
    }
    LAUNCHERROR("kCalculateAmoebaCudaElectrostaticPotential");

    kReducePotential( amoebaGpu->gpuContext );

   // ---------------------------------------------------------------------------------------
}
