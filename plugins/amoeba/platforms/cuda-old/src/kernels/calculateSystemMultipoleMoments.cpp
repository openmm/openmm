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

#include "amoebaCudaKernels.h"
#include "openmm/OpenMMException.h"

#include <stdio.h>
using namespace std; 

void kCalculateAmoebaSystemMultipoleMoments( amoebaGpuContext amoebaGpu, std::vector< double >& outputMultipoleMoments ) 
{

    // setup

    kSetupAmoebaMultipoleForces(amoebaGpu, false ); 

    gpuContext gpu         = amoebaGpu->gpuContext;

    gpu->psPosq4->Download();
    gpu->psVelm4->Download();
    float4* posq           = gpu->psPosq4->_pSysData;    
    float4* velm           = gpu->psVelm4->_pSysData;    
    float totalMass        = 0.0f;
    float centerOfMass[3]  = { 0.0f, 0.0f, 0.0f };
    for( unsigned int ii  = 0; ii < gpu->natoms; ii++ ){
        float mass;
        if( velm->w > 0.0f ){
            mass        = 1.0f/velm[ii].w;
        } else {
            mass        = 0.0f;
        }
        totalMass        += mass;
        centerOfMass[0]  += mass*posq[ii].x;
        centerOfMass[1]  += mass*posq[ii].y;
        centerOfMass[2]  += mass*posq[ii].z;
    }

    std::vector<float4>  posqLocal(gpu->natoms);
    if( totalMass > 0.0f ){
        centerOfMass[0]  /= totalMass;
        centerOfMass[1]  /= totalMass;
        centerOfMass[2]  /= totalMass;
    }
    for( unsigned int ii  = 0; ii < gpu->natoms; ii++ ){
        posqLocal[ii].x = posq[ii].x  - centerOfMass[0];
        posqLocal[ii].y = posq[ii].y  - centerOfMass[1];
        posqLocal[ii].z = posq[ii].z  - centerOfMass[2];
        posqLocal[ii].w = posq[ii].w;
    }

    float netchg  = 0.0f;

    float xdpl    = 0.0f;
    float ydpl    = 0.0f;
    float zdpl    = 0.0f;

    float xxqdp   = 0.0f;
    float xyqdp   = 0.0f;
    float xzqdp   = 0.0f;

    float yxqdp   = 0.0f;
    float yyqdp   = 0.0f;
    float yzqdp   = 0.0f;

    float zxqdp   = 0.0f;
    float zyqdp   = 0.0f;
    float zzqdp   = 0.0f;

    amoebaGpu->psLabFrameDipole->Download();
    float* labFrameDipole      = amoebaGpu->psLabFrameDipole->_pSysData;    

    amoebaGpu->psInducedDipole->Download();
    float* inducedDipole       = amoebaGpu->psInducedDipole->_pSysData;    

    amoebaGpu->psLabFrameQuadrupole->Download();
    float* labFrameQuadrupole  = amoebaGpu->psLabFrameQuadrupole->_pSysData;    
    for( unsigned int ii  = 0; ii < gpu->natoms; ii++ ){

        netchg              += posqLocal[ii].w;

        float netDipoleX     = (labFrameDipole[3*ii]    + inducedDipole[3*ii]);
        float netDipoleY     = (labFrameDipole[3*ii+1]  + inducedDipole[3*ii+1]);
        float netDipoleZ     = (labFrameDipole[3*ii+2]  + inducedDipole[3*ii+2]);

        xdpl    += posqLocal[ii].x*posqLocal[ii].w + netDipoleX;
        ydpl    += posqLocal[ii].y*posqLocal[ii].w + netDipoleY;
        zdpl    += posqLocal[ii].z*posqLocal[ii].w + netDipoleZ;

        xxqdp   += posqLocal[ii].x*posqLocal[ii].x*posqLocal[ii].w + 2.0f*posqLocal[ii].x*netDipoleX;
        xyqdp   += posqLocal[ii].x*posqLocal[ii].y*posqLocal[ii].w + posqLocal[ii].x*netDipoleY + posqLocal[ii].y*netDipoleX;
        xzqdp   += posqLocal[ii].x*posqLocal[ii].z*posqLocal[ii].w + posqLocal[ii].x*netDipoleZ + posqLocal[ii].z*netDipoleX;

        yxqdp   += posqLocal[ii].y*posqLocal[ii].x*posqLocal[ii].w + posqLocal[ii].y*netDipoleX + posqLocal[ii].x*netDipoleY;
        yyqdp   += posqLocal[ii].y*posqLocal[ii].y*posqLocal[ii].w + 2.0f*posqLocal[ii].y*netDipoleY;
        yzqdp   += posqLocal[ii].y*posqLocal[ii].z*posqLocal[ii].w + posqLocal[ii].y*netDipoleZ + posqLocal[ii].z*netDipoleY;

        zxqdp   += posqLocal[ii].z*posqLocal[ii].x*posqLocal[ii].w + posqLocal[ii].z*netDipoleX + posqLocal[ii].x*netDipoleZ;
        zyqdp   += posqLocal[ii].z*posqLocal[ii].y*posqLocal[ii].w + posqLocal[ii].z*netDipoleY + posqLocal[ii].y*netDipoleZ;
        zzqdp   += posqLocal[ii].z*posqLocal[ii].z*posqLocal[ii].w + 2.0f*posqLocal[ii].z*netDipoleZ;

    }

//  convert the quadrupole from traced to traceless form
 
    float qave   = (xxqdp + yyqdp + zzqdp)/3.0f;
          xxqdp  = 1.5f*(xxqdp-qave);
          xyqdp  = 1.5f*xyqdp;
          xzqdp  = 1.5f*xzqdp;
          yxqdp  = 1.5f*yxqdp;
          yyqdp  = 1.5f*(yyqdp-qave);
          yzqdp  = 1.5f*yzqdp;
          zxqdp  = 1.5f*zxqdp;
          zyqdp  = 1.5f*zyqdp;
          zzqdp  = 1.5f*(zzqdp-qave);

//  add the traceless atomic quadrupoles to total quadrupole

    for( unsigned int ii  = 0; ii < gpu->natoms; ii++ ){
        xxqdp  = xxqdp + 3.0f*labFrameQuadrupole[9*ii];
        xyqdp  = xyqdp + 3.0f*labFrameQuadrupole[9*ii+1];
        xzqdp  = xzqdp + 3.0f*labFrameQuadrupole[9*ii+2];
        yxqdp  = yxqdp + 3.0f*labFrameQuadrupole[9*ii+3];
        yyqdp  = yyqdp + 3.0f*labFrameQuadrupole[9*ii+4];
        yzqdp  = yzqdp + 3.0f*labFrameQuadrupole[9*ii+5];
        zxqdp  = zxqdp + 3.0f*labFrameQuadrupole[9*ii+6];
        zyqdp  = zyqdp + 3.0f*labFrameQuadrupole[9*ii+7];
        zzqdp  = zzqdp + 3.0f*labFrameQuadrupole[9*ii+8];
    }
 
    float debye                = 4.80321f;
    outputMultipoleMoments.resize( 13 );

    outputMultipoleMoments[0]  = netchg;

    outputMultipoleMoments[1]  = xdpl*debye;
    outputMultipoleMoments[2]  = ydpl*debye;
    outputMultipoleMoments[3]  = zdpl*debye;
    
    outputMultipoleMoments[4]  = xxqdp*debye;
    outputMultipoleMoments[5]  = xyqdp*debye;
    outputMultipoleMoments[6]  = xzqdp*debye;

    outputMultipoleMoments[7]  = yxqdp*debye;
    outputMultipoleMoments[8]  = yyqdp*debye;
    outputMultipoleMoments[9]  = yzqdp*debye;

    outputMultipoleMoments[10] = zxqdp*debye;
    outputMultipoleMoments[11] = zyqdp*debye;
    outputMultipoleMoments[12] = zzqdp*debye;

}
