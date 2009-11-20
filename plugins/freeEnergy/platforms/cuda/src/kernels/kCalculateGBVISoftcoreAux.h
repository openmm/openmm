/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Mark Friedrichs                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#ifndef __Gpu_GBVI_SOFTCORE_AUX_H__
#define __Gpu_GBVI_SOFTCORE_AUX_H__

/**
 * This file contains subroutines used in evaluating quantities associated w/ the GB/VI function
 */

__device__ float getGBVI_L( float r, float x, float S )
{

   float rInv   = 1.0f/r;

   float xInv   = 1.0f/x;
   float xInv2  = xInv*xInv;
   float diff2  = (r + S)*(r - S);

   return (1.5f*xInv2)*( (0.25f*rInv) - (xInv/3.0f) + (0.125f*diff2*xInv2*rInv) );
}

__device__ float getGBVI_Volume( float r_ij, float R, float S )
{

     float upperBound        = r_ij + S; 
     float rdiffS            = r_ij - S; 
     float lowerBound        = R > rdiffS ? R : rdiffS;
     float L_upper           = getGBVI_L( r_ij, upperBound, S );
     float L_lower           = getGBVI_L( r_ij, lowerBound, S );
     float mask              = r_ij < (R - S) ? 0.0f : 1.0f;  
     float addOn             = r_ij < (S - R) ? (1.0f/(R*R*R)) : 0.0f;  

     return (mask*( L_upper - L_lower ) + addOn);

}

__device__ float getGBVI_dL_dr( float r, float x, float S )
{

   float rInv   = 1.0f/r;
   float rInv2  = rInv*rInv;

   float xInv   = 1.0f/x;
   float xInv2  = xInv*xInv;
   float xInv3  = xInv2*xInv;
   
   float diff2  = (r + S)*(r - S);

   return ( (-1.5f*xInv2*rInv2)*( 0.25f + 0.125f*diff2*xInv2 ) + 0.375f*xInv3*xInv );
   //return 0.0f;

}

__device__ float getGBVI_dL_dx( float r, float x, float S )
{

   float rInv   = 1.0f/r;

   float xInv   = 1.0f/x;
   float xInv2  = xInv*xInv;
   float xInv3  = xInv2*xInv;

   float diff   = (r + S)*(r - S);

   return ( (-1.5f*xInv3)*( (0.5f*rInv) - xInv + (0.5f*diff*xInv2*rInv) ));


}

__device__ float getGBVI_dE2( float r, float R, float S, float bornForce )
{

    float diff              = S - R;
    float absDiff           = fabsf( S - R );
    float dE                = getGBVI_dL_dr( r, r+S, S ) + getGBVI_dL_dx( r, r+S, S );
    float mask;
    float lowerBound;
    if( (R > (r - S)) && (absDiff < r) ){
        mask       = 0.0f;
        lowerBound = R;
    } else {
        mask       = 1.0f;
        lowerBound = (r - S);  
    }   
    dE                    -= getGBVI_dL_dr( r, lowerBound, S ) + mask*getGBVI_dL_dx( r, lowerBound, S );
    dE                     = (absDiff >= r) && r >= diff ? 0.0f : dE; 
    dE                    *= ( (r > 1.0e-08f) ? (bornForce/r) : 0.0f);

    return (-dE);

}

__device__ float getGBVIBornForce2( float bornRadius, float R, float bornForce, float gamma )
{ 
    float ratio                     = (R/bornRadius);
    float returnBornForce           = bornForce + (3.0f*gamma*ratio*ratio*ratio)/bornRadius; // 'cavity' term
    float br2                       = bornRadius*bornRadius;
          returnBornForce          *= (1.0f/3.0f)*br2*br2;

   return returnBornForce;

}

#endif // __Gpu_GBVI_SOFTCORE_AUX_H__
