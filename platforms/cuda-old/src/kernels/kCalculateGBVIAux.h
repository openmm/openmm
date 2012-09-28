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

#ifndef __GpuGBVIAUX_H__
#define __GpuGBVIAUX_H__

/**
 * This file contains subroutines used in evaluating quantities associated w/ the GB/VI function
 */

static __device__ float getGBVI_L( float r, float x, float S )
{

   float rInv   = 1.0f/r;
   float xInv   = 1.0f/x;

   float xInv2  = xInv*xInv;
   float diff2  = (r + S)*(r - S);

   return (1.5f*xInv2)*( (0.25f*rInv) - (xInv/3.0f) + (0.125f*diff2*xInv2*rInv) );
}

static __device__ float getGBVI_Volume( float r, float R, float S )
{

    float addOn      = 0.0f;
    int mask         = 1;
    float lowerBound = (r - S); 

    float diff       = (S - R); 
    if( fabsf( diff ) < r ){
        lowerBound = R > lowerBound ? R : lowerBound; 
    } else if( r <= diff ){
        addOn      = (1.0f/(R*R*R));
    } else {
        mask       = 0;
    }   
    float s2         = getGBVI_L( r, lowerBound, S );
    float s1         = getGBVI_L( r, (r + S),    S );
    s1               = mask ? (s1 - s2 + addOn) : 0.0f;
    return s1;
}

static __device__ float getGBVI_dL_dr( float r, float x, float S )
{

   float rInv   = 1.0f/r;
   float rInv2  = rInv*rInv;

   float xInv   = 1.0f/x;
   float xInv2  = xInv*xInv;
   float xInv3  = xInv2*xInv;
   
   float diff2  = (r + S)*(r - S);

   return ( (-1.5f*xInv2*rInv2)*( 0.25f + 0.125f*diff2*xInv2 ) + 0.375f*xInv3*xInv );

}

static __device__ float getGBVI_dL_drNew( float r, float x, float S )
{

   float rInv   = 1.0f/r;
   float rInv2  = rInv*rInv;

   float xInv   = 1.0f/x;
   float xInv2  = xInv*xInv;
   
   float t1     = (S*rInv);
         t1     = 1.0f + t1*t1;

   return (-0.375f*xInv2)*( rInv2 - 0.5f*xInv2*t1 );

}

static __device__ float getGBVI_dL_dx( float r, float x, float S )
{

   float rInv   = 1.0f/r;

   float xInv   = 1.0f/x;
   float xInv2  = xInv*xInv;
   float xInv3  = xInv2*xInv;

   float diff   = (r + S)*(r - S);

   return ( (-1.5f*xInv3)*( (0.5f*rInv) - xInv + (0.5f*diff*xInv2*rInv) ));

}

static __device__ float getGBVI_dE2Old( float r, float R, float S, float bornForce )
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
    float dE2              = getGBVI_dL_dr( r, lowerBound, S ) + mask*getGBVI_dL_dx( r, lowerBound, S );
    dE                    -= (absDiff >= r) && r >= diff ? 0.0f : dE2;
    dE                     = r < -diff ? 0.0f : dE;

    dE                    *= ( (r > 1.0e-08f) ? (bornForce/r) : 0.0f);

    return (-dE);

}

static __device__ float getGBVI_dE2( float r, float R, float S, float bornForce )
{
    float diff  = S - R;
    float dE    = 0.0f;

    if( fabsf( diff ) < r ){
        dE = getGBVI_dL_dr( r, r+S, S ) + getGBVI_dL_dx( r, r+S, S );   
        float lowerBound;
        float mask;
        if( R > (r - S) ){
            lowerBound = R;
            mask       = 0.0f;
        } else {
            lowerBound = r - S;
            mask       = 1.0f;
        }
        dE -= getGBVI_dL_dr( r, lowerBound, S ) + mask*getGBVI_dL_dx( r, lowerBound, S );
    } else if( r < (S - R) ){
        dE  =   getGBVI_dL_dr( r, r+S, S ) + getGBVI_dL_dx( r, r+S, S );   
        dE -= ( getGBVI_dL_dr( r, r-S, S ) + getGBVI_dL_dx( r, r-S, S ) );   
    }   
   
    dE         *= ( (r > 1.0e-08f) ? (bornForce/r) : 0.0f);
    return (-dE);

}

static __device__ float getGBVIBornForce2( float bornRadius, float R, float bornForce, float gamma )
{ 
    float ratio                     = (R/bornRadius);
    float returnBornForce           = bornForce + (3.0f*gamma*ratio*ratio*ratio)/bornRadius; // 'cavity' term
    float br2                       = bornRadius*bornRadius;
          returnBornForce          *= (1.0f/3.0f)*br2*br2;

   return returnBornForce;

}

#endif // __GpuGBVIAUX_H__
