#ifndef __KGBSA_H__
#define __KGBSA_H__

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Peter Eastman, Mark Friedrichs, Chris Bruns                       *
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



void kMergeFloat4( 
      const float repfac, 
      const float atomStrWidth, 
      const float pforceStrWidth,
      const float natoms,
      const float iUnroll,
      ::brook::stream stream1, 
      ::brook::stream stream2, 
      ::brook::stream outputStream 
                 );

void kMergeFloat4_4( 
      const float repfac, 
      const float atomStrWidth, 
      const float pforceStrWidth,
      const float natoms,
      const float roundNatoms,
      const float iUnroll,
      ::brook::stream stream1, 
      ::brook::stream stream2, 
      ::brook::stream stream3, 
      ::brook::stream stream4, 
      ::brook::stream outputStream 
                 );

void kMergeFloat4_4X( 
      const float repfac, 
      const float atomStrWidth, 
      const float pforceStrWidth,
      const float natoms,
      const float roundNatoms,
      const float iUnroll,
      ::brook::stream stream1, 
      ::brook::stream stream2, 
      ::brook::stream stream3, 
      ::brook::stream stream4, 
      ::brook::stream outputStream 
                 );

void kCheck( 
      const float natoms,
      const float atomStrWidth, 
      const float pforceStrWidth,
      const float unroll,
      ::brook::stream stream1
                 );

void kAddAndMergeFloat4( 
      const float repfac, 
      const float atomStrWidth, 
      const float pforceStrWidth,
      const float natoms,
      const float iUnroll,
      ::brook::stream inStream, 
      ::brook::stream stream1, 
      ::brook::stream stream2, 
      ::brook::stream outputStream 
                 );

void kAddAndMergeFloat4_4( 
      const float repfac, 
      const float atomStrWidth, 
      const float pforceStrWidth,
      const float natoms,
      const float roundNatoms,
      const float iUnroll,
      ::brook::stream inStream, 
      ::brook::stream stream1, 
      ::brook::stream stream2, 
      ::brook::stream stream3, 
      ::brook::stream stream4, 
      ::brook::stream outputStream 
                 );

void kPostObcLoop2_nobranch( 
      const float repfac, 
      const float atomStrWidth, 
      const float pforceStrWidth,
      const float natoms,
      const float roundNatoms,
      const float iUnroll,
      const float conversion,
      const float mergeNonObcForces,
      ::brook::stream inObcForces, 
      ::brook::stream nonObcForces, 
      ::brook::stream stream1, 
      ::brook::stream stream2, 
      ::brook::stream stream3, 
      ::brook::stream stream4, 
      ::brook::stream atomicRadii,
      ::brook::stream bornRadii,
      ::brook::stream obcChain,
      ::brook::stream outputStream 
                 );

void kPreGbsaForce2( 
      ::brook::stream intermediateForceIn, 
      ::brook::stream bornRadii, 
      ::brook::stream bornRadii2Force 
                 );

void kPreObcForce2( 
      ::brook::stream intermediateForceIn, 
      ::brook::stream obcChainForceIn, 
      ::brook::stream bornRadii, 
      ::brook::stream bornRadii2Force 
                 );

void kPostObcLoop1( 
      const float repfac, 
      const float atomStrWidth, 
      const float pforceStrWidth,
      const float natoms,
      const float roundNatoms,
      const float iUnroll,
      ::brook::stream stream1, 
      ::brook::stream stream2, 
      ::brook::stream stream3, 
      ::brook::stream stream4, 
      ::brook::stream obcChainForceIn, 
      ::brook::stream bornRadii, 
      ::brook::stream outputStream, 
      ::brook::stream bornRadii2Force 
                 );


void kPostObcLoop1_nobranch( 
      const float repfac, 
      const float atomStrWidth, 
      const float pforceStrWidth,
      const float natoms,
      const float roundNatoms,
      const float iUnroll,
      ::brook::stream stream1, 
      ::brook::stream stream2, 
      ::brook::stream stream3, 
      ::brook::stream stream4, 
      ::brook::stream obcChainForceIn, 
      ::brook::stream bornRadii, 
      ::brook::stream outputStream, 
      ::brook::stream bornRadii2Force 
                 );

void kCopyFloat4( 
      ::brook::stream stream1, 
      ::brook::stream stream2
                 );

void kObcLoop1(
      const float natoms, 
      const float roundedUpAtoms, 
      const float dupfac, 
      const float strwidth, 
      const float fstrwidth,
      const float soluteDielectric,
      const float solventDielectric,
      const float includeAce,
      ::brook::stream posq, 
      ::brook::stream bornRadii, 
      ::brook::stream atomicRadii, 
      ::brook::stream bornForce1,
      ::brook::stream bornForce2,
      ::brook::stream bornForce3,
      ::brook::stream bornForce4
      );

void kCalculateGbsaForce1_4 (
      const float natoms, 
      const float roundedUpAtoms, 
      const float dupfac, 
      const float strwidth, 
      const float fstrwidth,
      const float soluteDielectric,
      const float solventDielectric,
      ::brook::stream posq, 
      ::brook::stream bornRadii, 
      ::brook::stream nonpolarForce, 
      ::brook::stream bornForce1,
      ::brook::stream bornForce2,
      ::brook::stream bornForce3,
      ::brook::stream bornForce4
      );

void kCalculateGbsaForce1_4X (
      const float natoms, 
      const float roundedUpAtoms, 
      const float dupfac, 
      const float strwidth, 
      const float fstrwidth,
      const float soluteDielectric,
      const float solventDielectric,
      ::brook::stream posq, 
      ::brook::stream bornRadii, 
      ::brook::stream nonpolarForce, 
      ::brook::stream bornForce1,
      ::brook::stream bornForce2,
      ::brook::stream bornForce3,
      ::brook::stream bornForce4
      );

void kCalculateGbsaForce2_4 (
      float natoms, 
      float dupfac, 
      float strwidth, 
      float fstrwidth,
      ::brook::stream posq, 
      ::brook::stream  bornRadii,
      ::brook::stream vdwRadii, 
      ::brook::stream volume, 
      ::brook::stream excludeBlockIndices,
      ::brook::stream inBornForce,
      ::brook::stream bornForce1,
      ::brook::stream bornForce2
      );

void kCalculateGbsaForce2_4Debug (
      float natoms, 
      float dupfac, 
      float strwidth, 
      float fstrwidth,
      const float debugJAtom,
      const float debugReport,
      ::brook::stream posq, 
      ::brook::stream  bornRadii,
      ::brook::stream vdwRadii, 
      ::brook::stream volume, 
      ::brook::stream excludeBlockIndices,
      ::brook::stream inBornForce,
      ::brook::stream bornForce1,
      ::brook::stream bornForce2,
		::brook::stream debugStream1,
		::brook::stream debugStream2,
		::brook::stream debugStream3,
		::brook::stream debugStream4
      );

void kCalculateGbsaForce2_4Excl (
      float natoms, 
      float dupfac, 
      float strwidth, 
      float fstrwidth,
      ::brook::stream posq, 
      ::brook::stream  bornRadii,
      ::brook::stream vdwRadii, 
      ::brook::stream volume, 
      ::brook::stream excl,
      ::brook::stream excludeBlockIndices,
      ::brook::stream inBornForce,
      ::brook::stream bornForce1,
      ::brook::stream bornForce2
      );

void kCalculateGbsaForce2_4ExclDebug (
      float natoms, 
      float dupfac, 
      float strwidth, 
      float fstrwidth,
      const float debugJAtom,
      const float debugReport,
      ::brook::stream posq, 
      ::brook::stream  bornRadii,
      ::brook::stream vdwRadii, 
      ::brook::stream volume, 
      ::brook::stream excl,
      ::brook::stream excludeBlockIndices,
      ::brook::stream inBornForce,
      ::brook::stream bornForce1,
      ::brook::stream bornForce2,
		::brook::stream debugStream1,
		::brook::stream debugStream2,
		::brook::stream debugStream3,
		::brook::stream debugStream4
      );

void kCalculateGbsaForce2_2 (
      float natoms, 
      float dupfac, 
      float strwidth, 
      float fstrwidth,
      ::brook::stream posq, 
      ::brook::stream  bornRadii,
      ::brook::stream vdwRadii, 
      ::brook::stream volume, 
      ::brook::stream excludeBlockIndices,
      ::brook::stream inBornForce,
      ::brook::stream bornForce1,
      ::brook::stream bornForce2
      );

void kCalculateGbsaForce2_2Debug (
      float natoms, 
      float dupfac, 
      float strwidth, 
      float fstrwidth,
      const float debugJAtom,
      const float debugReport,
      ::brook::stream posq, 
      ::brook::stream  bornRadii,
      ::brook::stream vdwRadii, 
      ::brook::stream volume, 
      ::brook::stream excludeBlockIndices,
      ::brook::stream inBornForce,
      ::brook::stream bornForce1,
      ::brook::stream bornForce2,
		::brook::stream debugStream1,
		::brook::stream debugStream2,
		::brook::stream debugStream3,
		::brook::stream debugStream4
      );

void kCalculateGbsaForce2_2Excl (
      float natoms, 
      float dupfac, 
      float strwidth, 
      float fstrwidth,
      ::brook::stream posq, 
      ::brook::stream  bornRadii,
      ::brook::stream vdwRadii, 
      ::brook::stream volume, 
      ::brook::stream excl,
      ::brook::stream excludeBlockIndices,
      ::brook::stream inBornForce,
      ::brook::stream bornForce1,
      ::brook::stream bornForce2
      );

void kCalculateGbsaForce2_2ExclDebug (
      float natoms, 
      float dupfac, 
      float strwidth, 
      float fstrwidth,
      const float debugJAtom,
      const float debugReport,
      ::brook::stream posq, 
      ::brook::stream  bornRadii,
      ::brook::stream vdwRadii, 
      ::brook::stream volume, 
      ::brook::stream excl,
      ::brook::stream excludeBlockIndices,
      ::brook::stream inBornForce,
      ::brook::stream bornForce1,
      ::brook::stream bornForce2,
		::brook::stream debugStream1,
		::brook::stream debugStream2,
		::brook::stream debugStream3,
		::brook::stream debugStream4
      );

void kBornRadii(
      ::brook::stream gpolNonBonded, 
      ::brook::stream gPolFixed,
      ::brook::stream bornRadii
      );

void kPostGbsaForce2Debug (
      float natoms, 
      float dupfac, 
      float strwidth, 
      float fstrwidth,
      const float debugJAtom,
      const float debugReport,
      ::brook::stream inBornForce2, 
		::brook::stream gPolFixed,
      ::brook::stream  bornRadii,
		::brook::stream debugStream1,
		::brook::stream debugStream2,
		::brook::stream debugStream3,
		::brook::stream debugStream4
      );

void kCalculateGbsaForce2_1_4 (
      float natoms, 
      float dupfac, 
      float strwidth, 
      float fstrwidth,
      ::brook::stream posq, 
      ::brook::stream  bornRadii,
      ::brook::stream vdwRadii, 
      ::brook::stream volume, 
      ::brook::stream excludeBlockIndices,
      ::brook::stream inBornForce,
      ::brook::stream bornForce1
      );

void kCalculateGbsaForce2_1_4Debug (
      float natoms, 
      float dupfac, 
      float strwidth, 
      float fstrwidth,
      const float debugJAtom,
      const float debugReport,
      ::brook::stream posq, 
      ::brook::stream  bornRadii,
      ::brook::stream vdwRadii, 
      ::brook::stream volume, 
      ::brook::stream excludeBlockIndices,
      ::brook::stream inBornForce,
      ::brook::stream bornForce1,
		::brook::stream debugStream1,
		::brook::stream debugStream2,
		::brook::stream debugStream3,
		::brook::stream debugStream4
      );

void kCalculateGbsaForce2_1_4Excl (
      float natoms, 
      float dupfac, 
      float strwidth, 
      float fstrwidth,
      ::brook::stream posq, 
      ::brook::stream  bornRadii,
      ::brook::stream vdwRadii, 
      ::brook::stream volume, 
      ::brook::stream excl,
      ::brook::stream excludeBlockIndices,
      ::brook::stream inBornForce,
      ::brook::stream bornForce1
      );

void kCalculateGbsaForce2_1_4ExclDebug (
      float natoms, 
      float dupfac, 
      float strwidth, 
      float fstrwidth,
      const float debugIAtom,
      const float debugJAtom,
      const float debugReport,
      ::brook::stream posq, 
      ::brook::stream  bornRadii,
      ::brook::stream vdwRadii, 
      ::brook::stream volume, 
      ::brook::stream excl,
      ::brook::stream excludeBlockIndices,
      ::brook::stream inBornForce,
      ::brook::stream bornForce1,
		::brook::stream debugStream1,
		::brook::stream debugStream2,
		::brook::stream debugStream3,
		::brook::stream debugStream4
      );

void kObcLoop2(
      float numberOfAtoms, 
      float roundNatoms,
      float duplicationFactor, 
      float streamWidth, 
      float fstreamWidth,
      ::brook::stream posq, 
      ::brook::stream atomicRadii, 
      ::brook::stream bornForceFactor, 
      ::brook::stream bornForce1,
      ::brook::stream bornForce2,
      ::brook::stream bornForce3,
      ::brook::stream bornForce4
      );

void kObcBornRadii(
      float numberOfAtoms, 
      float streamWidth, 
      ::brook::stream bornForce, 
      ::brook::stream atomicRadii, 
      ::brook::stream bornRadii,
      ::brook::stream obcChain
      );

void kAceNonPolar(
      float soluteDielectric,
      float solventDielectric,
      ::brook::stream posq, 
      ::brook::stream bornRadii, 
      ::brook::stream vdwRadii, 
      ::brook::stream bornForce
      );

void kZeroFloat4( ::brook::stream bornForce );
void kSetValue4( float value, ::brook::stream bornForce );
void kSetValue3( float value, ::brook::stream bornForce );
void kSetValue2( float value, ::brook::stream bornForce );
void kSetValue1( float value, ::brook::stream bornForce );

void kCopyFloat3To4( ::brook::stream inForce, ::brook::stream outForce );

void kObcLoop2Cdlj4( float numberOfAtoms, float roundedUpAtoms, float duplicationFactor, 
                     float streamWidth, float fstreamWidth,
                     float jstreamWidth, float epsfac, 
                     ::brook::stream posq, 
                     ::brook::stream atomicRadii,
                     ::brook::stream scaledAtomicRadii, 
                     ::brook::stream bornForceFactor,
                     ::brook::stream isigeps,
                     ::brook::stream sigma,
                     ::brook::stream epsilon,
                     ::brook::stream excl,
                     ::brook::stream bornForce1,
                     ::brook::stream bornForce2,
                     ::brook::stream bornForce3,
                     ::brook::stream bornForce4 );

#endif //__KGBSA_H__
