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

#include <stdio.h>
#include <string.h>
#include <cuda.h>
#include <vector_functions.h>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cmath>
#include <map>
#include <algorithm>
#ifdef WIN32
  #include <windows.h>
#else
  #include <stdint.h>
#endif
using namespace std;

#include "gputypes.h"
#include "cudaKernels.h"
#include "OpenMMException.h"

using OpenMM::OpenMMException;

#ifdef WIN32
  typedef unsigned __int64 u64;
  typedef signed __int64 s64;
#else
  typedef uint64_t u64;
  typedef int64_t s64;
#endif
typedef unsigned int u32;
typedef float f32;
typedef double f64;
typedef char ascii;
typedef char utf8;
typedef unsigned char u8;
typedef signed char s8;
typedef unsigned short u16;
typedef signed short s16;
typedef struct
{
  u8 type[4];
  f32 charge;
  f32 radius;
} FAH_ATOM;

typedef struct
{
  u32 a; /* rule: a < b */
  u32 b;
} FAH_BOND;

typedef struct
{
  f32 x;
  f32 y;
  f32 z;
} FAH_XYZ;

typedef struct
{
  u32    magic;
  u32    version;
  utf8   name[64];
  s64    timestamp;
  u64    iterations;
  u32    frames;
  u32    atom_count;
  u32    bond_count;
  /* v2 */
  utf8   user_name[64];
  utf8   user_team[16];
  utf8   user_done[16];
} FAH_INFO;

typedef struct
{
  u32    magic;
  u32    version;
  s64    timestamp;
  u64    iterations_done;
  u32    frames_done;
  f32    energy;
  f32    temperature;
} FAH_CURRENT;

typedef struct
{
  FAH_INFO info;
  FAH_CURRENT current;
  FAH_ATOM * atoms;
  FAH_BOND * bonds;
  FAH_XYZ  * xyz;
} PROTEIN;

struct ShakeCluster {
    int centralID;
    int peripheralID[3];
    int size;
    float distance;
    float centralInvMass, peripheralInvMass;
    ShakeCluster() {
    }
    ShakeCluster(int centralID, float invMass) : centralID(centralID), centralInvMass(invMass), size(0) {
    }
    void addAtom(int id, float dist, float invMass) {
        if (size == 3)
            throw OpenMMException("A single atom may only have three constraints");
        if (size > 0 && dist != distance)
            throw OpenMMException("All constraints for a central atom must have the same distance");
        if (size > 0 && invMass != peripheralInvMass)
            throw OpenMMException("All constraints for a central atom must have the same mass");
        peripheralID[size++] = id;
        distance = dist;
        peripheralInvMass = invMass;
    }
};

struct Constraint
{
    Constraint(int atom1, int atom2, float distance2) : atom1(atom1), atom2(atom2), distance2(distance2) {
    }
    int atom1, atom2;
    float distance2;
};

struct Molecule {
    vector<int> atoms;
    vector<int> bonds;
    vector<int> angles;
    vector<int> periodicTorsions;
    vector<int> rbTorsions;
    vector<int> constraints;
};

static const float dielectricOffset         =    0.009f;
static const float PI                       =    3.1415926535f;
static const float probeRadius              =    0.14f;
static const float forceConversionFactor    =    0.4184f;

//static const float surfaceAreaFactor        =   -6.0f * 0.06786f * forceConversionFactor * 1000.0f;  // PI * 4.0f * 0.0049f * 1000.0f;
//static const float surfaceAreaFactor        =   -6.0f * PI * 4.0f * 0.0049f * 1000.0f;
static const float surfaceAreaFactor        = -6.0f*PI*0.0216f*1000.0f*0.4184f;
//static const float surfaceAreaFactor        = -1.7035573959e+001;
//static const float surfaceAreaFactor        = -166.02691f;
//static const float surfaceAreaFactor        = 1.0f;

static const float alphaOBC                 =    1.0f;
static const float betaOBC                  =    0.8f;
static const float gammaOBC                 =    4.85f;
static const float kcalMolTokJNM            =   -0.4184f;
static const float electricConstant         = -166.02691f;
static const float defaultInnerDielectric   =    1.0f;
static const float defaultSolventDielectric =   78.3f;
static const float KILO                     =    1e3;                      // Thousand
static const float BOLTZMANN                =    1.380658e-23f;            // (J/K)    
static const float AVOGADRO                 =    6.0221367e23f;            // ()        
static const float RGAS                     =    BOLTZMANN * AVOGADRO;     // (J/(mol K))
static const float BOLTZ                    =    (RGAS / KILO);            // (kJ/(mol K)) 

#define DUMP_PARAMETERS 0

#define DeltaShake


extern "C"
int gpuReadBondParameters(gpuContext gpu, char* fname)
{
    ifstream infile(fname);
    
    if (!infile.fail())
    {
        char buff[512];
        int bonds;
        infile >> bonds;
        infile.getline(buff, 512);
        vector<int> atom1(bonds);
        vector<int> atom2(bonds);
        vector<float> length(bonds);
        vector<float> k(bonds);
        for (int i = 0; i < bonds; i++)
        {
            int junk;
            infile >> 
                junk >> 
                atom1[i] >> 
                atom2[i] >> 
                length[i] >> 
                k[i];
        }
        gpuSetBondParameters(gpu, atom1, atom2, length, k);
        return bonds;
    }
    else
    {
        cout << "Error opening harmonic bond parameter file " << fname << endl;
        exit(-1);
    }
    return 0;
}

extern "C"
void gpuSetBondParameters(gpuContext gpu, const vector<int>& atom1, const vector<int>& atom2, const vector<float>& length, const vector<float>& k)
{
    int bonds = atom1.size();
    gpu->sim.bonds                              = bonds;
    CUDAStream<int4>* psBondID                  = new CUDAStream<int4>(bonds, 1);
    gpu->psBondID                               = psBondID;
    gpu->sim.pBondID                            = psBondID->_pDevStream[0];
    CUDAStream<float2>* psBondParameter         = new CUDAStream<float2>(bonds, 1);
    gpu->psBondParameter                        = psBondParameter;
    gpu->sim.pBondParameter                     = psBondParameter->_pDevStream[0];
    for (int i = 0; i < bonds; i++)
    {
        psBondID->_pSysStream[0][i].x = atom1[i];
        psBondID->_pSysStream[0][i].y = atom2[i];
        psBondParameter->_pSysStream[0][i].x = length[i];
        psBondParameter->_pSysStream[0][i].y = k[i];
        psBondID->_pSysStream[0][i].z = gpu->pOutputBufferCounter[psBondID->_pSysStream[0][i].x]++;
        psBondID->_pSysStream[0][i].w = gpu->pOutputBufferCounter[psBondID->_pSysStream[0][i].y]++;
#if (DUMP_PARAMETERS == 1)                
        cout << 
            i << " " << 
            psBondID->_pSysStream[0][i].x << " " << 
            psBondID->_pSysStream[0][i].y << " " << 
            psBondID->_pSysStream[0][i].z << " " << 
            psBondID->_pSysStream[0][i].w << " " << 
            psBondParameter->_pSysStream[0][i].x << " " << 
            psBondParameter->_pSysStream[0][i].y << 
            endl;
#endif
    }
    psBondID->Upload();
    psBondParameter->Upload();
}

extern "C"
int gpuReadBondAngleParameters(gpuContext gpu, char* fname)
{
    ifstream infile(fname);
    
    if (!infile.fail())
    {
        char buff[512];
        int bond_angles;
        infile >> bond_angles;
        infile.getline(buff, 512);
        vector<int> atom1(bond_angles);
        vector<int> atom2(bond_angles);
        vector<int> atom3(bond_angles);
        vector<float> angle(bond_angles);
        vector<float> k(bond_angles);
     
        for (int i = 0; i < bond_angles; i++)
        {
            int junk;
            infile >> 
                junk >> 
                atom1[i] >> 
                atom2[i] >> 
                atom3[i] >> 
                angle[i] >> 
                k[i];
        }
        gpuSetBondAngleParameters(gpu, atom1, atom2, atom3, angle, k);
        return bond_angles;
    }
    else
    {
        cout << "Error opening harmonic bond angle parameter file " << fname << endl;
        exit(-1);
    }
    return 0;
}

extern "C"
void gpuSetBondAngleParameters(gpuContext gpu, const vector<int>& atom1, const vector<int>& atom2, const vector<int>& atom3,
        const vector<float>& angle, const vector<float>& k)
{
    int bond_angles = atom1.size();
    gpu->sim.bond_angles                        = bond_angles;
    CUDAStream<int4>* psBondAngleID1            = new CUDAStream<int4>(bond_angles, 1);
    gpu->psBondAngleID1                         = psBondAngleID1;
    gpu->sim.pBondAngleID1                      = psBondAngleID1->_pDevStream[0];
    CUDAStream<int2>* psBondAngleID2            = new CUDAStream<int2>(bond_angles, 1);
    gpu->psBondAngleID2                         = psBondAngleID2;
    gpu->sim.pBondAngleID2                      = psBondAngleID2->_pDevStream[0];
    CUDAStream<float2>* psBondAngleParameter    = new CUDAStream<float2>(bond_angles, 1);
    gpu->psBondAngleParameter                   = psBondAngleParameter;
    gpu->sim.pBondAngleParameter                = psBondAngleParameter->_pDevStream[0];        

    for (int i = 0; i < bond_angles; i++)
    {
        psBondAngleID1->_pSysStream[0][i].x = atom1[i];
        psBondAngleID1->_pSysStream[0][i].y = atom2[i];
        psBondAngleID1->_pSysStream[0][i].z = atom3[i];
        psBondAngleParameter->_pSysStream[0][i].x = angle[i];
        psBondAngleParameter->_pSysStream[0][i].y = k[i];
        psBondAngleID1->_pSysStream[0][i].w = gpu->pOutputBufferCounter[psBondAngleID1->_pSysStream[0][i].x]++;
        psBondAngleID2->_pSysStream[0][i].x = gpu->pOutputBufferCounter[psBondAngleID1->_pSysStream[0][i].y]++;
        psBondAngleID2->_pSysStream[0][i].y = gpu->pOutputBufferCounter[psBondAngleID1->_pSysStream[0][i].z]++;
#if (DUMP_PARAMETERS == 1)
         cout << 
            i << " " << 
            psBondAngleID1->_pSysStream[0][i].x << " " << 
            psBondAngleID1->_pSysStream[0][i].y << " " << 
            psBondAngleID1->_pSysStream[0][i].z << " " << 
            psBondAngleID1->_pSysStream[0][i].w << " " << 
            psBondAngleID2->_pSysStream[0][i].x << " " << 
            psBondAngleID2->_pSysStream[0][i].y << " " << 
            psBondAngleParameter->_pSysStream[0][i].x << " " << 
            psBondAngleParameter->_pSysStream[0][i].y << 
            endl;
#endif
    }
    psBondAngleID1->Upload();
    psBondAngleID2->Upload();
    psBondAngleParameter->Upload();
}

extern "C"
int gpuReadDihedralParameters(gpuContext gpu, char* fname)
{
    ifstream infile(fname);
    
    if (!infile.fail())
    {
        char buff[512];
        int dihedrals;
        infile >> dihedrals;
        infile.getline(buff, 512);
        vector<int> atom1(dihedrals);
        vector<int> atom2(dihedrals);
        vector<int> atom3(dihedrals);
        vector<int> atom4(dihedrals);
        vector<float> k(dihedrals);
        vector<float> phase(dihedrals);
        vector<int> periodicity(dihedrals);
        for (int i = 0; i < dihedrals; i++)
        {
            int junk;
            infile >> 
                junk >> 
                atom1[i] >> 
                atom2[i] >> 
                atom3[i] >>
                atom4[i] >> 
                k[i] >> 
                phase[i] >>
                periodicity[i];
        }
        gpuSetDihedralParameters(gpu, atom1, atom2, atom3, atom4, k, phase, periodicity);
        return dihedrals;
    }
    else
    {
        cout << "Error opening dihedral parameter file " << fname << endl;
        exit(-1);
    }
    return 0;
}

extern "C"
void gpuSetDihedralParameters(gpuContext gpu, const vector<int>& atom1, const vector<int>& atom2, const vector<int>& atom3, const vector<int>& atom4,
        const vector<float>& k, const vector<float>& phase, const vector<int>& periodicity)
{
        int dihedrals = atom1.size();
        gpu->sim.dihedrals = dihedrals;
        CUDAStream<int4>* psDihedralID1             = new CUDAStream<int4>(dihedrals, 1);
        gpu->psDihedralID1                          = psDihedralID1;
        gpu->sim.pDihedralID1                       = psDihedralID1->_pDevStream[0];
        CUDAStream<int4>* psDihedralID2             = new CUDAStream<int4>(dihedrals, 1);
        gpu->psDihedralID2                          = psDihedralID2;
        gpu->sim.pDihedralID2                       = psDihedralID2->_pDevStream[0];
        CUDAStream<float4>* psDihedralParameter     = new CUDAStream<float4>(dihedrals, 1);
        gpu->psDihedralParameter                    = psDihedralParameter;
        gpu->sim.pDihedralParameter                 = psDihedralParameter->_pDevStream[0];
        for (int i = 0; i < dihedrals; i++)
        {
            psDihedralID1->_pSysStream[0][i].x = atom1[i];
            psDihedralID1->_pSysStream[0][i].y = atom2[i];
            psDihedralID1->_pSysStream[0][i].z = atom3[i];
            psDihedralID1->_pSysStream[0][i].w = atom4[i];
            psDihedralParameter->_pSysStream[0][i].x = k[i];
            psDihedralParameter->_pSysStream[0][i].y = phase[i];
            psDihedralParameter->_pSysStream[0][i].z = (float) periodicity[i];
            psDihedralID2->_pSysStream[0][i].x = gpu->pOutputBufferCounter[psDihedralID1->_pSysStream[0][i].x]++;
            psDihedralID2->_pSysStream[0][i].y = gpu->pOutputBufferCounter[psDihedralID1->_pSysStream[0][i].y]++;
            psDihedralID2->_pSysStream[0][i].z = gpu->pOutputBufferCounter[psDihedralID1->_pSysStream[0][i].z]++;
            psDihedralID2->_pSysStream[0][i].w = gpu->pOutputBufferCounter[psDihedralID1->_pSysStream[0][i].w]++;
#if (DUMP_PARAMETERS == 1)
            cout << 
                i << " " << 
                psDihedralID1->_pSysStream[0][i].x << " " << 
                psDihedralID1->_pSysStream[0][i].y << " " << 
                psDihedralID1->_pSysStream[0][i].z << " " << 
                psDihedralID1->_pSysStream[0][i].w << " " << 
                psDihedralID2->_pSysStream[0][i].x << " " << 
                psDihedralID2->_pSysStream[0][i].y << " " << 
                psDihedralID2->_pSysStream[0][i].z << " " << 
                psDihedralID2->_pSysStream[0][i].w << " " << 
                psDihedralParameter->_pSysStream[0][i].x << " " << 
                psDihedralParameter->_pSysStream[0][i].y << " " << 
                psDihedralParameter->_pSysStream[0][i].z << endl;
#endif
        }
        psDihedralID1->Upload();
        psDihedralID2->Upload();
        psDihedralParameter->Upload();
}

extern "C"
int gpuReadRbDihedralParameters(gpuContext gpu, char* fname)
{
    ifstream infile(fname);
    
    if (!infile.fail())
    {
        char buff[512];
        int rb_dihedrals;
        infile >> rb_dihedrals;
        infile.getline(buff, 512);
        vector<int> atom1(rb_dihedrals);
        vector<int> atom2(rb_dihedrals);
        vector<int> atom3(rb_dihedrals);
        vector<int> atom4(rb_dihedrals);
        vector<float> c0(rb_dihedrals);
        vector<float> c1(rb_dihedrals);
        vector<float> c2(rb_dihedrals);
        vector<float> c3(rb_dihedrals);
        vector<float> c4(rb_dihedrals);
        vector<float> c5(rb_dihedrals);
        gpu->sim.rb_dihedrals = rb_dihedrals;
        CUDAStream<int4>* psRbDihedralID1           = new CUDAStream<int4>(rb_dihedrals, 1);
        gpu->psRbDihedralID1                        = psRbDihedralID1;
        gpu->sim.pRbDihedralID1                     = psRbDihedralID1->_pDevStream[0];
        CUDAStream<int4>* psRbDihedralID2           = new CUDAStream<int4>(rb_dihedrals, 1);
        gpu->psRbDihedralID2                        = psRbDihedralID2;
        gpu->sim.pRbDihedralID2                     = psRbDihedralID2->_pDevStream[0];
        CUDAStream<float4>* psRbDihedralParameter1  = new CUDAStream<float4>(rb_dihedrals, 1);
        gpu->psRbDihedralParameter1                 = psRbDihedralParameter1;
        gpu->sim.pRbDihedralParameter1              = psRbDihedralParameter1->_pDevStream[0];
        CUDAStream<float2>* psRbDihedralParameter2  = new CUDAStream<float2>(rb_dihedrals, 1);    
        gpu->psRbDihedralParameter2                 = psRbDihedralParameter2;
        gpu->sim.pRbDihedralParameter2              = psRbDihedralParameter2->_pDevStream[0];
        
        for (int i = 0; i < rb_dihedrals; i++)
        {
            int junk;
            infile >> 
                junk >> 
                atom1[i] >> 
                atom2[i] >> 
                atom3[i] >>
                atom4[i] >> 
                c0[i] >> 
                c1[i] >> 
                c2[i] >> 
                c3[i] >> 
                c4[i] >> 
                c5[i];
        }
        gpuSetRbDihedralParameters(gpu, atom1, atom2, atom3, atom4, c0, c1, c2, c3, c4, c5);
        return rb_dihedrals;
    }
    else
    {
        cout << "Error opening Ryckaert-Bellemans dihedral parameter file " << fname << endl;
        exit(-1);
    }
    return 0;
}

extern "C"
void gpuSetRbDihedralParameters(gpuContext gpu, const vector<int>& atom1, const vector<int>& atom2, const vector<int>& atom3, const vector<int>& atom4,
        const vector<float>& c0, const vector<float>& c1, const vector<float>& c2, const vector<float>& c3, const vector<float>& c4, const vector<float>& c5)
{
    int rb_dihedrals = atom1.size();
    gpu->sim.rb_dihedrals = rb_dihedrals;
    CUDAStream<int4>* psRbDihedralID1           = new CUDAStream<int4>(rb_dihedrals, 1);
    gpu->psRbDihedralID1                        = psRbDihedralID1;
    gpu->sim.pRbDihedralID1                     = psRbDihedralID1->_pDevStream[0];
    CUDAStream<int4>* psRbDihedralID2           = new CUDAStream<int4>(rb_dihedrals, 1);
    gpu->psRbDihedralID2                        = psRbDihedralID2;
    gpu->sim.pRbDihedralID2                     = psRbDihedralID2->_pDevStream[0];
    CUDAStream<float4>* psRbDihedralParameter1  = new CUDAStream<float4>(rb_dihedrals, 1);
    gpu->psRbDihedralParameter1                 = psRbDihedralParameter1;
    gpu->sim.pRbDihedralParameter1              = psRbDihedralParameter1->_pDevStream[0];
    CUDAStream<float2>* psRbDihedralParameter2  = new CUDAStream<float2>(rb_dihedrals, 1);    
    gpu->psRbDihedralParameter2                 = psRbDihedralParameter2;
    gpu->sim.pRbDihedralParameter2              = psRbDihedralParameter2->_pDevStream[0];

    for (int i = 0; i < rb_dihedrals; i++)
    {
        psRbDihedralID1->_pSysStream[0][i].x = atom1[i];
        psRbDihedralID1->_pSysStream[0][i].y = atom2[i];
        psRbDihedralID1->_pSysStream[0][i].z = atom3[i];
        psRbDihedralID1->_pSysStream[0][i].w = atom4[i];
        psRbDihedralParameter1->_pSysStream[0][i].x = c0[i];
        psRbDihedralParameter1->_pSysStream[0][i].y = c1[i];
        psRbDihedralParameter1->_pSysStream[0][i].z = c2[i];
        psRbDihedralParameter1->_pSysStream[0][i].w = c3[i];
        psRbDihedralParameter2->_pSysStream[0][i].x = c4[i];
        psRbDihedralParameter2->_pSysStream[0][i].y = c5[i];
        psRbDihedralID2->_pSysStream[0][i].x = gpu->pOutputBufferCounter[psRbDihedralID1->_pSysStream[0][i].x]++;
        psRbDihedralID2->_pSysStream[0][i].y = gpu->pOutputBufferCounter[psRbDihedralID1->_pSysStream[0][i].y]++;
        psRbDihedralID2->_pSysStream[0][i].z = gpu->pOutputBufferCounter[psRbDihedralID1->_pSysStream[0][i].z]++;
        psRbDihedralID2->_pSysStream[0][i].w = gpu->pOutputBufferCounter[psRbDihedralID1->_pSysStream[0][i].w]++;
#if (DUMP_PARAMETERS == 1)
        cout << 
            i << " " << 
            psRbDihedralID1->_pSysStream[0][i].x << " " << 
            psRbDihedralID1->_pSysStream[0][i].y << " " << 
            psRbDihedralID1->_pSysStream[0][i].z << " " << 
            psRbDihedralID1->_pSysStream[0][i].w <<" " << 
            psRbDihedralID2->_pSysStream[0][i].x << " " << 
            psRbDihedralID2->_pSysStream[0][i].y << " " << 
            psRbDihedralID2->_pSysStream[0][i].z << " " << 
            psRbDihedralID2->_pSysStream[0][i].w <<" " <<                 
            psRbDihedralParameter1->_pSysStream[0][i].x << " " << 
            psRbDihedralParameter1->_pSysStream[0][i].y << " " << 
            psRbDihedralParameter1->_pSysStream[0][i].z << " " << 
            psRbDihedralParameter1->_pSysStream[0][i].w << " " << 
            psRbDihedralParameter2->_pSysStream[0][i].x << " " << 
            psRbDihedralParameter2->_pSysStream[0][i].y << 
            endl;
#endif
    }
    psRbDihedralID1->Upload();
    psRbDihedralID2->Upload();
    psRbDihedralParameter1->Upload();
    psRbDihedralParameter2->Upload();
}

extern "C"
int gpuReadLJ14Parameters(gpuContext gpu, char* fname)
{
    ifstream infile(fname);
    
    if (!infile.fail())
    {
        char buff[1024];
        float epsfac = 0.0f;
        float fudge = 0.0f;
        int LJ14s;
        infile >> LJ14s;
        infile.get(buff, 61);
        // cout << buff << endl;
        infile >> epsfac;
        infile.get(buff, 8);
        infile >> fudge;
        infile.getline(buff, 512);
        // cout << buff << endl;
        
        vector<int> atom1(LJ14s);
        vector<int> atom2(LJ14s);
        vector<float> c6(LJ14s);
        vector<float> c12(LJ14s);
        vector<float> q1(LJ14s);
        vector<float> q2(LJ14s);
 
        for (int i = 0; i < LJ14s; i++)
        {
            int junk;
            infile >> 
                junk >> 
                atom1[i] >> 
                atom2[i] >> 
                c6[i] >> 
                c12[i] >>
                q1[i] >>
                q2[i];
        }
        gpuSetLJ14Parameters(gpu, epsfac, fudge, atom1, atom2, c6, c12, q1, q2);
        return LJ14s;
    }
    else
    {
        cout << "Error opening Lennard-Jones 1-4 parameter file " << fname << endl;
        exit(-1);
    }
    return 0;
}

extern "C"
void gpuSetLJ14Parameters(gpuContext gpu, float epsfac, float fudge, const vector<int>& atom1, const vector<int>& atom2,
        const vector<float>& c6, const vector<float>& c12, const vector<float>& q1, const vector<float>& q2)
{
    int LJ14s = atom1.size();
    float scale = epsfac * fudge;

    gpu->sim.LJ14s                              = LJ14s;
    CUDAStream<int4>* psLJ14ID                  = new CUDAStream<int4>(LJ14s, 1);
    gpu->psLJ14ID                               = psLJ14ID;
    gpu->sim.pLJ14ID                            = psLJ14ID->_pDevStream[0];
    CUDAStream<float4>* psLJ14Parameter         = new CUDAStream<float4>(LJ14s, 1);
    gpu->psLJ14Parameter                        = psLJ14Parameter;
    gpu->sim.pLJ14Parameter                     = psLJ14Parameter->_pDevStream[0];

    for (int i = 0; i < LJ14s; i++)
    {
        psLJ14ID->_pSysStream[0][i].x = atom1[i];
        psLJ14ID->_pSysStream[0][i].y = atom2[i];
        psLJ14ID->_pSysStream[0][i].z = gpu->pOutputBufferCounter[psLJ14ID->_pSysStream[0][i].x]++;
        psLJ14ID->_pSysStream[0][i].w = gpu->pOutputBufferCounter[psLJ14ID->_pSysStream[0][i].y]++;
        float p0, p1, p2;
        if (c12[i] == 0.0f)
        {
            p0 = 0.0f;
            p1 = 1.0f;
        }
        else
        {
            p0 = c6[i] * c6[i] / c12[i];
            p1 = pow(c12[i] / c6[i], 1.0f / 6.0f);
        }
        p2 = scale * q1[i] * q2[i];
        psLJ14Parameter->_pSysStream[0][i].x = p0;
        psLJ14Parameter->_pSysStream[0][i].y = p1;
        psLJ14Parameter->_pSysStream[0][i].z = p2;
    }
#if (DUMP_PARAMETERS == 1)
        cout << 
            i << " " <<
            psLJ14ID->_pSysStream[0][i].x << " " << 
            psLJ14ID->_pSysStream[0][i].y << " " << 
            psLJ14ID->_pSysStream[0][i].z << " " << 
            psLJ14ID->_pSysStream[0][i].w << " " << 
            psLJ14Parameter->_pSysStream[0][i].x << " " << 
            psLJ14Parameter->_pSysStream[0][i].y << " " <<
            psLJ14Parameter->_pSysStream[0][i].z << " " << 
            p0 << " " << 
            p1 << " " << 
            p2 << " " << 
            endl;
#endif
    psLJ14ID->Upload();
    psLJ14Parameter->Upload();
}

extern "C"
float gpuGetAtomicRadius(gpuContext gpu, string s)
{
    for (int i = 0; i < gpu->gAtomTypes; i++)
    {
        if (s == gpu->gpAtomTable[i].name)
        {
            return gpu->gpAtomTable[i].r;
        }
    }
    
    return 0.0f;
}

extern "C"
unsigned char gpuGetAtomicSymbol(gpuContext gpu, string s)
{
    for (int i = 0; i < gpu->gAtomTypes; i++)
    {
        if (s == gpu->gpAtomTable[i].name)
        {
            return gpu->gpAtomTable[i].symbol;
        }
    }
    
    return ' ';
}

extern "C"
int gpuReadAtomicParameters(gpuContext gpu, char* fname)
{
    gpu->gAtomTypes = 0;
    if (gpu->gpAtomTable)
        delete[] gpu->gpAtomTable;
    
    // Read file once to count atom types
    ifstream infile(fname);
    
    if (!infile.fail())
    {
        char buff[1024];
        int skips = 0;
        bool skipflag = true;
        while (infile.getline(buff, 512))
        {
            if (buff[0] == ' ')
            {
                skipflag = false;
                gpu->gAtomTypes++;
            }
            else if (skipflag)
                skips++;
        }
        infile.close();
        
        gpu->gpAtomTable = new gpuAtomType[gpu->gAtomTypes];
        ifstream infile1(fname);
        for (int i = 0; i < skips; i++)
        {
            infile1.getline(buff, 512);
        }
        for (int i = 0; i < gpu->gAtomTypes; i++)
        {
            infile1 >> gpu->gpAtomTable[i].name >> gpu->gpAtomTable[i].r;
            infile1.getline(buff, 512);
        
            // Determine symbol
            if (gpu->gpAtomTable[i].r < 1.3f)
                gpu->gpAtomTable[i].symbol = 'H';
            else if (gpu->gpAtomTable[i].r < 1.6f)
                gpu->gpAtomTable[i].symbol = 'O';
            else if (gpu->gpAtomTable[i].r < 1.7f)
                gpu->gpAtomTable[i].symbol = 'N';
            else
                gpu->gpAtomTable[i].symbol = 'C';

#if (DUMP_PARAMETERS == 1)            
            cout << i << " " << gpu->gpAtomTable[i].name << " " << gpu->gpAtomTable[i].symbol << " " << gpu->gpAtomTable[i].r << endl; 
#endif
        }
        return gpu->gAtomTypes;
    }
    else
    {
        cout << "Error opening atom parameter file " << fname << endl;
        exit(-1);
    }
    return 0;   

}

extern "C"
int gpuReadCoulombParameters(gpuContext gpu, char* fname)
{
    ifstream infile(fname);
    
    if (!infile.fail())
    {
        char buff[1024];
        unsigned int coulombs;
        float fudge = 0.0f;
        float epsfac = 1.0f;
        infile >> coulombs;
        infile.get(buff, 9);
        infile >> epsfac;
        infile.get(buff, 8);
        infile >> fudge;
        infile.getline(buff, 512);
        vector<int> atom(coulombs);
        vector<float> c6(coulombs);
        vector<float> c12(coulombs);
        vector<float> q(coulombs);
        vector<float> radius(coulombs);
        vector<float> scale(coulombs);
        vector<char> symbol(coulombs);
        vector<vector<int> > exclusions(coulombs);
        unsigned int total_exclusions = 0;

        for (unsigned int i = 0; i < coulombs; i++)
        {
            int junk, numExclusions;
            char atype[512];
            infile >> 
                junk >> 
                c6[i] >>
                c12[i] >>
                q[i] >>
                atype >>
                scale[i] >>
                numExclusions;
                radius[i] = gpuGetAtomicRadius(gpu, atype);
                symbol[i] = gpuGetAtomicSymbol(gpu, atype);
                for (int j = 0; j < numExclusions; j++)
                {
                    int exclusion;
                    infile >> exclusion;
                    exclusions[i].push_back(exclusion);
                }
        }
        cout << total_exclusions << " total exclusions.\n";
        gpuSetCoulombParameters(gpu, epsfac, atom, c6, c12, q, symbol, exclusions, NO_CUTOFF);
        gpuSetObcParameters(gpu, defaultInnerDielectric, defaultSolventDielectric, atom, radius, scale);
        return coulombs;
    }
    else
    {
        cout << "Error opening Coulomb parameter file " << fname << endl;
        exit(-1);
    }
    return 0;
}

extern "C"
void gpuSetCoulombParameters(gpuContext gpu, float epsfac, const vector<int>& atom, const vector<float>& c6, const vector<float>& c12, const vector<float>& q,
        const vector<char>& symbol, const vector<vector<int> >& exclusions, CudaNonbondedMethod method)
{
    unsigned int coulombs = atom.size();
    gpu->sim.epsfac = epsfac;
    gpu->sim.nonbondedMethod = method;
    gpu->exclusions = exclusions;

    for (unsigned int i = 0; i < coulombs; i++)
    {
            float p0 = q[i];
            float p1 = 0.5f, p2 = 0.0f;               
            if ((c6[i] > 0.0f) && (c12[i] > 0.0f))
            {
                p1 = 0.5f * pow(c12[i] / c6[i], 1.0f / 6.0f);
                p2 = c6[i] * sqrt(1.0f / c12[i]);
            }
            if (symbol.size() > 0)
                gpu->pAtomSymbol[i] = symbol[i];
            gpu->psPosq4->_pSysStream[0][i].w = p0;
            gpu->psSigEps2->_pSysStream[0][i].x = p1;
            gpu->psSigEps2->_pSysStream[0][i].y = p2;
    }

    // Dummy out extra atom data
    for (unsigned int i = coulombs; i < gpu->sim.paddedNumberOfAtoms; i++)
    {
        gpu->psPosq4->_pSysStream[0][i].x       = 100000.0f + i * 10.0f;
        gpu->psPosq4->_pSysStream[0][i].y       = 100000.0f + i * 10.0f;
        gpu->psPosq4->_pSysStream[0][i].z       = 100000.0f + i * 10.0f;
        gpu->psPosq4->_pSysStream[0][i].w       = 0.0f;
        gpu->psSigEps2->_pSysStream[0][i].x     = 0.0f;
        gpu->psSigEps2->_pSysStream[0][i].y     = 0.0f;   
    }

    gpu->psPosq4->Upload();
    gpu->psSigEps2->Upload();
}

extern "C"
void gpuSetNonbondedCutoff(gpuContext gpu, float cutoffDistance, float solventDielectric)
{
    gpu->sim.nonbondedCutoffSqr = cutoffDistance*cutoffDistance;
    gpu->sim.reactionFieldK = pow(cutoffDistance, -3.0f)*(solventDielectric-1.0f)/(2.0f*solventDielectric+1.0f);
}

extern "C"
void gpuSetPeriodicBoxSize(gpuContext gpu, float xsize, float ysize, float zsize)
{
    gpu->sim.periodicBoxSizeX = xsize;
    gpu->sim.periodicBoxSizeY = ysize;
    gpu->sim.periodicBoxSizeZ = zsize;
}

extern "C"
void gpuSetObcParameters(gpuContext gpu, float innerDielectric, float solventDielectric, const vector<int>& atom, const vector<float>& radius, const vector<float>& scale)
{
    unsigned int atoms = atom.size();

    gpu->bIncludeGBSA = true;
    for (unsigned int i = 0; i < atoms; i++)
    {
            gpu->psObcData->_pSysStream[0][i].x = radius[i] - dielectricOffset;
            gpu->psObcData->_pSysStream[0][i].y = scale[i] * gpu->psObcData->_pSysStream[0][i].x;

#if (DUMP_PARAMETERS == 1)
        cout << 
            i << " " << 
            gpu->psObcData->_pSysStream[0][i].x << " " <<
            gpu->psObcData->_pSysStream[0][i].y;
#endif
    }

    // Dummy out extra atom data
    for (unsigned int i = atoms; i < gpu->sim.paddedNumberOfAtoms; i++)
    {
        gpu->psBornRadii->_pSysStream[0][i]     = 0.2f;
        gpu->psObcData->_pSysStream[0][i].x     = 0.01f;
        gpu->psObcData->_pSysStream[0][i].y     = 0.01f;
    }

    gpu->psBornRadii->Upload();
    gpu->psObcData->Upload();
    gpu->sim.preFactor = 2.0f*electricConstant*((1.0f/innerDielectric)-(1.0f/solventDielectric))*gpu->sim.forceConversionFactor;
}

extern "C"
int gpuReadShakeParameters(gpuContext gpu, char* fname)
{
    ifstream infile(fname);
    if (!infile.fail())
    {
        char buff[512];
        int shake_constraints;
        infile >> buff >> shake_constraints;
        infile.getline(buff, 512);
        vector<int> atom1(shake_constraints);
        vector<int> atom2(shake_constraints);
        vector<float> distance(shake_constraints);
        vector<float> invMass1(shake_constraints);
        vector<float> invMass2(shake_constraints);

        for (int i = 0; i < shake_constraints; i++)
        {
            int junk;
            infile >> 
                junk >> 
                atom1[i] >> 
                atom2[i] >> 
                distance[i] >> 
                invMass1[i] >> 
                invMass2[i];
        }
        gpuSetShakeParameters(gpu, atom1, atom2, distance, invMass1, invMass2, 1e-4f);
        return gpu->sim.ShakeConstraints;
    }
    else
    {
        cout << "Error opening Shake parameter file " << fname << endl;
        exit(-1);
    }
    return 0;
}

extern "C"
void gpuSetShakeParameters(gpuContext gpu, const vector<int>& atom1, const vector<int>& atom2, const vector<float>& distance,
        const vector<float>& invMass1, const vector<float>& invMass2, float tolerance)
{
    // Find how many constraints each atom is involved in.
    
    vector<int> constraintCount(gpu->natoms, 0);
    for (int i = 0; i < atom1.size(); i++) {
        constraintCount[atom1[i]]++;
        constraintCount[atom2[i]]++;
    }

    // Identify clusters of three atoms that can be treated with SETTLE.  First, for every
    // atom that might be part of such a cluster, make a list of the two other atoms it is
    // connected to.

    vector<map<int, float> > settleConstraints(gpu->natoms);
    for (int i = 0; i < atom1.size(); i++) {
        if (constraintCount[atom1[i]] == 2 && constraintCount[atom2[i]] == 2) {
            settleConstraints[atom1[i]][atom2[i]] = distance[i];
            settleConstraints[atom2[i]][atom1[i]] = distance[i];
        }
    }

    // Now remove the ones that don't actually form closed loops of three atoms.

    vector<int> settleClusters;
    for (int i = 0; i < settleConstraints.size(); i++) {
        if (settleConstraints[i].size() == 2) {
            int partner1 = settleConstraints[i].begin()->first;
            int partner2 = (++settleConstraints[i].begin())->first;
            if (settleConstraints[partner1].size() != 2 || settleConstraints[partner2].size() != 2 ||
                    settleConstraints[partner1].find(partner2) == settleConstraints[partner1].end())
                settleConstraints[i].clear();
            else if (i < partner1 && i < partner2)
                settleClusters.push_back(i);
        }
        else
            settleConstraints[i].clear();
    }

    // Record the actual SETTLE clusters.

    CUDAStream<int4>* psSettleID          = new CUDAStream<int4>((int) settleClusters.size(), 1);
    gpu->psSettleID                       = psSettleID;
    gpu->sim.pSettleID                    = psSettleID->_pDevStream[0];
    CUDAStream<float2>* psSettleParameter = new CUDAStream<float2>((int) settleClusters.size(), 1);
    gpu->psSettleParameter                = psSettleParameter;
    gpu->sim.pSettleParameter             = psSettleParameter->_pDevStream[0];
    gpu->sim.settleConstraints            = settleClusters.size();
    for (int i = 0; i < settleClusters.size(); i++) {
        int atom1 = settleClusters[i];
        int atom2 = settleConstraints[atom1].begin()->first;
        int atom3 = (++settleConstraints[atom1].begin())->first;
        float dist12 = settleConstraints[atom1].find(atom2)->second;
        float dist13 = settleConstraints[atom1].find(atom3)->second;
        float dist23 = settleConstraints[atom2].find(atom3)->second;
        if (dist12 == dist13) { // atom1 is the central atom
            psSettleID->_pSysData[i].x = atom1;
            psSettleID->_pSysData[i].y = atom2;
            psSettleID->_pSysData[i].z = atom3;
            psSettleParameter->_pSysData[i].x = dist12;
            psSettleParameter->_pSysData[i].y = dist23;
        }
        else if (dist12 == dist23) { // atom2 is the central atom
            psSettleID->_pSysData[i].x = atom2;
            psSettleID->_pSysData[i].y = atom1;
            psSettleID->_pSysData[i].z = atom3;
            psSettleParameter->_pSysData[i].x = dist12;
            psSettleParameter->_pSysData[i].y = dist13;
        }
        else if (dist13 == dist23) { // atom3 is the central atom
            psSettleID->_pSysData[i].x = atom3;
            psSettleID->_pSysData[i].y = atom1;
            psSettleID->_pSysData[i].z = atom2;
            psSettleParameter->_pSysData[i].x = dist13;
            psSettleParameter->_pSysData[i].y = dist12;
        }
        else
            throw OpenMMException("Two of the three distances constrained with SETTLE must be the same.");
    }
    psSettleID->Upload();
    psSettleParameter->Upload();
    gpu->sim.settle_threads_per_block     = (gpu->sim.settleConstraints + gpu->sim.blocks - 1) / gpu->sim.blocks;
    if (gpu->sim.settle_threads_per_block > gpu->sim.max_shake_threads_per_block)
        gpu->sim.settle_threads_per_block = gpu->sim.max_shake_threads_per_block;
    if (gpu->sim.settle_threads_per_block < 1)
        gpu->sim.settle_threads_per_block = 1;

    // Find clusters consisting of a central atom with up to three peripheral atoms.
    
    map<int, ShakeCluster> clusters;
    for (int i = 0; i < atom1.size(); i++) {
        if (settleConstraints[atom1[i]].size() == 2)
            continue; // This is being taken care of with SETTLE.

        // Determine which is the central atom.
        
        bool firstIsCentral;
        if (constraintCount[atom1[i]] > 1)
            firstIsCentral = true;
        else if (constraintCount[atom2[i]] > 1)
            firstIsCentral = false;
        else if (atom1[i] < atom2[i])
            firstIsCentral = true;
        else
            firstIsCentral = false;
        int centralID, peripheralID;
        float centralInvMass, peripheralInvMass;
        if (firstIsCentral) {
            centralID = atom1[i];
            peripheralID = atom2[i];
            centralInvMass = invMass1[i];
            peripheralInvMass = invMass2[i];
        }
        else {
            centralID = atom2[i];
            peripheralID = atom1[i];
            centralInvMass = invMass2[i];
            peripheralInvMass = invMass1[i];
        }
        if (constraintCount[peripheralID] != 1)
            throw OpenMMException("Only bonds to hydrogens may be constrained");
        
        // Add it to the cluster.
        
        if (clusters.find(centralID) == clusters.end()) {
            clusters[centralID] = ShakeCluster(centralID, centralInvMass);
        }
        clusters[centralID].addAtom(peripheralID, distance[i], peripheralInvMass);
    }
    
    // Fill in the Cuda streams.
    
    CUDAStream<int4>* psShakeID             = new CUDAStream<int4>((int) clusters.size(), 1);
    gpu->psShakeID                          = psShakeID;
    gpu->sim.pShakeID                       = psShakeID->_pDevStream[0]; 
    CUDAStream<float4>* psShakeParameter    = new CUDAStream<float4>((int) clusters.size(), 1);
    gpu->psShakeParameter                   = psShakeParameter;
    gpu->sim.pShakeParameter                = psShakeParameter->_pDevStream[0];
    gpu->sim.ShakeConstraints               = clusters.size();
    int index = 0;
    for (map<int, ShakeCluster>::const_iterator iter = clusters.begin(); iter != clusters.end(); ++iter) {
        const ShakeCluster& cluster = iter->second;
        psShakeID->_pSysStream[0][index].x = cluster.centralID;
        psShakeID->_pSysStream[0][index].y = cluster.peripheralID[0];
        psShakeID->_pSysStream[0][index].z = cluster.size > 1 ? cluster.peripheralID[1] : -1;
        psShakeID->_pSysStream[0][index].w = cluster.size > 2 ? cluster.peripheralID[2] : -1;
        psShakeParameter->_pSysStream[0][index].x = cluster.centralInvMass;
        psShakeParameter->_pSysStream[0][index].y = 0.5f/(cluster.centralInvMass+cluster.peripheralInvMass);
        psShakeParameter->_pSysStream[0][index].z = cluster.distance*cluster.distance;
        psShakeParameter->_pSysStream[0][index].w = cluster.peripheralInvMass;
        ++index;
    }
    psShakeID->Upload();
    psShakeParameter->Upload();
    gpu->sim.shakeTolerance = tolerance;

    gpu->sim.shake_threads_per_block     = (gpu->sim.ShakeConstraints + gpu->sim.blocks - 1) / gpu->sim.blocks; 
    if (gpu->sim.shake_threads_per_block > gpu->sim.max_shake_threads_per_block)
        gpu->sim.shake_threads_per_block = gpu->sim.max_shake_threads_per_block;
    if (gpu->sim.shake_threads_per_block < 1)
        gpu->sim.shake_threads_per_block = 1;

#ifdef DeltaShake

    // count number of atoms w/o constraint

    int count = 0;
    for (int i = 0; i < gpu->natoms; i++)
       if (constraintCount[i] == 0)
          count++;

    // Allocate NonShake parameters

    gpu->sim.NonShakeConstraints                  = count;
    if( count || true ){

       CUDAStream<int>* psNonShakeID              = new CUDAStream<int>(count, 1);
       gpu->psNonShakeID                          = psNonShakeID;
       gpu->sim.pNonShakeID                       = psNonShakeID->_pDevStream[0];

       gpu->sim.nonshake_threads_per_block        = (count + gpu->sim.blocks - 1) / gpu->sim.blocks;

       if (gpu->sim.nonshake_threads_per_block > gpu->sim.max_shake_threads_per_block)
           gpu->sim.nonshake_threads_per_block = gpu->sim.max_shake_threads_per_block;

       if (gpu->sim.nonshake_threads_per_block < 1)
               gpu->sim.nonshake_threads_per_block = 1;

       // load indices

       count = 0;
       for (int i = 0; i < gpu->natoms; i++){
          if (constraintCount[i] == 0){
             psNonShakeID->_pSysStream[0][count++] = i;
          }
       }
       psNonShakeID->Upload();

    } else {
       gpu->sim.nonshake_threads_per_block           = 0;
    }
#endif
}

extern "C"
int gpuAllocateInitialBuffers(gpuContext gpu)
{
    gpu->sim.atoms                      = gpu->natoms;
    gpu->sim.paddedNumberOfAtoms        = ((gpu->sim.atoms + GRID - 1) >> GRIDBITS) << GRIDBITS;
    gpu->sim.degreesOfFreedom           = 3 * gpu->sim.atoms - 6;
    gpu->gpAtomTable                    = NULL;
    gpu->gAtomTypes                     = 0;
    gpu->psPosq4                        = new CUDAStream<float4>(gpu->sim.paddedNumberOfAtoms, 1);
    gpu->sim.stride                     = gpu->psPosq4->_stride;
    gpu->sim.stride2                    = gpu->sim.stride * 2;
    gpu->sim.stride3                    = gpu->sim.stride * 3;
    gpu->sim.stride4                    = gpu->sim.stride * 4;
    gpu->sim.pPosq                      = gpu->psPosq4->_pDevStream[0];
    gpu->sim.stride                     = gpu->psPosq4->_stride;
    gpu->sim.stride2                    = 2 * gpu->sim.stride;
    gpu->sim.stride3                    = 3 * gpu->sim.stride;
    gpu->sim.stride4                    = 4 * gpu->sim.stride;
    gpu->sim.exclusionStride            = gpu->sim.stride / GRID;
    gpu->psPosqP4                       = new CUDAStream<float4>(gpu->sim.paddedNumberOfAtoms, 1);
    gpu->sim.pPosqP                     = gpu->psPosqP4->_pDevStream[0];
    gpu->psOldPosq4                     = new CUDAStream<float4>(gpu->sim.paddedNumberOfAtoms, 1);
    gpu->sim.pOldPosq                   = gpu->psOldPosq4->_pDevStream[0];
    gpu->psVelm4                        = new CUDAStream<float4>(gpu->sim.paddedNumberOfAtoms, 1);
    gpu->sim.pVelm4                     = gpu->psVelm4->_pDevStream[0];
    gpu->psvVector4                     = new CUDAStream<float4>(gpu->sim.paddedNumberOfAtoms, 1);
    gpu->sim.pvVector4                  = gpu->psvVector4->_pDevStream[0];
    gpu->psxVector4                     = new CUDAStream<float4>(gpu->sim.paddedNumberOfAtoms, 1);
    gpu->sim.pxVector4                  = gpu->psxVector4->_pDevStream[0];
    gpu->psBornRadii                    = new CUDAStream<float>(gpu->sim.paddedNumberOfAtoms, 1);
    gpu->sim.pBornRadii                 = gpu->psBornRadii->_pDevStream[0];
    gpu->psObcChain                     = new CUDAStream<float>(gpu->sim.paddedNumberOfAtoms, 1);
    gpu->sim.pObcChain                  = gpu->psObcChain->_pDevStream[0];
    gpu->psSigEps2                      = new CUDAStream<float2>(gpu->sim.paddedNumberOfAtoms, 1);
    gpu->sim.pAttr                      = gpu->psSigEps2->_pDevStream[0];
    gpu->psObcData                      = new CUDAStream<float2>(gpu->sim.paddedNumberOfAtoms, 1);
    gpu->sim.pObcData                   = gpu->psObcData->_pDevStream[0];
    gpu->pAtomSymbol                    = new unsigned char[gpu->natoms];
    gpu->psAtomIndex                    = new CUDAStream<int>(gpu->sim.paddedNumberOfAtoms, 1);
    gpu->sim.pAtomIndex                 = gpu->psAtomIndex->_pDevStream[0];
    for (int i = 0; i < (int) gpu->sim.paddedNumberOfAtoms; i++)
        gpu->psAtomIndex->_pSysStream[0][i] = i;
    gpu->psAtomIndex->Upload();
    // Determine randoms
    gpu->seed                           = (unsigned long)time(NULL) & 0x000fffff;
    gpu->sim.randomFrames               = 95;
    gpu->sim.randomIterations           = gpu->sim.randomFrames;
    gpu->sim.randoms                    = gpu->sim.randomFrames * gpu->sim.paddedNumberOfAtoms - 5 * GRID;
    gpu->sim.totalRandoms               = gpu->sim.randoms + gpu->sim.paddedNumberOfAtoms;
    gpu->sim.totalRandomsTimesTwo       = gpu->sim.totalRandoms * 2;
    gpu->psRandom4                      = new CUDAStream<float4>(gpu->sim.totalRandomsTimesTwo, 1);
    gpu->psRandom2                      = new CUDAStream<float2>(gpu->sim.totalRandomsTimesTwo, 1);
    gpu->psRandomPosition               = new CUDAStream<int>(gpu->sim.blocks, 1);
    gpu->psRandomSeed                   = new CUDAStream<uint4>(gpu->sim.blocks * gpu->sim.random_threads_per_block, 1);
    gpu->sim.pRandom4a                  = gpu->psRandom4->_pDevStream[0];
    gpu->sim.pRandom2a                  = gpu->psRandom2->_pDevStream[0];
    gpu->sim.pRandom4b                  = gpu->psRandom4->_pDevStream[0] + gpu->sim.totalRandoms;
    gpu->sim.pRandom2b                  = gpu->psRandom2->_pDevStream[0] + gpu->sim.totalRandoms;
    gpu->sim.pRandomPosition            = gpu->psRandomPosition->_pDevStream[0];
    gpu->sim.pRandomSeed                = gpu->psRandomSeed->_pDevStream[0];
    for (int i = 0; i < (int) gpu->sim.blocks; i++)
    {
        gpu->psRandomPosition->_pSysStream[0][i] = 0;
    }
    int seed = gpu->seed | ((gpu->seed ^ 0xffffffff) << 16);
    srand(seed);
    for (int i = 0; i < (int) (gpu->sim.blocks * gpu->sim.random_threads_per_block); i++)
    {
        gpu->psRandomSeed->_pSysStream[0][i].x = rand();
        gpu->psRandomSeed->_pSysStream[0][i].y = rand();
        gpu->psRandomSeed->_pSysStream[0][i].z = rand();
        gpu->psRandomSeed->_pSysStream[0][i].w = rand();
    }

    float randomValue = 0.0f;
    for (int i = 0; i < (int) gpu->sim.totalRandomsTimesTwo; i++)
    {
        gpu->psRandom4->_pSysStream[0][i].x         = randomValue;
        gpu->psRandom4->_pSysStream[0][i].y         = randomValue;
        gpu->psRandom4->_pSysStream[0][i].z         = randomValue;
        gpu->psRandom4->_pSysStream[0][i].w         = randomValue;
        gpu->psRandom2->_pSysStream[0][i].x         = randomValue;
        gpu->psRandom2->_pSysStream[0][i].y         = randomValue;
    }

    gpu->psRandomSeed->Upload();
    gpu->psRandom4->Upload();
    gpu->psRandom2->Upload();
    gpu->psRandomPosition->Upload();

    // Allocate and clear linear momentum buffer
    gpu->psLinearMomentum = new CUDAStream<float4>(gpu->sim.blocks, 1);
    gpu->sim.pLinearMomentum = gpu->psLinearMomentum->_pDevStream[0];
    for (int i = 0; i < (int) gpu->sim.blocks; i++)
    {
        gpu->psLinearMomentum->_pSysStream[0][i].x = 0.0f;
        gpu->psLinearMomentum->_pSysStream[0][i].y = 0.0f;
        gpu->psLinearMomentum->_pSysStream[0][i].z = 0.0f;
        gpu->psLinearMomentum->_pSysStream[0][i].w = 0.0f;
    }
    gpu->psLinearMomentum->Upload();

    return 1;
}

extern "C"
void gpuReadCoordinates(gpuContext gpu, char* fname)
{
    ifstream infile(fname);
    gpu->natoms = 0;
    char buff[512];
    infile >> buff >> gpu->natoms;
    infile.getline(buff, 511);
    float totalMass = 0.0f;

    gpuAllocateInitialBuffers(gpu);
    
    for (int i = 0; i < gpu->natoms; i++)
    {
        int junk;
        infile >> junk >> 
            gpu->psPosq4->_pSysStream[0][i].x >> 
            gpu->psPosq4->_pSysStream[0][i].y >> 
            gpu->psPosq4->_pSysStream[0][i].z >>
            gpu->psPosq4->_pSysStream[0][i].w >>
            gpu->psVelm4->_pSysStream[0][i].x >> 
            gpu->psVelm4->_pSysStream[0][i].y >> 
            gpu->psVelm4->_pSysStream[0][i].z >>
            gpu->psVelm4->_pSysStream[0][i].w;
        gpu->psxVector4->_pSysStream[0][i].x = 0.0f;
        gpu->psxVector4->_pSysStream[0][i].y = 0.0f;
        gpu->psxVector4->_pSysStream[0][i].z = 0.0f;
        gpu->psxVector4->_pSysStream[0][i].w = 0.0f;

        // Accumulate mass
        totalMass += 1.0f / gpu->psVelm4->_pSysStream[0][i].w;
    }
    
    gpu->sim.inverseTotalMass = 1.0f / totalMass;
    gpu->psPosq4->Upload();
    gpu->psVelm4->Upload();
    gpu->psxVector4->Upload();
}

extern "C"
void gpuSetPositions(gpuContext gpu, const vector<float>& x, const vector<float>& y, const vector<float>& z)
{
    for (int i = 0; i < gpu->natoms; i++)
    {
        gpu->psPosq4->_pSysStream[0][i].x = x[i];
        gpu->psPosq4->_pSysStream[0][i].y = y[i];
        gpu->psPosq4->_pSysStream[0][i].z = z[i];
    }
    gpu->psPosq4->Upload();

	 // set flag to recalculate Born radii

	 gpu->bRecalculateBornRadii = true;
} 

extern "C"
void gpuSetVelocities(gpuContext gpu, const vector<float>& x, const vector<float>& y, const vector<float>& z)
{
    for (int i = 0; i < gpu->natoms; i++)
    {
        gpu->psVelm4->_pSysStream[0][i].x = x[i];
        gpu->psVelm4->_pSysStream[0][i].y = y[i];
        gpu->psVelm4->_pSysStream[0][i].z = z[i];
    }
    gpu->psVelm4->Upload();
} 

extern "C"
void gpuSetMass(gpuContext gpu, const vector<float>& mass)
{
    float totalMass = 0.0f;
    for (int i = 0; i < gpu->natoms; i++)
    {
        gpu->psVelm4->_pSysStream[0][i].w = 1.0f/mass[i];
        totalMass += mass[i];
    }
    gpu->sim.inverseTotalMass = 1.0f / totalMass;
    gpu->psVelm4->Upload();
} 

extern "C"
void gpuInitializeRandoms(gpuContext gpu)
{
    for (int i = 0; i < (int) gpu->sim.blocks; i++)
    {
        gpu->psRandomPosition->_pSysStream[0][i] = 0;
    }
    int seed = gpu->seed | ((gpu->seed ^ 0xffffffff) << 16);
    srand(seed);
    for (int i = 0; i < (int) (gpu->sim.blocks * gpu->sim.random_threads_per_block); i++)
    {
        gpu->psRandomSeed->_pSysStream[0][i].x = rand();
        gpu->psRandomSeed->_pSysStream[0][i].y = rand();
        gpu->psRandomSeed->_pSysStream[0][i].z = rand();
        gpu->psRandomSeed->_pSysStream[0][i].w = rand();
    }
    gpu->psRandomPosition->Upload();
    gpu->psRandomSeed->Upload();
    gpuSetConstants(gpu);
    kGenerateRandoms(gpu);
    return;
}

extern "C"
bool gpuIsAvailable()
{
    int deviceCount;
    cudaGetDeviceCount(&deviceCount);
    return (deviceCount > 0);
}

extern "C"
void* gpuInitFromFile(char* fname) 
{
    ifstream infile(fname);
    int numAtoms = 0;
    char buff[512];
    infile >> buff >> numAtoms;
    gpuContext gpu = (gpuContext) gpuInit(numAtoms);
    vector<float> x(numAtoms), y(numAtoms), z(numAtoms), charge(numAtoms), vx(numAtoms), vy(numAtoms), vz(numAtoms), mass(numAtoms);
    infile.getline(buff, 511);
    float totalMass = 0.0f;
    for (int i = 0; i < gpu->natoms; i++)
    {
        int junk;
        infile >> junk >> 
            x[i] >> 
            y[i] >> 
            z[i] >> 
            charge[i] >> 
            vx[i] >> 
            vy[i] >> 
            vz[i] >> 
            mass[i];
        mass[i] = 1.0f/mass[i];
    }
    gpuSetPositions(gpu, x, y, z);
    gpuSetVelocities(gpu, vx, vy, vz);
    gpuSetMass(gpu, mass);
    return (void*)gpu;
}

extern "C"
void* gpuInit(int numAtoms)
{
    gpuContext gpu = new _gpuContext;
    int LRFSize = 0;
    int SMCount = 0;
    int SMMajor = 0;
    int SMMinor = 0;

    // Get adapter
    unsigned int device = 0;
    char * pAdapter;
    pAdapter = getenv ("NV_FAH_DEVICE");
    if (pAdapter != NULL)
    {
        sscanf(pAdapter, "%d", &device);
    }
//    cudaError_t status = cudaSetDevice(device);
//    RTERROR(status, "Error setting CUDA device")

    // Determine which core to run on
#if 0
    SYSTEM_INFO info;
    GetSystemInfo(&info);
    unsigned int cores = info.dwNumberOfProcessors;
    if (cores > 1)
    {
        HANDLE hproc = GetCurrentProcess();
        unsigned int core = (cores - 1) - (device % (cores - 1)); 
        unsigned int mask = 1 << core;
        SetProcessAffinityMask(hproc, mask);
    }
#endif

    // Determine kernel call configuration
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, 0);

    // Determine SM version
    if (deviceProp.major == 1)
    {
        switch (deviceProp.minor)
        {
        case 0:
        case 1:
            gpu->sm_version = SM_10;
            gpu->sim.workUnitsPerSM = G8X_NONBOND_WORKUNITS_PER_SM;
            break;

        default:
            gpu->sm_version = SM_12;
            gpu->sim.workUnitsPerSM = GT2XX_NONBOND_WORKUNITS_PER_SM;
            break;
        }
    }

    gpu->sim.nonbond_blocks = deviceProp.multiProcessorCount;
    gpu->sim.bornForce2_blocks = deviceProp.multiProcessorCount;
    gpu->sim.blocks = deviceProp.multiProcessorCount;
    if (deviceProp.regsPerBlock == 8192)
    {
        gpu->sim.nonbond_threads_per_block          = G8X_NONBOND_THREADS_PER_BLOCK;
        gpu->sim.bornForce2_threads_per_block       = G8X_BORNFORCE2_THREADS_PER_BLOCK;
        gpu->sim.max_shake_threads_per_block        = G8X_SHAKE_THREADS_PER_BLOCK;
        gpu->sim.max_update_threads_per_block       = G8X_UPDATE_THREADS_PER_BLOCK;
        gpu->sim.max_localForces_threads_per_block  = G8X_LOCALFORCES_THREADS_PER_BLOCK;
        gpu->sim.threads_per_block                  = G8X_THREADS_PER_BLOCK;
        gpu->sim.random_threads_per_block           = G8X_RANDOM_THREADS_PER_BLOCK;
    }
    else
    {
        gpu->sim.nonbond_threads_per_block          = GT2XX_NONBOND_THREADS_PER_BLOCK;
        gpu->sim.bornForce2_threads_per_block       = GT2XX_BORNFORCE2_THREADS_PER_BLOCK;
        gpu->sim.max_shake_threads_per_block        = GT2XX_SHAKE_THREADS_PER_BLOCK;
        gpu->sim.max_update_threads_per_block       = GT2XX_UPDATE_THREADS_PER_BLOCK;
        gpu->sim.max_localForces_threads_per_block  = GT2XX_LOCALFORCES_THREADS_PER_BLOCK;
        gpu->sim.threads_per_block                  = GT2XX_NONBOND_THREADS_PER_BLOCK;
        gpu->sim.random_threads_per_block           = GT2XX_RANDOM_THREADS_PER_BLOCK;
    }
    gpu->sim.shake_threads_per_block                = gpu->sim.max_shake_threads_per_block;
    gpu->sim.localForces_threads_per_block          = gpu->sim.max_localForces_threads_per_block;

    gpu->natoms = numAtoms;
    gpuAllocateInitialBuffers(gpu);
    for (int i = 0; i < gpu->natoms; i++)
    {
        gpu->psxVector4->_pSysStream[0][i].x = 0.0f;
        gpu->psxVector4->_pSysStream[0][i].y = 0.0f;
        gpu->psxVector4->_pSysStream[0][i].z = 0.0f;
        gpu->psxVector4->_pSysStream[0][i].w = 0.0f;
    }
    gpu->psxVector4->Upload();

    gpu->iterations = 0;
    gpu->sim.update_threads_per_block               = (gpu->natoms + gpu->sim.blocks - 1) / gpu->sim.blocks;
    if (gpu->sim.update_threads_per_block > gpu->sim.max_update_threads_per_block)
        gpu->sim.update_threads_per_block = gpu->sim.max_update_threads_per_block;
    if (gpu->sim.update_threads_per_block < 1)
            gpu->sim.update_threads_per_block = 1;
    gpu->sim.bf_reduce_threads_per_block = gpu->sim.update_threads_per_block;
    gpu->sim.bsf_reduce_threads_per_block = (gpu->sim.stride4 + gpu->natoms + gpu->sim.blocks - 1) / gpu->sim.blocks;
    gpu->sim.bsf_reduce_threads_per_block = ((gpu->sim.bsf_reduce_threads_per_block + (GRID - 1)) / GRID) * GRID;
    if (gpu->sim.bsf_reduce_threads_per_block > gpu->sim.threads_per_block)
        gpu->sim.bsf_reduce_threads_per_block = gpu->sim.threads_per_block;
    if (gpu->sim.bsf_reduce_threads_per_block < 1)
        gpu->sim.bsf_reduce_threads_per_block = 1;

    // Initialize constants to reasonable values
    gpu->sim.probeRadius            = probeRadius;
    gpu->sim.surfaceAreaFactor      = surfaceAreaFactor;
    gpu->sim.electricConstant       = electricConstant;
    gpu->sim.nonbondedMethod        = NO_CUTOFF;
    gpu->sim.nonbondedCutoffSqr     = 0.0f;

    gpu->sim.bigFloat               = 99999999.0f;
    gpu->sim.forceConversionFactor  = forceConversionFactor;
    gpu->sim.preFactor              = 2.0f*electricConstant*((1.0f/defaultInnerDielectric)-(1.0f/defaultSolventDielectric))*gpu->sim.forceConversionFactor;
    gpu->sim.dielectricOffset       = dielectricOffset;
    gpu->sim.alphaOBC               = alphaOBC;
    gpu->sim.betaOBC                = betaOBC;
    gpu->sim.gammaOBC               = gammaOBC;
    gpuSetIntegrationParameters(gpu, 1.0f, 2.0e-3f, 300.0f);
    gpu->sim.maxShakeIterations     = 15;
    gpu->sim.shakeTolerance         = 1.0e-04f * 2.0f;
    gpu->sim.InvMassJ               = 9.920635e-001f;
    gpu->grid                       = GRID;
    gpu->bCalculateCM               = false;
    gpu->bRemoveCM                  = false;
    gpu->bRecalculateBornRadii      = true;
    gpu->bIncludeGBSA               = false;
    gpuInitializeRandoms(gpu);

    // To be determined later
    gpu->psLJ14ID                   = NULL;
    gpu->psForce4                   = NULL;
    gpu->sim.pForce4                = NULL;
    gpu->sim.pForce4a               = NULL;
    gpu->sim.pForce4b               = NULL;
    gpu->psBornForce                = NULL;
    gpu->sim.pBornForce             = NULL;
    gpu->psBornSum                  = NULL;
    gpu->sim.pBornSum               = NULL;
    gpu->psBondID                   = NULL;
    gpu->psBondParameter            = NULL;
    gpu->psBondAngleID1             = NULL;
    gpu->psBondAngleID2             = NULL;
    gpu->psBondAngleParameter       = NULL;
    gpu->psDihedralID1              = NULL;
    gpu->psDihedralID2              = NULL;
    gpu->psDihedralParameter        = NULL;
    gpu->psRbDihedralID1            = NULL;
    gpu->psRbDihedralID2            = NULL;
    gpu->psRbDihedralParameter1     = NULL;
    gpu->psRbDihedralParameter2     = NULL;
    gpu->psLJ14ID                   = NULL;
    gpu->psLJ14Parameter            = NULL;
    gpu->psShakeID                  = NULL;
    gpu->psShakeParameter           = NULL;
    gpu->psSettleID                 = NULL;
    gpu->psSettleParameter          = NULL;
    gpu->psExclusion                = NULL;
    gpu->psWorkUnit                 = NULL;
    gpu->psInteractingWorkUnit      = NULL;
    gpu->psInteractionFlag          = NULL;
    gpu->psInteractionCount         = NULL;
    gpu->psGridBoundingBox          = NULL;
    gpu->psGridCenter               = NULL;


    // Initialize output buffer before reading parameters
    gpu->pOutputBufferCounter       = new unsigned int[gpu->sim.paddedNumberOfAtoms];
    memset(gpu->pOutputBufferCounter, 0, gpu->sim.paddedNumberOfAtoms * sizeof(unsigned int));

    return (void*)gpu;
}

extern "C"
void gpuSetIntegrationParameters(gpuContext gpu, float tau, float deltaT, float temperature) {
    gpu->sim.deltaT                 = deltaT;
    gpu->sim.oneOverDeltaT          = 1.0f/deltaT;
    gpu->sim.tau                    = tau;
    gpu->sim.GDT                    = gpu->sim.deltaT / gpu->sim.tau;
    gpu->sim.EPH                    = exp(0.5f * gpu->sim.GDT);
    gpu->sim.EMH                    = exp(-0.5f * gpu->sim.GDT);
    gpu->sim.EP                     = exp(gpu->sim.GDT);
    gpu->sim.EM                     = exp(-gpu->sim.GDT);
    gpu->sim.OneMinusEM             = 1.0f - gpu->sim.EM;
    gpu->sim.TauOneMinusEM          = gpu->sim.tau * gpu->sim.OneMinusEM;
    if (gpu->sim.GDT >= 0.1f)
    {
        float term1                 = gpu->sim.EPH - 1.0f;
        term1                      *= term1;
        gpu->sim.B                  = gpu->sim.GDT * (gpu->sim.EP - 1.0f) - 4.0f * term1;
        gpu->sim.C                  = gpu->sim.GDT - 3.0f + 4.0f * gpu->sim.EMH - gpu->sim.EM;
        gpu->sim.D                  = 2.0f - gpu->sim.EPH - gpu->sim.EMH;
    }
    else
    {
        float term1                 = 0.5f * gpu->sim.GDT;
        float term2                 = term1 * term1;
        float term4                 = term2 * term2;

        float third                 = 1.0f / 3.0f;
        float o7_9                  = 7.0f / 9.0f;
        float o1_12                 = 1.0f / 12.0f;
        float o17_90                = 17.0f / 90.0f;
        float o7_30                 = 7.0f / 30.0f;
        float o31_1260              = 31.0f / 1260.0f;
        float o_360                 = 1.0f / 360.0f;

        gpu->sim.B                  = term4 * (third + term1 * (third + term1 * (o17_90 + term1 * o7_9)));
        gpu->sim.C                  = term2 * term1 * (2.0f * third + term1 * (-0.5f + term1 * (o7_30 + term1 * (-o1_12 + term1 * o31_1260))));
        gpu->sim.D                  = term2 * (-1.0f + term2 * (-o1_12 - term2 * o_360));   
    }
    gpu->sim.TauDOverEMMinusOne     = gpu->sim.tau * gpu->sim.D / (gpu->sim.EM - 1.0f);
    gpu->sim.DOverTauC              = gpu->sim.D / (gpu->sim.tau * gpu->sim.C);
    gpu->sim.fix1                   = gpu->sim.tau * (gpu->sim.EPH - gpu->sim.EMH);
    gpu->sim.oneOverFix1            = 1.0f / (gpu->sim.tau * (gpu->sim.EPH - gpu->sim.EMH));
    gpu->sim.T                      = temperature;
    gpu->sim.kT                     = BOLTZ * gpu->sim.T;
    gpu->sim.V                      = sqrt(gpu->sim.kT * (1.0f - gpu->sim.EM));
    gpu->sim.X                      = gpu->sim.tau * sqrt(gpu->sim.kT * gpu->sim.C);
    gpu->sim.Yv                     = sqrt(gpu->sim.kT * gpu->sim.B / gpu->sim.C);
    gpu->sim.Yx                     = gpu->sim.tau * sqrt(gpu->sim.kT * gpu->sim.B / (1.0f - gpu->sim.EM));
}

extern "C"
void gpuSetVerletIntegrationParameters(gpuContext gpu, float deltaT) {
    gpu->sim.deltaT                 = deltaT;
    gpu->sim.oneOverDeltaT          = 1.0f/deltaT;
}

extern "C"
void gpuSetBrownianIntegrationParameters(gpuContext gpu, float tau, float deltaT, float temperature) {
    gpu->sim.deltaT                 = deltaT;
    gpu->sim.oneOverDeltaT          = 1.0f/deltaT;
    gpu->sim.tau                    = tau;
    gpu->sim.GDT                    = gpu->sim.deltaT * gpu->sim.tau;
    gpu->sim.T                      = temperature;
    gpu->sim.kT                     = BOLTZ * gpu->sim.T;
    gpu->sim.Yv = gpu->sim.Yx       = sqrt(2.0f*gpu->sim.kT*deltaT*tau);
}

extern "C"
void gpuSetAndersenThermostatParameters(gpuContext gpu, float temperature, float collisionProbability) {
    gpu->sim.T                      = temperature;
    gpu->sim.kT                     = BOLTZ * gpu->sim.T;
    gpu->sim.collisionProbability   = collisionProbability;
    gpu->sim.Yv = gpu->sim.Yx       = 1.0f;
    gpu->sim.V = gpu->sim.X         = 1.0f;
}

extern "C"
void gpuShutDown(gpuContext gpu)
{
    // Delete sysmem pointers
    delete[] gpu->pOutputBufferCounter;
    delete[] gpu->gpAtomTable;
    delete[] gpu->pAtomSymbol;

    // Delete device pointers
    delete gpu->psPosq4;
    delete gpu->psPosqP4;
    delete gpu->psOldPosq4;
    delete gpu->psVelm4;
    delete gpu->psForce4;
    delete gpu->psxVector4;
    delete gpu->psvVector4;
    delete gpu->psSigEps2; 
    delete gpu->psObcData; 
    delete gpu->psObcChain;
    delete gpu->psBornForce;
    delete gpu->psBornRadii;
    delete gpu->psBornSum;
    delete gpu->psBondID;
    delete gpu->psBondParameter;
    delete gpu->psBondAngleID1;
    delete gpu->psBondAngleID2;
    delete gpu->psBondAngleParameter;
    delete gpu->psDihedralID1;
    delete gpu->psDihedralID2;
    delete gpu->psDihedralParameter;
    delete gpu->psRbDihedralID1;
    delete gpu->psRbDihedralID2;
    delete gpu->psRbDihedralParameter1;
    delete gpu->psRbDihedralParameter2;
    delete gpu->psLJ14ID;
    delete gpu->psLJ14Parameter;
    delete gpu->psShakeID;
    delete gpu->psShakeParameter;
    delete gpu->psSettleID;
    delete gpu->psSettleParameter;
    delete gpu->psExclusion;
    delete gpu->psWorkUnit;
    delete gpu->psInteractingWorkUnit;
    delete gpu->psInteractionFlag;
    delete gpu->psInteractionCount;
    delete gpu->psRandom4;
    delete gpu->psRandom2;
    delete gpu->psRandomPosition;    
    delete gpu->psRandomSeed;
    delete gpu->psLinearMomentum;
    delete gpu->psAtomIndex;
    delete gpu->psGridBoundingBox;
    delete gpu->psGridCenter;
    if (gpu->cudpp != 0)
        cudppDestroyPlan(gpu->cudpp);

    // Wrap up
    delete gpu;
    return;
}

extern "C"
int gpuBuildOutputBuffers(gpuContext gpu)
{
    // Select the number of output buffer to use.
    gpu->bOutputBufferPerWarp           = true;
    gpu->sim.nonbondOutputBuffers       = gpu->sim.nonbond_blocks * gpu->sim.nonbond_threads_per_block / GRID;
    if (gpu->sim.nonbondOutputBuffers >= gpu->sim.paddedNumberOfAtoms/GRID)
    {
        // For small systems, it is more efficient to have one output buffer per block of 32 atoms instead of one per warp.
        gpu->bOutputBufferPerWarp           = false;
        gpu->sim.nonbondOutputBuffers       = gpu->sim.paddedNumberOfAtoms / GRID;
    }
    gpu->sim.totalNonbondOutputBuffers  = (gpu->bIncludeGBSA ? 2 * gpu->sim.nonbondOutputBuffers : gpu->sim.nonbondOutputBuffers);
    gpu->sim.outputBuffers              = gpu->sim.totalNonbondOutputBuffers;


    unsigned int outputBuffers = gpu->sim.totalNonbondOutputBuffers;
    for (unsigned int i = 0; i < gpu->sim.paddedNumberOfAtoms; i++)
    {
        if (outputBuffers < gpu->pOutputBufferCounter[i])
        {
            outputBuffers = gpu->pOutputBufferCounter[i];
        }
    }    
    gpu->sim.outputBuffers      = outputBuffers;
    gpu->psForce4               = new CUDAStream<float4>(gpu->sim.paddedNumberOfAtoms, outputBuffers);
    gpu->psBornForce            = new CUDAStream<float>(gpu->sim.paddedNumberOfAtoms, gpu->sim.nonbondOutputBuffers);
    gpu->psBornSum              = new CUDAStream<float>(gpu->sim.paddedNumberOfAtoms, gpu->sim.nonbondOutputBuffers);
    gpu->sim.pForce4            = gpu->psForce4->_pDevStream[0];
    gpu->sim.pForce4a           = gpu->sim.pForce4;
    gpu->sim.pForce4b           = gpu->sim.pForce4 + 1 * gpu->sim.nonbondOutputBuffers * gpu->sim.stride;
    gpu->sim.pBornForce         = gpu->psBornForce->_pDevStream[0];
    gpu->sim.pBornSum           = gpu->psBornSum->_pDevStream[0];

    // Determine local energy paramter offsets for bonded interactions
    gpu->sim.bond_offset        =                                  gpu->psBondParameter->_stride;
    gpu->sim.bond_angle_offset  = gpu->sim.bond_offset           + gpu->psBondAngleParameter->_stride;
    gpu->sim.dihedral_offset    = gpu->sim.bond_angle_offset     + gpu->psDihedralParameter->_stride;
    gpu->sim.rb_dihedral_offset = gpu->sim.dihedral_offset       + gpu->psRbDihedralParameter1->_stride;
    gpu->sim.LJ14_offset        = gpu->sim.rb_dihedral_offset    + gpu->psLJ14Parameter->_stride;
    gpu->sim.localForces_threads_per_block  = (gpu->sim.LJ14_offset / gpu->sim.blocks + 15) & 0xfffffff0;
    if (gpu->sim.localForces_threads_per_block > gpu->sim.max_localForces_threads_per_block)
        gpu->sim.localForces_threads_per_block = gpu->sim.max_localForces_threads_per_block;
    if (gpu->sim.localForces_threads_per_block < 1)
        gpu->sim.localForces_threads_per_block = 1;

    // Flip local force output buffers
    int flip = outputBuffers - 1;
    for (int i = 0; i < (int) gpu->sim.bonds; i++)
    {
        gpu->psBondID->_pSysStream[0][i].z = flip - gpu->psBondID->_pSysStream[0][i].z;
        gpu->psBondID->_pSysStream[0][i].w = flip - gpu->psBondID->_pSysStream[0][i].w;
    }
    for (int i = 0; i < (int) gpu->sim.bond_angles; i++)
    {
        gpu->psBondAngleID1->_pSysStream[0][i].w = flip - gpu->psBondAngleID1->_pSysStream[0][i].w;
        gpu->psBondAngleID2->_pSysStream[0][i].x = flip - gpu->psBondAngleID2->_pSysStream[0][i].x;
        gpu->psBondAngleID2->_pSysStream[0][i].y = flip - gpu->psBondAngleID2->_pSysStream[0][i].y;
    }
    for (int i = 0; i < (int) gpu->sim.dihedrals; i++)
    {
        gpu->psDihedralID2->_pSysStream[0][i].x = flip - gpu->psDihedralID2->_pSysStream[0][i].x;
        gpu->psDihedralID2->_pSysStream[0][i].y = flip - gpu->psDihedralID2->_pSysStream[0][i].y;
        gpu->psDihedralID2->_pSysStream[0][i].z = flip - gpu->psDihedralID2->_pSysStream[0][i].z;
        gpu->psDihedralID2->_pSysStream[0][i].w = flip - gpu->psDihedralID2->_pSysStream[0][i].w;
    }
    for (int i = 0; i < (int) gpu->sim.rb_dihedrals; i++)
    {
        gpu->psRbDihedralID2->_pSysStream[0][i].x = flip - gpu->psRbDihedralID2->_pSysStream[0][i].x;
        gpu->psRbDihedralID2->_pSysStream[0][i].y = flip - gpu->psRbDihedralID2->_pSysStream[0][i].y;
        gpu->psRbDihedralID2->_pSysStream[0][i].z = flip - gpu->psRbDihedralID2->_pSysStream[0][i].z;
        gpu->psRbDihedralID2->_pSysStream[0][i].w = flip - gpu->psRbDihedralID2->_pSysStream[0][i].w;
    }
    for (int i = 0; i < (int) gpu->sim.LJ14s; i++)
    {
        gpu->psLJ14ID->_pSysStream[0][i].z = flip - gpu->psLJ14ID->_pSysStream[0][i].z;
        gpu->psLJ14ID->_pSysStream[0][i].w = flip - gpu->psLJ14ID->_pSysStream[0][i].w;
    }
    gpu->psBondID->Upload();
    gpu->psBondAngleID1->Upload();
    gpu->psBondAngleID2->Upload();
    gpu->psDihedralID2->Upload();
    gpu->psRbDihedralID2->Upload();
    gpu->psLJ14ID->Upload();

    return 1;
}

extern "C"
int gpuBuildThreadBlockWorkList(gpuContext gpu)
{
    const unsigned int atoms = gpu->sim.paddedNumberOfAtoms;
    const unsigned int grid = gpu->grid;
    const unsigned int dim = (atoms + (grid - 1)) / grid;
    const unsigned int cells = dim * (dim + 1) / 2;
    CUDAStream<unsigned int>* psWorkUnit = new CUDAStream<unsigned int>(cells, 1u);
    unsigned int* pWorkList = psWorkUnit->_pSysStream[0];
    gpu->psWorkUnit = psWorkUnit;
    gpu->sim.pWorkUnit = psWorkUnit->_pDevStream[0];
    CUDAStream<unsigned int>* psInteractingWorkUnit = new CUDAStream<unsigned int>(cells, 1u);
    gpu->psInteractingWorkUnit = psInteractingWorkUnit;
    gpu->sim.pInteractingWorkUnit = psInteractingWorkUnit->_pDevStream[0];
    CUDAStream<unsigned int>* psInteractionFlag = new CUDAStream<unsigned int>(cells, 1u);
    gpu->psInteractionFlag = psInteractionFlag;
    gpu->sim.pInteractionFlag = psInteractionFlag->_pDevStream[0];
    CUDAStream<size_t>* psInteractionCount = new CUDAStream<size_t>(1, 1u);
    gpu->psInteractionCount = psInteractionCount;
    gpu->sim.pInteractionCount = psInteractionCount->_pDevStream[0];
    CUDAStream<float4>* psGridBoundingBox = new CUDAStream<float4>(dim, 1u);
    gpu->psGridBoundingBox = psGridBoundingBox;
    gpu->sim.pGridBoundingBox = psGridBoundingBox->_pDevStream[0];
    CUDAStream<float4>* psGridCenter = new CUDAStream<float4>(dim, 1u);
    gpu->psGridCenter = psGridCenter;
    gpu->sim.pGridCenter = psGridCenter->_pDevStream[0];
    gpu->sim.nonbond_workBlock      = gpu->sim.nonbond_threads_per_block / GRID;
    gpu->sim.bornForce2_workBlock   = gpu->sim.bornForce2_threads_per_block / GRID;
    gpu->sim.workUnits = cells;

    // Initialize the CUDPP workspace.
    gpu->cudpp = 0;
    CUDPPConfiguration config;
    config.datatype = CUDPP_UINT;
    config.algorithm = CUDPP_COMPACT;
    config.options = CUDPP_OPTION_FORWARD;
    CUDPPResult result = cudppPlan(&gpu->cudpp, config, cells, 1, 0);
    if (CUDPP_SUCCESS != result)
    {
        printf("Error initializing CUDPP: %d\n", result);
        exit(-1);
    }

    // Increase block count if necessary for extra large molecules that would
    // otherwise overflow the SM workunit buffers
//    int minimumBlocks = (cells + gpu->sim.workUnitsPerSM - 1) / gpu->sim.workUnitsPerSM;
//    if ((int) gpu->sim.nonbond_blocks < minimumBlocks)
//    {
//        gpu->sim.nonbond_blocks = gpu->sim.nonbond_blocks * ((minimumBlocks + gpu->sim.nonbond_blocks - 1) / gpu->sim.nonbond_blocks);
//    }
//    if ((int) gpu->sim.bornForce2_blocks < minimumBlocks)
//    {
//        gpu->sim.bornForce2_blocks = gpu->sim.bornForce2_blocks * ((minimumBlocks + gpu->sim.bornForce2_blocks - 1) / gpu->sim.bornForce2_blocks);
//    }
    gpu->sim.nbWorkUnitsPerBlock            = cells / gpu->sim.nonbond_blocks;
    gpu->sim.nbWorkUnitsPerBlockRemainder   = cells - gpu->sim.nonbond_blocks * gpu->sim.nbWorkUnitsPerBlock;
    gpu->sim.bf2WorkUnitsPerBlock           = cells / gpu->sim.bornForce2_blocks;
    gpu->sim.bf2WorkUnitsPerBlockRemainder  = cells - gpu->sim.bornForce2_blocks * gpu->sim.bf2WorkUnitsPerBlock;
    gpu->sim.interaction_threads_per_block = 64;
    gpu->sim.interaction_blocks = (gpu->sim.workUnits + gpu->sim.interaction_threads_per_block - 1) / gpu->sim.interaction_threads_per_block;

    // Decrease thread count for extra small molecules to spread computation
    // across entire chip
    int activeWorkUnits = gpu->sim.nonbond_blocks * gpu->sim.nonbond_workBlock;
    if (activeWorkUnits > (int) cells)
    {
        int balancedWorkBlock                   = (cells + gpu->sim.nonbond_blocks - 1) / gpu->sim.nonbond_blocks;
        gpu->sim.nonbond_threads_per_block      = balancedWorkBlock * GRID;
        gpu->sim.nonbond_workBlock              = balancedWorkBlock;
    }
    activeWorkUnits = gpu->sim.bornForce2_blocks * gpu->sim.bornForce2_workBlock;
    if (activeWorkUnits > (int) cells)
    {
        int balancedWorkBlock                   = (cells + gpu->sim.bornForce2_blocks - 1) / gpu->sim.bornForce2_blocks;
        gpu->sim.bornForce2_threads_per_block   = balancedWorkBlock * GRID;
        gpu->sim.bornForce2_workBlock           = balancedWorkBlock;
    }

    unsigned int count = 0;
    for (unsigned int y = 0; y < dim; y++)
    {
        for (unsigned int x = y; x < dim; x++)
        {
            pWorkList[count] = (x << 17) | (y << 2);
            count++;
        }
    }

    psWorkUnit->Upload();
    gpuSetConstants(gpu);
    return cells;
}

extern "C"
void gpuBuildExclusionList(gpuContext gpu)
{
    const unsigned int atoms = gpu->sim.paddedNumberOfAtoms;
    const unsigned int grid = gpu->grid;
    const unsigned int dim = (atoms+(grid-1))/grid;
    CUDAStream<unsigned int>* psExclusion = new CUDAStream<unsigned int>((atoms*atoms+grid-1) / grid, 1u);
    gpu->psExclusion = psExclusion;
    gpu->sim.pExclusion = psExclusion->_pDevStream[0];
    unsigned int* pExList = psExclusion->_pSysStream[0];
    unsigned int* pWorkList = gpu->psWorkUnit->_pSysStream[0];
    for (int i = 0; i < psExclusion->_length; ++i)
        pExList[i] = 0xFFFFFFFF;

    // Fill in the exclusions.

    for (int atom1 = 0; atom1 < gpu->exclusions.size(); ++atom1)
    {
        int x = atom1/grid;
        int offset = atom1-x*grid;
        for (int j = 0; j < gpu->exclusions[atom1].size(); ++j)
        {
            int atom2 = gpu->exclusions[atom1][j];
            int y = atom2/grid;
            int index = x*atoms+y*grid+offset;
            pExList[index] &= 0xFFFFFFFF-(1<<(atom2-y*grid));
            int cell = (x > y ? x+y*dim-y*(y+1)/2 : y+x*dim-x*(x+1)/2);
            pWorkList[cell] |= 1;
        }
    }

    // Mark all interactions that involve a padding atom as being excluded.

    for (int atom1 = gpu->natoms; atom1 < atoms; ++atom1)
    {
        int x = atom1/grid;
        int offset1 = atom1-x*grid;
        for (int atom2 = 0; atom2 < atoms; ++atom2)
        {
            int y = atom2/grid;
            int index = x*atoms+y*grid+offset1;
            pExList[index] &= 0xFFFFFFFF-(1<<(atom2-y*grid));
            int offset2 = atom2-y*grid;
            index = y*atoms+x*grid+offset2;
            pExList[index] &= 0xFFFFFFFF-(1<<(atom1-x*grid));
            int cell = (x > y ? x+y*dim-y*(y+1)/2 : y+x*dim-x*(x+1)/2);
            pWorkList[cell] |= 1;
        }
    }
    
    psExclusion->Upload();
    gpu->psWorkUnit->Upload();
    gpuSetConstants(gpu);
}

extern "C"
int gpuSetConstants(gpuContext gpu)
{
    SetCalculateCDLJForcesSim(gpu);
    SetCalculateCDLJObcGbsaForces1Sim(gpu);
    SetCalculateLocalForcesSim(gpu);
    SetCalculateObcGbsaBornSumSim(gpu);
    SetCalculateObcGbsaForces1Sim(gpu);
    SetCalculateObcGbsaForces2Sim(gpu);
    SetCalculateAndersenThermostatSim(gpu);
    SetForcesSim(gpu);
    SetUpdateShakeHSim(gpu);
    SetVerletUpdateSim(gpu);
    SetBrownianUpdateSim(gpu);
    SetSettleSim(gpu);
    SetRandomSim(gpu);

    if (gpu->sm_version >= SM_12)
    {
        SetCalculateObcGbsaForces1_12Sim(gpu);
    }

    return 1;
}

extern "C"
void gpuDumpCoordinates(gpuContext gpu)
{
    gpu->psPosq4->Download();
    gpu->psVelm4->Download();
    (void) printf( "\n\nCoordinates and velocities\n" );
    for (int i = 0; i < gpu->natoms; i++)
    {
        printf("%4d: %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f\n", i, 
            gpu->psPosq4->_pSysStream[0][i].x,
            gpu->psPosq4->_pSysStream[0][i].y,
            gpu->psPosq4->_pSysStream[0][i].z,
            gpu->psPosq4->_pSysStream[0][i].w,
            gpu->psVelm4->_pSysStream[0][i].x,
            gpu->psVelm4->_pSysStream[0][i].y,
            gpu->psVelm4->_pSysStream[0][i].z,
            gpu->psVelm4->_pSysStream[0][i].w
        );
    }
}

bool ISNAN(float f)
{
    return !(f == f);
}

extern "C"
bool gpuCheckData(gpuContext gpu)
{
    gpu->psPosq4->Download();
    gpu->psVelm4->Download();
    gpu->psForce4->Download();
    gpu->psBornForce->Download();
    int violations = 0;
    for (int i = 0; i < gpu->natoms; i++)
    {
        if (ISNAN( gpu->psPosq4->_pSysStream[0][i].x) ||
            ISNAN( gpu->psPosq4->_pSysStream[0][i].y) ||
            ISNAN( gpu->psPosq4->_pSysStream[0][i].z) ||
            ISNAN( gpu->psVelm4->_pSysStream[0][i].x) ||
            ISNAN( gpu->psVelm4->_pSysStream[0][i].y) ||
            ISNAN( gpu->psVelm4->_pSysStream[0][i].z) ||
            ISNAN( gpu->psForce4->_pSysStream[0][i].x) ||
            ISNAN( gpu->psForce4->_pSysStream[0][i].y) ||
            ISNAN( gpu->psForce4->_pSysStream[0][i].z) ||
            ISNAN( gpu->psBornForce->_pSysStream[0][i]))
        {
            printf("%4d: %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f\n", i, 
                gpu->psPosq4->_pSysStream[0][i].x,
                gpu->psPosq4->_pSysStream[0][i].y,
                gpu->psPosq4->_pSysStream[0][i].z,
                gpu->psVelm4->_pSysStream[0][i].x,
                gpu->psVelm4->_pSysStream[0][i].y,
                gpu->psVelm4->_pSysStream[0][i].z,
                gpu->psForce4->_pSysStream[0][i].x,
                gpu->psForce4->_pSysStream[0][i].y,
                gpu->psForce4->_pSysStream[0][i].z,
                gpu->psBornForce->_pSysStream[0][i]
            );
            violations++;
        }
    }


    if (violations > 0)
    {
        printf("%d total violations\n", violations);
        for (int i = 0; i < gpu->natoms; i++)
        {
            float dmin = 99999999.0f;
            int closest = -9999;
            float x = gpu->psPosq4->_pSysStream[0][i].x;
            float y = gpu->psPosq4->_pSysStream[0][i].y;
            float z = gpu->psPosq4->_pSysStream[0][i].z;
            for (int j = 0; j < gpu->natoms; j++)
            {
                if (j != i)
                {
                    float dx = gpu->psPosq4->_pSysStream[0][j].x - x;
                    float dy = gpu->psPosq4->_pSysStream[0][j].y - y;
                    float dz = gpu->psPosq4->_pSysStream[0][j].z - z;
                    float r = sqrt(dx * dx + dy * dy + dz * dz);
                    if (r < dmin)
                    {
                        dmin = r;
                        closest = j;
                    }
                }
            }
            printf("Atom %4d: Closest neighbor is Atom %4d, %11.5e\n", i, closest, dmin);
        }

        gpuDumpAtomData(gpu);

        kClearBornForces(gpu);
        kClearForces(gpu);
        kCPUCalculateLocalForces(gpu);


        // Determine which forces have gone awry
        kClearBornForces(gpu);
        kClearForces(gpu);
        kCalculateCDLJForces(gpu);
        kReduceForces(gpu);
        printf("Nonbond Forces\n");
        gpuDumpForces(gpu);

        kClearBornForces(gpu);
        kClearForces(gpu);
        kCalculateObcGbsaForces1(gpu);
        kReduceObcGbsaBornForces(gpu);
        kCalculateObcGbsaForces2(gpu); 
        kReduceForces(gpu);
        printf("OBC Forces\n");
        gpuDumpForces(gpu);

        kClearBornForces(gpu);
        kClearForces(gpu);
        kCalculateLocalForces(gpu);
        kReduceForces(gpu);
        printf("Local Forces\n");
        gpuDumpForces(gpu);
        kClearBornForces(gpu);
        kClearForces(gpu);
        kReduceForces(gpu);
        printf("Cleared Forces\n");
        gpuDumpForces(gpu);
        
        return false;
    }
    return true;
}

extern "C"
void kCPUCalculate14(gpuContext gpu)
{
    gpu->psPosq4->Download();
    gpu->psForce4->Download();
 //   gpu->psLJ14ID->Download();
 //   gpu->psLJ14Parameter->Download();
    for (int pos = 0; pos < (int) gpu->sim.LJ14s; pos++)
    {
        int4 atom               = gpu->psLJ14ID->_pSysStream[0][pos];
        float4 LJ14             = gpu->psLJ14Parameter->_pSysStream[0][pos];
        float4 a1               = gpu->psPosq4->_pSysStream[0][atom.x];
        float4 a2               = gpu->psPosq4->_pSysStream[0][atom.y];
        float3 d;
        d.x                     = a1.x - a2.x;
        d.y                     = a1.y - a2.y;
        d.z                     = a1.z - a2.z;
        float r2                = d.x * d.x + d.y * d.y + d.z * d.z;
        float inverseR          = 1.0f / sqrt(r2);
        float sig2              = inverseR * LJ14.y;
        sig2                   *= sig2;
        float sig6              = sig2 * sig2 * sig2;
        float dEdR              = LJ14.x * (12.0f * sig6 - 6.0f) * sig6;
        dEdR                   += LJ14.z * inverseR;
        dEdR                   *= inverseR * inverseR;
        unsigned int offsetA    = atom.x + atom.z * gpu->sim.stride;
        unsigned int offsetB    = atom.y + atom.w * gpu->sim.stride;
        float4 forceA           = gpu->psForce4->_pSysStream[0][offsetA];
        float4 forceB           = gpu->psForce4->_pSysStream[0][offsetB];
        d.x                    *= dEdR;
        d.y                    *= dEdR;
        d.z                    *= dEdR;
        forceA.x               += d.x;
        forceA.y               += d.y;
        forceA.z               += d.z;
        forceB.x               -= d.x;
        forceB.y               -= d.y;
        forceB.z               -= d.z;        
        gpu->psForce4->_pSysStream[0][offsetA]   = forceA;
        gpu->psForce4->_pSysStream[0][offsetB]   = forceB;
        printf("%4d: %4d - %4d: %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f\n", pos, atom.x, atom.y, r2, dEdR, sig2, sig6, LJ14.x, LJ14.z); 
    }        
}


extern "C"
void gpuDumpPrimeCoordinates(gpuContext gpu)
{
    gpu->psPosqP4->Download();
    for (int i = 0; i < gpu->natoms; i++)
    {
        printf("%4d: %11.5f %11.5f %11.5f %11.5f\n", i, 
            gpu->psPosqP4->_pSysStream[0][i].x,
            gpu->psPosqP4->_pSysStream[0][i].y,
            gpu->psPosqP4->_pSysStream[0][i].z,
            gpu->psPosqP4->_pSysStream[0][i].w
        );
    }
}

extern "C"
void gpuDumpForces(gpuContext gpu)
{
    gpu->psForce4->Download();
    gpu->psBornForce->Download();
    for (int i = 0; i < gpu->natoms; i++)
    {
        char buff[512];
        sprintf(buff, "%4d: %11.5f %11.5f %11.5f %11.5f\n", i, 
            gpu->psForce4->_pSysStream[0][i].x,
            gpu->psForce4->_pSysStream[0][i].y,
            gpu->psForce4->_pSysStream[0][i].z,
            gpu->psBornForce->_pSysStream[0][i]
        );
//        OutputDebugString(buff);
    }
}

extern "C"
void gpuDumpAtomData(gpuContext gpu)
{
    gpu->psPosq4->Download();
    gpu->psSigEps2->Download();
    gpu->psBornRadii->Download();
    gpu->psObcChain->Download();
    for (int i = 0; i < gpu->natoms; i++)
    {
        char buff[512];
        sprintf(buff, "%4d: %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f\n", i, 
            gpu->psPosq4->_pSysStream[0][i].x,
            gpu->psPosq4->_pSysStream[0][i].y,
            gpu->psPosq4->_pSysStream[0][i].z,
            gpu->psPosq4->_pSysStream[0][i].w,
            gpu->psSigEps2->_pSysStream[0][i].x,
            gpu->psSigEps2->_pSysStream[0][i].y,
            gpu->psBornRadii->_pSysStream[0][i],
            gpu->psObcChain->_pSysStream[0][i]
        );
//        OutputDebugString((LPCWSTR)buff);
    }
}


extern "C"
void gpuSetup(void* pVoid)
{            
    gpuContext gpu = (gpuContext)pVoid;
    // Read parameters
    cout << gpuReadAtomicParameters(gpu, "Data/atomicradii.txt") << " atom types\n";
    cout << gpuReadBondParameters(gpu, "Data/GromacsHarmonicBondParameter.txt") << " bond parameters.\n";
    cout << gpuReadBondAngleParameters(gpu, "Data/GromacsAngleBondParameter.txt") << " bond angle parameters.\n";
    cout << gpuReadDihedralParameters(gpu, "Data/GromacsProperDihedralParameter.txt") << " proper dihedral parameters.\n";
    cout << gpuReadRbDihedralParameters(gpu, "Data/GromacsRbDihedralParameter.txt") << " Ryckaert-Bellemans dihedral parameters.\n";
    cout << gpuReadLJ14Parameters(gpu, "Data/GromacsLJ14Parameter.txt") << " Lennard-Jones 1-4 parameters.\n";
    cout << gpuReadCoulombParameters(gpu, "Data/GromacsLJCoulombParameter.txt") << " Coulomb parameters.\n";
    cout << gpuReadShakeParameters(gpu, "Data/GromacsShakeParameters.txt") << " shake parameters.\n";

    // Build thread block work list
    gpuBuildThreadBlockWorkList(gpu);

    // Build exclusion list
    gpuBuildExclusionList(gpu);
    
    // Create output buffers
    gpuBuildOutputBuffers(gpu);

    // Set constant blocks
    gpuSetConstants(gpu);

    // Initialize randoms
    gpuInitializeRandoms(gpu);

    // Initialize Born Radii;
    kCalculateObcGbsaBornSum(gpu);
    kReduceObcGbsaBornSum(gpu);
    kClearForces(gpu);
    kClearBornForces(gpu);
    return;
}


#define DOT3(v1, v2) (v1.x * v2.x + v1.y * v2.y + v1.z * v2.z)

#define GETNORMEDDOTPRODUCT(v1, v2, dp) \
{ \
    dp          = DOT3(v1, v2); \
    float norm1 = DOT3(v1, v1); \
    float norm2 = DOT3(v2, v2); \
    dp /= sqrt(norm1 * norm2); \
    dp = min(dp, 1.0f); \
    dp = max(dp, -1.0f); \
}

#define CROSS_PRODUCT(v1, v2, c) \
    c.x = v1.y * v2.z - v1.z * v2.y; \
    c.y = v1.z * v2.x - v1.x * v2.z; \
    c.z = v1.x * v2.y - v1.y * v2.x;

#define GETPREFACTORSGIVENANGLECOSINE(cosine, param, dEdR) \
{ \
   float angle          = acos(cosine); \
   float deltaIdeal     = angle - (param.x * (3.14159265f / 180.0f)); \
   dEdR                 = param.y * deltaIdeal; \
}

#define GETANGLEBETWEENTWOVECTORS(v1, v2, angle) \
{ \
    float dp; \
    GETNORMEDDOTPRODUCT(v1, v2, dp); \
    angle = acos(dp); \
}

#define GETANGLECOSINEBETWEENTWOVECTORS(v1, v2, angle, cosine) \
{ \
    GETNORMEDDOTPRODUCT(v1, v2, cosine); \
    angle = acos(cosine); \
}

#define GETDIHEDRALANGLEBETWEENTHREEVECTORS(vector1, vector2, vector3, signVector, cp0, cp1, angle) \
{ \
    CROSS_PRODUCT(vector1, vector2, cp0); \
    CROSS_PRODUCT(vector2, vector3, cp1); \
    GETANGLEBETWEENTWOVECTORS(cp0, cp1, angle); \
    float dp = DOT3(signVector, cp1); \
    angle = (dp >= 0) ? angle : -angle; \
}                                                          

#define GETDIHEDRALANGLECOSINEBETWEENTHREEVECTORS(vector1, vector2, vector3, signVector, cp0, cp1, angle, cosine) \
{ \
    CROSS_PRODUCT(vector1, vector2, cp0); \
    CROSS_PRODUCT(vector2, vector3, cp1); \
    GETANGLECOSINEBETWEENTWOVECTORS(cp0, cp1, angle, cosine); \
    float dp = DOT3(signVector, cp1); \
    angle = (dp >= 0) ? angle : -angle; \
}    

// Calculate Local forces on CPU
extern "C"
void kCPUCalculateLocalForces(gpuContext gpu)
{
    gpu->psPosq4->Download();
    gpu->psForce4->Download();
    gpu->psBondID->Download();
    gpu->psBondParameter->Download();
    gpu->psBondAngleID1->Download();
    gpu->psBondAngleID2->Download();
    gpu->psBondAngleParameter->Download();
    gpu->psDihedralID1->Download();
    gpu->psDihedralID2->Download();
    gpu->psDihedralParameter->Download();
    gpu->psRbDihedralID1->Download();
    gpu->psRbDihedralID2->Download();
    gpu->psRbDihedralParameter1->Download();
    gpu->psRbDihedralParameter2->Download();
    gpu->psLJ14ID->Download();
    gpu->psLJ14Parameter->Download();

    unsigned int pos = 0;
    Vectors V;
    Vectors* A = &V;
    int violations = 0;

    while (pos < gpu->sim.bond_offset)
    {
        if (pos < gpu->sim.bonds)
        {
            int4   atom         = gpu->psBondID->_pSysStream[0][pos];
            float4 atomA        = gpu->psPosq4->_pSysStream[0][atom.x];
            float4 atomB        = gpu->psPosq4->_pSysStream[0][atom.y];
            float2 bond         = gpu->psBondParameter->_pSysStream[0][pos];
            float dx            = atomB.x - atomA.x;
            float dy            = atomB.y - atomA.y;
            float dz            = atomB.z - atomA.z;
            float r2            = dx * dx + dy * dy + dz * dz;
            float r             = sqrt(r2);
            float deltaIdeal    = r - bond.x;
            float dEdR          = bond.y * deltaIdeal;
            dEdR                = (r > 0.0f) ? (dEdR / r) : 0.0f;
            if (fabs(deltaIdeal) > 1.0f)
            {
                printf("Bond %4d: %11.4f %11.4f %11.4f %11.4f %11.4f %11.4f\n", pos, dx, dy, dz, r, deltaIdeal, dEdR);
                violations++;
            }
            dx                 *= dEdR;
            dy                 *= dEdR;
            dz                 *= dEdR;
            unsigned int offsetA                = atom.x + atom.z * gpu->sim.stride;
            unsigned int offsetB                = atom.y + atom.w * gpu->sim.stride;
            float4 forceA                       = gpu->psForce4->_pSysStream[0][offsetA];
            float4 forceB                       = gpu->psForce4->_pSysStream[0][offsetB];
            forceA.x                           += dx;
            forceA.y                           += dy;
            forceA.z                           += dz;
            forceB.x                           -= dx;
            forceB.y                           -= dy;
            forceB.z                           -= dz;
            gpu->psForce4->_pSysStream[0][offsetA]               = forceA;
            gpu->psForce4->_pSysStream[0][offsetB]               = forceB;    
        }
        pos++;
    }
#if 0  
    while (pos < gpu->sim.bond_angle_offset)
    {
        unsigned int pos1   = pos - gpu->sim.bond_offset;
        if (pos1 < gpu->sim.bond_angles)
        {
            int4   atom1        = gpu->psBondAngleID1->_pSysStream[0][pos1];  
            float2 bond_angle   = gpu->psBondAngleParameter->_pSysStream[0][pos1];
            float4 a1           = gpu->psPosq4->_pSysStream[0][atom1.x];
            float4 a2           = gpu->psPosq4->_pSysStream[0][atom1.y];
            float4 a3           = gpu->psPosq4->_pSysStream[0][atom1.z];
            A->v0.x = a2.x - a1.x;
            A->v0.y = a2.y - a1.y;
            A->v0.z = a2.z - a1.z;
            A->v1.x = a2.x - a3.x;
            A->v1.y = a2.y - a3.y;
            A->v1.z = a2.z - a3.z;
            float3 cp;
            CROSS_PRODUCT(A->v0, A->v1, cp);
            float rp = DOT3(cp, cp); //cx * cx + cy * cy + cz * cz;
            rp = max(sqrt(rp), 1.0e-06f);
            float r21       = DOT3(A->v0, A->v0); // dx1 * dx1 + dy1 * dy1 + dz1 * dz1;
            float r23       = DOT3(A->v1, A->v1); // dx2 * dx2 + dy2 * dy2 + dz2 * dz2;
            float dot       = DOT3(A->v0, A->v1); // dx1 * dx2 + dy1 * dy2 + dz1 * dz2;
            float cosine    = dot / sqrt(r21 * r23);
            float dEdR;
            GETPREFACTORSGIVENANGLECOSINE(cosine, bond_angle, dEdR);
            printf("Bond angle %4d %11.4f %11.4f\n", pos1, cosine, dEdR);
            float termA =  dEdR / (r21 * rp);
            float termC = -dEdR / (r23 * rp);
            float3 c21;
            float3 c23;
            CROSS_PRODUCT(A->v0, cp, c21);
            CROSS_PRODUCT(A->v1, cp, c23);
            c21.x *= termA;
            c21.y *= termA;
            c21.z *= termA;
            c23.x *= termC;
            c23.y *= termC;
            c23.z *= termC;
            int2 atom2 = gpu->psBondAngleID2->_pSysStream[0][pos1];
            unsigned int offset = atom1.x + atom1.w * gpu->sim.stride;
            float4 force = gpu->psForce4->_pSysStream[0][offset]; 
            force.x += c21.x;
            force.y += c21.y;
            force.z += c21.z;
            gpu->psForce4->_pSysStream[0][offset] = force;
            offset = atom1.y + atom2.x * gpu->sim.stride;
            force = gpu->psForce4->_pSysStream[0][offset];
            force.x -= (c21.x + c23.x);
            force.y -= (c21.y + c23.y);
            force.z -= (c21.z + c23.z);
            gpu->psForce4->_pSysStream[0][offset] = force;
            offset = atom1.z + atom2.y * gpu->sim.stride;
            force = gpu->psForce4->_pSysStream[0][offset];
            force.x += c23.x;
            force.y += c23.y;
            force.z += c23.z;
            gpu->psForce4->_pSysStream[0][offset] = force;
        }
        pos++;
    }
            
    while (pos < gpu->sim.dihedral_offset)
    {
        unsigned int pos1 = pos - gpu->sim.bond_angle_offset;
        if (pos1 < gpu->sim.dihedrals)
        {
            int4   atom1        = gpu->psDihedralID1->_pSysStream[0][pos1];  
            float4 atomA        = gpu->psPosq4->_pSysStream[0][atom1.x];
            float4 atomB        = gpu->psPosq4->_pSysStream[0][atom1.y];
            float4 atomC        = gpu->psPosq4->_pSysStream[0][atom1.z];
            float4 atomD        = gpu->psPosq4->_pSysStream[0][atom1.w];            
            A->v0.x             = atomA.x - atomB.x;
            A->v0.y             = atomA.y - atomB.y;
            A->v0.z             = atomA.z - atomB.z;
            A->v1.x             = atomC.x - atomB.x;
            A->v1.y             = atomC.y - atomB.y;
            A->v1.z             = atomC.z - atomB.z;
            A->v2.x             = atomC.x - atomD.x;
            A->v2.y             = atomC.y - atomD.y;
            A->v2.z             = atomC.z - atomD.z; 
            float3 cp0, cp1;
            float dihedralAngle;
            GETDIHEDRALANGLEBETWEENTHREEVECTORS(A->v0, A->v1, A->v2, A->v0, cp0, cp1, dihedralAngle);
            float4 dihedral         = gpu->psDihedralParameter->_pSysStream[0][pos1];
            float deltaAngle        = dihedral.z * dihedralAngle - (dihedral.y * 3.14159265f / 180.0f);
            float sinDeltaAngle     = sin(deltaAngle);
            float dEdAngle          = -dihedral.x * dihedral.z * sinDeltaAngle;
            float normCross1        = DOT3(cp0, cp0);
            float normBC            = sqrt(DOT3(A->v1, A->v1));
            float4 ff;
            ff.x                    = (-dEdAngle * normBC) / normCross1;
            float normCross2        = DOT3(cp1, cp1);
            ff.w                    = (dEdAngle * normBC) / normCross2;
            float dp                = 1.0f / DOT3(A->v1, A->v1);
            ff.y                    = DOT3(A->v0, A->v1) * dp;
            ff.z                    = DOT3(A->v2, A->v1) * dp;
            int4  atom2             = gpu->psDihedralID2->_pSysStream[0][pos1];   
            float3 internalF0;
            float3 internalF3;
            float3 s;
            
//            printf("%4d: %9.4f %9.4f %9.4f %9.4f\n", pos1, ff.x, ff.y, ff.z, ff.w);  
            unsigned int offset                 = atom1.x + atom2.x * gpu->sim.stride;
            float4 force                        = gpu->psForce4->_pSysStream[0][offset];
            internalF0.x                        = ff.x * cp0.x; 
            force.x                            += internalF0.x;
            internalF0.y                        = ff.x * cp0.y;
            force.y                            += internalF0.y;
            internalF0.z                        = ff.x * cp0.z;       
            force.z                            += internalF0.z;
            gpu->psForce4->_pSysStream[0][offset]                = force;
            
            printf("Dihedral %4d - 0: %9.4f %9.4f %9.4f\n", pos1, gpu->psForce4->_pSysStream[0][offset], gpu->psForce4->_pSysStream[0][offset + gpu->sim.stride], gpu->psForce4->_pSysStream[0][offset + gpu->sim.stride2]);
            offset                              = atom1.w + atom2.w * gpu->sim.stride;
            force                               = gpu->psForce4->_pSysStream[0][offset];
            internalF3.x                        = ff.w * cp1.x;
            force.x                            += internalF3.x;
            internalF3.y                        = ff.w * cp1.y;
            force.y                            += internalF3.y;
            internalF3.z                        = ff.w * cp1.z;
            force.z                            += internalF3.z;
            gpu->psForce4->_pSysStream[0][offset]                = force;
            
            printf("Dihedral %4d - 3: %9.4f %9.4f %9.4f\n", pos1, gpu->psForce4->_pSysStream[0][offset], gpu->psForce4->_pSysStream[0][offset + gpu->sim.stride], gpu->psForce4->_pSysStream[0][offset + gpu->sim.stride2]);
            s.x                                 = ff.y * internalF0.x - ff.z * internalF3.x;   
            s.y                                 = ff.y * internalF0.y - ff.z * internalF3.y;  
            s.z                                 = ff.y * internalF0.z - ff.z * internalF3.z;        
            offset                              = atom1.y + atom2.y * gpu->sim.stride;
            force                               = gpu->psForce4->_pSysStream[0][offset];
            force.x                            += -internalF0.x + s.x;
            force.y                            += -internalF0.y + s.y;
            force.z                            += -internalF0.z + s.z;
            gpu->psForce4->_pSysStream[0][offset]                = force;
            
            printf("Dihedral %4d - 1: %9.4f %9.4f %9.4f\n", pos1, gpu->psForce4->_pSysStream[0][offset], gpu->psForce4->_pSysStream[0][offset + gpu->sim.stride], gpu->psForce4->_pSysStream[0][offset + gpu->sim.stride2]);
            offset                              = atom1.z + atom2.z * gpu->sim.stride;
            force                               = gpu->psForce4->_pSysStream[0][offset];
            force.x                            += -internalF3.x - s.x;
            force.y                            += -internalF3.y - s.y;
            force.z                            += -internalF3.z - s.z;
            gpu->psForce4->_pSysStream[0][offset]                = force;
            printf("Dihedral %4d - 2: %9.4f %9.4f %9.4f\n", pos1, gpu->psForce4->_pSysStream[0][offset], gpu->psForce4->_pSysStream[0][offset + gpu->sim.stride], gpu->psForce4->_pSysStream[0][offset + gpu->sim.stride2]);
        }        
        pos++;
    }

    while (pos < gpu->sim.rb_dihedral_offset)
    {
        unsigned int pos1 = pos - gpu->sim.dihedral_offset;
        if (pos1 < gpu->sim.rb_dihedrals)
        {
            int4   atom1        = gpu->psRbDihedralID1->_pSysStream[0][pos1];  
            float4 atomA        = gpu->psPosq4->_pSysStream[0][atom1.x];
            float4 atomB        = gpu->psPosq4->_pSysStream[0][atom1.y];
            float4 atomC        = gpu->psPosq4->_pSysStream[0][atom1.z];
            float4 atomD        = gpu->psPosq4->_pSysStream[0][atom1.w];            
            A->v0.x             = atomA.x - atomB.x;
            A->v0.y             = atomA.y - atomB.y;
            A->v0.z             = atomA.z - atomB.z;
            A->v1.x             = atomC.x - atomB.x;
            A->v1.y             = atomC.y - atomB.y;
            A->v1.z             = atomC.z - atomB.z;
            A->v2.x             = atomC.x - atomD.x;
            A->v2.y             = atomC.y - atomD.y;
            A->v2.z             = atomC.z - atomD.z; 
            float3 cp0, cp1;
            float dihedralAngle, cosPhi;
      //      printf("%4d - 0 : %9.4f %9.4f %9.4f\n", pos1, A->v0.x, A->v0.y, A->v0.z); 
      //      printf("%4d - 1 : %9.4f %9.4f %9.4f\n", pos1, A->v1.x, A->v1.y, A->v1.z); 
      //      printf("%4d - 2 : %9.4f %9.4f %9.4f\n", pos1, A->v2.x, A->v2.y, A->v2.z);  
            GETDIHEDRALANGLECOSINEBETWEENTHREEVECTORS(A->v0, A->v1, A->v2, A->v0, cp0, cp1, dihedralAngle, cosPhi);
            if (dihedralAngle < 0.0f )
            {
                dihedralAngle += 3.14159265f;
            } 
            else 
            {
                dihedralAngle -= 3.14159265f;
            }
            cosPhi                  = -cosPhi;
         //   printf("%4d: %9.4f %9.4f\n", pos1, dihedralAngle, cosPhi);
            float4 dihedral1        = gpu->psRbDihedralParameter1->_pSysStream[0][pos1];
            float2 dihedral2        = gpu->psRbDihedralParameter2->_pSysStream[0][pos1];
            float cosFactor         = cosPhi;
            float dEdAngle          = -dihedral1.y;
        //    printf("%4d - 1: %9.4f %9.4f\n", pos1, dEdAngle, 1.0f);
            dEdAngle               -= 2.0f * dihedral1.z * cosFactor;
       //     printf("%4d - 2: %9.4f %9.4f\n", pos1, dEdAngle, cosFactor);
            cosFactor              *= cosPhi;
            dEdAngle               -= 3.0f * dihedral1.w * cosFactor;
     //       printf("%4d - 3: %9.4f %9.4f\n", pos1, dEdAngle, cosFactor);
            cosFactor              *= cosPhi;
            dEdAngle               -= 4.0f * dihedral2.x * cosFactor;
   //         printf("%4d - 4: %9.4f %9.4f\n", pos1, dEdAngle, cosFactor);
            cosFactor              *= cosPhi;
            dEdAngle               -= 5.0f * dihedral2.y * cosFactor;
 //           printf("%4d - 5: %9.4f %9.4f\n", pos1, dEdAngle, cosFactor);
            dEdAngle               *= sin(dihedralAngle);  
//            printf("%4d - f: %9.4f\n", pos1, dEdAngle);
            
            float normCross1        = DOT3(cp0, cp0);
            float normBC            = sqrt(DOT3(A->v1, A->v1));
            float4 ff;
            ff.x                    = (-dEdAngle * normBC) / normCross1;
            float normCross2        = DOT3(cp1, cp1);
            ff.w                    = (dEdAngle * normBC) / normCross2;
            float dp                = 1.0f / DOT3(A->v1, A->v1);
            ff.y                    = DOT3(A->v0, A->v1) * dp;
            ff.z                    = DOT3(A->v2, A->v1) * dp;
            int4  atom2             = gpu->psRbDihedralID2->_pSysStream[0][pos1];   
            float3 internalF0;
            float3 internalF3;
            float3 s;
            
            printf("RB Dihedral %4d: %9.4f %9.4f %9.4f %9.4f\n", pos1, ff.x, ff.y, ff.z, ff.w);  
            unsigned int offset                 = atom1.x + atom2.x * gpu->sim.stride;
            float4 force                        = gpu->psForce4->_pSysStream[0][offset];
            internalF0.x                        = ff.x * cp0.x; 
            force.x                            += internalF0.x;
            internalF0.y                        = ff.x * cp0.y;
            force.y                            += internalF0.y;
            internalF0.z                        = ff.x * cp0.z;       
            force.z                            += internalF0.z;
            gpu->psForce4->_pSysStream[0][offset]                = force;
            
            printf("RB Dihedral %4d - 0: %9.4f %9.4f %9.4f\n", pos1, gpu->psForce4->_pSysStream[0][offset], gpu->psForce4->_pSysStream[0][offset + gpu->sim.stride], gpu->psForce4->_pSysStream[0][offset + gpu->sim.stride2]);
            offset                              = atom1.w + atom2.w * gpu->sim.stride;
            force                               = gpu->psForce4->_pSysStream[0][offset];
            internalF3.x                        = ff.w * cp1.x;
            force.x                            += internalF3.x;
            internalF3.y                        = ff.w * cp1.y;
            force.y                            += internalF3.y;
            internalF3.z                        = ff.w * cp1.z;
            force.z                            += internalF3.z;
            gpu->psForce4->_pSysStream[0][offset]                = force;
            
            printf("RB Dihedral %4d - 3: %9.4f %9.4f %9.4f\n", pos1, gpu->psForce4->_pSysStream[0][offset], gpu->psForce4->_pSysStream[0][offset + gpu->sim.stride], gpu->psForce4->_pSysStream[0][offset + gpu->sim.stride2]);
            s.x                                 = ff.y * internalF0.x - ff.z * internalF3.x;   
            s.y                                 = ff.y * internalF0.y - ff.z * internalF3.y;  
            s.z                                 = ff.y * internalF0.z - ff.z * internalF3.z;        
            offset                              = atom1.y + atom2.y * gpu->sim.stride;
            force                               = gpu->psForce4->_pSysStream[0][offset];
            force.x                            += -internalF0.x + s.x;
            force.y                            += -internalF0.y + s.y;
            force.z                            += -internalF0.z + s.z;
            gpu->psForce4->_pSysStream[0][offset]                = force;
            printf("RB Dihedral %4d - 1: %9.4f %9.4f %9.4f\n", pos1, gpu->psForce4->_pSysStream[0][offset], gpu->psForce4->_pSysStream[0][offset + gpu->sim.stride], gpu->psForce4->_pSysStream[0][offset + gpu->sim.stride2]);
            offset                              = atom1.z + atom2.z * gpu->sim.stride;
            force                               = gpu->psForce4->_pSysStream[0][offset];
            force.x                            += -internalF3.x - s.x;
            force.y                            += -internalF3.y - s.y;
            force.z                            += -internalF3.z - s.z;
            gpu->psForce4->_pSysStream[0][offset]                = force;
     //       printf("%4d - 2: %9.4f %9.4f %9.4f\n", pos1, gpu->psForce4->_pSysStream[0][offset], gpu->psForce4->_pSysStream[0][offset + gpu->sim.stride], gpu->psForce4->_pSysStream[0][offset + gpu->sim.stride2]);
        }            
        pos++;
    }   

    while (pos < gpu->sim.LJ14_offset)
    {  
        unsigned int pos1       = pos - gpu->sim.rb_dihedral_offset;
        if (pos1 < gpu->sim.LJ14s)
        {
            int4 atom               = gpu->psLJ14ID->_pSysStream[0][pos1];
            float4 LJ14             = gpu->psLJ14Parameter->_pSysStream[0][pos1];
            float4 a1               = gpu->psPosq4->_pSysStream[0][atom.x];
            float4 a2               = gpu->psPosq4->_pSysStream[0][atom.y];
            float3 d;
            d.x                     = a1.x - a2.x;
            d.y                     = a1.y - a2.y;
            d.z                     = a1.z - a2.z;
            float r2                = DOT3(d, d);
            float inverseR          = 1.0f / sqrt(r2);
            float sig2              = inverseR * LJ14.y;
            sig2                   *= sig2;
            float sig6              = sig2 * sig2 * sig2;
            float dEdR              = LJ14.x * (12.0f * sig6 - 6.0f) * sig6;
            dEdR                   += LJ14.z * inverseR;
            dEdR                   *= inverseR * inverseR;
            unsigned int offsetA    = atom.x + atom.z * gpu->sim.stride;
            unsigned int offsetB    = atom.y + atom.w * gpu->sim.stride;
            float4 forceA           = gpu->psForce4->_pSysStream[0][offsetA];
            float4 forceB           = gpu->psForce4->_pSysStream[0][offsetB];
            d.x                    *= dEdR;
            d.y                    *= dEdR;
            d.z                    *= dEdR;
            forceA.x               += d.x;
            forceA.y               += d.y;
            forceA.z               += d.z;
            forceB.x               -= d.x;
            forceB.y               -= d.y;
            forceB.z               -= d.z;        
            printf("LJ14 %d: %11.4f %11.4f %11.4f\n", pos1, d.x, d.y, d.z);
            gpu->psForce4->_pSysStream[0][offsetA]   = forceA;
            gpu->psForce4->_pSysStream[0][offsetB]   = forceB;
        }        
        pos++;
    }
#endif

    if (violations > 0)
    {
        gpuDumpCoordinates(gpu);
        gpuDumpForces(gpu);
    }
}

static FILE* getWriteToFilePtr( char* fname, int step )
{
   std::stringstream fileName;
   fileName << fname << "_";
   fileName << step;
   fileName << ".txt";
   FILE* filePtr = fopen( fileName.str().c_str(), "w" );
   if( filePtr == NULL ){
      (void) fprintf( stderr, "Could not open file=<%s> for writitng.", fileName.str().c_str() );
      exit(-1);
   }
   return filePtr;
}

extern "C" {
static void printValues( FILE* filePtr, int index, int numberOfValues, float* values )
{
   int i;
   (void) fprintf( filePtr, "%5d ", index );
   for ( i = 0; i < numberOfValues; i++ ) { 
      (void) fprintf( filePtr, " %18.10e", values[i] );
   }
   (void) fprintf( filePtr, "\n" );
   (void) fflush( filePtr );
} 
}

extern "C"
void WriteArrayToFile1( gpuContext gpu, char* fname, int step, CUDAStream<float>* psPos, int numPrint )
{
   int i;
   static const int numberOfValues = 1;
   FILE* filePtr = getWriteToFilePtr( fname, step );
   float values[numberOfValues];
   psPos->Download();

   numPrint = (numPrint > 0 && (numPrint < gpu->natoms)) ? numPrint : gpu->natoms;
   for ( i = 0; i < numPrint; i++ ) { 
      values[0] = psPos->_pSysStream[0][i];
      printValues( filePtr, i, numberOfValues, values ); 
   }
   for ( i = gpu->natoms - numPrint; i < gpu->natoms; i++ ) { 
      values[0] = psPos->_pSysStream[0][i];
      printValues( filePtr, i, numberOfValues, values ); 
   }
   (void) fclose( filePtr );
}

extern "C"
void WriteArrayToFile2( gpuContext gpu, char* fname, int step, CUDAStream<float2>* psPos, int numPrint )
{
   int i;
   static const int numberOfValues = 2;
   FILE* filePtr = getWriteToFilePtr( fname, step );
   float values[numberOfValues];
   psPos->Download();

   numPrint = (numPrint > 0 && (numPrint < gpu->natoms)) ? numPrint : gpu->natoms;
   for ( i = 0; i < numPrint; i++ ) { 
      values[0] = psPos->_pSysStream[0][i].x;
      values[1] = psPos->_pSysStream[0][i].y;
      printValues( filePtr, i, numberOfValues, values ); 
   }
   for ( i = gpu->natoms - numPrint; i < gpu->natoms; i++ ) { 
      values[0] = psPos->_pSysStream[0][i].x;
      values[1] = psPos->_pSysStream[0][i].y;
      printValues( filePtr, i, numberOfValues, values ); 
   }
   (void) fclose( filePtr );
}

extern "C"
void WriteArrayToFile4( gpuContext gpu, char* fname, int step, CUDAStream<float4>* psPos, int numPrint )
{
   int i;
   static const int numberOfValues = 4;
   FILE* filePtr = getWriteToFilePtr( fname, step );
   float values[numberOfValues];
   psPos->Download();

   numPrint = (numPrint > 0 && (numPrint < gpu->natoms)) ? numPrint : gpu->natoms;
   for ( i = 0; i < numPrint; i++ ) { 
      values[0] = psPos->_pSysStream[0][i].x;
      values[1] = psPos->_pSysStream[0][i].y;
      values[2] = psPos->_pSysStream[0][i].z;
      values[3] = psPos->_pSysStream[0][i].w;
      printValues( filePtr, i, numberOfValues, values ); 
   }
   for ( i = gpu->natoms - numPrint; i < gpu->natoms; i++ ) { 
      values[0] = psPos->_pSysStream[0][i].x;
      values[1] = psPos->_pSysStream[0][i].y;
      values[2] = psPos->_pSysStream[0][i].z;
      values[3] = psPos->_pSysStream[0][i].w;
      printValues( filePtr, i, numberOfValues, values ); 
   }
   (void) fclose( filePtr );
}

extern "C"
void gpuDumpObcInfo(gpuContext gpu)
{
    gpu->psPosq4->Download();
    gpu->psBornRadii->Download();
    gpu->psObcData->Download();
    gpu->psBornSum->Download();
    printf( "\n\nObc Info xyzw Brad atomR scaledAtomR\n" );
    for (int i = 0; i < gpu->natoms; i++)
    {
        printf("%4d: %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f\n", i, 
            gpu->psPosq4->_pSysStream[0][i].x,
            gpu->psPosq4->_pSysStream[0][i].y,
            gpu->psPosq4->_pSysStream[0][i].z,
            gpu->psPosq4->_pSysStream[0][i].w,
            gpu->psBornRadii->_pSysStream[0][i],
            gpu->psBornSum->_pSysStream[0][i],
            gpu->psObcData->_pSysStream[0][i].x,
            gpu->psObcData->_pSysStream[0][i].y
        );
    }
}

extern "C"
void gpuDumpObcLoop1(gpuContext gpu)
{
    float compF;
    gpu->psForce4->Download();
    gpu->psBornRadii->Download();
    gpu->psBornForce->Download();
    gpu->psObcChain->Download();
    gpu->psBornSum->Download();
    printf( "\n\nObc F3 BrnR BrnF Chn\n" );
    for (int i = 0; i < gpu->natoms; i++)
    {
	     compF = gpu->psBornForce->_pSysStream[0][i]/(gpu->psBornRadii->_pSysStream[0][i]*gpu->psBornRadii->_pSysStream[0][i]*gpu->psObcChain->_pSysStream[0][i]);
        printf("%4d: %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f\n", i, 
            gpu->psForce4->_pSysStream[0][i].x,
            gpu->psForce4->_pSysStream[0][i].y,
            gpu->psForce4->_pSysStream[0][i].z,
//            gpu->psForce4->_pSysStream[0][i].w,
            gpu->psBornRadii->_pSysStream[0][i],
				compF,
            gpu->psBornForce->_pSysStream[0][i],
//            gpu->psBornSum->_pSysStream[0][i],
            gpu->psObcChain->_pSysStream[0][i]
        );
    }
}

static void tagAtomsInMolecule(int atom, int molecule, vector<int>& atomMolecule, vector<vector<int> >& atomBonds)
{
    // Recursively tag atoms as belonging to a particular molecule.

    atomMolecule[atom] = molecule;
    for (int i = 0; i < atomBonds[atom].size(); i++)
        if (atomMolecule[atomBonds[atom][i]] == -1)
            tagAtomsInMolecule(atomBonds[atom][i], molecule, atomMolecule, atomBonds);
}

static void findMoleculeGroups(gpuContext gpu)
{
    // First make a list of constraints for future use.

    vector<Constraint> constraints;
    for (int i = 0; i < gpu->sim.ShakeConstraints; i++)
    {
        int atom1 = gpu->psShakeID->_pSysData[i].x;
        int atom2 = gpu->psShakeID->_pSysData[i].y;
        int atom3 = gpu->psShakeID->_pSysData[i].z;
        int atom4 = gpu->psShakeID->_pSysData[i].w;
        float distance2 = gpu->psShakeParameter->_pSysData[i].z;
        constraints.push_back(Constraint(atom1, atom2, distance2));
        if (atom3 != -1)
            constraints.push_back(Constraint(atom1, atom3, distance2));
        if (atom4 != -1)
            constraints.push_back(Constraint(atom1, atom3, distance2));
    }
    for (int i = 0; i < gpu->sim.settleConstraints; i++)
    {
        int atom1 = gpu->psSettleID->_pSysData[i].x;
        int atom2 = gpu->psSettleID->_pSysData[i].y;
        int atom3 = gpu->psSettleID->_pSysData[i].z;
        float distance12 = gpu->psSettleParameter->_pSysData[i].x;
        float distance23 = gpu->psSettleParameter->_pSysData[i].y;
        constraints.push_back(Constraint(atom1, atom2, distance12*distance12));
        constraints.push_back(Constraint(atom1, atom3, distance12*distance12));
        constraints.push_back(Constraint(atom2, atom3, distance23*distance23));
    }

    // First make a list of every other atom to which each atom is connect by a bond or constraint.

    int numAtoms = gpu->natoms;
    vector<vector<int> > atomBonds(numAtoms);
    for (int i = 0; i < gpu->sim.bonds; i++)
    {
        int atom1 = gpu->psBondID->_pSysData[i].x;
        int atom2 = gpu->psBondID->_pSysData[i].y;
        atomBonds[atom1].push_back(atom2);
        atomBonds[atom2].push_back(atom1);
    }
    for (int i = 0; i < constraints.size(); i++)
    {
        int atom1 = constraints[i].atom1;
        int atom2 = constraints[i].atom2;
        atomBonds[atom1].push_back(atom2);
        atomBonds[atom2].push_back(atom1);
    }

    // Now tag atoms by which molecule they belong to.

    vector<int> atomMolecule(numAtoms, -1);
    int numMolecules = 0;
    for (int i = 0; i < numAtoms; i++)
        if (atomMolecule[i] == -1)
            tagAtomsInMolecule(i, numMolecules++, atomMolecule, atomBonds);
    vector<vector<int> > atomIndices(numMolecules);
    for (int i = 0; i < numAtoms; i++)
        atomIndices[atomMolecule[i]].push_back(i);

    // Construct a description of each molecule.

    vector<Molecule> molecules(numMolecules);
    for (int i = 0; i < numMolecules; i++)
        molecules[i].atoms = atomIndices[i];
    for (int i = 0; i < gpu->sim.bonds; i++)
    {
        int atom1 = gpu->psBondID->_pSysData[i].x;
        molecules[atomMolecule[atom1]].bonds.push_back(i);
    }
    for (int i = 0; i < gpu->sim.bond_angles; i++)
    {
        int atom1 = gpu->psBondAngleID1->_pSysData[i].x;
        molecules[atomMolecule[atom1]].angles.push_back(i);
    }
    for (int i = 0; i < gpu->sim.dihedrals; i++)
    {
        int atom1 = gpu->psDihedralID1->_pSysData[i].x;
        molecules[atomMolecule[atom1]].periodicTorsions.push_back(i);
    }
    for (int i = 0; i < gpu->sim.rb_dihedrals; i++)
    {
        int atom1 = gpu->psRbDihedralID1->_pSysData[i].x;
        molecules[atomMolecule[atom1]].rbTorsions.push_back(i);
    }
    for (int i = 0; i < constraints.size(); i++)
    {
        molecules[atomMolecule[constraints[i].atom1]].constraints.push_back(i);
    }

    // Sort them into groups of identical molecules.

    vector<Molecule> uniqueMolecules;
    vector<vector<int> > moleculeInstances;
    for (int molIndex = 0; molIndex < molecules.size(); molIndex++)
    {
        Molecule& mol = molecules[molIndex];

        // See if it is identical to another molecule.

        bool isNew = true;
        for (int j = 0; j < uniqueMolecules.size() && isNew; j++)
        {
            Molecule& mol2 = uniqueMolecules[j];
            bool identical = true;
            if (mol.atoms.size() != mol2.atoms.size() || mol.bonds.size() != mol2.bonds.size()
                    || mol.angles.size() != mol2.angles.size() || mol.periodicTorsions.size() != mol2.periodicTorsions.size()
                    || mol.rbTorsions.size() != mol2.rbTorsions.size() || mol.constraints.size() != mol2.constraints.size())
                identical = false;
            int atomOffset = mol2.atoms[0]-mol.atoms[0];
            float4* posq = gpu->psPosq4->_pSysData;
            float4* velm = gpu->psVelm4->_pSysData;
            float2* sigeps = gpu->psSigEps2->_pSysData;
            for (int i = 0; i < mol.atoms.size() && identical; i++)
                if (mol.atoms[i] != mol2.atoms[i]-atomOffset || posq[mol.atoms[i]].w != posq[mol2.atoms[i]].w ||
                        velm[mol.atoms[i]].w != velm[mol2.atoms[i]].w || sigeps[mol.atoms[i]].x != sigeps[mol2.atoms[i]].x ||
                        sigeps[mol.atoms[i]].y != sigeps[mol2.atoms[i]].y)
                    identical = false;
            int4* bondID = gpu->psBondID->_pSysData;
            float2* bondParam = gpu->psBondParameter->_pSysData;
            for (int i = 0; i < mol.bonds.size() && identical; i++)
                if (bondID[mol.bonds[i]].x != bondID[mol2.bonds[i]].x-atomOffset || bondID[mol.bonds[i]].y != bondID[mol2.bonds[i]].y-atomOffset ||
                        bondParam[mol.bonds[i]].x != bondParam[mol2.bonds[i]].x || bondParam[mol.bonds[i]].y != bondParam[mol2.bonds[i]].y)
                    identical = false;
            int4* angleID = gpu->psBondAngleID1->_pSysData;
            float2* angleParam = gpu->psBondAngleParameter->_pSysData;
            for (int i = 0; i < mol.angles.size() && identical; i++)
                if (angleID[mol.angles[i]].x != angleID[mol2.angles[i]].x-atomOffset ||
                        angleID[mol.angles[i]].y != angleID[mol2.angles[i]].y-atomOffset ||
                        angleID[mol.angles[i]].z != angleID[mol2.angles[i]].z-atomOffset ||
                        angleParam[mol.angles[i]].x != angleParam[mol2.angles[i]].x ||
                        angleParam[mol.angles[i]].y != angleParam[mol2.angles[i]].y)
                    identical = false;
            int4* periodicID = gpu->psDihedralID1->_pSysData;
            float4* periodicParam = gpu->psDihedralParameter->_pSysData;
            for (int i = 0; i < mol.periodicTorsions.size() && identical; i++)
                if (periodicID[mol.periodicTorsions[i]].x != periodicID[mol2.periodicTorsions[i]].x-atomOffset ||
                        periodicID[mol.periodicTorsions[i]].y != periodicID[mol2.periodicTorsions[i]].y-atomOffset ||
                        periodicID[mol.periodicTorsions[i]].z != periodicID[mol2.periodicTorsions[i]].z-atomOffset ||
                        periodicID[mol.periodicTorsions[i]].w != periodicID[mol2.periodicTorsions[i]].w-atomOffset ||
                        periodicParam[mol.periodicTorsions[i]].x != periodicParam[mol2.periodicTorsions[i]].x ||
                        periodicParam[mol.periodicTorsions[i]].y != periodicParam[mol2.periodicTorsions[i]].y ||
                        periodicParam[mol.periodicTorsions[i]].z != periodicParam[mol2.periodicTorsions[i]].z)
                    identical = false;
            int4* rbID = gpu->psRbDihedralID1->_pSysData;
            float4* rbParam1 = gpu->psRbDihedralParameter1->_pSysData;
            float2* rbParam2 = gpu->psRbDihedralParameter2->_pSysData;
            for (int i = 0; i < mol.rbTorsions.size() && identical; i++)
                if (rbID[mol.rbTorsions[i]].x != rbID[mol2.rbTorsions[i]].x-atomOffset ||
                        rbID[mol.rbTorsions[i]].y != rbID[mol2.rbTorsions[i]].y-atomOffset ||
                        rbID[mol.rbTorsions[i]].z != rbID[mol2.rbTorsions[i]].z-atomOffset ||
                        rbID[mol.rbTorsions[i]].w != rbID[mol2.rbTorsions[i]].w-atomOffset ||
                        rbParam1[mol.rbTorsions[i]].x != rbParam1[mol2.rbTorsions[i]].x ||
                        rbParam1[mol.rbTorsions[i]].y != rbParam1[mol2.rbTorsions[i]].y ||
                        rbParam1[mol.rbTorsions[i]].z != rbParam1[mol2.rbTorsions[i]].z ||
                        rbParam1[mol.rbTorsions[i]].w != rbParam1[mol2.rbTorsions[i]].w ||
                        rbParam2[mol.rbTorsions[i]].x != rbParam2[mol2.rbTorsions[i]].x ||
                        rbParam2[mol.rbTorsions[i]].y != rbParam2[mol2.rbTorsions[i]].y)
                    identical = false;
            for (int i = 0; i < mol.constraints.size() && identical; i++)
                if (constraints[mol.constraints[i]].atom1 != constraints[mol2.constraints[i]].atom1-atomOffset ||
                        constraints[mol.constraints[i]].atom2 != constraints[mol2.constraints[i]].atom2-atomOffset ||
                        constraints[mol.constraints[i]].distance2 != constraints[mol2.constraints[i]].distance2)
                    identical = false;
            if (identical)
            {
                moleculeInstances[j].push_back(mol.atoms[0]);
                isNew = false;
            }
        }
        if (isNew)
        {
            uniqueMolecules.push_back(mol);
            moleculeInstances.push_back(vector<int>());
            moleculeInstances[moleculeInstances.size()-1].push_back(mol.atoms[0]);
        }
    }
    gpu->moleculeGroups.resize(moleculeInstances.size());
    for (int i = 0; i < moleculeInstances.size(); i++)
    {
        gpu->moleculeGroups[i].instances = moleculeInstances[i];
        vector<int>& atoms = uniqueMolecules[i].atoms;
        gpu->moleculeGroups[i].atoms.resize(atoms.size());
        for (int j = 0; j < atoms.size(); j++)
            gpu->moleculeGroups[i].atoms[j] = atoms[j]-atoms[0];
    }
}

extern "C"
void gpuReorderAtoms(gpuContext gpu)
{
    if (gpu->natoms == 0 || gpu->sim.nonbondedCutoffSqr == 0.0)
        return;
    if (gpu->moleculeGroups.size() == 0)
        findMoleculeGroups(gpu);

    // Find the range of positions and the number of bins along each axis.

    int numAtoms = gpu->natoms;
    gpu->psPosq4->Download();
    gpu->psVelm4->Download();
    float4* posq = gpu->psPosq4->_pSysData;
    float4* velm = gpu->psVelm4->_pSysData;
    float minx = posq[0].x, maxx = posq[0].x;
    float miny = posq[0].y, maxy = posq[0].y;
    float minz = posq[0].z, maxz = posq[0].z;
    if (gpu->sim.nonbondedMethod == PERIODIC)
    {
        minx = miny = minz = 0.0;
        maxx = gpu->sim.periodicBoxSizeX;
        maxy = gpu->sim.periodicBoxSizeY;
        maxz = gpu->sim.periodicBoxSizeZ;
    }
    else
    {
        for (int i = 1; i < numAtoms; i++)
        {
            minx = min(minx, posq[i].x);
            maxx = max(maxx, posq[i].x);
            miny = min(miny, posq[i].y);
            maxy = max(maxy, posq[i].y);
            minz = min(minz, posq[i].z);
            maxz = max(maxz, posq[i].z);
        }
    }
    float binWidth = 0.2*sqrt(gpu->sim.nonbondedCutoffSqr);
    int xbins = 1 + (int) ((maxx-minx)/binWidth);
    int ybins = 1 + (int) ((maxy-miny)/binWidth);
    int zbins = 1 + (int) ((maxz-minz)/binWidth);

    // Loop over each group of identical molecules and reorder them.

    vector<int> originalIndex(numAtoms);
    vector<float4> newPosq(numAtoms);
    vector<float4> newVelm(numAtoms);
    for (int group = 0; group < gpu->moleculeGroups.size(); group++)
    {
        // Find the center of each molecule.

        gpuMoleculeGroup& mol = gpu->moleculeGroups[group];
        int numMolecules = mol.instances.size();
        vector<int>& atoms = mol.atoms;
        vector<float3> molPos(numMolecules);
        for (int i = 0; i < numMolecules; i++)
        {
            molPos[i].x = 0.0f;
            molPos[i].y = 0.0f;
            molPos[i].z = 0.0f;
            for (int j = 0; j < atoms.size(); j++)
            {
                int atom = atoms[j]+mol.instances[i];
                molPos[i].x += posq[atom].x;
                molPos[i].y += posq[atom].y;
                molPos[i].z += posq[atom].z;
            }
            molPos[i].x /= atoms.size();
            molPos[i].y /= atoms.size();
            molPos[i].z /= atoms.size();
        }
        if (gpu->sim.nonbondedMethod == PERIODIC)
        {
            // Move each molecule position into the same box.

            for (int i = 0; i < numMolecules; i++)
            {
                molPos[i].x -= floor(molPos[i].x/gpu->sim.periodicBoxSizeX)*gpu->sim.periodicBoxSizeX;
                molPos[i].y -= floor(molPos[i].y/gpu->sim.periodicBoxSizeY)*gpu->sim.periodicBoxSizeY;
                molPos[i].z -= floor(molPos[i].z/gpu->sim.periodicBoxSizeZ)*gpu->sim.periodicBoxSizeZ;
            }
        }

        // Select a bin for each molecule, then sort them by bin.

        vector<pair<int, int> > molBins(numMolecules);
        for (int i = 0; i < numMolecules; i++)
        {
            int x = (int) ((molPos[i].x-minx)/binWidth);
            int y = (int) ((molPos[i].y-miny)/binWidth);
            int z = (int) ((molPos[i].z-minz)/binWidth);
            int yodd = y&1;
            int zodd = z&1;
            int bin = z*xbins*ybins;
            bin += (zodd ? ybins-y : y)*xbins;
            bin += (yodd ? xbins-x : x);
            molBins[i] = pair<int, int>(bin, i);
        }
        sort(molBins.begin(), molBins.end());

        // Reorder the atoms.

        for (int i = 0; i < numMolecules; i++)
        {
            for (int j = 0; j < atoms.size(); j++)
            {
                int oldIndex = mol.instances[molBins[i].second]+atoms[j];
                int newIndex = mol.instances[i]+atoms[j];
                originalIndex[newIndex] = gpu->psAtomIndex->_pSysStream[0][oldIndex];
                newPosq[newIndex] = posq[oldIndex];
                newVelm[newIndex] = velm[oldIndex];
            }
        }
    }

    // Update the streams.

    for (int i = 0; i < numAtoms; i++)
        posq[i] = newPosq[i];
    gpu->psPosq4->Upload();
    for (int i = 0; i < numAtoms; i++)
        velm[i] = newVelm[i];
    gpu->psVelm4->Upload();
    for (int i = 0; i < numAtoms; i++)
        gpu->psAtomIndex->_pSysData[i] = originalIndex[i];
    gpu->psAtomIndex->Upload();
}
