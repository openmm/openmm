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
#include "hilbert.h"
#include "OpenMMException.h"

using OpenMM::OpenMMException;

struct ShakeCluster {
    int centralID;
    int peripheralID[3];
    int size;
    bool valid;
    float distance;
    float centralInvMass, peripheralInvMass;
    ShakeCluster() : valid(true) {
    }
    ShakeCluster(int centralID, float invMass) : centralID(centralID), centralInvMass(invMass), size(0), valid(true) {
    }
    void addAtom(int id, float dist, float invMass) {
        if (size == 3 || (size > 0 && dist != distance) || (size > 0 && invMass != peripheralInvMass))
            valid = false;
        else {
            peripheralID[size++] = id;
            distance = dist;
            peripheralInvMass = invMass;
        }
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
void gpuSetBondParameters(gpuContext gpu, const vector<int>& atom1, const vector<int>& atom2, const vector<float>& length, const vector<float>& k)
{
    int bonds = atom1.size();
    gpu->sim.bonds                              = bonds;
    CUDAStream<int4>* psBondID                  = new CUDAStream<int4>(bonds, 1, "BondID");
    gpu->psBondID                               = psBondID;
    gpu->sim.pBondID                            = psBondID->_pDevStream[0];
    CUDAStream<float2>* psBondParameter         = new CUDAStream<float2>(bonds, 1, "BondParameter");
    gpu->psBondParameter                        = psBondParameter;
    gpu->sim.pBondParameter                     = psBondParameter->_pDevStream[0];
    for (int i = 0; i < bonds; i++)
    {
        (*psBondID)[i].x = atom1[i];
        (*psBondID)[i].y = atom2[i];
        (*psBondParameter)[i].x = length[i];
        (*psBondParameter)[i].y = k[i];
        psBondID->_pSysData[i].z = gpu->pOutputBufferCounter[psBondID->_pSysData[i].x]++;
        psBondID->_pSysData[i].w = gpu->pOutputBufferCounter[psBondID->_pSysData[i].y]++;
#if (DUMP_PARAMETERS == 1)                
        cout << 
            i << " " << 
            (*psBondID)[i].x << " " <<
            (*psBondID)[i].y << " " <<
            (*psBondID)[i].z << " " <<
            (*psBondID)[i].w << " " <<
            (*psBondParameter)[i].x << " " <<
            (*psBondParameter)[i].y <<
            endl;
#endif
    }
    psBondID->Upload();
    psBondParameter->Upload();
}

extern "C"
void gpuSetBondAngleParameters(gpuContext gpu, const vector<int>& atom1, const vector<int>& atom2, const vector<int>& atom3,
        const vector<float>& angle, const vector<float>& k)
{
    int bond_angles = atom1.size();
    gpu->sim.bond_angles                        = bond_angles;
    CUDAStream<int4>* psBondAngleID1            = new CUDAStream<int4>(bond_angles, 1, "BondAngleID1");
    gpu->psBondAngleID1                         = psBondAngleID1;
    gpu->sim.pBondAngleID1                      = psBondAngleID1->_pDevStream[0];
    CUDAStream<int2>* psBondAngleID2            = new CUDAStream<int2>(bond_angles, 1, "BondAngleID2");
    gpu->psBondAngleID2                         = psBondAngleID2;
    gpu->sim.pBondAngleID2                      = psBondAngleID2->_pDevStream[0];
    CUDAStream<float2>* psBondAngleParameter    = new CUDAStream<float2>(bond_angles, 1, "BondAngleParameter");
    gpu->psBondAngleParameter                   = psBondAngleParameter;
    gpu->sim.pBondAngleParameter                = psBondAngleParameter->_pDevStream[0];        

    for (int i = 0; i < bond_angles; i++)
    {
        (*psBondAngleID1)[i].x = atom1[i];
        (*psBondAngleID1)[i].y = atom2[i];
        (*psBondAngleID1)[i].z = atom3[i];
        (*psBondAngleParameter)[i].x = angle[i];
        (*psBondAngleParameter)[i].y = k[i];
        psBondAngleID1->_pSysData[i].w = gpu->pOutputBufferCounter[psBondAngleID1->_pSysData[i].x]++;
        psBondAngleID2->_pSysData[i].x = gpu->pOutputBufferCounter[psBondAngleID1->_pSysData[i].y]++;
        psBondAngleID2->_pSysData[i].y = gpu->pOutputBufferCounter[psBondAngleID1->_pSysData[i].z]++;
#if (DUMP_PARAMETERS == 1)
         cout << 
            i << " " << 
            (*psBondAngleID1)[i].x << " " <<
            (*psBondAngleID1)[i].y << " " <<
            (*psBondAngleID1)[i].z << " " <<
            (*psBondAngleID1)[i].w << " " <<
            (*psBondAngleID2)[i].x << " " <<
            (*psBondAngleID2)[i].y << " " <<
            (*psBondAngleParameter)[i].x << " " <<
            (*psBondAngleParameter)[i].y <<
            endl;
#endif
    }
    psBondAngleID1->Upload();
    psBondAngleID2->Upload();
    psBondAngleParameter->Upload();
}

extern "C"
void gpuSetDihedralParameters(gpuContext gpu, const vector<int>& atom1, const vector<int>& atom2, const vector<int>& atom3, const vector<int>& atom4,
        const vector<float>& k, const vector<float>& phase, const vector<int>& periodicity)
{
        int dihedrals = atom1.size();
        gpu->sim.dihedrals = dihedrals;
        CUDAStream<int4>* psDihedralID1             = new CUDAStream<int4>(dihedrals, 1, "DihedralID1");
        gpu->psDihedralID1                          = psDihedralID1;
        gpu->sim.pDihedralID1                       = psDihedralID1->_pDevStream[0];
        CUDAStream<int4>* psDihedralID2             = new CUDAStream<int4>(dihedrals, 1, "DihedralID2");
        gpu->psDihedralID2                          = psDihedralID2;
        gpu->sim.pDihedralID2                       = psDihedralID2->_pDevStream[0];
        CUDAStream<float4>* psDihedralParameter     = new CUDAStream<float4>(dihedrals, 1, "DihedralParameter");
        gpu->psDihedralParameter                    = psDihedralParameter;
        gpu->sim.pDihedralParameter                 = psDihedralParameter->_pDevStream[0];
        for (int i = 0; i < dihedrals; i++)
        {
            (*psDihedralID1)[i].x = atom1[i];
            (*psDihedralID1)[i].y = atom2[i];
            (*psDihedralID1)[i].z = atom3[i];
            (*psDihedralID1)[i].w = atom4[i];
            (*psDihedralParameter)[i].x = k[i];
            (*psDihedralParameter)[i].y = phase[i];
            (*psDihedralParameter)[i].z = (float) periodicity[i];
            psDihedralID2->_pSysData[i].x = gpu->pOutputBufferCounter[psDihedralID1->_pSysData[i].x]++;
            psDihedralID2->_pSysData[i].y = gpu->pOutputBufferCounter[psDihedralID1->_pSysData[i].y]++;
            psDihedralID2->_pSysData[i].z = gpu->pOutputBufferCounter[psDihedralID1->_pSysData[i].z]++;
            psDihedralID2->_pSysData[i].w = gpu->pOutputBufferCounter[psDihedralID1->_pSysData[i].w]++;
#if (DUMP_PARAMETERS == 1)
            cout << 
                i << " " << 
                (*psDihedralID1)[i].x << " " <<
                (*psDihedralID1)[i].y << " " <<
                (*psDihedralID1)[i].z << " " <<
                (*psDihedralID1)[i].w << " " <<
                (*psDihedralID2)[i].x << " " <<
                (*psDihedralID2)[i].y << " " <<
                (*psDihedralID2)[i].z << " " <<
                (*psDihedralID2)[i].w << " " <<
                (*psDihedralParameter)[i].x << " " <<
                (*psDihedralParameter)[i].y << " " <<
                (*psDihedralParameter)[i].z << endl;
#endif
        }
        psDihedralID1->Upload();
        psDihedralID2->Upload();
        psDihedralParameter->Upload();
}

extern "C"
void gpuSetRbDihedralParameters(gpuContext gpu, const vector<int>& atom1, const vector<int>& atom2, const vector<int>& atom3, const vector<int>& atom4,
        const vector<float>& c0, const vector<float>& c1, const vector<float>& c2, const vector<float>& c3, const vector<float>& c4, const vector<float>& c5)
{
    int rb_dihedrals = atom1.size();
    gpu->sim.rb_dihedrals = rb_dihedrals;
    CUDAStream<int4>* psRbDihedralID1           = new CUDAStream<int4>(rb_dihedrals, 1, "RbDihedralID1");
    gpu->psRbDihedralID1                        = psRbDihedralID1;
    gpu->sim.pRbDihedralID1                     = psRbDihedralID1->_pDevStream[0];
    CUDAStream<int4>* psRbDihedralID2           = new CUDAStream<int4>(rb_dihedrals, 1, "RbDihedralID2");
    gpu->psRbDihedralID2                        = psRbDihedralID2;
    gpu->sim.pRbDihedralID2                     = psRbDihedralID2->_pDevStream[0];
    CUDAStream<float4>* psRbDihedralParameter1  = new CUDAStream<float4>(rb_dihedrals, 1, "RbDihedralParameter1");
    gpu->psRbDihedralParameter1                 = psRbDihedralParameter1;
    gpu->sim.pRbDihedralParameter1              = psRbDihedralParameter1->_pDevStream[0];
    CUDAStream<float2>* psRbDihedralParameter2  = new CUDAStream<float2>(rb_dihedrals, 1, "RbDihedralParameter2");
    gpu->psRbDihedralParameter2                 = psRbDihedralParameter2;
    gpu->sim.pRbDihedralParameter2              = psRbDihedralParameter2->_pDevStream[0];

    for (int i = 0; i < rb_dihedrals; i++)
    {
        (*psRbDihedralID1)[i].x = atom1[i];
        (*psRbDihedralID1)[i].y = atom2[i];
        (*psRbDihedralID1)[i].z = atom3[i];
        (*psRbDihedralID1)[i].w = atom4[i];
        (*psRbDihedralParameter1)[i].x = c0[i];
        (*psRbDihedralParameter1)[i].y = c1[i];
        (*psRbDihedralParameter1)[i].z = c2[i];
        (*psRbDihedralParameter1)[i].w = c3[i];
        (*psRbDihedralParameter2)[i].x = c4[i];
        (*psRbDihedralParameter2)[i].y = c5[i];
        psRbDihedralID2->_pSysData[i].x = gpu->pOutputBufferCounter[psRbDihedralID1->_pSysData[i].x]++;
        psRbDihedralID2->_pSysData[i].y = gpu->pOutputBufferCounter[psRbDihedralID1->_pSysData[i].y]++;
        psRbDihedralID2->_pSysData[i].z = gpu->pOutputBufferCounter[psRbDihedralID1->_pSysData[i].z]++;
        psRbDihedralID2->_pSysData[i].w = gpu->pOutputBufferCounter[psRbDihedralID1->_pSysData[i].w]++;
#if (DUMP_PARAMETERS == 1)
        cout << 
            i << " " << 
            (*psRbDihedralID1)[i].x << " " <<
            (*psRbDihedralID1)[i].y << " " <<
            (*psRbDihedralID1)[i].z << " " <<
            (*psRbDihedralID1)[i].w <<" " <<
            (*psRbDihedralID2)[i].x << " " <<
            (*psRbDihedralID2)[i].y << " " <<
            (*psRbDihedralID2)[i].z << " " <<
            (*psRbDihedralID2)[i].w <<" " <<
            (*psRbDihedralParameter1)[i].x << " " <<
            (*psRbDihedralParameter1)[i].y << " " <<
            (*psRbDihedralParameter1)[i].z << " " <<
            (*psRbDihedralParameter1)[i].w << " " <<
            (*psRbDihedralParameter2)[i].x << " " <<
            (*psRbDihedralParameter2)[i].y <<
            endl;
#endif
    }
    psRbDihedralID1->Upload();
    psRbDihedralID2->Upload();
    psRbDihedralParameter1->Upload();
    psRbDihedralParameter2->Upload();
}

extern "C"
void gpuSetLJ14Parameters(gpuContext gpu, float epsfac, float fudge, const vector<int>& atom1, const vector<int>& atom2,
        const vector<float>& c6, const vector<float>& c12, const vector<float>& q1, const vector<float>& q2)
{
    int LJ14s = atom1.size();
    float scale = epsfac * fudge;

    gpu->sim.LJ14s                              = LJ14s;
    CUDAStream<int4>* psLJ14ID                  = new CUDAStream<int4>(LJ14s, 1, "LJ14ID");
    gpu->psLJ14ID                               = psLJ14ID;
    gpu->sim.pLJ14ID                            = psLJ14ID->_pDevStream[0];
    CUDAStream<float4>* psLJ14Parameter         = new CUDAStream<float4>(LJ14s, 1, "LJ14Parameter");
    gpu->psLJ14Parameter                        = psLJ14Parameter;
    gpu->sim.pLJ14Parameter                     = psLJ14Parameter->_pDevStream[0];

    for (int i = 0; i < LJ14s; i++)
    {
        (*psLJ14ID)[i].x = atom1[i];
        (*psLJ14ID)[i].y = atom2[i];
        psLJ14ID->_pSysData[i].z = gpu->pOutputBufferCounter[psLJ14ID->_pSysData[i].x]++;
        psLJ14ID->_pSysData[i].w = gpu->pOutputBufferCounter[psLJ14ID->_pSysData[i].y]++;
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
        (*psLJ14Parameter)[i].x = p0;
        (*psLJ14Parameter)[i].y = p1;
        (*psLJ14Parameter)[i].z = p2;
    }
#if (DUMP_PARAMETERS == 1)
        cout << 
            i << " " <<
            (*psLJ14ID)[i].x << " " <<
            (*psLJ14ID)[i].y << " " <<
            (*psLJ14ID)[i].z << " " <<
            (*psLJ14ID)[i].w << " " <<
            (*psLJ14Parameter)[i].x << " " <<
            (*psLJ14Parameter)[i].y << " " <<
            (*psLJ14Parameter)[i].z << " " <<
            p0 << " " << 
            p1 << " " << 
            p2 << " " << 
            endl;
#endif
    psLJ14ID->Upload();
    psLJ14Parameter->Upload();
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
            (*gpu->psPosq4)[i].w = p0;
            (*gpu->psSigEps2)[i].x = p1;
            (*gpu->psSigEps2)[i].y = p2;
    }

    // Dummy out extra atom data
    for (unsigned int i = coulombs; i < gpu->sim.paddedNumberOfAtoms; i++)
    {
        (*gpu->psPosq4)[i].x       = 100000.0f + i * 10.0f;
        (*gpu->psPosq4)[i].y       = 100000.0f + i * 10.0f;
        (*gpu->psPosq4)[i].z       = 100000.0f + i * 10.0f;
        (*gpu->psPosq4)[i].w       = 0.0f;
        (*gpu->psSigEps2)[i].x     = 0.0f;
        (*gpu->psSigEps2)[i].y     = 0.0f;
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
            (*gpu->psObcData)[i].x = radius[i] - dielectricOffset;
            (*gpu->psObcData)[i].y = scale[i] * (*gpu->psObcData)[i].x;

#if (DUMP_PARAMETERS == 1)
        cout << 
            i << " " << 
            (*gpu->psObcData)[i].x << " " <<
            (*gpu->psObcData)[i].y;
#endif
    }

    // Dummy out extra atom data
    for (unsigned int i = atoms; i < gpu->sim.paddedNumberOfAtoms; i++)
    {
        (*gpu->psBornRadii)[i]     = 0.2f;
        (*gpu->psObcData)[i].x     = 0.01f;
        (*gpu->psObcData)[i].y     = 0.01f;
    }

    gpu->psBornRadii->Upload();
    gpu->psObcData->Upload();
    gpu->sim.preFactor = 2.0f*electricConstant*((1.0f/innerDielectric)-(1.0f/solventDielectric))*gpu->sim.forceConversionFactor;
}

void markShakeClusterInvalid(ShakeCluster& cluster, map<int, ShakeCluster>& allClusters, vector<bool>& invalidForShake)
{
    cluster.valid = false;
    invalidForShake[cluster.centralID] = true;
    for (int i = 0; i < cluster.size; i++) {
        invalidForShake[cluster.peripheralID[i]] = true;
        map<int, ShakeCluster>::iterator otherCluster = allClusters.find(cluster.peripheralID[i]);
        if (otherCluster != allClusters.end() && otherCluster->second.valid)
            markShakeClusterInvalid(otherCluster->second, allClusters, invalidForShake);
    }
}

extern "C"
void gpuSetConstraintParameters(gpuContext gpu, const vector<int>& atom1, const vector<int>& atom2, const vector<float>& distance,
        const vector<float>& invMass1, const vector<float>& invMass2, float shakeTolerance, unsigned int lincsTerms)
{
    // Create a vector for recording which atoms are handled by SHAKE (or SETTLE).

    vector<bool> isShakeAtom(gpu->natoms, false);

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

    CUDAStream<int4>* psSettleID          = new CUDAStream<int4>((int) settleClusters.size(), 1, "SettleID");
    gpu->psSettleID                       = psSettleID;
    gpu->sim.pSettleID                    = psSettleID->_pDevStream[0];
    CUDAStream<float2>* psSettleParameter = new CUDAStream<float2>((int) settleClusters.size(), 1, "SettleParameter");
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
            (*psSettleID)[i].x = atom1;
            (*psSettleID)[i].y = atom2;
            (*psSettleID)[i].z = atom3;
            (*psSettleParameter)[i].x = dist12;
            (*psSettleParameter)[i].y = dist23;
        }
        else if (dist12 == dist23) { // atom2 is the central atom
            (*psSettleID)[i].x = atom2;
            (*psSettleID)[i].y = atom1;
            (*psSettleID)[i].z = atom3;
            (*psSettleParameter)[i].x = dist12;
            (*psSettleParameter)[i].y = dist13;
        }
        else if (dist13 == dist23) { // atom3 is the central atom
            (*psSettleID)[i].x = atom3;
            (*psSettleID)[i].y = atom1;
            (*psSettleID)[i].z = atom2;
            (*psSettleParameter)[i].x = dist13;
            (*psSettleParameter)[i].y = dist12;
        }
        else
            throw OpenMMException("Two of the three distances constrained with SETTLE must be the same.");
        isShakeAtom[atom1] = true;
        isShakeAtom[atom2] = true;
        isShakeAtom[atom3] = true;
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
    vector<bool> invalidForShake(gpu->natoms, false);
    for (int i = 0; i < atom1.size(); i++) {
        if (isShakeAtom[atom1[i]])
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

        // Add it to the cluster.

        if (clusters.find(centralID) == clusters.end()) {
            clusters[centralID] = ShakeCluster(centralID, centralInvMass);
        }
        ShakeCluster& cluster = clusters[centralID];
        cluster.addAtom(peripheralID, distance[i], peripheralInvMass);
        if (constraintCount[peripheralID] != 1 || invalidForShake[atom1[i]] || invalidForShake[atom2[i]]) {
            markShakeClusterInvalid(cluster, clusters, invalidForShake);
            map<int, ShakeCluster>::iterator otherCluster = clusters.find(peripheralID);
            if (otherCluster != clusters.end() && otherCluster->second.valid)
                markShakeClusterInvalid(otherCluster->second, clusters, invalidForShake);
        }
    }
    int validShakeClusters = 0;
    for (map<int, ShakeCluster>::iterator iter = clusters.begin(); iter != clusters.end(); ++iter) {
        ShakeCluster& cluster = iter->second;
        if (cluster.valid) {
            cluster.valid = !invalidForShake[cluster.centralID];
            for (int i = 0; i < cluster.size; i++)
                if (invalidForShake[cluster.peripheralID[i]])
                    cluster.valid = false;
            if (cluster.valid)
                ++validShakeClusters;
        }
    }

    // Fill in the Cuda streams.

    CUDAStream<int4>* psShakeID             = new CUDAStream<int4>(validShakeClusters, 1, "ShakeID");
    gpu->psShakeID                          = psShakeID;
    gpu->sim.pShakeID                       = psShakeID->_pDevStream[0];
    CUDAStream<float4>* psShakeParameter    = new CUDAStream<float4>(validShakeClusters, 1, "ShakeParameter");
    gpu->psShakeParameter                   = psShakeParameter;
    gpu->sim.pShakeParameter                = psShakeParameter->_pDevStream[0];
    gpu->sim.ShakeConstraints               = validShakeClusters;
    int index = 0;
    for (map<int, ShakeCluster>::const_iterator iter = clusters.begin(); iter != clusters.end(); ++iter) {
        const ShakeCluster& cluster = iter->second;
        if (!cluster.valid)
            continue;
        (*psShakeID)[index].x = cluster.centralID;
        (*psShakeID)[index].y = cluster.peripheralID[0];
        (*psShakeID)[index].z = cluster.size > 1 ? cluster.peripheralID[1] : -1;
        (*psShakeID)[index].w = cluster.size > 2 ? cluster.peripheralID[2] : -1;
        (*psShakeParameter)[index].x = cluster.centralInvMass;
        (*psShakeParameter)[index].y = 0.5f/(cluster.centralInvMass+cluster.peripheralInvMass);
        (*psShakeParameter)[index].z = cluster.distance*cluster.distance;
        (*psShakeParameter)[index].w = cluster.peripheralInvMass;
        isShakeAtom[cluster.centralID] = true;
        isShakeAtom[cluster.peripheralID[0]] = true;
        if (cluster.size > 1)
            isShakeAtom[cluster.peripheralID[1]] = true;
        if (cluster.size > 2)
            isShakeAtom[cluster.peripheralID[2]] = true;
        ++index;
    }
    psShakeID->Upload();
    psShakeParameter->Upload();
    gpu->sim.shakeTolerance = shakeTolerance;
    gpu->sim.shake_threads_per_block     = (gpu->sim.ShakeConstraints + gpu->sim.blocks - 1) / gpu->sim.blocks;
    if (gpu->sim.shake_threads_per_block > gpu->sim.max_shake_threads_per_block)
        gpu->sim.shake_threads_per_block = gpu->sim.max_shake_threads_per_block;
    if (gpu->sim.shake_threads_per_block < 1)
        gpu->sim.shake_threads_per_block = 1;

    // Find connected constraints for LINCS.

    vector<int> lincsConstraints;
    for (unsigned int i = 0; i < atom1.size(); i++)
        if (!isShakeAtom[atom1[i]])
            lincsConstraints.push_back(i);
    vector<vector<int> > atomConstraints(gpu->natoms);
    for (unsigned int i = 0; i < lincsConstraints.size(); i++) {
        atomConstraints[atom1[lincsConstraints[i]]].push_back(i);
        atomConstraints[atom2[lincsConstraints[i]]].push_back(i);
    }
    vector<vector<int> > linkedConstraints(lincsConstraints.size());
    for (unsigned int atom = 0; atom < atomConstraints.size(); atom++) {
        for (unsigned int i = 0; i < atomConstraints[atom].size(); i++)
            for (int j = 0; j < i; j++) {
                int c1 = atomConstraints[atom][i];
                int c2 = atomConstraints[atom][j];
                linkedConstraints[c1].push_back(c2);
                linkedConstraints[c2].push_back(c1);
            }
    }
    int totalLinks = 0;
    for (unsigned int i = 0; i < linkedConstraints.size(); i++)
        totalLinks += linkedConstraints[i].size();

    // Fill in the CUDA streams.

    CUDAStream<int2>* psLincsAtoms = new CUDAStream<int2>((int) lincsConstraints.size(), 1, "LincsAtoms");
    gpu->psLincsAtoms              = psLincsAtoms;
    gpu->sim.pLincsAtoms           = psLincsAtoms->_pDevData;
    CUDAStream<float4>* psLincsDistance = new CUDAStream<float4>((int) lincsConstraints.size(), 1, "LincsDistance");
    gpu->psLincsDistance                = psLincsDistance;
    gpu->sim.pLincsDistance             = psLincsDistance->_pDevData;
    CUDAStream<int>* psLincsConnections = new CUDAStream<int>(totalLinks, 1, "LincsConnections");
    gpu->psLincsConnections             = psLincsConnections;
    gpu->sim.pLincsConnections          = psLincsConnections->_pDevData;
    CUDAStream<int>* psLincsConnectionsIndex = new CUDAStream<int>((int) lincsConstraints.size()+1, 1, "LincsConnectionsIndex");
    gpu->psLincsConnectionsIndex             = psLincsConnectionsIndex;
    gpu->sim.pLincsConnectionsIndex          = psLincsConnectionsIndex->_pDevData;
    CUDAStream<int>* psLincsAtomConstraints = new CUDAStream<int>((int) lincsConstraints.size()*2, 1, "LincsAtomConstraints");
    gpu->psLincsAtomConstraints             = psLincsAtomConstraints;
    gpu->sim.pLincsAtomConstraints          = psLincsAtomConstraints->_pDevData;
    CUDAStream<int>* psLincsAtomConstraintsIndex = new CUDAStream<int>(gpu->natoms+1, 1, "LincsAtomConstraintsIndex");
    gpu->psLincsAtomConstraintsIndex             = psLincsAtomConstraintsIndex;
    gpu->sim.pLincsAtomConstraintsIndex          = psLincsAtomConstraintsIndex->_pDevData;
    CUDAStream<float>* psLincsS = new CUDAStream<float>((int) lincsConstraints.size(), 1, "LincsS");
    gpu->psLincsS             = psLincsS;
    gpu->sim.pLincsS          = psLincsS->_pDevData;
    CUDAStream<float>* psLincsCoupling = new CUDAStream<float>(totalLinks, 1, "LincsCoupling");
    gpu->psLincsCoupling               = psLincsCoupling;
    gpu->sim.pLincsCoupling            = psLincsCoupling->_pDevData;
    CUDAStream<float>* psLincsRhs1 = new CUDAStream<float>((int) lincsConstraints.size(), 1, "LincsRhs1");
    gpu->psLincsRhs1             = psLincsRhs1;
    gpu->sim.pLincsRhs1          = psLincsRhs1->_pDevData;
    CUDAStream<float>* psLincsRhs2 = new CUDAStream<float>((int) lincsConstraints.size(), 1, "LincsRhs2");
    gpu->psLincsRhs2             = psLincsRhs2;
    gpu->sim.pLincsRhs2          = psLincsRhs2->_pDevData;
    CUDAStream<float>* psLincsSolution = new CUDAStream<float>((int) lincsConstraints.size(), 1, "LincsSolution");
    gpu->psLincsSolution             = psLincsSolution;
    gpu->sim.pLincsSolution          = psLincsSolution->_pDevData;
    CUDAStream<unsigned int>* psSyncCounter = new CUDAStream<unsigned int>(2*lincsTerms+2, 1, "SyncCounter");
    gpu->psSyncCounter                      = psSyncCounter;
    gpu->sim.pSyncCounter                   = psSyncCounter->_pDevData;
    gpu->sim.lincsConstraints = lincsConstraints.size();
    index = 0;
    for (unsigned int i = 0; i < lincsConstraints.size(); i++) {
        int c = lincsConstraints[i];
        (*psLincsAtoms)[i].x = atom1[c];
        (*psLincsAtoms)[i].y = atom2[c];
        (*psLincsDistance)[i].w = distance[c];
        (*psLincsS)[i] = 1.0f/sqrt(invMass1[c]+invMass2[c]);
        (*psLincsConnectionsIndex)[i] = index;
        for (unsigned int j = 0; j < linkedConstraints[i].size(); j++)
            (*psLincsConnections)[index++] = linkedConstraints[i][j];
    }
    (*psLincsConnectionsIndex)[lincsConstraints.size()] = index;
    for (unsigned int i = 0; i < psSyncCounter->_length; i++)
        (*psSyncCounter)[i] = 0;
    index = 0;
    for (unsigned int i = 0; i < atomConstraints.size(); i++) {
        (*psLincsAtomConstraintsIndex)[i] = index;
        for (unsigned int j = 0; j < atomConstraints[i].size(); j++)
            (*psLincsAtomConstraints)[index++] = atomConstraints[i][j];
    }
    (*psLincsAtomConstraintsIndex)[atomConstraints.size()] = index;
    psLincsAtoms->Upload();
    psLincsDistance->Upload();
    psLincsS->Upload();
    psLincsConnections->Upload();
    psLincsConnectionsIndex->Upload();
    psLincsAtomConstraints->Upload();
    psLincsAtomConstraintsIndex->Upload();
    psSyncCounter->Upload();
    gpu->sim.lincsTerms = lincsTerms;
    gpu->sim.lincs_threads_per_block = (gpu->sim.lincsConstraints + gpu->sim.blocks - 1) / gpu->sim.blocks;
    if (gpu->sim.lincs_threads_per_block > gpu->sim.max_shake_threads_per_block)
        gpu->sim.lincs_threads_per_block = gpu->sim.max_shake_threads_per_block;
    if (gpu->sim.lincs_threads_per_block < 1)
        gpu->sim.lincs_threads_per_block = 1;


    gpu->psLincsConnectionsIndex->Download();
    gpu->psLincsConnections->Download();
    gpu->psLincsCoupling->Download();



#ifdef DeltaShake

    // count number of atoms w/o constraint

    int count = 0;
    for (int i = 0; i < gpu->natoms; i++)
       if (!isShakeAtom[i])
          count++;

    // Allocate NonShake parameters

    gpu->sim.NonShakeConstraints                  = count;
    if( count || true ){

       CUDAStream<int>* psNonShakeID              = new CUDAStream<int>(count, 1, "NonShakeID");
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
          if (!isShakeAtom[i]){
             (*psNonShakeID)[count++] = i;
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
    gpu->psPosq4                        = new CUDAStream<float4>(gpu->sim.paddedNumberOfAtoms, 1, "Posq");
    gpu->sim.stride                     = gpu->psPosq4->_stride;
    gpu->sim.stride2                    = gpu->sim.stride * 2;
    gpu->sim.stride3                    = gpu->sim.stride * 3;
    gpu->sim.stride4                    = gpu->sim.stride * 4;
    gpu->sim.pPosq                      = gpu->psPosq4->_pDevStream[0];
    gpu->sim.stride                     = gpu->psPosq4->_stride;
    gpu->sim.stride2                    = 2 * gpu->sim.stride;
    gpu->sim.stride3                    = 3 * gpu->sim.stride;
    gpu->sim.stride4                    = 4 * gpu->sim.stride;
    gpu->psPosqP4                       = new CUDAStream<float4>(gpu->sim.paddedNumberOfAtoms, 1, "PosqP");
    gpu->sim.pPosqP                     = gpu->psPosqP4->_pDevStream[0];
    gpu->psOldPosq4                     = new CUDAStream<float4>(gpu->sim.paddedNumberOfAtoms, 1, "OldPosq");
    gpu->sim.pOldPosq                   = gpu->psOldPosq4->_pDevStream[0];
    gpu->psVelm4                        = new CUDAStream<float4>(gpu->sim.paddedNumberOfAtoms, 1, "Velm");
    gpu->sim.pVelm4                     = gpu->psVelm4->_pDevStream[0];
    gpu->psvVector4                     = new CUDAStream<float4>(gpu->sim.paddedNumberOfAtoms, 1, "vVector");
    gpu->sim.pvVector4                  = gpu->psvVector4->_pDevStream[0];
    gpu->psxVector4                     = new CUDAStream<float4>(gpu->sim.paddedNumberOfAtoms, 1, "xVector");
    gpu->sim.pxVector4                  = gpu->psxVector4->_pDevStream[0];
    gpu->psBornRadii                    = new CUDAStream<float>(gpu->sim.paddedNumberOfAtoms, 1, "BornRadii");
    gpu->sim.pBornRadii                 = gpu->psBornRadii->_pDevStream[0];
    gpu->psObcChain                     = new CUDAStream<float>(gpu->sim.paddedNumberOfAtoms, 1, "ObcChain");
    gpu->sim.pObcChain                  = gpu->psObcChain->_pDevStream[0];
    gpu->psSigEps2                      = new CUDAStream<float2>(gpu->sim.paddedNumberOfAtoms, 1, "SigEps2");
    gpu->sim.pAttr                      = gpu->psSigEps2->_pDevStream[0];
    gpu->psObcData                      = new CUDAStream<float2>(gpu->sim.paddedNumberOfAtoms, 1, "ObcData");
    gpu->sim.pObcData                   = gpu->psObcData->_pDevStream[0];
    gpu->pAtomSymbol                    = new unsigned char[gpu->natoms];
    gpu->psAtomIndex                    = new CUDAStream<int>(gpu->sim.paddedNumberOfAtoms, 1, "AtomIndex");
    gpu->sim.pAtomIndex                 = gpu->psAtomIndex->_pDevStream[0];
    for (int i = 0; i < (int) gpu->sim.paddedNumberOfAtoms; i++)
        (*gpu->psAtomIndex)[i] = i;
    gpu->psAtomIndex->Upload();
    // Determine randoms
    gpu->seed                           = 1;
    gpu->sim.randomFrames               = 20;
    gpu->sim.randomIterations           = gpu->sim.randomFrames;
    gpu->sim.randoms                    = gpu->sim.randomFrames * gpu->sim.paddedNumberOfAtoms - 5 * GRID;
    gpu->sim.totalRandoms               = gpu->sim.randoms + gpu->sim.paddedNumberOfAtoms;
    gpu->sim.totalRandomsTimesTwo       = gpu->sim.totalRandoms * 2;
    gpu->psRandom4                      = new CUDAStream<float4>(gpu->sim.totalRandomsTimesTwo, 1, "Random4");
    gpu->psRandom2                      = new CUDAStream<float2>(gpu->sim.totalRandomsTimesTwo, 1, "Random2");
    gpu->psRandomPosition               = new CUDAStream<int>(gpu->sim.blocks, 1, "RandomPosition");
    gpu->psRandomSeed                   = new CUDAStream<uint4>(gpu->sim.blocks * gpu->sim.random_threads_per_block, 1, "RandomSeed");
    gpu->sim.pRandom4a                  = gpu->psRandom4->_pDevStream[0];
    gpu->sim.pRandom2a                  = gpu->psRandom2->_pDevStream[0];
    gpu->sim.pRandom4b                  = gpu->psRandom4->_pDevStream[0] + gpu->sim.totalRandoms;
    gpu->sim.pRandom2b                  = gpu->psRandom2->_pDevStream[0] + gpu->sim.totalRandoms;
    gpu->sim.pRandomPosition            = gpu->psRandomPosition->_pDevStream[0];
    gpu->sim.pRandomSeed                = gpu->psRandomSeed->_pDevStream[0];

    // Allocate and clear linear momentum buffer
    gpu->psLinearMomentum = new CUDAStream<float4>(gpu->sim.blocks, 1, "LinearMomentum");
    gpu->sim.pLinearMomentum = gpu->psLinearMomentum->_pDevStream[0];
    for (int i = 0; i < (int) gpu->sim.blocks; i++)
    {
        (*gpu->psLinearMomentum)[i].x = 0.0f;
        (*gpu->psLinearMomentum)[i].y = 0.0f;
        (*gpu->psLinearMomentum)[i].z = 0.0f;
        (*gpu->psLinearMomentum)[i].w = 0.0f;
    }
    gpu->psLinearMomentum->Upload();

    return 1;
}

extern "C"
void gpuSetPositions(gpuContext gpu, const vector<float>& x, const vector<float>& y, const vector<float>& z)
{
    for (int i = 0; i < gpu->natoms; i++)
    {
        (*gpu->psPosq4)[i].x = x[i];
        (*gpu->psPosq4)[i].y = y[i];
        (*gpu->psPosq4)[i].z = z[i];
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
        (*gpu->psVelm4)[i].x = x[i];
        (*gpu->psVelm4)[i].y = y[i];
        (*gpu->psVelm4)[i].z = z[i];
    }
    gpu->psVelm4->Upload();
} 

extern "C"
void gpuSetMass(gpuContext gpu, const vector<float>& mass)
{
    float totalMass = 0.0f;
    for (int i = 0; i < gpu->natoms; i++)
    {
        (*gpu->psVelm4)[i].w = 1.0f/mass[i];
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
        (*gpu->psRandomPosition)[i] = 0;
    }
    int seed = gpu->seed | ((gpu->seed ^ 0xffffffff) << 16);
    srand(seed);
    for (int i = 0; i < (int) (gpu->sim.blocks * gpu->sim.random_threads_per_block); i++)
    {
        (*gpu->psRandomSeed)[i].x = rand();
        (*gpu->psRandomSeed)[i].y = rand();
        (*gpu->psRandomSeed)[i].z = rand();
        (*gpu->psRandomSeed)[i].w = rand();
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
        (*gpu->psxVector4)[i].x = 0.0f;
        (*gpu->psxVector4)[i].y = 0.0f;
        (*gpu->psxVector4)[i].z = 0.0f;
        (*gpu->psxVector4)[i].w = 0.0f;
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
    gpu->psExclusionIndex           = NULL;
    gpu->psWorkUnit                 = NULL;
    gpu->psInteractingWorkUnit      = NULL;
    gpu->psInteractionFlag          = NULL;
    gpu->psInteractionCount         = NULL;
    gpu->psGridBoundingBox          = NULL;
    gpu->psGridCenter               = NULL;
    gpu->psLincsAtoms               = NULL;
    gpu->psLincsDistance            = NULL;
    gpu->psLincsConnections         = NULL;
    gpu->psLincsConnectionsIndex    = NULL;
    gpu->psLincsAtomConstraints     = NULL;
    gpu->psLincsAtomConstraintsIndex= NULL;
    gpu->psLincsS                   = NULL;
    gpu->psLincsCoupling            = NULL;
    gpu->psLincsRhs1                = NULL;
    gpu->psLincsRhs2                = NULL;
    gpu->psLincsSolution            = NULL;
    gpu->psSyncCounter              = NULL;

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
    delete gpu->psExclusionIndex;
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
    delete gpu->psLincsAtoms;
    delete gpu->psLincsDistance;
    delete gpu->psLincsConnections;
    delete gpu->psLincsConnectionsIndex;
    delete gpu->psLincsAtomConstraints;
    delete gpu->psLincsAtomConstraintsIndex;
    delete gpu->psLincsS;
    delete gpu->psLincsCoupling;
    delete gpu->psLincsRhs1;
    delete gpu->psLincsRhs2;
    delete gpu->psLincsSolution;
    delete gpu->psSyncCounter;
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
    gpu->psForce4               = new CUDAStream<float4>(gpu->sim.paddedNumberOfAtoms, outputBuffers, "Force");
    gpu->psBornForce            = new CUDAStream<float>(gpu->sim.paddedNumberOfAtoms, gpu->sim.nonbondOutputBuffers, "BornForce");
    gpu->psBornSum              = new CUDAStream<float>(gpu->sim.paddedNumberOfAtoms, gpu->sim.nonbondOutputBuffers, "BornSum");
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
        (*gpu->psBondID)[i].z = flip - (*gpu->psBondID)[i].z;
        (*gpu->psBondID)[i].w = flip - (*gpu->psBondID)[i].w;
    }
    for (int i = 0; i < (int) gpu->sim.bond_angles; i++)
    {
        (*gpu->psBondAngleID1)[i].w = flip - (*gpu->psBondAngleID1)[i].w;
        (*gpu->psBondAngleID2)[i].x = flip - (*gpu->psBondAngleID2)[i].x;
        (*gpu->psBondAngleID2)[i].y = flip - (*gpu->psBondAngleID2)[i].y;
    }
    for (int i = 0; i < (int) gpu->sim.dihedrals; i++)
    {
        (*gpu->psDihedralID2)[i].x = flip - (*gpu->psDihedralID2)[i].x;
        (*gpu->psDihedralID2)[i].y = flip - (*gpu->psDihedralID2)[i].y;
        (*gpu->psDihedralID2)[i].z = flip - (*gpu->psDihedralID2)[i].z;
        (*gpu->psDihedralID2)[i].w = flip - (*gpu->psDihedralID2)[i].w;
    }
    for (int i = 0; i < (int) gpu->sim.rb_dihedrals; i++)
    {
        (*gpu->psRbDihedralID2)[i].x = flip - (*gpu->psRbDihedralID2)[i].x;
        (*gpu->psRbDihedralID2)[i].y = flip - (*gpu->psRbDihedralID2)[i].y;
        (*gpu->psRbDihedralID2)[i].z = flip - (*gpu->psRbDihedralID2)[i].z;
        (*gpu->psRbDihedralID2)[i].w = flip - (*gpu->psRbDihedralID2)[i].w;
    }
    for (int i = 0; i < (int) gpu->sim.LJ14s; i++)
    {
        (*gpu->psLJ14ID)[i].z = flip - (*gpu->psLJ14ID)[i].z;
        (*gpu->psLJ14ID)[i].w = flip - (*gpu->psLJ14ID)[i].w;
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
    CUDAStream<unsigned int>* psWorkUnit = new CUDAStream<unsigned int>(cells, 1u, "WorkUnit");
    unsigned int* pWorkList = psWorkUnit->_pSysData;
    gpu->psWorkUnit = psWorkUnit;
    gpu->sim.pWorkUnit = psWorkUnit->_pDevStream[0];
    CUDAStream<unsigned int>* psInteractingWorkUnit = new CUDAStream<unsigned int>(cells, 1u, "InteractingWorkUnit");
    gpu->psInteractingWorkUnit = psInteractingWorkUnit;
    gpu->sim.pInteractingWorkUnit = psInteractingWorkUnit->_pDevStream[0];
    CUDAStream<unsigned int>* psInteractionFlag = new CUDAStream<unsigned int>(cells, 1u, "InteractionFlag");
    gpu->psInteractionFlag = psInteractionFlag;
    gpu->sim.pInteractionFlag = psInteractionFlag->_pDevStream[0];
    CUDAStream<size_t>* psInteractionCount = new CUDAStream<size_t>(1, 1u, "InteractionCount");
    gpu->psInteractionCount = psInteractionCount;
    gpu->sim.pInteractionCount = psInteractionCount->_pDevStream[0];
    CUDAStream<float4>* psGridBoundingBox = new CUDAStream<float4>(dim, 1u, "GridBoundingBox");
    gpu->psGridBoundingBox = psGridBoundingBox;
    gpu->sim.pGridBoundingBox = psGridBoundingBox->_pDevStream[0];
    CUDAStream<float4>* psGridCenter = new CUDAStream<float4>(dim, 1u, "GridCenter");
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
    const unsigned int dim = atoms/grid;
    unsigned int* pWorkList = gpu->psWorkUnit->_pSysData;

    // Mark which work units have exclusions.

    for (int atom1 = 0; atom1 < gpu->exclusions.size(); ++atom1)
    {
        int x = atom1/grid;
        for (int j = 0; j < gpu->exclusions[atom1].size(); ++j)
        {
            int atom2 = gpu->exclusions[atom1][j];
            int y = atom2/grid;
            int cell = (x > y ? x+y*dim-y*(y+1)/2 : y+x*dim-x*(x+1)/2);
            pWorkList[cell] |= 1;
        }
    }
    if (gpu->sim.paddedNumberOfAtoms > gpu->natoms)
    {
        int lastBlock = gpu->natoms/grid;
        for (int i = 0; i < gpu->sim.workUnits; ++i)
        {
            int x = pWorkList[i]>>17;
            int y = (pWorkList[i]>>2)&0x7FFF;
            if (x == lastBlock || y == lastBlock)
                pWorkList[i] |= 1;
        }
    }

    // Build a list of indexes for the work units with exclusions.

    CUDAStream<unsigned int>* psExclusionIndex = new CUDAStream<unsigned int>(gpu->sim.workUnits, 1u, "ExclusionIndex");
    gpu->psExclusionIndex = psExclusionIndex;
    unsigned int* pExclusionIndex = psExclusionIndex->_pSysData;
    gpu->sim.pExclusionIndex = psExclusionIndex->_pDevData;
    int numWithExclusions = 0;
    for (int i = 0; i < psExclusionIndex->_length; ++i)
        if ((pWorkList[i]&1) == 1)
            pExclusionIndex[i] = (numWithExclusions++)*grid;

    // Record the exclusion data.

    CUDAStream<unsigned int>* psExclusion = new CUDAStream<unsigned int>(numWithExclusions*grid, 1u, "Exclusion");
    gpu->psExclusion = psExclusion;
    unsigned int* pExclusion = psExclusion->_pSysData;
    gpu->sim.pExclusion = psExclusion->_pDevData;
    for (int i = 0; i < psExclusion->_length; ++i)
        pExclusion[i] = 0xFFFFFFFF;
    for (int atom1 = 0; atom1 < gpu->exclusions.size(); ++atom1)
    {
        int x = atom1/grid;
        int offset1 = atom1-x*grid;
        for (int j = 0; j < gpu->exclusions[atom1].size(); ++j)
        {
            int atom2 = gpu->exclusions[atom1][j];
            int y = atom2/grid;
            int offset2 = atom2-y*grid;
            if (x > y)
            {
                int cell = x+y*dim-y*(y+1)/2;
                pExclusion[pExclusionIndex[cell]+offset1] &= 0xFFFFFFFF-(1<<offset2);
            }
            else
            {
                int cell = y+x*dim-x*(x+1)/2;
                pExclusion[pExclusionIndex[cell]+offset2] &= 0xFFFFFFFF-(1<<offset1);
            }
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
            int offset2 = atom2-y*grid;
            if (x >= y)
            {
                int cell = x+y*dim-y*(y+1)/2;
                pExclusion[pExclusionIndex[cell]+offset1] &= 0xFFFFFFFF-(1<<offset2);
            }
            if (y >= x)
            {
                int cell = y+x*dim-x*(x+1)/2;
                pExclusion[pExclusionIndex[cell]+offset2] &= 0xFFFFFFFF-(1<<offset1);
            }
        }
    }
    
    psExclusion->Upload();
    psExclusionIndex->Upload();
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
    SetCalculateObcGbsaForces2Sim(gpu);
    SetCalculateAndersenThermostatSim(gpu);
    SetForcesSim(gpu);
    SetUpdateShakeHSim(gpu);
    SetVerletUpdateSim(gpu);
    SetBrownianUpdateSim(gpu);
    SetSettleSim(gpu);
    SetLincsSim(gpu);
    SetRandomSim(gpu);
    return 1;
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
        int atom1 = (*gpu->psShakeID)[i].x;
        int atom2 = (*gpu->psShakeID)[i].y;
        int atom3 = (*gpu->psShakeID)[i].z;
        int atom4 = (*gpu->psShakeID)[i].w;
        float distance2 = (*gpu->psShakeParameter)[i].z;
        constraints.push_back(Constraint(atom1, atom2, distance2));
        if (atom3 != -1)
            constraints.push_back(Constraint(atom1, atom3, distance2));
        if (atom4 != -1)
            constraints.push_back(Constraint(atom1, atom4, distance2));
    }
    for (int i = 0; i < gpu->sim.settleConstraints; i++)
    {
        int atom1 = (*gpu->psSettleID)[i].x;
        int atom2 = (*gpu->psSettleID)[i].y;
        int atom3 = (*gpu->psSettleID)[i].z;
        float distance12 = (*gpu->psSettleParameter)[i].x;
        float distance23 = (*gpu->psSettleParameter)[i].y;
        constraints.push_back(Constraint(atom1, atom2, distance12*distance12));
        constraints.push_back(Constraint(atom1, atom3, distance12*distance12));
        constraints.push_back(Constraint(atom2, atom3, distance23*distance23));
    }
    for (int i = 0; i < gpu->sim.lincsConstraints; i++)
    {
        int atom1 = (*gpu->psLincsAtoms)[i].x;
        int atom2 = (*gpu->psLincsAtoms)[i].y;
        float distance2 = (*gpu->psLincsDistance)[i].w;
        constraints.push_back(Constraint(atom1, atom2, distance2));
    }

    // First make a list of every other atom to which each atom is connect by a bond or constraint.

    int numAtoms = gpu->natoms;
    vector<vector<int> > atomBonds(numAtoms);
    for (int i = 0; i < gpu->sim.bonds; i++)
    {
        int atom1 = (*gpu->psBondID)[i].x;
        int atom2 = (*gpu->psBondID)[i].y;
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
        int atom1 = (*gpu->psBondID)[i].x;
        molecules[atomMolecule[atom1]].bonds.push_back(i);
    }
    for (int i = 0; i < gpu->sim.bond_angles; i++)
    {
        int atom1 = (*gpu->psBondAngleID1)[i].x;
        molecules[atomMolecule[atom1]].angles.push_back(i);
    }
    for (int i = 0; i < gpu->sim.dihedrals; i++)
    {
        int atom1 = (*gpu->psDihedralID1)[i].x;
        molecules[atomMolecule[atom1]].periodicTorsions.push_back(i);
    }
    for (int i = 0; i < gpu->sim.rb_dihedrals; i++)
    {
        int atom1 = (*gpu->psRbDihedralID1)[i].x;
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
                float dx = floor(molPos[i].x/gpu->sim.periodicBoxSizeX)*gpu->sim.periodicBoxSizeX;
                float dy = floor(molPos[i].y/gpu->sim.periodicBoxSizeY)*gpu->sim.periodicBoxSizeY;
                float dz = floor(molPos[i].z/gpu->sim.periodicBoxSizeZ)*gpu->sim.periodicBoxSizeZ;
                if (dx != 0.0f || dy != 0.0f || dz != 0.0f)
                {
                    molPos[i].x -= dx;
                    molPos[i].y -= dy;
                    molPos[i].z -= dz;
                    for (int j = 0; j < atoms.size(); j++)
                    {
                        int atom = atoms[j]+mol.instances[i];
                        posq[atom].x -= dx;
                        posq[atom].y -= dy;
                        posq[atom].z -= dz;
                    }
                }
            }
        }

        // Select a bin for each molecule, then sort them by bin.

        bool useHilbert = (numMolecules > 5000 || atoms.size() > 8); // For small systems, a simple zigzag curve works better than a Hilbert curve.
        float binWidth;
        if (useHilbert)
            binWidth = max(max(maxx-minx, maxy-miny), maxz-minz)/255.0;
        else
            binWidth = 0.2*sqrt(gpu->sim.nonbondedCutoffSqr);
        int xbins = 1 + (int) ((maxx-minx)/binWidth);
        int ybins = 1 + (int) ((maxy-miny)/binWidth);
        vector<pair<int, int> > molBins(numMolecules);
        bitmask_t coords[3];
        for (int i = 0; i < numMolecules; i++)
        {
            int x = (int) ((molPos[i].x-minx)/binWidth);
            int y = (int) ((molPos[i].y-miny)/binWidth);
            int z = (int) ((molPos[i].z-minz)/binWidth);
            int bin;
            if (useHilbert)
            {
                coords[0] = x;
                coords[1] = y;
                coords[2] = z;
                bin = (int) hilbert_c2i(3, 8, coords);
            }
            else
            {
                int yodd = y&1;
                int zodd = z&1;
                bin = z*xbins*ybins;
                bin += (zodd ? ybins-y : y)*xbins;
                bin += (yodd ? xbins-x : x);
            }
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
                originalIndex[newIndex] = (*gpu->psAtomIndex)[oldIndex];
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
        (*gpu->psAtomIndex)[i] = originalIndex[i];
    gpu->psAtomIndex->Upload();
}
