/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
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

/**
 * This tests the CPU implementation of PME.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "openmm/NonbondedForce.h"
#include "openmm/internal/NonbondedForceImpl.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/Units.h"
#include "../src/CpuPmeKernels.h"
#include "SimTKOpenMMRealType.h"
#include "sfmt/SFMT.h"
#include <iostream>
#include <vector>

using namespace OpenMM;
using namespace std;

class IO : public CalcPmeReciprocalForceKernel::IO {
public:
    vector<float> posq;
    float* force;
    float* getPosq() {
        return &posq[0];
    }
    void setForce(float* force) {
        this->force = force;
    }
};

void make_waterbox(int natoms, double boxEdgeLength, NonbondedForce *forceField,  vector<Vec3> &positions, vector<double>& eps, vector<double>& sig,
                   vector<pair<int, int> >& bonds, System &system, bool do_electrostatics) {
    const int RESSIZE = 3;
    const double masses[RESSIZE]    = {     8.0,    1.0,    1.0 };
    const double charges[RESSIZE]   = {  -0.834,  0.417,  0.417 };
    // Values from the CHARMM force field, in AKMA units
    const double epsilons[RESSIZE]  = { -0.1521, -0.046, -0.046 };
    const double halfrmins[RESSIZE] = {  1.7682, 0.2245, 0.2245 };
    positions.clear();
    if(natoms == 6){
        const double coords[6][3] = {
            {  2.000000, 2.000000, 2.000000},
            {  2.500000, 2.000000, 3.000000},
            {  1.500000, 2.000000, 3.000000},
            {  0.000000, 0.000000, 0.000000},
            {  0.500000, 0.000000, 1.000000},
            { -0.500000, 0.000000, 1.000000}
        };
        for(int atom = 0; atom < natoms; ++atom)
            positions.push_back(Vec3(coords[atom][0], coords[atom][1], coords[atom][2])*OpenMM::NmPerAngstrom);
    }else if(natoms == 375){
        const double coords[375][3] = {
            { -6.22, -6.25, -6.24 },
            { -5.32, -6.03, -6.00 },
            { -6.75, -5.56, -5.84 },
            { -3.04, -6.23, -6.19 },
            { -3.52, -5.55, -5.71 },
            { -3.59, -6.43, -6.94 },
            {  0.02, -6.23, -6.14 },
            { -0.87, -5.97, -6.37 },
            {  0.53, -6.03, -6.93 },
            {  3.10, -6.20, -6.27 },
            {  3.87, -6.35, -5.72 },
            {  2.37, -6.11, -5.64 },
            {  6.18, -6.14, -6.20 },
            {  6.46, -6.66, -5.44 },
            {  6.26, -6.74, -6.94 },
            { -6.21, -3.15, -6.24 },
            { -6.23, -3.07, -5.28 },
            { -6.02, -2.26, -6.55 },
            { -3.14, -3.07, -6.16 },
            { -3.38, -3.63, -6.90 },
            { -2.18, -3.05, -6.17 },
            { -0.00, -3.16, -6.23 },
            { -0.03, -2.30, -6.67 },
            {  0.05, -2.95, -5.29 },
            {  3.08, -3.11, -6.14 },
            {  2.65, -2.55, -6.79 },
            {  3.80, -3.53, -6.62 },
            {  6.16, -3.14, -6.16 },
            {  7.04, -3.32, -6.51 },
            {  5.95, -2.27, -6.51 },
            { -6.20, -0.04, -6.15 },
            { -5.43,  0.32, -6.59 },
            { -6.95,  0.33, -6.62 },
            { -3.10, -0.06, -6.19 },
            { -3.75,  0.42, -6.69 },
            { -2.46,  0.60, -5.93 },
            {  0.05, -0.01, -6.17 },
            { -0.10,  0.02, -7.12 },
            { -0.79,  0.16, -5.77 },
            {  3.03,  0.00, -6.19 },
            {  3.54,  0.08, -7.01 },
            {  3.69, -0.22, -5.53 },
            {  6.17,  0.05, -6.19 },
            {  5.78, -0.73, -6.57 },
            {  7.09, -0.17, -6.04 },
            { -6.20,  3.15, -6.25 },
            { -6.59,  3.18, -5.37 },
            { -5.87,  2.25, -6.33 },
            { -3.09,  3.04, -6.17 },
            { -3.88,  3.58, -6.26 },
            { -2.41,  3.54, -6.63 },
            {  0.00,  3.06, -6.26 },
            { -0.71,  3.64, -6.00 },
            {  0.65,  3.15, -5.55 },
            {  3.14,  3.06, -6.23 },
            {  3.11,  3.31, -5.30 },
            {  2.38,  3.49, -6.63 },
            {  6.19,  3.14, -6.25 },
            {  6.82,  3.25, -5.54 },
            {  5.76,  2.30, -6.07 },
            { -6.22,  6.26, -6.19 },
            { -6.22,  5.74, -7.00 },
            { -5.89,  5.67, -5.52 },
            { -3.04,  6.24, -6.20 },
            { -3.08,  5.28, -6.17 },
            { -3.96,  6.52, -6.25 },
            { -0.05,  6.21, -6.16 },
            {  0.82,  6.58, -6.06 },
            {  0.01,  5.64, -6.93 },
            {  3.10,  6.25, -6.15 },
            {  3.64,  5.47, -6.31 },
            {  2.46,  6.24, -6.87 },
            {  6.22,  6.20, -6.27 },
            {  5.37,  6.42, -5.88 },
            {  6.80,  6.07, -5.51 },
            { -6.19, -6.15, -3.13 },
            { -6.37, -7.01, -3.51 },
            { -6.25, -6.29, -2.18 },
            { -3.10, -6.27, -3.11 },
            { -2.29, -5.77, -2.99 },
            { -3.80, -5.62, -2.98 },
            { -0.03, -6.18, -3.15 },
            { -0.07, -7.05, -2.75 },
            {  0.68, -5.74, -2.70 },
            {  3.10, -6.14, -3.07 },
            {  2.35, -6.72, -3.23 },
            {  3.86, -6.65, -3.37 },
            {  6.22, -6.20, -3.16 },
            {  6.82, -6.36, -2.43 },
            {  5.35, -6.13, -2.75 },
            { -6.26, -3.13, -3.12 },
            { -6.16, -2.27, -2.70 },
            { -5.36, -3.47, -3.18 },
            { -3.11, -3.05, -3.14 },
            { -3.31, -3.96, -3.34 },
            { -2.77, -3.06, -2.24 },
            {  0.00, -3.13, -3.16 },
            {  0.48, -2.37, -2.81 },
            { -0.57, -3.40, -2.44 },
            {  3.09, -3.09, -3.16 },
            {  2.41, -3.19, -2.49 },
            {  3.91, -3.07, -2.67 },
            {  6.19, -3.04, -3.08 },
            {  5.64, -3.61, -3.61 },
            {  6.93, -3.58, -2.82 },
            { -6.18, -0.00, -3.04 },
            { -6.00, -0.59, -3.78 },
            { -6.79,  0.64, -3.39 },
            { -3.05, -0.03, -3.07 },
            { -2.95,  0.80, -3.52 },
            { -4.00, -0.20, -3.07 },
            { -0.03,  0.03, -3.06 },
            { -0.33, -0.37, -3.87 },
            {  0.89, -0.21, -2.99 },
            {  3.13, -0.05, -3.10 },
            {  3.44,  0.81, -3.34 },
            {  2.21,  0.07, -2.86 },
            {  6.20, -0.05, -3.13 },
            {  6.89,  0.60, -3.20 },
            {  5.58,  0.30, -2.49 },
            { -6.23,  3.09, -3.16 },
            { -5.62,  3.79, -2.94 },
            { -6.33,  2.60, -2.33 },
            { -3.10,  3.08, -3.04 },
            { -3.84,  3.47, -3.51 },
            { -2.40,  3.01, -3.69 },
            {  0.01,  3.04, -3.11 },
            { -0.56,  3.59, -3.64 },
            {  0.28,  3.60, -2.38 },
            {  3.04,  3.11, -3.09 },
            {  3.49,  2.30, -2.87 },
            {  3.70,  3.66, -3.51 },
            {  6.15,  3.14, -3.11 },
            {  6.52,  2.52, -3.74 },
            {  6.72,  3.06, -2.34 },
            { -6.22,  6.15, -3.13 },
            { -5.49,  6.21, -2.51 },
            { -6.56,  7.04, -3.18 },
            { -3.11,  6.24, -3.05 },
            { -3.76,  5.83, -3.62 },
            { -2.26,  5.92, -3.37 },
            {  0.03,  6.25, -3.07 },
            {  0.34,  5.63, -3.73 },
            { -0.87,  6.00, -2.91 },
            {  3.07,  6.15, -3.08 },
            {  3.29,  6.92, -2.56 },
            {  3.39,  6.35, -3.96 },
            {  6.22,  6.14, -3.12 },
            {  5.79,  6.38, -2.29 },
            {  6.25,  6.96, -3.62 },
            { -6.21, -6.20, -0.06 },
            { -5.79, -6.87,  0.48 },
            { -6.43, -5.50,  0.54 },
            { -3.16, -6.21, -0.02 },
            { -2.50, -6.87,  0.20 },
            { -2.77, -5.37,  0.23 },
            { -0.00, -6.14, -0.00 },
            {  0.68, -6.72, -0.33 },
            { -0.64, -6.73,  0.38 },
            {  3.03, -6.20, -0.01 },
            {  3.77, -6.56, -0.51 },
            {  3.43, -5.85,  0.78 },
            {  6.25, -6.16, -0.00 },
            {  5.36, -6.09, -0.36 },
            {  6.24, -6.97,  0.49 },
            { -6.24, -3.05, -0.01 },
            { -6.35, -3.64,  0.73 },
            { -5.42, -3.33, -0.42 },
            { -3.09, -3.06,  0.05 },
            { -2.44, -3.62, -0.38 },
            { -3.90, -3.21, -0.43 },
            {  0.05, -3.10,  0.02 },
            { -0.31, -2.35, -0.43 },
            { -0.63, -3.77,  0.01 },
            {  3.05, -3.09, -0.04 },
            {  3.28, -3.90,  0.41 },
            {  3.65, -2.43,  0.30 },
            {  6.20, -3.04, -0.03 },
            {  5.66, -3.31,  0.71 },
            {  6.78, -3.79, -0.19 },
            { -6.18,  0.04, -0.04 },
            { -6.73, -0.73, -0.15 },
            { -5.98,  0.06,  0.89 },
            { -3.11, -0.04, -0.04 },
            { -3.36, -0.08,  0.87 },
            { -2.70,  0.81, -0.14 },
            { -0.02, -0.02, -0.05 },
            { -0.45,  0.28,  0.75 },
            {  0.90,  0.15,  0.07 },
            {  3.04,  0.02, -0.01 },
            {  3.26, -0.82,  0.38 },
            {  3.89,  0.45, -0.13 },
            {  6.19,  0.05, -0.03 },
            {  5.52, -0.56,  0.25 },
            {  7.01, -0.29,  0.32 },
            { -6.14,  3.08,  0.00 },
            { -6.83,  2.82,  0.61 },
            { -6.59,  3.64, -0.64 },
            { -3.05,  3.09, -0.04 },
            { -3.79,  2.50,  0.09 },
            { -3.18,  3.80,  0.59 },
            {  0.02,  3.14,  0.04 },
            { -0.89,  3.04, -0.19 },
            {  0.49,  2.57, -0.57 },
            {  3.14,  3.15,  0.00 },
            {  3.28,  2.28,  0.37 },
            {  2.30,  3.08, -0.45 },
            {  6.27,  3.08, -0.00 },
            {  5.55,  2.54, -0.33 },
            {  5.83,  3.87,  0.34 },
            { -6.18,  6.15, -0.03 },
            { -6.45,  6.21,  0.88 },
            { -6.26,  7.05, -0.36 },
            { -3.06,  6.19, -0.05 },
            { -2.84,  6.64,  0.76 },
            { -3.99,  5.96,  0.03 },
            { -0.00,  6.20,  0.06 },
            { -0.67,  5.99, -0.59 },
            {  0.76,  6.46, -0.44 },
            {  3.10,  6.26, -0.03 },
            {  3.57,  6.09,  0.78 },
            {  2.57,  5.47, -0.18 },
            {  6.26,  6.18,  0.02 },
            {  5.53,  5.64, -0.29 },
            {  5.95,  7.08, -0.06 },
            { -6.26, -6.21,  3.07 },
            { -5.98, -6.38,  3.97 },
            { -5.46, -5.94,  2.62 },
            { -3.10, -6.24,  3.04 },
            { -2.69, -6.51,  3.87 },
            { -3.43, -5.35,  3.21 },
            { -0.03, -6.16,  3.06 },
            {  0.83, -6.00,  3.42 },
            { -0.30, -6.99,  3.45 },
            {  3.15, -6.25,  3.11 },
            {  2.77, -5.60,  3.72 },
            {  2.68, -6.10,  2.28 },
            {  6.20, -6.21,  3.16 },
            {  5.75, -6.73,  2.50 },
            {  6.69, -5.56,  2.66 },
            { -6.17, -3.10,  3.04 },
            { -6.82, -2.44,  3.28 },
            { -6.12, -3.69,  3.80 },
            { -3.08, -3.04,  3.11 },
            { -3.59, -3.56,  3.72 },
            { -2.97, -3.61,  2.34 },
            {  0.01, -3.04,  3.11 },
            { -0.86, -3.41,  3.20 },
            {  0.56, -3.78,  2.86 },
            {  3.07, -3.07,  3.15 },
            {  3.81, -3.68,  3.13 },
            {  2.80, -2.98,  2.23 },
            {  6.20, -3.04,  3.13 },
            {  5.48, -3.64,  2.92 },
            {  6.98, -3.49,  2.81 },
            { -6.18, -0.05,  3.12 },
            { -6.41,  0.66,  3.69 },
            { -6.33,  0.28,  2.23 },
            { -3.05,  0.03,  3.10 },
            { -3.46, -0.42,  3.83 },
            { -3.57, -0.19,  2.33 },
            {  0.03, -0.02,  3.15 },
            {  0.23, -0.08,  2.21 },
            { -0.81,  0.41,  3.18 },
            {  3.09,  0.00,  3.03 },
            {  2.48, -0.29,  3.71 },
            {  3.91,  0.16,  3.51 },
            {  6.19, -0.06,  3.11 },
            {  6.05,  0.47,  2.33 },
            {  6.59,  0.52,  3.74 },
            { -6.20,  3.05,  3.05 },
            { -6.87,  3.73,  3.17 },
            { -5.55,  3.24,  3.73 },
            { -3.11,  3.06,  3.15 },
            { -3.64,  3.74,  2.71 },
            { -2.32,  3.00,  2.62 },
            {  0.02,  3.05,  3.06 },
            { -0.87,  3.14,  3.38 },
            {  0.48,  3.82,  3.42 },
            {  3.07,  3.10,  3.16 },
            {  3.95,  3.44,  2.97 },
            {  2.76,  2.73,  2.32 },
            {  6.19,  3.07,  3.16 },
            {  7.02,  3.30,  2.72 },
            {  5.52,  3.27,  2.51 },
            { -6.19,  6.24,  3.15 },
            { -5.56,  5.88,  2.52 },
            { -7.05,  5.96,  2.83 },
            { -3.10,  6.14,  3.08 },
            { -2.34,  6.69,  3.27 },
            { -3.86,  6.69,  3.29 },
            { -0.04,  6.24,  3.13 },
            {  0.63,  6.54,  2.53 },
            {  0.08,  5.29,  3.18 },
            {  3.12,  6.24,  3.14 },
            {  3.57,  5.82,  2.40 },
            {  2.23,  5.90,  3.12 },
            {  6.25,  6.19,  3.06 },
            {  5.55,  5.59,  3.32 },
            {  6.08,  6.99,  3.55 },
            { -6.20, -6.16,  6.15 },
            { -6.29, -5.99,  7.09 },
            { -6.09, -7.11,  6.09 },
            { -3.09, -6.19,  6.27 },
            { -2.56, -5.90,  5.52 },
            { -3.80, -6.69,  5.87 },
            {  0.02, -6.25,  6.24 },
            { -0.70, -5.70,  6.51 },
            {  0.25, -5.93,  5.36 },
            {  3.11, -6.18,  6.14 },
            {  3.76, -6.54,  6.74 },
            {  2.29, -6.20,  6.64 },
            {  6.22, -6.17,  6.15 },
            {  6.61, -6.98,  6.47 },
            {  5.56, -5.94,  6.81 },
            { -6.21, -3.10,  6.14 },
            { -6.76, -2.66,  6.78 },
            { -5.51, -3.50,  6.65 },
            { -3.13, -3.05,  6.18 },
            { -2.19, -3.14,  6.34 },
            { -3.50, -3.89,  6.43 },
            {  0.01, -3.06,  6.15 },
            { -0.06, -2.81,  7.07 },
            { -0.25, -3.98,  6.13 },
            {  3.04, -3.09,  6.17 },
            {  3.84, -3.51,  5.84 },
            {  3.25, -2.85,  7.08 },
            {  6.26, -3.13,  6.19 },
            {  6.01, -2.20,  6.09 },
            {  5.47, -3.55,  6.54 },
            { -6.20,  0.01,  6.27 },
            { -5.79, -0.70,  5.78 },
            { -6.67,  0.51,  5.60 },
            { -3.13,  0.01,  6.14 },
            { -3.53, -0.35,  6.94 },
            { -2.21,  0.17,  6.39 },
            { -0.04, -0.04,  6.20 },
            {  0.26,  0.47,  5.46 },
            {  0.51,  0.22,  6.93 },
            {  3.10, -0.05,  6.23 },
            {  2.33,  0.44,  5.95 },
            {  3.85,  0.45,  5.92 },
            {  6.19, -0.01,  6.26 },
            {  7.05,  0.16,  5.88 },
            {  5.58,  0.02,  5.52 },
            { -6.22,  3.04,  6.17 },
            { -5.45,  3.57,  5.95 },
            { -6.62,  3.50,  6.92 },
            { -3.09,  3.16,  6.21 },
            { -3.71,  2.75,  5.61 },
            { -2.60,  2.43,  6.59 },
            { -0.02,  3.10,  6.26 },
            {  0.89,  3.27,  6.05 },
            { -0.44,  2.94,  5.41 },
            {  3.12,  3.04,  6.23 },
            {  2.31,  3.53,  6.43 },
            {  3.59,  3.60,  5.60 },
            {  6.23,  3.05,  6.24 },
            {  5.92,  3.91,  6.54 },
            {  6.02,  3.03,  5.30 },
            { -6.15,  6.21,  6.24 },
            { -6.27,  6.46,  5.32 },
            { -7.00,  5.85,  6.51 },
            { -3.07,  6.15,  6.22 },
            { -3.98,  6.27,  5.94 },
            { -2.66,  7.01,  6.10 },
            {  0.04,  6.20,  6.25 },
            { -0.38,  5.50,  5.75 },
            { -0.36,  7.00,  5.93 },
            {  3.12,  6.15,  6.24 },
            {  3.66,  6.88,  5.93 },
            {  2.25,  6.33,  5.86 },
            {  6.20,  6.27,  6.19 },
            {  5.46,  5.65,  6.19 },
            {  6.97,  5.73,  6.39 }
        };
        for(int atom = 0; atom < natoms; ++atom)
            positions.push_back(Vec3(coords[atom][0], coords[atom][1], coords[atom][2])*OpenMM::NmPerAngstrom);
    }else{
        throw exception();
    }

    system.setDefaultPeriodicBoxVectors(Vec3(boxEdgeLength, 0, 0),
                                        Vec3(0, boxEdgeLength, 0),
                                        Vec3(0, 0, boxEdgeLength));

    sig.clear();
    eps.clear();
    bonds.clear();
    for(int atom = 0; atom < natoms; ++atom){
        system.addParticle(masses[atom%RESSIZE]);
        double sigma = 2.0*pow(2.0, -1.0/6.0)*halfrmins[atom%RESSIZE]*OpenMM::NmPerAngstrom;
        double epsilon = fabs(epsilons[atom%RESSIZE])*OpenMM::KJPerKcal;
        sig.push_back(0.5*sigma);
        eps.push_back(2.0*sqrt(epsilon));
        if(atom%RESSIZE == 0){
            bonds.push_back(pair<int, int>(atom, atom+1));
            bonds.push_back(pair<int, int>(atom, atom+2));
        }
        double charge = do_electrostatics ? charges[atom] : 0;
        forceField->addParticle(charge, sigma, epsilon);
    }
}


void test_water2_dpme_energies_forces_no_exclusions() {
    const double cutoff = 7.0*OpenMM::NmPerAngstrom;
    const double dalpha = 4.0124063605;
    const int grid = 32;
    NonbondedForce* forceField = new NonbondedForce();

    vector<Vec3> positions;
    vector<double> epsvals;
    vector<double> sigvals;
    vector<pair<int, int> > bonds;
    System system;

    const int NATOMS = 6;
    double boxEdgeLength = 25*OpenMM::NmPerAngstrom;

    make_waterbox(NATOMS, boxEdgeLength, forceField,  positions, epsvals, sigvals, bonds, system, false);
    forceField->setNonbondedMethod(OpenMM::NonbondedForce::LJPME);
    forceField->setPMEParameters(0.0f, grid, grid, grid);
    forceField->setReciprocalSpaceForceGroup(1);
    forceField->setLJPMEParameters(dalpha, grid, grid, grid);
    forceField->setCutoffDistance(cutoff);
    forceField->setReactionFieldDielectric(1.0);
    system.addForce(forceField);

    // Reference calculation
    VerletIntegrator integrator(0.01);
    Platform& platform = Platform::getPlatformByName("Reference");
    Context context(system, integrator, platform);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy, false, 1<<1);
    double refenergy = state.getPotentialEnergy();
    const vector<Vec3>& refforces = state.getForces();

    // Optimized CPU calculation
    CpuCalcDispersionPmeReciprocalForceKernel pme(CalcPmeReciprocalForceKernel::Name(), platform);
    IO io;
    double selfEwaldEnergy = 0;
    double dalpha6 = pow(dalpha, 6.0);
    for (int i = 0; i < NATOMS; i++) {
        io.posq.push_back((float)positions[i][0]);
        io.posq.push_back((float)positions[i][1]);
        io.posq.push_back((float)positions[i][2]);
        double c6 = 8.0f * pow(sigvals[i], 3) * epsvals[i];
        io.posq.push_back(c6);
        selfEwaldEnergy += dalpha6 * c6 * c6 / 12.0;
    }
    pme.initialize(grid, grid, grid, NATOMS, dalpha, false);
    Vec3 boxVectors[3];
    system.getDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
    pme.beginComputation(io, boxVectors, true);
    double recenergy = pme.finishComputation(io);

    ASSERT_EQUAL_TOL(recenergy, -2.179629087, 5e-3);
    ASSERT_EQUAL_TOL(selfEwaldEnergy, 1.731404285, 1e-5);

    std::vector<Vec3> knownforces(6);
    knownforces[0] = Vec3(    -1.890360546,    -1.890723915,    -1.879662698);
    knownforces[1] = Vec3( -0.003161352455, -0.000922244929, -0.005391616425);
    knownforces[2] = Vec3( 0.0009199035545, -0.001453894176, -0.006188087146);
    knownforces[3] = Vec3(     1.887108856,     1.887241446,      1.89644647);
    knownforces[4] = Vec3( 0.0008242336483,  0.003778910089, -0.002116131106);
    knownforces[5] = Vec3(  0.004912763044,  0.002324059399, -0.002844482646);
    for (int i = 0; i < NATOMS; i++)
        ASSERT_EQUAL_VEC(refforces[i], knownforces[i], 5e-3);

    recenergy += selfEwaldEnergy;

    // See if they match.

    ASSERT_EQUAL_TOL(refenergy, recenergy, 1e-3);
    for (int i = 0; i < NATOMS; i++)
        ASSERT_EQUAL_VEC(refforces[i], Vec3(io.force[4*i], io.force[4*i+1], io.force[4*i+2]), 5e-3);

}


void testPME(bool triclinic) {
    // Create a cloud of random point charges.

    const int numParticles = 51;
    const double boxWidth = 5.0;
    const double cutoff = 1.0;
    Vec3 boxVectors[3];
    if (triclinic) {
        boxVectors[0] = Vec3(boxWidth, 0, 0);
        boxVectors[1] = Vec3(0.2*boxWidth, boxWidth, 0);
        boxVectors[2] = Vec3(-0.3*boxWidth, -0.1*boxWidth, boxWidth);
    }
    else {
        boxVectors[0] = Vec3(boxWidth, 0, 0);
        boxVectors[1] = Vec3(0, boxWidth, 0);
        boxVectors[2] = Vec3(0, 0, boxWidth);
    }
    System system;
    system.setDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
    NonbondedForce* force = new NonbondedForce();
    system.addForce(force);
    vector<Vec3> positions(numParticles);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);

    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0);
        force->addParticle(-1.0+i*2.0/(numParticles-1), 1.0, 0.0);
        positions[i] = Vec3(boxWidth*genrand_real2(sfmt), boxWidth*genrand_real2(sfmt), boxWidth*genrand_real2(sfmt));
    }
    force->setNonbondedMethod(NonbondedForce::PME);
    force->setCutoffDistance(cutoff);
    force->setReciprocalSpaceForceGroup(1);
    force->setEwaldErrorTolerance(1e-4);
    
    // Compute the reciprocal space forces with the reference platform.
    
    Platform& platform = Platform::getPlatformByName("Reference");
    VerletIntegrator integrator(0.01);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    State refState = context.getState(State::Forces | State::Energy, false, 1<<1);
    
    // Now compute them with the optimized kernel.
    
    double alpha;
    int gridx, gridy, gridz;
    NonbondedForceImpl::calcPMEParameters(system, *force, alpha, gridx, gridy, gridz, false);
    CpuCalcPmeReciprocalForceKernel pme(CalcPmeReciprocalForceKernel::Name(), platform);
    IO io;
    double sumSquaredCharges = 0;
    for (int i = 0; i < numParticles; i++) {
        io.posq.push_back(positions[i][0]);
        io.posq.push_back(positions[i][1]);
        io.posq.push_back(positions[i][2]);
        double charge, sigma, epsilon;
        force->getParticleParameters(i, charge, sigma, epsilon);
        io.posq.push_back(charge);
        sumSquaredCharges += charge*charge;
    }
    double ewaldSelfEnergy = -ONE_4PI_EPS0*alpha*sumSquaredCharges/sqrt(M_PI);
    pme.initialize(gridx, gridy, gridz, numParticles, alpha, true);
    pme.beginComputation(io, boxVectors, true);
    double energy = pme.finishComputation(io);

    // See if they match.
    
    ASSERT_EQUAL_TOL(refState.getPotentialEnergy(), energy+ewaldSelfEnergy, 1e-3);
    for (int i = 0; i < numParticles; i++)
        ASSERT_EQUAL_VEC(refState.getForces()[i], Vec3(io.force[4*i], io.force[4*i+1], io.force[4*i+2]), 1e-3);
}

void testLJPME(bool triclinic) {
    // Create a cloud of random LJ particles.

    const int numParticles = 51;
    const double boxWidth = 5.0;
    const double cutoff = 1.0;
    const double alpha = 2.91842;
    Vec3 boxVectors[3];
    if (triclinic) {
        boxVectors[0] = Vec3(boxWidth, 0, 0);
        boxVectors[1] = Vec3(0.2*boxWidth, boxWidth, 0);
        boxVectors[2] = Vec3(-0.3*boxWidth, -0.1*boxWidth, boxWidth);
    }
    else {
        boxVectors[0] = Vec3(boxWidth, 0, 0);
        boxVectors[1] = Vec3(0, boxWidth, 0);
        boxVectors[2] = Vec3(0, 0, boxWidth);
    }
    System system;
    system.setDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
    NonbondedForce* force = new NonbondedForce();
    system.addForce(force);
    vector<Vec3> positions(numParticles);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);

    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0);
        force->addParticle(0, 0.5, 1.0);
        positions[i] = Vec3(boxWidth*genrand_real2(sfmt), boxWidth*genrand_real2(sfmt), boxWidth*genrand_real2(sfmt));
    }
    force->setNonbondedMethod(NonbondedForce::LJPME);
    force->setCutoffDistance(cutoff);
    force->setReciprocalSpaceForceGroup(1);
    force->setLJPMEParameters(alpha, 64, 64, 64);
    
    // Compute the reciprocal space forces with the reference platform.
    
    Platform& platform = Platform::getPlatformByName("Reference");
    VerletIntegrator integrator(0.01);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    State refState = context.getState(State::Forces | State::Energy, false, 1<<1);
    
    // Now compute them with the optimized kernel.
    
    CpuCalcDispersionPmeReciprocalForceKernel pme(CalcDispersionPmeReciprocalForceKernel::Name(), platform);
    IO io;
    double ewaldSelfEnergy = 0;
    for (int i = 0; i < numParticles; i++) {
        io.posq.push_back(positions[i][0]);
        io.posq.push_back(positions[i][1]);
        io.posq.push_back(positions[i][2]);
        double charge, sigma, epsilon;
        force->getParticleParameters(i, charge, sigma, epsilon);
        io.posq.push_back(pow(sigma, 3.0) * 2.0*sqrt(epsilon));
        ewaldSelfEnergy += pow(alpha*sigma, 6.0) * epsilon / 3.0;
    }
    pme.initialize(64, 64, 64, numParticles, alpha, true);
    pme.beginComputation(io, boxVectors, true);
    double energy = pme.finishComputation(io);

    // See if they match.
    
    ASSERT_EQUAL_TOL(refState.getPotentialEnergy(), energy+ewaldSelfEnergy, 1e-3);
    for (int i = 0; i < numParticles; i++)
        ASSERT_EQUAL_VEC(refState.getForces()[i], Vec3(io.force[4*i], io.force[4*i+1], io.force[4*i+2]), 1e-3);
}

int main(int argc, char* argv[]) {
    try {
        if (!CpuCalcPmeReciprocalForceKernel::isProcessorSupported()) {
            cout << "CPU is not supported.  Exiting." << endl;
            return 0;
        }
        testPME(false);
        testPME(true);
        testLJPME(false);
        testLJPME(true);
        test_water2_dpme_energies_forces_no_exclusions();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
