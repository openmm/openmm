/* -------------------------------------------------------------------------- *
 *                                AmoebaOpenMM                                *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2009 Stanford University and the Authors.      *
 * Authors:                                                                   *
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

#include "openmm/Force.h"
#include "openmm/OpenMMException.h"
#include "AmoebaTorsionTorsionForce.h"
#include "internal/AmoebaTorsionTorsionForceImpl.h"

using namespace OpenMM;

AmoebaTorsionTorsionForce::AmoebaTorsionTorsionForce() {
}

int AmoebaTorsionTorsionForce::addTorsionTorsion(int particle1, int particle2, int particle3,
                                                 int particle4, int particle5, int chiralCheckAtomIndex,
                                                 int gridIndex) {
    torsionTorsions.push_back(TorsionTorsionInfo(particle1, particle2, particle3, particle4, particle5, chiralCheckAtomIndex, gridIndex));
    return torsionTorsions.size()-1;
}

void AmoebaTorsionTorsionForce::getTorsionTorsionParameters(int index, int& particle1, int& particle2, int& particle3, int& particle4,
                                                            int& particle5, int& chiralCheckAtomIndex, int& gridIndex) const {
    particle1                = torsionTorsions[index].particle1;
    particle2                = torsionTorsions[index].particle2;
    particle3                = torsionTorsions[index].particle3;
    particle4                = torsionTorsions[index].particle4;
    particle5                = torsionTorsions[index].particle5;
    chiralCheckAtomIndex     = torsionTorsions[index].chiralCheckAtomIndex;
    gridIndex                = torsionTorsions[index].gridIndex;
}

void AmoebaTorsionTorsionForce::setTorsionTorsionParameters(int index, int particle1, int particle2, int particle3, int particle4,
                                                            int particle5, int chiralCheckAtomIndex, int gridIndex) {
    torsionTorsions[index].particle1              = particle1;
    torsionTorsions[index].particle2              = particle2;
    torsionTorsions[index].particle3              = particle3;
    torsionTorsions[index].particle4              = particle4;
    torsionTorsions[index].particle5              = particle5;
    torsionTorsions[index].chiralCheckAtomIndex   = chiralCheckAtomIndex;
    torsionTorsions[index].gridIndex              = gridIndex;
}

const TorsionTorsionGrid& AmoebaTorsionTorsionForce::getTorsionTorsionGrid(int index ) const {
   return torsionTorsionGrids[index].getTorsionTorsionGrid();
}

void AmoebaTorsionTorsionForce::setTorsionTorsionGrid(int index, TorsionTorsionGrid& grid ) {
   if( index >= torsionTorsionGrids.size() ){
      torsionTorsionGrids.resize( index + 1);
   }
   torsionTorsionGrids[index] = grid;
}

ForceImpl* AmoebaTorsionTorsionForce::createImpl() {
    return new AmoebaTorsionTorsionForceImpl(*this);
}
