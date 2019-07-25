/* -------------------------------------------------------------------------- *
 *                               OpenMMAmoeba                                 *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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

#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/AmoebaTorsionTorsionForceImpl.h"
#include "openmm/amoebaKernels.h"
#include <cstdio>

using namespace OpenMM;

using std::pair;
using std::vector;
using std::set;

AmoebaTorsionTorsionForceImpl::AmoebaTorsionTorsionForceImpl(const AmoebaTorsionTorsionForce& owner) : owner(owner) {
}

AmoebaTorsionTorsionForceImpl::~AmoebaTorsionTorsionForceImpl() {
}

void AmoebaTorsionTorsionForceImpl::initialize(ContextImpl& context) {
    kernel = context.getPlatform().createKernel(CalcAmoebaTorsionTorsionForceKernel::Name(), context);
    kernel.getAs<CalcAmoebaTorsionTorsionForceKernel>().initialize(context.getSystem(), owner);
}

double AmoebaTorsionTorsionForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups&(1<<owner.getForceGroup())) != 0)
        return kernel.getAs<CalcAmoebaTorsionTorsionForceKernel>().execute(context, includeForces, includeEnergy);
    return 0.0;
}

struct IntPair {
    unsigned int index1;
    unsigned int index2;
};

typedef std::map< double, struct IntPair > Map_Double_IntPair;
typedef Map_Double_IntPair::iterator Map_Double_IntPairI;
typedef Map_Double_IntPair::const_iterator Map_Double_IntPairCI;

typedef std::map< double, Map_Double_IntPair > Map_Double_MapDoubleIntPair;
typedef Map_Double_MapDoubleIntPair::iterator Map_Double_MapDoubleIntPairI;
typedef Map_Double_MapDoubleIntPair::const_iterator Map_Double_MapDoubleIntPairCI;

void AmoebaTorsionTorsionForceImpl::reorderGrid(const TorsionTorsionGrid& grid, TorsionTorsionGrid& reorderedGrid) {

    reorderedGrid.resize(grid.size());
    std::vector<Map_Double_IntPair> map_Double_IntPair_Vector(grid.size());
    Map_Double_MapDoubleIntPair mapAngles;

    // (1) set dimensions for reorderd grid
    // (2) build map:
    //         map[angleX][angleY] = <ii, jj> indices
    //         assume map keys are sorted from least to greatest

    for (unsigned int ii = 0; ii < grid.size(); ii++) {
    
        reorderedGrid[ii].resize(grid[ii].size());
        for (unsigned int jj = 0; jj < grid[ii].size(); jj++) {
            reorderedGrid[ii][jj].resize(grid[ii][jj].size());

            double angleX =  grid[ii][jj][0]; 
            double angleY =  grid[ii][jj][1]; 

            if (mapAngles.find(angleX) == mapAngles.end()) {
                if (map_Double_IntPair_Vector[ii].size() > 0) {
                    char buffer[1024];
                    (void) sprintf(buffer, "TorsionTorsion grid reorder: x-angle not set correctly: x=%15.7e y=%15.7e size=%u should be zero; ii/jj indies=%u %u.\n",
                                    angleX, angleY, static_cast<unsigned int>(map_Double_IntPair_Vector[ii].size()), ii, jj);
                    throw OpenMMException(buffer);
                 }
                 mapAngles[angleX] = map_Double_IntPair_Vector[ii];
            }

            Map_Double_IntPair& map_Double_IntPair  = mapAngles[angleX];
            if (map_Double_IntPair.find(angleY) != map_Double_IntPair.end()) {
                char buffer[1024];
                (void) sprintf(buffer, "TorsionTorsion grid reorder: angle pair found twice: %15.7e %15.7e %u\n", angleX, angleY, static_cast<unsigned int>(map_Double_IntPair.size()));
                throw OpenMMException(buffer);
            }
            struct IntPair pair; 
            pair.index1 = ii;
            pair.index2 = jj;
            map_Double_IntPair[angleY] = pair; 
        }
    }

    // load reordered entries

    Map_Double_MapDoubleIntPairCI mapII    = mapAngles.begin();
    Map_Double_IntPair map_Double_IntPair  = mapII->second;
    Map_Double_IntPairCI mapJJ             = map_Double_IntPair.begin();

    for (unsigned int ii = 0; ii < grid.size(); ii++) {
        for (unsigned int jj = 0; jj < grid[ii].size(); jj++) {

            struct IntPair pair  = mapJJ->second;
            int index1           = pair.index1;
            int index2           = pair.index2;

            for (unsigned int kk = 0; kk < grid[ii][jj].size(); kk++) {
                reorderedGrid[ii][jj][kk] = static_cast<float>(grid[index1][index2][kk]);
            }

            // increment map iterators

            mapJJ++;
            if (mapJJ == map_Double_IntPair.end()) {
                mapII++;
                if (mapII == mapAngles.end()) {
                    if ((jj != (grid[ii].size()-1)) && (ii != (grid.size()-1))) {
                        char buffer[1024];
                        (void) sprintf(buffer, "AmoebaTorsionTorsionForceImpl::reorderGrid: error detected with map iterators.\n");
                        throw OpenMMException(buffer);
                    }
                } else {
                    map_Double_IntPair  = mapII->second;
                    mapJJ               = map_Double_IntPair.begin();
                }
            }
        }
    }

    return;
}

std::vector<std::string> AmoebaTorsionTorsionForceImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(CalcAmoebaTorsionTorsionForceKernel::Name());
    return names;
}

