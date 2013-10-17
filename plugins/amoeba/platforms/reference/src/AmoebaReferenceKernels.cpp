/* -------------------------------------------------------------------------- *
 *                               OpenMMAmoeba                                 *
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

#include "AmoebaReferenceKernels.h"
#include "AmoebaReferenceBondForce.h"
#include "AmoebaReferenceAngleForce.h"
#include "AmoebaReferenceInPlaneAngleForce.h"
#include "AmoebaReferencePiTorsionForce.h"
#include "AmoebaReferenceStretchBendForce.h"
#include "AmoebaReferenceOutOfPlaneBendForce.h"
#include "AmoebaReferenceTorsionTorsionForce.h"
#include "AmoebaReferenceVdwForce.h"
#include "AmoebaReferenceWcaDispersionForce.h"
#include "AmoebaReferenceGeneralizedKirkwoodForce.h"
#include "openmm/internal/AmoebaTorsionTorsionForceImpl.h"
#include "openmm/internal/AmoebaWcaDispersionForceImpl.h"
#include "ReferencePlatform.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/AmoebaMultipoleForce.h"
#include "openmm/internal/AmoebaMultipoleForceImpl.h"
#include "openmm/internal/AmoebaVdwForceImpl.h"
#include "openmm/internal/AmoebaGeneralizedKirkwoodForceImpl.h"
#include "openmm/NonbondedForce.h"
#include "openmm/internal/NonbondedForceImpl.h"

#include <cmath>
#ifdef _MSC_VER
#include <windows.h>
#endif

using namespace OpenMM;
using namespace std;

static vector<RealVec>& extractPositions(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->positions);
}

static vector<RealVec>& extractVelocities(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->velocities);
}

static vector<RealVec>& extractForces(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->forces);
}

static RealVec& extractBoxSize(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *(RealVec*) data->periodicBoxSize;
}

// ***************************************************************************

ReferenceCalcAmoebaBondForceKernel::ReferenceCalcAmoebaBondForceKernel(std::string name, const Platform& platform, const System& system) : 
                CalcAmoebaBondForceKernel(name, platform), system(system) {
}

ReferenceCalcAmoebaBondForceKernel::~ReferenceCalcAmoebaBondForceKernel() {
}

void ReferenceCalcAmoebaBondForceKernel::initialize(const System& system, const AmoebaBondForce& force) {

    numBonds = force.getNumBonds();
    for( int ii = 0; ii < numBonds; ii++) {

        int particle1Index, particle2Index;
        double lengthValue, kValue;
        force.getBondParameters(ii, particle1Index, particle2Index, lengthValue, kValue );

        particle1.push_back( particle1Index ); 
        particle2.push_back( particle2Index ); 
        length.push_back(    static_cast<RealOpenMM>( lengthValue ) );
        kQuadratic.push_back( static_cast<RealOpenMM>( kValue ) );
    } 
    globalBondCubic   = static_cast<RealOpenMM>(force.getAmoebaGlobalBondCubic());
    globalBondQuartic = static_cast<RealOpenMM>(force.getAmoebaGlobalBondQuartic());
}

double ReferenceCalcAmoebaBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& posData   = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    AmoebaReferenceBondForce amoebaReferenceBondForce;
    RealOpenMM energy      = amoebaReferenceBondForce.calculateForceAndEnergy( numBonds, posData, particle1, particle2, length, kQuadratic,
                                                                                       globalBondCubic, globalBondQuartic,
                                                                                       forceData );
    return static_cast<double>(energy);
}

void ReferenceCalcAmoebaBondForceKernel::copyParametersToContext(ContextImpl& context, const AmoebaBondForce& force) {
    if (numBonds != force.getNumBonds())
        throw OpenMMException("updateParametersInContext: The number of bonds has changed");

    // Record the values.

    for (int i = 0; i < numBonds; ++i) {
        int particle1Index, particle2Index;
        double lengthValue, kValue;
        force.getBondParameters(i, particle1Index, particle2Index, lengthValue, kValue);
        if (particle1Index != particle1[i] || particle2Index != particle2[i])
            throw OpenMMException("updateParametersInContext: The set of particles in a bond has changed");
        length[i] = (RealOpenMM) lengthValue;
        kQuadratic[i] = (RealOpenMM) kValue;
    }
}

// ***************************************************************************

ReferenceCalcAmoebaAngleForceKernel::ReferenceCalcAmoebaAngleForceKernel(std::string name, const Platform& platform, const System& system) :
            CalcAmoebaAngleForceKernel(name, platform), system(system) {
}

ReferenceCalcAmoebaAngleForceKernel::~ReferenceCalcAmoebaAngleForceKernel() {
}

void ReferenceCalcAmoebaAngleForceKernel::initialize(const System& system, const AmoebaAngleForce& force) {

    numAngles = force.getNumAngles();

    for (int ii = 0; ii < numAngles; ii++) {
        int particle1Index, particle2Index, particle3Index;
        double angleValue, k;
        force.getAngleParameters(ii, particle1Index, particle2Index, particle3Index, angleValue, k);
        particle1.push_back( particle1Index ); 
        particle2.push_back( particle2Index ); 
        particle3.push_back( particle3Index ); 
        angle.push_back(  static_cast<RealOpenMM>( angleValue ) );
        kQuadratic.push_back( static_cast<RealOpenMM>( k) );
    }
    globalAngleCubic    = static_cast<RealOpenMM>(force.getAmoebaGlobalAngleCubic());
    globalAngleQuartic  = static_cast<RealOpenMM>(force.getAmoebaGlobalAngleQuartic());
    globalAnglePentic   = static_cast<RealOpenMM>(force.getAmoebaGlobalAnglePentic());
    globalAngleSextic   = static_cast<RealOpenMM>(force.getAmoebaGlobalAngleSextic());
}

double ReferenceCalcAmoebaAngleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& posData   = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    AmoebaReferenceAngleForce amoebaReferenceAngleForce;
    RealOpenMM energy      = amoebaReferenceAngleForce.calculateForceAndEnergy( numAngles, 
                                       posData, particle1, particle2, particle3, angle, kQuadratic, globalAngleCubic, globalAngleQuartic, globalAnglePentic, globalAngleSextic, forceData );
    return static_cast<double>(energy);
}

void ReferenceCalcAmoebaAngleForceKernel::copyParametersToContext(ContextImpl& context, const AmoebaAngleForce& force) {
    if (numAngles != force.getNumAngles())
        throw OpenMMException("updateParametersInContext: The number of angles has changed");

    // Record the values.

    for (int i = 0; i < numAngles; ++i) {
        int particle1Index, particle2Index, particle3Index;
        double angleValue, k;
        force.getAngleParameters(i, particle1Index, particle2Index, particle3Index, angleValue, k);
        if (particle1Index != particle1[i] || particle2Index != particle2[i] || particle3Index != particle3[i])
            throw OpenMMException("updateParametersInContext: The set of particles in an angle has changed");
        angle[i] = (RealOpenMM) angleValue;
        kQuadratic[i] = (RealOpenMM) k;
    }
}

ReferenceCalcAmoebaInPlaneAngleForceKernel::ReferenceCalcAmoebaInPlaneAngleForceKernel(std::string name, const Platform& platform, const System& system) : 
          CalcAmoebaInPlaneAngleForceKernel(name, platform), system(system) {
}

ReferenceCalcAmoebaInPlaneAngleForceKernel::~ReferenceCalcAmoebaInPlaneAngleForceKernel() {
}

void ReferenceCalcAmoebaInPlaneAngleForceKernel::initialize(const System& system, const AmoebaInPlaneAngleForce& force) {

    numAngles = force.getNumAngles();
    for (int ii = 0; ii < numAngles; ii++) {
        int particle1Index, particle2Index, particle3Index, particle4Index;
        double angleValue, k;
        force.getAngleParameters(ii, particle1Index, particle2Index, particle3Index, particle4Index, angleValue, k);
        particle1.push_back( particle1Index ); 
        particle2.push_back( particle2Index ); 
        particle3.push_back( particle3Index ); 
        particle4.push_back( particle4Index ); 
        angle.push_back(       static_cast<RealOpenMM>( angleValue ) );
        kQuadratic.push_back(  static_cast<RealOpenMM>( k ) );
    }
    globalInPlaneAngleCubic    = static_cast<RealOpenMM>(force.getAmoebaGlobalInPlaneAngleCubic());
    globalInPlaneAngleQuartic  = static_cast<RealOpenMM>(force.getAmoebaGlobalInPlaneAngleQuartic());
    globalInPlaneAnglePentic   = static_cast<RealOpenMM>(force.getAmoebaGlobalInPlaneAnglePentic());
    globalInPlaneAngleSextic   = static_cast<RealOpenMM>(force.getAmoebaGlobalInPlaneAngleSextic());
}

double ReferenceCalcAmoebaInPlaneAngleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {

    vector<RealVec>& posData   = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    AmoebaReferenceInPlaneAngleForce amoebaReferenceInPlaneAngleForce;
    RealOpenMM energy      = amoebaReferenceInPlaneAngleForce.calculateForceAndEnergy( numAngles, posData, particle1, particle2, particle3, particle4, 
                                                                                               angle, kQuadratic, globalInPlaneAngleCubic, globalInPlaneAngleQuartic,
                                                                                               globalInPlaneAnglePentic, globalInPlaneAngleSextic, forceData );
    return static_cast<double>(energy);
}

void ReferenceCalcAmoebaInPlaneAngleForceKernel::copyParametersToContext(ContextImpl& context, const AmoebaInPlaneAngleForce& force) {
    if (numAngles != force.getNumAngles())
        throw OpenMMException("updateParametersInContext: The number of angles has changed");

    // Record the values.

    for (int i = 0; i < numAngles; ++i) {
        int particle1Index, particle2Index, particle3Index, particle4Index;
        double angleValue, k;
        force.getAngleParameters(i, particle1Index, particle2Index, particle3Index, particle4Index, angleValue, k);
        if (particle1Index != particle1[i] || particle2Index != particle2[i] || particle3Index != particle3[i] || particle4Index != particle4[i])
            throw OpenMMException("updateParametersInContext: The set of particles in an angle has changed");
        angle[i] = (RealOpenMM) angleValue;
        kQuadratic[i] = (RealOpenMM) k;
    }
}

ReferenceCalcAmoebaPiTorsionForceKernel::ReferenceCalcAmoebaPiTorsionForceKernel(std::string name, const Platform& platform, const System& system) :
         CalcAmoebaPiTorsionForceKernel(name, platform), system(system) {
}

ReferenceCalcAmoebaPiTorsionForceKernel::~ReferenceCalcAmoebaPiTorsionForceKernel() {
}

void ReferenceCalcAmoebaPiTorsionForceKernel::initialize(const System& system, const AmoebaPiTorsionForce& force) {

    numPiTorsions                     = force.getNumPiTorsions();
    for (int ii = 0; ii < numPiTorsions; ii++) {

        int particle1Index, particle2Index, particle3Index, particle4Index, particle5Index, particle6Index;
        double kTorsionParameter;
        force.getPiTorsionParameters(ii, particle1Index, particle2Index, particle3Index, particle4Index, particle5Index, particle6Index, kTorsionParameter );
        particle1.push_back( particle1Index ); 
        particle2.push_back( particle2Index ); 
        particle3.push_back( particle3Index ); 
        particle4.push_back( particle4Index ); 
        particle5.push_back( particle5Index ); 
        particle6.push_back( particle6Index ); 
        kTorsion.push_back( static_cast<RealOpenMM>(kTorsionParameter) );
    }
}

double ReferenceCalcAmoebaPiTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& posData   = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    AmoebaReferencePiTorsionForce amoebaReferencePiTorsionForce;
    RealOpenMM energy      = amoebaReferencePiTorsionForce.calculateForceAndEnergy( numPiTorsions, posData, particle1, particle2,
                                                                                    particle3, particle4, particle5, particle6,
                                                                                    kTorsion, forceData );
    return static_cast<double>(energy);
}

void ReferenceCalcAmoebaPiTorsionForceKernel::copyParametersToContext(ContextImpl& context, const AmoebaPiTorsionForce& force) {
    if (numPiTorsions != force.getNumPiTorsions())
        throw OpenMMException("updateParametersInContext: The number of torsions has changed");

    // Record the values.

    for (int i = 0; i < numPiTorsions; ++i) {
        int particle1Index, particle2Index, particle3Index, particle4Index, particle5Index, particle6Index;
        double kTorsionParameter;
        force.getPiTorsionParameters(i, particle1Index, particle2Index, particle3Index, particle4Index, particle5Index, particle6Index, kTorsionParameter);
        if (particle1Index != particle1[i] || particle2Index != particle2[i] || particle3Index != particle3[i] ||
            particle4Index != particle4[i] || particle5Index != particle5[i] || particle6Index != particle6[i])
            throw OpenMMException("updateParametersInContext: The set of particles in a torsion has changed");
        kTorsion[i] = (RealOpenMM) kTorsionParameter;
    }
}

ReferenceCalcAmoebaStretchBendForceKernel::ReferenceCalcAmoebaStretchBendForceKernel(std::string name, const Platform& platform, const System& system) :
                   CalcAmoebaStretchBendForceKernel(name, platform), system(system) {
}

ReferenceCalcAmoebaStretchBendForceKernel::~ReferenceCalcAmoebaStretchBendForceKernel() {
}

void ReferenceCalcAmoebaStretchBendForceKernel::initialize(const System& system, const AmoebaStretchBendForce& force) {

    numStretchBends = force.getNumStretchBends();
    for ( int ii = 0; ii < numStretchBends; ii++) {
        int particle1Index, particle2Index, particle3Index;
        double lengthAB, lengthCB, angle, k;
        force.getStretchBendParameters(ii, particle1Index, particle2Index, particle3Index, lengthAB, lengthCB, angle, k);
        particle1.push_back( particle1Index ); 
        particle2.push_back( particle2Index ); 
        particle3.push_back( particle3Index ); 
        lengthABParameters.push_back( static_cast<RealOpenMM>(lengthAB) );
        lengthCBParameters.push_back( static_cast<RealOpenMM>(lengthCB) );
        angleParameters.push_back(    static_cast<RealOpenMM>(angle) );
        kParameters.push_back(        static_cast<RealOpenMM>(k) );
    }
}

double ReferenceCalcAmoebaStretchBendForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& posData   = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    AmoebaReferenceStretchBendForce amoebaReferenceStretchBendForce;
    RealOpenMM energy      = amoebaReferenceStretchBendForce.calculateForceAndEnergy( numStretchBends, posData, particle1, particle2, particle3,
                                                                                      lengthABParameters, lengthCBParameters, angleParameters, kParameters, forceData );
    return static_cast<double>(energy);
}

void ReferenceCalcAmoebaStretchBendForceKernel::copyParametersToContext(ContextImpl& context, const AmoebaStretchBendForce& force) {
    if (numStretchBends != force.getNumStretchBends())
        throw OpenMMException("updateParametersInContext: The number of stretch-bends has changed");

    // Record the values.

    for (int i = 0; i < numStretchBends; ++i) {
        int particle1Index, particle2Index, particle3Index;
        double lengthAB, lengthCB, angle, k;
        force.getStretchBendParameters(i, particle1Index, particle2Index, particle3Index, lengthAB, lengthCB, angle, k);
        if (particle1Index != particle1[i] || particle2Index != particle2[i] || particle3Index != particle3[i])
            throw OpenMMException("updateParametersInContext: The set of particles in a stretch-bend has changed");
        lengthABParameters[i] = (RealOpenMM) lengthAB;
        lengthCBParameters[i] = (RealOpenMM) lengthCB;
        angleParameters[i] = (RealOpenMM) angle;
        kParameters[i] = (RealOpenMM) k;
    }
}

ReferenceCalcAmoebaOutOfPlaneBendForceKernel::ReferenceCalcAmoebaOutOfPlaneBendForceKernel(std::string name, const Platform& platform, const System& system) :
          CalcAmoebaOutOfPlaneBendForceKernel(name, platform), system(system) {
}

ReferenceCalcAmoebaOutOfPlaneBendForceKernel::~ReferenceCalcAmoebaOutOfPlaneBendForceKernel() {
}

void ReferenceCalcAmoebaOutOfPlaneBendForceKernel::initialize(const System& system, const AmoebaOutOfPlaneBendForce& force) {

    numOutOfPlaneBends = force.getNumOutOfPlaneBends();
    for (int ii = 0; ii < numOutOfPlaneBends; ii++) {

        int particle1Index, particle2Index, particle3Index, particle4Index;
        double k;

        force.getOutOfPlaneBendParameters(ii, particle1Index, particle2Index, particle3Index, particle4Index, k);
        particle1.push_back( particle1Index ); 
        particle2.push_back( particle2Index ); 
        particle3.push_back( particle3Index ); 
        particle4.push_back( particle4Index ); 
        kParameters.push_back( static_cast<RealOpenMM>(k) );
    }
    globalOutOfPlaneBendAngleCubic      = static_cast<RealOpenMM>( force.getAmoebaGlobalOutOfPlaneBendCubic());
    globalOutOfPlaneBendAngleQuartic    = static_cast<RealOpenMM>( force.getAmoebaGlobalOutOfPlaneBendQuartic());
    globalOutOfPlaneBendAnglePentic     = static_cast<RealOpenMM>( force.getAmoebaGlobalOutOfPlaneBendPentic());
    globalOutOfPlaneBendAngleSextic     = static_cast<RealOpenMM>( force.getAmoebaGlobalOutOfPlaneBendSextic());

}

double ReferenceCalcAmoebaOutOfPlaneBendForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& posData   = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    AmoebaReferenceOutOfPlaneBendForce amoebaReferenceOutOfPlaneBendForce;
    RealOpenMM energy      = amoebaReferenceOutOfPlaneBendForce.calculateForceAndEnergy( numOutOfPlaneBends, posData,
                                                                                         particle1, particle2, particle3, particle4,
                                                                                         kParameters, 
                                                                                         globalOutOfPlaneBendAngleCubic,
                                                                                         globalOutOfPlaneBendAngleQuartic,
                                                                                         globalOutOfPlaneBendAnglePentic,
                                                                                         globalOutOfPlaneBendAngleSextic, forceData ); 
    return static_cast<double>(energy);
}

void ReferenceCalcAmoebaOutOfPlaneBendForceKernel::copyParametersToContext(ContextImpl& context, const AmoebaOutOfPlaneBendForce& force) {
    if (numOutOfPlaneBends != force.getNumOutOfPlaneBends())
        throw OpenMMException("updateParametersInContext: The number of out-of-plane bends has changed");

    // Record the values.

    for (int i = 0; i < numOutOfPlaneBends; ++i) {
        int particle1Index, particle2Index, particle3Index, particle4Index;
        double k;
        force.getOutOfPlaneBendParameters(i, particle1Index, particle2Index, particle3Index, particle4Index, k);
        if (particle1Index != particle1[i] || particle2Index != particle2[i] || particle3Index != particle3[i] || particle4Index != particle4[i])
            throw OpenMMException("updateParametersInContext: The set of particles in an out-of-plane bend has changed");
        kParameters[i] = (RealOpenMM) k;
    }
}

ReferenceCalcAmoebaTorsionTorsionForceKernel::ReferenceCalcAmoebaTorsionTorsionForceKernel(std::string name, const Platform& platform, const System& system) :
                CalcAmoebaTorsionTorsionForceKernel(name, platform), system(system) {
}

ReferenceCalcAmoebaTorsionTorsionForceKernel::~ReferenceCalcAmoebaTorsionTorsionForceKernel() {
}

void ReferenceCalcAmoebaTorsionTorsionForceKernel::initialize(const System& system, const AmoebaTorsionTorsionForce& force) {

    numTorsionTorsions = force.getNumTorsionTorsions();

    // torsion-torsion parameters

    for (int ii = 0; ii < numTorsionTorsions; ii++) {
        int particle1Index, particle2Index, particle3Index, particle4Index, particle5Index, chiralCheckAtomIndex, gridIndex;
        force.getTorsionTorsionParameters(ii, particle1Index, particle2Index, particle3Index,
                                          particle4Index, particle5Index, chiralCheckAtomIndex, gridIndex);
        particle1.push_back( particle1Index ); 
        particle2.push_back( particle2Index ); 
        particle3.push_back( particle3Index ); 
        particle4.push_back( particle4Index ); 
        particle5.push_back( particle5Index ); 
        chiralCheckAtom.push_back( chiralCheckAtomIndex ); 
        gridIndices.push_back( gridIndex ); 
    }

    // torsion-torsion grids

    numTorsionTorsionGrids = force.getNumTorsionTorsionGrids();
    torsionTorsionGrids.resize(numTorsionTorsionGrids);
    for (int ii = 0; ii < numTorsionTorsionGrids; ii++) {

        const TorsionTorsionGrid grid = force.getTorsionTorsionGrid( ii );
        torsionTorsionGrids[ii].resize( grid.size() );

        // check if grid needs to be reordered: x-angle should be 'slow' index

        TorsionTorsionGrid reorderedGrid;
        int reorder = 0; 
        if( grid[0][0][0] != grid[0][1][0] ){
            AmoebaTorsionTorsionForceImpl::reorderGrid( grid, reorderedGrid );
            reorder = 1; 
        }    

        for (unsigned int kk = 0; kk < grid.size(); kk++) {

            torsionTorsionGrids[ii][kk].resize( grid[kk].size() );
            for (unsigned int jj = 0; jj < grid[kk].size(); jj++) {

                torsionTorsionGrids[ii][kk][jj].resize( grid[kk][jj].size() );
                if( reorder ){
                    for (unsigned int ll = 0; ll < grid[ll][jj].size(); ll++) {
                        torsionTorsionGrids[ii][kk][jj][ll] = static_cast<RealOpenMM>(reorderedGrid[kk][jj][ll]);
                    }
                } else {
                    for (unsigned int ll = 0; ll < grid[ll][jj].size(); ll++) {
                        torsionTorsionGrids[ii][kk][jj][ll] = static_cast<RealOpenMM>(grid[kk][jj][ll]);
                    }
                }
            }
        }
    }
}

double ReferenceCalcAmoebaTorsionTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {

    vector<RealVec>& posData   = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    AmoebaReferenceTorsionTorsionForce amoebaReferenceTorsionTorsionForce;
    RealOpenMM energy      = amoebaReferenceTorsionTorsionForce.calculateForceAndEnergy( numTorsionTorsions, posData,
                                                                                         particle1, particle2, particle3, particle4, particle5,
                                                                                         chiralCheckAtom, gridIndices, torsionTorsionGrids, forceData );
    return static_cast<double>(energy);
}

/* -------------------------------------------------------------------------- *
 *                             AmoebaMultipole                                *
 * -------------------------------------------------------------------------- */

ReferenceCalcAmoebaMultipoleForceKernel::ReferenceCalcAmoebaMultipoleForceKernel(std::string name, const Platform& platform, const System& system) : 
         CalcAmoebaMultipoleForceKernel(name, platform), system(system), numMultipoles(0), mutualInducedMaxIterations(60), mutualInducedTargetEpsilon(1.0e-03),
                                                         usePme(false),alphaEwald(0.0), cutoffDistance(1.0) {  

}

ReferenceCalcAmoebaMultipoleForceKernel::~ReferenceCalcAmoebaMultipoleForceKernel() {
}

void ReferenceCalcAmoebaMultipoleForceKernel::initialize(const System& system, const AmoebaMultipoleForce& force) {

    numMultipoles   = force.getNumMultipoles();

    charges.resize(numMultipoles);
    dipoles.resize(3*numMultipoles);
    quadrupoles.resize(9*numMultipoles);
    tholes.resize(numMultipoles);
    dampingFactors.resize(numMultipoles);
    polarity.resize(numMultipoles);
    axisTypes.resize(numMultipoles);
    multipoleAtomZs.resize(numMultipoles);
    multipoleAtomXs.resize(numMultipoles);
    multipoleAtomYs.resize(numMultipoles);
    multipoleAtomCovalentInfo.resize(numMultipoles);

    int dipoleIndex      = 0;
    int quadrupoleIndex  = 0;
    int maxCovalentRange = 0;
    double totalCharge   = 0.0;
    for( int ii = 0; ii < numMultipoles; ii++ ){

        // multipoles

        int axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY;
        double charge, tholeD, dampingFactorD, polarityD;
        std::vector<double> dipolesD;
        std::vector<double> quadrupolesD;
        force.getMultipoleParameters(ii, charge, dipolesD, quadrupolesD, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY,
                                     tholeD, dampingFactorD, polarityD );

        totalCharge                       += charge;
        axisTypes[ii]                      = axisType;
        multipoleAtomZs[ii]                = multipoleAtomZ;
        multipoleAtomXs[ii]                = multipoleAtomX;
        multipoleAtomYs[ii]                = multipoleAtomY;

        charges[ii]                        = static_cast<RealOpenMM>(charge);
        tholes[ii]                         = static_cast<RealOpenMM>(tholeD);
        dampingFactors[ii]                 = static_cast<RealOpenMM>(dampingFactorD);
        polarity[ii]                       = static_cast<RealOpenMM>(polarityD);

        dipoles[dipoleIndex++]             = static_cast<RealOpenMM>(dipolesD[0]);
        dipoles[dipoleIndex++]             = static_cast<RealOpenMM>(dipolesD[1]);
        dipoles[dipoleIndex++]             = static_cast<RealOpenMM>(dipolesD[2]);
        
        quadrupoles[quadrupoleIndex++]     = static_cast<RealOpenMM>(quadrupolesD[0]);
        quadrupoles[quadrupoleIndex++]     = static_cast<RealOpenMM>(quadrupolesD[1]);
        quadrupoles[quadrupoleIndex++]     = static_cast<RealOpenMM>(quadrupolesD[2]);
        quadrupoles[quadrupoleIndex++]     = static_cast<RealOpenMM>(quadrupolesD[3]);
        quadrupoles[quadrupoleIndex++]     = static_cast<RealOpenMM>(quadrupolesD[4]);
        quadrupoles[quadrupoleIndex++]     = static_cast<RealOpenMM>(quadrupolesD[5]);
        quadrupoles[quadrupoleIndex++]     = static_cast<RealOpenMM>(quadrupolesD[6]);
        quadrupoles[quadrupoleIndex++]     = static_cast<RealOpenMM>(quadrupolesD[7]);
        quadrupoles[quadrupoleIndex++]     = static_cast<RealOpenMM>(quadrupolesD[8]);

        // covalent info

        std::vector< std::vector<int> > covalentLists;
        force.getCovalentMaps(ii, covalentLists );
        multipoleAtomCovalentInfo[ii] = covalentLists;

    }

    polarizationType = force.getPolarizationType();
    if( polarizationType == AmoebaMultipoleForce::Mutual ){
        mutualInducedMaxIterations = force.getMutualInducedMaxIterations();
        mutualInducedTargetEpsilon = force.getMutualInducedTargetEpsilon();
    }

    // PME

    nonbondedMethod  = force.getNonbondedMethod();
    if( nonbondedMethod == AmoebaMultipoleForce::PME ){
        usePme     = true;
        alphaEwald = force.getAEwald();
        cutoffDistance = force.getCutoffDistance();
        force.getPmeGridDimensions(pmeGridDimension);
        if (pmeGridDimension[0] == 0 || alphaEwald == 0.0) {
            NonbondedForce nb;
            nb.setEwaldErrorTolerance(force.getEwaldErrorTolerance());
            nb.setCutoffDistance(force.getCutoffDistance());
            int gridSizeX, gridSizeY, gridSizeZ;
            NonbondedForceImpl::calcPMEParameters(system, nb, alphaEwald, gridSizeX, gridSizeY, gridSizeZ);
            pmeGridDimension[0] = gridSizeX;
            pmeGridDimension[1] = gridSizeY;
            pmeGridDimension[2] = gridSizeZ;
        }    
    } else {
        usePme = false;
    }
    return;
}

AmoebaReferenceMultipoleForce* ReferenceCalcAmoebaMultipoleForceKernel::setupAmoebaReferenceMultipoleForce(ContextImpl& context )
{

    // amoebaReferenceMultipoleForce is set to AmoebaReferenceGeneralizedKirkwoodForce if AmoebaGeneralizedKirkwoodForce is present
    // amoebaReferenceMultipoleForce is set to AmoebaReferencePmeMultipoleForce if 'usePme' is set
    // amoebaReferenceMultipoleForce is set to AmoebaReferenceMultipoleForce otherwise

    // check if AmoebaGeneralizedKirkwoodForce is present 

    ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel* gkKernel = NULL;
    for (unsigned int ii = 0; ii < context.getForceImpls().size() && gkKernel == NULL; ii++) {
        AmoebaGeneralizedKirkwoodForceImpl* gkImpl = dynamic_cast<AmoebaGeneralizedKirkwoodForceImpl*>(context.getForceImpls()[ii]);
        if (gkImpl != NULL) {
            gkKernel = dynamic_cast<ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel*>(&gkImpl->getKernel().getImpl());
        }
    }    

    AmoebaReferenceMultipoleForce* amoebaReferenceMultipoleForce = NULL;
    if( gkKernel ){

        // amoebaReferenceGeneralizedKirkwoodForce is deleted in AmoebaReferenceGeneralizedKirkwoodMultipoleForce
        // destructor

        AmoebaReferenceGeneralizedKirkwoodForce* amoebaReferenceGeneralizedKirkwoodForce = new AmoebaReferenceGeneralizedKirkwoodForce();
        amoebaReferenceGeneralizedKirkwoodForce->setNumParticles( gkKernel->getNumParticles() );
        amoebaReferenceGeneralizedKirkwoodForce->setSoluteDielectric( gkKernel->getSoluteDielectric() );
        amoebaReferenceGeneralizedKirkwoodForce->setSolventDielectric( gkKernel->getSolventDielectric() );
        amoebaReferenceGeneralizedKirkwoodForce->setDielectricOffset( gkKernel->getDielectricOffset() );
        amoebaReferenceGeneralizedKirkwoodForce->setProbeRadius( gkKernel->getProbeRadius() );
        amoebaReferenceGeneralizedKirkwoodForce->setSurfaceAreaFactor( gkKernel->getSurfaceAreaFactor() );
        amoebaReferenceGeneralizedKirkwoodForce->setIncludeCavityTerm( gkKernel->getIncludeCavityTerm() );
        amoebaReferenceGeneralizedKirkwoodForce->setDirectPolarization( gkKernel->getDirectPolarization() );

        vector<RealOpenMM> parameters; 
        gkKernel->getAtomicRadii( parameters );
        amoebaReferenceGeneralizedKirkwoodForce->setAtomicRadii( parameters );

        gkKernel->getScaleFactors( parameters );
        amoebaReferenceGeneralizedKirkwoodForce->setScaleFactors( parameters );

        gkKernel->getCharges( parameters );
        amoebaReferenceGeneralizedKirkwoodForce->setCharges( parameters );

        // calculate Grycuk Born radii

        vector<RealVec>& posData   = extractPositions(context);
        amoebaReferenceGeneralizedKirkwoodForce->calculateGrycukBornRadii( posData );

        amoebaReferenceMultipoleForce = new AmoebaReferenceGeneralizedKirkwoodMultipoleForce( amoebaReferenceGeneralizedKirkwoodForce );

    } else if( usePme ) {

         AmoebaReferencePmeMultipoleForce* amoebaReferencePmeMultipoleForce = new AmoebaReferencePmeMultipoleForce( );
         amoebaReferencePmeMultipoleForce->setAlphaEwald( alphaEwald );
         amoebaReferencePmeMultipoleForce->setCutoffDistance( cutoffDistance );
         amoebaReferencePmeMultipoleForce->setPmeGridDimensions( pmeGridDimension );
         RealVec& box = extractBoxSize(context);
         double minAllowedSize = 1.999999*cutoffDistance;
         if (box[0] < minAllowedSize || box[1] < minAllowedSize || box[2] < minAllowedSize){
            throw OpenMMException("The periodic box size has decreased to less than twice the nonbonded cutoff.");
         }
         amoebaReferencePmeMultipoleForce->setPeriodicBoxSize(box);
         amoebaReferenceMultipoleForce = static_cast<AmoebaReferenceMultipoleForce*>(amoebaReferencePmeMultipoleForce);

    } else {
         amoebaReferenceMultipoleForce = new AmoebaReferenceMultipoleForce( AmoebaReferenceMultipoleForce::NoCutoff );
    }

    // set polarization type

    if( polarizationType == AmoebaMultipoleForce::Mutual ){
        amoebaReferenceMultipoleForce->setPolarizationType( AmoebaReferenceMultipoleForce::Mutual );
        amoebaReferenceMultipoleForce->setMutualInducedDipoleTargetEpsilon( mutualInducedTargetEpsilon );
        amoebaReferenceMultipoleForce->setMaximumMutualInducedDipoleIterations( mutualInducedMaxIterations );
    } else if( polarizationType == AmoebaMultipoleForce::Direct ){
        amoebaReferenceMultipoleForce->setPolarizationType( AmoebaReferenceMultipoleForce::Direct );
    } else {
        throw OpenMMException("Polarization type not recognzied." );
    }

    return amoebaReferenceMultipoleForce;

}

double ReferenceCalcAmoebaMultipoleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {

    AmoebaReferenceMultipoleForce* amoebaReferenceMultipoleForce = setupAmoebaReferenceMultipoleForce( context );

    vector<RealVec>& posData   = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    RealOpenMM energy          = amoebaReferenceMultipoleForce->calculateForceAndEnergy( posData, charges, dipoles, quadrupoles, tholes,
                                                                                         dampingFactors, polarity, axisTypes, 
                                                                                         multipoleAtomZs, multipoleAtomXs, multipoleAtomYs,
                                                                                         multipoleAtomCovalentInfo, forceData);

    delete amoebaReferenceMultipoleForce;

    return static_cast<double>(energy);
}

void ReferenceCalcAmoebaMultipoleForceKernel::getInducedDipoles(ContextImpl& context, vector<Vec3>& outputDipoles) {
    int numParticles = context.getSystem().getNumParticles();
    outputDipoles.resize(numParticles);

    // Create an AmoebaReferenceMultipoleForce to do the calculation.
    
    AmoebaReferenceMultipoleForce* amoebaReferenceMultipoleForce = setupAmoebaReferenceMultipoleForce( context );
    vector<RealVec>& posData = extractPositions(context);
    
    // Retrieve the induced dipoles.
    
    vector<RealVec> inducedDipoles;
    amoebaReferenceMultipoleForce->calculateInducedDipoles(posData, charges, dipoles, quadrupoles, tholes,
            dampingFactors, polarity, axisTypes, multipoleAtomZs, multipoleAtomXs, multipoleAtomYs, multipoleAtomCovalentInfo, inducedDipoles);
    for (int i = 0; i < numParticles; i++)
        outputDipoles[i] = inducedDipoles[i];
    delete amoebaReferenceMultipoleForce;
}

void ReferenceCalcAmoebaMultipoleForceKernel::getElectrostaticPotential(ContextImpl& context, const std::vector< Vec3 >& inputGrid,
                                                                        std::vector< double >& outputElectrostaticPotential ){

    AmoebaReferenceMultipoleForce* amoebaReferenceMultipoleForce = setupAmoebaReferenceMultipoleForce( context );
    vector<RealVec>& posData                                     = extractPositions(context);
    vector<RealVec> grid( inputGrid.size() );
    vector<RealOpenMM> potential( inputGrid.size() );
    for( unsigned int ii = 0; ii < inputGrid.size(); ii++ ){
        grid[ii] = inputGrid[ii];
    }
    amoebaReferenceMultipoleForce->calculateElectrostaticPotential( posData, charges, dipoles, quadrupoles, tholes,
                                                                    dampingFactors, polarity, axisTypes, 
                                                                    multipoleAtomZs, multipoleAtomXs, multipoleAtomYs,
                                                                    multipoleAtomCovalentInfo, grid, potential );

    outputElectrostaticPotential.resize( inputGrid.size() );
    for( unsigned int ii = 0; ii < inputGrid.size(); ii++ ){
        outputElectrostaticPotential[ii] = potential[ii];
    }

    delete amoebaReferenceMultipoleForce;
}

void ReferenceCalcAmoebaMultipoleForceKernel::getSystemMultipoleMoments(ContextImpl& context, std::vector< double >& outputMultipoleMoments){

    // retrieve masses

    const System& system             = context.getSystem();
    vector<RealOpenMM> masses;
    for (int i = 0; i <  system.getNumParticles(); ++i) {
        masses.push_back( static_cast<RealOpenMM>(system.getParticleMass(i)) );
    }    

    AmoebaReferenceMultipoleForce* amoebaReferenceMultipoleForce = setupAmoebaReferenceMultipoleForce( context );
    vector<RealVec>& posData                                     = extractPositions(context);
    amoebaReferenceMultipoleForce->calculateAmoebaSystemMultipoleMoments( masses, posData, charges, dipoles, quadrupoles, tholes,
                                                                          dampingFactors, polarity, axisTypes, 
                                                                          multipoleAtomZs, multipoleAtomXs, multipoleAtomYs,
                                                                          multipoleAtomCovalentInfo, outputMultipoleMoments );

    delete amoebaReferenceMultipoleForce;
}

void ReferenceCalcAmoebaMultipoleForceKernel::copyParametersToContext(ContextImpl& context, const AmoebaMultipoleForce& force) {
    if (numMultipoles != force.getNumMultipoles())
        throw OpenMMException("updateParametersInContext: The number of multipoles has changed");

    // Record the values.

    int dipoleIndex = 0;
    int quadrupoleIndex = 0;
    for (int i = 0; i < numMultipoles; ++i) {
        int axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY;
        double charge, tholeD, dampingFactorD, polarityD;
        std::vector<double> dipolesD;
        std::vector<double> quadrupolesD;
        force.getMultipoleParameters(i, charge, dipolesD, quadrupolesD, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY, tholeD, dampingFactorD, polarityD);
        axisTypes[i] = axisType;
        multipoleAtomZs[i] = multipoleAtomZ;
        multipoleAtomXs[i] = multipoleAtomX;
        multipoleAtomYs[i] = multipoleAtomY;
        charges[i] = (RealOpenMM) charge;
        tholes[i] = (RealOpenMM) tholeD;
        dampingFactors[i] = (RealOpenMM) dampingFactorD;
        polarity[i] = (RealOpenMM) polarityD;
        dipoles[dipoleIndex++] = (RealOpenMM) dipolesD[0];
        dipoles[dipoleIndex++] = (RealOpenMM) dipolesD[1];
        dipoles[dipoleIndex++] = (RealOpenMM) dipolesD[2];
        quadrupoles[quadrupoleIndex++] = (RealOpenMM) quadrupolesD[0];
        quadrupoles[quadrupoleIndex++] = (RealOpenMM) quadrupolesD[1];
        quadrupoles[quadrupoleIndex++] = (RealOpenMM) quadrupolesD[2];
        quadrupoles[quadrupoleIndex++] = (RealOpenMM) quadrupolesD[3];
        quadrupoles[quadrupoleIndex++] = (RealOpenMM) quadrupolesD[4];
        quadrupoles[quadrupoleIndex++] = (RealOpenMM) quadrupolesD[5];
        quadrupoles[quadrupoleIndex++] = (RealOpenMM) quadrupolesD[6];
        quadrupoles[quadrupoleIndex++] = (RealOpenMM) quadrupolesD[7];
        quadrupoles[quadrupoleIndex++] = (RealOpenMM) quadrupolesD[8];
    }
}

/* -------------------------------------------------------------------------- *
 *                       AmoebaGeneralizedKirkwood                            *
 * -------------------------------------------------------------------------- */

ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel::ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel(std::string name, const Platform& platform, const System& system) : 
           CalcAmoebaGeneralizedKirkwoodForceKernel(name, platform), system(system) {
}

ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel::~ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel() {
}

int ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel::getNumParticles( void ) const {
    return numParticles;
}

int ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel::getIncludeCavityTerm( void ) const {
    return includeCavityTerm;
}

int ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel::getDirectPolarization( void ) const {
    return directPolarization;
}

RealOpenMM ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel::getSoluteDielectric( void ) const {
    return soluteDielectric;
}

RealOpenMM ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel::getSolventDielectric( void ) const {
    return solventDielectric;
}

RealOpenMM ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel::getDielectricOffset( void ) const {
    return dielectricOffset;
}

RealOpenMM ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel::getProbeRadius( void ) const {
    return probeRadius;
}

RealOpenMM ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel::getSurfaceAreaFactor( void ) const {
    return surfaceAreaFactor;
}

void ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel::getAtomicRadii( vector<RealOpenMM>& outputAtomicRadii ) const {
    outputAtomicRadii.resize( atomicRadii.size() );
    copy( atomicRadii.begin(), atomicRadii.end(), outputAtomicRadii.begin() );
}

void ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel::getScaleFactors( vector<RealOpenMM>& outputScaleFactors ) const {
    outputScaleFactors.resize( scaleFactors.size() );
    copy( scaleFactors.begin(), scaleFactors.end(), outputScaleFactors.begin() );
}

void ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel::getCharges( vector<RealOpenMM>& outputCharges ) const {
    outputCharges.resize( charges.size() );
    copy( charges.begin(), charges.end(), outputCharges.begin() );
}

void ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel::initialize(const System& system, const AmoebaGeneralizedKirkwoodForce& force) {

    // check that AmoebaMultipoleForce is present

    const AmoebaMultipoleForce* amoebaMultipoleForce = NULL;
    for (int ii = 0; ii < system.getNumForces() && amoebaMultipoleForce == NULL; ii++) {
        amoebaMultipoleForce = dynamic_cast<const AmoebaMultipoleForce*>(&system.getForce(ii));
    }

    if (amoebaMultipoleForce == NULL) {
        throw OpenMMException("AmoebaGeneralizedKirkwoodForce requires the System to also contain an AmoebaMultipoleForce.");
    }

    if (amoebaMultipoleForce->getNonbondedMethod() != AmoebaMultipoleForce::NoCutoff ) {
        throw OpenMMException("AmoebaGeneralizedKirkwoodForce requires the AmoebaMultipoleForce use the NoCutoff nonbonded method.");
    }

    numParticles = system.getNumParticles();

    for( int ii = 0; ii < numParticles; ii++ ){

        double particleCharge, particleRadius, scalingFactor;
        force.getParticleParameters(ii, particleCharge, particleRadius, scalingFactor);
        atomicRadii.push_back( static_cast<RealOpenMM>( particleRadius ) );
        scaleFactors.push_back( static_cast<RealOpenMM>( scalingFactor ) );
        charges.push_back( static_cast<RealOpenMM>( particleCharge ) );

        // Make sure the charge matches the one specified by the AmoebaMultipoleForce.

        double charge2, thole, damping, polarity;
        int axisType, atomX, atomY, atomZ;
        vector<double> dipole, quadrupole;
        amoebaMultipoleForce->getMultipoleParameters( ii, charge2, dipole, quadrupole, axisType, atomZ, atomX, atomY, thole, damping, polarity);
        if ( particleCharge != charge2 ){
            throw OpenMMException("AmoebaGeneralizedKirkwoodForce and AmoebaMultipoleForce must specify the same charge for every atom.");
        }

    }   
    includeCavityTerm  = force.getIncludeCavityTerm();
    soluteDielectric   = static_cast<RealOpenMM>( force.getSoluteDielectric() );
    solventDielectric  = static_cast<RealOpenMM>( force.getSolventDielectric() );
    dielectricOffset   = static_cast<RealOpenMM>( 0.009 );
    probeRadius        = static_cast<RealOpenMM>( force.getProbeRadius() ), 
    surfaceAreaFactor  = static_cast<RealOpenMM>( force.getSurfaceAreaFactor() ); 
    directPolarization = amoebaMultipoleForce->getPolarizationType() == AmoebaMultipoleForce::Direct ? 1 : 0;
}

double ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    // handled in AmoebaReferenceGeneralizedKirkwoodMultipoleForce, a derived class of the class AmoebaReferenceMultipoleForce
    return 0.0;
}

void ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel::copyParametersToContext(ContextImpl& context, const AmoebaGeneralizedKirkwoodForce& force) {
    if (numParticles != force.getNumParticles())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");

    // Record the values.

    for (int i = 0; i < numParticles; ++i) {
        double particleCharge, particleRadius, scalingFactor;
        force.getParticleParameters(i, particleCharge, particleRadius, scalingFactor);
        atomicRadii[i] = particleRadius;
        scaleFactors[i] = scalingFactor;
        charges[i] = particleCharge;
    }
}

ReferenceCalcAmoebaVdwForceKernel::ReferenceCalcAmoebaVdwForceKernel(std::string name, const Platform& platform, const System& system) :
       CalcAmoebaVdwForceKernel(name, platform), system(system) {
    useCutoff = 0;
    usePBC = 0;
    cutoff = 1.0e+10;
    neighborList = NULL;
}

ReferenceCalcAmoebaVdwForceKernel::~ReferenceCalcAmoebaVdwForceKernel() {
    if( neighborList ){
        delete neighborList;
    } 
}

void ReferenceCalcAmoebaVdwForceKernel::initialize(const System& system, const AmoebaVdwForce& force) {

    // per-particle parameters

    numParticles = system.getNumParticles();

    indexIVs.resize( numParticles );
    allExclusions.resize( numParticles );
    sigmas.resize( numParticles );
    epsilons.resize( numParticles );
    reductions.resize( numParticles );

    for( int ii = 0; ii < numParticles; ii++ ){

        int indexIV;
        double sigma, epsilon, reduction;
        std::vector<int> exclusions;

        force.getParticleParameters( ii, indexIV, sigma, epsilon, reduction );
        force.getParticleExclusions( ii, exclusions );
        for( unsigned int jj = 0; jj < exclusions.size(); jj++ ){
           allExclusions[ii].insert( exclusions[jj] );
        }

        indexIVs[ii]      = indexIV;
        sigmas[ii]        = static_cast<RealOpenMM>( sigma );
        epsilons[ii]      = static_cast<RealOpenMM>( epsilon );
        reductions[ii]    = static_cast<RealOpenMM>( reduction );
    }   
    sigmaCombiningRule     = force.getSigmaCombiningRule();
    epsilonCombiningRule   = force.getEpsilonCombiningRule();
    useCutoff              = (force.getNonbondedMethod() != AmoebaVdwForce::NoCutoff);
    usePBC                 = (force.getNonbondedMethod() == AmoebaVdwForce::CutoffPeriodic);
    cutoff                 = force.getCutoff();
    neighborList           = useCutoff ? new NeighborList() : NULL;
    dispersionCoefficient  = force.getUseDispersionCorrection() ?  AmoebaVdwForceImpl::calcDispersionCorrection(system, force) : 0.0;

}

double ReferenceCalcAmoebaVdwForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {

    vector<RealVec>& posData   = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    AmoebaReferenceVdwForce vdwForce( sigmaCombiningRule, epsilonCombiningRule );
    RealOpenMM energy;
    if( useCutoff ){
        vdwForce.setCutoff( cutoff );
        computeNeighborListVoxelHash( *neighborList, numParticles, posData, allExclusions, extractBoxSize(context), usePBC, cutoff, 0.0);
        if( usePBC ){
            vdwForce.setNonbondedMethod( AmoebaReferenceVdwForce::CutoffPeriodic);
            RealVec& box = extractBoxSize(context);
            double minAllowedSize = 1.999999*cutoff;
            if (box[0] < minAllowedSize || box[1] < minAllowedSize || box[2] < minAllowedSize){
                throw OpenMMException("The periodic box size has decreased to less than twice the cutoff.");
            }
            vdwForce.setPeriodicBox(box);
            energy  = vdwForce.calculateForceAndEnergy( numParticles, posData, indexIVs, sigmas, epsilons, reductions, *neighborList, forceData);
            energy += dispersionCoefficient/(box[0]*box[1]*box[2]);
        } else {
            vdwForce.setNonbondedMethod( AmoebaReferenceVdwForce::CutoffNonPeriodic);
        }
    } else {
        vdwForce.setNonbondedMethod( AmoebaReferenceVdwForce::NoCutoff );
        energy = vdwForce.calculateForceAndEnergy( numParticles, posData, indexIVs, sigmas, epsilons, reductions, allExclusions, forceData);
    }
    return static_cast<double>(energy);
}

void ReferenceCalcAmoebaVdwForceKernel::copyParametersToContext(ContextImpl& context, const AmoebaVdwForce& force) {
    if (numParticles != force.getNumParticles())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");

    // Record the values.

    for (int i = 0; i < numParticles; ++i) {
        int indexIV;
        double sigma, epsilon, reduction;
        force.getParticleParameters(i, indexIV, sigma, epsilon, reduction);
        indexIVs[i] = indexIV;
        sigmas[i] = (RealOpenMM) sigma;
        epsilons[i] = (RealOpenMM) epsilon;
        reductions[i]= (RealOpenMM) reduction;
    }
}

/* -------------------------------------------------------------------------- *
 *                           AmoebaWcaDispersion                              *
 * -------------------------------------------------------------------------- */

ReferenceCalcAmoebaWcaDispersionForceKernel::ReferenceCalcAmoebaWcaDispersionForceKernel(std::string name, const Platform& platform, const System& system) : 
           CalcAmoebaWcaDispersionForceKernel(name, platform), system(system) {
}

ReferenceCalcAmoebaWcaDispersionForceKernel::~ReferenceCalcAmoebaWcaDispersionForceKernel() {
}

void ReferenceCalcAmoebaWcaDispersionForceKernel::initialize(const System& system, const AmoebaWcaDispersionForce& force) {

    // per-particle parameters

    numParticles = system.getNumParticles();
    radii.resize(numParticles);
    epsilons.resize(numParticles);
    for( int ii = 0; ii < numParticles; ii++ ){

        double radius, epsilon;
        force.getParticleParameters( ii, radius, epsilon );

        radii[ii]         = static_cast<RealOpenMM>( radius );
        epsilons[ii]      = static_cast<RealOpenMM>( epsilon );
    }   

    totalMaximumDispersionEnergy = static_cast<RealOpenMM>( AmoebaWcaDispersionForceImpl::getTotalMaximumDispersionEnergy( force ) );

    epso                         = static_cast<RealOpenMM>( force.getEpso()   );
    epsh                         = static_cast<RealOpenMM>( force.getEpsh()   );
    rmino                        = static_cast<RealOpenMM>( force.getRmino()  );
    rminh                        = static_cast<RealOpenMM>( force.getRminh()  );
    awater                       = static_cast<RealOpenMM>( force.getAwater() );
    shctd                        = static_cast<RealOpenMM>( force.getShctd()  );
    dispoff                      = static_cast<RealOpenMM>( force.getDispoff());
    slevy                        = static_cast<RealOpenMM>( force.getSlevy()  );
}

double ReferenceCalcAmoebaWcaDispersionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& posData   = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    AmoebaReferenceWcaDispersionForce amoebaReferenceWcaDispersionForce( epso, epsh, rmino, rminh, awater, shctd, dispoff, slevy );
    RealOpenMM energy      = amoebaReferenceWcaDispersionForce.calculateForceAndEnergy( numParticles, posData, radii, epsilons, totalMaximumDispersionEnergy, forceData);
    return static_cast<double>(energy);
}

void ReferenceCalcAmoebaWcaDispersionForceKernel::copyParametersToContext(ContextImpl& context, const AmoebaWcaDispersionForce& force) {
    if (numParticles != force.getNumParticles())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");

    // Record the values.

    for (int i = 0; i < numParticles; ++i) {
        double radius, epsilon;
        force.getParticleParameters(i, radius, epsilon);
        radii[i] = (RealOpenMM) radius;
        epsilons[i] = (RealOpenMM) epsilon;
    }
    totalMaximumDispersionEnergy = (RealOpenMM) AmoebaWcaDispersionForceImpl::getTotalMaximumDispersionEnergy(force);
}
