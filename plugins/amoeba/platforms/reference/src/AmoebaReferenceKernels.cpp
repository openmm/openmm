/* -------------------------------------------------------------------------- *
 *                               AmoebaOpenMM                                 *
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
#include "AmoebaReferenceHarmonicBondForce.h"
#include "AmoebaReferenceHarmonicAngleForce.h"
#include "AmoebaReferenceHarmonicInPlaneAngleForce.h"
#include "AmoebaReferenceTorsionForce.h"
#include "AmoebaReferencePiTorsionForce.h"
#include "AmoebaReferenceStretchBendForce.h"
#include "AmoebaReferenceOutOfPlaneBendForce.h"
#include "AmoebaReferenceTorsionTorsionForce.h"
#include "AmoebaReferenceMultipoleForce.h"
#include "AmoebaReferenceVdwForce.h"
#include "AmoebaReferenceWcaDispersionForce.h"
#include "internal/AmoebaWcaDispersionForceImpl.h"
#include "ReferencePlatform.h"
#include "openmm/internal/ContextImpl.h"
#include "AmoebaMultipoleForce.h"
#include "internal/AmoebaMultipoleForceImpl.h"

#include <cmath>
#ifdef _MSC_VER
#include <windows.h>
#endif

using namespace OpenMM;
using namespace std;

static RealOpenMM** extractPositions(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return (RealOpenMM**) data->positions;
}
 
static RealOpenMM** extractVelocities(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return (RealOpenMM**) data->velocities;
}
 
static RealOpenMM** extractForces(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return (RealOpenMM**) data->forces;
}
 
static RealOpenMM* extractBoxSize(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return (RealOpenMM*) data->periodicBoxSize;
}

// ***************************************************************************

ReferenceCalcAmoebaHarmonicBondForceKernel::ReferenceCalcAmoebaHarmonicBondForceKernel(std::string name, const Platform& platform, System& system) : 
                CalcAmoebaHarmonicBondForceKernel(name, platform), system(system) {
}

ReferenceCalcAmoebaHarmonicBondForceKernel::~ReferenceCalcAmoebaHarmonicBondForceKernel() {
}

void ReferenceCalcAmoebaHarmonicBondForceKernel::initialize(const System& system, const AmoebaHarmonicBondForce& force) {

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
    globalHarmonicBondCubic   = static_cast<RealOpenMM>(force.getAmoebaGlobalHarmonicBondCubic());
    globalHarmonicBondQuartic = static_cast<RealOpenMM>(force.getAmoebaGlobalHarmonicBondQuartic());
}

double ReferenceCalcAmoebaHarmonicBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    RealOpenMM** posData   = extractPositions(context);
    RealOpenMM** forceData = extractForces(context);
    AmoebaReferenceHarmonicBondForce amoebaReferenceHarmonicBondForce;
    RealOpenMM energy      = amoebaReferenceHarmonicBondForce.calculateForceAndEnergy( numBonds, posData, particle1, particle2, length, kQuadratic,
                                                                                       globalHarmonicBondCubic, globalHarmonicBondQuartic,
                                                                                       forceData );
    return static_cast<double>(energy);
}

ReferenceCalcAmoebaHarmonicAngleForceKernel::ReferenceCalcAmoebaHarmonicAngleForceKernel(std::string name, const Platform& platform, System& system) :
            CalcAmoebaHarmonicAngleForceKernel(name, platform), system(system) {
}

ReferenceCalcAmoebaHarmonicAngleForceKernel::~ReferenceCalcAmoebaHarmonicAngleForceKernel() {
}

void ReferenceCalcAmoebaHarmonicAngleForceKernel::initialize(const System& system, const AmoebaHarmonicAngleForce& force) {

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
    globalHarmonicAngleCubic    = static_cast<RealOpenMM>(force.getAmoebaGlobalHarmonicAngleCubic());
    globalHarmonicAngleQuartic  = static_cast<RealOpenMM>(force.getAmoebaGlobalHarmonicAngleQuartic());
    globalHarmonicAnglePentic   = static_cast<RealOpenMM>(force.getAmoebaGlobalHarmonicAnglePentic());
    globalHarmonicAngleSextic   = static_cast<RealOpenMM>(force.getAmoebaGlobalHarmonicAngleSextic());
}

double ReferenceCalcAmoebaHarmonicAngleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    RealOpenMM** posData   = extractPositions(context);
    RealOpenMM** forceData = extractForces(context);
    AmoebaReferenceHarmonicAngleForce amoebaReferenceHarmonicAngleForce;
    RealOpenMM energy      = amoebaReferenceHarmonicAngleForce.calculateForceAndEnergy( numAngles, 
                                       posData, particle1, particle2, particle3, angle, kQuadratic, globalHarmonicAngleCubic, globalHarmonicAngleQuartic, globalHarmonicAnglePentic, globalHarmonicAngleSextic, forceData );
    return static_cast<double>(energy);
}

ReferenceCalcAmoebaHarmonicInPlaneAngleForceKernel::ReferenceCalcAmoebaHarmonicInPlaneAngleForceKernel(std::string name, const Platform& platform, System& system) : 
          CalcAmoebaHarmonicInPlaneAngleForceKernel(name, platform), system(system) {
}

ReferenceCalcAmoebaHarmonicInPlaneAngleForceKernel::~ReferenceCalcAmoebaHarmonicInPlaneAngleForceKernel() {
}

void ReferenceCalcAmoebaHarmonicInPlaneAngleForceKernel::initialize(const System& system, const AmoebaHarmonicInPlaneAngleForce& force) {

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
    globalHarmonicInPlaneAngleCubic    = static_cast<RealOpenMM>(force.getAmoebaGlobalHarmonicInPlaneAngleCubic());
    globalHarmonicInPlaneAngleQuartic  = static_cast<RealOpenMM>(force.getAmoebaGlobalHarmonicInPlaneAngleQuartic());
    globalHarmonicInPlaneAnglePentic   = static_cast<RealOpenMM>(force.getAmoebaGlobalHarmonicInPlaneAnglePentic());
    globalHarmonicInPlaneAngleSextic   = static_cast<RealOpenMM>(force.getAmoebaGlobalHarmonicInPlaneAngleSextic());
}

double ReferenceCalcAmoebaHarmonicInPlaneAngleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {

    RealOpenMM** posData   = extractPositions(context);
    RealOpenMM** forceData = extractForces(context);
    AmoebaReferenceHarmonicInPlaneAngleForce amoebaReferenceHarmonicInPlaneAngleForce;
    RealOpenMM energy      = amoebaReferenceHarmonicInPlaneAngleForce.calculateForceAndEnergy( numAngles, posData, particle1, particle2, particle3, particle4, 
                                                                                               angle, kQuadratic, globalHarmonicInPlaneAngleCubic, globalHarmonicInPlaneAngleQuartic,
                                                                                               globalHarmonicInPlaneAnglePentic, globalHarmonicInPlaneAngleSextic, forceData );
    return static_cast<double>(energy);
}

ReferenceCalcAmoebaTorsionForceKernel::ReferenceCalcAmoebaTorsionForceKernel(std::string name, const Platform& platform, System& system) :
             CalcAmoebaTorsionForceKernel(name, platform), system(system) {
}

ReferenceCalcAmoebaTorsionForceKernel::~ReferenceCalcAmoebaTorsionForceKernel() {
}

void ReferenceCalcAmoebaTorsionForceKernel::initialize(const System& system, const AmoebaTorsionForce& force) {

    numTorsions = force.getNumTorsions();
    torsionParameters1.resize( numTorsions );
    torsionParameters2.resize( numTorsions );
    torsionParameters3.resize( numTorsions );
    for (int ii = 0; ii < numTorsions; ii++) {

        int particle1Index, particle2Index, particle3Index, particle4Index;
        std::vector<double> torsionParameter1;
        std::vector<double> torsionParameter2;
        std::vector<double> torsionParameter3;

        std::vector<RealOpenMM> torsionParameters1F(3);
        std::vector<RealOpenMM> torsionParameters2F(3);
        std::vector<RealOpenMM> torsionParameters3F(3);

        force.getTorsionParameters(ii, particle1Index, particle2Index, particle3Index, particle4Index, torsionParameter1, torsionParameter2, torsionParameter3 );
        particle1.push_back( particle1Index ); 
        particle2.push_back( particle2Index ); 
        particle3.push_back( particle3Index ); 
        particle4.push_back( particle4Index ); 
        torsionParameters1[ii].resize( torsionParameter1.size() );
        torsionParameters2[ii].resize( torsionParameter2.size() );
        torsionParameters3[ii].resize( torsionParameter3.size() );
        for ( unsigned int jj = 0; jj < torsionParameter1.size(); jj++) {
            torsionParameters1[ii][jj] = static_cast<RealOpenMM>(torsionParameter1[jj]);
            torsionParameters2[ii][jj] = static_cast<RealOpenMM>(torsionParameter2[jj]);
            torsionParameters3[ii][jj] = static_cast<RealOpenMM>(torsionParameter3[jj]);
        }
    }
}

double ReferenceCalcAmoebaTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    RealOpenMM** posData   = extractPositions(context);
    RealOpenMM** forceData = extractForces(context);
    AmoebaReferenceTorsionForce amoebaReferenceTorsionForce;
    RealOpenMM energy      = amoebaReferenceTorsionForce.calculateForceAndEnergy( numTorsions, posData, particle1, particle2, particle3, particle4,
                                                                                  torsionParameters1, torsionParameters2, torsionParameters3, forceData );
    return static_cast<double>(energy);
}

ReferenceCalcAmoebaPiTorsionForceKernel::ReferenceCalcAmoebaPiTorsionForceKernel(std::string name, const Platform& platform, System& system) :
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
    RealOpenMM** posData   = extractPositions(context);
    RealOpenMM** forceData = extractForces(context);
    AmoebaReferencePiTorsionForce amoebaReferencePiTorsionForce;
    RealOpenMM energy      = amoebaReferencePiTorsionForce.calculateForceAndEnergy( numPiTorsions, posData, particle1, particle2,
                                                                                    particle3, particle4, particle5, particle6,
                                                                                    kTorsion, forceData );
    return static_cast<double>(energy);
}

ReferenceCalcAmoebaStretchBendForceKernel::ReferenceCalcAmoebaStretchBendForceKernel(std::string name, const Platform& platform, System& system) :
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
    RealOpenMM** posData   = extractPositions(context);
    RealOpenMM** forceData = extractForces(context);
    AmoebaReferenceStretchBendForce amoebaReferenceStretchBendForce;
    RealOpenMM energy      = amoebaReferenceStretchBendForce.calculateForceAndEnergy( numStretchBends, posData, particle1, particle2, particle3,
                                                                                      lengthABParameters, lengthCBParameters, angleParameters, kParameters, forceData );
    return static_cast<double>(energy);
}

ReferenceCalcAmoebaOutOfPlaneBendForceKernel::ReferenceCalcAmoebaOutOfPlaneBendForceKernel(std::string name, const Platform& platform, System& system) :
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
    RealOpenMM** posData   = extractPositions(context);
    RealOpenMM** forceData = extractForces(context);
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

ReferenceCalcAmoebaTorsionTorsionForceKernel::ReferenceCalcAmoebaTorsionTorsionForceKernel(std::string name, const Platform& platform, System& system) :
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
        for (unsigned int kk = 0; kk < grid.size(); kk++) {

            torsionTorsionGrids[ii][kk].resize( grid[kk].size() );
            for (unsigned int jj = 0; jj < grid[kk].size(); jj++) {

                torsionTorsionGrids[ii][kk][jj].resize( grid[kk][jj].size() );
                for (unsigned int ll = 0; ll < grid[ll][jj].size(); ll++) {
                    torsionTorsionGrids[ii][kk][jj][ll] = static_cast<RealOpenMM>(grid[kk][jj][ll]);
                }
            }
        }
    }
}

double ReferenceCalcAmoebaTorsionTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {

    RealOpenMM** posData   = extractPositions(context);
    RealOpenMM** forceData = extractForces(context);
    AmoebaReferenceTorsionTorsionForce amoebaReferenceTorsionTorsionForce;
    RealOpenMM energy      = amoebaReferenceTorsionTorsionForce.calculateForceAndEnergy( numTorsionTorsions, posData,
                                                                                         particle1, particle2, particle3, particle4, particle5,
                                                                                         chiralCheckAtom, gridIndices, torsionTorsionGrids, forceData );
    return static_cast<double>(energy);
}

/* -------------------------------------------------------------------------- *
 *                             AmoebaMultipole                                *
 * -------------------------------------------------------------------------- */

ReferenceCalcAmoebaMultipoleForceKernel::ReferenceCalcAmoebaMultipoleForceKernel(std::string name, const Platform& platform, System& system) : 
         CalcAmoebaMultipoleForceKernel(name, platform), system(system) {
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
        multipoleAtomZs[ii]              = multipoleAtomZ;
        multipoleAtomXs[ii]              = multipoleAtomX;
        multipoleAtomYs[ii]              = multipoleAtomY;

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

    mutualInducedMaxIterations = force.getMutualInducedMaxIterations();
    mutualInducedTargetEpsilon = force.getMutualInducedTargetEpsilon();

    nonbondedMethod = static_cast<int>(force.getNonbondedMethod());
    if( nonbondedMethod != 0 && nonbondedMethod != 1 ){
         throw OpenMMException("AmoebaMultipoleForce nonbonded method not recognized.\n");
    }
}

double ReferenceCalcAmoebaMultipoleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    RealOpenMM** posData   = extractPositions(context);
    RealOpenMM** forceData = extractForces(context);

    AmoebaReferenceMultipoleForce amoebaReferenceMultipoleForce( AmoebaReferenceMultipoleForce::NoCutoff );
    amoebaReferenceMultipoleForce.setMutualInducedDipoleTargetEpsilon( mutualInducedTargetEpsilon );
    amoebaReferenceMultipoleForce.setMaximumMutualInducedDipoleIterations( mutualInducedMaxIterations );

    RealOpenMM energy      = amoebaReferenceMultipoleForce.calculateForceAndEnergy( numMultipoles, posData, 
                                                                                    charges, dipoles, quadrupoles, tholes,
                                                                                    dampingFactors, polarity, axisTypes, 
                                                                                    multipoleAtomZs, multipoleAtomXs, multipoleAtomYs,
                                                                                    multipoleAtomCovalentInfo, forceData);

    return static_cast<double>(energy);
}

///* -------------------------------------------------------------------------- *
// *                       AmoebaGeneralizedKirkwood                            *
// * -------------------------------------------------------------------------- */
//
//ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel::ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel(std::string name, const Platform& platform, System& system) : 
//           CalcAmoebaGeneralizedKirkwoodForceKernel(name, platform), system(system) {
//    data.incrementKernelCount();
//}
//
//ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel::~ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel() {
//    data.decrementKernelCount();
//}
//
//void ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel::initialize(const System& system, const AmoebaGeneralizedKirkwoodForce& force) {
//
//    data.setHasAmoebaGeneralizedKirkwood( true );
//
//    int numParticles = system.getNumParticles();
//
//    std::vector<RealOpenMM> radius(numParticles);
//    std::vector<RealOpenMM> scale(numParticles);
//    std::vector<RealOpenMM> charge(numParticles);
//
//    for( int ii = 0; ii < numParticles; ii++ ){
//        double particleCharge, particleRadius, scalingFactor;
//        force.getParticleParameters(ii, particleCharge, particleRadius, scalingFactor);
//        radius[ii]  = static_cast<RealOpenMM>( particleRadius );
//        scale[ii]   = static_cast<RealOpenMM>( scalingFactor );
//        charge[ii]  = static_cast<RealOpenMM>( particleCharge );
//    }   
//    gpuSetAmoebaObcParameters( data.getAmoebaGpu(), static_cast<RealOpenMM>(force.getSoluteDielectric() ), 
//                               static_cast<RealOpenMM>( force.getSolventDielectric() ), 
//                               static_cast<RealOpenMM>( force.getDielectricOffset() ), radius, scale, charge,
//                               force.getIncludeCavityTerm(),
//                               static_cast<RealOpenMM>( force.getProbeRadius() ), 
//                               static_cast<RealOpenMM>( force.getSurfaceAreaFactor() ) ); 
//}
//
//double ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
//    // handled in computeAmoebaMultipoleForce()
//    return 0.0;
//}

ReferenceCalcAmoebaVdwForceKernel::ReferenceCalcAmoebaVdwForceKernel(std::string name, const Platform& platform, System& system) :
       CalcAmoebaVdwForceKernel(name, platform), system(system) {
}

ReferenceCalcAmoebaVdwForceKernel::~ReferenceCalcAmoebaVdwForceKernel() {
}

void ReferenceCalcAmoebaVdwForceKernel::initialize(const System& system, const AmoebaVdwForce& force) {

    // per-particle parameters

    numParticles = system.getNumParticles();

    indexIVs.resize( numParticles );
    indexClasses.resize( numParticles );
    allExclusions.resize( numParticles );
    sigmas.resize( numParticles );
    epsilons.resize( numParticles );
    reductions.resize( numParticles );

    for( int ii = 0; ii < numParticles; ii++ ){

        int indexIV, indexClass;
        double sigma, epsilon, reduction;
        std::vector<int> exclusions;

        force.getParticleParameters( ii, indexIV, indexClass, sigma, epsilon, reduction );
        force.getParticleExclusions( ii, exclusions );
        for( unsigned int jj = 0; jj < exclusions.size(); jj++ ){
           allExclusions[ii].push_back( exclusions[jj] );
        }

        indexIVs[ii]      = indexIV;
        indexClasses[ii]  = indexClass;
        sigmas[ii]        = static_cast<RealOpenMM>( sigma );
        epsilons[ii]      = static_cast<RealOpenMM>( epsilon );
        reductions[ii]    = static_cast<RealOpenMM>( reduction );
    }   
    sigmaCombiningRule   = force.getSigmaCombiningRule();
    epsilonCombiningRule = force.getEpsilonCombiningRule();
}

double ReferenceCalcAmoebaVdwForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {

    RealOpenMM** posData   = extractPositions(context);
    RealOpenMM** forceData = extractForces(context);
    AmoebaReferenceVdwForce vdwForce( sigmaCombiningRule, epsilonCombiningRule, AmoebaReferenceVdwForce::NoCutoff );
    RealOpenMM energy      = vdwForce.calculateForceAndEnergy( numParticles, posData, indexIVs, indexClasses, sigmas, epsilons, reductions, allExclusions, forceData);
    return static_cast<double>(energy);
}

/* -------------------------------------------------------------------------- *
 *                           AmoebaWcaDispersion                              *
 * -------------------------------------------------------------------------- */

ReferenceCalcAmoebaWcaDispersionForceKernel::ReferenceCalcAmoebaWcaDispersionForceKernel(std::string name, const Platform& platform, System& system) : 
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
    RealOpenMM** posData   = extractPositions(context);
    RealOpenMM** forceData = extractForces(context);
    AmoebaReferenceWcaDispersionForce amoebaReferenceWcaDispersionForce( epso, epsh, rmino, rminh, awater, shctd, dispoff, slevy );
    RealOpenMM energy      = amoebaReferenceWcaDispersionForce.calculateForceAndEnergy( numParticles, posData, radii, epsilons, totalMaximumDispersionEnergy, forceData);
    return static_cast<double>(energy);
}
