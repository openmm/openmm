/* -------------------------------------------------------------------------- *
 *                               OpenMMAmoeba                                 *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2020 Stanford University and the Authors.      *
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
#include "AmoebaReferenceTorsionTorsionForce.h"
#include "AmoebaReferenceWcaDispersionForce.h"
#include "AmoebaReferenceGeneralizedKirkwoodForce.h"
#include "openmm/internal/AmoebaTorsionTorsionForceImpl.h"
#include "openmm/internal/AmoebaWcaDispersionForceImpl.h"
#include "ReferencePlatform.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/AmoebaMultipoleForce.h"
#include "openmm/HippoNonbondedForce.h"
#include "openmm/internal/AmoebaMultipoleForceImpl.h"
#include "openmm/internal/AmoebaVdwForceImpl.h"
#include "openmm/internal/AmoebaGeneralizedKirkwoodForceImpl.h"
#include "openmm/NonbondedForce.h"
#include "openmm/internal/NonbondedForceImpl.h"
#include "SimTKReference/AmoebaReferenceHippoNonbondedForce.h"

#include <cmath>
#ifdef _MSC_VER
#include <windows.h>
#endif

using namespace OpenMM;
using namespace std;

static vector<Vec3>& extractPositions(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *data->positions;
}

static vector<Vec3>& extractVelocities(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *data->velocities;
}

static vector<Vec3>& extractForces(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *data->forces;
}

static Vec3& extractBoxSize(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *data->periodicBoxSize;
}

static Vec3* extractBoxVectors(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return data->periodicBoxVectors;
}

// ***************************************************************************

ReferenceCalcAmoebaTorsionTorsionForceKernel::ReferenceCalcAmoebaTorsionTorsionForceKernel(const std::string& name, const Platform& platform, const System& system) :
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
        particle1.push_back(particle1Index); 
        particle2.push_back(particle2Index); 
        particle3.push_back(particle3Index); 
        particle4.push_back(particle4Index); 
        particle5.push_back(particle5Index); 
        chiralCheckAtom.push_back(chiralCheckAtomIndex); 
        gridIndices.push_back(gridIndex); 
    }
    usePeriodic = force.usesPeriodicBoundaryConditions();

    // torsion-torsion grids

    numTorsionTorsionGrids = force.getNumTorsionTorsionGrids();
    torsionTorsionGrids.resize(numTorsionTorsionGrids);
    for (int ii = 0; ii < numTorsionTorsionGrids; ii++) {

        const TorsionTorsionGrid grid = force.getTorsionTorsionGrid(ii);
        torsionTorsionGrids[ii].resize(grid.size());

        // check if grid needs to be reordered: x-angle should be 'slow' index

        TorsionTorsionGrid reorderedGrid;
        int reorder = 0; 
        if (grid[0][0][0] != grid[0][1][0]) {
            AmoebaTorsionTorsionForceImpl::reorderGrid(grid, reorderedGrid);
            reorder = 1; 
        }    

        for (unsigned int kk = 0; kk < grid.size(); kk++) {

            torsionTorsionGrids[ii][kk].resize(grid[kk].size());
            for (unsigned int jj = 0; jj < grid[kk].size(); jj++) {

                torsionTorsionGrids[ii][kk][jj].resize(grid[kk][jj].size());
                if (reorder) {
                    for (unsigned int ll = 0; ll < grid[ll][jj].size(); ll++) {
                        torsionTorsionGrids[ii][kk][jj][ll] = reorderedGrid[kk][jj][ll];
                    }
                } else {
                    for (unsigned int ll = 0; ll < grid[ll][jj].size(); ll++) {
                        torsionTorsionGrids[ii][kk][jj][ll] = grid[kk][jj][ll];
                    }
                }
            }
        }
    }
}

double ReferenceCalcAmoebaTorsionTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {

    vector<Vec3>& posData   = extractPositions(context);
    vector<Vec3>& forceData = extractForces(context);
    AmoebaReferenceTorsionTorsionForce amoebaReferenceTorsionTorsionForce;
    if (usePeriodic)
        amoebaReferenceTorsionTorsionForce.setPeriodic(extractBoxVectors(context));
    double energy = amoebaReferenceTorsionTorsionForce.calculateForceAndEnergy(numTorsionTorsions, posData,
                                                                                         particle1, particle2, particle3, particle4, particle5,
                                                                                         chiralCheckAtom, gridIndices, torsionTorsionGrids, forceData);
    return static_cast<double>(energy);
}

/* -------------------------------------------------------------------------- *
 *                             AmoebaMultipole                                *
 * -------------------------------------------------------------------------- */

ReferenceCalcAmoebaMultipoleForceKernel::ReferenceCalcAmoebaMultipoleForceKernel(const std::string& name, const Platform& platform, const System& system) :
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
    double totalCharge   = 0.0;
    for (int ii = 0; ii < numMultipoles; ii++) {

        // multipoles

        int axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY;
        double charge, tholeD, dampingFactorD, polarityD;
        std::vector<double> dipolesD;
        std::vector<double> quadrupolesD;
        force.getMultipoleParameters(ii, charge, dipolesD, quadrupolesD, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY,
                                     tholeD, dampingFactorD, polarityD);

        totalCharge                       += charge;
        axisTypes[ii]                      = axisType;
        multipoleAtomZs[ii]                = multipoleAtomZ;
        multipoleAtomXs[ii]                = multipoleAtomX;
        multipoleAtomYs[ii]                = multipoleAtomY;

        charges[ii]                        = charge;
        tholes[ii]                         = tholeD;
        dampingFactors[ii]                 = dampingFactorD;
        polarity[ii]                       = polarityD;

        dipoles[dipoleIndex++]             = dipolesD[0];
        dipoles[dipoleIndex++]             = dipolesD[1];
        dipoles[dipoleIndex++]             = dipolesD[2];
        
        quadrupoles[quadrupoleIndex++]     = quadrupolesD[0];
        quadrupoles[quadrupoleIndex++]     = quadrupolesD[1];
        quadrupoles[quadrupoleIndex++]     = quadrupolesD[2];
        quadrupoles[quadrupoleIndex++]     = quadrupolesD[3];
        quadrupoles[quadrupoleIndex++]     = quadrupolesD[4];
        quadrupoles[quadrupoleIndex++]     = quadrupolesD[5];
        quadrupoles[quadrupoleIndex++]     = quadrupolesD[6];
        quadrupoles[quadrupoleIndex++]     = quadrupolesD[7];
        quadrupoles[quadrupoleIndex++]     = quadrupolesD[8];

        // covalent info

        std::vector< std::vector<int> > covalentLists;
        force.getCovalentMaps(ii, covalentLists);
        multipoleAtomCovalentInfo[ii] = covalentLists;

    }

    polarizationType = force.getPolarizationType();
    if (polarizationType == AmoebaMultipoleForce::Mutual) {
        mutualInducedMaxIterations = force.getMutualInducedMaxIterations();
        mutualInducedTargetEpsilon = force.getMutualInducedTargetEpsilon();
    } else if (polarizationType == AmoebaMultipoleForce::Extrapolated) {
        extrapolationCoefficients = force.getExtrapolationCoefficients();
    }

    // PME

    nonbondedMethod  = force.getNonbondedMethod();
    if (nonbondedMethod == AmoebaMultipoleForce::PME) {
        usePme     = true;
        pmeGridDimension.resize(3);
        force.getPMEParameters(alphaEwald, pmeGridDimension[0], pmeGridDimension[1], pmeGridDimension[2]);
        cutoffDistance = force.getCutoffDistance();
        if (pmeGridDimension[0] == 0 || alphaEwald == 0.0) {
            NonbondedForce nb;
            nb.setEwaldErrorTolerance(force.getEwaldErrorTolerance());
            nb.setCutoffDistance(force.getCutoffDistance());
            int gridSizeX, gridSizeY, gridSizeZ;
            NonbondedForceImpl::calcPMEParameters(system, nb, alphaEwald, gridSizeX, gridSizeY, gridSizeZ, false);
            pmeGridDimension[0] = gridSizeX;
            pmeGridDimension[1] = gridSizeY;
            pmeGridDimension[2] = gridSizeZ;
        }    
    } else {
        usePme = false;
    }
    return;
}

AmoebaReferenceMultipoleForce* ReferenceCalcAmoebaMultipoleForceKernel::setupAmoebaReferenceMultipoleForce(ContextImpl& context)
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
    if (gkKernel) {

        // amoebaReferenceGeneralizedKirkwoodForce is deleted in AmoebaReferenceGeneralizedKirkwoodMultipoleForce
        // destructor

        AmoebaReferenceGeneralizedKirkwoodForce* amoebaReferenceGeneralizedKirkwoodForce = new AmoebaReferenceGeneralizedKirkwoodForce();
        amoebaReferenceGeneralizedKirkwoodForce->setNumParticles(gkKernel->getNumParticles());
        amoebaReferenceGeneralizedKirkwoodForce->setSoluteDielectric(gkKernel->getSoluteDielectric());
        amoebaReferenceGeneralizedKirkwoodForce->setSolventDielectric(gkKernel->getSolventDielectric());
        amoebaReferenceGeneralizedKirkwoodForce->setDielectricOffset(gkKernel->getDielectricOffset());
        amoebaReferenceGeneralizedKirkwoodForce->setProbeRadius(gkKernel->getProbeRadius());
        amoebaReferenceGeneralizedKirkwoodForce->setSurfaceAreaFactor(gkKernel->getSurfaceAreaFactor());
        amoebaReferenceGeneralizedKirkwoodForce->setIncludeCavityTerm(gkKernel->getIncludeCavityTerm());
        amoebaReferenceGeneralizedKirkwoodForce->setDirectPolarization(gkKernel->getDirectPolarization());

        vector<double> parameters; 
        gkKernel->getAtomicRadii(parameters);
        amoebaReferenceGeneralizedKirkwoodForce->setAtomicRadii(parameters);

        gkKernel->getScaleFactors(parameters);
        amoebaReferenceGeneralizedKirkwoodForce->setScaleFactors(parameters);

        gkKernel->getCharges(parameters);
        amoebaReferenceGeneralizedKirkwoodForce->setCharges(parameters);

        // calculate Grycuk Born radii

        vector<Vec3>& posData   = extractPositions(context);
        amoebaReferenceGeneralizedKirkwoodForce->calculateGrycukBornRadii(posData);

        amoebaReferenceMultipoleForce = new AmoebaReferenceGeneralizedKirkwoodMultipoleForce(amoebaReferenceGeneralizedKirkwoodForce);

    } else if (usePme) {

        AmoebaReferencePmeMultipoleForce* amoebaReferencePmeMultipoleForce = new AmoebaReferencePmeMultipoleForce();
        amoebaReferencePmeMultipoleForce->setAlphaEwald(alphaEwald);
        amoebaReferencePmeMultipoleForce->setCutoffDistance(cutoffDistance);
        amoebaReferencePmeMultipoleForce->setPmeGridDimensions(pmeGridDimension);
        Vec3* boxVectors = extractBoxVectors(context);
        double minAllowedSize = 1.999999*cutoffDistance;
        if (boxVectors[0][0] < minAllowedSize || boxVectors[1][1] < minAllowedSize || boxVectors[2][2] < minAllowedSize) {
            throw OpenMMException("The periodic box size has decreased to less than twice the nonbonded cutoff.");
        }
        amoebaReferencePmeMultipoleForce->setPeriodicBoxSize(boxVectors);
        amoebaReferenceMultipoleForce = static_cast<AmoebaReferenceMultipoleForce*>(amoebaReferencePmeMultipoleForce);

    } else {
         amoebaReferenceMultipoleForce = new AmoebaReferenceMultipoleForce(AmoebaReferenceMultipoleForce::NoCutoff);
    }

    // set polarization type

    if (polarizationType == AmoebaMultipoleForce::Mutual) {
        amoebaReferenceMultipoleForce->setPolarizationType(AmoebaReferenceMultipoleForce::Mutual);
        amoebaReferenceMultipoleForce->setMutualInducedDipoleTargetEpsilon(mutualInducedTargetEpsilon);
        amoebaReferenceMultipoleForce->setMaximumMutualInducedDipoleIterations(mutualInducedMaxIterations);
    } else if (polarizationType == AmoebaMultipoleForce::Direct) {
        amoebaReferenceMultipoleForce->setPolarizationType(AmoebaReferenceMultipoleForce::Direct);
    } else if (polarizationType == AmoebaMultipoleForce::Extrapolated) {
        amoebaReferenceMultipoleForce->setPolarizationType(AmoebaReferenceMultipoleForce::Extrapolated);
        amoebaReferenceMultipoleForce->setExtrapolationCoefficients(extrapolationCoefficients);
    } else {
        throw OpenMMException("Polarization type not recognzied.");
    }

    return amoebaReferenceMultipoleForce;

}

double ReferenceCalcAmoebaMultipoleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {

    AmoebaReferenceMultipoleForce* amoebaReferenceMultipoleForce = setupAmoebaReferenceMultipoleForce(context);

    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& forceData = extractForces(context);
    double energy = amoebaReferenceMultipoleForce->calculateForceAndEnergy(posData, charges, dipoles, quadrupoles, tholes,
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
    
    AmoebaReferenceMultipoleForce* amoebaReferenceMultipoleForce = setupAmoebaReferenceMultipoleForce(context);
    vector<Vec3>& posData = extractPositions(context);
    
    // Retrieve the induced dipoles.
    
    vector<Vec3> inducedDipoles;
    amoebaReferenceMultipoleForce->calculateInducedDipoles(posData, charges, dipoles, quadrupoles, tholes,
            dampingFactors, polarity, axisTypes, multipoleAtomZs, multipoleAtomXs, multipoleAtomYs, multipoleAtomCovalentInfo, inducedDipoles);
    for (int i = 0; i < numParticles; i++)
        outputDipoles[i] = inducedDipoles[i];
    delete amoebaReferenceMultipoleForce;
}

void ReferenceCalcAmoebaMultipoleForceKernel::getLabFramePermanentDipoles(ContextImpl& context, vector<Vec3>& outputDipoles) {
    int numParticles = context.getSystem().getNumParticles();
    outputDipoles.resize(numParticles);

    // Create an AmoebaReferenceMultipoleForce to do the calculation.
    
    AmoebaReferenceMultipoleForce* amoebaReferenceMultipoleForce = setupAmoebaReferenceMultipoleForce(context);
    vector<Vec3>& posData = extractPositions(context);
    
    // Retrieve the permanent dipoles in the lab frame.
    
    vector<Vec3> labFramePermanentDipoles;
    amoebaReferenceMultipoleForce->calculateLabFramePermanentDipoles(posData, charges, dipoles, quadrupoles, tholes, 
            dampingFactors, polarity, axisTypes, multipoleAtomZs, multipoleAtomXs, multipoleAtomYs, multipoleAtomCovalentInfo, labFramePermanentDipoles);
    for (int i = 0; i < numParticles; i++)
        outputDipoles[i] = labFramePermanentDipoles[i];
    delete amoebaReferenceMultipoleForce;
}


void ReferenceCalcAmoebaMultipoleForceKernel::getTotalDipoles(ContextImpl& context, vector<Vec3>& outputDipoles) {
    int numParticles = context.getSystem().getNumParticles();
    outputDipoles.resize(numParticles);

    // Create an AmoebaReferenceMultipoleForce to do the calculation.
    
    AmoebaReferenceMultipoleForce* amoebaReferenceMultipoleForce = setupAmoebaReferenceMultipoleForce(context);
    vector<Vec3>& posData = extractPositions(context);
    
    // Retrieve the permanent dipoles in the lab frame.
    
    vector<Vec3> totalDipoles;
    amoebaReferenceMultipoleForce->calculateTotalDipoles(posData, charges, dipoles, quadrupoles, tholes,
            dampingFactors, polarity, axisTypes, multipoleAtomZs, multipoleAtomXs, multipoleAtomYs, multipoleAtomCovalentInfo, totalDipoles);

    for (int i = 0; i < numParticles; i++)
        outputDipoles[i] = totalDipoles[i];
    delete amoebaReferenceMultipoleForce;
}



void ReferenceCalcAmoebaMultipoleForceKernel::getElectrostaticPotential(ContextImpl& context, const std::vector< Vec3 >& inputGrid,
                                                                        std::vector< double >& outputElectrostaticPotential) {

    AmoebaReferenceMultipoleForce* amoebaReferenceMultipoleForce = setupAmoebaReferenceMultipoleForce(context);
    vector<Vec3>& posData                                     = extractPositions(context);
    vector<Vec3> grid(inputGrid.size());
    vector<double> potential(inputGrid.size());
    for (unsigned int ii = 0; ii < inputGrid.size(); ii++) {
        grid[ii] = inputGrid[ii];
    }
    amoebaReferenceMultipoleForce->calculateElectrostaticPotential(posData, charges, dipoles, quadrupoles, tholes,
                                                                   dampingFactors, polarity, axisTypes, 
                                                                   multipoleAtomZs, multipoleAtomXs, multipoleAtomYs,
                                                                   multipoleAtomCovalentInfo, grid, potential);

    outputElectrostaticPotential.resize(inputGrid.size());
    for (unsigned int ii = 0; ii < inputGrid.size(); ii++) {
        outputElectrostaticPotential[ii] = potential[ii];
    }

    delete amoebaReferenceMultipoleForce;
}

void ReferenceCalcAmoebaMultipoleForceKernel::getSystemMultipoleMoments(ContextImpl& context, std::vector< double >& outputMultipoleMoments) {

    // retrieve masses

    const System& system             = context.getSystem();
    vector<double> masses;
    for (int i = 0; i <  system.getNumParticles(); ++i) {
        masses.push_back(system.getParticleMass(i));
    }    

    AmoebaReferenceMultipoleForce* amoebaReferenceMultipoleForce = setupAmoebaReferenceMultipoleForce(context);
    vector<Vec3>& posData                                     = extractPositions(context);
    amoebaReferenceMultipoleForce->calculateAmoebaSystemMultipoleMoments(masses, posData, charges, dipoles, quadrupoles, tholes,
                                                                         dampingFactors, polarity, axisTypes, 
                                                                         multipoleAtomZs, multipoleAtomXs, multipoleAtomYs,
                                                                         multipoleAtomCovalentInfo, outputMultipoleMoments);

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
        charges[i] = charge;
        tholes[i] = tholeD;
        dampingFactors[i] = dampingFactorD;
        polarity[i] = polarityD;
        dipoles[dipoleIndex++] = dipolesD[0];
        dipoles[dipoleIndex++] = dipolesD[1];
        dipoles[dipoleIndex++] = dipolesD[2];
        quadrupoles[quadrupoleIndex++] = quadrupolesD[0];
        quadrupoles[quadrupoleIndex++] = quadrupolesD[1];
        quadrupoles[quadrupoleIndex++] = quadrupolesD[2];
        quadrupoles[quadrupoleIndex++] = quadrupolesD[3];
        quadrupoles[quadrupoleIndex++] = quadrupolesD[4];
        quadrupoles[quadrupoleIndex++] = quadrupolesD[5];
        quadrupoles[quadrupoleIndex++] = quadrupolesD[6];
        quadrupoles[quadrupoleIndex++] = quadrupolesD[7];
        quadrupoles[quadrupoleIndex++] = quadrupolesD[8];
    }
}

void ReferenceCalcAmoebaMultipoleForceKernel::getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    if (!usePme)
        throw OpenMMException("getPMEParametersInContext: This Context is not using PME");
    alpha = alphaEwald;
    nx = pmeGridDimension[0];
    ny = pmeGridDimension[1];
    nz = pmeGridDimension[2];
}

/* -------------------------------------------------------------------------- *
 *                       AmoebaGeneralizedKirkwood                            *
 * -------------------------------------------------------------------------- */

ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel::ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel(const std::string& name, const Platform& platform, const System& system) :
           CalcAmoebaGeneralizedKirkwoodForceKernel(name, platform), system(system) {
}

ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel::~ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel() {
}

int ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel::getNumParticles() const {
    return numParticles;
}

int ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel::getIncludeCavityTerm() const {
    return includeCavityTerm;
}

int ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel::getDirectPolarization() const {
    return directPolarization;
}

double ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel::getSoluteDielectric() const {
    return soluteDielectric;
}

double ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel::getSolventDielectric() const {
    return solventDielectric;
}

double ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel::getDielectricOffset() const {
    return dielectricOffset;
}

double ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel::getProbeRadius() const {
    return probeRadius;
}

double ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel::getSurfaceAreaFactor() const {
    return surfaceAreaFactor;
}

void ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel::getAtomicRadii(vector<double>& outputAtomicRadii) const {
    outputAtomicRadii.resize(atomicRadii.size());
    copy(atomicRadii.begin(), atomicRadii.end(), outputAtomicRadii.begin());
}

void ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel::getScaleFactors(vector<double>& outputScaleFactors) const {
    outputScaleFactors.resize(scaleFactors.size());
    copy(scaleFactors.begin(), scaleFactors.end(), outputScaleFactors.begin());
}

void ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel::getCharges(vector<double>& outputCharges) const {
    outputCharges.resize(charges.size());
    copy(charges.begin(), charges.end(), outputCharges.begin());
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

    if (amoebaMultipoleForce->getNonbondedMethod() != AmoebaMultipoleForce::NoCutoff) {
        throw OpenMMException("AmoebaGeneralizedKirkwoodForce requires the AmoebaMultipoleForce use the NoCutoff nonbonded method.");
    }

    numParticles = system.getNumParticles();

    for (int ii = 0; ii < numParticles; ii++) {

        double particleCharge, particleRadius, scalingFactor;
        force.getParticleParameters(ii, particleCharge, particleRadius, scalingFactor);
        atomicRadii.push_back(particleRadius);
        scaleFactors.push_back(scalingFactor);
        charges.push_back(particleCharge);

        // Make sure the charge matches the one specified by the AmoebaMultipoleForce.

        double charge2, thole, damping, polarity;
        int axisType, atomX, atomY, atomZ;
        vector<double> dipole, quadrupole;
        amoebaMultipoleForce->getMultipoleParameters(ii, charge2, dipole, quadrupole, axisType, atomZ, atomX, atomY, thole, damping, polarity);
        if (particleCharge != charge2) {
            throw OpenMMException("AmoebaGeneralizedKirkwoodForce and AmoebaMultipoleForce must specify the same charge for every atom.");
        }

    }   
    includeCavityTerm  = force.getIncludeCavityTerm();
    soluteDielectric   = force.getSoluteDielectric();
    solventDielectric  = force.getSolventDielectric();
    dielectricOffset   = 0.009;
    probeRadius        = force.getProbeRadius(), 
    surfaceAreaFactor  = force.getSurfaceAreaFactor(); 
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

ReferenceCalcAmoebaVdwForceKernel::ReferenceCalcAmoebaVdwForceKernel(const std::string& name, const Platform& platform, const System& system) :
       CalcAmoebaVdwForceKernel(name, platform), system(system) {
    useCutoff = 0;
    usePBC = 0;
    cutoff = 1.0e+10;
    neighborList = NULL;
}

ReferenceCalcAmoebaVdwForceKernel::~ReferenceCalcAmoebaVdwForceKernel() {
    if (neighborList)
        delete neighborList;
}

void ReferenceCalcAmoebaVdwForceKernel::initialize(const System& system, const AmoebaVdwForce& force) {

    // per-particle parameters

    numParticles = system.getNumParticles();
    useCutoff              = (force.getNonbondedMethod() != AmoebaVdwForce::NoCutoff);
    usePBC                 = (force.getNonbondedMethod() == AmoebaVdwForce::CutoffPeriodic);
    cutoff                 = force.getCutoffDistance();
    neighborList           = useCutoff ? new NeighborList() : NULL;
    dispersionCoefficient  = force.getUseDispersionCorrection() ?  AmoebaVdwForceImpl::calcDispersionCorrection(system, force) : 0.0;
    vdwForce.initialize(force);
}

double ReferenceCalcAmoebaVdwForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {

    vector<Vec3>& posData   = extractPositions(context);
    vector<Vec3>& forceData = extractForces(context);
    double energy;
    double lambda = context.getParameter(AmoebaVdwForce::Lambda());
    if (useCutoff) {
        computeNeighborListVoxelHash(*neighborList, numParticles, posData, vdwForce.getExclusions(), extractBoxVectors(context), usePBC, cutoff, 0.0);
        if (usePBC) {
            Vec3* boxVectors = extractBoxVectors(context);
            double minAllowedSize = 1.999999*cutoff;
            if (boxVectors[0][0] < minAllowedSize || boxVectors[1][1] < minAllowedSize || boxVectors[2][2] < minAllowedSize)
                throw OpenMMException("The periodic box size has decreased to less than twice the cutoff.");
            vdwForce.setPeriodicBox(boxVectors);
            energy  = vdwForce.calculateForceAndEnergy(numParticles, lambda, posData, *neighborList, forceData);
            energy += dispersionCoefficient/(boxVectors[0][0]*boxVectors[1][1]*boxVectors[2][2]);
        }
    }
    else
        energy = vdwForce.calculateForceAndEnergy(numParticles, lambda, posData, forceData);
    return static_cast<double>(energy);
}

void ReferenceCalcAmoebaVdwForceKernel::copyParametersToContext(ContextImpl& context, const AmoebaVdwForce& force) {
    if (numParticles != force.getNumParticles())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");
    vdwForce.initialize(force);
}

/* -------------------------------------------------------------------------- *
 *                           AmoebaWcaDispersion                              *
 * -------------------------------------------------------------------------- */

ReferenceCalcAmoebaWcaDispersionForceKernel::ReferenceCalcAmoebaWcaDispersionForceKernel(const std::string& name, const Platform& platform, const System& system) :
           CalcAmoebaWcaDispersionForceKernel(name, platform), system(system) {
}

ReferenceCalcAmoebaWcaDispersionForceKernel::~ReferenceCalcAmoebaWcaDispersionForceKernel() {
}

void ReferenceCalcAmoebaWcaDispersionForceKernel::initialize(const System& system, const AmoebaWcaDispersionForce& force) {

    // per-particle parameters

    numParticles = system.getNumParticles();
    radii.resize(numParticles);
    epsilons.resize(numParticles);
    for (int ii = 0; ii < numParticles; ii++) {

        double radius, epsilon;
        force.getParticleParameters(ii, radius, epsilon);

        radii[ii] = radius;
        epsilons[ii] = epsilon;
    }   

    totalMaximumDispersionEnergy = AmoebaWcaDispersionForceImpl::getTotalMaximumDispersionEnergy(force);

    epso    = force.getEpso();
    epsh    = force.getEpsh();
    rmino   = force.getRmino();
    rminh   = force.getRminh();
    awater  = force.getAwater();
    shctd   = force.getShctd();
    dispoff = force.getDispoff();
    slevy   = force.getSlevy();
}

double ReferenceCalcAmoebaWcaDispersionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& forceData = extractForces(context);
    AmoebaReferenceWcaDispersionForce amoebaReferenceWcaDispersionForce(epso, epsh, rmino, rminh, awater, shctd, dispoff, slevy);
    double energy = amoebaReferenceWcaDispersionForce.calculateForceAndEnergy(numParticles, posData, radii, epsilons, totalMaximumDispersionEnergy, forceData);
    return static_cast<double>(energy);
}

void ReferenceCalcAmoebaWcaDispersionForceKernel::copyParametersToContext(ContextImpl& context, const AmoebaWcaDispersionForce& force) {
    if (numParticles != force.getNumParticles())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");

    // Record the values.

    for (int i = 0; i < numParticles; ++i) {
        double radius, epsilon;
        force.getParticleParameters(i, radius, epsilon);
        radii[i] = radius;
        epsilons[i] = epsilon;
    }
    totalMaximumDispersionEnergy = AmoebaWcaDispersionForceImpl::getTotalMaximumDispersionEnergy(force);
}


/* -------------------------------------------------------------------------- *
 *                              HippoNonbonded                                *
 * -------------------------------------------------------------------------- */

ReferenceCalcHippoNonbondedForceKernel::ReferenceCalcHippoNonbondedForceKernel(const std::string& name, const Platform& platform, const System& system) :
         CalcHippoNonbondedForceKernel(name, platform), ixn(NULL) {
}

ReferenceCalcHippoNonbondedForceKernel::~ReferenceCalcHippoNonbondedForceKernel() {
    if (ixn != NULL)
        delete ixn;
}

void ReferenceCalcHippoNonbondedForceKernel::initialize(const System& system, const HippoNonbondedForce& force) {
    numParticles = force.getNumParticles();
    if (force.getNonbondedMethod() == HippoNonbondedForce::PME)
        ixn = new AmoebaReferencePmeHippoNonbondedForce(force, system);
    else
        ixn = new AmoebaReferenceHippoNonbondedForce(force);
}

void ReferenceCalcHippoNonbondedForceKernel::setupAmoebaReferenceHippoNonbondedForce(ContextImpl& context) {
    if (ixn->getNonbondedMethod() == HippoNonbondedForce::PME) {
        AmoebaReferencePmeHippoNonbondedForce* force = dynamic_cast<AmoebaReferencePmeHippoNonbondedForce*>(ixn);
        Vec3* boxVectors = extractBoxVectors(context);
        double minAllowedSize = 1.999999*force->getCutoffDistance();
        if (boxVectors[0][0] < minAllowedSize || boxVectors[1][1] < minAllowedSize || boxVectors[2][2] < minAllowedSize)
            throw OpenMMException("The periodic box size has decreased to less than twice the nonbonded cutoff.");
        force->setPeriodicBoxSize(boxVectors);
    }
}

double ReferenceCalcHippoNonbondedForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {

    setupAmoebaReferenceHippoNonbondedForce(context);

    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& forceData = extractForces(context);
    return ixn->calculateForceAndEnergy(posData, forceData);
}

void ReferenceCalcHippoNonbondedForceKernel::getInducedDipoles(ContextImpl& context, vector<Vec3>& outputDipoles) {
    outputDipoles.resize(numParticles);
    setupAmoebaReferenceHippoNonbondedForce(context);
    vector<Vec3>& posData = extractPositions(context);
    
    // Retrieve the induced dipoles.
    
    vector<Vec3> inducedDipoles;
    ixn->calculateInducedDipoles(posData, inducedDipoles);
    for (int i = 0; i < numParticles; i++)
        outputDipoles[i] = inducedDipoles[i];
}

void ReferenceCalcHippoNonbondedForceKernel::getLabFramePermanentDipoles(ContextImpl& context, vector<Vec3>& outputDipoles) {
    outputDipoles.resize(numParticles);
    setupAmoebaReferenceHippoNonbondedForce(context);
    vector<Vec3>& posData = extractPositions(context);
    
    // Retrieve the permanent dipoles in the lab frame.
    
    vector<Vec3> labFramePermanentDipoles;
    ixn->calculateLabFramePermanentDipoles(posData, labFramePermanentDipoles);
    for (int i = 0; i < numParticles; i++)
        outputDipoles[i] = labFramePermanentDipoles[i];
}

void ReferenceCalcHippoNonbondedForceKernel::copyParametersToContext(ContextImpl& context, const HippoNonbondedForce& force) {
    if (numParticles != force.getNumParticles())
        throw OpenMMException("updateParametersInContext: The number of multipoles has changed");
    delete ixn;
    ixn = NULL;
    if (force.getNonbondedMethod() == HippoNonbondedForce::PME)
        ixn = new AmoebaReferencePmeHippoNonbondedForce(force, context.getSystem());
    else
        ixn = new AmoebaReferenceHippoNonbondedForce(force);
}

void ReferenceCalcHippoNonbondedForceKernel::getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    if (ixn->getNonbondedMethod() != HippoNonbondedForce::PME)
        throw OpenMMException("getPMEParametersInContext: This Context is not using PME");
    AmoebaReferencePmeHippoNonbondedForce* force = dynamic_cast<AmoebaReferencePmeHippoNonbondedForce*>(ixn);
    alpha = force->getAlphaEwald();
    vector<int> dim;
    force->getPmeGridDimensions(dim);
    nx = dim[0];
    ny = dim[1];
    nz = dim[2];
}

void ReferenceCalcHippoNonbondedForceKernel::getDPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    if (ixn->getNonbondedMethod() != HippoNonbondedForce::PME)
        throw OpenMMException("getDPMEParametersInContext: This Context is not using PME");
    AmoebaReferencePmeHippoNonbondedForce* force = dynamic_cast<AmoebaReferencePmeHippoNonbondedForce*>(ixn);
    alpha = force->getDispersionAlphaEwald();
    vector<int> dim;
    force->getDispersionPmeGridDimensions(dim);
    nx = dim[0];
    ny = dim[1];
    nz = dim[2];
}
