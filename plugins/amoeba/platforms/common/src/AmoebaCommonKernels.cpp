/* -------------------------------------------------------------------------- *
 *                               OpenMMAmoeba                                 *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2021 Stanford University and the Authors.      *
 * Authors: Peter Eastman, Mark Friedrichs                                    *
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

#ifdef WIN32
  #define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include "AmoebaCommonKernels.h"
#include "CommonAmoebaKernelSources.h"
#include "openmm/common/ContextSelector.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/AmoebaGeneralizedKirkwoodForceImpl.h"
#include "openmm/internal/AmoebaMultipoleForceImpl.h"
#include "openmm/internal/AmoebaWcaDispersionForceImpl.h"
#include "openmm/internal/AmoebaTorsionTorsionForceImpl.h"
#include "openmm/internal/AmoebaVdwForceImpl.h"
#include "openmm/internal/NonbondedForceImpl.h"
#include "CommonKernelSources.h"
#include "SimTKOpenMMRealType.h"
#include "jama_lu.h"

#include <algorithm>
#include <cmath>
#ifdef _MSC_VER
#include <windows.h>
#endif

using namespace OpenMM;
using namespace std;

static void setPeriodicBoxArgs(ComputeContext& cc, ComputeKernel kernel, int index) {
    Vec3 a, b, c;
    cc.getPeriodicBoxVectors(a, b, c);
    if (cc.getUseDoublePrecision()) {
        kernel->setArg(index++, mm_double4(a[0], b[1], c[2], 0.0));
        kernel->setArg(index++, mm_double4(1.0/a[0], 1.0/b[1], 1.0/c[2], 0.0));
        kernel->setArg(index++, mm_double4(a[0], a[1], a[2], 0.0));
        kernel->setArg(index++, mm_double4(b[0], b[1], b[2], 0.0));
        kernel->setArg(index, mm_double4(c[0], c[1], c[2], 0.0));
    }
    else {
        kernel->setArg(index++, mm_float4((float) a[0], (float) b[1], (float) c[2], 0.0f));
        kernel->setArg(index++, mm_float4(1.0f/(float) a[0], 1.0f/(float) b[1], 1.0f/(float) c[2], 0.0f));
        kernel->setArg(index++, mm_float4((float) a[0], (float) a[1], (float) a[2], 0.0f));
        kernel->setArg(index++, mm_float4((float) b[0], (float) b[1], (float) b[2], 0.0f));
        kernel->setArg(index, mm_float4((float) c[0], (float) c[1], (float) c[2], 0.0f));
    }
}

/* -------------------------------------------------------------------------- *
 *                           AmoebaTorsionTorsion                             *
 * -------------------------------------------------------------------------- */

class CommonCalcAmoebaTorsionTorsionForceKernel::ForceInfo : public ComputeForceInfo {
public:
    ForceInfo(const AmoebaTorsionTorsionForce& force) : force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumTorsionTorsions();
    }
    void getParticlesInGroup(int index, std::vector<int>& particles) {
        int particle1, particle2, particle3, particle4, particle5, chiralCheckAtomIndex, gridIndex;
        force.getTorsionTorsionParameters(index, particle1, particle2, particle3, particle4, particle5, chiralCheckAtomIndex, gridIndex);
        particles.resize(5);
        particles[0] = particle1;
        particles[1] = particle2;
        particles[2] = particle3;
        particles[3] = particle4;
        particles[4] = particle5;
    }
    bool areGroupsIdentical(int group1, int group2) {
        int particle1, particle2, particle3, particle4, particle5;
        int chiral1, chiral2, grid1, grid2;
        force.getTorsionTorsionParameters(group1, particle1, particle2, particle3, particle4, particle5, chiral1, grid1);
        force.getTorsionTorsionParameters(group2, particle1, particle2, particle3, particle4, particle5, chiral2, grid2);
        return (grid1 == grid2);
    }
private:
    const AmoebaTorsionTorsionForce& force;
};

CommonCalcAmoebaTorsionTorsionForceKernel::CommonCalcAmoebaTorsionTorsionForceKernel(const std::string& name, const Platform& platform, ComputeContext& cc, const System& system) :
                CalcAmoebaTorsionTorsionForceKernel(name, platform), cc(cc), system(system) {
}

void CommonCalcAmoebaTorsionTorsionForceKernel::initialize(const System& system, const AmoebaTorsionTorsionForce& force) {
    ContextSelector selector(cc);
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumTorsionTorsions()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumTorsionTorsions()/numContexts;
    numTorsionTorsions = endIndex-startIndex;
    if (numTorsionTorsions == 0)
        return;
    
    // Record torsion parameters.
    
    vector<vector<int> > atoms(numTorsionTorsions, vector<int>(5));
    vector<mm_int2> torsionParamsVec(numTorsionTorsions);
    torsionParams.initialize<mm_int2>(cc, numTorsionTorsions, "torsionTorsionParams");
    for (int i = 0; i < numTorsionTorsions; i++)
        force.getTorsionTorsionParameters(startIndex+i, atoms[i][0], atoms[i][1], atoms[i][2], atoms[i][3], atoms[i][4], torsionParamsVec[i].x, torsionParamsVec[i].y);
    torsionParams.upload(torsionParamsVec);
    
    // Record the grids.
    
    vector<mm_float4> gridValuesVec;
    vector<mm_float4> gridParamsVec;
    for (int i = 0; i < force.getNumTorsionTorsionGrids(); i++) {
        const TorsionTorsionGrid& initialGrid = force.getTorsionTorsionGrid(i);

        // check if grid needs to be reordered: x-angle should be 'slow' index

        bool reordered = false;
        TorsionTorsionGrid reorderedGrid;
        if (initialGrid[0][0][0] != initialGrid[0][1][0]) {
            AmoebaTorsionTorsionForceImpl::reorderGrid(initialGrid, reorderedGrid);
            reordered = true;
        }
        const TorsionTorsionGrid& grid = (reordered ? reorderedGrid : initialGrid);
        float range = grid[0][grid[0].size()-1][1] - grid[0][0][1];
        gridParamsVec.push_back(mm_float4(gridValuesVec.size(), grid[0][0][0], range/(grid.size()-1), grid.size()));
        for (int j = 0; j < grid.size(); j++)
            for (int k = 0; k < grid[j].size(); k++)
                gridValuesVec.push_back(mm_float4((float) grid[j][k][2], (float) grid[j][k][3], (float) grid[j][k][4], (float) grid[j][k][5]));
    }
    gridValues.initialize<mm_float4>(cc, gridValuesVec.size(), "torsionTorsionGridValues");
    gridParams.initialize<mm_float4>(cc, gridParamsVec.size(), "torsionTorsionGridParams");
    gridValues.upload(gridValuesVec);
    gridParams.upload(gridParamsVec);
    map<string, string> replacements;
    replacements["APPLY_PERIODIC"] = (force.usesPeriodicBoundaryConditions() ? "1" : "0");
    replacements["GRID_VALUES"] = cc.getBondedUtilities().addArgument(gridValues, "float4");
    replacements["GRID_PARAMS"] = cc.getBondedUtilities().addArgument(gridParams, "float4");
    replacements["TORSION_PARAMS"] = cc.getBondedUtilities().addArgument(torsionParams, "int2");
    replacements["RAD_TO_DEG"] = cc.doubleToString(180/M_PI);
    cc.getBondedUtilities().addInteraction(atoms, cc.replaceStrings(CommonAmoebaKernelSources::amoebaTorsionTorsionForce, replacements), force.getForceGroup());
    cc.getBondedUtilities().addPrefixCode(CommonAmoebaKernelSources::bicubic);
    cc.addForce(new ForceInfo(force));
}

double CommonCalcAmoebaTorsionTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

/* -------------------------------------------------------------------------- *
 *                             AmoebaMultipole                                *
 * -------------------------------------------------------------------------- */

class CommonCalcAmoebaMultipoleForceKernel::ForceInfo : public ComputeForceInfo {
public:
    ForceInfo(const AmoebaMultipoleForce& force) : force(force) {
    }
    bool areParticlesIdentical(int particle1, int particle2) {
        double charge1, charge2, thole1, thole2, damping1, damping2, polarity1, polarity2;
        int axis1, axis2, multipole11, multipole12, multipole21, multipole22, multipole31, multipole32;
        vector<double> dipole1, dipole2, quadrupole1, quadrupole2;
        force.getMultipoleParameters(particle1, charge1, dipole1, quadrupole1, axis1, multipole11, multipole21, multipole31, thole1, damping1, polarity1);
        force.getMultipoleParameters(particle2, charge2, dipole2, quadrupole2, axis2, multipole12, multipole22, multipole32, thole2, damping2, polarity2);
        if (charge1 != charge2 || thole1 != thole2 || damping1 != damping2 || polarity1 != polarity2 || axis1 != axis2) {
            return false;
        }
        for (int i = 0; i < (int) dipole1.size(); ++i) {
            if (dipole1[i] != dipole2[i]) {
                return false;
            }
        }
        for (int i = 0; i < (int) quadrupole1.size(); ++i) {
            if (quadrupole1[i] != quadrupole2[i]) {
                return false;
            }
        }
        return true;
    }
    int getNumParticleGroups() {
        return 7*force.getNumMultipoles();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        int particle = index/7;
        int type = index-7*particle;
        force.getCovalentMap(particle, AmoebaMultipoleForce::CovalentType(type), particles);
    }
    bool areGroupsIdentical(int group1, int group2) {
        return ((group1%7) == (group2%7));
    }
private:
    const AmoebaMultipoleForce& force;
};

CommonCalcAmoebaMultipoleForceKernel::CommonCalcAmoebaMultipoleForceKernel(const std::string& name, const Platform& platform, ComputeContext& cc, const System& system) :
        CalcAmoebaMultipoleForceKernel(name, platform), cc(cc), system(system), hasInitializedScaleFactors(false), multipolesAreValid(false), hasCreatedEvent(false),
        gkKernel(NULL) {
}

CommonCalcAmoebaMultipoleForceKernel::~CommonCalcAmoebaMultipoleForceKernel() {
}

void CommonCalcAmoebaMultipoleForceKernel::initialize(const System& system, const AmoebaMultipoleForce& force) {
    ContextSelector selector(cc);

    // Initialize multipole parameters.

    numMultipoles = force.getNumMultipoles();
    ArrayInterface& posq = cc.getPosq();
    vector<mm_double4> temp(posq.getSize());
    mm_float4* posqf = (mm_float4*) &temp[0];
    mm_double4* posqd = (mm_double4*) &temp[0];
    vector<mm_float2> dampingAndTholeVec;
    vector<float> polarizabilityVec;
    vector<float> localDipolesVec;
    vector<float> localQuadrupolesVec;
    vector<mm_int4> multipoleParticlesVec;
    for (int i = 0; i < numMultipoles; i++) {
        double charge, thole, damping, polarity;
        int axisType, atomX, atomY, atomZ;
        vector<double> dipole, quadrupole;
        force.getMultipoleParameters(i, charge, dipole, quadrupole, axisType, atomZ, atomX, atomY, thole, damping, polarity);
        if (cc.getUseDoublePrecision())
            posqd[i] = mm_double4(0, 0, 0, charge);
        else
            posqf[i] = mm_float4(0, 0, 0, (float) charge);
        dampingAndTholeVec.push_back(mm_float2((float) damping, (float) thole));
        polarizabilityVec.push_back((float) polarity);
        multipoleParticlesVec.push_back(mm_int4(atomX, atomY, atomZ, axisType));
        for (int j = 0; j < 3; j++)
            localDipolesVec.push_back((float) dipole[j]);
        localQuadrupolesVec.push_back((float) quadrupole[0]);
        localQuadrupolesVec.push_back((float) quadrupole[1]);
        localQuadrupolesVec.push_back((float) quadrupole[2]);
        localQuadrupolesVec.push_back((float) quadrupole[4]);
        localQuadrupolesVec.push_back((float) quadrupole[5]);
    }
    hasQuadrupoles = false;
    for (auto q : localQuadrupolesVec)
        if (q != 0.0)
            hasQuadrupoles = true;
    int paddedNumAtoms = cc.getPaddedNumAtoms();
    for (int i = numMultipoles; i < paddedNumAtoms; i++) {
        dampingAndTholeVec.push_back(mm_float2(0, 0));
        polarizabilityVec.push_back(0);
        multipoleParticlesVec.push_back(mm_int4(0, 0, 0, 0));
        for (int j = 0; j < 3; j++)
            localDipolesVec.push_back(0);
        for (int j = 0; j < 5; j++)
            localQuadrupolesVec.push_back(0);
    }
    dampingAndThole.initialize<mm_float2>(cc, paddedNumAtoms, "dampingAndThole");
    polarizability.initialize<float>(cc, paddedNumAtoms, "polarizability");
    multipoleParticles.initialize<mm_int4>(cc, paddedNumAtoms, "multipoleParticles");
    localDipoles.initialize<float>(cc, 3*paddedNumAtoms, "localDipoles");
    localQuadrupoles.initialize<float>(cc, 5*paddedNumAtoms, "localQuadrupoles");
    lastPositions.initialize(cc, cc.getPosq().getSize(), cc.getPosq().getElementSize(), "lastPositions");
    dampingAndThole.upload(dampingAndTholeVec);
    polarizability.upload(polarizabilityVec);
    multipoleParticles.upload(multipoleParticlesVec);
    localDipoles.upload(localDipolesVec);
    localQuadrupoles.upload(localQuadrupolesVec);
    posq.upload(&temp[0]);
    
    // Create workspace arrays.
    
    polarizationType = force.getPolarizationType();
    int elementSize = (cc.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
    labDipoles.initialize(cc, 3*paddedNumAtoms, elementSize, "labDipoles");
    labQuadrupoles.initialize(cc, 5*paddedNumAtoms, elementSize, "labQuadrupoles");
    sphericalDipoles.initialize(cc, 3*paddedNumAtoms, elementSize, "sphericalDipoles");
    sphericalQuadrupoles.initialize(cc, 5*paddedNumAtoms, elementSize, "sphericalQuadrupoles");
    fracDipoles.initialize(cc, 3*paddedNumAtoms, elementSize, "fracDipoles");
    fracQuadrupoles.initialize(cc, 6*paddedNumAtoms, elementSize, "fracQuadrupoles");
    field.initialize(cc, 3*paddedNumAtoms, sizeof(long long), "field");
    fieldPolar.initialize(cc, 3*paddedNumAtoms, sizeof(long long), "fieldPolar");
    torque.initialize(cc, 3*paddedNumAtoms, sizeof(long long), "torque");
    inducedDipole.initialize(cc, 3*paddedNumAtoms, elementSize, "inducedDipole");
    inducedDipolePolar.initialize(cc, 3*paddedNumAtoms, elementSize, "inducedDipolePolar");
    if (polarizationType == AmoebaMultipoleForce::Mutual) {
        inducedDipoleErrors.initialize(cc, cc.getNumThreadBlocks(), sizeof(mm_float2), "inducedDipoleErrors");
        prevDipoles.initialize(cc, 3*numMultipoles*MaxPrevDIISDipoles, elementSize, "prevDipoles");
        prevDipolesPolar.initialize(cc, 3*numMultipoles*MaxPrevDIISDipoles, elementSize, "prevDipolesPolar");
        prevErrors.initialize(cc, 3*numMultipoles*MaxPrevDIISDipoles, elementSize, "prevErrors");
        diisMatrix.initialize(cc, MaxPrevDIISDipoles*MaxPrevDIISDipoles, elementSize, "diisMatrix");
        diisCoefficients.initialize(cc, MaxPrevDIISDipoles+1, sizeof(float), "diisCoefficients");
        syncEvent = cc.createEvent();
        hasCreatedEvent = true;
    }
    else if (polarizationType == AmoebaMultipoleForce::Extrapolated) {
        int numOrders = force.getExtrapolationCoefficients().size();
        extrapolatedDipole.initialize(cc, 3*numMultipoles*numOrders, elementSize, "extrapolatedDipole");
        extrapolatedDipolePolar.initialize(cc, 3*numMultipoles*numOrders, elementSize, "extrapolatedDipolePolar");
        inducedDipoleFieldGradient.initialize(cc, 6*paddedNumAtoms, sizeof(long long), "inducedDipoleFieldGradient");
        inducedDipoleFieldGradientPolar.initialize(cc, 6*paddedNumAtoms, sizeof(long long), "inducedDipoleFieldGradientPolar");
        extrapolatedDipoleFieldGradient.initialize(cc, 6*paddedNumAtoms*(numOrders-1), elementSize, "extrapolatedDipoleFieldGradient");
        extrapolatedDipoleFieldGradientPolar.initialize(cc, 6*paddedNumAtoms*(numOrders-1), elementSize, "extrapolatedDipoleFieldGradientPolar");
    }
    // The next two arrays will get resized once we know the number of exclusions.
    covalentFlags.initialize<mm_int2>(cc, 1, "covalentFlags");
    polarizationGroupFlags.initialize<unsigned int>(cc, 1, "polarizationGroupFlags");
    cc.addAutoclearBuffer(field);
    cc.addAutoclearBuffer(fieldPolar);
    cc.addAutoclearBuffer(torque);
    
    // Record which atoms should be flagged as exclusions based on covalent groups, and determine
    // the values for the covalent group flags.
    
    vector<vector<int> > exclusions(numMultipoles);
    for (int i = 0; i < numMultipoles; i++) {
        vector<int> atoms;
        set<int> allAtoms;
        allAtoms.insert(i);
        force.getCovalentMap(i, AmoebaMultipoleForce::Covalent12, atoms);
        allAtoms.insert(atoms.begin(), atoms.end());
        force.getCovalentMap(i, AmoebaMultipoleForce::Covalent13, atoms);
        allAtoms.insert(atoms.begin(), atoms.end());
        for (int atom : allAtoms)
            covalentFlagValues.push_back(mm_int4(i, atom, 0, 0));
        force.getCovalentMap(i, AmoebaMultipoleForce::Covalent14, atoms);
        allAtoms.insert(atoms.begin(), atoms.end());
        for (int atom : atoms)
            covalentFlagValues.push_back(mm_int4(i, atom, 1, 0));
        force.getCovalentMap(i, AmoebaMultipoleForce::Covalent15, atoms);
        for (int atom : atoms)
            covalentFlagValues.push_back(mm_int4(i, atom, 2, 0));
        allAtoms.insert(atoms.begin(), atoms.end());
        force.getCovalentMap(i, AmoebaMultipoleForce::PolarizationCovalent11, atoms);
        allAtoms.insert(atoms.begin(), atoms.end());
        exclusions[i].insert(exclusions[i].end(), allAtoms.begin(), allAtoms.end());

        // Workaround for bug in TINKER: if an atom is listed in both the PolarizationCovalent11
        // and PolarizationCovalent12 maps, the latter takes precedence.

        vector<int> atoms12;
        force.getCovalentMap(i, AmoebaMultipoleForce::PolarizationCovalent12, atoms12);
        for (int atom : atoms)
            if (find(atoms12.begin(), atoms12.end(), atom) == atoms12.end())
                polarizationFlagValues.push_back(mm_int2(i, atom));
    }
    set<pair<int, int> > tilesWithExclusions;
    for (int atom1 = 0; atom1 < (int) exclusions.size(); ++atom1) {
        int x = atom1/ComputeContext::TileSize;
        for (int atom2 : exclusions[atom1]) {
            int y = atom2/ComputeContext::TileSize;
            tilesWithExclusions.insert(make_pair(max(x, y), min(x, y)));
        }
    }
    
    // Record other options.
    
    if (polarizationType == AmoebaMultipoleForce::Mutual) {
        maxInducedIterations = force.getMutualInducedMaxIterations();
        inducedEpsilon = force.getMutualInducedTargetEpsilon();
    }
    else
        maxInducedIterations = 0;
    if (polarizationType != AmoebaMultipoleForce::Direct) {
        inducedField.initialize(cc, 3*paddedNumAtoms, sizeof(long long), "inducedField");
        inducedFieldPolar.initialize(cc, 3*paddedNumAtoms, sizeof(long long), "inducedFieldPolar");
    }
    usePME = (force.getNonbondedMethod() == AmoebaMultipoleForce::PME);
    
    // See whether there's an AmoebaGeneralizedKirkwoodForce in the System.

    const AmoebaGeneralizedKirkwoodForce* gk = NULL;
    for (int i = 0; i < system.getNumForces() && gk == NULL; i++)
        gk = dynamic_cast<const AmoebaGeneralizedKirkwoodForce*>(&system.getForce(i));
    double innerDielectric = (gk == NULL ? 1.0 : gk->getSoluteDielectric());
    
    // Create the kernels.

    bool useShuffle = false;//(cc.getComputeCapability() >= 3.0 && !cc.getUseDoublePrecision());
    double fixedThreadMemory = 19*elementSize+2*sizeof(float)+3*sizeof(int)/(double) cc.TileSize;
    double inducedThreadMemory = 15*elementSize+2*sizeof(float);
    if (polarizationType == AmoebaMultipoleForce::Extrapolated)
        inducedThreadMemory += 12*elementSize;
    double electrostaticsThreadMemory = 0;
    if (!useShuffle)
        fixedThreadMemory += 3*elementSize;
    map<string, string> defines;
    defines["NUM_ATOMS"] = cc.intToString(numMultipoles);
    defines["PADDED_NUM_ATOMS"] = cc.intToString(cc.getPaddedNumAtoms());
    defines["NUM_BLOCKS"] = cc.intToString(cc.getNumAtomBlocks());
    defines["ENERGY_SCALE_FACTOR"] = cc.doubleToString(ONE_4PI_EPS0/innerDielectric);
    if (polarizationType == AmoebaMultipoleForce::Direct)
        defines["DIRECT_POLARIZATION"] = "";
    else if (polarizationType == AmoebaMultipoleForce::Mutual)
        defines["MUTUAL_POLARIZATION"] = "";
    else if (polarizationType == AmoebaMultipoleForce::Extrapolated)
        defines["EXTRAPOLATED_POLARIZATION"] = "";
    if (useShuffle)
        defines["USE_SHUFFLE"] = "";
    if (hasQuadrupoles)
        defines["INCLUDE_QUADRUPOLES"] = "";
    defines["TILE_SIZE"] = cc.intToString(ComputeContext::TileSize);
    int numExclusionTiles = tilesWithExclusions.size();
    defines["NUM_TILES_WITH_EXCLUSIONS"] = cc.intToString(numExclusionTiles);
    int numContexts = cc.getNumContexts();
    int startExclusionIndex = cc.getContextIndex()*numExclusionTiles/numContexts;
    int endExclusionIndex = (cc.getContextIndex()+1)*numExclusionTiles/numContexts;
    defines["FIRST_EXCLUSION_TILE"] = cc.intToString(startExclusionIndex);
    defines["LAST_EXCLUSION_TILE"] = cc.intToString(endExclusionIndex);
    maxExtrapolationOrder = force.getExtrapolationCoefficients().size();
    defines["MAX_EXTRAPOLATION_ORDER"] = cc.intToString(maxExtrapolationOrder);
    stringstream coefficients;
    for (int i = 0; i < maxExtrapolationOrder; i++) {
        if (i > 0)
            coefficients << ",";
        double sum = 0;
        for (int j = i; j < maxExtrapolationOrder; j++)
            sum += force.getExtrapolationCoefficients()[j];
        coefficients << cc.doubleToString(sum);
    }
    defines["EXTRAPOLATION_COEFFICIENTS_SUM"] = coefficients.str();
    if (usePME) {
        int nx, ny, nz;
        force.getPMEParameters(pmeAlpha, nx, ny, nz);
        if (nx == 0 || pmeAlpha == 0) {
            NonbondedForce nb;
            nb.setEwaldErrorTolerance(force.getEwaldErrorTolerance());
            nb.setCutoffDistance(force.getCutoffDistance());
            NonbondedForceImpl::calcPMEParameters(system, nb, pmeAlpha, gridSizeX, gridSizeY, gridSizeZ, false);
            gridSizeX = cc.findLegalFFTDimension(gridSizeX);
            gridSizeY = cc.findLegalFFTDimension(gridSizeY);
            gridSizeZ = cc.findLegalFFTDimension(gridSizeZ);
        } else {
            gridSizeX = cc.findLegalFFTDimension(nx);
            gridSizeY = cc.findLegalFFTDimension(ny);
            gridSizeZ = cc.findLegalFFTDimension(nz);
        }
        defines["EWALD_ALPHA"] = cc.doubleToString(pmeAlpha);
        defines["SQRT_PI"] = cc.doubleToString(sqrt(M_PI));
        defines["USE_EWALD"] = "";
        defines["USE_CUTOFF"] = "";
        defines["USE_PERIODIC"] = "";
        defines["CUTOFF_SQUARED"] = cc.doubleToString(force.getCutoffDistance()*force.getCutoffDistance());
    }
    if (gk != NULL) {
        defines["USE_GK"] = "";
        defines["GK_C"] = cc.doubleToString(2.455);
        double solventDielectric = gk->getSolventDielectric();
        defines["GK_FC"] = cc.doubleToString(1*(1-solventDielectric)/(0+1*solventDielectric));
        defines["GK_FD"] = cc.doubleToString(2*(1-solventDielectric)/(1+2*solventDielectric));
        defines["GK_FQ"] = cc.doubleToString(3*(1-solventDielectric)/(2+3*solventDielectric));
        fixedThreadMemory += 4*elementSize;
        inducedThreadMemory += 13*elementSize;
        if (polarizationType == AmoebaMultipoleForce::Mutual) {
            prevDipolesGk.initialize(cc, 3*numMultipoles*MaxPrevDIISDipoles, elementSize, "prevDipolesGk");
            prevDipolesGkPolar.initialize(cc, 3*numMultipoles*MaxPrevDIISDipoles, elementSize, "prevDipolesGkPolar");
        }
        else if (polarizationType == AmoebaMultipoleForce::Extrapolated) {
            inducedThreadMemory += 12*elementSize;
            int numOrders = force.getExtrapolationCoefficients().size();
            extrapolatedDipoleGk.initialize(cc, 3*numMultipoles*numOrders, elementSize, "extrapolatedDipoleGk");
            extrapolatedDipoleGkPolar.initialize(cc, 3*numMultipoles*numOrders, elementSize, "extrapolatedDipoleGkPolar");
            inducedDipoleFieldGradientGk.initialize(cc, 6*paddedNumAtoms, sizeof(long long), "inducedDipoleFieldGradientGk");
            inducedDipoleFieldGradientGkPolar.initialize(cc, 6*paddedNumAtoms, sizeof(long long), "inducedDipoleFieldGradientGkPolar");
            extrapolatedDipoleFieldGradientGk.initialize(cc, 6*paddedNumAtoms*(numOrders-1), elementSize, "extrapolatedDipoleFieldGradientGk");
            extrapolatedDipoleFieldGradientGkPolar.initialize(cc, 6*paddedNumAtoms*(numOrders-1), elementSize, "extrapolatedDipoleFieldGradientGkPolar");
        }
    }
    NonbondedUtilities& nb = cc.getNonbondedUtilities();
    int maxThreads = max(32, nb.getForceThreadBlockSize());
    fixedFieldThreads = min(maxThreads, cc.computeThreadBlockSize(fixedThreadMemory));
    inducedFieldThreads = min(maxThreads, cc.computeThreadBlockSize(inducedThreadMemory));
    ComputeProgram program = cc.compileProgram(CommonAmoebaKernelSources::multipoles, defines);
    computeMomentsKernel = program->createKernel("computeLabFrameMoments");
    computeMomentsKernel->addArg(cc.getPosq());
    computeMomentsKernel->addArg(multipoleParticles);
    computeMomentsKernel->addArg(localDipoles);
    computeMomentsKernel->addArg(localQuadrupoles);
    computeMomentsKernel->addArg(labDipoles);
    computeMomentsKernel->addArg(labQuadrupoles);
    computeMomentsKernel->addArg(sphericalDipoles);
    computeMomentsKernel->addArg(sphericalQuadrupoles);
    recordInducedDipolesKernel = program->createKernel("recordInducedDipoles");
    recordInducedDipolesKernel->addArg(field);
    recordInducedDipolesKernel->addArg(fieldPolar);
    if (gk != NULL)
        for (int i = 0; i < 3; i++)
            recordInducedDipolesKernel->addArg();
    recordInducedDipolesKernel->addArg(inducedDipole);
    recordInducedDipolesKernel->addArg(inducedDipolePolar);
    recordInducedDipolesKernel->addArg(polarizability);
    mapTorqueKernel = program->createKernel("mapTorqueToForce");
    mapTorqueKernel->addArg(cc.getLongForceBuffer());
    mapTorqueKernel->addArg(torque);
    mapTorqueKernel->addArg(cc.getPosq());
    mapTorqueKernel->addArg(multipoleParticles);
    computePotentialKernel = program->createKernel("computePotentialAtPoints");
    computePotentialKernel->addArg(cc.getPosq());
    computePotentialKernel->addArg(labDipoles);
    computePotentialKernel->addArg(labQuadrupoles);
    computePotentialKernel->addArg(inducedDipole);
    for (int i = 0; i < 8; i++)
        computePotentialKernel->addArg();
    defines["THREAD_BLOCK_SIZE"] = cc.intToString(fixedFieldThreads);
    program = cc.compileProgram(CommonAmoebaKernelSources::multipoleFixedField, defines);
    computeFixedFieldKernel = program->createKernel("computeFixedField");
    computeFixedFieldKernel->addArg(field);
    computeFixedFieldKernel->addArg(fieldPolar);
    computeFixedFieldKernel->addArg(cc.getPosq());
    computeFixedFieldKernel->addArg(covalentFlags);
    computeFixedFieldKernel->addArg(polarizationGroupFlags);
    computeFixedFieldKernel->addArg(nb.getExclusionTiles());
    computeFixedFieldKernel->addArg();
    computeFixedFieldKernel->addArg();
    if (gk != NULL) {
        computeFixedFieldKernel->addArg();
        computeFixedFieldKernel->addArg();
    }
    else if (usePME) {
        computeFixedFieldKernel->addArg(nb.getInteractingTiles());
        computeFixedFieldKernel->addArg(nb.getInteractionCount());
        for (int i = 0; i < 6; i++)
            computeFixedFieldKernel->addArg();
        computeFixedFieldKernel->addArg(nb.getBlockCenters());
        computeFixedFieldKernel->addArg(nb.getInteractingAtoms());
    }
    computeFixedFieldKernel->addArg(labDipoles);
    computeFixedFieldKernel->addArg(labQuadrupoles);
    computeFixedFieldKernel->addArg(dampingAndThole);
    if (polarizationType != AmoebaMultipoleForce::Direct) {
        defines["THREAD_BLOCK_SIZE"] = cc.intToString(inducedFieldThreads);
        defines["MAX_PREV_DIIS_DIPOLES"] = cc.intToString(MaxPrevDIISDipoles);
        program = cc.compileProgram(CommonAmoebaKernelSources::multipoleInducedField, defines);
        computeInducedFieldKernel = program->createKernel("computeInducedField");
        computeInducedFieldKernel->addArg(inducedField);
        computeInducedFieldKernel->addArg(inducedFieldPolar);
        computeInducedFieldKernel->addArg(cc.getPosq());
        computeInducedFieldKernel->addArg(nb.getExclusionTiles());
        computeInducedFieldKernel->addArg(inducedDipole);
        computeInducedFieldKernel->addArg(inducedDipolePolar);
        computeInducedFieldKernel->addArg();
        computeInducedFieldKernel->addArg();
        if (usePME) {
            computeInducedFieldKernel->addArg(nb.getInteractingTiles());
            computeInducedFieldKernel->addArg(nb.getInteractionCount());
            for (int i = 0; i < 6; i++)
                computeInducedFieldKernel->addArg();
            computeInducedFieldKernel->addArg(nb.getBlockCenters());
            computeInducedFieldKernel->addArg(nb.getInteractingAtoms());
        }
        if (gk != NULL) {
            for (int i = 0; i < 5; i++)
                computeInducedFieldKernel->addArg();
        }
        if (polarizationType == AmoebaMultipoleForce::Extrapolated) {
            computeInducedFieldKernel->addArg(inducedDipoleFieldGradient);
            computeInducedFieldKernel->addArg(inducedDipoleFieldGradientPolar);
            if (gk != NULL) {
                computeInducedFieldKernel->addArg(inducedDipoleFieldGradientGk);
                computeInducedFieldKernel->addArg(inducedDipoleFieldGradientGkPolar);
            }
        }
        computeInducedFieldKernel->addArg(dampingAndThole);
        if (polarizationType == AmoebaMultipoleForce::Mutual) {
            updateInducedFieldKernel = program->createKernel("updateInducedFieldByDIIS");
            for (int i = 0; i < 4; i++)
                updateInducedFieldKernel->addArg();
            updateInducedFieldKernel->addArg(diisCoefficients);
            updateInducedFieldKernel->addArg();
            recordDIISDipolesKernel = program->createKernel("recordInducedDipolesForDIIS");
            recordDIISDipolesKernel->addArg(field);
            recordDIISDipolesKernel->addArg(fieldPolar);
            recordDIISDipolesKernel->addArg(polarizability);
            recordDIISDipolesKernel->addArg(inducedDipoleErrors);
            recordDIISDipolesKernel->addArg(prevErrors);
            recordDIISDipolesKernel->addArg(diisMatrix);
            for (int i = 0; i < 9; i++)
                recordDIISDipolesKernel->addArg();
            buildMatrixKernel = program->createKernel("computeDIISMatrix");
            buildMatrixKernel->addArg(prevErrors);
            buildMatrixKernel->addArg();
            buildMatrixKernel->addArg(diisMatrix);
            solveMatrixKernel = program->createKernel("solveDIISMatrix");
            solveMatrixKernel->addArg();
            solveMatrixKernel->addArg(diisMatrix);
            solveMatrixKernel->addArg(diisCoefficients);
        }
        if (polarizationType == AmoebaMultipoleForce::Extrapolated) {
            initExtrapolatedKernel = program->createKernel("initExtrapolatedDipoles");
            initExtrapolatedKernel->addArg(inducedDipole);
            initExtrapolatedKernel->addArg(extrapolatedDipole);
            initExtrapolatedKernel->addArg(inducedDipolePolar);
            initExtrapolatedKernel->addArg(extrapolatedDipolePolar);
            initExtrapolatedKernel->addArg(inducedDipoleFieldGradient);
            initExtrapolatedKernel->addArg(inducedDipoleFieldGradientPolar);
            if (gk != NULL) {
                initExtrapolatedKernel->addArg();
                initExtrapolatedKernel->addArg();
                initExtrapolatedKernel->addArg(extrapolatedDipoleGk);
                initExtrapolatedKernel->addArg(extrapolatedDipoleGkPolar);
                initExtrapolatedKernel->addArg(inducedDipoleFieldGradientGk);
                initExtrapolatedKernel->addArg(inducedDipoleFieldGradientGkPolar);
            }
            iterateExtrapolatedKernel = program->createKernel("iterateExtrapolatedDipoles");
            iterateExtrapolatedKernel->addArg();
            iterateExtrapolatedKernel->addArg(inducedDipole);
            iterateExtrapolatedKernel->addArg(extrapolatedDipole);
            iterateExtrapolatedKernel->addArg(inducedField);
            iterateExtrapolatedKernel->addArg(inducedDipolePolar);
            iterateExtrapolatedKernel->addArg(extrapolatedDipolePolar);
            iterateExtrapolatedKernel->addArg(inducedFieldPolar);
            iterateExtrapolatedKernel->addArg(inducedDipoleFieldGradient);
            iterateExtrapolatedKernel->addArg(inducedDipoleFieldGradientPolar);
            iterateExtrapolatedKernel->addArg(extrapolatedDipoleFieldGradient);
            iterateExtrapolatedKernel->addArg(extrapolatedDipoleFieldGradientPolar);
            if (gk != NULL) {
                iterateExtrapolatedKernel->addArg();
                iterateExtrapolatedKernel->addArg();
                iterateExtrapolatedKernel->addArg(extrapolatedDipoleGk);
                iterateExtrapolatedKernel->addArg(extrapolatedDipoleGkPolar);
                iterateExtrapolatedKernel->addArg(inducedDipoleFieldGradientGk);
                iterateExtrapolatedKernel->addArg(inducedDipoleFieldGradientGkPolar);
                iterateExtrapolatedKernel->addArg();
                iterateExtrapolatedKernel->addArg();
                iterateExtrapolatedKernel->addArg(extrapolatedDipoleFieldGradientGk);
                iterateExtrapolatedKernel->addArg(extrapolatedDipoleFieldGradientGkPolar);
            }
            iterateExtrapolatedKernel->addArg(polarizability);
            computeExtrapolatedKernel = program->createKernel("computeExtrapolatedDipoles");
            computeExtrapolatedKernel->addArg(inducedDipole);
            computeExtrapolatedKernel->addArg(extrapolatedDipole);
            computeExtrapolatedKernel->addArg(inducedDipolePolar);
            computeExtrapolatedKernel->addArg(extrapolatedDipolePolar);
            if (gk != NULL) {
                computeExtrapolatedKernel->addArg();
                computeExtrapolatedKernel->addArg();
                computeExtrapolatedKernel->addArg(extrapolatedDipoleGk);
                computeExtrapolatedKernel->addArg(extrapolatedDipoleGkPolar);
            }
            addExtrapolatedGradientKernel = program->createKernel("addExtrapolatedFieldGradientToForce");
            addExtrapolatedGradientKernel->addArg(cc.getLongForceBuffer());
            addExtrapolatedGradientKernel->addArg(extrapolatedDipole);
            addExtrapolatedGradientKernel->addArg(extrapolatedDipolePolar);
            addExtrapolatedGradientKernel->addArg(extrapolatedDipoleFieldGradient);
            addExtrapolatedGradientKernel->addArg(extrapolatedDipoleFieldGradientPolar);
            if (gk != NULL) {
                addExtrapolatedGradientKernel->addArg(extrapolatedDipoleGk);
                addExtrapolatedGradientKernel->addArg(extrapolatedDipoleGkPolar);
                addExtrapolatedGradientKernel->addArg(extrapolatedDipoleFieldGradientGk);
                addExtrapolatedGradientKernel->addArg(extrapolatedDipoleFieldGradientGkPolar);
            }
        }
    }
    stringstream electrostaticsSource;
    electrostaticsSource << CommonAmoebaKernelSources::sphericalMultipoles;
    if (usePME)
        electrostaticsSource << CommonAmoebaKernelSources::pmeMultipoleElectrostatics;
    else
        electrostaticsSource << CommonAmoebaKernelSources::multipoleElectrostatics;
    electrostaticsThreadMemory = 24*elementSize+3*sizeof(float)+3*sizeof(int)/(double) cc.TileSize;
    electrostaticsThreads = min(maxThreads, cc.computeThreadBlockSize(electrostaticsThreadMemory));
    defines["THREAD_BLOCK_SIZE"] = cc.intToString(electrostaticsThreads);
    program = cc.compileProgram(electrostaticsSource.str(), defines);
    electrostaticsKernel = program->createKernel("computeElectrostatics");
    electrostaticsKernel->addArg(cc.getLongForceBuffer());
    electrostaticsKernel->addArg(torque);
    electrostaticsKernel->addArg(cc.getEnergyBuffer());
    electrostaticsKernel->addArg(cc.getPosq());
    electrostaticsKernel->addArg(covalentFlags);
    electrostaticsKernel->addArg(polarizationGroupFlags);
    electrostaticsKernel->addArg(nb.getExclusionTiles());
    electrostaticsKernel->addArg();
    electrostaticsKernel->addArg();
    if (usePME) {
        electrostaticsKernel->addArg(nb.getInteractingTiles());
        electrostaticsKernel->addArg(nb.getInteractionCount());
        for (int i = 0; i < 6; i++)
            electrostaticsKernel->addArg();
        electrostaticsKernel->addArg(nb.getBlockCenters());
        electrostaticsKernel->addArg(nb.getInteractingAtoms());
    }
    electrostaticsKernel->addArg(sphericalDipoles);
    electrostaticsKernel->addArg(sphericalQuadrupoles);
    electrostaticsKernel->addArg(inducedDipole);
    electrostaticsKernel->addArg(inducedDipolePolar);
    electrostaticsKernel->addArg(dampingAndThole);

    // Set up PME.
    
    if (usePME) {
        // Create required data structures.

        int elementSize = (cc.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
        pmeGrid1.initialize(cc, gridSizeX*gridSizeY*gridSizeZ, 2*elementSize, "pmeGrid1");
        pmeGrid2.initialize(cc, gridSizeX*gridSizeY*gridSizeZ, 2*elementSize, "pmeGrid2");
        if (useFixedPointChargeSpreading()) {
            pmeGridLong.initialize(cc, 2*gridSizeX*gridSizeY*gridSizeZ, sizeof(long long), "pmeGridLong");
            cc.addAutoclearBuffer(pmeGridLong);
        }
        else
            cc.addAutoclearBuffer(pmeGrid1);
        pmeBsplineModuliX.initialize(cc, gridSizeX, elementSize, "pmeBsplineModuliX");
        pmeBsplineModuliY.initialize(cc, gridSizeY, elementSize, "pmeBsplineModuliY");
        pmeBsplineModuliZ.initialize(cc, gridSizeZ, elementSize, "pmeBsplineModuliZ");
        pmePhi.initialize(cc, 20*numMultipoles, elementSize, "pmePhi");
        pmePhid.initialize(cc, 10*numMultipoles, elementSize, "pmePhid");
        pmePhip.initialize(cc, 10*numMultipoles, elementSize, "pmePhip");
        pmePhidp.initialize(cc, 20*numMultipoles, elementSize, "pmePhidp");
        pmeCphi.initialize(cc, 10*numMultipoles, elementSize, "pmeCphi");

        // Create the PME kernels.

        map<string, string> pmeDefines;
        pmeDefines["EWALD_ALPHA"] = cc.doubleToString(pmeAlpha);
        pmeDefines["PME_ORDER"] = cc.intToString(PmeOrder);
        pmeDefines["NUM_ATOMS"] = cc.intToString(numMultipoles);
        pmeDefines["PADDED_NUM_ATOMS"] = cc.intToString(cc.getPaddedNumAtoms());
        pmeDefines["EPSILON_FACTOR"] = cc.doubleToString(ONE_4PI_EPS0);
        pmeDefines["GRID_SIZE_X"] = cc.intToString(gridSizeX);
        pmeDefines["GRID_SIZE_Y"] = cc.intToString(gridSizeY);
        pmeDefines["GRID_SIZE_Z"] = cc.intToString(gridSizeZ);
        pmeDefines["M_PI"] = cc.doubleToString(M_PI);
        pmeDefines["SQRT_PI"] = cc.doubleToString(sqrt(M_PI));
        if (polarizationType == AmoebaMultipoleForce::Direct)
            pmeDefines["DIRECT_POLARIZATION"] = "";
        else if (polarizationType == AmoebaMultipoleForce::Mutual)
            pmeDefines["MUTUAL_POLARIZATION"] = "";
        else if (polarizationType == AmoebaMultipoleForce::Extrapolated)
            pmeDefines["EXTRAPOLATED_POLARIZATION"] = "";
        if (useFixedPointChargeSpreading())
            pmeDefines["USE_FIXED_POINT_CHARGE_SPREADING"] = "";
        program = cc.compileProgram(CommonAmoebaKernelSources::multipolePme, pmeDefines);
        pmeTransformMultipolesKernel = program->createKernel("transformMultipolesToFractionalCoordinates");
        pmeTransformMultipolesKernel->addArg(labDipoles);
        pmeTransformMultipolesKernel->addArg(labQuadrupoles);
        pmeTransformMultipolesKernel->addArg(fracDipoles);
        pmeTransformMultipolesKernel->addArg(fracQuadrupoles);
        for (int i = 0; i < 3; i++)
            pmeTransformMultipolesKernel->addArg();
        pmeTransformPotentialKernel = program->createKernel("transformPotentialToCartesianCoordinates");
        pmeTransformPotentialKernel->addArg();
        pmeTransformPotentialKernel->addArg(pmeCphi);
        for (int i = 0; i < 3; i++)
            pmeTransformPotentialKernel->addArg();
        pmeSpreadFixedMultipolesKernel = program->createKernel("gridSpreadFixedMultipoles");
        pmeSpreadFixedMultipolesKernel->addArg(cc.getPosq());
        pmeSpreadFixedMultipolesKernel->addArg(fracDipoles);
        pmeSpreadFixedMultipolesKernel->addArg(fracQuadrupoles);
        if (useFixedPointChargeSpreading())
            pmeSpreadFixedMultipolesKernel->addArg(pmeGridLong);
        else
            pmeSpreadFixedMultipolesKernel->addArg(pmeGrid1);
        for (int i = 0; i < 6; i++)
            pmeSpreadFixedMultipolesKernel->addArg();
        pmeSpreadInducedDipolesKernel = program->createKernel("gridSpreadInducedDipoles");
        pmeSpreadInducedDipolesKernel->addArg(cc.getPosq());
        pmeSpreadInducedDipolesKernel->addArg(inducedDipole);
        pmeSpreadInducedDipolesKernel->addArg(inducedDipolePolar);
        if (useFixedPointChargeSpreading())
            pmeSpreadInducedDipolesKernel->addArg(pmeGridLong);
        else
            pmeSpreadInducedDipolesKernel->addArg(pmeGrid1);
        for (int i = 0; i < 6; i++)
            pmeSpreadInducedDipolesKernel->addArg();
        if (useFixedPointChargeSpreading()) {
            pmeFinishSpreadChargeKernel = program->createKernel("finishSpreadCharge");
            pmeFinishSpreadChargeKernel->addArg(pmeGridLong);
            pmeFinishSpreadChargeKernel->addArg(pmeGrid1);
        }
        pmeConvolutionKernel = program->createKernel("reciprocalConvolution");
        pmeConvolutionKernel->addArg(pmeGrid2);
        pmeConvolutionKernel->addArg(pmeBsplineModuliX);
        pmeConvolutionKernel->addArg(pmeBsplineModuliY);
        pmeConvolutionKernel->addArg(pmeBsplineModuliZ);
        for (int i = 0; i < 4; i++)
            pmeConvolutionKernel->addArg();
        pmeFixedPotentialKernel = program->createKernel("computeFixedPotentialFromGrid");
        pmeFixedPotentialKernel->addArg(pmeGrid1);
        pmeFixedPotentialKernel->addArg(pmePhi);
        pmeFixedPotentialKernel->addArg(field);
        pmeFixedPotentialKernel->addArg(fieldPolar);
        pmeFixedPotentialKernel->addArg(cc.getPosq());
        pmeFixedPotentialKernel->addArg(labDipoles);
        for (int i = 0; i < 6; i++)
            pmeFixedPotentialKernel->addArg();
        pmeInducedPotentialKernel = program->createKernel("computeInducedPotentialFromGrid");
        pmeInducedPotentialKernel->addArg(pmeGrid1);
        pmeInducedPotentialKernel->addArg(pmePhid);
        pmeInducedPotentialKernel->addArg(pmePhip);
        pmeInducedPotentialKernel->addArg(pmePhidp);
        pmeInducedPotentialKernel->addArg(cc.getPosq());
        for (int i = 0; i < 6; i++)
            pmeInducedPotentialKernel->addArg();
        pmeFixedForceKernel = program->createKernel("computeFixedMultipoleForceAndEnergy");
        pmeFixedForceKernel->addArg(cc.getPosq());
        pmeFixedForceKernel->addArg(cc.getLongForceBuffer());
        pmeFixedForceKernel->addArg(torque);
        pmeFixedForceKernel->addArg(cc.getEnergyBuffer());
        pmeFixedForceKernel->addArg(labDipoles);
        pmeFixedForceKernel->addArg(labQuadrupoles);
        pmeFixedForceKernel->addArg(fracDipoles);
        pmeFixedForceKernel->addArg(fracQuadrupoles);
        pmeFixedForceKernel->addArg(pmePhi);
        pmeFixedForceKernel->addArg(pmeCphi);
        for (int i = 0; i < 3; i++)
            pmeFixedForceKernel->addArg();
        pmeInducedForceKernel = program->createKernel("computeInducedDipoleForceAndEnergy");
        pmeInducedForceKernel->addArg(cc.getPosq());
        pmeInducedForceKernel->addArg(cc.getLongForceBuffer());
        pmeInducedForceKernel->addArg(torque);
        pmeInducedForceKernel->addArg(cc.getEnergyBuffer());
        pmeInducedForceKernel->addArg(labDipoles);
        pmeInducedForceKernel->addArg(labQuadrupoles);
        pmeInducedForceKernel->addArg(fracDipoles);
        pmeInducedForceKernel->addArg(fracQuadrupoles);
        pmeInducedForceKernel->addArg(inducedDipole);
        pmeInducedForceKernel->addArg(inducedDipolePolar);
        pmeInducedForceKernel->addArg(pmePhi);
        pmeInducedForceKernel->addArg(pmePhid);
        pmeInducedForceKernel->addArg(pmePhip);
        pmeInducedForceKernel->addArg(pmePhidp);
        pmeInducedForceKernel->addArg(pmeCphi);
        for (int i = 0; i < 3; i++)
            pmeInducedForceKernel->addArg();
        if (polarizationType != AmoebaMultipoleForce::Direct) {
            pmeRecordInducedFieldDipolesKernel = program->createKernel("recordInducedFieldDipoles");
            pmeRecordInducedFieldDipolesKernel->addArg(pmePhid);
            pmeRecordInducedFieldDipolesKernel->addArg(pmePhip);
            pmeRecordInducedFieldDipolesKernel->addArg(inducedField);
            pmeRecordInducedFieldDipolesKernel->addArg(inducedFieldPolar);
            pmeRecordInducedFieldDipolesKernel->addArg(inducedDipole);
            pmeRecordInducedFieldDipolesKernel->addArg(inducedDipolePolar);
            for (int i = 0; i < 3; i++)
                pmeRecordInducedFieldDipolesKernel->addArg();
            if (polarizationType == AmoebaMultipoleForce::Extrapolated) {
                pmeRecordInducedFieldDipolesKernel->addArg(inducedDipoleFieldGradient);
                pmeRecordInducedFieldDipolesKernel->addArg(inducedDipoleFieldGradientPolar);
            }
        }

        // Initialize the B-spline moduli.

        double data[PmeOrder];
        double x = 0.0;
        data[0] = 1.0 - x;
        data[1] = x;
        for (int i = 2; i < PmeOrder; i++) {
            double denom = 1.0/i;
            data[i] = x*data[i-1]*denom;
            for (int j = 1; j < i; j++)
                data[i-j] = ((x+j)*data[i-j-1] + ((i-j+1)-x)*data[i-j])*denom;
            data[0] = (1.0-x)*data[0]*denom;
        }
        int maxSize = max(max(gridSizeX, gridSizeY), gridSizeZ);
        vector<double> bsplines_data(maxSize+1, 0.0);
        for (int i = 2; i <= PmeOrder+1; i++)
            bsplines_data[i] = data[i-2];
        for (int dim = 0; dim < 3; dim++) {
            int ndata = (dim == 0 ? gridSizeX : dim == 1 ? gridSizeY : gridSizeZ);
            vector<double> moduli(ndata);

            // get the modulus of the discrete Fourier transform

            double factor = 2.0*M_PI/ndata;
            for (int i = 0; i < ndata; i++) {
                double sc = 0.0;
                double ss = 0.0;
                for (int j = 1; j <= ndata; j++) {
                    double arg = factor*i*(j-1);
                    sc += bsplines_data[j]*cos(arg);
                    ss += bsplines_data[j]*sin(arg);
                }
                moduli[i] = sc*sc+ss*ss;
            }

            // Fix for exponential Euler spline interpolation failure.

            double eps = 1.0e-7;
            if (moduli[0] < eps)
                moduli[0] = 0.9*moduli[1];
            for (int i = 1; i < ndata-1; i++)
                if (moduli[i] < eps)
                    moduli[i] = 0.9*(moduli[i-1]+moduli[i+1]);
            if (moduli[ndata-1] < eps)
                moduli[ndata-1] = 0.9*moduli[ndata-2];

            // Compute and apply the optimal zeta coefficient.

            int jcut = 50;
            for (int i = 1; i <= ndata; i++) {
                int k = i - 1;
                if (i > ndata/2)
                    k = k - ndata;
                double zeta;
                if (k == 0)
                    zeta = 1.0;
                else {
                    double sum1 = 1.0;
                    double sum2 = 1.0;
                    factor = M_PI*k/ndata;
                    for (int j = 1; j <= jcut; j++) {
                        double arg = factor/(factor+M_PI*j);
                        sum1 += pow(arg, PmeOrder);
                        sum2 += pow(arg, 2*PmeOrder);
                    }
                    for (int j = 1; j <= jcut; j++) {
                        double arg = factor/(factor-M_PI*j);
                        sum1 += pow(arg, PmeOrder);
                        sum2 += pow(arg, 2*PmeOrder);
                    }
                    zeta = sum2/sum1;
                }
                moduli[i-1] = moduli[i-1]*zeta*zeta;
            }
            if (cc.getUseDoublePrecision()) {
                if (dim == 0)
                    pmeBsplineModuliX.upload(moduli);
                else if (dim == 1)
                    pmeBsplineModuliY.upload(moduli);
                else
                    pmeBsplineModuliZ.upload(moduli);
            }
            else {
                vector<float> modulif(ndata);
                for (int i = 0; i < ndata; i++)
                    modulif[i] = (float) moduli[i];
                if (dim == 0)
                    pmeBsplineModuliX.upload(modulif);
                else if (dim == 1)
                    pmeBsplineModuliY.upload(modulif);
                else
                    pmeBsplineModuliZ.upload(modulif);
            }
        }
    }

    // Add an interaction to the default nonbonded kernel.  This doesn't actually do any calculations.  It's
    // just so that NonbondedUtilities will build the exclusion flags and maintain the neighbor list.
    
    cc.getNonbondedUtilities().addInteraction(usePME, usePME, true, force.getCutoffDistance(), exclusions, "", force.getForceGroup());
    cc.getNonbondedUtilities().setUsePadding(false);
    cc.addForce(new ForceInfo(force));
}

void CommonCalcAmoebaMultipoleForceKernel::initializeScaleFactors() {
    hasInitializedScaleFactors = true;
    NonbondedUtilities& nb = cc.getNonbondedUtilities();
    
    // Figure out the covalent flag values to use for each atom pair.

    vector<mm_int2> exclusionTiles;
    nb.getExclusionTiles().download(exclusionTiles);
    map<pair<int, int>, int> exclusionTileMap;
    for (int i = 0; i < (int) exclusionTiles.size(); i++) {
        mm_int2 tile = exclusionTiles[i];
        exclusionTileMap[make_pair(tile.x, tile.y)] = i;
    }
    covalentFlags.resize(nb.getExclusions().getSize());
    vector<mm_int2> covalentFlagsVec(nb.getExclusions().getSize(), mm_int2(0, 0));
    for (mm_int4 values : covalentFlagValues) {
        int atom1 = values.x;
        int atom2 = values.y;
        int value = values.z;
        int x = atom1/ComputeContext::TileSize;
        int offset1 = atom1-x*ComputeContext::TileSize;
        int y = atom2/ComputeContext::TileSize;
        int offset2 = atom2-y*ComputeContext::TileSize;
        int f1 = (value == 0 || value == 1 ? 1 : 0);
        int f2 = (value == 0 || value == 2 ? 1 : 0);
        if (x == y) {
            int index = exclusionTileMap[make_pair(x, y)]*ComputeContext::TileSize;
            covalentFlagsVec[index+offset1].x |= f1<<offset2;
            covalentFlagsVec[index+offset1].y |= f2<<offset2;
            covalentFlagsVec[index+offset2].x |= f1<<offset1;
            covalentFlagsVec[index+offset2].y |= f2<<offset1;
        }
        else if (x > y) {
            int index = exclusionTileMap[make_pair(x, y)]*ComputeContext::TileSize;
            covalentFlagsVec[index+offset1].x |= f1<<offset2;
            covalentFlagsVec[index+offset1].y |= f2<<offset2;
        }
        else {
            int index = exclusionTileMap[make_pair(y, x)]*ComputeContext::TileSize;
            covalentFlagsVec[index+offset2].x |= f1<<offset1;
            covalentFlagsVec[index+offset2].y |= f2<<offset1;
        }
    }
    covalentFlags.upload(covalentFlagsVec);
    
    // Do the same for the polarization flags.
    
    polarizationGroupFlags.resize(nb.getExclusions().getSize());
    vector<unsigned int> polarizationGroupFlagsVec(nb.getExclusions().getSize(), 0);
    for (mm_int2 values : polarizationFlagValues) {
        int atom1 = values.x;
        int atom2 = values.y;
        int x = atom1/ComputeContext::TileSize;
        int offset1 = atom1-x*ComputeContext::TileSize;
        int y = atom2/ComputeContext::TileSize;
        int offset2 = atom2-y*ComputeContext::TileSize;
        if (x == y) {
            int index = exclusionTileMap[make_pair(x, y)]*ComputeContext::TileSize;
            polarizationGroupFlagsVec[index+offset1] |= 1<<offset2;
            polarizationGroupFlagsVec[index+offset2] |= 1<<offset1;
        }
        else if (x > y) {
            int index = exclusionTileMap[make_pair(x, y)]*ComputeContext::TileSize;
            polarizationGroupFlagsVec[index+offset1] |= 1<<offset2;
        }
        else {
            int index = exclusionTileMap[make_pair(y, x)]*ComputeContext::TileSize;
            polarizationGroupFlagsVec[index+offset2] |= 1<<offset1;
        }
    }
    polarizationGroupFlags.upload(polarizationGroupFlagsVec);
}

double CommonCalcAmoebaMultipoleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    ContextSelector selector(cc);
    if (!hasInitializedScaleFactors) {
        initializeScaleFactors();
        for (auto impl : context.getForceImpls()) {
            AmoebaGeneralizedKirkwoodForceImpl* gkImpl = dynamic_cast<AmoebaGeneralizedKirkwoodForceImpl*>(impl);
            if (gkImpl != NULL) {
                gkKernel = dynamic_cast<CommonCalcAmoebaGeneralizedKirkwoodForceKernel*>(&gkImpl->getKernel().getImpl());
                recordInducedDipolesKernel->setArg(2, gkKernel->getField());
                recordInducedDipolesKernel->setArg(3, gkKernel->getInducedDipoles());
                recordInducedDipolesKernel->setArg(4, gkKernel->getInducedDipolesPolar());
                computeFixedFieldKernel->setArg(8, gkKernel->getBornRadii());
                computeFixedFieldKernel->setArg(9, gkKernel->getField());
                if (polarizationType != AmoebaMultipoleForce::Direct) {
                    computeInducedFieldKernel->setArg(8, gkKernel->getInducedField());
                    computeInducedFieldKernel->setArg(9,  gkKernel->getInducedFieldPolar());
                    computeInducedFieldKernel->setArg(10, gkKernel->getInducedDipoles());
                    computeInducedFieldKernel->setArg(11, gkKernel->getInducedDipolesPolar());
                    computeInducedFieldKernel->setArg(12, gkKernel->getBornRadii());
                }
                if (polarizationType == AmoebaMultipoleForce::Extrapolated) {
                    initExtrapolatedKernel->setArg(6, gkKernel->getInducedDipoles());
                    initExtrapolatedKernel->setArg(7, gkKernel->getInducedDipolesPolar());
                    iterateExtrapolatedKernel->setArg(11, gkKernel->getInducedDipoles());
                    iterateExtrapolatedKernel->setArg(12, gkKernel->getInducedDipolesPolar());
                    iterateExtrapolatedKernel->setArg(17, gkKernel->getInducedField());
                    iterateExtrapolatedKernel->setArg(18, gkKernel->getInducedFieldPolar());
                    computeExtrapolatedKernel->setArg(4, gkKernel->getInducedDipoles());
                    computeExtrapolatedKernel->setArg(5, gkKernel->getInducedDipolesPolar());
                }
                break;
            }
        }
    }
    NonbondedUtilities& nb = cc.getNonbondedUtilities();
    
    // Compute the lab frame moments.

    computeMomentsKernel->execute(cc.getNumAtoms());
    int startTileIndex = nb.getStartTileIndex();
    int numTileIndices = nb.getNumTiles();
    int numForceThreadBlocks = nb.getNumForceThreadBlocks();
    electrostaticsKernel->setArg(7, startTileIndex);
    electrostaticsKernel->setArg(8, numTileIndices);
    computeFixedFieldKernel->setArg(6, startTileIndex);
    computeFixedFieldKernel->setArg(7, numTileIndices);
    if (!pmeGrid1.isInitialized()) {
        // Compute induced dipoles.
        
        if (gkKernel != NULL)
            gkKernel->computeBornRadii(torque, labDipoles, labQuadrupoles, inducedDipole, inducedDipolePolar, dampingAndThole, covalentFlags, polarizationGroupFlags);
        computeFixedFieldKernel->execute(numForceThreadBlocks*fixedFieldThreads, fixedFieldThreads);
        recordInducedDipolesKernel->execute(cc.getNumAtoms());
        
        // Iterate until the dipoles converge.
        
        if (polarizationType == AmoebaMultipoleForce::Extrapolated)
            computeExtrapolatedDipoles();
        for (int i = 0; i < maxInducedIterations; i++) {
            computeInducedField();
            bool converged = iterateDipolesByDIIS(i);
            if (converged)
                break;
        }
        
        // Compute electrostatic force.
        
        electrostaticsKernel->execute(numForceThreadBlocks*electrostaticsThreads, electrostaticsThreads);
        if (gkKernel != NULL)
            gkKernel->finishComputation();
    }
    else {
        // Compute reciprocal box vectors.
        
        Vec3 a, b, c;
        cc.getPeriodicBoxVectors(a, b, c);
        double determinant = a[0]*b[1]*c[2];
        double scale = 1.0/determinant;
        mm_double4 recipBoxVectors[3];
        recipBoxVectors[0] = mm_double4(b[1]*c[2]*scale, 0, 0, 0);
        recipBoxVectors[1] = mm_double4(-b[0]*c[2]*scale, a[0]*c[2]*scale, 0, 0);
        recipBoxVectors[2] = mm_double4((b[0]*c[1]-b[1]*c[0])*scale, -a[0]*c[1]*scale, a[0]*b[1]*scale, 0);
        if (cc.getUseDoublePrecision()) {
            mm_double4 boxVectors[] = {mm_double4(a[0], a[1], a[2], 0), mm_double4(b[0], b[1], b[2], 0), mm_double4(c[0], c[1], c[2], 0)};
            pmeConvolutionKernel->setArg(4, mm_double4(a[0], b[1], c[2], 0));
            for (int i = 0; i < 3; i++) {
                pmeTransformMultipolesKernel->setArg(4+i, recipBoxVectors[i]);
                pmeTransformPotentialKernel->setArg(2+i, recipBoxVectors[i]);
                pmeSpreadFixedMultipolesKernel->setArg(4+i, boxVectors[i]);
                pmeSpreadFixedMultipolesKernel->setArg(7+i, recipBoxVectors[i]);
                pmeSpreadInducedDipolesKernel->setArg(4+i, boxVectors[i]);
                pmeSpreadInducedDipolesKernel->setArg(7+i, recipBoxVectors[i]);
                pmeConvolutionKernel->setArg(5+i, recipBoxVectors[i]);
                pmeFixedPotentialKernel->setArg(6+i, boxVectors[i]);
                pmeFixedPotentialKernel->setArg(9+i, recipBoxVectors[i]);
                pmeInducedPotentialKernel->setArg(5+i, boxVectors[i]);
                pmeInducedPotentialKernel->setArg(8+i, recipBoxVectors[i]);
                pmeFixedForceKernel->setArg(10+i, recipBoxVectors[i]);
                pmeInducedForceKernel->setArg(15+i, recipBoxVectors[i]);
                if (polarizationType != AmoebaMultipoleForce::Direct)
                    pmeRecordInducedFieldDipolesKernel->setArg(6+i, recipBoxVectors[i]);
            }
        }
        else {
            mm_float4 recipBoxVectorsFloat[3];
            recipBoxVectorsFloat[0] = mm_float4((float) recipBoxVectors[0].x, 0, 0, 0);
            recipBoxVectorsFloat[1] = mm_float4((float) recipBoxVectors[1].x, (float) recipBoxVectors[1].y, 0, 0);
            recipBoxVectorsFloat[2] = mm_float4((float) recipBoxVectors[2].x, (float) recipBoxVectors[2].y, (float) recipBoxVectors[2].z, 0);
            mm_float4 boxVectors[] = {mm_float4(a[0], a[1], a[2], 0), mm_float4(b[0], b[1], b[2], 0), mm_float4(c[0], c[1], c[2], 0)};
            pmeConvolutionKernel->setArg(4, mm_float4(a[0], b[1], c[2], 0));
            for (int i = 0; i < 3; i++) {
                pmeTransformMultipolesKernel->setArg(4+i, recipBoxVectorsFloat[i]);
                pmeTransformPotentialKernel->setArg(2+i, recipBoxVectorsFloat[i]);
                pmeSpreadFixedMultipolesKernel->setArg(4+i, boxVectors[i]);
                pmeSpreadFixedMultipolesKernel->setArg(7+i, recipBoxVectorsFloat[i]);
                pmeSpreadInducedDipolesKernel->setArg(4+i, boxVectors[i]);
                pmeSpreadInducedDipolesKernel->setArg(7+i, recipBoxVectorsFloat[i]);
                pmeConvolutionKernel->setArg(5+i, recipBoxVectorsFloat[i]);
                pmeFixedPotentialKernel->setArg(6+i, boxVectors[i]);
                pmeFixedPotentialKernel->setArg(9+i, recipBoxVectorsFloat[i]);
                pmeInducedPotentialKernel->setArg(5+i, boxVectors[i]);
                pmeInducedPotentialKernel->setArg(8+i, recipBoxVectorsFloat[i]);
                pmeFixedForceKernel->setArg(10+i, recipBoxVectorsFloat[i]);
                pmeInducedForceKernel->setArg(15+i, recipBoxVectorsFloat[i]);
                if (polarizationType != AmoebaMultipoleForce::Direct)
                    pmeRecordInducedFieldDipolesKernel->setArg(6+i, recipBoxVectorsFloat[i]);
            }
        }

        // Reciprocal space calculation.
        
        unsigned int maxTiles = nb.getInteractingTiles().getSize();
        pmeTransformMultipolesKernel->execute(cc.getNumAtoms());
        pmeSpreadFixedMultipolesKernel->execute(cc.getNumAtoms());
        if (useFixedPointChargeSpreading())
            pmeFinishSpreadChargeKernel->execute(pmeGrid1.getSize());
        computeFFT(true);
        pmeConvolutionKernel->execute(gridSizeX*gridSizeY*gridSizeZ, 256);
        computeFFT(false);
        pmeFixedPotentialKernel->execute(cc.getNumAtoms());
        pmeTransformPotentialKernel->setArg(0, pmePhi);
        pmeTransformPotentialKernel->execute(cc.getNumAtoms());
        pmeFixedForceKernel->execute(cc.getNumAtoms());

        // Direct space calculation.

        setPeriodicBoxArgs(cc, computeFixedFieldKernel, 10);
        computeFixedFieldKernel->setArg(15, maxTiles);
        computeFixedFieldKernel->execute(numForceThreadBlocks*fixedFieldThreads, fixedFieldThreads);
        recordInducedDipolesKernel->execute(cc.getNumAtoms());

        // Reciprocal space calculation for the induced dipoles.

        if (useFixedPointChargeSpreading())
            cc.clearBuffer(pmeGridLong);
        else
            cc.clearBuffer(pmeGrid1);
        pmeSpreadInducedDipolesKernel->execute(cc.getNumAtoms());
        if (useFixedPointChargeSpreading())
            pmeFinishSpreadChargeKernel->execute(pmeGrid1.getSize());
        computeFFT(true);
        pmeConvolutionKernel->execute(gridSizeX*gridSizeY*gridSizeZ, 256);
        computeFFT(false);
        pmeInducedPotentialKernel->execute(cc.getNumAtoms());
        
        // Iterate until the dipoles converge.
        
        if (polarizationType == AmoebaMultipoleForce::Extrapolated)
            computeExtrapolatedDipoles();
        for (int i = 0; i < maxInducedIterations; i++) {
            computeInducedField();
            bool converged = iterateDipolesByDIIS(i);
            if (converged)
                break;
        }
        
        // Compute electrostatic force.
        
        setPeriodicBoxArgs(cc, electrostaticsKernel, 11);
        electrostaticsKernel->setArg(16, maxTiles);
        electrostaticsKernel->execute(numForceThreadBlocks*electrostaticsThreads, electrostaticsThreads);
        pmeTransformPotentialKernel->setArg(0, pmePhidp);
        pmeTransformPotentialKernel->execute(cc.getNumAtoms());
        pmeInducedForceKernel->execute(cc.getNumAtoms());
    }
    
    // If using extrapolated polarization, add in force contributions from µ(m) T µ(n).
    
    if (polarizationType == AmoebaMultipoleForce::Extrapolated)
        addExtrapolatedGradientKernel->execute(numMultipoles);

    // Map torques to force.

    mapTorqueKernel->execute(cc.getNumAtoms());
    
    // Record the current atom positions so we can tell later if they have changed.
    
    cc.getPosq().copyTo(lastPositions);
    multipolesAreValid = true;
    return 0.0;
}

void CommonCalcAmoebaMultipoleForceKernel::computeInducedField() {
    NonbondedUtilities& nb = cc.getNonbondedUtilities();
    int startTileIndex = nb.getStartTileIndex();
    int numTileIndices = nb.getNumTiles();
    int numForceThreadBlocks = nb.getNumForceThreadBlocks();
    computeInducedFieldKernel->setArg(6, startTileIndex);
    computeInducedFieldKernel->setArg(7, numTileIndices);
    if (usePME) {
        setPeriodicBoxArgs(cc, computeInducedFieldKernel, 10);
        computeInducedFieldKernel->setArg(15, (int) nb.getInteractingTiles().getSize());
    }
    cc.clearBuffer(inducedField);
    cc.clearBuffer(inducedFieldPolar);
    if (polarizationType == AmoebaMultipoleForce::Extrapolated) {
        cc.clearBuffer(inducedDipoleFieldGradient);
        cc.clearBuffer(inducedDipoleFieldGradientPolar);
    }
    if (gkKernel != NULL) {
        cc.clearBuffer(gkKernel->getInducedField());
        cc.clearBuffer(gkKernel->getInducedFieldPolar());
        if (polarizationType == AmoebaMultipoleForce::Extrapolated) {
            cc.clearBuffer(inducedDipoleFieldGradientGk);
            cc.clearBuffer(inducedDipoleFieldGradientGkPolar);
        }
    }
    computeInducedFieldKernel->execute(numForceThreadBlocks*inducedFieldThreads, inducedFieldThreads);
    if (pmeGrid1.isInitialized()) {
        if (useFixedPointChargeSpreading())
            cc.clearBuffer(pmeGridLong);
        else
            cc.clearBuffer(pmeGrid1);
        pmeSpreadInducedDipolesKernel->execute(cc.getNumAtoms());
        if (useFixedPointChargeSpreading())
            pmeFinishSpreadChargeKernel->execute(pmeGrid1.getSize());
        computeFFT(true);
        pmeConvolutionKernel->execute(gridSizeX*gridSizeY*gridSizeZ, 256);
        computeFFT(false);
        pmeInducedPotentialKernel->execute(cc.getNumAtoms());
        if (polarizationType == AmoebaMultipoleForce::Extrapolated) {
            pmeRecordInducedFieldDipolesKernel->execute(cc.getNumAtoms());
        }
        else {
            pmeRecordInducedFieldDipolesKernel->execute(cc.getNumAtoms());
        }
    }
}

bool CommonCalcAmoebaMultipoleForceKernel::iterateDipolesByDIIS(int iteration) {
    void* npt = NULL;

    // Record the dipoles and errors into the lists of previous dipoles.

    recordDIISDipolesKernel->setArg(13, iteration);
    if (gkKernel != NULL) {
        recordDIISDipolesKernel->setArg(6, gkKernel->getField());
        recordDIISDipolesKernel->setArg(7, gkKernel->getInducedField());
        recordDIISDipolesKernel->setArg(8, gkKernel->getInducedFieldPolar());
        recordDIISDipolesKernel->setArg(9, gkKernel->getInducedDipoles());
        recordDIISDipolesKernel->setArg(10, gkKernel->getInducedDipolesPolar());
        recordDIISDipolesKernel->setArg(11, prevDipolesGk);
        recordDIISDipolesKernel->setArg(12, prevDipolesGkPolar);
        recordDIISDipolesKernel->setArg(14, 1);
        recordDIISDipolesKernel->execute(cc.getNumThreadBlocks()*64, 64);
    }
    recordDIISDipolesKernel->setArg(6, npt);
    recordDIISDipolesKernel->setArg(7, inducedField);
    recordDIISDipolesKernel->setArg(8, inducedFieldPolar);
    recordDIISDipolesKernel->setArg(9, inducedDipole);
    recordDIISDipolesKernel->setArg(10, inducedDipolePolar);
    recordDIISDipolesKernel->setArg(11, prevDipoles);
    recordDIISDipolesKernel->setArg(12, prevDipolesPolar);
    recordDIISDipolesKernel->setArg(14, 0);
    recordDIISDipolesKernel->execute(cc.getNumThreadBlocks()*64, 64);
    mm_float2* errors = (mm_float2*) cc.getPinnedBuffer();
    inducedDipoleErrors.download(errors, false);
    syncEvent->enqueue();
    
    // Build the DIIS matrix.
    
    int numPrev = (iteration+1 < MaxPrevDIISDipoles ? iteration+1 : MaxPrevDIISDipoles);
    int threadBlocks = min(numPrev, cc.getNumThreadBlocks());
    int blockSize = min(512, buildMatrixKernel->getMaxBlockSize());
    buildMatrixKernel->setArg(1, iteration);
    buildMatrixKernel->execute(threadBlocks*blockSize, blockSize);
    
    // Solve the matrix.

    solveMatrixKernel->setArg(0, iteration);
    solveMatrixKernel->execute(32, 32);
    
    // Determine whether the iteration has converged.
    
    syncEvent->wait();
    double total1 = 0.0, total2 = 0.0;
    for (int j = 0; j < inducedDipoleErrors.getSize(); j++) {
        total1 += errors[j].x;
        total2 += errors[j].y;
    }
    if (48.033324*sqrt(max(total1, total2)/cc.getNumAtoms()) < inducedEpsilon)
        return true;
    
    // Compute the dipoles.
    
    updateInducedFieldKernel->setArg(0, inducedDipole);
    updateInducedFieldKernel->setArg(1, inducedDipolePolar);
    updateInducedFieldKernel->setArg(2, prevDipoles);
    updateInducedFieldKernel->setArg(3, prevDipolesPolar);
    updateInducedFieldKernel->setArg(5, numPrev);
    updateInducedFieldKernel->execute(3*cc.getNumAtoms(), 256);
    if (gkKernel != NULL) {
        updateInducedFieldKernel->setArg(0, gkKernel->getInducedDipoles());
        updateInducedFieldKernel->setArg(1, gkKernel->getInducedDipolesPolar());
        updateInducedFieldKernel->setArg(2, prevDipolesGk);
        updateInducedFieldKernel->setArg(3, prevDipolesGkPolar);
        updateInducedFieldKernel->execute(3*cc.getNumAtoms(), 256);
    }
    return false;
}

void CommonCalcAmoebaMultipoleForceKernel::computeExtrapolatedDipoles() {
    // Start by storing the direct dipoles as PT0

    initExtrapolatedKernel->execute(extrapolatedDipole.getSize());

    // Recursively apply alpha.Tau to the µ_(n) components to generate µ_(n+1), and store the result

    for (int order = 1; order < maxExtrapolationOrder; ++order) {
        computeInducedField();
        iterateExtrapolatedKernel->setArg(0, order);
        iterateExtrapolatedKernel->execute(extrapolatedDipole.getSize());
    }
    
    // Take a linear combination of the µ_(n) components to form the total dipole

    computeExtrapolatedKernel->execute(extrapolatedDipole.getSize());
    computeInducedField();
}

void CommonCalcAmoebaMultipoleForceKernel::ensureMultipolesValid(ContextImpl& context) {
    if (multipolesAreValid) {
        int numParticles = cc.getNumAtoms();
        if (cc.getUseDoublePrecision()) {
            vector<mm_double4> pos1, pos2;
            cc.getPosq().download(pos1);
            lastPositions.download(pos2);
            for (int i = 0; i < numParticles; i++)
                if (pos1[i].x != pos2[i].x || pos1[i].y != pos2[i].y || pos1[i].z != pos2[i].z) {
                    multipolesAreValid = false;
                    break;
                }
        }
        else {
            vector<mm_float4> pos1, pos2;
            cc.getPosq().download(pos1);
            lastPositions.download(pos2);
            for (int i = 0; i < numParticles; i++)
                if (pos1[i].x != pos2[i].x || pos1[i].y != pos2[i].y || pos1[i].z != pos2[i].z) {
                    multipolesAreValid = false;
                    break;
                }
        }
    }
    if (!multipolesAreValid)
        context.calcForcesAndEnergy(false, false, context.getIntegrator().getIntegrationForceGroups());
}

void CommonCalcAmoebaMultipoleForceKernel::getLabFramePermanentDipoles(ContextImpl& context, vector<Vec3>& dipoles) {
    ContextSelector selector(cc);
    ensureMultipolesValid(context);
    int numParticles = cc.getNumAtoms();
    dipoles.resize(numParticles);
    const vector<int>& order = cc.getAtomIndex();
    if (cc.getUseDoublePrecision()) {
        vector<double> labDipoleVec;
        labDipoles.download(labDipoleVec);
        for (int i = 0; i < numParticles; i++)
            dipoles[order[i]] = Vec3(labDipoleVec[3*i], labDipoleVec[3*i+1], labDipoleVec[3*i+2]);
    }
    else {
        vector<float> labDipoleVec;
        labDipoles.download(labDipoleVec);
        for (int i = 0; i < numParticles; i++)
            dipoles[order[i]] = Vec3(labDipoleVec[3*i], labDipoleVec[3*i+1], labDipoleVec[3*i+2]);
    }
}


void CommonCalcAmoebaMultipoleForceKernel::getInducedDipoles(ContextImpl& context, vector<Vec3>& dipoles) {
    ContextSelector selector(cc);
    ensureMultipolesValid(context);
    int numParticles = cc.getNumAtoms();
    dipoles.resize(numParticles);
    const vector<int>& order = cc.getAtomIndex();
    if (cc.getUseDoublePrecision()) {
        vector<double> d;
        inducedDipole.download(d);
        for (int i = 0; i < numParticles; i++)
            dipoles[order[i]] = Vec3(d[3*i], d[3*i+1], d[3*i+2]);
    }
    else {
        vector<float> d;
        inducedDipole.download(d);
        for (int i = 0; i < numParticles; i++)
            dipoles[order[i]] = Vec3(d[3*i], d[3*i+1], d[3*i+2]);
    }
}


void CommonCalcAmoebaMultipoleForceKernel::getTotalDipoles(ContextImpl& context, vector<Vec3>& dipoles) {
    ContextSelector selector(cc);
    ensureMultipolesValid(context);
    int numParticles = cc.getNumAtoms();
    dipoles.resize(numParticles);
    const vector<int>& order = cc.getAtomIndex();
    if (cc.getUseDoublePrecision()) {
        vector<mm_double4> posqVec;
        vector<double> labDipoleVec;
        vector<double> inducedDipoleVec;
        double totalDipoleVecX;
        double totalDipoleVecY;
        double totalDipoleVecZ;
        inducedDipole.download(inducedDipoleVec);
        labDipoles.download(labDipoleVec);
        cc.getPosq().download(posqVec);
        for (int i = 0; i < numParticles; i++) {
            totalDipoleVecX = labDipoleVec[3*i] + inducedDipoleVec[3*i];
            totalDipoleVecY = labDipoleVec[3*i+1] + inducedDipoleVec[3*i+1];
            totalDipoleVecZ = labDipoleVec[3*i+2] + inducedDipoleVec[3*i+2];
            dipoles[order[i]] = Vec3(totalDipoleVecX, totalDipoleVecY, totalDipoleVecZ);
        }
    }
    else {
        vector<mm_float4> posqVec;
        vector<float> labDipoleVec;
        vector<float> inducedDipoleVec;
        float totalDipoleVecX;
        float totalDipoleVecY;
        float totalDipoleVecZ;
        inducedDipole.download(inducedDipoleVec);
        labDipoles.download(labDipoleVec);
        cc.getPosq().download(posqVec);
        for (int i = 0; i < numParticles; i++) {
            totalDipoleVecX = labDipoleVec[3*i] + inducedDipoleVec[3*i];
            totalDipoleVecY = labDipoleVec[3*i+1] + inducedDipoleVec[3*i+1];
            totalDipoleVecZ = labDipoleVec[3*i+2] + inducedDipoleVec[3*i+2];
            dipoles[order[i]] = Vec3(totalDipoleVecX, totalDipoleVecY, totalDipoleVecZ);
        }
    }
}

void CommonCalcAmoebaMultipoleForceKernel::getElectrostaticPotential(ContextImpl& context, const vector<Vec3>& inputGrid, vector<double>& outputElectrostaticPotential) {
    ContextSelector selector(cc);
    ensureMultipolesValid(context);
    int numPoints = inputGrid.size();
    int elementSize = (cc.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
    ComputeArray points, potential;
    points.initialize(cc, numPoints, 4*elementSize, "points");
    potential.initialize(cc, numPoints, elementSize, "potential");
    
    // Copy the grid points to the GPU.
    
    if (cc.getUseDoublePrecision()) {
        vector<mm_double4> p(numPoints);
        for (int i = 0; i < numPoints; i++)
            p[i] = mm_double4(inputGrid[i][0], inputGrid[i][1], inputGrid[i][2], 0);
        points.upload(p);
    }
    else {
        vector<mm_float4> p(numPoints);
        for (int i = 0; i < numPoints; i++)
            p[i] = mm_float4((float) inputGrid[i][0], (float) inputGrid[i][1], (float) inputGrid[i][2], 0);
        points.upload(p);
    }
    
    // Compute the potential.
    
    computePotentialKernel->setArg(4, points);
    computePotentialKernel->setArg(5, potential);
    computePotentialKernel->setArg(6, numPoints);
    setPeriodicBoxArgs(cc, computePotentialKernel, 7);
    computePotentialKernel->execute(numPoints, 128);
    outputElectrostaticPotential.resize(numPoints);
    if (cc.getUseDoublePrecision())
        potential.download(outputElectrostaticPotential);
    else {
        vector<float> p(numPoints);
        potential.download(p);
        for (int i = 0; i < numPoints; i++)
            outputElectrostaticPotential[i] = p[i];
    }
}

template <class T, class T4, class M4>
void CommonCalcAmoebaMultipoleForceKernel::computeSystemMultipoleMoments(ContextImpl& context, vector<double>& outputMultipoleMoments) {
    // Compute the local coordinates relative to the center of mass.
    int numAtoms = cc.getNumAtoms();
    vector<T4> posq;
    vector<M4> velm;
    cc.getPosq().download(posq);
    cc.getVelm().download(velm);
    double totalMass = 0.0;
    Vec3 centerOfMass(0, 0, 0);
    for (int i = 0; i < numAtoms; i++) {
        double mass = (velm[i].w > 0 ? 1.0/velm[i].w : 0.0);
        totalMass += mass;
        centerOfMass[0] += mass*posq[i].x;
        centerOfMass[1] += mass*posq[i].y;
        centerOfMass[2] += mass*posq[i].z;
    }
    if (totalMass > 0.0) {
        centerOfMass[0] /= totalMass;
        centerOfMass[1] /= totalMass;
        centerOfMass[2] /= totalMass;
    }
    vector<mm_double4> posqLocal(numAtoms);
    for (int i = 0; i < numAtoms; i++) {
        posqLocal[i].x = posq[i].x - centerOfMass[0];
        posqLocal[i].y = posq[i].y - centerOfMass[1];
        posqLocal[i].z = posq[i].z - centerOfMass[2];
        posqLocal[i].w = posq[i].w;
    }

    // Compute the multipole moments.
    
    double totalCharge = 0.0;
    double xdpl = 0.0;
    double ydpl = 0.0;
    double zdpl = 0.0;
    double xxqdp = 0.0;
    double xyqdp = 0.0;
    double xzqdp = 0.0;
    double yxqdp = 0.0;
    double yyqdp = 0.0;
    double yzqdp = 0.0;
    double zxqdp = 0.0;
    double zyqdp = 0.0;
    double zzqdp = 0.0;
    vector<T> labDipoleVec, inducedDipoleVec, quadrupoleVec;
    labDipoles.download(labDipoleVec);
    inducedDipole.download(inducedDipoleVec);
    labQuadrupoles.download(quadrupoleVec);
    for (int i = 0; i < numAtoms; i++) {
        totalCharge += posqLocal[i].w;
        double netDipoleX = (labDipoleVec[3*i] + inducedDipoleVec[3*i]);
        double netDipoleY = (labDipoleVec[3*i+1] + inducedDipoleVec[3*i+1]);
        double netDipoleZ = (labDipoleVec[3*i+2] + inducedDipoleVec[3*i+2]);
        xdpl += posqLocal[i].x*posqLocal[i].w + netDipoleX;
        ydpl += posqLocal[i].y*posqLocal[i].w + netDipoleY;
        zdpl += posqLocal[i].z*posqLocal[i].w + netDipoleZ;
        xxqdp += posqLocal[i].x*posqLocal[i].x*posqLocal[i].w + 2*posqLocal[i].x*netDipoleX;
        xyqdp += posqLocal[i].x*posqLocal[i].y*posqLocal[i].w + posqLocal[i].x*netDipoleY + posqLocal[i].y*netDipoleX;
        xzqdp += posqLocal[i].x*posqLocal[i].z*posqLocal[i].w + posqLocal[i].x*netDipoleZ + posqLocal[i].z*netDipoleX;
        yxqdp += posqLocal[i].y*posqLocal[i].x*posqLocal[i].w + posqLocal[i].y*netDipoleX + posqLocal[i].x*netDipoleY;
        yyqdp += posqLocal[i].y*posqLocal[i].y*posqLocal[i].w + 2*posqLocal[i].y*netDipoleY;
        yzqdp += posqLocal[i].y*posqLocal[i].z*posqLocal[i].w + posqLocal[i].y*netDipoleZ + posqLocal[i].z*netDipoleY;
        zxqdp += posqLocal[i].z*posqLocal[i].x*posqLocal[i].w + posqLocal[i].z*netDipoleX + posqLocal[i].x*netDipoleZ;
        zyqdp += posqLocal[i].z*posqLocal[i].y*posqLocal[i].w + posqLocal[i].z*netDipoleY + posqLocal[i].y*netDipoleZ;
        zzqdp += posqLocal[i].z*posqLocal[i].z*posqLocal[i].w + 2*posqLocal[i].z*netDipoleZ;
    }

    // Convert the quadrupole from traced to traceless form.
 
    double qave = (xxqdp + yyqdp + zzqdp)/3;
    xxqdp = 1.5*(xxqdp-qave);
    xyqdp = 1.5*xyqdp;
    xzqdp = 1.5*xzqdp;
    yxqdp = 1.5*yxqdp;
    yyqdp = 1.5*(yyqdp-qave);
    yzqdp = 1.5*yzqdp;
    zxqdp = 1.5*zxqdp;
    zyqdp = 1.5*zyqdp;
    zzqdp = 1.5*(zzqdp-qave);

    // Add the traceless atomic quadrupoles to the total quadrupole moment.

    for (int i = 0; i < numAtoms; i++) {
        xxqdp = xxqdp + 3*quadrupoleVec[5*i];
        xyqdp = xyqdp + 3*quadrupoleVec[5*i+1];
        xzqdp = xzqdp + 3*quadrupoleVec[5*i+2];
        yxqdp = yxqdp + 3*quadrupoleVec[5*i+1];
        yyqdp = yyqdp + 3*quadrupoleVec[5*i+3];
        yzqdp = yzqdp + 3*quadrupoleVec[5*i+4];
        zxqdp = zxqdp + 3*quadrupoleVec[5*i+2];
        zyqdp = zyqdp + 3*quadrupoleVec[5*i+4];
        zzqdp = zzqdp + -3*(quadrupoleVec[5*i]+quadrupoleVec[5*i+3]);
    }
 
    double debye = 4.80321;
    outputMultipoleMoments.resize(13);
    outputMultipoleMoments[0] = totalCharge;
    outputMultipoleMoments[1] = 10.0*xdpl*debye;
    outputMultipoleMoments[2] = 10.0*ydpl*debye;
    outputMultipoleMoments[3] = 10.0*zdpl*debye;
    outputMultipoleMoments[4] = 100.0*xxqdp*debye;
    outputMultipoleMoments[5] = 100.0*xyqdp*debye;
    outputMultipoleMoments[6] = 100.0*xzqdp*debye;
    outputMultipoleMoments[7] = 100.0*yxqdp*debye;
    outputMultipoleMoments[8] = 100.0*yyqdp*debye;
    outputMultipoleMoments[9] = 100.0*yzqdp*debye;
    outputMultipoleMoments[10] = 100.0*zxqdp*debye;
    outputMultipoleMoments[11] = 100.0*zyqdp*debye;
    outputMultipoleMoments[12] = 100.0*zzqdp*debye;
}


void CommonCalcAmoebaMultipoleForceKernel::getSystemMultipoleMoments(ContextImpl& context, vector<double>& outputMultipoleMoments) {
    ContextSelector selector(cc);
    ensureMultipolesValid(context);
    if (cc.getUseDoublePrecision())
        computeSystemMultipoleMoments<double, mm_double4, mm_double4>(context, outputMultipoleMoments);
    else if (cc.getUseMixedPrecision())
        computeSystemMultipoleMoments<float, mm_float4, mm_double4>(context, outputMultipoleMoments);
    else
        computeSystemMultipoleMoments<float, mm_float4, mm_float4>(context, outputMultipoleMoments);
}

void CommonCalcAmoebaMultipoleForceKernel::copyParametersToContext(ContextImpl& context, const AmoebaMultipoleForce& force) {
    // Make sure the new parameters are acceptable.
    
    ContextSelector selector(cc);
    if (force.getNumMultipoles() != cc.getNumAtoms())
        throw OpenMMException("updateParametersInContext: The number of multipoles has changed");
    
    // Record the per-multipole parameters.
    
    cc.getPosq().download(cc.getPinnedBuffer());
    mm_float4* posqf = (mm_float4*) cc.getPinnedBuffer();
    mm_double4* posqd = (mm_double4*) cc.getPinnedBuffer();
    vector<mm_float2> dampingAndTholeVec;
    vector<float> polarizabilityVec;
    vector<float> localDipolesVec;
    vector<float> localQuadrupolesVec;
    vector<mm_int4> multipoleParticlesVec;
    for (int i = 0; i < force.getNumMultipoles(); i++) {
        double charge, thole, damping, polarity;
        int axisType, atomX, atomY, atomZ;
        vector<double> dipole, quadrupole;
        force.getMultipoleParameters(i, charge, dipole, quadrupole, axisType, atomZ, atomX, atomY, thole, damping, polarity);
        if (cc.getUseDoublePrecision())
            posqd[i].w = charge;
        else
            posqf[i].w = (float) charge;
        dampingAndTholeVec.push_back(mm_float2((float) damping, (float) thole));
        polarizabilityVec.push_back((float) polarity);
        multipoleParticlesVec.push_back(mm_int4(atomX, atomY, atomZ, axisType));
        for (int j = 0; j < 3; j++)
            localDipolesVec.push_back((float) dipole[j]);
        localQuadrupolesVec.push_back((float) quadrupole[0]);
        localQuadrupolesVec.push_back((float) quadrupole[1]);
        localQuadrupolesVec.push_back((float) quadrupole[2]);
        localQuadrupolesVec.push_back((float) quadrupole[4]);
        localQuadrupolesVec.push_back((float) quadrupole[5]);
    }
    if (!hasQuadrupoles) {
        for (auto q : localQuadrupolesVec)
            if (q != 0.0)
                throw OpenMMException("updateParametersInContext: Cannot set a non-zero quadrupole moment, because quadrupoles were excluded from the kernel");
    }
    for (int i = force.getNumMultipoles(); i < cc.getPaddedNumAtoms(); i++) {
        dampingAndTholeVec.push_back(mm_float2(0, 0));
        polarizabilityVec.push_back(0);
        multipoleParticlesVec.push_back(mm_int4(0, 0, 0, 0));
        for (int j = 0; j < 3; j++)
            localDipolesVec.push_back(0);
        for (int j = 0; j < 5; j++)
            localQuadrupolesVec.push_back(0);
    }
    dampingAndThole.upload(dampingAndTholeVec);
    polarizability.upload(polarizabilityVec);
    multipoleParticles.upload(multipoleParticlesVec);
    localDipoles.upload(localDipolesVec);
    localQuadrupoles.upload(localQuadrupolesVec);
    cc.getPosq().upload(cc.getPinnedBuffer());
    cc.invalidateMolecules();
    multipolesAreValid = false;
}

void CommonCalcAmoebaMultipoleForceKernel::getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    if (!usePME)
        throw OpenMMException("getPMEParametersInContext: This Context is not using PME");
    alpha = pmeAlpha;
    nx = gridSizeX;
    ny = gridSizeY;
    nz = gridSizeZ;
}

/* -------------------------------------------------------------------------- *
 *                       AmoebaGeneralizedKirkwood                            *
 * -------------------------------------------------------------------------- */

class CommonCalcAmoebaGeneralizedKirkwoodForceKernel::ForceInfo : public ComputeForceInfo {
public:
    ForceInfo(const AmoebaGeneralizedKirkwoodForce& force) : force(force) {
    }
    bool areParticlesIdentical(int particle1, int particle2) {
        double charge1, charge2, radius1, radius2, scale1, scale2;
        force.getParticleParameters(particle1, charge1, radius1, scale1);
        force.getParticleParameters(particle2, charge2, radius2, scale2);
        return (charge1 == charge2 && radius1 == radius2 && scale1 == scale2);
    }
private:
    const AmoebaGeneralizedKirkwoodForce& force;
};

CommonCalcAmoebaGeneralizedKirkwoodForceKernel::CommonCalcAmoebaGeneralizedKirkwoodForceKernel(const std::string& name, const Platform& platform, ComputeContext& cc, const System& system) :
           CalcAmoebaGeneralizedKirkwoodForceKernel(name, platform), cc(cc), system(system), hasInitializedKernels(false) {
}

void CommonCalcAmoebaGeneralizedKirkwoodForceKernel::initialize(const System& system, const AmoebaGeneralizedKirkwoodForce& force) {
    ContextSelector selector(cc);
    if (cc.getNumContexts() > 1)
        throw OpenMMException("AmoebaGeneralizedKirkwoodForce does not support using multiple devices");
    const AmoebaMultipoleForce* multipoles = NULL;
    for (int i = 0; i < system.getNumForces() && multipoles == NULL; i++)
        multipoles = dynamic_cast<const AmoebaMultipoleForce*>(&system.getForce(i));
    if (multipoles == NULL)
        throw OpenMMException("AmoebaGeneralizedKirkwoodForce requires the System to also contain an AmoebaMultipoleForce");
    NonbondedUtilities& nb = cc.getNonbondedUtilities();
    int paddedNumAtoms = cc.getPaddedNumAtoms();
    int elementSize = (cc.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
    params.initialize<mm_float2>(cc, paddedNumAtoms, "amoebaGkParams");
    bornRadii.initialize(cc, paddedNumAtoms, elementSize, "bornRadii");
    field.initialize(cc, 3*paddedNumAtoms, sizeof(long long), "gkField");
    bornSum.initialize<long long>(cc, paddedNumAtoms, "bornSum");
    bornForce.initialize<long long>(cc, paddedNumAtoms, "bornForce");
    inducedDipoleS.initialize(cc, 3*paddedNumAtoms, elementSize, "inducedDipoleS");
    inducedDipolePolarS.initialize(cc, 3*paddedNumAtoms, elementSize, "inducedDipolePolarS");
    polarizationType = multipoles->getPolarizationType();
    if (polarizationType != AmoebaMultipoleForce::Direct) {
        inducedField.initialize(cc, 3*paddedNumAtoms, sizeof(long long), "gkInducedField");
        inducedFieldPolar.initialize(cc, 3*paddedNumAtoms, sizeof(long long), "gkInducedFieldPolar");
    }
    cc.addAutoclearBuffer(field);
    cc.addAutoclearBuffer(bornSum);
    cc.addAutoclearBuffer(bornForce);
    vector<mm_float2> paramsVector(paddedNumAtoms);
    for (int i = 0; i < force.getNumParticles(); i++) {
        double charge, radius, scalingFactor;
        force.getParticleParameters(i, charge, radius, scalingFactor);
        paramsVector[i] = mm_float2((float) radius, (float) (scalingFactor*radius));
        
        // Make sure the charge matches the one specified by the AmoebaMultipoleForce.
        
        double charge2, thole, damping, polarity;
        int axisType, atomX, atomY, atomZ;
        vector<double> dipole, quadrupole;
        multipoles->getMultipoleParameters(i, charge2, dipole, quadrupole, axisType, atomZ, atomX, atomY, thole, damping, polarity);
        if (charge != charge2)
            throw OpenMMException("AmoebaGeneralizedKirkwoodForce and AmoebaMultipoleForce must specify the same charge for every atom");
    }
    params.upload(paramsVector);
    
    // Select the number of threads for each kernel.
    
    double computeBornSumThreadMemory = 4*elementSize+3*sizeof(float);
    double gkForceThreadMemory = 24*elementSize;
    double chainRuleThreadMemory = 10*elementSize;
    double ediffThreadMemory = 28*elementSize+2*sizeof(float)+3*sizeof(int)/(double) cc.TileSize;
    int maxThreads = max(32, cc.getNonbondedUtilities().getForceThreadBlockSize());
    computeBornSumThreads = min(maxThreads, cc.computeThreadBlockSize(computeBornSumThreadMemory));
    gkForceThreads = min(maxThreads, cc.computeThreadBlockSize(gkForceThreadMemory));
    chainRuleThreads = min(maxThreads, cc.computeThreadBlockSize(chainRuleThreadMemory));
    ediffThreads = min(maxThreads, cc.computeThreadBlockSize(ediffThreadMemory));
    
    // Set preprocessor macros we will use when we create the kernels.
    
    defines["NUM_ATOMS"] = cc.intToString(cc.getNumAtoms());
    defines["PADDED_NUM_ATOMS"] = cc.intToString(paddedNumAtoms);
    defines["BORN_SUM_THREAD_BLOCK_SIZE"] = cc.intToString(computeBornSumThreads);
    defines["GK_FORCE_THREAD_BLOCK_SIZE"] = cc.intToString(gkForceThreads);
    defines["CHAIN_RULE_THREAD_BLOCK_SIZE"] = cc.intToString(chainRuleThreads);
    defines["EDIFF_THREAD_BLOCK_SIZE"] = cc.intToString(ediffThreads);
    defines["NUM_BLOCKS"] = cc.intToString(cc.getNumAtomBlocks());
    defines["GK_C"] = cc.doubleToString(2.455);
    double solventDielectric = force.getSolventDielectric();
    defines["GK_FC"] = cc.doubleToString(1*(1-solventDielectric)/(0+1*solventDielectric));
    defines["GK_FD"] = cc.doubleToString(2*(1-solventDielectric)/(1+2*solventDielectric));
    defines["GK_FQ"] = cc.doubleToString(3*(1-solventDielectric)/(2+3*solventDielectric));
    defines["EPSILON_FACTOR"] = cc.doubleToString(ONE_4PI_EPS0);
    defines["M_PI"] = cc.doubleToString(M_PI);
    defines["ENERGY_SCALE_FACTOR"] = cc.doubleToString(ONE_4PI_EPS0/force.getSoluteDielectric());
    if (polarizationType == AmoebaMultipoleForce::Direct)
        defines["DIRECT_POLARIZATION"] = "";
    else if (polarizationType == AmoebaMultipoleForce::Mutual)
        defines["MUTUAL_POLARIZATION"] = "";
    else if (polarizationType == AmoebaMultipoleForce::Extrapolated)
        defines["EXTRAPOLATED_POLARIZATION"] = "";
    includeSurfaceArea = force.getIncludeCavityTerm();
    if (includeSurfaceArea) {
        defines["SURFACE_AREA_FACTOR"] = cc.doubleToString(force.getSurfaceAreaFactor());
        defines["PROBE_RADIUS"] = cc.doubleToString(force.getProbeRadius());
        defines["DIELECTRIC_OFFSET"] = cc.doubleToString(0.009);
    }
    cc.addForce(new ForceInfo(force));
}

double CommonCalcAmoebaGeneralizedKirkwoodForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    // Since GK is so tightly entwined with the electrostatics, this method does nothing, and the force calculation
    // is driven by AmoebaMultipoleForce.
    return 0.0;
}

void CommonCalcAmoebaGeneralizedKirkwoodForceKernel::computeBornRadii(ComputeArray& torque, ComputeArray& labDipoles, ComputeArray& labQuadrupoles,
            ComputeArray& inducedDipole, ComputeArray& inducedDipolePolar, ComputeArray& dampingAndThole, ComputeArray& covalentFlags, ComputeArray& polarizationGroupFlags) {
    NonbondedUtilities& nb = cc.getNonbondedUtilities();
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
        
        // Create the kernels.
        
        int numExclusionTiles = nb.getExclusionTiles().getSize();
        defines["NUM_TILES_WITH_EXCLUSIONS"] = cc.intToString(numExclusionTiles);
        int numContexts = cc.getNumContexts();
        int startExclusionIndex = cc.getContextIndex()*numExclusionTiles/numContexts;
        int endExclusionIndex = (cc.getContextIndex()+1)*numExclusionTiles/numContexts;
        defines["FIRST_EXCLUSION_TILE"] = cc.intToString(startExclusionIndex);
        defines["LAST_EXCLUSION_TILE"] = cc.intToString(endExclusionIndex);
        stringstream forceSource;
        forceSource << CommonAmoebaKernelSources::amoebaGk;
        forceSource << "#define F1\n";
        forceSource << CommonAmoebaKernelSources::gkPairForce1;
        forceSource << CommonAmoebaKernelSources::gkPairForce2;
        forceSource << CommonAmoebaKernelSources::gkEDiffPairForce;
        forceSource << "#undef F1\n";
        forceSource << "#define F2\n";
        forceSource << CommonAmoebaKernelSources::gkPairForce1;
        forceSource << CommonAmoebaKernelSources::gkPairForce2;
        forceSource << "#undef F2\n";
        forceSource << "#define T1\n";
        forceSource << CommonAmoebaKernelSources::gkPairForce1;
        forceSource << CommonAmoebaKernelSources::gkPairForce2;
        forceSource << CommonAmoebaKernelSources::gkEDiffPairForce;
        forceSource << "#undef T1\n";
        forceSource << "#define T2\n";
        forceSource << CommonAmoebaKernelSources::gkPairForce1;
        forceSource << CommonAmoebaKernelSources::gkPairForce2;
        forceSource << "#undef T2\n";
        forceSource << "#define T3\n";
        forceSource << CommonAmoebaKernelSources::gkEDiffPairForce;
        forceSource << "#undef T3\n";
        forceSource << "#define B1\n";
        forceSource << "#define B2\n";
        forceSource << CommonAmoebaKernelSources::gkPairForce1;
        forceSource << CommonAmoebaKernelSources::gkPairForce2;
        ComputeProgram program = cc.compileProgram(forceSource.str(), defines);
        computeBornSumKernel = program->createKernel("computeBornSum");
        computeBornSumKernel->addArg(bornSum);
        computeBornSumKernel->addArg(cc.getPosq());
        computeBornSumKernel->addArg(params);
        computeBornSumKernel->addArg();
        reduceBornSumKernel = program->createKernel("reduceBornSum");
        reduceBornSumKernel->addArg(bornSum);
        reduceBornSumKernel->addArg(params);
        reduceBornSumKernel->addArg(bornRadii);
        gkForceKernel = program->createKernel("computeGKForces");
        gkForceKernel->addArg(cc.getLongForceBuffer());
        gkForceKernel->addArg(torque);
        gkForceKernel->addArg(cc.getEnergyBuffer());
        gkForceKernel->addArg(cc.getPosq());
        gkForceKernel->addArg();
        gkForceKernel->addArg();
        gkForceKernel->addArg(labDipoles);
        gkForceKernel->addArg(labQuadrupoles);
        gkForceKernel->addArg(inducedDipoleS);
        gkForceKernel->addArg(inducedDipolePolarS);
        gkForceKernel->addArg(bornRadii);
        gkForceKernel->addArg(bornForce);
        chainRuleKernel = program->createKernel("computeChainRuleForce");
        chainRuleKernel->addArg(cc.getLongForceBuffer());
        chainRuleKernel->addArg(cc.getPosq());
        chainRuleKernel->addArg();
        chainRuleKernel->addArg();
        chainRuleKernel->addArg(params);
        chainRuleKernel->addArg(bornRadii);
        chainRuleKernel->addArg(bornForce);
        ediffKernel = program->createKernel("computeEDiffForce");
        ediffKernel->addArg(cc.getLongForceBuffer());
        ediffKernel->addArg(torque);
        ediffKernel->addArg(cc.getEnergyBuffer());
        ediffKernel->addArg(cc.getPosq());
        ediffKernel->addArg(covalentFlags);
        ediffKernel->addArg(polarizationGroupFlags);
        ediffKernel->addArg(nb.getExclusionTiles());
        ediffKernel->addArg();
        ediffKernel->addArg();
        ediffKernel->addArg(labDipoles);
        ediffKernel->addArg(labQuadrupoles);
        ediffKernel->addArg(inducedDipole);
        ediffKernel->addArg(inducedDipolePolar);
        ediffKernel->addArg(inducedDipoleS);
        ediffKernel->addArg(inducedDipolePolarS);
        ediffKernel->addArg(dampingAndThole);
        if (includeSurfaceArea) {
            surfaceAreaKernel = program->createKernel("computeSurfaceAreaForce");
            surfaceAreaKernel->addArg(bornForce);
            surfaceAreaKernel->addArg(cc.getEnergyBuffer());
            surfaceAreaKernel->addArg(params);
            surfaceAreaKernel->addArg(bornRadii);
        }
    }
    computeBornSumKernel->setArg(3, nb.getNumTiles());
    int numForceThreadBlocks = nb.getNumForceThreadBlocks();
    computeBornSumKernel->execute(numForceThreadBlocks*computeBornSumThreads, computeBornSumThreads);
    reduceBornSumKernel->execute(cc.getNumAtoms());
}

void CommonCalcAmoebaGeneralizedKirkwoodForceKernel::finishComputation() {
    NonbondedUtilities& nb = cc.getNonbondedUtilities();
    int startTileIndex = nb.getStartTileIndex();
    int numTileIndices = nb.getNumTiles();
    int numForceThreadBlocks = nb.getNumForceThreadBlocks();
    
    // Compute the GK force.
    
    gkForceKernel->setArg(4, startTileIndex);
    gkForceKernel->setArg(5, numTileIndices);
    gkForceKernel->execute(numForceThreadBlocks*gkForceThreads, gkForceThreads);

    // Compute the surface area force.
    
    if (includeSurfaceArea)
        surfaceAreaKernel->execute(cc.getNumAtoms());
    
    // Apply the remaining terms.
    
    chainRuleKernel->setArg(2, startTileIndex);
    chainRuleKernel->setArg(3, numTileIndices);
    chainRuleKernel->execute(numForceThreadBlocks*chainRuleThreads, chainRuleThreads);    
    ediffKernel->setArg(7, startTileIndex);
    ediffKernel->setArg(8, numTileIndices);
    ediffKernel->execute(numForceThreadBlocks*ediffThreads, ediffThreads);
}

void CommonCalcAmoebaGeneralizedKirkwoodForceKernel::copyParametersToContext(ContextImpl& context, const AmoebaGeneralizedKirkwoodForce& force) {
    // Make sure the new parameters are acceptable.
    
    ContextSelector selector(cc);
    if (force.getNumParticles() != cc.getNumAtoms())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");
    
    // Record the per-particle parameters.
    
    vector<mm_float2> paramsVector(cc.getPaddedNumAtoms());
    for (int i = 0; i < force.getNumParticles(); i++) {
        double charge, radius, scalingFactor;
        force.getParticleParameters(i, charge, radius, scalingFactor);
        paramsVector[i] = mm_float2((float) radius, (float) (scalingFactor*radius));
    }
    params.upload(paramsVector);
    cc.invalidateMolecules();
}

/* -------------------------------------------------------------------------- *
 *                           AmoebaVdw                                        *
 * -------------------------------------------------------------------------- */

class CommonCalcAmoebaVdwForceKernel::ForceInfo : public ComputeForceInfo {
public:
    ForceInfo(const AmoebaVdwForce& force) : force(force) {
    }
    bool areParticlesIdentical(int particle1, int particle2) {
        int iv1, iv2, type1, type2;
        double sigma1, sigma2, epsilon1, epsilon2, reduction1, reduction2;
        bool isAlchemical1, isAlchemical2;
        force.getParticleParameters(particle1, iv1, sigma1, epsilon1, reduction1, isAlchemical1, type1);
        force.getParticleParameters(particle2, iv2, sigma2, epsilon2, reduction2, isAlchemical2, type2);
        return (sigma1 == sigma2 && epsilon1 == epsilon2 && reduction1 == reduction2 && isAlchemical1 == isAlchemical2 && type1 == type2);
    }
private:
    const AmoebaVdwForce& force;
};

CommonCalcAmoebaVdwForceKernel::CommonCalcAmoebaVdwForceKernel(const std::string& name, const Platform& platform, ComputeContext& cc, const System& system) :
        CalcAmoebaVdwForceKernel(name, platform), cc(cc), system(system), hasInitializedNonbonded(false), nonbonded(NULL) {
}

CommonCalcAmoebaVdwForceKernel::~CommonCalcAmoebaVdwForceKernel() {
    ContextSelector selector(cc);
    if (nonbonded != NULL)
        delete nonbonded;
}

void CommonCalcAmoebaVdwForceKernel::initialize(const System& system, const AmoebaVdwForce& force) {
    ContextSelector selector(cc);
    int paddedNumAtoms = cc.getPaddedNumAtoms();
    bondReductionAtoms.initialize<int>(cc, paddedNumAtoms, "bondReductionAtoms");
    bondReductionFactors.initialize<float>(cc, paddedNumAtoms, "bondReductionFactors");
    tempPosq.initialize(cc, paddedNumAtoms, cc.getUseDoublePrecision() ? sizeof(mm_double4) : sizeof(mm_float4), "tempPosq");
    tempForces.initialize<long long>(cc, 3*paddedNumAtoms, "tempForces");
    
    // Record atom parameters.
    vector<int> atomTypeVec;
    vector<vector<double> > sigmaMatrix, epsilonMatrix;
    AmoebaVdwForceImpl::createParameterMatrix(force, atomTypeVec, sigmaMatrix, epsilonMatrix);
    atomTypeVec.resize(paddedNumAtoms, 0);
    int numTypes = sigmaMatrix.size();
    atomType.initialize<int>(cc, paddedNumAtoms, "atomType");
    sigmaEpsilon.initialize<mm_float2>(cc, numTypes*numTypes, "sigmaEpsilon");
    vector<mm_float2> sigmaEpsilonVec(sigmaEpsilon.getSize());
    for (int i = 0; i < numTypes; i++)
        for (int j = 0; j < numTypes; j++)
            sigmaEpsilonVec[i*numTypes+j] = mm_float2((float) sigmaMatrix[i][j], (float) epsilonMatrix[i][j]);
    atomType.upload(atomTypeVec);
    sigmaEpsilon.upload(sigmaEpsilonVec);
    
    vector<float> isAlchemicalVec(paddedNumAtoms, 0);
    vector<int> bondReductionAtomsVec(paddedNumAtoms, 0);
    vector<float> bondReductionFactorsVec(paddedNumAtoms, 0);
    vector<vector<int> > exclusions(cc.getNumAtoms());

    // Handle Alchemical parameters.
    hasAlchemical = force.getAlchemicalMethod() != AmoebaVdwForce::None;
    if (hasAlchemical) {
       isAlchemical.initialize<float>(cc, paddedNumAtoms, "isAlchemical");
       vdwLambda.initialize<float>(cc, 1, "vdwLambda");
    }

    for (int i = 0; i < force.getNumParticles(); i++) {
        int ivIndex, type;
        double sigma, epsilon, reductionFactor;
        bool alchemical;
        force.getParticleParameters(i, ivIndex, sigma, epsilon, reductionFactor, alchemical, type);
        isAlchemicalVec[i] = (alchemical) ? 1.0f : 0.0f;
        bondReductionAtomsVec[i] = ivIndex;
        bondReductionFactorsVec[i] = (float) reductionFactor;
        force.getParticleExclusions(i, exclusions[i]);
        exclusions[i].push_back(i);
    }
    bondReductionAtoms.upload(bondReductionAtomsVec);
    bondReductionFactors.upload(bondReductionFactorsVec);
    if (force.getUseDispersionCorrection())
        dispersionCoefficient = AmoebaVdwForceImpl::calcDispersionCorrection(system, force);
    else
        dispersionCoefficient = 0.0;               

    // This force is applied based on modified atom positions, where hydrogens have been moved slightly
    // closer to their parent atoms.  We therefore create a separate NonbondedUtilities just for
    // this force, so it will have its own neighbor list and interaction kernel.
    
    nonbonded = cc.createNonbondedUtilities();
    nonbonded->addParameter(ComputeParameterInfo(atomType, "atomType", "int", 1));
    nonbonded->addArgument(ComputeParameterInfo(sigmaEpsilon, "sigmaEpsilon", "float", 2));

    if (hasAlchemical) {
       isAlchemical.upload(isAlchemicalVec);
       currentVdwLambda = 1.0f;
       vdwLambda.upload(&currentVdwLambda);
       nonbonded->addParameter(ComputeParameterInfo(isAlchemical, "isAlchemical", "float", 1));
       nonbonded->addArgument(ComputeParameterInfo(vdwLambda, "vdwLambda", "float", 1));
    }
    
    // Create the interaction kernel.
    
    map<string, string> replacements;
    replacements["VDW_ALCHEMICAL_METHOD"] = cc.intToString(force.getAlchemicalMethod()); 
    replacements["VDW_SOFTCORE_POWER"] = cc.intToString(force.getSoftcorePower());
    replacements["VDW_SOFTCORE_ALPHA"] = cc.doubleToString(force.getSoftcoreAlpha()); 
    replacements["POTENTIAL_FUNCTION"] = cc.intToString(force.getPotentialFunction());
    replacements["NUM_TYPES"] = cc.intToString(numTypes);

    double cutoff = force.getCutoffDistance();
    double taperCutoff = cutoff*0.9;
    replacements["CUTOFF_DISTANCE"] = cc.doubleToString(force.getCutoffDistance());
    replacements["TAPER_CUTOFF"] = cc.doubleToString(taperCutoff);
    replacements["TAPER_C3"] = cc.doubleToString(10/pow(taperCutoff-cutoff, 3.0));
    replacements["TAPER_C4"] = cc.doubleToString(15/pow(taperCutoff-cutoff, 4.0));
    replacements["TAPER_C5"] = cc.doubleToString(6/pow(taperCutoff-cutoff, 5.0));
    bool useCutoff = (force.getNonbondedMethod() != AmoebaVdwForce::NoCutoff);
    nonbonded->addInteraction(useCutoff, useCutoff, true, force.getCutoffDistance(), exclusions,
        cc.replaceStrings(CommonAmoebaKernelSources::amoebaVdwForce2, replacements), 0);
    
    // Create the other kernels.
    
    map<string, string> defines;
    defines["PADDED_NUM_ATOMS"] = cc.intToString(paddedNumAtoms);
    ComputeProgram program = cc.compileProgram(CommonAmoebaKernelSources::amoebaVdwForce1, defines);
    prepareKernel = program->createKernel("prepareToComputeForce");
    prepareKernel->addArg(cc.getLongForceBuffer());
    prepareKernel->addArg(cc.getPosq());
    prepareKernel->addArg(tempPosq);
    prepareKernel->addArg(bondReductionAtoms);
    prepareKernel->addArg(bondReductionFactors);
    spreadKernel = program->createKernel("spreadForces");
    spreadKernel->addArg(cc.getLongForceBuffer());
    spreadKernel->addArg(tempForces);
    spreadKernel->addArg(bondReductionAtoms);
    spreadKernel->addArg(bondReductionFactors);
    cc.addForce(new ForceInfo(force));
}

double CommonCalcAmoebaVdwForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    ContextSelector selector(cc);
    if (!hasInitializedNonbonded) {
        hasInitializedNonbonded = true;
        nonbonded->initialize(system);
    }

    if (hasAlchemical) {
       float contextLambda = context.getParameter(AmoebaVdwForce::Lambda());
       if (contextLambda != currentVdwLambda) {
          vdwLambda.upload(&contextLambda);
          currentVdwLambda = contextLambda;
       }
    }

    cc.getPosq().copyTo(tempPosq);
    cc.getLongForceBuffer().copyTo(tempForces);
    prepareKernel->execute(cc.getPaddedNumAtoms());
    nonbonded->prepareInteractions(1);
    nonbonded->computeInteractions(1, includeForces, includeEnergy);
    spreadKernel->execute(cc.getPaddedNumAtoms());
    tempPosq.copyTo(cc.getPosq());
    tempForces.copyTo(cc.getLongForceBuffer());
    Vec3 a, b, c;
    cc.getPeriodicBoxVectors(a, b, c);
    return dispersionCoefficient/(a[0]*b[1]*c[2]);
}

void CommonCalcAmoebaVdwForceKernel::copyParametersToContext(ContextImpl& context, const AmoebaVdwForce& force) {
    // Make sure the new parameters are acceptable.
    
    ContextSelector selector(cc);
    if (force.getNumParticles() != cc.getNumAtoms())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");
    
    vector<int> atomTypeVec;
    vector<vector<double> > sigmaMatrix, epsilonMatrix;
    AmoebaVdwForceImpl::createParameterMatrix(force, atomTypeVec, sigmaMatrix, epsilonMatrix);
    atomTypeVec.resize(cc.getPaddedNumAtoms(), 0);
    int numTypes = sigmaMatrix.size();
    if (sigmaEpsilon.getSize() != numTypes*numTypes)
        throw OpenMMException("updateParametersInContext: The number of particle types has changed");
    vector<mm_float2> sigmaEpsilonVec(sigmaEpsilon.getSize());
    for (int i = 0; i < numTypes; i++)
        for (int j = 0; j < numTypes; j++)
            sigmaEpsilonVec[i*numTypes+j] = mm_float2((float) sigmaMatrix[i][j], (float) epsilonMatrix[i][j]);
    atomType.upload(atomTypeVec);
    sigmaEpsilon.upload(sigmaEpsilonVec);

    // Record the per-particle parameters.
    vector<float> isAlchemicalVec(cc.getPaddedNumAtoms(), 0);
    vector<int> bondReductionAtomsVec(cc.getPaddedNumAtoms(), 0);
    vector<float> bondReductionFactorsVec(cc.getPaddedNumAtoms(), 0);
    for (int i = 0; i < force.getNumParticles(); i++) {
        int ivIndex, type;
        double sigma, epsilon, reductionFactor;
        bool alchemical;
        force.getParticleParameters(i, ivIndex, sigma, epsilon, reductionFactor, alchemical, type);
        isAlchemicalVec[i] = (alchemical) ? 1.0f : 0.0f;
        bondReductionAtomsVec[i] = ivIndex;
        bondReductionFactorsVec[i] = (float) reductionFactor;
    }
    if (hasAlchemical) isAlchemical.upload(isAlchemicalVec);
    bondReductionAtoms.upload(bondReductionAtomsVec);
    bondReductionFactors.upload(bondReductionFactorsVec);
    if (force.getUseDispersionCorrection())
        dispersionCoefficient = AmoebaVdwForceImpl::calcDispersionCorrection(system, force);
    else
        dispersionCoefficient = 0.0;               
    cc.invalidateMolecules();
}

/* -------------------------------------------------------------------------- *
 *                           AmoebaWcaDispersion                              *
 * -------------------------------------------------------------------------- */

class CommonCalcAmoebaWcaDispersionForceKernel::ForceInfo : public ComputeForceInfo {
public:
    ForceInfo(const AmoebaWcaDispersionForce& force) : force(force) {
    }
    bool areParticlesIdentical(int particle1, int particle2) {
        double radius1, radius2, epsilon1, epsilon2;
        force.getParticleParameters(particle1, radius1, epsilon1);
        force.getParticleParameters(particle2, radius2, epsilon2);
        return (radius1 == radius2 && epsilon1 == epsilon2);
    }
private:
    const AmoebaWcaDispersionForce& force;
};

CommonCalcAmoebaWcaDispersionForceKernel::CommonCalcAmoebaWcaDispersionForceKernel(const std::string& name, const Platform& platform, ComputeContext& cc, const System& system) :
           CalcAmoebaWcaDispersionForceKernel(name, platform), cc(cc), system(system) {
}

void CommonCalcAmoebaWcaDispersionForceKernel::initialize(const System& system, const AmoebaWcaDispersionForce& force) {
    int numParticles = system.getNumParticles();
    int paddedNumAtoms = cc.getPaddedNumAtoms();
    
    // Record parameters.
    
    ContextSelector selector(cc);
    vector<mm_float2> radiusEpsilonVec(paddedNumAtoms, mm_float2(0, 0));
    for (int i = 0; i < numParticles; i++) {
        double radius, epsilon;
        force.getParticleParameters(i, radius, epsilon);
        radiusEpsilonVec[i] = mm_float2((float) radius, (float) epsilon);
    }
    radiusEpsilon.initialize<mm_float2>(cc, paddedNumAtoms, "radiusEpsilon");
    radiusEpsilon.upload(radiusEpsilonVec);
    
    // Create the kernel.
    
    forceThreadBlockSize = max(32, cc.getNonbondedUtilities().getForceThreadBlockSize());
    map<string, string> defines;
    defines["NUM_ATOMS"] = cc.intToString(numParticles);
    defines["PADDED_NUM_ATOMS"] = cc.intToString(cc.getPaddedNumAtoms());
    defines["THREAD_BLOCK_SIZE"] = cc.intToString(forceThreadBlockSize);
    defines["NUM_BLOCKS"] = cc.intToString(cc.getNumAtomBlocks());
    defines["EPSO"] = cc.doubleToString(force.getEpso());
    defines["EPSH"] = cc.doubleToString(force.getEpsh());
    defines["RMINO"] = cc.doubleToString(force.getRmino());
    defines["RMINH"] = cc.doubleToString(force.getRminh());
    defines["AWATER"] = cc.doubleToString(force.getAwater());
    defines["SHCTD"] = cc.doubleToString(force.getShctd());
    defines["M_PI"] = cc.doubleToString(M_PI);
    ComputeProgram program = cc.compileProgram(CommonAmoebaKernelSources::amoebaWcaForce, defines);
    forceKernel = program->createKernel("computeWCAForce");
    forceKernel->addArg(cc.getLongForceBuffer());
    forceKernel->addArg(cc.getEnergyBuffer());
    forceKernel->addArg(cc.getPosq());
    forceKernel->addArg();
    forceKernel->addArg();
    forceKernel->addArg(radiusEpsilon);
    totalMaximumDispersionEnergy = AmoebaWcaDispersionForceImpl::getTotalMaximumDispersionEnergy(force);

    // Add an interaction to the default nonbonded kernel.  This doesn't actually do any calculations.  It's
    // just so that NonbondedUtilities will keep track of the tiles.
    
    vector<vector<int> > exclusions;
    cc.getNonbondedUtilities().addInteraction(false, false, false, 1.0, exclusions, "", force.getForceGroup());
    cc.addForce(new ForceInfo(force));
}

double CommonCalcAmoebaWcaDispersionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    ContextSelector selector(cc);
    NonbondedUtilities& nb = cc.getNonbondedUtilities();
    int startTileIndex = nb.getStartTileIndex();
    int numTileIndices = nb.getNumTiles();
    int numForceThreadBlocks = nb.getNumForceThreadBlocks();
    forceKernel->setArg(3, startTileIndex);
    forceKernel->setArg(4, numTileIndices);
    forceKernel->execute(numForceThreadBlocks*forceThreadBlockSize, forceThreadBlockSize);
    return totalMaximumDispersionEnergy;
}

void CommonCalcAmoebaWcaDispersionForceKernel::copyParametersToContext(ContextImpl& context, const AmoebaWcaDispersionForce& force) {
    // Make sure the new parameters are acceptable.
    
    ContextSelector selector(cc);
    if (force.getNumParticles() != cc.getNumAtoms())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");
    
    // Record the per-particle parameters.
    
    vector<mm_float2> radiusEpsilonVec(cc.getPaddedNumAtoms(), mm_float2(0, 0));
    for (int i = 0; i < cc.getNumAtoms(); i++) {
        double radius, epsilon;
        force.getParticleParameters(i, radius, epsilon);
        radiusEpsilonVec[i] = mm_float2((float) radius, (float) epsilon);
    }
    radiusEpsilon.upload(radiusEpsilonVec);
    totalMaximumDispersionEnergy = AmoebaWcaDispersionForceImpl::getTotalMaximumDispersionEnergy(force);
    cc.invalidateMolecules();
}


/* -------------------------------------------------------------------------- *
 *                           HippoNonbondedForce                              *
 * -------------------------------------------------------------------------- */

class CommonCalcHippoNonbondedForceKernel::ForceInfo : public ComputeForceInfo {
public:
    ForceInfo(const HippoNonbondedForce& force) : force(force) {
    }
    bool areParticlesIdentical(int particle1, int particle2) {
        double charge1, coreCharge1, alpha1, epsilon1, damping1, c61, pauliK1, pauliQ1, pauliAlpha1, polarizability1;
        double charge2, coreCharge2, alpha2, epsilon2, damping2, c62, pauliK2, pauliQ2, pauliAlpha2, polarizability2;
        int axisType1, multipoleZ1, multipoleX1, multipoleY1;
        int axisType2, multipoleZ2, multipoleX2, multipoleY2;
        vector<double> dipole1, dipole2, quadrupole1, quadrupole2;
        force.getParticleParameters(particle1, charge1, dipole1, quadrupole1, coreCharge1, alpha1, epsilon1, damping1, c61, pauliK1, pauliQ1, pauliAlpha1,
                                    polarizability1, axisType1, multipoleZ1, multipoleX1, multipoleY1);
        force.getParticleParameters(particle2, charge2, dipole2, quadrupole2, coreCharge2, alpha2, epsilon2, damping2, c62, pauliK2, pauliQ2, pauliAlpha2,
                                    polarizability2, axisType2, multipoleZ2, multipoleX2, multipoleY2);
        if (charge1 != charge2 || coreCharge1 != coreCharge2 || alpha1 != alpha2 || epsilon1 != epsilon1 || damping1 != damping2 || c61 != c62 ||
                pauliK1 != pauliK2 || pauliQ1 != pauliQ2 || pauliAlpha1 != pauliAlpha2 || polarizability1 != polarizability2 || axisType1 != axisType2) {
            return false;
        }
        for (int i = 0; i < dipole1.size(); ++i)
            if (dipole1[i] != dipole2[i])
                return false;
        for (int i = 0; i < quadrupole1.size(); ++i)
            if (quadrupole1[i] != quadrupole2[i])
                return false;
        return true;
    }
    int getNumParticleGroups() {
        return force.getNumExceptions();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        int particle1, particle2;
        double multipoleMultipoleScale, dipoleMultipoleScale, dipoleDipoleScale, dispersionScale, repulsionScale, chargeTransferScale;
        force.getExceptionParameters(index, particle1, particle2, multipoleMultipoleScale, dipoleMultipoleScale, dipoleDipoleScale, dispersionScale, repulsionScale, chargeTransferScale);
        particles.resize(2);
        particles[0] = particle1;
        particles[1] = particle2;
    }
    bool areGroupsIdentical(int group1, int group2) {
        int particle1, particle2;
        double multipoleMultipoleScale1, dipoleMultipoleScale1, dipoleDipoleScale1, dispersionScale1, repulsionScale1, chargeTransferScale1;
        double multipoleMultipoleScale2, dipoleMultipoleScale2, dipoleDipoleScale2, dispersionScale2, repulsionScale2, chargeTransferScale2;
        force.getExceptionParameters(group1, particle1, particle2, multipoleMultipoleScale1, dipoleMultipoleScale1, dipoleDipoleScale1, dispersionScale1, repulsionScale1, chargeTransferScale1);
        force.getExceptionParameters(group2, particle1, particle2, multipoleMultipoleScale2, dipoleMultipoleScale2, dipoleDipoleScale2, dispersionScale2, repulsionScale2, chargeTransferScale2);
        return (multipoleMultipoleScale1 == multipoleMultipoleScale2 && dipoleMultipoleScale1 == dipoleMultipoleScale2 &&
                dipoleDipoleScale1 == dipoleDipoleScale2 && dispersionScale1 == dispersionScale2 && repulsionScale1 == repulsionScale2 && chargeTransferScale1 == chargeTransferScale2);
    }
private:
    const HippoNonbondedForce& force;
};

class CommonCalcHippoNonbondedForceKernel::TorquePostComputation : public ComputeContext::ForcePostComputation {
public:
    TorquePostComputation(CommonCalcHippoNonbondedForceKernel& owner) : owner(owner) {
    }
    double computeForceAndEnergy(bool includeForces, bool includeEnergy, int groups) {
        owner.addTorquesToForces();
        return 0.0;
    }
private:
    CommonCalcHippoNonbondedForceKernel& owner;
};

CommonCalcHippoNonbondedForceKernel::CommonCalcHippoNonbondedForceKernel(const std::string& name, const Platform& platform, ComputeContext& cc, const System& system) :
        CalcHippoNonbondedForceKernel(name, platform), cc(cc), system(system), hasInitializedKernels(false), multipolesAreValid(false) {
}

void CommonCalcHippoNonbondedForceKernel::initialize(const System& system, const HippoNonbondedForce& force) {
    ContextSelector selector(cc);
    extrapolationCoefficients = force.getExtrapolationCoefficients();
    usePME = (force.getNonbondedMethod() == HippoNonbondedForce::PME);

    // Initialize particle parameters.

    numParticles = force.getNumParticles();
    vector<double> coreChargeVec, valenceChargeVec, alphaVec, epsilonVec, dampingVec, c6Vec, pauliKVec, pauliQVec, pauliAlphaVec, polarizabilityVec;
    vector<double> localDipolesVec, localQuadrupolesVec;
    vector<mm_int4> multipoleParticlesVec;
    vector<vector<int> > exclusions(numParticles);
    for (int i = 0; i < numParticles; i++) {
        double charge, coreCharge, alpha, epsilon, damping, c6, pauliK, pauliQ, pauliAlpha, polarizability;
        int axisType, atomX, atomY, atomZ;
        vector<double> dipole, quadrupole;
        force.getParticleParameters(i, charge, dipole, quadrupole, coreCharge, alpha, epsilon, damping, c6, pauliK, pauliQ, pauliAlpha,
                                    polarizability, axisType, atomZ, atomX, atomY);
        coreChargeVec.push_back(coreCharge);
        valenceChargeVec.push_back(charge-coreCharge);
        alphaVec.push_back(alpha);
        epsilonVec.push_back(epsilon);
        dampingVec.push_back(damping);
        c6Vec.push_back(c6);
        pauliKVec.push_back(pauliK);
        pauliQVec.push_back(pauliQ);
        pauliAlphaVec.push_back(pauliAlpha);
        polarizabilityVec.push_back(polarizability);
        multipoleParticlesVec.push_back(mm_int4(atomX, atomY, atomZ, axisType));
        for (int j = 0; j < 3; j++)
            localDipolesVec.push_back(dipole[j]);
        localQuadrupolesVec.push_back(quadrupole[0]);
        localQuadrupolesVec.push_back(quadrupole[1]);
        localQuadrupolesVec.push_back(quadrupole[2]);
        localQuadrupolesVec.push_back(quadrupole[4]);
        localQuadrupolesVec.push_back(quadrupole[5]);
        exclusions[i].push_back(i);
    }
    int paddedNumAtoms = cc.getPaddedNumAtoms();
    for (int i = numParticles; i < paddedNumAtoms; i++) {
        coreChargeVec.push_back(0);
        valenceChargeVec.push_back(0);
        alphaVec.push_back(0);
        epsilonVec.push_back(0);
        dampingVec.push_back(0);
        c6Vec.push_back(0);
        pauliKVec.push_back(0);
        pauliQVec.push_back(0);
        pauliAlphaVec.push_back(0);
        polarizabilityVec.push_back(0);
        multipoleParticlesVec.push_back(mm_int4(0, 0, 0, 0));
        for (int j = 0; j < 3; j++)
            localDipolesVec.push_back(0);
        for (int j = 0; j < 5; j++)
            localQuadrupolesVec.push_back(0);
    }
    int elementSize = (cc.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
    coreCharge.initialize(cc, paddedNumAtoms, elementSize, "coreCharge");
    valenceCharge.initialize(cc, paddedNumAtoms, elementSize, "valenceCharge");
    alpha.initialize(cc, paddedNumAtoms, elementSize, "alpha");
    epsilon.initialize(cc, paddedNumAtoms, elementSize, "epsilon");
    damping.initialize(cc, paddedNumAtoms, elementSize, "damping");
    c6.initialize(cc, paddedNumAtoms, elementSize, "c6");
    pauliK.initialize(cc, paddedNumAtoms, elementSize, "pauliK");
    pauliQ.initialize(cc, paddedNumAtoms, elementSize, "pauliQ");
    pauliAlpha.initialize(cc, paddedNumAtoms, elementSize, "pauliAlpha");
    polarizability.initialize(cc, paddedNumAtoms, elementSize, "polarizability");
    multipoleParticles.initialize<mm_int4>(cc, paddedNumAtoms, "multipoleParticles");
    localDipoles.initialize(cc, 3*paddedNumAtoms, elementSize, "localDipoles");
    localQuadrupoles.initialize(cc, 5*paddedNumAtoms, elementSize, "localQuadrupoles");
    lastPositions.initialize(cc, cc.getPosq().getSize(), cc.getPosq().getElementSize(), "lastPositions");
    coreCharge.upload(coreChargeVec, true);
    valenceCharge.upload(valenceChargeVec, true);
    alpha.upload(alphaVec, true);
    epsilon.upload(epsilonVec, true);
    damping.upload(dampingVec, true);
    c6.upload(c6Vec, true);
    pauliK.upload(pauliKVec, true);
    pauliQ.upload(pauliQVec, true);
    pauliAlpha.upload(pauliAlphaVec, true);
    polarizability.upload(polarizabilityVec, true);
    multipoleParticles.upload(multipoleParticlesVec);
    localDipoles.upload(localDipolesVec, true);
    localQuadrupoles.upload(localQuadrupolesVec, true);
    
    // Create workspace arrays.
    
    labDipoles.initialize(cc, 3*paddedNumAtoms, elementSize, "dipole");
    labQuadrupoles[0].initialize(cc, paddedNumAtoms, elementSize, "qXX");
    labQuadrupoles[1].initialize(cc, paddedNumAtoms, elementSize, "qXY");
    labQuadrupoles[2].initialize(cc, paddedNumAtoms, elementSize, "qXZ");
    labQuadrupoles[3].initialize(cc, paddedNumAtoms, elementSize, "qYY");
    labQuadrupoles[4].initialize(cc, paddedNumAtoms, elementSize, "qYZ");
    fracDipoles.initialize(cc, 3*paddedNumAtoms, elementSize, "fracDipoles");
    fracQuadrupoles.initialize(cc, 6*paddedNumAtoms, elementSize, "fracQuadrupoles");
    field.initialize(cc, 3*paddedNumAtoms, sizeof(long long), "field");
    inducedField.initialize(cc, 3*paddedNumAtoms, sizeof(long long), "inducedField");
    torque.initialize(cc, 3*paddedNumAtoms, sizeof(long long), "torque");
    inducedDipole.initialize(cc, 3*paddedNumAtoms, elementSize, "inducedDipole");
    int numOrders = extrapolationCoefficients.size();
    extrapolatedDipole.initialize(cc, 3*numParticles*numOrders, elementSize, "extrapolatedDipole");
    extrapolatedPhi.initialize(cc, 10*numParticles*numOrders, elementSize, "extrapolatedPhi");
    cc.addAutoclearBuffer(field);
    cc.addAutoclearBuffer(torque);
    
    // Record exceptions and exclusions.
    
    vector<double> exceptionScaleVec[6];
    vector<mm_int2> exceptionAtomsVec;
    for (int i = 0; i < force.getNumExceptions(); i++) {
        int particle1, particle2;
        double multipoleMultipoleScale, dipoleMultipoleScale, dipoleDipoleScale, dispersionScale, repulsionScale, chargeTransferScale;
        force.getExceptionParameters(i, particle1, particle2, multipoleMultipoleScale, dipoleMultipoleScale, dipoleDipoleScale, dispersionScale, repulsionScale, chargeTransferScale);
        exclusions[particle1].push_back(particle2);
        exclusions[particle2].push_back(particle1);
        if (usePME || multipoleMultipoleScale != 0 || dipoleMultipoleScale != 0 || dipoleDipoleScale != 0 || dispersionScale != 0 || repulsionScale != 0 || chargeTransferScale != 0) {
            exceptionAtomsVec.push_back(mm_int2(particle1, particle2));
            exceptionScaleVec[0].push_back(multipoleMultipoleScale);
            exceptionScaleVec[1].push_back(dipoleMultipoleScale);
            exceptionScaleVec[2].push_back(dipoleDipoleScale);
            exceptionScaleVec[3].push_back(dispersionScale);
            exceptionScaleVec[4].push_back(repulsionScale);
            exceptionScaleVec[5].push_back(chargeTransferScale);
        }
    }
    if (exceptionAtomsVec.size() > 0) {
        exceptionAtoms.initialize<mm_int2>(cc, exceptionAtomsVec.size(), "exceptionAtoms");
        exceptionAtoms.upload(exceptionAtomsVec);
        for (int i = 0; i < 6; i++) {
            exceptionScales[i].initialize(cc, exceptionAtomsVec.size(), elementSize, "exceptionScales");
            exceptionScales[i].upload(exceptionScaleVec[i], true);
        }
    }
    
    // Create the kernels.

    bool useShuffle = false;//(cc.getComputeCapability() >= 3.0 && !cc.getUseDoublePrecision());
    map<string, string> defines;
    defines["HIPPO"] = "1";
    defines["NUM_ATOMS"] = cc.intToString(numParticles);
    defines["PADDED_NUM_ATOMS"] = cc.intToString(cc.getPaddedNumAtoms());
    defines["NUM_BLOCKS"] = cc.intToString(cc.getNumAtomBlocks());
    defines["ENERGY_SCALE_FACTOR"] = cc.doubleToString(ONE_4PI_EPS0);
    if (useShuffle)
        defines["USE_SHUFFLE"] = "";
    maxExtrapolationOrder = extrapolationCoefficients.size();
    defines["MAX_EXTRAPOLATION_ORDER"] = cc.intToString(maxExtrapolationOrder);
    stringstream coefficients;
    for (int i = 0; i < maxExtrapolationOrder; i++) {
        if (i > 0)
            coefficients << ",";
        double sum = 0;
        for (int j = i; j < maxExtrapolationOrder; j++)
            sum += extrapolationCoefficients[j];
        coefficients << cc.doubleToString(sum);
    }
    defines["EXTRAPOLATION_COEFFICIENTS_SUM"] = coefficients.str();
    cutoff = force.getCutoffDistance();
    if (usePME) {
        int nx, ny, nz;
        force.getPMEParameters(pmeAlpha, nx, ny, nz);
        if (nx == 0 || pmeAlpha == 0) {
            NonbondedForce nb;
            nb.setEwaldErrorTolerance(force.getEwaldErrorTolerance());
            nb.setCutoffDistance(force.getCutoffDistance());
            NonbondedForceImpl::calcPMEParameters(system, nb, pmeAlpha, gridSizeX, gridSizeY, gridSizeZ, false);
            gridSizeX = cc.findLegalFFTDimension(gridSizeX);
            gridSizeY = cc.findLegalFFTDimension(gridSizeY);
            gridSizeZ = cc.findLegalFFTDimension(gridSizeZ);
        } else {
            gridSizeX = cc.findLegalFFTDimension(nx);
            gridSizeY = cc.findLegalFFTDimension(ny);
            gridSizeZ = cc.findLegalFFTDimension(nz);
        }
        force.getDPMEParameters(dpmeAlpha, nx, ny, nz);
        if (nx == 0 || dpmeAlpha == 0) {
            NonbondedForce nb;
            nb.setEwaldErrorTolerance(force.getEwaldErrorTolerance());
            nb.setCutoffDistance(force.getCutoffDistance());
            NonbondedForceImpl::calcPMEParameters(system, nb, dpmeAlpha, dispersionGridSizeX, dispersionGridSizeY, dispersionGridSizeZ, true);
            dispersionGridSizeX = cc.findLegalFFTDimension(dispersionGridSizeX);
            dispersionGridSizeY = cc.findLegalFFTDimension(dispersionGridSizeY);
            dispersionGridSizeZ = cc.findLegalFFTDimension(dispersionGridSizeZ);
        } else {
            dispersionGridSizeX = cc.findLegalFFTDimension(nx);
            dispersionGridSizeY = cc.findLegalFFTDimension(ny);
            dispersionGridSizeZ = cc.findLegalFFTDimension(nz);
        }
        defines["EWALD_ALPHA"] = cc.doubleToString(pmeAlpha);
        defines["SQRT_PI"] = cc.doubleToString(sqrt(M_PI));
        defines["USE_EWALD"] = "";
        defines["USE_CUTOFF"] = "";
        defines["USE_PERIODIC"] = "";
        defines["CUTOFF_SQUARED"] = cc.doubleToString(force.getCutoffDistance()*force.getCutoffDistance());
    }
    ComputeProgram program = cc.compileProgram(CommonAmoebaKernelSources::hippoMultipoles, defines);
    computeMomentsKernel = program->createKernel("computeLabFrameMoments");
    computeMomentsKernel->addArg(cc.getPosq());
    computeMomentsKernel->addArg(multipoleParticles);
    computeMomentsKernel->addArg(localDipoles);
    computeMomentsKernel->addArg(localQuadrupoles);
    computeMomentsKernel->addArg(labDipoles);
    computeMomentsKernel->addArg(labQuadrupoles[0]);
    computeMomentsKernel->addArg(labQuadrupoles[1]);
    computeMomentsKernel->addArg(labQuadrupoles[2]);
    computeMomentsKernel->addArg(labQuadrupoles[3]);
    computeMomentsKernel->addArg(labQuadrupoles[4]);
    recordInducedDipolesKernel = program->createKernel("recordInducedDipoles");
    recordInducedDipolesKernel->addArg(field);
    recordInducedDipolesKernel->addArg(inducedDipole);
    recordInducedDipolesKernel->addArg(polarizability);
    mapTorqueKernel = program->createKernel("mapTorqueToForce");
    mapTorqueKernel->addArg(cc.getLongForceBuffer());
    mapTorqueKernel->addArg(torque);
    mapTorqueKernel->addArg(cc.getPosq());
    mapTorqueKernel->addArg(multipoleParticles);
    program = cc.compileProgram(CommonAmoebaKernelSources::multipoleInducedField, defines);
    initExtrapolatedKernel = program->createKernel("initExtrapolatedDipoles");
    initExtrapolatedKernel->addArg(inducedDipole);
    initExtrapolatedKernel->addArg(extrapolatedDipole);
    iterateExtrapolatedKernel = program->createKernel("iterateExtrapolatedDipoles");
    iterateExtrapolatedKernel->addArg();
    iterateExtrapolatedKernel->addArg(inducedDipole);
    iterateExtrapolatedKernel->addArg(extrapolatedDipole);
    iterateExtrapolatedKernel->addArg(inducedField);
    iterateExtrapolatedKernel->addArg(polarizability);
    computeExtrapolatedKernel = program->createKernel("computeExtrapolatedDipoles");
    computeExtrapolatedKernel->addArg(inducedDipole);
    computeExtrapolatedKernel->addArg(extrapolatedDipole);
    polarizationEnergyKernel = program->createKernel("computePolarizationEnergy");
    polarizationEnergyKernel->addArg(cc.getEnergyBuffer());
    polarizationEnergyKernel->addArg(inducedDipole);
    polarizationEnergyKernel->addArg(extrapolatedDipole);
    polarizationEnergyKernel->addArg(polarizability);

    // Set up PME.
    
    if (usePME) {

        // Create required data structures.

        int gridElements = gridSizeX*gridSizeY*gridSizeZ;
        gridElements = max(gridElements, dispersionGridSizeX*dispersionGridSizeY*dispersionGridSizeZ);
        pmeGrid1.initialize(cc, gridElements, elementSize, "pmeGrid1");
        pmeGrid2.initialize(cc, gridElements, 2*elementSize, "pmeGrid2");
        if (useFixedPointChargeSpreading()) {
            pmeGridLong.initialize(cc, 2*gridElements, sizeof(long long), "pmeGridLong");
            cc.addAutoclearBuffer(pmeGridLong);
        }
        else
            cc.addAutoclearBuffer(pmeGrid1);
        pmeBsplineModuliX.initialize(cc, gridSizeX, elementSize, "pmeBsplineModuliX");
        pmeBsplineModuliY.initialize(cc, gridSizeY, elementSize, "pmeBsplineModuliY");
        pmeBsplineModuliZ.initialize(cc, gridSizeZ, elementSize, "pmeBsplineModuliZ");
        dpmeBsplineModuliX.initialize(cc, dispersionGridSizeX, elementSize, "dpmeBsplineModuliX");
        dpmeBsplineModuliY.initialize(cc, dispersionGridSizeY, elementSize, "dpmeBsplineModuliY");
        dpmeBsplineModuliZ.initialize(cc, dispersionGridSizeZ, elementSize, "dpmeBsplineModuliZ");
        pmePhi.initialize(cc, 20*numParticles, elementSize, "pmePhi");
        pmePhidp.initialize(cc, 20*numParticles, elementSize, "pmePhidp");
        pmeCphi.initialize(cc, 10*numParticles, elementSize, "pmeCphi");
        pmeAtomGridIndex.initialize<mm_int2>(cc, numParticles, "pmeAtomGridIndex");

        // Create the PME kernels.

        map<string, string> pmeDefines;
        pmeDefines["HIPPO"] = "1";
        pmeDefines["EWALD_ALPHA"] = cc.doubleToString(pmeAlpha);
        pmeDefines["DISPERSION_EWALD_ALPHA"] = cc.doubleToString(dpmeAlpha);
        pmeDefines["PME_ORDER"] = cc.intToString(PmeOrder);
        pmeDefines["NUM_ATOMS"] = cc.intToString(numParticles);
        pmeDefines["PADDED_NUM_ATOMS"] = cc.intToString(cc.getPaddedNumAtoms());
        pmeDefines["EPSILON_FACTOR"] = cc.doubleToString(ONE_4PI_EPS0);
        pmeDefines["GRID_SIZE_X"] = cc.intToString(gridSizeX);
        pmeDefines["GRID_SIZE_Y"] = cc.intToString(gridSizeY);
        pmeDefines["GRID_SIZE_Z"] = cc.intToString(gridSizeZ);
        pmeDefines["M_PI"] = cc.doubleToString(M_PI);
        pmeDefines["SQRT_PI"] = cc.doubleToString(sqrt(M_PI));
        pmeDefines["EXTRAPOLATION_COEFFICIENTS_SUM"] = coefficients.str();
        pmeDefines["MAX_EXTRAPOLATION_ORDER"] = cc.intToString(maxExtrapolationOrder);
        if (useFixedPointChargeSpreading())
            pmeDefines["USE_FIXED_POINT_CHARGE_SPREADING"] = "";
        program = cc.compileProgram(CommonAmoebaKernelSources::multipolePme, pmeDefines);
        pmeTransformMultipolesKernel = program->createKernel("transformMultipolesToFractionalCoordinates");
        pmeTransformMultipolesKernel->addArg(labDipoles);
        pmeTransformMultipolesKernel->addArg(labQuadrupoles[0]);
        pmeTransformMultipolesKernel->addArg(labQuadrupoles[1]);
        pmeTransformMultipolesKernel->addArg(labQuadrupoles[2]);
        pmeTransformMultipolesKernel->addArg(labQuadrupoles[3]);
        pmeTransformMultipolesKernel->addArg(labQuadrupoles[4]);
        pmeTransformMultipolesKernel->addArg(fracDipoles);
        pmeTransformMultipolesKernel->addArg(fracQuadrupoles);
        for (int i = 0; i < 3; i++)
            pmeTransformMultipolesKernel->addArg();
        pmeTransformPotentialKernel = program->createKernel("transformPotentialToCartesianCoordinates");
        pmeTransformPotentialKernel->addArg();
        pmeTransformPotentialKernel->addArg(pmeCphi);
        for (int i = 0; i < 3; i++)
            pmeTransformPotentialKernel->addArg();
        pmeSpreadFixedMultipolesKernel = program->createKernel("gridSpreadFixedMultipoles");
        pmeSpreadFixedMultipolesKernel->addArg(cc.getPosq());
        pmeSpreadFixedMultipolesKernel->addArg(fracDipoles);
        pmeSpreadFixedMultipolesKernel->addArg(fracQuadrupoles);
        if (useFixedPointChargeSpreading())
            pmeSpreadFixedMultipolesKernel->addArg(pmeGridLong);
        else
            pmeSpreadFixedMultipolesKernel->addArg(pmeGrid1);
        pmeSpreadFixedMultipolesKernel->addArg(coreCharge);
        pmeSpreadFixedMultipolesKernel->addArg(valenceCharge);
        for (int i = 0; i < 6; i++)
            pmeSpreadFixedMultipolesKernel->addArg();
        pmeSpreadInducedDipolesKernel = program->createKernel("gridSpreadInducedDipoles");
        pmeSpreadInducedDipolesKernel->addArg(cc.getPosq());
        pmeSpreadInducedDipolesKernel->addArg(inducedDipole);
        if (useFixedPointChargeSpreading())
            pmeSpreadInducedDipolesKernel->addArg(pmeGridLong);
        else
            pmeSpreadInducedDipolesKernel->addArg(pmeGrid1);
        for (int i = 0; i < 6; i++)
            pmeSpreadInducedDipolesKernel->addArg();
        if (useFixedPointChargeSpreading()) {
            pmeFinishSpreadChargeKernel = program->createKernel("finishSpreadCharge");
            pmeFinishSpreadChargeKernel->addArg(pmeGridLong);
            pmeFinishSpreadChargeKernel->addArg(pmeGrid1);
        }
        pmeConvolutionKernel = program->createKernel("reciprocalConvolution");
        pmeConvolutionKernel->addArg(pmeGrid2);
        pmeConvolutionKernel->addArg(pmeBsplineModuliX);
        pmeConvolutionKernel->addArg(pmeBsplineModuliY);
        pmeConvolutionKernel->addArg(pmeBsplineModuliX);
        for (int i = 0; i < 4; i++)
            pmeConvolutionKernel->addArg();
        pmeFixedPotentialKernel = program->createKernel("computeFixedPotentialFromGrid");
        pmeFixedPotentialKernel->addArg(pmeGrid1);
        pmeFixedPotentialKernel->addArg(pmePhi);
        pmeFixedPotentialKernel->addArg(field);
        pmeFixedPotentialKernel->addArg(cc.getPosq());
        pmeFixedPotentialKernel->addArg(labDipoles);
        for (int i = 0; i < 6; i++)
            pmeFixedPotentialKernel->addArg();
        pmeInducedPotentialKernel = program->createKernel("computeInducedPotentialFromGrid");
        pmeInducedPotentialKernel->addArg(pmeGrid1);
        pmeInducedPotentialKernel->addArg(extrapolatedPhi);
        pmeInducedPotentialKernel->addArg();
        pmeInducedPotentialKernel->addArg(pmePhidp);
        pmeInducedPotentialKernel->addArg(cc.getPosq());
        for (int i = 0; i < 6; i++)
            pmeInducedPotentialKernel->addArg();
        pmeFixedForceKernel = program->createKernel("computeFixedMultipoleForceAndEnergy");
        pmeFixedForceKernel->addArg(cc.getPosq());
        pmeFixedForceKernel->addArg(cc.getLongForceBuffer());
        pmeFixedForceKernel->addArg(torque);
        pmeFixedForceKernel->addArg(cc.getEnergyBuffer());
        pmeFixedForceKernel->addArg(labDipoles);
        pmeFixedForceKernel->addArg(coreCharge);
        pmeFixedForceKernel->addArg(valenceCharge);
        pmeFixedForceKernel->addArg(labQuadrupoles[0]);
        pmeFixedForceKernel->addArg(labQuadrupoles[1]);
        pmeFixedForceKernel->addArg(labQuadrupoles[2]);
        pmeFixedForceKernel->addArg(labQuadrupoles[3]);
        pmeFixedForceKernel->addArg(labQuadrupoles[4]);
        pmeFixedForceKernel->addArg(fracDipoles);
        pmeFixedForceKernel->addArg(fracQuadrupoles);
        pmeFixedForceKernel->addArg(pmePhi);
        pmeFixedForceKernel->addArg(pmeCphi);
        for (int i = 0; i < 3; i++)
            pmeFixedForceKernel->addArg();
        pmeInducedForceKernel = program->createKernel("computeInducedDipoleForceAndEnergy");
        pmeInducedForceKernel->addArg(cc.getPosq());
        pmeInducedForceKernel->addArg(cc.getLongForceBuffer());
        pmeInducedForceKernel->addArg(torque);
        pmeInducedForceKernel->addArg(cc.getEnergyBuffer());
        pmeInducedForceKernel->addArg(labDipoles);
        pmeInducedForceKernel->addArg(coreCharge);
        pmeInducedForceKernel->addArg(valenceCharge);
        pmeInducedForceKernel->addArg(extrapolatedDipole);
        pmeInducedForceKernel->addArg(extrapolatedPhi);
        pmeInducedForceKernel->addArg(labQuadrupoles[0]);
        pmeInducedForceKernel->addArg(labQuadrupoles[1]);
        pmeInducedForceKernel->addArg(labQuadrupoles[2]);
        pmeInducedForceKernel->addArg(labQuadrupoles[3]);
        pmeInducedForceKernel->addArg(labQuadrupoles[4]);
        pmeInducedForceKernel->addArg(fracDipoles);
        pmeInducedForceKernel->addArg(fracQuadrupoles);
        pmeInducedForceKernel->addArg(inducedDipole);
        pmeInducedForceKernel->addArg(pmePhi);
        pmeInducedForceKernel->addArg(pmePhidp);
        pmeInducedForceKernel->addArg(pmeCphi);
        for (int i = 0; i < 3; i++)
            pmeInducedForceKernel->addArg();
        pmeRecordInducedFieldDipolesKernel = program->createKernel("recordInducedFieldDipoles");
        pmeRecordInducedFieldDipolesKernel->addArg(pmePhidp);
        pmeRecordInducedFieldDipolesKernel->addArg(inducedField);
        pmeRecordInducedFieldDipolesKernel->addArg(inducedDipole);
        for (int i = 0; i < 3; i++)
            pmeRecordInducedFieldDipolesKernel->addArg();
        pmeSelfEnergyKernel = program->createKernel("calculateSelfEnergyAndTorque");
        pmeSelfEnergyKernel->addArg(torque);
        pmeSelfEnergyKernel->addArg(cc.getEnergyBuffer());
        pmeSelfEnergyKernel->addArg(labDipoles);
        pmeSelfEnergyKernel->addArg(coreCharge);
        pmeSelfEnergyKernel->addArg(valenceCharge);
        pmeSelfEnergyKernel->addArg(c6);
        pmeSelfEnergyKernel->addArg(inducedDipole);
        pmeSelfEnergyKernel->addArg(labQuadrupoles[0]);
        pmeSelfEnergyKernel->addArg(labQuadrupoles[1]);
        pmeSelfEnergyKernel->addArg(labQuadrupoles[2]);
        pmeSelfEnergyKernel->addArg(labQuadrupoles[3]);
        pmeSelfEnergyKernel->addArg(labQuadrupoles[4]);

        // Create the dispersion PME kernels.

        pmeDefines["EWALD_ALPHA"] = cc.doubleToString(dpmeAlpha);
        pmeDefines["EPSILON_FACTOR"] = "1";
        pmeDefines["GRID_SIZE_X"] = cc.intToString(dispersionGridSizeX);
        pmeDefines["GRID_SIZE_Y"] = cc.intToString(dispersionGridSizeY);
        pmeDefines["GRID_SIZE_Z"] = cc.intToString(dispersionGridSizeZ);
        pmeDefines["RECIP_EXP_FACTOR"] = cc.doubleToString(M_PI*M_PI/(dpmeAlpha*dpmeAlpha));
        pmeDefines["CHARGE"] = "charges[atom]";
        pmeDefines["USE_LJPME"] = "1";
        program = cc.compileProgram(CommonKernelSources::pme, pmeDefines);
        if (useFixedPointChargeSpreading()) {
            dpmeFinishSpreadChargeKernel = program->createKernel("finishSpreadCharge");
            dpmeFinishSpreadChargeKernel->addArg(pmeGrid2);
            dpmeFinishSpreadChargeKernel->addArg(pmeGrid1);
        }
        dpmeGridIndexKernel = program->createKernel("findAtomGridIndex");
        dpmeGridIndexKernel->addArg(cc.getPosq());
        dpmeGridIndexKernel->addArg(pmeAtomGridIndex);
        for (int i = 0; i < 8; i++)
            dpmeGridIndexKernel->addArg();
        dpmeSpreadChargeKernel = program->createKernel("gridSpreadCharge");
        dpmeSpreadChargeKernel->addArg(cc.getPosq());
        if (useFixedPointChargeSpreading())
            dpmeSpreadChargeKernel->addArg(pmeGrid2);
        else
            dpmeSpreadChargeKernel->addArg(pmeGrid1);
        for (int i = 0; i < 8; i++)
            dpmeSpreadChargeKernel->addArg();
        dpmeSpreadChargeKernel->addArg(pmeAtomGridIndex);
        dpmeSpreadChargeKernel->addArg(c6);
        dpmeConvolutionKernel = program->createKernel("reciprocalConvolution");
        dpmeConvolutionKernel->addArg(pmeGrid2);
        dpmeConvolutionKernel->addArg(dpmeBsplineModuliX);
        dpmeConvolutionKernel->addArg(dpmeBsplineModuliY);
        dpmeConvolutionKernel->addArg(dpmeBsplineModuliZ);
        for (int i = 0; i < 3; i++)
            dpmeConvolutionKernel->addArg();
        dpmeEvalEnergyKernel = program->createKernel("gridEvaluateEnergy");
        dpmeEvalEnergyKernel->addArg(pmeGrid2);
        dpmeEvalEnergyKernel->addArg(cc.getEnergyBuffer());
        dpmeEvalEnergyKernel->addArg(dpmeBsplineModuliX);
        dpmeEvalEnergyKernel->addArg(dpmeBsplineModuliY);
        dpmeEvalEnergyKernel->addArg(dpmeBsplineModuliZ);
        for (int i = 0; i < 3; i++)
            dpmeEvalEnergyKernel->addArg();
        dpmeInterpolateForceKernel = program->createKernel("gridInterpolateForce");
        dpmeInterpolateForceKernel->addArg(cc.getPosq());
        dpmeInterpolateForceKernel->addArg(cc.getLongForceBuffer());
        dpmeInterpolateForceKernel->addArg(pmeGrid1);
        for (int i = 0; i < 8; i++)
            dpmeInterpolateForceKernel->addArg();
        dpmeInterpolateForceKernel->addArg(pmeAtomGridIndex);
        dpmeInterpolateForceKernel->addArg(c6);

        // Initialize the B-spline moduli.

        double data[PmeOrder];
        double x = 0.0;
        data[0] = 1.0 - x;
        data[1] = x;
        for (int i = 2; i < PmeOrder; i++) {
            double denom = 1.0/i;
            data[i] = x*data[i-1]*denom;
            for (int j = 1; j < i; j++)
                data[i-j] = ((x+j)*data[i-j-1] + ((i-j+1)-x)*data[i-j])*denom;
            data[0] = (1.0-x)*data[0]*denom;
        }
        int maxSize = max(max(gridSizeX, gridSizeY), gridSizeZ);
        vector<double> bsplines_data(maxSize+1, 0.0);
        for (int i = 2; i <= PmeOrder+1; i++)
            bsplines_data[i] = data[i-2];
        for (int dim = 0; dim < 3; dim++) {
            int ndata = (dim == 0 ? gridSizeX : dim == 1 ? gridSizeY : gridSizeZ);
            vector<double> moduli(ndata);

            // get the modulus of the discrete Fourier transform

            double factor = 2.0*M_PI/ndata;
            for (int i = 0; i < ndata; i++) {
                double sc = 0.0;
                double ss = 0.0;
                for (int j = 1; j <= ndata; j++) {
                    double arg = factor*i*(j-1);
                    sc += bsplines_data[j]*cos(arg);
                    ss += bsplines_data[j]*sin(arg);
                }
                moduli[i] = sc*sc+ss*ss;
            }

            // Fix for exponential Euler spline interpolation failure.

            double eps = 1.0e-7;
            if (moduli[0] < eps)
                moduli[0] = 0.9*moduli[1];
            for (int i = 1; i < ndata-1; i++)
                if (moduli[i] < eps)
                    moduli[i] = 0.9*(moduli[i-1]+moduli[i+1]);
            if (moduli[ndata-1] < eps)
                moduli[ndata-1] = 0.9*moduli[ndata-2];

            // Compute and apply the optimal zeta coefficient.

            int jcut = 50;
            for (int i = 1; i <= ndata; i++) {
                int k = i - 1;
                if (i > ndata/2)
                    k = k - ndata;
                double zeta;
                if (k == 0)
                    zeta = 1.0;
                else {
                    double sum1 = 1.0;
                    double sum2 = 1.0;
                    factor = M_PI*k/ndata;
                    for (int j = 1; j <= jcut; j++) {
                        double arg = factor/(factor+M_PI*j);
                        sum1 += pow(arg, PmeOrder);
                        sum2 += pow(arg, 2*PmeOrder);
                    }
                    for (int j = 1; j <= jcut; j++) {
                        double arg = factor/(factor-M_PI*j);
                        sum1 += pow(arg, PmeOrder);
                        sum2 += pow(arg, 2*PmeOrder);
                    }
                    zeta = sum2/sum1;
                }
                moduli[i-1] = moduli[i-1]*zeta*zeta;
            }
            if (cc.getUseDoublePrecision()) {
                if (dim == 0)
                    pmeBsplineModuliX.upload(moduli);
                else if (dim == 1)
                    pmeBsplineModuliY.upload(moduli);
                else
                    pmeBsplineModuliZ.upload(moduli);
            }
            else {
                vector<float> modulif(ndata);
                for (int i = 0; i < ndata; i++)
                    modulif[i] = (float) moduli[i];
                if (dim == 0)
                    pmeBsplineModuliX.upload(modulif);
                else if (dim == 1)
                    pmeBsplineModuliY.upload(modulif);
                else
                    pmeBsplineModuliZ.upload(modulif);
            }
        }

        // Initialize the b-spline moduli for dispersion PME.

        maxSize = max(max(dispersionGridSizeX, dispersionGridSizeY), dispersionGridSizeZ);
        vector<double> ddata(PmeOrder);
        bsplines_data.resize(maxSize);
        data[PmeOrder-1] = 0.0;
        data[1] = 0.0;
        data[0] = 1.0;
        for (int i = 3; i < PmeOrder; i++) {
            double div = 1.0/(i-1.0);
            data[i-1] = 0.0;
            for (int j = 1; j < (i-1); j++)
                data[i-j-1] = div*(j*data[i-j-2]+(i-j)*data[i-j-1]);
            data[0] = div*data[0];
        }

        // Differentiate.

        ddata[0] = -data[0];
        for (int i = 1; i < PmeOrder; i++)
            ddata[i] = data[i-1]-data[i];
        double div = 1.0/(PmeOrder-1);
        data[PmeOrder-1] = 0.0;
        for (int i = 1; i < (PmeOrder-1); i++)
            data[PmeOrder-i-1] = div*(i*data[PmeOrder-i-2]+(PmeOrder-i)*data[PmeOrder-i-1]);
        data[0] = div*data[0];
        for (int i = 0; i < maxSize; i++)
            bsplines_data[i] = 0.0;
        for (int i = 1; i <= PmeOrder; i++)
            bsplines_data[i] = data[i-1];

        // Evaluate the actual bspline moduli for X/Y/Z.

        for(int dim = 0; dim < 3; dim++) {
            int ndata = (dim == 0 ? dispersionGridSizeX : dim == 1 ? dispersionGridSizeY : dispersionGridSizeZ);
            vector<double> moduli(ndata);
            for (int i = 0; i < ndata; i++) {
                double sc = 0.0;
                double ss = 0.0;
                for (int j = 0; j < ndata; j++) {
                    double arg = (2.0*M_PI*i*j)/ndata;
                    sc += bsplines_data[j]*cos(arg);
                    ss += bsplines_data[j]*sin(arg);
                }
                moduli[i] = sc*sc+ss*ss;
            }
            for (int i = 0; i < ndata; i++)
                if (moduli[i] < 1.0e-7)
                    moduli[i] = (moduli[i-1]+moduli[i+1])*0.5;
            if (dim == 0)
                dpmeBsplineModuliX.upload(moduli, true);
            else if (dim == 1)
                dpmeBsplineModuliY.upload(moduli, true);
            else
                dpmeBsplineModuliZ.upload(moduli, true);
        }
    }

    // Add the interaction to the default nonbonded kernel.
    
    NonbondedUtilities& nb = cc.getNonbondedUtilities();
    nb.setKernelSource(CommonAmoebaKernelSources::hippoInteractionHeader+CommonAmoebaKernelSources::hippoNonbonded);
    nb.addArgument(ComputeParameterInfo(torque, "torqueBuffers", "mm_ulong", 1, false));
    nb.addArgument(ComputeParameterInfo(extrapolatedDipole, "extrapolatedDipole", "real", 1));
    nb.addParameter(ComputeParameterInfo(coreCharge, "coreCharge", "real", 1));
    nb.addParameter(ComputeParameterInfo(valenceCharge, "valenceCharge", "real", 1));
    nb.addParameter(ComputeParameterInfo(alpha, "alpha", "real", 1));
    nb.addParameter(ComputeParameterInfo(epsilon, "epsilon", "real", 1));
    nb.addParameter(ComputeParameterInfo(damping, "damping", "real", 1));
    nb.addParameter(ComputeParameterInfo(c6, "c6", "real", 1));
    nb.addParameter(ComputeParameterInfo(pauliK, "pauliK", "real", 1));
    nb.addParameter(ComputeParameterInfo(pauliQ, "pauliQ", "real", 1));
    nb.addParameter(ComputeParameterInfo(pauliAlpha, "pauliAlpha", "real", 1));
    nb.addParameter(ComputeParameterInfo(labDipoles,"dipole", "real", 3));
    nb.addParameter(ComputeParameterInfo(inducedDipole, "inducedDipole", "real", 3));
    nb.addParameter(ComputeParameterInfo(labQuadrupoles[0], "qXX", "real", 1));
    nb.addParameter(ComputeParameterInfo(labQuadrupoles[1], "qXY", "real", 1));
    nb.addParameter(ComputeParameterInfo(labQuadrupoles[2], "qXZ", "real", 1));
    nb.addParameter(ComputeParameterInfo(labQuadrupoles[3], "qYY", "real", 1));
    nb.addParameter(ComputeParameterInfo(labQuadrupoles[4], "qYZ", "real", 1));
    map<string, string> replacements;
    replacements["ENERGY_SCALE_FACTOR"] = cc.doubleToString(ONE_4PI_EPS0);
    replacements["SWITCH_CUTOFF"] = cc.doubleToString(force.getSwitchingDistance());
    replacements["SWITCH_C3"] = cc.doubleToString(10/pow(force.getSwitchingDistance()-force.getCutoffDistance(), 3.0));
    replacements["SWITCH_C4"] = cc.doubleToString(15/pow(force.getSwitchingDistance()-force.getCutoffDistance(), 4.0));
    replacements["SWITCH_C5"] = cc.doubleToString(6/pow(force.getSwitchingDistance()-force.getCutoffDistance(), 5.0));
    replacements["MAX_EXTRAPOLATION_ORDER"] = cc.intToString(maxExtrapolationOrder);
    replacements["EXTRAPOLATION_COEFFICIENTS_SUM"] = coefficients.str();
    replacements["USE_EWALD"] = (usePME ? "1" : "0");
    replacements["PME_ALPHA"] = (usePME ? cc.doubleToString(pmeAlpha) : "0");
    replacements["DPME_ALPHA"] = (usePME ? cc.doubleToString(dpmeAlpha) : "0");
    replacements["SQRT_PI"] = cc.doubleToString(sqrt(M_PI));
    string interactionSource = cc.replaceStrings(CommonAmoebaKernelSources::hippoInteraction, replacements);
    nb.addInteraction(usePME, usePME, true, force.getCutoffDistance(), exclusions, interactionSource, force.getForceGroup());
    nb.setUsePadding(false);
    
    // Create the kernel for computing exceptions.
    
    if (exceptionAtoms.isInitialized()) {
        replacements["COMPUTE_INTERACTION"] = interactionSource;
        string exceptionsSrc = CommonAmoebaKernelSources::hippoInteractionHeader+CommonAmoebaKernelSources::hippoNonbondedExceptions;
        exceptionsSrc = cc.replaceStrings(exceptionsSrc, replacements);
        defines["NUM_EXCEPTIONS"] = cc.intToString(exceptionAtoms.getSize());
        program = cc.compileProgram(exceptionsSrc, defines);
        computeExceptionsKernel = program->createKernel("computeNonbondedExceptions");
        computeExceptionsKernel->addArg(cc.getLongForceBuffer());
        computeExceptionsKernel->addArg(cc.getEnergyBuffer());
        computeExceptionsKernel->addArg(torque);
        computeExceptionsKernel->addArg(cc.getPosq());
        computeExceptionsKernel->addArg(exceptionAtoms);
        computeExceptionsKernel->addArg(exceptionScales[0]);
        computeExceptionsKernel->addArg(exceptionScales[1]);
        computeExceptionsKernel->addArg(exceptionScales[2]);
        computeExceptionsKernel->addArg(exceptionScales[3]);
        computeExceptionsKernel->addArg(exceptionScales[4]);
        computeExceptionsKernel->addArg(exceptionScales[5]);
        computeExceptionsKernel->addArg(coreCharge);
        computeExceptionsKernel->addArg(valenceCharge);
        computeExceptionsKernel->addArg(alpha);
        computeExceptionsKernel->addArg(epsilon);
        computeExceptionsKernel->addArg(damping);
        computeExceptionsKernel->addArg(c6);
        computeExceptionsKernel->addArg(pauliK);
        computeExceptionsKernel->addArg(pauliQ);
        computeExceptionsKernel->addArg(pauliAlpha);
        computeExceptionsKernel->addArg(labDipoles);
        computeExceptionsKernel->addArg(inducedDipole);
        computeExceptionsKernel->addArg(labQuadrupoles[0]);
        computeExceptionsKernel->addArg(labQuadrupoles[1]);
        computeExceptionsKernel->addArg(labQuadrupoles[2]);
        computeExceptionsKernel->addArg(labQuadrupoles[3]);
        computeExceptionsKernel->addArg(labQuadrupoles[4]);
        computeExceptionsKernel->addArg(extrapolatedDipole);
        if (nb.getUseCutoff())
            for (int i = 0; i < 5; i++)
                computeExceptionsKernel->addArg();
    }
    cc.addForce(new ForceInfo(force));
    cc.addPostComputation(new TorquePostComputation(*this));
    fieldThreadBlockSize = max(32, cc.getNonbondedUtilities().getForceThreadBlockSize());
}

void CommonCalcHippoNonbondedForceKernel::createFieldKernel(const string& interactionSrc, vector<ComputeArray*> params,
            ComputeArray& fieldBuffer, ComputeKernel& kernel, ComputeKernel& exceptionKernel, ComputeArray& exceptionScale) {
    // Create the kernel source.

    map<string, string> replacements;
    replacements["COMPUTE_FIELD"] = interactionSrc;
    stringstream extraArgs, atomParams, loadLocal1, loadLocal2, load1, load2, load3;
    for (auto param : params) {
        string name = param->getName();
        bool isReal3 = (param->getSize() == cc.getPaddedNumAtoms()*3);
        string type = (isReal3 ? "real3" : "real");
        extraArgs << ", GLOBAL const real* RESTRICT " << name;
        atomParams << type << " " << name << ";\n";
        loadLocal1 << "localData[localAtomIndex]." << name << " = " << name << "1;\n";
        if (isReal3) {
            loadLocal2 << "localData[localAtomIndex]." << name << " = make_real3(" << name << "[3*j], " << name << "[3*j+1], " << name << "[3*j+2]);\n";
            load1 << type << " " << name << "1 = make_real3(" << name << "[3*atom1], " << name << "[3*atom1+1], " << name << "[3*atom1+2]);\n";
            load3 << type << " " << name << "2 = make_real3(" << name << "[3*atom2], " << name << "[3*atom2+1], " << name << "[3*atom2+2]);\n";
        }
        else {
            loadLocal2 << "localData[localAtomIndex]." << name << " = " << name << "[j];\n";
            load1 << type << " " << name << "1 = " << name << "[atom1];\n";
            load3 << type << " " << name << "2 = " << name << "[atom2];\n";
        }
        load2 << type << " " << name << "2 = localData[atom2]." << name << ";\n";
    }
    replacements["PARAMETER_ARGUMENTS"] = extraArgs.str();
    replacements["ATOM_PARAMETER_DATA"] = atomParams.str();
    replacements["LOAD_LOCAL_PARAMETERS_FROM_1"] = loadLocal1.str();
    replacements["LOAD_LOCAL_PARAMETERS_FROM_GLOBAL"] = loadLocal2.str();
    replacements["LOAD_ATOM1_PARAMETERS"] = load1.str();
    replacements["LOAD_ATOM2_PARAMETERS"] = load2.str();
    replacements["LOAD_ATOM2_PARAMETERS_FROM_GLOBAL"] = load3.str();
    string src = cc.replaceStrings(CommonAmoebaKernelSources::hippoComputeField, replacements);

    // Set defines and create the kernel.

    map<string, string> defines;
    if (usePME) {
        defines["USE_CUTOFF"] = "1";
        defines["USE_PERIODIC"] = "1";
        defines["USE_EWALD"] = "1";
        defines["PME_ALPHA"] = cc.doubleToString(pmeAlpha);
        defines["SQRT_PI"] = cc.doubleToString(sqrt(M_PI));
    }
    defines["WARPS_PER_GROUP"] = cc.intToString(fieldThreadBlockSize/ComputeContext::TileSize);
    defines["THREAD_BLOCK_SIZE"] = cc.intToString(fieldThreadBlockSize);
    defines["CUTOFF"] = cc.doubleToString(cutoff);
    defines["CUTOFF_SQUARED"] = cc.doubleToString(cutoff*cutoff);
    defines["NUM_ATOMS"] = cc.intToString(cc.getNumAtoms());
    defines["PADDED_NUM_ATOMS"] = cc.intToString(cc.getPaddedNumAtoms());
    defines["NUM_BLOCKS"] = cc.intToString(cc.getNumAtomBlocks());
    defines["TILE_SIZE"] = cc.intToString(ComputeContext::TileSize);
    defines["NUM_TILES_WITH_EXCLUSIONS"] = cc.intToString(cc.getNonbondedUtilities().getExclusionTiles().getSize());
    defines["NUM_EXCEPTIONS"] = cc.intToString(exceptionAtoms.isInitialized() ? exceptionAtoms.getSize() : 0);
    ComputeProgram program = cc.compileProgram(src, defines);
    kernel = program->createKernel("computeField");

    // Build the list of arguments.

    NonbondedUtilities& nb = cc.getNonbondedUtilities();
    kernel->addArg(cc.getPosq());
    kernel->addArg(cc.getNonbondedUtilities().getExclusions());
    kernel->addArg(cc.getNonbondedUtilities().getExclusionTiles());
    kernel->addArg(fieldBuffer);
    if (nb.getUseCutoff()) {
        kernel->addArg(nb.getInteractingTiles());
        kernel->addArg(nb.getInteractionCount());
        for (int i = 0; i < 5; i++)
            kernel->addArg();
        kernel->addArg(maxTiles);
        kernel->addArg(nb.getBlockCenters());
        kernel->addArg(nb.getBlockBoundingBoxes());
        kernel->addArg(nb.getInteractingAtoms());
    }
    else
        kernel->addArg(maxTiles);
    for (auto param : params)
        kernel->addArg(*param);
    
    // If there are any exceptions, build the kernel and arguments to compute them.
    
    if (exceptionAtoms.isInitialized()) {
        exceptionKernel = program->createKernel("computeFieldExceptions");
        exceptionKernel->addArg(cc.getPosq());
        exceptionKernel->addArg(fieldBuffer);
        exceptionKernel->addArg(exceptionAtoms);
        exceptionKernel->addArg(exceptionScale);
        if (nb.getUseCutoff())
            for (int i = 0; i < 5; i++)
                exceptionKernel->addArg();
        for (auto param : params)
            exceptionKernel->addArg(*param);
    }
}

double CommonCalcHippoNonbondedForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    ContextSelector selector(cc);
    NonbondedUtilities& nb = cc.getNonbondedUtilities();
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
        
        // These kernels can't be compiled in initialize(), because the nonbonded utilities object
        // has not yet been initialized then.

        maxTiles = (nb.getUseCutoff() ? nb.getInteractingTiles().getSize() : cc.getNumAtomBlocks()*(cc.getNumAtomBlocks()+1)/2);
        createFieldKernel(CommonAmoebaKernelSources::hippoFixedField, {&coreCharge, &valenceCharge, &alpha, &labDipoles, &labQuadrupoles[0],
                &labQuadrupoles[1], &labQuadrupoles[2], &labQuadrupoles[3], &labQuadrupoles[4]}, field, fixedFieldKernel,
                fixedFieldExceptionKernel, exceptionScales[1]);
        createFieldKernel(CommonAmoebaKernelSources::hippoMutualField, {&alpha, &inducedDipole}, inducedField, mutualFieldKernel,
                mutualFieldExceptionKernel, exceptionScales[2]);
    }

    // Compute the lab frame moments.

    computeMomentsKernel->execute(cc.getNumAtoms());

    if (usePME) {
        setPeriodicBoxArgs(cc, dpmeGridIndexKernel, 2);
        setPeriodicBoxArgs(cc, dpmeSpreadChargeKernel, 2);
        setPeriodicBoxArgs(cc, dpmeInterpolateForceKernel, 3);
        
        // Compute reciprocal box vectors.
        
        Vec3 a, b, c;
        cc.getPeriodicBoxVectors(a, b, c);
        double determinant = a[0]*b[1]*c[2];
        double scale = 1.0/determinant;
        mm_double4 recipBoxVectors[3];
        recipBoxVectors[0] = mm_double4(b[1]*c[2]*scale, 0, 0, 0);
        recipBoxVectors[1] = mm_double4(-b[0]*c[2]*scale, a[0]*c[2]*scale, 0, 0);
        recipBoxVectors[2] = mm_double4((b[0]*c[1]-b[1]*c[0])*scale, -a[0]*c[1]*scale, a[0]*b[1]*scale, 0);
        if (cc.getUseDoublePrecision()) {
            mm_double4 boxVectors[] = {mm_double4(a[0], a[1], a[2], 0), mm_double4(b[0], b[1], b[2], 0), mm_double4(c[0], c[1], c[2], 0)};
            pmeConvolutionKernel->setArg(4, mm_double4(a[0], b[1], c[2], 0));
            for (int i = 0; i < 3; i++) {
                pmeTransformMultipolesKernel->setArg(8+i, recipBoxVectors[i]);
                pmeTransformPotentialKernel->setArg(2+i, recipBoxVectors[i]);
                pmeSpreadFixedMultipolesKernel->setArg(6+i, boxVectors[i]);
                pmeSpreadFixedMultipolesKernel->setArg(9+i, recipBoxVectors[i]);
                pmeSpreadInducedDipolesKernel->setArg(3+i, boxVectors[i]);
                pmeSpreadInducedDipolesKernel->setArg(6+i, recipBoxVectors[i]);
                pmeConvolutionKernel->setArg(5+i, recipBoxVectors[i]);
                pmeFixedPotentialKernel->setArg(5+i, boxVectors[i]);
                pmeFixedPotentialKernel->setArg(8+i, recipBoxVectors[i]);
                pmeInducedPotentialKernel->setArg(5+i, boxVectors[i]);
                pmeInducedPotentialKernel->setArg(8+i, recipBoxVectors[i]);
                pmeFixedForceKernel->setArg(16+i, recipBoxVectors[i]);
                pmeInducedForceKernel->setArg(20+i, recipBoxVectors[i]);
                pmeRecordInducedFieldDipolesKernel->setArg(3+i, recipBoxVectors[i]);
                dpmeGridIndexKernel->setArg(7+i, recipBoxVectors[i]);
                dpmeSpreadChargeKernel->setArg(7+i, recipBoxVectors[i]);
                dpmeConvolutionKernel->setArg(4+i, recipBoxVectors[i]);
                dpmeEvalEnergyKernel->setArg(5+i, recipBoxVectors[i]);
                dpmeInterpolateForceKernel->setArg(8+i, recipBoxVectors[i]);
            }
        }
        else {
            mm_float4 boxVectors[] = {mm_float4(a[0], a[1], a[2], 0), mm_float4(b[0], b[1], b[2], 0), mm_float4(c[0], c[1], c[2], 0)};
            pmeConvolutionKernel->setArg(4, mm_float4(a[0], b[1], c[2], 0));
            mm_float4 recipBoxVectorsFloat[3];
            recipBoxVectorsFloat[0] = mm_float4((float) recipBoxVectors[0].x, 0, 0, 0);
            recipBoxVectorsFloat[1] = mm_float4((float) recipBoxVectors[1].x, (float) recipBoxVectors[1].y, 0, 0);
            recipBoxVectorsFloat[2] = mm_float4((float) recipBoxVectors[2].x, (float) recipBoxVectors[2].y, (float) recipBoxVectors[2].z, 0);
            for (int i = 0; i < 3; i++) {
                pmeTransformMultipolesKernel->setArg(8+i, recipBoxVectorsFloat[i]);
                pmeTransformPotentialKernel->setArg(2+i, recipBoxVectorsFloat[i]);
                pmeSpreadFixedMultipolesKernel->setArg(6+i, boxVectors[i]);
                pmeSpreadFixedMultipolesKernel->setArg(9+i, recipBoxVectorsFloat[i]);
                pmeSpreadInducedDipolesKernel->setArg(3+i, boxVectors[i]);
                pmeSpreadInducedDipolesKernel->setArg(6+i, recipBoxVectorsFloat[i]);
                pmeConvolutionKernel->setArg(5+i, recipBoxVectorsFloat[i]);
                pmeFixedPotentialKernel->setArg(5+i, boxVectors[i]);
                pmeFixedPotentialKernel->setArg(8+i, recipBoxVectorsFloat[i]);
                pmeInducedPotentialKernel->setArg(5+i, boxVectors[i]);
                pmeInducedPotentialKernel->setArg(8+i, recipBoxVectorsFloat[i]);
                pmeFixedForceKernel->setArg(16+i, recipBoxVectorsFloat[i]);
                pmeInducedForceKernel->setArg(20+i, recipBoxVectorsFloat[i]);
                pmeRecordInducedFieldDipolesKernel->setArg(3+i, recipBoxVectorsFloat[i]);
                dpmeGridIndexKernel->setArg(7+i, recipBoxVectorsFloat[i]);
                dpmeSpreadChargeKernel->setArg(7+i, recipBoxVectorsFloat[i]);
                dpmeConvolutionKernel->setArg(4+i, recipBoxVectorsFloat[i]);
                dpmeEvalEnergyKernel->setArg(5+i, recipBoxVectorsFloat[i]);
                dpmeInterpolateForceKernel->setArg(8+i, recipBoxVectorsFloat[i]);
            }
        }

        // Reciprocal space calculation for electrostatics.
        
        pmeTransformMultipolesKernel->execute(cc.getNumAtoms());
        pmeSpreadFixedMultipolesKernel->execute(cc.getNumAtoms());
        if (useFixedPointChargeSpreading())
            pmeFinishSpreadChargeKernel->execute(pmeGrid1.getSize());
        computeFFT(true, false);
        pmeConvolutionKernel->execute(gridSizeX*gridSizeY*gridSizeZ, 256);
        computeFFT(false, false);
        pmeFixedPotentialKernel->execute(cc.getNumAtoms());
        pmeTransformPotentialKernel->setArg(0, pmePhi);
        pmeTransformPotentialKernel->execute(cc.getNumAtoms());
        pmeFixedForceKernel->execute(cc.getNumAtoms());

        // Reciprocal space calculation for dispersion.

        dpmeGridIndexKernel->execute(cc.getNumAtoms());
        sortGridIndex();
        if (useFixedPointChargeSpreading())
            cc.clearBuffer(pmeGrid2);
        else
            cc.clearBuffer(pmeGrid1);
        dpmeSpreadChargeKernel->execute(PmeOrder*cc.getNumAtoms(), 128);
        if (useFixedPointChargeSpreading())
            dpmeFinishSpreadChargeKernel->execute(dispersionGridSizeX*dispersionGridSizeY*dispersionGridSizeZ, 256);
        computeFFT(true, true);
        if (includeEnergy)
            dpmeEvalEnergyKernel->execute(dispersionGridSizeX*dispersionGridSizeY*dispersionGridSizeZ);
        dpmeConvolutionKernel->execute(dispersionGridSizeX*dispersionGridSizeY*dispersionGridSizeZ, 256);
        computeFFT(false, true);
        dpmeInterpolateForceKernel->execute(cc.getNumAtoms(), 128);
    }

    // Compute the field from fixed multipoles.

    if (nb.getUseCutoff())
        setPeriodicBoxArgs(cc, fixedFieldKernel, 6);
    fixedFieldKernel->execute(nb.getNumForceThreadBlocks()*fieldThreadBlockSize, fieldThreadBlockSize);
    if (exceptionAtoms.isInitialized()) {
        if (nb.getUseCutoff())
            setPeriodicBoxArgs(cc, fixedFieldExceptionKernel, 4);
        fixedFieldExceptionKernel->execute(exceptionAtoms.getSize());
    }

    // Iterate the induced dipoles.

    computeExtrapolatedDipoles();

    // Add the polarization energy.

    if (includeEnergy)
        polarizationEnergyKernel->execute(cc.getNumAtoms());

    // Compute the forces due to the reciprocal space PME calculation for induced dipoles.

    if (usePME) {
        pmeTransformPotentialKernel->setArg(0, pmePhidp);
        pmeTransformPotentialKernel->execute(cc.getNumAtoms());
        pmeInducedForceKernel->execute(cc.getNumAtoms());
        pmeSelfEnergyKernel->execute(cc.getNumAtoms());
    }

    // Compute nonbonded exceptions.

    if (exceptionAtoms.isInitialized()) {
        if (nb.getUseCutoff())
            setPeriodicBoxArgs(cc, computeExceptionsKernel, 28);
        computeExceptionsKernel->execute(exceptionAtoms.getSize());
    }

    // Record the current atom positions so we can tell later if they have changed.
    
    cc.getPosq().copyTo(lastPositions);
    multipolesAreValid = true;
    return 0.0;
}

void CommonCalcHippoNonbondedForceKernel::computeInducedField(int optOrder) {
    NonbondedUtilities& nb = cc.getNonbondedUtilities();
    cc.clearBuffer(inducedField);
    if (nb.getUseCutoff())
        setPeriodicBoxArgs(cc, mutualFieldKernel, 6);
    mutualFieldKernel->execute(nb.getNumForceThreadBlocks()*fieldThreadBlockSize, fieldThreadBlockSize);
    if (exceptionAtoms.isInitialized()) {
        if (nb.getUseCutoff())
            setPeriodicBoxArgs(cc, mutualFieldExceptionKernel, 4);
        mutualFieldExceptionKernel->execute(exceptionAtoms.getSize());
    }
    if (usePME) {
        if (useFixedPointChargeSpreading())
            cc.clearBuffer(pmeGridLong);
        else
            cc.clearBuffer(pmeGrid1);
        pmeSpreadInducedDipolesKernel->execute(cc.getNumAtoms());
        if (useFixedPointChargeSpreading())
            pmeFinishSpreadChargeKernel->execute(pmeGrid1.getSize());
        computeFFT(true, false);
        pmeConvolutionKernel->execute(gridSizeX*gridSizeY*gridSizeZ, 256);
        computeFFT(false, false);
        pmeInducedPotentialKernel->setArg(2, optOrder);
        pmeInducedPotentialKernel->execute(cc.getNumAtoms());
        pmeRecordInducedFieldDipolesKernel->execute(cc.getNumAtoms());
    }
}

void CommonCalcHippoNonbondedForceKernel::computeExtrapolatedDipoles() {
    // Start by storing the direct dipoles as PT0

    recordInducedDipolesKernel->execute(cc.getNumAtoms());
    initExtrapolatedKernel->execute(extrapolatedDipole.getSize());

    // Recursively apply alpha.Tau to the µ_(n) components to generate µ_(n+1), and store the result

    for (int order = 1; order < maxExtrapolationOrder; ++order) {
        computeInducedField(order-1);
        iterateExtrapolatedKernel->setArg(0, order);
        iterateExtrapolatedKernel->execute(extrapolatedDipole.getSize());
    }
    
    // Take a linear combination of the µ_(n) components to form the total dipole

    computeExtrapolatedKernel->execute(extrapolatedDipole.getSize());
    computeInducedField(maxExtrapolationOrder-1);
}

void CommonCalcHippoNonbondedForceKernel::addTorquesToForces() {
    mapTorqueKernel->execute(cc.getNumAtoms());
}

void CommonCalcHippoNonbondedForceKernel::getInducedDipoles(ContextImpl& context, vector<Vec3>& dipoles) {
    ContextSelector selector(cc);
    ensureMultipolesValid(context);
    int numParticles = cc.getNumAtoms();
    dipoles.resize(numParticles);
    const vector<int>& order = cc.getAtomIndex();
    if (cc.getUseDoublePrecision()) {
        vector<double> d;
        inducedDipole.download(d);
        for (int i = 0; i < numParticles; i++)
            dipoles[order[i]] = Vec3(d[3*i], d[3*i+1], d[3*i+2]);
    }
    else {
        vector<float> d;
        inducedDipole.download(d);
        for (int i = 0; i < numParticles; i++)
            dipoles[order[i]] = Vec3(d[3*i], d[3*i+1], d[3*i+2]);
    }
}

void CommonCalcHippoNonbondedForceKernel::ensureMultipolesValid(ContextImpl& context) {
    if (multipolesAreValid) {
        int numParticles = cc.getNumAtoms();
        if (cc.getUseDoublePrecision()) {
            vector<mm_double4> pos1, pos2;
            cc.getPosq().download(pos1);
            lastPositions.download(pos2);
            for (int i = 0; i < numParticles; i++)
                if (pos1[i].x != pos2[i].x || pos1[i].y != pos2[i].y || pos1[i].z != pos2[i].z) {
                    multipolesAreValid = false;
                    break;
                }
        }
        else {
            vector<mm_float4> pos1, pos2;
            cc.getPosq().download(pos1);
            lastPositions.download(pos2);
            for (int i = 0; i < numParticles; i++)
                if (pos1[i].x != pos2[i].x || pos1[i].y != pos2[i].y || pos1[i].z != pos2[i].z) {
                    multipolesAreValid = false;
                    break;
                }
        }
    }
    if (!multipolesAreValid)
        context.calcForcesAndEnergy(false, false, context.getIntegrator().getIntegrationForceGroups());
}

void CommonCalcHippoNonbondedForceKernel::getLabFramePermanentDipoles(ContextImpl& context, vector<Vec3>& dipoles) {
    ContextSelector selector(cc);
    ensureMultipolesValid(context);
    int numParticles = cc.getNumAtoms();
    dipoles.resize(numParticles);
    const vector<int>& order = cc.getAtomIndex();
    if (cc.getUseDoublePrecision()) {
        vector<double> labDipoleVec;
        labDipoles.download(labDipoleVec);
        for (int i = 0; i < numParticles; i++)
            dipoles[order[i]] = Vec3(labDipoleVec[3*i], labDipoleVec[3*i+1], labDipoleVec[3*i+2]);
    }
    else {
        vector<float> labDipoleVec;
        labDipoles.download(labDipoleVec);
        for (int i = 0; i < numParticles; i++)
            dipoles[order[i]] = Vec3(labDipoleVec[3*i], labDipoleVec[3*i+1], labDipoleVec[3*i+2]);
    }
}

void CommonCalcHippoNonbondedForceKernel::copyParametersToContext(ContextImpl& context, const HippoNonbondedForce& force) {
    // Make sure the new parameters are acceptable.
    
    ContextSelector selector(cc);
    if (force.getNumParticles() != cc.getNumAtoms())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");
    
    // Record the per-particle parameters.
    
    vector<double> coreChargeVec, valenceChargeVec, alphaVec, epsilonVec, dampingVec, c6Vec, pauliKVec, pauliQVec, pauliAlphaVec, polarizabilityVec;
    vector<double> localDipolesVec, localQuadrupolesVec;
    vector<mm_int4> multipoleParticlesVec;
    for (int i = 0; i < numParticles; i++) {
        double charge, coreCharge, alpha, epsilon, damping, c6, pauliK, pauliQ, pauliAlpha, polarizability;
        int axisType, atomX, atomY, atomZ;
        vector<double> dipole, quadrupole;
        force.getParticleParameters(i, charge, dipole, quadrupole, coreCharge, alpha, epsilon, damping, c6, pauliK, pauliQ, pauliAlpha,
                                    polarizability, axisType, atomZ, atomX, atomY);
        coreChargeVec.push_back(coreCharge);
        valenceChargeVec.push_back(charge-coreCharge);
        alphaVec.push_back(alpha);
        epsilonVec.push_back(epsilon);
        dampingVec.push_back(damping);
        c6Vec.push_back(c6);
        pauliKVec.push_back(pauliK);
        pauliQVec.push_back(pauliQ);
        pauliAlphaVec.push_back(pauliAlpha);
        polarizabilityVec.push_back(polarizability);
        multipoleParticlesVec.push_back(mm_int4(atomX, atomY, atomZ, axisType));
        for (int j = 0; j < 3; j++)
            localDipolesVec.push_back(dipole[j]);
        localQuadrupolesVec.push_back(quadrupole[0]);
        localQuadrupolesVec.push_back(quadrupole[1]);
        localQuadrupolesVec.push_back(quadrupole[2]);
        localQuadrupolesVec.push_back(quadrupole[4]);
        localQuadrupolesVec.push_back(quadrupole[5]);
    }
    int paddedNumAtoms = cc.getPaddedNumAtoms();
    for (int i = numParticles; i < paddedNumAtoms; i++) {
        coreChargeVec.push_back(0);
        valenceChargeVec.push_back(0);
        alphaVec.push_back(0);
        epsilonVec.push_back(0);
        dampingVec.push_back(0);
        c6Vec.push_back(0);
        pauliKVec.push_back(0);
        pauliQVec.push_back(0);
        pauliAlphaVec.push_back(0);
        polarizabilityVec.push_back(0);
        multipoleParticlesVec.push_back(mm_int4(0, 0, 0, 0));
        for (int j = 0; j < 3; j++)
            localDipolesVec.push_back(0);
        for (int j = 0; j < 5; j++)
            localQuadrupolesVec.push_back(0);
    }
    coreCharge.upload(coreChargeVec, true);
    valenceCharge.upload(valenceChargeVec, true);
    alpha.upload(alphaVec, true);
    epsilon.upload(epsilonVec, true);
    damping.upload(dampingVec, true);
    c6.upload(c6Vec, true);
    pauliK.upload(pauliKVec, true);
    pauliQ.upload(pauliQVec, true);
    pauliAlpha.upload(pauliAlphaVec, true);
    polarizability.upload(polarizabilityVec, true);
    multipoleParticles.upload(multipoleParticlesVec);
    localDipoles.upload(localDipolesVec, true);
    localQuadrupoles.upload(localQuadrupolesVec, true);
    
    // Record the per-exception parameters.

    vector<double> exceptionScaleVec[6];
    vector<mm_int2> exceptionAtomsVec;
    for (int i = 0; i < force.getNumExceptions(); i++) {
        int particle1, particle2;
        double multipoleMultipoleScale, dipoleMultipoleScale, dipoleDipoleScale, dispersionScale, repulsionScale, chargeTransferScale;
        force.getExceptionParameters(i, particle1, particle2, multipoleMultipoleScale, dipoleMultipoleScale, dipoleDipoleScale, dispersionScale, repulsionScale, chargeTransferScale);
        if (usePME || multipoleMultipoleScale != 0 || dipoleMultipoleScale != 0 || dipoleDipoleScale != 0 || dispersionScale != 0 || repulsionScale != 0 || chargeTransferScale != 0) {
            exceptionAtomsVec.push_back(mm_int2(particle1, particle2));
            exceptionScaleVec[0].push_back(multipoleMultipoleScale);
            exceptionScaleVec[1].push_back(dipoleMultipoleScale);
            exceptionScaleVec[2].push_back(dipoleDipoleScale);
            exceptionScaleVec[3].push_back(dispersionScale);
            exceptionScaleVec[4].push_back(repulsionScale);
            exceptionScaleVec[5].push_back(chargeTransferScale);
        }
    }
    if (exceptionAtomsVec.size() > 0) {
        if (!exceptionAtoms.isInitialized() || exceptionAtoms.getSize() != exceptionAtomsVec.size())
            throw OpenMMException("updateParametersInContext: The number of exceptions has changed");
        exceptionAtoms.upload(exceptionAtomsVec);
        for (int i = 0; i < 6; i++)
            exceptionScales[i].upload(exceptionScaleVec[i], true);
    }
    else if (exceptionAtoms.isInitialized())
        throw OpenMMException("updateParametersInContext: The number of exceptions has changed");
    cc.invalidateMolecules();
    multipolesAreValid = false;
}

void CommonCalcHippoNonbondedForceKernel::getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    alpha = pmeAlpha;
    nx = gridSizeX;
    ny = gridSizeY;
    nz = gridSizeZ;
}

void CommonCalcHippoNonbondedForceKernel::getDPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    alpha = dpmeAlpha;
    nx = dispersionGridSizeX;
    ny = dispersionGridSizeY;
    nz = dispersionGridSizeZ;
}
