/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2019 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
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

#include "openmm/common/CommonKernels.h"
#include "openmm/Context.h"
#include "openmm/internal/AndersenThermostatImpl.h"
#include "openmm/internal/CMAPTorsionForceImpl.h"
#include "CommonKernelSources.h"
#include "SimTKOpenMMRealType.h"
#include <algorithm>
#include <cmath>
#include <iterator>
#include <set>

using namespace OpenMM;
using namespace std;

class CommonCalcHarmonicBondForceKernel::ForceInfo : public ComputeForceInfo {
public:
    ForceInfo(const HarmonicBondForce& force) : force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumBonds();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        int particle1, particle2;
        double length, k;
        force.getBondParameters(index, particle1, particle2, length, k);
        particles.resize(2);
        particles[0] = particle1;
        particles[1] = particle2;
    }
    bool areGroupsIdentical(int group1, int group2) {
        int particle1, particle2;
        double length1, length2, k1, k2;
        force.getBondParameters(group1, particle1, particle2, length1, k1);
        force.getBondParameters(group2, particle1, particle2, length2, k2);
        return (length1 == length2 && k1 == k2);
    }
private:
    const HarmonicBondForce& force;
};

void CommonCalcHarmonicBondForceKernel::initialize(const System& system, const HarmonicBondForce& force) {
    cc.setAsCurrent();
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumBonds()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumBonds()/numContexts;
    numBonds = endIndex-startIndex;
    if (numBonds == 0)
        return;
    vector<vector<int> > atoms(numBonds, vector<int>(2));
    params.initialize<mm_float2>(cc, numBonds, "bondParams");
    vector<mm_float2> paramVector(numBonds);
    for (int i = 0; i < numBonds; i++) {
        double length, k;
        force.getBondParameters(startIndex+i, atoms[i][0], atoms[i][1], length, k);
        paramVector[i] = mm_float2((float) length, (float) k);
    }
    params.upload(paramVector);
    map<string, string> replacements;
    replacements["APPLY_PERIODIC"] = (force.usesPeriodicBoundaryConditions() ? "1" : "0");
    replacements["COMPUTE_FORCE"] = CommonKernelSources::harmonicBondForce;
    replacements["PARAMS"] = cc.getBondedUtilities().addArgument(params, "float2");
    cc.getBondedUtilities().addInteraction(atoms, cc.replaceStrings(CommonKernelSources::bondForce, replacements), force.getForceGroup());
    info = new ForceInfo(force);
    cc.addForce(info);
}

double CommonCalcHarmonicBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

void CommonCalcHarmonicBondForceKernel::copyParametersToContext(ContextImpl& context, const HarmonicBondForce& force) {
    cc.setAsCurrent();
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumBonds()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumBonds()/numContexts;
    if (numBonds != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of bonds has changed");
    if (numBonds == 0)
        return;
    
    // Record the per-bond parameters.
    
    vector<mm_float2> paramVector(numBonds);
    for (int i = 0; i < numBonds; i++) {
        int atom1, atom2;
        double length, k;
        force.getBondParameters(startIndex+i, atom1, atom2, length, k);
        paramVector[i] = mm_float2((float) length, (float) k);
    }
    params.upload(paramVector);
    
    // Mark that the current reordering may be invalid.
    
    cc.invalidateMolecules(info);
}

class CommonCalcHarmonicAngleForceKernel::ForceInfo : public ComputeForceInfo {
public:
    ForceInfo(const HarmonicAngleForce& force) : force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumAngles();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        int particle1, particle2, particle3;
        double angle, k;
        force.getAngleParameters(index, particle1, particle2, particle3, angle, k);
        particles.resize(3);
        particles[0] = particle1;
        particles[1] = particle2;
        particles[2] = particle3;
    }
    bool areGroupsIdentical(int group1, int group2) {
        int particle1, particle2, particle3;
        double angle1, angle2, k1, k2;
        force.getAngleParameters(group1, particle1, particle2, particle3, angle1, k1);
        force.getAngleParameters(group2, particle1, particle2, particle3, angle2, k2);
        return (angle1 == angle2 && k1 == k2);
    }
private:
    const HarmonicAngleForce& force;
};

void CommonCalcHarmonicAngleForceKernel::initialize(const System& system, const HarmonicAngleForce& force) {
    cc.setAsCurrent();
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumAngles()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumAngles()/numContexts;
    numAngles = endIndex-startIndex;
    if (numAngles == 0)
        return;
    vector<vector<int> > atoms(numAngles, vector<int>(3));
    params.initialize<mm_float2>(cc, numAngles, "angleParams");
    vector<mm_float2> paramVector(numAngles);
    for (int i = 0; i < numAngles; i++) {
        double angle, k;
        force.getAngleParameters(startIndex+i, atoms[i][0], atoms[i][1], atoms[i][2], angle, k);
        paramVector[i] = mm_float2((float) angle, (float) k);

    }
    params.upload(paramVector);
    map<string, string> replacements;
    replacements["APPLY_PERIODIC"] = (force.usesPeriodicBoundaryConditions() ? "1" : "0");
    replacements["COMPUTE_FORCE"] = CommonKernelSources::harmonicAngleForce;
    replacements["PARAMS"] = cc.getBondedUtilities().addArgument(params, "float2");
    cc.getBondedUtilities().addInteraction(atoms, cc.replaceStrings(CommonKernelSources::angleForce, replacements), force.getForceGroup());
    info = new ForceInfo(force);
    cc.addForce(info);
}

double CommonCalcHarmonicAngleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

void CommonCalcHarmonicAngleForceKernel::copyParametersToContext(ContextImpl& context, const HarmonicAngleForce& force) {
    cc.setAsCurrent();
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumAngles()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumAngles()/numContexts;
    if (numAngles != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of angles has changed");
    if (numAngles == 0)
        return;
    
    // Record the per-angle parameters.
    
    vector<mm_float2> paramVector(numAngles);
    for (int i = 0; i < numAngles; i++) {
        int atom1, atom2, atom3;
        double angle, k;
        force.getAngleParameters(startIndex+i, atom1, atom2, atom3, angle, k);
        paramVector[i] = mm_float2((float) angle, (float) k);
    }
    params.upload(paramVector);
    
    // Mark that the current reordering may be invalid.
    
    cc.invalidateMolecules();
}

class CommonCalcPeriodicTorsionForceKernel::ForceInfo : public ComputeForceInfo {
public:
    ForceInfo(const PeriodicTorsionForce& force) : force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumTorsions();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        int particle1, particle2, particle3, particle4, periodicity;
        double phase, k;
        force.getTorsionParameters(index, particle1, particle2, particle3, particle4, periodicity, phase, k);
        particles.resize(4);
        particles[0] = particle1;
        particles[1] = particle2;
        particles[2] = particle3;
        particles[3] = particle4;
    }
    bool areGroupsIdentical(int group1, int group2) {
        int particle1, particle2, particle3, particle4, periodicity1, periodicity2;
        double phase1, phase2, k1, k2;
        force.getTorsionParameters(group1, particle1, particle2, particle3, particle4, periodicity1, phase1, k1);
        force.getTorsionParameters(group2, particle1, particle2, particle3, particle4, periodicity2, phase2, k2);
        return (periodicity1 == periodicity2 && phase1 == phase2 && k1 == k2);
    }
private:
    const PeriodicTorsionForce& force;
};

void CommonCalcPeriodicTorsionForceKernel::initialize(const System& system, const PeriodicTorsionForce& force) {
    cc.setAsCurrent();
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumTorsions()/numContexts;
    numTorsions = endIndex-startIndex;
    if (numTorsions == 0)
        return;
    vector<vector<int> > atoms(numTorsions, vector<int>(4));
    params.initialize<mm_float4>(cc, numTorsions, "periodicTorsionParams");
    vector<mm_float4> paramVector(numTorsions);
    for (int i = 0; i < numTorsions; i++) {
        int periodicity;
        double phase, k;
        force.getTorsionParameters(startIndex+i, atoms[i][0], atoms[i][1], atoms[i][2], atoms[i][3], periodicity, phase, k);
        paramVector[i] = mm_float4((float) k, (float) phase, (float) periodicity, 0.0f);
    }
    params.upload(paramVector);
    map<string, string> replacements;
    replacements["APPLY_PERIODIC"] = (force.usesPeriodicBoundaryConditions() ? "1" : "0");
    replacements["COMPUTE_FORCE"] = CommonKernelSources::periodicTorsionForce;
    replacements["PARAMS"] = cc.getBondedUtilities().addArgument(params, "float4");
    cc.getBondedUtilities().addInteraction(atoms, cc.replaceStrings(CommonKernelSources::torsionForce, replacements), force.getForceGroup());
    info = new ForceInfo(force);
    cc.addForce(info);
}

double CommonCalcPeriodicTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

void CommonCalcPeriodicTorsionForceKernel::copyParametersToContext(ContextImpl& context, const PeriodicTorsionForce& force) {
    cc.setAsCurrent();
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumTorsions()/numContexts;
    if (numTorsions != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of torsions has changed");
    if (numTorsions == 0)
        return;
    
    // Record the per-torsion parameters.
    
    vector<mm_float4> paramVector(numTorsions);
    for (int i = 0; i < numTorsions; i++) {
        int atom1, atom2, atom3, atom4, periodicity;
        double phase, k;
        force.getTorsionParameters(startIndex+i, atom1, atom2, atom3, atom4, periodicity, phase, k);
        paramVector[i] = mm_float4((float) k, (float) phase, (float) periodicity, 0.0f);
    }
    params.upload(paramVector);
    
    // Mark that the current reordering may be invalid.
    
    cc.invalidateMolecules();
}

class CommonCalcRBTorsionForceKernel::ForceInfo : public ComputeForceInfo {
public:
    ForceInfo(const RBTorsionForce& force) : force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumTorsions();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        int particle1, particle2, particle3, particle4;
        double c0, c1, c2, c3, c4, c5;
        force.getTorsionParameters(index, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5);
        particles.resize(4);
        particles[0] = particle1;
        particles[1] = particle2;
        particles[2] = particle3;
        particles[3] = particle4;
    }
    bool areGroupsIdentical(int group1, int group2) {
        int particle1, particle2, particle3, particle4;
        double c0a, c0b, c1a, c1b, c2a, c2b, c3a, c3b, c4a, c4b, c5a, c5b;
        force.getTorsionParameters(group1, particle1, particle2, particle3, particle4, c0a, c1a, c2a, c3a, c4a, c5a);
        force.getTorsionParameters(group2, particle1, particle2, particle3, particle4, c0b, c1b, c2b, c3b, c4b, c5b);
        return (c0a == c0b && c1a == c1b && c2a == c2b && c3a == c3b && c4a == c4b && c5a == c5b);
    }
private:
    const RBTorsionForce& force;
};

void CommonCalcRBTorsionForceKernel::initialize(const System& system, const RBTorsionForce& force) {
    cc.setAsCurrent();
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumTorsions()/numContexts;
    numTorsions = endIndex-startIndex;
    if (numTorsions == 0)
        return;
    vector<vector<int> > atoms(numTorsions, vector<int>(4));
    params1.initialize<mm_float4>(cc, numTorsions, "rbTorsionParams1");
    params2.initialize<mm_float2>(cc, numTorsions, "rbTorsionParams2");
    vector<mm_float4> paramVector1(numTorsions);
    vector<mm_float2> paramVector2(numTorsions);
    for (int i = 0; i < numTorsions; i++) {
        double c0, c1, c2, c3, c4, c5;
        force.getTorsionParameters(startIndex+i, atoms[i][0], atoms[i][1], atoms[i][2], atoms[i][3], c0, c1, c2, c3, c4, c5);
        paramVector1[i] = mm_float4((float) c0, (float) c1, (float) c2, (float) c3);
        paramVector2[i] = mm_float2((float) c4, (float) c5);

    }
    params1.upload(paramVector1);
    params2.upload(paramVector2);
    map<string, string> replacements;
    replacements["APPLY_PERIODIC"] = (force.usesPeriodicBoundaryConditions() ? "1" : "0");
    replacements["COMPUTE_FORCE"] = CommonKernelSources::rbTorsionForce;
    replacements["PARAMS1"] = cc.getBondedUtilities().addArgument(params1, "float4");
    replacements["PARAMS2"] = cc.getBondedUtilities().addArgument(params2, "float2");
    cc.getBondedUtilities().addInteraction(atoms, cc.replaceStrings(CommonKernelSources::torsionForce, replacements), force.getForceGroup());
    info = new ForceInfo(force);
    cc.addForce(info);
}

double CommonCalcRBTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

void CommonCalcRBTorsionForceKernel::copyParametersToContext(ContextImpl& context, const RBTorsionForce& force) {
    cc.setAsCurrent();
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumTorsions()/numContexts;
    if (numTorsions != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of torsions has changed");
    if (numTorsions == 0)
        return;
    
    // Record the per-torsion parameters.
    
    vector<mm_float4> paramVector1(numTorsions);
    vector<mm_float2> paramVector2(numTorsions);
    for (int i = 0; i < numTorsions; i++) {
        int atom1, atom2, atom3, atom4;
        double c0, c1, c2, c3, c4, c5;
        force.getTorsionParameters(startIndex+i, atom1, atom2, atom3, atom4, c0, c1, c2, c3, c4, c5);
        paramVector1[i] = mm_float4((float) c0, (float) c1, (float) c2, (float) c3);
        paramVector2[i] = mm_float2((float) c4, (float) c5);
    }
    params1.upload(paramVector1);
    params2.upload(paramVector2);
    
    // Mark that the current reordering may be invalid.
    
    cc.invalidateMolecules();
}

class CommonCalcCMAPTorsionForceKernel::ForceInfo : public ComputeForceInfo {
public:
    ForceInfo(const CMAPTorsionForce& force) : force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumTorsions();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        int map, a1, a2, a3, a4, b1, b2, b3, b4;
        force.getTorsionParameters(index, map, a1, a2, a3, a4, b1, b2, b3, b4);
        particles.resize(8);
        particles[0] = a1;
        particles[1] = a2;
        particles[2] = a3;
        particles[3] = a4;
        particles[4] = b1;
        particles[5] = b2;
        particles[6] = b3;
        particles[7] = b4;
    }
    bool areGroupsIdentical(int group1, int group2) {
        int map1, map2, a1, a2, a3, a4, b1, b2, b3, b4;
        force.getTorsionParameters(group1, map1, a1, a2, a3, a4, b1, b2, b3, b4);
        force.getTorsionParameters(group2, map2, a1, a2, a3, a4, b1, b2, b3, b4);
        return (map1 == map2);
    }
private:
    const CMAPTorsionForce& force;
};

void CommonCalcCMAPTorsionForceKernel::initialize(const System& system, const CMAPTorsionForce& force) {
    cc.setAsCurrent();
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumTorsions()/numContexts;
    numTorsions = endIndex-startIndex;
    if (numTorsions == 0)
        return;
    int numMaps = force.getNumMaps();
    vector<mm_float4> coeffVec;
    mapPositionsVec.resize(numMaps);
    vector<double> energy;
    vector<vector<double> > c;
    int currentPosition = 0;
    for (int i = 0; i < numMaps; i++) {
        int size;
        force.getMapParameters(i, size, energy);
        CMAPTorsionForceImpl::calcMapDerivatives(size, energy, c);
        mapPositionsVec[i] = mm_int2(currentPosition, size);
        currentPosition += 4*size*size;
        for (int j = 0; j < size*size; j++) {
            coeffVec.push_back(mm_float4((float) c[j][0], (float) c[j][1], (float) c[j][2], (float) c[j][3]));
            coeffVec.push_back(mm_float4((float) c[j][4], (float) c[j][5], (float) c[j][6], (float) c[j][7]));
            coeffVec.push_back(mm_float4((float) c[j][8], (float) c[j][9], (float) c[j][10], (float) c[j][11]));
            coeffVec.push_back(mm_float4((float) c[j][12], (float) c[j][13], (float) c[j][14], (float) c[j][15]));
        }
    }
    vector<vector<int> > atoms(numTorsions, vector<int>(8));
    vector<int> torsionMapsVec(numTorsions);
    for (int i = 0; i < numTorsions; i++)
        force.getTorsionParameters(startIndex+i, torsionMapsVec[i], atoms[i][0], atoms[i][1], atoms[i][2], atoms[i][3], atoms[i][4], atoms[i][5], atoms[i][6], atoms[i][7]);
    coefficients.initialize<mm_float4>(cc, coeffVec.size(), "cmapTorsionCoefficients");
    mapPositions.initialize<mm_int2>(cc, numMaps, "cmapTorsionMapPositions");
    torsionMaps.initialize<int>(cc, numTorsions, "cmapTorsionMaps");
    coefficients.upload(coeffVec);
    mapPositions.upload(mapPositionsVec);
    torsionMaps.upload(torsionMapsVec);
    map<string, string> replacements;
    replacements["APPLY_PERIODIC"] = (force.usesPeriodicBoundaryConditions() ? "1" : "0");
    replacements["COEFF"] = cc.getBondedUtilities().addArgument(coefficients, "float4");
    replacements["MAP_POS"] = cc.getBondedUtilities().addArgument(mapPositions, "int2");
    replacements["MAPS"] = cc.getBondedUtilities().addArgument(torsionMaps, "int");
    cc.getBondedUtilities().addInteraction(atoms, cc.replaceStrings(CommonKernelSources::cmapTorsionForce, replacements), force.getForceGroup());
    info = new ForceInfo(force);
    cc.addForce(info);
}

double CommonCalcCMAPTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

void CommonCalcCMAPTorsionForceKernel::copyParametersToContext(ContextImpl& context, const CMAPTorsionForce& force) {
    int numMaps = force.getNumMaps();
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumTorsions()/numContexts;
    numTorsions = endIndex-startIndex;
    if (mapPositions.getSize() != numMaps)
        throw OpenMMException("updateParametersInContext: The number of maps has changed");
    if (torsionMaps.getSize() != numTorsions)
        throw OpenMMException("updateParametersInContext: The number of CMAP torsions has changed");

    // Update the maps.

    vector<mm_float4> coeffVec;
    vector<double> energy;
    vector<vector<double> > c;
    int currentPosition = 0;
    for (int i = 0; i < numMaps; i++) {
        int size;
        force.getMapParameters(i, size, energy);
        if (size != mapPositionsVec[i].y)
            throw OpenMMException("updateParametersInContext: The size of a map has changed");
        CMAPTorsionForceImpl::calcMapDerivatives(size, energy, c);
        currentPosition += 4*size*size;
        for (int j = 0; j < size*size; j++) {
            coeffVec.push_back(mm_float4((float) c[j][0], (float) c[j][1], (float) c[j][2], (float) c[j][3]));
            coeffVec.push_back(mm_float4((float) c[j][4], (float) c[j][5], (float) c[j][6], (float) c[j][7]));
            coeffVec.push_back(mm_float4((float) c[j][8], (float) c[j][9], (float) c[j][10], (float) c[j][11]));
            coeffVec.push_back(mm_float4((float) c[j][12], (float) c[j][13], (float) c[j][14], (float) c[j][15]));
        }
    }
    coefficients.upload(coeffVec);

    // Update the indices.

    vector<int> torsionMapsVec(numTorsions);
    for (int i = 0; i < numTorsions; i++) {
        int index[8];
        force.getTorsionParameters(i, torsionMapsVec[i], index[0], index[1], index[2], index[3], index[4], index[5], index[6], index[7]);
    }
    torsionMaps.upload(torsionMapsVec);
}

void CommonRemoveCMMotionKernel::initialize(const System& system, const CMMotionRemover& force) {
    cc.setAsCurrent();
    frequency = force.getFrequency();
    int numAtoms = cc.getNumAtoms();
    cmMomentum.initialize<mm_float3>(cc, cc.getPaddedNumAtoms(), "cmMomentum");
    double totalMass = 0.0;
    for (int i = 0; i < numAtoms; i++)
        totalMass += system.getParticleMass(i);
    map<string, string> defines;
    defines["INVERSE_TOTAL_MASS"] = cc.doubleToString(totalMass == 0 ? 0.0 : 1.0/totalMass);
    ComputeProgram program = cc.compileProgram(CommonKernelSources::removeCM, defines);
    kernel1 = program->createKernel("calcCenterOfMassMomentum");
    kernel1->addArg(numAtoms);
    kernel1->addArg(cc.getVelm());
    kernel1->addArg(cmMomentum);
    kernel2 = program->createKernel("removeCenterOfMassMomentum");
    kernel2->addArg(numAtoms);
    kernel2->addArg(cc.getVelm());
    kernel2->addArg(cmMomentum);
}

void CommonRemoveCMMotionKernel::execute(ContextImpl& context) {
    cc.setAsCurrent();
    kernel1->execute(cc.getNumAtoms(), 64);
    kernel2->execute(cc.getNumAtoms(), 64);
}
