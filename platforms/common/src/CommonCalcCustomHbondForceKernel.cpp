/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2025 Stanford University and the Authors.      *
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

#include "openmm/common/CommonCalcCustomHbondForceKernel.h"
#include "openmm/common/CommonKernelUtilities.h"
#include "openmm/common/ContextSelector.h"
#include "openmm/common/ExpressionUtilities.h"
#include "openmm/Context.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/CustomHbondForceImpl.h"
#include "CommonKernelSources.h"
#include "SimTKOpenMMRealType.h"
#include "lepton/CustomFunction.h"
#include "lepton/ExpressionTreeNode.h"
#include "lepton/Operation.h"
#include "lepton/Parser.h"
#include "lepton/ParsedExpression.h"

using namespace OpenMM;
using namespace std;
using namespace Lepton;

class CommonCalcCustomHbondForceKernel::ForceInfo : public ComputeForceInfo {
public:
    ForceInfo(const CustomHbondForce& force) : force(force) {
    }
    bool areParticlesIdentical(int particle1, int particle2) {
        return true;
    }
    int getNumParticleGroups() {
        return force.getNumDonors()+force.getNumAcceptors()+force.getNumExclusions();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        int p1, p2, p3;
        thread_local static vector<double> parameters;
        if (index < force.getNumDonors()) {
            force.getDonorParameters(index, p1, p2, p3, parameters);
            particles.clear();
            particles.push_back(p1);
            if (p2 > -1)
                particles.push_back(p2);
            if (p3 > -1)
                particles.push_back(p3);
            return;
        }
        index -= force.getNumDonors();
        if (index < force.getNumAcceptors()) {
            force.getAcceptorParameters(index, p1, p2, p3, parameters);
            particles.clear();
            particles.push_back(p1);
            if (p2 > -1)
                particles.push_back(p2);
            if (p3 > -1)
                particles.push_back(p3);
            return;
        }
        index -= force.getNumAcceptors();
        int donor, acceptor;
        force.getExclusionParticles(index, donor, acceptor);
        particles.clear();
        force.getDonorParameters(donor, p1, p2, p3, parameters);
        particles.push_back(p1);
        if (p2 > -1)
            particles.push_back(p2);
        if (p3 > -1)
            particles.push_back(p3);
        force.getAcceptorParameters(acceptor, p1, p2, p3, parameters);
        particles.push_back(p1);
        if (p2 > -1)
            particles.push_back(p2);
        if (p3 > -1)
            particles.push_back(p3);
    }
    bool areGroupsIdentical(int group1, int group2) {
        int p1, p2, p3;
        thread_local static vector<double> params1, params2;
        if (group1 < force.getNumDonors() && group2 < force.getNumDonors()) {
            force.getDonorParameters(group1, p1, p2, p3, params1);
            force.getDonorParameters(group2, p1, p2, p3, params2);
            return (params1 == params2);
        }
        if (group1 < force.getNumDonors() || group2 < force.getNumDonors())
            return false;
        group1 -= force.getNumDonors();
        group2 -= force.getNumDonors();
        if (group1 < force.getNumAcceptors() && group2 < force.getNumAcceptors()) {
            force.getAcceptorParameters(group1, p1, p2, p3, params1);
            force.getAcceptorParameters(group2, p1, p2, p3, params2);
            return (params1 == params2);
        }
        if (group1 < force.getNumAcceptors() || group2 < force.getNumAcceptors())
            return false;
        return true;
    }
private:
    const CustomHbondForce& force;
};

CommonCalcCustomHbondForceKernel::~CommonCalcCustomHbondForceKernel() {
    ContextSelector selector(cc);
    if (donorParams != NULL)
        delete donorParams;
    if (acceptorParams != NULL)
        delete acceptorParams;
}

static void applyDonorAndAcceptorForces(stringstream& apply, int atom, const string& value, bool trim=true) {
    string forceNames[] = {"f1", "f2", "f3"};
    string toAdd = (trim ? "trimTo3("+value+")" : value);
    if (atom < 3)
        apply << "localData[tbx+index]." << forceNames[atom]<<" += "<<toAdd<<";\n";
    else
        apply << forceNames[atom-3]<<" += "<<toAdd<<";\n";
}

void CommonCalcCustomHbondForceKernel::initialize(const System& system, const CustomHbondForce& force) {
    // Record the lists of donors and acceptors, and the parameters for each one.

    ContextSelector selector(cc);
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumDonors()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumDonors()/numContexts;
    numDonors = endIndex-startIndex;
    numAcceptors = force.getNumAcceptors();
    if (numDonors == 0 || numAcceptors == 0)
        return;
    int numParticles = system.getNumParticles();
    donors.initialize<mm_int4>(cc, numDonors, "customHbondDonors");
    acceptors.initialize<mm_int4>(cc, numAcceptors, "customHbondAcceptors");
    donorParams = new ComputeParameterSet(cc, force.getNumPerDonorParameters(), numDonors, "customHbondDonorParameters");
    acceptorParams = new ComputeParameterSet(cc, force.getNumPerAcceptorParameters(), numAcceptors, "customHbondAcceptorParameters");
    if (force.getNumGlobalParameters() > 0)
        globals.initialize<float>(cc, force.getNumGlobalParameters(), "customHbondGlobals");
    vector<vector<float> > donorParamVector(numDonors);
    vector<mm_int4> donorVector(numDonors);
    for (int i = 0; i < numDonors; i++) {
        vector<double> parameters;
        force.getDonorParameters(startIndex+i, donorVector[i].x, donorVector[i].y, donorVector[i].z, parameters);
        donorParamVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            donorParamVector[i][j] = (float) parameters[j];
    }
    donors.upload(donorVector);
    donorParams->setParameterValues(donorParamVector);
    vector<vector<float> > acceptorParamVector(numAcceptors);
    vector<mm_int4> acceptorVector(numAcceptors);
    for (int i = 0; i < numAcceptors; i++) {
        vector<double> parameters;
        force.getAcceptorParameters(i, acceptorVector[i].x, acceptorVector[i].y, acceptorVector[i].z, parameters);
        acceptorParamVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            acceptorParamVector[i][j] = (float) parameters[j];
    }
    acceptors.upload(acceptorVector);
    acceptorParams->setParameterValues(acceptorParamVector);
    info = new ForceInfo(force);
    cc.addForce(info);

    // Decide whether to use bounding boxes to accelerate the calculation.

    int numDonorBlocks = (numDonors+31)/32;
    int numAcceptorBlocks = (numAcceptors+31)/32;
    useBoundingBoxes = (force.getNonbondedMethod() != CustomHbondForce::NoCutoff && numDonorBlocks*numAcceptorBlocks > cc.getNumThreadBlocks());
    if (useBoundingBoxes) {
        int elementSize = (cc.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
        donorBlockCenter.initialize(cc, numDonorBlocks, 4*elementSize, "donorBlockCenter");
        donorBlockSize.initialize(cc, numDonorBlocks, 4*elementSize, "donorBlockSize");
        acceptorBlockCenter.initialize(cc, numAcceptorBlocks, 4*elementSize, "acceptorBlockCenter");
        acceptorBlockSize.initialize(cc, numAcceptorBlocks, 4*elementSize, "acceptorBlockSize");
    }

    // Record exclusions.

    vector<mm_int4> donorExclusionVector(numDonors, mm_int4(-1, -1, -1, -1));
    vector<mm_int4> acceptorExclusionVector(numAcceptors, mm_int4(-1, -1, -1, -1));
    for (int i = 0; i < force.getNumExclusions(); i++) {
        int donor, acceptor;
        force.getExclusionParticles(i, donor, acceptor);
        if (donor < startIndex || donor >= endIndex)
            continue;
        donor -= startIndex;
        if (donorExclusionVector[donor].x == -1)
            donorExclusionVector[donor].x = acceptor;
        else if (donorExclusionVector[donor].y == -1)
            donorExclusionVector[donor].y = acceptor;
        else if (donorExclusionVector[donor].z == -1)
            donorExclusionVector[donor].z = acceptor;
        else if (donorExclusionVector[donor].w == -1)
            donorExclusionVector[donor].w = acceptor;
        else
            throw OpenMMException("CustomHbondForce: this platform does not support more than four exclusions per donor");
        if (acceptorExclusionVector[acceptor].x == -1)
            acceptorExclusionVector[acceptor].x = donor;
        else if (acceptorExclusionVector[acceptor].y == -1)
            acceptorExclusionVector[acceptor].y = donor;
        else if (acceptorExclusionVector[acceptor].z == -1)
            acceptorExclusionVector[acceptor].z = donor;
        else if (acceptorExclusionVector[acceptor].w == -1)
            acceptorExclusionVector[acceptor].w = donor;
        else
            throw OpenMMException("CustomHbondForce: this platform does not support more than four exclusions per acceptor");
    }
    donorExclusions.initialize<mm_int4>(cc, numDonors, "customHbondDonorExclusions");
    acceptorExclusions.initialize<mm_int4>(cc, numAcceptors, "customHbondAcceptorExclusions");
    donorExclusions.upload(donorExclusionVector);
    acceptorExclusions.upload(acceptorExclusionVector);

    // Record the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    vector<pair<string, string> > functionDefinitions;
    vector<const TabulatedFunction*> functionList;
    stringstream tableArgs;
    tabulatedFunctionArrays.resize(force.getNumTabulatedFunctions());
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        functionList.push_back(&force.getTabulatedFunction(i));
        string name = force.getTabulatedFunctionName(i);
        tabulatedFunctionUpdateCount[name] = force.getTabulatedFunction(i).getUpdateCount();
        string arrayName = "table"+cc.intToString(i);
        functionDefinitions.push_back(make_pair(name, arrayName));
        functions[name] = cc.getExpressionUtilities().getFunctionPlaceholder(force.getTabulatedFunction(i));
        int width;
        vector<float> f = cc.getExpressionUtilities().computeFunctionCoefficients(force.getTabulatedFunction(i), width);
        tabulatedFunctionArrays[i].initialize<float>(cc, f.size(), "TabulatedFunction");
        tabulatedFunctionArrays[i].upload(f);
        tableArgs << ", GLOBAL const float";
        if (width > 1)
            tableArgs << width;
        tableArgs << "* RESTRICT " << arrayName;
    }

    // Record information about parameters.

    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (float) force.getGlobalParameterDefaultValue(i);
    }
    if (globals.isInitialized())
        globals.upload(globalParamValues);
    map<string, string> variables;
    for (int i = 0; i < force.getNumPerDonorParameters(); i++) {
        const string& name = force.getPerDonorParameterName(i);
        variables[name] = "donorParams"+donorParams->getParameterSuffix(i);
    }
    for (int i = 0; i < force.getNumPerAcceptorParameters(); i++) {
        const string& name = force.getPerAcceptorParameterName(i);
        variables[name] = "acceptorParams"+acceptorParams->getParameterSuffix(i);
    }
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        const string& name = force.getGlobalParameterName(i);
        variables[name] = "globals["+cc.intToString(i)+"]";
    }

    // Now to generate the kernel.  First, it needs to calculate all distances, angles,
    // and dihedrals the expression depends on.

    map<string, vector<int> > distances;
    map<string, vector<int> > angles;
    map<string, vector<int> > dihedrals;
    Lepton::ParsedExpression energyExpression = CustomHbondForceImpl::prepareExpression(force, functions, distances, angles, dihedrals);
    map<string, Lepton::ParsedExpression> forceExpressions;
    set<string> computedDeltas;
    computedDeltas.insert("D1A1");
    string atomNames[] = {"A1", "A2", "A3", "D1", "D2", "D3"};
    string atomNamesLower[] = {"a1", "a2", "a3", "d1", "d2", "d3"};
    stringstream compute, extraArgs;
    int index = 0;
    for (auto& distance : distances) {
        const vector<int>& atoms = distance.second;
        string deltaName = atomNames[atoms[0]]+atomNames[atoms[1]];
        if (computedDeltas.count(deltaName) == 0) {
            compute << "real4 delta"+deltaName+" = delta("+atomNamesLower[atoms[0]]+", "+atomNamesLower[atoms[1]]+", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);\n";
            computedDeltas.insert(deltaName);
        }
        compute << "real r_"+deltaName+" = SQRT(delta"+deltaName+".w);\n";
        variables[distance.first] = "r_"+deltaName;
        forceExpressions["real dEdDistance"+cc.intToString(index)+" = "] = energyExpression.differentiate(distance.first).optimize();
        index++;
    }
    index = 0;
    for (auto& angle : angles) {
        const vector<int>& atoms = angle.second;
        string deltaName1 = atomNames[atoms[1]]+atomNames[atoms[0]];
        string deltaName2 = atomNames[atoms[1]]+atomNames[atoms[2]];
        string angleName = "angle_"+atomNames[atoms[0]]+atomNames[atoms[1]]+atomNames[atoms[2]];
        if (computedDeltas.count(deltaName1) == 0) {
            compute << "real4 delta"+deltaName1+" = delta("+atomNamesLower[atoms[1]]+", "+atomNamesLower[atoms[0]]+", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);\n";
            computedDeltas.insert(deltaName1);
        }
        if (computedDeltas.count(deltaName2) == 0) {
            compute << "real4 delta"+deltaName2+" = delta("+atomNamesLower[atoms[1]]+", "+atomNamesLower[atoms[2]]+", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);\n";
            computedDeltas.insert(deltaName2);
        }
        compute << "real "+angleName+" = computeAngle(delta"+deltaName1+", delta"+deltaName2+");\n";
        variables[angle.first] = angleName;
        forceExpressions["real dEdAngle"+cc.intToString(index)+" = "] = energyExpression.differentiate(angle.first).optimize();
        index++;
    }
    index = 0;
    for (auto& dihedral : dihedrals) {
        const vector<int>& atoms = dihedral.second;
        string deltaName1 = atomNames[atoms[0]]+atomNames[atoms[1]];
        string deltaName2 = atomNames[atoms[2]]+atomNames[atoms[1]];
        string deltaName3 = atomNames[atoms[2]]+atomNames[atoms[3]];
        string crossName1 = "cross_"+deltaName1+"_"+deltaName2;
        string crossName2 = "cross_"+deltaName2+"_"+deltaName3;
        string dihedralName = "dihedral_"+atomNames[atoms[0]]+atomNames[atoms[1]]+atomNames[atoms[2]]+atomNames[atoms[3]];
        if (computedDeltas.count(deltaName1) == 0) {
            compute << "real4 delta"+deltaName1+" = delta("+atomNamesLower[atoms[0]]+", "+atomNamesLower[atoms[1]]+", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);\n";
            computedDeltas.insert(deltaName1);
        }
        if (computedDeltas.count(deltaName2) == 0) {
            compute << "real4 delta"+deltaName2+" = delta("+atomNamesLower[atoms[2]]+", "+atomNamesLower[atoms[1]]+", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);\n";
            computedDeltas.insert(deltaName2);
        }
        if (computedDeltas.count(deltaName3) == 0) {
            compute << "real4 delta"+deltaName3+" = delta("+atomNamesLower[atoms[2]]+", "+atomNamesLower[atoms[3]]+", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);\n";
            computedDeltas.insert(deltaName3);
        }
        compute << "real4 "+crossName1+" = computeCross(delta"+deltaName1+", delta"+deltaName2+");\n";
        compute << "real4 "+crossName2+" = computeCross(delta"+deltaName2+", delta"+deltaName3+");\n";
        compute << "real "+dihedralName+" = computeAngle("+crossName1+", "+crossName2+");\n";
        compute << dihedralName+" *= (delta"+deltaName1+".x*"+crossName2+".x + delta"+deltaName1+".y*"+crossName2+".y + delta"+deltaName1+".z*"+crossName2+".z < 0 ? -1 : 1);\n";
        variables[dihedral.first] = dihedralName;
        forceExpressions["real dEdDihedral"+cc.intToString(index)+" = "] = energyExpression.differentiate(dihedral.first).optimize();
        index++;
    }

    // Next it needs to load parameters from global memory.

    if (force.getNumGlobalParameters() > 0)
        extraArgs << ", GLOBAL const float* RESTRICT globals";
    for (int i = 0; i < (int) donorParams->getParameterInfos().size(); i++) {
        ComputeParameterInfo& parameter = donorParams->getParameterInfos()[i];
        extraArgs << ", GLOBAL const "+parameter.getType()+"* RESTRICT donor"+parameter.getName();
        compute << parameter.getType()+" donorParams"+cc.intToString(i+1)+" = donor"+parameter.getName()+"[donorIndex];\n";
    }
    for (int i = 0; i < (int) acceptorParams->getParameterInfos().size(); i++) {
        ComputeParameterInfo& parameter = acceptorParams->getParameterInfos()[i];
        extraArgs << ", GLOBAL const "+parameter.getType()+"* RESTRICT acceptor"+parameter.getName();
        compute << parameter.getType()+" acceptorParams"+cc.intToString(i+1)+" = acceptor"+parameter.getName()+"[acceptorIndex];\n";
    }

    // Now evaluate the expressions.

    forceExpressions["energy += "] = energyExpression;
    compute << cc.getExpressionUtilities().createExpressions(forceExpressions, variables, functionList, functionDefinitions, "temp");

    // Finally, apply forces to atoms.

    index = 0;
    for (auto& distance : distances) {
        const vector<int>& atoms = distance.second;
        string deltaName = atomNames[atoms[0]]+atomNames[atoms[1]];
        string value = "(dEdDistance"+cc.intToString(index)+"/r_"+deltaName+")*delta"+deltaName;
        applyDonorAndAcceptorForces(compute, atoms[0], "-"+value);
        applyDonorAndAcceptorForces(compute, atoms[1], value);
        index++;
    }
    index = 0;
    for (auto& angle : angles) {
        const vector<int>& atoms = angle.second;
        string deltaName1 = atomNames[atoms[1]]+atomNames[atoms[0]];
        string deltaName2 = atomNames[atoms[1]]+atomNames[atoms[2]];
        compute << "{\n";
        compute << "real3 crossProd = trimTo3(cross(delta"+deltaName2+", delta"+deltaName1+"));\n";
        compute << "real lengthCross = max(SQRT(dot(crossProd,crossProd)), (real) 1e-6f);\n";
        compute << "real3 deltaCross0 = -cross(trimTo3(delta"+deltaName1+"), crossProd)*dEdAngle"+cc.intToString(index)+"/(delta"+deltaName1+".w*lengthCross);\n";
        compute << "real3 deltaCross2 = cross(trimTo3(delta"+deltaName2+"), crossProd)*dEdAngle"+cc.intToString(index)+"/(delta"+deltaName2+".w*lengthCross);\n";
        compute << "real3 deltaCross1 = -(deltaCross0+deltaCross2);\n";
        applyDonorAndAcceptorForces(compute, atoms[0], "deltaCross0", false);
        applyDonorAndAcceptorForces(compute, atoms[1], "deltaCross1", false);
        applyDonorAndAcceptorForces(compute, atoms[2], "deltaCross2", false);
        compute << "}\n";
        index++;
    }
    index = 0;
    for (auto& dihedral : dihedrals) {
        const vector<int>& atoms = dihedral.second;
        string deltaName1 = atomNames[atoms[0]]+atomNames[atoms[1]];
        string deltaName2 = atomNames[atoms[2]]+atomNames[atoms[1]];
        string deltaName3 = atomNames[atoms[2]]+atomNames[atoms[3]];
        string crossName1 = "cross_"+deltaName1+"_"+deltaName2;
        string crossName2 = "cross_"+deltaName2+"_"+deltaName3;
        compute << "{\n";
        compute << "real r = SQRT(delta"+deltaName2+".w);\n";
        compute << "real4 ff;\n";
        compute << "ff.x = (-dEdDihedral"+cc.intToString(index)+"*r)/"+crossName1+".w;\n";
        compute << "ff.y = (delta"+deltaName1+".x*delta"+deltaName2+".x + delta"+deltaName1+".y*delta"+deltaName2+".y + delta"+deltaName1+".z*delta"+deltaName2+".z)/delta"+deltaName2+".w;\n";
        compute << "ff.z = (delta"+deltaName3+".x*delta"+deltaName2+".x + delta"+deltaName3+".y*delta"+deltaName2+".y + delta"+deltaName3+".z*delta"+deltaName2+".z)/delta"+deltaName2+".w;\n";
        compute << "ff.w = (dEdDihedral"+cc.intToString(index)+"*r)/"+crossName2+".w;\n";
        compute << "real4 internalF0 = ff.x*"+crossName1+";\n";
        compute << "real4 internalF3 = ff.w*"+crossName2+";\n";
        compute << "real4 s = ff.y*internalF0 - ff.z*internalF3;\n";
        applyDonorAndAcceptorForces(compute, atoms[0], "internalF0");
        applyDonorAndAcceptorForces(compute, atoms[1], "s-internalF0");
        applyDonorAndAcceptorForces(compute, atoms[2], "-s-internalF3");
        applyDonorAndAcceptorForces(compute, atoms[3], "internalF3");
        compute << "}\n";
        index++;
    }

    // Generate the kernels.

    map<string, string> replacements;
    replacements["COMPUTE_FORCE"] = compute.str();
    replacements["PARAMETER_ARGUMENTS"] = extraArgs.str()+tableArgs.str();
    map<string, string> defines;
    defines["PADDED_NUM_ATOMS"] = cc.intToString(cc.getPaddedNumAtoms());
    defines["NUM_DONORS"] = cc.intToString(numDonors);
    defines["NUM_ACCEPTORS"] = cc.intToString(numAcceptors);
    defines["NUM_DONOR_BLOCKS"] = cc.intToString(numDonorBlocks);
    defines["NUM_ACCEPTOR_BLOCKS"] = cc.intToString(numAcceptorBlocks);
    defines["M_PI"] = cc.doubleToString(M_PI);
    defines["THREAD_BLOCK_SIZE"] = "128";
    if (force.getNonbondedMethod() != CustomHbondForce::NoCutoff) {
        defines["USE_CUTOFF"] = "1";
        defines["CUTOFF_SQUARED"] = cc.doubleToString(force.getCutoffDistance()*force.getCutoffDistance());
    }
    if (force.getNonbondedMethod() != CustomHbondForce::NoCutoff && force.getNonbondedMethod() != CustomHbondForce::CutoffNonPeriodic)
        defines["USE_PERIODIC"] = "1";
    if (force.getNumExclusions() > 0)
        defines["USE_EXCLUSIONS"] = "1";
    if (useBoundingBoxes)
        defines["USE_BOUNDING_BOXES"] = "1";
    ComputeProgram program = cc.compileProgram(cc.replaceStrings(CommonKernelSources::customHbondForce, replacements), defines);
    blockBoundsKernel = program->createKernel("findBlockBounds");
    forceKernel = program->createKernel("computeHbondForces");
}

double CommonCalcCustomHbondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if (numDonors == 0 || numAcceptors == 0)
        return 0.0;
    ContextSelector selector(cc);
    if (globals.isInitialized()) {
        bool changed = false;
        for (int i = 0; i < (int) globalParamNames.size(); i++) {
            float value = (float) context.getParameter(globalParamNames[i]);
            if (value != globalParamValues[i])
                changed = true;
            globalParamValues[i] = value;
        }
        if (changed)
            globals.upload(globalParamValues);
    }
    if (!hasInitializedKernel) {
        hasInitializedKernel = true;
        if (useBoundingBoxes) {
            blockBoundsKernel->addArg(donors);
            blockBoundsKernel->addArg(acceptors);
            for (int i = 0; i < 5; i++)
                blockBoundsKernel->addArg(); // Periodic box size arguments are set when the kernel is executed.
            blockBoundsKernel->addArg(cc.getPosq());
            blockBoundsKernel->addArg(donorBlockCenter);
            blockBoundsKernel->addArg(donorBlockSize);
            blockBoundsKernel->addArg(acceptorBlockCenter);
            blockBoundsKernel->addArg(acceptorBlockSize);
        }
        forceKernel->addArg(cc.getLongForceBuffer());
        forceKernel->addArg(cc.getEnergyBuffer());
        forceKernel->addArg(cc.getPosq());
        forceKernel->addArg(donorExclusions);
        forceKernel->addArg(donors);
        forceKernel->addArg(acceptors);
        for (int i = 0; i < 5; i++)
            forceKernel->addArg(); // Periodic box size arguments are set when the kernel is executed.
        if (useBoundingBoxes) {
            forceKernel->addArg(donorBlockCenter);
            forceKernel->addArg(donorBlockSize);
            forceKernel->addArg(acceptorBlockCenter);
            forceKernel->addArg(acceptorBlockSize);
        }
        if (globals.isInitialized())
            forceKernel->addArg(globals);
        for (auto& parameter : donorParams->getParameterInfos())
            forceKernel->addArg(parameter.getArray());
        for (auto& parameter : acceptorParams->getParameterInfos())
            forceKernel->addArg(parameter.getArray());
        for (auto& function : tabulatedFunctionArrays)
            forceKernel->addArg(function);
    }
    if (useBoundingBoxes) {
        setPeriodicBoxArgs(cc, blockBoundsKernel, 2);
        blockBoundsKernel->execute(max(numDonors, numAcceptors));
    }
    setPeriodicBoxArgs(cc, forceKernel, 6);
    int numDonorBlocks = (numDonors+31)/32;
    int numAcceptorBlocks = (numAcceptors+31)/32;
    forceKernel->execute(numDonorBlocks*numAcceptorBlocks*32, cc.getSIMDWidth() < 32 ? 32 : 128);
    return 0.0;
}

void CommonCalcCustomHbondForceKernel::copyParametersToContext(ContextImpl& context, const CustomHbondForce& force) {
    ContextSelector selector(cc);
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumDonors()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumDonors()/numContexts;
    if (numDonors != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of donors has changed");
    if (numAcceptors != force.getNumAcceptors())
        throw OpenMMException("updateParametersInContext: The number of acceptors has changed");

    // Record the per-donor parameters.

    if (numDonors > 0) {
        vector<vector<float> > donorParamVector(numDonors);
        vector<double> parameters;
        for (int i = 0; i < numDonors; i++) {
            int d1, d2, d3;
            force.getDonorParameters(startIndex+i, d1, d2, d3, parameters);
            donorParamVector[i].resize(parameters.size());
            for (int j = 0; j < (int) parameters.size(); j++)
                donorParamVector[i][j] = (float) parameters[j];
        }
        donorParams->setParameterValues(donorParamVector);
    }

    // Record the per-acceptor parameters.

    if (numAcceptors > 0) {
        vector<vector<float> > acceptorParamVector(numAcceptors);
        vector<double> parameters;
        for (int i = 0; i < numAcceptors; i++) {
            int a1, a2, a3;
            force.getAcceptorParameters(i, a1, a2, a3, parameters);
            acceptorParamVector[i].resize(parameters.size());
            for (int j = 0; j < (int) parameters.size(); j++)
                acceptorParamVector[i][j] = (float) parameters[j];
        }
        acceptorParams->setParameterValues(acceptorParamVector);
    }

    // See if any tabulated functions have changed.

    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        string name = force.getTabulatedFunctionName(i);
        if (force.getTabulatedFunction(i).getUpdateCount() != tabulatedFunctionUpdateCount[name]) {
            tabulatedFunctionUpdateCount[name] = force.getTabulatedFunction(i).getUpdateCount();
            int width;
            vector<float> f = cc.getExpressionUtilities().computeFunctionCoefficients(force.getTabulatedFunction(i), width);
            tabulatedFunctionArrays[i].upload(f);
        }
    }

    // Mark that the current reordering may be invalid.

    cc.invalidateMolecules(info, false, true);
}
