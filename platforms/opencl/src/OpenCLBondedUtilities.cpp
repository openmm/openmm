/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2011-2022 Stanford University and the Authors.      *
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

#include "OpenCLBondedUtilities.h"
#include "OpenCLContext.h"
#include "OpenCLExpressionUtilities.h"
#include "openmm/OpenMMException.h"
#include "OpenCLNonbondedUtilities.h"
#include <iostream>

using namespace OpenMM;
using namespace std;

OpenCLBondedUtilities::OpenCLBondedUtilities(OpenCLContext& context) : context(context), maxBonds(0), allGroups(0), hasInitializedKernels(false) {
}

void OpenCLBondedUtilities::addInteraction(const vector<vector<int> >& atoms, const string& source, int group) {
    if (atoms.size() > 0) {
        forceAtoms.push_back(atoms);
        forceSource.push_back(source);
        forceGroup.push_back(group);
        allGroups |= 1<<group;
        int width = 1;
        while (width < (int) atoms[0].size())
            width *= 2;
        indexWidth.push_back(width);
    }
}

string OpenCLBondedUtilities::addArgument(cl::Memory& data, const string& type) {
    arguments.push_back(&data);
    argTypes.push_back(type);
    return "customArg"+context.intToString(arguments.size());
}

string OpenCLBondedUtilities::addArgument(ArrayInterface& data, const string& type) {
    return addArgument(context.unwrap(data).getDeviceBuffer(), type);
}

string OpenCLBondedUtilities::addEnergyParameterDerivative(const string& param) {
    // See if the parameter has already been added.
    
    int index;
    for (index = 0; index < energyParameterDerivatives.size(); index++)
        if (param == energyParameterDerivatives[index])
            break;
    if (index == energyParameterDerivatives.size())
        energyParameterDerivatives.push_back(param);
    context.addEnergyParameterDerivative(param);
    return string("energyParamDeriv")+context.intToString(index);
}

void OpenCLBondedUtilities::addPrefixCode(const string& source) {
    for (int i = 0; i < (int) prefixCode.size(); i++)
        if (prefixCode[i] == source)
            return;
    prefixCode.push_back(source);
}

void OpenCLBondedUtilities::initialize(const System& system) {
    int numForces = forceAtoms.size();
    if (numForces == 0)
        return;
    
    // Build the lists of atom indices.
    
    atomIndices.resize(numForces);
    for (int i = 0; i < numForces; i++) {
        int numBonds = forceAtoms[i].size();
        int numAtoms = forceAtoms[i][0].size();
        int width = indexWidth[i];
        vector<cl_uint> indexVec(width*numBonds);
        for (int bond = 0; bond < numBonds; bond++) {
            for (int atom = 0; atom < numAtoms; atom++)
                indexVec[bond*width+atom] = forceAtoms[i][bond][atom];
        }
        atomIndices[i].initialize<cl_uint>(context, indexVec.size(), "bondedIndices");
        atomIndices[i].upload(indexVec);
    }

    // Create the kernel.

    stringstream s;
    for (int i = 0; i < (int) prefixCode.size(); i++)
        s<<prefixCode[i];
    s<<"__kernel void computeBondedForces(__global unsigned long* restrict forceBuffers, __global mixed* restrict energyBuffer, __global const real4* restrict posq, int groups, real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ";
    for (int force = 0; force < numForces; force++) {
        string indexType = "uint"+(indexWidth[force] == 1 ? "" : context.intToString(indexWidth[force]));
        s<<", __global const "<<indexType<<"* restrict atomIndices"<<force;
    }
    for (int i = 0; i < (int) arguments.size(); i++)
        s<<", __global "<<argTypes[i]<<"* customArg"<<(i+1);
    if (energyParameterDerivatives.size() > 0)
        s<<", __global mixed* restrict energyParamDerivs";
    s<<") {\n";
    s<<"mixed energy = 0;\n";
    for (int i = 0; i < energyParameterDerivatives.size(); i++)
        s<<"mixed energyParamDeriv"<<i<<" = 0;\n";
    for (int force = 0; force < numForces; force++)
        s<<createForceSource(force, forceAtoms[force].size(), forceAtoms[force][0].size(), forceGroup[force], forceSource[force]);
    s<<"energyBuffer[get_global_id(0)] += energy;\n";
    const vector<string>& allParamDerivNames = context.getEnergyParamDerivNames();
    int numDerivs = allParamDerivNames.size();
    for (int i = 0; i < energyParameterDerivatives.size(); i++)
        for (int index = 0; index < numDerivs; index++)
            if (allParamDerivNames[index] == energyParameterDerivatives[i])
                s<<"energyParamDerivs[get_global_id(0)*"<<numDerivs<<"+"<<index<<"] += energyParamDeriv"<<i<<";\n";
    s<<"}\n";
    map<string, string> defines;
    defines["PADDED_NUM_ATOMS"] = context.intToString(context.getPaddedNumAtoms());
    cl::Program program = context.createProgram(s.str(), defines);
    kernel = cl::Kernel(program, "computeBondedForces");
    forceAtoms.clear();
    forceSource.clear();
}

string OpenCLBondedUtilities::createForceSource(int forceIndex, int numBonds, int numAtoms, int group, const string& computeForce) {
    maxBonds = max(maxBonds, numBonds);
    int width = 1;
    while (width < numAtoms)
        width *= 2;
    string suffix1[] = {""};
    string suffix4[] = {".x", ".y", ".z", ".w"};
    string suffix16[] = {".s0", ".s1", ".s2", ".s3", ".s4", ".s5", ".s6", ".s7",
        ".s8", ".s9", ".s10", ".s11", ".s12", ".s13", ".s14", ".s15"};
    string* suffix;
    if (width == 1)
        suffix = suffix1;
    else if (width <= 4)
        suffix = suffix4;
    else
        suffix = suffix16;
    string indexType = "uint"+(width == 1 ? "" : context.intToString(width));
    stringstream s;
    s<<"if ((groups&"<<(1<<group)<<") != 0)\n";
    s<<"for (unsigned int index = get_global_id(0); index < "<<numBonds<<"; index += get_global_size(0)) {\n";
    s<<"    "<<indexType<<" atoms = atomIndices"<<forceIndex<<"[index];\n";
    for (int i = 0; i < numAtoms; i++) {
        s<<"    unsigned int atom"<<(i+1)<<" = atoms"<<suffix[i]<<";\n";
        s<<"    real4 pos"<<(i+1)<<" = posq[atom"<<(i+1)<<"];\n";
    }
    s<<computeForce<<"\n";
    for (int i = 0; i < numAtoms; i++) {
        s<<"    {\n";
        s<<"    ATOMIC_ADD(&forceBuffers[atom"<<(i+1)<<"], (mm_ulong) realToFixedPoint(force"<<(i+1)<<".x));\n";
        s<<"    ATOMIC_ADD(&forceBuffers[atom"<<(i+1)<<"+PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(force"<<(i+1)<<".y));\n";
        s<<"    ATOMIC_ADD(&forceBuffers[atom"<<(i+1)<<"+2*PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(force"<<(i+1)<<".z));\n";
        s<<"    }\n";
    }
    s<<"}\n";
    return s.str();
}

void OpenCLBondedUtilities::computeInteractions(int groups) {
    if ((groups&allGroups) == 0)
        return;
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
        int index = 0;
        kernel.setArg<cl::Buffer>(index++, context.getLongForceBuffer().getDeviceBuffer());
        kernel.setArg<cl::Buffer>(index++, context.getEnergyBuffer().getDeviceBuffer());
        kernel.setArg<cl::Buffer>(index++, context.getPosq().getDeviceBuffer());
        index += 6;
        for (int j = 0; j < (int) atomIndices.size(); j++)
            kernel.setArg<cl::Buffer>(index++, atomIndices[j].getDeviceBuffer());
        for (int j = 0; j < (int) arguments.size(); j++)
            kernel.setArg<cl::Memory>(index++, *arguments[j]);
        if (energyParameterDerivatives.size() > 0)
            kernel.setArg<cl::Memory>(index++, context.getEnergyParamDerivBuffer().getDeviceBuffer());
    }
    kernel.setArg<cl_int>(3, groups);
    if (context.getUseDoublePrecision()) {
        kernel.setArg<mm_double4>(4, context.getPeriodicBoxSizeDouble());
        kernel.setArg<mm_double4>(5, context.getInvPeriodicBoxSizeDouble());
        kernel.setArg<mm_double4>(6, context.getPeriodicBoxVecXDouble());
        kernel.setArg<mm_double4>(7, context.getPeriodicBoxVecYDouble());
        kernel.setArg<mm_double4>(8, context.getPeriodicBoxVecZDouble());
    }
    else {
        kernel.setArg<mm_float4>(4, context.getPeriodicBoxSize());
        kernel.setArg<mm_float4>(5, context.getInvPeriodicBoxSize());
        kernel.setArg<mm_float4>(6, context.getPeriodicBoxVecX());
        kernel.setArg<mm_float4>(7, context.getPeriodicBoxVecY());
        kernel.setArg<mm_float4>(8, context.getPeriodicBoxVecZ());
    }
    context.executeKernel(kernel, maxBonds);
}
