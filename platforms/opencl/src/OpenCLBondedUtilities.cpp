/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2011-2012 Stanford University and the Authors.      *
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
#include "OpenCLExpressionUtilities.h"
#include "openmm/OpenMMException.h"
#include "OpenCLNonbondedUtilities.h"
#include <iostream>

using namespace OpenMM;
using namespace std;

OpenCLBondedUtilities::OpenCLBondedUtilities(OpenCLContext& context) : context(context), numForceBuffers(0), maxBonds(0), hasInitializedKernels(false) {
}

OpenCLBondedUtilities::~OpenCLBondedUtilities() {
    for (int i = 0; i < (int) atomIndices.size(); i++)
        delete atomIndices[i];
    for (int i = 0; i < (int) bufferIndices.size(); i++)
        delete bufferIndices[i];
}

void OpenCLBondedUtilities::addInteraction(const vector<vector<int> >& atoms, const string& source, int group) {
    if (atoms.size() > 0) {
        forceAtoms.push_back(atoms);
        forceSource.push_back(source);
        forceGroup.push_back(group);
        int width = 1;
        while (width < (int) atoms[0].size())
            width *= 2;
        indexWidth.push_back(width);
    }
}

std::string OpenCLBondedUtilities::addArgument(cl::Memory& data, const string& type) {
    arguments.push_back(&data);
    argTypes.push_back(type);
    return "customArg"+context.intToString(arguments.size());
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
    
    // Build the lists of atom indicse and buffer indices.
    
    vector<vector<cl_uint> > bufferVec(numForces);
    vector<vector<int> > bufferCounter(numForces, vector<int>(system.getNumParticles(), 0));
    vector<int> numBuffers(numForces, 0);
    for (int i = 0; i < numForces; i++) {
        int numBonds = forceAtoms[i].size();
        int numAtoms = forceAtoms[i][0].size();
        int width = indexWidth[i];
        vector<cl_uint> indexVec(width*numBonds);
        for (int bond = 0; bond < numBonds; bond++) {
            for (int atom = 0; atom < numAtoms; atom++)
                indexVec[bond*width+atom] = forceAtoms[i][bond][atom];
        }
        OpenCLArray* indices = OpenCLArray::create<cl_uint>(context, indexVec.size(), "bondedIndices");
        indices->upload(indexVec);
        atomIndices.push_back(indices);
        bufferVec[i].resize(width*numBonds, 0);
        for (int bond = 0; bond < numBonds; bond++) {
            for (int atom = 0; atom < numAtoms; atom++)
                bufferVec[i][bond*width+atom] = bufferCounter[i][forceAtoms[i][bond][atom]]++;
        }
        for (int j = 0; j < (int) bufferCounter[i].size(); j++)
            numBuffers[i] = max(numBuffers[i], bufferCounter[i][j]);
    }
    
    // For efficiency, we want to merge multiple forces into a single kernel - but only if that
    // won't increase the number of force buffers.
    
    if (context.getSupports64BitGlobalAtomics()) {
        // Put all the forces in the same set.
        
        numForceBuffers = 1;
        forceSets.push_back(vector<int>());
        for (int i = 0; i < numForces; i++)
            forceSets[0].push_back(i);
    }
    else {
        // Figure out how many force buffers will be required.
    
        for (int i = 0; i < numForces; i++)
            numForceBuffers = max(numForceBuffers, numBuffers[i]);
        int bufferLimit = max(numForceBuffers, (int) context.getPlatformData().contexts.size());
        if (context.getNonbondedUtilities().getHasInteractions())
            bufferLimit = max(bufferLimit, context.getNonbondedUtilities().getNumForceBuffers());
        
        // Figure out sets of forces that can be merged.
        
        vector<int> unmerged(numForces);
        for (int i = 0; i < numForces; i++)
            unmerged[i] = i;
        for (int i = 0; i < numForces; i++)
            for (int j = i-1; j >= 0; j--) {
                if (numBuffers[unmerged[j]] <= numBuffers[unmerged[j+1]])
                    break;
                int temp = unmerged[j+1];
                unmerged[j+1] = unmerged[j];
                unmerged[j] = temp;
            }
        while (unmerged.size() > 0) {
            int sum = numBuffers[unmerged.back()];
            int i;
            for (i = 0; i < (int) unmerged.size()-1; i++) {
                if (sum+numBuffers[unmerged[i]] > bufferLimit)
                    break;
                sum += numBuffers[unmerged[i]];
            }
            forceSets.push_back(vector<int>());
            for (int j = 0; j < i; j++)
                forceSets.back().push_back(unmerged[j]);
            forceSets.back().push_back(unmerged.back());
            for (int j = 0; j < i; j++)
                unmerged.erase(unmerged.begin());
            unmerged.pop_back();
        }
    }

    // Update the buffer indices based on merged sets.
    
    bufferIndices.resize(numForces);
    for (int i = 0; i < (int) forceSets.size(); i++)
        for (int j = 0; j < (int) forceSets[i].size(); j++) {
            int force = forceSets[i][j];
            int numBonds = forceAtoms[force].size();
            int numAtoms = forceAtoms[force][0].size();
            int width = indexWidth[force];
            for (int k = 0; k < j; k++)
                for (int bond = 0; bond < numBonds; bond++)
                    for (int atom = 0; atom < numAtoms; atom++)
                        bufferVec[force][bond*width+atom] += bufferCounter[forceSets[i][k]][forceAtoms[force][bond][atom]];
            OpenCLArray* buffers = OpenCLArray::create<cl_uint>(context, bufferVec[force].size(), "bondedBufferIndices");
            buffers->upload(bufferVec[force]);
            bufferIndices[force] = buffers;
        }

    // Create the kernels.

    for (vector<vector<int> >::const_iterator iter = forceSets.begin(); iter != forceSets.end(); ++iter) {
        const vector<int>& set = *iter;
        int setSize = set.size();
        stringstream s;
        s<<"#ifdef SUPPORTS_64_BIT_ATOMICS\n";
        s<<"#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable\n";
        s<<"#endif\n";
        for (int i = 0; i < (int) prefixCode.size(); i++)
            s<<prefixCode[i];
        string bufferType = (context.getSupports64BitGlobalAtomics() ? "long" : "real4");
        s<<"__kernel void computeBondedForces(__global "<<bufferType<<"* restrict forceBuffers, __global real* restrict energyBuffer, __global const real4* restrict posq, int groups";
        for (int i = 0; i < setSize; i++) {
            int force = set[i];
            string indexType = "uint"+(indexWidth[force] == 1 ? "" : context.intToString(indexWidth[force]));
            s<<", __global const "<<indexType<<"* restrict atomIndices"<<i;
            s<<", __global const "<<indexType<<"* restrict bufferIndices"<<i;
        }
        for (int i = 0; i < (int) arguments.size(); i++)
            s<<", __global "<<argTypes[i]<<"* customArg"<<(i+1);
        s<<") {\n";
        s<<"real energy = 0.0f;\n";
        for (int i = 0; i < setSize; i++) {
            int force = set[i];
            s<<createForceSource(i, forceAtoms[force].size(), forceAtoms[force][0].size(), forceGroup[force], forceSource[force]);
        }
        s<<"energyBuffer[get_global_id(0)] += energy;\n";
        s<<"}\n";
        map<string, string> defines;
        defines["PADDED_NUM_ATOMS"] = context.intToString(context.getPaddedNumAtoms());
        cl::Program program = context.createProgram(s.str(), defines);
        kernels.push_back(cl::Kernel(program, "computeBondedForces"));
    }
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
    s<<"    "<<indexType<<" buffers = bufferIndices"<<forceIndex<<"[index];\n";
    for (int i = 0; i < numAtoms; i++) {
        s<<"    unsigned int atom"<<(i+1)<<" = atoms"<<suffix[i]<<";\n";
        s<<"    real4 pos"<<(i+1)<<" = posq[atom"<<(i+1)<<"];\n";
    }
    s<<computeForce<<"\n";
    for (int i = 0; i < numAtoms; i++) {
        s<<"    {\n";
        if (context.getSupports64BitGlobalAtomics()) {
            s<<"    atom_add(&forceBuffers[atom"<<(i+1)<<"], (long) (force"<<(i+1)<<".x*0x100000000));\n";
            s<<"    atom_add(&forceBuffers[atom"<<(i+1)<<"+PADDED_NUM_ATOMS], (long) (force"<<(i+1)<<".y*0x100000000));\n";
            s<<"    atom_add(&forceBuffers[atom"<<(i+1)<<"+2*PADDED_NUM_ATOMS], (long) (force"<<(i+1)<<".z*0x100000000));\n";
        }
        else {
            s<<"    unsigned int offset = atom"<<(i+1)<<"+buffers"<<suffix[i]<<"*PADDED_NUM_ATOMS;\n";
            s<<"    real4 force = forceBuffers[offset];\n";
            s<<"    force.xyz += force"<<(i+1)<<".xyz;\n";
            s<<"    forceBuffers[offset] = force;\n";
        }
        s<<"    }\n";
    }
    s<<"}\n";
    return s.str();
}

void OpenCLBondedUtilities::computeInteractions(int groups) {
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
        for (int i = 0; i < (int) forceSets.size(); i++) {
            int index = 0;
            cl::Kernel& kernel = kernels[i];
            if (context.getSupports64BitGlobalAtomics())
                kernel.setArg<cl::Buffer>(index++, context.getLongForceBuffer().getDeviceBuffer());
            else
                kernel.setArg<cl::Buffer>(index++, context.getForceBuffers().getDeviceBuffer());
            kernel.setArg<cl::Buffer>(index++, context.getEnergyBuffer().getDeviceBuffer());
            kernel.setArg<cl::Buffer>(index++, context.getPosq().getDeviceBuffer());
            index++;
            for (int j = 0; j < (int) forceSets[i].size(); j++) {
                kernel.setArg<cl::Buffer>(index++, atomIndices[forceSets[i][j]]->getDeviceBuffer());
                kernel.setArg<cl::Buffer>(index++, bufferIndices[forceSets[i][j]]->getDeviceBuffer());
            }
            for (int j = 0; j < (int) arguments.size(); j++)
                kernel.setArg<cl::Memory>(index++, *arguments[j]);
        }
    }
    for (int i = 0; i < (int) kernels.size(); i++) {
        kernels[i].setArg<cl_int>(3, groups);
        context.executeKernel(kernels[i], maxBonds);
    }
}
