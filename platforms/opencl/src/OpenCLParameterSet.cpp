/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
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

#include "OpenCLParameterSet.h"
#include "openmm/OpenMMException.h"
#include <cmath>
#include <sstream>

using namespace OpenMM;
using namespace std;

OpenCLParameterSet::OpenCLParameterSet(OpenCLContext& context, int numParameters, int numObjects, const string& name) :
            context(context), numParameters(numParameters), numObjects(numObjects), name(name) {
    int params = numParameters;
    int bufferCount = 0;
    try {
        while (params > 2) {
            cl::Buffer* buf = new cl::Buffer(context.getContext(), CL_MEM_READ_WRITE, numObjects*sizeof(mm_float4));
            std::stringstream name;
            name << "param" << (++bufferCount);
            buffers.push_back(OpenCLNonbondedUtilities::ParameterInfo(name.str(), "float", 4, sizeof(mm_float4), *buf));
            params -= 4;
        }
        if (params > 1) {
            cl::Buffer* buf = new cl::Buffer(context.getContext(), CL_MEM_READ_WRITE, numObjects*sizeof(mm_float2));
            std::stringstream name;
            name << "param" << (++bufferCount);
            buffers.push_back(OpenCLNonbondedUtilities::ParameterInfo(name.str(), "float", 2, sizeof(mm_float2), *buf));
            params -= 2;
        }
        if (params > 0) {
            cl::Buffer* buf = new cl::Buffer(context.getContext(), CL_MEM_READ_WRITE, numObjects*sizeof(cl_float));
            std::stringstream name;
            name << "param" << (++bufferCount);
            buffers.push_back(OpenCLNonbondedUtilities::ParameterInfo(name.str(), "float", 1, sizeof(cl_float), *buf));
        }
    }
    catch (cl::Error err) {
        stringstream str;
        str<<"Error creating parameter set "<<name<<": "<<err.what()<<" ("<<err.err()<<")";
        throw OpenMMException(str.str());
    }
}

OpenCLParameterSet::~OpenCLParameterSet() {
    for (int i = 0; i < (int) buffers.size(); i++)
        delete &buffers[i].getMemory();
}

void OpenCLParameterSet::getParameterValues(vector<vector<cl_float> >& values) const {
    values.resize(numObjects);
    for (int i = 0; i < numObjects; i++)
        values[i].resize(numParameters);
    try {
        int base = 0;
        for (int i = 0; i < (int) buffers.size(); i++) {
            if (buffers[i].getType() == "float4") {
                vector<mm_float4> data(numObjects);
                context.getQueue().enqueueReadBuffer(reinterpret_cast<cl::Buffer&>(buffers[i].getMemory()), CL_TRUE, 0, numObjects*buffers[i].getSize(), &data[0]);
                for (int j = 0; j < numObjects; j++) {
                    values[j][base] = data[j].x;
                    if (base+1 < numParameters)
                        values[j][base+1] = data[j].y;
                    if (base+2 < numParameters)
                        values[j][base+2] = data[j].z;
                    if (base+3 < numParameters)
                        values[j][base+3] = data[j].w;
                }
                base += 4;
            }
            else if (buffers[i].getType() == "float2") {
                vector<mm_float2> data(numObjects);
                context.getQueue().enqueueReadBuffer(reinterpret_cast<cl::Buffer&>(buffers[i].getMemory()), CL_TRUE, 0, numObjects*buffers[i].getSize(), &data[0]);
                for (int j = 0; j < numObjects; j++) {
                    values[j][base] = data[j].x;
                    if (base+1 < numParameters)
                        values[j][base+1] = data[j].y;
                }
                base += 2;
            }
            else if (buffers[i].getType() == "float") {
                vector<cl_float> data(numObjects);
                context.getQueue().enqueueReadBuffer(reinterpret_cast<cl::Buffer&>(buffers[i].getMemory()), CL_TRUE, 0, numObjects*buffers[i].getSize(), &data[0]);
                for (int j = 0; j < numObjects; j++)
                    values[j][base] = data[j];
            }
            else
                throw OpenMMException("Internal error: Unknown buffer type in OpenCLParameterSet");
        }
    }
    catch (cl::Error err) {
        stringstream str;
        str<<"Error downloading parameter set "<<name<<": "<<err.what()<<" ("<<err.err()<<")";
        throw OpenMMException(str.str());
    }
}

void OpenCLParameterSet::setParameterValues(const vector<vector<cl_float> >& values) {
    try {
        int base = 0;
        for (int i = 0; i < (int) buffers.size(); i++) {
            if (buffers[i].getType() == "float4") {
                vector<mm_float4> data(numObjects);
                for (int j = 0; j < numObjects; j++) {
                    data[j].x = values[j][base];
                    if (base+1 < numParameters)
                        data[j].y = values[j][base+1];
                    if (base+2 < numParameters)
                        data[j].z = values[j][base+2];
                    if (base+3 < numParameters)
                        data[j].w = values[j][base+3];
                }
                context.getQueue().enqueueWriteBuffer(reinterpret_cast<cl::Buffer&>(buffers[i].getMemory()), CL_TRUE, 0, numObjects*buffers[i].getSize(), &data[0]);
                base += 4;
            }
            else if (buffers[i].getType() == "float2") {
                vector<mm_float2> data(numObjects);
                for (int j = 0; j < numObjects; j++) {
                    data[j].x = values[j][base];
                    if (base+1 < numParameters)
                        data[j].y = values[j][base+1];
                }
                context.getQueue().enqueueWriteBuffer(reinterpret_cast<cl::Buffer&>(buffers[i].getMemory()), CL_TRUE, 0, numObjects*buffers[i].getSize(), &data[0]);
                base += 2;
            }
            else if (buffers[i].getType() == "float") {
                vector<cl_float> data(numObjects);
                for (int j = 0; j < numObjects; j++)
                    data[j] = values[j][base];
                context.getQueue().enqueueWriteBuffer(reinterpret_cast<cl::Buffer&>(buffers[i].getMemory()), CL_TRUE, 0, numObjects*buffers[i].getSize(), &data[0]);
            }
            else
                throw OpenMMException("Internal error: Unknown buffer type in OpenCLParameterSet");
        }
    }
    catch (cl::Error err) {
        stringstream str;
        str<<"Error uploading parameter set "<<name<<": "<<err.what()<<" ("<<err.err()<<")";
        throw OpenMMException(str.str());
    }
}

string OpenCLParameterSet::getParameterSuffix(int index, const std::string& extraSuffix) const {
    const string suffixes[] = {".x", ".y", ".z", ".w"};
    int buffer = -1;
    for (int i = 0; buffer == -1 && i < (int) buffers.size(); i++) {
        if (index*sizeof(cl_float) < buffers[i].getSize())
            buffer = i;
        else
            index -= buffers[i].getSize()/sizeof(cl_float);
    }
    if (buffer == -1)
        throw OpenMMException("Internal error: Illegal argument to OpenCLParameterSet::getParameterSuffix() ("+name+")");
    stringstream suffix;
    suffix << (buffer+1) << extraSuffix;
    if (buffers[buffer].getType() != "float")
        suffix << suffixes[index];
    return suffix.str();
}
