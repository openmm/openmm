/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2012 Stanford University and the Authors.      *
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

#include "CudaParameterSet.h"
#include "openmm/OpenMMException.h"
#include <cmath>
#include <sstream>

using namespace OpenMM;
using namespace std;

#define CHECK_RESULT(result) \
    if (result != CUDA_SUCCESS) { \
        std::stringstream m; \
        m<<errorMessage<<": "<<context.getErrorString(result)<<" ("<<result<<")"; \
        throw OpenMMException(m.str());\
    }

CudaParameterSet::CudaParameterSet(CudaContext& context, int numParameters, int numObjects, const string& name, bool bufferPerParameter, bool useDoublePrecision) :
            context(context), numParameters(numParameters), numObjects(numObjects), name(name) {
    int params = numParameters;
    int bufferCount = 0;
    elementSize = (useDoublePrecision ? sizeof(double) : sizeof(float));
    string elementType = (useDoublePrecision ? "double" : "float");
    CUdeviceptr pointer;
    string errorMessage = "Error creating parameter set "+name;
    if (!bufferPerParameter) {
        while (params > 2) {
            CHECK_RESULT(cuMemAlloc(&pointer, numObjects*elementSize*4));
            std::stringstream name;
            name << "param" << (++bufferCount);
            buffers.push_back(CudaNonbondedUtilities::ParameterInfo(name.str(), elementType, 4, elementSize*4, pointer));
            params -= 4;
        }
        if (params > 1) {
            CHECK_RESULT(cuMemAlloc(&pointer, numObjects*elementSize*2));
            std::stringstream name;
            name << "param" << (++bufferCount);
            buffers.push_back(CudaNonbondedUtilities::ParameterInfo(name.str(), elementType, 2, elementSize*2, pointer));
            params -= 2;
        }
    }
    while (params > 0) {
            CHECK_RESULT(cuMemAlloc(&pointer, numObjects*elementSize));
        std::stringstream name;
        name << "param" << (++bufferCount);
        buffers.push_back(CudaNonbondedUtilities::ParameterInfo(name.str(), elementType, 1, elementSize, pointer));
        params--;
    }
}

CudaParameterSet::~CudaParameterSet() {
    if (context.getContextIsValid()) {
        string errorMessage = "Error freeing device memory";
        for (int i = 0; i < (int) buffers.size(); i++)
            CHECK_RESULT(cuMemFree(buffers[i].getMemory()));
    }
}

template <class T>
void CudaParameterSet::getParameterValues(vector<vector<T> >& values) {
    if (sizeof(T) != elementSize)
        throw OpenMMException("Called getParameterValues() with vector of wrong type");
    values.resize(numObjects);
    for (int i = 0; i < numObjects; i++)
        values[i].resize(numParameters);
    int base = 0;
    string errorMessage = "Error downloading parameter set "+name;
    for (int i = 0; i < (int) buffers.size(); i++) {
        if (buffers[i].getSize() == 4*elementSize) {
            vector<T> data(4*numObjects);
            CHECK_RESULT(cuMemcpyDtoH(&data[0], buffers[i].getMemory(), numObjects*buffers[i].getSize()));
            for (int j = 0; j < numObjects; j++) {
                values[j][base] = data[4*j];
                if (base+1 < numParameters)
                    values[j][base+1] = data[4*j+1];
                if (base+2 < numParameters)
                    values[j][base+2] = data[4*j+2];
                if (base+3 < numParameters)
                    values[j][base+3] = data[4*j+3];
            }
            base += 4;
        }
        else if (buffers[i].getSize() == 2*elementSize) {
            vector<T> data(2*numObjects);
            CHECK_RESULT(cuMemcpyDtoH(&data[0], buffers[i].getMemory(), numObjects*buffers[i].getSize()));
            for (int j = 0; j < numObjects; j++) {
                values[j][base] = data[2*j];
                if (base+1 < numParameters)
                    values[j][base+1] = data[2*j+1];
            }
            base += 2;
        }
        else if (buffers[i].getSize() == elementSize) {
            vector<T> data(numObjects);
            CHECK_RESULT(cuMemcpyDtoH(&data[0], buffers[i].getMemory(), numObjects*buffers[i].getSize()));
            for (int j = 0; j < numObjects; j++)
                values[j][base] = data[j];
            base++;
        }
        else
            throw OpenMMException("Internal error: Unknown buffer type in CudaParameterSet");
    }
}

template <class T>
void CudaParameterSet::setParameterValues(const vector<vector<T> >& values) {
    if (sizeof(T) != elementSize)
        throw OpenMMException("Called setParameterValues() with vector of wrong type");
    int base = 0;
    string errorMessage = "Error uploading parameter set "+name;
    for (int i = 0; i < (int) buffers.size(); i++) {
        if (buffers[i].getSize() == 4*elementSize) {
            vector<T> data(4*numObjects);
            for (int j = 0; j < numObjects; j++) {
                data[4*j] = values[j][base];
                if (base+1 < numParameters)
                    data[4*j+1] = values[j][base+1];
                if (base+2 < numParameters)
                    data[4*j+2] = values[j][base+2];
                if (base+3 < numParameters)
                    data[4*j+3] = values[j][base+3];
            }
            CHECK_RESULT(cuMemcpyHtoD(buffers[i].getMemory(), &data[0], numObjects*buffers[i].getSize()));
            base += 4;
        }
        else if (buffers[i].getSize() == 2*elementSize) {
            vector<T> data(2*numObjects);
            for (int j = 0; j < numObjects; j++) {
                data[2*j] = values[j][base];
                if (base+1 < numParameters)
                    data[2*j+1] = values[j][base+1];
            }
            CHECK_RESULT(cuMemcpyHtoD(buffers[i].getMemory(), &data[0], numObjects*buffers[i].getSize()));
            base += 2;
        }
        else if (buffers[i].getSize() == elementSize) {
            vector<T> data(numObjects);
            for (int j = 0; j < numObjects; j++)
                data[j] = values[j][base];
            CHECK_RESULT(cuMemcpyHtoD(buffers[i].getMemory(), &data[0], numObjects*buffers[i].getSize()));
            base++;
        }
        else
            throw OpenMMException("Internal error: Unknown buffer type in CudaParameterSet");
    }
}

string CudaParameterSet::getParameterSuffix(int index, const std::string& extraSuffix) const {
    const string suffixes[] = {".x", ".y", ".z", ".w"};
    int buffer = -1;
    for (int i = 0; buffer == -1 && i < (int) buffers.size(); i++) {
        if (index*elementSize < buffers[i].getSize())
            buffer = i;
        else
            index -= buffers[i].getSize()/elementSize;
    }
    if (buffer == -1)
        throw OpenMMException("Internal error: Illegal argument to CudaParameterSet::getParameterSuffix() ("+name+")");
    stringstream suffix;
    suffix << (buffer+1) << extraSuffix;
    if (buffers[buffer].getSize() != elementSize)
        suffix << suffixes[index];
    return suffix.str();
}

/**
 * Define template instantiations for float and double versions of getParameterValues() and setParameterValues().
 */
namespace OpenMM {
template OPENMM_EXPORT_CUDA void CudaParameterSet::getParameterValues<float>(vector<vector<float> >& values);
template OPENMM_EXPORT_CUDA void CudaParameterSet::setParameterValues<float>(const vector<vector<float> >& values);
template OPENMM_EXPORT_CUDA void CudaParameterSet::getParameterValues<double>(vector<vector<double> >& values);
template OPENMM_EXPORT_CUDA void CudaParameterSet::setParameterValues<double>(const vector<vector<double> >& values);
}