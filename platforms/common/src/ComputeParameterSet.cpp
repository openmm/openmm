/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2024 Stanford University and the Authors.      *
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

#include "openmm/common/ComputeParameterSet.h"
#include "openmm/OpenMMException.h"
#include <cmath>
#include <sstream>

using namespace OpenMM;
using namespace std;

ComputeParameterSet::ComputeParameterSet(ComputeContext& context, int numParameters, int numObjects, const string& name, bool arrayPerParameter, bool useDoublePrecision) :
            context(context), numParameters(numParameters), numObjects(numObjects), name(name) {
    int params = numParameters;
    int bufferCount = 0;
    elementSize = (useDoublePrecision ? sizeof(double) : sizeof(float));
    string elementType = (useDoublePrecision ? "double" : "float");
    if (!arrayPerParameter) {
        while (params > 2) {
            std::stringstream arrayName;
            arrayName << "param" << (++bufferCount);
            arrays.push_back(context.createArray());
            arrays.back()->initialize(context, numObjects, elementSize*4, arrayName.str());
            params -= 4;
        }
        if (params > 1) {
            std::stringstream arrayName;
            arrayName << "param" << (++bufferCount);
            arrays.push_back(context.createArray());
            arrays.back()->initialize(context, numObjects, elementSize*2, arrayName.str());
            params -= 2;
        }
    }
    while (params > 0) {
        std::stringstream arrayName;
        arrayName << "param" << (++bufferCount);
            arrays.push_back(context.createArray());
        arrays.back()->initialize(context, numObjects, elementSize, arrayName.str());
        params--;
    }
    for (ArrayInterface* array : arrays)
        parameters.push_back(ComputeParameterInfo(*array, array->getName(), elementType, array->getElementSize()/elementSize));
}

ComputeParameterSet::~ComputeParameterSet() {
    for (ArrayInterface* array : arrays)
        delete array;
}

template <class T>
void ComputeParameterSet::getParameterValues(vector<vector<T> >& values) {
    if (sizeof(T) != elementSize)
        throw OpenMMException("Called getParameterValues() with vector of wrong type");
    values.resize(numObjects);
    for (int i = 0; i < numObjects; i++)
        values[i].resize(numParameters);
    int base = 0;
    for (int i = 0; i < (int) arrays.size(); i++) {
        if (arrays[i]->getElementSize() == 4*elementSize) {
            vector<T> data(4*numObjects);
            arrays[i]->download(data.data());
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
        else if (arrays[i]->getElementSize() == 2*elementSize) {
            vector<T> data(2*numObjects);
            arrays[i]->download(data.data());
            for (int j = 0; j < numObjects; j++) {
                values[j][base] = data[2*j];
                if (base+1 < numParameters)
                    values[j][base+1] = data[2*j+1];
            }
            base += 2;
        }
        else if (arrays[i]->getElementSize() == elementSize) {
            vector<T> data(numObjects);
            arrays[i]->download(data.data());
            for (int j = 0; j < numObjects; j++)
                values[j][base] = data[j];
            base++;
        }
        else
            throw OpenMMException("Internal error: Unknown buffer type in ComputeParameterSet");
    }
}

template <class T>
void ComputeParameterSet::setParameterValues(const vector<vector<T> >& values, bool convert) {
    setParameterValuesSubset(0, values, convert);
}

template <class T>
void ComputeParameterSet::setParameterValuesSubset(int first, const vector<vector<T> >& values, bool convert) {
    if (convert && sizeof(T) == sizeof(float) && elementSize == sizeof(double)) {
        vector<vector<double> > values2(values.size());
        for (int i = 0; i < values.size(); i++) {
            values2[i].resize(values[i].size(), 0.0);
            for (int j = 0; j < values[i].size(); j++)
                values2[i][j] = (double) values[i][j];
        }
        setParameterValuesSubset(first, values2);
        return;
    }
    if (convert && sizeof(T) == sizeof(double) && elementSize == sizeof(float)) {
        vector<vector<float> > values2(values.size());
        for (int i = 0; i < values.size(); i++) {
            values2[i].resize(values[i].size(), 0.0);
            for (int j = 0; j < values[i].size(); j++)
                values2[i][j] = (float) values[i][j];
        }
        setParameterValuesSubset(first, values2);
        return;
    }
    if (sizeof(T) != elementSize)
        throw OpenMMException("Called setParameterValues() with vector of wrong type");
    if (first < 0 || first+values.size() > numObjects)
        throw OpenMMException("Called setParameterValuesSubset() for an illegal range of objects");
    int numToSet = values.size();
    int base = 0;
    for (int i = 0; i < (int) arrays.size(); i++) {
        if (arrays[i]->getElementSize() == 4*elementSize) {
            vector<T> data(4*numToSet);
            for (int j = 0; j < numToSet; j++) {
                data[4*j] = values[j][base];
                if (base+1 < numParameters)
                    data[4*j+1] = values[j][base+1];
                if (base+2 < numParameters)
                    data[4*j+2] = values[j][base+2];
                if (base+3 < numParameters)
                    data[4*j+3] = values[j][base+3];
            }
            arrays[i]->uploadSubArray(data.data(), first, numToSet);
            base += 4;
        }
        else if (arrays[i]->getElementSize() == 2*elementSize) {
            vector<T> data(2*numToSet);
            for (int j = 0; j < numToSet; j++) {
                data[2*j] = values[j][base];
                if (base+1 < numParameters)
                    data[2*j+1] = values[j][base+1];
            }
            arrays[i]->uploadSubArray(data.data(), first, numToSet);
            base += 2;
        }
        else if (arrays[i]->getElementSize() == elementSize) {
            vector<T> data(numToSet);
            for (int j = 0; j < numToSet; j++)
                data[j] = values[j][base];
            arrays[i]->uploadSubArray(data.data(), first, numToSet);
            base++;
        }
        else
            throw OpenMMException("Internal error: Unknown buffer type in ComputeParameterSet");
    }
}

string ComputeParameterSet::getParameterSuffix(int index, const std::string& extraSuffix) const {
    const string suffixes[] = {".x", ".y", ".z", ".w"};
    int buffer = -1;
    for (int i = 0; buffer == -1 && i < (int) parameters.size(); i++) {
        if (index*elementSize < parameters[i].getSize())
            buffer = i;
        else
            index -= parameters[i].getSize()/elementSize;
    }
    if (buffer == -1)
        throw OpenMMException("Internal error: Illegal argument to ComputeParameterSet::getParameterSuffix() ("+name+")");
    stringstream suffix;
    suffix << (buffer+1) << extraSuffix;
    if (parameters[buffer].getSize() != elementSize)
        suffix << suffixes[index];
    return suffix.str();
}

/**
 * Define template instantiations for float and double versions of getParameterValues() and setParameterValues().
 */
namespace OpenMM {
template void ComputeParameterSet::getParameterValues<float>(vector<vector<float> >& values);
template void ComputeParameterSet::setParameterValues<float>(const vector<vector<float> >& values, bool convert);
template void ComputeParameterSet::setParameterValuesSubset<float>(int first, const vector<vector<float> >& values, bool convert);
template void ComputeParameterSet::getParameterValues<double>(vector<vector<double> >& values);
template void ComputeParameterSet::setParameterValues<double>(const vector<vector<double> >& values, bool convert);
template void ComputeParameterSet::setParameterValuesSubset<double>(int first, const vector<vector<double> >& values, bool convert);
}