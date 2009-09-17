#ifndef OPENMM_OPENCLSTREAMIMPL_H_
#define OPENMM_OPENCLSTREAMIMPL_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2009 Stanford University and the Authors.      *
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

#include "openmm/StreamImpl.h"
#include "OpenCLArray.h"

namespace OpenMM {

/**
 * This is the implementation of streams in the OpenCL Platform.
 */

template <class T>
class OpenCLStreamImpl : public StreamImpl {
public:
    OpenCLStreamImpl(std::string name, int size, Stream::DataType type, const Platform& platform, OpenCLContext& context);
    OpenCLStreamImpl(std::string name, int size, Stream::DataType type, const Platform& platform, OpenCLArray<T>* clarray, int rowOffset, OpenCLContext& context);
    ~OpenCLStreamImpl();
    void loadFromArray(const void* array);
    void saveToArray(void* array);
    void fillWithValue(void* value);
    const OpenCLArray<T>& getArray() const;
    OpenCLArray<T>& getArray();
private:
    void initType();
    OpenCLArray<T>* clarray;
    OpenCLContext& context;
    bool ownArray;
    int width, rowOffset;
    Stream::DataType baseType;
};

template <class T>
OpenCLStreamImpl<T>::OpenCLStreamImpl(std::string name, int size, Stream::DataType type, const Platform& platform, OpenCLContext& context) :
        StreamImpl(name, size, type, platform), ownArray(true), context(context) {
    initType();
    clarray = new OpenCLArray<T>(context, size*width, name, true);
    rowOffset = width;
};

template <class T>
OpenCLStreamImpl<T>::OpenCLStreamImpl(std::string name, int size, Stream::DataType type, const Platform& platform, OpenCLArray<T>* clarray, int rowOffset, OpenCLContext& context) :
        StreamImpl(name, size, type, platform), clarray(clarray), rowOffset(rowOffset), ownArray(false), context(context) {
    initType();
};

template <class T>
void OpenCLStreamImpl<T>::initType() {
    switch (getDataType()) {
    case Stream::Float:
    case Stream::Float2:
    case Stream::Float3:
    case Stream::Float4:
        baseType = Stream::Float;
        break;
    case Stream::Double:
    case Stream::Double2:
    case Stream::Double3:
    case Stream::Double4:
        baseType = Stream::Double;
        break;
    case Stream::Integer:
    case Stream::Integer2:
    case Stream::Integer3:
    case Stream::Integer4:
        baseType = Stream::Integer;
        break;
    }
    switch (getDataType()) {
    case Stream::Float:
    case Stream::Double:
    case Stream::Integer:
        width = 1;
        break;
    case Stream::Float2:
    case Stream::Double2:
    case Stream::Integer2:
        width = 2;
        break;
    case Stream::Float3:
    case Stream::Double3:
    case Stream::Integer3:
        width = 3;
        break;
    case Stream::Float4:
    case Stream::Double4:
    case Stream::Integer4:
        width = 4;
        break;
    }
}

template <class T>
OpenCLStreamImpl<T>::~OpenCLStreamImpl() {
    if (ownArray)
        delete clarray;
}

template <class T>
void OpenCLStreamImpl<T>::loadFromArray(const void* array) {
    float* data = reinterpret_cast<float*>(clarray->getHostBuffer());
    OpenCLArray<cl_int>& order = context.getAtomIndex();
    if (baseType == Stream::Float) {
        float* arrayData = (float*) array;
        for (int i = 0; i < getSize(); ++i)
            for (int j = 0; j < width; ++j)
                data[i*rowOffset+j] = arrayData[order[i]*width+j];
    }
    else if (baseType == Stream::Double) {
        double* arrayData = (double*) array;
        for (int i = 0; i < getSize(); ++i)
            for (int j = 0; j < width; ++j)
                data[i*rowOffset+j] = (float) arrayData[order[i]*width+j];
    }
    else {
        int* arrayData = (int*) array;
        for (int i = 0; i < getSize(); ++i)
            for (int j = 0; j < width; ++j)
                data[i*rowOffset+j] = (float) arrayData[order[i]*width+j];
    }
    clarray->upload();
//    if (getName() == "particlePositions") {
//        for (int i = 0; i < context->posCellOffsets.size(); i++)
//            context.posCellOffsets[i] = make_int3(0, 0, 0);
//    }
}

template <class T>
void OpenCLStreamImpl<T>::saveToArray(void* array) {
    clarray->download();
    float* data = reinterpret_cast<float*>(clarray->getHostBuffer());
    OpenCLArray<cl_int>& order = context.getAtomIndex();
    if (baseType == Stream::Float) {
        float* arrayData = (float*) array;
        for (int i = 0; i < getSize(); ++i)
            for (int j = 0; j < width; ++j)
                arrayData[order[i]*width+j] = data[i*rowOffset+j];
    }
    else if (baseType == Stream::Double) {
        double* arrayData = (double*) array;
        for (int i = 0; i < getSize(); ++i)
            for (int j = 0; j < width; ++j)
                arrayData[order[i]*width+j] = data[i*rowOffset+j];
//        if (context && getName() == "particlePositions") {
//            for (int i = 0; i < getSize(); i++) {
//                int3 offset = context->posCellOffsets[i];
//                arrayData[order[i]*width] -= offset.x*context->sim.periodicBoxSizeX;
//                arrayData[order[i]*width+1] -= offset.y*context->sim.periodicBoxSizeY;
//                arrayData[order[i]*width+2] -= offset.z*context->sim.periodicBoxSizeZ;
//            }
//        }
    }
    else {
        int* arrayData = (int*) array;
        for (int i = 0; i < getSize(); ++i)
            for (int j = 0; j < width; ++j)
                arrayData[order[i]*width+j] = (int) data[i*rowOffset+j];
    }
}

template <class T>
void OpenCLStreamImpl<T>::fillWithValue(void* value) {
    float* data = reinterpret_cast<float*>(clarray->getHostBuffer());
    if (baseType == Stream::Float) {
        float valueData = *((float*) value);
        for (int i = 0; i < getSize(); ++i)
            for (int j = 0; j < width; ++j)
                data[i*rowOffset+j] = valueData;
    }
    else if (baseType == Stream::Double) {
        double valueData = *((double*) value);
        for (int i = 0; i < getSize(); ++i)
            for (int j = 0; j < width; ++j)
                data[i*rowOffset+j] = (float) valueData;
    }
    else {
        int valueData = *((int*) value);
        for (int i = 0; i < getSize(); ++i)
            for (int j = 0; j < width; ++j)
                data[i*rowOffset+j] = (float) valueData;
    }
    clarray->upload();
//    if (context && getName() == "particlePositions") {
//        context->bRecalculateBornRadii = true;
//        for (int i = 0; i < context->posCellOffsets.size(); i++)
//            context->posCellOffsets[i] = make_int3(0, 0, 0);
//    }
}

template <class T>
const OpenCLArray<T>& OpenCLStreamImpl<T>::getArray() const {
    return clarray;
}

template <class T>
OpenCLArray<T>& OpenCLStreamImpl<T>::getArray() {
    return clarray;
}

} // namespace OpenMM

#endif /*OPENMM_OPENCLSTREAMIMPL_H_*/
