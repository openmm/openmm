#ifndef OPENMM_CUDASTREAMIMPL_H_
#define OPENMM_CUDASTREAMIMPL_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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
#include "kernels/gputypes.h"

namespace OpenMM {

/**
 * This is the implementation of streams in the CUDA Platform.
 */

template <class T>
class CudaStreamImpl : public StreamImpl {
public:
    CudaStreamImpl(std::string name, int size, Stream::DataType type, const Platform& platform, int substreams, _gpuContext* gpu);
    CudaStreamImpl(std::string name, int size, Stream::DataType type, const Platform& platform, CUDAStream<T>* stream, int rowOffset, float* padding, _gpuContext* gpu);
    ~CudaStreamImpl();
    void loadFromArray(const void* array);
    void saveToArray(void* array);
    void fillWithValue(void* value);
    const CUDAStream<T>& getStream() const;
    CUDAStream<T>& getStream();
private:
    void initType();
    CUDAStream<T>* stream;
    _gpuContext* gpu;
    bool ownStream;
    int width, rowOffset;
    float paddingValues[4];
    Stream::DataType baseType;
};

template <class T>
CudaStreamImpl<T>::CudaStreamImpl(std::string name, int size, Stream::DataType type, const Platform& platform, int substreams, _gpuContext* gpu) :
        StreamImpl(name, size, type, platform), stream(new CUDAStream<T>(size, substreams)), ownStream(true), gpu(gpu) {
    initType();
    rowOffset = width;
};

template <class T>
CudaStreamImpl<T>::CudaStreamImpl(std::string name, int size, Stream::DataType type, const Platform& platform, CUDAStream<T>* stream, int rowOffset, float* padding, _gpuContext* gpu) :
        StreamImpl(name, size, type, platform), stream(stream), rowOffset(rowOffset), ownStream(false), gpu(gpu) {
    initType();
    for (int i = 0; i < 4; ++i)
        paddingValues[i] = padding[i];
};
    
template <class T>
void CudaStreamImpl<T>::initType() {
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
CudaStreamImpl<T>::~CudaStreamImpl() {
    if (ownStream)
        delete stream;
}

template <class T>
void CudaStreamImpl<T>::loadFromArray(const void* array) {
    float* data = reinterpret_cast<float*>(stream->_pSysData);
    int* order = gpu->psAtomIndex->_pSysData;
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
    for (int i = getSize(); i < (int) stream->_length; ++i)
        for (int j = 0; j < rowOffset; ++j)
            data[i*rowOffset+j] = paddingValues[j];
    stream->Upload();

    // VisualStudio compiler did not like stream == gpu->psPosq4
    //if( gpu && stream == gpu->psPosq4 ){

    if( gpu && getName() == "particlePositions" ){
        gpu->bRecalculateBornRadii = true;
    }
}

template <class T>
void CudaStreamImpl<T>::saveToArray(void* array) {
    stream->Download();
    float* data = reinterpret_cast<float*>(stream->_pSysData);
    int* order = gpu->psAtomIndex->_pSysData;
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
    }
    else {
        int* arrayData = (int*) array;
        for (int i = 0; i < getSize(); ++i)
            for (int j = 0; j < width; ++j)
                arrayData[order[i]*width+j] = (int) data[i*rowOffset+j];
    }
}

template <class T>
void CudaStreamImpl<T>::fillWithValue(void* value) {
    float* data = reinterpret_cast<float*>(stream->_pSysData);
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
    for (int i = getSize(); i < (int) stream->_length; ++i)
        for (int j = 0; j < rowOffset; ++j)
            data[i*rowOffset+j] = paddingValues[j];
    stream->Upload();
}

template <class T>
const CUDAStream<T>& CudaStreamImpl<T>::getStream() const {
    return stream;
}

template <class T>
CUDAStream<T>& CudaStreamImpl<T>::getStream() {
    return stream;
}

} // namespace OpenMM

#endif /*OPENMM_CUDASTREAMIMPL_H_*/
