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
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "StreamImpl.h"
#include "kernels/cudatypes.h"

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
    if (baseType == Stream::Float) {
        float* arrayData = (float*) array;
        for (int i = 0; i < getSize(); ++i)
            for (int j = 0; j < width; ++j)
                data[i*rowOffset+j] = arrayData[i*width+j];
    }
    else if (baseType == Stream::Double) {
        double* arrayData = (double*) array;
        for (int i = 0; i < getSize(); ++i)
            for (int j = 0; j < width; ++j)
                data[i*rowOffset+j] = (float) arrayData[i*width+j];
    }
    else {
        int* arrayData = (int*) array;
        for (int i = 0; i < getSize(); ++i)
            for (int j = 0; j < width; ++j)
                data[i*rowOffset+j] = (float) arrayData[i*width+j];
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
    if (baseType == Stream::Float) {
        float* arrayData = (float*) array;
        for (int i = 0; i < getSize(); ++i)
            for (int j = 0; j < width; ++j)
                arrayData[i*width+j] = data[i*rowOffset+j];
    }
    else if (baseType == Stream::Double) {
        double* arrayData = (double*) array;
        for (int i = 0; i < getSize(); ++i)
            for (int j = 0; j < width; ++j)
                arrayData[i*width+j] = data[i*rowOffset+j];
    }
    else {
        int* arrayData = (int*) array;
        for (int i = 0; i < getSize(); ++i)
            for (int j = 0; j < width; ++j)
                arrayData[i*width+j] = (int) data[i*rowOffset+j];
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
