/*******************************************************************************
 * Copyright (c) 2008-2012 The Khronos Group Inc.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and/or associated documentation files (the
 * "Materials"), to deal in the Materials without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Materials, and to
 * permit persons to whom the Materials are furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Materials.
 *
 * THE MATERIALS ARE PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * MATERIALS OR THE USE OR OTHER DEALINGS IN THE MATERIALS.
 ******************************************************************************/

/*! \file
 *
 *   \brief C++ bindings for OpenCL 1.0 (rev 48), OpenCL 1.1 (rev 33) and 
 *       OpenCL 1.2 (rev 15)    
 *   \author Benedict R. Gaster and Laurent Morichetti
 *   
 *   Additions and fixes from:
 *       Brian Cole, March 3rd 2010 
 *       Lee Howes, October 2011, March 2012.
 *   
 *   \version 1.2
 *   \date April 2012
 *
 *   Optional extension support
 *
 *         cl
 *         cl_ext_device_fission
 *				#define USE_CL_DEVICE_FISSION
 */

/*! \mainpage
 * \section intro Introduction
 * For many large applications C++ is the language of choice and so it seems
 * reasonable to define C++ bindings for OpenCL.
 *
 *
 * The interface is contained with a single C++ header file \em cl.hpp and all
 * definitions are contained within the namespace \em cl. There is no additional
 * requirement to include \em cl.h and to use either the C++ or original C
 * bindings it is enough to simply include \em cl.hpp.
 *
 * The bindings themselves are lightweight and correspond closely to the
 * underlying C API. Using the C++ bindings introduces no additional execution
 * overhead.
 *
 * For detail documentation on the bindings see:
 *
 * The OpenCL C++ Wrapper API 1.1 (revision 04)
 *  http://www.khronos.org/registry/cl/specs/opencl-cplusplus-1.1.pdf
 *
 * \section example Example
 *
 * The following example shows a general use case for the C++
 * bindings, including support for the optional exception feature and
 * also the supplied vector and string classes, see following sections for
 * decriptions of these features.
 *
 * \code
 * #define __CL_ENABLE_EXCEPTIONS
 * 
 * #if defined(__APPLE__) || defined(__MACOSX)
 * #include <OpenCL/cl.hpp>
 * #else
 * #include <CL/cl.hpp>
 * #endif
 * #include <cstdio>
 * #include <cstdlib>
 * #include <iostream>
 * 
 *  const char * helloStr  = "__kernel void "
 *                           "hello(void) "
 *                           "{ "
 *                           "  "
 *                           "} ";
 * 
 *  int
 *  main(void)
 *  {
 *     cl_int err = CL_SUCCESS;
 *     try {
 *
 *       std::vector<cl::Platform> platforms;
 *       cl::Platform::get(&platforms);
 *       if (platforms.size() == 0) {
 *           std::cout << "Platform size 0\n";
 *           return -1;
 *       }
 *
 *       cl_context_properties properties[] = 
 *          { CL_CONTEXT_PLATFORM, (cl_context_properties)(platforms[0])(), 0};
 *       cl::Context context(CL_DEVICE_TYPE_CPU, properties); 
 * 
 *       std::vector<cl::Device> devices = context.getInfo<CL_CONTEXT_DEVICES>();
 * 
 *       cl::Program::Sources source(1,
 *           std::make_pair(helloStr,strlen(helloStr)));
 *       cl::Program program_ = cl::Program(context, source);
 *       program_.build(devices);
 * 
 *       cl::Kernel kernel(program_, "hello", &err);
 * 
 *       cl::Event event;
 *       cl::CommandQueue queue(context, devices[0], 0, &err);
 *       queue.enqueueNDRangeKernel(
 *           kernel, 
 *           cl::NullRange, 
 *           cl::NDRange(4,4),
 *           cl::NullRange,
 *           NULL,
 *           &event); 
 * 
 *       event.wait();
 *     }
 *     catch (cl::Error err) {
 *        std::cerr 
 *           << "ERROR: "
 *           << err.what()
 *           << "("
 *           << err.err()
 *           << ")"
 *           << std::endl;
 *     }
 * 
 *    return EXIT_SUCCESS;
 *  }
 * 
 * \endcode
 *
 */
#ifndef CL_HPP_
#define CL_HPP_

#ifdef _WIN32
 
#include <windows.h>
#include <malloc.h>
#include <iterator>

#if defined(__CL_ENABLE_EXCEPTIONS)
#include <exception>
#endif // #if defined(__CL_ENABLE_EXCEPTIONS)

#pragma push_macro("max")
#undef max
#if defined(USE_DX_INTEROP)
#include <CL/cl_d3d10.h>
#endif
#endif // _WIN32

// 
#if defined(USE_CL_DEVICE_FISSION)
#include <CL/cl_ext.h>
#endif

#define CL_USE_DEPRECATED_OPENCL_1_1_APIS 
#if defined(__APPLE__) || defined(__MACOSX)
#include <OpenGL/OpenGL.h>
#include <OpenCL/opencl.h>
#include <libkern/OSAtomic.h>
#else
#include <GL/gl.h>
#include <CL/opencl.h>
#endif // !__APPLE__

#if !defined(CL_CALLBACK)
#define CL_CALLBACK
#endif //CL_CALLBACK

#include <utility>
#include <limits>

#if !defined(__NO_STD_VECTOR)
#include <vector>
#endif

#if !defined(__NO_STD_STRING)
#include <string>
#endif 

#if defined(linux) || defined(__APPLE__) || defined(__MACOSX)
#include <alloca.h>
#endif // linux

#include <cstring>

#if (_MSC_VER >= 1600)
#include <functional>
#define __CL_FUNCTION_TYPE typename std::function
#define CL_USE_CPP_FUNCTORS
#elif (__GNUC__ == 4 && __GNUC_MINOR__ >= 1 && !defined(STLPORT)) || defined(__APPLE__) || defined(__MACOSX)
#include <tr1/functional>
#define __CL_FUNCTION_TYPE typename std::tr1::function
#define CL_USE_CPP_FUNCTORS
#endif

/*! \namespace cl
 *
 * \brief The OpenCL C++ bindings are defined within this namespace.
 *
 */
namespace cl {

/**
 * Deprecated APIs for 1.2
 */
#if defined(CL_VERSION_1_1)
#define __INIT_CL_EXT_FCN_PTR(name) \
    if(!pfn_##name) { \
        pfn_##name = (PFN_##name) \
            clGetExtensionFunctionAddress(#name); \
        if(!pfn_##name) { \
        } \
    }
#endif // #if defined(CL_VERSION_1_1)

#if defined(CL_VERSION_1_2)
#define __INIT_CL_EXT_FCN_PTR_PLATFORM(platform, name) \
    if(!pfn_##name) { \
        pfn_##name = (PFN_##name) \
            clGetExtensionFunctionAddressForPlatform(platform, #name); \
        if(!pfn_##name) { \
        } \
    }
#endif // #if defined(CL_VERSION_1_1)

class Program;
class Device;
class Context;
class CommandQueue;
class Memory;

#if defined(__CL_ENABLE_EXCEPTIONS)
/*! \class Error
 * \brief Exception class
 */
class Error : public std::exception
{
private:
    cl_int err_;
    const char * errStr_;
public:
    /*! Create a new CL error exception for a given error code
     *  and corresponding message.
     */
    Error(cl_int err, const char * errStr = NULL) : err_(err), errStr_(errStr)
    {}

    ~Error() throw() {}

    /*! \brief Get error string associated with exception
     *
     * \return A memory pointer to the error message string.
     */
    virtual const char * what() const throw ()
    {
        if (errStr_ == NULL) {
            return "empty";
        }
        else {
            return errStr_;
        }
    }

    /*! \brief Get error code associated with exception
     *
     *  \return The error code.
     */
    cl_int err(void) const { return err_; }
};

#define __ERR_STR(x) #x
#else
#define __ERR_STR(x) NULL
#endif // __CL_ENABLE_EXCEPTIONS


namespace detail
{
#if defined(__CL_ENABLE_EXCEPTIONS)
static inline cl_int errHandler (
    cl_int err,
    const char * errStr = NULL)
{
    if (err != CL_SUCCESS) {
        throw Error(err, errStr);
    }
    return err;
}
#else
static inline cl_int errHandler (cl_int err, const char * errStr = NULL)
{
    return err;
}
#endif // __CL_ENABLE_EXCEPTIONS
}



//! \cond DOXYGEN_DETAIL
#if !defined(__CL_USER_OVERRIDE_ERROR_STRINGS)
#define __GET_DEVICE_INFO_ERR               __ERR_STR(clGetDeviceInfo)
#define __GET_PLATFORM_INFO_ERR             __ERR_STR(clGetPlatformInfo)
#define __GET_DEVICE_IDS_ERR                __ERR_STR(clGetDeviceIDs)
#define __GET_PLATFORM_IDS_ERR              __ERR_STR(clGetPlatformIDs)
#define __GET_CONTEXT_INFO_ERR              __ERR_STR(clGetContextInfo)
#define __GET_EVENT_INFO_ERR                __ERR_STR(clGetEventInfo)
#define __GET_EVENT_PROFILE_INFO_ERR        __ERR_STR(clGetEventProfileInfo)
#define __GET_MEM_OBJECT_INFO_ERR           __ERR_STR(clGetMemObjectInfo)
#define __GET_IMAGE_INFO_ERR                __ERR_STR(clGetImageInfo)
#define __GET_SAMPLER_INFO_ERR              __ERR_STR(clGetSamplerInfo)
#define __GET_KERNEL_INFO_ERR               __ERR_STR(clGetKernelInfo)
#if defined(CL_VERSION_1_2)
#define __GET_KERNEL_ARG_INFO_ERR               __ERR_STR(clGetKernelArgInfo)
#endif // #if defined(CL_VERSION_1_2)
#define __GET_KERNEL_WORK_GROUP_INFO_ERR    __ERR_STR(clGetKernelWorkGroupInfo)
#define __GET_PROGRAM_INFO_ERR              __ERR_STR(clGetProgramInfo)
#define __GET_PROGRAM_BUILD_INFO_ERR        __ERR_STR(clGetProgramBuildInfo)
#define __GET_COMMAND_QUEUE_INFO_ERR        __ERR_STR(clGetCommandQueueInfo)

#define __CREATE_CONTEXT_ERR                __ERR_STR(clCreateContext)
#define __CREATE_CONTEXT_FROM_TYPE_ERR      __ERR_STR(clCreateContextFromType)
#define __GET_SUPPORTED_IMAGE_FORMATS_ERR   __ERR_STR(clGetSupportedImageFormats)

#define __CREATE_BUFFER_ERR                 __ERR_STR(clCreateBuffer)
#define __COPY_ERR                          __ERR_STR(cl::copy)
#define __CREATE_SUBBUFFER_ERR              __ERR_STR(clCreateSubBuffer)
#define __CREATE_GL_BUFFER_ERR              __ERR_STR(clCreateFromGLBuffer)
#define __CREATE_GL_RENDER_BUFFER_ERR       __ERR_STR(clCreateFromGLBuffer)
#define __GET_GL_OBJECT_INFO_ERR            __ERR_STR(clGetGLObjectInfo)
#if defined(CL_VERSION_1_2)
#define __CREATE_IMAGE_ERR                  __ERR_STR(clCreateImage)
#define __CREATE_GL_TEXTURE_ERR             __ERR_STR(clCreateFromGLTexture)
#define __IMAGE_DIMENSION_ERR               __ERR_STR(Incorrect image dimensions)
#endif // #if defined(CL_VERSION_1_2)
#define __CREATE_SAMPLER_ERR                __ERR_STR(clCreateSampler)
#define __SET_MEM_OBJECT_DESTRUCTOR_CALLBACK_ERR __ERR_STR(clSetMemObjectDestructorCallback)

#define __CREATE_USER_EVENT_ERR             __ERR_STR(clCreateUserEvent)
#define __SET_USER_EVENT_STATUS_ERR         __ERR_STR(clSetUserEventStatus)
#define __SET_EVENT_CALLBACK_ERR            __ERR_STR(clSetEventCallback)
#define __WAIT_FOR_EVENTS_ERR               __ERR_STR(clWaitForEvents)

#define __CREATE_KERNEL_ERR                 __ERR_STR(clCreateKernel)
#define __SET_KERNEL_ARGS_ERR               __ERR_STR(clSetKernelArg)
#define __CREATE_PROGRAM_WITH_SOURCE_ERR    __ERR_STR(clCreateProgramWithSource)
#define __CREATE_PROGRAM_WITH_BINARY_ERR    __ERR_STR(clCreateProgramWithBinary)
#if defined(CL_VERSION_1_2)
#define __CREATE_PROGRAM_WITH_BUILT_IN_KERNELS_ERR    __ERR_STR(clCreateProgramWithBuiltInKernels)
#endif // #if defined(CL_VERSION_1_2)
#define __BUILD_PROGRAM_ERR                 __ERR_STR(clBuildProgram)
#if defined(CL_VERSION_1_2)
#define __COMPILE_PROGRAM_ERR                  __ERR_STR(clCompileProgram)
#define __SET_PRINTF_CALLBACK_ERR           __ERR_STR(clSetPrintfCallback)
#endif // #if defined(CL_VERSION_1_2)
#define __CREATE_KERNELS_IN_PROGRAM_ERR     __ERR_STR(clCreateKernelsInProgram)

#define __CREATE_COMMAND_QUEUE_ERR          __ERR_STR(clCreateCommandQueue)
#define __SET_COMMAND_QUEUE_PROPERTY_ERR    __ERR_STR(clSetCommandQueueProperty)
#define __ENQUEUE_READ_BUFFER_ERR           __ERR_STR(clEnqueueReadBuffer)
#define __ENQUEUE_READ_BUFFER_RECT_ERR      __ERR_STR(clEnqueueReadBufferRect)
#define __ENQUEUE_WRITE_BUFFER_ERR          __ERR_STR(clEnqueueWriteBuffer)
#define __ENQUEUE_WRITE_BUFFER_RECT_ERR     __ERR_STR(clEnqueueWriteBufferRect)
#define __ENQEUE_COPY_BUFFER_ERR            __ERR_STR(clEnqueueCopyBuffer)
#define __ENQEUE_COPY_BUFFER_RECT_ERR       __ERR_STR(clEnqueueCopyBufferRect)
#define __ENQUEUE_FILL_BUFFER_ERR           __ERR_STR(clEnqueueFillBuffer)
#define __ENQUEUE_READ_IMAGE_ERR            __ERR_STR(clEnqueueReadImage)
#define __ENQUEUE_WRITE_IMAGE_ERR           __ERR_STR(clEnqueueWriteImage)
#define __ENQUEUE_COPY_IMAGE_ERR            __ERR_STR(clEnqueueCopyImage)
#define __ENQUEUE_FILL_IMAGE_ERR           __ERR_STR(clEnqueueFillImage)
#define __ENQUEUE_COPY_IMAGE_TO_BUFFER_ERR  __ERR_STR(clEnqueueCopyImageToBuffer)
#define __ENQUEUE_COPY_BUFFER_TO_IMAGE_ERR  __ERR_STR(clEnqueueCopyBufferToImage)
#define __ENQUEUE_MAP_BUFFER_ERR            __ERR_STR(clEnqueueMapBuffer)
#define __ENQUEUE_MAP_IMAGE_ERR             __ERR_STR(clEnqueueMapImage)
#define __ENQUEUE_UNMAP_MEM_OBJECT_ERR      __ERR_STR(clEnqueueUnMapMemObject)
#define __ENQUEUE_NDRANGE_KERNEL_ERR        __ERR_STR(clEnqueueNDRangeKernel)
#define __ENQUEUE_TASK_ERR                  __ERR_STR(clEnqueueTask)
#define __ENQUEUE_NATIVE_KERNEL             __ERR_STR(clEnqueueNativeKernel)
#if defined(CL_VERSION_1_2)
#define __ENQUEUE_MIGRATE_MEM_OBJECTS_ERR   __ERR_STR(clEnqueueMigrateMemObjects)
#endif // #if defined(CL_VERSION_1_2)

#define __ENQUEUE_ACQUIRE_GL_ERR            __ERR_STR(clEnqueueAcquireGLObjects)
#define __ENQUEUE_RELEASE_GL_ERR            __ERR_STR(clEnqueueReleaseGLObjects)


#define __RETAIN_ERR                        __ERR_STR(Retain Object)
#define __RELEASE_ERR                       __ERR_STR(Release Object)
#define __FLUSH_ERR                         __ERR_STR(clFlush)
#define __FINISH_ERR                        __ERR_STR(clFinish)
#define __VECTOR_CAPACITY_ERR               __ERR_STR(Vector capacity error)

/**
 * CL 1.1 version that uses device fission.
 */
#if defined(CL_VERSION_1_2)
#define __CREATE_SUB_DEVICES                __ERR_STR(clCreateSubDevices)
#else
#define __CREATE_SUB_DEVICES                __ERR_STR(clCreateSubDevicesEXT)
#endif // #if defined(CL_VERSION_1_2)

/**
 * Deprecated APIs for 1.2
 */
#if defined(CL_VERSION_1_1)
#define __ENQUEUE_MARKER_ERR                __ERR_STR(clEnqueueMarker)
#define __ENQUEUE_WAIT_FOR_EVENTS_ERR       __ERR_STR(clEnqueueWaitForEvents)
#define __ENQUEUE_BARRIER_ERR               __ERR_STR(clEnqueueBarrier)
#define __UNLOAD_COMPILER_ERR               __ERR_STR(clUnloadCompiler)
#define __CREATE_GL_TEXTURE_2D_ERR          __ERR_STR(clCreateFromGLTexture2D)
#define __CREATE_GL_TEXTURE_3D_ERR          __ERR_STR(clCreateFromGLTexture3D)
#define __CREATE_IMAGE2D_ERR                __ERR_STR(clCreateImage2D)
#define __CREATE_IMAGE3D_ERR                __ERR_STR(clCreateImage3D)
#endif // #if defined(CL_VERSION_1_1)

#endif // __CL_USER_OVERRIDE_ERROR_STRINGS
//! \endcond

/*! \class string
 * \brief Simple string class, that provides a limited subset of std::string
 * functionality but avoids many of the issues that come with that class.
 */
class string
{
private:
    ::size_t size_;
    char * str_;
public:
    string(void) : size_(0), str_(NULL)
    {
    }

    string(const char * str, ::size_t size) :
        size_(size),
        str_(NULL)
    {
        str_ = new char[size_+1];
        if (str_ != NULL) {
            memcpy(str_, str, size_  * sizeof(char));
            str_[size_] = '\0';
        }
        else {
            size_ = 0;
        }
    }

    string(const char * str) :
        str_(NULL)
    {
        size_= ::strlen(str);
        str_ = new char[size_ + 1];
        if (str_ != NULL) {
            memcpy(str_, str, (size_ + 1) * sizeof(char));
        }
        else {
            size_ = 0;
        }
    }

    string& operator=(const string& rhs)
    {
        if (this == &rhs) {
            return *this;
        }

        if (rhs.size_ == 0 || rhs.str_ == NULL) {
            size_ = 0;
            str_  = NULL;
        } 
        else {
            size_ = rhs.size_;
            str_ = new char[size_ + 1];
            if (str_ != NULL) {
                memcpy(str_, rhs.str_, (size_ + 1) * sizeof(char));
            }
            else {
                size_ = 0;
            }
        }

        return *this;
    }

    string(const string& rhs) :
        size_(0),
        str_(NULL)
    {
        *this = rhs;
    }

    ~string()
    {
        if (str_ != NULL) {
            delete[] str_;
            str_ = NULL;
        }
    }

    ::size_t size(void) const   { return size_; }
    ::size_t length(void) const { return size(); }

    const char * c_str(void) const { return (str_) ? str_ : "";}
};

#if !defined(__USE_DEV_STRING) && !defined(__NO_STD_STRING)
typedef std::string STRING_CLASS;
#elif !defined(__USE_DEV_STRING) 
typedef cl::string STRING_CLASS;
#endif

#if !defined(__USE_DEV_VECTOR) && !defined(__NO_STD_VECTOR)
#define VECTOR_CLASS std::vector
#elif !defined(__USE_DEV_VECTOR) 
#define VECTOR_CLASS cl::vector 
#endif

#if !defined(__MAX_DEFAULT_VECTOR_SIZE)
#define __MAX_DEFAULT_VECTOR_SIZE 10
#endif

/*! \class vector
 * \brief Fixed sized vector implementation that mirrors 
 * std::vector functionality.
 */
template <typename T, unsigned int N = __MAX_DEFAULT_VECTOR_SIZE>
class vector
{
private:
    T data_[N];

    unsigned int size_;
public:
    vector() :  
        size_(static_cast<unsigned int>(0))
    {}

    ~vector() 
    {
        clear();
    }

    unsigned int size(void) const
    {
        return size_;
    }

    void clear()
    {
        while(!empty()) {
            pop_back();
        }
    }

    void push_back (const T& x)
    { 
        if (size() < N) {    
            new (&data_[size_]) T(x);
            size_++;
        } else {
            detail::errHandler(CL_MEM_OBJECT_ALLOCATION_FAILURE, __VECTOR_CAPACITY_ERR);
        }
    }

    void pop_back(void)
    {
        if (size_ != 0) {
            --size_;
            data_[size_].~T();
        } else {
            detail::errHandler(CL_MEM_OBJECT_ALLOCATION_FAILURE, __VECTOR_CAPACITY_ERR);
        }
    }
  
    vector(const vector<T, N>& vec) : 
        size_(vec.size_)
    {
        if (size_ != 0) {	
            assign(vec.begin(), vec.end());
        }
    } 

    vector(unsigned int size, const T& val = T()) :
        size_(0)
    {
        for (unsigned int i = 0; i < size; i++) {
            push_back(val);
        }
    }

    vector<T, N>& operator=(const vector<T, N>& rhs)
    {
        if (this == &rhs) {
            return *this;
        }

        if (rhs.size_ != 0) {	
            assign(rhs.begin(), rhs.end());
        } else {
            clear();
        }
    
        return *this;
    }

    bool operator==(vector<T,N> &vec)
    {
        if (size() != vec.size()) {
            return false;
        }

        for( unsigned int i = 0; i < size(); ++i ) {
            if( operator[](i) != vec[i] ) {
                return false;
            }
        }
        return true;
    }
  
    operator T* ()             { return data_; }
    operator const T* () const { return data_; }
   
    bool empty (void) const
    {
        return size_==0;
    }
  
    unsigned int max_size (void) const
    {
        return N;
    }

    unsigned int capacity () const
    {
        return N;
    }

    T& operator[](int index)
    {
        return data_[index];
    }
  
    const T& operator[](int index) const
    {
        return data_[index];
    }
  
    template<class I>
    void assign(I start, I end)
    {
        clear();   
        while(start != end) {
            push_back(*start);
            start++;
        }
    }

    /*! \class iterator
     * \brief Const iterator class for vectors
     */
    class iterator
    {
    private:
        const vector<T,N> *vec_;
        int index_;

        /**
         * Internal iterator constructor to capture reference
         * to the vector it iterates over rather than taking 
         * the vector by copy.
         */
        iterator (const vector<T,N> &vec, int index) :
            vec_(&vec)
        {            
            if( !vec.empty() ) {
                index_ = index;
            } else {
                index_ = -1;
            }
        }

    public:
        iterator(void) : 
            index_(-1),
            vec_(NULL)
        {
        }

        iterator(const iterator& rhs) :
            vec_(rhs.vec_),
            index_(rhs.index_)
        {
        }

        ~iterator(void) {}

        static iterator begin(const cl::vector<T,N> &vec)
        {
            iterator i(vec, 0);

            return i;
        }

        static iterator end(const cl::vector<T,N> &vec)
        {
            iterator i(vec, vec.size());

            return i;
        }
    
        bool operator==(iterator i)
        {
            return ((vec_ == i.vec_) && 
                    (index_ == i.index_));
        }

        bool operator!=(iterator i)
        {
            return (!(*this==i));
        }

        iterator& operator++()
        {
            ++index_;
            return *this;
        }

        iterator operator++(int)
        {
            iterator retVal(*this);
            ++index_;
            return retVal;
        }

        iterator& operator--()
        {
            --index_;
            return *this;
        }

        iterator operator--(int)
        {
            iterator retVal(*this);
            --index_;
            return retVal;
        }

        const T& operator *() const
        {
            return (*vec_)[index_];
        }
    };

    iterator begin(void)
    {
        return iterator::begin(*this);
    }

    iterator begin(void) const
    {
        return iterator::begin(*this);
    }

    iterator end(void)
    {
        return iterator::end(*this);
    }

    iterator end(void) const
    {
        return iterator::end(*this);
    }

    T& front(void)
    {
        return data_[0];
    }

    T& back(void)
    {
        return data_[size_];
    }

    const T& front(void) const
    {
        return data_[0];
    }

    const T& back(void) const
    {
        return data_[size_];
    }
};  

namespace detail {
#define __DEFAULT_NOT_INITIALIZED 1 
#define __DEFAULT_BEING_INITIALIZED 2
#define __DEFAULT_INITIALIZED 4

    /*
     * Compare and exchange primitives are needed for handling of defaults
    */
    inline int compare_exchange(volatile int * dest, int exchange, int comparand)
    {
#ifdef _WIN32
        return (int)(InterlockedCompareExchange(
           (volatile long*)dest, 
           (long)exchange, 
           (long)comparand));
#elif defined(__APPLE__) || defined(__MACOSX)
		return OSAtomicOr32Orig((uint32_t)exchange, (volatile uint32_t*)dest);
#else // !_WIN32 || defined(__APPLE__) || defined(__MACOSX)
        return (__sync_val_compare_and_swap(
            dest, 
            comparand, 
            exchange));
#endif // !_WIN32
    }
}; // namespace detail

    
/*!
 * \brief size_t class used to interface between C++ and
 * OpenCL C calls that require arrays of size_t values, who's
 * size is known statically.
 */
template <int N>
struct size_t : public cl::vector< ::size_t, N> { };

namespace detail {

// GetInfo help struct
template <typename Functor, typename T>
struct GetInfoHelper
{
    static cl_int
    get(Functor f, cl_uint name, T* param)
    {
        return f(name, sizeof(T), param, NULL);
    }
};

// Specialized GetInfoHelper for VECTOR_CLASS params
template <typename Func, typename T>
struct GetInfoHelper<Func, VECTOR_CLASS<T> >
{
    static cl_int get(Func f, cl_uint name, VECTOR_CLASS<T>* param)
    {
        ::size_t required;
        cl_int err = f(name, 0, NULL, &required);
        if (err != CL_SUCCESS) {
            return err;
        }

        T* value = (T*) alloca(required);
        err = f(name, required, value, NULL);
        if (err != CL_SUCCESS) {
            return err;
        }

        param->assign(&value[0], &value[required/sizeof(T)]);
        return CL_SUCCESS;
    }
};

// Specialized for getInfo<CL_PROGRAM_BINARIES>
template <typename Func>
struct GetInfoHelper<Func, VECTOR_CLASS<char *> >
{
    static cl_int
    get(Func f, cl_uint name, VECTOR_CLASS<char *>* param)
    {
      cl_uint err = f(name, param->size() * sizeof(char *), &(*param)[0], NULL);

      if (err != CL_SUCCESS) {
        return err;
      }

      return CL_SUCCESS;
    }
};

// Specialized GetInfoHelper for STRING_CLASS params
template <typename Func>
struct GetInfoHelper<Func, STRING_CLASS>
{
    static cl_int get(Func f, cl_uint name, STRING_CLASS* param)
    {
        ::size_t required;
        cl_int err = f(name, 0, NULL, &required);
        if (err != CL_SUCCESS) {
            return err;
        }

        char* value = (char*) alloca(required);
        err = f(name, required, value, NULL);
        if (err != CL_SUCCESS) {
            return err;
        }

        *param = value;
        return CL_SUCCESS;
    }
};

// Specialized GetInfoHelper for cl::size_t params
template <typename Func, ::size_t N>
struct GetInfoHelper<Func, size_t<N> >
{
    static cl_int get(Func f, cl_uint name, size_t<N>* param)
    {
        ::size_t required;
        cl_int err = f(name, 0, NULL, &required);
        if (err != CL_SUCCESS) {
            return err;
        }

        ::size_t* value = (::size_t*) alloca(required);
        err = f(name, required, value, NULL);
        if (err != CL_SUCCESS) {
            return err;
        }

        param->assign(&value[0], &value[required/sizeof(::size_t)]);
        return CL_SUCCESS;
    }
};

#define __GET_INFO_HELPER_WITH_RETAIN(CPP_TYPE) \
namespace detail { \
template <typename Func> \
struct GetInfoHelper<Func, CPP_TYPE> \
{ \
    static cl_int get(Func f, cl_uint name, CPP_TYPE* param) \
    { \
      cl_uint err = f(name, sizeof(CPP_TYPE), param, NULL); \
      if (err != CL_SUCCESS) { \
        return err; \
      } \
      \
      return ReferenceHandler<CPP_TYPE::cl_type>::retain((*param)()); \
    } \
}; \
} 


#define __PARAM_NAME_INFO_1_0(F) \
    F(cl_platform_info, CL_PLATFORM_PROFILE, STRING_CLASS) \
    F(cl_platform_info, CL_PLATFORM_VERSION, STRING_CLASS) \
    F(cl_platform_info, CL_PLATFORM_NAME, STRING_CLASS) \
    F(cl_platform_info, CL_PLATFORM_VENDOR, STRING_CLASS) \
    F(cl_platform_info, CL_PLATFORM_EXTENSIONS, STRING_CLASS) \
    \
    F(cl_device_info, CL_DEVICE_TYPE, cl_device_type) \
    F(cl_device_info, CL_DEVICE_VENDOR_ID, cl_uint) \
    F(cl_device_info, CL_DEVICE_MAX_COMPUTE_UNITS, cl_uint) \
    F(cl_device_info, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, cl_uint) \
    F(cl_device_info, CL_DEVICE_MAX_WORK_GROUP_SIZE, ::size_t) \
    F(cl_device_info, CL_DEVICE_MAX_WORK_ITEM_SIZES, VECTOR_CLASS< ::size_t>) \
    F(cl_device_info, CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR, cl_uint) \
    F(cl_device_info, CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT, cl_uint) \
    F(cl_device_info, CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT, cl_uint) \
    F(cl_device_info, CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG, cl_uint) \
    F(cl_device_info, CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT, cl_uint) \
    F(cl_device_info, CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE, cl_uint) \
    F(cl_device_info, CL_DEVICE_MAX_CLOCK_FREQUENCY, cl_uint) \
    F(cl_device_info, CL_DEVICE_ADDRESS_BITS, cl_uint) \
    F(cl_device_info, CL_DEVICE_MAX_READ_IMAGE_ARGS, cl_uint) \
    F(cl_device_info, CL_DEVICE_MAX_WRITE_IMAGE_ARGS, cl_uint) \
    F(cl_device_info, CL_DEVICE_MAX_MEM_ALLOC_SIZE, cl_ulong) \
    F(cl_device_info, CL_DEVICE_IMAGE2D_MAX_WIDTH, ::size_t) \
    F(cl_device_info, CL_DEVICE_IMAGE2D_MAX_HEIGHT, ::size_t) \
    F(cl_device_info, CL_DEVICE_IMAGE3D_MAX_WIDTH, ::size_t) \
    F(cl_device_info, CL_DEVICE_IMAGE3D_MAX_HEIGHT, ::size_t) \
    F(cl_device_info, CL_DEVICE_IMAGE3D_MAX_DEPTH, ::size_t) \
    F(cl_device_info, CL_DEVICE_IMAGE_SUPPORT, cl_bool) \
    F(cl_device_info, CL_DEVICE_MAX_PARAMETER_SIZE, ::size_t) \
    F(cl_device_info, CL_DEVICE_MAX_SAMPLERS, cl_uint) \
    F(cl_device_info, CL_DEVICE_MEM_BASE_ADDR_ALIGN, cl_uint) \
    F(cl_device_info, CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE, cl_uint) \
    F(cl_device_info, CL_DEVICE_SINGLE_FP_CONFIG, cl_device_fp_config) \
    F(cl_device_info, CL_DEVICE_GLOBAL_MEM_CACHE_TYPE, cl_device_mem_cache_type) \
    F(cl_device_info, CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE, cl_uint)\
    F(cl_device_info, CL_DEVICE_GLOBAL_MEM_CACHE_SIZE, cl_ulong) \
    F(cl_device_info, CL_DEVICE_GLOBAL_MEM_SIZE, cl_ulong) \
    F(cl_device_info, CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, cl_ulong) \
    F(cl_device_info, CL_DEVICE_MAX_CONSTANT_ARGS, cl_uint) \
    F(cl_device_info, CL_DEVICE_LOCAL_MEM_TYPE, cl_device_local_mem_type) \
    F(cl_device_info, CL_DEVICE_LOCAL_MEM_SIZE, cl_ulong) \
    F(cl_device_info, CL_DEVICE_ERROR_CORRECTION_SUPPORT, cl_bool) \
    F(cl_device_info, CL_DEVICE_PROFILING_TIMER_RESOLUTION, ::size_t) \
    F(cl_device_info, CL_DEVICE_ENDIAN_LITTLE, cl_bool) \
    F(cl_device_info, CL_DEVICE_AVAILABLE, cl_bool) \
    F(cl_device_info, CL_DEVICE_COMPILER_AVAILABLE, cl_bool) \
    F(cl_device_info, CL_DEVICE_EXECUTION_CAPABILITIES, cl_device_exec_capabilities) \
    F(cl_device_info, CL_DEVICE_QUEUE_PROPERTIES, cl_command_queue_properties) \
    F(cl_device_info, CL_DEVICE_PLATFORM, cl_platform_id) \
    F(cl_device_info, CL_DEVICE_NAME, STRING_CLASS) \
    F(cl_device_info, CL_DEVICE_VENDOR, STRING_CLASS) \
    F(cl_device_info, CL_DRIVER_VERSION, STRING_CLASS) \
    F(cl_device_info, CL_DEVICE_PROFILE, STRING_CLASS) \
    F(cl_device_info, CL_DEVICE_VERSION, STRING_CLASS) \
    F(cl_device_info, CL_DEVICE_EXTENSIONS, STRING_CLASS) \
    \
    F(cl_context_info, CL_CONTEXT_REFERENCE_COUNT, cl_uint) \
    F(cl_context_info, CL_CONTEXT_DEVICES, VECTOR_CLASS<Device>) \
    F(cl_context_info, CL_CONTEXT_PROPERTIES, VECTOR_CLASS<cl_context_properties>) \
    \
    F(cl_event_info, CL_EVENT_COMMAND_QUEUE, cl::CommandQueue) \
    F(cl_event_info, CL_EVENT_COMMAND_TYPE, cl_command_type) \
    F(cl_event_info, CL_EVENT_REFERENCE_COUNT, cl_uint) \
    F(cl_event_info, CL_EVENT_COMMAND_EXECUTION_STATUS, cl_uint) \
    \
    F(cl_profiling_info, CL_PROFILING_COMMAND_QUEUED, cl_ulong) \
    F(cl_profiling_info, CL_PROFILING_COMMAND_SUBMIT, cl_ulong) \
    F(cl_profiling_info, CL_PROFILING_COMMAND_START, cl_ulong) \
    F(cl_profiling_info, CL_PROFILING_COMMAND_END, cl_ulong) \
    \
    F(cl_mem_info, CL_MEM_TYPE, cl_mem_object_type) \
    F(cl_mem_info, CL_MEM_FLAGS, cl_mem_flags) \
    F(cl_mem_info, CL_MEM_SIZE, ::size_t) \
    F(cl_mem_info, CL_MEM_HOST_PTR, void*) \
    F(cl_mem_info, CL_MEM_MAP_COUNT, cl_uint) \
    F(cl_mem_info, CL_MEM_REFERENCE_COUNT, cl_uint) \
    F(cl_mem_info, CL_MEM_CONTEXT, cl::Context) \
    \
    F(cl_image_info, CL_IMAGE_FORMAT, cl_image_format) \
    F(cl_image_info, CL_IMAGE_ELEMENT_SIZE, ::size_t) \
    F(cl_image_info, CL_IMAGE_ROW_PITCH, ::size_t) \
    F(cl_image_info, CL_IMAGE_SLICE_PITCH, ::size_t) \
    F(cl_image_info, CL_IMAGE_WIDTH, ::size_t) \
    F(cl_image_info, CL_IMAGE_HEIGHT, ::size_t) \
    F(cl_image_info, CL_IMAGE_DEPTH, ::size_t) \
    \
    F(cl_sampler_info, CL_SAMPLER_REFERENCE_COUNT, cl_uint) \
    F(cl_sampler_info, CL_SAMPLER_CONTEXT, cl::Context) \
    F(cl_sampler_info, CL_SAMPLER_NORMALIZED_COORDS, cl_addressing_mode) \
    F(cl_sampler_info, CL_SAMPLER_ADDRESSING_MODE, cl_filter_mode) \
    F(cl_sampler_info, CL_SAMPLER_FILTER_MODE, cl_bool) \
    \
    F(cl_program_info, CL_PROGRAM_REFERENCE_COUNT, cl_uint) \
    F(cl_program_info, CL_PROGRAM_CONTEXT, cl::Context) \
    F(cl_program_info, CL_PROGRAM_NUM_DEVICES, cl_uint) \
    F(cl_program_info, CL_PROGRAM_DEVICES, VECTOR_CLASS<cl_device_id>) \
    F(cl_program_info, CL_PROGRAM_SOURCE, STRING_CLASS) \
    F(cl_program_info, CL_PROGRAM_BINARY_SIZES, VECTOR_CLASS< ::size_t>) \
    F(cl_program_info, CL_PROGRAM_BINARIES, VECTOR_CLASS<char *>) \
    \
    F(cl_program_build_info, CL_PROGRAM_BUILD_STATUS, cl_build_status) \
    F(cl_program_build_info, CL_PROGRAM_BUILD_OPTIONS, STRING_CLASS) \
    F(cl_program_build_info, CL_PROGRAM_BUILD_LOG, STRING_CLASS) \
    \
    F(cl_kernel_info, CL_KERNEL_FUNCTION_NAME, STRING_CLASS) \
    F(cl_kernel_info, CL_KERNEL_NUM_ARGS, cl_uint) \
    F(cl_kernel_info, CL_KERNEL_REFERENCE_COUNT, cl_uint) \
    F(cl_kernel_info, CL_KERNEL_CONTEXT, cl::Context) \
    F(cl_kernel_info, CL_KERNEL_PROGRAM, cl::Program) \
    \
    F(cl_kernel_work_group_info, CL_KERNEL_WORK_GROUP_SIZE, ::size_t) \
    F(cl_kernel_work_group_info, CL_KERNEL_COMPILE_WORK_GROUP_SIZE, cl::size_t<3>) \
    F(cl_kernel_work_group_info, CL_KERNEL_LOCAL_MEM_SIZE, cl_ulong) \
    \
    F(cl_command_queue_info, CL_QUEUE_CONTEXT, cl::Context) \
    F(cl_command_queue_info, CL_QUEUE_DEVICE, cl::Device) \
    F(cl_command_queue_info, CL_QUEUE_REFERENCE_COUNT, cl_uint) \
    F(cl_command_queue_info, CL_QUEUE_PROPERTIES, cl_command_queue_properties)

#if defined(CL_VERSION_1_1)
#define __PARAM_NAME_INFO_1_1(F) \
    F(cl_context_info, CL_CONTEXT_NUM_DEVICES, cl_uint)\
    F(cl_device_info, CL_DEVICE_PREFERRED_VECTOR_WIDTH_HALF, cl_uint) \
    F(cl_device_info, CL_DEVICE_NATIVE_VECTOR_WIDTH_CHAR, cl_uint) \
    F(cl_device_info, CL_DEVICE_NATIVE_VECTOR_WIDTH_SHORT, cl_uint) \
    F(cl_device_info, CL_DEVICE_NATIVE_VECTOR_WIDTH_INT, cl_uint) \
    F(cl_device_info, CL_DEVICE_NATIVE_VECTOR_WIDTH_LONG, cl_uint) \
    F(cl_device_info, CL_DEVICE_NATIVE_VECTOR_WIDTH_FLOAT, cl_uint) \
    F(cl_device_info, CL_DEVICE_NATIVE_VECTOR_WIDTH_DOUBLE, cl_uint) \
    F(cl_device_info, CL_DEVICE_NATIVE_VECTOR_WIDTH_HALF, cl_uint) \
    F(cl_device_info, CL_DEVICE_DOUBLE_FP_CONFIG, cl_device_fp_config) \
    F(cl_device_info, CL_DEVICE_HALF_FP_CONFIG, cl_device_fp_config) \
    F(cl_device_info, CL_DEVICE_HOST_UNIFIED_MEMORY, cl_bool) \
    F(cl_device_info, CL_DEVICE_OPENCL_C_VERSION, STRING_CLASS) \
    \
    F(cl_mem_info, CL_MEM_ASSOCIATED_MEMOBJECT, cl::Memory) \
    F(cl_mem_info, CL_MEM_OFFSET, ::size_t) \
    \
    F(cl_kernel_work_group_info, CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, ::size_t) \
    F(cl_kernel_work_group_info, CL_KERNEL_PRIVATE_MEM_SIZE, cl_ulong) \
    \
    F(cl_event_info, CL_EVENT_CONTEXT, cl::Context)
#endif // CL_VERSION_1_1

    
#if defined(CL_VERSION_1_2)
#define __PARAM_NAME_INFO_1_2(F) \
    F(cl_program_info, CL_PROGRAM_NUM_KERNELS, ::size_t) \
    F(cl_program_info, CL_PROGRAM_KERNEL_NAMES, STRING_CLASS) \
    \
    F(cl_program_build_info, CL_PROGRAM_BINARY_TYPE, cl_program_binary_type) \
    \
    F(cl_kernel_info, CL_KERNEL_ATTRIBUTES, STRING_CLASS) \
    \
    F(cl_kernel_arg_info, CL_KERNEL_ARG_ADDRESS_QUALIFIER, cl_kernel_arg_address_qualifier) \
    F(cl_kernel_arg_info, CL_KERNEL_ARG_ACCESS_QUALIFIER, cl_kernel_arg_access_qualifier) \
    F(cl_kernel_arg_info, CL_KERNEL_ARG_TYPE_NAME, STRING_CLASS) \
    F(cl_kernel_arg_info, CL_KERNEL_ARG_NAME, STRING_CLASS) \
    \
    F(cl_device_info, CL_DEVICE_PARENT_DEVICE, cl_device_id) \
    F(cl_device_info, CL_DEVICE_PARTITION_PROPERTIES, VECTOR_CLASS<cl_device_partition_property>) \
    F(cl_device_info, CL_DEVICE_PARTITION_TYPE, VECTOR_CLASS<cl_device_partition_property>)  \
    F(cl_device_info, CL_DEVICE_REFERENCE_COUNT, cl_uint) \
    F(cl_device_info, CL_DEVICE_PREFERRED_INTEROP_USER_SYNC, ::size_t) \
    F(cl_device_info, CL_DEVICE_PARTITION_AFFINITY_DOMAIN, cl_device_affinity_domain) \
    F(cl_device_info, CL_DEVICE_BUILT_IN_KERNELS, STRING_CLASS)
#endif // #if defined(CL_VERSION_1_2)

#if defined(USE_CL_DEVICE_FISSION)
#define __PARAM_NAME_DEVICE_FISSION(F) \
    F(cl_device_info, CL_DEVICE_PARENT_DEVICE_EXT, cl_device_id) \
    F(cl_device_info, CL_DEVICE_PARTITION_TYPES_EXT, VECTOR_CLASS<cl_device_partition_property_ext>) \
    F(cl_device_info, CL_DEVICE_AFFINITY_DOMAINS_EXT, VECTOR_CLASS<cl_device_partition_property_ext>) \
    F(cl_device_info, CL_DEVICE_REFERENCE_COUNT_EXT , cl_uint) \
    F(cl_device_info, CL_DEVICE_PARTITION_STYLE_EXT, VECTOR_CLASS<cl_device_partition_property_ext>)
#endif // USE_CL_DEVICE_FISSION

template <typename enum_type, cl_int Name>
struct param_traits {};

#define __CL_DECLARE_PARAM_TRAITS(token, param_name, T) \
struct token;                                        \
template<>                                           \
struct param_traits<detail:: token,param_name>       \
{                                                    \
    enum { value = param_name };                     \
    typedef T param_type;                            \
};

__PARAM_NAME_INFO_1_0(__CL_DECLARE_PARAM_TRAITS)
#if defined(CL_VERSION_1_1)
__PARAM_NAME_INFO_1_1(__CL_DECLARE_PARAM_TRAITS)
#endif // CL_VERSION_1_1
#if defined(CL_VERSION_1_2)
__PARAM_NAME_INFO_1_2(__CL_DECLARE_PARAM_TRAITS)
#endif // CL_VERSION_1_1

#if defined(USE_CL_DEVICE_FISSION)
__PARAM_NAME_DEVICE_FISSION(__CL_DECLARE_PARAM_TRAITS);
#endif // USE_CL_DEVICE_FISSION

#ifdef CL_PLATFORM_ICD_SUFFIX_KHR
__CL_DECLARE_PARAM_TRAITS(cl_platform_info, CL_PLATFORM_ICD_SUFFIX_KHR, STRING_CLASS)
#endif

#ifdef CL_DEVICE_PROFILING_TIMER_OFFSET_AMD
__CL_DECLARE_PARAM_TRAITS(cl_device_info, CL_DEVICE_PROFILING_TIMER_OFFSET_AMD, cl_ulong)
#endif

#ifdef CL_DEVICE_GLOBAL_FREE_MEMORY_AMD
__CL_DECLARE_PARAM_TRAITS(cl_device_info, CL_DEVICE_GLOBAL_FREE_MEMORY_AMD, VECTOR_CLASS< ::size_t>)
#endif
#ifdef CL_DEVICE_SIMD_PER_COMPUTE_UNIT_AMD
__CL_DECLARE_PARAM_TRAITS(cl_device_info, CL_DEVICE_SIMD_PER_COMPUTE_UNIT_AMD, cl_uint)
#endif
#ifdef CL_DEVICE_SIMD_WIDTH_AMD
__CL_DECLARE_PARAM_TRAITS(cl_device_info, CL_DEVICE_SIMD_WIDTH_AMD, cl_uint)
#endif
#ifdef CL_DEVICE_SIMD_INSTRUCTION_WIDTH_AMD
__CL_DECLARE_PARAM_TRAITS(cl_device_info, CL_DEVICE_SIMD_INSTRUCTION_WIDTH_AMD, cl_uint)
#endif
#ifdef CL_DEVICE_WAVEFRONT_WIDTH_AMD
__CL_DECLARE_PARAM_TRAITS(cl_device_info, CL_DEVICE_WAVEFRONT_WIDTH_AMD, cl_uint)
#endif
#ifdef CL_DEVICE_GLOBAL_MEM_CHANNELS_AMD
__CL_DECLARE_PARAM_TRAITS(cl_device_info, CL_DEVICE_GLOBAL_MEM_CHANNELS_AMD, cl_uint)
#endif
#ifdef CL_DEVICE_GLOBAL_MEM_CHANNEL_BANKS_AMD
__CL_DECLARE_PARAM_TRAITS(cl_device_info, CL_DEVICE_GLOBAL_MEM_CHANNEL_BANKS_AMD, cl_uint)
#endif
#ifdef CL_DEVICE_GLOBAL_MEM_CHANNEL_BANK_WIDTH_AMD
__CL_DECLARE_PARAM_TRAITS(cl_device_info, CL_DEVICE_GLOBAL_MEM_CHANNEL_BANK_WIDTH_AMD, cl_uint)
#endif
#ifdef CL_DEVICE_LOCAL_MEM_SIZE_PER_COMPUTE_UNIT_AMD
__CL_DECLARE_PARAM_TRAITS(cl_device_info, CL_DEVICE_LOCAL_MEM_SIZE_PER_COMPUTE_UNIT_AMD, cl_uint)
#endif
#ifdef CL_DEVICE_LOCAL_MEM_BANKS_AMD
__CL_DECLARE_PARAM_TRAITS(cl_device_info, CL_DEVICE_LOCAL_MEM_BANKS_AMD, cl_uint)
#endif

#ifdef CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV
__CL_DECLARE_PARAM_TRAITS(cl_device_info, CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV, cl_uint)
#endif
#ifdef CL_DEVICE_COMPUTE_CAPABILITY_MINOR_NV
__CL_DECLARE_PARAM_TRAITS(cl_device_info, CL_DEVICE_COMPUTE_CAPABILITY_MINOR_NV, cl_uint)
#endif
#ifdef CL_DEVICE_REGISTERS_PER_BLOCK_NV
__CL_DECLARE_PARAM_TRAITS(cl_device_info, CL_DEVICE_REGISTERS_PER_BLOCK_NV, cl_uint)
#endif
#ifdef CL_DEVICE_WARP_SIZE_NV
__CL_DECLARE_PARAM_TRAITS(cl_device_info, CL_DEVICE_WARP_SIZE_NV, cl_uint)
#endif
#ifdef CL_DEVICE_GPU_OVERLAP_NV
__CL_DECLARE_PARAM_TRAITS(cl_device_info, CL_DEVICE_GPU_OVERLAP_NV, cl_bool)
#endif
#ifdef CL_DEVICE_KERNEL_EXEC_TIMEOUT_NV
__CL_DECLARE_PARAM_TRAITS(cl_device_info, CL_DEVICE_KERNEL_EXEC_TIMEOUT_NV, cl_bool)
#endif
#ifdef CL_DEVICE_INTEGRATED_MEMORY_NV
__CL_DECLARE_PARAM_TRAITS(cl_device_info, CL_DEVICE_INTEGRATED_MEMORY_NV, cl_bool)
#endif

// Convenience functions

template <typename Func, typename T>
inline cl_int
getInfo(Func f, cl_uint name, T* param)
{
    return GetInfoHelper<Func, T>::get(f, name, param);
}

template <typename Func, typename Arg0>
struct GetInfoFunctor0
{
    Func f_; const Arg0& arg0_;
    cl_int operator ()(
        cl_uint param, ::size_t size, void* value, ::size_t* size_ret)
    { return f_(arg0_, param, size, value, size_ret); }
};

template <typename Func, typename Arg0, typename Arg1>
struct GetInfoFunctor1
{
    Func f_; const Arg0& arg0_; const Arg1& arg1_;
    cl_int operator ()(
        cl_uint param, ::size_t size, void* value, ::size_t* size_ret)
    { return f_(arg0_, arg1_, param, size, value, size_ret); }
};

template <typename Func, typename Arg0, typename T>
inline cl_int
getInfo(Func f, const Arg0& arg0, cl_uint name, T* param)
{
    GetInfoFunctor0<Func, Arg0> f0 = { f, arg0 };
    return GetInfoHelper<GetInfoFunctor0<Func, Arg0>, T>
        ::get(f0, name, param);
}

template <typename Func, typename Arg0, typename Arg1, typename T>
inline cl_int
getInfo(Func f, const Arg0& arg0, const Arg1& arg1, cl_uint name, T* param)
{
    GetInfoFunctor1<Func, Arg0, Arg1> f0 = { f, arg0, arg1 };
    return GetInfoHelper<GetInfoFunctor1<Func, Arg0, Arg1>, T>
        ::get(f0, name, param);
}

template<typename T>
struct ReferenceHandler
{ };

#if defined(CL_VERSION_1_2)
/**
 * OpenCL 1.2 devices do have retain/release.
 */
template <>
struct ReferenceHandler<cl_device_id>
{
    /**
     * Retain the device.
     * \param device A valid device created using createSubDevices
     * \return 
     *   CL_SUCCESS if the function executed successfully.
     *   CL_INVALID_DEVICE if device was not a valid subdevice
     *   CL_OUT_OF_RESOURCES
     *   CL_OUT_OF_HOST_MEMORY
     */
    static cl_int retain(cl_device_id device)
    { return ::clRetainDevice(device); }
    /**
     * Retain the device.
     * \param device A valid device created using createSubDevices
     * \return 
     *   CL_SUCCESS if the function executed successfully.
     *   CL_INVALID_DEVICE if device was not a valid subdevice
     *   CL_OUT_OF_RESOURCES
     *   CL_OUT_OF_HOST_MEMORY
     */
    static cl_int release(cl_device_id device)
    { return ::clReleaseDevice(device); }
};
#else // #if defined(CL_VERSION_1_2)
/**
 * OpenCL 1.1 devices do not have retain/release.
 */
template <>
struct ReferenceHandler<cl_device_id>
{
    // cl_device_id does not have retain().
    static cl_int retain(cl_device_id)
    { return CL_SUCCESS; }
    // cl_device_id does not have release().
    static cl_int release(cl_device_id)
    { return CL_SUCCESS; }
};
#endif // #if defined(CL_VERSION_1_2)

template <>
struct ReferenceHandler<cl_platform_id>
{
    // cl_platform_id does not have retain().
    static cl_int retain(cl_platform_id)
    { return CL_SUCCESS; }
    // cl_platform_id does not have release().
    static cl_int release(cl_platform_id)
    { return CL_SUCCESS; }
};

template <>
struct ReferenceHandler<cl_context>
{
    static cl_int retain(cl_context context)
    { return ::clRetainContext(context); }
    static cl_int release(cl_context context)
    { return ::clReleaseContext(context); }
};

template <>
struct ReferenceHandler<cl_command_queue>
{
    static cl_int retain(cl_command_queue queue)
    { return ::clRetainCommandQueue(queue); }
    static cl_int release(cl_command_queue queue)
    { return ::clReleaseCommandQueue(queue); }
};

template <>
struct ReferenceHandler<cl_mem>
{
    static cl_int retain(cl_mem memory)
    { return ::clRetainMemObject(memory); }
    static cl_int release(cl_mem memory)
    { return ::clReleaseMemObject(memory); }
};

template <>
struct ReferenceHandler<cl_sampler>
{
    static cl_int retain(cl_sampler sampler)
    { return ::clRetainSampler(sampler); }
    static cl_int release(cl_sampler sampler)
    { return ::clReleaseSampler(sampler); }
};

template <>
struct ReferenceHandler<cl_program>
{
    static cl_int retain(cl_program program)
    { return ::clRetainProgram(program); }
    static cl_int release(cl_program program)
    { return ::clReleaseProgram(program); }
};

template <>
struct ReferenceHandler<cl_kernel>
{
    static cl_int retain(cl_kernel kernel)
    { return ::clRetainKernel(kernel); }
    static cl_int release(cl_kernel kernel)
    { return ::clReleaseKernel(kernel); }
};

template <>
struct ReferenceHandler<cl_event>
{
    static cl_int retain(cl_event event)
    { return ::clRetainEvent(event); }
    static cl_int release(cl_event event)
    { return ::clReleaseEvent(event); }
};

template <typename T>
class Wrapper
{
public:
    typedef T cl_type;

protected:
    cl_type object_;

public:
    Wrapper() : object_(NULL) { }

    Wrapper(const cl_type &obj) : object_(obj) { }

    ~Wrapper()
    {
        if (object_ != NULL) { detail::errHandler(release(), __RELEASE_ERR); }
    }

    Wrapper(const Wrapper<cl_type>& rhs)
    {
        object_ = rhs.object_;
        if (object_ != NULL) { detail::errHandler(retain(), __RETAIN_ERR); }
    }

    Wrapper<cl_type>& operator = (const Wrapper<cl_type>& rhs)
    {
        if (object_ != NULL) { detail::errHandler(release(), __RELEASE_ERR); }
        object_ = rhs.object_;
        if (object_ != NULL) { detail::errHandler(retain(), __RETAIN_ERR); }
        return *this;
    }

    Wrapper<cl_type>& operator = (const cl_type &rhs)
    {
        if (object_ != NULL) { detail::errHandler(release(), __RELEASE_ERR); }
        object_ = rhs;
        return *this;
    }

    cl_type operator ()() const { return object_; }

    cl_type& operator ()() { return object_; }

protected:

    cl_int retain() const
    {
        return ReferenceHandler<cl_type>::retain(object_);
    }

    cl_int release() const
    {
        return ReferenceHandler<cl_type>::release(object_);
    }
};

} // namespace detail
//! \endcond

/*! \stuct ImageFormat
 * \brief ImageFormat interface fro cl_image_format.
 */
struct ImageFormat : public cl_image_format
{
    ImageFormat(){}

    ImageFormat(cl_channel_order order, cl_channel_type type)
    {
        image_channel_order = order;
        image_channel_data_type = type;
    }

    ImageFormat& operator = (const ImageFormat& rhs)
    {
        if (this != &rhs) {
            this->image_channel_data_type = rhs.image_channel_data_type;
            this->image_channel_order     = rhs.image_channel_order;
        }
        return *this;
    }
};

/*! \class Device
 * \brief Device interface for cl_device_id.
 */
class Device : public detail::Wrapper<cl_device_id>
{
public:
    Device() : detail::Wrapper<cl_type>() { }

    Device(const Device& device) : detail::Wrapper<cl_type>(device) { }

    Device(const cl_device_id &device) : detail::Wrapper<cl_type>(device) { }

    static Device getDefault(cl_int * err = NULL);

    Device& operator = (const Device& rhs)
    {
        if (this != &rhs) {
            detail::Wrapper<cl_type>::operator=(rhs);
        }
        return *this;
    }

    Device& operator = (const cl_device_id& rhs)
    {
        detail::Wrapper<cl_type>::operator=(rhs);
        return *this;
    }

    template <typename T>
    cl_int getInfo(cl_device_info name, T* param) const
    {
        return detail::errHandler(
            detail::getInfo(&::clGetDeviceInfo, object_, name, param),
            __GET_DEVICE_INFO_ERR);
    }

    template <cl_int name> typename
    detail::param_traits<detail::cl_device_info, name>::param_type
    getInfo(cl_int* err = NULL) const
    {
        typename detail::param_traits<
            detail::cl_device_info, name>::param_type param;
        cl_int result = getInfo(name, &param);
        if (err != NULL) {
            *err = result;
        }
        return param;
    }

    /**
     * CL 1.2 version
     */
#if defined(CL_VERSION_1_2)
    cl_int createSubDevices(
        const cl_device_partition_property * properties,
        VECTOR_CLASS<Device>* devices)
    {
        cl_uint n = 0;
        cl_int err = clCreateSubDevices(object_, properties, 0, NULL, &n);
        if (err != CL_SUCCESS) {
            return detail::errHandler(err, __CREATE_SUB_DEVICES);
        }

        cl_device_id* ids = (cl_device_id*) alloca(n * sizeof(cl_device_id));
        err = clCreateSubDevices(object_, properties, n, ids, NULL);
        if (err != CL_SUCCESS) {
            return detail::errHandler(err, __CREATE_SUB_DEVICES);
        }

        devices->assign(&ids[0], &ids[n]);
        return CL_SUCCESS;
    }
#endif // #if defined(CL_VERSION_1_2)

/**
 * CL 1.1 version that uses device fission.
 */
#if defined(CL_VERSION_1_1)
#if defined(USE_CL_DEVICE_FISSION)
    cl_int createSubDevices(
        const cl_device_partition_property_ext * properties,
        VECTOR_CLASS<Device>* devices)
    {
        typedef CL_API_ENTRY cl_int 
            ( CL_API_CALL * PFN_clCreateSubDevicesEXT)(
                cl_device_id /*in_device*/,
                const cl_device_partition_property_ext * /* properties */,
                cl_uint /*num_entries*/,
                cl_device_id * /*out_devices*/,
                cl_uint * /*num_devices*/ ) CL_EXT_SUFFIX__VERSION_1_1;

        static PFN_clCreateSubDevicesEXT pfn_clCreateSubDevicesEXT = NULL;
        __INIT_CL_EXT_FCN_PTR(clCreateSubDevicesEXT);

        cl_uint n = 0;
        cl_int err = pfn_clCreateSubDevicesEXT(object_, properties, 0, NULL, &n);
        if (err != CL_SUCCESS) {
            return detail::errHandler(err, __CREATE_SUB_DEVICES);
        }

        cl_device_id* ids = (cl_device_id*) alloca(n * sizeof(cl_device_id));
        err = pfn_clCreateSubDevicesEXT(object_, properties, n, ids, NULL);
        if (err != CL_SUCCESS) {
            return detail::errHandler(err, __CREATE_SUB_DEVICES);
        }

        devices->assign(&ids[0], &ids[n]);
        return CL_SUCCESS;
    }
#endif // #if defined(USE_CL_DEVICE_FISSION)
#endif // #if defined(CL_VERSION_1_1)
};

/*! \class Platform
 *  \brief Platform interface.
 */
class Platform : public detail::Wrapper<cl_platform_id>
{
public:
    static const Platform null();

    Platform() : detail::Wrapper<cl_type>()  { }

    Platform(const Platform& platform) : detail::Wrapper<cl_type>(platform) { }

    Platform(const cl_platform_id &platform) : detail::Wrapper<cl_type>(platform) { }

    Platform& operator = (const Platform& rhs)
    {
        if (this != &rhs) {
            detail::Wrapper<cl_type>::operator=(rhs);
        }
        return *this;
    }

    Platform& operator = (const cl_platform_id& rhs)
    {
        detail::Wrapper<cl_type>::operator=(rhs);
        return *this;
    }

    cl_int getInfo(cl_platform_info name, STRING_CLASS* param) const
    {
        return detail::errHandler(
            detail::getInfo(&::clGetPlatformInfo, object_, name, param),
            __GET_PLATFORM_INFO_ERR);
    }

    template <cl_int name> typename
    detail::param_traits<detail::cl_platform_info, name>::param_type
    getInfo(cl_int* err = NULL) const
    {
        typename detail::param_traits<
            detail::cl_platform_info, name>::param_type param;
        cl_int result = getInfo(name, &param);
        if (err != NULL) {
            *err = result;
        }
        return param;
    }

    cl_int getDevices(
        cl_device_type type,
        VECTOR_CLASS<Device>* devices) const
    {
        cl_uint n = 0;
        cl_int err = ::clGetDeviceIDs(object_, type, 0, NULL, &n);
        if (err != CL_SUCCESS) {
            return detail::errHandler(err, __GET_DEVICE_IDS_ERR);
        }

        cl_device_id* ids = (cl_device_id*) alloca(n * sizeof(cl_device_id));
        err = ::clGetDeviceIDs(object_, type, n, ids, NULL);
        if (err != CL_SUCCESS) {
            return detail::errHandler(err, __GET_DEVICE_IDS_ERR);
        }

        devices->assign(&ids[0], &ids[n]);
        return CL_SUCCESS;
    }

#if defined(USE_DX_INTEROP)
   /*! \brief Get the list of available D3D10 devices.
     *
     *  \param d3d_device_source.
     *
     *  \param d3d_object.
     *
     *  \param d3d_device_set.
     *
     *  \param devices returns a vector of OpenCL D3D10 devices found. The cl::Device
     *  values returned in devices can be used to identify a specific OpenCL
     *  device. If \a devices argument is NULL, this argument is ignored.
     *
     *  \return One of the following values:
     *    - CL_SUCCESS if the function is executed successfully.
     *
     *  The application can query specific capabilities of the OpenCL device(s)
     *  returned by cl::getDevices. This can be used by the application to
     *  determine which device(s) to use.
     *
     * \note In the case that exceptions are enabled and a return value
     * other than CL_SUCCESS is generated, then cl::Error exception is
     * generated.
     */
    cl_int getDevices(
        cl_d3d10_device_source_khr d3d_device_source,
        void *                     d3d_object,
        cl_d3d10_device_set_khr    d3d_device_set,
        VECTOR_CLASS<Device>* devices) const
    {
        typedef CL_API_ENTRY cl_int (CL_API_CALL *PFN_clGetDeviceIDsFromD3D10KHR)(
            cl_platform_id platform, 
            cl_d3d10_device_source_khr d3d_device_source, 
            void * d3d_object,
            cl_d3d10_device_set_khr d3d_device_set,
            cl_uint num_entries,
            cl_device_id * devices,
            cl_uint* num_devices);

        static PFN_clGetDeviceIDsFromD3D10KHR pfn_clGetDeviceIDsFromD3D10KHR = NULL;
        __INIT_CL_EXT_FCN_PTR_PLATFORM(object_, clGetDeviceIDsFromD3D10KHR);

        cl_uint n = 0;
        cl_int err = pfn_clGetDeviceIDsFromD3D10KHR(
            object_, 
            d3d_device_source, 
            d3d_object,
            d3d_device_set, 
            0, 
            NULL, 
            &n);
        if (err != CL_SUCCESS) {
            return detail::errHandler(err, __GET_DEVICE_IDS_ERR);
        }

        cl_device_id* ids = (cl_device_id*) alloca(n * sizeof(cl_device_id));
        err = pfn_clGetDeviceIDsFromD3D10KHR(
            object_, 
            d3d_device_source, 
            d3d_object,
            d3d_device_set,
            n, 
            ids, 
            NULL);
        if (err != CL_SUCCESS) {
            return detail::errHandler(err, __GET_DEVICE_IDS_ERR);
        }

        devices->assign(&ids[0], &ids[n]);
        return CL_SUCCESS;
    }
#endif

    static cl_int get(
        VECTOR_CLASS<Platform>* platforms)
    {
        cl_uint n = 0;
        cl_int err = ::clGetPlatformIDs(0, NULL, &n);
        if (err != CL_SUCCESS) {
            return detail::errHandler(err, __GET_PLATFORM_IDS_ERR);
        }

        cl_platform_id* ids = (cl_platform_id*) alloca(
            n * sizeof(cl_platform_id));
        err = ::clGetPlatformIDs(n, ids, NULL);
        if (err != CL_SUCCESS) {
            return detail::errHandler(err, __GET_PLATFORM_IDS_ERR);
        }

        platforms->assign(&ids[0], &ids[n]);
        return CL_SUCCESS;
    }

    static cl_int get(
        Platform * platform)
    {
        cl_uint n = 0;
        cl_int err = ::clGetPlatformIDs(0, NULL, &n);
        if (err != CL_SUCCESS) {
            return detail::errHandler(err, __GET_PLATFORM_IDS_ERR);
        }

        cl_platform_id* ids = (cl_platform_id*) alloca(
            n * sizeof(cl_platform_id));
        err = ::clGetPlatformIDs(n, ids, NULL);
        if (err != CL_SUCCESS) {
            return detail::errHandler(err, __GET_PLATFORM_IDS_ERR);
        }

        *platform = ids[0];
        return CL_SUCCESS;
    }

    static Platform get(
        cl_int * errResult = NULL)
    {
        Platform platform;
        cl_uint n = 0;
        cl_int err = ::clGetPlatformIDs(0, NULL, &n);
        if (err != CL_SUCCESS) {
            detail::errHandler(err, __GET_PLATFORM_IDS_ERR);
            if (errResult != NULL) {
                *errResult = err;
            }
        }

        cl_platform_id* ids = (cl_platform_id*) alloca(
            n * sizeof(cl_platform_id));
        err = ::clGetPlatformIDs(n, ids, NULL);

        if (err != CL_SUCCESS) {
            detail::errHandler(err, __GET_PLATFORM_IDS_ERR);
        }

        if (errResult != NULL) {
            *errResult = err;
        }
        
        return ids[0];
    }

    static Platform getDefault( 
        cl_int *errResult = NULL )
    {
        return get(errResult);
    }

    
#if defined(CL_VERSION_1_2)
    /**
     * Unload the OpenCL compiler associated with this platform.
     */
    cl_int
    unloadCompiler()
    {
        return ::clUnloadPlatformCompiler(object_);
    }
#endif // #if defined(CL_VERSION_1_2)
};

/**
 * Deprecated APIs for 1.2
 */
#if defined(CL_VERSION_1_1)
/**
 * Unload the OpenCL compiler.
 * \note Deprecated for OpenCL 1.2. Use Platform::unloadCompiler instead.
 */
inline cl_int
UnloadCompiler()
{
    return ::clUnloadCompiler();
}
#endif // #if defined(CL_VERSION_1_1)

class Context : public detail::Wrapper<cl_context>
{
private:
    static volatile int default_initialized_;
    static Context default_;
    static volatile cl_int default_error_;
public:
    Context(
        const VECTOR_CLASS<Device>& devices,
        cl_context_properties* properties = NULL,
        void (CL_CALLBACK * notifyFptr)(
            const char *,
            const void *,
            ::size_t,
            void *) = NULL,
        void* data = NULL,
        cl_int* err = NULL)
    {
        cl_int error;
        object_ = ::clCreateContext(
            properties, (cl_uint) devices.size(),
            (cl_device_id*) &devices.front(),
            notifyFptr, data, &error);

        detail::errHandler(error, __CREATE_CONTEXT_ERR);
        if (err != NULL) {
            *err = error;
        }
    }

    Context(
        const Device& device,
        cl_context_properties* properties = NULL,
        void (CL_CALLBACK * notifyFptr)(
            const char *,
            const void *,
            ::size_t,
            void *) = NULL,
        void* data = NULL,
        cl_int* err = NULL)
    {
        cl_int error;
        object_ = ::clCreateContext(
            properties, 1,
            (cl_device_id*) &device,
            notifyFptr, data, &error);

        detail::errHandler(error, __CREATE_CONTEXT_ERR);
        if (err != NULL) {
            *err = error;
        }
    }

    Context(
        cl_device_type type,
        cl_context_properties* properties = NULL,
        void (CL_CALLBACK * notifyFptr)(
            const char *,
            const void *,
            ::size_t,
            void *) = NULL,
        void* data = NULL,
        cl_int* err = NULL)
    {
        cl_int error;

#if !defined(__APPLE__) || !defined(__MACOS)
        cl_context_properties prop[4] = {CL_CONTEXT_PLATFORM, 0, 0, 0 };	
        if (properties == NULL) {
            prop[1] = (cl_context_properties)Platform::get(&error)();
            if (error != CL_SUCCESS) {
                detail::errHandler(error, __CREATE_CONTEXT_FROM_TYPE_ERR);
                if (err != NULL) {
                    *err = error;
                    return;
                }
            }

            properties = &prop[0];
        }
#endif
        object_ = ::clCreateContextFromType(
            properties, type, notifyFptr, data, &error);

        detail::errHandler(error, __CREATE_CONTEXT_FROM_TYPE_ERR);
        if (err != NULL) {
            *err = error;
        }
    }

    static Context getDefault(cl_int * err = NULL) 
    {
        int state = detail::compare_exchange(
            &default_initialized_, 
            __DEFAULT_BEING_INITIALIZED, __DEFAULT_NOT_INITIALIZED);
        
        if (state & __DEFAULT_INITIALIZED) {
            if (err != NULL) {
                *err = default_error_;
            }
            return default_;
        }

        if (state & __DEFAULT_BEING_INITIALIZED) {
              // Assume writes will propagate eventually...
              while(default_initialized_ != __DEFAULT_INITIALIZED) {
                                ;
              }

            if (err != NULL) {
                *err = default_error_;
            }
            return default_;
        }

        cl_int error;
        default_ = Context(
            CL_DEVICE_TYPE_DEFAULT,
            NULL,
            NULL,
            NULL,
            &error);

        default_error_ = error;
        // Assume writes will propagate eventually...
        default_initialized_ = __DEFAULT_INITIALIZED;
        if (err != NULL) {
            *err = default_error_;
        }
        return default_;

    }

    Context() : detail::Wrapper<cl_type>() { }

    Context(const Context& context) : detail::Wrapper<cl_type>(context) { }

    Context(const cl_context& context) : detail::Wrapper<cl_type>(context) { }

    Context& operator = (const Context& rhs)
    {
        if (this != &rhs) {
            detail::Wrapper<cl_type>::operator=(rhs);
        }
        return *this;
    }

    Context& operator = (const cl_context& rhs)
    {
        detail::Wrapper<cl_type>::operator=(rhs);
        return *this;
    }

    template <typename T>
    cl_int getInfo(cl_context_info name, T* param) const
    {
        return detail::errHandler(
            detail::getInfo(&::clGetContextInfo, object_, name, param),
            __GET_CONTEXT_INFO_ERR);
    }

    template <cl_int name> typename
    detail::param_traits<detail::cl_context_info, name>::param_type
    getInfo(cl_int* err = NULL) const
    {
        typename detail::param_traits<
            detail::cl_context_info, name>::param_type param;
        cl_int result = getInfo(name, &param);
        if (err != NULL) {
            *err = result;
        }
        return param;
    }

    cl_int getSupportedImageFormats(
        cl_mem_flags flags,
        cl_mem_object_type type,
        VECTOR_CLASS<ImageFormat>* formats) const
    {
        cl_uint numEntries;
        cl_int err = ::clGetSupportedImageFormats(
           object_, 
           flags,
           type, 
           0, 
           NULL, 
           &numEntries);
        if (err != CL_SUCCESS) {
            return detail::errHandler(err, __GET_SUPPORTED_IMAGE_FORMATS_ERR);
        }

        ImageFormat* value = (ImageFormat*)
            alloca(numEntries * sizeof(ImageFormat));
        err = ::clGetSupportedImageFormats(
            object_, 
            flags, 
            type, 
            numEntries,
            (cl_image_format*) value, 
            NULL);
        if (err != CL_SUCCESS) {
            return detail::errHandler(err, __GET_SUPPORTED_IMAGE_FORMATS_ERR);
        }

        formats->assign(&value[0], &value[numEntries]);
        return CL_SUCCESS;
    }


#if defined(CL_VERSION_1_2)
    cl_int setPrintfCallback(
        void (CL_CALLBACK * pfn_notify)(
            cl_context /* program */, 
            cl_uint /*printf_data_len */, 
            char * /* printf_data_ptr */, 
            void * /* user_data */),
        void * user_data )
    {
        return detail::errHandler(
            ::clSetPrintfCallback(
                object_,
                pfn_notify,
                user_data), 
            __SET_PRINTF_CALLBACK_ERR);
    }
#endif // #if defined(CL_VERSION_1_2)
};

inline Device Device::getDefault(cl_int * err)
{
    cl_int error;
    Device device;

    Context context = Context::getDefault(&error);
    detail::errHandler(error, __CREATE_COMMAND_QUEUE_ERR);

    if (error != CL_SUCCESS) {
        if (err != NULL) {
            *err = error;
        }
    }
    else {
        device = context.getInfo<CL_CONTEXT_DEVICES>()[0];
        if (err != NULL) {
            *err = CL_SUCCESS;
        }
    }

    return device;
}


#ifdef _WIN32
__declspec(selectany) volatile int Context::default_initialized_ = __DEFAULT_NOT_INITIALIZED;
__declspec(selectany) Context Context::default_;
__declspec(selectany) volatile cl_int Context::default_error_ = CL_SUCCESS;
#else
__attribute__((weak)) volatile int Context::default_initialized_ = __DEFAULT_NOT_INITIALIZED;
__attribute__((weak)) Context Context::default_;
__attribute__((weak)) volatile cl_int Context::default_error_ = CL_SUCCESS;
#endif

__GET_INFO_HELPER_WITH_RETAIN(cl::Context)

/*! \class Event
 * \brief Event interface for cl_event.
 */
class Event : public detail::Wrapper<cl_event>
{
public:
    Event() : detail::Wrapper<cl_type>() { }

    Event(const Event& event) : detail::Wrapper<cl_type>(event) { }

    Event(const cl_event& event) : detail::Wrapper<cl_type>(event) { }

    Event& operator = (const Event& rhs)
    {
        if (this != &rhs) {
            detail::Wrapper<cl_type>::operator=(rhs);
        }
        return *this;
    }

    Event& operator = (const cl_event& rhs)
    {
        detail::Wrapper<cl_type>::operator=(rhs);
        return *this;
    }

    template <typename T>
    cl_int getInfo(cl_event_info name, T* param) const
    {
        return detail::errHandler(
            detail::getInfo(&::clGetEventInfo, object_, name, param),
            __GET_EVENT_INFO_ERR);
    }

    template <cl_int name> typename
    detail::param_traits<detail::cl_event_info, name>::param_type
    getInfo(cl_int* err = NULL) const
    {
        typename detail::param_traits<
            detail::cl_event_info, name>::param_type param;
        cl_int result = getInfo(name, &param);
        if (err != NULL) {
            *err = result;
        }
        return param;
    }

    template <typename T>
    cl_int getProfilingInfo(cl_profiling_info name, T* param) const
    {
        return detail::errHandler(detail::getInfo(
            &::clGetEventProfilingInfo, object_, name, param),
            __GET_EVENT_PROFILE_INFO_ERR);
    }

    template <cl_int name> typename
    detail::param_traits<detail::cl_profiling_info, name>::param_type
    getProfilingInfo(cl_int* err = NULL) const
    {
        typename detail::param_traits<
            detail::cl_profiling_info, name>::param_type param;
        cl_int result = getProfilingInfo(name, &param);
        if (err != NULL) {
            *err = result;
        }
        return param;
    }

    cl_int wait() const
    {
        return detail::errHandler(
            ::clWaitForEvents(1, &object_),
            __WAIT_FOR_EVENTS_ERR);
    }

#if defined(CL_VERSION_1_1)
    cl_int setCallback(
        cl_int type,
        void (CL_CALLBACK * pfn_notify)(cl_event, cl_int, void *),		
        void * user_data = NULL)
    {
        return detail::errHandler(
            ::clSetEventCallback(
                object_,
                type,
                pfn_notify,
                user_data), 
            __SET_EVENT_CALLBACK_ERR);
    }
#endif

    static cl_int
    waitForEvents(const VECTOR_CLASS<Event>& events)
    {
        return detail::errHandler(
            ::clWaitForEvents(
                (cl_uint) events.size(), (cl_event*)&events.front()),
            __WAIT_FOR_EVENTS_ERR);
    }
};

__GET_INFO_HELPER_WITH_RETAIN(cl::Event)

#if defined(CL_VERSION_1_1)
/*! \class UserEvent
 * \brief User event interface for cl_event.
 */
class UserEvent : public Event
{
public:
    UserEvent(
        const Context& context,
        cl_int * err = NULL)
    {
        cl_int error;
        object_ = ::clCreateUserEvent(
            context(),
            &error);

        detail::errHandler(error, __CREATE_USER_EVENT_ERR);
        if (err != NULL) {
            *err = error;
        }
    }

    UserEvent() : Event() { }

    UserEvent(const UserEvent& event) : Event(event) { }

    UserEvent& operator = (const UserEvent& rhs)
    {
        if (this != &rhs) {
            Event::operator=(rhs);
        }
        return *this;
    }

    cl_int setStatus(cl_int status)
    {
        return detail::errHandler(
            ::clSetUserEventStatus(object_,status), 
            __SET_USER_EVENT_STATUS_ERR);
    }
};
#endif

inline static cl_int
WaitForEvents(const VECTOR_CLASS<Event>& events)
{
    return detail::errHandler(
        ::clWaitForEvents(
            (cl_uint) events.size(), (cl_event*)&events.front()),
        __WAIT_FOR_EVENTS_ERR);
}

/*! \class Memory
 * \brief Memory interface for cl_mem.
 */
class Memory : public detail::Wrapper<cl_mem>
{
public:
    Memory() : detail::Wrapper<cl_type>() { }

    Memory(const Memory& memory) : detail::Wrapper<cl_type>(memory) { }

    Memory(const cl_mem& memory) : detail::Wrapper<cl_type>(memory) { }

    Memory& operator = (const Memory& rhs)
    {
        if (this != &rhs) {
            detail::Wrapper<cl_type>::operator=(rhs);
        }
        return *this;
    }

    Memory& operator = (const cl_mem& rhs)
    {
        detail::Wrapper<cl_type>::operator=(rhs);
        return *this;
    }

    template <typename T>
    cl_int getInfo(cl_mem_info name, T* param) const
    {
        return detail::errHandler(
            detail::getInfo(&::clGetMemObjectInfo, object_, name, param),
            __GET_MEM_OBJECT_INFO_ERR);
    }

    template <cl_int name> typename
    detail::param_traits<detail::cl_mem_info, name>::param_type
    getInfo(cl_int* err = NULL) const
    {
        typename detail::param_traits<
            detail::cl_mem_info, name>::param_type param;
        cl_int result = getInfo(name, &param);
        if (err != NULL) {
            *err = result;
        }
        return param;
    }

#if defined(CL_VERSION_1_1)
    cl_int setDestructorCallback(
        void (CL_CALLBACK * pfn_notify)(cl_mem, void *),		
        void * user_data = NULL)
    {
        return detail::errHandler(
            ::clSetMemObjectDestructorCallback(
                object_,
                pfn_notify,
                user_data), 
            __SET_MEM_OBJECT_DESTRUCTOR_CALLBACK_ERR);
    }
#endif

};

__GET_INFO_HELPER_WITH_RETAIN(cl::Memory)

// Pre-declare copy functions
class Buffer;
template< typename IteratorType >
cl_int copy( IteratorType startIterator, IteratorType endIterator, cl::Buffer &buffer );
template< typename IteratorType >
cl_int copy( cl::Buffer &buffer, IteratorType startIterator, IteratorType endIterator );

/*! \class Buffer
 * \brief Memory buffer interface.
 */
class Buffer : public Memory
{
public:
    Buffer(
        const Context& context,
        cl_mem_flags flags,
        ::size_t size,
        void* host_ptr = NULL,
        cl_int* err = NULL)
    {
        cl_int error;
        object_ = ::clCreateBuffer(context(), flags, size, host_ptr, &error);

        detail::errHandler(error, __CREATE_BUFFER_ERR);
        if (err != NULL) {
            *err = error;
        }
    }

    Buffer(
         cl_mem_flags flags,
        ::size_t size,
        void* host_ptr = NULL,
        cl_int* err = NULL)
    {
        cl_int error;

        Context context = Context::getDefault(err);

        object_ = ::clCreateBuffer(context(), flags, size, host_ptr, &error);

        detail::errHandler(error, __CREATE_BUFFER_ERR);
        if (err != NULL) {
            *err = error;
        }
    }

    /**
     * Construct a Buffer from a host container via iterators.
     */
    template< typename IteratorType >
    Buffer(
        IteratorType startIterator,
        IteratorType endIterator,
        bool readOnly = false,
        bool useHostPtr = true,
        cl_int* err = NULL)
    {
        typedef typename std::iterator_traits<IteratorType>::value_type DataType;

        cl_int error;
        cl_mem_flags flags = 0;
        if( readOnly ) {
            flags |= CL_MEM_READ_ONLY;
        }
        else {
            flags |= CL_MEM_READ_WRITE;
        }
        if( useHostPtr ) {
            flags |= CL_MEM_USE_HOST_PTR;
        }
        
        ::size_t size = sizeof(DataType)*(endIterator - startIterator);

        Context context = Context::getDefault(err);

        object_ = ::clCreateBuffer(context(), flags, size, static_cast<DataType*>(&*startIterator), &error);

        detail::errHandler(error, __CREATE_BUFFER_ERR);
        if (err != NULL) {
            *err = error;
        }

        if( !useHostPtr ) {
            error = cl::copy(startIterator, endIterator, *this);
            detail::errHandler(error, __CREATE_BUFFER_ERR);
            if (err != NULL) {
                *err = error;
            }
        }
    }


    Buffer() : Memory() { }

    Buffer(const Buffer& buffer) : Memory(buffer) { }

    Buffer(const cl_mem& buffer) : Memory(buffer) { }

    Buffer& operator = (const Buffer& rhs)
    {
        if (this != &rhs) {
            Memory::operator=(rhs);
        }
        return *this;
    }

    Buffer& operator = (const cl_mem& rhs)
    {
        Memory::operator=(rhs);
        return *this;
    }

#if defined(CL_VERSION_1_1)
    Buffer createSubBuffer(
        cl_mem_flags flags,
        cl_buffer_create_type buffer_create_type,
        const void * buffer_create_info,
        cl_int * err = NULL)
    {
        Buffer result;
        cl_int error;
        result.object_ = ::clCreateSubBuffer(
            object_, 
            flags, 
            buffer_create_type, 
            buffer_create_info, 
            &error);

        detail::errHandler(error, __CREATE_SUBBUFFER_ERR);
        if (err != NULL) {
            *err = error;
        }

        return result;
    }		
#endif
};

#if defined (USE_DX_INTEROP)
class BufferD3D10 : public Buffer
{
public:
    typedef CL_API_ENTRY cl_mem (CL_API_CALL *PFN_clCreateFromD3D10BufferKHR)(
    cl_context context, cl_mem_flags flags, ID3D10Buffer*  buffer,
    cl_int* errcode_ret);

    BufferD3D10(
        const Context& context,
        cl_mem_flags flags,
        ID3D10Buffer* bufobj,
        cl_int * err = NULL)
    {
        static PFN_clCreateFromD3D10BufferKHR pfn_clCreateFromD3D10BufferKHR = NULL;

#if defined(CL_VERSION_1_2)
        vector<cl_context_properties> props = context.getInfo<CL_CONTEXT_PROPERTIES>();
        cl_platform platform = -1;
        for( int i = 0; i < props.size(); ++i ) {
            if( props[i] == CL_CONTEXT_PLATFORM ) {
                platform = props[i+1];
            }
        }
        __INIT_CL_EXT_FCN_PTR_PLATFORM(platform, clCreateFromD3D10BufferKHR);
#endif
#if defined(CL_VERSION_1_1)
        __INIT_CL_EXT_FCN_PTR(clCreateFromD3D10BufferKHR);
#endif

        cl_int error;
        object_ = pfn_clCreateFromD3D10BufferKHR(
            context(),
            flags,
            bufobj,
            &error);

        detail::errHandler(error, __CREATE_GL_BUFFER_ERR);
        if (err != NULL) {
            *err = error;
        }
    }

    BufferD3D10() : Buffer() { }

    BufferD3D10(const BufferD3D10& buffer) : Buffer(buffer) { }

    BufferD3D10(const cl_mem& buffer) : Buffer(buffer) { }

    BufferD3D10& operator = (const BufferD3D10& rhs)
    {
        if (this != &rhs) {
            Buffer::operator=(rhs);
        }
        return *this;
    }

    BufferD3D10& operator = (const cl_mem& rhs)
    {
        Buffer::operator=(rhs);
        return *this;
    }
};
#endif

/*! \class BufferGL
 * \brief Memory buffer interface for GL interop.
 */
class BufferGL : public Buffer
{
public:
    BufferGL(
        const Context& context,
        cl_mem_flags flags,
        GLuint bufobj,
        cl_int * err = NULL)
    {
        cl_int error;
        object_ = ::clCreateFromGLBuffer(
            context(),
            flags,
            bufobj,
            &error);

        detail::errHandler(error, __CREATE_GL_BUFFER_ERR);
        if (err != NULL) {
            *err = error;
        }
    }

    BufferGL() : Buffer() { }

    BufferGL(const BufferGL& buffer) : Buffer(buffer) { }

    BufferGL(const cl_mem& buffer) : Buffer(buffer) { }

    BufferGL& operator = (const BufferGL& rhs)
    {
        if (this != &rhs) {
            Buffer::operator=(rhs);
        }
        return *this;
    }

    BufferGL& operator = (const cl_mem& rhs)
    {
        Buffer::operator=(rhs);
        return *this;
    }

    cl_int getObjectInfo(
        cl_gl_object_type *type,
        GLuint * gl_object_name)
    {
        return detail::errHandler(
            ::clGetGLObjectInfo(object_,type,gl_object_name),
            __GET_GL_OBJECT_INFO_ERR);
    }
};

/*! \class BufferRenderGL
 * \brief Memory buffer interface for GL interop with renderbuffer.
 */
class BufferRenderGL : public Buffer
{
public:
    BufferRenderGL(
        const Context& context,
        cl_mem_flags flags,
        GLuint bufobj,
        cl_int * err = NULL)
    {
        cl_int error;
        object_ = ::clCreateFromGLRenderbuffer(
            context(),
            flags,
            bufobj,
            &error);

        detail::errHandler(error, __CREATE_GL_RENDER_BUFFER_ERR);
        if (err != NULL) {
            *err = error;
        }
    }

    BufferRenderGL() : Buffer() { }

    BufferRenderGL(const BufferGL& buffer) : Buffer(buffer) { }

    BufferRenderGL(const cl_mem& buffer) : Buffer(buffer) { }

    BufferRenderGL& operator = (const BufferRenderGL& rhs)
    {
        if (this != &rhs) {
            Buffer::operator=(rhs);
        }
        return *this;
    }

    BufferRenderGL& operator = (const cl_mem& rhs)
    {
        Buffer::operator=(rhs);
        return *this;
    }

    cl_int getObjectInfo(
        cl_gl_object_type *type,
        GLuint * gl_object_name)
    {
        return detail::errHandler(
            ::clGetGLObjectInfo(object_,type,gl_object_name),
            __GET_GL_OBJECT_INFO_ERR);
    }
};

/*! \class Image
 * \brief Base class  interface for all images.
 */
class Image : public Memory
{
protected:
    Image() : Memory() { }

    Image(const Image& image) : Memory(image) { }

    Image(const cl_mem& image) : Memory(image) { }

    Image& operator = (const Image& rhs)
    {
        if (this != &rhs) {
            Memory::operator=(rhs);
        }
        return *this;
    }

    Image& operator = (const cl_mem& rhs)
    {
        Memory::operator=(rhs);
        return *this;
    }

public:
    template <typename T>
    cl_int getImageInfo(cl_image_info name, T* param) const
    {
        return detail::errHandler(
            detail::getInfo(&::clGetImageInfo, object_, name, param),
            __GET_IMAGE_INFO_ERR);
    }

    template <cl_int name> typename
    detail::param_traits<detail::cl_image_info, name>::param_type
    getImageInfo(cl_int* err = NULL) const
    {
        typename detail::param_traits<
            detail::cl_image_info, name>::param_type param;
        cl_int result = getImageInfo(name, &param);
        if (err != NULL) {
            *err = result;
        }
        return param;
    }
};

#if defined(CL_VERSION_1_2)
/*! \class Image1D
 * \brief Image interface for 1D images.
 */
class Image1D : public Image
{
public:
    Image1D(
        const Context& context,
        cl_mem_flags flags,
        ImageFormat format,
        ::size_t width,
        void* host_ptr = NULL,
        cl_int* err = NULL)
    {
        cl_int error;
        cl_image_desc desc;
        desc.image_type = CL_MEM_OBJECT_IMAGE1D;
        desc.image_width = width;
        desc.image_row_pitch = 0;
        desc.num_mip_levels = 0;
        desc.num_samples = 0;
        desc.buffer = 0;
        object_ = ::clCreateImage(
            context(), 
            flags, 
            &format, 
            &desc, 
            host_ptr, 
            &error);

        detail::errHandler(error, __CREATE_IMAGE_ERR);
        if (err != NULL) {
            *err = error;
        }
    }

    Image1D() { }

    Image1D(const Image1D& image1D) : Image(image1D) { }

    Image1D(const cl_mem& image1D) : Image(image1D) { }

    Image1D& operator = (const Image1D& rhs)
    {
        if (this != &rhs) {
            Image::operator=(rhs);
        }
        return *this;
    }

    Image1D& operator = (const cl_mem& rhs)
    {
        Image::operator=(rhs);
        return *this;
    }
};

/*! \class Image1DBuffer
 * \brief Image interface for 1D buffer images.
 */
class Image1DBuffer : public Image
{
public:
    Image1DBuffer(
        const Context& context,
        cl_mem_flags flags,
        ImageFormat format,
        ::size_t width,
        Buffer &buffer,
        void* host_ptr = NULL,
        cl_int* err = NULL)
    {
        cl_int error;
        cl_image_desc desc;
        desc.image_type = CL_MEM_OBJECT_IMAGE1D_BUFFER;
        desc.image_width = width;
        desc.image_row_pitch = 0;
        desc.num_mip_levels = 0;
        desc.num_samples = 0;
        desc.buffer = buffer();
        object_ = ::clCreateImage(
            context(), 
            flags, 
            &format, 
            &desc, 
            host_ptr, 
            &error);

        detail::errHandler(error, __CREATE_IMAGE_ERR);
        if (err != NULL) {
            *err = error;
        }
    }

    Image1DBuffer() { }

    Image1DBuffer(const Image1DBuffer& image1D) : Image(image1D) { }

    Image1DBuffer(const cl_mem& image1D) : Image(image1D) { }

    Image1DBuffer& operator = (const Image1DBuffer& rhs)
    {
        if (this != &rhs) {
            Image::operator=(rhs);
        }
        return *this;
    }

    Image1DBuffer& operator = (const cl_mem& rhs)
    {
        Image::operator=(rhs);
        return *this;
    }
};

/*! \class Image1DArray
 * \brief Image interface for arrays of 1D images.
 */
class Image1DArray : public Image
{
public:
    Image1DArray(
        const Context& context,
        cl_mem_flags flags,
        ImageFormat format,
        ::size_t arraySize,
        ::size_t width,
        ::size_t rowPitch,
        void* host_ptr = NULL,
        cl_int* err = NULL)
    {
        cl_int error;
        cl_image_desc desc;
        desc.image_type = CL_MEM_OBJECT_IMAGE1D_ARRAY;
        desc.image_array_size = arraySize;
        desc.image_width = width;
        desc.image_row_pitch = rowPitch;
        desc.num_mip_levels = 0;
        desc.num_samples = 0;
        desc.buffer = 0;
        object_ = ::clCreateImage(
            context(), 
            flags, 
            &format, 
            &desc, 
            host_ptr, 
            &error);

        detail::errHandler(error, __CREATE_IMAGE_ERR);
        if (err != NULL) {
            *err = error;
        }
    }

    Image1DArray() { }

    Image1DArray(const Image1DArray& imageArray) : Image(imageArray) { }

    Image1DArray(const cl_mem& imageArray) : Image(imageArray) { }

    Image1DArray& operator = (const Image1DArray& rhs)
    {
        if (this != &rhs) {
            Image::operator=(rhs);
        }
        return *this;
    }

    Image1DArray& operator = (const cl_mem& rhs)
    {
        Image::operator=(rhs);
        return *this;
    }
};
#endif // #if defined(CL_VERSION_1_2)


/*! \class Image2D
 * \brief Image interface for 2D images.
 */
class Image2D : public Image
{
public:
    Image2D(
        const Context& context,
        cl_mem_flags flags,
        ImageFormat format,
        ::size_t width,
        ::size_t height,
        ::size_t row_pitch = 0,
        void* host_ptr = NULL,
        cl_int* err = NULL)
    {
        cl_int error;
#if defined(CL_VERSION_1_2)
        cl_image_desc desc;
        desc.image_type = CL_MEM_OBJECT_IMAGE2D;
        desc.image_width = width;
        desc.image_height = height;
        desc.image_row_pitch = row_pitch;
        desc.num_mip_levels = 0;
        desc.num_samples = 0;
        desc.buffer = 0;
        object_ = ::clCreateImage(
            context(), 
            flags, 
            &format, 
            &desc, 
            host_ptr, 
            &error);

        detail::errHandler(error, __CREATE_IMAGE_ERR);
        if (err != NULL) {
            *err = error;
        }
#else
        object_ = ::clCreateImage2D(
            context(), flags,&format, width, height, row_pitch, host_ptr, &error);

        detail::errHandler(error, __CREATE_IMAGE2D_ERR);
        if (err != NULL) {
            *err = error;
        }
#endif // #if defined(CL_VERSION_1_2)
    }

    Image2D() { }

    Image2D(const Image2D& image2D) : Image(image2D) { }

    Image2D(const cl_mem& image2D) : Image(image2D) { }

    Image2D& operator = (const Image2D& rhs)
    {
        if (this != &rhs) {
            Image::operator=(rhs);
        }
        return *this;
    }

    Image2D& operator = (const cl_mem& rhs)
    {
        Image::operator=(rhs);
        return *this;
    }
};


#if !defined(CL_VERSION_1_2)
/*! \class Image2DGL
 * \brief 2D image interface for GL interop.
 * \note Deprecated for OpenCL 1.2. Please use ImageGL instead.
 */
class Image2DGL : public Image2D
{
public:
    Image2DGL(
        const Context& context,
        cl_mem_flags flags,
        GLenum target,
        GLint  miplevel,
        GLuint texobj,
        cl_int * err = NULL)
    {
        cl_int error;
        object_ = ::clCreateFromGLTexture2D(
            context(),
            flags,
            target,
            miplevel,
            texobj,
            &error);

        detail::errHandler(error, __CREATE_GL_TEXTURE_2D_ERR);
        if (err != NULL) {
            *err = error;
        }

    }

    Image2DGL() : Image2D() { }

    Image2DGL(const Image2DGL& image) : Image2D(image) { }

    Image2DGL(const cl_mem& image) : Image2D(image) { }

    Image2DGL& operator = (const Image2DGL& rhs)
    {
        if (this != &rhs) {
            Image2D::operator=(rhs);
        }
        return *this;
    }

    Image2DGL& operator = (const cl_mem& rhs)
    {
        Image2D::operator=(rhs);
        return *this;
    }
};
#endif // #if !defined(CL_VERSION_1_2)

#if defined(CL_VERSION_1_2)
/*! \class Image2DArray
 * \brief Image interface for arrays of 2D images.
 */
class Image2DArray : public Image
{
public:
    Image2DArray(
        const Context& context,
        cl_mem_flags flags,
        ImageFormat format,
        ::size_t arraySize,
        ::size_t width,
        ::size_t height,
        ::size_t rowPitch,
        ::size_t slicePitch,
        void* host_ptr = NULL,
        cl_int* err = NULL)
    {
        cl_int error;
        cl_image_desc desc;
        desc.image_type = CL_MEM_OBJECT_IMAGE2D_ARRAY;
        desc.image_array_size = arraySize;
        desc.image_width = width;
        desc.image_height = height;
        desc.image_row_pitch = rowPitch;
        desc.image_slice_pitch = slicePitch;
        desc.num_mip_levels = 0;
        desc.num_samples = 0;
        desc.buffer = 0;
        object_ = ::clCreateImage(
            context(), 
            flags, 
            &format, 
            &desc, 
            host_ptr, 
            &error);

        detail::errHandler(error, __CREATE_IMAGE_ERR);
        if (err != NULL) {
            *err = error;
        }
    }

    Image2DArray() { }

    Image2DArray(const Image2DArray& imageArray) : Image(imageArray) { }

    Image2DArray(const cl_mem& imageArray) : Image(imageArray) { }

    Image2DArray& operator = (const Image2DArray& rhs)
    {
        if (this != &rhs) {
            Image::operator=(rhs);
        }
        return *this;
    }

    Image2DArray& operator = (const cl_mem& rhs)
    {
        Image::operator=(rhs);
        return *this;
    }
};
#endif // #if defined(CL_VERSION_1_2)

/*! \class Image3D
 * \brief Image interface for 3D images.
 */
class Image3D : public Image
{
public:
    Image3D(
        const Context& context,
        cl_mem_flags flags,
        ImageFormat format,
        ::size_t width,
        ::size_t height,
        ::size_t depth,
        ::size_t row_pitch = 0,
        ::size_t slice_pitch = 0,
        void* host_ptr = NULL,
        cl_int* err = NULL)
    {
        cl_int error;
#if defined(CL_VERSION_1_2)
        cl_image_desc desc;
        desc.image_type = CL_MEM_OBJECT_IMAGE2D;
        desc.image_width = width;
        desc.image_height = height;
        desc.image_depth = depth;
        desc.image_row_pitch = row_pitch;
        desc.image_slice_pitch = slice_pitch;
        desc.num_mip_levels = 0;
        desc.num_samples = 0;
        desc.buffer = 0;
        object_ = ::clCreateImage(
            context(), 
            flags, 
            &format, 
            &desc, 
            host_ptr, 
            &error);

        detail::errHandler(error, __CREATE_IMAGE_ERR);
        if (err != NULL) {
            *err = error;
        }
#else  // #if defined(CL_VERSION_1_2)
        object_ = ::clCreateImage3D(
            context(), flags, &format, width, height, depth, row_pitch,
            slice_pitch, host_ptr, &error);

        detail::errHandler(error, __CREATE_IMAGE3D_ERR);
        if (err != NULL) {
            *err = error;
        }
#endif // #if defined(CL_VERSION_1_2)
    }

    Image3D() { }

    Image3D(const Image3D& image3D) : Image(image3D) { }

    Image3D(const cl_mem& image3D) : Image(image3D) { }

    Image3D& operator = (const Image3D& rhs)
    {
        if (this != &rhs) {
            Image::operator=(rhs);
        }
        return *this;
    }

    Image3D& operator = (const cl_mem& rhs)
    {
        Image::operator=(rhs);
        return *this;
    }
};

#if !defined(CL_VERSION_1_2)
/*! \class Image2DGL
 * \brief 2D image interface for GL interop.
 * \note Deprecated for OpenCL 1.2. Please use ImageGL instead.
 */
class Image3DGL : public Image3D
{
public:
    Image3DGL(
        const Context& context,
        cl_mem_flags flags,
        GLenum target,
        GLint  miplevel,
        GLuint texobj,
        cl_int * err = NULL)
    {
        cl_int error;
        object_ = ::clCreateFromGLTexture3D(
            context(),
            flags,
            target,
            miplevel,
            texobj,
            &error);

        detail::errHandler(error, __CREATE_GL_TEXTURE_3D_ERR);
        if (err != NULL) {
            *err = error;
        }
    }

    Image3DGL() : Image3D() { }

    Image3DGL(const Image3DGL& image) : Image3D(image) { }

    Image3DGL(const cl_mem& image) : Image3D(image) { }

    Image3DGL& operator = (const Image3DGL& rhs)
    {
        if (this != &rhs) {
            Image3D::operator=(rhs);
        }
        return *this;
    }

    Image3DGL& operator = (const cl_mem& rhs)
    {
        Image3D::operator=(rhs);
        return *this;
    }
};
#endif // #if !defined(CL_VERSION_1_2)

#if defined(CL_VERSION_1_2)
/*! \class ImageGL
 * \brief general image interface for GL interop.
 * We abstract the 2D and 3D GL images into a single instance here
 * that wraps all GL sourced images on the grounds that setup information
 * was performed by OpenCL anyway.
 */
class ImageGL : public Image
{
public:
    ImageGL(
        const Context& context,
        cl_mem_flags flags,
        GLenum target,
        GLint  miplevel,
        GLuint texobj,
        cl_int * err = NULL)
    {
        cl_int error;
        object_ = ::clCreateFromGLTexture(
            context(), 
            flags, 
            target,
            miplevel,
            texobj,
            &error);

        detail::errHandler(error, __CREATE_GL_TEXTURE_ERR);
        if (err != NULL) {
            *err = error;
        }
    }

    ImageGL() : Image() { }

    ImageGL(const ImageGL& image) : Image(image) { }

    ImageGL(const cl_mem& image) : Image(image) { }

    ImageGL& operator = (const ImageGL& rhs)
    {
        if (this != &rhs) {
            ImageGL::operator=(rhs);
        }
        return *this;
    }

    ImageGL& operator = (const cl_mem& rhs)
    {
        ImageGL::operator=(rhs);
        return *this;
    }
};
#endif // #if defined(CL_VERSION_1_2)

/*! \class Sampler
 * \brief Sampler interface for cl_sampler.
 */
class Sampler : public detail::Wrapper<cl_sampler>
{
public:
    Sampler() { }

    Sampler(
        const Context& context,
        cl_bool normalized_coords,
        cl_addressing_mode addressing_mode,
        cl_filter_mode filter_mode,
        cl_int* err = NULL)
    {
        cl_int error;
        object_ = ::clCreateSampler(
            context(), 
            normalized_coords,
            addressing_mode,
            filter_mode,
            &error);

        detail::errHandler(error, __CREATE_SAMPLER_ERR);
        if (err != NULL) {
            *err = error;
        }
    }

    Sampler(const Sampler& sampler) : detail::Wrapper<cl_type>(sampler) { }

    Sampler(const cl_sampler& sampler) : detail::Wrapper<cl_type>(sampler) { }

    Sampler& operator = (const Sampler& rhs)
    {
        if (this != &rhs) {
            detail::Wrapper<cl_type>::operator=(rhs);
        }
        return *this;
    }

    Sampler& operator = (const cl_sampler& rhs)
    {
        detail::Wrapper<cl_type>::operator=(rhs);
        return *this;
    }

    template <typename T>
    cl_int getInfo(cl_sampler_info name, T* param) const
    {
        return detail::errHandler(
            detail::getInfo(&::clGetSamplerInfo, object_, name, param),
            __GET_SAMPLER_INFO_ERR);
    }

    template <cl_int name> typename
    detail::param_traits<detail::cl_sampler_info, name>::param_type
    getInfo(cl_int* err = NULL) const
    {
        typename detail::param_traits<
            detail::cl_sampler_info, name>::param_type param;
        cl_int result = getInfo(name, &param);
        if (err != NULL) {
            *err = result;
        }
        return param;
    }
};

__GET_INFO_HELPER_WITH_RETAIN(cl::Sampler)

class Program;
class CommandQueue;
class Kernel;

/*! \class NDRange
 * \brief NDRange interface
 */
class NDRange
{
private:
    size_t<3> sizes_;
    cl_uint dimensions_;

public:
    NDRange()
        : dimensions_(0)
    { }

    NDRange(::size_t size0)
        : dimensions_(1)
    {
        sizes_.push_back(size0);
    }

    NDRange(::size_t size0, ::size_t size1)
        : dimensions_(2)
    {
        sizes_.push_back(size0);
        sizes_.push_back(size1);
    }

    NDRange(::size_t size0, ::size_t size1, ::size_t size2)
        : dimensions_(3)
    {
        sizes_.push_back(size0);
        sizes_.push_back(size1);
        sizes_.push_back(size2);
    }

    operator const ::size_t*() const { return (const ::size_t*) sizes_; }
    ::size_t dimensions() const { return dimensions_; }
};

static const NDRange NullRange;

/*!
 * \struct LocalSpaceArg
 * \brief Local address raper for use with Kernel::setArg
 */
struct LocalSpaceArg
{
    ::size_t size_;
};

namespace detail {

template <typename T>
struct KernelArgumentHandler
{
    static ::size_t size(const T&) { return sizeof(T); }
    static T* ptr(T& value) { return &value; }
};

template <>
struct KernelArgumentHandler<LocalSpaceArg>
{
    static ::size_t size(const LocalSpaceArg& value) { return value.size_; }
    static void* ptr(LocalSpaceArg&) { return NULL; }
};

} 
//! \endcond

/*! __local
 * \brief Returns a Local type to support setArg overloading.
 * Deprecated. Replaced with Local.
 */
inline LocalSpaceArg
__local(::size_t size)
{
    LocalSpaceArg ret = { size };
    return ret;
}

/*! Local
 * \brief Returns a Local type to support setArg overloading.
 */
inline LocalSpaceArg
Local(::size_t size)
{
    LocalSpaceArg ret = { size };
    return ret;
}

//class KernelFunctor;

/*! \class Kernel
 * \brief Kernel interface that implements cl_kernel
 */
class Kernel : public detail::Wrapper<cl_kernel>
{
public:
    inline Kernel(const Program& program, const char* name, cl_int* err = NULL);

    Kernel() { }

    Kernel(const Kernel& kernel) : detail::Wrapper<cl_type>(kernel) { }

    Kernel(const cl_kernel& kernel) : detail::Wrapper<cl_type>(kernel) { }

    Kernel& operator = (const Kernel& rhs)
    {
        if (this != &rhs) {
            detail::Wrapper<cl_type>::operator=(rhs);
        }
        return *this;
    }

    Kernel& operator = (const cl_kernel& rhs)
    {
        detail::Wrapper<cl_type>::operator=(rhs);
        return *this;
    }

    template <typename T>
    cl_int getInfo(cl_kernel_info name, T* param) const
    {
        return detail::errHandler(
            detail::getInfo(&::clGetKernelInfo, object_, name, param),
            __GET_KERNEL_INFO_ERR);
    }

    template <cl_int name> typename
    detail::param_traits<detail::cl_kernel_info, name>::param_type
    getInfo(cl_int* err = NULL) const
    {
        typename detail::param_traits<
            detail::cl_kernel_info, name>::param_type param;
        cl_int result = getInfo(name, &param);
        if (err != NULL) {
            *err = result;
        }
        return param;
    }

#if defined(CL_VERSION_1_2)
    template <typename T>
    cl_int getArgInfo(cl_uint argIndex, cl_kernel_arg_info name, T* param) const
    {
        return detail::errHandler(
            detail::getInfo(&::clGetKernelArgInfo, object_, argIndex, name, param),
            __GET_KERNEL_ARG_INFO_ERR);
    }

    template <cl_int name> typename
    detail::param_traits<detail::cl_kernel_arg_info, name>::param_type
    getArgInfo(cl_uint argIndex, cl_int* err = NULL) const
    {
        typename detail::param_traits<
            detail::cl_kernel_arg_info, name>::param_type param;
        cl_int result = getArgInfo(name, argIndex, &param);
        if (err != NULL) {
            *err = result;
        }
        return param;
    }
#endif // #if defined(CL_VERSION_1_2)

    template <typename T>
    cl_int getWorkGroupInfo(
        const Device& device, cl_kernel_work_group_info name, T* param) const
    {
        return detail::errHandler(
            detail::getInfo(
                &::clGetKernelWorkGroupInfo, object_, device(), name, param),
                __GET_KERNEL_WORK_GROUP_INFO_ERR);
    }

    template <cl_int name> typename
    detail::param_traits<detail::cl_kernel_work_group_info, name>::param_type
        getWorkGroupInfo(const Device& device, cl_int* err = NULL) const
    {
        typename detail::param_traits<
        detail::cl_kernel_work_group_info, name>::param_type param;
        cl_int result = getWorkGroupInfo(device, name, &param);
        if (err != NULL) {
            *err = result;
        }
        return param;
    }

    template <typename T>
    cl_int setArg(cl_uint index, T value)
    {
        return detail::errHandler(
            ::clSetKernelArg(
                object_,
                index,
                detail::KernelArgumentHandler<T>::size(value),
                detail::KernelArgumentHandler<T>::ptr(value)),
            __SET_KERNEL_ARGS_ERR);
    }

    cl_int setArg(cl_uint index, ::size_t size, void* argPtr)
    {
        return detail::errHandler(
            ::clSetKernelArg(object_, index, size, argPtr),
            __SET_KERNEL_ARGS_ERR);
    }
};

__GET_INFO_HELPER_WITH_RETAIN(cl::Kernel)

/*! \class Program
 * \brief Program interface that implements cl_program.
 */
class Program : public detail::Wrapper<cl_program>
{
public:
    typedef VECTOR_CLASS<std::pair<const void*, ::size_t> > Binaries;
    typedef VECTOR_CLASS<std::pair<const char*, ::size_t> > Sources;

    Program(
        const STRING_CLASS& source,
        cl_int* err = NULL)
    {
        cl_int error;

        const char * strings = source.c_str();
        const ::size_t length  = source.size();

        Context context = Context::getDefault(err);

        object_ = ::clCreateProgramWithSource(
            context(), (cl_uint)1, &strings, &length, &error);

        detail::errHandler(error, __CREATE_PROGRAM_WITH_SOURCE_ERR);

        if (error == CL_SUCCESS) {

            error = ::clBuildProgram(
                object_,
                0,
                NULL,
                "",
                NULL,
                NULL);

            detail::errHandler(error, __BUILD_PROGRAM_ERR);
        }

        if (err != NULL) {
            *err = error;
        }
    }

    Program(
        const STRING_CLASS& source,
		bool build,
        cl_int* err = NULL)
    {
        cl_int error;

        const char * strings = source.c_str();
        const ::size_t length  = source.size();

        Context context = Context::getDefault(err);

        object_ = ::clCreateProgramWithSource(
            context(), (cl_uint)1, &strings, &length, &error);

        detail::errHandler(error, __CREATE_PROGRAM_WITH_SOURCE_ERR);

        if (error == CL_SUCCESS && build) {

            error = ::clBuildProgram(
                object_,
                0,
                NULL,
                "",
                NULL,
                NULL);

            detail::errHandler(error, __BUILD_PROGRAM_ERR);
        }

        if (err != NULL) {
            *err = error;
        }
    }

    Program(
        const Context& context,
        const STRING_CLASS& source,
        bool build = false,
        cl_int* err = NULL)
    {
        cl_int error;

        const char * strings = source.c_str();
        const ::size_t length  = source.size();

        object_ = ::clCreateProgramWithSource(
            context(), (cl_uint)1, &strings, &length, &error);

        detail::errHandler(error, __CREATE_PROGRAM_WITH_SOURCE_ERR);

        if (error == CL_SUCCESS && build) {

            error = ::clBuildProgram(
                object_,
                0,
                NULL,
                "",
                NULL,
                NULL);

            detail::errHandler(error, __BUILD_PROGRAM_ERR);
        }

        if (err != NULL) {
            *err = error;
        }
    }

    Program(
        const Context& context,
        const Sources& sources,
        cl_int* err = NULL)
    {
        cl_int error;

        const ::size_t n = (::size_t)sources.size();
        ::size_t* lengths = (::size_t*) alloca(n * sizeof(::size_t));
        const char** strings = (const char**) alloca(n * sizeof(const char*));

        for (::size_t i = 0; i < n; ++i) {
            strings[i] = sources[(int)i].first;
            lengths[i] = sources[(int)i].second;
        }

        object_ = ::clCreateProgramWithSource(
            context(), (cl_uint)n, strings, lengths, &error);

        detail::errHandler(error, __CREATE_PROGRAM_WITH_SOURCE_ERR);
        if (err != NULL) {
            *err = error;
        }
    }

    Program(
        const Context& context,
        const VECTOR_CLASS<Device>& devices,
        const Binaries& binaries,
        VECTOR_CLASS<cl_int>* binaryStatus = NULL,
        cl_int* err = NULL)
    {
        cl_int error;
        const ::size_t n = binaries.size();
        ::size_t* lengths = (::size_t*) alloca(n * sizeof(::size_t));
        const unsigned char** images = (const unsigned char**) alloca(n * sizeof(const void*));

        for (::size_t i = 0; i < n; ++i) {
            images[i] = (const unsigned char*)binaries[(int)i].first;
            lengths[i] = binaries[(int)i].second;
        }

        object_ = ::clCreateProgramWithBinary(
            context(), (cl_uint) devices.size(),
            (cl_device_id*)&devices.front(),
            lengths, images, binaryStatus != NULL
               ? (cl_int*) &binaryStatus->front()
               : NULL, &error);

        detail::errHandler(error, __CREATE_PROGRAM_WITH_BINARY_ERR);
        if (err != NULL) {
            *err = error;
        }
    }

    
#if defined(CL_VERSION_1_2)
    /**
     * Create program using builtin kernels.
     * \param kernelNames Semi-colon separated list of builtin kernel names
     */
    Program(
        const Context& context,
        const VECTOR_CLASS<Device>& devices,
        const STRING_CLASS& kernelNames,
        cl_int* err = NULL)
    {
        cl_int error;
        
        object_ = ::clCreateProgramWithBuiltInKernels(
            context(), 
            (cl_uint) devices.size(),
            (cl_device_id*)&devices.front(),
            kernelNames.c_str(), 
            &error);

        detail::errHandler(error, __CREATE_PROGRAM_WITH_BUILT_IN_KERNELS_ERR);
        if (err != NULL) {
            *err = error;
        }
    }
#endif // #if defined(CL_VERSION_1_2)

    Program() { }

    Program(const Program& program) : detail::Wrapper<cl_type>(program) { }

    Program(const cl_program& program) : detail::Wrapper<cl_type>(program) { }

    Program& operator = (const Program& rhs)
    {
        if (this != &rhs) {
            detail::Wrapper<cl_type>::operator=(rhs);
        }
        return *this;
    }

    Program& operator = (const cl_program& rhs)
    {
        detail::Wrapper<cl_type>::operator=(rhs);
        return *this;
    }

    cl_int build(
        const VECTOR_CLASS<Device>& devices,
        const char* options = NULL,
        void (CL_CALLBACK * notifyFptr)(cl_program, void *) = NULL,
        void* data = NULL) const
    {
        return detail::errHandler(
            ::clBuildProgram(
                object_,
                (cl_uint)
                devices.size(),
                (cl_device_id*)&devices.front(),
                options,
                notifyFptr,
                data),
                __BUILD_PROGRAM_ERR);
    }

    cl_int build(
        const char* options = NULL,
        void (CL_CALLBACK * notifyFptr)(cl_program, void *) = NULL,
        void* data = NULL) const
    {
        return detail::errHandler(
            ::clBuildProgram(
                object_,
                0,
                NULL,
                options,
                notifyFptr,
                data),
                __BUILD_PROGRAM_ERR);
    }

#if defined(CL_VERSION_1_2)
	cl_int compile(
        const char* options = NULL,
        void (CL_CALLBACK * notifyFptr)(cl_program, void *) = NULL,
        void* data = NULL) const
    {
        return detail::errHandler(
            ::clCompileProgram(
                object_,
                0,
                NULL,
                options,
				0,
				NULL,
				NULL,
                notifyFptr,
                data),
                __COMPILE_PROGRAM_ERR);
    }
#endif

    template <typename T>
    cl_int getInfo(cl_program_info name, T* param) const
    {
        return detail::errHandler(
            detail::getInfo(&::clGetProgramInfo, object_, name, param),
            __GET_PROGRAM_INFO_ERR);
    }

    template <cl_int name> typename
    detail::param_traits<detail::cl_program_info, name>::param_type
    getInfo(cl_int* err = NULL) const
    {
        typename detail::param_traits<
            detail::cl_program_info, name>::param_type param;
        cl_int result = getInfo(name, &param);
        if (err != NULL) {
            *err = result;
        }
        return param;
    }

    template <typename T>
    cl_int getBuildInfo(
        const Device& device, cl_program_build_info name, T* param) const
    {
        return detail::errHandler(
            detail::getInfo(
                &::clGetProgramBuildInfo, object_, device(), name, param),
                __GET_PROGRAM_BUILD_INFO_ERR);
    }

    template <cl_int name> typename
    detail::param_traits<detail::cl_program_build_info, name>::param_type
    getBuildInfo(const Device& device, cl_int* err = NULL) const
    {
        typename detail::param_traits<
            detail::cl_program_build_info, name>::param_type param;
        cl_int result = getBuildInfo(device, name, &param);
        if (err != NULL) {
            *err = result;
        }
        return param;
    }

    cl_int createKernels(VECTOR_CLASS<Kernel>* kernels)
    {
        cl_uint numKernels;
        cl_int err = ::clCreateKernelsInProgram(object_, 0, NULL, &numKernels);
        if (err != CL_SUCCESS) {
            return detail::errHandler(err, __CREATE_KERNELS_IN_PROGRAM_ERR);
        }

        Kernel* value = (Kernel*) alloca(numKernels * sizeof(Kernel));
        err = ::clCreateKernelsInProgram(
            object_, numKernels, (cl_kernel*) value, NULL);
        if (err != CL_SUCCESS) {
            return detail::errHandler(err, __CREATE_KERNELS_IN_PROGRAM_ERR);
        }

        kernels->assign(&value[0], &value[numKernels]);
        return CL_SUCCESS;
    }
};

#if defined(CL_VERSION_1_2)
inline Program linkProgram(
    Program input1,
    Program input2,
    const char* options = NULL,
    void (CL_CALLBACK * notifyFptr)(cl_program, void *) = NULL,
    void* data = NULL,
    cl_int* err = NULL) 
{
    cl_int err_local;

    cl_program programs[2] = { input1(), input2() };

    cl_program prog = ::clLinkProgram(
        Context::getDefault()(),
        0,
        NULL,
        options,
        2,
        programs,
        notifyFptr,
        data,
        &err_local);

    if (err_local != CL_SUCCESS) {
        err_local = detail::errHandler(err_local,__COMPILE_PROGRAM_ERR);
        if (err != NULL) {
            *err = err_local;
        }
    }

    return Program(prog);
}

inline Program linkProgram(
    VECTOR_CLASS<Program> inputPrograms,
    const char* options = NULL,
    void (CL_CALLBACK * notifyFptr)(cl_program, void *) = NULL,
    void* data = NULL,
    cl_int* err = NULL) 
{
    cl_int err_local;

    if (inputPrograms.size() == 0) {
        err_local = detail::errHandler(CL_INVALID_VALUE,__COMPILE_PROGRAM_ERR);
        if (err != NULL) {
            *err = err_local;
        }
        return Program();
    }

    cl_program * programs = (cl_program*) alloca(inputPrograms.size() * sizeof(cl_program));

    if (programs != NULL) {
        for (unsigned int i = 0; i < inputPrograms.size(); i++) {
          programs[i] = inputPrograms[i]();
        }
    } 
    else {
        err_local = detail::errHandler(CL_OUT_OF_HOST_MEMORY,__COMPILE_PROGRAM_ERR);
        if (err != NULL) {
            *err = err_local;
        }
        return Program();
    }

    cl_program prog = ::clLinkProgram(
        Context::getDefault()(),
        0,
        NULL,
        options,
        (cl_uint)inputPrograms.size(),
        programs,
        notifyFptr,
        data,
        &err_local);

    if (err_local != CL_SUCCESS) {
        err_local = detail::errHandler(err_local,__COMPILE_PROGRAM_ERR);
        if (err != NULL) {
            *err = err_local;
        }
    }

    return Program(prog);
}
#endif

template<>
inline VECTOR_CLASS<char *> cl::Program::getInfo<CL_PROGRAM_BINARIES>(cl_int* err) const
{
    VECTOR_CLASS< ::size_t> sizes = getInfo<CL_PROGRAM_BINARY_SIZES>();
    VECTOR_CLASS<char *> binaries;
    for (VECTOR_CLASS< ::size_t>::iterator s = sizes.begin(); s != sizes.end(); ++s) 
    {
        char *ptr = NULL;
        if (*s != 0) 
            ptr = new char[*s];
        binaries.push_back(ptr);
    }
    
    cl_int result = getInfo(CL_PROGRAM_BINARIES, &binaries);
    if (err != NULL) {
        *err = result;
    }
    return binaries;
}

__GET_INFO_HELPER_WITH_RETAIN(cl::Program)

inline Kernel::Kernel(const Program& program, const char* name, cl_int* err)
{
    cl_int error;

    object_ = ::clCreateKernel(program(), name, &error);
    detail::errHandler(error, __CREATE_KERNEL_ERR);

    if (err != NULL) {
        *err = error;
    }

}

/*! \class CommandQueue
 * \brief CommandQueue interface for cl_command_queue.
 */
class CommandQueue : public detail::Wrapper<cl_command_queue>
{
private:
    static volatile int default_initialized_;
    static CommandQueue default_;
    static volatile cl_int default_error_;
public:
   CommandQueue(
        cl_command_queue_properties properties,
        cl_int* err = NULL)
    {
        cl_int error;

        Context context = Context::getDefault(&error);
        detail::errHandler(error, __CREATE_COMMAND_QUEUE_ERR);

        if (error != CL_SUCCESS) {
            if (err != NULL) {
                *err = error;
            }
        }
        else {
            Device device = context.getInfo<CL_CONTEXT_DEVICES>()[0];

            object_ = ::clCreateCommandQueue(
                context(), device(), properties, &error);

            detail::errHandler(error, __CREATE_COMMAND_QUEUE_ERR);
            if (err != NULL) {
                *err = error;
            }
        }
    }

    CommandQueue(
        const Context& context,
        const Device& device,
        cl_command_queue_properties properties = 0,
        cl_int* err = NULL)
    {
        cl_int error;
        object_ = ::clCreateCommandQueue(
            context(), device(), properties, &error);

        detail::errHandler(error, __CREATE_COMMAND_QUEUE_ERR);
        if (err != NULL) {
            *err = error;
        }
    }

    static CommandQueue getDefault(cl_int * err = NULL) 
    {
        int state = detail::compare_exchange(
            &default_initialized_, 
            __DEFAULT_BEING_INITIALIZED, __DEFAULT_NOT_INITIALIZED);
        
        if (state & __DEFAULT_INITIALIZED) {
            if (err != NULL) {
                *err = default_error_;
            }
            return default_;
        }

        if (state & __DEFAULT_BEING_INITIALIZED) {
              // Assume writes will propagate eventually...
              while(default_initialized_ != __DEFAULT_INITIALIZED) {
                                ;
              }

            if (err != NULL) {
                *err = default_error_;
            }
            return default_;
        }

        cl_int error;

        Context context = Context::getDefault(&error);
        detail::errHandler(error, __CREATE_COMMAND_QUEUE_ERR);

        if (error != CL_SUCCESS) {
            if (err != NULL) {
                *err = error;
            }
        }
        else {
            Device device = context.getInfo<CL_CONTEXT_DEVICES>()[0];

            default_ = CommandQueue(context, device, 0, &error);

            detail::errHandler(error, __CREATE_COMMAND_QUEUE_ERR);
            if (err != NULL) {
                *err = error;
            }
        }

        default_error_ = error;
        // Assume writes will propagate eventually...
        default_initialized_ = __DEFAULT_INITIALIZED;
        if (err != NULL) {
            *err = default_error_;
        }
        return default_;

    }

    CommandQueue() { }

    CommandQueue(const CommandQueue& commandQueue) : detail::Wrapper<cl_type>(commandQueue) { }

    CommandQueue(const cl_command_queue& commandQueue) : detail::Wrapper<cl_type>(commandQueue) { }

    CommandQueue& operator = (const CommandQueue& rhs)
    {
        if (this != &rhs) {
            detail::Wrapper<cl_type>::operator=(rhs);
        }
        return *this;
    }

    CommandQueue& operator = (const cl_command_queue& rhs)
    {
        detail::Wrapper<cl_type>::operator=(rhs);
        return *this;
    }

    template <typename T>
    cl_int getInfo(cl_command_queue_info name, T* param) const
    {
        return detail::errHandler(
            detail::getInfo(
                &::clGetCommandQueueInfo, object_, name, param),
                __GET_COMMAND_QUEUE_INFO_ERR);
    }

    template <cl_int name> typename
    detail::param_traits<detail::cl_command_queue_info, name>::param_type
    getInfo(cl_int* err = NULL) const
    {
        typename detail::param_traits<
            detail::cl_command_queue_info, name>::param_type param;
        cl_int result = getInfo(name, &param);
        if (err != NULL) {
            *err = result;
        }
        return param;
    }

    cl_int enqueueReadBuffer(
        const Buffer& buffer,
        cl_bool blocking,
        ::size_t offset,
        ::size_t size,
        void* ptr,
        const VECTOR_CLASS<Event>* events = NULL,
        Event* event = NULL) const
    {
        cl_event tmp;
        cl_int err = detail::errHandler(
            ::clEnqueueReadBuffer(
                object_, buffer(), blocking, offset, size,
                ptr,
                (events != NULL) ? (cl_uint) events->size() : 0,
                (events != NULL && events->size() > 0) ? (cl_event*) &events->front() : NULL,
                (event != NULL) ? &tmp : NULL),
            __ENQUEUE_READ_BUFFER_ERR);

        if (event != NULL && err == CL_SUCCESS)
            *event = tmp;

        return err;
    }

    cl_int enqueueWriteBuffer(
        const Buffer& buffer,
        cl_bool blocking,
        ::size_t offset,
        ::size_t size,
        const void* ptr,
        const VECTOR_CLASS<Event>* events = NULL,
        Event* event = NULL) const
    {
        cl_event tmp;
        cl_int err = detail::errHandler(
            ::clEnqueueWriteBuffer(
                object_, buffer(), blocking, offset, size,
                ptr,
                (events != NULL) ? (cl_uint) events->size() : 0,
                (events != NULL && events->size() > 0) ? (cl_event*) &events->front() : NULL,
                (event != NULL) ? &tmp : NULL),
                __ENQUEUE_WRITE_BUFFER_ERR);

        if (event != NULL && err == CL_SUCCESS)
            *event = tmp;

        return err;
    }

    cl_int enqueueCopyBuffer(
        const Buffer& src,
        const Buffer& dst,
        ::size_t src_offset,
        ::size_t dst_offset,
        ::size_t size,
        const VECTOR_CLASS<Event>* events = NULL,
        Event* event = NULL) const
    {
        cl_event tmp;
        cl_int err = detail::errHandler(
            ::clEnqueueCopyBuffer(
                object_, src(), dst(), src_offset, dst_offset, size,
                (events != NULL) ? (cl_uint) events->size() : 0,
                (events != NULL && events->size() > 0) ? (cl_event*) &events->front() : NULL,
                (event != NULL) ? &tmp : NULL),
            __ENQEUE_COPY_BUFFER_ERR);

        if (event != NULL && err == CL_SUCCESS)
            *event = tmp;

        return err;
    }

    cl_int enqueueReadBufferRect(
        const Buffer& buffer,
        cl_bool blocking,
        const size_t<3>& buffer_offset,
        const size_t<3>& host_offset,
        const size_t<3>& region,
        ::size_t buffer_row_pitch,
        ::size_t buffer_slice_pitch,
        ::size_t host_row_pitch,
        ::size_t host_slice_pitch,
        void *ptr,
        const VECTOR_CLASS<Event>* events = NULL,
        Event* event = NULL) const
    {
        cl_event tmp;
        cl_int err = detail::errHandler(
            ::clEnqueueReadBufferRect(
                object_, 
                buffer(), 
                blocking, 
                (const ::size_t *)buffer_offset,
                (const ::size_t *)host_offset,
                (const ::size_t *)region,
                buffer_row_pitch,
                buffer_slice_pitch,
                host_row_pitch,
                host_slice_pitch,
                ptr,
                (events != NULL) ? (cl_uint) events->size() : 0,
                (events != NULL && events->size() > 0) ? (cl_event*) &events->front() : NULL,
                (event != NULL) ? &tmp : NULL),
                __ENQUEUE_READ_BUFFER_RECT_ERR);

        if (event != NULL && err == CL_SUCCESS)
            *event = tmp;

        return err;
    }

    cl_int enqueueWriteBufferRect(
        const Buffer& buffer,
        cl_bool blocking,
        const size_t<3>& buffer_offset,
        const size_t<3>& host_offset,
        const size_t<3>& region,
        ::size_t buffer_row_pitch,
        ::size_t buffer_slice_pitch,
        ::size_t host_row_pitch,
        ::size_t host_slice_pitch,
        void *ptr,
        const VECTOR_CLASS<Event>* events = NULL,
        Event* event = NULL) const
    {
        cl_event tmp;
        cl_int err = detail::errHandler(
            ::clEnqueueWriteBufferRect(
                object_, 
                buffer(), 
                blocking, 
                (const ::size_t *)buffer_offset,
                (const ::size_t *)host_offset,
                (const ::size_t *)region,
                buffer_row_pitch,
                buffer_slice_pitch,
                host_row_pitch,
                host_slice_pitch,
                ptr,
                (events != NULL) ? (cl_uint) events->size() : 0,
                (events != NULL && events->size() > 0) ? (cl_event*) &events->front() : NULL,
                (event != NULL) ? &tmp : NULL),
                __ENQUEUE_WRITE_BUFFER_RECT_ERR);

        if (event != NULL && err == CL_SUCCESS)
            *event = tmp;

        return err;
    }

    cl_int enqueueCopyBufferRect(
        const Buffer& src,
        const Buffer& dst,
        const size_t<3>& src_origin,
        const size_t<3>& dst_origin,
        const size_t<3>& region,
        ::size_t src_row_pitch,
        ::size_t src_slice_pitch,
        ::size_t dst_row_pitch,
        ::size_t dst_slice_pitch,
        const VECTOR_CLASS<Event>* events = NULL,
        Event* event = NULL) const
    {
        cl_event tmp;
        cl_int err = detail::errHandler(
            ::clEnqueueCopyBufferRect(
                object_, 
                src(), 
                dst(), 
                (const ::size_t *)src_origin, 
                (const ::size_t *)dst_origin, 
                (const ::size_t *)region,
                src_row_pitch,
                src_slice_pitch,
                dst_row_pitch,
                dst_slice_pitch,
                (events != NULL) ? (cl_uint) events->size() : 0,
                (events != NULL && events->size() > 0) ? (cl_event*) &events->front() : NULL,
                (event != NULL) ? &tmp : NULL),
            __ENQEUE_COPY_BUFFER_RECT_ERR);

        if (event != NULL && err == CL_SUCCESS)
            *event = tmp;

        return err;
    }

#if defined(CL_VERSION_1_2)
    /**
     * Enqueue a command to fill a buffer object with a pattern
     * of a given size. The pattern is specified a as vector.
     * \tparam PatternType The datatype of the pattern field. 
     *     The pattern type must be an accepted OpenCL data type.
     */
    template<typename PatternType>
    cl_int enqueueFillBuffer(
        const Buffer& buffer,
        PatternType pattern,
        ::size_t offset,
        ::size_t size,
        const VECTOR_CLASS<Event>* events = NULL,
        Event* event = NULL) const
    {
        cl_event tmp;
        cl_int err = detail::errHandler(
            ::clEnqueueFillBuffer(
                object_, 
                buffer(),
                static_cast<void*>(&pattern),
                sizeof(PatternType), 
                offset, 
                size,
                (events != NULL) ? (cl_uint) events->size() : 0,
                (events != NULL && events->size() > 0) ? (cl_event*) &events->front() : NULL,
                (event != NULL) ? &tmp : NULL),
                __ENQUEUE_FILL_BUFFER_ERR);

        if (event != NULL && err == CL_SUCCESS)
            *event = tmp;

        return err;
    }
#endif // #if defined(CL_VERSION_1_2)

    cl_int enqueueReadImage(
        const Image& image,
        cl_bool blocking,
        const size_t<3>& origin,
        const size_t<3>& region,
        ::size_t row_pitch,
        ::size_t slice_pitch,
        void* ptr,
        const VECTOR_CLASS<Event>* events = NULL,
        Event* event = NULL) const
    {
        cl_event tmp;
        cl_int err = detail::errHandler(
            ::clEnqueueReadImage(
                object_, image(), blocking, (const ::size_t *) origin,
                (const ::size_t *) region, row_pitch, slice_pitch, ptr,
                (events != NULL) ? (cl_uint) events->size() : 0,
                (events != NULL && events->size() > 0) ? (cl_event*) &events->front() : NULL,
                (event != NULL) ? &tmp : NULL),
            __ENQUEUE_READ_IMAGE_ERR);

        if (event != NULL && err == CL_SUCCESS)
            *event = tmp;

        return err;
    }

    cl_int enqueueWriteImage(
        const Image& image,
        cl_bool blocking,
        const size_t<3>& origin,
        const size_t<3>& region,
        ::size_t row_pitch,
        ::size_t slice_pitch,
        void* ptr,
        const VECTOR_CLASS<Event>* events = NULL,
        Event* event = NULL) const
    {
        cl_event tmp;
        cl_int err = detail::errHandler(
            ::clEnqueueWriteImage(
                object_, image(), blocking, (const ::size_t *) origin,
                (const ::size_t *) region, row_pitch, slice_pitch, ptr,
                (events != NULL) ? (cl_uint) events->size() : 0,
                (events != NULL && events->size() > 0) ? (cl_event*) &events->front() : NULL,
                (event != NULL) ? &tmp : NULL),
            __ENQUEUE_WRITE_IMAGE_ERR);

        if (event != NULL && err == CL_SUCCESS)
            *event = tmp;

        return err;
    }

    cl_int enqueueCopyImage(
        const Image& src,
        const Image& dst,
        const size_t<3>& src_origin,
        const size_t<3>& dst_origin,
        const size_t<3>& region,
        const VECTOR_CLASS<Event>* events = NULL,
        Event* event = NULL) const
    {
        cl_event tmp;
        cl_int err = detail::errHandler(
            ::clEnqueueCopyImage(
                object_, src(), dst(), (const ::size_t *) src_origin,
                (const ::size_t *)dst_origin, (const ::size_t *) region,
                (events != NULL) ? (cl_uint) events->size() : 0,
                (events != NULL && events->size() > 0) ? (cl_event*) &events->front() : NULL,
                (event != NULL) ? &tmp : NULL),
            __ENQUEUE_COPY_IMAGE_ERR);

        if (event != NULL && err == CL_SUCCESS)
            *event = tmp;

        return err;
    }

#if defined(CL_VERSION_1_2)
    /**
     * Enqueue a command to fill an image object with a specified color.
     * \param fillColor is the color to use to fill the image.
     *     This is a four component RGBA floating-point color value if
     *     the image channel data type is not an unnormalized signed or
     *     unsigned data type.
     */
    cl_int enqueueFillImage(
        const Image& image,
        cl_float4 fillColor,
        const size_t<3>& origin,
        const size_t<3>& region,
        const VECTOR_CLASS<Event>* events = NULL,
        Event* event = NULL) const
    {
        cl_event tmp;
        cl_int err = detail::errHandler(
            ::clEnqueueFillImage(
                object_, 
                image(),
                static_cast<void*>(&fillColor), 
                (const ::size_t *) origin, 
                (const ::size_t *) region,
                (events != NULL) ? (cl_uint) events->size() : 0,
                (events != NULL && events->size() > 0) ? (cl_event*) &events->front() : NULL,
                (event != NULL) ? &tmp : NULL),
                __ENQUEUE_FILL_IMAGE_ERR);

        if (event != NULL && err == CL_SUCCESS)
            *event = tmp;

        return err;
    }

    /**
     * Enqueue a command to fill an image object with a specified color.
     * \param fillColor is the color to use to fill the image.
     *     This is a four component RGBA signed integer color value if
     *     the image channel data type is an unnormalized signed integer
     *     type.
     */
    cl_int enqueueFillImage(
        const Image& image,
        cl_int4 fillColor,
        const size_t<3>& origin,
        const size_t<3>& region,
        const VECTOR_CLASS<Event>* events = NULL,
        Event* event = NULL) const
    {
        cl_event tmp;
        cl_int err = detail::errHandler(
            ::clEnqueueFillImage(
                object_, 
                image(),
                static_cast<void*>(&fillColor), 
                (const ::size_t *) origin, 
                (const ::size_t *) region,
                (events != NULL) ? (cl_uint) events->size() : 0,
                (events != NULL && events->size() > 0) ? (cl_event*) &events->front() : NULL,
                (event != NULL) ? &tmp : NULL),
                __ENQUEUE_FILL_IMAGE_ERR);

        if (event != NULL && err == CL_SUCCESS)
            *event = tmp;

        return err;
    }

    /**
     * Enqueue a command to fill an image object with a specified color.
     * \param fillColor is the color to use to fill the image.
     *     This is a four component RGBA unsigned integer color value if
     *     the image channel data type is an unnormalized unsigned integer
     *     type.
     */
    cl_int enqueueFillImage(
        const Image& image,
        cl_uint4 fillColor,
        const size_t<3>& origin,
        const size_t<3>& region,
        const VECTOR_CLASS<Event>* events = NULL,
        Event* event = NULL) const
    {
        cl_event tmp;
        cl_int err = detail::errHandler(
            ::clEnqueueFillImage(
                object_, 
                image(),
                static_cast<void*>(&fillColor), 
                (const ::size_t *) origin, 
                (const ::size_t *) region,
                (events != NULL) ? (cl_uint) events->size() : 0,
                (events != NULL && events->size() > 0) ? (cl_event*) &events->front() : NULL,
                (event != NULL) ? &tmp : NULL),
                __ENQUEUE_FILL_IMAGE_ERR);

        if (event != NULL && err == CL_SUCCESS)
            *event = tmp;

        return err;
    }
#endif // #if defined(CL_VERSION_1_2)

    cl_int enqueueCopyImageToBuffer(
        const Image& src,
        const Buffer& dst,
        const size_t<3>& src_origin,
        const size_t<3>& region,
        ::size_t dst_offset,
        const VECTOR_CLASS<Event>* events = NULL,
        Event* event = NULL) const
    {
        cl_event tmp;
        cl_int err = detail::errHandler(
            ::clEnqueueCopyImageToBuffer(
                object_, src(), dst(), (const ::size_t *) src_origin,
                (const ::size_t *) region, dst_offset,
                (events != NULL) ? (cl_uint) events->size() : 0,
                (events != NULL && events->size() > 0) ? (cl_event*) &events->front() : NULL,
                (event != NULL) ? &tmp : NULL),
            __ENQUEUE_COPY_IMAGE_TO_BUFFER_ERR);

        if (event != NULL && err == CL_SUCCESS)
            *event = tmp;

        return err;
    }

    cl_int enqueueCopyBufferToImage(
        const Buffer& src,
        const Image& dst,
        ::size_t src_offset,
        const size_t<3>& dst_origin,
        const size_t<3>& region,
        const VECTOR_CLASS<Event>* events = NULL,
        Event* event = NULL) const
    {
        cl_event tmp;
        cl_int err = detail::errHandler(
            ::clEnqueueCopyBufferToImage(
                object_, src(), dst(), src_offset,
                (const ::size_t *) dst_origin, (const ::size_t *) region,
                (events != NULL) ? (cl_uint) events->size() : 0,
                (events != NULL && events->size() > 0) ? (cl_event*) &events->front() : NULL,
                (event != NULL) ? &tmp : NULL),
            __ENQUEUE_COPY_BUFFER_TO_IMAGE_ERR);

        if (event != NULL && err == CL_SUCCESS)
            *event = tmp;

        return err;
    }

    void* enqueueMapBuffer(
        const Buffer& buffer,
        cl_bool blocking,
        cl_map_flags flags,
        ::size_t offset,
        ::size_t size,
        const VECTOR_CLASS<Event>* events = NULL,
        Event* event = NULL,
        cl_int* err = NULL) const
    {
        cl_int error;
        void * result = ::clEnqueueMapBuffer(
            object_, buffer(), blocking, flags, offset, size,
            (events != NULL) ? (cl_uint) events->size() : 0,
            (events != NULL && events->size() > 0) ? (cl_event*) &events->front() : NULL,
            (cl_event*) event,
            &error);

        detail::errHandler(error, __ENQUEUE_MAP_BUFFER_ERR);
        if (err != NULL) {
            *err = error;
        }
        return result;
    }

    void* enqueueMapImage(
        const Image& buffer,
        cl_bool blocking,
        cl_map_flags flags,
        const size_t<3>& origin,
        const size_t<3>& region,
        ::size_t * row_pitch,
        ::size_t * slice_pitch,
        const VECTOR_CLASS<Event>* events = NULL,
        Event* event = NULL,
        cl_int* err = NULL) const
    {
        cl_int error;
        void * result = ::clEnqueueMapImage(
            object_, buffer(), blocking, flags,
            (const ::size_t *) origin, (const ::size_t *) region,
            row_pitch, slice_pitch,
            (events != NULL) ? (cl_uint) events->size() : 0,
            (events != NULL && events->size() > 0) ? (cl_event*) &events->front() : NULL,
            (cl_event*) event,
            &error);

        detail::errHandler(error, __ENQUEUE_MAP_IMAGE_ERR);
        if (err != NULL) {
              *err = error;
        }
        return result;
    }

    cl_int enqueueUnmapMemObject(
        const Memory& memory,
        void* mapped_ptr,
        const VECTOR_CLASS<Event>* events = NULL,
        Event* event = NULL) const
    {
        cl_event tmp;
        cl_int err = detail::errHandler(
            ::clEnqueueUnmapMemObject(
                object_, memory(), mapped_ptr,
                (events != NULL) ? (cl_uint) events->size() : 0,
                (events != NULL && events->size() > 0) ? (cl_event*) &events->front() : NULL,
                (event != NULL) ? &tmp : NULL),
            __ENQUEUE_UNMAP_MEM_OBJECT_ERR);

        if (event != NULL && err == CL_SUCCESS)
            *event = tmp;

        return err;
    }

#if defined(CL_VERSION_1_2)
    /**
     * Enqueues a marker command which waits for either a list of events to complete, 
     * or all previously enqueued commands to complete.
     *
     * Enqueues a marker command which waits for either a list of events to complete, 
     * or if the list is empty it waits for all commands previously enqueued in command_queue 
     * to complete before it completes. This command returns an event which can be waited on, 
     * i.e. this event can be waited on to insure that all events either in the event_wait_list 
     * or all previously enqueued commands, queued before this command to command_queue, 
     * have completed.
     */
    cl_int enqueueMarkerWithWaitList(
        const VECTOR_CLASS<Event> *events = 0,
        Event *event = 0)
    {
        cl_event tmp;
        cl_int err = detail::errHandler(
            ::clEnqueueMarkerWithWaitList(
                object_,
                (events != NULL) ? (cl_uint) events->size() : 0,
                (events != NULL && events->size() > 0) ? (cl_event*) &events->front() : NULL,
                (event != NULL) ? &tmp : NULL),
            __ENQUEUE_UNMAP_MEM_OBJECT_ERR);

        if (event != NULL && err == CL_SUCCESS)
            *event = tmp;

        return err;
    }

    /**
     * A synchronization point that enqueues a barrier operation.
     *
     * Enqueues a barrier command which waits for either a list of events to complete, 
     * or if the list is empty it waits for all commands previously enqueued in command_queue 
     * to complete before it completes. This command blocks command execution, that is, any 
     * following commands enqueued after it do not execute until it completes. This command 
     * returns an event which can be waited on, i.e. this event can be waited on to insure that 
     * all events either in the event_wait_list or all previously enqueued commands, queued 
     * before this command to command_queue, have completed.
     */
    cl_int clEnqueueBarrierWithWaitList(
        const VECTOR_CLASS<Event> *events = 0,
        Event *event = 0)
    {
        cl_event tmp;
        cl_int err = detail::errHandler(
            ::clEnqueueMarkerWithWaitList(
                object_,
                (events != NULL) ? (cl_uint) events->size() : 0,
                (events != NULL && events->size() > 0) ? (cl_event*) &events->front() : NULL,
                (event != NULL) ? &tmp : NULL),
            __ENQUEUE_UNMAP_MEM_OBJECT_ERR);

        if (event != NULL && err == CL_SUCCESS)
            *event = tmp;

        return err;
    }
    
    /**
     * Enqueues a command to indicate with which device a set of memory objects
     * should be associated.
     */
    cl_int enqueueMigrateMemObjects(
        const VECTOR_CLASS<Memory> &memObjects,
        cl_mem_migration_flags flags,
        const VECTOR_CLASS<Event>* events = NULL,
        Event* event = NULL
        )
    {
        cl_event tmp;
        
        cl_mem* localMemObjects = static_cast<cl_mem*>(alloca(memObjects.size() * sizeof(cl_mem)));
        for( int i = 0; i < (int)memObjects.size(); ++i ) {
            localMemObjects[i] = memObjects[i]();
        }


        cl_int err = detail::errHandler(
            ::clEnqueueMigrateMemObjects(
                object_, 
                (cl_uint)memObjects.size(), 
                static_cast<const cl_mem*>(localMemObjects),
                flags,
                (events != NULL) ? (cl_uint) events->size() : 0,
                (events != NULL && events->size() > 0) ? (cl_event*) &events->front() : NULL,
                (event != NULL) ? &tmp : NULL),
            __ENQUEUE_UNMAP_MEM_OBJECT_ERR);

        if (event != NULL && err == CL_SUCCESS)
            *event = tmp;

        return err;
    }
#endif // #if defined(CL_VERSION_1_2)

    cl_int enqueueNDRangeKernel(
        const Kernel& kernel,
        const NDRange& offset,
        const NDRange& global,
        const NDRange& local = NullRange,
        const VECTOR_CLASS<Event>* events = NULL,
        Event* event = NULL) const
    {
        cl_event tmp;
        cl_int err = detail::errHandler(
            ::clEnqueueNDRangeKernel(
                object_, kernel(), (cl_uint) global.dimensions(),
                offset.dimensions() != 0 ? (const ::size_t*) offset : NULL,
                (const ::size_t*) global,
                local.dimensions() != 0 ? (const ::size_t*) local : NULL,
                (events != NULL) ? (cl_uint) events->size() : 0,
                (events != NULL && events->size() > 0) ? (cl_event*) &events->front() : NULL,
                (event != NULL) ? &tmp : NULL),
            __ENQUEUE_NDRANGE_KERNEL_ERR);

        if (event != NULL && err == CL_SUCCESS)
            *event = tmp;

        return err;
    }

    cl_int enqueueTask(
        const Kernel& kernel,
        const VECTOR_CLASS<Event>* events = NULL,
        Event* event = NULL) const
    {
        cl_event tmp;
        cl_int err = detail::errHandler(
            ::clEnqueueTask(
                object_, kernel(),
                (events != NULL) ? (cl_uint) events->size() : 0,
                (events != NULL && events->size() > 0) ? (cl_event*) &events->front() : NULL,
                (event != NULL) ? &tmp : NULL),
            __ENQUEUE_TASK_ERR);

        if (event != NULL && err == CL_SUCCESS)
            *event = tmp;

        return err;
    }

    cl_int enqueueNativeKernel(
        void (CL_CALLBACK *userFptr)(void *),
        std::pair<void*, ::size_t> args,
        const VECTOR_CLASS<Memory>* mem_objects = NULL,
        const VECTOR_CLASS<const void*>* mem_locs = NULL,
        const VECTOR_CLASS<Event>* events = NULL,
        Event* event = NULL) const
    {
        cl_mem * mems = (mem_objects != NULL && mem_objects->size() > 0) 
            ? (cl_mem*) alloca(mem_objects->size() * sizeof(cl_mem))
            : NULL;

        if (mems != NULL) {
            for (unsigned int i = 0; i < mem_objects->size(); i++) {
                mems[i] = ((*mem_objects)[i])();
            }
        }

        cl_event tmp;
        cl_int err = detail::errHandler(
            ::clEnqueueNativeKernel(
                object_, userFptr, args.first, args.second,
                (mem_objects != NULL) ? (cl_uint) mem_objects->size() : 0,
                mems,
                (mem_locs != NULL) ? (const void **) &mem_locs->front() : NULL,
                (events != NULL) ? (cl_uint) events->size() : 0,
                (events != NULL && events->size() > 0) ? (cl_event*) &events->front() : NULL,
                (event != NULL) ? &tmp : NULL),
            __ENQUEUE_NATIVE_KERNEL);

        if (event != NULL && err == CL_SUCCESS)
            *event = tmp;

        return err;
    }

/**
 * Deprecated APIs for 1.2
 */
#if defined(CL_VERSION_1_1)
    cl_int enqueueMarker(Event* event = NULL) const
    {
        return detail::errHandler(
            ::clEnqueueMarker(object_, (cl_event*) event),
            __ENQUEUE_MARKER_ERR);
    }

    cl_int enqueueWaitForEvents(const VECTOR_CLASS<Event>& events) const
    {
        return detail::errHandler(
            ::clEnqueueWaitForEvents(
                object_,
                (cl_uint) events.size(),
                (const cl_event*) &events.front()),
            __ENQUEUE_WAIT_FOR_EVENTS_ERR);
    }
#endif // #if defined(CL_VERSION_1_1)

    cl_int enqueueAcquireGLObjects(
         const VECTOR_CLASS<Memory>* mem_objects = NULL,
         const VECTOR_CLASS<Event>* events = NULL,
         Event* event = NULL) const
     {
        cl_event tmp;
        cl_int err = detail::errHandler(
             ::clEnqueueAcquireGLObjects(
                 object_,
                 (mem_objects != NULL) ? (cl_uint) mem_objects->size() : 0,
                 (mem_objects != NULL) ? (const cl_mem *) &mem_objects->front(): NULL,
                 (events != NULL) ? (cl_uint) events->size() : 0,
                 (events != NULL && events->size() > 0) ? (cl_event*) &events->front() : NULL,
                 (event != NULL) ? &tmp : NULL),
             __ENQUEUE_ACQUIRE_GL_ERR);

        if (event != NULL && err == CL_SUCCESS)
            *event = tmp;

        return err;
     }

    cl_int enqueueReleaseGLObjects(
         const VECTOR_CLASS<Memory>* mem_objects = NULL,
         const VECTOR_CLASS<Event>* events = NULL,
         Event* event = NULL) const
     {
        cl_event tmp;
        cl_int err = detail::errHandler(
             ::clEnqueueReleaseGLObjects(
                 object_,
                 (mem_objects != NULL) ? (cl_uint) mem_objects->size() : 0,
                 (mem_objects != NULL) ? (const cl_mem *) &mem_objects->front(): NULL,
                 (events != NULL) ? (cl_uint) events->size() : 0,
                 (events != NULL && events->size() > 0) ? (cl_event*) &events->front() : NULL,
                 (event != NULL) ? &tmp : NULL),
             __ENQUEUE_RELEASE_GL_ERR);

        if (event != NULL && err == CL_SUCCESS)
            *event = tmp;

        return err;
     }

#if defined (USE_DX_INTEROP)
typedef CL_API_ENTRY cl_int (CL_API_CALL *PFN_clEnqueueAcquireD3D10ObjectsKHR)(
    cl_command_queue command_queue, cl_uint num_objects,
    const cl_mem* mem_objects, cl_uint num_events_in_wait_list,
    const cl_event* event_wait_list, cl_event* event);
typedef CL_API_ENTRY cl_int (CL_API_CALL *PFN_clEnqueueReleaseD3D10ObjectsKHR)(
    cl_command_queue command_queue, cl_uint num_objects,
    const cl_mem* mem_objects,  cl_uint num_events_in_wait_list,
    const cl_event* event_wait_list, cl_event* event);

    cl_int enqueueAcquireD3D10Objects(
         const VECTOR_CLASS<Memory>* mem_objects = NULL,
         const VECTOR_CLASS<Event>* events = NULL,
         Event* event = NULL) const
    {
        static PFN_clEnqueueAcquireD3D10ObjectsKHR pfn_clEnqueueAcquireD3D10ObjectsKHR = NULL;
#if defined(CL_VERSION_1_2)
        cl_context context = getInfo<CL_QUEUE_CONTEXT>();
        cl::Device device(getInfo<CL_QUEUE_DEVICE>());
        cl_platform_id platform = device.getInfo<CL_DEVICE_PLATFORM>();
        __INIT_CL_EXT_FCN_PTR_PLATFORM(platform, clEnqueueAcquireD3D10ObjectsKHR);
#endif
#if defined(CL_VERSION_1_1)
        __INIT_CL_EXT_FCN_PTR(clEnqueueAcquireD3D10ObjectsKHR);
#endif
        
        cl_event tmp;
        cl_int err = detail::errHandler(
             pfn_clEnqueueAcquireD3D10ObjectsKHR(
                 object_,
                 (mem_objects != NULL) ? (cl_uint) mem_objects->size() : 0,
                 (mem_objects != NULL) ? (const cl_mem *) &mem_objects->front(): NULL,
                 (events != NULL) ? (cl_uint) events->size() : 0,
                 (events != NULL) ? (cl_event*) &events->front() : NULL,
                 (event != NULL) ? &tmp : NULL),
             __ENQUEUE_ACQUIRE_GL_ERR);

        if (event != NULL && err == CL_SUCCESS)
            *event = tmp;

        return err;
     }

    cl_int enqueueReleaseD3D10Objects(
         const VECTOR_CLASS<Memory>* mem_objects = NULL,
         const VECTOR_CLASS<Event>* events = NULL,
         Event* event = NULL) const
    {
        static PFN_clEnqueueReleaseD3D10ObjectsKHR pfn_clEnqueueReleaseD3D10ObjectsKHR = NULL;
#if defined(CL_VERSION_1_2)
        cl_context context = getInfo<CL_QUEUE_CONTEXT>();
        cl::Device device(getInfo<CL_QUEUE_DEVICE>());
        cl_platform_id platform = device.getInfo<CL_DEVICE_PLATFORM>();
        __INIT_CL_EXT_FCN_PTR_PLATFORM(platform, clEnqueueReleaseD3D10ObjectsKHR);
#endif // #if defined(CL_VERSION_1_2)
#if defined(CL_VERSION_1_1)
        __INIT_CL_EXT_FCN_PTR(clEnqueueReleaseD3D10ObjectsKHR);
#endif // #if defined(CL_VERSION_1_1)

        cl_event tmp;
        cl_int err = detail::errHandler(
            pfn_clEnqueueReleaseD3D10ObjectsKHR(
                object_,
                (mem_objects != NULL) ? (cl_uint) mem_objects->size() : 0,
                (mem_objects != NULL) ? (const cl_mem *) &mem_objects->front(): NULL,
                (events != NULL) ? (cl_uint) events->size() : 0,
                (events != NULL) ? (cl_event*) &events->front() : NULL,
                (event != NULL) ? &tmp : NULL),
            __ENQUEUE_RELEASE_GL_ERR);

        if (event != NULL && err == CL_SUCCESS)
            *event = tmp;

        return err;
    }
#endif

/**
 * Deprecated APIs for 1.2
 */
#if defined(CL_VERSION_1_1)
    cl_int enqueueBarrier() const
    {
        return detail::errHandler(
            ::clEnqueueBarrier(object_),
            __ENQUEUE_BARRIER_ERR);
    }
#endif // #if defined(CL_VERSION_1_1)

    cl_int flush() const
    {
        return detail::errHandler(::clFlush(object_), __FLUSH_ERR);
    }

    cl_int finish() const
    {
        return detail::errHandler(::clFinish(object_), __FINISH_ERR);
    }
};

__GET_INFO_HELPER_WITH_RETAIN(cl::CommandQueue)

#ifdef _WIN32
__declspec(selectany) volatile int CommandQueue::default_initialized_ = __DEFAULT_NOT_INITIALIZED;
__declspec(selectany) CommandQueue CommandQueue::default_;
__declspec(selectany) volatile cl_int CommandQueue::default_error_ = CL_SUCCESS;
#else
__attribute__((weak)) volatile int CommandQueue::default_initialized_ = __DEFAULT_NOT_INITIALIZED;
__attribute__((weak)) CommandQueue CommandQueue::default_;
__attribute__((weak)) volatile cl_int CommandQueue::default_error_ = CL_SUCCESS;
#endif

inline cl_int enqueueReadBuffer(
    const Buffer& buffer,
    cl_bool blocking,
    ::size_t offset,
    ::size_t size,
    void* ptr,
    const VECTOR_CLASS<Event>* events = NULL,
    Event* event = NULL)
{
    cl_int error;
    CommandQueue queue = CommandQueue::getDefault(&error);

    if (error != CL_SUCCESS) {
        return error;
    }

    return queue.enqueueReadBuffer(buffer, blocking, offset, size, ptr, events, event);
}

inline cl_int enqueueWriteBuffer(
        const Buffer& buffer,
        cl_bool blocking,
        ::size_t offset,
        ::size_t size,
        const void* ptr,
        const VECTOR_CLASS<Event>* events = NULL,
        Event* event = NULL)
{
    cl_int error;
    CommandQueue queue = CommandQueue::getDefault(&error);

    if (error != CL_SUCCESS) {
        return error;
    }

    return queue.enqueueWriteBuffer(buffer, blocking, offset, size, ptr, events, event);
}

inline void* enqueueMapBuffer(
        const Buffer& buffer,
        cl_bool blocking,
        cl_map_flags flags,
        ::size_t offset,
        ::size_t size,
        const VECTOR_CLASS<Event>* events = NULL,
        Event* event = NULL,
        cl_int* err = NULL)
{
    cl_int error;
    CommandQueue queue = CommandQueue::getDefault(&error);
    detail::errHandler(error, __ENQUEUE_MAP_BUFFER_ERR);
    if (err != NULL) {
        *err = error;
    }

    void * result = ::clEnqueueMapBuffer(
            queue(), buffer(), blocking, flags, offset, size,
            (events != NULL) ? (cl_uint) events->size() : 0,
            (events != NULL && events->size() > 0) ? (cl_event*) &events->front() : NULL,
            (cl_event*) event,
            &error);

    detail::errHandler(error, __ENQUEUE_MAP_BUFFER_ERR);
    if (err != NULL) {
        *err = error;
    }
    return result;
}

inline cl_int enqueueUnmapMemObject(
    const Memory& memory,
    void* mapped_ptr,
    const VECTOR_CLASS<Event>* events = NULL,
    Event* event = NULL)
{
    cl_int error;
    CommandQueue queue = CommandQueue::getDefault(&error);
    detail::errHandler(error, __ENQUEUE_MAP_BUFFER_ERR);
    if (error != CL_SUCCESS) {
        return error;
    }

    cl_event tmp;
    cl_int err = detail::errHandler(
        ::clEnqueueUnmapMemObject(
            queue(), memory(), mapped_ptr,
            (events != NULL) ? (cl_uint) events->size() : 0,
            (events != NULL && events->size() > 0) ? (cl_event*) &events->front() : NULL,
            (event != NULL) ? &tmp : NULL),
        __ENQUEUE_UNMAP_MEM_OBJECT_ERR);

    if (event != NULL && err == CL_SUCCESS)
        *event = tmp;

    return err;
}

inline cl_int enqueueCopyBuffer(
        const Buffer& src,
        const Buffer& dst,
        ::size_t src_offset,
        ::size_t dst_offset,
        ::size_t size,
        const VECTOR_CLASS<Event>* events = NULL,
        Event* event = NULL)
{
    cl_int error;
    CommandQueue queue = CommandQueue::getDefault(&error);

    if (error != CL_SUCCESS) {
        return error;
    }

    return queue.enqueueCopyBuffer(src, dst, src_offset, dst_offset, size, events, event);
}

/**
 * Blocking copy operation between iterators and a buffer.
 */
template< typename IteratorType >
inline cl_int copy( IteratorType startIterator, IteratorType endIterator, cl::Buffer &buffer )
{
    typedef typename std::iterator_traits<IteratorType>::value_type DataType;
    cl_int error;
    CommandQueue queue = CommandQueue::getDefault(&error);
        
    detail::errHandler(error, __COPY_ERR);
    if( error != CL_SUCCESS ) {
        return error;
    }
    
    ::size_t length = endIterator-startIterator;
    ::size_t byteLength = length*sizeof(DataType);

    DataType *pointer = 
        static_cast<DataType*>(enqueueMapBuffer(buffer, CL_TRUE, CL_MEM_WRITE_ONLY, 0, byteLength, 0, 0, &error));
    detail::errHandler(error, __COPY_ERR);
    if( error != CL_SUCCESS ) {
        return error;
    }
#if (_WIN32)
    std::copy(
        startIterator, 
        endIterator, 
        stdext::checked_array_iterator<DataType*>(
            pointer, length));
#else
    std::copy(startIterator, endIterator, pointer);
#endif
    enqueueUnmapMemObject(buffer, pointer, 0, 0);
    detail::errHandler(error, __COPY_ERR);
    if( error != CL_SUCCESS ) { 
        return error;
    }
    return CL_SUCCESS;
}

/**
 * Blocking copy operation between iterators and a buffer.
 */
template< typename IteratorType >
inline cl_int copy( cl::Buffer &buffer, IteratorType startIterator, IteratorType endIterator )
{
    typedef typename std::iterator_traits<IteratorType>::value_type DataType;
    cl_int error;
    CommandQueue queue = CommandQueue::getDefault(&error);
        
    detail::errHandler(error, __COPY_ERR);
    if( error != CL_SUCCESS ) {
        return error;
    }

    ::size_t length = endIterator-startIterator;
    ::size_t byteLength = length*sizeof(DataType);

    DataType *pointer = 
        static_cast<DataType*>(enqueueMapBuffer(buffer, CL_TRUE, CL_MEM_READ_ONLY, 0, byteLength, 0, 0, &error));
    detail::errHandler(error, __COPY_ERR);
    if( error != CL_SUCCESS ) {
        return error;
    }
    std::copy(pointer, pointer + length, startIterator);
    enqueueUnmapMemObject(buffer, pointer, 0, 0);
    detail::errHandler(error, __COPY_ERR);
    if( error != CL_SUCCESS ) {
        return error;
    }
    return CL_SUCCESS;
}

#if defined(CL_VERSION_1_1)
inline cl_int enqueueReadBufferRect(
    const Buffer& buffer,
    cl_bool blocking,
    const size_t<3>& buffer_offset,
    const size_t<3>& host_offset,
    const size_t<3>& region,
    ::size_t buffer_row_pitch,
    ::size_t buffer_slice_pitch,
    ::size_t host_row_pitch,
    ::size_t host_slice_pitch,
    void *ptr,
    const VECTOR_CLASS<Event>* events = NULL,
    Event* event = NULL)
{
    cl_int error;
    CommandQueue queue = CommandQueue::getDefault(&error);

    if (error != CL_SUCCESS) {
        return error;
    }

    return queue.enqueueReadBufferRect(
        buffer, 
        blocking, 
        buffer_offset, 
        host_offset,
        region,
        buffer_row_pitch,
        buffer_slice_pitch,
        host_row_pitch,
        host_slice_pitch,
        ptr, 
        events, 
        event);
}

inline cl_int enqueueWriteBufferRect(
    const Buffer& buffer,
    cl_bool blocking,
    const size_t<3>& buffer_offset,
    const size_t<3>& host_offset,
    const size_t<3>& region,
    ::size_t buffer_row_pitch,
    ::size_t buffer_slice_pitch,
    ::size_t host_row_pitch,
    ::size_t host_slice_pitch,
    void *ptr,
    const VECTOR_CLASS<Event>* events = NULL,
    Event* event = NULL)
{
    cl_int error;
    CommandQueue queue = CommandQueue::getDefault(&error);

    if (error != CL_SUCCESS) {
        return error;
    }

    return queue.enqueueWriteBufferRect(
        buffer, 
        blocking, 
        buffer_offset, 
        host_offset,
        region,
        buffer_row_pitch,
        buffer_slice_pitch,
        host_row_pitch,
        host_slice_pitch,
        ptr, 
        events, 
        event);
}

inline cl_int enqueueCopyBufferRect(
    const Buffer& src,
    const Buffer& dst,
    const size_t<3>& src_origin,
    const size_t<3>& dst_origin,
    const size_t<3>& region,
    ::size_t src_row_pitch,
    ::size_t src_slice_pitch,
    ::size_t dst_row_pitch,
    ::size_t dst_slice_pitch,
    const VECTOR_CLASS<Event>* events = NULL,
    Event* event = NULL)
{
    cl_int error;
    CommandQueue queue = CommandQueue::getDefault(&error);

    if (error != CL_SUCCESS) {
        return error;
    }

    return queue.enqueueCopyBufferRect(
        src,
        dst,
        src_origin,
        dst_origin,
        region,
        src_row_pitch,
        src_slice_pitch,
        dst_row_pitch,
        dst_slice_pitch,
        events, 
        event);
}
#endif

inline cl_int enqueueReadImage(
    const Image& image,
    cl_bool blocking,
    const size_t<3>& origin,
    const size_t<3>& region,
    ::size_t row_pitch,
    ::size_t slice_pitch,
    void* ptr,
    const VECTOR_CLASS<Event>* events = NULL,
    Event* event = NULL) 
{
    cl_int error;
    CommandQueue queue = CommandQueue::getDefault(&error);

    if (error != CL_SUCCESS) {
        return error;
    }

    return queue.enqueueReadImage(
        image,
        blocking,
        origin,
        region,
        row_pitch,
        slice_pitch,
        ptr,
        events, 
        event);
}

inline cl_int enqueueWriteImage(
    const Image& image,
    cl_bool blocking,
    const size_t<3>& origin,
    const size_t<3>& region,
    ::size_t row_pitch,
    ::size_t slice_pitch,
    void* ptr,
    const VECTOR_CLASS<Event>* events = NULL,
    Event* event = NULL)
{
    cl_int error;
    CommandQueue queue = CommandQueue::getDefault(&error);

    if (error != CL_SUCCESS) {
        return error;
    }

    return queue.enqueueWriteImage(
        image,
        blocking,
        origin,
        region,
        row_pitch,
        slice_pitch,
        ptr,
        events, 
        event);
}

inline cl_int enqueueCopyImage(
    const Image& src,
    const Image& dst,
    const size_t<3>& src_origin,
    const size_t<3>& dst_origin,
    const size_t<3>& region,
    const VECTOR_CLASS<Event>* events = NULL,
    Event* event = NULL)
{
    cl_int error;
    CommandQueue queue = CommandQueue::getDefault(&error);

    if (error != CL_SUCCESS) {
        return error;
    }

    return queue.enqueueCopyImage(
        src,
        dst,
        src_origin,
        dst_origin,
        region,
        events,
        event);
}

inline cl_int enqueueCopyImageToBuffer(
    const Image& src,
    const Buffer& dst,
    const size_t<3>& src_origin,
    const size_t<3>& region,
    ::size_t dst_offset,
    const VECTOR_CLASS<Event>* events = NULL,
    Event* event = NULL)
{
    cl_int error;
    CommandQueue queue = CommandQueue::getDefault(&error);

    if (error != CL_SUCCESS) {
        return error;
    }

    return queue.enqueueCopyImageToBuffer(
        src,
        dst,
        src_origin,
        region,
        dst_offset,
        events,
        event);
}

inline cl_int enqueueCopyBufferToImage(
    const Buffer& src,
    const Image& dst,
    ::size_t src_offset,
    const size_t<3>& dst_origin,
    const size_t<3>& region,
    const VECTOR_CLASS<Event>* events = NULL,
    Event* event = NULL)
{
    cl_int error;
    CommandQueue queue = CommandQueue::getDefault(&error);

    if (error != CL_SUCCESS) {
        return error;
    }

    return queue.enqueueCopyBufferToImage(
        src,
        dst,
        src_offset,
        dst_origin,
        region,
        events,
        event);
}


inline cl_int flush(void)
{
    cl_int error;
    CommandQueue queue = CommandQueue::getDefault(&error);

    if (error != CL_SUCCESS) {
        return error;
    }

    return queue.flush();
}

inline cl_int finish(void)
{
    cl_int error;
    CommandQueue queue = CommandQueue::getDefault(&error);

    if (error != CL_SUCCESS) {
        return error;
    } 


    return queue.finish();
}


// Kernel Functor support
// New interface as of September 2011
// Requires the C++11 std::tr1::function (note do not support TR1)
// Visual Studio 2010 and GCC 4.2

#if defined(CL_USE_CPP_FUNCTORS)
struct EnqueueArgs
{
    CommandQueue queue_;
    const NDRange offset_;
    const NDRange global_;
    const NDRange local_;

    EnqueueArgs(NDRange global) : 
      queue_(CommandQueue::getDefault()),
      offset_(NullRange), 
      global_(global),
      local_(NullRange)
    {

    }

    EnqueueArgs(NDRange global, NDRange local) : 
      queue_(CommandQueue::getDefault()),
      offset_(NullRange), 
      global_(global),
      local_(local)
    {

    }

    EnqueueArgs(NDRange offset, NDRange global, NDRange local) : 
      queue_(CommandQueue::getDefault()),
      offset_(offset), 
      global_(global),
      local_(local)
    {

    }

    EnqueueArgs(CommandQueue &queue, NDRange global) : 
      queue_(queue),
      offset_(NullRange), 
      global_(global),
      local_(NullRange)
    {

    }

    EnqueueArgs(CommandQueue &queue, NDRange global, NDRange local) : 
      queue_(queue),
      offset_(NullRange), 
      global_(global),
      local_(local)
    {

    }

    EnqueueArgs(CommandQueue &queue, NDRange offset, NDRange global, NDRange local) : 
      queue_(queue),
      offset_(offset), 
      global_(global),
      local_(local)
    {

    }
};

namespace detail {

class NullType {};

template<int index, typename T0>
struct SetArg
{
    static void set (Kernel kernel, T0 arg)
    {
        kernel.setArg(index, arg);
    }
};  


template<int index>
struct SetArg<index, NullType>
{
    static void set (Kernel, NullType)
    { 
    }
};


template <
   typename T0,   typename T1,   typename T2,   typename T3,
   typename T4,   typename T5,   typename T6,   typename T7,
   typename T8,   typename T9,   typename T10,   typename T11,
   typename T12,   typename T13,   typename T14,   typename T15,
   typename T16,   typename T17,   typename T18,   typename T19,
   typename T20,   typename T21,   typename T22,   typename T23,
   typename T24,   typename T25,   typename T26,   typename T27,
   typename T28,   typename T29,   typename T30,   typename T31
>
class KernelFunctorGlobal
{
private:
    Kernel kernel_;

public:
   KernelFunctorGlobal(
        Kernel kernel) :
            kernel_(kernel)
    {}

   KernelFunctorGlobal(
        const Program& program,
        const STRING_CLASS name,
        cl_int * err = NULL) :
            kernel_(program, name.c_str(), err)
    {}

    Event operator() (
        EnqueueArgs& args,
        T0 t0,
        T1 t1 = NullType(),
        T2 t2 = NullType(),
        T3 t3 = NullType(),
        T4 t4 = NullType(),
        T5 t5 = NullType(),
        T6 t6 = NullType(),
        T7 t7 = NullType(),
        T8 t8 = NullType(),
        T9 t9 = NullType(),
        T10 t10 = NullType(),
        T11 t11 = NullType(),
        T12 t12 = NullType(),
        T13 t13 = NullType(),
        T14 t14 = NullType(),
        T15 t15 = NullType(),
        T16 t16 = NullType(),
        T17 t17 = NullType(),
        T18 t18 = NullType(),
        T19 t19 = NullType(),
        T20 t20 = NullType(),
        T21 t21 = NullType(),
        T22 t22 = NullType(),
        T23 t23 = NullType(),
        T24 t24 = NullType(),
        T25 t25 = NullType(),
        T26 t26 = NullType(),
        T27 t27 = NullType(),
        T28 t28 = NullType(),
        T29 t29 = NullType(),
        T30 t30 = NullType(),
        T31 t31 = NullType()
        )
    {
        Event event;
        SetArg<0, T0>::set(kernel_, t0);
        SetArg<1, T1>::set(kernel_, t1);
        SetArg<2, T2>::set(kernel_, t2);
        SetArg<3, T3>::set(kernel_, t3);
        SetArg<4, T4>::set(kernel_, t4);
        SetArg<5, T5>::set(kernel_, t5);
        SetArg<6, T6>::set(kernel_, t6);
        SetArg<7, T7>::set(kernel_, t7);
        SetArg<8, T8>::set(kernel_, t8);
        SetArg<9, T9>::set(kernel_, t9);
        SetArg<10, T10>::set(kernel_, t10);
        SetArg<11, T11>::set(kernel_, t11);
        SetArg<12, T12>::set(kernel_, t12);
        SetArg<13, T13>::set(kernel_, t13);
        SetArg<14, T14>::set(kernel_, t14);
        SetArg<15, T15>::set(kernel_, t15);
        SetArg<16, T16>::set(kernel_, t16);
        SetArg<17, T17>::set(kernel_, t17);
        SetArg<18, T18>::set(kernel_, t18);
        SetArg<19, T19>::set(kernel_, t19);
        SetArg<20, T20>::set(kernel_, t20);
        SetArg<21, T21>::set(kernel_, t21);
        SetArg<22, T22>::set(kernel_, t22);
        SetArg<23, T23>::set(kernel_, t23);
        SetArg<24, T24>::set(kernel_, t24);
        SetArg<25, T25>::set(kernel_, t25);
        SetArg<26, T26>::set(kernel_, t26);
        SetArg<27, T27>::set(kernel_, t27);
        SetArg<28, T28>::set(kernel_, t28);
        SetArg<29, T29>::set(kernel_, t29);
        SetArg<30, T30>::set(kernel_, t30);
        SetArg<31, T31>::set(kernel_, t31);
        

        args.queue_.enqueueNDRangeKernel(
            kernel_,
            args.offset_,
            args.global_,
            args.local_,
            NULL,
            &event);
        

        return event;
    }

    Event operator() (
        EnqueueArgs& args,
        const Event& waitEvent,
        T0 t0,
        T1 t1 = NullType(),
        T2 t2 = NullType(),
        T3 t3 = NullType(),
        T4 t4 = NullType(),
        T5 t5 = NullType(),
        T6 t6 = NullType(),
        T7 t7 = NullType(),
        T8 t8 = NullType(),
        T9 t9 = NullType(),
        T10 t10 = NullType(),
        T11 t11 = NullType(),
        T12 t12 = NullType(),
        T13 t13 = NullType(),
        T14 t14 = NullType(),
        T15 t15 = NullType(),
        T16 t16 = NullType(),
        T17 t17 = NullType(),
        T18 t18 = NullType(),
        T19 t19 = NullType(),
        T20 t20 = NullType(),
        T21 t21 = NullType(),
        T22 t22 = NullType(),
        T23 t23 = NullType(),
        T24 t24 = NullType(),
        T25 t25 = NullType(),
        T26 t26 = NullType(),
        T27 t27 = NullType(),
        T28 t28 = NullType(),
        T29 t29 = NullType(),
        T30 t30 = NullType(),
        T31 t31 = NullType()
        )
    {
        Event event;
        SetArg<0, T0>::set(kernel_, t0);
        SetArg<1, T1>::set(kernel_, t1);
        SetArg<2, T2>::set(kernel_, t2);
        SetArg<3, T3>::set(kernel_, t3);
        SetArg<4, T4>::set(kernel_, t4);
        SetArg<5, T5>::set(kernel_, t5);
        SetArg<6, T6>::set(kernel_, t6);
        SetArg<7, T7>::set(kernel_, t7);
        SetArg<8, T8>::set(kernel_, t8);
        SetArg<9, T9>::set(kernel_, t9);
        SetArg<10, T10>::set(kernel_, t10);
        SetArg<11, T11>::set(kernel_, t11);
        SetArg<12, T12>::set(kernel_, t12);
        SetArg<13, T13>::set(kernel_, t13);
        SetArg<14, T14>::set(kernel_, t14);
        SetArg<15, T15>::set(kernel_, t15);
        SetArg<16, T16>::set(kernel_, t16);
        SetArg<17, T17>::set(kernel_, t17);
        SetArg<18, T18>::set(kernel_, t18);
        SetArg<19, T19>::set(kernel_, t19);
        SetArg<20, T20>::set(kernel_, t20);
        SetArg<21, T21>::set(kernel_, t21);
        SetArg<22, T22>::set(kernel_, t22);
        SetArg<23, T23>::set(kernel_, t23);
        SetArg<24, T24>::set(kernel_, t24);
        SetArg<25, T25>::set(kernel_, t25);
        SetArg<26, T26>::set(kernel_, t26);
        SetArg<27, T27>::set(kernel_, t27);
        SetArg<28, T28>::set(kernel_, t28);
        SetArg<29, T29>::set(kernel_, t29);
        SetArg<30, T30>::set(kernel_, t30);
        SetArg<31, T31>::set(kernel_, t31);


		VECTOR_CLASS<Event> events(&waitEvent, &waitEvent);

        args.queue_.enqueueNDRangeKernel(
            kernel_,
            args.offset_,
            args.global_,
            args.local_,
            &events,
            &event);


        return event;
    }

    Event operator() (
        EnqueueArgs& args,
        const VECTOR_CLASS<Event>& waitEvents,
        T0 t0,
        T1 t1 = NullType(),
        T2 t2 = NullType(),
        T3 t3 = NullType(),
        T4 t4 = NullType(),
        T5 t5 = NullType(),
        T6 t6 = NullType(),
        T7 t7 = NullType(),
        T8 t8 = NullType(),
        T9 t9 = NullType(),
        T10 t10 = NullType(),
        T11 t11 = NullType(),
        T12 t12 = NullType(),
        T13 t13 = NullType(),
        T14 t14 = NullType(),
        T15 t15 = NullType(),
        T16 t16 = NullType(),
        T17 t17 = NullType(),
        T18 t18 = NullType(),
        T19 t19 = NullType(),
        T20 t20 = NullType(),
        T21 t21 = NullType(),
        T22 t22 = NullType(),
        T23 t23 = NullType(),
        T24 t24 = NullType(),
        T25 t25 = NullType(),
        T26 t26 = NullType(),
        T27 t27 = NullType(),
        T28 t28 = NullType(),
        T29 t29 = NullType(),
        T30 t30 = NullType(),
        T31 t31 = NullType()
        )
    {
        Event event;
        SetArg<0, T0>::set(kernel_, t0);
        SetArg<1, T1>::set(kernel_, t1);
        SetArg<2, T2>::set(kernel_, t2);
        SetArg<3, T3>::set(kernel_, t3);
        SetArg<4, T4>::set(kernel_, t4);
        SetArg<5, T5>::set(kernel_, t5);
        SetArg<6, T6>::set(kernel_, t6);
        SetArg<7, T7>::set(kernel_, t7);
        SetArg<8, T8>::set(kernel_, t8);
        SetArg<9, T9>::set(kernel_, t9);
        SetArg<10, T10>::set(kernel_, t10);
        SetArg<11, T11>::set(kernel_, t11);
        SetArg<12, T12>::set(kernel_, t12);
        SetArg<13, T13>::set(kernel_, t13);
        SetArg<14, T14>::set(kernel_, t14);
        SetArg<15, T15>::set(kernel_, t15);
        SetArg<16, T16>::set(kernel_, t16);
        SetArg<17, T17>::set(kernel_, t17);
        SetArg<18, T18>::set(kernel_, t18);
        SetArg<19, T19>::set(kernel_, t19);
        SetArg<20, T20>::set(kernel_, t20);
        SetArg<21, T21>::set(kernel_, t21);
        SetArg<22, T22>::set(kernel_, t22);
        SetArg<23, T23>::set(kernel_, t23);
        SetArg<24, T24>::set(kernel_, t24);
        SetArg<25, T25>::set(kernel_, t25);
        SetArg<26, T26>::set(kernel_, t26);
        SetArg<27, T27>::set(kernel_, t27);
        SetArg<28, T28>::set(kernel_, t28);
        SetArg<29, T29>::set(kernel_, t29);
        SetArg<30, T30>::set(kernel_, t30);
        SetArg<31, T31>::set(kernel_, t31);


        args.queue_.enqueueNDRangeKernel(
            kernel_,
            args.offset_,
            args.global_,
            args.local_,
            &waitEvents,
            &event);


        return event;
    }
};

//------------------------------------------------------------------------------------------------------


template<
	typename T0,
	typename T1,
	typename T2,
	typename T3,
	typename T4,
	typename T5,
	typename T6,
	typename T7,
	typename T8,
	typename T9,
	typename T10,
	typename T11,
	typename T12,
	typename T13,
	typename T14,
	typename T15,
	typename T16,
	typename T17,
	typename T18,
	typename T19,
	typename T20,
	typename T21,
	typename T22,
	typename T23,
	typename T24,
	typename T25,
	typename T26,
	typename T27,
	typename T28,
	typename T29,
	typename T30,
	typename T31>
struct functionImplementation_
{
	typedef detail::KernelFunctorGlobal<
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13,
		T14,
		T15,
		T16,
		T17,
		T18,
		T19,
		T20,
		T21,
		T22,
		T23,
		T24,
		T25,
		T26,
		T27,
		T28,
		T29,
		T30,
		T31> FunctorType;

    FunctorType functor_;

    functionImplementation_(const FunctorType &functor) :
        functor_(functor)
    {
    
        #if (defined(_WIN32) && (_VARIADIC_MAX < 32))
        // Fail variadic expansion for dev11
        static_assert(0, "Visual Studio has a hard limit of argument count for a std::function expansion. Please define _VARIADIC_MAX to be 10. If you need more arguments than that VC12 and below cannot support it.");
        #endif
            
    }

    typedef __CL_FUNCTION_TYPE<Event (
        EnqueueArgs&,
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13,
		T14,
		T15,
		T16,
		T17,
		T18,
		T19,
		T20,
		T21,
		T22,
		T23,
		T24,
		T25,
		T26,
		T27,
		T28,
		T29,
		T30,
		T31)> type_;

	Event operator()(
		EnqueueArgs& enqueueArgs,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19,
		T20 arg20,
		T21 arg21,
		T22 arg22,
		T23 arg23,
		T24 arg24,
		T25 arg25,
		T26 arg26,
		T27 arg27,
		T28 arg28,
		T29 arg29,
		T30 arg30,
		T31 arg31)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19,
			arg20,
			arg21,
			arg22,
			arg23,
			arg24,
			arg25,
			arg26,
			arg27,
			arg28,
			arg29,
			arg30,
			arg31);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const Event &event,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19,
		T20 arg20,
		T21 arg21,
		T22 arg22,
		T23 arg23,
		T24 arg24,
		T25 arg25,
		T26 arg26,
		T27 arg27,
		T28 arg28,
		T29 arg29,
		T30 arg30,
		T31 arg31)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19,
			arg20,
			arg21,
			arg22,
			arg23,
			arg24,
			arg25,
			arg26,
			arg27,
			arg28,
			arg29,
			arg30,
			arg31);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const VECTOR_CLASS<Event>& waitEvents,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19,
		T20 arg20,
		T21 arg21,
		T22 arg22,
		T23 arg23,
		T24 arg24,
		T25 arg25,
		T26 arg26,
		T27 arg27,
		T28 arg28,
		T29 arg29,
		T30 arg30,
		T31 arg31)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19,
			arg20,
			arg21,
			arg22,
			arg23,
			arg24,
			arg25,
			arg26,
			arg27,
			arg28,
			arg29,
			arg30,
			arg31);
	}

    operator type_ ()    
    {
        return type_(*this);
    }

};

template<
	typename T0,
	typename T1,
	typename T2,
	typename T3,
	typename T4,
	typename T5,
	typename T6,
	typename T7,
	typename T8,
	typename T9,
	typename T10,
	typename T11,
	typename T12,
	typename T13,
	typename T14,
	typename T15,
	typename T16,
	typename T17,
	typename T18,
	typename T19,
	typename T20,
	typename T21,
	typename T22,
	typename T23,
	typename T24,
	typename T25,
	typename T26,
	typename T27,
	typename T28,
	typename T29,
	typename T30>
struct functionImplementation_
<	T0,
	T1,
	T2,
	T3,
	T4,
	T5,
	T6,
	T7,
	T8,
	T9,
	T10,
	T11,
	T12,
	T13,
	T14,
	T15,
	T16,
	T17,
	T18,
	T19,
	T20,
	T21,
	T22,
	T23,
	T24,
	T25,
	T26,
	T27,
	T28,
	T29,
	T30,
	NullType>
{
	typedef detail::KernelFunctorGlobal<
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13,
		T14,
		T15,
		T16,
		T17,
		T18,
		T19,
		T20,
		T21,
		T22,
		T23,
		T24,
		T25,
		T26,
		T27,
		T28,
		T29,
		T30,
		NullType> FunctorType;

    FunctorType functor_;

    functionImplementation_(const FunctorType &functor) :
        functor_(functor)
    {
    
        #if (defined(_WIN32) && (_VARIADIC_MAX < 31))
        // Fail variadic expansion for dev11
        static_assert(0, "Visual Studio has a hard limit of argument count for a std::function expansion. Please define _VARIADIC_MAX to be 10. If you need more arguments than that VC12 and below cannot support it.");
        #endif
            
    }

    typedef __CL_FUNCTION_TYPE<Event (
        EnqueueArgs&,
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13,
		T14,
		T15,
		T16,
		T17,
		T18,
		T19,
		T20,
		T21,
		T22,
		T23,
		T24,
		T25,
		T26,
		T27,
		T28,
		T29,
		T30)> type_;

	Event operator()(
		EnqueueArgs& enqueueArgs,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19,
		T20 arg20,
		T21 arg21,
		T22 arg22,
		T23 arg23,
		T24 arg24,
		T25 arg25,
		T26 arg26,
		T27 arg27,
		T28 arg28,
		T29 arg29,
		T30 arg30)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19,
			arg20,
			arg21,
			arg22,
			arg23,
			arg24,
			arg25,
			arg26,
			arg27,
			arg28,
			arg29,
			arg30);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const Event &event,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19,
		T20 arg20,
		T21 arg21,
		T22 arg22,
		T23 arg23,
		T24 arg24,
		T25 arg25,
		T26 arg26,
		T27 arg27,
		T28 arg28,
		T29 arg29,
		T30 arg30)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19,
			arg20,
			arg21,
			arg22,
			arg23,
			arg24,
			arg25,
			arg26,
			arg27,
			arg28,
			arg29,
			arg30);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const VECTOR_CLASS<Event>& waitEvents,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19,
		T20 arg20,
		T21 arg21,
		T22 arg22,
		T23 arg23,
		T24 arg24,
		T25 arg25,
		T26 arg26,
		T27 arg27,
		T28 arg28,
		T29 arg29,
		T30 arg30)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19,
			arg20,
			arg21,
			arg22,
			arg23,
			arg24,
			arg25,
			arg26,
			arg27,
			arg28,
			arg29,
			arg30);
	}

    operator type_ ()    
    {
        return type_(*this);
    }

};

template<
	typename T0,
	typename T1,
	typename T2,
	typename T3,
	typename T4,
	typename T5,
	typename T6,
	typename T7,
	typename T8,
	typename T9,
	typename T10,
	typename T11,
	typename T12,
	typename T13,
	typename T14,
	typename T15,
	typename T16,
	typename T17,
	typename T18,
	typename T19,
	typename T20,
	typename T21,
	typename T22,
	typename T23,
	typename T24,
	typename T25,
	typename T26,
	typename T27,
	typename T28,
	typename T29>
struct functionImplementation_
<	T0,
	T1,
	T2,
	T3,
	T4,
	T5,
	T6,
	T7,
	T8,
	T9,
	T10,
	T11,
	T12,
	T13,
	T14,
	T15,
	T16,
	T17,
	T18,
	T19,
	T20,
	T21,
	T22,
	T23,
	T24,
	T25,
	T26,
	T27,
	T28,
	T29,
	NullType,
	NullType>
{
	typedef detail::KernelFunctorGlobal<
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13,
		T14,
		T15,
		T16,
		T17,
		T18,
		T19,
		T20,
		T21,
		T22,
		T23,
		T24,
		T25,
		T26,
		T27,
		T28,
		T29,
		NullType,
		NullType> FunctorType;

    FunctorType functor_;

    functionImplementation_(const FunctorType &functor) :
        functor_(functor)
    {
    
        #if (defined(_WIN32) && (_VARIADIC_MAX < 30))
        // Fail variadic expansion for dev11
        static_assert(0, "Visual Studio has a hard limit of argument count for a std::function expansion. Please define _VARIADIC_MAX to be 10. If you need more arguments than that VC12 and below cannot support it.");
        #endif
            
    }

    typedef __CL_FUNCTION_TYPE<Event (
        EnqueueArgs&,
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13,
		T14,
		T15,
		T16,
		T17,
		T18,
		T19,
		T20,
		T21,
		T22,
		T23,
		T24,
		T25,
		T26,
		T27,
		T28,
		T29)> type_;

	Event operator()(
		EnqueueArgs& enqueueArgs,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19,
		T20 arg20,
		T21 arg21,
		T22 arg22,
		T23 arg23,
		T24 arg24,
		T25 arg25,
		T26 arg26,
		T27 arg27,
		T28 arg28,
		T29 arg29)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19,
			arg20,
			arg21,
			arg22,
			arg23,
			arg24,
			arg25,
			arg26,
			arg27,
			arg28,
			arg29);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const Event &event,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19,
		T20 arg20,
		T21 arg21,
		T22 arg22,
		T23 arg23,
		T24 arg24,
		T25 arg25,
		T26 arg26,
		T27 arg27,
		T28 arg28,
		T29 arg29)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19,
			arg20,
			arg21,
			arg22,
			arg23,
			arg24,
			arg25,
			arg26,
			arg27,
			arg28,
			arg29);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const VECTOR_CLASS<Event>& waitEvents,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19,
		T20 arg20,
		T21 arg21,
		T22 arg22,
		T23 arg23,
		T24 arg24,
		T25 arg25,
		T26 arg26,
		T27 arg27,
		T28 arg28,
		T29 arg29)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19,
			arg20,
			arg21,
			arg22,
			arg23,
			arg24,
			arg25,
			arg26,
			arg27,
			arg28,
			arg29);
	}

    operator type_ ()    
    {
        return type_(*this);
    }

};

template<
	typename T0,
	typename T1,
	typename T2,
	typename T3,
	typename T4,
	typename T5,
	typename T6,
	typename T7,
	typename T8,
	typename T9,
	typename T10,
	typename T11,
	typename T12,
	typename T13,
	typename T14,
	typename T15,
	typename T16,
	typename T17,
	typename T18,
	typename T19,
	typename T20,
	typename T21,
	typename T22,
	typename T23,
	typename T24,
	typename T25,
	typename T26,
	typename T27,
	typename T28>
struct functionImplementation_
<	T0,
	T1,
	T2,
	T3,
	T4,
	T5,
	T6,
	T7,
	T8,
	T9,
	T10,
	T11,
	T12,
	T13,
	T14,
	T15,
	T16,
	T17,
	T18,
	T19,
	T20,
	T21,
	T22,
	T23,
	T24,
	T25,
	T26,
	T27,
	T28,
	NullType,
	NullType,
	NullType>
{
	typedef detail::KernelFunctorGlobal<
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13,
		T14,
		T15,
		T16,
		T17,
		T18,
		T19,
		T20,
		T21,
		T22,
		T23,
		T24,
		T25,
		T26,
		T27,
		T28,
		NullType,
		NullType,
		NullType> FunctorType;

    FunctorType functor_;

    functionImplementation_(const FunctorType &functor) :
        functor_(functor)
    {
    
        #if (defined(_WIN32) && (_VARIADIC_MAX < 29))
        // Fail variadic expansion for dev11
        static_assert(0, "Visual Studio has a hard limit of argument count for a std::function expansion. Please define _VARIADIC_MAX to be 10. If you need more arguments than that VC12 and below cannot support it.");
        #endif
            
    }

    typedef __CL_FUNCTION_TYPE<Event (
        EnqueueArgs&,
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13,
		T14,
		T15,
		T16,
		T17,
		T18,
		T19,
		T20,
		T21,
		T22,
		T23,
		T24,
		T25,
		T26,
		T27,
		T28)> type_;

	Event operator()(
		EnqueueArgs& enqueueArgs,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19,
		T20 arg20,
		T21 arg21,
		T22 arg22,
		T23 arg23,
		T24 arg24,
		T25 arg25,
		T26 arg26,
		T27 arg27,
		T28 arg28)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19,
			arg20,
			arg21,
			arg22,
			arg23,
			arg24,
			arg25,
			arg26,
			arg27,
			arg28);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const Event &event,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19,
		T20 arg20,
		T21 arg21,
		T22 arg22,
		T23 arg23,
		T24 arg24,
		T25 arg25,
		T26 arg26,
		T27 arg27,
		T28 arg28)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19,
			arg20,
			arg21,
			arg22,
			arg23,
			arg24,
			arg25,
			arg26,
			arg27,
			arg28);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const VECTOR_CLASS<Event>& waitEvents,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19,
		T20 arg20,
		T21 arg21,
		T22 arg22,
		T23 arg23,
		T24 arg24,
		T25 arg25,
		T26 arg26,
		T27 arg27,
		T28 arg28)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19,
			arg20,
			arg21,
			arg22,
			arg23,
			arg24,
			arg25,
			arg26,
			arg27,
			arg28);
	}

    operator type_ ()    
    {
        return type_(*this);
    }

};

template<
	typename T0,
	typename T1,
	typename T2,
	typename T3,
	typename T4,
	typename T5,
	typename T6,
	typename T7,
	typename T8,
	typename T9,
	typename T10,
	typename T11,
	typename T12,
	typename T13,
	typename T14,
	typename T15,
	typename T16,
	typename T17,
	typename T18,
	typename T19,
	typename T20,
	typename T21,
	typename T22,
	typename T23,
	typename T24,
	typename T25,
	typename T26,
	typename T27>
struct functionImplementation_
<	T0,
	T1,
	T2,
	T3,
	T4,
	T5,
	T6,
	T7,
	T8,
	T9,
	T10,
	T11,
	T12,
	T13,
	T14,
	T15,
	T16,
	T17,
	T18,
	T19,
	T20,
	T21,
	T22,
	T23,
	T24,
	T25,
	T26,
	T27,
	NullType,
	NullType,
	NullType,
	NullType>
{
	typedef detail::KernelFunctorGlobal<
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13,
		T14,
		T15,
		T16,
		T17,
		T18,
		T19,
		T20,
		T21,
		T22,
		T23,
		T24,
		T25,
		T26,
		T27,
		NullType,
		NullType,
		NullType,
		NullType> FunctorType;

    FunctorType functor_;

    functionImplementation_(const FunctorType &functor) :
        functor_(functor)
    {
    
        #if (defined(_WIN32) && (_VARIADIC_MAX < 28))
        // Fail variadic expansion for dev11
        static_assert(0, "Visual Studio has a hard limit of argument count for a std::function expansion. Please define _VARIADIC_MAX to be 10. If you need more arguments than that VC12 and below cannot support it.");
        #endif
            
    }

    typedef __CL_FUNCTION_TYPE<Event (
        EnqueueArgs&,
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13,
		T14,
		T15,
		T16,
		T17,
		T18,
		T19,
		T20,
		T21,
		T22,
		T23,
		T24,
		T25,
		T26,
		T27)> type_;

	Event operator()(
		EnqueueArgs& enqueueArgs,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19,
		T20 arg20,
		T21 arg21,
		T22 arg22,
		T23 arg23,
		T24 arg24,
		T25 arg25,
		T26 arg26,
		T27 arg27)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19,
			arg20,
			arg21,
			arg22,
			arg23,
			arg24,
			arg25,
			arg26,
			arg27);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const Event &event,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19,
		T20 arg20,
		T21 arg21,
		T22 arg22,
		T23 arg23,
		T24 arg24,
		T25 arg25,
		T26 arg26,
		T27 arg27)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19,
			arg20,
			arg21,
			arg22,
			arg23,
			arg24,
			arg25,
			arg26,
			arg27);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const VECTOR_CLASS<Event>& waitEvents,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19,
		T20 arg20,
		T21 arg21,
		T22 arg22,
		T23 arg23,
		T24 arg24,
		T25 arg25,
		T26 arg26,
		T27 arg27)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19,
			arg20,
			arg21,
			arg22,
			arg23,
			arg24,
			arg25,
			arg26,
			arg27);
	}

    operator type_ ()    
    {
        return type_(*this);
    }

};

template<
	typename T0,
	typename T1,
	typename T2,
	typename T3,
	typename T4,
	typename T5,
	typename T6,
	typename T7,
	typename T8,
	typename T9,
	typename T10,
	typename T11,
	typename T12,
	typename T13,
	typename T14,
	typename T15,
	typename T16,
	typename T17,
	typename T18,
	typename T19,
	typename T20,
	typename T21,
	typename T22,
	typename T23,
	typename T24,
	typename T25,
	typename T26>
struct functionImplementation_
<	T0,
	T1,
	T2,
	T3,
	T4,
	T5,
	T6,
	T7,
	T8,
	T9,
	T10,
	T11,
	T12,
	T13,
	T14,
	T15,
	T16,
	T17,
	T18,
	T19,
	T20,
	T21,
	T22,
	T23,
	T24,
	T25,
	T26,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType>
{
	typedef detail::KernelFunctorGlobal<
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13,
		T14,
		T15,
		T16,
		T17,
		T18,
		T19,
		T20,
		T21,
		T22,
		T23,
		T24,
		T25,
		T26,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType> FunctorType;

    FunctorType functor_;

    functionImplementation_(const FunctorType &functor) :
        functor_(functor)
    {
    
        #if (defined(_WIN32) && (_VARIADIC_MAX < 27))
        // Fail variadic expansion for dev11
        static_assert(0, "Visual Studio has a hard limit of argument count for a std::function expansion. Please define _VARIADIC_MAX to be 10. If you need more arguments than that VC12 and below cannot support it.");
        #endif
            
    }

    typedef __CL_FUNCTION_TYPE<Event (
        EnqueueArgs&,
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13,
		T14,
		T15,
		T16,
		T17,
		T18,
		T19,
		T20,
		T21,
		T22,
		T23,
		T24,
		T25,
		T26)> type_;

	Event operator()(
		EnqueueArgs& enqueueArgs,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19,
		T20 arg20,
		T21 arg21,
		T22 arg22,
		T23 arg23,
		T24 arg24,
		T25 arg25,
		T26 arg26)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19,
			arg20,
			arg21,
			arg22,
			arg23,
			arg24,
			arg25,
			arg26);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const Event &event,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19,
		T20 arg20,
		T21 arg21,
		T22 arg22,
		T23 arg23,
		T24 arg24,
		T25 arg25,
		T26 arg26)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19,
			arg20,
			arg21,
			arg22,
			arg23,
			arg24,
			arg25,
			arg26);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const VECTOR_CLASS<Event>& waitEvents,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19,
		T20 arg20,
		T21 arg21,
		T22 arg22,
		T23 arg23,
		T24 arg24,
		T25 arg25,
		T26 arg26)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19,
			arg20,
			arg21,
			arg22,
			arg23,
			arg24,
			arg25,
			arg26);
	}

    operator type_ ()    
    {
        return type_(*this);
    }

};

template<
	typename T0,
	typename T1,
	typename T2,
	typename T3,
	typename T4,
	typename T5,
	typename T6,
	typename T7,
	typename T8,
	typename T9,
	typename T10,
	typename T11,
	typename T12,
	typename T13,
	typename T14,
	typename T15,
	typename T16,
	typename T17,
	typename T18,
	typename T19,
	typename T20,
	typename T21,
	typename T22,
	typename T23,
	typename T24,
	typename T25>
struct functionImplementation_
<	T0,
	T1,
	T2,
	T3,
	T4,
	T5,
	T6,
	T7,
	T8,
	T9,
	T10,
	T11,
	T12,
	T13,
	T14,
	T15,
	T16,
	T17,
	T18,
	T19,
	T20,
	T21,
	T22,
	T23,
	T24,
	T25,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType>
{
	typedef detail::KernelFunctorGlobal<
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13,
		T14,
		T15,
		T16,
		T17,
		T18,
		T19,
		T20,
		T21,
		T22,
		T23,
		T24,
		T25,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType> FunctorType;

    FunctorType functor_;

    functionImplementation_(const FunctorType &functor) :
        functor_(functor)
    {
    
        #if (defined(_WIN32) && (_VARIADIC_MAX < 26))
        // Fail variadic expansion for dev11
        static_assert(0, "Visual Studio has a hard limit of argument count for a std::function expansion. Please define _VARIADIC_MAX to be 10. If you need more arguments than that VC12 and below cannot support it.");
        #endif
            
    }

    typedef __CL_FUNCTION_TYPE<Event (
        EnqueueArgs&,
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13,
		T14,
		T15,
		T16,
		T17,
		T18,
		T19,
		T20,
		T21,
		T22,
		T23,
		T24,
		T25)> type_;

	Event operator()(
		EnqueueArgs& enqueueArgs,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19,
		T20 arg20,
		T21 arg21,
		T22 arg22,
		T23 arg23,
		T24 arg24,
		T25 arg25)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19,
			arg20,
			arg21,
			arg22,
			arg23,
			arg24,
			arg25);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const Event &event,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19,
		T20 arg20,
		T21 arg21,
		T22 arg22,
		T23 arg23,
		T24 arg24,
		T25 arg25)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19,
			arg20,
			arg21,
			arg22,
			arg23,
			arg24,
			arg25);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const VECTOR_CLASS<Event>& waitEvents,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19,
		T20 arg20,
		T21 arg21,
		T22 arg22,
		T23 arg23,
		T24 arg24,
		T25 arg25)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19,
			arg20,
			arg21,
			arg22,
			arg23,
			arg24,
			arg25);
	}

    operator type_ ()    
    {
        return type_(*this);
    }

};

template<
	typename T0,
	typename T1,
	typename T2,
	typename T3,
	typename T4,
	typename T5,
	typename T6,
	typename T7,
	typename T8,
	typename T9,
	typename T10,
	typename T11,
	typename T12,
	typename T13,
	typename T14,
	typename T15,
	typename T16,
	typename T17,
	typename T18,
	typename T19,
	typename T20,
	typename T21,
	typename T22,
	typename T23,
	typename T24>
struct functionImplementation_
<	T0,
	T1,
	T2,
	T3,
	T4,
	T5,
	T6,
	T7,
	T8,
	T9,
	T10,
	T11,
	T12,
	T13,
	T14,
	T15,
	T16,
	T17,
	T18,
	T19,
	T20,
	T21,
	T22,
	T23,
	T24,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType>
{
	typedef detail::KernelFunctorGlobal<
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13,
		T14,
		T15,
		T16,
		T17,
		T18,
		T19,
		T20,
		T21,
		T22,
		T23,
		T24,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType> FunctorType;

    FunctorType functor_;

    functionImplementation_(const FunctorType &functor) :
        functor_(functor)
    {
    
        #if (defined(_WIN32) && (_VARIADIC_MAX < 25))
        // Fail variadic expansion for dev11
        static_assert(0, "Visual Studio has a hard limit of argument count for a std::function expansion. Please define _VARIADIC_MAX to be 10. If you need more arguments than that VC12 and below cannot support it.");
        #endif
            
    }

    typedef __CL_FUNCTION_TYPE<Event (
        EnqueueArgs&,
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13,
		T14,
		T15,
		T16,
		T17,
		T18,
		T19,
		T20,
		T21,
		T22,
		T23,
		T24)> type_;

	Event operator()(
		EnqueueArgs& enqueueArgs,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19,
		T20 arg20,
		T21 arg21,
		T22 arg22,
		T23 arg23,
		T24 arg24)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19,
			arg20,
			arg21,
			arg22,
			arg23,
			arg24);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const Event &event,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19,
		T20 arg20,
		T21 arg21,
		T22 arg22,
		T23 arg23,
		T24 arg24)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19,
			arg20,
			arg21,
			arg22,
			arg23,
			arg24);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const VECTOR_CLASS<Event>& waitEvents,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19,
		T20 arg20,
		T21 arg21,
		T22 arg22,
		T23 arg23,
		T24 arg24)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19,
			arg20,
			arg21,
			arg22,
			arg23,
			arg24);
	}

    operator type_ ()    
    {
        return type_(*this);
    }

};

template<
	typename T0,
	typename T1,
	typename T2,
	typename T3,
	typename T4,
	typename T5,
	typename T6,
	typename T7,
	typename T8,
	typename T9,
	typename T10,
	typename T11,
	typename T12,
	typename T13,
	typename T14,
	typename T15,
	typename T16,
	typename T17,
	typename T18,
	typename T19,
	typename T20,
	typename T21,
	typename T22,
	typename T23>
struct functionImplementation_
<	T0,
	T1,
	T2,
	T3,
	T4,
	T5,
	T6,
	T7,
	T8,
	T9,
	T10,
	T11,
	T12,
	T13,
	T14,
	T15,
	T16,
	T17,
	T18,
	T19,
	T20,
	T21,
	T22,
	T23,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType>
{
	typedef detail::KernelFunctorGlobal<
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13,
		T14,
		T15,
		T16,
		T17,
		T18,
		T19,
		T20,
		T21,
		T22,
		T23,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType> FunctorType;

    FunctorType functor_;

    functionImplementation_(const FunctorType &functor) :
        functor_(functor)
    {
    
        #if (defined(_WIN32) && (_VARIADIC_MAX < 24))
        // Fail variadic expansion for dev11
        static_assert(0, "Visual Studio has a hard limit of argument count for a std::function expansion. Please define _VARIADIC_MAX to be 10. If you need more arguments than that VC12 and below cannot support it.");
        #endif
            
    }

    typedef __CL_FUNCTION_TYPE<Event (
        EnqueueArgs&,
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13,
		T14,
		T15,
		T16,
		T17,
		T18,
		T19,
		T20,
		T21,
		T22,
		T23)> type_;

	Event operator()(
		EnqueueArgs& enqueueArgs,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19,
		T20 arg20,
		T21 arg21,
		T22 arg22,
		T23 arg23)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19,
			arg20,
			arg21,
			arg22,
			arg23);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const Event &event,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19,
		T20 arg20,
		T21 arg21,
		T22 arg22,
		T23 arg23)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19,
			arg20,
			arg21,
			arg22,
			arg23);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const VECTOR_CLASS<Event>& waitEvents,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19,
		T20 arg20,
		T21 arg21,
		T22 arg22,
		T23 arg23)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19,
			arg20,
			arg21,
			arg22,
			arg23);
	}

    operator type_ ()    
    {
        return type_(*this);
    }

};

template<
	typename T0,
	typename T1,
	typename T2,
	typename T3,
	typename T4,
	typename T5,
	typename T6,
	typename T7,
	typename T8,
	typename T9,
	typename T10,
	typename T11,
	typename T12,
	typename T13,
	typename T14,
	typename T15,
	typename T16,
	typename T17,
	typename T18,
	typename T19,
	typename T20,
	typename T21,
	typename T22>
struct functionImplementation_
<	T0,
	T1,
	T2,
	T3,
	T4,
	T5,
	T6,
	T7,
	T8,
	T9,
	T10,
	T11,
	T12,
	T13,
	T14,
	T15,
	T16,
	T17,
	T18,
	T19,
	T20,
	T21,
	T22,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType>
{
	typedef detail::KernelFunctorGlobal<
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13,
		T14,
		T15,
		T16,
		T17,
		T18,
		T19,
		T20,
		T21,
		T22,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType> FunctorType;

    FunctorType functor_;

    functionImplementation_(const FunctorType &functor) :
        functor_(functor)
    {
    
        #if (defined(_WIN32) && (_VARIADIC_MAX < 23))
        // Fail variadic expansion for dev11
        static_assert(0, "Visual Studio has a hard limit of argument count for a std::function expansion. Please define _VARIADIC_MAX to be 10. If you need more arguments than that VC12 and below cannot support it.");
        #endif
            
    }

    typedef __CL_FUNCTION_TYPE<Event (
        EnqueueArgs&,
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13,
		T14,
		T15,
		T16,
		T17,
		T18,
		T19,
		T20,
		T21,
		T22)> type_;

	Event operator()(
		EnqueueArgs& enqueueArgs,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19,
		T20 arg20,
		T21 arg21,
		T22 arg22)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19,
			arg20,
			arg21,
			arg22);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const Event &event,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19,
		T20 arg20,
		T21 arg21,
		T22 arg22)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19,
			arg20,
			arg21,
			arg22);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const VECTOR_CLASS<Event>& waitEvents,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19,
		T20 arg20,
		T21 arg21,
		T22 arg22)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19,
			arg20,
			arg21,
			arg22);
	}

    operator type_ ()    
    {
        return type_(*this);
    }

};

template<
	typename T0,
	typename T1,
	typename T2,
	typename T3,
	typename T4,
	typename T5,
	typename T6,
	typename T7,
	typename T8,
	typename T9,
	typename T10,
	typename T11,
	typename T12,
	typename T13,
	typename T14,
	typename T15,
	typename T16,
	typename T17,
	typename T18,
	typename T19,
	typename T20,
	typename T21>
struct functionImplementation_
<	T0,
	T1,
	T2,
	T3,
	T4,
	T5,
	T6,
	T7,
	T8,
	T9,
	T10,
	T11,
	T12,
	T13,
	T14,
	T15,
	T16,
	T17,
	T18,
	T19,
	T20,
	T21,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType>
{
	typedef detail::KernelFunctorGlobal<
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13,
		T14,
		T15,
		T16,
		T17,
		T18,
		T19,
		T20,
		T21,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType> FunctorType;

    FunctorType functor_;

    functionImplementation_(const FunctorType &functor) :
        functor_(functor)
    {
    
        #if (defined(_WIN32) && (_VARIADIC_MAX < 22))
        // Fail variadic expansion for dev11
        static_assert(0, "Visual Studio has a hard limit of argument count for a std::function expansion. Please define _VARIADIC_MAX to be 10. If you need more arguments than that VC12 and below cannot support it.");
        #endif
            
    }

    typedef __CL_FUNCTION_TYPE<Event (
        EnqueueArgs&,
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13,
		T14,
		T15,
		T16,
		T17,
		T18,
		T19,
		T20,
		T21)> type_;

	Event operator()(
		EnqueueArgs& enqueueArgs,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19,
		T20 arg20,
		T21 arg21)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19,
			arg20,
			arg21);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const Event &event,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19,
		T20 arg20,
		T21 arg21)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19,
			arg20,
			arg21);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const VECTOR_CLASS<Event>& waitEvents,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19,
		T20 arg20,
		T21 arg21)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19,
			arg20,
			arg21);
	}

    operator type_ ()    
    {
        return type_(*this);
    }

};

template<
	typename T0,
	typename T1,
	typename T2,
	typename T3,
	typename T4,
	typename T5,
	typename T6,
	typename T7,
	typename T8,
	typename T9,
	typename T10,
	typename T11,
	typename T12,
	typename T13,
	typename T14,
	typename T15,
	typename T16,
	typename T17,
	typename T18,
	typename T19,
	typename T20>
struct functionImplementation_
<	T0,
	T1,
	T2,
	T3,
	T4,
	T5,
	T6,
	T7,
	T8,
	T9,
	T10,
	T11,
	T12,
	T13,
	T14,
	T15,
	T16,
	T17,
	T18,
	T19,
	T20,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType>
{
	typedef detail::KernelFunctorGlobal<
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13,
		T14,
		T15,
		T16,
		T17,
		T18,
		T19,
		T20,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType> FunctorType;

    FunctorType functor_;

    functionImplementation_(const FunctorType &functor) :
        functor_(functor)
    {
    
        #if (defined(_WIN32) && (_VARIADIC_MAX < 21))
        // Fail variadic expansion for dev11
        static_assert(0, "Visual Studio has a hard limit of argument count for a std::function expansion. Please define _VARIADIC_MAX to be 10. If you need more arguments than that VC12 and below cannot support it.");
        #endif
            
    }

    typedef __CL_FUNCTION_TYPE<Event (
        EnqueueArgs&,
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13,
		T14,
		T15,
		T16,
		T17,
		T18,
		T19,
		T20)> type_;

	Event operator()(
		EnqueueArgs& enqueueArgs,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19,
		T20 arg20)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19,
			arg20);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const Event &event,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19,
		T20 arg20)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19,
			arg20);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const VECTOR_CLASS<Event>& waitEvents,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19,
		T20 arg20)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19,
			arg20);
	}

    operator type_ ()    
    {
        return type_(*this);
    }

};

template<
	typename T0,
	typename T1,
	typename T2,
	typename T3,
	typename T4,
	typename T5,
	typename T6,
	typename T7,
	typename T8,
	typename T9,
	typename T10,
	typename T11,
	typename T12,
	typename T13,
	typename T14,
	typename T15,
	typename T16,
	typename T17,
	typename T18,
	typename T19>
struct functionImplementation_
<	T0,
	T1,
	T2,
	T3,
	T4,
	T5,
	T6,
	T7,
	T8,
	T9,
	T10,
	T11,
	T12,
	T13,
	T14,
	T15,
	T16,
	T17,
	T18,
	T19,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType>
{
	typedef detail::KernelFunctorGlobal<
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13,
		T14,
		T15,
		T16,
		T17,
		T18,
		T19,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType> FunctorType;

    FunctorType functor_;

    functionImplementation_(const FunctorType &functor) :
        functor_(functor)
    {
    
        #if (defined(_WIN32) && (_VARIADIC_MAX < 20))
        // Fail variadic expansion for dev11
        static_assert(0, "Visual Studio has a hard limit of argument count for a std::function expansion. Please define _VARIADIC_MAX to be 10. If you need more arguments than that VC12 and below cannot support it.");
        #endif
            
    }

    typedef __CL_FUNCTION_TYPE<Event (
        EnqueueArgs&,
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13,
		T14,
		T15,
		T16,
		T17,
		T18,
		T19)> type_;

	Event operator()(
		EnqueueArgs& enqueueArgs,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const Event &event,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const VECTOR_CLASS<Event>& waitEvents,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18,
		T19 arg19)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18,
			arg19);
	}

    operator type_ ()    
    {
        return type_(*this);
    }

};

template<
	typename T0,
	typename T1,
	typename T2,
	typename T3,
	typename T4,
	typename T5,
	typename T6,
	typename T7,
	typename T8,
	typename T9,
	typename T10,
	typename T11,
	typename T12,
	typename T13,
	typename T14,
	typename T15,
	typename T16,
	typename T17,
	typename T18>
struct functionImplementation_
<	T0,
	T1,
	T2,
	T3,
	T4,
	T5,
	T6,
	T7,
	T8,
	T9,
	T10,
	T11,
	T12,
	T13,
	T14,
	T15,
	T16,
	T17,
	T18,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType>
{
	typedef detail::KernelFunctorGlobal<
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13,
		T14,
		T15,
		T16,
		T17,
		T18,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType> FunctorType;

    FunctorType functor_;

    functionImplementation_(const FunctorType &functor) :
        functor_(functor)
    {
    
        #if (defined(_WIN32) && (_VARIADIC_MAX < 19))
        // Fail variadic expansion for dev11
        static_assert(0, "Visual Studio has a hard limit of argument count for a std::function expansion. Please define _VARIADIC_MAX to be 10. If you need more arguments than that VC12 and below cannot support it.");
        #endif
            
    }

    typedef __CL_FUNCTION_TYPE<Event (
        EnqueueArgs&,
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13,
		T14,
		T15,
		T16,
		T17,
		T18)> type_;

	Event operator()(
		EnqueueArgs& enqueueArgs,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const Event &event,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const VECTOR_CLASS<Event>& waitEvents,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17,
		T18 arg18)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17,
			arg18);
	}

    operator type_ ()    
    {
        return type_(*this);
    }

};

template<
	typename T0,
	typename T1,
	typename T2,
	typename T3,
	typename T4,
	typename T5,
	typename T6,
	typename T7,
	typename T8,
	typename T9,
	typename T10,
	typename T11,
	typename T12,
	typename T13,
	typename T14,
	typename T15,
	typename T16,
	typename T17>
struct functionImplementation_
<	T0,
	T1,
	T2,
	T3,
	T4,
	T5,
	T6,
	T7,
	T8,
	T9,
	T10,
	T11,
	T12,
	T13,
	T14,
	T15,
	T16,
	T17,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType>
{
	typedef detail::KernelFunctorGlobal<
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13,
		T14,
		T15,
		T16,
		T17,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType> FunctorType;

    FunctorType functor_;

    functionImplementation_(const FunctorType &functor) :
        functor_(functor)
    {
    
        #if (defined(_WIN32) && (_VARIADIC_MAX < 18))
        // Fail variadic expansion for dev11
        static_assert(0, "Visual Studio has a hard limit of argument count for a std::function expansion. Please define _VARIADIC_MAX to be 10. If you need more arguments than that VC12 and below cannot support it.");
        #endif
            
    }

    typedef __CL_FUNCTION_TYPE<Event (
        EnqueueArgs&,
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13,
		T14,
		T15,
		T16,
		T17)> type_;

	Event operator()(
		EnqueueArgs& enqueueArgs,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const Event &event,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const VECTOR_CLASS<Event>& waitEvents,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16,
		T17 arg17)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16,
			arg17);
	}

    operator type_ ()    
    {
        return type_(*this);
    }

};

template<
	typename T0,
	typename T1,
	typename T2,
	typename T3,
	typename T4,
	typename T5,
	typename T6,
	typename T7,
	typename T8,
	typename T9,
	typename T10,
	typename T11,
	typename T12,
	typename T13,
	typename T14,
	typename T15,
	typename T16>
struct functionImplementation_
<	T0,
	T1,
	T2,
	T3,
	T4,
	T5,
	T6,
	T7,
	T8,
	T9,
	T10,
	T11,
	T12,
	T13,
	T14,
	T15,
	T16,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType>
{
	typedef detail::KernelFunctorGlobal<
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13,
		T14,
		T15,
		T16,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType> FunctorType;

    FunctorType functor_;

    functionImplementation_(const FunctorType &functor) :
        functor_(functor)
    {
    
        #if (defined(_WIN32) && (_VARIADIC_MAX < 17))
        // Fail variadic expansion for dev11
        static_assert(0, "Visual Studio has a hard limit of argument count for a std::function expansion. Please define _VARIADIC_MAX to be 10. If you need more arguments than that VC12 and below cannot support it.");
        #endif
            
    }

    typedef __CL_FUNCTION_TYPE<Event (
        EnqueueArgs&,
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13,
		T14,
		T15,
		T16)> type_;

	Event operator()(
		EnqueueArgs& enqueueArgs,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const Event &event,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const VECTOR_CLASS<Event>& waitEvents,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15,
		T16 arg16)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15,
			arg16);
	}

    operator type_ ()    
    {
        return type_(*this);
    }

};

template<
	typename T0,
	typename T1,
	typename T2,
	typename T3,
	typename T4,
	typename T5,
	typename T6,
	typename T7,
	typename T8,
	typename T9,
	typename T10,
	typename T11,
	typename T12,
	typename T13,
	typename T14,
	typename T15>
struct functionImplementation_
<	T0,
	T1,
	T2,
	T3,
	T4,
	T5,
	T6,
	T7,
	T8,
	T9,
	T10,
	T11,
	T12,
	T13,
	T14,
	T15,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType>
{
	typedef detail::KernelFunctorGlobal<
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13,
		T14,
		T15,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType> FunctorType;

    FunctorType functor_;

    functionImplementation_(const FunctorType &functor) :
        functor_(functor)
    {
    
        #if (defined(_WIN32) && (_VARIADIC_MAX < 16))
        // Fail variadic expansion for dev11
        static_assert(0, "Visual Studio has a hard limit of argument count for a std::function expansion. Please define _VARIADIC_MAX to be 10. If you need more arguments than that VC12 and below cannot support it.");
        #endif
            
    }

    typedef __CL_FUNCTION_TYPE<Event (
        EnqueueArgs&,
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13,
		T14,
		T15)> type_;

	Event operator()(
		EnqueueArgs& enqueueArgs,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const Event &event,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const VECTOR_CLASS<Event>& waitEvents,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14,
		T15 arg15)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14,
			arg15);
	}

    operator type_ ()    
    {
        return type_(*this);
    }

};

template<
	typename T0,
	typename T1,
	typename T2,
	typename T3,
	typename T4,
	typename T5,
	typename T6,
	typename T7,
	typename T8,
	typename T9,
	typename T10,
	typename T11,
	typename T12,
	typename T13,
	typename T14>
struct functionImplementation_
<	T0,
	T1,
	T2,
	T3,
	T4,
	T5,
	T6,
	T7,
	T8,
	T9,
	T10,
	T11,
	T12,
	T13,
	T14,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType>
{
	typedef detail::KernelFunctorGlobal<
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13,
		T14,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType> FunctorType;

    FunctorType functor_;

    functionImplementation_(const FunctorType &functor) :
        functor_(functor)
    {
    
        #if (defined(_WIN32) && (_VARIADIC_MAX < 15))
        // Fail variadic expansion for dev11
        static_assert(0, "Visual Studio has a hard limit of argument count for a std::function expansion. Please define _VARIADIC_MAX to be 10. If you need more arguments than that VC12 and below cannot support it.");
        #endif
            
    }

    typedef __CL_FUNCTION_TYPE<Event (
        EnqueueArgs&,
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13,
		T14)> type_;

	Event operator()(
		EnqueueArgs& enqueueArgs,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const Event &event,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const VECTOR_CLASS<Event>& waitEvents,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13,
		T14 arg14)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13,
			arg14);
	}

    operator type_ ()    
    {
        return type_(*this);
    }

};

template<
	typename T0,
	typename T1,
	typename T2,
	typename T3,
	typename T4,
	typename T5,
	typename T6,
	typename T7,
	typename T8,
	typename T9,
	typename T10,
	typename T11,
	typename T12,
	typename T13>
struct functionImplementation_
<	T0,
	T1,
	T2,
	T3,
	T4,
	T5,
	T6,
	T7,
	T8,
	T9,
	T10,
	T11,
	T12,
	T13,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType>
{
	typedef detail::KernelFunctorGlobal<
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType> FunctorType;

    FunctorType functor_;

    functionImplementation_(const FunctorType &functor) :
        functor_(functor)
    {
    
        #if (defined(_WIN32) && (_VARIADIC_MAX < 14))
        // Fail variadic expansion for dev11
        static_assert(0, "Visual Studio has a hard limit of argument count for a std::function expansion. Please define _VARIADIC_MAX to be 10. If you need more arguments than that VC12 and below cannot support it.");
        #endif
            
    }

    typedef __CL_FUNCTION_TYPE<Event (
        EnqueueArgs&,
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		T13)> type_;

	Event operator()(
		EnqueueArgs& enqueueArgs,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const Event &event,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const VECTOR_CLASS<Event>& waitEvents,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12,
		T13 arg13)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12,
			arg13);
	}

    operator type_ ()    
    {
        return type_(*this);
    }

};

template<
	typename T0,
	typename T1,
	typename T2,
	typename T3,
	typename T4,
	typename T5,
	typename T6,
	typename T7,
	typename T8,
	typename T9,
	typename T10,
	typename T11,
	typename T12>
struct functionImplementation_
<	T0,
	T1,
	T2,
	T3,
	T4,
	T5,
	T6,
	T7,
	T8,
	T9,
	T10,
	T11,
	T12,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType>
{
	typedef detail::KernelFunctorGlobal<
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType> FunctorType;

    FunctorType functor_;

    functionImplementation_(const FunctorType &functor) :
        functor_(functor)
    {
    
        #if (defined(_WIN32) && (_VARIADIC_MAX < 13))
        // Fail variadic expansion for dev11
        static_assert(0, "Visual Studio has a hard limit of argument count for a std::function expansion. Please define _VARIADIC_MAX to be 10. If you need more arguments than that VC12 and below cannot support it.");
        #endif
            
    }

    typedef __CL_FUNCTION_TYPE<Event (
        EnqueueArgs&,
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		T12)> type_;

	Event operator()(
		EnqueueArgs& enqueueArgs,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const Event &event,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const VECTOR_CLASS<Event>& waitEvents,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11,
		T12 arg12)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11,
			arg12);
	}

    operator type_ ()    
    {
        return type_(*this);
    }

};

template<
	typename T0,
	typename T1,
	typename T2,
	typename T3,
	typename T4,
	typename T5,
	typename T6,
	typename T7,
	typename T8,
	typename T9,
	typename T10,
	typename T11>
struct functionImplementation_
<	T0,
	T1,
	T2,
	T3,
	T4,
	T5,
	T6,
	T7,
	T8,
	T9,
	T10,
	T11,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType>
{
	typedef detail::KernelFunctorGlobal<
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType> FunctorType;

    FunctorType functor_;

    functionImplementation_(const FunctorType &functor) :
        functor_(functor)
    {
    
        #if (defined(_WIN32) && (_VARIADIC_MAX < 12))
        // Fail variadic expansion for dev11
        static_assert(0, "Visual Studio has a hard limit of argument count for a std::function expansion. Please define _VARIADIC_MAX to be 10. If you need more arguments than that VC12 and below cannot support it.");
        #endif
            
    }

    typedef __CL_FUNCTION_TYPE<Event (
        EnqueueArgs&,
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		T11)> type_;

	Event operator()(
		EnqueueArgs& enqueueArgs,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const Event &event,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const VECTOR_CLASS<Event>& waitEvents,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10,
		T11 arg11)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10,
			arg11);
	}

    operator type_ ()    
    {
        return type_(*this);
    }

};

template<
	typename T0,
	typename T1,
	typename T2,
	typename T3,
	typename T4,
	typename T5,
	typename T6,
	typename T7,
	typename T8,
	typename T9,
	typename T10>
struct functionImplementation_
<	T0,
	T1,
	T2,
	T3,
	T4,
	T5,
	T6,
	T7,
	T8,
	T9,
	T10,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType>
{
	typedef detail::KernelFunctorGlobal<
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType> FunctorType;

    FunctorType functor_;

    functionImplementation_(const FunctorType &functor) :
        functor_(functor)
    {
    
        #if (defined(_WIN32) && (_VARIADIC_MAX < 11))
        // Fail variadic expansion for dev11
        static_assert(0, "Visual Studio has a hard limit of argument count for a std::function expansion. Please define _VARIADIC_MAX to be 10. If you need more arguments than that VC12 and below cannot support it.");
        #endif
            
    }

    typedef __CL_FUNCTION_TYPE<Event (
        EnqueueArgs&,
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		T10)> type_;

	Event operator()(
		EnqueueArgs& enqueueArgs,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const Event &event,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const VECTOR_CLASS<Event>& waitEvents,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9,
		T10 arg10)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9,
			arg10);
	}

    operator type_ ()    
    {
        return type_(*this);
    }

};

template<
	typename T0,
	typename T1,
	typename T2,
	typename T3,
	typename T4,
	typename T5,
	typename T6,
	typename T7,
	typename T8,
	typename T9>
struct functionImplementation_
<	T0,
	T1,
	T2,
	T3,
	T4,
	T5,
	T6,
	T7,
	T8,
	T9,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType>
{
	typedef detail::KernelFunctorGlobal<
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType> FunctorType;

    FunctorType functor_;

    functionImplementation_(const FunctorType &functor) :
        functor_(functor)
    {
    
        #if (defined(_WIN32) && (_VARIADIC_MAX < 10))
        // Fail variadic expansion for dev11
        static_assert(0, "Visual Studio has a hard limit of argument count for a std::function expansion. Please define _VARIADIC_MAX to be 10. If you need more arguments than that VC12 and below cannot support it.");
        #endif
            
    }

    typedef __CL_FUNCTION_TYPE<Event (
        EnqueueArgs&,
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		T9)> type_;

	Event operator()(
		EnqueueArgs& enqueueArgs,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const Event &event,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const VECTOR_CLASS<Event>& waitEvents,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8,
		T9 arg9)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8,
			arg9);
	}

    operator type_ ()    
    {
        return type_(*this);
    }

};

template<
	typename T0,
	typename T1,
	typename T2,
	typename T3,
	typename T4,
	typename T5,
	typename T6,
	typename T7,
	typename T8>
struct functionImplementation_
<	T0,
	T1,
	T2,
	T3,
	T4,
	T5,
	T6,
	T7,
	T8,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType>
{
	typedef detail::KernelFunctorGlobal<
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType> FunctorType;

    FunctorType functor_;

    functionImplementation_(const FunctorType &functor) :
        functor_(functor)
    {
    
        #if (defined(_WIN32) && (_VARIADIC_MAX < 9))
        // Fail variadic expansion for dev11
        static_assert(0, "Visual Studio has a hard limit of argument count for a std::function expansion. Please define _VARIADIC_MAX to be 10. If you need more arguments than that VC12 and below cannot support it.");
        #endif
            
    }

    typedef __CL_FUNCTION_TYPE<Event (
        EnqueueArgs&,
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		T8)> type_;

	Event operator()(
		EnqueueArgs& enqueueArgs,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const Event &event,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const VECTOR_CLASS<Event>& waitEvents,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7,
		T8 arg8)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7,
			arg8);
	}

    operator type_ ()    
    {
        return type_(*this);
    }

};

template<
	typename T0,
	typename T1,
	typename T2,
	typename T3,
	typename T4,
	typename T5,
	typename T6,
	typename T7>
struct functionImplementation_
<	T0,
	T1,
	T2,
	T3,
	T4,
	T5,
	T6,
	T7,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType>
{
	typedef detail::KernelFunctorGlobal<
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType> FunctorType;

    FunctorType functor_;

    functionImplementation_(const FunctorType &functor) :
        functor_(functor)
    {
    
        #if (defined(_WIN32) && (_VARIADIC_MAX < 8))
        // Fail variadic expansion for dev11
        static_assert(0, "Visual Studio has a hard limit of argument count for a std::function expansion. Please define _VARIADIC_MAX to be 10. If you need more arguments than that VC12 and below cannot support it.");
        #endif
            
    }

    typedef __CL_FUNCTION_TYPE<Event (
        EnqueueArgs&,
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		T7)> type_;

	Event operator()(
		EnqueueArgs& enqueueArgs,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const Event &event,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const VECTOR_CLASS<Event>& waitEvents,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6,
		T7 arg7)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6,
			arg7);
	}

    operator type_ ()    
    {
        return type_(*this);
    }

};

template<
	typename T0,
	typename T1,
	typename T2,
	typename T3,
	typename T4,
	typename T5,
	typename T6>
struct functionImplementation_
<	T0,
	T1,
	T2,
	T3,
	T4,
	T5,
	T6,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType>
{
	typedef detail::KernelFunctorGlobal<
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType> FunctorType;

    FunctorType functor_;

    functionImplementation_(const FunctorType &functor) :
        functor_(functor)
    {
    
        #if (defined(_WIN32) && (_VARIADIC_MAX < 7))
        // Fail variadic expansion for dev11
        static_assert(0, "Visual Studio has a hard limit of argument count for a std::function expansion. Please define _VARIADIC_MAX to be 10. If you need more arguments than that VC12 and below cannot support it.");
        #endif
            
    }

    typedef __CL_FUNCTION_TYPE<Event (
        EnqueueArgs&,
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		T6)> type_;

	Event operator()(
		EnqueueArgs& enqueueArgs,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const Event &event,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const VECTOR_CLASS<Event>& waitEvents,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5,
		T6 arg6)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5,
			arg6);
	}

    operator type_ ()    
    {
        return type_(*this);
    }

};

template<
	typename T0,
	typename T1,
	typename T2,
	typename T3,
	typename T4,
	typename T5>
struct functionImplementation_
<	T0,
	T1,
	T2,
	T3,
	T4,
	T5,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType>
{
	typedef detail::KernelFunctorGlobal<
		T0,
		T1,
		T2,
		T3,
		T4,
		T5,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType> FunctorType;

    FunctorType functor_;

    functionImplementation_(const FunctorType &functor) :
        functor_(functor)
    {
    
        #if (defined(_WIN32) && (_VARIADIC_MAX < 6))
        // Fail variadic expansion for dev11
        static_assert(0, "Visual Studio has a hard limit of argument count for a std::function expansion. Please define _VARIADIC_MAX to be 10. If you need more arguments than that VC12 and below cannot support it.");
        #endif
            
    }

    typedef __CL_FUNCTION_TYPE<Event (
        EnqueueArgs&,
		T0,
		T1,
		T2,
		T3,
		T4,
		T5)> type_;

	Event operator()(
		EnqueueArgs& enqueueArgs,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const Event &event,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const VECTOR_CLASS<Event>& waitEvents,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4,
		T5 arg5)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4,
			arg5);
	}

    operator type_ ()    
    {
        return type_(*this);
    }

};

template<
	typename T0,
	typename T1,
	typename T2,
	typename T3,
	typename T4>
struct functionImplementation_
<	T0,
	T1,
	T2,
	T3,
	T4,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType>
{
	typedef detail::KernelFunctorGlobal<
		T0,
		T1,
		T2,
		T3,
		T4,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType> FunctorType;

    FunctorType functor_;

    functionImplementation_(const FunctorType &functor) :
        functor_(functor)
    {
    
        #if (defined(_WIN32) && (_VARIADIC_MAX < 5))
        // Fail variadic expansion for dev11
        static_assert(0, "Visual Studio has a hard limit of argument count for a std::function expansion. Please define _VARIADIC_MAX to be 10. If you need more arguments than that VC12 and below cannot support it.");
        #endif
            
    }

    typedef __CL_FUNCTION_TYPE<Event (
        EnqueueArgs&,
		T0,
		T1,
		T2,
		T3,
		T4)> type_;

	Event operator()(
		EnqueueArgs& enqueueArgs,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const Event &event,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const VECTOR_CLASS<Event>& waitEvents,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3,
		T4 arg4)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3,
			arg4);
	}

    operator type_ ()    
    {
        return type_(*this);
    }

};

template<
	typename T0,
	typename T1,
	typename T2,
	typename T3>
struct functionImplementation_
<	T0,
	T1,
	T2,
	T3,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType>
{
	typedef detail::KernelFunctorGlobal<
		T0,
		T1,
		T2,
		T3,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType> FunctorType;

    FunctorType functor_;

    functionImplementation_(const FunctorType &functor) :
        functor_(functor)
    {
    
        #if (defined(_WIN32) && (_VARIADIC_MAX < 4))
        // Fail variadic expansion for dev11
        static_assert(0, "Visual Studio has a hard limit of argument count for a std::function expansion. Please define _VARIADIC_MAX to be 10. If you need more arguments than that VC12 and below cannot support it.");
        #endif
            
    }

    typedef __CL_FUNCTION_TYPE<Event (
        EnqueueArgs&,
		T0,
		T1,
		T2,
		T3)> type_;

	Event operator()(
		EnqueueArgs& enqueueArgs,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const Event &event,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const VECTOR_CLASS<Event>& waitEvents,
		T0 arg0,
		T1 arg1,
		T2 arg2,
		T3 arg3)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2,
			arg3);
	}

    operator type_ ()    
    {
        return type_(*this);
    }

};

template<
	typename T0,
	typename T1,
	typename T2>
struct functionImplementation_
<	T0,
	T1,
	T2,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType>
{
	typedef detail::KernelFunctorGlobal<
		T0,
		T1,
		T2,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType> FunctorType;

    FunctorType functor_;

    functionImplementation_(const FunctorType &functor) :
        functor_(functor)
    {
    
        #if (defined(_WIN32) && (_VARIADIC_MAX < 3))
        // Fail variadic expansion for dev11
        static_assert(0, "Visual Studio has a hard limit of argument count for a std::function expansion. Please define _VARIADIC_MAX to be 10. If you need more arguments than that VC12 and below cannot support it.");
        #endif
            
    }

    typedef __CL_FUNCTION_TYPE<Event (
        EnqueueArgs&,
		T0,
		T1,
		T2)> type_;

	Event operator()(
		EnqueueArgs& enqueueArgs,
		T0 arg0,
		T1 arg1,
		T2 arg2)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const Event &event,
		T0 arg0,
		T1 arg1,
		T2 arg2)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const VECTOR_CLASS<Event>& waitEvents,
		T0 arg0,
		T1 arg1,
		T2 arg2)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1,
			arg2);
	}

    operator type_ ()    
    {
        return type_(*this);
    }

};

template<
	typename T0,
	typename T1>
struct functionImplementation_
<	T0,
	T1,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType>
{
	typedef detail::KernelFunctorGlobal<
		T0,
		T1,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType> FunctorType;

    FunctorType functor_;

    functionImplementation_(const FunctorType &functor) :
        functor_(functor)
    {
    
        #if (defined(_WIN32) && (_VARIADIC_MAX < 2))
        // Fail variadic expansion for dev11
        static_assert(0, "Visual Studio has a hard limit of argument count for a std::function expansion. Please define _VARIADIC_MAX to be 10. If you need more arguments than that VC12 and below cannot support it.");
        #endif
            
    }

    typedef __CL_FUNCTION_TYPE<Event (
        EnqueueArgs&,
		T0,
		T1)> type_;

	Event operator()(
		EnqueueArgs& enqueueArgs,
		T0 arg0,
		T1 arg1)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const Event &event,
		T0 arg0,
		T1 arg1)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const VECTOR_CLASS<Event>& waitEvents,
		T0 arg0,
		T1 arg1)
	{
		return functor_(
			enqueueArgs,
			arg0,
			arg1);
	}

    operator type_ ()    
    {
        return type_(*this);
    }

};

template<
	typename T0>
struct functionImplementation_
<	T0,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType,
	NullType>
{
	typedef detail::KernelFunctorGlobal<
		T0,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType,
		NullType> FunctorType;

    FunctorType functor_;

    functionImplementation_(const FunctorType &functor) :
        functor_(functor)
    {
    
        #if (defined(_WIN32) && (_VARIADIC_MAX < 1))
        // Fail variadic expansion for dev11
        static_assert(0, "Visual Studio has a hard limit of argument count for a std::function expansion. Please define _VARIADIC_MAX to be 10. If you need more arguments than that VC12 and below cannot support it.");
        #endif
            
    }

    typedef __CL_FUNCTION_TYPE<Event (
        EnqueueArgs&,
		T0)> type_;

	Event operator()(
		EnqueueArgs& enqueueArgs,
		T0 arg0)
	{
		return functor_(
			enqueueArgs,
			arg0);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const Event &event,
		T0 arg0)
	{
		return functor_(
			enqueueArgs,
			arg0);
	}
	Event operator()(
		EnqueueArgs& enqueueArgs, 
		const VECTOR_CLASS<Event>& waitEvents,
		T0 arg0)
	{
		return functor_(
			enqueueArgs,
			arg0);
	}

    operator type_ ()    
    {
        return type_(*this);
    }

};





} // namespace detail

//----------------------------------------------------------------------------------------------

template <
   typename T0,   typename T1 = detail::NullType,   typename T2 = detail::NullType,
   typename T3 = detail::NullType,   typename T4 = detail::NullType,
   typename T5 = detail::NullType,   typename T6 = detail::NullType,
   typename T7 = detail::NullType,   typename T8 = detail::NullType,
   typename T9 = detail::NullType,   typename T10 = detail::NullType,
   typename T11 = detail::NullType,   typename T12 = detail::NullType,
   typename T13 = detail::NullType,   typename T14 = detail::NullType,
   typename T15 = detail::NullType,   typename T16 = detail::NullType,
   typename T17 = detail::NullType,   typename T18 = detail::NullType,
   typename T19 = detail::NullType,   typename T20 = detail::NullType,
   typename T21 = detail::NullType,   typename T22 = detail::NullType,
   typename T23 = detail::NullType,   typename T24 = detail::NullType,
   typename T25 = detail::NullType,   typename T26 = detail::NullType,
   typename T27 = detail::NullType,   typename T28 = detail::NullType,
   typename T29 = detail::NullType,   typename T30 = detail::NullType,
   typename T31 = detail::NullType
>
struct make_kernel :
    public detail::functionImplementation_<
               T0,   T1,   T2,   T3,
               T4,   T5,   T6,   T7,
               T8,   T9,   T10,   T11,
               T12,   T13,   T14,   T15,
               T16,   T17,   T18,   T19,
               T20,   T21,   T22,   T23,
               T24,   T25,   T26,   T27,
               T28,   T29,   T30,   T31
    >
{
public:
	typedef detail::KernelFunctorGlobal<             
		       T0,   T1,   T2,   T3,
               T4,   T5,   T6,   T7,
               T8,   T9,   T10,   T11,
               T12,   T13,   T14,   T15,
               T16,   T17,   T18,   T19,
               T20,   T21,   T22,   T23,
               T24,   T25,   T26,   T27,
               T28,   T29,   T30,   T31
    > FunctorType;

    make_kernel(
        const Program& program,
        const STRING_CLASS name,
        cl_int * err = NULL) :
           detail::functionImplementation_<
                    T0,   T1,   T2,   T3,
                       T4,   T5,   T6,   T7,
                       T8,   T9,   T10,   T11,
                       T12,   T13,   T14,   T15,
                       T16,   T17,   T18,   T19,
                       T20,   T21,   T22,   T23,
                       T24,   T25,   T26,   T27,
                       T28,   T29,   T30,   T31
           >(
            FunctorType(program, name, err)) 
    {}

    make_kernel(
        const Kernel kernel) :
           detail::functionImplementation_<
                    T0,   T1,   T2,   T3,
                       T4,   T5,   T6,   T7,
                       T8,   T9,   T10,   T11,
                       T12,   T13,   T14,   T15,
                       T16,   T17,   T18,   T19,
                       T20,   T21,   T22,   T23,
                       T24,   T25,   T26,   T27,
                       T28,   T29,   T30,   T31
           >(
            FunctorType(kernel)) 
    {}    
};

#endif // #if (defined(CL_USE_CPP_FUNCTORS)


//----------------------------------------------------------------------------------------------------------------------

#undef __ERR_STR
#if !defined(__CL_USER_OVERRIDE_ERROR_STRINGS)
#undef __GET_DEVICE_INFO_ERR
#undef __GET_PLATFORM_INFO_ERR
#undef __GET_DEVICE_IDS_ERR
#undef __GET_CONTEXT_INFO_ERR
#undef __GET_EVENT_INFO_ERR
#undef __GET_EVENT_PROFILE_INFO_ERR
#undef __GET_MEM_OBJECT_INFO_ERR
#undef __GET_IMAGE_INFO_ERR
#undef __GET_SAMPLER_INFO_ERR
#undef __GET_KERNEL_INFO_ERR
#undef __GET_KERNEL_ARG_INFO_ERR
#undef __GET_KERNEL_WORK_GROUP_INFO_ERR
#undef __GET_PROGRAM_INFO_ERR
#undef __GET_PROGRAM_BUILD_INFO_ERR
#undef __GET_COMMAND_QUEUE_INFO_ERR

#undef __CREATE_CONTEXT_ERR
#undef __CREATE_CONTEXT_FROM_TYPE_ERR
#undef __GET_SUPPORTED_IMAGE_FORMATS_ERR

#undef __CREATE_BUFFER_ERR
#undef __CREATE_SUBBUFFER_ERR
#undef __CREATE_IMAGE2D_ERR
#undef __CREATE_IMAGE3D_ERR
#undef __CREATE_SAMPLER_ERR
#undef __SET_MEM_OBJECT_DESTRUCTOR_CALLBACK_ERR

#undef __CREATE_USER_EVENT_ERR
#undef __SET_USER_EVENT_STATUS_ERR
#undef __SET_EVENT_CALLBACK_ERR
#undef __SET_PRINTF_CALLBACK_ERR

#undef __WAIT_FOR_EVENTS_ERR

#undef __CREATE_KERNEL_ERR
#undef __SET_KERNEL_ARGS_ERR
#undef __CREATE_PROGRAM_WITH_SOURCE_ERR
#undef __CREATE_PROGRAM_WITH_BINARY_ERR
#undef __CREATE_PROGRAM_WITH_BUILT_IN_KERNELS_ERR
#undef __BUILD_PROGRAM_ERR
#undef __CREATE_KERNELS_IN_PROGRAM_ERR

#undef __CREATE_COMMAND_QUEUE_ERR
#undef __SET_COMMAND_QUEUE_PROPERTY_ERR
#undef __ENQUEUE_READ_BUFFER_ERR
#undef __ENQUEUE_WRITE_BUFFER_ERR
#undef __ENQUEUE_READ_BUFFER_RECT_ERR
#undef __ENQUEUE_WRITE_BUFFER_RECT_ERR
#undef __ENQEUE_COPY_BUFFER_ERR
#undef __ENQEUE_COPY_BUFFER_RECT_ERR
#undef __ENQUEUE_READ_IMAGE_ERR
#undef __ENQUEUE_WRITE_IMAGE_ERR
#undef __ENQUEUE_COPY_IMAGE_ERR
#undef __ENQUEUE_COPY_IMAGE_TO_BUFFER_ERR
#undef __ENQUEUE_COPY_BUFFER_TO_IMAGE_ERR
#undef __ENQUEUE_MAP_BUFFER_ERR
#undef __ENQUEUE_MAP_IMAGE_ERR
#undef __ENQUEUE_UNMAP_MEM_OBJECT_ERR
#undef __ENQUEUE_NDRANGE_KERNEL_ERR
#undef __ENQUEUE_TASK_ERR
#undef __ENQUEUE_NATIVE_KERNEL

#undef __UNLOAD_COMPILER_ERR
#endif //__CL_USER_OVERRIDE_ERROR_STRINGS

#undef __GET_INFO_HELPER_WITH_RETAIN

#undef __CL_FUNCTION_TYPE

// Extensions
/**
 * Deprecated APIs for 1.2
 */
#if defined(CL_VERSION_1_1)
#undef __INIT_CL_EXT_FCN_PTR
#endif // #if defined(CL_VERSION_1_1)
#undef __CREATE_SUB_DEVICES

#if defined(USE_CL_DEVICE_FISSION)
#undef __PARAM_NAME_DEVICE_FISSION
#endif // USE_CL_DEVICE_FISSION

#undef __DEFAULT_NOT_INITIALIZED 
#undef __DEFAULT_BEING_INITIALIZED 
#undef __DEFAULT_INITIALIZED

} // namespace cl

#ifdef _WIN32
#pragma pop_macro("max")
#endif // _WIN32

#endif // CL_HPP_
