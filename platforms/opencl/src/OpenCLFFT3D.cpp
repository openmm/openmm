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

#include "OpenCLFFT3D.h"
#include "OpenCLExpressionUtilities.h"
#include "OpenCLKernelSources.h"
#include "../src/SimTKUtilities/SimTKOpenMMRealType.h"
#include <map>
#include <sstream>
#include <string>

using namespace OpenMM;
using namespace std;

OpenCLFFT3D::OpenCLFFT3D(OpenCLContext& context, int xsize, int ysize, int zsize) : context(context), xsize(xsize), ysize(ysize), zsize(zsize) {
    xkernel = createKernel(xsize, ysize, zsize, ysize*zsize, zsize, 1);
    ykernel = createKernel(ysize, zsize, xsize, zsize, 1, ysize*zsize);
    zkernel = createKernel(zsize, xsize, ysize, 1, ysize*zsize, zsize);
}

void OpenCLFFT3D::execFFT(OpenCLArray<mm_float2>& data, bool forward) {
    xkernel.setArg<cl::Buffer>(0, data.getDeviceBuffer());
    xkernel.setArg<cl_float>(1, forward ? 1.0f : -1.0f);
    context.executeKernel(xkernel, xsize*ysize*zsize, xsize);
    ykernel.setArg<cl::Buffer>(0, data.getDeviceBuffer());
    ykernel.setArg<cl_float>(1, forward ? 1.0f : -1.0f);
    context.executeKernel(ykernel, xsize*ysize*zsize, ysize);
    zkernel.setArg<cl::Buffer>(0, data.getDeviceBuffer());
    zkernel.setArg<cl_float>(1, forward ? 1.0f : -1.0f);
    context.executeKernel(zkernel, xsize*ysize*zsize, zsize);
}

int OpenCLFFT3D::findLegalDimension(int minimum) {
    if (minimum < 1)
        return 1;
    while (true) {
        // Attempt to factor the current value.

        int unfactored = minimum;
        for (int factor = 2; factor < 6; factor++) {
            while (unfactored > 1 && unfactored%factor == 0)
                unfactored /= factor;
        }
        if (unfactored == 1)
            return minimum;
        minimum++;
    }
}

cl::Kernel OpenCLFFT3D::createKernel(int xsize, int ysize, int zsize, int xmult, int ymult, int zmult) {
    stringstream source;
    int unfactored = xsize;
    int stage = 0;
    int L = xsize;
    int m = 1;

    // Factor xsize, generating an appropriate block of code for each factor.

    while (unfactored > 1) {
        int input = stage%2;
        int output = 1-input;
        source<<"{\n";
        if (unfactored%5 == 0) {
            L = L/5;
            source<<"// Pass "<<(stage+1)<<" (radix 5)\n";
            source<<"int j = i/"<<m<<";\n";
            source<<"float2 c0 = data"<<input<<"[i];\n";
            source<<"float2 c1 = data"<<input<<"[i+"<<(L*m)<<"];\n";
            source<<"float2 c2 = data"<<input<<"[i+"<<(2*L*m)<<"];\n";
            source<<"float2 c3 = data"<<input<<"[i+"<<(3*L*m)<<"];\n";
            source<<"float2 c4 = data"<<input<<"[i+"<<(4*L*m)<<"];\n";
            source<<"float2 d0 = c1+c4;\n";
            source<<"float2 d1 = c2+c3;\n";
            source<<"float2 d2 = "<<OpenCLExpressionUtilities::doubleToString(sin(0.4*M_PI))<<"*(c1-c4);\n";
            source<<"float2 d3 = "<<OpenCLExpressionUtilities::doubleToString(sin(0.4*M_PI))<<"*(c2-c3);\n";
            source<<"float2 d4 = d0+d1;\n";
            source<<"float2 d5 = "<<OpenCLExpressionUtilities::doubleToString(0.25*sqrt(5.0))<<"*(d0-d1);\n";
            source<<"float2 d6 = c0-0.25f*d4;\n";
            source<<"float2 d7 = d6+d5;\n";
            source<<"float2 d8 = d6-d5;\n";
            string coeff = OpenCLExpressionUtilities::doubleToString(sin(0.2*M_PI)/sin(0.4*M_PI));
            source<<"float2 d9 = sign*(float2) (d2.y+"<<coeff<<"*d3.y, -d2.x-"<<coeff<<"*d3.x);\n";
            source<<"float2 d10 = sign*(float2) ("<<coeff<<"*d2.y-d3.y, d3.x-"<<coeff<<"*d2.x);\n";
            source<<"data"<<output<<"[i+4*j*"<<m<<"] = c0+d4;\n";
            source<<"data"<<output<<"[i+(4*j+1)*"<<m<<"] = multiplyComplex(w[j*"<<xsize<<"/"<<(5*L)<<"], d7+d9);\n";
            source<<"data"<<output<<"[i+(4*j+2)*"<<m<<"] = multiplyComplex(w[j*"<<(2*xsize)<<"/"<<(5*L)<<"], d8+d10);\n";
            source<<"data"<<output<<"[i+(4*j+3)*"<<m<<"] = multiplyComplex(w[j*"<<(3*xsize)<<"/"<<(5*L)<<"], d8-d10);\n";
            source<<"data"<<output<<"[i+(4*j+4)*"<<m<<"] = multiplyComplex(w[j*"<<(4*xsize)<<"/"<<(5*L)<<"], d7-d9);\n";
            m = m*5;
            unfactored /= 5;
        }
        else if (unfactored%4 == 0) {
            L = L/4;
            source<<"// Pass "<<(stage+1)<<" (radix 4)\n";
            source<<"int j = i/"<<m<<";\n";
            source<<"float2 c0 = data"<<input<<"[i];\n";
            source<<"float2 c1 = data"<<input<<"[i+"<<(L*m)<<"];\n";
            source<<"float2 c2 = data"<<input<<"[i+"<<(2*L*m)<<"];\n";
            source<<"float2 c3 = data"<<input<<"[i+"<<(3*L*m)<<"];\n";
            source<<"float2 d0 = c0+c2;\n";
            source<<"float2 d1 = c0-c2;\n";
            source<<"float2 d2 = c1+c3;\n";
            source<<"float2 d3 = sign*(float2) (c1.y-c3.y, c3.x-c1.x);\n";
            source<<"data"<<output<<"[i+3*j*"<<m<<"] = d0+d2;\n";
            source<<"data"<<output<<"[i+(3*j+1)*"<<m<<"] = multiplyComplex(w[j*"<<xsize<<"/"<<(4*L)<<"], d1+d3);\n";
            source<<"data"<<output<<"[i+(3*j+2)*"<<m<<"] = multiplyComplex(w[j*"<<(2*xsize)<<"/"<<(4*L)<<"], d0-d2);\n";
            source<<"data"<<output<<"[i+(3*j+3)*"<<m<<"] = multiplyComplex(w[j*"<<(3*xsize)<<"/"<<(4*L)<<"], d1-d3);\n";
            m = m*4;
            unfactored /= 4;
        }
        else if (unfactored%3 == 0) {
            L = L/3;
            source<<"// Pass "<<(stage+1)<<" (radix 3)\n";
            source<<"int j = i/"<<m<<";\n";
            source<<"float2 c0 = data"<<input<<"[i];\n";
            source<<"float2 c1 = data"<<input<<"[i+"<<(L*m)<<"];\n";
            source<<"float2 c2 = data"<<input<<"[i+"<<(2*L*m)<<"];\n";
            source<<"float2 d0 = c1+c2;\n";
            source<<"float2 d1 = c0-0.5f*d0;\n";
            source<<"float2 d2 = sign*"<<OpenCLExpressionUtilities::doubleToString(sin(M_PI/3.0))<<"*(float2) (c1.y-c2.y, c2.x-c1.x);\n";
            source<<"data"<<output<<"[i+2*j*"<<m<<"] = c0+d0;\n";
            source<<"data"<<output<<"[i+(2*j+1)*"<<m<<"] = multiplyComplex(w[j*"<<xsize<<"/"<<(3*L)<<"], d1+d2);\n";
            source<<"data"<<output<<"[i+(2*j+2)*"<<m<<"] = multiplyComplex(w[j*"<<(2*xsize)<<"/"<<(3*L)<<"], d1-d2);\n";
            m = m*3;
            unfactored /= 3;
        }
        else if (unfactored%2 == 0) {
            L = L/2;
            source<<"// Pass "<<(stage+1)<<" (radix 2)\n";
            source<<"int j = i/"<<m<<";\n";
            source<<"float2 c0 = data"<<input<<"[i];\n";
            source<<"float2 c1 = data"<<input<<"[i+"<<(L*m)<<"];\n";
            source<<"data"<<output<<"[i+j*"<<m<<"] = c0+c1;\n";
            source<<"data"<<output<<"[i+(j+1)*"<<m<<"] = multiplyComplex(w[j*"<<xsize<<"/"<<(2*L)<<"], c0-c1);\n";
            m = m*2;
            unfactored /= 2;
        }
        else
            throw OpenMMException("Illegal size for FFT: "+OpenCLExpressionUtilities::intToString(xsize));
        source<<"barrier(CLK_LOCAL_MEM_FENCE);\n";
        source<<"}\n";
        ++stage;
    }

    // Create the kernel.

    source<<"matrix[element] = data"<<(stage%2)<<"[i];";
    map<string, string> replacements;
    replacements["XSIZE"] = OpenCLExpressionUtilities::intToString(xsize);
    replacements["YSIZE"] = OpenCLExpressionUtilities::intToString(ysize);
    replacements["ZSIZE"] = OpenCLExpressionUtilities::intToString(zsize);
    replacements["XMULT"] = OpenCLExpressionUtilities::intToString(xmult);
    replacements["YMULT"] = OpenCLExpressionUtilities::intToString(ymult);
    replacements["ZMULT"] = OpenCLExpressionUtilities::intToString(zmult);
    replacements["M_PI"] = OpenCLExpressionUtilities::doubleToString(M_PI);
    replacements["COMPUTE_FFT"] = source.str();
    cl::Program program = context.createProgram(context.replaceStrings(OpenCLKernelSources::fft, replacements));
    cl::Kernel kernel(program, "execFFT");
    kernel.setArg(2, xsize*sizeof(mm_float2), NULL);
    kernel.setArg(3, xsize*sizeof(mm_float2), NULL);
    kernel.setArg(4, xsize*sizeof(mm_float2), NULL);
    return kernel;
}
