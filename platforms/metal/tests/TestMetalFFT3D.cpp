/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 *                                                                            *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Metal Platform code:                                                       *
 * Portions copyright (c) 2026 Chun-Chi Hung.                                 *
 * Authors: Chun-Chi Hung                                                     *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the               *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program. If not, see <http://www.gnu.org/licenses/>.       *
 * -------------------------------------------------------------------------- */

#include "MetalContext.h"
#include "MetalFFT3D.h"
#include "openmm/System.h"
#include "openmm/OpenMMException.h"
#include "openmm/common/ComputeArray.h"
#include "openmm/common/ComputeVectorTypes.h"
#include "openmm/internal/AssertionUtilities.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

using namespace OpenMM;
using namespace std;

void expectException(const function<void()>& operation) {
    try {
        operation();
    }
    catch (const OpenMMException&) {
        return;
    }
    throw OpenMMException("Expected an invalid Metal FFT request to fail");
}

vector<complex<double>> referenceFFT(const vector<complex<double>>& input, int nx, int ny, int nz) {
    vector<complex<double>> output(input.size());
    const double scale = -2*acos(-1.0);
    for (int kx = 0; kx < nx; kx++)
        for (int ky = 0; ky < ny; ky++)
            for (int kz = 0; kz < nz; kz++) {
                complex<double> sum(0, 0);
                for (int x = 0; x < nx; x++)
                    for (int y = 0; y < ny; y++)
                        for (int z = 0; z < nz; z++) {
                            const double angle = scale*(double(kx*x)/nx+double(ky*y)/ny+double(kz*z)/nz);
                            sum += input[(x*ny+y)*nz+z]*complex<double>(cos(angle), sin(angle));
                        }
                output[(kx*ny+ky)*nz+kz] = sum;
            }
    return output;
}

void checkValue(double expected, double actual) {
    ASSERT(isfinite(actual));
    // Diagnostic single-precision threshold, not a settled platform-wide accuracy policy.
    ASSERT_EQUAL_TOL(expected, actual, 2e-4);
}

void testTransform(MetalContext& context, int nx, int ny, int nz, bool realToComplex) {
    const int count = nx*ny*nz;
    ComputeArray input, output;
    input.initialize<float>(context, 2*count, "fftInput");
    output.initialize<float>(context, 2*count, "fftOutput");
    vector<float> values(2*count, 0);
    vector<complex<double>> source(count);
    for (int i = 0; i < count; i++) {
        float real = 0.125f+sin(0.17*i)+float(i%7)/11;
        float imag = realToComplex ? 0 : cos(0.31*i)-float(i%5)/13;
        source[i] = complex<double>(real, imag);
        if (realToComplex)
            values[i] = real;
        else {
            values[2*i] = real;
            values[2*i+1] = imag;
        }
    }
    if (realToComplex && count == 1)
        values[1] = 37; // Padding is not the imaginary component of real input.
    vector<complex<double>> expected = referenceFFT(source, nx, ny, nz);
    MetalFFT3D fft(context, nx, ny, nz, realToComplex);
    input.upload(values);
    fft.execFFT(input, output);
    vector<float> actual;
    output.download(actual);
    const int outputZ = realToComplex ? nz/2+1 : nz;
    for (int x = 0; x < nx; x++)
        for (int y = 0; y < ny; y++)
            for (int z = 0; z < outputZ; z++) {
                const int index = (x*ny+y)*outputZ+z;
                const complex<double> reference = expected[(x*ny+y)*nz+z];
                checkValue(reference.real(), actual[2*index]);
                checkValue(reference.imag(), actual[2*index+1]);
            }
    fft.execFFT(output, input, false);
    input.download(actual);
    for (int i = 0; i < (realToComplex ? count : 2*count); i++)
        checkValue(values[i], actual[i]/count);

    // The same plan must use current allocations after resize and a queue change.
    context.setCurrentQueue(context.createQueue());
    input.resize(2*count);
    output.resize(2*count);
    input.upload(values);
    fft.execFFT(input, output);
    fft.execFFT(output, input, false);
    input.download(actual);
    for (int i = 0; i < (realToComplex ? count : 2*count); i++)
        checkValue(values[i], actual[i]/count);
    context.restoreDefaultQueue();
}

void testBatchedPlanReuse(MetalContext& context) {
    const int nx = 3, ny = 5, nz = 8, count = nx*ny*nz;
    MetalFFT3D fft(context, nx, ny, nz, false);
    ComputeArray original, input, output, saved;
    original.initialize<mm_float2>(context, count, "fftOriginal");
    input.initialize<mm_float2>(context, count, "fftWork");
    output.initialize<mm_float2>(context, count, "fftSpectrum");
    saved.initialize<mm_float2>(context, count, "fftSaved");
    vector<mm_float2> values(count);
    for (int i = 0; i < count; i++)
        values[i] = mm_float2(float(i%11)/7, float(i%5)/9);
    original.upload(values);
    for (int i = 0; i < 130; i++) {
        original.copyTo(input);
        fft.execFFT(input, output);
        fft.execFFT(output, input, false);
        if (i == 3)
            input.copyTo(saved);
    }
    vector<mm_float2> actual;
    input.download(actual);
    for (int i = 0; i < count; i++) {
        checkValue(values[i].x, actual[i].x/count);
        checkValue(values[i].y, actual[i].y/count);
    }
    saved.download(actual);
    for (int i = 0; i < count; i++) {
        checkValue(values[i].x, actual[i].x/count);
        checkValue(values[i].y, actual[i].y/count);
    }
}

void testLargeRoundTrip(MetalContext& context, int nx, int ny, int nz, bool realToComplex) {
    const int count = nx*ny*nz;
    ComputeArray input, output;
    input.initialize<float>(context, 2*count, "largeFFTInput");
    output.initialize<float>(context, 2*count, "largeFFTOutput");
    vector<float> values(2*count, 0), actual;
    for (int i = 0; i < (realToComplex ? count : 2*count); i++)
        values[i] = float((i*17)%101-50)/37;
    input.upload(values);
    MetalFFT3D fft(context, nx, ny, nz, realToComplex);
    fft.execFFT(input, output);
    fft.execFFT(output, input, false);
    input.download(actual);
    for (int i = 0; i < (realToComplex ? count : 2*count); i++)
        checkValue(values[i], actual[i]/count);
}

void testValidation(MetalContext& context) {
    ASSERT_EQUAL(1, MetalFFT3D::findLegalDimension(0));
    ASSERT_EQUAL(13, MetalFFT3D::findLegalDimension(13));
    ASSERT_EQUAL(18, MetalFFT3D::findLegalDimension(17));
    expectException([&] { MetalFFT3D invalid(context, 0, 2, 3); });
    expectException([&] { MetalFFT3D invalid(context, -1, 2, 3); });
    expectException([&] { MetalFFT3D invalid(context, numeric_limits<int>::max(), 2, 3); });
    expectException([&] { MetalFFT3D::findLegalDimension(numeric_limits<int>::max()); });
    ComputeArray input, output;
    input.initialize<mm_float2>(context, 8, "fftValid");
    output.initialize<float>(context, 8, "fftTooSmall");
    MetalFFT3D fft(context, 2, 2, 2);
    expectException([&] { fft.execFFT(input, input); });
    expectException([&] { fft.execFFT(input, output); });
}

int main() {
    try {
        System system;
        system.addParticle(1);
        unique_ptr<MetalContext> context;
        try {
            context.reset(new MetalContext(system));
        }
        catch (const OpenMMException& error) {
            if (string(error.what()).find("No Metal device") != string::npos) {
                cout << error.what() << endl;
                return 77;
            }
            throw;
        }
        testValidation(*context);
        const int sizes[][3] = {{1, 1, 1}, {2, 3, 4}, {3, 4, 5}, {3, 5, 7},
                {2, 3, 11}, {2, 3, 13}, {1, 1, 17}, {2, 1, 5}, {1, 3, 5}, {2, 3, 1}};
        for (const auto& size : sizes)
            for (bool realToComplex : {false, true}) {
                cout << "FFT " << size[0] << "x" << size[1] << "x" << size[2]
                     << (realToComplex ? " R2C" : " C2C") << endl;
                testTransform(*context, size[0], size[1], size[2], realToComplex);
            }
        testBatchedPlanReuse(*context);
        for (bool realToComplex : {false, true}) {
            testLargeRoundTrip(*context, 32, 48, 50, realToComplex);
            testLargeRoundTrip(*context, 1, 1, 32768, realToComplex);
        }
        cout << "Metal FFT tests passed" << endl;
    }
    catch (const exception& error) {
        cerr << error.what() << endl;
        return 1;
    }
    return 0;
}
