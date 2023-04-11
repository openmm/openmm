/*
MIT License

Copyright (c) 2023 Accellera

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Contributors: Stefan Doerr, Raul P. Pelaez
*/
#include "xtc.h"
#include "xdrfile_xtc.h"
#include <vector>
#include <stdexcept>

// Helper functions to convert between atom/frame and x/y/z indices
static size_t Xf(size_t atom, size_t frame, size_t nframes) {
    return atom * 3 * nframes + frame;
}

static size_t Yf(size_t xidx, size_t nframes) {
    return xidx + nframes;
}

static size_t Zf(size_t yidx, size_t nframes) {
    return yidx + nframes;
}

// Helper struct to manage the XDRFILE pointer
struct XDRFILE_RAII {
    XDRFILE* xd;
    XDRFILE_RAII(const char* filename, const char* mode) : xd(xdrfile_open(filename, mode)) {
        if (!xd) {
            throw std::runtime_error("xtc file: Could not open file");
        }
    }

    ~XDRFILE_RAII() {
        if (xd)
            xdrfile_close(xd);
    }

    operator XDRFILE*() const {
        return xd;
    }
};

int xtc_natoms(const char* filename) {
    int natoms = 0;
    char* filename_c = const_cast<char*>(filename);
    if (exdrOK != read_xtc_natoms(filename_c, &natoms)) {
        throw std::runtime_error("xtc_read(): could not get natoms\n");
    }
    return natoms;
}

int xtc_nframes(const char* filename) {
    int nframes = 0;
    float time;
    int step;
    float prec;
    matrix box;
    int natoms = xtc_natoms(filename);
    if (!natoms) {
        throw std::runtime_error("xtc_read(): natoms is 0\n");
    }
    XDRFILE_RAII xd(filename, "r");
    std::vector<float> p(3 * natoms);
    auto* p_ptr = reinterpret_cast<rvec*>(p.data());
    int retval = 0;
    while (exdrOK == (retval = read_xtc(xd, natoms, &step, &time, box, p_ptr, &prec))) {
        nframes++;
    }
    if (retval == exdr3DX) {
        throw std::runtime_error("xtc_read(): XTC file is corrupt\n");
    }
    return nframes;
}

void xtc_read(const char* filename, float* coords_arr, float* box_arr, float* time_arr, int* step_arr, int natoms, int nframes) {
    if (natoms == 0) {
        throw std::runtime_error("xtc_read(): natoms is 0\n");
    }
    XDRFILE_RAII xd(filename, "r");
    std::vector<float> p(3 * natoms);
    float time;
    int step;
    float prec;
    matrix box;
    int retval = 0;
    int fidx = 0;
    int natoms_garbage = 0;
    auto* p_ptr = reinterpret_cast<rvec*>(p.data());
    while (exdrOK == (retval = read_xtc(xd, natoms_garbage, &step, &time, box, p_ptr, &prec))) {
        time_arr[fidx] = time;
        step_arr[fidx] = step;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                box_arr[fidx + (3 * i + j) * nframes] = box[i][j];
            }
        }
        for (int aidx = 0; aidx < natoms; aidx++) {
            int xidx = Xf(aidx, fidx, nframes);
            int yidx = Yf(xidx, nframes);
            int zidx = Zf(yidx, nframes);
            coords_arr[xidx] = p[3 * aidx + 0];
            coords_arr[yidx] = p[3 * aidx + 1];
            coords_arr[zidx] = p[3 * aidx + 2];
        }
        fidx++;
    }
    if (retval == exdr3DX) {
        throw std::runtime_error("xtc_read(): XTC file is corrupt\n");
    }
}

int xtc_write(const char* filename, int natoms, int nframes, int* step, float* timex, float* pos, float* box) {
    XDRFILE_RAII xd(filename, "a");
    std::vector<float> p(natoms * 3);
    matrix b;
    int err = 0;
    constexpr float prec = 1000;
    for (int f = 0; f < nframes; f++) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                b[i][j] = box[(3 * i + j) * nframes + f];
            }
        }
        for (int i = 0; i < natoms; i++) {
            int xidx = Xf(i, f, nframes);
            int yidx = Yf(xidx, nframes);
            int zidx = Zf(yidx, nframes);
            p[3 * i + 0] = pos[xidx];
            p[3 * i + 1] = pos[yidx];
            p[3 * i + 2] = pos[zidx];
        }
        auto* p_ptr = reinterpret_cast<rvec*>(p.data());
        int err = write_xtc(xd, natoms, (unsigned int)step[f], (float)timex[f], b, p_ptr, prec);
        if (err != 0) {
            throw std::runtime_error("Error writing frame to xtc file\n");
        }
    }
    return err;
}
