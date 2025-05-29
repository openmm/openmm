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
#include <string>
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
    XDRFILE_RAII(std::string filename, std::string mode) : xd(xdrfile_open(filename.c_str(), mode.c_str())) {
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

int xtc_natoms(std::string filename) {
    int natoms = 0;
    if (exdrOK != read_xtc_natoms(const_cast<char*>(filename.c_str()), &natoms)) {
        throw std::runtime_error("xtc_read(): could not get natoms\n");
    }
    return natoms;
}

struct XTCFrame {
    int step;
    float time;
    matrix box;
    std::vector<float> positions;
    int natoms;
    const float prec = 1000.0;
    XTCFrame(int natoms) : positions(3 * natoms), natoms(natoms) {
    }

    // Read the next frame from the XTC file and store it in this object
    int readNextFrame(XDRFILE* xd) {
        // Preinitialize in_prec for the precision check below since it may not
        // be modified by read_xtc if the coordinates to read are not compressed
        float in_prec = prec;
        auto* p_ptr = reinterpret_cast<rvec*>(positions.data());
        int status = read_xtc(xd, natoms, &step, &time, box, p_ptr, &in_prec);
        if (status == exdrOK && prec != in_prec) {
            throw std::runtime_error("xtc_read(): precision mismatch\n");
        }
        if (status == exdr3DX) {
            throw std::runtime_error("xtc_read(): XTC file is corrupt\n");
        }
        return status;
    }

    // Write the current frame to the XTC file
    int appendFrameToFile(XDRFILE* xd) {
        auto* p_ptr = reinterpret_cast<rvec*>(positions.data());
        int err = write_xtc(xd, natoms, step, time, box, p_ptr, prec);
        if (err != exdrOK) {
            throw std::runtime_error("xtc_write(): could not write frame\n");
        }
        return err;
    }
};

int xtc_nframes(std::string filename) {
    int nframes = 0;
    int natoms = xtc_natoms(filename);
    if (!natoms) {
        throw std::runtime_error("xtc_read(): natoms is 0\n");
    }
    XDRFILE_RAII xd(filename, "r");
    XTCFrame frame(natoms);
    while (exdrOK == frame.readNextFrame(xd)) {
        nframes++;
    }
    return nframes;
}

void xtc_read(std::string filename, float* coords_arr, float* box_arr, float* time_arr, int* step_arr, int natoms, int nframes) {
    if (natoms == 0) {
        throw std::runtime_error("xtc_read(): natoms is 0\n");
    }
    XDRFILE_RAII xd(filename, "r");
    size_t fidx = 0;
    size_t nframes_long = nframes;
    XTCFrame frame(natoms);
    while (exdrOK == frame.readNextFrame(xd)) {
        time_arr[fidx] = frame.time;
        step_arr[fidx] = frame.step;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                box_arr[fidx + (3 * i + j) * nframes_long] = frame.box[i][j];
            }
        }
        for (size_t aidx = 0; aidx < natoms; aidx++) {
	    size_t xidx = Xf(aidx, fidx, nframes_long);
            size_t yidx = Yf(xidx, nframes_long);
            size_t zidx = Zf(yidx, nframes_long);
            coords_arr[xidx] = frame.positions[3 * aidx + 0];
            coords_arr[yidx] = frame.positions[3 * aidx + 1];
            coords_arr[zidx] = frame.positions[3 * aidx + 2];
        }
        fidx++;
    }
}

static void box_from_array(matrix& matrix_box, float* box, size_t frame, size_t nframes) {
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            matrix_box[i][j] = box[(3 * i + j) * nframes + frame];
        }
    }
}

void xtc_write(std::string filename, int natoms, int nframes, int* step, float* timex, float* pos, float* box) {
    XDRFILE_RAII xd(filename, "a");
    XTCFrame frame(natoms);
    size_t nframes_long = nframes;
    for (size_t f = 0; f < nframes; f++) {
        box_from_array(frame.box, box, f, nframes_long);
        for (size_t i = 0; i < natoms; i++) {
            size_t xidx = Xf(i, f, nframes_long);
            size_t yidx = Yf(xidx, nframes_long);
            size_t zidx = Zf(yidx, nframes_long);
            frame.positions[3 * i + 0] = pos[xidx];
            frame.positions[3 * i + 1] = pos[yidx];
            frame.positions[3 * i + 2] = pos[zidx];
        }
        frame.step = step[f];
        frame.time = timex[f];
        frame.appendFrameToFile(xd);
    }
}

void xtc_rewrite_with_new_timestep(std::string filename_in, std::string filename_out, int first_step, int interval, float dt) {
    int natoms = xtc_natoms(filename_in);
    if (natoms == 0) {
        throw std::runtime_error("xtc_read(): natoms is 0\n");
    }
    XDRFILE_RAII xd_in(filename_in, "r");
    XDRFILE_RAII xd_out(filename_out, "a");
    XTCFrame frame(natoms);
    int i = 0;
    while (exdrOK == frame.readNextFrame(xd_in)) {
        frame.step = first_step + i * interval;
        frame.time = frame.step * dt;
        frame.appendFrameToFile(xd_out);
        i++;
    }
}
