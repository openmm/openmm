#define NOMINMAX
#include "ReorderHilbert.h"
#include <climits>
#include <cstdint>
#include <iostream>
#include <stdexcept>

extern "C" bitmask_t referenceHilbertC2I(unsigned, unsigned, const bitmask_t[]);
static std::uint64_t fallbackCalls = 0;
extern "C" bitmask_t hilbert_c2i(unsigned dimensions, unsigned bits, const bitmask_t xyz[]) {
    ++fallbackCalls;
    return referenceHilbertC2I(dimensions, bits, xyz);
}
static int reference(int x, int y, int z) {
    const bitmask_t xyz[3] = {static_cast<bitmask_t>(x), static_cast<bitmask_t>(y), static_cast<bitmask_t>(z)};
    return static_cast<int>(referenceHilbertC2I(3, 8, xyz));
}
static void require(bool condition, const char* reason) {
    if (!condition) throw std::runtime_error(reason);
}
int main() {
    try {
        std::uint64_t points = 0, checksum = 0;
        for (int z=0; z<256; ++z)
            for (int y=0; y<256; ++y)
                for (int x=0; x<256; ++x) {
                    const int original = reference(x,y,z);
                    const int actual = OpenMM::computeReorderHilbert3D8(x,y,z,true);
                    if (original != actual) {
                        std::cerr << "mismatch " << x << ' ' << y << ' ' << z << ' ' << original << ' ' << actual << '\n';
                        return 1;
                    }
                    require(actual >= 0 && actual < (1<<24), "key range");
                    checksum += static_cast<unsigned>(actual);
                    ++points;
                }
        require(fallbackCalls == 0, "fast domain called generic path");
        require(checksum == ((std::uint64_t(1)<<24)*((std::uint64_t(1)<<24)-1))/2, "full cube checksum");
        std::uint64_t outside = 0, disabled = 0;
        const int values[] = {INT_MIN, -1000, -1, 0, 1, 127, 255, 256, 257, 1000, INT_MAX};
        for (int x : values) for (int y : values) for (int z : values) {
            const int original = reference(x,y,z);
            const auto before = fallbackCalls;
            require(OpenMM::computeReorderHilbert3D8(x,y,z,true) == original, "outside fallback key mismatch");
            const bool eligible = static_cast<unsigned>(x)<256 && static_cast<unsigned>(y)<256 && static_cast<unsigned>(z)<256;
            require(fallbackCalls == before+(eligible ? 0 : 1), "outside fallback not selected exactly once");
            if (!eligible) ++outside;
            const auto beforeOff = fallbackCalls;
            require(OpenMM::computeReorderHilbert3D8(x,y,z,false) == original, "disabled result mismatch");
            require(fallbackCalls == beforeOff+1, "disabled path skipped generic");
            ++disabled;
        }
        std::cout << "{\"pass\":true,\"exhaustive_points\":" << points
                  << ",\"outside_fallback_cases\":" << outside
                  << ",\"disabled_cases\":" << disabled
                  << ",\"checksum\":" << checksum
                  << ",\"comparison\":\"bundled original hilbert.cpp compiled as separate translation unit\",\"GPU_executed\":false,\"timing_benchmark\":false}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
