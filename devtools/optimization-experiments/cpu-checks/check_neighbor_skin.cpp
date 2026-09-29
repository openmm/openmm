#include "NeighborSkinPolicy.h"
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

static unsigned long long checks = 0;
static void require(bool condition, const char* message) {
    ++checks;
    if (!condition)
        throw std::runtime_error(message);
}
static bool sameBits(double a, double b) {
    return std::memcmp(&a, &b, sizeof(double)) == 0;
}

int main() {
    try {
        using OpenMM::NeighborSkinPolicy::selectFraction;
        using OpenMM::NeighborSkinPolicy::padCutoff;
        struct Choice {const char* text; double expected;};
        const Choice choices[] = {{nullptr, 0.08}, {"0", 0.08}, {"0.08", 0.08}, {"0.10", 0.10}, {"0.12", 0.12}};
        const char* rejected[] = {"", "1", "0.1", "0.120", "0.080", "0.07", "0.16", "-0.10", "nan", "NaN", "inf", "INF",
            "+inf", "-inf", " 0.10", "0.10 ", "0.10junk", "0,10", "0.10\n", "\t0.10", "+0.10", "1e-1", "00.10", "0x1p-3"};
        const double max = std::numeric_limits<double>::max();
        std::vector<double> cutoffs = {0.0, -0.0, std::numeric_limits<double>::denorm_min(),
            std::numeric_limits<double>::min(), std::nextafter(std::numeric_limits<double>::min(), 0.0),
            1e-200, 0.8, 1.0, 2.4, 1e200, max/2, max/1.12, std::nextafter(max/1.12, 0.0),
            std::nextafter(max/1.12, max), max};
        unsigned acceptedCases = 0, rejectedCases = 0, paddingCases = 0, overflowCases = 0;
        for (const Choice& choice : choices) {
            double fraction = -1;
            require(selectFraction(choice.text, fraction), "valid selection rejected");
            require(sameBits(fraction, choice.expected), "selected fraction differs");
            ++acceptedCases;
            for (double cutoff : cutoffs) {
                require(std::isfinite(cutoff), "test input must be finite");
                for (bool padding : {false, true}) {
                    volatile double referencePadding = padding ? choice.expected*cutoff : 0.0;
                    const double expected = cutoff+referencePadding;
                    const double actual = padCutoff(cutoff, padding, fraction);
                    require(sameBits(actual, expected), "padding expression differs by representation");
                    require(actual >= cutoff, "padding decreased a nonnegative cutoff");
                    if (!padding)
                        require(std::isfinite(actual), "disabled padding overflowed");
                    if (choice.expected == 0.08) {
                        volatile double oldPadding = padding ? 0.08*cutoff : 0.0;
                        require(sameBits(actual, cutoff+oldPadding), "baseline does not match original expression");
                    }
                    if (std::isinf(actual))
                        ++overflowCases;
                    ++paddingCases;
                }
            }
        }
        for (const char* text : rejected) {
            double fraction = 9.0;
            require(!selectFraction(text, fraction), "invalid selection accepted");
            require(sameBits(fraction, 0.08), "invalid selection did not reset baseline");
            ++rejectedCases;
        }
        std::cout << "{\"pass\":true,\"accepted_cases\":" << acceptedCases
                  << ",\"rejected_cases\":" << rejectedCases << ",\"finite_cutoffs\":" << cutoffs.size()
                  << ",\"padding_cases\":" << paddingCases << ",\"expected_overflow_cases\":" << overflowCases
                  << ",\"checks\":" << checks << ",\"GPU_executed\":false,\"timing_benchmark\":false}\n";
        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
