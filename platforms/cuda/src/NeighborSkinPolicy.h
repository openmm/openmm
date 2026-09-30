#ifndef OPENMM_CUDA_NEIGHBOR_SKIN_POLICY_H_
#define OPENMM_CUDA_NEIGHBOR_SKIN_POLICY_H_

#include <cstring>

namespace OpenMM {
namespace NeighborSkinPolicy {

// Exact opt-in values avoid locale-dependent parsing, NaN and unbounded skins.
// False indicates an invalid request; the caller supplies its normal exception.
inline bool selectFraction(const char* value, double& fraction) {
    fraction = 0.08;
    if (value == nullptr || std::strcmp(value, "0") == 0 || std::strcmp(value, "0.08") == 0)
        return true;
    if (std::strcmp(value, "0.10") == 0) {
        fraction = 0.10;
        return true;
    }
    if (std::strcmp(value, "0.12") == 0) {
        fraction = 0.12;
        return true;
    }
    return false;
}

inline double padCutoff(double cutoff, bool usePadding, double fraction) {
    double padding = (usePadding ? fraction*cutoff : 0.0);
    return cutoff+padding;
}

} // namespace NeighborSkinPolicy
} // namespace OpenMM

#endif
