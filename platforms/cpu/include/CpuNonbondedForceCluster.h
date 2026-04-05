#ifndef OPENMM_CPU_NONBONDED_FORCE_CLUSTER_H_
#define OPENMM_CPU_NONBONDED_FORCE_CLUSTER_H_

#include "CpuClusterPairList.h"
#include "openmm/internal/ThreadPool.h"
#include "AlignedArray.h"
#include <vector>

namespace OpenMM {

class CpuNonbondedForceCluster {
public:
    /**
     * Check if the current CPU supports AVX2 (required for the cluster kernel).
     * This is a runtime check via CPUID, compiled in the AVX2 translation unit.
     */
    static bool isSupported();

    static void computeClusterPairIxn(
        const AtomCluster& ci, const AtomCluster& cj,
        uint64_t exclMask, float* forces, double& energy,
        float cutoff2, const fvec4& boxSize, const fvec4& invBoxSize,
        bool usePeriodic, bool useCutoff, float krf, float crf,
        bool useEwald, float alphaEwald,
        const float* erfcTable, const float* ewaldScaleTable,
        float ewaldDXInv, float erfcDXInv, int numTablePoints,
        uint8_t cp_jAtomMask);

    static void calculateDirectIxn(
        const CpuClusterPairList& clusterPairList,
        const float* posq,
        std::vector<AlignedArray<float>>& threadForce,
        double* totalEnergy,
        ThreadPool& threads,
        float cutoffDistance,
        bool usePeriodic,
        const Vec3* periodicBoxVectors,
        bool useCutoff,
        float krf,
        float crf,
        bool useEwald,
        float alphaEwald,
        const float* erfcTable,
        const float* ewaldScaleTable,
        float ewaldDXInv,
        float erfcDXInv,
        int numTablePoints,
        bool useSwitch,
        float switchingDistance);
};

} // namespace OpenMM

#endif // OPENMM_CPU_NONBONDED_FORCE_CLUSTER_H_
