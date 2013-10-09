#ifndef OPENMM_CPU_NEIGHBORLIST_H_
#define OPENMM_CPU_NEIGHBORLIST_H_

#include "windowsExportCpu.h"
#include <set>
#include <utility>
#include <vector>

namespace OpenMM {
    
class OPENMM_EXPORT_CPU CpuNeighborList {
public:
    void computeNeighborList(int nAtoms,
                             const std::vector<float>& atomLocations, 
                             const std::vector<std::set<int> >& exclusions,
                             const float* periodicBoxSize,
                             bool usePeriodic,
                             float maxDistance,
                             float minDistance = 0.0f,
                             bool reportSymmetricPairs = false);
    const std::vector<std::pair<int, int> >& getNeighbors();
private:
    std::vector<std::pair<int, int> > neighbors;
};

} // namespace OpenMM

#endif // OPENMM_REFERENCE_NEIGHBORLIST_H_
