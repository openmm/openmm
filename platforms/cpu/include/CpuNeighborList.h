#ifndef OPENMM_CPU_NEIGHBORLIST_H_
#define OPENMM_CPU_NEIGHBORLIST_H_

#include "windowsExportCpu.h"
#include <pthread.h>
#include <set>
#include <utility>
#include <vector>

namespace OpenMM {
    
class OPENMM_EXPORT_CPU CpuNeighborList {
public:
    class ThreadData;
    class VoxelHash;
    CpuNeighborList();
    ~CpuNeighborList();
    void computeNeighborList(int numAtoms, const std::vector<float>& atomLocations, const std::vector<std::set<int> >& exclusions,
            const float* periodicBoxSize, bool usePeriodic, float maxDistance);
    const std::vector<std::pair<int, int> >& getNeighbors();
    /**
     * This routine contains the code executed by each thread.
     */
    void runThread(int index, std::vector<std::pair<int, int> >& threadNeighbors);
private:
    bool isDeleted;
    int numThreads, waitCount;
    std::vector<std::pair<int, int> > neighbors;
    std::vector<pthread_t> thread;
    std::vector<ThreadData*> threadData;
    pthread_cond_t startCondition, endCondition;
    pthread_mutex_t lock;
    // The following variables are used to make information accessible to the individual threads.
    VoxelHash* voxelHash;
    const std::vector<std::set<int> >* exclusions;
    const float* atomLocations;
    const float* periodicBoxSize;
    int numAtoms;
    bool usePeriodic;
    float maxDistance;
};

} // namespace OpenMM

#endif // OPENMM_REFERENCE_NEIGHBORLIST_H_
