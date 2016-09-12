#ifdef SUPPORTS_64_BIT_ATOMICS
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable
#define STORE_DERIVATIVE_1(INDEX) atom_add(&derivBuffers[offset+(INDEX-1)*PADDED_NUM_ATOMS], (long) (deriv##INDEX##_1*0x100000000));
#define STORE_DERIVATIVE_2(INDEX) atom_add(&derivBuffers[offset+(INDEX-1)*PADDED_NUM_ATOMS], (long) (local_deriv##INDEX[get_local_id(0)]*0x100000000));
#else
#define STORE_DERIVATIVE_1(INDEX) derivBuffers##INDEX[offset] += deriv##INDEX##_1;
#define STORE_DERIVATIVE_2(INDEX) derivBuffers##INDEX[offset] += local_deriv##INDEX[get_local_id(0)];
#endif

/**
 * Compute a force based on pair interactions.
 */
__kernel void computeN2Energy(
#ifdef SUPPORTS_64_BIT_ATOMICS
        __global long* restrict forceBuffers,
#else
        __global real4* restrict forceBuffers,
#endif
        __global mixed* restrict energyBuffer, __local real4* restrict local_force,
        __global const real4* restrict posq, __local real4* restrict local_posq, __global const unsigned int* restrict exclusions,
        __global const ushort2* exclusionTiles, int needEnergy,
#ifdef USE_CUTOFF
        __global const int* restrict tiles, __global const unsigned int* restrict interactionCount, real4 periodicBoxSize, real4 invPeriodicBoxSize,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, unsigned int maxTiles, __global const real4* restrict blockCenter,
        __global const real4* restrict blockSize, __global const int* restrict interactingAtoms
#else
        unsigned int numTiles
#endif
        PARAMETER_ARGUMENTS) {
    const unsigned int totalWarps = get_global_size(0)/TILE_SIZE;
    const unsigned int warp = get_global_id(0)/TILE_SIZE;
    const unsigned int tgx = get_local_id(0) & (TILE_SIZE-1);
    const unsigned int tbx = get_local_id(0) - tgx;
    mixed energy = 0;
    INIT_PARAM_DERIVS

    // First loop: process tiles that contain exclusions.
    
    const unsigned int firstExclusionTile = FIRST_EXCLUSION_TILE+warp*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    const unsigned int lastExclusionTile = FIRST_EXCLUSION_TILE+(warp+1)*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    for (int pos = firstExclusionTile; pos < lastExclusionTile; pos++) {
        const ushort2 tileIndices = exclusionTiles[pos];
        const unsigned int x = tileIndices.x;
        const unsigned int y = tileIndices.y;
        real4 force = 0;
        DECLARE_ATOM1_DERIVATIVES
        unsigned int atom1 = x*TILE_SIZE + tgx;
        real4 posq1 = posq[atom1];
        LOAD_ATOM1_PARAMETERS
#ifdef USE_EXCLUSIONS
        unsigned int excl = exclusions[pos*TILE_SIZE+tgx];
#endif
        if (x == y) {
            // This tile is on the diagonal.

            const unsigned int localAtomIndex = get_local_id(0);
            local_posq[localAtomIndex] = posq1;
            LOAD_LOCAL_PARAMETERS_FROM_1
            SYNC_WARPS;
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                int atom2 = tbx+j;
                real4 posq2 = local_posq[atom2];
                real4 delta = (real4) (posq2.xyz - posq1.xyz, 0);
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
                if (r2 < CUTOFF_SQUARED) {
#endif
                    real invR = RSQRT(r2);
                    real r = r2*invR;
                    LOAD_ATOM2_PARAMETERS
                    atom2 = y*TILE_SIZE+j;
                    real dEdR = 0;
                    real tempEnergy = 0;
                    const real interactionScale = 0.5f;
#ifdef USE_EXCLUSIONS
                    bool isExcluded = !(excl & 0x1);
#endif
                    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && atom1 != atom2) {
                        COMPUTE_INTERACTION
                        dEdR /= -r;
                    }
                    if (needEnergy)
                        energy += 0.5f*tempEnergy;
                    delta.xyz *= dEdR;
                    force.xyz -= delta.xyz;
#ifdef USE_CUTOFF
                }
#endif
#ifdef USE_EXCLUSIONS
                excl >>= 1;
#endif
                SYNC_WARPS;
            }
        }
        else {
            // This is an off-diagonal tile.

            const unsigned int localAtomIndex = get_local_id(0);
            unsigned int j = y*TILE_SIZE + tgx;
            local_posq[localAtomIndex] = posq[j];
            LOAD_LOCAL_PARAMETERS_FROM_GLOBAL
            local_force[localAtomIndex] = 0;
            CLEAR_LOCAL_DERIVATIVES
            SYNC_WARPS;
#ifdef USE_EXCLUSIONS
            excl = (excl >> tgx) | (excl << (TILE_SIZE - tgx));
#endif
            unsigned int tj = tgx;
            for (j = 0; j < TILE_SIZE; j++) {
                int atom2 = tbx+tj;
                real4 posq2 = local_posq[atom2];
                real4 delta = (real4) (posq2.xyz - posq1.xyz, 0);
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
                if (r2 < CUTOFF_SQUARED) {
#endif
                    real invR = RSQRT(r2);
                    real r = r2*invR;
                    LOAD_ATOM2_PARAMETERS
                    atom2 = y*TILE_SIZE+tj;
                    real dEdR = 0;
                    real tempEnergy = 0;
                    const real interactionScale = 1.0f;
#ifdef USE_EXCLUSIONS
                    bool isExcluded = !(excl & 0x1);
#endif
                    if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                        COMPUTE_INTERACTION
                        dEdR /= -r;
                    }
                    if (needEnergy)
                        energy += tempEnergy;
                    delta.xyz *= dEdR;
                    force.xyz -= delta.xyz;
                    atom2 = tbx+tj;
                    local_force[atom2].xyz += delta.xyz;
                    RECORD_DERIVATIVE_2
#ifdef USE_CUTOFF
                }
#endif
#ifdef USE_EXCLUSIONS
                excl >>= 1;
#endif
                tj = (tj + 1) & (TILE_SIZE - 1);
                SYNC_WARPS;
            }
        }
        
        // Write results.
        
#ifdef SUPPORTS_64_BIT_ATOMICS
        unsigned int offset = x*TILE_SIZE + tgx;
        atom_add(&forceBuffers[offset], (long) (force.x*0x100000000));
        atom_add(&forceBuffers[offset+PADDED_NUM_ATOMS], (long) (force.y*0x100000000));
        atom_add(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (long) (force.z*0x100000000));
        STORE_DERIVATIVES_1
        if (x != y) {
            offset = y*TILE_SIZE + tgx;
            atom_add(&forceBuffers[offset], (long) (local_force[get_local_id(0)].x*0x100000000));
            atom_add(&forceBuffers[offset+PADDED_NUM_ATOMS], (long) (local_force[get_local_id(0)].y*0x100000000));
            atom_add(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (long) (local_force[get_local_id(0)].z*0x100000000));
            STORE_DERIVATIVES_2
        }
#else
        unsigned int offset1 = x*TILE_SIZE + tgx + warp*PADDED_NUM_ATOMS;
        unsigned int offset2 = y*TILE_SIZE + tgx + warp*PADDED_NUM_ATOMS;
        unsigned int offset = offset1;
        forceBuffers[offset1].xyz += force.xyz;
        STORE_DERIVATIVES_1
        if (x != y) {
            offset = offset2;
            forceBuffers[offset2] += (real4) (local_force[get_local_id(0)].x, local_force[get_local_id(0)].y, local_force[get_local_id(0)].z, 0.0f);
            STORE_DERIVATIVES_2
        }
#endif
    }

    // Second loop: tiles without exclusions, either from the neighbor list (with cutoff) or just enumerating all
    // of them (no cutoff).

#ifdef USE_CUTOFF
    unsigned int numTiles = interactionCount[0];
    if (numTiles > maxTiles)
        return; // There wasn't enough memory for the neighbor list.
    int pos = (int) (warp*(numTiles > maxTiles ? NUM_BLOCKS*((long)NUM_BLOCKS+1)/2 : (long)numTiles)/totalWarps);
    int end = (int) ((warp+1)*(numTiles > maxTiles ? NUM_BLOCKS*((long)NUM_BLOCKS+1)/2 : (long)numTiles)/totalWarps);
#else
    int pos = (int) (warp*(long)numTiles/totalWarps);
    int end = (int) ((warp+1)*(long)numTiles/totalWarps);
#endif
    int skipBase = 0;
    int currentSkipIndex = tbx;
    __local int atomIndices[FORCE_WORK_GROUP_SIZE];
    __local volatile int skipTiles[FORCE_WORK_GROUP_SIZE];
    skipTiles[get_local_id(0)] = -1;

    while (pos < end) {
        const bool isExcluded = false;
        real4 force = 0;
        DECLARE_ATOM1_DERIVATIVES
        bool includeTile = true;

        // Extract the coordinates of this tile.
        
        int x, y;
        bool singlePeriodicCopy = false;
#ifdef USE_CUTOFF
        x = tiles[pos];
        real4 blockSizeX = blockSize[x];
        singlePeriodicCopy = (0.5f*periodicBoxSize.x-blockSizeX.x >= CUTOFF &&
                              0.5f*periodicBoxSize.y-blockSizeX.y >= CUTOFF &&
                              0.5f*periodicBoxSize.z-blockSizeX.z >= CUTOFF);
#else
        y = (int) floor(NUM_BLOCKS+0.5f-SQRT((NUM_BLOCKS+0.5f)*(NUM_BLOCKS+0.5f)-2*pos));
        x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
        if (x < y || x >= NUM_BLOCKS) { // Occasionally happens due to roundoff error.
            y += (x < y ? -1 : 1);
            x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
        }

        // Skip over tiles that have exclusions, since they were already processed.

        SYNC_WARPS;
        while (skipTiles[tbx+TILE_SIZE-1] < pos) {
            SYNC_WARPS;
            if (skipBase+tgx < NUM_TILES_WITH_EXCLUSIONS) {
                ushort2 tile = exclusionTiles[skipBase+tgx];
                skipTiles[get_local_id(0)] = tile.x + tile.y*NUM_BLOCKS - tile.y*(tile.y+1)/2;
            }
            else
                skipTiles[get_local_id(0)] = end;
            skipBase += TILE_SIZE;            
            currentSkipIndex = tbx;
            SYNC_WARPS;
        }
        while (skipTiles[currentSkipIndex] < pos)
            currentSkipIndex++;
        includeTile = (skipTiles[currentSkipIndex] != pos);
#endif
        if (includeTile) {
            unsigned int atom1 = x*TILE_SIZE + tgx;

            // Load atom data for this tile.
            
            real4 posq1 = posq[atom1];
            LOAD_ATOM1_PARAMETERS
            const unsigned int localAtomIndex = get_local_id(0);
#ifdef USE_CUTOFF
            unsigned int j = interactingAtoms[pos*TILE_SIZE+tgx];
#else
            unsigned int j = y*TILE_SIZE + tgx;
#endif
            atomIndices[get_local_id(0)] = j;
            if (j < PADDED_NUM_ATOMS) {
                local_posq[localAtomIndex] = posq[j];
                LOAD_LOCAL_PARAMETERS_FROM_GLOBAL
                local_force[localAtomIndex] = 0;
                CLEAR_LOCAL_DERIVATIVES
            }
            SYNC_WARPS;
#ifdef USE_PERIODIC
            if (singlePeriodicCopy) {
                // The box is small enough that we can just translate all the atoms into a single periodic
                // box, then skip having to apply periodic boundary conditions later.

                real4 blockCenterX = blockCenter[x];
                APPLY_PERIODIC_TO_POS_WITH_CENTER(posq1, blockCenterX)
                APPLY_PERIODIC_TO_POS_WITH_CENTER(local_posq[get_local_id(0)], blockCenterX)
                SYNC_WARPS;
                unsigned int tj = tgx;
                for (j = 0; j < TILE_SIZE; j++) {
                    int atom2 = tbx+tj;
                    real4 posq2 = local_posq[atom2];
                    real4 delta = (real4) (posq2.xyz - posq1.xyz, 0);
                    real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                    if (r2 < CUTOFF_SQUARED) {
                        real invR = RSQRT(r2);
                        real r = r2*invR;
                        LOAD_ATOM2_PARAMETERS
                        atom2 = atomIndices[tbx+tj];
                        real dEdR = 0;
                        real tempEnergy = 0;
                        const real interactionScale = 1.0f;
                        if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                            COMPUTE_INTERACTION
                            dEdR /= -r;
                        }
                        if (needEnergy)
                            energy += tempEnergy;
                        delta.xyz *= dEdR;
                        force.xyz -= delta.xyz;
                        atom2 = tbx+tj;
                        local_force[atom2].xyz += delta.xyz;
                        RECORD_DERIVATIVE_2
                    }
                    tj = (tj + 1) & (TILE_SIZE - 1);
                    SYNC_WARPS;
                }
            }
            else
#endif
            {
                // We need to apply periodic boundary conditions separately for each interaction.

                unsigned int tj = tgx;
                for (j = 0; j < TILE_SIZE; j++) {
                    int atom2 = tbx+tj;
                    real4 posq2 = local_posq[atom2];
                    real4 delta = (real4) (posq2.xyz - posq1.xyz, 0);
#ifdef USE_PERIODIC
                    APPLY_PERIODIC_TO_DELTA(delta)
#endif
                    real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
                    if (r2 < CUTOFF_SQUARED) {
#endif
                        real invR = RSQRT(r2);
                        real r = r2*invR;
                        LOAD_ATOM2_PARAMETERS
                        atom2 = atomIndices[tbx+tj];
                        real dEdR = 0;
                        real tempEnergy = 0;
                        const real interactionScale = 1.0f;
                        if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                            COMPUTE_INTERACTION
                            dEdR /= -r;
                        }
                        if (needEnergy)
                            energy += tempEnergy;
                        delta.xyz *= dEdR;
                        force.xyz -= delta.xyz;
                        atom2 = tbx+tj;
                        local_force[atom2].xyz += delta.xyz;
                        RECORD_DERIVATIVE_2
#ifdef USE_CUTOFF
                    }
#endif
                    tj = (tj + 1) & (TILE_SIZE - 1);
                    SYNC_WARPS;
                }
            }
        
            // Write results.
        
#ifdef USE_CUTOFF
            unsigned int atom2 = atomIndices[get_local_id(0)];
#else
            unsigned int atom2 = y*TILE_SIZE + tgx;
#endif
#ifdef SUPPORTS_64_BIT_ATOMICS
            atom_add(&forceBuffers[atom1], (long) (force.x*0x100000000));
            atom_add(&forceBuffers[atom1+PADDED_NUM_ATOMS], (long) (force.y*0x100000000));
            atom_add(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], (long) (force.z*0x100000000));
            unsigned int offset = atom1;
            STORE_DERIVATIVES_1
            if (atom2 < PADDED_NUM_ATOMS) {
                atom_add(&forceBuffers[atom2], (long) (local_force[get_local_id(0)].x*0x100000000));
                atom_add(&forceBuffers[atom2+PADDED_NUM_ATOMS], (long) (local_force[get_local_id(0)].y*0x100000000));
                atom_add(&forceBuffers[atom2+2*PADDED_NUM_ATOMS], (long) (local_force[get_local_id(0)].z*0x100000000));
                offset = atom2;
                STORE_DERIVATIVES_2
            }
#else
            unsigned int offset1 = atom1 + warp*PADDED_NUM_ATOMS;
            unsigned int offset2 = atom2 + warp*PADDED_NUM_ATOMS;
            forceBuffers[offset1].xyz += force.xyz;
            unsigned int offset = offset1;
            STORE_DERIVATIVES_1
            if (atom2 < PADDED_NUM_ATOMS) {
                forceBuffers[offset2] += (real4) (local_force[get_local_id(0)].x, local_force[get_local_id(0)].y, local_force[get_local_id(0)].z, 0.0f);
                offset = offset2;
                STORE_DERIVATIVES_2
            }
#endif
        }
        pos++;
    }
    energyBuffer[get_global_id(0)] += energy;
    SAVE_PARAM_DERIVS
}
