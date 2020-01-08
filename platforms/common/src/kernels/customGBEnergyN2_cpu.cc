#ifdef SUPPORTS_64_BIT_ATOMICS
#define STORE_DERIVATIVE_1(INDEX) ATOMIC_ADD(&derivBuffers[offset+(INDEX-1)*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (deriv##INDEX##_1*0x100000000)));
#define STORE_DERIVATIVE_2(INDEX) ATOMIC_ADD(&derivBuffers[offset+(INDEX-1)*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (local_deriv##INDEX[tgx]*0x100000000)));
#else
#define STORE_DERIVATIVE_1(INDEX) derivBuffers##INDEX[offset] += deriv##INDEX##_1;
#define STORE_DERIVATIVE_2(INDEX) derivBuffers##INDEX[offset] += local_deriv##INDEX[tgx];
#endif

/**
 * Compute a force based on pair interactions.
 */
KERNEL void computeN2Energy(
#ifdef SUPPORTS_64_BIT_ATOMICS
        GLOBAL mm_ulong* RESTRICT forceBuffers,
#else
        GLOBAL real4* RESTRICT forceBuffers,
#endif
        GLOBAL mixed* RESTRICT energyBuffer,
        GLOBAL const real4* RESTRICT posq, GLOBAL const unsigned int* RESTRICT exclusions,
        GLOBAL const ushort2* exclusionTiles, int needEnergy,
#ifdef USE_CUTOFF
        GLOBAL const int* RESTRICT tiles, GLOBAL const unsigned int* RESTRICT interactionCount, real4 periodicBoxSize, real4 invPeriodicBoxSize,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, unsigned int maxTiles, GLOBAL const real4* RESTRICT blockCenter,
        GLOBAL const real4* RESTRICT blockSize, GLOBAL const int* RESTRICT interactingAtoms
#else
        unsigned int numTiles
#endif
        PARAMETER_ARGUMENTS) {
    mixed energy = 0;
    INIT_PARAM_DERIVS
    LOCAL real3 local_pos[LOCAL_BUFFER_SIZE];
    LOCAL real3 local_force[LOCAL_BUFFER_SIZE];
    ATOM_PARAMETER_DATA

    // First loop: process tiles that contain exclusions.
    
    const int firstExclusionTile = FIRST_EXCLUSION_TILE+GROUP_ID*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/NUM_GROUPS;
    const int lastExclusionTile = FIRST_EXCLUSION_TILE+(GROUP_ID+1)*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/NUM_GROUPS;
    for (int pos = firstExclusionTile; pos < lastExclusionTile; pos++) {
        const ushort2 tileIndices = exclusionTiles[pos];
        const unsigned int x = tileIndices.x;
        const unsigned int y = tileIndices.y;

        // Load the data for this tile.

        for (int localAtomIndex = 0; localAtomIndex < TILE_SIZE; localAtomIndex++) {
            unsigned int j = y*TILE_SIZE + localAtomIndex;
            local_pos[localAtomIndex] = trimTo3(posq[j]);
            LOAD_LOCAL_PARAMETERS_FROM_GLOBAL
        }
        if (x == y) {
            // This tile is on the diagonal.

            for (unsigned int tgx = 0; tgx < TILE_SIZE; tgx++) {
#ifdef USE_EXCLUSIONS
                unsigned int excl = exclusions[pos*TILE_SIZE+tgx];
#endif
                unsigned int atom1 = x*TILE_SIZE+tgx;
                real4 force = 0;
                DECLARE_ATOM1_DERIVATIVES
                real3 pos1 = trimTo3(posq[atom1]);
                LOAD_ATOM1_PARAMETERS
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    real3 pos2 = local_pos[j];
                    real3 delta = make_real3(pos2.x-pos1.x, pos2.y-pos1.y, pos2.z-pos1.z);
#ifdef USE_PERIODIC
                    APPLY_PERIODIC_TO_DELTA(delta)
#endif
                    real r2 = dot(delta, delta);
#ifdef USE_CUTOFF
                    if (r2 < CUTOFF_SQUARED) {
#endif
                        real invR = RSQRT(r2);
                        real r = r2*invR;
                        unsigned int atom2 = j;
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
                        energy += 0.5f*tempEnergy;
                        delta *= dEdR;
                        force.x -= delta.x;
                        force.y -= delta.y;
                        force.z -= delta.z;
#ifdef USE_CUTOFF
                    }
#endif
#ifdef USE_EXCLUSIONS
                    excl >>= 1;
#endif
                }

                // Write results.

#ifdef SUPPORTS_64_BIT_ATOMICS
                unsigned int offset = atom1;
                ATOMIC_ADD(&forceBuffers[offset], (mm_ulong) ((mm_long) (force.x*0x100000000)));
                ATOMIC_ADD(&forceBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (force.y*0x100000000)));
                ATOMIC_ADD(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (force.z*0x100000000)));
                STORE_DERIVATIVES_1
#else
                unsigned int offset = atom1 + GROUP_ID*PADDED_NUM_ATOMS;
                forceBuffers[offset].xyz += force.xyz;
                STORE_DERIVATIVES_1
#endif
            }
        }
        else {
            // This is an off-diagonal tile.

            for (int localAtomIndex = 0; localAtomIndex < TILE_SIZE; localAtomIndex++) {
                local_force[localAtomIndex] = 0;
                CLEAR_LOCAL_DERIVATIVES
            }
            for (unsigned int tgx = 0; tgx < TILE_SIZE; tgx++) {
#ifdef USE_EXCLUSIONS
                unsigned int excl = exclusions[pos*TILE_SIZE+tgx];
#endif
                unsigned int atom1 = x*TILE_SIZE+tgx;
                real4 force = 0;
                DECLARE_ATOM1_DERIVATIVES
                real3 pos1 = trimTo3(posq[atom1]);
                LOAD_ATOM1_PARAMETERS
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    real3 pos2 = local_pos[j];
                    real3 delta = make_real3(pos2.x-pos1.x, pos2.y-pos1.y, pos2.z-pos1.z);
#ifdef USE_PERIODIC
                    APPLY_PERIODIC_TO_DELTA(delta)
#endif
                    real r2 = dot(delta.xyz, delta.xyz);
#ifdef USE_CUTOFF
                    if (r2 < CUTOFF_SQUARED) {
#endif
                        real invR = RSQRT(r2);
                        real r = r2*invR;
                        unsigned int atom2 = j;
                        LOAD_ATOM2_PARAMETERS
                        atom2 = y*TILE_SIZE+j;
                        real dEdR = 0;
                        real tempEnergy = 0;
                        const real interactionScale = 1.0f;
#ifdef USE_EXCLUSIONS
                        bool isExcluded = (atom1 >= NUM_ATOMS || atom2 >= NUM_ATOMS || !(excl & 0x1));
                        if (!isExcluded) {
#else
                        if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
#endif
                            COMPUTE_INTERACTION
                            dEdR /= -r;
                        }
                        energy += tempEnergy;
                        delta *= dEdR;
                        force.x -= delta.x;
                        force.y -= delta.y;
                        force.z -= delta.z;
                        atom2 = j;
                        local_force[atom2].xyz += delta.xyz;
                        RECORD_DERIVATIVE_2
#ifdef USE_CUTOFF
                    }
#endif
#ifdef USE_EXCLUSIONS
                    excl >>= 1;
#endif
                }

                // Write results for atom1.

#ifdef SUPPORTS_64_BIT_ATOMICS
                unsigned int offset = atom1;
                ATOMIC_ADD(&forceBuffers[offset], (mm_ulong) ((mm_long) (force.x*0x100000000)));
                ATOMIC_ADD(&forceBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (force.y*0x100000000)));
                ATOMIC_ADD(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (force.z*0x100000000)));
                STORE_DERIVATIVES_1
#else
                unsigned int offset = atom1 + GROUP_ID*PADDED_NUM_ATOMS;
                forceBuffers[offset].xyz += force.xyz;
                STORE_DERIVATIVES_1
#endif
            }

            // Write results.

            for (int tgx = 0; tgx < TILE_SIZE; tgx++) {
#ifdef SUPPORTS_64_BIT_ATOMICS
                unsigned int offset = y*TILE_SIZE+tgx;
                ATOMIC_ADD(&forceBuffers[offset], (mm_ulong) ((mm_long) (local_force[tgx].x*0x100000000)));
                ATOMIC_ADD(&forceBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (local_force[tgx].y*0x100000000)));
                ATOMIC_ADD(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (local_force[tgx].z*0x100000000)));
                STORE_DERIVATIVES_2
#else
                unsigned int offset = y*TILE_SIZE+tgx + GROUP_ID*PADDED_NUM_ATOMS;
                forceBuffers[offset].xyz += local_force[tgx].xyz;
                STORE_DERIVATIVES_2
#endif
            }
        }
    }

    // Second loop: tiles without exclusions, either from the neighbor list (with cutoff) or just enumerating all
    // of them (no cutoff).

#ifdef USE_CUTOFF
    const unsigned int numTiles = interactionCount[0];
    if (numTiles > maxTiles)
        return; // There wasn't enough memory for the neighbor list.
    int pos = (int) (GROUP_ID*(numTiles > maxTiles ? NUM_BLOCKS*((mm_long)NUM_BLOCKS+1)/2 : numTiles)/NUM_GROUPS);
    int end = (int) ((GROUP_ID+1)*(numTiles > maxTiles ? NUM_BLOCKS*((mm_long)NUM_BLOCKS+1)/2 : numTiles)/NUM_GROUPS);
#else
    int pos = (int) (GROUP_ID*(mm_long)numTiles/NUM_GROUPS);
    int end = (int) ((GROUP_ID+1)*(mm_long)numTiles/NUM_GROUPS);
#endif
    int nextToSkip = -1;
    int currentSkipIndex = 0;
    LOCAL int atomIndices[TILE_SIZE];

    while (pos < end) {
        const bool isExcluded = false;
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

        while (nextToSkip < pos) {
            if (currentSkipIndex < NUM_TILES_WITH_EXCLUSIONS) {
                ushort2 tile = exclusionTiles[currentSkipIndex++];
                nextToSkip = tile.x + tile.y*NUM_BLOCKS - tile.y*(tile.y+1)/2;
            }
            else
                nextToSkip = end;
        }
        includeTile = (nextToSkip != pos);
#endif
        if (includeTile) {
            // Load the data for this tile.

            for (int localAtomIndex = 0; localAtomIndex < TILE_SIZE; localAtomIndex++) {
#ifdef USE_CUTOFF
                unsigned int j = interactingAtoms[pos*TILE_SIZE+localAtomIndex];
#else
                unsigned int j = y*TILE_SIZE+localAtomIndex;
#endif
                atomIndices[localAtomIndex] = j;
                if (j < PADDED_NUM_ATOMS) {
                    local_pos[localAtomIndex] = trimTo3(posq[j]);
                    LOAD_LOCAL_PARAMETERS_FROM_GLOBAL
                    local_force[localAtomIndex] = 0;
                    CLEAR_LOCAL_DERIVATIVES
                }
            }
#ifdef USE_PERIODIC
            if (singlePeriodicCopy) {
                // The box is small enough that we can just translate all the atoms into a single periodic
                // box, then skip having to apply periodic boundary conditions later.

                real4 blockCenterX = blockCenter[x];
                for (unsigned int tgx = 0; tgx < TILE_SIZE; tgx++)
                    APPLY_PERIODIC_TO_POS_WITH_CENTER(local_pos[tgx], blockCenterX)
                for (unsigned int tgx = 0; tgx < TILE_SIZE; tgx++) {
                    unsigned int atom1 = x*TILE_SIZE+tgx;
                    real4 force = 0;
                    DECLARE_ATOM1_DERIVATIVES
                    real3 pos1 = trimTo3(posq[atom1]);
                    APPLY_PERIODIC_TO_POS_WITH_CENTER(pos1, blockCenterX)
                    LOAD_ATOM1_PARAMETERS
                    for (unsigned int j = 0; j < TILE_SIZE; j++) {
                        real3 pos2 = local_pos[j];
                        real3 delta = make_real3(pos2.x-pos1.x, pos2.y-pos1.y, pos2.z-pos1.z);
                        real r2 = dot(delta.xyz, delta.xyz);
                        if (atom1 < NUM_ATOMS && atomIndices[j] < NUM_ATOMS && r2 < CUTOFF_SQUARED) {
                            real invR = RSQRT(r2);
                            real r = r2*invR;
                            unsigned int atom2 = j;
                            LOAD_ATOM2_PARAMETERS
                            atom2 = atomIndices[j];
                            real dEdR = 0;
                            real tempEnergy = 0;
                            const real interactionScale = 1.0f;
                            COMPUTE_INTERACTION
                            dEdR /= -r;
                            energy += tempEnergy;
                            delta *= dEdR;
                            force.x -= delta.x;
                            force.y -= delta.y;
                            force.z -= delta.z;
                            atom2 = j;
                            local_force[atom2].xyz += delta.xyz;
                            RECORD_DERIVATIVE_2
                        }
                    }

                    // Write results for atom1.

#ifdef SUPPORTS_64_BIT_ATOMICS
                    unsigned int offset = atom1;
                    ATOMIC_ADD(&forceBuffers[offset], (mm_ulong) ((mm_long) (force.x*0x100000000)));
                    ATOMIC_ADD(&forceBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (force.y*0x100000000)));
                    ATOMIC_ADD(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (force.z*0x100000000)));
                    STORE_DERIVATIVES_1
#else
                    unsigned int offset = atom1 + GROUP_ID*PADDED_NUM_ATOMS;
                    forceBuffers[offset].xyz += force.xyz;
                    STORE_DERIVATIVES_1
#endif
                }
            }
            else
#endif
            {
                // We need to apply periodic boundary conditions separately for each interaction.

                for (unsigned int tgx = 0; tgx < TILE_SIZE; tgx++) {
                    unsigned int atom1 = x*TILE_SIZE+tgx;
                    real4 force = 0;
                    DECLARE_ATOM1_DERIVATIVES
                    real3 pos1 = trimTo3(posq[atom1]);
                    LOAD_ATOM1_PARAMETERS
                    for (unsigned int j = 0; j < TILE_SIZE; j++) {
                        real3 pos2 = local_pos[j];
                        real3 delta = make_real3(pos2.x-pos1.x, pos2.y-pos1.y, pos2.z-pos1.z);
#ifdef USE_PERIODIC
                        APPLY_PERIODIC_TO_DELTA(delta)
#endif
                        real r2 = dot(delta.xyz, delta.xyz);
#ifdef USE_CUTOFF
                        if (atom1 < NUM_ATOMS && atomIndices[j] < NUM_ATOMS && r2 < CUTOFF_SQUARED) {
#else
                        if (atom1 < NUM_ATOMS && atomIndices[j] < NUM_ATOMS) {
#endif
                            real invR = RSQRT(r2);
                            real r = r2*invR;
                            unsigned int atom2 = j;
                            LOAD_ATOM2_PARAMETERS
                            atom2 = atomIndices[j];
                            real dEdR = 0;
                            real tempEnergy = 0;
                            const real interactionScale = 1.0f;
                            COMPUTE_INTERACTION
                            dEdR /= -r;
                            energy += tempEnergy;
                            delta *= dEdR;
                            force.x -= delta.x;
                            force.y -= delta.y;
                            force.z -= delta.z;
                            atom2 = j;
                            local_force[atom2] += delta;
                            RECORD_DERIVATIVE_2
                        }
                    }

                    // Write results for atom1.

#ifdef SUPPORTS_64_BIT_ATOMICS
                    unsigned int offset = atom1;
                    ATOMIC_ADD(&forceBuffers[offset], (mm_ulong) ((mm_long) (force.x*0x100000000)));
                    ATOMIC_ADD(&forceBuffers[offset+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (force.y*0x100000000)));
                    ATOMIC_ADD(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (force.z*0x100000000)));
                    STORE_DERIVATIVES_1
#else
                    unsigned int offset = atom1 + GROUP_ID*PADDED_NUM_ATOMS;
                    forceBuffers[offset].xyz += force.xyz;
                    STORE_DERIVATIVES_1
#endif
                }
            }

            // Write results.

            for (int tgx = 0; tgx < TILE_SIZE; tgx++) {
#ifdef USE_CUTOFF
                unsigned int atom2 = atomIndices[tgx];
#else
                unsigned int atom2 = y*TILE_SIZE + tgx;
#endif
                if (atom2 < PADDED_NUM_ATOMS) {
#ifdef SUPPORTS_64_BIT_ATOMICS
                    ATOMIC_ADD(&forceBuffers[atom2], (mm_ulong) ((mm_long) (local_force[tgx].x*0x100000000)));
                    ATOMIC_ADD(&forceBuffers[atom2+PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (local_force[tgx].y*0x100000000)));
                    ATOMIC_ADD(&forceBuffers[atom2+2*PADDED_NUM_ATOMS], (mm_ulong) ((mm_long) (local_force[tgx].z*0x100000000)));
                    unsigned int offset = atom2;
                    STORE_DERIVATIVES_2
#else
                    unsigned int offset = atom2 + GROUP_ID*PADDED_NUM_ATOMS;
                    forceBuffers[offset].xyz += local_force[tgx].xyz;
                    STORE_DERIVATIVES_2
#endif
                }
            }
        }
        pos++;
    }
    energyBuffer[GLOBAL_ID] += energy;
    SAVE_PARAM_DERIVS
}
