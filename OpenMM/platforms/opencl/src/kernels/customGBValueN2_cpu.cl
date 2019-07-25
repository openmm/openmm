#ifdef SUPPORTS_64_BIT_ATOMICS
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable
#endif

/**
 * Compute a value based on pair interactions.
 */
__kernel void computeN2Value(__global const real4* restrict posq, __local real4* restrict local_posq, __global const unsigned int* restrict exclusions,
        __global const ushort2* exclusionTiles,
#ifdef SUPPORTS_64_BIT_ATOMICS
        __global long* restrict global_value,
#else
        __global real* restrict global_value,
#endif
        __local real* restrict local_value,
#ifdef USE_CUTOFF
        __global const int* restrict tiles, __global const unsigned int* restrict interactionCount, real4 periodicBoxSize, real4 invPeriodicBoxSize,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, unsigned int maxTiles, __global const real4* restrict blockCenter,
        __global const real4* restrict blockSize, __global const int* restrict interactingAtoms
#else
        unsigned int numTiles
#endif
        PARAMETER_ARGUMENTS) {

    // First loop: process tiles that contain exclusions.
    
    const unsigned int firstExclusionTile = FIRST_EXCLUSION_TILE+get_group_id(0)*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/get_num_groups(0);
    const unsigned int lastExclusionTile = FIRST_EXCLUSION_TILE+(get_group_id(0)+1)*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/get_num_groups(0);
    for (int pos = firstExclusionTile; pos < lastExclusionTile; pos++) {
        const ushort2 tileIndices = exclusionTiles[pos];
        const unsigned int x = tileIndices.x;
        const unsigned int y = tileIndices.y;

        // Load the data for this tile.

        for (int localAtomIndex = 0; localAtomIndex < TILE_SIZE; localAtomIndex++) {
            unsigned int j = y*TILE_SIZE + localAtomIndex;
            local_posq[localAtomIndex] = posq[j];
            LOAD_LOCAL_PARAMETERS_FROM_GLOBAL
        }
        if (x == y) {
            // This tile is on the diagonal.

            for (unsigned int tgx = 0; tgx < TILE_SIZE; tgx++) {
#ifdef USE_EXCLUSIONS
                unsigned int excl = exclusions[pos*TILE_SIZE+tgx];
#endif
                unsigned int atom1 = x*TILE_SIZE+tgx;
                real value = 0;
                real4 posq1 = posq[atom1];
                LOAD_ATOM1_PARAMETERS
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    real4 posq2 = local_posq[j];
                    real4 delta = (real4) (posq2.xyz - posq1.xyz, 0);
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
                        real tempValue1 = 0;
                        real tempValue2 = 0;
#ifdef USE_EXCLUSIONS
                        bool isExcluded = (atom1 >= NUM_ATOMS || atom2 >= NUM_ATOMS || !(excl & 0x1));
                        if (!isExcluded && atom1 != atom2) {
#else
                        if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && atom1 != atom2) {
#endif
                            COMPUTE_VALUE
                        }
                        value += tempValue1;
                        ADD_TEMP_DERIVS1
#ifdef USE_CUTOFF
                    }
#endif
#ifdef USE_EXCLUSIONS
                    excl >>= 1;
#endif
                }

                // Write results.

#ifdef SUPPORTS_64_BIT_ATOMICS
                unsigned int offset1 = atom1;
                atom_add(&global_value[offset1], (long) (value*0x100000000));
#else
                unsigned int offset1 = atom1 + get_group_id(0)*PADDED_NUM_ATOMS;
                global_value[offset1] += value;
#endif
                STORE_PARAM_DERIVS1
            }
        }
        else {
            // This is an off-diagonal tile.

            for (int tgx = 0; tgx < TILE_SIZE; tgx++)
                local_value[tgx] = 0;
            for (unsigned int tgx = 0; tgx < TILE_SIZE; tgx++) {
#ifdef USE_EXCLUSIONS
                unsigned int excl = exclusions[pos*TILE_SIZE+tgx];
#endif
                unsigned int atom1 = x*TILE_SIZE+tgx;
                real value = 0;
                real4 posq1 = posq[atom1];
                LOAD_ATOM1_PARAMETERS
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    real4 posq2 = local_posq[j];
                    real4 delta = (real4) (posq2.xyz - posq1.xyz, 0);
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
                        real tempValue1 = 0;
                        real tempValue2 = 0;
#ifdef USE_EXCLUSIONS
                        bool isExcluded = (atom1 >= NUM_ATOMS || atom2 >= NUM_ATOMS || !(excl & 0x1));
                        if (!isExcluded) {
#else
                        if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS && atom1 != atom2) {
#endif
                            COMPUTE_VALUE
                        }
                        value += tempValue1;
                        local_value[j] += tempValue2;
                        ADD_TEMP_DERIVS1
                        ADD_TEMP_DERIVS2
#ifdef USE_CUTOFF
                    }
#endif
#ifdef USE_EXCLUSIONS
                    excl >>= 1;
#endif
                }

                // Write results for atom1.

#ifdef SUPPORTS_64_BIT_ATOMICS
                unsigned int offset1 = atom1;
                atom_add(&global_value[offset1], (long) (value*0x100000000));
#else
                unsigned int offset1 = atom1 + get_group_id(0)*PADDED_NUM_ATOMS;
                global_value[offset1] += value;
#endif
                STORE_PARAM_DERIVS1
            }

            // Write results.

            for (int tgx = 0; tgx < TILE_SIZE; tgx++) {
#ifdef SUPPORTS_64_BIT_ATOMICS
                unsigned int offset2 = y*TILE_SIZE+tgx;
                atom_add(&global_value[offset2], (long) (local_value[tgx]*0x100000000));
#else
                unsigned int offset2 = y*TILE_SIZE+tgx + get_group_id(0)*PADDED_NUM_ATOMS;
                global_value[offset2] += local_value[tgx];
#endif
                STORE_PARAM_DERIVS2
            }
        }
    }

    // Second loop: tiles without exclusions, either from the neighbor list (with cutoff) or just enumerating all
    // of them (no cutoff).

#ifdef USE_CUTOFF
    const unsigned int numTiles = interactionCount[0];
    if (numTiles > maxTiles)
        return; // There wasn't enough memory for the neighbor list.
    int pos = (int) (get_group_id(0)*(numTiles > maxTiles ? NUM_BLOCKS*((long)NUM_BLOCKS+1)/2 : numTiles)/get_num_groups(0));
    int end = (int) ((get_group_id(0)+1)*(numTiles > maxTiles ? NUM_BLOCKS*((long)NUM_BLOCKS+1)/2 : numTiles)/get_num_groups(0));
#else
    int pos = (int) (get_group_id(0)*(long)numTiles/get_num_groups(0));
    int end = (int) ((get_group_id(0)+1)*(long)numTiles/get_num_groups(0));
#endif
    int nextToSkip = -1;
    int currentSkipIndex = 0;
    __local int atomIndices[TILE_SIZE];

    while (pos < end) {
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
                    local_posq[localAtomIndex] = posq[j];
                    LOAD_LOCAL_PARAMETERS_FROM_GLOBAL
                    local_value[localAtomIndex] = 0;
                }
            }
#ifdef USE_PERIODIC
            if (singlePeriodicCopy) {
                // The box is small enough that we can just translate all the atoms into a single periodic
                // box, then skip having to apply periodic boundary conditions later.

                real4 blockCenterX = blockCenter[x];
                for (unsigned int tgx = 0; tgx < TILE_SIZE; tgx++)
                    APPLY_PERIODIC_TO_POS_WITH_CENTER(local_posq[tgx], blockCenterX)
                for (unsigned int tgx = 0; tgx < TILE_SIZE; tgx++) {
                    unsigned int atom1 = x*TILE_SIZE+tgx;
                    real value = 0;
                    real4 posq1 = posq[atom1];
                    APPLY_PERIODIC_TO_POS_WITH_CENTER(posq1, blockCenterX)
                    LOAD_ATOM1_PARAMETERS
                    for (unsigned int j = 0; j < TILE_SIZE; j++) {
                        real4 posq2 = local_posq[j];
                        real4 delta = (real4) (posq2.xyz - posq1.xyz, 0);
                        real r2 = dot(delta.xyz, delta.xyz);
                        if (atom1 < NUM_ATOMS && atomIndices[j] < NUM_ATOMS && r2 < CUTOFF_SQUARED) {
                            real invR = RSQRT(r2);
                            real r = r2*invR;
                            unsigned int atom2 = j;
                            LOAD_ATOM2_PARAMETERS
                            atom2 = atomIndices[j];
                            real tempValue1 = 0;
                            real tempValue2 = 0;
                            COMPUTE_VALUE
                            value += tempValue1;
                            local_value[j] += tempValue2;
                            ADD_TEMP_DERIVS1
                            ADD_TEMP_DERIVS2
                        }
                    }

                    // Write results for atom1.

#ifdef SUPPORTS_64_BIT_ATOMICS
                    unsigned int offset1 = atom1;
                    atom_add(&global_value[offset1], (long) (value*0x100000000));
#else
                    unsigned int offset1 = atom1 + get_group_id(0)*PADDED_NUM_ATOMS;
                    global_value[offset1] += value;
#endif
                    STORE_PARAM_DERIVS1
                }
            }
            else
#endif
            {
                // We need to apply periodic boundary conditions separately for each interaction.

                for (unsigned int tgx = 0; tgx < TILE_SIZE; tgx++) {
                    unsigned int atom1 = x*TILE_SIZE+tgx;
                    real value = 0;
                    real4 posq1 = posq[atom1];
                    LOAD_ATOM1_PARAMETERS
                    for (unsigned int j = 0; j < TILE_SIZE; j++) {
                        real4 posq2 = local_posq[j];
                        real4 delta = (real4) (posq2.xyz - posq1.xyz, 0);
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
                            real tempValue1 = 0;
                            real tempValue2 = 0;
                            COMPUTE_VALUE
                            value += tempValue1;
                            local_value[j] += tempValue2;
                            ADD_TEMP_DERIVS1
                            ADD_TEMP_DERIVS2
                        }
                    }

                    // Write results for atom1.

#ifdef SUPPORTS_64_BIT_ATOMICS
                    unsigned int offset1 = atom1;
                    atom_add(&global_value[offset1], (long) (value*0x100000000));
#else
                    unsigned int offset1 = atom1 + get_group_id(0)*PADDED_NUM_ATOMS;
                    global_value[offset1] += value;
#endif
                    STORE_PARAM_DERIVS1
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
                    unsigned int offset2 = atom2;
                    atom_add(&global_value[offset2], (long) (local_value[tgx]*0x100000000));
#else
                    unsigned int offset2 = atom2 + get_group_id(0)*PADDED_NUM_ATOMS;
                    global_value[offset2] += local_value[tgx];
#endif
                    STORE_PARAM_DERIVS2
                }
            }
        }
        pos++;
    }
}
