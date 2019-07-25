#ifdef SUPPORTS_64_BIT_ATOMICS
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable
#endif

typedef struct {
    real x, y, z;
    real q;
    real fx, fy, fz;
    ATOM_PARAMETER_DATA
} AtomData;

/**
 * Compute nonbonded interactions.
 */

__kernel void computeNonbonded(
#ifdef SUPPORTS_64_BIT_ATOMICS
        __global long* restrict forceBuffers,
#else
        __global real4* restrict forceBuffers,
#endif
        __global mixed* restrict energyBuffer, __global const real4* restrict posq, __global const unsigned int* restrict exclusions,
        __global const ushort2* restrict exclusionTiles, unsigned int startTileIndex, unsigned int numTileIndices
#ifdef USE_CUTOFF
        , __global const int* restrict tiles, __global const unsigned int* restrict interactionCount, real4 periodicBoxSize, real4 invPeriodicBoxSize,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, unsigned int maxTiles, __global const real4* restrict blockCenter,
        __global const real4* restrict blockSize, __global const int* restrict interactingAtoms
#endif
        PARAMETER_ARGUMENTS) {
    mixed energy = 0;
    INIT_DERIVATIVES
    __local AtomData localData[TILE_SIZE];

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
            real4 tempPosq = posq[j];
            localData[localAtomIndex].x = tempPosq.x;
            localData[localAtomIndex].y = tempPosq.y;
            localData[localAtomIndex].z = tempPosq.z;
            localData[localAtomIndex].q = tempPosq.w;
            LOAD_LOCAL_PARAMETERS_FROM_GLOBAL
        }
        const bool hasExclusions = true;
        if (x == y) {
            // This tile is on the diagonal.

            for (unsigned int tgx = 0; tgx < TILE_SIZE; tgx++) {
#ifdef USE_EXCLUSIONS
                unsigned int excl = exclusions[pos*TILE_SIZE+tgx];
#endif
                unsigned int atom1 = x*TILE_SIZE+tgx;
                real4 force = 0;
                real4 posq1 = posq[atom1];
                LOAD_ATOM1_PARAMETERS
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    real4 posq2 = (real4) (localData[j].x, localData[j].y, localData[j].z, localData[j].q);
                    real4 delta = (real4) (posq2.xyz - posq1.xyz, 0);
#ifdef USE_PERIODIC
                    APPLY_PERIODIC_TO_DELTA(delta)
#endif
                    real r2 = dot(delta.xyz, delta.xyz);
#ifdef USE_CUTOFF
                    if (r2 < MAX_CUTOFF*MAX_CUTOFF) {
#endif
                        real invR = RSQRT(r2);
                        real r = r2*invR;
                        unsigned int atom2 = j;
                        LOAD_ATOM2_PARAMETERS
                        atom2 = y*TILE_SIZE+j;
#ifdef USE_SYMMETRIC
                        real dEdR = 0;
#else
                        real4 dEdR1 = (real4) 0;
                        real4 dEdR2 = (real4) 0;
#endif
#ifdef USE_EXCLUSIONS
                        bool isExcluded = (atom1 >= NUM_ATOMS || atom2 >= NUM_ATOMS || !(excl & 0x1));
#endif
                        real tempEnergy = 0;
                        const real interactionScale = 0.5f;
                        COMPUTE_INTERACTION
                        energy += 0.5f*tempEnergy;
#ifdef USE_SYMMETRIC
                        force.xyz -= delta.xyz*dEdR;
#else
                        force.xyz -= dEdR1.xyz;
#endif
#ifdef USE_CUTOFF
                    }
#endif
#ifdef USE_EXCLUSIONS
                    excl >>= 1;
#endif
                }

                // Write results.

#ifdef SUPPORTS_64_BIT_ATOMICS
                atom_add(&forceBuffers[atom1], (long) (force.x*0x100000000));
                atom_add(&forceBuffers[atom1+PADDED_NUM_ATOMS], (long) (force.y*0x100000000));
                atom_add(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], (long) (force.z*0x100000000));
#else
                unsigned int offset = atom1 + get_group_id(0)*PADDED_NUM_ATOMS;
                forceBuffers[offset].xyz = forceBuffers[offset].xyz+force.xyz;
#endif
            }
        }
        else {
            // This is an off-diagonal tile.

            for (int tgx = 0; tgx < TILE_SIZE; tgx++) {
                localData[tgx].fx = 0;
                localData[tgx].fy = 0;
                localData[tgx].fz = 0;
            }
            for (unsigned int tgx = 0; tgx < TILE_SIZE; tgx++) {
                unsigned int atom1 = x*TILE_SIZE+tgx;
                real4 force = 0;
                real4 posq1 = posq[atom1];
                LOAD_ATOM1_PARAMETERS
#ifdef USE_EXCLUSIONS
                unsigned int excl = exclusions[pos*TILE_SIZE+tgx];
#endif
                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    real4 posq2 = (real4) (localData[j].x, localData[j].y, localData[j].z, localData[j].q);
                    real4 delta = (real4) (posq2.xyz - posq1.xyz, 0);
#ifdef USE_PERIODIC
                    APPLY_PERIODIC_TO_DELTA(delta)
#endif
                    real r2 = dot(delta.xyz, delta.xyz);
#ifdef USE_CUTOFF
                    if (r2 < MAX_CUTOFF*MAX_CUTOFF) {
#endif
                        real invR = RSQRT(r2);
                        real r = r2*invR;
                        unsigned int atom2 = j;
                        LOAD_ATOM2_PARAMETERS
                        atom2 = y*TILE_SIZE+j;
#ifdef USE_SYMMETRIC
                        real dEdR = 0;
#else
                        real4 dEdR1 = (real4) 0;
                        real4 dEdR2 = (real4) 0;
#endif
#ifdef USE_EXCLUSIONS
                        bool isExcluded = (atom1 >= NUM_ATOMS || atom2 >= NUM_ATOMS || !(excl & 0x1));
#endif
                        real tempEnergy = 0;
                        const real interactionScale = 1.0f;
                        COMPUTE_INTERACTION
                        energy += tempEnergy;
#ifdef USE_SYMMETRIC
                        delta.xyz *= dEdR;
                        force.xyz -= delta.xyz;
                        localData[j].fx += delta.x;
                        localData[j].fy += delta.y;
                        localData[j].fz += delta.z;
#else
                        force.xyz -= dEdR1.xyz;
                        localData[j].fx += dEdR2.x;
                        localData[j].fy += dEdR2.y;
                        localData[j].fz += dEdR2.z;
#endif
#ifdef USE_CUTOFF
                    }
#endif
#ifdef USE_EXCLUSIONS
                    excl >>= 1;
#endif
                }

               // Write results for atom1.

#ifdef SUPPORTS_64_BIT_ATOMICS
                atom_add(&forceBuffers[atom1], (long) (force.x*0x100000000));
                atom_add(&forceBuffers[atom1+PADDED_NUM_ATOMS], (long) (force.y*0x100000000));
                atom_add(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], (long) (force.z*0x100000000));
#else
                unsigned int offset = atom1 + get_group_id(0)*PADDED_NUM_ATOMS;
                forceBuffers[offset].xyz = forceBuffers[offset].xyz+force.xyz;
#endif
            }

            // Write results.

            for (int tgx = 0; tgx < TILE_SIZE; tgx++) {
#ifdef SUPPORTS_64_BIT_ATOMICS
                unsigned int offset = y*TILE_SIZE + tgx;
                atom_add(&forceBuffers[offset], (long) (localData[tgx].fx*0x100000000));
                atom_add(&forceBuffers[offset+PADDED_NUM_ATOMS], (long) (localData[tgx].fy*0x100000000));
                atom_add(&forceBuffers[offset+2*PADDED_NUM_ATOMS], (long) (localData[tgx].fz*0x100000000));
#else
                unsigned int offset = y*TILE_SIZE+tgx + get_group_id(0)*PADDED_NUM_ATOMS;
                real4 f = forceBuffers[offset];
                f.x += localData[tgx].fx;
                f.y += localData[tgx].fy;
                f.z += localData[tgx].fz;
                forceBuffers[offset] = f;
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
    int pos = (int) (numTiles > maxTiles ? (unsigned int) (startTileIndex+get_group_id(0)*(long)numTileIndices/get_num_groups(0)) : get_group_id(0)*(long)numTiles/get_num_groups(0));
    int end = (int) (numTiles > maxTiles ? (unsigned int) (startTileIndex+(get_group_id(0)+1)*(long)numTileIndices/get_num_groups(0)) : (get_group_id(0)+1)*(long)numTiles/get_num_groups(0));
#else
    const unsigned int numTiles = numTileIndices;
    int pos = (int) (startTileIndex+get_group_id(0)*(long)numTiles/get_num_groups(0));
    int end = (int) (startTileIndex+(get_group_id(0)+1)*(long)numTiles/get_num_groups(0));
#endif
    int nextToSkip = -1;
    int currentSkipIndex = 0;
    __local int atomIndices[TILE_SIZE];

    while (pos < end) {
        const bool hasExclusions = false;
        bool includeTile = true;

        // Extract the coordinates of this tile.

        int x, y;
        bool singlePeriodicCopy = false;
#ifdef USE_CUTOFF
        x = tiles[pos];
        real4 blockSizeX = blockSize[x];
        singlePeriodicCopy = (0.5f*periodicBoxSize.x-blockSizeX.x >= MAX_CUTOFF &&
                              0.5f*periodicBoxSize.y-blockSizeX.y >= MAX_CUTOFF &&
                              0.5f*periodicBoxSize.z-blockSizeX.z >= MAX_CUTOFF);
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
                    real4 tempPosq = posq[j];
                    localData[localAtomIndex].x = tempPosq.x;
                    localData[localAtomIndex].y = tempPosq.y;
                    localData[localAtomIndex].z = tempPosq.z;
                    localData[localAtomIndex].q = tempPosq.w;
                    LOAD_LOCAL_PARAMETERS_FROM_GLOBAL
                    localData[localAtomIndex].fx = 0;
                    localData[localAtomIndex].fy = 0;
                    localData[localAtomIndex].fz = 0;
                }
            }
#ifdef USE_PERIODIC
            if (singlePeriodicCopy) {
                // The box is small enough that we can just translate all the atoms into a single periodic
                // box, then skip having to apply periodic boundary conditions later.

                real4 blockCenterX = blockCenter[x];
                for (unsigned int tgx = 0; tgx < TILE_SIZE; tgx++) {
                    APPLY_PERIODIC_TO_POS_WITH_CENTER(localData[tgx], blockCenterX)
                }
                for (unsigned int tgx = 0; tgx < TILE_SIZE; tgx++) {
                    unsigned int atom1 = x*TILE_SIZE+tgx;
                    real4 force = 0;
                    real4 posq1 = posq[atom1];
                    APPLY_PERIODIC_TO_POS_WITH_CENTER(posq1, blockCenterX)
                    LOAD_ATOM1_PARAMETERS
                    for (unsigned int j = 0; j < TILE_SIZE; j++) {
                        real4 posq2 = (real4) (localData[j].x, localData[j].y, localData[j].z, localData[j].q);
                        real4 delta = (real4) (posq2.xyz - posq1.xyz, 0);
                        real r2 = dot(delta.xyz, delta.xyz);
                        if (r2 < MAX_CUTOFF*MAX_CUTOFF) {
                            real invR = RSQRT(r2);
                            real r = r2*invR;
                            unsigned int atom2 = j;
                            LOAD_ATOM2_PARAMETERS
                            atom2 = atomIndices[j];
#ifdef USE_SYMMETRIC
                            real dEdR = 0;
#else
                            real4 dEdR1 = (real4) 0;
                            real4 dEdR2 = (real4) 0;
#endif
#ifdef USE_EXCLUSIONS
                            bool isExcluded = (atom1 >= NUM_ATOMS || atom2 >= NUM_ATOMS);
#endif
                            real tempEnergy = 0;
                            const real interactionScale = 1.0f;
                            COMPUTE_INTERACTION
                            energy += tempEnergy;
#ifdef USE_SYMMETRIC
                            delta.xyz *= dEdR;
                            force.xyz -= delta.xyz;
                            localData[j].fx += delta.x;
                            localData[j].fy += delta.y;
                            localData[j].fz += delta.z;
#else
                            force.xyz -= dEdR1.xyz;
                            localData[j].fx += dEdR2.x;
                            localData[j].fy += dEdR2.y;
                            localData[j].fz += dEdR2.z;
#endif
                        }
                    }

                   // Write results for atom1.

#ifdef SUPPORTS_64_BIT_ATOMICS
                    atom_add(&forceBuffers[atom1], (long) (force.x*0x100000000));
                    atom_add(&forceBuffers[atom1+PADDED_NUM_ATOMS], (long) (force.y*0x100000000));
                    atom_add(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], (long) (force.z*0x100000000));
#else
                    unsigned int offset = atom1 + get_group_id(0)*PADDED_NUM_ATOMS;
                    forceBuffers[offset].xyz = forceBuffers[offset].xyz+force.xyz;
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
                    real4 posq1 = posq[atom1];
                    LOAD_ATOM1_PARAMETERS
                    for (unsigned int j = 0; j < TILE_SIZE; j++) {
                        real4 posq2 = (real4) (localData[j].x, localData[j].y, localData[j].z, localData[j].q);
                        real4 delta = (real4) (posq2.xyz - posq1.xyz, 0);
#ifdef USE_PERIODIC
                        APPLY_PERIODIC_TO_DELTA(delta)
#endif
                        real r2 = dot(delta.xyz, delta.xyz);
#ifdef USE_CUTOFF
                        if (r2 < MAX_CUTOFF*MAX_CUTOFF) {
#endif
                            real invR = RSQRT(r2);
                            real r = r2*invR;
                            unsigned int atom2 = j;
                            LOAD_ATOM2_PARAMETERS
                            atom2 = atomIndices[j];
#ifdef USE_SYMMETRIC
                            real dEdR = 0;
#else
                            real4 dEdR1 = (real4) 0;
                            real4 dEdR2 = (real4) 0;
#endif
#ifdef USE_EXCLUSIONS
                            bool isExcluded = (atom1 >= NUM_ATOMS || atom2 >= NUM_ATOMS);
#endif
                            real tempEnergy = 0;
                            const real interactionScale = 1.0f;
                            COMPUTE_INTERACTION
                            energy += tempEnergy;
#ifdef USE_SYMMETRIC
                            delta.xyz *= dEdR;
                            force.xyz -= delta.xyz;
                            localData[j].fx += delta.x;
                            localData[j].fy += delta.y;
                            localData[j].fz += delta.z;
#else
                            force.xyz -= dEdR1.xyz;
                            localData[j].fx += dEdR2.x;
                            localData[j].fy += dEdR2.y;
                            localData[j].fz += dEdR2.z;
#endif
#ifdef USE_CUTOFF
                        }
#endif
                    }

                    // Write results for atom1.

#ifdef SUPPORTS_64_BIT_ATOMICS
                    atom_add(&forceBuffers[atom1], (long) (force.x*0x100000000));
                    atom_add(&forceBuffers[atom1+PADDED_NUM_ATOMS], (long) (force.y*0x100000000));
                    atom_add(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], (long) (force.z*0x100000000));
#else
                    unsigned int offset = atom1 + get_group_id(0)*PADDED_NUM_ATOMS;
                    forceBuffers[offset].xyz = forceBuffers[offset].xyz+force.xyz;
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
                    atom_add(&forceBuffers[atom2], (long) (localData[tgx].fx*0x100000000));
                    atom_add(&forceBuffers[atom2+PADDED_NUM_ATOMS], (long) (localData[tgx].fy*0x100000000));
                    atom_add(&forceBuffers[atom2+2*PADDED_NUM_ATOMS], (long) (localData[tgx].fz*0x100000000));
#else
                    unsigned int offset = atom2 + get_group_id(0)*PADDED_NUM_ATOMS;
                    real4 f = forceBuffers[offset];
                    f.x += localData[tgx].fx;
                    f.y += localData[tgx].fy;
                    f.z += localData[tgx].fz;
                    forceBuffers[offset] = f;
#endif
                }
            }
        }
        pos++;
    }
    energyBuffer[get_global_id(0)] += energy;
    SAVE_DERIVATIVES
}
