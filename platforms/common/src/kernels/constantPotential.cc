#define WARP_SIZE 32

#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ >= 700
    #define WARP_SHUFFLE_DOWN(local, offset) __shfl_down_sync(0xffffffff, local, offset)
#elif defined(USE_HIP)
    #define WARP_SHUFFLE_DOWN(local, offset) __shfl_down(local, offset)
#endif

#ifdef WARP_SHUFFLE_DOWN
    #define TEMP_SIZE WARP_SIZE
#else
    #define TEMP_SIZE THREAD_BLOCK_SIZE
#endif

#ifdef USE_HIP
    #define ALIGN alignas(16)
#else
    #define ALIGN
#endif

typedef struct ALIGN {
    real x, y, z, q, width, derivative;
    int ii;
} AtomData;

DEVICE real reduceValue(real value, LOCAL_ARG volatile real* temp) {
    const int thread = LOCAL_ID;
    SYNC_THREADS;
#ifdef WARP_SHUFFLE_DOWN
    const int warpCount = LOCAL_SIZE / WARP_SIZE;
    const int warp = thread / WARP_SIZE;
    const int lane = thread % WARP_SIZE;
    for (int step = WARP_SIZE / 2; step > 0; step >>= 1) {
        value += WARP_SHUFFLE_DOWN(value, step);
    }
    if (!lane) {
        temp[warp] = value;
    }
    SYNC_THREADS;
    if (!warp) {
        value = lane < warpCount ? temp[lane] : 0;
        for (int step = WARP_SIZE / 2; step > 0; step >>= 1) {
            value += WARP_SHUFFLE_DOWN(value, step);
        }
        if (!lane) {
            temp[0] = value;
        }
    }
    SYNC_THREADS;
#else
    temp[thread] = value;
    SYNC_THREADS;
    for (int step = 1; step < WARP_SIZE / 2; step <<= 1) {
        if(thread + step < LOCAL_SIZE && thread % (2 * step) == 0) {
            temp[thread] += temp[thread + step];
        }
        SYNC_WARPS;
    }
    for (int step = WARP_SIZE / 2; step < LOCAL_SIZE; step <<= 1) {
        if(thread + step < LOCAL_SIZE && thread % (2 * step) == 0) {
            temp[thread] += temp[thread + step];
        }
        SYNC_THREADS;
    }
#endif
    return temp[0];
}

KERNEL void updateNonElectrodeCharges(
    GLOBAL real4* RESTRICT posq,
    GLOBAL real* RESTRICT charges,
    GLOBAL real* RESTRICT nonElectrodeCharges,
    GLOBAL int* RESTRICT sysElec
) {
    for (int i = GLOBAL_ID; i < NUM_PARTICLES; i += GLOBAL_SIZE) {
        if (sysElec[i] == -1) {
#ifdef USE_POSQ_CHARGES
            posq[i].w = nonElectrodeCharges[i];
#else
            charges[i] = nonElectrodeCharges[i];
#endif
        }
    }
}

KERNEL void updateElectrodeCharges(
    GLOBAL real4* RESTRICT posq,
    GLOBAL real* RESTRICT charges,
    GLOBAL real* RESTRICT electrodeCharges,
    GLOBAL int* RESTRICT elecToSys
) {
    for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
#ifdef USE_POSQ_CHARGES
        posq[elecToSys[ii]].w = electrodeCharges[ii];
#else
        charges[elecToSys[ii]] = electrodeCharges[ii];
#endif
    }
}

KERNEL void getTotalCharge(
    GLOBAL real4* RESTRICT posq,
    GLOBAL real* RESTRICT charges,
    GLOBAL real* RESTRICT totalChargeResult
) {
    // This kernel expects to be executed in a single thread block.

    LOCAL volatile real temp[TEMP_SIZE];

    real totalCharge = 0;
    for (int i = LOCAL_ID; i < NUM_PARTICLES; i += LOCAL_SIZE) {
#ifdef USE_POSQ_CHARGES
        totalCharge += posq[i].w;
#else
        totalCharge += charges[i];
#endif
    }
    totalCharge = reduceValue(totalCharge, temp);

    if (LOCAL_ID == 0) {
        totalChargeResult[0] = totalCharge;
    }
}

KERNEL void evaluateSelfEnergyForces(
    GLOBAL real4* RESTRICT posq,
    GLOBAL real* RESTRICT charges,
    GLOBAL int* RESTRICT sysElec,
    GLOBAL real4* RESTRICT electrodeParams,
    GLOBAL real* RESTRICT totalCharge,
    GLOBAL int4* RESTRICT posCellOffsets,
    real4 periodicBoxVecX,
    real4 periodicBoxVecY,
    real4 periodicBoxVecZ,
    real4 externalField,
    GLOBAL mixed* RESTRICT energyBuffer,
    GLOBAL mm_ulong* RESTRICT forceBuffers

) {
    if (GLOBAL_ID == 0) {
        energyBuffer[0] -= PLASMA_SCALE * totalCharge[0] * totalCharge[0] / (periodicBoxVecX.x * periodicBoxVecY.y * periodicBoxVecZ.z * EWALD_ALPHA * EWALD_ALPHA);
    }

    for (int i = GLOBAL_ID; i < NUM_PARTICLES; i += GLOBAL_SIZE) {
        const real4 pos = posq[i];
        const int4 offset = posCellOffsets[i];
#ifdef USE_POSQ_CHARGES
        const real charge = posq[i].w;
#else
        const real charge = charges[i];
#endif
        const real4 params = electrodeParams[sysElec[i] + 1];

        const real4 posOffset = pos - offset.x * periodicBoxVecX - offset.y * periodicBoxVecY - offset.z * periodicBoxVecZ;
        const real fieldTerm = posOffset.x * externalField.x + posOffset.y * externalField.y + posOffset.z * externalField.z;
        energyBuffer[GLOBAL_ID] += charge * (charge * params.w - params.x - fieldTerm);
        forceBuffers[i] += (mm_ulong) realToFixedPoint(charge * externalField.x);
        forceBuffers[i + PADDED_NUM_ATOMS] += (mm_ulong) realToFixedPoint(charge * externalField.y);
        forceBuffers[i + 2 * PADDED_NUM_ATOMS] += (mm_ulong) realToFixedPoint(charge * externalField.z);
    }
}

DEVICE real evaluateDirectDerivative(real r2, real width1, real width2) {
    const real r = SQRT(r2);
    const real alphaR = EWALD_ALPHA * r;
    const real etaR = r / SQRT(width1 * width1 + width2 * width2);
#ifdef USE_DOUBLE_PRECISION
    const real erfcAlphaR = erfc(alphaR);
    const real erfcEtaR = erfc(etaR);
#else
    const real expAlphaRSqr = EXP(-alphaR * alphaR);
    const real expEtaRSqr = EXP(-etaR * etaR);
    const real tAlpha = RECIP(1.0f+0.3275911f*alphaR);
    const real tEta = RECIP(1.0f+0.3275911f*etaR);
    const real erfcAlphaR = (0.254829592f+(-0.284496736f+(1.421413741f+(-1.453152027f+1.061405429f*tAlpha)*tAlpha)*tAlpha)*tAlpha)*tAlpha*expAlphaRSqr;
    const real erfcEtaR = (0.254829592f+(-0.284496736f+(1.421413741f+(-1.453152027f+1.061405429f*tEta)*tEta)*tEta)*tEta)*tEta*expEtaRSqr;
#endif
    return ONE_4PI_EPS0 * (erfcAlphaR - erfcEtaR) / r;
}

KERNEL void evaluateDirectDerivatives(
    GLOBAL const real4* RESTRICT posq,
    GLOBAL real* RESTRICT charges,
    GLOBAL int* RESTRICT sysToElec,
    GLOBAL int* RESTRICT sysElec,
    GLOBAL real4* RESTRICT electrodeParams,
    GLOBAL mm_ulong* RESTRICT chargeDerivativesFixed,
    real4 periodicBoxSize,
    real4 invPeriodicBoxSize,
    real4 periodicBoxVecX,
    real4 periodicBoxVecY,
    real4 periodicBoxVecZ,
    GLOBAL const int2* RESTRICT exclusionTiles,
    GLOBAL const int* RESTRICT tiles,
    GLOBAL const unsigned int* RESTRICT interactionCount,
    GLOBAL const real4* RESTRICT blockCenter,
    GLOBAL const real4* RESTRICT blockSize,
    GLOBAL const int* RESTRICT interactingAtoms,
    unsigned int maxTiles
) {
#ifndef DEVICE_IS_CPU
    // GPU-specific direct derivative calculation kernel.

    const unsigned int totalWarps = GLOBAL_SIZE / TILE_SIZE;
    const unsigned int warp = GLOBAL_ID / TILE_SIZE;
    const unsigned int tgx = LOCAL_ID & (TILE_SIZE - 1);
    const unsigned int tbx = LOCAL_ID - tgx;
    LOCAL AtomData localData[WORK_GROUP_SIZE];

    // First loop: process tiles that contain exclusions.  Exclusions cannot
    // involve electrode particles, so only self-interactions need be excluded.

    const unsigned int firstExclusionTile = warp * NUM_EXCLUSION_TILES / totalWarps;
    const unsigned int lastExclusionTile = (warp + 1) * NUM_EXCLUSION_TILES / totalWarps;
    for (int pos = firstExclusionTile; pos < lastExclusionTile; pos++) {
        const int2 tileIndices = exclusionTiles[pos];
        const unsigned int x = tileIndices.x;
        const unsigned int y = tileIndices.y;

        const unsigned int atom1 = x * TILE_SIZE + tgx;
        const real4 posq1 = posq[atom1];
#ifdef USE_POSQ_CHARGES
        const real charge1 = posq1.w;
#else
        const real charge1 = charges[atom1];
#endif
        const real width1 = electrodeParams[sysElec[atom1] + 1].y;
        const int ii1 = sysToElec[atom1];

        real derivative = 0;

        if (x == y) {
            // This tile is on the diagonal.

            localData[LOCAL_ID].x = posq1.x;
            localData[LOCAL_ID].y = posq1.y;
            localData[LOCAL_ID].z = posq1.z;
            localData[LOCAL_ID].q = charge1;
            localData[LOCAL_ID].width = width1;
            SYNC_WARPS;

            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                if (ii1 != -1 && atom1 < NUM_PARTICLES && y * TILE_SIZE + j < NUM_PARTICLES && atom1 != y * TILE_SIZE + j) {
                    const real3 pos2 = make_real3(localData[tbx + j].x, localData[tbx + j].y, localData[tbx + j].z);
                    real3 delta = make_real3(pos2.x - posq1.x, pos2.y - posq1.y, pos2.z - posq1.z);
                    APPLY_PERIODIC_TO_DELTA(delta)
                    const real r2 = delta.x * delta.x + delta.y * delta.y + delta.z * delta.z;
                    if (r2 < CUTOFF_SQUARED) {
                        derivative += localData[tbx + j].q * evaluateDirectDerivative(r2, width1, localData[tbx + j].width);
                    }
                }
                SYNC_WARPS;
            }
        }
        else {
            // This is an off-diagonal tile.

            unsigned int j = y * TILE_SIZE + tgx;
            const real4 posq2 = posq[j];
            localData[LOCAL_ID].x = posq2.x;
            localData[LOCAL_ID].y = posq2.y;
            localData[LOCAL_ID].z = posq2.z;
#ifdef USE_POSQ_CHARGES
            localData[LOCAL_ID].q = posq2.w;
#else
            localData[LOCAL_ID].q = charges[j];
#endif
            localData[LOCAL_ID].width = electrodeParams[sysElec[j] + 1].y;
            localData[LOCAL_ID].derivative = 0;
            localData[LOCAL_ID].ii = sysToElec[j];
            SYNC_WARPS;

            unsigned int tj = tgx;
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                if (atom1 < NUM_PARTICLES && y * TILE_SIZE + tj < NUM_PARTICLES && (ii1 != -1 || localData[tbx + tj].ii != -1)) {
                    const real3 pos2 = make_real3(localData[tbx + tj].x, localData[tbx + tj].y, localData[tbx + tj].z);
                    real3 delta = make_real3(pos2.x - posq1.x, pos2.y - posq1.y, pos2.z - posq1.z);
                    APPLY_PERIODIC_TO_DELTA(delta)
                    const real r2 = delta.x * delta.x + delta.y * delta.y + delta.z * delta.z;
                    if (r2 < CUTOFF_SQUARED) {
                        const real derivativeScale = evaluateDirectDerivative(r2, width1, localData[tbx + tj].width);
                        derivative += localData[tbx + tj].q * derivativeScale;
                        localData[tbx + tj].derivative += charge1 * derivativeScale;
                    }
                }
                tj = (tj + 1) & (TILE_SIZE - 1);
                SYNC_WARPS;
            }
        }

        // Write results.
        if (ii1 != -1) {
            ATOMIC_ADD(&chargeDerivativesFixed[ii1], (mm_ulong) realToFixedPoint(derivative));
        }
        if (x != y) {
            const int ii2 = localData[LOCAL_ID].ii;
            if (ii2 != -1) {
                ATOMIC_ADD(&chargeDerivativesFixed[ii2], (mm_ulong) realToFixedPoint(localData[LOCAL_ID].derivative));
            }
        }
    }

    // Second loop: process tiles without exclusions.

    const unsigned int numTiles = interactionCount[0];
    if (numTiles > maxTiles) {
        // There wasn't enough memory for the neighbor list.
        return;
    }
    int pos = (int) (warp * ((mm_long) numTiles) / totalWarps);
    const int end = (int) ((warp + 1) * ((mm_long) numTiles) / totalWarps);
    LOCAL int atomIndices[WORK_GROUP_SIZE];

    while (pos < end) {
        real derivative = 0;

        // Extract the coordinates of this tile.
        const int x = tiles[pos];
        const real4 blockSizeX = blockSize[x];
        const bool singlePeriodicCopy = (0.5f * periodicBoxSize.x - blockSizeX.x >= CUTOFF &&
                                         0.5f * periodicBoxSize.y - blockSizeX.y >= CUTOFF &&
                                         0.5f * periodicBoxSize.z - blockSizeX.z >= CUTOFF);

        const unsigned int atom1 = x * TILE_SIZE + tgx;

        // Load atom data for this tile.
        real4 posq1 = posq[atom1];
#ifdef USE_POSQ_CHARGES
        const real charge1 = posq1.w;
#else
        const real charge1 = charges[atom1];
#endif
        const real width1 = electrodeParams[sysElec[atom1] + 1].y;
        const int ii1 = sysToElec[atom1];

        unsigned int j = interactingAtoms[pos * TILE_SIZE + tgx];
        atomIndices[LOCAL_ID] = j;
        if (j < PADDED_NUM_ATOMS) {
            const real4 posq2 = posq[j];
            localData[LOCAL_ID].x = posq2.x;
            localData[LOCAL_ID].y = posq2.y;
            localData[LOCAL_ID].z = posq2.z;
#ifdef USE_POSQ_CHARGES
            localData[LOCAL_ID].q = posq2.w;
#else
            localData[LOCAL_ID].q = charges[j];
#endif
            localData[LOCAL_ID].width = electrodeParams[sysElec[j] + 1].y;
            localData[LOCAL_ID].derivative = 0;
            localData[LOCAL_ID].ii = sysToElec[j];
        }
        SYNC_WARPS;

        if (singlePeriodicCopy) {
            // The box is small enough; we can translate atoms and avoid having
            // to apply periodic boundary conditions to every interaction later.

            real4 blockCenterX = blockCenter[x];
            APPLY_PERIODIC_TO_POS_WITH_CENTER(posq1, blockCenterX)
            APPLY_PERIODIC_TO_POS_WITH_CENTER(localData[LOCAL_ID], blockCenterX)
            SYNC_WARPS;

            unsigned int tj = tgx;
            for (j = 0; j < TILE_SIZE; j++) {
                const int atom2 = atomIndices[tbx + tj];
                if (atom1 < NUM_PARTICLES && atom2 < NUM_PARTICLES && (ii1 != -1 || localData[tbx + tj].ii != -1)) {
                    const real3 pos2 = make_real3(localData[tbx + tj].x, localData[tbx + tj].y, localData[tbx + tj].z);
                    real3 delta = make_real3(pos2.x - posq1.x, pos2.y - posq1.y, pos2.z - posq1.z);
                    const real r2 = delta.x * delta.x + delta.y * delta.y + delta.z * delta.z;
                    if (r2 < CUTOFF_SQUARED) {
                        const real derivativeScale = evaluateDirectDerivative(r2, width1, localData[tbx + tj].width);
                        derivative += localData[tbx + tj].q * derivativeScale;
                        localData[tbx + tj].derivative += charge1 * derivativeScale;
                    }
                }
                tj = (tj + 1) & (TILE_SIZE - 1);
                SYNC_WARPS;
            }
        }
        else {
            // We must apply periodic boundary conditions to every interaction.

            unsigned int tj = tgx;
            for (j = 0; j < TILE_SIZE; j++) {
                const int atom2 = atomIndices[tbx + tj];
                if (atom1 < NUM_PARTICLES && atom2 < NUM_PARTICLES && (ii1 != -1 || localData[tbx + tj].ii != -1)) {
                    const real3 pos2 = make_real3(localData[tbx + tj].x, localData[tbx + tj].y, localData[tbx + tj].z);
                    real3 delta = make_real3(pos2.x - posq1.x, pos2.y - posq1.y, pos2.z - posq1.z);
                    APPLY_PERIODIC_TO_DELTA(delta)
                    const real r2 = delta.x * delta.x + delta.y * delta.y + delta.z * delta.z;
                    if (r2 < CUTOFF_SQUARED) {
                        const real derivativeScale = evaluateDirectDerivative(r2, width1, localData[tbx + tj].width);
                        derivative += localData[tbx + tj].q * derivativeScale;
                        localData[tbx + tj].derivative += charge1 * derivativeScale;
                    }
                }
                tj = (tj + 1) & (TILE_SIZE - 1);
                SYNC_WARPS;
            }
        }

        // Write results.
        if (ii1 != -1) {
            ATOMIC_ADD(&chargeDerivativesFixed[ii1], (mm_ulong) realToFixedPoint(derivative));
        }
        if (atomIndices[LOCAL_ID] < PADDED_NUM_ATOMS) {
            const int ii2 = localData[LOCAL_ID].ii;
            if (ii2 != -1) {
                ATOMIC_ADD(&chargeDerivativesFixed[ii2], (mm_ulong) realToFixedPoint(localData[LOCAL_ID].derivative));
            }
        }
        pos++;
    }
#else
    // CPU-specific direct derivative calculation kernel.

    LOCAL AtomData localData[TILE_SIZE];

    // First loop: process tiles that contain exclusions.  Exclusions cannot
    // involve electrode particles, so only self-interactions need be excluded.

    const unsigned int firstExclusionTile = GROUP_ID * NUM_EXCLUSION_TILES / NUM_GROUPS;
    const unsigned int lastExclusionTile = (GROUP_ID + 1) * NUM_EXCLUSION_TILES / NUM_GROUPS;
    for (int pos = firstExclusionTile; pos < lastExclusionTile; pos++) {
        const int2 tileIndices = exclusionTiles[pos];
        const unsigned int x = tileIndices.x;
        const unsigned int y = tileIndices.y;

        // Load atom data for this tile.
        for (unsigned int tgx = 0; tgx < TILE_SIZE; tgx++) {
            const unsigned int j = y * TILE_SIZE + tgx;
            const real4 posq2 = posq[j];
            localData[tgx].x = posq2.x;
            localData[tgx].y = posq2.y;
            localData[tgx].z = posq2.z;
#ifdef USE_POSQ_CHARGES
            localData[tgx].q = posq2.w;
#else
            localData[tgx].q = charges[j];
#endif
            localData[tgx].width = electrodeParams[sysElec[j] + 1].y;
            localData[tgx].ii = sysToElec[j];
        }

        if (x == y) {
            // This tile is on the diagonal.

            for (unsigned int tgx = 0; tgx < TILE_SIZE; tgx++) {
                const unsigned int atom1 = x * TILE_SIZE + tgx;
                const int ii1 = sysToElec[atom1];
                if (ii1 == -1) {
                    continue;
                }
                const real4 posq1 = posq[atom1];
#ifdef USE_POSQ_CHARGES
                const real charge1 = posq1.w;
#else
                const real charge1 = charges[atom1];
#endif
                const real width1 = electrodeParams[sysElec[atom1] + 1].y;

                real derivative = 0;

                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    if (atom1 < NUM_PARTICLES && y * TILE_SIZE + j < NUM_PARTICLES && atom1 != y * TILE_SIZE + j) {
                        const real3 pos2 = make_real3(localData[j].x, localData[j].y, localData[j].z);
                        real3 delta = make_real3(pos2.x - posq1.x, pos2.y - posq1.y, pos2.z - posq1.z);
                        APPLY_PERIODIC_TO_DELTA(delta)
                        const real r2 = delta.x * delta.x + delta.y * delta.y + delta.z * delta.z;
                        if (r2 < CUTOFF_SQUARED) {
                            derivative += localData[j].q * evaluateDirectDerivative(r2, width1, localData[j].width);
                        }
                    }
                }

                // Write results.
                ATOMIC_ADD(&chargeDerivativesFixed[ii1], (mm_ulong) realToFixedPoint(derivative));
            }
        }
        else {
            // This is an off-diagonal tile.

            for (unsigned int tgx = 0; tgx < TILE_SIZE; tgx++) {
                localData[tgx].derivative = 0;
            }
            for (unsigned int tgx = 0; tgx < TILE_SIZE; tgx++) {
                const unsigned int atom1 = x * TILE_SIZE + tgx;
                const real4 posq1 = posq[atom1];
#ifdef USE_POSQ_CHARGES
                const real charge1 = posq1.w;
#else
                const real charge1 = charges[atom1];
#endif
                const real width1 = electrodeParams[sysElec[atom1] + 1].y;
                const int ii1 = sysToElec[atom1];

                real derivative = 0;

                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    if (atom1 < NUM_PARTICLES && y * TILE_SIZE + j < NUM_PARTICLES && (ii1 != -1 || localData[j].ii != -1)) {
                        const real3 pos2 = make_real3(localData[j].x, localData[j].y, localData[j].z);
                        real3 delta = make_real3(pos2.x - posq1.x, pos2.y - posq1.y, pos2.z - posq1.z);
                        APPLY_PERIODIC_TO_DELTA(delta)
                        const real r2 = delta.x * delta.x + delta.y * delta.y + delta.z * delta.z;
                        if (r2 < CUTOFF_SQUARED) {
                            const real derivativeScale = evaluateDirectDerivative(r2, width1, localData[j].width);
                            derivative += localData[j].q * derivativeScale;
                            localData[j].derivative += charge1 * derivativeScale;
                        }
                    }
                }

                // Write results for "1" atoms.
                if (ii1 != -1) {
                    ATOMIC_ADD(&chargeDerivativesFixed[ii1], (mm_ulong) realToFixedPoint(derivative));
                }
            }

            // Write results for "2" atoms.
            for (unsigned int tgx = 0; tgx < TILE_SIZE; tgx++) {
                const int ii2 = localData[tgx].ii;
                if (ii2 != -1) {
                    ATOMIC_ADD(&chargeDerivativesFixed[ii2], (mm_ulong) realToFixedPoint(localData[tgx].derivative));
                }
            }
        }
    }

    // Second loop: process tiles without exclusions.

    const unsigned int numTiles = interactionCount[0];
    if (numTiles > maxTiles) {
        // There wasn't enough memory for the neighbor list.
        return;
    }
    int pos = (int) (GROUP_ID * ((mm_long) numTiles) / NUM_GROUPS);
    const int end = (int) ((GROUP_ID + 1) * ((mm_long) numTiles) / NUM_GROUPS);
    LOCAL int atomIndices[TILE_SIZE];

    while (pos < end) {
        real derivative = 0;

        // Extract the coordinates of this tile.
        const int x = tiles[pos];
        const real4 blockSizeX = blockSize[x];
        const bool singlePeriodicCopy = (0.5f * periodicBoxSize.x - blockSizeX.x >= CUTOFF &&
                                         0.5f * periodicBoxSize.y - blockSizeX.y >= CUTOFF &&
                                         0.5f * periodicBoxSize.z - blockSizeX.z >= CUTOFF);

        // Load atom data for this tile.
        for (unsigned int tgx = 0; tgx < TILE_SIZE; tgx++) {
            const unsigned int j = interactingAtoms[pos * TILE_SIZE + tgx];
            atomIndices[tgx] = j;
            if (j < PADDED_NUM_ATOMS) {
                const real4 posq2 = posq[j];
                localData[tgx].x = posq2.x;
                localData[tgx].y = posq2.y;
                localData[tgx].z = posq2.z;
    #ifdef USE_POSQ_CHARGES
                localData[tgx].q = posq2.w;
    #else
                localData[tgx].q = charges[j];
    #endif
                localData[tgx].width = electrodeParams[sysElec[j] + 1].y;
                localData[tgx].derivative = 0;
                localData[tgx].ii = sysToElec[j];
            }
        }

        if (singlePeriodicCopy) {
            // The box is small enough; we can translate atoms and avoid having
            // to apply periodic boundary conditions to every interaction later.

            real4 blockCenterX = blockCenter[x];
            for (unsigned int tgx = 0; tgx < TILE_SIZE; tgx++) {
                APPLY_PERIODIC_TO_POS_WITH_CENTER(localData[tgx], blockCenterX)
            }
            for (unsigned int tgx = 0; tgx < TILE_SIZE; tgx++) {
                const unsigned int atom1 = x * TILE_SIZE + tgx;
                real4 posq1 = posq[atom1];
                APPLY_PERIODIC_TO_POS_WITH_CENTER(posq1, blockCenterX)
#ifdef USE_POSQ_CHARGES
                const real charge1 = posq1.w;
#else
                const real charge1 = charges[atom1];
#endif
                const real width1 = electrodeParams[sysElec[atom1] + 1].y;
                const int ii1 = sysToElec[atom1];

                real derivative = 0;

                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    const int atom2 = atomIndices[j];
                    if (atom1 < NUM_PARTICLES && atom2 < NUM_PARTICLES && (ii1 != -1 || localData[j].ii != -1)) {
                        const real3 pos2 = make_real3(localData[j].x, localData[j].y, localData[j].z);
                        real3 delta = make_real3(pos2.x - posq1.x, pos2.y - posq1.y, pos2.z - posq1.z);
                        const real r2 = delta.x * delta.x + delta.y * delta.y + delta.z * delta.z;
                        if (r2 < CUTOFF_SQUARED) {
                            const real derivativeScale = evaluateDirectDerivative(r2, width1, localData[j].width);
                            derivative += localData[j].q * derivativeScale;
                            localData[j].derivative += charge1 * derivativeScale;
                        }
                    }
                }

                // Write results for "1" atoms.
                if (ii1 != -1) {
                    ATOMIC_ADD(&chargeDerivativesFixed[ii1], (mm_ulong) realToFixedPoint(derivative));
                }
            }
        }
        else {
            // We must apply periodic boundary conditions to every interaction.

            for (unsigned int tgx = 0; tgx < TILE_SIZE; tgx++) {
                const unsigned int atom1 = x * TILE_SIZE + tgx;
                const real4 posq1 = posq[atom1];
#ifdef USE_POSQ_CHARGES
                const real charge1 = posq1.w;
#else
                const real charge1 = charges[atom1];
#endif
                const real width1 = electrodeParams[sysElec[atom1] + 1].y;
                const int ii1 = sysToElec[atom1];

                real derivative = 0;

                for (unsigned int j = 0; j < TILE_SIZE; j++) {
                    const int atom2 = atomIndices[j];
                    if (atom1 < NUM_PARTICLES && atom2 < NUM_PARTICLES && (ii1 != -1 || localData[j].ii != -1)) {
                        const real3 pos2 = make_real3(localData[j].x, localData[j].y, localData[j].z);
                        real3 delta = make_real3(pos2.x - posq1.x, pos2.y - posq1.y, pos2.z - posq1.z);
                        APPLY_PERIODIC_TO_DELTA(delta)
                        const real r2 = delta.x * delta.x + delta.y * delta.y + delta.z * delta.z;
                        if (r2 < CUTOFF_SQUARED) {
                            const real derivativeScale = evaluateDirectDerivative(r2, width1, localData[j].width);
                            derivative += localData[j].q * derivativeScale;
                            localData[j].derivative += charge1 * derivativeScale;
                        }
                    }
                }

                // Write results for "1" atoms.
                if (ii1 != -1) {
                    ATOMIC_ADD(&chargeDerivativesFixed[ii1], (mm_ulong) realToFixedPoint(derivative));
                }
            }
        }

        // Write results for "2" atoms.
        for (unsigned int tgx = 0; tgx < TILE_SIZE; tgx++) {
            if (atomIndices[tgx] < PADDED_NUM_ATOMS) {
                const int ii2 = localData[tgx].ii;
                if (ii2 != -1) {
                    ATOMIC_ADD(&chargeDerivativesFixed[ii2], (mm_ulong) realToFixedPoint(localData[tgx].derivative));
                }
            }
        }
        pos++;
    }
#endif
}

KERNEL void finishDerivatives(
    GLOBAL real4* RESTRICT posq,
    GLOBAL real* RESTRICT charges,
    GLOBAL int* RESTRICT elecToSys,
    GLOBAL int* RESTRICT elecElec,
    GLOBAL real4* RESTRICT electrodeParams,
    GLOBAL real* RESTRICT totalCharge,
    GLOBAL int4* RESTRICT posCellOffsets,
    real4 periodicBoxVecX,
    real4 periodicBoxVecY,
    real4 periodicBoxVecZ,
    real4 externalField,
    GLOBAL real* RESTRICT chargeDerivatives,
    GLOBAL mm_long* RESTRICT chargeDerivativesFixed
) {
    const real fixedScale = 1 / (real) 0x100000000;
    const real plasmaScale = PLASMA_SCALE * totalCharge[0] / (periodicBoxVecX.x * periodicBoxVecY.y * periodicBoxVecZ.z * EWALD_ALPHA * EWALD_ALPHA);

    for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
        int i = elecToSys[ii];
        const real4 pos = posq[i];
        const int4 offset = posCellOffsets[i];
#ifdef USE_POSQ_CHARGES
        const real charge = pos.w;
#else
        const real charge = charges[i];
#endif
        const real4 params = electrodeParams[elecElec[ii] + 1];

        const real4 posOffset = pos - offset.x * periodicBoxVecX - offset.y * periodicBoxVecY - offset.z * periodicBoxVecZ;
        const real fieldTerm = posOffset.x * externalField.x + posOffset.y * externalField.y + posOffset.z * externalField.z;
        chargeDerivatives[ii] += 2 * (charge * params.w - plasmaScale) - params.x - fieldTerm + fixedScale * chargeDerivativesFixed[ii];
    }
}
