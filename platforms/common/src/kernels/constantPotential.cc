DEVICE real reduceValue(real value, LOCAL_ARG volatile real* temp) {
    const int thread = LOCAL_ID;
    SYNC_THREADS;
    temp[thread] = value;
    SYNC_THREADS;
    for (int step = 1; step < 32; step *= 2) {
        if (thread + step < LOCAL_SIZE && thread % (2 * step) == 0) {
            temp[thread] = temp[thread] + temp[thread + step];
        }
        SYNC_WARPS;
    }
    for (int step = 32; step < LOCAL_SIZE; step *= 2) {
        if (thread + step < LOCAL_SIZE && thread % (2 * step) == 0) {
            temp[thread] = temp[thread] + temp[thread + step];
        }
        SYNC_THREADS;
    }
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
    GLOBAL real* RESTRICT output
) {
    // This kernel expects to be executed in a single thread block.

    LOCAL volatile real temp[THREAD_BLOCK_SIZE];

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
        output[0] = totalCharge;
    }
}

KERNEL void evaluateSelfEnergyForces(
    GLOBAL real4* RESTRICT posq,
    GLOBAL real* RESTRICT charges,
    GLOBAL int* RESTRICT sysElec,
    GLOBAL real4* RESTRICT electrodeParams,
    GLOBAL int4* RESTRICT posCellOffsets,
    real4 periodicBoxVecX,
    real4 periodicBoxVecY,
    real4 periodicBoxVecZ,
    real4 externalField,
    GLOBAL mixed* RESTRICT energyBuffer,
    GLOBAL mm_ulong* RESTRICT forceBuffers

) {
    for (int i = GLOBAL_ID; i < NUM_PARTICLES; i += GLOBAL_SIZE) {
        const real4 pos = posq[i];
        const int4 offset = posCellOffsets[i];
#ifdef USE_POSQ_CHARGES
        const real charge = posq[i].w;
#else
        const real charge = charges[i];
#endif
        const real4 params = electrodeParams[sysElec[i] + 1];

        const real4 posOffset = pos + offset.x * periodicBoxVecX + offset.y * periodicBoxVecY + offset.z * periodicBoxVecZ;
        const real fieldTerm = posOffset.x * externalField.x + posOffset.y * externalField.y + posOffset.z * externalField.z;
        energyBuffer[GLOBAL_ID] += charge * (charge * params.w - params.x - fieldTerm);
        forceBuffers[i] += (mm_ulong) realToFixedPoint(charge * externalField.x);
        forceBuffers[i + PADDED_NUM_ATOMS] += (mm_ulong) realToFixedPoint(charge * externalField.y);
        forceBuffers[i + 2 * PADDED_NUM_ATOMS] += (mm_ulong) realToFixedPoint(charge * externalField.z);
    }
}

DEVICE void evaluateDirectDerivativesPair(GLOBAL mm_ulong* RESTRICT chargeDerivativesFixed, int ii, int jj, real3 delta, real charge1, real charge2, real width1, real width2) {
    if (ii == -1 && jj == -1) {
        return;
    }

    const real r2 = delta.x * delta.x + delta.y * delta.y + delta.z * delta.z;
    if (r2 >= CUTOFF_SQUARED) {
        return;
    }
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
    const real derivative = ONE_4PI_EPS0 * (erfcAlphaR - erfcEtaR) / r;
    if (ii != -1) {
        ATOMIC_ADD(&chargeDerivativesFixed[ii], (mm_ulong) realToFixedPoint(charge2 * derivative));
    }
    if (jj != -1) {
        ATOMIC_ADD(&chargeDerivativesFixed[jj], (mm_ulong) realToFixedPoint(charge1 * derivative));
    }
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
    int numExclusionTiles,
    GLOBAL const int* RESTRICT tiles,
    GLOBAL const unsigned int* RESTRICT interactionCount,
    GLOBAL const real4* RESTRICT blockCenter,
    GLOBAL const real4* RESTRICT blockSize,
    GLOBAL const int* RESTRICT interactingAtoms,
    GLOBAL int* RESTRICT tileCounter
) {
    if (GLOBAL_ID == 0) {
        *tileCounter = 0;
    }
    SYNC_THREADS;

    const int totalWarps = GLOBAL_SIZE/TILE_SIZE;
    const int warp = GLOBAL_ID/TILE_SIZE;
    const int tgx = LOCAL_ID & (TILE_SIZE-1);
    const int tbx = LOCAL_ID - tgx;

    LOCAL volatile int nextTile[WORK_GROUP_SIZE/TILE_SIZE];
    LOCAL volatile int atomIndices[WORK_GROUP_SIZE];
    LOCAL volatile real3 localPos[WORK_GROUP_SIZE];
    LOCAL volatile real localCharge[WORK_GROUP_SIZE];
    LOCAL volatile int localIndexElec[WORK_GROUP_SIZE];
    LOCAL volatile real localWidth[WORK_GROUP_SIZE];

    // First loop: process fixed tiles (ones that contain exclusions).  None of
    // the pairs involving electrode particles should actually be excluded, so
    // no special handling should be required for them.

    for (int tile = warp; tile < numExclusionTiles; tile += totalWarps) {
        const int2 tileIndices = exclusionTiles[tile];
        const unsigned int x = tileIndices.x;
        const unsigned int y = tileIndices.y;

        int i = x*TILE_SIZE + tgx;
        real4 posq1 = posq[i];
        real3 pos1 = trimTo3(posq1);
#ifdef USE_POSQ_CHARGES
        real charge1 = posq1.w;
#else
        real charge1 = charges[i];
#endif
        int ii = sysToElec[i];
        real width1 = electrodeParams[sysElec[i] + 1].y;

        if (x == y) {
            localPos[LOCAL_ID] = pos1;
            localCharge[LOCAL_ID] = charge1;
            localIndexElec[LOCAL_ID] = ii;
            localWidth[LOCAL_ID] = width1;
        }
        else {
            int j = y*TILE_SIZE + tgx;
            real4 posq2 = posq[j];
            localPos[LOCAL_ID] = trimTo3(posq2);
    #ifdef USE_POSQ_CHARGES
            localCharge[LOCAL_ID]  = posq2.w;
    #else
            localCharge[LOCAL_ID] = charges[j];
    #endif
            localIndexElec[LOCAL_ID] = sysToElec[j];
            localWidth[LOCAL_ID] = electrodeParams[sysElec[j] + 1].y;
        }
        SYNC_WARPS;

        if (i < NUM_PARTICLES) {
            for (int k = 0; k < TILE_SIZE; k++) {
#ifdef INTEL_WORKAROUND
                // Workaround for bug in Intel's OpenCL for CPUs.
                MEM_FENCE;
#endif
                int j = y*TILE_SIZE+k;
                if ((x != y || i < j) && j < NUM_PARTICLES) {
                    real3 pos2 = localPos[tbx + k];
                    real3 delta = make_real3(pos2.x-pos1.x, pos2.y-pos1.y, pos2.z-pos1.z);
                    APPLY_PERIODIC_TO_DELTA(delta)
                    evaluateDirectDerivativesPair(chargeDerivativesFixed, ii, localIndexElec[tbx + k], delta, charge1, localCharge[tbx + k], width1, localWidth[tbx + k]);
                }
            }
        }
#ifdef NVIDIA_WORKAROUND
        SYNC_THREADS;
#else
       SYNC_WARPS;
#endif
    }

    // Second loop: process tiles from the neighbor list.

    unsigned int numTiles = interactionCount[0];
    for (int tile = warp; tile < numTiles; tile += totalWarps) {
        if (tgx == 0)
            nextTile[tbx/TILE_SIZE] = ATOMIC_ADD(tileCounter, 1);
        SYNC_WARPS;
        int tileIndex = nextTile[tbx/TILE_SIZE];
        int x = tiles[tileIndex];
        real4 blockSizeX = blockSize[x];
        bool singlePeriodicCopy = (0.5f*periodicBoxSize.x-blockSizeX.x >= CUTOFF &&
                                   0.5f*periodicBoxSize.y-blockSizeX.y >= CUTOFF &&
                                   0.5f*periodicBoxSize.z-blockSizeX.z >= CUTOFF);

        int i = x*TILE_SIZE + tgx;
        real4 posq1 = posq[i];
        real3 pos1 = trimTo3(posq1);
#ifdef USE_POSQ_CHARGES
        real charge1 = posq1.w;
#else
        real charge1 = charges[i];
#endif
        int ii = sysToElec[i];
        real width1 = electrodeParams[sysElec[i] + 1].y;

        int j = interactingAtoms[tileIndex*TILE_SIZE+tgx];
        real4 posq2 = posq[j];
        real3 pos2 = trimTo3(posq2);
#ifdef USE_POSQ_CHARGES
        real charge2 = posq2.w;
#else
        real charge2 = charges[j];
#endif
        int jj = sysToElec[j];
        real width2 = electrodeParams[sysElec[j] + 1].y;

        atomIndices[LOCAL_ID] = j;
        localPos[LOCAL_ID] = pos2;
        localCharge[LOCAL_ID] = charge2;
        localIndexElec[LOCAL_ID] = jj;
        localWidth[LOCAL_ID] = width2;

        if (singlePeriodicCopy) {
            // The box is small enough that we can just translate all the atoms into a single periodic
            // box, then skip having to apply periodic boundary conditions later.

            real4 blockCenterX = blockCenter[x];
            APPLY_PERIODIC_TO_POS_WITH_CENTER(pos1, blockCenterX)
            APPLY_PERIODIC_TO_POS_WITH_CENTER(localPos[LOCAL_ID], blockCenterX)
            SYNC_WARPS;
            if (i < NUM_PARTICLES) {
                for (int k = 0; k < TILE_SIZE; k++) {
#ifdef INTEL_WORKAROUND
                    // Workaround for bug in Intel's OpenCL for CPUs.
                    MEM_FENCE;
#endif
                    if (atomIndices[tbx + k] < NUM_PARTICLES) {
                        pos2 = localPos[tbx + k];
                        real3 delta = make_real3(pos2.x-pos1.x, pos2.y-pos1.y, pos2.z-pos1.z);
                        evaluateDirectDerivativesPair(chargeDerivativesFixed, ii, localIndexElec[tbx + k], delta, charge1, localCharge[tbx + k], width1, localWidth[tbx + k]);
                    }
                }
            }
        }
        else
        {
            // We need to apply periodic boundary conditions separately for each interaction.

            SYNC_WARPS;
            if (i < NUM_PARTICLES) {
                for (int k = 0; k < TILE_SIZE; k++) {
#ifdef INTEL_WORKAROUND
                    // Workaround for bug in Intel's OpenCL for CPUs.
                    MEM_FENCE;
#endif
                    if (atomIndices[tbx + k] < NUM_PARTICLES) {
                        pos2 = localPos[tbx + k];
                        real3 delta = make_real3(pos2.x-pos1.x, pos2.y-pos1.y, pos2.z-pos1.z);
                        APPLY_PERIODIC_TO_DELTA(delta)
                        evaluateDirectDerivativesPair(chargeDerivativesFixed, ii, localIndexElec[tbx + k], delta, charge1, localCharge[tbx + k], width1, localWidth[tbx + k]);
                    }
                }
            }
        }
        SYNC_WARPS;
    }
}

KERNEL void finishDerivatives(
    GLOBAL real4* RESTRICT posq,
    GLOBAL real* RESTRICT charges,
    GLOBAL int* RESTRICT elecToSys,
    GLOBAL int* RESTRICT elecElec,
    GLOBAL real4* RESTRICT electrodeParams,
    real plasmaScale,
    GLOBAL int4* RESTRICT posCellOffsets,
    real4 periodicBoxVecX,
    real4 periodicBoxVecY,
    real4 periodicBoxVecZ,
    real4 externalField,
    GLOBAL real* RESTRICT chargeDerivatives,
    GLOBAL mm_long* RESTRICT chargeDerivativesFixed
) {
    const real fixed_scale = 1 / (real) 0x100000000;

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

        const real4 posOffset = pos + offset.x * periodicBoxVecX + offset.y * periodicBoxVecY + offset.z * periodicBoxVecZ;
        const real fieldTerm = posOffset.x * externalField.x + posOffset.y * externalField.y + posOffset.z * externalField.z;
        chargeDerivatives[ii] += 2 * (charge * params.w - plasmaScale) - params.x - fieldTerm + fixed_scale * chargeDerivativesFixed[ii];
    }
}
