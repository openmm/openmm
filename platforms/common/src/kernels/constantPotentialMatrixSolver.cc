#define WARP_SIZE 32

#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ >= 700
    #define WARP_SHUFFLE(local, index) __shfl_sync(0xffffffff, local, index)
    #define WARP_SHUFFLE_DOWN(local, offset) __shfl_down_sync(0xffffffff, local, offset)
#elif defined(USE_HIP)
    #define WARP_SHUFFLE(local, index) __shfl(local, index)
    #define WARP_SHUFFLE_DOWN(local, offset) __shfl_down(local, offset)
#endif

#ifdef WARP_SHUFFLE_DOWN
    #define TEMP_SIZE WARP_SIZE
#else
    #define TEMP_SIZE THREAD_BLOCK_SIZE
#endif

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

KERNEL void checkSavedElectrodePositions(GLOBAL real4* RESTRICT posq, GLOBAL real4* RESTRICT electrodePosData, GLOBAL int* RESTRICT elecToSys, GLOBAL int* RESTRICT result) {
    for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
        real4 posqPosition = posq[elecToSys[ii]];
        real4 savedPosition = electrodePosData[ii];
        if (posqPosition.x != savedPosition.x || posqPosition.y != savedPosition.y || posqPosition.z != savedPosition.z) {
            *result = 1;
            break;
        }
    }
}

KERNEL void saveElectrodePositions(GLOBAL real4* RESTRICT posq, GLOBAL real4* RESTRICT electrodePosData, GLOBAL int* RESTRICT elecToSys) {
    for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
        electrodePosData[ii] = posq[elecToSys[ii]];
    }
}

KERNEL void solve(GLOBAL real* RESTRICT electrodeCharges, GLOBAL real* RESTRICT chargeDerivatives, GLOBAL real* RESTRICT capacitance
#ifdef USE_CHARGE_CONSTRAINT
    , GLOBAL real* RESTRICT constraintVector, real chargeTarget
#endif
) {
    // This kernel expects to be executed in a single thread block.

#if CHUNK_SIZE > 1
    LOCAL volatile real chunkCharges[CHUNK_SIZE];
#endif

    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        electrodeCharges[ii] = -chargeDerivatives[ii];
    }
    SYNC_THREADS;

    // Cholesky solve step 1 (outer loop over chunks of rows).

    for (int jj = 0; jj < PADDED_PROBLEM_SIZE; jj += CHUNK_SIZE) {
        if (LOCAL_ID < CHUNK_SIZE) {
#if CHUNK_SIZE > 1
    #ifdef WARP_SHUFFLE
            real threadCharge = electrodeCharges[jj + LOCAL_ID];
            for (int k = 0; k < CHUNK_SIZE - 1; k++) {
                const real chargeShuffled = WARP_SHUFFLE(threadCharge, k);
                if (LOCAL_ID > k) {
                    threadCharge -= chargeShuffled * capacitance[(mm_long) (jj + k) * PADDED_PROBLEM_SIZE + (jj + LOCAL_ID)];
                }
            }
            SYNC_WARPS;
            electrodeCharges[jj + LOCAL_ID] = chunkCharges[LOCAL_ID] = threadCharge;
    #else
            chunkCharges[LOCAL_ID] = electrodeCharges[jj + LOCAL_ID];
            for (int k = 0; k < CHUNK_SIZE - 1; k++) {
                SYNC_WARPS;
                if (LOCAL_ID > k) {
                    chunkCharges[LOCAL_ID] -= chunkCharges[k] * capacitance[(mm_long) (jj + k) * PADDED_PROBLEM_SIZE + (jj + LOCAL_ID)];
                }
            }
            SYNC_WARPS;
            electrodeCharges[jj + LOCAL_ID] = chunkCharges[LOCAL_ID];
    #endif
#endif
        }
        SYNC_THREADS;
        for (int ii = jj + CHUNK_SIZE + LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
#if CHUNK_SIZE > 1
            real chargeOffset = 0;
            for (int k = 0; k < CHUNK_SIZE; k++) {
                chargeOffset += chunkCharges[k] * capacitance[(mm_long) (jj + k) * PADDED_PROBLEM_SIZE + ii];
            }
            electrodeCharges[ii] -= chargeOffset;
#else
            electrodeCharges[ii] -= electrodeCharges[jj] * capacitance[(mm_long) jj * PADDED_PROBLEM_SIZE + ii];
#endif
        }
        SYNC_THREADS;
    }

    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        electrodeCharges[ii] *= capacitance[(mm_long) ii * PADDED_PROBLEM_SIZE + ii];
    }
    SYNC_THREADS;

    // Cholesky solve step 2 (outer loop over chunks of columns).

    for (int jj = PADDED_PROBLEM_SIZE - CHUNK_SIZE; jj >= 0; jj -= CHUNK_SIZE) {
        if (LOCAL_ID < CHUNK_SIZE) {
#if CHUNK_SIZE > 1
    #ifdef WARP_SHUFFLE
            real threadCharge = electrodeCharges[jj + LOCAL_ID];
            for (int k = CHUNK_SIZE - 1; k >= 0; k--) {
                const real chargeShuffled = WARP_SHUFFLE(threadCharge, k);
                if (LOCAL_ID < k) {
                    threadCharge -= chargeShuffled * capacitance[(mm_long) (jj + k) * PADDED_PROBLEM_SIZE + (jj + LOCAL_ID)];
                }
            }
            SYNC_WARPS;
            electrodeCharges[jj + LOCAL_ID] = chunkCharges[LOCAL_ID] = threadCharge;
    #else
            chunkCharges[LOCAL_ID] = electrodeCharges[jj + LOCAL_ID];
            for (int k = CHUNK_SIZE - 1; k >= 0; k--) {
                SYNC_WARPS;
                if (LOCAL_ID < k) {
                    chunkCharges[LOCAL_ID] -= chunkCharges[k] * capacitance[(mm_long) (jj + k) * PADDED_PROBLEM_SIZE + (jj + LOCAL_ID)];
                }
            }
            SYNC_WARPS;
            electrodeCharges[jj + LOCAL_ID] = chunkCharges[LOCAL_ID];
    #endif
#endif
        }
        SYNC_THREADS;
        for (int ii = LOCAL_ID; ii < jj; ii += LOCAL_SIZE) {
#if CHUNK_SIZE > 1
            real chargeOffset = 0;
            for (int k = 0; k < CHUNK_SIZE; k++) {
                chargeOffset += chunkCharges[k] * capacitance[(mm_long) (jj + k) * PADDED_PROBLEM_SIZE + ii];
            }
            electrodeCharges[ii] -= chargeOffset;
#else
            electrodeCharges[ii] -= electrodeCharges[jj] * capacitance[(mm_long) jj * PADDED_PROBLEM_SIZE + ii];
#endif
        }
        SYNC_THREADS;
    }

    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        electrodeCharges[ii] *= capacitance[(mm_long) ii * PADDED_PROBLEM_SIZE + ii];
    }
    SYNC_THREADS;

#ifdef USE_CHARGE_CONSTRAINT
    LOCAL volatile real temp[TEMP_SIZE];

    real chargeOffset = 0;
    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        chargeOffset -= electrodeCharges[ii];
    }
    chargeOffset = chargeTarget + reduceValue(chargeOffset, temp);
    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        electrodeCharges[ii] += chargeOffset * constraintVector[ii];
    }
#endif
}
