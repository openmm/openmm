#define CHUNK_SIZE 32

DEVICE real reduceValue(real value, LOCAL_ARG volatile real* temp) {
    const int thread = LOCAL_ID;
    SYNC_THREADS;
    temp[thread] = value;
    SYNC_THREADS;
    for (int step = 1; step < 16; step *= 2) {
        if (thread + step < LOCAL_SIZE && thread % (2 * step) == 0) {
            temp[thread] = temp[thread] + temp[thread + step];
        }
        SYNC_WARPS;
    }
    for (int step = 16; step < LOCAL_SIZE; step *= 2) {
        if (thread + step < LOCAL_SIZE && thread % (2 * step) == 0) {
            temp[thread] = temp[thread] + temp[thread + step];
        }
        SYNC_THREADS;
    }
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

    const int NUM_CHUNKS = NUM_ELECTRODE_PARTICLES / CHUNK_SIZE;

    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        electrodeCharges[ii] = -chargeDerivatives[ii];
    }
    SYNC_THREADS;

    LOCAL volatile real chunkCharges[CHUNK_SIZE];
    for (int jj = 0; jj < NUM_CHUNKS * CHUNK_SIZE; jj += CHUNK_SIZE) {
        if (LOCAL_ID < CHUNK_SIZE) {
            chunkCharges[LOCAL_ID] = electrodeCharges[jj + LOCAL_ID];
            SYNC_WARPS;
            for (int k = 0; k < CHUNK_SIZE - 1; k++) {
                if (LOCAL_ID > k) {
                    chunkCharges[LOCAL_ID] -= chunkCharges[k] * capacitance[(mm_long) (jj + k) * NUM_ELECTRODE_PARTICLES + (LOCAL_ID + jj)];
                }
                SYNC_WARPS;
            }
            electrodeCharges[jj + LOCAL_ID] = chunkCharges[LOCAL_ID];
        }
        SYNC_THREADS;
        for (int ii = LOCAL_ID + jj + CHUNK_SIZE; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
            real chargeOffset = 0;
            for (int k = 0; k < CHUNK_SIZE; k++) {
                chargeOffset += chunkCharges[k] * capacitance[(mm_long) (jj + k) * NUM_ELECTRODE_PARTICLES + ii];
            }
            electrodeCharges[ii] -= chargeOffset;
        }
        SYNC_THREADS;
    }
    for (int jj = NUM_CHUNKS * CHUNK_SIZE; jj < NUM_ELECTRODE_PARTICLES; jj++) {
        const mm_long offset = (mm_long) jj * NUM_ELECTRODE_PARTICLES;
        for (int ii = LOCAL_ID + jj + 1; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
            electrodeCharges[ii] -= electrodeCharges[jj] * capacitance[offset + ii];
        }
        SYNC_THREADS;
    }

    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        const mm_long offset = (mm_long) ii * NUM_ELECTRODE_PARTICLES;
        electrodeCharges[ii] *= capacitance[offset + ii];
    }
    SYNC_THREADS;

    // Cholesky solve step 2 (outer loop over columns).

    for (int jj = NUM_ELECTRODE_PARTICLES - 1; jj >= 0; jj--) {
        const mm_long offset = (mm_long) jj * NUM_ELECTRODE_PARTICLES;

        for (int ii = LOCAL_ID; ii < jj; ii += LOCAL_SIZE) {
            // Retrieve capacitance[ii, jj] by retrieving capacitance[jj, ii] from the upper triangle.
            electrodeCharges[ii] -= electrodeCharges[jj] * capacitance[offset + ii];
        }
        SYNC_THREADS;
    }

    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        const mm_long offset = (mm_long) ii * NUM_ELECTRODE_PARTICLES;
        electrodeCharges[ii] *= capacitance[offset + ii];
    }
    SYNC_THREADS;

#ifdef USE_CHARGE_CONSTRAINT
    LOCAL volatile real temp[THREAD_BLOCK_SIZE];

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
