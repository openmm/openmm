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

    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        electrodeCharges[ii] = -chargeDerivatives[ii];
    }

    SYNC_THREADS;

    // Cholesky solve step 1 (outer loop over rows).
    for (int jj = 0; jj < NUM_ELECTRODE_PARTICLES; jj++) {
        const mm_long offset = (mm_long) jj * NUM_ELECTRODE_PARTICLES;

        if (LOCAL_ID == 0) {
            // Retrieve 1 / capacitance[jj, jj] (reciprocal already taken on host).
            electrodeCharges[jj] *= capacitance[offset + jj];
        }
        SYNC_THREADS;

        for (int ii = LOCAL_ID + jj + 1; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
            // Retrieve capacitance[jj, ii] from the lower triangle.
            electrodeCharges[ii] -= electrodeCharges[jj] * capacitance[offset + ii];
        }
        SYNC_THREADS;
    }

    // Cholesky solve step 2 (outer loop over columns).
    for (int jj = NUM_ELECTRODE_PARTICLES - 1; jj >= 0; jj--) {
        const mm_long offset = (mm_long) jj * NUM_ELECTRODE_PARTICLES;

        if (LOCAL_ID == 0) {
            // Retrieve 1 / capacitance[jj, jj] (reciprocal already taken on host).
            electrodeCharges[jj] *= capacitance[offset + jj];
        }
        SYNC_THREADS;

        for (int ii = LOCAL_ID; ii < jj; ii += LOCAL_SIZE) {
            // Retrieve capacitance[ii, jj] by retrieving capacitance[jj, ii] from the upper triangle.
            electrodeCharges[ii] -= electrodeCharges[jj] * capacitance[offset + ii];
        }
        SYNC_THREADS;
    }

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
