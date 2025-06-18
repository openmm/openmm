DEVICE real reduceLocalBuffer(LOCAL_ARG volatile real* temp) {
    const int thread = LOCAL_ID;
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

    LOCAL volatile real temp[THREAD_BLOCK_SIZE];
    const int thread = LOCAL_ID;

    for (int i = GLOBAL_ID; i < NUM_ELECTRODE_PARTICLES; i += GLOBAL_SIZE) {
        electrodeCharges[i] = -chargeDerivatives[i];
    }

    // Cholesky solve step 1 (outer loop over rows).
    mm_long offset = 0;
    for (int j = 0; j < NUM_ELECTRODE_PARTICLES; j++) {
        offset += j;
        temp[thread] = 0;
        for (int i = GLOBAL_ID; i < j; i += GLOBAL_SIZE) {
            // Retrieve capacitance[j, i].
            temp[thread] -= electrodeCharges[i] * capacitance[offset + i];
        }
        const real total = reduceLocalBuffer(temp);

        if (thread == 0) {
            // Retrieve 1 / capacitance[j, j] (reciprocal already taken on host).
            electrodeCharges[j] = (electrodeCharges[j] + total) * capacitance[offset + j];
        }
        SYNC_THREADS;
    }

    // Cholesky solve step 2 (outer loop over columns).
    offset = (mm_long)(NUM_ELECTRODE_PARTICLES - 1) * (NUM_ELECTRODE_PARTICLES + 2) / 2;
    for (int j = NUM_ELECTRODE_PARTICLES - 1; j >= 0; j--) {
        temp[thread] = 0;
        for (int i = GLOBAL_ID; i < NUM_ELECTRODE_PARTICLES; i += GLOBAL_SIZE) {
            if (i <= j) {
                continue;
            }
            // Retrieve capacitance[i, j].
            temp[thread] -= electrodeCharges[i] * capacitance[(mm_long)i * (i + 1) / 2 + j];
        }
        const real total = reduceLocalBuffer(temp);

        if (thread == 0) {
            // Retrieve 1 / capacitance[j, j] (reciprocal already taken on host).
            electrodeCharges[j] = (electrodeCharges[j] + total) * capacitance[offset];
        }
        SYNC_THREADS;
        offset -= j + 1;
    }

#ifdef USE_CHARGE_CONSTRAINT
    temp[thread] = thread == 0 ? chargeTarget : 0;
    for (int i = GLOBAL_ID; i < NUM_ELECTRODE_PARTICLES; i += GLOBAL_SIZE) {
        temp[thread] -= electrodeCharges[i];
    }
    const real chargeOffset = reduceLocalBuffer(temp);
    for (int i = GLOBAL_ID; i < NUM_ELECTRODE_PARTICLES; i += GLOBAL_SIZE) {
        electrodeCharges[i] += chargeOffset * constraintVector[i];
    }
#endif
}
