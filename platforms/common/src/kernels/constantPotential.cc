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

KERNEL void evaluateDirectDerivatives(
    GLOBAL real4* RESTRICT posq,
    GLOBAL real* RESTRICT charges,
    GLOBAL int* RESTRICT elecToSys,
    GLOBAL int* RESTRICT sysElec,
    GLOBAL real4* RESTRICT electrodeParams,
    GLOBAL mm_ulong* RESTRICT chargeDerivativesFixed,
    real4 periodicBoxSize,
    real4 invPeriodicBoxSize,
    real4 periodicBoxVecX,
    real4 periodicBoxVecY,
    real4 periodicBoxVecZ
) {
    for (int index = GLOBAL_ID; index < NUM_PARTICLES * NUM_ELECTRODE_PARTICLES; index += GLOBAL_SIZE) {
        const int ii = index / NUM_PARTICLES;
        const int j = index - ii * NUM_PARTICLES;
        const int i = elecToSys[ii];
        if (i == j) {
            continue;
        }

        const real4 posq1 = posq[i];
        const real4 posq2 = posq[j];
        real4 delta = posq2 - posq1;
        APPLY_PERIODIC_TO_DELTA(delta)
        const real r2 = delta.x * delta.x + delta.y * delta.y + delta.z * delta.z;
        if (r2 >= CUTOFF_SQUARED) {
            continue;
        }
        const real r = SQRT(r2);

#ifdef USE_POSQ_CHARGES
        const real charge2 = posq2.w;
#else
        const real charge2 = charges[j];
#endif
        const real4 params1 = electrodeParams[sysElec[i] + 1];
        const real4 params2 = electrodeParams[sysElec[j] + 1];

        const real alphaR = EWALD_ALPHA * r;
        const real etaR = r / SQRT(params1.y * params1.y + params2.y * params2.y);
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
        const real derivative = ONE_4PI_EPS0 * charge2 * (erfcAlphaR - erfcEtaR) / r;
        ATOMIC_ADD(&chargeDerivativesFixed[ii], (mm_ulong) realToFixedPoint(derivative));
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
