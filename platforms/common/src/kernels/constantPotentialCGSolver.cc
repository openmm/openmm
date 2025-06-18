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

KERNEL void solveInitializeStep1(GLOBAL real* RESTRICT electrodeCharges, GLOBAL real* RESTRICT qLast
#ifdef USE_CHARGE_CONSTRAINT
    , real chargeTarget
#endif
) {
    // This kernel expects to be executed in a single thread block.

#ifdef USE_CHARGE_CONSTRAINT
    LOCAL volatile real temp[THREAD_BLOCK_SIZE];
    real offset;
    const int thread = LOCAL_ID;
#endif

    // Set initial guess charges as linear extrapolations from the current and
    // previous charges fed through the solver, and save the current charges as
    // the previous charges.
    for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
        const real qGuess = electrodeCharges[ii];
        electrodeCharges[ii] = 2 * qGuess - qLast[ii];
        qLast[ii] = qGuess;
    }

#ifdef USE_CHARGE_CONSTRAINT
    // Ensure that initial guess charges satisfy the constraint.
    temp[thread] = 0;
    for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
        temp[thread] -= electrodeCharges[ii];
    }
    offset = (chargeTarget + reduceLocalBuffer(temp)) / NUM_ELECTRODE_PARTICLES;
    for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
        electrodeCharges[ii] += offset;
    }
#endif
}

KERNEL void solveInitializeStep2(GLOBAL real* RESTRICT chargeDerivatives, GLOBAL real* RESTRICT grad, GLOBAL real* RESTRICT projGrad, GLOBAL real* RESTRICT error) {
    // This kernel expects to be executed in a single thread block.

    LOCAL volatile real temp[THREAD_BLOCK_SIZE];
    real offset;
    const int thread = LOCAL_ID;

    for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
        grad[ii] = chargeDerivatives[ii];
    }

#ifdef USE_CHARGE_CONSTRAINT
    // Project the initial gradient without preconditioning.
    temp[thread] = 0;
    for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
        temp[thread] += grad[ii];
    }
    offset = reduceLocalBuffer(temp) / NUM_ELECTRODE_PARTICLES;
    for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
        projGrad[ii] = grad[ii] - offset;
    }
#else
    for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
        projGrad[ii] = grad[ii];
    }
#endif

    // Check for convergence at the initial guess charges.
    temp[thread] = 0;
    for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
        temp[thread] += projGrad[ii] * projGrad[ii];
    }
    offset = reduceLocalBuffer(temp);
    if (thread == 0) {
        error[0] = offset;
    }
}

KERNEL void solveInitializeStep3(GLOBAL real* RESTRICT chargeDerivatives, GLOBAL real* RESTRICT grad, GLOBAL real* RESTRICT projGrad, GLOBAL real* RESTRICT precGrad, GLOBAL real* RESTRICT qStep, GLOBAL real* RESTRICT grad0
#ifdef PRECOND_REQUESTED
    , GLOBAL real* RESTRICT precondVector, int precondActivated
#endif
) {
    // This kernel expects to be executed in a single thread block.

#ifdef USE_CHARGE_CONSTRAINT
    LOCAL volatile real temp[THREAD_BLOCK_SIZE];
    real offset;
    const int thread = LOCAL_ID;
#endif

    for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
        grad0[ii] = chargeDerivatives[ii];
    }

#ifdef PRECOND_REQUESTED
    // Project the initial gradient with preconditioning.
    if (precondActivated) {
#ifdef USE_CHARGE_CONSTRAINT
        temp[thread] = 0;
        for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
            temp[thread] += precondVector[ii] * grad[ii];
        }
        offset = reduceLocalBuffer(temp);
        for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
            precGrad[ii] = precondVector[ii] * (grad[ii] - offset);
        }
#else
        for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
            precGrad[ii] = precondVector[ii] * grad[ii];
        }
#endif
    }
    else {
#endif
        for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
            precGrad[ii] = projGrad[ii];
        }
#ifdef PRECOND_REQUESTED
    }
#endif

    // Initialize step vector for conjugate gradient iterations.
    for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
        qStep[ii] = -precGrad[ii];
    }
}

KERNEL void solveLoopStep1(GLOBAL real* RESTRICT chargeDerivatives, GLOBAL real* RESTRICT gradStep, GLOBAL real* RESTRICT grad0, GLOBAL real* RESTRICT error) {
    // This kernel expects to be executed in a single thread block.

    LOCAL volatile real temp[THREAD_BLOCK_SIZE];
    real offset;
    const int thread = LOCAL_ID;

    for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
        gradStep[ii] = chargeDerivatives[ii] - grad0[ii];
    }

    // If A qStep is small enough, stop to prevent, e.g., division by
    // zero in the calculation of alpha, or too large step sizes.
    temp[thread] = 0;
    for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
        temp[thread] += gradStep[ii] * gradStep[ii];
    }
    offset = reduceLocalBuffer(temp);
    if (thread == 0) {
        error[0] = offset;
    }
}

KERNEL void solveLoopStep2(GLOBAL real* RESTRICT electrodeCharges, GLOBAL real* RESTRICT q, GLOBAL real* RESTRICT grad, GLOBAL real* RESTRICT qStep, GLOBAL real* RESTRICT gradStep, GLOBAL real* RESTRICT paramScaleResult, int recomputeGradient
#ifdef USE_CHARGE_CONSTRAINT
    , real chargeTarget
#endif
) {
    // This kernel expects to be executed in a single thread block.

    LOCAL volatile real temp[THREAD_BLOCK_SIZE];
    real offset;
    const int thread = LOCAL_ID;

    // Evaluate the scalar 1 / (qStep^T A qStep).
    temp[thread] = 0;
    for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
        temp[thread] += qStep[ii] * gradStep[ii];
    }
    real paramScale = 1 / reduceLocalBuffer(temp);
    if (thread == 0) {
        paramScaleResult[0] = paramScale;
    }

    // Evaluate the conjugate gradient parameter alpha.
    temp[thread] = 0;
    for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
        temp[thread] -= qStep[ii] * grad[ii];
    }
    real alpha = reduceLocalBuffer(temp) * paramScale;

    // Update the charge vector.
    for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
        q[ii] += alpha * qStep[ii];
    }

#ifdef USE_CHARGE_CONSTRAINT
    // Remove any accumulated drift from the charge vector.  This
    // would be zero in exact arithmetic, but error can accumulate
    // over time in finite precision.
    temp[thread] = 0;
    for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
        temp[thread] -= q[ii];
    }
    offset = (chargeTarget + reduceLocalBuffer(temp)) / NUM_ELECTRODE_PARTICLES;
    for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
        q[ii] += offset;
    }
#endif

    if (recomputeGradient) {
        for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
            electrodeCharges[ii] = q[ii];
        }
    }
    else {
        for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
            grad[ii] += alpha * gradStep[ii];
        }
    }
}

KERNEL void solveLoopStep3(GLOBAL real* RESTRICT chargeDerivatives, GLOBAL real* RESTRICT grad, GLOBAL real* RESTRICT projGrad, int recomputeGradient, GLOBAL real* RESTRICT error) {
    // This kernel expects to be executed in a single thread block.

    LOCAL volatile real temp[THREAD_BLOCK_SIZE];
    real offset;
    const int thread = LOCAL_ID;

    if (recomputeGradient) {
        for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
            grad[ii] = chargeDerivatives[ii];
        }
    }

    // Project the current gradient without preconditioning.
#ifdef USE_CHARGE_CONSTRAINT
    temp[thread] = 0;
    for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
        temp[thread] += grad[ii];
    }
    offset = reduceLocalBuffer(temp) / NUM_ELECTRODE_PARTICLES;
    for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
        projGrad[ii] = grad[ii] - offset;
    }
#else
    for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
        projGrad[ii] = grad[ii];
    }
#endif
    
    // Check for convergence.
    temp[thread] = 0;
    for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
        temp[thread] += projGrad[ii] * projGrad[ii];
    }
    offset = reduceLocalBuffer(temp);
    if (thread == 0) {
        error[0] = offset;
    }
}

KERNEL void solveLoopStep4(GLOBAL real* RESTRICT grad, GLOBAL real* RESTRICT projGrad, GLOBAL real* RESTRICT precGrad, GLOBAL real* RESTRICT qStep, GLOBAL real* RESTRICT gradStep, GLOBAL real* RESTRICT paramScale
#ifdef PRECOND_REQUESTED
    , GLOBAL real* RESTRICT precondVector, int precondActivated
#endif
) {
    // This kernel expects to be executed in a single thread block.

    LOCAL volatile real temp[THREAD_BLOCK_SIZE];
    real offset;
    const int thread = LOCAL_ID;

    // Project the current gradient with preconditioning.
#ifdef PRECOND_REQUESTED
    if (precondActivated) {
#ifdef USE_CHARGE_CONSTRAINT
        temp[thread] = 0;
        for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
            temp[thread] += precondVector[ii] * grad[ii];
        }
        offset = reduceLocalBuffer(temp);
        for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
            precGrad[ii] = precondVector[ii] * (grad[ii] - offset);
        }
#else
        for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
            precGrad[ii] = precondVector[ii] * grad[ii];
        }
#endif
    }
    else {
#endif
        for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
            precGrad[ii] = projGrad[ii];
        }
#ifdef PRECOND_REQUESTED
    }
#endif

    // Evaluate the conjugate gradient parameter beta.
    temp[thread] = 0;
    for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
        temp[thread] += precGrad[ii] * gradStep[ii];
    }
    real beta = reduceLocalBuffer(temp) * paramScale[0];

    // Update the step vector.
    for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
        qStep[ii] = beta * qStep[ii] - precGrad[ii];
    }

#ifdef USE_CHARGE_CONSTRAINT
    // Project out any deviation off of the constraint plane from
    // the step vector.  This would be zero in exact arithmetic, but
    // error can accumulate over time in finite precision.

    temp[thread] = 0;
    for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
        temp[thread] += qStep[ii];
    }
    offset = reduceLocalBuffer(temp) / NUM_ELECTRODE_PARTICLES;
    for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
        qStep[ii] -= offset;
    }
#endif
}
