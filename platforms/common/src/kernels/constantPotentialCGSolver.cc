#define REDUCE_BODY_1 \
    const int thread = LOCAL_ID; \
    SYNC_THREADS; \
    temp[thread] = value; \
    SYNC_THREADS; \
    for (int step = 1; step < 32; step *= 2) { \
        if (thread + step < LOCAL_SIZE && thread % (2 * step) == 0) {

#define REDUCE_BODY_2 \
        } \
        SYNC_WARPS; \
    } \
    for (int step = 32; step < LOCAL_SIZE; step *= 2) { \
        if (thread + step < LOCAL_SIZE && thread % (2 * step) == 0) {

#define REDUCE_BODY_3 \
        } \
        SYNC_THREADS; \
    }

// Sum value from each thread (using temp).  Use real type variables (float on
// single and mixed precision modes, double on double precision mode).
DEVICE real reduceReal(real value, LOCAL_ARG volatile real* temp) {
    REDUCE_BODY_1
    temp[thread] = temp[thread] + temp[thread + step];
    REDUCE_BODY_2
    temp[thread] = temp[thread] + temp[thread + step];
    REDUCE_BODY_3
    return temp[0];
}

// We need more than single precision for accumulation regardless of the mode
// selected, so use double if double precision is supported, and otherwise use
// double-float arithmetic.  In the latter case, use float2 with .x storing the
// "high" part and .y storing the "low" part of each double-float number.
#ifdef SUPPORTS_DOUBLE_PRECISION

#define ACCUM double
#define ACCUM_ZERO 0.0

// Perform accum = accum + real.
#define ACCUM_ADD(x, y) ((x) + (y))
// Perform real = real + accum.
#define ACCUM_APPLY(x, y) ((real) ((x) + (y)))

// Sum value from each thread (using temp) and return (sum + offset) * scale.
DEVICE ACCUM reduceAccum(ACCUM value, LOCAL_ARG volatile ACCUM* temp, real offset, real scale) {
    REDUCE_BODY_1
    temp[thread] = temp[thread] + temp[thread + step];
    REDUCE_BODY_2
    temp[thread] = temp[thread] + temp[thread + step];
    REDUCE_BODY_3
    return (temp[0] + offset) * scale;
}

#else

#define ACCUM float2
#define ACCUM_ZERO make_float2(0.0f, 0.0f)

#define ACCUM_ADD(x, y) compensatedAdd2(x, y)
#define ACCUM_APPLY(x, y) compensatedAdd1(y, x)

// For details of the compensated summation implemented, see Joldes et al.,
// ACM Trans. Math. Softw. 2017, 44, 15res (DOI: 10.1145/3121432).

// float + float -> float2, only valid if the floating-point exponent of a is
// not less than that of b.
DEVICE inline float2 compensatedAddKernel1(float a, float b) {
    float s = a + b;
    return make_float2(s, b - (s - a));
}

// float + float -> float2, valid for any inputs.
DEVICE inline float2 compensatedAddKernel2(float a, float b) {
    float s = a + b;
    float c = s - b;
    float d = s - c;
    return make_float2(s, (a - c) + (b - d));
}

// float * float -> float2.
DEVICE inline float2 compensatedMultiplyKernel(float a, float b) {
    float c = a * b;
    return make_float2(c, FMA(a, b, -c));
}

// float2 + float -> float.  Like compensatedAdd2, but only computes the high
// part of the result.
DEVICE inline float compensatedAdd1(float2 x, float y) {
    float2 s = compensatedAddKernel2(x.x, y);
    return s.x + (x.y + s.y);
}

// float2 + float -> float2, with a relative error of 2^-47.
DEVICE inline float2 compensatedAdd2(float2 x, float y) {
    float2 s = compensatedAddKernel2(x.x, y);
    float v = x.y + s.y;
    return compensatedAddKernel1(s.x, v);
}

// float2 + float2 -> float2, with a relative error of 2^-46.
DEVICE inline float2 compensatedAdd3(float2 x, float2 y) {
    float2 s = compensatedAddKernel2(x.x, y.x);
    float2 t = compensatedAddKernel2(x.y, y.y);
    float c = s.y + t.x;
    float2 v = compensatedAddKernel1(s.x, c);
    float w = t.y + v.y;
    return compensatedAddKernel1(v.x, w);
}

// float2 * float -> float2, with a relative error of 2^-47.
DEVICE inline float2 compensatedMultiply(float2 x, float y) {
    float c = x.x * y;
    return compensatedAddKernel1(c, FMA(x.x, y, -c) + x.y * y);
}

// Sum value from each thread (using temp) and return (sum + offset) * scale.
DEVICE ACCUM reduceAccum(ACCUM value, LOCAL_ARG volatile ACCUM* temp, real offset, real scale) {
    REDUCE_BODY_1
    temp[thread] = compensatedAdd3(temp[thread], temp[thread + step]);
    REDUCE_BODY_2
    temp[thread] = compensatedAdd3(temp[thread], temp[thread + step]);
    REDUCE_BODY_3
    return compensatedMultiply(compensatedAdd2(temp[0], offset), scale);
}

#endif

KERNEL void solveInitializeStep1(GLOBAL real* RESTRICT electrodeCharges, GLOBAL real* RESTRICT qLast
#ifdef USE_CHARGE_CONSTRAINT
    , real chargeTarget
#endif
) {
    // This kernel expects to be executed in a single thread block.

#ifdef USE_CHARGE_CONSTRAINT
    LOCAL volatile ACCUM tempAccum[THREAD_BLOCK_SIZE];
#endif

    // Set initial guess charges as linear extrapolations from the current and
    // previous charges fed through the solver, and save the current charges as
    // the previous charges.
    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        const real qGuess = electrodeCharges[ii];
        electrodeCharges[ii] = 2 * qGuess - qLast[ii];
        qLast[ii] = qGuess;
    }

#ifdef USE_CHARGE_CONSTRAINT
    // Ensure that initial guess charges satisfy the constraint.
    ACCUM offsetAccum = ACCUM_ZERO;
    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        offsetAccum = ACCUM_ADD(offsetAccum, -electrodeCharges[ii]);
    }
    const ACCUM offset = reduceAccum(offsetAccum, tempAccum, chargeTarget, 1 / (real) NUM_ELECTRODE_PARTICLES);
    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        electrodeCharges[ii] = ACCUM_APPLY(electrodeCharges[ii], offset);
    }
#endif
}

KERNEL void solveInitializeStep2(GLOBAL real* RESTRICT chargeDerivatives, GLOBAL real* RESTRICT grad, GLOBAL real* RESTRICT projGrad, GLOBAL real* RESTRICT errorResult) {
    // This kernel expects to be executed in a single thread block.

    LOCAL volatile ACCUM tempAccum[THREAD_BLOCK_SIZE];
    LOCAL volatile real temp[THREAD_BLOCK_SIZE];

    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        grad[ii] = chargeDerivatives[ii];
    }

#ifdef USE_CHARGE_CONSTRAINT
    // Project the initial gradient without preconditioning.
    ACCUM offsetAccum = ACCUM_ZERO;
    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        offsetAccum = ACCUM_ADD(offsetAccum, grad[ii]);
    }
    const ACCUM offset = reduceAccum(offsetAccum, tempAccum, 0, -1 / (real) NUM_ELECTRODE_PARTICLES);
    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        projGrad[ii] = ACCUM_APPLY(grad[ii], offset);
    }
#else
    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        projGrad[ii] = grad[ii];
    }
#endif

    // Check for convergence at the initial guess charges.
    real error = 0;
    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        error += projGrad[ii] * projGrad[ii];
    }
    error = reduceReal(error, temp);
    if (LOCAL_ID == 0) {
        errorResult[0] = error;
    }
}

KERNEL void solveInitializeStep3(GLOBAL real* RESTRICT chargeDerivatives, GLOBAL real* RESTRICT grad, GLOBAL real* RESTRICT projGrad, GLOBAL real* RESTRICT precGrad, GLOBAL real* RESTRICT qStep, GLOBAL real* RESTRICT grad0
#ifdef PRECOND_REQUESTED
    , GLOBAL real* RESTRICT precondVector, int precondActivated
#endif
) {
    // This kernel expects to be executed in a single thread block.

#if defined(PRECOND_REQUESTED) && defined(USE_CHARGE_CONSTRAINT)
    LOCAL volatile ACCUM tempAccum[THREAD_BLOCK_SIZE];
#endif

    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        grad0[ii] = chargeDerivatives[ii];
    }

#ifdef PRECOND_REQUESTED
    // Project the initial gradient with preconditioning.
    if (precondActivated) {
#ifdef USE_CHARGE_CONSTRAINT
        ACCUM offsetAccum = ACCUM_ZERO;
        for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
            offsetAccum = ACCUM_ADD(offsetAccum, precondVector[ii] * grad[ii]);
        }
        const ACCUM offset = reduceAccum(offsetAccum, tempAccum, 0, -1);
        for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
            precGrad[ii] = precondVector[ii] * ACCUM_APPLY(grad[ii], offset);
        }
#else
        for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
            precGrad[ii] = precondVector[ii] * grad[ii];
        }
#endif
    }
    else {
#endif
        for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
            precGrad[ii] = projGrad[ii];
        }
#ifdef PRECOND_REQUESTED
    }
#endif

    // Initialize step vector for conjugate gradient iterations.
    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        qStep[ii] = -precGrad[ii];
    }
}

KERNEL void solveLoopStep1(GLOBAL real* RESTRICT chargeDerivatives, GLOBAL real* RESTRICT gradStep, GLOBAL real* RESTRICT grad0, GLOBAL real* RESTRICT errorResult) {
    // This kernel expects to be executed in a single thread block.

    LOCAL volatile real temp[THREAD_BLOCK_SIZE];

    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        gradStep[ii] = chargeDerivatives[ii] - grad0[ii];
    }

    // If A qStep is small enough, stop to prevent, e.g., division by zero in
    // the calculation of alpha, or too large step sizes.
    real error = 0;
    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        error += gradStep[ii] * gradStep[ii];
    }
    error = reduceReal(error, temp);
    if (LOCAL_ID == 0) {
        errorResult[0] = error;
    }
}

KERNEL void solveLoopStep2(GLOBAL real* RESTRICT electrodeCharges, GLOBAL real* RESTRICT q, GLOBAL real* RESTRICT grad, GLOBAL real* RESTRICT qStep, GLOBAL real* RESTRICT gradStep, GLOBAL real* RESTRICT paramScaleResult, int recomputeGradient
#ifdef USE_CHARGE_CONSTRAINT
    , real chargeTarget
#endif
) {
    // This kernel expects to be executed in a single thread block.

    LOCAL volatile ACCUM tempAccum[THREAD_BLOCK_SIZE];
    LOCAL volatile real temp[THREAD_BLOCK_SIZE];

    // Evaluate the scalar 1 / (qStep^T A qStep).
    real paramScale = 0;
    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        paramScale += qStep[ii] * gradStep[ii];
    }
    paramScale = 1 / reduceReal(paramScale, temp);
    if (LOCAL_ID == 0) {
        paramScaleResult[0] = paramScale;
    }

    // Evaluate the conjugate gradient parameter alpha.
    real alpha = 0;
    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        alpha -= qStep[ii] * grad[ii];
    }
    alpha = reduceReal(alpha, temp) * paramScale;

    // Update the charge vector.
    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        q[ii] += alpha * qStep[ii];
    }

#ifdef USE_CHARGE_CONSTRAINT
    // Remove any accumulated drift from the charge vector.  This would be zero
    // in exact arithmetic, but error can accumulate over time in finite
    // precision.
    ACCUM offsetAccum = ACCUM_ZERO;
    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        offsetAccum = ACCUM_ADD(offsetAccum, -q[ii]);
    }
    const ACCUM offset = reduceAccum(offsetAccum, tempAccum, chargeTarget, 1 / (real) NUM_ELECTRODE_PARTICLES);
    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        q[ii] = ACCUM_APPLY(q[ii], offset);
    }
#endif

    if (recomputeGradient) {
        for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
            electrodeCharges[ii] = q[ii];
        }
    }
    else {
        for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
            grad[ii] += alpha * gradStep[ii];
        }
    }
}

KERNEL void solveLoopStep3(GLOBAL real* RESTRICT chargeDerivatives, GLOBAL real* RESTRICT grad, GLOBAL real* RESTRICT projGrad, int recomputeGradient, GLOBAL real* RESTRICT errorResult) {
    // This kernel expects to be executed in a single thread block.

    LOCAL volatile ACCUM tempAccum[THREAD_BLOCK_SIZE];
    LOCAL volatile real temp[THREAD_BLOCK_SIZE];

    if (recomputeGradient) {
        for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
            grad[ii] = chargeDerivatives[ii];
        }
    }

    // Project the current gradient without preconditioning.
#ifdef USE_CHARGE_CONSTRAINT
    ACCUM offsetAccum = ACCUM_ZERO;
    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        offsetAccum = ACCUM_ADD(offsetAccum, grad[ii]);
    }
    const ACCUM offset = reduceAccum(offsetAccum, tempAccum, 0, -1 / (real) NUM_ELECTRODE_PARTICLES);
    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        projGrad[ii] = ACCUM_APPLY(grad[ii], offset);
    }
#else
    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        projGrad[ii] = grad[ii];
    }
#endif
    
    // Check for convergence.
    real error = 0;
    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        error += projGrad[ii] * projGrad[ii];
    }
    error = reduceReal(error, temp);
    if (LOCAL_ID == 0) {
        errorResult[0] = error;
    }
}

KERNEL void solveLoopStep4(GLOBAL real* RESTRICT grad, GLOBAL real* RESTRICT projGrad, GLOBAL real* RESTRICT precGrad, GLOBAL real* RESTRICT qStep, GLOBAL real* RESTRICT gradStep, GLOBAL real* RESTRICT paramScale
#ifdef PRECOND_REQUESTED
    , GLOBAL real* RESTRICT precondVector, int precondActivated
#endif
) {
    // This kernel expects to be executed in a single thread block.

    LOCAL volatile ACCUM tempAccum[THREAD_BLOCK_SIZE];
    LOCAL volatile real temp[THREAD_BLOCK_SIZE];

    // Project the current gradient with preconditioning.
#ifdef PRECOND_REQUESTED
    if (precondActivated) {
#ifdef USE_CHARGE_CONSTRAINT
        ACCUM offsetAccum = ACCUM_ZERO;
        for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
            offsetAccum = ACCUM_ADD(offsetAccum, precondVector[ii] * grad[ii]);
        }
        const ACCUM offset = reduceAccum(offsetAccum, tempAccum, 0, -1);
        for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
            precGrad[ii] = precondVector[ii] * ACCUM_APPLY(grad[ii], offset);
        }
#else
        for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
            precGrad[ii] = precondVector[ii] * grad[ii];
        }
#endif
    }
    else {
#endif
        for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
            precGrad[ii] = projGrad[ii];
        }
#ifdef PRECOND_REQUESTED
    }
#endif

    // Evaluate the conjugate gradient parameter beta.
    real beta = 0;
    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        beta += precGrad[ii] * gradStep[ii];
    }
    beta = reduceReal(beta, temp) * paramScale[0];

    // Update the step vector.
    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        qStep[ii] = beta * qStep[ii] - precGrad[ii];
    }

#ifdef USE_CHARGE_CONSTRAINT
    // Project out any deviation off of the constraint plane from the step
    // vector.  This would be zero in exact arithmetic, but error can accumulate
    // over time in finite precision.

    ACCUM offsetAccum = ACCUM_ZERO;
    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        offsetAccum = ACCUM_ADD(offsetAccum, qStep[ii]);
    }
    const ACCUM offset = reduceAccum(offsetAccum, tempAccum, 0, -1 / (real) NUM_ELECTRODE_PARTICLES);
    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        qStep[ii] = ACCUM_APPLY(qStep[ii], offset);
    }
#endif
}
