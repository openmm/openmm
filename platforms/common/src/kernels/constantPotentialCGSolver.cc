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

typedef struct {
    real gradStepSq, qStepGradStep, qStepGrad, q, qStep;
} BlockSums1;

typedef struct {
    real projGradSq, precGradStep, precGrad;
} BlockSums2;

// Sum value from each thread (using temp).  Use real type variables (float on
// single and mixed precision modes, double on double precision mode).
DEVICE real reduceReal(real value, LOCAL_ARG volatile real* temp) {
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

// Performs the equivalent of reduceReal() on 5 values simultaneously.
DEVICE BlockSums1 reduceBlockSums1(BlockSums1 value, LOCAL_ARG BlockSums1* temp) {
    const int thread = LOCAL_ID;
    SYNC_THREADS;
#ifdef WARP_SHUFFLE_DOWN
    const int warpCount = LOCAL_SIZE / WARP_SIZE;
    const int warp = thread / WARP_SIZE;
    const int lane = thread % WARP_SIZE;
    for (int step = WARP_SIZE / 2; step > 0; step >>= 1) {
        value.gradStepSq += WARP_SHUFFLE_DOWN(value.gradStepSq, step);
        value.qStepGradStep += WARP_SHUFFLE_DOWN(value.qStepGradStep, step);
        value.qStepGrad += WARP_SHUFFLE_DOWN(value.qStepGrad, step);
        value.q += WARP_SHUFFLE_DOWN(value.q, step);
        value.qStep += WARP_SHUFFLE_DOWN(value.qStep, step);
    }
    if (!lane) {
        temp[warp] = value;
    }
    SYNC_THREADS;
    if (!warp) {
        value.gradStepSq = value.qStepGradStep = value.qStepGrad = value.q = value.qStep = 0;
        if (lane < warpCount) {
            value = temp[lane];
        }
        for (int step = WARP_SIZE / 2; step > 0; step >>= 1) {
            value.gradStepSq += WARP_SHUFFLE_DOWN(value.gradStepSq, step);
            value.qStepGradStep += WARP_SHUFFLE_DOWN(value.qStepGradStep, step);
            value.qStepGrad += WARP_SHUFFLE_DOWN(value.qStepGrad, step);
            value.q += WARP_SHUFFLE_DOWN(value.q, step);
            value.qStep += WARP_SHUFFLE_DOWN(value.qStep, step);
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
            temp[thread].gradStepSq += temp[thread + step].gradStepSq;
            temp[thread].qStepGradStep += temp[thread + step].qStepGradStep;
            temp[thread].qStepGrad += temp[thread + step].qStepGrad;
            temp[thread].q += temp[thread + step].q;
            temp[thread].qStep += temp[thread + step].qStep;
        }
        SYNC_WARPS;
    }
    for (int step = WARP_SIZE / 2; step < LOCAL_SIZE; step <<= 1) {
        if(thread + step < LOCAL_SIZE && thread % (2 * step) == 0) {
            temp[thread].gradStepSq += temp[thread + step].gradStepSq;
            temp[thread].qStepGradStep += temp[thread + step].qStepGradStep;
            temp[thread].qStepGrad += temp[thread + step].qStepGrad;
            temp[thread].q += temp[thread + step].q;
            temp[thread].qStep += temp[thread + step].qStep;
        }
        SYNC_THREADS;
    }
#endif
    return temp[0];
}

// Performs the equivalent of reduceReal() on 3 values simultaneously.
DEVICE BlockSums2 reduceBlockSums2(BlockSums2 value, LOCAL_ARG BlockSums2* temp) {
    const int thread = LOCAL_ID;
    SYNC_THREADS;
#ifdef WARP_SHUFFLE_DOWN
    const int warpCount = LOCAL_SIZE / WARP_SIZE;
    const int warp = thread / WARP_SIZE;
    const int lane = thread % WARP_SIZE;
    for (int step = WARP_SIZE / 2; step > 0; step >>= 1) {
        value.projGradSq += WARP_SHUFFLE_DOWN(value.projGradSq, step);
        value.precGradStep += WARP_SHUFFLE_DOWN(value.precGradStep, step);
        value.precGrad += WARP_SHUFFLE_DOWN(value.precGrad, step);
    }
    if (!lane) {
        temp[warp] = value;
    }
    SYNC_THREADS;
    if (!warp) {
        value.projGradSq = value.precGradStep = value.precGrad = 0;
        if (lane < warpCount) {
            value = temp[lane];
        }
        for (int step = WARP_SIZE / 2; step > 0; step >>= 1) {
            value.projGradSq += WARP_SHUFFLE_DOWN(value.projGradSq, step);
            value.precGradStep += WARP_SHUFFLE_DOWN(value.precGradStep, step);
            value.precGrad += WARP_SHUFFLE_DOWN(value.precGrad, step);
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
            temp[thread].projGradSq += temp[thread + step].projGradSq;
            temp[thread].precGradStep += temp[thread + step].precGradStep;
            temp[thread].precGrad += temp[thread + step].precGrad;
        }
        SYNC_WARPS;
    }
    for (int step = WARP_SIZE / 2; step < LOCAL_SIZE; step <<= 1) {
        if(thread + step < LOCAL_SIZE && thread % (2 * step) == 0) {
            temp[thread].projGradSq += temp[thread + step].projGradSq;
            temp[thread].precGradStep += temp[thread + step].precGradStep;
            temp[thread].precGrad += temp[thread + step].precGrad;
        }
        SYNC_THREADS;
    }
#endif
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
#define ACCUM_ADD(x, y) ((x) + (ACCUM) (y))
// Perform real = real + accum.
#define ACCUM_APPLY(x, y) ((real) ((ACCUM) (x) + (y)))
// Perform real = accum * real.
#define ACCUM_MUL(x, y) ((real) ((x) * (ACCUM) (y)))
// Perform accum = (accum * real) + accum.
#define ACCUM_MUL_ADD(x, y, z) (((x) * (ACCUM) (y)) + (z))
// Perform real = (real + accum) * accum.
#define ACCUM_ADD_MUL(x, y, z) ((real) (((ACCUM) (x) + (y)) * (z)))

// Sum value from each thread (using temp) and return (sum + offset) * scale.
DEVICE ACCUM reduceAccum(ACCUM value, LOCAL_ARG volatile ACCUM* temp, real offset, real scale) {
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
    return (temp[0] + offset) * scale;
}

#else

#define ACCUM float2
#define ACCUM_ZERO make_float2(0.0f, 0.0f)

#define ACCUM_ADD(x, y) compensatedAdd2(x, y)
#define ACCUM_APPLY(x, y) compensatedAdd1(y, x)
#define ACCUM_MUL(x, y) compensatedMultiply1(x, y)
#define ACCUM_MUL_ADD(x, y, z) compensatedAdd3(compensatedMultiply2(x, y), z)
#define ACCUM_ADD_MUL(x, y, z) compensatedMultiply3(compensatedAdd2(y, x), z)

// For details of the compensated addition and multiplication implemented, see
// Joldes et al., ACM Trans. Math. Softw. 2017, 44, 15res (DOI: 10.1145/3121432).

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

// float2 * float -> float.  Like compensatedMultiply2, but only computes the
// high part of the result.
DEVICE inline float compensatedMultiply1(float2 x, float y) {
    float c = x.x * y;
    return c + (FMA(x.x, y, -c) + x.y * y);
}

// float2 * float -> float2, with a relative error of 2^-47.
DEVICE inline float2 compensatedMultiply2(float2 x, float y) {
    float c = x.x * y;
    return compensatedAddKernel1(c, FMA(x.x, y, -c) + x.y * y);
}

// float2 * float2 -> float.
DEVICE inline float compensatedMultiply3(float2 x, float2 y) {
    float2 c = compensatedMultiplyKernel(x.x, y.x);
    return c.x + (c.y + FMA(x.y, y.x, FMA(x.x, y.y, x.y * y.y)));
}

// Sum value from each thread (using temp) and return (sum + offset) * scale.
DEVICE ACCUM reduceAccum(ACCUM value, LOCAL_ARG volatile ACCUM* temp, real offset, real scale) {
    const int thread = LOCAL_ID;
    SYNC_THREADS;
#ifdef WARP_SHUFFLE_DOWN
    const int warpCount = LOCAL_SIZE / WARP_SIZE;
    const int warp = thread / WARP_SIZE;
    const int lane = thread % WARP_SIZE;
    for (int step = WARP_SIZE / 2; step > 0; step >>= 1) {
        value = compensatedAdd3(value, WARP_SHUFFLE_DOWN(value, step));
    }
    if (!lane) {
        temp[warp] = value;
    }
    SYNC_THREADS;
    if (!warp) {
        value = lane < warpCount ? temp[lane] : ACCUM_ZERO;
        for (int step = WARP_SIZE / 2; step > 0; step >>= 1) {
            value = compensatedAdd3(value, WARP_SHUFFLE_DOWN(value, step));
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
            temp[thread] = compensatedAdd3(temp[thread], temp[thread + step]);
        }
        SYNC_WARPS;
    }
    for (int step = WARP_SIZE / 2; step < LOCAL_SIZE; step <<= 1) {
        if(thread + step < LOCAL_SIZE && thread % (2 * step) == 0) {
            temp[thread] = compensatedAdd3(temp[thread], temp[thread + step]);
        }
        SYNC_THREADS;
    }
#endif
    return compensatedMultiply2(compensatedAdd2(temp[0], offset), scale);
}

#endif

KERNEL void solveInitializeStep1(GLOBAL real* RESTRICT electrodeCharges, GLOBAL real* RESTRICT qLast
#ifdef USE_CHARGE_CONSTRAINT
    , real chargeTarget
#endif
) {
    // This kernel expects to be executed in a single thread block.

#ifdef USE_CHARGE_CONSTRAINT
    LOCAL volatile ACCUM tempAccum[TEMP_SIZE];
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

KERNEL void solveInitializeStep2(GLOBAL real* RESTRICT chargeDerivatives, GLOBAL real* RESTRICT grad, GLOBAL real* RESTRICT projGrad, GLOBAL int* RESTRICT convergedResult) {
    // This kernel expects to be executed in a single thread block.

    LOCAL volatile real temp[TEMP_SIZE];
#ifdef USE_CHARGE_CONSTRAINT
    LOCAL volatile ACCUM tempAccum[TEMP_SIZE];
#endif

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
        convergedResult[0] = (int) (error <= ERROR_TARGET);
    }
}

KERNEL void solveInitializeStep3(GLOBAL real* RESTRICT electrodeCharges, GLOBAL real* RESTRICT chargeDerivatives, GLOBAL real* RESTRICT grad, GLOBAL real* RESTRICT projGrad,
    GLOBAL real* RESTRICT precGrad, GLOBAL real* RESTRICT qStep, GLOBAL real* RESTRICT grad0
#ifdef PRECOND_REQUESTED
    , GLOBAL ACCUM* RESTRICT precondVector, int precondActivated
#endif
) {
    // This kernel expects to be executed in a single thread block.

#if defined(PRECOND_REQUESTED) && defined(USE_CHARGE_CONSTRAINT)
    LOCAL volatile ACCUM tempAccum[TEMP_SIZE];
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
            offsetAccum = ACCUM_MUL_ADD(precondVector[ii], grad[ii], offsetAccum);
        }
        const ACCUM offset = reduceAccum(offsetAccum, tempAccum, 0, -1);
        for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
            precGrad[ii] = ACCUM_ADD_MUL(grad[ii], offset, precondVector[ii]);
        }
#else
        for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
            precGrad[ii] = ACCUM_MUL(precondVector[ii], grad[ii]);
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
        electrodeCharges[ii] = qStep[ii] = -precGrad[ii];
    }
}

KERNEL void solveLoopStep1(
    GLOBAL real* RESTRICT chargeDerivatives,
    GLOBAL real* RESTRICT q,
    GLOBAL real* RESTRICT grad,
    GLOBAL real* RESTRICT qStep,
    GLOBAL real* RESTRICT gradStep,
    GLOBAL real* RESTRICT grad0,
    GLOBAL BlockSums1* RESTRICT blockSums1Buffer
) {
    // This kernel can be executed across multiple thread blocks.

    LOCAL BlockSums1 temp[TEMP_SIZE];

    for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
        gradStep[ii] = chargeDerivatives[ii] - grad0[ii];
    }

    // Reduce values within each block and store results.
    BlockSums1 blockSums1 = {0, 0, 0, 0, 0};
    for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
        blockSums1.gradStepSq += gradStep[ii] * gradStep[ii];
        blockSums1.qStepGradStep += qStep[ii] * gradStep[ii];
        blockSums1.qStepGrad += qStep[ii] * grad[ii];
        blockSums1.q += q[ii];
        blockSums1.qStep += qStep[ii];
    }
    blockSums1 = reduceBlockSums1(blockSums1, temp);

    if (LOCAL_ID == 0) {
        blockSums1Buffer[GROUP_ID] = blockSums1;
    }
}

KERNEL void solveLoopStep2(
    GLOBAL BlockSums1* RESTRICT blockSums1Buffer,
    GLOBAL int* RESTRICT convergedResult
) {
    // This kernel expects to be executed in a single thread block.

    LOCAL BlockSums1 tempSums[TEMP_SIZE];

    // Reduce values from all blocks.
    BlockSums1 blockSums1 = {0, 0, 0, 0, 0};
    for (int ii = LOCAL_ID; ii < THREAD_BLOCK_COUNT; ii += LOCAL_SIZE) {
        blockSums1.gradStepSq += blockSums1Buffer[ii].gradStepSq;
        blockSums1.qStepGradStep += blockSums1Buffer[ii].qStepGradStep;
        blockSums1.qStepGrad += blockSums1Buffer[ii].qStepGrad;
        blockSums1.q += blockSums1Buffer[ii].q;
        blockSums1.qStep += blockSums1Buffer[ii].qStep;
    }
    blockSums1 = reduceBlockSums1(blockSums1, tempSums);

    if (LOCAL_ID == 0) {
        blockSums1Buffer[0] = blockSums1;
        // If A qStep is small enough, stop to prevent, e.g., division by zero
        // in the calculation of alpha, or too large step sizes.
        convergedResult[0] = (int) (blockSums1.gradStepSq <= ERROR_TARGET);
    }
}

KERNEL void solveLoopStep3(
    GLOBAL real* RESTRICT q,
    GLOBAL real* RESTRICT grad,
    GLOBAL real* RESTRICT qStep,
    GLOBAL real* RESTRICT gradStep,
    GLOBAL BlockSums1* RESTRICT blockSums1Buffer,
    GLOBAL int* RESTRICT convergedResult
#ifdef USE_CHARGE_CONSTRAINT
    , real chargeTarget
#endif
) {
    // This kernel can be executed across multiple thread blocks.

    if (convergedResult[0] != 0) {
        return;
    }

    const BlockSums1 blockSums1 = blockSums1Buffer[0];
    const real alpha = -blockSums1.qStepGrad / blockSums1.qStepGradStep;

    // Update the charge vector.
    for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
        q[ii] += alpha * qStep[ii];
    }

#ifdef USE_CHARGE_CONSTRAINT
    // Remove any accumulated drift from the charge vector.  This would be zero
    // in exact arithmetic, but error can accumulate over time in finite
    // precision.
    const real offset = (chargeTarget - (blockSums1.q + alpha * blockSums1.qStep)) / NUM_ELECTRODE_PARTICLES;
    for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
        q[ii] += offset;
    }
#endif

    // Update the gradient vector.  If on this iteration, the gradient is to be
    // recomputed, the contents of grad will be overwritten.
    for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
        grad[ii] += alpha * gradStep[ii];
    }
}

KERNEL void solveLoopStep4(
    GLOBAL real* RESTRICT grad,
    GLOBAL real* RESTRICT projGrad,
    GLOBAL real* RESTRICT precGrad,
    GLOBAL real* RESTRICT gradStep,
    GLOBAL BlockSums2* RESTRICT blockSums2Buffer,
    GLOBAL int* RESTRICT convergedResult
#ifdef PRECOND_REQUESTED
    , GLOBAL ACCUM* RESTRICT precondVector,
    int precondActivated
#endif
) {
    // This kernel expects to be executed in a single thread block.

    if (convergedResult[0] != 0) {
        return;
    }

    LOCAL volatile real temp[TEMP_SIZE];
    LOCAL BlockSums2 tempSums[TEMP_SIZE];
#ifdef USE_CHARGE_CONSTRAINT
    LOCAL volatile ACCUM tempAccum[TEMP_SIZE];
#endif

    // Project the current gradient without preconditioning.
#ifdef USE_CHARGE_CONSTRAINT
    ACCUM offsetAccum = ACCUM_ZERO;
    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        offsetAccum = ACCUM_ADD(offsetAccum, grad[ii]);
    }
    ACCUM offset = reduceAccum(offsetAccum, tempAccum, 0, -1 / (real) NUM_ELECTRODE_PARTICLES);
    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        projGrad[ii] = ACCUM_APPLY(grad[ii], offset);
    }
#else
    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        projGrad[ii] = grad[ii];
    }
#endif

    // Project the current gradient with preconditioning.
#ifdef PRECOND_REQUESTED
    if (precondActivated) {
#ifdef USE_CHARGE_CONSTRAINT
        offsetAccum = ACCUM_ZERO;
        for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
            offsetAccum = ACCUM_MUL_ADD(precondVector[ii], grad[ii], offsetAccum);
        }
        offset = reduceAccum(offsetAccum, tempAccum, 0, -1);
        for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
            precGrad[ii] = ACCUM_ADD_MUL(grad[ii], offset, precondVector[ii]);
        }
#else
        for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
            precGrad[ii] = ACCUM_MUL(precondVector[ii], grad[ii]);
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

    // Reduce values to be used by all blocks in the final kernel.
    BlockSums2 blockSums2 = {0, 0, 0};
    for (int ii = LOCAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += LOCAL_SIZE) {
        blockSums2.projGradSq += projGrad[ii] * projGrad[ii];
        blockSums2.precGradStep += precGrad[ii] * gradStep[ii];
        blockSums2.precGrad += precGrad[ii];
    }
    blockSums2 = reduceBlockSums2(blockSums2, tempSums);

    if (LOCAL_ID == 0) {
        blockSums2Buffer[0] = blockSums2;
        convergedResult[0] = (int) (blockSums2.projGradSq <= ERROR_TARGET);
    }
}

KERNEL void solveLoopStep5(
    GLOBAL real* RESTRICT electrodeCharges,
    GLOBAL real* RESTRICT precGrad,
    GLOBAL real* RESTRICT qStep,
    GLOBAL BlockSums1* RESTRICT blockSums1Buffer,
    GLOBAL BlockSums2* RESTRICT blockSums2Buffer,
    GLOBAL int* RESTRICT convergedResult
) {
    // This kernel can be executed across multiple thread blocks.

    if (convergedResult[0] != 0) {
        return;
    }

    const BlockSums1 blockSums1 = blockSums1Buffer[0];
    const BlockSums2 blockSums2 = blockSums2Buffer[0];

    // Evaluate the conjugate gradient parameter beta.
    const real beta = blockSums2.precGradStep / blockSums1.qStepGradStep;

    // Update the step vector.
    for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
        qStep[ii] = beta * qStep[ii] - precGrad[ii];
    }

#ifdef USE_CHARGE_CONSTRAINT
    // Project out any deviation off of the constraint plane from the step
    // vector.  This would be zero in exact arithmetic, but error can accumulate
    // over time in finite precision.
    const real offset = (beta * blockSums1.qStep - blockSums2.precGrad) / NUM_ELECTRODE_PARTICLES;
    for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
        qStep[ii] -= offset;
    }
#endif

    // Prepare for the next derivative calculation.
    for (int ii = GLOBAL_ID; ii < NUM_ELECTRODE_PARTICLES; ii += GLOBAL_SIZE) {
        electrodeCharges[ii] = qStep[ii];
    }
}
