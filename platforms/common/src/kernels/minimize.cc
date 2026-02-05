#define LS_DOT_START 0
#define LS_DOT 1
#define LS_ENERGY 2
#define LS_STEP 3

#define LS_FAIL 0
#define LS_SUCCEED 1
#define LS_CONTINUE 2

#define WARP_SIZE 32

#ifdef USE_MIXED_PRECISION
    // When mixed is double precision, use double precision sqrt and fabs to
    // avoid overflow.
    #define SQRT_MIXED sqrt
    #define FABS_MIXED fabs
#else
    #define SQRT_MIXED SQRT
    #define FABS_MIXED FABS
#endif

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

DEVICE void resetLineSearchData(GLOBAL volatile mixed* lineSearchData) {
    if (GLOBAL_ID == 0) {
        lineSearchData[LS_DOT_START] = 0;
        lineSearchData[LS_STEP] = 1;
    }
}

DEVICE mixed reduceAdd(mixed value, LOCAL_ARG volatile mixed* temp) {
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
        if (thread + step < LOCAL_SIZE && thread % (2 * step) == 0) {
            temp[thread] += temp[thread + step];
        }
        SYNC_WARPS;
    }
    for (int step = WARP_SIZE / 2; step < LOCAL_SIZE; step <<= 1) {
        if (thread + step < LOCAL_SIZE && thread % (2 * step) == 0) {
            temp[thread] += temp[thread + step];
        }
        SYNC_THREADS;
    }
#endif
    return temp[0];
}

DEVICE mixed reduceMax(mixed value, LOCAL_ARG volatile mixed* temp) {
    const int thread = LOCAL_ID;
    SYNC_THREADS;
#ifdef WARP_SHUFFLE_DOWN
    const int warpCount = LOCAL_SIZE / WARP_SIZE;
    const int warp = thread / WARP_SIZE;
    const int lane = thread % WARP_SIZE;
    for (int step = WARP_SIZE / 2; step > 0; step >>= 1) {
        value = max(value, WARP_SHUFFLE_DOWN(value, step));
    }
    if (!lane) {
        temp[warp] = value;
    }
    SYNC_THREADS;
    if (!warp) {
        value = lane < warpCount ? temp[lane] : 0;
        for (int step = WARP_SIZE / 2; step > 0; step >>= 1) {
            value = max(value, WARP_SHUFFLE_DOWN(value, step));
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
        if (thread + step < LOCAL_SIZE && thread % (2 * step) == 0) {
            temp[thread] = max(temp[thread], temp[thread + step]);
        }
        SYNC_WARPS;
    }
    for (int step = WARP_SIZE / 2; step < LOCAL_SIZE; step <<= 1) {
        if (thread + step < LOCAL_SIZE && thread % (2 * step) == 0) {
            temp[thread] = max(temp[thread], temp[thread + step]);
        }
        SYNC_THREADS;
    }
#endif
    return temp[0];
}

DEVICE void atomicAddMixed(GLOBAL mixed* RESTRICT target, const mixed value) {
#if defined(__CUDA_ARCH__) || defined(USE_HIP)
    atomicAdd(target, value);
#elif defined(USE_MIXED_PRECISION) || defined(USE_DOUBLE_PRECISION)
    unsigned long old = as_ulong(*target);
    unsigned long check;
    do {
        check = old;
        old = atom_cmpxchg((GLOBAL unsigned long*)target, check, as_ulong(value + as_double(check)));
    } while(check != old);
#else
    unsigned int old = as_uint(*target);
    unsigned int check;
    do {
        check = old;
        old = atomic_cmpxchg((GLOBAL unsigned int*)target, check, as_uint(value + as_float(check)));
    } while(check != old);
#endif
}

KERNEL void recordInitialPos(
    GLOBAL const real4* RESTRICT posq,
    GLOBAL const int* RESTRICT order,
    GLOBAL mixed* RESTRICT xInit,
    GLOBAL mixed* RESTRICT x,
    const int numParticles
#ifdef USE_MIXED_PRECISION
    , GLOBAL const real4* RESTRICT posqCorrection
#endif
) {
    for (int i = GLOBAL_ID; i < numParticles; i += GLOBAL_SIZE) {
        int offset = 3 * order[i];
#ifdef USE_MIXED_PRECISION
        real4 pos1 = posq[i];
        real4 pos2 = posqCorrection[i];
        mixed3 pos = make_mixed3(pos1.x + (mixed) pos2.x, pos1.y + (mixed) pos2.y, pos1.z + (mixed) pos2.z);
#else
        real4 pos = posq[i];
#endif
        xInit[offset] = x[offset] = pos.x;
        xInit[offset + 1] = x[offset + 1] = pos.y;
        xInit[offset + 2] = x[offset + 2] = pos.z;
    }
}

KERNEL void restorePos(
    GLOBAL real4* RESTRICT posq,
    GLOBAL const int* RESTRICT order,
    GLOBAL const mixed* RESTRICT x,
    GLOBAL mixed* RESTRICT returnValue,
    const int numParticles
#ifdef USE_MIXED_PRECISION
    , GLOBAL real4* RESTRICT posqCorrection
#endif
) {
    for (int i = GLOBAL_ID; i < numParticles; i += GLOBAL_SIZE) {
        int offset = 3 * order[i];
        posq[i].x = x[offset];
        posq[i].y = x[offset + 1];
        posq[i].z = x[offset + 2];
#ifdef USE_MIXED_PRECISION
        posqCorrection[i].x = x[offset] - (real) x[offset];
        posqCorrection[i].y = x[offset + 1] - (real) x[offset + 1];
        posqCorrection[i].z = x[offset + 2] - (real) x[offset + 2];
#endif
    }

    // Reset in case we will accumulate constraint energies in returnValue.

    if (GLOBAL_ID == 0) {
        *returnValue = 0;
    }
}

KERNEL void convertForces(
    GLOBAL const mixed4* RESTRICT velm,
    GLOBAL const mm_long* RESTRICT forceBuffer,
    GLOBAL const int* RESTRICT order,
    GLOBAL mixed* RESTRICT grad,
    GLOBAL mixed* RESTRICT returnValue,
    const int numParticles,
    const int numPadded
) {
    const mm_long limit = 0x4000000000000000;
    const mixed scale = -1 / (mixed) 0x100000000;
    const mixed nanValue = NAN;

    for (int i = GLOBAL_ID; i < numParticles; i += GLOBAL_SIZE) {
        int offset = 3 * order[i];
        if (velm[i].w == 0) {
            grad[offset] = 0;
            grad[offset + 1] = 0;
            grad[offset + 2] = 0;
        }
        else {
            mm_long fx = forceBuffer[i];
            mm_long fy = forceBuffer[i + numPadded];
            mm_long fz = forceBuffer[i + 2 * numPadded];
            if (fx < -limit || fx > limit || fy < -limit || fy > limit || fz < -limit || fz > limit) {
                *returnValue = nanValue;
            }
            grad[offset] = scale * fx;
            grad[offset + 1] = scale * fy;
            grad[offset + 2] = scale * fz;
        }
    }
}

KERNEL void getConstraintEnergyForces(
    GLOBAL mm_ulong* RESTRICT forceBuffer,
    GLOBAL const int* RESTRICT order,
    GLOBAL const int2* RESTRICT constraintIndices,
    GLOBAL const mixed* RESTRICT constraintDistances,
    GLOBAL const mixed* RESTRICT x,
    GLOBAL mixed* RESTRICT returnValue,
    const int numPadded,
    const int numConstraints,
    const mixed kRestraint
) {
    LOCAL volatile mixed temp[TEMP_SIZE];

    const mixed scale = 0x100000000;

    mixed energy = 0;
    for (int i = GLOBAL_ID; i < numConstraints; i += GLOBAL_SIZE) {
        int2 indices = constraintIndices[i];
        mixed distance = constraintDistances[i];
        int offset1 = 3 * order[indices.x];
        int offset2 = 3 * order[indices.y];
        mixed3 delta = make_mixed3(x[offset2] - x[offset1], x[offset2 + 1] - x[offset1 + 1], x[offset2 + 2] - x[offset1 + 2]);
        mixed r2 = dot(delta, delta);
        mixed r = SQRT_MIXED(r2);
        delta *= 1 / r;
        mixed dr = r - distance;
        mixed kdr = kRestraint * dr;
        energy += (mixed) 0.5 * kdr * dr;
        mm_long fx = (mm_long) (kdr * scale * delta.x);
        mm_long fy = (mm_long) (kdr * scale * delta.y);
        mm_long fz = (mm_long) (kdr * scale * delta.z);
        ATOMIC_ADD(&forceBuffer[indices.x], (mm_ulong) fx);
        ATOMIC_ADD(&forceBuffer[indices.x + numPadded], (mm_ulong) fy);
        ATOMIC_ADD(&forceBuffer[indices.x + 2 * numPadded], (mm_ulong) fz);
        ATOMIC_ADD(&forceBuffer[indices.y], (mm_ulong) -fx);
        ATOMIC_ADD(&forceBuffer[indices.y + numPadded], (mm_ulong) -fy);
        ATOMIC_ADD(&forceBuffer[indices.y + 2 * numPadded], (mm_ulong) -fz);
    }
    energy = reduceAdd(energy, temp);

    if (LOCAL_ID == 0) {
        atomicAddMixed(returnValue, energy);
    }
}

KERNEL void getConstraintError(
    GLOBAL const int2* RESTRICT constraintIndices,
    GLOBAL const mixed* RESTRICT constraintDistances,
    GLOBAL const mixed* RESTRICT x,
    GLOBAL mixed* RESTRICT returnValue,
    const int numConstraints
) {
    // This kernel is expected to be executed in a single thread block.

    LOCAL volatile mixed temp[TEMP_SIZE];

    mixed maxError = 0;
    for (int i = LOCAL_ID; i < numConstraints; i += LOCAL_SIZE) {
        int2 indices = constraintIndices[i];
        mixed distance = constraintDistances[i];
        mixed3 delta = make_mixed3(x[3 * indices.y] - x[3 * indices.x], x[3 * indices.y + 1] - x[3 * indices.x + 1], x[3 * indices.y + 2] - x[3 * indices.x + 2]);
        mixed r = SQRT_MIXED(dot(delta, delta));
        maxError = max(maxError, FABS_MIXED(r - distance) / distance);
    }
    maxError = reduceMax(maxError, temp);

    if (LOCAL_ID == 0) {
        returnValue[0] = maxError;
    }
}

KERNEL void initializeDir(
    GLOBAL const mixed* RESTRICT grad,
    GLOBAL mixed* RESTRICT dir,
    GLOBAL const mixed* RESTRICT gradNorm,
    GLOBAL mixed* RESTRICT lineSearchData,
    const int numVariables
) {
    const real scale = -1 / gradNorm[0];

    for (int i = GLOBAL_ID; i < numVariables; i += GLOBAL_SIZE) {
        dir[i] = scale * grad[i];
    }

    // Prepare for the first line search of the optimization.

    resetLineSearchData(lineSearchData);
}

KERNEL void gradNormPart1(
    GLOBAL const mixed* RESTRICT grad,
    GLOBAL mixed* RESTRICT reduceBuffer,
    const int numVariables
) {
    LOCAL volatile mixed temp[TEMP_SIZE];

    mixed norm = 0;
    for (int i = GLOBAL_ID; i < numVariables; i += GLOBAL_SIZE) {
        norm += grad[i] * grad[i];
    }
    norm = reduceAdd(norm, temp);

    if (LOCAL_ID == 0) {
        reduceBuffer[GROUP_ID] = norm;
    }
}


KERNEL void gradNormPart2(
    GLOBAL const mixed* RESTRICT grad,
    GLOBAL const mixed* RESTRICT reduceBuffer,
    GLOBAL mixed* RESTRICT gradNorm,
    const int numVariables,
    const int numVariableBlocks
) {
    // This kernel is expected to be executed in a single thread block.

    LOCAL volatile mixed temp[TEMP_SIZE];

    mixed norm = 0;
    for (int i = LOCAL_ID; i < numVariableBlocks; i += LOCAL_SIZE) {
        norm += reduceBuffer[i];
    }
    norm = reduceAdd(norm, temp);

    if (isfinite(norm)) {
        if (LOCAL_ID == 0) {
            gradNorm[0] = SQRT_MIXED(norm);
        }
        return;
    }

    // The square norm overflowed, so take a slow fallback path.

    mixed gradScale = 1;
    for (int i = LOCAL_ID; i < numVariables; i += LOCAL_SIZE) {
        gradScale = max(gradScale, FABS_MIXED(grad[i]));
    }
    gradScale = reduceMax(gradScale, temp);
    mixed gradScaleInv = 1 / gradScale;

    mixed scaledNorm = 0;
    for (int i = LOCAL_ID; i < numVariables; i += LOCAL_SIZE) {
        mixed gradScaled = gradScaleInv * grad[i];
        scaledNorm += gradScaled * gradScaled;
    }
    scaledNorm = reduceAdd(scaledNorm, temp);

    if (LOCAL_ID == 0) {
        gradNorm[0] = gradScale * SQRT_MIXED(scaledNorm);
    }
}

KERNEL void getDiff(
    GLOBAL const mixed* RESTRICT x,
    GLOBAL const mixed* RESTRICT xPrev,
    GLOBAL const mixed* RESTRICT grad,
    GLOBAL const mixed* RESTRICT gradPrev,
    GLOBAL mixed* RESTRICT xDiff,
    GLOBAL mixed* RESTRICT gradDiff,
    GLOBAL mixed* RESTRICT reduceBuffer,
    GLOBAL int* RESTRICT returnFlag,
    GLOBAL const mixed* RESTRICT gradNorm,
    const int numVariables,
    const int numVariableBlocks,
    const mixed tolerance,
    const int end
) {
    LOCAL volatile mixed temp[TEMP_SIZE];

    // If the convergence condition is satisfied, update returnFlag.  After this
    // kernel has run, it is safe to start downloading returnFlag.

    if (*gradNorm <= tolerance) {
        if (GLOBAL_ID == 0) {
            *returnFlag = 1;
        }
        return;
    }
    else {
        if (GLOBAL_ID == 0) {
            *returnFlag = 0;
        }
    }

    const int endOffset = numVariables * end;

    mixed xGrad = 0;
    mixed gradGrad = 0;
    for (int i = GLOBAL_ID; i < numVariables; i += GLOBAL_SIZE) {
        xDiff[endOffset + i] = x[i] - xPrev[i];
        gradDiff[endOffset + i] = grad[i] - gradPrev[i];
        xGrad += xDiff[endOffset + i] * gradDiff[endOffset + i];
        gradGrad += gradDiff[endOffset + i] * gradDiff[endOffset + i];
    }
    xGrad = reduceAdd(xGrad, temp);
    gradGrad = reduceAdd(gradGrad, temp);

    if (LOCAL_ID == 0) {
        reduceBuffer[GROUP_ID] = xGrad;
        reduceBuffer[GROUP_ID + numVariableBlocks] = gradGrad;
    }
}

KERNEL void getScale(
    GLOBAL mixed* RESTRICT alpha,
    GLOBAL mixed* RESTRICT scale,
    GLOBAL const mixed* RESTRICT xDiff,
    GLOBAL const mixed* RESTRICT gradDiff,
    GLOBAL const mixed* RESTRICT reduceBuffer,
    GLOBAL const int* RESTRICT returnFlag,
    GLOBAL mixed* RESTRICT returnValue,
    const int numVariables,
    const int numVariableBlocks,
    const int end
) {
    // This kernel is expected to be executed in a single thread block.

    LOCAL volatile mixed temp[TEMP_SIZE];

    if (*returnFlag) {
        return;
    }

    const int endOffset = numVariables * end;

    mixed xGrad = 0;
    mixed gradGrad = 0;
    for (int i = LOCAL_ID; i < numVariableBlocks; i += LOCAL_SIZE) {
        xGrad += reduceBuffer[i];
        gradGrad += reduceBuffer[i + numVariableBlocks];
    }
    xGrad = reduceAdd(xGrad, temp);
    gradGrad = reduceAdd(gradGrad, temp);

    // Clear all values of alpha (plus the extra slot at alpha[NUM_VECTORS])
    // to prepare them for use in reductions by subsequent kernels.
    for (int i = LOCAL_ID; i <= NUM_VECTORS; i += LOCAL_SIZE) {
        alpha[i] = 0;
    }

    if (isfinite(gradGrad)) {
        if (LOCAL_ID == 0) {
            scale[end] = xGrad;
            returnValue[0] = xGrad / gradGrad;
        }
        return;
    }

    // The square norm overflowed, so take a slow fallback path.

    mixed gradScale = 1;
    for (int i = LOCAL_ID; i < numVariables; i += LOCAL_SIZE) {
        gradScale = max(gradScale, FABS_MIXED(gradDiff[endOffset + i]));
    }
    gradScale = reduceMax(gradScale, temp);
    mixed gradScaleInv = 1 / gradScale;

    mixed xGradScaled = 0;
    mixed gradGradScaled = 0;
    for (int i = LOCAL_ID; i < numVariables; i += LOCAL_SIZE) {
        mixed gradDiffScaled = gradScaleInv * gradDiff[endOffset + i];
        xGradScaled += xDiff[endOffset + i] * gradDiffScaled;
        gradGradScaled += gradDiffScaled * gradDiffScaled;
    }
    xGradScaled = reduceAdd(xGradScaled, temp);
    gradGradScaled = reduceAdd(gradGradScaled, temp);

    if (LOCAL_ID == 0) {
        scale[end] = xGrad;
        returnValue[0] = xGradScaled / (gradScale * gradGradScaled);
    }
}

KERNEL void reinitializeDir(
    GLOBAL const mixed* RESTRICT grad,
    GLOBAL mixed* RESTRICT dir,
    GLOBAL mixed* RESTRICT alpha,
    GLOBAL const mixed* RESTRICT scale,
    GLOBAL const mixed* RESTRICT xDiff,
    GLOBAL const int* RESTRICT returnFlag,
    const int numVariables,
    const int vectorIndex
) {
    LOCAL volatile mixed temp[TEMP_SIZE];

    if (*returnFlag) {
        return;
    }

    const int indexOffset = numVariables * vectorIndex;

    mixed vectorAlpha = 0;
    for (int i = GLOBAL_ID; i < numVariables; i += GLOBAL_SIZE) {
        dir[i] = -grad[i];
        vectorAlpha -= grad[i] * xDiff[indexOffset + i];
    }
    vectorAlpha = reduceAdd(vectorAlpha, temp);

    if (LOCAL_ID == 0) {
        atomicAddMixed(&alpha[vectorIndex], vectorAlpha / scale[vectorIndex]);
    }
}

KERNEL void updateDirAlpha(
    GLOBAL mixed* dir,
    GLOBAL mixed* RESTRICT alpha,
    GLOBAL const mixed* RESTRICT scale,
    GLOBAL const mixed* RESTRICT xDiff,
    GLOBAL const mixed* RESTRICT gradDiff,
    GLOBAL const int* RESTRICT returnFlag,
    const int numVariables,
    const int vectorIndex1
) {
    LOCAL volatile mixed temp[TEMP_SIZE];

    if (*returnFlag) {
        return;
    }

    const int vectorIndex2 = (vectorIndex1 ? vectorIndex1 : NUM_VECTORS) - 1;
    const int indexOffset1 = numVariables * vectorIndex1;
    const int indexOffset2 = numVariables * vectorIndex2;
    const mixed vectorScale = alpha[vectorIndex1];

    mixed vectorAlpha = 0;
    for (int i = GLOBAL_ID; i < numVariables; i += GLOBAL_SIZE) {
        dir[i] -= vectorScale * gradDiff[indexOffset1 + i];
        vectorAlpha += dir[i] * xDiff[indexOffset2 + i];
    }
    vectorAlpha = reduceAdd(vectorAlpha, temp);

    if (LOCAL_ID == 0) {
        atomicAddMixed(&alpha[vectorIndex2], vectorAlpha / scale[vectorIndex2]);
    }
}

KERNEL void scaleDir(
    GLOBAL mixed* dir,
    GLOBAL mixed* RESTRICT alpha,
    GLOBAL const mixed* RESTRICT scale,
    GLOBAL const mixed* RESTRICT gradDiff,
    GLOBAL const int* RESTRICT returnFlag,
    GLOBAL const mixed* RESTRICT returnValue,
    const int numVariables,
    const int vectorIndex
) {
    LOCAL volatile mixed temp[TEMP_SIZE];

    if (*returnFlag) {
        return;
    }

    const int indexOffset = numVariables * vectorIndex;
    const mixed innerScale = alpha[vectorIndex];
    const mixed outerScale = returnValue[0];

    mixed vectorBeta = 0;
    for (int i = GLOBAL_ID; i < numVariables; i += GLOBAL_SIZE) {
        dir[i] = (dir[i] - innerScale * gradDiff[indexOffset + i]) * outerScale;
        vectorBeta += dir[i] * gradDiff[indexOffset + i];
    }
    vectorBeta = reduceAdd(vectorBeta, temp);

    if (LOCAL_ID == 0) {
        // Store the result in an extra slot alpha[NUM_VECTORS] instead of using
        // alpha[vectorIndex] since other blocks might still need to read it.

        atomicAddMixed(&alpha[NUM_VECTORS], (GLOBAL_ID == 0 ? innerScale : 0) - vectorBeta / scale[vectorIndex]);
    }
}

KERNEL void updateDirBeta(
    GLOBAL mixed* dir,
    GLOBAL mixed* RESTRICT alpha,
    GLOBAL const mixed* RESTRICT scale,
    GLOBAL const mixed* RESTRICT xDiff,
    GLOBAL const mixed* RESTRICT gradDiff,
    GLOBAL const int* RESTRICT returnFlag,
    const int numVariables,
    const int vectorIndex1,
    const int vectorIndexAlpha
) {
    LOCAL volatile mixed temp[TEMP_SIZE];

    if (*returnFlag) {
        return;
    }

    const int vectorIndex2 = vectorIndex1 == NUM_VECTORS - 1 ? 0 : vectorIndex1 + 1;
    const int indexOffset1 = numVariables * vectorIndex1;
    const int indexOffset2 = numVariables * vectorIndex2;
    const mixed vectorScale = alpha[vectorIndexAlpha];

    mixed vectorBeta = 0;
    for (int i = GLOBAL_ID; i < numVariables; i += GLOBAL_SIZE) {
        dir[i] += vectorScale * xDiff[indexOffset1 + i];
        vectorBeta += dir[i] * gradDiff[indexOffset2 + i];
    }
    vectorBeta = reduceAdd(vectorBeta, temp);

    if (LOCAL_ID == 0) {
        atomicAddMixed(&alpha[vectorIndex2], -vectorBeta / scale[vectorIndex2]);
    }
}

KERNEL void updateDirFinal(
    GLOBAL mixed* dir,
    GLOBAL const mixed* RESTRICT alpha,
    GLOBAL const mixed* RESTRICT xDiff,
    GLOBAL const int* RESTRICT returnFlag,
    GLOBAL mixed* RESTRICT lineSearchData,
    const int numVariables,
    const int vectorIndex,
    const int vectorIndexAlpha
) {
    LOCAL volatile mixed temp[TEMP_SIZE];

    if (*returnFlag) {
        return;
    }

    const int indexOffset = numVariables * vectorIndex;
    const mixed vectorScale = alpha[vectorIndexAlpha];

    for (int i = GLOBAL_ID; i < numVariables; i += GLOBAL_SIZE) {
        dir[i] += vectorScale * xDiff[indexOffset + i];
    }

    // Prepare for a line search in the next iteration.

    resetLineSearchData(lineSearchData);
}

KERNEL void lineSearchSetup(
    GLOBAL const mixed* RESTRICT x,
    GLOBAL mixed* RESTRICT xPrev,
    GLOBAL const mixed* RESTRICT grad,
    GLOBAL mixed* RESTRICT gradPrev,
    GLOBAL const mixed* RESTRICT dir,
    GLOBAL int* RESTRICT returnFlag,
    GLOBAL mixed* RESTRICT lineSearchData,
    const int numVariables,
    const mixed energyStart
) {
    LOCAL volatile mixed temp[TEMP_SIZE];

    for (int i = GLOBAL_ID; i < numVariables; i += GLOBAL_SIZE) {
        xPrev[i] = x[i];
    }

    for (int i = GLOBAL_ID; i < numVariables; i += GLOBAL_SIZE) {
        gradPrev[i] = grad[i];
    }

    mixed result = 0;
    for (int i = GLOBAL_ID; i < numVariables; i += GLOBAL_SIZE) {
        result += grad[i] * dir[i];
    }
    result = reduceAdd(result, temp);

    if (LOCAL_ID == 0) {
        atomicAddMixed(&lineSearchData[LS_DOT_START], result);
    }

    if (GLOBAL_ID == 0) {
        *returnFlag = LS_CONTINUE;
        lineSearchData[LS_ENERGY] = energyStart;
    }
}

KERNEL void lineSearchStep(
    GLOBAL mixed* RESTRICT x,
    GLOBAL const mixed* RESTRICT xPrev,
    GLOBAL mixed* RESTRICT grad,
    GLOBAL const mixed* RESTRICT gradPrev,
    GLOBAL const mixed* RESTRICT dir,
    GLOBAL mixed* RESTRICT reduceBuffer,
    GLOBAL int* RESTRICT returnFlag,
    GLOBAL mixed* RESTRICT lineSearchData,
    const int numVariables
) {
    LOCAL volatile mixed temp[TEMP_SIZE];

    if (*returnFlag == LS_FAIL) {
        // The line search failed.  Put back the old coordinates and quit.

        for (int i = GLOBAL_ID; i < numVariables; i += GLOBAL_SIZE) {
            x[i] = xPrev[i];
        }

        for (int i = GLOBAL_ID; i < numVariables; i += GLOBAL_SIZE) {
            grad[i] = gradPrev[i];
        }

        return;
    }

    if (*returnFlag == LS_SUCCEED) {
        // The strong Wolfe condition was satisfied on the last iteration.
        // Instead of taking another step, start evaluating the gradient.

        mixed norm = 0;
        for (int i = GLOBAL_ID; i < numVariables; i += GLOBAL_SIZE) {
            norm += grad[i] * grad[i];
        }
        norm = reduceAdd(norm, temp);

        if (LOCAL_ID == 0) {
            reduceBuffer[GROUP_ID] = norm;
        }

        return;
    }

    // We are either getting started or going around for another normal line
    // search iteration, so check to see whether the initial search direction is
    // bad, or the updated step is bad.

    if (lineSearchData[LS_DOT_START] > 0 || lineSearchData[LS_STEP] < LBFGS_MIN_STEP || lineSearchData[LS_STEP] > LBFGS_MAX_STEP) {
        if (GLOBAL_ID == 0) {
            *returnFlag = LS_FAIL;
        }
        return;
    }

    const mixed step = lineSearchData[LS_STEP];

    for (int i = GLOBAL_ID; i < numVariables; i += GLOBAL_SIZE) {
        x[i] = xPrev[i] + step * dir[i];
    }

    if (GLOBAL_ID == 0) {
        lineSearchData[LS_DOT] = 0;
    }
}

KERNEL void lineSearchDot(
    GLOBAL const mixed* RESTRICT grad,
    GLOBAL const mixed* RESTRICT dir,
    GLOBAL mixed* RESTRICT lineSearchData,
    GLOBAL int* RESTRICT returnFlag,
    const int numVariables,
    const mixed energy
) {
    LOCAL volatile mixed temp[TEMP_SIZE];

    if (*returnFlag == LS_FAIL) {
        return;
    }

    // The energy may be such that we don't need to do a dot product and can
    // immediately decide to scale the step, so mark this case with LS_SUCCEED.
    // This will be checked in the following kernel.

    if (!isfinite(energy) || energy > lineSearchData[LS_ENERGY] + lineSearchData[LS_STEP] * LBFGS_FTOL * lineSearchData[LS_DOT_START]) {
        if (GLOBAL_ID == 0) {
            *returnFlag = LS_SUCCEED;
        }
        return;
    }

    mixed result = 0;
    for (int i = GLOBAL_ID; i < numVariables; i += GLOBAL_SIZE) {
        result += grad[i] * dir[i];
    }
    result = reduceAdd(result, temp);

    if (LOCAL_ID == 0) {
        atomicAddMixed(&lineSearchData[LS_DOT], result);
    }
}

KERNEL void lineSearchContinue(
    GLOBAL int* RESTRICT returnFlag,
    GLOBAL mixed* RESTRICT lineSearchData
) {
    // This kernel should run a single thread and just lets us update the step
    // and decide how to continue the line search without returning to the CPU.

    if (GLOBAL_ID != 0 || *returnFlag == LS_FAIL) {
        return;
    }

    // If the last kernel set LS_SUCCEED, the energy was out of range and we
    // should scale the step back and go around again.

    if (*returnFlag == LS_SUCCEED) {
        lineSearchData[LS_STEP] *= LBFGS_SCALE_DOWN;
        *returnFlag = LS_CONTINUE;
        return;
    }

    // If we are here, we should be reading LS_CONTINUE that was set by
    // lineSearchSetup, so either scale the step or set LS_SUCCEED if the line
    // search is successful and we should exit.

    if (lineSearchData[LS_DOT] < LBFGS_WOLFE * lineSearchData[LS_DOT_START]) {
        lineSearchData[LS_STEP] *= LBFGS_SCALE_UP;
    }
    else if(lineSearchData[LS_DOT] > -LBFGS_WOLFE * lineSearchData[LS_DOT_START]) {
        lineSearchData[LS_STEP] *= LBFGS_SCALE_DOWN;
    }
    else {
        *returnFlag = LS_SUCCEED;
    }
}
