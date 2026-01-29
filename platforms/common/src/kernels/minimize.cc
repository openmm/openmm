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
        if(thread + step < LOCAL_SIZE && thread % (2 * step) == 0) {
            temp[thread] = max(temp[thread], temp[thread + step]);
        }
        SYNC_WARPS;
    }
    for (int step = WARP_SIZE / 2; step < LOCAL_SIZE; step <<= 1) {
        if(thread + step < LOCAL_SIZE && thread % (2 * step) == 0) {
            temp[thread] = max(temp[thread], temp[thread + step]);
        }
        SYNC_THREADS;
    }
#endif
    return temp[0];
}

KERNEL void recordInitialPos(
    GLOBAL const real4* RESTRICT posq,
    GLOBAL const int* RESTRICT order,
    GLOBAL mixed* RESTRICT xInit,
    GLOBAL mixed* RESTRICT x
#ifdef USE_MIXED_PRECISION
    , GLOBAL const real4* RESTRICT posqCorrection
#endif
) {
    for (int i = GLOBAL_ID; i < NUM_PARTICLES; i += GLOBAL_SIZE) {
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
    GLOBAL const mixed* RESTRICT x
#ifdef USE_MIXED_PRECISION
    , GLOBAL real4* RESTRICT posqCorrection
#endif
) {
    for (int i = GLOBAL_ID; i < NUM_PARTICLES; i += GLOBAL_SIZE) {
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
}

KERNEL void convertForces(
    GLOBAL const mixed4* RESTRICT velm,
    GLOBAL const mm_long* RESTRICT forceBuffer,
    GLOBAL const int* RESTRICT order,
    GLOBAL mixed* RESTRICT grad,
    GLOBAL int* RESTRICT returnFlag
) {
    const mm_long limit = 0x4000000000000000;
    const mixed scale = -1 / (mixed) 0x100000000;

    for (int i = GLOBAL_ID; i < NUM_PARTICLES; i += GLOBAL_SIZE) {
        if(velm[i].w == 0) {
            continue;
        }

        mm_long fx = forceBuffer[i];
        mm_long fy = forceBuffer[i + NUM_PADDED];
        mm_long fz = forceBuffer[i + 2 * NUM_PADDED];
        if (fx < -limit || fx > limit || fy < -limit || fy > limit || fz < -limit || fz > limit) {
            *returnFlag = 1;
        }
        int offset = 3 * order[i];
        grad[offset] = scale * fx;
        grad[offset + 1] = scale * fy;
        grad[offset + 2] = scale * fz;
    }
}

KERNEL void getConstraintEnergyForces(
    GLOBAL mm_ulong* RESTRICT forceBuffer,
    GLOBAL const int* RESTRICT order,
    GLOBAL const int2* RESTRICT constraintIndices,
    GLOBAL const mixed* RESTRICT constraintDistances,
    GLOBAL const mixed* RESTRICT x,
    GLOBAL mixed* RESTRICT returnValue,
    const mixed kRestraint
) {
    // This kernel is expected to be executed in a single thread block.

    LOCAL volatile mixed temp[TEMP_SIZE];

    const mixed scale = 0x100000000;

    mixed energy = 0;
    for (int i = LOCAL_ID; i < NUM_CONSTRAINTS; i += LOCAL_SIZE) {
        int2 indices = constraintIndices[i];
        mixed distance = constraintDistances[i];
        int offset1 = 3 * order[indices.x];
        int offset2 = 3 * order[indices.y];
        mixed3 delta = make_mixed3(x[offset2] - x[offset1], x[offset2 + 1] - x[offset1 + 1], x[offset2 + 2] - x[offset1 + 2]);
        mixed r2 = dot(delta, delta);
        mixed r = SQRT(r2);
        delta *= 1 / r;
        mixed dr = r - distance;
        mixed kdr = kRestraint * dr;
        energy += (mixed) 0.5 * kdr * dr;
        mm_long fx = (mm_long) (kdr * scale * delta.x);
        mm_long fy = (mm_long) (kdr * scale * delta.y);
        mm_long fz = (mm_long) (kdr * scale * delta.z);
        ATOMIC_ADD(&forceBuffer[indices.x], (mm_ulong) fx);
        ATOMIC_ADD(&forceBuffer[indices.x + NUM_PADDED], (mm_ulong) fy);
        ATOMIC_ADD(&forceBuffer[indices.x + 2 * NUM_PADDED], (mm_ulong) fz);
        ATOMIC_ADD(&forceBuffer[indices.y], (mm_ulong) -fx);
        ATOMIC_ADD(&forceBuffer[indices.y + NUM_PADDED], (mm_ulong) -fy);
        ATOMIC_ADD(&forceBuffer[indices.y + 2 * NUM_PADDED], (mm_ulong) -fz);
    }
    energy = reduceAdd(energy, temp);

    if (LOCAL_ID == 0) {
        returnValue[0] = energy;
    }
}

KERNEL void getConstraintError(
    GLOBAL const int2* RESTRICT constraintIndices,
    GLOBAL const mixed* RESTRICT constraintDistances,
    GLOBAL const mixed* RESTRICT x,
    GLOBAL mixed* RESTRICT returnValue
) {
    // This kernel is expected to be executed in a single thread block.

    LOCAL volatile mixed temp[TEMP_SIZE];

    mixed maxError = 0;
    for (int i = LOCAL_ID; i < NUM_CONSTRAINTS; i += LOCAL_SIZE) {
        int2 indices = constraintIndices[i];
        mixed distance = constraintDistances[i];
        mixed3 delta = make_mixed3(x[3 * indices.y] - x[3 * indices.x], x[3 * indices.y + 1] - x[3 * indices.x + 1], x[3 * indices.y + 2] - x[3 * indices.x + 2]);
        mixed r = SQRT(dot(delta, delta));
        maxError = max(maxError, FABS(r - distance) / distance);
    }
    maxError = reduceMax(maxError, temp);

    if (LOCAL_ID == 0) {
        returnValue[0] = maxError;
    }
}

KERNEL void initializeDir(
    GLOBAL const mixed* RESTRICT grad,
    GLOBAL mixed* RESTRICT dir
) {
    for (int i = GLOBAL_ID; i < NUM_VARIABLES; i += GLOBAL_SIZE) {
        dir[i] = -grad[i];
    }
}

KERNEL void gradNorm(
    GLOBAL const mixed* RESTRICT grad,
    GLOBAL mixed* RESTRICT returnValue
) {
    // This kernel is expected to be executed in a single thread block.

    LOCAL volatile mixed temp[TEMP_SIZE];

    mixed norm = 0;
    for (int i = LOCAL_ID; i < NUM_VARIABLES; i += LOCAL_SIZE) {
        norm += grad[i] * grad[i];
    }
    norm = reduceAdd(norm, temp);

    if (LOCAL_ID == 0) {
        returnValue[0] = SQRT(norm);
    }
}

KERNEL void getDiff(
    GLOBAL const mixed* RESTRICT x,
    GLOBAL const mixed* RESTRICT xPrev,
    GLOBAL const mixed* RESTRICT grad,
    GLOBAL const mixed* RESTRICT gradPrev,
    GLOBAL mixed* RESTRICT xDiff,
    GLOBAL mixed* RESTRICT gradDiff,
    const int end
) {
    const int endOffset = NUM_VARIABLES * end;

    for (int i = GLOBAL_ID; i < NUM_VARIABLES; i += GLOBAL_SIZE) {
        xDiff[endOffset + i] = x[i] - xPrev[i];
    }

    for (int i = GLOBAL_ID; i < NUM_VARIABLES; i += GLOBAL_SIZE) {
        gradDiff[endOffset + i] = grad[i] - gradPrev[i];
    }
}

KERNEL void getScale(
    GLOBAL const mixed* RESTRICT xDiff,
    GLOBAL const mixed* RESTRICT gradDiff,
    GLOBAL mixed* RESTRICT scale,
    GLOBAL mixed* RESTRICT returnValue,
    const int end
) {
    // This kernel is expected to be executed in a single thread block.

    LOCAL volatile mixed temp[TEMP_SIZE];

    const int endOffset = NUM_VARIABLES * end;

    mixed xGrad = 0;
    mixed gradGrad = 0;
    for (int i = LOCAL_ID; i < NUM_VARIABLES; i += LOCAL_SIZE) {
        xGrad += xDiff[endOffset + i] * gradDiff[endOffset + i];
        gradGrad += gradDiff[endOffset + i] * gradDiff[endOffset + i];
    }
    xGrad = reduceAdd(xGrad, temp);
    gradGrad = reduceAdd(gradGrad, temp);

    if (LOCAL_ID == 0) {
        scale[end] = xGrad;
        returnValue[0] = gradGrad;
    }
}

KERNEL void getAlpha(
    GLOBAL const mixed* RESTRICT dir,
    GLOBAL mixed* RESTRICT alpha,
    GLOBAL const mixed* RESTRICT scale,
    GLOBAL const mixed* RESTRICT xDiff,
    const int vectorIndex
) {
    // This kernel is expected to be executed in a single thread block.

    LOCAL volatile mixed temp[TEMP_SIZE];

    const int indexOffset = NUM_VARIABLES * vectorIndex;

    mixed vectorAlpha = 0;
    for (int i = LOCAL_ID; i < NUM_VARIABLES; i += LOCAL_SIZE) {
        vectorAlpha += dir[i] * xDiff[indexOffset + i];
    }
    vectorAlpha = reduceAdd(vectorAlpha, temp);

    if (LOCAL_ID == 0) {
        alpha[vectorIndex] = vectorAlpha / scale[vectorIndex];
    }
}

KERNEL void getBeta(
    GLOBAL const mixed* RESTRICT dir,
    GLOBAL const mixed* RESTRICT scale,
    GLOBAL const mixed* RESTRICT gradDiff,
    GLOBAL mixed* RESTRICT returnValue,
    const int vectorIndex
) {
    // This kernel is expected to be executed in a single thread block.

    LOCAL volatile mixed temp[TEMP_SIZE];

    const int indexOffset = NUM_VARIABLES * vectorIndex;

    mixed vectorBeta = 0;
    for (int i = LOCAL_ID; i < NUM_VARIABLES; i += LOCAL_SIZE) {
        vectorBeta += dir[i] * gradDiff[indexOffset + i];
    }
    vectorBeta = reduceAdd(vectorBeta, temp);

    if (LOCAL_ID == 0) {
        returnValue[0] = vectorBeta / scale[vectorIndex];
    }
}

KERNEL void updateDirAlpha(
    GLOBAL mixed* dir,
    GLOBAL const mixed* RESTRICT alpha,
    GLOBAL const mixed* RESTRICT gradDiff,
    const int vectorIndex
) {
    const int indexOffset = NUM_VARIABLES * vectorIndex;
    const mixed vectorScale = -alpha[vectorIndex];

    for (int i = GLOBAL_ID; i < NUM_VARIABLES; i += GLOBAL_SIZE) {
        dir[i] += vectorScale * gradDiff[indexOffset + i];
    }
}

KERNEL void updateDirBeta(
    GLOBAL mixed* dir,
    GLOBAL const mixed* RESTRICT alpha,
    GLOBAL const mixed* RESTRICT xDiff,
    GLOBAL const mixed* RESTRICT returnValue,
    const int vectorIndex
) {
    const int indexOffset = NUM_VARIABLES * vectorIndex;
    const mixed vectorScale = alpha[vectorIndex] - returnValue[0];

    for (int i = GLOBAL_ID; i < NUM_VARIABLES; i += GLOBAL_SIZE) {
        dir[i] += vectorScale * xDiff[indexOffset + i];
    }
}


KERNEL void scaleDir(
    GLOBAL mixed* dir,
    GLOBAL const mixed* RESTRICT scale,
    GLOBAL const mixed* RESTRICT returnValue,
    const int vectorIndex
) {
    const mixed vectorScale = scale[vectorIndex] / returnValue[0];

    for (int i = GLOBAL_ID; i < NUM_VARIABLES; i += GLOBAL_SIZE) {
        dir[i] *= vectorScale;
    }
}

KERNEL void lineSearchStep(
    GLOBAL mixed* RESTRICT x,
    GLOBAL const mixed* RESTRICT xPrev,
    GLOBAL const mixed* RESTRICT dir,
    const mixed step
) {
    for (int i = GLOBAL_ID; i < NUM_VARIABLES; i += GLOBAL_SIZE) {
        x[i] = xPrev[i] + step * dir[i];
    }
}

KERNEL void lineSearchDot(
    GLOBAL const mixed* RESTRICT grad,
    GLOBAL const mixed* RESTRICT dir,
    GLOBAL mixed* RESTRICT returnValue
) {
    // This kernel is expected to be executed in a single thread block.

    LOCAL volatile mixed temp[TEMP_SIZE];

    mixed result = 0;
    for (int i = LOCAL_ID; i < NUM_VARIABLES; i += LOCAL_SIZE) {
        result += grad[i] * dir[i];
    }
    result = reduceAdd(result, temp);

    if (LOCAL_ID == 0) {
        returnValue[0] = result;
    }
}
