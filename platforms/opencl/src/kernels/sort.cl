#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable

float getValue(TYPE value) {
    return SORT_KEY;
}

/**
 * Calculate the minimum and maximum value in the array to be sorted.  This kernel
 * is executed as a single work group.
 */
__kernel void computeRange(__global const TYPE* restrict data, int length, __global float2* restrict range, __local float* restrict buffer) {
    float minimum = MAXFLOAT;
    float maximum = -MAXFLOAT;

    // Each thread calculates the range of a subset of values.

    for (int index = get_local_id(0); index < length; index += get_local_size(0)) {
        float value = getValue(data[index]);
        minimum = min(minimum, value);
        maximum = max(maximum, value);
    }

    // Now reduce them.

    buffer[get_local_id(0)] = minimum;
    barrier(CLK_LOCAL_MEM_FENCE);
    for (int step = 1; step < get_local_size(0); step *= 2) {
        if (get_local_id(0)+step < get_local_size(0) && get_local_id(0)%(2*step) == 0)
            buffer[get_local_id(0)] = min(buffer[get_local_id(0)], buffer[get_local_id(0)+step]);
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    minimum = buffer[0];
    buffer[get_local_id(0)] = maximum;
    barrier(CLK_LOCAL_MEM_FENCE);
    for (int step = 1; step < get_local_size(0); step *= 2) {
        if (get_local_id(0)+step < get_local_size(0) && get_local_id(0)%(2*step) == 0)
            buffer[get_local_id(0)] = max(buffer[get_local_id(0)], buffer[get_local_id(0)+step]);
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    maximum = buffer[0];
    if (get_local_id(0) == 0)
        range[0] = (float2) (minimum, maximum);
}

/**
 * Assign elements to buckets.
 */
__kernel void assignElementsToBuckets(__global const TYPE* restrict data, int length, int numBuckets, __global const float2* restrict range,
        __global int* bucketOffset, __global int* restrict bucketOfElement, __global int* restrict offsetInBucket) {
#ifdef AMD_ATOMIC_WORK_AROUND
    // Do a byte write to force all memory accesses to interactionCount to use the complete path.
    // This avoids the atomic access from causing all word accesses to other buffers from using the slow complete path.
    // The IF actually causes the write to never be executed, its presence is all that is needed.
    // AMD APP SDK 2.4 has this problem.
    if (get_global_id(0) == get_local_id(0)+1)
        ((__global char*)bucketOffset)[sizeof(int)*numBuckets+1] = 0;
#endif
    float2 dataRange = range[0];
    float minValue = dataRange.x;
    float maxValue = dataRange.y;
    float bucketWidth = (maxValue-minValue)/numBuckets;
    for (int index = get_global_id(0); index < length; index += get_global_size(0)) {
        TYPE element = data[index];
        float value = getValue(element);
        int bucketIndex = min((int) ((value-minValue)/bucketWidth), numBuckets-1);
        offsetInBucket[index] = atom_inc(&bucketOffset[bucketIndex]);
        bucketOfElement[index] = bucketIndex;
    }
}

/**
 * Sum the bucket sizes to compute the start position of each bucket.  This kernel
 * is executed as a single work group.
 */
__kernel void computeBucketPositions(int numBuckets, __global int* restrict bucketOffset, __local int* restrict buffer) {
    int globalOffset = 0;
    for (int startBucket = 0; startBucket < numBuckets; startBucket += get_local_size(0)) {
        // Load the bucket sizes into local memory.

        int globalIndex = startBucket+get_local_id(0);
        buffer[get_local_id(0)] = (globalIndex < numBuckets ? bucketOffset[globalIndex] : 0);
        barrier(CLK_LOCAL_MEM_FENCE);

        // Perform a parallel prefix sum.

        for (int step = 1; step < get_local_size(0); step *= 2) {
            int add = (get_local_id(0) >= step ? buffer[get_local_id(0)-step] : 0);
            barrier(CLK_LOCAL_MEM_FENCE);
            buffer[get_local_id(0)] += add;
            barrier(CLK_LOCAL_MEM_FENCE);
        }

        // Write the results back to global memory.

        if (globalIndex < numBuckets)
            bucketOffset[globalIndex] = buffer[get_local_id(0)]+globalOffset;
        globalOffset += buffer[get_local_size(0)-1];
    }
}

/**
 * Copy the input data into the buckets for sorting.
 */
__kernel void copyDataToBuckets(__global const TYPE* restrict data, __global TYPE* restrict buckets, int length, __global const int* restrict bucketOffset, __global const int* restrict bucketOfElement, __global const int* restrict offsetInBucket) {
    for (int index = get_global_id(0); index < length; index += get_global_size(0)) {
        TYPE element = data[index];
        int bucketIndex = bucketOfElement[index];
        int offset = (bucketIndex == 0 ? 0 : bucketOffset[bucketIndex-1]);
        buckets[offset+offsetInBucket[index]] = element;
    }
}

/**
 * Sort the data in each bucket.
 */
__kernel void sortBuckets(__global TYPE* restrict data, __global const TYPE* restrict buckets, int numBuckets, __global const int* restrict bucketOffset, __local TYPE* restrict buffer) {
    for (int index = get_group_id(0); index < numBuckets; index += get_num_groups(0)) {
        int startIndex = (index == 0 ? 0 : bucketOffset[index-1]);
        int endIndex = bucketOffset[index];
        int length = endIndex-startIndex;
        if (length <= get_local_size(0)) {
            // Load the data into local memory.

            buffer[get_local_id(0)] = (get_local_id(0) < length ? buckets[startIndex+get_local_id(0)] : (TYPE) MAXFLOAT);
            barrier(CLK_LOCAL_MEM_FENCE);

            // Perform a bitonic sort in local memory.

            for (int k = 2; k <= get_local_size(0); k *= 2) {
                for (int j = k/2; j > 0; j /= 2) {
                    int ixj = get_local_id(0)^j;
                    if (ixj > get_local_id(0)) {
                        if (((get_local_id(0)&k) == 0 && getValue(buffer[get_local_id(0)]) > getValue(buffer[ixj])) ||
                            ((get_local_id(0)&k) != 0 && getValue(buffer[get_local_id(0)]) < getValue(buffer[ixj]))) {
                            TYPE temp = buffer[get_local_id(0)];
                            buffer[get_local_id(0)] = buffer[ixj];
                            buffer[ixj] = temp;
                        }
                    }
                    barrier(CLK_LOCAL_MEM_FENCE);
                }
            }

            // Write the data to the sorted array.

            if (get_local_id(0) < length)
                data[startIndex+get_local_id(0)] = buffer[get_local_id(0)];
        }
        else {
            // Copy the bucket data over to the output array.

            for (int i = get_local_id(0); i < length; i += get_local_size(0))
                data[startIndex+i] = buckets[startIndex+i];
            barrier(CLK_GLOBAL_MEM_FENCE);

            // Perform a bitonic sort in global memory.

            for (int k = 2; k < 2*length; k *= 2) {
                for (int j = k/2; j > 0; j /= 2) {
                    for (int i = get_local_id(0); i < length; i += get_local_size(0)) {
                        int ixj = i^j;
                        if (ixj > i && ixj < length) {
                            TYPE value1 = data[startIndex+i];
                            TYPE value2 = data[startIndex+ixj];
                            bool ascending = ((i&k) == 0);
                            for (int mask = k*2; mask < 2*length; mask *= 2)
                                ascending = ((i&mask) == 0 ? !ascending : ascending);
                            if ((ascending && getValue(value1) > getValue(value2)) ||
                                (!ascending && getValue(value1) < getValue(value2))) {
                                data[startIndex+i] = value2;
                                data[startIndex+ixj] = value1;
                            }
                        }
                    }
                    barrier(CLK_GLOBAL_MEM_FENCE);
                }
            }
        }
    }
}
