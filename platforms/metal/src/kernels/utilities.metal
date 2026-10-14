/**
 * This is called to determine the accuracy of the fast versions of various functions.
 */
KERNEL void determineNativeAccuracy(GLOBAL float* RESTRICT values, int numValues) {
    for (int i = GLOBAL_ID; i < numValues; i += GLOBAL_SIZE) {
        float v = values[6*i];
        values[6*i+1] = fast::sqrt(v);
        values[6*i+2] = fast::rsqrt(v);
        values[6*i+3] = fast::divide(1.0, v);
        values[6*i+4] = fast::exp(v);
        values[6*i+5] = fast::log(v);
    }
}
