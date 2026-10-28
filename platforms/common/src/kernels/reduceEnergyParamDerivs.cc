/**
 * Sum the buffer of energy derivatives with respect to parameters.  The buffer stores the
 * derivatives with respect to all parameters for each thread, and each work group sums the
 * derivatives with respect to one parameter.
 */
KERNEL void reduceEnergyParamDerivs(GLOBAL const mixed* RESTRICT derivBuffer, GLOBAL mixed* RESTRICT result, int bufferSize, int numDerivs, int workGroupSize) {
    LOCAL mixed tempBuffer[512];
    mixed sum = 0;
    for (int index = GROUP_ID+numDerivs*LOCAL_ID; index < bufferSize; index += numDerivs*workGroupSize)
        sum += derivBuffer[index];
    tempBuffer[LOCAL_ID] = sum;
    for (int i = 1; i < workGroupSize; i *= 2) {
        SYNC_THREADS;
        if (LOCAL_ID%(i*2) == 0 && LOCAL_ID+i < workGroupSize)
            tempBuffer[LOCAL_ID] += tempBuffer[LOCAL_ID+i];
    }
    if (LOCAL_ID == 0)
        result[GROUP_ID] = tempBuffer[0];
}
