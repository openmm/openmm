/**
 * Sum the buffer of energy derivatives with respect to parameters.  The buffer stores the
 * derivatives with respect to all parameters for each thread, and each work group sums the
 * derivatives with respect to one parameter.
 */
KERNEL void reduceEnergyParamDerivs(GLOBAL const mixed* RESTRICT derivBuffer, GLOBAL mixed* RESTRICT result, int bufferSize, int numDerivs) {
    LOCAL mixed tempBuffer[512];
    for (int param = GROUP_ID; param < numDerivs; param += NUM_GROUPS) {
        mixed sum = 0;
        for (int index = param+numDerivs*LOCAL_ID; index < bufferSize; index += numDerivs*LOCAL_SIZE)
            sum += derivBuffer[index];
        SYNC_THREADS;
        tempBuffer[LOCAL_ID] = sum;
        for (int i = 1; i < LOCAL_SIZE; i *= 2) {
            SYNC_THREADS;
            if (LOCAL_ID%(i*2) == 0 && LOCAL_ID+i < LOCAL_SIZE)
                tempBuffer[LOCAL_ID] += tempBuffer[LOCAL_ID+i];
        }
        if (LOCAL_ID == 0)
            result[param] = tempBuffer[0];
    }
}
