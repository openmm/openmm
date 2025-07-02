KERNEL void checkSavedPositions(GLOBAL real4* RESTRICT posq, GLOBAL real4* RESTRICT savedPositions, GLOBAL int* RESTRICT result) {
    for (int i = GLOBAL_ID; i < NUM_PARTICLES; i += GLOBAL_SIZE) {
        real4 posqPosition = posq[i];
        real4 savedPosition = savedPositions[i];
        if (posqPosition.x != savedPosition.x || posqPosition.y != savedPosition.y || posqPosition.z != savedPosition.z) {
            *result = 1;
            break;
        }
    }
}
