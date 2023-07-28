KERNEL void hybridForce(int numParticles,
                        int paddedNumParticles,
                        GLOBAL mm_long* RESTRICT force,
                        GLOBAL mm_long* RESTRICT force1,
                        GLOBAL mm_long* RESTRICT force2,
                        real dEdu0,
                        real dEdu1) {
    for (int i = GLOBAL_ID; i < numParticles; i += GLOBAL_SIZE) {
        force[i] += (mm_long) (dEdu0*force1[i] + dEdu1*force2[i]);
        force[i+paddedNumParticles] += (mm_long) (dEdu0*force1[i+paddedNumParticles] + dEdu1*force2[i+paddedNumParticles]);
        force[i+paddedNumParticles*2] += (mm_long) (dEdu0*force1[i+paddedNumParticles*2] + dEdu1*force2[i+paddedNumParticles*2]);
    }
}

KERNEL void copyState(int numParticles,
                      GLOBAL real4* RESTRICT posq,
                      GLOBAL real4* RESTRICT posq1,
                      GLOBAL real4* RESTRICT posq2,
                      GLOBAL float4* RESTRICT displ
#ifdef USE_MIXED_PRECISION
                      ,
                      GLOBAL real4* RESTRICT posqCorrection,
                      GLOBAL real4* RESTRICT posq1Correction,
                      GLOBAL real4* RESTRICT posq2Correction
#endif
                    ) {
    for (int i = GLOBAL_ID; i < numParticles; i += GLOBAL_SIZE) {
        real4 d = make_real4((real) displ[i].x, (real) displ[i].y, (real) displ[i].z, 0);
        posq1[i] = posq[i];
        posq2[i] = posq[i] + d;
#ifdef USE_MIXED_PRECISION
        posq1Correction[i] = posqCorrection[i];
        posq2Correction[i] = posqCorrection[i];
#endif
    }
}


