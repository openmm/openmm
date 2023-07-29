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
                      GLOBAL float4* RESTRICT displ1,
		      GLOBAL float4* RESTRICT displ0
#ifdef USE_MIXED_PRECISION
                      ,
                      GLOBAL real4* RESTRICT posqCorrection,
                      GLOBAL real4* RESTRICT posq1Correction,
                      GLOBAL real4* RESTRICT posq2Correction
#endif
                    ) {
    for (int i = GLOBAL_ID; i < numParticles; i += GLOBAL_SIZE) {
        //displacement for initial state
        real4 d0 = make_real4((real) displ0[i].x, (real) displ0[i].y, (real) displ0[i].z, 0);
	//displacement for target state
        real4 d1 = make_real4((real) displ1[i].x, (real) displ1[i].y, (real) displ1[i].z, 0);
        posq1[i] = posq[i] + d0;
        posq2[i] = posq[i] + d1;
#ifdef USE_MIXED_PRECISION
        posq1Correction[i] = posqCorrection[i];
        posq2Correction[i] = posqCorrection[i];
#endif
    }
}


