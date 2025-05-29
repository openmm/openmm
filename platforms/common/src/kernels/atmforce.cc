KERNEL void hybridForce(int numParticles,
                        int paddedNumParticles,
                        GLOBAL mm_long* RESTRICT force,
                        GLOBAL mm_long* RESTRICT force0,
                        GLOBAL mm_long* RESTRICT force1,
                        GLOBAL int* RESTRICT invAtomOrder,
                        GLOBAL int* RESTRICT inner0InvAtomOrder,
                        GLOBAL int* RESTRICT inner1InvAtomOrder,
                        real dEdu0,
                        real dEdu1) {
    for (int i = GLOBAL_ID; i < numParticles; i += GLOBAL_SIZE) {
        int index = invAtomOrder[i];
        int index0 = inner0InvAtomOrder[i];
        int index1 = inner1InvAtomOrder[i];
        force[index] += (mm_long) (dEdu0*force0[index0] + dEdu1*force1[index1]);
        force[index+paddedNumParticles] += (mm_long) (dEdu0*force0[index0+paddedNumParticles] + dEdu1*force1[index1+paddedNumParticles]);
        force[index+paddedNumParticles*2] += (mm_long) (dEdu0*force0[index0+paddedNumParticles*2] + dEdu1*force1[index1+paddedNumParticles*2]);
    }
}

KERNEL void copyState(int numParticles,
                      GLOBAL real4* RESTRICT posq,
                      GLOBAL real4* RESTRICT posq0,
                      GLOBAL real4* RESTRICT posq1,
                      GLOBAL float4* RESTRICT displ0,
		      GLOBAL float4* RESTRICT displ1,
                      GLOBAL int* RESTRICT atomOrder,
                      GLOBAL int* RESTRICT inner0InvAtomOrder,
                      GLOBAL int* RESTRICT inner1InvAtomOrder
#ifdef USE_MIXED_PRECISION
                      ,
                      GLOBAL real4* RESTRICT posqCorrection,
                      GLOBAL real4* RESTRICT posq0Correction,
                      GLOBAL real4* RESTRICT posq1Correction
#endif
                    ) {
    for (int i = GLOBAL_ID; i < numParticles; i += GLOBAL_SIZE) {
        int atom = atomOrder[i];
        int index0 = inner0InvAtomOrder[atom];
        int index1 = inner1InvAtomOrder[atom];
        real4 p0 = posq[i] + make_real4((real) displ0[atom].x, (real) displ0[atom].y, (real) displ0[atom].z, 0);
        real4 p1 = posq[i] + make_real4((real) displ1[atom].x, (real) displ1[atom].y, (real) displ1[atom].z, 0);
        p0.w = posq0[i].w;
        p1.w = posq1[i].w;
        posq0[index0] = p0;
        posq1[index1] = p1;
#ifdef USE_MIXED_PRECISION
        posq0Correction[index0] = posqCorrection[i];
        posq1Correction[index1] = posqCorrection[i];
#endif
    }
}


