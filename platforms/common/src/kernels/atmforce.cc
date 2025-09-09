KERNEL void hybridForce(int numParticles,
                        int paddedNumParticles,
                        GLOBAL mm_long* RESTRICT force,
                        GLOBAL mm_long* RESTRICT force0,
                        GLOBAL mm_long* RESTRICT force1,
                        GLOBAL mm_long* RESTRICT dforce0,
                        GLOBAL mm_long* RESTRICT dforce1,
                        GLOBAL int* RESTRICT invAtomOrder,
                        GLOBAL int* RESTRICT inner0InvAtomOrder,
                        GLOBAL int* RESTRICT inner1InvAtomOrder,
                        real dEdu0,
                        real dEdu1) {
    for (int i = GLOBAL_ID; i < numParticles; i += GLOBAL_SIZE) {
        int index = invAtomOrder[i];
        int index0 = inner0InvAtomOrder[i];
        int index1 = inner1InvAtomOrder[i];
        mm_long fx0 = force0[index0]+dforce0[index];
        mm_long fy0 = force0[index0+paddedNumParticles]+dforce0[index+paddedNumParticles];
        mm_long fz0 = force0[index0+paddedNumParticles*2]+dforce0[index+paddedNumParticles*2];
        mm_long fx1 = force1[index1]+dforce1[index];
        mm_long fy1 = force1[index1+paddedNumParticles]+dforce1[index+paddedNumParticles];
        mm_long fz1 = force1[index1+paddedNumParticles*2]+dforce1[index+paddedNumParticles*2];
        force[index]                      += (mm_long) (dEdu0*fx0 + dEdu1*fx1);
        force[index+paddedNumParticles]   += (mm_long) (dEdu0*fy0 + dEdu1*fy1);
        force[index+paddedNumParticles*2] += (mm_long) (dEdu0*fz0 + dEdu1*fz1);
    }
}

//reset variable displacement forces
KERNEL void resetDisplForce(int numParticles,
                            int paddedNumParticles,
                            GLOBAL mm_ulong* RESTRICT dforce0,
                            GLOBAL mm_ulong* RESTRICT dforce1) {
    mm_ulong zero = 0;
    for (int index = GLOBAL_ID; index < numParticles; index += GLOBAL_SIZE) {
        dforce0[index]                      = zero;
        dforce0[index+paddedNumParticles]   = zero;
        dforce0[index+paddedNumParticles*2] = zero;
        dforce1[index]                      = zero;
        dforce1[index+paddedNumParticles]   = zero;
        dforce1[index+paddedNumParticles*2] = zero;
    }
}

//add forces due to variable displacements
KERNEL void displForce(int numParticles,
                       int paddedNumParticles,
                       GLOBAL mm_long* RESTRICT force0,
                       GLOBAL mm_long* RESTRICT force1,
                       GLOBAL mm_long* RESTRICT dforce0,
                       GLOBAL mm_long* RESTRICT dforce1,
                       GLOBAL int4* displParticles,
                       GLOBAL int* RESTRICT atomOrder,
                       GLOBAL int* RESTRICT invAtomOrder,
                       GLOBAL int* RESTRICT inner0InvAtomOrder,
                       GLOBAL int* RESTRICT inner1InvAtomOrder) {
    GLOBAL mm_ulong* df0 = (GLOBAL mm_ulong*) dforce0;
    GLOBAL mm_ulong* df1 = (GLOBAL mm_ulong*) dforce1;
    for (int index = GLOBAL_ID; index < numParticles; index += GLOBAL_SIZE) {
        int atom = atomOrder[index];
        int pj1 = displParticles[atom].x;
        int pi1 = displParticles[atom].y;
        int pj0 = displParticles[atom].z;
        int pi0 = displParticles[atom].w;
        int index0 = inner0InvAtomOrder[atom];
        int index1 = inner1InvAtomOrder[atom];
        if (pj1 >= 0 && pi1 >= 0) {
            int j1 = invAtomOrder[pj1];
            int i1 = invAtomOrder[pi1];
            ATOMIC_ADD(&df1[j1], (mm_ulong) force1[index1]);
            ATOMIC_ADD(&df1[j1+paddedNumParticles], (mm_ulong) force1[index1+paddedNumParticles]);
            ATOMIC_ADD(&df1[j1+paddedNumParticles*2], (mm_ulong) force1[index1+paddedNumParticles*2]);
            ATOMIC_ADD(&df1[i1], (mm_ulong) -force1[index1]);
            ATOMIC_ADD(&df1[i1+paddedNumParticles], (mm_ulong) -force1[index1+paddedNumParticles]);
            ATOMIC_ADD(&df1[i1+paddedNumParticles*2],(mm_ulong)  -force1[index1+paddedNumParticles*2]);
        }
        if (pj0 >= 0 && pi0 >= 0) {
            int j0 = invAtomOrder[pj0];
            int i0 = invAtomOrder[pi0];
            ATOMIC_ADD(&df0[j0], (mm_ulong) force0[index0]);
            ATOMIC_ADD(&df0[j0+paddedNumParticles], (mm_ulong) force0[index0+paddedNumParticles]);
            ATOMIC_ADD(&df0[j0+paddedNumParticles*2],(mm_ulong)  force0[index0+paddedNumParticles*2]);
            ATOMIC_ADD(&df0[i0], (mm_ulong) -force0[index0]);
            ATOMIC_ADD(&df0[i0+paddedNumParticles], (mm_ulong) -force0[index0+paddedNumParticles]);
            ATOMIC_ADD(&df0[i0+paddedNumParticles*2], (mm_ulong) -force0[index0+paddedNumParticles*2]);
        }
    }
}



KERNEL void copyState(int numParticles,
                      GLOBAL real4* RESTRICT posq,
                      GLOBAL real4* RESTRICT posq0,
                      GLOBAL real4* RESTRICT posq1,
                      GLOBAL real4* RESTRICT displacement0,
                      GLOBAL real4* RESTRICT displacement1,
                      GLOBAL int4* displParticles,
                      GLOBAL int* RESTRICT atomOrder,
                      GLOBAL int* RESTRICT invAtomOrder,
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

        //default fixed lab frame displacement
        real4 displ0 = displacement0[atom];
        real4 displ1 = displacement1[atom];
        //override with variable displacements if set
        int pj1 = displParticles[atom].x;
        int pi1 = displParticles[atom].y;
        int pj0 = displParticles[atom].z;
        int pi0 = displParticles[atom].w;
        if (pj1 >= 0 && pi1 >= 0) {
            // variable system coordinate displacements
            int indexj1 = invAtomOrder[pj1];
            int indexi1 = invAtomOrder[pi1];
            displ1 = make_real4((real) posq[indexj1].x- posq[indexi1].x,
                                (real) posq[indexj1].y- posq[indexi1].y,
                                (real) posq[indexj1].z- posq[indexi1].z, (real) 0);
            if (pj0 >= 0 && pi0 >= 0) {
                int indexj0 = invAtomOrder[pj0];
                int indexi0 = invAtomOrder[pi0];
                displ0 = make_real4((real) posq[indexj0].x - posq[indexi0].x,
                                    (real) posq[indexj0].y - posq[indexi0].y,
                                    (real) posq[indexj0].z - posq[indexi0].z, (real) 0);
            }
            else {
                displ0 = make_real4((real) 0, (real) 0, (real) 0, (real) 0);
            }
        }

        int index0 = inner0InvAtomOrder[atom];
        int index1 = inner1InvAtomOrder[atom];
        real4 p0 = posq[i] + make_real4((real) displ0.x, (real) displ0.y, (real) displ0.z, 0);
        real4 p1 = posq[i] + make_real4((real) displ1.x, (real) displ1.y, (real) displ1.z, 0);
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
