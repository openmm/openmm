KERNEL void minimizeDrudePositions(int numDrude, int paddedNumAtoms, float tolerance, GLOBAL real4* RESTRICT posq,
        GLOBAL const mm_long* RESTRICT force, GLOBAL float4* RESTRICT drudeParams, GLOBAL int* RESTRICT drudeIndex,
        GLOBAL int4* RESTRICT drudeParents) {
    const real scale = 1/(real) 0x100000000;
    for (int i = GLOBAL_ID; i < numDrude; i += GLOBAL_SIZE) {
        int index = drudeIndex[i];
        int4 parents = drudeParents[i];
        float4 params = drudeParams[i];
        real3 fscale = make_real3(params.z, params.z, params.z);
        if (parents.y != -1) {
            real3 dir = trimTo3(posq[parents.x]-posq[parents.y]);
            dir *= RSQRT(dot(dir, dir));
            fscale += params.x*dir;
        }
        if (parents.z != -1 && parents.w != -1) {
            real3 dir = trimTo3(posq[parents.z]-posq[parents.w]);
            dir *= RSQRT(dot(dir, dir));
            fscale += params.y*dir;
        }
        real4 pos = posq[index];
        real4 f = make_real4(scale*force[index], scale*force[index+paddedNumAtoms], scale*force[index+paddedNumAtoms*2], 0);
        real damping = (SQRT(f.x*f.x + f.y*f.y + f.z*f.z) > 10*tolerance ? 0.5f : 1.0f);
        pos.x += damping*f.x/fscale.x;
        pos.y += damping*f.y/fscale.y;
        pos.z += damping*f.z/fscale.z;
        posq[index] = pos;
    }
}

