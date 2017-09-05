/**
 * Load the position of a particle.
 */
mixed4 loadPos(__global const real4* restrict posq, __global const real4* restrict posqCorrection, int index) {
#ifdef USE_MIXED_PRECISION
    real4 pos1 = posq[index];
    real4 pos2 = posqCorrection[index];
    return (mixed4) (pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
#else
    return posq[index];
#endif
}

/**
 * Store the position of a particle.
 */
void storePos(__global real4* restrict posq, __global real4* restrict posqCorrection, int index, mixed4 pos) {
#ifdef USE_MIXED_PRECISION
    posq[index] = (real4) ((real) pos.x, (real) pos.y, (real) pos.z, (real) pos.w);
    posqCorrection[index] = (real4) (pos.x-(real) pos.x, pos.y-(real) pos.y, pos.z-(real) pos.z, 0);
#else
    posq[index] = pos;
#endif
}

/**
 * Compute the positions of virtual sites
 */
__kernel void computeVirtualSites(__global real4* restrict posq,
#ifdef USE_MIXED_PRECISION
        __global real4* restrict posqCorrection,
#endif
        __global const int4* restrict avg2Atoms, __global const real2* restrict avg2Weights,
        __global const int4* restrict avg3Atoms, __global const real4* restrict avg3Weights,
        __global const int4* restrict outOfPlaneAtoms, __global const real4* restrict outOfPlaneWeights,
        __global const int* restrict localCoordsIndex, __global const int* restrict localCoordsAtoms,
        __global const real* restrict localCoordsWeights, __global const real4* restrict localCoordsPos,
        __global const int* restrict localCoordsStartIndex) {
#ifndef USE_MIXED_PRECISION
        __global real4* posqCorrection = 0;
#endif
    
    // Two particle average sites.
    
    for (int index = get_global_id(0); index < NUM_2_AVERAGE; index += get_global_size(0)) {
        int4 atoms = avg2Atoms[index];
        real2 weights = avg2Weights[index];
        mixed4 pos = loadPos(posq, posqCorrection, atoms.x);
        mixed4 pos1 = loadPos(posq, posqCorrection, atoms.y);
        mixed4 pos2 = loadPos(posq, posqCorrection, atoms.z);
        pos.xyz = pos1.xyz*weights.x + pos2.xyz*weights.y;
        storePos(posq, posqCorrection, atoms.x, pos);
    }
    
    // Three particle average sites.
    
    for (int index = get_global_id(0); index < NUM_3_AVERAGE; index += get_global_size(0)) {
        int4 atoms = avg3Atoms[index];
        real4 weights = avg3Weights[index];
        mixed4 pos = loadPos(posq, posqCorrection, atoms.x);
        mixed4 pos1 = loadPos(posq, posqCorrection, atoms.y);
        mixed4 pos2 = loadPos(posq, posqCorrection, atoms.z);
        mixed4 pos3 = loadPos(posq, posqCorrection, atoms.w);
        pos.xyz = pos1.xyz*weights.x + pos2.xyz*weights.y + pos3.xyz*weights.z;
        storePos(posq, posqCorrection, atoms.x, pos);
    }
    
    // Out of plane sites.
    
    for (int index = get_global_id(0); index < NUM_OUT_OF_PLANE; index += get_global_size(0)) {
        int4 atoms = outOfPlaneAtoms[index];
        real4 weights = outOfPlaneWeights[index];
        mixed4 pos = loadPos(posq, posqCorrection, atoms.x);
        mixed4 pos1 = loadPos(posq, posqCorrection, atoms.y);
        mixed4 pos2 = loadPos(posq, posqCorrection, atoms.z);
        mixed4 pos3 = loadPos(posq, posqCorrection, atoms.w);
        mixed4 v12 = pos2-pos1;
        mixed4 v13 = pos3-pos1;
        pos.xyz = pos1.xyz + v12.xyz*weights.x + v13.xyz*weights.y + cross(v12, v13).xyz*weights.z;
        storePos(posq, posqCorrection, atoms.x, pos);
    }
    
    // Local coordinates sites.
    
    for (int index = get_global_id(0); index < NUM_LOCAL_COORDS; index += get_global_size(0)) {
        int siteAtomIndex = localCoordsIndex[index];
        int start = localCoordsStartIndex[index];
        int end = localCoordsStartIndex[index+1];
        mixed3 origin = 0, xdir = 0, ydir = 0;
        for (int j = start; j < end; j++) {
            mixed3 pos = loadPos(posq, posqCorrection, localCoordsAtoms[j]).xyz;
            origin += pos*localCoordsWeights[3*j];
            xdir += pos*localCoordsWeights[3*j+1];
            ydir += pos*localCoordsWeights[3*j+2];
        }
        mixed3 zdir = cross(xdir, ydir);
        mixed normXdir = sqrt(xdir.x*xdir.x+xdir.y*xdir.y+xdir.z*xdir.z);
        mixed normZdir = sqrt(zdir.x*zdir.x+zdir.y*zdir.y+zdir.z*zdir.z);
        mixed invNormXdir = (normXdir > 0 ? 1/normXdir : 0);
        mixed invNormZdir = (normZdir > 0 ? 1/normZdir : 0);
        xdir *= invNormXdir;
        zdir *= invNormZdir;
        ydir = cross(zdir, xdir);
        mixed3 localPosition = convert_mixed4(localCoordsPos[index]).xyz;
        mixed4 pos = loadPos(posq, posqCorrection, siteAtomIndex);
        pos.x = origin.x + xdir.x*localPosition.x + ydir.x*localPosition.y + zdir.x*localPosition.z;
        pos.y = origin.y + xdir.y*localPosition.x + ydir.y*localPosition.y + zdir.y*localPosition.z;
        pos.z = origin.z + xdir.z*localPosition.x + ydir.z*localPosition.y + zdir.z*localPosition.z;
        storePos(posq, posqCorrection, siteAtomIndex, pos);
    }
}

#ifdef HAS_OVERLAPPING_VSITES
    #ifdef SUPPORTS_64_BIT_ATOMICS
        // We will use 64 bit atomics for force redistribution.

        #define ADD_FORCE(index, f) addForce(index, f, longForce);

        #pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable

        void addForce(int index, real4 f,  __global long* longForce) {
            atom_add(&longForce[index], (long) (f.x*0x100000000));
            atom_add(&longForce[index+PADDED_NUM_ATOMS], (long) (f.y*0x100000000));
            atom_add(&longForce[index+2*PADDED_NUM_ATOMS], (long) (f.z*0x100000000));
        }

        __kernel void addDistributedForces(__global const long* restrict longForces, __global real4* restrict forces) {
            real scale = 1/(real) 0x100000000;
            for (int index = get_global_id(0); index < NUM_ATOMS; index += get_global_size(0)) {
                real4 f = (real4) (scale*longForces[index], scale*longForces[index+PADDED_NUM_ATOMS], scale*longForces[index+2*PADDED_NUM_ATOMS], 0);
                forces[index] += f;
            }
        }
    #else
        // 64 bit atomics aren't supported, so we have to use atomic_cmpxchg() for force redistribution.
        
        #define ADD_FORCE(index, f) addForce(index, f, force);

        void atomicAddFloat(__global float* p, float v) {
            __global int* ip = (__global int*) p;
            while (true) {
                union {
                    float f;
                    int i;
                } oldval, newval;
                oldval.f = *p;
                newval.f = oldval.f+v;
                if (atomic_cmpxchg(ip, oldval.i, newval.i) == oldval.i)
                    return;
            }
        }

        void addForce(int index, float4 f, __global float4* force) {
            __global float* components = (__global float*) force;
            atomicAddFloat(&components[4*index], f.x);
            atomicAddFloat(&components[4*index+1], f.y);
            atomicAddFloat(&components[4*index+2], f.z);
        }
    #endif
#else
    // There are no overlapping virtual sites, so we can just store forces directly.

    #define ADD_FORCE(index, f) force[index].xyz += (f).xyz;
#endif

/**
 * Distribute forces from virtual sites to the atoms they are based on.
 */
__kernel void distributeForces(__global const real4* restrict posq, __global real4* restrict force,
#ifdef SUPPORTS_64_BIT_ATOMICS
        __global long* restrict longForce,
#endif
#ifdef USE_MIXED_PRECISION
        __global real4* restrict posqCorrection,
#endif
        __global const int4* restrict avg2Atoms, __global const real2* restrict avg2Weights,
        __global const int4* restrict avg3Atoms, __global const real4* restrict avg3Weights,
        __global const int4* restrict outOfPlaneAtoms, __global const real4* restrict outOfPlaneWeights,
        __global const int* restrict localCoordsIndex, __global const int* restrict localCoordsAtoms,
        __global const real* restrict localCoordsWeights, __global const real4* restrict localCoordsPos,
        __global const int* restrict localCoordsStartIndex) {
#ifndef USE_MIXED_PRECISION
        __global real4* posqCorrection = 0;
#endif

    // Two particle average sites.
    
    for (int index = get_global_id(0); index < NUM_2_AVERAGE; index += get_global_size(0)) {
        int4 atoms = avg2Atoms[index];
        real2 weights = avg2Weights[index];
        real4 f = force[atoms.x];
        ADD_FORCE(atoms.y, f*weights.x);
        ADD_FORCE(atoms.z, f*weights.y);
    }
    
    // Three particle average sites.
    
    for (int index = get_global_id(0); index < NUM_3_AVERAGE; index += get_global_size(0)) {
        int4 atoms = avg3Atoms[index];
        real4 weights = avg3Weights[index];
        real4 f = force[atoms.x];
        ADD_FORCE(atoms.y, f*weights.x);
        ADD_FORCE(atoms.z, f*weights.y);
        ADD_FORCE(atoms.w, f*weights.z);
    }
    
    // Out of plane sites.
    
    for (int index = get_global_id(0); index < NUM_OUT_OF_PLANE; index += get_global_size(0)) {
        int4 atoms = outOfPlaneAtoms[index];
        real4 weights = outOfPlaneWeights[index];
        mixed4 pos1 = loadPos(posq, posqCorrection, atoms.y);
        mixed4 pos2 = loadPos(posq, posqCorrection, atoms.z);
        mixed4 pos3 = loadPos(posq, posqCorrection, atoms.w);
        mixed4 v12 = pos2-pos1;
        mixed4 v13 = pos3-pos1;
        real4 f = force[atoms.x];
        real4 fp2 = (real4) (weights.x*f.x - weights.z*v13.z*f.y + weights.z*v13.y*f.z,
                   weights.z*v13.z*f.x + weights.x*f.y - weights.z*v13.x*f.z,
                  -weights.z*v13.y*f.x + weights.z*v13.x*f.y + weights.x*f.z, 0.0f);
        real4 fp3 = (real4) (weights.y*f.x + weights.z*v12.z*f.y - weights.z*v12.y*f.z,
                  -weights.z*v12.z*f.x + weights.y*f.y + weights.z*v12.x*f.z,
                   weights.z*v12.y*f.x - weights.z*v12.x*f.y + weights.y*f.z, 0.0f);
        ADD_FORCE(atoms.y, f-fp2-fp3);
        ADD_FORCE(atoms.z, fp2);
        ADD_FORCE(atoms.w, fp3);
    }
    
    // Local coordinates sites.
    
    for (int index = get_global_id(0); index < NUM_LOCAL_COORDS; index += get_global_size(0)) {
        int siteAtomIndex = localCoordsIndex[index];
        int start = localCoordsStartIndex[index];
        int end = localCoordsStartIndex[index+1];
        mixed3 origin = 0, xdir = 0, ydir = 0;
        for (int j = start; j < end; j++) {
            mixed3 pos = loadPos(posq, posqCorrection, localCoordsAtoms[j]).xyz;
            origin += pos*localCoordsWeights[3*j];
            xdir += pos*localCoordsWeights[3*j+1];
            ydir += pos*localCoordsWeights[3*j+2];
        }
        mixed3 zdir = cross(xdir, ydir);
        mixed normXdir = sqrt(xdir.x*xdir.x+xdir.y*xdir.y+xdir.z*xdir.z);
        mixed normZdir = sqrt(zdir.x*zdir.x+zdir.y*zdir.y+zdir.z*zdir.z);
        mixed invNormXdir = (normXdir > 0 ? 1/normXdir : 0);
        mixed invNormZdir = (normZdir > 0 ? 1/normZdir : 0);
        mixed3 dx = xdir*invNormXdir;
        mixed3 dz = zdir*invNormZdir;
        mixed3 dy = cross(dz, dx);
        mixed3 localPosition = convert_mixed4(localCoordsPos[index]).xyz;

        // The derivatives for this case are very complicated.  They were computed with SymPy then simplified by hand.

        real4 f = force[siteAtomIndex];
        mixed3 fp1 = localPosition*f.x;
        mixed3 fp2 = localPosition*f.y;
        mixed3 fp3 = localPosition*f.z;
        for (int j = start; j < end; j++) {
            real originWeight = localCoordsWeights[3*j];
            real wx = localCoordsWeights[3*j+1];
            real wy = localCoordsWeights[3*j+2];
            mixed wxScaled = wx*invNormXdir;
            mixed t1 = (wx*ydir.x-wy*xdir.x)*invNormZdir;
            mixed t2 = (wx*ydir.y-wy*xdir.y)*invNormZdir;
            mixed t3 = (wx*ydir.z-wy*xdir.z)*invNormZdir;
            mixed sx = t3*dz.y-t2*dz.z;
            mixed sy = t1*dz.z-t3*dz.x;
            mixed sz = t2*dz.x-t1*dz.y;
            real4 fresult = 0;
            fresult.x += fp1.x*wxScaled*(1-dx.x*dx.x) + fp1.z*(dz.x*sx   ) + fp1.y*((-dx.x*dy.x     )*wxScaled + dy.x*sx - dx.y*t2 - dx.z*t3) + f.x*originWeight;
            fresult.y += fp1.x*wxScaled*( -dx.x*dx.y) + fp1.z*(dz.x*sy+t3) + fp1.y*((-dx.y*dy.x-dz.z)*wxScaled + dy.x*sy + dx.y*t1);
            fresult.z += fp1.x*wxScaled*( -dx.x*dx.z) + fp1.z*(dz.x*sz-t2) + fp1.y*((-dx.z*dy.x+dz.y)*wxScaled + dy.x*sz + dx.z*t1);
            fresult.x += fp2.x*wxScaled*( -dx.y*dx.x) + fp2.z*(dz.y*sx-t3) - fp2.y*(( dx.x*dy.y-dz.z)*wxScaled - dy.y*sx - dx.x*t2);
            fresult.y += fp2.x*wxScaled*(1-dx.y*dx.y) + fp2.z*(dz.y*sy   ) - fp2.y*(( dx.y*dy.y     )*wxScaled - dy.y*sy + dx.x*t1 + dx.z*t3) + f.y*originWeight;
            fresult.z += fp2.x*wxScaled*( -dx.y*dx.z) + fp2.z*(dz.y*sz+t1) - fp2.y*(( dx.z*dy.y+dz.x)*wxScaled - dy.y*sz - dx.z*t2);
            fresult.x += fp3.x*wxScaled*( -dx.z*dx.x) + fp3.z*(dz.z*sx+t2) + fp3.y*((-dx.x*dy.z-dz.y)*wxScaled + dy.z*sx + dx.x*t3);
            fresult.y += fp3.x*wxScaled*( -dx.z*dx.y) + fp3.z*(dz.z*sy-t1) + fp3.y*((-dx.y*dy.z+dz.x)*wxScaled + dy.z*sy + dx.y*t3);
            fresult.z += fp3.x*wxScaled*(1-dx.z*dx.z) + fp3.z*(dz.z*sz   ) + fp3.y*((-dx.z*dy.z     )*wxScaled + dy.z*sz - dx.x*t1 - dx.y*t2) + f.z*originWeight;
            ADD_FORCE(localCoordsAtoms[j], fresult);
        }
    }
}
