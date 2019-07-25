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
        __global const int4* restrict localCoordsAtoms, __global const real* restrict localCoordsParams) {
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
        int4 atoms = localCoordsAtoms[index];
        __global const real* params = &localCoordsParams[12*index];
        mixed4 pos = loadPos(posq, posqCorrection, atoms.x);
        mixed4 pos1_4 = loadPos(posq, posqCorrection, atoms.y);
        mixed4 pos2_4 = loadPos(posq, posqCorrection, atoms.z);
        mixed4 pos3_4 = loadPos(posq, posqCorrection, atoms.w);
        mixed4 pos1 = (mixed4) (pos1_4.x, pos1_4.y, pos1_4.z, 0);
        mixed4 pos2 = (mixed4) (pos2_4.x, pos2_4.y, pos2_4.z, 0);
        mixed4 pos3 = (mixed4) (pos3_4.x, pos3_4.y, pos3_4.z, 0);
        mixed4 originWeights = (mixed4) (params[0], params[1], params[2], 0);
        mixed4 xWeights = (mixed4) (params[3], params[4], params[5], 0);
        mixed4 yWeights = (mixed4) (params[6], params[7], params[8], 0);
        mixed4 localPosition = (mixed4) (params[9], params[10], params[11], 0);
        mixed4 origin = pos1*originWeights.x + pos2*originWeights.y + pos3*originWeights.z;
        mixed4 xdir = pos1*xWeights.x + pos2*xWeights.y + pos3*xWeights.z;
        mixed4 ydir = pos1*yWeights.x + pos2*yWeights.y + pos3*yWeights.z;
        mixed4 zdir = cross(xdir, ydir);
        xdir *= rsqrt(xdir.x*xdir.x+xdir.y*xdir.y+xdir.z*xdir.z);
        zdir *= rsqrt(zdir.x*zdir.x+zdir.y*zdir.y+zdir.z*zdir.z);
        ydir = cross(zdir, xdir);
        pos.x = origin.x + xdir.x*localPosition.x + ydir.x*localPosition.y + zdir.x*localPosition.z;
        pos.y = origin.y + xdir.y*localPosition.x + ydir.y*localPosition.y + zdir.y*localPosition.z;
        pos.z = origin.z + xdir.z*localPosition.x + ydir.z*localPosition.y + zdir.z*localPosition.z;
        storePos(posq, posqCorrection, atoms.x, pos);
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
        __global const int4* restrict localCoordsAtoms, __global const real* restrict localCoordsParams) {
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
        int4 atoms = localCoordsAtoms[index];
        __global const real* params = &localCoordsParams[12*index];
        mixed4 pos = loadPos(posq, posqCorrection, atoms.x);
        mixed4 pos1_4 = loadPos(posq, posqCorrection, atoms.y);
        mixed4 pos2_4 = loadPos(posq, posqCorrection, atoms.z);
        mixed4 pos3_4 = loadPos(posq, posqCorrection, atoms.w);
        mixed4 pos1 = (mixed4) (pos1_4.x, pos1_4.y, pos1_4.z, 0);
        mixed4 pos2 = (mixed4) (pos2_4.x, pos2_4.y, pos2_4.z, 0);
        mixed4 pos3 = (mixed4) (pos3_4.x, pos3_4.y, pos3_4.z, 0);
        mixed4 originWeights = (mixed4) (params[0], params[1], params[2], 0);
        mixed4 wx = (mixed4) (params[3], params[4], params[5], 0);
        mixed4 wy = (mixed4) (params[6], params[7], params[8], 0);
        mixed4 localPosition = (mixed4) (params[9], params[10], params[11], 0);
        mixed4 origin = pos1*originWeights.x + pos2*originWeights.y + pos3*originWeights.z;
        mixed4 xdir = pos1*wx.x + pos2*wx.y + pos3*wx.z;
        mixed4 ydir = pos1*wy.x + pos2*wy.y + pos3*wy.z;
        mixed4 zdir = cross(xdir, ydir);
        mixed invNormXdir = rsqrt(xdir.x*xdir.x+xdir.y*xdir.y+xdir.z*xdir.z);
        mixed invNormZdir = rsqrt(zdir.x*zdir.x+zdir.y*zdir.y+zdir.z*zdir.z);
        mixed4 dx = xdir*invNormXdir;
        mixed4 dz = zdir*invNormZdir;
        mixed4 dy = cross(dz, dx);

        // The derivatives for this case are very complicated.  They were computed with SymPy then simplified by hand.

        mixed t11 = (wx.x*ydir.x-wy.x*xdir.x)*invNormZdir;
        mixed t12 = (wx.x*ydir.y-wy.x*xdir.y)*invNormZdir;
        mixed t13 = (wx.x*ydir.z-wy.x*xdir.z)*invNormZdir;
        mixed t21 = (wx.y*ydir.x-wy.y*xdir.x)*invNormZdir;
        mixed t22 = (wx.y*ydir.y-wy.y*xdir.y)*invNormZdir;
        mixed t23 = (wx.y*ydir.z-wy.y*xdir.z)*invNormZdir;
        mixed t31 = (wx.z*ydir.x-wy.z*xdir.x)*invNormZdir;
        mixed t32 = (wx.z*ydir.y-wy.z*xdir.y)*invNormZdir;
        mixed t33 = (wx.z*ydir.z-wy.z*xdir.z)*invNormZdir;
        mixed sx1 = t13*dz.y-t12*dz.z;
        mixed sy1 = t11*dz.z-t13*dz.x;
        mixed sz1 = t12*dz.x-t11*dz.y;
        mixed sx2 = t23*dz.y-t22*dz.z;
        mixed sy2 = t21*dz.z-t23*dz.x;
        mixed sz2 = t22*dz.x-t21*dz.y;
        mixed sx3 = t33*dz.y-t32*dz.z;
        mixed sy3 = t31*dz.z-t33*dz.x;
        mixed sz3 = t32*dz.x-t31*dz.y;
        mixed4 wxScaled = wx*invNormXdir;
        real4 f = force[atoms.x];
        real4 f1 = 0;
        real4 f2 = 0;
        real4 f3 = 0;
        mixed4 fp1 = localPosition*f.x;
        mixed4 fp2 = localPosition*f.y;
        mixed4 fp3 = localPosition*f.z;
        f1.x += fp1.x*wxScaled.x*(1-dx.x*dx.x) + fp1.z*(dz.x*sx1    ) + fp1.y*((-dx.x*dy.x     )*wxScaled.x + dy.x*sx1 - dx.y*t12 - dx.z*t13) + f.x*originWeights.x;
        f1.y += fp1.x*wxScaled.x*( -dx.x*dx.y) + fp1.z*(dz.x*sy1+t13) + fp1.y*((-dx.y*dy.x-dz.z)*wxScaled.x + dy.x*sy1 + dx.y*t11);
        f1.z += fp1.x*wxScaled.x*( -dx.x*dx.z) + fp1.z*(dz.x*sz1-t12) + fp1.y*((-dx.z*dy.x+dz.y)*wxScaled.x + dy.x*sz1 + dx.z*t11);
        f2.x += fp1.x*wxScaled.y*(1-dx.x*dx.x) + fp1.z*(dz.x*sx2    ) + fp1.y*((-dx.x*dy.x     )*wxScaled.y + dy.x*sx2 - dx.y*t22 - dx.z*t23) + f.x*originWeights.y;
        f2.y += fp1.x*wxScaled.y*( -dx.x*dx.y) + fp1.z*(dz.x*sy2+t23) + fp1.y*((-dx.y*dy.x-dz.z)*wxScaled.y + dy.x*sy2 + dx.y*t21);
        f2.z += fp1.x*wxScaled.y*( -dx.x*dx.z) + fp1.z*(dz.x*sz2-t22) + fp1.y*((-dx.z*dy.x+dz.y)*wxScaled.y + dy.x*sz2 + dx.z*t21);
        f3.x += fp1.x*wxScaled.z*(1-dx.x*dx.x) + fp1.z*(dz.x*sx3    ) + fp1.y*((-dx.x*dy.x     )*wxScaled.z + dy.x*sx3 - dx.y*t32 - dx.z*t33) + f.x*originWeights.z;
        f3.y += fp1.x*wxScaled.z*( -dx.x*dx.y) + fp1.z*(dz.x*sy3+t33) + fp1.y*((-dx.y*dy.x-dz.z)*wxScaled.z + dy.x*sy3 + dx.y*t31);
        f3.z += fp1.x*wxScaled.z*( -dx.x*dx.z) + fp1.z*(dz.x*sz3-t32) + fp1.y*((-dx.z*dy.x+dz.y)*wxScaled.z + dy.x*sz3 + dx.z*t31);
        f1.x += fp2.x*wxScaled.x*( -dx.y*dx.x) + fp2.z*(dz.y*sx1-t13) - fp2.y*(( dx.x*dy.y-dz.z)*wxScaled.x - dy.y*sx1 - dx.x*t12);
        f1.y += fp2.x*wxScaled.x*(1-dx.y*dx.y) + fp2.z*(dz.y*sy1    ) - fp2.y*(( dx.y*dy.y     )*wxScaled.x - dy.y*sy1 + dx.x*t11 + dx.z*t13) + f.y*originWeights.x;
        f1.z += fp2.x*wxScaled.x*( -dx.y*dx.z) + fp2.z*(dz.y*sz1+t11) - fp2.y*(( dx.z*dy.y+dz.x)*wxScaled.x - dy.y*sz1 - dx.z*t12);
        f2.x += fp2.x*wxScaled.y*( -dx.y*dx.x) + fp2.z*(dz.y*sx2-t23) - fp2.y*(( dx.x*dy.y-dz.z)*wxScaled.y - dy.y*sx2 - dx.x*t22);
        f2.y += fp2.x*wxScaled.y*(1-dx.y*dx.y) + fp2.z*(dz.y*sy2    ) - fp2.y*(( dx.y*dy.y     )*wxScaled.y - dy.y*sy2 + dx.x*t21 + dx.z*t23) + f.y*originWeights.y;
        f2.z += fp2.x*wxScaled.y*( -dx.y*dx.z) + fp2.z*(dz.y*sz2+t21) - fp2.y*(( dx.z*dy.y+dz.x)*wxScaled.y - dy.y*sz2 - dx.z*t22);
        f3.x += fp2.x*wxScaled.z*( -dx.y*dx.x) + fp2.z*(dz.y*sx3-t33) - fp2.y*(( dx.x*dy.y-dz.z)*wxScaled.z - dy.y*sx3 - dx.x*t32);
        f3.y += fp2.x*wxScaled.z*(1-dx.y*dx.y) + fp2.z*(dz.y*sy3    ) - fp2.y*(( dx.y*dy.y     )*wxScaled.z - dy.y*sy3 + dx.x*t31 + dx.z*t33) + f.y*originWeights.z;
        f3.z += fp2.x*wxScaled.z*( -dx.y*dx.z) + fp2.z*(dz.y*sz3+t31) - fp2.y*(( dx.z*dy.y+dz.x)*wxScaled.z - dy.y*sz3 - dx.z*t32);
        f1.x += fp3.x*wxScaled.x*( -dx.z*dx.x) + fp3.z*(dz.z*sx1+t12) + fp3.y*((-dx.x*dy.z-dz.y)*wxScaled.x + dy.z*sx1 + dx.x*t13);
        f1.y += fp3.x*wxScaled.x*( -dx.z*dx.y) + fp3.z*(dz.z*sy1-t11) + fp3.y*((-dx.y*dy.z+dz.x)*wxScaled.x + dy.z*sy1 + dx.y*t13);
        f1.z += fp3.x*wxScaled.x*(1-dx.z*dx.z) + fp3.z*(dz.z*sz1    ) + fp3.y*((-dx.z*dy.z     )*wxScaled.x + dy.z*sz1 - dx.x*t11 - dx.y*t12) + f.z*originWeights.x;
        f2.x += fp3.x*wxScaled.y*( -dx.z*dx.x) + fp3.z*(dz.z*sx2+t22) + fp3.y*((-dx.x*dy.z-dz.y)*wxScaled.y + dy.z*sx2 + dx.x*t23);
        f2.y += fp3.x*wxScaled.y*( -dx.z*dx.y) + fp3.z*(dz.z*sy2-t21) + fp3.y*((-dx.y*dy.z+dz.x)*wxScaled.y + dy.z*sy2 + dx.y*t23);
        f2.z += fp3.x*wxScaled.y*(1-dx.z*dx.z) + fp3.z*(dz.z*sz2    ) + fp3.y*((-dx.z*dy.z     )*wxScaled.y + dy.z*sz2 - dx.x*t21 - dx.y*t22) + f.z*originWeights.y;
        f3.x += fp3.x*wxScaled.z*( -dx.z*dx.x) + fp3.z*(dz.z*sx3+t32) + fp3.y*((-dx.x*dy.z-dz.y)*wxScaled.z + dy.z*sx3 + dx.x*t33);
        f3.y += fp3.x*wxScaled.z*( -dx.z*dx.y) + fp3.z*(dz.z*sy3-t31) + fp3.y*((-dx.y*dy.z+dz.x)*wxScaled.z + dy.z*sy3 + dx.y*t33);
        f3.z += fp3.x*wxScaled.z*(1-dx.z*dx.z) + fp3.z*(dz.z*sz3    ) + fp3.y*((-dx.z*dy.z     )*wxScaled.z + dy.z*sz3 - dx.x*t31 - dx.y*t32) + f.z*originWeights.z;
        ADD_FORCE(atoms.y, f1);
        ADD_FORCE(atoms.z, f2);
        ADD_FORCE(atoms.w, f3);
    }
}
