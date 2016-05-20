#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable

#define TILE_SIZE 32

/**
 * Calculate the ellipsoid coordinate frames and associated matrices.
 */
__kernel void computeEllipsoidFrames(int numParticles, __global const real4* restrict posq, __global int2* const restrict axisParticleIndices,
        __global const float4* restrict sigParams, __global const float4* restrict scale, __global real* restrict aMatrix,
        __global real* restrict bMatrix, __global real* restrict gMatrix, __global const int* sortedParticles) {
    for (int sortedIndex = get_global_id(0); sortedIndex < numParticles; sortedIndex += get_global_size(0)) {
        // Compute the local coordinate system of the ellipsoid;

        int originalIndex = sortedParticles[sortedIndex];
        real3 pos = posq[originalIndex].xyz;
        int2 axisParticles = axisParticleIndices[originalIndex];
        real3 xdir, ydir, zdir;
        if (axisParticles.x == -1) {
            xdir = (real3) (1, 0, 0);
            ydir = (real3) (0, 1, 0);
        }
        else {
            xdir = pos-posq[axisParticles.x].xyz;
            xdir = normalize(xdir);
            if (axisParticles.y == -1) {
                if (xdir.y > -0.5f && xdir.y < 0.5f)
                    ydir = (real3) (0, 1, 0);
                else
                    ydir = (real3) (1, 0, 0);
            }
            else
                ydir = pos-posq[axisParticles.y].xyz;
            ydir -= xdir*dot(xdir, ydir);
            ydir = normalize(ydir);
        }
        zdir = cross(xdir, ydir);

        // Compute matrices we will need later.

        __global real (*a)[3] = (__global real (*)[3]) (aMatrix+sortedIndex*9);
        __global real (*b)[3] = (__global real (*)[3]) (bMatrix+sortedIndex*9);
        __global real (*g)[3] = (__global real (*)[3]) (gMatrix+sortedIndex*9);
        a[0][0] = xdir.x;
        a[0][1] = xdir.y;
        a[0][2] = xdir.z;
        a[1][0] = ydir.x;
        a[1][1] = ydir.y;
        a[1][2] = ydir.z;
        a[2][0] = zdir.x;
        a[2][1] = zdir.y;
        a[2][2] = zdir.z;
        float4 sig = sigParams[originalIndex];
        float3 r2 = sig.yzw;
        float3 e2 = scale[originalIndex].xyz;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) {
                real belem = 0;
                real gelem = 0;
                for (int k = 0; k < 3; k++) {
                    belem += a[k][i]*e2[k]*a[k][j];
                    gelem += a[k][i]*r2[k]*a[k][j];
                }
                b[i][j] = belem;
                g[i][j] = gelem;
            }
    }
}

/**
 * Find a bounding box for the atoms in each block.
 */
__kernel void findBlockBounds(int numAtoms, real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        __global const int* sortedAtoms, __global const real4* restrict posq, __global real4* restrict sortedPos, __global real4* restrict blockCenter,
        __global real4* restrict blockBoundingBox, __global int* restrict neighborBlockCount) {
    int index = get_global_id(0);
    int base = index*TILE_SIZE;
    while (base < numAtoms) {
        real4 pos = posq[sortedAtoms[base]];
        sortedPos[base] = pos;
#ifdef USE_PERIODIC
        APPLY_PERIODIC_TO_POS(pos)
#endif
        real4 minPos = pos;
        real4 maxPos = pos;
        int last = min(base+TILE_SIZE, numAtoms);
        for (int i = base+1; i < last; i++) {
            pos = posq[sortedAtoms[i]];
            sortedPos[i] = pos;
#ifdef USE_PERIODIC
            real4 center = 0.5f*(maxPos+minPos);
            APPLY_PERIODIC_TO_POS_WITH_CENTER(pos, center)
#endif
            minPos = min(minPos, pos);
            maxPos = max(maxPos, pos);
        }
        real4 blockSize = 0.5f*(maxPos-minPos);
        blockBoundingBox[index] = blockSize;
        blockCenter[index] = 0.5f*(maxPos+minPos);
        index += get_global_size(0);
        base = index*TILE_SIZE;
    }
    if (get_global_id(0) == 0)
        *neighborBlockCount = 0;
}

/**
 * Build a list of neighbors for each atom.
 */
__kernel void findNeighbors(int numAtoms, int maxNeighborBlocks, real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        __global real4* restrict sortedPos, __global real4* restrict blockCenter, __global real4* restrict blockBoundingBox, __global int* restrict neighbors,
        __global int2* restrict neighborIndex, __global int* restrict neighborBlockCount) {
}

typedef struct {
    float4 sig;
    float2 eps;
    real3 pos;
    real a[3][3], b[3][3], g[3][3];
} AtomData;

void loadAtomData(AtomData* data, int sortedIndex, int originalIndex, __global const real4* restrict pos, __global const float4* restrict sigParams,
        __global const float2* restrict epsParams, __global const real* restrict aMatrix, __global const real* restrict bMatrix, __global const real* restrict gMatrix) {
    data->sig = sigParams[originalIndex];
    data->eps = epsParams[originalIndex];
    data->pos = pos[sortedIndex].xyz;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            int k = 9*sortedIndex+3*i+j;
            data->a[i][j] = aMatrix[k];
            data->b[i][j] = bMatrix[k];
            data->g[i][j] = gMatrix[k];
        }
}

real3 matrixVectorProduct(real (*m)[3], real3 v) {
    return (real3) (m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2],
                    m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2],
                    m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2]);
}

real3 vectorMatrixProduct(real3 v, real (*m)[3]) {
    return (real3) (m[0][0]*v[0] + m[1][0]*v[1] + m[2][0]*v[2],
                    m[0][1]*v[0] + m[1][1]*v[1] + m[2][1]*v[2],
                    m[0][2]*v[0] + m[1][2]*v[1] + m[2][2]*v[2]);
}


void matrixSum(real (*result)[3], real (*a)[3], real (*b)[3]) {
    result[0][0] = a[0][0]+b[0][0];
    result[0][1] = a[0][1]+b[0][1];
    result[0][2] = a[0][2]+b[0][2];
    result[1][0] = a[1][0]+b[1][0];
    result[1][1] = a[1][1]+b[1][1];
    result[1][2] = a[1][2]+b[1][2];
    result[2][0] = a[2][0]+b[2][0];
    result[2][1] = a[2][1]+b[2][1];
    result[2][2] = a[2][2]+b[2][2];
}

real determinant(real (*m)[3]) {
    return (m[0][0]*m[1][1]*m[2][2] + m[0][1]*m[1][2]*m[2][0] + m[0][2]*m[1][0]*m[2][1] -
            m[0][0]*m[1][2]*m[2][1] - m[0][1]*m[1][0]*m[2][2] - m[0][2]*m[1][1]*m[2][0]);
}


void matrixInverse(real (*result)[3], real (*m)[3]) {
    real invDet = RECIP(determinant(m));
    result[0][0] = invDet*(m[1][1]*m[2][2] - m[1][2]*m[2][1]);
    result[1][0] = -invDet*(m[1][0]*m[2][2] - m[1][2]*m[2][0]);
    result[2][0] = invDet*(m[1][0]*m[2][1] - m[1][1]*m[2][0]);
    result[0][1] = -invDet*(m[0][1]*m[2][2] - m[0][2]*m[2][1]);
    result[1][1] = invDet*(m[0][0]*m[2][2] - m[0][2]*m[2][0]);
    result[2][1] = -invDet*(m[0][0]*m[2][1] - m[0][1]*m[2][0]);
    result[0][2] = invDet*(m[0][1]*m[1][2] - m[0][2]*m[1][1]);
    result[1][2] = -invDet*(m[0][0]*m[1][2] - m[0][2]*m[1][0]);
    result[2][2] = invDet*(m[0][0]*m[1][1] - m[0][1]*m[1][0]);
}

void computeOneInteraction(AtomData* data1, AtomData* data2, real sigma, real epsilon, real3 dr, real r2, real3* force1, real3* force2, real3* torque1, real3* torque2, real *totalEnergy) {
    real rInv = RSQRT(r2);
    real r = r2*rInv;
    real3 drUnit = dr*rInv;
    
    // Compute the switching function.

    real switchValue = 1, switchDeriv = 0;
    #if USE_LJ_SWITCH
    if (r > LJ_SWITCH_CUTOFF) {
        real x = r-LJ_SWITCH_CUTOFF;
        switchValue = 1+x*x*x*(LJ_SWITCH_C3+x*(LJ_SWITCH_C4+x*LJ_SWITCH_C5));
        switchDeriv = x*x*(3*LJ_SWITCH_C3+x*(4*LJ_SWITCH_C4+x*5*LJ_SWITCH_C5));
    }
    #endif

    // Compute vectors and matrices we'll be needing.

    real B12[3][3], G12[3][3], B12inv[3][3], G12inv[3][3];
    matrixSum(B12, data1->b, data2->b);
    matrixSum(G12, data1->g, data2->g);
    matrixInverse(B12inv, B12);
    matrixInverse(G12inv, G12);
    real detG12 = determinant(G12);

    // Estimate the distance between the ellipsoids and compute the first terms needed for the energy.

    real sigma12 = 1/SQRT(0.5f*dot(drUnit, matrixVectorProduct(G12inv, drUnit)));
    real h12 = r - sigma12;
    real rho = sigma/(h12+sigma);
    real rho2 = rho*rho;
    real rho6 = rho2*rho2*rho2;
    real u = 4*epsilon*(rho6*rho6-rho6);
    real eta = SQRT(2*data1->eps.y*data2->eps.y/detG12);
    real chi = 2*dot(drUnit, matrixVectorProduct(B12inv, drUnit));
    chi *= chi;
    real energy = u*eta*chi;
    
    // Compute the terms needed for the force.

    real3 kappa = matrixVectorProduct(G12inv, dr);
    real3 iota = matrixVectorProduct(B12inv, dr);
    real rInv2 = rInv*rInv;
    real dUSLJdr = 24*epsilon*(2*rho6-1)*rho6*rho/sigma;
    real temp = 0.5f*sigma12*sigma12*sigma12*rInv2;
    real3 dudr = (drUnit + (kappa-drUnit*dot(kappa, drUnit))*temp)*dUSLJdr;
    real3 dchidr = (iota-drUnit*dot(iota, drUnit))*(-8*rInv2*SQRT(chi));
    real3 force = (dchidr*u + dudr*chi)*(eta*switchValue) - drUnit*(energy*switchDeriv);
    *force1 += force;
    *force2 -= force;

    // Compute the terms needed for the torque.

    for (int j = 0; j < 2; j++) {
        real (*a)[3] = (j == 0 ? data1->a : data2->a);
        real (*b)[3] = (j == 0 ? data1->b : data2->b);
        real (*g)[3] = (j == 0 ? data1->g : data2->g);
        float4 sig = (j == 0 ? data1->sig : data2->sig);
        real3 dudq = cross(vectorMatrixProduct(kappa, g), kappa*(temp*dUSLJdr));
        real3 dchidq = cross(vectorMatrixProduct(iota, b), iota)*(-4*rInv2);
        real3 scale = (real3) (sig.y, sig.z, sig.w)*(-0.5f*eta/detG12);
        real d[3][3];
        d[0][0] = scale[0]*(2*a[0][0]*(G12[1][1]*G12[2][2] - G12[1][2]*G12[2][1]) +
                              a[0][2]*(G12[1][2]*G12[0][1] + G12[1][0]*G12[2][1] - G12[1][1]*(G12[0][2] + G12[2][0])) +
                              a[0][1]*(G12[0][2]*G12[2][1] + G12[2][0]*G12[1][2] - G12[2][2]*(G12[0][1] + G12[1][0])));
        d[0][1] = scale[0]*(  a[0][0]*(G12[0][2]*G12[2][1] + G12[2][0]*G12[1][2] - G12[2][2]*(G12[0][1] + G12[1][0])) +
                            2*a[0][1]*(G12[0][0]*G12[2][2] - G12[2][0]*G12[0][2]) +
                              a[0][2]*(G12[1][0]*G12[0][2] + G12[2][0]*G12[0][1] - G12[0][0]*(G12[1][2] + G12[2][1])));
        d[0][2] = scale[0]*(  a[0][0]*(G12[0][1]*G12[1][2] + G12[1][0]*G12[2][1] - G12[1][1]*(G12[0][2] + G12[2][0])) +
                              a[0][1]*(G12[1][0]*G12[0][2] + G12[2][0]*G12[0][1] - G12[0][0]*(G12[1][2] + G12[2][1])) +
                            2*a[0][2]*(G12[1][1]*G12[0][0] - G12[1][0]*G12[0][1]));
        d[1][0] = scale[1]*(2*a[1][0]*(G12[1][1]*G12[2][2] - G12[1][2]*G12[2][1]) +
                              a[1][1]*(G12[0][2]*G12[2][1] + G12[2][0]*G12[1][2] - G12[2][2]*(G12[0][1] + G12[1][0])) +
                              a[1][2]*(G12[1][2]*G12[0][1] + G12[1][0]*G12[2][1] - G12[1][1]*(G12[0][2] + G12[2][0])));
        d[1][1] = scale[1]*(  a[1][0]*(G12[0][2]*G12[2][1] + G12[2][0]*G12[1][2] - G12[2][2]*(G12[0][1] + G12[1][0])) +
                            2*a[1][1]*(G12[2][2]*G12[0][0] - G12[2][0]*G12[0][2]) +
                              a[1][2]*(G12[1][0]*G12[0][2] + G12[0][1]*G12[2][0] - G12[0][0]*(G12[1][2] + G12[2][1])));
        d[1][2] = scale[1]*(  a[1][0]*(G12[0][1]*G12[1][2] + G12[1][0]*G12[2][1] - G12[1][1]*(G12[0][2] + G12[2][0])) +
                              a[1][1]*(G12[1][0]*G12[0][2] + G12[0][1]*G12[2][0] - G12[0][0]*(G12[1][2] + G12[2][1])) +
                            2*a[1][2]*(G12[1][1]*G12[0][0] - G12[1][0]*G12[0][1]));
        d[2][0] = scale[2]*(2*a[2][0]*(G12[1][1]*G12[2][2] - G12[2][1]*G12[1][2]) +
                              a[2][1]*(G12[0][2]*G12[2][1] + G12[1][2]*G12[2][0] - G12[2][2]*(G12[0][1] + G12[1][0])) +
                              a[2][2]*(G12[0][1]*G12[1][2] + G12[2][1]*G12[1][0] - G12[1][1]*(G12[0][2] + G12[2][0])));
        d[2][1] = scale[2]*(  a[2][0]*(G12[0][2]*G12[2][1] + G12[1][2]*G12[2][0] - G12[2][2]*(G12[0][1] + G12[1][0])) +
                            2*a[2][1]*(G12[0][0]*G12[2][2] - G12[0][2]*G12[2][0]) +
                              a[2][2]*(G12[1][0]*G12[0][2] + G12[0][1]*G12[2][0] - G12[0][0]*(G12[1][2] + G12[2][1])));
        d[2][2] = scale[2]*(  a[2][0]*(G12[0][1]*G12[1][2] + G12[2][1]*G12[1][0] - G12[1][1]*(G12[0][2] + G12[2][0])) +
                              a[2][1]*(G12[1][0]*G12[0][2] + G12[2][0]*G12[0][1] - G12[0][0]*(G12[1][2] + G12[2][1])) +
                            2*a[2][2]*(G12[1][1]*G12[0][0] - G12[1][0]*G12[0][1]));
        real3 detadq;
        for (int i = 0; i < 3; i++)
            detadq += cross((real3) (a[i][0], a[i][1], a[i][2]), (real3) (d[i][0], d[i][1], d[i][2]));
        real3 torque = (dchidq*(u*eta) + detadq*(u*chi) + dudq*(eta*chi))*switchValue;
        *(j == 0 ? torque1 : torque2) -= torque;
    }
    *totalEnergy += switchValue*energy;
}

/**
 * Compute the interactions.
 */
__kernel void computeForce(
#ifdef SUPPORTS_64_BIT_ATOMICS
        __global long* restrict forceBuffers, __global long* restrict torqueBuffers,
#else
        __global real4* restrict forceBuffers, __global real4* restrict torqueBuffers,
#endif
        int numAtoms, int numExceptions, __global mixed* restrict energyBuffer, __global const real4* restrict pos,
        __global const float4* restrict sigParams, __global const float2* restrict epsParams, __global const int* restrict sortedAtoms,
        __global const real* restrict aMatrix, __global const real* restrict bMatrix, __global const real* restrict gMatrix,
        __global const int* restrict exclusions, __global const int* restrict exclusionStartIndex,
        __global const int4* restrict exceptionParticles, __global const float2* restrict exceptionParams
#ifdef USE_CUTOFF
        __global const unsigned int* restrict interactionCount, real4 periodicBoxSize, real4 invPeriodicBoxSize, 
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, unsigned int maxTiles, __global const real4* restrict blockCenter,
        __global const real4* restrict blockSize, __global const int* restrict interactingAtoms,
#endif
        ) {
    const unsigned int warp = get_global_id(0)/TILE_SIZE;
    mixed energy = 0;
    for (int atom1 = get_global_id(0); atom1 < numAtoms; atom1 += get_global_size(0)) {
        // Load parameters for atom1.
        
        int index1 = sortedAtoms[atom1];
        AtomData data1;
        loadAtomData(&data1, atom1, index1, pos, sigParams, epsParams, aMatrix, bMatrix, gMatrix);
        real3 force1 = 0.0f;
        real3 torque1 = 0.0f;
        int nextExclusion = exclusionStartIndex[atom1];
        int lastExclusion = exclusionStartIndex[atom1+1];
        for (int atom2 = atom1+1; atom2 < numAtoms; atom2++) {
            // Skip over excluded interactions.
            
            int index2 = sortedAtoms[atom2];
            if (nextExclusion < lastExclusion && exclusions[nextExclusion] == index2) {
                nextExclusion++;
                continue;
            }
            
            // Load parameters for atom2.
            
            AtomData data2;
            loadAtomData(&data2, atom2, index2, pos, sigParams, epsParams, aMatrix, bMatrix, gMatrix);
            real3 force2 = 0.0f;
            real3 torque2 = 0.0f;
            
            // Compute the interaction.
            
            real3 delta = data1.pos-data2.pos;
#ifdef USE_PERIODIC
            APPLY_PERIODIC_TO_DELTA(delta)
#endif
            real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
            if (r2 < CUTOFF_SQUARED) {
#endif
                real sigma = data1.sig.x+data2.sig.x;
                real epsilon = data1.eps.x*data2.eps.x;
                computeOneInteraction(&data1, &data2, sigma, epsilon, delta, r2, &force1, &force2, &torque1, &torque2, &energy);
#ifdef SUPPORTS_64_BIT_ATOMICS
                atom_add(&forceBuffers[index2], (long) (force2.x*0x100000000));
                atom_add(&forceBuffers[index2+PADDED_NUM_ATOMS], (long) (force2.y*0x100000000));
                atom_add(&forceBuffers[index2+2*PADDED_NUM_ATOMS], (long) (force2.z*0x100000000));
                atom_add(&torqueBuffers[index2], (long) (torque2.x*0x100000000));
                atom_add(&torqueBuffers[index2+PADDED_NUM_ATOMS], (long) (torque2.y*0x100000000));
                atom_add(&torqueBuffers[index2+2*PADDED_NUM_ATOMS], (long) (torque2.z*0x100000000));
#else
                unsigned int offset = index2 + warp*PADDED_NUM_ATOMS;
                forceBuffers[offset].xyz += force2.xyz;
                torqueBuffers[offset].xyz += torque2.xyz;
#endif
#ifdef USE_CUTOFF
            }
#endif
        }
#ifdef SUPPORTS_64_BIT_ATOMICS
        atom_add(&forceBuffers[index1], (long) (force1.x*0x100000000));
        atom_add(&forceBuffers[index1+PADDED_NUM_ATOMS], (long) (force1.y*0x100000000));
        atom_add(&forceBuffers[index1+2*PADDED_NUM_ATOMS], (long) (force1.z*0x100000000));
        atom_add(&torqueBuffers[index1], (long) (torque1.x*0x100000000));
        atom_add(&torqueBuffers[index1+PADDED_NUM_ATOMS], (long) (torque1.y*0x100000000));
        atom_add(&torqueBuffers[index1+2*PADDED_NUM_ATOMS], (long) (torque1.z*0x100000000));
#else
        unsigned int offset = index1 + warp*PADDED_NUM_ATOMS;
        forceBuffers[offset].xyz += force1.xyz;
        torqueBuffers[offset].xyz += torque1.xyz;
#endif
    }
    
    // Now compute exceptions.
    
    for (int index = get_global_id(0); index < numExceptions; index += get_global_size(0)) {
        int4 atomIndices = exceptionParticles[index];
        real2 params = exceptionParams[index];
        int index1 = atomIndices.x, index2 = atomIndices.y;
        int atom1 = atomIndices.z, atom2 = atomIndices.w;
        AtomData data1, data2;
        loadAtomData(&data1, atom1, index1, pos, sigParams, epsParams, aMatrix, bMatrix, gMatrix);
        loadAtomData(&data2, atom2, index2, pos, sigParams, epsParams, aMatrix, bMatrix, gMatrix);
        real3 force1 = 0, force2 = 0;
        real3 torque1 = 0, torque2 = 0;
        real3 delta = data1.pos-data2.pos;
        real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
#ifdef USE_CUTOFF
        if (r2 < CUTOFF_SQUARED) {
#endif
            computeOneInteraction(&data1, &data2, params.x, params.y, delta, r2, &force1, &force2, &torque1, &torque2, &energy);
#ifdef SUPPORTS_64_BIT_ATOMICS
            atom_add(&forceBuffers[index1], (long) (force1.x*0x100000000));
            atom_add(&forceBuffers[index1+PADDED_NUM_ATOMS], (long) (force1.y*0x100000000));
            atom_add(&forceBuffers[index1+2*PADDED_NUM_ATOMS], (long) (force1.z*0x100000000));
            atom_add(&forceBuffers[index2], (long) (force2.x*0x100000000));
            atom_add(&forceBuffers[index2+PADDED_NUM_ATOMS], (long) (force2.y*0x100000000));
            atom_add(&forceBuffers[index2+2*PADDED_NUM_ATOMS], (long) (force2.z*0x100000000));
            atom_add(&torqueBuffers[index1], (long) (torque1.x*0x100000000));
            atom_add(&torqueBuffers[index1+PADDED_NUM_ATOMS], (long) (torque1.y*0x100000000));
            atom_add(&torqueBuffers[index1+2*PADDED_NUM_ATOMS], (long) (torque1.z*0x100000000));
            atom_add(&torqueBuffers[index2], (long) (torque2.x*0x100000000));
            atom_add(&torqueBuffers[index2+PADDED_NUM_ATOMS], (long) (torque2.y*0x100000000));
            atom_add(&torqueBuffers[index2+2*PADDED_NUM_ATOMS], (long) (torque2.z*0x100000000));
#else
            int offset = index1 + warp*PADDED_NUM_ATOMS;
            forceBuffers[offset].xyz += force1.xyz;
            torqueBuffers[offset].xyz += torque1.xyz;
            offset = index2 + warp*PADDED_NUM_ATOMS;
            forceBuffers[offset].xyz += force2.xyz;
            torqueBuffers[offset].xyz += torque2.xyz;
#endif
#ifdef USE_CUTOFF
        }
#endif
    }
    energyBuffer[get_global_id(0)] += energy;
}

/**
 * Convert the torques to forces on the connected particles.
 */
__kernel void applyTorques(
#ifdef SUPPORTS_64_BIT_ATOMICS
        __global long* restrict forceBuffers, __global long* restrict torqueBuffers,
#else
        __global real4* restrict forceBuffers, __global real4* restrict torqueBuffers, int numBuffers,
#endif
        int numParticles, __global const real4* restrict posq, __global int2* const restrict axisParticleIndices,
        __global const int* sortedParticles) {
    const unsigned int warp = get_global_id(0)/TILE_SIZE;
    for (int sortedIndex = get_global_id(0); sortedIndex < numParticles; sortedIndex += get_global_size(0)) {
        int originalIndex = sortedParticles[sortedIndex];
        real3 pos = posq[originalIndex].xyz;
        int2 axisParticles = axisParticleIndices[originalIndex];
        if (axisParticles.x != -1) {
            // Load the torque.

#ifdef SUPPORTS_64_BIT_ATOMICS
            real scale = 1/(real) 0x100000000;
            real3 torque = (real3) (scale*torqueBuffers[originalIndex], scale*torqueBuffers[originalIndex+PADDED_NUM_ATOMS], scale*torqueBuffers[originalIndex+2*PADDED_NUM_ATOMS]);
#else
            real3 torque = 0;
            for (int i = 0; i < numBuffers; i++)
                torque += torqueBuffers[originalIndex+i*PADDED_NUM_ATOMS].xyz;
#endif
            real3 force = 0, xforce = 0, yforce = 0;

            // Apply a force to the x particle.
            
            real3 dx = posq[axisParticles.x].xyz-pos;
            real dx2 = dot(dx, dx);
            real3 f = cross(torque, dx)/dx2;
            xforce += f;
            force -= f;
            if (axisParticles.y != -1) {
                // Apply a force to the y particle.  This is based on the component of the torque
                // that was not already applied to the x particle.
                
                real3 dy = posq[axisParticles.y].xyz-pos;
                real dy2 = dot(dy, dy);
                real3 torque2 = dx*dot(torque, dx)/dx2;
                f = cross(torque2, dy)/dy2;
                yforce += f;
                force -= f;
            }
#ifdef SUPPORTS_64_BIT_ATOMICS
            atom_add(&forceBuffers[originalIndex], (long) (force.x*0x100000000));
            atom_add(&forceBuffers[originalIndex+PADDED_NUM_ATOMS], (long) (force.y*0x100000000));
            atom_add(&forceBuffers[originalIndex+2*PADDED_NUM_ATOMS], (long) (force.z*0x100000000));
            atom_add(&forceBuffers[axisParticles.x], (long) (xforce.x*0x100000000));
            atom_add(&forceBuffers[axisParticles.x+PADDED_NUM_ATOMS], (long) (xforce.y*0x100000000));
            atom_add(&forceBuffers[axisParticles.x+2*PADDED_NUM_ATOMS], (long) (xforce.z*0x100000000));
            if (axisParticles.y != -1) {
                atom_add(&forceBuffers[axisParticles.y], (long) (yforce.x*0x100000000));
                atom_add(&forceBuffers[axisParticles.y+PADDED_NUM_ATOMS], (long) (yforce.y*0x100000000));
                atom_add(&forceBuffers[axisParticles.y+2*PADDED_NUM_ATOMS], (long) (yforce.z*0x100000000));
            }
#else
            forceBuffers[originalIndex + warp*PADDED_NUM_ATOMS].xyz += force.xyz;
            forceBuffers[axisParticles.x + warp*PADDED_NUM_ATOMS].xyz += xforce.xyz;
            if (axisParticles.y != -1)
                forceBuffers[axisParticles.y + warp*PADDED_NUM_ATOMS].xyz += yforce.xyz;
#endif
        }
    }
}