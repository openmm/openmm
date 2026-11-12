/**
 * Generate random numbers
 */
KERNEL void generateRandomNumbers(int numValues, GLOBAL float4* RESTRICT random, GLOBAL uint4* RESTRICT seed) {
    int index = GLOBAL_ID;
    uint4 state = seed[index];
    unsigned int carry = 0;
    while (index < numValues) {
        float4 value;

        // Generate first two values.

        state.x = state.x * 69069 + 1;
        state.y ^= state.y << 13;
        state.y ^= state.y >> 17;
        state.y ^= state.y << 5;
        unsigned int k = (state.z >> 2) + (state.w >> 3) + (carry >> 2);
        unsigned int m = state.w + state.w + state.z + carry;
        state.z = state.w;
        state.w = m;
        carry = k >> 30;
        float x1 = (float)max(state.x + state.y + state.w, 0x00000001u) / (float)0xffffffff;
        state.x = state.x * 69069 + 1;
        state.y ^= state.y << 13;
        state.y ^= state.y >> 17;
        state.y ^= state.y << 5;
        x1 = SQRT(-2.0f * LOG(x1));
        k = (state.z >> 2) + (state.w >> 3) + (carry >> 2);
        m = state.w + state.w + state.z + carry;
        state.z = state.w;
        state.w = m;
        carry = k >> 30;
        float x2 = (float)(state.x + state.y + state.w) / (float)0xffffffff;
        value.x = x1 * COS(2.0f * 3.14159265f * x2);
        value.y = x1 * SIN(2.0f * 3.14159265f * x2);

        // Generate next two values.

        state.x = state.x * 69069 + 1;
        state.y ^= state.y << 13;
        state.y ^= state.y >> 17;
        state.y ^= state.y << 5;
        k = (state.z >> 2) + (state.w >> 3) + (carry >> 2);
        m = state.w + state.w + state.z + carry;
        state.z = state.w;
        state.w = m;
        carry = k >> 30;
        float x3 = (float)max(state.x + state.y + state.w, 0x00000001u) / (float)0xffffffff;
        state.x = state.x * 69069 + 1;
        state.y ^= state.y << 13;
        state.y ^= state.y >> 17;
        state.y ^= state.y << 5;
        x3 = SQRT(-2.0f * LOG(x3));
        k = (state.z >> 2) + (state.w >> 3) + (carry >> 2);
        m = state.w + state.w + state.z + carry;
        state.z = state.w;
        state.w = m;
        carry = k >> 30;
        float x4 = (float)(state.x + state.y + state.w) / (float)0xffffffff;
        value.z = x3 * COS(2.0f * 3.14159265f * x4);
        value.w = x3 * SIN(2.0f * 3.14159265f * x4);

        // Record the values.

        random[index] = value;
        index += GLOBAL_SIZE;
    }
    seed[GLOBAL_ID] = state;
}

/**
 * Enforce constraints on SHAKE clusters
 */
KERNEL void applyShakeToPositions(int numClusters, mixed tol, GLOBAL const real4* RESTRICT oldPos,
        GLOBAL mixed4* RESTRICT posDelta, GLOBAL const int4* RESTRICT clusterAtoms, GLOBAL const float4* RESTRICT clusterParams
#ifdef USE_MIXED_PRECISION
        , GLOBAL const real4* RESTRICT posqCorrection
#endif
        ) {
#ifndef USE_MIXED_PRECISION
        GLOBAL real4* posqCorrection = 0;
#endif
    int index = GLOBAL_ID;
    while (index < numClusters) {
        // Load the data for this cluster.

        int4 atoms = clusterAtoms[index];
        float4 params = clusterParams[index];
        mixed4 pos = loadPos(oldPos, posqCorrection, atoms.x);
        mixed4 xpi = posDelta[atoms.x];
        mixed4 pos1 = loadPos(oldPos, posqCorrection, atoms.y);
        mixed4 xpj1 = posDelta[atoms.y];
        mixed4 pos2 = make_mixed4(0);
        mixed4 xpj2 = make_mixed4(0);
        if (atoms.z != -1) {
            pos2 = loadPos(oldPos, posqCorrection, atoms.z);
            xpj2 = posDelta[atoms.z];
        }
        mixed4 pos3 = make_mixed4(0);
        mixed4 xpj3 = make_mixed4(0);
        if (atoms.w != -1) {
            pos3 = loadPos(oldPos, posqCorrection, atoms.w);
            xpj3 = posDelta[atoms.w];
        }

        applyShakeToPositionsCluster(atoms, params, tol, pos, pos1, pos2, pos3, &xpi, &xpj1, &xpj2, &xpj3);

        // Record the new positions.

        posDelta[atoms.x] = xpi;
        posDelta[atoms.y] = xpj1;
        if (atoms.z != -1)
            posDelta[atoms.z] = xpj2;
        if (atoms.w != -1)
            posDelta[atoms.w] = xpj3;
        index += GLOBAL_SIZE;
    }
}

/**
 * Enforce velocity constraints on SHAKE clusters
 */
KERNEL void applyShakeToVelocities(int numClusters, mixed tol, GLOBAL const real4* RESTRICT oldPos,
        GLOBAL mixed4* RESTRICT posDelta, GLOBAL const int4* RESTRICT clusterAtoms, GLOBAL const float4* RESTRICT clusterParams
#ifdef USE_MIXED_PRECISION
        , GLOBAL const real4* RESTRICT posqCorrection
#endif
        ) {
#ifndef USE_MIXED_PRECISION
        GLOBAL real4* posqCorrection = 0;
#endif
    int index = GLOBAL_ID;
    while (index < numClusters) {
        // Load the data for this cluster.

        int4 atoms = clusterAtoms[index];
        float4 params = clusterParams[index];
        mixed4 pos = loadPos(oldPos, posqCorrection, atoms.x);
        mixed4 xpi = posDelta[atoms.x];
        mixed4 pos1 = loadPos(oldPos, posqCorrection, atoms.y);
        mixed4 xpj1 = posDelta[atoms.y];
        mixed4 pos2 = make_mixed4(0);
        mixed4 xpj2 = make_mixed4(0);
        if (atoms.z != -1) {
            pos2 = loadPos(oldPos, posqCorrection, atoms.z);
            xpj2 = posDelta[atoms.z];
        }
        mixed4 pos3 = make_mixed4(0);
        mixed4 xpj3 = make_mixed4(0);
        if (atoms.w != -1) {
            pos3 = loadPos(oldPos, posqCorrection, atoms.w);
            xpj3 = posDelta[atoms.w];
        }

        applyShakeToVelocitiesCluster(atoms, params, tol, pos, pos1, pos2, pos3, &xpi, &xpj1, &xpj2, &xpj3);

        // Record the new positions.

        posDelta[atoms.x] = xpi;
        posDelta[atoms.y] = xpj1;
        if (atoms.z != -1)
            posDelta[atoms.z] = xpj2;
        if (atoms.w != -1)
            posDelta[atoms.w] = xpj3;
        index += GLOBAL_SIZE;
    }
}

/**
 * Enforce constraints on SETTLE clusters
 */
KERNEL void applySettleToPositions(int numClusters, mixed tol, GLOBAL const real4* RESTRICT oldPos,
        GLOBAL mixed4* RESTRICT posDelta, GLOBAL const mixed4* RESTRICT velm, GLOBAL const int4* RESTRICT clusterAtoms,
        GLOBAL const float2* RESTRICT clusterParams
#ifdef USE_MIXED_PRECISION
        , GLOBAL const real4* RESTRICT posqCorrection
#endif
        ) {
#ifndef USE_MIXED_PRECISION
        GLOBAL real4* posqCorrection = 0;
#endif
    int index = GLOBAL_ID;
    while (index < numClusters) {
        // Load the data for this cluster.

        int4 atoms = clusterAtoms[index];
        float2 params = clusterParams[index];
        mixed4 apos0 = loadPos(oldPos, posqCorrection, atoms.x);
        mixed4 xp0 = posDelta[atoms.x];
        mixed4 apos1 = loadPos(oldPos, posqCorrection, atoms.y);
        mixed4 xp1 = posDelta[atoms.y];
        mixed4 apos2 = loadPos(oldPos, posqCorrection, atoms.z);
        mixed4 xp2 = posDelta[atoms.z];
        mixed m0 = 1/velm[atoms.x].w;
        mixed m1 = 1/velm[atoms.y].w;
        mixed m2 = 1/velm[atoms.z].w;

        applySettleToPositionsCluster(params, apos0, apos1, apos2, m0, m1, m2, &xp0, &xp1, &xp2);

        // Record the new positions.

        posDelta[atoms.x] = xp0;
        posDelta[atoms.y] = xp1;
        posDelta[atoms.z] = xp2;
        index += GLOBAL_SIZE;
    }
}

/**
 * Enforce velocity constraints on SETTLE clusters
 */
KERNEL void applySettleToVelocities(int numClusters, mixed tol, GLOBAL const real4* RESTRICT oldPos,
        GLOBAL mixed4* RESTRICT posDelta, GLOBAL mixed4* RESTRICT velm, GLOBAL const int4* RESTRICT clusterAtoms,
        GLOBAL const float2* RESTRICT clusterParams
#ifdef USE_MIXED_PRECISION
        , GLOBAL const real4* RESTRICT posqCorrection
#endif
        ) {
#ifndef USE_MIXED_PRECISION
        GLOBAL real4* posqCorrection = 0;
#endif
    for (int index = GLOBAL_ID; index < numClusters; index += GLOBAL_SIZE) {
        // Load the data for this cluster.

        int4 atoms = clusterAtoms[index];
        mixed4 apos0 = loadPos(oldPos, posqCorrection, atoms.x);
        mixed4 apos1 = loadPos(oldPos, posqCorrection, atoms.y);
        mixed4 apos2 = loadPos(oldPos, posqCorrection, atoms.z);
        mixed4 v0 = velm[atoms.x];
        mixed4 v1 = velm[atoms.y];
        mixed4 v2 = velm[atoms.z];
        
        applySettleToVelocitiesCluster(apos0, apos1, apos2, &v0, &v1, &v2);
        velm[atoms.x] = v0;
        velm[atoms.y] = v1;
        velm[atoms.z] = v2;
    }
}

/**
 * Compute the direction each CCMA constraint is pointing in.  This is called once at the beginning of constraint evaluation.
 */
DEVICE void computeCCMAConstraintDirections(GLOBAL const int2* RESTRICT constraintAtoms, GLOBAL mixed4* RESTRICT constraintDistance,
        GLOBAL const real4* RESTRICT atomPositions
#ifdef USE_MIXED_PRECISION
        , GLOBAL const real4* RESTRICT posqCorrection
#endif
        ) {
#ifndef USE_MIXED_PRECISION
        GLOBAL real4* posqCorrection = 0;
#endif
    for (int index = GLOBAL_ID; index < NUM_CCMA_CONSTRAINTS; index += GLOBAL_SIZE) {
        // Compute the direction for this constraint.

        int2 atoms = constraintAtoms[index];
        mixed4 dir = constraintDistance[index];
        mixed4 oldPos1 = loadPos(atomPositions, posqCorrection, atoms.x);
        mixed4 oldPos2 = loadPos(atomPositions, posqCorrection, atoms.y);
        dir.x = oldPos1.x-oldPos2.x;
        dir.y = oldPos1.y-oldPos2.y;
        dir.z = oldPos1.z-oldPos2.z;
        constraintDistance[index] = dir;
    }
}

KERNEL void computeCCMAConstraintDirectionsKernel(GLOBAL const int2* RESTRICT constraintAtoms, GLOBAL mixed4* RESTRICT constraintDistance,
        GLOBAL const real4* RESTRICT atomPositions, GLOBAL int* RESTRICT converged
#ifdef USE_MIXED_PRECISION
        , GLOBAL const real4* RESTRICT posqCorrection
#endif
        ) {
#ifdef USE_MIXED_PRECISION
    computeCCMAConstraintDirections(constraintAtoms, constraintDistance, atomPositions, posqCorrection);
#else
    computeCCMAConstraintDirections(constraintAtoms, constraintDistance, atomPositions);
#endif
    if (GLOBAL_ID == 0) {
        converged[0] = 1;
        converged[1] = 0;
    }
}

/**
 * Compute the force applied by each CCMA position constraint.
 */
DEVICE void computeCCMAPositionConstraintForce(GLOBAL const int2* RESTRICT constraintAtoms, GLOBAL const mixed4* RESTRICT constraintDistance,
        GLOBAL const mixed4* RESTRICT atomPositions, GLOBAL const mixed* RESTRICT reducedMass, GLOBAL mixed* RESTRICT delta1,
        mixed tol, int iteration, LOCAL_ARG int* groupConverged) {
    if (LOCAL_ID == 0)
        *groupConverged = 1;
    SYNC_THREADS;
    mixed lowerTol = 1-2*tol+tol*tol;
    mixed upperTol = 1+2*tol+tol*tol;
    bool threadConverged = true;
    for (int index = GLOBAL_ID; index < NUM_CCMA_CONSTRAINTS; index += GLOBAL_SIZE) {
        // Compute the force due to this constraint.

        int2 atoms = constraintAtoms[index];
        mixed4 dir = constraintDistance[index];
        mixed4 rp_ij = atomPositions[atoms.x]-atomPositions[atoms.y];
        rp_ij.x += dir.x;
        rp_ij.y += dir.y;
        rp_ij.z += dir.z;
        mixed rrpr = rp_ij.x*dir.x + rp_ij.y*dir.y + rp_ij.z*dir.z;
        mixed d_ij2 = dir.x*dir.x + dir.y*dir.y + dir.z*dir.z;
        mixed rp2 = rp_ij.x*rp_ij.x + rp_ij.y*rp_ij.y + rp_ij.z*rp_ij.z;
        mixed dist2 = dir.w*dir.w;
        mixed diff = dist2 - rp2;
        delta1[index] = (rrpr > d_ij2*1e-6f ? reducedMass[index]*diff/rrpr : 0);
        threadConverged &= (rp2 > lowerTol*dist2 && rp2 < upperTol*dist2);
    }
    if (*groupConverged && !threadConverged)
        *groupConverged = 0;
}

KERNEL void computeCCMAPositionConstraintForceKernel(GLOBAL const int2* RESTRICT constraintAtoms, GLOBAL const mixed4* RESTRICT constraintDistance,
        GLOBAL const mixed4* RESTRICT atomPositions, GLOBAL const mixed* RESTRICT reducedMass, GLOBAL mixed* RESTRICT delta1,
        GLOBAL int* RESTRICT converged, GLOBAL int* RESTRICT hostConvergedFlag, mixed tol, int iteration) {
    LOCAL int groupConverged;
    if (converged[1-iteration%2]) {
        if (GLOBAL_ID == 0) {
            converged[iteration%2] = 1;
            hostConvergedFlag[0] = 1;
        }
        return; // The constraint iteration has already converged.
    }
    computeCCMAPositionConstraintForce(constraintAtoms, constraintDistance, atomPositions, reducedMass,
            delta1, tol, iteration, &groupConverged);
    SYNC_THREADS;
    if (LOCAL_ID == 0 && !groupConverged)
        converged[iteration%2] = 0;
}

/**
 * Compute the force applied by each CCMA velocity constraint.
 */
DEVICE void computeCCMAVelocityConstraintForce(GLOBAL const int2* RESTRICT constraintAtoms, GLOBAL const mixed4* RESTRICT constraintDistance,
        GLOBAL const mixed4* RESTRICT atomPositions, GLOBAL const mixed* RESTRICT reducedMass, GLOBAL mixed* RESTRICT delta1,
        mixed tol, int iteration, LOCAL_ARG int* groupConverged) {
    if (LOCAL_ID == 0)
        *groupConverged = 1;
    SYNC_THREADS;
    bool threadConverged = true;
    for (int index = GLOBAL_ID; index < NUM_CCMA_CONSTRAINTS; index += GLOBAL_SIZE) {
        // Compute the force due to this constraint.

        int2 atoms = constraintAtoms[index];
        mixed4 dir = constraintDistance[index];
        mixed4 rp_ij = atomPositions[atoms.x]-atomPositions[atoms.y];
        mixed rrpr = rp_ij.x*dir.x + rp_ij.y*dir.y + rp_ij.z*dir.z;
        mixed d_ij2 = dir.x*dir.x + dir.y*dir.y + dir.z*dir.z;
        delta1[index] = -2*reducedMass[index]*rrpr/d_ij2;
        threadConverged &= (fabs(delta1[index]) <= tol);
    }
    if (*groupConverged && !threadConverged)
        *groupConverged = 0;
}

KERNEL void computeCCMAVelocityConstraintForceKernel(GLOBAL const int2* RESTRICT constraintAtoms, GLOBAL const mixed4* RESTRICT constraintDistance,
        GLOBAL const mixed4* RESTRICT atomPositions, GLOBAL const mixed* RESTRICT reducedMass, GLOBAL mixed* RESTRICT delta1,
        GLOBAL int* RESTRICT converged, GLOBAL int* RESTRICT hostConvergedFlag, mixed tol, int iteration) {
    LOCAL int groupConverged;
    if (converged[1-iteration%2]) {
        if (GROUP_ID == 0 && LOCAL_ID == 0) {
            converged[iteration%2] = 1;
            hostConvergedFlag[0] = 1;
        }
        return; // The constraint iteration has already converged.
    }
    computeCCMAVelocityConstraintForce(constraintAtoms, constraintDistance, atomPositions, reducedMass,
            delta1, tol, iteration, &groupConverged);
    if (LOCAL_ID == 0 && !groupConverged)
        converged[iteration%2] = 0;
}

/**
 * Multiply the vector of CCMA constraint forces by the constraint matrix.
 */
DEVICE void multiplyByCCMAConstraintMatrix(GLOBAL const mixed* RESTRICT delta1, GLOBAL mixed* RESTRICT delta2, GLOBAL const int* RESTRICT constraintMatrixColumn,
        GLOBAL const mixed* RESTRICT constraintMatrixValue, int iteration) {
    // Multiply by the inverse constraint matrix.

    for (int index = GLOBAL_ID; index < NUM_CCMA_CONSTRAINTS; index += GLOBAL_SIZE) {
        mixed sum = 0;
        for (int i = 0; ; i++) {
            int element = index+i*NUM_CCMA_CONSTRAINTS;
            int column = constraintMatrixColumn[element];
            if (column >= NUM_CCMA_CONSTRAINTS)
                break;
            sum += delta1[column]*constraintMatrixValue[element];
        }
        delta2[index] = sum;
    }
}

KERNEL void multiplyByCCMAConstraintMatrixKernel(GLOBAL const mixed* RESTRICT delta1, GLOBAL mixed* RESTRICT delta2, GLOBAL const int* RESTRICT constraintMatrixColumn,
        GLOBAL const mixed* RESTRICT constraintMatrixValue, GLOBAL const int* RESTRICT converged, int iteration) {
    if (converged[iteration%2])
        return; // The constraint iteration has already converged.
    multiplyByCCMAConstraintMatrix(delta1, delta2, constraintMatrixColumn, constraintMatrixValue, iteration);
}

/**
 * Update the atom positions based on CCMA constraint forces.
 */
DEVICE void updateCCMAAtomPositions(GLOBAL const int* RESTRICT atoms, GLOBAL const int* RESTRICT numAtomConstraints, GLOBAL const int* RESTRICT atomConstraints,
        GLOBAL const mixed4* RESTRICT constraintDistance, GLOBAL mixed4* RESTRICT atomPositions, GLOBAL const mixed4* RESTRICT velm,
        GLOBAL const mixed* RESTRICT delta1, GLOBAL const mixed* RESTRICT delta2, int iteration) {
    mixed damping = (iteration < 2 ? 0.5f : 1.0f);
    for (int i = GLOBAL_ID; i < NUM_CCMA_ATOMS; i += GLOBAL_SIZE) {
        // Compute the new position of this atom.

        int index = atoms[i];
        mixed4 atomPos = atomPositions[index];
        mixed invMass = velm[index].w;
        int num = numAtomConstraints[index];
        for (int i = 0; i < num; i++) {
            int constraint = atomConstraints[index+i*NUM_ATOMS];
            bool forward = (constraint > 0);
            constraint = (forward ? constraint-1 : -constraint-1);
            mixed constraintForce = damping*invMass*delta2[constraint];
            constraintForce = (forward ? constraintForce : -constraintForce);
            mixed4 dir = constraintDistance[constraint];
            atomPos.x += constraintForce*dir.x;
            atomPos.y += constraintForce*dir.y;
            atomPos.z += constraintForce*dir.z;
        }
        atomPositions[index] = atomPos;
    }
}

KERNEL void updateCCMAAtomPositionsKernel(GLOBAL const int* RESTRICT atoms, GLOBAL const int* RESTRICT numAtomConstraints, GLOBAL const int* RESTRICT atomConstraints,
        GLOBAL const mixed4* RESTRICT constraintDistance, GLOBAL mixed4* RESTRICT atomPositions, GLOBAL const mixed4* RESTRICT velm,
        GLOBAL const mixed* RESTRICT delta1, GLOBAL const mixed* RESTRICT delta2, GLOBAL int* RESTRICT converged, int iteration) {
    if (GROUP_ID == 0 && LOCAL_ID == 0)
        converged[1-iteration%2] = 1;
    if (converged[iteration%2])
        return; // The constraint iteration has already converged.
    updateCCMAAtomPositions(atoms, numAtomConstraints, atomConstraints, constraintDistance, atomPositions, velm,
            delta1, delta2, iteration);
}

/**
 * Run the entire CCMA iteration within a single kernel.  This has far less overhead than
 * using multiple kernels, but requires the calculation to use only a single workgroup.
 * That makes it faster for small numbers of constraints, but slower for large numbers.
 */
KERNEL void runCCMA(int constrainVelocities, GLOBAL const int* RESTRICT atoms, GLOBAL const int* RESTRICT numAtomConstraints, GLOBAL const int* RESTRICT atomConstraints,
        GLOBAL const int2* RESTRICT constraintAtoms, GLOBAL mixed4* RESTRICT constraintDistance, GLOBAL const real4* RESTRICT atomPositions,
        GLOBAL mixed4* RESTRICT velm, GLOBAL mixed4* RESTRICT posDelta, GLOBAL const mixed* RESTRICT reducedMass,
        GLOBAL mixed* RESTRICT delta1, GLOBAL mixed* RESTRICT delta2, GLOBAL const int* RESTRICT constraintMatrixColumn,
        GLOBAL const mixed* RESTRICT constraintMatrixValue, mixed tol
#ifdef USE_MIXED_PRECISION
        , GLOBAL const real4* RESTRICT posqCorrection
#endif
        ) {
    LOCAL int groupConverged;
#ifdef USE_MIXED_PRECISION
    computeCCMAConstraintDirections(constraintAtoms, constraintDistance, atomPositions, posqCorrection);
#else
    computeCCMAConstraintDirections(constraintAtoms, constraintDistance, atomPositions);
#endif
    for (int iteration = 0; iteration < 150; iteration++) {
        SYNC_THREADS
        if (constrainVelocities)
            computeCCMAVelocityConstraintForce(constraintAtoms, constraintDistance, velm, reducedMass,
                    delta1, tol, iteration, &groupConverged);
        else
            computeCCMAPositionConstraintForce(constraintAtoms, constraintDistance, posDelta, reducedMass,
                    delta1, tol, iteration, &groupConverged);
        SYNC_THREADS
        multiplyByCCMAConstraintMatrix(delta1, delta2, constraintMatrixColumn, constraintMatrixValue, iteration);
        SYNC_THREADS
        if (constrainVelocities)
            updateCCMAAtomPositions(atoms, numAtomConstraints, atomConstraints, constraintDistance, velm, velm,
                    delta1, delta2, iteration);
        else
            updateCCMAAtomPositions(atoms, numAtomConstraints, atomConstraints, constraintDistance, posDelta, velm,
                    delta1, delta2, iteration);
        SYNC_THREADS
        if (groupConverged)
            return;
    }
}

/**
 * Compute the positions of virtual sites
 */
KERNEL void computeVirtualSites(GLOBAL real4* RESTRICT posq, GLOBAL real4* RESTRICT posqCorrection,
        GLOBAL const int4* RESTRICT avg2Atoms, GLOBAL const real2* RESTRICT avg2Weights,
        GLOBAL const int4* RESTRICT avg3Atoms, GLOBAL const real4* RESTRICT avg3Weights,
        GLOBAL const int4* RESTRICT outOfPlaneAtoms, GLOBAL const real4* RESTRICT outOfPlaneWeights,
        GLOBAL const int* RESTRICT localCoordsIndex, GLOBAL const int* RESTRICT localCoordsAtoms,
        GLOBAL const real* RESTRICT localCoordsWeights, GLOBAL const real4* RESTRICT localCoordsPos,
        GLOBAL const int* RESTRICT localCoordsStartIndex, GLOBAL const int2* RESTRICT symmetryAtoms,
        GLOBAL const real4* RESTRICT symmetryMatrix, GLOBAL const real4* RESTRICT symmetryOffset,
        GLOBAL const int* RESTRICT symmetryUseBox, real4 periodicBoxVecX, real4 periodicBoxVecY,
        real4 periodicBoxVecZ, real4 recipBoxVecX, real4 recipBoxVecY, real4 recipBoxVecZ,
        GLOBAL const int* RESTRICT vsiteStage, int currentStage) {
    
    // Two particle average sites.
    
    for (int index = GLOBAL_ID; index < NUM_2_AVERAGE; index += GLOBAL_SIZE) {
        int4 atoms = avg2Atoms[index];
#ifdef MULTIPLE_VSITE_STAGES
        if (vsiteStage[atoms.x] != currentStage)
            continue;
#endif
        real2 weights = avg2Weights[index];
        mixed4 pos = loadPos(posq, posqCorrection, atoms.x);
        mixed4 pos1 = loadPos(posq, posqCorrection, atoms.y);
        mixed4 pos2 = loadPos(posq, posqCorrection, atoms.z);
        pos.x = pos1.x*weights.x + pos2.x*weights.y;
        pos.y = pos1.y*weights.x + pos2.y*weights.y;
        pos.z = pos1.z*weights.x + pos2.z*weights.y;
        storePos(posq, posqCorrection, atoms.x, pos);
    }
    
    // Three particle average sites.
    
    for (int index = GLOBAL_ID; index < NUM_3_AVERAGE; index += GLOBAL_SIZE) {
        int4 atoms = avg3Atoms[index];
#ifdef MULTIPLE_VSITE_STAGES
        if (vsiteStage[atoms.x] != currentStage)
            continue;
#endif
        real4 weights = avg3Weights[index];
        mixed4 pos = loadPos(posq, posqCorrection, atoms.x);
        mixed4 pos1 = loadPos(posq, posqCorrection, atoms.y);
        mixed4 pos2 = loadPos(posq, posqCorrection, atoms.z);
        mixed4 pos3 = loadPos(posq, posqCorrection, atoms.w);
        pos.x = pos1.x*weights.x + pos2.x*weights.y + pos3.x*weights.z;
        pos.y = pos1.y*weights.x + pos2.y*weights.y + pos3.y*weights.z;
        pos.z = pos1.z*weights.x + pos2.z*weights.y + pos3.z*weights.z;
        storePos(posq, posqCorrection, atoms.x, pos);
    }
    
    // Out of plane sites.
    
    for (int index = GLOBAL_ID; index < NUM_OUT_OF_PLANE; index += GLOBAL_SIZE) {
        int4 atoms = outOfPlaneAtoms[index];
#ifdef MULTIPLE_VSITE_STAGES
        if (vsiteStage[atoms.x] != currentStage)
            continue;
#endif
        real4 weights = outOfPlaneWeights[index];
        mixed4 pos = loadPos(posq, posqCorrection, atoms.x);
        mixed4 pos1 = loadPos(posq, posqCorrection, atoms.y);
        mixed4 pos2 = loadPos(posq, posqCorrection, atoms.z);
        mixed4 pos3 = loadPos(posq, posqCorrection, atoms.w);
        mixed4 v12 = pos2-pos1;
        mixed4 v13 = pos3-pos1;
        mixed4 cr = cross(v12, v13);
        pos.x = pos1.x + v12.x*weights.x + v13.x*weights.y + cr.x*weights.z;
        pos.y = pos1.y + v12.y*weights.x + v13.y*weights.y + cr.y*weights.z;
        pos.z = pos1.z + v12.z*weights.x + v13.z*weights.y + cr.z*weights.z;
        storePos(posq, posqCorrection, atoms.x, pos);
    }
    
    // Local coordinates sites.
    
    for (int index = GLOBAL_ID; index < NUM_LOCAL_COORDS; index += GLOBAL_SIZE) {
        int siteAtomIndex = localCoordsIndex[index];
#ifdef MULTIPLE_VSITE_STAGES
        if (vsiteStage[siteAtomIndex] != currentStage)
            continue;
#endif
        int start = localCoordsStartIndex[index];
        int end = localCoordsStartIndex[index+1];
        mixed3 origin = make_mixed3(0), xdir = make_mixed3(0), ydir = make_mixed3(0);
        for (int j = start; j < end; j++) {
            mixed3 pos = trimTo3(loadPos(posq, posqCorrection, localCoordsAtoms[j]));
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
        real4 localPosition_4 = localCoordsPos[index];
        mixed3 localPosition = make_mixed3(localPosition_4.x, localPosition_4.y, localPosition_4.z);
        mixed4 pos = loadPos(posq, posqCorrection, siteAtomIndex);
        pos.x = origin.x + xdir.x*localPosition.x + ydir.x*localPosition.y + zdir.x*localPosition.z;
        pos.y = origin.y + xdir.y*localPosition.x + ydir.y*localPosition.y + zdir.y*localPosition.z;
        pos.z = origin.z + xdir.z*localPosition.x + ydir.z*localPosition.y + zdir.z*localPosition.z;
        storePos(posq, posqCorrection, siteAtomIndex, pos);
    }

    // Symmetry sites.

    for (int index = GLOBAL_ID; index < NUM_SYMMETRY; index += GLOBAL_SIZE) {
        int2 atoms = symmetryAtoms[index];
        real4 Rx = symmetryMatrix[3*index];
        real4 Ry = symmetryMatrix[3*index+1];
        real4 Rz = symmetryMatrix[3*index+2];
        real4 v = symmetryOffset[index];
        int useBoxVectors = symmetryUseBox[index];
        mixed4 pos = loadPos(posq, posqCorrection, atoms.x);
        mixed4 r = loadPos(posq, posqCorrection, atoms.y);
        if (useBoxVectors) {
            r.x = r.x*recipBoxVecX.x + r.y*recipBoxVecY.x + r.z*recipBoxVecZ.x;
            r.y = r.y*recipBoxVecY.y + r.z*recipBoxVecZ.y,
            r.z = r.z*recipBoxVecZ.z;
        }
        pos.x = Rx.x*r.x + Rx.y*r.y + Rx.z*r.z + v.x;
        pos.y = Ry.x*r.x + Ry.y*r.y + Ry.z*r.z + v.y;
        pos.z = Rz.x*r.x + Rz.y*r.y + Rz.z*r.z + v.z;
        if (useBoxVectors) {
            pos.x = pos.x*periodicBoxVecX.x + pos.y*periodicBoxVecY.x + pos.z*periodicBoxVecZ.x;
            pos.y = pos.y*periodicBoxVecY.y + pos.z*periodicBoxVecZ.y,
            pos.z = pos.z*periodicBoxVecZ.z;
        }
        storePos(posq, posqCorrection, atoms.x, pos);
    }
}

inline DEVICE real3 loadForce(int index, GLOBAL const mm_long* RESTRICT force) {
    real scale = 1/((real) 0x100000000);
    return make_real3(scale*force[index], scale*force[index+PADDED_NUM_ATOMS], scale*force[index+PADDED_NUM_ATOMS*2]);
}

inline DEVICE void addForce(int index, GLOBAL mm_long* RESTRICT force, real3 value) {
    GLOBAL mm_ulong* f = (GLOBAL mm_ulong*) force;
#ifdef HAS_OVERLAPPING_VSITES
    ATOMIC_ADD(&f[index], (mm_ulong) realToFixedPoint(value.x));
    ATOMIC_ADD(&f[index+PADDED_NUM_ATOMS], (mm_ulong) realToFixedPoint(value.y));
    ATOMIC_ADD(&f[index+PADDED_NUM_ATOMS*2], (mm_ulong) realToFixedPoint(value.z));
#else
    f[index] += (mm_ulong) realToFixedPoint(value.x);
    f[index+PADDED_NUM_ATOMS] += (mm_ulong) realToFixedPoint(value.y);
    f[index+PADDED_NUM_ATOMS*2] += (mm_ulong) realToFixedPoint(value.z);
#endif
}

/**
 * Distribute forces from virtual sites to the atoms they are based on.
 */
KERNEL void distributeVirtualSiteForces(GLOBAL const real4* RESTRICT posq, GLOBAL const real4* RESTRICT posqCorrection, GLOBAL mm_long* RESTRICT force,
        GLOBAL const int4* RESTRICT avg2Atoms, GLOBAL const real2* RESTRICT avg2Weights,
        GLOBAL const int4* RESTRICT avg3Atoms, GLOBAL const real4* RESTRICT avg3Weights,
        GLOBAL const int4* RESTRICT outOfPlaneAtoms, GLOBAL const real4* RESTRICT outOfPlaneWeights,
        GLOBAL const int* RESTRICT localCoordsIndex, GLOBAL const int* RESTRICT localCoordsAtoms,
        GLOBAL const real* RESTRICT localCoordsWeights, GLOBAL const real4* RESTRICT localCoordsPos,
        GLOBAL const int* RESTRICT localCoordsStartIndex, GLOBAL const int2* RESTRICT symmetryAtoms,
        GLOBAL const real4* RESTRICT symmetryMatrix, GLOBAL const real4* RESTRICT symmetryOffset,
        GLOBAL const int* RESTRICT symmetryUseBox, real4 periodicBoxVecX, real4 periodicBoxVecY,
        real4 periodicBoxVecZ, real4 recipBoxVecX, real4 recipBoxVecY, real4 recipBoxVecZ,
        GLOBAL const int* RESTRICT vsiteStage, int currentStage) {
    
    // Two particle average sites.
    
    for (int index = GLOBAL_ID; index < NUM_2_AVERAGE; index += GLOBAL_SIZE) {
        int4 atoms = avg2Atoms[index];
#ifdef MULTIPLE_VSITE_STAGES
        if (vsiteStage[atoms.x] != currentStage)
            continue;
#endif
        real2 weights = avg2Weights[index];
        real3 f = loadForce(atoms.x, force);
        addForce(atoms.y, force, f*weights.x);
        addForce(atoms.z, force, f*weights.y);
    }
    
    // Three particle average sites.
    
    for (int index = GLOBAL_ID; index < NUM_3_AVERAGE; index += GLOBAL_SIZE) {
        int4 atoms = avg3Atoms[index];
#ifdef MULTIPLE_VSITE_STAGES
        if (vsiteStage[atoms.x] != currentStage)
            continue;
#endif
        real4 weights = avg3Weights[index];
        real3 f = loadForce(atoms.x, force);
        addForce(atoms.y, force, f*weights.x);
        addForce(atoms.z, force, f*weights.y);
        addForce(atoms.w, force, f*weights.z);
    }
    
    // Out of plane sites.
    
    for (int index = GLOBAL_ID; index < NUM_OUT_OF_PLANE; index += GLOBAL_SIZE) {
        int4 atoms = outOfPlaneAtoms[index];
#ifdef MULTIPLE_VSITE_STAGES
        if (vsiteStage[atoms.x] != currentStage)
            continue;
#endif
        real4 weights = outOfPlaneWeights[index];
        mixed4 pos1 = loadPos(posq, posqCorrection, atoms.y);
        mixed4 pos2 = loadPos(posq, posqCorrection, atoms.z);
        mixed4 pos3 = loadPos(posq, posqCorrection, atoms.w);
        mixed4 v12 = pos2-pos1;
        mixed4 v13 = pos3-pos1;
        real3 f = loadForce(atoms.x, force);
        real3 fp2 = make_real3((real) (weights.x*f.x - weights.z*v13.z*f.y + weights.z*v13.y*f.z),
                   (real) (weights.z*v13.z*f.x + weights.x*f.y - weights.z*v13.x*f.z),
                   (real) (-weights.z*v13.y*f.x + weights.z*v13.x*f.y + weights.x*f.z));
        real3 fp3 = make_real3((real) (weights.y*f.x + weights.z*v12.z*f.y - weights.z*v12.y*f.z),
                   (real) (-weights.z*v12.z*f.x + weights.y*f.y + weights.z*v12.x*f.z),
                   (real) (weights.z*v12.y*f.x - weights.z*v12.x*f.y + weights.y*f.z));
        addForce(atoms.y, force, f-fp2-fp3);
        addForce(atoms.z, force, fp2);
        addForce(atoms.w, force, fp3);
    }
    
    // Local coordinates sites.
    
    for (int index = GLOBAL_ID; index < NUM_LOCAL_COORDS; index += GLOBAL_SIZE) {
        int siteAtomIndex = localCoordsIndex[index];
#ifdef MULTIPLE_VSITE_STAGES
        if (vsiteStage[siteAtomIndex] != currentStage)
            continue;
#endif
        int start = localCoordsStartIndex[index];
        int end = localCoordsStartIndex[index+1];
        mixed3 origin = make_mixed3(0), xdir = make_mixed3(0), ydir = make_mixed3(0);
        for (int j = start; j < end; j++) {
            mixed3 pos = trimTo3(loadPos(posq, posqCorrection, localCoordsAtoms[j]));
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
        real4 localPosition_4 = localCoordsPos[index];
        mixed3 localPosition = make_mixed3(localPosition_4.x, localPosition_4.y, localPosition_4.z);

        // The derivatives for this case are very complicated.  They were computed with SymPy then simplified by hand.

        real3 f = loadForce(siteAtomIndex, force);
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
            real3 fresult = make_real3(0);
            fresult.x += fp1.x*wxScaled*(1-dx.x*dx.x) + fp1.z*(dz.x*sx   ) + fp1.y*((-dx.x*dy.x     )*wxScaled + dy.x*sx - dx.y*t2 - dx.z*t3) + f.x*originWeight;
            fresult.y += fp1.x*wxScaled*( -dx.x*dx.y) + fp1.z*(dz.x*sy+t3) + fp1.y*((-dx.y*dy.x-dz.z)*wxScaled + dy.x*sy + dx.y*t1);
            fresult.z += fp1.x*wxScaled*( -dx.x*dx.z) + fp1.z*(dz.x*sz-t2) + fp1.y*((-dx.z*dy.x+dz.y)*wxScaled + dy.x*sz + dx.z*t1);
            fresult.x += fp2.x*wxScaled*( -dx.y*dx.x) + fp2.z*(dz.y*sx-t3) - fp2.y*(( dx.x*dy.y-dz.z)*wxScaled - dy.y*sx - dx.x*t2);
            fresult.y += fp2.x*wxScaled*(1-dx.y*dx.y) + fp2.z*(dz.y*sy   ) - fp2.y*(( dx.y*dy.y     )*wxScaled - dy.y*sy + dx.x*t1 + dx.z*t3) + f.y*originWeight;
            fresult.z += fp2.x*wxScaled*( -dx.y*dx.z) + fp2.z*(dz.y*sz+t1) - fp2.y*(( dx.z*dy.y+dz.x)*wxScaled - dy.y*sz - dx.z*t2);
            fresult.x += fp3.x*wxScaled*( -dx.z*dx.x) + fp3.z*(dz.z*sx+t2) + fp3.y*((-dx.x*dy.z-dz.y)*wxScaled + dy.z*sx + dx.x*t3);
            fresult.y += fp3.x*wxScaled*( -dx.z*dx.y) + fp3.z*(dz.z*sy-t1) + fp3.y*((-dx.y*dy.z+dz.x)*wxScaled + dy.z*sy + dx.y*t3);
            fresult.z += fp3.x*wxScaled*(1-dx.z*dx.z) + fp3.z*(dz.z*sz   ) + fp3.y*((-dx.z*dy.z     )*wxScaled + dy.z*sz - dx.x*t1 - dx.y*t2) + f.z*originWeight;
            addForce(localCoordsAtoms[j], force, fresult);
        }
    }

    // Symmetry sites.

    for (int index = GLOBAL_ID; index < NUM_SYMMETRY; index += GLOBAL_SIZE) {
        int2 atoms = symmetryAtoms[index];
        real4 Rx = symmetryMatrix[3*index];
        real4 Ry = symmetryMatrix[3*index+1];
        real4 Rz = symmetryMatrix[3*index+2];
        real4 v = symmetryOffset[index];
        int useBoxVectors = symmetryUseBox[index];
        real3 f = loadForce(atoms.x, force);
        if (useBoxVectors) {
            f.x = f.x*periodicBoxVecX.x + f.y*periodicBoxVecY.x + f.z*periodicBoxVecZ.x;
            f.y = f.y*periodicBoxVecY.y + f.z*periodicBoxVecZ.y,
            f.z = f.z*periodicBoxVecZ.z;
        }
        real3 f1 = make_real3(Rx.x*f.x + Ry.x*f.y + Rz.x*f.z,
                              Rx.y*f.x + Ry.y*f.y + Rz.y*f.z,
                              Rx.z*f.x + Ry.z*f.y + Rz.z*f.z);
        if (useBoxVectors) {
            f1.x = f1.x*recipBoxVecX.x + f1.y*recipBoxVecY.x + f1.z*recipBoxVecZ.x;
            f1.y = f1.y*recipBoxVecY.y + f1.z*recipBoxVecZ.y,
            f1.z = f1.z*recipBoxVecZ.z;
        }
        addForce(atoms.y, force, f1);
    }
}

/**
 * Copy the distributed forces from the long buffer back to the float buffer.
 */
KERNEL void saveDistributedForces(GLOBAL const mm_long* RESTRICT longForces, GLOBAL real4* RESTRICT forces) {
    for (int index = GLOBAL_ID; index < NUM_ATOMS; index += GLOBAL_SIZE) {
        real3 f = loadForce(index, longForces);
        forces[index] = make_real4(f.x, f.y, f.z, 0);
    }
}

/**
 * Apply a time shift to the velocities before computing kinetic energy.
 */
KERNEL void timeShiftVelocities(GLOBAL mixed4* RESTRICT velm, GLOBAL const mm_long* RESTRICT force, real timeShift) {
    const mixed scale = timeShift/(mixed) 0x100000000;
    for (int index = GLOBAL_ID; index < NUM_ATOMS; index += GLOBAL_SIZE) {
        mixed4 velocity = velm[index];
        if (velocity.w != 0.0) {
            velocity.x += scale*force[index]*velocity.w;
            velocity.y += scale*force[index+PADDED_NUM_ATOMS]*velocity.w;
            velocity.z += scale*force[index+PADDED_NUM_ATOMS*2]*velocity.w;
            velm[index] = velocity;
        }
    }
}

/**
 * Compute the total kinetic energy.
 */
KERNEL void computeKineticEnergy(GLOBAL mixed4* RESTRICT velm, GLOBAL mixed* result) {
    LOCAL mixed tempBuffer[KE_WORK_GROUP_SIZE];
    mixed sum = 0;
    for (unsigned int index = LOCAL_ID; index < NUM_ATOMS; index += LOCAL_SIZE) {
        mixed4 v = velm[index];
        if (v.w != 0)
            sum += (v.x*v.x+v.y*v.y+v.z*v.z)/v.w;
    }
    tempBuffer[LOCAL_ID] = sum;
    for (int i = 1; i < KE_WORK_GROUP_SIZE; i *= 2) {
        SYNC_THREADS;
        if (LOCAL_ID%(i*2) == 0 && LOCAL_ID+i < KE_WORK_GROUP_SIZE)
            tempBuffer[LOCAL_ID] += tempBuffer[LOCAL_ID+i];
    }
    if (LOCAL_ID == 0)
        *result = 0.5f*tempBuffer[0];
}
