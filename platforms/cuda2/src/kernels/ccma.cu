/**
 * Compute the direction each constraint is pointing in.  This is called once at the beginning of constraint evaluation.
 */
extern "C" __global__ void computeConstraintDirections(const int2* __restrict__ constraintAtoms, real4* __restrict__ constraintDistance, const real4* __restrict__ atomPositions) {
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_CONSTRAINTS; index += blockDim.x*gridDim.x) {
        // Compute the direction for this constraint.

        int2 atoms = constraintAtoms[index];
        real4 dir = constraintDistance[index];
        real4 oldPos1 = atomPositions[atoms.x];
        real4 oldPos2 = atomPositions[atoms.y];
        dir.x = oldPos1.x-oldPos2.x;
        dir.y = oldPos1.y-oldPos2.y;
        dir.z = oldPos1.z-oldPos2.z;
        constraintDistance[index] = dir;
    }
}

/**
 * Compute the force applied by each constraint.
 */
extern "C" __global__ void computeConstraintForce(const int2* __restrict__ constraintAtoms, const real4* __restrict__ constraintDistance, const real4* __restrict__ atomPositions,
        const real* __restrict__ reducedMass, real* __restrict__ delta1, int* __restrict__ converged, float tol, int iteration) {
    __shared__ int groupConverged;
    if (converged[1-iteration%2]) {
        if (blockIdx.x == 0 && threadIdx.x == 0)
            converged[iteration%2] = 1;
        return; // The constraint iteration has already converged.
    }
    if (threadIdx.x == 0)
        groupConverged = 1;
    __syncthreads();
    real lowerTol = 1.0f-2.0f*tol+tol*tol;
    real upperTol = 1.0f+2.0f*tol+tol*tol;
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_CONSTRAINTS; index += blockDim.x*gridDim.x) {
        // Compute the force due to this constraint.

        int2 atoms = constraintAtoms[index];
        real4 dir = constraintDistance[index];
        real4 rp_ij = atomPositions[atoms.x]-atomPositions[atoms.y];
#ifndef CONSTRAIN_VELOCITIES
        rp_ij.x += dir.x;
        rp_ij.y += dir.y;
        rp_ij.z += dir.z;
#endif
        real rrpr = rp_ij.x*dir.x + rp_ij.y*dir.y + rp_ij.z*dir.z;
        real d_ij2 = dir.x*dir.x + dir.y*dir.y + dir.z*dir.z;
#ifdef CONSTRAIN_VELOCITIES
        delta1[index] = -2.0f*reducedMass[index]*rrpr/d_ij2;

        // See whether it has converged.

        if (groupConverged && fabs(delta1[index]) > tol) {
            groupConverged = 0;
            converged[iteration%2] = 0;
        }
#else
        real rp2 = rp_ij.x*rp_ij.x + rp_ij.y*rp_ij.y + rp_ij.z*rp_ij.z;
        real dist2 = dir.w*dir.w;
        real diff = dist2 - rp2;
        delta1[index] = (rrpr > d_ij2*1e-6f ? reducedMass[index]*diff/rrpr : 0.0f);

        // See whether it has converged.

        if (groupConverged && (rp2 < lowerTol*dist2 || rp2 > upperTol*dist2)) {
            groupConverged = 0;
            converged[iteration%2] = 0;
        }
#endif
    }
}

/**
 * Multiply the vector of constraint forces by the constraint matrix.
 */
extern "C" __global__ void multiplyByConstraintMatrix(const real* __restrict__ delta1, real* __restrict__ delta2, const int* __restrict__ constraintMatrixColumn,
        const real* __restrict__ constraintMatrixValue, const int* __restrict__ converged, int iteration) {
    if (converged[iteration%2])
        return; // The constraint iteration has already converged.

    // Multiply by the inverse constraint matrix.

    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_CONSTRAINTS; index += blockDim.x*gridDim.x) {
        real sum = 0.0f;
        for (int i = 0; ; i++) {
            int element = index+i*NUM_CONSTRAINTS;
            int column = constraintMatrixColumn[element];
            if (column >= NUM_CONSTRAINTS)
                break;
            sum += delta1[column]*constraintMatrixValue[element];
        }
        delta2[index] = sum;
    }
}

/**
 * Update the atom positions based on constraint forces.
 */
extern "C" __global__ void updateAtomPositions(const int* __restrict__ numAtomConstraints, const int* __restrict__ atomConstraints, const real4* __restrict__ constraintDistance,
        real4* __restrict__ atomPositions, const real4* __restrict__ velm, const real* __restrict__ delta1, const real* __restrict__ delta2, int* __restrict__ converged, int iteration) {
    if (blockIdx.x == 0 && threadIdx.x == 0)
        converged[1-iteration%2] = 1;
    if (converged[iteration%2])
        return; // The constraint iteration has already converged.
    real damping = (iteration < 2 ? 0.5f : 1.0f);
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_ATOMS; index += blockDim.x*gridDim.x) {
        // Compute the new position of this atom.

        real4 atomPos = atomPositions[index];
        real invMass = velm[index].w;
        int num = numAtomConstraints[index];
        for (int i = 0; i < num; i++) {
            int constraint = atomConstraints[index+i*NUM_ATOMS];
            bool forward = (constraint > 0);
            constraint = (forward ? constraint-1 : -constraint-1);
            real constraintForce = damping*invMass*delta2[constraint];
            constraintForce = (forward ? constraintForce : -constraintForce);
            real4 dir = constraintDistance[constraint];
            atomPos.x += constraintForce*dir.x;
            atomPos.y += constraintForce*dir.y;
            atomPos.z += constraintForce*dir.z;
        }
        atomPositions[index] = atomPos;
    }
}
