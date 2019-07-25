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
 * Compute the direction each constraint is pointing in.  This is called once at the beginning of constraint evaluation.
 */
__kernel void computeConstraintDirections(__global const int2* restrict constraintAtoms, __global mixed4* restrict constraintDistance,
        __global const real4* restrict atomPositions, __global const real4* restrict posCorrection, __global int* restrict converged) {
    for (int index = get_global_id(0); index < NUM_CONSTRAINTS; index += get_global_size(0)) {
        // Compute the direction for this constraint.

        int2 atoms = constraintAtoms[index];
        mixed4 dir = constraintDistance[index];
        mixed4 oldPos1 = loadPos(atomPositions, posCorrection, atoms.x);
        mixed4 oldPos2 = loadPos(atomPositions, posCorrection, atoms.y);
        dir.x = oldPos1.x-oldPos2.x;
        dir.y = oldPos1.y-oldPos2.y;
        dir.z = oldPos1.z-oldPos2.z;
        constraintDistance[index] = dir;
    }
    if (get_global_id(0) == 0) {
        converged[0] = 1;
        converged[1] = 0;
    }
}

/**
 * Compute the force applied by each constraint.
 */
__kernel void computeConstraintForce(__global const int2* restrict constraintAtoms, __global const mixed4* restrict constraintDistance, __global const mixed4* restrict atomPositions,
        __global const mixed* restrict reducedMass, __global mixed* restrict delta1, __global int* restrict converged, __global int* restrict hostConvergedFlag, mixed tol, int iteration) {
    __local int groupConverged;
    if (converged[1-iteration%2]) {
        if (get_global_id(0) == 0) {
            converged[iteration%2] = 1;
            hostConvergedFlag[0] = 1;
        }
        return; // The constraint iteration has already converged.
    }
    if (get_local_id(0) == 0)
        groupConverged = 1;
    barrier(CLK_LOCAL_MEM_FENCE);
    mixed lowerTol = 1-2*tol+tol*tol;
    mixed upperTol = 1+2*tol+tol*tol;
    for (int index = get_global_id(0); index < NUM_CONSTRAINTS; index += get_global_size(0)) {
        // Compute the force due to this constraint.

        int2 atoms = constraintAtoms[index];
        mixed4 dir = constraintDistance[index];
        mixed4 rp_ij = atomPositions[atoms.x]-atomPositions[atoms.y];
#ifndef CONSTRAIN_VELOCITIES
        rp_ij.xyz += dir.xyz;
#endif
        mixed rrpr = rp_ij.x*dir.x + rp_ij.y*dir.y + rp_ij.z*dir.z;
        mixed d_ij2 = dir.x*dir.x + dir.y*dir.y + dir.z*dir.z;
#ifdef CONSTRAIN_VELOCITIES
        delta1[index] = -2*reducedMass[index]*rrpr/d_ij2;

        // See whether it has converged.

        if (groupConverged && fabs(delta1[index]) > tol) {
            groupConverged = 0;
            converged[iteration%2] = 0;
        }
#else
        mixed rp2 = rp_ij.x*rp_ij.x + rp_ij.y*rp_ij.y + rp_ij.z*rp_ij.z;
        mixed dist2 = dir.w*dir.w;
        mixed diff = dist2 - rp2;
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
__kernel void multiplyByConstraintMatrix(__global const mixed* restrict delta1, __global mixed* restrict delta2, __global const int* restrict constraintMatrixColumn,
        __global const mixed* restrict constraintMatrixValue, __global const int* restrict converged, int iteration) {
    if (converged[iteration%2])
        return; // The constraint iteration has already converged.

    // Multiply by the inverse constraint matrix.

    for (int index = get_global_id(0); index < NUM_CONSTRAINTS; index += get_global_size(0)) {
        mixed sum = 0;
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
__kernel void updateAtomPositions(__global const int* restrict numAtomConstraints, __global const int* restrict atomConstraints, __global const mixed4* restrict constraintDistance,
        __global mixed4* restrict atomPositions, __global const mixed4* restrict velm, __global const mixed* restrict delta1, __global const mixed* restrict delta2, __global int* restrict converged, int iteration) {
    if (get_global_id(0) == 0)
        converged[1-iteration%2] = 1;
    if (converged[iteration%2])
        return; // The constraint iteration has already converged.
    mixed damping = (iteration < 2 ? 0.5f : 1.0f);
    for (int index = get_global_id(0); index < NUM_ATOMS; index += get_global_size(0)) {
        // Compute the new position of this atom.

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
