/**
 * Compute the direction each constraint is pointing in.  This is called once at the beginning of constraint evaluation.
 */
__kernel void computeConstraintDirections(__global const int2* restrict constraintAtoms, __global float4* restrict constraintDistance, __global const float4* restrict atomPositions) {
    for (int index = get_global_id(0); index < NUM_CONSTRAINTS; index += get_global_size(0)) {
        // Compute the direction for this constraint.

        int2 atoms = constraintAtoms[index];
        float4 dir = constraintDistance[index];
        float4 oldPos1 = atomPositions[atoms.x];
        float4 oldPos2 = atomPositions[atoms.y];
        dir.x = oldPos1.x-oldPos2.x;
        dir.y = oldPos1.y-oldPos2.y;
        dir.z = oldPos1.z-oldPos2.z;
        constraintDistance[index] = dir;
    }
}

/**
 * Compute the force applied by each constraint.
 */
__kernel void computeConstraintForce(__global const int2* restrict constraintAtoms, __global const float4* restrict constraintDistance, __global const float4* restrict atomPositions,
        __global const float* restrict reducedMass, __global float* restrict delta1, __global int* restrict converged, float tol, int iteration) {
    __local int groupConverged;
    if (converged[1-iteration%2]) {
        if (get_global_id(0) == 0)
            converged[iteration%2] = 1;
        return; // The constraint iteration has already converged.
    }
    if (get_local_id(0) == 0)
        groupConverged = 1;
    barrier(CLK_LOCAL_MEM_FENCE);
    float lowerTol = 1.0f-2.0f*tol+tol*tol;
    float upperTol = 1.0f+2.0f*tol+tol*tol;
    for (int index = get_global_id(0); index < NUM_CONSTRAINTS; index += get_global_size(0)) {
        // Compute the force due to this constraint.

        int2 atoms = constraintAtoms[index];
        float4 dir = constraintDistance[index];
        float4 rp_ij = atomPositions[atoms.x]-atomPositions[atoms.y];
#ifndef CONSTRAIN_VELOCITIES
        rp_ij.xyz += dir.xyz;
#endif
        float rrpr = rp_ij.x*dir.x + rp_ij.y*dir.y + rp_ij.z*dir.z;
        float d_ij2 = dir.x*dir.x + dir.y*dir.y + dir.z*dir.z;
#ifdef CONSTRAIN_VELOCITIES
        delta1[index] = -2.0f*reducedMass[index]*rrpr/d_ij2;

        // See whether it has converged.

        if (groupConverged && fabs(delta1[index]) > tol) {
            groupConverged = 0;
            converged[iteration%2] = 0;
        }
#else
        float rp2 = rp_ij.x*rp_ij.x + rp_ij.y*rp_ij.y + rp_ij.z*rp_ij.z;
        float dist2 = dir.w*dir.w;
        float diff = dist2 - rp2;
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
__kernel void multiplyByConstraintMatrix(__global const float* restrict delta1, __global float* restrict delta2, __global const int* restrict constraintMatrixColumn,
        __global const float* restrict constraintMatrixValue, __global const int* restrict converged, int iteration) {
    if (converged[iteration%2])
        return; // The constraint iteration has already converged.

    // Multiply by the inverse constraint matrix.

    for (int index = get_global_id(0); index < NUM_CONSTRAINTS; index += get_global_size(0)) {
        float sum = 0.0f;
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
__kernel void updateAtomPositions(__global const int* restrict numAtomConstraints, __global const int* restrict atomConstraints, __global const float4* restrict constraintDistance,
        __global float4* restrict atomPositions, __global const float4* restrict velm, __global const float* restrict delta1, __global const float* restrict delta2, __global int* restrict converged, int iteration) {
    if (get_global_id(0) == 0)
        converged[1-iteration%2] = 1;
    if (converged[iteration%2])
        return; // The constraint iteration has already converged.
    float damping = (iteration < 2 ? 0.5f : 1.0f);
    for (int index = get_global_id(0); index < NUM_ATOMS; index += get_global_size(0)) {
        // Compute the new position of this atom.

        float4 atomPos = atomPositions[index];
        float invMass = velm[index].w;
        int num = numAtomConstraints[index];
        for (int i = 0; i < num; i++) {
            int constraint = atomConstraints[index+i*NUM_ATOMS];
            bool forward = (constraint > 0);
            constraint = (forward ? constraint-1 : -constraint-1);
            float constraintForce = damping*invMass*delta2[constraint];
            constraintForce = (forward ? constraintForce : -constraintForce);
            float4 dir = constraintDistance[constraint];
            atomPos.x += constraintForce*dir.x;
            atomPos.y += constraintForce*dir.y;
            atomPos.z += constraintForce*dir.z;
        }
        atomPositions[index] = atomPos;
    }
}
