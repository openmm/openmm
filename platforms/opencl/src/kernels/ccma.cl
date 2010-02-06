/**
 * Compute the direction each constraint is pointing in.  This is called once at the beginning of constraint evaluation.
 */
__kernel void computeConstraintDirections(__global int2* constraintAtoms, __global float4* constraintDistance, __global float4* atomPositions,
        __global int* converged) {
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

    // Mark that no work groups have converged yet.

    for (int index = get_global_id(0); index < get_num_groups(0); index += get_global_size(0))
        converged[index] = false;
}

/**
 * Compute the force applied by each constraint.
 */
__kernel void computeConstraintForce(__global int2* constraintAtoms, __global float4* constraintDistance, __global float4* atomPositions,
        __global float* reducedMass, __global float* delta1, __global int* converged, __local int* convergedBuffer, float tol) {
    if (converged[get_group_id(0)])
        return; // The constraint iteration has already converged.
    float lowerTol = 1.0f-2.0f*tol+tol*tol;
    float upperTol = 1.0f+2.0f*tol+tol*tol;
    int threadConverged = 1;
    for (int index = get_global_id(0); index < NUM_CONSTRAINTS; index += get_global_size(0)) {
        // Compute the force due to this constraint.

        int2 atoms = constraintAtoms[index];
        float4 dir = constraintDistance[index];
        float4 rp_ij = atomPositions[atoms.x]-atomPositions[atoms.y];
        rp_ij.xyz += dir.xyz;
        float rp2 = rp_ij.x*rp_ij.x + rp_ij.y*rp_ij.y + rp_ij.z*rp_ij.z;
        float dist2 = dir.w*dir.w;
        float diff = dist2 - rp2;
        float rrpr = rp_ij.x*dir.x + rp_ij.y*dir.y + rp_ij.z*dir.z;
        float d_ij2 = dir.x*dir.x + dir.y*dir.y + dir.z*dir.z;
        delta1[index] = (rrpr > d_ij2*1e-6f ? reducedMass[index]*diff/rrpr : 0.0f);

        // See whether it has converged.

        threadConverged &= (rp2 >= lowerTol*dist2 && rp2 <= upperTol*dist2);
    }

    // Perform a parallel reduction to see if all constraints handled by this work group have converged.

    convergedBuffer[get_local_id(0)] = threadConverged;
    barrier(CLK_LOCAL_MEM_FENCE);
    for (int step = 1; step < get_local_size(0); step *= 2) {
        if (get_local_id(0)%(2*step) == 0)
            convergedBuffer[get_local_id(0)] &= convergedBuffer[get_local_id(0)+step];
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    if (get_local_id(0) == 0)
        converged[get_group_id(0)] = convergedBuffer[0];
}

/**
 * Multiply the vector of constraint forces by the constraint matrix.
 */
__kernel void multiplyByConstraintMatrix(__global float* delta1, __global float* delta2, __global int* constraintMatrixColumn,
        __global float* constraintMatrixValue, __global int* converged, __local int* convergedBuffer) {
    // First see whether all work groups have converged.

    convergedBuffer[get_local_id(0)] = true;
    for (int index = get_local_id(0); index < get_num_groups(0); index += get_local_size(0))
        convergedBuffer[get_local_id(0)] &= converged[index];
    barrier(CLK_LOCAL_MEM_FENCE);
    for (int step = 1; step < get_local_size(0); step *= 2) {
        if (get_local_id(0)%(2*step) == 0)
            convergedBuffer[get_local_id(0)] &= convergedBuffer[get_local_id(0)+step];
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    if (get_local_id(0) == 0)
        converged[get_group_id(0)] = convergedBuffer[0];
    if (converged[0])
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
__kernel void updateAtomPositions(__global int* numAtomConstraints, __global int* atomConstraints, __global float4* constraintDistance,
        __global float4* atomPositions, __global float4* velm, __global float* delta1, __global float* delta2, __global int* converged, int iteration) {
    if (converged[get_group_id(0)])
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
