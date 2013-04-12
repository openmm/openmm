/**
 * Apply a time shift to the velocities before computing kinetic energy.
 */
__kernel void timeShiftVelocities(__global mixed4* restrict velm, __global const real4* restrict force, real timeShift) {
    for (int index = get_global_id(0); index < NUM_ATOMS; index += get_global_size(0)) {
        mixed4 velocity = velm[index];
        if (velocity.w != 0.0) {
            mixed4 f = convert_mixed4(force[index]);
            velocity.xyz += timeShift*f.xyz*velocity.w;
            velm[index] = velocity;
        }
    }
}