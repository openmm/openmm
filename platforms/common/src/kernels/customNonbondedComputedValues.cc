/**
 * Calculate per-particle computed values for a CustomNonbondedForce.
 */

KERNEL void computePerParticleValues(PARAMETER_ARGUMENTS) {
    for (int index = GLOBAL_ID; index < NUM_ATOMS; index += GLOBAL_SIZE) {
        COMPUTE_VALUES
    }
}
