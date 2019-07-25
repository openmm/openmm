/**
 * Scale the particle positions with each axis independent.
 */

__kernel void scalePositions(float scaleX, float scaleY, float scaleZ, int numMolecules, real4 periodicBoxSize, real4 invPeriodicBoxSize,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, __global real4* restrict posq,
        __global const int* restrict moleculeAtoms, __global const int* restrict moleculeStartIndex) {
    for (int index = get_global_id(0); index < numMolecules; index += get_global_size(0)) {
        int first = moleculeStartIndex[index];
        int last = moleculeStartIndex[index+1];
        int numAtoms = last-first;

        // Find the center of each molecule.

        real3 center = (real3) 0;
        for (int atom = first; atom < last; atom++)
            center += posq[moleculeAtoms[atom]].xyz;
        center /= (real) numAtoms;

        // Move it into the first periodic box.

        real3 oldCenter = center;
        APPLY_PERIODIC_TO_POS(center)
        real3 delta = oldCenter-center;;
        real3 scaleXYZ = (real3) (scaleX, scaleY, scaleZ);

        // Now scale the position of the molecule center.

        delta = center*(scaleXYZ-1)-delta;
        for (int atom = first; atom < last; atom++) {
            real4 pos = posq[moleculeAtoms[atom]];
            pos.xyz += delta.xyz;
            posq[moleculeAtoms[atom]] = pos;
        }
    }
}
