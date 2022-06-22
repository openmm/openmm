/**
 * Scale the particle positions with each axis independent.
 */

KERNEL void scalePositions(float scaleX, float scaleY, float scaleZ, int numMolecules, real4 periodicBoxSize,
        real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, GLOBAL real4* RESTRICT posq,
        GLOBAL const int* RESTRICT moleculeAtoms, GLOBAL const int* RESTRICT moleculeStartIndex) {
    for (int index = GLOBAL_ID; index < numMolecules; index += GLOBAL_SIZE) {
        int first = moleculeStartIndex[index];
        int last = moleculeStartIndex[index+1];
        int numAtoms = last-first;

        // Find the center of each molecule.

        real3 center = make_real3(0, 0, 0);
        for (int atom = first; atom < last; atom++) {
            real4 pos = posq[moleculeAtoms[atom]];
            center.x += pos.x;
            center.y += pos.y;
            center.z += pos.z;
        }
        real invNumAtoms = RECIP((real) numAtoms);
        center.x *= invNumAtoms;
        center.y *= invNumAtoms;
        center.z *= invNumAtoms;

        // Now scale the position of the molecule center.

        real3 delta;
        delta.x = center.x*(scaleX-1);
        delta.y = center.y*(scaleY-1);
        delta.z = center.z*(scaleZ-1);
        for (int atom = first; atom < last; atom++) {
            real4 pos = posq[moleculeAtoms[atom]];
            pos.x += delta.x;
            pos.y += delta.y;
            pos.z += delta.z;
            posq[moleculeAtoms[atom]] = pos;
        }
    }
}
