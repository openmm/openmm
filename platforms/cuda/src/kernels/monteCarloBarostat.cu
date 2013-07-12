/**
 * Scale the particle positions with each axis independent
 */

extern "C" __global__ void scalePositions(float scaleX, float scaleY, float scaleZ, int numMolecules, real4 periodicBoxSize, real4 invPeriodicBoxSize, real4* __restrict__ posq,
        const int* __restrict__ moleculeAtoms, const int* __restrict__ moleculeStartIndex) {
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < numMolecules; index += blockDim.x*gridDim.x) {
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
        real invNumAtoms = RECIP(numAtoms);
        center.x *= invNumAtoms;
        center.y *= invNumAtoms;
        center.z *= invNumAtoms;

        // Move it into the first periodic box.

        int xcell = (int) floor(center.x*invPeriodicBoxSize.x);
        int ycell = (int) floor(center.y*invPeriodicBoxSize.y);
        int zcell = (int) floor(center.z*invPeriodicBoxSize.z);
        real3 delta = make_real3(xcell*periodicBoxSize.x, ycell*periodicBoxSize.y, zcell*periodicBoxSize.z);
        center.x -= delta.x;
        center.y -= delta.y;
        center.z -= delta.z;

        // Now scale the position of the molecule center.

        delta.x = center.x*(scaleX-1)-delta.x;
        delta.y = center.y*(scaleY-1)-delta.y;
        delta.z = center.z*(scaleZ-1)-delta.z;
        for (int atom = first; atom < last; atom++) {
            real4 pos = posq[moleculeAtoms[atom]];
            pos.x += delta.x;
            pos.y += delta.y;
            pos.z += delta.z;
            posq[moleculeAtoms[atom]] = pos;
        }
    }
}
