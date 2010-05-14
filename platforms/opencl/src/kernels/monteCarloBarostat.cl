/**
 * Scale the particle positions.
 */

__kernel void scalePositions(float scale, int numMolecules, float4 periodicBoxSize, float4 invPeriodicBoxSize, __global float4* posq,
        __global int* moleculeAtoms, __global int* moleculeStartIndex) {
    for (int index = get_global_id(0); index < numMolecules; index += get_global_size(0)) {
        int first = moleculeStartIndex[index];
        int last = moleculeStartIndex[index+1];
        int numAtoms = last-first;

        // Find the center of each molecule.

        float4 center = (float4) 0.0f;
        for (int atom = first; atom < last; atom++)
            center += posq[moleculeAtoms[atom]];
        center /= (float) numAtoms;

        // Move it into the first periodic box.

        int xcell = (int) floor(center.x*invPeriodicBoxSize.x);
        int ycell = (int) floor(center.y*invPeriodicBoxSize.y);
        int zcell = (int) floor(center.z*invPeriodicBoxSize.z);
        float4 delta = (float) (xcell*periodicBoxSize.x, ycell*periodicBoxSize.y, zcell*periodicBoxSize.z, 0);
        center -= delta;

        // Now scale the position of the molecule center.

        delta = center*(scale-1)-delta;
        for (int atom = first; atom < last; atom++) {
            float4 pos = posq[moleculeAtoms[atom]];
            pos.xyz += delta.xyz;
            posq[moleculeAtoms[atom]] = pos;
        }
    }
}
