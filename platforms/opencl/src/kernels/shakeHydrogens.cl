/**
 * Enforce constraints on SHAKE clusters
 */

__kernel void applyShakeToHydrogens(int numClusters, float tol, __global float4* oldPos, __global float4* posDelta, __global float4* newDelta, __global int4* clusterAtoms, __global float4* clusterParams) {
    int index = get_global_id(0);
    while (index < numClusters) {
        // Load the data for this cluster.

        int4 atoms = clusterAtoms[index];
        float4 params = clusterParams[index];
        float4 pos = oldPos[atoms.x];
        float4 xpi = posDelta[atoms.x];
        float4 pos1 = oldPos[atoms.y];
        float4 xpj1 = posDelta[atoms.y];
        float4 pos2 = {0.0f, 0.0f, 0.0f, 0.0f};
        float4 xpj2 = {0.0f, 0.0f, 0.0f, 0.0f};
        float invMassCentral = params.x;
        float avgMass = params.y;
        float d2 = params.z;
        float invMassPeripheral = params.w;
        if (atoms.z != -1) {
            pos2 = oldPos[atoms.z];
            xpj2 = posDelta[atoms.z];
        }
        float4 pos3 = {0.0f, 0.0f, 0.0f, 0.0f};
        float4 xpj3 = {0.0f, 0.0f, 0.0f, 0.0f};
        if (atoms.w != -1) {
            pos3 = oldPos[atoms.w];
            xpj3 = posDelta[atoms.w];
        }

        // Precompute quantities.

        float4 rij1 = pos-pos1;
        float4 rij2 = pos-pos2;
        float4 rij3 = pos-pos3;
        float rij1sq = rij1.x*rij1.x + rij1.y*rij1.y + rij1.z*rij1.z;
        float rij2sq = rij2.x*rij2.x + rij2.y*rij2.y + rij2.z*rij2.z;
        float rij3sq = rij3.x*rij3.x + rij3.y*rij3.y + rij3.z*rij3.z;
        float ld1 = d2-rij1sq;
        float ld2 = d2-rij2sq;
        float ld3 = d2-rij3sq;

        // Iterate until convergence.

        bool converged = false;
        int iteration = 0;
        while (iteration < 15 && !converged) {
            converged = true;
            float4 rpij = xpi-xpj1;
            float rpsqij = rpij.x*rpij.x + rpij.y*rpij.y + rpij.z*rpij.z;
            float rrpr = rij1.x*rpij.x + rij1.y*rpij.y + rij1.z*rpij.z;
            float diff = fabs(ld1-2.0f*rrpr-rpsqij) / (d2*tol);
            if (diff >= 1.0f) {
                float acor  = (ld1-2.0f*rrpr-rpsqij)*avgMass / (rrpr+rij1sq);
                float4 dr = rij1*acor;
                xpi.xyz += dr.xyz*invMassCentral;
                xpj1.xyz -= dr.xyz*invMassPeripheral;
                converged = false;
            }
            if (atoms.z != -1) {
                rpij.xyz = xpi.xyz-xpj2.xyz;
                rpsqij = rpij.x*rpij.x + rpij.y*rpij.y + rpij.z*rpij.z;
                rrpr = rij2.x*rpij.x + rij2.y*rpij.y + rij2.z*rpij.z;
                diff = fabs(ld2-2.0f*rrpr-rpsqij) / (d2*tol);
                if (diff >= 1.0f) {
                    float acor  = (ld2 - 2.0f*rrpr - rpsqij)*avgMass / (rrpr + rij2sq);
                    float4 dr = rij2*acor;
                    xpi.xyz += dr.xyz*invMassCentral;
                    xpj2.xyz -= dr.xyz*invMassPeripheral;
                    converged = false;
                }
            }
            if (atoms.w != -1) {
                rpij.xyz = xpi.xyz-xpj3.xyz;
                rpsqij = rpij.x*rpij.x + rpij.y*rpij.y + rpij.z*rpij.z;
                rrpr = rij3.x*rpij.x + rij3.y*rpij.y + rij3.z*rpij.z;
                diff = fabs(ld3 - 2.0f*rrpr - rpsqij) / (d2*tol);
                if (diff >= 1.0f) {
                    float acor  = (ld3-2.0f*rrpr-rpsqij)*avgMass / (rrpr+rij3sq);
                    float4 dr = rij3*acor;
                    xpi.xyz += dr.xyz*invMassCentral;
                    xpj3.xyz -= dr.xyz*invMassPeripheral;
                    converged = false;
                }
            }
            iteration++;
        }

        // Record the new positions.

        newDelta[atoms.x] = xpi;
        newDelta[atoms.y] = xpj1;
        if (atoms.z != -1)
            newDelta[atoms.z] = xpj2;
        if (atoms.w != -1)
            newDelta[atoms.w] = xpj3;
        index += get_global_size(0);
    }
}
