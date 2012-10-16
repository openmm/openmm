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
 * Enforce constraints on SHAKE clusters
 */

__kernel void applyShakeToHydrogens(int numClusters, mixed tol, __global const real4* restrict oldPos, __global const real4* restrict posCorrection, __global mixed4* restrict posDelta, __global const int4* restrict clusterAtoms, __global const float4* restrict clusterParams) {
    int index = get_global_id(0);
    while (index < numClusters) {
        // Load the data for this cluster.

        int4 atoms = clusterAtoms[index];
        float4 params = clusterParams[index];
        mixed4 pos = loadPos(oldPos, posCorrection, atoms.x);
        mixed4 xpi = posDelta[atoms.x];
        mixed4 pos1 = loadPos(oldPos, posCorrection, atoms.y);
        mixed4 xpj1 = posDelta[atoms.y];
        mixed4 pos2 = {0.0f, 0.0f, 0.0f, 0.0f};
        mixed4 xpj2 = {0.0f, 0.0f, 0.0f, 0.0f};
        float invMassCentral = params.x;
        float avgMass = params.y;
        float d2 = params.z;
        float invMassPeripheral = params.w;
        if (atoms.z != -1) {
            pos2 = loadPos(oldPos, posCorrection, atoms.z);
            xpj2 = posDelta[atoms.z];
        }
        mixed4 pos3 = {0.0f, 0.0f, 0.0f, 0.0f};
        mixed4 xpj3 = {0.0f, 0.0f, 0.0f, 0.0f};
        if (atoms.w != -1) {
            pos3 = loadPos(oldPos, posCorrection, atoms.w);
            xpj3 = posDelta[atoms.w];
        }

        // Precompute quantities.

        mixed4 rij1 = pos-pos1;
        mixed4 rij2 = pos-pos2;
        mixed4 rij3 = pos-pos3;
        mixed rij1sq = rij1.x*rij1.x + rij1.y*rij1.y + rij1.z*rij1.z;
        mixed rij2sq = rij2.x*rij2.x + rij2.y*rij2.y + rij2.z*rij2.z;
        mixed rij3sq = rij3.x*rij3.x + rij3.y*rij3.y + rij3.z*rij3.z;
        mixed ld1 = d2-rij1sq;
        mixed ld2 = d2-rij2sq;
        mixed ld3 = d2-rij3sq;

        // Iterate until convergence.

        int converged = false;
        int iteration = 0;
        while (iteration < 15 && !converged) {
            converged = true;
#ifdef CONSTRAIN_VELOCITIES
            mixed4 rpij = xpi-xpj1;
            mixed rrpr = rpij.x*rij1.x + rpij.y*rij1.y + rpij.z*rij1.z;
            mixed delta = -2.0f*avgMass*rrpr/rij1sq;
            mixed4 dr = rij1*delta;
            xpi.xyz += dr.xyz*invMassCentral;
            xpj1.xyz -= dr.xyz*invMassPeripheral;
            if (fabs(delta) > tol)
                converged = false;
            if (atoms.z != -1) {
                rpij = xpi-xpj2;
                rrpr = rpij.x*rij2.x + rpij.y*rij2.y + rpij.z*rij2.z;
                delta = -2.0f*avgMass*rrpr/rij2sq;
                dr = rij2*delta;
                xpi.xyz += dr.xyz*invMassCentral;
                xpj2.xyz -= dr.xyz*invMassPeripheral;
                if (fabs(delta) > tol)
                    converged = false;
            }
            if (atoms.w != -1) {
                rpij = xpi-xpj3;
                rrpr = rpij.x*rij3.x + rpij.y*rij3.y + rpij.z*rij3.z;
                delta = -2.0f*avgMass*rrpr/rij3sq;
                dr = rij3*delta;
                xpi.xyz += dr.xyz*invMassCentral;
                xpj3.xyz -= dr.xyz*invMassPeripheral;
                if (fabs(delta) > tol)
                    converged = false;
            }
#else
            mixed4 rpij = xpi-xpj1;
            mixed rpsqij = rpij.x*rpij.x + rpij.y*rpij.y + rpij.z*rpij.z;
            mixed rrpr = rij1.x*rpij.x + rij1.y*rpij.y + rij1.z*rpij.z;
            mixed diff = fabs(ld1-2.0f*rrpr-rpsqij) / (d2*tol);
            if (diff >= 1.0f) {
                mixed acor  = (ld1-2.0f*rrpr-rpsqij)*avgMass / (rrpr+rij1sq);
                mixed4 dr = rij1*acor;
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
                    mixed acor  = (ld2 - 2.0f*rrpr - rpsqij)*avgMass / (rrpr + rij2sq);
                    mixed4 dr = rij2*acor;
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
                    mixed acor  = (ld3-2.0f*rrpr-rpsqij)*avgMass / (rrpr+rij3sq);
                    mixed4 dr = rij3*acor;
                    xpi.xyz += dr.xyz*invMassCentral;
                    xpj3.xyz -= dr.xyz*invMassPeripheral;
                    converged = false;
                }
            }
#endif
            iteration++;
        }

        // Record the new positions.

        posDelta[atoms.x] = xpi;
        posDelta[atoms.y] = xpj1;
        if (atoms.z != -1)
            posDelta[atoms.z] = xpj2;
        if (atoms.w != -1)
            posDelta[atoms.w] = xpj3;
        index += get_global_size(0);
    }
}
