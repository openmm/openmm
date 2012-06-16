/**
 * Enforce constraints on SHAKE clusters
 */

extern "C" __global__ void applyShakeToHydrogens(int numClusters, real tol, const real4* __restrict__ oldPos, real4* __restrict__ posDelta, const int4* __restrict__ clusterAtoms, const float4* __restrict__ clusterParams) {
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    while (index < numClusters) {
        // Load the data for this cluster.

        int4 atoms = clusterAtoms[index];
        float4 params = clusterParams[index];
        real4 pos = oldPos[atoms.x];
        real4 xpi = posDelta[atoms.x];
        real4 pos1 = oldPos[atoms.y];
        real4 xpj1 = posDelta[atoms.y];
        real4 pos2 = make_real4(0);
        real4 xpj2 = make_real4(0);
        real invMassCentral = params.x;
        real avgMass = params.y;
        float d2 = params.z;
        float invMassPeripheral = params.w;
        if (atoms.z != -1) {
            pos2 = oldPos[atoms.z];
            xpj2 = posDelta[atoms.z];
        }
        real4 pos3 = make_real4(0);
        real4 xpj3 = make_real4(0);
        if (atoms.w != -1) {
            pos3 = oldPos[atoms.w];
            xpj3 = posDelta[atoms.w];
        }

        // Precompute quantities.

        real3 rij1 = make_real3(pos.x-pos1.x, pos.y-pos1.y, pos.z-pos1.z);
        real3 rij2 = make_real3(pos.x-pos2.x, pos.y-pos2.y, pos.z-pos2.z);
        real3 rij3 = make_real3(pos.x-pos3.x, pos.y-pos3.y, pos.z-pos3.z);
        real rij1sq = rij1.x*rij1.x + rij1.y*rij1.y + rij1.z*rij1.z;
        real rij2sq = rij2.x*rij2.x + rij2.y*rij2.y + rij2.z*rij2.z;
        real rij3sq = rij3.x*rij3.x + rij3.y*rij3.y + rij3.z*rij3.z;
        real ld1 = d2-rij1sq;
        real ld2 = d2-rij2sq;
        real ld3 = d2-rij3sq;

        // Iterate until convergence.

        bool converged = false;
        int iteration = 0;
        while (iteration < 15 && !converged) {
            converged = true;
#ifdef CONSTRAIN_VELOCITIES
            real3 rpij = make_real3(xpi.x-xpj1.x, xpi.y-xpj1.y, xpi.z-xpj1.z);
            real rrpr = rpij.x*rij1.x + rpij.y*rij1.y + rpij.z*rij1.z;
            real delta = -2.0f*avgMass*rrpr/rij1sq;
            real3 dr = rij1*delta;
            xpi.x += dr.x*invMassCentral;
            xpi.y += dr.y*invMassCentral;
            xpi.z += dr.z*invMassCentral;
            xpj1.x -= dr.x*invMassPeripheral;
            xpj1.y -= dr.y*invMassPeripheral;
            xpj1.z -= dr.z*invMassPeripheral;
            if (fabs(delta) > tol)
                converged = false;
            if (atoms.z != -1) {
                rpij = make_real3(xpi.x-xpj2.x, xpi.y-xpj2.y, xpi.z-xpj2.z);
                rrpr = rpij.x*rij2.x + rpij.y*rij2.y + rpij.z*rij2.z;
                delta = -2.0f*avgMass*rrpr/rij2sq;
                dr = rij2*delta;
                xpi.x += dr.x*invMassCentral;
                xpi.y += dr.y*invMassCentral;
                xpi.z += dr.z*invMassCentral;
                xpj2.x -= dr.x*invMassPeripheral;
                xpj2.y -= dr.y*invMassPeripheral;
                xpj2.z -= dr.z*invMassPeripheral;
                if (fabs(delta) > tol)
                    converged = false;
            }
            if (atoms.w != -1) {
                rpij = make_real3(xpi.x-xpj3.x, xpi.y-xpj3.y, xpi.z-xpj3.z);
                rrpr = rpij.x*rij3.x + rpij.y*rij3.y + rpij.z*rij3.z;
                delta = -2.0f*avgMass*rrpr/rij3sq;
                dr = rij3*delta;
                xpi.x += dr.x*invMassCentral;
                xpi.y += dr.y*invMassCentral;
                xpi.z += dr.z*invMassCentral;
                xpj3.x -= dr.x*invMassPeripheral;
                xpj3.y -= dr.y*invMassPeripheral;
                xpj3.z -= dr.z*invMassPeripheral;
                if (fabs(delta) > tol)
                    converged = false;
            }
#else
            real3 rpij = make_real3(xpi.x-xpj1.x, xpi.y-xpj1.y, xpi.z-xpj1.z);
            real rpsqij = rpij.x*rpij.x + rpij.y*rpij.y + rpij.z*rpij.z;
            real rrpr = rij1.x*rpij.x + rij1.y*rpij.y + rij1.z*rpij.z;
            real diff = fabs(ld1-2.0f*rrpr-rpsqij) / (d2*tol);
            if (diff >= 1.0f) {
                real acor  = (ld1-2.0f*rrpr-rpsqij)*avgMass / (rrpr+rij1sq);
                real3 dr = rij1*acor;
                xpi.x += dr.x*invMassCentral;
                xpi.y += dr.y*invMassCentral;
                xpi.z += dr.z*invMassCentral;
                xpj1.x -= dr.x*invMassPeripheral;
                xpj1.y -= dr.y*invMassPeripheral;
                xpj1.z -= dr.z*invMassPeripheral;
                converged = false;
            }
            if (atoms.z != -1) {
                rpij = make_real3(xpi.x-xpj2.x, xpi.y-xpj2.y, xpi.z-xpj2.z);
                rpsqij = rpij.x*rpij.x + rpij.y*rpij.y + rpij.z*rpij.z;
                rrpr = rij2.x*rpij.x + rij2.y*rpij.y + rij2.z*rpij.z;
                diff = fabs(ld2-2.0f*rrpr-rpsqij) / (d2*tol);
                if (diff >= 1.0f) {
                    real acor  = (ld2 - 2.0f*rrpr - rpsqij)*avgMass / (rrpr + rij2sq);
                    real3 dr = rij2*acor;
                    xpi.x += dr.x*invMassCentral;
                    xpi.y += dr.y*invMassCentral;
                    xpi.z += dr.z*invMassCentral;
                    xpj2.x -= dr.x*invMassPeripheral;
                    xpj2.y -= dr.y*invMassPeripheral;
                    xpj2.z -= dr.z*invMassPeripheral;
                    converged = false;
                }
            }
            if (atoms.w != -1) {
                rpij = make_real3(xpi.x-xpj3.x, xpi.y-xpj3.y, xpi.z-xpj3.z);
                rpsqij = rpij.x*rpij.x + rpij.y*rpij.y + rpij.z*rpij.z;
                rrpr = rij3.x*rpij.x + rij3.y*rpij.y + rij3.z*rpij.z;
                diff = fabs(ld3 - 2.0f*rrpr - rpsqij) / (d2*tol);
                if (diff >= 1.0f) {
                    real acor  = (ld3-2.0f*rrpr-rpsqij)*avgMass / (rrpr+rij3sq);
                    real3 dr = rij3*acor;
                    xpi.x += dr.x*invMassCentral;
                    xpi.y += dr.y*invMassCentral;
                    xpi.z += dr.z*invMassCentral;
                    xpj3.x -= dr.x*invMassPeripheral;
                    xpj3.y -= dr.y*invMassPeripheral;
                    xpj3.z -= dr.z*invMassPeripheral;
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
        index += blockDim.x*gridDim.x;
    }
}
