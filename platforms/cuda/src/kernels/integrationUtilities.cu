/**
 * Generate random numbers
 */
extern "C" __global__ void generateRandomNumbers(int numValues, float4* __restrict__ random, uint4* __restrict__ seed) {
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    uint4 state = seed[index];
    unsigned int carry = 0;
    while (index < numValues) {
        float4 value;

        // Generate first two values.

        state.x = state.x * 69069 + 1;
        state.y ^= state.y << 13;
        state.y ^= state.y >> 17;
        state.y ^= state.y << 5;
        unsigned int k = (state.z >> 2) + (state.w >> 3) + (carry >> 2);
        unsigned int m = state.w + state.w + state.z + carry;
        state.z = state.w;
        state.w = m;
        carry = k >> 30;
        float x1 = (float)max(state.x + state.y + state.w, 0x00000001u) / (float)0xffffffff;
        state.x = state.x * 69069 + 1;
        state.y ^= state.y << 13;
        state.y ^= state.y >> 17;
        state.y ^= state.y << 5;
        x1 = SQRT(-2.0f * LOG(x1));
        k = (state.z >> 2) + (state.w >> 3) + (carry >> 2);
        m = state.w + state.w + state.z + carry;
        state.z = state.w;
        state.w = m;
        carry = k >> 30;
        float x2 = (float)(state.x + state.y + state.w) / (float)0xffffffff;
        value.x = x1 * COS(2.0f * 3.14159265f * x2);
        value.y = x1 * SIN(2.0f * 3.14159265f * x2);

        // Generate next two values.

        state.x = state.x * 69069 + 1;
        state.y ^= state.y << 13;
        state.y ^= state.y >> 17;
        state.y ^= state.y << 5;
        k = (state.z >> 2) + (state.w >> 3) + (carry >> 2);
        m = state.w + state.w + state.z + carry;
        state.z = state.w;
        state.w = m;
        carry = k >> 30;
        float x3 = (float)max(state.x + state.y + state.w, 0x00000001u) / (float)0xffffffff;
        state.x = state.x * 69069 + 1;
        state.y ^= state.y << 13;
        state.y ^= state.y >> 17;
        state.y ^= state.y << 5;
        x3 = SQRT(-2.0f * LOG(x3));
        k = (state.z >> 2) + (state.w >> 3) + (carry >> 2);
        m = state.w + state.w + state.z + carry;
        state.z = state.w;
        state.w = m;
        carry = k >> 30;
        float x4 = (float)(state.x + state.y + state.w) / (float)0xffffffff;
        value.z = x3 * COS(2.0f * 3.14159265f * x4);
        value.w = x3 * SIN(2.0f * 3.14159265f * x4);

        // Record the values.

        random[index] = value;
        index += blockDim.x*gridDim.x;
    }
    seed[blockIdx.x*blockDim.x+threadIdx.x] = state;
}

/**
 * Load the position of a particle.
 */
inline __device__ mixed4 loadPos(const real4* __restrict__ posq, const real4* __restrict__ posqCorrection, int index) {
#ifdef USE_MIXED_PRECISION
    real4 pos1 = posq[index];
    real4 pos2 = posqCorrection[index];
    return make_mixed4(pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
#else
    return posq[index];
#endif
}

/**
 * Store the position of a particle.
 */
inline __device__ void storePos(real4* __restrict__ posq, real4* __restrict__ posqCorrection, int index, mixed4 pos) {
#ifdef USE_MIXED_PRECISION
    posq[index] = make_real4((real) pos.x, (real) pos.y, (real) pos.z, (real) pos.w);
    posqCorrection[index] = make_real4(pos.x-(real) pos.x, pos.y-(real) pos.y, pos.z-(real) pos.z, 0);
#else
    posq[index] = pos;
#endif
}

/**
 * Enforce constraints on SHAKE clusters
 */
extern "C" __global__ void applyShakeToPositions(int numClusters, mixed tol, const real4* __restrict__ oldPos, real4* __restrict__ posCorrection, mixed4* __restrict__ posDelta, const int4* __restrict__ clusterAtoms, const float4* __restrict__ clusterParams) {
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    while (index < numClusters) {
        // Load the data for this cluster.

        int4 atoms = clusterAtoms[index];
        float4 params = clusterParams[index];
        mixed4 pos = loadPos(oldPos, posCorrection, atoms.x);
        mixed4 xpi = posDelta[atoms.x];
        mixed4 pos1 = loadPos(oldPos, posCorrection, atoms.y);
        mixed4 xpj1 = posDelta[atoms.y];
        mixed4 pos2 = make_mixed4(0);
        mixed4 xpj2 = make_mixed4(0);
        float invMassCentral = params.x;
        float avgMass = params.y;
        float d2 = params.z;
        float invMassPeripheral = params.w;
        if (atoms.z != -1) {
            pos2 = loadPos(oldPos, posCorrection, atoms.z);
            xpj2 = posDelta[atoms.z];
        }
        mixed4 pos3 = make_mixed4(0);
        mixed4 xpj3 = make_mixed4(0);
        if (atoms.w != -1) {
            pos3 = loadPos(oldPos, posCorrection, atoms.w);
            xpj3 = posDelta[atoms.w];
        }

        // Precompute quantities.

        mixed3 rij1 = make_mixed3(pos.x-pos1.x, pos.y-pos1.y, pos.z-pos1.z);
        mixed3 rij2 = make_mixed3(pos.x-pos2.x, pos.y-pos2.y, pos.z-pos2.z);
        mixed3 rij3 = make_mixed3(pos.x-pos3.x, pos.y-pos3.y, pos.z-pos3.z);
        mixed rij1sq = rij1.x*rij1.x + rij1.y*rij1.y + rij1.z*rij1.z;
        mixed rij2sq = rij2.x*rij2.x + rij2.y*rij2.y + rij2.z*rij2.z;
        mixed rij3sq = rij3.x*rij3.x + rij3.y*rij3.y + rij3.z*rij3.z;
        mixed ld1 = d2-rij1sq;
        mixed ld2 = d2-rij2sq;
        mixed ld3 = d2-rij3sq;

        // Iterate until convergence.

        bool converged = false;
        int iteration = 0;
        while (iteration < 15 && !converged) {
            converged = true;
            mixed3 rpij = make_mixed3(xpi.x-xpj1.x, xpi.y-xpj1.y, xpi.z-xpj1.z);
            mixed rpsqij = rpij.x*rpij.x + rpij.y*rpij.y + rpij.z*rpij.z;
            mixed rrpr = rij1.x*rpij.x + rij1.y*rpij.y + rij1.z*rpij.z;
            mixed diff = fabs(ld1-2.0f*rrpr-rpsqij) / (d2*tol);
            if (diff >= 1.0f) {
                mixed acor  = (ld1-2.0f*rrpr-rpsqij)*avgMass / (rrpr+rij1sq);
                mixed3 dr = rij1*acor;
                xpi.x += dr.x*invMassCentral;
                xpi.y += dr.y*invMassCentral;
                xpi.z += dr.z*invMassCentral;
                xpj1.x -= dr.x*invMassPeripheral;
                xpj1.y -= dr.y*invMassPeripheral;
                xpj1.z -= dr.z*invMassPeripheral;
                converged = false;
            }
            if (atoms.z != -1) {
                rpij = make_mixed3(xpi.x-xpj2.x, xpi.y-xpj2.y, xpi.z-xpj2.z);
                rpsqij = rpij.x*rpij.x + rpij.y*rpij.y + rpij.z*rpij.z;
                rrpr = rij2.x*rpij.x + rij2.y*rpij.y + rij2.z*rpij.z;
                diff = fabs(ld2-2.0f*rrpr-rpsqij) / (d2*tol);
                if (diff >= 1.0f) {
                    mixed acor  = (ld2 - 2.0f*rrpr - rpsqij)*avgMass / (rrpr + rij2sq);
                    mixed3 dr = rij2*acor;
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
                rpij = make_mixed3(xpi.x-xpj3.x, xpi.y-xpj3.y, xpi.z-xpj3.z);
                rpsqij = rpij.x*rpij.x + rpij.y*rpij.y + rpij.z*rpij.z;
                rrpr = rij3.x*rpij.x + rij3.y*rpij.y + rij3.z*rpij.z;
                diff = fabs(ld3 - 2.0f*rrpr - rpsqij) / (d2*tol);
                if (diff >= 1.0f) {
                    mixed acor  = (ld3-2.0f*rrpr-rpsqij)*avgMass / (rrpr+rij3sq);
                    mixed3 dr = rij3*acor;
                    xpi.x += dr.x*invMassCentral;
                    xpi.y += dr.y*invMassCentral;
                    xpi.z += dr.z*invMassCentral;
                    xpj3.x -= dr.x*invMassPeripheral;
                    xpj3.y -= dr.y*invMassPeripheral;
                    xpj3.z -= dr.z*invMassPeripheral;
                    converged = false;
                }
            }
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

/**
 * Enforce velocity constraints on SHAKE clusters
 */
extern "C" __global__ void applyShakeToVelocities(int numClusters, mixed tol, const real4* __restrict__ oldPos, real4* __restrict__ posCorrection, mixed4* __restrict__ posDelta, const int4* __restrict__ clusterAtoms, const float4* __restrict__ clusterParams) {
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    while (index < numClusters) {
        // Load the data for this cluster.

        int4 atoms = clusterAtoms[index];
        float4 params = clusterParams[index];
        mixed4 pos = loadPos(oldPos, posCorrection, atoms.x);
        mixed4 xpi = posDelta[atoms.x];
        mixed4 pos1 = loadPos(oldPos, posCorrection, atoms.y);
        mixed4 xpj1 = posDelta[atoms.y];
        mixed4 pos2 = make_mixed4(0);
        mixed4 xpj2 = make_mixed4(0);
        float invMassCentral = params.x;
        float avgMass = params.y;
        float invMassPeripheral = params.w;
        if (atoms.z != -1) {
            pos2 = loadPos(oldPos, posCorrection, atoms.z);
            xpj2 = posDelta[atoms.z];
        }
        mixed4 pos3 = make_mixed4(0);
        mixed4 xpj3 = make_mixed4(0);
        if (atoms.w != -1) {
            pos3 = loadPos(oldPos, posCorrection, atoms.w);
            xpj3 = posDelta[atoms.w];
        }

        // Precompute quantities.

        mixed3 rij1 = make_mixed3(pos.x-pos1.x, pos.y-pos1.y, pos.z-pos1.z);
        mixed3 rij2 = make_mixed3(pos.x-pos2.x, pos.y-pos2.y, pos.z-pos2.z);
        mixed3 rij3 = make_mixed3(pos.x-pos3.x, pos.y-pos3.y, pos.z-pos3.z);
        mixed rij1sq = rij1.x*rij1.x + rij1.y*rij1.y + rij1.z*rij1.z;
        mixed rij2sq = rij2.x*rij2.x + rij2.y*rij2.y + rij2.z*rij2.z;
        mixed rij3sq = rij3.x*rij3.x + rij3.y*rij3.y + rij3.z*rij3.z;

        // Iterate until convergence.

        bool converged = false;
        int iteration = 0;
        while (iteration < 15 && !converged) {
            converged = true;
            mixed3 rpij = make_mixed3(xpi.x-xpj1.x, xpi.y-xpj1.y, xpi.z-xpj1.z);
            mixed rrpr = rpij.x*rij1.x + rpij.y*rij1.y + rpij.z*rij1.z;
            mixed delta = -2.0f*avgMass*rrpr/rij1sq;
            mixed3 dr = rij1*delta;
            xpi.x += dr.x*invMassCentral;
            xpi.y += dr.y*invMassCentral;
            xpi.z += dr.z*invMassCentral;
            xpj1.x -= dr.x*invMassPeripheral;
            xpj1.y -= dr.y*invMassPeripheral;
            xpj1.z -= dr.z*invMassPeripheral;
            if (fabs(delta) > tol)
                converged = false;
            if (atoms.z != -1) {
                rpij = make_mixed3(xpi.x-xpj2.x, xpi.y-xpj2.y, xpi.z-xpj2.z);
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
                rpij = make_mixed3(xpi.x-xpj3.x, xpi.y-xpj3.y, xpi.z-xpj3.z);
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

/**
 * Enforce constraints on SETTLE clusters
 */
extern "C" __global__ void applySettleToPositions(int numClusters, mixed tol, const real4* __restrict__ oldPos, real4* __restrict__ posCorrection, mixed4* __restrict__ posDelta, const mixed4* __restrict__ velm, const int4* __restrict__ clusterAtoms, const float2* __restrict__ clusterParams) {
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    while (index < numClusters) {
        // Load the data for this cluster.

        int4 atoms = clusterAtoms[index];
        float2 params = clusterParams[index];
        mixed4 apos0 = loadPos(oldPos, posCorrection, atoms.x);
        mixed4 xp0 = posDelta[atoms.x];
        mixed4 apos1 = loadPos(oldPos, posCorrection, atoms.y);
        mixed4 xp1 = posDelta[atoms.y];
        mixed4 apos2 = loadPos(oldPos, posCorrection, atoms.z);
        mixed4 xp2 = posDelta[atoms.z];
        mixed m0 = 1/velm[atoms.x].w;
        mixed m1 = 1/velm[atoms.y].w;
        mixed m2 = 1/velm[atoms.z].w;

        // Apply the SETTLE algorithm.

        mixed xb0 = apos1.x-apos0.x;
        mixed yb0 = apos1.y-apos0.y;
        mixed zb0 = apos1.z-apos0.z;
        mixed xc0 = apos2.x-apos0.x;
        mixed yc0 = apos2.y-apos0.y;
        mixed zc0 = apos2.z-apos0.z;

        mixed invTotalMass = 1/(m0+m1+m2);
        mixed xcom = (xp0.x*m0 + (xb0+xp1.x)*m1 + (xc0+xp2.x)*m2) * invTotalMass;
        mixed ycom = (xp0.y*m0 + (yb0+xp1.y)*m1 + (yc0+xp2.y)*m2) * invTotalMass;
        mixed zcom = (xp0.z*m0 + (zb0+xp1.z)*m1 + (zc0+xp2.z)*m2) * invTotalMass;

        mixed xa1 = xp0.x - xcom;
        mixed ya1 = xp0.y - ycom;
        mixed za1 = xp0.z - zcom;
        mixed xb1 = xb0 + xp1.x - xcom;
        mixed yb1 = yb0 + xp1.y - ycom;
        mixed zb1 = zb0 + xp1.z - zcom;
        mixed xc1 = xc0 + xp2.x - xcom;
        mixed yc1 = yc0 + xp2.y - ycom;
        mixed zc1 = zc0 + xp2.z - zcom;

        mixed xaksZd = yb0*zc0 - zb0*yc0;
        mixed yaksZd = zb0*xc0 - xb0*zc0;
        mixed zaksZd = xb0*yc0 - yb0*xc0;
        mixed xaksXd = ya1*zaksZd - za1*yaksZd;
        mixed yaksXd = za1*xaksZd - xa1*zaksZd;
        mixed zaksXd = xa1*yaksZd - ya1*xaksZd;
        mixed xaksYd = yaksZd*zaksXd - zaksZd*yaksXd;
        mixed yaksYd = zaksZd*xaksXd - xaksZd*zaksXd;
        mixed zaksYd = xaksZd*yaksXd - yaksZd*xaksXd;

        mixed axlng = sqrt(xaksXd*xaksXd + yaksXd*yaksXd + zaksXd*zaksXd);
        mixed aylng = sqrt(xaksYd*xaksYd + yaksYd*yaksYd + zaksYd*zaksYd);
        mixed azlng = sqrt(xaksZd*xaksZd + yaksZd*yaksZd + zaksZd*zaksZd);
        mixed trns11 = xaksXd / axlng;
        mixed trns21 = yaksXd / axlng;
        mixed trns31 = zaksXd / axlng;
        mixed trns12 = xaksYd / aylng;
        mixed trns22 = yaksYd / aylng;
        mixed trns32 = zaksYd / aylng;
        mixed trns13 = xaksZd / azlng;
        mixed trns23 = yaksZd / azlng;
        mixed trns33 = zaksZd / azlng;

        mixed xb0d = trns11*xb0 + trns21*yb0 + trns31*zb0;
        mixed yb0d = trns12*xb0 + trns22*yb0 + trns32*zb0;
        mixed xc0d = trns11*xc0 + trns21*yc0 + trns31*zc0;
        mixed yc0d = trns12*xc0 + trns22*yc0 + trns32*zc0;
        mixed za1d = trns13*xa1 + trns23*ya1 + trns33*za1;
        mixed xb1d = trns11*xb1 + trns21*yb1 + trns31*zb1;
        mixed yb1d = trns12*xb1 + trns22*yb1 + trns32*zb1;
        mixed zb1d = trns13*xb1 + trns23*yb1 + trns33*zb1;
        mixed xc1d = trns11*xc1 + trns21*yc1 + trns31*zc1;
        mixed yc1d = trns12*xc1 + trns22*yc1 + trns32*zc1;
        mixed zc1d = trns13*xc1 + trns23*yc1 + trns33*zc1;

        //                                        --- Step2  A2' ---

        float rc = 0.5f*params.y;
        mixed rb = sqrt(params.x*params.x-rc*rc);
        mixed ra = rb*(m1+m2)*invTotalMass;
        rb -= ra;
        mixed sinphi = za1d/ra;
        mixed cosphi = sqrt(1-sinphi*sinphi);
        mixed sinpsi = (zb1d-zc1d) / (2*rc*cosphi);
        mixed cospsi = sqrt(1-sinpsi*sinpsi);

        mixed ya2d =   ra*cosphi;
        mixed xb2d = - rc*cospsi;
        mixed yb2d = - rb*cosphi - rc*sinpsi*sinphi;
        mixed yc2d = - rb*cosphi + rc*sinpsi*sinphi;
        mixed xb2d2 = xb2d*xb2d;
        mixed hh2 = 4.0f*xb2d2 + (yb2d-yc2d)*(yb2d-yc2d) + (zb1d-zc1d)*(zb1d-zc1d);
        mixed deltx = 2.0f*xb2d + sqrt(4.0f*xb2d2 - hh2 + params.y*params.y);
        xb2d -= deltx*0.5f;

        //                                        --- Step3  al,be,ga ---

        mixed alpha = (xb2d*(xb0d-xc0d) + yb0d*yb2d + yc0d*yc2d);
        mixed beta = (xb2d*(yc0d-yb0d) + xb0d*yb2d + xc0d*yc2d);
        mixed gamma = xb0d*yb1d - xb1d*yb0d + xc0d*yc1d - xc1d*yc0d;

        mixed al2be2 = alpha*alpha + beta*beta;
        mixed sintheta = (alpha*gamma - beta*sqrt(al2be2 - gamma*gamma)) / al2be2;

        //                                        --- Step4  A3' ---

        mixed costheta = sqrt(1-sintheta*sintheta);
        mixed xa3d = - ya2d*sintheta;
        mixed ya3d =   ya2d*costheta;
        mixed za3d = za1d;
        mixed xb3d =   xb2d*costheta - yb2d*sintheta;
        mixed yb3d =   xb2d*sintheta + yb2d*costheta;
        mixed zb3d = zb1d;
        mixed xc3d = - xb2d*costheta - yc2d*sintheta;
        mixed yc3d = - xb2d*sintheta + yc2d*costheta;
        mixed zc3d = zc1d;

        //                                        --- Step5  A3 ---

        mixed xa3 = trns11*xa3d + trns12*ya3d + trns13*za3d;
        mixed ya3 = trns21*xa3d + trns22*ya3d + trns23*za3d;
        mixed za3 = trns31*xa3d + trns32*ya3d + trns33*za3d;
        mixed xb3 = trns11*xb3d + trns12*yb3d + trns13*zb3d;
        mixed yb3 = trns21*xb3d + trns22*yb3d + trns23*zb3d;
        mixed zb3 = trns31*xb3d + trns32*yb3d + trns33*zb3d;
        mixed xc3 = trns11*xc3d + trns12*yc3d + trns13*zc3d;
        mixed yc3 = trns21*xc3d + trns22*yc3d + trns23*zc3d;
        mixed zc3 = trns31*xc3d + trns32*yc3d + trns33*zc3d;

        xp0.x = xcom + xa3;
        xp0.y = ycom + ya3;
        xp0.z = zcom + za3;
        xp1.x = xcom + xb3 - xb0;
        xp1.y = ycom + yb3 - yb0;
        xp1.z = zcom + zb3 - zb0;
        xp2.x = xcom + xc3 - xc0;
        xp2.y = ycom + yc3 - yc0;
        xp2.z = zcom + zc3 - zc0;

        // Record the new positions.

        posDelta[atoms.x] = xp0;
        posDelta[atoms.y] = xp1;
        posDelta[atoms.z] = xp2;
        index += blockDim.x*gridDim.x;
    }
}

/**
 * Enforce velocity constraints on SETTLE clusters
 */
extern "C" __global__ void applySettleToVelocities(int numClusters, mixed tol, const real4* __restrict__ oldPos, real4* __restrict__ posCorrection, mixed4* __restrict__ posDelta, mixed4* __restrict__ velm, const int4* __restrict__ clusterAtoms, const float2* __restrict__ clusterParams) {
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < numClusters; index += blockDim.x*gridDim.x) {
        // Load the data for this cluster.

        int4 atoms = clusterAtoms[index];
        mixed4 apos0 = loadPos(oldPos, posCorrection, atoms.x);
        mixed4 apos1 = loadPos(oldPos, posCorrection, atoms.y);
        mixed4 apos2 = loadPos(oldPos, posCorrection, atoms.z);
        mixed4 v0 = velm[atoms.x];
        mixed4 v1 = velm[atoms.y];
        mixed4 v2 = velm[atoms.z];
        
        // Compute intermediate quantities: the atom masses, the bond directions, the relative velocities,
        // and the angle cosines and sines.
        
        mixed mA = 1/v0.w;
        mixed mB = 1/v1.w;
        mixed mC = 1/v2.w;
        mixed3 eAB = make_mixed3(apos1.x-apos0.x, apos1.y-apos0.y, apos1.z-apos0.z);
        mixed3 eBC = make_mixed3(apos2.x-apos1.x, apos2.y-apos1.y, apos2.z-apos1.z);
        mixed3 eCA = make_mixed3(apos0.x-apos2.x, apos0.y-apos2.y, apos0.z-apos2.z);
        eAB *= RSQRT(eAB.x*eAB.x + eAB.y*eAB.y + eAB.z*eAB.z);
        eBC *= RSQRT(eBC.x*eBC.x + eBC.y*eBC.y + eBC.z*eBC.z);
        eCA *= RSQRT(eCA.x*eCA.x + eCA.y*eCA.y + eCA.z*eCA.z);
        mixed vAB = (v1.x-v0.x)*eAB.x + (v1.y-v0.y)*eAB.y + (v1.z-v0.z)*eAB.z;
        mixed vBC = (v2.x-v1.x)*eBC.x + (v2.y-v1.y)*eBC.y + (v2.z-v1.z)*eBC.z;
        mixed vCA = (v0.x-v2.x)*eCA.x + (v0.y-v2.y)*eCA.y + (v0.z-v2.z)*eCA.z;
        mixed cA = -(eAB.x*eCA.x + eAB.y*eCA.y + eAB.z*eCA.z);
        mixed cB = -(eAB.x*eBC.x + eAB.y*eBC.y + eAB.z*eBC.z);
        mixed cC = -(eBC.x*eCA.x + eBC.y*eCA.y + eBC.z*eCA.z);
        mixed s2A = 1-cA*cA;
        mixed s2B = 1-cB*cB;
        mixed s2C = 1-cC*cC;
        
        // Solve the equations.  These are different from those in the SETTLE paper (JCC 13(8), pp. 952-962, 1992), because
        // in going from equations B1 to B2, they make the assumption that mB=mC (but don't bother to mention they're
        // making that assumption).  We allow all three atoms to have different masses.
        
        mixed mABCinv = 1/(mA*mB*mC);
        mixed denom = (((s2A*mB+s2B*mA)*mC+(s2A*mB*mB+2*(cA*cB*cC+1)*mA*mB+s2B*mA*mA))*mC+s2C*mA*mB*(mA+mB))*mABCinv;
        mixed tab = ((cB*cC*mA-cA*mB-cA*mC)*vCA + (cA*cC*mB-cB*mC-cB*mA)*vBC + (s2C*mA*mA*mB*mB*mABCinv+(mA+mB+mC))*vAB)/denom;
        mixed tbc = ((cA*cB*mC-cC*mB-cC*mA)*vCA + (s2A*mB*mB*mC*mC*mABCinv+(mA+mB+mC))*vBC + (cA*cC*mB-cB*mA-cB*mC)*vAB)/denom;
        mixed tca = ((s2B*mA*mA*mC*mC*mABCinv+(mA+mB+mC))*vCA + (cA*cB*mC-cC*mB-cC*mA)*vBC + (cB*cC*mA-cA*mB-cA*mC)*vAB)/denom;
        v0.x += (tab*eAB.x - tca*eCA.x)*v0.w;
        v0.y += (tab*eAB.y - tca*eCA.y)*v0.w;
        v0.z += (tab*eAB.z - tca*eCA.z)*v0.w;
        v1.x += (tbc*eBC.x - tab*eAB.x)*v1.w;
        v1.y += (tbc*eBC.y - tab*eAB.y)*v1.w;
        v1.z += (tbc*eBC.z - tab*eAB.z)*v1.w;
        v2.x += (tca*eCA.x - tbc*eBC.x)*v2.w;
        v2.y += (tca*eCA.y - tbc*eBC.y)*v2.w;
        v2.z += (tca*eCA.z - tbc*eBC.z)*v2.w;
        velm[atoms.x] = v0;
        velm[atoms.y] = v1;
        velm[atoms.z] = v2;
    }
}

/**
 * Compute the direction each CCMA constraint is pointing in.  This is called once at the beginning of constraint evaluation.
 */
extern "C" __global__ void computeCCMAConstraintDirections(const int2* __restrict__ constraintAtoms, mixed4* __restrict__ constraintDistance,
        const real4* __restrict__ atomPositions, const real4* __restrict__ posqCorrection, int* __restrict__ converged) {
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_CCMA_CONSTRAINTS; index += blockDim.x*gridDim.x) {
        // Compute the direction for this constraint.

        int2 atoms = constraintAtoms[index];
        mixed4 dir = constraintDistance[index];
        mixed4 oldPos1 = loadPos(atomPositions, posqCorrection, atoms.x);
        mixed4 oldPos2 = loadPos(atomPositions, posqCorrection, atoms.y);
        dir.x = oldPos1.x-oldPos2.x;
        dir.y = oldPos1.y-oldPos2.y;
        dir.z = oldPos1.z-oldPos2.z;
        constraintDistance[index] = dir;
    }
    if (threadIdx.x == 0 && blockIdx.x == 0) {
        converged[0] = 1;
        converged[1] = 0;
    }
}

/**
 * Compute the force applied by each CCMA position constraint.
 */
extern "C" __global__ void computeCCMAPositionConstraintForce(const int2* __restrict__ constraintAtoms, const mixed4* __restrict__ constraintDistance, const mixed4* __restrict__ atomPositions,
        const mixed* __restrict__ reducedMass, mixed* __restrict__ delta1, int* __restrict__ converged, int* __restrict__ hostConvergedFlag, mixed tol, int iteration) {
    __shared__ int groupConverged;
    if (converged[1-iteration%2]) {
        if (blockIdx.x == 0 && threadIdx.x == 0) {
            converged[iteration%2] = 1;
            hostConvergedFlag[0] = 1;
        }
        return; // The constraint iteration has already converged.
    }
    if (threadIdx.x == 0)
        groupConverged = 1;
    __syncthreads();
    mixed lowerTol = 1-2*tol+tol*tol;
    mixed upperTol = 1+2*tol+tol*tol;
    bool threadConverged = true;
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_CCMA_CONSTRAINTS; index += blockDim.x*gridDim.x) {
        // Compute the force due to this constraint.

        int2 atoms = constraintAtoms[index];
        mixed4 dir = constraintDistance[index];
        mixed4 rp_ij = atomPositions[atoms.x]-atomPositions[atoms.y];
        rp_ij.x += dir.x;
        rp_ij.y += dir.y;
        rp_ij.z += dir.z;
        mixed rrpr = rp_ij.x*dir.x + rp_ij.y*dir.y + rp_ij.z*dir.z;
        mixed d_ij2 = dir.x*dir.x + dir.y*dir.y + dir.z*dir.z;
        mixed rp2 = rp_ij.x*rp_ij.x + rp_ij.y*rp_ij.y + rp_ij.z*rp_ij.z;
        mixed dist2 = dir.w*dir.w;
        mixed diff = dist2 - rp2;
        delta1[index] = (rrpr > d_ij2*1e-6f ? reducedMass[index]*diff/rrpr : 0.0f);
        threadConverged &= (rp2 > lowerTol*dist2 && rp2 < upperTol*dist2);
    }
    if (groupConverged && !threadConverged)
        groupConverged = 0;
    __syncthreads();
    if (threadIdx.x == 0 && !groupConverged)
        converged[iteration%2] = 0;
}

/**
 * Compute the force applied by each CCMA velocity constraint.
 */
extern "C" __global__ void computeCCMAVelocityConstraintForce(const int2* __restrict__ constraintAtoms, const mixed4* __restrict__ constraintDistance, const mixed4* __restrict__ atomPositions,
        const mixed* __restrict__ reducedMass, mixed* __restrict__ delta1, int* __restrict__ converged, int* __restrict__ hostConvergedFlag, mixed tol, int iteration) {
    __shared__ int groupConverged;
    if (converged[1-iteration%2]) {
        if (blockIdx.x == 0 && threadIdx.x == 0) {
            converged[iteration%2] = 1;
            hostConvergedFlag[0] = 1;
        }
        return; // The constraint iteration has already converged.
    }
    if (threadIdx.x == 0)
        groupConverged = 1;
    __syncthreads();
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_CCMA_CONSTRAINTS; index += blockDim.x*gridDim.x) {
        // Compute the force due to this constraint.

        int2 atoms = constraintAtoms[index];
        mixed4 dir = constraintDistance[index];
        mixed4 rp_ij = atomPositions[atoms.x]-atomPositions[atoms.y];
        mixed rrpr = rp_ij.x*dir.x + rp_ij.y*dir.y + rp_ij.z*dir.z;
        mixed d_ij2 = dir.x*dir.x + dir.y*dir.y + dir.z*dir.z;
        delta1[index] = -2*reducedMass[index]*rrpr/d_ij2;

        // See whether it has converged.

        if (groupConverged && fabs(delta1[index]) > tol) {
            groupConverged = 0;
            converged[iteration%2] = 0;
        }
    }
}

/**
 * Multiply the vector of CCMA constraint forces by the constraint matrix.
 */
extern "C" __global__ void multiplyByCCMAConstraintMatrix(const mixed* __restrict__ delta1, mixed* __restrict__ delta2, const int* __restrict__ constraintMatrixColumn,
        const mixed* __restrict__ constraintMatrixValue, const int* __restrict__ converged, int iteration) {
    if (converged[iteration%2])
        return; // The constraint iteration has already converged.

    // Multiply by the inverse constraint matrix.

    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_CCMA_CONSTRAINTS; index += blockDim.x*gridDim.x) {
        mixed sum = 0;
        for (int i = 0; ; i++) {
            int element = index+i*NUM_CCMA_CONSTRAINTS;
            int column = constraintMatrixColumn[element];
            if (column >= NUM_CCMA_CONSTRAINTS)
                break;
            sum += delta1[column]*constraintMatrixValue[element];
        }
        delta2[index] = sum;
    }
}

/**
 * Update the atom positions based on CCMA constraint forces.
 */
extern "C" __global__ void updateCCMAAtomPositions(const int* __restrict__ numAtomConstraints, const int* __restrict__ atomConstraints, const mixed4* __restrict__ constraintDistance,
        mixed4* __restrict__ atomPositions, const mixed4* __restrict__ velm, const mixed* __restrict__ delta1, const mixed* __restrict__ delta2, int* __restrict__ converged, int iteration) {
    if (blockIdx.x == 0 && threadIdx.x == 0)
        converged[1-iteration%2] = 1;
    if (converged[iteration%2])
        return; // The constraint iteration has already converged.
    mixed damping = (iteration < 2 ? 0.5f : 1.0f);
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_ATOMS; index += blockDim.x*gridDim.x) {
        // Compute the new position of this atom.

        mixed4 atomPos = atomPositions[index];
        mixed invMass = velm[index].w;
        int num = numAtomConstraints[index];
        for (int i = 0; i < num; i++) {
            int constraint = atomConstraints[index+i*NUM_ATOMS];
            bool forward = (constraint > 0);
            constraint = (forward ? constraint-1 : -constraint-1);
            mixed constraintForce = damping*invMass*delta2[constraint];
            constraintForce = (forward ? constraintForce : -constraintForce);
            mixed4 dir = constraintDistance[constraint];
            atomPos.x += constraintForce*dir.x;
            atomPos.y += constraintForce*dir.y;
            atomPos.z += constraintForce*dir.z;
        }
        atomPositions[index] = atomPos;
    }
}

/**
 * Compute the positions of virtual sites
 */
extern "C" __global__ void computeVirtualSites(real4* __restrict__ posq, real4* __restrict__ posqCorrection, const int4* __restrict__ avg2Atoms, const real2* __restrict__ avg2Weights,
        const int4* __restrict__ avg3Atoms, const real4* __restrict__ avg3Weights,
        const int4* __restrict__ outOfPlaneAtoms, const real4* __restrict__ outOfPlaneWeights,
        const int* __restrict__ localCoordsIndex, const int* __restrict__ localCoordsAtoms,
        const real* __restrict__ localCoordsWeights, const real4* __restrict__ localCoordsPos,
        const int* __restrict__ localCoordsStartIndex) {
    
    // Two particle average sites.
    
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_2_AVERAGE; index += blockDim.x*gridDim.x) {
        int4 atoms = avg2Atoms[index];
        real2 weights = avg2Weights[index];
        mixed4 pos = loadPos(posq, posqCorrection, atoms.x);
        mixed4 pos1 = loadPos(posq, posqCorrection, atoms.y);
        mixed4 pos2 = loadPos(posq, posqCorrection, atoms.z);
        pos.x = pos1.x*weights.x + pos2.x*weights.y;
        pos.y = pos1.y*weights.x + pos2.y*weights.y;
        pos.z = pos1.z*weights.x + pos2.z*weights.y;
        storePos(posq, posqCorrection, atoms.x, pos);
    }
    
    // Three particle average sites.
    
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_3_AVERAGE; index += blockDim.x*gridDim.x) {
        int4 atoms = avg3Atoms[index];
        real4 weights = avg3Weights[index];
        mixed4 pos = loadPos(posq, posqCorrection, atoms.x);
        mixed4 pos1 = loadPos(posq, posqCorrection, atoms.y);
        mixed4 pos2 = loadPos(posq, posqCorrection, atoms.z);
        mixed4 pos3 = loadPos(posq, posqCorrection, atoms.w);
        pos.x = pos1.x*weights.x + pos2.x*weights.y + pos3.x*weights.z;
        pos.y = pos1.y*weights.x + pos2.y*weights.y + pos3.y*weights.z;
        pos.z = pos1.z*weights.x + pos2.z*weights.y + pos3.z*weights.z;
        storePos(posq, posqCorrection, atoms.x, pos);
    }
    
    // Out of plane sites.
    
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_OUT_OF_PLANE; index += blockDim.x*gridDim.x) {
        int4 atoms = outOfPlaneAtoms[index];
        real4 weights = outOfPlaneWeights[index];
        mixed4 pos = loadPos(posq, posqCorrection, atoms.x);
        mixed4 pos1 = loadPos(posq, posqCorrection, atoms.y);
        mixed4 pos2 = loadPos(posq, posqCorrection, atoms.z);
        mixed4 pos3 = loadPos(posq, posqCorrection, atoms.w);
        mixed4 v12 = pos2-pos1;
        mixed4 v13 = pos3-pos1;
        mixed3 cr = cross(v12, v13);
        pos.x = pos1.x + v12.x*weights.x + v13.x*weights.y + cr.x*weights.z;
        pos.y = pos1.y + v12.y*weights.x + v13.y*weights.y + cr.y*weights.z;
        pos.z = pos1.z + v12.z*weights.x + v13.z*weights.y + cr.z*weights.z;
        storePos(posq, posqCorrection, atoms.x, pos);
    }
    
    // Local coordinates sites.
    
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_LOCAL_COORDS; index += blockDim.x*gridDim.x) {
        int siteAtomIndex = localCoordsIndex[index];
        int start = localCoordsStartIndex[index];
        int end = localCoordsStartIndex[index+1];
        mixed3 origin = make_mixed3(0), xdir = make_mixed3(0), ydir = make_mixed3(0);
        for (int j = start; j < end; j++) {
            mixed3 pos = trimTo3(loadPos(posq, posqCorrection, localCoordsAtoms[j]));
            origin += pos*localCoordsWeights[3*j];
            xdir += pos*localCoordsWeights[3*j+1];
            ydir += pos*localCoordsWeights[3*j+2];
        }
        mixed3 zdir = cross(xdir, ydir);
        mixed normXdir = sqrt(xdir.x*xdir.x+xdir.y*xdir.y+xdir.z*xdir.z);
        mixed normZdir = sqrt(zdir.x*zdir.x+zdir.y*zdir.y+zdir.z*zdir.z);
        mixed invNormXdir = (normXdir > 0 ? 1/normXdir : 0);
        mixed invNormZdir = (normZdir > 0 ? 1/normZdir : 0);
        xdir *= invNormXdir;
        zdir *= invNormZdir;
        ydir = cross(zdir, xdir);
        real4 localPosition_4 = localCoordsPos[index];
        mixed3 localPosition = make_mixed3(localPosition_4.x, localPosition_4.y, localPosition_4.z);
        mixed4 pos = loadPos(posq, posqCorrection, siteAtomIndex);
        pos.x = origin.x + xdir.x*localPosition.x + ydir.x*localPosition.y + zdir.x*localPosition.z;
        pos.y = origin.y + xdir.y*localPosition.x + ydir.y*localPosition.y + zdir.y*localPosition.z;
        pos.z = origin.z + xdir.z*localPosition.x + ydir.z*localPosition.y + zdir.z*localPosition.z;
        storePos(posq, posqCorrection, siteAtomIndex, pos);
    }
}

inline __device__ real3 loadForce(int index, long long* __restrict__ force) {
    real scale = 1/((real) 0x100000000);
    return make_real3(scale*force[index], scale*force[index+PADDED_NUM_ATOMS], scale*force[index+PADDED_NUM_ATOMS*2]);
}

inline __device__ void addForce(int index, long long* __restrict__ force, real3 value) {
    unsigned long long* f = (unsigned long long*) force;
    atomicAdd(&f[index], static_cast<unsigned long long>((long long) (value.x*0x100000000)));
    atomicAdd(&f[index+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (value.y*0x100000000)));
    atomicAdd(&f[index+PADDED_NUM_ATOMS*2], static_cast<unsigned long long>((long long) (value.z*0x100000000)));
}

/**
 * Distribute forces from virtual sites to the atoms they are based on.
 */
extern "C" __global__ void distributeVirtualSiteForces(const real4* __restrict__ posq, const real4* __restrict__ posqCorrection, long long* __restrict__ force,
        const int4* __restrict__ avg2Atoms, const real2* __restrict__ avg2Weights,
        const int4* __restrict__ avg3Atoms, const real4* __restrict__ avg3Weights,
        const int4* __restrict__ outOfPlaneAtoms, const real4* __restrict__ outOfPlaneWeights,
        const int* __restrict__ localCoordsIndex, const int* __restrict__ localCoordsAtoms,
        const real* __restrict__ localCoordsWeights, const real4* __restrict__ localCoordsPos,
        const int* __restrict__ localCoordsStartIndex) {
    
    // Two particle average sites.
    
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_2_AVERAGE; index += blockDim.x*gridDim.x) {
        int4 atoms = avg2Atoms[index];
        real2 weights = avg2Weights[index];
        real3 f = loadForce(atoms.x, force);
        addForce(atoms.y, force, f*weights.x);
        addForce(atoms.z, force, f*weights.y);
    }
    
    // Three particle average sites.
    
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_3_AVERAGE; index += blockDim.x*gridDim.x) {
        int4 atoms = avg3Atoms[index];
        real4 weights = avg3Weights[index];
        real3 f = loadForce(atoms.x, force);
        addForce(atoms.y, force, f*weights.x);
        addForce(atoms.z, force, f*weights.y);
        addForce(atoms.w, force, f*weights.z);
    }
    
    // Out of plane sites.
    
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_OUT_OF_PLANE; index += blockDim.x*gridDim.x) {
        int4 atoms = outOfPlaneAtoms[index];
        real4 weights = outOfPlaneWeights[index];
        mixed4 pos1 = loadPos(posq, posqCorrection, atoms.y);
        mixed4 pos2 = loadPos(posq, posqCorrection, atoms.z);
        mixed4 pos3 = loadPos(posq, posqCorrection, atoms.w);
        mixed4 v12 = pos2-pos1;
        mixed4 v13 = pos3-pos1;
        real3 f = loadForce(atoms.x, force);
        real3 fp2 = make_real3((real) (weights.x*f.x - weights.z*v13.z*f.y + weights.z*v13.y*f.z),
                   (real) (weights.z*v13.z*f.x + weights.x*f.y - weights.z*v13.x*f.z),
                   (real) (-weights.z*v13.y*f.x + weights.z*v13.x*f.y + weights.x*f.z));
        real3 fp3 = make_real3((real) (weights.y*f.x + weights.z*v12.z*f.y - weights.z*v12.y*f.z),
                   (real) (-weights.z*v12.z*f.x + weights.y*f.y + weights.z*v12.x*f.z),
                   (real) (weights.z*v12.y*f.x - weights.z*v12.x*f.y + weights.y*f.z));
        addForce(atoms.y, force, f-fp2-fp3);
        addForce(atoms.z, force, fp2);
        addForce(atoms.w, force, fp3);
    }
    
    // Local coordinates sites.
    
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_LOCAL_COORDS; index += blockDim.x*gridDim.x) {
        int siteAtomIndex = localCoordsIndex[index];
        int start = localCoordsStartIndex[index];
        int end = localCoordsStartIndex[index+1];
        mixed3 origin = make_mixed3(0), xdir = make_mixed3(0), ydir = make_mixed3(0);
        for (int j = start; j < end; j++) {
            mixed3 pos = trimTo3(loadPos(posq, posqCorrection, localCoordsAtoms[j]));
            origin += pos*localCoordsWeights[3*j];
            xdir += pos*localCoordsWeights[3*j+1];
            ydir += pos*localCoordsWeights[3*j+2];
        }
        mixed3 zdir = cross(xdir, ydir);
        mixed normXdir = sqrt(xdir.x*xdir.x+xdir.y*xdir.y+xdir.z*xdir.z);
        mixed normZdir = sqrt(zdir.x*zdir.x+zdir.y*zdir.y+zdir.z*zdir.z);
        mixed invNormXdir = (normXdir > 0 ? 1/normXdir : 0);
        mixed invNormZdir = (normZdir > 0 ? 1/normZdir : 0);
        mixed3 dx = xdir*invNormXdir;
        mixed3 dz = zdir*invNormZdir;
        mixed3 dy = cross(dz, dx);
        real4 localPosition_4 = localCoordsPos[index];
        mixed3 localPosition = make_mixed3(localPosition_4.x, localPosition_4.y, localPosition_4.z);

        // The derivatives for this case are very complicated.  They were computed with SymPy then simplified by hand.

        real3 f = loadForce(siteAtomIndex, force);
        mixed3 fp1 = localPosition*f.x;
        mixed3 fp2 = localPosition*f.y;
        mixed3 fp3 = localPosition*f.z;
        for (int j = start; j < end; j++) {
            real originWeight = localCoordsWeights[3*j];
            real wx = localCoordsWeights[3*j+1];
            real wy = localCoordsWeights[3*j+2];
            mixed wxScaled = wx*invNormXdir;
            mixed t1 = (wx*ydir.x-wy*xdir.x)*invNormZdir;
            mixed t2 = (wx*ydir.y-wy*xdir.y)*invNormZdir;
            mixed t3 = (wx*ydir.z-wy*xdir.z)*invNormZdir;
            mixed sx = t3*dz.y-t2*dz.z;
            mixed sy = t1*dz.z-t3*dz.x;
            mixed sz = t2*dz.x-t1*dz.y;
            real3 fresult = make_real3(0);
            fresult.x += fp1.x*wxScaled*(1-dx.x*dx.x) + fp1.z*(dz.x*sx   ) + fp1.y*((-dx.x*dy.x     )*wxScaled + dy.x*sx - dx.y*t2 - dx.z*t3) + f.x*originWeight;
            fresult.y += fp1.x*wxScaled*( -dx.x*dx.y) + fp1.z*(dz.x*sy+t3) + fp1.y*((-dx.y*dy.x-dz.z)*wxScaled + dy.x*sy + dx.y*t1);
            fresult.z += fp1.x*wxScaled*( -dx.x*dx.z) + fp1.z*(dz.x*sz-t2) + fp1.y*((-dx.z*dy.x+dz.y)*wxScaled + dy.x*sz + dx.z*t1);
            fresult.x += fp2.x*wxScaled*( -dx.y*dx.x) + fp2.z*(dz.y*sx-t3) - fp2.y*(( dx.x*dy.y-dz.z)*wxScaled - dy.y*sx - dx.x*t2);
            fresult.y += fp2.x*wxScaled*(1-dx.y*dx.y) + fp2.z*(dz.y*sy   ) - fp2.y*(( dx.y*dy.y     )*wxScaled - dy.y*sy + dx.x*t1 + dx.z*t3) + f.y*originWeight;
            fresult.z += fp2.x*wxScaled*( -dx.y*dx.z) + fp2.z*(dz.y*sz+t1) - fp2.y*(( dx.z*dy.y+dz.x)*wxScaled - dy.y*sz - dx.z*t2);
            fresult.x += fp3.x*wxScaled*( -dx.z*dx.x) + fp3.z*(dz.z*sx+t2) + fp3.y*((-dx.x*dy.z-dz.y)*wxScaled + dy.z*sx + dx.x*t3);
            fresult.y += fp3.x*wxScaled*( -dx.z*dx.y) + fp3.z*(dz.z*sy-t1) + fp3.y*((-dx.y*dy.z+dz.x)*wxScaled + dy.z*sy + dx.y*t3);
            fresult.z += fp3.x*wxScaled*(1-dx.z*dx.z) + fp3.z*(dz.z*sz   ) + fp3.y*((-dx.z*dy.z     )*wxScaled + dy.z*sz - dx.x*t1 - dx.y*t2) + f.z*originWeight;
            addForce(localCoordsAtoms[j], force, fresult);
        }
    }
}

/**
 * Apply a time shift to the velocities before computing kinetic energy.
 */
extern "C" __global__ void timeShiftVelocities(mixed4* __restrict__ velm, const long long* __restrict__ force, real timeShift) {
    const mixed scale = timeShift/(mixed) 0x100000000;
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_ATOMS; index += blockDim.x*gridDim.x) {
        mixed4 velocity = velm[index];
        if (velocity.w != 0.0) {
            velocity.x += scale*force[index]*velocity.w;
            velocity.y += scale*force[index+PADDED_NUM_ATOMS]*velocity.w;
            velocity.z += scale*force[index+PADDED_NUM_ATOMS*2]*velocity.w;
            velm[index] = velocity;
        }
    }
}
