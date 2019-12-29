/**
 * Perform the first step of Velocity Verlet integration.
 */

__kernel void integrateVelocityVerletPart1(int numAtoms, int numPairs, int paddedNumAtoms, __global const mixed2* restrict dt, __global const real4* restrict posq,
        __global const real4* restrict posqCorrection, __global mixed4* restrict velm, __global const real4* restrict force, __global mixed4* restrict posDelta,
        __global const int* restrict atomList, __global const int2* restrict pairList) {
    const mixed2 stepSize = dt[0];
    const mixed dtPos = stepSize.y;
    const mixed dtVel = 0.5f*(stepSize.x+stepSize.y);
    int index = get_global_id(0);
    while (index < numAtoms) {
        int atom = atomList[index];
        mixed4 velocity = velm[atom];
        if (velocity.w != 0.0) {
#ifdef USE_MIXED_PRECISION
            real4 pos1 = posq[atom];
            real4 pos2 = posqCorrection[atom];
            mixed4 pos = (mixed4) (pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
#else
            real4 pos = posq[atom];
#endif
            velocity.x += 0.5f * dtVel*force[atom].x*velocity.w;
            velocity.y += 0.5f * dtVel*force[atom].y*velocity.w;
            velocity.z += 0.5f * dtVel*force[atom].z*velocity.w;
            pos.x = velocity.x*dtPos;
            pos.y = velocity.y*dtPos;
            pos.z = velocity.z*dtPos;
            posDelta[atom] = pos;
            velm[atom] = velocity;
        }
        index += get_global_size(0);
    }
    index = get_global_id(0);
    while (index < numPairs){
        int atom1 = pairList[index].x;
        int atom2 = pairList[index].y;
        mixed4 v1 = velm[atom1];
        mixed4 v2 = velm[atom2];
        mixed m1 = v1.w == 0.0f ? 0.0f : 1.0f / v1.w;
        mixed m2 = v2.w == 0.0f ? 0.0f : 1.0f / v2.w;
        mixed mass1fract = m1 / (m1 + m2);
        mixed mass2fract = m2 / (m1 + m2);
        mixed invRedMass = (m1 * m2 != 0.0f) ? (m1 + m2)/(m1 * m2) : 0.0f;
        mixed invTotMass = (m1 + m2 != 0.0f) ? 1.0f /(m1 + m2) : 0.0f;
        mixed3 comVel;
        comVel.x= v1.x*mass1fract + v2.x*mass2fract;
        comVel.y= v1.y*mass1fract + v2.y*mass2fract;
        comVel.z= v1.z*mass1fract + v2.z*mass2fract;
        mixed3 relVel;
        relVel.x= v2.x - v1.x;
        relVel.y= v2.y - v1.y;
        relVel.z= v2.z - v1.z;

        mixed3 comFrc;
        comFrc.x = force[atom1].x + force[atom2].x;
        comFrc.y = force[atom1].y + force[atom2].y;
        comFrc.z = force[atom1].z + force[atom2].z;
        mixed3 relFrc;
        relFrc.x = mass1fract*force[atom2].x - mass2fract*force[atom1].x;
        relFrc.y = mass1fract*force[atom2].y - mass2fract*force[atom1].y;
        relFrc.z = mass1fract*force[atom2].z - mass2fract*force[atom1].z;
        comVel.x += comFrc.x * 0.5f * dtVel * invTotMass;
        comVel.y += comFrc.y * 0.5f * dtVel * invTotMass;
        comVel.z += comFrc.z * 0.5f * dtVel * invTotMass;
        relVel.x += relFrc.x * 0.5f * dtVel * invRedMass;
        relVel.y += relFrc.y * 0.5f * dtVel * invRedMass;
        relVel.z += relFrc.z * 0.5f * dtVel * invRedMass;
#ifdef USE_MIXED_PRECISION
        real4 posv1 = posq[atom1];
        real4 posv2 = posq[atom2];
        real4 posc1 = posqCorrection[atom1];
        real4 posc2 = posqCorrection[atom2];
        mixed4 pos1 = (mixed4) (posv1.x+(mixed)posc1.x, posv1.y+(mixed)posc1.y, posv1.z+(mixed)posc1.z, posv1.w);
        mixed4 pos2 = (mixed4) (posv2.x+(mixed)posc2.x, posv2.y+(mixed)posc2.y, posv2.z+(mixed)posc2.z, posv2.w);
#else
        real4 pos1 = posq[atom1];
        real4 pos2 = posq[atom2];
#endif
        if (v1.w != 0.0f) {
            v1.x = comVel.x - relVel.x*mass2fract;
            v1.y = comVel.y - relVel.y*mass2fract;
            v1.z = comVel.z - relVel.z*mass2fract;
            pos1.x = v1.x*dtPos;
            pos1.y = v1.y*dtPos;
            pos1.z = v1.z*dtPos;
            posDelta[atom1] = pos1;
            velm[atom1] = v1;
        }
        if (v2.w != 0.0f) {
            v2.x = comVel.x + relVel.x*mass1fract;
            v2.y = comVel.y + relVel.y*mass1fract;
            v2.z = comVel.z + relVel.z*mass1fract;
            pos2.x = v2.x*dtPos;
            pos2.y = v2.y*dtPos;
            pos2.z = v2.z*dtPos;
            posDelta[atom2] = pos2;
            velm[atom2] = v2;
        }
        index += get_global_size(0);
     }
}

/**
 * Perform the second step of Velocity Verlet integration.
 */

__kernel void integrateVelocityVerletPart2(int numAtoms, __global mixed2* restrict dt, __global real4* restrict posq,
        __global real4* restrict posqCorrection, __global mixed4* restrict velm, __global const mixed4* restrict posDelta) {
    mixed2 stepSize = dt[0];
    int index = get_global_id(0);
    if (index == 0)
        dt[0].x = stepSize.y;
    while(index < numAtoms) {
        mixed4 velocity = velm[index];
        if (velocity.w != 0.0) {
#ifdef USE_MIXED_PRECISION
            real4 pos1 = posq[index];
            real4 pos2 = posqCorrection[index];
            mixed4 pos = (mixed4) (pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
#else
            real4 pos = posq[index];
#endif
            mixed4 delta = posDelta[index];
            pos.xyz += delta.xyz;
#ifdef USE_MIXED_PRECISION
            posq[index] = (real4) ((real) pos.x, (real) pos.y, (real) pos.z, (real) pos.w);
            posqCorrection[index] = (real4) (pos.x-(real) pos.x, pos.y-(real) pos.y, pos.z-(real) pos.z, 0);
#else
            posq[index] = pos;
#endif
        }
        index += get_global_size(0);
    }
}

/**
 * Perform the third step of Velocity Verlet integration.
 */

__kernel void integrateVelocityVerletPart3(int numAtoms, int numPairs, int paddedNumAtoms, __global mixed2* restrict dt, __global real4* restrict posq,
        __global real4* restrict posqCorrection, __global mixed4* restrict velm, __global const real4* restrict force, __global const mixed4* restrict posDelta,
        __global const int* restrict atomList, __global const int2* __restrict__ pairList) {
    mixed2 stepSize = dt[0];
#ifndef SUPPORTS_DOUBLE_PRECISION
    double oneOverDt = 1.0/stepSize.y;
#else
    float oneOverDt = 1.0f/stepSize.y;
    float correction = (1.0f-oneOverDt*stepSize.y)/stepSize.y;
#endif
    const mixed dtVel = 0.5f*(stepSize.x+stepSize.y);
    int index = get_global_id(0);
    if (index == 0)
        dt[0].x = stepSize.y;
    while(index < numAtoms) {
        int atom = atomList[index];
        mixed4 velocity = velm[atom];
        if (velocity.w != 0.0) {
            mixed4 deltaXconstrained = posDelta[atom];
            velocity.x += 0.5f * dtVel*force[atom].x*velocity.w + (deltaXconstrained.x - velocity.x*stepSize.y)*oneOverDt;
            velocity.y += 0.5f * dtVel*force[atom].y*velocity.w + (deltaXconstrained.y - velocity.y*stepSize.y)*oneOverDt;
            velocity.z += 0.5f * dtVel*force[atom].z*velocity.w + (deltaXconstrained.z - velocity.z*stepSize.y)*oneOverDt;
#ifdef SUPPORTS_DOUBLE_PRECISION
            velocity.x += (deltaXconstrained.x - velocity.x*stepSize.y)*correction;
            velocity.y += (deltaXconstrained.y - velocity.y*stepSize.y)*correction;
            velocity.z += (deltaXconstrained.z - velocity.z*stepSize.y)*correction;
#endif
            velm[atom] = velocity;
        }
        index += get_global_size(0);
    }
    index = get_global_id(0);
    while(index < numPairs) {
        int atom1 = pairList[index].x;
        int atom2 = pairList[index].y;
        mixed4 v1 = velm[atom1];
        mixed4 v2 = velm[atom2];
        mixed m1 = v1.w == 0.0f ? 0.0f : 1.0f / v1.w;
        mixed m2 = v2.w == 0.0f ? 0.0f : 1.0f / v2.w;
        mixed mass1fract = m1 / (m1 + m2);
        mixed mass2fract = m2 / (m1 + m2);
        mixed invRedMass = (m1 * m2 != 0.0f) ? (m1 + m2)/(m1 * m2) : 0.0f;
        mixed invTotMass = (m1 + m2 != 0.0f) ? 1.0f /(m1 + m2) : 0.0f;
        mixed3 comVel;
        comVel.x= v1.x*mass1fract + v2.x*mass2fract;
        comVel.y= v1.y*mass1fract + v2.y*mass2fract;
        comVel.z= v1.z*mass1fract + v2.z*mass2fract;
        mixed3 relVel;
        relVel.x= v2.x - v1.x;
        relVel.y= v2.y - v1.y;
        relVel.z= v2.z - v1.z;

        mixed3 comFrc;
        comFrc.x = force[atom1].x + force[atom2].x;
        comFrc.y = force[atom1].y + force[atom2].y;
        comFrc.z = force[atom1].z + force[atom2].z;
        mixed3 relFrc;
        relFrc.x = mass1fract*force[atom2].x - mass2fract*force[atom1].x;
        relFrc.y = mass1fract*force[atom2].y - mass2fract*force[atom1].y;
        relFrc.z = mass1fract*force[atom2].z - mass2fract*force[atom1].z;
        comVel.x += comFrc.x * 0.5f * dtVel * invTotMass;
        comVel.y += comFrc.y * 0.5f * dtVel * invTotMass;
        comVel.z += comFrc.z * 0.5f * dtVel * invTotMass;
        relVel.x += relFrc.x * 0.5f * dtVel * invRedMass;
        relVel.y += relFrc.y * 0.5f * dtVel * invRedMass;
        relVel.z += relFrc.z * 0.5f * dtVel * invRedMass;
        if (v1.w != 0.0f) {
            mixed4 deltaXconstrained = posDelta[atom1];
            v1.x = comVel.x - relVel.x*mass2fract + (deltaXconstrained.x - v1.x*stepSize.y)*oneOverDt;
            v1.y = comVel.y - relVel.y*mass2fract + (deltaXconstrained.y - v1.y*stepSize.y)*oneOverDt;
            v1.z = comVel.z - relVel.z*mass2fract + (deltaXconstrained.z - v1.z*stepSize.y)*oneOverDt;
#ifdef SUPPORTS_DOUBLE_PRECISION
            v1.x += (deltaXconstrained.x - v1.x*stepSize.y)*correction;
            v1.y += (deltaXconstrained.y - v1.y*stepSize.y)*correction;
            v1.z += (deltaXconstrained.z - v1.z*stepSize.y)*correction;
#endif
            velm[atom1] = v1;
        }
        if (v2.w != 0.0f) {
            mixed4 deltaXconstrained = posDelta[atom2];
            v2.x = comVel.x + relVel.x*mass1fract + (deltaXconstrained.x - v2.x*stepSize.y)*oneOverDt;
            v2.y = comVel.y + relVel.y*mass1fract + (deltaXconstrained.y - v2.y*stepSize.y)*oneOverDt;
            v2.z = comVel.z + relVel.z*mass1fract + (deltaXconstrained.z - v2.z*stepSize.y)*oneOverDt;
#ifdef SUPPORTS_DOUBLE_PRECISION
            v2.x += (deltaXconstrained.x - v2.x*stepSize.y)*correction;
            v2.y += (deltaXconstrained.y - v2.y*stepSize.y)*correction;
            v2.z += (deltaXconstrained.z - v2.z*stepSize.y)*correction;
#endif
            velm[atom2] = v2;
        }
        index += get_global_size(0);
    }
}

__kernel void integrateVelocityVerletHardWall(int numPairs, __global const float* restrict maxPairDistance, 
        __global mixed2* restrict dt, __global real4* restrict posq,
        __global real4* restrict posqCorrection, __global mixed4* restrict velm, 
        __global const int2* restrict pairList, __global const float* __restrict__ pairTemperature) {

    mixed dtPos = dt[0].y;
    mixed maxDelta = (mixed) maxPairDistance[0];
    if (maxDelta > 0){
        int index = get_global_id(0);
        while(index < numPairs) {

            const mixed hardWallScale = sqrt( ((mixed) pairTemperature[index]) * ((mixed) BOLTZ));
            int2 atom = (int2) (pairList[index].x, pairList[index].y);
#ifdef USE_MIXED_PRECISION
            real4 posv1 = posq[atom.x];
            real4 posc1 = posqCorrection[atom.x];
            mixed4 pos1 = (mixed4) (posv1.x+(mixed)posc1.x, posv1.y+(mixed)posc1.y, posv1.z+(mixed)posc1.z, posv1.w);
            real4 posv2 = posq[atom.y];
            real4 posc2 = posqCorrection[atom.y];
            mixed4 pos2 = (mixed4) (posv2.x+(mixed)posc2.x, posv2.y+(mixed)posc2.y, posv2.z+(mixed)posc2.z, posv2.w);
#else
            real4 pos1 = posq[atom.x];
            real4 pos2 = posq[atom.y];
#endif
            mixed3 delta = (mixed3) (
                (mixed) (pos1.x - pos2.x),
                (mixed) (pos1.y - pos2.y),
                (mixed) (pos1.z - pos2.z)
            );
            mixed r = sqrt(delta.x*delta.x + delta.y*delta.y + delta.z*delta.z);
            mixed rInv = 1/r;
            if (rInv*maxDelta < 1.0) {
                // The constraint has been violated, so make the inter-particle distance "bounce"
                // off the hard wall.
                mixed3 bondDir = (mixed3) (delta.x * rInv, delta.y * rInv, delta.z * rInv);
                mixed3 vel1 = (mixed3) (velm[atom.x].x, velm[atom.x].y, velm[atom.x].z);
                mixed3 vel2 = (mixed3) (velm[atom.y].x, velm[atom.y].y, velm[atom.y].z);
                mixed m1 = velm[atom.x].w != 0.0 ? 1.0/velm[atom.x].w : 0.0;
                mixed m2 = velm[atom.y].w != 0.0 ? 1.0/velm[atom.y].w : 0.0;
                mixed invTotMass = (m1 + m2 != 0.0) ? 1.0 /(m1 + m2) : 0.0;
                mixed deltaR = r-maxDelta;
                mixed deltaT = dtPos;
                mixed dt = dtPos;

                mixed dotvr1 = vel1.x*bondDir.x + vel1.y*bondDir.y + vel1.z*bondDir.z;
                mixed3 vb1 = (mixed3) (bondDir.x*dotvr1, bondDir.y*dotvr1, bondDir.z*dotvr1);
                mixed3 vp1 = (mixed3) (vel1.x-vb1.x, vel1.y-vb1.y, vel1.z-vb1.z);
                if (m2 == 0) {
                    // The parent particle is massless, so move only the Drude particle.

                    if (dotvr1 != 0.0)
                        deltaT = deltaR/fabs(dotvr1);
                    if (deltaT > dtPos)
                        deltaT = dtPos;
                    dotvr1 = -dotvr1*hardWallScale/(fabs(dotvr1)*sqrt(m1));
                    mixed dr = -deltaR + deltaT*dotvr1;
                    pos1.x += bondDir.x*dr;
                    pos1.y += bondDir.y*dr;
                    pos1.z += bondDir.z*dr;
                    velm[atom.x] = (mixed4) (vp1.x + bondDir.x*dotvr1, vp1.y + bondDir.y*dotvr1, vp1.z + bondDir.z*dotvr1, velm[atom.x].w);
#ifdef USE_MIXED_PRECISION
                    posq[atom.x] = (real4) ((real) pos1.x, (real) pos1.y, (real) pos1.z, (real) pos1.w);
                    posqCorrection[atom.x] = (real4) (pos1.x-(real) pos1.x, pos1.y-(real) pos1.y, pos1.z-(real) pos1.z, 0);
#else
                    posq[atom.x] = pos1;
#endif
                }
                else {
                    // Move both particles.
                    mixed dotvr2 = vel2.x*bondDir.x + vel2.y*bondDir.y + vel2.z*bondDir.z;
                    mixed3 vb2 = (mixed3) (bondDir.x*dotvr2, bondDir.y*dotvr2, bondDir.z*dotvr2);
                    mixed3 vp2 = (mixed3) (vel2.x-vb2.x, vel2.y-vb2.y, vel2.z-vb2.z);
                    mixed vbCMass = (m1*dotvr1 + m2*dotvr2)*invTotMass;
                    dotvr1 -= vbCMass;
                    dotvr2 -= vbCMass;
                    if (dotvr1 != dotvr2)
                        deltaT = deltaR/fabs(dotvr1-dotvr2);
                    if (deltaT > dt)
                        deltaT = dt;
                    mixed vBond = hardWallScale/sqrt(m1);
                    dotvr1 = -dotvr1*vBond*m2*invTotMass/fabs(dotvr1);
                    dotvr2 = -dotvr2*vBond*m1*invTotMass/fabs(dotvr2);
                    mixed dr1 = -deltaR*m2*invTotMass + deltaT*dotvr1;
                    mixed dr2 = deltaR*m1*invTotMass + deltaT*dotvr2;
                    dotvr1 += vbCMass;
                    dotvr2 += vbCMass;
                    pos1.x += bondDir.x*dr1;
                    pos1.y += bondDir.y*dr1;
                    pos1.z += bondDir.z*dr1;
                    pos2.x += bondDir.x*dr2;
                    pos2.y += bondDir.y*dr2;
                    pos2.z += bondDir.z*dr2;
                    velm[atom.x] = (mixed4) (vp1.x + bondDir.x*dotvr1, vp1.y + bondDir.y*dotvr1, vp1.z + bondDir.z*dotvr1, velm[atom.x].w);
                    velm[atom.y] = (mixed4) (vp2.x + bondDir.x*dotvr2, vp2.y + bondDir.y*dotvr2, vp2.z + bondDir.z*dotvr2, velm[atom.y].w);
#ifdef USE_MIXED_PRECISION
                    posq[atom.x] = (real4) ((real) pos1.x, (real) pos1.y, (real) pos1.z, (real) pos1.w);
                    posq[atom.y] = (real4) ((real) pos2.x, (real) pos2.y, (real) pos2.z, (real) pos2.w);
                    posqCorrection[atom.x] = (real4) (pos1.x-(real) pos1.x, pos1.y-(real) pos1.y, pos1.z-(real) pos1.z, 0);
                    posqCorrection[atom.y] = (real4) (pos2.x-(real) pos2.x, pos2.y-(real) pos2.y, pos2.z-(real) pos2.z, 0);
#else
                    posq[atom.x] = pos1;
                    posq[atom.y] = pos2;
#endif
                }
            }
            index += get_global_size(0);
        }
    }
}

