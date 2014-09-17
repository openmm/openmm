#define WARPS_PER_GROUP (THREAD_BLOCK_SIZE/TILE_SIZE)

typedef struct {
    real3 pos, force, torque, dipole, inducedDipole, inducedDipolePolar;
    real q;
    float thole, damp;
#ifdef INCLUDE_QUADRUPOLES
    real quadrupoleXX, quadrupoleXY, quadrupoleXZ, quadrupoleYY, quadrupoleYZ;
    float padding;
#endif
} AtomData;

__device__ void computeOneInteractionF1(AtomData& atom1, volatile AtomData& atom2, real4 delta, real4 bn, real bn5, float forceFactor, float dScale, float pScale, float mScale, real3& force, real& energy);
__device__ void computeOneInteractionF2(AtomData& atom1, volatile AtomData& atom2, real4 delta, real4 bn, float forceFactor, float dScale, float pScale, float mScale, real3& force, real& energy);
__device__ void computeOneInteractionT1(AtomData& atom1, volatile AtomData& atom2, const real4 delta, const real4 bn, float dScale, float pScale, float mScale);
__device__ void computeOneInteractionT2(AtomData& atom1, volatile AtomData& atom2, const real4 delta, const real4 bn, float dScale, float pScale, float mScale);
__device__ void computeOneInteractionF1NoScale(AtomData& atom1, volatile AtomData& atom2, real4 delta, real4 bn, real bn5, float forceFactor, real3& force, real& energy);
__device__ void computeOneInteractionF2NoScale(AtomData& atom1, volatile AtomData& atom2, real4 delta, real4 bn, float forceFactor, real3& force, real& energy);
__device__ void computeOneInteractionT1NoScale(AtomData& atom1, volatile AtomData& atom2, const real4 delta, const real4 bn);
__device__ void computeOneInteractionT2NoScale(AtomData& atom1, volatile AtomData& atom2, const real4 delta, const real4 bn);

inline __device__ void loadAtomData(AtomData& data, int atom, const real4* __restrict__ posq, const real* __restrict__ labFrameDipole,
        const real* __restrict__ labFrameQuadrupole, const real* __restrict__ inducedDipole, const real* __restrict__ inducedDipolePolar, const float2* __restrict__ dampingAndThole) {
    real4 atomPosq = posq[atom];
    data.pos = make_real3(atomPosq.x, atomPosq.y, atomPosq.z);
    data.q = atomPosq.w;
    data.dipole.x = labFrameDipole[atom*3];
    data.dipole.y = labFrameDipole[atom*3+1];
    data.dipole.z = labFrameDipole[atom*3+2];
#ifdef INCLUDE_QUADRUPOLES
    data.quadrupoleXX = labFrameQuadrupole[atom*5];
    data.quadrupoleXY = labFrameQuadrupole[atom*5+1];
    data.quadrupoleXZ = labFrameQuadrupole[atom*5+2];
    data.quadrupoleYY = labFrameQuadrupole[atom*5+3];
    data.quadrupoleYZ = labFrameQuadrupole[atom*5+4];
#endif
    data.inducedDipole.x = inducedDipole[atom*3];
    data.inducedDipole.y = inducedDipole[atom*3+1];
    data.inducedDipole.z = inducedDipole[atom*3+2];
    data.inducedDipolePolar.x = inducedDipolePolar[atom*3];
    data.inducedDipolePolar.y = inducedDipolePolar[atom*3+1];
    data.inducedDipolePolar.z = inducedDipolePolar[atom*3+2];
    float2 temp = dampingAndThole[atom];
    data.damp = temp.x;
    data.thole = temp.y;
}

__device__ real computeDScaleFactor(unsigned int polarizationGroup, int index) {
    return (polarizationGroup & 1<<index ? 0 : 1);
}

__device__ float computeMScaleFactor(uint2 covalent, int index) {
    int mask = 1<<index;
    bool x = (covalent.x & mask);
    bool y = (covalent.y & mask);
    return (x ? (y ? 0.0f : 0.4f) : (y ? 0.8f : 1.0f));
}

__device__ float computePScaleFactor(uint2 covalent, unsigned int polarizationGroup, int index) {
    int mask = 1<<index;
    bool x = (covalent.x & mask);
    bool y = (covalent.y & mask);
    bool p = (polarizationGroup & mask);
    return (x && y ? 0.0f : (x && p ? 0.5f : 1.0f));
}

__device__ void computeOneInteraction(AtomData& atom1, AtomData& atom2, bool hasExclusions, float dScale, float pScale, float mScale, float forceFactor,
                                      real& energy, real4 periodicBoxSize, real4 invPeriodicBoxSize) {
    real4 delta;
    delta.x = atom2.pos.x - atom1.pos.x;
    delta.y = atom2.pos.y - atom1.pos.y;
    delta.z = atom2.pos.z - atom1.pos.z;

    // periodic box

    delta.x -= floor(delta.x*invPeriodicBoxSize.x+0.5f)*periodicBoxSize.x;
    delta.y -= floor(delta.y*invPeriodicBoxSize.y+0.5f)*periodicBoxSize.y;
    delta.z -= floor(delta.z*invPeriodicBoxSize.z+0.5f)*periodicBoxSize.z;

    delta.w = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
    if (delta.w > CUTOFF_SQUARED)
        return;

    real r = SQRT(delta.w);
    real ralpha = EWALD_ALPHA*r;

    real alsq2 = 2*EWALD_ALPHA*EWALD_ALPHA;
    real alsq2n = 0;
    if (EWALD_ALPHA > 0)
        alsq2n = RECIP(SQRT_PI*EWALD_ALPHA);
    real exp2a = EXP(-(ralpha*ralpha));

    real rr1 = RECIP(r);
    delta.w = rr1;
#ifdef USE_DOUBLE_PRECISION
    const real erfcAlphaR = erfc(ralpha);
#else
    // This approximation for erfc is from Abramowitz and Stegun (1964) p. 299.  They cite the following as
    // the original source: C. Hastings, Jr., Approximations for Digital Computers (1955).  It has a maximum
    // error of 1.5e-7.

    const real t = RECIP(1.0f+0.3275911f*ralpha);
    const real erfcAlphaR = (0.254829592f+(-0.284496736f+(1.421413741f+(-1.453152027f+1.061405429f*t)*t)*t)*t)*t*exp2a;
#endif
    real bn0 = erfcAlphaR*rr1;
    energy += forceFactor*atom1.q*atom2.q*bn0;
    real rr2 = rr1*rr1;
    alsq2n *= alsq2;

    real4 bn;
    bn.x = (bn0+alsq2n*exp2a)*rr2;

    alsq2n *= alsq2;
    bn.y = (3*bn.x+alsq2n*exp2a)*rr2;

    alsq2n *= alsq2;
    bn.z = (5*bn.y+alsq2n*exp2a)*rr2;

    alsq2n *= alsq2;
    bn.w = (7*bn.z+alsq2n*exp2a)*rr2;

    alsq2n *= alsq2;
    real bn5 = (9*bn.w+alsq2n*exp2a)*rr2;

    real3 force;

    if (hasExclusions) {
        computeOneInteractionF1(atom1, atom2, delta, bn, bn5, forceFactor, dScale, pScale, mScale, force, energy);
        computeOneInteractionF2(atom1, atom2, delta, bn, forceFactor, dScale, pScale, mScale, force, energy);
    }
    else {
        computeOneInteractionF1NoScale(atom1, atom2, delta, bn, bn5, forceFactor, force, energy);
        computeOneInteractionF2NoScale(atom1, atom2, delta, bn, forceFactor, force, energy);
    }

    atom1.force += force;
    if (forceFactor == 1)
        atom2.force -= force;
    
    if (hasExclusions) {
        computeOneInteractionT1(atom1, atom2, delta, bn, dScale, pScale, mScale);
        computeOneInteractionT2(atom1, atom2, delta, bn, dScale, pScale, mScale);
    }
    else {
        computeOneInteractionT1NoScale(atom1, atom2, delta, bn);
        computeOneInteractionT2NoScale(atom1, atom2, delta, bn);
    }


    if (forceFactor == 1) {
        // T3 == T1 w/ particles I and J reversed
        // T4 == T2 w/ particles I and J reversed

        delta.x = -delta.x;
        delta.y = -delta.y;
        delta.z = -delta.z;
        if (hasExclusions) {
            computeOneInteractionT1(atom2, atom1, delta, bn, dScale, pScale, mScale);
            computeOneInteractionT2(atom2, atom1, delta, bn, dScale, pScale, mScale);
        }
        else {
            computeOneInteractionT1NoScale(atom2, atom1, delta, bn);
            computeOneInteractionT2NoScale(atom2, atom1, delta, bn);
        }
    }
}

/**
 * Compute the self energy and self torque.
 */
__device__ void computeSelfEnergyAndTorque(AtomData& atom1, real& energy) {
    real term = 2*EWALD_ALPHA*EWALD_ALPHA;
    real fterm = -EWALD_ALPHA/SQRT_PI;
    real cii = atom1.q*atom1.q;
    real dii = dot(atom1.dipole, atom1.dipole);
#ifdef INCLUDE_QUADRUPOLES
    real qii = 2*(atom1.quadrupoleXX*atom1.quadrupoleXX +
                  atom1.quadrupoleYY*atom1.quadrupoleYY +
                  atom1.quadrupoleXX*atom1.quadrupoleYY +
                  atom1.quadrupoleXY*atom1.quadrupoleXY +
                  atom1.quadrupoleXZ*atom1.quadrupoleXZ +
                  atom1.quadrupoleYZ*atom1.quadrupoleYZ);
#else
    real qii = 0;
#endif
    real uii = dot(atom1.dipole, atom1.inducedDipole);
    real selfEnergy = (cii + term*(dii/3 + 2*term*qii/5));
    selfEnergy += term*uii/3;
    selfEnergy *= fterm;
    energy += selfEnergy;

    // self-torque for PME

    real3 ui = atom1.inducedDipole+atom1.inducedDipolePolar;
    atom1.torque += ((2/(real) 3)*(EWALD_ALPHA*EWALD_ALPHA*EWALD_ALPHA)/SQRT_PI)*cross(atom1.dipole, ui);
}

/**
 * Compute electrostatic interactions.
 */
extern "C" __global__ void computeElectrostatics(
        unsigned long long* __restrict__ forceBuffers, unsigned long long* __restrict__ torqueBuffers, real* __restrict__ energyBuffer,
        const real4* __restrict__ posq, const uint2* __restrict__ covalentFlags, const unsigned int* __restrict__ polarizationGroupFlags,
        const ushort2* __restrict__ exclusionTiles, unsigned int startTileIndex, unsigned int numTileIndices,
#ifdef USE_CUTOFF
        const int* __restrict__ tiles, const unsigned int* __restrict__ interactionCount, real4 periodicBoxSize, real4 invPeriodicBoxSize, unsigned int maxTiles, const real4* __restrict__ blockCenter, const unsigned int* __restrict__ interactingAtoms,
#endif
        const real* __restrict__ labFrameDipole, const real* __restrict__ labFrameQuadrupole, const real* __restrict__ inducedDipole,
        const real* __restrict__ inducedDipolePolar, const float2* __restrict__ dampingAndThole) {
    const unsigned int totalWarps = (blockDim.x*gridDim.x)/TILE_SIZE;
    const unsigned int warp = (blockIdx.x*blockDim.x+threadIdx.x)/TILE_SIZE;
    const unsigned int tgx = threadIdx.x & (TILE_SIZE-1);
    const unsigned int tbx = threadIdx.x - tgx;
    real energy = 0;
    __shared__ AtomData localData[THREAD_BLOCK_SIZE];

    // First loop: process tiles that contain exclusions.
    
    const unsigned int firstExclusionTile = FIRST_EXCLUSION_TILE+warp*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    const unsigned int lastExclusionTile = FIRST_EXCLUSION_TILE+(warp+1)*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    for (int pos = firstExclusionTile; pos < lastExclusionTile; pos++) {
        const ushort2 tileIndices = exclusionTiles[pos];
        const unsigned int x = tileIndices.x;
        const unsigned int y = tileIndices.y;
        AtomData data;
        unsigned int atom1 = x*TILE_SIZE + tgx;
        loadAtomData(data, atom1, posq, labFrameDipole, labFrameQuadrupole, inducedDipole, inducedDipolePolar, dampingAndThole);
        data.force = make_real3(0);
        data.torque = make_real3(0);
        uint2 covalent = covalentFlags[pos*TILE_SIZE+tgx];
        unsigned int polarizationGroup = polarizationGroupFlags[pos*TILE_SIZE+tgx];
        if (x == y) {
            // This tile is on the diagonal.

            localData[threadIdx.x].pos = data.pos;
            localData[threadIdx.x].q = data.q;
            localData[threadIdx.x].dipole = data.dipole;
#ifdef INCLUDE_QUADRUPOLES
            localData[threadIdx.x].quadrupoleXX = data.quadrupoleXX;
            localData[threadIdx.x].quadrupoleXY = data.quadrupoleXY;
            localData[threadIdx.x].quadrupoleXZ = data.quadrupoleXZ;
            localData[threadIdx.x].quadrupoleYY = data.quadrupoleYY;
            localData[threadIdx.x].quadrupoleYZ = data.quadrupoleYZ;
#endif
            localData[threadIdx.x].inducedDipole = data.inducedDipole;
            localData[threadIdx.x].inducedDipolePolar = data.inducedDipolePolar;
            localData[threadIdx.x].thole = data.thole;
            localData[threadIdx.x].damp = data.damp;

            // Compute forces.

            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                int atom2 = y*TILE_SIZE+j;
                if (atom1 != atom2 && atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                    float d = computeDScaleFactor(polarizationGroup, j);
                    float p = computePScaleFactor(covalent, polarizationGroup, j);
                    float m = computeMScaleFactor(covalent, j);
                    computeOneInteraction(data, localData[tbx+j], true, d, p, m, 0.5f, energy, periodicBoxSize, invPeriodicBoxSize);
                }
            }
            if (atom1 < NUM_ATOMS)
                computeSelfEnergyAndTorque(data, energy);
            data.force *= -ENERGY_SCALE_FACTOR;
            data.torque *= ENERGY_SCALE_FACTOR;
            atomicAdd(&forceBuffers[atom1], static_cast<unsigned long long>((long long) (data.force.x*0x100000000)));
            atomicAdd(&forceBuffers[atom1+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.y*0x100000000)));
            atomicAdd(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.z*0x100000000)));
            atomicAdd(&torqueBuffers[atom1], static_cast<unsigned long long>((long long) (data.torque.x*0x100000000)));
            atomicAdd(&torqueBuffers[atom1+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.torque.y*0x100000000)));
            atomicAdd(&torqueBuffers[atom1+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.torque.z*0x100000000)));
        }
        else {
            // This is an off-diagonal tile.

            unsigned int j = y*TILE_SIZE + tgx;
            loadAtomData(localData[threadIdx.x], j, posq, labFrameDipole, labFrameQuadrupole, inducedDipole, inducedDipolePolar, dampingAndThole);
            localData[threadIdx.x].force = make_real3(0);
            localData[threadIdx.x].torque = make_real3(0);
            unsigned int tj = tgx;
            for (j = 0; j < TILE_SIZE; j++) {
                int atom2 = y*TILE_SIZE+tj;
                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                    float d = computeDScaleFactor(polarizationGroup, tj);
                    float p = computePScaleFactor(covalent, polarizationGroup, tj);
                    float m = computeMScaleFactor(covalent, tj);
                    computeOneInteraction(data, localData[tbx+tj], true, d, p, m, 1, energy, periodicBoxSize, invPeriodicBoxSize);
                }
                tj = (tj + 1) & (TILE_SIZE - 1);
            }
            data.force *= -ENERGY_SCALE_FACTOR;
            data.torque *= ENERGY_SCALE_FACTOR;
            localData[threadIdx.x].force *= -ENERGY_SCALE_FACTOR;
            localData[threadIdx.x].torque *= ENERGY_SCALE_FACTOR;
            unsigned int offset = x*TILE_SIZE + tgx;
            atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long) (data.force.x*0x100000000)));
            atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.y*0x100000000)));
            atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.z*0x100000000)));
            atomicAdd(&torqueBuffers[offset], static_cast<unsigned long long>((long long) (data.torque.x*0x100000000)));
            atomicAdd(&torqueBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.torque.y*0x100000000)));
            atomicAdd(&torqueBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.torque.z*0x100000000)));
            offset = y*TILE_SIZE + tgx;
            atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.x*0x100000000)));
            atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.y*0x100000000)));
            atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.z*0x100000000)));
            atomicAdd(&torqueBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].torque.x*0x100000000)));
            atomicAdd(&torqueBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].torque.y*0x100000000)));
            atomicAdd(&torqueBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].torque.z*0x100000000)));
        }
    }

    // Second loop: tiles without exclusions, either from the neighbor list (with cutoff) or just enumerating all
    // of them (no cutoff).

#ifdef USE_CUTOFF
    const unsigned int numTiles = interactionCount[0];
    int pos = (int) (numTiles > maxTiles ? startTileIndex+warp*(long long)numTileIndices/totalWarps : warp*(long long)numTiles/totalWarps);
    int end = (int) (numTiles > maxTiles ? startTileIndex+(warp+1)*(long long)numTileIndices/totalWarps : (warp+1)*(long long)numTiles/totalWarps);
#else
    const unsigned int numTiles = numTileIndices;
    int pos = (int) (startTileIndex+warp*(long long)numTiles/totalWarps);
    int end = (int) (startTileIndex+(warp+1)*(long long)numTiles/totalWarps);
#endif
    int skipBase = 0;
    int currentSkipIndex = tbx;
    __shared__ int atomIndices[THREAD_BLOCK_SIZE];
    __shared__ volatile int skipTiles[THREAD_BLOCK_SIZE];
    skipTiles[threadIdx.x] = -1;
    
    while (pos < end) {
        bool includeTile = true;

        // Extract the coordinates of this tile.
        
        int x, y;
#ifdef USE_CUTOFF
        if (numTiles <= maxTiles)
            x = tiles[pos];
        else
#endif
        {
            y = (int) floor(NUM_BLOCKS+0.5f-SQRT((NUM_BLOCKS+0.5f)*(NUM_BLOCKS+0.5f)-2*pos));
            x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
            if (x < y || x >= NUM_BLOCKS) { // Occasionally happens due to roundoff error.
                y += (x < y ? -1 : 1);
                x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
            }

            // Skip over tiles that have exclusions, since they were already processed.

            while (skipTiles[tbx+TILE_SIZE-1] < pos) {
                if (skipBase+tgx < NUM_TILES_WITH_EXCLUSIONS) {
                    ushort2 tile = exclusionTiles[skipBase+tgx];
                    skipTiles[threadIdx.x] = tile.x + tile.y*NUM_BLOCKS - tile.y*(tile.y+1)/2;
                }
                else
                    skipTiles[threadIdx.x] = end;
                skipBase += TILE_SIZE;            
                currentSkipIndex = tbx;
            }
            while (skipTiles[currentSkipIndex] < pos)
                currentSkipIndex++;
            includeTile = (skipTiles[currentSkipIndex] != pos);
        }
        if (includeTile) {
            unsigned int atom1 = x*TILE_SIZE + tgx;

            // Load atom data for this tile.

            AtomData data;
            loadAtomData(data, atom1, posq, labFrameDipole, labFrameQuadrupole, inducedDipole, inducedDipolePolar, dampingAndThole);
            data.force = make_real3(0);
            data.torque = make_real3(0);
#ifdef USE_CUTOFF
            unsigned int j = (numTiles <= maxTiles ? interactingAtoms[pos*TILE_SIZE+tgx] : y*TILE_SIZE + tgx);
#else
            unsigned int j = y*TILE_SIZE + tgx;
#endif
            atomIndices[threadIdx.x] = j;
            loadAtomData(localData[threadIdx.x], j, posq, labFrameDipole, labFrameQuadrupole, inducedDipole, inducedDipolePolar, dampingAndThole);
            localData[threadIdx.x].force = make_real3(0);
            localData[threadIdx.x].torque = make_real3(0);

            // Compute forces.

            unsigned int tj = tgx;
            for (j = 0; j < TILE_SIZE; j++) {
                int atom2 = atomIndices[tbx+tj];
                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                    computeOneInteraction(data, localData[tbx+tj], false, 1, 1, 1, 1, energy, periodicBoxSize, invPeriodicBoxSize);
                }
                tj = (tj + 1) & (TILE_SIZE - 1);
            }
            data.force *= -ENERGY_SCALE_FACTOR;
            data.torque *= ENERGY_SCALE_FACTOR;
            localData[threadIdx.x].force *= -ENERGY_SCALE_FACTOR;
            localData[threadIdx.x].torque *= ENERGY_SCALE_FACTOR;

            // Write results.

            unsigned int offset = x*TILE_SIZE + tgx;
            atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long) (data.force.x*0x100000000)));
            atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.y*0x100000000)));
            atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.z*0x100000000)));
            atomicAdd(&torqueBuffers[offset], static_cast<unsigned long long>((long long) (data.torque.x*0x100000000)));
            atomicAdd(&torqueBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.torque.y*0x100000000)));
            atomicAdd(&torqueBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.torque.z*0x100000000)));
#ifdef USE_CUTOFF
            offset = atomIndices[threadIdx.x];
#else
            offset = y*TILE_SIZE + tgx;
#endif
            atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.x*0x100000000)));
            atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.y*0x100000000)));
            atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.z*0x100000000)));
            atomicAdd(&torqueBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].torque.x*0x100000000)));
            atomicAdd(&torqueBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].torque.y*0x100000000)));
            atomicAdd(&torqueBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].torque.z*0x100000000)));
        }
        pos++;
    }
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x] += energy*ENERGY_SCALE_FACTOR;
}
