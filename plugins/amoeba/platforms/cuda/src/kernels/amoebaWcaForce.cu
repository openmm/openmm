#define TILE_SIZE 32
#define WARPS_PER_GROUP (THREAD_BLOCK_SIZE/TILE_SIZE)

typedef struct {
    real3 pos, force;
    float radius, epsilon, padding;
} AtomData;

inline __device__ void loadAtomData(AtomData& data, int atom, const real4* __restrict__ posq, const float2* __restrict__ radiusEpsilon) {
    real4 atomPosq = posq[atom];
    data.pos = make_real3(atomPosq.x, atomPosq.y, atomPosq.z);
    float2 temp = radiusEpsilon[atom];
    data.radius = temp.x;
    data.epsilon = temp.y;
}

__device__ void initParticleParameters(float radius, float epsilon, real& rmixo, real& rmixh, real& emixo, real& emixh) {
    real sqrtEps = SQRT(epsilon);
    real denominator = SQRT(EPSO) + sqrtEps;
    emixo = 4*EPSO*epsilon / (denominator*denominator);
    denominator = SQRT(EPSH) + sqrtEps;
    emixh = 4*EPSH*epsilon / (denominator*denominator);
    real radius2 = radius*radius;
    real rmino2 = RMINO*RMINO; 
    rmixo = 2*(rmino2*RMINO + radius2*radius) / (rmino2 + radius2);
    real rminh2 = RMINH*RMINH;
    rmixh = 2*(rminh2*RMINH + radius2*radius) / (rminh2+radius2);
}

__device__ void computeOneInteraction(AtomData& atom1, AtomData& atom2, real rmixo, real rmixh, real emixo, real emixh, real3& force, real& energy) {
    // get deltaR and r between 2 atoms
    
    force = atom2.pos - atom1.pos;
    real r2 = dot(force, force);
    if (r2 <= 0) {
        force = make_real3(0);
        energy = 0;
        return;
    }
    real rI = RSQRT(r2);
    real r = RECIP(rI);

    real sk = atom2.radius*SHCTD;
    real sk2 = sk*sk;
    if (atom1.radius >= (r+sk)) {
        force = make_real3(0);
        energy = 0;
        return;
    }

    real rmax = atom1.radius > (r - sk) ? atom1.radius : (r - sk);
    real lik = rmax;
    real lik2 = lik*lik;
    real lik3 = lik2*lik;
    real lik4 = lik2*lik2;
 
    real uik = (r+sk) < rmixo ? (r+sk) : rmixo;
    real uik2 = uik*uik;
    real uik3 = uik2*uik;
    real uik4 = uik2*uik2;

    real term = 4*M_PI/(48*r)*(3*(lik4-uik4) - 8*r*(lik3-uik3) + 6*(r2-sk2)*(lik2-uik2));

    real r3 = r2*r;
    real dl1 = lik2*(-lik2 + 2*(r2 + sk2));
    real dl2 = lik*(-lik3 + 4*lik2*r - 6*lik*r2 + 2*lik*sk2 + 4*r3 - 4*r*sk2);
    real dl = atom1.radius > (r-sk)? dl1 : dl2;

    real du1 = uik2*(-uik2 + 2*(r2 + sk2));
    real du2 = uik*(-uik3 + 4*uik2*r - 2*uik*(3*r2 - sk2) + 4*r*(r2 - sk2));
    real du = (r+sk) > rmixo ? -du1 : -du2;

    real mask2 = lik < rmixo ? 1 : 0;
    real sum = -mask2*(emixo*term);
    real de = -mask2*emixo*M_PI*(dl+du)/(4*r2);

    uik = (r+sk) < rmixh ? (r+sk) : rmixh;
    uik2 = uik*uik;
    uik3 = uik2*uik;
    uik4 = uik2*uik2;

    term = (M_PI)/ (12*r) * (3*(lik4-uik4) - 8*r*(lik3-uik3) + 6*(r2-sk2)*(lik2-uik2));

    dl1 = lik2*(-lik2 + 2*r2 + 2*sk2);
    dl2 = lik*(-lik3 + 4*lik2*r - 6*lik*r2 + 2*lik*sk2 + 4*r3 - 4*r*sk2);
    dl = atom1.radius > (r-sk) ? dl1 : dl2;

    du1 = -uik2*(-uik2 + 2*r2 + 2*sk2);
    du2 = -uik*(-uik3 + 4*uik2*r - 6*uik*r2 + 2*uik*sk2 + 4*r3 - 4*r*sk2);
    du = (r+sk) > rmixh ? du1 : du2;

    mask2 = lik < rmixh ? 1 : 0;
    sum -= mask2*(2*emixh*term);
    de -= mask2*(2*emixh*M_PI*(dl+du)/(4*r2));

    uik = r + sk;
    uik2 = uik*uik;
    uik3 = uik2*uik;
    uik4 = uik2*uik2;
    real uik5 = uik4*uik;
    real uik6 = uik3*uik3;
    real uik10 = uik5*uik5;
    real uik11 = uik10*uik;
    real uik12 = uik6*uik6;
    real uik13 = uik12*uik;

    lik = rmax > rmixo ? rmax : rmixo;
    lik2 = lik*lik;
    lik3 = lik2*lik;
    lik4 = lik2*lik2;
    real lik5 = lik4*lik;
    real lik6 = lik3*lik3;
    real lik10 = lik5*lik5;
    real lik11 = lik10*lik;
    real lik12 = lik6*lik6;
    real lik13 = lik12*lik;

    term = 4*M_PI/(120*r*lik5*uik5)*(15*uik*lik*r*(uik4-lik4) - 10*uik2*lik2*(uik3-lik3) + 6*(sk2-r2)*(uik5-lik5));
    dl1 = (-5*lik2 + 3*r2 + 3*sk2)/lik5;
    dl2 = (5*lik3 - 33*lik*r2 - 3*lik*sk2 + 15*(lik2*r+r3-r*sk2))/lik6;
    dl = (atom1.radius > (r-sk)) || (rmax < rmixo) ? -dl1 : dl2;

    du = (-5*uik3 + 33*uik*r2 + 3*uik*sk2 - 15*(uik2*r+r3-r*sk2))/uik6;

    real rmixo7 = rmixo*rmixo*rmixo;
    rmixo7 = rmixo7*rmixo7*rmixo;
    real ao = emixo*rmixo7;

    real idisp = -2*ao*term;
    mask2 = uik > rmixo ? 1 : 0;

    de -= mask2*(2*ao*M_PI*(dl + du)/(15*r2));

    term = 4*M_PI/(2640*r*lik12*uik12) * (120*uik*lik*r*(uik11-lik11) - 66*uik2*lik2*(uik10-lik10) + 55*(sk2-r2)*(uik12-lik12));

    dl1 = (6*lik2 - 5*r2 - 5*sk2)/lik12;
    dl2 = (6*lik3 - 125*lik*r2 - 5*lik*sk2 + 60*(lik2*r+r3-r*sk2))/lik13;
    dl = (atom1.radius > (r-sk)) || (rmax < rmixo) ? dl1 : dl2;

    du = (-6*uik3 + 125*uik*r2 + 5*uik*sk2 - 60*(uik2*r+r3-r*sk2))/uik13;

    de += mask2*(ao*rmixo7*M_PI*(dl + du)/(60*r2));
    real irep = ao*rmixo7*term;
    sum += mask2*(irep + idisp);

    lik = rmax > rmixh ? rmax : rmixh;
    lik2 = lik*lik;
    lik3 = lik2*lik;
    lik4 = lik2*lik2;
    lik5 = lik4*lik;
    lik6 = lik3*lik3;
    lik10 = lik5*lik5;
    lik11 = lik10*lik;
    lik12 = lik6*lik6;
    lik13 = lik12*lik;

    term = 4*M_PI / (120*r*lik5*uik5) * (15*uik*lik*r*(uik4-lik4) - 10*uik2*lik2*(uik3-lik3) + 6*(sk2-r2)*(uik5-lik5));

    dl1 = (-5*lik2 + 3*r2 + 3*sk2)/lik5;
    dl2 = (5*lik3 - 33*lik*r2 - 3*lik*sk2+ 15*(lik2*r+r3-r*sk2))/lik6;
    dl = (atom1.radius > (r-sk)) || (rmax < rmixh) ? -dl1 : dl2;

    du = -(5*uik3 - 33*uik*r2 - 3*uik*sk2 + 15*(uik2*r+r3-r*sk2))/uik6;

    real rmixh7 = rmixh*rmixh*rmixh;
    rmixh7 = rmixh7*rmixh7*rmixh;
    real ah = emixh*rmixh7;

    idisp = -4*ah*term;

    mask2 = uik > rmixh ? 1 : 0;
    de -= mask2*(4*ah*M_PI*(dl + du)/(15*r2));

    term = 4*M_PI / (2640*r*lik12*uik12) * (120*uik*lik*r*(uik11-lik11) - 66*uik2*lik2*(uik10-lik10) + 55*(sk2-r2)*(uik12-lik12));

    dl1 = -(-6*lik2 + 5*r2 + 5*sk2)/lik12;
    dl2 = (6*lik3 - 125*lik*r2 - 5*lik*sk2 + 60*(lik2*r+r3-r*sk2))/lik13;
    dl = ((atom1.radius > (r-sk)) || (rmax < rmixh)) ? dl1 : dl2;

    du = -(6*uik3 - 125*uik*r2 -5*uik*sk2 + 60*(uik2*r+r3-r*sk2))/uik13;
    irep = 2*ah*rmixh7*term;

    de += mask2*(ah*rmixh7*M_PI*(dl+du)/(30*r2));
    sum += mask2*(irep+idisp);

    energy = sum;
    de *= -AWATER*rI;
    force *= de;
}

/**
 * Compute WCA interaction.
 */
extern "C" __global__ void computeWCAForce(unsigned long long* __restrict__ forceBuffers, mixed* __restrict__ energyBuffer,
        const real4* __restrict__ posq, unsigned int startTileIndex, unsigned int numTileIndices, const float2* __restrict__ radiusEpsilon) {
    unsigned int totalWarps = (blockDim.x*gridDim.x)/TILE_SIZE;
    unsigned int warp = (blockIdx.x*blockDim.x+threadIdx.x)/TILE_SIZE;
    const unsigned int numTiles = numTileIndices;
    unsigned int pos = (unsigned int) (startTileIndex+warp*(long long)numTiles/totalWarps);
    unsigned int end = (unsigned int) (startTileIndex+(warp+1)*(long long)numTiles/totalWarps);
    mixed energy = 0;
    __shared__ AtomData localData[THREAD_BLOCK_SIZE];
    
    do {
        // Extract the coordinates of this tile
        
        const unsigned int tgx = threadIdx.x & (TILE_SIZE-1);
        const unsigned int tbx = threadIdx.x - tgx;
        const unsigned int localGroupIndex = threadIdx.x/TILE_SIZE;
        int x, y;
        AtomData data;
        if (pos < end) {
            y = (int) floor(NUM_BLOCKS+0.5f-SQRT((NUM_BLOCKS+0.5f)*(NUM_BLOCKS+0.5f)-2*pos));
            x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
            if (x < y || x >= NUM_BLOCKS) { // Occasionally happens due to roundoff error.
                y += (x < y ? -1 : 1);
                x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
            }
            unsigned int atom1 = x*TILE_SIZE + tgx;
            loadAtomData(data, atom1, posq, radiusEpsilon);
            loadAtomData(localData[threadIdx.x], y*TILE_SIZE+tgx, posq, radiusEpsilon);
            real emixo, emixh, rmixo, rmixh;
            initParticleParameters(data.radius, data.epsilon, rmixo, rmixh, emixo, emixh);
            data.force = make_real3(0);
            localData[threadIdx.x].force = make_real3(0);

            // Compute forces.

            unsigned int tj = tgx;
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                int atom2 = y*TILE_SIZE+tj;
                if (atom1 != atom2 && atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                    real3 tempForce;
                    real tempEnergy;
                    computeOneInteraction(data, localData[tbx+tj], rmixo, rmixh, emixo, emixh, tempForce, tempEnergy);
                    data.force += tempForce;
                    localData[tbx+tj].force -= tempForce;
                    energy += (x == y ? 0.5f*tempEnergy : tempEnergy);
                    real emjxo, emjxh, rmjxo, rmjxh;
                    initParticleParameters(localData[tbx+tj].radius, localData[tbx+tj].epsilon, rmjxo, rmjxh, emjxo, emjxh);
                    computeOneInteraction(localData[tbx+tj], data, rmjxo, rmjxh, emjxo, emjxh, tempForce, tempEnergy);
                    data.force -= tempForce;
                    localData[tbx+tj].force += tempForce;
                    energy += (x == y ? 0.5f*tempEnergy : tempEnergy);
                }
                tj = (tj+1) & (TILE_SIZE-1);
            }
            unsigned int offset = x*TILE_SIZE + tgx;
            atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long) (data.force.x*0x100000000)));
            atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.y*0x100000000)));
            atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.z*0x100000000)));
            if (x != y) {
                offset = y*TILE_SIZE + tgx;
                atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.x*0x100000000)));
                atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.y*0x100000000)));
                atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.z*0x100000000)));
            }
        }
        pos++;
    } while (pos < end);
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x] -= AWATER*energy;
}
