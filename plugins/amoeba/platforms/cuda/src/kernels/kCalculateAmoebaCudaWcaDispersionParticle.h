#ifndef AMOEBA_CUDA_WCA_DISPERSION_PARTICLE_H
#define AMOEBA_CUDA_WCA_DISPERSION_PARTICLE_H

struct WcaDispersionParticle {

    // coordinates

    float x;
    float y;
    float z;

    // radius &  epsilon

    float radius;
    float epsilon;

    float force[3];
    float padding;

};

#endif

