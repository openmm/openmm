#ifndef AMOEBA_CUDA_VDW_PARTICLE_H
#define AMOEBA_CUDA_VDW_PARTICLE_H

struct Vdw14_7Particle {

    // coordinates, sigma, epsilon

    float x;
    float y;
    float z;

    float sigma;
    float epsilon;

    float force[3];
    float tempForce[3];

};

#endif

