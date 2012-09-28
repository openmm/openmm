#ifndef AMOEBA_CUDA_KIRKWOOD_PARTICLE_H
#define AMOEBA_CUDA_KIRKWOOD_PARTICLE_H

struct KirkwoodParticle {

    // coordinates charge

    float x;
    float y;
    float z;
    float q;

    // lab frame dipole

    float labFrameDipole[3]; 

    // lab frame quadrupole

    float labFrameQuadrupole_XX;
    float labFrameQuadrupole_XY;
    float labFrameQuadrupole_XZ;
    float labFrameQuadrupole_YY;
    float labFrameQuadrupole_YZ;
    float labFrameQuadrupole_ZZ;

    // induced dipole

    float inducedDipole[3]; 

    // polar induced dipole

    float inducedDipoleP[3];

    // Born radii

    float bornRadius;

    float force[3];
#ifdef INCLUDE_TORQUE
    float torque[3];
#endif

    float dBornRadius;
    float dBornRadiusPolar;
    //float padding;

};

#endif
