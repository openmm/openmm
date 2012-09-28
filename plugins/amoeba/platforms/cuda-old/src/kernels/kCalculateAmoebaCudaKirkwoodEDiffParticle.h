#ifndef AMOEBA_CUDA_KIRKWOOD_PARTICLE_H
#define AMOEBA_CUDA_KIRKWOOD_PARTICLE_H

struct KirkwoodEDiffParticle {

    // coordinates charge

    float x;
    float y;
    float z;
    float q;

    // scaling factor

    float thole;
    float damp;
    
    // lab frame dipole

    float labFrameDipole[3]; 

    // lab frame quadrupole

    float labFrameQuadrupole_XX;
    float labFrameQuadrupole_XY;
    float labFrameQuadrupole_XZ;
    float labFrameQuadrupole_YY;
    float labFrameQuadrupole_YZ;
    float labFrameQuadrupole_ZZ;

    // induced dipole and polar counterpart

    float inducedDipole[3]; 
    float inducedDipoleP[3]; 

    // solvent induced dipole and polar counterpart

    float inducedDipoleS[3]; 
    float inducedDipolePS[3]; 

    // Born radii

    float force[3];
//    float torque[3];

};

#endif
