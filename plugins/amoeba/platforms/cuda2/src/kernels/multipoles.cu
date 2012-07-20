extern "C" __global__ void computeLabFrameMoments(real4* __restrict__ posq, int4* __restrict__ multipoleParticles, float* __restrict__ molecularDipoles,
        float* __restrict__ molecularQuadrupoles, real* __restrict__ labFrameDipoles, real* __restrict__ labFrameQuadrupoles) {
    // get coordinates of this atom and the z & x axis atoms
    // compute the vector between the atoms and 1/sqrt(d2), d2 is distance between
    // this atom and the axis atom
 
    // this atom is referred to as the k-atom in notes below
 
    // code common to ZThenX and Bisector
    
    for (int particleIndex = blockIdx.x*blockDim.x+threadIdx.x; particleIndex < NUM_ATOMS; particleIndex += gridDim.x*blockDim.x) {
        int4 particles = multipoleParticles[particleIndex];
        if (particles.x >= 0 && particles.z >= 0) {
            real4 thisParticlePos = posq[particleIndex];
            real4 posZ = posq[particles.z];
            real3 vectorZ = make_real3(posZ.x-thisParticlePos.x, posZ.y-thisParticlePos.y, posZ.z-thisParticlePos.z);
            real4 posX = posq[particles.x];
            real3 vectorX = make_real3(posX.x-thisParticlePos.x, posX.y-thisParticlePos.y, posX.z-thisParticlePos.z);
            int axisType = particles.w; 
    
            /*
                z-only
                   (1) norm z
                   (2) select random x
                   (3) x = x - (x.z)z
                   (4) norm x
        
                z-then-x
                   (1) norm z
                   (2) norm x (not needed)
                   (3) x = x - (x.z)z
                   (4) norm x
        
                bisector
                   (1) norm z
                   (2) norm x 
                   (3) z = x + z
                   (4) norm z
                   (5) x = x - (x.z)z 
                   (6) norm x 
        
                z-bisect
                   (1) norm z
                   (2) norm x 
                   (3) norm y 
                   (3) x = x + y
                   (4) norm x
                   (5) x = x - (x.z)z 
                   (6) norm x 
        
                3-fold
                   (1) norm z
                   (2) norm x 
                   (3) norm y 
                   (4) z = x + y + z
                   (5) norm z
                   (6) x = x - (x.z)z 
                   (7) norm x 
        
            */
        
            // branch based on axis type
             
            vectorZ = normalize(vectorZ);
        
            if (axisType == 1) {
        
                // bisector
                
                vectorX = normalize(vectorX);
                vectorZ += vectorX;
                vectorZ = normalize(vectorZ);
            }
            else if (axisType == 2 || axisType == 3) { 
         
                // z-bisect
        
                if (particles.y >= 0 && particles.y < NUM_ATOMS) {
                    real4 posY = posq[particles.y];
                    real3 vectorY = make_real3(posY.x-thisParticlePos.x, posY.y-thisParticlePos.y, posY.z-thisParticlePos.z);
                    vectorY = normalize(vectorY);
                    vectorX = normalize(vectorX);
                    if (axisType == 2) {
                        vectorX += vectorY;
                        vectorX = normalize(vectorX);
                    }
                    else { 
             
                        // 3-fold
                
                        vectorZ += vectorX + vectorY;
                        vectorZ = normalize(vectorZ);
                    }
                }
         
            }
            else if (axisType >= 4)
                vectorX = make_real3((real) 0.1f);
            
            // x = x - (x.z)z
        
            vectorX -= dot(vectorZ, vectorX)*vectorZ;
            vectorX = normalize(vectorX);
            real3 vectorY = cross(vectorZ, vectorX);
         
            // use identity rotation matrix for unrecognized axis types
        
            if (axisType < 0 || axisType > 4) {
        
                vectorX.x = 1;
                vectorX.y = 0;
                vectorX.z = 0;
        
                vectorY.x = 0;
                vectorY.y = 1;
                vectorY.z = 0;
        
                vectorZ.x = 0;
                vectorZ.y = 0;
                vectorZ.z = 1;
            }
            
            // Check the chirality and see whether it needs to be reversed
            
            bool reverse = false;
            if (axisType != 0 && particles.x >= 0 && particles.y >=0 && particles.z >= 0) {
                real4 posY = posq[particles.y];
                real delta[4][3];

                delta[0][0] = thisParticlePos.x - posY.x;
                delta[0][1] = thisParticlePos.y - posY.y;
                delta[0][2] = thisParticlePos.z - posY.z;

                delta[1][0] = posZ.x - posY.x;
                delta[1][1] = posZ.y - posY.y;
                delta[1][2] = posZ.z - posY.z;

                delta[2][0] = posX.x - posY.x;
                delta[2][1] = posX.y - posY.y;
                delta[2][2] = posX.z - posY.z;

                delta[3][0] = delta[1][1]*delta[2][2] - delta[1][2]*delta[2][1];
                delta[3][1] = delta[2][1]*delta[0][2] - delta[2][2]*delta[0][1];
                delta[3][2] = delta[0][1]*delta[1][2] - delta[0][2]*delta[1][1];

                real volume = delta[3][0]*delta[0][0] + delta[3][1]*delta[1][0] + delta[3][2]*delta[2][0];
                reverse = (volume < 0);
            }
        
            // Transform the dipole
            
            unsigned int offset = 3*particleIndex;
            real molDipole[3];
            molDipole[0] = molecularDipoles[offset];
            molDipole[1] = molecularDipoles[offset+1];
            molDipole[2] = molecularDipoles[offset+2];
            if (reverse)
                molDipole[1] *= -1;
            labFrameDipoles[offset] = molDipole[0]*vectorX.x + molDipole[1]*vectorY.x + molDipole[2]*vectorZ.x;
            labFrameDipoles[offset+1] = molDipole[0]*vectorX.y + molDipole[1]*vectorY.y + molDipole[2]*vectorZ.y;
            labFrameDipoles[offset+2] = molDipole[0]*vectorX.z + molDipole[1]*vectorY.z + molDipole[2]*vectorZ.z;
            
            // ---------------------------------------------------------------------------------------
            
            // Transform the quadrupole
            
            real mPole[3][3];
            offset = 9*particleIndex;
            mPole[0][0] = molecularQuadrupoles[offset];
            mPole[0][1] = molecularQuadrupoles[offset+1];
            mPole[0][2] = molecularQuadrupoles[offset+2];
        
            mPole[1][0] = molecularQuadrupoles[offset+3];
            mPole[1][1] = molecularQuadrupoles[offset+4];
            mPole[1][2] = molecularQuadrupoles[offset+5];
        
            mPole[2][0] = molecularQuadrupoles[offset+6];
            mPole[2][1] = molecularQuadrupoles[offset+7];
            mPole[2][2] = molecularQuadrupoles[offset+8];
        
            if (reverse) {
                mPole[0][1] *= -1;
                mPole[1][0] *= -1;
                mPole[1][2] *= -1;
                mPole[2][1] *= -1;
            }
            
            labFrameQuadrupoles[offset+8] = vectorX.z*(vectorX.z*mPole[0][0] + vectorY.z*mPole[0][1] + vectorZ.z*mPole[0][2]);
            labFrameQuadrupoles[offset+8] += vectorY.z*(vectorX.z*mPole[1][0] + vectorY.z*mPole[1][1] + vectorZ.z*mPole[1][2]);
            labFrameQuadrupoles[offset+8] += vectorZ.z*(vectorX.z*mPole[2][0] + vectorY.z*mPole[2][1] + vectorZ.z*mPole[2][2]);
    
            labFrameQuadrupoles[offset+4] = vectorX.y*(vectorX.y*mPole[0][0] + vectorY.y*mPole[0][1] + vectorZ.y*mPole[0][2]);
            labFrameQuadrupoles[offset+4] += vectorY.y*(vectorX.y*mPole[1][0] + vectorY.y*mPole[1][1] + vectorZ.y*mPole[1][2]);
            labFrameQuadrupoles[offset+4] += vectorZ.y*(vectorX.y*mPole[2][0] + vectorY.y*mPole[2][1] + vectorZ.y*mPole[2][2]);
    
            labFrameQuadrupoles[offset+5] = vectorX.y*(vectorX.z*mPole[0][0] + vectorY.z*mPole[0][1] + vectorZ.z*mPole[0][2]);
            labFrameQuadrupoles[offset+5] += vectorY.y*(vectorX.z*mPole[1][0] + vectorY.z*mPole[1][1] + vectorZ.z*mPole[1][2]);
            labFrameQuadrupoles[offset+5] += vectorZ.y*(vectorX.z*mPole[2][0] + vectorY.z*mPole[2][1] + vectorZ.z*mPole[2][2]);
    
            labFrameQuadrupoles[offset] = vectorX.x*(vectorX.x*mPole[0][0] + vectorY.x*mPole[0][1] + vectorZ.x*mPole[0][2]);
            labFrameQuadrupoles[offset] += vectorY.x*(vectorX.x*mPole[1][0] + vectorY.x*mPole[1][1] + vectorZ.x*mPole[1][2]);
            labFrameQuadrupoles[offset] += vectorZ.x*(vectorX.x*mPole[2][0] + vectorY.x*mPole[2][1] + vectorZ.x*mPole[2][2]);
    
            labFrameQuadrupoles[offset+1] = vectorX.x*(vectorX.y*mPole[0][0] + vectorY.y*mPole[0][1] + vectorZ.y*mPole[0][2]);
            labFrameQuadrupoles[offset+1] += vectorY.x*(vectorX.y*mPole[1][0] + vectorY.y*mPole[1][1] + vectorZ.y*mPole[1][2]);
            labFrameQuadrupoles[offset+1] += vectorZ.x*(vectorX.y*mPole[2][0] + vectorY.y*mPole[2][1] + vectorZ.y*mPole[2][2]);
    
            labFrameQuadrupoles[offset+2] = vectorX.x*(vectorX.z*mPole[0][0] + vectorY.z*mPole[0][1] + vectorZ.z*mPole[0][2]);
            labFrameQuadrupoles[offset+2] += vectorY.x*(vectorX.z*mPole[1][0] + vectorY.z*mPole[1][1] + vectorZ.z*mPole[1][2]);
            labFrameQuadrupoles[offset+2] += vectorZ.x*(vectorX.z*mPole[2][0] + vectorY.z*mPole[2][1] + vectorZ.z*mPole[2][2]);
     
            labFrameQuadrupoles[offset+3] = labFrameQuadrupoles[offset+1];
            labFrameQuadrupoles[offset+6] = labFrameQuadrupoles[offset+2];
            labFrameQuadrupoles[offset+7] = labFrameQuadrupoles[offset+5];
        }
    }
}

