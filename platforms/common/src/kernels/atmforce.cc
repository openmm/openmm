KERNEL void HybridForce(int numParticles,
			int paddedNumParticles,
			GLOBAL mm_long* RESTRICT force,
			GLOBAL mm_long* RESTRICT force_state1,
			GLOBAL mm_long* RESTRICT force_state2,
			float sp
			){
  int i = GLOBAL_ID;
  real lmb  = sp;
  real lmb1 = 1.0f - sp;
  while (i < numParticles) {
    force[i] += (mm_long) (lmb*force_state2[i] + lmb1*force_state1[i]);
    force[i+paddedNumParticles  ] += (mm_long) (lmb*force_state2[i+paddedNumParticles  ]+ lmb1*force_state1[i+paddedNumParticles  ]);
    force[i+paddedNumParticles*2] += (mm_long) (lmb*force_state2[i+paddedNumParticles*2]+ lmb1*force_state1[i+paddedNumParticles*2]);
    i += GLOBAL_SIZE;
  }
}

KERNEL void CopyState(int numParticles,
			GLOBAL real4* RESTRICT posq,
			GLOBAL real4* RESTRICT posq1,
			GLOBAL real4* RESTRICT posq2,
			GLOBAL float4* RESTRICT displ
#ifdef USE_MIXED_PRECISION
			,
			GLOBAL real4* RESTRICT posqCorrection,
			GLOBAL real4* RESTRICT posq1Correction,
			GLOBAL real4* RESTRICT posq2Correction
#endif
			){

  //set the coordinates of the context for state 1
  int i = GLOBAL_ID;
  while (i < numParticles) {
    posq1[i] = posq[i];
#ifdef USE_MIXED_PRECISION
    posq1Correction[i] = posqCorrection[i];
#endif
    i += GLOBAL_SIZE;
  }
  
  //set the coordinates of the context for state 2
  i = GLOBAL_ID;
  while (i < numParticles) {
    real4 d = make_real4((real)displ[i].x, (real)displ[i].y, (real)displ[i].z, 0);
    posq2[i] = posq[i] + d;
#ifdef USE_MIXED_PRECISION
    posq2Correction[i] = posqCorrection[i];
#endif
    i += GLOBAL_SIZE;
  }
}


