{
    float sig = sigmaEpsilon1.x + sigmaEpsilon2.x;
    float sig2 = invR*sig;
    sig2 *= sig2;
    float sig6 = sig2*sig2*sig2;
    float eps = sigmaEpsilon1.y*sigmaEpsilon2.y;
    float tempForce = eps*(12.0f*sig6 - 6.0f)*sig6;
    tempEnergy += eps*(sig6 - 1.0f)*sig6;
    const float EpsilonFactor = 138.935485f;
#ifdef USE_CUTOFF
    tempForce += EpsilonFactor*posq1.w*posq2.w*(invR - 2.0f*REACTION_FIELD_K*r2);
    tempEnergy += EpsilonFactor*posq1.w*posq2.w*(invR + REACTION_FIELD_K*r2 - REACTION_FIELD_C);
#else
    tempForce += EpsilonFactor*posq1.w*posq2.w*invR;
    tempEnergy += EpsilonFactor*posq1.w*posq2.w*invR;
#endif
    dEdR += tempForce*invR*invR;
}
