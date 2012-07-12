{
#ifdef USE_CUTOFF
    unsigned int includeInteraction = (!isExcluded && r2 < CUTOFF_SQUARED);
#else
    unsigned int includeInteraction = (!isExcluded);
#endif
    real tempForce = 0.0f;
#if SIGMA_COMBINING_RULE == 1
    real sigma = sigmaEpsilon1.x + sigmaEpsilon2.x;
#elif SIGMA_COMBINING_RULE == 2
    real sigma = 2*SQRT(sigmaEpsilon1.x*sigmaEpsilon2.x);
#else
    real sigma1_2 = sigmaEpsilon1.x*sigmaEpsilon1.x;
    real sigma2_2 = sigmaEpsilon2.x*sigmaEpsilon2.x;
    real sigma = 2*(sigmaEpsilon1.x*sigma1_2 + sigmaEpsilon2.x*sigma2_2)/(sigma1_2+sigma2_2);
#endif
#if EPSILON_COMBINING_RULE == 1
    real epsilon = sigmaEpsilon1.y + sigmaEpsilon2.y;
#elif EPSILON_COMBINING_RULE == 2
    real epsilon = 2*SQRT(sigmaEpsilon1.y*sigmaEpsilon2.y);
#elif EPSILON_COMBINING_RULE == 3
    real epsilon1_2 = sigmaEpsilon1.x*sigmaEpsilon1.x;
    real epsilon2_2 = sigmaEpsilon2.x*sigmaEpsilon2.x;
    real epsilon = 2*(sigmaEpsilon1.x*epsilon1_2 + sigmaEpsilon2.x*epsilon2_2)/(epsilon1_2+epsilon2_2);
#else
    real epsilon_s = SQRT(sigmaEpsilon1.y) + SQRT(sigmaEpsilon2.y);
    real epsilon = 4*sigmaEpsilon1.y*sigmaEpsilon2.y/(epsilon_s*epsilon_s);
#endif
    real r6 = r2*r2*r2;
    real r7 = r6*r;
    real sigma7 = sigma*sigma;
    sigma7 = sigma7*sigma7*sigma7*sigma;
    real rho = r7 + sigma7*0.12f;
    real invRho = RECIP(rho);
    real tau = 1.07f/(r + 0.07f*sigma);
    real tau7 = tau*tau*tau;
    tau7 = tau7*tau7*tau;
    real dTau = tau/1.07f;
    real tmp = sigma7*invRho;
    real gTau = epsilon*tau7*r6*1.12f*tmp*tmp;
    real termEnergy = epsilon*sigma7*tau7*((sigma7*1.12f*invRho)-2.0f);
    real deltaE = (-7.0f*(dTau*termEnergy+gTau))*invR;
    if (r > TAPER_CUTOFF) {
        real taper = TAPER_C0+r*(TAPER_C1+r*(TAPER_C2+r*(TAPER_C3+r*(TAPER_C4+r*TAPER_C5))));
        real dtaper = TAPER_C1+r*(2*TAPER_C2+r*(3*TAPER_C3+r*(4*TAPER_C4+r*5*TAPER_C5)));
        deltaE = termEnergy*dtaper + deltaE*taper;
        termEnergy *= taper;
    }
    tempEnergy += (includeInteraction ? termEnergy : 0);
    dEdR -= (includeInteraction ? deltaE : 0);
}
