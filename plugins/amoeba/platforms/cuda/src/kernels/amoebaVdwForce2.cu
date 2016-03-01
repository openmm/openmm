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
    real sigmasum = sigma1_2+sigma2_2;
    real sigma = (sigmasum == 0.0f ? (real) 0 : 2*(sigmaEpsilon1.x*sigma1_2 + sigmaEpsilon2.x*sigma2_2)/(sigma1_2+sigma2_2));
#endif
#if EPSILON_COMBINING_RULE == 1
    real epsilon = 0.5f*(sigmaEpsilon1.y + sigmaEpsilon2.y);
#elif EPSILON_COMBINING_RULE == 2
    real epsilon = SQRT(sigmaEpsilon1.y*sigmaEpsilon2.y);
#elif EPSILON_COMBINING_RULE == 3
    real epssum = sigmaEpsilon1.y+sigmaEpsilon2.y;
    real epsilon = (epssum == 0.0f ? (real) 0 : 2*(sigmaEpsilon1.y*sigmaEpsilon2.y)/(sigmaEpsilon1.y+sigmaEpsilon2.y));
#else
    real epsilon_s = SQRT(sigmaEpsilon1.y) + SQRT(sigmaEpsilon2.y);
    real epsilon = (epsilon_s == 0.0f ? (real) 0 : 4*sigmaEpsilon1.y*sigmaEpsilon2.y/(epsilon_s*epsilon_s));
#endif

    real lambdaI = sigmaEpsilon1.z;
    real lambdaJ = sigmaEpsilon2.z;
    real combindedLambda = 1.0f ;
if (lambdaI != lambdaJ){
      combindedLambda = min(lambdaI,lambdaJ);
}
     real comblambda2= combindedLambda*combindedLambda;
     epsilon = epsilon * comblambda2*comblambda2*combindedLambda;
     real rho = r *  1.0f / (sigma);
     real rho2= rho*rho;
     real rho6= rho2*rho2*rho2;
     real rhoplus= rho+0.07f;
     real rhodec2=(rhoplus)*(rhoplus);
     real rhodec = rhodec2*rhodec2*rhodec2;
     real  scal = 0.7f * (1.0f - combindedLambda)*(1.0f - combindedLambda);
     real s1  = 1.0f / (scal+rhodec*rhoplus);
     real s2 =  1.0f / (scal+rho6*rho + 0.12f );
     real point72= 1.60578147648f ;
     real t1 = point72*s1;
     real t2= 1.12f * s2;
     real t2min= t2 - 2.0f ;
     real dt1= -7.0f * rhodec*t1*s1;
     real dt2 = -7.0f * rho6*t2*s2;
     real termEnergy       = epsilon*t1*(t2min);
     real deltaE         = epsilon*(dt1*(t2- 2.0f )+t1*dt2)* 1.0f /(sigma);
#ifdef USE_CUTOFF
    if (r > TAPER_CUTOFF) {
        real x = r-TAPER_CUTOFF;
        real taper = 1+x*x*x*(TAPER_C3+x*(TAPER_C4+x*TAPER_C5));
        real dtaper = x*x*(3*TAPER_C3+x*(4*TAPER_C4+x*5*TAPER_C5));
        deltaE = termEnergy*dtaper + deltaE*taper;
        termEnergy *= taper;
    }
#endif
    tempEnergy += (includeInteraction ? termEnergy : 0);
    dEdR -= (includeInteraction ? deltaE*invR : 0);

     
