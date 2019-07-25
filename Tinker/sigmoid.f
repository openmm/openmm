c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  function sigmoid  --  general sigmoidal functional form  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "sigmoid" implements a normalized sigmoidal function on the
c     interval [0,1]; the curves connect (0,0) to (1,1) and have
c     a cooperativity controlled by beta, they approach a straight
c     line as beta -> 0 and get more nonlinear as beta increases
c
c
      function sigmoid (beta,x)
      implicit none
      real*8 beta,x
      real*8 sigmoid
      real*8 expmax
      real*8 expmin
      real*8 expterm
c
c
c     compute the value of the normalized sigmoidal function
c
      if (beta .eq. 0.0d0) then
         sigmoid = x
      else
         expmax = 1.0d0 / (exp(-beta) + 1.0d0)
         expmin = 1.0d0 / (exp(beta) + 1.0d0)
         expterm = 1.0d0 / (exp(beta*(2.0d0*x-1.0d0)) + 1.0d0)
         sigmoid = (expmax - expterm) / (expmax - expmin)
      end if
      return
      end
