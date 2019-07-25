c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  function maxwell  --  Maxwell-Boltzmann distribution value  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "maxwell" returns a speed in Angstroms/picosecond randomly
c     selected from a 3-D Maxwell-Boltzmann distribution for the
c     specified particle mass and system temperature
c
c     literature reference:
c
c     P. W. Atkins, "Physical Chemistry, 4th Edition", W. H. Freeman,
c     New York, 1990; see section 24.2 for general discussion
c
c
      function maxwell (mass,temper)
      use units
      implicit none
      real*8 maxwell
      real*8 mass,temper
      real*8 rho,beta
      real*8 random,erfinv
      real*8 xspeed,yspeed
      real*8 zspeed
      external random,erfinv
c
c
c     initialize the speed magnitude of the particle to zero
c
      maxwell = 0.0d0
c
c     set normalization factor for cumulative velocity distribution
c
      if (mass.gt.0.0d0 .and. temper.gt.0.0d0) then
         beta = sqrt(mass / (2.0d0*boltzmann*temper))
c
c     pick a randomly distributed velocity along each of three axes
c
         rho = random ()
         xspeed = erfinv(rho) / beta
         rho = random ()
         yspeed = erfinv(rho) / beta
         rho = random ()
         zspeed = erfinv(rho) / beta
c
c     set the final value of the particle speed in 3-dimensions
c
         maxwell = sqrt(xspeed**2 + yspeed**2 + zspeed**2)
      end if
      return
      end
