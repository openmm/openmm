c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine ebond  --  bond stretch potential energy  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "ebond" calculates the bond stretching energy
c
c
      subroutine ebond
      use sizes
      use atoms
      use bndpot
      use bndstr
      use bound
      use energi
      use group
      use usage
      implicit none
      integer i,ia,ib
      real*8 e,ideal,force
      real*8 expterm,bde
      real*8 dt,dt2,fgrp
      real*8 xab,yab,zab,rab
      logical proceed
c
c
c     zero out the bond stretching energy
c
      eb = 0.0d0
      if (nbond .eq. 0)  return
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(nbond,ibnd,bl,bk,use,
!$OMP& x,y,z,cbnd,qbnd,bndtyp,bndunit,use_group,use_polymer)
!$OMP& shared(eb)
!$OMP DO reduction(+:eb) schedule(guided)
c
c     calculate the bond stretching energy term
c
      do i = 1, nbond
         ia = ibnd(1,i)
         ib = ibnd(2,i)
         ideal = bl(i)
         force = bk(i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,ia,ib,0,0,0,0)
         if (proceed)  proceed = (use(ia) .or. use(ib))
c
c     compute the value of the bond length deviation
c
         if (proceed) then
            xab = x(ia) - x(ib)
            yab = y(ia) - y(ib)
            zab = z(ia) - z(ib)
            if (use_polymer)  call image (xab,yab,zab)
            rab = sqrt(xab*xab + yab*yab + zab*zab)
            dt = rab - ideal
c
c     harmonic potential uses Taylor expansion of Morse potential
c     through the fourth power of the bond length deviation
c
            if (bndtyp .eq. 'HARMONIC') then
               dt2 = dt * dt
               e = bndunit * force * dt2 * (1.0d0+cbnd*dt+qbnd*dt2)
c
c     Morse potential uses energy = BDE * (1 - e**(-alpha*dt))**2)
c     with the approximations alpha = sqrt(ForceConst/BDE) = -2
c     and BDE = Bond Dissociation Energy = ForceConst/alpha**2
c
            else if (bndtyp .eq. 'MORSE') then
               expterm = exp(-2.0d0*dt)
               bde = 0.25d0 * bndunit * force
               e = bde * (1.0d0-expterm)**2
            end if
c
c     scale the interaction based on its group membership
c
            if (use_group)  e = e * fgrp
c
c     increment the total bond stretching energy
c
            eb = eb + e
         end if
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end
