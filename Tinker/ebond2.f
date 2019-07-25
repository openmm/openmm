c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ebond2  --  atom-by-atom bond stretch Hessian  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ebond2" calculates second derivatives of the bond
c     stretching energy for a single atom at a time
c
c
      subroutine ebond2 (i)
      use sizes
      use atmlst
      use atoms
      use bndpot
      use bndstr
      use bound
      use couple
      use group
      use hessn
      implicit none
      integer i,j,k,ia,ib
      real*8 ideal,force,fgrp
      real*8 xab,yab,zab
      real*8 rab,rab2
      real*8 expterm,bde
      real*8 dt,dt2,term
      real*8 termx,termy,termz
      real*8 de,deddt,d2eddt2
      real*8 d2e(3,3)
      logical proceed
c
c
c     calculate the bond stretch interaction Hessian elements
c
      ia = i
      do k = 1, n12(ia)
         j = bndlist(k,ia)
         if (ibnd(1,j) .eq. ia) then
            ib = ibnd(2,j)
         else
            ib = ibnd(1,j)
         end if
         ideal = bl(j)
         force = bk(j)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,ia,ib,0,0,0,0)
c
c     compute the value of the bond length deviation
c
         if (proceed) then
            xab = x(ia) - x(ib)
            yab = y(ia) - y(ib)
            zab = z(ia) - z(ib)
            if (use_polymer)  call image (xab,yab,zab)
            rab2 = xab*xab + yab*yab + zab*zab
            rab = sqrt(rab2)
            dt = rab - ideal
c
c     harmonic potential uses Taylor expansion of Morse potential
c     through the fourth power of the bond length deviation
c
            if (bndtyp .eq. 'HARMONIC') then
               dt2 = dt * dt
               deddt = 2.0d0 * bndunit * force * dt
     &                    * (1.0d0+1.5d0*cbnd*dt+2.0d0*qbnd*dt2)
               d2eddt2 = 2.0d0 * bndunit * force
     &                      * (1.0d0+3.0d0*cbnd*dt+6.0d0*qbnd*dt2)
c
c     Morse potential uses energy = BDE * (1 - e**(-alpha*dt))**2)
c     with the approximations alpha = sqrt(ForceConst/BDE) = -2
c     and BDE = Bond Dissociation Energy = ForceConst/alpha**2
c
            else if (bndtyp .eq. 'MORSE') then
               expterm = exp(-2.0d0*dt)
               bde = 0.25d0 * bndunit * force
               deddt = 4.0d0 * bde * (1.0d0-expterm) * expterm
               d2eddt2 = -8.0d0 * bde * (1.0d0-2.0d0*expterm) * expterm
            end if
c
c     scale the interaction based on its group membership
c
            if (use_group) then
               deddt = deddt * fgrp
               d2eddt2 = d2eddt2 * fgrp
            end if
c
c     set the chain rule terms for the Hessian elements
c
            if (rab2 .eq. 0.0d0) then
               de = 0.0d0
               term = 0.0d0
            else
               de = deddt / rab
               term = (d2eddt2-de) / rab2
            end if
            termx = term * xab
            termy = term * yab
            termz = term * zab
            d2e(1,1) = termx*xab + de
            d2e(1,2) = termx*yab
            d2e(1,3) = termx*zab
            d2e(2,1) = d2e(1,2)
            d2e(2,2) = termy*yab + de
            d2e(2,3) = termy*zab
            d2e(3,1) = d2e(1,3)
            d2e(3,2) = d2e(2,3)
            d2e(3,3) = termz*zab + de
c
c     increment diagonal and non-diagonal Hessian elements
c
            do j = 1, 3
               hessx(j,ia) = hessx(j,ia) + d2e(1,j)
               hessy(j,ia) = hessy(j,ia) + d2e(2,j)
               hessz(j,ia) = hessz(j,ia) + d2e(3,j)
               hessx(j,ib) = hessx(j,ib) - d2e(1,j)
               hessy(j,ib) = hessy(j,ib) - d2e(2,j)
               hessz(j,ib) = hessz(j,ib) - d2e(3,j)
            end do
         end if
      end do
      return
      end
