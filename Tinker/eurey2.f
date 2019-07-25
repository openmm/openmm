c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine eurey2  --  atom-by-atom Urey-Bradley Hessian  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "eurey2" calculates second derivatives of the Urey-Bradley
c     interaction energy for a single atom at a time
c
c
      subroutine eurey2 (i)
      use sizes
      use atoms
      use bound
      use couple
      use group
      use hessn
      use urey
      use urypot
      implicit none
      integer i,j,ia,ic,iurey
      real*8 ideal,force,fgrp
      real*8 xac,yac,zac
      real*8 rac,rac2
      real*8 dt,dt2,term
      real*8 termx,termy,termz
      real*8 de,deddt,d2eddt2
      real*8 d2e(3,3)
      logical proceed
c
c
c     calculate the Urey-Bradley interaction Hessian elements
c
      do iurey = 1, nurey
         ia = iury(1,iurey)
         ic = iury(3,iurey)
         ideal = ul(iurey)
         force = uk(iurey)
c
c     decide whether to compute the current interaction
c
         proceed = (i.eq.ia .or. i.eq.ic)
         if (proceed .and. use_group)
     &      call groups (proceed,fgrp,ia,ic,0,0,0,0)
c
c     compute the value of the 1-3 distance deviation
c
         if (proceed) then
            if (i .eq. ic) then
               ic = ia
               ia = i
            end if
            xac = x(ia) - x(ic)
            yac = y(ia) - y(ic)
            zac = z(ia) - z(ic)
            if (use_polymer)  call image (xac,yac,zac)
            rac2 = xac*xac + yac*yac + zac*zac
            rac = sqrt(rac2)
            dt = rac - ideal
            dt2 = dt * dt
            deddt = 2.0d0 * ureyunit * force * dt
     &                 * (1.0d0+1.5d0*cury*dt+2.0d0*qury*dt2)
            d2eddt2 = 2.0d0 * ureyunit * force
     &                   * (1.0d0+3.0d0*cury*dt+6.0d0*qury*dt2)
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
            de = deddt / rac
            term = (d2eddt2-de) / rac2
            termx = term * xac
            termy = term * yac
            termz = term * zac
            d2e(1,1) = termx*xac + de
            d2e(1,2) = termx*yac
            d2e(1,3) = termx*zac
            d2e(2,1) = d2e(1,2)
            d2e(2,2) = termy*yac + de
            d2e(2,3) = termy*zac
            d2e(3,1) = d2e(1,3)
            d2e(3,2) = d2e(2,3)
            d2e(3,3) = termz*zac + de
c
c     increment diagonal and non-diagonal Hessian elements
c
            do j = 1, 3
               hessx(j,ia) = hessx(j,ia) + d2e(1,j)
               hessy(j,ia) = hessy(j,ia) + d2e(2,j)
               hessz(j,ia) = hessz(j,ia) + d2e(3,j)
               hessx(j,ic) = hessx(j,ic) - d2e(1,j)
               hessy(j,ic) = hessy(j,ic) - d2e(2,j)
               hessz(j,ic) = hessz(j,ic) - d2e(3,j)
            end do
         end if
      end do
      return
      end
