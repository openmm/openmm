c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine hessrot  --  torsional Hessian elements  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "hessrot" computes numerical Hessian elements with respect
c     to torsional angles; either the diagonal or the full matrix
c     can be calculated; the full matrix needs nomega+1 gradient
c     evaluations while the diagonal needs just two evaluations
c
c
      subroutine hessrot (mode,hrot)
      use sizes
      use omega
      use math
      use zcoord
      implicit none
      integer i,j,line
      real*8 e,eps
      real*8, allocatable :: g(:)
      real*8, allocatable :: g0(:)
      real*8, allocatable :: old(:)
      real*8 hrot(nomega,*)
      character*4 mode
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (g(nomega))
      allocate (g0(nomega))
      allocate (old(nomega))
c
c     calculate base values for the torsional gradient
c
      eps = 0.0001d0
      call gradrot (e,g0)
c
c     compute one-sided numerical Hessian from gradient values;
c     set off-diagonal elements to the average symmetric value
c
      if (mode .eq. 'FULL') then
         do i = 1, nomega
            line = zline(i)
            old(i) = ztors(line)
            ztors(line) = ztors(line) + radian*eps
            call makexyz
            call gradrot (e,g)
            ztors(line) = old(i)
            do j = 1, nomega
               hrot(j,i) = (g(j)-g0(j)) / eps
            end do
            do j = 1, i-1
               hrot(j,i) = 0.5d0 * (hrot(j,i)+hrot(i,j))
               hrot(i,j) = hrot(j,i)
            end do
         end do
c
c     compute numerical Hessian diagonal from gradient values
c
      else if (mode .eq. 'DIAG') then
         do i = 1, nomega
            line = zline(i)
            old(i) = ztors(line)
            ztors(line) = ztors(line) + radian*eps
         end do
         call makexyz
         call gradrot (e,g)
         do i = 1, nomega
            hrot(i,i) = (g(i)-g0(i)) / eps
            line = zline(i)
            ztors(line) = old(i)
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (g)
      deallocate (g0)
      deallocate (old)
c
c     restore the Cartesian coordinates to original values
c
      call makexyz
      return
      end
