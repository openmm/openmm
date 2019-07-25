c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine hessrgd  --  rigid body Hessian elements  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "hessrgd" computes the numerical Hessian elements with
c     respect to rigid body coordinates via 6*ngroup+1 gradient
c     evaluations
c
c
      subroutine hessrgd (hrigid)
      use sizes
      use atoms
      use group
      use rigid
      implicit none
      integer i,j,k,m,nvar
      real*8 e,eps,old
      real*8 hrigid(6*ngrp,*)
      real*8, allocatable :: g(:,:)
      real*8, allocatable :: g0(:,:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (g(6,ngrp))
      allocate (g0(6,ngrp))
c
c     calculate base values for the rigid body gradient
c
      eps = 0.00001d0
      call gradrgd (e,g0)
c
c     compute one-sided numerical Hessian from gradient values;
c     set off-diagonal elements to the average symmetric value
c
      nvar = 6 * ngrp
      do i = 1, nvar
         j = (i-1)/6 + 1
         k = mod(i-1,6) + 1
         old = rbc(k,j)
         rbc(k,j) = rbc(k,j) + eps
         call rigidxyz
         call gradrgd (e,g)
         rbc(k,j) = old
         do m = 1, nvar
            j = (m-1)/6 + 1
            k = mod(m-1,6) + 1
            hrigid(m,i) = (g(k,j)-g0(k,j)) / eps
         end do
         do m = 1, i-1
            hrigid(m,i) = 0.5d0 * (hrigid(m,i)+hrigid(i,m))
            hrigid(i,m) = hrigid(m,i)
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (g)
      deallocate (g0)
c
c     restore the Cartesian coordinates to original values
c
      call rigidxyz
      return
      end
