c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1996  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine deflate  --  eigenvalues by method of deflation  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "deflate" uses the power method with deflation to compute the
c     few largest eigenvalues and eigenvectors of a symmetric matrix
c
c     n      dimension of the matrix to be diagonalized
c     nv     number of largest eigenvalues to be extracted
c     a      input with the matrix to be diagonalized; only
c               the lower triangle and diagonal are required
c     ev     returned with the eigenvalues in descending order
c     vec    returned with the eigenvectors of the matrix
c     work   local vector containing temporary work space
c
c
      subroutine deflate (n,nv,a,ev,vec)
      use iounit
      implicit none
      integer i,j,k,n,nv
      integer iter,maxiter
      real*8 random,eps
      real*8 dot1,dot2,ratio
      real*8 ev(*)
      real*8, allocatable :: work(:)
      real*8 a(n,*)
      real*8 vec(n,*)
      external random
c
c
c     initialize number of iterations and convergence criteria
c
      maxiter = 500
      eps = 1.0d-6
c
c     use identity vector as initial guess for eigenvectors
c
      do j = 1, nv
         do i = 1, n
            vec(i,j) = 1.0d0
         end do
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (work(n))
c
c     find the few largest eigenvalues and eigenvectors
c
      do k = 1, nv
         ev(k) = 0.0d0
         dot1 = 0.0d0
         do i = 1, n
            work(i) = 0.0d0
            do j = 1, i-1
               work(i) = work(i) + a(i,j)*vec(j,k)
            end do
            do j = i, n
               work(i) = work(i) + a(j,i)*vec(j,k)
            end do
            dot1 = dot1 + work(i)**2
         end do
c
c     if in or near null space, use random guess as eigenvector
c
         if (dot1 .le. 100.0d0*eps*dble(n)) then
            do i = 1, n
               work(i) = random ()
            end do
         end if
c
c     find the current eigenvalue by iterating to convergence;
c     first multiply vector by matrix and compute dot products
c
         do iter = 1, maxiter
            dot1 = 0.0d0
            dot2 = 0.0d0
            do i = 1, n
               vec(i,k) = 0.0d0
               do j = 1, i-1
                  vec(i,k) = vec(i,k) + a(i,j)*work(j)
               end do
               do j = i, n
                  vec(i,k) = vec(i,k) + a(j,i)*work(j)
               end do
               dot1 = dot1 + vec(i,k)**2
               dot2 = dot2 + vec(i,k)*work(i)
            end do
c
c     normalize new eigenvector and substitute for old one
c
            ratio = abs((ev(k)-dot2) / dot2)
            ev(k) = dot2
            dot1 = sqrt(dot1)
            do i = 1, n
               vec(i,k) = vec(i,k) / dot1
               work(i) = vec(i,k)
            end do
            if (ratio .lt. eps)  goto 20
         end do
         write (iout,10)  k
   10    format (/,' DEFLATE  --  Eigenvalue',i3,' not Fully Converged')
c
c     eliminate the current eigenvalue from the matrix
c
   20    continue
         do i = 1, n
            do j = i, n
               a(j,i) = a(j,i) - ev(k)*vec(i,k)*vec(j,k)
            end do
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (work)
      return
      end
