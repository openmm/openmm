c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1991  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine orthog  --  Gram-Schmidt orthogonalization  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "orthog" performs an orthogonalization of an input matrix
c     via the modified Gram-Schmidt algorithm
c
c     variables and parameters:
c
c     m     first dimension of the matrix to orthogonalize
c     n     second dimension of the matrix to orthogonalize
c     a     matrix to orthogonalize; contains result on exit
c
c
      subroutine orthog (m,n,a)
      implicit none
      integer i,j,k
      integer m,n
      real*8 rkk,rkj
      real*8 a(m,*)
c
c
c     compute the modified Gram-Schmidt orthogonalization
c
      do k = 1, n
         rkk = 0.0d0
         do i = 1, m
            rkk = rkk + a(i,k)**2
         end do
         rkk = sqrt(rkk)
         do i = 1, m
            a(i,k) = a(i,k) / rkk
         end do
         do j = k+1, n
            rkj = 0.0d0
            do i = 1, m
               rkj = rkj + a(i,k)*a(i,j)
            end do
            do i = 1, m
               a(i,j) = a(i,j) - a(i,k)*rkj
            end do
         end do
      end do
      return
      end
