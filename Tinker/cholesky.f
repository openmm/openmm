c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine cholesky  --  modified Cholesky linear solver  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "cholesky" uses a modified Cholesky method to solve the linear
c     system Ax = b, returning "x" in "b"; "A" is a real symmetric
c     positive definite matrix with its upper triangle (including the
c     diagonal) stored by rows
c
c     literature reference:
c
c     R. S. Martin, G. Peters and J. H. Wilkinson, "Symmetric
c     Decomposition of a Positive Definite Matrix", Numerische
c     Mathematik, 7, 362-383 (1965)
c
c
      subroutine cholesky (nvar,a,b)
      implicit none
      integer i,j,k,nvar
      integer ii,ij,ik,ki,kk
      integer im,jk,jm
      real*8 r,s,t
      real*8 a(*)
      real*8 b(*)
c
c
c     Cholesky factorization to reduce "A" to (L)(D)(L transpose)
c     "L" has a unit diagonal; store 1.0/D on the diagonal of "A"
c
      ii = 1
      do i = 1, nvar
         im = i - 1
         if (i .ne. 1) then
            ij = i
            do j = 1, im
               r = a(ij)
               if (j .ne. 1) then
                  ik = i
                  jk = j
                  jm = j - 1
                  do k = 1, jm
                     r = r - a(ik)*a(jk)
                     ik = nvar - k + ik
                     jk = nvar - k + jk
                  end do
               end if
               a(ij) = r
               ij = nvar - j + ij
            end do
         end if
         r = a(ii)
         if (i .ne. 1) then
            kk = 1
            ik = i
            do k = 1, im
               s = a(ik)
               t = s * a(kk)
               a(ik) = t
               r = r - s*t
               ik = nvar - k + ik
               kk = nvar - k + 1 + kk
            end do
         end if
         a(ii) = 1.0d0 / r
         ii = nvar - i + 1 + ii
      end do
c
c     solve linear equations; first solve Ly = b for y
c
      do i = 1, nvar
         if (i .ne. 1) then
            ik = i
            im = i - 1
            r = b(i)
            do k = 1, im
               r = r - b(k)*a(ik)
               ik = nvar - k + ik
            end do
            b(i) = r
         end if
      end do
c
c     finally, solve (D)(L transpose)(x) = y for x
c
      ii = nvar*(nvar+1)/2
      do j = 1, nvar
         i = nvar + 1 - j
         r = b(i) * a(ii)
         if (j .ne. 1) then
            im = i + 1
            ki = ii + 1
            do k = im, nvar
               r = r - a(ki)*b(k)
               ki = ki + 1
            end do
         end if
         b(i) = r
         ii = ii - j - 1
      end do
      return
      end
