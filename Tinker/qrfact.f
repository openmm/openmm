c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine qrfact  --  rectangular matrix QR factorization  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "qrfact" computes the QR factorization of an m by n matrix a
c     via Householder transformations with optional column pivoting;
c     the routine determines an orthogonal matrix q, a permutation
c     matrix p, and an upper trapezoidal matrix r with diagonal
c     elements of nonincreasing magnitude, such that a*p = q*r; the
c     Householder transformation for column k, k = 1,2,...,min(m,n),
c     is of the form:
c
c               i - (1/u(k))*u*u(transpose)
c
c     where u has zeros in the first k-1 positions
c
c     arguments and variables :
c
c     n        number of columns in the "a" matrix
c     m        number of rows in the "a" matrix
c     a        on input contains the m by n matrix for which the QR
c                factorization is to be computed; on output the
c                strict upper trapezoidal part contains the strict
c                upper trapezoidal part of r, the lower trapezoidal
c                part contains a factored form of q
c     pivot    logical flag governing use of pivoting
c     ipvt     integer output array which defines the permutation
c                matrix p such that a*p = q*r; column j of p is
c                column ipvt(j) of the identity matrix
c     rdiag    output vector of length n with diagonal elements of r
c
c
      subroutine qrfact (n,m,a,pivot,ipvt,rdiag)
      implicit none
      integer i,j,k
      integer m,n,minmn
      integer jmax,itemp
      integer ipvt(*)
      real*8 aknorm,temp
      real*8 rdiag(*)
      real*8, allocatable :: work(:)
      real*8 a(m,*)
      logical pivot
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (work(n))
c
c     initialize variables, and find the initial column norms
c
      do j = 1, n
         temp = 0.0d0
         do i = 1, m
            temp = temp + a(i,j)**2
         end do
         rdiag(j) = sqrt(temp)
         work(j) = rdiag(j)
         if (pivot)  ipvt(j) = j
      end do
c
c     bring the column of largest norm into the pivot position
c
      minmn = min(m,n)
      do k = 1, minmn
         if (pivot) then
            jmax = k
            do j = k, n
               if (rdiag(j) .gt. rdiag(jmax))  jmax = j
            end do
            if (jmax .ne. k) then
               do i = 1, m
                  temp = a(i,k)
                  a(i,k) = a(i,jmax)
                  a(i,jmax) = temp
               end do
               rdiag(jmax) = rdiag(k)
               work(jmax) = work(k)
               itemp = ipvt(k)
               ipvt(k) = ipvt(jmax)
               ipvt(jmax) = itemp
            end if
         end if
c
c     compute the Householder transformation to reduce the
c     k-th column of "a" to a multiple of the k-th unit vector
c
         aknorm = 0.0d0
         do i = k, m
            aknorm = aknorm + a(i,k)**2
         end do
         aknorm = sqrt(aknorm)
         if (aknorm .ne. 0.0d0) then
            if (a(k,k) .lt. 0.0d0)  aknorm = -aknorm
            do i = k, m
               a(i,k) = a(i,k) / aknorm
            end do
            a(k,k) = a(k,k) + 1.0d0
c
c     apply transform to remaining columns and update column norms
c
            if (n .ge. k+1) then
               do j = k+1, n
                  temp = 0.0d0
                  do i = k, m
                     temp = temp + a(i,k)*a(i,j)
                  end do
                  temp = temp / a(k,k)
                  do i = k, m
                     a(i,j) = a(i,j) - temp*a(i,k)
                  end do
                  if (pivot .and. rdiag(j).ne.0.0d0) then
                     temp = a(k,j) / rdiag(j)
                     if (abs(temp) .lt. 1.0d0) then
                        rdiag(j) = rdiag(j) * sqrt(1.0d0-temp**2)
                     else
                        temp = 0.0d0
                        do i = k+1, m
                           temp = temp + a(i,j)**2
                        end do
                        rdiag(j) = sqrt(temp)
                        work(j) = rdiag(j)
                     end if
                  end if
               end do
            end if
         end if
         rdiag(k) = -aknorm
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (work)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine qrsolve  --  triangular least squares solution  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "qrsolve" solves a*x = b and d*x = 0 in the least squares sense;
c     used with routine "qrfact" to solve least squares problems
c
c     arguments and variables :
c
c     n        number of rows and columns in the matrix r
c     np       leading physical dimension of r in the calling program
c     r        on input, an n by n array with the upper triangular
c                matrix r; on output the full triangle is unaltered,
c                and the strict lower triangle contains the transpose
c                of the strict upper triangular matrix s
c     ipvt     vector of length n which defines the permutation
c                matrix p such that a*p = q*r; column j of p is
c                column ipvt(j) of the identity matrix
c     diag     vector of length n containing the diagonal elements
c                of the matrix d
c     qtb      vector of length n containing the first n elements
c                of the vector q(transpose)*b
c     x        vector of length n containing the least squares
c                solution of the systems a*x = b and d*x = 0
c     sdiag    vector of length n containing the diagonal elements
c                of the upper triangular matrix s
c     xpvt     vector of length n containing permuted (pivoted)
c                solution of the systems
c
c
      subroutine qrsolve (n,np,r,ipvt,diag,qtb,x,sdiag,xpvt)
      implicit none
      integer i,j,k,jj
      integer n,np,nsing
      integer ipvt(*)
      real*8 sine,cosine
      real*8 tangent
      real*8 cotangent
      real*8 qtbpj,temp
      real*8 diag(*)
      real*8 qtb(*)
      real*8 x(*)
      real*8 sdiag(*)
      real*8 xpvt(*)
      real*8 r(np,*)
c
c
c     copy r and (q transpose)*b to preserve input and initialize s;
c     in particular, save the diagonal elements of r in x
c
      do j = 1, n-1
         do k = j+1, n
            r(k,j) = r(j,k)
         end do
      end do
      do j = 1, n
         x(j) = r(j,j)
         xpvt(j) = qtb(j)
      end do
c
c     eliminate the diagonal matrix d using a Givens rotation
c
      do j = 1, n
c
c     prepare the row of d to be eliminated, locating the
c     diagonal element using p from the QR factorization
c
         jj = ipvt(j)
         if (diag(jj) .ne. 0.0d0) then
            do k = j, n
               sdiag(k) = 0.0d0
            end do
            sdiag(j) = diag(jj)
c
c     transform to eliminate the row of d modify only one element
c     of (q transpose)*b beyond the first n, which is initially zero
c
            qtbpj = 0.0d0
            do k = j, n
c
c     determine a Givens rotation which eliminates the
c     appropriate element in the current row of d
c
               if (sdiag(k) .ne. 0.0d0) then
                  if (abs(r(k,k)) .lt. abs(sdiag(k))) then
                     cotangent = r(k,k) / sdiag(k)
                     sine = 0.5d0 / sqrt(0.25d0+0.25d0*cotangent**2)
                     cosine = sine * cotangent
                  else
                     tangent = sdiag(k) / r(k,k)
                     cosine = 0.5d0 / sqrt(0.25d0+0.25d0*tangent**2)
                     sine = cosine * tangent
                  end if
c
c     compute the modified diagonal element of r
c     and the modified element of ((q transpose)*b,0)
c
                  r(k,k) = cosine*r(k,k) + sine*sdiag(k)
                  temp = cosine*xpvt(k) + sine*qtbpj
                  qtbpj = -sine*xpvt(k) + cosine*qtbpj
                  xpvt(k) = temp
c
c     accumulate the tranformation in the row of s
c
                  if (n .ge. k+1) then
                     do i = k+1, n
                        temp = cosine*r(i,k) + sine*sdiag(i)
                        sdiag(i) = -sine*r(i,k) + cosine*sdiag(i)
                        r(i,k) = temp
                     end do
                  end if
               end if
            end do
         end if
c
c     store the diagonal element of s and restore
c     the corresponding diagonal elements of r
c
         sdiag(j) = r(j,j)
         r(j,j) = x(j)
      end do
c
c     solve the triangular system for xpvt; if the system
c     is singular, then obtain a least squares solution
c
      nsing = n
      do j = 1, n
         if (sdiag(j).eq.0.0d0 .and. nsing.eq.n)  nsing = j - 1
         if (nsing .lt. n)  xpvt(j) = 0.0d0
      end do
      if (nsing .ge. 1) then
         do k = 1, nsing
            j = nsing - k + 1
            temp = 0.0d0
            if (nsing .ge. j+1) then
               do i = j+1, nsing
                  temp = temp + r(i,j)*xpvt(i)
               end do
            end if
            xpvt(j) = (xpvt(j)-temp) / sdiag(j)
         end do
      end if
c
c     permute the components of xpvt back to components of x
c
      do j = 1, n
         k = ipvt(j)
         x(k) = xpvt(j)
      end do
      return
      end
