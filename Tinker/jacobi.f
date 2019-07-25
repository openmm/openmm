c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine jacobi  --  jacobi matrix diagonalization  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "jacobi" performs a matrix diagonalization of a real
c     symmetric matrix by the method of Jacobi rotations
c
c     variables and parameters:
c
c     n     dimension of the matrix to be diagonalized
c     a     input with the matrix to be diagonalized; only
c              the upper triangle and diagonal are required
c     d     returned with the eigenvalues in ascending order
c     v     returned with the eigenvectors of the matrix
c     b     temporary work vector
c     z     temporary work vector
c
c
      subroutine jacobi (n,a,d,v)
      use iounit
      implicit none
      integer i,j,k
      integer n,ip,iq
      integer nrot,maxrot
      real*8 sm,tresh,s,c,t
      real*8 theta,tau,h,g,p
      real*8 d(*)
      real*8, allocatable :: b(:)
      real*8, allocatable :: z(:)
      real*8 a(n,*)
      real*8 v(n,*)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (b(n))
      allocate (z(n))
c
c     setup and initialization
c
      maxrot = 100
      nrot = 0
      do ip = 1, n
         do iq = 1, n
            v(ip,iq) = 0.0d0
         end do
         v(ip,ip) = 1.0d0
      end do
      do ip = 1, n
         b(ip) = a(ip,ip)
         d(ip) = b(ip)
         z(ip) = 0.0d0
      end do
c
c     perform the jacobi rotations
c
      do i = 1, maxrot
         sm = 0.0d0
         do ip = 1, n-1
            do iq = ip+1, n
               sm = sm + abs(a(ip,iq))
            end do
         end do
         if (sm .eq. 0.0d0)  goto 10
         if (i .lt. 4) then
            tresh = 0.2d0*sm / n**2
         else
            tresh = 0.0d0
         end if
         do ip = 1, n-1
            do iq = ip+1, n
               g = 100.0d0 * abs(a(ip,iq))
               if (i.gt.4 .and. abs(d(ip))+g.eq.abs(d(ip))
     &                    .and. abs(d(iq))+g.eq.abs(d(iq))) then
                  a(ip,iq) = 0.0d0
               else if (abs(a(ip,iq)) .gt. tresh) then
                  h = d(iq) - d(ip)
                  if (abs(h)+g .eq. abs(h)) then
                     t = a(ip,iq) / h
                  else
                     theta = 0.5d0*h / a(ip,iq)
                     t = 1.0d0 / (abs(theta)+sqrt(1.0d0+theta**2))
                     if (theta .lt. 0.0d0)  t = -t
                  end if
                  c = 1.0d0 / sqrt(1.0d0+t**2)
                  s = t * c
                  tau = s / (1.0d0+c)
                  h = t * a(ip,iq)
                  z(ip) = z(ip) - h
                  z(iq) = z(iq) + h
                  d(ip) = d(ip) - h
                  d(iq) = d(iq) + h
                  a(ip,iq) = 0.0d0
                  do j = 1, ip-1
                     g = a(j,ip)
                     h = a(j,iq)
                     a(j,ip) = g - s*(h+g*tau)
                     a(j,iq) = h + s*(g-h*tau)
                  end do
                  do j = ip+1, iq-1
                     g = a(ip,j)
                     h = a(j,iq)
                     a(ip,j) = g - s*(h+g*tau)
                     a(j,iq) = h + s*(g-h*tau)
                  end do
                  do j = iq+1, n
                     g = a(ip,j)
                     h = a(iq,j)
                     a(ip,j) = g - s*(h+g*tau)
                     a(iq,j) = h + s*(g-h*tau)
                  end do
                  do j = 1, n
                     g = v(j,ip)
                     h = v(j,iq)
                     v(j,ip) = g - s*(h+g*tau)
                     v(j,iq) = h + s*(g-h*tau)
                  end do
                  nrot = nrot + 1
               end if
            end do
         end do
         do ip = 1, n
            b(ip) = b(ip) + z(ip)
            d(ip) = b(ip)
            z(ip) = 0.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (b)
      deallocate (z)
c
c     print warning if not converged
c
   10 continue
      if (nrot .eq. maxrot) then
         write (iout,20)
   20    format (/,' JACOBI  --  Matrix Diagonalization not Converged')
      end if
c
c     sort the eigenvalues and vectors
c
      do i = 1, n-1
         k = i
         p = d(i)
         do j = i+1, n
            if (d(j) .lt. p) then
               k = j
               p = d(j)
            end if
         end do
         if (k .ne. i) then
            d(k) = d(i)
            d(i) = p
            do j = 1, n
               p = v(j,i)
               v(j,i) = v(j,k)
               v(j,k) = p
            end do
         end if
      end do
      return
      end
