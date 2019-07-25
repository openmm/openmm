c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine invert  --  gauss-jordan matrix inversion  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "invert" inverts a matrix using the Gauss-Jordan method
c
c     variables and parameters:
c
c     n     dimension of the matrix to be inverted
c     a     matrix to invert; contains inverse on exit
c
c
      subroutine invert (n,a)
      use iounit
      implicit none
      integer i,j,k,n
      integer icol,irow
      integer, allocatable :: ipivot(:)
      integer, allocatable :: indxc(:)
      integer, allocatable :: indxr(:)
      real*8 big,temp
      real*8 pivot
      real*8 a(n,*)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (ipivot(n))
      allocate (indxc(n))
      allocate (indxr(n))
c
c     perform matrix inversion via the Gauss-Jordan algorithm
c
      do i = 1, n
         ipivot(i) = 0
      end do
      do i = 1, n
         big = 0.0d0
         do j = 1, n
            if (ipivot(j) .ne. 1) then
               do k = 1, n
                  if (ipivot(k) .eq. 0) then
                     if (abs(a(j,k)) .ge. big) then
                        big = abs(a(j,k))
                        irow = j
                        icol = k
                     end if
                  else if (ipivot(k) .gt. 1) then
                     write (iout,10)
   10                format (/,' INVERT  --  Cannot Invert',
     &                          ' a Singular Matrix')
                     call fatal
                  end if
               end do
            end if
         end do
         ipivot(icol) = ipivot(icol) + 1
         if (irow .ne. icol) then
            do j = 1, n
               temp = a(irow,j)
               a(irow,j) = a(icol,j)
               a(icol,j) = temp
            end do
         end if
         indxr(i) = irow
         indxc(i) = icol
         if (a(icol,icol) .eq. 0.0d0) then
            write (iout,20)
   20       format (/,' INVERT  --  Cannot Invert a Singular Matrix')
            call fatal
         end if
         pivot = a(icol,icol)
         a(icol,icol) = 1.0d0
         do j = 1, n
            a(icol,j) = a(icol,j) / pivot
         end do
         do j = 1, n
            if (j .ne. icol) then
               temp = a(j,icol)
               a(j,icol) = 0.0d0
               do k = 1, n
                  a(j,k) = a(j,k) - a(icol,k)*temp
               end do
            end if
         end do
      end do
      do i = n, 1, -1
         if (indxr(i) .ne. indxc(i)) then
            do k = 1, n
               temp = a(k,indxr(i))
               a(k,indxr(i)) = a(k,indxc(i))
               a(k,indxc(i)) = temp
            end do
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (ipivot)
      deallocate (indxc)
      deallocate (indxr)
      return
      end
