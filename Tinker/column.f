c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine column  --  access Hessian elements by column  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "column" takes the off-diagonal Hessian elements stored
c     as sparse rows and sets up indices to allow column access
c
c
      subroutine column (nvar,hinit,hstop,hindex,
     &                   cinit,cstop,cindex,cvalue)
      implicit none
      integer i,j,k
      integer m,nvar
      integer hinit(*)
      integer hstop(*)
      integer cinit(*)
      integer cstop(*)
      integer hindex(*)
      integer cindex(*)
      integer cvalue(*)
c
c
c     zero out the start and end marker for each column
c
      do i = 1, nvar
         cinit(i) = 0
         cstop(i) = 0
      end do
c
c     count the number of elements in each column
c
      do i = 1, nvar
         do j = hinit(i), hstop(i)
            k = hindex(j)
            cstop(k) = cstop(k) + 1
         end do
      end do
c
c     set each start marker just past last element for its column
c
      cinit(1) = cstop(1) + 1
      do i = 2, nvar
         cinit(i) = cinit(i-1) + cstop(i)
      end do
c
c     set column index by scanning rows in reverse order
c
      do i = nvar, 1, -1
         do j = hinit(i), hstop(i)
            k = hindex(j)
            m = cinit(k) - 1
            cinit(k) = m
            cindex(m) = i
            cvalue(m) = j
         end do
      end do
c
c     convert from number of elements to end marker for column
c
      do i = 1, nvar
         cstop(i) = cinit(i) + cstop(i) - 1
      end do
      return
      end
