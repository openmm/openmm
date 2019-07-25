c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine sort  --  heapsort of an integer array  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "sort" takes an input list of integers and sorts it
c     into ascending order using the Heapsort algorithm
c
c
      subroutine sort (n,list)
      implicit none
      integer i,j,k,n
      integer index,lists
      integer list(*)
c
c
c     perform the heapsort of the input list
c
      k = n/2 + 1
      index = n
      do while (n .gt. 1)
         if (k .gt. 1) then
            k = k - 1
            lists = list(k)
         else
            lists = list(index)
            list(index) = list(1)
            index = index - 1
            if (index .le. 1) then
               list(1) = lists
               return
            end if
         end if
         i = k
         j = k + k
         do while (j .le. index)
            if (j .lt. index) then
               if (list(j) .lt. list(j+1))  j = j + 1
            end if
            if (lists .lt. list(j)) then
               list(i) = list(j)
               i = j
               j = j + j
            else
               j = index + 1
            end if
         end do
         list(i) = lists
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine sort2  --  heapsort of real array with keys  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "sort2" takes an input list of reals and sorts it
c     into ascending order using the Heapsort algorithm;
c     it also returns a key into the original ordering
c
c
      subroutine sort2 (n,list,key)
      implicit none
      integer i,j,k,n
      integer index,keys
      integer key(*)
      real*8 lists
      real*8 list(*)
c
c
c     initialize index into the original ordering
c
      do i = 1, n
         key(i) = i
      end do
c
c     perform the heapsort of the input list
c
      k = n/2 + 1
      index = n
      do while (n .gt. 1)
         if (k .gt. 1) then
            k = k - 1
            lists = list(k)
            keys = key(k)
         else
            lists = list(index)
            keys = key(index)
            list(index) = list(1)
            key(index) = key(1)
            index = index - 1
            if (index .le. 1) then
               list(1) = lists
               key(1) = keys
               return
            end if
         end if
         i = k
         j = k + k
         do while (j .le. index)
            if (j .lt. index) then
               if (list(j) .lt. list(j+1))  j = j + 1
            end if
            if (lists .lt. list(j)) then
               list(i) = list(j)
               key(i) = key(j)
               i = j
               j = j + j
            else
               j = index + 1
            end if
         end do
         list(i) = lists
         key(i) = keys
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine sort3  --  heapsort of integer array with keys  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "sort3" takes an input list of integers and sorts it
c     into ascending order using the Heapsort algorithm;
c     it also returns a key into the original ordering
c
c
      subroutine sort3 (n,list,key)
      implicit none
      integer i,j,k,n
      integer index
      integer lists
      integer keys
      integer list(*)
      integer key(*)
c
c
c     initialize index into the original ordering
c
      do i = 1, n
         key(i) = i
      end do
c
c     perform the heapsort of the input list
c
      k = n/2 + 1
      index = n
      do while (n .gt. 1)
         if (k .gt. 1) then
            k = k - 1
            lists = list(k)
            keys = key(k)
         else
            lists = list(index)
            keys = key(index)
            list(index) = list(1)
            key(index) = key(1)
            index = index - 1
            if (index .le. 1) then
               list(1) = lists
               key(1) = keys
               return
            end if
         end if
         i = k
         j = k + k
         do while (j .le. index)
            if (j .lt. index) then
               if (list(j) .lt. list(j+1))  j = j + 1
            end if
            if (lists .lt. list(j)) then
               list(i) = list(j)
               key(i) = key(j)
               i = j
               j = j + j
            else
               j = index + 1
            end if
         end do
         list(i) = lists
         key(i) = keys
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine sort4  --  heapsort of integer absolute values  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "sort4" takes an input list of integers and sorts it into
c     ascending absolute value using the Heapsort algorithm
c
c
      subroutine sort4 (n,list)
      implicit none
      integer i,j,k,n
      integer index
      integer lists
      integer list(*)
c
c
c     perform the heapsort of the input list
c
      k = n/2 + 1
      index = n
      do while (n .gt. 1)
         if (k .gt. 1) then
            k = k - 1
            lists = list(k)
         else
            lists = list(index)
            list(index) = list(1)
            index = index - 1
            if (index .le. 1) then
               list(1) = lists
               return
            end if
         end if
         i = k
         j = k + k
         do while (j .le. index)
            if (j .lt. index) then
               if (abs(list(j)) .lt. abs(list(j+1)))  j = j + 1
            end if
            if (abs(lists) .lt. abs(list(j))) then
               list(i) = list(j)
               i = j
               j = j + j
            else
               j = index + 1
            end if
         end do
         list(i) = lists
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine sort5  --  heapsort of integer array modulo m  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "sort5" takes an input list of integers and sorts it
c     into ascending order based on each value modulo "m"
c
c
      subroutine sort5 (n,list,m)
      implicit none
      integer i,j,k,m,n
      integer index,smod
      integer jmod,j1mod
      integer lists
      integer list(*)
c
c
c     perform the heapsort of the input list
c
      k = n/2 + 1
      index = n
      do while (n .gt. 1)
         if (k .gt. 1) then
            k = k - 1
            lists = list(k)
         else
            lists = list(index)
            list(index) = list(1)
            index = index - 1
            if (index .le. 1) then
               list(1) = lists
               return
            end if
         end if
         i = k
         j = k + k
         do while (j .le. index)
            if (j .lt. index) then
               jmod = mod(list(j),m)
               j1mod = mod(list(j+1),m)
               if (jmod .lt. j1mod) then
                  j = j + 1
               else if (jmod.eq.j1mod .and. list(j).lt.list(j+1)) then
                  j = j + 1
               end if
            end if
            smod = mod(lists,m)
            jmod = mod(list(j),m)
            if (smod .lt. jmod) then
               list(i) = list(j)
               i = j
               j = j + j
            else if (smod.eq.jmod .and. lists.lt.list(j)) then
               list(i) = list(j)
               i = j
               j = j + j
            else
               j = index + 1
            end if
         end do
         list(i) = lists
      end do
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine sort6  --  heapsort of a text string array  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "sort6" takes an input list of character strings and sorts
c     it into alphabetical order using the Heapsort algorithm
c
c
      subroutine sort6 (n,list)
      implicit none
      integer i,j,k,n
      integer index
      character*256 lists
      character*(*) list(*)
c
c
c     perform the heapsort of the input list
c
      k = n/2 + 1
      index = n
      do while (n .gt. 1)
         if (k .gt. 1) then
            k = k - 1
            lists = list(k)
         else
            lists = list(index)
            list(index) = list(1)
            index = index - 1
            if (index .le. 1) then
               list(1) = lists
               return
            end if
         end if
         i = k
         j = k + k
         do while (j .le. index)
            if (j .lt. index) then
               if (list(j) .lt. list(j+1))  j = j + 1
            end if
            if (lists .lt. list(j)) then
               list(i) = list(j)
               i = j
               j = j + j
            else
               j = index + 1
            end if
         end do
         list(i) = lists
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine sort7  --  heapsort of text strings with keys  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "sort7" takes an input list of character strings and sorts it
c     into alphabetical order using the Heapsort algorithm; it also
c     returns a key into the original ordering
c
c
      subroutine sort7 (n,list,key)
      implicit none
      integer i,j,k,n
      integer index
      integer keys
      integer key(*)
      character*256 lists
      character*(*) list(*)
c
c
c     initialize index into the original ordering
c
      do i = 1, n
         key(i) = i
      end do
c
c     perform the heapsort of the input list
c
      k = n/2 + 1
      index = n
      do while (n .gt. 1)
         if (k .gt. 1) then
            k = k - 1
            lists = list(k)
            keys = key(k)
         else
            lists = list(index)
            keys = key(index)
            list(index) = list(1)
            key(index) = key(1)
            index = index - 1
            if (index .le. 1) then
               list(1) = lists
               key(1) = keys
               return
            end if
         end if
         i = k
         j = k + k
         do while (j .le. index)
            if (j .lt. index) then
               if (list(j) .lt. list(j+1))  j = j + 1
            end if
            if (lists .lt. list(j)) then
               list(i) = list(j)
               key(i) = key(j)
               i = j
               j = j + j
            else
               j = index + 1
            end if
         end do
         list(i) = lists
         key(i) = keys
      end do
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine sort8  --  heapsort to unique integers  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "sort8" takes an input list of integers and sorts it into
c     ascending order using the Heapsort algorithm, duplicate
c     values are removed from the final sorted list
c
c
      subroutine sort8 (n,list)
      implicit none
      integer i,j,k,n
      integer index
      integer lists
      integer list(*)
c
c
c     perform the heapsort of the input list
c
      k = n/2 + 1
      index = n
      do while (n .gt. 1)
         if (k .gt. 1) then
            k = k - 1
            lists = list(k)
         else
            lists = list(index)
            list(index) = list(1)
            index = index - 1
            if (index .le. 1) then
               list(1) = lists
c
c     remove duplicate values from final list
c
               j = 1
               do i = 2, n
                  if (list(i-1) .ne. list(i)) then
                     j = j + 1
                     list(j) = list(i)
                  end if
               end do
               if (j .lt. n)  n = j
               return
            end if
         end if
         i = k
         j = k + k
         do while (j .le. index)
            if (j .lt. index) then
               if (list(j) .lt. list(j+1))  j = j + 1
            end if
            if (lists .lt. list(j)) then
               list(i) = list(j)
               i = j
               j = j + j
            else
               j = index + 1
            end if
         end do
         list(i) = lists
      end do
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine sort9  --  heapsort to unique real values  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "sort9" takes an input list of reals and sorts it into
c     ascending order using the Heapsort algorithm, duplicate
c     values are removed from the final sorted list
c
c
      subroutine sort9 (n,list)
      implicit none
      integer i,j,k,n
      integer index
      real*8 lists
      real*8 list(*)
c
c
c     perform the heapsort of the input list
c
      k = n/2 + 1
      index = n
      do while (n .gt. 1)
         if (k .gt. 1) then
            k = k - 1
            lists = list(k)
         else
            lists = list(index)
            list(index) = list(1)
            index = index - 1
            if (index .le. 1) then
               list(1) = lists
c
c     remove duplicate values from final list
c
               j = 1
               do i = 2, n
                  if (list(i-1) .ne. list(i)) then
                     j = j + 1
                     list(j) = list(i)
                  end if
               end do
               if (j .lt. n)  n = j
               return
            end if
         end if
         i = k
         j = k + k
         do while (j .le. index)
            if (j .lt. index) then
               if (list(j) .lt. list(j+1))  j = j + 1
            end if
            if (lists .lt. list(j)) then
               list(i) = list(j)
               i = j
               j = j + j
            else
               j = index + 1
            end if
         end do
         list(i) = lists
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine sort10  --  heapsort to unique text strings  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "sort10" takes an input list of character strings and sorts
c     it into alphabetical order using the Heapsort algorithm,
c     duplicate values are removed from the final sorted list
c
c
      subroutine sort10 (n,list)
      implicit none
      integer i,j,k,n
      integer index
      character*256 lists
      character*(*) list(*)
c
c
c     perform the heapsort of the input list
c
      k = n/2 + 1
      index = n
      do while (n .gt. 1)
         if (k .gt. 1) then
            k = k - 1
            lists = list(k)
         else
            lists = list(index)
            list(index) = list(1)
            index = index - 1
            if (index .le. 1) then
               list(1) = lists
c
c     remove duplicate values from final list
c
               j = 1
               do i = 2, n
                  if (list(i-1) .ne. list(i)) then
                     j = j + 1
                     list(j) = list(i)
                  end if
               end do
               if (j .lt. n)  n = j
               return
            end if
         end if
         i = k
         j = k + k
         do while (j .le. index)
            if (j .lt. index) then
               if (list(j) .lt. list(j+1))  j = j + 1
            end if
            if (lists .lt. list(j)) then
               list(i) = list(j)
               i = j
               j = j + j
            else
               j = index + 1
            end if
         end do
         list(i) = lists
      end do
      return
      end
