c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine cluster  --  set user-defined groups of atoms  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "cluster" gets the partitioning of the system into groups
c     and stores a list of the group to which each atom belongs
c
c
      subroutine cluster
      use sizes
      use atomid
      use atoms
      use bound
      use group
      use inform
      use iounit
      use keys
      use limits
      use molcul
      implicit none
      integer i,j,k
      integer next,size
      integer gnum,ga,gb
      integer, allocatable :: list(:)
      real*8 wg
      logical header
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     perform dynamic allocation of some global arrays
c
      if (.not. allocated(igrp))  allocate (igrp(2,0:maxgrp))
      if (.not. allocated(grpmass))  allocate (grpmass(0:maxgrp))
      if (.not. allocated(wgrp))  allocate (wgrp(0:maxgrp,0:maxgrp))
      if (allocated(kgrp))  deallocate (kgrp)
      if (allocated(grplist))  deallocate (grplist)
      allocate (kgrp(n))
      allocate (grplist(n))
c
c     set defaults for the group atom list and weight options
c
      use_group = .false.
      use_intra = .false.
      use_inter = .false.
      ngrp = 0
      do i = 1, n
         kgrp(i) = 0
         grplist(i) = 0
      end do
      do i = 1, maxgrp
         igrp(1,i) = 1
         igrp(2,i) = 0
      end do
      do i = 0, maxgrp
         do j = 0, maxgrp
            wgrp(j,i) = 1.0d0
         end do
      end do
c
c     perform dynamic allocation of some local arrays
c
      size = max(20,n)
      allocate (list(size))
c
c     get any keywords containing atom group definitions
c
      do j = 1, nkey
         next = 1
         record = keyline(j)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:6) .eq. 'GROUP ') then
            use_group = .true.
            gnum = 0
            do i = 1, 20
               list(i) = 0
            end do
            call getnumb (record,gnum,next)
            if (gnum .gt. maxgrp) then
               write (iout,10)
   10          format (/,' CLUSTER  --  Too many Atom Groups;',
     &                    ' Increase MAXGRP')
               call fatal
            end if
            string = record(next:240)
            read (string,*,err=20,end=20)  (list(i),i=1,20)
   20       continue
            i = 1
            do while (list(i) .ne. 0)
               if (list(i) .gt. 0) then
                  grplist(list(i)) = gnum
                  i = i + 1
               else
                  do k = abs(list(i)), abs(list(i+1))
                     grplist(k) = gnum
                  end do
                  i = i + 2
               end if
            end do
c
c     get any keywords with weights for group interactions
c
         else if (keyword(1:15) .eq. 'GROUP-MOLECULE ') then
            use_group = .true.
            use_inter = .true.
            use_intra = .false.
            if (nmol .gt. maxgrp) then
               write (iout,30)
   30          format (/,' CLUSTER  --  Too many Atom Groups;',
     &                    ' Increase MAXGRP')
               call fatal
            end if
            do i = 1, nmol
               do k = imol(1,i), imol(2,i)
                  grplist(kmol(k)) = i
               end do
            end do
c
c     get any keywords with weights for group interactions
c
         else if (keyword(1:13) .eq. 'GROUP-SELECT ') then
            ga = 0
            gb = 0
            wg = -1.0d0
            string = record(next:240)
            read (string,*,err=40,end=40)  ga,gb,wg
   40       continue
            if (wg .lt. 0.0d0)  wg = 1.0d0
            wgrp(ga,gb) = wg
            wgrp(gb,ga) = wg
            use_inter = .false.
c
c     get keywords to select common sets of group interactions
c
         else if (keyword(1:12) .eq. 'GROUP-INTRA ') then
            use_intra = .true.
            use_inter = .false.
         else if (keyword(1:12) .eq. 'GROUP-INTER ') then
            use_inter = .true.
            use_intra = .false.
         end if
      end do
c
c     pack atoms of each group into a contiguous indexed list
c
      if (use_group) then
         do i = 1, n
            list(i) = grplist(i)
         end do
         call sort3 (n,list,kgrp)
c
c     find the first and last atom in each of the groups
c
         k = list(1)
         igrp(1,k) = 1
         do i = 1, n
            j = list(i)
            if (j .ne. k) then
               igrp(2,k) = i - 1
               igrp(1,j) = i
               k = j
            end if
            ngrp = max(j,ngrp)
         end do
         igrp(2,j) = n
c
c     sort the list of atoms in each group by atom number
c
         do i = 0, ngrp
            size = igrp(2,i) - igrp(1,i) + 1
            if (igrp(1,i) .ne. 0)
     &         call sort (size,kgrp(igrp(1,i)))
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (list)
c
c     use only intragroup or intergroup interactions if selected
c
      if (use_intra) then
         do i = 0, ngrp
            do j = 0, ngrp
               wgrp(j,i) = 0.0d0
            end do
            wgrp(i,i) = 1.0d0
         end do
      end if
      if (use_inter) then
         do i = 0, ngrp
            do j = 0, ngrp
               wgrp(j,i) = 1.0d0
            end do
            wgrp(i,i) = 0.0d0
         end do
      end if
c
c     disable consideration of interactions with any empty groups
c
      do i = 0, ngrp
         size = igrp(2,i) - igrp(1,i) + 1
         if (size .eq. 0) then
            do j = 0, ngrp
               wgrp(j,i) = 0.0d0
               wgrp(i,j) = 0.0d0
            end do
         end if
      end do
c
c     turn off bounds and replicas for intragroup calculations
c
      if (use_intra) then
         use_bounds = .false.
         use_replica = .false.
         call cutoffs
      end if
c
c     compute the total mass of all atoms in each group
c
      do i = 1, ngrp
         grpmass(i) = 0.0d0
         do j = igrp(1,i), igrp(2,i)
            grpmass(i) = grpmass(i) + mass(kgrp(j))
         end do
      end do
c
c     output the final list of atoms in each group
c
      if (debug .and. use_group) then
         do i = 1, ngrp
            size = igrp(2,i) - igrp(1,i) + 1
            if (size .ne. 0) then
               write (iout,50)  i
   50          format (/,' List of Atoms in Group',i3,' :',/)
               write (iout,60)  (kgrp(j),j=igrp(1,i),igrp(2,i))
   60          format (3x,10i7)
            end if
         end do
      end if
c
c     output the weights for intragroup and intergroup interactions
c
      if (debug .and. use_group) then
         header = .true.
         do i = 0, ngrp
            do j = i, ngrp
               if (wgrp(j,i) .ne. 0.0d0) then
                  if (header) then
                     header = .false.
                     write (iout,70)
   70                format (/,' Active Sets of Intra- and InterGroup',
     &                          ' Interactions :',
     &                       //,11x,'Groups',15x,'Type',14x,'Weight',/)
                  end if
                  if (i .eq. j) then
                     write (iout,80)  i,j,wgrp(j,i)
   80                format (5x,2i6,12x,'IntraGroup',5x,f12.4)
                  else
                     write (iout,90)  i,j,wgrp(j,i)
   90                format (5x,2i6,12x,'InterGroup',5x,f12.4)
                  end if
               end if
            end do
         end do
      end if
      return
      end
