c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine kangang  --  angle-angle parameter assignment  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "kangang" assigns the parameters for angle-angle cross term
c     interactions and processes new or changed parameter values
c
c
      subroutine kangang
      use sizes
      use angang
      use angbnd
      use atmlst
      use atomid
      use atoms
      use couple
      use inform
      use iounit
      use kanang
      use keys
      use potent
      use tors
      implicit none
      integer i,j,k,m,next
      integer it,ia,ic
      integer nang,jang,kang
      integer maxaa
      real*8 fa,faa,aak(3)
      logical header
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     process keywords containing angle-angle parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'ANGANG ') then
            it = 0
            do j = 1, 3
               aak(j) = 0.0d0
            end do
            string = record(next:240)
            read (string,*,err=10,end=10)  it,(aak(j),j=1,3)
   10       continue
            if (.not. silent) then
               if (header) then
                  header = .false.
                  write (iout,20)
   20             format (/,' Additional Angle-Angle Parameters :',
     &                    //,5x,'Atom Class',8x,'K(AA) 1',4x,'K(AA) 2',
     &                       4x,'K(AA) 3',/)
               end if
               write (iout,30)  it,(aak(j),j=1,3)
   30          format (9x,i3,7x,3f11.3)
            end if
            do j = 1, 3
               anan(j,it) = aak(j)
            end do
         end if
      end do
c
c     perform dynamic allocation of some global arrays
c
      maxaa = 15 * n
      if (allocated(iaa))  deallocate (iaa)
      if (allocated(kaa))  deallocate (kaa)
      allocate (iaa(2,maxaa))
      allocate (kaa(maxaa))
c
c     assign the angle-angle parameters for each angle pair
c
      nangang = 0
      do i = 1, n
         nang = n12(i) * (n12(i)-1) / 2
         it = class(i)
         do j = 1, nang-1
            jang = anglist(j,i)
            ia = iang(1,jang)
            ic = iang(3,jang)
            m = 1
            if (atomic(ia) .le. 1)  m = m + 1
            if (atomic(ic) .le. 1)  m = m + 1
            fa = anan(m,it)
            do k = j+1, nang
               kang = anglist(k,i)
               ia = iang(1,kang)
               ic = iang(3,kang)
               m = 1
               if (atomic(ia) .le. 1)  m = m + 1
               if (atomic(ic) .le. 1)  m = m + 1
               faa = fa * anan(m,it)
               if (faa .ne. 0.0d0) then
                  nangang = nangang + 1
                  iaa(1,nangang) = jang
                  iaa(2,nangang) = kang
                  kaa(nangang) = faa
               end if
            end do
         end do
      end do
c
c     turn off the angle-angle potential if it is not used
c
      if (nangang .eq. 0)  use_angang = .false.
      return
      end
