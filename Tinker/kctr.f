c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kctr  --  Charge transfer parameter assignment##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kctr" assigns the parameters to be used in computing the
c     van der Waals interactions and processes any new or changed
c     values for these parameters
c
c
      subroutine kctr
      use sizes
      use atomid
      use atoms
      use couple
      use fields
      use inform
      use iounit
      use keys
      use math
      use potent
      use ctran
      implicit none
      integer i,k
      integer next
      real*8 apr,bex 
      real*8 pen 
      logical header
      character*8 blank
      character*20 keyword
      character*240 record
      character*240 string
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(ict))  deallocate (ict)
      if (allocated(jct))  deallocate (jct)
      allocate (ict(n))
      allocate (jct(n))

      do i = 1, n
         jct(i) = type(i)
      end do
c
c     process keywords containing charge transfer parameters
c
      blank = '        '
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:3) .eq. 'CT ') then
            call getnumb (record,k,next)
            if (k.ge.1 .and. k.le.maxtyp) then
               apr = apre(k)
               bex = bexp(k)
               string = record(next:120)
               read (string,*,err=10,end=10) apr,bex 
  10           continue
               apre(k) = apr 
               bexp(k) = bex
            end if
         end if
         if ((keyword(1:9)) .eq. 'APRERULE ') then
            call getnumb (record,k,next)
            string = record(next:120) 
            read (string,*,err=20,end=20)aprerule 
  20        continue
         end if 
         if ((keyword(1:9)) .eq. 'BEXPRULE ') then
            call getnumb (record,k,next)
            string = record(next:120) 
            read (string,*,err=30,end=30)bexprule 
  30        continue
         end if 
      end do

c
c     remove zero-sized atoms from the list of CT sites
c
      nct = 0 
      do i = 1, n
         if (apre(jct(i)) .ne. 0.0d0) then
            nct = nct + 1
            ict(nct) = i
         end if
      end do
c
c     turn off the CT potential if it is not used
c
      if (nct .eq. 0)  use_ct = .false.
      return
      end
