c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine kopdist  --  out-of-plane distance parameters  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "kopdist" assigns the force constants for out-of-plane
c     distance at trigonal centers via the central atom height;
c     also processes any new or changed parameter values
c
c
      subroutine kopdist
      use sizes
      use angbnd
      use angpot
      use atmlst
      use atomid
      use atoms
      use couple
      use inform
      use iounit
      use keys
      use kopdst
      use opdist
      use potent
      implicit none
      integer i,j,k,nopd
      integer ia,ib,ic,id
      integer ita,itb,itc,itd
      integer imin,itmin
      integer size,next
      real*8 fopd
      logical header
      character*4 pa,pb,pc,pd
      character*12 zeros
      character*16 blank
      character*16 pt,pt0
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     process keywords containing out-of-plane distance parameters
c
      blank = '                '
      zeros = '000000000000'
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'OPDIST ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            fopd = 0.0d0
            string = record(next:240)
            read (string,*,err=10,end=10)  ia,ib,ic,id,fopd
   10       continue
            size = 4
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            imin = min(ib,ic,id)
            if (ib .eq. imin) then
               if (ic .le. id) then
                  pt = pa//pb//pc//pd
               else
                  pt = pa//pb//pd//pc
               end if
            else if (ic .eq. imin) then
               if (ib .le. id) then
                  pt = pa//pc//pb//pd
               else
                  pt = pa//pc//pd//pb
               end if
            else if (id .eq. imin) then
               if (ib .le. ic) then
                  pt = pa//pd//pb//pc
               else
                  pt = pa//pd//pc//pb
               end if
            end if
            if (.not. silent) then
               if (header) then
                  header = .false.
                  write (iout,20)
   20             format (/,' Additional Out-of-Plane Distance',
     &                       ' Parameters :',
     &                    //,5x,'Atom Classes',19x,'K(OPD)',/)
               end if
               write (iout,30)  ia,ib,ic,id,fopd
   30          format (4x,4i4,10x,2f12.3)
            end if
            do j = 1, maxnopd
               if (kopd(j).eq.blank .or. kopd(j).eq.pt) then
                  kopd(j) = pt
                  opds(j) = fopd
                  goto 50
               end if
            end do
            write (iout,40)
   40       format (/,' KOPDIST  --  Too many Out-of-Plane Distance',
     &                 ' Parameters')
            abort = .true.
   50       continue
         end if
      end do
c
c     determine the total number of forcefield parameters
c
      nopd = maxnopd
      do i = maxnopd, 1, -1
         if (kopd(i) .eq. blank)  nopd = i - 1
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (.not. allocated(angtyp))  allocate (angtyp(nangle))
      if (allocated(iopd))  deallocate (iopd)
      if (allocated(opdk))  deallocate (opdk)
      allocate (iopd(4,n))
      allocate (opdk(n))
c
c     assign out-of-plane distance parameters for trigonal sites
c
      nopdist = 0
      if (nopd .ne. 0) then
         do i = 1, n
            if (n12(i) .eq. 3) then
               ia = i
               ib = i12(1,i)
               ic = i12(2,i)
               id = i12(3,i)
               ita = class(ia)
               itb = class(ib)
               itc = class(ic)
               itd = class(id)
               size = 4
               call numeral (ita,pa,size)
               call numeral (itb,pb,size)
               call numeral (itc,pc,size)
               call numeral (itd,pd,size)
               itmin = min(itb,itc,itd)
               if (itb .eq. itmin) then
                  if (itc .le. itd) then
                     pt = pa//pb//pc//pd
                  else
                     pt = pa//pb//pd//pc
                  end if
               else if (itc .eq. itmin) then
                  if (itb .le. itd) then
                     pt = pa//pc//pb//pd
                  else
                     pt = pa//pc//pd//pb
                  end if
               else if (itd .eq. itmin) then
                  if (itb .le. itc) then
                     pt = pa//pd//pb//pc
                  else
                     pt = pa//pd//pc//pb
                  end if
               end if
               pt0 = pa//zeros
               do j = 1, nopd
                  if (kopd(j) .eq. pt) then
                     nopdist = nopdist + 1
                     iopd(1,nopdist) = ia
                     iopd(2,nopdist) = ib
                     iopd(3,nopdist) = ic
                     iopd(4,nopdist) = id
                     opdk(nopdist) = opds(j)
                     goto 60
                  end if
               end do
               do j = 1, nopd
                  if (kopd(j) .eq. pt0) then
                     nopdist = nopdist + 1
                     iopd(1,nopdist) = ia
                     iopd(2,nopdist) = ib
                     iopd(3,nopdist) = ic
                     iopd(4,nopdist) = id
                     opdk(nopdist) = opds(j)
                     goto 60
                  end if
               end do
   60          continue
            end if
         end do
      end if
c
c     mark angles at trigonal sites to use projected in-plane values
c
      do i = 1, nopdist
         ia = iopd(1,i)
         do j = 1, 3
            k = anglist(j,ia)
            if (angtyp(k) .eq. 'HARMONIC')  angtyp(k) = 'IN-PLANE'
         end do
      end do
c
c     turn off out-of-plane distance potential if it is not used
c
      if (nopdist .eq. 0)  use_opdist = .false.
      return
      end
