c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine kpitors  --  find pi-system torsion parameters  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "kpitors" assigns pi-system torsion parameters to torsions
c     needing them, and processes any new or changed values
c
c
      subroutine kpitors
      use sizes
      use atomid
      use atoms
      use bndstr
      use couple
      use inform
      use iounit
      use keys
      use kpitor
      use pitors
      use potent
      use tors
      implicit none
      integer i,j,npt
      integer ia,ib
      integer ita,itb
      integer size,next
      real*8 tp
      logical header
      character*4 pa,pb
      character*8 blank,pt
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     process keywords containing pi-system torsion parameters
c
      blank = '        '
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'PITORS ') then
            ia = 0
            ib = 0
            tp = 0.0d0
            string = record(next:240)
            read (string,*,err=10,end=10)  ia,ib,tp
   10       continue
            if (.not. silent) then
               if (header) then
                  header = .false.
                  write (iout,20)
   20             format (/,' Additional Pi-Orbital Torsion',
     &                       ' Parameters :',
     &                    //,5x,'Atom Classes',7x,'2-Fold',/)
               end if
               write (iout,30)  ia,ib,tp
   30          format (6x,2i4,4x,f12.3)
            end if
            size = 4
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            if (ia .le. ib) then
               pt = pa//pb
            else
               pt = pb//pa
            end if
            do j = 1, maxnpt
               if (kpt(j).eq.blank .or. kpt(j).eq.pt) then
                  kpt(j) = pt
                  ptcon(j) = tp
                  goto 50
               end if
            end do
            write (iout,40)
   40       format (/,' KPITORS  --  Too many Pi-Orbital Torsion',
     &                 ' Parameters')
            abort = .true.
   50       continue
         end if
      end do
c
c     determine the total number of forcefield parameters
c
      npt = maxnpt
      do i = maxnpt, 1, -1
         if (kpt(i) .eq. blank)  npt = i - 1
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(ipit))  deallocate (ipit)
      if (allocated(kpit))  deallocate (kpit)
      allocate (ipit(6,ntors))
      allocate (kpit(ntors))
c
c     assign pi-system torsion parameters as required
c
      npitors = 0
      if (npt .ne. 0) then
         do i = 1, nbond
            ia = ibnd(1,i)
            ib = ibnd(2,i)
            if (n12(ia).eq.3 .and. n12(ib).eq.3) then
               ita = class(ia)
               itb = class(ib)
               size = 4
               call numeral (ita,pa,size)
               call numeral (itb,pb,size)
               if (ita .le. itb) then
                  pt = pa//pb
               else
                  pt = pb//pa
               end if
               do j = 1, npt
                  if (kpt(j) .eq. pt) then
                     npitors = npitors + 1
                     kpit(npitors) = ptcon(j)
                     ipit(1,npitors) = i12(1,ia)
                     ipit(2,npitors) = i12(2,ia)
                     ipit(3,npitors) = ia
                     ipit(4,npitors) = ib
                     ipit(5,npitors) = i12(1,ib)
                     ipit(6,npitors) = i12(2,ib)
                     if (i12(1,ia) .eq. ib)
     &                  ipit(1,npitors) = i12(3,ia)
                     if (i12(2,ia) .eq. ib)
     &                  ipit(2,npitors) = i12(3,ia)
                     if (i12(1,ib) .eq. ia)
     &                  ipit(5,npitors) = i12(3,ib)
                     if (i12(2,ib) .eq. ia)
     &                  ipit(6,npitors) = i12(3,ib)
                     goto 60
                  end if
               end do
            end if
   60       continue
         end do
      end if
c
c     turn off the pi-system torsion potential if it is not used
c
      if (npitors .eq. 0)  use_pitors = .false.
      return
      end
