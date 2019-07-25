c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1991  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kstrtor  --  find stretch-torsion parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kstrtor" assigns stretch-torsion parameters to torsions
c     needing them, and processes any new or changed values
c
c
      subroutine kstrtor
      use sizes
      use atmlst
      use atomid
      use atoms
      use couple
      use inform
      use iounit
      use keys
      use ksttor
      use potent
      use strtor
      use tors
      implicit none
      integer i,j,k,nbt
      integer ia,ib,ic,id
      integer ita,itb,itc,itd
      integer size,next
      real*8 bt1,bt2,bt3
      real*8 bt4,bt5,bt6
      real*8 bt7,bt8,bt9
      logical header,swap
      character*4 pa,pb,pc,pd
      character*4 zeros
      character*16 blank
      character*16 pt,pt0
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     process keywords containing stretch-torsion parameters
c
      blank = '                '
      zeros = '0000'
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:8) .eq. 'STRTORS ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            bt1 = 0.0d0
            bt2 = 0.0d0
            bt3 = 0.0d0
            bt4 = 0.0d0
            bt5 = 0.0d0
            bt6 = 0.0d0
            bt7 = 0.0d0
            bt8 = 0.0d0
            bt9 = 0.0d0
            string = record(next:240)
            read (string,*,err=10,end=10)  ia,ib,ic,id,bt1,bt2,bt3,
     &                                     bt4,bt5,bt6,bt7,bt8,bt9
   10       continue
            size = 4
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            if (ib .lt. ic) then
               pt = pa//pb//pc//pd
               swap = .false.
            else if (ic .lt. ib) then
               pt = pd//pc//pb//pa
               swap = .true.
            else if (ia .le. id) then
               pt = pa//pb//pc//pd
               swap = .false.
            else if (id .lt. ia) then
               pt = pd//pc//pb//pa
               swap = .true.
            end if
            if (.not. silent) then
               if (header) then
                  header = .false.
                  write (iout,20)
   20             format (/,' Additional Stretch-Torsion Parameters :',
     &                    //,5x,'Atom Classes',9x,'Stretch',7x,'1-Fold',
     &                       6x,'2-Fold',6x,'3-Fold',/)
               end if
               write (iout,30)  ia,ib,ic,id,bt1,bt2,bt3,
     &                          bt4,bt5,bt6,bt7,bt8,bt9
   30          format (1x,4i4,8x,'1st Bond',1x,3f12.3,
     &                 /,25x,'2nd Bond',1x,3f12.3,
     &                 /,25x,'3rd Bond',1x,3f12.3)
            end if
            do j = 1, maxnbt
               if (kbt(j).eq.blank .or. kbt(j).eq.pt) then
                  kbt(j) = pt
                  btcon(4,j) = bt4
                  btcon(5,j) = bt5
                  btcon(6,j) = bt6
                  if (swap) then
                     btcon(1,j) = bt7
                     btcon(2,j) = bt8
                     btcon(3,j) = bt9
                     btcon(7,j) = bt1
                     btcon(8,j) = bt2
                     btcon(9,j) = bt3
                  else
                     btcon(1,j) = bt1
                     btcon(2,j) = bt2
                     btcon(3,j) = bt3
                     btcon(7,j) = bt7
                     btcon(8,j) = bt8
                     btcon(9,j) = bt9
                  end if
                  goto 50
               end if
            end do
            write (iout,40)
   40       format (/,' KSTRTOR  --  Too many Stretch-Torsion',
     &                 ' Parameters')
            abort = .true.
   50       continue
         end if
      end do
c
c     determine the total number of forcefield parameters
c
      nbt = maxnbt
      do i = maxnbt, 1, -1
         if (kbt(i) .eq. blank)  nbt = i - 1
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(ist))  deallocate (ist)
      if (allocated(kst))  deallocate (kst)
      allocate (ist(4,ntors))
      allocate (kst(9,ntors))
c
c     assign the stretch-torsion parameters for each torsion
c
      nstrtor = 0
      if (nbt .ne. 0) then
         do i = 1, ntors
            ia = itors(1,i)
            ib = itors(2,i)
            ic = itors(3,i)
            id = itors(4,i)
            ita = class(ia)
            itb = class(ib)
            itc = class(ic)
            itd = class(id)
            size = 4
            call numeral (ita,pa,size)
            call numeral (itb,pb,size)
            call numeral (itc,pc,size)
            call numeral (itd,pd,size)
            if (itb .lt. itc) then
               pt = pa//pb//pc//pd
               swap = .false.
            else if (itc .lt. itb) then
               pt = pd//pc//pb//pa
               swap = .true.
            else if (ita .le. itd) then
               pt = pa//pb//pc//pd
               swap = .false.
            else if (itd .lt. ita) then
               pt = pd//pc//pb//pa
               swap = .true.
            end if
            pt0 = zeros//pt(5:12)//zeros
            do j = 1, nbt
               if (kbt(j) .eq. pt) then
                  nstrtor = nstrtor + 1
                  kst(4,nstrtor) = btcon(4,j)
                  kst(5,nstrtor) = btcon(5,j)
                  kst(6,nstrtor) = btcon(6,j)
                  if (swap) then
                     kst(1,nstrtor) = btcon(7,j)
                     kst(2,nstrtor) = btcon(8,j)
                     kst(3,nstrtor) = btcon(9,j)
                     kst(7,nstrtor) = btcon(1,j)
                     kst(8,nstrtor) = btcon(2,j)
                     kst(9,nstrtor) = btcon(3,j)
                  else
                     kst(1,nstrtor) = btcon(1,j)
                     kst(2,nstrtor) = btcon(2,j)
                     kst(3,nstrtor) = btcon(3,j)
                     kst(7,nstrtor) = btcon(7,j)
                     kst(8,nstrtor) = btcon(8,j)
                     kst(9,nstrtor) = btcon(9,j)
                  end if
                  ist(1,nstrtor) = i
                  do k = 1, n12(ia)
                     if (i12(k,ia) .eq. ib) then
                        ist(2,nstrtor) = bndlist(k,ia)
                        goto 60
                     endif
                  end do
   60             continue
                  do k = 1, n12(ib)
                     if (i12(k,ib) .eq. ic) then
                        ist(3,nstrtor) = bndlist(k,ib)
                        goto 70
                     end if
                  end do
   70             continue
                  do k = 1, n12(ic)
                     if (i12(k,ic) .eq. id) then
                        ist(4,nstrtor) = bndlist(k,ic)
                        goto 100
                     end if
                  end do
               end if
            end do
            do j = 1, nbt
               if (kbt(j) .eq. pt0) then
                  nstrtor = nstrtor + 1
                  kst(4,nstrtor) = btcon(4,j)
                  kst(5,nstrtor) = btcon(5,j)
                  kst(6,nstrtor) = btcon(6,j)
                  if (swap) then
                     kst(1,nstrtor) = btcon(7,j)
                     kst(2,nstrtor) = btcon(8,j)
                     kst(3,nstrtor) = btcon(9,j)
                     kst(7,nstrtor) = btcon(1,j)
                     kst(8,nstrtor) = btcon(2,j)
                     kst(9,nstrtor) = btcon(3,j)
                  else
                     kst(1,nstrtor) = btcon(1,j)
                     kst(2,nstrtor) = btcon(2,j)
                     kst(3,nstrtor) = btcon(3,j)
                     kst(7,nstrtor) = btcon(7,j)
                     kst(8,nstrtor) = btcon(8,j)
                     kst(9,nstrtor) = btcon(9,j)
                  end if
                  ist(1,nstrtor) = i
                  do k = 1, n12(ia)
                     if (i12(k,ia) .eq. ib) then
                        ist(2,nstrtor) = bndlist(k,ia)
                        goto 80
                     endif
                  end do
   80             continue
                  do k = 1, n12(ib)
                     if (i12(k,ib) .eq. ic) then
                        ist(3,nstrtor) = bndlist(k,ib)
                        goto 90
                     end if
                  end do
   90             continue
                  do k = 1, n12(ic)
                     if (i12(k,ic) .eq. id) then
                        ist(4,nstrtor) = bndlist(k,ic)
                        goto 100
                     end if
                  end do
               end if
            end do
  100       continue
         end do
      end if
c
c     turn off the stretch-torsion potential if it is not used
c
      if (nstrtor .eq. 0)  use_strtor = .false.
      return
      end
