c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1994  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine kmpole  --  multipole parameter assignment  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "kmpole" assigns atomic multipole moments to the atoms of
c     the structure and processes any new or changed values
c
c
      subroutine kmpole
      use sizes
      use atoms
      use couple
      use inform
      use iounit
      use keys
      use kmulti
      use mpole
      use polar
      use polgrp
      use potent
      use units
      implicit none
      integer i,j,k,l,m
      integer ji,ki,li
      integer it,jt,kt,lt
      integer imp,nmp
      integer size,next
      integer number
      integer kz,kx,ky
      integer ztyp,xtyp,ytyp
      integer, allocatable :: mpt(:)
      integer, allocatable :: mpz(:)
      integer, allocatable :: mpx(:)
      integer, allocatable :: mpy(:)
      real*8 mpl(13)
      logical header,path
      character*4 pa,pb,pc,pd
      character*8 axt
      character*16 blank,pt
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     count the number of existing multipole parameters
c
      blank = '                '
      nmp = maxnmp
      do i = maxnmp, 1, -1
         if (kmp(i) .eq. blank)  nmp = i - 1
      end do
c
c     find and count new multipole parameters in the keyfile
c
      imp = 0
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:10) .eq. 'MULTIPOLE ') then
            k = 0
            string = record(next:240)
            read (string,*,err=10,end=10)  k,kz,kx,ky,mpl(1)
            goto 40
   10       continue
            read (string,*,err=20,end=20)  k,kz,kx,mpl(1)
            goto 40
   20       continue
            read (string,*,err=30,end=30)  k,kz,mpl(1)
            goto 40
   30       continue
            read (string,*,err=50,end=50)  k,mpl(1)
   40       continue
            if (k .gt. 0) then
               record = keyline(i+1)
               read (record,*,err=50,end=50)  mpl(2),mpl(3),mpl(4)
               record = keyline(i+2)
               read (record,*,err=50,end=50)  mpl(5)
               record = keyline(i+3)
               read (record,*,err=50,end=50)  mpl(8),mpl(9)
               record = keyline(i+4)
               read (record,*,err=50,end=50)  mpl(11),mpl(12),mpl(13)
               imp = imp + 1
            end if
   50       continue
         end if
      end do
c
c     check for too many combined parameter values
c
      nmp = nmp + imp
      if (nmp .gt. maxnmp) then
         write (iout,60)
   60    format (/,' KMPOLE  --  Too many Atomic Multipole',
     &              ' Parameters')
         abort = .true.
      end if
c
c     move existing parameters to make room for new values
c
      if (imp .ne. 0) then
         do j = nmp, imp+1, -1
            k = j - imp
            kmp(j) = kmp(k)
            mpaxis(j) = mpaxis(k)
            do m = 1, 13
               multip(m,j) = multip(m,k)
            end do
         end do
      end if
c
c     process keywords containing atomic multipole parameters
c
      imp = 0
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:10) .eq. 'MULTIPOLE ') then
            k = 0
            kz = 0
            kx = 0
            ky = 0
            axt = 'Z-then-X'
            do j = 1, 13
               mpl(j) = 0.0d0
            end do
            string = record(next:240)
            read (string,*,err=70,end=70)  k,kz,kx,ky,mpl(1)
            goto 100
   70       continue
            ky = 0
            read (string,*,err=80,end=80)  k,kz,kx,mpl(1)
            goto 100
   80       continue
            kx = 0
            read (string,*,err=90,end=90)  k,kz,mpl(1)
            goto 100
   90       continue
            kz = 0
            read (string,*,err=130,end=130)  k,mpl(1)
  100       continue
            if (k .gt. 0) then
               if (kz .eq. 0)  axt = 'None'
               if (kz.ne.0 .and. kx.eq.0)  axt = 'Z-Only'
               if (kz.lt.0 .or. kx.lt.0)  axt = 'Bisector'
               if (kx.lt.0 .and. ky.lt.0)  axt = 'Z-Bisect'
               if (max(kz,kx,ky) .lt. 0)  axt = '3-Fold'
               kz = abs(kz)
               kx = abs(kx)
               ky = abs(ky)
               record = keyline(i+1)
               read (record,*,err=130,end=130)  mpl(2),mpl(3),mpl(4)
               record = keyline(i+2)
               read (record,*,err=130,end=130)  mpl(5)
               record = keyline(i+3)
               read (record,*,err=130,end=130)  mpl(8),mpl(9)
               record = keyline(i+4)
               read (record,*,err=130,end=130)  mpl(11),mpl(12),mpl(13)
               mpl(6) = mpl(8)
               mpl(7) = mpl(11)
               mpl(10) = mpl(12)
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,110)
  110             format (/,' Additional Atomic Multipole Parameters :',
     &                    //,5x,'Atom Type',5x,'Coordinate Frame',
     &                       ' Definition',8x,'Multipole Moments')
               end if
               if (.not. silent) then
                  write (iout,120)  k,kz,kx,ky,axt,(mpl(j),j=1,5),
     &                             mpl(8),mpl(9),(mpl(j),j=11,13)
  120             format (/,4x,i6,5x,i6,1x,i6,1x,i6,3x,a8,2x,f9.5,
     &                       /,48x,3f9.5,/,48x,f9.5,
     &                       /,48x,2f9.5,/,48x,3f9.5)
               end if
               size = 4
               call numeral (k,pa,size)
               call numeral (kz,pb,size)
               call numeral (kx,pc,size)
               call numeral (ky,pd,size)
               pt = pa//pb//pc//pd
               imp = imp + 1
               kmp(imp) = pt
               mpaxis(imp) = axt
               do j = 1, 13
                  multip(j,imp) = mpl(j)
               end do
            end if
  130       continue
         end if
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(ipole))  deallocate (ipole)
      if (allocated(polsiz))  deallocate (polsiz)
      if (allocated(pollist))  deallocate (pollist)
      if (allocated(zaxis))  deallocate (zaxis)
      if (allocated(xaxis))  deallocate (xaxis)
      if (allocated(yaxis))  deallocate (yaxis)
      if (allocated(pole))  deallocate (pole)
      if (allocated(rpole))  deallocate (rpole)
      if (allocated(polaxe))  deallocate (polaxe)
      if (allocated(np11))  deallocate (np11)
      if (allocated(np12))  deallocate (np12)
      if (allocated(np13))  deallocate (np13)
      if (allocated(np14))  deallocate (np14)
      allocate (ipole(n))
      allocate (polsiz(n))
      allocate (pollist(n))
      allocate (zaxis(n))
      allocate (xaxis(n))
      allocate (yaxis(n))
      allocate (pole(maxpole,n))
      allocate (rpole(maxpole,n))
      allocate (polaxe(n))
      allocate (np11(n))
      allocate (np12(n))
      allocate (np13(n))
      allocate (np14(n))
c
c     zero out local axes, multipoles and polarization attachments
c
      npole = n
      do i = 1, n
         ipole(i) = 0
         polsiz(i) = 0
         pollist(i) = 0
         zaxis(i) = 0
         xaxis(i) = 0
         yaxis(i) = 0
         polaxe(i) = '        '
         do j = 1, 13
            pole(j,i) = 0.0d0
         end do
         np11(i) = 0
         np12(i) = 0
         np13(i) = 0
         np14(i) = 0
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (mpt(maxnmp))
      allocate (mpz(maxnmp))
      allocate (mpx(maxnmp))
      allocate (mpy(maxnmp))
c
c     store the atom types associated with each parameter
c
      do i = 1, nmp
         mpt(i) = number(kmp(i)(1:4))
         mpz(i) = number(kmp(i)(5:8))
         mpx(i) = number(kmp(i)(9:12))
         mpy(i) = number(kmp(i)(13:16))
      end do
c
c     assign multipole parameters via only 1-2 connected atoms
c
      do i = 1, n
         it = type(i)
         do imp = 1, nmp
            if (it .eq. mpt(imp)) then
               ztyp = mpz(imp)
               xtyp = mpx(imp)
               ytyp = mpy(imp)
               do j = 1, n12(i)
                  ji = i12(j,i)
                  jt = type(ji)
                  if (jt .eq. ztyp) then
                     do k = 1, n12(i)
                        ki = i12(k,i)
                        kt = type(ki)
                        if (kt.eq.xtyp .and. ki.ne.ji) then
                           if (ytyp .eq. 0) then
                              zaxis(i) = ji
                              xaxis(i) = ki
                              polaxe(i) = mpaxis(imp)
                              do m = 1, 13
                                 pole(m,i) = multip(m,imp)
                              end do
                              goto 140
                           end if
                           do l = 1, n12(i)
                              li = i12(l,i)
                              lt = type(li)
                              if (lt.eq.ytyp .and. li.ne.ji
     &                               .and. li.ne.ki) then
                                 zaxis(i) = ji
                                 xaxis(i) = ki
                                 yaxis(i) = li
                                 polaxe(i) = mpaxis(imp)
                                 do m = 1, 13
                                    pole(m,i) = multip(m,imp)
                                 end do
                                 goto 140
                              end if
                           end do
                        end if
                     end do
                  end if
               end do
            end if
         end do
c
c     assign multipole parameters via 1-2 and 1-3 connected atoms
c
         do imp = 1, nmp
            if (it .eq. mpt(imp)) then
               ztyp = mpz(imp)
               xtyp = mpx(imp)
               ytyp = mpy(imp)
               do j = 1, n12(i)
                  ji = i12(j,i)
                  jt = type(ji)
                  if (jt .eq. ztyp) then
                     do k = 1, n13(i)
                        ki = i13(k,i)
                        kt = type(ki)
                        path = .false.
                        do m = 1, n12(ki)
                           if (i12(m,ki) .eq. ji)  path = .true.
                        end do
                        if (kt.eq.xtyp .and. path) then
                           if (ytyp .eq. 0) then
                              zaxis(i) = ji
                              xaxis(i) = ki
                              polaxe(i) = mpaxis(imp)
                              do m = 1, 13
                                 pole(m,i) = multip(m,imp)
                              end do
                              goto 140
                           end if
                           do l = 1, n13(i)
                              li = i13(l,i)
                              lt = type(li)
                              path = .false.
                              do m = 1, n12(li)
                                 if (i12(m,li) .eq. ji)  path = .true.
                              end do
                              if (lt.eq.ytyp .and. li.ne.ki
     &                               .and. path) then
                                 zaxis(i) = ji
                                 xaxis(i) = ki
                                 yaxis(i) = li
                                 polaxe(i) = mpaxis(imp)
                                 do m = 1, 13
                                    pole(m,i) = multip(m,imp)
                                 end do
                                 goto 140
                              end if
                           end do
                        end if
                     end do
                  end if
               end do
            end if
         end do
c
c     assign multipole parameters via only a z-defining atom
c
         do imp = 1, nmp
            if (it .eq. mpt(imp)) then
               ztyp = mpz(imp)
               xtyp = mpx(imp)
               ytyp = mpy(imp)
               do j = 1, n12(i)
                  ji = i12(j,i)
                  jt = type(ji)
                  if (jt .eq. ztyp) then
                     if (xtyp .eq. 0) then
                        zaxis(i) = ji
                        polaxe(i) = mpaxis(imp)
                        do m = 1, 13
                           pole(m,i) = multip(m,imp)
                        end do
                        goto 140
                     end if
                  end if
               end do
            end if
         end do
c
c     assign multipole parameters via no connected atoms
c
         do imp = 1, nmp
            if (it .eq. mpt(imp)) then
               ztyp = mpz(imp)
               xtyp = mpx(imp)
               ytyp = mpy(imp)
               if (ztyp .eq. 0) then
                  polaxe(i) = mpaxis(imp)
                  do m = 1, 13
                     pole(m,i) = multip(m,imp)
                  end do
                  goto 140
               end if
            end if
         end do
  140    continue
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (mpt)
      deallocate (mpz)
      deallocate (mpx)
      deallocate (mpy)
c
c     process keywords with multipole parameters for specific atoms
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:10) .eq. 'MULTIPOLE ') then
            k = 0
            kz = 0
            kx = 0
            ky = 0
            axt = 'Z-then-X'
            do j = 1, 13
               mpl(j) = 0.0d0
            end do
            string = record(next:240)
            read (string,*,err=150,end=150)  k,kz,kx,ky,mpl(1)
            goto 180
  150       continue
            ky = 0
            read (string,*,err=160,end=160)  k,kz,kx,mpl(1)
            goto 180
  160       continue
            kx = 0
            read (string,*,err=170,end=170)  k,kz,mpl(1)
            goto 180
  170       continue
            kz = 0
            read (string,*,err=210,end=210)  k,mpl(1)
  180       continue
            if (k.lt.0 .and. k.ge.-n) then
               k = -k
               if (kz .eq. 0)  axt = 'None'
               if (kz.ne.0 .and. kx.eq.0)  axt = 'Z-Only'
               if (kz.lt.0 .or. kx.lt.0)  axt = 'Bisector'
               if (kx.lt.0 .and. ky.lt.0)  axt = 'Z-Bisect'
               if (max(kz,kx,ky) .lt. 0)  axt = '3-Fold'
               kz = abs(kz)
               kx = abs(kx)
               ky = abs(ky)
               record = keyline(i+1)
               read (record,*,err=210,end=210)  mpl(2),mpl(3),mpl(4)
               record = keyline(i+2)
               read (record,*,err=210,end=210)  mpl(5)
               record = keyline(i+3)
               read (record,*,err=210,end=210)  mpl(8),mpl(9)
               record = keyline(i+4)
               read (record,*,err=210,end=210)  mpl(11),mpl(12),mpl(13)
               mpl(6) = mpl(8)
               mpl(7) = mpl(11)
               mpl(10) = mpl(12)
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,190)
  190             format (/,' Additional Atomic Multipoles',
     &                       ' for Specific Atoms :',
     &                    //,6x,'Atom',9x,'Coordinate Frame',
     &                       ' Definition',8x,'Multipole Moments')
               end if
               if (.not. silent) then
                  write (iout,200)  k,kz,kx,ky,axt,(mpl(j),j=1,5),
     &                              mpl(8),mpl(9),(mpl(j),j=11,13)
  200             format (/,4x,i6,5x,i6,1x,i6,1x,i6,3x,a8,2x,f9.5,
     &                       /,48x,3f9.5,/,48x,f9.5,
     &                       /,48x,2f9.5,/,48x,3f9.5)
               end if
               zaxis(k) = kz
               xaxis(k) = kx
               yaxis(k) = ky
               polaxe(k) = axt
               do j = 1, 13
                  pole(j,k) = mpl(j)
               end do
            end if
  210       continue
         end if
      end do
c
c     convert the dipole and quadrupole moments to Angstroms,
c     quadrupole divided by 3 for use as traceless values
c
      do i = 1, n
         do k = 2, 4
            pole(k,i) = pole(k,i) * bohr
         end do
         do k = 5, 13
            pole(k,i) = pole(k,i) * bohr**2 / 3.0d0
         end do
      end do
c
c     get the order of the multipole expansion at each site
c
      do i = 1, n
         size = 0
         do k = 1, maxpole
            if (pole(k,i) .ne. 0.0d0)  size = max(k,size)
         end do
         if (size .gt. 4) then
            size = 13
         else if (size .gt. 1) then
            size = 4
         end if
         polsiz(i) = size
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (.not. use_polar) then
         if (allocated(uind))  deallocate (uind)
         if (allocated(uinp))  deallocate (uinp)
         if (allocated(uinds))  deallocate (uinds)
         if (allocated(uinps))  deallocate (uinps)
         allocate (uind(3,n))
         allocate (uinp(3,n))
         allocate (uinds(3,n))
         allocate (uinps(3,n))
c
c     if polarization not used, zero out induced dipoles
c
         do i = 1, n
            do j = 1, 3
               uind(j,i) = 0.0d0
               uinp(j,i) = 0.0d0
               uinds(j,i) = 0.0d0
               uinps(j,i) = 0.0d0
            end do
         end do
      end if
c
c     remove any zero or undefined atomic multipoles
c
      if (.not. use_polar) then
         npole = 0
         do i = 1, n
            if (polsiz(i) .ne. 0) then
               npole = npole + 1
               ipole(npole) = i
               pollist(i) = npole
               zaxis(npole) = zaxis(i)
               xaxis(npole) = xaxis(i)
               yaxis(npole) = yaxis(i)
               polaxe(npole) = polaxe(i)
               do j = 1, maxpole
                  pole(j,npole) = pole(j,i)
               end do
            end if
         end do
c
c     test multipoles at chiral sites and invert if necessary
c
         call chkpole
c
c     turn off the atomic multipole potential if it is not used
c
         if (npole .eq. 0)  use_mpole = .false.
      end if
      return
      end
