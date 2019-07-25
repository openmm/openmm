c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kvdw  --  van der Waals parameter assignment  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kvdw" assigns the parameters to be used in computing the
c     van der Waals interactions and processes any new or changed
c     values for these parameters
c
c
      subroutine kvdw
      use sizes
      use atomid
      use atoms
      use couple
      use fields
      use inform
      use iounit
      use keys
      use khbond
      use kvdws
      use kvdwpr
      use math
      use merck
      use potent
      use vdw
      use vdwpot
      implicit none
      integer i,k,ia,ib
      integer next,size
      integer number
      real*8 rd,ep,rdn,gik
      real*8, allocatable :: srad(:)
      real*8, allocatable :: srad4(:)
      real*8, allocatable :: seps(:)
      real*8, allocatable :: seps4(:)
      logical header
      character*4 pa,pb
      character*8 blank,pt
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     process keywords containing van der Waals parameters
c
      blank = '        '
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:4) .eq. 'VDW ') then
            call getnumb (record,k,next)
            if (k.ge.1 .and. k.le.maxclass) then
               rd = rad(k)
               ep = eps(k)
               rdn = reduct(k)
               string = record(next:240)
               read (string,*,err=10,end=10)  rd,ep,rdn
   10          continue
               if (header .and. .not.silent) then
                  header = .false.
                  if (vdwindex .eq. 'CLASS') then
                     write (iout,20)
   20                format (/,' Additional van der Waals Parameters :',
     &                       //,5x,'Atom Class',10x,'Size',6x,
     &                          'Epsilon',5x,'Reduction',/)
                  else
                     write (iout,30)
   30                format (/,' Additional van der Waals Parameters :',
     &                       //,5x,'Atom Type',11x,'Size',6x,
     &                          'Epsilon',5x,'Reduction',/)
                  end if
               end if
               rad(k) = rd
               eps(k) = ep
               reduct(k) = rdn
               if (.not. silent) then
                  write (iout,40)  k,rd,ep,rdn
   40             format (4x,i6,8x,2f12.4,f12.3)
               end if
            else if (k .gt. maxclass) then
               write (iout,50)  maxclass
   50          format (/,' KVDW  --  Only Atom Classes through',i4,
     &                    ' are Allowed')
               abort = .true.
            end if
         end if
      end do
c
c     process keywords containing 1-4 van der Waals parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:6) .eq. 'VDW14 ') then
            call getnumb (record,k,next)
            if (k.ge.1 .and. k.le.maxclass) then
               rd = rad4(k)
               ep = eps4(k)
               string = record(next:240)
               read (string,*,err=60,end=60)  rd,ep
   60          continue
               if (header .and. .not.silent) then
                  header = .false.
                  if (vdwindex .eq. 'CLASS') then
                     write (iout,70)
   70                format (/,' Additional 1-4 van der Waals',
     &                          ' Parameters :',
     &                       //,5x,'Atom Class',10x,'Size',6x,
     &                          'Epsilon',/)
                  else
                     write (iout,80)
   80                format (/,' Additional 1-4 van der Waals',
     &                          ' Parameters :',
     &                       //,5x,'Atom Type',11x,'Size',6x,
     &                          'Epsilon',/)
                  end if
               end if
               rad4(k) = rd
               eps4(k) = ep
               if (.not. silent) then
                  write (iout,90)  k,rd,ep
   90             format (4x,i6,8x,2f12.4)
               end if
            else if (k .gt. maxclass) then
               write (iout,100)  maxclass
  100          format (/,' KVDW  --  Only Atom Classes through',i4,
     &                    ' are Allowed')
               abort = .true.
            end if
         end if
      end do
c
c     process keywords containing specific pair vdw parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:6) .eq. 'VDWPR ') then
            ia = 0
            ib = 0
            rd = 0.0d0
            ep = 0.0d0
            string = record(next:240)
            read (string,*,err=150,end=150)  ia,ib,rd,ep
            if (header .and. .not.silent) then
               header = .false.
               if (vdwindex .eq. 'CLASS') then
                  write (iout,110)
  110             format (/,' Additional van der Waals Parameters',
     &                       ' for Specific Pairs :',
     &                    //,5x,'Atom Classes',6x,'Size Sum',
     &                       4x,'Epsilon',/)
               else
                  write (iout,120)
  120             format (/,' Additional van der Waals Parameters',
     &                       ' for Specific Pairs :',
     &                    //,5x,'Atom Types',8x,'Size Sum',
     &                       4x,'Epsilon',/)
               end if
            end if
            if (.not. silent) then
               write (iout,130)  ia,ib,rd,ep
  130          format (6x,2i4,4x,2f12.4)
            end if
            size = 4
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            if (ia .le. ib) then
               pt = pa//pb
            else
               pt = pb//pa
            end if
            do k = 1, maxnvp
               if (kvpr(k).eq.blank .or. kvpr(k).eq.pt) then
                  kvpr(k) = pt
                  radpr(k) = rd
                  epspr(k) = ep
                  goto 150
               end if
            end do
            write (iout,140)
  140       format (/,' KVDW  --  Too many Special VDW Pair',
     &                 ' Parameters')
            abort = .true.
  150       continue
         end if
      end do
c
c     process keywords containing hydrogen bonding vdw parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:6) .eq. 'HBOND ') then
            ia = 0
            ib = 0
            rd = 0.0d0
            ep = 0.0d0
            string = record(next:240)
            read (string,*,err=200,end=200)  ia,ib,rd,ep
            if (header .and. .not.silent) then
               header = .false.
               if (vdwindex .eq. 'CLASS') then
                  write (iout,160)
  160             format (/,' Additional van der Waals Hydrogen',
     &                       ' Bonding Parameters :',
     &                    //,5x,'Atom Classes',6x,'Size Sum',
     &                       4x,'Epsilon',/)
               else
                  write (iout,170)
  170             format (/,' Additional van der Waals Hydrogen',
     &                       ' Bonding Parameters :',
     &                    //,5x,'Atom Types',8x,'Size Sum',
     &                       4x,'Epsilon',/)
               end if
            end if
            if (.not. silent) then
               write (iout,180)  ia,ib,rd,ep
  180          format (6x,2i4,4x,2f12.4)
            end if
            size = 4
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            if (ia .le. ib) then
               pt = pa//pb
            else
               pt = pb//pa
            end if
            do k = 1, maxnvp
               if (khb(k).eq.blank .or. khb(k).eq.pt) then
                  khb(k) = pt
                  radhb(k) = rd
                  epshb(k) = ep
                  goto 200
               end if
            end do
            write (iout,190)
  190       format (/,' KVDW  --  Too many Hydrogen Bonding Pair',
     &                 ' Parameters')
            abort = .true.
  200       continue
         end if
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(ivdw))  deallocate (ivdw)
      if (allocated(jvdw))  deallocate (jvdw)
      if (allocated(ired))  deallocate (ired)
      if (allocated(kred))  deallocate (kred)
      if (allocated(radmin))  deallocate (radmin)
      if (allocated(epsilon))  deallocate (epsilon)
      if (allocated(radmin4))  deallocate (radmin4)
      if (allocated(epsilon4))  deallocate (epsilon4)
      if (allocated(radhbnd))  deallocate (radhbnd)
      if (allocated(epshbnd))  deallocate (epshbnd)
      allocate (ivdw(n))
      allocate (jvdw(n))
      allocate (ired(n))
      allocate (kred(n))
      allocate (radmin(maxclass,maxclass))
      allocate (epsilon(maxclass,maxclass))
      allocate (radmin4(maxclass,maxclass))
      allocate (epsilon4(maxclass,maxclass))
      allocate (radhbnd(maxclass,maxclass))
      allocate (epshbnd(maxclass,maxclass))
c
c     use atom class or type as index into vdw parameters
c
      k = 0
      do i = 1, n
         jvdw(i) = class(i)
         if (vdwindex .eq. 'TYPE')  jvdw(i) = type(i)
         k = max(k,jvdw(i))
      end do
      if (vdwindex.eq.'CLASS' .and. k.gt.maxclass) then
         write (iout,210)
  210    format (/,' KVDW  --  Unable to Index VDW Parameters;',
     &              ' Increase MAXCLASS')
         abort = .true.
      else if (vdwindex.eq.'TYPE' .and. k.gt.maxtyp) then
         write (iout,220)
  220    format (/,' KVDW  --  Unable to Index VDW Parameters;',
     &              ' Increase MAXTYP')
         abort = .true.
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (srad(maxtyp))
      allocate (srad4(maxtyp))
      allocate (seps(maxtyp))
      allocate (seps4(maxtyp))
c
c     get the vdw radii and well depths for each atom type
c
      do i = 1, maxtyp
         if (rad4(i) .eq. 0.0d0)  rad4(i) = rad(i)
         if (eps4(i) .eq. 0.0d0)  eps4(i) = eps(i)
         if (radtyp .eq. 'SIGMA') then
            rad(i) = twosix * rad(i)
            rad4(i) = twosix * rad4(i)
         end if
         if (radsiz .eq. 'DIAMETER') then
            rad(i) = 0.5d0 * rad(i)
            rad4(i) = 0.5d0 * rad4(i)
         end if
         srad(i) = sqrt(rad(i))
         eps(i) = abs(eps(i))
         seps(i) = sqrt(eps(i))
         srad4(i) = sqrt(rad4(i))
         eps4(i) = abs(eps4(i))
         seps4(i) = sqrt(eps4(i))
      end do
c
c     use combination rules to set pairwise vdw radii sums
c
      do i = 1, maxclass
         do k = i, maxclass
            if (radrule(1:6) .eq. 'MMFF94') then
               if (i .ne. k) then
                  if (DA(i).eq.'D' .or. DA(k).eq.'D') then
                     rd = 0.5d0 * (rad(i)+rad(k))
                  else
                     gik = (rad(i)-rad(k))/(rad(i)+rad(k))
                     rd = 0.5d0 * (rad(i)+rad(k))
     &                     * (1.0d0+0.2d0*(1.0d0-exp(-12.0d0*gik*gik)))
                  end if
               else
                  rd = rad(i)
               end if
            else if (rad(i).eq.0.0d0 .and. rad(k).eq.0.0d0) then
               rd = 0.0d0
            else if (radrule(1:10) .eq. 'ARITHMETIC') then
               rd = rad(i) + rad(k)
            else if (radrule(1:9) .eq. 'GEOMETRIC') then
               rd = 2.0d0 * (srad(i) * srad(k))
            else if (radrule(1:10) .eq. 'CUBIC-MEAN') then
               rd = 2.0d0 * (rad(i)**3+rad(k)**3)/(rad(i)**2+rad(k)**2)
            else
               rd = rad(i) + rad(k)
            end if
            radmin(i,k) = rd
            radmin(k,i) = rd
         end do
      end do
c
c     use combination rules to set pairwise well depths
c
      do i = 1, maxclass
         do k = i, maxclass
            if (epsrule(1:6) .eq. 'MMFF94') then
               ep = 181.16d0*G(i)*G(k)*alph(i)*alph(k)
     &                 / ((sqrt(alph(i)/Nn(i))+sqrt(alph(k)/Nn(k)))
     &                              *radmin(i,k)**6)
               if (i .eq. k)  eps(i) = ep
            else if (eps(i).eq.0.0d0 .and. eps(k).eq.0.0d0) then
               ep = 0.0d0
            else if (epsrule(1:10) .eq. 'ARITHMETIC') then
               ep = 0.5d0 * (eps(i) + eps(k))
            else if (epsrule(1:9) .eq. 'GEOMETRIC') then
               ep = seps(i) * seps(k)
            else if (epsrule(1:8) .eq. 'HARMONIC') then
               ep = 2.0d0 * (eps(i)*eps(k)) / (eps(i)+eps(k))
            else if (epsrule(1:3) .eq. 'HHG') then
               ep = 4.0d0 * (eps(i)*eps(k)) / (seps(i)+seps(k))**2
            else if (epsrule(1:3) .eq. 'W-H') then
               ep = 2.0D0 * seps(i) * seps(k) * (rad(i)*rad(k))**3 /
     &              (rad(i)**6 + rad(k)**6)
            else
               ep = seps(i) * seps(k)
            end if
            epsilon(i,k) = ep
            epsilon(k,i) = ep
         end do
      end do
c
c     use combination rules to set pairwise 1-4 vdw radii sums
c
      do i = 1, maxclass
         do k = i, maxclass
            if (radrule(1:6) .eq. 'MMFF94') then
               if (i .ne. k) then
                  if (DA(i).eq.'D' .or. DA(k).eq.'D') then
                     rd = 0.5d0 * (rad(i)+rad(k))
                  else
                     gik = (rad(i)-rad(k))/(rad(i)+rad(k))
                     rd = 0.5d0 * (rad(i)+rad(k))
     &                     * (1.0d0+0.2d0*(1.0d0-exp(-12.0d0*gik*gik)))
                  end if
               else
                  rd = rad(i)
               end if
            else if (rad4(i).eq.0.0d0 .and. rad4(k).eq.0.0d0) then
               rd = 0.0d0
            else if (radrule(1:10) .eq. 'ARITHMETIC') then
               rd = rad4(i) + rad4(k)
            else if (radrule(1:9) .eq. 'GEOMETRIC') then
               rd = 2.0d0 * (srad4(i) * srad4(k))
            else if (radrule(1:10) .eq. 'CUBIC-MEAN') then
               rd = 2.0d0 * (rad4(i)**3+rad4(k)**3)
     &                         / (rad4(i)**2+rad4(k)**2)
            else
               rd = rad4(i) + rad4(k)
            end if
            radmin4(i,k) = rd
            radmin4(k,i) = rd
         end do
      end do
c
c     use combination rules to set pairwise 1-4 well depths
c
      do i = 1, maxclass
         do k = i, maxclass
            if (epsrule(1:6) .eq. 'MMFF94') then
               ep = 181.16d0*G(i)*G(k)*alph(i)*alph(k)
     &                 / ((sqrt(alph(i)/Nn(i))+sqrt(alph(k)/Nn(k)))
     &                              *radmin(i,k)**6)
               if (i .eq. k)  eps4(i) = ep
            else if (eps4(i).eq.0.0d0 .and. eps4(k).eq.0.0d0) then
               ep = 0.0d0
            else if (epsrule(1:10) .eq. 'ARITHMETIC') then
               ep = 0.5d0 * (eps4(i) + eps4(k))
            else if (epsrule(1:9) .eq. 'GEOMETRIC') then
               ep = seps4(i) * seps4(k)
            else if (epsrule(1:8) .eq. 'HARMONIC') then
               ep = 2.0d0 * (eps4(i)*eps4(k)) / (eps4(i)+eps4(k))
            else if (epsrule(1:3) .eq. 'HHG') then
               ep = 4.0d0 * (eps4(i)*eps4(k)) / (seps4(i)+seps4(k))**2
            else if (epsrule(1:3) .eq. 'W-H') then
               ep = 2.0d0 * seps(i) * seps(k) * (rad(i)*rad(k))**3 /
     &              (rad(i)**6 + rad(k)**6)
            else
               ep = seps4(i) * seps4(k)
            end if
            epsilon4(i,k) = ep
            epsilon4(k,i) = ep
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (srad)
      deallocate (srad4)
      deallocate (seps)
      deallocate (seps4)
c
c     use reduced values for MMFF donor-acceptor pairs
c
      if (forcefield .eq. 'MMFF94') then
         do i = 1, maxclass
            do k = i, maxclass
               if ((da(i).eq.'D' .and. da(k).eq.'A') .or.
     &             (da(i).eq.'A' .and. da(k).eq.'D')) then
                  epsilon(i,k) = epsilon(i,k) * 0.5d0
                  epsilon(k,i) = epsilon(k,i) * 0.5d0
                  radmin(i,k) = radmin(i,k) * 0.8d0
                  radmin(k,i) = radmin(k,i) * 0.8d0
                  epsilon4(i,k) = epsilon4(i,k) * 0.5d0
                  epsilon4(k,i) = epsilon4(k,i) * 0.5d0
                  radmin4(i,k) = radmin4(i,k) * 0.8d0
                  radmin4(k,i) = radmin4(k,i) * 0.8d0
               end if
            end do
         end do
      end if
c
c     vdw reduction factor information for each individual atom
c
      do i = 1, n
         kred(i) = reduct(jvdw(i))
         if (n12(i).ne.1 .or. kred(i).eq.0.0d0) then
            ired(i) = i
         else
            ired(i) = i12(1,i)
         end if
      end do
c
c     radii and well depths for special atom class pairs
c
      do i = 1, maxnvp
         if (kvpr(i) .eq. blank)  goto 230
         ia = number(kvpr(i)(1:4))
         ib = number(kvpr(i)(5:8))
         if (rad(ia) .eq. 0.0d0)  rad(ia) = 0.001d0
         if (rad(ib) .eq. 0.0d0)  rad(ib) = 0.001d0
         if (radtyp .eq. 'SIGMA')  radpr(i) = twosix * radpr(i)
         radmin(ia,ib) = radpr(i)
         radmin(ib,ia) = radpr(i)
         epsilon(ia,ib) = abs(epspr(i))
         epsilon(ib,ia) = abs(epspr(i))
         radmin4(ia,ib) = radpr(i)
         radmin4(ib,ia) = radpr(i)
         epsilon4(ia,ib) = abs(epspr(i))
         epsilon4(ib,ia) = abs(epspr(i))
      end do
  230 continue
c
c     radii and well depths for hydrogen bonding pairs
c
      if (vdwtyp .eq. 'MM3-HBOND') then
         do i = 1, maxclass
            do k = 1, maxclass
               radhbnd(k,i) = 0.0d0
               epshbnd(k,i) = 0.0d0
            end do
         end do
         do i = 1, maxnhb
            if (khb(i) .eq. blank)  goto 240
            ia = number(khb(i)(1:4))
            ib = number(khb(i)(5:8))
            if (rad(ia) .eq. 0.0d0)  rad(ia) = 0.001d0
            if (rad(ib) .eq. 0.0d0)  rad(ib) = 0.001d0
            if (radtyp .eq. 'SIGMA')  radhb(i) = twosix * radhb(i)
            radhbnd(ia,ib) = radhb(i)
            radhbnd(ib,ia) = radhb(i)
            epshbnd(ia,ib) = abs(epshb(i))
            epshbnd(ib,ia) = abs(epshb(i))
         end do
  240    continue
      end if
c
c     set coefficients for Gaussian fit to eps=1 and radmin=1
c
      if (vdwtyp .eq. 'GAUSSIAN') then
         if (gausstyp .eq. 'LJ-4') then
            ngauss = 4
            igauss(1,1) = 846706.7d0
            igauss(2,1) = 15.464405d0 * twosix**2
            igauss(1,2) = 2713.651d0
            igauss(2,2) = 7.346875d0 * twosix**2
            igauss(1,3) = -9.699172d0
            igauss(2,3) = 1.8503725d0 * twosix**2
            igauss(1,4) = -0.7154420d0
            igauss(2,4) = 0.639621d0 * twosix**2
         else if (gausstyp .eq. 'LJ-2') then
            ngauss = 2
            igauss(1,1) = 14487.1d0
            igauss(2,1) = 9.05148d0 * twosix**2
            igauss(1,2) = -5.55338d0
            igauss(2,2) = 1.22536d0 * twosix**2
         else if (gausstyp .eq. 'MM3-2') then
            ngauss = 2
            igauss(1,1) = 2438.886d0
            igauss(2,1) = 9.342616d0
            igauss(1,2) = -6.197368d0
            igauss(2,2) = 1.564486d0
         else if (gausstyp .eq. 'MM2-2') then
            ngauss = 2
            igauss(1,1) = 3423.562d0
            igauss(2,1) = 9.692821d0
            igauss(1,2) = -6.503760d0
            igauss(2,2) = 1.585344d0
         else if (gausstyp .eq. 'IN-PLACE') then
            ngauss = 2
            igauss(1,1) = 500.0d0
            igauss(2,1) = 6.143d0
            igauss(1,2) = -18.831d0
            igauss(2,2) = 2.209d0
         end if
      end if
c
c     remove zero-sized atoms from the list of vdw sites
c
      nvdw = 0
      do i = 1, n
         if (rad(jvdw(i)) .ne. 0.0d0) then
            nvdw = nvdw + 1
            ivdw(nvdw) = i
         end if
      end do
c
c     turn off the van der Waals potential if it is not used
c
      if (nvdw .eq. 0)  use_vdw = .false.
      return
      end
