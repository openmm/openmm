c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine prtprm  --  output of force field parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "prtprm" writes out a formatted listing of the default
c     set of potential energy parameters for a force field
c
c
      subroutine prtprm (itxt)
      use sizes
      use angpot
      use bndpot
      use chgpot
      use fields
      use kanang
      use kangs
      use kantor
      use katoms
      use kbonds
      use kchrge
      use kdipol
      use khbond
      use kiprop
      use kitors
      use kmulti
      use kopbnd
      use kopdst
      use korbs
      use kpitor
      use kpolr
      use kstbnd
      use ksttor
      use ktorsn
      use ktrtor
      use kurybr
      use kvdws
      use kvdwpr
      use mplpot
      use polpot
      use urypot
      use vdwpot
      use ctran
      implicit none
      integer i,j,k,itxt
      integer number,npg
      integer k1,k2,k3
      integer k4,k5
      integer fold(6)
      real*8 ampli(6)
      real*8 phase(6)
      logical exist
      character*1 formfeed
      character*3 blank3
      character*8 blank8
      character*12 blank12
      character*16 blank16
      character*20 blank20
c
c
c     define blank character strings of various lengths
c
      blank3 = '   '
      blank8 = '        '
      blank12 = '            '
      blank16 = '                '
      blank20 = '                    '
c
c     set the string value of the formfeed character (Ctrl-L)
c
      formfeed = char(12)
c
c     force field atom type definitions
c
      exist = .false.
      do i = 1, maxtyp
         if (symbol(i) .ne. blank3)  exist = .true.
      end do
      if (exist) then
         write (itxt,10)  forcefield
   10    format (//,15x,'TINKER Force Field Parameters for ',a20)
         write (itxt,20)
   20    format (//,15x,'Force Field Atom Definitions',
     &           //,54x,'Atomic',4x,'Atomic',
     &           /,5x,'Type',3x,'Class',3x,'Symbol',3x,'Description',
     &              14x,'Number',4x,'Weight',3x,'Valence',/)
         do i = 1, maxtyp
            if (symbol(i) .ne. blank3) then
               write (itxt,30)  i,atmcls(i),symbol(i),describe(i),
     &                          atmnum(i),weight(i),ligand(i)
   30          format (3x,i5,3x,i5,5x,a3,5x,a24,i5,f12.3,i7)
            end if
         end do
      end if
c
c     van der Waals parameters for atom types
c
      exist = .false.
      do i = 1, maxtyp
         if (rad(i) .ne. 0.0d0)  exist = .true.
      end do
      if (exist) then
         write (itxt,40)  formfeed,forcefield
   40    format (a1,//,15x,'TINKER Force Field Parameters for ',a20)
         if (vdwindex .eq. 'CLASS') then
            write (itxt,50)
   50       format (//,15x,'Van der Waals Parameters',
     &              ///,20x,'Class',7x,'Radius',6x,'Epsilon',
     &                    4x,'Reduction',/)
         else
            write (itxt,60)
   60       format (//,15x,'Van der Waals Parameters',
     &              ///,20x,'Type',8x,'Radius',6x,'Epsilon',
     &                    4x,'Reduction',/)
         end if
         k = 0
         do i = 1, maxtyp
            if (rad(i) .ne. 0.0d0) then
               k = k + 1
               write (itxt,70)  k,i,rad(i),eps(i),reduct(i)
   70          format (10x,i5,5x,i4,2x,3f12.3)
            end if
         end do
c        Charge transfer parameters
         k = 0
         do i = 1, maxtyp
            if (apre(i) .ne. 0.0d0) then
               k = k + 1
               write (itxt,75)  k,i,apre(i),bexp(i)
   75          format (10x,i5,5x,i4,2x,2f12.3)
            end if
         end do
c
c     van der Waals scaling parameters
c
         write (itxt,80)  v2scale,v3scale,v4scale,v5scale
   80    format (//,15x,'Van der Waals Scaling Factors',
     &           ///,20x,'1-2 Atoms',f17.3,/,20x,'1-3 Atoms',f17.3,
     &           /,20x,'1-4 Atoms',f17.3,/,20x,'1-5 Atoms',f17.3)
      end if
c
c     van der Waals 1-4 parameters for atom types
c
      exist = .false.
      do i = 1, maxtyp
         if (rad4(i) .ne. 0.0d0)  exist = .true.
      end do
      if (exist) then
         if (vdwindex .eq. 'CLASS') then
            write (itxt,90)
   90       format (//,15x,'Van der Waals Parameters for 1-4',
     &                 ' Interactions',
     &              ///,20x,'Class',7x,'Radius',6x,'Epsilon',/)
         else
            write (itxt,100)
  100       format (//,15x,'Van der Waals Parameters for 1-4',
     &                 ' Interactions',
     &              ///,20x,'Type',8x,'Radius',6x,'Epsilon',/)
         end if
         k = 0
         do i = 1, maxtyp
            if (rad4(i) .ne. 0.0d0) then
               k = k + 1
               write (itxt,110)  k,i,rad4(i),eps4(i)
  110          format (10x,i5,5x,i4,2x,2f12.3)
            end if
         end do
      end if
c
c     van der Waals parameters for specific atom pairs
c
      if (kvpr(1) .ne. blank8) then
         if (vdwindex .eq. 'CLASS') then
            write (itxt,120)
  120       format (//,15x,'Van der Waals Parameters for Atom Pairs',
     &              ///,22x,'Classes',7x,'Radii Sum',4x,'Epsilon',/)
         else
            write (itxt,130)
  130       format (//,15x,'Van der Waals Parameters for Atom Pairs',
     &              ///,23x,'Types',8x,'Radii Sum',4x,'Epsilon',/)
         end if
         do i = 1, maxnvp
            if (kvpr(i) .eq. blank8)  goto 150
            k1 = number(kvpr(i)(1:4))
            k2 = number(kvpr(i)(5:8))
            write (itxt,140)  i,k1,k2,radpr(i),epspr(i)
  140       format (10x,i5,5x,i4,'-',i4,2x,2f12.3)
         end do
  150    continue
      end if
c
c     hydrogen bonding parameters for specific atom pairs
c
      if (khb(1) .ne. blank8) then
         if (vdwindex .eq. 'CLASS') then
            write (itxt,160)
  160       format (//,15x,'Hydrogen Bonding Parameters for Atom Pairs',
     &              ///,22x,'Classes',7x,'Radii Sum',4x,'Epsilon',/)
         else
            write (itxt,170)
  170       format (//,15x,'Hydrogen Bonding Parameters for Atom Pairs',
     &              ///,23x,'Types',8x,'Radii Sum',4x,'Epsilon',/)
         end if
         do i = 1, maxnhb
            if (khb(i) .eq. blank8)  goto 190
            k1 = number(khb(i)(1:4))
            k2 = number(khb(i)(5:8))
            write (itxt,180)  i,k1,k2,radhb(i),epshb(i)
  180       format (10x,i5,5x,i4,'-',i4,2x,2f12.3)
         end do
  190    continue
      end if
c
c     bond stretching parameters
c
      if (kb(1) .ne. blank8) then
         write (itxt,200)  formfeed,forcefield
  200    format (a1,//,15x,'TINKER Force Field Parameters for ',a20)
         write (itxt,210)
  210    format (//,15x,'Bond Stretching Parameters',
     &           ///,22x,'Classes',15x,'KS',7x,'Length',/)
         do i = 1, maxnb
            if (kb(i) .eq. blank8)  goto 230
            k1 = number(kb(i)(1:4))
            k2 = number(kb(i)(5:8))
            write (itxt,220)  i,k1,k2,bcon(i),blen(i)
  220       format (10x,i5,5x,i4,'-',i4,6x,f12.3,f12.4)
         end do
  230    continue
      end if
c
c     bond stretching parameters for 5-membered rings
c
      if (kb5(1) .ne. blank8) then
         write (itxt,240)
  240    format (//,15x,'5-Membered Ring Stretch Parameters',
     &           ///,22x,'Classes',15x,'KS',7x,'Length',/)
         do i = 1, maxnb5
            if (kb5(i) .eq. blank8)  goto 260
            k1 = number(kb5(i)(1:4))
            k2 = number(kb5(i)(5:8))
            write (itxt,250)  i,k1,k2,bcon5(i),blen5(i)
  250       format (10x,i5,5x,i4,'-',i4,6x,f12.3,f12.4)
         end do
  260    continue
      end if
c
c     bond stretching parameters for 4-membered rings
c
      if (kb4(1) .ne. blank8) then
         write (itxt,270)
  270    format (//,15x,'4-Membered Ring Stretch Parameters',
     &           ///,22x,'Classes',15x,'KS',7x,'Length',/)
         do i = 1, maxnb4
            if (kb4(i) .eq. blank8)  goto 290
            k1 = number(kb4(i)(1:4))
            k2 = number(kb4(i)(5:8))
            write (itxt,280)  i,k1,k2,bcon4(i),blen4(i)
  280       format (10x,i5,5x,i4,'-',i4,6x,f12.3,f12.4)
         end do
  290    continue
      end if
c
c     bond stretching parameters for 3-membered rings
c
      if (kb3(1) .ne. blank8) then
         write (itxt,300)
  300    format (//,15x,'3-Membered Ring Stretch Parameters',
     &           ///,22x,'Classes',15x,'KS',7x,'Length',/)
         do i = 1, maxnb3
            if (kb3(i) .eq. blank8)  goto 320
            k1 = number(kb3(i)(1:4))
            k2 = number(kb3(i)(5:8))
            write (itxt,310)  i,k1,k2,bcon3(i),blen3(i)
  310       format (10x,i5,5x,i4,'-',i4,6x,f12.3,f12.4)
         end do
  320    continue
      end if
c
c     cubic and quartic bond stretching parameters
c
      if (cbnd.ne.0.0d0 .or. qbnd.ne.0.0d0) then
         write (itxt,330)  cbnd,qbnd
  330    format (//,15x,'Higher Order Stretching Constants',
     &           ///,20x,'Cubic',f17.3,/,20x,'Quartic',f15.3)
      end if
c
c     electronegativity bond length correction parameters
c
      if (kel(1) .ne. blank12) then
         write (itxt,340)
  340    format (//,15x,'Electronegativity Bond Length Parameters',
     &           ///,25x,'Classes',21x,'dLength',/)
         do i = 1, maxnel
            if (kel(i) .eq. blank12)  goto 360
            k1 = number(kel(i)(1:4))
            k2 = number(kel(i)(5:8))
            k3 = number(kel(i)(9:12))
            write (itxt,350)  i,k1,k2,k3,dlen(i)
  350       format (10x,i5,5x,i4,'-',i4,'-',i4,14x,f12.4)
         end do
  360    continue
      end if
c
c     bond angle bending parameters
c
      if (ka(1) .ne. blank12) then
         write (itxt,370)  formfeed,forcefield
  370    format (a1,//,15x,'TINKER Force Field Parameters for ',a20)
         write (itxt,380)
  380    format (//,15x,'Angle Bending Parameters',
     &           ///,18x,'Classes',11x,'KB',6x,'Value 1',
     &              5x,'Value 2',5x,'Value 3',
     &           /,44x,'(R-X-R)',5x,'(R-X-H)',5x,'(H-X-H)',/)
         do i = 1, maxna
            if (ka(i) .eq. blank12)  goto 410
            k1 = number(ka(i)(1:4))
            k2 = number(ka(i)(5:8))
            k3 = number(ka(i)(9:12))
            if (ang(2,i).eq.0.0d0 .and. ang(3,i).eq.0.0d0) then
               write (itxt,390)  i,k1,k2,k3,acon(i),ang(1,i)
  390          format (3x,i5,5x,i4,'-',i4,'-',i4,2f12.3)
            else
               write (itxt,400)  i,k1,k2,k3,acon(i),(ang(j,i),j=1,3)
  400          format (3x,i5,5x,i4,'-',i4,'-',i4,4f12.3)
            end if
         end do
  410    continue
      end if
c
c     bond angle bending parameters for 5-membered rings
c
      if (ka5(1) .ne. blank12) then
         write (itxt,420)
  420    format (//,17x,'5-Membered Ring Bend Parameters',
     &           ///,18x,'Classes',11x,'KB',6x,'Value 1',
     &              5x,'Value 2',5x,'Value 3',
     &           /,44x,'(R-X-R)',5x,'(R-X-H)',5x,'(H-X-H)',/)
         do i = 1, maxna5
            if (ka5(i) .eq. blank12)  goto 450
            k1 = number(ka5(i)(1:4))
            k2 = number(ka5(i)(5:8))
            k3 = number(ka5(i)(9:12))
            if (ang5(2,i).eq.0.0d0 .and. ang5(3,i).eq.0.0d0) then
               write (itxt,430)  i,k1,k2,k3,acon5(i),ang5(1,i)
  430          format (3x,i5,5x,i4,'-',i4,'-',i4,2f12.3)
            else
               write (itxt,440)  i,k1,k2,k3,acon5(i),(ang5(j,i),j=1,3)
  440          format (3x,i5,5x,i4,'-',i4,'-',i4,4f12.3)
            end if
         end do
  450    continue
      end if
c
c     bond angle bending parameters for 4-membered rings
c
      if (ka4(1) .ne. blank12) then
         write (itxt,460)
  460    format (//,15x,'4-Membered Ring Bend Parameters',
     &           ///,18x,'Classes',11x,'KB',6x,'Value 1',
     &              5x,'Value 2',5x,'Value 3',
     &           /,44x,'(R-X-R)',5x,'(R-X-H)',5x,'(H-X-H)',/)
         do i = 1, maxna4
            if (ka4(i) .eq. blank12)  goto 490
            k1 = number(ka4(i)(1:4))
            k2 = number(ka4(i)(5:8))
            k3 = number(ka4(i)(9:12))
            if (ang4(2,i).eq.0.0d0 .and. ang4(3,i).eq.0.0d0) then
               write (itxt,470)  i,k1,k2,k3,acon4(i),ang4(1,i)
  470          format (3x,i5,5x,i4,'-',i4,'-',i4,2f12.3)
            else
               write (itxt,480)  i,k1,k2,k3,acon4(i),(ang4(j,i),j=1,3)
  480          format (3x,i5,5x,i4,'-',i4,'-',i4,4f12.3)
            end if
         end do
  490    continue
      end if
c
c     bond angle bending parameters for 3-membered rings
c
      if (ka3(1) .ne. blank12) then
         write (itxt,500)
  500    format (//,15x,'3-Membered Ring Bend Parameters',
     &           ///,18x,'Classes',11x,'KB',6x,'Value 1',
     &              5x,'Value 2',5x,'Value 3',
     &           /,44x,'(R-X-R)',5x,'(R-X-H)',5x,'(H-X-H)',/)
         do  i = 1, maxna3
            if (ka3(i) .eq. blank12)  goto 530
            k1 = number(ka3(i)(1:4))
            k2 = number(ka3(i)(5:8))
            k3 = number(ka3(i)(9:12))
            if (ang3(2,i).eq.0.0d0 .and. ang3(3,i).eq.0.0d0) then
               write (itxt,510)  i,k1,k2,k3,acon3(i),ang3(1,i)
  510          format (3x,i5,5x,i4,'-',i4,'-',i4,2f12.3)
            else
               write (itxt,520)  i,k1,k2,k3,acon3(i),(ang3(j,i),j=1,3)
  520          format (3x,i5,5x,i4,'-',i4,'-',i4,4f12.3)
            end if
         end do
  530    continue
      end if
c
c     Fourier bond angle bending parameters
c
      if (kaf(1) .ne. blank12) then
         write (itxt,540)
  540    format (//,15x,'Fourier Angle Bending Parameters',
     &           ///,18x,'Classes',11x,'KB',8x,'Shift',6x,'Period',/)
         do  i = 1, maxnaf
            if (kaf(i) .eq. blank12)  goto 560
            k1 = number(kaf(i)(1:4))
            k2 = number(kaf(i)(5:8))
            k3 = number(kaf(i)(9:12))
            write (itxt,550)  i,k1,k2,k3,aconf(i),(angf(j,i),j=1,2)
  550       format (3x,i5,5x,i4,'-',i4,'-',i4,3f12.3)
         end do
  560    continue
      end if
c
c     cubic through sextic bond angle bending parameters
c
      if (cang.ne.0.0d0 .or. qang.ne.0.0d0 .or.
     &    pang.ne.0.0d0 .or. sang.ne.0.0d0) then
         write (itxt,570)  cang,qang,pang,sang
  570    format (//,15x,'Higher Order Bending Constants',
     &           ///,20x,'Cubic',d17.3,/,20x,'Quartic',d15.3,
     &           /,20x,'Pentic',d16.3,/,20x,'Sextic',d16.3)
      end if
c
c     stretch-bend parameters
c
      if (ksb(1) .ne. blank12) then
         write (itxt,580)  formfeed,forcefield
  580    format (a1,//,15x,'TINKER Force Field Parameters for ',a20)
         write (itxt,590)
  590    format (//,15x,'Stretch-Bend Parameters',
     &           ///,18x,'Classes',10x,'KSB1',8x,'KSB2',/)
         do i = 1, maxnsb
            if (ksb(i) .eq. blank12)  goto 610
            k1 = number(ksb(i)(1:4))
            k2 = number(ksb(i)(5:8))
            k3 = number(ksb(i)(9:12))
            write (itxt,600)  i,k1,k2,k3,stbn(1,i),stbn(2,i)
  600       format (3x,i5,5x,i4,'-',i4,'-',i4,2f12.3)
         end do
  610    continue
      end if
c
c     Urey-Bradley parameters
c
      if (ku(1) .ne. blank12) then
         write (itxt,620)  formfeed,forcefield
  620    format (a1,//,15x,'TINKER Force Field Parameters for ',a20)
         write (itxt,630)
  630    format (//,15x,'Urey-Bradley Parameters',
     &           ///,18x,'Classes',11x,'KB',6x,'Distance',/)
         do i = 1, maxnu
            if (ku(i) .eq. blank12)  goto 650
            k1 = number(ku(i)(1:4))
            k2 = number(ku(i)(5:8))
            k3 = number(ku(i)(9:12))
            write (itxt,640)  i,k1,k2,k3,ucon(i),dst13(i)
  640       format (3x,i5,5x,i4,'-',i4,'-',i4,f12.3,f12.4)
         end do
  650    continue
      end if
c
c     cubic and quartic Urey-Bradley parameters
c
      if (cury.ne.0.0d0 .or. qury.ne.0.0d0) then
         write (itxt,660)  cury,qury
  660    format (//,15x,'Higher Order Urey-Bradley Constants',
     &           ///,20x,'Cubic',f17.3,/,20x,'Quartic',f15.3)
      end if
c
c     angle-angle parameters
c
      exist = .false.
      do i = 1, maxclass
         do k = 1, 3
            if (anan(k,i) .ne. 0.0d0)  exist = .true.
         end do
      end do
      if (exist) then
         write (itxt,670)  formfeed,forcefield
  670    format (a1,//,15x,'TINKER Force Field Parameters for ',a20)
         write (itxt,680)
  680    format (//,15x,'Angle-Angle Parameters',
     &           ///,20x,'Class',9x,'KAA 1',7x,'KAA 2',7x,'KAA 3',
     &           /,33x,'(R-X-R)',5x,'(R-X-H)',5x,'(H-X-H)',/)
         k = 0
         do i = 1, maxclass
            if (anan(1,i).ne.0.0d0 .or. anan(2,i).ne.0.0d0
     &               .or. anan(3,i).ne.0.0d0) then
               k = k + 1
               write (itxt,690)  k,i,(anan(j,i),j=1,3)
  690          format (8x,i5,7x,i4,3x,3f12.3)
            end if
         end do
      end if
c
c     out-of-plane bending parameters
c
      if (kopb(1) .ne. blank16) then
         write (itxt,700)  formfeed,forcefield
  700    format (a1,//,15x,'TINKER Force Field Parameters for ',a20)
         write (itxt,710)
  710    format (//,15x,'Out-of-Plane Bend Parameters',
     &           ///,26x,'Classes',11x,'KOPB',/)
         do i = 1, maxnopb
            if (kopb(i) .eq. blank16)  goto 730
            k1 = number(kopb(i)(1:4))
            k2 = number(kopb(i)(5:8))
            k3 = number(kopb(i)(9:12))
            k4 = number(kopb(i)(13:16))
            write (itxt,720)  i,k1,k2,k3,k4,opbn(i)
  720       format (8x,i5,5x,i4,'-',i4,'-',i4,'-',i4,f12.3)
         end do
  730    continue
      end if
c
c     out-of-plane distance parameters
c
      if (kopd(1) .ne. blank16) then
         write (itxt,740)  formfeed,forcefield
  740    format (a1,//,15x,'TINKER Force Field Parameters for ',a20)
         write (itxt,750)
  750    format (//,15x,'Out-of-Plane Distance Parameters',
     &           ///,26x,'Classes',11x,'KOPD',/)
         do i = 1, maxnopd
            if (kopd(i) .eq. blank16)  goto 770
            k1 = number(kopd(i)(1:4))
            k2 = number(kopd(i)(5:8))
            k3 = number(kopd(i)(9:12))
            k4 = number(kopd(i)(13:16))
            write (itxt,760)  i,k1,k2,k3,k4,opds(i)
  760       format (8x,i5,5x,i4,'-',i4,'-',i4,'-',i4,f12.3)
         end do
  770    continue
      end if
c
c     improper dihedral parameters
c
      if (kdi(1) .ne. blank16) then
         write (itxt,780)  formfeed,forcefield
  780    format (a1,//,15x,'TINKER Force Field Parameters for ',a20)
         write (itxt,790)
  790    format (//,15x,'Improper Dihedral Parameters',
     &           ///,20x,'Classes',12x,'KID',7x,'Target',/)
         do i = 1, maxndi
            if (kdi(i) .eq. blank16)  goto 810
            k1 = number(kdi(i)(1:4))
            k2 = number(kdi(i)(5:8))
            k3 = number(kdi(i)(9:12))
            k4 = number(kdi(i)(13:16))
            write (itxt,800)  i,k1,k2,k3,k4,dcon(i),tdi(i)
  800       format (2x,i5,5x,i4,'-',i4,'-',i4,'-',i4,f12.3,f12.4)
         end do
  810    continue
      end if
c
c     improper torsional parameters
c
      if (kti(1) .ne. blank16) then
         write (itxt,820)  formfeed,forcefield
  820    format (a1,//,15x,'TINKER Force Field Parameters for ',a20)
         write (itxt,830)
  830    format (//,15x,'Improper Torsion Parameters',
     &           ///,17x,'Classes',15x,'KTI Values',/)
         do i = 1, maxnti
            if (kti(i) .eq. blank16)  goto 850
            k1 = number(kti(i)(1:4))
            k2 = number(kti(i)(5:8))
            k3 = number(kti(i)(9:12))
            k4 = number(kti(i)(13:16))
            j = 0
            if (ti1(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 1
               ampli(j) = ti1(1,i)
               phase(j) = ti1(2,i)
            end if
            if (ti2(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 2
               ampli(j) = ti2(1,i)
               phase(j) = ti2(2,i)
            end if
            if (ti3(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 3
               ampli(j) = ti3(1,i)
               phase(j) = ti3(2,i)
            end if
            write (itxt,840)  i,k1,k2,k3,k4,(ampli(k),
     &                        phase(k),fold(k),k=1,j)
  840       format (2x,i5,2x,i4,'-',i4,'-',i4,'-',i4,2x,3(f8.3,f6.1,i2))
         end do
  850    continue
      end if
c
c     torsional angle parameters
c
      if (kt(1) .ne. blank16) then
         write (itxt,860)  formfeed,forcefield
  860    format (a1,//,15x,'TINKER Force Field Parameters for ',a20)
         write (itxt,870)
  870    format (//,15x,'Torsional Parameters',
     &           ///,17x,'Classes',15x,'KT Values',/)
         do i = 1, maxnt
            if (kt(i) .eq. blank16)  goto 890
            k1 = number(kt(i)(1:4))
            k2 = number(kt(i)(5:8))
            k3 = number(kt(i)(9:12))
            k4 = number(kt(i)(13:16))
            j = 0
            if (t1(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 1
               ampli(j) = t1(1,i)
               phase(j) = t1(2,i)
            end if
            if (t2(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 2
               ampli(j) = t2(1,i)
               phase(j) = t2(2,i)
            end if
            if (t3(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 3
               ampli(j) = t3(1,i)
               phase(j) = t3(2,i)
            end if
            if (t4(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 4
               ampli(j) = t4(1,i)
               phase(j) = t4(2,i)
            end if
            if (t5(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 5
               ampli(j) = t5(1,i)
               phase(j) = t5(2,i)
            end if
            if (t6(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 6
               ampli(j) = t6(1,i)
               phase(j) = t6(2,i)
            end if
            write (itxt,880)  i,k1,k2,k3,k4,(ampli(k),
     &                        phase(k),fold(k),k=1,j)
  880       format (2x,i5,2x,i4,'-',i4,'-',i4,'-',i4,2x,6(f8.3,f6.1,i2))
         end do
  890    continue
      end if
c
c     torsional angle parameters for 5-membered rings
c
      if (kt5(1) .ne. blank16) then
         write (itxt,900)
  900    format (//,15x,'5-Membered Ring Torsion Parameters',
     &           ///,17x,'Classes',15x,'KT Values',/)
         do i = 1, maxnt5
            if (kt5(i) .eq. blank16)  goto 920
            k1 = number(kt5(i)(1:4))
            k2 = number(kt5(i)(5:8))
            k3 = number(kt5(i)(9:12))
            k4 = number(kt5(i)(13:16))
            j = 0
            if (t15(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 1
               ampli(j) = t15(1,i)
               phase(j) = t15(2,i)
            end if
            if (t25(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 2
               ampli(j) = t25(1,i)
               phase(j) = t25(2,i)
            end if
            if (t35(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 3
               ampli(j) = t35(1,i)
               phase(j) = t35(2,i)
            end if
            if (t45(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 4
               ampli(j) = t45(1,i)
               phase(j) = t45(2,i)
            end if
            if (t55(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 5
               ampli(j) = t55(1,i)
               phase(j) = t55(2,i)
            end if
            if (t65(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 6
               ampli(j) = t65(1,i)
               phase(j) = t65(2,i)
            end if
            write (itxt,910)  i,k1,k2,k3,k4,(ampli(k),
     &                        phase(k),fold(k),k=1,j)
  910       format (2x,i5,2x,i4,'-',i4,'-',i4,'-',i4,2x,6(f8.3,f6.1,i2))
         end do
  920    continue
      end if
c
c     torsional angle parameters for 4-membered rings
c
      if (kt4(1) .ne. blank16) then
         write (itxt,930)
  930    format (//,15x,'4-Membered Ring Torsion Parameters',
     &           ///,17x,'Classes',15x,'KT Values',/)
         do i = 1, maxnt4
            if (kt4(i) .eq. blank16)  goto 950
            k1 = number(kt4(i)(1:4))
            k2 = number(kt4(i)(5:8))
            k3 = number(kt4(i)(9:12))
            k4 = number(kt4(i)(13:16))
            j = 0
            if (t14(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 1
               ampli(j) = t14(1,i)
               phase(j) = t14(2,i)
            end if
            if (t24(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 2
               ampli(j) = t24(1,i)
               phase(j) = t24(2,i)
            end if
            if (t34(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 3
               ampli(j) = t34(1,i)
               phase(j) = t34(2,i)
            end if
            if (t44(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 4
               ampli(j) = t44(1,i)
               phase(j) = t44(2,i)
            end if
            if (t54(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 5
               ampli(j) = t54(1,i)
               phase(j) = t54(2,i)
            end if
            if (t64(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 6
               ampli(j) = t64(1,i)
               phase(j) = t64(2,i)
            end if
            write (itxt,940)  i,k1,k2,k3,k4,(ampli(k),
     &                        phase(k),fold(k),k=1,j)
  940       format (2x,i5,2x,i4,'-',i4,'-',i4,'-',i4,2x,6(f8.3,f6.1,i2))
         end do
  950    continue
      end if
c
c     pi-system torsion parameters
c
      if (kpt(1) .ne. blank8) then
         write (itxt,960)  formfeed,forcefield
  960    format (a1,//,15x,'TINKER Force Field Parameters for ',a20)
         write (itxt,970)
  970    format (//,15x,'Pi-Orbital Torsion Parameters',
     &           ///,18x,'Classes',15x,'KPT',/)
         do i = 1, maxnpt
            if (kpt(i) .eq. blank8)  goto 990
            k1 = number(kpt(i)(1:4))
            k2 = number(kpt(i)(5:8))
            write (itxt,980)  i,k1,k2,ptcon(i)
  980       format (6x,i5,5x,i4,'-',i4,6x,f12.3)
         end do
  990    continue
      end if
c
c     stretch-torsion parameters
c
      if (kbt(1) .ne. blank16) then
         write (itxt,1000)  formfeed,forcefield
 1000    format (a1,//,15x,'TINKER Force Field Parameters for ',a20)
         write (itxt,1010)
 1010    format (//,15x,'Stretch-Torsion Parameters',
     &           ///,17x,'Classes',12x,'Bond',8x,'KST1',
     &              8x,'KST2',8x,'KST3',/)
         do i = 1, maxnbt
            if (kbt(i) .eq. blank16)  goto 1030
            k1 = number(kbt(i)(1:4))
            k2 = number(kbt(i)(5:8))
            k3 = number(kbt(i)(9:12))
            k4 = number(kbt(i)(13:16))
            write (itxt,1020)  i,k1,k2,k3,k4,(btcon(j,i),j=1,9)
 1020       format (2x,i5,2x,i4,'-',i4,'-',i4,'-',i4,9x,'1st',3f12.3,
     &              /,37x,'2nd',3f12.3,/,37x,'3rd',3f12.3)
         end do
 1030    continue
      end if
c
c     angle-torsion parameters
c
      if (kat(1) .ne. blank16) then
         write (itxt,1040)  formfeed,forcefield
 1040    format (a1,//,15x,'TINKER Force Field Parameters for ',a20)
         write (itxt,1050)
 1050    format (//,15x,'Angle-Torsion Parameters',
     &           ///,17x,'Classes',12x,'Angle',7x,'KAT1',
     &              8x,'KAT2',8x,'KAT3',/)
         do i = 1, maxnat
            if (kat(i) .eq. blank16)  goto 1070
            k1 = number(kat(i)(1:4))
            k2 = number(kat(i)(5:8))
            k3 = number(kat(i)(9:12))
            k4 = number(kat(i)(13:16))
            write (itxt,1060)  i,k1,k2,k3,k4,(atcon(j,i),j=1,6)
 1060       format (2x,i5,2x,i4,'-',i4,'-',i4,'-',i4,9x,'1st',3f12.3
     &              /,37x,'2nd',3f12.3)
         end do
 1070    continue
      end if
c
c     torsion-torsion parameters
c
      if (ktt(1) .ne. blank20) then
         write (itxt,1080)  formfeed,forcefield
 1080    format (a1,//,15x,'TINKER Force Field Parameters for ',a20)
         write (itxt,1090)
 1090    format (//,15x,'Torsion-Torsion Parameters',
     &           ///,19x,'Classes',18x,'KNX',9x,'KNY')
         do i = 1, maxntt
            if (ktt(i) .eq. blank20)  goto 1120
            k1 = number(ktt(i)(1:4))
            k2 = number(ktt(i)(5:8))
            k3 = number(ktt(i)(9:12))
            k4 = number(ktt(i)(13:16))
            k5 = number(ktt(i)(17:20))
            write (itxt,1100)  i,k1,k2,k3,k4,k5,tnx(i),tny(i)
 1100       format (/,2x,i5,2x,i4,'-',i4,'-',i4,'-',i4,'-',i4,2x,2i12,/)
            k = tnx(i) * tny(i)
            write (itxt,1110)  (tbf(j,i),j=1,k)
 1110       format (3x,6f12.4)
         end do
 1120    continue
      end if
c
c     atomic partial charge parameters
c
      exist = .false.
      do i = 1, maxtyp
         if (chg(i) .ne. 0.0d0)  exist = .true.
      end do
      if (exist) then
         write (itxt,1130)  formfeed,forcefield
 1130    format (a1,//,15x,'TINKER Force Field Parameters for ',a20)
         write (itxt,1140)
 1140    format (//,15x,'Atomic Partial Charge Parameters',
     &           ///,24x,'Type',9x,'Partial Chg',/)
         k = 0
         do i = 1, maxtyp
            if (chg(i) .ne. 0.0d0) then
               k = k + 1
               write (itxt,1150)  k,i,chg(i)
 1150          format (12x,i5,7x,i3,6x,f12.3)
            end if
         end do
c
c     atomic partial charge scaling parameters
c
         write (itxt,1160)  c2scale,c3scale,c4scale,c5scale
 1160    format (//,15x,'Atomic Partial Charge Scaling Factors',
     &           ///,20x,'1-2 Atoms',f17.3,/,20x,'1-3 Atoms',f17.3,
     &           /,20x,'1-4 Atoms',f17.3,/,20x,'1-5 Atoms',f17.3)
      end if
c
c     bond dipole moment parameters
c
      if (kd(1) .ne. blank8) then
         write (itxt,1170)  formfeed,forcefield
 1170    format (a1,//,15x,'TINKER Force Field Parameters for ',a20)
         write (itxt,1180)
 1180    format (//,15x,'Bond Dipole Moment Parameters',
     &           ///,25x,'Types',10x,'Bond Dipole',4x,'Position',/)
         do i = 1, maxnd
            if (kd(i) .eq. blank8)  goto 1200
            k1 = number(kd(i)(1:4))
            k2 = number(kd(i)(5:8))
            write (itxt,1190)  i,k1,k2,dpl(i),pos(i)
 1190       format (12x,i5,5x,i4,'-',i4,6x,2f12.3)
         end do
 1200    continue
      end if
c
c     bond dipole moment parameters for 5-membered rings
c
      if (kd5(1) .ne. blank8) then
         write (itxt,1210)
 1210    format (//,15x,'5-Membered Ring Bond Dipole Parameters',
     &           ///,25x,'Types',10x,'Bond Dipole',4x,'Position',/)
         do i = 1, maxnd5
            if (kd5(i) .eq. blank8)  goto 1230
            k1 = number(kd5(i)(1:4))
            k2 = number(kd5(i)(5:8))
            write (itxt,1220)  i,k1,k2,dpl5(i),pos5(i)
 1220       format (12x,i5,5x,i4,'-',i4,6x,2f12.3)
         end do
 1230    continue
      end if
c
c     bond dipole moment parameters for 4-membered rings
c
      if (kd4(1) .ne. blank8) then
         write (itxt,1240)
 1240    format (//,15x,'4-Membered Ring Bond Dipole Parameters',
     &           ///,25x,'Types',10x,'Bond Dipole',4x,'Position',/)
         do i = 1, maxnd4
            if (kd4(i) .eq. blank8)  goto 1260
            k1 = number(kd4(i)(1:4))
            k2 = number(kd4(i)(5:8))
            write (itxt,1250)  i,k1,k2,dpl4(i),pos4(i)
 1250       format (12x,i5,5x,i4,'-',i4,6x,2f12.3)
         end do
 1260    continue
      end if
c
c     bond dipole moment parameters for 3-membered rings
c
      if (kd3(1) .ne. blank8) then
         write (itxt,1270)
 1270    format (//,15x,'3-Membered Ring Bond Dipole Parameters',
     &           ///,25x,'Types',10x,'Bond Dipole',4x,'Position',/)
         do i = 1, maxnd3
            if (kd3(i) .eq. blank8)  goto 1290
            k1 = number(kd3(i)(1:4))
            k2 = number(kd3(i)(5:8))
            write (itxt,1280)  i,k1,k2,dpl3(i),pos3(i)
 1280       format (12x,i5,5x,i4,'-',i4,6x,2f12.3)
         end do
 1290    continue
      end if
c
c     atomic multipole electrostatic parameters
c
      if (kmp(1) .ne. blank16) then
         write (itxt,1300)  formfeed,forcefield
 1300    format (a1,//,15x,'TINKER Force Field Parameters for ',a20)
         write (itxt,1310)
 1310    format (//,17x,'Atomic Multipole Parameters',
     &           ///,11x,'Type',7x,'Axis Types',8x,'Frame',
     &              9x,'Multipoles (M-D-Q)',/)
         do i = 1, maxnmp
            if (kmp(i) .eq. blank16)  goto 1330
            k1 = number(kmp(i)(1:4))
            k2 = number(kmp(i)(5:8))
            k3 = number(kmp(i)(9:12))
            k4 = number(kmp(i)(13:16))
            write (itxt,1320)  i,k1,k2,k3,k4,mpaxis(i),multip(1,i),
     &                         multip(2,i),multip(3,i),multip(4,i),
     &                         multip(5,i),multip(8,i),multip(9,i),
     &                         multip(11,i),multip(12,i),multip(13,i)
 1320       format (2x,i5,3x,i4,3x,i4,2x,i4,2x,i4,5x,a8,2x,f10.5,
     &                 /,48x,3f10.5,/,48x,f10.5,
     &                 /,48x,2f10.5,/,48x,3f10.5)
         end do
 1330    continue
c
c     atomic multipole scaling parameters
c
         write (itxt,1340)  m2scale,m3scale,m4scale,m5scale
 1340    format (//,15x,'Atomic Multipole Scaling Factors',
     &           ///,20x,'1-2 Atoms',f17.3,/,20x,'1-3 Atoms',f17.3,
     &           /,20x,'1-4 Atoms',f17.3,/,20x,'1-5 Atoms',f17.3)
      end if
c
c     atomic dipole polarizability parameters
c
      exist = .false.
      do i = 1, maxtyp
         if (polr(i) .ne. 0.0d0)  exist = .true.
      end do
      if (exist) then
         write (itxt,1350)  formfeed,forcefield
 1350    format (a1,//,15x,'TINKER Force Field Parameters for ',a20)
         write (itxt,1360)
 1360    format (//,15x,'Dipole Polarizability Parameters',
     &           ///,23x,'Type',7x,'Alpha',6x,'Damp',
     &              6x,'Group Atom Types',/)
         k = 0
         do i = 1, maxtyp
            if (polr(i) .ne. 0.0d0) then
               k = k + 1
               npg = 0
               do j = 1, maxval
                  if (pgrp(j,i) .ne. 0)  npg = npg + 1
               end do
               if (npg .eq. 0) then
                  write (itxt,1370)  k,i,polr(i),athl(i),adird(i)
 1370             format (10x,i5,7x,i4,3x,3f10.3)
               else
                  write (itxt,1380)  k,i,polr(i),athl(i),adird(i),
     &                               (pgrp(j,i),j=1,npg)
 1380             format (10x,i5,7x,i4,3x,3f10.3,4x,6i5)
               end if
            end if
         end do
c
c     dipole polarizability scaling parameters
c
         write (itxt,1390)  d1scale,d2scale,d3scale,d4scale
 1390    format (//,15x,'Direct Induction Scaling Factors',
     &           ///,20x,'1-1 Groups',f15.3,/,20x,'1-2 Groups',f15.3,
     &           /,20x,'1-3 Groups',f15.3,/,20x,'1-4 Groups',f15.3)
         write (itxt,1400)  u1scale,u2scale,u3scale,u4scale
 1400    format (//,15x,'Mutual Induction Scaling Factors',
     &           ///,20x,'1-1 Groups',f15.3,/,20x,'1-2 Groups',f15.3,
     &           /,20x,'1-3 Groups',f15.3,/,20x,'1-4 Groups',f15.3)
         write (itxt,1410)  p2scale,p3scale,p4scale,p5scale,p41scale
 1410    format (//,15x,'Polarizability Energy Scaling Factors',
     &           ///,20x,'1-2 Atoms',f16.3,/,20x,'1-3 Atoms',f16.3,
     &           /,20x,'1-4 Atoms',f16.3,/,20x,'1-5 Atoms',f16.3,
     &           /,20x,'1-4 Intra',f16.3)
      end if
c
c     conjugated pisystem atom parameters
c
      exist = .false.
      do i = 1, maxclass
         if (ionize(i) .ne. 0.0d0)  exist = .true.
      end do
      if (exist) then
         write (itxt,1420)  formfeed,forcefield
 1420    format (a1,//,15x,'TINKER Force Field Parameters for ',a20)
         write (itxt,1430)
 1430    format (//,15x,'Conjugated Pisystem Atom Parameters',
     &           ///,20x,'Class',3x,'Electron',
     &              3x,'Ionization',3x,'Repulsion',/)
         k = 0
         do i = 1, maxclass
            if (ionize(i) .ne. 0.0d0) then
               k = k + 1
               write (itxt,1440)  k,i,electron(i),ionize(i),repulse(i)
 1440          format (8x,i5,7x,i4,f10.1,2x,2f12.3)
            end if
         end do
      end if
c
c     conjugated pisystem bond parameters
c
      if (kpi(1) .ne. blank8) then
         write (itxt,1450)
 1450    format (//,15x,'Conjugated Pisystem Bond Parameters',
     &           ///,20x,'Classes',8x,'d Force',4x,'d Length',/)
         do i = 1, maxnpi
            if (kpi(i) .eq. blank8)  goto 1470
            k1 = number(kpi(i)(1:4))
            k2 = number(kpi(i)(5:8))
            write (itxt,1460)  i,k1,k2,sslope(i),tslope(i)
 1460       format (8x,i5,5x,i4,'-',i4,3x,f12.3,f12.3)
         end do
 1470    continue
      end if
c
c     conjugated pisystem bond parameters for 5-membered rings
c
      if (kpi5(1) .ne. blank8) then
         write (itxt,1480)
 1480    format (//,15x,'5-Membered Ring Pisystem Bond Parameters',
     &           ///,20x,'Classes',8x,'d Force',4x,'d Length',/)
         do i = 1, maxnpi5
            if (kpi5(i) .eq. blank8)  goto 1500
            k1 = number(kpi5(i)(1:4))
            k2 = number(kpi5(i)(5:8))
            write (itxt,1490)  i,k1,k2,sslope5(i),tslope5(i)
 1490       format (8x,i5,5x,i4,'-',i4,3x,f12.3,f12.3)
         end do
 1500    continue
      end if
c
c     conjugated pisystem bond parameters for 4-membered rings
c
      if (kpi4(1) .ne. blank8) then
         write (itxt,1510)
 1510    format (//,15x,'4-Membered Ring Pisystem Bond Parameters',
     &           ///,20x,'Classes',8x,'d Force',4x,'d Length',/)
         do i = 1, maxnpi4
            if (kpi4(i) .eq. blank8)  goto 1530
            k1 = number(kpi4(i)(1:4))
            k2 = number(kpi4(i)(5:8))
            write (itxt,1520)  i,k1,k2,sslope4(i),tslope4(i)
 1520       format (8x,i5,5x,i4,'-',i4,3x,f12.3,f12.3)
         end do
 1530    continue
      end if
      return
      end
