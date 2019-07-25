c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2004  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  program prmedit  --  edit and renumber parameter files  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "prmedit" reformats an existing parameter file, and revises
c     type and class numbers based on the "atom" parameter ordering
c
c
      program prmedit
      use iounit
      implicit none
      integer iprm,mode
      integer freeunit
      integer trimtext
      logical dotype,doclass
      logical exist,query
      character*240 prmfile
      character*240 string
c
c
c     read and store the original force field parameter file
c
      call initial
      call getprm
c
c     get the desired type of parameter file modification
c
      mode = 0
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  mode
         query = .false.
      end if
   10 continue
      if (query) then
         write (iout,20)
   20    format (/,' The Parameter Editing Facility can Provide :',
     &           //,4x,'(1) Format Individual Parameter Records',
     &           /,4x,'(2) Reorder Individual Parameter Records',
     &           /,4x,'(3) Renumber the Atom Types, and Reorder',
     &           /,4x,'(4) Renumber the Atom Classes, and Reorder',
     &           /,4x,'(5) Renumber Types and Classes, and Reorder',
     &           /,4x,'(6) Sort and Format Multipole Parameters',
     &           /,4x,'(7) Renumber and Format Biotype Parameters')
         do while (mode.lt.1 .or. mode.gt.7)
            mode = 0
            write (iout,30)
   30       format (/,' Enter the Number of the Desired Choice :  ',$)
            read (input,40,err=50,end=50)  mode
   40       format (i10)
   50       continue
         end do
      end if
c
c     set the renumbering operations to be performed
c
      dotype = .false.
      doclass = .false.
      if (mode .eq. 3)  dotype = .true.
      if (mode .eq. 4)  doclass = .true.
      if (mode .eq. 5) then
         dotype = .true.
         doclass = .true.
      end if
c
c     format records in the original parameter file
c
      if (mode .eq. 1) then
         iprm = freeunit ()
         prmfile = 'parameter.prm'
         call version (prmfile,'new')
         open (unit=iprm,file=prmfile,status='new')
         call prmform (iprm)
         write (iout,60)  prmfile(1:trimtext(prmfile))
   60    format (/,' Reformated Parameter File Written To:  ',a)
         close (unit=iprm)
      end if
c
c     reorder and renumber the original parameter file
c
      if (mode.ge.2 .and. mode.le.5) then
         iprm = freeunit ()
         prmfile = 'parameter.prm'
         call version (prmfile,'new')
         open (unit=iprm,file=prmfile,status='new')
         call prmorder (iprm,dotype,doclass)
         write (iout,70)  prmfile(1:trimtext(prmfile))
   70    format (/,' Renumbered Parameter File Written To:  ',a)
         close (unit=iprm)
      end if
c
c     sort the atomic multipole parameters by atom type
c
      if (mode .eq. 6) then
         iprm = freeunit ()
         prmfile = 'multipole.prm'
         call version (prmfile,'new')
         open (unit=iprm,file=prmfile,status='new')
         call polesort (iprm)
         write (iout,80)  prmfile(1:trimtext(prmfile))
   80    format (/,' Sorted Multipole Values Written To:  ',a)
         close (unit=iprm)
      end if
c
c     renumber and format any biotype parameter values
c
      if (mode .eq. 7) then
         iprm = freeunit ()
         prmfile = 'biotype.prm'
         call version (prmfile,'new')
         open (unit=iprm,file=prmfile,status='new')
         call biosort (iprm)
         write (iout,90)  prmfile(1:trimtext(prmfile))
   90    format (/,' Renumbered Biotype Values Written To:  ',a)
         close (unit=iprm)
      end if
      call final
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine prmform  --  reformat individual parameters  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "prmform" formats each individual parameter record to conform
c     to a consistent text layout
c
c
      subroutine prmform (iprm)
      use sizes
      use angpot
      use bndpot
      use math
      use params
      use urypot
      implicit none
      integer i,j,iprm
      integer ia,ib,ic
      integer id,ie
      integer length,next
      integer trimtext
      integer atn,lig
      integer kg,kt
      integer nx,ny,nxy
      integer ft(6)
      integer ig(20)
      real*8 wght
      real*8 rd,ep,rdn
      real*8 dl,fc,bd
      real*8 an1,an2,an3
      real*8 an,pr
      real*8 ba1,ba2
      real*8 ds,dk,vd,pt
      real*8 aa1,aa2,aa3
      real*8 bt1,bt2,bt3
      real*8 bt4,bt5,bt6
      real*8 bt7,bt8,bt9
      real*8 at1,at2,at3
      real*8 at4,at5,at6
      real*8 tx,ty,tf
      real*8 cg,dp,ps,pl
      real*8 pl1,pl2,pl3
      real*8 pol,thl
      real*8 el,iz,rp
      real*8 ss,ts
      real*8 vt(6),st(6)
      character*3 sym
      character*20 keyword
      character*24 note
      character*30 blank
      character*240 record
      character*240 string
c
c
c     reformat and print the various parameters
c
      i = 0
      blank = '                              '
      do while (i .lt. nprm)
         i = i + 1
         record = prmline(i)
         length = trimtext (record)
         next = 1
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:5) .eq. 'ATOM ') then
            ia = -1
            ib = -1
            sym = '   '
            note = '                        '
            atn = 0
            wght = 0.0d0
            lig = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            call gettext (record,sym,next)
            call getstring (record,note,next)
            string = record(next:240)
            read (string,*,err=10,end=10)  atn,wght,lig
   10       continue
            length = trimtext(note)
            string = '"'//note(1:length)//'"'//blank
            if (ib .ge. 0) then
               write (iprm,20)  ia,ib,sym,string(1:26),atn,wght,lig
   20          format ('atom',6x,2i5,4x,a3,3x,a26,1x,i5,f10.3,i5)
            else if (ia .ge. 0) then
               write (iprm,30)  ia,sym,string(1:26),atn,wght,lig
   30          format ('atom',6x,i5,4x,a3,3x,a26,1x,i5,f10.3,i5)
            else
               write (iprm,40)  record(1:length)
   40          format (a)
            end if
         else if (keyword(1:4) .eq. 'VDW ') then
            ia = 0
            rd = 0.0d0
            ep = 0.0d0
            rdn = 0.0d0
            read (string,*,err=50,end=50)  ia,rd,ep,rdn
   50       continue
            if (rdn .eq. 0.0d0) then
               write (iprm,60)  ia,rd,ep
   60          format ('vdw',7x,i5,10x,2f11.4)
            else
               write (iprm,70)  ia,rd,ep,rdn
   70          format ('vdw',7x,i5,10x,2f11.4,f11.3)
            end if
         else if (keyword(1:6) .eq. 'VDW14 ') then
            ia = 0
            rd = 0.0d0
            ep = 0.0d0
            read (string,*,err=80,end=80)  ia,rd,ep
   80       continue
            write (iprm,90)  ia,rd,ep
   90       format ('vdw14',5x,i5,10x,2f11.4)
         else if (keyword(1:6) .eq. 'VDWPR ') then
            ia = 0
            ib = 0
            rd = 0.0d0
            ep = 0.0d0
            read (string,*,err=100,end=100)  ia,ib,rd,ep
  100       continue
            write (iprm,110)  ia,ib,rd,ep
  110       format ('vdwpr',5x,2i5,5x,2f11.4)
         else if (keyword(1:6) .eq. 'HBOND ') then
            ia = 0
            ib = 0
            rd = 0.0d0
            ep = 0.0d0
            read (string,*,err=120,end=120)  ia,ib,rd,ep
  120       continue
            write (iprm,130)  ia,ib,rd,ep
  130       format ('hbond',5x,2i5,5x,2f11.4)
         else if (keyword(1:5) .eq. 'BOND ') then
            ia = 0
            ib = 0
            fc = 0.0d0
            bd = 0.0d0
            read (string,*,err=140,end=140)  ia,ib,fc,bd
  140       continue
            if (bndunit .lt. 10.0d0) then
               write (iprm,150)  ia,ib,fc,bd
  150          format ('bond',6x,2i5,5x,f11.2,f11.4)
            else
               write (iprm,160)  ia,ib,fc,bd
  160          format ('bond',6x,2i5,5x,f11.3,f11.4)
            end if
         else if (keyword(1:6) .eq. 'BOND5 ') then
            ia = 0
            ib = 0
            fc = 0.0d0
            bd = 0.0d0
            read (string,*,err=170,end=170)  ia,ib,fc,bd
  170       continue
            if (bndunit .lt. 10.0d0) then
               write (iprm,180)  ia,ib,fc,bd
  180          format ('bond5',5x,2i5,5x,f11.2,f11.4)
            else
               write (iprm,190)  ia,ib,fc,bd
  190          format ('bond5',5x,2i5,5x,f11.3,f11.4)
            end if
         else if (keyword(1:6) .eq. 'BOND4 ') then
            ia = 0
            ib = 0
            fc = 0.0d0
            bd = 0.0d0
            read (string,*,err=200,end=200)  ia,ib,fc,bd
  200       continue
            if (bndunit .lt. 10.0d0) then
               write (iprm,210)  ia,ib,fc,bd
  210          format ('bond4',5x,2i5,5x,f11.2,f11.4)
            else
               write (iprm,220)  ia,ib,fc,bd
  220          format ('bond4',5x,2i5,5x,f11.3,f11.4)
            end if
         else if (keyword(1:6) .eq. 'BOND3 ') then
            ia = 0
            ib = 0
            fc = 0.0d0
            bd = 0.0d0
            read (string,*,err=230,end=230)  ia,ib,fc,bd
  230       continue
            if (bndunit .lt. 10.0d0) then
               write (iprm,240)  ia,ib,fc,bd
  240          format ('bond3',5x,2i5,5x,f11.2,f11.4)
            else
               write (iprm,250)  ia,ib,fc,bd
  250          format ('bond3',5x,2i5,5x,f11.3,f11.4)
            end if
         else if (keyword(1:9) .eq. 'ELECTNEG ') then
            ia = 0
            ib = 0
            ic = 0
            dl = 0.0d0
            read (string,*,err=260,end=260)  ia,ib,ic,dl
  260       continue
            write (iprm,270)  ia,ib,ic,dl
  270       format ('electneg',2x,3i5,11x,f11.4)
         else if (keyword(1:6) .eq. 'ANGLE ') then
            ia = 0
            ib = 0
            ic = 0
            fc = 0.0d0
            an1 = 0.0d0
            an2 = 0.0d0
            an3 = 0.0d0
            read (string,*,err=280,end=280)  ia,ib,ic,fc,an1,an2,an3
  280       continue
            if (an2.eq.0.0d0 .and. an3.eq.0.0d0) then
               if (angunit .lt. 10.0d0/radian**2) then
                  write (iprm,290)  ia,ib,ic,fc,an1
  290             format ('angle',5x,3i5,f11.2,f11.2)
               else
                  write (iprm,300)  ia,ib,ic,fc,an1
  300             format ('angle',5x,3i5,f11.3,f11.2)
               end if
            else
               if (angunit .lt. 10.0d0/radian**2) then
                  write (iprm,310)  ia,ib,ic,fc,an1,an2,an3
  310             format ('angle',5x,3i5,f11.2,3f11.2)
               else
                  write (iprm,320)  ia,ib,ic,fc,an1,an2,an3
  320             format ('angle',5x,3i5,f11.3,3f11.2)
               end if
            end if
         else if (keyword(1:7) .eq. 'ANGLE5 ') then
            ia = 0
            ib = 0
            ic = 0
            fc = 0.0d0
            an1 = 0.0d0
            an2 = 0.0d0
            an3 = 0.0d0
            read (string,*,err=330,end=330)  ia,ib,ic,fc,an1,an2,an3
  330       continue
            if (an2.eq.0.0d0 .and. an3.eq.0.0d0) then
               if (angunit .lt. 10.0d0/radian**2) then
                  write (iprm,340)  ia,ib,ic,fc,an1
  340             format ('angle5',4x,3i5,f11.2,f11.2)
               else
                  write (iprm,350)  ia,ib,ic,fc,an1
  350             format ('angle5',4x,3i5,f11.3,f11.2)
               end if
            else
               if (angunit .lt. 10.0d0/radian**2) then
                  write (iprm,360)  ia,ib,ic,fc,an1,an2,an3
  360             format ('angle5',4x,3i5,f11.2,3f11.2)
               else
                  write (iprm,370)  ia,ib,ic,fc,an1,an2,an3
  370             format ('angle5',4x,3i5,f11.3,3f11.2)
               end if
            end if
         else if (keyword(1:7) .eq. 'ANGLE4 ') then
            ia = 0
            ib = 0
            ic = 0
            fc = 0.0d0
            an1 = 0.0d0
            an2 = 0.0d0
            an3 = 0.0d0
            read (string,*,err=380,end=380)  ia,ib,ic,fc,an1,an2,an3
  380       continue
            if (an2.eq.0.0d0 .and. an3.eq.0.0d0) then
               if (angunit .lt. 10.0d0/radian**2) then
                  write (iprm,390)  ia,ib,ic,fc,an1
  390             format ('angle4',4x,3i5,f11.2,f11.2)
               else
                  write (iprm,400)  ia,ib,ic,fc,an1
  400             format ('angle4',4x,3i5,f11.3,f11.2)
               end if
            else
               if (angunit .lt. 10.0d0/radian**2) then
                  write (iprm,410)  ia,ib,ic,fc,an1,an2,an3
  410             format ('angle4',4x,3i5,f11.2,3f11.2)
               else
                  write (iprm,420)  ia,ib,ic,fc,an1,an2,an3
  420             format ('angle4',4x,3i5,f11.3,3f11.2)
               end if
            end if
         else if (keyword(1:7) .eq. 'ANGLE3 ') then
            ia = 0
            ib = 0
            ic = 0
            fc = 0.0d0
            an1 = 0.0d0
            an2 = 0.0d0
            an3 = 0.0d0
            read (string,*,err=430,end=430)  ia,ib,ic,fc,an1,an2,an3
  430       continue
            if (an2.eq.0.0d0 .and. an3.eq.0.0d0) then
               if (angunit .lt. 10.0d0/radian**2) then
                  write (iprm,440)  ia,ib,ic,fc,an1
  440             format ('angle3',4x,3i5,f11.2,f11.2)
               else
                  write (iprm,450)  ia,ib,ic,fc,an1
  450             format ('angle3',4x,3i5,f11.3,f11.2)
               end if
            else
               if (angunit .lt. 10.0d0/radian**2) then
                  write (iprm,460)  ia,ib,ic,fc,an1,an2,an3
  460             format ('angle3',4x,3i5,f11.2,3f11.2)
               else
                  write (iprm,470)  ia,ib,ic,fc,an1,an2,an3
  470             format ('angle3',4x,3i5,f11.3,3f11.2)
               end if
            end if
         else if (keyword(1:7) .eq. 'ANGLEF ') then
            ia = 0
            ib = 0
            ic = 0
            fc = 0.0d0
            an = 0.0d0
            pr = 0.0d0
            read (string,*,err=480,end=480)  ia,ib,ic,fc,an,pr
  480       continue
            if (angunit .lt. 10.0d0/radian**2) then
               write (iprm,490)  ia,ib,ic,fc,an,pr
  490          format ('anglef',4x,3i5,f11.2,f11.2,f11.1)
            else
               write (iprm,500)  ia,ib,ic,fc,an,pr
  500          format ('anglef',4x,3i5,f11.3,f11.2,f11.1)
            end if
         else if (keyword(1:7) .eq. 'STRBND ') then
            ia = 0
            ib = 0
            ic = 0
            ba1 = 0.0d0
            ba2 = 0.0d0
            read (string,*,err=510,end=510)  ia,ib,ic,ba1,ba2
  510       continue
            if (stbnunit .lt. 10.0d0/radian) then
               write (iprm,520)  ia,ib,ic,ba1,ba2
  520          format ('strbnd',4x,3i5,2f11.2)
            else
               write (iprm,530)  ia,ib,ic,ba1,ba2
  530          format ('strbnd',4x,3i5,2f11.3)
            end if
         else if (keyword(1:9) .eq. 'UREYBRAD ') then
            ia = 0
            ib = 0
            ic = 0
            fc = 0.0d0
            ds = 0.0d0
            read (string,*,err=540,end=540)  ia,ib,ic,fc,ds
  540       continue
            if (ureyunit .lt. 10.0d0) then
               write (iprm,550)  ia,ib,ic,fc,ds
  550          format ('ureybrad',2x,3i5,f11.2,f11.4)
            else
               write (iprm,560)  ia,ib,ic,fc,ds
  560          format ('ureybrad',2x,3i5,f11.3,f11.4)
            end if
         else if (keyword(1:7) .eq. 'ANGANG ') then
            ia = 0
            aa1 = 0.0d0
            aa2 = 0.0d0
            aa3 = 0.0d0
            read (string,*,err=570,end=570)  ia,aa1,aa2,aa3
  570       continue
            if (abs(aaunit) .lt. 10.0d0/radian**2) then
               write (iprm,580)  ia,aa1,aa2,aa3
  580          format ('angang',4x,i5,10x,3f11.2)
            else
               write (iprm,590)  ia,aa1,aa2,aa3
  590          format ('angang',4x,i5,10x,3f11.3)
            end if
         else if (keyword(1:7) .eq. 'OPBEND ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            fc = 0.0d0
            read (string,*,err=600,end=600)  ia,ib,ic,id,fc
  600       continue
            if (opbunit .lt. 10.0d0/radian**2) then
               write (iprm,610)  ia,ib,ic,id,fc
  610          format ('opbend',4x,4i5,6x,f11.2)
            else
               write (iprm,620)  ia,ib,ic,id,fc
  620          format ('opbend',4x,4i5,6x,f11.3)
            end if
         else if (keyword(1:7) .eq. 'OPDIST ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            fc = 0.0d0
            read (string,*,err=630,end=630)  ia,ib,ic,id,fc
  630       continue
            if (opdunit .lt. 10.0d0) then
               write (iprm,640)  ia,ib,ic,id,fc
  640          format ('opdist',4x,4i5,6x,f11.2)
            else
               write (iprm,650)  ia,ib,ic,id,fc
  650          format ('opdist',4x,4i5,6x,f11.3)
            end if
         else if (keyword(1:9) .eq. 'IMPROPER ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            dk = 0.0d0
            vd = 0.0d0
            read (string,*,err=660,end=660)  ia,ib,ic,id,dk,vd
  660       continue
            write (iprm,670)  ia,ib,ic,id,dk,vd
  670       format ('improper',2x,4i5,6x,2f11.2)
         else if (keyword(1:8) .eq. 'IMPTORS ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            do j = 1, 6
               vt(j) = 0.0d0
               st(j) = 0.0d0
               ft(j) = 0
            end do
            read (string,*,err=680,end=680)  ia,ib,ic,id,
     &                                       (vt(j),st(j),ft(j),j=1,6)
  680       continue
            kt = 0
            do j = 1, 6
               if (ft(j) .ne. 0) then
                  kt = j
               end if
            end do
            write (iprm,690)  ia,ib,ic,id,(vt(j),st(j),ft(j),j=1,kt)
  690       format ('imptors',3x,4i5,6x,6(f11.3,f7.1,i3))
         else if (keyword(1:8) .eq. 'TORSION ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            do j = 1, 6
               vt(j) = 0.0d0
               st(j) = 0.0d0
               ft(j) = 0
            end do
            read (string,*,err=700,end=700)  ia,ib,ic,id,
     &                                       (vt(j),st(j),ft(j),j=1,6)
  700       continue
            kt = 0
            do j = 1, 6
               if (ft(j) .ne. 0) then
                  kt = j
               end if
            end do
            if (kt.eq.3 .and. st(1).eq.0.0d0 .and. st(2).eq.180.0d0
     &                  .and. st(3).eq.0.0d0) then
               write (iprm,710)  ia,ib,ic,id,(vt(j),st(j),ft(j),j=1,kt)
  710          format ('torsion',3x,4i5,3x,f8.3,f4.1,i2,
     &                    f8.3,f6.1,i2,f8.3,f4.1,i2)
            else if (kt .le. 2) then
               write (iprm,720)  ia,ib,ic,id,(vt(j),st(j),ft(j),j=1,kt)
  720          format ('torsion',3x,4i5,6x,2(f11.3,f7.1,i3))
            else
               write (iprm,730)  ia,ib,ic,id,(vt(j),st(j),ft(j),j=1,kt)
  730          format ('torsion',3x,4i5,3x,6(f8.3,f6.1,i2))
            end if
         else if (keyword(1:9) .eq. 'TORSION5 ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            do j = 1, 6
               vt(j) = 0.0d0
               st(j) = 0.0d0
               ft(j) = 0
            end do
            read (string,*,err=740,end=740)  ia,ib,ic,id,
     &                                       (vt(j),st(j),ft(j),j=1,6)
  740       continue
            kt = 0
            do j = 1, 6
               if (ft(j) .ne. 0) then
                  kt = j
               end if
            end do
            if (kt.eq.3 .and. st(1).eq.0.0d0 .and. st(2).eq.180.0d0
     &                  .and. st(3).eq.0.0d0) then
               write (iprm,750)  ia,ib,ic,id,(vt(j),st(j),ft(j),j=1,kt)
  750          format ('torsion5',2x,4i5,3x,f8.3,f4.1,i2,
     &                    f8.3,f6.1,i2,f8.3,f4.1,i2)
            else if (kt .le. 2) then
               write (iprm,760)  ia,ib,ic,id,(vt(j),st(j),ft(j),j=1,kt)
  760          format ('torsion5',2x,4i5,6x,2(f11.3,f7.1,i3))
            else
               write (iprm,770)  ia,ib,ic,id,(vt(j),st(j),ft(j),j=1,kt)
  770          format ('torsion5',2x,4i5,3x,6(f8.3,f6.1,i2))
            end if
         else if (keyword(1:9) .eq. 'TORSION4 ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            do j = 1, 6
               vt(j) = 0.0d0
               st(j) = 0.0d0
               ft(j) = 0
            end do
            read (string,*,err=780,end=780)  ia,ib,ic,id,
     &                                       (vt(j),st(j),ft(j),j=1,6)
  780       continue
            kt = 0
            do j = 1, 6
               if (ft(j) .ne. 0) then
                  kt = j
               end if
            end do
            if (kt.eq.3 .and. st(1).eq.0.0d0 .and. st(2).eq.180.0d0
     &                  .and. st(3).eq.0.0d0) then
               write (iprm,790)  ia,ib,ic,id,(vt(j),st(j),ft(j),j=1,kt)
  790          format ('torsion4',2x,4i5,3x,f8.3,f4.1,i2,
     &                    f8.3,f6.1,i2,f8.3,f4.1,i2)
            else if (kt .le. 2) then
               write (iprm,800)  ia,ib,ic,id,(vt(j),st(j),ft(j),j=1,kt)
  800          format ('torsion4',2x,4i5,6x,2(f11.3,f7.1,i3))
            else
               write (iprm,810)  ia,ib,ic,id,(vt(j),st(j),ft(j),j=1,kt)
  810          format ('torsion4',2x,4i5,3x,6(f8.3,f6.1,i2))
            end if
         else if (keyword(1:7) .eq. 'PITORS ') then
            ia = 0
            ib = 0
            pt = 0.0d0
            read (string,*,err=820,end=820)  ia,ib,pt
  820       continue
            write (iprm,830)  ia,ib,pt
  830       format ('pitors',4x,2i5,5x,f11.2)
         else if (keyword(1:8) .eq. 'STRTORS ') then
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
            read (string,*,err=840,end=840)  ia,ib,ic,id,bt1,bt2,bt3,
     &                                       bt4,bt5,bt6,bt7,bt8,bt9
  840       continue
            write (iprm,850)  ia,ib,ic,id,bt1,bt2,bt3,
     &                        bt4,bt5,bt6,bt7,bt8,bt9
  850       format ('strtors',3x,4i5,1x,9f8.3)
         else if (keyword(1:8) .eq. 'ANGTORS ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            at1 = 0.0d0
            at2 = 0.0d0
            at3 = 0.0d0
            at4 = 0.0d0
            at5 = 0.0d0
            at6 = 0.0d0
            read (string,*,err=860,end=860)  ia,ib,ic,id,at1,at2,
     &                                       at3,at4,at5,at6
  860       continue
            write (iprm,870)  ia,ib,ic,id,at1,at2,at3,at4,at5,at6
  870       format ('angtors',3x,4i5,1x,6f8.3)
         else if (keyword(1:8) .eq. 'TORTORS ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            ie = 0
            nx = 0
            ny = 0
            read (string,*,err=880,end=880)  ia,ib,ic,id,ie,nx,ny
  880       continue
            write (iprm,890)  ia,ib,ic,id,ie,nx,ny
  890       format ('tortors',3x,5i5,5x,2i5)
            nxy = nx * ny
            do j = 1, nxy
               i = i + 1
               record = prmline(i)
               read (record,*,err=900,end=900)  tx,ty,tf
  900          continue
               write (iprm,910)  tx,ty,tf
  910          format (f8.1,1x,f8.1,1x,f11.5)
            end do
         else if (keyword(1:7) .eq. 'CHARGE ') then
            ia = 0
            cg = 0.0d0
            read (string,*,err=920,end=920)  ia,cg
  920       continue
            write (iprm,930)  ia,cg
  930       format ('charge',4x,i5,10x,f11.4)
         else if (keyword(1:7) .eq. 'DIPOLE ') then
            ia = 0
            ib = 0
            dp = 0.0d0
            ps = 0.5d0
            read (string,*,err=940,end=940)  ia,ib,dp,ps
  940       continue
            write (iprm,950)  ia,ib,dp,ps
  950       format ('dipole',4x,2i5,5x,f11.4,f11.3)
         else if (keyword(1:8) .eq. 'DIPOLE5 ') then
            ia = 0
            ib = 0
            dp = 0.0d0
            ps = 0.5d0
            read (string,*,err=960,end=960)  ia,ib,dp,ps
  960       continue
            write (iprm,970)  ia,ib,dp,ps
  970       format ('dipole5',3x,2i5,5x,f11.4,f11.3)
         else if (keyword(1:8) .eq. 'DIPOLE4 ') then
            ia = 0
            ib = 0
            dp = 0.0d0
            ps = 0.5d0
            read (string,*,err=980,end=980)  ia,ib,dp,ps
  980       continue
            write (iprm,990)  ia,ib,dp,ps
  990       format ('dipole4',3x,2i5,5x,f11.4,f11.3)
         else if (keyword(1:8) .eq. 'DIPOLE3 ') then
            ia = 0
            ib = 0
            dp = 0.0d0
            ps = 0.5d0
            read (string,*,err=1000,end=1000)  ia,ib,dp,ps
 1000       continue
            write (iprm,1010)  ia,ib,dp,ps
 1010       format ('dipole3',3x,2i5,5x,f11.4,f11.3)
         else if (keyword(1:10) .eq. 'MULTIPOLE ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            pl = 0.0d0
            read (string,*,err=1020,end=1020)  ia,ib,ic,id,pl
            goto 1030
 1020       continue
            id = 0
            read (string,*,err=1030,end=1030)  ia,ib,ic,pl
 1030       continue
            if (id .eq. 0) then
               write (iprm,1040)  ia,ib,ic,pl
 1040          format ('multipole',1x,3i5,11x,f11.5)
            else
               write (iprm,1050)  ia,ib,ic,id,pl
 1050          format ('multipole',1x,4i5,6x,f11.5)
            end if
            i = i + 1
            record = prmline(i)
            read (record,*,err=1060,end=1060)  pl1,pl2,pl3
 1060       continue
            write (iprm,1070)  pl1,pl2,pl3
 1070       format (36x,3f11.5)
            i = i + 1
            record = prmline(i)
            read (record,*,err=1080,end=1080)  pl1
 1080       continue
            write (iprm,1090)  pl1
 1090       format (36x,f11.5)
            i = i + 1
            record = prmline(i)
            read (record,*,err=1100,end=1100)  pl1,pl2
 1100       continue
            write (iprm,1110)  pl1,pl2
 1110       format (36x,2f11.5)
            i = i + 1
            record = prmline(i)
            read (record,*,err=1120,end=1120)  pl1,pl2,pl3
 1120       continue
            write (iprm,1130)  pl1,pl2,pl3
 1130       format (36x,3f11.5)
         else if (keyword(1:9) .eq. 'POLARIZE ') then
            ia = 0
            pol = 0.0d0
            thl = 0.0d0
            do j = 1, 20
               ig(j) = 0
            end do
            read (string,*,err=1140,end=1140)  ia,pol,thl,
     &                                         (ig(j),j=1,20)
 1140       continue
            kg = 0
            do j = 1, maxval
               if (ig(j) .ne. 0) then
                  kg = j
               end if
            end do
            call sort (kg,ig)
            write (iprm,1150)  ia,pol,thl,(ig(j),j=1,kg)
 1150       format ('polarize',2x,i5,5x,2f11.4,2x,20i5)
         else if (keyword(1:7) .eq. 'PIATOM ') then
            ia = 0
            el = 0.0d0
            iz = 0.0d0
            rp = 0.0d0
            read (string,*,err=1160,end=1160)  ia,el,iz,rp
 1160       continue
            write (iprm,1170)  ia,el,iz,rp
 1170       format ('piatom',4x,i5,10x,f11.1,2f11.3)
         else if (keyword(1:7) .eq. 'PIBOND ') then
            ia = 0
            ib = 0
            ss = 0.0d0
            ts = 0.0d0
            read (string,*,err=1180,end=1180)  ia,ib,ss,ts
 1180       continue
            write (iprm,1190)  ia,ib,ss,ts
 1190       format ('pibond',4x,2i5,5x,f11.3,f11.4)
         else if (keyword(1:8) .eq. 'PIBOND5 ') then
            ia = 0
            ib = 0
            ss = 0.0d0
            ts = 0.0d0
            read (string,*,err=1200,end=1200)  ia,ib,ss,ts
 1200       continue
            write (iprm,1210)  ia,ib,ss,ts
 1210       format ('pibond5',3x,2i5,5x,f11.3,f11.4)
         else if (keyword(1:8) .eq. 'PIBOND4 ') then
            ia = 0
            ib = 0
            ss = 0.0d0
            ts = 0.0d0
            read (string,*,err=1220,end=1220)  ia,ib,ss,ts
 1220       continue
            write (iprm,1230)  ia,ib,ss,ts
 1230       format ('pibond4',3x,2i5,5x,f11.3,f11.4)
         else if (keyword(1:6) .eq. 'METAL ') then
            ia = 0
            call getnumb (record,ia,next)
            write (iprm,1240)  ia,record(next:length)
 1240       format ('metal',5x,i5,a)
         else if (keyword(1:8) .eq. 'BIOTYPE ') then
            ia = 0
            ib = 0
            sym = '   '
            note = '                        '
            read (string,*,err=1250,end=1250)  ia
            call getword (record,sym,next)
            call getstring (record,note,next)
            string = record(next:240)
            read (string,*,err=1250,end=1250)  ib
 1250       continue
            length = trimtext(note)
            string = '"'//note(1:length)//'"'//blank
            write (iprm,1260)  ia,sym,string(1:30),ib
 1260       format ('biotype',3x,i5,4x,a3,5x,a30,2x,i5)
         else if (length .eq. 0) then
            write (iprm,1270)
 1270       format ()
         else
            write (iprm,1280)  record(1:length)
 1280       format (a)
         end if
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine prmorder  --  reorder atom types and classes  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "prmorder" places a list of atom type or class numbers into
c     canonical order for potential energy parameter definitions
c
c
      subroutine prmorder (iprm,dotype,doclass)
      use sizes
      use iounit
      use params
      use vdwpot
      implicit none
      integer i,j,iprm
      integer it,ic,kt,kc
      integer ia,ib,id,ie
      integer offset,next
      integer length
      integer trimtext
      integer kg,ig(20)
      integer itype(0:maxtyp)
      integer iclass(0:maxclass)
      real*8 pol,thl
      logical dotype,doclass
      logical prtclass
      character*20 keyword
      character*30 blank
      character*240 record
      character*240 string
c
c
c     zero out the storage for atom types and classes
c
      ia = 0
      ib = 0
      ic = 0
      id = 0
      ie = 0
      kt = 0
      kc = 0
      do i = 0, maxtyp
         itype(i) = 0
      end do
      do i = 0, maxclass
         iclass(i) = 0
      end do
      blank = '                              '
c
c     get the starting numbers for atom types and classes
c
      if (dotype) then
         write (iout,10)
   10    format (/,' Enter Starting Number for Atom Types [1] :  ',$)
         read (input,20)  offset
   20    format (i10)
         if (offset .gt. 0)  kt = offset - 1
      end if
      if (doclass) then
         write (iout,30)
   30    format (/,' Enter Starting Number for Atom Classes [1] :  ',$)
         read (input,40)  offset
   40    format (i10)
         if (offset .gt. 0)  kc = offset - 1
      end if
c
c     count, order and test equivalence of atom types and classes
c
      prtclass = .false.
      do i = 1, nprm
         record = prmline(i)
         next = 1
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:5) .eq. 'ATOM ') then
            it = 0
            ic = 0
            call getnumb (record,it,next)
            call getnumb (record,ic,next)
            if (ic .eq. 0)  ic = it
            if (it .ne. ic)  prtclass = .true.
            if (itype(it) .eq. 0) then
               kt = kt + 1
               if (dotype) then
                  itype(it) = kt
               else
                  itype(it) = it
               end if
            end if
            if (iclass(ic) .eq. 0) then
               kc = kc + 1
               if (doclass) then
                  iclass(ic) = kc
               else
                  iclass(ic) = ic
               end if
            end if
         end if
      end do
c
c     reorder, renumber and print the various parameters
c
      do i = 1, nprm
         record = prmline(i)
         length = trimtext (record)
         next = 1
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:5) .eq. 'ATOM ') then
            ia = 0
            ib = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            if (ib .eq. 0)  ib = ia
            ia = itype(ia)
            ib = iclass(ib)
            if (prtclass) then
               write (iprm,50)  ia,ib,record(next:length)
   50          format ('atom',6x,2i5,a)
            else
               write (iprm,60)  ia,record(next:length)
   60          format ('atom',6x,i5,a)
            end if
         else if (keyword(1:4) .eq. 'VDW ') then
            ia = 0
            call getnumb (record,ia,next)
            if (vdwindex .eq. 'TYPE') then
               ia = itype(ia)
            else
               ia = iclass(ia)
            end if
            write (iprm,70)  ia,record(next:length)
   70       format ('vdw',7x,i5,a)
         else if (keyword(1:6) .eq. 'VDW14 ') then
            ia = 0
            call getnumb (record,ia,next)
            if (vdwindex .eq. 'TYPE') then
               ia = itype(ia)
            else
               ia = iclass(ia)
            end if
            write (iprm,80)  ia,record(next:length)
   80       format ('vdw14',5x,i5,a)
         else if (keyword(1:6) .eq. 'VDWPR ') then
            ia = 0
            ib = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            if (vdwindex .eq. 'TYPE') then
               ia = itype(ia)
               ib = itype(ib)
            else
               ia = iclass(ia)
               ib = iclass(ib)
            end if
            call prmsort (2,ia,ib,0,0,0)
            write (iprm,90)  ia,ib,record(next:length)
   90       format ('vdwpr',5x,2i5,a)
         else if (keyword(1:6) .eq. 'HBOND ') then
            ia = 0
            ib = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            if (vdwindex .eq. 'TYPE') then
               ia = itype(ia)
               ib = itype(ib)
            else
               ia = iclass(ia)
               ib = iclass(ib)
            end if
            call prmsort (2,ia,ib,0,0,0)
            write (iprm,100)  ia,ib,record(next:length)
  100       format ('hbond',5x,2i5,a)
         else if (keyword(1:5) .eq. 'BOND ') then
            ia = 0
            ib = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            ia = iclass(ia)
            ib = iclass(ib)
            call prmsort (2,ia,ib,0,0,0)
            write (iprm,110)  ia,ib,record(next:length)
  110       format ('bond',6x,2i5,a)
         else if (keyword(1:6) .eq. 'BOND5 ') then
            ia = 0
            ib = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            ia = iclass(ia)
            ib = iclass(ib)
            call prmsort (2,ia,ib,0,0,0)
            write (iprm,120)  ia,ib,record(next:length)
  120       format ('bond5',5x,2i5,a)
         else if (keyword(1:6) .eq. 'BOND4 ') then
            ia = 0
            ib = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            ia = iclass(ia)
            ib = iclass(ib)
            call prmsort (2,ia,ib,0,0,0)
            write (iprm,130)  ia,ib,record(next:length)
  130       format ('bond4',5x,2i5,a)
         else if (keyword(1:6) .eq. 'BOND3 ') then
            ia = 0
            ib = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            ia = iclass(ia)
            ib = iclass(ib)
            call prmsort (2,ia,ib,0,0,0)
            write (iprm,140)  ia,ib,record(next:length)
  140       format ('bond3',5x,2i5,a)
         else if (keyword(1:9) .eq. 'ELECTNEG ') then
            ia = 0
            ib = 0
            ic = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            call getnumb (record,ic,next)
            ia = iclass(ia)
            ib = iclass(ib)
            ic = iclass(ic)
            write (iprm,150)  ia,ib,ic,record(next:length)
  150       format ('electneg',2x,3i5,a)
         else if (keyword(1:6) .eq. 'ANGLE ') then
            ia = 0
            ib = 0
            ic = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            call getnumb (record,ic,next)
            ia = iclass(ia)
            ib = iclass(ib)
            ic = iclass(ic)
            call prmsort (3,ia,ib,ic,0,0)
            write (iprm,160)  ia,ib,ic,record(next:length)
  160       format ('angle',5x,3i5,a)
         else if (keyword(1:7) .eq. 'ANGLE5 ') then
            ia = 0
            ib = 0
            ic = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            call getnumb (record,ic,next)
            ia = iclass(ia)
            ib = iclass(ib)
            ic = iclass(ic)
            call prmsort (3,ia,ib,ic,0,0)
            write (iprm,170)  ia,ib,ic,record(next:length)
  170       format ('angle5',4x,3i5,a)
         else if (keyword(1:7) .eq. 'ANGLE4 ') then
            ia = 0
            ib = 0
            ic = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            call getnumb (record,ic,next)
            ia = iclass(ia)
            ib = iclass(ib)
            ic = iclass(ic)
            call prmsort (3,ia,ib,ic,0,0)
            write (iprm,180)  ia,ib,ic,record(next:length)
  180       format ('angle4',4x,3i5,a)
         else if (keyword(1:7) .eq. 'ANGLE3 ') then
            ia = 0
            ib = 0
            ic = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            call getnumb (record,ic,next)
            ia = iclass(ia)
            ib = iclass(ib)
            ic = iclass(ic)
            call prmsort (3,ia,ib,ic,0,0)
            write (iprm,190)  ia,ib,ic,record(next:length)
  190       format ('angle3',4x,3i5,a)
         else if (keyword(1:7) .eq. 'ANGLEF ') then
            ia = 0
            ib = 0
            ic = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            call getnumb (record,ic,next)
            ia = iclass(ia)
            ib = iclass(ib)
            ic = iclass(ic)
            call prmsort (3,ia,ib,ic,0,0)
            write (iprm,200)  ia,ib,ic,record(next:length)
  200       format ('anglef',4x,3i5,a)
         else if (keyword(1:7) .eq. 'STRBND ') then
            ia = 0
            ib = 0
            ic = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            call getnumb (record,ic,next)
            ia = iclass(ia)
            ib = iclass(ib)
            ic = iclass(ic)
            write (iprm,210)  ia,ib,ic,record(next:length)
  210       format ('strbnd',4x,3i5,a)
         else if (keyword(1:9) .eq. 'UREYBRAD ') then
            ia = 0
            ib = 0
            ic = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            call getnumb (record,ic,next)
            ia = iclass(ia)
            ib = iclass(ib)
            ic = iclass(ic)
            call prmsort (3,ia,ib,ic,0,0)
            write (iprm,220)  ia,ib,ic,record(next:length)
  220       format ('ureybrad',2x,3i5,a)
         else if (keyword(1:7) .eq. 'ANGANG ') then
            ia = 0
            call getnumb (record,ia,next)
            ia = iclass(ia)
            write (iprm,230)  ia,record(next:length)
  230       format ('angang',4x,i5,a)
         else if (keyword(1:7) .eq. 'OPBEND ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            call getnumb (record,ic,next)
            call getnumb (record,id,next)
            ia = iclass(ia)
            ib = iclass(ib)
            ic = iclass(ic)
            id = iclass(id)
            call prmsort (2,ic,id,0,0,0)
            write (iprm,240)  ia,ib,ic,id,record(next:length)
  240       format ('opbend',4x,4i5,a)
         else if (keyword(1:7) .eq. 'OPDIST ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            call getnumb (record,ic,next)
            call getnumb (record,id,next)
            ia = iclass(ia)
            ib = iclass(ib)
            ic = iclass(ic)
            id = iclass(id)
            call prmsort (2,ib,ic,0,0,0)
            call prmsort (2,ib,id,0,0,0)
            call prmsort (2,ic,id,0,0,0)
            write (iprm,250)  ia,ib,ic,id,record(next:length)
  250       format ('opdist',4x,4i5,a)
         else if (keyword(1:9) .eq. 'IMPROPER ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            call getnumb (record,ic,next)
            call getnumb (record,id,next)
            ia = iclass(ia)
            ib = iclass(ib)
            ic = iclass(ic)
            id = iclass(id)
            write (iprm,260)  ia,ib,ic,id,record(next:length)
  260       format ('improper',2x,4i5,a)
         else if (keyword(1:8) .eq. 'IMPTORS ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            call getnumb (record,ic,next)
            call getnumb (record,id,next)
            ia = iclass(ia)
            ib = iclass(ib)
            ic = iclass(ic)
            id = iclass(id)
            write (iprm,270)  ia,ib,ic,id,record(next:length)
  270       format ('imptors',3x,4i5,a)
         else if (keyword(1:8) .eq. 'TORSION ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            call getnumb (record,ic,next)
            call getnumb (record,id,next)
            ia = iclass(ia)
            ib = iclass(ib)
            ic = iclass(ic)
            id = iclass(id)
            call prmsort (4,ia,ib,ic,id,0)
            write (iprm,280)  ia,ib,ic,id,record(next:length)
  280       format ('torsion',3x,4i5,a)
         else if (keyword(1:9) .eq. 'TORSION5 ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            call getnumb (record,ic,next)
            call getnumb (record,id,next)
            ia = iclass(ia)
            ib = iclass(ib)
            ic = iclass(ic)
            id = iclass(id)
            call prmsort (4,ia,ib,ic,id,0)
            write (iprm,290)  ia,ib,ic,id,record(next:length)
  290       format ('torsion5',2x,4i5,a)
         else if (keyword(1:9) .eq. 'TORSION4 ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            call getnumb (record,ic,next)
            call getnumb (record,id,next)
            ia = iclass(ia)
            ib = iclass(ib)
            ic = iclass(ic)
            id = iclass(id)
            call prmsort (4,ia,ib,ic,id,0)
            write (iprm,300)  ia,ib,ic,id,record(next:length)
  300       format ('torsion4',2x,4i5,a)
         else if (keyword(1:7) .eq. 'PITORS ') then
            ia = 0
            ib = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            ia = iclass(ia)
            ib = iclass(ib)
            call prmsort (2,ia,ib,0,0,0)
            write (iprm,310)  ia,ib,record(next:length)
  310       format ('pitors',4x,2i5,a)
         else if (keyword(1:8) .eq. 'STRTORS ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            call getnumb (record,ic,next)
            call getnumb (record,id,next)
            ia = iclass(ia)
            ib = iclass(ib)
            ic = iclass(ic)
            id = iclass(id)
            call prmsort (4,ia,ib,ic,id,0)
            write (iprm,320)  ia,ib,ic,id,record(next:length)
  320       format ('strtors',3x,4i5,a)
         else if (keyword(1:8) .eq. 'ANGTORS ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            call getnumb (record,ic,next)
            call getnumb (record,id,next)
            ia = iclass(ia)
            ib = iclass(ib)
            ic = iclass(ic)
            id = iclass(id)
            call prmsort (4,ia,ib,ic,id,0)
            write (iprm,330)  ia,ib,ic,id,record(next:length)
  330       format ('angtors',3x,4i5,a)
         else if (keyword(1:8) .eq. 'TORTORS ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            ie = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            call getnumb (record,ic,next)
            call getnumb (record,id,next)
            call getnumb (record,ie,next)
            ia = iclass(ia)
            ib = iclass(ib)
            ic = iclass(ic)
            id = iclass(id)
            ie = iclass(ie)
            write (iprm,340)  ia,ib,ic,id,ie,record(next:length)
  340       format ('tortors',3x,5i5,a)
         else if (keyword(1:7) .eq. 'CHARGE ') then
            ia = 0
            call getnumb (record,ia,next)
            ia = itype(ia)
            write (iprm,350)  ia,record(next:length)
  350       format ('charge',4x,i5,a)
         else if (keyword(1:7) .eq. 'DIPOLE ') then
            ia = 0
            ib = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            ia = itype(ia)
            ib = itype(ib)
            write (iprm,360)  ia,ib,record(next:length)
  360       format ('dipole',4x,2i5,a)
         else if (keyword(1:8) .eq. 'DIPOLE5 ') then
            ia = 0
            ib = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            ia = itype(ia)
            ib = itype(ib)
            write (iprm,370)  ia,ib,record(next:length)
  370       format ('dipole5',3x,2i5,a)
         else if (keyword(1:8) .eq. 'DIPOLE4 ') then
            ia = 0
            ib = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            ia = itype(ia)
            ib = itype(ib)
            write (iprm,380)  ia,ib,record(next:length)
  380       format ('dipole4',3x,2i5,a)
         else if (keyword(1:8) .eq. 'DIPOLE3 ') then
            ia = 0
            ib = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            ia = itype(ia)
            ib = itype(ib)
            write (iprm,390)  ia,ib,record(next:length)
  390       format ('dipole3',3x,2i5,a)
         else if (keyword(1:10) .eq. 'MULTIPOLE ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            call getnumb (record,ic,next)
            call getnumb (record,id,next)
            ia = itype(ia)
            ib = isign(1,ib) * itype(abs(ib))
            ic = isign(1,ic) * itype(abs(ic))
            id = isign(1,id) * itype(abs(id))
            if (id .eq. 0) then
               write (iprm,400)  ia,ib,ic,record(next:length)
  400          format ('multipole',1x,3i5,a)
            else
               write (iprm,410)  ia,ib,ic,id,record(next:length)
  410          format ('multipole',1x,4i5,a)
            end if
         else if (keyword(1:9) .eq. 'POLARIZE ') then
            ia = 0
            pol = 0.0d0
            thl = 0.0d0
            kg = 0
            do j = 1, 20
               ig(j) = 0
            end do
            string = record(next:240)
            read (string,*,err=420,end=420)  ia,pol,thl,
     &                                       (ig(j),j=1,20)
  420       continue
            ia = itype(ia)
            do j = 1, maxval
               if (ig(j) .ne. 0) then
                  kg = j
                  ig(j) = itype(ig(j))
               end if
            end do
            call sort (kg,ig)
            write (iprm,430)  ia,pol,thl,(ig(j),j=1,kg)
  430       format ('polarize',2x,i5,5x,2f11.4,2x,20i5)
         else if (keyword(1:7) .eq. 'PIATOM ') then
            ia = 0
            call getnumb (record,ia,next)
            ia = iclass(ia)
            write (iprm,440)  ia,record(next:length)
  440       format ('piatom',4x,i5,a)
         else if (keyword(1:7) .eq. 'PIBOND ') then
            ia = 0
            ib = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            ia = iclass(ia)
            ib = iclass(ib)
            call prmsort (2,ia,ib,0,0,0)
            write (iprm,450)  ia,ib,record(next:length)
  450       format ('pibond',4x,2i5,a)
         else if (keyword(1:8) .eq. 'PIBOND5 ') then
            ia = 0
            ib = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            ia = iclass(ia)
            ib = iclass(ib)
            call prmsort (2,ia,ib,0,0,0)
            write (iprm,460)  ia,ib,record(next:length)
  460       format ('pibond5',3x,2i5,a)
         else if (keyword(1:8) .eq. 'PIBOND4 ') then
            ia = 0
            ib = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            ia = iclass(ia)
            ib = iclass(ib)
            call prmsort (2,ia,ib,0,0,0)
            write (iprm,470)  ia,ib,record(next:length)
  470       format ('pibond4',3x,2i5,a)
         else if (keyword(1:6) .eq. 'METAL ') then
            ia = 0
            call getnumb (record,ia,next)
            ia = iclass(ia)
            write (iprm,480)  ia,record(next:length)
  480       format ('metal',5x,i5,a)
         else if (keyword(1:8) .eq. 'BIOTYPE ') then
            ia = 0
            ib = 0
            string = record(next:240)
            read (string,*,err=490,end=490)  ia
            call getword (record,string,next)
            call getstring (record,string,next)
            string = record(next:240)
            read (string,*,err=490,end=490)  ib
  490       continue
            if (ib .gt. 0)  ib = itype(ib)
            length = min(30,max(1,59-next))
            write (iprm,500)  record(8:next)//blank(1:length),ib
  500       format ('biotype',a,i5)
         else if (length .eq. 0) then
            write (iprm,510)
  510       format ()
         else
            write (iprm,520)  record(1:length)
  520       format (a)
         end if
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine prmsort  --  reorder atom types and classes  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "prmsort" places a list of atom type or class numbers into
c     canonical order for potential energy parameter definitions
c
c
      subroutine prmsort (index,ia,ib,ic,id,ie)
      implicit none
      integer ia,ib,ic,id,ie
      integer index,temp
c
c
c     put the atom type or class numbers into canonical order
c
      if (index .eq. 2) then
         if (ia .gt. ib) then
            temp = ia
            ia = ib
            ib = temp
         end if
      else if (index .eq. 3) then
         if (ia .gt. ic) then
            temp = ia
            ia = ic
            ic = temp
         end if
      else if (index .eq. 4) then
         if (ib.gt.ic .or. (ib.eq.ic.and.ia.gt.id)) then
            temp = ib
            ib = ic
            ic = temp
            temp = ia
            ia = id
            id = temp
         end if
      else if (index .eq. 5) then
         if (ib.gt.id .or. (ib.eq.id.and.ia.gt.ie)) then
            temp = ib
            ib = id
            id = temp
            temp = ia
            ia = ie
            ie = temp
         end if
      end if
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine polesort  --  sort multipoles by atom type  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "polesort" sorts a set of atomic multipole parameters based
c     on the atom types of centers involved
c
c
      subroutine polesort (iprm)
      use sizes
      use params
      implicit none
      integer i,j,n,iprm
      integer size,next
      integer ia,ib,ic,id
      integer, allocatable :: key(:)
      integer, allocatable :: line(:)
      real*8 v1,v2,v3
      character*4 pa,pb,pc,pd
      character*16, allocatable :: list(:)
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (key(nprm))
      allocate (line(nprm))
      allocate (list(nprm))
c
c     find and store atom types for the multipole parameters
c
      n = 0
      do i = 1, nprm
         record = prmline(i)
         next = 1
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:10) .eq. 'MULTIPOLE ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            call getnumb (record,ic,next)
            call getnumb (record,id,next)
            ia = abs(ia)
            ib = abs(ib)
            ic = abs(ic)
            id = abs(id)
            size = 4
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            n = n + 1
            line(n) = i
            list(n) = pa//pb//pc//pd
         end if
      end do
c
c     sort the parameters based on the atom type numbers
c
      call sort7 (n,list,key)
c
c     format and output the sorted multipole parameters
c
      do i = 1, n
         j = line(key(i))
         record = prmline(j)
         next = 1
         call gettext (record,keyword,next)
         ia = 0
         ib = 0
         ic = 0
         id = 0
         string = record(next:240)
         read (string,*,err=20,end=20)  ia,ib,ic,id,v1
         write (iprm,10)  ia,ib,ic,id,v1
   10    format ('multipole ',4i5,6x,f11.5)
         goto 40
   20    continue
         read (string,*,err=90,end=90)  ia,ib,ic,v1
         write (iprm,30)  ia,ib,ic,v1
   30    format ('multipole ',3i5,11x,f11.5)
   40    continue
         j = j + 1
         record = prmline(j)
         read (record,*,err=90,end=90)  v1,v2,v3
         write (iprm,50)  v1,v2,v3
   50    format (36x,3f11.5)
         j = j + 1
         record = prmline(j)
         read (record,*,err=90,end=90)  v1
         write (iprm,60)  v1
   60    format (36x,f11.5)
         j = j + 1
         record = prmline(j)
         read (record,*,err=90,end=90)  v1,v2
         write (iprm,70)  v1,v2
   70    format (36x,2f11.5)
         j = j + 1
         record = prmline(j)
         read (record,*,err=90,end=90)  v1,v2,v3
         write (iprm,80)  v1,v2,v3
   80    format (36x,3f11.5)
   90    continue
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (key)
      deallocate (line)
      deallocate (list)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine biosort  --  renumber and format biotypes  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "biosort" renumbers and formats biotype parameters used to
c     convert biomolecular structure into force field atom types
c
c
      subroutine biosort (iprm)
      use sizes
      use params
      implicit none
      integer i,n,iprm
      integer next
      integer length
      integer trimtext
      integer ia,ib
      character*3 sym
      character*20 keyword
      character*30 blank
      character*240 record
      character*240 string
c
c
c     find, renumber and format the biotype parameters
c
      blank = '                              '
      n = 0
      do i = 1, nprm
         record = prmline(i)
         next = 1
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:8) .eq. 'BIOTYPE ') then
            n = n + 1
            call getnumb (record,ia,next)
            call getword (record,sym,next)
            call getstring (record,string,next)
            call getnumb (record,ib,next)
            if (ia .gt. n)  n = ia
            length = trimtext (string)
            string = '"'//string(1:length)//'"'//blank(1:30-length)
            write (iprm,10)  n,sym,string(1:32),ib
   10       format ('biotype',i8,4x,a3,5x,a32,i5)
         end if
      end do
      return
      end
