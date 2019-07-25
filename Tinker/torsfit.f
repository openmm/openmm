c
c
c     ##############################################################
c     ##  COPYRIGHT (C) 2010 by Chuanjie Wu & Jay William Ponder  ##
c     ##                   All Rights Reserved                    ##
c     ##############################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  program torsfit  --  fit torsional force field parameters  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "torsfit" refines torsional force field parameters based on
c     a quantum mechanical potential surface and analytical gradient
c
c
      program torsfit
      use sizes
      use files
      use inform
      use iounit
      use keys
      implicit none
      integer i,length
      integer torbnd(10)
      logical exist,query
      character*240 record
      character*240 string
      character*240 xyzfile
c
c
c     get the Cartesian coordinates and connectivity info
c
      call initial
      call getxyz
      xyzfile = filename
      length = leng
c
c     find keyword options and setup force field parameters
c
      call getkey
      call mechanic
c
c     choose the first torsion based on the center bond atoms
c
      do i = 1, 10
         torbnd(i) = 0
      end do
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  torbnd(1),torbnd(2)
         query = .false.
      end if
   10 continue
      if (query) then
         do while (torbnd(1).eq.0 .or. torbnd(2).eq.0)
            write (iout,20)
   20       format (/,' Enter Central Atoms of the 1st Torsion : ',$)
            read (input,*,err=30,end=30)  torbnd(1),torbnd(2)
   30       continue
         end do
      end if
c
c     choose the second torsion based on the center bond atoms
c
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=40,end=40)  torbnd(3),torbnd(4)
         query = .false.
      end if
   40 continue
      if (query) then
         write(iout,50)
   50    format (/,' Enter Central Atoms of the 2nd Torsion',
     &              ' [Optional, <CR>=None] : ',$)
         read (input,60,err=70,end=70)  record
   60    format (a240)
         read (record,*,err=70,end=70)  torbnd(3),torbnd(4)
   70    continue
      end if
c
c     fit the torsional parameters based on potential surface
c
      call fittors (torbnd)
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine fittors  --  torsional parameter refinement  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "fittors" refines torsion parameters based on a quantum
c     mechanical optimized energy surface
c
c
      subroutine fittors (torbnd)
      use sizes
      use atoms
      use atomid
      use files
      use inform
      use iounit
      use keys
      use ktorsn
      use math
      use output
      use potent
      use qmstuf
      use restrn
      use scales
      use tors
      use usage
      implicit none
      integer maxfit,maxconf
      parameter (maxfit=12)
      parameter (maxconf=500)
      integer i,j,k,ii,jj,kk
      integer ia,ib,ic,id
      integer ita,itb,itc,itd
      integer itmpa,itmpb,otfix
      integer ntorfit,ntorcrs
      integer nconf,size
      integer oldleng,oldnkey
      integer istep,maxstep
      integer nvxx,ivxx
      integer ikey,nvar
      integer freeunit
      integer trimtext
      integer torbnd(*)
      integer ctorid(maxfit)
      integer ftorid(maxfit)
      integer tflg(maxfit)
      integer torcrs(4,maxfit)
      integer cflg(9*maxfit)
      integer refconf(maxconf)
      real*8 tmpa,tmpb,tv,vcon
      real*8 eqmmin,emmmin
      real*8 rms,zrms,avedl
      real*8 minimum,grdmin
      real*8 energy,torfit1
      real*8 geometry
      real*8 vxx(6*maxfit)
      real*8 vxxl(6*maxfit)
      real*8 eqm(maxconf)
      real*8 emm(maxconf)
      real*8 erqm(maxconf)
      real*8 ermm(maxconf)
      real*8 delte(maxconf)
      real*8 fwt(maxconf)
      real*8 torf(maxconf)
      real*8, allocatable :: xx(:)
      real*8 tord(6*maxfit,6*maxfit)
      real*8 mata(6*maxfit,6*maxfit)
      real*8 ftv(maxconf,maxfit)
      real*8 rftv(maxconf,maxfit)
      real*8 coeff(maxconf,6*maxfit)
      real*8 ctv(maxconf,9*maxfit)
      logical done
      logical vflg(6,maxfit)
      logical confvisited(maxconf)
      character*4 pa,pb,pc,pd
      character*16 kft(maxfit)
      character*16 kct(9*maxfit)
      character*240 record
      character*240 keyfile
      character*240 oldfilename
      character*240, allocatable :: oldkeyline(:)
      external torfit1
      external optsave
c
c
c     set initial values
c
      ntorfit = 0
      ntorcrs = 0
      otfix = ntfix
      istep = 0
      tv = 0.0d0
      vcon = 0.5d0
      do i = 1, maxfit
         ftorid(i) = 0
         tflg(i) = 0
         do j = 1, 6
            vflg(j,i) = .false.
         end do
      end do
      do i = 1, 6*maxfit
         vxx(i) = 0.0d0
         vxxl(i) = 0.1d0
         avedl = 0.0d0
         do j = 1, 6*maxfit
            tord(i,j) = 0.0d0
         end do
      end do
      do i = 1, 9*maxfit
         cflg(i) = 0
      end do
      do i = 1, maxconf
         fwt(i) = 1.0d0
         torf(i) = 0.0d0
         confvisited(i) = .false.
      end do
      do i = 1, maxconf
         do j = 1, 6*maxfit
            coeff(i,j) = 0.0d0
         end do
         refconf = 0
      end do
      grdmin = 0.01
      if (torbnd(1) .gt. torbnd(2)) then
         itmpa = torbnd(1)
         torbnd(1) = torbnd(2)
         torbnd(2) = itmpa
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (oldkeyline(nkey))
c
c     store the information from the keyfile
c
      oldfilename = filename
      oldleng = leng
      oldnkey = nkey
      do i = 1, nkey
         oldkeyline(i) = keyline(i)
      end do
c
c     check all the torsions cross the two center bond atoms
c
      write (iout,10)
   10 format (/,' Torsions Crossing the Central Bond :',/)
      do i = 1, ntors
         ia = itors(1,i)
         ib = itors(2,i)
         ic = itors(3,i)
         id = itors(4,i)
         itmpa = ib
         itmpb = ic
         if (itmpa .gt. itmpb) then
            j = itmpa
            itmpa = itmpb
            itmpb = j
         end if
         if ((torbnd(1).eq.itmpa .and. torbnd(2).eq.itmpb) .or.
     &       (torbnd(3).eq.itmpa .and. torbnd(4).eq.itmpb)) then
            ntorcrs = ntorcrs + 1
            torcrs(1,ntorcrs) = ia
            torcrs(2,ntorcrs) = ib
            torcrs(3,ntorcrs) = ic
            torcrs(4,ntorcrs) = id
            ctorid(ntorcrs) = i
            write (iout,20)  ntorcrs,ia,name(ia),ib,name(ib),
     &                       ic,name(ic),id,name(id)
   20       format (' Torsion',i5,' :',3x,4(i6,'-',a3))
         end if
      end do
c
c     choose the specific torsions for fitting
c
      write (iout,30)
   30 format (/,' Choose Torsions for Fitting from Above List :  ',$)
      read (input,40,err=50,end=50)  record
   40 format (a240)
   50 continue
      read (record,*,err=60,end=60)  (ftorid(i),i=1,ntorcrs)
   60 continue
c
c     count the torsions to be fitted
c
      do i = 1, ntorcrs
         if (ftorid(i) .gt. 0)  ntorfit = ntorfit + 1
      end do
c
c     get the number of conformations for fitting
c
      write (iout,70)
   70 format (/,' Enter Total Number of Conformations :  ',$)
      read (input,*,err=80,end=80)  nconf
   80 continue
c
c     read the QM coordinates and conformations energies
c
      do i = 1, nconf
         call readgau
         write (iout,90)  i
   90    format (/ ,' Finished Reading Conformation',i4)
         do j = 1, n
            x(j) = gx(j)
            y(j) = gy(j)
            z(j) = gz(j)
         end do
         call makeref (i)
         eqm(i) = egau
      end do
c
c     calculate the relative QM conformational energies
c
      eqmmin = eqm(1)
      do i = 2, nconf
         if (eqm(i) .lt. eqmmin)  eqmmin = eqm(i)
      end do
      write (iout,100)
  100 format ()
      do i = 1, nconf
         erqm(i) = eqm(i) - eqmmin
         write (iout,110)  i,erqm(i)
  110    format (' Relative Conformational Energy (QM)',i8,f12.4,
     &              ' Kcal/mole')
      end do
c
c     get fitting torsion type (atom classes)
c
      do i = 1, ntorfit
         j = ftorid(i)
         k = ctorid(j)
         ia = itors(1,k)
         ib = itors(2,k)
         ic = itors(3,k)
         id = itors(4,k)
         ita = class(ia)
         itb = class(ib)
         itc = class(ic)
         itd = class(id)
         size = 4
         call numeral (ita,pa,size)
         call numeral (itb,pb,size)
         call numeral (itc,pc,size)
         call numeral (itd,pd,size)
         if (itb .le. itc) then
            kft(i) = pa//pb//pc//pd
         else
            kft(i) = pd//pc//pb//pa
         end if
      end do
c
c     get all the cross torsion types
c
      do i = 1, ntorcrs
         k = ctorid(i)
         ia = itors(1,k)
         ib = itors(2,k)
         ic = itors(3,k)
         id = itors(4,k)
         ita = class(ia)
         itb = class(ib)
         itc = class(ic)
         itd = class(id)
         size = 4
         call numeral (ita,pa,size)
         call numeral (itb,pb,size)
         call numeral (itc,pc,size)
         call numeral (itd,pd,size)
         if (itb .le. itc) then
            kct(i) = pa//pb//pc//pd
         else
            kct(i) = pd//pc//pb//pa
         end if
      end do
c
c     initialize the torsion and geometry restrain parameters
c
      write (iout,120)
  120 format (/,' Initial Torsional Parameters:',/)
      nvxx = 0
      do i = 1, ntorfit
         j = ftorid(i)
         k = ctorid(j)
         done = .false.
         tflg(i) = 0
         ia = itors(1,k)
         ib = itors(2,k)
         ic = itors(3,k)
         id = itors(4,k)
         ita = class(ia)
         itb = class(ib)
         itc = class(ic)
         itd = class(id)
         write (iout,130)  ita,itb,itc,itd,tors1(1,k),tors2(1,k),
     &                     tors3(1,k),tors4(1,k),tors5(1,k),tors6(1,k)
  130    format (' torsion ',4i4,6f8.3)
         do ii = 1, i-1
            jj = ftorid(ii)
            kk = ctorid(jj)
            if (kft(i) .eq. kft(ii)) then
               done = .true.
               tflg(i) = ii
               goto 150
            end if
         end do
         do ii = 1, ntorcrs
            if (kct(ii).eq.kft(i) .and. ii.ne.j)  cflg(ii) = j
         end do
         if (abs(tors1(1,k)) .gt. 0.0d0) then
            nvxx = nvxx +1
            vxx(nvxx) = tors1(1,k)
            vflg(1,i) = .true.
         end if
         tors1(1,k) = 0.0d0
         tors1(2,k) = 0.0d0
         if (abs(tors2(1,k)) .gt. 0.0d0) then
            nvxx = nvxx +1
            vxx(nvxx) = tors2(1,k)
            vflg(2,i) = .true.
         end if
         tors2(1,k) = 0.0d0
         tors2(2,k) = 180.0d0
         if (abs(tors3(1,k)) .gt. 0.0d0) then
            nvxx = nvxx +1
            vxx(nvxx) = tors3(1,k)
            vflg(3,i) = .true.
         end if
         tors3(1,k) = 0.0d0
         tors3(2,k) = 0.0d0
         if (abs(tors4(1,k)) .gt. 0.0d0) then
            nvxx = nvxx +1
            vxx(nvxx) = tors4(1,k)
            vflg(4,i) = .true.
         end if
         tors4(1,k) = 0.0d0
         tors4(2,k) = 180.0d0
         if (abs(tors5(1,k)) .gt. 0.0d0) then
            nvxx = nvxx +1
            vxx(nvxx) = tors5(1,k)
            vflg(5,i) = .true.
         end if
         tors5(1,k) = 0.0d0
         tors5(2,k) = 0.0d0
         if (abs(tors6(1,k)) .gt. 0.0d0) then
            nvxx = nvxx +1
            vxx(nvxx) = tors6(1,k)
            vflg(6,i) = .true.
         end if
         tors6(1,k) = 0.0d0
         tors6(2,k) = 180.0d0
         ntfix = ntfix+1
         itfix(1,ntfix) = ia
         itfix(2,ntfix) = ib
         itfix(3,ntfix) = ic
         itfix(4,ntfix) = id
         tfix(1,ntfix) = 5.0d0
         write (iout,140)  ia,ib,ic,id
  140    format (' Fixed Torsion',3x,4i6)
  150    continue
      end do
c
c     print torsion flags (check duplicated torsion types)
c
      do i = 1, ntorfit
         write (iout,160)  i,tflg(i)
  160    format (/,' Fitting Torsion Number',i5,5x,'Flag',i5)
         do j = 1, 6
            write (iout,170)  i,j,vflg(j,i)
  170       format (' Variable',2i4,5x,'Variable Flag',l5)
         end do
      end do
c
c     print torsion flags for all the torsions across the bond
c
      write (iout,180)
  180 format (/,' All the Torsions Across the Bond :')
      do i = 1, ntorcrs
         k = ctorid(i)
         if (cflg(i) .gt. 0) then
            tors1(1,k) = 0.0d0
            tors2(1,k) = 0.0d0
            tors3(1,k) = 0.0d0
            tors4(1,k) = 0.0d0
            tors5(1,k) = 0.0d0
            tors6(1,k) = 0.0d0
         end if
         write (iout,190)  i,cflg(i)
  190    format (' Fitting Torsion Number',i5,5x,'Flag',i5)
      end do
c
c     add one constant variable
c
      nvxx = nvxx + 1
      vxx(nvxx) = vcon
c
c     get initial energy difference
c
      do i = 1, nconf
         call getref (i)
         kk = 0
         do j = 1, ntorfit
            k = ftorid(j)
            ia = torcrs(1,k)
            ib = torcrs(2,k)
            ic = torcrs(3,k)
            id = torcrs(4,k)
            ftv(i,j) = geometry (ia,ib,ic,id)
            write (iout,200)  i,j,ftv(i,j)
  200       format (' Fitting Torsion Value',2i5,f12.4)
            if (tflg(j) .eq. 0) then
               kk = kk+1
               tfix(2,otfix+kk) = ftv(i,j)
               tfix(3,otfix+kk) = tfix(2,otfix+kk)
            end if
         end do
         do k = 1, ntorcrs
            ia = torcrs(1,k)
            ib = torcrs(2,k)
            ic = torcrs(3,k)
            id = torcrs(4,k)
            ctv(i,k) = geometry (ia,ib,ic,id)
         end do
c
c     perform dynamic allocation of some local arrays
c
         allocate (xx(3*nuse))
c
c     scale the coordinates of each active atom
c
         nvar = 0
         do j = 1, n
            if (use(j)) then
               nvar = nvar + 1
               scale(nvar) = 12.0d0
               xx(nvar) = x(j) * scale(nvar)
               nvar = nvar + 1
               scale(nvar) = 12.0d0
               xx(nvar) = y(j) * scale(nvar)
               nvar = nvar + 1
               scale(nvar) = 12.0d0
               xx(nvar) = z(j) * scale(nvar)
            end if
         end do
c
c     make the call to the optimization routine
c
         write (iout,210)  i
  210    format (/,' Minimizing Structure',i6)
         coordtype = 'CARTESIAN'
         use_geom = .true.
         grdmin = 0.01d0
         iwrite = 0
         iprint = 0
         call lbfgs (nvar,xx,minimum,grdmin,torfit1,optsave)
c
c     unscale the final coordinates for active atoms
c
         nvar = 0
         do j = 1, n
            if (use(j)) then
               nvar = nvar + 1
               x(j) = xx(nvar) / scale(nvar)
               nvar = nvar + 1
               y(j) = xx(nvar) / scale(nvar)
               nvar = nvar + 1
               z(j) = xx(nvar) / scale(nvar)
            end if
         end do
c
c     perform deallocation of some local arrays
c
         deallocate (xx)
c
c     set the energy value for the current minimum
c
         emm(i) = energy ()
      end do
c
c     calculate relative value for each torsional angle
c
      do i = 1, nconf
         do j = 1, ntorfit
            rftv(i,j) = ftv(i,j) - ftv(i,1)
         end do
      end do
c
c     calculate the relative MM energies
c
      emmmin = emm(1)
      do i = 2, nconf
         if (emm(i) .lt. emmmin)  emmmin = emm(i)
      end do
c
c     calculate the energy difference and RMS
c
      rms = 0.0d0
      zrms = 0.0d0
      write (iout,220)
  220 format ()
      do i = 1, nconf
         ermm (i) = emm(i) - emmmin
         delte (i) = erqm (i) - ermm(i)
         rms = rms + delte(i)*delte(i)
         write (iout,230)  i,ermm(i)
  230    format (' Relative Conformational Energy (MM)',i8,f12.4,
     &              ' Kcal/mole')
      end do
      rms = sqrt(rms/dble(nconf))
      zrms = rms
      write (iout,240)  rms
  240 format (/,' Energy RMS Difference :',8x,f12.4)
c
c     calculate the weights
c
c     do i = 1, nconf
c        do j = 1, nconf
c           if (.not. confvisited(j)) then
c              tmpa = ftv(j,1)
c              itmpa = j
c              confvisited(j) = .true.
c              goto 241
c           end if
c        end do
c 241    continue
c        do j = 1, nconf
c           if (ftv(j,1).lt.tmpa .and. .not.confvisited(j)) then
c              confvisited(itmpa) = .false.
c              itmpa = j
c              tmpa = ftv(j,1)
c           end if
c        end do
c        refconf(itmpa) = i
c        confvisited(itmpa) = .true.
c        write (iout,242)  itmpa,refconf(itmpa)
c 242    format (i8,' <===> ',i8)
c     end do
      if (nconf .gt. 1 .and. torbnd(3) .eq. 0) then
         if ((ftv(nconf,1)+180.0d0) .lt. 1.0d0
     &           .and. ftv(nconf-1,1) .gt. 0.0d0)
     &       ftv(nconf,1) = 180.0d0
         tmpa = erqm(2) - erqm(1)
         tmpb = (ftv(2,1)-ftv(1,1)) / radian
         fwt(1) = 1.0d0 / sqrt(1.0d0+(tmpa/tmpb)**2)
         if (nconf .gt. 2) then
            do i = 2, nconf-1
               tmpa = erqm(i+1) - erqm(i-1)
               tmpb = (ftv(i+1,1) - ftv(i-1,1))/radian
               fwt(i) = 1.0d0 / sqrt(1.0d0+(tmpa/tmpb)**2)
            end do
         end if
         tmpa = erqm(nconf) - erqm(nconf-1)
         tmpb = (ftv(nconf,1) - ftv(nconf-1,1))/radian
         fwt(nconf) = 1.0d0 / sqrt(1.0d0+(tmpa/tmpb)**2)
      end if
      write (iout,250)
  250 format ()
      do i = 1, nconf
         write (iout,260)  i,fwt(i)
  260    format (' Conformation',i5,5x,'Weight',f8.4)
      end do
c
c     set initial values for torsions to be fitted
c
      ivxx = 0
      do i = 1, ntorfit
         j = ftorid(i)
         k = ctorid(j)
         do ii = 1, 6
            if (vflg(ii,i) .and. tflg(i).eq.0) then
               ivxx = ivxx + 1
               if (ii .eq. 1) then
                  tors1(1,k) = vxx(ivxx)
               else if (ii .eq. 2) then
                  tors2(1,k) = vxx(ivxx)
               else if (ii .eq. 3) then
                  tors3(1,k) = vxx(ivxx)
               else if (ii .eq. 4) then
                  tors4(1,k) = vxx(ivxx)
               else if (ii .eq. 5) then
                  tors5(1,k) = vxx(ivxx)
               else if (ii .eq. 6) then
                  tors6(1,k) = vxx(ivxx)
               end if
               do jj = 1, ntorfit
                  kk = ctorid(ftorid(jj))
                  if (tflg(jj) .eq. i) then
                     if (ii .eq. 1) then
                        tors1(1,kk) = vxx(ivxx)
                     else if (ii .eq. 2) then
                        tors2(1,kk) = vxx(ivxx)
                     else if (ii .eq. 3) then
                        tors3(1,kk) = vxx(ivxx)
                     else if (ii .eq. 4) then
                        tors4(1,kk) = vxx(ivxx)
                     else if (ii .eq. 5) then
                        tors5(1,kk) = vxx(ivxx)
                     else if (ii .eq. 6) then
                        tors6(1,kk) = vxx(ivxx)
                     end if
                  end if
               end do
            end if
         end do
         ivxx = ivxx + 1
         vcon = vxx(ivxx)
      end do
c
c     fitting the torsion parameters
c
      write (iout,270)
  270 format ()
      maxstep = 1
      avedl = 0.5d0
      do while (avedl.gt.0.1d0 .and. istep.lt.maxstep)
         do i = 1, nconf
            ivxx = 0
            torf(i) = 0.0d0
            do j = 1, ntorfit
               jj = ftorid(j)
               kk = ctorid(jj)
               ia = itors(1,kk)
               ib = itors(2,kk)
               ic = itors(3,kk)
               id = itors(4,kk)
               ita = class(ia)
               itb = class(ib)
               itc = class(ic)
               itd = class(id)
               tv = ftv(i,j) / radian
               tmpa = tors1(1,kk)*(1+cos(tv))
     &                    + tors2(1,kk)*(1-cos(2*tv))
     &                    + tors3(1,kk)*(1+cos(3*tv))
     &                    + tors4(1,kk)*(1-cos(4*tv))
     &                    + tors5(1,kk)*(1+cos(5*tv))
     &                    + tors6(1,kk)*(1-cos(6*tv))
               torf(i) = torf(i) + 0.5*tmpa
               do ii = 1, 6
                  if (vflg(ii,j) .and. tflg(j).eq.0) then
                     ivxx = ivxx +1
                     coeff(i,ivxx) = 0.5*(1+(-1)**(ii+1)
     &                                  *cos(dble(ii)*tv))
                     do k = 1, ntorcrs
                        if (cflg(k).gt.0 .and. cflg(k).eq.jj) then
                           coeff(i,ivxx) = coeff(i,ivxx)
     &                                     +0.5*(1+(-1)**(ii+1)
     &                               *cos(dble(ii)*ctv(i,k)/radian))
                        end if
                     end do
                     write (iout,280)  i,ivxx,coeff(i,ivxx)
  280                format (' Derivative :',5x,2i4,f8.4)
                  end if
               end do
            end do
            torf(i) = torf(i) + vcon - delte(i)
            ivxx = ivxx + 1
            coeff(i,ivxx) = 1.0d0
            write (iout,290)  i,torf(i)
  290       format (' Energy Difference :',i8,f12.4)
         end do
c
c     set matrix elements for matrix A
c
         do i = 1, nvxx
            do j = 1, nvxx
               tord(i,j) = 0.0d0
               do k = 1, nconf
                  tord(i,j) = tord(i,j) + coeff(k,i)*coeff(k,j)*fwt(k)
               end do
            end do
         end do
c
c     print the matrix A elements
c
         write (iout,300)  nvxx
  300    format (/,' Total Variable Number ',i8)
         write (iout,310)
  310    format (/,' Matrix A Elements :')
         do i = 1, nvxx
            do j = 1, nvxx
               mata(i,j) = tord(i,j)
            end do
            write (iout,320)  (mata(i,j),j=1,nvxx)
  320       format (1x,5f12.4)
         end do
c
c     multiply vector: Yi * Coeff * Weight
c
         do i = 1, nvxx
            torf(i) = 0.0d0
            do j = 1, nconf
               torf(i) = torf(i) + delte(j)*fwt(j)*coeff(j,i)
            end do
         end do
         do i = 1, nvxx
            mata(i,nvxx+1) = torf(i)
         end do
c
c     solve the linear equations via Gauss-Jordan elimination
c
         call gaussjordan (nvxx,mata)
c
c     get new torsion force constants
c
         do i = 1, nvxx
            vxx(i) = mata(i,nvxx+1)
         end do
         ivxx = 0
         do i = 1, ntorfit
            j = ftorid(i)
            k = ctorid(j)
            do ii = 1, 6
               if (vflg(ii,i) .and. tflg(i).eq.0) then
                  ivxx = ivxx + 1
                  if (ii .eq. 1) then
                     tors1(1,k) = vxx(ivxx)
                  else if (ii .eq. 2) then
                     tors2(1,k) = vxx(ivxx)
                  else if (ii .eq. 3) then
                     tors3(1,k) = vxx(ivxx)
                  else if (ii .eq. 4) then
                     tors4(1,k) = vxx(ivxx)
                  else if (ii .eq. 5) then
                     tors5(1,k) = vxx(ivxx)
                  else if (ii .eq. 6) then
                     tors6(1,k) = vxx(ivxx)
                  end if
                  do jj = 1, ntorcrs
                     kk = ctorid(jj)
                     if (cflg(j).gt.0 .and. cflg(jj).eq.j) then
                        if (ii .eq. 1) then
                           tors1(1,kk) = vxx(ivxx)
                        else if (ii .eq. 2) then
                           tors2(1,kk) = vxx(ivxx)
                        else if (ii .eq. 3) then
                           tors3(1,kk) = vxx(ivxx)
                        else if (ii .eq. 4) then
                           tors4(1,kk) = vxx(ivxx)
                        else if (ii .eq. 5) then
                           tors5(1,kk) = vxx(ivxx)
                        else if (ii .eq. 6) then
                           tors6(1,kk) = vxx(ivxx)
                        end if
                     end if
                  end do
               end if
            end do
            ivxx = ivxx + 1
            vcon = vxx(ivxx)
         end do
         istep = istep + 1
      end do
c
c     validate the fitted results
c
      write (iout,330)
  330 format ()
      do i = 1, nconf
         call getref (i)
         kk = 0
         do j = 1, ntorfit
            k = ftorid(j)
            ia = torcrs(1,k)
            ib = torcrs(2,k)
            ic = torcrs(3,k)
            id = torcrs(4,k)
            ftv(i,j) = geometry (ia,ib,ic,id)
            if (tflg(j) .eq. 0) then
               kk = kk + 1
               tfix(2,otfix+kk) = ftv(i,j)
               tfix(3,otfix+kk) = tfix(2,otfix+kk)
            end if
         end do
c
c     perform dynamic allocation of some local arrays
c
         allocate (xx(3*nuse))
c
c     scale the coordinates of each active atom
c
         nvar = 0
         do j = 1, n
            if (use(j)) then
               nvar = nvar + 1
               xx(nvar) = x(j) * scale(nvar)
               nvar = nvar + 1
               xx(nvar) = y(j) * scale(nvar)
               nvar = nvar + 1
               xx(nvar) = z(j) * scale(nvar)
            end if
         end do
c
c     make the call to the optimization routine
c
         write (iout,340)  i
  340    format (' Minimizing Structure',i5,2x,'with New Parameters')
         coordtype = 'CARTESIAN'
         call lbfgs (nvar,xx,minimum,grdmin,torfit1,optsave)
c
c     unscale the final coordinates for active atoms
c
         nvar = 0
         do j = 1, n
            if (use(j)) then
               nvar = nvar + 1
               x(j) = xx(nvar) / scale(nvar)
               nvar = nvar + 1
               y(j) = xx(nvar) / scale(nvar)
               nvar = nvar + 1
               z(j) = xx(nvar) / scale(nvar)
            end if
         end do
c
c     perform deallocation of some local arrays
c
         deallocate (xx)
c
c     set the energy value for the current minimum
c
         emm(i) = energy ()
      end do
c
c     calculate the relative MM energies
c
      emmmin = emm(1)
      do i = 2, nconf
         if (emm(i) .lt. emmmin)  emmmin = emm(i)
      end do
c
c     calculate the energy difference and RMS
c
      rms = 0.0d0
      write (iout,350)
  350 format ()
      do i = 1, nconf
         ermm (i) = emm(i) - emmmin
         delte (i) = erqm (i) - ermm(i)
         rms = rms + delte(i)*delte(i)
         write (iout,360)  i,ermm(i)
  360    format (' Relative Conformational Energy (MM)',i8,f12.4,
     &              ' Kcal/mole')
      end do
      rms = sqrt(rms/dble(nconf))
      write (iout,370)  rms
  370 format (/,' Energy RMS With Fitting Parmeters :',8x,f12.4)
      if (rms .gt. zrms ) then
         write (iout,380)  zrms
  380    format (/,' Annihilating the Torsions is Preferable',
     &           /,' Final RMS :',f12.6,' Kcal/mole',/)
      end if
c
c     output keyfile information with the fitted parameters
c
      filename = oldfilename
      leng = oldleng
      nkey = oldnkey
      do i = 1, nkey
         keyline(i) = oldkeyline(i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (oldkeyline)
c
c     output some definitions and parameters to a keyfile
c
      ikey = freeunit ()
      keyfile = filename(1:leng)//'.key'
      call version (keyfile,'new')
      open (unit=ikey,file=keyfile,status='new')
c
c     copy the contents of any previously existing keyfile
c
      do i = 1, nkey
         record = keyline(i)
         size = trimtext (record)
         write (ikey,390)  record(1:size)
  390    format (a)
      end do
c
c     list the valence parameters
c
      write (ikey,400)
  400 format (/,'#',/,'# Results of Valence Parameter Fitting',
     &        /,'#',/)
      write (iout,410)
  410 format (/,' Optimized Torsional Parameters:',/)
      do i = 1, ntorfit
         if (tflg(i) .eq. 0) then
            j = ftorid(i)
            k = ctorid(j)
            ia = itors(1,k)
            ib = itors(2,k)
            ic = itors(3,k)
            id = itors(4,k)
            ita = class(ia)
            itb = class(ib)
            itc = class(ic)
            itd = class(id)
            if (rms .gt. zrms) then
               tors1(1,k) = 0.0d0
               tors2(1,k) = 0.0d0
               tors3(1,k) = 0.0d0
            end if
            write (iout,420)  ita,itb,itc,itd,tors1(1,k),
     &                        tors2(1,k),tors3(1,k)
  420       format (' torsion ',4i4,f8.3,' 0.0 1 ',f8.3,
     &                 ' 180.0 2 ',f8.3,' 0.0 3')
            write (ikey,430)  ita,itb,itc,itd,tors1(1,k),
     &                        tors2(1,k),tors3(1,k)
  430       format (' torsion ',4i4,f8.3,' 0.0 1 ',f8.3,
     &                 ' 180.0 2 ',f8.3,' 0.0 3')
         end if
      end do
      close (unit=ikey)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine gaussjordan  --  Gauss-Jordan elimination  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "gaussjordan" solves a system of linear equations by using
c     the method of Gaussian elimination with partial pivoting
c
c
      subroutine gaussjordan (n,a)
      use iounit
      implicit none
      integer maxfit
      parameter (maxfit=12)
      integer i,j,k,l,n
      real*8 t,av
      real*8 a(6*maxfit,*)
c
c
c     perform the Gauss-Jordan elimination procedure
c
      do k = 1, n-1
         av = 0.0d0
         do i = k, n
            if (abs(a(i,k)) .gt. abs(av)) then
               av = a(i,k)
               l = i
            end if
         end do
         if (abs(av) .lt. 1.0d-8) then
            write (iout,10)
   10       format (/,' GAUSSJORDAN  --  Singular Coefficient Matrix')
            call fatal
         end if
         if (l .ne. k) then
            do j = k, n+1
               t = a(k,j)
               a(k,j) = a(l,j)
               a(l,j) = t
            end do
         end if
         av = 1.0d0 / av
         do j = k+1, n+1
            a(k,j) = a(k,j) * av
            do i = k+1, n
               a(i,j) = a(i,j) - a(i,k)*a(k,j)
            end do
         end do
      end do
      a(n,n+1) = a(n,n+1) / a(n,n)
      do k = 1, n-1
         i = n - k
         av = 0.0d0
         do j = i+1, n
            av = av + a(i,j)*a(j,n+1)
         end do
         a(i,n+1) = a(i,n+1) - av
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  function torfit1  --  energy and gradient for minimize  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "torfit1" is a service routine that computes the energy and
c     gradient for a low storage BFGS optimization in Cartesian
c     coordinate space
c
c
      function torfit1 (xx,g)
      use sizes
      use atoms
      use scales
      use usage
      implicit none
      integer i,nvar
      real*8 torfit1,e
      real*8 energy,eps
      real*8 xx(*)
      real*8 g(*)
      real*8, allocatable :: derivs(:,:)
      logical analytic
      external energy
c
c
c     use either analytical or numerical gradients
c
      analytic = .true.
      eps = 0.00001d0
c
c     translate optimization parameters to atomic coordinates
c
      nvar = 0
      do i = 1, n
         if (use(i)) then
            nvar = nvar + 1
            x(i) = xx(nvar) / scale(nvar)
            nvar = nvar + 1
            y(i) = xx(nvar) / scale(nvar)
            nvar = nvar + 1
            z(i) = xx(nvar) / scale(nvar)
         end if
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (derivs(3,n))
c
c     compute and store the energy and gradient
c
      if (analytic) then
         call gradient (e,derivs)
      else
         e = energy ()
         call numgrad (energy,derivs,eps)
      end if
      torfit1 = e
c
c     store Cartesian gradient as optimization gradient
c
      nvar = 0
      do i = 1, n
         if (use(i)) then
            nvar = nvar + 1
            g(nvar) = derivs(1,i) / scale(nvar)
            nvar = nvar + 1
            g(nvar) = derivs(2,i) / scale(nvar)
            nvar = nvar + 1
            g(nvar) = derivs(3,i) / scale(nvar)
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (derivs)
      return
      end
