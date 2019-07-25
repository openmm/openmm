c
c
c     ##############################################################
c     ##  COPYRIGHT (C) 2008 by Chuanjie Wu & Jay William Ponder  ##
c     ##                   All Rights Reserved                    ##
c     ##############################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  program potential  --  electrostatic potential utility  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "potential" calculates the electrostatic potential for a
c     molecule at a set of grid points; optionally compares to a
c     target potential or optimizes electrostatic parameters
c
c
      program potential
      use sizes
      use atoms
      use charge
      use files
      use inform
      use iounit
      use keys
      use minima
      use mpole
      use neigh
      use output
      use potent
      use potfit
      use titles
      use units
      implicit none
      integer i,j,k
      integer ixyz,ipot
      integer igrd,icub
      integer next,mode
      integer nvar,nmodel
      integer maxpgrd
      integer nglist,nflist
      integer freeunit
      integer trimtext
      integer, allocatable :: glist(:)
      integer, allocatable :: flist(:)
      real*8 xi,yi,zi,pot
      real*8 x0,y0,z0
      real*8 xx0,xy0,xz0
      real*8 yy0,yz0,zz0
      real*8 minimum,grdmin
      real*8 potfit1
      real*8, allocatable :: xx(:)
      logical exist,query
      logical dogrid,docube
      logical domodel,dopair
      logical dotarget,dofit
      logical dofull
      logical, allocatable :: tmpchg(:)
      logical, allocatable :: tmppol(:)
      character*1 answer
      character*20 keyword
      character*240 record
      character*240 string
      character*240 xyzfile
      character*240 xyz2file
      character*240 potfile
      character*240 gridfile
      character*240 cubefile
      external potfit1
      external optsave
c
c
c     setup the computation and assign some default values
c
      call initial
      nmodel = 1
      dogrid = .false.
      docube = .false.
      domodel = .false.
      dopair = .false.
      dotarget = .false.
      dofit = .false.
c
c     perform dynamic allocation of some global arrays
c
      maxpgrd = 100000
      allocate (ipgrid(maxpgrd,maxref))
      allocate (pgrid(3,maxpgrd,maxref))
      allocate (epot(2,maxpgrd,maxref))
      allocate (gatm(maxatm))
      allocate (fatm(maxatm))
c
c     initialize target molecular dipole and quadrupole values
c
      use_dpl = .false.
      use_qdp = .false.
      fit_mpl = .true.
      fit_dpl = .true.
      fit_qdp = .true.
      do i = 1, maxref
         xdpl0(i) = 0.0d0
         ydpl0(i) = 0.0d0
         zdpl0(i) = 0.0d0
         xxqdp0(i) = 0.0d0
         xyqdp0(i) = 0.0d0
         xzqdp0(i) = 0.0d0
         yyqdp0(i) = 0.0d0
         yzqdp0(i) = 0.0d0
         zzqdp0(i) = 0.0d0
      end do
c
c     find electrostatic potential manipulation to perform
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
   20    format (/,' The TINKER Electrostatic Potential Utility Can :',
     &           //,4x,'(1) Create an Input File for Gaussian CUBEGEN',
     &           /,4x,'(2) Get QM Potential from a Gaussian CUBE File',
     &           /,4x,'(3) Calculate the Model Potential for a System',
     &           /,4x,'(4) Compare Two Model Potentials for a System',
     &           /,4x,'(5) Compare a Model Potential to a Target Grid',
     &           /,4x,'(6) Fit Electrostatic Parameters to Target Grid')
         do while (mode.lt.1 .or. mode.gt.6)
            mode = 0
            write (iout,30)
   30       format (/,' Enter the Number of the Desired Choice :  ',$)
            read (input,40,err=50,end=50)  mode
   40       format (i10)
   50       continue
         end do
      end if
      if (mode .eq. 1) then
         dogrid = .true.
      else if (mode .eq. 2) then
         docube = .true.
      else if (mode .eq. 3) then
         domodel = .true.
      else if (mode .eq. 4) then
         nmodel = 2
         dopair = .true.
      else if (mode .eq. 5) then
         dotarget = .true.
      else if (mode .eq. 6) then
         dotarget = .true.
         dofit = .true.
      end if
c
c     read electrostatic potential from a Gaussian CUBE file
c
      if (docube) then
         call nextarg (cubefile,exist)
         if (exist) then
            call basefile (cubefile)
            call suffix (cubefile,'cube','old')
            inquire (file=cubefile,exist=exist)
         end if
         do while (.not. exist)
            write (iout,60)
   60       format (/,' Enter the Gaussian CUBE File Name :  ',$)
            read (input,70)  cubefile
   70       format (a240)
            call basefile (cubefile)
            call suffix (cubefile,'cube','old')
            inquire (file=cubefile,exist=exist)
         end do
         icub = freeunit ()
         open (unit=icub,file=cubefile,status ='old')
         rewind (unit=icub)
         read (icub,80)  title
   80    format (1x,a240)
         ltitle = trimtext (title)
         read (icub,90)
   90    format ()
         read (icub,100)  n
  100    format (i5)
         read (icub,110)  npgrid(1)
  110    format (i5)
         do i = 1, n+2
            read (icub,120)
  120       format ()
         end do
         do i = 1, npgrid(1)
            read (icub,130)  record
  130       format (a240)
            read (record,*)  xi,yi,zi,pot
            pgrid(1,i,1) = xi
            pgrid(2,i,1) = yi
            pgrid(3,i,1) = zi
            epot(1,i,1) = hartree * pot
         end do
         close (unit=icub)
c
c     write the electrostatic potential to a TINKER POT file
c
         potfile = filename(1:leng)
         call suffix (potfile,'pot','new')
         open (unit=ipot,file=potfile,status ='new')
         rewind (unit=ipot)
         write (ipot,140)  npgrid(1),title(1:ltitle)
  140    format (i8,2x,a)
         do i = 1, npgrid(1)
            xi = pgrid(1,i,1)
            yi = pgrid(2,i,1)
            zi = pgrid(3,i,1)
            pot = epot(1,i,1)
            write (ipot,150)  i,xi,yi,zi,pot
  150       format (i8,3x,3f12.6,2x,f12.4)
         end do
         close (unit=ipot)
         write (iout,160)  potfile(1:trimtext(potfile))
  160    format (/,' Electrostatic Potential Written To :  ',a)
         goto 380
      end if
c
c     read first structure and get electrostatic parameters
c
      call getxyz
      call field
      call setelect
c
c     reopen the structure file and read all the structures
c
      ixyz = freeunit ()
      xyzfile = filename
      call suffix (xyzfile,'xyz','old')
      open (unit=ixyz,file=xyzfile,status ='old')
      rewind (unit=ixyz)
      call readxyz (ixyz)
      nconf = 0
      namax = n
      do while (.not. abort)
         nconf = nconf + 1
         call makeref (nconf)
         call readxyz (ixyz)
         namax = max(namax,n)
      end do
      close (unit=ixyz)
      if (nconf .gt. 1) then
         write (iout,170)  nconf
  170    format (/,' Structures Used for Potential Analysis :',i6)
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (glist(namax))
      allocate (flist(namax))
c
c     set defaults for the active grid atoms and fit atoms
c
      nglist = 0
      nflist = 0
      ngatm = namax
      nfatm = namax
      do i = 1, namax
         glist(i) = 0
         flist(i) = 0
         gatm(i) = .true.
         fatm(i) = .true.
      end do
c
c     get control parameters and target values from keyfile
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:16) .eq. 'POTENTIAL-ATOMS ') then
            read (string,*,err=180,end=180)  (glist(k),k=nglist+1,namax)
  180       continue
            do while (glist(nglist+1) .ne. 0)
               nglist = nglist + 1
               glist(nglist) = max(-namax,min(namax,glist(nglist)))
            end do
         else if (keyword(1:14) .eq. 'POTENTIAL-FIT ') then
            read (string,*,err=190,end=190)  (flist(k),k=nflist+1,namax)
  190       continue
            do while (flist(nflist+1) .ne. 0)
               nflist = nflist + 1
               flist(nflist) = max(-namax,min(namax,flist(nflist)))
            end do
         else if (keyword(1:13) .eq. 'FIX-MONOPOLE ') then
            fit_mpl = .false.
         else if (keyword(1:11) .eq. 'FIX-DIPOLE ') then
            fit_dpl = .false.
         else if (keyword(1:15) .eq. 'FIX-QUADRUPOLE ') then
            fit_qdp = .false.
         else if (keyword(1:14) .eq. 'TARGET-DIPOLE ') then
            use_dpl = .true.
            k = 1
            read (string,*,err=200,end=200)  x0,y0,z0,k
  200       continue
            xdpl0(k) = x0
            ydpl0(k) = y0
            zdpl0(k) = z0
         else if (keyword(1:18) .eq. 'TARGET-QUADRUPOLE ') then
            use_qdp = .true.
            k = 1
            read (string,*,err=210,end=210)  xx0,xy0,xz0,yy0,yz0,zz0,k
  210       continue
            xxqdp0(k) = xx0
            xyqdp0(k) = xy0
            xzqdp0(k) = xz0
            yyqdp0(k) = yy0
            yzqdp0(k) = yz0
            zzqdp0(k) = zz0
         end if
      end do
c
c     set active grid atoms to only those marked for use
c
      i = 1
      do while (glist(i) .ne. 0)
         if (i .eq. 1) then
            ngatm = 0
            do k = 1, namax
               gatm(k) = .false.
            end do
         end if
         if (glist(i) .gt. 0) then
            k = glist(i)
            if (.not. gatm(k)) then
               gatm(k) = .true.
               ngatm = ngatm + 1
            end if
            i = i + 1
         else
            do k = abs(glist(i)), abs(glist(i+1))
               if (.not. gatm(k)) then
                  gatm(k) = .true.
                  ngatm = ngatm + 1
               end if
            end do
            i = i + 2
         end if
      end do
c
c     set active fitting atoms to only those marked for use
c
      i = 1
      do while (flist(i) .ne. 0)
         if (i .eq. 1) then
            nfatm = 0
            do k = 1, namax
               fatm(k) = .false.
            end do
         end if
         if (flist(i) .gt. 0) then
            k = flist(i)
            if (.not. fatm(k)) then
               fatm(k) = .true.
               nfatm = nfatm + 1
            end if
            i = i + 1
         else
            do k = abs(flist(i)), abs(flist(i+1))
               if (.not. fatm(k)) then
                  fatm(k) = .true.
                  nfatm = nfatm + 1
               end if
            end do
            i = i + 2
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (glist)
      deallocate (flist)
c
c     generate potential grid based on the molecular surface
c
      if (.not. dotarget) then
         do i = 1, nconf
            call getref (i)
            call potgrid (i)
         end do
      end if
c
c     get name of optional second structure for comparison
c
      if (dopair) then
         call nextarg (xyz2file,exist)
         if (exist) then
            call basefile (xyz2file)
            call suffix (xyz2file,'xyz','old')
            inquire (file=xyz2file,exist=exist)
         end if
         do while (.not. exist)
            write (iout,220)
  220       format (/,' Enter Name of Second Coordinate File :  ',$)
            read (input,230)  xyz2file
  230       format (a240)
            call basefile (xyz2file)
            call suffix (xyz2file,'xyz','old')
            inquire (file=xyz2file,exist=exist)
         end do
      end if
c
c     get optional file with grid points and target potential
c
      if (dotarget) then
         call nextarg (potfile,exist)
         if (exist) then
            call basefile (potfile)
            call suffix (potfile,'pot','old')
            inquire (file=potfile,exist=exist)
         end if
         do while (.not. exist)
            write (iout,240)
  240       format (/,' Enter Target Grid/Potential File Name :  ',$)
            read (input,250)  potfile
  250       format (a240)
            call basefile (potfile)
            call suffix (potfile,'pot','old')
            inquire (file=potfile,exist=exist)
         end do
      end if
c
c     decide whether to output potential at each grid point
c
      dofull = .false.
      if (domodel .or. dopair .or. dotarget) then
         call nextarg (answer,exist)
         if (.not. exist) then
            write (iout,260)
  260       format (/,' Output Potential Value at Each Grid Point',
     &                 ' [N] :  ',$)
            read (input,270)  record
  270       format (a240)
            next = 1
            call gettext (record,answer,next)
         end if
         call upcase (answer)
         if (answer .eq. 'Y')  dofull = .true.
      end if
c
c     read grid points where potential will be computed
c
      if (dotarget) then
         ipot = freeunit ()
         open (unit=ipot,file=potfile,status='old')
         rewind (unit=ipot)
         do i = 1, nconf
            call getref (i)
            call readpot (ipot,i)
         end do
         close (unit=ipot)
      end if
c
c     output the number of potential grid points to be used
c
      do i = 1, nconf
         if (i .eq. 1) then
            write (iout,280)
  280       format ()
         end if
         if (npgrid(i) .gt. maxpgrd) then
            write (iout,290)
  290       format (' POTENTIAL  --  Too many Grid Points;',
     &                 ' Increase MAXGRID')
            call fatal
         else if (nconf .eq. 1) then
            write (iout,300)  npgrid(1)
  300       format (' Electrostatic Potential Grid Points :',6x,i10)
         else
            write (iout,310)  i,npgrid(i)
  310       format (' Potential Grid Points for Structure',i4,' :',
     &                 2x,i10)
         end if
      end do
c
c     output grid points at which to compute QM potential
c
      if (dogrid) then
         igrd = freeunit ()
         gridfile = filename(1:leng)
         call suffix (gridfile,'grid','new')
         open (unit=igrd,file=gridfile,status='new')
         do j = 1, nconf
            do i = 1, npgrid(j)
               xi = pgrid(1,i,j)
               yi = pgrid(2,i,j)
               zi = pgrid(3,i,j)
               write (igrd,320)  xi,yi,zi
  320          format (3f15.8)
            end do
         end do
         close (unit=igrd)
         write (iout,330)  gridfile(1:trimtext(gridfile))
  330    format (/,' Gaussian CUBEGEN Input Written To :   ',a)
         write (iout,340)
  340    format (/,' Next, run the Gaussian CUBEGEN program; for',
     &              ' example:',
     &           /,' cubegen 0 potential=MP2 xxx.fchk',
     &              ' xxx.cube -5 h < xxx.grid',
     &           //,' See the Gaussian documentation for additional',
     &              ' details;',
     &           /,' After CUBEGEN, rerun TINKER POTENTIAL using',
     &              ' Option 2')
      end if
c
c     get termination criterion for fitting as RMS gradient
c
      if (dofit) then
         grdmin = -1.0d0
         call nextarg (string,exist)
         if (exist)  read (string,*,err=350,end=350)  grdmin
  350    continue
         if (grdmin .le. 0.0d0) then
            write (iout,360)
  360       format (/,' Enter RMS Gradient Termination Criterion',
     &                 ' [0.1] :  ',$)
            read (input,370)  grdmin
  370       format (f20.0)
         end if
         if (grdmin .le. 0.0d0)  grdmin = 0.1d0
      end if
c
c     setup the potential computation for alternative models
c
      if (.not. dogrid) then
         do k = 1, nmodel
            ixyz = freeunit ()
            if (k .eq. 1) then
               call basefile (xyzfile)
               open (unit=ixyz,file=xyzfile,status='old')
            else
               call basefile (xyz2file)
               open (unit=ixyz,file=xyz2file,status='old')
            end if
            rewind (unit=ixyz)
            do j = 1, nconf
               call readxyz (ixyz)
               call makeref (j)
            end do
            close (unit=ixyz)
c
c     get potential for each structure and print statistics
c
            do j = 1, nconf
               call getref (j)
               call field
               call setelect
               if (use_mpole) then
                  call chkpole
                  call rotpole
               end if
               if (use_polar) then
                  domlst = .true.
                  doulst = .true.
                  call nblist
                  call induce
               end if
               do i = 1, npgrid(j)
                  xi = pgrid(1,i,j)
                  yi = pgrid(2,i,j)
                  zi = pgrid(3,i,j)
                  call potpoint (xi,yi,zi,pot)
                  epot(k,i,j) = pot
               end do
            end do
         end do
         call potstat (dofull,domodel,dopair,dotarget)
      end if
c
c     perform dynamic allocation of some global arrays
c
      if (dofit) then
         allocate (fchg(maxtyp))
         allocate (fpol(13,maxtyp))
         allocate (fitchg(maxtyp))
         allocate (fitpol(maxtyp))
c
c     perform dynamic allocation of some local arrays
c
         allocate (xx(12*nconf*namax))
         allocate (tmpchg(maxtyp))
         allocate (tmppol(maxtyp))
c
c     set parameters, run optimization, get final parameters
c
         nvar = 0
         do j = 1, maxtyp
            fitchg(j) = .false.
            fitpol(j) = .false.
         end do
         do j = 1, nconf
            call getref (j)
            call setelect
            call prmvar (nvar,xx)
         end do
         hguess = 1.0d-4
         coordtype = 'NONE'
         call ocvm (nvar,xx,minimum,grdmin,potfit1,optsave)
         nvar = 0
         do j = 1, maxtyp
            fitchg(j) = .false.
            fitpol(j) = .false.
         end do
         do j = 1, nconf
            call getref (j)
            call setelect
            next = nvar
            do k = 1, maxtyp
               tmpchg(k) = fitchg(k)
               tmppol(k) = fitpol(k)
            end do
            call varprm (nvar,xx,0,0.0d0)
            nvar = next
            do k = 1, maxtyp
               fitchg(k) = tmpchg(k)
               fitpol(k) = tmppol(k)
            end do
            call prmvar (nvar,xx)
         end do
c
c     get potential for each structure and print statistics
c
         nvar = 0
         do j = 1, maxtyp
            fitchg(j) = .false.
            fitpol(j) = .false.
         end do
         do j = 1, nconf
            call getref (j)
            call setelect
            call varprm (nvar,xx,0,0.0d0)
            if (use_mpole) then
               call chkpole
               call rotpole
            end if
            if (use_polar) then
               domlst = .true.
               doulst = .true.
               call nblist
               call induce
            end if
            do i = 1, npgrid(j)
               xi = pgrid(1,i,j)
               yi = pgrid(2,i,j)
               zi = pgrid(3,i,j)
               call potpoint (xi,yi,zi,pot)
               epot(1,i,j) = pot
            end do
         end do
         call prtfit
         call potstat (dofull,domodel,dopair,dotarget)
c
c     perform deallocation of some local arrays
c
         deallocate (xx)
         deallocate (tmpchg)
         deallocate (tmppol)
      end if
c
c     perform any final tasks before program exit
c
  380 continue
      call final
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine readpot  --  get and assign potential grid  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "readpot" gets a set of grid points and target electrostatic
c     potential values from an external disk file
c
c
      subroutine readpot (ipot,iconf)
      use sizes
      use atoms
      use katoms
      use potfit
      use ptable
      implicit none
      integer i,j,k
      integer ipot,iconf
      integer npoint,atn
      real*8 xi,yi,zi
      real*8 big,small
      real*8 r2,dist
      real*8, allocatable :: rad(:)
      character*240 record
c
c
c     read the grid points and target potential from a file
c
      npoint = 0
      read (ipot,10,err=20,end=20)  record
   10 format (a240)
      read (record,*,err=20,end=20)  npoint
   20 continue
      do i = 1, npoint
         pgrid(1,i,iconf) = 0.0d0
         pgrid(2,i,iconf) = 0.0d0
         pgrid(3,i,iconf) = 0.0d0
         epot(2,i,iconf) = 0.0d0
         read (ipot,30,err=40,end=40)  record
   30    format (a240)
         read (record,*,err=40,end=40)  k,(pgrid(j,i,iconf),j=1,3),
     &                                  epot(2,i,iconf)
   40    continue
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (rad(n))
c
c     set base atomic radii from consensus vdw values
c
      do i = 1, n
         atn = atmnum(type(i))
         rad(i) = vdwrad(atn)
         if (rad(i) .eq. 0.0d0)  rad(i) = 1.7d0
      end do
c
c     assign each grid point to atom on molecular surface
c
      big = 1000.0d0
      do i = 1, npoint
         small = big
         xi = pgrid(1,i,iconf)
         yi = pgrid(2,i,iconf)
         zi = pgrid(3,i,iconf)
         do k = 1, n
            r2 = (xi-x(k))**2 + (yi-y(k))**2 + (zi-z(k))**2
            dist = sqrt(r2) - rad(k)
            if (dist .lt. small) then
               small = dist
               ipgrid(i,iconf) = k
            end if
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (rad)
c
c     use potential grid points only for active grid atoms
c
      k = npoint
      npoint = 0
      do i = 1, k
         if (gatm(ipgrid(i,iconf))) then
            npoint = npoint + 1
            ipgrid(npoint,iconf) = ipgrid(i,iconf)
            pgrid(1,npoint,iconf) = pgrid(1,i,iconf)
            pgrid(2,npoint,iconf) = pgrid(2,i,iconf)
            pgrid(3,npoint,iconf) = pgrid(3,i,iconf)
            epot(2,npoint,iconf) = epot(2,i,iconf)
         end if
      end do
      npgrid(iconf) = npoint
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine potgrid  --  generate shells of grid points  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "potgrid" generates electrostatic potential grid points in
c     radially distributed shells based on the molecular surface
c
c
      subroutine potgrid (iconf)
      use sizes
      use atoms
      use iounit
      use katoms
      use keys
      use math
      use potfit
      use ptable
      implicit none
      integer i,j,k,m
      integer iconf,next
      integer npoint,nshell
      integer maxdot
      integer ndot,atn
      real*8 r2,rfactor
      real*8 roffset
      real*8 spacing
      real*8 density
      real*8 round
      real*8 xi,yi,zi
      real*8 xj,yj,zj
      real*8 xr,yr,zr
      real*8, allocatable :: rad(:)
      real*8, allocatable :: rad2(:)
      real*8, allocatable :: dot(:,:)
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     set default values for grid point generation parameters
c
      npoint = 0
      nshell = 4
      maxdot = 50000
      spacing = 0.35d0
      density = 4.0d0 * pi / spacing**2
      rfactor = 1.0d0
      roffset = 1.0d0
      round = 0.000001d0
c
c     check for keywords containing any altered parameters
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:17) .eq. 'POTENTIAL-SHELLS ') then
            read (string,*,err=10,end=10)  nshell
         else if (keyword(1:18) .eq. 'POTENTIAL-SPACING ') then
            read (string,*,err=10,end=10)  spacing
            density = 4.0d0 * pi / spacing**2
         else if (keyword(1:17) .eq. 'POTENTIAL-FACTOR ') then
            read (string,*,err=10,end=10)  rfactor
         else if (keyword(1:17) .eq. 'POTENTIAL-OFFSET ') then
            read (string,*,err=10,end=10)  roffset
         end if
   10    continue
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (rad(n))
      allocate (rad2(n))
c
c     get modified atomic radii from consensus vdw values
c
      do i = 1, n
         atn = atmnum(type(i))
         rad(i) = vdwrad(atn)
         if (rad(i) .eq. 0.0d0)  rad(i) = 1.7d0
         rad(i) = rfactor*rad(i) + roffset
         rad2(i) = rad(i) * rad(i)
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (dot(3,maxdot))
c
c     find points on each of the molecular surface shells
c
      do m = 1, nshell
         if (m .ne. 1) then
            do i = 1, n
               rad(i) = rad(i) + spacing
               rad2(i) = rad(i) * rad(i)
            end do
         end if
         do i = 1, n
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ndot = int(density*rad2(i))
            if (ndot .gt. maxdot) then
               write (iout,20)
   20          format (/,' POTGRID  --  Too many Surface Grid Points;',
     &                    ' Increase MAXDOT')
               call fatal
            end if
            call sphere (ndot,dot)
            do j = 1, ndot
               xj = xi + rad(i)*dot(1,j)
               yj = yi + rad(i)*dot(2,j)
               zj = zi + rad(i)*dot(3,j)
               xj = dble(nint(xj/round)) * round
               yj = dble(nint(yj/round)) * round
               zj = dble(nint(zj/round)) * round
               do k = 1, i-1
                  xr = xj - x(k)
                  yr = yj - y(k)
                  zr = zj - z(k)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .lt. rad2(k))  goto 30
               end do
               do k = i+1, n
                  xr = xj - x(k)
                  yr = yj - y(k)
                  zr = zj - z(k)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .lt. rad2(k))  goto 30
               end do
               npoint = npoint + 1
               ipgrid(npoint,iconf) = i
               pgrid(1,npoint,iconf) = xj
               pgrid(2,npoint,iconf) = yj
               pgrid(3,npoint,iconf) = zj
   30          continue
            end do
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (rad)
      deallocate (rad2)
      deallocate (dot)
c
c     use potential grid points only for active grid atoms
c
      k = npoint
      npoint = 0
      do i = 1, k
         if (gatm(ipgrid(i,iconf))) then
            npoint = npoint + 1
            ipgrid(npoint,iconf) = ipgrid(i,iconf)
            pgrid(1,npoint,iconf) = pgrid(1,i,iconf)
            pgrid(2,npoint,iconf) = pgrid(2,i,iconf)
            pgrid(3,npoint,iconf) = pgrid(3,i,iconf)
         end if
      end do
      npgrid(iconf) = npoint
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine setelect  --  assign electrostatic parameters  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "setelect" assigns partial charge, bond dipole and atomic
c     multipole parameters for the current structure, as needed
c     for computation of the electrostatic potential
c
c
      subroutine setelect
      implicit none
c
c
c     get connectivity info and make parameter assignments
c
      call attach
      call active
      call bonds
      call angles
      call torsions
      call bitors
      call rings
      call cutoffs
      call katom
      call kcharge
      call kdipole
      call kmpole
      call kpolar
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine potpoint  --  electrostatic potential at point  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "potpoint" calculates the electrostatic potential at a grid
c     point "i" as the total electrostatic interaction energy of
c     the system with a positive charge located at the grid point
c
c
      subroutine potpoint (xi,yi,zi,pot)
      use sizes
      use atoms
      use charge
      use chgpot
      use dipole
      use mpole
      use polar
      use potent
      use units
      implicit none
      integer k,kk,k1,k2
      real*8 e,ei,pot
      real*8 ec,ed,em,ep
      real*8 xi,yi,zi
      real*8 xk,yk,zk
      real*8 xr,yr,zr
      real*8 r,r2,dotk
      real*8 rk2,rkr3
      real*8 rr1,rr3,rr5
      real*8 f,fi,ci,ck
      real*8 dkx,dky,dkz
      real*8 ukx,uky,ukz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 qkx,qky,qkz
      real*8 scd,scq,scu
c
c
c     zero out charge, dipole and multipole potential terms
c
      ec = 0.0d0
      ed = 0.0d0
      em = 0.0d0
      ep = 0.0d0
c
c     set charge of probe site and electrostatic constants
c
      f = electric / dielec
      ci = 1.0d0
      fi = f * ci
c
c     calculate the charge contribution to the potential
c
      do k = 1, nion
         kk = iion(k)
         xr = x(kk) - xi
         yr = y(kk) - yi
         zr = z(kk) - zi
         r2 = xr*xr + yr* yr + zr*zr
         r = sqrt(r2)
         e = fi * pchg(k) / r
         ec = ec + e
      end do
c
c     calculate the bond dipole contribution to the potential
c
      do k = 1, ndipole
         k1 = idpl(1,k)
         k2 = idpl(2,k)
         xk = x(k2) - x(k1)
         yk = y(k2) - y(k1)
         zk = z(k2) - z(k1)
         xr = x(k1) + xk*sdpl(k) - xi
         yr = y(k1) + yk*sdpl(k) - yi
         zr = z(k1) + zk*sdpl(k) - zi
         r2 = xr*xr + yr* yr + zr*zr
         rk2 = xk*xk + yk*yk + zk*zk
         rkr3 = sqrt(rk2*r2) * r2
         dotk = xk*xr + yk*yr + zk*zr
         e = (fi/debye) * bdpl(k) * dotk / rkr3
         ed = ed + e
      end do
c
c     calculate the multipole contribution to the potential
c
      do k = 1, npole
         kk = ipole(k)
         xr = x(kk) - xi
         yr = y(kk) - yi
         zr = z(kk) - zi
         r2 = xr*xr + yr* yr + zr*zr
         r = sqrt(r2)
         ck = rpole(1,k)
         dkx = rpole(2,k)
         dky = rpole(3,k)
         dkz = rpole(4,k)
         qkxx = rpole(5,k)
         qkxy = rpole(6,k)
         qkxz = rpole(7,k)
         qkyy = rpole(9,k)
         qkyz = rpole(10,k)
         qkzz = rpole(13,k)
         if (use_polar) then
            ukx = uind(1,k)
            uky = uind(2,k)
            ukz = uind(3,k)
         else
            uind(1,k) = 0.0d0
            uind(2,k) = 0.0d0
            uind(3,k) = 0.0d0
         end if
c
c     construct some intermediate quadrupole values
c
         qkx = qkxx*xr + qkxy*yr + qkxz*zr
         qky = qkxy*xr + qkyy*yr + qkyz*zr
         qkz = qkxz*xr + qkyz*yr + qkzz*zr
c
c     calculate scalar products for permanent and induced
c
         scd = dkx*xr + dky*yr + dkz*zr
         scq = qkx*xr + qky*yr + qkz*zr
         scu = ukx*xr + uky*yr + ukz*zr
c
c     compute the potential contributions for this interaction
c
         rr1 = 1.0d0 / r
         rr3 = rr1 / r2
         rr5 = 3.0d0 * rr3 / r2
         e = ck*rr1 - scd*rr3 + scq*rr5
         ei = -scu * rr3
c
c     increment the overall multipole and polarization terms
c
         e = fi * e
         ei = fi * ei
         em = em + e
         ep = ep + ei
      end do
c
c     potential is sum of all interactions with probe site
c
      pot = ec + ed + em + ep
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  function potfit1  --  potential fit error and gradient  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "potfit1" is a service routine that computes the RMS error
c     and gradient for electrostatic parameters fit to a potential
c
c
      function potfit1 (xx,g)
      use sizes
      use atoms
      use moment
      use neigh
      use potent
      use potfit
      implicit none
      integer i,j,k,m
      integer nvar,npoint
      real*8 potfit1
      real*8 pot,eps
      real*8 e,e0
      real*8 er,ec,et
      real*8 xi,yi,zi
      real*8 cscale,tscale
      real*8 xx(*)
      real*8 g(*)
c
c
c     initialize scaling factors for error and gradient
c
      npoint = 0
      do j = 1, nconf
         npoint = npoint + npgrid(j)
      end do
      cscale = 100000000.0d0 / dble(nconf)
      tscale = 10000.0d0 / dble(nconf)
      eps = 0.000001d0
c
c     find total error by cycling over all conformations
c
      er = 0.0d0
      ec = 0.0d0
      et = 0.0d0
      nvar = 0
      do j = 1, maxtyp
         fitchg(j) = .false.
         fitpol(j) = .false.
      end do
      do j = 1, nconf
         call getref (j)
         call setelect
         call varprm (nvar,xx,0,0.0d0)
         if (use_mpole)  call rotpole
         if (use_polar) then
            domlst = .true.
            doulst = .true.
            call nblist
            call induce
         end if
c
c     get the RMS potential error summed over grid points
c
         do i = 1, npgrid(j)
            xi = pgrid(1,i,j)
            yi = pgrid(2,i,j)
            zi = pgrid(3,i,j)
            call potpoint (xi,yi,zi,pot)
            epot(1,i,j) = pot
            er = er + (epot(1,i,j)-epot(2,i,j))**2
         end do
c
c     get deviation from integral net molecular charge
c
         call momfull
         ec = ec + cscale*(netchg-dble(nint(netchg)))**2
c
c     get deviation from dipole and quadrupole targets
c
         if (use_dpl) then
            et = et + tscale*(xdpl-xdpl0(j))**2
            et = et + tscale*(ydpl-ydpl0(j))**2
            et = et + tscale*(zdpl-zdpl0(j))**2
         end if
         if (use_qdp) then
            et = et + tscale*(xxqdp-xxqdp0(j))**2
            et = et + tscale*(xyqdp-xyqdp0(j))**2
            et = et + tscale*(xzqdp-xzqdp0(j))**2
            et = et + tscale*(yyqdp-yyqdp0(j))**2
            et = et + tscale*(yzqdp-yzqdp0(j))**2
            et = et + tscale*(zzqdp-zzqdp0(j))**2
         end if
      end do
      er = sqrt(er/dble(npoint))
      potfit1 = er + ec + et
c
c     compute numerical gradient for electrostatic parameters
c
      m = nvar
      do k = 1, m
         er = 0.0d0
         ec = 0.0d0
         et = 0.0d0
         nvar = 0
         do j = 1, maxtyp
            fitchg(j) = .false.
            fitpol(j) = .false.
         end do
         do j = 1, nconf
            call getref (j)
            call setelect
            call varprm (nvar,xx,k,-0.5d0*eps)
            if (use_mpole)  call rotpole
            if (use_polar) then
               domlst = .true.
               doulst = .true.
               call nblist
               call induce
            end if
            do i = 1, npgrid(j)
               xi = pgrid(1,i,j)
               yi = pgrid(2,i,j)
               zi = pgrid(3,i,j)
               call potpoint (xi,yi,zi,pot)
               epot(1,i,j) = pot
               er = er + (epot(1,i,j)-epot(2,i,j))**2
            end do
            call momfull
            ec = ec + cscale*(netchg-dble(nint(netchg)))**2
            if (use_dpl) then
               et = et + tscale*(xdpl-xdpl0(j))**2
               et = et + tscale*(ydpl-ydpl0(j))**2
               et = et + tscale*(zdpl-zdpl0(j))**2
            end if
            if (use_qdp) then
               et = et + tscale*(xxqdp-xxqdp0(j))**2
               et = et + tscale*(xyqdp-xyqdp0(j))**2
               et = et + tscale*(xzqdp-xzqdp0(j))**2
               et = et + tscale*(yyqdp-yyqdp0(j))**2
               et = et + tscale*(yzqdp-yzqdp0(j))**2
               et = et + tscale*(zzqdp-zzqdp0(j))**2
            end if
         end do
         er = sqrt(er/dble(npoint))
         e0 = er + ec + et
         er = 0.0d0
         ec = 0.0d0
         et = 0.0d0
         nvar = 0
         do j = 1, maxtyp
            fitchg(j) = .false.
            fitpol(j) = .false.
         end do
         do j = 1, nconf
            call getref (j)
            call setelect
            call varprm (nvar,xx,k,0.5d0*eps)
            if (use_mpole)  call rotpole
            if (use_polar) then
               domlst = .true.
               doulst = .true.
               call nblist
               call induce
            end if
            do i = 1, npgrid(j)
               xi = pgrid(1,i,j)
               yi = pgrid(2,i,j)
               zi = pgrid(3,i,j)
               call potpoint (xi,yi,zi,pot)
               epot(1,i,j) = pot
               er = er + (epot(1,i,j)-epot(2,i,j))**2
            end do
            call momfull
            ec = ec + cscale*(netchg-dble(nint(netchg)))**2
            if (use_dpl) then
               et = et + tscale*(xdpl-xdpl0(j))**2
               et = et + tscale*(ydpl-ydpl0(j))**2
               et = et + tscale*(zdpl-zdpl0(j))**2
            end if
            if (use_qdp) then
               et = et + tscale*(xxqdp-xxqdp0(j))**2
               et = et + tscale*(xyqdp-xyqdp0(j))**2
               et = et + tscale*(xzqdp-xzqdp0(j))**2
               et = et + tscale*(yyqdp-yyqdp0(j))**2
               et = et + tscale*(yzqdp-yzqdp0(j))**2
               et = et + tscale*(zzqdp-zzqdp0(j))**2
            end if
         end do
         er = sqrt(er/dble(npoint))
         e = er + ec + et
         g(k) = (e-e0) / eps
      end do
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine prmvar  --  electrostatics to optimization  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "prmvar" determines the optimization values from the
c     corresponding electrostatic potential energy parameters
c
c
      subroutine prmvar (nvar,xx)
      use sizes
      use atomid
      use atoms
      use charge
      use iounit
      use mpole
      use potfit
      use units
      implicit none
      integer i,j,k
      integer ii,it
      integer nvar
      real*8 dterm,qterm
      real*8 ci,cj
      real*8 eps,big,sum
      real*8 xx(*)
      logical done
      character*17 prmtyp
c
c
c     convert dipole and quadrupole moments to atomic units
c
      dterm = 1.0d0 / bohr
      qterm = 3.0d0 / bohr**2
      do i = 1, npole
         do j = 2, 4
            pole(j,i) = dterm * pole(j,i)
         end do
         do j = 5, 13
            pole(j,i) = qterm * pole(j,i)
         end do
      end do
      dterm = bohr
      qterm = bohr**2 / 3.0d0
c
c     regularize charges and multipoles to desired precision
c
      eps = 0.00001d0
      do i = 1, nion
         pchg(i) = dble(nint(pchg(i)/eps)) * eps
      end do
      do i = 1, npole
         do j = 1, 13
            pole(j,i) = dble(nint(pole(j,i)/eps)) * eps
         end do
      end do
c
c     maintain integer net charge over partial charges
c
      k = 0
      big = 0.0d0
      sum = 0.0d0
      do i = 1, nion
         sum = sum + pchg(i)
         ci = abs(pchg(i))
         if (ci .gt. big) then
            do j = 1, n
               cj = abs(pchg(j))
               if (i.ne.j .and. ci.eq.cj)  goto 10
            end do
            k = i
            big = ci
   10       continue
         end if
      end do
      sum = sum - dble(nint(sum))
      if (k .ne. 0)  pchg(k) = pchg(k) - sum
c
c     maintain integer net charge over atomic multipoles
c
      k = 0
      big = 0.0d0
      sum = 0.0d0
      do i = 1, npole
         sum = sum + pole(1,i)
         ci = abs(pole(1,i))
         if (ci .gt. big) then
            do j = 1, n
               cj = abs(pole(1,j))
               if (i.ne.j .and. ci.eq.cj)  goto 20
            end do
            k = i
            big = ci
   20       continue
         end if
      end do
      sum = sum - dble(nint(sum))
      if (k .ne. 0)  pole(1,k) = pole(1,k) - sum
c
c     maintain traceless quadrupole at each multipole site
c
      do i = 1, npole
         sum = pole(5,i) + pole(9,i) + pole(13,i)
         big = max(abs(pole(5,i)),abs(pole(9,i)),abs(pole(13,i)))
         k = 0
         if (big .eq. abs(pole(5,i)))  k = 5
         if (big .eq. abs(pole(9,i)))  k = 9
         if (big .eq. abs(pole(13,i)))  k = 13
         if (k .ne. 0)  pole(k,i) = pole(k,i) - sum
      end do
c
c     list active atoms when not all are used in optimization
c
      if (nconf.eq.1 .and. nfatm.ne.n) then
         write (iout,30)
   30    format (/,' Atomic Parameters Included in Potential Fitting :',
     &           //,3x,'Atom',10x,'Atom Name',9x,'Atom Type',
     &              9x,'Parameters',/)
         do i = 1, nion
            ii = iion(i)
            if (fatm(ii)) then
               it = type(ii)
               prmtyp = 'Partial Charge'
               write (iout,40)  ii,name(ii),it,prmtyp
   40          format (i6,15x,a3,10x,i6,13x,a)
            end if
         end do
         do i = 1, npole
            ii = ipole(i)
            if (fatm(ii)) then
               it = type(ii)
               prmtyp = 'Atomic Multipoles'
               write (iout,50)  ii,name(ii),it,prmtyp
   50          format (i6,15x,a3,10x,i6,13x,a)
            end if
         end do
      end if
c
c     print header information for electrostatic parameters
c
      if (nvar .eq. 0) then
         write (iout,60)
   60    format (/,' Potential Fitting of Electrostatic Parameters :',
     &           //,1x,'Parameter',6x,'Atom Type',9x,'Category',
     &              12x,'Value',8x,'Fixed',/)
      end if
c
c     get optimization parameters from partial charge values
c
      do i = 1, nion
         done = .true.
         ii = iion(i)
         it = type(ii)
         if (fatm(ii))  done = .false.
         if (.not. done) then
            if (fitchg(it))  done = .true.
            fitchg(it) = .true.
         end if
         if (.not. done) then
            if (pchg(i) .ne. 0.0d0) then
               nvar = nvar + 1
               xx(nvar) = pchg(i)
               write (iout,70)  nvar,it,'Charge  ',pchg(i)
   70          format (i6,7x,i8,13x,a8,5x,f12.5)
            else
               write (iout,80)  it,'Charge  ',pchg(i)
   80          format (4x,'--',7x,i8,13x,a8,5x,f12.5,10x,'X')
            end if
         end if
      end do
c
c     get optimization parameters from atomic multipole values
c
      do i = 1, npole
         done = .true.
         ii = ipole(i)
         it = type(ii)
         if (fatm(ii))  done = .false.
         if (.not. done) then
            if (fitpol(it))  done = .true.
            fitpol(it) = .true.
         end if
         if (.not. done) then
            if (fit_mpl .and. pole(1,i).ne.0.0d0) then
               nvar = nvar + 1
               xx(nvar) = pole(1,i)
               write (iout,90)  nvar,it,'Monopole',pole(1,i)
   90          format (i6,7x,i8,13x,a8,5x,f12.5)
            else
               write (iout,100)  it,'Monopole',pole(1,i)
  100          format (4x,'--',7x,i8,13x,a8,5x,f12.5,10x,'X')
            end if
            if (fit_dpl .and. pole(2,i).ne.0.0d0) then
               nvar = nvar + 1
               xx(nvar) = dterm * pole(2,i)
               write (iout,110)  nvar,it,'X-Dipole',pole(2,i)
  110          format (i6,7x,i8,13x,a8,5x,f12.5)
            else
               write (iout,120)  it,'X-Dipole',pole(2,i)
  120          format (4x,'--',7x,i8,13x,a8,5x,f12.5,10x,'X')
            end if
            if (fit_dpl .and. pole(3,i).ne.0.0d0) then
               nvar = nvar + 1
               xx(nvar) = dterm * pole(3,i)
               write (iout,130)  nvar,it,'Y-Dipole',pole(3,i)
  130          format (i6,7x,i8,13x,a8,5x,f12.5)
            else
               write (iout,140)  it,'Y-Dipole',pole(3,i)
  140          format (4x,'--',7x,i8,13x,a8,5x,f12.5,10x,'X')
            end if
            if (fit_dpl .and. pole(4,i).ne.0.0d0) then
               nvar = nvar + 1
               xx(nvar) = dterm * pole(4,i)
               write (iout,150)  nvar,it,'Z-Dipole',pole(4,i)
  150          format (i6,7x,i8,13x,a8,5x,f12.5)
            else
               write (iout,160)  it,'Z-Dipole',pole(4,i)
  160          format (4x,'--',7x,i8,13x,a8,5x,f12.5,10x,'X')
            end if
            if (fit_qdp .and. pole(5,i).ne.0.0d0) then
               if (pole(5,i).ne.pole(9,i) .and.
     &             pole(5,i).ne.pole(13,i)) then
                  nvar = nvar + 1
                  xx(nvar) = qterm * pole(5,i)
                  write (iout,170)  nvar,it,'XX-Quad ',pole(5,i)
  170             format (i6,7x,i8,13x,a8,5x,f12.5)
               else
                  write (iout,180)    it,'XX-Quad ',pole(5,i)
  180             format (4x,'--',7x,i8,13x,a8,5x,f12.5)
               end if
            else
               write (iout,190)  it,'XX-Quad ',pole(5,i)
  190          format (4x,'--',7x,i8,13x,a8,5x,f12.5,10x,'X')
            end if
            if (fit_qdp .and. pole(6,i).ne.0.0d0) then
               nvar = nvar + 1
               xx(nvar) = qterm * pole(6,i)
               write (iout,200)  nvar,it,'XY-Quad ',pole(6,i)
  200          format (i6,7x,i8,13x,a8,5x,f12.5)
            else
               write (iout,210)  it,'XY-Quad ',pole(6,i)
  210          format (4x,'--',7x,i8,13x,a8,5x,f12.5,10x,'X')
            end if
            if (fit_qdp .and. pole(7,i).ne.0.0d0) then
               nvar = nvar + 1
               xx(nvar) = qterm * pole(7,i)
               write (iout,220)  nvar,it,'XZ-Quad ',pole(7,i)
  220          format (i6,7x,i8,13x,a8,5x,f12.5)
            else
               write (iout,230)  it,'XZ-Quad ',pole(7,i)
  230          format (4x,'--',7x,i8,13x,a8,5x,f12.5,10x,'X')
            end if
            if (fit_qdp .and. pole(9,i).ne.0.0d0) then
               if (pole(9,i).ne.pole(5,i) .and.
     &             pole(9,i).ne.pole(13,i)) then
                  nvar = nvar + 1
                  xx(nvar) = qterm * pole(9,i)
                  write (iout,240)  nvar,it,'YY-Quad ',pole(9,i)
  240             format (i6,7x,i8,13x,a8,5x,f12.5)
               else
                  write (iout,250)  it,'YY-Quad ',pole(9,i)
  250             format (4x,'--',7x,i8,13x,a8,5x,f12.5)
               end if
            else
               write (iout,260)  it,'YY-Quad ',pole(9,i)
  260          format (4x,'--',7x,i8,13x,a8,5x,f12.5,10x,'X')
            end if
            if (fit_qdp .and. pole(10,i).ne.0.0d0) then
               nvar = nvar + 1
               xx(nvar) = qterm * pole(10,i)
               write (iout,270)  nvar,it,'YZ-Quad ',pole(10,i)
  270          format (i6,7x,i8,13x,a8,5x,f12.5)
            else
               write (iout,280)  it,'YZ-Quad ',pole(10,i)
  280          format (4x,'--',7x,i8,13x,a8,5x,f12.5,10x,'X')
            end if
            if (fit_qdp .and. pole(13,i).ne.0.0d0) then
               if (pole(5,i) .eq. pole(9,i)) then
                  nvar = nvar + 1
                  xx(nvar) = qterm * pole(13,i)
                  write (iout,290)  nvar,it,'ZZ-Quad ',pole(13,i)
  290             format (i6,7x,i8,13x,a8,5x,f12.5)
               else
                  write (iout,300)  it,'ZZ-Quad ',pole(13,i)
  300             format (4x,'--',7x,i8,13x,a8,5x,f12.5)
               end if
            else
               write (iout,310)  it,'ZZ-Quad ',pole(13,i)
  310          format (4x,'--',7x,i8,13x,a8,5x,f12.5,10x,'X')
            end if
         end if
      end do
c
c     reconvert dipole and quadrupole moments to original units
c
      do i = 1, npole
         do j = 2, 4
            pole(j,i) = dterm * pole(j,i)
         end do
         do j = 5, 13
            pole(j,i) = qterm * pole(j,i)
         end do
      end do
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine varprm  --  optimization to electrostatics  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "varprm" copies the current optimization values into the
c     corresponding electrostatic potential energy parameters
c
c
      subroutine varprm (nvar,xx,ivar,eps)
      use sizes
      use atoms
      use charge
      use mpole
      use potent
      use potfit
      implicit none
      integer i,j
      integer ii,it
      integer nvar,ivar
      real*8 eps
      real*8 xx(*)
      logical done
c
c
c     translate optimization values back to partial charges
c
      do i = 1, nion
         done = .true.
         ii = iion(i)
         it = type(ii)
         if (fatm(ii))  done = .false.
         if (.not. done) then
            if (fitchg(it)) then
               done = .true.
               pchg(i) = fchg(it)
            end if
         end if
         if (.not. done) then
            if (pchg(i) .ne. 0.0d0) then
               nvar = nvar + 1
               pchg(i) = xx(nvar)
               if (ivar .eq. nvar)  pchg(i) = pchg(i) + eps
            end if
            fitchg(it) = .true.
            fchg(it) = pchg(i)
         end if
      end do
c
c     translate optimization values back to atomic multipoles
c
      do i = 1, npole
         done = .true.
         ii = ipole(i)
         it = type(ii)
         if (fatm(ii))  done = .false.
         if (.not. done) then
            if (fitpol(it)) then
               done = .true.
               do j = 1, 13
                  pole(j,i) = fpol(j,it)
               end do
            end if
         end if
         if (.not. done) then
            if (fit_mpl .and. pole(1,i).ne.0.0d0) then
               nvar = nvar + 1
               pole(1,i) = xx(nvar)
               if (ivar .eq. nvar)  pole(1,i) = pole(1,i) + eps
            end if
            if (fit_dpl .and. pole(2,i).ne.0.0d0) then
               nvar = nvar + 1
               pole(2,i) = xx(nvar)
               if (ivar .eq. nvar)  pole(2,i) = pole(2,i) + eps
            end if
            if (fit_dpl .and. pole(3,i).ne.0.0d0) then
               nvar = nvar + 1
               pole(3,i) = xx(nvar)
               if (ivar .eq. nvar)  pole(3,i) = pole(3,i) + eps
            end if
            if (fit_dpl .and. pole(4,i).ne.0.0d0) then
               nvar = nvar + 1
               pole(4,i) = xx(nvar)
               if (ivar .eq. nvar)  pole(4,i) = pole(4,i) + eps
            end if
            if (fit_qdp .and. pole(5,i).ne.0.0d0) then
               if (pole(5,i).ne.pole(9,i) .and.
     &             pole(5,i).ne.pole(13,i)) then
                  nvar = nvar + 1
                  pole(5,i) = xx(nvar)
                  if (ivar .eq. nvar)  pole(5,i) = pole(5,i) + eps
                  if (pole(9,i) .eq. pole(13,i)) then
                     pole(9,i) = -0.5d0 * pole(5,i)
                     pole(13,i) = pole(9,i)
                  end if
               end if
            end if
            if (fit_qdp .and. pole(6,i).ne.0.0d0) then
               nvar = nvar + 1
               pole(6,i) = xx(nvar)
               if (ivar .eq. nvar)  pole(6,i) = pole(6,i) + eps
               pole(8,i) = xx(nvar)
               if (ivar .eq. nvar)  pole(8,i) = pole(8,i) + eps
            end if
            if (fit_qdp .and. pole(7,i).ne.0.0d0) then
               nvar = nvar + 1
               pole(7,i) = xx(nvar)
               if (ivar .eq. nvar)  pole(7,i) = pole(7,i) + eps
               pole(11,i) = xx(nvar)
               if (ivar .eq. nvar)  pole(11,i) = pole(11,i) + eps
            end if
            if (fit_qdp .and. pole(9,i).ne.0.0d0) then
               if (pole(9,i).ne.pole(5,i) .and.
     &             pole(9,i).ne.pole(13,i)) then
                  nvar = nvar + 1
                  pole(9,i) = xx(nvar)
                  if (ivar .eq. nvar)  pole(9,i) = pole(9,i) + eps
                  if (pole(5,i) .eq. pole(13,i)) then
                     pole(5,i) = -0.5d0 * pole(9,i)
                     pole(13,i) = pole(5,i)
                  end if
               end if
            end if
            if (fit_qdp .and. pole(10,i).ne.0.0d0) then
               nvar = nvar + 1
               pole(10,i) = xx(nvar)
               if (ivar .eq. nvar)  pole(10,i) = pole(10,i) + eps
               pole(12,i) = xx(nvar)
               if (ivar .eq. nvar)  pole(12,i) = pole(12,i) + eps
            end if
            if (fit_qdp .and. pole(13,i).ne.0.0d0) then
               if (pole(5,i) .eq. pole(9,i)) then
                  nvar = nvar + 1
                  pole(13,i) = xx(nvar)
                  if (ivar .eq. nvar)  pole(13,i) = pole(13,i) + eps
                  pole(5,i) = -0.5d0 * pole(13,i)
                  pole(9,i) = pole(5,i)
               else
                  pole(13,i) = -pole(5,i) - pole(9,i)
               end if
            end if
            fitpol(it) = .true.
            do j = 1, 13
               fpol(j,it) = pole(j,i)
            end do
         end if
      end do
c
c     check chiral multipoles and rotate into global frame
c
      if (use_mpole) then
         call chkpole
         call rotpole
      end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine potstat  --  electrostatic potential statistics  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "potstat" computes and prints statistics for the electrostatic
c     potential over a set of grid points
c
c
      subroutine potstat (dofull,domodel,dopair,dotarget)
      use sizes
      use atoms
      use files
      use iounit
      use potfit
      use refer
      use titles
      implicit none
      integer i,j,k
      integer ipot,npoint
      integer freeunit
      integer trimtext
      integer, allocatable :: natm(:)
      real*8 xi,yi,zi
      real*8 pave1,pave2
      real*8 mave1,mave2
      real*8 tave,uave,rmsd
      real*8, allocatable :: patm1(:)
      real*8, allocatable :: patm2(:)
      real*8, allocatable :: rmsa(:)
      logical dofull,domodel
      logical dopair,dotarget
      character*240 potfile
c
c
c     output potential values for each model at each point
c
      if (dofull) then
         if (domodel) then
            ipot = freeunit ()
            potfile = filename(1:leng)//'.pot'
            call version (potfile,'new')
            open (unit=ipot,file=potfile,status='new')
         end if
         do j = 1, nconf
            if (nconf .eq. 1) then
               write (iout,10)
   10          format (/,' Electrostatic Potential at Each Grid',
     &                    ' Point :',
     &                 /,8x,'(Kcal/mole per unit charge)')
            else
               write (iout,20)  j
   20          format (/,' Electrostatic Potential at Grid Points',
     &                    ' for Structure',i4,' :',
     &                 /,12x,'(Kcal/mole per unit charge)')
            end if
            if (dotarget) then
               write (iout,30)
   30          format (/,3x,'Point',15x,'XYZ-Coordinates',15x,
     &                    'Potential',5x,'Target',/)
            else if (dopair) then
               write (iout,40)
   40          format (/,3x,'Point',15x,'XYZ-Coordinates',13x,
     &                    'Potential 1',3x,'Potential 2',/)
            else if (domodel) then
               write (iout,50)
   50          format (/,3x,'Point',15x,'XYZ-Coordinates',14x,
     &                    'Potential',/)
               write (ipot,60)  npgrid(j),title(1:ltitle)
   60          format (i8,2x,a)
            end if
            do i = 1, npgrid(j)
               xi = pgrid(1,i,j)
               yi = pgrid(2,i,j)
               zi = pgrid(3,i,j)
               if (dotarget .or. dopair) then
                  write (iout,70)  i,xi,yi,zi,epot(1,i,j),epot(2,i,j)
   70             format (i8,3x,3f12.6,2x,2f12.4)
               else if (domodel) then
                  write (iout,80)  i,xi,yi,zi,epot(1,i,j)
   80             format (i8,3x,3f12.6,2x,f12.4)
                  write (ipot,90)  i,xi,yi,zi,epot(1,i,j)
   90             format (i8,3x,3f12.6,2x,f12.4)
               end if
            end do
         end do
         if (domodel) then
            close (unit=ipot)
            write (iout,100)  potfile(1:trimtext(potfile))
  100       format (/,' Electrostatic Potential Written To :  ',a)
         end if
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (natm(namax))
      allocate (patm1(namax))
      allocate (patm2(namax))
      allocate (rmsa(namax))
c
c     find average electrostatic potential around each atom
c
      write (iout,110)
  110 format (/,' Average Electrostatic Potential over Atoms :',
     &        /,6x,'(Kcal/mole per unit charge)')
      if (dotarget) then
         write (iout,120)
  120    format (/,3x,'Structure',3x,'Atom',6x,'Points',
     &              6x,'Potential',8x,'Target',8x,'RMS Diff',/)
      else if (dopair) then
         write (iout,130)
  130    format (/,3x,'Structure',3x,'Atom',6x,'Points',
     &              5x,'Potential 1',4x,'Potential 2',6x,'RMS Diff',/)
      else if (domodel) then
         write (iout,140)
  140    format (/,3x,'Structure',3x,'Atom',5x,'Points',
     &              6x,'Potential',/)
      end if
      do j = 1, nconf
         call getref (j)
         do i = 1, n
            natm(i) = 0
            patm1(i) = 0.0d0
            patm2(i) = 0.0d0
            rmsa(i) = 0.0d0
         end do
         do i = 1, npgrid(j)
            k = ipgrid(i,j)
            natm(k) = natm(k) + 1
            patm1(k) = patm1(k) + epot(1,i,j)
            patm2(k) = patm2(k) + epot(2,i,j)
            rmsa(k) = rmsa(k) + (epot(1,i,j)-epot(2,i,j))**2
         end do
         do i = 1, n
            if (natm(i) .ne. 0) then
               patm1(i) = patm1(i) / dble(natm(i))
               patm2(i) = patm2(i) / dble(natm(i))
               rmsa(i) = sqrt(rmsa(i)/dble(natm(i)))
            end if
            if (gatm(i)) then
               if (dotarget .or. dopair) then
                  write (iout,150)  j,i,natm(i),patm1(i),
     &                              patm2(i),rmsa(i)
  150             format (2i9,3x,i9,3x,f12.4,3x,f12.4,3x,f12.4)
               else if (domodel) then
                  write (iout,160)  j,i,natm(i),patm1(i)
  160             format (2i9,3x,i9,3x,f12.4)
               end if
            end if
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (natm)
      deallocate (patm1)
      deallocate (patm2)
      deallocate (rmsa)
c
c     overall averages for the sets of electrostatic potentials
c
      npoint = 0
      pave1 = 0.0d0
      pave2 = 0.0d0
      mave1 = 0.0d0
      mave2 = 0.0d0
      tave = 0.0d0
      uave = 0.0d0
      rmsd = 0.0d0
      do j = 1, nconf
         npoint = npoint + npgrid(j)
         do i = 1, npgrid(j)
            pave1 = pave1 + epot(1,i,j)
            pave2 = pave2 + epot(2,i,j)
            mave1 = mave1 + abs(epot(1,i,j))
            mave2 = mave2 + abs(epot(2,i,j))
            tave = tave + epot(1,i,j) - epot(2,i,j)
            uave = uave + abs(epot(1,i,j)-epot(2,i,j))
            rmsd = rmsd + (epot(1,i,j)-epot(2,i,j))**2
         end do
      end do
      pave1 = pave1 / dble(npoint)
      pave2 = pave2 / dble(npoint)
      mave1 = mave1 / dble(npoint)
      mave2 = mave2 / dble(npoint)
      tave = tave / dble(npoint)
      uave = uave / dble(npoint)
      rmsd = sqrt(rmsd/dble(npoint))
      if (dopair) then
         write (iout,170)  pave1,mave1
  170    format (/,' Electrostatic Potential over all Grid Points :',
     &           //,' Average Potential Value for Model 1 :',10x,f12.4,
     &           /,' Average Potential Magnitude for Model 1 :',
     &              6x,f12.4)
      else
         write (iout,180)  pave1,mave1
  180    format (/,' Electrostatic Potential over all Grid Points :',
     &           //,' Average Potential Value for Model :',12x,f12.4,
     &           /,' Average Potential Magnitude for Model :',8x,f12.4)
      end if
      if (dotarget) then
         write (iout,190)  pave2,mave2,tave,uave,rmsd
  190    format (' Average Potential Value for Target :',11x,f12.4,
     &           /,' Average Potential Magnitude for Target :',7x,f12.4,
     &           /,' Average Signed Potential Difference :',10x,f12.4,
     &           /,' Average Unsigned Potential Difference :',8x,f12.4,
     &           /,' Root Mean Square Potential Difference :',8x,f12.4)
      else if (dopair) then
         write (iout,200)  pave2,mave2,tave,uave,rmsd
  200    format (' Average Potential Value for Model 2 :',10x,f12.4,
     &           /,' Average Potential Magnitude for Model 2 :',6x,f12.4,
     &           /,' Average Signed Potential Difference :',10x,f12.4,
     &           /,' Average Unsigned Potential Difference :',8x,f12.4,
     &           /,' Root Mean Square Potential Difference :',8x,f12.4)
      end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine prtfit  --  create file with optimal parameters  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "prtfit" makes a key file containing results from fitting a
c     charge or multipole model to an electrostatic potential grid
c
c
      subroutine prtfit
      use sizes
      use atoms
      use charge
      use files
      use keys
      use mpole
      use potfit
      use units
      implicit none
      integer i,j,k
      integer ii,it
      integer ix,iy,iz
      integer ikey,size
      integer freeunit
      integer trimtext
      real*8 dterm,qterm
      logical done,header
      character*240 keyfile
      character*240 record
c
c
c     open a new keyfile to contain the optimized parameters
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
         write (ikey,10)  record(1:size)
   10    format (a)
      end do
c
c     output the optimized partial charge values to the keyfile
c
      header = .true.
      do i = 1, maxtyp
         fitchg(i) = .false.
      end do
      do k = 1, nconf
         call getref (k)
         call setelect
         do i = 1, nion
            done = .true.
            ii = iion(i)
            it = type(ii)
            if (fatm(ii))  done = .false.
            if (.not. done) then
               if (fitchg(it))  done = .true.
               fitchg(it) = .true.
            end if
            if (.not. done) then
               pchg(i) = fchg(it)
               if (header) then
                  header = .false.
                  write (ikey,20)
   20             format (/,'#',/,'# Charges from Electrostatic',
     &                       ' Potential Fitting',/,'#',/)
               end if
               write (ikey,30)  it,pchg(i)
   30          format ('charge',4x,i5,10x,f11.4)
            end if
         end do
      end do
c
c     output the optimized atomic multipole values to the keyfile
c
      dterm = 1.0d0 / bohr
      qterm = 3.0d0 / bohr**2
      header = .true.
      do i = 1, maxtyp
         fitpol(i) = .false.
      end do
      do k = 1, nconf
         call getref (k)
         call setelect
         do i = 1, npole
            done = .true.
            ii = ipole(i)
            it = type(ii)
            if (fatm(ii))  done = .false.
            if (.not. done) then
               if (fitpol(it))  done = .true.
               fitpol(it) = .true.
            end if
            if (.not. done) then
               if (header) then
                  header = .false.
                  write (ikey,40)
   40             format (/,'#',/,'# Multipoles from Electrostatic',
     &                       ' Potential Fitting',/,'#',/)
               end if
               pole(1,i) = fpol(1,it)
               do j = 2, 4
                  pole(j,i) = dterm * fpol(j,it)
               end do
               do j = 5, 13
                  pole(j,i) = qterm * fpol(j,it)
               end do
               iz = zaxis(i)
               ix = xaxis(i)
               iy = yaxis(i)
               if (iy .lt. 0) then
                  yaxis(i) = -yaxis(i)
                  iy = yaxis(i)
                  pole(3,i) = -pole(3,i)
                  pole(6,i) = -pole(6,i)
                  pole(8,i) = -pole(8,i)
                  pole(10,i) = -pole(10,i)
                  pole(12,i) = -pole(12,i)
               end if
               if (iz .ne. 0)  iz = type(iz)
               if (ix .ne. 0)  ix = type(ix)
               if (iy .ne. 0)  iy = type(iy)
               if (polaxe(i) .eq. 'None') then
                  write (ikey,50)  it,pole(1,i)
   50             format ('multipole',27x,f11.5)
               else if (polaxe(i) .eq. 'Z-Only') then
                  write (ikey,60)  it,iz,pole(1,i)
   60             format ('multipole',1x,2i5,16x,f11.5)
               else if (polaxe(i) .eq. 'Z-then-X') then
                  if (yaxis(i) .eq. 0) then
                     write (ikey,70)  it,iz,ix,pole(1,i)
   70                format ('multipole',1x,3i5,11x,f11.5)
                  else
                     write (ikey,80)  it,iz,ix,iy,pole(1,i)
   80                format ('multipole',1x,4i5,6x,f11.5)
                  end if
               else if (polaxe(i) .eq. 'Bisector') then
                  if (yaxis(i) .eq. 0) then
                     write (ikey,90)  it,-iz,-ix,pole(1,i)
   90                format ('multipole',1x,3i5,11x,f11.5)
                  else
                     write (ikey,100)  it,-iz,-ix,iy,pole(1,i)
  100                format ('multipole',1x,4i5,6x,f11.5)
                  end if
               else if (polaxe(i) .eq. 'Z-Bisect') then
                  write (ikey,110)  it,iz,-ix,-iy,pole(1,i)
  110             format ('multipole',1x,4i5,6x,f11.5)
               else if (polaxe(i) .eq. '3-Fold') then
                  write (ikey,120)  it,-iz,-ix,-iy,pole(1,i)
  120             format ('multipole',1x,4i5,6x,f11.5)
               end if
               write (ikey,130)  pole(2,i),pole(3,i),pole(4,i)
  130          format (36x,3f11.5)
               write (ikey,140)  pole(5,i)
  140          format (36x,f11.5)
               write (ikey,150)  pole(8,i),pole(9,i)
  150          format (36x,2f11.5)
               write (ikey,160)  pole(11,i),pole(12,i),pole(13,i)
  160          format (36x,3f11.5)
            end if
         end do
      end do
      close (unit=ikey)
      return
      end
