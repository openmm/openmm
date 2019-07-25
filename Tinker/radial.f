c
c
c     ##############################################################
c     ##  COPYRIGHT (C) 1995 by Yong Kong and Jay William Ponder  ##
c     ##                   All Rights Reserved                    ##
c     ##############################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  program radial  --  compute radial distribution function  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "radial" finds the radial distribution function for a specified
c     pair of atom types via analysis of a set of coordinate frames
c
c
      program radial
      use sizes
      use argue
      use atomid
      use atoms
      use bound
      use boxes
      use files
      use inform
      use iounit
      use limits
      use math
      use molcul
      use potent
      implicit none
      integer i,j,k,iarc
      integer nframe,iframe
      integer freeunit,next
      integer molj,molk
      integer numj,numk
      integer typej,typek
      integer start,stop
      integer step,skip
      integer nbin,bin
      integer, allocatable :: hist(:)
      real*8 xj,yj,zj
      real*8 dx,dy,dz
      real*8 rjk,rmax,width
      real*8 rlower,rupper
      real*8 factor,pairs
      real*8 volume,expect
      real*8, allocatable :: gr(:)
      real*8, allocatable :: gs(:)
      logical exist,query
      logical intramol
      character*1 answer
      character*3 namej,namek
      character*6 labelj,labelk
      character*240 arcfile
      character*240 record
      character*240 string
c
c
c     perform the standard initialization functions
c
      call initial
c
c     try to get a filename from the command line arguments
c
      call nextarg (arcfile,exist)
      if (exist) then
         call basefile (arcfile)
         call suffix (arcfile,'arc','old')
         inquire (file=arcfile,exist=exist)
      end if
c
c     ask for the user specified input structure filename
c
      do while (.not. exist)
         write (iout,10)
   10    format (/,' Enter Coordinate Archive File Name :  ',$)
         read (input,20)  arcfile
   20    format (a240)
         call basefile (arcfile)
         call suffix (arcfile,'arc','old')
         inquire (file=arcfile,exist=exist)
      end do
c
c     read the first coordinate set in the archive
c
      iarc = freeunit ()
      open (unit=iarc,file=arcfile,status='old')
      call readxyz (iarc)
c
c     get the unitcell parameters and number of molecules
c
      call unitcell
      call molecule
c
c     set cutoffs small to enforce use of minimum images
c
      use_vdw = .true.
      use_ct = .true.
      use_charge = .false.
      use_dipole = .false.
      use_mpole = .false.
      use_ewald = .false.
      vdwcut = 0.01d0
      call lattice
c
c     get numbers of the coordinate frames to be processed
c
      start = 1
      stop = 100000
      step = 1
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=30,end=30)  start
         query = .false.
      end if
      call nextarg (string,exist)
      if (exist)  read (string,*,err=30,end=30)  stop
      call nextarg (string,exist)
      if (exist)  read (string,*,err=30,end=30)  step
   30 continue
      if (query) then
         write (iout,40)
   40    format (/,' Numbers of First & Last Frame and Step',
     &              ' Increment :  ',$)
         read (input,50)  record
   50    format (a240)
         read (record,*,err=60,end=60)  start,stop,step
   60    continue
      end if
c
c     get the names of the atoms to be used in rdf computation
c
      call nextarg (labelj,exist)
      call nextarg (labelk,exist)
      if (.not. exist) then
         write (iout,70)
   70    format (/,' Enter 1st & 2nd Atom Names or Type Numbers :  ',$)
         read (input,80)  record
   80    format (a240)
         next = 1
         call gettext (record,labelj,next)
         call gettext (record,labelk,next)
      end if
c
c     convert the labels to either atom names or type numbers
c
      namej = '   '
      typej = -1
      read (labelj,*,err=90,end=90)  typej
   90 continue
      if (typej .le. 0) then
         next = 1
         call gettext (labelj,namej,next)
      end if
      namek = '   '
      typek = -1
      read (labelk,*,err=100,end=100)  typek
  100 continue
      if (typek .le. 0) then
         next = 1
         call gettext (labelk,namek,next)
      end if
c
c     get maximum distance from input or minimum image convention
c
      if (.not. use_bounds) then
         rmax = -1.0d0
         query = .true.
         call nextarg (string,exist)
         if (exist) then
            read (string,*,err=110,end=110)  rmax
            query = .false.
         end if
  110    continue
         if (query) then
            write (iout,120)
  120       format (/,' Enter Maximum Distance to Accumulate',
     &                 ' [10.0 Ang] :  ',$)
            read (input,130)  rmax
  130       format (f20.0)
         end if
         if (rmax .le. 0.0d0)  rmax = 10.0d0
      else if (octahedron) then
         rmax = (sqrt(3.0d0)/4.0d0) * xbox
         rmax = 0.95d0 * rmax
      else
         rmax = min(xbox2*beta_sin*gamma_sin,ybox2*gamma_sin,
     &                         zbox2*beta_sin)
         rmax = 0.95d0 * rmax
      end if
c
c     get the desired width of the radial distance bins
c
      width = -1.0d0
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=140,end=140)  width
         query = .false.
      end if
  140 continue
      if (query) then
         write (iout,150)
  150    format (/,' Enter Width of Distance Bins [0.01 Ang] :  ',$)
         read (input,160)  width
  160    format (f20.0)
      end if
      if (width .le. 0.0d0)  width = 0.01d0
c
c     decide whether to restrict to intermolecular atom pairs
c
      intramol = .false.
      call nextarg (answer,exist)
      if (.not. exist) then
         write (iout,170)
  170    format (/,' Include Intramolecular Pairs in Distribution',
     &              ' [N] :  ',$)
         read (input,180)  record
  180    format (a240)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .eq. 'Y')  intramol = .true.
c
c     count the number of coordinate frames in the archive file
c
      abort = .false.
      rewind (unit=iarc)
      nframe = 0
      do while (.not. abort)
         call readxyz (iarc)
         nframe = nframe + 1
      end do
      nframe = nframe - 1
      rewind (unit=iarc)
      stop = min(nframe,stop)
      nframe = (stop-start)/step + 1
      write (iout,190)  nframe
  190 format (/,' Number of Coordinate Frames :',i14)
c
c     set the number of distance bins to be accumulated
c
      nbin = int(rmax/width)
      write (iout,200)  nbin
  200 format (' Number of Distance Bins :',i18)
c
c     perform dynamic allocation of some local arrays
c
      allocate (hist(nbin))
      allocate (gr(nbin))
      allocate (gs(nbin))
c
c     zero out the distance bins and distribution functions
c
      do i = 1, nbin
         hist(i) = 0
         gr(i) = 0.0d0
         gs(i) = 0.0d0
      end do
c
c     get the archived coordinates for each frame in turn
c
      write (iout,210)
  210 format (/,' Reading the Coordinates Archive File :',/)
      nframe = 0
      iframe = start
      skip = start
      do while (iframe.ge.start .and. iframe.le.stop)
         do j = 1, skip-1
            call readxyz (iarc)
         end do
         iframe = iframe + step
         skip = step
         call readxyz (iarc)
         if (.not. abort) then
            nframe = nframe + 1
            if (mod(nframe,100) .eq. 0) then
               write (iout,220)  nframe
  220          format (4x,'Processing Coordinate Frame',i13)
            end if
            do j = 1, n
               if (name(j).eq.namej .or. type(j).eq.typej) then
                  xj = x(j)
                  yj = y(j)
                  zj = z(j)
                  molj = molcule(j)
                  do k = 1, n
                     if (name(k).eq.namek .or. type(k).eq.typek) then
                        if (j .ne. k) then
                           molk = molcule(k)
                           if (intramol .or. molj.ne.molk) then
                              dx = x(k) - xj
                              dy = y(k) - yj
                              dz = z(k) - zj
                              call image (dx,dy,dz)
                              rjk = sqrt(dx*dx + dy*dy + dz*dz)
                              bin = int(rjk/width) + 1
                              if (bin .le. nbin)
     &                           hist(bin) = hist(bin) + 1
                           end if
                        end if
                     end if
                  end do
               end if
            end do
         end if
      end do
c
c     ensure a valid frame is loaded and report total frames
c
      if (abort) then
         rewind (unit=iarc)
         call readxyz (iarc)
      end if
      close (unit=iarc)
      if (mod(nframe,100) .ne. 0) then
         write (iout,230)  nframe
  230    format (4x,'Processing Coordinate Frame',i13)
      end if
c
c     count the number of occurrences of each atom type
c
      numj = 0
      numk = 0
      do i = 1, n
         if (name(i).eq.namej .or. type(i).eq.typej)  numj = numj + 1
         if (name(i).eq.namek .or. type(i).eq.typek)  numk = numk + 1
      end do
c
c     normalize the distance bins to give radial distribution
c
      if (numj.ne.0 .and. numk.ne.0) then
         factor = (4.0d0/3.0d0) * pi * dble(nframe)
         if (use_bounds) then
            pairs = dble(numj) * dble(numk)
            volume = (gamma_sin*gamma_term) * xbox * ybox * zbox
            if (octahedron)  volume = 0.5d0 * volume
            factor = factor * pairs / volume
         end if
         do i = 1, nbin
            rupper = dble(i) * width
            rlower = rupper - width
            expect = factor * (rupper**3 - rlower**3)
            gr(i) = dble(hist(i)) / expect
         end do
      end if
c
c     find the 5th degree polynomial smoothed distribution function
c
      if (nbin .ge. 5) then
         gs(1) = (69.0d0*gr(1) + 4.0d0*gr(2) - 6.0d0*gr(3)
     &             + 4.0d0*gr(4) - gr(5)) / 70.0d0
         gs(2) = (2.0d0*gr(1) + 27.0d0*gr(2) + 12.0d0*gr(3)
     &             - 8.0d0*gr(4) + 2.0d0*gr(5)) / 35.0d0
         do i = 3, nbin-2
            gs(i) = (-3.0d0*gr(i-2) + 12.0d0*gr(i-1) + 17.0d0*gr(i)
     &                + 12.0d0*gr(i+1) - 3.0d0*gr(i+2)) / 35.0d0
         end do
         gs(nbin-1) = (2.0d0*gr(nbin-4) - 8.0d0*gr(nbin-3)
     &                  + 12.0d0*gr(nbin-2) + 27.0d0*gr(nbin-1)
     &                       + 2.0d0*gr(nbin)) / 35.0d0
         gs(nbin) = (-gr(nbin-4) + 4.0d0*gr(nbin-3) - 6.0d0*gr(nbin-2)
     &                + 4.0d0*gr(nbin-1) + 69.0d0*gr(nbin)) / 70.0d0
         do i = 1, nbin
            gs(i) = max(0.0d0,gs(i))
         end do
      end if
c
c     output the final radial distribution function results
c
      write (iout,240)  labelj,labelk
  240 format (/,' Pairwise Radial Distribution Function :'
     &        //,7x,'First Name or Type :  ',a6,
     &           5x,'Second Name or Type :  ',a6)
      write (iout,250)
  250 format (/,5x,'Bin',9x,'Counts',7x,'Distance',7x,'Raw g(r)',
     &           4x,'Smooth g(r)',/)
      do i = 1, nbin
         write (iout,260)  i,hist(i),(dble(i)-0.5d0)*width,gr(i),gs(i)
  260    format (i8,i15,3x,f12.4,3x,f12.4,3x,f12.4)
      end do
c
c     perform deallocation of some local arrays
c
c     deallocate (hist)
c     deallocate (gr)
c     deallocate (gs)
c
c     perform any final tasks before program exit
c
      call final
      end
