c
c
c     ############################################################
c     ##  COPYRIGHT (C) 1995 by Yong Kong & Jay William Ponder  ##
c     ##                  All Rights Reserved                   ##
c     ############################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  program diffuse  --  find liquid self-diffusion constant  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "diffuse" finds the self-diffusion constant for a homogeneous
c     liquid via the Einstein relation from a set of stored molecular
c     dynamics frames; molecular centers of mass are unfolded and mean
c     squared displacements are computed versus time separation
c
c     the estimate for the self-diffusion constant in 10-5 cm**2/sec
c     is printed in the far right column of output and can be checked
c     by plotting mean squared displacements as a function of the time
c     separation; values for very large time separation are inaccurate
c     due to the small amount of data
c
c
      program diffuse
      use sizes
      use atomid
      use atoms
      use bound
      use inform
      use iounit
      use molcul
      use usage
      implicit none
      integer i,j,k,m
      integer nframe,iframe
      integer iarc,freeunit
      integer start,stop
      integer step,skip
      integer list(20)
      integer, allocatable :: ntime(:)
      real*8 xmid,ymid,zmid
      real*8 xold,yold,zold
      real*8 xdiff,ydiff,zdiff
      real*8 xr,yr,zr,weigh
      real*8 tstep,dunits,delta
      real*8 xvalue,yvalue,zvalue
      real*8 rvalue,dvalue,counts
      real*8, allocatable :: xmsd(:)
      real*8, allocatable :: ymsd(:)
      real*8, allocatable :: zmsd(:)
      real*8, allocatable :: xcm(:,:)
      real*8, allocatable :: ycm(:,:)
      real*8, allocatable :: zcm(:,:)
      logical exist,query
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
c     get the time increment between frames in picoseconds
c
      tstep = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=70,end=70)  tstep
   70 continue
      if (tstep .le. 0.0d0) then
         write (iout,80)
   80    format (/,' Enter the Time Increment in Picoseconds',
     &              ' [1.0] :  ',$)
         read (input,90)  tstep
   90    format (f20.0)
      end if
      if (tstep .le. 0.0d0)  tstep = 1.0d0
c
c     get the atom parameters, lattice type and molecule count
c
      call field
      call unitcell
      call lattice
      call katom
      call molecule
c
c     find atoms and molecules to be excluded from consideration
c
      call active
      if (nuse .eq. n) then
         do i = 1, 20
            list(i) = 0
         end do
         i = 0
         do while (exist)
            call nextarg (string,exist)
            if (exist) then
               read (string,*,err=100,end=100)  list(i+1)
               i = i + 1
            end if
         end do
  100    continue
         if (i .eq. 0) then
            write (iout,110)
  110       format (/,' Numbers of any Atoms to be Removed :  ',$)
            read (input,120)  record
  120       format (a240)
            read (record,*,err=130,end=130)  (list(i),i=1,20)
  130       continue
         end if
         i = 1
         do while (list(i) .ne. 0)
            list(i) = max(-n,min(n,list(i)))
            if (list(i) .gt. 0) then
               k = list(i)
               if (use(k)) then
                  use(k) = .false.
                  nuse = nuse - 1
               end if
               i = i + 1
            else
               list(i+1) = max(-n,min(n,list(i+1)))
               do k = abs(list(i)), abs(list(i+1))
                  if (use(k)) then
                     use(k) = .false.
                     nuse = nuse - 1
                  end if
               end do
               i = i + 2
            end if
         end do
      end if
c
c     alter the molecule list to include only active molecules
c
      do i = 1, nmol
         do j = imol(1,i), imol(2,i)
            k = kmol(j)
            if (.not. use(k))  imol(1,i) = 0
         end do
      end do
      k = 0
      do i = 1, nmol
         if (imol(1,i) .ne. 0) then
            k = k + 1
            imol(1,k) = imol(1,i)
            imol(2,k) = imol(2,i)
            molmass(k) = molmass(i)
         end if
      end do
      nmol = k
      write (iout,140)  nmol
  140 format (/,' Total Number of Molecules :',i16)
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
      write (iout,150)  nframe
  150 format (/,' Number of Coordinate Frames :',i14)
c
c     perform dynamic allocation of some local arrays
c
      allocate (ntime(nframe))
      allocate (xmsd(nframe))
      allocate (ymsd(nframe))
      allocate (zmsd(nframe))
      allocate (xcm(nmol,nframe))
      allocate (ycm(nmol,nframe))
      allocate (zcm(nmol,nframe))
c
c     get the archived coordinates for each frame in turn
c
      write (iout,160)
  160 format (/,' Reading the Coordinates Archive File :',/)
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
         if (n .eq. 0)  goto 180
         nframe = nframe + 1
         if (mod(nframe,100) .eq. 0) then
            write (iout,170)  nframe
  170       format (4x,'Processing Coordinate Frame',i13)
         end if
c
c     unfold each molecule to get its corrected center of mass
c
         do i = 1, nmol
            xmid = 0.0d0
            ymid = 0.0d0
            zmid = 0.0d0
            do j = imol(1,i), imol(2,i)
               k = kmol(j)
               weigh = mass(k)
               xmid = xmid + x(k)*weigh
               ymid = ymid + y(k)*weigh
               zmid = zmid + z(k)*weigh
            end do
            weigh = molmass(i)
            xmid = xmid / weigh
            ymid = ymid / weigh
            zmid = zmid / weigh
            if (nframe .eq. 1) then
               xold = xmid
               yold = ymid
               zold = zmid
            else
               xold = xcm(i,nframe-1)
               yold = ycm(i,nframe-1)
               zold = zcm(i,nframe-1)
            end if
            xr = xmid - xold
            yr = ymid - yold
            zr = zmid - zold
            if (use_bounds)  call image (xr,yr,zr)
            xcm(i,nframe) = xold + xr
            ycm(i,nframe) = yold + yr
            zcm(i,nframe) = zold + zr
         end do
      end do
  180 continue
      close (unit=iarc)
      if (mod(nframe,100) .ne. 0) then
         write (iout,190)  nframe
  190    format (4x,'Processing Coordinate Frame',i13)
      end if
c
c     increment the squared displacements for each frame pair
c
      do i = 1, nframe
         ntime(i) = 0
         xmsd(i) = 0.0d0
         ymsd(i) = 0.0d0
         zmsd(i) = 0.0d0
      end do
      do i = 1, nframe-1
         do j = i+1, nframe
            m = j - i
            ntime(m) = ntime(m) + 1
            do k = 1, nmol
               xdiff = xcm(k,j) - xcm(k,i)
               ydiff = ycm(k,j) - ycm(k,i)
               zdiff = zcm(k,j) - zcm(k,i)
               xmsd(m) = xmsd(m) + xdiff*xdiff
               ymsd(m) = ymsd(m) + ydiff*ydiff
               zmsd(m) = zmsd(m) + zdiff*zdiff
            end do
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xcm)
      deallocate (ycm)
      deallocate (zcm)
c
c     get mean squared displacements and convert units;
c     conversion is from sq. Ang/ps to 10-5 sq. cm/sec
c
      dunits = 10.0d0
      do i = 1, nframe-1
         counts = dble(nmol) * dble(ntime(i))
         xmsd(i) = xmsd(i) * (dunits/counts)
         ymsd(i) = ymsd(i) * (dunits/counts)
         zmsd(i) = zmsd(i) * (dunits/counts)
      end do
c
c     estimate the diffusion constant via the Einstein relation
c
      write (iout,200)
  200 format (/,' Mean Squared Displacements and Self-Diffusion',
     &           ' Constant :',
     &        //,5x,'Time Gap',6x,'X MSD',7x,'Y MSD',7x,'Z MSD',
     &           7x,'R MSD',4x,'Diff Const',
     &        /,7x,'(ps)',9x,'(/2)',8x,'(/2)',8x,'(/2)',8x,'(/6)',
     &           5x,'(x 10^5)',/)
      do i = 1, nframe-1
         delta = tstep * dble(i)
         xvalue = xmsd(i) / 2.0d0
         yvalue = ymsd(i) / 2.0d0
         zvalue = zmsd(i) / 2.0d0
         rvalue = (xmsd(i) + ymsd(i) + zmsd(i)) / 6.0d0
         dvalue = rvalue / delta
         write (iout,210)  delta,xvalue,yvalue,zvalue,rvalue,dvalue
  210    format (f12.2,4f12.2,f12.4)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (ntime)
      deallocate (xmsd)
      deallocate (ymsd)
      deallocate (zmsd)
c
c     perform any final tasks before program exit
c
      call final
      end
