c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  program correlate  --  time correlation of a property  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "correlate" computes the time correlation function of some
c     user-supplied property from individual snapshot frames taken
c     from a molecular dynamics or other trajectory
c
c
      program correlate
      use sizes
      use ascii
      use atoms
      use files
      use inform
      use iounit
      implicit none
      integer maxsite,maxblock
      parameter (maxsite=1000)
      parameter (maxblock=1000)
      integer i,j,k,m
      integer n1,n2,dt
      integer first,last,mode
      integer start,stop,step
      integer nframe,nblock,maxgap
      integer blksize,blkgap,blkdiff
      integer, allocatable :: t1(:)
      integer, allocatable :: t2(:)
      integer, allocatable :: icorr(:)
      real*8 value,property
      real*8, allocatable :: vcorr(:)
      real*8, allocatable :: x1(:,:)
      real*8, allocatable :: y1(:,:)
      real*8, allocatable :: z1(:,:)
      real*8, allocatable :: x2(:,:)
      real*8, allocatable :: y2(:,:)
      real*8, allocatable :: z2(:,:)
      logical exist,query,normal
      character*240 string
c
c
c     determine the desired type of time correlation function
c
      call initial
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
   20    format (/,' The TINKER Correlation Function Utility Can :',
     &           //,4x,'(1) Find Velocity Autocorrelation Function',
     &           /,4x,'(2) Find Superposition Correlation Function')
         do while (mode.lt.1 .or. mode.gt.2)
            mode = 0
            write (iout,30)
   30       format (/,' Enter the Number of the Desired Choice :  ',$)
            read (input,40,err=50,end=50)  mode
   40       format (i10)
   50       continue
         end do
      end if
c
c     get the base name of user specified input structures
c
      call nextarg (filename,exist)
      if (.not. exist) then
         write (iout,60)
   60    format (/,' Enter Base Name of Coordinate Cycle Files :  ',$)
         read (input,70)  filename
   70    format (a240)
      end if
      call basefile (filename)
c
c     set first and last snapshot frames and step increment
c
      first = 0
      last = 0
      step = 0
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=80,end=80)  first
         query = .false.
      end if
      call nextarg (string,exist)
      if (exist)  read (string,*,err=80,end=80)  last
      call nextarg (string,exist)
      if (exist)  read (string,*,err=80,end=80)  step
   80 continue
      if (query) then
         write (iout,90)
   90    format (/,' Numbers of First & Last File and Step',
     &              ' Increment :  ',$)
         read (input,100)  string
  100    format (a240)
         read (string,*,err=110,end=110)  first,last,step
  110    continue
      end if
      if (last .eq. 0)  last = first
      if (step .eq. 0)  step = 1
c
c     set the maximum frame separation to be used for correlation
c
      maxgap = last - first
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=120,end=120)  maxgap
         query = .false.
      end if
  120 continue
      if (query) then
         write (iout,130)
  130    format (/,' Maximum Frame Separation to be Used in',
     &              ' Correlation [ALL] :  ',$)
         read (input,140)  string
  140    format (a240)
         read (string,*,err=150,end=150)  maxgap
  150    continue
      end if
      if (maxgap .eq. 0)  maxgap = last - first
c
c     get the number of frame blocks from the total frames
c
      nframe = 1 + (last-first)/step
      nblock = 1 + (nframe-1)/maxblock
      blksize = maxblock * step
      blkgap = 1 + (maxgap-1)/blksize
      write (iout,160)  nblock,min(nframe,maxblock)
  160 format (/,' Correlation Function Computed using',i5,
     &           ' Blocks of',i6,' Frames')
c
c     perform dynamic allocation of some local arrays
c
      allocate (t1(maxblock))
      allocate (t2(maxblock))
      allocate (icorr(0:maxgap))
      allocate (vcorr(0:maxgap))
      allocate (x1(maxsite,maxblock))
      allocate (y1(maxsite,maxblock))
      allocate (z1(maxsite,maxblock))
      allocate (x2(maxsite,maxblock))
      allocate (y2(maxsite,maxblock))
      allocate (z2(maxsite,maxblock))
c
c     zero out the time correlation function cumulative values
c
      do i = 0, maxgap
         icorr(i) = 0
         vcorr(i) = 0.0d0
      end do
c
c     cycle over all pairs of snapshot frame blocks
c
      do i = 1, nblock
         start = first + (i-1) * blksize
         stop = start + blksize - step
         stop = min(last,stop)
         call readblk (mode,start,stop,step,n1,t1,x1,y1,z1)
         write (iout,170)  i
  170    format (/,3x,'Correlation within Frame Block :    ',i8)
c
c     compute time correlation for frames within single block
c
         do k = 1, n1
            do m = k, n1
               dt = t1(m) - t1(k)
               if (dt .le. maxgap) then
                  value = property (mode,k,x1,y1,z1,m,x1,y1,z1)
                  icorr(dt) = icorr(dt) + 1
                  vcorr(dt) = vcorr(dt) + value
               end if
            end do
         end do
c
c     compute time correlation for frames between two blocks
c
         do j = i+1, min(i+blkgap,nblock)
            start = first + (j-1) * blksize
            stop = start + blksize - step
            stop = min(last,stop)
            blkdiff = (j-i) * maxblock
            call readblk (mode,start,stop,step,n2,t2,x2,y2,z2)
            write (iout,180)  i,j
  180       format (3x,'Correlation between Frame Blocks :  ',2i8)
            do k = 1, n1
               do m = 1, n2
                  dt = t2(m) - t1(k) + blkdiff
                  if (dt .le. maxgap) then
                     value = property (mode,k,x1,y1,z1,m,x2,y2,z2)
                     icorr(dt) = icorr(dt) + 1
                     vcorr(dt) = vcorr(dt) + value
                  end if
               end do
            end do
         end do
      end do
c
c     compute the average correlation function values
c
      do i = 0, maxgap
         if (icorr(i) .ne. 0)  vcorr(i) = vcorr(i)/dble(icorr(i))
      end do
      normal = .false.
      if (vcorr(0) .ne. 0.0d0)  normal = .true.
c
c     print the final values of the correlation function
c
      if (normal) then
         write (iout,190)
  190    format (/,3x,'Separation',7x,'Samples',8x,'Average Value',
     &              7x,'Normalized',/)
         do i = 0, maxgap
            if (icorr(i) .ne. 0) then
               write (iout,200)  i*step,icorr(i),vcorr(i),
     &                           vcorr(i)/vcorr(0)
  200          format (i9,6x,i10,6x,2f17.6)
            end if
         end do
      else
         write (iout,210)
  210    format (/,3x,'Separation',7x,'Samples',8x,'Average Value',/)
         do i = 0, maxgap
            if (icorr(i) .ne. 0) then
               write (iout,220)  i*step,icorr(i),vcorr(i)
  220          format (i9,6x,i10,6x,f17.6)
            end if
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (t1)
      deallocate (t2)
      deallocate (icorr)
      deallocate (vcorr)
      deallocate (x1)
      deallocate (y1)
      deallocate (z1)
      deallocate (x2)
      deallocate (y2)
      deallocate (z2)
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine readblk  --  read a block of snapshot frames  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "readblk" reads in a set of snapshot frames and transfers
c     the values to internal arrays for use in the computation
c     of time correlation functions
c
c
      subroutine readblk (mode,start,stop,step,nb,tb,xb,yb,zb)
      use sizes
      use atomid
      use atoms
      use files
      use iounit
      implicit none
      integer maxsite
      parameter (maxsite=1000)
      integer i,k,ixyz
      integer nt,nb,mode
      integer start,stop
      integer next,label
      integer step,lext
      integer freeunit
      integer tb(*)
      real*8 xb(maxsite,*)
      real*8 yb(maxsite,*)
      real*8 zb(maxsite,*)
      logical exist
      character*7 ext
      character*240 record
      character*240 string
      character*240 xyzfile
c
c
c     initialize the number of files and the numeral size
c
      nt = 0
      nb = 0
      lext = 3
c
c     cycle over all snapshot frames in the block of files
c
      do i = start, stop, step
         nt = nt + 1
         call numeral (i,ext,lext)
         if (mode .eq. 1) then
            xyzfile = filename(1:leng)//'.'//ext(1:lext)//'v'
         else if (mode .eq. 2) then
            xyzfile = filename(1:leng)//'.'//ext(1:lext)
         end if
         inquire (file=xyzfile,exist=exist)
c
c     add file to the current block and get number of atoms
c
         if (exist) then
            nb = nb + 1
            tb(nb) = nt
            ixyz = freeunit ()
            open (unit=ixyz,file=xyzfile,status='old')
            read (ixyz,10)  record
   10       format (a240)
            read (record,*)  n
c
c     check for too many correlation sites in the frame
c
            if (n .gt. maxsite) then
               write (iout,20)
   20          format (/,' READBLK  --  Too many Correlation Sites;',
     &                    ' Increase MAXSITE')
               call fatal
            end if
c
c     read the frame in the TINKER-generated coordinate format;
c     this is fast, but assumes the fixed format shown below
c
c           do k = 1, n
c              read (ixyz,30)  name(k),xb(k,nb),yb(k,nb),zb(k,nb)
c  30          format (8x,a3,3f12.6)
c           end do
c
c     alternatively, get each frame from a free formated file;
c     this is slow, but correctly handles any valid TINKER file
c
            do k = 1, n
               next = 1
               read (ixyz,30)  record
   30          format (a240)
               read (record,*)  label
               call getword (record,name(k),next)
               string = record(next:240)
               read (string,*)  xb(k,nb),yb(k,nb),zb(k,nb)
            end do
            close (unit=ixyz)
         end if
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  function property  --  compute correlation property value  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "property" takes two input snapshot frames and computes the
c     value of the property for which the correlation function is
c     being accumulated
c
c     this version of "property" finds the velocity autocorrelation
c     or the rms fit as a function of time, and is merely provided
c     as an example; the user will need to write a similar custom
c     function to compute other properties to be correlated
c
c
      function property (mode,i,xi,yi,zi,k,xk,yk,zk)
      use sizes
      use atoms
      implicit none
      integer maxsite
      parameter (maxsite=1000)
      integer i,j,k,mode
      real*8 property,value
      real*8, allocatable :: x1(:)
      real*8, allocatable :: y1(:)
      real*8, allocatable :: z1(:)
      real*8, allocatable :: x2(:)
      real*8, allocatable :: y2(:)
      real*8, allocatable :: z2(:)
      real*8 xi(maxsite,*)
      real*8 yi(maxsite,*)
      real*8 zi(maxsite,*)
      real*8 xk(maxsite,*)
      real*8 yk(maxsite,*)
      real*8 zk(maxsite,*)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (x1(maxsite))
      allocate (y1(maxsite))
      allocate (z1(maxsite))
      allocate (x2(maxsite))
      allocate (y2(maxsite))
      allocate (z2(maxsite))
c
c     transfer the input trajectory frames to local vectors
c
      value = 0.0d0
      do j = 1, n
         x1(j) = xi(j,i)
         y1(j) = yi(j,i)
         z1(j) = zi(j,i)
         x2(j) = xk(j,k)
         y2(j) = yk(j,k)
         z2(j) = zk(j,k)
      end do
c
c     sample code to find the velocity autocorrelation function
c
      if (mode .eq. 1) then
         do j = 1, n
            value = value + x1(j)*x2(j) + y1(j)*y2(j) + z1(j)*z2(j)
         end do
      end if
c
c     sample code to find the rms deviation upon superposition
c
      if (mode .eq. 2) then
         call impose (n,x1,y1,z1,n,x2,y2,z2,value)
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (x1)
      deallocate (y1)
      deallocate (z1)
      deallocate (x2)
      deallocate (y2)
      deallocate (z2)
c
c     set property value to be returned for this frame pair
c
      property = value
      return
      end
