c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  program archive  --  create or extract from an archive  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "archive" is a utility program for coordinate files which
c     concatenates multiple coordinate sets into a new archive or
c     performs any of several manipulations on an existing archive
c
c
      program archive
      use sizes
      use atoms
      use bound
      use files
      use inform
      use iounit
      use usage
      implicit none
      integer i,j,k
      integer iarc,ixyz
      integer start,stop
      integer step,now,next
      integer lext,lengb
      integer nmode,mode
      integer leng1,leng2
      integer freeunit
      integer list(20)
      real*8 xr,yr,zr
      real*8, allocatable :: xold(:)
      real*8, allocatable :: yold(:)
      real*8, allocatable :: zold(:)
      logical exist,query
      character*1 answer
      character*7 ext,modtyp
      character*240 arcfile
      character*240 basename
      character*240 xyzfile
      character*240 record
      character*240 string
c
c
c     present a list of possible archive modifications
c
      call initial
      write (iout,10)
   10 format (/,' The TINKER Archive File Utility Can :',
     &        //,4x,'(1) Create an Archive from Individual Frames',
     &        /,4x,'(2) Extract Individual Frames from an Archive',
     &        /,4x,'(3) Trim an Archive to Remove Atoms or Frames',
     &        /,4x,'(4) Unfold Periodic Boundaries for a Trajectory',
     &        /,4x,'(5) Remove Periodic Box Size from a Trajectory')
c
c     get the desired type of archive file modification
c
      next = 1
      nmode = 5
      mode = 0
      call nextarg (answer,exist)
      call getnumb (answer,mode,next)
   20 continue
      if (mode .eq. 0)  mode = -1
      do while (mode.lt.0 .or. mode.gt.nmode)
         mode = 0
         write (iout,30)
   30    format (/,' Number of the Desired Choice [<CR>=Exit] :  ',$)
         read (input,40,err=20,end=50)  mode
   40    format (i10)
   50    continue
      end do
      if (mode .eq. 0)  modtyp = 'EXIT'
      if (mode .eq. 1)  modtyp = 'CREATE'
      if (mode .eq. 2)  modtyp = 'EXTRACT'
      if (mode .eq. 3)  modtyp = 'TRIM'
      if (mode .eq. 4)  modtyp = 'UNFOLD'
      if (mode .eq. 5)  modtyp = 'UNBOUND'
c
c     get the name to use for the coordinate archive file
c
      call nextarg (arcfile,exist)
      if (.not. exist) then
         if (mode .eq. 1) then
            write (iout,60)
   60       format (/,' Enter Base Name of the Individual',
     &                 ' Frames :  ',$)
         else
            write (iout,70)
   70       format (/,' Enter Name of the Coordinate Archive',
     &                 ' File :  ',$)
         end if
         read (input,80)  arcfile
   80    format (a240)
      end if
c
c     create a new archive file or open an existing one
c
      iarc = freeunit ()
      call basefile (arcfile)
      basename = arcfile
      lengb = leng
      if (modtyp .eq. 'CREATE') then
         call suffix (arcfile,'arc','new')
         open (unit=iarc,file=arcfile,status='new')
      else if (modtyp .ne. 'EXIT') then
         call suffix (arcfile,'arc','old')
         inquire (file=arcfile,exist=exist)
         do while (.not. exist)
            write (iout,90)
   90       format (/,' Enter Name of the Coordinate Archive',
     &                 ' File :  ',$)
            read (input,100)  arcfile
  100       format (a240)
            call basefile (arcfile)
            basename = arcfile
            lengb = leng
            call suffix (arcfile,'arc','old')
            inquire (file=arcfile,exist=exist)
         end do
         open (unit=iarc,file=arcfile,status='old')
         rewind (unit=iarc)
         call readxyz (iarc)
         rewind (unit=iarc)
         call active
      end if
c
c     combine individual files into a single archive file
c
      if (modtyp .eq. 'CREATE') then
         modtyp = 'EXIT'
         start = 0
         stop = 0
         step = 0
         query = .true.
         call nextarg (string,exist)
         if (exist) then
            read (string,*,err=110,end=110)  start
            query = .false.
         end if
         call nextarg (string,exist)
         if (exist)  read (string,*,err=110,end=110)  stop
         call nextarg (string,exist)
         if (exist)  read (string,*,err=110,end=110)  step
  110    continue
         if (query) then
            write (iout,120)
  120       format (/,' Numbers of First & Last File and Step',
     &                 ' Increment :  ',$)
            read (input,130)  record
  130       format (a240)
            read (record,*,err=140,end=140)  start,stop,step
  140       continue
         end if
         if (stop .eq. 0)  stop = start
         if (step .eq. 0)  step = 1
c
c     cycle over the user specified coordinate files
c
         i = start
         do while (i.ge.start .and. i.le.stop)
            ixyz = freeunit ()
            lext = 3
            call numeral (i,ext,lext)
            xyzfile = basename(1:lengb)//'.'//ext(1:lext)
            call version (xyzfile,'old')
            inquire (file=xyzfile,exist=exist)
            if (.not.exist .and. i.lt.100) then
               lext = 2
               call numeral (i,ext,lext)
               xyzfile = basename(1:lengb)//'.'//ext(1:lext)
               call version (xyzfile,'old')
               inquire (file=xyzfile,exist=exist)
            end if
            if (.not.exist .and. i.lt.10) then
               lext = 1
               call numeral (i,ext,lext)
               xyzfile = basename(1:lengb)//'.'//ext(1:lext)
               call version (xyzfile,'old')
               inquire (file=xyzfile,exist=exist)
            end if
            if (exist) then
               open (unit=ixyz,file=xyzfile,status='old')
               rewind (unit=ixyz)
               call readxyz (ixyz)
               close (unit=ixyz)
               if (i .eq. start)  call active
               nuse = n
               do j = 1, n
                  use(j) = .true.
               end do
               call prtarc (iarc)
            end if
            i = i + step
         end do
      end if
c
c     decide whether atoms are to be removed from each frame
c
      if (modtyp .eq. 'TRIM') then
         call active
         if (nuse .eq. n) then
            do i = 1, 20
               list(i) = 0
            end do
            write (iout,150)
  150       format (/,' Numbers of the Atoms to be Removed :  ',$)
            read (input,160)  record
  160       format (a240)
            read (record,*,err=170,end=170)  (list(i),i=1,20)
  170       continue
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
      end if
c
c     store index to use in renumbering the untrimmed atoms
c
      k = 0
      do i = 1, n
         iuse(i) = 0
         if (use(i)) then
            k = k + 1
            iuse(i) = k
         end if
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (xold(n))
      allocate (yold(n))
      allocate (zold(n))
c
c     get the initial and final coordinate frames to process
c
      if (modtyp .ne. 'EXIT') then
         now = 1
         leng1 = 1
         leng2 = leng
         do i = 1, leng
            if (filename(i:i) .eq. '/')  leng1 = i+1
            if (filename(i:i) .eq. ']')  leng1 = i+1
            if (filename(i:i) .eq. ':')  leng1 = i+1
         end do
         do i = leng, leng1, -1
            if (filename(i:i) .eq. '.')  leng2 = i-1
         end do
         leng = leng2 - leng1 + 1
         filename(1:leng) = filename(leng1:leng2)
         start = 0
         stop = 0
         step = 0
         query = .true.
         call nextarg (string,exist)
         if (exist) then
            read (string,*,err=180,end=180)  start
            query = .false.
         end if
         call nextarg (string,exist)
         if (exist)  read (string,*,err=180,end=180)  stop
         call nextarg (string,exist)
         if (exist)  read (string,*,err=180,end=180)  step
  180    continue
         if (query) then
            write (iout,190)
  190       format (/,' Numbers of First & Last File and Step',
     &                 ' [<CR>=Exit] :  ',$)
            read (input,200)  record
  200       format (a240)
            read (record,*,err=210,end=210)  start,stop,step
  210       continue
         end if
         if (stop .eq. 0)  stop = start
         if (step .eq. 0)  step = 1
c
c     loop over the individual coordinate files to be extracted
c
         do while (start .ne. 0)
            if (start .le. now) then
               now = 1
               rewind (unit=iarc)
            end if
            do k = 1, start-now
               call readxyz (iarc)
            end do
            i = start
            if (modtyp .eq. 'EXTRACT') then
               do while (i.ge.start .and. i.le.stop)
                  lext = 3
                  call numeral (i,ext,lext)
                  call readxyz (iarc)
                  if (abort)  goto 220
                  nuse = n
                  do j = 1, n
                     use(j) = .true.
                  end do
                  ixyz = freeunit ()
                  xyzfile = filename(1:leng)//'.'//ext(1:lext)
                  call version (xyzfile,'new')
                  open (unit=ixyz,file=xyzfile,status='new')
                  call prtarc (ixyz)
                  close (unit=ixyz)
                  i = i + step
                  do k = 1, step-1
                     call readxyz (iarc)
                  end do
               end do
            else
               ixyz = freeunit ()
               xyzfile = basename
               call suffix (xyzfile,'arc','new')
               open (unit=ixyz,file=xyzfile,status='new')
               do while (i.ge.start .and. i.le.stop)
                  call readxyz (iarc)
                  if (abort)  goto 220
                  if (modtyp .eq. 'UNFOLD') then
                     nuse = n
                     do j = 1, n
                        use(j) = .true.
                     end do
                     if (i .eq. start) then
                        call unitcell
                        do j = 1, n
                           xold(j) = x(j)
                           yold(j) = y(j)
                           zold(j) = z(j)
                        end do
                     end if
                     call lattice
                     do j = 1, n
                        xr = x(j) - xold(j)
                        yr = y(j) - yold(j)
                        zr = z(j) - zold(j)
                        if (use_bounds)  call image (xr,yr,zr)
                        x(j) = xold(j) + xr
                        y(j) = yold(j) + yr
                        z(j) = zold(j) + zr
                        xold(j) = x(j)
                        yold(j) = y(j)
                        zold(j) = z(j)
                     end do
                  else if (modtyp .eq. 'UNBOUND') then
                     use_bounds = .false.
                  end if
                  call prtarc (ixyz)
                  i = i + step
                  do k = 1, step-1
                     call readxyz (iarc)
                  end do
               end do
               close (unit=ixyz)
            end if
  220       continue
            now = stop
            start = 0
            stop = 0
            step = 0
            query = .true.
            call nextarg (string,exist)
            if (exist) then
               read (string,*,err=230,end=230)  start
               query = .false.
            end if
            call nextarg (string,exist)
            if (exist)  read (string,*,err=230,end=230)  stop
            call nextarg (string,exist)
            if (exist)  read (string,*,err=230,end=230)  step
  230       continue
            if (query) then
               write (iout,240)
  240          format (/,' Numbers of First & Last File and Step',
     &                    ' [<CR>=Exit] :  ',$)
               read (input,250)  record
  250          format (a240)
               read (record,*,err=260,end=260)  start,stop,step
  260          continue
            end if
            if (stop .eq. 0)  stop = start
            if (step .eq. 0)  step = 1
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (xold)
      deallocate (yold)
      deallocate (zold)
c
c     perform any final tasks before program exit
c
      close (unit=iarc)
      call final
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine prtarc  --  output of a TINKER archive file  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "prtarc" writes out a set of Cartesian coordinates for
c     all active atoms in the TINKER XYZ archive format
c
c
      subroutine prtarc (iarc)
      use sizes
      use atomid
      use atoms
      use bound
      use boxes
      use couple
      use files
      use inform
      use titles
      use usage
      implicit none
      integer i,k,iarc
      integer size,crdsiz
      real*8 crdmin,crdmax
      logical opened
      character*2 atmc
      character*2 crdc
      character*2 digc
      character*25 fstr
      character*240 arcfile
c
c
c     open output unit if not already done
c
      inquire (unit=iarc,opened=opened)
      if (.not. opened) then
         arcfile = filename(1:leng)//'.arc'
         call version (arcfile,'new')
         open (unit=iarc,file=arcfile,status='new')
      end if
c
c     check for large systems needing extended formatting
c
      atmc = 'i6'
      if (n .ge. 100000)  atmc = 'i7'
      if (n .ge. 1000000)  atmc = 'i8'
      crdmin = 0.0d0
      crdmax = 0.0d0
      do i = 1, n
         crdmin = min(crdmin,x(i),y(i),z(i))
         crdmax = max(crdmax,x(i),y(i),z(i))
      end do
      crdsiz = 6
      if (crdmin .le. -1000.0d0)  crdsiz = 7
      if (crdmax .ge. 10000.0d0)  crdsiz = 7
      if (crdmin .le. -10000.0d0)  crdsiz = 8
      if (crdmax .ge. 100000.0d0)  crdsiz = 8
      crdsiz = crdsiz + max(6,digits)
      size = 0
      call numeral (crdsiz,crdc,size)
      if (digits .le. 6) then
         digc = '6 '
      else if (digits .le. 8) then
         digc = '8'
      else
         digc = '10'
      end if
c
c     write out the number of atoms and the title
c
      if (ltitle .eq. 0) then
         fstr = '('//atmc//')'
         write (iarc,fstr(1:4))  nuse
      else
         fstr = '('//atmc//',2x,a)'
         write (iarc,fstr(1:9))  nuse,title(1:ltitle)
      end if
c
c     write out the periodic cell lengths and angles
c
      if (use_bounds) then
         fstr = '(1x,6f'//crdc//'.'//digc//')'
         write (iarc,fstr)  xbox,ybox,zbox,alpha,beta,gamma
      end if
c
c     write out the coordinate line for each atom
c
      fstr = '('//atmc//',2x,a3,3f'//crdc//
     &          '.'//digc//',i6,8'//atmc//')'
      do i = 1, n
         if (use(i)) then
            write (iarc,fstr)  iuse(i),name(i),x(i),y(i),z(i),type(i),
     &                         (iuse(i12(k,i)),k=1,n12(i))
         end if
      end do
c
c     close the output unit if opened by this routine
c
      if (.not. opened)  close (unit=iarc)
      return
      end
