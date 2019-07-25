c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2000 by P. Bagossi, P. Ren & Jay W. Ponder  ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  program poledit  --  manipulate atomic multipole values  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "poledit" provides for modification and manipulation of
c     the atomic multipole electrostatic models used in TINKER
c
c
      program poledit
      use iounit
      use potent
      implicit none
      integer mode
      logical exist,query
      character*240 string
c
c
c     get the desired type of coordinate file modification
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
   20    format (/,' The TINKER Multipole Editing Utility Can :',
     &           //,4x,'(1) Multipole Parameters from GDMA Output',
     &           /,4x,'(2) Alter Local Coordinate Frame Definitions',
     &           /,4x,'(3) Removal of Intramolecular Polarization')
         do while (mode.lt.1 .or. mode.gt.3)
            mode = 0
            write (iout,30)
   30       format (/,' Enter the Number of the Desired Choice :  ',$)
            read (input,40,err=50,end=50)  mode
   40       format (i10)
   50       continue
         end do
      end if
c
c     perform the desired multipole manipulation operation
c
      if (mode .eq. 1) then
         use_mpole = .true.
         use_polar = .true.
         call readgdma
         call initprm
         call molsetup
         call setframe
         call rotframe
         call setpolar
         call alterpol
         call fixpolar
         call prtpolar
      else if (mode .eq. 2) then
         call getxyz
         call attach
         call field
         call katom
         call kmpole
         call kpolar
         call fixframe
         call prtpolar
      else if (mode .eq. 3) then
         call getxyz
         call attach
         call field
         call katom
         call kmpole
         call kpolar
         call alterpol
         call fixpolar
         call prtpolar
      end if
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine readgdma  --  get information from GDMA output  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "readgdma" takes the DMA output in spherical harmonics from
c     the GDMA program and converts to Cartesian multipoles in
c     the global coordinate frame
c
c     this version is compatible with the formatted output from
c     GDMA 2.2.04 released by Anthony Stone in Fall 2008
c
c
      subroutine readgdma
      use sizes
      use atomid
      use atoms
      use dma
      use iounit
      use mpole
      use units
      implicit none
      integer i,j,k
      integer idma,next
      integer freeunit
      real*8 term
      logical exist,done
      logical use_bohr
      character*3 atmnam
      character*240 record
      character*240 dmafile
c
c
c     perform dynamic allocation of some global arrays
c
      if (.not. allocated(mp))  allocate (mp(maxatm))
      if (.not. allocated(dpx))  allocate (dpx(maxatm))
      if (.not. allocated(dpy))  allocate (dpy(maxatm))
      if (.not. allocated(dpz))  allocate (dpz(maxatm))
      if (.not. allocated(q20))  allocate (q20(maxatm))
      if (.not. allocated(q21c))  allocate (q21c(maxatm))
      if (.not. allocated(q21s))  allocate (q21s(maxatm))
      if (.not. allocated(q22c))  allocate (q22c(maxatm))
      if (.not. allocated(q22s))  allocate (q22s(maxatm))
c
c     zero out the atomic coordinates and DMA values
c
      do i = 1, maxatm
         x(i) = 0.0d0
         y(i) = 0.0d0
         z(i) = 0.0d0
         mp(i) = 0.0d0
         dpx(i) = 0.0d0
         dpy(i) = 0.0d0
         dpz(i) = 0.0d0
         q20(i) = 0.0d0
         q21c(i) = 0.0d0
         q21s(i) = 0.0d0
         q22c(i) = 0.0d0
         q22s(i) = 0.0d0
      end do
c
c     try to get a filename from the command line arguments
c
      call nextarg (dmafile,exist)
      if (exist) then
         call basefile (dmafile)
         call suffix (dmafile,'dma','old')
         inquire (file=dmafile,exist=exist)
      end if
c
c     ask for the user specified GDMA output filename
c
      do while (.not. exist)
         write (iout,10)
   10    format (/,' Enter GDMA Output File Name :  ',$)
         read (input,20)  dmafile
   20    format (a240)
         call basefile (dmafile)
         call suffix (dmafile,'dma','old')
         inquire (file=dmafile,exist=exist)
      end do
c
c     first open and then read the GDMA output file
c
      idma = freeunit ()
      open (unit=idma,file=dmafile,status='old')
c
c     get coordinates and multipoles from GDMA output file
c
      i = 0
      rewind (unit=idma)
      do while (.true.)
         read (idma,30,err=50,end=50)  record
   30    format (a240)
         if (i .ne. 0)  call match1 (i,record)
         if (record(12:14) .eq. 'x =') then
            i = i + 1
            next = 1
            call gettext (record,name(i),next)
            read (record(15:24),*)  x(i)
            read (record(30:39),*)  y(i)
            read (record(45:54),*)  z(i)
            read (idma,40,err=50,end=50)
   40       format ()
         else if (record(1:16) .eq. 'Total multipoles') then
            goto 50
         end if
      end do
   50 continue
      n = i
c
c     perform dynamic allocation of some global arrays
c
      if (.not. allocated(rpole))  allocate (rpole(maxpole,n))
c
c     convert quadrupole from spherical harmonic to Cartesian
c
      term = sqrt(0.75d0)
      do i = 1, n
         rpole(1,i) = mp(i)
         rpole(2,i) = dpx(i)
         rpole(3,i) = dpy(i)
         rpole(4,i) = dpz(i)
         rpole(5,i) = -0.5d0*q20(i) + term*q22c(i)
         rpole(6,i) = term*q22s(i)
         rpole(7,i) = term*q21c(i)
         rpole(8,i) = rpole(6,i)
         rpole(9,i) = -0.5d0*q20(i) - term*q22c(i)
         rpole(10,i) = term*q21s(i)
         rpole(11,i) = rpole(7,i)
         rpole(12,i) = rpole(10,i)
         rpole(13,i) = q20(i)
      end do
c
c     check for GDMA coordinate values in atomic units
c
      use_bohr = .false.
      rewind (unit=idma)
      do while (.true.)
         read (idma,60,err=70,end=70)  record
   60    format (a240)
         if (record(1:27) .eq. 'Positions and radii in bohr') then
            use_bohr = .true.
            goto 70
         end if
      end do
   70 continue
c
c     convert coordinates from Bohrs to Angstroms if needed
c
      if (use_bohr) then
         do i = 1, n
            x(i) = x(i) * bohr
            y(i) = y(i) * bohr
            z(i) = z(i) * bohr
         end do
      end if
c
c     find atomic numbers in verbose GDMA output if available
c
      done = .false.
      rewind (unit=idma)
      do while (.true.)
         read (idma,80,err=100,end=100)  record
   80    format (a240)
         if (record(1:16) .eq. 'Nuclear charges:') then
            k = min(n,20)
            read (record(17:240),*,err=100,end=100)  (atomic(i),i=1,k)
            do while (k .ne. n)
               j = k + 1
               k = min(n,k+20)
               read (idma,90,err=100,end=100)  record
   90          format (a240)
               read (record,*,err=100,end=100)  (atomic(i),i=j,k)
            end do
            done = .true.
         end if
      end do
  100 continue
      close (unit=idma)
c
c     attempt to get atomic numbers from GDMA atom names
c
      if (.not. done) then
         do i = 1, n
            atomic(i) = 0
            atmnam = name(i)
            call upcase (atmnam)
            if (atmnam(1:2) .eq. 'SI') then
               atomic(i) = 14
            else if (atmnam(1:2) .eq. 'CL') then
               atomic(i) = 17
            else if (atmnam(1:2) .eq. 'BR') then
               atomic(i) = 35
            else if (atmnam(1:1) .eq. 'H') then
               atomic(i) = 1
            else if (atmnam(1:1) .eq. 'B') then
               atomic(i) = 5
            else if (atmnam(1:1) .eq. 'C') then
               atomic(i) = 6
            else if (atmnam(1:1) .eq. 'N') then
               atomic(i) = 7
            else if (atmnam(1:1) .eq. 'O') then
               atomic(i) = 8
            else if (atmnam(1:1) .eq. 'F') then
               atomic(i) = 9
            else if (atmnam(1:1) .eq. 'P') then
               atomic(i) = 15
            else if (atmnam(1:1) .eq. 'S') then
               atomic(i) = 16
            else if (atmnam(1:1) .eq. 'I') then
               atomic(i) = 53
            else
               read (atmnam,*,err=110,end=110)  atomic(i)
  110          continue
            end if
         end do
      end if
c
c     print the global frame Cartesian atomic multipoles
c
      write (iout,120)
  120 format (/,' Global Frame Cartesian Multipole Moments :')
      do i = 1, n
         write (iout,130)  i,name(i),atomic(i)
  130    format (/,' Site:',i8,9x,'Name:',3x,a3,7x,'Atomic Number:',i8)
         write (iout,140)  x(i),y(i),z(i)
  140    format (/,' Coordinates:',5x,3f15.6)
         write (iout,150)  rpole(1,i)
  150    format (/,' Charge:',10x,f15.5)
         write (iout,160)  rpole(2,i),rpole(3,i),rpole(4,i)
  160    format (' Dipole:',10x,3f15.5)
         write (iout,170)  rpole(5,i)
  170    format (' Quadrupole:',6x,f15.5)
         write (iout,180)  rpole(8,i),rpole(9,i)
  180    format (18x,2f15.5)
         write (iout,190)  rpole(11,i),rpole(12,i),rpole(13,i)
  190    format (18x,3f15.5)
      end do
c
c     convert the dipole and quadrupole moments to Angstroms,
c     quadrupole divided by 3 for use as traceless values
c
      do i = 1, n
         do k = 2, 4
            rpole(k,i) = rpole(k,i) * bohr
         end do
         do k = 5, 13
            rpole(k,i) = rpole(k,i) * bohr**2 / 3.0d0
         end do
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine match1  --  match first value from GDMA output  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "match1" finds and stores the first multipole component found
c     on a line of output from Stone's GDMA program
c
c
      subroutine match1 (i,record)
      use sizes
      use dma
      implicit none
      integer i
      character*240 record
c
c
c     store first multipole components on a line of GDMA output
c
      if (record(6:8) .eq. 'Q0 ') then
         read (record(13:23),*)  mp(i)
         call match2 (i,record)
      else if (record(20:23) .eq. 'Q00 ') then
         read (record(26:36),*)  mp(i)
      else if (record(20:23) .eq. 'Q10 ') then
         read (record(26:36),*)  dpz(i)
         call match2 (i,record)
      else if (record(20:23) .eq. 'Q11c') then
         read (record(26:36),*)  dpx(i)
         call match2 (i,record)
      else if (record(20:23) .eq. 'Q11s') then
         read (record(26:36),*)  dpy(i)
         call match2 (i,record)
      else if (record(20:23) .eq. 'Q20 ') then
         read (record(26:36),*)  q20(i)
         call match2 (i,record)
      else if (record(20:23) .eq. 'Q21c') then
         read (record(26:36),*)  q21c(i)
         call match2 (i,record)
      else if (record(20:23) .eq. 'Q21s') then
         read (record(26:36),*)  q21s(i)
         call match2 (i,record)
      else if (record(20:23) .eq. 'Q22c') then
         read (record(26:36),*)  q22c(i)
         call match2 (i,record)
      else if (record(20:23) .eq. 'Q22s') then
         read (record(26:36),*)  q22s(i)
         call match2 (i,record)
      end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine match2  --  match second value from GDMA output  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "match2" finds and stores the second multipole component found
c     on a line of output from Stone's GDMA program
c
c
      subroutine match2 (i,record)
      use sizes
      use dma
      implicit none
      integer i
      character*240 record
c
c
c     store second multipole component on a line of GDMA output
c
      if (record(29:31) .eq. 'Q1 ') then
         read (record(36:46),*)  dpz(i)
         call match3 (i,record)
      else if (record(39:42) .eq. 'Q11c') then
         read (record(45:55),*)  dpx(i)
         call match3 (i,record)
      else if (record(39:42) .eq. 'Q11s') then
         read (record(45:55),*)  dpy(i)
         call match3 (i,record)
      else if (record(39:42) .eq. 'Q21c') then
         read (record(45:55),*)  q21c(i)
         call match3 (i,record)
      else if (record(39:42) .eq. 'Q21s') then
         read (record(45:55),*)  q21s(i)
         call match3 (i,record)
      else if (record(39:42) .eq. 'Q22c') then
         read (record(45:55),*)  q22c(i)
         call match3 (i,record)
      else if (record(39:42) .eq. 'Q22s') then
         read (record(45:55),*)  q22s(i)
         call match3 (i,record)
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine match3  --  match third value from GDMA output  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "match3" finds and stores the third multipole component found
c     on a line of output from Stone's GDMA program
c
c
      subroutine match3 (i,record)
      use sizes
      use dma
      implicit none
      integer i
      character*240 record
c
c
c     store third multipole component on a line of GDMA output
c
      if (record(52:54) .eq. 'Q2 ') then
         read (record(59:69),*)  q20(i)
      else if (record(58:61) .eq. 'Q11s') then
         read (record(64:74),*)  dpy(i)
      else if (record(58:61) .eq. 'Q21s') then
         read (record(64:74),*)  q21s(i)
      else if (record(58:61) .eq. 'Q22c') then
         read (record(64:74),*)  q22c(i)
      else if (record(58:61) .eq. 'Q22s') then
         read (record(64:74),*)  q22s(i)
      end if
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine molsetup  --  set molecule for polarization  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "molsetup" generates trial parameters needed to perform
c     polarizable multipole calculations on a structure read
c     from a GDMA output file
c
c
      subroutine molsetup
      use sizes
      use atomid
      use atoms
      use couple
      use files
      use mpole
      use polar
      implicit none
      integer i,j,k,m,ixyz
      integer atmnum,size
      integer freeunit
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 ri,rij,dij
      real*8 sixth
      real*8, allocatable :: rad(:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (rad(n))
c
c     set base atomic radii from covalent radius values
c
      do i = 1, n
         rad(i) = 0.77d0
         atmnum = atomic(i)
         if (atmnum .eq. 0)  rad(i) = 0.00d0
         if (atmnum .eq. 1)  rad(i) = 0.37d0
         if (atmnum .eq. 2)  rad(i) = 0.32d0
         if (atmnum .eq. 6)  rad(i) = 0.77d0
         if (atmnum .eq. 7)  rad(i) = 0.75d0
         if (atmnum .eq. 8)  rad(i) = 0.73d0
         if (atmnum .eq. 9)  rad(i) = 0.71d0
         if (atmnum .eq. 10)  rad(i) = 0.69d0
         if (atmnum .eq. 14)  rad(i) = 1.11d0
         if (atmnum .eq. 15)  rad(i) = 1.06d0
         if (atmnum .eq. 16)  rad(i) = 1.02d0
         if (atmnum .eq. 17)  rad(i) = 0.99d0
         if (atmnum .eq. 18)  rad(i) = 0.97d0
         if (atmnum .eq. 35)  rad(i) = 1.14d0
         if (atmnum .eq. 36)  rad(i) = 1.10d0
         if (atmnum .eq. 53)  rad(i) = 1.33d0
         if (atmnum .eq. 54)  rad(i) = 1.30d0
         rad(i) = 1.1d0 * rad(i)
      end do
c
c     assign atom connectivities based on interatomic distances
c
      do i = 1, n
         n12(i) = 0
      end do
      do i = 1, n-1
         xi = x(i)
         yi = y(i)
         zi = z(i)
         ri = rad(i)
         do j = i+1, n
            xr = x(j) - xi
            yr = y(j) - yi
            zr = z(j) - zi
            rij = ri + rad(j)
            dij = sqrt(xr*xr + yr*yr + zr*zr)
            if (dij .lt. rij) then
               n12(i) = n12(i) + 1
               i12(n12(i),i) = j
               n12(j) = n12(j) + 1
               i12(n12(j),j) = i
            end if
         end do
      end do
      do i = 1, n
         call sort (n12(i),i12(1,i))
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (rad)
c
c     assign unique atom types and set the valence values
c
      size = min(20,leng)
      do i = 1, n
         type(i) = i
         valence(i) = n12(i)
         story(i) = filename(1:size)
      end do
c
c     create a file with coordinates and connectivities
c
      ixyz = freeunit ()
      call prtxyz (ixyz)
c
c     perform dynamic allocation of some global arrays
c
      if (.not. allocated(polarity))  allocate (polarity(n))
      if (.not. allocated(thole))  allocate (thole(n))
      if (.not. allocated(pdamp))  allocate (pdamp(n))
      if (.not. allocated(ipole))  allocate (ipole(n))
      if (.not. allocated(polsiz))  allocate (polsiz(n))
      if (.not. allocated(pollist))  allocate (pollist(n))
c
c     assign atomic mass and polarizability by atomic number
c
      do i = 1, n
         atmnum = atomic(i)
         mass(i) = 1.0d0
         polarity(i) = 0.0d0
         thole(i) = 0.39d0
         if (atmnum .eq. 1) then
            mass(i) = 1.008d0
            polarity(i) = 0.496d0
         else if (atmnum .eq. 5) then
            mass(i) = 10.810d0
            polarity(i) = 1.600d0
         else if (atmnum .eq. 6) then
            mass(i) = 12.011d0
            polarity(i) = 1.334d0
         else if (atmnum .eq. 7) then
            mass(i) = 14.007d0
            polarity(i) = 1.073d0
         else if (atmnum .eq. 8) then
            mass(i) = 15.999d0
            polarity(i) = 0.837d0
         else if (atmnum .eq. 9) then
            mass(i) = 18.998d0
            polarity(i) = 0.507d0
         else if (atmnum .eq. 14) then
            mass(i) = 28.086d0
         else if (atmnum .eq. 15) then
            mass(i) = 30.974d0
            polarity(i) = 1.828d0
         else if (atmnum .eq. 16) then
            mass(i) = 32.066d0
            polarity(i) = 3.300d0
         else if (atmnum .eq. 17) then
            mass(i) = 35.453d0
            polarity(i) = 2.500d0
         else if (atmnum .eq. 35) then
            mass(i) = 79.904d0
            polarity(i) = 3.595d0
         else if (atmnum .eq. 53) then
            mass(i) = 126.904d0
            polarity(i) = 5.705d0
         end if
      end do
c
c     alter polarizabilities for aromatic carbon and hydrogen
c
      do i = 1, n
         atmnum = atomic(i)
         if (atmnum .eq. 1) then
            j = i12(1,i)
            if (atomic(j).eq.6 .and. n12(j).eq.3) then
               polarity(i) = 0.696d0
               do k = 1, n12(j)
                  m = i12(k,j)
                  if (atomic(m).eq.8 .and. n12(m).eq.1) then
                     polarity(i) = 0.494d0
                  end if
               end do
            end if
         else if (atmnum .eq. 6) then
            if (n12(i) .eq. 3) then
               polarity(i) = 1.75d0
               do j = 1, n12(i)
                  k = i12(j,i)
                  if (atomic(k).eq.8 .and. n12(k).eq.1) then
                     polarity(i) = 1.334d0
                  end if
               end do
            end if
         end if
      end do
c
c     set atomic multipole and polarizability scaling values
c
      npole = n
      npolar = n
      sixth = 1.0d0 / 6.0d0
      do i = 1, n
         ipole(i) = i
         polsiz(i) = 13
         pollist(i) = i
         if (thole(i) .eq. 0.0d0) then
            pdamp(i) = 0.0d0
         else
            pdamp(i) = polarity(i)**sixth
         end if
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine setframe  --  define local coordinate frames  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "setframe" assigns a local coordinate frame at each atomic
c     multipole site using high priority connected atoms along axes
c
c
      subroutine setframe
      use sizes
      use atomid
      use atoms
      use couple
      use iounit
      use mpole
      implicit none
      integer i,j,k,m,kb
      integer ia,ib,ic,id
      integer kab,kac,kad
      integer kbc,kbd,kcd
      integer priority
      logical exist,query
      logical change
      character*240 record
      character*240 string
c
c
c     perform dynamic allocation of some global arrays
c
      if (.not. allocated(zaxis))  allocate (zaxis(npole))
      if (.not. allocated(xaxis))  allocate (xaxis(npole))
      if (.not. allocated(yaxis))  allocate (yaxis(npole))
      if (.not. allocated(polaxe))  allocate (polaxe(npole))
c
c     initialize the local frame type and defining atoms
c
      do i = 1, npole
         polaxe(i) = 'None'
         zaxis(i) = 0
         xaxis(i) = 0
         yaxis(i) = 0
      end do
c
c     assign the local frame definition for an isolated atom
c
      do i = 1, npole
         j = n12(i)
         if (j .eq. 0) then
            polaxe(i) = 'None'
            zaxis(i) = 0
            xaxis(i) = 0
            yaxis(i) = 0
c
c     assign the local frame definition for a monovalent atom
c
         else if (j .eq. 1) then
            ia = i12(1,i)
            if (n12(ia) .eq. 1) then
               polaxe(i) = 'Z-Only'
               zaxis(i) = ia
               xaxis(i) = 0
               yaxis(i) = 0
            else
               polaxe(i) = 'Z-then-X'
               zaxis(i) = ia
               m = 0
               do k = 1, n12(ia)
                  ib = i12(k,ia)
                  kb = atomic(ib)
                  if (kb.gt.m .and. ib.ne.i) then
                     xaxis(i) = ib
                     m = kb
                  end if
               end do
               yaxis(i) = 0
            end if
c
c     assign the local frame definition for a divalent atom
c
         else if (j .eq. 2) then
            ia = i12(1,i)
            ib = i12(2,i)
            kab = priority (i,ia,ib)
            if (kab .eq. ia) then
               polaxe(i) = 'Z-then-X'
               zaxis(i) = ia
               xaxis(i) = ib
               yaxis(i) = 0
            else if (kab .eq. ib) then
               polaxe(i) = 'Z-then-X'
               zaxis(i) = ib
               xaxis(i) = ia
               yaxis(i) = 0
            else
               polaxe(i) = 'Bisector'
               zaxis(i) = ia
               xaxis(i) = ib
               yaxis(i) = 0
            end if
c
c     assign the local frame definition for a trivalent atom
c
         else if (j .eq. 3) then
            ia = i12(1,i)
            ib = i12(2,i)
            ic = i12(3,i)
            kab = priority (i,ia,ib)
            kac = priority (i,ia,ic)
            kbc = priority (i,ib,ic)
            if (kab.eq.0 .and. kac.eq.0) then
               polaxe(i) = 'Z-then-X'
               zaxis(i) = ia
               xaxis(i) = ib
               yaxis(i) = 0
            else if (kab.eq.ia .and. kac.eq.ia) then
               polaxe(i) = 'Z-then-X'
               zaxis(i) = ia
               xaxis(i) = ib
               if (kbc .eq. ic)  xaxis(i) = ic
               yaxis(i) = 0
            else if (kab.eq.ib .and. kbc.eq.ib) then
               polaxe(i) = 'Z-then-X'
               zaxis(i) = ib
               xaxis(i) = ia
               if (kac .eq. ic)  xaxis(i) = ic
               yaxis(i) = 0
            else if (kac.eq.ic .and. kbc.eq.ic) then
               polaxe(i) = 'Z-then-X'
               zaxis(i) = ic
               xaxis(i) = ia
               if (kab .eq. ib)  xaxis(i) = ib
               yaxis(i) = 0
            else if (kab .eq. 0) then
               polaxe(i) = 'Bisector'
               zaxis(i) = ia
               xaxis(i) = ib
               yaxis(i) = 0
            else if (kac .eq. 0) then
               polaxe(i) = 'Bisector'
               zaxis(i) = ia
               xaxis(i) = ic
               yaxis(i) = 0
            else if (kbc .eq. 0) then
               polaxe(i) = 'Bisector'
               zaxis(i) = ib
               xaxis(i) = ic
               yaxis(i) = 0
            end if
c
c     assign the local frame definition for a tetravalent atom
c
         else if (j .eq. 4) then
            ia = i12(1,i)
            ib = i12(2,i)
            ic = i12(3,i)
            id = i12(4,i)
            kab = priority (i,ia,ib)
            kac = priority (i,ia,ic)
            kad = priority (i,ia,id)
            kbc = priority (i,ib,ic)
            kbd = priority (i,ib,id)
            kcd = priority (i,ic,id)
            if (kab.eq.0 .and. kac.eq.0 .and. kad.eq.0) then
               polaxe(i) = 'Z-then-X'
               zaxis(i) = ia
               xaxis(i) = ib
               yaxis(i) = 0
            else if (kab.eq.ia .and. kac.eq.ia .and. kad.eq.ia) then
               polaxe(i) = 'Z-then-X'
               zaxis(i) = ia
               xaxis(i) = ib
               if (kbc.eq.ic .and. kcd.eq.ic)  xaxis(i) = ic
               if (kbd.eq.id .and. kcd.eq.id)  xaxis(i) = id
               if (kbc.eq.ic .and. kcd.eq.0)  xaxis(i) = ic
               yaxis(i) = 0
            else if (kab.eq.ib .and. kbc.eq.ib .and. kbd.eq.ib) then
               polaxe(i) = 'Z-then-X'
               zaxis(i) = ib
               xaxis(i) = ia
               if (kac.eq.ic .and. kcd.eq.ic)  xaxis(i) = ic
               if (kad.eq.id .and. kcd.eq.id)  xaxis(i) = id
               if (kac.eq.ic .and. kcd.eq.0)  xaxis(i) = ic
               yaxis(i) = 0
            else if (kac.eq.ic .and. kbc.eq.ic .and. kcd.eq.ic) then
               polaxe(i) = 'Z-then-X'
               zaxis(i) = ic
               xaxis(i) = ia
               if (kab.eq.ib .and. kbd.eq.ib)  xaxis(i) = ib
               if (kad.eq.id .and. kbd.eq.id)  xaxis(i) = id
               if (kab.eq.ib .and. kbd.eq.0)  xaxis(i) = ib
               yaxis(i) = 0
            else if (kad.eq.id .and. kbd.eq.id .and. kcd.eq.id) then
               polaxe(i) = 'Z-then-X'
               zaxis(i) = id
               xaxis(i) = ia
               if (kab.eq.ib .and. kbc.eq.ib)  xaxis(i) = ib
               if (kac.eq.ic .and. kbc.eq.ic)  xaxis(i) = ic
               if (kab.eq.ib .and. kbc.eq.0)  xaxis(i) = ib
               yaxis(i) = 0
            else if (kab.eq.0 .and. kac.eq.0 .and. kbc.eq.0) then
               polaxe(i) = 'Z-Bisect'
               zaxis(i) = ia
               xaxis(i) = ib
               yaxis(i) = ic
            else if (kab.eq.0 .and. kad.eq.0 .and. kbd.eq.0) then
               polaxe(i) = 'Z-Bisect'
               zaxis(i) = ia
               xaxis(i) = ib
               yaxis(i) = id
            else if (kac.eq.0 .and. kad.eq.0 .and. kcd.eq.0) then
               polaxe(i) = 'Z-Bisect'
               zaxis(i) = ia
               xaxis(i) = ic
               yaxis(i) = id
            else if (kbc.eq.0 .and. kbd.eq.0 .and. kcd.eq.0) then
               polaxe(i) = 'Z-Bisect'
               zaxis(i) = ib
               xaxis(i) = ic
               yaxis(i) = id
            else if (kab.eq.0 .and. kac.eq.ia .and. kad.eq.ia) then
               polaxe(i) = 'Bisector'
               zaxis(i) = ia
               xaxis(i) = ib
               yaxis(i) = 0
            else if (kac.eq.0 .and. kab.eq.ia .and. kad.eq.ia) then
               polaxe(i) = 'Bisector'
               zaxis(i) = ia
               xaxis(i) = ic
               yaxis(i) = 0
            else if (kad.eq.0 .and. kab.eq.ia .and. kac.eq.ia) then
               polaxe(i) = 'Bisector'
               zaxis(i) = ia
               xaxis(i) = id
               yaxis(i) = 0
            else if (kbc.eq.0 .and. kab.eq.ib .and. kbd.eq.ib) then
               polaxe(i) = 'Bisector'
               zaxis(i) = ib
               xaxis(i) = ic
               yaxis(i) = 0
            else if (kbd.eq.0 .and. kab.eq.ib .and. kbc.eq.ib) then
               polaxe(i) = 'Bisector'
               zaxis(i) = ib
               xaxis(i) = id
               yaxis(i) = 0
            else if (kcd.eq.0 .and. kac.eq.ic .and. kbc.eq.ic) then
               polaxe(i) = 'Bisector'
               zaxis(i) = ic
               xaxis(i) = id
               yaxis(i) = 0
            end if
         end if
      end do
c
c     list the local frame definition for each multipole site
c
      write (iout,10)
   10 format (/,' Local Frame Definition for Multipole Sites :')
      write (iout,20)
   20 format (/,5x,'Atom',5x,'Name',6x,'Axis Type',5x,'Z Axis',
     &           2x,'X Axis',2x,'Y Axis',/)
      do i = 1, npole
         write (iout,30)  i,name(i),polaxe(i),zaxis(i),
     &                    xaxis(i),yaxis(i)
   30    format (i8,6x,a3,7x,a8,2x,3i8)
      end do
c
c     allow the user to manually alter local coordinate frames
c
      change = .false.
      query = .true.
      i = -1
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=40,end=40)  i
         if (i .eq. 0)  query = .false.
      end if
   40 continue
      do while (query)
         i = 0
         ia = 0
         ib = 0
         ic = 0
         write (iout,50)
   50    format (/,' Enter Altered Local Frame Definition',
     &              ' [<CR>=Exit] :  ',$)
         read (input,60)  record
   60    format (a240)
         read (record,*,err=70,end=70)  i,ia,ib,ic
   70    continue
         if (i .eq. 0) then
            query = .false.
         else
            change = .true.
            if (ia .eq. 0)  polaxe(i)= 'None'
            if (ia.ne.0 .and. ib.eq.0)  polaxe(i) = 'Z-Only'
            if (ia.gt.0 .and. ib.gt.0)  polaxe(i) = 'Z-then-X'
            if (ia.lt.0 .or. ib.lt.0)  polaxe(i) = 'Bisector'
            if (ib.lt.0 .and. ic.lt.0)  polaxe(i) = 'Z-Bisect'
            if (max(ia,ib,ic) .lt. 0)  polaxe(i) = '3-Fold'
            zaxis(i) = abs(ia)
            xaxis(i) = abs(ib)
            yaxis(i) = abs(ic)
         end if
      end do
c
c     repeat local frame list if definitions were altered
c
      if (change) then
         write (iout,80)
   80    format (/,' Local Frame Definition for Multipole Sites :')
         write (iout,90)
   90    format (/,5x,'Atom',5x,'Name',6x,'Axis Type',5x,'Z Axis',
     &              2x,'X Axis',2x,'Y Axis',/)
         do i = 1, npole
            write (iout,100)  i,name(i),polaxe(i),zaxis(i),
     &                        xaxis(i),yaxis(i)
  100       format (i8,6x,a3,7x,a8,2x,3i8)
         end do
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  function priority  --  atom priority for axis assignment  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "priority" decides which of two connected atoms should be
c     preferred during construction of a local coordinate frame
c     and returns that atom number; if the two atoms have equal
c     priority then a zero is returned
c
c
      function priority (i,ia,ib)
      use sizes
      use atomid
      use couple
      implicit none
      integer i,k,m
      integer ia,ib
      integer ka,kb
      integer ka1,kb1
      integer ka2,kb2
      integer priority
c
c
c     get priority based on atomic number and connected atoms
c
      ka = atomic(ia)
      kb = atomic(ib)
      if (ka .gt. kb) then
         priority = ia
      else if (kb .gt. ka) then
         priority = ib
      else
         ka1 = 0
         ka2 = 0
         do k = 1, n12(ia)
            m = i12(k,ia)
            if (i .ne. m) then
               m = atomic(m)
               if (m .ge. ka1) then
                  ka2 = ka1
                  ka1 = m
               else if (m .gt. ka2) then
                  ka2 = m
               end if
            end if
         end do
         kb1 = 0
         kb2 = 0
         do k = 1, n12(ib)
            m = i12(k,ib)
            if (i .ne. m) then
               m = atomic(m)
               if (m .gt. kb1) then
                  kb2 = kb1
                  kb1 = m
               else if (m .gt. kb2) then
                  kb2 = m
               end if
            end if
         end do
         if (n12(ia) .lt. n12(ib)) then
            priority = ia
         else if (n12(ib) .lt. n12(ia)) then
            priority = ib
         else if (ka1 .gt. kb1) then
            priority = ia
         else if (kb1 .gt. ka1) then
            priority = ib
         else if (ka2 .gt. kb2) then
            priority = ia
         else if (kb2 .gt. ka2) then
            priority = ib
         else
            priority = 0
         end if
      end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine rotframe  --  convert multipoles to local frame  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "rotframe" takes the global multipole moments and rotates them
c     into the local coordinate frame defined at each atomic site
c
c
      subroutine rotframe
      use sizes
      use atomid
      use atoms
      use iounit
      use mpole
      use units
      implicit none
      integer i,j,ii
      integer ixaxe
      integer iyaxe
      integer izaxe
      real*8 a(3,3)
c
c
c     perform dynamic allocation of some global arrays
c
      if (.not. allocated(pole))  allocate (pole(maxpole,npole))
c
c     store the global multipoles in the local frame array
c
      do i = 1, npole
         do j = 1, 13
            pole(j,i) = rpole(j,i)
         end do
      end do
c
c     rotate the multipoles from global frame to local frame
c
      do i = 1, npole
         call rotmat (i,a)
         call invert (3,a)
         call rotsite (i,a)
      end do
c
c     copy the rotated multipoles back to local frame array
c
      do i = 1, npole
         do j = 1, 13
            pole(j,i) = rpole(j,i)
         end do
      end do
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     convert dipole and quadrupole moments back to atomic units
c
      do i = 1, npole
         rpole(1,i) = pole(1,i)
         do j = 2, 4
            rpole(j,i) = pole(j,i) / bohr
         end do
         do j = 5, 13
            rpole(j,i) = 3.0d0 * pole(j,i) / bohr**2
         end do
      end do
c
c     print the local frame Cartesian atomic multipoles
c
      write (iout,10)
   10 format (/,' Local Frame Cartesian Multipole Moments :')
      do ii = 1, npole
         i = pollist(ii)
         if (i .eq. 0) then
            write (iout,20)  ii,name(ii),atomic(ii)
   20       format (/,' Atom:',i8,9x,'Name:',3x,a3,7x,
     &                 'Atomic Number:',i8)
            write (iout,30)
   30       format (/,' No Atomic Multipole Moments for this Site')
         else
            izaxe = zaxis(i)
            ixaxe = xaxis(i)
            iyaxe = yaxis(i)
            if (iyaxe .lt. 0)  iyaxe = -iyaxe
            write (iout,40)  ii,name(ii),atomic(ii)
   40       format (/,' Atom:',i8,9x,'Name:',3x,a3,
     &                 7x,'Atomic Number:',i8)
            write (iout,50)  polaxe(i),izaxe,ixaxe,iyaxe
   50       format (/,' Local Frame:',12x,a8,6x,3i8)
            write (iout,60)  rpole(1,i)
   60       format (/,' Charge:',10x,f15.5)
            write (iout,70)  rpole(2,i),rpole(3,i),rpole(4,i)
   70       format (' Dipole:',10x,3f15.5)
            write (iout,80)  rpole(5,i)
   80       format (' Quadrupole:',6x,f15.5)
            write (iout,90)  rpole(8,i),rpole(9,i)
   90       format (18x,2f15.5)
            write (iout,100)  rpole(11,i),rpole(12,i),rpole(13,i)
  100       format (18x,3f15.5)
         end if
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine fixframe  --  alter the local frame definition  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "fixframe" is a service routine that alters the local frame
c     definition for specified atoms
c
c
      subroutine fixframe
      use sizes
      use atomid
      use atoms
      use couple
      use files
      use keys
      use kpolr
      use iounit
      use mpole
      use polar
      use units
      implicit none
      integer i,j,k,ii
      integer ia,ib,ic
      integer ixaxe
      integer iyaxe
      integer izaxe
      real*8 eps,ci,cj
      real*8 big,sum
      real*8 a(3,3)
      logical query,change
      character*240 record
c
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     list the local frame definition for each multipole site
c
      write (iout,10)
   10 format (/,' Local Frame Definition for Multipole Sites :')
      write (iout,20)
   20 format (/,5x,'Atom',5x,'Name',6x,'Axis Type',5x,'Z Axis',2x,
     &          'X Axis',2x,'Y Axis',/)
      do ii = 1, n
         i = pollist(ii)
         if (i .eq. 0) then
            write (iout,30)  ii,name(ii)
   30       format (i8,6x,a3,10x,'--',11x,'--',6x,'--',6x,'--')
         else
            izaxe = zaxis(i)
            ixaxe = xaxis(i)
            iyaxe = yaxis(i)
            if (iyaxe .lt. 0)  iyaxe = -iyaxe
            write (iout,40)  ii,name(ii),polaxe(i),izaxe,ixaxe,iyaxe
   40       format (i8,6x,a3,7x,a8,2x,3i8)
         end if
      end do
c
c     allow the user to manually alter local coordinate frames
c
      query = .true.
      change = .false.
      do while (query)
         i = 0
         ii = 0
         ia = 0
         ib = 0
         ic = 0
         write (iout,50)
   50    format (/,' Enter Altered Local Frame Definition',
     &              ' [<CR>=Exit] :  ',$)
         read (input,60)  record
   60    format (a240)
         read (record,*,err=70,end=70)  ii,ia,ib,ic
   70    continue
         if (ii .eq. 0) then
            query = .false.
         else
            i = pollist(ii)
         end if
         if (i .ne. 0) then
            change = .true.
            if (ia .eq. 0)  polaxe(i) = 'None'
            if (ia.ne.0 .and. ib.eq.0)  polaxe(i) = 'Z-Only'
            if (ia.gt.0 .and. ib.gt.0)  polaxe(i) = 'Z-then-X'
            if (ia.lt.0  .or. ib.lt.0)  polaxe(i) = 'Bisector'
            if (ib.lt.0 .and. ic.lt.0)  polaxe(i) = 'Z-Bisect'
            if (max(ia,ib,ic)  .lt. 0)  polaxe(i) = '3-Fold'
            zaxis(i) = abs(ia)
            xaxis(i) = abs(ib)
            yaxis(i) = abs(ic)
         end if
      end do
c
c     repeat local frame list if definitions were altered
c
      if (change) then
         write (iout,80)
   80    format (/,' Local Frame Definition for Multipole Sites :')
         write (iout,90)
   90    format (/,5x,'Atom',5x,'Name',6x,'Axis Type',5x,'Z Axis',2x,
     &             'X Axis',2x,'Y Axis',/)
         do ii = 1, npole
            i = pollist(ii)
            if (i .eq. 0) then
               write (iout,100)  ii,name(ii)
  100          format (i8,6x,a3,10x,'--',11x,'--',6x,'--',6x,'--')
            else
               izaxe = zaxis(i)
               ixaxe = xaxis(i)
               iyaxe = yaxis(i)
               if (iyaxe .lt. 0)  iyaxe = -iyaxe
               write (iout,110)  ii,name(ii),polaxe(i),izaxe,ixaxe,iyaxe
  110          format (i8,6x,a3,7x,a8,2x,3i8)
            end if
         end do
      end if
c
c     store the global multipoles in the local frame array
c
      do i = 1, npole
         do j = 1, 13
            pole(j,i) = rpole(j,i)
         end do
      end do
c
c     rotate the multipoles from global frame to local frame
c
      do i = 1, npole
         call rotmat (i,a)
         call invert (3,a)
         call rotsite (i,a)
      end do
c
c     copy the rotated multipoles back to local frame array
c
      do i = 1, npole
         do j = 1, 13
            pole(j,i) = rpole(j,i)
         end do
      end do
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     convert dipole and quadrupole moments back to atomic units
c
      do i = 1, npole
         pole(1,i) = pole(1,i)
         do j = 2, 4
            pole(j,i) = pole(j,i) / bohr
         end do
         do j = 5, 13
            pole(j,i) = 3.0d0 * pole(j,i) / bohr**2
         end do
      end do
c
c     regularize the multipole moments to desired precision
c
      eps = 0.00001d0
      do i = 1, npole
         do j = 1, 13
            pole(j,i) = dble(nint(pole(j,i)/eps)) * eps
         end do
      end do
c
c     maintain integer net charge for the whole system
c
      k = 0
      big = 0.0d0
      sum = 0.0d0
      do i = 1, n
         sum = sum + pole(1,i)
         ci = abs(pole(1,i))
         if (ci .gt. big) then
            do j = 1, n
               cj = abs(pole(1,j))
               if (i.ne.j .and. ci.eq.cj)  goto 120
            end do
            k = i
            big = ci
  120       continue
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
c     print the altered local frame atomic multipole values
c
      write (iout,130)
  130 format (/,' Multipoles With Altered Local Frame Definition :')
      do ii = 1, n
         i = pollist(ii)
         if (i .eq. 0) then
            write (iout,140)  ii,name(ii),atomic(ii)
  140       format (/,' Atom:',i8,9x,'Name:',3x,a3,7x,
     &                 'Atomic Number:',i8)
            write (iout,150)
  150       format (/,' No Atomic Multipole Moments for this Site')
         else
            izaxe = zaxis(i)
            ixaxe = xaxis(i)
            iyaxe = yaxis(i)
            if (iyaxe .lt. 0)  iyaxe = -iyaxe
            write (iout,160)  ii,name(ii),atomic(ii)
  160       format (/,' Atom:',i8,9x,'Name:',3x,a3,
     &                 7x,'Atomic Number:',i8)
            write (iout,170)  polaxe(i),izaxe,ixaxe,iyaxe
  170       format (/,' Local Frame:',12x,a8,6x,3i8)
            write (iout,180)  pole(1,i)
  180       format (/,' Charge:',10x,f15.5)
            write (iout,190)  pole(2,i),pole(3,i),pole(4,i)
  190       format (' Dipole:',10x,3f15.5)
            write (iout,200)  pole(5,i)
  200       format (' Quadrupole:',6x,f15.5)
            write (iout,210)  pole(8,i),pole(9,i)
  210       format (18x,2f15.5)
            write (iout,220)  pole(11,i),pole(12,i),pole(13,i)
  220       format (18x,3f15.5)
         end if
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine setpolar  --  define polarization and groups  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "setpolar" assigns atomic polarizability and Thole damping
c     parameters and allows user alteration of these values
c
c
      subroutine setpolar
      use sizes
      use atomid
      use atoms
      use couple
      use iounit
      use kpolr
      use mpole
      use polar
      use polgrp
      implicit none
      integer i,j,k,ii
      integer ia,ib
      real*8 pol,thl
      logical exist,query
      logical change
      character*240 record
      character*240 string
c
c
c     list the polariability values for each multipole site
c
      write (iout,10)
   10 format (/,' Atomic Polarizabilities for Multipole Sites :')
      write (iout,20)
   20 format (/,5x,'Atom',5x,'Name',7x,'Polarize',10x,'Thole',/)
      do ii = 1, n
         i = pollist(ii)
         if (i .eq. 0) then
            write (iout,30)  ii,name(ii)
   30       format (i8,6x,a3,12x,'--',13x,'--')
         else
            write (iout,40)  ii,name(ii),polarity(i),thole(i)
   40       format (i8,6x,a3,4x,f12.4,3x,f12.4)
         end if
      end do
c
c     allow the user to manually alter polarizability values
c
      change = .false.
      query = .true.
      i = -1
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=50,end=50)  i
         if (i .eq. 0)  query = .false.
      end if
   50 continue
      do while (query)
         i = 0
         ii = 0
         pol = 0.0d0
         thl = 0.39d0
         write (iout,60)
   60    format (/,' Enter Atom Number & Polarizability Values',
     &              ' [<CR>=Exit] :  ',$)
         read (input,70)  record
   70    format (a240)
         read (record,*,err=80,end=80)  ii,pol,thl
   80    continue
         if (ii .eq. 0) then
            query = .false.
         else
            i = pollist(ii)
         end if
         if (i .ne. 0) then
            change = .true.
            polarity(i) = pol
            thole(i) = thl
         end if
      end do
c
c     repeat polarizability values if parameters were altered
c
      if (change) then
         write (iout,90)
   90    format (/,' Atomic Polarizabilities for Multipole Sites :')
         write (iout,100)
  100    format (/,5x,'Atom',5x,'Name',7x,'Polarize',10x,'Thole',/)
         do ii = 1, n
            i = pollist(ii)
            if (i .eq. 0) then
               write (iout,110)  ii,name(ii)
  110          format (i8,6x,a3,12x,'--',13x,'--')
            else
               write (iout,120)  ii,name(ii),polarity(i),thole(i)
  120          format (i8,6x,a3,4x,f12.4,3x,f12.4)
            end if
         end do
      end if
c
c     use bonded atoms as initial guess at polarization groups
c
      write (iout,130)
  130 format (/,' The default is to place all Atoms into one',
     &           ' Polarization Group;',
     &        /,' This can be altered by entering a series of',
     &           ' Bonded Atom Pairs',
     &        /,' that separate the Molecule into distinct',
     &           ' Polarization Groups')
      do i = 1, n
         do j = 1, n12(i)
            pgrp(j,i) = i12(j,i)
         end do
      end do
c
c     get the bonds that separate the polarization groups
c
      query = .true.
      i = -1
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=140,end=140)  i
         if (i .eq. 0)  query = .false.
      end if
  140 continue
      do while (query)
         ia = 0
         ib = 0
         write (iout,150)
  150    format (/,' Enter a Bond between Polarization Groups',
     &              ' [<CR>=Exit] :  ',$)
         read (input,160)  record
  160    format (a240)
         read (record,*,err=170,end=170)  ia,ib
  170    continue
         if (ia.eq.0 .or. ib.eq.0) then
            query = .false.
         else
            do i = 1, n12(ia)
               if (pgrp(i,ia) .eq. ib) then
                  do j = i+1, n12(ia)
                     pgrp(j-1,ia) = pgrp(j,ia)
                  end do
                  pgrp(n12(ia),ia) = 0
               end if
            end do
            do i = 1, n12(ib)
               if (pgrp(i,ib) .eq. ia) then
                  do j = i+1, n12(ib)
                     pgrp(j-1,ib) = pgrp(j,ib)
                  end do
                  pgrp(n12(ib),ib) = 0
               end if
            end do
         end if
      end do
      call polargrp
c
c     list the polarization group for each multipole site
c
      write (iout,180)
  180 format (/,' Polarization Groups for Multipole Sites :')
      write (iout,190)
  190 format (/,5x,'Atom',5x,'Name',7x,'Polarization Group',
     &           ' Definition',/)
      do i = 1, n
         k = 0
         do j = 1, maxval
            if (pgrp(j,i) .ne. 0)  k = j
         end do
         write (iout,200)  i,name(i),(pgrp(j,i),j=1,k)
  200    format (i8,6x,a3,8x,20i6)
      end do
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine alterpol  --  alter multipoles for polarization  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "alterpol" finds an output set of TINKER multipole parameters
c     which when used with an intergroup polarization model will
c     give the same electrostatic potential around the molecule as
c     the input set of multipole parameters with all atoms in one
c     polarization group
c
c     for example, the input parameters could be from a distributed
c     multipole analysis of a molecular wavefunction and the output
c     will be the parameter values that achieve the same potential
c     in the presence of intergroup (intramolecular) polarization
c
c
      subroutine alterpol
      use sizes
      use atomid
      use atoms
      use iounit
      use mpole
      use polar
      use units
      implicit none
      integer i,j,ii
      integer ixaxe
      integer iyaxe
      integer izaxe
      real*8 a(3,3)
c
c
c     rotate the multipole components into the global frame
c
      do i = 1, npole
         call rotpole
      end do
c
c     compute induced dipoles to be removed from QM multipoles
c
      call interpol
c
c     remove induced dipole from global frame multipoles
c
      do i = 1, npole
         rpole(2,i) = rpole(2,i) - uind(1,i)
         rpole(3,i) = rpole(3,i) - uind(2,i)
         rpole(4,i) = rpole(4,i) - uind(3,i)
         do j = 1, 13
            pole(j,i) = rpole(j,i)
         end do
      end do
c
c     rotate the multipoles from global frame to local frame
c
      do i = 1, npole
         call rotmat (i,a)
         call invert (3,a)
         call rotsite (i,a)
      end do
c
c     copy the rotated multipoles back to local frame array
c
      do i = 1, npole
         do j = 1, 13
            pole(j,i) = rpole(j,i)
         end do
      end do
c
c     convert dipole and quadrupole moments back to atomic units
c
      do i = 1, npole
         rpole(1,i) = pole(1,i)
         do j = 2, 4
            rpole(j,i) = pole(j,i) / bohr
         end do
         do j = 5, 13
            rpole(j,i) = 3.0d0 * pole(j,i) / bohr**2
         end do
      end do
c
c     print multipoles with intergroup polarization removed
c
      write (iout,10)
   10 format (/,' Multipoles after Removal of Intergroup',
     &           ' Polarization :')
      do ii = 1, n
         i = pollist(ii)
         if (i .eq. 0) then
            write (iout,20)  ii,name(ii),atomic(ii)
   20       format (/,' Atom:',i8,9x,'Name:',3x,a3,
     &                 7x,'Atomic Number:',i8)
            write (iout,30)
   30       format (/,' No Atomic Multipole Moments for this Site')
         else
            izaxe = zaxis(i)
            ixaxe = xaxis(i)
            iyaxe = yaxis(i)
            if (iyaxe .lt. 0)  iyaxe = -iyaxe
            write (iout,40)  ii,name(ii),atomic(ii)
   40       format (/,' Atom:',i8,9x,'Name:',3x,a3,
     &                 7x,'Atomic Number:',i8)
            write (iout,50)  polaxe(i),izaxe,ixaxe,iyaxe
   50       format (/,' Local Frame:',12x,a8,6x,3i8)
            write (iout,60)  rpole(1,i)
   60       format (/,' Charge:',10x,f15.5)
            write (iout,70)  rpole(2,i),rpole(3,i),rpole(4,i)
   70       format (' Dipole:',10x,3f15.5)
            write (iout,80)  rpole(5,i)
   80       format (' Quadrupole:',6x,f15.5)
            write (iout,90)  rpole(8,i),rpole(9,i)
   90       format (18x,2f15.5)
            write (iout,100)  rpole(11,i),rpole(12,i),rpole(13,i)
  100       format (18x,3f15.5)
         end if
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine interpol  --  get intergroup induced dipoles  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "interpol" is a service routine that computes induced dipole
c     moments for use during removal of intergroup polarization
c
c
      subroutine interpol
      use sizes
      use atoms
      use iounit
      use mpole
      use polar
      use polpot
      use units
      implicit none
      integer i,j,k,iter
      integer maxiter
      real*8 eps,epsold
      real*8 polmin,norm
      real*8 a,b,sum
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: poli(:)
      real*8, allocatable :: field(:,:)
      real*8, allocatable :: rsd(:,:)
      real*8, allocatable :: zrsd(:,:)
      real*8, allocatable :: conj(:,:)
      real*8, allocatable :: vec(:,:)
      logical done
c
c
c     perform dynamic allocation of some global arrays
c
      if (.not. allocated(uind))  allocate (uind(3,npole))
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
      allocate (poli(npole))
      allocate (field(3,npole))
      allocate (rsd(3,npole))
      allocate (zrsd(3,npole))
      allocate (conj(3,npole))
      allocate (vec(3,npole))
c
c     compute induced dipoles as polarizability times field
c
      call dfieldi (field,pscale)
      do i = 1, npole
         do j = 1, 3
            uind(j,i) = polarity(i) * field(j,i)
         end do
      end do
c
c     for direct-only models set mutual scale factors to zero
c
      if (poltyp .eq. 'DIRECT') then
         u1scale = 0.0d0
         u2scale = 0.0d0
         u3scale = 0.0d0
         u4scale = 0.0d0
      end if
c
c     compute intergroup induced dipole moments via CG algorithm
c
      done = .false.
      maxiter = 500
      iter = 0
      eps = 100.0d0
      polmin = 0.00000001d0
      call ufieldi (field,pscale)
      do i = 1, npole
         poli(i) = max(polmin,polarity(i))
         do j = 1, 3
            rsd(j,i) = field(j,i)
            zrsd(j,i) = rsd(j,i) * poli(i)
            conj(j,i) = zrsd(j,i)
         end do
      end do
c
c     iterate the intergroup induced dipoles and check convergence
c
      do while (.not. done)
         iter = iter + 1
         do i = 1, npole
            do j = 1, 3
               vec(j,i) = uind(j,i)
               uind(j,i) = conj(j,i)
            end do
         end do
         call ufieldi (field,pscale)
         do i = 1, npole
            do j = 1, 3
               uind(j,i) = vec(j,i)
               vec(j,i) = conj(j,i)/poli(i) - field(j,i)
            end do
         end do
         a = 0.0d0
         sum = 0.0d0
         do i = 1, npole
            do j = 1, 3
               a = a + conj(j,i)*vec(j,i)
               sum = sum + rsd(j,i)*zrsd(j,i)
            end do
         end do
         if (a .ne. 0.0d0)  a = sum / a
         do i = 1, npole
            do j = 1, 3
               uind(j,i) = uind(j,i) + a*conj(j,i)
               rsd(j,i) = rsd(j,i) - a*vec(j,i)
            end do
         end do
         b = 0.0d0
         do i = 1, npole
            do j = 1, 3
               zrsd(j,i) = rsd(j,i) * poli(i)
               b = b + rsd(j,i)*zrsd(j,i)
            end do
         end do
         if (sum .ne. 0.0d0)  b = b / sum
         eps = 0.0d0
         do i = 1, npole
            do j = 1, 3
               conj(j,i) = zrsd(j,i) + b*conj(j,i)
               eps = eps + rsd(j,i)*rsd(j,i)
            end do
         end do
         eps = debye * sqrt(eps/dble(npolar))
         epsold = eps
         if (iter .eq. 1) then
            write (iout,10)
   10       format (/,' Determination of Intergroup Induced',
     &                 ' Dipoles :',
     &              //,4x,'Iter',8x,'RMS Change (Debyes)',/)
         end if
         write (iout,20)  iter,eps
   20    format (i8,7x,f16.10)
         if (eps .lt. poleps)  done = .true.
         if (eps .gt. epsold)  done = .true.
         if (iter .ge. maxiter)  done = .true.
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (pscale)
      deallocate (poli)
      deallocate (field)
      deallocate (rsd)
      deallocate (zrsd)
      deallocate (conj)
      deallocate (vec)
c
c     terminate the calculation if dipoles failed to converge
c
      if (eps .gt. poleps) then
         write (iout,30)
   30    format (/,' INTERPOL  --  Warning, Induced Dipoles',
     &              ' are not Converged')
         call prterr
         call fatal
      end if
c
c     print out a list of the final induced dipole moments
c
      write (iout,40)
   40 format (/,' Intergroup Induced Dipoles to be Removed',
     &           ' (Debyes) :')
      write (iout,50)
   50 format (/,4x,'Atom',14x,'X',11x,'Y',11x,'Z',
     &           9x,'Total'/)
      do i = 1, npole
         k = ipole(i)
         norm = sqrt(uind(1,i)**2+uind(2,i)**2+uind(3,i)**2)
         write (iout,60)  k,(debye*uind(j,i),j=1,3),debye*norm
   60    format (i8,5x,4f12.4)
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine dfieldi  --  find permanent multipole field  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "dfieldi" computes the electrostatic field due to permanent
c     multipole moments
c
c
      subroutine dfieldi (field,pscale)
      use sizes
      use atoms
      use mpole
      use polar
      use polgrp
      use polpot
      implicit none
      integer i,j,k,ii,kk
      real*8 r,r2,xr,yr,zr
      real*8 rr3,rr5,rr7
      real*8 pdi,pti,pgamma
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dir,dkr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 scale3,scale5
      real*8 scale7,damp
      real*8 fi(3),fk(3)
      real*8 pscale(*)
      real*8 field(3,*)
c
c
c     zero out the induced dipole and the field at each site
c
      do i = 1, npole
         do j = 1, 3
            uind(j,i) = 0.0d0
            field(j,i) = 0.0d0
         end do
      end do
c
c     compute the electrostatic field due to permanent multipoles
c
      do i = 1, npole-1
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         do j = i+1, npole
            pscale(ipole(j)) = 1.0d0
         end do
         do j = 1, np11(ii)
            pscale(ip11(j,ii)) = d1scale
         end do
         do j = 1, np12(ii)
            pscale(ip12(j,ii)) = d2scale
         end do
         do j = 1, np13(ii)
            pscale(ip13(j,ii)) = d3scale
         end do
         do j = 1, np14(ii)
            pscale(ip14(j,ii)) = d4scale
         end do
         do k = i+1, npole
            kk = ipole(k)
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
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
            scale3 = pscale(kk)
            scale5 = pscale(kk)
            scale7 = pscale(kk)
            damp = pdi * pdamp(k)
            if (damp .ne. 0.0d0) then
               pgamma = min(pti,thole(k))
               damp = -pgamma * (r/damp)**3
               if (damp .gt. -50.0d0) then
                  scale3 = scale3 * (1.0d0-exp(damp))
                  scale5 = scale5 * (1.0d0-(1.0d0-damp)*exp(damp))
                  scale7 = scale7 * (1.0d0-(1.0d0-damp+0.6d0*damp**2)
     &                                               *exp(damp))
               end if
            end if
            rr3 = scale3 / (r*r2)
            rr5 = 3.0d0 * scale5 / (r*r2*r2)
            rr7 = 15.0d0 * scale7 / (r*r2*r2*r2)
            dir = dix*xr + diy*yr + diz*zr
            qix = qixx*xr + qixy*yr + qixz*zr
            qiy = qixy*xr + qiyy*yr + qiyz*zr
            qiz = qixz*xr + qiyz*yr + qizz*zr
            qir = qix*xr + qiy*yr + qiz*zr
            dkr = dkx*xr + dky*yr + dkz*zr
            qkx = qkxx*xr + qkxy*yr + qkxz*zr
            qky = qkxy*xr + qkyy*yr + qkyz*zr
            qkz = qkxz*xr + qkyz*yr + qkzz*zr
            qkr = qkx*xr + qky*yr + qkz*zr
            fi(1) = -xr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                 - rr3*dkx + 2.0d0*rr5*qkx
            fi(2) = -yr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                 - rr3*dky + 2.0d0*rr5*qky
            fi(3) = -zr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                 - rr3*dkz + 2.0d0*rr5*qkz
            fk(1) = xr*(rr3*ci+rr5*dir+rr7*qir)
     &                 - rr3*dix - 2.0d0*rr5*qix
            fk(2) = yr*(rr3*ci+rr5*dir+rr7*qir)
     &                 - rr3*diy - 2.0d0*rr5*qiy
            fk(3) = zr*(rr3*ci+rr5*dir+rr7*qir)
     &                 - rr3*diz - 2.0d0*rr5*qiz
            do j = 1, 3
               field(j,i) = field(j,i) + fi(j)
               field(j,k) = field(j,k) + fk(j)
            end do
         end do
      end do
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine ufieldi  --  find induced intergroup field  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "ufieldi" computes the electrostatic field due to intergroup
c     induced dipole moments
c
c
      subroutine ufieldi (field,pscale)
      use sizes
      use atoms
      use mpole
      use polar
      use polgrp
      use polpot
      implicit none
      integer i,j,k,ii,kk
      real*8 xr,yr,zr
      real*8 r,r2,rr3,rr5
      real*8 pdi,pti,pgamma
      real*8 uix,uiy,uiz
      real*8 ukx,uky,ukz
      real*8 uir,ukr,damp
      real*8 scale3,scale5
      real*8 fi(3),fk(3)
      real*8 pscale(*)
      real*8 field(3,*)
c
c
c     zero out the value of the field at each site
c
      do i = 1, npole
         do j = 1, 3
            field(j,i) = 0.0d0
         end do
      end do
c
c     find the electrostatic field due to induced dipoles
c
      do i = 1, npole-1
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         uix = uind(1,i)
         uiy = uind(2,i)
         uiz = uind(3,i)
         do j = i+1, npole
            pscale(ipole(j)) = 0.0d0
         end do
         do j = 1, np11(ii)
            pscale(ip11(j,ii)) = u1scale - d1scale
         end do
         do j = 1, np12(ii)
            pscale(ip12(j,ii)) = u2scale - d2scale
         end do
         do j = 1, np13(ii)
            pscale(ip13(j,ii)) = u3scale - d3scale
         end do
         do j = 1, np14(ii)
            pscale(ip14(j,ii)) = u4scale - d4scale
         end do
         do k = i+1, npole
            kk = ipole(k)
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            r2 = xr*xr + yr* yr + zr*zr
            r = sqrt(r2)
            ukx = uind(1,k)
            uky = uind(2,k)
            ukz = uind(3,k)
            scale3 = pscale(kk)
            scale5 = pscale(kk)
            damp = pdi * pdamp(k)
            if (damp .ne. 0.0d0) then
               pgamma = min(pti,thole(k))
               damp = -pgamma * (r/damp)**3
               if (damp .gt. -50.0d0) then
                  scale3 = scale3*(1.0d0-exp(damp))
                  scale5 = scale5*(1.0d0-(1.0d0-damp)*exp(damp))
               end if
            end if
            rr3 = scale3 / (r*r2)
            rr5 = 3.0d0 * scale5 / (r*r2*r2)
            uir = xr*uix + yr*uiy + zr*uiz
            ukr = xr*ukx + yr*uky + zr*ukz
            fi(1) = -rr3*ukx + rr5*ukr*xr
            fi(2) = -rr3*uky + rr5*ukr*yr
            fi(3) = -rr3*ukz + rr5*ukr*zr
            fk(1) = -rr3*uix + rr5*uir*xr
            fk(2) = -rr3*uiy + rr5*uir*yr
            fk(3) = -rr3*uiz + rr5*uir*zr
            do j = 1, 3
               field(j,i) = field(j,i) + fi(j)
               field(j,k) = field(j,k) + fk(j)
            end do
         end do
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine fixpolar  --  postprocess multipole moments  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "fixpolar" averages multipoles over equivalent sites, sets
c     symmetric components at achiral atoms to zero, and maintains
c     an integer net charge and traceless quadrupoles
c
c
      subroutine fixpolar
      use sizes
      use atomid
      use atoms
      use couple
      use iounit
      use mpole
      use units
      implicit none
      integer i,j,k,m
      integer ii,kk
      integer ixaxe
      integer iyaxe
      integer izaxe
      integer next,nx
      integer nlist
      integer list(maxval)
      real*8 eps,ci,cj
      real*8 big,sum
      real*8 pave(13)
      logical exist,query
      logical yzero
      character*1 answer
      character*240 record
      character*240 string
c
c
c     optionally average multipoles for equivalent atoms
c
      query = .true.
      answer = ' '
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  answer
         call upcase (answer)
         if (answer.eq.'N' .or. answer.eq.'Y')  query = .false.
      end if
   10 continue
      if (query) then
         write (iout,20)
   20    format (/,' Average the Multipole Moments of Equivalent',
     &              ' Atoms [N] :  ',$)
         read (input,30)  record
   30    format (a240)
         next = 1
         call gettext (record,answer,next)
         call upcase (answer)
      end if
c
c     perform averaging for equivalent monovalent atoms
c
      if (answer .eq. 'Y') then
         do i = 1, n
            nlist = 0
            do j = 1, n12(i)
               k = i12(j,i)
               if (n12(k) .eq. 1) then
                  do m = 1, nlist
                     if (list(m) .eq. atomic(k))  goto 40
                  end do
                  nlist = nlist + 1
                  list(nlist) = atomic(k)
   40             continue
               end if
            end do
            do ii = 1, nlist
               kk = list(ii)
               do j = 1, 13
                  pave(j) = 0.0d0
               end do
               nx = 0
               do j = 1, n12(i)
                  k = i12(j,i)
                  if (pollist(k) .ne. 0) then
                     if (atomic(k).eq.kk .and. n12(k).eq.1) then
                        nx = nx + 1
                        do m = 1, 13
                           pave(m) = pave(m) + pole(m,k)
                        end do
                     end if
                  end if
               end do
               if (nx .ge. 2) then
                  do j = 1, 13
                     pave(j) = pave(j) / dble(nx)
                  end do
                  do j = 1, n12(i)
                     k = i12(j,i)
                     if (pollist(k) .ne. 0) then
                        if (atomic(k).eq.kk .and. n12(k).eq.1) then
                           do m = 1, 13
                              pole(m,k) = pave(m)
                           end do
                        end if
                     end if
                  end do
               end if
            end do
         end do
      end if
c
c     optionally set symmetric multipole components to zero
c
      query = .true.
      answer = ' '
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=50,end=50)  answer
         call upcase (answer)
         if (answer.eq.'N' .or. answer.eq.'Y')  query = .false.
      end if
   50 continue
      if (query) then
         write (iout,60)
   60    format (/,' Remove Multipole Components Zeroed by',
     &              ' Symmetry [N] :  ',$)
         read (input,70)  record
   70    format (a240)
         next = 1
         call gettext (record,answer,next)
         call upcase (answer)
      end if
c
c     remove multipole components that are zero by symmetry
c
      if (answer .eq. 'Y') then
         do i = 1, npole
            yzero = .false.
            if (yaxis(i) .eq. 0)  yzero = .true.
            if (polaxe(i) .eq. 'Bisector')  yzero = .true.
            if (polaxe(i) .eq. 'Z-Bisect')  yzero = .true.
            if (zaxis(i).eq.0 .or. zaxis(i).gt.n) then
               pole(13,i) = 0.0d0
            end if
            if (xaxis(i).eq.0 .or. xaxis(i).gt.n) then
               pole(2,i) = 0.0d0
               pole(5,i) = -0.5d0 * pole(13,i)
               pole(7,i) = 0.0d0
               pole(9,i) = pole(5,i)
               pole(11,i) = 0.0d0
            end if
            if (yzero) then
               pole(3,i) = 0.0d0
               pole(6,i) = 0.0d0
               pole(8,i) = 0.0d0
               pole(10,i) = 0.0d0
               pole(12,i) = 0.0d0
            end if
         end do
      end if
c
c     convert dipole and quadrupole moments back to atomic units
c
      do i = 1, npole
         pole(1,i) = pole(1,i)
         do j = 2, 4
            pole(j,i) = pole(j,i) / bohr
         end do
         do j = 5, 13
            pole(j,i) = 3.0d0 * pole(j,i) / bohr**2
         end do
      end do
c
c     regularize the multipole moments to desired precision
c
      eps = 0.00001d0
      do i = 1, npole
         do j = 1, 13
            pole(j,i) = dble(nint(pole(j,i)/eps)) * eps
         end do
      end do
c
c     maintain integer net charge for the whole system
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
               if (i.ne.j .and. ci.eq.cj)  goto 80
            end do
            k = i
            big = ci
   80       continue
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
c     print the final post-processed multipoles for AMOEBA
c
      write (iout,90)
   90 format (/,' Final Multipole Moments for the AMOEBA Force',
     &           ' Field :')
      do ii = 1, n
         i = pollist(ii)
         if (i .eq. 0) then
            write (iout,100)  ii,name(ii),atomic(ii)
  100       format (/,' Atom:',i8,9x,'Name:',3x,a3,
     &                 7x,'Atomic Number:',i8)
            write (iout,110)
  110       format (/,' No Atomic Multipole Moments for this Site')
         else
            izaxe = zaxis(i)
            ixaxe = xaxis(i)
            iyaxe = yaxis(i)
            if (iyaxe .lt. 0)  iyaxe = -iyaxe
            write (iout,120)  ii,name(ii),atomic(ii)
  120       format (/,' Atom:',i8,9x,'Name:',3x,a3,
     &                 7x,'Atomic Number:',i8)
            write (iout,130)  polaxe(i),izaxe,ixaxe,iyaxe
  130       format (/,' Local Frame:',12x,a8,6x,3i8)
            write (iout,140)  pole(1,i)
  140       format (/,' Charge:',10x,f15.5)
            write (iout,150)  pole(2,i),pole(3,i),pole(4,i)
  150       format (' Dipole:',10x,3f15.5)
            write (iout,160)  pole(5,i)
  160       format (' Quadrupole:',6x,f15.5)
            write (iout,170)  pole(8,i),pole(9,i)
  170       format (18x,2f15.5)
            write (iout,180)  pole(11,i),pole(12,i),pole(13,i)
  180       format (18x,3f15.5)
         end if
      end do
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine prtpolar  --  create file with final multipoles  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "prtpolar" makes a key file containing results from distributed
c     multipole analysis or removal of intramolecular polarization
c
c
      subroutine prtpolar
      use sizes
      use atoms
      use atomid
      use files
      use keys
      use kpolr
      use mpole
      use polar
      use units
      implicit none
      integer i,j,k,it
      integer ikey,size
      integer ixaxe
      integer iyaxe
      integer izaxe
      integer freeunit
      integer trimtext
      logical dofull
      character*240 keyfile
      character*240 record
c
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
         write (ikey,10)  record(1:size)
   10    format (a)
      end do
      if (nkey .ne. 0) then
         write (ikey,20)
   20    format ()
      end if
c
c     output the atom definitions to the keyfile as appropriate
c
      dofull = .true.
      do i = 1, n
         if (type(i) .ne. i)  dofull = .false.
      end do
      if (dofull) then
         do i = 1, n
            write (ikey,30)  i,i,name(i),story(i),atomic(i),
     &                       mass(i),valence(i)
   30       format ('atom',6x,2i5,4x,a3,3x,'"',a20,'"',i10,f10.3,i5)
         end do
         if (n .ne. 0) then
            write (ikey,40)
   40       format ()
         end if
      end if
c
c     output the local frame multipole values to the keyfile
c
      do i = 1, npole
         it = ipole(i)
         if (.not. dofull)  it = -it
         izaxe = zaxis(i)
         ixaxe = xaxis(i)
         iyaxe = yaxis(i)
         if (iyaxe .lt. 0)  iyaxe = -iyaxe
         if (polaxe(i) .eq. 'None') then
            write (ikey,50)  it,pole(1,i)
   50       format ('multipole',1x,i5,21x,f11.5)
         else if (polaxe(i) .eq. 'Z-Only') then
            write (ikey,60)  it,izaxe,pole(1,i)
   60       format ('multipole',1x,2i5,16x,f11.5)
         else if (polaxe(i) .eq. 'Z-then-X') then
            if (yaxis(i) .eq. 0) then
               write (ikey,70)  it,izaxe,ixaxe,pole(1,i)
   70          format ('multipole',1x,3i5,11x,f11.5)
            else
               write (ikey,80)  it,izaxe,ixaxe,iyaxe,pole(1,i)
   80          format ('multipole',1x,4i5,6x,f11.5)
            end if
         else if (polaxe(i) .eq. 'Bisector') then
            if (yaxis(i) .eq. 0) then
               write (ikey,90)  it,-izaxe,-ixaxe,pole(1,i)
   90          format ('multipole',1x,3i5,11x,f11.5)
            else
               write (ikey,100)  it,-izaxe,-ixaxe,iyaxe,pole(1,i)
  100          format ('multipole',1x,4i5,6x,f11.5)
            end if
         else if (polaxe(i) .eq. 'Z-Bisect') then
            write (ikey,110)  it,izaxe,-ixaxe,-iyaxe,pole(1,i)
  110       format ('multipole',1x,4i5,6x,f11.5)
         else if (polaxe(i) .eq. '3-Fold') then
            write (ikey,120)  it,-izaxe,-ixaxe,-iyaxe,pole(1,i)
  120       format ('multipole',1x,4i5,6x,f11.5)
         end if
         write (ikey,130)  pole(2,i),pole(3,i),pole(4,i)
  130    format (36x,3f11.5)
         write (ikey,140)  pole(5,i)
  140    format (36x,f11.5)
         write (ikey,150)  pole(8,i),pole(9,i)
  150    format (36x,2f11.5)
         write (ikey,160)  pole(11,i),pole(12,i),pole(13,i)
  160    format (36x,3f11.5)
      end do
c
c     output the polarizability parameters to the keyfile
c
      if (dofull) then
         if (n .ne. 0) then
            write (ikey,170)
  170       format ()
         end if
         do i = 1, npole
            k = 0
            do j = 1, maxval
               if (pgrp(j,i) .ne. 0)  k = j
            end do
            write (ikey,180)  ipole(i),polarity(i),thole(i),
     &                        (pgrp(j,i),j=1,k)
  180       format ('polarize',2x,i5,5x,2f11.4,2x,20i5)
         end do
      end if
      close (unit=ikey)
      return
      end
