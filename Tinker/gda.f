c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  program gda  --  simulated annealing on gaussian density  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "gda" implements Gaussian Density Annealing (GDA) algorithm
c     for global optimization via simulated annealing
c
c     literature reference:
c
c     J. Ma and J. E. Straub, "Simulated Annealing using the
c     Classical Density Distribution", Journal of Chemical Physics,
c     101, 533-541 (1994)
c
c
      program gda
      use sizes
      use atoms
      use files
      use iounit
      use minima
      use potent
      use vdwpot
      use warp
      implicit none
      integer i,igda,itrial,ntrial
      integer nstep,nvar,nok,nbad
      integer lext,next,freeunit
      real*8 bstart,bstop
      real*8 random,boxsize
      real*8 eps,h1,hmin,gda2
      real*8 minimum,grdmin
      real*8 xcm,ycm,zcm
      real*8, allocatable :: m2init(:)
      real*8, allocatable :: xx(:)
      logical exist,randomize
      character*1 answer
      character*6 mode,method
      character*7 ext,status
      character*240 gdafile
      character*240 record
      character*240 string
      external gda1,gda2,gda3
      external random,optsave
c
c
c     set up the structure, mechanics calculation and smoothing
c
      call initial
      call getxyz
      use_smooth = .true.
      use_gda = .true.
      call mechanic
c
c     get the number of optimized structures to be constructed
c
      ntrial = 0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=10,end=10)  ntrial
   10 continue
      if (ntrial .le. 0) then
         write (iout,20)
   20    format (/,' Enter Number of Annealing Trials [1] :  ',$)
         read (input,30)  ntrial
   30    format (i10)
      end if
      if (ntrial .le. 0)  ntrial = 1
c
c     see if random coordinates are desired as starting structures
c
      randomize = .true.
      if (ntrial .eq. 1) then
         randomize = .false.
         call nextarg (answer,exist)
         if (.not. exist) then
            write (iout,40)
   40       format (/,' Use Randomized Initial Coordinates [N] :  ',$)
            read (input,50)  record
   50       format (a240)
            next = 1
            call gettext (record,answer,next)
         end if
         call upcase (answer)
         if (answer .eq. 'Y')  randomize = .true.
      end if
      if (randomize)  boxsize = 10.0d0 * (dble(n))**(1.0d0/3.0d0)
c
c     get initial and final values of inverse temperature
c
      bstart = -1.0d0
      bstop = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=60,end=60)  bstart
      call nextarg (string,exist)
      if (exist)  read (string,*,err=60,end=60)  bstop
   60 continue
      if (bstart.le.0.0d0 .or. bstop.le.0.0d0) then
         write (iout,70)
   70    format (/,' Enter Initial and Final Beta [0.01, 10**10] :  ',$)
         read (input,80)  record
   80    format (a240)
         read (record,*,err=90,end=90)  bstart,bstop
   90    continue
      end if
      if (bstart .le. 0.0d0)  bstart = 0.01d0
      if (bstop .le. 0.0d0)  bstop = 1.0d10
c
c     perform dynamic allocation of some local arrays
c
      allocate (m2init(n))
      allocate (xx(4*n))
c
c     store the initial values of the squared mean Gaussian width
c
      do i = 1, n
         m2init(i) = m2(1)
      end do
c
c     write out a copy of coordinates for later update
c
      do itrial = 1, ntrial
         lext = 3
         call numeral (itrial,ext,lext)
         gdafile = filename(1:leng)//'.'//ext(1:lext)
         call version (gdafile,'new')
         igda = freeunit ()
         open (unit=igda,file=gdafile,status='new')
         call prtxyz (igda)
         close (unit=igda)
         outfile = gdafile
c
c     set an initial box size and generate random coordinates
c
         if (randomize) then
            do i = 1, n
               x(i) = boxsize * random ()
               y(i) = boxsize * random ()
               z(i) = boxsize * random ()
            end do
         end if
c
c     translate coordinates and M2's to optimization parameters
c
         nvar = 0
         do i = 1, n
            nvar = nvar + 1
            xx(nvar) = x(i)
            nvar = nvar + 1
            xx(nvar) = y(i)
            nvar = nvar + 1
            xx(nvar) = z(i)
         end do
         do i = 1, n
            nvar = nvar + 1
            xx(nvar) = m2init(i)
         end do
c
c     make changes to the potential to use potential smoothing
c
         use_smooth = .true.
         use_geom = .true.
         vdwtyp = 'GAUSSIAN'
c
c     make the call to the Bulirsch-Stoer integration routine
c
         nstep = 0
         status = '       '
         eps = 1.0d-8
         h1 = 0.01d0
         hmin = 0.0d0
         write (iout,100)
  100    format (//,' Gaussian Density Annealing Global Optimization :',
     &           //,' BS Step',5x,'Log(Beta)',6x,'Energy',
     &              9x,'Rg',8x,'Log(M2)',7x,'Status',/)
         call gdastat (nstep,bstart,xx,status)
         call diffeq (nvar,xx,bstart,bstop,eps,h1,hmin,nok,nbad,gda1)
         nstep = nok + nbad
c
c     make changes to the potential for standard optimization
c
         use_smooth = .false.
         use_geom = .false.
         vdwtyp = 'LENNARD-JONES'
c
c     make the call to the energy minimization routine
c
         mode = 'DTNCG'
         method = 'AUTO'
         nvar = 3 * n
         grdmin = 0.0001d0
         nextiter = nstep + 1
         call tncg (mode,method,nvar,xx,minimum,grdmin,
     &                     gda2,gda3,optsave)
c        call lbfgs (nvar,xx,minimum,grdmin,gda2,optsave)
         write (iout,110)  itrial,minimum
  110    format (/,' Global Energy Minimum for Trial',i4,' :',f15.4)
c
c     translate optimization parameters into atomic coordinates
c
         nvar = 0
         do i = 1, n
            nvar = nvar + 1
            x(i) = xx(nvar)
            nvar = nvar + 1
            y(i) = xx(nvar)
            nvar = nvar + 1
            z(i) = xx(nvar)
         end do
c
c     move the center of mass to the origin
c
         xcm = 0.0d0
         ycm = 0.0d0
         zcm = 0.0d0
         do i = 1, n
            xcm = xcm + x(i)
            ycm = ycm + y(i)
            zcm = zcm + z(i)
         end do
         xcm = xcm / dble(n)
         ycm = ycm / dble(n)
         zcm = zcm / dble(n)
         do i = 1, n
            x(i) = x(i) - xcm
            y(i) = y(i) - xcm
            z(i) = z(i) - xcm
         end do
c
c     write the final coordinates into a file
c
         igda = freeunit ()
         open (unit=igda,file=gdafile,status='old')
         rewind (unit=igda)
         call prtxyz (igda)
         close (unit=igda)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (m2init)
      deallocate (xx)
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine gda1  --  gaussian density annealing gradient  ##
c     ##                                                            ##
c     ################################################################
c
c
      subroutine gda1 (beta,xx,g)
      use sizes
      use atoms
      use iounit
      use warp
      implicit none
      integer i,nvar
      integer, allocatable :: hinit(:)
      integer, allocatable :: hstop(:)
      integer, allocatable :: hindex(:)
      real*8 beta,e,sum
      real*8, allocatable :: hdiag(:)
      real*8, allocatable :: h(:)
      real*8 xx(*)
      real*8 g(*)
      real*8, allocatable :: derivs(:,:)
c
c
c     translate optimization parameters to coordinates and M2's
c
      nvar = 0
      do i = 1, n
         nvar = nvar + 1
         x(i) = xx(nvar)
         nvar = nvar + 1
         y(i) = xx(nvar)
         nvar = nvar + 1
         z(i) = xx(nvar)
      end do
      do i = 1, n
         nvar = nvar + 1
         m2(i) = xx(nvar)
         if (m2(i) .lt. 0.0d0) then
            write (iout,10)  i,m2(i)
   10       format (' GDA1  --  Warning, Negative M2 at Atom',i6,
     &                  ' with Value',d12.4)
            m2(i) = -m2(i)
         end if
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (derivs(3,n))
c
c     compute and store the Cartesian energy gradient vector
c
      call gradient (e,derivs)
c
c     translate the energy gradient into a dr/dbeta vector
c
      nvar = 0
      do i = 1, n
         nvar = nvar + 1
         g(nvar) = -(m2(i)/3.0d0) * derivs(1,i)
         nvar = nvar + 1
         g(nvar) = -(m2(i)/3.0d0) * derivs(2,i)
         nvar = nvar + 1
         g(nvar) = -(m2(i)/3.0d0) * derivs(3,i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (derivs)
c
c     perform dynamic allocation of some local arrays
c
      allocate (hinit(3*n))
      allocate (hstop(3*n))
      allocate (hindex((3*n*(3*n-1))/2))
      allocate (hdiag(3*n))
      allocate (h((3*n*(3*n-1))/2))
c
c     compute and store the Hessian elements
c
      call hessian (h,hinit,hstop,hindex,hdiag)
c
c     translate the Hessian diagonal into a dM2/dbeta vector
c
      do i = 1, n
         nvar = nvar + 1
         sum = hdiag(3*i-2) + hdiag(3*i-1) + hdiag(3*i)
         g(nvar) = -(m2(i)/3.0d0)**2 * sum
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (hinit)
      deallocate (hstop)
      deallocate (hindex)
      deallocate (hdiag)
      deallocate (h)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  function gda2  --  energy/gradient for TNCG optimization  ##
c     ##                                                            ##
c     ################################################################
c
c
      function gda2 (xx,g)
      use sizes
      use atoms
      implicit none
      integer i,nvar
      real*8 gda2,e
      real*8 xx(*)
      real*8 g(*)
      real*8, allocatable :: derivs(:,:)
c
c
c     translate optimization parameters to atomic coordinates
c
      nvar = 0
      do i = 1, n
         nvar = nvar + 1
         x(i) = xx(nvar)
         nvar = nvar + 1
         y(i) = xx(nvar)
         nvar = nvar + 1
         z(i) = xx(nvar)
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (derivs(3,n))
c
c     compute and store the energy and gradient
c
      call gradient (e,derivs)
      gda2 = e
c
c     store Cartesian gradient as optimization gradient
c
      nvar = 0
      do i = 1, n
         nvar = nvar + 1
         g(nvar) = derivs(1,i)
         nvar = nvar + 1
         g(nvar) = derivs(2,i)
         nvar = nvar + 1
         g(nvar) = derivs(3,i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (derivs)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine gda3  --  Hessian values for TNCG optimization  ##
c     ##                                                             ##
c     #################################################################
c
c
      subroutine gda3 (mode,xx,h,hinit,hstop,hindex,hdiag)
      use sizes
      use atoms
      implicit none
      integer i,nvar
      integer hinit(*)
      integer hstop(*)
      integer hindex(*)
      real*8 xx(*)
      real*8 hdiag(*)
      real*8 h(*)
      character*4 mode
c
c
c     translate optimization parameters to atomic coordinates
c
      if (mode .eq. 'NONE')  return
      nvar = 0
      do i = 1, n
         nvar = nvar + 1
         x(i) = xx(nvar)
         nvar = nvar + 1
         y(i) = xx(nvar)
         nvar = nvar + 1
         z(i) = xx(nvar)
      end do
c
c     compute and store the Hessian elements
c
      call hessian (h,hinit,hstop,hindex,hdiag)
      return
      end
