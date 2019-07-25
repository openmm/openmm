c
c
c     ##############################################################
c     ##  COPYRIGHT (C) 1997 by Rohit Pappu & Jay William Ponder  ##
c     ##                   All Rights Reserved                    ##
c     ##############################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  program pss  --  Cartesian potential smoothing & search  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "pss" implements the potential smoothing plus search method
c     for global optimization in Cartesian coordinate space with
c     local searches performed in Cartesian or torsional space
c
c     literature reference:
c
c     J. Kostrowicki and H. A. Scheraga, "Application of the Diffusion
c     Equation Method for Global Optimization to Oligopeptides", Journal
c     of Physical Chemistry, 96, 7442-7449 (1992)
c
c     S. Nakamura, H. Hirose, M. Ikeguchi and J. Doi, "Conformational
c     Energy Minimization Using a Two-Stage Method", Journal of Physical
c     Chemistry, 99, 8374-8378 (1995)
c
c
      program pss
      use sizes
      use atoms
      use inform
      use iounit
      use omega
      use refer
      use tree
      use warp
      implicit none
      integer i,next,range
      integer start,stop
      real*8 minimum,grdmin
      real*8 srchmax,rms
      real*8 ratio,sigmoid
      logical exist,check
      logical use_forward
      logical use_cart
      logical use_tors
      character*1 answer
      character*1 formtyp
      character*240 record
      character*240 string
c
c
c     set up the structure, mechanics calculation and smoothing
c
      call initial
      call getxyz
      use_smooth = .true.
      use_dem = .true.
      call mechanic
      iwrite = 0
c
c     get the number of points along the deformation schedule
c
      nlevel = -1
      call nextarg (string,exist)
      if (exist)  read (string,*,err=10,end=10)  nlevel
   10 continue
      if (nlevel .lt. 0) then
         write (iout,20)
   20    format (/,' Enter the Number of Steps for Smoothing Schedule',
     &              ' [100] :  ',$)
         read (input,30)  nlevel
   30    format (i10)
         if (nlevel .le. 0)  nlevel = 100
      end if
c
c     decide whether to use forward smoothing of initial structure
c
      use_forward = .true.
      call nextarg (answer,exist)
      if (.not. exist) then
         write (iout,40)
   40    format (/,' Perform Forward Smoothing from Input Structure',
     &              ' [Y] :  ',$)
         read (input,50)  record
   50    format (a240)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .eq. 'N')  use_forward = .false.
c
c     get the functional form for the deformation schedule
c
      formtyp = 'C'
      call nextarg (answer,exist)
      if (.not. exist) then
         write (iout,60)
   60    format (/,' Use Quadratic, Cubic or Sigmoidal Schedule',
     &              ' (Q [C] or S) :  ',$)
         read (input,70)  record
   70    format (a240)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .eq. 'Q')  formtyp = answer
      if (answer .eq. 'S')  formtyp = answer
c
c     decide which type of local search procedure to use
c
      use_cart = .false.
      use_tors = .false.
      call nextarg (answer,exist)
      if (.not. exist) then
         write (iout,80)
   80    format (/,' Local Search Type - Cartesian, Torsional or None',
     &              ' (C T or [N]) :  ',$)
         read (input,90)  record
   90    format (a240)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .eq. 'C')  use_cart = .true.
      if (answer .eq. 'T')  use_tors = .true.
c
c     get the rotatable bonds for torsional local search
c
      if (use_tors) then
         call makeint (0)
         call initrot
         call active
      end if
c
c     get the number of eigenvectors to use for local search
c
      if (use_cart .or. use_tors) then
         start = -1
         stop = -1
         call nextarg (string,exist)
         if (exist)  read (string,*,err=100,end=100)  start
         call nextarg (string,exist)
         if (exist)  read (string,*,err=100,end=100)  stop
  100    continue
         if (stop .le. 0) then
            write (iout,110)
  110       format (/,' Enter the Range of Local Search Directions',
     &                 ' (1=Highest Freq) :  ',$)
            read (input,120)  record
  120       format (a240)
            read (record,*)  start,stop
            range = abs(stop-start)
            start = min(start,stop)
            stop = start + range
         end if
         if (use_cart)  stop = min(stop,3*n-6)
         if (use_tors)  stop = min(stop,nomega)
      end if
c
c     get the maximal smoothing level for use of local search
c
      if (use_cart .or. use_tors) then
         srchmax = -1.0d0
         call nextarg (string,exist)
         if (exist)  read (string,*,err=130,end=130)  srchmax
  130    continue
         if (srchmax .lt. 0.0d0) then
            write (iout,140)
  140       format (/,' Enter the Largest Smoothing Level for',
     &                 ' Local Search [5.0] :  ',$)
            read (input,150)  srchmax
  150       format (f20.0)
            if (srchmax .lt. 0.0d0)  srchmax = 5.0d0
         end if
      end if
c
c     decide whether to use forward smoothing of initial structure
c
      check = .false.
      if ((use_cart .or. use_tors) .and. .not.use_forward) then
         call nextarg (answer,exist)
         if (.not. exist) then
            write (iout,160)
  160       format (/,' Restrict Local Search to Children of Input',
     &                 ' Structure [N] :  ',$)
            read (input,170)  record
  170       format (a240)
            next = 1
            call gettext (record,answer,next)
         end if
         call upcase (answer)
         if (answer .eq. 'Y')  check = .true.
      end if
c
c     get the termination criterion as RMS gradient per atom
c
      grdmin = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=180,end=180)  grdmin
  180 continue
      if (grdmin .le. 0.0d0) then
         write (iout,190)
  190    format (/,' Enter RMS Gradient per Atom Criterion',
     &              ' [0.0001] :  ',$)
         read (input,200)  grdmin
  200    format (f20.0)
      end if
      if (grdmin .le. 0.0d0)  grdmin = 0.0001d0
c
c     compute the smoothing levels for the desired protocol
c
      do i = 0, nlevel
         ratio = 1.0d0 - dble(nlevel-i)/dble(nlevel)
         if (formtyp .eq. 'Q') then
            ilevel(i) = deform * ratio**2
         else if (formtyp .eq. 'C') then
            ilevel(i) = deform * ratio**3
         else if (formtyp .eq. 'S') then
            ilevel(i) = deform * sigmoid (12.0d0,ratio)
         end if
      end do
c
c     perform forward PSS by looping over smoothed surfaces
c
      if (use_forward) then
         do i = 0, nlevel-1
            deform = ilevel(i)
            call makeref (1)
            iprint = 1
            call localxyz (minimum,grdmin)
            call impose (n,xref,yref,zref,n,x,y,z,rms)
            call psswrite (i)
            write (iout,210)  minimum,deform
  210       format (/,' Final Function Value and Deformation :',2f15.4)
         end do
      end if
c
c     perform PSS reversal by looping over smoothed surfaces
c
      do i = nlevel, 0, -1
         deform = ilevel(i)
         call makeref (1)
         iprint = 1
         call localxyz (minimum,grdmin)
         call impose (n,xref,yref,zref,n,x,y,z,rms)
         if (i .eq. nlevel)  etree = minimum
         if (deform .le. srchmax) then
            if (use_cart) then
               call modecart (start,stop,minimum,grdmin,check)
            else if (use_tors) then
               call modetors (start,stop,minimum,grdmin,check)
            end if
         end if
         if (use_forward) then
            call psswrite (2*nlevel-i)
         else
            call psswrite (nlevel-i)
         end if
         write (iout,220)  minimum,deform
  220    format (/,' Final Function Value and Deformation :',2f15.4)
      end do
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  function pss1  --  energy and gradient values for PSS  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "pss1" is a service routine that computes the energy
c     and gradient during PSS global optimization in Cartesian
c     coordinate space
c
c
      function pss1 (xx,g)
      use sizes
      use atoms
      implicit none
      integer i,nvar
      real*8 pss1,e
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
      pss1 = e
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
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine pss2  --  Hessian matrix values for PSS  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "pss2" is a service routine that computes the sparse
c     matrix Hessian elements during PSS global optimization
c     in Cartesian coordinate space
c
c
      subroutine pss2 (mode,xx,h,hinit,hstop,hindex,hdiag)
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
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine modecart  --  Cartesian local search for PSS  ##
c     ##                                                           ##
c     ###############################################################
c
c
      subroutine modecart (start,stop,minimum,grdmin,check)
      use sizes
      use atoms
      use iounit
      use omega
      use refer
      implicit none
      integer i,j,k,nfreq
      integer start,stop
      integer niter,nsearch
      real*8 minimum,grdmin
      real*8 minref,minbest
      real*8 eps,rms,size
      real*8, allocatable :: xbest(:)
      real*8, allocatable :: ybest(:)
      real*8, allocatable :: zbest(:)
      real*8, allocatable :: eigen(:)
      real*8, allocatable :: step(:,:)
      real*8, allocatable :: vects(:,:)
      logical done,check
c
c
c     store the current coordinates as the reference set
c
      call makeref (1)
c
c     set parameters related to the local search procedure
c
      done = .false.
      eps = 0.0001d0
      minref = minimum
      minbest = minimum
      niter = 0
c
c     perform dynamic allocation of some local arrays
c
      nfreq = 3 * n
      allocate (xbest(n))
      allocate (ybest(n))
      allocate (zbest(n))
      allocate (eigen(nfreq))
      allocate (step(3,nfreq))
      allocate (vects(nfreq,nfreq))
c
c     find local minimum along each of the steepest directions
c
      do while (.not. done)
         niter = niter + 1
         write (iout,10)  niter,minref
   10    format (/,' Cartesian Mode Search :',5x,'Iteration',i4,
     &              6x,'Energy',f12.4,/)
         call eigenxyz (eigen,vects)
c
c     search both directions along each eigenvector in turn
c
         nsearch = 0
         do i = start, stop
            do k = 1, n
               j = 3*(k-1)
               size = 1.0d0 / sqrt(abs(eigen(3*n-i+1)))
               step(1,k) = size * vects(j+1,3*n-i+1)
               step(2,k) = size * vects(j+2,3*n-i+1)
               step(3,k) = size * vects(j+3,3*n-i+1)
            end do
            nsearch = nsearch + 1
            call getref (1)
            call climbxyz (nsearch,minimum,step,grdmin,check)
            if (minimum .lt. minbest) then
               minbest = minimum
               do k = 1, n
                  xbest(k) = x(k)
                  ybest(k) = y(k)
                  zbest(k) = z(k)
               end do
            end if
            do k = 1, n
               step(1,k) = -step(1,k)
               step(2,k) = -step(2,k)
               step(3,k) = -step(3,k)
            end do
            nsearch = nsearch + 1
            call getref (1)
            call climbxyz (nsearch,minimum,step,grdmin,check)
            if (minimum .lt. minbest) then
               minbest = minimum
               do k = 1, n
                  xbest(k) = x(k)
                  ybest(k) = y(k)
                  zbest(k) = z(k)
               end do
            end if
         end do
c
c     check for convergence of the local search procedure
c
         if (minbest .lt. minref-eps) then
            done = .false.
            minref = minbest
            call impose (n,xref,yref,zref,n,xbest,ybest,zbest,rms)
            do k = 1, n
               x(k) = xbest(k)
               y(k) = ybest(k)
               z(k) = zbest(k)
            end do
            call makeref (1)
         else
            done = .true.
            minimum = minref
            call getref (1)
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xbest)
      deallocate (ybest)
      deallocate (zbest)
      deallocate (eigen)
      deallocate (step)
      deallocate (vects)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine modetors  --  torsional local search for PSS  ##
c     ##                                                           ##
c     ###############################################################
c
c
      subroutine modetors (start,stop,minimum,grdmin,check)
      use sizes
      use atoms
      use iounit
      use omega
      use refer
      implicit none
      integer i,k
      integer start,stop
      integer niter,nsearch
      real*8 minimum,grdmin
      real*8 minref,minbest
      real*8 eps,rms
      real*8, allocatable :: xbest(:)
      real*8, allocatable :: ybest(:)
      real*8, allocatable :: zbest(:)
      real*8, allocatable :: step(:)
      real*8, allocatable :: eigen(:)
      real*8, allocatable :: vects(:,:)
      logical done,check
c
c
c     store the current coordinates as the reference set
c
      call makeref (1)
c
c     set parameters related to the local search procedure
c
      done = .false.
      eps = 0.0001d0
      minref = minimum
      minbest = minimum
      niter = 0
c
c     perform dynamic allocation of some local arrays
c
      allocate (xbest(n))
      allocate (ybest(n))
      allocate (zbest(n))
      allocate (step(nomega))
      allocate (eigen(nomega))
      allocate (vects(nomega,nomega))
c
c     find local minimum along each of the steepest directions
c
      do while (.not. done)
         niter = niter + 1
         write (iout,10)  niter,minref
   10    format (/,' Torsional Mode Search :',5x,'Iteration',i4,
     &              6x,'Energy',f12.4,/)
         call makeint (0)
         call eigentor (eigen,vects)
c
c     search both directions along each eigenvector in turn
c
         nsearch = 0
         do i = start, stop
            do k = 1, nomega
               step(k) = vects(k,nomega-i+1)
            end do
            nsearch = nsearch + 1
            call climbtor (nsearch,minimum,step,grdmin,check)
            if (minimum .lt. minbest) then
               minbest = minimum
               do k = 1, n
                  xbest(k) = x(k)
                  ybest(k) = y(k)
                  zbest(k) = z(k)
               end do
            end if
            do k = 1, nomega
               step(k) = -step(k)
            end do
            nsearch = nsearch + 1
            call climbtor (nsearch,minimum,step,grdmin,check)
            if (minimum .lt. minbest) then
               minbest = minimum
               do k = 1, n
                  xbest(k) = x(k)
                  ybest(k) = y(k)
                  zbest(k) = z(k)
               end do
            end if
         end do
c
c     check for convergence of the local search procedure
c
         if (minbest .lt. minref-eps) then
            done = .false.
            minref = minbest
            call impose (n,xref,yref,zref,n,xbest,ybest,zbest,rms)
            do k = 1, n
               x(k) = xbest(k)
               y(k) = ybest(k)
               z(k) = zbest(k)
            end do
            call makeref (1)
         else
            done = .true.
            minimum = minref
            call getref (1)
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xbest)
      deallocate (ybest)
      deallocate (zbest)
      deallocate (step)
      deallocate (eigen)
      deallocate (vects)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine eigenxyz  --  Cartesian Hessian eigenvectors  ##
c     ##                                                           ##
c     ###############################################################
c
c
      subroutine eigenxyz (eigen,vects)
      use sizes
      use atoms
      use hescut
      implicit none
      integer i,j,k,nfreq,ihess
      integer, allocatable :: hindex(:)
      integer, allocatable :: hinit(:,:)
      integer, allocatable :: hstop(:,:)
      real*8 eigen(*)
      real*8, allocatable :: matrix(:)
      real*8, allocatable :: h(:)
      real*8 vects(3*n,*)
      real*8, allocatable :: hdiag(:,:)
c
c
c     perform dynamic allocation of some local arrays
c
      nfreq = 3 * n
      allocate (hindex((nfreq*(nfreq-1))/2))
      allocate (hinit(3,n))
      allocate (hstop(3,n))
      allocate (matrix((nfreq*(nfreq+1))/2))
      allocate (h((nfreq*(nfreq-1))/2))
      allocate (hdiag(3,n))
c
c     compute the Hessian matrix in Cartesian space
c
      hesscut = 0.0d0
      call hessian (h,hinit,hstop,hindex,hdiag)
c
c     place Hessian elements into triangular form
c
      ihess = 0
      do i = 1, n
         do j = 1, 3
            ihess = ihess + 1
            matrix(ihess) = hdiag(j,i)
            do k = hinit(j,i), hstop(j,i)
               ihess = ihess + 1
               matrix(ihess) = h(k)
            end do
         end do
      end do
c
c     diagonalize the Hessian to obtain eigenvalues
c
      call diagq (nfreq,nfreq,matrix,eigen,vects)
c
c     perform deallocation of some local arrays
c
      deallocate (hindex)
      deallocate (hinit)
      deallocate (hstop)
      deallocate (matrix)
      deallocate (h)
      deallocate (hdiag)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine eigentor  --  torsional Hessian eigenvectors  ##
c     ##                                                           ##
c     ###############################################################
c
c
      subroutine eigentor (eigen,vects)
      use sizes
      use atoms
      use omega
      implicit none
      integer i,j,ihess
      real*8 eigen(*)
      real*8, allocatable :: matrix(:)
      real*8 vects(nomega,*)
      real*8, allocatable :: hrot(:,:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (matrix(nomega*(nomega+1)/2))
      allocate (hrot(nomega,nomega))
c
c     compute the Hessian in torsional space
c
      call hessrot ('FULL',hrot)
c
c     place Hessian elements into triangular form
c
      ihess = 0
      do i = 1, nomega
         do j = i, nomega
            ihess = ihess + 1
            matrix(ihess) = hrot(i,j)
         end do
      end do
c
c     diagonalize the Hessian to obtain eigenvalues
c
      call diagq (nomega,nomega,matrix,eigen,vects)
c
c     perform deallocation of some local arrays
c
      deallocate (matrix)
      deallocate (hrot)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine climbxyz  --  Cartesian local search direction  ##
c     ##                                                             ##
c     #################################################################
c
c
      subroutine climbxyz (nsearch,minimum,step,grdmin,check)
      use sizes
      use atoms
      use inform
      use iounit
      use refer
      implicit none
      integer maxstep
      parameter (maxstep=500)
      integer i,kstep
      integer nstep,nsearch
      real*8 minimum,grdmin
      real*8 parent
      real*8 energy,big
      real*8 step(3,*)
      real*8 estep(0:maxstep)
      logical done,check,keep
c
c
c     convert current reference coordinates to a Z-matrix
c
      call getref (1)
c
c     set the maximum number of steps and the step size
c
      done = .false.
      keep = .true.
      iprint = 0
      big = 1000000.0d0
      minimum = big
      kstep = 0
      nstep = 65
c
c     scan the search direction for a minimization candidate
c
      do while (.not. done)
         if (kstep .ne. 0) then
            do i = 1, n
               x(i) = x(i) + step(1,i)
               y(i) = y(i) + step(2,i)
               z(i) = z(i) + step(3,i)
            end do
         end if
         estep(kstep) = energy ()
         if (kstep .ge. 2) then
            if (estep(kstep) .lt. estep(kstep-2) .and.
     &          estep(kstep-1) .lt. estep(kstep-2)) then
               done = .true.
               do i = 1, n
                  x(i) = x(i) - step(1,i)
                  y(i) = y(i) - step(2,i)
                  z(i) = z(i) - step(3,i)
               end do
               call localxyz (minimum,grdmin)
               parent = minimum
               if (check)  call chktree (parent,grdmin,keep)
               if (minimum .ge. -big) then
                  if (check) then
                     write (iout,10)  nsearch,kstep-1,minimum,parent
   10                format (4x,'Search Direction',i4,10x,'Step',
     &                          i6,10x,2f12.4)
                  else
                     write (iout,20)  nsearch,kstep-1,minimum
   20                format (4x,'Search Direction',i4,10x,'Step',
     &                          i6,10x,f12.4)
                  end if
               else
                  minimum = big
                  write (iout,30)  nsearch
   30             format (4x,'Search Direction',i4,36x,'------')
               end if
               if (.not. keep)  minimum = big
            end if
         end if
         if (kstep.ge.nstep .and. .not.done) then
            done = .true.
            write (iout,40)  nsearch
   40       format (4x,'Search Direction',i4,36x,'------')
         end if
         kstep = kstep + 1
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine climbtor  --  torsional local search direction  ##
c     ##                                                             ##
c     #################################################################
c
c
      subroutine climbtor (nsearch,minimum,step,grdmin,check)
      use sizes
      use inform
      use iounit
      use math
      use omega
      use zcoord
      implicit none
      integer maxstep
      parameter (maxstep=500)
      integer i,kstep
      integer nstep,nsearch
      real*8 minimum,grdmin
      real*8 parent
      real*8 energy,big,size
      real*8 step(*)
      real*8 estep(0:maxstep)
      logical done,check,keep
c
c
c     convert current reference coordinates to a Z-matrix
c
      call getref (1)
      call makeint (0)
c
c     set the maximum number of steps and the step size
c
      done = .false.
      keep = .true.
      iprint = 0
      big = 1000000.0d0
      minimum = big
      kstep = 0
      nstep = 65
      size = 0.1d0 * radian
      do i = 1, nomega
         step(i) = size * step(i)
      end do
c
c     scan the search direction for a minimization candidate
c
      do while (.not. done)
         if (kstep .ne. 0) then
            do i = 1, nomega
               ztors(zline(i)) = ztors(zline(i)) + step(i)
            end do
         end if
         call makexyz
         estep(kstep) = energy ()
         if (kstep .ge. 2) then
            if (estep(kstep) .lt. estep(kstep-2) .and.
     &          estep(kstep-1) .lt. estep(kstep-2)) then
               done = .true.
               do i = 1, nomega
                  ztors(zline(i)) = ztors(zline(i)) - step(i)
               end do
               call makexyz
               call localxyz (minimum,grdmin)
               parent = minimum
               if (check)  call chktree (parent,grdmin,keep)
               if (minimum .ge. -big) then
                  if (check) then
                     write (iout,10)  nsearch,kstep-1,minimum,parent
   10                format (4x,'Search Direction',i4,10x,'Step',
     &                          i6,10x,2f12.4)
                  else
                     write (iout,20)  nsearch,kstep-1,minimum
   20                format (4x,'Search Direction',i4,10x,'Step',
     &                          i6,10x,f12.4)
                  end if
               else
                  minimum = big
                  write (iout,30)  nsearch
   30             format (4x,'Search Direction',i4,36x,'------')
               end if
               if (.not. keep)  minimum = big
            end if
         end if
         if (kstep.ge.nstep .and. .not.done) then
            done = .true.
            write (iout,40)  nsearch
   40       format (4x,'Search Direction',i4,36x,'------')
         end if
         kstep = kstep + 1
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine localxyz  --  PSS local search optimization  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "localxyz" is used during the potential smoothing and search
c     procedure to perform a local optimization at the current
c     smoothing level
c
c
      subroutine localxyz (minimum,grdmin)
      use sizes
      use atoms
      use inform
      implicit none
      integer i,nvar
      real*8 minimum
      real*8 grdmin
      real*8 pss1
      real*8, allocatable :: xx(:)
      logical oldverb
      character*6 mode
      character*6 method
      external pss1,pss2
      external optsave
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (xx(3*n))
c
c     translate the coordinates of each atom
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
c
c     make the call to the optimization routine
c
      oldverb = verbose
      verbose = .false.
      mode = 'AUTO'
      method = 'AUTO'
      call tncg (mode,method,nvar,xx,minimum,grdmin,
     &                  pss1,pss2,optsave)
      verbose = oldverb
c
c     untranslate the final coordinates for each atom
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
c     perform deallocation of some local arrays
c
      deallocate (xx)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine chktree  --  check for legitimacy of branch  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "chktree" tests a minimum energy structure to see if it
c     belongs to the correct progenitor in the existing map
c
c
      subroutine chktree (parent,grdmin,keep)
      use sizes
      use atoms
      use tree
      use warp
      implicit none
      integer i
      real*8 parent,grdmin
      real*8 deform0,eps
      real*8, allocatable :: x0(:)
      real*8, allocatable :: y0(:)
      real*8, allocatable :: z0(:)
      logical keep
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (x0(n))
      allocate (y0(n))
      allocate (z0(n))
c
c     store the current smoothing level and coordinates
c
      deform0 = deform
      do i = 1, n
         x0(i) = x(i)
         y0(i) = y(i)
         z0(i) = z(i)
      end do
c
c     forward smoothing optimizations back to highest level
c
      do i = 1, nlevel
         if (deform .lt. ilevel(i)) then
            deform = ilevel(i)
            call localxyz (parent,grdmin)
         end if
      end do
c
c     compare energy to reference value for this tree branch
c
      eps = 1.0d-4
      keep = .false.
      if (abs(parent-etree) .lt. eps)  keep = .true.
c
c     restore the original smoothing level and coordinates
c
      deform = deform0
      do i = 1, n
         x(i) = x0(i)
         y(i) = y0(i)
         z(i) = z0(i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (x0)
      deallocate (y0)
      deallocate (z0)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine psswrite  --  output structures on PSS path  ##
c     ##                                                          ##
c     ##############################################################
c
c
      subroutine psswrite (i)
      use files
      implicit none
      integer i,ixyz
      integer lext,freeunit
      character*7 ext
      character*240 xyzfile
c
c
c     write the coordinates of the current minimum to a file
c
      lext = 3
      call numeral (i,ext,lext)
      ixyz = freeunit ()
      xyzfile = filename(1:leng)//'.'//ext(1:lext)
      call version (xyzfile,'new')
      open (unit=ixyz,file=xyzfile,status='new')
      call prtxyz (ixyz)
      close (unit=ixyz)
      return
      end
