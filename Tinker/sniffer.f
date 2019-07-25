c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  program sniffer  --  discrete generalized descent search  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "sniffer" performs a global energy minimization using a
c     discrete version of Griewank's global search trajectory
c
c     literature references:
c
c     A. O. Griewank, "Generalized Descent for Global Optimization",
c     Journal of Optimization Theory and Applications, 34, 11-39 (1981)
c
c     R. A. R. Butler and E. E. Slaminka, "An Evaluation of the Sniffer
c     Global Optimization Algorithm Using Standard Test Functions",
c     Journal of Computational Physics, 99, 28-32 (1992)
c
c     J. W. Rogers, Jr. and R. A. Donnelly, "Potential Transformation
c     Methods for Large-Scale Global Optimization", SIAM Journal of
c     Optimization, 5, 871-891 (1995)
c
c
      program sniffer
      use sizes
      use atoms
      use files
      use inform
      use iounit
      use linmin
      use math
      use minima
      use output
      use scales
      use usage
      implicit none
      integer i,j,imin
      integer nvar,niter
      integer start,stop
      integer freeunit
      integer istep,maxstep
      real*8 sniffer1,gnorm
      real*8 grms,grdmin
      real*8 f,eps,mu
      real*8 scaler,angle
      real*8 rms,size
      real*8 fmin,gmin
      real*8 dnorm,dot
      real*8 alpha,cosine
      real*8 epsfac,mufac
      real*8 stepfac
      real*8 minimum
      real*8, allocatable :: xx(:)
      real*8, allocatable :: g(:)
      real*8, allocatable :: d(:)
      real*8, allocatable :: derivs(:,:)
      logical exist,done
      character*9 status
      character*240 minfile
      character*240 string
c
c
c     set up the structure and mechanics calculation
c
      call initial
      call getxyz
      call mechanic
      call makeref (1)
c
c     get the number of steps in the initial block
c
      maxstep = -1
      call nextarg (string,exist)
      if (exist)  read (string,*,err=10,end=10)  maxstep
   10 continue
      if (maxstep .le. 0) then
         write (iout,20)
   20    format (/,' Enter Number of Steps in the Initial Set',
     &              ' [100] :  ',$)
         read (input,30)  maxstep
   30    format (i10)
      end if
      if (maxstep .le. 0)  maxstep = 100
c
c     get the target value for the global energy minimum
c
      fctmin = 1000000.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=40,end=40)  fctmin
   40 continue
      if (fctmin .ge. 1000000.0d0) then
         write (iout,50)
   50    format (/,' Enter Target Energy for the Global Minimum',
     &              ' [0.0] :  ',$)
         read (input,60)  fctmin
   60    format (f20.0)
      end if
      if (fctmin .ge. 1000000.0d0)  fctmin = 0.0d0
c
c     get termination criterion as RMS gradient per atom
c
      grdmin = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=70,end=70)  grdmin
   70 continue
      if (grdmin .le. 0.0d0) then
         write (iout,80)
   80    format (/,' Enter RMS Gradient per Atom Criterion [1.0] :  ',$)
         read (input,90)  grdmin
   90    format (f20.0)
      end if
      if (grdmin .le. 0.0d0)  grdmin = 1.0d0
c
c     write out a copy of coordinates for later update
c
      imin = freeunit ()
      minfile = filename(1:leng)//'.xyz'
      call version (minfile,'new')
      open (unit=imin,file=minfile,status='new')
      call prtxyz (imin)
      close (unit=imin)
      outfile = minfile
      coordtype = 'CARTESIAN'
c
c     set scaling parameter for function and derivative values;
c     use square root of median eigenvalue of a typical Hessian
c
      set_scale = .true.
      scaler = 1.0d0
      nvar = 0
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               nvar = nvar + 1
               scale(nvar) = scaler
            end do
         end if
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (xx(nvar))
      allocate (g(nvar))
      allocate (d(nvar))
      allocate (derivs(3,n))
c
c     scale the coordinates of each active atom
c
      nvar = 0
      do i = 1, n
         if (use(i)) then
            nvar = nvar + 1
            xx(nvar) = x(i) * scale(nvar)
            nvar = nvar + 1
            xx(nvar) = y(i) * scale(nvar)
            nvar = nvar + 1
            xx(nvar) = z(i) * scale(nvar)
         end if
      end do
c
c     set initial values for the control parameters
c
      epsfac = 1.1d0
      mufac = 1.7d0
      stepfac = 1.1d0
c
c     set initial values for optimization parameters
c
      iprint = 1
      iwrite = 100
      rms = sqrt(dble(n))
      start = 0
      stop = start + maxstep
      eps = 1.0d0
      mu = 1.0d0
      stpmax = 0.1d0 * rms
      stpmin = 0.001d0
c
c     initialize unit direction vector along negative gradient
c
      f = sniffer1 (xx,g)
      gnorm = 0.0d0
      do i = 1, nvar
         gnorm = gnorm + g(i)**2
      end do
      gnorm = sqrt(gnorm)
      grms = gnorm / rms
      do i = 1, nvar
         d(i) = -g(i) / gnorm
      end do
      fmin = f
      gmin = grms
c
c     tests of the successful termination criteria
c
      if (fmin .le. fctmin) then
         status = 'TargetVal'
         done = .true.
      else if (gmin .le. grdmin) then
         status = 'SmallGrad'
         done = .true.
      else
         done = .false.
      end if
c
c     print header information prior to iterations
c
      if (iprint .gt. 0) then
         write (iout,100)
  100    format (/,' Discrete Generalized Descent Global',
     &              ' Optimization :')
      end if
c
c     perform a set of basic sniffer search steps
c
      niter = 0
      do while (.not. done)
         write (iout,110)
  110    format (/,4x,'Iter',11x,'F Value',13x,'G RMS',
     &             8x,'X Move',9x,'Angle',/)
         do istep = start, stop
c
c     get the current energy and gradient values
c
            f = sniffer1 (xx,g)
c
c     if current energy is lowest yet, save the coordinates
c
            if (f .lt. fmin) then
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
               call makeref (1)
               imin = freeunit ()
               open (unit=imin,file=minfile,status='old')
               rewind (unit=imin)
               call prtxyz (imin)
               close (unit=imin)
            end if
c
c     get rms gradient and dot product with search direction
c
            gnorm = 0.0d0
            dot = 0.0d0
            do i = 1, nvar
               gnorm = gnorm + g(i)*g(i)
               dot = dot + d(i)*g(i)
            end do
            gnorm = sqrt(gnorm)
            grms = gnorm / (scaler*rms)
c
c     compute the next direction vector and its length
c
            alpha = max(0.0d0,1.0d0+(1.0d0+eps)*dot)
            dnorm = 0.0d0
            do i = 1, nvar
               d(i) = -eps*g(i) + alpha*d(i)
               dnorm = dnorm + d(i)*d(i)
            end do
            dnorm = sqrt(dnorm)
c
c     normalize direction and get angle with negative gradient
c
            dot = 0.0d0
            do i = 1, nvar
               d(i) = d(i) / dnorm
               dot = dot + d(i)*g(i)
            end do
            cosine = -dot / gnorm
            cosine = min(1.0d0,max(-1.0d0,cosine))
            angle = radian * acos(cosine)
c
c     move atomic positions along the direction vector
c
            size = min(stpmax,mu*(f-fctmin))
            do i = 1, nvar
               xx(i) = xx(i) + size*d(i)
            end do
c
c     compute the size of the step taken
c
            size = min(stpmax,mu*(f-fctmin))
            size = size / rms
c
c     update the best value and gradient found so far
c
            fmin = min(fmin,f)
            gmin = min(gmin,grms)
c
c     print intermediate results every few iterations
c
            if (iprint .gt. 0) then
               if (done .or. mod(niter,iprint).eq.0) then
                  if (f.lt.1.0d12 .and. f.gt.-1.0d11
     &                .and. grms.lt.1.0d12) then
                     write (iout,120)  istep,f,grms,size,angle
  120                format (i8,2f18.4,2f14.4)
                  else
                     write (iout,130)  istep,f,grms,size,angle
  130                format (i8,2d18.4,2f14.4)
                  end if
               end if
            end if
         end do
c
c     tests of the various termination and error criteria
c
         if (fmin .le. fctmin) then
            status = 'TargetVal'
            done = .true.
         else if (gmin .le. grdmin) then
            status = 'SmallGrad'
            done = .true.
         else if (size .le. stpmin) then
            status = 'SmallMove'
            done = .true.
         end if
c
c     write the final coordinates for this set of steps
c
         niter = niter + 1
         if (cyclesave)  call optsave (niter,fmin,xx)
c
c     update the optimization parameters for the next set
c
         eps = eps * epsfac
         mu = mu / mufac
         maxstep = nint(dble(maxstep)*stepfac)
         start = stop + 1
         stop = start + maxstep
      end do
c
c     write message about satisfaction of termination criteria
c
      if (status .eq. 'SmallMove') then
         write (iout,140)  status
  140    format (/,' SNIFFER  --  Incomplete Convergence due to ',a9)
      else
         write (iout,150)  status
  150    format (/,' SNIFFER  --  Normal Termination due to ',a9)
      end if
c
c     use lowest energy structure as global minimum estimate
c
      call getref (1)
c
c     compute the final function and RMS gradient values
c
      call gradient (minimum,derivs)
      gnorm = 0.0d0
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               gnorm = gnorm + derivs(j,i)**2
            end do
         end if
      end do
      gnorm = sqrt(gnorm)
      grms = gnorm / rms
c
c     perform deallocation of some local arrays
c
      deallocate (xx)
      deallocate (g)
      deallocate (d)
      deallocate (derivs)
c
c     write out the final function and gradient values
c
      if (grms .gt. 0.0001d0) then
         write (iout,160)  minimum,grms,gnorm
  160    format (/,' Final Function Value :',f15.4,
     &           /,' Final RMS Gradient :  ',f15.4,
     &           /,' Final Gradient Norm : ',f15.4)
      else
         write (iout,170)  minimum,grms,gnorm
  170    format (/,' Final Function Value :',f15.4,
     &           /,' Final RMS Gradient :  ',d15.4,
     &           /,' Final Gradient Norm : ',d15.4)
      end if
c
c     write the final coordinates into a file
c
      imin = freeunit ()
      open (unit=imin,file=minfile,status='old')
      rewind (unit=imin)
      call prtxyz (imin)
      close (unit=imin)
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  function sniffer1  --  energy and gradient for sniffer  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "sniffer1" is a service routine that computes the energy
c     and gradient for the Sniffer global optimization method
c
c
      function sniffer1 (xx,g)
      use sizes
      use atoms
      use scales
      use usage
      implicit none
      integer i,nvar
      real*8 sniffer1,e
      real*8 xx(*)
      real*8 g(*)
      real*8, allocatable :: derivs(:,:)
c
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
      call gradient (e,derivs)
      sniffer1 = e
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
