c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine lbfgs  --  limited memory BFGS optimization  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "lbfgs" is a limited memory BFGS quasi-newton nonlinear
c     optimization routine
c
c     literature references:
c
c     J. Nocedal, "Updating Quasi-Newton Matrices with Limited
c     Storage", Mathematics of Computation, 35, 773-782 (1980)
c
c     D. C. Lui and J. Nocedal, "On the Limited Memory BFGS Method
c     for Large Scale Optimization", Mathematical Programming,
c     45, 503-528 (1989)
c
c     J. Nocedal and S. J. Wright, "Numerical Optimization",
c     Springer-Verlag, New York, 1999, Section 9.1
c
c     variables and parameters:
c
c     nvar      number of parameters in the objective function
c     x0        contains starting point upon input, upon return
c                 contains the best point found
c     minimum   during optimization contains best current function
c                 value; returns final best function value
c     grdmin    normal exit if rms gradient gets below this value
c     ncalls    total number of function/gradient evaluations
c
c     required external routines:
c
c     fgvalue    function to evaluate function and gradient values
c     optsave    subroutine to write out info about current status
c
c
      subroutine lbfgs (nvar,x0,minimum,grdmin,fgvalue,optsave)
      use inform
      use iounit
      use keys
      use linmin
      use math
      use minima
      use output
      use scales
      implicit none
      integer i,j,k,m
      integer nvar,next
      integer msav,muse
      integer niter,ncalls
      integer nerr,maxerr
      real*8 f,f_old,fgvalue
      real*8 f_move,x_move
      real*8 g_norm,g_rms
      real*8 minimum,grdmin
      real*8 angle,rms,beta
      real*8 ys,yy,gamma
      real*8 x0(*)
      real*8, allocatable :: rho(:)
      real*8, allocatable :: alpha(:)
      real*8, allocatable :: x_old(:)
      real*8, allocatable :: g(:)
      real*8, allocatable :: g_old(:)
      real*8, allocatable :: p(:)
      real*8, allocatable :: q(:)
      real*8, allocatable :: r(:)
      real*8, allocatable :: h0(:)
      real*8, allocatable :: s(:,:)
      real*8, allocatable :: y(:,:)
      logical done
      character*9 blank,status
      character*20 keyword
      character*240 record
      character*240 string
      external fgvalue,optsave
c
c
c     initialize some values to be used below
c
      ncalls = 0
      rms = sqrt(dble(nvar))
      if (coordtype .eq. 'CARTESIAN') then
         rms = rms / sqrt(3.0d0)
      else if (coordtype .eq. 'RIGIDBODY') then
         rms = rms / sqrt(6.0d0)
      end if
      blank = '         '
      done = .false.
      nerr = 0
      maxerr = 2
c
c     set default values for variable scale factors
c
      if (.not. set_scale) then
         do i = 1, nvar
            if (scale(i) .eq. 0.0d0)  scale(i) = 1.0d0
         end do
      end if
c
c     set default parameters for the optimization
c
      msav = min(nvar,20)
      if (fctmin .eq. 0.0d0)  fctmin = -100000000.0d0
      if (maxiter .eq. 0)  maxiter = 1000000
      if (nextiter .eq. 0)  nextiter = 1
      if (iprint .lt. 0)  iprint = 1
      if (iwrite .lt. 0)  iwrite = 1
c
c     set default parameters for the line search
c
      if (stpmax .eq. 0.0d0)  stpmax = 5.0d0
      stpmin = 1.0d-16
      cappa = 0.9d0
      slpmax = 10000.0d0
      angmax = 180.0d0
      intmax = 5
c
c     search the keywords for optimization parameters
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:14) .eq. 'LBFGS-VECTORS ') then
            read (string,*,err=10,end=10)  msav
            msav = max(0,min(msav,nvar))
         else if (keyword(1:17) .eq. 'STEEPEST-DESCENT ') then
            msav = 0
         else if (keyword(1:7) .eq. 'FCTMIN ') then
            read (string,*,err=10,end=10)  fctmin
         else if (keyword(1:8) .eq. 'MAXITER ') then
            read (string,*,err=10,end=10)  maxiter
         else if (keyword(1:8) .eq. 'STEPMAX ') then
            read (string,*,err=10,end=10)  stpmax
         else if (keyword(1:8) .eq. 'STEPMIN ') then
            read (string,*,err=10,end=10)  stpmin
         else if (keyword(1:6) .eq. 'CAPPA ') then
            read (string,*,err=10,end=10)  cappa
         else if (keyword(1:9) .eq. 'SLOPEMAX ') then
            read (string,*,err=10,end=10)  slpmax
         else if (keyword(1:7) .eq. 'ANGMAX ') then
            read (string,*,err=10,end=10)  angmax
         else if (keyword(1:7) .eq. 'INTMAX ') then
            read (string,*,err=10,end=10)  intmax
         end if
   10    continue
      end do
c
c     print header information about the optimization method
c
      if (iprint .gt. 0) then
         if (msav .eq. 0) then
            write (iout,20)
   20       format (/,' Steepest Descent Gradient Optimization :')
            write (iout,30)
   30       format (/,' SD Iter     F Value      G RMS      F Move',
     &                 '   X Move   Angle  FG Call  Comment',/)
         else
            write (iout,40)
   40       format (/,' Limited Memory BFGS Quasi-Newton',
     &                 ' Optimization :')
            write (iout,50)
   50       format (/,' QN Iter     F Value      G RMS      F Move',
     &                 '   X Move   Angle  FG Call  Comment',/)
         end if
         flush (iout)
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (x_old(nvar))
      allocate (g(nvar))
      allocate (g_old(nvar))
      allocate (p(nvar))
      allocate (q(nvar))
      allocate (r(nvar))
      allocate (h0(nvar))
      if (msav .ne. 0) then
         allocate (rho(msav))
         allocate (alpha(msav))
         allocate (s(nvar,msav))
         allocate (y(nvar,msav))
      end if
c
c     evaluate the function and get the initial gradient
c
      niter = nextiter - 1
      maxiter = niter + maxiter
      ncalls = ncalls + 1
      f = fgvalue (x0,g)
      f_old = f
      m = 0
      gamma = 1.0d0
      g_norm = 0.0d0
      g_rms = 0.0d0
      do i = 1, nvar
         g_norm = g_norm + g(i)*g(i)
         g_rms = g_rms + (g(i)*scale(i))**2
      end do
      g_norm = sqrt(g_norm)
      g_rms = sqrt(g_rms) / rms
      f_move = 0.5d0 * stpmax * g_norm
c
c     print initial information prior to first iteration
c
      if (iprint .gt. 0) then
         if (f.lt.1.0d8 .and. f.gt.-1.0d7 .and. g_rms.lt.1.0d5) then
            write (iout,60)  niter,f,g_rms,ncalls
   60       format (i6,f14.4,f11.4,29x,i7)
         else
            write (iout,70)  niter,f,g_rms,ncalls
   70       format (i6,d14.4,d11.4,29x,i7)
         end if
         flush (iout)
      end if
c
c     write initial intermediate prior to first iteration
c
      if (iwrite .gt. 0)  call optsave (niter,f,x0)
c
c     tests of the various termination criteria
c
      if (niter .ge. maxiter) then
         status = 'IterLimit'
         done = .true.
      end if
      if (f .le. fctmin) then
         status = 'SmallFct '
         done = .true.
      end if
      if (g_rms .le. grdmin) then
         status = 'SmallGrad'
         done = .true.
      end if
c
c     start of a new limited memory BFGS iteration
c
      do while (.not. done)
         niter = niter + 1
         muse = min(niter-1,msav)
         m = m + 1
         if (m .gt. msav)  m = 1
c
c     estimate Hessian diagonal and compute the Hg product
c
         do i = 1, nvar
            h0(i) = gamma
            q(i) = g(i)
         end do
         k = m
         do j = 1, muse
            k = k - 1
            if (k .eq. 0)  k = msav
            alpha(k) = 0.0d0
            do i = 1, nvar
               alpha(k) = alpha(k) + s(i,k)*q(i)
            end do
            alpha(k) = alpha(k) * rho(k)
            do i = 1, nvar
               q(i) = q(i) - alpha(k)*y(i,k)
            end do
         end do
         do i = 1, nvar
            r(i) = h0(i) * q(i)
         end do
         do j = 1, muse
            beta = 0.0d0
            do i = 1, nvar
               beta = beta + y(i,k)*r(i)
            end do
            beta = beta * rho(k)
            do i = 1, nvar
               r(i) = r(i) + s(i,k)*(alpha(k)-beta)
            end do
            k = k + 1
            if (k .gt. msav)  k = 1
         end do
c
c     set search direction and store current point and gradient
c
         do i = 1, nvar
            p(i) = -r(i)
            x_old(i) = x0(i)
            g_old(i) = g(i)
         end do
c
c     perform line search along the new conjugate direction
c
         status = blank
         call search (nvar,f,g,x0,p,f_move,angle,ncalls,fgvalue,status)
c
c     update variables based on results of this iteration
c
         if (msav .ne. 0) then
            ys = 0.0d0
            yy = 0.0d0
            do i = 1, nvar
               s(i,m) = x0(i) - x_old(i)
               y(i,m) = g(i) - g_old(i)
               ys = ys + y(i,m)*s(i,m)
               yy = yy + y(i,m)*y(i,m)
            end do
            gamma = abs(ys/yy)
            rho(m) = 1.0d0 / ys
         end if
c
c     get the sizes of the moves made during this iteration
c
         f_move = f_old - f
         f_old = f
         x_move = 0.0d0
         do i = 1, nvar
            x_move = x_move + ((x0(i)-x_old(i))/scale(i))**2
         end do
         x_move = sqrt(x_move) / rms
         if (coordtype .eq. 'INTERNAL') then
            x_move = radian * x_move
         end if
c
c     compute the rms gradient per optimization parameter
c
         g_rms = 0.0d0
         do i = 1, nvar
            g_rms = g_rms + (g(i)*scale(i))**2
         end do
         g_rms = sqrt(g_rms) / rms
c
c     test for error due to line search problems
c
         if (status.eq.'BadIntpln' .or. status.eq.'IntplnErr') then
            nerr = nerr + 1
            if (nerr .ge. maxerr)  done = .true.
         else
            nerr = 0
         end if
c
c     test for too many total iterations
c
         if (niter .ge. maxiter) then
            status = 'IterLimit'
            done = .true.
         end if
c
c     test the normal termination criteria
c
         if (f .le. fctmin) then
            status = 'SmallFct '
            done = .true.
         end if
         if (g_rms .le. grdmin) then
            status = 'SmallGrad'
            done = .true.
         end if
c
c     print intermediate results for the current iteration
c
         if (iprint .gt. 0) then
            if (done .or. mod(niter,iprint).eq.0) then
               if (f.lt.1.0d8 .and. f.gt.-1.0d7 .and.
     &             g_rms.lt.1.0d5 .and. f_move.lt.1.0d6 .and.
     &             f_move.gt.-1.0d5) then
                  write (iout,80)  niter,f,g_rms,f_move,x_move,
     &                             angle,ncalls,status
   80             format (i6,f14.4,f11.4,f12.4,f9.4,f8.2,i7,3x,a9)
               else
                  write (iout,90)  niter,f,g_rms,f_move,x_move,
     &                             angle,ncalls,status
   90             format (i6,d14.4,d11.4,d12.4,f9.4,f8.2,i7,3x,a9)
               end if
            end if
            flush (iout)
         end if
c
c     write intermediate results for the current iteration
c
         if (iwrite .gt. 0) then
            if (done .or. mod(niter,iwrite).eq.0) then
               call optsave (niter,f,x0)
            end if
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (x_old)
      deallocate (g)
      deallocate (g_old)
      deallocate (p)
      deallocate (q)
      deallocate (r)
      deallocate (h0)
      if (msav .ne. 0) then
         deallocate (rho)
         deallocate (alpha)
         deallocate (s)
         deallocate (y)
      end if
c
c     set final value of the objective function
c
      minimum = f
      if (iprint .gt. 0) then
         if (status.eq.'SmallGrad' .or. status.eq.'SmallFct ') then
            write (iout,100)  status
  100       format (/,' LBFGS  --  Normal Termination due to ',a9)
         else
            write (iout,110)  status
  110       format (/,' LBFGS  --  Incomplete Convergence due to ',a9)
         end if
         flush (iout)
      end if
      return
      end
