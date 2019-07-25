c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine tncg  --  truncated Newton optimization method  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "tncg" implements a truncated Newton optimization algorithm
c     in which a preconditioned linear conjugate gradient method is
c     used to approximately solve Newton's equations; special features
c     include use of an explicit sparse Hessian or finite-difference
c     gradient-Hessian products within the PCG iteration; the exact
c     Newton search directions can be used optionally; by default the
c     algorithm checks for negative curvature to prevent convergence
c     to a stationary point having negative eigenvalues; if a saddle
c     point is desired this test can be removed by disabling "negtest"
c
c     literature references:
c
c     J. W. Ponder and F. M Richards, "An Efficient Newton-like
c     Method for Molecular Mechanics Energy Minimization of
c     Large Molecules", Journal of Computational Chemistry,
c     8, 1016-1024 (1987)
c
c     R. S. Dembo and T. Steihaug, "Truncated-Newton Algorithms
c     for Large-Scale Unconstrained Optimization", Mathematical
c     Programming, 26, 190-212 (1983)
c
c     variables and parameters:
c
c     mode       determines optimization method; choice of
c                  Newton's method, truncated Newton, or
c                  truncated Newton with finite differencing
c     method     determines which type of preconditioning will
c                  be used on the Newton equations; choice
c                  of none, diagonal, 3x3 block diagonal,
c                  SSOR or incomplete Cholesky preconditioning
c     nvar       number of parameters in the objective function
c     minimum    upon return contains the best value of the
c                  function found during the optimization
c     f          contains current best value of function
c     x0         contains starting point upon input, upon
c                  return contains the best point found
c     g          contains gradient of current best point
c     h          contains the Hessian matrix values in an
c                  indexed linear array
c     h_mode     controls amount of Hessian matrix computed;
c                  either the full matrix, diagonal or none
c     h_init     points to the first Hessian matrix element
c                  associated with each parameter
c     h_stop     points to the last Hessian matrix element
c                  associated with each parameter
c     h_index    contains second parameter involved in each
c                  element of the Hessian array
c     h_diag     contains diagonal of the Hessian matrix
c     p          search direction resulting from pcg iteration
c     f_move     function decrease over last tncg iteration
c     f_old      function value at end of last iteration
c     x_move     rms movement per atom over last tn iteration
c     x_old      parameters value at end of last tn iteration
c     g_norm     Euclidian norm of the gradient vector
c     g_rms      root mean square gradient value
c     fg_call    cumulative number of function/gradient calls
c     grdmin     termination criterion based on RMS gradient
c     iprint     print iteration results every iprint iterations
c     iwrite     call user-supplied output every iwrite iterations
c     newhess    number of iterations between the computation
c                  of new Hessian matrix values
c     negtest    determines whether test for negative curvature
c                  is performed during the PCG iterations
c     maxiter    maximum number of tncg iterations to attempt
c
c     parameters used in the line search:
c
c     cappa      accuarcy of line search control  (0 < cappa < 1)
c     stpmin     minimum allowed line search step size
c     stpmax     maximum allowed line search step size
c     angmax     maximum angle between search and -grad directions
c     intmax     maximum number of interpolations in line search
c
c     required external routines:
c
c     fgvalue    function to evaluate function and gradient values
c     hmatrix    subroutine which evaluates Hessian diagonal
c                  and large off-diagonal matrix elements
c     optsave    subroutine to write out info about current status
c
c
      subroutine tncg (mode,method,nvar,x0,minimum,grdmin,
     &                       fgvalue,hmatrix,optsave)
      use sizes
      use atoms
      use hescut
      use inform
      use iounit
      use keys
      use linmin
      use math
      use minima
      use output
      use piorbs
      use potent
      implicit none
      integer i,fg_call
      integer nvar,nmax
      integer iter_tn,iter_cg
      integer next,newhess
      integer nerr,maxerr
      integer, allocatable :: h_init(:)
      integer, allocatable :: h_stop(:)
      integer, allocatable :: h_index(:)
      real*8 f,fgvalue,grdmin
      real*8 minimum,angle,rms
      real*8 x_move,f_move,f_old
      real*8 g_norm,g_rms
      real*8 x0(*)
      real*8, allocatable :: x_old(:)
      real*8, allocatable :: g(:)
      real*8, allocatable :: p(:)
      real*8, allocatable :: h_diag(:)
      real*8, allocatable :: h(:)
      logical done,negtest
      logical automode,automatic
      character*4 h_mode
      character*6 mode,method
      character*9 status
      character*9 info_solve
      character*9 info_search
      character*20 keyword
      character*240 record
      character*240 string
      save h_index,h
      external fgvalue
      external hmatrix
      external optsave
c
c
c     check number of variables and get type of optimization
c
      rms = sqrt(dble(nvar))
      if (coordtype .eq. 'CARTESIAN') then
         rms = rms / sqrt(3.0d0)
      else if (coordtype .eq. 'RIGIDBODY') then
         rms = rms / sqrt(6.0d0)
      end if
c
c     set default parameters for the optimization
c
      if (fctmin .eq. 0.0d0)  fctmin = -100000000.0d0
      if (iwrite .lt. 0)  iwrite = 1
      if (iprint .lt. 0)  iprint = 1
      if (maxiter .eq. 0)  maxiter = 1000
      if (nextiter .eq. 0)  nextiter = 1
      newhess = 1
      maxerr = 3
      done = .false.
      status = '         '
      negtest = .true.
      automode = .false.
      automatic = .false.
      if (mode .eq. 'AUTO')  automode = .true.
      if (method .eq. 'AUTO')  automatic = .true.
c
c     set default parameters for the line search
c
      if (stpmax .eq. 0.0d0)  stpmax = 5.0d0
      stpmin = 1.0d-16
      cappa = 0.1d0
      slpmax = 10000.0d0
      angmax = 180.0d0
      intmax = 8
c
c     search each line of the keyword file for options
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:7) .eq. 'FCTMIN ') then
            read (string,*,err=10,end=10)  fctmin
         else if (keyword(1:8) .eq. 'MAXITER ') then
            read (string,*,err=10,end=10)  maxiter
         else if (keyword(1:9) .eq. 'NEXTITER ') then
            read (string,*,err=10,end=10)  nextiter
         else if (keyword(1:8) .eq. 'NEWHESS ') then
            read (string,*,err=10,end=10)  newhess
         else if (keyword(1:12) .eq. 'SADDLEPOINT ') then
            negtest = .false.
         else if (keyword(1:8) .eq. 'STEPMIN ') then
            read (string,*,err=10,end=10)  stpmin
         else if (keyword(1:8) .eq. 'STEPMAX ') then
            read (string,*,err=10,end=10)  stpmax
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
c     initialize the function call and iteration counters
c
      fg_call = 0
      nerr = 0
      iter_tn = nextiter - 1
      maxiter = iter_tn + maxiter
c
c     print header information about the method used
c
      if (iprint .gt. 0) then
         if (mode .eq. 'NEWTON') then
            write (iout,20)
   20       format (/,' Full-Newton Conjugate-Gradient',
     &                 ' Optimization :')
         else if (mode .eq. 'TNCG') then
            write (iout,30)
   30       format (/,' Truncated-Newton Conjugate-Gradient',
     &                 ' Optimization :')
         else if (mode .eq. 'DTNCG') then
            write (iout,40)
   40       format (/,' Finite-Difference Truncated-Newton',
     &                 ' Conjugate-Gradient Optimization :')
         else if (mode .eq. 'AUTO') then
            write (iout,50)
   50       format (/,' Variable-Mode Truncated-Newton',
     &                 ' Conjugate-Gradient Optimization :')
         end if
         write (iout,60)  mode,method,grdmin
   60    format (/,' Algorithm : ',a6,5x,'Preconditioning : ',a6,5x,
     &              ' RMS Grad :',d9.2)
         write (iout,70)
   70    format (/,' TN Iter     F Value      G RMS      F Move',
     &              '   X Move   CG Iter  Solve   FG Call',/)
         flush (iout)
      end if
c
c     perform dynamic allocation of some local arrays
c
      nmax = 3 * n
      allocate (h_init(nmax))
      allocate (h_stop(nmax))
      allocate (x_old(nmax))
      allocate (g(nmax))
      allocate (p(nmax))
      allocate (h_diag(nmax))
      allocate (h_index((nmax*(nmax-1))/2))
      allocate (h((nmax*(nmax-1))/2))
c
c     evaluate the function and get the initial gradient
c
      iter_cg = 0
      fg_call = fg_call + 1
      f = fgvalue (x0,g)
      f_old = f
      g_norm = 0.0d0
      do i = 1, nvar
         x_old(i) = x0(i)
         g_norm = g_norm + g(i)**2
      end do
      g_norm = sqrt(g_norm)
      f_move = 0.5d0 * stpmax * g_norm
      g_rms = g_norm / rms
c
c     print initial information prior to first iteration
c
      if (iprint .gt. 0) then
         if (f.lt.1.0d8 .and. f.gt.-1.0d7 .and. g_rms.lt.1.0d5) then
            write (iout,80)  iter_tn,f,g_rms,fg_call
   80       format (i6,f14.4,f11.4,41x,i7)
         else
            write (iout,90)  iter_tn,f,g_rms,fg_call
   90       format (i6,d14.4,d11.4,41x,i7)
         end if
         flush (iout)
      end if
c
c     write initial intermediate prior to first iteration
c
      if (iwrite .gt. 0)  call optsave (iter_tn,f,x0)
c
c     check for termination criteria met by initial point
c
      if (g_rms .le. grdmin) then
         done = .true.
         minimum = f
         if (iprint .gt. 0) then
            write (iout,100)
  100       format (/,' TNCG  --  Normal Termination due to SmallGrad')
         end if
      else if (f .le. fctmin) then
         done = .true.
         minimum = f
         if (iprint .gt. 0) then
            write (iout,110)
  110       format (/,' TNCG  --  Normal Termination due to SmallFct')
         end if
      else if (iter_tn .ge. maxiter) then
         done = .true.
         minimum = f
         if (iprint .gt. 0) then
            write (iout,120)
  120       format (/,' TNCG  --  Incomplete Convergence',
     &                 ' due to IterLimit')
         end if
      end if
c
c     beginning of the outer truncated Newton iteration
c
      do while (.not. done)
         iter_tn = iter_tn + 1
c
c     if pisystem is present, update the molecular orbitals
c
         if (use_orbit) then
            reorbit = 1
            call picalc
            fg_call = fg_call + 1
            f = fgvalue (x0,g)
            reorbit = 0
         end if
c
c     choose the optimization mode based on the gradient value
c
         if (automode) then
            if (g_rms .ge. 3.0d0) then
               mode = 'TNCG'
            else
               mode = 'DTNCG'
            end if
         end if
c
c     decide on an optimal preconditioning based on the gradient
c
         if (automatic) then
            if (nvar .lt. 10) then
               method = 'DIAG'
               hesscut = 0.0d0
            else if (g_rms .ge. 10.0d0) then
               method = 'DIAG'
               hesscut = 1.0d0
            else if (g_rms .ge. 1.0d0) then
               method = 'ICCG'
               hesscut = 0.001d0 * nvar
               if (hesscut .gt. 1.0d0)  hesscut = 1.0d0
            else
               method = 'ICCG'
               hesscut = 0.001d0 * nvar
               if (hesscut .gt. 0.1d0)  hesscut = 0.1d0
            end if
         end if
c
c     compute needed portions of the Hessian matrix
c
         h_mode = 'FULL'
         if (mod(iter_tn-1,newhess) .ne. 0)  h_mode = 'NONE'
         if (mode.eq.'DTNCG' .and. method.eq.'NONE')  h_mode = 'NONE'
         if (mode.eq.'DTNCG' .and. method.eq.'DIAG')  h_mode = 'DIAG'
         call hmatrix (h_mode,x0,h,h_init,h_stop,h_index,h_diag)
c
c     find the next approximate Newton search direction
c
         call tnsolve (mode,method,negtest,nvar,p,x0,g,h,
     &                 h_init,h_stop,h_index,h_diag,iter_tn,
     &                 iter_cg,fg_call,fgvalue,info_solve)
c
c     perform a line search in the chosen direction
c
         info_search = '         '
         call search (nvar,f,g,x0,p,f_move,angle,fg_call,
     &                fgvalue,info_search)
         if (info_search .ne. ' Success ') then
            info_solve = info_search
         end if
c
c     update variables to reflect this iteration
c
         f_move = f_old - f
         f_old = f
         x_move = 0.0d0
         g_norm = 0.0d0
         do i = 1, nvar
            x_move = x_move + (x0(i)-x_old(i))**2
            x_old(i) = x0(i)
            g_norm = g_norm + g(i)**2
         end do
         x_move = sqrt(x_move)
         x_move = x_move / rms
         if (coordtype .eq. 'INTERNAL') then
            x_move = x_move * radian
         end if
         g_norm = sqrt(g_norm)
         g_rms = g_norm / rms
c
c     quit if the maximum number of iterations is exceeded
c
         if (iter_tn .ge. maxiter) then
            done = .true.
            status = 'IterLimit'
         end if
c
c     quit if the function value did not change
c
         if (f_move .eq. 0.0d0) then
            done = .true.
            status = 'NoMotion '
         end if
c
c     quit if either of the normal termination tests are met
c
         if (g_rms .le. grdmin) then
            done = .true.
            status = 'SmallGrad'
         else if (f .le. fctmin) then
            done = .true.
            status = 'SmallFct '
         end if
c
c     quit if the line search encounters successive problems
c
         if (info_search.eq.'BadIntpln' .or.
     &       info_search.eq.'IntplnErr') then
            nerr = nerr + 1
            if (nerr .ge. maxerr) then
               done = .true.
               status = info_search
            end if
         else
            nerr = 0
         end if
c
c     print intermediate results for the current iteration
c
         if (iprint .gt. 0) then
            if (done .or. mod(iter_tn,iprint).eq.0) then
               if (f.lt.1.0d8 .and. f.gt.-1.0d7 .and.
     &             g_rms.lt.1.0d5 .and. f_move.lt.1.0d6 .and.
     &             f_move.gt.-1.0d5) then
                  write (iout,130)  iter_tn,f,g_rms,f_move,x_move,
     &                              iter_cg,info_solve,fg_call
  130             format (i6,f14.4,f11.4,f12.4,f9.4,i8,3x,a9,i7)
               else
                  write (iout,140)  iter_tn,f,g_rms,f_move,x_move,
     &                              iter_cg,info_solve,fg_call
  140             format (i6,d14.4,d11.4,d12.4,f9.4,i8,3x,a9,i7)
               end if
               flush (iout)
            end if
         end if
c
c     write intermediate results for the current iteration
c
         if (iwrite .gt. 0) then
            if (done .or. mod(iter_tn,iwrite).eq.0) then
               call optsave (iter_tn,f,x0)
            end if
         end if
c
c     print the reason for terminating the optimization
c
         if (done) then
            minimum = f
            if (iprint .gt. 0) then
               if (g_rms.le.grdmin .or. f.le.fctmin) then
                  write (iout,150)  status
  150             format (/,' TNCG  --  Normal Termination due to ',a9)
               else
                  write (iout,160)  status
  160             format (/,' TNCG  --  Incomplete Convergence',
     &                       ' due to ',a9)
               end if
               flush (iout)
            end if
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (h_init)
      deallocate (h_stop)
      deallocate (x_old)
      deallocate (g)
      deallocate (p)
      deallocate (h_diag)
      deallocate (h_index)
      deallocate (h)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine tnsolve  --  approx linear equation solution   ##
c     ##                                                            ##
c     ################################################################
c
c
c     "tnsolve" uses a linear conjugate gradient method to find
c     an approximate solution to the set of linear equations
c     represented in matrix form by Hp = -g (Newton's equations)
c
c     status codes upon return:
c
c     TruncNewt    convergence to (truncated) Newton criterion
c     NegCurve     termination upon detecting negative curvature
c     OverLimit    maximum number of CG iterations exceeded
c
c
      subroutine tnsolve (mode,method,negtest,nvar,p,x0,g,h,
     &                    h_init,h_stop,h_index,h_diag,cycle,
     &                    iter_cg,fg_call,fgvalue,status)
      use sizes
      use output
      implicit none
      integer i,j,k,nvar,cycle
      integer iter,iter_cg
      integer fg_call,maxiter
      integer h_init(*)
      integer h_stop(*)
      integer h_index(*)
      real*8 alpha,beta,delta
      real*8 sigma,f_sigma
      real*8 fgvalue,eps
      real*8 g_norm,g_rms
      real*8 hj,gg,dq,rr,dd
      real*8 rs,rs_new,r_norm
      real*8 converge
      real*8 x0(*)
      real*8 g(*)
      real*8 p(*)
      real*8 h_diag(*)
      real*8 h(*)
      real*8, allocatable :: m(:)
      real*8, allocatable :: r(:)
      real*8, allocatable :: s(:)
      real*8, allocatable :: d(:)
      real*8, allocatable :: q(:)
      real*8, allocatable :: x_sigma(:)
      real*8, allocatable :: g_sigma(:)
      logical negtest
      character*6 mode,method
      character*9 status
      external fgvalue
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (m(nvar))
      allocate (r(nvar))
      allocate (s(nvar))
      allocate (d(nvar))
      allocate (q(nvar))
      allocate (x_sigma(nvar))
      allocate (g_sigma(nvar))
c
c     transformation using exact Hessian diagonal
c
      if (mode.ne.'DTNCG' .and. method.ne.'NONE') then
         do i = 1, nvar
            m(i) = 1.0d0 / sqrt(abs(h_diag(i)))
         end do
         do i = 1, nvar
            g(i) = g(i) * m(i)
            h_diag(i) = h_diag(i) * m(i) * m(i)
            do j = h_init(i), h_stop(i)
               k = h_index(j)
               h(j) = h(j) * m(i) * m(k)
            end do
         end do
      end if
c
c     setup prior to linear conjugate gradient iterations
c
      iter = 0
      gg = 0.0d0
      do i = 1, nvar
         p(i) = 0.0d0
         r(i) = -g(i)
         gg = gg + g(i)*g(i)
      end do
      g_norm = sqrt(gg)
      call precond (method,iter,nvar,s,r,h,h_init,
     &                  h_stop,h_index,h_diag)
      rs = 0.0d0
      do i = 1, nvar
         d(i) = s(i)
         rs = rs + r(i)*s(i)
      end do
      if (mode .eq. 'NEWTON') then
         eps = 1.0d-10
         maxiter = nvar
      else if (mode.eq.'TNCG' .or. mode.eq.'DTNCG') then
         delta = 1.0d0
         eps = delta / dble(cycle)
         g_rms = g_norm / sqrt(dble(nvar))
         eps = min(eps,g_rms)
         converge = 1.0d0
         eps = eps**converge
         maxiter = nint(10.0d0*sqrt(dble(nvar)))
      end if
      iter = 1
c
c     evaluate or estimate the matrix-vector product
c
      do while (.true.)
         if (mode.eq.'TNCG' .or. mode.eq.'NEWTON') then
            do i = 1, nvar
               q(i) = 0.0d0
            end do
            do i = 1, nvar
               q(i) = q(i) + h_diag(i)*d(i)
               do j = h_init(i), h_stop(i)
                  k = h_index(j)
                  hj = h(j)
                  q(i) = q(i) + hj*d(k)
                  q(k) = q(k) + hj*d(i)
               end do
            end do
         else if (mode .eq. 'DTNCG') then
            dd = 0.0d0
            do i = 1, nvar
               dd = dd + d(i)*d(i)
            end do
            sigma = 1.0d-7 / sqrt(dd)
            if (coordtype .eq. 'INTERNAL') then
               sigma = 1.0d-4 / sqrt(dd)
            end if
            do i = 1, nvar
               x_sigma(i) = x0(i) + sigma*d(i)
            end do
            fg_call = fg_call + 1
            f_sigma = fgvalue (x_sigma,g_sigma)
            do i = 1, nvar
               q(i) = (g_sigma(i)-g(i)) / sigma
            end do
         end if
c
c     check for a direction of negative curvature
c
         dq = 0.0d0
         do i = 1, nvar
            dq = dq + d(i)*q(i)
         end do
         if (negtest) then
            if (dq .le. 0.0d0) then
               if (iter .eq. 1) then
                  do i = 1, nvar
                     p(i) = d(i)
                  end do
               end if
               status = ' NegCurve'
               goto 10
            end if
         end if
c
c     test the truncated Newton termination criterion
c
         alpha = rs / dq
         rr = 0.0d0
         do i = 1, nvar
            p(i) = p(i) + alpha*d(i)
            r(i) = r(i) - alpha*q(i)
            rr = rr + r(i)*r(i)
         end do
         r_norm = sqrt(rr)
         if (r_norm/g_norm .le. eps) then
            status = 'TruncNewt'
            goto 10
         end if
c
c     solve the preconditioning equations
c
         call precond (method,iter,nvar,s,r,h,h_init,
     &                     h_stop,h_index,h_diag)
c
c     update the truncated Newton direction
c
         rs_new = 0.0d0
         do i = 1, nvar
            rs_new = rs_new + r(i)*s(i)
         end do
         beta = rs_new / rs
         rs = rs_new
         do i = 1, nvar
            d(i) = s(i) + beta*d(i)
         end do
c
c     check for overlimit, then begin next iteration
c
         if (iter .ge. maxiter) then
            status = 'OverLimit'
            goto 10
         end if
         iter = iter + 1
      end do
c
c     retransform and increment total iterations, then terminate
c
   10 continue
      if (mode.ne.'DTNCG' .and. method.ne.'NONE') then
         do i = 1, nvar
            p(i) = p(i) * m(i)
            g(i) = g(i) / m(i)
         end do
      end if
      iter_cg = iter_cg + iter
c
c     perform deallocation of some local arrays
c
      deallocate (m)
      deallocate (r)
      deallocate (s)
      deallocate (d)
      deallocate (q)
      deallocate (x_sigma)
      deallocate (g_sigma)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine precond  --  precondition linear CG method  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "precond" solves a simplified version of the Newton equations
c     Ms = r, and uses the result to precondition linear conjugate
c     gradient iterations on the full Newton equations in "tnsolve"
c
c     reference for incomplete Cholesky factorization :
c
c     T. A. Manteuffel, "An Incomplete Factorization Technique
c     for Positive Definite Linear Systems", Mathematics of
c     Computation, 34, 473-497 (1980); the present method is
c     based upon the SICCG(0) method described in this paper
c
c     types of preconditioning methods :
c
c     none     use no preconditioning at all
c     diag     exact Hessian diagonal preconditioning
c     block    3x3 block diagonal preconditioning
c     ssor     symmetric successive over-relaxation
c     iccg     shifted incomplete Cholesky factorization
c
c
      subroutine precond (method,iter,nvar,s,r,h,h_init,
     &                       h_stop,h_index,h_diag)
      use sizes
      use inform
      use iounit
      implicit none
      integer i,j,k,ii,kk
      integer iii,kkk,iter
      integer nvar,nblock
      integer ix,iy,iz,icount
      integer h_init(*)
      integer h_stop(*)
      integer h_index(*)
      integer, allocatable :: c_init(:)
      integer, allocatable :: c_stop(:)
      integer, allocatable :: c_index(:)
      integer, allocatable :: c_value(:)
      real*8 f_i,f_k
      real*8 omega,factor
      real*8 maxalpha,alpha
      real*8 a(6),b(3)
      real*8 h_diag(*)
      real*8 h(*)
      real*8 s(*)
      real*8 r(*)
      real*8, allocatable :: diag(:)
      real*8, allocatable :: f_diag(:)
      real*8, allocatable :: f(:)
      logical stable
      character*6 method
      save f,f_diag,stable
c
c
c     perform dynamic allocation of some local arrays
c
      if (method .eq. 'SSOR')  allocate (diag(nvar))
      if (method .eq. 'ICCG') then
         if (iter .eq. 0) then
            allocate (c_init(nvar))
            allocate (c_stop(nvar))
            allocate (c_index((nvar*(nvar-1))/2))
            allocate (c_value((nvar*(nvar-1))/2))
         end if
         if (.not. allocated(f_diag))  allocate (f_diag(nvar))
         if (.not. allocated(f))  allocate (f((nvar*(nvar-1))/2))
      end if
c
c     use no preconditioning, using M = identity matrix
c
      if (method .eq. 'NONE') then
         do i = 1, nvar
            s(i) = r(i)
         end do
      end if
c
c     diagonal preconditioning, using M = abs(Hessian diagonal)
c
      if (method .eq. 'DIAG') then
         do i = 1, nvar
            s(i) = r(i) / abs(h_diag(i))
         end do
      end if
c
c     block diagonal preconditioning with exact atom blocks
c     (using M = 3x3 blocks from diagonal of full Hessian)
c
      if (method .eq. 'BLOCK') then
         nblock = 3
         do i = 1, nvar/3
            iz = 3 * i
            iy = iz - 1
            ix = iz - 2
            a(1) = h_diag(ix)
            if (h_index(h_init(ix)) .eq. iy) then
               a(2) = h(h_init(ix))
            else
               a(2) = 0.0d0
            end if
            if (h_index(h_init(ix)+1) .eq. iz) then
               a(3) = h(h_init(ix)+1)
            else
               a(3) = 0.0d0
            end if
            a(4) = h_diag(iy)
            if (h_index(h_init(iy)) .eq. iz) then
               a(5) = h(h_init(iy))
            else
               a(5) = 0.0d0
            end if
            a(6) = h_diag(iz)
            b(1) = r(ix)
            b(2) = r(iy)
            b(3) = r(iz)
            call cholesky (nblock,a,b)
            s(ix) = b(1)
            s(iy) = b(2)
            s(iz) = b(3)
         end do
      end if
c
c     symmetric successive over-relaxation (SSOR) preconditioning
c     (using M = (D/w+U)T * (D/w)-1 * (D/w+U) with 0 < w < 2)
c
      if (method .eq. 'SSOR') then
         omega = 1.0d0
         factor = 2.0d0 - omega
         do i = 1, nvar
            s(i) = r(i) * factor
            diag(i) = h_diag(i) / omega
         end do
         do i = 1, nvar
            s(i) = s(i) / diag(i)
            do j = h_init(i), h_stop(i)
               k = h_index(j)
               s(k) = s(k) - h(j)*s(i)
            end do
         end do
         do i = nvar, 1, -1
            s(i) = s(i) * diag(i)
            do j = h_init(i), h_stop(i)
               k = h_index(j)
               s(i) = s(i) - h(j)*s(k)
            end do
            s(i) = s(i) / diag(i)
         end do
      end if
c
c     factorization phase of incomplete cholesky preconditioning
c
      if (method.eq.'ICCG' .and. iter.eq.0) then
         call column (nvar,h_init,h_stop,h_index,
     &                c_init,c_stop,c_index,c_value)
         stable = .true.
         icount = 0
         maxalpha = 2.1d0
         alpha = -0.001d0
   10    continue
         if (alpha .le. 0.0d0) then
            alpha = alpha + 0.001d0
         else
            alpha = 2.0d0 * alpha
         end if
         if (alpha .gt. maxalpha) then
            stable = .false.
            if (verbose) then
               write (iout,20)
   20          format (' PRECOND  --  Incomplete Cholesky is',
     &                    ' Unstable, using Diagonal Method')
            end if
         else
            factor = 1.0d0 + alpha
            do i = 1, nvar
               f_diag(i) = factor * h_diag(i)
               do j = c_init(i), c_stop(i)
                  k = c_index(j)
                  f_i = f(c_value(j))
                  f_diag(i) = f_diag(i) - f_i*f_i*f_diag(k)
                  icount = icount + 1
               end do
               if (f_diag(i) .le. 0.0d0)  goto 10
               if (f_diag(i) .lt. 1.0d-7)  f_diag(i) = 1.0d-7
               f_diag(i) = 1.0d0 / f_diag(i)
               do j = h_init(i), h_stop(i)
                  k = h_index(j)
                  f(j) = h(j)
                  ii = c_init(i)
                  kk = c_init(k)
                  do while (ii.le.c_stop(i) .and. kk.le.c_stop(k))
                     iii = c_index(ii)
                     kkk = c_index(kk)
                     if (iii .lt. kkk) then
                        ii = ii + 1
                     else if (kkk .lt. iii) then
                        kk = kk + 1
                     else
                        f_i = f(c_value(ii))
                        f_k = f(c_value(kk))
                        f(j) = f(j) - f_i*f_k*f_diag(iii)
                        ii = ii + 1
                        kk = kk + 1
                        icount = icount + 1
                     end if
                  end do
               end do
            end do
            if (verbose) then
               write (iout,30)  icount,alpha
   30          format (' PRECOND  --  Incomplete Cholesky',i12,
     &                    ' Operations',f8.3,' Alpha Value')
            end if
         end if
      end if
c
c     solution phase of incomplete cholesky preconditioning
c
      if (method .eq. 'ICCG') then
         if (stable) then
            do i = 1, nvar
               s(i) = r(i)
            end do
            do i = 1, nvar
               s(i) = s(i) * f_diag(i)
               do j = h_init(i), h_stop(i)
                  k = h_index(j)
                  s(k) = s(k) - f(j)*s(i)
               end do
            end do
            do i = nvar, 1, -1
               s(i) = s(i) / f_diag(i)
               do j = h_init(i), h_stop(i)
                  k = h_index(j)
                  s(i) = s(i) - f(j)*s(k)
               end do
               s(i) = s(i) * f_diag(i)
            end do
         else
            do i = 1, nvar
               s(i) = r(i) / abs(h_diag(i))
            end do
         end if
      end if
c
c     perform deallocation of some local arrays
c
      if (method .eq. 'SSOR')  deallocate (diag)
      if (method .eq. 'ICCG') then
         if (iter .eq. 0) then
            deallocate (c_init)
            deallocate (c_stop)
            deallocate (c_index)
            deallocate (c_value)
         end if
      end if
      return
      end
