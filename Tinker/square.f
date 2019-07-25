c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine square  --  nonlinear least squares with bounds  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "square" is a nonlinear least squares routine derived from
c     the IMSL routine BCLSF and More's Minpack routine LMDER; the
c     Jacobian is estimated by finite differences and bounds can
c     be specified for the variables to be refined
c
c     arguments and variables :
c
c     n         number of least squares variables
c     m         number of residual functions
c     xlo       vector with the lower bounds for the variables
c     xhi       vector with the upper bounds for the variables
c     xscale    vector with the diagonal scaling matrix for variables
c     xc        vector with variable values at the approximate solution
c     fc        vector with the residuals at the approximate solution
c     fp        vector containing the updated residuals
c     xp        vector containing the updated point
c     sc        vector containing the last step taken
c     gc        vector with gradient estimate at approximate solution
c     fjac      matrix with estimate of Jacobian at approximate solution
c     iactive   vector showing if variable is at upper or lower bound
c     ipvt      vector with permutation matrix used in QR factorization
c                 of the Jacobian at the approximate solution
c     stpmax    scalar containing maximum allowed step size
c     delta     scalar containing the trust region radius
c
c     required external routines :
c
c     rsdvalue   subroutine to evaluate residual function values
c     lsqwrite   subroutine to write out info about current status
c
c
      subroutine square (n,m,xlo,xhi,xc,fc,gc,fjac,grdmin,
     &                          rsdvalue,lsqwrite)
      use sizes
      use inform
      use iounit
      use keys
      use minima
      implicit none
      integer i,j,k,m,n
      integer icode,next
      integer niter,ncalls
      integer nactive
      integer nbigstp,ndigit
      integer, allocatable :: iactive(:)
      integer, allocatable :: ipvt(:)
      real*8 amu,delta,epsfcn
      real*8 fcnorm,fpnorm
      real*8 gcnorm,ganorm
      real*8 precise,eps
      real*8 grdmin,stpnorm
      real*8 stpmax,stpmin
      real*8 rftol,faketol
      real*8 xtemp,stepsz
      real*8 sum,temp
      real*8 xc(*)
      real*8 xlo(*)
      real*8 xhi(*)
      real*8 fc(*)
      real*8 gc(*)
      real*8, allocatable :: xp(:)
      real*8, allocatable :: ga(:)
      real*8, allocatable :: gs(:)
      real*8, allocatable :: sc(:)
      real*8, allocatable :: sa(:)
      real*8, allocatable :: xsa(:)
      real*8, allocatable :: xscale(:)
      real*8, allocatable :: rdiag(:)
      real*8, allocatable :: fp(:)
      real*8, allocatable :: ftemp(:)
      real*8, allocatable :: qtf(:)
      real*8 fjac(m,*)
      logical done,first
      logical gauss,bigstp
      logical pivot
      character*20 keyword
      character*240 record
      character*240 string
      external rsdvalue
      external lsqwrite
      external precise
c
c
c     initialize various counters and status code
c
      niter = 0
      ncalls = 0
      nbigstp = 0
      done = .false.
c
c     setup the default tolerances and parameter values
c
      ndigit = 10
      eps = precise (2)
      eps = max(eps,10.0d0**(-ndigit))
      if (maxiter .eq. 0)  maxiter = 100
      if (iprint .lt. 0)  iprint = 1
      if (iwrite .lt. 0)  iwrite = 1
      if (fctmin .eq. 0.0d0)  fctmin = eps
      if (grdmin .eq. 0.0d0)  grdmin = eps**(1.0d0/3.0d0)
      epsfcn = sqrt(eps)
      delta = 0.0d0
      stpmax = 1000.0d0 * sqrt(dble(n))
      stpmin = eps**(2.0d0/3.0d0)
      rftol = eps**(2.0d0/3.0d0)
      faketol = 100.0d0 * eps
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
         else if (keyword(1:9) .eq. 'PRINTOUT ') then
            read (string,*,err=10,end=10)  iprint
         else if (keyword(1:9) .eq. 'WRITEOUT ') then
            read (string,*,err=10,end=10)  iwrite
         end if
   10    continue
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (iactive(n))
      allocate (ipvt(n))
      allocate (xp(n))
      allocate (ga(n))
      allocate (gs(n))
      allocate (sc(n))
      allocate (sa(n))
      allocate (xsa(n))
      allocate (xscale(n))
      allocate (rdiag(n))
      allocate (fp(m))
      allocate (ftemp(m))
      allocate (qtf(m))
c
c     check feasibility of variables and use bounds if needed
c
      nactive = 0
      do j = 1, n
         if (xc(j) .lt. xlo(j)) then
            xc(j) = xlo(j)
            iactive(j) = -1
         else if (xc(j) .gt. xhi(j)) then
            xc(j) = xhi(j)
            iactive(j) = 1
         else
            nactive = nactive + 1
            iactive(j) = 0
         end if
      end do
c
c     evaluate the function at the initial point
c
      ncalls = ncalls + 1
      call rsdvalue (n,m,xc,fc)
      fcnorm = 0.0d0
      do i = 1, m
         fcnorm = fcnorm + fc(i)**2
      end do
      fcnorm = 0.5d0 * fcnorm
c
c     evaluate the Jacobian at the initial point by finite
c     differences; replace loop with user routine if desired
c
      do j = 1, n
         stepsz = epsfcn * abs(xc(j))
         if (stepsz .lt. epsfcn)  stepsz = epsfcn
         if (xc(j) .lt. 0.0d0)  stepsz = -stepsz
         xtemp = xc(j)
         xc(j) = xtemp + stepsz
         ncalls = ncalls + 1
         call rsdvalue (n,m,xc,ftemp)
         xc(j) = xtemp
         do i = 1, m
            fjac(i,j) = (ftemp(i)-fc(i)) / stepsz
         end do
      end do
c
c     compute More's adaptive variable scale factors
c
      do j = 1, n
         temp = 0.0d0
         do i = 1, m
            temp = temp + fjac(i,j)**2
         end do
         xscale(j) = sqrt(temp)
         if (xscale(j) .eq. 0.0d0)  xscale(j) = 1.0d0
      end do
c
c     compute the total gradient vector for all variables
c
      do j = 1, n
         gc(j) = 0.0d0
         do i = 1, m
            gc(j) = gc(j) + fjac(i,j)*fc(i)
         end do
      end do
c
c     compute the norm of the scaled total gradient
c     and the scaled gradient for active variables
c
      gcnorm = 0.0d0
      ganorm = 0.0d0
      do j = 1, n
         gs(j) = gc(j) * max(abs(xc(j)),1.0d0/xscale(j))
         gcnorm = gcnorm + gs(j)**2
         if (iactive(j) .eq. 0) then
            ganorm = ganorm + gs(j)**2
         end if
      end do
      gcnorm = sqrt(gcnorm/dble(n))
      if (nactive .ne. 0)  ganorm = sqrt(ganorm/dble(nactive))
c
c     print out information about initial conditions
c
      if (iprint .gt. 0) then
         write (iout,20)
   20    format (/,' Levenberg-Marquardt Nonlinear Least Squares :')
         write (iout,30)
   30    format (/,' LS Iter     F Value      Total G     Active G',
     &              '    N Active   F Calls',/)
         if (max(fcnorm,gcnorm) .lt. 10000000.0d0) then
            write (iout,40)  niter,fcnorm,gcnorm,ganorm,nactive,ncalls
   40       format (i6,f14.4,2f13.4,2i10)
         else
            write (iout,50)  niter,fcnorm,gcnorm,ganorm,nactive,ncalls
   50       format (i6,f14.4,2d13.4,2i10)
         end if
      end if
c
c     write out the parameters, derivatives and residuals
c
      if (iwrite .ne. 0)  call lsqwrite (niter,m,xc,gs,fc)
c
c     check stopping criteria at the initial point; test the
c     absolute function value and gradient norm for termination
c
      if (fcnorm .le. fctmin)  return
      if (ganorm .le. grdmin)  return
c
c     start of the main body of least squares iteration
c
   60 continue
      niter = niter + 1
c
c     repack the Jacobian to include only active variables
c
      if (nactive .ne. n) then
         k = 0
         do j = 1, n
            if (iactive(j) .ne. 0) then
               if (k .eq. 0)  k = j
            else
               if (k .ne. 0) then
                  do i = 1, m
                     fjac(i,k) = fjac(i,j)
                  end do
                  k = k + 1
               end if
            end if
         end do
      end if
c
c     repack scale factors and gradient for active variables
c
      k = 0
      do j = 1, n
         if (iactive(j) .eq. 0) then
            k = k + 1
            xsa(k) = xscale(j)
            ga(k) = gc(j)
         end if
      end do
c
c     compute the QR factorization of the Jacobian
c
      pivot = .true.
      call qrfact (nactive,m,fjac,pivot,ipvt,rdiag)
c
c     compute the vector Q(transpose) * residuals
c
      do i = 1, m
         qtf(i) = fc(i)
      end do
      do j = 1, nactive
         if (fjac(j,j) .ne. 0.0d0) then
            sum = 0.0d0
            do i = j, m
               sum = sum + fjac(i,j)*qtf(i)
            end do
            temp = -sum / fjac(j,j)
            do i = j, m
               qtf(i) = qtf(i) + fjac(i,j)*temp
            end do
         end if
         fjac(j,j) = rdiag(j)
      end do
c
c     compute the Levenberg-Marquardt step
c
      icode = 6
      first = .true.
      do while (icode .ge. 4)
         call lmstep (nactive,m,ga,fjac,ipvt,xsa,qtf,stpmax,
     &                      delta,amu,first,sa,gauss)
c
c     unpack the step vector to include all variables
c
         k = 0
         do i = 1, n
            if (iactive(i) .ne. 0) then
               sc(i) = 0.0d0
            else
               k = k + 1
               sc(i) = sa(k)
            end if
         end do
c
c     check new point and update the trust region
c
         call trust (n,m,xc,fcnorm,gc,fjac,ipvt,sc,sa,xscale,gauss,
     &               stpmax,delta,icode,xp,fc,fp,fpnorm,bigstp,ncalls,
     &               xlo,xhi,nactive,stpmin,rftol,faketol,rsdvalue)
      end do
      if (icode .eq. 1)  done = .true.
c
c     update to the new variables and residuals
c
      do j = 1, n
         xc(j) = xp(j)
      end do
      do i = 1, m
         fc(i) = fp(i)
      end do
      fcnorm = fpnorm
c
c     update the active vs inactive status of the variables;
c     in a true active set strategy, at most one constraint
c     is added to the active set per iteration
c
      do j = 1, n
         if (iactive(j) .eq. 0) then
            if (abs(xc(j)-xlo(j)) .le. eps) then
               nactive = nactive - 1
               iactive(j) = -1
c              goto 99
            else if (abs(xc(j)-xhi(j)) .le. eps) then
               nactive = nactive - 1
               iactive(j) = 1
c              goto 99
            end if
         end if
      end do
c  99 continue
c
c     evaluate the Jacobian at the new point using finite
c     differences; replace loop with user routine if desired
c
      do j = 1, n
         stepsz = epsfcn * max(abs(xc(j)),1.0d0/xscale(j))
         if (xc(j) .lt. 0.0d0)  stepsz = -stepsz
         xtemp = xc(j)
         xc(j) = xtemp + stepsz
         ncalls = ncalls + 1
         call rsdvalue (n,m,xc,ftemp)
         xc(j) = xtemp
         do i = 1, m
            fjac(i,j) = (ftemp(i)-fc(i)) / stepsz
         end do
      end do
c
c     compute More's adaptive variable scale factors
c
      do j = 1, n
         temp = 0.0d0
         do i = 1, m
            temp = temp + fjac(i,j)**2
         end do
         xscale(j) = max(xscale(j),sqrt(temp))
      end do
c
c     compute the total gradient vector for all variables
c
      do j = 1, n
         gc(j) = 0.0d0
         do i = 1, m
            gc(j) = gc(j) + fjac(i,j)*fc(i)
         end do
      end do
c
c     compute the norm of the scaled total gradient
c     and the scaled gradient for active variables
c
      gcnorm = 0.0d0
      ganorm = 0.0d0
      do j = 1, n
         gs(j) = gc(j) * max(abs(xc(j)),1.0d0/xscale(j))
         gcnorm = gcnorm + gs(j)**2
         if (iactive(j) .eq. 0) then
            ganorm = ganorm + gs(j)**2
         end if
      end do
      gcnorm = sqrt(gcnorm/dble(n))
      if (nactive .ne. 0)  ganorm = sqrt(ganorm/dble(nactive))
c
c     print out information about current iteration
c
      if (iprint.ne.0 .and. mod(niter,iprint).eq.0) then
         if (max(fcnorm,gcnorm) .lt. 10000000.0d0) then
            write (iout,70)  niter,fcnorm,gcnorm,ganorm,nactive,ncalls
   70       format (i6,f14.4,2f13.4,2i10)
         else
            write (iout,80)  niter,fcnorm,gcnorm,ganorm,nactive,ncalls
   80       format (i6,f14.4,2d13.4,2i10)
         end if
      end if
c
c     check stopping criteria at the new point; test the absolute
c     function value, the gradient norm and step for termination
c
      if (fcnorm .le. fctmin)  done = .true.
      if (ganorm .le. grdmin)  done = .true.
      stpnorm = 0.0d0
      do j = 1, n
         temp = max(abs(xc(j)),1.0d0/xscale(j))
         stpnorm = stpnorm + (sc(j)/temp)**2
      end do
      stpnorm = sqrt(stpnorm/n)
      if (stpnorm .le. stpmin)  done = .true.
c
c     check for inactive variables that can be made active;
c     in a true active set strategy, variables are released
c     one at a time at a minimum of the current active set
c
c     if (done) then
      if (nactive .ne. n) then
         do j = 1, n
            if (iactive(j).eq.-1 .and. gc(j).lt.0.0d0) then
               nactive = nactive + 1
               iactive(j) = 0
               done = .false.
c              goto 99
            else if (iactive(j).eq.1 .and. gc(j).gt.0.0d0) then
               nactive = nactive + 1
               iactive(j) = 0
               done = .false.
c              goto 99
            end if
         end do
c  99    continue
      end if
c     end if
c
c     if still done, then normal termination has been achieved
c
      if (done) then
         write (iout,90)
   90    format (/,' SQUARE  --  Normal Termination of Least Squares')
c
c     check the limit on the number of iterations
c
      else if (niter .ge. maxiter) then
         done = .true.
         write (iout,100)
  100    format (/,' SQUARE  --  Maximum Number of Allowed Iterations')
c
c     check for termination due to relative function convergence
c
      else if (icode .eq. 2) then
         done = .true.
         write (iout,110)
  110    format (/,' SQUARE  --  Relative Function Convergence',
     &           //,' Both the scaled actual and predicted',
     &              ' reductions in the function',
     &           /,' are less than or equal to the relative',
     &              ' convergence tolerance')
c
c     check for termination due to false convergence
c
      else if (icode .eq. 3) then
         done = .true.
         write (iout,120)
  120    format (/,' SQUARE  --  Possible False Convergence',
     &           //,' The iterates appear to be converging to',
     &              ' a noncritical point due',
     &           /,' to bad gradient information, discontinuous',
     &              ' function, or stopping',
     &           /,' tolerances being too tight')
c
c     check for several consecutive maximum steps taken
c
      else if (bigstp) then
         nbigstp = nbigstp + 1
         if (nbigstp .eq. 5) then
            done = .true.
            write (iout,130)
  130       format (/,' SQUARE  --  Five Consecutive Maximum',
     &                 ' Length Steps',
     &              //,' Either the function is unbounded below,',
     &                 ' or has a finite',
     &              /,' asymptote in some direction, or STEPMAX',
     &                 ' is too small')
         end if
c
c     no reason to quit, so prepare to take another step
c
      else
         nbigstp = 0
      end if
c
c     write out the parameters, derivatives and residuals
c
      if (iwrite.ne.0 .and. mod(niter,iwrite).eq.0) then
         if (.not. done)  call lsqwrite (niter,m,xc,gs,fc)
      end if
c
c     continue with the next iteration if not finished
c
      if (.not. done)  goto 60
c
c     perform deallocation of some local arrays
c
      deallocate (iactive)
      deallocate (ipvt)
      deallocate (xp)
      deallocate (ga)
      deallocate (gs)
      deallocate (sc)
      deallocate (sa)
      deallocate (xsa)
      deallocate (xscale)
      deallocate (rdiag)
      deallocate (fp)
      deallocate (ftemp)
      deallocate (qtf)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine lmstep  --  computes Levenberg-Marquardt step  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "lmstep" computes the Levenberg-Marquardt step during a
c     nonlinear least squares calculation; based on ideas from
c     the Minpack routine LMPAR together with the internal doubling
c     strategy of Dennis and Schnabel
c
c     arguments and variables :
c
c     n        number of least squares variables
c     m        number of residual functions
c     ga       vector with the gradient of the residual vector
c     a        array of size n by n which on input contains in the full
c                upper triangle of the matrix r resulting from the QR
c                factorization of the Jacobian; on output the full upper
c                triangle is unaltered, and the strict lower triangle
c                contains the strict lower triangle of the matrix l
c                which is the Cholesky factor of (j**t)*j + amu*xscale
c     ipvt     vector with pivoting information from QR factorization
c     xscale   vector with the diagonal scaling matrix for variables
c     qtf      vector with first n elements of Q(transpose)
c                * (scaled residual)
c     amu      scalar with initial estimate of the Levenberg-Marquardt
c                parameter on input, and the final estimate of the
c                parameter on output
c     first    logical flag set true only if this is the first
c                call to this routine in this iteration
c     sa       vector with the Levenberg-Marquardt step
c     gnstep   vector with the Gauss-Newton step
c     gauss    logical flag set true if the Gauss-Newton step
c                is acceptable, and false otherwise
c     diag     vector with the diagonal elements of the Cholesky
c                factor of (j**t)*j + amu*xscale
c
c
      subroutine lmstep (n,m,ga,a,ipvt,xscale,qtf,stpmax,
     &                      delta,amu,first,sa,gauss)
      implicit none
      integer i,j,k
      integer m,n,nsing
      integer ipvt(*)
      real*8 stpmax,delta,amu
      real*8 alow,alpha
      real*8 amulow,amuhi
      real*8 beta,small,sum
      real*8 deltap,gnleng
      real*8 phi,phip,phipi
      real*8 sgnorm,precise
      real*8 high,tiny
      real*8 stplen,temp
      real*8 ga(*)
      real*8 xscale(*)
      real*8 qtf(*)
      real*8 sa(*)
      real*8, allocatable :: gnstep(:)
      real*8, allocatable :: diag(:)
      real*8, allocatable :: work1(:)
      real*8, allocatable :: work2(:)
      real*8 a(m,*)
      logical first,gauss
      logical done
      save deltap,nsing
      save phi,phip
      external precise
c
c
c     set smallest floating point magnitude and spacing
c
      tiny = precise (1)
      small = precise (2)
c
c     perform dynamic allocation of some local arrays
c
      allocate (gnstep(n))
      allocate (diag(n))
      allocate (work1(n))
      allocate (work2(n))
c
c     if initial trust region is not provided by the user,
c     compute and use the length of the Cauchy step given
c     by beta = norm2(r*trans(p)*d**(-2)*g)**2
c
      if (delta .eq. 0.0d0) then
         amu = 0.0d0
         do i = 1, n
            work1(i) = ga(i) / xscale(i)
         end do
         alpha = 0.0d0
         do i = 1, n
            alpha = alpha + work1(i)**2
         end do
         beta = 0.0d0
         do i = 1, n
            temp = 0.0d0
            do j = i, n
               k = ipvt(j)
               temp = temp + a(i,j)*ga(k)/xscale(k)**2
            end do
            beta = beta + temp**2
         end do
         if (beta .le. tiny) then
            delta = alpha * sqrt(alpha)
         else
            delta = alpha * sqrt(alpha)/beta
         end if
         delta = min(delta,stpmax)
      end if
c
c     the following is done only on the first time through
c     this iteration: (1) compute the Gauss-Newton step;
c     if the Jacobian is rank-deficient, obtain a least
c     squares solution, (2) compute the length of the scaled
c     Gauss-Newton step, (3) compute the norm of the scaled
c     gradient used in computing an upper bound for "amu"
c
      if (first) then
         nsing = n
         do j = 1, n
            if (a(j,j).eq.0.0d0 .and. nsing.eq.n)  nsing = j - 1
            if (nsing .lt. n)  work1(j) = 0.0d0
         end do
         work1(nsing) = qtf(nsing) / a(nsing,nsing)
         do j = nsing-1, 1, -1
            sum = 0.0d0
            do i = j+1, nsing
               sum = sum + a(j,i)*work1(i)
            end do
            work1(j) = (qtf(j)-sum) / a(j,j)
         end do
         do j = 1, n
            gnstep(ipvt(j)) = -work1(j)
         end do
c
c     find the length of scaled Gauss-Newton step
c
         do j = 1, n
            work1(j) = xscale(j) * gnstep(j)
         end do
         gnleng = 0.0d0
         do j = 1, n
            gnleng = gnleng + work1(j)**2
         end do
         gnleng = sqrt(gnleng)
c
c     find the length of the scaled gradient
c
         do j = 1, n
            work1(j) = ga(j) / xscale(j)
         end do
         sgnorm = 0.0d0
         do j = 1, n
            sgnorm = sgnorm + work1(j)**2
         end do
         sgnorm = sqrt(sgnorm)
      end if
c
c     set the bounds on the computed step
c
      high = 1.5d0
      alow = 0.75d0
c
c     check to see if the Gauss-Newton step is acceptable
c
      if (gnleng .le. high*delta) then
         gauss = .true.
         do j = 1, n
            sa(j) = gnstep(j)
         end do
         amu = 0.0d0
         delta = min(delta,gnleng)
c
c     the Gauss-Newton step is rejected, find a nontrivial step;
c     first compute a starting value of "amu" if previous step
c     was not a Gauss-Newton step
c
      else
         gauss = .false.
         if (amu .gt. 0.0d0)
     &      amu = amu - ((phi+deltap)/delta)*(((deltap-delta)+phi)/phip)
         phi = gnleng - delta
c
c     if the Jacobian is not rank deficient, the Newton step
c     provides a lower bound for "amu"; else set bound to zero
c
         if (nsing .eq. n) then
            if (first) then
               first = .false.
               do j = 1, n
                  k = ipvt(j)
                  work1(j) = gnstep(k) * xscale(k)**2
               end do
c
c     obtain trans(r**-1)*(trans(p)*s) by solving the
c     system of equations trans(r)*work1 = work1
c
               work1(n) = work1(n) / a(n,n)
               do j = n-1, 1, -1
                  sum = 0.0d0
                  do i = j+1, n
                     sum = sum + a(j,i)*work1(i)
                  end do
                  work1(j) = (work1(j)-sum) / a(j,j)
               end do
               phipi = 0.0d0
               do j = 1, n
                  phipi = phipi - work1(j)**2
               end do
               phipi = phipi / gnleng
            end if
            amulow = -phi / phipi
         else
            first = .false.
            amulow = 0.0d0
         end if
         amuhi = sgnorm / delta
c
c     iterate until a satisfactory "amu" is generated
c
         done = .false.
         do while (.not. done)
            if (amu.lt.amulow .or. amu.gt.amuhi) then
               amu = max(sqrt(amulow*amuhi),0.001d0*amuhi)
            end if
            temp = sqrt(amu)
            do j = 1, n
               work1(j) = temp * xscale(j)
            end do
c
c     solve the damped least squares system for the value of the
c     Levenberg-Marquardt step using More's Minpack technique
c
            call qrsolve (n,m,a,ipvt,work1,qtf,sa,diag,work2)
            do j = 1, n
               sa(j) = -sa(j)
            end do
            do j = 1, n
               work2(j) = xscale(j) * sa(j)
            end do
            stplen = 0.0d0
            do j = 1, n
               stplen = stplen + work2(j)**2
            end do
            stplen = sqrt(stplen)
            phi = stplen - delta
            do j = 1, n
               k = ipvt(j)
               work1(j) = xscale(k) * work2(k)
            end do
            do j = 1, n
               if (abs(diag(j)) .ge. tiny) then
                  work1(j) = work1(j) / diag(j)
               end if
               if (j .lt. n) then
                  do i = j+1, n
                     work1(i) = work1(i) - work1(j)*a(i,j)
                  end do
               end if
            end do
            phip = 0.0d0
            do j = 1, n
               phip = phip - work1(j)**2
            end do
            phip = phip / stplen
c
c     check to see if the step is acceptable; if not,
c     update amulow, amuhi and amu for next iteration
c
            if ((stplen.ge.alow*delta.and.stplen.le.high*delta)
     &            .or. (amuhi-amulow).le.small) then
               done = .true.
            else
               amulow = max(amulow,amu-(phi/phip))
               if (phi .lt. 0.0d0)  amuhi = amu
               amu = amu - (stplen/delta)*(phi/phip)
            end if
         end do
      end if
      deltap = delta
c
c     perform deallocation of some local arrays
c
      deallocate (gnstep)
      deallocate (diag)
      deallocate (work1)
      deallocate (work2)
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine trust  --  update the model trust region  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "trust" updates the model trust region for a nonlinear
c     least squares calculation; based on ideas found in NL2SOL
c     and in Dennis and Schnabel's book
c
c     arguments and variables :
c
c     n         number of least squares variables
c     m         number of residual functions
c     xc        vector with the current iterate
c     fcnorm    scalar containing the norm of f(xc)
c     gc        vector with the gradient at xc
c     a         real m by n matrix containing the upper triangular
c                 matrix r from the QR factorization of the current
c                 Jacobian in the upper triangle
c     ipvt      vector of length n containing the permutation matrix
c                 from QR factorization of the Jacobian
c     sc        vector containing the Newton step
c     sa        vector containing current step
c     xscale    vector containing the diagonal scaling matrix for x
c     gauss     flag set to true when the Gauss-Newton step is taken
c     stpmax    maximum allowable step size
c     delta     trust region radius with value retained between calls
c     icode     return code values, set upon exit
c                 0  means xp accepted as next iterate, delta
c                      is trust region for next iteration
c                 1  means the algorithm was unable to find a
c                      satisfactory xp sufficiently distinct from xc
c                 2  means both the scaled actual and predicted
c                      function reductions are smaller than rftol
c                 3  means that false convergence is detected
c                 4  means fpnorm is too large, current iteration is
c                      continued with a new, reduced trust region
c                 5  means fpnorm is sufficiently small, but the
c                      chance of taking a longer successful step
c                      seems good that the current iteration is to
c                      be continued with a new, doubled trust region
c     xpprev    vector with the value of xp at the  previous call
c                 within this iteration
c     fpprev    vector of length m containing f(xpprev)
c     xp        vector of length n containing the new iterate
c     fp        vector of length m containing the functions at xp
c     fpnorm    scalar containing the norm of f(xp)
c     bigstp    flag set to true if maximum step length was taken
c     ncalls    number of function evaluations used
c     xlo       vector of length n containing the lower bounds
c     xhi       vector of length n containing the upper bounds
c     nactive   number of columns in the active Jacobian
c
c     required external routines :
c
c     rsdvalue   subroutine to evaluate residual function values
c
c
      subroutine trust (n,m,xc,fcnorm,gc,a,ipvt,sc,sa,xscale,gauss,
     &                  stpmax,delta,icode,xp,fc,fp,fpnorm,bigstp,
     &                  ncalls,xlo,xhi,nactive,stpmin,rftol,faketol,
     &                  rsdvalue)
      implicit none
      integer maxlsq,maxrsd
      parameter (maxlsq=1000)
      parameter (maxrsd=1000)
      integer i,j,k
      integer m,n,icode
      integer ncalls,nactive
      integer ipvt(*)
      real*8 fcnorm,stpmax
      real*8 delta,fpnorm,fpnrmp
      real*8 reduce,predict
      real*8 rellen,slope,tiny
      real*8 stplen,stpmin
      real*8 rftol,faketol
      real*8 alpha,temp,precise
      real*8 xpprev(maxlsq)
      real*8 fpprev(maxrsd)
      real*8 xc(*)
      real*8 gc(*)
      real*8 sc(*)
      real*8 sa(*)
      real*8 xp(*)
      real*8 fc(*)
      real*8 fp(*)
      real*8 xlo(*)
      real*8 xhi(*)
      real*8 xscale(*)
      real*8 a(m,*)
      logical gauss,bigstp
      logical feas,ltemp
      save xpprev,fpprev
      save fpnrmp
      external rsdvalue
      external precise
c
c
c     set value of alpha, logicals and step length
c
      alpha = 0.0001d0
      bigstp = .false.
      feas = .true.
      stplen = 0.0d0
      do i = 1, n
         stplen = stplen + (xscale(i)*sc(i))**2
      end do
      stplen = sqrt(stplen)
c
c     compute new trial point and new function values
c
      do i = 1, n
         xp(i) = xc(i) + sc(i)
         if (xp(i) .gt. xhi(i)) then
            sc(i) = xhi(i) - xc(i)
            xp(i) = xhi(i)
            feas = .false.
         else if (xp(i) .lt. xlo(i)) then
            sc(i) = xlo(i) - xc(i)
            xp(i) = xlo(i)
            feas = .false.
         end if
      end do
      ncalls = ncalls + 1
      call rsdvalue (n,m,xp,fp)
      fpnorm = 0.0d0
      do i = 1, m
         fpnorm = fpnorm + fp(i)**2
      end do
      fpnorm = 0.5d0 * fpnorm
      reduce = fpnorm - fcnorm
      slope = 0.0d0
      do i = 1, n
         slope = slope + gc(i)*sc(i)
      end do
      if (icode .ne. 5)  fpnrmp = 0.0d0
c
c     internal doubling no good; reset to previous and quit
c
      if (icode.eq.5 .and.
     &     ((fpnorm.ge.fpnrmp).or.(reduce.gt.alpha*slope))) then
         icode = 0
         do i = 1, n
            xp(i) = xpprev(i)
         end do
         do i = 1, m
            fp(i) = fpprev(i)
         end do
         fpnorm = fpnrmp
         delta = 0.5d0 * delta
c
c     fpnorm is too large; the step is unacceptable
c
      else if (reduce .ge. alpha*slope) then
         rellen = 0.0d0
         do i = 1, n
            temp = abs(sc(i))/max(abs(xp(i)),1.0d0/xscale(i))
            rellen = max(rellen,temp)
         end do
c
c     magnitude of (xp-xc) is too small, end the global step
c
         if (rellen .lt. stpmin) then
            icode = 1
            do i = 1, n
               xp(i) = xc(i)
            end do
            do i = 1, m
               fp(i) = fc(i)
            end do
c
c     quadratic interpolation step; reduce delta and continue
c
         else
            icode = 4
            tiny = precise (1)
            if (abs(reduce-slope) .gt. tiny) then
               temp = -slope*stplen / (2.0d0*(reduce-slope))
            else
               temp = -slope*stplen / 2.0d0
            end if
            if (temp .lt. 0.1d0*delta) then
               delta = 0.1d0 * delta
            else if (temp .gt. 0.5d0*delta) then
               delta = 0.5d0 * delta
            else
               delta = temp
            end if
         end if
c
c     fpnorm is sufficiently small; the step is acceptable compute
c     the predicted reduction as predict = g(T)*s + (1/2)*s(T)*h*s
c     with h = p * r**t * r * p**t
c
      else
         predict = slope
         do i = 1, nactive
            k = ipvt(i)
            temp = 0.0d0
            do j = i, nactive
               temp = temp + sa(k)*a(i,j)
            end do
            predict = predict + 0.5d0*temp**2
         end do
         ltemp = (abs(predict-reduce) .le. 0.1d0*abs(reduce))
c
c     if reduce and predict agree to within relative error of 0.1
c     or if negative curvature is indicated, and a longer step is
c     possible and delta has not been decreased this iteration,
c     then double trust region and continue global step
c
         if (icode.ne.4 .and. (ltemp.or.(reduce.le.slope)) .and. feas
     &        .and. .not.gauss .and. (delta.le.0.99d0*stpmax)) then
            icode = 5
            do i = 1, n
               xpprev(i) = xp(i)
            end do
            do i = 1, m
               fpprev(i) = fp(i)
            end do
            fpnrmp = fpnorm
            delta = min(2.0d0*delta,stpmax)
c
c     accept the point; choose new trust region for next iteration
c
         else
            icode = 0
            if (stplen .gt. 0.99d0*stpmax)  bigstp = .true.
            if (reduce .ge. 0.1d0*predict) then
               delta = 0.5d0 * delta
            else if (reduce .le. 0.75d0*predict) then
               delta = min(2.0d0*delta,stpmax)
            end if
         end if
c
c     check relative function convergence and false convergence
c
         if (reduce .le. 2.0d0*predict) then
            if (abs(reduce).le.rftol*abs(fcnorm) .and.
     &          abs(predict).le.rftol*abs(fcnorm)) then
               icode = 2
            end if
         else
            rellen = 0.0d0
            do i = 1, n
               temp = abs(sc(i))/max(abs(xp(i)),1.0d0/xscale(i))
               rellen = max(rellen,temp)
            end do
            if (rellen .lt. faketol)  icode = 3
         end if
      end if
      return
      end
