c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1998  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  function invbeta  --  inverse Beta distribution function  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "invbeta" computes the inverse Beta distribution function
c     via a combination of Newton iteration and bisection search
c
c     literature reference:
c
c     K. L. Majumder and G. P. Bhattacharjee, "Inverse of the
c     Incomplete Beta Function Ratio", Applied Statistics, 22,
c     411-414 (1973)
c
c
      function invbeta (a,b,y)
      implicit none
      real*8 eps
      parameter (eps=1.0d-5)
      real*8 invbeta,a,b,y
      real*8 x,x0,x1
      real*8 aexp,bexp,beta
      real*8 mean,stdev
      real*8 betai,gammln
      real*8 slope,error
      logical done
      external betai
c
c
c     use limiting values when input argument is out of range
c
      done = .false.
      if (y .le. 0.0d0) then
         x = 0.0d0
         done = .true.
      else if (y .ge. 1.0d0) then
         x = 1.0d0
         done = .true.
      end if
c
c     initial guess from mean and variance of probability function
c
      if (.not. done) then
         aexp = a - 1.0d0
         bexp = b - 1.0d0
         beta = exp(gammln(a) + gammln(b) - gammln(a+b))
         mean = a / (a+b)
         stdev = sqrt(a*b/((a+b+1.0d0)*(a+b)**2))
         if (y.gt.0.0d0 .and. y.le.0.167d0) then
            x = mean + (y/0.167d0-2.0d0)*stdev
         else if (y.gt.0.167d0 .and. y.lt.0.833d0) then
            x = mean + (y/0.333d0-1.5d0)*stdev
         else if (y.ge.0.833d0 .and. y.lt.1.0d0) then
            x = mean + (y/0.167d0-4.0d0)*stdev
         end if
         x = max(eps,min(1.0d0-eps,x))
      end if
c
c     refine inverse distribution value via Newton iteration
c
      do while (.not. done)
         slope = (x**aexp * (1.0d0-x)**bexp) / beta
         error = betai(a,b,x) - y
         x = x - error/slope
         if (abs(error) .lt. eps)  done = .true.
         if (x.lt.0.0d0 .or. x.gt.1.0d0)  done = .true.
      end do
c
c     try bisection search if Newton iteration moved out of range
c
      if (x.lt.0.0d0 .or. x.gt.1.0d0) then
         x0 = 0.0d0
         x1 = 1.0d0
         done = .false.
      end if
c
c     refine inverse distribution value via bisection search
c
      do while (.not. done)
         x = 0.5d0 * (x0+x1)
         error = betai(a,b,x) - y
         if (error .gt. 0.0d0)  x1 = x
         if (error .lt. 0.0d0)  x0 = x
         if (abs(error) .lt. eps)  done = .true.
      end do
c
c     return best estimate of the inverse beta distribution value
c
      invbeta = x
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  function betai  --  cumulative Beta distribution function  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "betai" evaluates the cumulative Beta distribution function
c     as the probability that a random variable from a distribution
c     with Beta parameters "a" and "b" will be less than "x"
c
c
      function betai (a,b,x)
      implicit none
      real*8 betai,a,b,x
      real*8 bt,gammln
      real*8 betacf
      external betacf
c
c
c     get cumulative distribution directly or via reflection
c
      if (x .le. 0.0d0) then
         betai = 0.0d0
      else if (x .ge. 1.0d0) then
         betai = 1.0d0
      else
         bt = exp(gammln(a+b) - gammln(a) - gammln(b)
     &               + a*log(x) + b*log(1.0d0-x))
         if (x .lt. (a+1.0d0)/(a+b+2.0d0)) then
            betai = (bt/a) * betacf (a,b,x)
         else
            betai = 1.0d0 - (bt/b) * betacf (b,a,1.0d0-x)
         end if
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  function betacf  --  continued fraction routine for betai  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "betacf" computes a rapidly convergent continued fraction needed
c     by routine "betai" to evaluate the cumulative Beta distribution
c
c
      function betacf (a,b,x)
      implicit none
      integer maxiter
      real*8 eps,delta
      parameter (maxiter=100)
      parameter (eps=1.0d-10)
      parameter (delta=1.0d-30)
      integer i
      real*8 betacf,a,b,x
      real*8 m,m2,aa
      real*8 c,d,del,h
      real*8 qab,qam,qap
c
c
c     establish an initial guess for the Beta continued fraction
c
      qab = a + b
      qap = a + 1.0d0
      qam = a - 1.0d0
      c = 1.0d0
      d = 1.0d0 - qab*x/qap
      if (abs(d) .lt. delta)  d = delta
      d = 1.0d0 / d
      h = d
c
c     iteratively improve the continued fraction to convergence
c
      do i = 1, maxiter
         m = dble(i)
         m2 = 2.0d0 * m
         aa = m * (b-m) * x / ((qam+m2)*(a+m2))
         d = 1.0d0 + aa*d
         if (abs(d) .lt. delta)  d = delta
         c = 1.0d0 + aa/c
         if (abs(c) .lt. delta)  c = delta
         d = 1.0d0 / d
         h = h * d * c
         aa = -(a+m) * (qab+m) * x / ((a+m2)*(qap+m2))
         d = 1.0d0 + aa*d
         if (abs(d) .lt. delta)  d = delta
         c = 1.0d0 + aa/c
         if (abs(c) .lt. delta)  c = delta
         d = 1.0d0 / d
         del = d * c
         h = h * del
         if (abs(del-1.0d0) .lt. eps)  goto 10
      end do
   10 continue
      betacf = h
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  function gammln  --  natural log of the Gamma function  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "gammln" uses a series expansion due to Lanczos to compute
c     the natural logarithm of the Gamma function at "x" in [0,1]
c
c
      function gammln (x)
      implicit none
      real*8 step,c0,c1,c2,c3,c4,c5,c6
      parameter (step=2.5066282746310005d0)
      parameter (c0=1.000000000190015d0)
      parameter (c1=7.618009172947146d1)
      parameter (c2=-8.650532032941677d1)
      parameter (c3=2.401409824083091d1)
      parameter (c4=-1.231739572450155d0)
      parameter (c5=1.208650973866179d-3)
      parameter (c6=-5.395239384953d-6)
      real*8 gammln,x
      real*8 series,temp
c
c
c     get the natural log of Gamma via a series expansion
c
      temp = x + 5.5d0
      temp = (x+0.5d0)*log(temp) - temp
      series = c0 + c1/(x+1.0d0) + c2/(x+2.0d0) + c3/(x+3.0d0)
     &            + c4/(x+4.0d0) + c5/(x+5.0d0) + c6/(x+6.0d0)
      gammln = temp + log(step*series/x)
      return
      end
