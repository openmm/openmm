c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  function random  --  portable random number generator  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "random" generates a random number on [0,1] via a long
c     period generator due to L'Ecuyer with Bays-Durham shuffle
c
c     literature references:
c
c     P. L'Ecuyer, Communications of the ACM, 31, 742-774 (1988)
c
c     W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P.
c     Flannery, Numerical Recipes (Fortran), 2nd Ed., Cambridge
c     University Press, 1992, Section 7.1
c
c
      function random ()
      use inform
      use iounit
      use keys
      implicit none
      integer im1,ia1,iq1,ir1
      integer im2,ia2,iq2,ir2
      integer big,nshuffle
      integer imm1,ndiv
      real*8 factor
      parameter (im1=2147483563)
      parameter (ia1=40014)
      parameter (iq1=53668)
      parameter (ir1=12211)
      parameter (im2=2147483399)
      parameter (ia2=40692)
      parameter (iq2=52774)
      parameter (ir2=3791)
      parameter (big=141803398)
      parameter (nshuffle=32)
      parameter (imm1=im1-1)
      parameter (ndiv=1+imm1/nshuffle)
      parameter (factor=1.0d0/im1)
      integer i,k,iy,next
      integer seed,seed2
      integer year,month,day
      integer hour,minute,second
      integer ishuffle(nshuffle)
      real*8 random
      logical first
      character*20 keyword
      character*240 record
      character*240 string
      save first
      save seed,seed2
      save iy,ishuffle
      data first  / .true. /
c
c
c     random number seed is first set to a big number,
c     then incremented by the seconds elapsed this decade
c
      if (first) then
         first = .false.
         seed = big
         call calendar (year,month,day,hour,minute,second)
         year = mod(year,10)
         seed = seed + 32140800*year + 2678400*(month-1)
         seed = seed + 86400*(day-1) + 3600*hour
         seed = seed + 60*minute + second
c
c     search the keywords for a random number seed
c
         do i = 1, nkey
            next = 1
            record = keyline(i)
            call gettext (record,keyword,next)
            call upcase (keyword)
            if (keyword(1:11) .eq. 'RANDOMSEED ') then
               string = record(next:240)
               read (string,*,err=10)  seed
               seed = max(1,seed)
            end if
   10       continue
         end do
c
c     print the value used for the random number seed
c
         if (verbose) then
            write (iout,20)  seed
   20       format (/,' Random Number Generator Initialized',
     &                 ' with SEED :',3x,i12)
         end if
c
c     warm up and then load the shuffling table
c
         seed2 = seed
         do i = nshuffle+8, 1, -1
            k = seed / iq1
            seed = ia1 * (seed-k*iq1) - k*ir1
            if (seed .lt. 0)  seed = seed + im1
            if (i .le. nshuffle)  ishuffle(i) = seed
         end do
         iy = ishuffle(1)
      end if
c
c     get a new random number value each call
c
      k = seed / iq1
      seed = ia1*(seed-k*iq1) - k*ir1
      if (seed .lt. 0)  seed = seed + im1
      k = seed2 / iq2
      seed2 = ia2*(seed2-k*iq2) - k*ir2
      if (seed2 .lt. 0)  seed2 = seed2 + im2
      i = 1 + iy/ndiv
      iy = ishuffle(i) - seed2
      ishuffle(i) = seed
      if (iy .lt. 1)  iy = iy + imm1
      random = factor * iy
c
c     print the value of the current random number
c
c     if (debug) then
c        write (iout,30)  random
c  30    format (' RANDOM  --  The Random Number Value is',f12.8)
c     end if
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  function normal  --  random number from normal curve  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "normal" generates a random number from a normal Gaussian
c     distribution with a mean of zero and a variance of one
c
c
      function normal ()
      use inform
      use iounit
      implicit none
      real*8 v1,v2,rsq
      real*8 factor,store
      real*8 normal,random
      logical compute
      external random
      save compute,store
      data compute  / .true. /
c
c
c     get a pair of random values from the distribution
c
      if (compute) then
   10    continue
         v1 = 2.0d0*random() - 1.0d0
         v2 = 2.0d0*random() - 1.0d0
         rsq = v1**2 + v2**2
         if (rsq .ge. 1.0d0)  goto 10
         factor = sqrt(-2.0d0*log(rsq)/rsq)
         store = v1 * factor
         normal = v2 * factor
         compute = .false.
c
c     use the second random value computed at the last call
c
      else
         normal = store
         compute = .true.
      end if
c
c     print the value of the current random number
c
c     if (debug) then
c        write (iout,20)  normal
c  20    format (' NORMAL  --  The Random Number Value is',f12.8)
c     end if
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine ranvec  --  unit vector in random direction  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "ranvec" generates a unit vector in 3-dimensional
c     space with uniformly distributed random orientation
c
c     literature references:
c
c     G. Marsaglia, Ann. Math. Stat., 43, 645 (1972)
c
c     R. C. Rapaport, The Art of Molecular Dynamics Simulation,
c     2nd Edition, Cambridge University Press, 2004, Section 18.4
c
c
      subroutine ranvec (vector)
      use inform
      use iounit
      implicit none
      real*8 x,y,s
      real*8 random
      real*8 vector(3)
      external random
c
c
c     get a pair of appropriate components in the plane
c
      s = 2.0d0
      do while (s .ge. 1.0d0)
         x = 2.0d0*random() - 1.0d0
         y = 2.0d0*random() - 1.0d0
         s = x**2 + y**2
      end do
c
c     construct the 3-dimensional random unit vector
c
      vector(3) = 1.0d0 - 2.0d0*s
      s = 2.0d0 * sqrt(1.0d0 - s)
      vector(2) = s * y
      vector(1) = s * x
c
c     print the components of the random unit vector
c
c     if (debug) then
c        write (iout,10)  vector(1),vector(2),vector(3)
c  10    format (' RANVEC  --  The Random Vector is',3f10.4)
c     end if
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine sphere  --  uniform set of points on sphere  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "sphere" finds a specified number of uniformly distributed
c     points on a sphere of unit radius centered at the origin
c
c     literature reference:
c
c     E. B. Saff and A. B. J. Kuijlaars, "Distributing Many
c     Points on a Sphere", The Mathematical Intelligencer,
c     19, 5-11 (1997)
c
c
      subroutine sphere (ndot,dot)
      use math
      implicit none
      integer i,ndot
      real*8 theta,phi
      real*8 h,phiold
      real*8 tot,tot1
      real*8 dot(3,*)
c
c
c     find spherical coordinates then convert to Cartesian
c
      tot = dble(ndot)
      tot1 = dble(ndot-1)
      do i = 1, ndot
         h = -1.0d0 + 2.0d0*dble(i-1)/tot1
         theta = acos(h)
         if (i.eq.1 .or. i.eq.ndot) then
            phi = 0.0d0
         else
            phi = mod(phiold+3.6d0/sqrt(tot*(1.0d0-h*h)),2.0d0*pi)
         end if
         dot(1,i) = sin(theta) * cos(phi)
         dot(2,i) = sin(theta) * sin(phi)
         dot(3,i) = cos(theta)
         phiold = phi
      end do
      return
      end
