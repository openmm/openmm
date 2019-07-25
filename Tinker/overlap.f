c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine overlap  --  p-orbital overlap for pisystem  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "overlap" computes the overlap for two parallel p-orbitals
c     given the atomic numbers and distance of separation
c
c
      subroutine overlap (atmnum1,atmnum2,rang,ovlap)
      use units
      implicit none
      integer atmnum1
      integer atmnum2
      integer na,nb,la,lb
      real*8 ovlap
      real*8 rbohr,rang
      real*8 za,zb,s(3)
      real*8 zeta(18)
      save zeta
c
c     Slater orbital exponents for hydrogen through argon
c
      data zeta  / 1.000, 1.700, 0.650, 0.975, 1.300, 1.625,
     &             1.950, 2.275, 2.600, 2.925, 0.733, 0.950,
     &             1.167, 1.383, 1.600, 1.817, 2.033, 2.250 /
c
c
c     principal quantum number from atomic number
c
      na = 2
      nb = 2
      if (atmnum1 .gt. 10)  na = 3
      if (atmnum2 .gt. 10)  nb = 3
c
c     azimuthal quantum number for p-orbitals
c
      la = 1
      lb = 1
c
c     orbital exponent from stored ideal values
c
      za = zeta(atmnum1)
      zb = zeta(atmnum2)
c
c     convert interatomic distance to bohrs
c
      rbohr = rang / bohr
c
c     get pi-overlap via generic overlap integral routine
c
      call slater (na,la,za,nb,lb,zb,rbohr,s)
      ovlap = s(2)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine slater  --  find overlap integrals for STO's  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "slater" is a general routine for computing the overlap
c     integrals between two Slater-type orbitals
c
c     literature reference:
c
c     D. B. Cook, "Structures and Approximations for Electrons in
c     Molecules", Ellis Horwood Limited, Sussex, England, 1978,
c     adapted from the code in Chapter 7
c
c     variables and parameters:
c
c     na   principle quantum number for first orbital
c     la   azimuthal quantum number for first orbital
c     za   orbital exponent for the first orbital
c     nb   principle quantum number for second orbital
c     lb   azimuthal quantum number for second orbital
c     zb   orbital exponent for the second orbital
c     r    interatomic distance in atomic units
c     s    vector containing the sigma-sigma, pi-pi
c            and delta-delta overlaps upon output
c
c
      subroutine slater (na,la,za,nb,lb,zb,r,s)
      implicit none
      real*8 rmin,eps
      parameter (rmin=0.000001d0)
      parameter (eps=0.00000001d0)
      integer j,k,m,na,nb,la,lb,ja,jb
      integer nn,max,maxx,novi
      integer idsga(5),idsgb(5)
      integer icosa(2),icosb(2)
      integer isina(4),isinb(4)
      integer ia(200),ib(200)
      real*8 an,ana,anb,anr
      real*8 rhalf,coef,p,pt
      real*8 r,za,zb,cjkm
      real*8 s(3),fact(15)
      real*8 cbase(20),theta(6)
      real*8 cosa(2),cosb(2)
      real*8 sinab(4)
      real*8 dsiga(5),dsigb(5)
      real*8 a(20),b(20),c(200)
      logical done
      save icosa,icosb
      save cosa,cosb
      save idsga,idsgb
      save dsiga,dsigb
      save isina,isinb
      save sinab,theta,fact
      external cjkm
      data icosa  / 0, 1 /
      data icosb  / 0, 1 /
      data cosa   /  1.0d0, 1.0d0 /
      data cosb   / -1.0d0, 1.0d0 /
      data idsga  / 0, 1, 2, 2, 0 /
      data idsgb  / 0, 1, 2, 0, 2 /
      data dsiga  / 3.0d0, 4.0d0, 3.0d0, -1.0d0, -1.0d0 /
      data dsigb  / 3.0d0,-4.0d0, 3.0d0, -1.0d0, -1.0d0 /
      data isina  / 0, 2, 0, 2 /
      data isinb  / 0, 0, 2, 2 /
      data sinab  / -1.0d0, 1.0d0, 1.0d0, -1.0d0 /
      data theta  / 0.7071068d0, 1.2247450d0, 0.8660254d0,
     &              0.7905694d0, 1.9364916d0, 0.9682458d0 /
      data fact   / 1.0d0, 1.0d0, 2.0d0, 6.0d0, 24.0d0, 120.0d0,
     &              720.0d0, 5040.0d0, 40320.0d0, 362880.0d0,
     &              3628800.0d0, 39916800.0d0, 479001600.0d0,
     &              6227020800.0d0, 87178291200.0d0 /
c
c
c     zero out the overlap integrals
c
      done = .false.
      s(1) = 0.0d0
      s(2) = 0.0d0
      s(3) = 0.0d0
      ana = (2.0d0*za)**(2*na+1) / fact(2*na+1)
      anb = (2.0d0*zb)**(2*nb+1) / fact(2*nb+1)
c
c     orbitals are on the same atomic center
c
      if (r .lt. rmin) then
         anr = 1.0d0
         j = na + nb + 1
         s(1) = fact(j) / ((za+zb)**j)
         an = sqrt(ana*anb)
         do novi = 1, 3
            s(novi) = s(novi) * an * anr
         end do
         return
      end if
c
c     compute overlap integrals for general case
c
      rhalf = 0.5d0 * r
      p = rhalf * (za+zb)
      pt = rhalf * (za-zb)
      nn = na + nb
      call aset (p,nn,a)
      call bset (pt,nn,b)
      k = na - la
      m = nb - lb
      max = k + m + 1
      do j = 1, max
         ia(j) = j - 1
         ib(j) = max - j
         cbase(j) = cjkm(j-1,k,m)
         c(j) = cbase(j)
      end do
      maxx = max
      if (la .eq. 1) then
         call polyp (c,ia,ib,maxx,cosa,icosa,icosb,2)
      else if (la .eq. 2) then
         call polyp (c,ia,ib,maxx,dsiga,idsga,idsgb,5)
      end if
      if (lb .eq. 1) then
         call polyp (c,ia,ib,maxx,cosb,icosa,icosb,2)
      else if (lb .eq. 2) then
         call polyp (c,ia,ib,maxx,dsigb,idsga,idsgb,5)
      end if
      novi = 1
      do while (.not. done)
         do j = 1, maxx
            ja = ia(j) + 1
            jb = ib(j) + 1
            coef = c(j)
            if (abs(coef) .ge. eps) then
               s(novi) = s(novi) + coef*a(ja)*b(jb)
            end if
         end do
         ja = la*(la+1)/2 + novi
         jb = lb*(lb+1)/2 + novi
         s(novi) = s(novi) * theta(ja) * theta(jb)
         if (novi.eq.1 .and. la.ne.0 .and. lb.ne.0) then
            maxx = max
            do j = 1, maxx
               c(j) = cbase(j)
            end do
            call polyp (c,ia,ib,maxx,sinab,isina,isinb,4)
            if (la .eq. 2) then
               call polyp (c,ia,ib,maxx,cosa,icosa,icosb,2)
            end if
            if (lb .eq. 2) then
               call polyp (c,ia,ib,maxx,cosb,icosa,icosb,2)
            end if
            novi = 2
         else if (novi.eq.2 .and. la.eq.2 .and. lb.eq.2) then
            maxx = max
            do j = 1, maxx
               c(j) = cbase(j)
            end do
            call polyp (c,ia,ib,maxx,sinab,isina,isinb,4)
            call polyp (c,ia,ib,maxx,sinab,isina,isinb,4)
            novi = 3
         else
            anr = rhalf**(na+nb+1)
            an = sqrt(ana*anb)
            do novi = 1, 3
               s(novi) = s(novi) * an * anr
            end do
            done = .true.
         end if
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine polyp  --  polynomial product for STO overlap  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "polyp" is a polynomial product routine that multiplies two
c     algebraic forms
c
c
      subroutine polyp (c,ia,ib,max,d,iaa,ibb,n)
      implicit none
      integer i,j,k,m,max,n
      integer ia(200),ib(200)
      integer iaa(*),ibb(*)
      real*8 c(200),d(*)
c
c
      do j = 1, max
         do k = 1, n
            i = n - k + 1
            m = (i-1)*max + j
            c(m) = c(j) * d(i)
            ia(m) = ia(j) + iaa(i)
            ib(m) = ib(j) + ibb(i)
         end do
      end do
      max = n * max
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  function cjkm  --  coefficients of spherical harmonics  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "cjkm" computes the coefficients of spherical harmonics
c     expressed in prolate spheroidal coordinates
c
c
      function cjkm (j,k,m)
      implicit none
      integer i,j,k,m
      integer min,max
      integer id,idd,ip1
      real*8 cjkm,b1,b2,sum
      real*8 fact(15)
      save fact
      data fact  / 1.0d0, 1.0d0, 2.0d0, 6.0d0, 24.0d0, 120.0d0,
     &             720.0d0, 5040.0d0, 40320.0d0, 362880.0d0,
     &             3628800.0d0, 39916800.0d0, 479001600.0d0,
     &             6227020800.0d0, 87178291200.0d0 /
c
c
      min = 1
      if (j .gt. m)  min = j - m + 1
      max = j + 1
      if (k .lt. j)  max = k + 1
      sum = 0.0d0
      do ip1 = min, max
         i = ip1 - 1
         id = k - i + 1
         b1 = fact(k+1) / (fact(i+1)*fact(id))
         if (j .lt. i) then
            b2 = 1.0d0
         else
            id = m - (j-i) + 1
            idd = j - i + 1
            b2 = fact(m+1) / (fact(idd)*fact(id))
         end if
         sum = sum + b1*b2*(-1.0d0)**i
      end do
      cjkm = sum * (-1.0d0)**(m-j)
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine aset  --  get "A" functions by recursion  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "aset" computes by recursion the A functions used in the
c     evaluation of Slater-type (STO) overlap integrals
c
c
      subroutine aset (alpha,n,a)
      implicit none
      integer i,n
      real*8 alpha,alp
      real*8 a(20)
c
c
      alp = 1.0d0 / alpha
      a(1) = exp(-alpha) * alp
      do i = 1, n
         a(i+1) = a(1) + dble(i)*a(i)*alp
      end do
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine bset  --  get "B" functions by recursion  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "bset" computes by downward recursion the B functions used
c     in the evaluation of Slater-type (STO) overlap integrals
c
c
      subroutine bset (beta,n,b)
      implicit none
      real*8 eps
      parameter (eps=0.000001d0)
      integer i,j,n
      real*8 beta,bmax
      real*8 betam,d1,d2
      real*8 b(20)
      external bmax
c
c
      if (abs(beta) .lt. eps) then
         do i = 1, n+1
            b(i) = 2.0d0 / dble(i)
            if ((i/2)*2 .eq. i)  b(i) = 0.0d0
         end do
      else if (abs(beta) .gt. (dble(n)/2.3d0)) then
         d1 = exp(beta)
         d2 = 1.0d0 / d1
         betam = 1.0d0 / beta
         b(1) = (d1-d2) * betam
         do i = 1, n
            d1 = -d1
            b(i+1) = (d1-d2+dble(i)*b(i)) * betam
         end do
      else
         b(n+1) = bmax(beta,n)
         d1 = exp(beta)
         d2 = 1.0d0 / d1
         if ((n/2)*2 .ne. n)  d1 = -d1
         do i = 1, n
            j = n - i + 1
            d1 = -d1
            b(j) = (d1+d2+beta*b(j+1)) / dble(j)
         end do
      end if
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  function bmax  --  find maximum order of "B" functions  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "bmax" computes the maximum order of the B functions needed
c     for evaluation of Slater-type (STO) overlap integrals
c
c
      function bmax (beta,n)
      implicit none
      real*8 eps
      parameter (eps=0.0000001d0)
      integer n
      real*8 bmax,beta
      real*8 b,top,bot
      real*8 sum,fi
      real*8 sign,term
      logical done
c
c
      done = .false.
      b = beta**2
      top = dble(n) + 1.0d0
      sum = 1.0d0 / top
      fi = 2.0d0
      sign = 2.0d0
      if ((n/2)*2 .ne. n) then
         top = top + 1.0d0
         sum = beta / top
         fi = fi + 1.0d0
         sign = -2.0d0
      end if
      term = sum
      do while (.not. done)
         bot = top + 2.0d0
         term = term * b * top / (fi*(fi-1.0d0)*bot)
         sum = sum + term
         if (abs(term) .le. eps) then
            done = .true.
         else
            fi = fi + 2.0d0
            top = bot
         end if
      end do
      bmax = sign * sum
      return
      end
