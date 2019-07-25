c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1996  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine born  --  Born radii for implicit solvation  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "born" computes the Born radius of each atom for use with
c     the various implicit solvation models
c
c     literature references:
c
c     W. C. Still, A. Tempczyk, R. C. Hawley and T. Hendrickson,
c     "A Semianalytical Treatment of Solvation for Molecular
c     Mechanics and Dynamics", J. Amer. Chem. Soc., 112, 6127-6129
c     (1990)  ("Onion" Method; see supplimentary material)
c
c     D. Qiu, P. S. Shenkin, F. P. Hollinger and W. C. Still, "The
c     GB/SA Continuum Model for Solvation. A Fast Analytical Method
c     for the Calculation of Approximate Radii", J. Phys. Chem. A,
c     101, 3005-3014 (1997)  (Analytical Still Method)
c
c     G. D. Hawkins, C. J. Cramer and D. G. Truhlar, "Parametrized
c     Models of Aqueous Free Energies of Solvation Based on Pairwise
c     Descreening of Solute Atomic Charges from a Dielectric Medium",
c     J. Phys. Chem., 100, 19824-19839 (1996)  (HCT Method)
c
c     A. Onufriev, D. Bashford and D. A. Case, "Exploring Protein
c     Native States and Large-Scale Conformational Changes with a
c     Modified Generalized Born Model", PROTEINS, 55, 383-394 (2004)
c     (OBC Method)
c
c     T. Grycuk, "Deficiency of the Coulomb-field Approximation
c     in the Generalized Born Model: An Improved Formula for Born
c     Radii Evaluation", J. Chem. Phys., 119, 4817-4826 (2003)
c     (Grycuk Method)
c
c     M. Schaefer, C. Bartels and M. Karplus, "Solution Conformations
c     and Thermodynamics of Structured Peptides: Molecular Dynamics
c     Simulation with an Implicit Solvation Model", J. Mol. Biol.,
c     284, 835-848 (1998)  (ACE Method)
c
c
      subroutine born
      use sizes
      use atomid
      use atoms
      use bath
      use chgpot
      use couple
      use inform
      use iounit
      use math
      use pbstuf
      use solute
      implicit none
      integer i,j,k,it,kt
      integer, allocatable :: skip(:)
      real*8 area,rold,t
      real*8 shell,fraction
      real*8 inner,outer,tinit
      real*8 ratio,total
      real*8 xi,yi,zi,ri
      real*8 rk,sk,sk2
      real*8 lik,lik2
      real*8 uik,uik2
      real*8 tsum,tchain
      real*8 sum,sum2,sum3
      real*8 alpha,beta,gamma
      real*8 xr,yr,zr,rvdw
      real*8 r,r2,r3,r4
      real*8 gpi,pip5,p5inv
      real*8 theta,term,ccf
      real*8 l2,l4,lr,l4r
      real*8 u2,u4,ur,u4r
      real*8 expterm,rmu
      real*8 b0,gself
      real*8 third,pi43
      real*8 bornmax
      real*8, allocatable :: roff(:)
      real*8, allocatable :: pos(:,:)
      real*8, allocatable :: pbpole(:,:)
      logical done
c
c
c     perform dynamic allocation of some local arrays
c
      if (borntyp .eq. 'STILL')  allocate (skip(n))
      allocate (roff(n))
      if (borntyp .eq. 'PERFECT') then
         allocate (pos(3,n))
         allocate (pbpole(13,n))
      end if
c
c     perform dynamic allocation of some global arrays
c
      if (borntyp .eq. 'PERFECT') then
         if (.not. allocated(apbe))  allocate (apbe(n))
         if (.not. allocated(pbep))  allocate (pbep(3,n))
         if (.not. allocated(pbfp))  allocate (pbfp(3,n))
         if (.not. allocated(pbtp))  allocate (pbtp(3,n))
      end if
c
c     set offset modified radii and OBC chain rule factor
c
      do i = 1, n
         roff(i) = rsolv(i) - doffset
         drobc(i) = 1.0d0
      end do
c
c     get the Born radii via the numerical "Onion" method
c
      if (borntyp .eq. 'ONION') then
         tinit = 0.1d0
         ratio = 1.5d0
         do i = 1, n
            t = tinit
            rold = roff(i)
            total = 0.0d0
            done = .false.
            do while (.not. done)
               roff(i) = roff(i) + 0.5d0*t
               call surfatom (i,area,roff)
               fraction = area / (4.0d0*pi*roff(i)**2)
               if (fraction .lt. 0.99d0) then
                  inner = roff(i) - 0.5d0*t
                  outer = inner + t
                  shell = 1.0d0/inner - 1.0d0/outer
                  total = total + fraction*shell
                  roff(i) = roff(i) + 0.5d0*t
                  t = ratio * t
               else
                  inner = roff(i) - 0.5d0*t
                  total = total + 1.0d0/inner
                  done = .true.
               end if
            end do
            rborn(i) = 1.0d0 / total
            roff(i) = rold
         end do
c
c     get the Born radii via the analytical Still method;
c     note this code only loops over the variable parts
c
      else if (borntyp .eq. 'STILL') then
         do i = 1, n
            skip(i) = 0
         end do
         p5inv = 1.0d0 / p5
         pip5 = pi * p5
         do i = 1, n
            xi = x(i)
            yi = y(i)
            zi = z(i)
            gpi = gpol(i)
            skip(i) = i
            do j = 1, n12(i)
               skip(i12(j,i)) = i
            end do
            do j = 1, n13(i)
               skip(i13(j,i)) = i
            end do
            do k = 1, n
               if (skip(k) .ne. i) then
                  xr = x(k) - xi
                  yr = y(k) - yi
                  zr = z(k) - zi
                  r2 = xr**2 + yr**2 + zr**2
                  r4 = r2 * r2
                  rvdw = rsolv(i) + rsolv(k)
                  ratio = r2 / (rvdw*rvdw)
                  if (ratio .gt. p5inv) then
                     ccf = 1.0d0
                  else
                     theta = ratio * pip5
                     term = 0.5d0 * (1.0d0-cos(theta))
                     ccf = term * term
                  end if
                  gpi = gpi + p4*ccf*vsolv(k)/r4
               end if
            end do
            rborn(i) = -0.5d0 * electric / gpi
         end do
c
c     get the Born radii via the Hawkins-Cramer-Truhlar method
c
      else if (borntyp .eq. 'HCT') then
         do i = 1, n
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ri = roff(i)
            sum = 1.0d0 / ri
            do k = 1, n
               if (i .ne. k) then
                  xr = x(k) - xi
                  yr = y(k) - yi
                  zr = z(k) - zi
                  r2 = xr**2 + yr**2 + zr**2
                  r = sqrt(r2)
                  rk = roff(k)
                  sk = rk * shct(k)
                  sk2 = sk * sk
                  if (ri .lt. r+sk) then
                     lik = 1.0d0 / max(ri,abs(r-sk))
                     uik = 1.0d0 / (r+sk)
                     lik2 = lik * lik
                     uik2 = uik * uik
                     term = lik - uik + 0.25d0*r*(uik2-lik2)
     &                         + (0.5d0/r)*log(uik/lik)
     &                         + (0.25d0*sk2/r)*(lik2-uik2)
                     if (ri .lt. sk-r) then
                        term = term + 2.0d0*(1.0d0/ri-lik)
                     end if
                     sum = sum - 0.5d0*term
                  end if
               end if
            end do
            rborn(i) = 1.0d0 / sum
         end do
c
c     get the Born radii via the Onufriev-Bashford-Case method
c
      else if (borntyp .eq. 'OBC') then
         do i = 1, n
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ri = roff(i)
            sum = 0.0d0
            do k = 1, n
               if (i .ne. k) then
                  xr = x(k) - xi
                  yr = y(k) - yi
                  zr = z(k) - zi
                  r2 = xr**2 + yr**2 + zr**2
                  r = sqrt(r2)
                  rk = roff(k)
                  sk = rk * shct(k)
                  sk2 = sk * sk
                  if (ri .lt. r+sk) then
                     lik = 1.0d0 / max(ri,abs(r-sk))
                     uik = 1.0d0 / (r+sk)
                     lik2 = lik * lik
                     uik2 = uik * uik
                     term = lik - uik + 0.25d0*r*(uik2-lik2)
     &                         + (0.5d0/r)*log(uik/lik)
     &                         + (0.25d0*sk2/r)*(lik2-uik2)
                     if (ri .lt. sk-r) then
                        term = term + 2.0d0*(1.0d0/ri-lik)
                     end if
                     sum = sum + 0.5d0*term
                  end if
               end if
            end do
            alpha = aobc(i)
            beta = bobc(i)
            gamma = gobc(i)
            sum = ri * sum
            sum2 = sum * sum
            sum3 = sum * sum2
            tsum = tanh(alpha*sum - beta*sum2 + gamma*sum3)
            rborn(i) = 1.0d0/ri - tsum/rsolv(i)
            rborn(i) = 1.0d0 / rborn(i)
            tchain = ri * (alpha-2.0d0*beta*sum+3.0d0*gamma*sum2)
            drobc(i) = (1.0d0-tsum*tsum) * tchain / rsolv(i)
         end do
c
c     get the Born radii via Grycuk's modified HCT method
c
      else if (borntyp .eq. 'GRYCUK') then
         third = 1.0d0 / 3.0d0
         pi43 = 4.0d0 * third * pi
         do i = 1, n
            rborn(i) = 0.0d0
            ri = rsolv(i)
            if (ri .gt. 0.0d0) then
               xi = x(i)
               yi = y(i)
               zi = z(i)
               sum = pi43 / ri**3
               do k = 1, n
                  rk = rsolv(k)
                  if (i.ne.k .and. rk.gt.0.0d0) then
                     xr = x(k) - xi
                     yr = y(k) - yi
                     zr = z(k) - zi
                     r2 = xr**2 + yr**2 + zr**2
                     r = sqrt(r2)
                     sk = rk * shct(k)
                     sk2 = sk * sk
                     if (ri+r .lt. sk) then
                        lik = ri
                        uik = sk - r
                        sum = sum + pi43*(1.0d0/uik**3-1.0d0/lik**3)
                     end if
                     uik = r + sk
                     if (ri+r .lt. sk) then
                        lik = sk - r
                     else if (r .lt. ri+sk) then
                        lik = ri
                     else
                        lik = r - sk
                     end if
                     l2 = lik * lik
                     l4 = l2 * l2
                     lr = lik * r
                     l4r = l4 * r
                     u2 = uik * uik
                     u4 = u2 * u2
                     ur = uik * r
                     u4r = u4 * r
                     term = (3.0d0*(r2-sk2)+6.0d0*u2-8.0d0*ur)/u4r
     &                         - (3.0d0*(r2-sk2)+6.0d0*l2-8.0d0*lr)/l4r
                     sum = sum - pi*term/12.0d0
                  end if
               end do
               rborn(i) = (sum/pi43)**third
               if (rborn(i) .le. 0.0d0)  rborn(i) = 0.0001d0
               rborn(i) = 1.0d0 / rborn(i)
            end if
         end do
c
c     get the Born radii via analytical continuum electrostatics
c
      else if (borntyp .eq. 'ACE') then
         third = 1.0d0 / 3.0d0
         b0 = 0.0d0
         do i = 1, n
            b0 = b0 + vsolv(i)
         end do
         b0 = (0.75d0*b0/pi)**third
         do i = 1, n
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ri = rsolv(i)
            it = class(i)
            gself = 1.0d0/ri + 2.0d0*wace(it,it)
            do k = 1, n
               if (k .ne. i) then
                  xr = x(k) - xi
                  yr = y(k) - yi
                  zr = z(k) - zi
                  kt = class(k)
                  r2 = xr**2 + yr**2 + zr**2
                  r3 = r2 * sqrt(r2)
                  r4 = r2 * r2
                  expterm = wace(it,kt) * exp(-r2/s2ace(it,kt))
                  rmu = r4 + uace(it,kt)**4
                  term = (vsolv(k)/(8.0d0*pi)) * (r3/rmu)**4
                  gself = gself - 2.0d0*(expterm+term)
               end if
            end do
            if (gself .ge. 0.5d0/b0) then
               rborn(i) = 1.0d0 / gself
            else
               rborn(i) = 2.0d0 * b0 * (1.0d0+b0*gself)
            end if
         end do
c
c     get the "perfect" Born radii via Poisson-Boltzmann
c
      else if (borntyp .eq. 'PERFECT') then
         do i = 1, n
            pos(1,i) = x(i)
            pos(2,i) = y(i)
            pos(3,i) = z(i)
            do j = 1, 13
               pbpole(j,i) = 0.0d0
            end do
         end do
         term = -0.5d0 * electric * (1.0d0-1.0d0/sdie)
         do i = 1, n
            pbpole(1,i) = 1.0d0
            call apbsempole (n,pos,rsolv,pbpole,pbe,apbe,pbep,pbfp,pbtp)
            pbpole(1,i) = 0.0d0
            rborn(i) = term / pbe
         end do
      end if
c
c     perform deallocation of some local arrays
c
      if (borntyp .eq. 'STILL')  deallocate (skip)
      deallocate (roff)
      if (borntyp .eq. 'PERFECT') then
         deallocate (pos)
         deallocate (pbpole)
      end if
c
c     make sure the final values are in a reasonable range
c
      bornmax = 500.0d0
      do i = 1, n
         if (rborn(i).lt.0.0d0 .or. rborn(i).gt.bornmax)
     &      rborn(i) = bornmax
      end do
c
c     write out the final Born radius value for each atom
c
      if (debug) then
         write (iout,10)
   10    format (/,' Born Radii for Individual Atoms :',/)
         k = 1
         do while (k .le. n)
            write (iout,20)  (i,rborn(i),i=k,min(k+4,n))
   20       format (1x,5(i7,f8.3))
            k = k + 5
         end do
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine born1  --  Born radii chain rule derivatives  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "born1" computes derivatives of the Born radii with respect
c     to atomic coordinates and increments total energy derivatives
c     and virial components for potentials involving Born radii
c
c
      subroutine born1
      use sizes
      use atomid
      use atoms
      use chgpot
      use couple
      use deriv
      use math
      use solute
      use virial
      implicit none
      integer i,j,k,it,kt
      integer, allocatable :: skip(:)
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 de,de1,de2
      real*8 r,r2,r3,r4,r6
      real*8 p5inv,pip5
      real*8 gpi,vk,ratio
      real*8 ccf,cosq,dccf
      real*8 sinq,theta
      real*8 factor,term
      real*8 rb2,ri,rk
      real*8 sk,sk2
      real*8 lik,lik2,lik3
      real*8 uik,uik2,uik3
      real*8 dlik,duik
      real*8 t1,t2,t3
      real*8 rbi,rbi2,vi
      real*8 ws2,s2ik,uik4
      real*8 third,pi43
      real*8 dbr,dborn
      real*8 expterm,rusum
      real*8 dedx,dedy,dedz
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8, allocatable :: roff(:)
      logical use_gk
c
c
c     perform dynamic allocation of some local arrays
c
      if (borntyp .eq. 'STILL')  allocate (skip(n))
      allocate (roff(n))
c
c     compute atomic radii modified by the dielectric offset
c
      do i = 1, n
         roff(i) = rsolv(i) - doffset
      end do
c
c     set flag for use of generalized Kirkwood with polarization
c
      use_gk = .false.
      if (solvtyp(1:2) .eq. 'GK')  use_gk = .true.
c
c     get Born radius chain rule components for the Still method
c
      if (borntyp .eq. 'STILL') then
         p5inv = 1.0d0 / p5
         pip5 = pi * p5
         do i = 1, n
            skip(i) = 0
         end do
         do i = 1, n
            xi = x(i)
            yi = y(i)
            zi = z(i)
            skip(i) = i
            do j = 1, n12(i)
               skip(i12(j,i)) = i
            end do
            do j = 1, n13(i)
               skip(i13(j,i)) = i
            end do
            gpi = 2.0d0 * rborn(i)**2 / electric
            do k = 1, n
               if (skip(k) .ne. i) then
                  xr = x(k) - xi
                  yr = y(k) - yi
                  zr = z(k) - zi
                  vk = vsolv(k)
                  r2 = xr**2 + yr**2 + zr**2
                  r =  sqrt(r2)
                  r6 = r2 * r2 * r2
                  ratio = r2 / (rsolv(i)+rsolv(k))**2
                  if (ratio .gt. p5inv) then
                     ccf = 1.0d0
                     dccf = 0.0d0
                  else
                     theta = ratio * pip5
                     cosq = cos(theta)
                     term = 0.5d0 * (1.0d0-cosq)
                     ccf = term * term
                     sinq = sin(theta)
                     dccf = 2.0d0 * term * sinq * pip5 * ratio
                  end if
                  dborn = drb(i)
                  if (use_gk)  dborn = dborn + drbp(i)
                  de = dborn * p4 * gpi * vk * (4.0d0*ccf-dccf)/r6
c
c     increment the overall implicit solvation derivatives
c
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
                  des(1,i) = des(1,i) + dedx
                  des(2,i) = des(2,i) + dedy
                  des(3,i) = des(3,i) + dedz
                  des(1,k) = des(1,k) - dedx
                  des(2,k) = des(2,k) - dedy
                  des(3,k) = des(3,k) - dedz
c
c     increment the internal virial tensor components
c
                  vxx = xr * dedx
                  vyx = yr * dedx
                  vzx = zr * dedx
                  vyy = yr * dedy
                  vzy = zr * dedy
                  vzz = zr * dedz
                  vir(1,1) = vir(1,1) + vxx
                  vir(2,1) = vir(2,1) + vyx
                  vir(3,1) = vir(3,1) + vzx
                  vir(1,2) = vir(1,2) + vyx
                  vir(2,2) = vir(2,2) + vyy
                  vir(3,2) = vir(3,2) + vzy
                  vir(1,3) = vir(1,3) + vzx
                  vir(2,3) = vir(2,3) + vzy
                  vir(3,3) = vir(3,3) + vzz
               end if
            end do
         end do
c
c     get Born radius chain rule components for the HCT method
c
      else if (borntyp .eq. 'HCT') then
         do i = 1, n
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ri = roff(i)
            rb2 = rborn(i) * rborn(i)
            do k = 1, n
               if (k .ne. i) then
                  xr = x(k) - xi
                  yr = y(k) - yi
                  zr = z(k) - zi
                  rk = roff(k)
                  sk = rk * shct(k)
                  sk2 = sk * sk
                  r2 = xr**2 + yr**2 + zr**2
                  r = sqrt(r2)
                  if (ri .lt. r+sk) then
                     lik = 1.0d0 / max(ri,abs(r-sk))
                     uik = 1.0d0 / (r+sk)
                     lik2 = lik * lik
                     uik2 = uik * uik
                     lik3 = lik * lik2
                     uik3 = uik * uik2
                     dlik = 1.0d0
                     if (ri .ge. r-sk)  dlik = 0.0d0
                     duik = 1.0d0
                     t1 = 0.5d0*lik2 + 0.25d0*sk2*lik3/r
     &                       - 0.25d0*(lik/r+lik3*r)
                     t2 = -0.5d0*uik2 - 0.25d0*sk2*uik3/r
     &                       + 0.25d0*(uik/r+uik3*r)
                     t3 = 0.125d0*(1.0d0+sk2/r2)*(lik2-uik2)
     &                       + 0.25d0*log(uik/lik)/r2
                     dborn = drb(i)
                     if (use_gk)  dborn = dborn + drbp(i)
                     de = dborn * rb2 * (dlik*t1+duik*t2+t3) / r
c
c     increment the overall implicit solvation derivatives
c
                     dedx = de * xr
                     dedy = de * yr
                     dedz = de * zr
                     des(1,i) = des(1,i) + dedx
                     des(2,i) = des(2,i) + dedy
                     des(3,i) = des(3,i) + dedz
                     des(1,k) = des(1,k) - dedx
                     des(2,k) = des(2,k) - dedy
                     des(3,k) = des(3,k) - dedz
c
c     increment the internal virial tensor components
c
                     vxx = xr * dedx
                     vyx = yr * dedx
                     vzx = zr * dedx
                     vyy = yr * dedy
                     vzy = zr * dedy
                     vzz = zr * dedz
                     vir(1,1) = vir(1,1) + vxx
                     vir(2,1) = vir(2,1) + vyx
                     vir(3,1) = vir(3,1) + vzx
                     vir(1,2) = vir(1,2) + vyx
                     vir(2,2) = vir(2,2) + vyy
                     vir(3,2) = vir(3,2) + vzy
                     vir(1,3) = vir(1,3) + vzx
                     vir(2,3) = vir(2,3) + vzy
                     vir(3,3) = vir(3,3) + vzz
                  end if
               end if
            end do
         end do
c
c     get Born radius chain rule components for the OBC method
c
      else if (borntyp .eq. 'OBC') then
         do i = 1, n
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ri = roff(i)
            rb2 = rborn(i) * rborn(i) * drobc(i)
            do k = 1, n
               if (k .ne. i) then
                  xr = x(k) - xi
                  yr = y(k) - yi
                  zr = z(k) - zi
                  rk = roff(k)
                  sk = rk * shct(k)
                  sk2 = sk * sk
                  r2 = xr**2 + yr**2 + zr**2
                  r = sqrt(r2)
                  if (ri .lt. r+sk) then
                     lik = 1.0d0 / max(ri,abs(r-sk))
                     uik = 1.0d0 / (r+sk)
                     lik2 = lik * lik
                     uik2 = uik * uik
                     lik3 = lik * lik2
                     uik3 = uik * uik2
                     dlik = 1.0d0
                     if (ri .ge. r-sk)  dlik = 0.0d0
                     duik = 1.0d0
                     t1 = 0.5d0*lik2 + 0.25d0*sk2*lik3/r
     &                       - 0.25d0*(lik/r+lik3*r)
                     t2 = -0.5d0*uik2 - 0.25d0*sk2*uik3/r
     &                       + 0.25d0*(uik/r+uik3*r)
                     t3 = 0.125d0*(1.0d0+sk2/r2)*(lik2-uik2)
     &                       + 0.25d0*log(uik/lik)/r2
                     dborn = drb(i)
                     if (use_gk)  dborn = dborn + drbp(i)
                     de = dborn * rb2 * (dlik*t1+duik*t2+t3) / r
c
c     increment the overall permanent solvation derivatives
c
                     dedx = de * xr
                     dedy = de * yr
                     dedz = de * zr
                     des(1,i) = des(1,i) + dedx
                     des(2,i) = des(2,i) + dedy
                     des(3,i) = des(3,i) + dedz
                     des(1,k) = des(1,k) - dedx
                     des(2,k) = des(2,k) - dedy
                     des(3,k) = des(3,k) - dedz
c
c     increment the internal virial tensor components
c
                     vxx = xr * dedx
                     vyx = yr * dedx
                     vzx = zr * dedx
                     vyy = yr * dedy
                     vzy = zr * dedy
                     vzz = zr * dedz
                     vir(1,1) = vir(1,1) + vxx
                     vir(2,1) = vir(2,1) + vyx
                     vir(3,1) = vir(3,1) + vzx
                     vir(1,2) = vir(1,2) + vyx
                     vir(2,2) = vir(2,2) + vyy
                     vir(3,2) = vir(3,2) + vzy
                     vir(1,3) = vir(1,3) + vzx
                     vir(2,3) = vir(2,3) + vzy
                     vir(3,3) = vir(3,3) + vzz
                  end if
               end if
            end do
         end do
c
c     get Born radius chain rule components for Grycuk's HCT method
c
      else if (borntyp .eq. 'GRYCUK') then
         third = 1.0d0 / 3.0d0
         pi43 = 4.0d0 * third * pi
         factor = -(pi**third) * 6.0d0**(2.0d0*third) / 9.0d0
         do i = 1, n
            ri = rsolv(i)
            if (ri .gt. 0.0d0) then
               xi = x(i)
               yi = y(i)
               zi = z(i)
               term = pi43 / rborn(i)**3.0d0
               term = factor / term**(4.0d0*third)
               do k = 1, n
                  rk = rsolv(k)
                  if (k.ne.i .and. rk.gt.0.0d0) then
                     xr = x(k) - xi
                     yr = y(k) - yi
                     zr = z(k) - zi
                     sk = rk * shct(k)
                     sk2 = sk * sk
                     r2 = xr**2 + yr**2 + zr**2
                     r = sqrt(r2)
                     de = 0.0d0
                     if (ri+r .lt. sk) then
                        uik = sk - r
                        de = -4.0d0 * pi / uik**4
                     end if
                     if (ri+r .lt. sk) then
                        lik = sk - r
                        de = de + 0.25d0*pi*(sk2-4.0d0*sk*r+17.0d0*r2)
     &                               / (r2*lik**4)
                     else if (r .lt. ri+sk) then
                        lik = ri
                        de = de + 0.25d0*pi*(2.0d0*ri*ri-sk2-r2)
     &                               / (r2*lik**4)
                     else
                        lik = r - sk
                        de = de + 0.25d0*pi*(sk2-4.0d0*sk*r+r2)
     &                               / (r2*lik**4)
                     end if
                     uik = r + sk
                     de = de - 0.25d0*pi*(sk2+4.0d0*sk*r+r2)
     &                            / (r2*uik**4)
                     dbr = term * de/r
                     dborn = drb(i)
                     if (use_gk)  dborn = dborn + drbp(i)
                     de = dbr * dborn
c
c     increment the overall permanent solvation derivatives
c
                     dedx = de * xr
                     dedy = de * yr
                     dedz = de * zr
                     des(1,i) = des(1,i) + dedx
                     des(2,i) = des(2,i) + dedy
                     des(3,i) = des(3,i) + dedz
                     des(1,k) = des(1,k) - dedx
                     des(2,k) = des(2,k) - dedy
                     des(3,k) = des(3,k) - dedz
c
c     increment the internal virial tensor components
c
                     vxx = xr * dedx
                     vyx = yr * dedx
                     vzx = zr * dedx
                     vyy = yr * dedy
                     vzy = zr * dedy
                     vzz = zr * dedz
                     vir(1,1) = vir(1,1) + vxx
                     vir(2,1) = vir(2,1) + vyx
                     vir(3,1) = vir(3,1) + vzx
                     vir(1,2) = vir(1,2) + vyx
                     vir(2,2) = vir(2,2) + vyy
                     vir(3,2) = vir(3,2) + vzy
                     vir(1,3) = vir(1,3) + vzx
                     vir(2,3) = vir(2,3) + vzy
                     vir(3,3) = vir(3,3) + vzz
                  end if
               end do
            end if
         end do
c
c     get Born radius chain rule components for the ACE method
c
      else if (borntyp .eq. 'ACE') then
         do i = 1, n
            xi = x(i)
            yi = y(i)
            zi = z(i)
            it = class(i)
            vi = vsolv(i)
            rbi = rborn(i)
            rbi2 = rbi * rbi
            do k = 1, n
               if (k .ne. i) then
                  xr = xi - x(k)
                  yr = yi - y(k)
                  zr = zi - z(k)
                  kt = class(k)
                  vk = vsolv(k)
                  s2ik = 1.0d0 / s2ace(it,kt)
                  ws2 = wace(it,kt) * s2ik
                  uik4 = uace(it,kt)**4
                  r2 = xr**2 + yr**2 + zr**2
                  r = sqrt(r2)
                  r3 = r2 * r
                  r4 = r2 * r2
                  r6 = r2 * r4
                  rusum = r4 + uik4
                  ratio = r3 / rusum
                  expterm = exp(-r2*s2ik)
                  de1 = -4.0d0 * r * ws2 * expterm
                  de2 = 3.0d0*r2/rusum - 4.0d0*r6/rusum**2
                  dborn = drb(i)
                  if (use_gk)  dborn = dborn + drbp(i)
                  de = dborn * rbi2 * (de1+vk*ratio**3*de2/pi) / r
c
c     increment the overall implicit solvation derivatives
c
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
                  des(1,i) = des(1,i) + dedx
                  des(2,i) = des(2,i) + dedy
                  des(3,i) = des(3,i) + dedz
                  des(1,k) = des(1,k) - dedx
                  des(2,k) = des(2,k) - dedy
                  des(3,k) = des(3,k) - dedz
c
c     increment the internal virial tensor components
c
                  vxx = xr * dedx
                  vyx = yr * dedx
                  vzx = zr * dedx
                  vyy = yr * dedy
                  vzy = zr * dedy
                  vzz = zr * dedz
                  vir(1,1) = vir(1,1) + vxx
                  vir(2,1) = vir(2,1) + vyx
                  vir(3,1) = vir(3,1) + vzx
                  vir(1,2) = vir(1,2) + vyx
                  vir(2,2) = vir(2,2) + vyy
                  vir(3,2) = vir(3,2) + vzy
                  vir(1,3) = vir(1,3) + vzx
                  vir(2,3) = vir(2,3) + vzy
                  vir(3,3) = vir(3,3) + vzz
               end if
            end do
         end do
      end if
c
c     perform deallocation of some local arrays
c
      if (borntyp .eq. 'STILL')  deallocate (skip)
      deallocate (roff)
      return
      end
