c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine echarge2  --  atomwise charge-charge Hessian  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "echarge2" calculates second derivatives of the
c     charge-charge interaction energy for a single atom
c
c
      subroutine echarge2 (i)
      use sizes
      use limits
      use warp
      implicit none
      integer i
c
c
c     choose the method for summing over pairwise interactions
c
      if (use_smooth) then
         call echarge2c (i)
      else if (use_ewald) then
         call echarge2b (i)
      else
         call echarge2a (i)
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine echarge2a  --  charge Hessian via double loop  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "echarge2a" calculates second derivatives of the charge-charge
c     interaction energy for a single atom using a pairwise double loop
c
c
      subroutine echarge2a (i)
      use sizes
      use atoms
      use bound
      use cell
      use charge
      use chgpot
      use couple
      use group
      use hessn
      use shunt
      implicit none
      integer i,j,k,kk
      integer in,kn,jcell
      real*8 e,de,d2e
      real*8 fi,fik,fgrp
      real*8 d2edx,d2edy,d2edz
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 shift,taper,trans
      real*8 dtaper,dtrans
      real*8 d2taper,d2trans
      real*8 r,rb,rb2
      real*8 r2,r3,r4
      real*8 r5,r6,r7
      real*8 term(3,3)
      real*8, allocatable :: cscale(:)
      logical proceed
      character*6 mode
c
c
c     first see if the atom of interest carries a charge
c
      do k = 1, nion
         if (iion(k) .eq. i) then
            fi = electric * pchg(k) / dielec
            in = jion(k)
            goto 10
         end if
      end do
      return
   10 continue
c
c     store the coordinates of the atom of interest
c
      xi = x(i)
      yi = y(i)
      zi = z(i)
c
c     perform dynamic allocation of some local arrays
c
      allocate (cscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do j = 1, nion
         cscale(iion(j)) = 1.0d0
      end do
      do j = 1, n12(in)
         cscale(i12(j,in)) = c2scale
      end do
      do j = 1, n13(in)
         cscale(i13(j,in)) = c3scale
      end do
      do j = 1, n14(in)
         cscale(i14(j,in)) = c4scale
      end do
      do j = 1, n15(in)
         cscale(i15(j,in)) = c5scale
      end do
c
c     set cutoff distances and switching function coefficients
c
      mode = 'CHARGE'
      call switch (mode)
c
c     calculate the charge interaction energy Hessian elements
c
      do kk = 1, nion
         k = iion(kk)
         kn = jion(kk)
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
         if (proceed)  proceed = (kn .ne. i)
c
c     compute the energy contribution for this interaction
c
         if (proceed) then
            xr = xi - x(k)
            yr = yi - y(k)
            zr = zi - z(k)
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
               rb = r + ebuffer
               rb2 = rb * rb
               fik = fi * pchg(kk) * cscale(kn)
c
c     compute chain rule terms for Hessian matrix elements
c
               de = -fik / rb2
               d2e = -2.0d0 * de/rb
c
c     use shifted energy switching if near the cutoff distance
c
               if (r2 .gt. cut2) then
                  e = fik / r
                  shift = fik / (0.5d0*(off+cut))
                  e = e - shift
                  r3 = r2 * r
                  r4 = r2 * r2
                  r5 = r2 * r3
                  r6 = r3 * r3
                  r7 = r3 * r4
                  taper = c5*r5 + c4*r4 + c3*r3 + c2*r2 + c1*r + c0
                  dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                        + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                  d2taper = 20.0d0*c5*r3 + 12.0d0*c4*r2
     &                         + 6.0d0*c3*r + 2.0d0*c2
                  trans = fik * (f7*r7 + f6*r6 + f5*r5 + f4*r4
     &                            + f3*r3 + f2*r2 + f1*r + f0)
                  dtrans = fik * (7.0d0*f7*r6 + 6.0d0*f6*r5
     &                            + 5.0d0*f5*r4 + 4.0d0*f4*r3
     &                            + 3.0d0*f3*r2 + 2.0d0*f2*r + f1)
                  d2trans = fik * (42.0d0*f7*r5 + 30.0d0*f6*r4
     &                             + 20.0d0*f5*r3 + 12.0d0*f4*r2
     &                             + 6.0d0*f3*r + 2.0d0*f2)
                  d2e = e*d2taper + 2.0d0*de*dtaper
     &                     + d2e*taper + d2trans
                  de = e*dtaper + de*taper + dtrans
               end if
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  de = de * fgrp
                  d2e = d2e * fgrp
               end if
c
c     form the individual Hessian element components
c
               de = de / r
               d2e = (d2e-de) / r2
               d2edx = d2e * xr
               d2edy = d2e * yr
               d2edz = d2e * zr
               term(1,1) = d2edx*xr + de
               term(1,2) = d2edx*yr
               term(1,3) = d2edx*zr
               term(2,1) = term(1,2)
               term(2,2) = d2edy*yr + de
               term(2,3) = d2edy*zr
               term(3,1) = term(1,3)
               term(3,2) = term(2,3)
               term(3,3) = d2edz*zr + de
c
c     increment diagonal and non-diagonal Hessian elements
c
               do j = 1, 3
                  hessx(j,i) = hessx(j,i) + term(1,j)
                  hessy(j,i) = hessy(j,i) + term(2,j)
                  hessz(j,i) = hessz(j,i) + term(3,j)
                  hessx(j,k) = hessx(j,k) - term(1,j)
                  hessy(j,k) = hessy(j,k) - term(2,j)
                  hessz(j,k) = hessz(j,k) - term(3,j)
               end do
            end if
         end if
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do kk = 1, nion
         k = iion(kk)
         kn = jion(kk)
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
c
c     compute the energy contribution for this interaction
c
         if (proceed) then
            do jcell = 1, ncell
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
               call imager (xr,yr,zr,jcell)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  rb = r + ebuffer
                  rb2 = rb * rb
                  fik = fi * pchg(kk)
                  if (use_polymer) then
                     if (r2 .le. polycut2)  fik = fik * cscale(kn)
                  end if
c
c     compute chain rule terms for Hessian matrix elements
c
                  de = -fik / rb2
                  d2e = -2.0d0 * de/rb
c
c     use shifted energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     e = fik / r
                     shift = fik / (0.5d0*(off+cut))
                     e = e - shift
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     r6 = r3 * r3
                     r7 = r3 * r4
                     taper = c5*r5 + c4*r4 + c3*r3 + c2*r2 + c1*r + c0
                     dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                           + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                     d2taper = 20.0d0*c5*r3 + 12.0d0*c4*r2
     &                            + 6.0d0*c3*r + 2.0d0*c2
                     trans = fik * (f7*r7 + f6*r6 + f5*r5 + f4*r4
     &                               + f3*r3 + f2*r2 + f1*r + f0)
                     dtrans = fik * (7.0d0*f7*r6 + 6.0d0*f6*r5
     &                               + 5.0d0*f5*r4 + 4.0d0*f4*r3
     &                               + 3.0d0*f3*r2 + 2.0d0*f2*r + f1)
                     d2trans = fik * (42.0d0*f7*r5 + 30.0d0*f6*r4
     &                                + 20.0d0*f5*r3 + 12.0d0*f4*r2
     &                                + 6.0d0*f3*r + 2.0d0*f2)
                     d2e = e*d2taper + 2.0d0*de*dtaper
     &                        + d2e*taper + d2trans
                     de = e*dtaper + de*taper + dtrans
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     de = de * fgrp
                     d2e = d2e * fgrp
                  end if
c
c     form the individual Hessian element components
c
                  de = de / r
                  d2e = (d2e-de) / r2
                  d2edx = d2e * xr
                  d2edy = d2e * yr
                  d2edz = d2e * zr
                  term(1,1) = d2edx*xr + de
                  term(1,2) = d2edx*yr
                  term(1,3) = d2edx*zr
                  term(2,1) = term(1,2)
                  term(2,2) = d2edy*yr + de
                  term(2,3) = d2edy*zr
                  term(3,1) = term(1,3)
                  term(3,2) = term(2,3)
                  term(3,3) = d2edz*zr + de
c
c     increment diagonal and non-diagonal Hessian elements
c
                  do j = 1, 3
                     hessx(j,i) = hessx(j,i) + term(1,j)
                     hessy(j,i) = hessy(j,i) + term(2,j)
                     hessz(j,i) = hessz(j,i) + term(3,j)
                     hessx(j,k) = hessx(j,k) - term(1,j)
                     hessy(j,k) = hessy(j,k) - term(2,j)
                     hessz(j,k) = hessz(j,k) - term(3,j)
                  end do
               end if
            end do
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (cscale)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine echarge2b  --  Ewald summation charge Hessian  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "echarge2b" calculates second derivatives of the charge-charge
c     interaction energy for a single atom using a particle mesh
c     Ewald summation
c
c     note this version does not incorporate the reciprocal space
c     contribution to the Hessian calculation
c
c
      subroutine echarge2b (i)
      use sizes
      use atoms
      use bound
      use cell
      use charge
      use chgpot
      use couple
      use ewald
      use group
      use hessn
      use limits
      use math
      implicit none
      integer i,j,k,kk
      integer in,kn,jcell
      real*8 fi,fik,fgrp
      real*8 r,r2,rb,rb2
      real*8 de,d2e
      real*8 d2edx,d2edy,d2edz
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 rew,erfc,cut2
      real*8 erfterm,expterm
      real*8 scale,scaleterm
      real*8 term(3,3)
      real*8, allocatable :: cscale(:)
      logical proceed
      external erfc
c
c
c     first see if the atom of interest carries a charge
c
      do k = 1, nion
         if (iion(k) .eq. i) then
            fi = electric * pchg(k) / dielec
            in = jion(k)
            cut2 = ewaldcut * ewaldcut
            goto 10
         end if
      end do
      return
   10 continue
c
c     store the coordinates of the atom of interest
c
      xi = x(i)
      yi = y(i)
      zi = z(i)
c
c     perform dynamic allocation of some local arrays
c
      allocate (cscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do j = 1, nion
         cscale(iion(j)) = 1.0d0
      end do
      do j = 1, n12(in)
         cscale(i12(j,in)) = c2scale
      end do
      do j = 1, n13(in)
         cscale(i13(j,in)) = c3scale
      end do
      do j = 1, n14(in)
         cscale(i14(j,in)) = c4scale
      end do
      do j = 1, n15(in)
         cscale(i15(j,in)) = c5scale
      end do
c
c     calculate the real space Ewald interaction Hessian elements
c
      do kk = 1, nion
         k = iion(kk)
         kn = jion(kk)
         if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
         proceed = .true.
         if (proceed)  proceed = (kn .ne. i)
c
c     compute the energy contribution for this interaction
c
         if (proceed) then
            xr = xi - x(k)
            yr = yi - y(k)
            zr = zi - z(k)
            call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. cut2) then
               r = sqrt(r2)
               rb = r + ebuffer
               rb2 = rb * rb
               fik = fi * pchg(kk) * cscale(kn)
               rew = aewald * r
               erfterm = erfc (rew)
               expterm = exp(-rew**2)
               scale = cscale(kn)
               if (use_group)  scale = scale * fgrp
               scaleterm = scale - 1.0d0
c
c     compute chain rule terms for Hessian matrix elements
c
               de = -fik * ((erfterm+scaleterm)/rb2
     &                        + (2.0d0*aewald/sqrtpi)*expterm/r)
               d2e = -2.0d0*de/rb + 2.0d0*(fik/(rb*rb2))*scaleterm
     &                  + (4.0d0*fik*aewald**3/sqrtpi)*expterm
     &                  + 2.0d0*(fik/(rb*rb2))*scaleterm
c
c     form the individual Hessian element components
c
               de = de / r
               d2e = (d2e-de) / r2
               d2edx = d2e * xr
               d2edy = d2e * yr
               d2edz = d2e * zr
               term(1,1) = d2edx*xr + de
               term(1,2) = d2edx*yr
               term(1,3) = d2edx*zr
               term(2,1) = term(1,2)
               term(2,2) = d2edy*yr + de
               term(2,3) = d2edy*zr
               term(3,1) = term(1,3)
               term(3,2) = term(2,3)
               term(3,3) = d2edz*zr + de
c
c     increment diagonal and non-diagonal Hessian elements
c
               do j = 1, 3
                  hessx(j,i) = hessx(j,i) + term(1,j)
                  hessy(j,i) = hessy(j,i) + term(2,j)
                  hessz(j,i) = hessz(j,i) + term(3,j)
                  hessx(j,k) = hessx(j,k) - term(1,j)
                  hessy(j,k) = hessy(j,k) - term(2,j)
                  hessz(j,k) = hessz(j,k) - term(3,j)
               end do
            end if
         end if
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do kk = 1, nion
         k = iion(kk)
         kn = jion(kk)
         if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
         proceed = .true.
c
c     compute the energy contribution for this interaction
c
         if (proceed) then
            do jcell = 1, ncell
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
               call imager (xr,yr,zr,jcell)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. cut2) then
                  r = sqrt(r2)
                  rb = r + ebuffer
                  rb2 = rb * rb
                  fik = fi * pchg(kk)
                  rew = aewald * r
                  erfterm = erfc (rew)
                  expterm = exp(-rew**2)
                  scale = 1.0d0
                  if (use_group)  scale = scale * fgrp
                  if (use_polymer) then
                     if (r2 .le. polycut2) then
                        scale = scale * cscale(kn)
                     end if
                  end if
                  scaleterm = scale - 1.0d0
c
c     compute chain rule terms for Hessian matrix elements
c
                  de = -fik * ((erfterm+scaleterm)/rb2
     &                    + (2.0d0*aewald/sqrtpi)*exp(-rew**2)/r)
                  d2e = -2.0d0*de/rb + 2.0d0*(fik/(rb*rb2))*scaleterm
     &                     + (4.0d0*fik*aewald**3/sqrtpi)*expterm
     &                     + 2.0d0*(fik/(rb*rb2))*scaleterm
c
c     form the individual Hessian element components
c
                  de = de / r
                  d2e = (d2e-de) / r2
                  d2edx = d2e * xr
                  d2edy = d2e * yr
                  d2edz = d2e * zr
                  term(1,1) = d2edx*xr + de
                  term(1,2) = d2edx*yr
                  term(1,3) = d2edx*zr
                  term(2,1) = term(1,2)
                  term(2,2) = d2edy*yr + de
                  term(2,3) = d2edy*zr
                  term(3,1) = term(1,3)
                  term(3,2) = term(2,3)
                  term(3,3) = d2edz*zr + de
c
c     increment diagonal and non-diagonal Hessian elements
c
                  do j = 1, 3
                     hessx(j,i) = hessx(j,i) + term(1,j)
                     hessy(j,i) = hessy(j,i) + term(2,j)
                     hessz(j,i) = hessz(j,i) + term(3,j)
                     hessx(j,k) = hessx(j,k) - term(1,j)
                     hessy(j,k) = hessy(j,k) - term(2,j)
                     hessz(j,k) = hessz(j,k) - term(3,j)
                  end do
               end if
            end do
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (cscale)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine echarge2c  --  charge Hessian for smoothing  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "echarge2c" calculates second derivatives of the charge-charge
c     interaction energy for a single atom for use with potential
c     smoothing methods
c
c
      subroutine echarge2c (i)
      use sizes
      use atoms
      use charge
      use chgpot
      use couple
      use group
      use hessn
      use math
      use warp
      implicit none
      integer i,j,k,kk
      integer in,kn
      real*8 fi,fik,fgrp
      real*8 r,r2,rb,rb2
      real*8 de,d2e
      real*8 d2edx,d2edy,d2edz
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 erf,erfterm
      real*8 expcut,expterm
      real*8 wterm,width
      real*8 width2,width3
      real*8 term(3,3)
      real*8, allocatable :: cscale(:)
      logical proceed
      external erf
c
c
c     first see if the atom of interest carries a charge
c
      do k = 1, nion
         if (iion(k) .eq. i) then
            fi = electric * pchg(k) / dielec
            in = jion(k)
            goto 10
         end if
      end do
      return
   10 continue
c
c     store the coordinates of the atom of interest
c
      xi = x(i)
      yi = y(i)
      zi = z(i)
c
c     perform dynamic allocation of some local arrays
c
      allocate (cscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do j = 1, nion
         cscale(iion(j)) = 1.0d0
      end do
      do j = 1, n12(in)
         cscale(i12(j,in)) = c2scale
      end do
      do j = 1, n13(in)
         cscale(i13(j,in)) = c3scale
      end do
      do j = 1, n14(in)
         cscale(i14(j,in)) = c4scale
      end do
      do j = 1, n15(in)
         cscale(i15(j,in)) = c5scale
      end do
c
c     set the smallest exponential terms to be calculated
c
      expcut = -50.0d0
c
c     set the extent of smoothing to be performed
c
      width = deform * diffc
      if (use_dem) then
         if (width .gt. 0.0d0)  width = 0.5d0 / sqrt(width)
      else if (use_gda) then
         wterm = sqrt(3.0d0/(2.0d0*diffc))
      end if
      width2 = width * width
      width3 = width * width2
c
c     calculate the charge interaction energy Hessian elements
c
      do kk = 1, nion
         k = iion(kk)
         kn = jion(kk)
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
         if (proceed)  proceed = (kn .ne. i)
c
c     compute the energy contribution for this interaction
c
         if (proceed) then
            xr = xi - x(k)
            yr = yi - y(k)
            zr = zi - z(k)
            r2 = xr*xr + yr*yr + zr*zr
            r = sqrt(r2)
            rb = r + ebuffer
            rb2 = rb * rb
            fik = fi * pchg(kk) * cscale(kn)
c
c     compute chain rule terms for Hessian matrix elements
c
            de = -fik / rb2
            d2e = -2.0d0 * de/rb
c
c     transform the potential function via smoothing
c
            if (use_dem) then
               if (width .gt. 0.0d0) then
                  erfterm = erf(width*rb)
                  expterm = -rb2 * width2
                  if (expterm .gt. expcut) then
                     expterm = 2.0d0*fik*width*exp(expterm)
     &                            / (sqrtpi*rb)
                  else
                     expterm = 0.0d0
                  end if
                  de = de*erfterm + expterm
                  d2e = -2.0d0 * (de/rb + expterm*rb*width2)
               end if
            else if (use_gda) then
               width = m2(i) + m2(k)
               if (width .gt. 0.0d0) then
                  width = wterm / sqrt(width)
                  width2 = width * width
                  erfterm = erf(width*rb)
                  expterm = -rb2 * width2
                  if (expterm .gt. expcut) then
                     expterm = 2.0d0*fik*width*exp(expterm)
     &                            / (sqrtpi*rb)
                  else
                     expterm = 0.0d0
                  end if
                  de = de*erfterm + expterm
                  d2e = -2.0d0 * (de/rb + expterm*r*width2)
               end if
            else if (use_tophat) then
               if (width .gt. rb) then
                  d2e = -fik / width3
                  de = d2e * rb
               end if
            else if (use_stophat) then
               wterm = rb + width
               de = -fik / (wterm*wterm)
               d2e = -2.0d0 * de / wterm
            end if
c
c     scale the interaction based on its group membership
c
            if (use_group) then
               de = de * fgrp
               d2e = d2e * fgrp
            end if
c
c     form the individual Hessian element components
c
            de = de / r
            d2e = (d2e-de) / r2
            d2edx = d2e * xr
            d2edy = d2e * yr
            d2edz = d2e * zr
            term(1,1) = d2edx*xr + de
            term(1,2) = d2edx*yr
            term(1,3) = d2edx*zr
            term(2,1) = term(1,2)
            term(2,2) = d2edy*yr + de
            term(2,3) = d2edy*zr
            term(3,1) = term(1,3)
            term(3,2) = term(2,3)
            term(3,3) = d2edz*zr + de
c
c     increment diagonal and non-diagonal Hessian elements
c
            do j = 1, 3
               hessx(j,i) = hessx(j,i) + term(1,j)
               hessy(j,i) = hessy(j,i) + term(2,j)
               hessz(j,i) = hessz(j,i) + term(3,j)
               hessx(j,k) = hessx(j,k) - term(1,j)
               hessy(j,k) = hessy(j,k) - term(2,j)
               hessz(j,k) = hessz(j,k) - term(3,j)
            end do
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (cscale)
      return
      end
