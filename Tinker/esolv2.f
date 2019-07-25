c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine esolv2  --  atom-by-atom solvation Hessian  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "esolv2" calculates second derivatives of the implicit
c     solvation energy for surface area, generalized Born,
c     generalized Kirkwood and Poisson-Boltzmann solvation models
c
c
      subroutine esolv2 (i)
      use sizes
      use potent
      use solute
      use warp
      implicit none
      integer i
      real*8 probe
c     real*8, allocatable :: aes(:)
c     real*8, allocatable :: des(:,:)
c
c
c     set a value for the solvent molecule probe radius
c
      probe = 1.4d0
c
c     perform dynamic allocation of some local arrays
c
c     allocate (aes(n))
c     allocate (des(3,n))
c
c     compute the surface area-based solvation energy term
c
c     call surface1 (es,aes,des,rsolv,asolv,probe)
c
c     perform deallocation of some local arrays
c
c     deallocate (aes)
c     deallocate (des)
c
c     setup for generalized Kirkwood solvation only calculation
c
      if (solvtyp(1:2) .eq. 'GK') then
         if (.not.use_mpole .and. .not.use_polar) then
            call chkpole
            call rotpole
            call induce
         end if
      end if
c
c     get the electrostatic Hessian for GB/SA solvation
c
      if (use_born .and. solvtyp(1:2).ne.'GK') then
         if (use_smooth) then
            call egb2b (i)
         else
            call egb2a (i)
         end if
c
c     get full finite difference Hessian for other models
c
      else
         call esolv2a (i)
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine esolv2a  --  implicit solvation Hessian matrix  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "esolv2a" calculates second derivatives of the implicit solvation
c     potential energy by finite differences
c
c
      subroutine esolv2a (i)
      use sizes
      use atoms
      use charge
      use deriv
      use hessn
      use mpole
      use potent
      use solute
      implicit none
      integer i,j,k
      integer nlist
      integer, allocatable :: list(:)
      real*8 eps,old
      real*8, allocatable :: d0(:,:)
      logical prior
      logical biglist
      logical reborn
      logical reinduce
      logical twosided
c
c
c     set the default stepsize and flag for induced dipoles
c
      eps = 1.0d-5
      biglist = .false.
      reborn = .false.
      reinduce = .false.
      twosided = .false.
      if (n .le. 300) then
         biglist = .true.
         if (use_born)  reborn = .true.
         if (use_mpole .or. use_polar)  reinduce = .true.
         reborn = .true.
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (list(n))
      allocate (d0(3,n))
c
c     perform dynamic allocation of some global arrays
c
      prior = .false.
      if (allocated(des)) then
         prior = .true.
         if (size(des) .lt. 3*n) then
            deallocate (des)
         end if
      end if
      if (.not. allocated(des))  allocate (des(3,n))
c
c     optionally restrict calculation to current atom and any
c     auxiliaries; results in a faster but approximate Hessian
c
      nlist = 0
      if (biglist) then
         nlist = n
         do k = 1, n
            list(k) = k
         end do
      else
         if (use_born .and. solvtyp(1:2).ne.'GK') then
            do k = 1, nion
               if (iion(k) .eq. i) then
                  nlist = nlist + 1
                  list(nlist) = k
               end if
            end do
         else if (solvtyp(1:2) .eq. 'GK') then
            do k = 1, npole
               if (ipole(k).eq.i .or. zaxis(k).eq.i .or.
     &             xaxis(k).eq.i .or. yaxis(k).eq.i) then
                  nlist = nlist + 1
                  list(nlist) = k
               end if
            end do
         else
            nlist = 1
            list(1) = i
         end if
      end if
c
c     get solvation first derivatives for the base structure
c
      if (.not. twosided) then
         call esolv2b (nlist,list,reborn,reinduce)
         do k = 1, n
            do j = 1, 3
               d0(j,k) = des(j,k)
            end do
         end do
      end if
c
c     find numerical x-components via perturbed structures
c
      old = x(i)
      if (twosided) then
         x(i) = x(i) - 0.5d0*eps
         call esolv2b (nlist,list,reborn,reinduce)
         do k = 1, n
            do j = 1, 3
               d0(j,k) = des(j,k)
            end do
         end do
      end if
      x(i) = x(i) + eps
      call esolv2b (nlist,list,reborn,reinduce)
      x(i) = old
      do k = 1, n
         do j = 1, 3
            hessx(j,k) = hessx(j,k) + (des(j,k)-d0(j,k))/eps
         end do
      end do
c
c     find numerical y-components via perturbed structures
c
      old = y(i)
      if (twosided) then
         y(i) = y(i) - 0.5d0*eps
         call esolv2b (nlist,list,reborn,reinduce)
         do k = 1, n
            do j = 1, 3
               d0(j,k) = des(j,k)
            end do
         end do
      end if
      y(i) = y(i) + eps
      call esolv2b (nlist,list,reborn,reinduce)
      y(i) = old
      do k = 1, n
         do j = 1, 3
            hessy(j,k) = hessy(j,k) + (des(j,k)-d0(j,k))/eps
         end do
      end do
c
c     find numerical z-components via perturbed structures
c
      old = z(i)
      if (twosided) then
         z(i) = z(i) - 0.5d0*eps
         call esolv2b (nlist,list,reborn,reinduce)
         do k = 1, n
            do j = 1, 3
               d0(j,k) = des(j,k)
            end do
         end do
      end if
      z(i) = z(i) + eps
      call esolv2b (nlist,list,reborn,reinduce)
      z(i) = old
      do k = 1, n
         do j = 1, 3
            hessz(j,k) = hessz(j,k) + (des(j,k)-d0(j,k))/eps
         end do
      end do
c
c     perform deallocation of some global arrays
c
      if (.not. prior) then
         deallocate (des)
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (list)
      deallocate (d0)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine esolv2b  --  finite diffs implicit solvation  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "esolv2b" finds implicit solvation gradients needed for
c     calculation of the Hessian matrix by finite differences
c
c
      subroutine esolv2b (nlist,list,reborn,reinduce)
      implicit none
      integer nlist
      integer list(*)
      logical reborn
      logical reinduce
c
c
c     get implicit solvation gradient for finite differences
c
      if (reborn)  call born
      if (reinduce) then
         call chkpole
         call rotpole
         call induce
      end if     
      call esolv1
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine egb2a  --  atom-by-atom GB solvation Hessian  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "egb2a" calculates second derivatives of the generalized
c     Born energy term for the GB/SA solvation models
c
c     note this version does not contain the chain rule terms
c     for derivatives of Born radii with respect to coordinates
c
c
      subroutine egb2a (i)
      use sizes
      use atoms
      use charge
      use chgpot
      use hessn
      use shunt
      use solute
      implicit none
      integer i,j,k,kk
      real*8 e,de,d2e
      real*8 fi,fik
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,r3,r4
      real*8 r5,r6,r7
      real*8 dwater,rb2,rm2
      real*8 expterm,shift
      real*8 d2edx,d2edy,d2edz
      real*8 taper,dtaper,d2taper
      real*8 trans,dtrans,d2trans
      real*8 fgb,fgb2,dfgb
      real*8 dfgb2,d2fgb
      real*8 term(3,3)
      character*6 mode
c
c
c     first see if the atom of interest carries a charge
c
      do k = 1, nion
         if (iion(k) .eq. i) then
            fi = pchg(k)
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
c     set the solvent dielectric and energy conversion factor
c
      dwater = 78.3d0
      fi = -electric * (1.0d0 - 1.0d0/dwater) * fi
c
c     set cutoff distances and switching function coefficients
c
      mode = 'CHARGE'
      call switch (mode)
c
c     calculate GB polarization energy Hessian elements
c
      do kk = 1, nion
         k = iion(kk)
         if (i .ne. k) then
            xr = xi - x(k)
            yr = yi - y(k)
            zr = zi - z(k)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
               fik = fi * pchg(kk)
c
c     compute chain rule terms for Hessian matrix elements
c
               rb2 = rborn(i) * rborn(k)
               expterm = exp(-0.25d0*r2/rb2)
               fgb2 = r2 + rb2*expterm
               fgb = sqrt(fgb2)
               dfgb = (1.0d0-0.25d0*expterm) * r / fgb
               dfgb2 = dfgb * dfgb
               d2fgb = -dfgb2/fgb + dfgb/r
     &                    + 0.125d0*(r2/rb2)*expterm/fgb
               de = -fik * dfgb / fgb2
               d2e = -fik * (d2fgb-2.0d0*dfgb2/fgb) / fgb2
c
c     use energy switching if near the cutoff distance
c
               if (r2 .gt. cut2) then
                  e = fik / fgb
                  rm2 = (0.5d0 * (off+cut))**2
                  shift = fik / sqrt(rm2 + rb2*exp(-0.25d0*rm2/rb2))
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
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine egb2b  --  GB solvation Hessian for smoothing  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "egb2b" calculates second derivatives of the generalized
c     Born energy term for the GB/SA solvation models for use with
c     potential smoothing methods
c
c     note this version does not contain the chain rule terms
c     for derivatives of Born radii with respect to coordinates
c
c
      subroutine egb2b (i)
      use sizes
      use atoms
      use charge
      use chgpot
      use hessn
      use math
      use solute
      use warp
      implicit none
      integer i,j,k,kk
      real*8 fi,fik,de,d2e
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 dwater,width
      real*8 r,r2,rb2
      real*8 fgb,fgb2
      real*8 dfgb,dfgb2,d2fgb
      real*8 d2edx,d2edy,d2edz
      real*8 sterm,expterm
      real*8 erf,erfterm
      real*8 term(3,3)
      external erf
c
c
c     first see if the atom of interest carries a charge
c
      do k = 1, nion
         if (iion(k) .eq. i) then
            fi = pchg(k)
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
c     set the solvent dielectric and energy conversion factor
c
      dwater = 78.3d0
      fi = -electric * (1.0d0 - 1.0d0/dwater) * fi
c
c     set the extent of smoothing to be performed
c
      sterm = 0.5d0 / sqrt(diffc)
c
c     calculate GB polarization energy Hessian elements
c
      do kk = 1, nion
         k = iion(kk)
         if (i .ne. k) then
            xr = xi - x(k)
            yr = yi - y(k)
            zr = zi - z(k)
            r2 = xr*xr + yr*yr + zr*zr
            r = sqrt(r2)
            fik = fi * pchg(kk)
c
c     compute chain rule terms for Hessian matrix elements
c
            rb2 = rborn(i) * rborn(k)
            expterm = exp(-0.25d0*r2/rb2)
            fgb2 = r2 + rb2*expterm
            fgb = sqrt(fgb2)
            dfgb = (1.0d0-0.25d0*expterm) * r / fgb
            dfgb2 = dfgb * dfgb
            d2fgb = -dfgb2/fgb + dfgb/r
     &                 + 0.125d0*(r2/rb2)*expterm/fgb
            de = -fik * dfgb / fgb2
            d2e = -fik * (d2fgb-2.0d0*dfgb2/fgb) / fgb2
c
c     use a smoothable GB analogous to the Coulomb solution
c
            if (deform .gt. 0.0d0) then
               width = deform + 0.15d0*rb2*exp(-0.006d0*rb2/deform)
               width = sterm / sqrt(width)
               erfterm = erf(width*fgb)
               expterm = width * exp(-(width*fgb)**2) / sqrtpi
               de = de * (erfterm-2.0d0*expterm*fgb)
               d2e = d2e*erfterm + 2.0d0*fik*expterm
     &                  * (d2fgb/fgb-2.0d0*dfgb2*(1.0d0/fgb2+width**2))
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
      return
      end
