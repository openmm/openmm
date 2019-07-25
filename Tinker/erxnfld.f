c
c
c     ############################################################
c     ##  COPYRIGHT (C) 1996 by Yong Kong & Jay William Ponder  ##
c     ##                  All Rights Reserved                   ##
c     ############################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine erxnfld  --  reaction field potential energy  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "erxnfld" calculates the macroscopic reaction field energy
c     arising from a set of atomic multipoles
c
c     literature reference:
c
c     Y. Kong and J. W. Ponder, "Reaction Field Methods for Off-Center
c     Multipoles", Journal of Chemical Physics, 107, 481-492 (1997)
c
c
      subroutine erxnfld
      use sizes
      use atoms
      use chgpot
      use energi
      use mpole
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer ii,kk
      integer ix,iy,iz
      integer kx,ky,kz
      real*8 eik,r2
      real*8 xr,yr,zr
      real*8 rpi(13)
      real*8 rpk(13)
      logical usei,usek
      character*6 mode
c
c
c     zero out the macroscopic reaction field energy
c
      er = 0.0d0
c
c     set the switching function coefficients
c
      mode = 'MPOLE'
      call switch (mode)
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     compute the indices used in reaction field calculations
c
      call ijkpts
c
c     calculate the reaction field interaction energy term
c
      do ii = 1, npole
         i = ipole(ii)
         iz = zaxis(ii)
         ix = xaxis(ii)
         iy = yaxis(ii)
         usei = (use(i) .or. use(iz) .or. use(ix) .or. use(iy))
         do j = 1, polsiz(ii)
            rpi(j) = rpole(j,ii)
         end do
         do kk = ii, npole
            k = ipole(kk)
            kz = zaxis(kk)
            kx = xaxis(kk)
            ky = yaxis(kk)
            usek = (use(k) .or. use(kz) .or. use(kx) .or. use(ky))
            if (usei .or. usek) then
               xr = x(k) - x(i)
               yr = y(k) - y(i)
               zr = z(k) - z(i)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  do j = 1, polsiz(kk)
                     rpk(j) = rpole(j,kk)
                  end do
                  call erfik (ii,kk,i,k,rpi,rpk,eik)
                  er = er + eik
               end if
            end if
         end do
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine erfik  --  reaction field energy of site pair   ##
c     ##                                                             ##
c     #################################################################
c
c
c     "erfik" compute the reaction field energy due to a single pair
c     of atomic multipoles
c
c
      subroutine erfik (ii,kk,i,k,rpi,rpk,eik)
      use sizes
      use atoms
      use chgpot
      use mpole
      use rxnfld
      use rxnpot
      implicit none
      integer i,k,ii,kk
      integer isiz,ksiz
      integer m,n1,n2,nn
      integer fii,fi,fj
      integer p_s1,p_s2,p_e1,p_e2
      integer ind1_x(13),ind2_x(13)
      integer ind1_y(13),ind2_y(13)
      integer ind1_z(13),ind2_z(13)
      real*8 eik,term,ratio,factor
      real*8 xi,yi,zi,xk,yk,zk
      real*8 size,size2,d,d2,ri2,rk2
      real*8 xi2,yi2,zi2,xk2,yk2,zk2
      real*8 rpi(13),rpk(13)
      real*8 rklij(13,13),d1d2
      real*8 m2t2(13)
c
c
c     get numbers of the atoms
c
      isiz = polsiz(ii)
      ksiz = polsiz(kk)
      xi = x(i)
      yi = y(i)
      zi = z(i)
      xk = x(k)
      yk = y(k)
      zk = z(k)
      d = xi*xk + yi*yk + zi*zk
      ri2 = xi*xi + yi*yi + zi*zi
      rk2 = xk*xk + yk*yk + zk*zk
c
c     set highest order of multipoles at each site (M=0, D=1, Q=2)
c
      eik = 0.0d0
      n1 = 2
      n2 = 2
      nn = rfterms
      ratio = rfbulkd / dielec
      factor = electric * (1.0d0-ratio)
      if (i .eq. k)  factor = 0.5d0 * factor
      size = 1.0d0 / rfsize
      size2 = size * size
c
c     get the values of the indices
c
      m = (3**(n1+1)-1)/2
      call rfindex (n1,m,ind1_x,ind1_y,ind1_z,p_s1,p_e1)
      m = (3**(n2+1)-1)/2
      call rfindex (n2,m,ind2_x,ind2_y,ind2_z,p_s2,p_e2)
c
c     initialize the stored matrix element arrays
c
      do fi = 1, p_e1
         do fj = 1, p_e2
            b1(fi,fj) = 0.0d0
            b2(fi,fj) = 0.0d0
         end do
      end do
c
c     explicit formula for the 0th summation term
c
      if (nn .ge. 0) then
         eik = size * rpi(1) * rpk(1) / ratio
         size = size * size2
      end if
c
c     explicit formula for the 1st summation term
c
      if (nn .ge. 1) then
         b2(1,1) = d
         b2(1,2) = xi
         b2(1,3) = yi
         b2(1,4) = zi
         b2(2,1) = xk
         b2(3,1) = yk
         b2(4,1) = zk
         b2(2,2) = 1.0d0
         b2(3,3) = 1.0d0
         b2(4,4) = 1.0d0
         do fi = 1, 4
            m2t2(fi) = 0.0d0
            do fj = 1, 4
               m2t2(fi) = m2t2(fi) + b2(fi,fj)*rpk(fj)
            end do
         end do
         term = 0.0d0
         do fi = 1, 4
            term = term + rpi(fi)*m2t2(fi)
         end do
         term = 2.0d0 * size * term / (2.0d0*ratio+1.0d0)
         eik = eik + term
         size = size * size2
      end if
c
c     explicit formula for the 2nd summation term
c
      if (nn .ge. 2) then
         b2(1,1) = (3.0d0*d*d-ri2*rk2) * 0.5d0
         b2(1,2) = 3.0d0*xi*d - xk*ri2
         b2(1,3) = 3.0d0*yi*d - yk*ri2
         b2(1,4) = 3.0d0*zi*d - zk*ri2
         b2(1,5) = 3.0d0*xi*xi - ri2
         b2(1,6) = 3.0d0*xi*yi
         b2(1,7) = 3.0d0*xi*zi
         b2(1,8) = b2(1,6)
         b2(1,9) = 3.0d0*yi*yi - ri2
         b2(1,10) = 3.0d0*yi*zi
         b2(1,11) = b2(1,7)
         b2(1,12) = b2(1,10)
         b2(1,13) = 3.0d0*zi*zi - ri2
         b2(2,1) = 3.0d0*xk*d - xi*rk2
         b2(2,2) = 3.0d0*d + xi*xk
         b2(2,3) = 3.0d0*xk*yi - 2.0d0*xi*yk
         b2(2,4) = 3.0d0*zi*xk - 2.0d0*xi*zk
         b2(2,5) = 4.0d0*xi
         b2(2,6) = 3.0d0*yi
         b2(2,7) = 3.0d0*zi
         b2(2,8) = b2(2,6)
         b2(2,9) = -2.0d0*xi
         b2(2,11) = b2(2,7)
         b2(2,13) = b2(2,9)
         b2(3,1) = 3.0d0*yk*d - yi*rk2
         b2(3,2) = 3.0d0*yk*xi - 2.0d0*yi*xk
         b2(3,3) = 3.0d0*d + yi*yk
         b2(3,4) = 3.0d0*yk*zi - 2.0d0*yi*zk
         b2(3,5) = -2.0d0*yi
         b2(3,6) = 3.0d0*xi
         b2(3,8) = b2(3,6)
         b2(3,9) = 4.0d0*yi
         b2(3,10) = 3.0d0*zi
         b2(3,12) = b2(3,10)
         b2(3,13) = b2(3,5)
         b2(4,1) = 3.0d0*zk*d - zi*rk2
         b2(4,2) = 3.0d0*zk*xi - 2.0d0*zi*xk
         b2(4,3) = 3.0d0*zk*yi - 2.0d0*zi*yk
         b2(4,4) = 3.0d0*d + zi*zk
         b2(4,5) = -2.0d0*zi
         b2(4,7) = 3.0d0*xi
         b2(4,9) = b2(4,5)
         b2(4,10) = 3.0d0*yi
         b2(4,11) = b2(4,7)
         b2(4,12) = b2(4,10)
         b2(4,13) = 4.0d0*zi
         b2(5,1) = 3.0d0*xk*xk - rk2
         b2(5,2) = 4.0d0*xk
         b2(5,3) = -2.0d0*yk
         b2(5,4) = -2.0d0*zk
         b2(5,5) = 4.0d0
         b2(5,9) = -2.0d0
         b2(5,13) = -2.0d0
         b2(6,1) = 3.0d0*xk*yk
         b2(6,2) = 3.0d0*yk
         b2(6,3) = 3.0d0*xk
         b2(6,6) = 3.0d0
         b2(6,8) = 3.0d0
         b2(7,1) = 3.0d0*xk*zk
         b2(7,2) = 3.0d0*zk
         b2(7,4) = 3.0d0*xk
         b2(7,7) = 3.0d0
         b2(7,11) = 3.0d0
         b2(8,1) = b2(6,1)
         b2(8,2) = b2(6,2)
         b2(8,3) = b2(6,3)
         b2(8,6) = 3.0d0
         b2(8,8) = 3.0d0
         b2(9,1) = 3.0d0*yk*yk - rk2
         b2(9,2) = -2.0d0*xk
         b2(9,3) = 4.0d0*yk
         b2(9,4) = -2.0d0*zk
         b2(9,5) = -2.0d0
         b2(9,9) = 4.0d0
         b2(9,13) = -2.0d0
         b2(10,1) = 3.0d0*yk*zk
         b2(10,3) = 3.0d0*zk
         b2(10,4) = 3.0d0*yk
         b2(10,10) = 3.0d0
         b2(10,12) = 3.0d0
         b2(11,1) = b2(7,1)
         b2(11,2) = b2(7,2)
         b2(11,4) = b2(7,4)
         b2(11,7) = 3.0d0
         b2(11,11) = 3.0d0
         b2(12,1) = b2(10,1)
         b2(12,3) = b2(10,3)
         b2(12,4) = b2(10,4)
         b2(12,10) = 3.0d0
         b2(12,12) = 3.0d0
         b2(13,1) = 3.0d0*zk*zk - rk2
         b2(13,2) = -2.0d0*xk
         b2(13,3) = -2.0d0*yk
         b2(13,4) = 4.0d0*zk
         b2(13,5) = -2.0d0
         b2(13,9) = -2.0d0
         b2(13,13) = 4.0d0
         do fi = 1, isiz
            m2t2(fi) = 0.0d0
            do fj = 1, ksiz
               m2t2(fi) = m2t2(fi) + b2(fi,fj)*rpk(fj)
            end do
         end do
         term = 0.0d0
         do fi = 1, isiz
            term = term + rpi(fi)*m2t2(fi)
         end do
         term = 3.0d0 * size * term / (3.0d0*ratio+2.0d0)
         eik = eik + term
         size = size * size2
      end if
c
c     explicit formula for the 3rd summation term
c
      if (nn .ge. 3) then
         d2 = d*d
         xi2 = xi*xi
         yi2 = yi*yi
         zi2 = zi*zi
         xk2 = xk*xk
         yk2 = yk*yk
         zk2 = zk*zk
         b1(1,1) = d*(2.5d0*d2-1.5d0*ri2*rk2)
         b1(1,2) = 7.5d0*d2*xi-3.0d0*xk*ri2*d-1.5d0*xi*ri2*rk2
         b1(1,3) = 7.5d0*d2*yi-3.0d0*yk*ri2*d-1.5d0*yi*ri2*rk2
         b1(1,4) = 7.5d0*d2*zi-3.0d0*zk*ri2*d-1.5d0*zi*ri2*rk2
         b1(1,5) = 15.0d0*d*xi2-3.0d0*ri2*(d+2.0d0*xi*xk)
         b1(1,6) = 15.0d0*xi*yi*d - 3.0d0*ri2*(xi*yk+xk*yi)
         b1(1,7) = 15.0d0*xi*zi*d - 3.0d0*ri2*(xi*zk+xk*zi)
         b1(1,8) = b1(1,6)
         b1(1,9) = 15.0d0*d*yi2-3.0d0*ri2*(d+2.0d0*yi*yk)
         b1(1,10) = 15.0d0*yi*zi*d - 3.0d0*ri2*(yi*zk+yk*zi)
         b1(1,11) = b1(1,7)
         b1(1,12) = b1(1,10)
         b1(1,13) = 15.0d0*d*zi2-3.0d0*ri2*(d+2.0d0*zi*zk)
         b1(2,1) = 7.5d0*d2*xk-3.0d0*xi*rk2*d-1.5d0*xk*ri2*rk2
         b1(2,2) = 7.5d0*d2+9.0d0*xi*xk*d-3.0d0*xi2*rk2-3.0d0*xk2*ri2
     &                -1.5d0*ri2*rk2
         b1(2,3) = 3.0d0*((5.0d0*xk*yi-2.0d0*xi*yk)*d
     &                -xi*yi*rk2-xk*yk*ri2)
         b1(2,4) = 3.0d0*((5.0d0*xk*zi-2.0d0*xi*zk)*d
     &                -xi*zi*rk2-xk*zk*ri2)
         b1(2,5) = 24.0d0*xi*yi*yk + 24.0d0*xi*zi*zk + 18.0d0*xi2*xk
     &                - 9.0d0*xk*yi2  - 9.0d0*xk*zi2
         b1(2,6) = (8.0d0*yi*xk*xi - 3.0d0*xi2*yk + 4.0d0*yi2*yk
     &                - yk*zi2  + 5.0d0*yi*zi*zk)*3.0d0
         b1(2,7) = 15.0d0*zi*yi*yk + 12.0d0*zi2*zk - 9.0d0*xi2*zk
     &                - 3.0d0*zk*yi2  + 24.0d0*zi*xk*xi
         b1(2,8) = b1(2,6)
         b1(2,9) = - 9.0d0*xi2*xk + 12.0d0*xk*yi2  - 3.0d0*xk*zi2
     &                - 18.0d0*xi*yi*yk - 6.0d0*xi*zi*zk
         b1(2,10) = 15.0d0*zi*xk*yi - 6.0d0*zi*xi*yk - 6.0d0*yi*xi*zk
         b1(2,11) = b1(2,7)
         b1(2,12) = b1(2,10)
         b1(2,13) = - 6.0d0*xi*yi*yk - 9.0d0*xi2*xk - 3.0d0*xk*yi2
     &                 + 12.0d0*xk*zi2  - 18.0d0*xi*zi*zk
         b1(3,1) = 7.5d0*d2*yk-3.0d0*yi*rk2*d-1.5d0*yk*ri2*rk2
         b1(3,2) = 3.0d0*((5.0d0*xi*yk-2.0d0*xk*yi)*d
     &                -xi*yi*rk2-xk*yk*ri2)
         b1(3,3) = 7.5d0*d2+9.0d0*yi*yk*d-3.0d0*yi2*rk2-3.0d0*yk2*ri2
     &                -1.5d0*ri2*rk2
         b1(3,4) = 3.0d0*((5.0d0*yk*zi-2.0d0*yi*zk)*d
     &                -yi*zi*rk2-yk*zk*ri2)
         b1(3,5) = - 9.0d0*yi2*yk - 6.0d0*yi*zi*zk - 18.0d0*yi*xk*xi
     &                + 12.0d0*xi2*yk - 3.0d0*yk*zi2
         b1(3,6) = 12.0d0*xi2*xk + 15.0d0*xi*zi*zk - 9.0d0*xk*yi2
     &                - 3.0d0*xk*zi2  + 24.0d0*xi*yi*yk
         b1(3,7) = 15.0d0*zi*xi*yk - 6.0d0*yi*xi*zk - 6.0d0*zi*xk*yi
         b1(3,8) = b1(3,6)
         b1(3,9) = - 9.0d0*xi2*yk + 18.0d0*yi2*yk - 9.0d0*yk*zi2
     &                + 24.0d0*yi*xk*xi + 24.0d0*yi*zi*zk
         b1(3,10) = 24.0d0*zi*yi*yk - 3.0d0*xi2*zk - 9.0d0*zk*yi2
     &                + 12.0d0*zi2*zk + 15.0d0*zi*xk*xi
         b1(3,11) = b1(3,7)
         b1(3,12) = b1(3,10)
         b1(3,13) = - 3.0d0*xi2*yk - 9.0d0*yi2*yk + 12.0d0*yk*zi2
     &                 - 18.0d0*yi*zi*zk - 6.0d0*yi*xk*xi
         b1(4,1) = 7.5d0*d2*zk-3.0d0*zi*rk2*d-1.5d0*zk*ri2*rk2
         b1(4,2) = 3.0d0*((5.0d0*xi*zk-2.0d0*xk*zi)*d
     &                -xi*zi*rk2-xk*zk*ri2)
         b1(4,3) = 3.0d0*((5.0d0*yi*zk-2.0d0*yk*zi)*d
     &                -yi*zi*rk2-yk*zk*ri2)
         b1(4,4) = 7.5d0*d2+9.0d0*zi*zk*d-3.0d0*zi2*rk2-3.0d0*zk2*ri2
     &                -1.5d0*ri2*rk2
         b1(4,5) = 12.0d0*xi2*zk - 3.0d0*zk*yi2 - 9.0d0*zi2*zk
     &                - 18.0d0*zi*xk*xi - 6.0d0*zi*yi*yk
         b1(4,6) = 15.0d0*yi*xi*zk - 6.0d0*zi*xi*yk - 6.0d0*zi*xk*yi
         b1(4,7) = 24.0d0*xi*zi*zk + 12.0d0*xi2*xk - 3.0d0*xk*yi2
     &                - 9.0d0*xk*zi2  + 15.0d0*xi*yi*yk
         b1(4,8) = b1(4,6)
         b1(4,9) = - 6.0d0*zi*xk*xi - 9.0d0*zi2*zk - 3.0d0*xi2*zk
     &                + 12.0d0*zk*yi2  - 18.0d0*zi*yi*yk
         b1(4,10) = 15.0d0*yi*xk*xi + 12.0d0*yi2*yk - 9.0d0*yk*zi2
     &                + 24.0d0*yi*zi*zk - 3.0d0*xi2*yk
         b1(4,11) = b1(4,7)
         b1(4,12) = b1(4,10)
         b1(4,13) = 24.0d0*zi*xk*xi + 18.0d0*zi2*zk - 9.0d0*xi2*zk
     &                - 9.0d0*zk*yi2  + 24.0d0*zi*yi*yk
         b1(5,1) = 15.0d0*d*xk2-3.0d0*rk2*(d+2.0d0*xi*xk)
         b1(5,2) = 18.0d0*xi*xk2  + 24.0d0*xk*yi*yk + 24.0d0*xk*zi*zk
     &                - 9.0d0*xi*yk2  - 9.0d0*xi*zk2
         b1(5,3) = 12.0d0*yi*xk2 - 9.0d0*yk2*yi - 3.0d0*yi*zk2
     &                - 18.0d0*xk*xi*yk - 6.0d0*yk*zi*zk
         b1(5,4) = - 9.0d0*zk2*zi - 6.0d0*zk*yi*yk - 18.0d0*xk*xi*zk
     &                + 12.0d0*zi*xk2  - 3.0d0*zi*yk2
         b1(5,5) = 24.0d0*zi*zk + 24.0d0*yi*yk + 36.0d0*xi*xk
         b1(5,6) = -18.0d0*xi*yk + 24.0d0*yi*xk
         b1(5,7) = -18.0d0*xi*zk + 24.0d0*zi*xk
         b1(5,8) = b1(5,6)
         b1(5,9) = -6.0d0*zi*zk - 18.0d0*yi*yk - 18.0d0*xi*xk
         b1(5,10) = -6.0d0*(yi*zk + zi*yk)
         b1(5,11) = b1(5,7)
         b1(5,12) = b1(5,10)
         b1(5,13) = -6.0d0*yi*yk - 18.0d0*xi*xk - 18.0d0*zi*zk
         b1(6,1) = 15.0d0*xk*yk*d - 3.0d0*rk2*(xi*yk+xk*yi)
         b1(6,2) = -9.0d0*yi*xk2 + 12.0d0*yk2*yi - 3.0d0*yi*zk2
     &                + 24.0d0*xk*xi*yk + 15.0d0*yk*zi*zk
         b1(6,3) = 12.0d0*xi*xk2 + 15.0d0*xk*zi*zk - 9.0d0*xi*yk2
     &                - 3.0d0*xi*zk2  + 24.0d0*xk*yi*yk
         b1(6,4) = -6.0d0*xk*yi*zk - 6.0d0*yk*xi*zk + 15.0d0*zi*xk*yk
         b1(6,5) = -18.0d0*yi*xk + 24.0d0*xi*yk
         b1(6,6) = 24.0d0*yi*yk + 24.0d0*xi*xk + 15.0d0*zi*zk
         b1(6,7) = -6.0d0*yi*zk + 15.0d0*zi*yk
         b1(6,8) = b1(6,6)
         b1(6,9) = -18.0d0*xi*yk + 24.0d0*yi*xk
         b1(6,10) = -6.0d0*xi*zk + 15.0d0*zi*xk
         b1(6,11) = b1(6,7)
         b1(6,12) = b1(6,10)
         b1(6,13) = -6.0d0*yi*xk - 6.0d0*xi*yk
         b1(7,1) = 15.0d0*xk*zk*d - 3.0d0*rk2*(xi*zk+xk*zi)
         b1(7,2) = 15.0d0*zk*yi*yk + 12.0d0*zk2*zi - 9.0d0*zi*xk2
     &                - 3.0d0*zi*yk2  + 24.0d0*xk*xi*zk
         b1(7,3) = - 6.0d0*zi*xk*yk - 6.0d0*yk*xi*zk + 15.0d0*xk*yi*zk
         b1(7,4) = 12.0d0*xi*xk2  - 3.0d0*xi*yk2  - 9.0d0*xi*zk2
     &                + 15.0d0*xk*yi*yk + 24.0d0*xk*zi*zk
         b1(7,5) = -18.0d0*zi*xk + 24.0d0*xi*zk
         b1(7,6) = -6.0d0*zi*yk + 15.0d0*yi*zk
         b1(7,7) = 24.0d0*xi*xk + 24.0d0*zi*zk + 15.0d0*yi*yk
         b1(7,8) = b1(7,6)
         b1(7,9) = -6.0d0*zi*xk - 6.0d0*xi*zk
         b1(7,10) = -6.0d0*xi*yk + 15.0d0*yi*xk
         b1(7,11) = b1(7,7)
         b1(7,12) = b1(7,10)
         b1(7,13) = -18.0d0*xi*zk + 24.0d0*zi*xk
         b1(9,1) = 15.0d0*d*yk2-3.0d0*rk2*(d+2.0d0*yi*yk)
         b1(9,2) = -9.0d0*xi*xk2 + 12.0d0*xi*yk2 - 3.0d0*xi*zk2
     &                - 18.0d0*xk*yi*yk - 6.0d0*xk*zi*zk
         b1(9,3) = -9.0d0*yi*xk2  + 18.0d0*yk2*yi - 9.0d0*yi*zk2
     &                + 24.0d0*yk*zi*zk + 24.0d0*xk*xi*yk
         b1(9,4) = 12.0d0*zi*yk2 - 18.0d0*zk*yi*yk - 3.0d0*zi*xk2
     &                - 9.0d0*zk2*zi - 6.0d0*xk*xi*zk
         b1(9,5) = -18.0d0*xi*xk - 6.0d0*zi*zk - 18.0d0*yi*yk
         b1(9,6) = -18.0d0*yi*xk + 24.0d0*xi*yk
         b1(9,7) = -6.0d0*zi*xk - 6.0d0*xi*zk
         b1(9,8) = b1(9,6)
         b1(9,9) = 24.0d0*xi*xk + 24.0d0*zi*zk + 36.0d0*yi*yk
         b1(9,10) = -18.0d0*yi*zk + 24.0d0*zi*yk
         b1(9,11) = b1(9,7)
         b1(9,12) = b1(9,10)
         b1(9,13) = -18.0d0*yi*yk - 6.0d0*xi*xk - 18.0d0*zi*zk
         b1(10,1) = 15.0d0*yk*zk*d - 3.0d0*rk2*(yi*zk+yk*zi)
         b1(10,2) = -6.0d0*zi*xk*yk -6.0d0*xk*yi*zk + 15.0d0*yk*xi*zk
         b1(10,3) = 12.0d0* zk2*zi + 15.0d0*xk*xi*zk - 3.0d0*zi*xk2
     &                 - 9.0d0*zi*yk2  + 24.0d0*zk*yi*yk
         b1(10,4) = 15.0d0*xk*xi*yk + 12.0d0*yk2*yi - 3.0d0*yi*xk2
     &                 - 9.0d0*yi*zk2  + 24.0d0*yk*zi*zk
         b1(10,5) = -6.0d0*yi*zk - 6.0d0*zi*yk
         b1(10,6) = -6.0d0*zi*xk + 15.0d0*xi*zk
         b1(10,7) = -6.0d0*yi*xk + 15.0d0*xi*yk
         b1(10,8) = b1(10,6)
         b1(10,9) = 24.0d0*yi*zk - 18.0d0*zi*yk
         b1(10,10) = 15.0d0*xi*xk + 24.0d0*zi*zk + 24.0d0*yi*yk
         b1(10,11) = b1(10,7)
         b1(10,12) = b1(10,10)
         b1(10,13) = -18.0d0*yi*zk + 24.0d0*zi*yk
         b1(13,1) = 15.0d0*d*zk2-3.0d0*rk2*(d+2.0d0*zi*zk)
         b1(13,2) = 12.0d0*xi*zk2 - 18.0d0*xk*zi*zk - 9.0d0*xi*xk2
     &                 - 3.0d0*xi*yk2 - 6.0d0*xk*yi*yk
         b1(13,3) = 12.0d0*yi*zk2 - 3.0d0*yi*xk2 - 9.0d0*yk2*yi
     &                 - 18.0d0*yk*zi*zk - 6.0d0*xk*xi*yk
         b1(13,4) = -9.0d0*zi*xk2 - 9.0d0*zi*yk2 + 18.0d0*zk2*zi
     &                 + 24.0d0*xk*xi*zk + 24.0d0*zk*yi*yk
         b1(13,5) = -6.0d0*yi*yk - 18.0d0*zi*zk - 18.0d0*xi*xk
         b1(13,6) = -6.0d0*yi*xk - 6.0d0*xi*yk
         b1(13,7) = 24.0d0*xi*zk - 18.0d0*zi*xk
         b1(13,8) = b1(13,6)
         b1(13,9) = -18.0d0*yi*yk - 6.0d0*xi*xk - 18.0d0*zi*zk
         b1(13,10) = 24.0d0*yi*zk - 18.0d0*zi*yk
         b1(13,11) = b1(13,7)
         b1(13,12) = b1(13,10)
         b1(13,13) = 36.0d0*zi*zk + 24.0d0*xi*xk + 24.0d0*yi*yk
         do fi = 1, isiz
            b1(8,fi) = b1(6,fi)
            b1(11,fi) = b1(7,fi)
            b1(12,fi) = b1(10,fi)
         end do
         do fi = 1, isiz
            m2t2(fi) = 0.0d0
            do fj = 1, ksiz
               m2t2(fi) = m2t2(fi) + b1(fi,fj)*rpk(fj)
            end do
         end do
         term = 0.0d0
         do fi = 1, isiz
            term = term + rpi(fi)*m2t2(fi)
         end do
         term = 4.0d0 * size * term / (4.0d0*ratio+3.0d0)
         eik = eik + term
         size = size * size2
      end if
c
c     recursive formulation of 4th through nth summation terms
c
      do fii = 4, nn
         do fi = 1, p_e1
            if (fi .eq. 8) then
               do fj = 1, p_e2
                  rklij(fi,fj) = rklij(6,fj)
               end do
            else if (fi .eq. 11) then
               do fj = 1, p_e2
                  rklij(fi,fj) = rklij(7,fj)
               end do
            else if (fi .eq. 12) then
               do fj = 1, p_e2
                  rklij(fi,fj) = rklij(10,fj)
               end do
            else
               do fj = 1, p_e2
                  if (fj .eq. 8) then
                     rklij(fi,fj) = rklij(fi,6)
                  else if (fj .eq. 11) then
                     rklij(fi,fj) = rklij(fi,7)
                  else if (fj .eq. 12) then
                     rklij(fi,fj) = rklij(fi,10)
                  else
                     rklij(fi,fj) = d1d2 (fii,xi,yi,zi,xk,yk,zk,
     &                                    d,ri2,rk2,ind1_x(fi),
     &                                    ind1_y(fi),ind1_z(fi),
     &                              ind2_x(fj),ind2_y(fj),ind2_z(fj))
                  end if
               end do
            end if
         end do
c
c     update storage of the last two sets of matrix elements
c
         do fi = 1, p_e1
           do fj = 1, p_e2
              b2(fj,fi) = b1(fj,fi)
              b1(fj,fi) = rklij(fj,fi)
           end do
         end do
c
c     compute interaction energy between the two multipole sites
c
         do fi = 1, isiz
            m2t2(fi) = 0.0d0
            do fj = 1, ksiz
               m2t2(fi) = m2t2(fi) + rklij(fi,fj)*rpk(fj)
            end do
         end do
         term = 0.0d0
         do fi = 1, isiz
            term = term + rpi(fi)*m2t2(fi)
         end do
         term = term * size * dble(fii+1)
     &             / (dble(fii+1)*ratio+dble(fii))
         eik = eik + term
         size = size * size2
      end do
      eik = factor * eik
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine rfindex  --  reaction field indices for sites  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "rfindex" finds indices for each multipole site for use
c     in computing reaction field energetics
c
c
      subroutine rfindex (n,m,ind_x,ind_y,ind_z,p_s,p_e)
      implicit none
      integer i,j,k,n,m
      integer p_s,p_e
      integer ind_x(*)
      integer ind_y(*)
      integer ind_z(*)
c
c
      p_s = 1
      p_e = 1
      do i = 1, m
         ind_x(i) = 0
         ind_y(i) = 0
         ind_z(i) = 0
      end do
      k = 1
      do i = 1, n
         do j = p_s, p_e
            k = k + 1
            ind_x(k) = ind_x(j) + 1
            ind_y(k) = ind_y(j)
            ind_z(k) = ind_z(j)
         end do
         do j = p_s, p_e
            k = k + 1
            ind_x(k) = ind_x(j)
            ind_y(k) = ind_y(j) + 1
            ind_z(k) = ind_z(j)
         end do
         do j = p_s, p_e
            k = k + 1
            ind_x(k) = ind_x(j)
            ind_y(k) = ind_y(j)
            ind_z(k) = ind_z(j) + 1
         end do
         p_s = p_e + 1
         p_e = k
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ijkpts  --  indices into "b1" and "b2" arrays  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ijkpts" stores a set of indices used during calculation
c     of macroscopic reaction field energetics
c
      subroutine ijkpts
      use rxnfld
      implicit none
      integer i,j,k
c
c
c     find and store indices into the "b1" and "b2" arrays
c
      do i = 0, 5
         do j = 0, 5
            do k = 0, 5
               ijk(i,j,k) = (3**(i+j+k) + 3**(j+k) + 3**k - 1) / 2
            end do
         end do
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  function d1d2  --  recursive summation element utility  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "d1d2" is a utility function used in computation of the
c     reaction field recursive summation elements
c
c
      function d1d2 (n,x1,y1,z1,x2,y2,z2,d,r1sq,r2sq,i,j,k,s,t,u)
      use rxnfld
      implicit none
      integer n,i,j,k
      integer s,t,u
      integer is,it,iu
      integer js,jt,ju
      integer ks,kt,ku
      real*8 x1,y1,z1
      real*8 x2,y2,z2
      real*8 d1d2,d,f,g
      real*8 r1sq,r2sq
c
c
      if (n.lt.i+j+k  .or. n.lt.s+t+u) then
         d1d2 = 0.0d0
         return
      end if
      is = i*s
      it = i*t
      iu = i*u
      js = j*s
      jt = j*t
      ju = j*u
      ks = k*s
      kt = k*t
      ku = k*u
      f = d*b1(ijk(i,j,k),ijk(s,t,u))
      g = r1sq*r2sq*b2(ijk(i,j,k),ijk(s,t,u))
      if (i .ne. 0) then
         f = f + i*x2*b1(ijk(i-1,j,k),ijk(s,t,u))
         g = g + 2.0d0*i*x1*r2sq*b2(ijk(i-1,j,k),ijk(s,t,u))
         if (i .ne. 1)  g = g + i*(i-1)*r2sq*b2(ijk(i-2,j,k),ijk(s,t,u))
      end if
      if (j .ne. 0) then
         f = f + j*y2*b1(ijk(i,j-1,k),ijk(s,t,u))
         g = g + 2.0d0*j*y1*r2sq*b2(ijk(i,j-1,k),ijk(s,t,u))
         if (j .ne. 1)  g = g + j*(j-1)*r2sq*b2(ijk(i,j-2,k),ijk(s,t,u))
      end if
      if (k .ne. 0) then
         f = f + k*z2*b1(ijk(i,j,k-1),ijk(s,t,u))
         g = g + 2.0d0*k*z1*r2sq*b2(ijk(i,j,k-1),ijk(s,t,u))
         if (k .ne. 1)  g = g + k*(k-1)*r2sq*b2(ijk(i,j,k-2),ijk(s,t,u))
      end if
      if (s .ne. 0) then
         f = f + s*x1*b1(ijk(i,j,k),ijk(s-1,t,u))
         g = g + 2.0d0*s*x2*r1sq*b2(ijk(i,j,k),ijk(s-1,t,u))
         if (s .ne. 1)  g = g + s*(s-1)*r1sq*b2(ijk(i,j,k),ijk(s-2,t,u))
      end if
      if (t .ne. 0) then
         f = f + t*y1*b1(ijk(i,j,k),ijk(s,t-1,u))
         g = g + 2.0d0*t*y2*r1sq*b2(ijk(i,j,k),ijk(s,t-1,u))
         if (t .ne. 1)  g = g + t*(t-1)*r1sq*b2(ijk(i,j,k),ijk(s,t-2,u))
      end if
      if (u .ne. 0) then
         f = f + u*z1*b1(ijk(i,j,k),ijk(s,t,u-1))
         g = g + 2.0d0*u*z2*r1sq*b2(ijk(i,j,k),ijk(s,t,u-1))
         if (u .ne. 1)  g = g + u*(u-1)*r1sq*b2(ijk(i,j,k),ijk(s,t,u-2))
      end if
      if (is .ne. 0) then
         f = f + is*b1(ijk(i-1,j,k),ijk(s-1,t,u))
         g = g + 4.0d0*is*x1*x2*b2(ijk(i-1,j,k),ijk(s-1,t,u))
         if (i .ne. 1) then
            g = g + 2.0d0*(i-1)*is*x2*b2(ijk(i-2,j,k),ijk(s-1,t,u))
            if (s .ne. 1)
     &         g = g + (i-1)*(s-1)*is*b2(ijk(i-2,j,k),ijk(s-2,t,u))
         end if
         if (s .ne. 1)
     &      g = g + 2.0d0*(s-1)*is*x1*b2(ijk(i-1,j,k),ijk(s-2,t,u))
      end if
      if (jt .ne. 0) then
         f = f + jt*b1(ijk(i,j-1,k),ijk(s,t-1,u))
         g = g + 4.0d0*jt*y1*y2*b2(ijk(i,j-1,k),ijk(s,t-1,u))
         if (j .ne. 1) then
            g = g + 2.0d0*(j-1)*jt*y2*b2(ijk(i,j-2,k),ijk(s,t-1,u))
            if (t .ne. 1)
     &         g = g + (j-1)*(t-1)*jt*b2(ijk(i,j-2,k),ijk(s,t-2,u))
         end if
         if (t .ne. 1)
     &      g = g + 2.0d0*(t-1)*jt*y1*b2(ijk(i,j-1,k),ijk(s,t-2,u))
      end if
      if (ku .ne. 0) then
         f = f + ku*b1(ijk(i,j,k-1),ijk(s,t,u-1))
         g = g + 4.0d0*ku*z1*z2*b2(ijk(i,j,k-1),ijk(s,t,u-1))
         if (k .ne. 1) then
            g = g + 2.0d0*(k-1)*ku*z2*b2(ijk(i,j,k-2),ijk(s,t,u-1))
            if (u .ne. 1)
     &         g = g + (k-1)*(u-1)*ku*b2(ijk(i,j,k-2),ijk(s,t,u-2))
         end if
         if (u .ne. 1)
     &      g = g + 2.0d0*(u-1)*ku*z1*b2(ijk(i,j,k-1),ijk(s,t,u-2))
      end if
      if (it .ne. 0) then
         g = g + 4.0d0*it*x1*y2*b2(ijk(i-1,j,k),ijk(s,t-1,u))
         if (i .ne. 1) then
            g = g + 2.0d0*(i-1)*it*y2*b2(ijk(i-2,j,k),ijk(s,t-1,u))
            if (t .ne. 1)
     &         g = g + (i-1)*(t-1)*it*b2(ijk(i-2,j,k),ijk(s,t-2,u))
         end if
         if (t .ne. 1)
     &      g = g + 2.0d0*(t-1)*it*x1*b2(ijk(i-1,j,k),ijk(s,t-2,u))
      end if
      if (iu .ne. 0) then
         g = g + 4.0d0*iu*x1*z2*b2(ijk(i-1,j,k),ijk(s,t,u-1))
         if (i .ne. 1) then
            g = g + 2.0d0*(i-1)*iu*z2*b2(ijk(i-2,j,k),ijk(s,t,u-1))
            if (u .ne. 1)
     &         g = g + (i-1)*(u-1)*iu*b2(ijk(i-2,j,k),ijk(s,t,u-2))
         end if
         if (u .ne. 1)
     &      g = g + 2.0d0*(u-1)*iu*x1*b2(ijk(i-1,j,k),ijk(s,t,u-2))
      end if
      if (js .ne. 0) then
         g = g + 4.0d0*js*y1*x2*b2(ijk(i,j-1,k),ijk(s-1,t,u))
         if (j .ne. 1) then
            g = g + 2.0d0*(j-1)*js*x2*b2(ijk(i,j-2,k),ijk(s-1,t,u))
            if (s .ne. 1)
     &         g = g + (j-1)*(s-1)*js*b2(ijk(i,j-2,k),ijk(s-2,t,u))
         end if
         if (s .ne. 1)
     &      g = g + 2.0d0*(s-1)*js*y1*b2(ijk(i,j-1,k),ijk(s-2,t,u))
      end if
      if (ju .ne. 0) then
         g = g + 4.0d0*ju*y1*z2*b2(ijk(i,j-1,k),ijk(s,t,u-1))
         if (j .ne. 1) then
            g = g + 2.0d0*(j-1)*ju*z2*b2(ijk(i,j-2,k),ijk(s,t,u-1))
            if (u .ne. 1)
     &         g = g + (j-1)*(u-1)*ju*b2(ijk(i,j-2,k),ijk(s,t,u-2))
         end if
         if (u .ne. 1)
     &      g = g + 2.0d0*(u-1)*ju*y1*b2(ijk(i,j-1,k),ijk(s,t,u-2))
      end if
      if (ks .ne. 0) then
         g = g + 4.0d0*ks*z1*x2*b2(ijk(i,j,k-1),ijk(s-1,t,u))
         if (k .ne. 1) then
            g = g + 2.0d0*(k-1)*ks*x2*b2(ijk(i,j,k-2),ijk(s-1,t,u))
            if (s .ne. 1)
     &         g = g + (k-1)*(s-1)*ks*b2(ijk(i,j,k-2),ijk(s-2,t,u))
         end if
         if (s .ne. 1)
     &      g = g + 2.0d0*(s-1)*ks*z1*b2(ijk(i,j,k-1),ijk(s-2,t,u))
      end if
      if (kt .ne. 0) then
         g = g + 4.0d0*kt*z1*y2*b2(ijk(i,j,k-1),ijk(s,t-1,u))
         if (k .ne. 1) then
            g = g + 2.0d0*(k-1)*kt*y2*b2(ijk(i,j,k-2),ijk(s,t-1,u))
            if (t .ne. 1)
     &         g = g + (k-1)*(t-1)*kt*b2(ijk(i,j,k-2),ijk(s,t-2,u))
         end if
         if (t .ne. 1)
     &      g = g + 2.0d0*(t-1)*kt*z1*b2(ijk(i,j,k-1),ijk(s,t-2,u))
      end if
      f = dble(2*n-1) * f
      g = dble(n-1) * g
      d1d2 = (f-g) / dble(n)
      return
      end
