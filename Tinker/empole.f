c
c
c     #############################################################
c     ##  COPYRIGHT (C) 1999 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine empole  --  atomic multipole moment energy  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "empole" calculates the electrostatic energy due to atomic
c     multipole interactions
c
c
      subroutine empole
      use limits
      implicit none
c
c
c     choose the method for summing over multipole interactions
c
      if (use_ewald) then
         if (use_mlist) then
            call empole0d
         else
            call empole0c
         end if
      else
         if (use_mlist) then
            call empole0b
         else
            call empole0a
         end if
      end if
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine empole0a  --  double loop multipole energy  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "empole0a" calculates the atomic multipole interaction energy
c     using a double loop
c
c
      subroutine empole0a
      use sizes
      use atomid 
      use atoms
      use bound
      use cell
      use polar 
      use chgpot
      use couple
      use energi
      use group
      use math
      use mplpot
      use mpole
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer ii,kk
      integer ix,iy,iz
      integer kx,ky,kz
      real*8 e,f,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 qrix,qriy,qriz
      real*8 qrkx,qrky,qrkz
      real*8 dri,drk,dik
      real*8 qrri,qrrk
      real*8 qrrik,qik
      real*8 diqrk,dkqri
      real*8 term1,term2,term3
      real*8 term4,term5
!     charge penetration varables
      real*8 alphai,alphak
      real*8 dampi,dampk
      real*8 expdampi,expdampk
      real*8 termi,termk
      real*8, allocatable ::  scalei(:),scalek(:)
      real*8, allocatable ::  scaleik(:)
      real*8 nuci,qi,nuck,qk
!     end 
      real*8, allocatable :: mscale(:)
      logical proceed,usei,usek
      character*6 mode
c
c
c     zero out the total atomic multipole energy
c
      em = 0.0d0
      if (npole .eq. 0)  return
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
c
c     maximum order to rr9 
c
      allocate (scalei(9))
      allocate (scalek(9))
      allocate (scaleik(9))
c
c     maximum order to rr9
c
      do i = 1, 9 
         scalei(i) = 1.0d0
         scalek(i) = 1.0d0
         scaleik(i) = 1.0d0
      end do
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         mscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'MPOLE'
      call switch (mode)
c
c     calculate the multipole interaction energy term
c
      do i = 1, npole-1
         ii = ipole(i)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = yaxis(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         usei = (use(ii) .or. use(iz) .or. use(ix) .or. use(iy))
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
         end do
c
c     evaluate all sites within the cutoff distance
c
         do k = i+1, npole
            kk = ipole(k)
            kz = zaxis(k)
            kx = xaxis(k)
            ky = yaxis(k)
            usek = (use(kk) .or. use(kz) .or. use(kx) .or. use(ky))
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (.not. use_intra)  proceed = .true.
            if (proceed)  proceed = (usei .or. usek)
            if (proceed) then
               xr = x(kk) - xi
               yr = y(kk) - yi
               zr = z(kk) - zi
               if (use_bounds)  call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  ck = rpole(1,k)
                  dkx = rpole(2,k)
                  dky = rpole(3,k)
                  dkz = rpole(4,k)
                  qkxx = rpole(5,k)
                  qkxy = rpole(6,k)
                  qkxz = rpole(7,k)
                  qkyy = rpole(9,k)
                  qkyz = rpole(10,k)
                  qkzz = rpole(13,k)
c
c     get reciprocal distance terms for this interaction
c
                  rr1 = f * mscale(kk) / r
                  rr3 = rr1 / r2
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr9 = 7.0d0 * rr7 / r2
c
c     intermediates involving moments and distance separation
c
                  dri = dix*xr + diy*yr + diz*zr
                  drk = dkx*xr + dky*yr + dkz*zr
                  dik = dix*dkx + diy*dky + diz*dkz
                  qrix = qixx*xr + qixy*yr + qixz*zr
                  qriy = qixy*xr + qiyy*yr + qiyz*zr
                  qriz = qixz*xr + qiyz*yr + qizz*zr
                  qrkx = qkxx*xr + qkxy*yr + qkxz*zr
                  qrky = qkxy*xr + qkyy*yr + qkyz*zr
                  qrkz = qkxz*xr + qkyz*yr + qkzz*zr
                  qrri = qrix*xr + qriy*yr + qriz*zr
                  qrrk = qrkx*xr + qrky*yr + qrkz*zr
                  qrrik = qrix*qrkx + qriy*qrky + qriz*qrkz
                  qik = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                     + qixx*qkxx + qiyy*qkyy + qizz*qkzz
                  diqrk = dix*qrkx + diy*qrky + diz*qrkz
                  dkqri = dkx*qrix + dky*qriy + dkz*qriz
c
c     read in charge penetration damping parameters
c
                  alphai = penalpha(type(ii))
                  alphak = penalpha(type(kk))
c
c     compute common factors for damping
c
                  dampi = alphai*r
                  dampk = alphak*r
                  expdampi = exp(-dampi)
                  expdampk = exp(-dampk)
c
c     calculate one-site scale factors
c
                  scalei(1) = 1.0d0 - expdampi
                  scalek(1) = 1.0d0 - expdampk
                  scalei(3) = 1.0d0 - (1.0d0 + dampi)*expdampi
                  scalek(3) = 1.0d0 - (1.0d0 +dampk)*expdampk
                  scalei(5) = 1.0d0 - (1.0d0 + dampi + 
     &                 (1.0d0/3.0d0)*dampi**2)*expdampi
                  scalek(5) = 1.0d0-(1.0d0 +dampk +
     &                 (1.0d0/3.0d0)*dampk**2)*expdampk
                  scalei(7) = 1.0d0-(1.0d0 + dampi + 0.4d0*dampi**2 +
     &                 (1.0d0/15.0d0)*dampi**3)*expdampi
                  scalek(7) = 1.0d0-(1.0d0 + dampk + 0.4d0*dampk**2 +
     &                 (1.0d0/15.0d0)*dampk**3)*expdampk
c
c     calculate two-site scale factors
c
                  if (alphai .ne. alphak) then
                     termi = alphak**2/(alphak**2 - alphai**2)
                     termk = alphai**2/(alphai**2 - alphak**2)
                     scaleik(1) =1.0d0-termi*expdampi -termk*expdampk
                     scaleik(3) =1.0d0-termi*(1.0d0 +dampi)*expdampi
     &                           - termk*(1.0d0 + dampk)*expdampk
                     scaleik(5) = 1.0d0 - termi*(1.0d0 + dampi +
     &                    (1.0d0/3.0d0)*dampi**2)*expdampi -
     &                    termk*(1.0d0 + dampk +
     &                    (1.0d0/3.0d0)*dampk**2)*expdampk
                     scaleik(7) = 1.0d0 - termi*(1.0d0 + dampi +
     &                    0.4d0*dampi**2 + (1.0d0/15.0d0)*dampi**3)*
     &                    expdampi -
     &                    termk*(1.0d0 + dampk +
     &                    0.4d0*dampk**2 + (1.0d0/15.0d0)*dampk**3)*
     &                    expdampk
                     scaleik(9) = 1.0d0 - termi*(1.0d0 + dampi +
     &                    (3.0d0/7.0d0)*dampi**2 +
     &                    (2.0d0/21.0d0)*dampi**3 +
     &                    (1.0d0/105.0d0)*dampi**4)*expdampi -
     &                    termk*(1.0d0 + dampk +
     &                    (3.0d0/7.0d0)*dampk**2 +
     &                    (2.0d0/21.0d0)*dampk**3 +
     &                    (1.0d0/105.0d0)*dampk**4)*expdampk
                  else
                     scaleik(1) = 1.0d0 -(1.0d0+0.5d0*dampi)*expdampi
                     scaleik(3) = 1.0d0 -(1.0d0+dampi+0.5d0*dampi**2)
     &                    *expdampi
                     scaleik(5) = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2
     &                    + (1.0d0/6.0d0)*dampi**3)*expdampi
                     scaleik(7) = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2
     &                    + (1.0d0/6.0d0)*dampi**3
     &                    + (1.0d0/30.0d0)*dampi**4)*expdampi
                     scaleik(9) = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2
     &                    + (1.0d0/6.0d0)*dampi**3
     &                    + (4.0d0/105.0d0)*dampi**4
     &                    + (1.0d0/210.0d0)*dampi**5)*expdampi
                  end if
            
                  nuci = atomic(i)
                  nuck = atomic(k)
                  qi = ci - nuci
                  qk = ck - nuck

                  term1 = nuci*nuck + nuci*qk*scalek(1) 
     &                    + nuck*qi*scalei(1) + qi*qk*scaleik(1)
                  term2 = nuck*dri*scalei(3) - nuci*drk*scalek(3)
     &                    + (dik + qk*dri - qi*drk)*scaleik(3)
                  term3 = nuci*qrrk*scalek(5) + nuck*qrri*scalei(5)+
     &                    (-dri*drk + 2.0d0*(dkqri-diqrk+qik) + 
     &                    qi*qrrk + qk*qrri)*scaleik(5) 
                  term4 = (dri*qrrk-drk*qrri-4.0d0*qrrik)*scaleik(7)
                  term5 = qrri*qrrk*scaleik(9)
c
c     compute the energy contribution for this interaction
c
                  e = term1*rr1 + term2*rr3 + term3*rr5
     &                   + term4*rr7 + term5*rr9
                  if (use_group)  e = e * fgrp
                  em = em + e
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (use_replica) then
c
c     calculate interaction energy with other unit cells
c
         do i = 1, npole
            ii = ipole(i)
            iz = zaxis(i)
            ix = xaxis(i)
            iy = yaxis(i)
            xi = x(ii)
            yi = y(ii)
            zi = z(ii)
            ci = rpole(1,i)
            dix = rpole(2,i)
            diy = rpole(3,i)
            diz = rpole(4,i)
            qixx = rpole(5,i)
            qixy = rpole(6,i)
            qixz = rpole(7,i)
            qiyy = rpole(9,i)
            qiyz = rpole(10,i)
            qizz = rpole(13,i)
            usei = (use(ii) .or. use(iz) .or. use(ix) .or. use(iy))
            do j = 1, n12(ii)
               mscale(i12(j,ii)) = m2scale
            end do
            do j = 1, n13(ii)
               mscale(i13(j,ii)) = m3scale
            end do
            do j = 1, n14(ii)
               mscale(i14(j,ii)) = m4scale
            end do
            do j = 1, n15(ii)
               mscale(i15(j,ii)) = m5scale
            end do
c
c     evaluate all sites within the cutoff distance
c
            do k = i, npole
               kk = ipole(k)
               kz = zaxis(k)
               kx = xaxis(k)
               ky = yaxis(k)
               usek = (use(kk) .or. use(kz) .or. use(kx) .or. use(ky))
               if (use_group)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
               proceed = .true.
               if (proceed)  proceed = (usei .or. usek)
               if (proceed) then
                  do j = 1, ncell
                     xr = x(kk) - xi
                     yr = y(kk) - yi
                     zr = z(kk) - zi
                     call imager (xr,yr,zr,j)
                     r2 = xr*xr + yr* yr + zr*zr
                     if (.not. (use_polymer .and. r2.le.polycut2))
     &                  mscale(kk) = 1.0d0
                     if (r2 .le. off2) then
                        r = sqrt(r2)
                        ck = rpole(1,k)
                        dkx = rpole(2,k)
                        dky = rpole(3,k)
                        dkz = rpole(4,k)
                        qkxx = rpole(5,k)
                        qkxy = rpole(6,k)
                        qkxz = rpole(7,k)
                        qkyy = rpole(9,k)
                        qkyz = rpole(10,k)
                        qkzz = rpole(13,k)
c
c     get reciprocal distance terms for this interaction
c
                        rr1 = f / r
                        rr3 = rr1 / r2
                        rr5 = 3.0d0 * rr3 / r2
                        rr7 = 5.0d0 * rr5 / r2
                        rr9 = 7.0d0 * rr7 / r2
c
c     intermediates involving moments and distance separation
c
                        dri = dix*xr + diy*yr + diz*zr
                        drk = dkx*xr + dky*yr + dkz*zr
                        dik = dix*dkx + diy*dky + diz*dkz
                        qrix = qixx*xr + qixy*yr + qixz*zr
                        qriy = qixy*xr + qiyy*yr + qiyz*zr
                        qriz = qixz*xr + qiyz*yr + qizz*zr
                        qrkx = qkxx*xr + qkxy*yr + qkxz*zr
                        qrky = qkxy*xr + qkyy*yr + qkyz*zr
                        qrkz = qkxz*xr + qkyz*yr + qkzz*zr
                        qrri = qrix*xr + qriy*yr + qriz*zr
                        qrrk = qrkx*xr + qrky*yr + qrkz*zr
                        qrrik = qrix*qrkx + qriy*qrky + qriz*qrkz
                        qik = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                           + qixx*qkxx + qiyy*qkyy + qizz*qkzz
                        diqrk = dix*qrkx + diy*qrky + diz*qrkz
                        dkqri = dkx*qrix + dky*qriy + dkz*qriz
c
c     read in charge penetration damping parameters
c
                        alphai = penalpha(type(ii))
                        alphak = penalpha(type(kk))
c
c     compute common factors for damping
c
                        dampi = alphai*r
                        dampk = alphak*r
                        expdampi = exp(-dampi)
                        expdampk = exp(-dampk)
c
c     calculate one-site scale factors
c
                        scalei(1) = 1.0d0 - expdampi
                        scalek(1) = 1.0d0 - expdampk
                        scalei(3) = 1.0d0 - (1.0d0 + dampi)*expdampi
                        scalek(3) = 1.0d0 - (1.0d0 +dampk)*expdampk
                        scalei(5) = 1.0d0 - (1.0d0 + dampi + 
     &                       (1.0d0/3.0d0)*dampi**2)*expdampi
                        scalek(5) = 1.0d0-(1.0d0 +dampk +
     &                       (1.0d0/3.0d0)*dampk**2)*expdampk
                        scalei(7) = 1.0d0-(1.0d0 + dampi + 
     &                        0.4d0*dampi**2 +
     &                       (1.0d0/15.0d0)*dampi**3)*expdampi
                        scalek(7) = 1.0d0-(1.0d0 + dampk + 
     &                        0.4d0*dampk**2 +
     &                       (1.0d0/15.0d0)*dampk**3)*expdampk
c
c     calculate two-site scale factors
c
                        if (alphai .ne. alphak) then
                           termi = alphak**2/(alphak**2 - alphai**2)
                           termk = alphai**2/(alphai**2 - alphak**2)
                           scaleik(1) = 1.0d0 - termi*expdampi - 
     &                               termk*expdampk
                           scaleik(3) = 1.0d0 - termi*(1.0d0+dampi)*
     &                               expdampi - termk*(1.0d0 + 
     &                               dampk)*expdampk
                           scaleik(5) = 1.0d0- termi*(1.0d0 + dampi +
     &                          (1.0d0/3.0d0)*dampi**2)*expdampi -
     &                          termk*(1.0d0 + dampk +
     &                          (1.0d0/3.0d0)*dampk**2)*expdampk
                           scaleik(7) = 1.0d0-termi*(1.0d0 + dampi +
     &                       0.4d0*dampi**2+(1.0d0/15.0d0)*dampi**3)*
     &                       expdampi - termk*(1.0d0 + dampk +
     &                       0.4d0*dampk**2+(1.0d0/15.0d0)*dampk**3)*
     &                       expdampk
                           scaleik(9) =1.0d0 - termi*(1.0d0 + dampi +
     &                          (3.0d0/7.0d0)*dampi**2 +
     &                          (2.0d0/21.0d0)*dampi**3 +
     &                          (1.0d0/105.0d0)*dampi**4)*expdampi -
     &                          termk*(1.0d0 + dampk +
     &                          (3.0d0/7.0d0)*dampk**2 +
     &                          (2.0d0/21.0d0)*dampk**3 +
     &                          (1.0d0/105.0d0)*dampk**4)*expdampk
                        else
                           scaleik(1)=1.0d0-(1.0d0+0.5d0*dampi)
     &                                *expdampi
                           scaleik(3)=1.0d0-(1.0d0+dampi+
     &                                0.5d0*dampi**2)*expdampi
                           scaleik(5) = 1.0d0 - (1.0d0+dampi + 
     &                          0.5d0*dampi**2+ (1.0d0/6.0d0)
     &                          *dampi**3)*expdampi 
                           scaleik(7) = 1.0d0 - (1.0d0+dampi + 
     &                          0.5d0*dampi**2
     &                          + (1.0d0/6.0d0)*dampi**3
     &                          + (1.0d0/30.0d0)*dampi**4)*expdampi
                           scaleik(9) = 1.0d0 - (1.0d0+dampi + 
     &                            0.5d0*dampi**2
     &                          + (1.0d0/6.0d0)*dampi**3
     &                          + (4.0d0/105.0d0)*dampi**4
     &                          + (1.0d0/210.0d0)*dampi**5)*expdampi
                        end if
            
                        nuci = atomic(i)
                        nuck = atomic(k)
                        qi = ci - nuci
                        qk = ck - nuck

                        term1 = nuci*nuck + nuci*qk*scalek(1) 
     &                        + nuck*qi*scalei(1) + qi*qk*scaleik(1)
                        term2 =nuck*dri*scalei(3)-nuci*drk*scalek(3)
     &                        + (dik + qk*dri - qi*drk)*scaleik(3)
                        term3=nuci*qrrk*scalek(5)+nuck*qrri*scalei(5)
     &                        +(-dri*drk + 2.0d0*(dkqri-diqrk+qik) + 
     &                         qi*qrrk + qk*qrri)*scaleik(5) 
                        term4 = (dri*qrrk-drk*qrri-4.0d0*qrrik)
     &                          *scaleik(7)
                        term5 = qrri*qrrk*scaleik(9)
c
c     compute the energy contribution for this interaction
c
                        e = term1*rr1 + term2*rr3 + term3*rr5
     &                         + term4*rr7 + term5*rr9
                        e = e * mscale(kk)
                        if (use_group)  e = e * fgrp
                        if (ii .eq. kk)  e = 0.5d0 * e
                        em = em + e
                     end if
                  end do
               end if
            end do
c
c     reset exclusion coefficients for connected atoms
c
            do j = 1, n12(ii)
               mscale(i12(j,ii)) = 1.0d0
            end do
            do j = 1, n13(ii)
               mscale(i13(j,ii)) = 1.0d0
            end do
            do j = 1, n14(ii)
               mscale(i14(j,ii)) = 1.0d0
            end do
            do j = 1, n15(ii)
               mscale(i15(j,ii)) = 1.0d0
            end do
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine empole0b  --  neighbor list multipole energy  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "empole0b" calculates the atomic multipole interaction energy
c     using a neighbor list
c
c
      subroutine empole0b
      use sizes
      use atomid 
      use atoms
      use bound
      use polar 
      use chgpot
      use couple
      use energi
      use group
      use math
      use mplpot
      use mpole
      use neigh
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      integer ix,iy,iz
      integer kx,ky,kz
      real*8 e,f,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 qrix,qriy,qriz
      real*8 qrkx,qrky,qrkz
      real*8 dri,drk,dik
      real*8 qrri,qrrk
      real*8 qrrik,qik
      real*8 diqrk,dkqri
      real*8 term1,term2,term3
      real*8 term4,term5
!     charge penetration varables
      real*8 alphai,alphak
      real*8 dampi,dampk
      real*8 expdampi,expdampk
      real*8 termi,termk
      real*8, allocatable ::  scalei(:),scalek(:)
      real*8, allocatable ::  scaleik(:)
      real*8 nuci,qi,nuck,qk
!     end 
      real*8, allocatable :: mscale(:)
      logical proceed,usei,usek
      character*6 mode
c
c
c     zero out the total atomic multipole energy
c
      em = 0.0d0
      if (npole .eq. 0)  return
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         mscale(i) = 1.0d0
      end do
c
c     maximum order to rr9 
c
      allocate (scalei(9))
      allocate (scalek(9))
      allocate (scaleik(9))
c
c     maximum order to rr9
c
      do i = 1, 9
         scalei(i) = 1.0d0
         scalek(i) = 1.0d0
         scaleik(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'MPOLE'
      call switch (mode)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private)
!$OMP& shared(npole,ipole,x,y,z,xaxis,yaxis,zaxis,rpole,use,
!$OMP& n12,i12,n13,i13,n14,i14,n15,i15,m2scale,m3scale,m4scale,
!$OMP& atomic,penalpha,type,
!$OMP& m5scale,nelst,elst,use_group,use_intra,use_bounds,off2,f)
!$OMP& firstprivate(mscale) shared (em)
!$OMP DO reduction(+:em) schedule(guided)
c
c     compute the real space portion of the Ewald summation
c
      do i = 1, npole
         ii = ipole(i)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = yaxis(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         usei = (use(ii) .or. use(iz) .or. use(ix) .or. use(iy))
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
         end do
c
c     evaluate all sites within the cutoff distance
c
         do kkk = 1, nelst(i)
            k = elst(kkk,i)
            kk = ipole(k)
            kz = zaxis(k)
            kx = xaxis(k)
            ky = yaxis(k)
            usek = (use(kk) .or. use(kz) .or. use(kx) .or. use(ky))
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (.not. use_intra)  proceed = .true.
            if (proceed)  proceed = (usei .or. usek)
            if (proceed) then
               xr = x(kk) - xi
               yr = y(kk) - yi
               zr = z(kk) - zi
               if (use_bounds)  call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  ck = rpole(1,k)
                  dkx = rpole(2,k)
                  dky = rpole(3,k)
                  dkz = rpole(4,k)
                  qkxx = rpole(5,k)
                  qkxy = rpole(6,k)
                  qkxz = rpole(7,k)
                  qkyy = rpole(9,k)
                  qkyz = rpole(10,k)
                  qkzz = rpole(13,k)
c
c     get reciprocal distance terms for this interaction
c
                  rr1 = f * mscale(kk) / r
                  rr3 = rr1 / r2
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr9 = 7.0d0 * rr7 / r2
c
c     intermediates involving moments and distance separation
c
                  dri = dix*xr + diy*yr + diz*zr
                  drk = dkx*xr + dky*yr + dkz*zr
                  dik = dix*dkx + diy*dky + diz*dkz
                  qrix = qixx*xr + qixy*yr + qixz*zr
                  qriy = qixy*xr + qiyy*yr + qiyz*zr
                  qriz = qixz*xr + qiyz*yr + qizz*zr
                  qrkx = qkxx*xr + qkxy*yr + qkxz*zr
                  qrky = qkxy*xr + qkyy*yr + qkyz*zr
                  qrkz = qkxz*xr + qkyz*yr + qkzz*zr
                  qrri = qrix*xr + qriy*yr + qriz*zr
                  qrrk = qrkx*xr + qrky*yr + qrkz*zr
                  qrrik = qrix*qrkx + qriy*qrky + qriz*qrkz
                  qik = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                     + qixx*qkxx + qiyy*qkyy + qizz*qkzz
                  diqrk = dix*qrkx + diy*qrky + diz*qrkz
                  dkqri = dkx*qrix + dky*qriy + dkz*qriz
c
c     read in charge penetration damping parameters
c
                  alphai = penalpha(type(ii))
                  alphak = penalpha(type(kk))
c
c     compute common factors for damping
c
                  dampi = alphai*r
                  dampk = alphak*r
                  expdampi = exp(-dampi)
                  expdampk = exp(-dampk)
c
c     calculate one-site scale factors
c
                  scalei(1) = 1.0d0 - expdampi
                  scalek(1) = 1.0d0 - expdampk
                  scalei(3) = 1.0d0 - (1.0d0 + dampi)*expdampi
                  scalek(3) = 1.0d0 - (1.0d0 +dampk)*expdampk
                  scalei(5) = 1.0d0 - (1.0d0 + dampi + 
     &                 (1.0d0/3.0d0)*dampi**2)*expdampi
                  scalek(5) = 1.0d0-(1.0d0 +dampk +
     &                 (1.0d0/3.0d0)*dampk**2)*expdampk
                  scalei(7) = 1.0d0-(1.0d0 + dampi + 0.4d0*dampi**2 +
     &                 (1.0d0/15.0d0)*dampi**3)*expdampi
                  scalek(7) = 1.0d0-(1.0d0 + dampk + 0.4d0*dampk**2 +
     &                 (1.0d0/15.0d0)*dampk**3)*expdampk
c
c     calculate two-site scale factors
c
                  if (alphai .ne. alphak) then
                     termi = alphak**2/(alphak**2 - alphai**2)
                     termk = alphai**2/(alphai**2 - alphak**2)
                     scaleik(1) =1.0d0-termi*expdampi -termk*expdampk
                     scaleik(3) =1.0d0-termi*(1.0d0 +dampi)*expdampi
     &                           - termk*(1.0d0 + dampk)*expdampk
                     scaleik(5) = 1.0d0 - termi*(1.0d0 + dampi +
     &                    (1.0d0/3.0d0)*dampi**2)*expdampi -
     &                    termk*(1.0d0 + dampk +
     &                    (1.0d0/3.0d0)*dampk**2)*expdampk
                     scaleik(7) = 1.0d0 - termi*(1.0d0 + dampi +
     &                    0.4d0*dampi**2 + (1.0d0/15.0d0)*dampi**3)*
     &                    expdampi -
     &                    termk*(1.0d0 + dampk +
     &                    0.4d0*dampk**2 + (1.0d0/15.0d0)*dampk**3)*
     &                    expdampk
                     scaleik(9) = 1.0d0 - termi*(1.0d0 + dampi +
     &                    (3.0d0/7.0d0)*dampi**2 +
     &                    (2.0d0/21.0d0)*dampi**3 +
     &                    (1.0d0/105.0d0)*dampi**4)*expdampi -
     &                    termk*(1.0d0 + dampk +
     &                    (3.0d0/7.0d0)*dampk**2 +
     &                    (2.0d0/21.0d0)*dampk**3 +
     &                    (1.0d0/105.0d0)*dampk**4)*expdampk
                  else
                     scaleik(1) = 1.0d0 -(1.0d0+0.5d0*dampi)*expdampi
                     scaleik(3) = 1.0d0 -(1.0d0+dampi+0.5d0*dampi**2)
     &                    *expdampi
                     scaleik(5) = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2
     &                    + (1.0d0/6.0d0)*dampi**3)*expdampi
                     scaleik(7) = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2
     &                    + (1.0d0/6.0d0)*dampi**3
     &                    + (1.0d0/30.0d0)*dampi**4)*expdampi
                     scaleik(9) = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2
     &                    + (1.0d0/6.0d0)*dampi**3
     &                    + (4.0d0/105.0d0)*dampi**4
     &                    + (1.0d0/210.0d0)*dampi**5)*expdampi
                  end if
            
                  nuci = atomic(i)
                  nuck = atomic(k)
                  qi = ci - nuci
                  qk = ck - nuck

                  term1 = nuci*nuck + nuci*qk*scalek(1) 
     &                    + nuck*qi*scalei(1) + qi*qk*scaleik(1)
                  term2 = nuck*dri*scalei(3) - nuci*drk*scalek(3)
     &                    + (dik + qk*dri - qi*drk)*scaleik(3)
                  term3 = nuci*qrrk*scalek(5) + nuck*qrri*scalei(5)+
     &                    (-dri*drk + 2.0d0*(dkqri-diqrk+qik) + 
     &                    qi*qrrk + qk*qrri)*scaleik(5) 
                  term4 = (dri*qrrk-drk*qrri-4.0d0*qrrik)*scaleik(7)
                  term5 = qrri*qrrk*scaleik(9)
c
c     compute the energy contribution for this interaction
c
                  e = term1*rr1 + term2*rr3 + term3*rr5
     &                   + term4*rr7 + term5*rr9
                  if (use_group)  e = e * fgrp
                  em = em + e
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine empole0c  --  Ewald multipole energy via loop  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "empole0c" calculates the atomic multipole interaction energy
c     using particle mesh Ewald summation and a double loop
c
c
      subroutine empole0c
      use sizes
      use atoms
      use boxes
      use chgpot
      use energi
      use ewald
      use math
      use mpole
      implicit none
      integer i,ii
      real*8 e,f
      real*8 term,fterm
      real*8 cii,dii,qii
      real*8 xd,yd,zd
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
c
c
c     zero out the total atomic multipole energy
c
      em = 0.0d0
      if (npole .eq. 0)  return
c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     compute the real space portion of the Ewald summation
c
      call emreal0c
c
c     compute the reciprocal space part of the Ewald summation
c
      call emrecip
c
c     compute the self-energy portion of the Ewald summation
c
      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / sqrtpi
      do i = 1, npole
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         cii = ci*ci
         dii = dix*dix + diy*diy + diz*diz
         qii = 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
     &            + qixx*qixx + qiyy*qiyy + qizz*qizz
         e = fterm * (cii + term*(dii/3.0d0+2.0d0*term*qii/5.0d0))
         em = em + e
      end do
c
c     compute the cell dipole boundary correction term
c
      if (boundary .eq. 'VACUUM') then
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            dix = rpole(2,i)
            diy = rpole(3,i)
            diz = rpole(4,i)
            xd = xd + dix + rpole(1,i)*x(ii)
            yd = yd + diy + rpole(1,i)*y(ii)
            zd = zd + diz + rpole(1,i)*z(ii)
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         e = term * (xd*xd+yd*yd+zd*zd)
         em = em + e
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine emreal0c  --  real space mpole energy via loop  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "emreal0c" evaluates the real space portion of the Ewald sum
c     energy due to atomic multipoles using a double loop
c
c     literature reference:
c
c     W. Smith, "Point Multipoles in the Ewald Summation (Revisited)",
c     CCP5 Newsletter, 46, 18-30, 1998  (see http://www.ccp5.org/)
c
c
      subroutine emreal0c
      use sizes
      use atomid 
      use atoms
      use bound
      use cell
      use polar 
      use chgpot
      use couple
      use energi
      use ewald
      use math
      use mplpot
      use mpole
      use shunt
      implicit none
      integer i,j,k
      integer ii,kk
      integer jcell
      real*8 e,f,bfac,erfc
      real*8 alsq2,alsq2n
      real*8 exp2a,ralpha
      real*8 scalekk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 qrix,qriy,qriz
      real*8 qrkx,qrky,qrkz
      real*8 dri,drk,dik
      real*8 qrri,qrrk
      real*8 qrrik,qik
      real*8 diqrk,dkqri
      real*8 term1,term2,term3
      real*8 term4,term5
      real*8 bn(0:4)
!     charge penetration varables
      real*8 alphai,alphak
      real*8 dampi,dampk
      real*8 expdampi,expdampk
      real*8 termi,termk
      real*8, allocatable ::  scalei(:),scalek(:)
      real*8, allocatable ::  scaleik(:)
      real*8 nuci,qi,nuck,qk
!     end 
      real*8, allocatable :: mscale(:)
      character*6 mode
      external erfc
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         mscale(i) = 1.0d0
      end do
c
c     maximum order to rr9
c
      allocate (scalei(9))
      allocate (scalek(9))
      allocate (scaleik(9))
c
c     maximum order to rr9 
c
      do i = 1, 9 
         scalei(i) = 1.0d0
         scalek(i) = 1.0d0
         scaleik(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     compute the real space portion of the Ewald summation
c
      do i = 1, npole-1
         ii = ipole(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
         end do
c
c     evaluate all sites within the cutoff distance
c
         do k = i+1, npole
            kk = ipole(k)
            xr = x(kk) - xi
            yr = y(kk) - yi
            zr = z(kk) - zi
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
               ck = rpole(1,k)
               dkx = rpole(2,k)
               dky = rpole(3,k)
               dkz = rpole(4,k)
               qkxx = rpole(5,k)
               qkxy = rpole(6,k)
               qkxz = rpole(7,k)
               qkyy = rpole(9,k)
               qkyz = rpole(10,k)
               qkzz = rpole(13,k)
c
c     get reciprocal distance terms for this interaction
c
               rr1 = f / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
c
c     calculate the real space Ewald error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 4
                  bfac = dble(j+j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
               do j = 0, 4
                  bn(j) = f * bn(j)
               end do
c
c     intermediates involving moments and distance separation
c
               dri = dix*xr + diy*yr + diz*zr
               drk = dkx*xr + dky*yr + dkz*zr
               dik = dix*dkx + diy*dky + diz*dkz
               qrix = qixx*xr + qixy*yr + qixz*zr
               qriy = qixy*xr + qiyy*yr + qiyz*zr
               qriz = qixz*xr + qiyz*yr + qizz*zr
               qrkx = qkxx*xr + qkxy*yr + qkxz*zr
               qrky = qkxy*xr + qkyy*yr + qkyz*zr
               qrkz = qkxz*xr + qkyz*yr + qkzz*zr
               qrri = qrix*xr + qriy*yr + qriz*zr
               qrrk = qrkx*xr + qrky*yr + qrkz*zr
               qrrik = qrix*qrkx + qriy*qrky + qriz*qrkz
               qik = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                  + qixx*qkxx + qiyy*qkyy + qizz*qkzz
               diqrk = dix*qrkx + diy*qrky + diz*qrkz
               dkqri = dkx*qrix + dky*qriy + dkz*qriz
c
c     read in charge penetration damping parameters
c
               alphai = penalpha(type(ii))
               alphak = penalpha(type(kk))
c
c     compute common factors for damping
c
               dampi = alphai*r
               dampk = alphak*r
               expdampi = exp(-dampi)
               expdampk = exp(-dampk)
c
c     calculate one-site scale factors
c
               scalei(1) = 1.0d0 - expdampi
               scalek(1) = 1.0d0 - expdampk
               scalei(3) = 1.0d0 - (1.0d0 + dampi)*expdampi
               scalek(3) = 1.0d0 - (1.0d0 +dampk)*expdampk
               scalei(5) = 1.0d0 - (1.0d0 + dampi + 
     &              (1.0d0/3.0d0)*dampi**2)*expdampi
               scalek(5) = 1.0d0-(1.0d0 +dampk +
     &              (1.0d0/3.0d0)*dampk**2)*expdampk
               scalei(7) = 1.0d0-(1.0d0 + dampi + 0.4d0*dampi**2 +
     &              (1.0d0/15.0d0)*dampi**3)*expdampi
               scalek(7) = 1.0d0-(1.0d0 + dampk + 0.4d0*dampk**2 +
     &              (1.0d0/15.0d0)*dampk**3)*expdampk
c
c     calculate two-site scale factors
c
               if (alphai .ne. alphak) then
                  termi = alphak**2/(alphak**2 - alphai**2)
                  termk = alphai**2/(alphai**2 - alphak**2)
                  scaleik(1) =1.0d0-termi*expdampi -termk*expdampk
                  scaleik(3) =1.0d0-termi*(1.0d0 +dampi)*expdampi
     &                        - termk*(1.0d0 + dampk)*expdampk
                  scaleik(5) = 1.0d0 - termi*(1.0d0 + dampi +
     &                 (1.0d0/3.0d0)*dampi**2)*expdampi -
     &                 termk*(1.0d0 + dampk +
     &                 (1.0d0/3.0d0)*dampk**2)*expdampk
                  scaleik(7) = 1.0d0 - termi*(1.0d0 + dampi +
     &                 0.4d0*dampi**2 + (1.0d0/15.0d0)*dampi**3)*
     &                 expdampi -
     &                 termk*(1.0d0 + dampk +
     &                 0.4d0*dampk**2 + (1.0d0/15.0d0)*dampk**3)*
     &                 expdampk
                  scaleik(9) = 1.0d0 - termi*(1.0d0 + dampi +
     &                 (3.0d0/7.0d0)*dampi**2 +
     &                 (2.0d0/21.0d0)*dampi**3 +
     &                 (1.0d0/105.0d0)*dampi**4)*expdampi -
     &                 termk*(1.0d0 + dampk +
     &                 (3.0d0/7.0d0)*dampk**2 +
     &                 (2.0d0/21.0d0)*dampk**3 +
     &                 (1.0d0/105.0d0)*dampk**4)*expdampk
               else
                  scaleik(1) = 1.0d0 -(1.0d0+0.5d0*dampi)*expdampi
                  scaleik(3) = 1.0d0 -(1.0d0+dampi+0.5d0*dampi**2)
     &                 *expdampi
                  scaleik(5) = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2
     &                 + (1.0d0/6.0d0)*dampi**3)*expdampi
                  scaleik(7) = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2
     &                 + (1.0d0/6.0d0)*dampi**3
     &                 + (1.0d0/30.0d0)*dampi**4)*expdampi
                  scaleik(9) = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2
     &                 + (1.0d0/6.0d0)*dampi**3
     &                 + (4.0d0/105.0d0)*dampi**4
     &                 + (1.0d0/210.0d0)*dampi**5)*expdampi
               end if
            
               nuci = atomic(i)
               nuck = atomic(k)
               qi = ci - nuci
               qk = ck - nuck
c
c     modify distances to account for Ewald and exclusions
c
               bn(0) = bn(0)/rr1
               bn(1) = bn(1)/rr3
               bn(2) = bn(2)/rr5
               bn(3) = bn(3)/rr7
               bn(4) = bn(4)/rr9

               term1 = nuci*nuck*(bn(0)-(1.0d0-mscale(kk))) 
     &             + nuci*qk*(bn(0)-(1.0d0-scalek(1)*mscale(kk))) 
     &             + nuck*qi*(bn(0)-(1.0d0-scalei(1)*mscale(kk))) 
     &             + qi*qk*(bn(0)-(1.0d0-scaleik(1)*mscale(kk)))

               term2 = nuck*dri*(bn(1)-(1.0d0-scalei(3)*mscale(kk))) 
     &             - nuci*drk*(bn(1)-(1.0d0-scalek(3)*mscale(kk)))
     &             + (dik + qk*dri - qi*drk)
     &               *(bn(1)-(1.0d0-scaleik(3)*mscale(kk)))
               term3 = nuci*qrrk*(bn(2)-(1.0d0-scalek(5)*mscale(kk)))
     &             + nuck*qrri*(bn(2)-(1.0d0-scalei(5)*mscale(kk)))
     &             +(-dri*drk + 2.0d0*(dkqri-diqrk+qik) 
     &             + qi*qrrk + qk*qrri)
     &             * (bn(2)-(1.0d0-scaleik(5)*mscale(kk))) 

               term4 = (dri*qrrk-drk*qrri-4.0d0*qrrik)
     &             *(bn(3)-(1.0d0-scaleik(7)*mscale(kk)))

               term5=qrri*qrrk*(bn(4)-(1.0d0-scaleik(9)*mscale(kk)))
c
c     compute the energy contribution for this interaction
c
               e = term1*rr1 + term2*rr3 + term3*rr5
     &                + term4*rr7 + term5*rr9
               em = em + e
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (use_replica) then
c
c     calculate interaction energy with other unit cells
c
         do i = 1, npole
            ii = ipole(i)
            xi = x(ii)
            yi = y(ii)
            zi = z(ii)
            ci = rpole(1,i)
            dix = rpole(2,i)
            diy = rpole(3,i)
            diz = rpole(4,i)
            qixx = rpole(5,i)
            qixy = rpole(6,i)
            qixz = rpole(7,i)
            qiyy = rpole(9,i)
            qiyz = rpole(10,i)
            qizz = rpole(13,i)
            do j = 1, n12(ii)
               mscale(i12(j,ii)) = m2scale
            end do
            do j = 1, n13(ii)
               mscale(i13(j,ii)) = m3scale
            end do
            do j = 1, n14(ii)
               mscale(i14(j,ii)) = m4scale
            end do
            do j = 1, n15(ii)
               mscale(i15(j,ii)) = m5scale
            end do
c
c     evaluate all sites within the cutoff distance
c
            do k = i, npole
               kk = ipole(k)
               do jcell = 1, ncell
                  xr = x(kk) - xi
                  yr = y(kk) - yi
                  zr = z(kk) - zi
                  call imager (xr,yr,zr,jcell)
                  r2 = xr*xr + yr* yr + zr*zr
                  if (.not. (use_polymer .and. r2.le.polycut2))
     &               mscale(kk) = 1.0d0
                  if (r2 .le. off2) then
                     r = sqrt(r2)
                     ck = rpole(1,k)
                     dkx = rpole(2,k)
                     dky = rpole(3,k)
                     dkz = rpole(4,k)
                     qkxx = rpole(5,k)
                     qkxy = rpole(6,k)
                     qkxz = rpole(7,k)
                     qkyy = rpole(9,k)
                     qkyz = rpole(10,k)
                     qkzz = rpole(13,k)
c
c     get reciprocal distance terms for this interaction
c
                     rr1 = f / r
                     rr3 = rr1 / r2
                     rr5 = 3.0d0 * rr3 / r2
                     rr7 = 5.0d0 * rr5 / r2
                     rr9 = 7.0d0 * rr7 / r2
c
c     calculate the real space Ewald error function terms
c
                     ralpha = aewald * r
                     bn(0) = erfc(ralpha) / r
                     alsq2 = 2.0d0 * aewald**2
                     alsq2n = 0.0d0
                     if (aewald .gt. 0.0d0)
     &                  alsq2n = 1.0d0 / (sqrtpi*aewald)
                     exp2a = exp(-ralpha**2)
                     do j = 1, 4
                        bfac = dble(j+j-1)
                        alsq2n = alsq2 * alsq2n
                        bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
                     end do
                     do j = 0, 4
                        bn(j) = f * bn(j)
                     end do
c
c     intermediates involving moments and distance separation
c
                     dri = dix*xr + diy*yr + diz*zr
                     drk = dkx*xr + dky*yr + dkz*zr
                     dik = dix*dkx + diy*dky + diz*dkz
                     qrix = qixx*xr + qixy*yr + qixz*zr
                     qriy = qixy*xr + qiyy*yr + qiyz*zr
                     qriz = qixz*xr + qiyz*yr + qizz*zr
                     qrkx = qkxx*xr + qkxy*yr + qkxz*zr
                     qrky = qkxy*xr + qkyy*yr + qkyz*zr
                     qrkz = qkxz*xr + qkyz*yr + qkzz*zr
                     qrri = qrix*xr + qriy*yr + qriz*zr
                     qrrk = qrkx*xr + qrky*yr + qrkz*zr
                     qrrik = qrix*qrkx + qriy*qrky + qriz*qrkz
                     qik = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                        + qixx*qkxx + qiyy*qkyy + qizz*qkzz
                     diqrk = dix*qrkx + diy*qrky + diz*qrkz
                     dkqri = dkx*qrix + dky*qriy + dkz*qriz
c
c     read in charge penetration damping parameters
c
                     alphai = penalpha(type(ii))
                     alphak = penalpha(type(kk))
c
c     compute common factors for damping
c
                     dampi = alphai*r
                     dampk = alphak*r
                     expdampi = exp(-dampi)
                     expdampk = exp(-dampk)
c
c     calculate one-site scale factors
c
                     scalei(1) = 1.0d0 - expdampi
                     scalek(1) = 1.0d0 - expdampk
                     scalei(3) = 1.0d0 - (1.0d0 + dampi)*expdampi
                     scalek(3) = 1.0d0 - (1.0d0 +dampk)*expdampk
                     scalei(5) = 1.0d0 - (1.0d0 + dampi + 
     &                    (1.0d0/3.0d0)*dampi**2)*expdampi
                     scalek(5) = 1.0d0-(1.0d0 +dampk +
     &                    (1.0d0/3.0d0)*dampk**2)*expdampk
                     scalei(7) = 1.0d0-(1.0d0 + dampi + 0.4d0*dampi**2 +
     &                    (1.0d0/15.0d0)*dampi**3)*expdampi
                     scalek(7) = 1.0d0-(1.0d0 + dampk + 0.4d0*dampk**2 +
     &                    (1.0d0/15.0d0)*dampk**3)*expdampk
c
c     calculate two-site scale factors
c
                     if (alphai .ne. alphak) then
                        termi = alphak**2/(alphak**2 - alphai**2)
                        termk = alphai**2/(alphai**2 - alphak**2)
                        scaleik(1) =1.0d0-termi*expdampi -termk*expdampk
                        scaleik(3) =1.0d0-termi*(1.0d0 +dampi)*expdampi
     &                              - termk*(1.0d0 + dampk)*expdampk
                        scaleik(5) = 1.0d0 - termi*(1.0d0 + dampi +
     &                       (1.0d0/3.0d0)*dampi**2)*expdampi -
     &                       termk*(1.0d0 + dampk +
     &                       (1.0d0/3.0d0)*dampk**2)*expdampk
                        scaleik(7) = 1.0d0 - termi*(1.0d0 + dampi +
     &                       0.4d0*dampi**2 + (1.0d0/15.0d0)*dampi**3)*
     &                       expdampi -
     &                       termk*(1.0d0 + dampk +
     &                       0.4d0*dampk**2 + (1.0d0/15.0d0)*dampk**3)*
     &                       expdampk
                        scaleik(9) = 1.0d0 - termi*(1.0d0 + dampi +
     &                       (3.0d0/7.0d0)*dampi**2 +
     &                       (2.0d0/21.0d0)*dampi**3 +
     &                       (1.0d0/105.0d0)*dampi**4)*expdampi -
     &                       termk*(1.0d0 + dampk +
     &                       (3.0d0/7.0d0)*dampk**2 +
     &                       (2.0d0/21.0d0)*dampk**3 +
     &                       (1.0d0/105.0d0)*dampk**4)*expdampk
                     else
                        scaleik(1) = 1.0d0 -(1.0d0+0.5d0*dampi)*expdampi
                        scaleik(3) = 1.0d0 -(1.0d0+dampi+0.5d0*dampi**2)
     &                       *expdampi
                        scaleik(5) = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2
     &                       + (1.0d0/6.0d0)*dampi**3)*expdampi
                        scaleik(7) = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2
     &                       + (1.0d0/6.0d0)*dampi**3
     &                       + (1.0d0/30.0d0)*dampi**4)*expdampi
                        scaleik(9) = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2
     &                       + (1.0d0/6.0d0)*dampi**3
     &                       + (4.0d0/105.0d0)*dampi**4
     &                       + (1.0d0/210.0d0)*dampi**5)*expdampi
                     end if
            
                     nuci = atomic(i)
                     nuck = atomic(k)
                     qi = ci - nuci
                     qk = ck - nuck
c
c     modify distances to account for Ewald and exclusions
c
                     bn(0) = bn(0)/rr1
                     bn(1) = bn(1)/rr3
                     bn(2) = bn(2)/rr5
                     bn(3) = bn(3)/rr7
                     bn(4) = bn(4)/rr9

                     term1 = nuci*nuck*(bn(0)-(1.0d0-mscale(kk))) 
     &                   + nuci*qk*(bn(0)-(1.0d0-scalek(1)*mscale(kk))) 
     &                   + nuck*qi*(bn(0)-(1.0d0-scalei(1)*mscale(kk))) 
     &                   + qi*qk*(bn(0)-(1.0d0-scaleik(1)*mscale(kk)))

                     term2 = nuck*dri*(bn(1)-(1.0d0-scalei(3)
     &                   * mscale(kk))) 
     &                   - nuci*drk*(bn(1)-(1.0d0-scalek(3)*mscale(kk)))
     &                   + (dik + qk*dri - qi*drk)
     &                   * (bn(1)-(1.0d0-scaleik(3)*mscale(kk)))
                     term3 = nuci*qrrk*(bn(2)-(1.0d0-scalek(5)
     &                   * mscale(kk)))
     &                   + nuck*qrri*(bn(2)-(1.0d0-scalei(5)
     &                   * mscale(kk)))
     &                   +(-dri*drk + 2.0d0*(dkqri-diqrk+qik) 
     &                   + qi*qrrk + qk*qrri)
     &                   * (bn(2)-(1.0d0-scaleik(5)*mscale(kk))) 

                     term4 = (dri*qrrk-drk*qrri-4.0d0*qrrik)
     &                   * (bn(3)-(1.0d0-scaleik(7)*mscale(kk)))

                     term5=qrri*qrrk*(bn(4)-(1.0d0-scaleik(9)
     &                   * mscale(kk)))
c
c     compute the energy contribution for this interaction
c
                     e = term1*rr1 + term2*rr3 + term3*rr5
     &                      + term4*rr7 + term5*rr9
                     if (ii .eq. kk)  e = 0.5d0 * e
                     em = em + e
                  end if
               end do
            end do
c
c     reset exclusion coefficients for connected atoms
c
            do j = 1, n12(ii)
               mscale(i12(j,ii)) = 1.0d0
            end do
            do j = 1, n13(ii)
               mscale(i13(j,ii)) = 1.0d0
            end do
            do j = 1, n14(ii)
               mscale(i14(j,ii)) = 1.0d0
            end do
            do j = 1, n15(ii)
               mscale(i15(j,ii)) = 1.0d0
            end do
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine empole0d  --  Ewald multipole energy via list  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "empole0d" calculates the atomic multipole interaction energy
c     using particle mesh Ewald summation and a neighbor list
c
c
      subroutine empole0d
      use sizes
      use atoms
      use boxes
      use chgpot
      use energi
      use ewald
      use math
      use mpole
      implicit none
      integer i,ii
      real*8 e,f
      real*8 term,fterm
      real*8 cii,dii,qii
      real*8 xd,yd,zd
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
c
c
c     zero out the total atomic multipole energy
c
      em = 0.0d0
      if (npole .eq. 0)  return
c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     compute the real space portion of the Ewald summation
c
      call emreal0d
c
c     compute the reciprocal space part of the Ewald summation
c
      call emrecip
c
c     compute the self-energy portion of the Ewald summation
c
      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / sqrtpi
      do i = 1, npole
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         cii = ci*ci
         dii = dix*dix + diy*diy + diz*diz
         qii = 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
     &            + qixx*qixx + qiyy*qiyy + qizz*qizz
         e = fterm * (cii + term*(dii/3.0d0+2.0d0*term*qii/5.0d0))
         em = em + e
      end do
c
c     compute the cell dipole boundary correction term
c
      if (boundary .eq. 'VACUUM') then
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            dix = rpole(2,i)
            diy = rpole(3,i)
            diz = rpole(4,i)
            xd = xd + dix + rpole(1,i)*x(ii)
            yd = yd + diy + rpole(1,i)*y(ii)
            zd = zd + diz + rpole(1,i)*z(ii)
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         e = term * (xd*xd+yd*yd+zd*zd)
         em = em + e
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine emreal0d  --  real space mpole energy via list  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "emreal0d" evaluates the real space portion of the Ewald sum
c     energy due to atomic multipoles using a neighbor list
c
c     literature reference:
c
c     W. Smith, "Point Multipoles in the Ewald Summation (Revisited)",
c     CCP5 Newsletter, 46, 18-30, 1998  (see http://www.ccp5.org/)
c
c
      subroutine emreal0d
      use sizes
      use atomid
      use atoms
      use bound
      use polar 
      use chgpot
      use couple
      use energi
      use ewald
      use math
      use mplpot
      use mpole
      use neigh
      use shunt
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      real*8 e,f,bfac,erfc
      real*8 alsq2,alsq2n
      real*8 exp2a,ralpha
      real*8 scalekk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 qrix,qriy,qriz
      real*8 qrkx,qrky,qrkz
      real*8 dri,drk,dik
      real*8 qrri,qrrk
      real*8 qrrik,qik
      real*8 diqrk,dkqri
      real*8 term1,term2,term3
      real*8 term4,term5
      real*8 bn(0:4)
!     charge penetration varables
      real*8 alphai,alphak
      real*8 dampi,dampk
      real*8 expdampi,expdampk
      real*8 termi,termk
      real*8, allocatable ::  scalei(:),scalek(:)
      real*8, allocatable ::  scaleik(:)
      real*8 nuci,qi,nuck,qk
!     end 
      real*8, allocatable :: mscale(:)
      character*6 mode
      external erfc
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         mscale(i) = 1.0d0
      end do
c
c     maximum order to rr9 
c
      allocate (scalei(9))
      allocate (scalek(9))
      allocate (scaleik(9))
c
c     maximum order to rr9
c
      do i = 1, 9
         scalei(i) = 1.0d0
         scalek(i) = 1.0d0
         scaleik(i) = 1.0d0
      end do

c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private)
!$OMP& shared(npole,ipole,x,y,z,rpole,n12,i12,n13,i13,
!$OMP& n14,i14,n15,i15,m2scale,m3scale,m4scale,m5scale,
!$OMP& atomic,penalpha,type,
!$OMP& nelst,elst,use_bounds,f,off2,aewald)
!$OMP& firstprivate(mscale) shared (em)
!$OMP DO reduction(+:em) schedule(guided)
c
c     compute the real space portion of the Ewald summation
c
      do i = 1, npole
         ii = ipole(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
         end do
c
c     evaluate all sites within the cutoff distance
c
         do kkk = 1, nelst(i)
            k = elst(kkk,i)
            kk = ipole(k)
            xr = x(kk) - xi
            yr = y(kk) - yi
            zr = z(kk) - zi
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
               ck = rpole(1,k)
               dkx = rpole(2,k)
               dky = rpole(3,k)
               dkz = rpole(4,k)
               qkxx = rpole(5,k)
               qkxy = rpole(6,k)
               qkxz = rpole(7,k)
               qkyy = rpole(9,k)
               qkyz = rpole(10,k)
               qkzz = rpole(13,k)
c
c     get reciprocal distance terms for this interaction
c
               rr1 = f / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
c
c     calculate the real space Ewald error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 4
                  bfac = dble(j+j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
               do j = 0, 4
                  bn(j) = f * bn(j)
               end do
c
c     intermediates involving moments and distance separation
c
               dri = dix*xr + diy*yr + diz*zr
               drk = dkx*xr + dky*yr + dkz*zr
               dik = dix*dkx + diy*dky + diz*dkz
               qrix = qixx*xr + qixy*yr + qixz*zr
               qriy = qixy*xr + qiyy*yr + qiyz*zr
               qriz = qixz*xr + qiyz*yr + qizz*zr
               qrkx = qkxx*xr + qkxy*yr + qkxz*zr
               qrky = qkxy*xr + qkyy*yr + qkyz*zr
               qrkz = qkxz*xr + qkyz*yr + qkzz*zr
               qrri = qrix*xr + qriy*yr + qriz*zr
               qrrk = qrkx*xr + qrky*yr + qrkz*zr
               qrrik = qrix*qrkx + qriy*qrky + qriz*qrkz
               qik = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                  + qixx*qkxx + qiyy*qkyy + qizz*qkzz
               diqrk = dix*qrkx + diy*qrky + diz*qrkz
               dkqri = dkx*qrix + dky*qriy + dkz*qriz
c
c     read in charge penetration damping parameters
c
               alphai = penalpha(type(ii))
               alphak = penalpha(type(kk))
c
c     compute common factors for damping
c
               dampi = alphai*r
               dampk = alphak*r
               expdampi = exp(-dampi)
               expdampk = exp(-dampk)
c
c     calculate one-site scale factors
c
               scalei(1) = 1.0d0 - expdampi
               scalek(1) = 1.0d0 - expdampk
               scalei(3) = 1.0d0 - (1.0d0 + dampi)*expdampi
               scalek(3) = 1.0d0 - (1.0d0 +dampk)*expdampk
               scalei(5) = 1.0d0 - (1.0d0 + dampi + 
     &              (1.0d0/3.0d0)*dampi**2)*expdampi
               scalek(5) = 1.0d0-(1.0d0 +dampk +
     &              (1.0d0/3.0d0)*dampk**2)*expdampk
               scalei(7) = 1.0d0-(1.0d0 + dampi + 0.4d0*dampi**2 +
     &              (1.0d0/15.0d0)*dampi**3)*expdampi
               scalek(7) = 1.0d0-(1.0d0 + dampk + 0.4d0*dampk**2 +
     &              (1.0d0/15.0d0)*dampk**3)*expdampk
c
c     calculate two-site scale factors
c
               if (alphai .ne. alphak) then
                  termi = alphak**2/(alphak**2 - alphai**2)
                  termk = alphai**2/(alphai**2 - alphak**2)
                  scaleik(1) =1.0d0-termi*expdampi -termk*expdampk
                  scaleik(3) =1.0d0-termi*(1.0d0 +dampi)*expdampi
     &                        - termk*(1.0d0 + dampk)*expdampk
                  scaleik(5) = 1.0d0 - termi*(1.0d0 + dampi +
     &                 (1.0d0/3.0d0)*dampi**2)*expdampi -
     &                 termk*(1.0d0 + dampk +
     &                 (1.0d0/3.0d0)*dampk**2)*expdampk
                  scaleik(7) = 1.0d0 - termi*(1.0d0 + dampi +
     &                 0.4d0*dampi**2 + (1.0d0/15.0d0)*dampi**3)*
     &                 expdampi -
     &                 termk*(1.0d0 + dampk +
     &                 0.4d0*dampk**2 + (1.0d0/15.0d0)*dampk**3)*
     &                 expdampk
                  scaleik(9) = 1.0d0 - termi*(1.0d0 + dampi +
     &                 (3.0d0/7.0d0)*dampi**2 +
     &                 (2.0d0/21.0d0)*dampi**3 +
     &                 (1.0d0/105.0d0)*dampi**4)*expdampi -
     &                 termk*(1.0d0 + dampk +
     &                 (3.0d0/7.0d0)*dampk**2 +
     &                 (2.0d0/21.0d0)*dampk**3 +
     &                 (1.0d0/105.0d0)*dampk**4)*expdampk
               else
                  scaleik(1) = 1.0d0 -(1.0d0+0.5d0*dampi)*expdampi
                  scaleik(3) = 1.0d0 -(1.0d0+dampi+0.5d0*dampi**2)
     &                 *expdampi
                  scaleik(5) = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2
     &                 + (1.0d0/6.0d0)*dampi**3)*expdampi
                  scaleik(7) = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2
     &                 + (1.0d0/6.0d0)*dampi**3
     &                 + (1.0d0/30.0d0)*dampi**4)*expdampi
                  scaleik(9) = 1.0d0 - (1.0d0+dampi+0.5d0*dampi**2
     &                 + (1.0d0/6.0d0)*dampi**3
     &                 + (4.0d0/105.0d0)*dampi**4
     &                 + (1.0d0/210.0d0)*dampi**5)*expdampi
               end if
            
               nuci = atomic(i)
               nuck = atomic(k)
               qi = ci - nuci
               qk = ck - nuck

               term1 = nuci*nuck + nuci*qk*scalek(1) 
     &                 + nuck*qi*scalei(1) + qi*qk*scaleik(1)
               term2 = nuck*dri*scalei(3) - nuci*drk*scalek(3)
     &                 + (dik + qk*dri - qi*drk)*scaleik(3)
               term3 = nuci*qrrk*scalek(5) + nuck*qrri*scalei(5)+
     &                 (-dri*drk + 2.0d0*(dkqri-diqrk+qik) + 
     &                 qi*qrrk + qk*qrri)*scaleik(5) 
               term4 = (dri*qrrk-drk*qrri-4.0d0*qrrik)*scaleik(7)
               term5 = qrri*qrrk*scaleik(9)
c
c     modify distances to account for Ewald and exclusions
c
               bn(0) = bn(0)/rr1
               bn(1) = bn(1)/rr3
               bn(2) = bn(2)/rr5
               bn(3) = bn(3)/rr7
               bn(4) = bn(4)/rr9

               term1 = nuci*nuck*(bn(0)-(1.0d0-mscale(kk))) 
     &             + nuci*qk*(bn(0)-(1.0d0-scalek(1)*mscale(kk))) 
     &             + nuck*qi*(bn(0)-(1.0d0-scalei(1)*mscale(kk))) 
     &             + qi*qk*(bn(0)-(1.0d0-scaleik(1)*mscale(kk)))

               term2 = nuck*dri*(bn(1)-(1.0d0-scalei(3)*mscale(kk))) 
     &             - nuci*drk*(bn(1)-(1.0d0-scalek(3)*mscale(kk)))
     &             + (dik + qk*dri - qi*drk)
     &               *(bn(1)-(1.0d0-scaleik(3)*mscale(kk)))
               term3 = nuci*qrrk*(bn(2)-(1.0d0-scalek(5)*mscale(kk)))
     &             + nuck*qrri*(bn(2)-(1.0d0-scalei(5)*mscale(kk)))
     &             +(-dri*drk + 2.0d0*(dkqri-diqrk+qik) 
     &             + qi*qrrk + qk*qrri)
     &             * (bn(2)-(1.0d0-scaleik(5)*mscale(kk))) 

               term4 = (dri*qrrk-drk*qrri-4.0d0*qrrik)
     &             *(bn(3)-(1.0d0-scaleik(7)*mscale(kk)))

               term5=qrri*qrrk*(bn(4)-(1.0d0-scaleik(9)*mscale(kk)))

c
c     compute the energy contribution for this interaction
c
               e = term1*rr1 + term2*rr3 + term3*rr5
     &                + term4*rr7 + term5*rr9
               em = em + e
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine emrecip  --  PME recip space multipole energy  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "emrecip" evaluates the reciprocal space portion of the particle
c     mesh Ewald energy due to atomic multipole interactions
c
c     literature reference:
c
c     C. Sagui, L. G. Pedersen and T. A. Darden, "Towards an Accurate
c     Representation of Electrostatics in Classical Force Fields:
c     Efficient Implementation of Multipolar Interactions in
c     Biomolecular Simulations", Journal of Chemical Physics, 120,
c     73-87 (2004)
c
c     modifications for nonperiodic systems suggested by Tom Darden
c     during May 2007
c
c
      subroutine emrecip
      use sizes
      use bound
      use boxes
      use chgpot
      use energi
      use ewald
      use math
      use mpole
      use mrecip
      use pme
      use potent
      implicit none
      integer i,j,k
      integer k1,k2,k3
      integer m1,m2,m3
      integer ntot,nff
      integer nf1,nf2,nf3
      real*8 e,r1,r2,r3
      real*8 f,h1,h2,h3
      real*8 volterm,denom
      real*8 hsq,expterm
      real*8 term,pterm
      real*8 struc2
c
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
      f = electric / dielec
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(cmp)) then
         if (size(cmp) .lt. 10*npole) then
            deallocate (cmp)
            deallocate (fmp)
            deallocate (cphi)
            deallocate (fphi)
         end if
      end if
      if (.not. allocated(cmp)) then
         allocate (cmp(10,npole))
         allocate (fmp(10,npole))
         allocate (cphi(10,npole))
         allocate (fphi(20,npole))
      end if
c
c     copy the multipole moments into local storage areas
c
      do i = 1, npole
         cmp(1,i) = rpole(1,i)
         cmp(2,i) = rpole(2,i)
         cmp(3,i) = rpole(3,i)
         cmp(4,i) = rpole(4,i)
         cmp(5,i) = rpole(5,i)
         cmp(6,i) = rpole(9,i)
         cmp(7,i) = rpole(13,i)
         cmp(8,i) = 2.0d0 * rpole(6,i)
         cmp(9,i) = 2.0d0 * rpole(7,i)
         cmp(10,i) = 2.0d0 * rpole(10,i)
      end do
c
c     compute B-spline coefficients and spatial decomposition
c
      call bspline_fill
      call table_fill
c
c     convert Cartesian multipoles to fractional coordinates
c
      call cmp_to_fmp (cmp,fmp)
c
c     assign PME grid and perform 3-D FFT forward transform
c
      call grid_mpole (fmp)
      call fftfront
c
c     make the scalar summation over reciprocal lattice
c
      ntot = nfft1 * nfft2 * nfft3
      pterm = (pi/aewald)**2
      volterm = pi * volbox
      nff = nfft1 * nfft2
      nf1 = (nfft1+1) / 2
      nf2 = (nfft2+1) / 2
      nf3 = (nfft3+1) / 2
      do i = 1, ntot-1
         k3 = i/nff + 1
         j = i - (k3-1)*nff
         k2 = j/nfft1 + 1
         k1 = j - (k2-1)*nfft1 + 1
         m1 = k1 - 1
         m2 = k2 - 1
         m3 = k3 - 1
         if (k1 .gt. nf1)  m1 = m1 - nfft1
         if (k2 .gt. nf2)  m2 = m2 - nfft2
         if (k3 .gt. nf3)  m3 = m3 - nfft3
         r1 = dble(m1)
         r2 = dble(m2)
         r3 = dble(m3)
         h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
         h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
         h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
         hsq = h1*h1 + h2*h2 + h3*h3
         term = -pterm * hsq
         expterm = 0.0d0
         if (term .gt. -50.0d0) then
            denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
            expterm = exp(term) / denom
            if (.not. use_bounds) then
               expterm = expterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
            else if (octahedron) then
               if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
            end if
         end if
         qfac(k1,k2,k3) = expterm
      end do
c
c     account for zeroth grid point for nonperiodic system
c
      qfac(1,1,1) = 0.0d0
      if (.not. use_bounds) then
         expterm = 0.5d0 * pi / xbox
         struc2 = qgrid(1,1,1,1)**2 + qgrid(2,1,1,1)**2
         e = f * expterm * struc2
         em = em + e
         qfac(1,1,1) = expterm
      end if
c
c     complete the transformation of the charge grid
c
      do k = 1, nfft3
         do j = 1, nfft2
            do i = 1, nfft1
               term = qfac(i,j,k)
               qgrid(1,i,j,k) = term * qgrid(1,i,j,k)
               qgrid(2,i,j,k) = term * qgrid(2,i,j,k)
            end do
         end do
      end do
c
c     perform 3-D FFT backward transform and get potential
c
      call fftback
      call fphi_mpole (fphi)
c
c     sum over multipoles and increment total multipole energy
c
      e = 0.0d0
      do i = 1, npole
         do k = 1, 10
            e = e + fmp(k,i)*fphi(k,i)
         end do
      end do
      e = 0.5d0 * f * e
      em = em + e
      return
      end
