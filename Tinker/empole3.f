c
c
c     #############################################################
c     ##  COPYRIGHT (C) 1999 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine empole3  --  atomic multipole energy & analysis  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "empole3" calculates the electrostatic energy due to atomic
c     multipole interactions, and partitions the energy among atoms
c
c
      subroutine empole3
      use limits
      implicit none
c
c
c     choose the method for summing over multipole interactions
c
      if (use_ewald) then
         if (use_mlist) then
            call empole3d
         else
            call empole3c
         end if
      else
         if (use_mlist) then
            call empole3b
         else
            call empole3a
         end if
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine empole3a  --  double loop multipole analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "empole3a" calculates the atomic multipole interaction energy
c     using a double loop, and partitions the energy among atoms
c
c
      subroutine empole3a
      use sizes
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use cell
      use polar 
      use chgpot
      use couple
      use energi
      use group
      use inform
      use inter
      use iounit
      use math
      use molcul
      use mplpot
      use mpole
      use potent
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
      logical proceed
      logical header,huge
      logical usei,usek
      character*6 mode
c
c
c     zero out total atomic multipole energy and partitioning
c
      nem = 0
      em = 0.0d0
      do i = 1, n
         aem(i) = 0.0d0
      end do
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
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. npole.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Atomic Multipole Interactions :',
     &           //,' Type',14x,'Atom Names',15x,'Distance',
     &              8x,'Energy',/)
      end if
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
c
c     increment the overall multipole energy components
c
                  if (e .ne. 0.0d0) then
                     nem = nem + 1
                     em = em + e
                     aem(ii) = aem(ii) + 0.5d0*e
                     aem(kk) = aem(kk) + 0.5d0*e
                     if (molcule(ii) .ne. molcule(kk))
     &                  einter = einter + e
                  end if
c
c     print message if the energy of this interaction is large
c
                  huge = (abs(e) .gt. 100.0d0)
                  if ((debug.and.e.ne.0.0d0)
     &                  .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,20)
   20                   format (/,' Individual Atomic Multipole',
     &                             ' Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             15x,'Distance',8x,'Energy',/)
                     end if
                     write (iout,30)  ii,name(ii),kk,name(kk),r,e
   30                format (' M-Pole',4x,2(i7,'-',a3),9x,
     &                          f10.4,2x,f12.4)
                  end if
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
c
c     increment the overall multipole energy components
c
                        if (e .ne. 0.0d0) then
                           nem = nem + 1
                           em = em + e
                           aem(ii) = aem(ii) + 0.5d0*e
                           aem(kk) = aem(kk) + 0.5d0*e
                           einter = einter + e
                        end if
c
c     print message if the energy of this interaction is large
c
                        huge = (abs(e) .gt. 100.0d0)
                        if ((debug.and.e.ne.0.0d0)
     &                        .or. (verbose.and.huge)) then
                           if (header) then
                              header = .false.
                              write (iout,40)
   40                         format (/,' Individual Atomic Multipole',
     &                                   ' Interactions :',
     &                                //,' Type',14x,'Atom Names',
     &                                   15x,'Distance',8x,'Energy',/)
                           end if
                           write (iout,50)  ii,name(ii),kk,name(kk),r,e
   50                      format (' M-Pole',4x,2(i7,'-',a3),1x,
     &                                '(XTAL)',2x,f10.4,2x,f12.4)
                        end if
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
c     #################################################################
c     ##                                                             ##
c     ##  subroutine empole3b  --  neighbor list multipole analysis  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "empole3b" calculates the atomic multipole interaction energy
c     using a neighbor list, and partitions the energy among the atoms
c
c
      subroutine empole3b
      use sizes
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use polar 
      use chgpot
      use couple
      use energi
      use group
      use inform
      use inter
      use iounit
      use math
      use molcul
      use mplpot
      use mpole
      use neigh
      use potent
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
      logical proceed
      logical header,huge
      logical usei,usek
      character*6 mode
c
c
c     zero out total atomic multipole energy and partitioning
c
      nem = 0
      em = 0.0d0
      do i = 1, n
         aem(i) = 0.0d0
      end do
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. npole.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Atomic Multipole Interactions :',
     &           //,' Type',14x,'Atom Names',15x,'Distance',
     &              8x,'Energy',/)
      end if
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
!$OMP& m5scale,nelst,elst,use_group,use_intra,use_bounds,off2,
!$OMP& f,molcule,name,verbose,debug,header,iout)
!$OMP& firstprivate(mscale) shared (em,einter,nem,aem)
!$OMP DO reduction(+:em,einter,nem,aem) schedule(guided)
c
c     calculate the multipole interaction energy term
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
c
c     increment the overall multipole energy components
c
                  if (e .ne. 0.0d0) then
                     nem = nem + 1
                     em = em + e
                     aem(ii) = aem(ii) + 0.5d0*e
                     aem(kk) = aem(kk) + 0.5d0*e
                     if (molcule(ii) .ne. molcule(kk))
     &                  einter = einter + e
                  end if
c
c     print message if the energy of this interaction is large
c
                  huge = (abs(e) .gt. 100.0d0)
                  if ((debug.and.e.ne.0.0d0)
     &                  .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,20)
   20                   format (/,' Individual Atomic Multipole',
     &                             ' Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             15x,'Distance',8x,'Energy',/)
                     end if
                     write (iout,30)  ii,name(ii),kk,name(kk),r,e
   30                format (' M-Pole',4x,2(i7,'-',a3),9x,
     &                          f10.4,2x,f12.4)
                  end if
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
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine empole3c  --  Ewald multipole analysis via loop  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "empole3c" calculates the atomic multipole interaction energy
c     using a particle mesh Ewald summation and double loop, and
c     partitions the energy among the atoms
c
c
      subroutine empole3c
      use sizes
      use action
      use analyz
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
c     zero out the multipole and polarization energies
c
      nem = 0
      em = 0.0d0
      do i = 1, n
         aem(i) = 0.0d0
      end do
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
c     compute the real space part of the Ewald summation
c
      call emreal3c
c
c     compute the reciprocal space part of the Ewald summation
c
      call emrecip
c
c     compute the self-energy part of the Ewald summation
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
         nem = nem + 1
         aem(i) = aem(i) + e
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
         nem = nem + 1
         do i = 1, npole
            ii = ipole(i)
            aem(ii) = aem(ii) + e/dble(npole)
         end do
      end if
      return
      end
c
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine emreal3c  --  real space mpole analysis via loop  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "emreal3c" evaluates the real space portion of the Ewald sum
c     energy due to atomic multipole interactions and partitions
c     the energy among the atoms
c
c     literature reference:
c
c     W. Smith, "Point Multipoles in the Ewald Summation (Revisited)",
c     CCP5 Newsletter, 46, 18-30, 1998  (see http://www.ccp5.org/)
c
c
      subroutine emreal3c
      use sizes
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use cell
      use polar 
      use chgpot
      use couple
      use energi
      use ewald
      use inform
      use inter
      use iounit
      use math
      use molcul
      use mplpot
      use mpole
      use potent
      use shunt
      implicit none
      integer i,j,k
      integer ii,kk
      integer jcell
      real*8 e,efull,f
      real*8 bfac,erfc
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
      logical header,huge
      character*6 mode
      external erfc
c
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. npole.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Atomic Multipole Interactions :',
     &           //,' Type',14x,'Atom Names',15x,'Distance',
     &              8x,'Energy',/)
      end if
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
c     compute the full undamped energy for this interaction
c
               efull = term1*rr1 + term2*rr3 + term3*rr5
     &                    + term4*rr7 + term5*rr9
               efull = mscale(kk) * efull
               if (efull .ne. 0.0d0) then
                  nem = nem + 1
                  aem(ii) = aem(ii) + 0.5d0*efull
                  aem(kk) = aem(kk) + 0.5d0*efull
                  if (molcule(ii) .ne. molcule(kk))
     &               einter = einter + efull
               end if
c
c     modify error function terms to account for scaling
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
c
c     print a message if the energy of this interaction is large
c
               huge = (abs(efull) .gt. 100.0d0)
               if ((debug.and.efull.ne.0.0d0)
     &               .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,20)
   20                format (/,' Individual Atomic Multipole',
     &                          ' Interactions :',
     &                       //,' Type',14x,'Atom Names',
     &                          15x,'Distance',8x,'Energy',/)
                  end if
                  write (iout,30)  ii,name(ii),kk,name(kk),r,efull
   30             format (' M-Pole',4x,2(i7,'-',a3),9x,
     &                       f10.4,2x,f12.4)
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

                     term1 = nuci*nuck + nuci*qk*scalek(1) 
     &                       + nuck*qi*scalei(1) + qi*qk*scaleik(1)
                     term2 = nuck*dri*scalei(3) - nuci*drk*scalek(3)
     &                       + (dik + qk*dri - qi*drk)*scaleik(3)
                     term3 = nuci*qrrk*scalek(5) + nuck*qrri*scalei(5)+
     &                       (-dri*drk + 2.0d0*(dkqri-diqrk+qik) + 
     &                       qi*qrrk + qk*qrri)*scaleik(5) 
                     term4 = (dri*qrrk-drk*qrri-4.0d0*qrrik)*scaleik(7)
                     term5 = qrri*qrrk*scaleik(9)
c
c     compute the full undamped energy for this interaction
c
                     efull = term1*rr1 + term2*rr3 + term3*rr5
     &                          + term4*rr7 + term5*rr9
                     efull = mscale(kk) * efull
                     if (ii .eq. kk)  efull = 0.5d0 * efull
                     if (efull .ne. 0.0d0) then
                        nem = nem + 1
                        aem(ii) = aem(ii) + 0.5d0*efull
                        aem(kk) = aem(kk) + 0.5d0*efull
                        einter = einter + efull
                     end if
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
     &                     *(bn(1)-(1.0d0-scaleik(3)*mscale(kk)))
                     term3 = nuci*qrrk*(bn(2)-(1.0d0-scalek(5)
     &                   * mscale(kk)))
     &                   + nuck*qrri*(bn(2)-(1.0d0-scalei(5)
     &                   * mscale(kk)))
     &                   +(-dri*drk + 2.0d0*(dkqri-diqrk+qik) 
     &                   + qi*qrrk + qk*qrri)
     &                   * (bn(2)-(1.0d0-scaleik(5)*mscale(kk))) 
                     term4 = (dri*qrrk-drk*qrri-4.0d0*qrrik)
     &                   *(bn(3)-(1.0d0-scaleik(7)*mscale(kk)))

                     term5=qrri*qrrk*(bn(4)-(1.0d0-scaleik(9)
     &                   * mscale(kk)))

c
c     compute the energy contribution for this interaction
c
                     e = term1*rr1 + term2*rr3 + term3*rr5
     &                      + term4*rr7 + term5*rr9
                     if (ii .eq. kk)  e = 0.5d0 * e
                     em = em + e
c
c     print message if the energy of this interaction is large
c
                     huge = (abs(efull) .gt. 100.0d0)
                     if ((debug.and.efull.ne.0.0d0)
     &                     .or. (verbose.and.huge)) then
                        if (header) then
                           header = .false.
                           write (iout,40)
   40                      format (/,' Individual Atomic Multipole',
     &                                ' Interactions :',
     &                             //,' Type',14x,'Atom Names',
     &                                15x,'Distance',8x,'Energy',/)
                        end if
                        write (iout,50)  ii,name(ii),kk,name(kk),r,efull
   50                   format (' M-Pole',4x,2(i7,'-',a3),1x,
     &                             '(XTAL)',2x,f10.4,2x,f12.4)
                     end if
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
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine empole3d  --  Ewald multipole analysis via list  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "empole3d" calculates the atomic multipole interaction energy
c     using particle mesh Ewald summation and a neighbor list, and
c     partitions the energy among the atoms
c
c
      subroutine empole3d
      use sizes
      use action
      use analyz
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
c     zero out the multipole and polarization energies
c
      nem = 0
      em = 0.0d0
      do i = 1, n
         aem(i) = 0.0d0
      end do
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
c     compute the real space part of the Ewald summation
c
      call emreal3d
c
c     compute the reciprocal space part of the Ewald summation
c
      call emrecip
c
c     compute the self-energy part of the Ewald summation
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
         nem = nem + 1
         aem(i) = aem(i) + e
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
         nem = nem + 1
         do i = 1, npole
            ii = ipole(i)
            aem(ii) = aem(ii) + e/dble(npole)
         end do
      end if
      return
      end
c
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine emreal3d  --  real space mpole analysis via list  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "emreal3d" evaluates the real space portion of the Ewald sum
c     energy due to atomic multipole interactions, and partitions
c     the energy among the atoms using a pairwise neighbor list
c
c     literature reference:
c
c     W. Smith, "Point Multipoles in the Ewald Summation (Revisited)",
c     CCP5 Newsletter, 46, 18-30, 1998  (see http://www.ccp5.org/)
c
c
      subroutine emreal3d
      use sizes
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use polar 
      use chgpot
      use couple
      use energi
      use ewald
      use inform
      use inter
      use iounit
      use math
      use molcul
      use mplpot
      use mpole
      use neigh
      use potent
      use shunt
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      real*8 e,efull,f
      real*8 bfac,erfc
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
      logical header,huge
      character*6 mode
      external erfc
c
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. npole.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Atomic Multipole Interactions :',
     &           //,' Type',14x,'Atom Names',15x,'Distance',
     &              8x,'Energy',/)
      end if
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
!$OMP& shared(npole,ipole,x,y,z,rpole,n12,i12,n13,i13,n14,i14,n15,
!$OMP& i15,m2scale,m3scale,m4scale,m5scale,nelst,elst,use_bounds,
!$OMP& atomic,penalpha,type,
!$OMP& f,off2,aewald,molcule,name,verbose,debug,header,iout)
!$OMP& firstprivate(mscale) shared (em,einter,nem,aem)
!$OMP DO reduction(+:em,einter,nem,aem) schedule(guided)
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
c     compute the full undamped energy for this interaction
c
               efull = term1*rr1 + term2*rr3 + term3*rr5
     &                    + term4*rr7 + term5*rr9
               efull = mscale(kk) * efull
               if (efull .ne. 0.0d0) then
                  nem = nem + 1
                  aem(ii) = aem(ii) + 0.5d0*efull
                  aem(kk) = aem(kk) + 0.5d0*efull
                  if (molcule(ii) .ne. molcule(kk))
     &               einter = einter + efull
               end if
c
c     modify error function terms to account for scaling
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
c
c     print a message if the energy of this interaction is large
c
               huge = (abs(efull) .gt. 100.0d0)
               if ((debug.and.efull.ne.0.0d0)
     &               .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,20)
   20                format (/,' Individual Atomic Multipole',
     &                          ' Interactions :',
     &                       //,' Type',14x,'Atom Names',
     &                          15x,'Distance',8x,'Energy',/)
                  end if
                  write (iout,30)  ii,name(ii),kk,name(kk),r,efull
   30             format (' M-Pole',4x,2(i7,'-',a3),9x,
     &                       f10.4,2x,f12.4)
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
