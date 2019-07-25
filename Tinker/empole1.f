c
c
c     #############################################################
c     ##  COPYRIGHT (C) 1999 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine empole1  --  mpole/polar energy & derivatives  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "empole1" calculates the multipole and dipole polarization
c     energy and derivatives with respect to Cartesian coordinates
c
c
      subroutine empole1
      use limits
      implicit none
c
c
c     choose the method for summing over multipole interactions
c
      if (use_ewald) then
         if (use_mlist) then
            call empole1d
         else
            call empole1c
         end if
      else
         if (use_mlist) then
            call empole1b
         else
            call empole1a
         end if
      end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine empole1a  --  double loop multipole derivatives  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "empole1a" calculates the multipole energy and derivatives with 
c     respect to Cartesian coordinates using a pairwise double loop
c
c
      subroutine empole1a
      use sizes
      use atomid 
      use atoms
      use bound
      use cell
      use polar 
      use chgpot
      use couple
      use deriv
      use energi
      use group
      use inter
      use molcul
      use mplpot
      use mpole
      use shunt
      use usage
      use virial
      implicit none
      integer i,j,k
      integer ii,kk,jcell
      integer ix,iy,iz
      integer kx,ky,kz
      integer iax,iay,iaz
      real*8 e,de,f,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9,rr11
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dikx,diky,dikz
      real*8 dirx,diry,dirz
      real*8 dkrx,dkry,dkrz
      real*8 qrix,qriy,qriz
      real*8 qrkx,qrky,qrkz
      real*8 qrixr,qriyr,qrizr
      real*8 qrkxr,qrkyr,qrkzr
      real*8 qrrx,qrry,qrrz
      real*8 qikrx,qikry,qikrz
      real*8 qkirx,qkiry,qkirz
      real*8 qikrxr,qikryr,qikrzr
      real*8 qkirxr,qkiryr,qkirzr
      real*8 diqkx,diqky,diqkz
      real*8 dkqix,dkqiy,dkqiz
      real*8 diqkxr,diqkyr,diqkzr
      real*8 dkqixr,dkqiyr,dkqizr
      real*8 dqiqkx,dqiqky,dqiqkz
      real*8 dri,drk,qrri,qrrk
      real*8 diqrk,dkqri
      real*8 dik,qik,qrrik
      real*8 term1,term2,term3
      real*8 term4,term5,term6
      real*8 frcx,frcy,frcz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 ttmi(3),ttmk(3)
      real*8 fix(3),fiy(3),fiz(3)
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
      real*8, allocatable :: tem(:,:)
      logical proceed,usei,usek
      character*6 mode
c
c
c     zero out the atomic multipole energy and derivatives
c
      em = 0.0d0
      do i = 1, n
         do j = 1, 3
            dem(j,i) = 0.0d0
         end do
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
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (tem(3,n))
c
c     maximum order to rr11 (for force)
c
      allocate (scalei(11))
      allocate (scalek(11))
      allocate (scaleik(11))
c
c     maximum order to rr11 (for force)
c
      do i = 1, 11 
         scalei(i) = 1.0d0
         scalek(i) = 1.0d0
         scaleik(i) = 1.0d0
      end do
c
c     initialize connected atom scaling and torque arrays
c
      do i = 1, n
         mscale(i) = 1.0d0
         do j = 1, 3
            tem(j,i) = 0.0d0
         end do
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'MPOLE'
      call switch (mode)
c
c     compute the multipole interaction energy and gradient
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
            if (.not. proceed)  goto 10
            xr = x(kk) - xi
            yr = y(kk) - yi
            zr = z(kk) - zi
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
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
               rr11 = 9.0d0 * rr9 / r2
c
c     intermediates involving moments and distance separation
c
               dikx = diy*dkz - diz*dky
               diky = diz*dkx - dix*dkz
               dikz = dix*dky - diy*dkx
               dirx = diy*zr - diz*yr
               diry = diz*xr - dix*zr
               dirz = dix*yr - diy*xr
               dkrx = dky*zr - dkz*yr
               dkry = dkz*xr - dkx*zr
               dkrz = dkx*yr - dky*xr
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
               qrixr = qriz*yr - qriy*zr
               qriyr = qrix*zr - qriz*xr
               qrizr = qriy*xr - qrix*yr
               qrkxr = qrkz*yr - qrky*zr
               qrkyr = qrkx*zr - qrkz*xr
               qrkzr = qrky*xr - qrkx*yr
               qrrx = qrky*qriz - qrkz*qriy
               qrry = qrkz*qrix - qrkx*qriz
               qrrz = qrkx*qriy - qrky*qrix
               qikrx = qixx*qrkx + qixy*qrky + qixz*qrkz
               qikry = qixy*qrkx + qiyy*qrky + qiyz*qrkz
               qikrz = qixz*qrkx + qiyz*qrky + qizz*qrkz
               qkirx = qkxx*qrix + qkxy*qriy + qkxz*qriz
               qkiry = qkxy*qrix + qkyy*qriy + qkyz*qriz
               qkirz = qkxz*qrix + qkyz*qriy + qkzz*qriz
               qikrxr = qikrz*yr - qikry*zr
               qikryr = qikrx*zr - qikrz*xr
               qikrzr = qikry*xr - qikrx*yr
               qkirxr = qkirz*yr - qkiry*zr
               qkiryr = qkirx*zr - qkirz*xr
               qkirzr = qkiry*xr - qkirx*yr
               diqkx = dix*qkxx + diy*qkxy + diz*qkxz
               diqky = dix*qkxy + diy*qkyy + diz*qkyz
               diqkz = dix*qkxz + diy*qkyz + diz*qkzz
               dkqix = dkx*qixx + dky*qixy + dkz*qixz
               dkqiy = dkx*qixy + dky*qiyy + dkz*qiyz
               dkqiz = dkx*qixz + dky*qiyz + dkz*qizz
               diqrk = dix*qrkx + diy*qrky + diz*qrkz
               dkqri = dkx*qrix + dky*qriy + dkz*qriz
               diqkxr = diqkz*yr - diqky*zr
               diqkyr = diqkx*zr - diqkz*xr
               diqkzr = diqky*xr - diqkx*yr
               dkqixr = dkqiz*yr - dkqiy*zr
               dkqiyr = dkqix*zr - dkqiz*xr
               dkqizr = dkqiy*xr - dkqix*yr
               dqiqkx = diy*qrkz - diz*qrky + dky*qriz - dkz*qriy
     &                     - 2.0d0*(qixy*qkxz+qiyy*qkyz+qiyz*qkzz
     &                             -qixz*qkxy-qiyz*qkyy-qizz*qkyz)
               dqiqky = diz*qrkx - dix*qrkz + dkz*qrix - dkx*qriz
     &                     - 2.0d0*(qixz*qkxx+qiyz*qkxy+qizz*qkxz
     &                             -qixx*qkxz-qixy*qkyz-qixz*qkzz)
               dqiqkz = dix*qrky - diy*qrkx + dkx*qriy - dky*qrix
     &                     - 2.0d0*(qixx*qkxy+qixy*qkyy+qixz*qkyz
     &                             -qixy*qkxx-qiyy*qkxy-qiyz*qkxz)
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
                  scaleik(11) = 1.0d0 - termi*(1.0d0 + dampi +
     &                 (4.0d0/9.0d0)*dampi**2 +
     &                 (1.0d0/9.0d0)*dampi**3 +
     &                 (1.0d0/63.0d0)*dampi**4 + 
     &                 (1.0d0/945.0d0)*dampi**5)*expdampi - 
     &                 termk*(1.0d0 + dampk +
     &                 (4.0d0/9.0d0)*dampk**2+
     &                 (1.0d0/9.0d0)*dampk**3+
     &                 (1.0d0/63.0d0)*dampk**4 +
     &                 (1.0d0/945.0d0)*dampk**5)*expdampk
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
                  scaleik(11) = 1.0d0-(1.0d0+dampi+0.5d0*dampi**2
     &                 + (1.0d0/6.0d0)*dampi**3 
     &                 + (5.0d0/126.0d0)*dampi**4
     &                 + (2.0d0/315.0d0)*dampi**5
     &                 + (1.0d0/1890.0d0)*dampi**6)*expdampi
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
c     compute the energy contribution for this interaction
c
               e = term1*rr1 + term2*rr3 + term3*rr5
     &                + term4*rr7 + term5*rr9
               em = em + e
               if (molcule(ii) .ne. molcule(kk))
     &            einter = einter + e
c
c     calculate intermediate terms for force and torque
c
               de =  nuci*nuck*rr3 
     &             + nuci * qk*rr3*scalek(3)
     &             + nuck * qi*rr3*scalei(3) 
     &             + qi * qk  *rr3*scaleik(3)
     &             + nuck*dri *rr5*scalei(5) 
     &             - nuci*drk *rr5*scalek(5)
     &             + (dik + qk*dri - qi*drk)*rr5*scaleik(5)
     &             + nuci*qrrk*rr7*scalek(7) 
     &             + nuck*qrri*rr7*scalei(7)
     &             +(-dri*drk + 2.0d0*(dkqri-diqrk+qik)  
     &             + qi*qrrk + qk*qrri)*rr7*scaleik(7) 
     &             +(dri*qrrk-drk*qrri-4.0d0*qrrik)*rr9*scaleik(9)
     &             + qrri*qrrk*rr11*scaleik(11)

               term1 = -nuck*rr3*scalei(3) 
     &                 - qk*rr3*scaleik(3) 
     &                 + drk*rr5*scaleik(5) 
     &                 - qrrk*rr7*scaleik(7)
               term2 = nuci*rr3*scalek(3) 
     &                 + qi*rr3*scaleik(3)  
     &                 + dri*rr5*scaleik(5) 
     &                 + qrri*rr7*scaleik(7)
               term3 = 2.0d0 * rr5 * scaleik(5) 
               term4 = 2.0d0 *(-nuck*rr5*scalei(5)
     &                       -qk*rr5*scaleik(5)
     &                       +drk*rr7*scaleik(7)
     &                       -qrrk*rr9*scaleik(9))
               term5 = 2.0d0 *(-nuci*rr5*scalek(5)
     &                       -qi*rr5*scaleik(5)
     &                       -dri*rr7*scaleik(7)
     &                       -qrri*rr9*scaleik(9))
               term6 = 4.0d0 * rr7 *scaleik(7) 

c
c     compute the force components for this interaction
c
               frcx = de*xr + term1*dix + term2*dkx
     &                   + term3*(diqkx-dkqix) + term4*qrix
     &                   + term5*qrkx + term6*(qikrx+qkirx)
               frcy = de*yr + term1*diy + term2*dky
     &                   + term3*(diqky-dkqiy) + term4*qriy
     &                   + term5*qrky + term6*(qikry+qkiry)
               frcz = de*zr + term1*diz + term2*dkz
     &                   + term3*(diqkz-dkqiz) + term4*qriz
     &                   + term5*qrkz + term6*(qikrz+qkirz)

c
c     compute the torque components for this interaction
c
               rr3 = rr3*scaleik(3)
               ttmi(1) = -rr3*dikx
     &                    + term1*dirx + term3*(dqiqkx+dkqixr)
     &                    - term4*qrixr - term6*(qikrxr+qrrx)
               ttmi(2) = -rr3*diky
     &                    + term1*diry + term3*(dqiqky+dkqiyr)
     &                    - term4*qriyr - term6*(qikryr+qrry)
               ttmi(3) = -rr3*dikz
     &                    + term1*dirz + term3*(dqiqkz+dkqizr)
     &                    - term4*qrizr - term6*(qikrzr+qrrz)
               ttmk(1) = rr3*dikx
     &                    + term2*dkrx - term3*(dqiqkx+diqkxr)
     &                    - term5*qrkxr - term6*(qkirxr-qrrx)
               ttmk(2) = rr3*diky
     &                    + term2*dkry - term3*(dqiqky+diqkyr)
     &                    - term5*qrkyr - term6*(qkiryr-qrry)
               ttmk(3) = rr3*dikz
     &                    + term2*dkrz - term3*(dqiqkz+diqkzr)
     &                    - term5*qrkzr - term6*(qkirzr-qrrz)

c
c     force and torque components scaled by group membership
c
               if (use_group) then
                  frcx = fgrp * frcx
                  frcy = fgrp * frcy
                  frcz = fgrp * frcz
                  do j = 1, 3
                     ttmi(j) = fgrp * ttmi(j)
                     ttmk(j) = fgrp * ttmk(j)
                  end do
               end if
c
c     increment force-based gradient and torque on first site
c
               dem(1,ii) = dem(1,ii) + frcx
               dem(2,ii) = dem(2,ii) + frcy
               dem(3,ii) = dem(3,ii) + frcz
               tem(1,i) = tem(1,i) + ttmi(1)
               tem(2,i) = tem(2,i) + ttmi(2)
               tem(3,i) = tem(3,i) + ttmi(3)
c
c     increment force-based gradient and torque on second site
c
               dem(1,kk) = dem(1,kk) - frcx
               dem(2,kk) = dem(2,kk) - frcy
               dem(3,kk) = dem(3,kk) - frcz
               tem(1,k) = tem(1,k) + ttmk(1)
               tem(2,k) = tem(2,k) + ttmk(2)
               tem(3,k) = tem(3,k) + ttmk(3)
c
c     increment the virial due to pairwise Cartesian forces
c
               vxx = -xr * frcx
               vxy = -yr * frcx
               vxz = -zr * frcx
               vyy = -yr * frcy
               vyz = -zr * frcy
               vzz = -zr * frcz
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vxy
               vir(3,1) = vir(3,1) + vxz
               vir(1,2) = vir(1,2) + vxy
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vyz
               vir(1,3) = vir(1,3) + vxz
               vir(2,3) = vir(2,3) + vyz
               vir(3,3) = vir(3,3) + vzz
            end if
   10       continue
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
c     calculate interaction with other unit cells
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
            if (.not. proceed)  goto 20
            do jcell = 1, ncell
            xr = x(kk) - xi
            yr = y(kk) - yi
            zr = z(kk) - zi
            call imager (xr,yr,zr,jcell)
            r2 = xr*xr + yr*yr + zr*zr
            if (.not. (use_polymer .and. r2.le.polycut2)) then
               mscale(kk) = 1.0d0
            end if
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
               rr11 = 9.0d0 * rr9 / r2
c
c     intermediates involving moments and distance separation
c
               dikx = diy*dkz - diz*dky
               diky = diz*dkx - dix*dkz
               dikz = dix*dky - diy*dkx
               dirx = diy*zr - diz*yr
               diry = diz*xr - dix*zr
               dirz = dix*yr - diy*xr
               dkrx = dky*zr - dkz*yr
               dkry = dkz*xr - dkx*zr
               dkrz = dkx*yr - dky*xr
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
               qrixr = qriz*yr - qriy*zr
               qriyr = qrix*zr - qriz*xr
               qrizr = qriy*xr - qrix*yr
               qrkxr = qrkz*yr - qrky*zr
               qrkyr = qrkx*zr - qrkz*xr
               qrkzr = qrky*xr - qrkx*yr
               qrrx = qrky*qriz - qrkz*qriy
               qrry = qrkz*qrix - qrkx*qriz
               qrrz = qrkx*qriy - qrky*qrix
               qikrx = qixx*qrkx + qixy*qrky + qixz*qrkz
               qikry = qixy*qrkx + qiyy*qrky + qiyz*qrkz
               qikrz = qixz*qrkx + qiyz*qrky + qizz*qrkz
               qkirx = qkxx*qrix + qkxy*qriy + qkxz*qriz
               qkiry = qkxy*qrix + qkyy*qriy + qkyz*qriz
               qkirz = qkxz*qrix + qkyz*qriy + qkzz*qriz
               qikrxr = qikrz*yr - qikry*zr
               qikryr = qikrx*zr - qikrz*xr
               qikrzr = qikry*xr - qikrx*yr
               qkirxr = qkirz*yr - qkiry*zr
               qkiryr = qkirx*zr - qkirz*xr
               qkirzr = qkiry*xr - qkirx*yr
               diqkx = dix*qkxx + diy*qkxy + diz*qkxz
               diqky = dix*qkxy + diy*qkyy + diz*qkyz
               diqkz = dix*qkxz + diy*qkyz + diz*qkzz
               dkqix = dkx*qixx + dky*qixy + dkz*qixz
               dkqiy = dkx*qixy + dky*qiyy + dkz*qiyz
               dkqiz = dkx*qixz + dky*qiyz + dkz*qizz
               diqrk = dix*qrkx + diy*qrky + diz*qrkz
               dkqri = dkx*qrix + dky*qriy + dkz*qriz
               diqkxr = diqkz*yr - diqky*zr
               diqkyr = diqkx*zr - diqkz*xr
               diqkzr = diqky*xr - diqkx*yr
               dkqixr = dkqiz*yr - dkqiy*zr
               dkqiyr = dkqix*zr - dkqiz*xr
               dkqizr = dkqiy*xr - dkqix*yr
               dqiqkx = diy*qrkz - diz*qrky + dky*qriz - dkz*qriy
     &                     - 2.0d0*(qixy*qkxz+qiyy*qkyz+qiyz*qkzz
     &                             -qixz*qkxy-qiyz*qkyy-qizz*qkyz)
               dqiqky = diz*qrkx - dix*qrkz + dkz*qrix - dkx*qriz
     &                     - 2.0d0*(qixz*qkxx+qiyz*qkxy+qizz*qkxz
     &                             -qixx*qkxz-qixy*qkyz-qixz*qkzz)
               dqiqkz = dix*qrky - diy*qrkx + dkx*qriy - dky*qrix
     &                     - 2.0d0*(qixx*qkxy+qixy*qkyy+qixz*qkyz
     &                             -qixy*qkxx-qiyy*qkxy-qiyz*qkxz)
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
                  scaleik(11) = 1.0d0 - termi*(1.0d0 + dampi +
     &                 (4.0d0/9.0d0)*dampi**2 +
     &                 (1.0d0/9.0d0)*dampi**3 +
     &                 (1.0d0/63.0d0)*dampi**4 + 
     &                 (1.0d0/945.0d0)*dampi**5)*expdampi - 
     &                 termk*(1.0d0 + dampk +
     &                 (4.0d0/9.0d0)*dampk**2+
     &                 (1.0d0/9.0d0)*dampk**3+
     &                 (1.0d0/63.0d0)*dampk**4 +
     &                 (1.0d0/945.0d0)*dampk**5)*expdampk
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
                  scaleik(11) = 1.0d0-(1.0d0+dampi+0.5d0*dampi**2
     &                 + (1.0d0/6.0d0)*dampi**3 
     &                 + (5.0d0/126.0d0)*dampi**4
     &                 + (2.0d0/315.0d0)*dampi**5
     &                 + (1.0d0/1890.0d0)*dampi**6)*expdampi
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
c     compute the energy contribution for this interaction
c
               e = term1*rr1 + term2*rr3 + term3*rr5
     &                + term4*rr7 + term5*rr9
               if (ii .eq. kk)  e = 0.5d0 * e
               em = em + e
               einter = einter + e
c
c     calculate intermediate terms for force and torque
c
               de =  nuci*nuck*rr3 
     &             + nuci * qk*rr3*scalek(3)
     &             + nuck * qi*rr3*scalei(3) 
     &             + qi * qk  *rr3*scaleik(3)
     &             + nuck*dri *rr5*scalei(5) 
     &             - nuci*drk *rr5*scalek(5)
     &             + (dik + qk*dri - qi*drk)*rr5*scaleik(5)
     &             + nuci*qrrk*rr7*scalek(7) 
     &             + nuck*qrri*rr7*scalei(7)
     &             +(-dri*drk + 2.0d0*(dkqri-diqrk+qik)  
     &             + qi*qrrk + qk*qrri)*rr7*scaleik(7) 
     &             +(dri*qrrk-drk*qrri-4.0d0*qrrik)*rr9*scaleik(9)
     &             + qrri*qrrk*rr11*scaleik(11)

               term1 = -nuck*rr3*scalei(3) 
     &                 - qk*rr3*scaleik(3) 
     &                 + drk*rr5*scaleik(5) 
     &                 - qrrk*rr7*scaleik(7)
               term2 = nuci*rr3*scalek(3) 
     &                 + qi*rr3*scaleik(3)  
     &                 + dri*rr5*scaleik(5) 
     &                 + qrri*rr7*scaleik(7)
               term3 = 2.0d0 * rr5 * scaleik(5) 
               term4 = 2.0d0 *(-nuck*rr5*scalei(5)
     &                       -qk*rr5*scaleik(5)
     &                       +drk*rr7*scaleik(7)
     &                       -qrrk*rr9*scaleik(9))
               term5 = 2.0d0 *(-nuci*rr5*scalek(5)
     &                       -qi*rr5*scaleik(5)
     &                       -dri*rr7*scaleik(7)
     &                       -qrri*rr9*scaleik(9))
               term6 = 4.0d0 * rr7 *scaleik(7) 

c
c     compute the force components for this interaction
c
               frcx = de*xr + term1*dix + term2*dkx
     &                   + term3*(diqkx-dkqix) + term4*qrix
     &                   + term5*qrkx + term6*(qikrx+qkirx)
               frcy = de*yr + term1*diy + term2*dky
     &                   + term3*(diqky-dkqiy) + term4*qriy
     &                   + term5*qrky + term6*(qikry+qkiry)
               frcz = de*zr + term1*diz + term2*dkz
     &                   + term3*(diqkz-dkqiz) + term4*qriz
     &                   + term5*qrkz + term6*(qikrz+qkirz)
c
c     compute the torque components for this interaction
c
               rr3 = rr3*scaleik(3)
               ttmi(1) = -rr3*dikx + term1*dirx + term3*(dqiqkx+dkqixr)
     &                      - term4*qrixr - term6*(qikrxr+qrrx)
               ttmi(2) = -rr3*diky + term1*diry + term3*(dqiqky+dkqiyr)
     &                      - term4*qriyr - term6*(qikryr+qrry)
               ttmi(3) = -rr3*dikz + term1*dirz + term3*(dqiqkz+dkqizr)
     &                      - term4*qrizr - term6*(qikrzr+qrrz)
               ttmk(1) = rr3*dikx + term2*dkrx - term3*(dqiqkx+diqkxr)
     &                      - term5*qrkxr - term6*(qkirxr-qrrx)
               ttmk(2) = rr3*diky + term2*dkry - term3*(dqiqky+diqkyr)
     &                      - term5*qrkyr - term6*(qkiryr-qrry)
               ttmk(3) = rr3*dikz + term2*dkrz - term3*(dqiqkz+diqkzr)
     &                      - term5*qrkzr - term6*(qkirzr-qrrz)
c
c     force and torque scaled for self-interactions and groups
c
               if (ii .eq. kk) then
                  frcx = 0.5d0 * frcx
                  frcy = 0.5d0 * frcy
                  frcz = 0.5d0 * frcz
                  do j = 1, 3
                     ttmi(j) = 0.5d0 * ttmi(j)
                     ttmk(j) = 0.5d0 * ttmk(j)
                  end do
               end if
               if (use_group) then
                  frcx = fgrp * frcx
                  frcy = fgrp * frcy
                  frcz = fgrp * frcz
                  do j = 1, 3
                     ttmi(j) = fgrp * ttmi(j)
                     ttmk(j) = fgrp * ttmk(j)
                  end do
               end if
c
c     increment force-based gradient and torque on first site
c
               dem(1,ii) = dem(1,ii) + frcx
               dem(2,ii) = dem(2,ii) + frcy
               dem(3,ii) = dem(3,ii) + frcz
               tem(1,i) = tem(1,i) + ttmi(1)
               tem(2,i) = tem(2,i) + ttmi(2)
               tem(3,i) = tem(3,i) + ttmi(3)
c
c     increment force-based gradient and torque on second site
c
               dem(1,kk) = dem(1,kk) - frcx
               dem(2,kk) = dem(2,kk) - frcy
               dem(3,kk) = dem(3,kk) - frcz
               tem(1,k) = tem(1,k) + ttmk(1)
               tem(2,k) = tem(2,k) + ttmk(2)
               tem(3,k) = tem(3,k) + ttmk(3)
c
c     increment the virial due to pairwise Cartesian forces
c
               vxx = -xr * frcx
               vxy = -yr * frcx
               vxz = -zr * frcx
               vyy = -yr * frcy
               vyz = -zr * frcy
               vzz = -zr * frcz
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vxy
               vir(3,1) = vir(3,1) + vxz
               vir(1,2) = vir(1,2) + vxy
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vyz
               vir(1,3) = vir(1,3) + vxz
               vir(2,3) = vir(2,3) + vyz
               vir(3,3) = vir(3,3) + vzz
            end if
            end do
   20       continue
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
c     resolve site torques then increment forces and virial
c
      do i = 1, npole
         call torque (i,tem(1,i),fix,fiy,fiz,dem)
         ii = ipole(i)
         iaz = zaxis(i)
         iax = xaxis(i)
         iay = yaxis(i)
         if (iaz .eq. 0)  iaz = ii
         if (iax .eq. 0)  iax = ii
         if (iay .eq. 0)  iay = ii
         xiz = x(iaz) - x(ii)
         yiz = y(iaz) - y(ii)
         ziz = z(iaz) - z(ii)
         xix = x(iax) - x(ii)
         yix = y(iax) - y(ii)
         zix = z(iax) - z(ii)
         xiy = x(iay) - x(ii)
         yiy = y(iay) - y(ii)
         ziy = z(iay) - z(ii)
         vxx = xix*fix(1) + xiy*fiy(1) + xiz*fiz(1)
         vxy = 0.5d0 * (yix*fix(1) + yiy*fiy(1) + yiz*fiz(1)
     &                    + xix*fix(2) + xiy*fiy(2) + xiz*fiz(2))
         vxz = 0.5d0 * (zix*fix(1) + ziy*fiy(1) + ziz*fiz(1)
     &                    + xix*fix(3) + xiy*fiy(3) + xiz*fiz(3)) 
         vyy = yix*fix(2) + yiy*fiy(2) + yiz*fiz(2)
         vyz = 0.5d0 * (zix*fix(2) + ziy*fiy(2) + ziz*fiz(2)
     &                    + yix*fix(3) + yiy*fiy(3) + yiz*fiz(3))
         vzz = zix*fix(3) + ziy*fiy(3) + ziz*fiz(3)
         vir(1,1) = vir(1,1) + vxx
         vir(2,1) = vir(2,1) + vxy
         vir(3,1) = vir(3,1) + vxz
         vir(1,2) = vir(1,2) + vxy
         vir(2,2) = vir(2,2) + vyy
         vir(3,2) = vir(3,2) + vyz
         vir(1,3) = vir(1,3) + vxz
         vir(2,3) = vir(2,3) + vyz
         vir(3,3) = vir(3,3) + vzz
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      deallocate (tem)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine empole1b  --  neighbor list multipole derivs  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "empole1b" calculates the multipole energy and derivatives
c     with respect to Cartesian coordinates using a neighbor list
c
c
      subroutine empole1b
      use sizes
      use atomid 
      use atoms
      use bound
      use polar 
      use chgpot
      use couple
      use deriv
      use energi
      use group
      use inter
      use molcul
      use mplpot
      use mpole
      use neigh
      use shunt
      use usage
      use virial
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      integer ix,iy,iz
      integer kx,ky,kz
      integer iax,iay,iaz
      real*8 e,de,f,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9,rr11
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dikx,diky,dikz
      real*8 dirx,diry,dirz
      real*8 dkrx,dkry,dkrz
      real*8 qrix,qriy,qriz
      real*8 qrkx,qrky,qrkz
      real*8 qrixr,qriyr,qrizr
      real*8 qrkxr,qrkyr,qrkzr
      real*8 qrrx,qrry,qrrz
      real*8 qikrx,qikry,qikrz
      real*8 qkirx,qkiry,qkirz
      real*8 qikrxr,qikryr,qikrzr
      real*8 qkirxr,qkiryr,qkirzr
      real*8 diqkx,diqky,diqkz
      real*8 dkqix,dkqiy,dkqiz
      real*8 diqkxr,diqkyr,diqkzr
      real*8 dkqixr,dkqiyr,dkqizr
      real*8 dqiqkx,dqiqky,dqiqkz
      real*8 dri,drk,qrri,qrrk
      real*8 diqrk,dkqri
      real*8 dik,qik,qrrik
      real*8 term1,term2,term3
      real*8 term4,term5,term6
      real*8 frcx,frcy,frcz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 ttmi(3),ttmk(3)
      real*8 fix(3),fiy(3),fiz(3)
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
      real*8, allocatable :: tem(:,:)
      logical proceed,usei,usek
      character*6 mode
c
c
c     zero out the atomic multipole energy and derivatives
c
      em = 0.0d0
      do i = 1, n
         do j = 1, 3
            dem(j,i) = 0.0d0
         end do
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
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (tem(3,n))
c
c     initialize connected atom scaling and torque arrays
c
      do i = 1, n
         mscale(i) = 1.0d0
         do j = 1, 3
            tem(j,i) = 0.0d0
         end do
      end do
c
c     maximum order to rr11 (for force)
c
      allocate (scalei(11))
      allocate (scalek(11))
      allocate (scaleik(11))
c
c     maximum order to rr11 (for force)
c
      do i = 1, 11
         scalei(i) = 1.0d0
         scalek(i) = 1.0d0
         scaleik(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and scaling coefficients
c
      f = electric / dielec
      mode = 'MPOLE'
      call switch (mode)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private)
!$OMP& shared(npole,ipole,x,y,z,xaxis,yaxis,zaxis,rpole,use,n12,
!$OMP& i12,n13,i13,n14,i14,n15,i15,m2scale,m3scale,m4scale,m5scale,
!$OMP& atomic,penalpha,type,
!$OMP& nelst,elst,use_group,use_intra,use_bounds,off2,f,molcule)
!$OMP& firstprivate(mscale) shared (em,einter,dem,tem,vir)
!$OMP DO reduction(+:em,einter,dem,tem,vir) schedule(guided)
c
c     compute the multipole interaction energy and gradient
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
            if (.not. proceed)  goto 10
            xr = x(kk) - xi
            yr = y(kk) - yi
            zr = z(kk) - zi
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
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
               rr11 = 9.0d0 * rr9 / r2
c
c     intermediates involving moments and distance separation
c
               dikx = diy*dkz - diz*dky
               diky = diz*dkx - dix*dkz
               dikz = dix*dky - diy*dkx
               dirx = diy*zr - diz*yr
               diry = diz*xr - dix*zr
               dirz = dix*yr - diy*xr
               dkrx = dky*zr - dkz*yr
               dkry = dkz*xr - dkx*zr
               dkrz = dkx*yr - dky*xr
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
               qrixr = qriz*yr - qriy*zr
               qriyr = qrix*zr - qriz*xr
               qrizr = qriy*xr - qrix*yr
               qrkxr = qrkz*yr - qrky*zr
               qrkyr = qrkx*zr - qrkz*xr
               qrkzr = qrky*xr - qrkx*yr
               qrrx = qrky*qriz - qrkz*qriy
               qrry = qrkz*qrix - qrkx*qriz
               qrrz = qrkx*qriy - qrky*qrix
               qikrx = qixx*qrkx + qixy*qrky + qixz*qrkz
               qikry = qixy*qrkx + qiyy*qrky + qiyz*qrkz
               qikrz = qixz*qrkx + qiyz*qrky + qizz*qrkz
               qkirx = qkxx*qrix + qkxy*qriy + qkxz*qriz
               qkiry = qkxy*qrix + qkyy*qriy + qkyz*qriz
               qkirz = qkxz*qrix + qkyz*qriy + qkzz*qriz
               qikrxr = qikrz*yr - qikry*zr
               qikryr = qikrx*zr - qikrz*xr
               qikrzr = qikry*xr - qikrx*yr
               qkirxr = qkirz*yr - qkiry*zr
               qkiryr = qkirx*zr - qkirz*xr
               qkirzr = qkiry*xr - qkirx*yr
               diqkx = dix*qkxx + diy*qkxy + diz*qkxz
               diqky = dix*qkxy + diy*qkyy + diz*qkyz
               diqkz = dix*qkxz + diy*qkyz + diz*qkzz
               dkqix = dkx*qixx + dky*qixy + dkz*qixz
               dkqiy = dkx*qixy + dky*qiyy + dkz*qiyz
               dkqiz = dkx*qixz + dky*qiyz + dkz*qizz
               diqrk = dix*qrkx + diy*qrky + diz*qrkz
               dkqri = dkx*qrix + dky*qriy + dkz*qriz
               diqkxr = diqkz*yr - diqky*zr
               diqkyr = diqkx*zr - diqkz*xr
               diqkzr = diqky*xr - diqkx*yr
               dkqixr = dkqiz*yr - dkqiy*zr
               dkqiyr = dkqix*zr - dkqiz*xr
               dkqizr = dkqiy*xr - dkqix*yr
               dqiqkx = diy*qrkz - diz*qrky + dky*qriz - dkz*qriy
     &                     - 2.0d0*(qixy*qkxz+qiyy*qkyz+qiyz*qkzz
     &                             -qixz*qkxy-qiyz*qkyy-qizz*qkyz)
               dqiqky = diz*qrkx - dix*qrkz + dkz*qrix - dkx*qriz
     &                     - 2.0d0*(qixz*qkxx+qiyz*qkxy+qizz*qkxz
     &                             -qixx*qkxz-qixy*qkyz-qixz*qkzz)
               dqiqkz = dix*qrky - diy*qrkx + dkx*qriy - dky*qrix
     &                     - 2.0d0*(qixx*qkxy+qixy*qkyy+qixz*qkyz
     &                             -qixy*qkxx-qiyy*qkxy-qiyz*qkxz)
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
                  scaleik(11) = 1.0d0 - termi*(1.0d0 + dampi +
     &                 (4.0d0/9.0d0)*dampi**2 +
     &                 (1.0d0/9.0d0)*dampi**3 +
     &                 (1.0d0/63.0d0)*dampi**4 + 
     &                 (1.0d0/945.0d0)*dampi**5)*expdampi - 
     &                 termk*(1.0d0 + dampk +
     &                 (4.0d0/9.0d0)*dampk**2+
     &                 (1.0d0/9.0d0)*dampk**3+
     &                 (1.0d0/63.0d0)*dampk**4 +
     &                 (1.0d0/945.0d0)*dampk**5)*expdampk
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
                  scaleik(11) = 1.0d0-(1.0d0+dampi+0.5d0*dampi**2
     &                 + (1.0d0/6.0d0)*dampi**3 
     &                 + (5.0d0/126.0d0)*dampi**4
     &                 + (2.0d0/315.0d0)*dampi**5
     &                 + (1.0d0/1890.0d0)*dampi**6)*expdampi
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
c     compute the energy contribution for this interaction
c
               e = term1*rr1 + term2*rr3 + term3*rr5
     &                + term4*rr7 + term5*rr9
               if (use_group)  e = e * fgrp
               em = em + e
               if (molcule(ii) .ne. molcule(kk))
     &            einter = einter + e
c
c     calculate intermediate terms for force and torque
c
               de =  nuci*nuck*rr3 
     &             + nuci * qk*rr3*scalek(3)
     &             + nuck * qi*rr3*scalei(3) 
     &             + qi * qk  *rr3*scaleik(3)
     &             + nuck*dri *rr5*scalei(5) 
     &             - nuci*drk *rr5*scalek(5)
     &             + (dik + qk*dri - qi*drk)*rr5*scaleik(5)
     &             + nuci*qrrk*rr7*scalek(7) 
     &             + nuck*qrri*rr7*scalei(7)
     &             +(-dri*drk + 2.0d0*(dkqri-diqrk+qik)  
     &             + qi*qrrk + qk*qrri)*rr7*scaleik(7) 
     &             +(dri*qrrk-drk*qrri-4.0d0*qrrik)*rr9*scaleik(9)
     &             + qrri*qrrk*rr11*scaleik(11)

               term1 = -nuck*rr3*scalei(3) 
     &                 - qk*rr3*scaleik(3) 
     &                 + drk*rr5*scaleik(5) 
     &                 - qrrk*rr7*scaleik(7)
               term2 = nuci*rr3*scalek(3) 
     &                 + qi*rr3*scaleik(3)  
     &                 + dri*rr5*scaleik(5) 
     &                 + qrri*rr7*scaleik(7)
               term3 = 2.0d0 * rr5 * scaleik(5) 
               term4 = 2.0d0 *(-nuck*rr5*scalei(5)
     &                       -qk*rr5*scaleik(5)
     &                       +drk*rr7*scaleik(7)
     &                       -qrrk*rr9*scaleik(9))
               term5 = 2.0d0 *(-nuci*rr5*scalek(5)
     &                       -qi*rr5*scaleik(5)
     &                       -dri*rr7*scaleik(7)
     &                       -qrri*rr9*scaleik(9))
               term6 = 4.0d0 * rr7 *scaleik(7) 

c
c     compute the force components for this interaction
c
               frcx = de*xr + term1*dix + term2*dkx
     &                   + term3*(diqkx-dkqix) + term4*qrix
     &                   + term5*qrkx + term6*(qikrx+qkirx)
               frcy = de*yr + term1*diy + term2*dky
     &                   + term3*(diqky-dkqiy) + term4*qriy
     &                   + term5*qrky + term6*(qikry+qkiry)
               frcz = de*zr + term1*diz + term2*dkz
     &                   + term3*(diqkz-dkqiz) + term4*qriz
     &                   + term5*qrkz + term6*(qikrz+qkirz)
c
c     compute the torque components for this interaction
c
               rr3 = rr3*scaleik(3)

               ttmi(1) = -rr3*dikx + term1*dirx + term3*(dqiqkx+dkqixr)
     &                      - term4*qrixr - term6*(qikrxr+qrrx)
               ttmi(2) = -rr3*diky + term1*diry + term3*(dqiqky+dkqiyr)
     &                      - term4*qriyr - term6*(qikryr+qrry)
               ttmi(3) = -rr3*dikz + term1*dirz + term3*(dqiqkz+dkqizr)
     &                      - term4*qrizr - term6*(qikrzr+qrrz)
               ttmk(1) = rr3*dikx + term2*dkrx - term3*(dqiqkx+diqkxr)
     &                      - term5*qrkxr - term6*(qkirxr-qrrx)
               ttmk(2) = rr3*diky + term2*dkry - term3*(dqiqky+diqkyr)
     &                      - term5*qrkyr - term6*(qkiryr-qrry)
               ttmk(3) = rr3*dikz + term2*dkrz - term3*(dqiqkz+diqkzr)
     &                      - term5*qrkzr - term6*(qkirzr-qrrz)
c
c     increment force-based gradient and torque on first site
c
               dem(1,ii) = dem(1,ii) + frcx
               dem(2,ii) = dem(2,ii) + frcy
               dem(3,ii) = dem(3,ii) + frcz
               tem(1,i) = tem(1,i) + ttmi(1)
               tem(2,i) = tem(2,i) + ttmi(2)
               tem(3,i) = tem(3,i) + ttmi(3)
c
c     increment force-based gradient and torque on second site
c
               dem(1,kk) = dem(1,kk) - frcx
               dem(2,kk) = dem(2,kk) - frcy
               dem(3,kk) = dem(3,kk) - frcz
               tem(1,k) = tem(1,k) + ttmk(1)
               tem(2,k) = tem(2,k) + ttmk(2)
               tem(3,k) = tem(3,k) + ttmk(3)
c
c     increment the virial due to pairwise Cartesian forces
c
               vxx = -xr * frcx
               vxy = -yr * frcx
               vxz = -zr * frcx
               vyy = -yr * frcy
               vyz = -zr * frcy
               vzz = -zr * frcz
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vxy
               vir(3,1) = vir(3,1) + vxz
               vir(1,2) = vir(1,2) + vxy
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vyz
               vir(1,3) = vir(1,3) + vxz
               vir(2,3) = vir(2,3) + vyz
               vir(3,3) = vir(3,3) + vzz
            end if
   10       continue
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
!$OMP DO reduction(+:dem,vir) schedule(guided)
c
c     resolve site torques then increment forces and virial
c
      do i = 1, npole
         call torque (i,tem(1,i),fix,fiy,fiz,dem)
         ii = ipole(i)
         iaz = zaxis(i)
         iax = xaxis(i)
         iay = yaxis(i)
         if (iaz .eq. 0)  iaz = ii
         if (iax .eq. 0)  iax = ii
         if (iay .eq. 0)  iay = ii
         xiz = x(iaz) - x(ii)
         yiz = y(iaz) - y(ii)
         ziz = z(iaz) - z(ii)
         xix = x(iax) - x(ii)
         yix = y(iax) - y(ii)
         zix = z(iax) - z(ii)
         xiy = x(iay) - x(ii)
         yiy = y(iay) - y(ii)
         ziy = z(iay) - z(ii)
         vxx = xix*fix(1) + xiy*fiy(1) + xiz*fiz(1)
         vxy = 0.5d0 * (yix*fix(1) + yiy*fiy(1) + yiz*fiz(1)
     &                    + xix*fix(2) + xiy*fiy(2) + xiz*fiz(2))
         vxz = 0.5d0 * (zix*fix(1) + ziy*fiy(1) + ziz*fiz(1)
     &                    + xix*fix(3) + xiy*fiy(3) + xiz*fiz(3)) 
         vyy = yix*fix(2) + yiy*fiy(2) + yiz*fiz(2)
         vyz = 0.5d0 * (zix*fix(2) + ziy*fiy(2) + ziz*fiz(2)
     &                    + yix*fix(3) + yiy*fiy(3) + yiz*fiz(3))
         vzz = zix*fix(3) + ziy*fiy(3) + ziz*fiz(3)
         vir(1,1) = vir(1,1) + vxx
         vir(2,1) = vir(2,1) + vxy
         vir(3,1) = vir(3,1) + vxz
         vir(1,2) = vir(1,2) + vxy
         vir(2,2) = vir(2,2) + vyy
         vir(3,2) = vir(3,2) + vyz
         vir(1,3) = vir(1,3) + vxz
         vir(2,3) = vir(2,3) + vyz
         vir(3,3) = vir(3,3) + vzz
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
      deallocate (tem)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine empole1c  --  Ewald multipole derivs via loop  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "empole1c" calculates the multipole energy and derivatives
c     with respect to Cartesian coordinates using particle mesh
c     Ewald summation and a double loop
c
c
      subroutine empole1c
      use sizes
      use atoms
      use boxes
      use chgpot
      use deriv
      use energi
      use ewald
      use math
      use mpole
      use virial
      implicit none
      integer i,j,ii
      real*8 e,f
      real*8 term,fterm
      real*8 cii,dii,qii
      real*8 xd,yd,zd
      real*8 xq,yq,zq
      real*8 xv,yv,zv,vterm
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 xdfield,ydfield
      real*8 zdfield
      real*8 trq(3),frcx(3)
      real*8 frcy(3),frcz(3)
c
c
c     zero out the atomic multipole energy and derivatives
c
      em = 0.0d0
      do i = 1, n
         do j = 1, 3
            dem(j,i) = 0.0d0
         end do
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
      call emreal1c
c
c     compute the reciprocal space part of the Ewald summation
c
      call emrecip1
c
c     compute the Ewald self-energy term over all the atoms
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
            xd = xd + rpole(2,i) + rpole(1,i)*x(ii)
            yd = yd + rpole(3,i) + rpole(1,i)*y(ii)
            zd = zd + rpole(4,i) + rpole(1,i)*z(ii)
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         em = em + term*(xd*xd+yd*yd+zd*zd)
         do i = 1, npole
            ii = ipole(i)
            dem(1,ii) = dem(1,ii) + 2.0d0*term*rpole(1,i)*xd
            dem(2,ii) = dem(2,ii) + 2.0d0*term*rpole(1,i)*yd
            dem(3,ii) = dem(3,ii) + 2.0d0*term*rpole(1,i)*zd
         end do
         xdfield = -2.0d0 * term * xd
         ydfield = -2.0d0 * term * yd
         zdfield = -2.0d0 * term * zd
         do i = 1, npole
            trq(1) = rpole(3,i)*zdfield - rpole(4,i)*ydfield
            trq(2) = rpole(4,i)*xdfield - rpole(2,i)*zdfield
            trq(3) = rpole(2,i)*ydfield - rpole(3,i)*xdfield
            call torque (i,trq,frcx,frcy,frcz,dem)
         end do
c
c     boundary correction to virial due to overall cell dipole
c
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         xq = 0.0d0
         yq = 0.0d0
         zq = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            xd = xd + rpole(2,i)
            yd = yd + rpole(3,i)
            zd = zd + rpole(4,i)
            xq = xq + rpole(1,i)*x(ii)
            yq = yq + rpole(1,i)*y(ii)
            zq = zq + rpole(1,i)*z(ii)
         end do
         xv = xd * xq
         yv = yd * yq
         zv = zd * zq
         vterm = term * (xd*xd + yd*yd + zd*zd + 2.0d0*(xv+yv+zv)
     &                      + xq*xq + yq*yq + zq*zq)
         vir(1,1) = vir(1,1) + 2.0d0*term*(xq*xq+xv) + vterm
         vir(2,1) = vir(2,1) + 2.0d0*term*(xq*yq+xv)
         vir(3,1) = vir(3,1) + 2.0d0*term*(xq*zq+xv)
         vir(1,2) = vir(1,2) + 2.0d0*term*(yq*xq+yv)
         vir(2,2) = vir(2,2) + 2.0d0*term*(yq*yq+yv) + vterm
         vir(3,2) = vir(3,2) + 2.0d0*term*(yq*zq+yv)
         vir(1,3) = vir(1,3) + 2.0d0*term*(zq*xq+zv)
         vir(2,3) = vir(2,3) + 2.0d0*term*(zq*yq+zv)
         vir(3,3) = vir(3,3) + 2.0d0*term*(zq*zq+zv) + vterm
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine emreal1c  --  Ewald real space derivs via loop  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "emreal1c" evaluates the real space portion of the Ewald
c     summation energy and gradient due to multipole interactions
c     via a double loop
c
c
      subroutine emreal1c
      use sizes
      use atomid 
      use atoms
      use bound
      use cell
      use polar 
      use chgpot
      use couple
      use deriv
      use energi
      use ewald
      use inter
      use math
      use molcul
      use mplpot
      use mpole
      use shunt
      use virial
      implicit none
      integer i,j,k
      integer ii,kk,jcell
      integer iax,iay,iaz
      real*8 e,efull,de,f
      real*8 bfac,erfc
      real*8 alsq2,alsq2n
      real*8 exp2a,ralpha
      real*8 scalekk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9,rr11
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dikx,diky,dikz
      real*8 dirx,diry,dirz
      real*8 dkrx,dkry,dkrz
      real*8 qrix,qriy,qriz
      real*8 qrkx,qrky,qrkz
      real*8 qrixr,qriyr,qrizr
      real*8 qrkxr,qrkyr,qrkzr
      real*8 qrrx,qrry,qrrz
      real*8 qikrx,qikry,qikrz
      real*8 qkirx,qkiry,qkirz
      real*8 qikrxr,qikryr,qikrzr
      real*8 qkirxr,qkiryr,qkirzr
      real*8 diqkx,diqky,diqkz
      real*8 dkqix,dkqiy,dkqiz
      real*8 diqkxr,diqkyr,diqkzr
      real*8 dkqixr,dkqiyr,dkqizr
      real*8 dqiqkx,dqiqky,dqiqkz
      real*8 dri,drk,qrri,qrrk
      real*8 diqrk,dkqri
      real*8 dik,qik,qrrik
      real*8 term1,term2,term3
      real*8 term4,term5,term6
      real*8 frcx,frcy,frcz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 ttmi(3),ttmk(3)
      real*8 fix(3),fiy(3),fiz(3)
      real*8 bn(0:5)
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
      real*8, allocatable :: tem(:,:)
      character*6 mode
      external erfc
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (tem(3,n))
c
c     maximum order to rr11 (for force)
c
      allocate (scalei(11))
      allocate (scalek(11))
      allocate (scaleik(11))
c
c     maximum order to rr11 (for force)
c
      do i = 1, 11 
         scalei(i) = 1.0d0
         scalek(i) = 1.0d0
         scaleik(i) = 1.0d0
      end do
c
c     initialize connected atom scaling and torque arrays
c
      do i = 1, n
         mscale(i) = 1.0d0
         do j = 1, 3
            tem(j,i) = 0.0d0
         end do
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
            r2 = xr*xr + yr*yr + zr*zr
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
               rr11 = 9.0d0 * rr9 / r2
c
c     calculate the real space Ewald error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 5
                  bfac = dble(j+j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
               do j = 0, 5
                  bn(j) = f * bn(j)
               end do
c
c     intermediates involving moments and distance separation
c
               dikx = diy*dkz - diz*dky
               diky = diz*dkx - dix*dkz
               dikz = dix*dky - diy*dkx
               dirx = diy*zr - diz*yr
               diry = diz*xr - dix*zr
               dirz = dix*yr - diy*xr
               dkrx = dky*zr - dkz*yr
               dkry = dkz*xr - dkx*zr
               dkrz = dkx*yr - dky*xr
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
               qrixr = qriz*yr - qriy*zr
               qriyr = qrix*zr - qriz*xr
               qrizr = qriy*xr - qrix*yr
               qrkxr = qrkz*yr - qrky*zr
               qrkyr = qrkx*zr - qrkz*xr
               qrkzr = qrky*xr - qrkx*yr
               qrrx = qrky*qriz - qrkz*qriy
               qrry = qrkz*qrix - qrkx*qriz
               qrrz = qrkx*qriy - qrky*qrix
               qikrx = qixx*qrkx + qixy*qrky + qixz*qrkz
               qikry = qixy*qrkx + qiyy*qrky + qiyz*qrkz
               qikrz = qixz*qrkx + qiyz*qrky + qizz*qrkz
               qkirx = qkxx*qrix + qkxy*qriy + qkxz*qriz
               qkiry = qkxy*qrix + qkyy*qriy + qkyz*qriz
               qkirz = qkxz*qrix + qkyz*qriy + qkzz*qriz
               qikrxr = qikrz*yr - qikry*zr
               qikryr = qikrx*zr - qikrz*xr
               qikrzr = qikry*xr - qikrx*yr
               qkirxr = qkirz*yr - qkiry*zr
               qkiryr = qkirx*zr - qkirz*xr
               qkirzr = qkiry*xr - qkirx*yr
               diqkx = dix*qkxx + diy*qkxy + diz*qkxz
               diqky = dix*qkxy + diy*qkyy + diz*qkyz
               diqkz = dix*qkxz + diy*qkyz + diz*qkzz
               dkqix = dkx*qixx + dky*qixy + dkz*qixz
               dkqiy = dkx*qixy + dky*qiyy + dkz*qiyz
               dkqiz = dkx*qixz + dky*qiyz + dkz*qizz
               diqrk = dix*qrkx + diy*qrky + diz*qrkz
               dkqri = dkx*qrix + dky*qriy + dkz*qriz
               diqkxr = diqkz*yr - diqky*zr
               diqkyr = diqkx*zr - diqkz*xr
               diqkzr = diqky*xr - diqkx*yr
               dkqixr = dkqiz*yr - dkqiy*zr
               dkqiyr = dkqix*zr - dkqiz*xr
               dkqizr = dkqiy*xr - dkqix*yr
               dqiqkx = diy*qrkz - diz*qrky + dky*qriz - dkz*qriy
     &                     - 2.0d0*(qixy*qkxz+qiyy*qkyz+qiyz*qkzz
     &                             -qixz*qkxy-qiyz*qkyy-qizz*qkyz)
               dqiqky = diz*qrkx - dix*qrkz + dkz*qrix - dkx*qriz
     &                     - 2.0d0*(qixz*qkxx+qiyz*qkxy+qizz*qkxz
     &                             -qixx*qkxz-qixy*qkyz-qixz*qkzz)
               dqiqkz = dix*qrky - diy*qrkx + dkx*qriy - dky*qrix
     &                     - 2.0d0*(qixx*qkxy+qixy*qkyy+qixz*qkyz
     &                             -qixy*qkxx-qiyy*qkxy-qiyz*qkxz)
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
                  scaleik(11) = 1.0d0 - termi*(1.0d0 + dampi +
     &                 (4.0d0/9.0d0)*dampi**2 +
     &                 (1.0d0/9.0d0)*dampi**3 +
     &                 (1.0d0/63.0d0)*dampi**4 + 
     &                 (1.0d0/945.0d0)*dampi**5)*expdampi - 
     &                 termk*(1.0d0 + dampk +
     &                 (4.0d0/9.0d0)*dampk**2+
     &                 (1.0d0/9.0d0)*dampk**3+
     &                 (1.0d0/63.0d0)*dampk**4 +
     &                 (1.0d0/945.0d0)*dampk**5)*expdampk
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
                  scaleik(11) = 1.0d0-(1.0d0+dampi+0.5d0*dampi**2
     &                 + (1.0d0/6.0d0)*dampi**3 
     &                 + (5.0d0/126.0d0)*dampi**4
     &                 + (2.0d0/315.0d0)*dampi**5
     &                 + (1.0d0/1890.0d0)*dampi**6)*expdampi
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
c     compute the full energy without any Ewald scaling
c
               efull = term1*rr1 + term2*rr3 + term3*rr5
     &                + term4*rr7 + term5*rr9
               efull = efull * mscale(kk)
               if (molcule(ii) .ne. molcule(kk))
     &            einter = einter + efull
c
c     modify distances to account for Ewald and exclusions
c
               bn(0) = bn(0)/rr1
               bn(1) = bn(1)/rr3
               bn(2) = bn(2)/rr5
               bn(3) = bn(3)/rr7
               bn(4) = bn(4)/rr9
               bn(5) = bn(5)/rr11

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
c     calculate intermediate terms for force and torque
c
               de =  nuci*nuck*rr3*(bn(1)-(1.0d0-mscale(kk))) 
     &           + nuci*qk*rr3*(bn(1)-(1.0d0-scalek(3)*mscale(kk)))
     &           + nuck * qi*rr3*(bn(1)-(1.0d0-scalei(3)*mscale(kk)))
     &           + qi*qk*rr3*(bn(1)-(1.0d0-scaleik(3)*mscale(kk)))
     &           + nuck*dri*rr5*(bn(2)-(1.0d0-scalei(5)*mscale(kk))) 
     &           - nuci*drk*rr5*(bn(2)-(1.0d0-scalek(5)*mscale(kk)))
     &           + (dik + qk*dri - qi*drk)*rr5
     &            *(bn(2)-(1.0d0-scaleik(5)*mscale(kk)))
     &           + nuci*qrrk*rr7*(bn(3)-(1.0d0-scalek(7)*mscale(kk)))
     &           + nuck*qrri*rr7*(bn(3)-(1.0d0-scalei(7)*mscale(kk)))
     &           +(-dri*drk + 2.0d0*(dkqri-diqrk+qik)  
     &           + qi*qrrk + qk*qrri)*rr7
     &            *(bn(3)-(1.0d0-scaleik(7)*mscale(kk)))
     &           +(dri*qrrk-drk*qrri-4.0d0*qrrik)*rr9
     *            *(bn(4)-(1.0d0-scaleik(9)*mscale(kk)))
     &           + qrri*qrrk*rr11
     &            *(bn(5)-(1.0d0-scaleik(11)*mscale(kk)))

               term1 =-nuck*rr3*(bn(1)-(1.0d0-mscale(kk)*scalei(3))) 
     &            - qk*rr3*(bn(1)-(1.0d0-mscale(kk)*scaleik(3))) 
     &            + drk*rr5*(bn(2)-(1.0d0-mscale(kk)*scaleik(5))) 
     &            - qrrk*rr7*(bn(3)-(1.0d0-mscale(kk)*scaleik(7)))
               term2 = nuci*rr3*(bn(1)-(1.0d0-mscale(kk)*scalek(3))) 
     &            + qi*rr3*(bn(1)-(1.0d0-mscale(kk)*scaleik(3)))  
     &            + dri*rr5*(bn(2)-(1.0d0-mscale(kk)*scaleik(5))) 
     &            + qrri*rr7*(bn(3)-(1.0d0-mscale(kk)*scaleik(7)))
               term3=2.0d0*rr5*(bn(2)-(1.0d0-mscale(kk)*scaleik(5))) 
               term4=2.0d0*(-nuck*rr5
     &             *(bn(2)-(1.0d0-mscale(kk)*scalei(5)))
     &            -qk*rr5*(bn(2)-(1.0d0-mscale(kk)*scaleik(5)))
     &            +drk*rr7*(bn(3)-(1.0d0-mscale(kk)*scaleik(7)))
     &            -qrrk*rr9*(bn(4)-(1.0d0-mscale(kk)*scaleik(9))))
               term5 = 2.0d0 *(-nuci*rr5
     &             *(bn(2)-(1.0d0-mscale(kk)*scalek(5)))
     &            -qi*rr5*(bn(2)-(1.0d0-mscale(kk)*scaleik(5)))
     &            -dri*rr7*(bn(3)-(1.0d0-mscale(kk)*scaleik(7)))
     &            -qrri*rr9*(bn(4)-(1.0d0-mscale(kk)*scaleik(9))))
               term6 = 4.0d0 * rr7 
     &             *(bn(3)-(1.0d0-mscale(kk)*scaleik(7))) 

c
c     compute the force components for this interaction
c
               frcx = de*xr + term1*dix + term2*dkx
     &                   + term3*(diqkx-dkqix) + term4*qrix
     &                   + term5*qrkx + term6*(qikrx+qkirx)
               frcy = de*yr + term1*diy + term2*dky
     &                   + term3*(diqky-dkqiy) + term4*qriy
     &                   + term5*qrky + term6*(qikry+qkiry)
               frcz = de*zr + term1*diz + term2*dkz
     &                   + term3*(diqkz-dkqiz) + term4*qriz
     &                   + term5*qrkz + term6*(qikrz+qkirz)
c
c     compute the torque components for this interaction
c
               rr3 =rr3*(bn(1)-(1.0d0-scaleik(3)*mscale(kk)))

               ttmi(1) = -rr3*dikx + term1*dirx + term3*(dqiqkx+dkqixr)
     &                      - term4*qrixr - term6*(qikrxr+qrrx)
               ttmi(2) = -rr3*diky + term1*diry + term3*(dqiqky+dkqiyr)
     &                      - term4*qriyr - term6*(qikryr+qrry)
               ttmi(3) = -rr3*dikz + term1*dirz + term3*(dqiqkz+dkqizr)
     &                      - term4*qrizr - term6*(qikrzr+qrrz)
               ttmk(1) = rr3*dikx + term2*dkrx - term3*(dqiqkx+diqkxr)
     &                      - term5*qrkxr - term6*(qkirxr-qrrx)
               ttmk(2) = rr3*diky + term2*dkry - term3*(dqiqky+diqkyr)
     &                      - term5*qrkyr - term6*(qkiryr-qrry)
               ttmk(3) = rr3*dikz + term2*dkrz - term3*(dqiqkz+diqkzr)
     &                      - term5*qrkzr - term6*(qkirzr-qrrz)
c
c     increment force-based gradient and torque on first site
c
               dem(1,ii) = dem(1,ii) + frcx
               dem(2,ii) = dem(2,ii) + frcy
               dem(3,ii) = dem(3,ii) + frcz
               tem(1,i) = tem(1,i) + ttmi(1)
               tem(2,i) = tem(2,i) + ttmi(2)
               tem(3,i) = tem(3,i) + ttmi(3)
c
c     increment force-based gradient and torque on second site
c
               dem(1,kk) = dem(1,kk) - frcx
               dem(2,kk) = dem(2,kk) - frcy
               dem(3,kk) = dem(3,kk) - frcz
               tem(1,k) = tem(1,k) + ttmk(1)
               tem(2,k) = tem(2,k) + ttmk(2)
               tem(3,k) = tem(3,k) + ttmk(3)
c
c     increment the virial due to pairwise Cartesian forces
c
               vxx = -xr * frcx
               vxy = -yr * frcx
               vxz = -zr * frcx
               vyy = -yr * frcy
               vyz = -zr * frcy
               vzz = -zr * frcz
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vxy
               vir(3,1) = vir(3,1) + vxz
               vir(1,2) = vir(1,2) + vxy
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vyz
               vir(1,3) = vir(1,3) + vxz
               vir(2,3) = vir(2,3) + vyz
               vir(3,3) = vir(3,3) + vzz
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
c     calculate interaction with other unit cells
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
            r2 = xr*xr + yr*yr + zr*zr
            if (.not. (use_polymer .and. r2.le.polycut2)) then
               mscale(kk) = 1.0d0
            end if
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
               rr11 = 9.0d0 * rr9 / r2
c
c     calculate the real space Ewald error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 5
                  bfac = dble(j+j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
               do j = 0, 5
                  bn(j) = f * bn(j)
               end do
c
c     intermediates involving moments and distance separation
c
               dikx = diy*dkz - diz*dky
               diky = diz*dkx - dix*dkz
               dikz = dix*dky - diy*dkx
               dirx = diy*zr - diz*yr
               diry = diz*xr - dix*zr
               dirz = dix*yr - diy*xr
               dkrx = dky*zr - dkz*yr
               dkry = dkz*xr - dkx*zr
               dkrz = dkx*yr - dky*xr
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
               qrixr = qriz*yr - qriy*zr
               qriyr = qrix*zr - qriz*xr
               qrizr = qriy*xr - qrix*yr
               qrkxr = qrkz*yr - qrky*zr
               qrkyr = qrkx*zr - qrkz*xr
               qrkzr = qrky*xr - qrkx*yr
               qrrx = qrky*qriz - qrkz*qriy
               qrry = qrkz*qrix - qrkx*qriz
               qrrz = qrkx*qriy - qrky*qrix
               qikrx = qixx*qrkx + qixy*qrky + qixz*qrkz
               qikry = qixy*qrkx + qiyy*qrky + qiyz*qrkz
               qikrz = qixz*qrkx + qiyz*qrky + qizz*qrkz
               qkirx = qkxx*qrix + qkxy*qriy + qkxz*qriz
               qkiry = qkxy*qrix + qkyy*qriy + qkyz*qriz
               qkirz = qkxz*qrix + qkyz*qriy + qkzz*qriz
               qikrxr = qikrz*yr - qikry*zr
               qikryr = qikrx*zr - qikrz*xr
               qikrzr = qikry*xr - qikrx*yr
               qkirxr = qkirz*yr - qkiry*zr
               qkiryr = qkirx*zr - qkirz*xr
               qkirzr = qkiry*xr - qkirx*yr
               diqkx = dix*qkxx + diy*qkxy + diz*qkxz
               diqky = dix*qkxy + diy*qkyy + diz*qkyz
               diqkz = dix*qkxz + diy*qkyz + diz*qkzz
               dkqix = dkx*qixx + dky*qixy + dkz*qixz
               dkqiy = dkx*qixy + dky*qiyy + dkz*qiyz
               dkqiz = dkx*qixz + dky*qiyz + dkz*qizz
               diqrk = dix*qrkx + diy*qrky + diz*qrkz
               dkqri = dkx*qrix + dky*qriy + dkz*qriz
               diqkxr = diqkz*yr - diqky*zr
               diqkyr = diqkx*zr - diqkz*xr
               diqkzr = diqky*xr - diqkx*yr
               dkqixr = dkqiz*yr - dkqiy*zr
               dkqiyr = dkqix*zr - dkqiz*xr
               dkqizr = dkqiy*xr - dkqix*yr
               dqiqkx = diy*qrkz - diz*qrky + dky*qriz - dkz*qriy
     &                     - 2.0d0*(qixy*qkxz+qiyy*qkyz+qiyz*qkzz
     &                             -qixz*qkxy-qiyz*qkyy-qizz*qkyz)
               dqiqky = diz*qrkx - dix*qrkz + dkz*qrix - dkx*qriz
     &                     - 2.0d0*(qixz*qkxx+qiyz*qkxy+qizz*qkxz
     &                             -qixx*qkxz-qixy*qkyz-qixz*qkzz)
               dqiqkz = dix*qrky - diy*qrkx + dkx*qriy - dky*qrix
     &                     - 2.0d0*(qixx*qkxy+qixy*qkyy+qixz*qkyz
     &                             -qixy*qkxx-qiyy*qkxy-qiyz*qkxz)
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
                  scaleik(11) = 1.0d0 - termi*(1.0d0 + dampi +
     &                 (4.0d0/9.0d0)*dampi**2 +
     &                 (1.0d0/9.0d0)*dampi**3 +
     &                 (1.0d0/63.0d0)*dampi**4 + 
     &                 (1.0d0/945.0d0)*dampi**5)*expdampi - 
     &                 termk*(1.0d0 + dampk +
     &                 (4.0d0/9.0d0)*dampk**2+
     &                 (1.0d0/9.0d0)*dampk**3+
     &                 (1.0d0/63.0d0)*dampk**4 +
     &                 (1.0d0/945.0d0)*dampk**5)*expdampk
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
                  scaleik(11) = 1.0d0-(1.0d0+dampi+0.5d0*dampi**2
     &                 + (1.0d0/6.0d0)*dampi**3 
     &                 + (5.0d0/126.0d0)*dampi**4
     &                 + (2.0d0/315.0d0)*dampi**5
     &                 + (1.0d0/1890.0d0)*dampi**6)*expdampi
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
c     compute the full energy without any Ewald scaling
c
               efull = term1*rr1 + term2*rr3 + term3*rr5
     &                    + term4*rr7 + term5*rr9
               efull = efull * mscale(kk)
               if (ii .eq. kk)  efull = 0.5d0 * e
               einter = einter + efull
c
c     modify distances to account for Ewald and exclusions
c
               bn(0) = bn(0)/rr1
               bn(1) = bn(1)/rr3
               bn(2) = bn(2)/rr5
               bn(3) = bn(3)/rr7
               bn(4) = bn(4)/rr9
               bn(5) = bn(5)/rr11

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
               if (ii .eq. kk)  e = 0.5d0 * e
               em = em + e
c
c     calculate intermediate terms for force and torque
c
               de =  nuci*nuck*rr3*(bn(1)-(1.0d0-mscale(kk))) 
     &           + nuci*qk*rr3*(bn(1)-(1.0d0-scalek(3)*mscale(kk)))
     &           + nuck * qi*rr3*(bn(1)-(1.0d0-scalei(3)*mscale(kk)))
     &           + qi*qk*rr3*(bn(1)-(1.0d0-scaleik(3)*mscale(kk)))
     &           + nuck*dri*rr5*(bn(2)-(1.0d0-scalei(5)*mscale(kk))) 
     &           - nuci*drk*rr5*(bn(2)-(1.0d0-scalek(5)*mscale(kk)))
     &           + (dik + qk*dri - qi*drk)*rr5
     &            *(bn(2)-(1.0d0-scaleik(5)*mscale(kk)))
     &           + nuci*qrrk*rr7*(bn(3)-(1.0d0-scalek(7)*mscale(kk)))
     &           + nuck*qrri*rr7*(bn(3)-(1.0d0-scalei(7)*mscale(kk)))
     &           +(-dri*drk + 2.0d0*(dkqri-diqrk+qik)  
     &           + qi*qrrk + qk*qrri)*rr7
     &            *(bn(3)-(1.0d0-scaleik(7)*mscale(kk)))
     &           +(dri*qrrk-drk*qrri-4.0d0*qrrik)*rr9
     *            *(bn(4)-(1.0d0-scaleik(9)*mscale(kk)))
     &           + qrri*qrrk*rr11
     &            *(bn(5)-(1.0d0-scaleik(11)*mscale(kk)))

               term1 =-nuck*rr3*(bn(1)-(1.0d0-mscale(kk)*scalei(3))) 
     &            - qk*rr3*(bn(1)-(1.0d0-mscale(kk)*scaleik(3))) 
     &            + drk*rr5*(bn(2)-(1.0d0-mscale(kk)*scaleik(5))) 
     &            - qrrk*rr7*(bn(3)-(1.0d0-mscale(kk)*scaleik(7)))
               term2 = nuci*rr3*(bn(1)-(1.0d0-mscale(kk)*scalek(3))) 
     &            + qi*rr3*(bn(1)-(1.0d0-mscale(kk)*scaleik(3)))  
     &            + dri*rr5*(bn(2)-(1.0d0-mscale(kk)*scaleik(5))) 
     &            + qrri*rr7*(bn(3)-(1.0d0-mscale(kk)*scaleik(7)))
               term3=2.0d0*rr5*(bn(2)-(1.0d0-mscale(kk)*scaleik(5))) 
               term4=2.0d0*(-nuck*rr5
     &             *(bn(2)-(1.0d0-mscale(kk)*scalei(5)))
     &            -qk*rr5*(bn(2)-(1.0d0-mscale(kk)*scaleik(5)))
     &            +drk*rr7*(bn(3)-(1.0d0-mscale(kk)*scaleik(7)))
     &            -qrrk*rr9*(bn(4)-(1.0d0-mscale(kk)*scaleik(9))))
               term5 = 2.0d0 *(-nuci*rr5
     &             *(bn(2)-(1.0d0-mscale(kk)*scalek(5)))
     &            -qi*rr5*(bn(2)-(1.0d0-mscale(kk)*scaleik(5)))
     &            -dri*rr7*(bn(3)-(1.0d0-mscale(kk)*scaleik(7)))
     &            -qrri*rr9*(bn(4)-(1.0d0-mscale(kk)*scaleik(9))))
               term6 = 4.0d0 * rr7 
     &             *(bn(3)-(1.0d0-mscale(kk)*scaleik(7))) 

c
c     compute the force components for this interaction
c
               frcx = de*xr + term1*dix + term2*dkx
     &                   + term3*(diqkx-dkqix) + term4*qrix
     &                   + term5*qrkx + term6*(qikrx+qkirx)
               frcy = de*yr + term1*diy + term2*dky
     &                   + term3*(diqky-dkqiy) + term4*qriy
     &                   + term5*qrky + term6*(qikry+qkiry)
               frcz = de*zr + term1*diz + term2*dkz
     &                   + term3*(diqkz-dkqiz) + term4*qriz
     &                   + term5*qrkz + term6*(qikrz+qkirz)
c
c     compute the torque components for this interaction
c
               rr3 =rr3*(bn(1)-(1.0d0-scaleik(3)*mscale(kk)))
               ttmi(1) = -rr3*dikx + term1*dirx + term3*(dqiqkx+dkqixr)
     &                      - term4*qrixr - term6*(qikrxr+qrrx)
               ttmi(2) = -rr3*diky + term1*diry + term3*(dqiqky+dkqiyr)
     &                      - term4*qriyr - term6*(qikryr+qrry)
               ttmi(3) = -rr3*dikz + term1*dirz + term3*(dqiqkz+dkqizr)
     &                      - term4*qrizr - term6*(qikrzr+qrrz)
               ttmk(1) = rr3*dikx + term2*dkrx - term3*(dqiqkx+diqkxr)
     &                      - term5*qrkxr - term6*(qkirxr-qrrx)
               ttmk(2) = rr3*diky + term2*dkry - term3*(dqiqky+diqkyr)
     &                      - term5*qrkyr - term6*(qkiryr-qrry)
               ttmk(3) = rr3*dikz + term2*dkrz - term3*(dqiqkz+diqkzr)
     &                      - term5*qrkzr - term6*(qkirzr-qrrz)
c
c     force and torque components scaled for self-interactions
c
               if (ii .eq. kk) then
                  frcx = 0.5d0 * frcx
                  frcy = 0.5d0 * frcy
                  frcz = 0.5d0 * frcz
                  do j = 1, 3
                     ttmi(j) = 0.5d0 * ttmi(j)
                     ttmk(j) = 0.5d0 * ttmk(j)
                  end do
               end if
c
c     increment force-based gradient and torque on first site
c
               dem(1,ii) = dem(1,ii) + frcx
               dem(2,ii) = dem(2,ii) + frcy
               dem(3,ii) = dem(3,ii) + frcz
               tem(1,i) = tem(1,i) + ttmi(1)
               tem(2,i) = tem(2,i) + ttmi(2)
               tem(3,i) = tem(3,i) + ttmi(3)
c
c     increment force-based gradient and torque on second site
c
               dem(1,kk) = dem(1,kk) - frcx
               dem(2,kk) = dem(2,kk) - frcy
               dem(3,kk) = dem(3,kk) - frcz
               tem(1,k) = tem(1,k) + ttmk(1)
               tem(2,k) = tem(2,k) + ttmk(2)
               tem(3,k) = tem(3,k) + ttmk(3)
c
c     increment the virial due to pairwise Cartesian forces
c
               vxx = -xr * frcx
               vxy = -yr * frcx
               vxz = -zr * frcx
               vyy = -yr * frcy
               vyz = -zr * frcy
               vzz = -zr * frcz
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vxy
               vir(3,1) = vir(3,1) + vxz
               vir(1,2) = vir(1,2) + vxy
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vyz
               vir(1,3) = vir(1,3) + vxz
               vir(2,3) = vir(2,3) + vyz
               vir(3,3) = vir(3,3) + vzz
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
c     resolve site torques then increment forces and virial
c
      do i = 1, npole
         call torque (i,tem(1,i),fix,fiy,fiz,dem)
         ii = ipole(i)
         iaz = zaxis(i)
         iax = xaxis(i)
         iay = yaxis(i)
         if (iaz .eq. 0)  iaz = ii
         if (iax .eq. 0)  iax = ii
         if (iay .eq. 0)  iay = ii
         xiz = x(iaz) - x(ii)
         yiz = y(iaz) - y(ii)
         ziz = z(iaz) - z(ii)
         xix = x(iax) - x(ii)
         yix = y(iax) - y(ii)
         zix = z(iax) - z(ii)
         xiy = x(iay) - x(ii)
         yiy = y(iay) - y(ii)
         ziy = z(iay) - z(ii)
         vxx = xix*fix(1) + xiy*fiy(1) + xiz*fiz(1)
         vxy = 0.5d0 * (yix*fix(1) + yiy*fiy(1) + yiz*fiz(1)
     &                    + xix*fix(2) + xiy*fiy(2) + xiz*fiz(2))
         vxz = 0.5d0 * (zix*fix(1) + ziy*fiy(1) + ziz*fiz(1)
     &                    + xix*fix(3) + xiy*fiy(3) + xiz*fiz(3)) 
         vyy = yix*fix(2) + yiy*fiy(2) + yiz*fiz(2)
         vyz = 0.5d0 * (zix*fix(2) + ziy*fiy(2) + ziz*fiz(2)
     &                    + yix*fix(3) + yiy*fiy(3) + yiz*fiz(3))
         vzz = zix*fix(3) + ziy*fiy(3) + ziz*fiz(3)
         vir(1,1) = vir(1,1) + vxx
         vir(2,1) = vir(2,1) + vxy
         vir(3,1) = vir(3,1) + vxz
         vir(1,2) = vir(1,2) + vxy
         vir(2,2) = vir(2,2) + vyy
         vir(3,2) = vir(3,2) + vyz
         vir(1,3) = vir(1,3) + vxz
         vir(2,3) = vir(2,3) + vyz
         vir(3,3) = vir(3,3) + vzz
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      deallocate (tem)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine empole1d  --  Ewald multipole derivs via list  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "empole1d" calculates the multipole energy and derivatives
c     with respect to Cartesian coordinates using particle mesh Ewald
c     summation and a neighbor list
c
c
      subroutine empole1d
      use sizes
      use atoms
      use boxes
      use chgpot
      use deriv
      use energi
      use ewald
      use math
      use mpole
      use virial
      implicit none
      integer i,j,ii
      real*8 e,f
      real*8 term,fterm
      real*8 cii,dii,qii
      real*8 xd,yd,zd
      real*8 xq,yq,zq
      real*8 xv,yv,zv,vterm
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 xdfield,ydfield
      real*8 zdfield
      real*8 trq(3),frcx(3)
      real*8 frcy(3),frcz(3)
c
c
c     zero out the atomic multipole energy and derivatives
c
      em = 0.0d0
      do i = 1, n
         do j = 1, 3
            dem(j,i) = 0.0d0
         end do
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
      call emreal1d
c
c     compute the reciprocal space part of the Ewald summation
c
      call emrecip1
c
c     compute the Ewald self-energy term over all the atoms
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
            xd = xd + rpole(2,i) + rpole(1,i)*x(ii)
            yd = yd + rpole(3,i) + rpole(1,i)*y(ii)
            zd = zd + rpole(4,i) + rpole(1,i)*z(ii)
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         em = em + term*(xd*xd+yd*yd+zd*zd)
         do i = 1, npole
            ii = ipole(i)
            dem(1,ii) = dem(1,ii) + 2.0d0*term*rpole(1,i)*xd
            dem(2,ii) = dem(2,ii) + 2.0d0*term*rpole(1,i)*yd
            dem(3,ii) = dem(3,ii) + 2.0d0*term*rpole(1,i)*zd
         end do
         xdfield = -2.0d0 * term * xd
         ydfield = -2.0d0 * term * yd
         zdfield = -2.0d0 * term * zd
         do i = 1, npole
            trq(1) = rpole(3,i)*zdfield - rpole(4,i)*ydfield
            trq(2) = rpole(4,i)*xdfield - rpole(2,i)*zdfield
            trq(3) = rpole(2,i)*ydfield - rpole(3,i)*xdfield
            call torque (i,trq,frcx,frcy,frcz,dem)
         end do
c
c     boundary correction to virial due to overall cell dipole
c
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         xq = 0.0d0
         yq = 0.0d0
         zq = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            xd = xd + rpole(2,i)
            yd = yd + rpole(3,i)
            zd = zd + rpole(4,i)
            xq = xq + rpole(1,i)*x(ii)
            yq = yq + rpole(1,i)*y(ii)
            zq = zq + rpole(1,i)*z(ii)
         end do
         xv = xd * xq
         yv = yd * yq
         zv = zd * zq
         vterm = term * (xd*xd + yd*yd + zd*zd + 2.0d0*(xv+yv+zv)
     &                      + xq*xq + yq*yq + zq*zq)
         vir(1,1) = vir(1,1) + 2.0d0*term*(xq*xq+xv) + vterm
         vir(2,1) = vir(2,1) + 2.0d0*term*(xq*yq+xv)
         vir(3,1) = vir(3,1) + 2.0d0*term*(xq*zq+xv)
         vir(1,2) = vir(1,2) + 2.0d0*term*(yq*xq+yv)
         vir(2,2) = vir(2,2) + 2.0d0*term*(yq*yq+yv) + vterm
         vir(3,2) = vir(3,2) + 2.0d0*term*(yq*zq+yv)
         vir(1,3) = vir(1,3) + 2.0d0*term*(zq*xq+zv)
         vir(2,3) = vir(2,3) + 2.0d0*term*(zq*yq+zv)
         vir(3,3) = vir(3,3) + 2.0d0*term*(zq*zq+zv) + vterm
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine emreal1d  --  Ewald real space derivs via list  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "emreal1d" evaluates the real space portion of the Ewald
c     summation energy and gradient due to multipole interactions
c     via a neighbor list
c
c
      subroutine emreal1d
      use sizes
      use atomid 
      use atoms
      use bound
      use polar 
      use chgpot
      use couple
      use deriv
      use energi
      use ewald
      use inter
      use math
      use molcul
      use mplpot
      use mpole
      use neigh
      use shunt
      use virial
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      integer iax,iay,iaz
      real*8 e,efull,de,f
      real*8 bfac,erfc
      real*8 alsq2,alsq2n
      real*8 exp2a,ralpha
      real*8 scalekk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9,rr11
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dikx,diky,dikz
      real*8 dirx,diry,dirz
      real*8 dkrx,dkry,dkrz
      real*8 qrix,qriy,qriz
      real*8 qrkx,qrky,qrkz
      real*8 qrixr,qriyr,qrizr
      real*8 qrkxr,qrkyr,qrkzr
      real*8 qrrx,qrry,qrrz
      real*8 qikrx,qikry,qikrz
      real*8 qkirx,qkiry,qkirz
      real*8 qikrxr,qikryr,qikrzr
      real*8 qkirxr,qkiryr,qkirzr
      real*8 diqkx,diqky,diqkz
      real*8 dkqix,dkqiy,dkqiz
      real*8 diqkxr,diqkyr,diqkzr
      real*8 dkqixr,dkqiyr,dkqizr
      real*8 dqiqkx,dqiqky,dqiqkz
      real*8 dri,drk,qrri,qrrk
      real*8 diqrk,dkqri
      real*8 dik,qik,qrrik
      real*8 term1,term2,term3
      real*8 term4,term5,term6
      real*8 frcx,frcy,frcz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 ttmi(3),ttmk(3)
      real*8 fix(3),fiy(3),fiz(3)
      real*8 bn(0:5)
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
      real*8, allocatable :: tem(:,:)
      character*6 mode
      external erfc
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (tem(3,n))
c
c     initialize connected atom scaling and torque arrays
c
      do i = 1, n
         mscale(i) = 1.0d0
         do j = 1, 3
            tem(j,i) = 0.0d0
         end do
      end do
c
c     maximum order to rr11 (for force)
c
      allocate (scalei(11))
      allocate (scalek(11))
      allocate (scaleik(11))
c
c     maximum order to rr11 (for force)
c
      do i = 1, 11 
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
!$OMP& shared(npole,ipole,x,y,z,rpole,n12,i12,n13,i13,n14,i14,
!$OMP& n15,i15,m2scale,m3scale,m4scale,m5scale,nelst,elst,
!$OMP& atomic,penalpha,type,
!$OMP& use_bounds,f,off2,aewald,molcule,xaxis,yaxis,zaxis)
!$OMP& firstprivate(mscale) shared (em,einter,dem,tem,vir)
!$OMP DO reduction(+:em,einter,dem,tem,vir) schedule(guided)
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
            r2 = xr*xr + yr*yr + zr*zr
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
               rr11 = 9.0d0 * rr9 / r2
c
c     calculate the real space Ewald error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 5
                  bfac = dble(j+j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
               do j = 0, 5
                  bn(j) = f * bn(j)
               end do
c
c     intermediates involving moments and distance separation
c
               dikx = diy*dkz - diz*dky
               diky = diz*dkx - dix*dkz
               dikz = dix*dky - diy*dkx
               dirx = diy*zr - diz*yr
               diry = diz*xr - dix*zr
               dirz = dix*yr - diy*xr
               dkrx = dky*zr - dkz*yr
               dkry = dkz*xr - dkx*zr
               dkrz = dkx*yr - dky*xr
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
               qrixr = qriz*yr - qriy*zr
               qriyr = qrix*zr - qriz*xr
               qrizr = qriy*xr - qrix*yr
               qrkxr = qrkz*yr - qrky*zr
               qrkyr = qrkx*zr - qrkz*xr
               qrkzr = qrky*xr - qrkx*yr
               qrrx = qrky*qriz - qrkz*qriy
               qrry = qrkz*qrix - qrkx*qriz
               qrrz = qrkx*qriy - qrky*qrix
               qikrx = qixx*qrkx + qixy*qrky + qixz*qrkz
               qikry = qixy*qrkx + qiyy*qrky + qiyz*qrkz
               qikrz = qixz*qrkx + qiyz*qrky + qizz*qrkz
               qkirx = qkxx*qrix + qkxy*qriy + qkxz*qriz
               qkiry = qkxy*qrix + qkyy*qriy + qkyz*qriz
               qkirz = qkxz*qrix + qkyz*qriy + qkzz*qriz
               qikrxr = qikrz*yr - qikry*zr
               qikryr = qikrx*zr - qikrz*xr
               qikrzr = qikry*xr - qikrx*yr
               qkirxr = qkirz*yr - qkiry*zr
               qkiryr = qkirx*zr - qkirz*xr
               qkirzr = qkiry*xr - qkirx*yr
               diqkx = dix*qkxx + diy*qkxy + diz*qkxz
               diqky = dix*qkxy + diy*qkyy + diz*qkyz
               diqkz = dix*qkxz + diy*qkyz + diz*qkzz
               dkqix = dkx*qixx + dky*qixy + dkz*qixz
               dkqiy = dkx*qixy + dky*qiyy + dkz*qiyz
               dkqiz = dkx*qixz + dky*qiyz + dkz*qizz
               diqrk = dix*qrkx + diy*qrky + diz*qrkz
               dkqri = dkx*qrix + dky*qriy + dkz*qriz
               diqkxr = diqkz*yr - diqky*zr
               diqkyr = diqkx*zr - diqkz*xr
               diqkzr = diqky*xr - diqkx*yr
               dkqixr = dkqiz*yr - dkqiy*zr
               dkqiyr = dkqix*zr - dkqiz*xr
               dkqizr = dkqiy*xr - dkqix*yr
               dqiqkx = diy*qrkz - diz*qrky + dky*qriz - dkz*qriy
     &                     - 2.0d0*(qixy*qkxz+qiyy*qkyz+qiyz*qkzz
     &                             -qixz*qkxy-qiyz*qkyy-qizz*qkyz)
               dqiqky = diz*qrkx - dix*qrkz + dkz*qrix - dkx*qriz
     &                     - 2.0d0*(qixz*qkxx+qiyz*qkxy+qizz*qkxz
     &                             -qixx*qkxz-qixy*qkyz-qixz*qkzz)
               dqiqkz = dix*qrky - diy*qrkx + dkx*qriy - dky*qrix
     &                     - 2.0d0*(qixx*qkxy+qixy*qkyy+qixz*qkyz
     &                             -qixy*qkxx-qiyy*qkxy-qiyz*qkxz)
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
                  scaleik(11) = 1.0d0 - termi*(1.0d0 + dampi +
     &                 (4.0d0/9.0d0)*dampi**2 +
     &                 (1.0d0/9.0d0)*dampi**3 +
     &                 (1.0d0/63.0d0)*dampi**4 + 
     &                 (1.0d0/945.0d0)*dampi**5)*expdampi - 
     &                 termk*(1.0d0 + dampk +
     &                 (4.0d0/9.0d0)*dampk**2+
     &                 (1.0d0/9.0d0)*dampk**3+
     &                 (1.0d0/63.0d0)*dampk**4 +
     &                 (1.0d0/945.0d0)*dampk**5)*expdampk
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
                  scaleik(11) = 1.0d0-(1.0d0+dampi+0.5d0*dampi**2
     &                 + (1.0d0/6.0d0)*dampi**3 
     &                 + (5.0d0/126.0d0)*dampi**4
     &                 + (2.0d0/315.0d0)*dampi**5
     &                 + (1.0d0/1890.0d0)*dampi**6)*expdampi
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
c     compute the full energy without any Ewald scaling
c
               efull = term1*rr1 + term2*rr3 + term3*rr5
     &                    + term4*rr7 + term5*rr9
               efull = efull * mscale(kk)
               if (molcule(ii) .ne. molcule(kk))
     &            einter = einter + efull
c
c     modify distances to account for Ewald and exclusions
c
               bn(0) = bn(0)/rr1
               bn(1) = bn(1)/rr3
               bn(2) = bn(2)/rr5
               bn(3) = bn(3)/rr7
               bn(4) = bn(4)/rr9
               bn(5) = bn(5)/rr11

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
c     compute the energy contributions for this interaction
c
               e = term1*rr1 + term2*rr3 + term3*rr5
     &                + term4*rr7 + term5*rr9
               em = em + e
c
c     calculate intermediate terms for force and torque
c
               de =  nuci*nuck*rr3*(bn(1)-(1.0d0-mscale(kk))) 
     &           + nuci*qk*rr3*(bn(1)-(1.0d0-scalek(3)*mscale(kk)))
     &           + nuck * qi*rr3*(bn(1)-(1.0d0-scalei(3)*mscale(kk)))
     &           + qi*qk*rr3*(bn(1)-(1.0d0-scaleik(3)*mscale(kk)))
     &           + nuck*dri*rr5*(bn(2)-(1.0d0-scalei(5)*mscale(kk))) 
     &           - nuci*drk*rr5*(bn(2)-(1.0d0-scalek(5)*mscale(kk)))
     &           + (dik + qk*dri - qi*drk)*rr5
     &            *(bn(2)-(1.0d0-scaleik(5)*mscale(kk)))
     &           + nuci*qrrk*rr7*(bn(3)-(1.0d0-scalek(7)*mscale(kk)))
     &           + nuck*qrri*rr7*(bn(3)-(1.0d0-scalei(7)*mscale(kk)))
     &           +(-dri*drk + 2.0d0*(dkqri-diqrk+qik)  
     &           + qi*qrrk + qk*qrri)*rr7
     &            *(bn(3)-(1.0d0-scaleik(7)*mscale(kk)))
     &           +(dri*qrrk-drk*qrri-4.0d0*qrrik)*rr9
     *            *(bn(4)-(1.0d0-scaleik(9)*mscale(kk)))
     &           + qrri*qrrk*rr11
     &            *(bn(5)-(1.0d0-scaleik(11)*mscale(kk)))

               term1 =-nuck*rr3*(bn(1)-(1.0d0-mscale(kk)*scalei(3))) 
     &            - qk*rr3*(bn(1)-(1.0d0-mscale(kk)*scaleik(3))) 
     &            + drk*rr5*(bn(2)-(1.0d0-mscale(kk)*scaleik(5))) 
     &            - qrrk*rr7*(bn(3)-(1.0d0-mscale(kk)*scaleik(7)))
               term2 = nuci*rr3*(bn(1)-(1.0d0-mscale(kk)*scalek(3))) 
     &            + qi*rr3*(bn(1)-(1.0d0-mscale(kk)*scaleik(3)))  
     &            + dri*rr5*(bn(2)-(1.0d0-mscale(kk)*scaleik(5))) 
     &            + qrri*rr7*(bn(3)-(1.0d0-mscale(kk)*scaleik(7)))
               term3=2.0d0*rr5*(bn(2)-(1.0d0-mscale(kk)*scaleik(5))) 
               term4=2.0d0*(-nuck*rr5
     &             *(bn(2)-(1.0d0-mscale(kk)*scalei(5)))
     &            -qk*rr5*(bn(2)-(1.0d0-mscale(kk)*scaleik(5)))
     &            +drk*rr7*(bn(3)-(1.0d0-mscale(kk)*scaleik(7)))
     &            -qrrk*rr9*(bn(4)-(1.0d0-mscale(kk)*scaleik(9))))
               term5 = 2.0d0 *(-nuci*rr5
     &             *(bn(2)-(1.0d0-mscale(kk)*scalek(5)))
     &            -qi*rr5*(bn(2)-(1.0d0-mscale(kk)*scaleik(5)))
     &            -dri*rr7*(bn(3)-(1.0d0-mscale(kk)*scaleik(7)))
     &            -qrri*rr9*(bn(4)-(1.0d0-mscale(kk)*scaleik(9))))
               term6 = 4.0d0 * rr7 
     &             *(bn(3)-(1.0d0-mscale(kk)*scaleik(7))) 

c
c     compute the force components for this interaction
c
               frcx = de*xr + term1*dix + term2*dkx
     &                   + term3*(diqkx-dkqix) + term4*qrix
     &                   + term5*qrkx + term6*(qikrx+qkirx)
               frcy = de*yr + term1*diy + term2*dky
     &                   + term3*(diqky-dkqiy) + term4*qriy
     &                   + term5*qrky + term6*(qikry+qkiry)
               frcz = de*zr + term1*diz + term2*dkz
     &                   + term3*(diqkz-dkqiz) + term4*qriz
     &                   + term5*qrkz + term6*(qikrz+qkirz)
c
c     compute the torque components for this interaction
c
               rr3 =rr3*(bn(1)-(1.0d0-scaleik(3)*mscale(kk)))
 
               ttmi(1) = -rr3*dikx + term1*dirx + term3*(dqiqkx+dkqixr)
     &                      - term4*qrixr - term6*(qikrxr+qrrx)
               ttmi(2) = -rr3*diky + term1*diry + term3*(dqiqky+dkqiyr)
     &                      - term4*qriyr - term6*(qikryr+qrry)
               ttmi(3) = -rr3*dikz + term1*dirz + term3*(dqiqkz+dkqizr)
     &                      - term4*qrizr - term6*(qikrzr+qrrz)
               ttmk(1) = rr3*dikx + term2*dkrx - term3*(dqiqkx+diqkxr)
     &                      - term5*qrkxr - term6*(qkirxr-qrrx)
               ttmk(2) = rr3*diky + term2*dkry - term3*(dqiqky+diqkyr)
     &                      - term5*qrkyr - term6*(qkiryr-qrry)
               ttmk(3) = rr3*dikz + term2*dkrz - term3*(dqiqkz+diqkzr)
     &                      - term5*qrkzr - term6*(qkirzr-qrrz)
c
c     increment force-based gradient and torque on first site
c
               dem(1,ii) = dem(1,ii) + frcx
               dem(2,ii) = dem(2,ii) + frcy
               dem(3,ii) = dem(3,ii) + frcz
               tem(1,i) = tem(1,i) + ttmi(1)
               tem(2,i) = tem(2,i) + ttmi(2)
               tem(3,i) = tem(3,i) + ttmi(3)
c
c     increment force-based gradient and torque on second site
c
               dem(1,kk) = dem(1,kk) - frcx
               dem(2,kk) = dem(2,kk) - frcy
               dem(3,kk) = dem(3,kk) - frcz
               tem(1,k) = tem(1,k) + ttmk(1)
               tem(2,k) = tem(2,k) + ttmk(2)
               tem(3,k) = tem(3,k) + ttmk(3)
c
c     increment the virial due to pairwise Cartesian forces
c
               vxx = -xr * frcx
               vxy = -yr * frcx
               vxz = -zr * frcx
               vyy = -yr * frcy
               vyz = -zr * frcy
               vzz = -zr * frcz
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vxy
               vir(3,1) = vir(3,1) + vxz
               vir(1,2) = vir(1,2) + vxy
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vyz
               vir(1,3) = vir(1,3) + vxz
               vir(2,3) = vir(2,3) + vyz
               vir(3,3) = vir(3,3) + vzz
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
!$OMP DO reduction(+:dem,vir) schedule(guided)
c
c     resolve site torques then increment forces and virial
c
      do i = 1, npole
         call torque (i,tem(1,i),fix,fiy,fiz,dem)
         ii = ipole(i)
         iaz = zaxis(i)
         iax = xaxis(i)
         iay = yaxis(i)
         if (iaz .eq. 0)  iaz = ii
         if (iax .eq. 0)  iax = ii
         if (iay .eq. 0)  iay = ii
         xiz = x(iaz) - x(ii)
         yiz = y(iaz) - y(ii)
         ziz = z(iaz) - z(ii)
         xix = x(iax) - x(ii)
         yix = y(iax) - y(ii)
         zix = z(iax) - z(ii)
         xiy = x(iay) - x(ii)
         yiy = y(iay) - y(ii)
         ziy = z(iay) - z(ii)
         vxx = xix*fix(1) + xiy*fiy(1) + xiz*fiz(1)
         vxy = 0.5d0 * (yix*fix(1) + yiy*fiy(1) + yiz*fiz(1)
     &                    + xix*fix(2) + xiy*fiy(2) + xiz*fiz(2))
         vxz = 0.5d0 * (zix*fix(1) + ziy*fiy(1) + ziz*fiz(1)
     &                    + xix*fix(3) + xiy*fiy(3) + xiz*fiz(3)) 
         vyy = yix*fix(2) + yiy*fiy(2) + yiz*fiz(2)
         vyz = 0.5d0 * (zix*fix(2) + ziy*fiy(2) + ziz*fiz(2)
     &                    + yix*fix(3) + yiy*fiy(3) + yiz*fiz(3))
         vzz = zix*fix(3) + ziy*fiy(3) + ziz*fiz(3)
         vir(1,1) = vir(1,1) + vxx
         vir(2,1) = vir(2,1) + vxy
         vir(3,1) = vir(3,1) + vxz
         vir(1,2) = vir(1,2) + vxy
         vir(2,2) = vir(2,2) + vyy
         vir(3,2) = vir(3,2) + vyz
         vir(1,3) = vir(1,3) + vxz
         vir(2,3) = vir(2,3) + vyz
         vir(3,3) = vir(3,3) + vzz
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
      deallocate (tem)
      return
      end
c
c
c     ####################################################################
c     ##                                                                ##
c     ##  subroutine emrecip1  --  PME recip multipole energy & derivs  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     "emrecip1" evaluates the reciprocal space portion of the particle
c     mesh Ewald summation energy and gradient due to multipoles
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
      subroutine emrecip1
      use sizes
      use atoms
      use bound
      use boxes
      use chgpot
      use deriv
      use energi
      use ewald
      use math
      use mpole
      use mrecip
      use pme
      use virial
      implicit none
      integer i,j,k,ii
      integer k1,k2,k3
      integer m1,m2,m3
      integer iax,iay,iaz
      integer ntot,nff
      integer nf1,nf2,nf3
      integer deriv1(10)
      integer deriv2(10)
      integer deriv3(10)
      real*8 e,eterm,f
      real*8 r1,r2,r3
      real*8 h1,h2,h3
      real*8 f1,f2,f3
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 volterm,denom
      real*8 hsq,expterm
      real*8 term,pterm
      real*8 vterm,struc2
      real*8 trq(3),fix(3)
      real*8 fiy(3),fiz(3)
c
c     indices into the electrostatic field array
c
      data deriv1  / 2, 5,  8,  9, 11, 16, 18, 14, 15, 20 /
      data deriv2  / 3, 8,  6, 10, 14, 12, 19, 16, 20, 17 /
      data deriv3  / 4, 9, 10,  7, 15, 17, 13, 20, 18, 19 /
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
c     zero out the temporary virial accumulation variables
c
      vxx = 0.0d0
      vxy = 0.0d0
      vxz = 0.0d0
      vyy = 0.0d0
      vyz = 0.0d0
      vzz = 0.0d0
c
c     copy multipole moments and coordinates to local storage
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
c     compute the arrays of B-spline coefficients
c
      call bspline_fill
      call table_fill
c
c     assign permanent multipoles to PME grid and perform
c     the 3-D FFT forward transformation
c
      call cmp_to_fmp (cmp,fmp)
      call grid_mpole (fmp)
      call fftfront
c
c     initialize variables required for the scalar summation
c
      ntot = nfft1 * nfft2 * nfft3
      pterm = (pi/aewald)**2
      volterm = pi * volbox
      nff = nfft1 * nfft2
      nf1 = (nfft1+1) / 2
      nf2 = (nfft2+1) / 2
      nf3 = (nfft3+1) / 2
c
c     make the scalar summation over reciprocal lattice
c
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
            struc2 = qgrid(1,k1,k2,k3)**2 + qgrid(2,k1,k2,k3)**2
            eterm = 0.5d0 * electric * expterm * struc2
            vterm = (2.0d0/hsq) * (1.0d0-term) * eterm
            vxx = vxx + h1*h1*vterm - eterm
            vxy = vxy + h1*h2*vterm
            vxz = vxz + h1*h3*vterm
            vyy = vyy + h2*h2*vterm - eterm
            vyz = vyz + h2*h3*vterm
            vzz = vzz + h3*h3*vterm - eterm
         end if
         qfac(k1,k2,k3) = expterm
      end do
c
c     save the virial for use in polarization computation
c
      vmxx = vxx
      vmxy = vxy
      vmxz = vxz
      vmyy = vyy
      vmyz = vyz
      vmzz = vzz
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
c     complete the transformation of the PME grid
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
      do i = 1, npole
         do j = 1, 20
            fphi(j,i) = f * fphi(j,i)
         end do
      end do
      call fphi_to_cphi (fphi,cphi)
c
c     increment the permanent multipole energy and gradient
c
      e = 0.0d0
      do i = 1, npole
         f1 = 0.0d0
         f2 = 0.0d0
         f3 = 0.0d0
         do k = 1, 10
            e = e + fmp(k,i)*fphi(k,i)
            f1 = f1 + fmp(k,i)*fphi(deriv1(k),i)
            f2 = f2 + fmp(k,i)*fphi(deriv2(k),i)
            f3 = f3 + fmp(k,i)*fphi(deriv3(k),i)
         end do
         f1 = dble(nfft1) * f1
         f2 = dble(nfft2) * f2
         f3 = dble(nfft3) * f3
         h1 = recip(1,1)*f1 + recip(1,2)*f2 + recip(1,3)*f3
         h2 = recip(2,1)*f1 + recip(2,2)*f2 + recip(2,3)*f3
         h3 = recip(3,1)*f1 + recip(3,2)*f2 + recip(3,3)*f3
         ii = ipole(i)
         dem(1,ii) = dem(1,ii) + h1
         dem(2,ii) = dem(2,ii) + h2
         dem(3,ii) = dem(3,ii) + h3
      end do
      e = 0.5d0 * e
      em = em + e
c
c     increment the permanent multipole virial contributions
c
      do i = 1, npole
         vxx = vxx - cmp(2,i)*cphi(2,i) - 2.0d0*cmp(5,i)*cphi(5,i)
     &            - cmp(8,i)*cphi(8,i) - cmp(9,i)*cphi(9,i)
         vxy = vxy - 0.5d0*(cmp(3,i)*cphi(2,i)+cmp(2,i)*cphi(3,i))
     &            - (cmp(5,i)+cmp(6,i))*cphi(8,i)
     &            - 0.5d0*cmp(8,i)*(cphi(5,i)+cphi(6,i))
     &            - 0.5d0*(cmp(9,i)*cphi(10,i)+cmp(10,i)*cphi(9,i))
         vxz = vxz - 0.5d0*(cmp(4,i)*cphi(2,i)+cmp(2,i)*cphi(4,i))
     &            - (cmp(5,i)+cmp(7,i))*cphi(9,i)
     &            - 0.5d0*cmp(9,i)*(cphi(5,i)+cphi(7,i))
     &            - 0.5d0*(cmp(8,i)*cphi(10,i)+cmp(10,i)*cphi(8,i))
         vyy = vyy - cmp(3,i)*cphi(3,i) - 2.0d0*cmp(6,i)*cphi(6,i)
     &            - cmp(8,i)*cphi(8,i) - cmp(10,i)*cphi(10,i)
         vyz = vyz - 0.5d0*(cmp(4,i)*cphi(3,i)+cmp(3,i)*cphi(4,i))
     &            - (cmp(6,i)+cmp(7,i))*cphi(10,i)
     &            - 0.5d0*cmp(10,i)*(cphi(6,i)+cphi(7,i))
     &            - 0.5d0*(cmp(8,i)*cphi(9,i)+cmp(9,i)*cphi(8,i))
         vzz = vzz - cmp(4,i)*cphi(4,i) - 2.0d0*cmp(7,i)*cphi(7,i)
     &            - cmp(9,i)*cphi(9,i) - cmp(10,i)*cphi(10,i)
      end do
c
c     resolve site torques then increment forces and virial
c
      do i = 1, npole
         trq(1) = cmp(4,i)*cphi(3,i) - cmp(3,i)*cphi(4,i)
     &               + 2.0d0*(cmp(7,i)-cmp(6,i))*cphi(10,i)
     &               + cmp(9,i)*cphi(8,i) + cmp(10,i)*cphi(6,i)
     &               - cmp(8,i)*cphi(9,i) - cmp(10,i)*cphi(7,i)
         trq(2) = cmp(2,i)*cphi(4,i) - cmp(4,i)*cphi(2,i)
     &               + 2.0d0*(cmp(5,i)-cmp(7,i))*cphi(9,i)
     &               + cmp(8,i)*cphi(10,i) + cmp(9,i)*cphi(7,i)
     &               - cmp(9,i)*cphi(5,i) - cmp(10,i)*cphi(8,i)
         trq(3) = cmp(3,i)*cphi(2,i) - cmp(2,i)*cphi(3,i)
     &               + 2.0d0*(cmp(6,i)-cmp(5,i))*cphi(8,i)
     &               + cmp(8,i)*cphi(5,i) + cmp(10,i)*cphi(9,i)
     &               - cmp(8,i)*cphi(6,i) - cmp(9,i)*cphi(10,i)
         call torque (i,trq,fix,fiy,fiz,dem)
         ii = ipole(i)
         iaz = zaxis(i)
         iax = xaxis(i)
         iay = yaxis(i)
         if (iaz .eq. 0)  iaz = ii
         if (iax .eq. 0)  iax = ii
         if (iay .eq. 0)  iay = ii
         xiz = x(iaz) - x(ii)
         yiz = y(iaz) - y(ii)
         ziz = z(iaz) - z(ii)
         xix = x(iax) - x(ii)
         yix = y(iax) - y(ii)
         zix = z(iax) - z(ii)
         xiy = x(iay) - x(ii)
         yiy = y(iay) - y(ii)
         ziy = z(iay) - z(ii)
         vxx = vxx + xix*fix(1) + xiy*fiy(1) + xiz*fiz(1)
         vxy = vxy + 0.5d0*(yix*fix(1) + yiy*fiy(1) + yiz*fiz(1)
     &                    + xix*fix(2) + xiy*fiy(2) + xiz*fiz(2))
         vxz = vxz + 0.5d0*(zix*fix(1) + ziy*fiy(1) + ziz*fiz(1)
     &                    + xix*fix(3) + xiy*fiy(3) + xiz*fiz(3)) 
         vyy = vyy + yix*fix(2) + yiy*fiy(2) + yiz*fiz(2)
         vyz = vyz + 0.5d0*(zix*fix(2) + ziy*fiy(2) + ziz*fiz(2)
     &                    + yix*fix(3) + yiy*fiy(3) + yiz*fiz(3))
         vzz = vzz + zix*fix(3) + ziy*fiy(3) + ziz*fiz(3)
      end do
c
c     increment the total internal virial tensor components
c
      vir(1,1) = vir(1,1) + vxx
      vir(2,1) = vir(2,1) + vxy
      vir(3,1) = vir(3,1) + vxz
      vir(1,2) = vir(1,2) + vxy
      vir(2,2) = vir(2,2) + vyy
      vir(3,2) = vir(3,2) + vyz
      vir(1,3) = vir(1,3) + vxz
      vir(2,3) = vir(2,3) + vyz
      vir(3,3) = vir(3,3) + vzz
      return
      end
