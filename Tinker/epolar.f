c
c
c     ##################################################
c     ##  COPYRIGHT (C) 2015  by  Jay William Ponder  ##
c     ##              All Rights Reserved             ##
c     ##################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine epolar  --  induced dipole polarization energy  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "epolar" calculates the polarization energy due to induced
c     dipole interactions
c
c
      subroutine epolar
      use limits
      implicit none
      logical pairwise
c
c
c     choose the method for summing over polarization interactions
c
      pairwise = .false.
      if (pairwise) then
         if (use_ewald) then
            if (use_mlist) then
               call epolar0d
            else
               call epolar0c
            end if
         else
            if (use_mlist) then
               call epolar0b
            else
               call epolar0a
            end if
         end if
      else
         call epolar0e
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine epolar0a  --  double loop polarization energy  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "epolar0a" calculates the induced dipole polarization energy
c     using a double loop, and partitions the energy among atoms
c
c
      subroutine epolar0a
      use sizes
      use atoms
      use bound
      use cell
      use chgpot
      use couple
      use energi
      use mpole
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      implicit none
      integer i,j,k
      integer ii,kk,jcell
      real*8 e,f
      real*8 damp,expdamp
      real*8 pdi,pti,pgamma,pdiri
      real*8 sc3,sc5,sc7
      real*8 psr3,psr5,psr7
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr1
      real*8 rr3,rr5,rr7
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 uix,uiy,uiz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 ukx,uky,ukz
      real*8 dri,drk,uri,urk
      real*8 qrix,qriy,qriz
      real*8 qrkx,qrky,qrkz
      real*8 qrri,qrrk
      real*8 duik,quik
      real*8 term1,term2,term3
      real*8, allocatable :: pscale(:)
      character*6 mode
c
c
c     zero out the total induced dipole polarization energy
c
      ep = 0.0d0
      if (npole .eq. 0)  return
c
c     check the sign of multipole components at chiral sites
c
      if (.not. use_mpole)  call chkpole
c
c     rotate the multipole components into the global frame
c
      if (.not. use_mpole)  call rotpole
c
c     compute the induced dipoles at each polarizable atom
c
      call induce
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         pscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = 0.5d0 * electric / dielec
      mode = 'MPOLE'
      call switch (mode)
c
c     compute the dipole polarization energy component
c
      do i = 1, npole-1
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         pdiri = dirdamp(i)
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
         uix = uind(1,i)
         uiy = uind(2,i)
         uiz = uind(3,i)
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
                if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = p5scale
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
               ukx = uind(1,k)
               uky = uind(2,k)
               ukz = uind(3,k)
c
c     get reciprocal distance terms for this interaction
c
               rr1 = f / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               sc3 = 1.0d0
               sc5 = 1.0d0
               sc7 = 1.0d0
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                 pgamma = min(pdiri, dirdamp(k))
                 damp  = pgamma*sqrt((r/damp)**3)
                 if (damp .gt. -50.0d0) then
                   expdamp= exp(-damp) 
                   sc3 = 1.0d0 - expdamp 
                   sc5 = 1.0d0-(1.0d0 + 0.5d0*damp)*expdamp
                   sc7 = 1.0d0-(1.0d0 + 0.65d0*damp + 
     &                      0.15d0*damp**2)*expdamp
                 end if
               endif
c
c     intermediates involving Thole damping and scale factors
c
               psr3 = rr3 * sc3 * pscale(kk)
               psr5 = rr5 * sc5 * pscale(kk)
               psr7 = rr7 * sc7 * pscale(kk)
c
c     intermediates involving moments and distance separation
c
               dri = dix*xr + diy*yr + diz*zr
               drk = dkx*xr + dky*yr + dkz*zr
               qrix = qixx*xr + qixy*yr + qixz*zr
               qriy = qixy*xr + qiyy*yr + qiyz*zr
               qriz = qixz*xr + qiyz*yr + qizz*zr
               qrkx = qkxx*xr + qkxy*yr + qkxz*zr
               qrky = qkxy*xr + qkyy*yr + qkyz*zr
               qrkz = qkxz*xr + qkyz*yr + qkzz*zr
               qrri = qrix*xr + qriy*yr + qriz*zr
               qrrk = qrkx*xr + qrky*yr + qrkz*zr
               uri = uix*xr + uiy*yr + uiz*zr
               urk = ukx*xr + uky*yr + ukz*zr
               duik = dix*ukx + diy*uky + diz*ukz
     &                   + dkx*uix + dky*uiy + dkz*uiz
               quik = qrix*ukx + qriy*uky + qriz*ukz
     &                   - qrkx*uix - qrky*uiy - qrkz*uiz
c
c     calculate intermediate terms for polarization interaction
c
               term1 = ck*uri - ci*urk + duik
               term2 = 2.0d0*quik - uri*drk - dri*urk
               term3 = uri*qrrk - urk*qrri
c
c     compute the energy contribution for this interaction
c
               e = term1*psr3 + term2*psr5 + term3*psr7
c
c     increment the overall polarization energy components
c
               ep = ep + e
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = 1.0d0
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
            pdi = pdamp(i)
            pti = thole(i)
            pdiri = dirdamp(i)
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
            uix = uind(1,i)
            uiy = uind(2,i)
            uiz = uind(3,i)
            do j = 1, n12(ii)
               pscale(i12(j,ii)) = p2scale
            end do
            do j = 1, n13(ii)
               pscale(i13(j,ii)) = p3scale
            end do
            do j = 1, n14(ii)
               pscale(i14(j,ii)) = p4scale
               do k = 1, np11(ii)
                   if (i14(j,ii) .eq. ip11(k,ii))
     &               pscale(i14(j,ii)) = p4scale * p41scale
               end do
            end do
            do j = 1, n15(ii)
               pscale(i15(j,ii)) = p5scale
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
                  if (use_bounds)  call imager (xr,yr,zr,jcell)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (.not. (use_polymer .and. r2.le.polycut2))
     &               pscale(kk) = 1.0d0
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
                     ukx = uind(1,k)
                     uky = uind(2,k)
                     ukz = uind(3,k)
c
c     get reciprocal distance terms for this interaction
c
                     rr1 = f / r
                     rr3 = rr1 / r2
                     rr5 = 3.0d0 * rr3 / r2
                     rr7 = 5.0d0 * rr5 / r2
                     sc3 = 1.0d0
                     sc5 = 1.0d0
                     sc7 = 1.0d0

                     damp = pdi * pdamp(k)
                     if (damp .ne. 0.0d0) then
                       pgamma = min(pdiri, dirdamp(k))
                       damp  = pgamma*sqrt((r/damp)**3)
                       if (damp .gt. -50.0d0) then
                         expdamp= exp(-damp) 
                         sc3 = 1.0d0 - expdamp 
                         sc5 = 1.0d0-(1.0d0 +0.5d0*damp)*expdamp
                         sc7 = 1.0d0-(1.0d0 + 0.65d0*damp + 
     &                            0.15d0*damp**2)*expdamp
                       end if
                     endif
c
c     intermediates involving Thole damping and scale factors
c
                     psr3 = rr3 * sc3 * pscale(kk)
                     psr5 = rr5 * sc5 * pscale(kk)
                     psr7 = rr7 * sc7 * pscale(kk)
c
c     intermediates involving moments and distance separation
c
                     dri = dix*xr + diy*yr + diz*zr
                     drk = dkx*xr + dky*yr + dkz*zr
                     qrix = qixx*xr + qixy*yr + qixz*zr
                     qriy = qixy*xr + qiyy*yr + qiyz*zr
                     qriz = qixz*xr + qiyz*yr + qizz*zr
                     qrkx = qkxx*xr + qkxy*yr + qkxz*zr
                     qrky = qkxy*xr + qkyy*yr + qkyz*zr
                     qrkz = qkxz*xr + qkyz*yr + qkzz*zr
                     qrri = qrix*xr + qriy*yr + qriz*zr
                     qrrk = qrkx*xr + qrky*yr + qrkz*zr
                     uri = uix*xr + uiy*yr + uiz*zr
                     urk = ukx*xr + uky*yr + ukz*zr
                     duik = dix*ukx + diy*uky + diz*ukz
     &                         + dkx*uix + dky*uiy + dkz*uiz
                     quik = qrix*ukx + qriy*uky + qriz*ukz
     &                         - qrkx*uix - qrky*uiy - qrkz*uiz
c
c     calculate intermediate terms for polarization interaction
c
                     term1 = ck*uri - ci*urk + duik
                     term2 = 2.0d0*quik - uri*drk - dri*urk
                     term3 = uri*qrrk - urk*qrri
c
c     compute the energy contribution for this interaction
c
                     e = term1*psr3 + term2*psr5 + term3*psr7
                     if (ii .eq. kk)  e = 0.5d0 * e
c
c     increment the overall polarization energy components
c
                     ep = ep + e
                  end if
               end do
            end do
c
c     reset exclusion coefficients for connected atoms
c
            do j = 1, n12(ii)
               pscale(i12(j,ii)) = 1.0d0
            end do
            do j = 1, n13(ii)
               pscale(i13(j,ii)) = 1.0d0
            end do
            do j = 1, n14(ii)
               pscale(i14(j,ii)) = 1.0d0
            end do
            do j = 1, n15(ii)
               pscale(i15(j,ii)) = 1.0d0
            end do
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (pscale)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine epolar0b  --  neighbor list polarization energy  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "epolar0b" calculates the induced dipole polarization energy
c     using a neighbor list
c
c
      subroutine epolar0b
      use sizes
      use atoms
      use bound
      use chgpot
      use couple
      use energi
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      real*8 e,f
      real*8 damp,expdamp
      real*8 pdi,pti,pgamma,pdiri
      real*8 sc3,sc5,sc7
      real*8 psr3,psr5,psr7
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr1
      real*8 rr3,rr5,rr7
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 uix,uiy,uiz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 ukx,uky,ukz
      real*8 dri,drk,uri,urk
      real*8 qrix,qriy,qriz
      real*8 qrkx,qrky,qrkz
      real*8 qrri,qrrk
      real*8 duik,quik
      real*8 term1,term2,term3
      real*8, allocatable :: pscale(:)
      character*6 mode
c
c
c     zero out the total polarization energy and partitioning
c
      ep = 0.0d0
      if (npole .eq. 0)  return
c
c     check the sign of multipole components at chiral sites
c
      if (.not. use_mpole)  call chkpole
c
c     rotate the multipole components into the global frame
c
      if (.not. use_mpole)  call rotpole
c
c     compute the induced dipoles at each polarizable atom
c
      call induce
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         pscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = 0.5d0 * electric / dielec
      mode = 'MPOLE'
      call switch (mode)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private)
!$OMP& shared(npole,ipole,pdamp,thole,x,y,z,rpole,uind,n12,
!$OMP& i12,n13,i13,n14,i14,n15,i15,np11,ip11,p2scale,p3scale,
!$OMP& dirdamp, 
!$OMP& p4scale,p5scale,p41scale,nelst,elst,use_bounds,off2,f)
!$OMP& firstprivate(pscale) shared (ep)
!$OMP DO reduction(+:ep) schedule(guided)
c
c     compute the dipole polarization energy component
c
      do i = 1, npole
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         pdiri = dirdamp(i)
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
         uix = uind(1,i)
         uiy = uind(2,i)
         uiz = uind(3,i)
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
                if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = p5scale
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
               ukx = uind(1,k)
               uky = uind(2,k)
               ukz = uind(3,k)
c
c     get reciprocal distance terms for this interaction
c
               rr1 = f / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               sc3 = 1.0d0
               sc5 = 1.0d0
               sc7 = 1.0d0
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                 pgamma = min(pdiri, dirdamp(k))
                 damp  = pgamma*sqrt((r/damp)**3)
                 if (damp .gt. -50.0d0) then
                   expdamp= exp(-damp) 
                   sc3 = 1.0d0 - expdamp 
                   sc5 = 1.0d0-(1.0d0 + 0.5d0*damp)*expdamp
                   sc7 = 1.0d0-(1.0d0 + 0.65d0*damp + 
     &                      0.15d0*damp**2)*expdamp
                 end if
               endif
c
c     intermediates involving Thole damping and scale factors
c
               psr3 = rr3 * sc3 * pscale(kk)
               psr5 = rr5 * sc5 * pscale(kk)
               psr7 = rr7 * sc7 * pscale(kk)
c
c     intermediates involving moments and distance separation
c
               dri = dix*xr + diy*yr + diz*zr
               drk = dkx*xr + dky*yr + dkz*zr
               qrix = qixx*xr + qixy*yr + qixz*zr
               qriy = qixy*xr + qiyy*yr + qiyz*zr
               qriz = qixz*xr + qiyz*yr + qizz*zr
               qrkx = qkxx*xr + qkxy*yr + qkxz*zr
               qrky = qkxy*xr + qkyy*yr + qkyz*zr
               qrkz = qkxz*xr + qkyz*yr + qkzz*zr
               qrri = qrix*xr + qriy*yr + qriz*zr
               qrrk = qrkx*xr + qrky*yr + qrkz*zr
               uri = uix*xr + uiy*yr + uiz*zr
               urk = ukx*xr + uky*yr + ukz*zr
               duik = dix*ukx + diy*uky + diz*ukz
     &                   + dkx*uix + dky*uiy + dkz*uiz
               quik = qrix*ukx + qriy*uky + qriz*ukz
     &                   - qrkx*uix - qrky*uiy - qrkz*uiz
c
c     calculate intermediate terms for polarization interaction
c
               term1 = ck*uri - ci*urk + duik
               term2 = 2.0d0*quik - uri*drk - dri*urk
               term3 = uri*qrrk - urk*qrri
c
c     compute the energy contribution for this interaction
c
               e = term1*psr3 + term2*psr5 + term3*psr7
c
c     increment the overall polarization energy components
c
               ep = ep + e
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = 1.0d0
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
      deallocate (pscale)
      return
      end
c
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine epolar0c  --  Ewald polarization derivs via loop  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "epolar0c" calculates the dipole polarization energy with respect
c     to Cartesian coordinates using particle mesh Ewald summation and
c     a double loop
c
c
      subroutine epolar0c
      use sizes
      use atoms
      use boxes
      use chgpot
      use energi
      use ewald
      use math
      use mpole
      use polar
      use polpot
      use potent
      implicit none
      integer i,ii
      real*8 e,f,term,fterm
      real*8 dix,diy,diz
      real*8 uix,uiy,uiz,uii
      real*8 xd,yd,zd
      real*8 xu,yu,zu
c
c
c     zero out the polarization energy and derivatives
c
      ep = 0.0d0
      if (npole .eq. 0)  return
c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     check the sign of multipole components at chiral sites
c
      if (.not. use_mpole)  call chkpole
c
c     rotate the multipole components into the global frame
c
      if (.not. use_mpole)  call rotpole
c
c     compute the induced dipoles at each polarizable atom
c
      call induce
c
c     compute the real space part of the Ewald summation
c
      call epreal0c
c
c     compute the reciprocal space part of the Ewald summation
c
      call eprecip
c
c     compute the Ewald self-energy term over all the atoms
c
      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / sqrtpi
      do i = 1, npole
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         uix = uind(1,i)
         uiy = uind(2,i)
         uiz = uind(3,i)
         uii = dix*uix + diy*uiy + diz*uiz
         e = fterm * term * uii / 3.0d0
         ep = ep + e
      end do
c
c     compute the cell dipole boundary correction term
c
      if (boundary .eq. 'VACUUM') then
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         xu = 0.0d0
         yu = 0.0d0
         zu = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            xd = xd + rpole(2,i) + rpole(1,i)*x(ii)
            yd = yd + rpole(3,i) + rpole(1,i)*y(ii)
            zd = zd + rpole(4,i) + rpole(1,i)*z(ii)
            xu = xu + uind(1,i)
            yu = yu + uind(2,i)
            zu = zu + uind(3,i)
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         ep = ep + term*(xd*xu+yd*yu+zd*zu)
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine epreal0c  --  real space polar energy via loop  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "epreal0c" calculates the induced dipole polarization energy
c     using particle mesh Ewald summation and a double loop
c
c
      subroutine epreal0c
      use sizes
      use atoms
      use bound
      use cell
      use chgpot
      use couple
      use energi
      use ewald
      use math
      use mpole
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      implicit none
      integer i,j,k
      integer ii,kk,jcell
      real*8 e,f
      real*8 damp,expdamp
      real*8 erfc,bfac
      real*8 alsq2,alsq2n
      real*8 exp2a,ralpha
      real*8 pdi,pti,pgamma,pdiri
      real*8 sc3,sc5,sc7
      real*8 psc3,psc5,psc7
      real*8 psr3,psr5,psr7
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr1
      real*8 rr3,rr5,rr7
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 uix,uiy,uiz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 ukx,uky,ukz
      real*8 dri,drk,uri,urk
      real*8 qrix,qriy,qriz
      real*8 qrkx,qrky,qrkz
      real*8 qrri,qrrk
      real*8 duik,quik
      real*8 term1,term2,term3
      real*8 bn(0:3)
      real*8, allocatable :: pscale(:)
      character*6 mode
      external erfc
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         pscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = 0.5d0 * electric / dielec
      mode = 'MPOLE'
      call switch (mode)
c
c     compute the dipole polarization energy component
c
      do i = 1, npole-1
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         pdiri = dirdamp(i)
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
         uix = uind(1,i)
         uiy = uind(2,i)
         uiz = uind(3,i)
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
                if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = p5scale
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
               ukx = uind(1,k)
               uky = uind(2,k)
               ukz = uind(3,k)
c
c     get reciprocal distance terms for this interaction
c
               rr1 = f / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
c
c     calculate the real space Ewald error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 3
                  bfac = dble(j+j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
               do j = 0, 3
                  bn(j) = f * bn(j)
               end do
c
c    use different damping functions for direct induction !CW
c

               sc3 = 1.0d0
               sc5 = 1.0d0
               sc7 = 1.0d0
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                 pgamma = min(pdiri, dirdamp(k))
                 damp  = pgamma*sqrt((r/damp)**3)
                 if (damp .gt. -50.0d0) then
                   expdamp= exp(-damp) 
                   sc3 = 1.0d0 - expdamp 
                   sc5 = 1.0d0-(1.0d0 + 0.5d0*damp)*expdamp
                   sc7 = 1.0d0-(1.0d0 + 0.65d0*damp + 
     &                      0.15d0*damp**2)*expdamp
                 end if
               endif
c
c     intermediates involving Thole damping and scale factors
c
               psc3 = 1.0d0 - sc3*pscale(kk)
               psc5 = 1.0d0 - sc5*pscale(kk)
               psc7 = 1.0d0 - sc7*pscale(kk)
               psr3 = bn(1) - psc3*rr3
               psr5 = bn(2) - psc5*rr5
               psr7 = bn(3) - psc7*rr7
c
c     intermediates involving moments and distance separation
c
               dri = dix*xr + diy*yr + diz*zr
               drk = dkx*xr + dky*yr + dkz*zr
               qrix = qixx*xr + qixy*yr + qixz*zr
               qriy = qixy*xr + qiyy*yr + qiyz*zr
               qriz = qixz*xr + qiyz*yr + qizz*zr
               qrkx = qkxx*xr + qkxy*yr + qkxz*zr
               qrky = qkxy*xr + qkyy*yr + qkyz*zr
               qrkz = qkxz*xr + qkyz*yr + qkzz*zr
               qrri = qrix*xr + qriy*yr + qriz*zr
               qrrk = qrkx*xr + qrky*yr + qrkz*zr
               uri = uix*xr + uiy*yr + uiz*zr
               urk = ukx*xr + uky*yr + ukz*zr
               duik = dix*ukx + diy*uky + diz*ukz
     &                   + dkx*uix + dky*uiy + dkz*uiz
               quik = qrix*ukx + qriy*uky + qriz*ukz
     &                   - qrkx*uix - qrky*uiy - qrkz*uiz
c
c     calculate intermediate terms for polarization interaction
c
               term1 = ck*uri - ci*urk + duik
               term2 = 2.0d0*quik - uri*drk - dri*urk
               term3 = uri*qrrk - urk*qrri
c
c     compute the energy contribution for this interaction
c
               e = term1*psr3 + term2*psr5 + term3*psr7
c
c     increment the overall polarization energy components
c
               ep = ep + e
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = 1.0d0
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
            pdi = pdamp(i)
            pti = thole(i)
            pdiri = dirdamp(i)
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
            uix = uind(1,i)
            uiy = uind(2,i)
            uiz = uind(3,i)
            do j = 1, n12(ii)
               pscale(i12(j,ii)) = p2scale
            end do
            do j = 1, n13(ii)
               pscale(i13(j,ii)) = p3scale
            end do
            do j = 1, n14(ii)
               pscale(i14(j,ii)) = p4scale
               do k = 1, np11(ii)
                   if (i14(j,ii) .eq. ip11(k,ii))
     &               pscale(i14(j,ii)) = p4scale * p41scale
               end do
            end do
            do j = 1, n15(ii)
               pscale(i15(j,ii)) = p5scale
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
                  if (use_bounds)  call imager (xr,yr,zr,jcell)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (.not. (use_polymer .and. r2.le.polycut2))
     &               pscale(kk) = 1.0d0
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
                     ukx = uind(1,k)
                     uky = uind(2,k)
                     ukz = uind(3,k)
c
c     get reciprocal distance terms for this interaction
c
                     rr1 = f / r
                     rr3 = rr1 / r2
                     rr5 = 3.0d0 * rr3 / r2
                     rr7 = 5.0d0 * rr5 / r2
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
                     do j = 1, 3
                        bfac = dble(j+j-1)
                        alsq2n = alsq2 * alsq2n
                        bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
                     end do
                     do j = 0, 3
                        bn(j) = f * bn(j)
                     end do
                     sc3 = 1.0d0
                     sc5 = 1.0d0
                     sc7 = 1.0d0
                     damp = pdi * pdamp(k)
                     if (damp .ne. 0.0d0) then
                       pgamma = min(pdiri, dirdamp(k))
                       damp  = pgamma*sqrt((r/damp)**3)
                       if (damp .gt. -50.0d0) then
                         expdamp = exp(-damp) 
                         sc3 = 1.0d0 - expdamp 
                         sc5 = 1.0d0-(1.0d0 + 0.5d0*damp)*expdamp
                         sc7 = 1.0d0-(1.0d0 + 0.65d0*damp + 
     &                            0.15d0*damp**2)*expdamp
                       end if
                     endif
c
c     intermediates involving Thole damping and scale factors
c
                     psc3 = 1.0d0 - sc3*pscale(kk)
                     psc5 = 1.0d0 - sc5*pscale(kk)
                     psc7 = 1.0d0 - sc7*pscale(kk)
                     psr3 = bn(1) - psc3*rr3
                     psr5 = bn(2) - psc5*rr5
                     psr7 = bn(3) - psc7*rr7
c
c     intermediates involving moments and distance separation
c
                     dri = dix*xr + diy*yr + diz*zr
                     drk = dkx*xr + dky*yr + dkz*zr
                     qrix = qixx*xr + qixy*yr + qixz*zr
                     qriy = qixy*xr + qiyy*yr + qiyz*zr
                     qriz = qixz*xr + qiyz*yr + qizz*zr
                     qrkx = qkxx*xr + qkxy*yr + qkxz*zr
                     qrky = qkxy*xr + qkyy*yr + qkyz*zr
                     qrkz = qkxz*xr + qkyz*yr + qkzz*zr
                     qrri = qrix*xr + qriy*yr + qriz*zr
                     qrrk = qrkx*xr + qrky*yr + qrkz*zr
                     uri = uix*xr + uiy*yr + uiz*zr
                     urk = ukx*xr + uky*yr + ukz*zr
                     duik = dix*ukx + diy*uky + diz*ukz
     &                         + dkx*uix + dky*uiy + dkz*uiz
                     quik = qrix*ukx + qriy*uky + qriz*ukz
     &                         - qrkx*uix - qrky*uiy - qrkz*uiz
c
c     calculate intermediate terms for polarization interaction
c
                     term1 = ck*uri - ci*urk + duik
                     term2 = 2.0d0*quik - uri*drk - dri*urk
                     term3 = uri*qrrk - urk*qrri
c
c     compute the energy contribution for this interaction
c
                     e = term1*psr3 + term2*psr5 + term3*psr7
                     if (ii .eq. kk)  e = 0.5d0 * e
c
c     increment the overall polarization energy components
c
                     ep = ep + e
                  end if
               end do
            end do
c
c     reset exclusion coefficients for connected atoms
c
            do j = 1, n12(ii)
               pscale(i12(j,ii)) = 1.0d0
            end do
            do j = 1, n13(ii)
               pscale(i13(j,ii)) = 1.0d0
            end do
            do j = 1, n14(ii)
               pscale(i14(j,ii)) = 1.0d0
            end do
            do j = 1, n15(ii)
               pscale(i15(j,ii)) = 1.0d0
            end do
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (pscale)
      return
      end
c
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine epolar0d  --  Ewald polarization derivs via list  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "epolar0d" calculates the dipole polarization energy with respect
c     to Cartesian coordinates using particle mesh Ewald summation and
c     a neighbor list
c
c
      subroutine epolar0d
      use sizes
      use atoms
      use boxes
      use chgpot
      use energi
      use ewald
      use math
      use mpole
      use polar
      use polpot
      use potent
      implicit none
      integer i,ii
      real*8 e,f,term,fterm
      real*8 dix,diy,diz
      real*8 uix,uiy,uiz,uii
      real*8 xd,yd,zd
      real*8 xu,yu,zu
c
c
c     zero out the polarization energy and derivatives
c
      ep = 0.0d0
      if (npole .eq. 0)  return
c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     check the sign of multipole components at chiral sites
c
      if (.not. use_mpole)  call chkpole
c
c     rotate the multipole components into the global frame
c
      if (.not. use_mpole)  call rotpole
c
c     compute the induced dipoles at each polarizable atom
c
      call induce
c
c     compute the real space part of the Ewald summation
c
      call epreal0d
c
c     compute the reciprocal space part of the Ewald summation
c
      call eprecip
c
c     compute the Ewald self-energy term over all the atoms
c
      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / sqrtpi
      do i = 1, npole
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         uix = uind(1,i)
         uiy = uind(2,i)
         uiz = uind(3,i)
         uii = dix*uix + diy*uiy + diz*uiz
         e = fterm * term * uii / 3.0d0
         ep = ep + e
      end do
c
c     compute the cell dipole boundary correction term
c
      if (boundary .eq. 'VACUUM') then
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         xu = 0.0d0
         yu = 0.0d0
         zu = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            xd = xd + rpole(2,i) + rpole(1,i)*x(ii)
            yd = yd + rpole(3,i) + rpole(1,i)*y(ii)
            zd = zd + rpole(4,i) + rpole(1,i)*z(ii)
            xu = xu + uind(1,i)
            yu = yu + uind(2,i)
            zu = zu + uind(3,i)
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         ep = ep + term*(xd*xu+yd*yu+zd*zu)
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine epreal0d  --  real space polar energy via list  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "epreal0d" calculates the induced dipole polarization energy
c     using particle mesh Ewald summation and a neighbor list
c
c
      subroutine epreal0d
      use sizes
      use atoms
      use bound
      use chgpot
      use couple
      use energi
      use ewald
      use math
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      real*8 e,f
      real*8 damp,expdamp
      real*8 erfc,bfac
      real*8 alsq2,alsq2n
      real*8 exp2a,ralpha
      real*8 pdi,pti,pgamma,pdiri
      real*8 sc3,sc5,sc7
      real*8 psc3,psc5,psc7
      real*8 psr3,psr5,psr7
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr1
      real*8 rr3,rr5,rr7
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 uix,uiy,uiz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 ukx,uky,ukz
      real*8 dri,drk,uri,urk
      real*8 qrix,qriy,qriz
      real*8 qrkx,qrky,qrkz
      real*8 qrri,qrrk
      real*8 duik,quik
      real*8 term1,term2,term3
      real*8 bn(0:3)
      real*8, allocatable :: pscale(:)
      character*6 mode
      external erfc
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         pscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = 0.5d0 * electric / dielec
      mode = 'MPOLE'
      call switch (mode)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private)
!$OMP& shared(npole,ipole,pdamp,thole,x,y,z,rpole,uind,n12,i12,
!$OMP& n13,i13,n14,i14,n15,i15,np11,ip11,p2scale,p3scale,p4scale,
!$OMP& dirdamp,
!$OMP& p5scale,p41scale,nelst,elst,use_bounds,off2,f,aewald)
!$OMP& firstprivate(pscale) shared (ep)
!$OMP DO reduction(+:ep) schedule(guided)
c
c     compute the dipole polarization energy component
c
      do i = 1, npole
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         pdiri = dirdamp(i)
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
         uix = uind(1,i)
         uiy = uind(2,i)
         uiz = uind(3,i)
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
                if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = p5scale
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
               ukx = uind(1,k)
               uky = uind(2,k)
               ukz = uind(3,k)
c
c     get reciprocal distance terms for this interaction
c
               rr1 = f / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
c
c     calculate the real space Ewald error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 3
                  bfac = dble(j+j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
               do j = 0, 3
                  bn(j) = f * bn(j)
               end do
               sc3 = 1.0d0
               sc5 = 1.0d0
               sc7 = 1.0d0
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                 pgamma = min(pdiri, dirdamp(k))
                 damp  = pgamma*sqrt((r/damp)**3)
                 if (damp .gt. -50.0d0) then
                   expdamp= exp(-damp) 
                   sc3 = 1.0d0 - expdamp 
                   sc5 = 1.0d0-(1.0d0 + 0.5d0*damp)*expdamp
                   sc7 = 1.0d0-(1.0d0 + 0.65d0*damp + 
     &                      0.15d0*damp**2)*expdamp
                 end if
               endif
c
c     intermediates involving Thole damping and scale factors
c
               psc3 = 1.0d0 - sc3*pscale(kk)
               psc5 = 1.0d0 - sc5*pscale(kk)
               psc7 = 1.0d0 - sc7*pscale(kk)
               psr3 = bn(1) - psc3*rr3
               psr5 = bn(2) - psc5*rr5
               psr7 = bn(3) - psc7*rr7
c
c     intermediates involving moments and distance separation
c
               dri = dix*xr + diy*yr + diz*zr
               drk = dkx*xr + dky*yr + dkz*zr
               qrix = qixx*xr + qixy*yr + qixz*zr
               qriy = qixy*xr + qiyy*yr + qiyz*zr
               qriz = qixz*xr + qiyz*yr + qizz*zr
               qrkx = qkxx*xr + qkxy*yr + qkxz*zr
               qrky = qkxy*xr + qkyy*yr + qkyz*zr
               qrkz = qkxz*xr + qkyz*yr + qkzz*zr
               qrri = qrix*xr + qriy*yr + qriz*zr
               qrrk = qrkx*xr + qrky*yr + qrkz*zr
               uri = uix*xr + uiy*yr + uiz*zr
               urk = ukx*xr + uky*yr + ukz*zr
               duik = dix*ukx + diy*uky + diz*ukz
     &                   + dkx*uix + dky*uiy + dkz*uiz
               quik = qrix*ukx + qriy*uky + qriz*ukz
     &                   - qrkx*uix - qrky*uiy - qrkz*uiz
c
c     calculate intermediate terms for polarization interaction
c
               term1 = ck*uri - ci*urk + duik
               term2 = 2.0d0*quik - uri*drk - dri*urk
               term3 = uri*qrrk - urk*qrri
c
c     compute the energy contribution for this interaction
c
               e = term1*psr3 + term2*psr5 + term3*psr7
c
c     increment the overall polarization energy components
c
               ep = ep + e
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = 1.0d0
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
      deallocate (pscale)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine epolar0e  --  single-loop polarization energy  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "epolar0e" calculates the dipole polarizability interaction
c     from the induced dipoles times the electric field
c
c
      subroutine epolar0e
      use sizes
      use atoms
      use boxes
      use chgpot
      use energi
      use ewald
      use limits
      use math
      use mpole
      use polar
      use potent
      use units
      implicit none
      integer i,j,ii
      real*8 e,f,fi,term
      real*8 xd,yd,zd
      real*8 xu,yu,zu
      real*8 dix,diy,diz
      real*8 uix,uiy,uiz
c
c
c     zero out the total polarization energy
c
      ep = 0.0d0
      if (npole .eq. 0)  return
c
c     check the sign of multipole components at chiral sites
c
      if (.not. use_mpole)  call chkpole
c
c     rotate the multipole components into the global frame
c
      if (.not. use_mpole)  call rotpole
c
c     compute the induced dipoles at each polarizable atom
c
      call induce
c
c     set the energy conversion factor
c
      f = -0.5d0 * electric / dielec
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private)
!$OMP& shared(npole,polarity,f,uind,udirp,ep)
!$OMP DO reduction(+:ep) schedule(guided)
c
c     get polarization energy via induced dipoles times field
c
      do i = 1, npole
         if (polarity(i) .ne. 0.0d0) then
            fi = f / polarity(i)
            e = 0.0d0
            do j = 1, 3
               e = fi * uind(j,i) * udirp(j,i)
               ep = ep + e
            end do
         end if
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     compute the cell dipole boundary correction term
c
      if (use_ewald) then
         if (boundary .eq. 'VACUUM') then
            f = electric / dielec
            xd = 0.0d0
            yd = 0.0d0
            zd = 0.0d0
            xu = 0.0d0
            yu = 0.0d0
            zu = 0.0d0
            do i = 1, npole
               ii = ipole(i)
               dix = rpole(2,i)
               diy = rpole(3,i)
               diz = rpole(4,i)
               uix = uind(1,i)
               uiy = uind(2,i)
               uiz = uind(3,i)
               xd = xd + dix + rpole(1,i)*x(ii)
               yd = yd + diy + rpole(1,i)*y(ii)
               zd = zd + diz + rpole(1,i)*z(ii)
               xu = xu + uix
               yu = yu + uiy
               zu = zu + uiz
            end do
            term = (2.0d0/3.0d0) * f * (pi/volbox)
            e = term * (xd*xu+yd*yu+zd*zu)
            ep = ep + e
         end if
      end if
      return
      end
c
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine eprecip  --  PME recip space polarization energy  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "eprecip" evaluates the reciprocal space portion of particle
c     mesh Ewald summation energy due to dipole polarization
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
      subroutine eprecip
      use sizes
      use atoms
      use bound
      use boxes
      use chgpot
      use energi
      use ewald
      use math
      use mpole
      use mrecip
      use pme
      use polar
      use polpot
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
      real*8 a(3,3),ftc(10,10)
      real*8, allocatable :: fuind(:,:)
      real*8, allocatable :: fuinp(:,:)
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
c     get the fractional to Cartesian transformation matrix
c
      call frac_to_cart (ftc)
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
c     remove scalar sum virial from prior multipole 3-D FFT
c
      if (.not. use_mpole) then
         call bspline_fill
         call table_fill
c
c     assign only the permanent multipoles to the PME grid
c     and perform the 3-D FFT forward transformation
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
         call cmp_to_fmp (cmp,fmp)
         call grid_mpole (fmp)
         call fftfront
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
            end if
            qfac(k1,k2,k3) = expterm
         end do
c
c     account for zeroth grid point for nonperiodic system
c
         qfac(1,1,1) = 0.0d0
         if (.not. use_bounds) then
            expterm = 0.5d0 * pi / xbox
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
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (fuind(3,npole))
      allocate (fuinp(3,npole))
c
c     convert Cartesian induced dipoles to fractional coordinates
c
      do i = 1, 3
         a(1,i) = dble(nfft1) * recip(i,1)
         a(2,i) = dble(nfft2) * recip(i,2)
         a(3,i) = dble(nfft3) * recip(i,3)
      end do
      do i = 1, npole
         do j = 1, 3
            fuind(j,i) = a(j,1)*uind(1,i) + a(j,2)*uind(2,i)
     &                      + a(j,3)*uind(3,i)
            fuinp(j,i) = a(j,1)*uinp(1,i) + a(j,2)*uinp(2,i)
     &                      + a(j,3)*uinp(3,i)
         end do
      end do
c
c     assign PME grid and perform 3-D FFT forward transform
c
      call grid_uind (fuind,fuinp)
      call fftfront
c
c     account for zeroth grid point for nonperiodic system
c
      if (.not. use_bounds) then
         expterm = 0.5d0 * pi / xbox
         struc2 = qgrid(1,1,1,1)**2 + qgrid(2,1,1,1)**2
         e = f * expterm * struc2
         ep = ep + e
      end if
c
c     increment the induced dipole polarization energy
c
      e = 0.0d0
      do i = 1, npole
         do k = 1, 3
            e = e + fuind(k,i)*fphi(k+1,i)
         end do
      end do
      e = 0.5d0 * f * e
      ep = ep + e
c
c     perform deallocation of some local arrays
c
      deallocate (fuind)
      deallocate (fuinp)
      return
      end
