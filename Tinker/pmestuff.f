c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2010 by T. Darden, D. Gohara & Jay W. Ponder  ##
c     ##                      All Rights Reserved                     ##
c     ##################################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  routines below implement various B-spline and coordinate  ##
c     ##  manipulations for particle mesh Ewald summation; spatial  ##
c     ##  grid assignment by David Gohara; modified from original   ##
c     ##  PME code by Thomas Darden, NIEHS, Research Triangle, NC   ##
c     ##                                                            ##
c     ################################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine bspline_fill  --  get PME B-spline coefficients  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "bspline_fill" finds B-spline coefficients and derivatives
c     for PME atomic sites along the fractional coordinate axes
c
c
      subroutine bspline_fill
      use sizes
      use atoms
      use boxes
      use pme
      implicit none
      integer i,ifr
      real*8 xi,yi,zi
      real*8 w,fr,eps
      logical first
      save first
      data first  / .true. /
c
c
c     perform dynamic allocation of some global arrays
c
      if (first) then
         first = .false.
         if (.not. allocated(igrid))  allocate (igrid(3,n))
      end if
c
c     offset used to shift sites off exact lattice bounds
c
      eps = 1.0d-8
c
c     get the B-spline coefficients for each atomic site
c
      do i = 1, n
         xi = x(i)
         yi = y(i)
         zi = z(i)
         w = xi*recip(1,1) + yi*recip(2,1) + zi*recip(3,1)
         fr = dble(nfft1) * (w-dble(anint(w))+0.5d0)
         ifr = int(fr-eps)
         w = fr - dble(ifr)
         igrid(1,i) = ifr - bsorder
         call bsplgen (w,thetai1(1,1,i))
         w = xi*recip(1,2) + yi*recip(2,2) + zi*recip(3,2)
         fr = dble(nfft2) * (w-dble(anint(w))+0.5d0)
         ifr = int(fr-eps)
         w = fr - dble(ifr)
         igrid(2,i) = ifr - bsorder
         call bsplgen (w,thetai2(1,1,i))
         w = xi*recip(1,3) + yi*recip(2,3) + zi*recip(3,3)
         fr = dble(nfft3) * (w-dble(anint(w))+0.5d0)
         ifr = int(fr-eps)
         w = fr - dble(ifr)
         igrid(3,i) = ifr - bsorder
         call bsplgen (w,thetai3(1,1,i))
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine bsplgen  --  B-spline coefficients for an atom  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "bsplgen" gets B-spline coefficients and derivatives for
c     a single PME atomic site along a particular direction
c
c
      subroutine bsplgen (w,thetai)
      use pme
      use potent
      implicit none
      integer i,j,k
      integer level
      real*8 w,denom
      real*8 thetai(4,*)
c
c
c     set B-spline depth for partial charges or multipoles
c
      level = 2
      if (use_mpole .or. use_polar)  level = 4
c
c     initialization to get to 2nd order recursion
c
      bsbuild(2,2) = w
      bsbuild(2,1) = 1.0d0 - w
c
c     perform one pass to get to 3rd order recursion
c
      bsbuild(3,3) = 0.5d0 * w * bsbuild(2,2)
      bsbuild(3,2) = 0.5d0 * ((1.0d0+w)*bsbuild(2,1)
     &                       +(2.0d0-w)*bsbuild(2,2))
      bsbuild(3,1) = 0.5d0 * (1.0d0-w) * bsbuild(2,1)
c
c     compute standard B-spline recursion to desired order
c
      do i = 4, bsorder
         k = i - 1
         denom = 1.0d0 / dble(k)
         bsbuild(i,i) = denom * w * bsbuild(k,k)
         do j = 1, i-2
            bsbuild(i,i-j) = denom * ((w+dble(j))*bsbuild(k,i-j-1)
     &                               +(dble(i-j)-w)*bsbuild(k,i-j))
         end do
         bsbuild(i,1) = denom * (1.0d0-w) * bsbuild(k,1)
      end do
c
c     get coefficients for the B-spline first derivative
c
      k = bsorder - 1
      bsbuild(k,bsorder) = bsbuild(k,bsorder-1)
      do i = bsorder-1, 2, -1
         bsbuild(k,i) = bsbuild(k,i-1) - bsbuild(k,i)
      end do
      bsbuild(k,1) = -bsbuild(k,1)
c
c     get coefficients for the B-spline second derivative
c
      if (level .eq. 4) then
         k = bsorder - 2
         bsbuild(k,bsorder-1) = bsbuild(k,bsorder-2)
         do i = bsorder-2, 2, -1
            bsbuild(k,i) = bsbuild(k,i-1) - bsbuild(k,i)
         end do
         bsbuild(k,1) = -bsbuild(k,1)
         bsbuild(k,bsorder) = bsbuild(k,bsorder-1)
         do i = bsorder-1, 2, -1
            bsbuild(k,i) = bsbuild(k,i-1) - bsbuild(k,i)
         end do
         bsbuild(k,1) = -bsbuild(k,1)
c
c     get coefficients for the B-spline third derivative
c
         k = bsorder - 3
         bsbuild(k,bsorder-2) = bsbuild(k,bsorder-3)
         do i = bsorder-3, 2, -1
            bsbuild(k,i) = bsbuild(k,i-1) - bsbuild(k,i)
         end do
         bsbuild(k,1) = -bsbuild(k,1)
         bsbuild(k,bsorder-1) = bsbuild(k,bsorder-2)
         do i = bsorder-2, 2, -1
            bsbuild(k,i) = bsbuild(k,i-1) - bsbuild(k,i)
         end do
         bsbuild(k,1) = -bsbuild(k,1)
         bsbuild(k,bsorder) = bsbuild(k,bsorder-1)
         do i = bsorder-1, 2, -1
            bsbuild(k,i) = bsbuild(k,i-1) - bsbuild(k,i)
         end do
         bsbuild(k,1) = -bsbuild(k,1)
      end if
c
c     copy coefficients from temporary to permanent storage
c
      do i = 1, bsorder
         do j = 1, level
            thetai(j,i) = bsbuild(bsorder-j+1,i)
         end do
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine table_fill  --  spatial chunks for each site  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "table_fill" constructs an array which stores the spatial
c     regions of the particle mesh Ewald grid with contributions
c     from each electrostatic site
c
c
      subroutine table_fill
      use sizes
      use atoms
      use chunks
      use pme
      implicit none
      integer i,k
      integer cid(3)
      integer nearpt(3)
      integer abound(6)
      integer cbound(6)
      logical negx,negy,negz
      logical posx,posy,posz
      logical midx,midy,midz
c
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private(i,k,cid,nearpt,abound,
!$OMP& cbound,negx,negy,negz,posx,posy,posz,midx,midy,midz)
!$OMP DO
c
c     zero out the PME table marking chunks per site
c
      do k = 1, nchunk
         do i = 1, n
            pmetable(i,k) = 0
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP DO
c
c     loop over sites to find the spatial chunks for each
c
      do i = 1, n
         nearpt(1) = igrid(1,i) + grdoff
         nearpt(2) = igrid(2,i) + grdoff
         nearpt(3) = igrid(3,i) + grdoff
         if (nearpt(1) .lt. 1) then
            nearpt(1) = mod(nearpt(1),nfft1) + nfft1
         else if (nearpt(1) .gt. nfft1) then
            nearpt(1) = mod(nearpt(1),nfft1)
         end if
         if (nearpt(2) .lt. 1) then
            nearpt(2) = mod(nearpt(2),nfft2) + nfft2
         else if (nearpt(2) .gt. nfft2) then
            nearpt(2) = mod(nearpt(2),nfft2)
         end if
         if (nearpt(3) .lt. 1) then
            nearpt(3) = mod(nearpt(3),nfft3) + nfft3
         else if (nearpt(3) .gt. nfft3) then
            nearpt(3) = mod(nearpt(3),nfft3)
         end if
         abound(1) = nearpt(1) - nlpts
         abound(2) = nearpt(1) + nrpts
         abound(3) = nearpt(2) - nlpts
         abound(4) = nearpt(2) + nrpts
         abound(5) = nearpt(3) - nlpts
         abound(6) = nearpt(3) + nrpts
         cid(1) = (nearpt(1)-1)/ngrd1 + 1
         cid(2) = (nearpt(2)-1)/ngrd2 + 1
         cid(3) = (nearpt(3)-1)/ngrd3 + 1
         cbound(1) = (cid(1)-1)*ngrd1 + 1
         cbound(2) = cbound(1) + ngrd1 - 1
         cbound(3) = (cid(2)-1)*ngrd2 + 1
         cbound(4) = cbound(3) + ngrd2 - 1
         cbound(5) = (cid(3)-1)*ngrd3 + 1
         cbound(6) = cbound(5) + ngrd3 - 1
c
c     set and store central chunk where the site is located
c
         k = (cid(3)-1)*nchk1*nchk2 + (cid(2)-1)*nchk1 + cid(1)
         pmetable(i,k) = 1
c
c     flags for atom bounds to left or right of central chunk
c
         negx = (abound(1) .lt. cbound(1))
         negy = (abound(3) .lt. cbound(3))
         negz = (abound(5) .lt. cbound(5))
         posx = (abound(2) .gt. cbound(2))
         posy = (abound(4) .gt. cbound(4))
         posz = (abound(6) .gt. cbound(6))
c
c     flags for atom bounds fully inside the central chunk
c
         midx = (.not.negx .and. .not.posx)
         midy = (.not.negy .and. .not.posy)
         midz = (.not.negz .and. .not.posz)
         if (midx .and. midy .and. midz)  goto 10
c
c     flags for atom bounds that overlap the central chunk
c
         midx = (.not.negx .or. .not.posx)
         midy = (.not.negy .or. .not.posy)
         midz = (.not.negz .or. .not.posz)
c
c     check for overlap with any of the neighboring chunks
c
         if (midx .and. midy .and. negz)  call setchunk (i,cid,0,0,-1)
         if (midx .and. midy .and. posz)  call setchunk (i,cid,0,0,1)
         if (midx .and. negy .and. midz)  call setchunk (i,cid,0,-1,0)
         if (midx .and. posy .and. midz)  call setchunk (i,cid,0,1,0)
         if (negx .and. midy .and. midz)  call setchunk (i,cid,-1,0,0)
         if (posx .and. midy .and. midz)  call setchunk (i,cid,1,0,0)
         if (midx .and. negy .and. negz)  call setchunk (i,cid,0,-1,-1)
         if (midx .and. negy .and. posz)  call setchunk (i,cid,0,-1,1)
         if (midx .and. posy .and. negz)  call setchunk (i,cid,0,1,-1)
         if (midx .and. posy .and. posz)  call setchunk (i,cid,0,1,1)
         if (negx .and. midy .and. negz)  call setchunk (i,cid,-1,0,-1)
         if (negx .and. midy .and. posz)  call setchunk (i,cid,-1,0,1)
         if (posx .and. midy .and. negz)  call setchunk (i,cid,1,0,-1)
         if (posx .and. midy .and. posz)  call setchunk (i,cid,1,0,1)
         if (negx .and. negy .and. midz)  call setchunk (i,cid,-1,-1,0)
         if (negx .and. posy .and. midz)  call setchunk (i,cid,-1,1,0)
         if (posx .and. negy .and. midz)  call setchunk (i,cid,1,-1,0)
         if (posx .and. posy .and. midz)  call setchunk (i,cid,1,1,0)
         if (negx .and. negy .and. negz)  call setchunk (i,cid,-1,-1,-1)
         if (negx .and. negy .and. posz)  call setchunk (i,cid,-1,-1,1)
         if (negx .and. posy .and. negz)  call setchunk (i,cid,-1,1,-1)
         if (posx .and. negy .and. negz)  call setchunk (i,cid,1,-1,-1)
         if (negx .and. posy .and. posz)  call setchunk (i,cid,-1,1,1)
         if (posx .and. negy .and. posz)  call setchunk (i,cid,1,-1,1)
         if (posx .and. posy .and. negz)  call setchunk (i,cid,1,1,-1)
         if (posx .and. posy .and. posz)  call setchunk (i,cid,1,1,1)
   10    continue
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine setchunk  --  site overlaps neighboring chunk  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "setchunk" marks a chunk in the PME spatial table which is
c     overlapped by the B-splines for an electrostatic site
c
c
      subroutine setchunk (i,cid,off1,off2,off3)
      use sizes
      use chunks
      use pme
      implicit none
      integer i,k
      integer off1,off2,off3
      integer cid(3),temp(3)
c
c
c     mark neighboring chunk overlapped by an electrostatic site
c
      temp(1) = cid(1) + off1
      if (temp(1) .lt. 1)  temp(1) = nchk1
      if (temp(1) .gt. nchk1)  temp(1) = 1
      temp(2) = cid(2) + off2
      if (temp(2) .lt. 1)  temp(2) = nchk2
      if (temp(2) .gt. nchk2)  temp(2) = 1
      temp(3) = cid(3) + off3
      if (temp(3) .lt. 1)  temp(3) = nchk3
      if (temp(3) .gt. nchk3)  temp(3) = 1
      k = (temp(3)-1)*nchk1*nchk2 + (temp(2)-1)*nchk1 + temp(1)
      pmetable(i,k) = 1
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine grid_pchg  --  put partial charges on PME grid  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "grid_pchg" places the fractional atomic partial charges onto
c     the particle mesh Ewald grid
c
c
      subroutine grid_pchg
      use sizes
      use atoms
      use charge
      use chunks
      use pme
      implicit none
      integer i,j,k,m
      integer ii,jj,kk
      integer ichk,isite,iatm
      integer offsetx,offsety
      integer offsetz
      integer cid(3)
      integer nearpt(3)
      integer abound(6)
      integer cbound(6)
      real*8 v0,u0,t0
      real*8 term
c
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private(i,j,k,m,ii,jj,kk,ichk,
!$OMP& isite,iatm,cid,nearpt,cbound,abound,offsetx,offsety,
!$OMP& offsetz,v0,u0,term,t0)
!$OMP DO
c
c     zero out the particle mesh Ewald charge grid
c
      do k = 1, nfft3
         do j = 1, nfft2
            do i = 1, nfft1
               qgrid(1,i,j,k) = 0.0d0
               qgrid(2,i,j,k) = 0.0d0
            end do
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP DO
c
c     put the permanent multipole moments onto the grid
c
      do ichk = 1, nchunk
         cid(1) = mod(ichk-1,nchk1)
         cid(2) = mod(((ichk-1-cid(1))/nchk1),nchk2)
         cid(3) = mod((ichk-1)/(nchk1*nchk2),nchk3)
         cbound(1) = cid(1)*ngrd1 + 1
         cbound(2) = cbound(1) + ngrd1 - 1
         cbound(3) = cid(2)*ngrd2 + 1
         cbound(4) = cbound(3) + ngrd2 - 1
         cbound(5) = cid(3)*ngrd3 + 1
         cbound(6) = cbound(5) + ngrd3 - 1
         do isite = 1, nion
            iatm = iion(isite)
            if (pmetable(iatm,ichk) .eq. 1) then
               nearpt(1) = igrid(1,iatm) + grdoff
               nearpt(2) = igrid(2,iatm) + grdoff
               nearpt(3) = igrid(3,iatm) + grdoff
               abound(1) = nearpt(1) - nlpts
               abound(2) = nearpt(1) + nrpts
               abound(3) = nearpt(2) - nlpts
               abound(4) = nearpt(2) + nrpts
               abound(5) = nearpt(3) - nlpts
               abound(6) = nearpt(3) + nrpts
               call adjust (offsetx,nfft1,nchk1,abound(1),
     &                        abound(2),cbound(1),cbound(2))
               call adjust (offsety,nfft2,nchk2,abound(3),
     &                        abound(4),cbound(3),cbound(4))
               call adjust (offsetz,nfft3,nchk3,abound(5),
     &                        abound(6),cbound(5),cbound(6))
               do kk = abound(5), abound(6)
                  k = kk
                  m = k + offsetz
                  if (k .lt. 1)  k = k + nfft3
                  v0 = thetai3(1,m,iatm) * pchg(isite)
                  do jj = abound(3), abound(4)
                     j = jj
                     m = j + offsety
                     if (j .lt. 1)  j = j + nfft2
                     u0 = thetai2(1,m,iatm)
                     term = v0 * u0
                     do ii = abound(1), abound(2)
                        i = ii
                        m = i + offsetx
                        if (i .lt. 1)  i = i + nfft1
                        t0 = thetai1(1,m,iatm)
                        qgrid(1,i,j,k) = qgrid(1,i,j,k) + term*t0
                     end do
                  end do
               end do
            end if
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine grid_mpole  --  put multipoles on PME grid  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "grid_mpole" places the fractional atomic multipoles onto
c     the particle mesh Ewald grid
c
c
      subroutine grid_mpole (fmp)
      use sizes
      use atoms
      use chunks
      use mpole
      use pme
      implicit none
      integer i,j,k,m
      integer ii,jj,kk
      integer ichk,isite,iatm
      integer offsetx,offsety
      integer offsetz
      integer cid(3)
      integer nearpt(3)
      integer abound(6)
      integer cbound(6)
      real*8 v0,u0,t0
      real*8 v1,u1,t1
      real*8 v2,u2,t2
      real*8 term0,term1,term2
      real*8 fmp(10,*)
c
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private(i,j,k,m,ii,jj,kk,ichk,
!$OMP& isite,iatm,cid,nearpt,cbound,abound,offsetx,offsety,
!$OMP& offsetz,v0,v1,v2,u0,u1,u2,term0,term1,term2,t0,t1,t2)
!$OMP DO
c
c     zero out the particle mesh Ewald charge grid
c
      do k = 1, nfft3
         do j = 1, nfft2
            do i = 1, nfft1
               qgrid(1,i,j,k) = 0.0d0
               qgrid(2,i,j,k) = 0.0d0
            end do
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP DO
c
c     put the permanent multipole moments onto the grid
c
      do ichk = 1, nchunk
         cid(1) = mod(ichk-1,nchk1)
         cid(2) = mod(((ichk-1-cid(1))/nchk1),nchk2)
         cid(3) = mod((ichk-1)/(nchk1*nchk2),nchk3)
         cbound(1) = cid(1)*ngrd1 + 1
         cbound(2) = cbound(1) + ngrd1 - 1
         cbound(3) = cid(2)*ngrd2 + 1
         cbound(4) = cbound(3) + ngrd2 - 1
         cbound(5) = cid(3)*ngrd3 + 1
         cbound(6) = cbound(5) + ngrd3 - 1
         do isite = 1, npole
            iatm = ipole(isite)
            if (pmetable(iatm,ichk) .eq. 1) then
               nearpt(1) = igrid(1,iatm) + grdoff
               nearpt(2) = igrid(2,iatm) + grdoff
               nearpt(3) = igrid(3,iatm) + grdoff
               abound(1) = nearpt(1) - nlpts
               abound(2) = nearpt(1) + nrpts
               abound(3) = nearpt(2) - nlpts
               abound(4) = nearpt(2) + nrpts
               abound(5) = nearpt(3) - nlpts
               abound(6) = nearpt(3) + nrpts
               call adjust (offsetx,nfft1,nchk1,abound(1),
     &                        abound(2),cbound(1),cbound(2))
               call adjust (offsety,nfft2,nchk2,abound(3),
     &                        abound(4),cbound(3),cbound(4))
               call adjust (offsetz,nfft3,nchk3,abound(5),
     &                        abound(6),cbound(5),cbound(6))
               do kk = abound(5), abound(6)
                  k = kk
                  m = k + offsetz
                  if (k .lt. 1)  k = k + nfft3
                  v0 = thetai3(1,m,iatm)
                  v1 = thetai3(2,m,iatm)
                  v2 = thetai3(3,m,iatm)
                  do jj = abound(3), abound(4)
                     j = jj
                     m = j + offsety
                     if (j .lt. 1)  j = j + nfft2
                     u0 = thetai2(1,m,iatm)
                     u1 = thetai2(2,m,iatm)
                     u2 = thetai2(3,m,iatm)
                     term0 = fmp(1,isite)*u0*v0 + fmp(3,isite)*u1*v0
     &                     + fmp(4,isite)*u0*v1 + fmp(6,isite)*u2*v0
     &                     + fmp(7,isite)*u0*v2 + fmp(10,isite)*u1*v1
                     term1 = fmp(2,isite)*u0*v0 + fmp(8,isite)*u1*v0
     &                          + fmp(9,isite)*u0*v1
                     term2 = fmp(5,isite) * u0 * v0
                     do ii = abound(1), abound(2)
                        i = ii
                        m = i + offsetx
                        if (i .lt. 1)  i = i + nfft1
                        t0 = thetai1(1,m,iatm)
                        t1 = thetai1(2,m,iatm)
                        t2 = thetai1(3,m,iatm)
                        qgrid(1,i,j,k) = qgrid(1,i,j,k) + term0*t0
     &                                      + term1*t1 + term2*t2
                     end do
                  end do
               end do
            end if
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine grid_uind  --  put induced dipoles on PME grid  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "grid_uind" places the fractional induced dipoles onto the
c     particle mesh Ewald grid
c
c
      subroutine grid_uind (fuind,fuinp)
      use sizes
      use atoms
      use chunks
      use mpole
      use pme
      implicit none
      integer i,j,k,m
      integer ii,jj,kk
      integer ichk,isite,iatm
      integer offsetx,offsety
      integer offsetz
      integer cid(3)
      integer nearpt(3)
      integer abound(6)
      integer cbound(6)
      real*8 v0,u0,t0
      real*8 v1,u1,t1
      real*8 term01,term11
      real*8 term02,term12
      real*8 fuind(3,*)
      real*8 fuinp(3,*)
c
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private(i,j,k,m,ii,jj,kk,ichk,
!$OMP& isite,iatm,cid,nearpt,cbound,abound,offsetx,offsety,
!$OMP& offsetz,v0,v1,u0,u1,term01,term11,term02,term12,t0,t1)
!$OMP DO
c
c     zero out the particle mesh Ewald charge grid
c
      do k = 1, nfft3
         do j = 1, nfft2
            do i = 1, nfft1
               qgrid(1,i,j,k) = 0.0d0
               qgrid(2,i,j,k) = 0.0d0
            end do
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP DO
c
c     put the induced dipole moments onto the grid
c
      do ichk = 1, nchunk
         cid(1) = mod(ichk-1,nchk1)
         cid(2) = mod(((ichk-1-cid(1))/nchk1),nchk2)
         cid(3) = mod((ichk-1)/(nchk1*nchk2),nchk3)
         cbound(1) = cid(1)*ngrd1 + 1
         cbound(2) = cbound(1) + ngrd1 - 1
         cbound(3) = cid(2)*ngrd2 + 1
         cbound(4) = cbound(3) + ngrd2 - 1
         cbound(5) = cid(3)*ngrd3 + 1
         cbound(6) = cbound(5) + ngrd3 - 1
         do isite = 1, npole
            iatm = ipole(isite)
            if (pmetable(iatm,ichk) .eq. 1) then
               nearpt(1) = igrid(1,iatm) + grdoff
               nearpt(2) = igrid(2,iatm) + grdoff
               nearpt(3) = igrid(3,iatm) + grdoff
               abound(1) = nearpt(1) - nlpts
               abound(2) = nearpt(1) + nrpts
               abound(3) = nearpt(2) - nlpts
               abound(4) = nearpt(2) + nrpts
               abound(5) = nearpt(3) - nlpts
               abound(6) = nearpt(3) + nrpts
               call adjust (offsetx,nfft1,nchk1,abound(1),
     &                        abound(2),cbound(1),cbound(2))
               call adjust (offsety,nfft2,nchk2,abound(3),
     &                        abound(4),cbound(3),cbound(4))
               call adjust (offsetz,nfft3,nchk3,abound(5),
     &                        abound(6),cbound(5),cbound(6))
               do kk = abound(5), abound(6)
                  k = kk
                  m = k + offsetz
                  if (k .lt. 1)  k = k + nfft3
                  v0 = thetai3(1,m,iatm)
                  v1 = thetai3(2,m,iatm)
                  do jj = abound(3), abound(4)
                     j = jj
                     m = j + offsety
                     if (j .lt. 1)  j = j + nfft2
                     u0 = thetai2(1,m,iatm)
                     u1 = thetai2(2,m,iatm)
                     term01 = fuind(2,isite)*u1*v0
     &                           + fuind(3,isite)*u0*v1
                     term11 = fuind(1,isite)*u0*v0
                     term02 = fuinp(2,isite)*u1*v0
     &                           + fuinp(3,isite)*u0*v1
                     term12 = fuinp(1,isite)*u0*v0
                     do ii = abound(1), abound(2)
                        i = ii
                        m = i + offsetx
                        if (i .lt. 1)  i = i + nfft1
                        t0 = thetai1(1,m,iatm)
                        t1 = thetai1(2,m,iatm)
                        qgrid(1,i,j,k) = qgrid(1,i,j,k) + term01*t0
     &                                      + term11*t1
                        qgrid(2,i,j,k) = qgrid(2,i,j,k) + term02*t0
     &                                      + term12*t1
                     end do
                  end do
               end do
            end if
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine adjust  --  alter site bounds for the PME grid  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "adjust" modifies site bounds on the PME grid and returns
c     an offset into the B-spline coefficient arrays
c
c
      subroutine adjust (offset,nfft,nchk,amin,amax,cmin,cmax)
      implicit none
      integer offset
      integer nfft,nchk
      integer amin,amax
      integer cmin,cmax
c
c
c     modify grid offset and bounds for site at edge of chunk
c
      offset = 0
      if (nchk .ne. 1) then
         if (amin.lt.cmin .or. amax.gt.cmax) then
            if (amin.lt.1 .or. amax.gt.nfft) then
               if (cmin .eq. 1) then
                  offset = 1 - amin
                  amin = 1
               else if (cmax .eq. nfft) then
                  amax = nfft
                  amin = amin + nfft
               end if
            else
               if (cmin .gt. amin) then
                  offset = cmin - amin
                  amin = cmin
               else
                  amax = cmax
               end if
            end if
         end if
      end if
      offset = offset + 1 - amin
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine fphi_mpole  --  multipole potential from grid  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "fphi_mpole" extracts the permanent multipole potential from
c     the particle mesh Ewald grid
c
c
      subroutine fphi_mpole (fphi)
      use sizes
      use mpole
      use pme
      implicit none
      integer i,j,k
      integer isite,iatm
      integer i0,j0,k0
      integer it1,it2,it3
      integer igrd0,jgrd0,kgrd0
      real*8 v0,v1,v2,v3
      real*8 u0,u1,u2,u3
      real*8 t0,t1,t2,t3,tq
      real*8 tu00,tu10,tu01,tu20,tu11
      real*8 tu02,tu21,tu12,tu30,tu03
      real*8 tuv000,tuv100,tuv010,tuv001
      real*8 tuv200,tuv020,tuv002,tuv110
      real*8 tuv101,tuv011,tuv300,tuv030
      real*8 tuv003,tuv210,tuv201,tuv120
      real*8 tuv021,tuv102,tuv012,tuv111
      real*8 fphi(20,*)
c
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(npole,ipole,igrid,bsorder,
!$OMP& nfft3,thetai3,nfft2,thetai2,nfft1,thetai1,qgrid,fphi)
!$OMP DO
c
c     extract the permanent multipole field at each site
c
      do isite = 1, npole
         iatm = ipole(isite)
         igrd0 = igrid(1,iatm)
         jgrd0 = igrid(2,iatm)
         kgrd0 = igrid(3,iatm)
         tuv000 = 0.0d0
         tuv001 = 0.0d0
         tuv010 = 0.0d0
         tuv100 = 0.0d0
         tuv200 = 0.0d0
         tuv020 = 0.0d0
         tuv002 = 0.0d0
         tuv110 = 0.0d0
         tuv101 = 0.0d0
         tuv011 = 0.0d0
         tuv300 = 0.0d0
         tuv030 = 0.0d0
         tuv003 = 0.0d0
         tuv210 = 0.0d0
         tuv201 = 0.0d0
         tuv120 = 0.0d0
         tuv021 = 0.0d0
         tuv102 = 0.0d0
         tuv012 = 0.0d0
         tuv111 = 0.0d0
         k0 = kgrd0
         do it3 = 1, bsorder
            k0 = k0 + 1
            k = k0 + 1 + (nfft3-isign(nfft3,k0))/2
            v0 = thetai3(1,it3,iatm)
            v1 = thetai3(2,it3,iatm)
            v2 = thetai3(3,it3,iatm)
            v3 = thetai3(4,it3,iatm)
            tu00 = 0.0d0
            tu10 = 0.0d0
            tu01 = 0.0d0
            tu20 = 0.0d0
            tu11 = 0.0d0
            tu02 = 0.0d0
            tu30 = 0.0d0
            tu21 = 0.0d0
            tu12 = 0.0d0
            tu03 = 0.0d0
            j0 = jgrd0
            do it2 = 1, bsorder
               j0 = j0 + 1
               j = j0 + 1 + (nfft2-isign(nfft2,j0))/2
               u0 = thetai2(1,it2,iatm)
               u1 = thetai2(2,it2,iatm)
               u2 = thetai2(3,it2,iatm)
               u3 = thetai2(4,it2,iatm)
               t0 = 0.0d0
               t1 = 0.0d0
               t2 = 0.0d0
               t3 = 0.0d0
               i0 = igrd0
               do it1 = 1, bsorder
                  i0 = i0 + 1
                  i = i0 + 1 + (nfft1-isign(nfft1,i0))/2
                  tq = qgrid(1,i,j,k)
                  t0 = t0 + tq*thetai1(1,it1,iatm)
                  t1 = t1 + tq*thetai1(2,it1,iatm)
                  t2 = t2 + tq*thetai1(3,it1,iatm)
                  t3 = t3 + tq*thetai1(4,it1,iatm)
               end do
               tu00 = tu00 + t0*u0
               tu10 = tu10 + t1*u0
               tu01 = tu01 + t0*u1
               tu20 = tu20 + t2*u0
               tu11 = tu11 + t1*u1
               tu02 = tu02 + t0*u2
               tu30 = tu30 + t3*u0
               tu21 = tu21 + t2*u1
               tu12 = tu12 + t1*u2
               tu03 = tu03 + t0*u3
            end do
            tuv000 = tuv000 + tu00*v0
            tuv100 = tuv100 + tu10*v0
            tuv010 = tuv010 + tu01*v0
            tuv001 = tuv001 + tu00*v1
            tuv200 = tuv200 + tu20*v0
            tuv020 = tuv020 + tu02*v0
            tuv002 = tuv002 + tu00*v2
            tuv110 = tuv110 + tu11*v0
            tuv101 = tuv101 + tu10*v1
            tuv011 = tuv011 + tu01*v1
            tuv300 = tuv300 + tu30*v0
            tuv030 = tuv030 + tu03*v0
            tuv003 = tuv003 + tu00*v3
            tuv210 = tuv210 + tu21*v0
            tuv201 = tuv201 + tu20*v1
            tuv120 = tuv120 + tu12*v0
            tuv021 = tuv021 + tu02*v1
            tuv102 = tuv102 + tu10*v2
            tuv012 = tuv012 + tu01*v2
            tuv111 = tuv111 + tu11*v1
         end do
         fphi(1,isite) = tuv000
         fphi(2,isite) = tuv100
         fphi(3,isite) = tuv010
         fphi(4,isite) = tuv001
         fphi(5,isite) = tuv200
         fphi(6,isite) = tuv020
         fphi(7,isite) = tuv002
         fphi(8,isite) = tuv110
         fphi(9,isite) = tuv101
         fphi(10,isite) = tuv011
         fphi(11,isite) = tuv300
         fphi(12,isite) = tuv030
         fphi(13,isite) = tuv003
         fphi(14,isite) = tuv210
         fphi(15,isite) = tuv201
         fphi(16,isite) = tuv120
         fphi(17,isite) = tuv021
         fphi(18,isite) = tuv102
         fphi(19,isite) = tuv012
         fphi(20,isite) = tuv111
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine fphi_uind  --  induced potential from grid  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "fphi_uind" extracts the induced dipole potential from
c     the particle mesh Ewald grid
c
c
      subroutine fphi_uind (fdip_phi1,fdip_phi2,fdip_sum_phi)
      use sizes
      use mpole
      use pme
      implicit none
      integer i,j,k
      integer isite,iatm
      integer i0,j0,k0
      integer it1,it2,it3
      integer igrd0,jgrd0,kgrd0
      real*8 v0,v1,v2,v3
      real*8 u0,u1,u2,u3
      real*8 t0,t1,t2,t3
      real*8 t0_1,t0_2,t1_1,t1_2
      real*8 t2_1,t2_2,tq_1,tq_2
      real*8 tu00,tu10,tu01,tu20,tu11
      real*8 tu02,tu30,tu21,tu12,tu03
      real*8 tu00_1,tu01_1,tu10_1
      real*8 tu00_2,tu01_2,tu10_2
      real*8 tu20_1,tu11_1,tu02_1
      real*8 tu20_2,tu11_2,tu02_2
      real*8 tuv100_1,tuv010_1,tuv001_1
      real*8 tuv100_2,tuv010_2,tuv001_2
      real*8 tuv200_1,tuv020_1,tuv002_1
      real*8 tuv110_1,tuv101_1,tuv011_1
      real*8 tuv200_2,tuv020_2,tuv002_2
      real*8 tuv110_2,tuv101_2,tuv011_2
      real*8 tuv000,tuv100,tuv010,tuv001
      real*8 tuv200,tuv020,tuv002,tuv110
      real*8 tuv101,tuv011,tuv300,tuv030
      real*8 tuv003,tuv210,tuv201,tuv120
      real*8 tuv021,tuv102,tuv012,tuv111
      real*8 fdip_phi1(10,*)
      real*8 fdip_phi2(10,*)
      real*8 fdip_sum_phi(20,*)
c
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(npole,ipole,
!$OMP& igrid,bsorder,nfft3,thetai3,nfft2,thetai2,nfft1,
!$OMP& thetai1,qgrid,fdip_phi1,fdip_phi2,fdip_sum_phi)
!$OMP DO
c
c     extract the induced dipole field at each site
c
      do isite = 1, npole
         iatm = ipole(isite)
         igrd0 = igrid(1,iatm)
         jgrd0 = igrid(2,iatm)
         kgrd0 = igrid(3,iatm)
         tuv100_1 = 0.0d0
         tuv010_1 = 0.0d0
         tuv001_1 = 0.0d0
         tuv200_1 = 0.0d0
         tuv020_1 = 0.0d0
         tuv002_1 = 0.0d0
         tuv110_1 = 0.0d0
         tuv101_1 = 0.0d0
         tuv011_1 = 0.0d0
         tuv100_2 = 0.0d0
         tuv010_2 = 0.0d0
         tuv001_2 = 0.0d0
         tuv200_2 = 0.0d0
         tuv020_2 = 0.0d0
         tuv002_2 = 0.0d0
         tuv110_2 = 0.0d0
         tuv101_2 = 0.0d0
         tuv011_2 = 0.0d0
         tuv000 = 0.0d0
         tuv001 = 0.0d0
         tuv010 = 0.0d0
         tuv100 = 0.0d0
         tuv200 = 0.0d0
         tuv020 = 0.0d0
         tuv002 = 0.0d0
         tuv110 = 0.0d0
         tuv101 = 0.0d0
         tuv011 = 0.0d0
         tuv300 = 0.0d0
         tuv030 = 0.0d0
         tuv003 = 0.0d0
         tuv210 = 0.0d0
         tuv201 = 0.0d0
         tuv120 = 0.0d0
         tuv021 = 0.0d0
         tuv102 = 0.0d0
         tuv012 = 0.0d0
         tuv111 = 0.0d0
         k0 = kgrd0
         do it3 = 1, bsorder
            k0 = k0 + 1
            k = k0 + 1 + (nfft3-isign(nfft3,k0))/2
            v0 = thetai3(1,it3,iatm)
            v1 = thetai3(2,it3,iatm)
            v2 = thetai3(3,it3,iatm)
            v3 = thetai3(4,it3,iatm)
            tu00_1 = 0.0d0
            tu01_1 = 0.0d0
            tu10_1 = 0.0d0
            tu20_1 = 0.0d0
            tu11_1 = 0.0d0
            tu02_1 = 0.0d0
            tu00_2 = 0.0d0
            tu01_2 = 0.0d0
            tu10_2 = 0.0d0
            tu20_2 = 0.0d0
            tu11_2 = 0.0d0
            tu02_2 = 0.0d0
            tu00 = 0.0d0
            tu10 = 0.0d0
            tu01 = 0.0d0
            tu20 = 0.0d0
            tu11 = 0.0d0
            tu02 = 0.0d0
            tu30 = 0.0d0
            tu21 = 0.0d0
            tu12 = 0.0d0
            tu03 = 0.0d0
            j0 = jgrd0
            do it2 = 1, bsorder
               j0 = j0 + 1
               j = j0 + 1 + (nfft2-isign(nfft2,j0))/2
               u0 = thetai2(1,it2,iatm)
               u1 = thetai2(2,it2,iatm)
               u2 = thetai2(3,it2,iatm)
               u3 = thetai2(4,it2,iatm)
               t0_1 = 0.0d0
               t1_1 = 0.0d0
               t2_1 = 0.0d0
               t0_2 = 0.0d0
               t1_2 = 0.0d0
               t2_2 = 0.0d0
               t3 = 0.0d0
               i0 = igrd0
               do it1 = 1, bsorder
                  i0 = i0 + 1
                  i = i0 + 1 + (nfft1-isign(nfft1,i0))/2
                  tq_1 = qgrid(1,i,j,k)
                  tq_2 = qgrid(2,i,j,k)
                  t0_1 = t0_1 + tq_1*thetai1(1,it1,iatm)
                  t1_1 = t1_1 + tq_1*thetai1(2,it1,iatm)
                  t2_1 = t2_1 + tq_1*thetai1(3,it1,iatm)
                  t0_2 = t0_2 + tq_2*thetai1(1,it1,iatm)
                  t1_2 = t1_2 + tq_2*thetai1(2,it1,iatm)
                  t2_2 = t2_2 + tq_2*thetai1(3,it1,iatm)
                  t3 = t3 + (tq_1+tq_2)*thetai1(4,it1,iatm)
               end do
               tu00_1 = tu00_1 + t0_1*u0
               tu10_1 = tu10_1 + t1_1*u0
               tu01_1 = tu01_1 + t0_1*u1
               tu20_1 = tu20_1 + t2_1*u0
               tu11_1 = tu11_1 + t1_1*u1
               tu02_1 = tu02_1 + t0_1*u2
               tu00_2 = tu00_2 + t0_2*u0
               tu10_2 = tu10_2 + t1_2*u0
               tu01_2 = tu01_2 + t0_2*u1
               tu20_2 = tu20_2 + t2_2*u0
               tu11_2 = tu11_2 + t1_2*u1
               tu02_2 = tu02_2 + t0_2*u2
               t0 = t0_1 + t0_2
               t1 = t1_1 + t1_2
               t2 = t2_1 + t2_2
               tu00 = tu00 + t0*u0
               tu10 = tu10 + t1*u0
               tu01 = tu01 + t0*u1
               tu20 = tu20 + t2*u0
               tu11 = tu11 + t1*u1
               tu02 = tu02 + t0*u2
               tu30 = tu30 + t3*u0
               tu21 = tu21 + t2*u1
               tu12 = tu12 + t1*u2
               tu03 = tu03 + t0*u3
            end do
            tuv100_1 = tuv100_1 + tu10_1*v0
            tuv010_1 = tuv010_1 + tu01_1*v0
            tuv001_1 = tuv001_1 + tu00_1*v1
            tuv200_1 = tuv200_1 + tu20_1*v0
            tuv020_1 = tuv020_1 + tu02_1*v0
            tuv002_1 = tuv002_1 + tu00_1*v2
            tuv110_1 = tuv110_1 + tu11_1*v0
            tuv101_1 = tuv101_1 + tu10_1*v1
            tuv011_1 = tuv011_1 + tu01_1*v1
            tuv100_2 = tuv100_2 + tu10_2*v0
            tuv010_2 = tuv010_2 + tu01_2*v0
            tuv001_2 = tuv001_2 + tu00_2*v1
            tuv200_2 = tuv200_2 + tu20_2*v0
            tuv020_2 = tuv020_2 + tu02_2*v0
            tuv002_2 = tuv002_2 + tu00_2*v2
            tuv110_2 = tuv110_2 + tu11_2*v0
            tuv101_2 = tuv101_2 + tu10_2*v1
            tuv011_2 = tuv011_2 + tu01_2*v1
            tuv000 = tuv000 + tu00*v0
            tuv100 = tuv100 + tu10*v0
            tuv010 = tuv010 + tu01*v0
            tuv001 = tuv001 + tu00*v1
            tuv200 = tuv200 + tu20*v0
            tuv020 = tuv020 + tu02*v0
            tuv002 = tuv002 + tu00*v2
            tuv110 = tuv110 + tu11*v0
            tuv101 = tuv101 + tu10*v1
            tuv011 = tuv011 + tu01*v1
            tuv300 = tuv300 + tu30*v0
            tuv030 = tuv030 + tu03*v0
            tuv003 = tuv003 + tu00*v3
            tuv210 = tuv210 + tu21*v0
            tuv201 = tuv201 + tu20*v1
            tuv120 = tuv120 + tu12*v0
            tuv021 = tuv021 + tu02*v1
            tuv102 = tuv102 + tu10*v2
            tuv012 = tuv012 + tu01*v2
            tuv111 = tuv111 + tu11*v1
         end do
         fdip_phi1(2,isite) = tuv100_1
         fdip_phi1(3,isite) = tuv010_1
         fdip_phi1(4,isite) = tuv001_1
         fdip_phi1(5,isite) = tuv200_1
         fdip_phi1(6,isite) = tuv020_1
         fdip_phi1(7,isite) = tuv002_1
         fdip_phi1(8,isite) = tuv110_1
         fdip_phi1(9,isite) = tuv101_1
         fdip_phi1(10,isite) = tuv011_1
         fdip_phi2(2,isite) = tuv100_2
         fdip_phi2(3,isite) = tuv010_2
         fdip_phi2(4,isite) = tuv001_2
         fdip_phi2(5,isite) = tuv200_2
         fdip_phi2(6,isite) = tuv020_2
         fdip_phi2(7,isite) = tuv002_2
         fdip_phi2(8,isite) = tuv110_2
         fdip_phi2(9,isite) = tuv101_2
         fdip_phi2(10,isite) = tuv011_2
         fdip_sum_phi(1,isite) = tuv000
         fdip_sum_phi(2,isite) = tuv100
         fdip_sum_phi(3,isite) = tuv010
         fdip_sum_phi(4,isite) = tuv001
         fdip_sum_phi(5,isite) = tuv200
         fdip_sum_phi(6,isite) = tuv020
         fdip_sum_phi(7,isite) = tuv002
         fdip_sum_phi(8,isite) = tuv110
         fdip_sum_phi(9,isite) = tuv101
         fdip_sum_phi(10,isite) = tuv011
         fdip_sum_phi(11,isite) = tuv300
         fdip_sum_phi(12,isite) = tuv030
         fdip_sum_phi(13,isite) = tuv003
         fdip_sum_phi(14,isite) = tuv210
         fdip_sum_phi(15,isite) = tuv201
         fdip_sum_phi(16,isite) = tuv120
         fdip_sum_phi(17,isite) = tuv021
         fdip_sum_phi(18,isite) = tuv102
         fdip_sum_phi(19,isite) = tuv012
         fdip_sum_phi(20,isite) = tuv111
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine cmp_to_fmp  --  transformation of multipoles  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "cmp_to_fmp" transforms the atomic multipoles from Cartesian
c     to fractional coordinates
c
c
      subroutine cmp_to_fmp (cmp,fmp)
      use sizes
      use mpole
      implicit none
      integer i,j,k
      real*8 ctf(10,10)
      real*8 cmp(10,*)
      real*8 fmp(10,*)
c
c
c     find the matrix to convert Cartesian to fractional
c
      call cart_to_frac (ctf)
c
c     apply the transformation to get the fractional multipoles
c
      do i = 1, npole
         fmp(1,i) = ctf(1,1) * cmp(1,i)
         do j = 2, 4
            fmp(j,i) = 0.0d0
            do k = 2, 4
               fmp(j,i) = fmp(j,i) + ctf(j,k)*cmp(k,i)
            end do
         end do
         do j = 5, 10
            fmp(j,i) = 0.0d0
            do k = 5, 10
               fmp(j,i) = fmp(j,i) + ctf(j,k)*cmp(k,i)
            end do
         end do
      end do
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine cart_to_frac  --  Cartesian to fractional  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "cart_to_frac" computes a transformation matrix to convert
c     a multipole object in Cartesian coordinates to fractional
c
c     note the multipole components are stored in the condensed
c     order (m,dx,dy,dz,qxx,qyy,qzz,qxy,qxz,qyz)
c
c
      subroutine cart_to_frac (ctf)
      use sizes
      use boxes
      use pme
      implicit none
      integer i,j,k,m
      integer i1,i2
      integer qi1(6)
      integer qi2(6)
      real*8 a(3,3)
      real*8 ctf(10,10)
      data qi1  / 1, 2, 3, 1, 1, 2 /
      data qi2  / 1, 2, 3, 2, 3, 3 /
c
c
c     set the reciprocal vector transformation matrix
c
      do i = 1, 3
         a(1,i) = dble(nfft1) * recip(i,1)
         a(2,i) = dble(nfft2) * recip(i,2)
         a(3,i) = dble(nfft3) * recip(i,3)
      end do
c
c     get the Cartesian to fractional conversion matrix
c
      do i = 1, 10
         do j = 1, 10
            ctf(j,i) = 0.0d0
         end do
      end do
      ctf(1,1) = 1.0d0
      do i = 2, 4
         do j = 2, 4
            ctf(i,j) = a(i-1,j-1)
         end do
      end do
      do i1 = 1, 3
         k = qi1(i1)
         do i2 = 1, 6
            i = qi1(i2)
            j = qi2(i2)
            ctf(i1+4,i2+4) = a(k,i) * a(k,j)
         end do
      end do
      do i1 = 4, 6
         k = qi1(i1)
         m = qi2(i1)
         do i2 = 1, 6
            i = qi1(i2)
            j = qi2(i2)
            ctf(i1+4,i2+4) = a(k,i)*a(m,j) + a(k,j)*a(m,i)
         end do
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine fphi_to_cphi  --  transformation of potential  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "fphi_to_cphi" transforms the reciprocal space potential from
c     fractional to Cartesian coordinates
c
c
      subroutine fphi_to_cphi (fphi,cphi)
      use sizes
      use mpole
      implicit none
      integer i,j,k
      real*8 ftc(10,10)
      real*8 cphi(10,*)
      real*8 fphi(20,*)
c
c
c     find the matrix to convert fractional to Cartesian
c
      call frac_to_cart (ftc)
c
c     apply the transformation to get the Cartesian potential
c
      do i = 1, npole
         cphi(1,i) = ftc(1,1) * fphi(1,i)
         do j = 2, 4
            cphi(j,i) = 0.0d0
            do k = 2, 4
               cphi(j,i) = cphi(j,i) + ftc(j,k)*fphi(k,i)
            end do
         end do
         do j = 5, 10
            cphi(j,i) = 0.0d0
            do k = 5, 10
               cphi(j,i) = cphi(j,i) + ftc(j,k)*fphi(k,i)
            end do
         end do
      end do
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine frac_to_cart  --  fractional to Cartesian  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "frac_to_cart" computes a transformation matrix to convert
c     a multipole object in fraction coordinates to Cartesian
c
c     note the multipole components are stored in the condensed
c     order (m,dx,dy,dz,qxx,qyy,qzz,qxy,qxz,qyz)
c
c
      subroutine frac_to_cart (ftc)
      use sizes
      use boxes
      use pme
      implicit none
      integer i,j,k,m
      integer i1,i2
      integer qi1(6)
      integer qi2(6)
      real*8 a(3,3)
      real*8 ftc(10,10)
      data qi1  / 1, 2, 3, 1, 1, 2 /
      data qi2  / 1, 2, 3, 2, 3, 3 /
c
c
c     set the reciprocal vector transformation matrix
c
      do i = 1, 3
         a(i,1) = dble(nfft1) * recip(i,1)
         a(i,2) = dble(nfft2) * recip(i,2)
         a(i,3) = dble(nfft3) * recip(i,3)
      end do
c
c     get the fractional to Cartesian conversion matrix
c
      do i = 1, 10
         do j = 1, 10
            ftc(j,i) = 0.0d0
         end do
      end do
      ftc(1,1) = 1.0d0
      do i = 2, 4
         do j = 2, 4
            ftc(i,j) = a(i-1,j-1)
         end do
      end do
      do i1 = 1, 3
         k = qi1(i1)
         do i2 = 1, 3
            i = qi1(i2)
            ftc(i1+4,i2+4) = a(k,i) * a(k,i)
         end do
         do i2 = 4, 6
            i = qi1(i2)
            j = qi2(i2)
            ftc(i1+4,i2+4) = 2.0d0 * a(k,i) * a(k,j)
         end do
      end do
      do i1 = 4, 6
         k = qi1(i1)
         m = qi2(i1)
         do i2 = 1, 3
            i = qi1(i2)
            ftc(i1+4,i2+4) = a(k,i) * a(m,i)
         end do
         do i2 = 4, 6
            i = qi1(i2)
            j = qi2(i2)
            ftc(i1+4,i2+4) = a(k,i)*a(m,j) + a(m,i)*a(k,j)
         end do
      end do
      return
      end
