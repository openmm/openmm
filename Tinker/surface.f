c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine surface  --  find the accessible surface area  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "surface" performs an analytical computation of the weighted
c     solvent accessible surface area of each atom and the first
c     derivatives of the area with respect to Cartesian coordinates
c
c     literature references:
c
c     T. J. Richmond, "Solvent Accessible Surface Area and
c     Excluded Volume in Proteins", Journal of Molecular Biology,
c     178, 63-89 (1984)
c
c     L. Wesson and D. Eisenberg, "Atomic Solvation Parameters
c     Applied to Molecular Dynamics of Proteins in Solution",
c     Protein Science, 1, 227-235 (1992)
c
c     variables and parameters:
c
c     total    total surface area of the whole structure
c     area     accessible surface area of each atom
c     radius   radii of the individual atoms
c     weight   weight assigned to each atom's area; if set to
c                1.0, return is actual area in square Angstroms
c     probe    radius of the probe sphere
c     delta    tolerance used in the tests for sphere overlaps
c                and for colinearity
c     rmove    connectivity errors can usually be avoided if the
c                offending atom is shifted by this small amount
c
c
      subroutine surface (total,area,radius,weight,probe)
      use sizes
      use atoms
      use inform
      use iounit
      use math
      use usage
      implicit none
      integer maxarc
      parameter (maxarc=1000)
      integer i,j,k,l,m
      integer ii,ib,jb
      integer io,ir
      integer mi,ni,narc
      integer key(maxarc)
      integer intag(maxarc)
      integer intag1(maxarc)
      integer lt(maxarc)
      integer kent(maxarc)
      integer kout(maxarc)
      real*8 total,wght
      real*8 delta,delta2
      real*8 eps,rmove,dsql
      real*8 probe,arcsum,cosine
      real*8 axx,axy,axz
      real*8 ayx,ayy,azx
      real*8 azy,azz
      real*8 uxl,uyl,uzl
      real*8 tx,ty,tz
      real*8 txb,tyb,td
      real*8 tr2,tr,txr,tyr
      real*8 tk1,tk2
      real*8 thec,the
      real*8 t,tb,txk,tyk,tzk
      real*8 t1,ti,tf,tt
      real*8 txl,tyl,tzl
      real*8 arclen,exang
      real*8 xr,yr,zr
      real*8 rr,rrx2,rrsq
      real*8 rplus,rminus
      real*8 ccsq,cc,xysq
      real*8 bk,gi,bsqk
      real*8 pix2,pix4,pid2
      real*8 therk,dk,gk
      real*8 risqk,rik
      real*8 area(*)
      real*8 radius(*)
      real*8 weight(*)
      real*8 ri(maxarc),risq(maxarc)
      real*8 bsq(maxarc),bsq1(maxarc)
      real*8 dsq(maxarc),dsq1(maxarc)
      real*8 arci(maxarc),arcf(maxarc)
      real*8 ex(maxarc),gr(maxarc)
      real*8 b(maxarc),b1(maxarc)
      real*8 bg(maxarc),ther(maxarc)
      real*8 xc(maxarc),xc1(maxarc)
      real*8 yc(maxarc),yc1(maxarc)
      real*8 zc(maxarc),zc1(maxarc)
      real*8 ux(maxarc),uy(maxarc)
      real*8 uz(maxarc)
      real*8, allocatable :: r(:)
      logical moved,top,komit
      logical omit(maxarc)
      logical, allocatable :: skip(:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (r(n))
      allocate (skip(n))
c
c     zero the area and derivatives, and set the sphere radii
c
      total = 0.0d0
      do i = 1, n
         area(i) = 0.0d0
         r(i) = radius(i)
         if (r(i) .ne. 0.0d0)  r(i) = r(i) + probe
      end do
c
c     set pi multiples, overlap criterion and tolerances
c
      pix2 = 2.0d0 * pi
      pix4 = 4.0d0 * pi
      pid2 = pi / 2.0d0
      delta = 1.0d-8
      delta2 = delta**2
      eps = 1.0d-8
      rmove = 1.0d-8
c
c     set the "skip" array to exclude all inactive atoms
c     that do not overlap any of the current active atoms
c
      do i = 1, n
         skip(i) = .true.
      end do
      do i = 1, n
         if (use(i)) then
            xr = x(i)
            yr = y(i)
            zr = z(i)
            rr = r(i)
            do k = 1, n
               rplus = (rr + r(k))**2
               ccsq = (x(k)-xr)**2 + (y(k)-yr)**2 + (z(k)-zr)**2
               if (ccsq .le. rplus)  skip(k) = .false.
            end do
         end if
      end do
c
c     compute the area and derivatives of current "ir" sphere
c
      do ir = 1, n
         if (skip(ir))  goto 180
         xr = x(ir)
         yr = y(ir)
         zr = z(ir)
         rr = r(ir)
         rrx2 = 2.0d0 * rr
         rrsq = rr**2
         wght = weight(ir)
         moved = .false.
c
c     initialize some counters and sums for the "ir" sphere
c
   10    continue
         io = 0
         jb = 0
         ib = 0
         arclen = 0.0d0
         exang = 0.0d0
c
c     test each sphere to see if it overlaps the "ir" sphere
c
         do i = 1, n
            if (i .eq. ir)  goto 30
            rplus = rr + r(i)
            tx = x(i) - xr
            if (abs(tx) .ge. rplus)  goto 30
            ty = y(i) - yr
            if (abs(ty) .ge. rplus)  goto 30
            tz = z(i) - zr
            if (abs(tz) .ge. rplus)  goto 30
c
c     check for overlap of spheres by testing center to
c     center distance against sum and difference of radii
c
            xysq = tx**2 + ty**2
            if (xysq .lt. delta2) then
               tx = delta
               ty = 0.0d0
               xysq = delta2
            end if
            ccsq = xysq + tz**2
            cc = sqrt(ccsq)
            if (rplus-cc .le. delta)  goto 30
            rminus = rr - r(i)
c
c     check for a completely buried "ir" sphere
c
            if (cc-abs(rminus) .le. delta) then
               if (rminus .le. 0.0d0)  goto 180
               goto 30
            end if
c
c     calculate overlap parameters between "i" and "ir" sphere
c
            io = io + 1
            xc1(io) = tx
            yc1(io) = ty
            zc1(io) = tz
            dsq1(io) = xysq
            bsq1(io) = ccsq
            b1(io) = cc
            gr(io) = (ccsq+rplus*rminus) / (rrx2*b1(io))
            intag1(io) = i
            if (io .gt. maxarc) then
               write (iout,20)
   20          format (/,' SURFACE  --  Increase the Value of MAXARC')
               call fatal
            end if
   30       continue
         end do
c
c     case where no other spheres overlap the current sphere
c
         if (io .eq. 0) then
            area(ir) = pix4
            goto 160
         end if
c
c     case where only one sphere overlaps the current sphere
c
         if (io .eq. 1) then
            k = 1
            txk = xc1(1)
            tyk = yc1(1)
            tzk = zc1(1)
            bsqk = bsq1(1)
            bk = b1(1)
            intag(1) = intag1(1)
            arcsum = pix2
            ib = ib + 1
            arclen = arclen + gr(k)*arcsum
            goto 150
         end if
c
c     general case where more than one sphere intersects the
c     current sphere; sort intersecting spheres by their degree
c     of overlap with the current main sphere
c
         call sort2 (io,gr,key)
         do i = 1, io
            k = key(i)
            intag(i) = intag1(k)
            xc(i) = xc1(k)
            yc(i) = yc1(k)
            zc(i) = zc1(k)
            dsq(i) = dsq1(k)
            b(i) = b1(k)
            bsq(i) = bsq1(k)
            omit(i) = .false.
         end do
c
c     radius of the each circle on the surface of the "ir" sphere
c
         do i = 1, io
            gi = gr(i) * rr
            bg(i) = b(i) * gi
            risq(i) = rrsq - gi**2
            ri(i) = sqrt(risq(i))
            ther(i) = pid2 - asin(min(1.0d0,max(-1.0d0,gr(i))))
         end do
c
c     find boundary of inaccessible area on "ir" sphere
c
         do k = 1, io-1
            if (.not. omit(k)) then
               txk = xc(k)
               tyk = yc(k)
               tzk = zc(k)
               bk = b(k)
               therk = ther(k)
               do j = k+1, io
                  if (omit(j))  goto 60
c
c     check to see if J circle is intersecting K circle;
c     get distance between circle centers and sum of radii
c
                  cc = (txk*xc(j)+tyk*yc(j)+tzk*zc(j))/(bk*b(j))
                  cc = acos(min(1.0d0,max(-1.0d0,cc)))
                  td = therk + ther(j)
c
c     check to see if circles enclose separate regions
c
                  if (cc .ge. td)  goto 60
c
c     check for circle J completely inside circle K
c
                  if (cc+ther(j) .lt. therk)  goto 40
c
c     check for circles essentially parallel
c
                  if (cc .gt. delta)  goto 50
   40             continue
                  omit(j) = .true.
                  goto 60
c
c     check for "ir" sphere completely buried
c
   50             continue
                  if (pix2-cc .le. td)  goto 180
   60             continue
               end do
            end if
         end do
c
c     find T value of circle intersections
c
         do k = 1, io
            if (omit(k))  goto 110
            komit = omit(k)
            omit(k) = .true.
            narc = 0
            top = .false.
            txk = xc(k)
            tyk = yc(k)
            tzk = zc(k)
            dk = sqrt(dsq(k))
            bsqk = bsq(k)
            bk = b(k)
            gk = gr(k) * rr
            risqk = risq(k)
            rik = ri(k)
            therk = ther(k)
c
c     rotation matrix elements
c
            t1 = tzk / (bk*dk)
            axx = txk * t1
            axy = tyk * t1
            axz = dk / bk
            ayx = tyk / dk
            ayy = txk / dk
            azx = txk / bk
            azy = tyk / bk
            azz = tzk / bk
            do l = 1, io
               if (.not. omit(l)) then
                  txl = xc(l)
                  tyl = yc(l)
                  tzl = zc(l)
c
c     rotate spheres so K vector colinear with z-axis
c
                  uxl = txl*axx + tyl*axy - tzl*axz
                  uyl = tyl*ayy - txl*ayx
                  uzl = txl*azx + tyl*azy + tzl*azz
                  cosine = min(1.0d0,max(-1.0d0,uzl/b(l)))
                  if (acos(cosine) .lt. therk+ther(l)) then
                     dsql = uxl**2 + uyl**2
                     tb = uzl*gk - bg(l)
                     txb = uxl * tb
                     tyb = uyl * tb
                     td = rik * dsql
                     tr2 = risqk*dsql - tb**2
                     tr2 = max(eps,tr2)
                     tr = sqrt(tr2)
                     txr = uxl * tr
                     tyr = uyl * tr
c
c     get T values of intersection for K circle
c
                     tb = (txb+tyr) / td
                     tb = min(1.0d0,max(-1.0d0,tb))
                     tk1 = acos(tb)
                     if (tyb-txr .lt. 0.0d0)  tk1 = pix2 - tk1
                     tb = (txb-tyr) / td
                     tb = min(1.0d0,max(-1.0d0,tb))
                     tk2 = acos(tb)
                     if (tyb+txr .lt. 0.0d0)  tk2 = pix2 - tk2
                     thec = (rrsq*uzl-gk*bg(l)) / (rik*ri(l)*b(l))
                     if (abs(thec) .lt. 1.0d0) then
                        the = -acos(thec)
                     else if (thec .ge. 1.0d0) then
                        the = 0.0d0
                     else if (thec .le. -1.0d0) then
                        the = -pi
                     end if
c
c     see if "tk1" is entry or exit point; check t=0 point;
c     "ti" is exit point, "tf" is entry point
c
                     cosine = min(1.0d0,max(-1.0d0,
     &                               (uzl*gk-uxl*rik)/(b(l)*rr)))
                     if ((acos(cosine)-ther(l))*(tk2-tk1)
     &                          .le. 0.0d0) then
                        ti = tk2
                        tf = tk1
                     else
                        ti = tk1
                        tf = tk2
                     end if
                     narc = narc + 1
                     if (narc .ge. maxarc) then
                        write (iout,70)
   70                   format (/,' SURFACE  --  Increase the Value',
     &                             ' of MAXARC')
                        call fatal
                     end if
                     if (tf .le. ti) then
                        arcf(narc) = tf
                        arci(narc) = 0.0d0
                        tf = pix2
                        lt(narc) = l
                        ex(narc) = the
                        top = .true.
                        narc = narc + 1
                     end if
                     arcf(narc) = tf
                     arci(narc) = ti
                     lt(narc) = l
                     ex(narc) = the
                     ux(l) = uxl
                     uy(l) = uyl
                     uz(l) = uzl
                  end if
               end if
            end do
            omit(k) = komit
c
c     special case; K circle without intersections
c
            if (narc .le. 0)  goto 90
c
c     general case; sum up arclength and set connectivity code
c
            call sort2 (narc,arci,key)
            arcsum = arci(1)
            mi = key(1)
            t = arcf(mi)
            ni = mi
            if (narc .gt. 1) then
               do j = 2, narc
                  m = key(j)
                  if (t .lt. arci(j)) then
                     arcsum = arcsum + arci(j) - t
                     exang = exang + ex(ni)
                     jb = jb + 1
                     if (jb .ge. maxarc) then
                        write (iout,80)
   80                   format (/,' SURFACE  --  Increase the Value',
     &                             ' of MAXARC')
                        call fatal
                     end if
                     l = lt(ni)
                     kent(jb) = maxarc*l + k
                     l = lt(m)
                     kout(jb) = maxarc*k + l
                  end if
                  tt = arcf(m)
                  if (tt .ge. t) then
                     t = tt
                     ni = m
                  end if
               end do
            end if
            arcsum = arcsum + pix2 - t
            if (.not. top) then
               exang = exang + ex(ni)
               jb = jb + 1
               l = lt(ni)
               kent(jb) = maxarc*l + k
               l = lt(mi)
               kout(jb) = maxarc*k + l
            end if
            goto 100
   90       continue
            arcsum = pix2
            ib = ib + 1
  100       continue
            arclen = arclen + gr(k)*arcsum
  110       continue
         end do
         if (arclen .eq. 0.0d0)  goto 180
         if (jb .eq. 0)  goto 150
c
c     find number of independent boundaries and check connectivity
c
         j = 0
         do k = 1, jb
            if (kout(k) .ne. 0) then
               i = k
  120          continue
               m = kout(i)
               kout(i) = 0
               j = j + 1
               do ii = 1, jb
                  if (m .eq. kent(ii)) then
                     if (ii .eq. k) then
                        ib = ib + 1
                        if (j .eq. jb)  goto 150
                        goto 130
                     end if
                     i = ii
                     goto 120
                  end if
               end do
  130          continue
            end if
         end do
         ib = ib + 1
c
c     attempt to fix connectivity error by moving atom slightly
c
         if (moved) then
            write (iout,140)  ir
  140       format (/,' SURFACE  --  Connectivity Error at Atom',i6)
         else
            moved = .true.
            xr = xr + rmove
            yr = yr + rmove
            zr = zr + rmove
            goto 10
         end if
c
c     form the accessible area for the current atom
c
  150    continue
         area(ir) = ib*pix2 + exang + arclen
         area(ir) = mod(area(ir),pix4)
  160    continue
         area(ir) = area(ir) * rrsq
c
c     attempt to fix negative area by moving atom slightly
c
         if (area(ir) .lt. 0.0d0) then
            if (moved) then
               write (iout,170)  ir
  170          format (/,' SURFACE  --  Negative Area at Atom',i6)
            else
               moved = .true.
               xr = xr + rmove
               yr = yr + rmove
               zr = zr + rmove
               goto 10
            end if
         end if
c
c     weight the accessible area by the scale factor
c
         area(ir) = area(ir) * wght
         total = total + area(ir)
  180    continue
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (r)
      deallocate (skip)
c
c     print out the surface area values for each atom
c
c     if (debug) then
c        write (iout,190)
c 190    format (/,' Weighted Atomic Surface Areas Values :',
c    &           //,4x,'Atom',7x,'Area Term',6x,'Weight',/)
c        do i = 1, n
c           if (.not. skip(i)) then
c              write (iout,200)  i,area(i),weight(i)
c 200          format (i8,4x,2f12.4)
c           end if
c        end do
c        write (iout,210)  total
c 210    format (/,' Total Weighted Surface Area :',5x,f16.4)
c     end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine surface1  --  accessible surface area & derivs  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "surface1" performs an analytical computation of the weighted
c     solvent accessible surface area of each atom and the first
c     derivatives of the area with respect to Cartesian coordinates
c
c     literature references:
c
c     T. J. Richmond, "Solvent Accessible Surface Area and
c     Excluded Volume in Proteins", Journal of Molecular Biology,
c     178, 63-89 (1984)
c
c     L. Wesson and D. Eisenberg, "Atomic Solvation Parameters
c     Applied to Molecular Dynamics of Proteins in Solution",
c     Protein Science, 1, 227-235 (1992)
c
c     variables and parameters:
c
c     total    total surface area of the whole structure
c     area     accessible surface area of each atom
c     darea    x,y,z components of the gradient of the area of
c                the molecule with respect to atomic coordinates
c     radius   radii of the individual atoms
c     weight   weight assigned to each atom's area; if set to
c                1.0, return is actual area in square Angstroms
c     probe    radius of the probe sphere
c     delta    tolerance used in the tests for sphere overlaps
c                and for colinearity
c     rmove    connectivity errors can usually be avoided if the
c                offending atom is shifted by this small amount
c
c
      subroutine surface1 (total,area,darea,radius,weight,probe)
      use sizes
      use atoms
      use inform
      use iounit
      use math
      use usage
      implicit none
      integer maxarc
      parameter (maxarc=1000)
      integer i,j,k,l,m
      integer ii,ib,jb
      integer in,io,ir
      integer mi,ni,narc
      integer key(maxarc)
      integer intag(maxarc)
      integer intag1(maxarc)
      integer lt(maxarc)
      integer kent(maxarc)
      integer kout(maxarc)
      integer ider(maxarc)
      integer sign_yder(maxarc)
      real*8 total,wght
      real*8 delta,delta2
      real*8 eps,rmove
      real*8 probe,arcsum,cosine
      real*8 dsql,wxl,wxlsq
      real*8 p,s,v,rcn
      real*8 axx,axy,axz
      real*8 ayx,ayy,azx
      real*8 azy,azz
      real*8 uxl,uyl,uzl
      real*8 tx,ty,tz
      real*8 txb,tyb,t2,td
      real*8 tr2,tr,txr,tyr
      real*8 tk1,tk2
      real*8 thec,the
      real*8 t,tb,txk,tyk,tzk
      real*8 t1,ti,tf,tt
      real*8 txl,tyl,tzl
      real*8 arclen,exang
      real*8 xr,yr,zr
      real*8 rr,rrx2,rrsq
      real*8 rplus,rminus
      real*8 ccsq,cc,xysq
      real*8 bgl,bsqk,bsql
      real*8 bk,gi,gl
      real*8 pix2,pix4,pid2
      real*8 dax,day,daz
      real*8 deal,decl
      real*8 dtkal,dtkcl
      real*8 dtlal,dtlcl
      real*8 therk,dk,gk
      real*8 risqk,rik,risql
      real*8 faca,facb,facc
      real*8 gaca,gacb
      real*8 area(*)
      real*8 radius(*)
      real*8 weight(*)
      real*8 darea(3,*)
      real*8 ri(maxarc),risq(maxarc)
      real*8 bsq(maxarc),bsq1(maxarc)
      real*8 dsq(maxarc),dsq1(maxarc)
      real*8 arci(maxarc),arcf(maxarc)
      real*8 ex(maxarc),gr(maxarc)
      real*8 b(maxarc),b1(maxarc)
      real*8 bg(maxarc),ther(maxarc)
      real*8 xc(maxarc),xc1(maxarc)
      real*8 yc(maxarc),yc1(maxarc)
      real*8 zc(maxarc),zc1(maxarc)
      real*8 ux(maxarc),uy(maxarc)
      real*8 uz(maxarc)
      real*8, allocatable :: r(:)
      logical moved,top,komit
      logical omit(maxarc)
      logical, allocatable :: skip(:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (r(n))
      allocate (skip(n))
c
c     zero the area and derivatives, and set the sphere radii
c
      total = 0.0d0
      do i = 1, n
         area(i) = 0.0d0
         darea(1,i) = 0.0d0
         darea(2,i) = 0.0d0
         darea(3,i) = 0.0d0
         r(i) = radius(i)
         if (r(i) .ne. 0.0d0)  r(i) = r(i) + probe
      end do
c
c     set pi multiples, overlap criterion and tolerances
c
      pix2 = 2.0d0 * pi
      pix4 = 4.0d0 * pi
      pid2 = pi / 2.0d0
      delta = 1.0d-8
      delta2 = delta**2
      eps = 1.0d-8
      rmove = 1.0d-8
      do i = 1, maxarc
         ider(i) = 0
         sign_yder(i) = 0
      end do
c
c     set the "skip" array to exclude all inactive atoms
c     that do not overlap any of the current active atoms
c
      do i = 1, n
         skip(i) = .true.
      end do
      do i = 1, n
         if (use(i)) then
            xr = x(i)
            yr = y(i)
            zr = z(i)
            rr = r(i)
            do k = 1, n
               rplus = (rr + r(k))**2
               ccsq = (x(k)-xr)**2 + (y(k)-yr)**2 + (z(k)-zr)**2
               if (ccsq .le. rplus)  skip(k) = .false.
            end do
         end if
      end do
c
c     compute the area and derivatives of current "ir" sphere
c
      do ir = 1, n
         if (skip(ir))  goto 180
         xr = x(ir)
         yr = y(ir)
         zr = z(ir)
         rr = r(ir)
         rrx2 = 2.0d0 * rr
         rrsq = rr**2
         wght = weight(ir)
         moved = .false.
c
c     initialize some counters and sums for the "ir" sphere
c
   10    continue
         io = 0
         jb = 0
         ib = 0
         arclen = 0.0d0
         exang = 0.0d0
c
c     test each sphere to see if it overlaps the "ir" sphere
c
         do i = 1, n
            if (i .eq. ir)  goto 30
            rplus = rr + r(i)
            tx = x(i) - xr
            if (abs(tx) .ge. rplus)  goto 30
            ty = y(i) - yr
            if (abs(ty) .ge. rplus)  goto 30
            tz = z(i) - zr
            if (abs(tz) .ge. rplus)  goto 30
c
c     check for overlap of spheres by testing center to
c     center distance against sum and difference of radii
c
            xysq = tx**2 + ty**2
            if (xysq .lt. delta2) then
               tx = delta
               ty = 0.0d0
               xysq = delta2
            end if
            ccsq = xysq + tz**2
            cc = sqrt(ccsq)
            if (rplus-cc .le. delta)  goto 30
            rminus = rr - r(i)
c
c     check for a completely buried "ir" sphere
c
            if (cc-abs(rminus) .le. delta) then
               if (rminus .le. 0.0d0)  goto 180
               goto 30
            end if
c
c     calculate overlap parameters between "i" and "ir" sphere
c
            io = io + 1
            xc1(io) = tx
            yc1(io) = ty
            zc1(io) = tz
            dsq1(io) = xysq
            bsq1(io) = ccsq
            b1(io) = cc
            gr(io) = (ccsq+rplus*rminus) / (rrx2*b1(io))
            intag1(io) = i
            if (io .gt. maxarc) then
               write (iout,20)
   20          format (/,' SURFACE1  --  Increase the Value of MAXARC')
               call fatal
            end if
   30       continue
         end do
c
c     case where no other spheres overlap the current sphere
c
         if (io .eq. 0) then
            area(ir) = pix4
            goto 160
         end if
c
c     case where only one sphere overlaps the current sphere
c
         if (io .eq. 1) then
            k = 1
            txk = xc1(1)
            tyk = yc1(1)
            tzk = zc1(1)
            bsqk = bsq1(1)
            bk = b1(1)
            intag(1) = intag1(1)
            arcsum = pix2
            ib = ib + 1
            arclen = arclen + gr(k)*arcsum
            if (.not. moved) then
               in = intag(k)
               t1 = arcsum*rrsq*(bsqk-rrsq+r(in)**2) / (rrx2*bsqk*bk)
               darea(1,ir) = darea(1,ir) - txk*t1*wght
               darea(2,ir) = darea(2,ir) - tyk*t1*wght
               darea(3,ir) = darea(3,ir) - tzk*t1*wght
               darea(1,in) = darea(1,in) + txk*t1*wght
               darea(2,in) = darea(2,in) + tyk*t1*wght
               darea(3,in) = darea(3,in) + tzk*t1*wght
            end if
            goto 150
         end if
c
c     general case where more than one sphere intersects the
c     current sphere; sort intersecting spheres by their degree
c     of overlap with the current main sphere
c
         call sort2 (io,gr,key)
         do i = 1, io
            k = key(i)
            intag(i) = intag1(k)
            xc(i) = xc1(k)
            yc(i) = yc1(k)
            zc(i) = zc1(k)
            dsq(i) = dsq1(k)
            b(i) = b1(k)
            bsq(i) = bsq1(k)
            omit(i) = .false.
         end do
c
c     radius of the each circle on the surface of the "ir" sphere
c
         do i = 1, io
            gi = gr(i) * rr
            bg(i) = b(i) * gi
            risq(i) = rrsq - gi**2
            ri(i) = sqrt(risq(i))
            ther(i) = pid2 - asin(min(1.0d0,max(-1.0d0,gr(i))))
         end do
c
c     find boundary of inaccessible area on "ir" sphere
c
         do k = 1, io-1
            if (.not. omit(k)) then
               txk = xc(k)
               tyk = yc(k)
               tzk = zc(k)
               bk = b(k)
               therk = ther(k)
               do j = k+1, io
                  if (omit(j))  goto 60
c
c     check to see if J circle is intersecting K circle;
c     get distance between circle centers and sum of radii
c
                  cc = (txk*xc(j)+tyk*yc(j)+tzk*zc(j))/(bk*b(j))
                  cc = acos(min(1.0d0,max(-1.0d0,cc)))
                  td = therk + ther(j)
c
c     check to see if circles enclose separate regions
c
                  if (cc .ge. td)  goto 60
c
c     check for circle J completely inside circle K
c
                  if (cc+ther(j) .lt. therk)  goto 40
c
c     check for circles essentially parallel
c
                  if (cc .gt. delta)  goto 50
   40             continue
                  omit(j) = .true.
                  goto 60
c
c     check for "ir" sphere completely buried
c
   50             continue
                  if (pix2-cc .le. td)  goto 180
   60             continue
               end do
            end if
         end do
c
c     find T value of circle intersections
c
         do k = 1, io
            if (omit(k))  goto 110
            komit = omit(k)
            omit(k) = .true.
            narc = 0
            top = .false.
            txk = xc(k)
            tyk = yc(k)
            tzk = zc(k)
            dk = sqrt(dsq(k))
            bsqk = bsq(k)
            bk = b(k)
            gk = gr(k) * rr
            risqk = risq(k)
            rik = ri(k)
            therk = ther(k)
c
c     rotation matrix elements
c
            t1 = tzk / (bk*dk)
            axx = txk * t1
            axy = tyk * t1
            axz = dk / bk
            ayx = tyk / dk
            ayy = txk / dk
            azx = txk / bk
            azy = tyk / bk
            azz = tzk / bk
            do l = 1, io
               if (.not. omit(l)) then
                  txl = xc(l)
                  tyl = yc(l)
                  tzl = zc(l)
c
c     rotate spheres so K vector colinear with z-axis
c
                  uxl = txl*axx + tyl*axy - tzl*axz
                  uyl = tyl*ayy - txl*ayx
                  uzl = txl*azx + tyl*azy + tzl*azz
                  cosine = min(1.0d0,max(-1.0d0,uzl/b(l)))
                  if (acos(cosine) .lt. therk+ther(l)) then
                     dsql = uxl**2 + uyl**2
                     tb = uzl*gk - bg(l)
                     txb = uxl * tb
                     tyb = uyl * tb
                     td = rik * dsql
                     tr2 = risqk*dsql - tb**2
                     tr2 = max(eps,tr2)
                     tr = sqrt(tr2)
                     txr = uxl * tr
                     tyr = uyl * tr
c
c     get T values of intersection for K circle
c
                     tb = (txb+tyr) / td
                     tb = min(1.0d0,max(-1.0d0,tb))
                     tk1 = acos(tb)
                     if (tyb-txr .lt. 0.0d0)  tk1 = pix2 - tk1
                     tb = (txb-tyr) / td
                     tb = min(1.0d0,max(-1.0d0,tb))
                     tk2 = acos(tb)
                     if (tyb+txr .lt. 0.0d0)  tk2 = pix2 - tk2
                     thec = (rrsq*uzl-gk*bg(l)) / (rik*ri(l)*b(l))
                     if (abs(thec) .lt. 1.0d0) then
                        the = -acos(thec)
                     else if (thec .ge. 1.0d0) then
                        the = 0.0d0
                     else if (thec .le. -1.0d0) then
                        the = -pi
                     end if
c
c     see if "tk1" is entry or exit point; check t=0 point;
c     "ti" is exit point, "tf" is entry point
c
                     cosine = min(1.0d0,max(-1.0d0,
     &                               (uzl*gk-uxl*rik)/(b(l)*rr)))
                     if ((acos(cosine)-ther(l))*(tk2-tk1)
     &                          .le. 0.0d0) then
                        ti = tk2
                        tf = tk1
                     else
                        ti = tk1
                        tf = tk2
                     end if
                     narc = narc + 1
                     if (narc .ge. maxarc) then
                        write (iout,70)
   70                   format (/,' SURFACE1  --  Increase the Value',
     &                             ' of MAXARC')
                        call fatal
                     end if
                     if (tf .le. ti) then
                        arcf(narc) = tf
                        arci(narc) = 0.0d0
                        tf = pix2
                        lt(narc) = l
                        ex(narc) = the
                        top = .true.
                        narc = narc + 1
                     end if
                     arcf(narc) = tf
                     arci(narc) = ti
                     lt(narc) = l
                     ex(narc) = the
                     ux(l) = uxl
                     uy(l) = uyl
                     uz(l) = uzl
                  end if
               end if
            end do
            omit(k) = komit
c
c     special case; K circle without intersections
c
            if (narc .le. 0)  goto 90
c
c     general case; sum up arclength and set connectivity code
c
            call sort2 (narc,arci,key)
            arcsum = arci(1)
            mi = key(1)
            t = arcf(mi)
            ni = mi
            if (narc .gt. 1) then
               do j = 2, narc
                  m = key(j)
                  if (t .lt. arci(j)) then
                     arcsum = arcsum + arci(j) - t
                     exang = exang + ex(ni)
                     jb = jb + 1
                     if (jb .ge. maxarc) then
                        write (iout,80)
   80                   format (/,' SURFACE1  --  Increase the Value',
     &                             ' of MAXARC')
                        call fatal
                     end if
                     l = lt(ni)
                     ider(l) = ider(l) + 1
                     sign_yder(l) = sign_yder(l) + 1
                     kent(jb) = maxarc*l + k
                     l = lt(m)
                     ider(l) = ider(l) + 1
                     sign_yder(l) = sign_yder(l) - 1
                     kout(jb) = maxarc*k + l
                  end if
                  tt = arcf(m)
                  if (tt .ge. t) then
                     t = tt
                     ni = m
                  end if
               end do
            end if
            arcsum = arcsum + pix2 - t
            if (.not. top) then
               exang = exang + ex(ni)
               jb = jb + 1
               l = lt(ni)
               ider(l) = ider(l) + 1
               sign_yder(l) = sign_yder(l) + 1
               kent(jb) = maxarc*l + k
               l = lt(mi)
               ider(l) = ider(l) + 1
               sign_yder(l) = sign_yder(l) - 1
               kout(jb) = maxarc*k + l
            end if
c
c     calculate the surface area derivatives
c
            do l = 1, io
               if (ider(l) .ne. 0) then
                  rcn = ider(l) * rrsq
                  ider(l) = 0
                  uzl = uz(l)
                  gl = gr(l) * rr
                  bgl = bg(l)
                  bsql = bsq(l)
                  risql = risq(l)
                  wxlsq = bsql - uzl**2
                  wxl = sqrt(wxlsq)
                  p = bgl - gk*uzl
                  v = risqk*wxlsq - p**2
                  v = max(eps,v)
                  v = sqrt(v)
                  t1 = rr * (gk*(bgl-bsql)+uzl*(bgl-rrsq))
     &                             / (v*risql*bsql)
                  deal = -wxl*t1
                  decl = -uzl*t1 - rr/v
                  dtkal = (wxlsq-p) / (wxl*v)
                  dtkcl = (uzl-gk) / v
                  s = gk*b(l) - gl*uzl
                  t1 = 2.0d0*gk - uzl
                  t2 = rrsq - bgl
                  dtlal = -(risql*wxlsq*b(l)*t1
     &                         -s*(wxlsq*t2+risql*bsql))
     &                             / (risql*wxl*bsql*v)
                  dtlcl = -(risql*b(l)*(uzl*t1-bgl)-uzl*t2*s)
     &                              / (risql*bsql*v)
                  gaca = rcn * (deal-(gk*dtkal-gl*dtlal)/rr) / wxl
                  gacb = (gk-uzl*gl/b(l)) * sign_yder(l) * rr / wxlsq
                  sign_yder(l) = 0
                  if (.not. moved) then
                     faca = ux(l)*gaca - uy(l)*gacb
                     facb = uy(l)*gaca + ux(l)*gacb
                     facc = rcn * (decl-(gk*dtkcl-gl*dtlcl)/rr)
                     dax = axx*faca - ayx*facb + azx*facc
                     day = axy*faca + ayy*facb + azy*facc
                     daz = azz*facc - axz*faca
                     in = intag(l)
                     darea(1,ir) = darea(1,ir) + dax*wght
                     darea(2,ir) = darea(2,ir) + day*wght
                     darea(3,ir) = darea(3,ir) + daz*wght
                     darea(1,in) = darea(1,in) - dax*wght
                     darea(2,in) = darea(2,in) - day*wght
                     darea(3,in) = darea(3,in) - daz*wght
                  end if
               end if
            end do
            goto 100
   90       continue
            arcsum = pix2
            ib = ib + 1
  100       continue
            arclen = arclen + gr(k)*arcsum
            if (.not. moved) then
               in = intag(k)
               t1 = arcsum*rrsq*(bsqk-rrsq+r(in)**2) / (rrx2*bsqk*bk)
               darea(1,ir) = darea(1,ir) - txk*t1*wght
               darea(2,ir) = darea(2,ir) - tyk*t1*wght
               darea(3,ir) = darea(3,ir) - tzk*t1*wght
               darea(1,in) = darea(1,in) + txk*t1*wght
               darea(2,in) = darea(2,in) + tyk*t1*wght
               darea(3,in) = darea(3,in) + tzk*t1*wght
            end if
  110       continue
         end do
         if (arclen .eq. 0.0d0)  goto 180
         if (jb .eq. 0)  goto 150
c
c     find number of independent boundaries and check connectivity
c
         j = 0
         do k = 1, jb
            if (kout(k) .ne. 0) then
               i = k
  120          continue
               m = kout(i)
               kout(i) = 0
               j = j + 1
               do ii = 1, jb
                  if (m .eq. kent(ii)) then
                     if (ii .eq. k) then
                        ib = ib + 1
                        if (j .eq. jb)  goto 150
                        goto 130
                     end if
                     i = ii
                     goto 120
                  end if
               end do
  130          continue
            end if
         end do
         ib = ib + 1
c
c     attempt to fix connectivity error by moving atom slightly
c
         if (moved) then
            write (iout,140)  ir
  140       format (/,' SURFACE1  --  Connectivity Error at Atom',i6)
         else
            moved = .true.
            xr = xr + rmove
            yr = yr + rmove
            zr = zr + rmove
            goto 10
         end if
c
c     form the accessible area for the current atom
c
  150    continue
         area(ir) = ib*pix2 + exang + arclen
         area(ir) = mod(area(ir),pix4)
  160    continue
         area(ir) = area(ir) * rrsq
c
c     attempt to fix negative area by moving atom slightly
c
         if (area(ir) .lt. 0.0d0) then
            if (moved) then
               write (iout,170)  ir
  170          format (/,' SURFACE1  --  Negative Area at Atom',i6)
            else
               moved = .true.
               xr = xr + rmove
               yr = yr + rmove
               zr = zr + rmove
               goto 10
            end if
         end if
c
c     weight the accessible area by the scale factor
c
         area(ir) = area(ir) * wght
         total = total + area(ir)
  180    continue
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (r)
      deallocate (skip)
c
c     zero out the area derivatives for the inactive atoms
c
      do i = 1, n
         if (.not. use(i)) then
            darea(1,i) = 0.0d0
            darea(2,i) = 0.0d0
            darea(3,i) = 0.0d0
         end if
      end do
c
c     print out the surface area and derivatives for each atom
c
c     if (debug) then
c        write (iout,190)
c 190    format (/,' Weighted Atomic Surface Areas and Derivatives :',
c    &           //,4x,'Atom',7x,'Area Term',10x,'dA/dx',
c    &              7x,'dA/dy',7x,'dA/dz',6x,'Weight',/)
c        do i = 1, n
c           if (.not. skip(i)) then
c              write (iout,200)  i,area(i),(darea(j,i),j=1,3),weight(i)
c 200          format (i8,4x,f12.4,3x,3f12.4,f12.4)
c           end if
c        end do
c        write (iout,210)  total
c 210    format (/,' Total Weighted Surface Area :',5x,f16.4)
c     end if
      return
      end
