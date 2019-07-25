c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1996  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine surfatom  --  exposed surface area of an atom  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "surfatom" performs an analytical computation of the surface
c     area of a specified atom; a simplified version of "surface"
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
c     ir       number of atom for which area is desired
c     area     accessible surface area of the atom
c     radius   radii of each of the individual atoms
c
c
      subroutine surfatom (ir,area,radius)
      use sizes
      use atoms
      use iounit
      use math
      implicit none
      integer maxarc
      parameter (maxarc=1000)
      integer i,j,k,m
      integer ii,ib,jb
      integer io,ir
      integer mi,ni,narc
      integer key(maxarc)
      integer intag(maxarc)
      integer intag1(maxarc)
      integer lt(maxarc)
      integer kent(maxarc)
      integer kout(maxarc)
      real*8 area,arcsum
      real*8 arclen,exang
      real*8 delta,delta2
      real*8 eps,rmove
      real*8 xr,yr,zr
      real*8 rr,rrsq
      real*8 rplus,rminus
      real*8 axx,axy,axz
      real*8 ayx,ayy
      real*8 azx,azy,azz
      real*8 uxj,uyj,uzj
      real*8 tx,ty,tz
      real*8 txb,tyb,td
      real*8 tr2,tr,txr,tyr
      real*8 tk1,tk2
      real*8 thec,the,t,tb
      real*8 txk,tyk,tzk
      real*8 t1,ti,tf,tt
      real*8 txj,tyj,tzj
      real*8 ccsq,cc,xysq
      real*8 bsqk,bk,cosine
      real*8 dsqj,gi,pix2
      real*8 therk,dk,gk
      real*8 risqk,rik
      real*8 radius(*)
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
      logical moved,top
      logical omit(maxarc)
c
c
c     zero out the surface area for the sphere of interest
c
      area = 0.0d0
      if (radius(ir) .eq. 0.0d0)  return
c
c     set the overlap significance and connectivity shift
c
      pix2 = 2.0d0 * pi
      delta = 1.0d-8
      delta2 = delta * delta
      eps = 1.0d-8
      moved = .false.
      rmove = 1.0d-8
c
c     store coordinates and radius of the sphere of interest
c
      xr = x(ir)
      yr = y(ir)
      zr = z(ir)
      rr = radius(ir)
      rrsq = rr * rr
c
c     initialize values of some counters and summations
c
   10 continue
      io = 0
      jb = 0
      ib = 0
      arclen = 0.0d0
      exang = 0.0d0
c
c     test each sphere to see if it overlaps the sphere of interest
c
      do i = 1, n
         if (i.eq.ir .or. radius(i).eq.0.0d0)  goto 30
         rplus = rr + radius(i)
         tx = x(i) - xr
         if (abs(tx) .ge. rplus)  goto 30
         ty = y(i) - yr
         if (abs(ty) .ge. rplus)  goto 30
         tz = z(i) - zr
         if (abs(tz) .ge. rplus)  goto 30
c
c     check for sphere overlap by testing distance against radii
c
         xysq = tx*tx + ty*ty
         if (xysq .lt. delta2) then
            tx = delta
            ty = 0.0d0
            xysq = delta2
         end if
         ccsq = xysq + tz*tz
         cc = sqrt(ccsq)
         if (rplus-cc .le. delta)  goto 30
         rminus = rr - radius(i)
c
c     check to see if sphere of interest is completely buried
c
         if (cc-abs(rminus) .le. delta) then
            if (rminus .le. 0.0d0)  goto 170
            goto 30
         end if
c
c     check for too many overlaps with sphere of interest
c
         if (io .ge. maxarc) then
            write (iout,20)
   20       format (/,' SURFATOM  --  Increase the Value of MAXARC')
            call fatal
         end if
c
c     get overlap between current sphere and sphere of interest
c
         io = io + 1
         xc1(io) = tx
         yc1(io) = ty
         zc1(io) = tz
         dsq1(io) = xysq
         bsq1(io) = ccsq
         b1(io) = cc
         gr(io) = (ccsq+rplus*rminus) / (2.0d0*rr*b1(io))
         intag1(io) = i
         omit(io) = .false.
   30    continue
      end do
c
c     case where no other spheres overlap the sphere of interest
c
      if (io .eq. 0) then
         area = 4.0d0 * pi * rrsq
         return
      end if
c
c     case where only one sphere overlaps the sphere of interest
c
      if (io .eq. 1) then
         area = pix2 * (1.0d0 + gr(1))
         area = mod(area,4.0d0*pi) * rrsq
         return
      end if
c
c     case where many spheres intersect the sphere of interest;
c     sort the intersecting spheres by their degree of overlap
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
      end do
c
c     get radius of each overlap circle on surface of the sphere
c
      do i = 1, io
         gi = gr(i) * rr
         bg(i) = b(i) * gi
         risq(i) = rrsq - gi*gi
         ri(i) = sqrt(risq(i))
         ther(i) = 0.5d0*pi - asin(min(1.0d0,max(-1.0d0,gr(i))))
      end do
c
c     find boundary of inaccessible area on sphere of interest
c
      do k = 1, io-1
         if (.not. omit(k)) then
            txk = xc(k)
            tyk = yc(k)
            tzk = zc(k)
            bk = b(k)
            therk = ther(k)
c
c     check to see if J circle is intersecting K circle;
c     get distance between circle centers and sum of radii
c
            do j = k+1, io
               if (omit(j))  goto 60
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
c     check for circles that are essentially parallel
c
               if (cc .gt. delta)  goto 50
   40          continue
               omit(j) = .true.
               goto 60
c
c     check to see if sphere of interest is completely buried
c
   50          continue
               if (pix2-cc .le. td)  goto 170
   60          continue
            end do
         end if
      end do
c
c     find T value of circle intersections
c
      do k = 1, io
         if (omit(k))  goto 110
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
         do j = 1, io
            if (.not. omit(j)) then
               txj = xc(j)
               tyj = yc(j)
               tzj = zc(j)
c
c     rotate spheres so K vector colinear with z-axis
c
               uxj = txj*axx + tyj*axy - tzj*axz
               uyj = tyj*ayy - txj*ayx
               uzj = txj*azx + tyj*azy + tzj*azz
               cosine = min(1.0d0,max(-1.0d0,uzj/b(j)))
               if (acos(cosine) .lt. therk+ther(j)) then
                  dsqj = uxj*uxj + uyj*uyj
                  tb = uzj*gk - bg(j)
                  txb = uxj * tb
                  tyb = uyj * tb
                  td = rik * dsqj
                  tr2 = risqk*dsqj - tb*tb
                  tr2 = max(eps,tr2)
                  tr = sqrt(tr2)
                  txr = uxj * tr
                  tyr = uyj * tr
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
                  thec = (rrsq*uzj-gk*bg(j)) / (rik*ri(j)*b(j))
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
     &                            (uzj*gk-uxj*rik)/(b(j)*rr)))
                  if ((acos(cosine)-ther(j))*(tk2-tk1) .le. 0.0d0) then
                     ti = tk2
                     tf = tk1
                  else
                     ti = tk2
                     tf = tk1
                  end if
                  narc = narc + 1
                  if (narc .ge. maxarc) then
                     write (iout,70)
   70                format (/,' SURFATOM  --  Increase the Value',
     &                          ' of MAXARC')
                     call fatal
                  end if
                  if (tf .le. ti) then
                     arcf(narc) = tf
                     arci(narc) = 0.0d0
                     tf = pix2
                     lt(narc) = j
                     ex(narc) = the
                     top = .true.
                     narc = narc + 1
                  end if
                  arcf(narc) = tf
                  arci(narc) = ti
                  lt(narc) = j
                  ex(narc) = the
                  ux(j) = uxj
                  uy(j) = uyj
                  uz(j) = uzj
               end if
            end if
         end do
         omit(k) = .false.
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
   80                format (/,' SURFATOM  --  Increase the Value',
     &                          ' of MAXARC')
                     call fatal
                  end if
                  i = lt(ni)
                  kent(jb) = maxarc*i + k
                  i = lt(m)
                  kout(jb) = maxarc*k + i
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
            i = lt(ni)
            kent(jb) = maxarc*i + k
            i = lt(mi)
            kout(jb) = maxarc*k + i
         end if
         goto 100
   90    continue
         arcsum = pix2
         ib = ib + 1
  100    continue
         arclen = arclen + gr(k)*arcsum
  110    continue
      end do
      if (arclen .eq. 0.0d0)  goto 170
      if (jb .eq. 0)  goto 150
c
c     find number of independent boundaries and check connectivity
c
      j = 0
      do k = 1, jb
         if (kout(k) .ne. 0) then
            i = k
  120       continue
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
  130       continue
         end if
      end do
      ib = ib + 1
c
c     attempt to fix connectivity error by moving atom slightly
c
      if (moved) then
         write (iout,140)  ir
  140    format (/,' SURFATOM  --  Connectivity Error at Atom',i6)
      else
         moved = .true.
         xr = xr + rmove
         yr = yr + rmove
         zr = zr + rmove
         goto 10
      end if
c
c     compute the exposed surface area for the sphere of interest
c
  150 continue
      area = ib*pix2 + exang + arclen
      area = mod(area,4.0d0*pi) * rrsq
c
c     attempt to fix negative area by moving atom slightly
c
      if (area .lt. 0.0d0) then
         if (moved) then
            write (iout,160)  ir
  160       format (/,' SURFATOM  --  Negative Area at Atom',i6)
         else
            moved = .true.
            xr = xr + rmove
            yr = yr + rmove
            zr = zr + rmove
            goto 10
         end if
      end if
  170 continue
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine surfatom1  --  surface area and derivs of atom  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "surfatom1" performs an analytical computation of the surface
c     area and first derivatives with respect to Cartesian coordinates
c     of a specified atom
c
c
      subroutine surfatom1 (ir,area,darea,radius)
      use sizes
      use atoms
      use iounit
      use math
      implicit none
      integer maxarc
      parameter (maxarc=1000)
      integer i,j,k,m
      integer ii,ib,jb
      integer io,ir,in
      integer mi,ni,narc
      integer key(maxarc)
      integer intag(maxarc)
      integer intag1(maxarc)
      integer lt(maxarc)
      integer kent(maxarc)
      integer kout(maxarc)
      integer ider(maxarc)
      integer sign_yder(maxarc)
      real*8 area,arcsum
      real*8 arclen,exang
      real*8 delta,delta2
      real*8 wxl,wxlsq
      real*8 p,s,v,rcn
      real*8 eps,rmove
      real*8 xr,yr,zr
      real*8 rr,rin
      real*8 rrx2,rrsq
      real*8 rplus,rminus
      real*8 axx,axy,axz
      real*8 ayx,ayy
      real*8 azx,azy,azz
      real*8 uxj,uyj,uzj
      real*8 tx,ty,tz
      real*8 txb,tyb,td
      real*8 tr2,tr,txr,tyr
      real*8 tk1,tk2
      real*8 thec,the,t,tb
      real*8 txk,tyk,tzk
      real*8 t1,ti,tf,tt
      real*8 txj,tyj,tzj
      real*8 ccsq,cc,xysq
      real*8 bgl,bsqk,bsql
      real*8 bk,cosine
      real*8 gl,uzl,t2
      real*8 dsqj,gi,pix2
      real*8 dax,day,daz
      real*8 deal,decl
      real*8 dtkal,dtkcl
      real*8 dtlal,dtlcl
      real*8 therk,dk,gk
      real*8 risqk,rik,risql
      real*8 faca,facb,facc
      real*8 gaca,gacb
      real*8 radius(*)
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
      logical moved,top
      logical omit(maxarc)
c
c
c     zero out the area and derivatives for sphere of interest
c
      area = 0.0d0
      do i = 1, n
         darea(1,i) = 0.0d0
         darea(2,i) = 0.0d0
         darea(3,i) = 0.0d0
      end do
      if (radius(ir) .eq. 0.0d0)  return
c
c     set the overlap significance and connectivity shift
c
      pix2 = 2.0d0 * pi
      delta = 1.0d-8
      delta2 = delta * delta
      eps = 1.0d-8
      moved = .false.
      rmove = 1.0d-8
      do i = 1, maxarc
         ider(i) = 0
         sign_yder(i) = 0
      end do
c
c     store coordinates and radius of the sphere of interest
c
      xr = x(ir)
      yr = y(ir)
      zr = z(ir)
      rr = radius(ir)
      rrx2 = 2.0d0 * rr
      rrsq = rr * rr
c
c     initialize values of some counters and summations
c
   10 continue
      io = 0
      jb = 0
      ib = 0
      arclen = 0.0d0
      exang = 0.0d0
c
c     test each sphere to see if it overlaps the sphere of interest
c
      do i = 1, n
         if (i.eq.ir .or. radius(i).eq.0.0d0)  goto 30
         rplus = rr + radius(i)
         tx = x(i) - xr
         if (abs(tx) .ge. rplus)  goto 30
         ty = y(i) - yr
         if (abs(ty) .ge. rplus)  goto 30
         tz = z(i) - zr
         if (abs(tz) .ge. rplus)  goto 30
c
c     check for sphere overlap by testing distance against radii
c
         xysq = tx*tx + ty*ty
         if (xysq .lt. delta2) then
            tx = delta
            ty = 0.0d0
            xysq = delta2
         end if
         ccsq = xysq + tz*tz
         cc = sqrt(ccsq)
         if (rplus-cc .le. delta)  goto 30
         rminus = rr - radius(i)
c
c     check to see if sphere of interest is completely buried
c
         if (cc-abs(rminus) .le. delta) then
            if (rminus .le. 0.0d0)  goto 170
            goto 30
         end if
c
c     check for too many overlaps with sphere of interest
c
         if (io .ge. maxarc) then
            write (iout,20)
   20       format (/,' SURFATOM  --  Increase the Value of MAXARC')
            call fatal
         end if
c
c     get overlap between current sphere and sphere of interest
c
         io = io + 1
         xc1(io) = tx
         yc1(io) = ty
         zc1(io) = tz
         dsq1(io) = xysq
         bsq1(io) = ccsq
         b1(io) = cc
         gr(io) = (ccsq+rplus*rminus) / (2.0d0*rr*b1(io))
         intag1(io) = i
         omit(io) = .false.
   30    continue
      end do
c
c     case where no other spheres overlap the sphere of interest
c
      if (io .eq. 0) then
         area = 4.0d0 * pi * rrsq
         return
      end if
c
c     case where only one sphere overlaps the sphere of interest
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
            rin = radius(in)
            t1 = arcsum*rrsq*(bsqk-rrsq+rin*rin) / (rrx2*bsqk*bk)
            darea(1,ir) = darea(1,ir) - txk*t1
            darea(2,ir) = darea(2,ir) - tyk*t1
            darea(3,ir) = darea(3,ir) - tzk*t1
            darea(1,in) = darea(1,in) + txk*t1
            darea(2,in) = darea(2,in) + tyk*t1
            darea(3,in) = darea(3,in) + tzk*t1
         end if
         area = pix2 * (1.0d0 + gr(1))
         area = mod(area,4.0d0*pi) * rrsq
         return
      end if
c
c     case where many spheres intersect the sphere of interest;
c     sort the intersecting spheres by their degree of overlap
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
      end do
c
c     get radius of each overlap circle on surface of the sphere
c
      do i = 1, io
         gi = gr(i) * rr
         bg(i) = b(i) * gi
         risq(i) = rrsq - gi*gi
         ri(i) = sqrt(risq(i))
         ther(i) = 0.5d0*pi - asin(min(1.0d0,max(-1.0d0,gr(i))))
      end do
c
c     find boundary of inaccessible area on sphere of interest
c
      do k = 1, io-1
         if (.not. omit(k)) then
            txk = xc(k)
            tyk = yc(k)
            tzk = zc(k)
            bk = b(k)
            therk = ther(k)
c
c     check to see if J circle is intersecting K circle;
c     get distance between circle centers and sum of radii
c
            do j = k+1, io
               if (omit(j))  goto 60
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
c     check for circles that are essentially parallel
c
               if (cc .gt. delta)  goto 50
   40          continue
               omit(j) = .true.
               goto 60
c
c     check to see if sphere of interest is completely buried
c
   50          continue
               if (pix2-cc .le. td)  goto 170
   60          continue
            end do
         end if
      end do
c
c     find T value of circle intersections
c
      do k = 1, io
         if (omit(k))  goto 110
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
         do j = 1, io
            if (.not. omit(j)) then
               txj = xc(j)
               tyj = yc(j)
               tzj = zc(j)
c
c     rotate spheres so K vector colinear with z-axis
c
               uxj = txj*axx + tyj*axy - tzj*axz
               uyj = tyj*ayy - txj*ayx
               uzj = txj*azx + tyj*azy + tzj*azz
               cosine = min(1.0d0,max(-1.0d0,uzj/b(j)))
               if (acos(cosine) .lt. therk+ther(j)) then
                  dsqj = uxj*uxj + uyj*uyj
                  tb = uzj*gk - bg(j)
                  txb = uxj * tb
                  tyb = uyj * tb
                  td = rik * dsqj
                  tr2 = risqk*dsqj - tb*tb
                  tr2 = max(eps,tr2)
                  tr = sqrt(tr2)
                  txr = uxj * tr
                  tyr = uyj * tr
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
                  thec = (rrsq*uzj-gk*bg(j)) / (rik*ri(j)*b(j))
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
     &                            (uzj*gk-uxj*rik)/(b(j)*rr)))
                  if ((acos(cosine)-ther(j))*(tk2-tk1) .le. 0.0d0) then
                     ti = tk2
                     tf = tk1
                  else
                     ti = tk2
                     tf = tk1
                  end if
                  narc = narc + 1
                  if (narc .ge. maxarc) then
                     write (iout,70)
   70                format (/,' SURFATOM  --  Increase the Value',
     &                          ' of MAXARC')
                     call fatal
                  end if
                  if (tf .le. ti) then
                     arcf(narc) = tf
                     arci(narc) = 0.0d0
                     tf = pix2
                     lt(narc) = j
                     ex(narc) = the
                     top = .true.
                     narc = narc + 1
                  end if
                  arcf(narc) = tf
                  arci(narc) = ti
                  lt(narc) = j
                  ex(narc) = the
                  ux(j) = uxj
                  uy(j) = uyj
                  uz(j) = uzj
               end if
            end if
         end do
         omit(k) = .false.
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
   80                format (/,' SURFATOM  --  Increase the Value',
     &                          ' of MAXARC')
                     call fatal
                  end if
                  i = lt(ni)
                  ider(i) = ider(i) + 1
                  sign_yder(i) = sign_yder(i) + 1
                  kent(jb) = maxarc*i + k
                  i = lt(m)
                  ider(i) = ider(i) + 1
                  sign_yder(i) = sign_yder(i) - 1
                  kout(jb) = maxarc*k + i
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
            i = lt(ni)
            ider(i) = ider(i) + 1
            sign_yder(i) = sign_yder(i) + 1
            kent(jb) = maxarc*i + k
            i = lt(mi)
            ider(i) = ider(i) + 1
            sign_yder(i) = sign_yder(i) - 1
            kout(jb) = maxarc*k + i
         end if
c
c     calculate the surface area derivatives
c
         do j = 1, io
            if (ider(j) .ne. 0) then
               rcn = ider(j) * rrsq
               ider(j) = 0
               uzl = uz(j)
               gl = gr(j) * rr
               bgl = bg(j)
               bsql = bsq(j)
               risql = risq(j)
               wxlsq = bsql - uzl**2
               wxl = sqrt(wxlsq)
               p = bgl - gk*uzl
               v = risqk*wxlsq - p**2
               v = max(eps,v)
               v = sqrt(v)
               t1 = rr * (gk*(bgl-bsql)+uzl*(bgl-rrsq))
     &                          / (v*risql*bsql)
               deal = -wxl*t1
               decl = -uzl*t1 - rr/v
               dtkal = (wxlsq-p) / (wxl*v)
               dtkcl = (uzl-gk) / v
               s = gk*b(j) - gl*uzl
               t1 = 2.0d0*gk - uzl
               t2 = rrsq - bgl
               dtlal = -(risql*wxlsq*b(j)*t1
     &                      -s*(wxlsq*t2+risql*bsql))
     &                          / (risql*wxl*bsql*v)
               dtlcl = -(risql*b(j)*(uzl*t1-bgl)-uzl*t2*s)
     &                          / (risql*bsql*v)
               gaca = rcn * (deal-(gk*dtkal-gl*dtlal)/rr) / wxl
               gacb = (gk-uzl*gl/b(j)) * sign_yder(j) * rr / wxlsq
               sign_yder(j) = 0
               if (.not. moved) then
                  faca = ux(j)*gaca - uy(j)*gacb
                  facb = uy(j)*gaca + ux(j)*gacb
                  facc = rcn * (decl-(gk*dtkcl-gl*dtlcl)/rr)
                  dax = axx*faca - ayx*facb + azx*facc
                  day = axy*faca + ayy*facb + azy*facc
                  daz = azz*facc - axz*faca
                  in = intag(j)
                  darea(1,ir) = darea(1,ir) + dax
                  darea(2,ir) = darea(2,ir) + day
                  darea(3,ir) = darea(3,ir) + daz
                  darea(1,in) = darea(1,in) - dax
                  darea(2,in) = darea(2,in) - day
                  darea(3,in) = darea(3,in) - daz
               end if
            end if
         end do
         goto 100
   90    continue
         arcsum = pix2
         ib = ib + 1
  100    continue
         arclen = arclen + gr(k)*arcsum
         if (.not. moved) then
            in = intag(k)
            rin = radius(in)
            t1 = arcsum*rrsq*(bsqk-rrsq+rin*rin) / (rrx2*bsqk*bk)
            darea(1,ir) = darea(1,ir) - txk*t1
            darea(2,ir) = darea(2,ir) - tyk*t1
            darea(3,ir) = darea(3,ir) - tzk*t1
            darea(1,in) = darea(1,in) + txk*t1
            darea(2,in) = darea(2,in) + tyk*t1
            darea(3,in) = darea(3,in) + tzk*t1
         end if
  110    continue
      end do
      if (arclen .eq. 0.0d0)  goto 170
      if (jb .eq. 0)  goto 150
c
c     find number of independent boundaries and check connectivity
c
      j = 0
      do k = 1, jb
         if (kout(k) .ne. 0) then
            i = k
  120       continue
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
  130       continue
         end if
      end do
      ib = ib + 1
c
c     attempt to fix connectivity error by moving atom slightly
c
      if (moved) then
         write (iout,140)  ir
  140    format (/,' SURFATOM  --  Connectivity Error at Atom',i6)
      else
         moved = .true.
         xr = xr + rmove
         yr = yr + rmove
         zr = zr + rmove
         goto 10
      end if
c
c     compute the exposed surface area for the sphere of interest
c
  150 continue
      area = ib*pix2 + exang + arclen
      area = mod(area,4.0d0*pi) * rrsq
c
c     attempt to fix negative area by moving atom slightly
c
      if (area .lt. 0.0d0) then
         if (moved) then
            write (iout,160)  ir
  160       format (/,' SURFATOM  --  Negative Area at Atom',i6)
         else
            moved = .true.
            xr = xr + rmove
            yr = yr + rmove
            zr = zr + rmove
            goto 10
         end if
      end if
  170 continue
      return
      end
