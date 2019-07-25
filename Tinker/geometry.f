c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  function geometry  --  evaluate distance, angle, torsion  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "geometry" finds the value of the interatomic distance, angle
c     or dihedral angle defined by two to four input atoms
c
c
      function geometry (ia,ib,ic,id)
      use sizes
      use atoms
      use math
      implicit none
      integer ia,ib,ic,id
      real*8 xab,yab,zab
      real*8 xba,yba,zba
      real*8 xcb,ycb,zcb
      real*8 xdc,ydc,zdc
      real*8 xt,yt,zt
      real*8 xu,yu,zu
      real*8 rab2,rcb2,rabc
      real*8 rt2,ru2,rtru
      real*8 cosine,sign
      real*8 geometry
c
c
c     set default in case atoms are coincident or colinear
c
      geometry = 0.0d0
c
c     compute the value of the distance in angstroms
c
      if (ic .eq. 0) then
         xab = x(ia) - x(ib)
         yab = y(ia) - y(ib)
         zab = z(ia) - z(ib)
         geometry = sqrt(xab*xab + yab*yab + zab*zab)
c
c     compute the value of the angle in degrees
c
      else if (id .eq. 0) then
         xab = x(ia) - x(ib)
         yab = y(ia) - y(ib)
         zab = z(ia) - z(ib)
         xcb = x(ic) - x(ib)
         ycb = y(ic) - y(ib)
         zcb = z(ic) - z(ib)
         rab2 = xab*xab + yab*yab + zab*zab
         rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
         rabc = sqrt(rab2 * rcb2)
         if (rabc .ne. 0.0d0) then
            cosine = (xab*xcb + yab*ycb + zab*zcb) / rabc
            cosine = min(1.0d0,max(-1.0d0,cosine))
            geometry = radian * acos(cosine)
         end if
c
c     compute the value of the dihedral angle in degrees
c
      else
         xba = x(ib) - x(ia)
         yba = y(ib) - y(ia)
         zba = z(ib) - z(ia)
         xcb = x(ic) - x(ib)
         ycb = y(ic) - y(ib)
         zcb = z(ic) - z(ib)
         xdc = x(id) - x(ic)
         ydc = y(id) - y(ic)
         zdc = z(id) - z(ic)
         xt = yba*zcb - ycb*zba
         yt = xcb*zba - xba*zcb
         zt = xba*ycb - xcb*yba
         xu = ycb*zdc - ydc*zcb
         yu = xdc*zcb - xcb*zdc
         zu = xcb*ydc - xdc*ycb
         rt2 = xt*xt + yt*yt + zt*zt
         ru2 = xu*xu + yu*yu + zu*zu
         rtru = sqrt(rt2 * ru2)
         if (rtru .ne. 0.0d0) then
            cosine = (xt*xu + yt*yu + zt*zu) / rtru
            cosine = min(1.0d0,max(-1.0d0,cosine))
            geometry = radian * acos(cosine)
            sign = xba*xu + yba*yu + zba*zu
            if (sign .lt. 0.0d0)  geometry = -geometry
         end if
      end if
      return
      end
