c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine xyzatm  --  single atom internal to Cartesian  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "xyzatm" computes the Cartesian coordinates of a single
c     atom from its defining internal coordinate values
c
c
      subroutine xyzatm (i,ia,bond,ib,angle1,ic,angle2,chiral)
      use sizes
      use atoms
      use inform
      use iounit
      use math
      implicit none
      integer i,ia,ib,ic,chiral
      real*8 bond,angle1,angle2
      real*8 eps,rad1,rad2
      real*8 sin1,cos1,sin2,cos2
      real*8 cosine,sine,sine2
      real*8 xab,yab,zab,rab
      real*8 xba,yba,zba,rba
      real*8 xbc,ybc,zbc,rbc
      real*8 xac,yac,zac,rac
      real*8 xt,yt,zt,xu,yu,zu
      real*8 cosb,sinb,cosg,sing
      real*8 xtmp,ztmp,a,b,c
c
c
c     convert angles to radians, and get their sines and cosines
c
      eps = 0.00000001d0
      rad1 = angle1 / radian
      rad2 = angle2 / radian
      sin1 = sin(rad1)
      cos1 = cos(rad1)
      sin2 = sin(rad2)
      cos2 = cos(rad2)
c
c     if no second site given, place the atom at the origin
c
      if (ia .eq. 0) then
         x(i) = 0.0d0
         y(i) = 0.0d0
         z(i) = 0.0d0
c
c     if no third site given, place the atom along the z-axis
c
      else if (ib .eq. 0) then
         x(i) = x(ia)
         y(i) = y(ia)
         z(i) = z(ia) + bond
c
c     if no fourth site given, place the atom in the x,z-plane
c
      else if (ic .eq. 0) then
         xab = x(ia) - x(ib)
         yab = y(ia) - y(ib)
         zab = z(ia) - z(ib)
         rab = sqrt(xab**2 + yab**2 + zab**2)
         xab = xab / rab
         yab = yab / rab
         zab = zab / rab
         cosb = zab
         sinb = sqrt(xab**2 + yab**2)
         if (sinb .eq. 0.0d0) then
            cosg = 1.0d0
            sing = 0.0d0
         else
            cosg = yab / sinb
            sing = xab / sinb
         end if
         xtmp = bond*sin1
         ztmp = rab - bond*cos1
         x(i) = x(ib) + xtmp*cosg + ztmp*sing*sinb
         y(i) = y(ib) - xtmp*sing + ztmp*cosg*sinb
         z(i) = z(ib) + ztmp*cosb
c
c     general case where the second angle is a dihedral angle
c
      else if (chiral .eq. 0) then
         xab = x(ia) - x(ib)
         yab = y(ia) - y(ib)
         zab = z(ia) - z(ib)
         rab = sqrt(xab**2 + yab**2 + zab**2)
         xab = xab / rab
         yab = yab / rab
         zab = zab / rab
         xbc = x(ib) - x(ic)
         ybc = y(ib) - y(ic)
         zbc = z(ib) - z(ic)
         rbc = sqrt(xbc**2 + ybc**2 + zbc**2)
         xbc = xbc / rbc
         ybc = ybc / rbc
         zbc = zbc / rbc
         xt = zab*ybc - yab*zbc
         yt = xab*zbc - zab*xbc
         zt = yab*xbc - xab*ybc
         cosine = xab*xbc + yab*ybc + zab*zbc
         sine = sqrt(max(1.0d0-cosine**2,eps))
         xt = xt / sine
         yt = yt / sine
         zt = zt / sine
         xu = yt*zab - zt*yab
         yu = zt*xab - xt*zab
         zu = xt*yab - yt*xab
         x(i) = x(ia) + bond * (xu*sin1*cos2 + xt*sin1*sin2 - xab*cos1)
         y(i) = y(ia) + bond * (yu*sin1*cos2 + yt*sin1*sin2 - yab*cos1)
         z(i) = z(ia) + bond * (zu*sin1*cos2 + zt*sin1*sin2 - zab*cos1)
         if (abs(cosine) .ge. 1.0d0) then
            cosb = zab
            sinb = sqrt(xab**2 + yab**2)
            if (sinb .eq. 0.0d0) then
               cosg = 1.0d0
               sing = 0.0d0
            else
               cosg = yab / sinb
               sing = xab / sinb
            end if
            xtmp = bond*sin1
            ztmp = rab - bond*cos1
            x(i) = x(ib) + xtmp*cosg + ztmp*sing*sinb
            y(i) = y(ib) - xtmp*sing + ztmp*cosg*sinb
            z(i) = z(ib) + ztmp*cosb
            write (iout,10)  i
   10       format (/,' XYZATM  --  Warning, Undefined Dihedral',
     &                 ' Angle at Atom',i6)
         end if
c
c     general case where the second angle is a bond angle
c
      else if (abs(chiral) .eq. 1) then
         xba = x(ib) - x(ia)
         yba = y(ib) - y(ia)
         zba = z(ib) - z(ia)
         rba = sqrt(xba**2 + yba**2 + zba**2)
         xba = xba / rba
         yba = yba / rba
         zba = zba / rba
         xac = x(ia) - x(ic)
         yac = y(ia) - y(ic)
         zac = z(ia) - z(ic)
         rac = sqrt(xac**2 + yac**2 + zac**2)
         xac = xac / rac
         yac = yac / rac
         zac = zac / rac
         xt = zba*yac - yba*zac
         yt = xba*zac - zba*xac
         zt = yba*xac - xba*yac
         cosine = xba*xac + yba*yac + zba*zac
         sine2 = max(1.0d0-cosine**2,eps)
         if (abs(cosine) .ge. 1.0d0) then
            write (iout,20)  i
   20       format (/,' XYZATM  --  Warning, Collinear Defining',
     &                 ' Atoms at Atom',i6)
         end if
         a = (-cos2 - cosine*cos1) / sine2
         b = (cos1 + cosine*cos2) / sine2
         c = (1.0d0 + a*cos2 - b*cos1) / sine2
         if (c .gt. eps) then
            c = chiral * sqrt(c)
         else if (c .lt. -eps) then
            c = sqrt((a*xac+b*xba)**2 + (a*yac+b*yba)**2
     &                       + (a*zac+b*zba)**2)
            a = a / c
            b = b / c
            c = 0.0d0
            if (debug) then
               write (iout,30)  ia
   30          format (/,' XYZATM  --  Warning, Sum of Bond Angles',
     &                    ' Too Large at Atom',i6)
            end if
         else
            c = 0.0d0
         end if
         x(i) = x(ia) + bond * (a*xac + b*xba + c*xt)
         y(i) = y(ia) + bond * (a*yac + b*yba + c*yt)
         z(i) = z(ia) + bond * (a*zac + b*zba + c*zt)
      end if
      return
      end
