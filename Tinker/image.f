c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine image  --  compute the minimum image distance  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "image" takes the components of pairwise distance between
c     two points in a periodic box and converts to the components
c     of the minimum image distance
c
c
      subroutine image (xr,yr,zr)
      use boxes
      use cell
      implicit none
      real*8 xr,yr,zr
c
c
c     for orthogonal lattice, find the desired image directly
c
      if (orthogonal) then
         do while (abs(xr) .gt. xcell2)
            xr = xr - sign(xcell,xr)
         end do
         do while (abs(yr) .gt. ycell2)
            yr = yr - sign(ycell,yr)
         end do
         do while (abs(zr) .gt. zcell2)
            zr = zr - sign(zcell,zr)
         end do
c
c     for monoclinic lattice, convert "xr" and "zr" to
c     fractional coordinates, find desired image and then
c     translate fractional coordinates back to Cartesian
c
      else if (monoclinic) then
         zr = zr / beta_sin
         xr = xr - zr*beta_cos
         do while (abs(xr) .gt. xcell2)
            xr = xr - sign(xcell,xr)
         end do
         do while (abs(yr) .gt. ycell2)
            yr = yr - sign(ycell,yr)
         end do
         do while (abs(zr) .gt. zcell2)
            zr = zr - sign(zcell,zr)
         end do
         xr = xr + zr*beta_cos
         zr = zr * beta_sin
c
c     for triclinic lattice, convert pairwise components to
c     fractional coordinates, find desired image and then
c     translate fractional coordinates back to Cartesian
c
      else if (triclinic) then
         zr = zr / gamma_term
         yr = (yr - zr*beta_term) / gamma_sin
         xr = xr - yr*gamma_cos - zr*beta_cos
         do while (abs(xr) .gt. xcell2)
            xr = xr - sign(xcell,xr)
         end do
         do while (abs(yr) .gt. ycell2)
            yr = yr - sign(ycell,yr)
         end do
         do while (abs(zr) .gt. zcell2)
            zr = zr - sign(zcell,zr)
         end do
         xr = xr + yr*gamma_cos + zr*beta_cos
         yr = yr*gamma_sin + zr*beta_term
         zr = zr * gamma_term
c
c     for truncated octahedron, use orthogonal box equations,
c     then perform extra tests to remove corner pieces
c
      else if (octahedron) then
         do while (abs(xr) .gt. xbox2)
            xr = xr - sign(xbox,xr)
         end do
         do while (abs(yr) .gt. ybox2)
            yr = yr - sign(ybox,yr)
         end do
         do while (abs(zr) .gt. zbox2)
            zr = zr - sign(zbox,zr)
         end do
         if (abs(xr)+abs(yr)+abs(zr) .gt. box34) then
            xr = xr - sign(xbox2,xr)
            yr = yr - sign(ybox2,yr)
            zr = zr - sign(zbox2,zr)
         end if
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine imager  --  replicate minimum image distance  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "imager" takes the components of pairwise distance between
c     two points in the same or neighboring periodic boxes and
c     converts to the components of the minimum image distance
c
c
      subroutine imager (xr,yr,zr,i)
      use sizes
      use boxes
      use cell
      implicit none
      integer i
      real*8 xr,yr,zr
      real*8 xsize,ysize,zsize
      real*8 xsize2,ysize2,zsize2
      real*8 xmove,ymove,zmove
c
c
c     set dimensions for either single box or replicated cell
c
      if (i .ge. 0) then
         xsize = xcell
         ysize = ycell
         zsize = zcell
         xsize2 = xcell2
         ysize2 = ycell2
         zsize2 = zcell2
      else
         xsize = xbox
         ysize = ybox
         zsize = zbox
         xsize2 = xbox2
         ysize2 = ybox2
         zsize2 = zbox2
      end if
c
c     compute the distance to translate along each cell axis
c
      if (i .le. 0) then
         xmove = 0.0d0
         ymove = 0.0d0
         zmove = 0.0d0
      else
         xmove = icell(1,i) * xbox
         ymove = icell(2,i) * ybox
         zmove = icell(3,i) * zbox
      end if
c
c     for orthogonal lattice, find the desired image directly
c
      if (orthogonal) then
         xr = xr + xmove
         do while (abs(xr) .gt. xsize2)
            xr = xr - sign(xsize,xr)
         end do
         yr = yr + ymove
         do while (abs(yr) .gt. ysize2)
            yr = yr - sign(ysize,yr)
         end do
         zr = zr + zmove
         do while (abs(zr) .gt. zsize2)
            zr = zr - sign(zsize,zr)
         end do
c
c     for monoclinic lattice, convert "xr" and "zr" to
c     fractional coordinates, find desired image and then
c     translate fractional coordinates back to Cartesian
c
      else if (monoclinic) then
         zr = zr / beta_sin
         xr = xr - zr*beta_cos
         xr = xr + xmove
         do while (abs(xr) .gt. xsize2)
            xr = xr - sign(xsize,xr)
         end do
         yr = yr + ymove
         do while (abs(yr) .gt. ysize2)
            yr = yr - sign(ysize,yr)
         end do
         zr = zr + zmove
         do while (abs(zr) .gt. zsize2)
            zr = zr - sign(zsize,zr)
         end do
         xr = xr + zr*beta_cos
         zr = zr * beta_sin
c
c     for triclinic lattice, convert pairwise components to
c     fractional coordinates, find desired image and then
c     translate fractional coordinates back to Cartesian
c
      else if (triclinic) then
         zr = zr / gamma_term
         yr = (yr - zr*beta_term) / gamma_sin
         xr = xr - yr*gamma_cos - zr*beta_cos
         xr = xr + xmove
         do while (abs(xr) .gt. xsize2)
            xr = xr - sign(xsize,xr)
         end do
         yr = yr + ymove
         do while (abs(yr) .gt. ysize2)
            yr = yr - sign(ysize,yr)
         end do
         zr = zr + zmove
         do while (abs(zr) .gt. zsize2)
            zr = zr - sign(zsize,zr)
         end do
         xr = xr + yr*gamma_cos + zr*beta_cos
         yr = yr*gamma_sin + zr*beta_term
         zr = zr * gamma_term
c
c     for truncated octahedron, use orthogonal box equations,
c     then perform extra tests to remove corner pieces
c
      else if (octahedron) then
         do while (abs(xr) .gt. xbox2)
            xr = xr - sign(xbox,xr)
         end do
         do while (abs(yr) .gt. ybox2)
            yr = yr - sign(ybox,yr)
         end do
         do while (abs(zr) .gt. zbox2)
            zr = zr - sign(zbox,zr)
         end do
         if (abs(xr)+abs(yr)+abs(zr) .gt. box34) then
            xr = xr - sign(xbox2,xr)
            yr = yr - sign(ybox2,yr)
            zr = zr - sign(zbox2,zr)
         end if
      end if
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine imagen  --  fast minimum image magnitude  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "imagen" takes the components of pairwise distance between
c     two points and converts to the components of the minimum
c     image distance
c
c     note this is a fast version which only returns the correct
c     component magnitudes for use in computing the 3D distance
c
c
      subroutine imagen (xr,yr,zr)
      use boxes
      implicit none
      real*8 xr,yr,zr
c
c
c     for orthogonal lattice, find the desired image directly
c
      if (orthogonal) then
         xr = abs(xr)
         yr = abs(yr)
         zr = abs(zr)
         if (xr .gt. xbox2)  xr = xr - xbox
         if (yr .gt. ybox2)  yr = yr - ybox
         if (zr .gt. zbox2)  zr = zr - zbox
c
c     for monoclinic lattice, convert "xr" and "zr" specially
c
      else if (monoclinic) then
         zr = zr / beta_sin
         yr = abs(yr)
         xr = xr - zr*beta_cos
         if (abs(xr) .gt. xbox2)  xr = xr - sign(xbox,xr)
         if (yr .gt. ybox2)  yr = yr - ybox
         if (abs(zr) .gt. zbox2)  zr = zr - sign(zbox,zr)
         xr = xr + zr*beta_cos
         zr = zr * beta_sin
c
c     for triclinic lattice, use general conversion equations
c
      else if (triclinic) then
         zr = zr / gamma_term
         yr = (yr - zr*beta_term) / gamma_sin
         xr = xr - yr*gamma_cos - zr*beta_cos
         if (abs(xr) .gt. xbox2)  xr = xr - sign(xbox,xr)
         if (abs(yr) .gt. ybox2)  yr = yr - sign(ybox,yr)
         if (abs(zr) .gt. zbox2)  zr = zr - sign(zbox,zr)
         xr = xr + yr*gamma_cos + zr*beta_cos
         yr = yr*gamma_sin + zr*beta_term
         zr = zr * gamma_term
c
c     for truncated octahedron, remove the corner pieces
c
      else if (octahedron) then
         if (abs(xr) .gt. xbox2)  xr = xr - sign(xbox,xr)
         if (abs(yr) .gt. ybox2)  yr = yr - sign(ybox,yr)
         if (abs(zr) .gt. zbox2)  zr = zr - sign(zbox,zr)
         if (abs(xr)+abs(yr)+abs(zr) .gt. box34) then
            xr = xr - sign(xbox2,xr)
            yr = yr - sign(ybox2,yr)
            zr = zr - sign(zbox2,zr)
         end if
      end if
      return
      end
