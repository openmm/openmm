c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1996  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine getref  --  get structure from reference area  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "getref" copies structure information from the reference area
c     into the standard variables for the current system structure
c
c
      subroutine getref (iref)
      use sizes
      use atomid
      use atoms
      use boxes
      use couple
      use files
      use refer
      use titles
      implicit none
      integer i,j,iref
c
c
c     retrieve the filename and title line for the structure
c
      filename = reffile(iref)
      leng = refleng(iref)
      title = reftitle(iref)
      ltitle = refltitle(iref)
c
c     retrieve the coordinates, type and connectivity of each atom
c
      n = nref(iref)
      do i = 1, n
         name(i) = refnam(i,iref)
         x(i) = xref(i,iref)
         y(i) = yref(i,iref)
         z(i) = zref(i,iref)
         type(i) = reftyp(i,iref)
         n12(i) = n12ref(i,iref)
         do j = 1, n12(i)
            i12(j,i) = i12ref(j,i,iref)
         end do
      end do
c
c     retrieve any unitcell parameters defining a periodic box
c
      xbox = xboxref(iref)
      ybox = yboxref(iref)
      zbox = zboxref(iref)
      alpha = alpharef(iref)
      beta = betaref(iref)
      gamma = gammaref(iref)
      return
      end
