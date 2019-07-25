c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1996  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine makeref  --  copy structure to reference area  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "makeref" copies the information contained in the "xyz" file
c     of the current structure into corresponding reference areas
c
c
      subroutine makeref (iref)
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
      logical first
      save first
      data first  / .true. /
c
c
c     perform dynamic allocation of some global arrays
c
      if (first) then
         first = .false.
         if (.not. allocated(reftyp))  allocate (reftyp(maxatm,maxref))
         if (.not. allocated(n12ref))  allocate (n12ref(maxatm,maxref))
         if (.not. allocated(i12ref))
     &      allocate (i12ref(maxval,maxatm,maxref))
         if (.not. allocated(xref))  allocate (xref(maxatm,maxref))
         if (.not. allocated(yref))  allocate (yref(maxatm,maxref))
         if (.not. allocated(zref))  allocate (zref(maxatm,maxref))
         if (.not. allocated(refnam))  allocate (refnam(maxatm,maxref))
      end if
c
c     copy the filename and title line for the structure
c
      reffile(iref) = filename
      refleng(iref) = leng
      reftitle(iref) = title
      refltitle(iref) = ltitle
c
c     copy the coordinates, type and connectivity of each atom
c
      nref(iref) = n
      do i = 1, n
         refnam(i,iref) = name(i)
         xref(i,iref) = x(i)
         yref(i,iref) = y(i)
         zref(i,iref) = z(i)
         reftyp(i,iref) = type(i)
         n12ref(i,iref) = n12(i)
         do j = 1, n12(i)
            i12ref(j,i,iref) = i12(j,i)
         end do
      end do
c
c     copy any unitcell parameters from the coordinates file
c
      xboxref(iref) = xbox
      yboxref(iref) = ybox
      zboxref(iref) = zbox
      alpharef(iref) = alpha
      betaref(iref) = beta
      gammaref(iref) = gamma
      return
      end
