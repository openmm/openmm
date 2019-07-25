c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine lights  --  get neighbors via method of lights  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "lights" computes the set of nearest neighbor interactions
c     using the method of lights algorithm
c
c     note this routine can include each pair only once via setting
c     of the negative x-coordinate boundaries, or it can optionally
c     include each pair in both directions, ie, both (A,B) and (B,A);
c     inclusion of one vs both directions is controlled by "unique"
c
c     literature reference:
c
c     F. Sullivan, R. D. Mountain and J. O'Connell, "Molecular
c     Dynamics on Vector Computers", Journal of Computational
c     Physics, 61, 138-153 (1985)
c
c
      subroutine lights (cutoff,nsite,xsort,ysort,zsort,unique)
      use sizes
      use bound
      use boxes
      use cell
      use iounit
      use light
      implicit none
      integer i,j,k
      integer nsite
      integer extent
      real*8 cutoff,box
      real*8 xcut,ycut,zcut
      real*8 xmove,ymove,zmove
      real*8 xsort(*)
      real*8 ysort(*)
      real*8 zsort(*)
      real*8, allocatable :: xfrac(:)
      real*8, allocatable :: yfrac(:)
      real*8, allocatable :: zfrac(:)
      logical unique
c
c
c     check that maximum number of replicates is not exceeded
c
      if (use_replica) then
         if (xcell2.gt.xbox .or. ycell2.gt.ybox
     &           .or. zcell2.gt.zbox) then
            write (iout,10)
   10       format (/,' LIGHTS  --  Number of Replicas is Too',
     &                 ' Large for Method of Lights')
            call fatal
         end if
      end if
c
c     truncated octahedron periodicity is not handled at present
c
      if (use_bounds) then
         if (octahedron) then
            write (iout,20)
   20       format (/,' LIGHTS  --  Truncated Octahedron not',
     &                 ' Supported by Method of Lights')
            call fatal
         end if
      end if
c
c     set the light width based on input distance cutoff
c
      xcut = cutoff
      ycut = cutoff
      zcut = cutoff
      if (use_bounds) then
         if (monoclinic) then
            zcut = zcut / beta_sin
            xcut = xcut + zcut*abs(beta_cos)
         else if (triclinic) then
            zcut = zcut / gamma_term
            ycut = (ycut + zcut*abs(beta_term)) / gamma_sin
            xcut = xcut + ycut*abs(gamma_cos) + zcut*abs(beta_cos)
         end if
         xcut = min(xcut,xcell2)
         ycut = min(ycut,ycell2)
         zcut = min(zcut,zcell2)
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (xfrac(nsite))
      allocate (yfrac(nsite))
      allocate (zfrac(nsite))
c
c     find fractional coordinates for the unitcell atoms
c
      if (use_bounds) then
         if (orthogonal) then
            do i = 1, nsite
               zfrac(i) = zsort(i)
               yfrac(i) = ysort(i)
               xfrac(i) = xsort(i)
            end do
         else if (monoclinic) then
            do i = 1, nsite
               zfrac(i) = zsort(i) / beta_sin
               yfrac(i) = ysort(i)
               xfrac(i) = xsort(i) - zfrac(i)*beta_cos
            end do
         else if (triclinic) then
            do i = 1, nsite
               zfrac(i) = zsort(i) / gamma_term
               yfrac(i) = (ysort(i) - zfrac(i)*beta_term) / gamma_sin
               xfrac(i) = xsort(i) - yfrac(i)*gamma_cos
     &                       - zfrac(i)*beta_cos
            end do
         end if
      end if
c
c     use images to move coordinates into periodic cell
c
      if (use_bounds) then
         do i = 1, nsite
            xsort(i) = xfrac(i)
            ysort(i) = yfrac(i)
            zsort(i) = zfrac(i)
            do while (abs(xsort(i)) .gt. xcell2)
               xsort(i) = xsort(i) - sign(xcell,xsort(i))
            end do
            do while (abs(ysort(i)) .gt. ycell2)
               ysort(i) = ysort(i) - sign(ycell,ysort(i))
            end do
            do while (abs(zsort(i)) .gt. zcell2)
               zsort(i) = zsort(i) - sign(zcell,zsort(i))
            end do
         end do
      end if
c
c     generate the replica coordinates for the sort arrays
c
      if (use_replica) then
         k = nsite
         do j = 1, ncell
            xmove = icell(1,j) * xbox
            ymove = icell(2,j) * ybox
            zmove = icell(3,j) * zbox
            do i = 1, nsite
               k = k + 1
               xsort(k) = xfrac(i) + xmove
               ysort(k) = yfrac(i) + ymove
               zsort(k) = zfrac(i) + zmove
               do while (abs(xsort(k)) .gt. xcell2)
                  xsort(k) = xsort(k) - sign(xcell,xsort(k))
               end do
               do while (abs(ysort(k)) .gt. ycell2)
                  ysort(k) = ysort(k) - sign(ycell,ysort(k))
               end do
               do while (abs(zsort(k)) .gt. zcell2)
                  zsort(k) = zsort(k) - sign(zcell,zsort(k))
               end do
            end do
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (xfrac)
      deallocate (yfrac)
      deallocate (zfrac)
c
c     perform dynamic allocation of some global arrays
c
      nlight = (ncell+1) * nsite
      extent = 0
      if (allocated(rgx))  extent = size(rgx)
      if (extent .lt. nlight) then
         if (allocated(kbx))  deallocate (kbx)
         if (allocated(kby))  deallocate (kby)
         if (allocated(kbz))  deallocate (kbz)
         if (allocated(kex))  deallocate (kex)
         if (allocated(key))  deallocate (key)
         if (allocated(kez))  deallocate (kez)
         if (allocated(locx))  deallocate (locx)
         if (allocated(locy))  deallocate (locy)
         if (allocated(locz))  deallocate (locz)
         if (allocated(rgx))  deallocate (rgx)
         if (allocated(rgy))  deallocate (rgy)
         if (allocated(rgz))  deallocate (rgz)
         allocate (kbx(nsite))
         allocate (kby(nsite))
         allocate (kbz(nsite))
         allocate (kex(nsite))
         allocate (key(nsite))
         allocate (kez(nsite))
         allocate (locx(nlight))
         allocate (locy(nlight))
         allocate (locz(nlight))
         allocate (rgx(nlight))
         allocate (rgy(nlight))
         allocate (rgz(nlight))
      end if
c
c     sort the coordinate components into ascending order
c
      call sort2 (nlight,xsort,locx)
      call sort2 (nlight,ysort,locy)
      call sort2 (nlight,zsort,locz)
c
c     use of replicates requires secondary sorting along x-axis
c
      if (use_replica) then
         j = 1
         do i = 1, nlight-1
            if (xsort(i+1) .ne. xsort(i)) then
               call sort5 (i-j+1,locx(j),nsite)
               j = i + 1
            end if
         end do
         call sort5 (nlight-j+1,locx(j),nsite)
      end if
c
c     index the position of each atom in the sorted coordinates
c
      do i = 1, nlight
         rgx(locx(i)) = i
         rgy(locy(i)) = i
         rgz(locz(i)) = i
      end do
c
c     find the negative x-coordinate boundary for each atom
c
      if (unique) then
         do i = nlight, 1, -1
            k = locx(i)
            if (k .le. nsite) then
               kbx(k) = i
            end if
         end do
      else
         j = nlight
         box = 0.0d0
         do i = nlight, 1, -1
            k = locx(i)
            do while (xsort(i)-xsort(j)+box .le. xcut)
               if (j .eq. 1) then
                  if (use_bounds) then
                     j = nlight + 1
                     box = xcell
                  end if
               end if
               j = j - 1
               if (j .lt. 1)  goto 30
            end do
   30       continue
            j = j + 1
            if (j .gt. nlight) then
               j = 1
               box = 0.0d0
            end if
            kbx(k) = j
         end do
      end if
c
c     find the positive x-coordinate boundary for each atom
c
      j = 1
      box = 0.0d0
      do i = 1, nlight
         k = locx(i)
         if (k .le. nsite) then
            do while (xsort(j)-xsort(i)+box .lt. xcut)
               if (j .eq. nlight) then
                  if (use_bounds) then
                     j = 0
                     box = xcell
                  end if
               end if
               j = j + 1
               if (j .gt. nlight)  goto 40
            end do
   40       continue
            j = j - 1
            if (j .lt. 1) then
               j = nlight
               box = 0.0d0
            end if
            kex(k) = j
         end if
      end do
c
c     find the negative y-coordinate boundary for each atom
c
      j = nlight
      box = 0.0d0
      do i = nlight, 1, -1
         k = locy(i)
         if (k .le. nsite) then
            do while (ysort(i)-ysort(j)+box .le. ycut)
               if (j .eq. 1) then
                  if (use_bounds) then
                     j = nlight + 1
                     box = ycell
                  end if
               end if
               j = j - 1
               if (j .lt. 1)  goto 50
            end do
   50       continue
            j = j + 1
            if (j .gt. nlight) then
               j = 1
               box = 0.0d0
            end if
            kby(k) = j
         end if
      end do
c
c     find the positive y-coordinate boundary for each atom
c
      j = 1
      box = 0.0d0
      do i = 1, nlight
         k = locy(i)
         if (k .le. nsite) then
            do while (ysort(j)-ysort(i)+box .lt. ycut)
               if (j .eq. nlight) then
                  if (use_bounds) then
                     j = 0
                     box = ycell
                  end if
               end if
               j = j + 1
               if (j .gt. nlight)  goto 60
            end do
   60       continue
            j = j - 1
            if (j .lt. 1) then
               j = nlight
               box = 0.0d0
            end if
            key(k) = j
         end if
      end do
c
c     find the negative z-coordinate boundary for each atom
c
      j = nlight
      box = 0.0d0
      do i = nlight, 1, -1
         k = locz(i)
         if (k .le. nsite) then
            do while (zsort(i)-zsort(j)+box .le. zcut)
               if (j .eq. 1) then
                  if (use_bounds) then
                     j = nlight + 1
                     box = zcell
                  end if
               end if
               j = j - 1
               if (j .lt. 1)  goto 70
            end do
   70       continue
            j = j + 1
            if (j .gt. nlight) then
               j = 1
               box = 0.0d0
            end if
            kbz(k) = j
         end if
      end do
c
c     find the positive z-coordinate boundary for each atom
c
      j = 1
      box = 0.0d0
      do i = 1, nlight
         k = locz(i)
         if (k .le. nsite) then
            do while (zsort(j)-zsort(i)+box .lt. zcut)
               if (j .eq. nlight) then
                  if (use_bounds) then
                     j = 0
                     box = zcell
                  end if
               end if
               j = j + 1
               if (j .gt. nlight)  goto 80
            end do
   80       continue
            j = j - 1
            if (j .lt. 1) then
               j = nlight
               box = 0.0d0
            end if
            kez(k) = j
         end if
      end do
      return
      end
