c
c
c     ###############################################################
c     ##  COPYRIGHT (C) 2006 by David Gohara & Jay William Ponder  ##
c     ##                    All Rights Reserved                    ##
c     ###############################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine nblist  --  maintain pairwise neighbor lists  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "nblist" constructs and maintains nonbonded pair neighbor lists
c     for vdw, electrostatic and polarization interactions
c
c
      subroutine nblist
      use limits
      use potent
      implicit none
c
c
c     update the vdw, CT and electrostatic neighbor lists
c
      if (use_vdw .and. use_vlist)  call vlist
      if (use_ct .and. use_ctlist)  call ctlist 
      if ((use_charge.or.use_solv) .and. use_clist)  call clist
      if ((use_mpole.or.use_polar.or.use_solv) .and. use_mlist)
     &      call mlist
      if (use_polar .and. use_ulist)  call ulist
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine vlist  --  get van der Waals neighbor lists  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "vlist" performs an update or a complete rebuild of the
c     van der Waals neighbor list
c
c
      subroutine vlist
      use sizes
      use atoms
      use bound
      use boxes
      use iounit
      use neigh
      use vdw
      implicit none
      integer i,j,k
      integer ii,iv
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 radius
      real*8 rdn,r2
      real*8, allocatable :: xred(:)
      real*8, allocatable :: yred(:)
      real*8, allocatable :: zred(:)
      logical, allocatable :: update(:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (xred(n))
      allocate (yred(n))
      allocate (zred(n))
      allocate (update(n))
c
c     apply reduction factors to find coordinates for each site
c
      do i = 1, nvdw
         ii = ivdw(i)
         iv = ired(ii)
         rdn = kred(ii)
         xred(i) = rdn*(x(ii)-x(iv)) + x(iv)
         yred(i) = rdn*(y(ii)-y(iv)) + y(iv)
         zred(i) = rdn*(z(ii)-z(iv)) + z(iv)
      end do
c
c     neighbor list cannot be used with the replicates method
c
      radius = sqrt(vbuf2)
      call replica (radius)
      if (use_replica) then
         write (iout,10)
   10    format (/,' VLIST  --  Pairwise Neighbor List cannot',
     &              ' be used with Replicas')
         call fatal
      end if
c
c     perform a complete list build instead of an update
c
      if (dovlst) then
         dovlst = .false.
         if (octahedron) then
            call vbuild (xred,yred,zred)
         else
            call vlight (xred,yred,zred)
         end if
         return
      end if
c
c     test sites for displacement exceeding half the buffer, and
c     rebuild the higher numbered neighbors of updated sites
c
!$OMP PARALLEL default(shared) private(i,j,k,xi,yi,zi,xr,yr,zr,r2)
!$OMP DO schedule(guided)
      do i = 1, nvdw
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         xr = xi - xvold(i)
         yr = yi - yvold(i)
         zr = zi - zvold(i)
         call imagen (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr
         update(i) = .false.
         if (r2 .ge. lbuf2) then
            update(i) = .true.
            xvold(i) = xi
            yvold(i) = yi
            zvold(i) = zi
            nvlst(i) = 0
            do k = i+1, nvdw
               xr = xi - xred(k)
               yr = yi - yred(k)
               zr = zi - zred(k)
               call imagen (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. vbuf2) then
                  nvlst(i) = nvlst(i) + 1
                  vlst(nvlst(i),i) = k
               end if
            end do
         end if
      end do
!$OMP END DO
c
c     adjust lists of lower numbered neighbors of updated sites
c
!$OMP DO schedule(guided)
      do i = 1, nvdw
         if (update(i)) then
            xi = xred(i)
            yi = yred(i)
            zi = zred(i)
            do k = 1, i-1
               if (.not. update(k)) then
                  xr = xi - xvold(k)
                  yr = yi - yvold(k)
                  zr = zi - zvold(k)
                  call imagen (xr,yr,zr)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. vbuf2) then
!$OMP CRITICAL
                     do j = 1, nvlst(k)
                        if (vlst(j,k) .eq. i)  goto 20
                     end do
                     nvlst(k) = nvlst(k) + 1
                     vlst(nvlst(k),k) = i
   20                continue
!$OMP END CRITICAL
                  else if (r2 .le. vbufx) then
!$OMP CRITICAL
                     do j = 1, nvlst(k)
                        if (vlst(j,k) .eq. i) then
                           vlst(j,k) = vlst(nvlst(k),k)
                           nvlst(k) = nvlst(k) - 1
                           goto 30
                        end if
                     end do
   30                continue
!$OMP END CRITICAL
                  end if
               end if
            end do
         end if
      end do
!$OMP END DO
c
c     check to see if any neighbor lists are too long
c
!$OMP DO schedule(guided)
      do i = 1, nvdw
         if (nvlst(i) .ge. maxvlst) then
            write (iout,40)
   40       format (/,' VLIST  --  Too many Neighbors;',
     &                 ' Increase MAXVLST')
            call fatal
         end if
      end do
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (xred)
      deallocate (yred)
      deallocate (zred)
      deallocate (update)
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine vbuild  --  build vdw list for all sites  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "vbuild" performs a complete rebuild of the van der Waals
c     pair neighbor list for all sites
c
c
      subroutine vbuild (xred,yred,zred)
      use sizes
      use bound
      use iounit
      use neigh
      use vdw
      implicit none
      integer i,k
      real*8 xi,yi,zi
      real*8 xr,yr,zr,r2
      real*8 xred(*)
      real*8 yred(*)
      real*8 zred(*)
c
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private(i,k,xi,yi,zi,xr,yr,zr,r2)
!$OMP DO schedule(guided)
c
c     store coordinates to reflect update of the site
c
      do i = 1, nvdw
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         xvold(i) = xi
         yvold(i) = yi
         zvold(i) = zi
c
c     generate all neighbors for the site being rebuilt
c
         nvlst(i) = 0
         do k = i+1, nvdw
            xr = xi - xred(k)
            yr = yi - yred(k)
            zr = zi - zred(k)
            call imagen (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. vbuf2) then
               nvlst(i) = nvlst(i) + 1
               vlst(nvlst(i),i) = k
            end if
         end do
c
c     check to see if the neighbor list is too long
c
         if (nvlst(i) .ge. maxvlst) then
            write (iout,10)
   10       format (/,' VBUILD  --  Too many Neighbors;',
     &                 ' Increase MAXVLST')
            call fatal
         end if
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
c     ##  subroutine vlight  --  build vdw pair list via lights  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "vlight" performs a complete rebuild of the van der Waals
c     pair neighbor list for all sites using the method of lights
c
c
      subroutine vlight (xred,yred,zred)
      use sizes
      use atoms
      use bound
      use cell
      use iounit
      use light
      use neigh
      use vdw
      implicit none
      integer i,j,k
      integer kgy,kgz
      integer start,stop
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r2,off
      real*8 xred(*)
      real*8 yred(*)
      real*8 zred(*)
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      logical unique,repeat
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (xsort(nvdw))
      allocate (ysort(nvdw))
      allocate (zsort(nvdw))
c
c     transfer interaction site coordinates to sorting arrays
c
      do i = 1, nvdw
         nvlst(i) = 0
         xvold(i) = xred(i)
         yvold(i) = yred(i)
         zvold(i) = zred(i)
         xsort(i) = xred(i)
         ysort(i) = yred(i)
         zsort(i) = zred(i)
      end do
c
c     use the method of lights to generate neighbors
c
      unique = .false.
      off = sqrt(vbuf2)
      call lights (off,nvdw,xsort,ysort,zsort,unique)
c
c     perform deallocation of some local arrays
c
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private(i,j,k,xi,yi,zi,
!$OMP& xr,yr,zr,r2,kgy,kgz,start,stop,repeat)
!$OMP DO schedule(guided)
c
c     loop over all atoms computing the neighbor lists
c
      do i = 1, nvdw
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         if (kbx(i) .le. kex(i)) then
            repeat = .false.
            start = kbx(i)
            stop = kex(i)
         else
            repeat = .true.
            start = 1
            stop = kex(i)
         end if
   10    continue
         do j = start, stop
            k = locx(j)
            if (k .le. i)  goto 20
            kgy = rgy(k)
            if (kby(i) .le. key(i)) then
               if (kgy.lt.kby(i) .or. kgy.gt.key(i))  goto 20
            else
               if (kgy.lt.kby(i) .and. kgy.gt.key(i))  goto 20
            end if
            kgz = rgz(k)
            if (kbz(i) .le. kez(i)) then
               if (kgz.lt.kbz(i) .or. kgz.gt.kez(i))  goto 20
            else
               if (kgz.lt.kbz(i) .and. kgz.gt.kez(i))  goto 20
            end if
            xr = xi - xred(k)
            yr = yi - yred(k)
            zr = zi - zred(k)
            call imagen (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. vbuf2) then
               nvlst(i) = nvlst(i) + 1
               vlst(nvlst(i),i) = k
            end if
   20       continue
         end do
         if (repeat) then
            repeat = .false.
            start = kbx(i)
            stop = nvdw
            goto 10
         end if
c
c     check to see if the neighbor list is too long
c
         if (nvlst(i) .ge. maxvlst) then
            write (iout,30)
   30       format (/,' VLIGHT  --  Too many Neighbors;',
     &                 ' Increase MAXVLST')
            call fatal
         end if
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
c     ##  subroutine clist  --  get partial charge neighbor lists  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "clist" performs an update or a complete rebuild of the
c     electrostatic neighbor lists for partial charges
c
c
      subroutine clist
      use sizes
      use atoms
      use bound
      use boxes
      use charge
      use iounit
      use neigh
      implicit none
      integer i,j,k
      integer ii,kk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 radius,r2
      logical, allocatable :: update(:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (update(n))
c
c     neighbor list cannot be used with the replicates method
c
      radius = sqrt(cbuf2)
      call replica (radius)
      if (use_replica) then
         write (iout,10)
   10    format (/,' CLIST  --  Pairwise Neighbor List cannot',
     &              ' be used with Replicas')
         call fatal
      end if
c
c     perform a complete list build instead of an update
c
      if (doclst) then
         doclst = .false.
         if (octahedron) then
            call cbuild
         else
            call clight
         end if
         return
      end if
c
c     test sites for displacement exceeding half the buffer, and
c     rebuild the higher numbered neighbors of updated sites
c
!$OMP PARALLEL default(shared) private(i,j,k,ii,kk,xi,yi,zi,xr,yr,zr,r2)
!$OMP DO schedule(guided)
      do i = 1, nion
         ii = kion(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         xr = xi - xcold(i)
         yr = yi - ycold(i)
         zr = zi - zcold(i)
         call imagen (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr
         update(i) = .false.
         if (r2 .ge. lbuf2) then
            update(i) = .true.
            xcold(i) = xi
            ycold(i) = yi
            zcold(i) = zi
            nelst(i) = 0
            do k = i+1, nion
               kk = kion(k)
               xr = xi - x(kk)
               yr = yi - y(kk)
               zr = zi - z(kk)
               call imagen (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. cbuf2) then
                  nelst(i) = nelst(i) + 1
                  elst(nelst(i),i) = k
               end if
            end do
         end if
      end do
!$OMP END DO
c
c     adjust lists of lower numbered neighbors of updated sites
c
!$OMP DO schedule(guided)
      do i = 1, nion
         if (update(i)) then
            ii = kion(i)
            xi = x(ii)
            yi = y(ii)
            zi = z(ii)
            do k = 1, i-1
               if (.not. update(k)) then
                  xr = xi - xcold(k)
                  yr = yi - ycold(k)
                  zr = zi - zcold(k)
                  call imagen (xr,yr,zr)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. cbuf2) then
!$OMP CRITICAL
                     do j = 1, nelst(k)
                        if (elst(j,k) .eq. i)  goto 20
                     end do
                     nelst(k) = nelst(k) + 1
                     elst(nelst(k),k) = i
   20                continue
!$OMP END CRITICAL
                  else if (r2 .le. cbufx) then
!$OMP CRITICAL
                     do j = 1, nelst(k)
                        if (elst(j,k) .eq. i) then
                           elst(j,k) = elst(nelst(k),k)
                           nelst(k) = nelst(k) - 1
                           goto 30
                        end if
                     end do
   30                continue
!$OMP END CRITICAL
                  end if
               end if
            end do
         end if
      end do
!$OMP END DO
c
c     check to see if any neighbor lists are too long
c
!$OMP DO schedule(guided)
      do i = 1, nion
         if (nelst(i) .ge. maxelst) then
            write (iout,40)
   40       format (/,' CLIST  --  Too many Neighbors;',
     &                 ' Increase MAXELST')
            call fatal
         end if
      end do
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (update)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine cbuild  --  build charge list for all sites  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "cbuild" performs a complete rebuild of the partial charge
c     electrostatic neighbor list for all sites
c
c
      subroutine cbuild
      use sizes
      use atoms
      use bound
      use charge
      use iounit
      use neigh
      implicit none
      integer i,k
      integer ii,kk
      real*8 xi,yi,zi
      real*8 xr,yr,zr,r2
c
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private(i,k,ii,kk,xi,yi,zi,xr,yr,zr,r2)
!$OMP DO schedule(guided)
c
c     store new coordinates to reflect update of the site
c
      do i = 1, nion
         ii = kion(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         xcold(i) = xi
         ycold(i) = yi
         zcold(i) = zi
c
c     generate all neighbors for the site being rebuilt
c
         nelst(i) = 0
         do k = i+1, nion
            kk = kion(k)
            xr = xi - x(kk)
            yr = yi - y(kk)
            zr = zi - z(kk)
            call imagen (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. cbuf2) then
               nelst(i) = nelst(i) + 1
               elst(nelst(i),i) = k
            end if
         end do
c
c     check to see if the neighbor list is too long
c
         if (nelst(i) .ge. maxelst) then
            write (iout,10)
   10       format (/,' CBUILD  --  Too many Neighbors;',
     &                 ' Increase MAXELST')
            call fatal
         end if
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
c     ##  subroutine clight  --  get partial charge list via lights  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "clight" performs a complete rebuild of the partial charge
c     pair neighbor list for all sites using the method of lights
c
c
      subroutine clight
      use sizes
      use atoms
      use bound
      use cell
      use charge
      use iounit
      use light
      use neigh
      implicit none
      integer i,j,k
      integer ii,kk
      integer kgy,kgz
      integer start,stop
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r2,off
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      logical unique,repeat
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (xsort(nion))
      allocate (ysort(nion))
      allocate (zsort(nion))
c
c     transfer interaction site coordinates to sorting arrays
c
      do i = 1, nion
         nelst(i) = 0
         ii = kion(i)
         xcold(i) = x(ii)
         ycold(i) = y(ii)
         zcold(i) = z(ii)
         xsort(i) = x(ii)
         ysort(i) = y(ii)
         zsort(i) = z(ii)
      end do
c
c     use the method of lights to generate neighbors
c
      unique = .false.
      off = sqrt(cbuf2)
      call lights (off,nion,xsort,ysort,zsort,unique)
c
c     perform deallocation of some local arrays
c
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private(i,j,k,ii,kk,xi,yi,zi,
!$OMP& xr,yr,zr,r2,kgy,kgz,start,stop,repeat)
!$OMP DO schedule(guided)
c
c     loop over all atoms computing the neighbor lists
c
      do i = 1, nion
         ii = kion(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         if (kbx(i) .le. kex(i)) then
            repeat = .false.
            start = kbx(i)
            stop = kex(i)
         else
            repeat = .true.
            start = 1
            stop = kex(i)
         end if
   10    continue
         do j = start, stop
            k = locx(j)
            if (k .le. i)  goto 20
            kk = kion(k)
            kgy = rgy(k)
            if (kby(i) .le. key(i)) then
               if (kgy.lt.kby(i) .or. kgy.gt.key(i))  goto 20
            else
               if (kgy.lt.kby(i) .and. kgy.gt.key(i))  goto 20
            end if
            kgz = rgz(k)
            if (kbz(i) .le. kez(i)) then
               if (kgz.lt.kbz(i) .or. kgz.gt.kez(i))  goto 20
            else
               if (kgz.lt.kbz(i) .and. kgz.gt.kez(i))  goto 20
            end if
            xr = xi - x(kk)
            yr = yi - y(kk)
            zr = zi - z(kk)
            call imagen (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. cbuf2) then
               nelst(i) = nelst(i) + 1
               elst(nelst(i),i) = k
            end if
   20       continue
         end do
         if (repeat) then
            repeat = .false.
            start = kbx(i)
            stop = nion
            goto 10
         end if
c
c     check to see if the neighbor list is too long
c
         if (nelst(i) .ge. maxelst) then
            write (iout,30)
   30       format (/,' CLIGHT  --  Too many Neighbors;',
     &                 ' Increase MAXELST')
            call fatal
         end if
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
c     ##  subroutine mlist  --  get atomic multipole neighbor lists  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "mlist" performs an update or a complete rebuild of the
c     electrostatic neighbor lists for atomic multipoles
c
c
      subroutine mlist
      use sizes
      use atoms
      use bound
      use boxes
      use iounit
      use mpole
      use neigh
      implicit none
      integer i,j,k
      integer ii,kk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 radius,r2
      logical, allocatable :: update(:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (update(n))
c
c     neighbor list cannot be used with the replicates method
c
      radius = sqrt(mbuf2)
      call replica (radius)
      if (use_replica) then
         write (iout,10)
   10    format (/,' MLIST  --  Pairwise Neighbor List cannot',
     &              ' be used with Replicas')
         call fatal
      end if
c
c     perform a complete list build instead of an update
c
      if (domlst) then
         domlst = .false.
         if (octahedron) then
            call mbuild
         else
            call mlight
         end if
         return
      end if
c
c     test sites for displacement exceeding half the buffer, and
c     rebuild the higher numbered neighbors of updated sites
c
!$OMP PARALLEL default(shared) private(i,j,k,ii,kk,xi,yi,zi,xr,yr,zr,r2)
!$OMP DO schedule(guided)
      do i = 1, npole
         ii = ipole(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         xr = xi - xmold(i)
         yr = yi - ymold(i)
         zr = zi - zmold(i)
         call imagen (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr
         update(i) = .false.
         if (r2 .ge. lbuf2) then
            update(i) = .true.
            xmold(i) = xi
            ymold(i) = yi
            zmold(i) = zi
            nelst(i) = 0
            do k = i+1, npole
               kk = ipole(k)
               xr = xi - x(kk)
               yr = yi - y(kk)
               zr = zi - z(kk)
               call imagen (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. mbuf2) then
                  nelst(i) = nelst(i) + 1
                  elst(nelst(i),i) = k
               end if
            end do
         end if
      end do
!$OMP END DO
c
c     adjust lists of lower numbered neighbors of updated sites
c
!$OMP DO schedule (guided)
      do i = 1, npole
         if (update(i)) then
            ii = ipole(i)
            xi = x(ii)
            yi = y(ii)
            zi = z(ii)
            do k = 1, i-1
               if (.not. update(k)) then
                  xr = xi - xmold(k)
                  yr = yi - ymold(k)
                  zr = zi - zmold(k)
                  call imagen (xr,yr,zr)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. mbuf2) then
!$OMP CRITICAL
                     do j = 1, nelst(k)
                        if (elst(j,k) .eq. i)  goto 20
                     end do
                     nelst(k) = nelst(k) + 1
                     elst(nelst(k),k) = i
   20                continue
!$OMP END CRITICAL
                  else if (r2 .le. mbufx) then
!$OMP CRITICAL
                     do j = 1, nelst(k)
                        if (elst(j,k) .eq. i) then
                           elst(j,k) = elst(nelst(k),k)
                           nelst(k) = nelst(k) - 1
                           goto 30
                        end if
                     end do
   30                continue
!$OMP END CRITICAL
                  end if
               end if
            end do
         end if
      end do
!$OMP END DO
c
c     check to see if any neighbor lists are too long
c
!$OMP DO schedule(guided)
      do i = 1, npole
         if (nelst(i) .ge. maxelst) then
            write (iout,40)
   40       format (/,' MLIST  --  Too many Neighbors;',
     &                 ' Increase MAXELST')
            call fatal
         end if
      end do
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (update)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine mbuild  --  build mpole list for all sites  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "mbuild" performs a complete rebuild of the atomic multipole
c     electrostatic neighbor list for all sites
c
c
      subroutine mbuild
      use sizes
      use atoms
      use bound
      use iounit
      use mpole
      use neigh
      implicit none
      integer i,k
      integer ii,kk
      real*8 xi,yi,zi
      real*8 xr,yr,zr,r2
c
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private(i,k,ii,kk,xi,yi,zi,xr,yr,zr,r2)
!$OMP DO schedule(guided)
c
c     store new coordinates to reflect update of the site
c
      do i = 1, npole
         ii = ipole(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         xmold(i) = xi
         ymold(i) = yi
         zmold(i) = zi
c
c     generate all neighbors for the site being rebuilt
c
         nelst(i) = 0
         do k = i+1, npole
            kk = ipole(k)
            xr = xi - x(kk)
            yr = yi - y(kk)
            zr = zi - z(kk)
            call imagen (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. mbuf2) then
               nelst(i) = nelst(i) + 1
               elst(nelst(i),i) = k
            end if
         end do
c
c     check to see if the neighbor list is too long
c
         if (nelst(i) .ge. maxelst) then
            write (iout,10)
   10       format (/,' MBUILD  --  Too many Neighbors;',
     &                 ' Increase MAXELST')
            call fatal
         end if
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
c     ##  subroutine mlight  --  get multipole pair list via lights  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "mlight" performs a complete rebuild of the atomic multipole
c     pair neighbor list for all sites using the method of lights
c
c
      subroutine mlight
      use sizes
      use atoms
      use bound
      use cell
      use iounit
      use light
      use mpole
      use neigh
      implicit none
      integer i,j,k
      integer ii,kk
      integer kgy,kgz
      integer start,stop
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r2,off
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      logical unique,repeat
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (xsort(npole))
      allocate (ysort(npole))
      allocate (zsort(npole))
c
c     transfer interaction site coordinates to sorting arrays
c
      do i = 1, npole
         nelst(i) = 0
         ii = ipole(i)
         xmold(i) = x(ii)
         ymold(i) = y(ii)
         zmold(i) = z(ii)
         xsort(i) = x(ii)
         ysort(i) = y(ii)
         zsort(i) = z(ii)
      end do
c
c     use the method of lights to generate neighbors
c
      unique = .false.
      off = sqrt(mbuf2)
      call lights (off,npole,xsort,ysort,zsort,unique)
c
c     perform deallocation of some local arrays
c
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private(i,j,k,ii,kk,xi,yi,zi,
!$OMP& xr,yr,zr,r2,kgy,kgz,start,stop,repeat)
!$OMP DO schedule(guided)
c
c     loop over all atoms computing the neighbor lists
c
      do i = 1, npole
         ii = ipole(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         if (kbx(i) .le. kex(i)) then
            repeat = .false.
            start = kbx(i)
            stop = kex(i)
         else
            repeat = .true.
            start = 1
            stop = kex(i)
         end if
   10    continue
         do j = start, stop
            k = locx(j)
            if (k .le. i)  goto 20
            kk = ipole(k)
            kgy = rgy(k)
            if (kby(i) .le. key(i)) then
               if (kgy.lt.kby(i) .or. kgy.gt.key(i))  goto 20
            else
               if (kgy.lt.kby(i) .and. kgy.gt.key(i))  goto 20
            end if
            kgz = rgz(k)
            if (kbz(i) .le. kez(i)) then
               if (kgz.lt.kbz(i) .or. kgz.gt.kez(i))  goto 20
            else
               if (kgz.lt.kbz(i) .and. kgz.gt.kez(i))  goto 20
            end if
            xr = xi - x(kk)
            yr = yi - y(kk)
            zr = zi - z(kk)
            call imagen (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. mbuf2) then
               nelst(i) = nelst(i) + 1
               elst(nelst(i),i) = k
            end if
   20       continue
         end do
         if (repeat) then
            repeat = .false.
            start = kbx(i)
            stop = npole
            goto 10
         end if
c
c     check to see if the neighbor list is too long
c
         if (nelst(i) .ge. maxelst) then
            write (iout,30)
   30       format (/,' MLIGHT  --  Too many Neighbors;',
     &                 ' Increase MAXELST')
            call fatal
         end if
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
c     ##  subroutine ulist  --  get preconditioner neighbor lists  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "ulist" performs an update or a complete rebuild of the
c     neighbor lists for the polarization preconditioner
c
c
      subroutine ulist
      use sizes
      use atoms
      use bound
      use boxes
      use iounit
      use mpole
      use neigh
      implicit none
      integer i,j,k
      integer ii,kk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 radius,r2
      logical, allocatable :: update(:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (update(n))
c
c     neighbor list cannot be used with the replicates method
c
      radius = sqrt(ubuf2)
      call replica (radius)
      if (use_replica) then
         write (iout,10)
   10    format (/,' ULIST  --  Pairwise Neighbor List cannot',
     &              ' be used with Replicas')
         call fatal
      end if
c
c     perform a complete list build instead of an update
c
      if (doulst) then
         doulst = .false.
         if (octahedron) then
            call ubuild
         else
            call ulight
         end if
         return
      end if
c
c     test sites for displacement exceeding half the buffer, and
c     rebuild the higher numbered neighbors of updated sites
c
!$OMP PARALLEL default(shared) private(i,j,k,ii,kk,xi,yi,zi,xr,yr,zr,r2)
!$OMP DO schedule(guided)
      do i = 1, npole
         ii = ipole(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         xr = xi - xuold(i)
         yr = yi - yuold(i)
         zr = zi - zuold(i)
         call imagen (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr
         update(i) = .false.
         if (r2 .ge. pbuf2) then
            update(i) = .true.
            xuold(i) = xi
            yuold(i) = yi
            zuold(i) = zi
            nulst(i) = 0
            do k = i+1, npole
               kk = ipole(k)
               xr = xi - x(kk)
               yr = yi - y(kk)
               zr = zi - z(kk)
               call imagen (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. ubuf2) then
                  nulst(i) = nulst(i) + 1
                  ulst(nulst(i),i) = k
               end if
            end do
         end if
      end do
!$OMP END DO
c
c     adjust lists of lower numbered neighbors of updated sites
c
!$OMP DO schedule(guided)
      do i = 1, npole
         if (update(i)) then
            ii = ipole(i)
            xi = x(ii)
            yi = y(ii)
            zi = z(ii)
            do k = 1, i-1
               if (.not. update(k)) then
                  xr = xi - xuold(k)
                  yr = yi - yuold(k)
                  zr = zi - zuold(k)
                  call imagen (xr,yr,zr)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. ubuf2) then
!$OMP CRITICAL
                     do j = 1, nulst(k)
                        if (ulst(j,k) .eq. i)  goto 20
                     end do
                     nulst(k) = nulst(k) + 1
                     ulst(nulst(k),k) = i
   20                continue
!$OMP END CRITICAL
                  else if (r2 .le. ubufx) then
!$OMP CRITICAL
                     do j = 1, nulst(k)
                        if (ulst(j,k) .eq. i) then
                           ulst(j,k) = ulst(nulst(k),k)
                           nulst(k) = nulst(k) - 1
                           goto 30
                        end if
                     end do
   30                continue
!$OMP END CRITICAL
                  end if
               end if
            end do
         end if
      end do
!$OMP END DO
c
c     check to see if any neighbor lists are too long
c
!$OMP DO schedule(guided)
      do i = 1, npole
         if (nulst(i) .ge. maxulst) then
            write (iout,40)
   40       format (/,' ULIST  --  Too many Neighbors;',
     &                 ' Increase MAXULST')
            call fatal
         end if
      end do
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (update)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ubuild  --  preconditioner list for all sites  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ubuild" performs a complete rebuild of the polarization
c     preconditioner neighbor list for all sites
c
c
      subroutine ubuild
      use sizes
      use atoms
      use bound
      use iounit
      use mpole
      use neigh
      implicit none
      integer i,k
      integer ii,kk
      real*8 xi,yi,zi
      real*8 xr,yr,zr,r2
c
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private(i,k,ii,kk,xi,yi,zi,xr,yr,zr,r2)
!$OMP DO schedule(guided)
c
c     store new coordinates to reflect update of the site
c
      do i = 1, npole
         ii = ipole(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         xuold(i) = xi
         yuold(i) = yi
         zuold(i) = zi
c
c     generate all neighbors for the site being rebuilt
c
         nulst(i) = 0
         do k = i+1, npole
            kk = ipole(k)
            xr = xi - x(kk)
            yr = yi - y(kk)
            zr = zi - z(kk)
            call imagen (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. ubuf2) then
               nulst(i) = nulst(i) + 1
               ulst(nulst(i),i) = k
            end if
         end do
c
c     check to see if the neighbor list is too long
c
         if (nulst(i) .ge. maxulst) then
            write (iout,10)
   10       format (/,' UBUILD  --  Too many Neighbors;',
     &                 ' Increase MAXULST')
            call fatal
         end if
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
c     ##  subroutine ulight  --  get preconditioner list via lights  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ulight" performs a complete rebuild of the polarization
c     preconditioner pair neighbor list for all sites using the
c     method of lights
c
c
      subroutine ulight
      use sizes
      use atoms
      use bound
      use cell
      use iounit
      use light
      use mpole
      use neigh
      implicit none
      integer i,j,k
      integer ii,kk
      integer kgy,kgz
      integer start,stop
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r2,off
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      logical unique,repeat
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (xsort(npole))
      allocate (ysort(npole))
      allocate (zsort(npole))
c
c     transfer interaction site coordinates to sorting arrays
c
      do i = 1, npole
         nulst(i) = 0
         ii = ipole(i)
         xuold(i) = x(ii)
         yuold(i) = y(ii)
         zuold(i) = z(ii)
         xsort(i) = x(ii)
         ysort(i) = y(ii)
         zsort(i) = z(ii)
      end do
c
c     use the method of lights to generate neighbors
c
      unique = .false.
      off = sqrt(ubuf2)
      call lights (off,npole,xsort,ysort,zsort,unique)
c
c     perform deallocation of some local arrays
c
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private(i,j,k,ii,kk,xi,yi,zi,
!$OMP& xr,yr,zr,r2,kgy,kgz,start,stop,repeat)
!$OMP DO schedule(guided)
c
c     loop over all atoms computing the neighbor lists
c
      do i = 1, npole
         ii = ipole(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         if (kbx(i) .le. kex(i)) then
            repeat = .false.
            start = kbx(i)
            stop = kex(i)
         else
            repeat = .true.
            start = 1
            stop = kex(i)
         end if
   10    continue
         do j = start, stop
            k = locx(j)
            if (k .le. i)  goto 20
            kk = ipole(k)
            kgy = rgy(k)
            if (kby(i) .le. key(i)) then
               if (kgy.lt.kby(i) .or. kgy.gt.key(i))  goto 20
            else
               if (kgy.lt.kby(i) .and. kgy.gt.key(i))  goto 20
            end if
            kgz = rgz(k)
            if (kbz(i) .le. kez(i)) then
               if (kgz.lt.kbz(i) .or. kgz.gt.kez(i))  goto 20
            else
               if (kgz.lt.kbz(i) .and. kgz.gt.kez(i))  goto 20
            end if
            xr = xi - x(kk)
            yr = yi - y(kk)
            zr = zi - z(kk)
            call imagen (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. ubuf2) then
               nulst(i) = nulst(i) + 1
               ulst(nulst(i),i) = k
            end if
   20       continue
         end do
         if (repeat) then
            repeat = .false.
            start = kbx(i)
            stop = npole
            goto 10
         end if
c
c     check to see if the neighbor list is too long
c
         if (nulst(i) .ge. maxulst) then
            write (iout,30)
   30       format (/,' ULIGHT  --  Too many Neighbors;',
     &                 ' Increase MAXULST')
            call fatal
         end if
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end

c     ###############################################################
c     ##                                                           ##
c     ##  subroutine ctlist  --  get charge transfer neighbor lists##
c     ##                                                           ##
c     ###############################################################
c
c
c     "ctlist" performs an update or a complete rebuild of the
c     charge transfer neighbor lists 
c
c
      subroutine ctlist
      use sizes
      use atoms
      use bound
      use boxes
      use ctran 
      use iounit
      use neigh
      implicit none
      integer i,j,k
      integer ii,kk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 radius,r2
      logical, allocatable :: update(:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (update(n))
c
c     neighbor list cannot be used with the replicates method
c
      radius = sqrt(ctbuf2)
      call replica (radius)
      if (use_replica) then
         write (iout,10)
   10    format (/,' CTLIST  --  Pairwise Neighbor List cannot',
     &              ' be used with Replicas')
         call fatal
      end if
c
c     perform a complete list build instead of an update
c
      if (doctlst) then
         doctlst = .false.
         if (octahedron) then
            call ctbuild
         else
            call ctlight
         end if
         return
      end if
c
c     test sites for displacement exceeding half the buffer, and
c     rebuild the higher numbered neighbors of updated sites
c
!$OMP PARALLEL default(shared) private(i,j,k,ii,kk,xi,yi,zi,xr,yr,zr,r2)
!$OMP DO schedule(guided)
      do i = 1, nct 
         ii = ict(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         xr = xi - xctold(i)
         yr = yi - yctold(i)
         zr = zi - zctold(i)
         call imagen (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr
         update(i) = .false.
         if (r2 .ge. lbuf2) then
            update(i) = .true.
            xctold(i) = xi
            yctold(i) = yi
            zctold(i) = zi
            nctlst(i) = 0
            do k = i+1, nct
               kk = ict(k)
               xr = xi - x(kk)
               yr = yi - y(kk)
               zr = zi - z(kk)
               call imagen (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. ctbuf2) then
                  nctlst(i) = nctlst(i) + 1
                  ctlst(nctlst(i),i) = k
               end if
            end do
         end if
      end do
!$OMP END DO
c
c     adjust lists of lower numbered neighbors of updated sites
c
!$OMP DO schedule(guided)
      do i = 1, nct
         if (update(i)) then
            ii = ict(i)
            xi = x(ii)
            yi = y(ii)
            zi = z(ii)
            do k = 1, i-1
               if (.not. update(k)) then
                  xr = xi - xctold(k)
                  yr = yi - yctold(k)
                  zr = zi - zctold(k)
                  call imagen (xr,yr,zr)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. ctbuf2) then
!$OMP CRITICAL
                     do j = 1, nctlst(k)
                        if (ctlst(j,k) .eq. i)  goto 20
                     end do
                     nctlst(k) = nctlst(k) + 1
                     ctlst(nctlst(k),k) = i
   20                continue
!$OMP END CRITICAL
                  else if (r2 .le. ctbufx) then
!$OMP CRITICAL
                     do j = 1, nctlst(k)
                        if (ctlst(j,k) .eq. i) then
                           ctlst(j,k) = ctlst(nctlst(k),k)
                           nctlst(k) = nctlst(k) - 1
                           goto 30
                        end if
                     end do
   30                continue
!$OMP END CRITICAL
                  end if
               end if
            end do
         end if
      end do
!$OMP END DO
c
c     check to see if any neighbor lists are too long
c
!$OMP DO schedule(guided)
      do i = 1, nct
         if (nctlst(i) .ge. maxctlst) then

            write (iout,40)
   40       format (/,' CTLIST  --  Too many Neighbors;',
     &                 ' Increase MAXCTLST')
            call fatal
         end if
      end do
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (update)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine cbuild  --  build charge list for all sites  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "cbuild" performs a complete rebuild of the partial charge
c     electrostatic neighbor list for all sites
c
c
      subroutine ctbuild
      use sizes
      use atoms
      use bound
      use ctran
      use iounit
      use neigh
      implicit none
      integer i,k
      integer ii,kk
      real*8 xi,yi,zi
      real*8 xr,yr,zr,r2
c
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private(i,k,ii,kk,xi,yi,zi,xr,yr,zr,r2)
!$OMP DO schedule(guided)
c
c     store new coordinates to reflect update of the site
c
      do i = 1, nct
         ii = ict(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         xctold(i) = xi
         yctold(i) = yi
         zctold(i) = zi
c
c     generate all neighbors for the site being rebuilt
c
         nctlst(i) = 0
         do k = i+1, nct
            kk = ict(k)
            xr = xi - x(kk)
            yr = yi - y(kk)
            zr = zi - z(kk)
            call imagen (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. ctbuf2) then
               nctlst(i) = nctlst(i) + 1
               ctlst(nctlst(i),i) = k
            end if
         end do
c
c     check to see if the neighbor list is too long
c
         if (nctlst(i) .ge. maxctlst) then
            write (iout,10)
   10       format (/,' CTBUILD  --  Too many Neighbors;',
     &                 ' Increase MAXCTLST')
            call fatal
         end if
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
c     ##  subroutine ctlight  --  get CT list via lights             ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ctlight" performs a complete rebuild of the CT 
c     pair neighbor list for all sites using the method of lights
c
c
      subroutine ctlight
      use sizes
      use atoms
      use bound
      use cell
      use ctran
      use iounit
      use light
      use neigh
      implicit none
      integer i,j,k
      integer ii,kk
      integer kgy,kgz
      integer start,stop
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r2,off
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      logical unique,repeat
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (xsort(nct))
      allocate (ysort(nct))
      allocate (zsort(nct))
c
c     transfer interaction site coordinates to sorting arrays
c
      do i = 1, nct
         nctlst(i) = 0
         ii = ict(i)
         xctold(i) = x(ii)
         yctold(i) = y(ii)
         zctold(i) = z(ii)
         xsort(i) = x(ii)
         ysort(i) = y(ii)
         zsort(i) = z(ii)
      end do
c
c     use the method of lights to generate neighbors
c
      unique = .false.
      off = sqrt(ctbuf2)
      call lights (off,nct,xsort,ysort,zsort,unique)
c
c     perform deallocation of some local arrays
c
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private(i,j,k,ii,kk,xi,yi,zi,
!$OMP& xr,yr,zr,r2,kgy,kgz,start,stop,repeat)
!$OMP DO schedule(guided)
c
c     loop over all atoms computing the neighbor lists
c
      do i = 1, nct
         ii = ict(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         if (kbx(i) .le. kex(i)) then
            repeat = .false.
            start = kbx(i)
            stop = kex(i)
         else
            repeat = .true.
            start = 1
            stop = kex(i)
         end if
   10    continue
         do j = start, stop
            k = locx(j)
            if (k .le. i)  goto 20
            kk = ict(k)
            kgy = rgy(k)
            if (kby(i) .le. key(i)) then
               if (kgy.lt.kby(i) .or. kgy.gt.key(i))  goto 20
            else
               if (kgy.lt.kby(i) .and. kgy.gt.key(i))  goto 20
            end if
            kgz = rgz(k)
            if (kbz(i) .le. kez(i)) then
               if (kgz.lt.kbz(i) .or. kgz.gt.kez(i))  goto 20
            else
               if (kgz.lt.kbz(i) .and. kgz.gt.kez(i))  goto 20
            end if
            xr = xi - x(kk)
            yr = yi - y(kk)
            zr = zi - z(kk)
            call imagen (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. ctbuf2) then
               nctlst(i) = nctlst(i) + 1
               ctlst(nctlst(i),i) = k
            end if
   20       continue
         end do
         if (repeat) then
            repeat = .false.
            start = kbx(i)
            stop = nct
            goto 10
         end if
c
c     check to see if the neighbor list is too long
c
         if (nctlst(i) .ge. maxctlst) then
            write (iout,30)
   30       format (/,' CTLIGHT  --  Too many Neighbors;',
     &                 ' Increase MAXCTLST')
            call fatal
         end if
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end
