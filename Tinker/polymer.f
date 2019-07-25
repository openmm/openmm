c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2001  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine polymer  --  check for an infinite polymer  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "polymer" tests for the presence of an infinite polymer
c     extending across periodic boundaries
c
c
      subroutine polymer
      use sizes
      use atoms
      use bndstr
      use bound
      use boxes
      use iounit
      use keys
      implicit none
      integer i,j,next
      integer ia,ib
      real*8 xr,yr,zr
      real*8 xab,yab,zab
      real*8 eps,delta
      real*8 xlimit
      real*8 ylimit
      real*8 zlimit
      real*8 maximage
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     set defaults of infinite polymer usage and cutoff distance
c
      use_polymer = .false.
      polycut = 5.5d0
c
c     get any keywords containing infinite polymer cutoff parameters
c
      do j = 1, nkey
         next = 1
         record = keyline(j)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:15) .eq. 'POLYMER-CUTOFF ') then
            string = record(next:240)
            read (string,*,err=10,end=10)  polycut
   10       continue
         end if
      end do
c
c     see if any bond connections require a minimum image
c
      if (use_bounds) then
         eps = 0.0001d0
         do i = 1, nbond
            ia = ibnd(1,i)
            ib = ibnd(2,i)
            xab = x(ia) - x(ib)
            yab = y(ia) - y(ib)
            zab = z(ia) - z(ib)
            xr = xab
            yr = yab
            zr = zab
            call image (xr,yr,zr)
            delta = abs(xr-xab) + abs(yr-yab) + abs(zr-zab)
            if (delta .gt. eps) then
               use_polymer = .true.
               goto 20
            end if
         end do
   20    continue
      end if
c
c     find the maximum sphere radius inscribed in periodic box
c
      if (use_polymer) then
         if (orthogonal) then
            xlimit = xbox2
            ylimit = ybox2
            zlimit = zbox2
         else if (monoclinic) then
            xlimit = xbox2 * beta_sin
            ylimit = ybox2
            zlimit = zbox2 * beta_sin
         else if (triclinic) then
            xlimit = xbox2 * beta_sin * gamma_sin
            ylimit = ybox2 * gamma_sin
            zlimit = zbox2 * beta_sin
         else if (octahedron) then
            xlimit = (sqrt(3.0d0)/4.0d0) * xbox
            ylimit = xlimit
            zlimit = xlimit
         end if
         maximage = min(xlimit,ylimit,zlimit)
c
c     check for too large or small an infinite polymer cutoff
c
         if (polycut .gt. maximage) then
            write (iout,30)
   30       format (/,' POLYMER  --  Image Conflicts for Infinite',
     &                 ' Polymer in Small Cell')
            call fatal
         else if (polycut .lt. 5.5d0) then
            write (iout,40)
   40       format (/,' POLYMER  --  Warning, Infinite Polymer',
     &                 ' Cutoff may be Too Small')
         end if
      end if
c
c     set square of cutoff distance for use with nonbonded terms
c
      polycut2 = polycut * polycut
      return
      end
