c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine kgeom  --  restraint term parameter assignment  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "kgeom" asisgns parameters for geometric restraint terms
c     to be included in the potential energy calculation
c
c
      subroutine kgeom
      use sizes
      use atomid
      use atoms
      use bound
      use couple
      use group
      use iounit
      use keys
      use molcul
      use potent
      use restrn
      implicit none
      integer i,j,k,next
      integer ia,ib,ic,id
      real*8 p1,p2,p3,p4,p5
      real*8 d1,d2,d3
      real*8 a1,a2,a3
      real*8 t1,t2,t3
      real*8 g1,g2,g3
      real*8 xr,yr,zr
      real*8 xcm,ycm,zcm
      real*8 geometry,weigh
      real*8 xad,yad,zad
      real*8 xbd,ybd,zbd
      real*8 xcd,ycd,zcd
      real*8 c1,c2,c3
      real*8 vol,ratio
      logical exist,keep
      logical intermol
      logical first
      character*1 letter
      character*20 keyword
      character*240 record
      character*240 string
      save first
      data first  / .true. /
c
c
c     set the default values for the restraint variables
c
      npfix = 0
      ndfix = 0
      nafix = 0
      ntfix = 0
      ngfix = 0
      nchir = 0
      depth = 0.0d0
      width = 0.0d0
      rwall = 0.0d0
      use_basin = .false.
      use_wall = .false.
c
c     perform dynamic allocation of some global arrays
c
      if (first) then
         first = .false.
         if (.not. allocated(ipfix))  allocate (ipfix(maxfix))
         if (.not. allocated(kpfix))  allocate (kpfix(3,maxfix))
         if (.not. allocated(idfix))  allocate (idfix(2,maxfix))
         if (.not. allocated(iafix))  allocate (iafix(3,maxfix))
         if (.not. allocated(itfix))  allocate (itfix(4,maxfix))
         if (.not. allocated(igfix))  allocate (igfix(2,maxfix))
         if (.not. allocated(ichir))  allocate (ichir(4,maxfix))
         if (.not. allocated(xpfix))  allocate (xpfix(maxfix))
         if (.not. allocated(ypfix))  allocate (ypfix(maxfix))
         if (.not. allocated(zpfix))  allocate (zpfix(maxfix))
         if (.not. allocated(pfix))  allocate (pfix(2,maxfix))
         if (.not. allocated(dfix))  allocate (dfix(3,maxfix))
         if (.not. allocated(afix))  allocate (afix(3,maxfix))
         if (.not. allocated(tfix))  allocate (tfix(3,maxfix))
         if (.not. allocated(gfix))  allocate (gfix(3,maxfix))
         if (.not. allocated(chir))  allocate (chir(3,maxfix))
      end if
c
c     search the keywords for restraint parameters
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
c
c     get atom restrained to a specified position range
c
         if (keyword(1:18) .eq. 'RESTRAIN-POSITION ') then
            ia = 0
            ib = 0
            p1 = 0.0d0
            p2 = 0.0d0
            p3 = 0.0d0
            p4 = 0.0d0
            p5 = 0.0d0
            next = 1
            call getword (string,letter,next)
            if (letter .eq. ' ') then
               call getnumb (string,ia,next)
               if (ia.ge.1 .and. ia.le.n) then
                  p1 = x(ia)
                  p2 = y(ia)
                  p3 = z(ia)
                  string = string(next:240)
                  read (string,*,err=10,end=10)  p1,p2,p3,p4,p5
   10             continue
                  if (p4 .eq. 0.0d0)  p4 = 100.0d0
                  npfix = npfix + 1
                  ipfix(npfix) = ia
                  kpfix(1,npfix) = 1
                  kpfix(2,npfix) = 1
                  kpfix(3,npfix) = 1
                  xpfix(npfix) = p1
                  ypfix(npfix) = p2
                  zpfix(npfix) = p3
                  pfix(1,npfix) = p4
                  pfix(2,npfix) = p5
               else if (ia.ge.-n .and. ia.le.-1) then
                  ia = abs(ia)
                  call getnumb (string,ib,next)
                  ib = min(abs(ib),n)
                  string = string(next:240)
                  read (string,*,err=20,end=20)  p1,p2
   20             continue
                  if (p1 .eq. 0.0d0)  p1 = 100.0d0
                  do j = ia, ib
                     npfix = npfix + 1
                     ipfix(npfix) = j
                     kpfix(1,npfix) = 1
                     kpfix(2,npfix) = 1
                     kpfix(3,npfix) = 1
                     xpfix(npfix) = x(j)
                     ypfix(npfix) = y(j)
                     zpfix(npfix) = z(j)
                     pfix(1,npfix) = p1
                     pfix(2,npfix) = p2
                  end do
               end if
            else
               call upcase (letter)
               read (string,*,err=30,end=30)  ia
               string = string(next:240)
               read (string,*,err=30,end=30)  p1,p2,p3
   30          continue
               if (p2 .eq. 0.0d0)  p2 = 100.0d0
               npfix = npfix + 1
               ipfix(npfix) = ia
               kpfix(1,npfix) = 0
               kpfix(2,npfix) = 0
               kpfix(3,npfix) = 0
               if (letter .eq. 'X') then
                  kpfix(1,npfix) = 1
                  xpfix(npfix) = p1
               else if (letter .eq. 'Y') then
                  kpfix(2,npfix) = 1
                  ypfix(npfix) = p1
               else if (letter .eq. 'Z') then
                  kpfix(3,npfix) = 1
                  zpfix(npfix) = p1
               end if
               pfix(1,npfix) = p2
               pfix(2,npfix) = p3
            end if
            if (npfix .gt. maxfix) then
               write (iout,40)
   40          format (/,' KGEOM  --  Too many Position Restraints;',
     &                    ' Increase MAXFIX')
               call fatal
            end if
c
c     get atoms restrained to a specified distance range
c
         else if (keyword(1:18) .eq. 'RESTRAIN-DISTANCE ') then
            ia = 0
            ib = 0
            d1 = 100.0d0
            d2 = 0.0d0
            d3 = 0.0d0
            exist = .false.
            read (string,*,err=50,end=50)  ia,ib,d1,d2
            exist = .true.
   50       continue
            read (string,*,err=60,end=60)  ia,ib,d1,d2,d3
   60       continue
            if (.not. exist) then
               xr = x(ia) - x(ib)
               yr = y(ia) - y(ib)
               zr = z(ia) - z(ib)
               intermol = (molcule(ia) .ne. molcule(ib))
               if (use_bounds .and. intermol)  call image (xr,yr,zr)
               d2 = sqrt(xr*xr + yr*yr + zr*zr)
            end if
            if (d3 .eq. 0.0d0)  d3 = d2
            ndfix = ndfix + 1
            idfix(1,ndfix) = ia
            idfix(2,ndfix) = ib
            dfix(1,ndfix) = d1
            dfix(2,ndfix) = d2
            dfix(3,ndfix) = d3
            if (ndfix .gt. maxfix) then
               write (iout,70)
   70          format (/,' KGEOM  --  Too many Distance Restraints;',
     &                    ' Increase MAXFIX')
               call fatal
            end if
c
c     get atoms restrained to a specified angle range
c
         else if (keyword(1:15) .eq. 'RESTRAIN-ANGLE ') then
            ia = 0
            ib = 0
            ic = 0
            a1 = 10.0d0
            a2 = 0.0d0
            a3 = 0.0d0
            exist = .false.
            read (string,*,err=80,end=80)  ia,ib,ic,a1,a2
            exist = .true.
   80       continue
            read (string,*,err=90,end=90)  ia,ib,ic,a1,a2,a3
   90       continue
            if (.not. exist)  a2 = geometry (ia,ib,ic,0)
            if (a3 .eq. 0.0d0)  a3 = a2
            nafix = nafix + 1
            iafix(1,nafix) = ia
            iafix(2,nafix) = ib
            iafix(3,nafix) = ic
            afix(1,nafix) = a1
            afix(2,nafix) = a2
            afix(3,nafix) = a3
            if (nafix .gt. maxfix) then
               write (iout,100)
  100          format (/,' KGEOM  --  Too many Angle Restraints;',
     &                    ' Increase MAXFIX')
               call fatal
            end if
c
c     get atoms restrained to a specified torsion range
c
         else if (keyword(1:17).eq.'RESTRAIN-TORSION ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            t1 = 1.0d0
            t2 = 0.0d0
            t3 = 0.0d0
            exist = .false.
            read (string,*,err=110,end=110)  ia,ib,ic,id,t1,t2
            exist = .true.
  110       continue
            read (string,*,err=120,end=120)  ia,ib,ic,id,t1,t2,t3
            exist = .true.
  120       continue
            if (.not. exist)  t2 = geometry (ia,ib,ic,id)
            if (t3 .eq. 0.0d0)  t3 = t2
            do while (t2 .gt. 180.0d0)
               t2 = t2 - 360.0d0
            end do
            do while (t2 .lt. -180.0d0)
               t2 = t2 + 360.0d0
            end do
            do while (t3 .gt. 180.0d0)
               t3 = t3 - 360.0d0
            end do
            do while (t3 .lt. -180.0d0)
               t3 = t3 + 360.0d0
            end do
            ntfix = ntfix + 1
            itfix(1,ntfix) = ia
            itfix(2,ntfix) = ib
            itfix(3,ntfix) = ic
            itfix(4,ntfix) = id
            tfix(1,ntfix) = t1
            tfix(2,ntfix) = t2
            tfix(3,ntfix) = t3
            if (ntfix .gt. maxfix) then
               write (iout,130)
  130          format (/,' KGEOM  --  Too many Torsion Restraints;',
     &                    ' Increase MAXFIX')
               call fatal
            end if
c
c     get groups restrained to a specified distance range
c
         else if (keyword(1:16) .eq. 'RESTRAIN-GROUPS ') then
            ia = 0
            ib = 0
            g1 = 100.0d0
            g2 = 0.0d0
            g3 = 0.0d0
            exist = .false.
            read (string,*,err=140,end=140)  ia,ib,g1,g2
            exist = .true.
  140       continue
            read (string,*,err=150,end=150)  ia,ib,g1,g2,g3
  150       continue
            if (.not. exist) then
               xcm = 0.0d0
               ycm = 0.0d0
               zcm = 0.0d0
               do j = igrp(1,ia), igrp(2,ia)
                  k = kgrp(j)
                  weigh = mass(k)
                  xcm = xcm + x(k)*weigh
                  ycm = ycm + y(k)*weigh
                  zcm = zcm + z(k)*weigh
               end do
               weigh = max(1.0d0,grpmass(ia))
               xr = xcm / weigh
               yr = ycm / weigh
               zr = zcm / weigh
               xcm = 0.0d0
               ycm = 0.0d0
               zcm = 0.0d0
               do j = igrp(1,ib), igrp(2,ib)
                  k = kgrp(j)
                  weigh = mass(k)
                  xcm = xcm + x(k)*weigh
                  ycm = ycm + y(k)*weigh
                  zcm = zcm + z(k)*weigh
               end do
               weigh = max(1.0d0,grpmass(ib))
               xr = xr - xcm/weigh
               yr = yr - ycm/weigh
               zr = zr - zcm/weigh
               intermol = (molcule(kgrp(igrp(1,ia))) .ne.
     &                     molcule(kgrp(igrp(1,ib))))
               if (use_bounds .and. intermol)  call image (xr,yr,zr)
               g2 = sqrt(xr*xr + yr*yr + zr*zr)
            end if
            if (g3 .eq. 0.0d0)  g3 = g2
            ngfix = ngfix + 1
            igfix(1,ngfix) = ia
            igfix(2,ngfix) = ib
            gfix(1,ngfix) = g1
            gfix(2,ngfix) = g2
            gfix(3,ngfix) = g3
            if (ngfix .gt. maxfix) then
               write (iout,160)
  160          format (/,' KGEOM  --  Too many Group Restraints;',
     &                    ' Increase MAXFIX')
               call fatal
            end if
c
c     maintain chirality as found in the original input structure
c
         else if (keyword(1:18) .eq. 'ENFORCE-CHIRALITY ') then
            do j = 1, n
               if (n12(j) .eq. 4) then
                  ia = i12(1,j)
                  ib = i12(2,j)
                  ic = i12(3,j)
                  id = i12(4,j)
                  keep = .true.
                  if (n12(ia) .eq. 1) then
                     if (type(ia) .eq. type(ib))  keep = .false.
                     if (type(ia) .eq. type(ic))  keep = .false.
                     if (type(ia) .eq. type(id))  keep = .false.
                  else if (n12(ib) .eq. 1) then
                     if (type(ib) .eq. type(ic))  keep = .false.
                     if (type(ib) .eq. type(id))  keep = .false.
                  else if (n12(ic) .eq. 1) then
                     if (type(ic) .eq. type(id))  keep = .false.
                  end if
                  if (keep) then
                     nchir = nchir + 1
                     ichir(1,nchir) = ia
                     ichir(2,nchir) = ib
                     ichir(3,nchir) = ic
                     ichir(4,nchir) = id
                     xad = x(ia) - x(id)
                     yad = y(ia) - y(id)
                     zad = z(ia) - z(id)
                     xbd = x(ib) - x(id)
                     ybd = y(ib) - y(id)
                     zbd = z(ib) - z(id)
                     xcd = x(ic) - x(id)
                     ycd = y(ic) - y(id)
                     zcd = z(ic) - z(id)
                     c1 = ybd*zcd - zbd*ycd
                     c2 = ycd*zad - zcd*yad
                     c3 = yad*zbd - zad*ybd
                     vol = xad*c1 + xbd*c2 + xcd*c3
                     ratio = abs(vol/(xad*xbd*xcd))
                     chir(1,nchir) = 10.0d0
                     if (ratio .gt. 0.1d0) then
                        chir(2,nchir) = 0.5d0 * vol
                        chir(3,nchir) = 2.0d0 * vol
                     else
                        chir(2,nchir) = -2.0d0 * abs(vol)
                        chir(3,nchir) = 2.0d0 * abs(vol)
                     end if
                  end if
               end if
            end do
c
c     setup any shallow Gaussian basin restraint between atoms
c
         else if (keyword(1:6) .eq. 'BASIN ') then
            depth = 0.0d0
            width = 0.0d0
            read (string,*,err=170,end=170)  depth,width
  170       continue
            use_basin = .true.
            if (depth .eq. 0.0d0)  use_basin = .false.
            if (width .eq. 0.0d0)  use_basin = .false.
            if (depth .gt. 0.0d0)  depth = -depth
c
c     setup any spherical droplet restraint between atoms
c
         else if (keyword(1:5) .eq. 'WALL ') then
            rwall = 0.0d0
            read (string,*,err=180,end=180)  rwall
  180       continue
            if (rwall .gt. 0.0d0)  use_wall = .true.
         end if
      end do
c
c     turn on the geometric restraint potential if it is used
c
      use_geom = .false.
      if (npfix .ne. 0)  use_geom = .true.
      if (ndfix .ne. 0)  use_geom = .true.
      if (nafix .ne. 0)  use_geom = .true.
      if (ntfix .ne. 0)  use_geom = .true.
      if (ngfix .ne. 0)  use_geom = .true.
      if (nchir .ne. 0)  use_geom = .true.
      if (use_basin)  use_geom = .true.
      if (use_wall)  use_geom = .true.
      return
      end
