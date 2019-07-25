c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine egeom3  --  restraint energy terms & analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "egeom3" calculates the energy due to restraints on positions,
c     distances, angles and torsions as well as Gaussian basin and
c     droplet restraints; also partitions energy among the atoms
c
c
      subroutine egeom3
      use sizes
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use energi
      use group
      use inform
      use inter
      use iounit
      use math
      use molcul
      use restrn
      use usage
      implicit none
      integer i,j,k
      integer ia,ib,ic,id
      real*8 e,dt,dt2,fgrp
      real*8 xr,yr,zr
      real*8 r,r2,r6,r12
      real*8 angle,target
      real*8 dot,force
      real*8 cosine,sine
      real*8 rab2,rcb2
      real*8 xt,yt,zt
      real*8 xu,yu,zu
      real*8 xtu,ytu,ztu
      real*8 rt2,ru2
      real*8 rtru,rcb
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xab,yab,zab
      real*8 xba,yba,zba
      real*8 xcb,ycb,zcb
      real*8 xdc,ydc,zdc
      real*8 xad,yad,zad
      real*8 xbd,ybd,zbd
      real*8 xcd,ycd,zcd
      real*8 df1,df2
      real*8 af1,af2
      real*8 tf1,tf2,t1,t2
      real*8 gf1,gf2
      real*8 weigh,size
      real*8 weigha,weighb
      real*8 xcm,ycm,zcm
      real*8 cf1,cf2,vol
      real*8 c1,c2,c3
      real*8 xi,yi,zi,ri
      real*8 a,b,buffer,term
      logical proceed,intermol
      logical header,huge
c
c
c     zero out the restraint energy and partitioning terms
c
      neg = 0
      eg = 0.0d0
      do i = 1, n
         aeg(i) = 0.0d0
      end do
c
c     compute the energy for position restraint terms
c
      header = .true.
      do i = 1, npfix
         ia = ipfix(i)
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,ia,0,0,0,0,0)
         if (proceed)  proceed = (use(ia))
         if (proceed) then
            xr = 0.0d0
            yr = 0.0d0
            zr = 0.0d0
            if (kpfix(1,i) .ne. 0)  xr = x(ia) - xpfix(i)
            if (kpfix(2,i) .ne. 0)  yr = y(ia) - ypfix(i)
            if (kpfix(3,i) .ne. 0)  zr = z(ia) - zpfix(i)
            if (use_bounds)  call image (xr,yr,zr)
            r = sqrt(xr*xr + yr*yr + zr*zr)
            force = pfix(1,i)
            dt = max(0.0d0,r-pfix(2,i))
            dt2 = dt * dt
            e = force * dt2
            if (use_group)  e = e * fgrp
            neg = neg + 1
            eg = eg + e
            aeg(ia) = aeg(ia) + e
            huge = (e .gt. 10.0d0)
            if (debug .or. (verbose.and.huge)) then
               if (header) then
                  header = .false.
                  write (iout,10)
   10             format (/,' Individual Atomic Position Restraint',
     &                       ' Terms :',
     &                    //,' Type',9x,'Atom Name',13x,'Target',
     &                       ' Position',7x,'Distance',6x,'Energy',/)
               end if
               dt = sqrt(dt2)
               if (kpfix(2,i).eq.0 .and. kpfix(3,i).eq.0) then
                  write (iout,20)  ia,name(ia),xpfix(i),dt,e
   20             format (' Position',2x,i7,'-',a3,4x,f10.4,
     &                       5x,'----',6x,'----',1x,f10.4,f12.4)
               else if (kpfix(1,i).eq.0 .and. kpfix(3,i).eq.0) then
                  write (iout,30)  ia,name(ia),ypfix(i),dt,e
   30             format (' Position',2x,i7,'-',a3,9x,'----',1x,
     &                       f10.4,5x,'----',1x,f10.4,f12.4)
               else if (kpfix(1,i).eq.0 .and. kpfix(2,i).eq.0) then
                  write (iout,40)  ia,name(ia),zpfix(i),dt,e
   40             format (' Position',2x,i7,'-',a3,9x,'----',
     &                       6x,'----',1x,2f10.4,f12.4)
               else
                  write (iout,50)  ia,name(ia),xpfix(i),ypfix(i),
     &                             zpfix(i),dt,e
   50             format (' Position',2x,i7,'-',a3,4x,4f10.4,f12.4)
               end if
            end if
         end if
      end do
c
c     compute the energy for distance restraint terms
c
      header = .true.
      do i = 1, ndfix
         ia = idfix(1,i)
         ib = idfix(2,i)
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,ia,ib,0,0,0,0)
         if (proceed)  proceed = (use(ia) .or. use(ib))
         if (proceed) then
            xr = x(ia) - x(ib)
            yr = y(ia) - y(ib)
            zr = z(ia) - z(ib)
            intermol = (molcule(ia) .ne. molcule(ib))
            if (use_bounds .and. intermol)  call image (xr,yr,zr)
            r = sqrt(xr*xr + yr*yr + zr*zr)
            force = dfix(1,i)
            df1 = dfix(2,i)
            df2 = dfix(3,i)
            target = r
            if (r .lt. df1)  target = df1
            if (r .gt. df2)  target = df2
            dt = r - target
            dt2 = dt * dt
            e = force * dt2
            if (use_group)  e = e * fgrp
            neg = neg + 1
            eg = eg + e
            aeg(ia) = aeg(ia) + 0.5d0*e
            aeg(ib) = aeg(ib) + 0.5d0*e
            if (intermol) then
               einter = einter + e
            end if
            huge = (e .gt. 10.0d0)
            if (debug .or. (verbose.and.huge)) then
               if (header) then
                  header = .false.
                  write (iout,60)
   60             format (/,' Individual Interatomic Distance',
     &                       ' Restraint Terms :',
     &                    //,' Type',14x,'Atom Names',16x,'Ideal Range',
     &                       4x,'Actual',6x,'Energy',/)
               end if
               write (iout,70)  ia,name(ia),ib,name(ib),df1,df2,r,e
   70          format (' Distance',2x,2(i7,'-',a3),
     &                    7x,2f8.2,f10.4,f12.4)
            end if
         end if
      end do
c
c     compute the energy for angle restraint terms
c
      header = .true.
      do i = 1, nafix
         ia = iafix(1,i)
         ib = iafix(2,i)
         ic = iafix(3,i)
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,ia,ib,ic,0,0,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or. use(ic))
         if (proceed) then
            xia = x(ia)
            yia = y(ia)
            zia = z(ia)
            xib = x(ib)
            yib = y(ib)
            zib = z(ib)
            xic = x(ic)
            yic = y(ic)
            zic = z(ic)
            xab = xia - xib
            yab = yia - yib
            zab = zia - zib
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            rab2 = xab*xab + yab*yab + zab*zab
            rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
            if (rab2.ne.0.0d0 .and. rcb2.ne.0.0d0) then
               dot = xab*xcb + yab*ycb + zab*zcb
               cosine = dot / sqrt(rab2*rcb2)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
               force = afix(1,i)
               af1 = afix(2,i)
               af2 = afix(3,i)
               target = angle
               if (angle .lt. af1)  target = af1
               if (angle .gt. af2)  target = af2
               dt = angle - target
               dt2 = dt * dt
               e = force * dt2
               if (use_group)  e = e * fgrp
               neg = neg + 1
               eg = eg + e
               aeg(ib) = aeg(ib) + e
               if (intermol) then
                  einter = einter + e
               end if
               huge = (e .gt. 10.0d0)
               if (debug .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,80)
   80                format (/,' Individual Interatomic Angle',
     &                          ' Restraint Terms :',
     &                       //,' Type',14x,'Atom Numbers',14x,'Ideal',
     &                          ' Range',4x,'Actual',6x,'Energy',/)
                  end if
                  write (iout,90)  ia,ib,ic,af1,af2,angle,e
   90             format (' Angle',8x,3i6,8x,2f8.2,f10.4,f12.4)
               end if
            end if
         end if
      end do
c
c     compute the energy for torsional restraint terms
c
      header = .true.
      do i = 1, ntfix
         ia = itfix(1,i)
         ib = itfix(2,i)
         ic = itfix(3,i)
         id = itfix(4,i)
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,ia,ib,ic,id,0,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
         if (proceed) then
            xia = x(ia)
            yia = y(ia)
            zia = z(ia)
            xib = x(ib)
            yib = y(ib)
            zib = z(ib)
            xic = x(ic)
            yic = y(ic)
            zic = z(ic)
            xid = x(id)
            yid = y(id)
            zid = z(id)
            xba = xib - xia
            yba = yib - yia
            zba = zib - zia
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdc = xid - xic
            ydc = yid - yic
            zdc = zid - zic
            xt = yba*zcb - ycb*zba
            yt = zba*xcb - zcb*xba
            zt = xba*ycb - xcb*yba
            xu = ycb*zdc - ydc*zcb
            yu = zcb*xdc - zdc*xcb
            zu = xcb*ydc - xdc*ycb
            xtu = yt*zu - yu*zt
            ytu = zt*xu - zu*xt
            ztu = xt*yu - xu*yt
            rt2 = xt*xt + yt*yt + zt*zt
            ru2 = xu*xu + yu*yu + zu*zu
            rtru = sqrt(rt2 * ru2)
            if (rtru .ne. 0.0d0) then
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
               if (sine .lt. 0.0d0)  angle = -angle
               force = tfix(1,i)
               tf1 = tfix(2,i)
               tf2 = tfix(3,i)
               if (angle.gt.tf1 .and. angle.lt.tf2) then
                  target = angle
               else if (angle.gt.tf1 .and. tf1.gt.tf2) then
                  target = angle
               else if (angle.lt.tf2 .and. tf1.gt.tf2) then
                  target = angle
               else
                  t1 = angle - tf1
                  t2 = angle - tf2
                  if (t1 .gt. 180.0d0) then
                     t1 = t1 - 360.0d0
                  else if (t1 .lt. -180.0d0) then
                     t1 = t1 + 360.0d0
                  end if
                  if (t2 .gt. 180.0d0) then
                     t2 = t2 - 360.0d0
                  else if (t2 .lt. -180.0d0) then
                     t2 = t2 + 360.0d0
                  end if
                  if (abs(t1) .lt. abs(t2)) then
                     target = tf1
                  else
                     target = tf2
                  end if
               end if
               dt = angle - target
               if (dt .gt. 180.0d0) then
                  dt = dt - 360.0d0
               else if (dt .lt. -180.0d0) then
                  dt = dt + 360.0d0
               end if
               dt2 = dt * dt
               e = force * dt2
               if (use_group)  e = e * fgrp
               neg = neg + 1
               eg = eg + e
               aeg(ib) = aeg(ib) + 0.5d0*e
               aeg(ic) = aeg(ic) + 0.5d0*e
               if (molcule(ia).ne.molcule(ib) .or.
     &             molcule(ia).ne.molcule(ic) .or.
     &             molcule(ia).ne.molcule(id)) then
                  einter = einter + e
               end if
               huge = (e .gt. 10.0d0)
               if (debug .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,100)
  100                format (/,' Individual Torsional Angle Restraint
     &',                          ' Terms :',
     &                       //,' Type',14x,'Atom Numbers',14x,'Ideal',
     &                          ' Range',4x,'Actual',6x,'Energy',/)
                  end if
                  write (iout,110)  ia,ib,ic,id,tf1,tf2,angle,e
  110             format (' Torsion',4x,4i6,4x,2f8.2,f10.4,f12.4)
               end if
            end if
         end if
      end do
c
c     compute the energy for group distance restraint terms
c
      header = .true.
      do i = 1, ngfix
         ia = igfix(1,i)
         ib = igfix(2,i)
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
         weigha = max(1.0d0,grpmass(ia))
         xr = xcm / weigha
         yr = ycm / weigha
         zr = zcm / weigha
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
         weighb = max(1.0d0,grpmass(ib))
         xr = xr - xcm/weighb
         yr = yr - ycm/weighb
         zr = zr - zcm/weighb
         intermol = (molcule(kgrp(igrp(1,ia))) .ne.
     &               molcule(kgrp(igrp(1,ib))))
         if (use_bounds .and. intermol)  call image (xr,yr,zr)
         r = sqrt(xr*xr + yr*yr + zr*zr)
         force = gfix(1,i)
         gf1 = gfix(2,i)
         gf2 = gfix(3,i)
         target = r
         if (r .lt. gf1)  target = gf1
         if (r .gt. gf2)  target = gf2
         dt = r - target
         dt2 = dt * dt
         e = force * dt2
         neg = neg + 1
         eg = eg + e
         size = dble(igrp(2,ia) - igrp(1,ia) + 1)
         do j = igrp(1,ia), igrp(2,ia)
            k = kgrp(j)
            aeg(k) = aeg(k) + 0.5d0*e/size
         end do
         size = dble(igrp(2,ib) - igrp(1,ib) + 1)
         do j = igrp(1,ib), igrp(2,ib)
            k = kgrp(j)
            aeg(k) = aeg(k) + 0.5d0*e/size
         end do
         if (intermol) then
            einter = einter + e
         end if
         huge = (e .gt. 10.0d0)
         if (debug .or. (verbose.and.huge)) then
            if (header) then
               header = .false.
               write (iout,120)
  120          format (/,' Individual Intergroup Distance',
     &                    ' Restraint Terms :',
     &                 //,' Type',13x,'Group Numbers',14x,'Ideal Range',
     &                    4x,'Actual',6x,'Energy',/)
            end if
            write (iout,130)  ia,ib,gf1,gf2,r,e
  130       format (' Distance',7x,2i7,10x,2f8.2,f10.4,f12.4)
         end if
      end do
c
c     compute the energy for chirality restraint terms
c
      header = .true.
      do i = 1, nchir
         ia = ichir(1,i)
         ib = ichir(2,i)
         ic = ichir(3,i)
         id = ichir(4,i)
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,ia,ib,ic,id,0,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
         if (proceed) then
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
            force = chir(1,i)
            cf1 = chir(2,i)
            cf2 = chir(3,i)
            target = vol
            if (vol .lt. min(cf1,cf2))  target = min(cf1,cf2)
            if (vol .gt. max(cf1,cf2))  target = max(cf1,cf2)
            dt = vol - target
            dt2 = dt * dt
            e = force * dt2
            if (use_group)  e = e * fgrp
            neg = neg + 1
            eg = eg + e
            aeg(ia) = aeg(ia) + 0.25d0*e
            aeg(ib) = aeg(ib) + 0.25d0*e
            aeg(ic) = aeg(ic) + 0.25d0*e
            aeg(id) = aeg(id) + 0.25d0*e
            huge = (e .gt. 10.0d0)
            if (debug .or. (verbose.and.huge)) then
               if (header) then
                  header = .false.
                  write (iout,140)
  140             format (/,' Individual Chirality Restraint Terms :'
     &,                    //,' Type',14x,'Atom Numbers',14x,'Ideal',
     &                       ' Range',4x,'Actual',6x,'Energy',/)
               end if
               write (iout,150)  ia,ib,ic,id,cf1,cf2,vol,e
  150          format (' Chiral',5x,4i6,4x,2f8.2,f10.4,f12.4)
            end if
         end if
      end do
c
c     compute the energy for a Gaussian basin restraint
c
      if (use_basin) then
         header = .true.
         do i = 1, n-1
            xi = x(i)
            yi = y(i)
            zi = z(i)
            do k = i+1, n
               proceed = .true.
               if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
               if (proceed)  proceed = (use(i) .or. use(k))
               if (proceed) then
                  xr = xi - x(k)
                  yr = yi - y(k)
                  zr = zi - z(k)
                  r2 = xr*xr + yr*yr + zr*zr
                  term = -width * r2
                  e = 0.0d0
                  if (term .gt. -50.0d0)  e = depth * exp(term)
                  e = e - depth
                  if (use_group)  e = e * fgrp
                  neg = neg + 1
                  eg = eg + e
                  aeg(i) = aeg(i) + 0.5d0*e
                  aeg(k) = aeg(k) + 0.5d0*e
                  huge = (e .gt. 10.0d0)
                  if (debug .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,160)
  160                   format (/,' Individual Gaussian Basin',
     &                             ' Restraint Terms :',
     &                          //,' Type',14x,'Atom Names',22x,'Ideal',
     &                             4x,'Actual',6x,'Energy',/)
                     end if
                     r = sqrt(r2)
                     write (iout,170)  i,name(i),k,name(k),0.0d0,r,e
  170                format (' Distance',2x,2(i7,'-',a3),
     &                          13x,2f10.4,f12.4)
                  end if
               end if
            end do
         end do
      end if
c
c     compute the energy for a spherical droplet restraint
c
      if (use_wall) then
         header = .true.
         buffer = 2.5d0
         a = 2048.0d0
         b = 64.0d0
         do i = 1, n
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,0,0,0,0,0)
            if (proceed)  proceed = (use(i))
            if (proceed) then
               xi = x(i)
               yi = y(i)
               zi = z(i)
               ri = sqrt(xi**2 + yi**2 + zi**2)
               r = rwall + buffer - ri
               r2 = r * r
               r6 = r2 * r2 * r2
               r12 = r6 * r6
               e = a/r12 - b/r6
               if (use_group)  e = e * fgrp
               neg = neg + 1
               eg = eg + e
               aeg(i) = aeg(i) + e
               huge = (e .gt. 10.0d0)
               if (debug .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,180)
  180                format (/,' Individual Spherical Boundary',
     &                          ' Restraint Terms :',
     &                       //,' Type',14x,'Atom Name',30x,'Distance',
     &                          6x,'Energy',/)
                  end if
                  write (iout,190)  i,name(i),ri,e
  190             format (' Wall',11x,i7,'-',a3,29x,f10.4,f12.4)
               end if
            end if
         end do
      end if
      return
      end
