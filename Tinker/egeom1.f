c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine egeom1  --  restraint energy & derivatives  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "egeom1" calculates the energy and first derivatives
c     with respect to Cartesian coordinates due to restraints
c     on positions, distances, angles and torsions as well as
c     Gaussian basin and spherical droplet restraints
c
c
      subroutine egeom1
      use sizes
      use atomid
      use atoms
      use bound
      use deriv
      use energi
      use group
      use inter
      use molcul
      use math
      use restrn
      use usage
      use virial
      implicit none
      integer i,j,k
      integer ia,ib,ic,id
      real*8 e,xr,yr,zr,fgrp
      real*8 de,dt,dt2,deddt
      real*8 r,r2,r6,r12
      real*8 dedx,dedy,dedz
      real*8 angle,target
      real*8 dot,force
      real*8 cosine,sine
      real*8 terma,termc
      real*8 rab2,rcb2
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xab,yab,zab
      real*8 xba,yba,zba
      real*8 xcb,ycb,zcb
      real*8 xdc,ydc,zdc
      real*8 xca,yca,zca
      real*8 xdb,ydb,zdb
      real*8 xad,yad,zad
      real*8 xbd,ybd,zbd
      real*8 xcd,ycd,zcd
      real*8 xp,yp,zp,rp
      real*8 xt,yt,zt
      real*8 xu,yu,zu
      real*8 xtu,ytu,ztu
      real*8 rt2,ru2,rtru
      real*8 rcb,dedphi
      real*8 dedxt,dedyt,dedzt
      real*8 dedxu,dedyu,dedzu
      real*8 dedxia,dedyia,dedzia
      real*8 dedxib,dedyib,dedzib
      real*8 dedxic,dedyic,dedzic
      real*8 dedxid,dedyid,dedzid
      real*8 df1,df2
      real*8 af1,af2
      real*8 tf1,tf2,t1,t2
      real*8 gf1,gf2
      real*8 weigh,ratio
      real*8 weigha,weighb
      real*8 xcm,ycm,zcm
      real*8 cf1,cf2,vol
      real*8 c1,c2,c3
      real*8 xi,yi,zi,ri
      real*8 a,b,buffer,term
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      logical proceed,intermol
c
c
c     zero out the restraint energy term and first derivatives
c
      eg = 0.0d0
      do i = 1, n
         deg(1,i) = 0.0d0
         deg(2,i) = 0.0d0
         deg(3,i) = 0.0d0
      end do
c
c     get energy and derivatives for position restraint terms
c
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
            if (r .eq. 0.0d0)  r = 1.0d0
            de = 2.0d0 * force * dt / r
c
c     scale the interaction based on its group membership
c
            if (use_group) then
               e = e * fgrp
               de = de * fgrp
            end if
c
c     compute chain rule terms needed for derivatives
c
            dedx = de * xr
            dedy = de * yr
            dedz = de * zr
c
c     increment the total energy and first derivatives
c
            eg = eg + e
            deg(1,ia) = deg(1,ia) + dedx
            deg(2,ia) = deg(2,ia) + dedy
            deg(3,ia) = deg(3,ia) + dedz
c
c     increment the internal virial tensor components
c
            vxx = xr * dedx
            vyx = yr * dedx
            vzx = zr * dedx
            vyy = yr * dedy
            vzy = zr * dedy
            vzz = zr * dedz
            vir(1,1) = vir(1,1) + vxx
            vir(2,1) = vir(2,1) + vyx
            vir(3,1) = vir(3,1) + vzx
            vir(1,2) = vir(1,2) + vyx
            vir(2,2) = vir(2,2) + vyy
            vir(3,2) = vir(3,2) + vzy
            vir(1,3) = vir(1,3) + vzx
            vir(2,3) = vir(2,3) + vzy
            vir(3,3) = vir(3,3) + vzz
         end if
      end do
c
c     get energy and derivatives for distance restraint terms
c
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
            if (r .eq. 0.0d0)  r = 1.0d0
            de = 2.0d0 * force * dt / r
c
c     scale the interaction based on its group membership
c
            if (use_group) then
               e = e * fgrp
               de = de * fgrp
            end if
c
c     compute chain rule terms needed for derivatives
c
            dedx = de * xr
            dedy = de * yr
            dedz = de * zr
c
c     increment the total energy and first derivatives
c
            eg = eg + e
            deg(1,ia) = deg(1,ia) + dedx
            deg(2,ia) = deg(2,ia) + dedy
            deg(3,ia) = deg(3,ia) + dedz
            deg(1,ib) = deg(1,ib) - dedx
            deg(2,ib) = deg(2,ib) - dedy
            deg(3,ib) = deg(3,ib) - dedz
c
c     increment the internal virial tensor components
c
            vxx = xr * dedx
            vyx = yr * dedx
            vzx = zr * dedx
            vyy = yr * dedy
            vzy = zr * dedy
            vzz = zr * dedz
            vir(1,1) = vir(1,1) + vxx
            vir(2,1) = vir(2,1) + vyx
            vir(3,1) = vir(3,1) + vzx
            vir(1,2) = vir(1,2) + vyx
            vir(2,2) = vir(2,2) + vyy
            vir(3,2) = vir(3,2) + vzy
            vir(1,3) = vir(1,3) + vzx
            vir(2,3) = vir(2,3) + vzy
            vir(3,3) = vir(3,3) + vzz
c
c     increment the total intermolecular energy
c
            if (intermol) then
               einter = einter + e
            end if
         end if
      end do
c
c     get energy and derivatives for angle restraint terms
c
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
               xp = ycb*zab - zcb*yab
               yp = zcb*xab - xcb*zab
               zp = xcb*yab - ycb*xab
               rp = sqrt(xp*xp + yp*yp + zp*zp)
               rp = max(rp,0.000001d0)
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
               deddt = 2.0d0 * force * dt * radian
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  e = e * fgrp
                  deddt = deddt * fgrp
               end if
c
c     compute derivative components for this interaction
c
               terma = -deddt / (rab2*rp)
               termc = deddt / (rcb2*rp)
               dedxia = terma * (yab*zp-zab*yp)
               dedyia = terma * (zab*xp-xab*zp)
               dedzia = terma * (xab*yp-yab*xp)
               dedxic = termc * (ycb*zp-zcb*yp)
               dedyic = termc * (zcb*xp-xcb*zp)
               dedzic = termc * (xcb*yp-ycb*xp)
               dedxib = -dedxia - dedxic
               dedyib = -dedyia - dedyic
               dedzib = -dedzia - dedzic
c
c     increment the overall energy term and derivatives
c
               eg = eg + e
               deg(1,ia) = deg(1,ia) + dedxia
               deg(2,ia) = deg(2,ia) + dedyia
               deg(3,ia) = deg(3,ia) + dedzia
               deg(1,ib) = deg(1,ib) + dedxib
               deg(2,ib) = deg(2,ib) + dedyib
               deg(3,ib) = deg(3,ib) + dedzib
               deg(1,ic) = deg(1,ic) + dedxic
               deg(2,ic) = deg(2,ic) + dedyic
               deg(3,ic) = deg(3,ic) + dedzic
c
c     increment the internal virial tensor components
c
               vxx = xab*dedxia + xcb*dedxic
               vyx = yab*dedxia + ycb*dedxic
               vzx = zab*dedxia + zcb*dedxic
               vyy = yab*dedyia + ycb*dedyic
               vzy = zab*dedyia + zcb*dedyic
               vzz = zab*dedzia + zcb*dedzic
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vyx
               vir(3,1) = vir(3,1) + vzx
               vir(1,2) = vir(1,2) + vyx
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vzy
               vir(1,3) = vir(1,3) + vzx
               vir(2,3) = vir(2,3) + vzy
               vir(3,3) = vir(3,3) + vzz
c
c     increment the total intermolecular energy
c
               if (intermol) then
                  einter = einter + e
               end if
            end if
         end if
      end do
c
c     get energy and derivatives for torsion restraint terms
c
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
               dedphi = 2.0d0 * radian * force * dt
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  e = e * fgrp
                  dedphi = dedphi * fgrp
               end if
c
c     chain rule terms for first derivative components
c
               xca = xic - xia
               yca = yic - yia
               zca = zic - zia
               xdb = xid - xib
               ydb = yid - yib
               zdb = zid - zib
               dedxt = dedphi * (yt*zcb - ycb*zt) / (rt2*rcb)
               dedyt = dedphi * (zt*xcb - zcb*xt) / (rt2*rcb)
               dedzt = dedphi * (xt*ycb - xcb*yt) / (rt2*rcb)
               dedxu = -dedphi * (yu*zcb - ycb*zu) / (ru2*rcb)
               dedyu = -dedphi * (zu*xcb - zcb*xu) / (ru2*rcb)
               dedzu = -dedphi * (xu*ycb - xcb*yu) / (ru2*rcb)
c
c     compute derivative components for this interaction
c
               dedxia = zcb*dedyt - ycb*dedzt
               dedyia = xcb*dedzt - zcb*dedxt
               dedzia = ycb*dedxt - xcb*dedyt
               dedxib = yca*dedzt - zca*dedyt + zdc*dedyu - ydc*dedzu
               dedyib = zca*dedxt - xca*dedzt + xdc*dedzu - zdc*dedxu
               dedzib = xca*dedyt - yca*dedxt + ydc*dedxu - xdc*dedyu
               dedxic = zba*dedyt - yba*dedzt + ydb*dedzu - zdb*dedyu
               dedyic = xba*dedzt - zba*dedxt + zdb*dedxu - xdb*dedzu
               dedzic = yba*dedxt - xba*dedyt + xdb*dedyu - ydb*dedxu
               dedxid = zcb*dedyu - ycb*dedzu
               dedyid = xcb*dedzu - zcb*dedxu
               dedzid = ycb*dedxu - xcb*dedyu
c
c     increment the overall energy term and derivatives
c
               eg = eg + e
               deg(1,ia) = deg(1,ia) + dedxia
               deg(2,ia) = deg(2,ia) + dedyia
               deg(3,ia) = deg(3,ia) + dedzia
               deg(1,ib) = deg(1,ib) + dedxib
               deg(2,ib) = deg(2,ib) + dedyib
               deg(3,ib) = deg(3,ib) + dedzib
               deg(1,ic) = deg(1,ic) + dedxic
               deg(2,ic) = deg(2,ic) + dedyic
               deg(3,ic) = deg(3,ic) + dedzic
               deg(1,id) = deg(1,id) + dedxid
               deg(2,id) = deg(2,id) + dedyid
               deg(3,id) = deg(3,id) + dedzid
c
c     increment the internal virial tensor components
c
               vxx = xcb*(dedxic+dedxid) - xba*dedxia + xdc*dedxid
               vyx = ycb*(dedxic+dedxid) - yba*dedxia + ydc*dedxid
               vzx = zcb*(dedxic+dedxid) - zba*dedxia + zdc*dedxid
               vyy = ycb*(dedyic+dedyid) - yba*dedyia + ydc*dedyid
               vzy = zcb*(dedyic+dedyid) - zba*dedyia + zdc*dedyid
               vzz = zcb*(dedzic+dedzid) - zba*dedzia + zdc*dedzid
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vyx
               vir(3,1) = vir(3,1) + vzx
               vir(1,2) = vir(1,2) + vyx
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vzy
               vir(1,3) = vir(1,3) + vzx
               vir(2,3) = vir(2,3) + vzy
               vir(3,3) = vir(3,3) + vzz
c
c     increment the total intermolecular energy
c
               if (molcule(ia).ne.molcule(ib) .or.
     &             molcule(ia).ne.molcule(ic) .or.
     &             molcule(ia).ne.molcule(id)) then
                  einter = einter + e
               end if
            end if
         end if
      end do
c
c     get energy and derivatives for group distance restraint terms
c
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
         if (r .eq. 0.0d0)  r = 1.0d0
         de = 2.0d0 * force * dt / r
c
c     compute chain rule terms needed for derivatives
c
         dedx = de * xr
         dedy = de * yr
         dedz = de * zr
c
c     increment the total energy and first derivatives
c
         eg = eg + e
         do j = igrp(1,ia), igrp(2,ia)
            k = kgrp(j)
            ratio = mass(k) / weigha
            deg(1,k) = deg(1,k) + dedx*ratio
            deg(2,k) = deg(2,k) + dedy*ratio
            deg(3,k) = deg(3,k) + dedz*ratio
         end do
         do j = igrp(1,ib), igrp(2,ib)
            k = kgrp(j)
            ratio = mass(k) / weighb
            deg(1,k) = deg(1,k) - dedx*ratio
            deg(2,k) = deg(2,k) - dedy*ratio
            deg(3,k) = deg(3,k) - dedz*ratio
         end do
c
c     increment the internal virial tensor components
c
         vxx = xr * dedx
         vyx = yr * dedx
         vzx = zr * dedx
         vyy = yr * dedy
         vzy = zr * dedy
         vzz = zr * dedz
         vir(1,1) = vir(1,1) + vxx
         vir(2,1) = vir(2,1) + vyx
         vir(3,1) = vir(3,1) + vzx
         vir(1,2) = vir(1,2) + vyx
         vir(2,2) = vir(2,2) + vyy
         vir(3,2) = vir(3,2) + vzy
         vir(1,3) = vir(1,3) + vzx
         vir(2,3) = vir(2,3) + vzy
         vir(3,3) = vir(3,3) + vzz
c
c     increment the total intermolecular energy
c
         if (intermol) then
            einter = einter + e
         end if
      end do
c
c     get energy and derivatives for chirality restraint terms
c
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
            deddt = 2.0d0 * force * dt
c
c     scale the interaction based on its group membership
c
            if (use_group) then
               e = e * fgrp
               deddt = deddt * fgrp
            end if
c
c     compute derivative components for this interaction
c
            dedxia = deddt * (ybd*zcd - zbd*ycd)
            dedyia = deddt * (zbd*xcd - xbd*zcd)
            dedzia = deddt * (xbd*ycd - ybd*xcd)
            dedxib = deddt * (zad*ycd - yad*zcd)
            dedyib = deddt * (xad*zcd - zad*xcd)
            dedzib = deddt * (yad*xcd - xad*ycd)
            dedxic = deddt * (yad*zbd - zad*ybd)
            dedyic = deddt * (zad*xbd - xad*zbd)
            dedzic = deddt * (xad*ybd - yad*xbd)
            dedxid = -dedxia - dedxib - dedxic
            dedyid = -dedyia - dedyib - dedyic
            dedzid = -dedzia - dedzib - dedzic
c
c     increment the overall energy term and derivatives
c
            eg = eg + e
            deg(1,ia) = deg(1,ia) + dedxia
            deg(2,ia) = deg(2,ia) + dedyia
            deg(3,ia) = deg(3,ia) + dedzia
            deg(1,ib) = deg(1,ib) + dedxib
            deg(2,ib) = deg(2,ib) + dedyib
            deg(3,ib) = deg(3,ib) + dedzib
            deg(1,ic) = deg(1,ic) + dedxic
            deg(2,ic) = deg(2,ic) + dedyic
            deg(3,ic) = deg(3,ic) + dedzic
            deg(1,id) = deg(1,id) + dedxid
            deg(2,id) = deg(2,id) + dedyid
            deg(3,id) = deg(3,id) + dedzid
c
c     increment the internal virial tensor components
c
            vxx = xad*dedxia + xbd*dedxib + xcd*dedxic
            vyx = yad*dedxia + ybd*dedxib + ycd*dedxic
            vzx = zad*dedxia + zbd*dedxib + zcd*dedxic
            vyy = yad*dedyia + ybd*dedyib + ycd*dedyic
            vzy = zad*dedyia + zbd*dedyib + zcd*dedyic
            vzz = zad*dedzia + zbd*dedzib + zcd*dedzic
            vir(1,1) = vir(1,1) + vxx
            vir(2,1) = vir(2,1) + vyx
            vir(3,1) = vir(3,1) + vzx
            vir(1,2) = vir(1,2) + vyx
            vir(2,2) = vir(2,2) + vyy
            vir(3,2) = vir(3,2) + vzy
            vir(1,3) = vir(1,3) + vzx
            vir(2,3) = vir(2,3) + vzy
            vir(3,3) = vir(3,3) + vzz
c
c     increment the total intermolecular energy
c
            if (molcule(ia).ne.molcule(ib) .or.
     &          molcule(ia).ne.molcule(ic) .or.
     &          molcule(ia).ne.molcule(id)) then
               einter = einter + e
            end if
         end if
      end do
c
c     get energy and derivatives for a Gaussian basin restraint
c
      if (use_basin) then
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
                  de = -2.0d0 * width * e
                  e = e - depth
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     de = de * fgrp
                  end if
c
c     compute chain rule terms needed for derivatives
c
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
c
c     increment the overall energy term and derivatives
c
                  eg = eg + e
                  deg(1,i) = deg(1,i) + dedx
                  deg(2,i) = deg(2,i) + dedy
                  deg(3,i) = deg(3,i) + dedz
                  deg(1,k) = deg(1,k) - dedx
                  deg(2,k) = deg(2,k) - dedy
                  deg(3,k) = deg(3,k) - dedz
c
c     increment the internal virial tensor components
c
                  vxx = xr * dedx
                  vyx = yr * dedx
                  vzx = zr * dedx
                  vyy = yr * dedy
                  vzy = zr * dedy
                  vzz = zr * dedz
                  vir(1,1) = vir(1,1) + vxx
                  vir(2,1) = vir(2,1) + vyx
                  vir(3,1) = vir(3,1) + vzx
                  vir(1,2) = vir(1,2) + vyx
                  vir(2,2) = vir(2,2) + vyy
                  vir(3,2) = vir(3,2) + vzy
                  vir(1,3) = vir(1,3) + vzx
                  vir(2,3) = vir(2,3) + vzy
                  vir(3,3) = vir(3,3) + vzz
               end if
            end do
         end do
      end if
c
c     get energy and derivatives for a spherical droplet restraint
c
      if (use_wall) then
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
               if (ri .eq. 0.0d0)  ri = 1.0d0
               de = (12.0d0*a/r12 - 6.0d0*b/r6) / (r*ri)
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  e = e * fgrp
                  de = de * fgrp
               end if
c
c     compute chain rule terms needed for derivatives
c
               dedx = de * xi
               dedy = de * yi
               dedz = de * zi
c
c     increment the overall energy term and derivatives
c
               eg = eg + e
               deg(1,i) = deg(1,i) + dedx
               deg(2,i) = deg(2,i) + dedy
               deg(3,i) = deg(3,i) + dedz
c
c     increment the internal virial tensor components
c
               xr = r * xi/ri
               yr = r * yi/ri
               zr = r * zi/ri
               vxx = xr * dedx
               vyx = yr * dedx
               vzx = zr * dedx
               vyy = yr * dedy
               vzy = zr * dedy
               vzz = zr * dedz
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vyx
               vir(3,1) = vir(3,1) + vzx
               vir(1,2) = vir(1,2) + vyx
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vzy
               vir(1,3) = vir(1,3) + vzx
               vir(2,3) = vir(2,3) + vzy
               vir(3,3) = vir(3,3) + vzz
            end if
         end do
      end if
      return
      end
