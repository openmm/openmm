c
c
c     #############################################################
c     ##  COPYRIGHT (C) 2003 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine etortor1  --  torsion-torsion energy & derivs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "etortor1" calculates the torsion-torsion energy and first
c     derivatives with respect to Cartesian coordinates
c
c
      subroutine etortor1
      use sizes
      use atoms
      use bitor
      use bound
      use deriv
      use energi
      use group
      use ktrtor
      use math
      use torpot
      use tortor
      use usage
      use virial
      implicit none
      integer i,k,itortor
      integer pos1,pos2
      integer ia,ib,ic,id,ie
      integer nlo,nhi,nt
      integer xlo,ylo
      real*8 e,fgrp,sign
      real*8 angle1,angle2
      real*8 value1,value2
      real*8 cosine1,cosine2
      real*8 xt,yt,zt,rt2
      real*8 xu,yu,zu,ru2
      real*8 xv,yv,zv,rv2
      real*8 xtu,ytu,ztu,rtru
      real*8 xuv,yuv,zuv,rurv
      real*8 xh,yh,x1l,x1u
      real*8 y1l,y1u
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xie,yie,zie
      real*8 xba,yba,zba
      real*8 xdc,ydc,zdc
      real*8 xcb,ycb,zcb
      real*8 xed,yed,zed
      real*8 rcb,rdc
      real*8 xca,yca,zca
      real*8 xdb,ydb,zdb
      real*8 xec,yec,zec
      real*8 dedang1,dedang2
      real*8 dedxt,dedyt,dedzt
      real*8 dedxu,dedyu,dedzu
      real*8 dedxu2,dedyu2,dedzu2
      real*8 dedxv2,dedyv2,dedzv2
      real*8 dedxia,dedyia,dedzia
      real*8 dedxib,dedyib,dedzib
      real*8 dedxic,dedyic,dedzic
      real*8 dedxid,dedyid,dedzid
      real*8 dedxib2,dedyib2,dedzib2
      real*8 dedxic2,dedyic2,dedzic2
      real*8 dedxid2,dedyid2,dedzid2
      real*8 dedxie2,dedyie2,dedzie2
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 vxx2,vyy2,vzz2
      real*8 vyx2,vzx2,vzy2
      real*8 ftt(4),ft12(4)
      real*8 ft1(4),ft2(4)
      logical proceed
c
c
c     zero out the torsion-torsion energy and first derivatives
c
      ett = 0.0d0
      do i = 1, n
         dett(1,i) = 0.0d0
         dett(2,i) = 0.0d0
         dett(3,i) = 0.0d0
      end do
      if (ntortor .eq. 0)  return
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(ntortor,itt,ibitor,
!$OMP& use,x,y,z,tnx,ttx,tny,tty,tbf,tbx,tby,tbxy,ttorunit,
!$OMP& use_group,use_polymer)
!$OMP& shared(ett,dett,vir)
!$OMP DO reduction(+:ett,dett,vir) schedule(guided)
c
c     calculate the torsion-torsion interaction energy term
c
      do itortor = 1, ntortor
         i = itt(1,itortor)
         k = itt(2,itortor)
         if (itt(3,itortor) .eq. 1) then
            ia = ibitor(1,i)
            ib = ibitor(2,i)
            ic = ibitor(3,i)
            id = ibitor(4,i)
            ie = ibitor(5,i)
         else
            ia = ibitor(5,i)
            ib = ibitor(4,i)
            ic = ibitor(3,i)
            id = ibitor(2,i)
            ie = ibitor(1,i)
         end if
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,ia,ib,ic,id,ie,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or. use(ic)
     &                               .or. use(id) .or. use(ie))
c
c     compute the values of the torsional angles
c
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
            xie = x(ie)
            yie = y(ie)
            zie = z(ie)
            xba = xib - xia
            yba = yib - yia
            zba = zib - zia
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdc = xid - xic
            ydc = yid - yic
            zdc = zid - zic
            xed = xie - xid
            yed = yie - yid
            zed = zie - zid
            if (use_polymer) then
               call image (xba,yba,zba)
               call image (xcb,ycb,zcb)
               call image (xdc,ydc,zdc)
               call image (xed,yed,zed)
            end if
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
            xv = ydc*zed - yed*zdc
            yv = zdc*xed - zed*xdc
            zv = xdc*yed - xed*ydc
            xuv = yu*zv - yv*zu
            yuv = zu*xv - zv*xu
            zuv = xu*yv - xv*yu
            rv2 = xv*xv + yv*yv + zv*zv
            rurv = sqrt(ru2 * rv2)
            if (rtru.ne.0.0d0 .and. rurv.ne.0.0d0) then
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cosine1 = (xt*xu + yt*yu + zt*zu) / rtru
               cosine1 = min(1.0d0,max(-1.0d0,cosine1))
               angle1 = radian * acos(cosine1)
               sign = xba*xu + yba*yu + zba*zu
               if (sign .lt. 0.0d0)  angle1 = -angle1
               value1 = angle1
               rdc = sqrt(xdc*xdc + ydc*ydc + zdc*zdc)
               cosine2 = (xu*xv + yu*yv + zu*zv) / rurv
               cosine2 = min(1.0d0,max(-1.0d0,cosine2))
               angle2 = radian * acos(cosine2)
               sign = xcb*xv + ycb*yv + zcb*zv
               if (sign .lt. 0.0d0)  angle2 = -angle2
               value2 = angle2
c
c     check for inverted chirality at the central atom
c
               call chkttor (ib,ic,id,sign,value1,value2)
c
c     use bicubic interpolation to compute spline values
c
               nlo = 1
               nhi = tnx(k)
               do while (nhi-nlo .gt. 1)
                  nt = (nhi+nlo) / 2
                  if (ttx(nt,k) .gt. value1) then
                     nhi = nt
                  else
                     nlo = nt
                  end if
               end do
               xlo = nlo
               nlo = 1
               nhi = tny(k)
               do while (nhi-nlo .gt. 1)
                  nt = (nhi + nlo)/2
                  if (tty(nt,k) .gt. value2) then
                     nhi = nt
                  else
                     nlo = nt
                  end if
               end do
               ylo = nlo
               xh = ttx(xlo+1,k) - ttx(xlo,k)
               yh = tty(ylo+1,k) - tty(ylo,k)
               x1l = ttx(xlo,k)
               x1u = ttx(xlo+1,k)
               y1l = tty(ylo,k)
               y1u = tty(ylo+1,k)
               pos2 = ylo*tnx(k) + xlo
               pos1 = pos2 - tnx(k)
               ftt(1) = tbf(pos1,k)
               ftt(2) = tbf(pos1+1,k)
               ftt(3) = tbf(pos2+1,k)
               ftt(4) = tbf(pos2,k)
               ft1(1) = tbx(pos1,k)
               ft1(2) = tbx(pos1+1,k)
               ft1(3) = tbx(pos2+1,k)
               ft1(4) = tbx(pos2,k)
               ft2(1) = tby(pos1,k)
               ft2(2) = tby(pos1+1,k)
               ft2(3) = tby(pos2+1,k)
               ft2(4) = tby(pos2,k)
               ft12(1) = tbxy(pos1,k)
               ft12(2) = tbxy(pos1+1,k)
               ft12(3) = tbxy(pos2+1,k)
               ft12(4) = tbxy(pos2,k)
               call bcuint1 (ftt,ft1,ft2,ft12,x1l,x1u,y1l,y1u,
     &                       value1,value2,e,dedang1,dedang2)
               e = ttorunit * e
               dedang1 = sign * ttorunit * radian * dedang1
               dedang2 = sign * ttorunit * radian * dedang2
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  e = e * fgrp
                  dedang1 = dedang1 * fgrp
                  dedang2 = dedang2 * fgrp
               end if
c
c     chain rule terms for first angle derivative components
c
               xca = xic - xia
               yca = yic - yia
               zca = zic - zia
               xdb = xid - xib
               ydb = yid - yib
               zdb = zid - zib
               if (use_polymer) then
                  call image (xca,yca,zca)
                  call image (xdb,ydb,zdb)
               end if
               dedxt = dedang1 * (yt*zcb - ycb*zt) / (rt2*rcb)
               dedyt = dedang1 * (zt*xcb - zcb*xt) / (rt2*rcb)
               dedzt = dedang1 * (xt*ycb - xcb*yt) / (rt2*rcb)
               dedxu = -dedang1 * (yu*zcb - ycb*zu) / (ru2*rcb)
               dedyu = -dedang1 * (zu*xcb - zcb*xu) / (ru2*rcb)
               dedzu = -dedang1 * (xu*ycb - xcb*yu) / (ru2*rcb)
c
c     compute first derivative components for first angle
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
c     chain rule terms for second angle derivative components
c
               xec = xie - xic
               yec = yie - yic
               zec = zie - zic
               if (use_polymer) then
                  call image (xdb,ydb,zdb)
                  call image (xec,yec,zec)
               end if
               dedxu2 = dedang2 * (yu*zdc - ydc*zu) / (ru2*rdc)
               dedyu2 = dedang2 * (zu*xdc - zdc*xu) / (ru2*rdc)
               dedzu2 = dedang2 * (xu*ydc - xdc*yu) / (ru2*rdc)
               dedxv2 = -dedang2 * (yv*zdc - ydc*zv) / (rv2*rdc)
               dedyv2 = -dedang2 * (zv*xdc - zdc*xv) / (rv2*rdc)
               dedzv2 = -dedang2 * (xv*ydc - xdc*yv) / (rv2*rdc)
c
c     compute first derivative components for second angle
c
               dedxib2 = zdc*dedyu2 - ydc*dedzu2
               dedyib2 = xdc*dedzu2 - zdc*dedxu2
               dedzib2 = ydc*dedxu2 - xdc*dedyu2
               dedxic2 = ydb*dedzu2 - zdb*dedyu2
     &                      + zed*dedyv2 - yed*dedzv2
               dedyic2 = zdb*dedxu2 - xdb*dedzu2
     &                      + xed*dedzv2 - zed*dedxv2
               dedzic2 = xdb*dedyu2 - ydb*dedxu2
     &                      + yed*dedxv2 - xed*dedyv2
               dedxid2 = zcb*dedyu2 - ycb*dedzu2
     &                      + yec*dedzv2 - zec*dedyv2
               dedyid2 = xcb*dedzu2 - zcb*dedxu2
     &                      + zec*dedxv2 - xec*dedzv2
               dedzid2 = ycb*dedxu2 - xcb*dedyu2
     &                      + xec*dedyv2 - yec*dedxv2
               dedxie2 = zdc*dedyv2 - ydc*dedzv2
               dedyie2 = xdc*dedzv2 - zdc*dedxv2
               dedzie2 = ydc*dedxv2 - xdc*dedyv2
c
c     increment the torsion-torsion energy and gradient
c
               ett = ett + e
               dett(1,ia) = dett(1,ia) + dedxia
               dett(2,ia) = dett(2,ia) + dedyia
               dett(3,ia) = dett(3,ia) + dedzia
               dett(1,ib) = dett(1,ib) + dedxib + dedxib2
               dett(2,ib) = dett(2,ib) + dedyib + dedyib2
               dett(3,ib) = dett(3,ib) + dedzib + dedzib2
               dett(1,ic) = dett(1,ic) + dedxic + dedxic2
               dett(2,ic) = dett(2,ic) + dedyic + dedyic2
               dett(3,ic) = dett(3,ic) + dedzic + dedzic2
               dett(1,id) = dett(1,id) + dedxid + dedxid2
               dett(2,id) = dett(2,id) + dedyid + dedyid2
               dett(3,id) = dett(3,id) + dedzid + dedzid2
               dett(1,ie) = dett(1,ie) + dedxie2
               dett(2,ie) = dett(2,ie) + dedyie2
               dett(3,ie) = dett(3,ie) + dedzie2
c
c     increment the internal virial tensor components
c
               vxx = xcb*(dedxic+dedxid) - xba*dedxia + xdc*dedxid
               vyx = ycb*(dedxic+dedxid) - yba*dedxia + ydc*dedxid
               vzx = zcb*(dedxic+dedxid) - zba*dedxia + zdc*dedxid
               vyy = ycb*(dedyic+dedyid) - yba*dedyia + ydc*dedyid
               vzy = zcb*(dedyic+dedyid) - zba*dedyia + zdc*dedyid
               vzz = zcb*(dedzic+dedzid) - zba*dedzia + zdc*dedzid
               vxx2 = xdc*(dedxid2+dedxie2) - xcb*dedxib2 + xed*dedxie2
               vyx2 = ydc*(dedxid2+dedxie2) - ycb*dedxib2 + yed*dedxie2
               vzx2 = zdc*(dedxid2+dedxie2) - zcb*dedxib2 + zed*dedxie2
               vyy2 = ydc*(dedyid2+dedyie2) - ycb*dedyib2 + yed*dedyie2
               vzy2 = zdc*(dedyid2+dedyie2) - zcb*dedyib2 + zed*dedyie2
               vzz2 = zdc*(dedzid2+dedzie2) - zcb*dedzib2 + zed*dedzie2
               vir(1,1) = vir(1,1) + vxx + vxx2
               vir(2,1) = vir(2,1) + vyx + vyx2
               vir(3,1) = vir(3,1) + vzx + vzx2
               vir(1,2) = vir(1,2) + vyx + vyx2
               vir(2,2) = vir(2,2) + vyy + vyy2
               vir(3,2) = vir(3,2) + vzy + vzy2
               vir(1,3) = vir(1,3) + vzx + vzx2
               vir(2,3) = vir(2,3) + vzy + vzy2
               vir(3,3) = vir(3,3) + vzz + vzz2
            end if
         end if
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end
