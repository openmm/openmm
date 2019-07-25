c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine eangang1  --  angle-angle energy & derivatives  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "eangang1" calculates the angle-angle potential energy and
c     first derivatives with respect to Cartesian coordinates
c
c
      subroutine eangang1
      use sizes
      use angang
      use angbnd
      use angpot
      use atoms
      use bound
      use deriv
      use energi
      use group
      use math
      use usage
      use virial
      implicit none
      integer i,k,iangang
      integer ia,ib,ic,id,ie
      real*8 e,angle,fgrp
      real*8 dot,cosine
      real*8 dt1,deddt1
      real*8 dt2,deddt2
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xie,yie,zie
      real*8 xab,yab,zab
      real*8 xcb,ycb,zcb
      real*8 xdb,ydb,zdb
      real*8 xeb,yeb,zeb
      real*8 rab2,rcb2
      real*8 rdb2,reb2
      real*8 xp,yp,zp,rp
      real*8 xq,yq,zq,rq
      real*8 terma,termc
      real*8 termd,terme
      real*8 dedxia,dedyia,dedzia
      real*8 dedxib,dedyib,dedzib
      real*8 dedxic,dedyic,dedzic
      real*8 dedxid,dedyid,dedzid
      real*8 dedxie,dedyie,dedzie
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      logical proceed
c
c
c     zero out the angle-angle energy and first derivatives
c
      eaa = 0.0d0
      do i = 1, n
         deaa(1,i) = 0.0d0
         deaa(2,i) = 0.0d0
         deaa(3,i) = 0.0d0
      end do
      if (nangang .eq. 0)  return
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(nangang,iaa,iang,
!$OMP& use,x,y,z,anat,kaa,aaunit,use_group,use_polymer)
!$OMP& shared(eaa,deaa,vir)
!$OMP DO reduction(+:eaa,deaa,vir) schedule(guided)
c
c     find the energy of each angle-angle interaction
c
      do iangang = 1, nangang
         i = iaa(1,iangang)
         k = iaa(2,iangang)
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         id = iang(1,k)
         ie = iang(3,k)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,ia,ib,ic,id,ie,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or. use(ic)
     &                               .or. use(id) .or. use(ie))
c
c     get the coordinates of the atoms in the angle
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
c
c     compute the values of the two bond angles
c
            xab = xia - xib
            yab = yia - yib
            zab = zia - zib
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdb = xid - xib
            ydb = yid - yib
            zdb = zid - zib
            xeb = xie - xib
            yeb = yie - yib
            zeb = zie - zib
            if (use_polymer) then
               call image (xab,yab,zab)
               call image (xcb,ycb,zcb)
               call image (xdb,ydb,zdb)
               call image (xeb,yeb,zeb)
            end if
            rab2 = xab*xab + yab*yab + zab*zab
            rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
            rdb2 = xdb*xdb + ydb*ydb + zdb*zdb
            reb2 = xeb*xeb + yeb*yeb + zeb*zeb
            xp = ycb*zab - zcb*yab
            yp = zcb*xab - xcb*zab
            zp = xcb*yab - ycb*xab
            xq = yeb*zdb - zeb*ydb
            yq = zeb*xdb - xeb*zdb
            zq = xeb*ydb - yeb*xdb
            rp = sqrt(xp*xp + yp*yp + zp*zp)
            rq = sqrt(xq*xq + yq*yq + zq*zq)
            if (rp*rq .ne. 0.0d0) then
               dot = xab*xcb + yab*ycb + zab*zcb
               cosine = dot / sqrt(rab2*rcb2)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
               dt1 = angle - anat(i)
               dot = xdb*xeb + ydb*yeb + zdb*zeb
               cosine = dot / sqrt(rdb2*reb2)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
               dt2 = angle - anat(k)
c
c     get the energy and master chain rule terms for derivatives
c
               e = aaunit * kaa(iangang) * dt1 * dt2
               deddt1 = radian * e / dt1
               deddt2 = radian * e / dt2
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  e = e * fgrp
                  deddt1 = deddt1 * fgrp
                  deddt2 = deddt2 * fgrp
               end if
c
c     find chain rule terms for the first bond angle deviation
c
               terma = -deddt1 / (rab2*rp)
               termc = deddt1 / (rcb2*rp)
               dedxia = terma * (yab*zp-zab*yp)
               dedyia = terma * (zab*xp-xab*zp)
               dedzia = terma * (xab*yp-yab*xp)
               dedxic = termc * (ycb*zp-zcb*yp)
               dedyic = termc * (zcb*xp-xcb*zp)
               dedzic = termc * (xcb*yp-ycb*xp)
c
c     find chain rule terms for the second bond angle deviation
c
               termd = -deddt2 / (rdb2*rq)
               terme = deddt2 / (reb2*rq)
               dedxid = termd * (ydb*zq-zdb*yq)
               dedyid = termd * (zdb*xq-xdb*zq)
               dedzid = termd * (xdb*yq-ydb*xq)
               dedxie = terme * (yeb*zq-zeb*yq)
               dedyie = terme * (zeb*xq-xeb*zq)
               dedzie = terme * (xeb*yq-yeb*xq)
c
c     get the central atom derivative terms by difference
c
               dedxib = -dedxia - dedxic - dedxid - dedxie
               dedyib = -dedyia - dedyic - dedyid - dedyie
               dedzib = -dedzia - dedzic - dedzid - dedzie
c
c     increment the total angle-angle energy and derivatives
c
               eaa = eaa + e
               deaa(1,ia) = deaa(1,ia) + dedxia
               deaa(2,ia) = deaa(2,ia) + dedyia
               deaa(3,ia) = deaa(3,ia) + dedzia
               deaa(1,ib) = deaa(1,ib) + dedxib
               deaa(2,ib) = deaa(2,ib) + dedyib
               deaa(3,ib) = deaa(3,ib) + dedzib
               deaa(1,ic) = deaa(1,ic) + dedxic
               deaa(2,ic) = deaa(2,ic) + dedyic
               deaa(3,ic) = deaa(3,ic) + dedzic
               deaa(1,id) = deaa(1,id) + dedxid
               deaa(2,id) = deaa(2,id) + dedyid
               deaa(3,id) = deaa(3,id) + dedzid
               deaa(1,ie) = deaa(1,ie) + dedxie
               deaa(2,ie) = deaa(2,ie) + dedyie
               deaa(3,ie) = deaa(3,ie) + dedzie
c
c     increment the internal virial tensor components
c
               vxx = xab*dedxia + xcb*dedxic + xdb*dedxid + xeb*dedxie
               vyx = yab*dedxia + ycb*dedxic + ydb*dedxid + yeb*dedxie
               vzx = zab*dedxia + zcb*dedxic + zdb*dedxid + zeb*dedxie
               vyy = yab*dedyia + ycb*dedyic + ydb*dedyid + yeb*dedyie
               vzy = zab*dedyia + zcb*dedyic + zdb*dedyid + zeb*dedyie
               vzz = zab*dedzia + zcb*dedzic + zdb*dedzid + zeb*dedzie
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
         end if
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end
