c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine eopbend2  --  out-of-plane bend Hessian; numer  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "eopbend2" calculates second derivatives of the out-of-plane
c     bend energy via a Wilson-Decius-Cross or Allinger angle for
c     a single atom using finite difference methods
c
c
      subroutine eopbend2 (i)
      use sizes
      use angbnd
      use atoms
      use group
      use hessn
      use opbend
      implicit none
      integer i,j,k,iopbend
      integer ia,ib,ic,id
      real*8 eps,fgrp
      real*8 old,term
      real*8, allocatable :: de(:,:)
      real*8, allocatable :: d0(:,:)
      logical proceed
      logical twosided
c
c
c     set stepsize for derivatives and default group weight
c
      eps = 1.0d-5
      fgrp = 1.0d0
      twosided = .false.
      if (n .le. 50)  twosided = .true.
c
c     perform dynamic allocation of some local arrays
c
      allocate (de(3,n))
      allocate (d0(3,n))
c
c     compute numerical out-of-plane Hessian for current atom
c
      do iopbend = 1, nopbend
         k = iopb(iopbend)
         ia = iang(1,k)
         ib = iang(2,k)
         ic = iang(3,k)
         id = iang(4,k)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,ia,ib,ic,id,0,0)
         if (proceed)  proceed = (i.eq.ia .or. i.eq.ib .or.
     &                              i.eq.ic .or. i.eq.id)
         if (proceed) then
            term = fgrp / eps
c
c     find first derivatives for the base structure
c
            if (.not. twosided) then
               call eopbend2a (iopbend,de)
               do j = 1, 3
                  d0(j,ia) = de(j,ia)
                  d0(j,ib) = de(j,ib)
                  d0(j,ic) = de(j,ic)
                  d0(j,id) = de(j,id)
               end do
            end if
c
c     find numerical x-components via perturbed structures
c
            old = x(i)
            if (twosided) then
               x(i) = x(i) - 0.5d0*eps
               call eopbend2a (iopbend,de)
               do j = 1, 3
                  d0(j,ia) = de(j,ia)
                  d0(j,ib) = de(j,ib)
                  d0(j,ic) = de(j,ic)
                  d0(j,id) = de(j,id)
               end do
            end if
            x(i) = x(i) + eps
            call eopbend2a (iopbend,de)
            x(i) = old
            do j = 1, 3
               hessx(j,ia) = hessx(j,ia) + term*(de(j,ia)-d0(j,ia))
               hessx(j,ib) = hessx(j,ib) + term*(de(j,ib)-d0(j,ib))
               hessx(j,ic) = hessx(j,ic) + term*(de(j,ic)-d0(j,ic))
               hessx(j,id) = hessx(j,id) + term*(de(j,id)-d0(j,id))
            end do
c
c     find numerical y-components via perturbed structures
c
            old = y(i)
            if (twosided) then
               y(i) = y(i) - 0.5d0*eps
               call eopbend2a (iopbend,de)
               do j = 1, 3
                  d0(j,ia) = de(j,ia)
                  d0(j,ib) = de(j,ib)
                  d0(j,ic) = de(j,ic)
                  d0(j,id) = de(j,id)
               end do
            end if
            y(i) = y(i) + eps
            call eopbend2a (iopbend,de)
            y(i) = old
            do j = 1, 3
               hessy(j,ia) = hessy(j,ia) + term*(de(j,ia)-d0(j,ia))
               hessy(j,ib) = hessy(j,ib) + term*(de(j,ib)-d0(j,ib))
               hessy(j,ic) = hessy(j,ic) + term*(de(j,ic)-d0(j,ic))
               hessy(j,id) = hessy(j,id) + term*(de(j,id)-d0(j,id))
            end do
c
c     find numerical z-components via perturbed structures
c
            old = z(i)
            if (twosided) then
               z(i) = z(i) - 0.5d0*eps
               call eopbend2a (iopbend,de)
               do j = 1, 3
                  d0(j,ia) = de(j,ia)
                  d0(j,ib) = de(j,ib)
                  d0(j,ic) = de(j,ic)
                  d0(j,id) = de(j,id)
               end do
            end if
            z(i) = z(i) + eps
            call eopbend2a (iopbend,de)
            z(i) = old
            do j = 1, 3
               hessz(j,ia) = hessz(j,ia) + term*(de(j,ia)-d0(j,ia))
               hessz(j,ib) = hessz(j,ib) + term*(de(j,ib)-d0(j,ib))
               hessz(j,ic) = hessz(j,ic) + term*(de(j,ic)-d0(j,ic))
               hessz(j,id) = hessz(j,id) + term*(de(j,id)-d0(j,id))
            end do
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (de)
      deallocate (d0)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine eopbend2a  --  out-of-plane bend derivatives  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "eopbend2a" calculates out-of-plane bend first derivatives at
c     a trigonal center via a Wilson-Decius-Cross or Allinger angle;
c     used in computation of finite difference second derivatives
c
c
      subroutine eopbend2a (i,de)
      use sizes
      use angbnd
      use angpot
      use atoms
      use bound
      use math
      use opbend
      implicit none
      integer i,k
      integer ia,ib,ic,id
      real*8 angle,force
      real*8 dot,cosine
      real*8 cc,ee,bkk2,term
      real*8 deddt,dedcos
      real*8 dt,dt2,dt3,dt4
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xab,yab,zab
      real*8 xcb,ycb,zcb
      real*8 xdb,ydb,zdb
      real*8 xad,yad,zad
      real*8 xcd,ycd,zcd
      real*8 rdb2,rad2,rcd2
      real*8 rab2,rcb2
      real*8 dccdxia,dccdyia,dccdzia
      real*8 dccdxic,dccdyic,dccdzic
      real*8 dccdxid,dccdyid,dccdzid
      real*8 deedxia,deedyia,deedzia
      real*8 deedxic,deedyic,deedzic
      real*8 deedxid,deedyid,deedzid
      real*8 dedxia,dedyia,dedzia
      real*8 dedxib,dedyib,dedzib
      real*8 dedxic,dedyic,dedzic
      real*8 dedxid,dedyid,dedzid
      real*8 de(3,*)
c
c
c     set the atom numbers and parameters for this angle
c
      k = iopb(i)
      ia = iang(1,k)
      ib = iang(2,k)
      ic = iang(3,k)
      id = iang(4,k)
      force = opbk(i)
c
c     get the coordinates of the atoms in the angle
c
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
c
c     zero out the first derivative components
c
      de(1,ia) = 0.0d0
      de(2,ia) = 0.0d0
      de(3,ia) = 0.0d0
      de(1,ib) = 0.0d0
      de(2,ib) = 0.0d0
      de(3,ib) = 0.0d0
      de(1,ic) = 0.0d0
      de(2,ic) = 0.0d0
      de(3,ic) = 0.0d0
      de(1,id) = 0.0d0
      de(2,id) = 0.0d0
      de(3,id) = 0.0d0
c
c     compute the out-of-plane bending angle
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
      xad = xia - xid
      yad = yia - yid
      zad = zia - zid
      xcd = xic - xid
      ycd = yic - yid
      zcd = zic - zid
      if (use_polymer) then
         call image (xab,yab,zab)
         call image (xcb,ycb,zcb)
         call image (xdb,ydb,zdb)
         call image (xad,yad,zad)
         call image (xcd,ycd,zcd)
      end if
c
c     W-D-C angle between A-B-C plane and B-D vector for D-B<AC
c
      if (opbtyp .eq. 'W-D-C') then
         rab2 = xab*xab + yab*yab + zab*zab
         rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
         dot = xab*xcb+yab*ycb+zab*zcb
         cc = rab2*rcb2 - dot*dot
c
c     Allinger angle between A-C-D plane and D-B vector for D-B<AC
c
      else if (opbtyp .eq. 'ALLINGER') then
         rad2 = xad*xad + yad*yad + zad*zad
         rcd2 = xcd*xcd + ycd*ycd + zcd*zcd
         dot = xad*xcd + yad*ycd + zad*zcd
         cc = rad2*rcd2 - dot*dot
      end if
c
c     get the out-of-plane bending master chain rule terms
c
      ee = xdb*(yab*zcb-zab*ycb) + ydb*(zab*xcb-xab*zcb)
     &        + zdb*(xab*ycb-yab*xcb)
      rdb2 = xdb*xdb + ydb*ydb + zdb*zdb
      if (rdb2.ne.0.0d0 .and. cc.ne.0.0d0) then
         bkk2 = rdb2 - ee*ee/cc
         cosine = sqrt(bkk2/rdb2)
         cosine = min(1.0d0,max(-1.0d0,cosine))
         angle = radian * acos(cosine)
         dt = angle
         dt2 = dt * dt
         dt3 = dt2 * dt
         dt4 = dt2 * dt2
         deddt = opbunit * force * dt * radian
     &              * (2.0d0 + 3.0d0*copb*dt + 4.0d0*qopb*dt2
     &                  + 5.0d0*popb*dt3 + 6.0d0*sopb*dt4)
         dedcos = -deddt * sign(1.0d0,ee) / sqrt(cc*bkk2)
c
c     chain rule terms for first derivative components
c
         if (opbtyp .eq. 'W-D-C') then
            term = ee / cc
            dccdxia = (xab*rcb2-xcb*dot) * term
            dccdyia = (yab*rcb2-ycb*dot) * term
            dccdzia = (zab*rcb2-zcb*dot) * term
            dccdxic = (xcb*rab2-xab*dot) * term
            dccdyic = (ycb*rab2-yab*dot) * term
            dccdzic = (zcb*rab2-zab*dot) * term
            dccdxid = 0.0d0
            dccdyid = 0.0d0
            dccdzid = 0.0d0
         else if (opbtyp .eq. 'ALLINGER') then
            term = ee / cc
            dccdxia = (xad*rcd2-xcd*dot) * term
            dccdyia = (yad*rcd2-ycd*dot) * term
            dccdzia = (zad*rcd2-zcd*dot) * term
            dccdxic = (xcd*rad2-xad*dot) * term
            dccdyic = (ycd*rad2-yad*dot) * term
            dccdzic = (zcd*rad2-zad*dot) * term
            dccdxid = -dccdxia - dccdxic
            dccdyid = -dccdyia - dccdyic
            dccdzid = -dccdzia - dccdzic
         end if
         term = ee / rdb2
         deedxia = ydb*zcb - zdb*ycb
         deedyia = zdb*xcb - xdb*zcb
         deedzia = xdb*ycb - ydb*xcb
         deedxic = yab*zdb - zab*ydb
         deedyic = zab*xdb - xab*zdb
         deedzic = xab*ydb - yab*xdb
         deedxid = ycb*zab - zcb*yab + xdb*term
         deedyid = zcb*xab - xcb*zab + ydb*term
         deedzid = xcb*yab - ycb*xab + zdb*term
c
c     compute first derivative components for this angle
c
         dedxia = dedcos * (dccdxia+deedxia)
         dedyia = dedcos * (dccdyia+deedyia)
         dedzia = dedcos * (dccdzia+deedzia)
         dedxic = dedcos * (dccdxic+deedxic)
         dedyic = dedcos * (dccdyic+deedyic)
         dedzic = dedcos * (dccdzic+deedzic)
         dedxid = dedcos * (dccdxid+deedxid)
         dedyid = dedcos * (dccdyid+deedyid)
         dedzid = dedcos * (dccdzid+deedzid)
         dedxib = -dedxia - dedxic - dedxid
         dedyib = -dedyia - dedyic - dedyid
         dedzib = -dedzia - dedzic - dedzid
c
c     set the out-of-plane bending derivatives
c
         de(1,ia) = dedxia
         de(2,ia) = dedyia
         de(3,ia) = dedzia
         de(1,ib) = dedxib
         de(2,ib) = dedyib
         de(3,ib) = dedzib
         de(1,ic) = dedxic
         de(2,ic) = dedyic
         de(3,ic) = dedzic
         de(1,id) = dedxid
         de(2,id) = dedyid
         de(3,id) = dedzid
      end if
      return
      end
