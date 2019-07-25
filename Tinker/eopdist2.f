c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine eopdist2  --  atomwise out-plane dist Hessian  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "eopdist2" calculates second derivatives of the out-of-plane
c     distance energy for a single atom via the central atom height
c
c
      subroutine eopdist2 (i)
      use sizes
      use angpot
      use atoms
      use bound
      use group
      use hessn
      use opdist
      use usage
      implicit none
      integer i,kopdist
      integer ia,ib,ic,id
      real*8 force,fgrp
      real*8 dt,dt2,dt3,dt4
      real*8 deddt,drt2,ddrt
      real*8 dot,dotx,doty,dotz
      real*8 term,termd
      real*8 termx,termy,termz
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xad,yad,zad
      real*8 xbd,ybd,zbd
      real*8 xcd,ycd,zcd
      real*8 xt,yt,zt,rt2
      real*8 xtd,ytd,ztd
      real*8 xyt,xzt,yxt
      real*8 yzt,zxt,zyt
      real*8 dxiaxia,dyiayia,dziazia
      real*8 dxibxib,dyibyib,dzibzib
      real*8 dxicxic,dyicyic,dziczic
      real*8 dxidxid,dyidyid,dzidzid
      real*8 dxiayia,dxiazia,dyiazia
      real*8 dxibyib,dxibzib,dyibzib
      real*8 dxicyic,dxiczic,dyiczic
      real*8 dxidyid,dxidzid,dyidzid
      real*8 dxiaxib,dxiayib,dxiazib
      real*8 dyiaxib,dyiayib,dyiazib
      real*8 dziaxib,dziayib,dziazib
      real*8 dxiaxic,dxiayic,dxiazic
      real*8 dyiaxic,dyiayic,dyiazic
      real*8 dziaxic,dziayic,dziazic
      real*8 dxiaxid,dxiayid,dxiazid
      real*8 dyiaxid,dyiayid,dyiazid
      real*8 dziaxid,dziayid,dziazid
      real*8 dxibxic,dxibyic,dxibzic
      real*8 dyibxic,dyibyic,dyibzic
      real*8 dzibxic,dzibyic,dzibzic
      real*8 dxibxid,dxibyid,dxibzid
      real*8 dyibxid,dyibyid,dyibzid
      real*8 dzibxid,dzibyid,dzibzid
      real*8 dxicxid,dxicyid,dxiczid
      real*8 dyicxid,dyicyid,dyiczid
      real*8 dzicxid,dzicyid,dziczid
      logical proceed
c
c
c     compute the out-of-plane distance term Hessian elements
c
      do kopdist = 1, nopdist
         ia = iopd(1,kopdist)
         ib = iopd(2,kopdist)
         ic = iopd(3,kopdist)
         id = iopd(4,kopdist)
         force = opdk(kopdist)
c
c     decide whether to compute the current interaction
c
         proceed = (i.eq.ia .or. i.eq.ib .or. i.eq.ic .or. i.eq.id)
         if (proceed .and. use_group)
     &      call groups (proceed,fgrp,ia,ib,ic,id,0,0)
c
c     get the coordinates of the defining atoms
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
c
c     compute the out-of-plane distance for central atom
c
            xad = xia - xid
            yad = yia - yid
            zad = zia - zid
            xbd = xib - xid
            ybd = yib - yid
            zbd = zib - zid
            xcd = xic - xid
            ycd = yic - yid
            zcd = zic - zid
            if (use_polymer) then
               call image (xad,yad,zad)
               call image (xbd,ybd,zbd)
               call image (xcd,ycd,zcd)
            end if
            xt = ybd*zcd - zbd*ycd
            yt = zbd*xcd - xbd*zcd
            zt = xbd*ycd - ybd*xcd
            rt2 = xt*xt + yt*yt + zt*zt
            dot = xt*xad + yt*yad + zt*zad
            drt2 = dot / rt2
            dt2 = dot * drt2
            dt = sqrt(dt2)
            dt3 = dt2 * dt
            dt4 = dt2 * dt2
            deddt = opdunit * force
     &                 * (2.0d0 + 3.0d0*copd*dt + 4.0d0*qopd*dt2
     &                     + 5.0d0*popd*dt3 + 6.0d0*sopd*dt4)
c
c     scale the interaction based on its group membership
c
            if (use_group) then
               deddt = deddt * fgrp
            end if
c
c     abbreviations for second derivative chain rule terms
c
            term = deddt / rt2
            termx = term * xt
            termy = term * yt
            termz = term * zt
            termd = term * dot
            xtd = xad - 2.0d0*xt*drt2
            ytd = yad - 2.0d0*yt*drt2
            ztd = zad - 2.0d0*zt*drt2
            xyt = xcd*ytd - ycd*xtd
            xzt = xbd*ztd - zbd*xtd
            yxt = ybd*xtd - xbd*ytd
            yzt = ycd*ztd - zcd*ytd
            zxt = zcd*xtd - xcd*ztd
            zyt = zbd*ytd - ybd*ztd
            ddrt = dot * drt2
            dotx = dot * xad
            doty = dot * yad
            dotz = dot * zad
c
c     chain rule terms for second derivative components
c
            dxiaxia = termx*xt
            dxiayia = termx*yt
            dxiazia = termx*zt
            dxiaxib = termx*yzt
            dxiayib = termx*zxt + termd*zcd
            dxiazib = termx*xyt - termd*ycd
            dxiaxic = termx*zyt
            dxiayic = termx*xzt - termd*zbd
            dxiazic = termx*yxt + termd*ybd
            dyiayia = termy*yt
            dyiazia = termy*zt
            dyiaxib = termy*yzt - termd*zcd
            dyiayib = termy*zxt
            dyiazib = termy*xyt + termd*xcd
            dyiaxic = termy*zyt + termd*zbd
            dyiayic = termy*xzt
            dyiazic = termy*yxt - termd*xbd
            dziazia = termz*zt
            dziaxib = termz*yzt + termd*ycd
            dziayib = termz*zxt - termd*xcd
            dziazib = termz*xyt
            dziaxic = termz*zyt - termd*ybd
            dziayic = termz*xzt + termd*xbd
            dziazic = termz*yxt
            dxibxib = term * (yzt*yzt - ddrt*(ycd*ycd+zcd*zcd))
            dxibyib = term * (yzt*zxt + ddrt*xcd*ycd)
            dxibzib = term * (yzt*xyt + ddrt*xcd*zcd)
            dxibxic = term * (yzt*zyt + ddrt*(ybd*ycd+zbd*zcd))
            dxibyic = term * (xzt*yzt - ddrt*(xbd*ycd+zt) + dotz)
            dxibzic = term * (yxt*yzt - ddrt*(xbd*zcd-yt) - doty)
            dyibyib = term * (zxt*zxt - ddrt*(xcd*xcd+zcd*zcd))
            dyibzib = term * (zxt*xyt + ddrt*ycd*zcd)
            dyibxic = term * (zyt*zxt - ddrt*(ybd*xcd-zt) - dotz)
            dyibyic = term * (zxt*xzt + ddrt*(xbd*xcd+zbd*zcd))
            dyibzic = term * (yxt*zxt - ddrt*(ybd*zcd+xt) + dotx)
            dzibzib = term * (xyt*xyt - ddrt*(xcd*xcd+ycd*ycd))
            dzibxic = term * (zyt*xyt - ddrt*(zbd*xcd+yt) + doty)
            dzibyic = term * (xzt*xyt - ddrt*(zbd*ycd-xt) - dotx)
            dzibzic = term * (xyt*yxt + ddrt*(xbd*xcd+ybd*ycd))
            dxicxic = term * (zyt*zyt - ddrt*(ybd*ybd+zbd*zbd))
            dxicyic = term * (zyt*xzt + ddrt*xbd*ybd)
            dxiczic = term * (zyt*yxt + ddrt*xbd*zbd)
            dyicyic = term * (xzt*xzt - ddrt*(xbd*xbd+zbd*zbd))
            dyiczic = term * (xzt*yxt + ddrt*ybd*zbd)
            dziczic = term * (yxt*yxt - ddrt*(xbd*xbd+ybd*ybd))
c
c     get some second derivative chain rule terms by difference
c
            dxiaxid = -dxiaxia - dxiaxib - dxiaxic
            dxiayid = -dxiayia - dxiayib - dxiayic
            dxiazid = -dxiazia - dxiazib - dxiazic
            dyiaxid = -dxiayia - dyiaxib - dyiaxic
            dyiayid = -dyiayia - dyiayib - dyiayic
            dyiazid = -dyiazia - dyiazib - dyiazic
            dziaxid = -dxiazia - dziaxib - dziaxic
            dziayid = -dyiazia - dziayib - dziayic
            dziazid = -dziazia - dziazib - dziazic
            dxibxid = -dxiaxib - dxibxib - dxibxic
            dxibyid = -dyiaxib - dxibyib - dxibyic
            dxibzid = -dziaxib - dxibzib - dxibzic
            dyibxid = -dxiayib - dxibyib - dyibxic
            dyibyid = -dyiayib - dyibyib - dyibyic
            dyibzid = -dziayib - dyibzib - dyibzic
            dzibxid = -dxiazib - dxibzib - dzibxic
            dzibyid = -dyiazib - dyibzib - dzibyic
            dzibzid = -dziazib - dzibzib - dzibzic
            dxicxid = -dxiaxic - dxibxic - dxicxic
            dxicyid = -dyiaxic - dyibxic - dxicyic
            dxiczid = -dziaxic - dzibxic - dxiczic
            dyicxid = -dxiayic - dxibyic - dxicyic
            dyicyid = -dyiayic - dyibyic - dyicyic
            dyiczid = -dziayic - dzibyic - dyiczic
            dzicxid = -dxiazic - dxibzic - dxiczic
            dzicyid = -dyiazic - dyibzic - dyiczic
            dziczid = -dziazic - dzibzic - dziczic
            dxidxid = -dxiaxid - dxibxid - dxicxid
            dxidyid = -dxiayid - dxibyid - dxicyid
            dxidzid = -dxiazid - dxibzid - dxiczid
            dyidyid = -dyiayid - dyibyid - dyicyid
            dyidzid = -dyiazid - dyibzid - dyiczid
            dzidzid = -dziazid - dzibzid - dziczid
c
c     increment diagonal and off-diagonal Hessian elements
c
            if (i .eq. ia) then
               hessx(1,ia) = hessx(1,ia) + dxiaxia
               hessy(1,ia) = hessy(1,ia) + dxiayia
               hessz(1,ia) = hessz(1,ia) + dxiazia
               hessx(2,ia) = hessx(2,ia) + dxiayia
               hessy(2,ia) = hessy(2,ia) + dyiayia
               hessz(2,ia) = hessz(2,ia) + dyiazia
               hessx(3,ia) = hessx(3,ia) + dxiazia
               hessy(3,ia) = hessy(3,ia) + dyiazia
               hessz(3,ia) = hessz(3,ia) + dziazia
               hessx(1,ib) = hessx(1,ib) + dxiaxib
               hessy(1,ib) = hessy(1,ib) + dyiaxib
               hessz(1,ib) = hessz(1,ib) + dziaxib
               hessx(2,ib) = hessx(2,ib) + dxiayib
               hessy(2,ib) = hessy(2,ib) + dyiayib
               hessz(2,ib) = hessz(2,ib) + dziayib
               hessx(3,ib) = hessx(3,ib) + dxiazib
               hessy(3,ib) = hessy(3,ib) + dyiazib
               hessz(3,ib) = hessz(3,ib) + dziazib
               hessx(1,ic) = hessx(1,ic) + dxiaxic
               hessy(1,ic) = hessy(1,ic) + dyiaxic
               hessz(1,ic) = hessz(1,ic) + dziaxic
               hessx(2,ic) = hessx(2,ic) + dxiayic
               hessy(2,ic) = hessy(2,ic) + dyiayic
               hessz(2,ic) = hessz(2,ic) + dziayic
               hessx(3,ic) = hessx(3,ic) + dxiazic
               hessy(3,ic) = hessy(3,ic) + dyiazic
               hessz(3,ic) = hessz(3,ic) + dziazic
               hessx(1,id) = hessx(1,id) + dxiaxid
               hessy(1,id) = hessy(1,id) + dyiaxid
               hessz(1,id) = hessz(1,id) + dziaxid
               hessx(2,id) = hessx(2,id) + dxiayid
               hessy(2,id) = hessy(2,id) + dyiayid
               hessz(2,id) = hessz(2,id) + dziayid
               hessx(3,id) = hessx(3,id) + dxiazid
               hessy(3,id) = hessy(3,id) + dyiazid
               hessz(3,id) = hessz(3,id) + dziazid
            else if (i .eq. ib) then
               hessx(1,ib) = hessx(1,ib) + dxibxib
               hessy(1,ib) = hessy(1,ib) + dxibyib
               hessz(1,ib) = hessz(1,ib) + dxibzib
               hessx(2,ib) = hessx(2,ib) + dxibyib
               hessy(2,ib) = hessy(2,ib) + dyibyib
               hessz(2,ib) = hessz(2,ib) + dyibzib
               hessx(3,ib) = hessx(3,ib) + dxibzib
               hessy(3,ib) = hessy(3,ib) + dyibzib
               hessz(3,ib) = hessz(3,ib) + dzibzib
               hessx(1,ia) = hessx(1,ia) + dxiaxib
               hessy(1,ia) = hessy(1,ia) + dxiayib
               hessz(1,ia) = hessz(1,ia) + dxiazib
               hessx(2,ia) = hessx(2,ia) + dyiaxib
               hessy(2,ia) = hessy(2,ia) + dyiayib
               hessz(2,ia) = hessz(2,ia) + dyiazib
               hessx(3,ia) = hessx(3,ia) + dziaxib
               hessy(3,ia) = hessy(3,ia) + dziayib
               hessz(3,ia) = hessz(3,ia) + dziazib
               hessx(1,ic) = hessx(1,ic) + dxibxic
               hessy(1,ic) = hessy(1,ic) + dyibxic
               hessz(1,ic) = hessz(1,ic) + dzibxic
               hessx(2,ic) = hessx(2,ic) + dxibyic
               hessy(2,ic) = hessy(2,ic) + dyibyic
               hessz(2,ic) = hessz(2,ic) + dzibyic
               hessx(3,ic) = hessx(3,ic) + dxibzic
               hessy(3,ic) = hessy(3,ic) + dyibzic
               hessz(3,ic) = hessz(3,ic) + dzibzic
               hessx(1,id) = hessx(1,id) + dxibxid
               hessy(1,id) = hessy(1,id) + dyibxid
               hessz(1,id) = hessz(1,id) + dzibxid
               hessx(2,id) = hessx(2,id) + dxibyid
               hessy(2,id) = hessy(2,id) + dyibyid
               hessz(2,id) = hessz(2,id) + dzibyid
               hessx(3,id) = hessx(3,id) + dxibzid
               hessy(3,id) = hessy(3,id) + dyibzid
               hessz(3,id) = hessz(3,id) + dzibzid
            else if (i .eq. ic) then
               hessx(1,ic) = hessx(1,ic) + dxicxic
               hessy(1,ic) = hessy(1,ic) + dxicyic
               hessz(1,ic) = hessz(1,ic) + dxiczic
               hessx(2,ic) = hessx(2,ic) + dxicyic
               hessy(2,ic) = hessy(2,ic) + dyicyic
               hessz(2,ic) = hessz(2,ic) + dyiczic
               hessx(3,ic) = hessx(3,ic) + dxiczic
               hessy(3,ic) = hessy(3,ic) + dyiczic
               hessz(3,ic) = hessz(3,ic) + dziczic
               hessx(1,ia) = hessx(1,ia) + dxiaxic
               hessy(1,ia) = hessy(1,ia) + dxiayic
               hessz(1,ia) = hessz(1,ia) + dxiazic
               hessx(2,ia) = hessx(2,ia) + dyiaxic
               hessy(2,ia) = hessy(2,ia) + dyiayic
               hessz(2,ia) = hessz(2,ia) + dyiazic
               hessx(3,ia) = hessx(3,ia) + dziaxic
               hessy(3,ia) = hessy(3,ia) + dziayic
               hessz(3,ia) = hessz(3,ia) + dziazic
               hessx(1,ib) = hessx(1,ib) + dxibxic
               hessy(1,ib) = hessy(1,ib) + dxibyic
               hessz(1,ib) = hessz(1,ib) + dxibzic
               hessx(2,ib) = hessx(2,ib) + dyibxic
               hessy(2,ib) = hessy(2,ib) + dyibyic
               hessz(2,ib) = hessz(2,ib) + dyibzic
               hessx(3,ib) = hessx(3,ib) + dzibxic
               hessy(3,ib) = hessy(3,ib) + dzibyic
               hessz(3,ib) = hessz(3,ib) + dzibzic
               hessx(1,id) = hessx(1,id) + dxicxid
               hessy(1,id) = hessy(1,id) + dyicxid
               hessz(1,id) = hessz(1,id) + dzicxid
               hessx(2,id) = hessx(2,id) + dxicyid
               hessy(2,id) = hessy(2,id) + dyicyid
               hessz(2,id) = hessz(2,id) + dzicyid
               hessx(3,id) = hessx(3,id) + dxiczid
               hessy(3,id) = hessy(3,id) + dyiczid
               hessz(3,id) = hessz(3,id) + dziczid
            else if (i .eq. id) then
               hessx(1,id) = hessx(1,id) + dxidxid
               hessy(1,id) = hessy(1,id) + dxidyid
               hessz(1,id) = hessz(1,id) + dxidzid
               hessx(2,id) = hessx(2,id) + dxidyid
               hessy(2,id) = hessy(2,id) + dyidyid
               hessz(2,id) = hessz(2,id) + dyidzid
               hessx(3,id) = hessx(3,id) + dxidzid
               hessy(3,id) = hessy(3,id) + dyidzid
               hessz(3,id) = hessz(3,id) + dzidzid
               hessx(1,ia) = hessx(1,ia) + dxiaxid
               hessy(1,ia) = hessy(1,ia) + dxiayid
               hessz(1,ia) = hessz(1,ia) + dxiazid
               hessx(2,ia) = hessx(2,ia) + dyiaxid
               hessy(2,ia) = hessy(2,ia) + dyiayid
               hessz(2,ia) = hessz(2,ia) + dyiazid
               hessx(3,ia) = hessx(3,ia) + dziaxid
               hessy(3,ia) = hessy(3,ia) + dziayid
               hessz(3,ia) = hessz(3,ia) + dziazid
               hessx(1,ib) = hessx(1,ib) + dxibxid
               hessy(1,ib) = hessy(1,ib) + dxibyid
               hessz(1,ib) = hessz(1,ib) + dxibzid
               hessx(2,ib) = hessx(2,ib) + dyibxid
               hessy(2,ib) = hessy(2,ib) + dyibyid
               hessz(2,ib) = hessz(2,ib) + dyibzid
               hessx(3,ib) = hessx(3,ib) + dzibxid
               hessy(3,ib) = hessy(3,ib) + dzibyid
               hessz(3,ib) = hessz(3,ib) + dzibzid
               hessx(1,ic) = hessx(1,ic) + dxicxid
               hessy(1,ic) = hessy(1,ic) + dxicyid
               hessz(1,ic) = hessz(1,ic) + dxiczid
               hessx(2,ic) = hessx(2,ic) + dyicxid
               hessy(2,ic) = hessy(2,ic) + dyicyid
               hessz(2,ic) = hessz(2,ic) + dyiczid
               hessx(3,ic) = hessx(3,ic) + dzicxid
               hessy(3,ic) = hessy(3,ic) + dzicyid
               hessz(3,ic) = hessz(3,ic) + dziczid
            end if
         end if
      end do
      return
      end
