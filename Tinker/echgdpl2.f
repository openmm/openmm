c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine echgdpl2  --  atomwise charge-dipole Hessian  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "echgdpl2" calculates second derivatives of the
c     charge-dipole interaction energy for a single atom
c
c
      subroutine echgdpl2 (i)
      use sizes
      use atoms
      use bound
      use cell
      use charge
      use chgpot
      use couple
      use dipole
      use group
      use hessn
      use shunt
      use units
      implicit none
      integer i,k,jcell
      integer ii,i1,k1,k2
      integer, allocatable :: skip(:)
      integer, allocatable :: omit(:)
      real*8 f,fi,fk,fik
      real*8 fgrp,sk1,sk2
      real*8 xi,yi,zi,xk,yk,zk
      real*8 xr,yr,zr,xq,yq,zq
      real*8 e,r2,rk2,rkr3,dotk
      real*8 term,term2,part,part2
      real*8 termx,termy,termz
      real*8 termxk,termyk,termzk
      real*8 xrr2,yrr2,zrr2
      real*8 xkrk2,ykrk2,zkrk2
      real*8 dotk2,dotkr2,dotkrk2
      real*8 factor,factork
      real*8 dedxi1,dedyi1,dedzi1
      real*8 dedxk1,dedyk1,dedzk1
      real*8 dedxk2,dedyk2,dedzk2
      real*8 dtdxi1,dtdyi1,dtdzi1
      real*8 dtdxk1,dtdyk1,dtdzk1
      real*8 dtdxk2,dtdyk2,dtdzk2
      real*8 dtxdxi1,dtxkdxi1,dtxdxk1
      real*8 dtxkdxk1,dtxdxk2,dtxkdxk2
      real*8 dtydxi1,dtykdxi1,dtydxk1
      real*8 dtykdxk1,dtydxk2,dtykdxk2
      real*8 dtzdxi1,dtzkdxi1,dtzdxk1
      real*8 dtzkdxk1,dtzdxk2,dtzkdxk2
      real*8 dtxdyi1,dtxkdyi1,dtxdyk1
      real*8 dtxkdyk1,dtxdyk2,dtxkdyk2
      real*8 dtydyi1,dtykdyi1,dtydyk1
      real*8 dtykdyk1,dtydyk2,dtykdyk2
      real*8 dtzdyi1,dtzkdyi1,dtzdyk1
      real*8 dtzkdyk1,dtzdyk2,dtzkdyk2
      real*8 dtxdzi1,dtxkdzi1,dtxdzk1
      real*8 dtxkdzk1,dtxdzk2,dtxkdzk2
      real*8 dtydzi1,dtykdzi1,dtydzk1
      real*8 dtykdzk1,dtydzk2,dtykdzk2
      real*8 dtzdzi1,dtzkdzi1,dtzdzk1
      real*8 dtzkdzk1,dtzdzk2,dtzkdzk2
      real*8 r,r3,r4,r5
      real*8 taper,dtaper,d2taper
      real*8 dtaperx,dtapery,dtaperz
      real*8 d2taperxx,d2taperyy
      real*8 d2taperzz,d2taperxy
      real*8 d2taperxz,d2taperyz
      logical proceed
      character*6 mode
c
c
c     check for the presence of both charges and dipoles
c
      if (ndipole.eq.0 .or. nion.eq.0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (skip(n))
      allocate (omit(n))
c
c     zero out the lists of atoms to be skipped
c
      do k = 1, n
         skip(k) = 0
         omit(k) = 0
      end do
c
c     set conversion factor and switching function coefficients
c
      f = -electric / (debye * dielec)
      mode = 'CHGDPL'
      call switch (mode)
c
c     first see if the atom of interest carries a charge
c
      do ii = 1, nion
         i1 = iion(ii)
         if (i1 .ne. i)  goto 10
         skip(i1) = i1
         do k = 1, n12(i1)
            skip(i12(k,i1)) = i1
         end do
         xi = x(i1)
         yi = y(i1)
         zi = z(i1)
         fi = f * pchg(ii)
c
c     decide whether to compute the current interaction
c
         do k = 1, ndipole
            k1 = idpl(1,k)
            k2 = idpl(2,k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i1,k1,k2,0,0,0)
            if (proceed)  proceed = (skip(k1).ne.i1 .and.
     &                                 skip(k2).ne.i1)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               sk1 = 1.0d0 - sdpl(k)
               sk2 = sdpl(k)
               xk = x(k1) - x(k2)
               yk = y(k1) - y(k2)
               zk = z(k1) - z(k2)
               xr = xi - x(k1) + xk*sk2
               yr = yi - y(k1) + yk*sk2
               zr = zi - z(k1) + zk*sk2
               call image (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  rk2 = xk*xk + yk*yk + zk*zk
                  rkr3 = sqrt(rk2*r2) * r2
                  dotk = xk*xr + yk*yr + zk*zr
                  fik = -fi * bdpl(k)
c
c     scale the interaction based on its group membership
c
                  if (use_group)  fik = fik * fgrp
c
c     some abbreviations used in various chain rule terms
c
                  xrr2 = xr / r2
                  yrr2 = yr / r2
                  zrr2 = zr / r2
                  xkrk2 = xk / rk2
                  ykrk2 = yk / rk2
                  zkrk2 = zk / rk2
                  dotk2 = 2.0d0 * dotk
                  dotkr2 = dotk / r2
c
c     form the chain rule terms for first derivatives
c
                  term = fik / rkr3
                  term2 = -3.0d0 * dotk
                  termx = term * (xk+xrr2*term2)
                  termy = term * (yk+yrr2*term2)
                  termz = term * (zk+zrr2*term2)
                  termxk = term * (xr-dotk*xkrk2)
                  termyk = term * (yr-dotk*ykrk2)
                  termzk = term * (zr-dotk*zkrk2)
c
c     use energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     e = fik * dotk / rkr3
                     dedxi1 = termx
                     dedyi1 = termy
                     dedzi1 = termz
                     dedxk1 = -sk1*termx + termxk
                     dedyk1 = -sk1*termy + termyk
                     dedzk1 = -sk1*termz + termzk
                     dedxk2 = -sk2*termx - termxk
                     dedyk2 = -sk2*termy - termyk
                     dedzk2 = -sk2*termz - termzk
                     r = sqrt(r2)
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     taper = c5*r5 + c4*r4 + c3*r3
     &                          + c2*r2 + c1*r + c0
                     dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                           + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                     d2taper = 20.0d0*c5*r3 + 12.0d0*c4*r2
     &                            + 6.0d0*c3*r + 2.0d0*c2
                     dtaper = dtaper / r
                     dtaperx = xr * dtaper
                     dtapery = yr * dtaper
                     dtaperz = zr * dtaper
                     d2taper = e * (d2taper-dtaper)
                     dtaper = e * dtaper
                     d2taperxx = xr*xrr2*d2taper + dtaper
                     d2taperxy = xr*yrr2*d2taper
                     d2taperxz = xr*zrr2*d2taper
                     d2taperyy = yr*yrr2*d2taper + dtaper
                     d2taperyz = yr*zrr2*d2taper
                     d2taperzz = zr*zrr2*d2taper + dtaper
                     term = term * taper
                     termx = termx * taper
                     termy = termy * taper
                     termz = termz * taper
                     termxk = termxk * taper
                     termyk = termyk * taper
                     termzk = termzk * taper
                  end if
c
c     chain rule terms for second derivative components
c
                  dtdxi1 = -3.0d0 * xrr2
                  part = xk - dotk2*xrr2
                  factor = -3.0d0 * (dotkr2 + xrr2*part)
                  factork = 1.0d0 - xk*xkrk2
                  dtxdxi1 = dtdxi1*termx + term*factor
                  dtxkdxi1 = dtdxi1*termxk + term*factork
                  factor = -3.0d0 * yrr2 * part
                  factork = -yk * xkrk2
                  dtydxi1 = dtdxi1*termy + term*factor
                  dtykdxi1 = dtdxi1*termyk + term*factork
                  factor = -3.0d0 * zrr2 * part
                  factork = -zk * xkrk2
                  dtzdxi1 = dtdxi1*termz + term*factor
                  dtzkdxi1 = dtdxi1*termzk + term*factork
                  dtdyi1 = -3.0d0 * yrr2
                  part = yk - dotk2*yrr2
                  factor = -3.0d0 * xrr2 * part
                  factork = -xk * ykrk2
                  dtxdyi1 = dtdyi1*termx + term*factor
                  dtxkdyi1 = dtdyi1*termxk + term*factork
                  factor = -3.0d0 * (dotkr2 + yrr2*part)
                  factork = 1.0d0 - yk*ykrk2
                  dtydyi1 = dtdyi1*termy + term*factor
                  dtykdyi1 = dtdyi1*termyk + term*factork
                  factor = -3.0d0 * zrr2 * part
                  factork = -zk * ykrk2
                  dtzdyi1 = dtdyi1*termz + term*factor
                  dtzkdyi1 = dtdyi1*termzk + term*factork
                  dtdzi1 = -3.0d0 * zrr2
                  part = zk - dotk2*zrr2
                  factor = -3.0d0 * xrr2 * part
                  factork = -xk * zkrk2
                  dtxdzi1 = dtdzi1*termx + term*factor
                  dtxkdzi1 = dtdzi1*termxk + term*factork
                  factor = -3.0d0 * yrr2 * part
                  factork = -yk * zkrk2
                  dtydzi1 = dtdzi1*termy + term*factor
                  dtykdzi1 = dtdzi1*termyk + term*factork
                  factor = -3.0d0 * (dotkr2 + zrr2*part)
                  factork = 1.0d0 - zk*zkrk2
                  dtzdzi1 = dtdzi1*termz + term*factor
                  dtzkdzi1 = dtdzi1*termzk + term*factork
c
c     increment diagonal and off-diagonal Hessian elements
c
                  hessx(1,i1) = hessx(1,i1) + dtxdxi1
                  hessx(2,i1) = hessx(2,i1) + dtydxi1
                  hessx(3,i1) = hessx(3,i1) + dtzdxi1
                  hessx(1,k1) = hessx(1,k1) - sk1*dtxdxi1 + dtxkdxi1
                  hessx(2,k1) = hessx(2,k1) - sk1*dtydxi1 + dtykdxi1
                  hessx(3,k1) = hessx(3,k1) - sk1*dtzdxi1 + dtzkdxi1
                  hessx(1,k2) = hessx(1,k2) - sk2*dtxdxi1 - dtxkdxi1
                  hessx(2,k2) = hessx(2,k2) - sk2*dtydxi1 - dtykdxi1
                  hessx(3,k2) = hessx(3,k2) - sk2*dtzdxi1 - dtzkdxi1
                  hessy(1,i1) = hessy(1,i1) + dtxdyi1
                  hessy(2,i1) = hessy(2,i1) + dtydyi1
                  hessy(3,i1) = hessy(3,i1) + dtzdyi1
                  hessy(1,k1) = hessy(1,k1) - sk1*dtxdyi1 + dtxkdyi1
                  hessy(2,k1) = hessy(2,k1) - sk1*dtydyi1 + dtykdyi1
                  hessy(3,k1) = hessy(3,k1) - sk1*dtzdyi1 + dtzkdyi1
                  hessy(1,k2) = hessy(1,k2) - sk2*dtxdyi1 - dtxkdyi1
                  hessy(2,k2) = hessy(2,k2) - sk2*dtydyi1 - dtykdyi1
                  hessy(3,k2) = hessy(3,k2) - sk2*dtzdyi1 - dtzkdyi1
                  hessz(1,i1) = hessz(1,i1) + dtxdzi1
                  hessz(2,i1) = hessz(2,i1) + dtydzi1
                  hessz(3,i1) = hessz(3,i1) + dtzdzi1
                  hessz(1,k1) = hessz(1,k1) - sk1*dtxdzi1 + dtxkdzi1
                  hessz(2,k1) = hessz(2,k1) - sk1*dtydzi1 + dtykdzi1
                  hessz(3,k1) = hessz(3,k1) - sk1*dtzdzi1 + dtzkdzi1
                  hessz(1,k2) = hessz(1,k2) - sk2*dtxdzi1 - dtxkdzi1
                  hessz(2,k2) = hessz(2,k2) - sk2*dtydzi1 - dtykdzi1
                  hessz(3,k2) = hessz(3,k2) - sk2*dtzdzi1 - dtzkdzi1
c
c     more energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     hessx(1,i1) = hessx(1,i1) + dtaperx*dedxi1
     &                             + dtaperx*dedxi1 + d2taperxx
                     hessx(2,i1) = hessx(2,i1) + dtaperx*dedyi1
     &                             + dtapery*dedxi1 + d2taperxy
                     hessx(3,i1) = hessx(3,i1) + dtaperx*dedzi1
     &                             + dtaperz*dedxi1 + d2taperxz
                     hessx(1,k1) = hessx(1,k1) + dtaperx*dedxk1
     &                             - sk1*(dtaperx*dedxi1+d2taperxx)
                     hessx(2,k1) = hessx(2,k1) + dtaperx*dedyk1
     &                             - sk1*(dtapery*dedxi1+d2taperxy)
                     hessx(3,k1) = hessx(3,k1) + dtaperx*dedzk1
     &                             - sk1*(dtaperz*dedxi1+d2taperxz)
                     hessx(1,k2) = hessx(1,k2) + dtaperx*dedxk2
     &                             - sk2*(dtaperx*dedxi1+d2taperxx)
                     hessx(2,k2) = hessx(2,k2) + dtaperx*dedyk2
     &                             - sk2*(dtapery*dedxi1+d2taperxy)
                     hessx(3,k2) = hessx(3,k2) + dtaperx*dedzk2
     &                             - sk2*(dtaperz*dedxi1+d2taperxz)
                     hessy(1,i1) = hessy(1,i1) + dtapery*dedxi1
     &                             + dtaperx*dedyi1 + d2taperxy
                     hessy(2,i1) = hessy(2,i1) + dtapery*dedyi1
     &                             + dtapery*dedyi1 + d2taperyy
                     hessy(3,i1) = hessy(3,i1) + dtapery*dedzi1
     &                             + dtaperz*dedyi1 + d2taperyz
                     hessy(1,k1) = hessy(1,k1) + dtapery*dedxk1
     &                             - sk1*(dtaperx*dedyi1+d2taperxy)
                     hessy(2,k1) = hessy(2,k1) + dtapery*dedyk1
     &                             - sk1*(dtapery*dedyi1+d2taperyy)
                     hessy(3,k1) = hessy(3,k1) + dtapery*dedzk1
     &                             - sk1*(dtaperz*dedyi1+d2taperyz)
                     hessy(1,k2) = hessy(1,k2) + dtapery*dedxk2
     &                             - sk2*(dtaperx*dedyi1+d2taperxy)
                     hessy(2,k2) = hessy(2,k2) + dtapery*dedyk2
     &                             - sk2*(dtapery*dedyi1+d2taperyy)
                     hessy(3,k2) = hessy(3,k2) + dtapery*dedzk2
     &                             - sk2*(dtaperz*dedyi1+d2taperyz)
                     hessz(1,i1) = hessz(1,i1) + dtaperz*dedxi1
     &                             + dtaperx*dedzi1 + d2taperxz
                     hessz(2,i1) = hessz(2,i1) + dtaperz*dedyi1
     &                             + dtapery*dedzi1 + d2taperyz
                     hessz(3,i1) = hessz(3,i1) + dtaperz*dedzi1
     &                             + dtaperz*dedzi1 + d2taperzz
                     hessz(1,k1) = hessz(1,k1) + dtaperz*dedxk1
     &                             - sk1*(dtaperx*dedzi1+d2taperxz)
                     hessz(2,k1) = hessz(2,k1) + dtaperz*dedyk1
     &                             - sk1*(dtapery*dedzi1+d2taperyz)
                     hessz(3,k1) = hessz(3,k1) + dtaperz*dedzk1
     &                             - sk1*(dtaperz*dedzi1+d2taperzz)
                     hessz(1,k2) = hessz(1,k2) + dtaperz*dedxk2
     &                             - sk2*(dtaperx*dedzi1+d2taperxz)
                     hessz(2,k2) = hessz(2,k2) + dtaperz*dedyk2
     &                             - sk2*(dtapery*dedzi1+d2taperyz)
                     hessz(3,k2) = hessz(3,k2) + dtaperz*dedzk2
     &                             - sk2*(dtaperz*dedzi1+d2taperzz)
                  end if
               end if
            end if
         end do
   10    continue
      end do
c
c     see if the atom of interest is part of a dipole
c
      do k = 1, ndipole
         k1 = idpl(1,k)
         k2 = idpl(2,k)
         if (k1.ne.i .and. k2.ne.i)  goto 20
         do ii = 1, n12(k1)
            omit(i12(ii,k1)) = k
         end do
         do ii = 1, n12(k2)
            omit(i12(ii,k2)) = k
         end do
         sk1 = 1.0d0 - sdpl(k)
         sk2 = sdpl(k)
         xk = x(k1) - x(k2)
         yk = y(k1) - y(k2)
         zk = z(k1) - z(k2)
         rk2 = xk*xk + yk*yk + zk*zk
         xq = x(k1) - xk*sk2
         yq = y(k1) - yk*sk2
         zq = z(k1) - zk*sk2
         fk = -f * bdpl(k)
c
c     decide whether to compute the current interaction
c
         do ii = 1, nion
            i1 = iion(ii)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i1,k1,k2,0,0,0)
            if (proceed)  proceed = (omit(i1) .ne. k)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = x(i1) - xq
               yr = y(i1) - yq
               zr = z(i1) - zq
               call image (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  rkr3 = sqrt(rk2*r2) * r2
                  dotk = xk*xr + yk*yr + zk*zr
                  fik = fk * pchg(ii)
c
c     scale the interaction based on its group membership
c
                  if (use_group)  fik = fik * fgrp
c
c     some abbreviations used in various chain rule terms
c
                  xrr2 = xr / r2
                  yrr2 = yr / r2
                  zrr2 = zr / r2
                  xkrk2 = xk / rk2
                  ykrk2 = yk / rk2
                  zkrk2 = zk / rk2
                  dotk2 = 2.0d0 * dotk
                  dotkr2 = dotk / r2
                  dotkrk2 = dotk / rk2
c
c     form the chain rule terms for first derivatives
c
                  term = fik / rkr3
                  term2 = -3.0d0 * dotk
                  termx = term * (xk+xrr2*term2)
                  termy = term * (yk+yrr2*term2)
                  termz = term * (zk+zrr2*term2)
                  termxk = term * (xr-dotk*xkrk2)
                  termyk = term * (yr-dotk*ykrk2)
                  termzk = term * (zr-dotk*zkrk2)
c
c     use energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     e = fik * dotk / rkr3
                     dedxi1 = termx
                     dedyi1 = termy
                     dedzi1 = termz
                     dedxk1 = -sk1*termx + termxk
                     dedyk1 = -sk1*termy + termyk
                     dedzk1 = -sk1*termz + termzk
                     dedxk2 = -sk2*termx - termxk
                     dedyk2 = -sk2*termy - termyk
                     dedzk2 = -sk2*termz - termzk
                     r = sqrt(r2)
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     taper = c5*r5 + c4*r4 + c3*r3
     &                          + c2*r2 + c1*r + c0
                     dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                           + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                     d2taper = 20.0d0*c5*r3 + 12.0d0*c4*r2
     &                            + 6.0d0*c3*r + 2.0d0*c2
                     dtaper = dtaper / r
                     dtaperx = xr * dtaper
                     dtapery = yr * dtaper
                     dtaperz = zr * dtaper
                     d2taper = e * (d2taper-dtaper)
                     dtaper = e * dtaper
                     d2taperxx = xr*xrr2*d2taper + dtaper
                     d2taperxy = xr*yrr2*d2taper
                     d2taperxz = xr*zrr2*d2taper
                     d2taperyy = yr*yrr2*d2taper + dtaper
                     d2taperyz = yr*zrr2*d2taper
                     d2taperzz = zr*zrr2*d2taper + dtaper
                     term = term * taper
                     termx = termx * taper
                     termy = termy * taper
                     termz = termz * taper
                     termxk = termxk * taper
                     termyk = termyk * taper
                     termzk = termzk * taper
                  end if
c
c     chain rule terms for second derivative components
c
                  if (k1 .eq. i) then
                     dtdxk1 = 3.0d0*sk1*xrr2 - xkrk2
                     part = sk1*xk - xr
                     part2 = sk1*dotk2*xrr2 - part
                     factor = 1.0d0 - 3.0d0*xrr2*part2
     &                           + 3.0d0*sk1*dotkr2
                     factork = -sk1 + dotk2*xkrk2*xkrk2
     &                            + xkrk2*part - dotkrk2
                     dtxdxk1 = dtdxk1*termx + term*factor
                     dtxkdxk1 = dtdxk1*termxk + term*factork
                     factor = -3.0d0 * yrr2 * part2
                     factork = dotk2*ykrk2*xkrk2 + ykrk2*part
                     dtydxk1 = dtdxk1*termy + term*factor
                     dtykdxk1 = dtdxk1*termyk + term*factork
                     factor = -3.0d0 * zrr2 * part2
                     factork = dotk2*zkrk2*xkrk2 + zkrk2*part
                     dtzdxk1 = dtdxk1*termz + term*factor
                     dtzkdxk1 = dtdxk1*termzk + term*factork
                     dtdyk1 = 3.0d0*sk1*yrr2 - ykrk2
                     part = sk1*yk - yr
                     part2 = sk1*dotk2*yrr2 - part
                     factor = -3.0d0 * xrr2 * part2
                     factork = dotk2*xkrk2*ykrk2 + xkrk2*part
                     dtxdyk1 = dtdyk1*termx + term*factor
                     dtxkdyk1 = dtdyk1*termxk + term*factork
                     factor = 1.0d0 - 3.0d0*yrr2*part2
     &                           + 3.0d0*sk1*dotkr2
                     factork = -sk1 + dotk2*ykrk2*ykrk2
     &                            + ykrk2*part - dotkrk2
                     dtydyk1 = dtdyk1*termy + term*factor
                     dtykdyk1 = dtdyk1*termyk + term*factork
                     factor = -3.0d0 * zrr2 * part2
                     factork = dotk2*zkrk2*ykrk2 + zkrk2*part
                     dtzdyk1 = dtdyk1*termz + term*factor
                     dtzkdyk1 = dtdyk1*termzk + term*factork
                     dtdzk1 = 3.0d0*sk1*zrr2 - zkrk2
                     part = sk1*zk - zr
                     part2 = sk1*dotk2*zrr2 - part
                     factor = -3.0d0 * xrr2 * part2
                     factork = dotk2*xkrk2*zkrk2 + xkrk2*part
                     dtxdzk1 = dtdzk1*termx + term*factor
                     dtxkdzk1 = dtdzk1*termxk + term*factork
                     factor = -3.0d0 * yrr2 * part2
                     factork = dotk2*ykrk2*zkrk2 + ykrk2*part
                     dtydzk1 = dtdzk1*termy + term*factor
                     dtykdzk1 = dtdzk1*termyk + term*factork
                     factor = 1.0d0 - 3.0d0*zrr2*part2
     &                           + 3.0d0*sk1*dotkr2
                     factork = -sk1 + dotk2*zkrk2*zkrk2
     &                            + zkrk2*part - dotkrk2
                     dtzdzk1 = dtdzk1*termz + term*factor
                     dtzkdzk1 = dtdzk1*termzk + term*factork
                  else if (k2 .eq. i) then
                     dtdxk2 = 3.0d0*sk2*xrr2 + xkrk2
                     part = sk2*xk + xr
                     part2 = sk2*dotk2*xrr2 - part
                     factor = -1.0d0 - 3.0d0*xrr2*part2
     &                           + 3.0d0*sk2*dotkr2
                     factork = -sk2 - dotk2*xkrk2*xkrk2
     &                            + xkrk2*part + dotkrk2
                     dtxdxk2 = dtdxk2*termx + term*factor
                     dtxkdxk2 = dtdxk2*termxk + term*factork
                     factor = -3.0d0 * yrr2 * part2
                     factork = -dotk2*ykrk2*xkrk2 + ykrk2*part
                     dtydxk2 = dtdxk2*termy + term*factor
                     dtykdxk2 = dtdxk2*termyk + term*factork
                     factor = -3.0d0 * zrr2 * part2
                     factork = -dotk2*zkrk2*xkrk2 + zkrk2*part
                     dtzdxk2 = dtdxk2*termz + term*factor
                     dtzkdxk2 = dtdxk2*termzk + term*factork
                     dtdyk2 = 3.0d0*sk2*yrr2 + ykrk2
                     part = sk2*yk + yr
                     part2 = sk2*dotk2*yrr2 - part
                     factor = -3.0d0 * xrr2 * part2
                     factork = -dotk2*xkrk2*ykrk2 + xkrk2*part
                     dtxdyk2 = dtdyk2*termx + term*factor
                     dtxkdyk2 = dtdyk2*termxk + term*factork
                     factor = -1.0d0 - 3.0d0*yrr2*part2
     &                           + 3.0d0*sk2*dotkr2
                     factork = -sk2 - dotk2*ykrk2*ykrk2
     &                            + ykrk2*part + dotkrk2
                     dtydyk2 = dtdyk2*termy + term*factor
                     dtykdyk2 = dtdyk2*termyk + term*factork
                     factor = -3.0d0 * zrr2 * part2
                     factork = -dotk2*zkrk2*ykrk2 + zkrk2*part
                     dtzdyk2 = dtdyk2*termz + term*factor
                     dtzkdyk2 = dtdyk2*termzk + term*factork
                     dtdzk2 = 3.0d0*sk2*zrr2 + zkrk2
                     part = sk2*zk + zr
                     part2 = sk2*dotk2*zrr2 - part
                     factor = -3.0d0 * xrr2 * part2
                     factork = -dotk2*xkrk2*zkrk2 + xkrk2*part
                     dtxdzk2 = dtdzk2*termx + term*factor
                     dtxkdzk2 = dtdzk2*termxk + term*factork
                     factor = -3.0d0 * yrr2 * part2
                     factork = -dotk2*ykrk2*zkrk2 + ykrk2*part
                     dtydzk2 = dtdzk2*termy + term*factor
                     dtykdzk2 = dtdzk2*termyk + term*factork
                     factor = -1.0d0 - 3.0d0*zrr2*part2
     &                           + 3.0d0*sk2*dotkr2
                     factork = -sk2 - dotk2*zkrk2*zkrk2
     &                            + zkrk2*part + dotkrk2
                     dtzdzk2 = dtdzk2*termz + term*factor
                     dtzkdzk2 = dtdzk2*termzk + term*factork
                  end if
c
c     increment diagonal and off-diagonal Hessian elements
c
                  if (i .eq. k1) then
                     hessx(1,i1) = hessx(1,i1) + dtxdxk1
                     hessx(2,i1) = hessx(2,i1) + dtydxk1
                     hessx(3,i1) = hessx(3,i1) + dtzdxk1
                     hessx(1,k1) = hessx(1,k1) - sk1*dtxdxk1 + dtxkdxk1
                     hessx(2,k1) = hessx(2,k1) - sk1*dtydxk1 + dtykdxk1
                     hessx(3,k1) = hessx(3,k1) - sk1*dtzdxk1 + dtzkdxk1
                     hessx(1,k2) = hessx(1,k2) - sk2*dtxdxk1 - dtxkdxk1
                     hessx(2,k2) = hessx(2,k2) - sk2*dtydxk1 - dtykdxk1
                     hessx(3,k2) = hessx(3,k2) - sk2*dtzdxk1 - dtzkdxk1
                     hessy(1,i1) = hessy(1,i1) + dtxdyk1
                     hessy(2,i1) = hessy(2,i1) + dtydyk1
                     hessy(3,i1) = hessy(3,i1) + dtzdyk1
                     hessy(1,k1) = hessy(1,k1) - sk1*dtxdyk1 + dtxkdyk1
                     hessy(2,k1) = hessy(2,k1) - sk1*dtydyk1 + dtykdyk1
                     hessy(3,k1) = hessy(3,k1) - sk1*dtzdyk1 + dtzkdyk1
                     hessy(1,k2) = hessy(1,k2) - sk2*dtxdyk1 - dtxkdyk1
                     hessy(2,k2) = hessy(2,k2) - sk2*dtydyk1 - dtykdyk1
                     hessy(3,k2) = hessy(3,k2) - sk2*dtzdyk1 - dtzkdyk1
                     hessz(1,i1) = hessz(1,i1) + dtxdzk1
                     hessz(2,i1) = hessz(2,i1) + dtydzk1
                     hessz(3,i1) = hessz(3,i1) + dtzdzk1
                     hessz(1,k1) = hessz(1,k1) - sk1*dtxdzk1 + dtxkdzk1
                     hessz(2,k1) = hessz(2,k1) - sk1*dtydzk1 + dtykdzk1
                     hessz(3,k1) = hessz(3,k1) - sk1*dtzdzk1 + dtzkdzk1
                     hessz(1,k2) = hessz(1,k2) - sk2*dtxdzk1 - dtxkdzk1
                     hessz(2,k2) = hessz(2,k2) - sk2*dtydzk1 - dtykdzk1
                     hessz(3,k2) = hessz(3,k2) - sk2*dtzdzk1 - dtzkdzk1
                  else if (i .eq. k2) then
                     hessx(1,i1) = hessx(1,i1) + dtxdxk2
                     hessx(2,i1) = hessx(2,i1) + dtydxk2
                     hessx(3,i1) = hessx(3,i1) + dtzdxk2
                     hessx(1,k1) = hessx(1,k1) - sk1*dtxdxk2 + dtxkdxk2
                     hessx(2,k1) = hessx(2,k1) - sk1*dtydxk2 + dtykdxk2
                     hessx(3,k1) = hessx(3,k1) - sk1*dtzdxk2 + dtzkdxk2
                     hessx(1,k2) = hessx(1,k2) - sk2*dtxdxk2 - dtxkdxk2
                     hessx(2,k2) = hessx(2,k2) - sk2*dtydxk2 - dtykdxk2
                     hessx(3,k2) = hessx(3,k2) - sk2*dtzdxk2 - dtzkdxk2
                     hessy(1,i1) = hessy(1,i1) + dtxdyk2
                     hessy(2,i1) = hessy(2,i1) + dtydyk2
                     hessy(3,i1) = hessy(3,i1) + dtzdyk2
                     hessy(1,k1) = hessy(1,k1) - sk1*dtxdyk2 + dtxkdyk2
                     hessy(2,k1) = hessy(2,k1) - sk1*dtydyk2 + dtykdyk2
                     hessy(3,k1) = hessy(3,k1) - sk1*dtzdyk2 + dtzkdyk2
                     hessy(1,k2) = hessy(1,k2) - sk2*dtxdyk2 - dtxkdyk2
                     hessy(2,k2) = hessy(2,k2) - sk2*dtydyk2 - dtykdyk2
                     hessy(3,k2) = hessy(3,k2) - sk2*dtzdyk2 - dtzkdyk2
                     hessz(1,i1) = hessz(1,i1) + dtxdzk2
                     hessz(2,i1) = hessz(2,i1) + dtydzk2
                     hessz(3,i1) = hessz(3,i1) + dtzdzk2
                     hessz(1,k1) = hessz(1,k1) - sk1*dtxdzk2 + dtxkdzk2
                     hessz(2,k1) = hessz(2,k1) - sk1*dtydzk2 + dtykdzk2
                     hessz(3,k1) = hessz(3,k1) - sk1*dtzdzk2 + dtzkdzk2
                     hessz(1,k2) = hessz(1,k2) - sk2*dtxdzk2 - dtxkdzk2
                     hessz(2,k2) = hessz(2,k2) - sk2*dtydzk2 - dtykdzk2
                     hessz(3,k2) = hessz(3,k2) - sk2*dtzdzk2 - dtzkdzk2
                  end if
c
c     more energy switching if near the cutoff distance
c
                  if (r2.gt.cut2 .and. i.eq.k1) then
                     hessx(1,i1) = hessx(1,i1) - sk1*dtaperx*dedxi1
     &                      + dtaperx*dedxk1 - sk1*d2taperxx
                     hessx(2,i1) = hessx(2,i1) - sk1*dtaperx*dedyi1
     &                      + dtapery*dedxk1 - sk1*d2taperxy
                     hessx(3,i1) = hessx(3,i1) - sk1*dtaperx*dedzi1
     &                      + dtaperz*dedxk1 - sk1*d2taperxz
                     hessx(1,k1) = hessx(1,k1) - sk1*dtaperx*dedxk1
     &                      - sk1*dtaperx*dedxk1 + sk1*sk1*d2taperxx
                     hessx(2,k1) = hessx(2,k1) - sk1*dtaperx*dedyk1
     &                      - sk1*dtapery*dedxk1 + sk1*sk1*d2taperxy
                     hessx(3,k1) = hessx(3,k1) - sk1*dtaperx*dedzk1
     &                      - sk1*dtaperz*dedxk1 + sk1*sk1*d2taperxz
                     hessx(1,k2) = hessx(1,k2) - sk1*dtaperx*dedxk2
     &                      - sk2*dtaperx*dedxk1 + sk1*sk2*d2taperxx
                     hessx(2,k2) = hessx(2,k2) - sk1*dtaperx*dedyk2
     &                      - sk2*dtapery*dedxk1 + sk1*sk2*d2taperxy
                     hessx(3,k2) = hessx(3,k2) - sk1*dtaperx*dedzk2
     &                      - sk2*dtaperz*dedxk1 + sk1*sk2*d2taperxz
                     hessy(1,i1) = hessy(1,i1) - sk1*dtapery*dedxi1
     &                      + dtaperx*dedyk1 - sk1*d2taperxy
                     hessy(2,i1) = hessy(2,i1) - sk1*dtapery*dedyi1
     &                      + dtapery*dedyk1 - sk1*d2taperyy
                     hessy(3,i1) = hessy(3,i1) - sk1*dtapery*dedzi1
     &                      + dtaperz*dedyk1 - sk1*d2taperyz
                     hessy(1,k1) = hessy(1,k1) - sk1*dtapery*dedxk1
     &                      - sk1*dtaperx*dedyk1 + sk1*sk1*d2taperxy
                     hessy(2,k1) = hessy(2,k1) - sk1*dtapery*dedyk1
     &                      - sk1*dtapery*dedyk1 + sk1*sk1*d2taperyy
                     hessy(3,k1) = hessy(3,k1) - sk1*dtapery*dedzk1
     &                      - sk1*dtaperz*dedyk1 + sk1*sk1*d2taperyz
                     hessy(1,k2) = hessy(1,k2) - sk1*dtapery*dedxk2
     &                      - sk2*dtaperx*dedyk1 + sk1*sk2*d2taperxy
                     hessy(2,k2) = hessy(2,k2) - sk1*dtapery*dedyk2
     &                      - sk2*dtapery*dedyk1 + sk1*sk2*d2taperyy
                     hessy(3,k2) = hessy(3,k2) - sk1*dtapery*dedzk2
     &                      - sk2*dtaperz*dedyk1 + sk1*sk2*d2taperyz
                     hessz(1,i1) = hessz(1,i1) - sk1*dtaperz*dedxi1
     &                      + dtaperx*dedzk1 - sk1*d2taperxz
                     hessz(2,i1) = hessz(2,i1) - sk1*dtaperz*dedyi1
     &                      + dtapery*dedzk1 - sk1*d2taperyz
                     hessz(3,i1) = hessz(3,i1) - sk1*dtaperz*dedzi1
     &                      + dtaperz*dedzk1 - sk1*d2taperzz
                     hessz(1,k1) = hessz(1,k1) - sk1*dtaperz*dedxk1
     &                      - sk1*dtaperx*dedzk1 + sk1*sk1*d2taperxz
                     hessz(2,k1) = hessz(2,k1) - sk1*dtaperz*dedyk1
     &                      - sk1*dtapery*dedzk1 + sk1*sk1*d2taperyz
                     hessz(3,k1) = hessz(3,k1) - sk1*dtaperz*dedzk1
     &                      - sk1*dtaperz*dedzk1 + sk1*sk1*d2taperzz
                     hessz(1,k2) = hessz(1,k2) - sk1*dtaperz*dedxk2
     &                      - sk2*dtaperx*dedzk1 + sk1*sk2*d2taperxz
                     hessz(2,k2) = hessz(2,k2) - sk1*dtaperz*dedyk2
     &                      - sk2*dtapery*dedzk1 + sk1*sk2*d2taperyz
                     hessz(3,k2) = hessz(3,k2) - sk1*dtaperz*dedzk2
     &                      - sk2*dtaperz*dedzk1 + sk1*sk2*d2taperzz
                  else if (r2.gt.cut2 .and. i.eq.k2) then
                     hessx(1,i1) = hessx(1,i1) - sk2*dtaperx*dedxi1
     &                      + dtaperx*dedxk2 - sk2*d2taperxx
                     hessx(2,i1) = hessx(2,i1) - sk2*dtaperx*dedyi1
     &                      + dtapery*dedxk2 - sk2*d2taperxy
                     hessx(3,i1) = hessx(3,i1) - sk2*dtaperx*dedzi1
     &                      + dtaperz*dedxk2 - sk2*d2taperxz
                     hessx(1,k1) = hessx(1,k1) - sk2*dtaperx*dedxk1
     &                      - sk1*dtaperx*dedxk2 + sk1*sk2*d2taperxx
                     hessx(2,k1) = hessx(2,k1) - sk2*dtaperx*dedyk1
     &                      - sk1*dtapery*dedxk2 + sk1*sk2*d2taperxy
                     hessx(3,k1) = hessx(3,k1) - sk2*dtaperx*dedzk1
     &                      - sk1*dtaperz*dedxk2 + sk1*sk2*d2taperxz
                     hessx(1,k2) = hessx(1,k2) - sk2*dtaperx*dedxk2
     &                      - sk2*dtaperx*dedxk2 + sk2*sk2*d2taperxx
                     hessx(2,k2) = hessx(2,k2) - sk2*dtaperx*dedyk2
     &                      - sk2*dtapery*dedxk2 + sk2*sk2*d2taperxy
                     hessx(3,k2) = hessx(3,k2) - sk2*dtaperx*dedzk2
     &                      - sk2*dtaperz*dedxk2 + sk2*sk2*d2taperxz
                     hessy(1,i1) = hessy(1,i1) - sk2*dtapery*dedxi1
     &                      + dtaperx*dedyk2 - sk2*d2taperxy
                     hessy(2,i1) = hessy(2,i1) - sk2*dtapery*dedyi1
     &                      + dtapery*dedyk2 - sk2*d2taperyy
                     hessy(3,i1) = hessy(3,i1) - sk2*dtapery*dedzi1
     &                      + dtaperz*dedyk2 - sk2*d2taperyz
                     hessy(1,k1) = hessy(1,k1) - sk2*dtapery*dedxk1
     &                      - sk1*dtaperx*dedyk2 + sk1*sk2*d2taperxy
                     hessy(2,k1) = hessy(2,k1) - sk2*dtapery*dedyk1
     &                      - sk1*dtapery*dedyk2 + sk1*sk2*d2taperyy
                     hessy(3,k1) = hessy(3,k1) - sk2*dtapery*dedzk1
     &                      - sk1*dtaperz*dedyk2 + sk1*sk2*d2taperyz
                     hessy(1,k2) = hessy(1,k2) - sk2*dtapery*dedxk2
     &                      - sk2*dtaperx*dedyk2 + sk2*sk2*d2taperxy
                     hessy(2,k2) = hessy(2,k2) - sk2*dtapery*dedyk2
     &                      - sk2*dtapery*dedyk2 + sk2*sk2*d2taperyy
                     hessy(3,k2) = hessy(3,k2) - sk2*dtapery*dedzk2
     &                      - sk2*dtaperz*dedyk2 + sk2*sk2*d2taperyz
                     hessz(1,i1) = hessz(1,i1) - sk2*dtaperz*dedxi1
     &                      + dtaperx*dedzk2 - sk2*d2taperxz
                     hessz(2,i1) = hessz(2,i1) - sk2*dtaperz*dedyi1
     &                      + dtapery*dedzk2 - sk2*d2taperyz
                     hessz(3,i1) = hessz(3,i1) - sk2*dtaperz*dedzi1
     &                      + dtaperz*dedzk2 - sk2*d2taperzz
                     hessz(1,k1) = hessz(1,k1) - sk2*dtaperz*dedxk1
     &                      - sk1*dtaperx*dedzk2 + sk1*sk2*d2taperxz
                     hessz(2,k1) = hessz(2,k1) - sk2*dtaperz*dedyk1
     &                      - sk1*dtapery*dedzk2 + sk1*sk2*d2taperyz
                     hessz(3,k1) = hessz(3,k1) - sk2*dtaperz*dedzk1
     &                      - sk1*dtaperz*dedzk2 + sk1*sk2*d2taperzz
                     hessz(1,k2) = hessz(1,k2) - sk2*dtaperz*dedxk2
     &                      - sk2*dtaperx*dedzk2 + sk2*sk2*d2taperxz
                     hessz(2,k2) = hessz(2,k2) - sk2*dtaperz*dedyk2
     &                      - sk2*dtapery*dedzk2 + sk2*sk2*d2taperyz
                     hessz(3,k2) = hessz(3,k2) - sk2*dtaperz*dedzk2
     &                      - sk2*dtaperz*dedzk2 + sk2*sk2*d2taperzz
                  end if
               end if
            end if
         end do
   20    continue
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do ii = 1, nion
         i1 = iion(ii)
         if (i1 .ne. i)  goto 30
         skip(i1) = i1
         do k = 1, n12(i1)
            skip(i12(k,i1)) = i1
         end do
         xi = x(i1)
         yi = y(i1)
         zi = z(i1)
         fi = f * pchg(ii)
c
c     decide whether to compute the current interaction
c
         do k = 1, ndipole
            k1 = idpl(1,k)
            k2 = idpl(2,k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i1,k1,k2,0,0,0)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               sk1 = 1.0d0 - sdpl(k)
               sk2 = sdpl(k)
               do jcell = 1, ncell
                  xk = x(k1) - x(k2)
                  yk = y(k1) - y(k2)
                  zk = z(k1) - z(k2)
                  xr = xi - x(k1) + xk*sk2
                  yr = yi - y(k1) + yk*sk2
                  zr = zi - z(k1) + zk*sk2
                  call imager (xr,yr,zr,jcell)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. off2) then
                     rk2 = xk*xk + yk*yk + zk*zk
                     rkr3 = sqrt(rk2*r2) * r2
                     dotk = xk*xr + yk*yr + zk*zr
                     fik = -fi * bdpl(k)
                     if (use_polymer) then
                        if (r2 .lt. polycut2) then
                           if (skip(k1).eq.i1 .or. skip(k2).ne.i1)
     &                        fik = 0.0d0
                        end if
                     end if
c
c     scale the interaction based on its group membership
c
                     if (use_group)  fik = fik * fgrp
c
c     some abbreviations used in various chain rule terms
c
                     xrr2 = xr / r2
                     yrr2 = yr / r2
                     zrr2 = zr / r2
                     xkrk2 = xk / rk2
                     ykrk2 = yk / rk2
                     zkrk2 = zk / rk2
                     dotk2 = 2.0d0 * dotk
                     dotkr2 = dotk / r2
c
c     form the chain rule terms for first derivatives
c
                     term = fik / rkr3
                     term2 = -3.0d0 * dotk
                     termx = term * (xk+xrr2*term2)
                     termy = term * (yk+yrr2*term2)
                     termz = term * (zk+zrr2*term2)
                     termxk = term * (xr-dotk*xkrk2)
                     termyk = term * (yr-dotk*ykrk2)
                     termzk = term * (zr-dotk*zkrk2)
c
c     use energy switching if near the cutoff distance
c
                     if (r2 .gt. cut2) then
                        e = fik * dotk / rkr3
                        dedxi1 = termx
                        dedyi1 = termy
                        dedzi1 = termz
                        dedxk1 = -sk1*termx + termxk
                        dedyk1 = -sk1*termy + termyk
                        dedzk1 = -sk1*termz + termzk
                        dedxk2 = -sk2*termx - termxk
                        dedyk2 = -sk2*termy - termyk
                        dedzk2 = -sk2*termz - termzk
                        r = sqrt(r2)
                        r3 = r2 * r
                        r4 = r2 * r2
                        r5 = r2 * r3
                        taper = c5*r5 + c4*r4 + c3*r3
     &                             + c2*r2 + c1*r + c0
                        dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                              + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                        d2taper = 20.0d0*c5*r3 + 12.0d0*c4*r2
     &                               + 6.0d0*c3*r + 2.0d0*c2
                        dtaper = dtaper / r
                        dtaperx = xr * dtaper
                        dtapery = yr * dtaper
                        dtaperz = zr * dtaper
                        d2taper = e * (d2taper-dtaper)
                        dtaper = e * dtaper
                        d2taperxx = xr*xrr2*d2taper + dtaper
                        d2taperxy = xr*yrr2*d2taper
                        d2taperxz = xr*zrr2*d2taper
                        d2taperyy = yr*yrr2*d2taper + dtaper
                        d2taperyz = yr*zrr2*d2taper
                        d2taperzz = zr*zrr2*d2taper + dtaper
                        term = term * taper
                        termx = termx * taper
                        termy = termy * taper
                        termz = termz * taper
                        termxk = termxk * taper
                        termyk = termyk * taper
                        termzk = termzk * taper
                     end if
c
c     chain rule terms for second derivative components
c
                     dtdxi1 = -3.0d0 * xrr2
                     part = xk - dotk2*xrr2
                     factor = -3.0d0 * (dotkr2 + xrr2*part)
                     factork = 1.0d0 - xk*xkrk2
                     dtxdxi1 = dtdxi1*termx + term*factor
                     dtxkdxi1 = dtdxi1*termxk + term*factork
                     factor = -3.0d0 * yrr2 * part
                     factork = -yk * xkrk2
                     dtydxi1 = dtdxi1*termy + term*factor
                     dtykdxi1 = dtdxi1*termyk + term*factork
                     factor = -3.0d0 * zrr2 * part
                     factork = -zk * xkrk2
                     dtzdxi1 = dtdxi1*termz + term*factor
                     dtzkdxi1 = dtdxi1*termzk + term*factork
                     dtdyi1 = -3.0d0 * yrr2
                     part = yk - dotk2*yrr2
                     factor = -3.0d0 * xrr2 * part
                     factork = -xk * ykrk2
                     dtxdyi1 = dtdyi1*termx + term*factor
                     dtxkdyi1 = dtdyi1*termxk + term*factork
                     factor = -3.0d0 * (dotkr2 + yrr2*part)
                     factork = 1.0d0 - yk*ykrk2
                     dtydyi1 = dtdyi1*termy + term*factor
                     dtykdyi1 = dtdyi1*termyk + term*factork
                     factor = -3.0d0 * zrr2 * part
                     factork = -zk * ykrk2
                     dtzdyi1 = dtdyi1*termz + term*factor
                     dtzkdyi1 = dtdyi1*termzk + term*factork
                     dtdzi1 = -3.0d0 * zrr2
                     part = zk - dotk2*zrr2
                     factor = -3.0d0 * xrr2 * part
                     factork = -xk * zkrk2
                     dtxdzi1 = dtdzi1*termx + term*factor
                     dtxkdzi1 = dtdzi1*termxk + term*factork
                     factor = -3.0d0 * yrr2 * part
                     factork = -yk * zkrk2
                     dtydzi1 = dtdzi1*termy + term*factor
                     dtykdzi1 = dtdzi1*termyk + term*factork
                     factor = -3.0d0 * (dotkr2 + zrr2*part)
                     factork = 1.0d0 - zk*zkrk2
                     dtzdzi1 = dtdzi1*termz + term*factor
                     dtzkdzi1 = dtdzi1*termzk + term*factork
c
c     increment diagonal and off-diagonal Hessian elements
c
                     hessx(1,i1) = hessx(1,i1) + dtxdxi1
                     hessx(2,i1) = hessx(2,i1) + dtydxi1
                     hessx(3,i1) = hessx(3,i1) + dtzdxi1
                     hessx(1,k1) = hessx(1,k1) - sk1*dtxdxi1 + dtxkdxi1
                     hessx(2,k1) = hessx(2,k1) - sk1*dtydxi1 + dtykdxi1
                     hessx(3,k1) = hessx(3,k1) - sk1*dtzdxi1 + dtzkdxi1
                     hessx(1,k2) = hessx(1,k2) - sk2*dtxdxi1 - dtxkdxi1
                     hessx(2,k2) = hessx(2,k2) - sk2*dtydxi1 - dtykdxi1
                     hessx(3,k2) = hessx(3,k2) - sk2*dtzdxi1 - dtzkdxi1
                     hessy(1,i1) = hessy(1,i1) + dtxdyi1
                     hessy(2,i1) = hessy(2,i1) + dtydyi1
                     hessy(3,i1) = hessy(3,i1) + dtzdyi1
                     hessy(1,k1) = hessy(1,k1) - sk1*dtxdyi1 + dtxkdyi1
                     hessy(2,k1) = hessy(2,k1) - sk1*dtydyi1 + dtykdyi1
                     hessy(3,k1) = hessy(3,k1) - sk1*dtzdyi1 + dtzkdyi1
                     hessy(1,k2) = hessy(1,k2) - sk2*dtxdyi1 - dtxkdyi1
                     hessy(2,k2) = hessy(2,k2) - sk2*dtydyi1 - dtykdyi1
                     hessy(3,k2) = hessy(3,k2) - sk2*dtzdyi1 - dtzkdyi1
                     hessz(1,i1) = hessz(1,i1) + dtxdzi1
                     hessz(2,i1) = hessz(2,i1) + dtydzi1
                     hessz(3,i1) = hessz(3,i1) + dtzdzi1
                     hessz(1,k1) = hessz(1,k1) - sk1*dtxdzi1 + dtxkdzi1
                     hessz(2,k1) = hessz(2,k1) - sk1*dtydzi1 + dtykdzi1
                     hessz(3,k1) = hessz(3,k1) - sk1*dtzdzi1 + dtzkdzi1
                     hessz(1,k2) = hessz(1,k2) - sk2*dtxdzi1 - dtxkdzi1
                     hessz(2,k2) = hessz(2,k2) - sk2*dtydzi1 - dtykdzi1
                     hessz(3,k2) = hessz(3,k2) - sk2*dtzdzi1 - dtzkdzi1
c
c     more energy switching if near the cutoff distance
c
                     if (r2 .gt. cut2) then
                        hessx(1,i1) = hessx(1,i1) + dtaperx*dedxi1
     &                                + dtaperx*dedxi1 + d2taperxx
                        hessx(2,i1) = hessx(2,i1) + dtaperx*dedyi1
     &                                + dtapery*dedxi1 + d2taperxy
                        hessx(3,i1) = hessx(3,i1) + dtaperx*dedzi1
     &                                + dtaperz*dedxi1 + d2taperxz
                        hessx(1,k1) = hessx(1,k1) + dtaperx*dedxk1
     &                                - sk1*(dtaperx*dedxi1+d2taperxx)
                        hessx(2,k1) = hessx(2,k1) + dtaperx*dedyk1
     &                                - sk1*(dtapery*dedxi1+d2taperxy)
                        hessx(3,k1) = hessx(3,k1) + dtaperx*dedzk1
     &                                - sk1*(dtaperz*dedxi1+d2taperxz)
                        hessx(1,k2) = hessx(1,k2) + dtaperx*dedxk2
     &                                - sk2*(dtaperx*dedxi1+d2taperxx)
                        hessx(2,k2) = hessx(2,k2) + dtaperx*dedyk2
     &                                - sk2*(dtapery*dedxi1+d2taperxy)
                        hessx(3,k2) = hessx(3,k2) + dtaperx*dedzk2
     &                                - sk2*(dtaperz*dedxi1+d2taperxz)
                        hessy(1,i1) = hessy(1,i1) + dtapery*dedxi1
     &                                + dtaperx*dedyi1 + d2taperxy
                        hessy(2,i1) = hessy(2,i1) + dtapery*dedyi1
     &                                + dtapery*dedyi1 + d2taperyy
                        hessy(3,i1) = hessy(3,i1) + dtapery*dedzi1
     &                                + dtaperz*dedyi1 + d2taperyz
                        hessy(1,k1) = hessy(1,k1) + dtapery*dedxk1
     &                                - sk1*(dtaperx*dedyi1+d2taperxy)
                        hessy(2,k1) = hessy(2,k1) + dtapery*dedyk1
     &                                - sk1*(dtapery*dedyi1+d2taperyy)
                        hessy(3,k1) = hessy(3,k1) + dtapery*dedzk1
     &                                - sk1*(dtaperz*dedyi1+d2taperyz)
                        hessy(1,k2) = hessy(1,k2) + dtapery*dedxk2
     &                                - sk2*(dtaperx*dedyi1+d2taperxy)
                        hessy(2,k2) = hessy(2,k2) + dtapery*dedyk2
     &                                - sk2*(dtapery*dedyi1+d2taperyy)
                        hessy(3,k2) = hessy(3,k2) + dtapery*dedzk2
     &                                - sk2*(dtaperz*dedyi1+d2taperyz)
                        hessz(1,i1) = hessz(1,i1) + dtaperz*dedxi1
     &                                + dtaperx*dedzi1 + d2taperxz
                        hessz(2,i1) = hessz(2,i1) + dtaperz*dedyi1
     &                                + dtapery*dedzi1 + d2taperyz
                        hessz(3,i1) = hessz(3,i1) + dtaperz*dedzi1
     &                                + dtaperz*dedzi1 + d2taperzz
                        hessz(1,k1) = hessz(1,k1) + dtaperz*dedxk1
     &                                - sk1*(dtaperx*dedzi1+d2taperxz)
                        hessz(2,k1) = hessz(2,k1) + dtaperz*dedyk1
     &                                - sk1*(dtapery*dedzi1+d2taperyz)
                        hessz(3,k1) = hessz(3,k1) + dtaperz*dedzk1
     &                                - sk1*(dtaperz*dedzi1+d2taperzz)
                        hessz(1,k2) = hessz(1,k2) + dtaperz*dedxk2
     &                                - sk2*(dtaperx*dedzi1+d2taperxz)
                        hessz(2,k2) = hessz(2,k2) + dtaperz*dedyk2
     &                                - sk2*(dtapery*dedzi1+d2taperyz)
                        hessz(3,k2) = hessz(3,k2) + dtaperz*dedzk2
     &                                - sk2*(dtaperz*dedzi1+d2taperzz)
                     end if
                  end if
               end do
            end if
         end do
   30    continue
      end do
c
c     see if the atom of interest is part of a dipole
c
      do k = 1, ndipole
         k1 = idpl(1,k)
         k2 = idpl(2,k)
         if (k1.ne.i .and. k2.ne.i)  goto 40
         do ii = 1, n12(k1)
            omit(i12(ii,k1)) = k
         end do
         do ii = 1, n12(k2)
            omit(i12(ii,k2)) = k
         end do
         sk1 = 1.0d0 - sdpl(k)
         sk2 = sdpl(k)
         xk = x(k1) - x(k2)
         yk = y(k1) - y(k2)
         zk = z(k1) - z(k2)
         rk2 = xk*xk + yk*yk + zk*zk
         xq = x(k1) - xk*sk2
         yq = y(k1) - yk*sk2
         zq = z(k1) - zk*sk2
         fk = -f * bdpl(k)
c
c     decide whether to compute the current interaction
c
         do ii = 1, nion
            i1 = iion(ii)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i1,k1,k2,0,0,0)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do jcell = 1, ncell
                  xr = x(i1) - xq
                  yr = y(i1) - yq
                  zr = z(i1) - zq
                  call imager (xr,yr,zr,jcell)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. off2) then
                     rkr3 = sqrt(rk2*r2) * r2
                     dotk = xk*xr + yk*yr + zk*zr
                     fik = fk * pchg(ii)
                     if (use_polymer) then
                        if (r2 .lt. polycut2) then
                           if (omit(i1) .ne. k)  fik = 0.0d0
                        end if
                     end if
c
c     scale the interaction based on its group membership
c
                     if (use_group)  fik = fik * fgrp
c
c     some abbreviations used in various chain rule terms
c
                     xrr2 = xr / r2
                     yrr2 = yr / r2
                     zrr2 = zr / r2
                     xkrk2 = xk / rk2
                     ykrk2 = yk / rk2
                     zkrk2 = zk / rk2
                     dotk2 = 2.0d0 * dotk
                     dotkr2 = dotk / r2
                     dotkrk2 = dotk / rk2
c
c     form the chain rule terms for first derivatives
c
                     term = fik / rkr3
                     term2 = -3.0d0 * dotk
                     termx = term * (xk+xrr2*term2)
                     termy = term * (yk+yrr2*term2)
                     termz = term * (zk+zrr2*term2)
                     termxk = term * (xr-dotk*xkrk2)
                     termyk = term * (yr-dotk*ykrk2)
                     termzk = term * (zr-dotk*zkrk2)
c
c     use energy switching if near the cutoff distance
c
                     if (r2 .gt. cut2) then
                        e = fik * dotk / rkr3
                        dedxi1 = termx
                        dedyi1 = termy
                        dedzi1 = termz
                        dedxk1 = -sk1*termx + termxk
                        dedyk1 = -sk1*termy + termyk
                        dedzk1 = -sk1*termz + termzk
                        dedxk2 = -sk2*termx - termxk
                        dedyk2 = -sk2*termy - termyk
                        dedzk2 = -sk2*termz - termzk
                        r = sqrt(r2)
                        r3 = r2 * r
                        r4 = r2 * r2
                        r5 = r2 * r3
                        taper = c5*r5 + c4*r4 + c3*r3
     &                             + c2*r2 + c1*r + c0
                        dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                              + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                        d2taper = 20.0d0*c5*r3 + 12.0d0*c4*r2
     &                               + 6.0d0*c3*r + 2.0d0*c2
                        dtaper = dtaper / r
                        dtaperx = xr * dtaper
                        dtapery = yr * dtaper
                        dtaperz = zr * dtaper
                        d2taper = e * (d2taper-dtaper)
                        dtaper = e * dtaper
                        d2taperxx = xr*xrr2*d2taper + dtaper
                        d2taperxy = xr*yrr2*d2taper
                        d2taperxz = xr*zrr2*d2taper
                        d2taperyy = yr*yrr2*d2taper + dtaper
                        d2taperyz = yr*zrr2*d2taper
                        d2taperzz = zr*zrr2*d2taper + dtaper
                        term = term * taper
                        termx = termx * taper
                        termy = termy * taper
                        termz = termz * taper
                        termxk = termxk * taper
                        termyk = termyk * taper
                        termzk = termzk * taper
                     end if
c
c     chain rule terms for second derivative components
c
                     if (k1 .eq. i) then
                        dtdxk1 = 3.0d0*sk1*xrr2 - xkrk2
                        part = sk1*xk - xr
                        part2 = sk1*dotk2*xrr2 - part
                        factor = 1.0d0 - 3.0d0*xrr2*part2
     &                              + 3.0d0*sk1*dotkr2
                        factork = -sk1 + dotk2*xkrk2*xkrk2
     &                               + xkrk2*part - dotkrk2
                        dtxdxk1 = dtdxk1*termx + term*factor
                        dtxkdxk1 = dtdxk1*termxk + term*factork
                        factor = -3.0d0 * yrr2 * part2
                        factork = dotk2*ykrk2*xkrk2 + ykrk2*part
                        dtydxk1 = dtdxk1*termy + term*factor
                        dtykdxk1 = dtdxk1*termyk + term*factork
                        factor = -3.0d0 * zrr2 * part2
                        factork = dotk2*zkrk2*xkrk2 + zkrk2*part
                        dtzdxk1 = dtdxk1*termz + term*factor
                        dtzkdxk1 = dtdxk1*termzk + term*factork
                        dtdyk1 = 3.0d0*sk1*yrr2 - ykrk2
                        part = sk1*yk - yr
                        part2 = sk1*dotk2*yrr2 - part
                        factor = -3.0d0 * xrr2 * part2
                        factork = dotk2*xkrk2*ykrk2 + xkrk2*part
                        dtxdyk1 = dtdyk1*termx + term*factor
                        dtxkdyk1 = dtdyk1*termxk + term*factork
                        factor = 1.0d0 - 3.0d0*yrr2*part2
     &                              + 3.0d0*sk1*dotkr2
                        factork = -sk1 + dotk2*ykrk2*ykrk2
     &                               + ykrk2*part - dotkrk2
                        dtydyk1 = dtdyk1*termy + term*factor
                        dtykdyk1 = dtdyk1*termyk + term*factork
                        factor = -3.0d0 * zrr2 * part2
                        factork = dotk2*zkrk2*ykrk2 + zkrk2*part
                        dtzdyk1 = dtdyk1*termz + term*factor
                        dtzkdyk1 = dtdyk1*termzk + term*factork
                        dtdzk1 = 3.0d0*sk1*zrr2 - zkrk2
                        part = sk1*zk - zr
                        part2 = sk1*dotk2*zrr2 - part
                        factor = -3.0d0 * xrr2 * part2
                        factork = dotk2*xkrk2*zkrk2 + xkrk2*part
                        dtxdzk1 = dtdzk1*termx + term*factor
                        dtxkdzk1 = dtdzk1*termxk + term*factork
                        factor = -3.0d0 * yrr2 * part2
                        factork = dotk2*ykrk2*zkrk2 + ykrk2*part
                        dtydzk1 = dtdzk1*termy + term*factor
                        dtykdzk1 = dtdzk1*termyk + term*factork
                        factor = 1.0d0 - 3.0d0*zrr2*part2
     &                              + 3.0d0*sk1*dotkr2
                        factork = -sk1 + dotk2*zkrk2*zkrk2
     &                               + zkrk2*part - dotkrk2
                        dtzdzk1 = dtdzk1*termz + term*factor
                        dtzkdzk1 = dtdzk1*termzk + term*factork
                     else if (k2 .eq. i) then
                        dtdxk2 = 3.0d0*sk2*xrr2 + xkrk2
                        part = sk2*xk + xr
                        part2 = sk2*dotk2*xrr2 - part
                        factor = -1.0d0 - 3.0d0*xrr2*part2
     &                              + 3.0d0*sk2*dotkr2
                        factork = -sk2 - dotk2*xkrk2*xkrk2
     &                               + xkrk2*part + dotkrk2
                        dtxdxk2 = dtdxk2*termx + term*factor
                        dtxkdxk2 = dtdxk2*termxk + term*factork
                        factor = -3.0d0 * yrr2 * part2
                        factork = -dotk2*ykrk2*xkrk2 + ykrk2*part
                        dtydxk2 = dtdxk2*termy + term*factor
                        dtykdxk2 = dtdxk2*termyk + term*factork
                        factor = -3.0d0 * zrr2 * part2
                        factork = -dotk2*zkrk2*xkrk2 + zkrk2*part
                        dtzdxk2 = dtdxk2*termz + term*factor
                        dtzkdxk2 = dtdxk2*termzk + term*factork
                        dtdyk2 = 3.0d0*sk2*yrr2 + ykrk2
                        part = sk2*yk + yr
                        part2 = sk2*dotk2*yrr2 - part
                        factor = -3.0d0 * xrr2 * part2
                        factork = -dotk2*xkrk2*ykrk2 + xkrk2*part
                        dtxdyk2 = dtdyk2*termx + term*factor
                        dtxkdyk2 = dtdyk2*termxk + term*factork
                        factor = -1.0d0 - 3.0d0*yrr2*part2
     &                              + 3.0d0*sk2*dotkr2
                        factork = -sk2 - dotk2*ykrk2*ykrk2
     &                               + ykrk2*part + dotkrk2
                        dtydyk2 = dtdyk2*termy + term*factor
                        dtykdyk2 = dtdyk2*termyk + term*factork
                        factor = -3.0d0 * zrr2 * part2
                        factork = -dotk2*zkrk2*ykrk2 + zkrk2*part
                        dtzdyk2 = dtdyk2*termz + term*factor
                        dtzkdyk2 = dtdyk2*termzk + term*factork
                        dtdzk2 = 3.0d0*sk2*zrr2 + zkrk2
                        part = sk2*zk + zr
                        part2 = sk2*dotk2*zrr2 - part
                        factor = -3.0d0 * xrr2 * part2
                        factork = -dotk2*xkrk2*zkrk2 + xkrk2*part
                        dtxdzk2 = dtdzk2*termx + term*factor
                        dtxkdzk2 = dtdzk2*termxk + term*factork
                        factor = -3.0d0 * yrr2 * part2
                        factork = -dotk2*ykrk2*zkrk2 + ykrk2*part
                        dtydzk2 = dtdzk2*termy + term*factor
                        dtykdzk2 = dtdzk2*termyk + term*factork
                        factor = -1.0d0 - 3.0d0*zrr2*part2
     &                              + 3.0d0*sk2*dotkr2
                        factork = -sk2 - dotk2*zkrk2*zkrk2
     &                               + zkrk2*part + dotkrk2
                        dtzdzk2 = dtdzk2*termz + term*factor
                        dtzkdzk2 = dtdzk2*termzk + term*factork
                     end if
c
c     increment diagonal and off-diagonal Hessian elements
c
                     if (i .eq. k1) then
                        hessx(1,i1) = hessx(1,i1) + dtxdxk1
                        hessx(2,i1) = hessx(2,i1) + dtydxk1
                        hessx(3,i1) = hessx(3,i1) + dtzdxk1
                        hessx(1,k1) = hessx(1,k1) - sk1*dtxdxk1
     &                                   + dtxkdxk1
                        hessx(2,k1) = hessx(2,k1) - sk1*dtydxk1
     &                                   + dtykdxk1
                        hessx(3,k1) = hessx(3,k1) - sk1*dtzdxk1
     &                                   + dtzkdxk1
                        hessx(1,k2) = hessx(1,k2) - sk2*dtxdxk1
     &                                   - dtxkdxk1
                        hessx(2,k2) = hessx(2,k2) - sk2*dtydxk1
     &                                   - dtykdxk1
                        hessx(3,k2) = hessx(3,k2) - sk2*dtzdxk1
     &                                   - dtzkdxk1
                        hessy(1,i1) = hessy(1,i1) + dtxdyk1
                        hessy(2,i1) = hessy(2,i1) + dtydyk1
                        hessy(3,i1) = hessy(3,i1) + dtzdyk1
                        hessy(1,k1) = hessy(1,k1) - sk1*dtxdyk1
     &                                   + dtxkdyk1
                        hessy(2,k1) = hessy(2,k1) - sk1*dtydyk1
     &                                   + dtykdyk1
                        hessy(3,k1) = hessy(3,k1) - sk1*dtzdyk1
     &                                   + dtzkdyk1
                        hessy(1,k2) = hessy(1,k2) - sk2*dtxdyk1
     &                                   - dtxkdyk1
                        hessy(2,k2) = hessy(2,k2) - sk2*dtydyk1
     &                                   - dtykdyk1
                        hessy(3,k2) = hessy(3,k2) - sk2*dtzdyk1
     &                                   - dtzkdyk1
                        hessz(1,i1) = hessz(1,i1) + dtxdzk1
                        hessz(2,i1) = hessz(2,i1) + dtydzk1
                        hessz(3,i1) = hessz(3,i1) + dtzdzk1
                        hessz(1,k1) = hessz(1,k1) - sk1*dtxdzk1
     &                                   + dtxkdzk1
                        hessz(2,k1) = hessz(2,k1) - sk1*dtydzk1
     &                                   + dtykdzk1
                        hessz(3,k1) = hessz(3,k1) - sk1*dtzdzk1
     &                                   + dtzkdzk1
                        hessz(1,k2) = hessz(1,k2) - sk2*dtxdzk1
     &                                   - dtxkdzk1
                        hessz(2,k2) = hessz(2,k2) - sk2*dtydzk1
     &                                   - dtykdzk1
                        hessz(3,k2) = hessz(3,k2) - sk2*dtzdzk1
     &                                   - dtzkdzk1
                     else if (i .eq. k2) then
                        hessx(1,i1) = hessx(1,i1) + dtxdxk2
                        hessx(2,i1) = hessx(2,i1) + dtydxk2
                        hessx(3,i1) = hessx(3,i1) + dtzdxk2
                        hessx(1,k1) = hessx(1,k1) - sk1*dtxdxk2
     &                                   + dtxkdxk2
                        hessx(2,k1) = hessx(2,k1) - sk1*dtydxk2
     &                                   + dtykdxk2
                        hessx(3,k1) = hessx(3,k1) - sk1*dtzdxk2
     &                                   + dtzkdxk2
                        hessx(1,k2) = hessx(1,k2) - sk2*dtxdxk2
     &                                   - dtxkdxk2
                        hessx(2,k2) = hessx(2,k2) - sk2*dtydxk2
     &                                   - dtykdxk2
                        hessx(3,k2) = hessx(3,k2) - sk2*dtzdxk2
     &                                   - dtzkdxk2
                        hessy(1,i1) = hessy(1,i1) + dtxdyk2
                        hessy(2,i1) = hessy(2,i1) + dtydyk2
                        hessy(3,i1) = hessy(3,i1) + dtzdyk2
                        hessy(1,k1) = hessy(1,k1) - sk1*dtxdyk2
     &                                   + dtxkdyk2
                        hessy(2,k1) = hessy(2,k1) - sk1*dtydyk2
     &                                   + dtykdyk2
                        hessy(3,k1) = hessy(3,k1) - sk1*dtzdyk2
     &                                   + dtzkdyk2
                        hessy(1,k2) = hessy(1,k2) - sk2*dtxdyk2
     &                                   - dtxkdyk2
                        hessy(2,k2) = hessy(2,k2) - sk2*dtydyk2
     &                                   - dtykdyk2
                        hessy(3,k2) = hessy(3,k2) - sk2*dtzdyk2
     &                                   - dtzkdyk2
                        hessz(1,i1) = hessz(1,i1) + dtxdzk2
                        hessz(2,i1) = hessz(2,i1) + dtydzk2
                        hessz(3,i1) = hessz(3,i1) + dtzdzk2
                        hessz(1,k1) = hessz(1,k1) - sk1*dtxdzk2
     &                                   + dtxkdzk2
                        hessz(2,k1) = hessz(2,k1) - sk1*dtydzk2
     &                                   + dtykdzk2
                        hessz(3,k1) = hessz(3,k1) - sk1*dtzdzk2
     &                                   + dtzkdzk2
                        hessz(1,k2) = hessz(1,k2) - sk2*dtxdzk2
     &                                   - dtxkdzk2
                        hessz(2,k2) = hessz(2,k2) - sk2*dtydzk2
     &                                   - dtykdzk2
                        hessz(3,k2) = hessz(3,k2) - sk2*dtzdzk2
     &                                   - dtzkdzk2
                     end if
c
c     more energy switching if near the cutoff distance
c
                     if (r2.gt.cut2 .and. i.eq.k1) then
                        hessx(1,i1) = hessx(1,i1) - sk1*dtaperx*dedxi1
     &                         + dtaperx*dedxk1 - sk1*d2taperxx
                        hessx(2,i1) = hessx(2,i1) - sk1*dtaperx*dedyi1
     &                         + dtapery*dedxk1 - sk1*d2taperxy
                        hessx(3,i1) = hessx(3,i1) - sk1*dtaperx*dedzi1
     &                         + dtaperz*dedxk1 - sk1*d2taperxz
                        hessx(1,k1) = hessx(1,k1) - sk1*dtaperx*dedxk1
     &                         - sk1*dtaperx*dedxk1 + sk1*sk1*d2taperxx
                        hessx(2,k1) = hessx(2,k1) - sk1*dtaperx*dedyk1
     &                         - sk1*dtapery*dedxk1 + sk1*sk1*d2taperxy
                        hessx(3,k1) = hessx(3,k1) - sk1*dtaperx*dedzk1
     &                         - sk1*dtaperz*dedxk1 + sk1*sk1*d2taperxz
                        hessx(1,k2) = hessx(1,k2) - sk1*dtaperx*dedxk2
     &                         - sk2*dtaperx*dedxk1 + sk1*sk2*d2taperxx
                        hessx(2,k2) = hessx(2,k2) - sk1*dtaperx*dedyk2
     &                         - sk2*dtapery*dedxk1 + sk1*sk2*d2taperxy
                        hessx(3,k2) = hessx(3,k2) - sk1*dtaperx*dedzk2
     &                         - sk2*dtaperz*dedxk1 + sk1*sk2*d2taperxz
                        hessy(1,i1) = hessy(1,i1) - sk1*dtapery*dedxi1
     &                         + dtaperx*dedyk1 - sk1*d2taperxy
                        hessy(2,i1) = hessy(2,i1) - sk1*dtapery*dedyi1
     &                         + dtapery*dedyk1 - sk1*d2taperyy
                        hessy(3,i1) = hessy(3,i1) - sk1*dtapery*dedzi1
     &                         + dtaperz*dedyk1 - sk1*d2taperyz
                        hessy(1,k1) = hessy(1,k1) - sk1*dtapery*dedxk1
     &                         - sk1*dtaperx*dedyk1 + sk1*sk1*d2taperxy
                        hessy(2,k1) = hessy(2,k1) - sk1*dtapery*dedyk1
     &                         - sk1*dtapery*dedyk1 + sk1*sk1*d2taperyy
                        hessy(3,k1) = hessy(3,k1) - sk1*dtapery*dedzk1
     &                         - sk1*dtaperz*dedyk1 + sk1*sk1*d2taperyz
                        hessy(1,k2) = hessy(1,k2) - sk1*dtapery*dedxk2
     &                         - sk2*dtaperx*dedyk1 + sk1*sk2*d2taperxy
                        hessy(2,k2) = hessy(2,k2) - sk1*dtapery*dedyk2
     &                         - sk2*dtapery*dedyk1 + sk1*sk2*d2taperyy
                        hessy(3,k2) = hessy(3,k2) - sk1*dtapery*dedzk2
     &                         - sk2*dtaperz*dedyk1 + sk1*sk2*d2taperyz
                        hessz(1,i1) = hessz(1,i1) - sk1*dtaperz*dedxi1
     &                         + dtaperx*dedzk1 - sk1*d2taperxz
                        hessz(2,i1) = hessz(2,i1) - sk1*dtaperz*dedyi1
     &                         + dtapery*dedzk1 - sk1*d2taperyz
                        hessz(3,i1) = hessz(3,i1) - sk1*dtaperz*dedzi1
     &                         + dtaperz*dedzk1 - sk1*d2taperzz
                        hessz(1,k1) = hessz(1,k1) - sk1*dtaperz*dedxk1
     &                         - sk1*dtaperx*dedzk1 + sk1*sk1*d2taperxz
                        hessz(2,k1) = hessz(2,k1) - sk1*dtaperz*dedyk1
     &                         - sk1*dtapery*dedzk1 + sk1*sk1*d2taperyz
                        hessz(3,k1) = hessz(3,k1) - sk1*dtaperz*dedzk1
     &                         - sk1*dtaperz*dedzk1 + sk1*sk1*d2taperzz
                        hessz(1,k2) = hessz(1,k2) - sk1*dtaperz*dedxk2
     &                         - sk2*dtaperx*dedzk1 + sk1*sk2*d2taperxz
                        hessz(2,k2) = hessz(2,k2) - sk1*dtaperz*dedyk2
     &                         - sk2*dtapery*dedzk1 + sk1*sk2*d2taperyz
                        hessz(3,k2) = hessz(3,k2) - sk1*dtaperz*dedzk2
     &                         - sk2*dtaperz*dedzk1 + sk1*sk2*d2taperzz
                     else if (r2.gt.cut2 .and. i.eq.k2) then
                        hessx(1,i1) = hessx(1,i1) - sk2*dtaperx*dedxi1
     &                         + dtaperx*dedxk2 - sk2*d2taperxx
                        hessx(2,i1) = hessx(2,i1) - sk2*dtaperx*dedyi1
     &                         + dtapery*dedxk2 - sk2*d2taperxy
                        hessx(3,i1) = hessx(3,i1) - sk2*dtaperx*dedzi1
     &                         + dtaperz*dedxk2 - sk2*d2taperxz
                        hessx(1,k1) = hessx(1,k1) - sk2*dtaperx*dedxk1
     &                         - sk1*dtaperx*dedxk2 + sk1*sk2*d2taperxx
                        hessx(2,k1) = hessx(2,k1) - sk2*dtaperx*dedyk1
     &                         - sk1*dtapery*dedxk2 + sk1*sk2*d2taperxy
                        hessx(3,k1) = hessx(3,k1) - sk2*dtaperx*dedzk1
     &                         - sk1*dtaperz*dedxk2 + sk1*sk2*d2taperxz
                        hessx(1,k2) = hessx(1,k2) - sk2*dtaperx*dedxk2
     &                         - sk2*dtaperx*dedxk2 + sk2*sk2*d2taperxx
                        hessx(2,k2) = hessx(2,k2) - sk2*dtaperx*dedyk2
     &                         - sk2*dtapery*dedxk2 + sk2*sk2*d2taperxy
                        hessx(3,k2) = hessx(3,k2) - sk2*dtaperx*dedzk2
     &                         - sk2*dtaperz*dedxk2 + sk2*sk2*d2taperxz
                        hessy(1,i1) = hessy(1,i1) - sk2*dtapery*dedxi1
     &                         + dtaperx*dedyk2 - sk2*d2taperxy
                        hessy(2,i1) = hessy(2,i1) - sk2*dtapery*dedyi1
     &                         + dtapery*dedyk2 - sk2*d2taperyy
                        hessy(3,i1) = hessy(3,i1) - sk2*dtapery*dedzi1
     &                         + dtaperz*dedyk2 - sk2*d2taperyz
                        hessy(1,k1) = hessy(1,k1) - sk2*dtapery*dedxk1
     &                         - sk1*dtaperx*dedyk2 + sk1*sk2*d2taperxy
                        hessy(2,k1) = hessy(2,k1) - sk2*dtapery*dedyk1
     &                         - sk1*dtapery*dedyk2 + sk1*sk2*d2taperyy
                        hessy(3,k1) = hessy(3,k1) - sk2*dtapery*dedzk1
     &                         - sk1*dtaperz*dedyk2 + sk1*sk2*d2taperyz
                        hessy(1,k2) = hessy(1,k2) - sk2*dtapery*dedxk2
     &                         - sk2*dtaperx*dedyk2 + sk2*sk2*d2taperxy
                        hessy(2,k2) = hessy(2,k2) - sk2*dtapery*dedyk2
     &                         - sk2*dtapery*dedyk2 + sk2*sk2*d2taperyy
                        hessy(3,k2) = hessy(3,k2) - sk2*dtapery*dedzk2
     &                         - sk2*dtaperz*dedyk2 + sk2*sk2*d2taperyz
                        hessz(1,i1) = hessz(1,i1) - sk2*dtaperz*dedxi1
     &                         + dtaperx*dedzk2 - sk2*d2taperxz
                        hessz(2,i1) = hessz(2,i1) - sk2*dtaperz*dedyi1
     &                         + dtapery*dedzk2 - sk2*d2taperyz
                        hessz(3,i1) = hessz(3,i1) - sk2*dtaperz*dedzi1
     &                         + dtaperz*dedzk2 - sk2*d2taperzz
                        hessz(1,k1) = hessz(1,k1) - sk2*dtaperz*dedxk1
     &                         - sk1*dtaperx*dedzk2 + sk1*sk2*d2taperxz
                        hessz(2,k1) = hessz(2,k1) - sk2*dtaperz*dedyk1
     &                         - sk1*dtapery*dedzk2 + sk1*sk2*d2taperyz
                        hessz(3,k1) = hessz(3,k1) - sk2*dtaperz*dedzk1
     &                         - sk1*dtaperz*dedzk2 + sk1*sk2*d2taperzz
                        hessz(1,k2) = hessz(1,k2) - sk2*dtaperz*dedxk2
     &                         - sk2*dtaperx*dedzk2 + sk2*sk2*d2taperxz
                        hessz(2,k2) = hessz(2,k2) - sk2*dtaperz*dedyk2
     &                         - sk2*dtapery*dedzk2 + sk2*sk2*d2taperyz
                        hessz(3,k2) = hessz(3,k2) - sk2*dtaperz*dedzk2
     &                         - sk2*dtaperz*dedzk2 + sk2*sk2*d2taperzz
                     end if
                  end if
               end do
            end if
         end do
   40    continue
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (skip)
      deallocate (omit)
      return
      end
