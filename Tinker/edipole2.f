c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine edipole2  --  atomwise dipole-dipole Hessian  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "edipole2" calculates second derivatives of the
c     dipole-dipole interaction energy for a single atom
c
c
      subroutine edipole2 (i)
      use sizes
      use atoms
      use bound
      use cell
      use chgpot
      use dipole
      use group
      use hessn
      use units
      use shunt
      implicit none
      integer i,i1,i2,k1,k2
      integer jcell,idipole,kdipole
      real*8 f,fi,fik,fgrp
      real*8 xi,yi,zi,xk,yk,zk
      real*8 xq,yq,zq,xr,yr,zr
      real*8 e,r2,ri2,rk2,rirkr3
      real*8 doti,dotk,dotp
      real*8 si1,si2,sk1,sk2
      real*8 de,dedr,dedrirk
      real*8 deddoti,deddotk,deddotp
      real*8 termx,termy,termz
      real*8 termxi,termyi,termzi
      real*8 termxk,termyk,termzk
      real*8 enum,r2inv,ri2inv
      real*8 dotik,xrr2,yrr2,zrr2
      real*8 xiri2,yiri2,ziri2
      real*8 xkrk2,ykrk2,zkrk2
      real*8 xixr,xiyr,xizr
      real*8 yixr,yiyr,yizr
      real*8 zixr,ziyr,zizr
      real*8 xkxr,xkyr,xkzr
      real*8 ykxr,ykyr,ykzr
      real*8 zkxr,zkyr,zkzr
      real*8 xixk,xiyk,xizk
      real*8 yixk,yiyk,yizk
      real*8 zixk,ziyk,zizk
      real*8 xrxr,xryr,xrzr
      real*8 yryr,yrzr,zrzr
      real*8 xidotk,yidotk,zidotk
      real*8 xkdoti,ykdoti,zkdoti
      real*8 factor,factori,factork
      real*8 part,partik
      real*8 dedxi1,dedyi1,dedzi1
      real*8 dedxi2,dedyi2,dedzi2
      real*8 dedxk1,dedyk1,dedzk1
      real*8 dedxk2,dedyk2,dedzk2
      real*8 dtdxi1,dtdyi1,dtdzi1
      real*8 dtdxi2,dtdyi2,dtdzi2
      real*8 dtxdxi1,dtxidxi1,dtxkdxi1
      real*8 dtxdxi2,dtxidxi2,dtxkdxi2
      real*8 dtydxi1,dtyidxi1,dtykdxi1
      real*8 dtydxi2,dtyidxi2,dtykdxi2
      real*8 dtzdxi1,dtzidxi1,dtzkdxi1
      real*8 dtzdxi2,dtzidxi2,dtzkdxi2
      real*8 dtxdyi1,dtxidyi1,dtxkdyi1
      real*8 dtxdyi2,dtxidyi2,dtxkdyi2
      real*8 dtydyi1,dtyidyi1,dtykdyi1
      real*8 dtydyi2,dtyidyi2,dtykdyi2
      real*8 dtzdyi1,dtzidyi1,dtzkdyi1
      real*8 dtzdyi2,dtzidyi2,dtzkdyi2
      real*8 dtxdzi1,dtxidzi1,dtxkdzi1
      real*8 dtxdzi2,dtxidzi2,dtxkdzi2
      real*8 dtydzi1,dtyidzi1,dtykdzi1
      real*8 dtydzi2,dtyidzi2,dtykdzi2
      real*8 dtzdzi1,dtzidzi1,dtzkdzi1
      real*8 dtzdzi2,dtzidzi2,dtzkdzi2
      real*8 r,r3,r4,r5
      real*8 taper,dtaper,d2taper
      real*8 dtaperx,dtapery,dtaperz
      real*8 d2taperxx,d2taperyy,d2taperzz
      real*8 d2taperxy,d2taperxz,d2taperyz
      logical proceed
      character*6 mode
c
c
c     set conversion factor and switching function coefficients
c
      if (ndipole .eq. 0)  return
      f = electric / (debye**2 * dielec)
      mode = 'DIPOLE'
      call switch (mode)
c
c     calculate the dipole interaction energy Hessian elements
c
      do idipole = 1, ndipole
         i1 = idpl(1,idipole)
         i2 = idpl(2,idipole)
         si1 = 1.0d0 - sdpl(idipole)
         si2 = sdpl(idipole)
         if (i1.ne.i .and. i2.ne.i)  goto 10
         xi = x(i2) - x(i1)
         yi = y(i2) - y(i1)
         zi = z(i2) - z(i1)
         if (use_polymer)  call imager (xi,yi,zi,-1)
         ri2 = xi*xi + yi*yi + zi*zi
         xq = x(i1) + xi*si2
         yq = y(i1) + yi*si2
         zq = z(i1) + zi*si2
         fi = f * bdpl(idipole)
c
c     decide whether to compute the current interaction
c
         do kdipole = 1, ndipole
            k1 = idpl(1,kdipole)
            k2 = idpl(2,kdipole)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i1,i2,k1,k2,0,0)
            if (proceed)  proceed = (k1.ne.i1 .and. k1.ne.i2 .and.
     &                                 k2.ne.i1 .and. k2.ne.i2)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               sk1 = 1.0d0 - sdpl(kdipole)
               sk2 = sdpl(kdipole)
               xk = x(k2) - x(k1)
               yk = y(k2) - y(k1)
               zk = z(k2) - z(k1)
               if (use_polymer)  call imager (xk,yk,zk,-1)
               xr = xq - x(k1) - xk*sk2
               yr = yq - y(k1) - yk*sk2
               zr = zq - z(k1) - zk*sk2
               call image (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  rk2 = xk*xk + yk*yk + zk*zk
                  rirkr3 = sqrt(ri2*rk2*r2) * r2
                  dotp = xi*xk + yi*yk + zi*zk
                  doti = xi*xr + yi*yr + zi*zr
                  dotk = xk*xr + yk*yr + zk*zr
                  fik = fi * bdpl(kdipole)
c
c     some abbreviations used in various chain rule terms
c
                  dotik = doti * dotk
                  enum = dotp*r2 - 3.0d0*dotik
                  r2inv = 15.0d0 / r2
                  ri2inv = 1.0d0 / ri2
                  xrr2 = xr / r2
                  yrr2 = yr / r2
                  zrr2 = zr / r2
                  xiri2 = xi / ri2
                  yiri2 = yi / ri2
                  ziri2 = zi / ri2
                  xkrk2 = xk / rk2
                  ykrk2 = yk / rk2
                  zkrk2 = zk / rk2
                  xixr = xi * xr
                  xiyr = xi * yr
                  xizr = xi * zr
                  yixr = yi * xr
                  yiyr = yi * yr
                  yizr = yi * zr
                  zixr = zi * xr
                  ziyr = zi * yr
                  zizr = zi * zr
                  xkxr = xk * xr
                  xkyr = xk * yr
                  xkzr = xk * zr
                  ykxr = yk * xr
                  ykyr = yk * yr
                  ykzr = yk * zr
                  zkxr = zk * xr
                  zkyr = zk * yr
                  zkzr = zk * zr
                  xixk = xi * xk
                  xiyk = xi * yk
                  xizk = xi * zk
                  yixk = yi * xk
                  yiyk = yi * yk
                  yizk = yi * zk
                  zixk = zi * xk
                  ziyk = zi * yk
                  zizk = zi * zk
                  xrxr = 3.0d0 * xr * xr
                  xryr = 3.0d0 * xr * yr
                  xrzr = 3.0d0 * xr * zr
                  yryr = 3.0d0 * yr * yr
                  yrzr = 3.0d0 * yr * zr
                  zrzr = 3.0d0 * zr * zr
                  xidotk = xi * dotk
                  yidotk = yi * dotk
                  zidotk = zi * dotk
                  xkdoti = xk * doti
                  ykdoti = yk * doti
                  zkdoti = zk * doti
c
c     scale the interaction based on its group membership
c
                  if (use_group)  fik = fik * fgrp
c
c     form the master chain rule term for derivatives
c
                  de = -fik / (rirkr3*r2)
c
c     form the chain rule terms for first derivatives
c
                  deddotp = -de * r2
                  deddoti = de * 3.0d0*dotk
                  deddotk = de * 3.0d0*doti
                  dedr = de * (3.0d0*dotp-15.0d0*dotik/r2)
                  dedrirk = de * enum
c
c     more first derivative chain rule expressions
c
                  termx = dedr*xr + deddoti*xi + deddotk*xk
                  termy = dedr*yr + deddoti*yi + deddotk*yk
                  termz = dedr*zr + deddoti*zi + deddotk*zk
                  termxi = dedrirk*xiri2 + deddotp*xk + deddoti*xr
                  termyi = dedrirk*yiri2 + deddotp*yk + deddoti*yr
                  termzi = dedrirk*ziri2 + deddotp*zk + deddoti*zr
                  termxk = dedrirk*xkrk2 + deddotp*xi + deddotk*xr
                  termyk = dedrirk*ykrk2 + deddotp*yi + deddotk*yr
                  termzk = dedrirk*zkrk2 + deddotp*zi + deddotk*zr
c
c     use energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     e = fik * (dotp-3.0d0*doti*dotk/r2) / rirkr3
                     dedxi1 = si1*termx - termxi
                     dedyi1 = si1*termy - termyi
                     dedzi1 = si1*termz - termzi
                     dedxi2 = si2*termx + termxi
                     dedyi2 = si2*termy + termyi
                     dedzi2 = si2*termz + termzi
                     dedxk1 = -sk1*termx - termxk
                     dedyk1 = -sk1*termy - termyk
                     dedzk1 = -sk1*termz - termzk
                     dedxk2 = -sk2*termx + termxk
                     dedyk2 = -sk2*termy + termyk
                     dedzk2 = -sk2*termz + termzk
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
                     de = de * taper
                     termx = termx * taper
                     termy = termy * taper
                     termz = termz * taper
                     termxi = termxi * taper
                     termyi = termyi * taper
                     termzi = termzi * taper
                     termxk = termxk * taper
                     termyk = termyk * taper
                     termzk = termzk * taper
                  end if
c
c     chain rule terms for second derivative components
c
                  if (i .eq. i1) then
                     dtdxi1 = -5.0d0*si1*xrr2 + xiri2
                     part = si1*xkdoti - dotk*xr + si1*xidotk
     &                         - 2.0d0*si1*dotik*xrr2
                     partik = -xk*r2 + 2.0d0*si1*dotp*xr
     &                           - 3.0d0*si1*xkdoti + 3.0d0*xr*dotk
     &                           - 3.0d0*si1*xidotk
                     factor = 3.0d0*si1*dotp - 6.0d0*xkxr
     &                           + 6.0d0*si1*xixk - 3.0d0*dotk
     &                           - r2inv*(xr*part+si1*dotik)
                     factori = 3.0d0*si1*dotk + si1*xkxr + xiri2*partik
     &                            - enum*(ri2inv-2.0d0*xiri2*xiri2)
                     factork = r2 + 3.0d0*si1*doti + si1*xixr
     &                            - xrxr + xkrk2*partik
                     dtxdxi1 = dtdxi1*termx + de*factor
                     dtxidxi1 = dtdxi1*termxi + de*factori
                     dtxkdxi1 = dtdxi1*termxk + de*factork
                     factor = -3.0d0*xkyr - 3.0d0*ykxr + 3.0d0*si1*xiyk
     &                           + 3.0d0*si1*yixk - r2inv*yr*part
                     factori = -2.0d0*si1*ykxr + 3.0d0*si1*xkyr
     &                           + yiri2*partik + 2.0d0*enum*yiri2*xiri2
                     factork = -2.0d0*si1*yixr - xryr + 3.0d0*si1*xiyr
     &                            + ykrk2*partik
                     dtydxi1 = dtdxi1*termy + de*factor
                     dtyidxi1 = dtdxi1*termyi + de*factori
                     dtykdxi1 = dtdxi1*termyk + de*factork
                     factor = -3.0d0*xkzr - 3.0d0*zkxr + 3.0d0*si1*xizk
     &                           + 3.0d0*si1*zixk - r2inv*zr*part
                     factori = -2.0d0*si1*zkxr + 3.0d0*si1*xkzr
     &                           + ziri2*partik + 2.0d0*enum*ziri2*xiri2
                     factork = -2.0d0*si1*zixr - xrzr + 3.0d0*si1*xizr
     &                            + zkrk2*partik
                     dtzdxi1 = dtdxi1*termz + de*factor
                     dtzidxi1 = dtdxi1*termzi + de*factori
                     dtzkdxi1 = dtdxi1*termzk + de*factork
                     dtdyi1 = -5.0d0*si1*yrr2 + yiri2
                     part = si1*ykdoti - dotk*yr + si1*yidotk
     &                         - 2.0d0*si1*dotik*yrr2
                     partik = -yk*r2 + 2.0d0*si1*dotp*yr
     &                           - 3.0d0*si1*ykdoti + 3.0d0*yr*dotk
     &                           - 3.0d0*si1*yidotk
                     factor = -3.0d0*ykxr - 3.0d0*xkyr + 3.0d0*si1*yixk
     &                           + 3.0d0*si1*xiyk - r2inv*xr*part
                     factori = -2.0d0*si1*xkyr + 3.0d0*si1*ykxr
     &                           + xiri2*partik + 2.0d0*enum*xiri2*yiri2
                     factork = -2.0d0*si1*xiyr - xryr + 3.0d0*si1*yixr
     &                            + xkrk2*partik
                     dtxdyi1 = dtdyi1*termx + de*factor
                     dtxidyi1 = dtdyi1*termxi + de*factori
                     dtxkdyi1 = dtdyi1*termxk + de*factork
                     factor = 3.0d0*si1*dotp - 6.0d0*ykyr
     &                           + 6.0d0*si1*yiyk - 3.0d0*dotk
     &                           - r2inv*(yr*part+si1*dotik)
                     factori = 3.0d0*si1*dotk + si1*ykyr + yiri2*partik
     &                            - enum*(ri2inv-2.0d0*yiri2*yiri2)
                     factork = r2 + 3.0d0*si1*doti + si1*yiyr
     &                            - yryr + ykrk2*partik
                     dtydyi1 = dtdyi1*termy + de*factor
                     dtyidyi1 = dtdyi1*termyi + de*factori
                     dtykdyi1 = dtdyi1*termyk + de*factork
                     factor = -3.0d0*ykzr - 3.0d0*zkyr + 3.0d0*si1*yizk
     &                           + 3.0d0*si1*ziyk - r2inv*zr*part
                     factori = -2.0d0*si1*zkyr + 3.0d0*si1*ykzr
     &                           + ziri2*partik + 2.0d0*enum*ziri2*yiri2
                     factork = -2.0d0*si1*ziyr - yrzr + 3.0d0*si1*yizr
     &                            + zkrk2*partik
                     dtzdyi1 = dtdyi1*termz + de*factor
                     dtzidyi1 = dtdyi1*termzi + de*factori
                     dtzkdyi1 = dtdyi1*termzk + de*factork
                     dtdzi1 = -5.0d0*si1*zrr2 + ziri2
                     part = si1*zkdoti - dotk*zr + si1*zidotk
     &                         - 2.0d0*si1*dotik*zrr2
                     partik = -zk*r2 + 2.0d0*si1*dotp*zr
     &                           - 3.0d0*si1*zkdoti + 3.0d0*zr*dotk
     &                           - 3.0d0*si1*zidotk
                     factor = -3.0d0*zkxr - 3.0d0*xkzr + 3.0d0*si1*zixk
     &                           + 3.0d0*si1*xizk - r2inv*xr*part
                     factori = -2.0d0*si1*xkzr + 3.0d0*si1*zkxr
     &                           + xiri2*partik + 2.0d0*enum*xiri2*ziri2
                     factork = -2.0d0*si1*xizr - xrzr + 3.0d0*si1*zixr
     &                            + xkrk2*partik
                     dtxdzi1 = dtdzi1*termx + de*factor
                     dtxidzi1 = dtdzi1*termxi + de*factori
                     dtxkdzi1 = dtdzi1*termxk + de*factork
                     factor = -3.0d0*zkyr - 3.0d0*ykzr + 3.0d0*si1*ziyk
     &                           + 3.0d0*si1*yizk - r2inv*yr*part
                     factori = -2.0d0*si1*ykzr + 3.0d0*si1*zkyr
     &                           + yiri2*partik + 2.0d0*enum*yiri2*ziri2
                     factork = -2.0d0*si1*yizr - yrzr + 3.0d0*si1*ziyr
     &                            + ykrk2*partik
                     dtydzi1 = dtdzi1*termy + de*factor
                     dtyidzi1 = dtdzi1*termyi + de*factori
                     dtykdzi1 = dtdzi1*termyk + de*factork
                     factor = 3.0d0*si1*dotp - 6.0d0*zkzr
     &                           + 6.0d0*si1*zizk - 3.0d0*dotk
     &                           - r2inv*(zr*part+si1*dotik)
                     factori = 3.0d0*si1*dotk + si1*zkzr + ziri2*partik
     &                            - enum*(ri2inv-2.0d0*ziri2*ziri2)
                     factork = r2 + 3.0d0*si1*doti + si1*zizr
     &                            - zrzr + zkrk2*partik
                     dtzdzi1 = dtdzi1*termz + de*factor
                     dtzidzi1 = dtdzi1*termzi + de*factori
                     dtzkdzi1 = dtdzi1*termzk + de*factork
                  else if (i .eq. i2) then
                     dtdxi2 = -5.0d0*si2*xrr2 - xiri2
                     part = si2*xkdoti + dotk*xr + si2*xidotk
     &                         - 2.0d0*si2*dotik*xrr2
                     partik = xk*r2 + 2.0d0*si2*dotp*xr
     &                           - 3.0d0*si2*xkdoti - 3.0d0*xr*dotk
     &                           - 3.0d0*si2*xidotk
                     factor = 3.0d0*si2*dotp + 6.0d0*xkxr
     &                           + 6.0d0*si2*xixk + 3.0d0*dotk
     &                           - r2inv*(xr*part+si2*dotik)
                     factori = 3.0d0*si2*dotk + si2*xkxr + xiri2*partik
     &                            + enum*(ri2inv-2.0d0*xiri2*xiri2)
                     factork = -r2 + 3.0d0*si2*doti + si2*xixr
     &                            + xrxr + xkrk2*partik
                     dtxdxi2 = dtdxi2*termx + de*factor
                     dtxidxi2 = dtdxi2*termxi + de*factori
                     dtxkdxi2 = dtdxi2*termxk + de*factork
                     factor = 3.0d0*xkyr + 3.0d0*ykxr + 3.0d0*si2*xiyk
     &                           + 3.0d0*si2*yixk - r2inv*yr*part
                     factori = -2.0d0*si2*ykxr + 3.0d0*si2*xkyr
     &                           + yiri2*partik - 2.0d0*enum*yiri2*xiri2
                     factork = -2.0d0*si2*yixr + xryr + 3.0d0*si2*xiyr
     &                            + ykrk2*partik
                     dtydxi2 = dtdxi2*termy + de*factor
                     dtyidxi2 = dtdxi2*termyi + de*factori
                     dtykdxi2 = dtdxi2*termyk + de*factork
                     factor = 3.0d0*xkzr + 3.0d0*zkxr + 3.0d0*si2*xizk
     &                           + 3.0d0*si2*zixk - r2inv*zr*part
                     factori = -2.0d0*si2*zkxr + 3.0d0*si2*xkzr
     &                           + ziri2*partik - 2.0d0*enum*ziri2*xiri2
                     factork = -2.0d0*si2*zixr + xrzr + 3.0d0*si2*xizr
     &                            + zkrk2*partik
                     dtzdxi2 = dtdxi2*termz + de*factor
                     dtzidxi2 = dtdxi2*termzi + de*factori
                     dtzkdxi2 = dtdxi2*termzk + de*factork
                     dtdyi2 = -5.0d0*si2*yrr2 - yiri2
                     part = si2*ykdoti + dotk*yr + si2*yidotk
     &                         - 2.0d0*si2*dotik*yrr2
                     partik = yk*r2 + 2.0d0*si2*dotp*yr
     &                           - 3.0d0*si2*ykdoti - 3.0d0*yr*dotk
     &                           - 3.0d0*si2*yidotk
                     factor = 3.0d0*ykxr + 3.0d0*xkyr + 3.0d0*si2*yixk
     &                           + 3.0d0*si2*xiyk - r2inv*xr*part
                     factori = -2.0d0*si2*xkyr + 3.0d0*si2*ykxr
     &                           + xiri2*partik - 2.0d0*enum*xiri2*yiri2
                     factork = -2.0d0*si2*xiyr + xryr + 3.0d0*si2*yixr
     &                            + xkrk2*partik
                     dtxdyi2 = dtdyi2*termx + de*factor
                     dtxidyi2 = dtdyi2*termxi + de*factori
                     dtxkdyi2 = dtdyi2*termxk + de*factork
                     factor = 3.0d0*si2*dotp + 6.0d0*ykyr
     &                           + 6.0d0*si2*yiyk + 3.0d0*dotk
     &                           - r2inv*(yr*part+si2*dotik)
                     factori = 3.0d0*si2*dotk + si2*ykyr + yiri2*partik
     &                            + enum*(ri2inv-2.0d0*yiri2*yiri2)
                     factork = -r2 + 3.0d0*si2*doti + si2*yiyr
     &                            + yryr + ykrk2*partik
                     dtydyi2 = dtdyi2*termy + de*factor
                     dtyidyi2 = dtdyi2*termyi + de*factori
                     dtykdyi2 = dtdyi2*termyk + de*factork
                     factor = 3.0d0*ykzr + 3.0d0*zkyr + 3.0d0*si2*yizk
     &                           + 3.0d0*si2*ziyk - r2inv*zr*part
                     factori = -2.0d0*si2*zkyr + 3.0d0*si2*ykzr
     &                           + ziri2*partik - 2.0d0*enum*ziri2*yiri2
                     factork = -2.0d0*si2*ziyr + yrzr + 3.0d0*si2*yizr
     &                            + zkrk2*partik
                     dtzdyi2 = dtdyi2*termz + de*factor
                     dtzidyi2 = dtdyi2*termzi + de*factori
                     dtzkdyi2 = dtdyi2*termzk + de*factork
                     dtdzi2 = -5.0d0*si2*zrr2 - ziri2
                     part = si2*zkdoti + dotk*zr + si2*zidotk
     &                         - 2.0d0*si2*dotik*zrr2
                     partik = zk*r2 + 2.0d0*si2*dotp*zr
     &                           - 3.0d0*si2*zkdoti - 3.0d0*zr*dotk
     &                           - 3.0d0*si2*zidotk
                     factor = 3.0d0*zkxr + 3.0d0*xkzr + 3.0d0*si2*zixk
     &                           + 3.0d0*si2*xizk - r2inv*xr*part
                     factori = -2.0d0*si2*xkzr + 3.0d0*si2*zkxr
     &                           + xiri2*partik - 2.0d0*enum*xiri2*ziri2
                     factork = -2.0d0*si2*xizr + xrzr + 3.0d0*si2*zixr
     &                            + xkrk2*partik
                     dtxdzi2 = dtdzi2*termx + de*factor
                     dtxidzi2 = dtdzi2*termxi + de*factori
                     dtxkdzi2 = dtdzi2*termxk + de*factork
                     factor = 3.0d0*zkyr + 3.0d0*ykzr + 3.0d0*si2*ziyk
     &                           + 3.0d0*si2*yizk - r2inv*yr*part
                     factori = -2.0d0*si2*ykzr + 3.0d0*si2*zkyr
     &                           + yiri2*partik - 2.0d0*enum*yiri2*ziri2
                     factork = -2.0d0*si2*yizr + yrzr + 3.0d0*si2*ziyr
     &                            + ykrk2*partik
                     dtydzi2 = dtdzi2*termy + de*factor
                     dtyidzi2 = dtdzi2*termyi + de*factori
                     dtykdzi2 = dtdzi2*termyk + de*factork
                     factor = 3.0d0*si2*dotp + 6.0d0*zkzr
     &                           + 6.0d0*si2*zizk + 3.0d0*dotk
     &                           - r2inv*(zr*part+si2*dotik)
                     factori = 3.0d0*si2*dotk + si2*zkzr + ziri2*partik
     &                            + enum*(ri2inv-2.0d0*ziri2*ziri2)
                     factork = -r2 + 3.0d0*si2*doti + si2*zizr
     &                            + zrzr + zkrk2*partik
                     dtzdzi2 = dtdzi2*termz + de*factor
                     dtzidzi2 = dtdzi2*termzi + de*factori
                     dtzkdzi2 = dtdzi2*termzk + de*factork
                  end if
c
c     increment diagonal and off-diagonal Hessian elements
c
                  if (i .eq. i1) then
                     hessx(1,i1) = hessx(1,i1) + si1*dtxdxi1 - dtxidxi1
                     hessx(2,i1) = hessx(2,i1) + si1*dtydxi1 - dtyidxi1
                     hessx(3,i1) = hessx(3,i1) + si1*dtzdxi1 - dtzidxi1
                     hessx(1,i2) = hessx(1,i2) + si2*dtxdxi1 + dtxidxi1
                     hessx(2,i2) = hessx(2,i2) + si2*dtydxi1 + dtyidxi1
                     hessx(3,i2) = hessx(3,i2) + si2*dtzdxi1 + dtzidxi1
                     hessx(1,k1) = hessx(1,k1) - sk1*dtxdxi1 - dtxkdxi1
                     hessx(2,k1) = hessx(2,k1) - sk1*dtydxi1 - dtykdxi1
                     hessx(3,k1) = hessx(3,k1) - sk1*dtzdxi1 - dtzkdxi1
                     hessx(1,k2) = hessx(1,k2) - sk2*dtxdxi1 + dtxkdxi1
                     hessx(2,k2) = hessx(2,k2) - sk2*dtydxi1 + dtykdxi1
                     hessx(3,k2) = hessx(3,k2) - sk2*dtzdxi1 + dtzkdxi1
                     hessy(1,i1) = hessy(1,i1) + si1*dtxdyi1 - dtxidyi1
                     hessy(2,i1) = hessy(2,i1) + si1*dtydyi1 - dtyidyi1
                     hessy(3,i1) = hessy(3,i1) + si1*dtzdyi1 - dtzidyi1
                     hessy(1,i2) = hessy(1,i2) + si2*dtxdyi1 + dtxidyi1
                     hessy(3,i2) = hessy(3,i2) + si2*dtzdyi1 + dtzidyi1
                     hessy(2,i2) = hessy(2,i2) + si2*dtydyi1 + dtyidyi1
                     hessy(1,k1) = hessy(1,k1) - sk1*dtxdyi1 - dtxkdyi1
                     hessy(2,k1) = hessy(2,k1) - sk1*dtydyi1 - dtykdyi1
                     hessy(3,k1) = hessy(3,k1) - sk1*dtzdyi1 - dtzkdyi1
                     hessy(1,k2) = hessy(1,k2) - sk2*dtxdyi1 + dtxkdyi1
                     hessy(2,k2) = hessy(2,k2) - sk2*dtydyi1 + dtykdyi1
                     hessy(3,k2) = hessy(3,k2) - sk2*dtzdyi1 + dtzkdyi1
                     hessz(1,i1) = hessz(1,i1) + si1*dtxdzi1 - dtxidzi1
                     hessz(2,i1) = hessz(2,i1) + si1*dtydzi1 - dtyidzi1
                     hessz(3,i1) = hessz(3,i1) + si1*dtzdzi1 - dtzidzi1
                     hessz(1,i2) = hessz(1,i2) + si2*dtxdzi1 + dtxidzi1
                     hessz(2,i2) = hessz(2,i2) + si2*dtydzi1 + dtyidzi1
                     hessz(3,i2) = hessz(3,i2) + si2*dtzdzi1 + dtzidzi1
                     hessz(1,k1) = hessz(1,k1) - sk1*dtxdzi1 - dtxkdzi1
                     hessz(2,k1) = hessz(2,k1) - sk1*dtydzi1 - dtykdzi1
                     hessz(3,k1) = hessz(3,k1) - sk1*dtzdzi1 - dtzkdzi1
                     hessz(1,k2) = hessz(1,k2) - sk2*dtxdzi1 + dtxkdzi1
                     hessz(2,k2) = hessz(2,k2) - sk2*dtydzi1 + dtykdzi1
                     hessz(3,k2) = hessz(3,k2) - sk2*dtzdzi1 + dtzkdzi1
                  else if (i .eq. i2) then
                     hessx(1,i1) = hessx(1,i1) + si1*dtxdxi2 - dtxidxi2
                     hessx(2,i1) = hessx(2,i1) + si1*dtydxi2 - dtyidxi2
                     hessx(3,i1) = hessx(3,i1) + si1*dtzdxi2 - dtzidxi2
                     hessx(1,i2) = hessx(1,i2) + si2*dtxdxi2 + dtxidxi2
                     hessx(2,i2) = hessx(2,i2) + si2*dtydxi2 + dtyidxi2
                     hessx(3,i2) = hessx(3,i2) + si2*dtzdxi2 + dtzidxi2
                     hessx(1,k1) = hessx(1,k1) - sk1*dtxdxi2 - dtxkdxi2
                     hessx(2,k1) = hessx(2,k1) - sk1*dtydxi2 - dtykdxi2
                     hessx(3,k1) = hessx(3,k1) - sk1*dtzdxi2 - dtzkdxi2
                     hessx(1,k2) = hessx(1,k2) - sk2*dtxdxi2 + dtxkdxi2
                     hessx(2,k2) = hessx(2,k2) - sk2*dtydxi2 + dtykdxi2
                     hessx(3,k2) = hessx(3,k2) - sk2*dtzdxi2 + dtzkdxi2
                     hessy(1,i1) = hessy(1,i1) + si1*dtxdyi2 - dtxidyi2
                     hessy(2,i1) = hessy(2,i1) + si1*dtydyi2 - dtyidyi2
                     hessy(3,i1) = hessy(3,i1) + si1*dtzdyi2 - dtzidyi2
                     hessy(1,i2) = hessy(1,i2) + si2*dtxdyi2 + dtxidyi2
                     hessy(2,i2) = hessy(2,i2) + si2*dtydyi2 + dtyidyi2
                     hessy(3,i2) = hessy(3,i2) + si2*dtzdyi2 + dtzidyi2
                     hessy(1,k1) = hessy(1,k1) - sk1*dtxdyi2 - dtxkdyi2
                     hessy(2,k1) = hessy(2,k1) - sk1*dtydyi2 - dtykdyi2
                     hessy(3,k1) = hessy(3,k1) - sk1*dtzdyi2 - dtzkdyi2
                     hessy(1,k2) = hessy(1,k2) - sk2*dtxdyi2 + dtxkdyi2
                     hessy(2,k2) = hessy(2,k2) - sk2*dtydyi2 + dtykdyi2
                     hessy(3,k2) = hessy(3,k2) - sk2*dtzdyi2 + dtzkdyi2
                     hessz(1,i1) = hessz(1,i1) + si1*dtxdzi2 - dtxidzi2
                     hessz(2,i1) = hessz(2,i1) + si1*dtydzi2 - dtyidzi2
                     hessz(3,i1) = hessz(3,i1) + si1*dtzdzi2 - dtzidzi2
                     hessz(1,i2) = hessz(1,i2) + si2*dtxdzi2 + dtxidzi2
                     hessz(2,i2) = hessz(2,i2) + si2*dtydzi2 + dtyidzi2
                     hessz(3,i2) = hessz(3,i2) + si2*dtzdzi2 + dtzidzi2
                     hessz(1,k1) = hessz(1,k1) - sk1*dtxdzi2 - dtxkdzi2
                     hessz(2,k1) = hessz(2,k1) - sk1*dtydzi2 - dtykdzi2
                     hessz(3,k1) = hessz(3,k1) - sk1*dtzdzi2 - dtzkdzi2
                     hessz(1,k2) = hessz(1,k2) - sk2*dtxdzi2 + dtxkdzi2
                     hessz(2,k2) = hessz(2,k2) - sk2*dtydzi2 + dtykdzi2
                     hessz(3,k2) = hessz(3,k2) - sk2*dtzdzi2 + dtzkdzi2
                  end if
c
c     more energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     if (i .eq. i1) then
                        hessx(1,i1) = hessx(1,i1) + si1*dtaperx*dedxi1
     &                         + si1*dtaperx*dedxi1 + si1*si1*d2taperxx
                        hessx(2,i1) = hessx(2,i1) + si1*dtapery*dedxi1
     &                         + si1*dtaperx*dedyi1 + si1*si1*d2taperxy
                        hessx(3,i1) = hessx(3,i1) + si1*dtaperz*dedxi1
     &                         + si1*dtaperx*dedzi1 + si1*si1*d2taperxz
                        hessx(1,i2) = hessx(1,i2) + si2*dtaperx*dedxi1
     &                         + si1*dtaperx*dedxi2 + si2*si1*d2taperxx
                        hessx(2,i2) = hessx(2,i2) + si2*dtapery*dedxi1
     &                         + si1*dtaperx*dedyi2 + si2*si1*d2taperxy
                        hessx(3,i2) = hessx(3,i2) + si2*dtaperz*dedxi1
     &                         + si1*dtaperx*dedzi2 + si2*si1*d2taperxz
                        hessx(1,k1) = hessx(1,k1) - sk1*dtaperx*dedxi1
     &                         + si1*dtaperx*dedxk1 - sk1*si1*d2taperxx
                        hessx(2,k1) = hessx(2,k1) - sk1*dtapery*dedxi1
     &                         + si1*dtaperx*dedyk1 - sk1*si1*d2taperxy
                        hessx(3,k1) = hessx(3,k1) - sk1*dtaperz*dedxi1
     &                         + si1*dtaperx*dedzk1 - sk1*si1*d2taperxz
                        hessx(1,k2) = hessx(1,k2) - sk2*dtaperx*dedxi1
     &                         + si1*dtaperx*dedxk2 - sk2*si1*d2taperxx
                        hessx(2,k2) = hessx(2,k2) - sk2*dtapery*dedxi1
     &                         + si1*dtaperx*dedyk2 - sk2*si1*d2taperxy
                        hessx(3,k2) = hessx(3,k2) - sk2*dtaperz*dedxi1
     &                         + si1*dtaperx*dedzk2 - sk2*si1*d2taperxz
                        hessy(1,i1) = hessy(1,i1) + si1*dtaperx*dedyi1
     &                         + si1*dtapery*dedxi1 + si1*si1*d2taperxy
                        hessy(2,i1) = hessy(2,i1) + si1*dtapery*dedyi1
     &                         + si1*dtapery*dedyi1 + si1*si1*d2taperyy
                        hessy(3,i1) = hessy(3,i1) + si1*dtaperz*dedyi1
     &                         + si1*dtapery*dedzi1 + si1*si1*d2taperyz
                        hessy(1,i2) = hessy(1,i2) + si2*dtaperx*dedyi1
     &                         + si1*dtapery*dedxi2 + si2*si1*d2taperxy
                        hessy(2,i2) = hessy(2,i2) + si2*dtapery*dedyi1
     &                         + si1*dtapery*dedyi2 + si2*si1*d2taperyy
                        hessy(3,i2) = hessy(3,i2) + si2*dtaperz*dedyi1
     &                         + si1*dtapery*dedzi2 + si2*si1*d2taperyz
                        hessy(1,k1) = hessy(1,k1) - sk1*dtaperx*dedyi1
     &                         + si1*dtapery*dedxk1 - sk1*si1*d2taperxy
                        hessy(2,k1) = hessy(2,k1) - sk1*dtapery*dedyi1
     &                         + si1*dtapery*dedyk1 - sk1*si1*d2taperyy
                        hessy(3,k1) = hessy(3,k1) - sk1*dtaperz*dedyi1
     &                         + si1*dtapery*dedzk1 - sk1*si1*d2taperyz
                        hessy(1,k2) = hessy(1,k2) - sk2*dtaperx*dedyi1
     &                         + si1*dtapery*dedxk2 - sk2*si1*d2taperxy
                        hessy(2,k2) = hessy(2,k2) - sk2*dtapery*dedyi1
     &                         + si1*dtapery*dedyk2 - sk2*si1*d2taperyy
                        hessy(3,k2) = hessy(3,k2) - sk2*dtaperz*dedyi1
     &                         + si1*dtapery*dedzk2 - sk2*si1*d2taperyz
                        hessz(1,i1) = hessz(1,i1) + si1*dtaperx*dedzi1
     &                         + si1*dtaperz*dedxi1 + si1*si1*d2taperxz
                        hessz(2,i1) = hessz(2,i1) + si1*dtapery*dedzi1
     &                         + si1*dtaperz*dedyi1 + si1*si1*d2taperyz
                        hessz(3,i1) = hessz(3,i1) + si1*dtaperz*dedzi1
     &                         + si1*dtaperz*dedzi1 + si1*si1*d2taperzz
                        hessz(1,i2) = hessz(1,i2) + si2*dtaperx*dedzi1
     &                         + si1*dtaperz*dedxi2 + si2*si1*d2taperxz
                        hessz(2,i2) = hessz(2,i2) + si2*dtapery*dedzi1
     &                         + si1*dtaperz*dedyi2 + si2*si1*d2taperyz
                        hessz(3,i2) = hessz(3,i2) + si2*dtaperz*dedzi1
     &                         + si1*dtaperz*dedzi2 + si2*si1*d2taperzz
                        hessz(1,k1) = hessz(1,k1) - sk1*dtaperx*dedzi1
     &                         + si1*dtaperz*dedxk1 - sk1*si1*d2taperxz
                        hessz(2,k1) = hessz(2,k1) - sk1*dtapery*dedzi1
     &                         + si1*dtaperz*dedyk1 - sk1*si1*d2taperyz
                        hessz(3,k1) = hessz(3,k1) - sk1*dtaperz*dedzi1
     &                         + si1*dtaperz*dedzk1 - sk1*si1*d2taperzz
                        hessz(1,k2) = hessz(1,k2) - sk2*dtaperx*dedzi1
     &                         + si1*dtaperz*dedxk2 - sk2*si1*d2taperxz
                        hessz(2,k2) = hessz(2,k2) - sk2*dtapery*dedzi1
     &                         + si1*dtaperz*dedyk2 - sk2*si1*d2taperyz
                        hessz(3,k2) = hessz(3,k2) - sk2*dtaperz*dedzi1
     &                         + si1*dtaperz*dedzk2 - sk2*si1*d2taperzz
                     else if (i .eq. i2) then
                        hessx(1,i1) = hessx(1,i1) + si1*dtaperx*dedxi2
     &                         + si2*dtaperx*dedxi1 + si1*si2*d2taperxx
                        hessx(2,i1) = hessx(2,i1) + si1*dtapery*dedxi2
     &                         + si2*dtaperx*dedyi1 + si1*si2*d2taperxy
                        hessx(3,i1) = hessx(3,i1) + si1*dtaperz*dedxi2
     &                         + si2*dtaperx*dedzi1 + si1*si2*d2taperxz
                        hessx(1,i2) = hessx(1,i2) + si2*dtaperx*dedxi2
     &                         + si2*dtaperx*dedxi2 + si2*si2*d2taperxx
                        hessx(2,i2) = hessx(2,i2) + si2*dtapery*dedxi2
     &                         + si2*dtaperx*dedyi2 + si2*si2*d2taperxy
                        hessx(3,i2) = hessx(3,i2) + si2*dtaperz*dedxi2
     &                         + si2*dtaperx*dedzi2 + si2*si2*d2taperxz
                        hessx(1,k1) = hessx(1,k1) - sk1*dtaperx*dedxi2
     &                         + si2*dtaperx*dedxk1 - sk1*si2*d2taperxx
                        hessx(2,k1) = hessx(2,k1) - sk1*dtapery*dedxi2
     &                         + si2*dtaperx*dedyk1 - sk1*si2*d2taperxy
                        hessx(3,k1) = hessx(3,k1) - sk1*dtaperz*dedxi2
     &                         + si2*dtaperx*dedzk1 - sk1*si2*d2taperxz
                        hessx(1,k2) = hessx(1,k2) - sk2*dtaperx*dedxi2
     &                         + si2*dtaperx*dedxk2 - sk2*si2*d2taperxx
                        hessx(2,k2) = hessx(2,k2) - sk2*dtapery*dedxi2
     &                         + si2*dtaperx*dedyk2 - sk2*si2*d2taperxy
                        hessx(3,k2) = hessx(3,k2) - sk2*dtaperz*dedxi2
     &                         + si2*dtaperx*dedzk2 - sk2*si2*d2taperxz
                        hessy(1,i1) = hessy(1,i1) + si1*dtaperx*dedyi2
     &                         + si2*dtapery*dedxi1 + si1*si2*d2taperxy
                        hessy(2,i1) = hessy(2,i1) + si1*dtapery*dedyi2
     &                         + si2*dtapery*dedyi1 + si1*si2*d2taperyy
                        hessy(3,i1) = hessy(3,i1) + si1*dtaperz*dedyi2
     &                         + si2*dtapery*dedzi1 + si1*si2*d2taperyz
                        hessy(1,i2) = hessy(1,i2) + si2*dtaperx*dedyi2
     &                         + si2*dtapery*dedxi2 + si2*si2*d2taperxy
                        hessy(2,i2) = hessy(2,i2) + si2*dtapery*dedyi2
     &                         + si2*dtapery*dedyi2 + si2*si2*d2taperyy
                        hessy(3,i2) = hessy(3,i2) + si2*dtaperz*dedyi2
     &                         + si2*dtapery*dedzi2 + si2*si2*d2taperyz
                        hessy(1,k1) = hessy(1,k1) - sk1*dtaperx*dedyi2
     &                         + si2*dtapery*dedxk1 - sk1*si2*d2taperxy
                        hessy(2,k1) = hessy(2,k1) - sk1*dtapery*dedyi2
     &                         + si2*dtapery*dedyk1 - sk1*si2*d2taperyy
                        hessy(3,k1) = hessy(3,k1) - sk1*dtaperz*dedyi2
     &                         + si2*dtapery*dedzk1 - sk1*si2*d2taperyz
                        hessy(1,k2) = hessy(1,k2) - sk2*dtaperx*dedyi2
     &                         + si2*dtapery*dedxk2 - sk2*si2*d2taperxy
                        hessy(2,k2) = hessy(2,k2) - sk2*dtapery*dedyi2
     &                         + si2*dtapery*dedyk2 - sk2*si2*d2taperyy
                        hessy(3,k2) = hessy(3,k2) - sk2*dtaperz*dedyi2
     &                         + si2*dtapery*dedzk2 - sk2*si2*d2taperyz
                        hessz(1,i1) = hessz(1,i1) + si1*dtaperx*dedzi2
     &                         + si2*dtaperz*dedxi1 + si1*si2*d2taperxz
                        hessz(2,i1) = hessz(2,i1) + si1*dtapery*dedzi2
     &                         + si2*dtaperz*dedyi1 + si1*si2*d2taperyz
                        hessz(3,i1) = hessz(3,i1) + si1*dtaperz*dedzi2
     &                         + si2*dtaperz*dedzi1 + si1*si2*d2taperzz
                        hessz(1,i2) = hessz(1,i2) + si2*dtaperx*dedzi2
     &                         + si2*dtaperz*dedxi2 + si2*si2*d2taperxz
                        hessz(2,i2) = hessz(2,i2) + si2*dtapery*dedzi2
     &                         + si2*dtaperz*dedyi2 + si2*si2*d2taperyz
                        hessz(3,i2) = hessz(3,i2) + si2*dtaperz*dedzi2
     &                         + si2*dtaperz*dedzi2 + si2*si2*d2taperzz
                        hessz(1,k1) = hessz(1,k1) - sk1*dtaperx*dedzi2
     &                         + si2*dtaperz*dedxk1 - sk1*si2*d2taperxz
                        hessz(2,k1) = hessz(2,k1) - sk1*dtapery*dedzi2
     &                         + si2*dtaperz*dedyk1 - sk1*si2*d2taperyz
                        hessz(3,k1) = hessz(3,k1) - sk1*dtaperz*dedzi2
     &                         + si2*dtaperz*dedzk1 - sk1*si2*d2taperzz
                        hessz(1,k2) = hessz(1,k2) - sk2*dtaperx*dedzi2
     &                         + si2*dtaperz*dedxk2 - sk2*si2*d2taperxz
                        hessz(2,k2) = hessz(2,k2) - sk2*dtapery*dedzi2
     &                         + si2*dtaperz*dedyk2 - sk2*si2*d2taperyz
                        hessz(3,k2) = hessz(3,k2) - sk2*dtaperz*dedzi2
     &                         + si2*dtaperz*dedzk2 - sk2*si2*d2taperzz
                     end if
                  end if
               end if
            end if
         end do
   10    continue
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do idipole = 1, ndipole
         i1 = idpl(1,idipole)
         i2 = idpl(2,idipole)
         si1 = 1.0d0 - sdpl(idipole)
         si2 = sdpl(idipole)
         if (i1.ne.i .and. i2.ne.i)  goto 30
         xi = x(i2) - x(i1)
         yi = y(i2) - y(i1)
         zi = z(i2) - z(i1)
         if (use_polymer)  call imager (xi,yi,zi,-1)
         ri2 = xi*xi + yi*yi + zi*zi
         xq = x(i1) + xi*si2
         yq = y(i1) + yi*si2
         zq = z(i1) + zi*si2
         fi = f * bdpl(idipole)
c
c     decide whether to compute the current interaction
c
         do kdipole = 1, ndipole
            k1 = idpl(1,kdipole)
            k2 = idpl(2,kdipole)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i1,i2,k1,k2,0,0)
            if (.not. proceed)  goto 20
c
c     compute the energy contribution for this interaction
c
            sk1 = 1.0d0 - sdpl(kdipole)
            sk2 = sdpl(kdipole)
            do jcell = 1, ncell
               xk = x(k2) - x(k1)
               yk = y(k2) - y(k1)
               zk = z(k2) - z(k1)
               if (use_polymer)  call imager (xk,yk,zk,-1)
               xr = xq - x(k1) - xk*sk2
               yr = yq - y(k1) - yk*sk2
               zr = zq - z(k1) - zk*sk2
               call imager (xr,yr,zr,jcell)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  rk2 = xk*xk + yk*yk + zk*zk
                  rirkr3 = sqrt(ri2*rk2*r2) * r2
                  dotp = xi*xk + yi*yk + zi*zk
                  doti = xi*xr + yi*yr + zi*zr
                  dotk = xk*xr + yk*yr + zk*zr
                  fik = fi * bdpl(kdipole)
                  if (use_polymer) then
                     if (r2 .lt. polycut2) then
                        if (k1.eq.i1 .or. k1.eq.i2 .or.
     &                      k2.eq.i1 .or. k2.eq.i2)  fik = 0.0d0
                     end if
                  end if
c
c     some abbreviations used in various chain rule terms
c
                  dotik = doti * dotk
                  enum = dotp*r2 - 3.0d0*dotik
                  r2inv = 15.0d0 / r2
                  ri2inv = 1.0d0 / ri2
                  xrr2 = xr / r2
                  yrr2 = yr / r2
                  zrr2 = zr / r2
                  xiri2 = xi / ri2
                  yiri2 = yi / ri2
                  ziri2 = zi / ri2
                  xkrk2 = xk / rk2
                  ykrk2 = yk / rk2
                  zkrk2 = zk / rk2
                  xixr = xi * xr
                  xiyr = xi * yr
                  xizr = xi * zr
                  yixr = yi * xr
                  yiyr = yi * yr
                  yizr = yi * zr
                  zixr = zi * xr
                  ziyr = zi * yr
                  zizr = zi * zr
                  xkxr = xk * xr
                  xkyr = xk * yr
                  xkzr = xk * zr
                  ykxr = yk * xr
                  ykyr = yk * yr
                  ykzr = yk * zr
                  zkxr = zk * xr
                  zkyr = zk * yr
                  zkzr = zk * zr
                  xixk = xi * xk
                  xiyk = xi * yk
                  xizk = xi * zk
                  yixk = yi * xk
                  yiyk = yi * yk
                  yizk = yi * zk
                  zixk = zi * xk
                  ziyk = zi * yk
                  zizk = zi * zk
                  xrxr = 3.0d0 * xr * xr
                  xryr = 3.0d0 * xr * yr
                  xrzr = 3.0d0 * xr * zr
                  yryr = 3.0d0 * yr * yr
                  yrzr = 3.0d0 * yr * zr
                  zrzr = 3.0d0 * zr * zr
                  xidotk = xi * dotk
                  yidotk = yi * dotk
                  zidotk = zi * dotk
                  xkdoti = xk * doti
                  ykdoti = yk * doti
                  zkdoti = zk * doti
c
c     scale the interaction based on its group membership
c
                  if (use_group)  fik = fik * fgrp
c
c     form the master chain rule term for derivatives
c
                  de = -fik / (rirkr3*r2)
c
c     form the chain rule terms for first derivatives
c
                  deddotp = -de * r2
                  deddoti = de * 3.0d0*dotk
                  deddotk = de * 3.0d0*doti
                  dedr = de * (3.0d0*dotp-15.0d0*dotik/r2)
                  dedrirk = de * enum
c
c     more first derivative chain rule expressions
c
                  termx = dedr*xr + deddoti*xi + deddotk*xk
                  termy = dedr*yr + deddoti*yi + deddotk*yk
                  termz = dedr*zr + deddoti*zi + deddotk*zk
                  termxi = dedrirk*xiri2 + deddotp*xk + deddoti*xr
                  termyi = dedrirk*yiri2 + deddotp*yk + deddoti*yr
                  termzi = dedrirk*ziri2 + deddotp*zk + deddoti*zr
                  termxk = dedrirk*xkrk2 + deddotp*xi + deddotk*xr
                  termyk = dedrirk*ykrk2 + deddotp*yi + deddotk*yr
                  termzk = dedrirk*zkrk2 + deddotp*zi + deddotk*zr
c
c     use energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     e = fik * (dotp-3.0d0*doti*dotk/r2) / rirkr3
                     dedxi1 = si1*termx - termxi
                     dedyi1 = si1*termy - termyi
                     dedzi1 = si1*termz - termzi
                     dedxi2 = si2*termx + termxi
                     dedyi2 = si2*termy + termyi
                     dedzi2 = si2*termz + termzi
                     dedxk1 = -sk1*termx - termxk
                     dedyk1 = -sk1*termy - termyk
                     dedzk1 = -sk1*termz - termzk
                     dedxk2 = -sk2*termx + termxk
                     dedyk2 = -sk2*termy + termyk
                     dedzk2 = -sk2*termz + termzk
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
                     de = de * taper
                     termx = termx * taper
                     termy = termy * taper
                     termz = termz * taper
                     termxi = termxi * taper
                     termyi = termyi * taper
                     termzi = termzi * taper
                     termxk = termxk * taper
                     termyk = termyk * taper
                     termzk = termzk * taper
                  end if
c
c     chain rule terms for second derivative components
c
                  if (i .eq. i1) then
                     dtdxi1 = -5.0d0*si1*xrr2 + xiri2
                     part = si1*xkdoti - dotk*xr + si1*xidotk
     &                         - 2.0d0*si1*dotik*xrr2
                     partik = -xk*r2 + 2.0d0*si1*dotp*xr
     &                           - 3.0d0*si1*xkdoti + 3.0d0*xr*dotk
     &                           - 3.0d0*si1*xidotk
                     factor = 3.0d0*si1*dotp - 6.0d0*xkxr
     &                           + 6.0d0*si1*xixk - 3.0d0*dotk
     &                           - r2inv*(xr*part+si1*dotik)
                     factori = 3.0d0*si1*dotk + si1*xkxr + xiri2*partik
     &                            - enum*(ri2inv-2.0d0*xiri2*xiri2)
                     factork = r2 + 3.0d0*si1*doti + si1*xixr
     &                            - xrxr + xkrk2*partik
                     dtxdxi1 = dtdxi1*termx + de*factor
                     dtxidxi1 = dtdxi1*termxi + de*factori
                     dtxkdxi1 = dtdxi1*termxk + de*factork
                     factor = -3.0d0*xkyr - 3.0d0*ykxr + 3.0d0*si1*xiyk
     &                           + 3.0d0*si1*yixk - r2inv*yr*part
                     factori = -2.0d0*si1*ykxr + 3.0d0*si1*xkyr
     &                           + yiri2*partik + 2.0d0*enum*yiri2*xiri2
                     factork = -2.0d0*si1*yixr - xryr + 3.0d0*si1*xiyr
     &                            + ykrk2*partik
                     dtydxi1 = dtdxi1*termy + de*factor
                     dtyidxi1 = dtdxi1*termyi + de*factori
                     dtykdxi1 = dtdxi1*termyk + de*factork
                     factor = -3.0d0*xkzr - 3.0d0*zkxr + 3.0d0*si1*xizk
     &                           + 3.0d0*si1*zixk - r2inv*zr*part
                     factori = -2.0d0*si1*zkxr + 3.0d0*si1*xkzr
     &                           + ziri2*partik + 2.0d0*enum*ziri2*xiri2
                     factork = -2.0d0*si1*zixr - xrzr + 3.0d0*si1*xizr
     &                            + zkrk2*partik
                     dtzdxi1 = dtdxi1*termz + de*factor
                     dtzidxi1 = dtdxi1*termzi + de*factori
                     dtzkdxi1 = dtdxi1*termzk + de*factork
                     dtdyi1 = -5.0d0*si1*yrr2 + yiri2
                     part = si1*ykdoti - dotk*yr + si1*yidotk
     &                         - 2.0d0*si1*dotik*yrr2
                     partik = -yk*r2 + 2.0d0*si1*dotp*yr
     &                           - 3.0d0*si1*ykdoti + 3.0d0*yr*dotk
     &                           - 3.0d0*si1*yidotk
                     factor = -3.0d0*ykxr - 3.0d0*xkyr + 3.0d0*si1*yixk
     &                           + 3.0d0*si1*xiyk - r2inv*xr*part
                     factori = -2.0d0*si1*xkyr + 3.0d0*si1*ykxr
     &                           + xiri2*partik + 2.0d0*enum*xiri2*yiri2
                     factork = -2.0d0*si1*xiyr - xryr + 3.0d0*si1*yixr
     &                            + xkrk2*partik
                     dtxdyi1 = dtdyi1*termx + de*factor
                     dtxidyi1 = dtdyi1*termxi + de*factori
                     dtxkdyi1 = dtdyi1*termxk + de*factork
                     factor = 3.0d0*si1*dotp - 6.0d0*ykyr
     &                           + 6.0d0*si1*yiyk - 3.0d0*dotk
     &                           - r2inv*(yr*part+si1*dotik)
                     factori = 3.0d0*si1*dotk + si1*ykyr + yiri2*partik
     &                            - enum*(ri2inv-2.0d0*yiri2*yiri2)
                     factork = r2 + 3.0d0*si1*doti + si1*yiyr
     &                            - yryr + ykrk2*partik
                     dtydyi1 = dtdyi1*termy + de*factor
                     dtyidyi1 = dtdyi1*termyi + de*factori
                     dtykdyi1 = dtdyi1*termyk + de*factork
                     factor = -3.0d0*ykzr - 3.0d0*zkyr + 3.0d0*si1*yizk
     &                           + 3.0d0*si1*ziyk - r2inv*zr*part
                     factori = -2.0d0*si1*zkyr + 3.0d0*si1*ykzr
     &                           + ziri2*partik + 2.0d0*enum*ziri2*yiri2
                     factork = -2.0d0*si1*ziyr - yrzr + 3.0d0*si1*yizr
     &                            + zkrk2*partik
                     dtzdyi1 = dtdyi1*termz + de*factor
                     dtzidyi1 = dtdyi1*termzi + de*factori
                     dtzkdyi1 = dtdyi1*termzk + de*factork
                     dtdzi1 = -5.0d0*si1*zrr2 + ziri2
                     part = si1*zkdoti - dotk*zr + si1*zidotk
     &                         - 2.0d0*si1*dotik*zrr2
                     partik = -zk*r2 + 2.0d0*si1*dotp*zr
     &                           - 3.0d0*si1*zkdoti + 3.0d0*zr*dotk
     &                           - 3.0d0*si1*zidotk
                     factor = -3.0d0*zkxr - 3.0d0*xkzr + 3.0d0*si1*zixk
     &                           + 3.0d0*si1*xizk - r2inv*xr*part
                     factori = -2.0d0*si1*xkzr + 3.0d0*si1*zkxr
     &                           + xiri2*partik + 2.0d0*enum*xiri2*ziri2
                     factork = -2.0d0*si1*xizr - xrzr + 3.0d0*si1*zixr
     &                            + xkrk2*partik
                     dtxdzi1 = dtdzi1*termx + de*factor
                     dtxidzi1 = dtdzi1*termxi + de*factori
                     dtxkdzi1 = dtdzi1*termxk + de*factork
                     factor = -3.0d0*zkyr - 3.0d0*ykzr + 3.0d0*si1*ziyk
     &                           + 3.0d0*si1*yizk - r2inv*yr*part
                     factori = -2.0d0*si1*ykzr + 3.0d0*si1*zkyr
     &                           + yiri2*partik + 2.0d0*enum*yiri2*ziri2
                     factork = -2.0d0*si1*yizr - yrzr + 3.0d0*si1*ziyr
     &                            + ykrk2*partik
                     dtydzi1 = dtdzi1*termy + de*factor
                     dtyidzi1 = dtdzi1*termyi + de*factori
                     dtykdzi1 = dtdzi1*termyk + de*factork
                     factor = 3.0d0*si1*dotp - 6.0d0*zkzr
     &                           + 6.0d0*si1*zizk - 3.0d0*dotk
     &                           - r2inv*(zr*part+si1*dotik)
                     factori = 3.0d0*si1*dotk + si1*zkzr + ziri2*partik
     &                            - enum*(ri2inv-2.0d0*ziri2*ziri2)
                     factork = r2 + 3.0d0*si1*doti + si1*zizr
     &                            - zrzr + zkrk2*partik
                     dtzdzi1 = dtdzi1*termz + de*factor
                     dtzidzi1 = dtdzi1*termzi + de*factori
                     dtzkdzi1 = dtdzi1*termzk + de*factork
                  else if (i .eq. i2) then
                     dtdxi2 = -5.0d0*si2*xrr2 - xiri2
                     part = si2*xkdoti + dotk*xr + si2*xidotk
     &                         - 2.0d0*si2*dotik*xrr2
                     partik = xk*r2 + 2.0d0*si2*dotp*xr
     &                           - 3.0d0*si2*xkdoti - 3.0d0*xr*dotk
     &                           - 3.0d0*si2*xidotk
                     factor = 3.0d0*si2*dotp + 6.0d0*xkxr
     &                           + 6.0d0*si2*xixk + 3.0d0*dotk
     &                           - r2inv*(xr*part+si2*dotik)
                     factori = 3.0d0*si2*dotk + si2*xkxr + xiri2*partik
     &                            + enum*(ri2inv-2.0d0*xiri2*xiri2)
                     factork = -r2 + 3.0d0*si2*doti + si2*xixr
     &                            + xrxr + xkrk2*partik
                     dtxdxi2 = dtdxi2*termx + de*factor
                     dtxidxi2 = dtdxi2*termxi + de*factori
                     dtxkdxi2 = dtdxi2*termxk + de*factork
                     factor = 3.0d0*xkyr + 3.0d0*ykxr + 3.0d0*si2*xiyk
     &                           + 3.0d0*si2*yixk - r2inv*yr*part
                     factori = -2.0d0*si2*ykxr + 3.0d0*si2*xkyr
     &                           + yiri2*partik - 2.0d0*enum*yiri2*xiri2
                     factork = -2.0d0*si2*yixr + xryr + 3.0d0*si2*xiyr
     &                            + ykrk2*partik
                     dtydxi2 = dtdxi2*termy + de*factor
                     dtyidxi2 = dtdxi2*termyi + de*factori
                     dtykdxi2 = dtdxi2*termyk + de*factork
                     factor = 3.0d0*xkzr + 3.0d0*zkxr + 3.0d0*si2*xizk
     &                           + 3.0d0*si2*zixk - r2inv*zr*part
                     factori = -2.0d0*si2*zkxr + 3.0d0*si2*xkzr
     &                           + ziri2*partik - 2.0d0*enum*ziri2*xiri2
                     factork = -2.0d0*si2*zixr + xrzr + 3.0d0*si2*xizr
     &                            + zkrk2*partik
                     dtzdxi2 = dtdxi2*termz + de*factor
                     dtzidxi2 = dtdxi2*termzi + de*factori
                     dtzkdxi2 = dtdxi2*termzk + de*factork
                     dtdyi2 = -5.0d0*si2*yrr2 - yiri2
                     part = si2*ykdoti + dotk*yr + si2*yidotk
     &                         - 2.0d0*si2*dotik*yrr2
                     partik = yk*r2 + 2.0d0*si2*dotp*yr
     &                           - 3.0d0*si2*ykdoti - 3.0d0*yr*dotk
     &                           - 3.0d0*si2*yidotk
                     factor = 3.0d0*ykxr + 3.0d0*xkyr + 3.0d0*si2*yixk
     &                           + 3.0d0*si2*xiyk - r2inv*xr*part
                     factori = -2.0d0*si2*xkyr + 3.0d0*si2*ykxr
     &                           + xiri2*partik - 2.0d0*enum*xiri2*yiri2
                     factork = -2.0d0*si2*xiyr + xryr + 3.0d0*si2*yixr
     &                            + xkrk2*partik
                     dtxdyi2 = dtdyi2*termx + de*factor
                     dtxidyi2 = dtdyi2*termxi + de*factori
                     dtxkdyi2 = dtdyi2*termxk + de*factork
                     factor = 3.0d0*si2*dotp + 6.0d0*ykyr
     &                           + 6.0d0*si2*yiyk + 3.0d0*dotk
     &                           - r2inv*(yr*part+si2*dotik)
                     factori = 3.0d0*si2*dotk + si2*ykyr + yiri2*partik
     &                            + enum*(ri2inv-2.0d0*yiri2*yiri2)
                     factork = -r2 + 3.0d0*si2*doti + si2*yiyr
     &                            + yryr + ykrk2*partik
                     dtydyi2 = dtdyi2*termy + de*factor
                     dtyidyi2 = dtdyi2*termyi + de*factori
                     dtykdyi2 = dtdyi2*termyk + de*factork
                     factor = 3.0d0*ykzr + 3.0d0*zkyr + 3.0d0*si2*yizk
     &                           + 3.0d0*si2*ziyk - r2inv*zr*part
                     factori = -2.0d0*si2*zkyr + 3.0d0*si2*ykzr
     &                           + ziri2*partik - 2.0d0*enum*ziri2*yiri2
                     factork = -2.0d0*si2*ziyr + yrzr + 3.0d0*si2*yizr
     &                            + zkrk2*partik
                     dtzdyi2 = dtdyi2*termz + de*factor
                     dtzidyi2 = dtdyi2*termzi + de*factori
                     dtzkdyi2 = dtdyi2*termzk + de*factork
                     dtdzi2 = -5.0d0*si2*zrr2 - ziri2
                     part = si2*zkdoti + dotk*zr + si2*zidotk
     &                         - 2.0d0*si2*dotik*zrr2
                     partik = zk*r2 + 2.0d0*si2*dotp*zr
     &                           - 3.0d0*si2*zkdoti - 3.0d0*zr*dotk
     &                           - 3.0d0*si2*zidotk
                     factor = 3.0d0*zkxr + 3.0d0*xkzr + 3.0d0*si2*zixk
     &                           + 3.0d0*si2*xizk - r2inv*xr*part
                     factori = -2.0d0*si2*xkzr + 3.0d0*si2*zkxr
     &                           + xiri2*partik - 2.0d0*enum*xiri2*ziri2
                     factork = -2.0d0*si2*xizr + xrzr + 3.0d0*si2*zixr
     &                            + xkrk2*partik
                     dtxdzi2 = dtdzi2*termx + de*factor
                     dtxidzi2 = dtdzi2*termxi + de*factori
                     dtxkdzi2 = dtdzi2*termxk + de*factork
                     factor = 3.0d0*zkyr + 3.0d0*ykzr + 3.0d0*si2*ziyk
     &                           + 3.0d0*si2*yizk - r2inv*yr*part
                     factori = -2.0d0*si2*ykzr + 3.0d0*si2*zkyr
     &                           + yiri2*partik - 2.0d0*enum*yiri2*ziri2
                     factork = -2.0d0*si2*yizr + yrzr + 3.0d0*si2*ziyr
     &                            + ykrk2*partik
                     dtydzi2 = dtdzi2*termy + de*factor
                     dtyidzi2 = dtdzi2*termyi + de*factori
                     dtykdzi2 = dtdzi2*termyk + de*factork
                     factor = 3.0d0*si2*dotp + 6.0d0*zkzr
     &                           + 6.0d0*si2*zizk + 3.0d0*dotk
     &                           - r2inv*(zr*part+si2*dotik)
                     factori = 3.0d0*si2*dotk + si2*zkzr + ziri2*partik
     &                            + enum*(ri2inv-2.0d0*ziri2*ziri2)
                     factork = -r2 + 3.0d0*si2*doti + si2*zizr
     &                            + zrzr + zkrk2*partik
                     dtzdzi2 = dtdzi2*termz + de*factor
                     dtzidzi2 = dtdzi2*termzi + de*factori
                     dtzkdzi2 = dtdzi2*termzk + de*factork
                  end if
c
c     increment diagonal and off-diagonal Hessian elements
c
                  if (i .eq. i1) then
                     hessx(1,i1) = hessx(1,i1) + si1*dtxdxi1 - dtxidxi1
                     hessx(2,i1) = hessx(2,i1) + si1*dtydxi1 - dtyidxi1
                     hessx(3,i1) = hessx(3,i1) + si1*dtzdxi1 - dtzidxi1
                     hessx(1,i2) = hessx(1,i2) + si2*dtxdxi1 + dtxidxi1
                     hessx(2,i2) = hessx(2,i2) + si2*dtydxi1 + dtyidxi1
                     hessx(3,i2) = hessx(3,i2) + si2*dtzdxi1 + dtzidxi1
                     hessx(1,k1) = hessx(1,k1) - sk1*dtxdxi1 - dtxkdxi1
                     hessx(2,k1) = hessx(2,k1) - sk1*dtydxi1 - dtykdxi1
                     hessx(3,k1) = hessx(3,k1) - sk1*dtzdxi1 - dtzkdxi1
                     hessx(1,k2) = hessx(1,k2) - sk2*dtxdxi1 + dtxkdxi1
                     hessx(2,k2) = hessx(2,k2) - sk2*dtydxi1 + dtykdxi1
                     hessx(3,k2) = hessx(3,k2) - sk2*dtzdxi1 + dtzkdxi1
                     hessy(1,i1) = hessy(1,i1) + si1*dtxdyi1 - dtxidyi1
                     hessy(2,i1) = hessy(2,i1) + si1*dtydyi1 - dtyidyi1
                     hessy(3,i1) = hessy(3,i1) + si1*dtzdyi1 - dtzidyi1
                     hessy(1,i2) = hessy(1,i2) + si2*dtxdyi1 + dtxidyi1
                     hessy(3,i2) = hessy(3,i2) + si2*dtzdyi1 + dtzidyi1
                     hessy(2,i2) = hessy(2,i2) + si2*dtydyi1 + dtyidyi1
                     hessy(1,k1) = hessy(1,k1) - sk1*dtxdyi1 - dtxkdyi1
                     hessy(2,k1) = hessy(2,k1) - sk1*dtydyi1 - dtykdyi1
                     hessy(3,k1) = hessy(3,k1) - sk1*dtzdyi1 - dtzkdyi1
                     hessy(1,k2) = hessy(1,k2) - sk2*dtxdyi1 + dtxkdyi1
                     hessy(2,k2) = hessy(2,k2) - sk2*dtydyi1 + dtykdyi1
                     hessy(3,k2) = hessy(3,k2) - sk2*dtzdyi1 + dtzkdyi1
                     hessz(1,i1) = hessz(1,i1) + si1*dtxdzi1 - dtxidzi1
                     hessz(2,i1) = hessz(2,i1) + si1*dtydzi1 - dtyidzi1
                     hessz(3,i1) = hessz(3,i1) + si1*dtzdzi1 - dtzidzi1
                     hessz(1,i2) = hessz(1,i2) + si2*dtxdzi1 + dtxidzi1
                     hessz(2,i2) = hessz(2,i2) + si2*dtydzi1 + dtyidzi1
                     hessz(3,i2) = hessz(3,i2) + si2*dtzdzi1 + dtzidzi1
                     hessz(1,k1) = hessz(1,k1) - sk1*dtxdzi1 - dtxkdzi1
                     hessz(2,k1) = hessz(2,k1) - sk1*dtydzi1 - dtykdzi1
                     hessz(3,k1) = hessz(3,k1) - sk1*dtzdzi1 - dtzkdzi1
                     hessz(1,k2) = hessz(1,k2) - sk2*dtxdzi1 + dtxkdzi1
                     hessz(2,k2) = hessz(2,k2) - sk2*dtydzi1 + dtykdzi1
                     hessz(3,k2) = hessz(3,k2) - sk2*dtzdzi1 + dtzkdzi1
                  else if (i .eq. i2) then
                     hessx(1,i1) = hessx(1,i1) + si1*dtxdxi2 - dtxidxi2
                     hessx(2,i1) = hessx(2,i1) + si1*dtydxi2 - dtyidxi2
                     hessx(3,i1) = hessx(3,i1) + si1*dtzdxi2 - dtzidxi2
                     hessx(1,i2) = hessx(1,i2) + si2*dtxdxi2 + dtxidxi2
                     hessx(2,i2) = hessx(2,i2) + si2*dtydxi2 + dtyidxi2
                     hessx(3,i2) = hessx(3,i2) + si2*dtzdxi2 + dtzidxi2
                     hessx(1,k1) = hessx(1,k1) - sk1*dtxdxi2 - dtxkdxi2
                     hessx(2,k1) = hessx(2,k1) - sk1*dtydxi2 - dtykdxi2
                     hessx(3,k1) = hessx(3,k1) - sk1*dtzdxi2 - dtzkdxi2
                     hessx(1,k2) = hessx(1,k2) - sk2*dtxdxi2 + dtxkdxi2
                     hessx(2,k2) = hessx(2,k2) - sk2*dtydxi2 + dtykdxi2
                     hessx(3,k2) = hessx(3,k2) - sk2*dtzdxi2 + dtzkdxi2
                     hessy(1,i1) = hessy(1,i1) + si1*dtxdyi2 - dtxidyi2
                     hessy(2,i1) = hessy(2,i1) + si1*dtydyi2 - dtyidyi2
                     hessy(3,i1) = hessy(3,i1) + si1*dtzdyi2 - dtzidyi2
                     hessy(1,i2) = hessy(1,i2) + si2*dtxdyi2 + dtxidyi2
                     hessy(2,i2) = hessy(2,i2) + si2*dtydyi2 + dtyidyi2
                     hessy(3,i2) = hessy(3,i2) + si2*dtzdyi2 + dtzidyi2
                     hessy(1,k1) = hessy(1,k1) - sk1*dtxdyi2 - dtxkdyi2
                     hessy(2,k1) = hessy(2,k1) - sk1*dtydyi2 - dtykdyi2
                     hessy(3,k1) = hessy(3,k1) - sk1*dtzdyi2 - dtzkdyi2
                     hessy(1,k2) = hessy(1,k2) - sk2*dtxdyi2 + dtxkdyi2
                     hessy(2,k2) = hessy(2,k2) - sk2*dtydyi2 + dtykdyi2
                     hessy(3,k2) = hessy(3,k2) - sk2*dtzdyi2 + dtzkdyi2
                     hessz(1,i1) = hessz(1,i1) + si1*dtxdzi2 - dtxidzi2
                     hessz(2,i1) = hessz(2,i1) + si1*dtydzi2 - dtyidzi2
                     hessz(3,i1) = hessz(3,i1) + si1*dtzdzi2 - dtzidzi2
                     hessz(1,i2) = hessz(1,i2) + si2*dtxdzi2 + dtxidzi2
                     hessz(2,i2) = hessz(2,i2) + si2*dtydzi2 + dtyidzi2
                     hessz(3,i2) = hessz(3,i2) + si2*dtzdzi2 + dtzidzi2
                     hessz(1,k1) = hessz(1,k1) - sk1*dtxdzi2 - dtxkdzi2
                     hessz(2,k1) = hessz(2,k1) - sk1*dtydzi2 - dtykdzi2
                     hessz(3,k1) = hessz(3,k1) - sk1*dtzdzi2 - dtzkdzi2
                     hessz(1,k2) = hessz(1,k2) - sk2*dtxdzi2 + dtxkdzi2
                     hessz(2,k2) = hessz(2,k2) - sk2*dtydzi2 + dtykdzi2
                     hessz(3,k2) = hessz(3,k2) - sk2*dtzdzi2 + dtzkdzi2
                  end if
c
c     more energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     if (i .eq. i1) then
                        hessx(1,i1) = hessx(1,i1) + si1*dtaperx*dedxi1
     &                         + si1*dtaperx*dedxi1 + si1*si1*d2taperxx
                        hessx(2,i1) = hessx(2,i1) + si1*dtapery*dedxi1
     &                         + si1*dtaperx*dedyi1 + si1*si1*d2taperxy
                        hessx(3,i1) = hessx(3,i1) + si1*dtaperz*dedxi1
     &                         + si1*dtaperx*dedzi1 + si1*si1*d2taperxz
                        hessx(1,i2) = hessx(1,i2) + si2*dtaperx*dedxi1
     &                         + si1*dtaperx*dedxi2 + si2*si1*d2taperxx
                        hessx(2,i2) = hessx(2,i2) + si2*dtapery*dedxi1
     &                         + si1*dtaperx*dedyi2 + si2*si1*d2taperxy
                        hessx(3,i2) = hessx(3,i2) + si2*dtaperz*dedxi1
     &                         + si1*dtaperx*dedzi2 + si2*si1*d2taperxz
                        hessx(1,k1) = hessx(1,k1) - sk1*dtaperx*dedxi1
     &                         + si1*dtaperx*dedxk1 - sk1*si1*d2taperxx
                        hessx(2,k1) = hessx(2,k1) - sk1*dtapery*dedxi1
     &                         + si1*dtaperx*dedyk1 - sk1*si1*d2taperxy
                        hessx(3,k1) = hessx(3,k1) - sk1*dtaperz*dedxi1
     &                         + si1*dtaperx*dedzk1 - sk1*si1*d2taperxz
                        hessx(1,k2) = hessx(1,k2) - sk2*dtaperx*dedxi1
     &                         + si1*dtaperx*dedxk2 - sk2*si1*d2taperxx
                        hessx(2,k2) = hessx(2,k2) - sk2*dtapery*dedxi1
     &                         + si1*dtaperx*dedyk2 - sk2*si1*d2taperxy
                        hessx(3,k2) = hessx(3,k2) - sk2*dtaperz*dedxi1
     &                         + si1*dtaperx*dedzk2 - sk2*si1*d2taperxz
                        hessy(1,i1) = hessy(1,i1) + si1*dtaperx*dedyi1
     &                         + si1*dtapery*dedxi1 + si1*si1*d2taperxy
                        hessy(2,i1) = hessy(2,i1) + si1*dtapery*dedyi1
     &                         + si1*dtapery*dedyi1 + si1*si1*d2taperyy
                        hessy(3,i1) = hessy(3,i1) + si1*dtaperz*dedyi1
     &                         + si1*dtapery*dedzi1 + si1*si1*d2taperyz
                        hessy(1,i2) = hessy(1,i2) + si2*dtaperx*dedyi1
     &                         + si1*dtapery*dedxi2 + si2*si1*d2taperxy
                        hessy(2,i2) = hessy(2,i2) + si2*dtapery*dedyi1
     &                         + si1*dtapery*dedyi2 + si2*si1*d2taperyy
                        hessy(3,i2) = hessy(3,i2) + si2*dtaperz*dedyi1
     &                         + si1*dtapery*dedzi2 + si2*si1*d2taperyz
                        hessy(1,k1) = hessy(1,k1) - sk1*dtaperx*dedyi1
     &                         + si1*dtapery*dedxk1 - sk1*si1*d2taperxy
                        hessy(2,k1) = hessy(2,k1) - sk1*dtapery*dedyi1
     &                         + si1*dtapery*dedyk1 - sk1*si1*d2taperyy
                        hessy(3,k1) = hessy(3,k1) - sk1*dtaperz*dedyi1
     &                         + si1*dtapery*dedzk1 - sk1*si1*d2taperyz
                        hessy(1,k2) = hessy(1,k2) - sk2*dtaperx*dedyi1
     &                         + si1*dtapery*dedxk2 - sk2*si1*d2taperxy
                        hessy(2,k2) = hessy(2,k2) - sk2*dtapery*dedyi1
     &                         + si1*dtapery*dedyk2 - sk2*si1*d2taperyy
                        hessy(3,k2) = hessy(3,k2) - sk2*dtaperz*dedyi1
     &                         + si1*dtapery*dedzk2 - sk2*si1*d2taperyz
                        hessz(1,i1) = hessz(1,i1) + si1*dtaperx*dedzi1
     &                         + si1*dtaperz*dedxi1 + si1*si1*d2taperxz
                        hessz(2,i1) = hessz(2,i1) + si1*dtapery*dedzi1
     &                         + si1*dtaperz*dedyi1 + si1*si1*d2taperyz
                        hessz(3,i1) = hessz(3,i1) + si1*dtaperz*dedzi1
     &                         + si1*dtaperz*dedzi1 + si1*si1*d2taperzz
                        hessz(1,i2) = hessz(1,i2) + si2*dtaperx*dedzi1
     &                         + si1*dtaperz*dedxi2 + si2*si1*d2taperxz
                        hessz(2,i2) = hessz(2,i2) + si2*dtapery*dedzi1
     &                         + si1*dtaperz*dedyi2 + si2*si1*d2taperyz
                        hessz(3,i2) = hessz(3,i2) + si2*dtaperz*dedzi1
     &                         + si1*dtaperz*dedzi2 + si2*si1*d2taperzz
                        hessz(1,k1) = hessz(1,k1) - sk1*dtaperx*dedzi1
     &                         + si1*dtaperz*dedxk1 - sk1*si1*d2taperxz
                        hessz(2,k1) = hessz(2,k1) - sk1*dtapery*dedzi1
     &                         + si1*dtaperz*dedyk1 - sk1*si1*d2taperyz
                        hessz(3,k1) = hessz(3,k1) - sk1*dtaperz*dedzi1
     &                         + si1*dtaperz*dedzk1 - sk1*si1*d2taperzz
                        hessz(1,k2) = hessz(1,k2) - sk2*dtaperx*dedzi1
     &                         + si1*dtaperz*dedxk2 - sk2*si1*d2taperxz
                        hessz(2,k2) = hessz(2,k2) - sk2*dtapery*dedzi1
     &                         + si1*dtaperz*dedyk2 - sk2*si1*d2taperyz
                        hessz(3,k2) = hessz(3,k2) - sk2*dtaperz*dedzi1
     &                         + si1*dtaperz*dedzk2 - sk2*si1*d2taperzz
                     else if (i .eq. i2) then
                        hessx(1,i1) = hessx(1,i1) + si1*dtaperx*dedxi2
     &                         + si2*dtaperx*dedxi1 + si1*si2*d2taperxx
                        hessx(2,i1) = hessx(2,i1) + si1*dtapery*dedxi2
     &                         + si2*dtaperx*dedyi1 + si1*si2*d2taperxy
                        hessx(3,i1) = hessx(3,i1) + si1*dtaperz*dedxi2
     &                         + si2*dtaperx*dedzi1 + si1*si2*d2taperxz
                        hessx(1,i2) = hessx(1,i2) + si2*dtaperx*dedxi2
     &                         + si2*dtaperx*dedxi2 + si2*si2*d2taperxx
                        hessx(2,i2) = hessx(2,i2) + si2*dtapery*dedxi2
     &                         + si2*dtaperx*dedyi2 + si2*si2*d2taperxy
                        hessx(3,i2) = hessx(3,i2) + si2*dtaperz*dedxi2
     &                         + si2*dtaperx*dedzi2 + si2*si2*d2taperxz
                        hessx(1,k1) = hessx(1,k1) - sk1*dtaperx*dedxi2
     &                         + si2*dtaperx*dedxk1 - sk1*si2*d2taperxx
                        hessx(2,k1) = hessx(2,k1) - sk1*dtapery*dedxi2
     &                         + si2*dtaperx*dedyk1 - sk1*si2*d2taperxy
                        hessx(3,k1) = hessx(3,k1) - sk1*dtaperz*dedxi2
     &                         + si2*dtaperx*dedzk1 - sk1*si2*d2taperxz
                        hessx(1,k2) = hessx(1,k2) - sk2*dtaperx*dedxi2
     &                         + si2*dtaperx*dedxk2 - sk2*si2*d2taperxx
                        hessx(2,k2) = hessx(2,k2) - sk2*dtapery*dedxi2
     &                         + si2*dtaperx*dedyk2 - sk2*si2*d2taperxy
                        hessx(3,k2) = hessx(3,k2) - sk2*dtaperz*dedxi2
     &                         + si2*dtaperx*dedzk2 - sk2*si2*d2taperxz
                        hessy(1,i1) = hessy(1,i1) + si1*dtaperx*dedyi2
     &                         + si2*dtapery*dedxi1 + si1*si2*d2taperxy
                        hessy(2,i1) = hessy(2,i1) + si1*dtapery*dedyi2
     &                         + si2*dtapery*dedyi1 + si1*si2*d2taperyy
                        hessy(3,i1) = hessy(3,i1) + si1*dtaperz*dedyi2
     &                         + si2*dtapery*dedzi1 + si1*si2*d2taperyz
                        hessy(1,i2) = hessy(1,i2) + si2*dtaperx*dedyi2
     &                         + si2*dtapery*dedxi2 + si2*si2*d2taperxy
                        hessy(2,i2) = hessy(2,i2) + si2*dtapery*dedyi2
     &                         + si2*dtapery*dedyi2 + si2*si2*d2taperyy
                        hessy(3,i2) = hessy(3,i2) + si2*dtaperz*dedyi2
     &                         + si2*dtapery*dedzi2 + si2*si2*d2taperyz
                        hessy(1,k1) = hessy(1,k1) - sk1*dtaperx*dedyi2
     &                         + si2*dtapery*dedxk1 - sk1*si2*d2taperxy
                        hessy(2,k1) = hessy(2,k1) - sk1*dtapery*dedyi2
     &                         + si2*dtapery*dedyk1 - sk1*si2*d2taperyy
                        hessy(3,k1) = hessy(3,k1) - sk1*dtaperz*dedyi2
     &                         + si2*dtapery*dedzk1 - sk1*si2*d2taperyz
                        hessy(1,k2) = hessy(1,k2) - sk2*dtaperx*dedyi2
     &                         + si2*dtapery*dedxk2 - sk2*si2*d2taperxy
                        hessy(2,k2) = hessy(2,k2) - sk2*dtapery*dedyi2
     &                         + si2*dtapery*dedyk2 - sk2*si2*d2taperyy
                        hessy(3,k2) = hessy(3,k2) - sk2*dtaperz*dedyi2
     &                         + si2*dtapery*dedzk2 - sk2*si2*d2taperyz
                        hessz(1,i1) = hessz(1,i1) + si1*dtaperx*dedzi2
     &                         + si2*dtaperz*dedxi1 + si1*si2*d2taperxz
                        hessz(2,i1) = hessz(2,i1) + si1*dtapery*dedzi2
     &                         + si2*dtaperz*dedyi1 + si1*si2*d2taperyz
                        hessz(3,i1) = hessz(3,i1) + si1*dtaperz*dedzi2
     &                         + si2*dtaperz*dedzi1 + si1*si2*d2taperzz
                        hessz(1,i2) = hessz(1,i2) + si2*dtaperx*dedzi2
     &                         + si2*dtaperz*dedxi2 + si2*si2*d2taperxz
                        hessz(2,i2) = hessz(2,i2) + si2*dtapery*dedzi2
     &                         + si2*dtaperz*dedyi2 + si2*si2*d2taperyz
                        hessz(3,i2) = hessz(3,i2) + si2*dtaperz*dedzi2
     &                         + si2*dtaperz*dedzi2 + si2*si2*d2taperzz
                        hessz(1,k1) = hessz(1,k1) - sk1*dtaperx*dedzi2
     &                         + si2*dtaperz*dedxk1 - sk1*si2*d2taperxz
                        hessz(2,k1) = hessz(2,k1) - sk1*dtapery*dedzi2
     &                         + si2*dtaperz*dedyk1 - sk1*si2*d2taperyz
                        hessz(3,k1) = hessz(3,k1) - sk1*dtaperz*dedzi2
     &                         + si2*dtaperz*dedzk1 - sk1*si2*d2taperzz
                        hessz(1,k2) = hessz(1,k2) - sk2*dtaperx*dedzi2
     &                         + si2*dtaperz*dedxk2 - sk2*si2*d2taperxz
                        hessz(2,k2) = hessz(2,k2) - sk2*dtapery*dedzi2
     &                         + si2*dtaperz*dedyk2 - sk2*si2*d2taperyz
                        hessz(3,k2) = hessz(3,k2) - sk2*dtaperz*dedzi2
     &                         + si2*dtaperz*dedzk2 - sk2*si2*d2taperzz
                     end if
                  end if
               end if
            end do
   20       continue
         end do
   30    continue
      end do
      return
      end
