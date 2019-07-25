c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine edipole1  --  dipole-dipole energy & derivs  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "edipole1" calculates the dipole-dipole interaction energy
c     and first derivatives with respect to Cartesian coordinates
c
c
      subroutine edipole1
      use sizes
      use atoms
      use bound
      use cell
      use chgpot
      use deriv
      use dipole
      use energi
      use group
      use inter
      use molcul
      use shunt
      use units
      use usage
      use virial
      implicit none
      integer i,j,k
      integer i1,i2,k1,k2
      real*8 xi,yi,zi
      real*8 xk,yk,zk
      real*8 xq,yq,zq
      real*8 xr,yr,zr
      real*8 xq1,yq1,zq1
      real*8 xq2,yq2,zq2
      real*8 f,fi,fik,fgrp
      real*8 e,r2,ri2,rk2,rirkr3
      real*8 doti,dotk,dotp
      real*8 si1,si2,sk1,sk2
      real*8 de,dedr,dedrirk
      real*8 deddoti,deddotk,deddotp
      real*8 termx,termy,termz
      real*8 dedrirkri2,dedrirkrk2
      real*8 termxi,termyi,termzi
      real*8 termxk,termyk,termzk
      real*8 dedxi1,dedyi1,dedzi1
      real*8 dedxi2,dedyi2,dedzi2
      real*8 dedxk1,dedyk1,dedzk1
      real*8 dedxk2,dedyk2,dedzk2
      real*8 r,r3,r4,r5,taper,dtaper
      real*8 dtaperx,dtapery,dtaperz
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      logical proceed
      character*6 mode
c
c
c     zero out the overall dipole interaction energy and derivs,
c     then set up the constants for the calculation
c
      ed = 0.0d0
      do i = 1, n
         ded(1,i) = 0.0d0
         ded(2,i) = 0.0d0
         ded(3,i) = 0.0d0
      end do
      if (ndipole .eq. 0)  return
c
c     set conversion factor and switching function coefficients
c
      f = electric / (debye**2 * dielec)
      mode = 'DIPOLE'
      call switch (mode)
c
c     compute the dipole interaction energy and first derivatives
c
      do i = 1, ndipole-1
         i1 = idpl(1,i)
         i2 = idpl(2,i)
         si1 = 1.0d0 - sdpl(i)
         si2 = sdpl(i)
         xi = x(i2) - x(i1)
         yi = y(i2) - y(i1)
         zi = z(i2) - z(i1)
         if (use_polymer)  call imager (xi,yi,zi,-1)
         ri2 = xi*xi + yi*yi + zi*zi
         xq = x(i1) + xi*si2
         yq = y(i1) + yi*si2
         zq = z(i1) + zi*si2
         fi = f * bdpl(i)
c
c     decide whether to compute the current interaction
c
         do k = i+1, ndipole
            k1 = idpl(1,k)
            k2 = idpl(2,k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i1,i2,k1,k2,0,0)
            if (proceed)  proceed = (use(i1) .or. use(i2) .or.
     &                                 use(k1) .or. use(k2))
            if (proceed)  proceed = (k1.ne.i1 .and. k1.ne.i2 .and.
     &                                 k2.ne.i1 .and. k2.ne.i2)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               sk1 = 1.0d0 - sdpl(k)
               sk2 = sdpl(k)
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
                  fik = fi * bdpl(k)
c
c     form the energy and master chain rule term for derivatives
c
                  e = fik * (dotp-3.0d0*doti*dotk/r2) / rirkr3
                  de = -fik / (rirkr3*r2)
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     de = de * fgrp
                  end if
c
c     secondary chain rule terms for derivative expressions
c
                  deddotp = -de * r2
                  deddoti = de * 3.0d0*dotk
                  deddotk = de * 3.0d0*doti
                  dedr = de * (3.0d0*dotp-15.0d0*doti*dotk/r2)
                  dedrirk = -e
                  dedrirkri2 = dedrirk / ri2
                  dedrirkrk2 = dedrirk / rk2
c
c     more chain rule terms for derivative expressions
c
                  termx = dedr*xr + deddoti*xi + deddotk*xk
                  termy = dedr*yr + deddoti*yi + deddotk*yk
                  termz = dedr*zr + deddoti*zi + deddotk*zk
                  termxi = dedrirkri2*xi + deddotp*xk + deddoti*xr
                  termyi = dedrirkri2*yi + deddotp*yk + deddoti*yr
                  termzi = dedrirkri2*zi + deddotp*zk + deddoti*zr
                  termxk = dedrirkrk2*xk + deddotp*xi + deddotk*xr
                  termyk = dedrirkrk2*yk + deddotp*yi + deddotk*yr
                  termzk = dedrirkrk2*zk + deddotp*zi + deddotk*zr
c
c     finally, the individual first derivative components
c
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
c
c     use energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     r = sqrt(r2)
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     taper = c5*r5 + c4*r4 + c3*r3
     &                          + c2*r2 + c1*r + c0
                     dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                           + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                     dtaper = dtaper * e/r
                     dtaperx = xr * dtaper
                     dtapery = yr * dtaper
                     dtaperz = zr * dtaper
                     e = e * taper
                     dedxi1 = dedxi1*taper + si1*dtaperx
                     dedyi1 = dedyi1*taper + si1*dtapery
                     dedzi1 = dedzi1*taper + si1*dtaperz
                     dedxi2 = dedxi2*taper + si2*dtaperx
                     dedyi2 = dedyi2*taper + si2*dtapery
                     dedzi2 = dedzi2*taper + si2*dtaperz
                     dedxk1 = dedxk1*taper - sk1*dtaperx
                     dedyk1 = dedyk1*taper - sk1*dtapery
                     dedzk1 = dedzk1*taper - sk1*dtaperz
                     dedxk2 = dedxk2*taper - sk2*dtaperx
                     dedyk2 = dedyk2*taper - sk2*dtapery
                     dedzk2 = dedzk2*taper - sk2*dtaperz
                  end if
c
c     increment the overall energy and derivative expressions
c
                  ed = ed + e
                  ded(1,i1) = ded(1,i1) + dedxi1
                  ded(2,i1) = ded(2,i1) + dedyi1
                  ded(3,i1) = ded(3,i1) + dedzi1
                  ded(1,i2) = ded(1,i2) + dedxi2
                  ded(2,i2) = ded(2,i2) + dedyi2
                  ded(3,i2) = ded(3,i2) + dedzi2
                  ded(1,k1) = ded(1,k1) + dedxk1
                  ded(2,k1) = ded(2,k1) + dedyk1
                  ded(3,k1) = ded(3,k1) + dedzk1
                  ded(1,k2) = ded(1,k2) + dedxk2
                  ded(2,k2) = ded(2,k2) + dedyk2
                  ded(3,k2) = ded(3,k2) + dedzk2
c
c     increment the internal virial tensor components
c
                  xq1 = x(k1) - xq
                  yq1 = y(k1) - yq
                  zq1 = z(k1) - zq
                  xq2 = x(k2) - xq
                  yq2 = y(k2) - yq
                  zq2 = z(k2) - zq
                  vxx = xq1*dedxk1 + xq2*dedxk2
                  vyx = 0.5d0 * (yq1*dedxk1 + yq2*dedxk2
     &                              + xq1*dedyk1 + xq2*dedyk2)
                  vzx = 0.5d0 * (zq1*dedxk1 + zq2*dedxk2
     &                              + xq1*dedzk1 + xq2*dedzk2)
                  vyy = yq1*dedyk1 + yq2*dedyk2
                  vzy = 0.5d0 * (zq1*dedyk1 + zq2*dedyk2
     &                              + yq1*dedzk1 + yq2*dedzk2)
                  vzz = zq1*dedzk1 + zq2*dedzk2
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
                  if (molcule(i1) .ne. molcule(k1)) then
                     einter = einter + e
                  end if
               end if
            end if
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do i = 1, ndipole
         i1 = idpl(1,i)
         i2 = idpl(2,i)
         si1 = 1.0d0 - sdpl(i)
         si2 = sdpl(i)
         xi = x(i2) - x(i1)
         yi = y(i2) - y(i1)
         zi = z(i2) - z(i1)
         if (use_polymer)  call imager (xi,yi,zi,-1)
         ri2 = xi*xi + yi*yi + zi*zi
         xq = x(i1) + xi*si2
         yq = y(i1) + yi*si2
         zq = z(i1) + zi*si2
         fi = f * bdpl(i)
c
c     decide whether to compute the current interaction
c
         do k = i, ndipole
            k1 = idpl(1,k)
            k2 = idpl(2,k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i1,i2,k1,k2,0,0)
            if (proceed)  proceed = (use(i1) .or. use(i2) .or.
     &                                 use(k1) .or. use(k2))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               sk1 = 1.0d0 - sdpl(k)
               sk2 = sdpl(k)
               do j = 1, ncell
                  xk = x(k2) - x(k1)
                  yk = y(k2) - y(k1)
                  zk = z(k2) - z(k1)
                  if (use_polymer)  call imager (xk,yk,zk,-1)
                  xr = xq - x(k1) - xk*sk2
                  yr = yq - y(k1) - yk*sk2
                  zr = zq - z(k1) - zk*sk2
                  call imager (xr,yr,zr,j)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. off2) then
                     rk2 = xk*xk + yk*yk + zk*zk
                     rirkr3 = sqrt(ri2*rk2*r2) * r2
                     dotp = xi*xk + yi*yk + zi*zk
                     doti = xi*xr + yi*yr + zi*zr
                     dotk = xk*xr + yk*yr + zk*zr
                     fik = fi * bdpl(k)
                     if (use_polymer) then
                        if (r2 .lt. polycut2) then
                           if (k1.eq.i1 .or. k1.eq.i2 .or.
     &                         k2.eq.i1 .or. k2.eq.i2)  fik = 0.0d0
                        end if
                     end if
c
c     form the energy and master chain rule term for derivatives
c
                     e = fik * (dotp-3.0d0*doti*dotk/r2) / rirkr3
                     de = -fik / (rirkr3*r2)
c
c     scale the interaction based on its group membership
c
                     if (use_group) then
                        e = e * fgrp
                        de = de * fgrp
                     end if
c
c     secondary chain rule terms for derivative expressions
c
                     deddotp = -de * r2
                     deddoti = de * 3.0d0*dotk
                     deddotk = de * 3.0d0*doti
                     dedr = de * (3.0d0*dotp-15.0d0*doti*dotk/r2)
                     dedrirk = -e
                     dedrirkri2 = dedrirk / ri2
                     dedrirkrk2 = dedrirk / rk2
c
c     more chain rule terms for derivative expressions
c
                     termx = dedr*xr + deddoti*xi + deddotk*xk
                     termy = dedr*yr + deddoti*yi + deddotk*yk
                     termz = dedr*zr + deddoti*zi + deddotk*zk
                     termxi = dedrirkri2*xi + deddotp*xk + deddoti*xr
                     termyi = dedrirkri2*yi + deddotp*yk + deddoti*yr
                     termzi = dedrirkri2*zi + deddotp*zk + deddoti*zr
                     termxk = dedrirkrk2*xk + deddotp*xi + deddotk*xr
                     termyk = dedrirkrk2*yk + deddotp*yi + deddotk*yr
                     termzk = dedrirkrk2*zk + deddotp*zi + deddotk*zr
c
c     finally, the individual first derivative components
c
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
c
c     use energy switching if near the cutoff distance
c
                     if (r2 .gt. cut2) then
                        r = sqrt(r2)
                        r3 = r2 * r
                        r4 = r2 * r2
                        r5 = r2 * r3
                        taper = c5*r5 + c4*r4 + c3*r3
     &                             + c2*r2 + c1*r + c0
                        dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                              + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                        dtaper = dtaper * e/r
                        dtaperx = xr * dtaper
                        dtapery = yr * dtaper
                        dtaperz = zr * dtaper
                        e = e * taper
                        dedxi1 = dedxi1*taper + si1*dtaperx
                        dedyi1 = dedyi1*taper + si1*dtapery
                        dedzi1 = dedzi1*taper + si1*dtaperz
                        dedxi2 = dedxi2*taper + si2*dtaperx
                        dedyi2 = dedyi2*taper + si2*dtapery
                        dedzi2 = dedzi2*taper + si2*dtaperz
                        dedxk1 = dedxk1*taper - sk1*dtaperx
                        dedyk1 = dedyk1*taper - sk1*dtapery
                        dedzk1 = dedzk1*taper - sk1*dtaperz
                        dedxk2 = dedxk2*taper - sk2*dtaperx
                        dedyk2 = dedyk2*taper - sk2*dtapery
                        dedzk2 = dedzk2*taper - sk2*dtaperz
                     end if
c
c     increment the overall energy and derivative expressions
c
                     if (i .eq. k)  e = 0.5d0 * e
                     ed = ed + e
                     ded(1,i1) = ded(1,i1) + dedxi1
                     ded(2,i1) = ded(2,i1) + dedyi1
                     ded(3,i1) = ded(3,i1) + dedzi1
                     ded(1,i2) = ded(1,i2) + dedxi2
                     ded(2,i2) = ded(2,i2) + dedyi2
                     ded(3,i2) = ded(3,i2) + dedzi2
                     if (i .ne. k) then
                        ded(1,k1) = ded(1,k1) + dedxk1
                        ded(2,k1) = ded(2,k1) + dedyk1
                        ded(3,k1) = ded(3,k1) + dedzk1
                        ded(1,k2) = ded(1,k2) + dedxk2
                        ded(2,k2) = ded(2,k2) + dedyk2
                        ded(3,k2) = ded(3,k2) + dedzk2
                     end if
c
c     increment the internal virial tensor components
c
                     xq1 = x(k1) - xq
                     yq1 = y(k1) - yq
                     zq1 = z(k1) - zq
                     xq2 = x(k2) - xq
                     yq2 = y(k2) - yq
                     zq2 = z(k2) - zq
                     vxx = xq1*dedxk1 + xq2*dedxk2
                     vyx = 0.5d0 * (yq1*dedxk1 + yq2*dedxk2
     &                                 + xq1*dedyk1 + xq2*dedyk2)
                     vzx = 0.5d0 * (zq1*dedxk1 + zq2*dedxk2
     &                                 + xq1*dedzk1 + xq2*dedzk2)
                     vyy = yq1*dedyk1 + yq2*dedyk2
                     vzy = 0.5d0 * (zq1*dedyk1 + zq2*dedyk2
     &                                 + yq1*dedzk1 + yq2*dedzk2)
                     vzz = zq1*dedzk1 + zq2*dedzk2
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
                     einter = einter + e
                  end if
               end do
            end if
         end do
      end do
      return
      end
