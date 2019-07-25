c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine echgdpl1  --  charge-dipole energy & derivs  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "echgdpl1" calculates the charge-dipole interaction energy
c     and first derivatives with respect to Cartesian coordinates
c
c
      subroutine echgdpl1
      use sizes
      use atoms
      use bound
      use cell
      use charge
      use chgpot
      use couple
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
      integer i1,k1,k2
      integer, allocatable :: skip(:)
      real*8 e,r2,rk2,rkr3,dotk
      real*8 f,fi,fik,fgrp
      real*8 sk1,sk2
      real*8 xi,yi,zi
      real*8 xk,yk,zk
      real*8 xr,yr,zr
      real*8 xr1,yr1,zr1
      real*8 xr2,yr2,zr2
      real*8 term,term2,term3
      real*8 termx,termy,termz
      real*8 termxk,termyk,termzk
      real*8 dedxi1,dedyi1,dedzi1
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
c     zero out the overall charge-dipole interaction energy
c     and set up the constants for the calculation
c
      ecd = 0.0d0
      do i = 1, n
         decd(1,i) = 0.0d0
         decd(2,i) = 0.0d0
         decd(3,i) = 0.0d0
      end do
      if (nion.eq.0 .or. ndipole.eq.0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (skip(n))
c
c     zero out the list of atoms to be skipped
c
      do i = 1, n
         skip(i) = 0
      end do
c
c     set conversion factor and switching function coefficients
c
      f = electric / (debye * dielec)
      mode = 'CHGDPL'
      call switch (mode)
c
c     get energy and derivs by looping over each charge-dipole pair
c
      do i = 1, nion
         i1 = iion(i)
         skip(i1) = i1
         do k = 1, n12(i1)
            skip(i12(k,i1)) = i1
         end do
         xi = x(i1)
         yi = y(i1)
         zi = z(i1)
         fi = f * pchg(i)
c
c     decide whether to compute the current interaction
c
         do k = 1, ndipole
            k1 = idpl(1,k)
            k2 = idpl(2,k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i1,k1,k2,0,0,0)
            if (proceed)  proceed = (use(i1) .or. use(k1) .or. use(k2))
            if (proceed)  proceed = (skip(k1).ne.i1 .and.
     &                                 skip(k2).ne.i1)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               sk1 = 1.0d0 - sdpl(k)
               sk2 = sdpl(k)
               xk = x(k2) - x(k1)
               yk = y(k2) - y(k1)
               zk = z(k2) - z(k1)
               xr = x(k1) + xk*sk2 - xi
               yr = y(k1) + yk*sk2 - yi
               zr = z(k1) + zk*sk2 - zi
               call image (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  fik = fi * bdpl(k)
                  rk2 = xk*xk + yk*yk + zk*zk
                  rkr3 = sqrt(rk2*r2) * r2
                  dotk = xk*xr + yk*yr + zk*zr
c
c     form the energy and master chain rule term for derivatives
c
                  e = fik * dotk / rkr3
                  term = -fik / rkr3
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     term = term * fgrp
                  end if
c
c     secondary chain rule terms for derivative expressions
c
                  term2 = -3.0d0 * dotk / r2
                  term3 = -dotk / rk2
                  termx = term * (xk+xr*term2)
                  termy = term * (yk+yr*term2)
                  termz = term * (zk+zr*term2)
                  termxk = -term * (xr+xk*term3)
                  termyk = -term * (yr+yk*term3)
                  termzk = -term * (zr+zk*term3)
                  dedxi1 = termx
                  dedyi1 = termy
                  dedzi1 = termz
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
                     dedxi1 = dedxi1*taper - dtaperx
                     dedyi1 = dedyi1*taper - dtapery
                     dedzi1 = dedzi1*taper - dtaperz
                     dedxk1 = dedxk1*taper + sk1*dtaperx
                     dedyk1 = dedyk1*taper + sk1*dtapery
                     dedzk1 = dedzk1*taper + sk1*dtaperz
                     dedxk2 = dedxk2*taper + sk2*dtaperx
                     dedyk2 = dedyk2*taper + sk2*dtapery
                     dedzk2 = dedzk2*taper + sk2*dtaperz
                  end if
c
c     increment the overall energy and derivative expressions
c
                  ecd = ecd + e
                  decd(1,i1) = decd(1,i1) + dedxi1
                  decd(2,i1) = decd(2,i1) + dedyi1
                  decd(3,i1) = decd(3,i1) + dedzi1
                  decd(1,k1) = decd(1,k1) + dedxk1
                  decd(2,k1) = decd(2,k1) + dedyk1
                  decd(3,k1) = decd(3,k1) + dedzk1
                  decd(1,k2) = decd(1,k2) + dedxk2
                  decd(2,k2) = decd(2,k2) + dedyk2
                  decd(3,k2) = decd(3,k2) + dedzk2
c
c     increment the internal virial tensor components
c
                  xr1 = x(k1) - xi
                  yr1 = y(k1) - yi
                  zr1 = z(k1) - zi
                  xr2 = x(k2) - xi
                  yr2 = y(k2) - yi
                  zr2 = z(k2) - zi
                  vxx = xr1*dedxk1 + xr2*dedxk2
                  vyx = yr1*dedxk1 + yr2*dedxk2
                  vzx = zr1*dedxk1 + zr2*dedxk2
                  vyy = yr1*dedyk1 + yr2*dedyk2
                  vzy = zr1*dedyk1 + zr2*dedyk2
                  vzz = zr1*dedzk1 + zr2*dedzk2
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
      do i = 1, nion
         i1 = iion(i)
         skip(i1) = i1
         do k = 1, n12(i1)
            skip(i12(k,i1)) = i1
         end do
         xi = x(i1)
         yi = y(i1)
         zi = z(i1)
         fi = f * pchg(i)
c
c     decide whether to compute the current interaction
c
         do k = 1, ndipole
            k1 = idpl(1,k)
            k2 = idpl(2,k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i1,k1,k2,0,0,0)
            if (proceed)  proceed = (use(i1) .or. use(k1) .or. use(k2))
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
                  xr = x(k1) + xk*sk2 - xi
                  yr = y(k1) + yk*sk2 - yi
                  zr = z(k1) + zk*sk2 - zi
                  call imager (xr,yr,zr,j)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. off2) then
                     fik = fi * bdpl(k)
                     if (use_polymer) then
                        if (r2 .lt. polycut2) then
                           if (skip(k1).eq.i1 .or. skip(k2).ne.i1)
     &                        fik = 0.0d0
                        end if
                     end if
                     rk2 = xk*xk + yk*yk + zk*zk
                     rkr3 = sqrt(rk2*r2) * r2
                     dotk = xk*xr + yk*yr + zk*zr
c
c     form the energy and master chain rule term for derivatives
c
                     e = fik * dotk / rkr3
                     term = -fik / rkr3
c
c     scale the interaction based on its group membership
c
                     if (use_group) then
                        e = e * fgrp
                        term = term * fgrp
                     end if
c
c     secondary chain rule terms for derivative expressions
c
                     term2 = -3.0d0 * dotk / r2
                     term3 = -dotk / rk2
                     termx = term * (xk+xr*term2)
                     termy = term * (yk+yr*term2)
                     termz = term * (zk+zr*term2)
                     termxk = -term * (xr+xk*term3)
                     termyk = -term * (yr+yk*term3)
                     termzk = -term * (zr+zk*term3)
                     dedxi1 = termx
                     dedyi1 = termy
                     dedzi1 = termz
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
                        dedxi1 = dedxi1*taper - dtaperx
                        dedyi1 = dedyi1*taper - dtapery
                        dedzi1 = dedzi1*taper - dtaperz
                        dedxk1 = dedxk1*taper + sk1*dtaperx
                        dedyk1 = dedyk1*taper + sk1*dtapery
                        dedzk1 = dedzk1*taper + sk1*dtaperz
                        dedxk2 = dedxk2*taper + sk2*dtaperx
                        dedyk2 = dedyk2*taper + sk2*dtapery
                        dedzk2 = dedzk2*taper + sk2*dtaperz
                     end if
c
c     increment the overall energy and derivative expressions
c
                     ecd = ecd + e
                     decd(1,i1) = decd(1,i1) + dedxi1
                     decd(2,i1) = decd(2,i1) + dedyi1
                     decd(3,i1) = decd(3,i1) + dedzi1
                     decd(1,k1) = decd(1,k1) + dedxk1
                     decd(2,k1) = decd(2,k1) + dedyk1
                     decd(3,k1) = decd(3,k1) + dedzk1
                     decd(1,k2) = decd(1,k2) + dedxk2
                     decd(2,k2) = decd(2,k2) + dedyk2
                     decd(3,k2) = decd(3,k2) + dedzk2
                  end if
c
c     increment the internal virial tensor components
c
                  xr1 = x(k1) - xi
                  yr1 = y(k1) - yi
                  zr1 = z(k1) - zi
                  xr2 = x(k2) - xi
                  yr2 = y(k2) - yi
                  zr2 = z(k2) - zi
                  vxx = xr1*dedxk1 + xr2*dedxk2
                  vyx = yr1*dedxk1 + yr2*dedxk2
                  vzx = zr1*dedxk1 + zr2*dedxk2
                  vyy = yr1*dedyk1 + yr2*dedyk2
                  vzy = zr1*dedyk1 + zr2*dedyk2
                  vzz = zr1*dedzk1 + zr2*dedzk2
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
               end do
            end if
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (skip)
      return
      end
