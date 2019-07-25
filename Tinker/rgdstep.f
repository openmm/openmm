c
c
c     ###########################################################
c     ##                 COPYRIGHT (C) 2001 by                 ##
c     ##  Andrey Kutepov, Marina A. Vorobieva & Jay W. Ponder  ##
c     ##                  All Rights Reserved                  ##
c     ###########################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine rgdstep  --  rigid body molecular dynamics step  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "rgdstep" performs a single molecular dynamics time step
c     via a rigid body integration algorithm
c
c     literature reference:
c
c     W. Smith, "Hail Euler and Farewell: Rotational Motion in the
c     Laboratory Frame", CCP5 Newsletter, February 2005
c
c     based on an original algorithm developed by Andrey Kutapov
c     and Marina A. Vorobieva, VNIITF, Russian Federal Nuclear
c     Center, Chelyabinsk, Russia, February 2001
c
c
      subroutine rgdstep (istep,dt)
      use sizes
      use atomid
      use atoms
      use bound
      use group
      use iounit
      use rgddyn
      use units
      use virial
      implicit none
      integer i,j,k
      integer istep,size
      integer start,stop
      integer iter,maxiter
      real*8 dt,epot,etot
      real*8 eksum,weigh
      real*8 eps,delta
      real*8 temp,pres
      real*8 xr,yr,zr
      real*8 x2,y2,z2
      real*8 fx,fy,fz
      real*8 fc(3),tc(3)
      real*8 inert(6)
      real*8 rc(3),rcold(3)
      real*8 dfi(3),dfiold(3)
      real*8 vcp(3),wcp(3)
      real*8 ekin(3,3)
      real*8 stress(3,3)
      real*8 arot(3,3)
      real*8, allocatable :: xp(:)
      real*8, allocatable :: yp(:)
      real*8, allocatable :: zp(:)
      real*8, allocatable :: derivs(:,:)
c
c
c     set iteration limit and tolerance for angular momenta
c
      maxiter = 15
      eps = 1.0d-12
c
c     perform dynamic allocation of some local arrays
c
      allocate (xp(n))
      allocate (yp(n))
      allocate (zp(n))
      allocate (derivs(3,n))
c
c     get the energy and atomic forces prior to the step
c
      call gradient (epot,derivs)
c
c     perform the integration step for each rigid body
c
      do i = 1, ngrp
         start = igrp(1,i)
         stop = igrp(2,i)
         size = stop - start + 1
         do j = 1, 3
            rc(j) = 0.0d0
         end do
         do j = start, stop
            k = kgrp(j)
            weigh = mass(k)
            rc(1) = rc(1) + x(k)*weigh
            rc(2) = rc(2) + y(k)*weigh
            rc(3) = rc(3) + z(k)*weigh
         end do
         do j = 1, 3
            rc(j) = rc(j) / grpmass(i)
         end do
c
c     find center of mass offsets only for first step
c
         if (istep .eq. 1) then
            do j = start, stop
               k = kgrp(j)
               xcmo(k) = x(k) - rc(1)
               ycmo(k) = y(k) - rc(2)
               zcmo(k) = z(k) - rc(3)
            end do
         end if
c
c     compute the force and torque components for rigid body
c
         do j = 1, 3
            fc(j) = 0.0d0
            tc(j) = 0.0d0
         end do
         do j = start, stop
            k = kgrp(j)
            xr = x(k) - rc(1)
            yr = y(k) - rc(2)
            zr = z(k) - rc(3)
            fx = -convert * derivs(1,k)
            fy = -convert * derivs(2,k)
            fz = -convert * derivs(3,k)
            fc(1) = fc(1) + fx
            fc(2) = fc(2) + fy
            fc(3) = fc(3) + fz
            tc(1) = tc(1) + yr*fz - zr*fy
            tc(2) = tc(2) + zr*fx - xr*fz
            tc(3) = tc(3) + xr*fy - yr*fx
         end do
c
c     update the translational velocity of the center of mass
c
         do j = 1, 3
            vcp(j) = vcm(j,i) + dt*fc(j)/grpmass(i)
            vc(j,i) = 0.5d0 * (vcm(j,i)+vcp(j))
            vcm(j,i) = vcp(j)
         end do
c
c     update the coordinates of the group center of mass
c
         do j = 1, 3
            rcold(j) = rc(j)
            rc(j) = rc(j) + dt*vcp(j)
         end do
c
c     single atom groups are treated as a separate case
c
         if (size .eq. 1) then
            k = kgrp(igrp(1,i))
            x(k) = rc(1)
            y(k) = rc(2)
            z(k) = rc(3)
            do j = 1, 3
               wcm(j,i) = 0.0d0
               lm(j,i) = 0.0d0
            end do
c
c     get impulse moment in fixed space coordinate system
c
         else
            do j = 1, 3
               lm(j,i) = lm(j,i) + dt*tc(j)
               dfi(j) = dt * wcm(j,i)
               dfiold(j) = dfi(j)
            end do
c
c     use iterative scheme to converge the angular momenta
c
            iter = 0
            delta = 1.0d0
            do while (delta.gt.eps .and. iter.lt.maxiter)
               iter = iter + 1
               call rotrgd (dfi,arot)
c
c     calculate the inertia tensor from rotated coordinates
c
               do j = 1, 6
                  inert(j) = 0.0d0
               end do
               do j = start, stop
                  k = kgrp(j)
                  xr = arot(1,1)*xcmo(k) + arot(1,2)*ycmo(k)
     &                    + arot(1,3)*zcmo(k)
                  yr = arot(2,1)*xcmo(k) + arot(2,2)*ycmo(k)
     &                    + arot(2,3)*zcmo(k)
                  zr = arot(3,1)*xcmo(k) + arot(3,2)*ycmo(k)
     &                    + arot(3,3)*zcmo(k)
                  x2 = xr * xr
                  y2 = yr * yr
                  z2 = zr * zr
                  weigh = mass(k)
                  inert(1) = inert(1) + weigh*(y2+z2)
                  inert(2) = inert(2) - weigh*xr*yr
                  inert(3) = inert(3) - weigh*xr*zr
                  inert(4) = inert(4) + weigh*(x2+z2)
                  inert(5) = inert(5) - weigh*yr*zr
                  inert(6) = inert(6) + weigh*(x2+y2)
                  xp(k) = xr
                  yp(k) = yr
                  zp(k) = zr
               end do
c
c     compute the angular velocity from the relation L=Iw
c
               do j = 1, 3
                  wcp(j) = lm(j,i)
               end do
               if (linear(i)) then
                  call linbody (i,inert,wcp)
               else
                  call cholesky (3,inert,wcp)
               end if
               delta = 0.d0
               do j = 1, 3
                  dfi(j) = 0.5d0 * dt * (wcm(j,i)+wcp(j))
                  delta = delta + abs(dfi(j)-dfiold(j))
                  dfiold(j) = dfi(j)
               end do
            end do
c
c     check to make sure the angular momenta converged
c
            if (delta .gt. eps) then
               write (iout,10)
   10          format (/,' RGDSTEP  --  Angular Momentum Convergence',
     &                    ' Failure')
               call prterr
               call fatal
            end if
c
c     set the final angular velocity and atomic coordinates
c
            do j = 1, 3
               dfi(j) = dt * wcp(j)
            end do
            call rotrgd (dfi,arot)
            do j = start, stop
               k = kgrp(j)
               xr = x(k) - rcold(1)
               yr = y(k) - rcold(2)
               zr = z(k) - rcold(3)
               x(k) = arot(1,1)*xr + arot(1,2)*yr + arot(1,3)*zr + rc(1)
               y(k) = arot(2,1)*xr + arot(2,2)*yr + arot(2,3)*zr + rc(2)
               z(k) = arot(3,1)*xr + arot(3,2)*yr + arot(3,3)*zr + rc(3)
            end do
            do j = 1, 3
               wc(j,i) = 0.5d0 * (wcm(j,i)+wcp(j))
               wcm(j,i) = wcp(j)
            end do
         end if
      end do
c
c     update the distance to center of mass for each atom
c
      do i = 1, n
         xcmo(i) = xp(i)
         ycmo(i) = yp(i)
         zcmo(i) = zp(i)
      end do
c
c     make center of mass correction to virial for rigid body
c
      do i = 1, n
         vir(1,1) = vir(1,1) - xcmo(i)*derivs(1,i)
         vir(2,1) = vir(2,1) - ycmo(i)*derivs(1,i)
         vir(3,1) = vir(3,1) - zcmo(i)*derivs(1,i)
         vir(1,2) = vir(1,2) - xcmo(i)*derivs(2,i)
         vir(2,2) = vir(2,2) - ycmo(i)*derivs(2,i)
         vir(3,2) = vir(3,2) - zcmo(i)*derivs(2,i)
         vir(1,3) = vir(1,3) - xcmo(i)*derivs(3,i)
         vir(2,3) = vir(2,3) - ycmo(i)*derivs(3,i)
         vir(3,3) = vir(3,3) - zcmo(i)*derivs(3,i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xp)
      deallocate (yp)
      deallocate (zp)
      deallocate (derivs)
c
c     make any temperature and pressure corrections
c
      call temper2 (dt,temp)
      call pressure2 (epot,temp)
      call temper (dt,eksum,ekin,temp)
      call pressure (dt,epot,ekin,temp,pres,stress)
c
c     total energy is sum of kinetic and potential energies
c
      etot = eksum + epot
c
c     compute statistics and save trajectory for this step
c
      call mdstat (istep,dt,etot,epot,eksum,temp,pres)
      call mdsave (istep,dt,epot,eksum)
      call mdrest (istep)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine rotrgd  --  rigid dynamics rotation matrix  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "rotrgd" finds the rotation matrix for a rigid body due
c     to a single step of dynamics
c
c
      subroutine rotrgd (dfi,arot)
      implicit none
      real*8 x,xc,xs
      real*8 y,yc,ys
      real*8 z,zc,zs
      real*8 cosine,sine
      real*8 anorm,coterm
      real*8 dfi(3)
      real*8 arot(3,3)
c
c
c     construct rotation matrix from angular distance
c
      anorm = sqrt(dfi(1)**2 + dfi(2)**2 + dfi(3)**2)
      cosine = cos(anorm)
      sine = sin(anorm)
      coterm = 1.0d0 - cosine
      if (anorm .le. 0.0d0)  anorm = 1.0d0
      x = dfi(1) / anorm
      y = dfi(2) / anorm
      z = dfi(3) / anorm
      xc = x * coterm
      yc = y * coterm
      zc = z * coterm
      xs = x * sine
      ys = y * sine
      zs = z * sine
      arot(1,1) = xc*x + cosine
      arot(2,1) = xc*y + zs
      arot(3,1) = xc*z - ys
      arot(1,2) = yc*x - zs
      arot(2,2) = yc*y + cosine
      arot(3,2) = yc*z + xs
      arot(1,3) = zc*x + ys
      arot(2,3) = zc*y - xs
      arot(3,3) = zc*z + cosine
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine linbody  --  angular velocity of linear body  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "linbody" finds the angular velocity of a linear rigid body
c     given the inertia tensor and angular momentum
c
c
      subroutine linbody (i,inert,wcp)
      use sizes
      use atoms
      use group
      implicit none
      integer i,j,k
      real*8 rinv,rmin
      real*8 a11,a12,a22
      real*8 b1,b2,w1,w2
      real*8 wcp(3),rmol(3)
      real*8 r1(3),r2(3),r3(3)
      real*8 inert(6)
c
c
c     construct a normalized vector along the molecular axis
c
      j = kgrp(igrp(1,i))
      k = kgrp(igrp(2,i))
      rmol(1) = x(k) - x(j)
      rmol(2) = y(k) - y(j)
      rmol(3) = z(k) - z(j)
      rinv = 1.0d0 / sqrt(rmol(1)**2+rmol(2)**2+rmol(3)**2)
      do j = 1, 3
         rmol(j) = rmol(j) * rinv
      end do
c
c     find two orthogonal vectors to complete coordinate frame
c
      k = 1
      rmin = abs(rmol(1))
      do j = 2, 3
         if (abs(rmol(j)) .lt. rmin) then
            k = j
            rmin = abs(rmol(j))
         end if
      end do
      do j = 1, 3
         r1(j) = -rmol(k) * rmol(j)
      end do
      r1(k) = 1.0d0 + r1(k)
      rinv = 1.0d0 / sqrt(r1(1)**2+r1(2)**2+r1(3)**2)
      do j = 1, 3
         r1(j) = r1(j) * rinv
      end do
      r2(1) = r1(2)*rmol(3) - r1(3)*rmol(2)
      r2(2) = r1(3)*rmol(1) - r1(1)*rmol(3)
      r2(3) = r1(1)*rmol(2) - r1(2)*rmol(1)
c
c     solve the 2-by-2 linear system for angular velocity
c
      r3(1) = inert(1)*r1(1) + inert(2)*r1(2) + inert(3)*r1(3)
      r3(2) = inert(2)*r1(1) + inert(4)*r1(2) + inert(5)*r1(3)
      r3(3) = inert(3)*r1(1) + inert(5)*r1(2) + inert(6)*r1(3)
      a11 = r1(1)*r3(1) + r1(2)*r3(2) + r1(3)*r3(3)
      r3(1) = inert(1)*r2(1) + inert(2)*r2(2) + inert(3)*r2(3)
      r3(2) = inert(2)*r2(1) + inert(4)*r2(2) + inert(5)*r2(3)
      r3(3) = inert(3)*r2(1) + inert(5)*r2(2) + inert(6)*r2(3)
      a12 = r1(1)*r3(1) + r1(2)*r3(2) + r1(3)*r3(3)
      a22 = r2(1)*r3(1) + r2(2)*r3(2) + r2(3)*r3(3)
      b1 = r1(1)*wcp(1) + r1(2)*wcp(2) + r1(3)*wcp(3)
      b2 = r2(1)*wcp(1) + r2(2)*wcp(2) + r2(3)*wcp(3)
      w1 = (a12*b2-a22*b1) / (a12*a12-a11*a22)
      w2 = (b2-a12*w1) / a22
      do j = 1, 3
         wcp(j) = w1*r1(j) + w2*r2(j)
      end do
      return
      end
