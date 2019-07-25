c
c
c     ###############################################################
c     ##  COPYRIGHT (C) 2014 by Alex Albaugh & Jay William Ponder  ##
c     ##                    All Rights Reserved                    ##
c     ###############################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine kinetic  --  compute kinetic energy components  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "kinetic" computes the total kinetic energy and kinetic energy
c     contributions to the pressure tensor by summing over velocities
c
c
      subroutine kinetic (eksum,ekin,temp)
      use sizes
      use atomid
      use atoms
      use bath
      use group
      use mdstuf
      use moldyn
      use rgddyn
      use units
      use usage
      implicit none
      integer i,j,k
      integer start,stop
      real*8 eksum,temp
      real*8 weigh
      real*8 term,value
      real*8 xr,yr,zr
      real*8 x2,y2,z2
      real*8 xcm,ycm,zcm
      real*8 ekin(3,3)
      real*8 inert(3,3)
c
c
c     zero out the temperature and kinetic energy components
c
      temp = 0.0d0
      eksum = 0.0d0
      do i = 1, 3
         do j = 1, 3
            ekin(j,i) = 0.0d0
         end do
      end do
c
c     get the total kinetic energy and tensor for atomic sites
c
      if (integrate .ne. 'RIGIDBODY') then
         do i = 1, n
            if (use(i)) then
               term = 0.5d0 * mass(i) / convert
               do j = 1, 3
                  do k = 1, 3
                     value = term * v(j,i) * v(k,i)
                     ekin(k,j) = ekin(k,j) + value
                  end do
               end do
            end if
         end do
         eksum = ekin(1,1) + ekin(2,2) + ekin(3,3)
c
c     get the total kinetic energy and tensor for rigid bodies
c
      else
         do i = 1, ngrp
            start = igrp(1,i)
            stop = igrp(2,i)
            xcm = 0.0d0
            ycm = 0.0d0
            zcm = 0.0d0
            do j = start, stop
               k = kgrp(j)
               weigh = mass(k)
               xcm = xcm + x(k)*weigh
               ycm = ycm + y(k)*weigh
               zcm = zcm + z(k)*weigh
            end do
            xcm = xcm / grpmass(i)
            ycm = ycm / grpmass(i)
            zcm = zcm / grpmass(i)
c
c     find the inertial tensor relative to the center of mass
c
            do j = 1, 3
               do k = 1, 3
                  inert(k,j) = 0.0d0
               end do
            end do
            do j = start, stop
               k = kgrp(j)
               xr = x(k) - xcm
               yr = y(k) - ycm
               zr = z(k) - zcm
               x2 = xr * xr
               y2 = yr * yr
               z2 = zr * zr
               weigh = mass(k)
               inert(1,1) = inert(1,1) + weigh*(y2+z2)
               inert(2,1) = inert(2,1) - weigh*xr*yr
               inert(3,1) = inert(3,1) - weigh*xr*zr
               inert(2,2) = inert(2,2) + weigh*(x2+z2)
               inert(3,2) = inert(3,2) - weigh*yr*zr
               inert(3,3) = inert(3,3) + weigh*(x2+y2)
            end do
            inert(1,2) = inert(2,1)
            inert(1,3) = inert(3,1)
            inert(2,3) = inert(3,2)
c
c     increment the kinetic energy due to translational motion
c
            term = 0.5d0 * grpmass(i) / convert
            do j = 1, 3
               do k = 1, 3
                  value = term * vc(j,i) * vc(k,i)
                  ekin(k,j) = ekin(k,j) + value
                  if (j .eq. k)  eksum = eksum + value
               end do
            end do
c
c     increment the kinetic energy due to rotational motion
c
            term = 0.5d0 / convert
            do j = 1, 3
               do k = 1, 3
                  value = term * inert(k,j) * wc(j,i) * wc(k,i)
                  eksum = eksum + value
               end do
            end do
         end do
      end if
c
c     set the instantaneous temperature from total kinetic energy
c
      temp = 2.0d0 * eksum / (dble(nfree) * gasconst)
c
c     get the kinetic energy for Bussi-Parrinello barostat
c
      if (isobaric .and. barostat.eq.'BUSSI') then
         term = dble(nfree) * gasconst * kelvin * taupres * taupres
         value = 0.5d0 * term * eta * eta
         do j = 1, 3
            ekin(j,j) = ekin(j,j) + value/3.0d0
         end do
         eksum = eksum + value
      end if
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine kinaux -- compute iEL dipole kinetic energy  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "kinaux" computes the total kinetic energy and temperature
c     for auxiliary dipole variables used in iEL polarization
c
c
      subroutine kinaux (temp_aux,temp_auxp)
      use atoms
      use ielscf
      implicit none
      integer i,j,k
      real*8 term
      real*8 vj,vjp
      real*8 vk,vkp
      real*8 temp_aux
      real*8 temp_auxp
      real*8 eksum_aux
      real*8 eksum_auxp
      real*8 ekaux(3,3)
      real*8 ekauxp(3,3)
c
c
c     zero out the temperature and kinetic energy components
c
      temp_aux = 0.0d0
      temp_auxp = 0.0d0
      do i = 1, 3
         do j = 1, 3
            ekaux(j,i) = 0.0d0
            ekauxp(j,i) = 0.0d0
         end do
      end do
c
c     get the kinetic energy tensor for auxiliary variables
c
      do i = 1, n
         term = 0.5d0
         do j = 1, 3
            vj = vaux(j,i)
            vjp = vpaux(j,i)
            do k = 1, 3
               vk = vaux(k,i)
               vkp = vpaux(k,i)
               ekaux(k,j) = ekaux(k,j) + term*vj*vk
               ekauxp(k,j) = ekauxp(k,j) + term*vjp*vkp
            end do
         end do
      end do
c
c     find the total kinetic energy and auxiliary temperatures
c
      eksum_aux = ekaux(1,1) + ekaux(2,2) + ekaux(3,3)
      eksum_auxp = ekauxp(1,1) + ekauxp(2,2) + ekauxp(3,3)
      if (nfree_aux .ne. 0) then
         temp_aux = 2.0d0 * eksum_aux / dble(nfree_aux)
         temp_auxp = 2.0d0 * eksum_auxp / dble(nfree_aux)
      end if
      return
      end
