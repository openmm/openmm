c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2001 by Anders Carlsson & Jay William Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine emetal1  --  ligand field energy & derivatives  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "emetal1" calculates the transition metal ligand field energy
c     and its first derivatives with respect to Cartesian coordinates
c
c
      subroutine emetal1
      use sizes
      use atomid
      use atoms
      use couple
      use deriv
      use energi
      use kchrge
      implicit none
      integer maxneigh
      parameter (maxneigh=10)
      integer i,j,k,jj,k0
      integer nneigh,ineigh,jneigh
      integer neighnum(maxneigh)
      real*8 e,e0,elf0,esqtet
      real*8 rcut,r2,rback2
      real*8 cij,argument,dot
      real*8 xmet,ymet,zmet
      real*8 alpha,beta
      real*8 kappa,r0cu
      real*8 demet(3)
      real*8 rback(3,maxneigh)
      real*8 dedrback(3,maxneigh)
      real*8 facback(maxneigh)
      real*8 dfacback(3,maxneigh)
      real*8 dftback(3,maxneigh)
      real*8 detback(3,maxneigh)
      real*8 de(3,maxneigh)
      real*8 delfdh(maxneigh)
      real*8 xneigh(maxneigh)
      real*8 yneigh(maxneigh)
      real*8 zneigh(maxneigh)
      real*8 rneigh(maxneigh)
      real*8 expfac(maxneigh)
      real*8 delfrad(maxneigh)
      real*8 angfac(maxneigh,maxneigh)
      real*8 dangfac(maxneigh,maxneigh)
      real*8 cosmat(maxneigh,maxneigh)
c
c
c     zero out metal ligand field energy and first derivatives
c
      elf = 0.0d0
      do i = 1, n
         delf(1,i) = 0.0d0
         delf(2,i) = 0.0d0
         delf(3,i) = 0.0d0
      end do
c
c     begin ligand field calculation; for now, only for Cu+2
c
      do i = 1, n
         if (atomic(i).ne.29 .or. chg(type(i)).ne.2.0d0)  goto 30
         nneigh = 0
         xmet = x(i)
         ymet = y(i)
         zmet = z(i)
         do j = 1, n
            if (j .eq. i)  goto 10
            if (atomic(j).ne.7 .and. atomic(j).ne.8)  goto 10
c
c     next are standard distance, decay factor and splitting energy
c
            r0cu = 2.06d0
            kappa = 1.0d0
            esqtet = 1.64d0
c
c     semiclassical method obtains only 65% of sq-tet difference
c
            elf0 = esqtet * 1.78d0/0.78d0
            elf0 = elf0 * 0.65d0
            e0 = elf0 * 23.05d0/2.608d0
            rcut = 2.50d0
            alpha = 0.0d0
            beta = -1.0d0
            r2 = (x(i)-x(j))**2 + (y(i)-y(j))**2 + (z(i)-z(j))**2
            if (r2 .gt. rcut**2)  goto 10
            nneigh = nneigh + 1
            xneigh(nneigh) = x(j)
            yneigh(nneigh) = y(j)
            zneigh(nneigh) = z(j)
            neighnum(nneigh) = j
            if (n12(j) .le. 0)  call fatal
c
c     set average of back neighbors relative to j,
c     notice that it is defined as a difference
c
            rback(1,nneigh) = 0.0d0
            rback(2,nneigh) = 0.0d0
            rback(3,nneigh) = 0.0d0
            do k0 = 1, n12(j)
               k = i12(k0,j)
               rback(1,nneigh) = rback(1,nneigh)
     &                              + (x(k)-x(j))/dble(n12(j))
               rback(2,nneigh) = rback(2,nneigh)
     &                              + (y(k)-y(j))/dble(n12(j))
               rback(3,nneigh) = rback(3,nneigh)
     &                              + (z(k)-z(j))/dble(n12(j))
            end do
            facback(nneigh) = 0.0d0
            dot = rback(1,nneigh)*(x(i)-x(j))
     &               + rback(2,nneigh)*(y(i)-y(j))
     &               + rback(3,nneigh)*(z(i)-z(j))
            rback2 = rback(1,nneigh)**2 + rback(2,nneigh)**2
     &                  + rback(3,nneigh)**2
            facback(nneigh) = (alpha*sqrt(rback2)*sqrt(r2)+beta*dot)**2
     &                                      / (rback2*r2)
c
c     dfacback is derivative of back factor with respect to center
c     of gravity of back atoms; detback is corresponding energy
c
            dfacback(1,nneigh) = 2.0d0*(alpha*sqrt(r2)*rback(1,nneigh)
     &                              /sqrt(rback2)+beta*(x(i)-x(j)))
     &                   *(alpha*sqrt(rback2*r2)+beta*dot)/(rback2*r2)
     &                   -2.0d0*facback(nneigh)*rback(1,nneigh)/rback2
            dfacback(2,nneigh) = 2.0d0*(alpha*sqrt(r2)*rback(2,nneigh)
     &                              /sqrt(rback2)+beta*(y(i)-y(j)))
     &                   *(alpha*sqrt(rback2*r2)+beta*dot)/(rback2*r2)
     &                   -2.0d0*facback(nneigh)*rback(2,nneigh)/rback2
            dfacback(3,nneigh) = 2.0d0*(alpha*sqrt(r2)*rback(3,nneigh)
     &                              /sqrt(rback2)+beta*(z(i)-z(j)))
     &                   *(alpha*sqrt(rback2*r2)+beta*dot)/(rback2*r2)
     &                   -2.0d0*facback(nneigh)*rback(3,nneigh)/rback2
            dftback(1,nneigh) = 2.0d0*(alpha*sqrt(rback2)*(x(i)-x(j))
     &                             /sqrt(r2)+beta*(rback(1,nneigh)))
     &                   *(alpha*sqrt(rback2*r2)+beta*dot)/(rback2*r2)
     &                          -2.0d0*facback(nneigh)*(x(i)-x(j))/r2
            dftback(2,nneigh) = 2.0d0*(alpha*sqrt(rback2)*(y(i)-y(j))
     &                             /sqrt(r2)+beta*(rback(2,nneigh)))
     &                   *(alpha*sqrt(rback2*r2)+beta*dot)/(rback2*r2)
     &                          -2.0d0*facback(nneigh)*(y(i)-y(j))/r2
            dftback(3,nneigh) = 2.0d0*(alpha*sqrt(rback2)*(z(i)-z(j))
     &                             /sqrt(r2)+beta*(rback(3,nneigh)))
     &                   *(alpha*sqrt(rback2*r2)+beta*dot)/(rback2*r2)
     &                          -2.0d0*facback(nneigh)*(z(i)-z(j))/r2
   10       continue
         end do
c
c     calculate the energy and derivatives for current interaction
c
         do ineigh = 1, nneigh
            rneigh(ineigh) = (xneigh(ineigh)-xmet)**2
     &                          + (yneigh(ineigh)-ymet)**2
     &                          + (zneigh(ineigh)-zmet)**2
            rneigh(ineigh) = sqrt(rneigh(ineigh))
            expfac(ineigh) = exp(-kappa*(rneigh(ineigh)-r0cu))
            jj = neighnum(ineigh)
            if (atomic(jj).eq.8)  expfac(ineigh) = 0.4d0*expfac(ineigh)
         end do
         do ineigh = 1, nneigh
            jj = neighnum(ineigh)
         end do
         do ineigh = 1, nneigh
            do jneigh = 1, nneigh
               dot = (xneigh(ineigh)-xmet)*(xneigh(jneigh)-xmet)
     &                  + (yneigh(ineigh)-ymet)*(yneigh(jneigh)-ymet)
     &                  + (zneigh(ineigh)-zmet)*(zneigh(jneigh)-zmet)
               cij = dot / (rneigh(ineigh)*rneigh(jneigh))
               cosmat(ineigh,jneigh) = cij
               angfac(ineigh,jneigh) = 2.25d0*cij**4 - 1.5d0*cij**2
     &                                          + 0.005d0
               dangfac(ineigh,jneigh) = 9.0d0*cij**3 - 3.0d0*cij
            end do
         end do
         argument = 0.0d0
         do ineigh = 1, nneigh
            do jneigh = 1, nneigh
               argument = argument + expfac(ineigh)*expfac(jneigh)
     &                                   *angfac(ineigh,jneigh)
     &                              *facback(ineigh)*facback(jneigh)
            end do
         end do
         e = 0.0d0
         do ineigh = 1, nneigh
            de(1,ineigh) = 0.0d0
            de(2,ineigh) = 0.0d0
            de(3,ineigh) = 0.0d0
         end do
         if (argument .le. 0)  goto 20
         if (argument .gt. 0)  e = -e0 * sqrt(argument)
c
c     set up radial derivatives of energy
c
         do ineigh = 1, nneigh
            delfrad(ineigh) = 0.0d0
            do jneigh = 1,nneigh
               delfrad(ineigh) = delfrad(ineigh) + expfac(jneigh)
     &                                         *angfac(ineigh,jneigh)
     &                                            *facback(jneigh)
            end do
            delfdh(ineigh) = delfrad(ineigh) * (e0/e)
c
c     note two minus signs above cancel
c
            delfrad(ineigh) = -delfrad(ineigh) * (e0**2/e)
     &                           * kappa*expfac(ineigh)*facback(ineigh)
c
c     note cancelling factors of two from square root and product
c
         end do
c
c     below does angular derivatives
c
         do ineigh = 1, nneigh
            de(1,ineigh) = 0.0d0
            de(2,ineigh) = 0.0d0
            de(3,ineigh) = 0.0d0
            do jneigh = 1, nneigh
               de(1,ineigh) = de(1,ineigh) +
     &    expfac(jneigh)*facback(jneigh)*dangfac(ineigh,jneigh)*
     &    ((xneigh(jneigh)-xmet)/(rneigh(ineigh)*rneigh(jneigh))-
     &    (xneigh(ineigh)-xmet)*cosmat(ineigh,jneigh)/rneigh(ineigh)**2)
               de(2,ineigh) = de(2,ineigh) +
     &    expfac(jneigh)*facback(jneigh)*dangfac(ineigh,jneigh)*
     &    ((yneigh(jneigh)-ymet)/(rneigh(ineigh)*rneigh(jneigh))-
     &    (yneigh(ineigh)-ymet)*cosmat(ineigh,jneigh)/rneigh(ineigh)**2)
               de(3,ineigh) = de(3,ineigh) +
     &    expfac(jneigh)*facback(jneigh)*dangfac(ineigh,jneigh)*
     &    ((zneigh(jneigh)-zmet)/(rneigh(ineigh)*rneigh(jneigh))-
     &    (zneigh(ineigh)-zmet)*cosmat(ineigh,jneigh)/rneigh(ineigh)**2)
            end do
            de(1,ineigh) = de(1,ineigh)*e0*e0*expfac(ineigh)
     &                            *facback(ineigh)/e
            de(2,ineigh) = de(2,ineigh)*e0*e0*expfac(ineigh)
     &                            *facback(ineigh)/e
            de(3,ineigh) = de(3,ineigh)*e0*e0*expfac(ineigh)
     &                            *facback(ineigh)/e
         end do
         do ineigh = 1, nneigh
            de(1,ineigh) = de(1,ineigh) + delfrad(ineigh)
     &                        *(xneigh(ineigh)-xmet)/rneigh(ineigh)
            de(2,ineigh) = de(2,ineigh) + delfrad(ineigh)
     &                        *(yneigh(ineigh)-ymet)/rneigh(ineigh)
            de(3,ineigh) = de(3,ineigh) + delfrad(ineigh)
     &                        *(zneigh(ineigh)-zmet)/rneigh(ineigh)
         end do
         do ineigh = 1, nneigh
            dedrback(1,ineigh) = dfacback(1,ineigh)*e0
     &                              *expfac(ineigh)*delfdh(ineigh)
            dedrback(2,ineigh) = dfacback(2,ineigh)*e0
     &                              *expfac(ineigh)*delfdh(ineigh)
            dedrback(3,ineigh) = dfacback(3,ineigh)*e0
     &                              *expfac(ineigh)*delfdh(ineigh)
            detback(1,ineigh) = dftback(1,ineigh)*e0
     &                             *expfac(ineigh)*delfdh(ineigh)
            detback(2,ineigh) = dftback(2,ineigh)*e0
     &                             *expfac(ineigh)*delfdh(ineigh)
            detback(3,ineigh) = dftback(3,ineigh)*e0
     &                             *expfac(ineigh)*delfdh(ineigh)
         end do
   20    continue
         demet(1) = 0.0d0
         demet(2) = 0.0d0
         demet(3) = 0.0d0
         do ineigh = 1, nneigh
            demet(1) = demet(1) - de(1,ineigh) + detback(1,ineigh)
            demet(2) = demet(2) - de(2,ineigh) + detback(2,ineigh)
            demet(3) = demet(3) - de(3,ineigh) + detback(3,ineigh)
         end do
         elf = elf + e
         delf(1,i) = delf(1,i) + demet(1)
         delf(2,i) = delf(2,i) + demet(2)
         delf(3,i) = delf(3,i) + demet(3)
         do ineigh = 1, nneigh
            j = neighnum(ineigh)
            delf(1,j) = delf(1,j) + de(1,ineigh) - dedrback(1,ineigh)
     &                     - detback(1,ineigh)
            delf(2,j) = delf(2,j) + de(2,ineigh) - dedrback(2,ineigh)
     &                     - detback(2,ineigh)
            delf(3,j) = delf(3,j) + de(3,ineigh) - dedrback(3,ineigh)
     &                     - detback(3,ineigh)
         end do
         do ineigh = 1, nneigh
            j = neighnum(ineigh)
            if (n12(j) .le. 0)  call fatal
            do k0 = 1, n12(j)
               k = i12(k0,j)
               delf(1,k) = delf(1,k) + dedrback(1,ineigh)/dble(n12(j))
               delf(2,k) = delf(2,k) + dedrback(2,ineigh)/dble(n12(j))
               delf(3,k) = delf(3,k) + dedrback(3,ineigh)/dble(n12(j))
            end do
         end do
   30    continue
      end do
      return
      end
