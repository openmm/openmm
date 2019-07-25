c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2001 by Anders Carlsson & Jay William Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine emetal  --  metal ligand field potential energy  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "emetal" calculates the transition metal ligand field energy
c
c     literature reference:
c
c     A. E. Carlsson and S. Zapata, "The Functional Form of Angular
c     Forces around Transition Metal Ions in Biomolecules", Biophysical
c     Journal, 81, 1-10 (2001)
c
c
      subroutine emetal
      use sizes
      use atomid
      use atoms
      use couple
      use energi
      use kchrge
      implicit none
      integer maxneigh
      parameter (maxneigh=10)
      integer i,j,k,jj,k0
      integer nneigh,ineigh,jneigh
      integer neighnum(maxneigh)
      real*8 e,e0,elf0
      real*8 dot,cij
      real*8 rcut,r2
      real*8 argument
      real*8 rback2,esqtet
      real*8 xmet,ymet,zmet
      real*8 alpha,beta
      real*8 kappa,r0cu
      real*8 rback(3,maxneigh)
      real*8 facback(maxneigh)
      real*8 xneigh(maxneigh)
      real*8 yneigh(maxneigh)
      real*8 zneigh(maxneigh)
      real*8 rneigh(maxneigh)
      real*8 expfac(maxneigh)
      real*8 angfac(maxneigh,maxneigh)
      real*8 cosmat(maxneigh,maxneigh)
c
c
c     zero out the metal ligand field energy
c
      elf = 0.0d0
c
c     begin ligand field calculation; for now, only for Cu+2
c
      do i = 1, n
         if (atomic(i).ne.29 .or. chg(type(i)).ne.2.0d0)  goto 20
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
   10       continue
         end do
c
c     calculate the energy for the current interaction
c
         do ineigh = 1, nneigh
            rneigh(ineigh) = (xneigh(ineigh)-xmet)**2
     &                          + (yneigh(ineigh)-ymet)**2
     &                          + (zneigh(ineigh)-zmet)**2
            rneigh(ineigh) = sqrt(rneigh(ineigh))
            expfac(ineigh) = exp(-kappa*(rneigh(ineigh)-r0cu))
            jj = neighnum(ineigh)
            if (atomic(jj) .eq. 8)
     &         expfac(ineigh) = 0.4d0 * expfac(ineigh)
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
            end do
         end do
         argument = 0.0d0
         do ineigh = 1, nneigh
            do jneigh = 1, nneigh
               argument = argument + expfac(ineigh)*expfac(jneigh)
     &                                  *angfac(ineigh,jneigh)
     &                              *facback(ineigh)*facback(jneigh)
            end do
         end do
c
c     increment the total metal ligand field energy
c
         e = 0.0d0
         if (argument .gt. 0) then
            e = -e0 * sqrt(argument)
            elf = elf + e
         end if
   20    continue
      end do
      return
      end
