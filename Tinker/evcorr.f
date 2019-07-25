c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2010 by Chuanjie Wu and Jay William Ponder  ##
c     ##                    All Rights Reserved                     ##
c     ################################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine evcorr  --  long range vdw energy correction  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "evcorr" computes the long range van der Waals correction
c     to the energy via numerical integration
c
c     literature reference:
c
c     M. P. Allen and D. J. Tildesley, "Computer Simulation of
c     Liquids", Oxford University Press, 1987, Section 2.8
c
c
      subroutine evcorr (elrc)
      use sizes
      use atoms
      use bound
      use boxes
      use limits
      use math
      use mutant
      use shunt
      use vdw
      use vdwpot
      implicit none
      integer i,j,k,it,kt
      integer nstep,ndelta,nvt
      integer, allocatable :: ivt(:)
      integer, allocatable :: jvt(:)
      integer, allocatable :: mvt(:)
      real*8 elrc,etot
      real*8 range,rdelta
      real*8 fi,fk,fim,fkm,fik
      real*8 e,eps,vlam1
      real*8 offset,taper
      real*8 rv,rv2,rv6,rv7
      real*8 r,r2,r3,r4
      real*8 r5,r6,r7
      real*8 p,p6,p12
      real*8 rho,tau,tau7
      real*8 expterm
      character*6 mode
c
c
c     zero out the long range van der Waals correction
c
      elrc = 0.0d0
c
c     only applicable if periodic boundaries are in use
c
      if (.not. use_bounds)  return
c
c     set the coefficients for the switching function
c
      mode = 'VDW'
      call switch (mode)
c
c     set number of steps and range for numerical integration
c
      nstep = 2
      range = 100.0d0
      ndelta = int(dble(nstep)*(range-cut))
      rdelta = (range-cut) / dble(ndelta)
      offset = cut - 0.5d0*rdelta
      vlam1 = 1.0d0 - vlambda
c
c     perform dynamic allocation of some local arrays
c
      allocate (ivt(n))
      allocate (jvt(n))
      allocate (mvt(n))
c
c     count the number of vdw types and their frequencies
c
      nvt = 0
      do i = 1, n
         it = jvdw(i)
         do k = 1, nvt
            if (ivt(k) .eq. it) then
               jvt(k) = jvt(k) + 1
               if (mut(i))  mvt(k) = mvt(k) + 1
               goto 10
            end if
         end do
         nvt = nvt + 1
         ivt(nvt) = it
         jvt(nvt) = 1
         mvt(nvt) = 0
         if (mut(i))  mvt(nvt) = 1
   10    continue
      end do
c
c     find the van der Waals energy via double loop search
c
      do i = 1, nvt
         it = ivt(i)
         fi = 4.0d0 * pi * dble(jvt(i))
         fim = 4.0d0 * pi * dble(mvt(i))
         do k = i, nvt
            kt = ivt(k)
            fk = dble(jvt(k))
            fkm = dble(mvt(k))
            fik = fi*fk - vlam1*(fim*(fk-fkm)+(fi-fim)*fkm)
            if (k .eq. i)  fik = 0.5d0 * fik
            rv = radmin(kt,it)
            eps = epsilon(kt,it)
            rv2 = rv * rv
            rv6 = rv2 * rv2 * rv2
            rv7 = rv6 * rv
            etot = 0.0d0
            do j = 1, ndelta
               r = offset + dble(j)*rdelta
               r2 = r * r
               r3 = r2 * r
               r6 = r3 * r3
               r7 = r6 * r
               e = 0.0d0
               if (vdwtyp .eq. 'LENNARD-JONES') then
                  p6 = rv6 / r6
                  p12 = p6 * p6
                  e = eps * (p12 - 2.0d0*p6)
               else if (vdwtyp .eq. 'BUFFERED-14-7') then
                  rho = r7 + ghal*rv7
                  tau = (dhal+1.0d0) / (r+dhal*rv)
                  tau7 = tau**7
                  e = eps * rv7 * tau7 * ((ghal+1.0d0)*rv7/rho-2.0d0)
               else if (vdwtyp.eq.'BUCKINGHAM' .or.
     &                  vdwtyp.eq.'MM3-HBOND') then
                  p = sqrt(rv2/r2)
                  p6 = rv6 / r6
                  expterm = abuck * exp(-bbuck/p)
                  e = eps * (expterm - cbuck*p6)
               end if
               if (r .lt. off) then
                  r4 = r2 * r2
                  r5 = r2 * r3
                  taper = c5*r5 + c4*r4 + c3*r3 + c2*r2 + c1*r + c0
                  e = e * (1.0d0-taper)
               end if
               etot = etot + e*rdelta*r2
            end do
            elrc = elrc + fik*etot
         end do
      end do
      elrc = elrc / volbox
c
c     perform deallocation of some local arrays
c
      deallocate (ivt)
      deallocate (jvt)
      deallocate (mvt)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine evcorr1  --  long range vdw energy & virial  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "evcorr1" computes the long range van der Waals correction
c     to the energy and virial via numerical integration
c
c     literature reference:
c
c     M. P. Allen and D. J. Tildesley, "Computer Simulation of
c     Liquids", Oxford University Press, 1987, Section 2.8
c
c
      subroutine evcorr1 (elrc,vlrc)
      use sizes
      use atoms
      use bound
      use boxes
      use limits
      use math
      use mutant
      use shunt
      use vdw
      use vdwpot
      implicit none
      integer i,j,k,it,kt
      integer nstep,ndelta,nvt
      integer, allocatable :: ivt(:)
      integer, allocatable :: jvt(:)
      integer, allocatable :: mvt(:)
      real*8 elrc,vlrc
      real*8 etot,vtot
      real*8 range,rdelta
      real*8 fi,fk,fim,fkm,fik
      real*8 e,de,eps
      real*8 offset,vlam1
      real*8 taper,dtaper
      real*8 rv,rv2,rv6,rv7
      real*8 r,r2,r3,r4
      real*8 r5,r6,r7
      real*8 p,p6,p12
      real*8 rho,tau,tau7
      real*8 dtau,gtau
      real*8 rvterm,expterm
      character*6 mode
c
c
c     zero out the long range van der Waals corrections
c
      elrc = 0.0d0
      vlrc = 0.0d0
c
c     only applicable if periodic boundaries are in use
c
      if (.not. use_bounds)  return
c
c     set the coefficients for the switching function
c
      mode = 'VDW'
      call switch (mode)
c
c     set number of steps and range for numerical integration
c
      nstep = 2
      range = 100.0d0
      ndelta = int(dble(nstep)*(range-cut))
      rdelta = (range-cut) / dble(ndelta)
      offset = cut - 0.5d0*rdelta
      vlam1 = 1.0d0 - vlambda
c
c     perform dynamic allocation of some local arrays
c
      allocate (ivt(n))
      allocate (jvt(n))
      allocate (mvt(n))
c
c     count the number of vdw types and their frequencies
c
      nvt = 0
      do i = 1, n
         it = jvdw(i)
         do k = 1, nvt
            if (ivt(k) .eq. it) then
               jvt(k) = jvt(k) + 1
               if (mut(i))  mvt(k) = mvt(k) + 1
               goto 10
            end if
         end do
         nvt = nvt + 1
         ivt(nvt) = it
         jvt(nvt) = 1
         mvt(nvt) = 0
         if (mut(i))  mvt(nvt) = 1
   10    continue
      end do
c
c     find the van der Waals energy via double loop search
c
      do i = 1, nvt
         it = ivt(i)
         fi = 4.0d0 * pi * dble(jvt(i))
         fim = 4.0d0 * pi * dble(mvt(i))
         do k = i, nvt
            kt = ivt(k)
            fk = dble(jvt(k))
            fkm = dble(mvt(k))
            fik = fi*fk - vlam1*(fim*(fk-fkm)+(fi-fim)*fkm)
            if (k .eq. i)  fik = 0.5d0 * fik
            rv = radmin(kt,it)
            eps = epsilon(kt,it)
            rv2 = rv * rv
            rv6 = rv2 * rv2 * rv2
            rv7 = rv6 * rv
            etot = 0.0d0
            vtot = 0.0d0
            do j = 1, ndelta
               r = offset + dble(j)*rdelta
               r2 = r * r
               r3 = r2 * r
               r6 = r3 * r3
               r7 = r6 * r
               e = 0.0d0
               de = 0.0d0
               if (vdwtyp .eq. 'LENNARD-JONES') then
                  p6 = rv6 / r6
                  p12 = p6 * p6
                  e = eps * (p12 - 2.0d0*p6)
                  de = eps * (p12-p6) * (-12.0d0/r)
               else if (vdwtyp .eq. 'BUFFERED-14-7') then
                  rho = r7 + ghal*rv7
                  tau = (dhal+1.0d0) / (r+dhal*rv)
                  tau7 = tau**7
                  dtau = tau / (dhal+1.0d0)
                  gtau = eps*tau7*r6*(ghal+1.0d0)*(rv7/rho)**2
                  e = eps * rv7 * tau7 * ((ghal+1.0d0)*rv7/rho-2.0d0)
                  de = -7.0d0 * (dtau*e+gtau)
               else if (vdwtyp.eq.'BUCKINGHAM' .or.
     &                  vdwtyp.eq.'MM3-HBOND') then
                  p = sqrt(rv2/r2)
                  p6 = rv6 / r6
                  rvterm = -bbuck / rv
                  expterm = abuck * exp(-bbuck/p)
                  e = eps * (expterm - cbuck*p6)
                  de = eps * (rvterm*expterm+6.0d0*cbuck*p6/r)
               end if
               if (r .lt. off) then
                  r4 = r2 * r2
                  r5 = r2 * r3
                  taper = c5*r5 + c4*r4 + c3*r3 + c2*r2 + c1*r + c0
                  dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                        + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                  de = de*(1.0d0-taper) - e*dtaper
                  e = e*(1.0d0-taper)
               end if
               etot = etot + e*rdelta*r2
               vtot = vtot + de*rdelta*r3
            end do
            elrc = elrc + fik*etot
            vlrc = vlrc + fik*vtot
         end do
      end do
      elrc = elrc / volbox
      vlrc = vlrc / (3.0d0*volbox)
c
c     perform deallocation of some local arrays
c
      deallocate (ivt)
      deallocate (jvt)
      deallocate (mvt)
      return
      end
