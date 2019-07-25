c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine esolv1  --  solvation energy and derivatives  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "esolv1" calculates the implicit solvation energy and
c     first derivatives with respect to Cartesian coordinates
c     for surface area, generalized Born, generalized Kirkwood
c     and Poisson-Boltzmann solvation models
c
c
      subroutine esolv1
      use sizes
      use atoms
      use deriv
      use energi
      use limits
      use math
      use potent
      use solute
      use warp
      implicit none
      integer i
      real*8 e,ai,ri,rb
      real*8 term,probe
      real*8 esurf,ehp,eace
      real*8 ecav,edisp
      real*8, allocatable :: aes(:)
c
c
c     zero out the implicit solvation energy and derivatives
c
      es = 0.0d0
      esurf = 0.0d0
      ecav = 0.0d0
      edisp = 0.0d0
      ehp = 0.0d0
      eace = 0.0d0
      do i = 1, n
         drb(i) = 0.0d0
         des(1,i) = 0.0d0
         des(2,i) = 0.0d0
         des(3,i) = 0.0d0
      end do
      if (solvtyp(1:2) .eq. 'GK') then
         do i = 1, n
            drbp(i) = 0.0d0
         end do
      end if
c
c     set a value for the solvent molecule probe radius
c
      probe = 1.4d0
c
c     perform dynamic allocation of some local arrays
c
      allocate (aes(n))
c
c     solvation energy and derivs for surface area only models
c
      if (solvtyp.eq.'ASP' .or. solvtyp.eq.'SASA') then
         call surface1 (es,aes,des,rsolv,asolv,probe)
c
c     nonpolar energy and derivs for Onion method via exact area
c
      else if (solvtyp .eq. 'ONION') then
         call surface1 (esurf,aes,des,rsolv,asolv,probe)
         es = esurf
c
c     nonpolar energy and derivs as cavity plus dispersion
c
      else if (solvtyp.eq.'GK' .or. solvtyp.eq.'PB') then
         call enp1 (ecav,edisp)
         es = ecav + edisp
c
c     nonpolar energy and derivs as hydrophobic PMF term
c
      else if (solvtyp.eq.'GB-HPMF' .or. solvtyp.eq.'GK-HPMF'
     &            .or. solvtyp.eq.'PB-HPMF') then
         call ehpmf1 (ehp)
         es = ehp
c
c     nonpolar energy and derivs via ACE area approximation
c
      else
         term = 4.0d0 * pi
         do i = 1, n
            ai = asolv(i)
            ri = rsolv(i)
            rb = rborn(i)
            if (rb .ne. 0.0d0) then
               e = ai * term * (ri+probe)**2 * (ri/rb)**6
               eace = eace + e
               drb(i) = drb(i) - 6.0d0*e/rb
            end if
         end do
         es = eace
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (aes)
c
c     get polarization energy and derivatives for solvation methods
c
      if (solvtyp(1:2) .eq. 'GK') then
         if (.not.use_mpole .and. .not.use_polar) then
            call chkpole
            call rotpole
            call induce
         end if
         call egk1
      else if (solvtyp(1:2) .eq. 'PB') then
         call epb1
      else if (use_born) then
         if (use_smooth) then
            call egb1c
         else if (use_clist) then
            call egb1b
         else
            call egb1a
         end if
         call born1
      end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine egb1a  --  GB energy and derivs via double loop  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "egb1a" calculates the generalized Born electrostatic energy
c     and first derivatives of the GB/SA solvation models using a
c     double loop
c
c     note application of distance cutoff scaling directly to
c     the Born radii chain rule term "derb" is an approximation
c
c
      subroutine egb1a
      use sizes
      use atoms
      use charge
      use chgpot
      use deriv
      use energi
      use group
      use inter
      use molcul
      use shunt
      use solute
      use usage
      use virial
      implicit none
      integer i,k,ii,kk
      real*8 e,de,fgrp
      real*8 f,fi,fik
      real*8 fgb,fgb2,fgm
      real*8 rb2,rm2
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,r3,r4
      real*8 r5,r6,r7
      real*8 dwater,rbi,rbk
      real*8 dedx,dedy,dedz
      real*8 derb,drbi,drbk
      real*8 expterm,shift
      real*8 taper,dtaper
      real*8 trans,dtrans
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      logical proceed,usei
      character*6 mode
c
c
c     set the solvent dielectric and energy conversion factor
c
      if (nion .eq. 0)  return
      dwater = 78.3d0
      f = -electric * (1.0d0 - 1.0d0/dwater)
c
c     set cutoff distances and switching function coefficients
c
      mode = 'CHARGE'
      call switch (mode)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(nion,iion,use,x,y,z,f,
!$OMP& pchg,rborn,use_group,off,off2,cut,cut2,c0,c1,c2,c3,c4,c5,
!$OMP% f0,f1,f2,f3,f4,f5,f6,f7,molcule)
!$OMP& shared(es,einter,des,drb,vir)
!$OMP DO reduction(+:es,einter,des,drb,vir) schedule(guided)
c
c     calculate GB electrostatic polarization energy term
c
      do ii = 1, nion
         i = iion(ii)
         usei = use(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         fi = f * pchg(ii)
         rbi = rborn(i)
c
c     decide whether to compute the current interaction
c
         do kk = ii, nion
            k = iion(kk)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  rbk = rborn(k)
                  fik = fi * pchg(kk)
                  rb2 = rbi * rbk
                  expterm = exp(-0.25d0*r2/rb2)
                  fgb2 = r2 + rb2*expterm
                  fgb = sqrt(fgb2)
                  e = fik / fgb
                  de = -e * (r-0.25d0*r*expterm) / fgb2
                  derb = -e * expterm*(0.5d0+0.125d0*r2/rb2) / fgb2
c
c     use energy switching if near the cutoff distance
c
                  rm2 = (0.5d0 * (off+cut))**2
                  fgm = sqrt(rm2 + rb2*exp(-0.25d0*rm2/rb2))
                  shift = fik / fgm
                  e = e - shift
                  if (r2 .gt. cut2) then
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     r6 = r3 * r3
                     r7 = r3 * r4
                     taper = c5*r5 + c4*r4 + c3*r3
     &                          + c2*r2 + c1*r + c0
                     dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                           + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                     trans = fik * (f7*r7 + f6*r6 + f5*r5 + f4*r4
     &                               + f3*r3 + f2*r2 + f1*r + f0)
                     dtrans = fik * (7.0d0*f7*r6 + 6.0d0*f6*r5
     &                               + 5.0d0*f5*r4 + 4.0d0*f4*r3
     &                             + 3.0d0*f3*r2 + 2.0d0*f2*r + f1)
                     derb = derb * taper
                     de = e*dtaper + de*taper + dtrans
                     e = e*taper + trans
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     de = de * fgrp
                     derb = derb * fgrp
                  end if
c
c     increment the overall energy and derivative expressions
c
                  if (i .eq. k) then
                     e = 0.5d0 * e
                     es = es + e
                     drbi = derb * rbk
                     drb(i) = drb(i) + drbi
                  else
                     es = es + e
                     de = de / r
                     dedx = de * xr
                     dedy = de * yr
                     dedz = de * zr
                     des(1,i) = des(1,i) + dedx
                     des(2,i) = des(2,i) + dedy
                     des(3,i) = des(3,i) + dedz
                     des(1,k) = des(1,k) - dedx
                     des(2,k) = des(2,k) - dedy
                     des(3,k) = des(3,k) - dedz
                     drbi = derb * rbk
                     drbk = derb * rbi
                     drb(i) = drb(i) + drbi
                     drb(k) = drb(k) + drbk
c
c     increment the internal virial tensor components
c
                     vxx = xr * dedx
                     vyx = yr * dedx
                     vzx = zr * dedx
                     vyy = yr * dedy
                     vzy = zr * dedy
                     vzz = zr * dedz
                     vir(1,1) = vir(1,1) + vxx
                     vir(2,1) = vir(2,1) + vyx
                     vir(3,1) = vir(3,1) + vzx
                     vir(1,2) = vir(1,2) + vyx
                     vir(2,2) = vir(2,2) + vyy
                     vir(3,2) = vir(3,2) + vzy
                     vir(1,3) = vir(1,3) + vzx
                     vir(2,3) = vir(2,3) + vzy
                     vir(3,3) = vir(3,3) + vzz
                  end if
c
c     increment the total intermolecular energy
c
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + e
                  end if
               end if
            end if
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine egb1b  --  GB energy and derivs via pair list  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "egb1b" calculates the generalized Born electrostatic energy
c     and first derivatives of the GB/SA solvation models using a
c     neighbor list
c
c     note application of distance cutoff scaling directly to
c     the Born radii chain rule term "derb" is an approximation
c
c
      subroutine egb1b
      use sizes
      use atoms
      use charge
      use chgpot
      use deriv
      use energi
      use group
      use inter
      use molcul
      use neigh
      use shunt
      use solute
      use usage
      use virial
      implicit none
      integer i,k
      integer ii,kk,kkk
      real*8 e,de,fgrp
      real*8 f,fi,fik
      real*8 fgb,fgb2,fgm
      real*8 rb2,rm2
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,r3,r4
      real*8 r5,r6,r7
      real*8 dwater,rbi,rbk
      real*8 dedx,dedy,dedz
      real*8 derb,drbi,drbk
      real*8 expterm,shift
      real*8 taper,dtaper
      real*8 trans,dtrans
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      logical proceed,usei
      character*6 mode
c
c
c     set the solvent dielectric and energy conversion factor
c
      if (nion .eq. 0)  return
      dwater = 78.3d0
      f = -electric * (1.0d0 - 1.0d0/dwater)
c
c     set cutoff distances and switching function coefficients
c
      mode = 'CHARGE'
      call switch (mode)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(nion,iion,use,x,y,z,
!$OMP& f,pchg,rborn,nelst,elst,use_group,off,off2,cut,cut2,
!$OMP% c0,c1,c2,c3,c4,c5,f0,f1,f2,f3,f4,f5,f6,f7,molcule)
!$OMP& shared(es,einter,des,drb,vir)
!$OMP DO reduction(+:es,einter,des,drb,vir) schedule(guided)
c
c     calculate GB electrostatic polarization energy term
c
      do ii = 1, nion
         i = iion(ii)
         usei = use(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         fi = f * pchg(ii)
         rbi = rborn(i)
c
c     calculate the self-energy term for the current atom
c
         fik = fi * pchg(ii)
         rb2 = rbi * rbi
         e = fik / rbi
         derb = -0.5d0 * e / rb2
         rm2 = (0.5d0 * (off+cut))**2
         fgm = sqrt(rm2 + rb2*exp(-0.25d0*rm2/rb2))
         shift = fik / fgm
         e = e - shift
         es = es + 0.5d0*e
         drbi = derb * rbi
         drb(i) = drb(i) + drbi
c
c     decide whether to compute the current interaction
c
         do kkk = 1, nelst(ii)
            kk = elst(kkk,ii)
            k = iion(kk)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  rbk = rborn(k)
                  fik = fi * pchg(kk)
                  rb2 = rbi * rbk
                  expterm = exp(-0.25d0*r2/rb2)
                  fgb2 = r2 + rb2*expterm
                  fgb = sqrt(fgb2)
                  e = fik / fgb
                  de = -e * (r-0.25d0*r*expterm) / fgb2
                  derb = -e * expterm*(0.5d0+0.125d0*r2/rb2) / fgb2
c
c     use energy switching if near the cutoff distance
c
                  rm2 = (0.5d0 * (off+cut))**2
                  fgm = sqrt(rm2 + rb2*exp(-0.25d0*rm2/rb2))
                  shift = fik / fgm
                  e = e - shift
                  if (r2 .gt. cut2) then
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     r6 = r3 * r3
                     r7 = r3 * r4
                     taper = c5*r5 + c4*r4 + c3*r3
     &                          + c2*r2 + c1*r + c0
                     dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                           + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                     trans = fik * (f7*r7 + f6*r6 + f5*r5 + f4*r4
     &                               + f3*r3 + f2*r2 + f1*r + f0)
                     dtrans = fik * (7.0d0*f7*r6 + 6.0d0*f6*r5
     &                               + 5.0d0*f5*r4 + 4.0d0*f4*r3
     &                             + 3.0d0*f3*r2 + 2.0d0*f2*r + f1)
                     derb = derb * taper
                     de = e*dtaper + de*taper + dtrans
                     e = e*taper + trans
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     de = de * fgrp
                     derb = derb * fgrp
                  end if
c
c     increment the overall energy and derivative expressions
c
                  es = es + e
                  de = de / r
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
                  des(1,i) = des(1,i) + dedx
                  des(2,i) = des(2,i) + dedy
                  des(3,i) = des(3,i) + dedz
                  des(1,k) = des(1,k) - dedx
                  des(2,k) = des(2,k) - dedy
                  des(3,k) = des(3,k) - dedz
                  drbi = derb * rbk
                  drbk = derb * rbi
                  drb(i) = drb(i) + drbi
                  drb(k) = drb(k) + drbk
c
c     increment the internal virial tensor components
c
                  vxx = xr * dedx
                  vyx = yr * dedx
                  vzx = zr * dedx
                  vyy = yr * dedy
                  vzy = zr * dedy
                  vzz = zr * dedz
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
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + e
                  end if
               end if
            end if
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine egb1c  --  GB energy and derivs for smoothing  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "egb1c" calculates the generalized Born energy and first
c     derivatives of the GB/SA solvation models for use with
c     potential smoothing methods
c
c
      subroutine egb1c
      use sizes
      use atoms
      use charge
      use chgpot
      use deriv
      use energi
      use group
      use inter
      use math
      use molcul
      use solute
      use usage
      use virial
      use warp
      implicit none
      integer i,k,ii,kk
      real*8 e,de,fgrp
      real*8 f,fi,fik
      real*8 fgb,fgb2
      real*8 rb2,width
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,sterm
      real*8 expterm
      real*8 dwater,rbi,rbk
      real*8 dedx,dedy,dedz
      real*8 derb,drbi,drbk
      real*8 erf,erfterm,term
      real*8 wterm,rterm,bterm
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      logical proceed,usei
      external erf
c
c
c     set the solvent dielectric and energy conversion factor
c
      if (nion .eq. 0)  return
      dwater = 78.3d0
      f = -electric * (1.0d0 - 1.0d0/dwater)
c
c     set the extent of smoothing to be performed
c
      sterm = 0.5d0 / sqrt(diffc)
c
c     calculate GB electrostatic polarization energy term
c
      do ii = 1, nion
         i = iion(ii)
         usei = use(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         fi = f * pchg(ii)
         rbi = rborn(i)
c
c     decide whether to compute the current interaction
c
         do kk = ii, nion
            k = iion(kk)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
               r2 = xr*xr + yr*yr + zr*zr
               r = sqrt(r2)
               rbk = rborn(k)
               fik = fi * pchg(kk)
               rb2 = rbi * rbk
               expterm = exp(-0.25d0*r2/rb2)
               fgb2 = r2 + rb2*expterm
               fgb = sqrt(fgb2)
               e = fik / fgb
               de = -e * (r-0.25d0*r*expterm) / fgb2
               derb = -e * expterm*(0.5d0+0.125d0*r2/rb2) / fgb2
c
c     use a smoothable GB analogous to the Coulomb solution
c
               if (deform .gt. 0.0d0) then
                  wterm = exp(-0.006d0*rb2/deform)
                  width = sterm / sqrt(deform+0.15d0*rb2*wterm)
                  erfterm = erf(width*fgb)
                  term = width * exp(-(width*fgb)**2) / sqrtpi
                  rterm = term * (2.0d0*r-0.5d0*r*expterm)/fgb
                  bterm = term * ((expterm*(1.0d0+0.25d0*r2/rb2)/fgb)
     &                              - (fgb*(width/sterm)**2) * wterm
     &                                 * (0.15d0-0.0009d0*rb2/deform))
                  derb = derb*erfterm + e*bterm
                  de = de*erfterm + e*rterm
                  e = e * erfterm
               end if
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  e = e * fgrp
                  de = de * fgrp
                  derb = derb * fgrp
               end if
c
c     increment the overall energy and derivative expressions
c
               if (i .eq. k) then
                  e = 0.5d0 * e
                  es = es + e
                  drbi = derb * rbk
                  drb(i) = drb(i) + drbi
               else
                  es = es + e
                  de = de / r
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
                  des(1,i) = des(1,i) + dedx
                  des(2,i) = des(2,i) + dedy
                  des(3,i) = des(3,i) + dedz
                  des(1,k) = des(1,k) - dedx
                  des(2,k) = des(2,k) - dedy
                  des(3,k) = des(3,k) - dedz
                  drbi = derb * rbk
                  drbk = derb * rbi
                  drb(i) = drb(i) + drbi
                  drb(k) = drb(k) + drbk
c
c     increment the internal virial tensor components
c
                  vxx = xr * dedx
                  vyx = yr * dedx
                  vzx = zr * dedx
                  vyy = yr * dedy
                  vzy = zr * dedy
                  vzz = zr * dedz
                  vir(1,1) = vir(1,1) + vxx
                  vir(2,1) = vir(2,1) + vyx
                  vir(3,1) = vir(3,1) + vzx
                  vir(1,2) = vir(1,2) + vyx
                  vir(2,2) = vir(2,2) + vyy
                  vir(3,2) = vir(3,2) + vzy
                  vir(1,3) = vir(1,3) + vzx
                  vir(2,3) = vir(2,3) + vzy
                  vir(3,3) = vir(3,3) + vzz
               end if
c
c     increment the total intermolecular energy
c
               if (molcule(i) .ne. molcule(k)) then
                  einter = einter + e
               end if
            end if
         end do
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine egk1  --  generalized Kirkwood energy & derivs  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "egk1" calculates the implicit solvation energy and derivatives
c     via the generalized Kirkwood plus nonpolar implicit solvation
c
c
      subroutine egk1
      use sizes
      use energi
      use limits
      use potent
      implicit none
c
c
c     setup the multipoles for solvation only calculations
c
      if (.not.use_mpole .and. .not.use_polar) then
         call chkpole
         call rotpole
         call induce
      end if
c
c     compute the generalized Kirkwood energy and gradient
c
      call egk1a
      call born1
c
c     correct energy and derivatives for vacuum to polarized state
c
      if (use_mlist) then
         call ediff1b
      else
         call ediff1a
      end if
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine egk1a  --  find GK energy and derivatives  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "egk1a" calculates the electrostatic portion of the implicit
c     solvation energy and derivatives via the generalized Kirkwood
c     model
c
c
      subroutine egk1a
      use sizes
      use atoms
      use charge
      use chgpot
      use deriv
      use energi
      use gkstuf
      use group
      use inter
      use molcul
      use mpole
      use polar
      use polpot
      use shunt
      use solute
      use usage
      use virial
      implicit none
      integer i,j,k,ii,kk
      real*8 e,ei,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xr2,yr2,zr2
      real*8 ci,ck
      real*8 uxi,uyi,uzi
      real*8 uxk,uyk,uzk
      real*8 qxxi,qxyi,qxzi
      real*8 qyyi,qyzi,qzzi
      real*8 qxxk,qxyk,qxzk
      real*8 qyyk,qyzk,qzzk
      real*8 dxi,dyi,dzi
      real*8 dxk,dyk,dzk
      real*8 pxi,pyi,pzi
      real*8 pxk,pyk,pzk
      real*8 sxi,syi,szi
      real*8 sxk,syk,szk
      real*8 r2,rb2
      real*8 dedx,dedy,dedz
      real*8 drbi,drbk
      real*8 dpdx,dpdy,dpdz
      real*8 dpbi,dpbk
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 dwater
      real*8 fc,fd,fq
      real*8 rbi,rbk
      real*8 expterm
      real*8 gf,gf2,gf3,gf5
      real*8 gf7,gf9,gf11
      real*8 expc,dexpc
      real*8 expc1,expcdexpc
      real*8 expcr,dexpcr
      real*8 dgfdr
      real*8 esym,ewi,ewk
      real*8 desymdx,dewidx,dewkdx
      real*8 desymdy,dewidy,dewkdy
      real*8 desymdz,dewidz,dewkdz
      real*8 dsumdr,desymdr
      real*8 dewidr,dewkdr
      real*8 dsymdr
      real*8 esymi,ewii,ewki
      real*8 dpsymdx,dpwidx,dpwkdx
      real*8 dpsymdy,dpwidy,dpwkdy
      real*8 dpsymdz,dpwidz,dpwkdz
      real*8 dwipdr,dwkpdr,duvdr
      real*8 a(0:5,0:3)
      real*8 b(0:4,0:2)
      real*8 fid(3),fkd(3)
      real*8 fidg(3,3),fkdg(3,3)
      real*8 gc(30),gux(30)
      real*8 guy(30),guz(30)
      real*8 gqxx(30),gqxy(30)
      real*8 gqxz(30),gqyy(30)
      real*8 gqyz(30),gqzz(30)
      real*8, allocatable :: trq(:,:)
      real*8, allocatable :: trqi(:,:)
      logical proceed,usei
      character*6 mode
c
c
c     set the bulk dielectric constant to the water value
c
      if (npole .eq. 0)  return
      dwater = 78.3d0
      fc = electric * 1.0d0 * (1.0d0-dwater)/(0.0d0+1.0d0*dwater)
      fd = electric * 2.0d0 * (1.0d0-dwater)/(1.0d0+2.0d0*dwater)
      fq = electric * 3.0d0 * (1.0d0-dwater)/(2.0d0+3.0d0*dwater)
c
c     set cutoff distances and switching function coefficients
c
      mode = 'MPOLE'
      call switch (mode)
c
c     perform dynamic allocation of some local arrays
c
      allocate (trq(3,n))
      allocate (trqi(3,n))
c
c     initialize local variables for OpenMP calculation
c
      do i = 1, n
         do j = 1, 3
            trq(j,i) = 0.0d0
            trqi(j,i) = 0.0d0
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(npole,ipole,use,x,y,z,rborn,
!$OMP& rpole,uinds,uinps,use_group,off2,gkc,fc,fd,fq,poltyp,molcule)
!$OMP& shared(es,einter,des,drb,drbp,trq,trqi,vir)
!$OMP DO reduction(+:es,einter,des,drb,drbp,trq,trqi,vir)
!$OMP& schedule(guided)
c
c     calculate GK electrostatic solvation free energy
c
      do ii = 1, npole
         i = ipole(ii)
         usei = use(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         rbi = rborn(i)
         ci = rpole(1,ii)
         uxi = rpole(2,ii)
         uyi = rpole(3,ii)
         uzi = rpole(4,ii)
         qxxi = rpole(5,ii)
         qxyi = rpole(6,ii)
         qxzi = rpole(7,ii)
         qyyi = rpole(9,ii)
         qyzi = rpole(10,ii)
         qzzi = rpole(13,ii)
         dxi = uinds(1,ii)
         dyi = uinds(2,ii)
         dzi = uinds(3,ii)
         pxi = uinps(1,ii)
         pyi = uinps(2,ii)
         pzi = uinps(3,ii)
         sxi = dxi + pxi
         syi = dyi + pyi
         szi = dzi + pzi
c
c     decide whether to compute the current interaction
c
         do kk = ii, npole
            k = ipole(kk)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = x(k) - xi
               yr = y(k) - yi
               zr = z(k) - zi
               xr2 = xr*xr
               yr2 = yr*yr
               zr2 = zr*zr
               r2 = xr2 + yr2 + zr2
               if (r2 .le. off2) then
                  rbk = rborn(k)
                  ck = rpole(1,kk)
                  uxk = rpole(2,kk)
                  uyk = rpole(3,kk)
                  uzk = rpole(4,kk)
                  qxxk = rpole(5,kk)
                  qxyk = rpole(6,kk)
                  qxzk = rpole(7,kk)
                  qyyk = rpole(9,kk)
                  qyzk = rpole(10,kk)
                  qzzk = rpole(13,kk)
                  dxk = uinds(1,kk)
                  dyk = uinds(2,kk)
                  dzk = uinds(3,kk)
                  pxk = uinps(1,kk)
                  pyk = uinps(2,kk)
                  pzk = uinps(3,kk)
                  sxk = dxk + pxk
                  syk = dyk + pyk
                  szk = dzk + pzk
                  rb2 = rbi * rbk
                  expterm = exp(-r2/(gkc*rb2))
                  expc = expterm / gkc
                  expcr = r2*expterm / (gkc*gkc*rb2*rb2)
                  dexpc = -2.0d0 / (gkc*rb2)
                  dexpcr = 2.0d0 / (gkc*rb2*rb2)
                  dgfdr = 0.5d0 * expterm * (1.0d0+r2/(rb2*gkc))
                  gf2 = 1.0d0 / (r2+rb2*expterm)
                  gf = sqrt(gf2)
                  gf3 = gf2 * gf
                  gf5 = gf3 * gf2
                  gf7 = gf5 * gf2
                  gf9 = gf7 * gf2
                  gf11 = gf9 * gf2
c
c     reaction potential auxiliary terms
c
                  a(0,0) = gf
                  a(1,0) = -gf3
                  a(2,0) = 3.0d0 * gf5
                  a(3,0) = -15.0d0 * gf7
                  a(4,0) = 105.0d0 * gf9
                  a(5,0) = -945.0d0 * gf11
c
c     Born radii derivatives of reaction potential auxiliary terms
c
                  b(0,0) = dgfdr * a(1,0)
                  b(1,0) = dgfdr * a(2,0)
                  b(2,0) = dgfdr * a(3,0)
                  b(3,0) = dgfdr * a(4,0)
                  b(4,0) = dgfdr * a(5,0)
c
c     get reaction potential gradient auxiliary terms
c
                  expc1 = 1.0d0 - expc
                  a(0,1) = expc1 * a(1,0)
                  a(1,1) = expc1 * a(2,0)
                  a(2,1) = expc1 * a(3,0)
                  a(3,1) = expc1 * a(4,0)
                  a(4,1) = expc1 * a(5,0)
c
c     Born radii derivs of reaction potential gradient auxiliary terms
c
                  b(0,1) = b(1,0) - expcr*a(1,0) - expc*b(1,0)
                  b(1,1) = b(2,0) - expcr*a(2,0) - expc*b(2,0)
                  b(2,1) = b(3,0) - expcr*a(3,0) - expc*b(3,0)
                  b(3,1) = b(4,0) - expcr*a(4,0) - expc*b(4,0)
c
c     get 2nd reaction potential gradient auxiliary terms
c
                  expcdexpc = -expc * dexpc
                  a(0,2) = expc1*a(1,1) + expcdexpc*a(1,0)
                  a(1,2) = expc1*a(2,1) + expcdexpc*a(2,0)
                  a(2,2) = expc1*a(3,1) + expcdexpc*a(3,0)
                  a(3,2) = expc1*a(4,1) + expcdexpc*a(4,0)
c
c     Born radii derivatives of the 2nd reaction potential
c     gradient auxiliary terms
c
                  b(0,2) = b(1,1) - (expcr*(a(1,1) + dexpc*a(1,0))
     &                     + expc*(b(1,1)+dexpcr*a(1,0)+dexpc*b(1,0)))
                  b(1,2) = b(2,1) - (expcr*(a(2,1) + dexpc*a(2,0))
     &                     + expc*(b(2,1)+dexpcr*a(2,0)+dexpc*b(2,0)))
                  b(2,2) = b(3,1) - (expcr*(a(3,1) + dexpc*a(3,0))
     &                     + expc*(b(3,1)+dexpcr*a(3,0)+dexpc*b(3,0)))
c
c     get 3rd reaction potential gradient auxiliary terms
c
                  expcdexpc = 2.0d0 * expcdexpc
                  a(0,3) = expc1*a(1,2) + expcdexpc*a(1,1)
                  a(1,3) = expc1*a(2,2) + expcdexpc*a(2,1)
                  a(2,3) = expc1*a(3,2) + expcdexpc*a(3,1)
                  expcdexpc = -expc * dexpc * dexpc
                  a(0,3) = a(0,3) + expcdexpc*a(1,0)
                  a(1,3) = a(1,3) + expcdexpc*a(2,0)
                  a(2,3) = a(2,3) + expcdexpc*a(3,0)
c
c     multiply the auxillary terms by their dielectric functions
c
                  a(0,0) = fc * a(0,0)
                  a(0,1) = fc * a(0,1)
                  a(0,2) = fc * a(0,2)
                  a(0,3) = fc * a(0,3)
                  b(0,0) = fc * b(0,0)
                  b(0,1) = fc * b(0,1)
                  b(0,2) = fc * b(0,2)
                  a(1,0) = fd * a(1,0)
                  a(1,1) = fd * a(1,1)
                  a(1,2) = fd * a(1,2)
                  a(1,3) = fd * a(1,3)
                  b(1,0) = fd * b(1,0)
                  b(1,1) = fd * b(1,1)
                  b(1,2) = fd * b(1,2)
                  a(2,0) = fq * a(2,0)
                  a(2,1) = fq * a(2,1)
                  a(2,2) = fq * a(2,2)
                  a(2,3) = fq * a(2,3)
                  b(2,0) = fq * b(2,0)
                  b(2,1) = fq * b(2,1)
                  b(2,2) = fq * b(2,2)
c
c     unweighted reaction potential tensor
c
                  gc(1) = a(0,0)
                  gux(1) = xr * a(1,0)
                  guy(1) = yr * a(1,0)
                  guz(1) = zr * a(1,0)
                  gqxx(1) = xr2 * a(2,0)
                  gqyy(1) = yr2 * a(2,0)
                  gqzz(1) = zr2 * a(2,0)
                  gqxy(1) = xr * yr * a(2,0)
                  gqxz(1) = xr * zr * a(2,0)
                  gqyz(1) = yr * zr * a(2,0)
c
c     Born radii derivs of unweighted reaction potential tensor
c
                  gc(21) = b(0,0)
                  gux(21) = xr * b(1,0)
                  guy(21) = yr * b(1,0)
                  guz(21) = zr * b(1,0)
                  gqxx(21) = xr2 * b(2,0)
                  gqyy(21) = yr2 * b(2,0)
                  gqzz(21) = zr2 * b(2,0)
                  gqxy(21) = xr * yr * b(2,0)
                  gqxz(21) = xr * zr * b(2,0)
                  gqyz(21) = yr * zr * b(2,0)
c
c     unweighted reaction potential gradient tensor
c
                  gc(2) = xr * a(0,1)
                  gc(3) = yr * a(0,1)
                  gc(4) = zr * a(0,1)
                  gux(2) = a(1,0) + xr2*a(1,1)
                  gux(3) = xr * yr * a(1,1)
                  gux(4) = xr * zr * a(1,1)
                  guy(2) = gux(3)
                  guy(3) = a(1,0) + yr2*a(1,1)
                  guy(4) = yr * zr * a(1,1)
                  guz(2) = gux(4)
                  guz(3) = guy(4)
                  guz(4) = a(1,0) + zr2*a(1,1)
                  gqxx(2) = xr * (2.0d0*a(2,0)+xr2*a(2,1))
                  gqxx(3) = yr * xr2 * a(2,1)
                  gqxx(4) = zr * xr2 * a(2,1)
                  gqyy(2) = xr * yr2 * a(2,1)
                  gqyy(3) = yr * (2.0d0*a(2,0)+yr2*a(2,1))
                  gqyy(4) = zr * yr2 * a(2,1)
                  gqzz(2) = xr * zr2 * a(2,1)
                  gqzz(3) = yr * zr2 * a(2,1)
                  gqzz(4) = zr * (2.0d0*a(2,0)+zr2*a(2,1))
                  gqxy(2) = yr * (a(2,0)+xr2*a(2,1))
                  gqxy(3) = xr * (a(2,0)+yr2*a(2,1))
                  gqxy(4) = zr * xr * yr * a(2,1)
                  gqxz(2) = zr * (a(2,0)+xr2*a(2,1))
                  gqxz(3) = gqxy(4)
                  gqxz(4) = xr * (a(2,0)+zr2*a(2,1))
                  gqyz(2) = gqxy(4)
                  gqyz(3) = zr * (a(2,0)+yr2*a(2,1))
                  gqyz(4) = yr * (a(2,0)+zr2*a(2,1))
c
c     Born derivs of the unweighted reaction potential gradient tensor
c
                  gc(22) = xr * b(0,1)
                  gc(23) = yr * b(0,1)
                  gc(24) = zr * b(0,1)
                  gux(22) = b(1,0) + xr2*b(1,1)
                  gux(23) = xr * yr * b(1,1)
                  gux(24) = xr * zr * b(1,1)
                  guy(22) = gux(23)
                  guy(23) = b(1,0) + yr2*b(1,1)
                  guy(24) = yr * zr * b(1,1)
                  guz(22) = gux(24)
                  guz(23) = guy(24)
                  guz(24) = b(1,0) + zr2*b(1,1)
                  gqxx(22) = xr * (2.0d0*b(2,0)+xr2*b(2,1))
                  gqxx(23) = yr * xr2 * b(2,1)
                  gqxx(24) = zr * xr2 * b(2,1)
                  gqyy(22) = xr * yr2 * b(2,1)
                  gqyy(23) = yr * (2.0d0*b(2,0)+yr2*b(2,1))
                  gqyy(24) = zr * yr2 * b(2,1)
                  gqzz(22) = xr * zr2 * b(2,1)
                  gqzz(23) = yr * zr2 * b(2,1)
                  gqzz(24) = zr * (2.0d0*b(2,0) + zr2*b(2,1))
                  gqxy(22) = yr * (b(2,0)+xr2*b(2,1))
                  gqxy(23) = xr * (b(2,0)+yr2*b(2,1))
                  gqxy(24) = zr * xr * yr * b(2,1)
                  gqxz(22) = zr * (b(2,0)+xr2*b(2,1))
                  gqxz(23) = gqxy(24)
                  gqxz(24) = xr * (b(2,0)+zr2*b(2,1))
                  gqyz(22) = gqxy(24)
                  gqyz(23) = zr * (b(2,0)+yr2*b(2,1))
                  gqyz(24) = yr * (b(2,0)+zr2*b(2,1))
c
c     unweighted 2nd reaction potential gradient tensor
c
                  gc(5) = a(0,1) + xr2*a(0,2)
                  gc(6) = xr * yr * a(0,2)
                  gc(7) = xr * zr * a(0,2)
                  gc(8) = a(0,1) + yr2*a(0,2)
                  gc(9) = yr * zr * a(0,2)
                  gc(10) = a(0,1) + zr2*a(0,2)
                  gux(5) = xr * (3.0d0*a(1,1)+xr2*a(1,2))
                  gux(6) = yr * (a(1,1)+xr2*a(1,2))
                  gux(7) = zr * (a(1,1)+xr2*a(1,2))
                  gux(8) = xr * (a(1,1)+yr2*a(1,2))
                  gux(9) = zr * xr * yr * a(1,2)
                  gux(10) = xr * (a(1,1)+zr2*a(1,2))
                  guy(5) = yr * (a(1,1)+xr2*a(1,2))
                  guy(6) = xr * (a(1,1)+yr2*a(1,2))
                  guy(7) = gux(9)
                  guy(8) = yr * (3.0d0*a(1,1)+yr2*a(1,2))
                  guy(9) = zr * (a(1,1)+yr2*a(1,2))
                  guy(10) = yr * (a(1,1)+zr2*a(1,2))
                  guz(5) = zr * (a(1,1)+xr2*a(1,2))
                  guz(6) = gux(9)
                  guz(7) = xr * (a(1,1)+zr2*a(1,2))
                  guz(8) = zr * (a(1,1)+yr2*a(1,2))
                  guz(9) = yr * (a(1,1)+zr2*a(1,2))
                  guz(10) = zr * (3.0d0*a(1,1)+zr2*a(1,2))
                  gqxx(5) = 2.0d0*a(2,0) + xr2*(5.0d0*a(2,1)+xr2*a(2,2))
                  gqxx(6) = yr * xr * (2.0d0*a(2,1)+xr2*a(2,2))
                  gqxx(7) = zr * xr * (2.0d0*a(2,1)+xr2*a(2,2))
                  gqxx(8) = xr2 * (a(2,1)+yr2*a(2,2))
                  gqxx(9) = zr * yr * xr2 * a(2,2)
                  gqxx(10) = xr2 * (a(2,1)+zr2*a(2,2))
                  gqyy(5) = yr2 * (a(2,1)+xr2*a(2,2))
                  gqyy(6) = xr * yr * (2.0d0*a(2,1)+yr2*a(2,2))
                  gqyy(7) = xr * zr * yr2 * a(2,2)
                  gqyy(8) = 2.0d0*a(2,0) + yr2*(5.0d0*a(2,1)+yr2*a(2,2))
                  gqyy(9) = yr * zr * (2.0d0*a(2,1)+yr2*a(2,2))
                  gqyy(10) = yr2 * (a(2,1)+zr2*a(2,2))
                  gqzz(5) = zr2 * (a(2,1)+xr2*a(2,2))
                  gqzz(6) = xr * yr * zr2 * a(2,2)
                  gqzz(7) = xr * zr * (2.0d0*a(2,1)+zr2*a(2,2))
                  gqzz(8) = zr2 * (a(2,1)+yr2*a(2,2))
                  gqzz(9) = yr * zr * (2.0d0*a(2,1)+zr2*a(2,2))
                  gqzz(10) = 2.0d0*a(2,0)
     &                          + zr2*(5.0d0*a(2,1)+zr2*a(2,2))
                  gqxy(5) = xr * yr * (3.0d0*a(2,1)+xr2*a(2,2))
                  gqxy(6) = a(2,0) + (xr2+yr2)*a(2,1) + xr2*yr2*a(2,2)
                  gqxy(7) = zr * yr * (a(2,1)+xr2*a(2,2))
                  gqxy(8) = xr * yr * (3.0d0*a(2,1)+yr2*a(2,2))
                  gqxy(9) = zr * xr * (a(2,1)+yr2*a(2,2))
                  gqxy(10) = xr * yr * (a(2,1)+zr2*a(2,2))
                  gqxz(5) = xr * zr * (3.0d0*a(2,1)+xr2*a(2,2))
                  gqxz(6) = yr * zr * (a(2,1)+xr2*a(2,2))
                  gqxz(7) = a(2,0) + (xr2+zr2)*a(2,1) + xr2*zr2*a(2,2)
                  gqxz(8) = xr * zr * (a(2,1)+yr2*a(2,2))
                  gqxz(9) = xr * yr * (a(2,1)+zr2*a(2,2))
                  gqxz(10) = xr * zr * (3.0d0*a(2,1)+zr2*a(2,2))
                  gqyz(5) = zr * yr * (a(2,1)+xr2*a(2,2))
                  gqyz(6) = xr * zr * (a(2,1)+yr2*a(2,2))
                  gqyz(7) = xr * yr * (a(2,1)+zr2*a(2,2))
                  gqyz(8) = yr * zr * (3.0d0*a(2,1)+yr2*a(2,2))
                  gqyz(9) = a(2,0) + (yr2+zr2)*a(2,1) + yr2*zr2*a(2,2)
                  gqyz(10) = yr * zr * (3.0d0*a(2,1)+zr2*a(2,2))
c
c     Born radii derivatives of the unweighted 2nd reaction
c     potential gradient tensor
c
                  gc(25) = b(0,1) + xr2*b(0,2)
                  gc(26) = xr * yr * b(0,2)
                  gc(27) = xr * zr * b(0,2)
                  gc(28) = b(0,1) + yr2*b(0,2)
                  gc(29) = yr * zr * b(0,2)
                  gc(30) = b(0,1) + zr2*b(0,2)
                  gux(25) = xr * (3.0d0*b(1,1)+xr2*b(1,2))
                  gux(26) = yr * (b(1,1)+xr2*b(1,2))
                  gux(27) = zr * (b(1,1)+xr2*b(1,2))
                  gux(28) = xr * (b(1,1)+yr2*b(1,2))
                  gux(29) = zr * xr * yr * b(1,2)
                  gux(30) = xr * (b(1,1)+zr2*b(1,2))
                  guy(25) = yr * (b(1,1)+xr2*b(1,2))
                  guy(26) = xr * (b(1,1)+yr2*b(1,2))
                  guy(27) = gux(29)
                  guy(28) = yr * (3.0d0*b(1,1)+yr2*b(1,2))
                  guy(29) = zr * (b(1,1)+yr2*b(1,2))
                  guy(30) = yr * (b(1,1)+zr2*b(1,2))
                  guz(25) = zr * (b(1,1)+xr2*b(1,2))
                  guz(26) = gux(29)
                  guz(27) = xr * (b(1,1)+zr2*b(1,2))
                  guz(28) = zr * (b(1,1)+yr2*b(1,2))
                  guz(29) = yr * (b(1,1)+zr2*b(1,2))
                  guz(30) = zr * (3.0d0*b(1,1)+zr2*b(1,2))
                  gqxx(25) = 2.0d0*b(2,0)
     &                          + xr2*(5.0d0*b(2,1)+xr2*b(2,2))
                  gqxx(26) = yr * xr * (2.0d0*b(2,1)+xr2*b(2,2))
                  gqxx(27) = zr * xr * (2.0d0*b(2,1)+xr2*b(2,2))
                  gqxx(28) = xr2 * (b(2,1)+yr2*b(2,2))
                  gqxx(29) = zr * yr * xr2 * b(2,2)
                  gqxx(30) = xr2 * (b(2,1)+zr2*b(2,2))
                  gqyy(25) = yr2 * (b(2,1)+xr2*b(2,2))
                  gqyy(26) = xr * yr * (2.0d0*b(2,1)+yr2*b(2,2))
                  gqyy(27) = xr * zr * yr2 * b(2,2)
                  gqyy(28) = 2.0d0*b(2,0)
     &                          + yr2*(5.0d0*b(2,1)+yr2*b(2,2))
                  gqyy(29) = yr * zr * (2.0d0*b(2,1)+yr2*b(2,2))
                  gqyy(30) = yr2 * (b(2,1)+zr2*b(2,2))
                  gqzz(25) = zr2 * (b(2,1)+xr2*b(2,2))
                  gqzz(26) = xr * yr * zr2 * b(2,2)
                  gqzz(27) = xr * zr * (2.0d0*b(2,1)+zr2*b(2,2))
                  gqzz(28) = zr2 * (b(2,1)+yr2*b(2,2))
                  gqzz(29) = yr * zr * (2.0d0*b(2,1)+zr2*b(2,2))
                  gqzz(30) = 2.0d0*b(2,0)
     &                          + zr2*(5.0d0*b(2,1)+zr2*b(2,2))
                  gqxy(25) = xr * yr * (3.0d0*b(2,1) + xr2*b(2,2))
                  gqxy(26) = b(2,0) + (xr2+yr2)*b(2,1) + xr2*yr2*b(2,2)
                  gqxy(27) = zr * yr * (b(2,1)+xr2*b(2,2))
                  gqxy(28) = xr * yr * (3.0d0*b(2,1)+yr2*b(2,2))
                  gqxy(29) = zr * xr * (b(2,1)+yr2*b(2,2))
                  gqxy(30) = xr * yr * (b(2,1)+zr2*b(2,2))
                  gqxz(25) = xr * zr * (3.0d0*b(2,1)+xr2*b(2,2))
                  gqxz(26) = yr * zr * (b(2,1)+xr2*b(2,2))
                  gqxz(27) = b(2,0) + (xr2+zr2)*b(2,1) + xr2*zr2*b(2,2)
                  gqxz(28) = xr * zr * (b(2,1)+yr2*b(2,2))
                  gqxz(29) = xr * yr * (b(2,1)+zr2*b(2,2))
                  gqxz(30) = xr * zr * (3.0d0*b(2,1)+zr2*b(2,2))
                  gqyz(25) = zr * yr * (b(2,1)+xr2*b(2,2))
                  gqyz(26) = xr * zr * (b(2,1)+yr2*b(2,2))
                  gqyz(27) = xr * yr * (b(2,1)+zr2*b(2,2))
                  gqyz(28) = yr * zr * (3.0d0*b(2,1)+yr2*b(2,2))
                  gqyz(29) = b(2,0) + (yr2+zr2)*b(2,1) + yr2*zr2*b(2,2)
                  gqyz(30) = yr * zr * (3.0d0*b(2,1)+zr2*b(2,2))
c
c     unweighted 3rd reaction potential gradient tensor
c
                  gc(11) = xr * (3.0d0*a(0,2)+xr2*a(0,3))
                  gc(12) = yr * (a(0,2)+xr2*a(0,3))
                  gc(13) = zr * (a(0,2)+xr2*a(0,3))
                  gc(14) = xr * (a(0,2)+yr2*a(0,3))
                  gc(15) = xr * yr * zr * a(0,3)
                  gc(16) = xr * (a(0,2)+zr2*a(0,3))
                  gc(17) = yr * (3.0d0*a(0,2)+yr2*a(0,3))
                  gc(18) = zr * (a(0,2)+yr2*a(0,3))
                  gc(19) = yr * (a(0,2)+zr2*a(0,3))
                  gc(20) = zr * (3.0d0*a(0,2)+zr2*a(0,3))
                  gux(11) = 3.0d0*a(1,1) + xr2*(6.0d0*a(1,2)+xr2*a(1,3))
                  gux(12) = xr * yr * (3.0d0*a(1,2)+xr2*a(1,3))
                  gux(13) = xr * zr * (3.0d0*a(1,2)+xr2*a(1,3))
                  gux(14) = a(1,1) + (xr2+yr2)*a(1,2) + xr2*yr2*a(1,3)
                  gux(15) = yr * zr * (a(1,2)+xr2*a(1,3))
                  gux(16) = a(1,1) + (xr2+zr2)*a(1,2) + xr2*zr2*a(1,3)
                  gux(17) = xr * yr * (3.0d0*a(1,2)+yr2*a(1,3))
                  gux(18) = xr * zr * (a(1,2)+yr2*a(1,3))
                  gux(19) = xr * yr * (a(1,2)+zr2*a(1,3))
                  gux(20) = xr * zr * (3.0d0*a(1,2)+zr2*a(1,3))
                  guy(11) = gux(12)
                  guy(12) = gux(14)
                  guy(13) = gux(15)
                  guy(14) = gux(17)
                  guy(15) = gux(18)
                  guy(16) = gux(19)
                  guy(17) = 3.0d0*a(1,1) + yr2*(6.0d0*a(1,2)+yr2*a(1,3))
                  guy(18) = yr * zr * (3.0d0*a(1,2)+yr2*a(1,3))
                  guy(19) = a(1,1) + (yr2+zr2)*a(1,2) + yr2*zr2*a(1,3)
                  guy(20) = yr * zr * (3.0d0*a(1,2)+zr2*a(1,3))
                  guz(11) = gux(13)
                  guz(12) = gux(15)
                  guz(13) = gux(16)
                  guz(14) = gux(18)
                  guz(15) = gux(19)
                  guz(16) = gux(20)
                  guz(17) = guy(18)
                  guz(18) = guy(19)
                  guz(19) = guy(20)
                  guz(20) = 3.0d0*a(1,1) + zr2*(6.0d0*a(1,2)+zr2*a(1,3))
                  gqxx(11) = xr * (12.0d0*a(2,1)+xr2*(9.0d0*a(2,2)
     &                                +xr2*a(2,3)))
                  gqxx(12) = yr * (2.0d0*a(2,1)+xr2*(5.0d0*a(2,2)
     &                                +xr2*a(2,3)))
                  gqxx(13) = zr * (2.0d0*a(2,1)+xr2*(5.0d0*a(2,2)
     &                                +xr2*a(2,3)))
                  gqxx(14) = xr * (2.0d0*a(2,1)+yr2*2.0d0*a(2,2)
     &                                +xr2*(a(2,2)+yr2*a(2,3)))
                  gqxx(15) = xr * yr * zr * (2.0d0*a(2,2)+xr2*a(2,3))
                  gqxx(16) = xr * (2.0d0*a(2,1)+zr2*2.0d0*a(2,2)
     &                                +xr2*(a(2,2)+zr2*a(2,3)))
                  gqxx(17) = yr * xr2 * (3.0d0*a(2,2)+yr2*a(2,3))
                  gqxx(18) = zr * xr2 * (a(2,2)+yr2*a(2,3))
                  gqxx(19) = yr * xr2 * (a(2,2)+zr2*a(2,3))
                  gqxx(20) = zr * xr2 * (3.0d0*a(2,2)+zr2*a(2,3))
                  gqxy(11) = yr * (3.0d0*a(2,1)+xr2*(6.0d0*a(2,2)
     &                                +xr2*a(2,3)))
                  gqxy(12) = xr * (3.0d0*(a(2,1)+yr2*a(2,2))
     &                                +xr2*(a(2,2)+yr2*a(2,3)))
                  gqxy(13) = xr * yr * zr * (3.0d0*a(2,2)+xr2*a(2,3))
                  gqxy(14) = yr * (3.0d0*(a(2,1)+xr2*a(2,2))
     &                                +yr2*(a(2,2)+xr2*a(2,3)))
                  gqxy(15) = zr * (a(2,1)+(yr2+xr2)*a(2,2)
     &                                +yr2*xr2*a(2,3))
                  gqxy(16) = yr * (a(2,1)+(xr2+zr2)*a(2,2)
     &                                +xr2*zr2*a(2,3))
                  gqxy(17) = xr * (3.0d0*(a(2,1)+yr2*a(2,2))
     &                                +yr2*(3.0d0*a(2,2)+yr2*a(2,3)))
                  gqxy(18) = xr * yr * zr * (3.0d0*a(2,2)+yr2*a(2,3))
                  gqxy(19) = xr * (a(2,1)+(yr2+zr2)*a(2,2)
     &                                +yr2*zr2*a(2,3))
                  gqxy(20) = xr * yr * zr * (3.0d0*a(2,2)+zr2*a(2,3))
                  gqxz(11) = zr * (3.0d0*a(2,1)+xr2*(6.0d0*a(2,2)
     &                                +xr2*a(2,3)))
                  gqxz(12) = xr * yr * zr * (3.0d0*a(2,2)+xr2*a(2,3))
                  gqxz(13) = xr * (3.0d0*(a(2,1)+zr2*a(2,2))
     &                                +xr2*(a(2,2)+zr2*a(2,3)))
                  gqxz(14) = zr * (a(2,1)+(xr2+yr2)*a(2,2)
     &                                +xr2*yr2*a(2,3))
                  gqxz(15) = yr * (a(2,1)+(xr2+zr2)*a(2,2)
     &                                +zr2*xr2*a(2,3))
                  gqxz(16) = zr * (3.0d0*(a(2,1)+xr2*a(2,2))
     &                                +zr2*(a(2,2)+xr2*a(2,3)))
                  gqxz(17) = xr * yr * zr * (3.0d0*a(2,2)+yr2*a(2,3))
                  gqxz(18) = xr * (a(2,1)+(zr2+yr2)*a(2,2)
     &                                +zr2*yr2*a(2,3))
                  gqxz(19) = xr * yr * zr * (3.0d0*a(2,2)+zr2*a(2,3))
                  gqxz(20) = xr * (3.0d0*a(2,1)+zr2*(6.0d0*a(2,2)
     &                                +zr2*a(2,3)))
                  gqyy(11) = xr * yr2 * (3.0d0*a(2,2)+xr2*a(2,3))
                  gqyy(12) = yr * (2.0d0*a(2,1)+xr2*2.0d0*a(2,2)
     &                                +yr2*(a(2,2)+xr2*a(2,3)))
                  gqyy(13) = zr * yr2 * (a(2,2)+xr2*a(2,3))
                  gqyy(14) = xr * (2.0d0*a(2,1)+yr2*(5.0d0*a(2,2)
     &                                +yr2*a(2,3)))
                  gqyy(15) = xr * yr * zr * (2.0d0*a(2,2)+yr2*a(2,3))
                  gqyy(16) = xr * yr2 * (a(2,2)+zr2*a(2,3))
                  gqyy(17) = yr * (12.0d0*a(2,1)+yr2*(9.0d0*a(2,2)
     &                                +yr2*a(2,3)))
                  gqyy(18) = zr * (2.0d0*a(2,1)+yr2*(5.0d0*a(2,2)
     &                                +yr2*a(2,3)))
                  gqyy(19) = yr * (2.0d0*a(2,1)+zr2*2.0d0*a(2,2)
     &                                +yr2*(a(2,2)+zr2*a(2,3)))
                  gqyy(20) = zr * yr2 * (3.0d0*a(2,2)+zr2*a(2,3))
                  gqyz(11) = xr * yr * zr * (3.0d0*a(2,2)+xr2*a(2,3))
                  gqyz(12) = zr * (a(2,1)+(xr2+yr2)*a(2,2)
     &                                +xr2*yr2*a(2,3))
                  gqyz(13) = yr * (a(2,1)+(xr2+zr2)*a(2,2)
     &                                +xr2*zr2*a(2,3))
                  gqyz(14) = xr * yr * zr * (3.0d0*a(2,2)+yr2*a(2,3))
                  gqyz(15) = xr * (a(2,1)+(yr2+zr2)*a(2,2)
     &                                +yr2*zr2*a(2,3))
                  gqyz(16) = xr * yr * zr * (3.0d0*a(2,2)+zr2*a(2,3))
                  gqyz(17) = zr * (3.0d0*a(2,1)+yr2*(6.0d0*a(2,2)
     &                                +yr2*a(2,3)))
                  gqyz(18) = yr * (3.0d0*(a(2,1)+zr2*a(2,2))
     &                                +yr2*(a(2,2)+zr2*a(2,3)))
                  gqyz(19) = zr * (3.0d0*(a(2,1)+yr2*a(2,2))
     &                                +zr2*(a(2,2)+yr2*a(2,3)))
                  gqyz(20) = yr * (3.0d0*a(2,1)+zr2*(6.0d0*a(2,2)
     &                                +zr2*a(2,3)))
                  gqzz(11) = xr * zr2 * (3.0d0*a(2,2)+xr2*a(2,3))
                  gqzz(12) = yr * (zr2*a(2,2)+xr2*(zr2*a(2,3)))
                  gqzz(13) = zr * (2.0d0*a(2,1)+xr2*2.0d0*a(2,2)
     &                                +zr2*(a(2,2)+xr2*a(2,3)))
                  gqzz(14) = xr * zr2 * (a(2,2)+yr2*a(2,3))
                  gqzz(15) = xr * yr * zr * (2.0d0*a(2,2)+zr2*a(2,3))
                  gqzz(16) = xr * (2.0d0*a(2,1)+zr2*(5.0d0*a(2,2)
     &                                +zr2*a(2,3)))
                  gqzz(17) = yr * zr2 * (3.0d0*a(2,2)+yr2*a(2,3))
                  gqzz(18) = zr * (2.0d0*a(2,1)+yr2*2.0d0*a(2,2)
     &                                +zr2*(a(2,2)+yr2*a(2,3)))
                  gqzz(19) = yr * (2.0d0*a(2,1)+zr2*(5.0d0*a(2,2)
     &                                +zr2*a(2,3)))
                  gqzz(20) = zr * (12.0d0*a(2,1)+zr2*(9.0d0*a(2,2)
     &                                +zr2*a(2,3)))
c
c     electrostatic solvation energy of the permanent multipoles
c     in their own GK reaction potential
c
                  esym = ci * ck * gc(1)
     &                   - (uxi*(uxk*gux(2)+uyk*guy(2)+uzk*guz(2))
     &                     +uyi*(uxk*gux(3)+uyk*guy(3)+uzk*guz(3))
     &                     +uzi*(uxk*gux(4)+uyk*guy(4)+uzk*guz(4)))
                  ewi = ci*(uxk*gc(2)+uyk*gc(3)+uzk*gc(4))
     &                 -ck*(uxi*gux(1)+uyi*guy(1)+uzi*guz(1))
     &                 +ci*(qxxk*gc(5)+qyyk*gc(8)+qzzk*gc(10)
     &                 +2.0d0*(qxyk*gc(6)+qxzk*gc(7)+qyzk*gc(9)))
     &                 +ck*(qxxi*gqxx(1)+qyyi*gqyy(1)+qzzi*gqzz(1)
     &                 +2.0d0*(qxyi*gqxy(1)+qxzi*gqxz(1)+qyzi*gqyz(1)))
     &                 - uxi*(qxxk*gux(5)+qyyk*gux(8)+qzzk*gux(10)
     &                 +2.0d0*(qxyk*gux(6)+qxzk*gux(7)+qyzk*gux(9)))
     &                 - uyi*(qxxk*guy(5)+qyyk*guy(8)+qzzk*guy(10)
     &                 +2.0d0*(qxyk*guy(6)+qxzk*guy(7)+qyzk*guy(9)))
     &                 - uzi*(qxxk*guz(5)+qyyk*guz(8)+qzzk*guz(10)
     &                 +2.0d0*(qxyk*guz(6)+qxzk*guz(7)+qyzk*guz(9)))
     &                 + uxk*(qxxi*gqxx(2)+qyyi*gqyy(2)+qzzi*gqzz(2)
     &                 +2.0d0*(qxyi*gqxy(2)+qxzi*gqxz(2)+qyzi*gqyz(2)))
     &                 + uyk*(qxxi*gqxx(3)+qyyi*gqyy(3)+qzzi*gqzz(3)
     &                 +2.0d0*(qxyi*gqxy(3)+qxzi*gqxz(3)+qyzi*gqyz(3)))
     &                 + uzk*(qxxi*gqxx(4)+qyyi*gqyy(4)+qzzi*gqzz(4)
     &                 +2.0d0*(qxyi*gqxy(4)+qxzi*gqxz(4)+qyzi*gqyz(4)))
     &                 + qxxi*(qxxk*gqxx(5)+qyyk*gqxx(8)+qzzk*gqxx(10)
     &                 +2.0d0*(qxyk*gqxx(6)+qxzk*gqxx(7)+qyzk*gqxx(9)))
     &                 + qyyi*(qxxk*gqyy(5)+qyyk*gqyy(8)+qzzk*gqyy(10)
     &                 +2.0d0*(qxyk*gqyy(6)+qxzk*gqyy(7)+qyzk*gqyy(9)))
     &                 + qzzi*(qxxk*gqzz(5)+qyyk*gqzz(8)+qzzk*gqzz(10)
     &                 +2.0d0*(qxyk*gqzz(6)+qxzk*gqzz(7)+qyzk*gqzz(9)))
     &           + 2.0d0*(qxyi*(qxxk*gqxy(5)+qyyk*gqxy(8)+qzzk*gqxy(10)
     &                 +2.0d0*(qxyk*gqxy(6)+qxzk*gqxy(7)+qyzk*gqxy(9)))
     &                 + qxzi*(qxxk*gqxz(5)+qyyk*gqxz(8)+qzzk*gqxz(10)
     &                 +2.0d0*(qxyk*gqxz(6)+qxzk*gqxz(7)+qyzk*gqxz(9)))
     &                 + qyzi*(qxxk*gqyz(5)+qyyk*gqyz(8)+qzzk*gqyz(10)
     &                 +2.0d0*(qxyk*gqyz(6)+qxzk*gqyz(7)+qyzk*gqyz(9))))
                  ewk = ci*(uxk*gux(1)+uyk*guy(1)+uzk*guz(1))
     &                 -ck*(uxi*gc(2)+uyi*gc(3)+uzi*gc(4))
     &                 +ci*(qxxk*gqxx(1)+qyyk*gqyy(1)+qzzk*gqzz(1)
     &                 +2.0d0*(qxyk*gqxy(1)+qxzk*gqxz(1)+qyzk*gqyz(1)))
     &                 +ck*(qxxi*gc(5)+qyyi*gc(8)+qzzi*gc(10)
     &                 +2.0d0*(qxyi*gc(6)+qxzi*gc(7)+qyzi*gc(9)))
     &                 - uxi*(qxxk*gqxx(2)+qyyk*gqyy(2)+qzzk*gqzz(2)
     &                 +2.0d0*(qxyk*gqxy(2)+qxzk*gqxz(2)+qyzk*gqyz(2)))
     &                 - uyi*(qxxk*gqxx(3)+qyyk*gqyy(3)+qzzk*gqzz(3)
     &                 +2.0d0*(qxyk*gqxy(3)+qxzk*gqxz(3)+qyzk*gqyz(3)))
     &                 - uzi*(qxxk*gqxx(4)+qyyk*gqyy(4)+qzzk*gqzz(4)
     &                 +2.0d0*(qxyk*gqxy(4)+qxzk*gqxz(4)+qyzk*gqyz(4)))
     &                 + uxk*(qxxi*gux(5)+qyyi*gux(8)+qzzi*gux(10)
     &                 +2.0d0*(qxyi*gux(6)+qxzi*gux(7)+qyzi*gux(9)))
     &                 + uyk*(qxxi*guy(5)+qyyi*guy(8)+qzzi*guy(10)
     &                 +2.0d0*(qxyi*guy(6)+qxzi*guy(7)+qyzi*guy(9)))
     &                 + uzk*(qxxi*guz(5)+qyyi*guz(8)+qzzi*guz(10)
     &                 +2.0d0*(qxyi*guz(6)+qxzi*guz(7)+qyzi*guz(9)))
     &                 + qxxi*(qxxk*gqxx(5)+qyyk*gqyy(5)+qzzk*gqzz(5)
     &                 +2.0d0*(qxyk*gqxy(5)+qxzk*gqxz(5)+qyzk*gqyz(5)))
     &                 + qyyi*(qxxk*gqxx(8)+qyyk*gqyy(8)+qzzk*gqzz(8)
     &                 +2.0d0*(qxyk*gqxy(8)+qxzk*gqxz(8)+qyzk*gqyz(8)))
     &                 + qzzi*(qxxk*gqxx(10)+qyyk*gqyy(10)+qzzk*gqzz(10)
     &           +2.0d0*(qxyk*gqxy(10)+qxzk*gqxz(10)+qyzk*gqyz(10)))
     &           + 2.0d0*(qxyi*(qxxk*gqxx(6)+qyyk*gqyy(6)+qzzk*gqzz(6)
     &                 +2.0d0*(qxyk*gqxy(6)+qxzk*gqxz(6)+qyzk*gqyz(6)))
     &                 + qxzi*(qxxk*gqxx(7)+qyyk*gqyy(7)+qzzk*gqzz(7)
     &                 +2.0d0*(qxyk*gqxy(7)+qxzk*gqxz(7)+qyzk*gqyz(7)))
     &                 + qyzi*(qxxk*gqxx(9)+qyyk*gqyy(9)+qzzk*gqzz(9)
     &                 +2.0d0*(qxyk*gqxy(9)+qxzk*gqxz(9)+qyzk*gqyz(9))))
c
                  desymdx = ci * ck * gc(2)
     &                      - (uxi*(uxk*gux(5)+uyk*guy(5)+uzk*guz(5))
     &                        +uyi*(uxk*gux(6)+uyk*guy(6)+uzk*guz(6))
     &                        +uzi*(uxk*gux(7)+uyk*guy(7)+uzk*guz(7)))
                  dewidx = ci*(uxk*gc(5)+uyk*gc(6)+uzk*gc(7))
     &                    -ck*(uxi*gux(2)+uyi*guy(2)+uzi*guz(2))
     &                 +ci*(qxxk*gc(11)+qyyk*gc(14)+qzzk*gc(16)
     &              +2.0d0*(qxyk*gc(12)+qxzk*gc(13)+qyzk*gc(15)))
     &                 +ck*(qxxi*gqxx(2)+qyyi*gqyy(2)+qzzi*gqzz(2)
     &              +2.0d0*(qxyi*gqxy(2)+qxzi*gqxz(2)+qyzi*gqyz(2)))
     &                 - uxi*(qxxk*gux(11)+qyyk*gux(14)+qzzk*gux(16)
     &              +2.0d0*(qxyk*gux(12)+qxzk*gux(13)+qyzk*gux(15)))
     &                 - uyi*(qxxk*guy(11)+qyyk*guy(14)+qzzk*guy(16)
     &              +2.0d0*(qxyk*guy(12)+qxzk*guy(13)+qyzk*guy(15)))
     &                 - uzi*(qxxk*guz(11)+qyyk*guz(14)+qzzk*guz(16)
     &              +2.0d0*(qxyk*guz(12)+qxzk*guz(13)+qyzk*guz(15)))
     &                 + uxk*(qxxi*gqxx(5)+qyyi*gqyy(5)+qzzi*gqzz(5)
     &              +2.0d0*(qxyi*gqxy(5)+qxzi*gqxz(5)+qyzi*gqyz(5)))
     &                 + uyk*(qxxi*gqxx(6)+qyyi*gqyy(6)+qzzi*gqzz(6)
     &              +2.0d0*(qxyi*gqxy(6)+qxzi*gqxz(6)+qyzi*gqyz(6)))
     &                 + uzk*(qxxi*gqxx(7)+qyyi*gqyy(7)+qzzi*gqzz(7)
     &              +2.0d0*(qxyi*gqxy(7)+qxzi*gqxz(7)+qyzi*gqyz(7)))
     &                 + qxxi*(qxxk*gqxx(11)+qyyk*gqxx(14)+qzzk*gqxx(16)
     &              +2.0d0*(qxyk*gqxx(12)+qxzk*gqxx(13)+qyzk*gqxx(15)))
     &                 + qyyi*(qxxk*gqyy(11)+qyyk*gqyy(14)+qzzk*gqyy(16)
     &              +2.0d0*(qxyk*gqyy(12)+qxzk*gqyy(13)+qyzk*gqyy(15)))
     &                 + qzzi*(qxxk*gqzz(11)+qyyk*gqzz(14)+qzzk*gqzz(16)
     &              +2.0d0*(qxyk*gqzz(12)+qxzk*gqzz(13)+qyzk*gqzz(15)))
     &        + 2.0d0*(qxyi*(qxxk*gqxy(11)+qyyk*gqxy(14)+qzzk*gqxy(16)
     &              +2.0d0*(qxyk*gqxy(12)+qxzk*gqxy(13)+qyzk*gqxy(15)))
     &                 + qxzi*(qxxk*gqxz(11)+qyyk*gqxz(14)+qzzk*gqxz(16)
     &              +2.0d0*(qxyk*gqxz(12)+qxzk*gqxz(13)+qyzk*gqxz(15)))
     &                 + qyzi*(qxxk*gqyz(11)+qyyk*gqyz(14)+qzzk*gqyz(16)
     &              +2.0d0*(qxyk*gqyz(12)+qxzk*gqyz(13)+qyzk*gqyz(15))))
                  dewkdx = ci*(uxk*gux(2)+uyk*guy(2)+uzk*guz(2))
     &                    -ck*(uxi*gc(5)+uyi*gc(6)+uzi*gc(7))
     &                 +ci*(qxxk*gqxx(2)+qyyk*gqyy(2)+qzzk*gqzz(2)
     &              +2.0d0*(qxyk*gqxy(2)+qxzk*gqxz(2)+qyzk*gqyz(2)))
     &                 +ck*(qxxi*gc(11)+qyyi*gc(14)+qzzi*gc(16)
     &              +2.0d0*(qxyi*gc(12)+qxzi*gc(13)+qyzi*gc(15)))
     &                 - uxi*(qxxk*gqxx(5)+qyyk*gqyy(5)+qzzk*gqzz(5)
     &              +2.0d0*(qxyk*gqxy(5)+qxzk*gqxz(5)+qyzk*gqyz(5)))
     &                 - uyi*(qxxk*gqxx(6)+qyyk*gqyy(6)+qzzk*gqzz(6)
     &              +2.0d0*(qxyk*gqxy(6)+qxzk*gqxz(6)+qyzk*gqyz(6)))
     &                 - uzi*(qxxk*gqxx(7)+qyyk*gqyy(7)+qzzk*gqzz(7)
     &              +2.0d0*(qxyk*gqxy(7)+qxzk*gqxz(7)+qyzk*gqyz(7)))
     &                 + uxk*(qxxi*gux(11)+qyyi*gux(14)+qzzi*gux(16)
     &              +2.0d0*(qxyi*gux(12)+qxzi*gux(13)+qyzi*gux(15)))
     &                 + uyk*(qxxi*guy(11)+qyyi*guy(14)+qzzi*guy(16)
     &              +2.0d0*(qxyi*guy(12)+qxzi*guy(13)+qyzi*guy(15)))
     &                 + uzk*(qxxi*guz(11)+qyyi*guz(14)+qzzi*guz(16)
     &              +2.0d0*(qxyi*guz(12)+qxzi*guz(13)+qyzi*guz(15)))
     &                 + qxxi*(qxxk*gqxx(11)+qyyk*gqyy(11)+qzzk*gqzz(11)
     &              +2.0d0*(qxyk*gqxy(11)+qxzk*gqxz(11)+qyzk*gqyz(11)))
     &                 + qyyi*(qxxk*gqxx(14)+qyyk*gqyy(14)+qzzk*gqzz(14)
     &              +2.0d0*(qxyk*gqxy(14)+qxzk*gqxz(14)+qyzk*gqyz(14)))
     &                 + qzzi*(qxxk*gqxx(16)+qyyk*gqyy(16)+qzzk*gqzz(16)
     &              +2.0d0*(qxyk*gqxy(16)+qxzk*gqxz(16)+qyzk*gqyz(16)))
     &        + 2.0d0*(qxyi*(qxxk*gqxx(12)+qyyk*gqyy(12)+qzzk*gqzz(12)
     &              +2.0d0*(qxyk*gqxy(12)+qxzk*gqxz(12)+qyzk*gqyz(12)))
     &                 + qxzi*(qxxk*gqxx(13)+qyyk*gqyy(13)+qzzk*gqzz(13)
     &              +2.0d0*(qxyk*gqxy(13)+qxzk*gqxz(13)+qyzk*gqyz(13)))
     &                 + qyzi*(qxxk*gqxx(15)+qyyk*gqyy(15)+qzzk*gqzz(15)
     &              +2.0d0*(qxyk*gqxy(15)+qxzk*gqxz(15)+qyzk*gqyz(15))))
                  dedx = desymdx + 0.5d0*(dewidx+dewkdx)
c
                  desymdy = ci * ck * gc(3)
     &                      - (uxi*(uxk*gux(6)+uyk*guy(6)+uzk*guz(6))
     &                        +uyi*(uxk*gux(8)+uyk*guy(8)+uzk*guz(8))
     &                        +uzi*(uxk*gux(9)+uyk*guy(9)+uzk*guz(9)))
                  dewidy = ci*(uxk*gc(6)+uyk*gc(8)+uzk*gc(9))
     &                    -ck*(uxi*gux(3)+uyi*guy(3)+uzi*guz(3))
     &                 +ci*(qxxk*gc(12)+qyyk*gc(17)+qzzk*gc(19)
     &              +2.0d0*(qxyk*gc(14)+qxzk*gc(15)+qyzk*gc(18)))
     &                 +ck*(qxxi*gqxx(3)+qyyi*gqyy(3)+qzzi*gqzz(3)
     &              +2.0d0*(qxyi*gqxy(3)+qxzi*gqxz(3)+qyzi*gqyz(3)))
     &                 - uxi*(qxxk*gux(12)+qyyk*gux(17)+qzzk*gux(19)
     &              +2.0d0*(qxyk*gux(14)+qxzk*gux(15)+qyzk*gux(18)))
     &                 - uyi*(qxxk*guy(12)+qyyk*guy(17)+qzzk*guy(19)
     &              +2.0d0*(qxyk*guy(14)+qxzk*guy(15)+qyzk*guy(18)))
     &                 - uzi*(qxxk*guz(12)+qyyk*guz(17)+qzzk*guz(19)
     &              +2.0d0*(qxyk*guz(14)+qxzk*guz(15)+qyzk*guz(18)))
     &                 + uxk*(qxxi*gqxx(6)+qyyi*gqyy(6)+qzzi*gqzz(6)
     &              +2.0d0*(qxyi*gqxy(6)+qxzi*gqxz(6)+qyzi*gqyz(6)))
     &                 + uyk*(qxxi*gqxx(8)+qyyi*gqyy(8)+qzzi*gqzz(8)
     &              +2.0d0*(qxyi*gqxy(8)+qxzi*gqxz(8)+qyzi*gqyz(8)))
     &                 + uzk*(qxxi*gqxx(9)+qyyi*gqyy(9)+qzzi*gqzz(9)
     &              +2.0d0*(qxyi*gqxy(9)+qxzi*gqxz(9)+qyzi*gqyz(9)))
     &                 + qxxi*(qxxk*gqxx(12)+qyyk*gqxx(17)+qzzk*gqxx(19)
     &              +2.0d0*(qxyk*gqxx(14)+qxzk*gqxx(15)+qyzk*gqxx(18)))
     &                 + qyyi*(qxxk*gqyy(12)+qyyk*gqyy(17)+qzzk*gqyy(19)
     &              +2.0d0*(qxyk*gqyy(14)+qxzk*gqyy(15)+qyzk*gqyy(18)))
     &                 + qzzi*(qxxk*gqzz(12)+qyyk*gqzz(17)+qzzk*gqzz(19)
     &              +2.0d0*(qxyk*gqzz(14)+qxzk*gqzz(15)+qyzk*gqzz(18)))
     &        + 2.0d0*(qxyi*(qxxk*gqxy(12)+qyyk*gqxy(17)+qzzk*gqxy(19)
     &              +2.0d0*(qxyk*gqxy(14)+qxzk*gqxy(15)+qyzk*gqxy(18)))
     &                 + qxzi*(qxxk*gqxz(12)+qyyk*gqxz(17)+qzzk*gqxz(19)
     &              +2.0d0*(qxyk*gqxz(14)+qxzk*gqxz(15)+qyzk*gqxz(18)))
     &                 + qyzi*(qxxk*gqyz(12)+qyyk*gqyz(17)+qzzk*gqyz(19)
     &              +2.0d0*(qxyk*gqyz(14)+qxzk*gqyz(15)+qyzk*gqyz(18))))
                  dewkdy = ci*(uxk*gux(3)+uyk*guy(3)+uzk*guz(3))
     &                    -ck*(uxi*gc(6)+uyi*gc(8)+uzi*gc(9))
     &                 +ci*(qxxk*gqxx(3)+qyyk*gqyy(3)+qzzk*gqzz(3)
     &              +2.0d0*(qxyk*gqxy(3)+qxzk*gqxz(3)+qyzk*gqyz(3)))
     &                 +ck*(qxxi*gc(12)+qyyi*gc(17)+qzzi*gc(19)
     &              +2.0d0*(qxyi*gc(14)+qxzi*gc(15)+qyzi*gc(18)))
     &                 - uxi*(qxxk*gqxx(6)+qyyk*gqyy(6)+qzzk*gqzz(6)
     &              +2.0d0*(qxyk*gqxy(6)+qxzk*gqxz(6)+qyzk*gqyz(6)))
     &                 - uyi*(qxxk*gqxx(8)+qyyk*gqyy(8)+qzzk*gqzz(8)
     &              +2.0d0*(qxyk*gqxy(8)+qxzk*gqxz(8)+qyzk*gqyz(8)))
     &                 - uzi*(qxxk*gqxx(9)+qyyk*gqyy(9)+qzzk*gqzz(9)
     &              +2.0d0*(qxyk*gqxy(9)+qxzk*gqxz(9)+qyzk*gqyz(9)))
     &                 + uxk*(qxxi*gux(12)+qyyi*gux(17)+qzzi*gux(19)
     &              +2.0d0*(qxyi*gux(14)+qxzi*gux(15)+qyzi*gux(18)))
     &                 + uyk*(qxxi*guy(12)+qyyi*guy(17)+qzzi*guy(19)
     &              +2.0d0*(qxyi*guy(14)+qxzi*guy(15)+qyzi*guy(18)))
     &                 + uzk*(qxxi*guz(12)+qyyi*guz(17)+qzzi*guz(19)
     &              +2.0d0*(qxyi*guz(14)+qxzi*guz(15)+qyzi*guz(18)))
     &                 + qxxi*(qxxk*gqxx(12)+qyyk*gqyy(12)+qzzk*gqzz(12)
     &              +2.0d0*(qxyk*gqxy(12)+qxzk*gqxz(12)+qyzk*gqyz(12)))
     &                 + qyyi*(qxxk*gqxx(17)+qyyk*gqyy(17)+qzzk*gqzz(17)
     &              +2.0d0*(qxyk*gqxy(17)+qxzk*gqxz(17)+qyzk*gqyz(17)))
     &                 + qzzi*(qxxk*gqxx(19)+qyyk*gqyy(19)+qzzk*gqzz(19)
     &              +2.0d0*(qxyk*gqxy(19)+qxzk*gqxz(19)+qyzk*gqyz(19)))
     &        + 2.0d0*(qxyi*(qxxk*gqxx(14)+qyyk*gqyy(14)+qzzk*gqzz(14)
     &              +2.0d0*(qxyk*gqxy(14)+qxzk*gqxz(14)+qyzk*gqyz(14)))
     &                 + qxzi*(qxxk*gqxx(15)+qyyk*gqyy(15)+qzzk*gqzz(15)
     &              +2.0d0*(qxyk*gqxy(15)+qxzk*gqxz(15)+qyzk*gqyz(15)))
     &                 + qyzi*(qxxk*gqxx(18)+qyyk*gqyy(18)+qzzk*gqzz(18)
     &              +2.0d0*(qxyk*gqxy(18)+qxzk*gqxz(18)+qyzk*gqyz(18))))
                  dedy = desymdy + 0.5d0*(dewidy+dewkdy)
c
                  desymdz = ci * ck * gc(4)
     &                      - (uxi*(uxk*gux(7)+uyk*guy(7)+uzk*guz(7))
     &                        +uyi*(uxk*gux(9)+uyk*guy(9)+uzk*guz(9))
     &                       +uzi*(uxk*gux(10)+uyk*guy(10)+uzk*guz(10)))
                  dewidz = ci*(uxk*gc(7)+uyk*gc(9)+uzk*gc(10))
     &                    -ck*(uxi*gux(4)+uyi*guy(4)+uzi*guz(4))
     &                 +ci*(qxxk*gc(13)+qyyk*gc(18)+qzzk*gc(20)
     &              +2.0d0*(qxyk*gc(15)+qxzk*gc(16)+qyzk*gc(19)))
     &                 +ck*(qxxi*gqxx(4)+qyyi*gqyy(4)+qzzi*gqzz(4)
     &              +2.0d0*(qxyi*gqxy(4)+qxzi*gqxz(4)+qyzi*gqyz(4)))
     &                 - uxi*(qxxk*gux(13)+qyyk*gux(18)+qzzk*gux(20)
     &              +2.0d0*(qxyk*gux(15)+qxzk*gux(16)+qyzk*gux(19)))
     &                 - uyi*(qxxk*guy(13)+qyyk*guy(18)+qzzk*guy(20)
     &              +2.0d0*(qxyk*guy(15)+qxzk*guy(16)+qyzk*guy(19)))
     &                 - uzi*(qxxk*guz(13)+qyyk*guz(18)+qzzk*guz(20)
     &              +2.0d0*(qxyk*guz(15)+qxzk*guz(16)+qyzk*guz(19)))
     &                 + uxk*(qxxi*gqxx(7)+qyyi*gqyy(7)+qzzi*gqzz(7)
     &              +2.0d0*(qxyi*gqxy(7)+qxzi*gqxz(7)+qyzi*gqyz(7)))
     &                 + uyk*(qxxi*gqxx(9)+qyyi*gqyy(9)+qzzi*gqzz(9)
     &              +2.0d0*(qxyi*gqxy(9)+qxzi*gqxz(9)+qyzi*gqyz(9)))
     &                 + uzk*(qxxi*gqxx(10)+qyyi*gqyy(10)+qzzi*gqzz(10)
     &              +2.0d0*(qxyi*gqxy(10)+qxzi*gqxz(10)+qyzi*gqyz(10)))
     &                 + qxxi*(qxxk*gqxx(13)+qyyk*gqxx(18)+qzzk*gqxx(20)
     &              +2.0d0*(qxyk*gqxx(15)+qxzk*gqxx(16)+qyzk*gqxx(19)))
     &                 + qyyi*(qxxk*gqyy(13)+qyyk*gqyy(18)+qzzk*gqyy(20)
     &              +2.0d0*(qxyk*gqyy(15)+qxzk*gqyy(16)+qyzk*gqyy(19)))
     &                 + qzzi*(qxxk*gqzz(13)+qyyk*gqzz(18)+qzzk*gqzz(20)
     &              +2.0d0*(qxyk*gqzz(15)+qxzk*gqzz(16)+qyzk*gqzz(19)))
     &        + 2.0d0*(qxyi*(qxxk*gqxy(13)+qyyk*gqxy(18)+qzzk*gqxy(20)
     &              +2.0d0*(qxyk*gqxy(15)+qxzk*gqxy(16)+qyzk*gqxy(19)))
     &                 + qxzi*(qxxk*gqxz(13)+qyyk*gqxz(18)+qzzk*gqxz(20)
     &              +2.0d0*(qxyk*gqxz(15)+qxzk*gqxz(16)+qyzk*gqxz(19)))
     &                 + qyzi*(qxxk*gqyz(13)+qyyk*gqyz(18)+qzzk*gqyz(20)
     &              +2.0d0*(qxyk*gqyz(15)+qxzk*gqyz(16)+qyzk*gqyz(19))))
                  dewkdz = ci*(uxk*gux(4)+uyk*guy(4)+uzk*guz(4))
     &                    -ck*(uxi*gc(7)+uyi*gc(9)+uzi*gc(10))
     &                 +ci*(qxxk*gqxx(4)+qyyk*gqyy(4)+qzzk*gqzz(4)
     &              +2.0d0*(qxyk*gqxy(4)+qxzk*gqxz(4)+qyzk*gqyz(4)))
     &                 +ck*(qxxi*gc(13)+qyyi*gc(18)+qzzi*gc(20)
     &              +2.0d0*(qxyi*gc(15)+qxzi*gc(16)+qyzi*gc(19)))
     &                 - uxi*(qxxk*gqxx(7)+qyyk*gqyy(7)+qzzk*gqzz(7)
     &              +2.0d0*(qxyk*gqxy(7)+qxzk*gqxz(7)+qyzk*gqyz(7)))
     &                 - uyi*(qxxk*gqxx(9)+qyyk*gqyy(9)+qzzk*gqzz(9)
     &              +2.0d0*(qxyk*gqxy(9)+qxzk*gqxz(9)+qyzk*gqyz(9)))
     &                 - uzi*(qxxk*gqxx(10)+qyyk*gqyy(10)+qzzk*gqzz(10)
     &              +2.0d0*(qxyk*gqxy(10)+qxzk*gqxz(10)+qyzk*gqyz(10)))
     &                 + uxk*(qxxi*gux(13)+qyyi*gux(18)+qzzi*gux(20)
     &              +2.0d0*(qxyi*gux(15)+qxzi*gux(16)+qyzi*gux(19)))
     &                 + uyk*(qxxi*guy(13)+qyyi*guy(18)+qzzi*guy(20)
     &              +2.0d0*(qxyi*guy(15)+qxzi*guy(16)+qyzi*guy(19)))
     &                 + uzk*(qxxi*guz(13)+qyyi*guz(18)+qzzi*guz(20)
     &              +2.0d0*(qxyi*guz(15)+qxzi*guz(16)+qyzi*guz(19)))
     &                 + qxxi*(qxxk*gqxx(13)+qyyk*gqyy(13)+qzzk*gqzz(13)
     &              +2.0d0*(qxyk*gqxy(13)+qxzk*gqxz(13)+qyzk*gqyz(13)))
     &                 + qyyi*(qxxk*gqxx(18)+qyyk*gqyy(18)+qzzk*gqzz(18)
     &              +2.0d0*(qxyk*gqxy(18)+qxzk*gqxz(18)+qyzk*gqyz(18)))
     &                 + qzzi*(qxxk*gqxx(20)+qyyk*gqyy(20)+qzzk*gqzz(20)
     &              +2.0d0*(qxyk*gqxy(20)+qxzk*gqxz(20)+qyzk*gqyz(20)))
     &        + 2.0d0*(qxyi*(qxxk*gqxx(15)+qyyk*gqyy(15)+qzzk*gqzz(15)
     &              +2.0d0*(qxyk*gqxy(15)+qxzk*gqxz(15)+qyzk*gqyz(15)))
     &                 + qxzi*(qxxk*gqxx(16)+qyyk*gqyy(16)+qzzk*gqzz(16)
     &              +2.0d0*(qxyk*gqxy(16)+qxzk*gqxz(16)+qyzk*gqyz(16)))
     &                 + qyzi*(qxxk*gqxx(19)+qyyk*gqyy(19)+qzzk*gqzz(19)
     &              +2.0d0*(qxyk*gqxy(19)+qxzk*gqxz(19)+qyzk*gqyz(19))))
                  dedz = desymdz + 0.5d0*(dewidz+dewkdz)
c
                  desymdr = ci * ck * gc(21)
     &                      - (uxi*(uxk*gux(22)+uyk*guy(22)+uzk*guz(22))
     &                        +uyi*(uxk*gux(23)+uyk*guy(23)+uzk*guz(23))
     &                       +uzi*(uxk*gux(24)+uyk*guy(24)+uzk*guz(24)))
                  dewidr = ci*(uxk*gc(22)+uyk*gc(23)+uzk*gc(24))
     &                    -ck*(uxi*gux(21)+uyi*guy(21)+uzi*guz(21))
     &                 +ci*(qxxk*gc(25)+qyyk*gc(28)+qzzk*gc(30)
     &              +2.0d0*(qxyk*gc(26)+qxzk*gc(27)+qyzk*gc(29)))
     &                 +ck*(qxxi*gqxx(21)+qyyi*gqyy(21)+qzzi*gqzz(21)
     &              +2.0d0*(qxyi*gqxy(21)+qxzi*gqxz(21)+qyzi*gqyz(21)))
     &                 - uxi*(qxxk*gux(25)+qyyk*gux(28)+qzzk*gux(30)
     &              +2.0d0*(qxyk*gux(26)+qxzk*gux(27)+qyzk*gux(29)))
     &                 - uyi*(qxxk*guy(25)+qyyk*guy(28)+qzzk*guy(30)
     &              +2.0d0*(qxyk*guy(26)+qxzk*guy(27)+qyzk*guy(29)))
     &                 - uzi*(qxxk*guz(25)+qyyk*guz(28)+qzzk*guz(30)
     &              +2.0d0*(qxyk*guz(26)+qxzk*guz(27)+qyzk*guz(29)))
     &                 + uxk*(qxxi*gqxx(22)+qyyi*gqyy(22)+qzzi*gqzz(22)
     &              +2.0d0*(qxyi*gqxy(22)+qxzi*gqxz(22)+qyzi*gqyz(22)))
     &                 + uyk*(qxxi*gqxx(23)+qyyi*gqyy(23)+qzzi*gqzz(23)
     &              +2.0d0*(qxyi*gqxy(23)+qxzi*gqxz(23)+qyzi*gqyz(23)))
     &                 + uzk*(qxxi*gqxx(24)+qyyi*gqyy(24)+qzzi*gqzz(24)
     &              +2.0d0*(qxyi*gqxy(24)+qxzi*gqxz(24)+qyzi*gqyz(24)))
     &                 + qxxi*(qxxk*gqxx(25)+qyyk*gqxx(28)+qzzk*gqxx(30)
     &              +2.0d0*(qxyk*gqxx(26)+qxzk*gqxx(27)+qyzk*gqxx(29)))
     &                 + qyyi*(qxxk*gqyy(25)+qyyk*gqyy(28)+qzzk*gqyy(30)
     &              +2.0d0*(qxyk*gqyy(26)+qxzk*gqyy(27)+qyzk*gqyy(29)))
     &                 + qzzi*(qxxk*gqzz(25)+qyyk*gqzz(28)+qzzk*gqzz(30)
     &              +2.0d0*(qxyk*gqzz(26)+qxzk*gqzz(27)+qyzk*gqzz(29)))
     &        + 2.0d0*(qxyi*(qxxk*gqxy(25)+qyyk*gqxy(28)+qzzk*gqxy(30)
     &              +2.0d0*(qxyk*gqxy(26)+qxzk*gqxy(27)+qyzk*gqxy(29)))
     &                 + qxzi*(qxxk*gqxz(25)+qyyk*gqxz(28)+qzzk*gqxz(30)
     &              +2.0d0*(qxyk*gqxz(26)+qxzk*gqxz(27)+qyzk*gqxz(29)))
     &                 + qyzi*(qxxk*gqyz(25)+qyyk*gqyz(28)+qzzk*gqyz(30)
     &              +2.0d0*(qxyk*gqyz(26)+qxzk*gqyz(27)+qyzk*gqyz(29))))
                  dewkdr = ci*(uxk*gux(21)+uyk*guy(21)+uzk*guz(21))
     &                    -ck*(uxi*gc(22)+uyi*gc(23)+uzi*gc(24))
     &                 +ci*(qxxk*gqxx(21)+qyyk*gqyy(21)+qzzk*gqzz(21)
     &              +2.0d0*(qxyk*gqxy(21)+qxzk*gqxz(21)+qyzk*gqyz(21)))
     &                 +ck*(qxxi*gc(25)+qyyi*gc(28)+qzzi*gc(30)
     &              +2.0d0*(qxyi*gc(26)+qxzi*gc(27)+qyzi*gc(29)))
     &                 - uxi*(qxxk*gqxx(22)+qyyk*gqyy(22)+qzzk*gqzz(22)
     &              +2.0d0*(qxyk*gqxy(22)+qxzk*gqxz(22)+qyzk*gqyz(22)))
     &                 - uyi*(qxxk*gqxx(23)+qyyk*gqyy(23)+qzzk*gqzz(23)
     &              +2.0d0*(qxyk*gqxy(23)+qxzk*gqxz(23)+qyzk*gqyz(23)))
     &                 - uzi*(qxxk*gqxx(24)+qyyk*gqyy(24)+qzzk*gqzz(24)
     &              +2.0d0*(qxyk*gqxy(24)+qxzk*gqxz(24)+qyzk*gqyz(24)))
     &                 + uxk*(qxxi*gux(25)+qyyi*gux(28)+qzzi*gux(30)
     &              +2.0d0*(qxyi*gux(26)+qxzi*gux(27)+qyzi*gux(29)))
     &                 + uyk*(qxxi*guy(25)+qyyi*guy(28)+qzzi*guy(30)
     &              +2.0d0*(qxyi*guy(26)+qxzi*guy(27)+qyzi*guy(29)))
     &                 + uzk*(qxxi*guz(25)+qyyi*guz(28)+qzzi*guz(30)
     &              +2.0d0*(qxyi*guz(26)+qxzi*guz(27)+qyzi*guz(29)))
     &                 + qxxi*(qxxk*gqxx(25)+qyyk*gqyy(25)+qzzk*gqzz(25)
     &              +2.0d0*(qxyk*gqxy(25)+qxzk*gqxz(25)+qyzk*gqyz(25)))
     &                 + qyyi*(qxxk*gqxx(28)+qyyk*gqyy(28)+qzzk*gqzz(28)
     &              +2.0d0*(qxyk*gqxy(28)+qxzk*gqxz(28)+qyzk*gqyz(28)))
     &                 + qzzi*(qxxk*gqxx(30)+qyyk*gqyy(30)+qzzk*gqzz(30)
     &              +2.0d0*(qxyk*gqxy(30)+qxzk*gqxz(30)+qyzk*gqyz(30)))
     &        + 2.0d0*(qxyi*(qxxk*gqxx(26)+qyyk*gqyy(26)+qzzk*gqzz(26)
     &              +2.0d0*(qxyk*gqxy(26)+qxzk*gqxz(26)+qyzk*gqyz(26)))
     &                 + qxzi*(qxxk*gqxx(27)+qyyk*gqyy(27)+qzzk*gqzz(27)
     &              +2.0d0*(qxyk*gqxy(27)+qxzk*gqxz(27)+qyzk*gqyz(27)))
     &                 + qyzi*(qxxk*gqxx(29)+qyyk*gqyy(29)+qzzk*gqzz(29)
     &              +2.0d0*(qxyk*gqxy(29)+qxzk*gqxz(29)+qyzk*gqyz(29))))
                  dsumdr = desymdr + 0.5d0*(dewidr+dewkdr)
                  drbi = rbk*dsumdr
                  drbk = rbi*dsumdr
c
c     torque on permanent dipoles due to permanent reaction field
c
                  if (i .ne. k) then
                     fid(1) = uxk*gux(2) + uyk*gux(3) + uzk*gux(4)
     &          + 0.5d0*(ck*gux(1)+qxxk*gux(5)+qyyk*gux(8)+qzzk*gux(10)
     &                +2.0d0*(qxyk*gux(6)+qxzk*gux(7)+qyzk*gux(9))
     &                +ck*gc(2)+qxxk*gqxx(2)+qyyk*gqyy(2)+qzzk*gqzz(2)
     &                +2.0d0*(qxyk*gqxy(2)+qxzk*gqxz(2)+qyzk*gqyz(2)))
                     fid(2) = uxk*guy(2) + uyk*guy(3) + uzk*guy(4)
     &          + 0.5d0*(ck*guy(1)+qxxk*guy(5)+qyyk*guy(8)+qzzk*guy(10)
     &                +2.0d0*(qxyk*guy(6)+qxzk*guy(7)+qyzk*guy(9))
     &                +ck*gc(3)+qxxk*gqxx(3)+qyyk*gqyy(3)+qzzk*gqzz(3)
     &                +2.0d0*(qxyk*gqxy(3)+qxzk*gqxz(3)+qyzk*gqyz(3)))
                     fid(3) = uxk*guz(2) + uyk*guz(3) + uzk*guz(4)
     &          + 0.5d0*(ck*guz(1)+qxxk*guz(5)+qyyk*guz(8)+qzzk*guz(10)
     &                +2.0d0*(qxyk*guz(6)+qxzk*guz(7)+qyzk*guz(9))
     &                +ck*gc(4)+qxxk*gqxx(4)+qyyk*gqyy(4)+qzzk*gqzz(4)
     &                +2.0d0*(qxyk*gqxy(4)+qxzk*gqxz(4)+qyzk*gqyz(4)))
                     fkd(1) = uxi*gux(2) + uyi*gux(3) + uzi*gux(4)
     &          - 0.5d0*(ci*gux(1)+qxxi*gux(5)+qyyi*gux(8)+qzzi*gux(10)
     &                +2.0d0*(qxyi*gux(6)+qxzi*gux(7)+qyzi*gux(9))
     &                +ci*gc(2)+qxxi*gqxx(2)+qyyi*gqyy(2)+qzzi*gqzz(2)
     &                +2.0d0*(qxyi*gqxy(2)+qxzi*gqxz(2)+qyzi*gqyz(2)))
                     fkd(2) = uxi*guy(2) + uyi*guy(3) + uzi*guy(4)
     &          - 0.5d0*(ci*guy(1)+qxxi*guy(5)+qyyi*guy(8)+qzzi*guy(10)
     &                +2.0d0*(qxyi*guy(6)+qxzi*guy(7)+qyzi*guy(9))
     &                +ci*gc(3)+qxxi*gqxx(3)+qyyi*gqyy(3)+qzzi*gqzz(3)
     &                +2.0d0*(qxyi*gqxy(3)+qxzi*gqxz(3)+qyzi*gqyz(3)))
                     fkd(3) = uxi*guz(2) + uyi*guz(3) + uzi*guz(4)
     &          - 0.5d0*(ci*guz(1)+qxxi*guz(5)+qyyi*guz(8)+qzzi*guz(10)
     &                +2.0d0*(qxyi*guz(6)+qxzi*guz(7)+qyzi*guz(9))
     &                +ci*gc(4)+qxxi*gqxx(4)+qyyi*gqyy(4)+qzzi*gqzz(4)
     &                +2.0d0*(qxyi*gqxy(4)+qxzi*gqxz(4)+qyzi*gqyz(4)))
                     trq(1,ii) = trq(1,ii) + uyi*fid(3) - uzi*fid(2)
                     trq(2,ii) = trq(2,ii) + uzi*fid(1) - uxi*fid(3)
                     trq(3,ii) = trq(3,ii) + uxi*fid(2) - uyi*fid(1)
                     trq(1,kk) = trq(1,kk) + uyk*fkd(3) - uzk*fkd(2)
                     trq(2,kk) = trq(2,kk) + uzk*fkd(1) - uxk*fkd(3)
                     trq(3,kk) = trq(3,kk) + uxk*fkd(2) - uyk*fkd(1)
c
c     torque on quadrupoles due to permanent reaction field gradient
c
                     fidg(1,1) =
     &          - 0.5d0*(ck*gqxx(1)+uxk*gqxx(2)+uyk*gqxx(3)+uzk*gqxx(4)
     &                +qxxk*gqxx(5)+qyyk*gqxx(8)+qzzk*gqxx(10)
     &                +2.0d0*(qxyk*gqxx(6)+qxzk*gqxx(7)+qyzk*gqxx(9))
     &                +ck*gc(5)+uxk*gux(5)+uyk*guy(5)+uzk*guz(5)
     &                +qxxk*gqxx(5)+qyyk*gqyy(5)+qzzk*gqzz(5)
     &                +2.0d0*(qxyk*gqxy(5)+qxzk*gqxz(5)+qyzk*gqyz(5)))
                     fidg(1,2) =
     &          - 0.5d0*(ck*gqxy(1)+uxk*gqxy(2)+uyk*gqxy(3)+uzk*gqxy(4)
     &                +qxxk*gqxy(5)+qyyk*gqxy(8)+qzzk*gqxy(10)
     &                +2.0d0*(qxyk*gqxy(6)+qxzk*gqxy(7)+qyzk*gqxy(9))
     &                +ck*gc(6)+uxk*gux(6)+uyk*guy(6)+uzk*guz(6)
     &                +qxxk*gqxx(6)+qyyk*gqyy(6)+qzzk*gqzz(6)
     &                +2.0d0*(qxyk*gqxy(6)+qxzk*gqxz(6)+qyzk*gqyz(6)))
                     fidg(1,3) =
     &          - 0.5d0*(ck*gqxz(1)+uxk*gqxz(2)+uyk*gqxz(3)+uzk*gqxz(4)
     &                +qxxk*gqxz(5)+qyyk*gqxz(8)+qzzk*gqxz(10)
     &                +2.0d0*(qxyk*gqxz(6)+qxzk*gqxz(7)+qyzk*gqxz(9))
     &                +ck*gc(7)+uxk*gux(7)+uyk*guy(7)+uzk*guz(7)
     &                +qxxk*gqxx(7)+qyyk*gqyy(7)+qzzk*gqzz(7)
     &                +2.0d0*(qxyk*gqxy(7)+qxzk*gqxz(7)+qyzk*gqyz(7)))
                     fidg(2,2) =
     &          - 0.5d0*(ck*gqyy(1)+uxk*gqyy(2)+uyk*gqyy(3)+uzk*gqyy(4)
     &                +qxxk*gqyy(5)+qyyk*gqyy(8)+qzzk*gqyy(10)
     &                +2.0d0*(qxyk*gqyy(6)+qxzk*gqyy(7)+qyzk*gqyy(9))
     &                +ck*gc(8)+uxk*gux(8)+uyk*guy(8)+uzk*guz(8)
     &                +qxxk*gqxx(8)+qyyk*gqyy(8)+qzzk*gqzz(8)
     &                +2.0d0*(qxyk*gqxy(8)+qxzk*gqxz(8)+qyzk*gqyz(8)))
                     fidg(2,3) =
     &          - 0.5d0*(ck*gqyz(1)+uxk*gqyz(2)+uyk*gqyz(3)+uzk*gqyz(4)
     &                +qxxk*gqyz(5)+qyyk*gqyz(8)+qzzk*gqyz(10)
     &                +2.0d0*(qxyk*gqyz(6)+qxzk*gqyz(7)+qyzk*gqyz(9))
     &                +ck*gc(9)+uxk*gux(9)+uyk*guy(9)+uzk*guz(9)
     &                +qxxk*gqxx(9)+qyyk*gqyy(9)+qzzk*gqzz(9)
     &                +2.0d0*(qxyk*gqxy(9)+qxzk*gqxz(9)+qyzk*gqyz(9)))
                     fidg(3,3) =
     &          - 0.5d0*(ck*gqzz(1)+uxk*gqzz(2)+uyk*gqzz(3)+uzk*gqzz(4)
     &                +qxxk*gqzz(5)+qyyk*gqzz(8)+qzzk*gqzz(10)
     &                +2.0d0*(qxyk*gqzz(6)+qxzk*gqzz(7)+qyzk*gqzz(9))
     &                +ck*gc(10)+uxk*gux(10)+uyk*guy(10)+uzk*guz(10)
     &                +qxxk*gqxx(10)+qyyk*gqyy(10)+qzzk*gqzz(10)
     &             +2.0d0*(qxyk*gqxy(10)+qxzk*gqxz(10)+qyzk*gqyz(10)))
                     fidg(2,1) = fidg(1,2)
                     fidg(3,1) = fidg(1,3)
                     fidg(3,2) = fidg(2,3)
                     fkdg(1,1) =
     &          - 0.5d0*(ci*gqxx(1)-uxi*gqxx(2)-uyi*gqxx(3)-uzi *gqxx(4)
     &                +qxxi*gqxx(5)+qyyi*gqxx(8)+qzzi*gqxx(10)
     &                +2.0d0*(qxyi*gqxx(6)+qxzi*gqxx(7)+qyzi*gqxx(9))
     &                +ci*gc(5)-uxi*gux(5)-uyi*guy(5)-uzi*guz(5)
     &                +qxxi*gqxx(5)+qyyi*gqyy(5)+qzzi*gqzz(5)
     &                +2.0d0*(qxyi*gqxy(5)+qxzi*gqxz(5)+qyzi*gqyz(5)))
                     fkdg(1,2) =
     &          - 0.5d0*(ci*gqxy(1)-uxi*gqxy(2)-uyi*gqxy(3)-uzi*gqxy(4)
     &                +qxxi*gqxy(5)+qyyi*gqxy(8)+qzzi*gqxy(10)
     &                +2.0d0*(qxyi*gqxy(6)+qxzi*gqxy(7)+qyzi*gqxy(9))
     &                +ci*gc(6)-uxi*gux(6)-uyi*guy(6)-uzi*guz(6)
     &                +qxxi*gqxx(6)+qyyi*gqyy(6)+qzzi*gqzz(6)
     &                +2.0d0*(qxyi*gqxy(6)+qxzi*gqxz(6)+qyzi*gqyz(6)))
                     fkdg(1,3) =
     &          - 0.5d0*(ci*gqxz(1)-uxi*gqxz(2)-uyi*gqxz(3)-uzi*gqxz(4)
     &                +qxxi*gqxz(5)+qyyi*gqxz(8)+qzzi*gqxz(10)
     &                +2.0d0*(qxyi*gqxz(6)+qxzi*gqxz(7)+qyzi*gqxz(9))
     &                +ci*gc(7)-uxi*gux(7)-uyi*guy(7)-uzi*guz(7)
     &                +qxxi*gqxx(7)+qyyi*gqyy(7)+qzzi*gqzz(7)
     &                +2.0d0*(qxyi*gqxy(7)+qxzi*gqxz(7)+qyzi*gqyz(7)))
                     fkdg(2,2) =
     &          - 0.5d0*(ci*gqyy(1)-uxi*gqyy(2)-uyi*gqyy(3)-uzi*gqyy(4)
     &                +qxxi*gqyy(5)+qyyi*gqyy(8)+qzzi*gqyy(10)
     &                +2.0d0*(qxyi*gqyy(6)+qxzi*gqyy(7)+qyzi*gqyy(9))
     &                +ci*gc(8)-uxi*gux(8)-uyi*guy(8)-uzi*guz(8)
     &                +qxxi*gqxx(8)+qyyi*gqyy(8)+qzzi*gqzz(8)
     &                +2.0d0*(qxyi*gqxy(8)+qxzi*gqxz(8)+qyzi*gqyz(8)))
                     fkdg(2,3) =
     &          - 0.5d0*(ci*gqyz(1)-uxi*gqyz(2)-uyi*gqyz(3)-uzi*gqyz(4)
     &                +qxxi*gqyz(5)+qyyi*gqyz(8)+qzzi*gqyz(10)
     &                +2.0d0*(qxyi*gqyz(6)+qxzi*gqyz(7)+qyzi*gqyz(9))
     &                +ci*gc(9)-uxi*gux(9)-uyi*guy(9)-uzi*guz(9)
     &                +qxxi*gqxx(9)+qyyi*gqyy(9)+qzzi*gqzz(9)
     &                +2.0d0*(qxyi*gqxy(9)+qxzi*gqxz(9)+qyzi*gqyz(9)))
                     fkdg(3,3) =
     &          - 0.5d0*(ci*gqzz(1)-uxi*gqzz(2)-uyi*gqzz(3)-uzi*gqzz(4)
     &                +qxxi*gqzz(5)+qyyi*gqzz(8)+qzzi*gqzz(10)
     &                +2.0d0*(qxyi*gqzz(6)+qxzi*gqzz(7)+qyzi*gqzz(9))
     &                +ci*gc(10)-uxi*gux(10)-uyi*guy(10)-uzi*guz(10)
     &                +qxxi*gqxx(10)+qyyi*gqyy(10)+qzzi*gqzz(10)
     &              +2.0d0*(qxyi*gqxy(10)+qxzi*gqxz(10)+qyzi*gqyz(10)))
                     fkdg(2,1) = fkdg(1,2)
                     fkdg(3,1) = fkdg(1,3)
                     fkdg(3,2) = fkdg(2,3)
                     trq(1,ii) = trq(1,ii) + 2.0d0*
     &                    (qxyi*fidg(1,3)+qyyi*fidg(2,3)+qyzi*fidg(3,3)
     &                    -qxzi*fidg(1,2)-qyzi*fidg(2,2)-qzzi*fidg(3,2))
                     trq(2,ii) = trq(2,ii) + 2.0d0*
     &                    (qxzi*fidg(1,1)+qyzi*fidg(2,1)+qzzi*fidg(3,1)
     &                    -qxxi*fidg(1,3)-qxyi*fidg(2,3)-qxzi*fidg(3,3))
                     trq(3,ii) = trq(3,ii) + 2.0d0*
     &                    (qxxi*fidg(1,2)+qxyi*fidg(2,2)+qxzi*fidg(3,2)
     &                    -qxyi*fidg(1,1)-qyyi*fidg(2,1)-qyzi*fidg(3,1))
                     trq(1,kk) = trq(1,kk) + 2.0d0*
     &                    (qxyk*fkdg(1,3)+qyyk*fkdg(2,3)+qyzk*fkdg(3,3)
     &                    -qxzk*fkdg(1,2)-qyzk*fkdg(2,2)-qzzk*fkdg(3,2))
                     trq(2,kk) = trq(2,kk) + 2.0d0*
     &                    (qxzk*fkdg(1,1)+qyzk*fkdg(2,1)+qzzk*fkdg(3,1)
     &                    -qxxk*fkdg(1,3)-qxyk*fkdg(2,3)-qxzk*fkdg(3,3))
                     trq(3,kk) = trq(3,kk) + 2.0d0*
     &                    (qxxk*fkdg(1,2)+qxyk*fkdg(2,2)+qxzk*fkdg(3,2)
     &                    -qxyk*fkdg(1,1)-qyyk*fkdg(2,1)-qyzk*fkdg(3,1))
                  end if
c
c     electrostatic solvation energy of the permanent multipoles in
c     the GK reaction potential of the induced dipoles
c
                  esymi = -uxi*(dxk*gux(2)+dyk*guy(2)+dzk*guz(2))
     &                   - uyi*(dxk*gux(3)+dyk*guy(3)+dzk*guz(3))
     &                   - uzi*(dxk*gux(4)+dyk*guy(4)+dzk*guz(4))
     &                   - uxk*(dxi*gux(2)+dyi*guy(2)+dzi*guz(2))
     &                   - uyk*(dxi*gux(3)+dyi*guy(3)+dzi*guz(3))
     &                   - uzk*(dxi*gux(4)+dyi*guy(4)+dzi*guz(4))
                  ewii = ci*(dxk*gc(2)+dyk*gc(3)+dzk*gc(4))
     &                 - ck*(dxi*gux(1)+dyi*guy(1)+dzi*guz(1))
     &                 - dxi*(qxxk*gux(5)+qyyk*gux(8)+qzzk*gux(10)
     &                +2.0d0*(qxyk*gux(6)+qxzk*gux(7)+qyzk*gux(9)))
     &                 - dyi*(qxxk*guy(5)+qyyk*guy(8)+qzzk*guy(10)
     &                +2.0d0*(qxyk*guy(6)+qxzk*guy(7)+qyzk*guy(9)))
     &                 - dzi*(qxxk*guz(5)+qyyk*guz(8)+qzzk*guz(10)
     &                +2.0d0*(qxyk*guz(6)+qxzk*guz(7)+qyzk*guz(9)))
     &                 + dxk*(qxxi*gqxx(2)+qyyi*gqyy(2)+qzzi*gqzz(2)
     &                +2.0d0*(qxyi*gqxy(2)+qxzi*gqxz(2)+qyzi*gqyz(2)))
     &                 + dyk*(qxxi*gqxx(3)+qyyi*gqyy(3)+qzzi*gqzz(3)
     &                +2.0d0*(qxyi*gqxy(3)+qxzi*gqxz(3)+qyzi*gqyz(3)))
     &                 + dzk*(qxxi*gqxx(4)+qyyi*gqyy(4)+qzzi*gqzz(4)
     &                +2.0d0*(qxyi*gqxy(4)+qxzi*gqxz(4)+qyzi*gqyz(4)))
                  ewki = ci*(dxk*gux(1)+dyk*guy(1)+dzk*guz(1))
     &                 - ck*(dxi*gc(2)+dyi*gc(3)+dzi*gc(4))
     &                 - dxi*(qxxk*gqxx(2)+qyyk*gqyy(2)+qzzk*gqzz(2)
     &                +2.0d0*(qxyk*gqxy(2)+qxzk*gqxz(2)+qyzk*gqyz(2)))
     &                 - dyi*(qxxk*gqxx(3)+qyyk*gqyy(3)+qzzk*gqzz(3)
     &                +2.0d0*(qxyk*gqxy(3)+qxzk*gqxz(3)+qyzk*gqyz(3)))
     &                 - dzi*(qxxk*gqxx(4)+qyyk*gqyy(4)+qzzk*gqzz(4)
     &                +2.0d0*(qxyk*gqxy(4)+qxzk*gqxz(4)+qyzk*gqyz(4)))
     &                 + dxk*(qxxi*gux(5)+qyyi*gux(8)+qzzi*gux(10)
     &                +2.0d0*(qxyi*gux(6)+qxzi*gux(7)+qyzi*gux(9)))
     &                 + dyk*(qxxi*guy(5)+qyyi*guy(8)+qzzi*guy(10)
     &                +2.0d0*(qxyi*guy(6)+qxzi*guy(7)+qyzi*guy(9)))
     &                 + dzk*(qxxi*guz(5)+qyyi*guz(8)+qzzi*guz(10)
     &                +2.0d0*(qxyi*guz(6)+qxzi*guz(7)+qyzi*guz(9)))
c
c     electrostatic solvation free energy gradient of the permanent
c     multipoles in the reaction potential of the induced dipoles
c
                  dpsymdx = -uxi*(sxk*gux(5)+syk*guy(5)+szk*guz(5))
     &                     - uyi*(sxk*gux(6)+syk*guy(6)+szk*guz(6))
     &                     - uzi*(sxk*gux(7)+syk*guy(7)+szk*guz(7))
     &                     - uxk*(sxi*gux(5)+syi*guy(5)+szi*guz(5))
     &                     - uyk*(sxi*gux(6)+syi*guy(6)+szi*guz(6))
     &                     - uzk*(sxi*gux(7)+syi*guy(7)+szi*guz(7))
                  dpwidx = ci*(sxk*gc(5)+syk*gc(6)+szk*gc(7))
     &                   - ck*(sxi*gux(2)+syi*guy(2)+szi*guz(2))
     &                   - sxi*(qxxk*gux(11)+qyyk*gux(14)+qzzk*gux(16)
     &                +2.0d0*(qxyk*gux(12)+qxzk*gux(13)+qyzk*gux(15)))
     &                   - syi*(qxxk*guy(11)+qyyk*guy(14)+qzzk*guy(16)
     &                +2.0d0*(qxyk*guy(12)+qxzk*guy(13)+qyzk*guy(15)))
     &                   - szi*(qxxk*guz(11)+qyyk*guz(14)+qzzk*guz(16)
     &                +2.0d0*(qxyk*guz(12)+qxzk*guz(13)+qyzk*guz(15)))
     &                   + sxk*(qxxi*gqxx(5)+qyyi*gqyy(5)+qzzi*gqzz(5)
     &                +2.0d0*(qxyi*gqxy(5)+qxzi*gqxz(5)+qyzi*gqyz(5)))
     &                   + syk*(qxxi*gqxx(6)+qyyi*gqyy(6)+qzzi*gqzz(6)
     &                +2.0d0*(qxyi*gqxy(6)+qxzi*gqxz(6)+qyzi*gqyz(6)))
     &                   + szk*(qxxi*gqxx(7)+qyyi*gqyy(7)+qzzi*gqzz(7)
     &                +2.0d0*(qxyi*gqxy(7)+qxzi*gqxz(7)+qyzi*gqyz(7)))
                  dpwkdx = ci*(sxk*gux(2)+syk*guy(2)+szk*guz(2))
     &                   - ck*(sxi*gc(5)+syi*gc(6)+szi*gc(7))
     &                   - sxi*(qxxk*gqxx(5)+qyyk*gqyy(5)+qzzk*gqzz(5)
     &                +2.0d0*(qxyk*gqxy(5)+qxzk*gqxz(5)+qyzk*gqyz(5)))
     &                   - syi*(qxxk*gqxx(6)+qyyk*gqyy(6)+qzzk*gqzz(6)
     &                +2.0d0*(qxyk*gqxy(6)+qxzk*gqxz(6)+qyzk*gqyz(6)))
     &                   - szi*(qxxk*gqxx(7)+qyyk*gqyy(7)+qzzk*gqzz(7)
     &                +2.0d0*(qxyk*gqxy(7)+qxzk*gqxz(7)+qyzk*gqyz(7)))
     &                   + sxk*(qxxi*gux(11)+qyyi*gux(14)+qzzi*gux(16)
     &                +2.0d0*(qxyi*gux(12)+qxzi*gux(13)+qyzi*gux(15)))
     &                   + syk*(qxxi*guy(11)+qyyi*guy(14)+qzzi*guy(16)
     &                +2.0d0*(qxyi*guy(12)+qxzi*guy(13)+qyzi*guy(15)))
     &                   + szk*(qxxi*guz(11)+qyyi*guz(14)+qzzi*guz(16)
     &                +2.0d0*(qxyi*guz(12)+qxzi*guz(13)+qyzi*guz(15)))
                  dpdx = 0.5d0 * (dpsymdx + 0.5d0*(dpwidx + dpwkdx))
                  dpsymdy = -uxi*(sxk*gux(6)+syk*guy(6)+szk*guz(6))
     &                     - uyi*(sxk*gux(8)+syk*guy(8)+szk*guz(8))
     &                     - uzi*(sxk*gux(9)+syk*guy(9)+szk*guz(9))
     &                     - uxk*(sxi*gux(6)+syi*guy(6)+szi*guz(6))
     &                     - uyk*(sxi*gux(8)+syi*guy(8)+szi*guz(8))
     &                     - uzk*(sxi*gux(9)+syi*guy(9)+szi*guz(9))
                  dpwidy = ci*(sxk*gc(6)+syk*gc(8)+szk*gc(9))
     &                   - ck*(sxi*gux(3)+syi*guy(3)+szi*guz(3))
     &                   - sxi*(qxxk*gux(12)+qyyk*gux(17)+qzzk*gux(19)
     &                +2.0d0*(qxyk*gux(14)+qxzk*gux(15)+qyzk*gux(18)))
     &                   - syi*(qxxk*guy(12)+qyyk*guy(17)+qzzk*guy(19)
     &                +2.0d0*(qxyk*guy(14)+qxzk*guy(15)+qyzk*guy(18)))
     &                   - szi*(qxxk*guz(12)+qyyk*guz(17)+qzzk*guz(19)
     &                +2.0d0*(qxyk*guz(14)+qxzk*guz(15)+qyzk*guz(18)))
     &                   + sxk*(qxxi*gqxx(6)+qyyi*gqyy(6)+qzzi*gqzz(6)
     &                +2.0d0*(qxyi*gqxy(6)+qxzi*gqxz(6)+qyzi*gqyz(6)))
     &                   + syk*(qxxi*gqxx(8)+qyyi*gqyy(8)+qzzi*gqzz(8)
     &                +2.0d0*(qxyi*gqxy(8)+qxzi*gqxz(8)+qyzi*gqyz(8)))
     &                   + szk*(qxxi*gqxx(9)+qyyi*gqyy(9)+qzzi*gqzz(9)
     &                +2.0d0*(qxyi*gqxy(9)+qxzi*gqxz(9)+qyzi*gqyz(9)))
                  dpwkdy = ci*(sxk*gux(3)+syk*guy(3)+szk*guz(3))
     &                   - ck*(sxi*gc(6)+syi*gc(8)+szi*gc(9))
     &                   - sxi*(qxxk*gqxx(6)+qyyk*gqyy(6)+qzzk*gqzz(6)
     &                +2.0d0*(qxyk*gqxy(6)+qxzk*gqxz(6)+qyzk*gqyz(6)))
     &                   - syi*(qxxk*gqxx(8)+qyyk*gqyy(8)+qzzk*gqzz(8)
     &                +2.0d0*(qxyk*gqxy(8)+qxzk*gqxz(8)+qyzk*gqyz(8)))
     &                   - szi*(qxxk*gqxx(9)+qyyk*gqyy(9)+qzzk*gqzz(9)
     &                +2.0d0*(qxyk*gqxy(9)+qxzk*gqxz(9)+qyzk*gqyz(9)))
     &                   + sxk*(qxxi*gux(12)+qyyi*gux(17)+qzzi*gux(19)
     &                +2.0d0*(qxyi*gux(14)+qxzi*gux(15)+qyzi*gux(18)))
     &                   + syk*(qxxi*guy(12)+qyyi*guy(17)+qzzi*guy(19)
     &                +2.0d0*(qxyi*guy(14)+qxzi*guy(15)+qyzi*guy(18)))
     &                   + szk*(qxxi*guz(12)+qyyi*guz(17)+qzzi*guz(19)
     &                +2.0d0*(qxyi*guz(14)+qxzi*guz(15)+qyzi*guz(18)))
                  dpdy = 0.5d0 * (dpsymdy + 0.5d0*(dpwidy + dpwkdy))
                  dpsymdz = -uxi*(sxk*gux(7)+syk*guy(7)+szk*guz(7))
     &                     - uyi*(sxk*gux(9)+syk*guy(9)+szk*guz(9))
     &                     - uzi*(sxk*gux(10)+syk*guy(10)+szk*guz(10))
     &                     - uxk*(sxi*gux(7)+syi*guy(7)+szi*guz(7))
     &                     - uyk*(sxi*gux(9)+syi*guy(9)+szi*guz(9))
     &                     - uzk*(sxi*gux(10)+syi*guy(10)+szi*guz(10))
                  dpwidz = ci*(sxk*gc(7)+syk*gc(9)+szk*gc(10))
     &                   - ck*(sxi*gux(4)+syi*guy(4)+szi*guz(4))
     &                   - sxi*(qxxk*gux(13)+qyyk*gux(18)+qzzk*gux(20)
     &                +2.0d0*(qxyk*gux(15)+qxzk*gux(16)+qyzk*gux(19)))
     &                   - syi*(qxxk*guy(13)+qyyk*guy(18)+qzzk*guy(20)
     &                +2.0d0*(qxyk*guy(15)+qxzk*guy(16)+qyzk*guy(19)))
     &                   - szi*(qxxk*guz(13)+qyyk*guz(18)+qzzk*guz(20)
     &                +2.0d0*(qxyk*guz(15)+qxzk*guz(16)+qyzk*guz(19)))
     &                   + sxk*(qxxi*gqxx(7)+qyyi*gqyy(7)+qzzi*gqzz(7)
     &                +2.0d0*(qxyi*gqxy(7)+qxzi*gqxz(7)+qyzi*gqyz(7)))
     &                   + syk*(qxxi*gqxx(9)+qyyi*gqyy(9)+qzzi*gqzz(9)
     &                +2.0d0*(qxyi*gqxy(9)+qxzi*gqxz(9)+qyzi*gqyz(9)))
     &                  + szk*(qxxi*gqxx(10)+qyyi*gqyy(10)+qzzi*gqzz(10)
     &               +2.0d0*(qxyi*gqxy(10)+qxzi*gqxz(10)+qyzi*gqyz(10)))
                  dpwkdz = ci*(sxk*gux(4)+syk*guy(4)+szk*guz(4))
     &                   - ck*(sxi*gc(7)+syi*gc(9)+szi*gc(10))
     &                   - sxi*(qxxk*gqxx(7)+qyyk*gqyy(7)+qzzk*gqzz(7)
     &                +2.0d0*(qxyk*gqxy(7)+qxzk*gqxz(7)+qyzk*gqyz(7)))
     &                   - syi*(qxxk*gqxx(9)+qyyk*gqyy(9)+qzzk*gqzz(9)
     &                +2.0d0*(qxyk*gqxy(9)+qxzk*gqxz(9)+qyzk*gqyz(9)))
     &                  - szi*(qxxk*gqxx(10)+qyyk*gqyy(10)+qzzk*gqzz(10)
     &               +2.0d0*(qxyk*gqxy(10)+qxzk*gqxz(10)+qyzk*gqyz(10)))
     &                   + sxk*(qxxi*gux(13)+qyyi*gux(18)+qzzi*gux(20)
     &                +2.0d0*(qxyi*gux(15)+qxzi*gux(16)+qyzi*gux(19)))
     &                   + syk*(qxxi*guy(13)+qyyi*guy(18)+qzzi*guy(20)
     &                +2.0d0*(qxyi*guy(15)+qxzi*guy(16)+qyzi*guy(19)))
     &                   + szk*(qxxi*guz(13)+qyyi*guz(18)+qzzi*guz(20)
     &                +2.0d0*(qxyi*guz(15)+qxzi*guz(16)+qyzi*guz(19)))
                  dpdz = 0.5d0 * (dpsymdz + 0.5d0*(dpwidz+dpwkdz))
c
c     effective radii chain rule terms for the electrostatic solvation
c     free energy gradient of the permanent multipoles in the reaction
c     potential of the induced dipoles
c
                  dsymdr = -uxi*(sxk*gux(22)+syk*guy(22)+szk*guz(22))
     &                    - uyi*(sxk*gux(23)+syk*guy(23)+szk*guz(23))
     &                    - uzi*(sxk*gux(24)+syk*guy(24)+szk*guz(24))
     &                    - uxk*(sxi*gux(22)+syi*guy(22)+szi*guz(22))
     &                    - uyk*(sxi*gux(23)+syi*guy(23)+szi*guz(23))
     &                    - uzk*(sxi*gux(24)+syi*guy(24)+szi*guz(24))
                  dwipdr = ci*(sxk*gc(22)+syk*gc(23)+szk*gc(24))
     &                   - ck*(sxi*gux(21)+syi*guy(21)+szi*guz(21))
     &                - sxi*(qxxk*gux(25)+qyyk*gux(28)+qzzk*gux(30)
     &               +2.0d0*(qxyk*gux(26)+qxzk*gux(27)+qyzk*gux(29)))
     &                - syi*(qxxk*guy(25)+qyyk*guy(28)+qzzk*guy(30)
     &               +2.0d0*(qxyk*guy(26)+qxzk*guy(27)+qyzk*guy(29)))
     &                - szi*(qxxk*guz(25)+qyyk*guz(28)+qzzk*guz(30)
     &               +2.0d0*(qxyk*guz(26)+qxzk*guz(27)+qyzk*guz(29)))
     &                + sxk*(qxxi*gqxx(22)+qyyi*gqyy(22)+qzzi*gqzz(22)
     &               +2.0d0*(qxyi*gqxy(22)+qxzi*gqxz(22)+qyzi*gqyz(22)))
     &                + syk*(qxxi*gqxx(23)+qyyi*gqyy(23)+qzzi*gqzz(23)
     &               +2.0d0*(qxyi*gqxy(23)+qxzi*gqxz(23)+qyzi*gqyz(23)))
     &                + szk*(qxxi*gqxx(24)+qyyi*gqyy(24)+qzzi*gqzz(24)
     &               +2.0d0*(qxyi*gqxy(24)+qxzi*gqxz(24)+qyzi*gqyz(24)))
                  dwkpdr = ci*(sxk*gux(21)+syk*guy(21)+szk*guz(21))
     &                   - ck*(sxi*gc(22)+syi*gc(23)+szi*gc(24))
     &                - sxi*(qxxk*gqxx(22)+qyyk*gqyy(22)+qzzk*gqzz(22)
     &               +2.0d0*(qxyk*gqxy(22)+qxzk*gqxz(22)+qyzk*gqyz(22)))
     &                - syi*(qxxk*gqxx(23)+qyyk*gqyy(23)+qzzk*gqzz(23)
     &               +2.0d0*(qxyk*gqxy(23)+qxzk*gqxz(23)+qyzk*gqyz(23)))
     &                - szi*(qxxk*gqxx(24)+qyyk*gqyy(24)+qzzk*gqzz(24)
     &               +2.0d0*(qxyk*gqxy(24)+qxzk*gqxz(24)+qyzk*gqyz(24)))
     &                + sxk*(qxxi*gux(25)+qyyi*gux(28)+qzzi*gux(30)
     &               +2.0d0*(qxyi*gux(26)+qxzi*gux(27)+qyzi*gux(29)))
     &                + syk*(qxxi*guy(25)+qyyi*guy(28)+qzzi*guy(30)
     &               +2.0d0*(qxyi*guy(26)+qxzi*guy(27)+qyzi*guy(29)))
     &                + szk*(qxxi*guz(25)+qyyi*guz(28)+qzzi*guz(30)
     &               +2.0d0*(qxyi*guz(26)+qxzi*guz(27)+qyzi*guz(29)))
                  dsumdr = dsymdr + 0.5d0*(dwipdr+dwkpdr)
                  dpbi = 0.5d0*rbk*dsumdr
                  dpbk = 0.5d0*rbi*dsumdr
c
c     mutual polarization electrostatic solvation free energy gradient
c
                  if (poltyp .eq. 'MUTUAL') then
                     dpdx = dpdx - 0.5d0 *
     &                      (dxi*(pxk*gux(5)+pyk*gux(6)+pzk*gux(7))
     &                      + dyi*(pxk*guy(5)+pyk*guy(6)+pzk*guy(7))
     &                      + dzi*(pxk*guz(5)+pyk*guz(6)+pzk*guz(7))
     &                      + dxk*(pxi*gux(5)+pyi*gux(6)+pzi*gux(7))
     &                      + dyk*(pxi*guy(5)+pyi*guy(6)+pzi*guy(7))
     &                      + dzk*(pxi*guz(5)+pyi*guz(6)+pzi*guz(7)))
                     dpdy = dpdy - 0.5d0 *
     &                      (dxi*(pxk*gux(6)+pyk*gux(8)+pzk*gux(9))
     &                      + dyi*(pxk*guy(6)+pyk*guy(8)+pzk*guy(9))
     &                      + dzi*(pxk*guz(6)+pyk*guz(8)+pzk*guz(9))
     &                      + dxk*(pxi*gux(6)+pyi*gux(8)+pzi*gux(9))
     &                      + dyk*(pxi*guy(6)+pyi*guy(8)+pzi*guy(9))
     &                      + dzk*(pxi*guz(6)+pyi*guz(8)+pzi*guz(9)))
                     dpdz = dpdz - 0.5d0 *
     &                      (dxi*(pxk*gux(7)+pyk*gux(9)+pzk*gux(10))
     &                      + dyi*(pxk*guy(7)+pyk*guy(9)+pzk*guy(10))
     &                      + dzi*(pxk*guz(7)+pyk*guz(9)+pzk*guz(10))
     &                      + dxk*(pxi*gux(7)+pyi*gux(9)+pzi*gux(10))
     &                      + dyk*(pxi*guy(7)+pyi*guy(9)+pzi*guy(10))
     &                      + dzk*(pxi*guz(7)+pyi*guz(9)+pzi*guz(10)))
                     duvdr = dxi*(pxk*gux(22)+pyk*gux(23)+pzk*gux(24))
     &                      + dyi*(pxk*guy(22)+pyk*guy(23)+pzk*guy(24))
     &                      + dzi*(pxk*guz(22)+pyk*guz(23)+pzk*guz(24))
     &                      + dxk*(pxi*gux(22)+pyi*gux(23)+pzi*gux(24))
     &                      + dyk*(pxi*guy(22)+pyi*guy(23)+pzi*guy(24))
     &                      + dzk*(pxi*guz(22)+pyi*guz(23)+pzi*guz(24))
                     dpbi = dpbi - 0.5d0*rbk*duvdr
                     dpbk = dpbk - 0.5d0*rbi*duvdr
                  end if
c
c     torque due to induced reaction field on permanent dipoles
c
                  fid(1) = 0.5d0 * (sxk*gux(2)+syk*guy(2)+szk*guz(2))
                  fid(2) = 0.5d0 * (sxk*gux(3)+syk*guy(3)+szk*guz(3))
                  fid(3) = 0.5d0 * (sxk*gux(4)+syk*guy(4)+szk*guz(4))
                  fkd(1) = 0.5d0 * (sxi*gux(2)+syi*guy(2)+szi*guz(2))
                  fkd(2) = 0.5d0 * (sxi*gux(3)+syi*guy(3)+szi*guz(3))
                  fkd(3) = 0.5d0 * (sxi*gux(4)+syi*guy(4)+szi*guz(4))
                  if (i .eq. k) then
                     fid(1) = 0.5d0 * fid(1)
                     fid(2) = 0.5d0 * fid(2)
                     fid(3) = 0.5d0 * fid(3)
                     fkd(1) = 0.5d0 * fkd(1)
                     fkd(2) = 0.5d0 * fkd(2)
                     fkd(3) = 0.5d0 * fkd(3)
                  end if
                  trqi(1,ii) = trqi(1,ii) + uyi*fid(3) - uzi*fid(2)
                  trqi(2,ii) = trqi(2,ii) + uzi*fid(1) - uxi*fid(3)
                  trqi(3,ii) = trqi(3,ii) + uxi*fid(2) - uyi*fid(1)
                  trqi(1,kk) = trqi(1,kk) + uyk*fkd(3) - uzk*fkd(2)
                  trqi(2,kk) = trqi(2,kk) + uzk*fkd(1) - uxk*fkd(3)
                  trqi(3,kk) = trqi(3,kk) + uxk*fkd(2) - uyk*fkd(1)
c
c     torque due to induced reaction field gradient on quadrupoles
c
                  fidg(1,1) = -0.25d0 *
     &                           ((sxk*gqxx(2)+syk*gqxx(3)+szk*gqxx(4))
     &                          + (sxk*gux(5)+syk*guy(5)+szk*guz(5)))
                  fidg(1,2) = -0.25d0 *
     &                           ((sxk*gqxy(2)+syk*gqxy(3)+szk*gqxy(4))
     &                          + (sxk*gux(6)+syk*guy(6)+szk*guz(6)))
                  fidg(1,3) = -0.25d0 *
     &                           ((sxk*gqxz(2)+syk*gqxz(3)+szk*gqxz(4))
     &                          + (sxk*gux(7)+syk*guy(7)+szk*guz(7)))
                  fidg(2,2) = -0.25d0 *
     &                           ((sxk*gqyy(2)+syk*gqyy(3)+szk*gqyy(4))
     &                          + (sxk*gux(8)+syk*guy(8)+szk*guz(8)))
                  fidg(2,3) = -0.25d0 *
     &                           ((sxk*gqyz(2)+syk*gqyz(3)+szk*gqyz(4))
     &                          + (sxk*gux(9)+syk*guy(9)+szk*guz(9)))
                  fidg(3,3) = -0.25d0 *
     &                           ((sxk*gqzz(2)+syk*gqzz(3)+szk*gqzz(4))
     &                          + (sxk*gux(10)+syk*guy(10)+szk*guz(10)))
                  fidg(2,1) = fidg(1,2)
                  fidg(3,1) = fidg(1,3)
                  fidg(3,2) = fidg(2,3)
                  fkdg(1,1) = 0.25d0 *
     &                           ((sxi*gqxx(2)+syi*gqxx(3)+szi*gqxx(4))
     &                          + (sxi*gux(5)+syi*guy(5)+szi*guz(5)))
                  fkdg(1,2) = 0.25d0 *
     &                           ((sxi*gqxy(2)+syi*gqxy(3)+szi*gqxy(4))
     &                          + (sxi*gux(6)+syi*guy(6)+szi*guz(6)))
                  fkdg(1,3) = 0.25d0 *
     &                           ((sxi*gqxz(2)+syi*gqxz(3)+szi*gqxz(4))
     &                          + (sxi*gux(7)+syi*guy(7)+szi*guz(7)))
                  fkdg(2,2) = 0.25d0 *
     &                           ((sxi*gqyy(2)+syi*gqyy(3)+szi*gqyy(4))
     &                          + (sxi*gux(8)+syi*guy(8)+szi*guz(8)))
                  fkdg(2,3) = 0.25d0 *
     &                           ((sxi*gqyz(2)+syi*gqyz(3)+szi*gqyz(4))
     &                          + (sxi*gux(9)+syi*guy(9)+szi*guz(9)))
                  fkdg(3,3) = 0.25d0 *
     &                           ((sxi*gqzz(2)+syi*gqzz(3)+szi*gqzz(4))
     &                          + (sxi*gux(10)+syi*guy(10)+szi*guz(10)))
                  fkdg(2,1) = fkdg(1,2)
                  fkdg(3,1) = fkdg(1,3)
                  fkdg(3,2) = fkdg(2,3)
                  if (i .eq. k) then
                     fidg(1,1) = 0.5d0 * fidg(1,1)
                     fidg(1,2) = 0.5d0 * fidg(1,2)
                     fidg(1,3) = 0.5d0 * fidg(1,3)
                     fidg(2,1) = 0.5d0 * fidg(2,1)
                     fidg(2,2) = 0.5d0 * fidg(2,2)
                     fidg(2,3) = 0.5d0 * fidg(2,3)
                     fidg(3,1) = 0.5d0 * fidg(3,1)
                     fidg(3,2) = 0.5d0 * fidg(3,2)
                     fidg(3,3) = 0.5d0 * fidg(3,3)
                     fkdg(1,1) = 0.5d0 * fkdg(1,1)
                     fkdg(1,2) = 0.5d0 * fkdg(1,2)
                     fkdg(1,3) = 0.5d0 * fkdg(1,3)
                     fkdg(2,1) = 0.5d0 * fkdg(2,1)
                     fkdg(2,2) = 0.5d0 * fkdg(2,2)
                     fkdg(2,3) = 0.5d0 * fkdg(2,3)
                     fkdg(3,1) = 0.5d0 * fkdg(3,1)
                     fkdg(3,2) = 0.5d0 * fkdg(3,2)
                     fkdg(3,3) = 0.5d0 * fkdg(3,3)
                  end if
                  trqi(1,ii) = trqi(1,ii) + 2.0d0*
     &                  (qxyi*fidg(1,3)+qyyi*fidg(2,3)+qyzi*fidg(3,3)
     &                  -qxzi*fidg(1,2)-qyzi*fidg(2,2)-qzzi*fidg(3,2))
                  trqi(2,ii) = trqi(2,ii) + 2.0d0*
     &                  (qxzi*fidg(1,1)+qyzi*fidg(2,1)+qzzi*fidg(3,1)
     &                  -qxxi*fidg(1,3)-qxyi*fidg(2,3)-qxzi*fidg(3,3))
                  trqi(3,ii) = trqi(3,ii) + 2.0d0*
     &                  (qxxi*fidg(1,2)+qxyi*fidg(2,2)+qxzi*fidg(3,2)
     &                  -qxyi*fidg(1,1)-qyyi*fidg(2,1)-qyzi*fidg(3,1))
                  trqi(1,kk) = trqi(1,kk) + 2.0d0*
     &                  (qxyk*fkdg(1,3)+qyyk*fkdg(2,3)+qyzk*fkdg(3,3)
     &                  -qxzk*fkdg(1,2)-qyzk*fkdg(2,2)-qzzk*fkdg(3,2))
                  trqi(2,kk) = trqi(2,kk) + 2.0d0*
     &                  (qxzk*fkdg(1,1)+qyzk*fkdg(2,1)+qzzk*fkdg(3,1)
     &                  -qxxk*fkdg(1,3)-qxyk*fkdg(2,3)-qxzk*fkdg(3,3))
                  trqi(3,kk) = trqi(3,kk) + 2.0d0*
     &                  (qxxk*fkdg(1,2)+qxyk*fkdg(2,2)+qxzk*fkdg(3,2)
     &                  -qxyk*fkdg(1,1)-qyyk*fkdg(2,1)-qyzk*fkdg(3,1))
c
c     total permanent and induced energies for this interaction
c
                  e = esym + 0.5d0*(ewi+ewk)
                  ei = 0.5d0 * (esymi + 0.5d0*(ewii+ewki))
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     dedx = dedx * fgrp
                     dedy = dedy * fgrp
                     dedz = dedz * fgrp
                     drbi = drbi * fgrp
                     drbk = drbk * fgrp
                     ei = ei * fgrp
                     dpdx = dpdx * fgrp
                     dpdy = dpdy * fgrp
                     dpdz = dpdz * fgrp
                     dpbi = dpbi * fgrp
                     dpbk = dpbk * fgrp
                  end if
c
c     increment the overall energy and derivative expressions
c
                  if (i .eq. k) then
                     e = 0.5d0 * e
                     ei = 0.5d0 * ei
                     es = es + e + ei
                     drb(i) = drb(i) + drbi
                     drbp(i) = drbp(i) + dpbi
                  else
                     es = es + e + ei
                     des(1,i) = des(1,i) - dedx - dpdx
                     des(2,i) = des(2,i) - dedy - dpdy
                     des(3,i) = des(3,i) - dedz - dpdz
                     des(1,k) = des(1,k) + dedx + dpdx
                     des(2,k) = des(2,k) + dedy + dpdy
                     des(3,k) = des(3,k) + dedz + dpdz
                     drb(i) = drb(i) + drbi
                     drb(k) = drb(k) + drbk
                     drbp(i) = drbp(i) + dpbi
                     drbp(k) = drbp(k) + dpbk
c
c     increment the internal virial tensor components
c
                     vxx = xr * dedx
                     vyx = yr * dedx
                     vzx = zr * dedx
                     vyy = yr * dedy
                     vzy = zr * dedy
                     vzz = zr * dedz
                     vir(1,1) = vir(1,1) + vxx
                     vir(2,1) = vir(2,1) + vyx
                     vir(3,1) = vir(3,1) + vzx
                     vir(1,2) = vir(1,2) + vyx
                     vir(2,2) = vir(2,2) + vyy
                     vir(3,2) = vir(3,2) + vzy
                     vir(1,3) = vir(1,3) + vzx
                     vir(2,3) = vir(2,3) + vzy
                     vir(3,3) = vir(3,3) + vzz
                  end if
c
c     increment the total intermolecular energy
c
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + e + ei
                  end if
               end if
            end if
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     convert torque derivative components to Cartesian form
c
      call torque2 (trq,des)
      call torque2 (trqi,des)
c
c     perform deallocation of some local arrays
c
      deallocate (trq)
      deallocate (trqi)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine ediff1a  --  vacuum to SCRF via pair list  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "ediff1a" calculates the energy and derivatives of polarizing
c     the vacuum induced dipoles to their SCRF polarized values using
c     a double loop
c
c
      subroutine ediff1a
      use sizes
      use atoms
      use bound
      use boxes
      use chgpot
      use couple
      use deriv
      use energi
      use group
      use inter
      use limits
      use molcul
      use mplpot
      use mpole
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer ii,kk
      integer ix,iy,iz
      integer kx,ky,kz
      real*8 ei,f,fgrp
      real*8 damp,gfd
      real*8 scale3,scale5
      real*8 scale7,scale9
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 psc3,psc5,psc7,psc9
      real*8 dsc3,dsc5,dsc7,dsc9
      real*8 scale3i,scale5i
      real*8 scale7i
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9
      real*8 pdi,pti,pgamma
      real*8 ci,di(3),qi(9)
      real*8 ck,dk(3),qk(9)
      real*8 fridmp(3),findmp(3)
      real*8 ftm2i(3)
      real*8 ttm2i(3),ttm3i(3)
      real*8 dixdk(3),fdir(3)
      real*8 dixuk(3),dkxui(3)
      real*8 dixukp(3),dkxuip(3)
      real*8 uixqkr(3),ukxqir(3)
      real*8 uixqkrp(3),ukxqirp(3)
      real*8 qiuk(3),qkui(3)
      real*8 qiukp(3),qkuip(3)
      real*8 rxqiuk(3),rxqkui(3)
      real*8 rxqiukp(3),rxqkuip(3)
      real*8 qidk(3),qkdi(3)
      real*8 qir(3),qkr(3)
      real*8 qiqkr(3),qkqir(3)
      real*8 qixqk(3),rxqir(3)
      real*8 dixr(3),dkxr(3)
      real*8 dixqkr(3),dkxqir(3)
      real*8 rxqkr(3),qkrxqir(3)
      real*8 rxqikr(3),rxqkir(3)
      real*8 rxqidk(3),rxqkdi(3)
      real*8 ddsc3(3),ddsc5(3)
      real*8 ddsc7(3)
      real*8 gli(7),glip(7)
      real*8 sc(10)
      real*8 sci(8),scip(8)
      real*8 gfi(6),gti(6)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: uscale(:)
      real*8, allocatable :: trqi(:,:)
      logical proceed,usei,usek
      character*6 mode
c
c
c     set conversion factor, cutoff and scaling coefficients
c
      if (npole .eq. 0)  return
      f = electric / dielec
      mode = 'MPOLE'
      call switch (mode)
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
      allocate (dscale(n))
      allocate (uscale(n))
      allocate (trqi(3,n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
         uscale(i) = 1.0d0
      end do
c
c     initialize local variables for OpenMP calculation
c
      do i = 1, n
         do j = 1, 3
            trqi(j,i) = 0.0d0
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(npole,ipole,x,y,z,xaxis,yaxis,
!$OMP& zaxis,pdamp,thole,rpole,uind,uinp,uinds,uinps,use,n12,n13,n14,
!$OMP% n15,i12,i13,i14,i15,np11,ip11,np12,ip12,np13,ip13,np14,ip14,
!$OMP% p2scale,p3scale,p4scale,p41scale,p5scale,d1scale,d2scale,
!$OMP% d3scale,d4scale,u1scale,u2scale,u3scale,u4scale,use_group,
!$OMP% use_intra,off2,f)
!$OMP& firstprivate(pscale,dscale,uscale)
!$OMP% shared(es,des,trqi)
!$OMP DO reduction(+:es,des,trqi) schedule(guided)
c
c     calculate the multipole interaction energy and gradient
c
      do i = 1, npole-1
         ii = ipole(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = yaxis(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         di(1) = rpole(2,i)
         di(2) = rpole(3,i)
         di(3) = rpole(4,i)
         qi(1) = rpole(5,i)
         qi(2) = rpole(6,i)
         qi(3) = rpole(7,i)
         qi(4) = rpole(8,i)
         qi(5) = rpole(9,i)
         qi(6) = rpole(10,i)
         qi(7) = rpole(11,i)
         qi(8) = rpole(12,i)
         qi(9) = rpole(13,i)
         usei = (use(ii) .or. use(iz) .or. use(ix) .or. use(iy))
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
                if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
            uscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
            uscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
            uscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
            uscale(ip14(j,ii)) = u4scale
         end do
         do k = i+1, npole
            kk = ipole(k)
            kz = zaxis(k)
            kx = xaxis(k)
            ky = yaxis(k)
            usek = (use(kk) .or. use(kz) .or. use(kx) .or. use(ky))
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (.not. use_intra)  proceed = .true.
            if (proceed)  proceed = (usei .or. usek)
            if (.not. proceed)  goto 10
            ck = rpole(1,k)
            dk(1) = rpole(2,k)
            dk(2) = rpole(3,k)
            dk(3) = rpole(4,k)
            qk(1) = rpole(5,k)
            qk(2) = rpole(6,k)
            qk(3) = rpole(7,k)
            qk(4) = rpole(8,k)
            qk(5) = rpole(9,k)
            qk(6) = rpole(10,k)
            qk(7) = rpole(11,k)
            qk(8) = rpole(12,k)
            qk(9) = rpole(13,k)
            xr = x(kk) - xi
            yr = y(kk) - yi
            zr = z(kk) - zi
            call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
               rr1 = 1.0d0 / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               scale3 = 1.0d0
               scale5 = 1.0d0
               scale7 = 1.0d0
               scale9 = 1.0d0
               do j = 1, 3
                  ddsc3(j) = 0.0d0
                  ddsc5(j) = 0.0d0
                  ddsc7(j) = 0.0d0
               end do
c
c     apply Thole polarization damping to scale factors
c
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     scale3 = 1.0d0 - exp(damp)
                     scale5 = 1.0d0 - (1.0d0-damp)*exp(damp)
                     scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                       *exp(damp)
                     scale9 = 1.0d0 - (1.0d0-damp+(18.0d0*damp**2
     &                                 -9.0d0*damp**3)/35.0d0)*exp(damp)
                     ddsc3(1) = -3.0d0*damp*exp(damp) * xr/r2
                     ddsc3(2) = -3.0d0*damp*exp(damp) * yr/r2
                     ddsc3(3) = -3.0d0*damp*exp(damp) * zr/r2
                     ddsc5(1) = -damp * ddsc3(1)
                     ddsc5(2) = -damp * ddsc3(2)
                     ddsc5(3) = -damp * ddsc3(3)
                     ddsc7(1) = (-0.2d0-0.6d0*damp) * ddsc5(1)
                     ddsc7(2) = (-0.2d0-0.6d0*damp) * ddsc5(2)
                     ddsc7(3) = (-0.2d0-0.6d0*damp) * ddsc5(3)
                  end if
               end if
               scale3i = scale3 * uscale(kk)
               scale5i = scale5 * uscale(kk)
               scale7i = scale7 * uscale(kk)
               dsc3 = scale3 * dscale(kk)
               dsc5 = scale5 * dscale(kk)
               dsc7 = scale7 * dscale(kk)
               dsc9 = scale9 * dscale(kk)
               psc3 = scale3 * pscale(kk)
               psc5 = scale5 * pscale(kk)
               psc7 = scale7 * pscale(kk)
               psc9 = scale9 * pscale(kk)
c
c     construct auxiliary vectors for permanent terms
c
               dixdk(1) = di(2)*dk(3) - di(3)*dk(2)
               dixdk(2) = di(3)*dk(1) - di(1)*dk(3)
               dixdk(3) = di(1)*dk(2) - di(2)*dk(1)
               dixr(1) = di(2)*zr - di(3)*yr
               dixr(2) = di(3)*xr - di(1)*zr
               dixr(3) = di(1)*yr - di(2)*xr
               dkxr(1) = dk(2)*zr - dk(3)*yr
               dkxr(2) = dk(3)*xr - dk(1)*zr
               dkxr(3) = dk(1)*yr - dk(2)*xr
               qir(1) = qi(1)*xr + qi(4)*yr + qi(7)*zr
               qir(2) = qi(2)*xr + qi(5)*yr + qi(8)*zr
               qir(3) = qi(3)*xr + qi(6)*yr + qi(9)*zr
               qkr(1) = qk(1)*xr + qk(4)*yr + qk(7)*zr
               qkr(2) = qk(2)*xr + qk(5)*yr + qk(8)*zr
               qkr(3) = qk(3)*xr + qk(6)*yr + qk(9)*zr
               qiqkr(1) = qi(1)*qkr(1) + qi(4)*qkr(2) + qi(7)*qkr(3)
               qiqkr(2) = qi(2)*qkr(1) + qi(5)*qkr(2) + qi(8)*qkr(3)
               qiqkr(3) = qi(3)*qkr(1) + qi(6)*qkr(2) + qi(9)*qkr(3)
               qkqir(1) = qk(1)*qir(1) + qk(4)*qir(2) + qk(7)*qir(3)
               qkqir(2) = qk(2)*qir(1) + qk(5)*qir(2) + qk(8)*qir(3)
               qkqir(3) = qk(3)*qir(1) + qk(6)*qir(2) + qk(9)*qir(3)
               qixqk(1) = qi(2)*qk(3) + qi(5)*qk(6) + qi(8)*qk(9)
     &                       - qi(3)*qk(2) - qi(6)*qk(5) - qi(9)*qk(8)
               qixqk(2) = qi(3)*qk(1) + qi(6)*qk(4) + qi(9)*qk(7)
     &                       - qi(1)*qk(3) - qi(4)*qk(6) - qi(7)*qk(9)
               qixqk(3) = qi(1)*qk(2) + qi(4)*qk(5) + qi(7)*qk(8)
     &                       - qi(2)*qk(1) - qi(5)*qk(4) - qi(8)*qk(7)
               rxqir(1) = yr*qir(3) - zr*qir(2)
               rxqir(2) = zr*qir(1) - xr*qir(3)
               rxqir(3) = xr*qir(2) - yr*qir(1)
               rxqkr(1) = yr*qkr(3) - zr*qkr(2)
               rxqkr(2) = zr*qkr(1) - xr*qkr(3)
               rxqkr(3) = xr*qkr(2) - yr*qkr(1)
               rxqikr(1) = yr*qiqkr(3) - zr*qiqkr(2)
               rxqikr(2) = zr*qiqkr(1) - xr*qiqkr(3)
               rxqikr(3) = xr*qiqkr(2) - yr*qiqkr(1)
               rxqkir(1) = yr*qkqir(3) - zr*qkqir(2)
               rxqkir(2) = zr*qkqir(1) - xr*qkqir(3)
               rxqkir(3) = xr*qkqir(2) - yr*qkqir(1)
               qkrxqir(1) = qkr(2)*qir(3) - qkr(3)*qir(2)
               qkrxqir(2) = qkr(3)*qir(1) - qkr(1)*qir(3)
               qkrxqir(3) = qkr(1)*qir(2) - qkr(2)*qir(1)
               qidk(1) = qi(1)*dk(1) + qi(4)*dk(2) + qi(7)*dk(3)
               qidk(2) = qi(2)*dk(1) + qi(5)*dk(2) + qi(8)*dk(3)
               qidk(3) = qi(3)*dk(1) + qi(6)*dk(2) + qi(9)*dk(3)
               qkdi(1) = qk(1)*di(1) + qk(4)*di(2) + qk(7)*di(3)
               qkdi(2) = qk(2)*di(1) + qk(5)*di(2) + qk(8)*di(3)
               qkdi(3) = qk(3)*di(1) + qk(6)*di(2) + qk(9)*di(3)
               dixqkr(1) = di(2)*qkr(3) - di(3)*qkr(2)
               dixqkr(2) = di(3)*qkr(1) - di(1)*qkr(3)
               dixqkr(3) = di(1)*qkr(2) - di(2)*qkr(1)
               dkxqir(1) = dk(2)*qir(3) - dk(3)*qir(2)
               dkxqir(2) = dk(3)*qir(1) - dk(1)*qir(3)
               dkxqir(3) = dk(1)*qir(2) - dk(2)*qir(1)
               rxqidk(1) = yr*qidk(3) - zr*qidk(2)
               rxqidk(2) = zr*qidk(1) - xr*qidk(3)
               rxqidk(3) = xr*qidk(2) - yr*qidk(1)
               rxqkdi(1) = yr*qkdi(3) - zr*qkdi(2)
               rxqkdi(2) = zr*qkdi(1) - xr*qkdi(3)
               rxqkdi(3) = xr*qkdi(2) - yr*qkdi(1)
c
c     get intermediate variables for permanent energy terms
c
               sc(3) = di(1)*xr + di(2)*yr + di(3)*zr
               sc(4) = dk(1)*xr + dk(2)*yr + dk(3)*zr
               sc(5) = qir(1)*xr + qir(2)*yr + qir(3)*zr
               sc(6) = qkr(1)*xr + qkr(2)*yr + qkr(3)*zr
c
c     construct auxiliary vectors for induced terms
c
               dixuk(1) = di(2)*uinds(3,k) - di(3)*uinds(2,k)
               dixuk(2) = di(3)*uinds(1,k) - di(1)*uinds(3,k)
               dixuk(3) = di(1)*uinds(2,k) - di(2)*uinds(1,k)
               dkxui(1) = dk(2)*uinds(3,i) - dk(3)*uinds(2,i)
               dkxui(2) = dk(3)*uinds(1,i) - dk(1)*uinds(3,i)
               dkxui(3) = dk(1)*uinds(2,i) - dk(2)*uinds(1,i)
               dixukp(1) = di(2)*uinps(3,k) - di(3)*uinps(2,k)
               dixukp(2) = di(3)*uinps(1,k) - di(1)*uinps(3,k)
               dixukp(3) = di(1)*uinps(2,k) - di(2)*uinps(1,k)
               dkxuip(1) = dk(2)*uinps(3,i) - dk(3)*uinps(2,i)
               dkxuip(2) = dk(3)*uinps(1,i) - dk(1)*uinps(3,i)
               dkxuip(3) = dk(1)*uinps(2,i) - dk(2)*uinps(1,i)
               qiuk(1) = qi(1)*uinds(1,k) + qi(4)*uinds(2,k)
     &                      + qi(7)*uinds(3,k)
               qiuk(2) = qi(2)*uinds(1,k) + qi(5)*uinds(2,k)
     &                      + qi(8)*uinds(3,k)
               qiuk(3) = qi(3)*uinds(1,k) + qi(6)*uinds(2,k)
     &                      + qi(9)*uinds(3,k)
               qkui(1) = qk(1)*uinds(1,i) + qk(4)*uinds(2,i)
     &                      + qk(7)*uinds(3,i)
               qkui(2) = qk(2)*uinds(1,i) + qk(5)*uinds(2,i)
     &                      + qk(8)*uinds(3,i)
               qkui(3) = qk(3)*uinds(1,i) + qk(6)*uinds(2,i)
     &                      + qk(9)*uinds(3,i)
               qiukp(1) = qi(1)*uinps(1,k) + qi(4)*uinps(2,k)
     &                       + qi(7)*uinps(3,k)
               qiukp(2) = qi(2)*uinps(1,k) + qi(5)*uinps(2,k)
     &                       + qi(8)*uinps(3,k)
               qiukp(3) = qi(3)*uinps(1,k) + qi(6)*uinps(2,k)
     &                       + qi(9)*uinps(3,k)
               qkuip(1) = qk(1)*uinps(1,i) + qk(4)*uinps(2,i)
     &                       + qk(7)*uinps(3,i)
               qkuip(2) = qk(2)*uinps(1,i) + qk(5)*uinps(2,i)
     &                       + qk(8)*uinps(3,i)
               qkuip(3) = qk(3)*uinps(1,i) + qk(6)*uinps(2,i)
     &                       + qk(9)*uinps(3,i)
               uixqkr(1) = uinds(2,i)*qkr(3) - uinds(3,i)*qkr(2)
               uixqkr(2) = uinds(3,i)*qkr(1) - uinds(1,i)*qkr(3)
               uixqkr(3) = uinds(1,i)*qkr(2) - uinds(2,i)*qkr(1)
               ukxqir(1) = uinds(2,k)*qir(3) - uinds(3,k)*qir(2)
               ukxqir(2) = uinds(3,k)*qir(1) - uinds(1,k)*qir(3)
               ukxqir(3) = uinds(1,k)*qir(2) - uinds(2,k)*qir(1)
               uixqkrp(1) = uinps(2,i)*qkr(3) - uinps(3,i)*qkr(2)
               uixqkrp(2) = uinps(3,i)*qkr(1) - uinps(1,i)*qkr(3)
               uixqkrp(3) = uinps(1,i)*qkr(2) - uinps(2,i)*qkr(1)
               ukxqirp(1) = uinps(2,k)*qir(3) - uinps(3,k)*qir(2)
               ukxqirp(2) = uinps(3,k)*qir(1) - uinps(1,k)*qir(3)
               ukxqirp(3) = uinps(1,k)*qir(2) - uinps(2,k)*qir(1)
               rxqiuk(1) = yr*qiuk(3) - zr*qiuk(2)
               rxqiuk(2) = zr*qiuk(1) - xr*qiuk(3)
               rxqiuk(3) = xr*qiuk(2) - yr*qiuk(1)
               rxqkui(1) = yr*qkui(3) - zr*qkui(2)
               rxqkui(2) = zr*qkui(1) - xr*qkui(3)
               rxqkui(3) = xr*qkui(2) - yr*qkui(1)
               rxqiukp(1) = yr*qiukp(3) - zr*qiukp(2)
               rxqiukp(2) = zr*qiukp(1) - xr*qiukp(3)
               rxqiukp(3) = xr*qiukp(2) - yr*qiukp(1)
               rxqkuip(1) = yr*qkuip(3) - zr*qkuip(2)
               rxqkuip(2) = zr*qkuip(1) - xr*qkuip(3)
               rxqkuip(3) = xr*qkuip(2) - yr*qkuip(1)
c
c     get intermediate variables for induction energy terms
c
               sci(1) = uinds(1,i)*dk(1) + uinds(2,i)*dk(2)
     &                     + uinds(3,i)*dk(3) + di(1)*uinds(1,k)
     &                     + di(2)*uinds(2,k) + di(3)*uinds(3,k)
               sci(2) = uinds(1,i)*uinds(1,k) + uinds(2,i)*uinds(2,k)
     &                     + uinds(3,i)*uinds(3,k)
               sci(3) = uinds(1,i)*xr + uinds(2,i)*yr + uinds(3,i)*zr
               sci(4) = uinds(1,k)*xr + uinds(2,k)*yr + uinds(3,k)*zr
               sci(7) = qir(1)*uinds(1,k) + qir(2)*uinds(2,k)
     &                     + qir(3)*uinds(3,k)
               sci(8) = qkr(1)*uinds(1,i) + qkr(2)*uinds(2,i)
     &                     + qkr(3)*uinds(3,i)
               scip(1) = uinps(1,i)*dk(1) + uinps(2,i)*dk(2)
     &                      + uinps(3,i)*dk(3) + di(1)*uinps(1,k)
     &                      + di(2)*uinps(2,k) + di(3)*uinps(3,k)
               scip(2) = uinds(1,i)*uinps(1,k) + uinds(2,i)*uinps(2,k)
     &                 + uinds(3,i)*uinps(3,k) + uinps(1,i)*uinds(1,k)
     &                 + uinps(2,i)*uinds(2,k) + uinps(3,i)*uinds(3,k)
               scip(3) = uinps(1,i)*xr + uinps(2,i)*yr + uinps(3,i)*zr
               scip(4) = uinps(1,k)*xr + uinps(2,k)*yr + uinps(3,k)*zr
               scip(7) = qir(1)*uinps(1,k) + qir(2)*uinps(2,k)
     &                      + qir(3)*uinps(3,k)
               scip(8) = qkr(1)*uinps(1,i) + qkr(2)*uinps(2,i)
     &                      + qkr(3)*uinps(3,i)
c
c     calculate the gl functions for potential energy
c
               gli(1) = ck*sci(3) - ci*sci(4)
               gli(2) = -sc(3)*sci(4) - sci(3)*sc(4)
               gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
               gli(6) = sci(1)
               gli(7) = 2.0d0 * (sci(7)-sci(8))
               glip(1) = ck*scip(3) - ci*scip(4)
               glip(2) = -sc(3)*scip(4) - scip(3)*sc(4)
               glip(3) = scip(3)*sc(6) - scip(4)*sc(5)
               glip(6) = scip(1)
               glip(7) = 2.0d0 * (scip(7)-scip(8))
c
c     get the permanent multipole and induced energies
c
               ei = 0.5d0 * (rr3*(gli(1)+gli(6))*psc3
     &                          + rr5*(gli(2)+gli(7))*psc5
     &                          + rr7*gli(3)*psc7)
               ei = f * ei
               es = es + ei
c
c     intermediate variables for the induced-permanent terms
c
               gfi(1) = 0.5d0*rr5*((gli(1)+gli(6))*psc3
     &                     +(glip(1)+glip(6))*dsc3+scip(2)*scale3i)
     &                + 0.5d0*rr7*((gli(7)+gli(2))*psc5
     &                            +(glip(7)+glip(2))*dsc5
     &                     -(sci(3)*scip(4)+scip(3)*sci(4))*scale5i)
     &                + 0.5d0*rr9*(gli(3)*psc7+glip(3)*dsc7)
               gfi(2) = -rr3*ck + rr5*sc(4) - rr7*sc(6)
               gfi(3) = rr3*ci + rr5*sc(3) + rr7*sc(5)
               gfi(4) = 2.0d0 * rr5
               gfi(5) = rr7 * (sci(4)*psc7+scip(4)*dsc7)
               gfi(6) = -rr7 * (sci(3)*psc7+scip(3)*dsc7)
c
c     get the induced force
c
               ftm2i(1) = gfi(1)*xr + 0.5d0*
     &            (- rr3*ck*(uinds(1,i)*psc3+uinps(1,i)*dsc3)
     &             + rr5*sc(4)*(uinds(1,i)*psc5+uinps(1,i)*dsc5)
     &             - rr7*sc(6)*(uinds(1,i)*psc7+uinps(1,i)*dsc7))
     &             +(rr3*ci*(uinds(1,k)*psc3+uinps(1,k)*dsc3)
     &             + rr5*sc(3)*(uinds(1,k)*psc5+uinps(1,k)*dsc5)
     &             + rr7*sc(5)*(uinds(1,k)*psc7+uinps(1,k)*dsc7))*0.5d0
     &             + rr5*scale5i*(sci(4)*uinps(1,i)+scip(4)*uinds(1,i)
     &             + sci(3)*uinps(1,k)+scip(3)*uinds(1,k))*0.5d0
     &             + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(1)
     &             + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(1)
     &             + 0.5d0*gfi(4)*((qkui(1)-qiuk(1))*psc5
     &             + (qkuip(1)-qiukp(1))*dsc5)
     &             + gfi(5)*qir(1) + gfi(6)*qkr(1)
               ftm2i(2) = gfi(1)*yr + 0.5d0*
     &            (- rr3*ck*(uinds(2,i)*psc3+uinps(2,i)*dsc3)
     &             + rr5*sc(4)*(uinds(2,i)*psc5+uinps(2,i)*dsc5)
     &             - rr7*sc(6)*(uinds(2,i)*psc7+uinps(2,i)*dsc7))
     &             +(rr3*ci*(uinds(2,k)*psc3+uinps(2,k)*dsc3)
     &             + rr5*sc(3)*(uinds(2,k)*psc5+uinps(2,k)*dsc5)
     &             + rr7*sc(5)*(uinds(2,k)*psc7+uinps(2,k)*dsc7))*0.5d0
     &             + rr5*scale5i*(sci(4)*uinps(2,i)+scip(4)*uinds(2,i)
     &             + sci(3)*uinps(2,k)+scip(3)*uinds(2,k))*0.5d0
     &             + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(2)
     &             + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(2)
     &             + 0.5d0*gfi(4)*((qkui(2)-qiuk(2))*psc5
     &             + (qkuip(2)-qiukp(2))*dsc5)
     &             + gfi(5)*qir(2) + gfi(6)*qkr(2)
               ftm2i(3) = gfi(1)*zr  + 0.5d0*
     &            (- rr3*ck*(uinds(3,i)*psc3+uinps(3,i)*dsc3)
     &             + rr5*sc(4)*(uinds(3,i)*psc5+uinps(3,i)*dsc5)
     &             - rr7*sc(6)*(uinds(3,i)*psc7+uinps(3,i)*dsc7))
     &             +(rr3*ci*(uinds(3,k)*psc3+uinps(3,k)*dsc3)
     &             + rr5*sc(3)*(uinds(3,k)*psc5+uinps(3,k)*dsc5)
     &             + rr7*sc(5)*(uinds(3,k)*psc7+uinps(3,k)*dsc7))*0.5d0
     &             + rr5*scale5i*(sci(4)*uinps(3,i)+scip(4)*uinds(3,i)
     &             + sci(3)*uinps(3,k)+scip(3)*uinds(3,k))*0.5d0
     &             + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(3)
     &             + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(3)
     &             + 0.5d0*gfi(4)*((qkui(3)-qiuk(3))*psc5
     &             + (qkuip(3)-qiukp(3))*dsc5)
     &             + gfi(5)*qir(3) + gfi(6)*qkr(3)
c
c     intermediate values needed for partially excluded interactions
c
               fridmp(1) = 0.5d0 * (rr3*((gli(1)+gli(6))*pscale(kk)
     &                        +(glip(1)+glip(6))*dscale(kk))*ddsc3(1)
     &            + rr5*((gli(2)+gli(7))*pscale(kk)
     &                +(glip(2)+glip(7))*dscale(kk))*ddsc5(1)
     &            + rr7*(gli(3)*pscale(kk)+glip(3)*dscale(kk))*ddsc7(1))
               fridmp(2) = 0.5d0 * (rr3*((gli(1)+gli(6))*pscale(kk)
     &                        +(glip(1)+glip(6))*dscale(kk))*ddsc3(2)
     &            + rr5*((gli(2)+gli(7))*pscale(kk)
     &                +(glip(2)+glip(7))*dscale(kk))*ddsc5(2)
     &            + rr7*(gli(3)*pscale(kk)+glip(3)*dscale(kk))*ddsc7(2))
               fridmp(3) = 0.5d0 * (rr3*((gli(1)+gli(6))*pscale(kk)
     &                        +(glip(1)+glip(6))*dscale(kk))*ddsc3(3)
     &            + rr5*((gli(2)+gli(7))*pscale(kk)
     &                +(glip(2)+glip(7))*dscale(kk))*ddsc5(3)
     &            + rr7*(gli(3)*pscale(kk)+glip(3)*dscale(kk))*ddsc7(3))
c
c     get the induced-induced derivative terms
c
               findmp(1) = 0.5d0 * uscale(kk) * (scip(2)*rr3*ddsc3(1)
     &                   - rr5*ddsc5(1)*(sci(3)*scip(4)+scip(3)*sci(4)))
               findmp(2) = 0.5d0 * uscale(kk) * (scip(2)*rr3*ddsc3(2)
     &                   - rr5*ddsc5(2)*(sci(3)*scip(4)+scip(3)*sci(4)))
               findmp(3) = 0.5d0 * uscale(kk) * (scip(2)*rr3*ddsc3(3)
     &                   - rr5*ddsc5(3)*(sci(3)*scip(4)+scip(3)*sci(4)))
c
c     handle of scaling for partially excluded interactions
c
               ftm2i(1) = ftm2i(1) - fridmp(1) - findmp(1)
               ftm2i(2) = ftm2i(2) - fridmp(2) - findmp(2)
               ftm2i(3) = ftm2i(3) - fridmp(3) - findmp(3)
c
c     correction to convert mutual to direct polarization force
c
               if (poltyp .eq. 'DIRECT') then
                  gfd = 0.5d0 * (rr5*scip(2)*scale3i
     &                  - rr7*(scip(3)*sci(4)+sci(3)*scip(4))*scale5i)
                  fdir(1) = gfd*xr + 0.5d0*rr5*scale5i
     &                         * (sci(4)*uinps(1,i)+scip(4)*uinds(1,i)
     &                           +sci(3)*uinps(1,k)+scip(3)*uinds(1,k))
                  fdir(2) = gfd*yr + 0.5d0*rr5*scale5i
     &                         * (sci(4)*uinps(2,i)+scip(4)*uinds(2,i)
     &                           +sci(3)*uinps(2,k)+scip(3)*uinds(2,k))
                  fdir(3) = gfd*zr + 0.5d0*rr5*scale5i
     &                         * (sci(4)*uinps(3,i)+scip(4)*uinds(3,i)
     &                           +sci(3)*uinps(3,k)+scip(3)*uinds(3,k))
                  ftm2i(1) = ftm2i(1) - fdir(1) + findmp(1)
                  ftm2i(2) = ftm2i(2) - fdir(2) + findmp(2)
                  ftm2i(3) = ftm2i(3) - fdir(3) + findmp(3)
               end if
c
c     now perform the torque calculation
c     intermediate terms for torque between multipoles i and k
c
               gti(2) = 0.5d0 * (sci(4)*psc5+scip(4)*dsc5) * rr5
               gti(3) = 0.5d0 * (sci(3)*psc5+scip(3)*dsc5) * rr5
               gti(4) = gfi(4)
               gti(5) = gfi(5)
               gti(6) = gfi(6)
c
c     calculate the induced torque components
c
               ttm2i(1) = -rr3*(dixuk(1)*psc3+dixukp(1)*dsc3)*0.5d0
     &            + gti(2)*dixr(1) + gti(4)*((ukxqir(1)+rxqiuk(1))*psc5
     &            +(ukxqirp(1)+rxqiukp(1))*dsc5)*0.5d0 - gti(5)*rxqir(1)
               ttm2i(2) = -rr3*(dixuk(2)*psc3+dixukp(2)*dsc3)*0.5d0
     &            + gti(2)*dixr(2) + gti(4)*((ukxqir(2)+rxqiuk(2))*psc5
     &            +(ukxqirp(2)+rxqiukp(2))*dsc5)*0.5d0 - gti(5)*rxqir(2)
               ttm2i(3) = -rr3*(dixuk(3)*psc3+dixukp(3)*dsc3)*0.5d0
     &            + gti(2)*dixr(3) + gti(4)*((ukxqir(3)+rxqiuk(3))*psc5
     &            +(ukxqirp(3)+rxqiukp(3))*dsc5)*0.5d0 - gti(5)*rxqir(3)
               ttm3i(1) = -rr3*(dkxui(1)*psc3+dkxuip(1)*dsc3)*0.5d0
     &            + gti(3)*dkxr(1) - gti(4)*((uixqkr(1)+rxqkui(1))*psc5
     &            +(uixqkrp(1)+rxqkuip(1))*dsc5)*0.5d0 - gti(6)*rxqkr(1)
               ttm3i(2) = -rr3*(dkxui(2)*psc3+dkxuip(2)*dsc3)*0.5d0
     &            + gti(3)*dkxr(2) - gti(4)*((uixqkr(2)+rxqkui(2))*psc5
     &            +(uixqkrp(2)+rxqkuip(2))*dsc5)*0.5d0 - gti(6)*rxqkr(2)
               ttm3i(3) = -rr3*(dkxui(3)*psc3+dkxuip(3)*dsc3)*0.5d0
     &            + gti(3)*dkxr(3) - gti(4)*((uixqkr(3)+rxqkui(3))*psc5
     &            +(uixqkrp(3)+rxqkuip(3))*dsc5)*0.5d0 - gti(6)*rxqkr(3)
c
c     update the force components on sites i and k
c
               des(1,ii) = des(1,ii) + f*ftm2i(1)
               des(2,ii) = des(2,ii) + f*ftm2i(2)
               des(3,ii) = des(3,ii) + f*ftm2i(3)
               des(1,kk) = des(1,kk) - f*ftm2i(1)
               des(2,kk) = des(2,kk) - f*ftm2i(2)
               des(3,kk) = des(3,kk) - f*ftm2i(3)
c
c     update the torque components on sites i and k
c
               trqi(1,ii) = trqi(1,ii) + f*ttm2i(1)
               trqi(2,ii) = trqi(2,ii) + f*ttm2i(2)
               trqi(3,ii) = trqi(3,ii) + f*ttm2i(3)
               trqi(1,kk) = trqi(1,kk) + f*ttm3i(1)
               trqi(2,kk) = trqi(2,kk) + f*ttm3i(2)
               trqi(3,kk) = trqi(3,kk) + f*ttm3i(3)
c
c     construct auxiliary vectors for induced terms
c
               dixuk(1) = di(2)*uind(3,k) - di(3)*uind(2,k)
               dixuk(2) = di(3)*uind(1,k) - di(1)*uind(3,k)
               dixuk(3) = di(1)*uind(2,k) - di(2)*uind(1,k)
               dkxui(1) = dk(2)*uind(3,i) - dk(3)*uind(2,i)
               dkxui(2) = dk(3)*uind(1,i) - dk(1)*uind(3,i)
               dkxui(3) = dk(1)*uind(2,i) - dk(2)*uind(1,i)
               dixukp(1) = di(2)*uinp(3,k) - di(3)*uinp(2,k)
               dixukp(2) = di(3)*uinp(1,k) - di(1)*uinp(3,k)
               dixukp(3) = di(1)*uinp(2,k) - di(2)*uinp(1,k)
               dkxuip(1) = dk(2)*uinp(3,i) - dk(3)*uinp(2,i)
               dkxuip(2) = dk(3)*uinp(1,i) - dk(1)*uinp(3,i)
               dkxuip(3) = dk(1)*uinp(2,i) - dk(2)*uinp(1,i)
               qiuk(1) = qi(1)*uind(1,k) + qi(4)*uind(2,k)
     &                      + qi(7)*uind(3,k)
               qiuk(2) = qi(2)*uind(1,k) + qi(5)*uind(2,k)
     &                      + qi(8)*uind(3,k)
               qiuk(3) = qi(3)*uind(1,k) + qi(6)*uind(2,k)
     &                      + qi(9)*uind(3,k)
               qkui(1) = qk(1)*uind(1,i) + qk(4)*uind(2,i)
     &                      + qk(7)*uind(3,i)
               qkui(2) = qk(2)*uind(1,i) + qk(5)*uind(2,i)
     &                      + qk(8)*uind(3,i)
               qkui(3) = qk(3)*uind(1,i) + qk(6)*uind(2,i)
     &                      + qk(9)*uind(3,i)
               qiukp(1) = qi(1)*uinp(1,k) + qi(4)*uinp(2,k)
     &                       + qi(7)*uinp(3,k)
               qiukp(2) = qi(2)*uinp(1,k) + qi(5)*uinp(2,k)
     &                       + qi(8)*uinp(3,k)
               qiukp(3) = qi(3)*uinp(1,k) + qi(6)*uinp(2,k)
     &                       + qi(9)*uinp(3,k)
               qkuip(1) = qk(1)*uinp(1,i) + qk(4)*uinp(2,i)
     &                       + qk(7)*uinp(3,i)
               qkuip(2) = qk(2)*uinp(1,i) + qk(5)*uinp(2,i)
     &                       + qk(8)*uinp(3,i)
               qkuip(3) = qk(3)*uinp(1,i) + qk(6)*uinp(2,i)
     &                       + qk(9)*uinp(3,i)
               uixqkr(1) = uind(2,i)*qkr(3) - uind(3,i)*qkr(2)
               uixqkr(2) = uind(3,i)*qkr(1) - uind(1,i)*qkr(3)
               uixqkr(3) = uind(1,i)*qkr(2) - uind(2,i)*qkr(1)
               ukxqir(1) = uind(2,k)*qir(3) - uind(3,k)*qir(2)
               ukxqir(2) = uind(3,k)*qir(1) - uind(1,k)*qir(3)
               ukxqir(3) = uind(1,k)*qir(2) - uind(2,k)*qir(1)
               uixqkrp(1) = uinp(2,i)*qkr(3) - uinp(3,i)*qkr(2)
               uixqkrp(2) = uinp(3,i)*qkr(1) - uinp(1,i)*qkr(3)
               uixqkrp(3) = uinp(1,i)*qkr(2) - uinp(2,i)*qkr(1)
               ukxqirp(1) = uinp(2,k)*qir(3) - uinp(3,k)*qir(2)
               ukxqirp(2) = uinp(3,k)*qir(1) - uinp(1,k)*qir(3)
               ukxqirp(3) = uinp(1,k)*qir(2) - uinp(2,k)*qir(1)
               rxqiuk(1) = yr*qiuk(3) - zr*qiuk(2)
               rxqiuk(2) = zr*qiuk(1) - xr*qiuk(3)
               rxqiuk(3) = xr*qiuk(2) - yr*qiuk(1)
               rxqkui(1) = yr*qkui(3) - zr*qkui(2)
               rxqkui(2) = zr*qkui(1) - xr*qkui(3)
               rxqkui(3) = xr*qkui(2) - yr*qkui(1)
               rxqiukp(1) = yr*qiukp(3) - zr*qiukp(2)
               rxqiukp(2) = zr*qiukp(1) - xr*qiukp(3)
               rxqiukp(3) = xr*qiukp(2) - yr*qiukp(1)
               rxqkuip(1) = yr*qkuip(3) - zr*qkuip(2)
               rxqkuip(2) = zr*qkuip(1) - xr*qkuip(3)
               rxqkuip(3) = xr*qkuip(2) - yr*qkuip(1)
c
c     get intermediate variables for induction energy terms
c
               sci(1) = uind(1,i)*dk(1) + uind(2,i)*dk(2)
     &                     + uind(3,i)*dk(3) + di(1)*uind(1,k)
     &                     + di(2)*uind(2,k) + di(3)*uind(3,k)
               sci(2) = uind(1,i)*uind(1,k) + uind(2,i)*uind(2,k)
     &                     + uind(3,i)*uind(3,k)
               sci(3) = uind(1,i)*xr + uind(2,i)*yr + uind(3,i)*zr
               sci(4) = uind(1,k)*xr + uind(2,k)*yr + uind(3,k)*zr
               sci(7) = qir(1)*uind(1,k) + qir(2)*uind(2,k)
     &                     + qir(3)*uind(3,k)
               sci(8) = qkr(1)*uind(1,i) + qkr(2)*uind(2,i)
     &                     + qkr(3)*uind(3,i)
               scip(1) = uinp(1,i)*dk(1) + uinp(2,i)*dk(2)
     &                      + uinp(3,i)*dk(3) + di(1)*uinp(1,k)
     &                      + di(2)*uinp(2,k) + di(3)*uinp(3,k)
               scip(2) = uind(1,i)*uinp(1,k)+uind(2,i)*uinp(2,k)
     &                      + uind(3,i)*uinp(3,k)+uinp(1,i)*uind(1,k)
     &                      + uinp(2,i)*uind(2,k)+uinp(3,i)*uind(3,k)
               scip(3) = uinp(1,i)*xr + uinp(2,i)*yr + uinp(3,i)*zr
               scip(4) = uinp(1,k)*xr + uinp(2,k)*yr + uinp(3,k)*zr
               scip(7) = qir(1)*uinp(1,k) + qir(2)*uinp(2,k)
     &                      + qir(3)*uinp(3,k)
               scip(8) = qkr(1)*uinp(1,i) + qkr(2)*uinp(2,i)
     &                      + qkr(3)*uinp(3,i)
c
c     calculate the gl functions for potential energy
c
               gli(1) = ck*sci(3) - ci*sci(4)
               gli(2) = -sc(3)*sci(4) - sci(3)*sc(4)
               gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
               gli(6) = sci(1)
               gli(7) = 2.0d0 * (sci(7)-sci(8))
               glip(1) = ck*scip(3) - ci*scip(4)
               glip(2) = -sc(3)*scip(4) - scip(3)*sc(4)
               glip(3) = scip(3)*sc(6) - scip(4)*sc(5)
               glip(6) = scip(1)
               glip(7) = 2.0d0 * (scip(7)-scip(8))
c
c     get the permanent multipole and induced energies
c
               ei = -0.5d0 * (rr3*(gli(1)+gli(6))*psc3
     &                           + rr5*(gli(2)+gli(7))*psc5
     &                           + rr7*gli(3)*psc7)
               ei = f * ei
               es = es + ei
c
c     intermediate variables for the induced-permanent terms
c
               gfi(1) = 0.5d0*rr5*((gli(1)+gli(6))*psc3
     &                     +(glip(1)+glip(6))*dsc3+scip(2)*scale3i)
     &                + 0.5d0*rr7*((gli(7)+gli(2))*psc5
     &                            +(glip(7)+glip(2))*dsc5
     &                     -(sci(3)*scip(4)+scip(3)*sci(4))*scale5i)
     &                + 0.5d0*rr9*(gli(3)*psc7+glip(3)*dsc7)
               gfi(2) = -rr3*ck + rr5*sc(4) - rr7*sc(6)
               gfi(3) = rr3*ci + rr5*sc(3) + rr7*sc(5)
               gfi(4) = 2.0d0 * rr5
               gfi(5) = rr7 * (sci(4)*psc7+scip(4)*dsc7)
               gfi(6) = -rr7 * (sci(3)*psc7+scip(3)*dsc7)
c
c     get the induced force
c
               ftm2i(1) = gfi(1)*xr + 0.5d0*
     &            (- rr3*ck*(uind(1,i)*psc3+uinp(1,i)*dsc3)
     &             + rr5*sc(4)*(uind(1,i)*psc5+uinp(1,i)*dsc5)
     &             - rr7*sc(6)*(uind(1,i)*psc7+uinp(1,i)*dsc7))
     &             +(rr3*ci*(uind(1,k)*psc3+uinp(1,k)*dsc3)
     &             + rr5*sc(3)*(uind(1,k)*psc5+uinp(1,k)*dsc5)
     &             + rr7*sc(5)*(uind(1,k)*psc7+uinp(1,k)*dsc7))*0.5d0
     &             + rr5*scale5i*(sci(4)*uinp(1,i)+scip(4)*uind(1,i)
     &             + sci(3)*uinp(1,k)+scip(3)*uind(1,k))*0.5d0
     &             + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(1)
     &             + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(1)
     &             + 0.5d0*gfi(4)*((qkui(1)-qiuk(1))*psc5
     &             + (qkuip(1)-qiukp(1))*dsc5)
     &             + gfi(5)*qir(1) + gfi(6)*qkr(1)
               ftm2i(2) = gfi(1)*yr + 0.5d0*
     &            (- rr3*ck*(uind(2,i)*psc3+uinp(2,i)*dsc3)
     &             + rr5*sc(4)*(uind(2,i)*psc5+uinp(2,i)*dsc5)
     &             - rr7*sc(6)*(uind(2,i)*psc7+uinp(2,i)*dsc7))
     &             +(rr3*ci*(uind(2,k)*psc3+uinp(2,k)*dsc3)
     &             + rr5*sc(3)*(uind(2,k)*psc5+uinp(2,k)*dsc5)
     &             + rr7*sc(5)*(uind(2,k)*psc7+uinp(2,k)*dsc7))*0.5d0
     &             + rr5*scale5i*(sci(4)*uinp(2,i)+scip(4)*uind(2,i)
     &             + sci(3)*uinp(2,k)+scip(3)*uind(2,k))*0.5d0
     &             + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(2)
     &             + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(2)
     &             + 0.5d0*gfi(4)*((qkui(2)-qiuk(2))*psc5
     &             + (qkuip(2)-qiukp(2))*dsc5)
     &             + gfi(5)*qir(2) + gfi(6)*qkr(2)
               ftm2i(3) = gfi(1)*zr  + 0.5d0*
     &            (- rr3*ck*(uind(3,i)*psc3+uinp(3,i)*dsc3)
     &             + rr5*sc(4)*(uind(3,i)*psc5+uinp(3,i)*dsc5)
     &             - rr7*sc(6)*(uind(3,i)*psc7+uinp(3,i)*dsc7))
     &             +(rr3*ci*(uind(3,k)*psc3+uinp(3,k)*dsc3)
     &             + rr5*sc(3)*(uind(3,k)*psc5+uinp(3,k)*dsc5)
     &             + rr7*sc(5)*(uind(3,k)*psc7+uinp(3,k)*dsc7))*0.5d0
     &             + rr5*scale5i*(sci(4)*uinp(3,i)+scip(4)*uind(3,i)
     &             + sci(3)*uinp(3,k)+scip(3)*uind(3,k))*0.5d0
     &             + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(3)
     &             + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(3)
     &             + 0.5d0*gfi(4)*((qkui(3)-qiuk(3))*psc5
     &             + (qkuip(3)-qiukp(3))*dsc5)
     &             + gfi(5)*qir(3) + gfi(6)*qkr(3)
c
c     intermediate values needed for partially excluded interactions
c
               fridmp(1) = 0.5d0 * (rr3*((gli(1)+gli(6))*pscale(kk)
     &                        +(glip(1)+glip(6))*dscale(kk))*ddsc3(1)
     &            + rr5*((gli(2)+gli(7))*pscale(kk)
     &                +(glip(2)+glip(7))*dscale(kk))*ddsc5(1)
     &            + rr7*(gli(3)*pscale(kk)+glip(3)*dscale(kk))*ddsc7(1))
               fridmp(2) = 0.5d0 * (rr3*((gli(1)+gli(6))*pscale(kk)
     &                        +(glip(1)+glip(6))*dscale(kk))*ddsc3(2)
     &            + rr5*((gli(2)+gli(7))*pscale(kk)
     &                +(glip(2)+glip(7))*dscale(kk))*ddsc5(2)
     &            + rr7*(gli(3)*pscale(kk)+glip(3)*dscale(kk))*ddsc7(2))
               fridmp(3) = 0.5d0 * (rr3*((gli(1)+gli(6))*pscale(kk)
     &                        +(glip(1)+glip(6))*dscale(kk))*ddsc3(3)
     &            + rr5*((gli(2)+gli(7))*pscale(kk)
     &                +(glip(2)+glip(7))*dscale(kk))*ddsc5(3)
     &            + rr7*(gli(3)*pscale(kk)+glip(3)*dscale(kk))*ddsc7(3))
c
c     get the induced-induced derivative terms
c
               findmp(1) = 0.5d0 * uscale(kk) * (scip(2)*rr3*ddsc3(1)
     &                   - rr5*ddsc5(1)*(sci(3)*scip(4)+scip(3)*sci(4)))
               findmp(2) = 0.5d0 * uscale(kk) * (scip(2)*rr3*ddsc3(2)
     &                   - rr5*ddsc5(2)*(sci(3)*scip(4)+scip(3)*sci(4)))
               findmp(3) = 0.5d0 * uscale(kk) * (scip(2)*rr3*ddsc3(3)
     &                   - rr5*ddsc5(3)*(sci(3)*scip(4)+scip(3)*sci(4)))
c
c     handle of scaling for partially excluded interactions
c
               ftm2i(1) = ftm2i(1) - fridmp(1) - findmp(1)
               ftm2i(2) = ftm2i(2) - fridmp(2) - findmp(2)
               ftm2i(3) = ftm2i(3) - fridmp(3) - findmp(3)
c
c     correction to convert mutual to direct polarization force
c
               if (poltyp .eq. 'DIRECT') then
                  gfd = 0.5d0 * (rr5*scip(2)*scale3i
     &                  - rr7*(scip(3)*sci(4)+sci(3)*scip(4))*scale5i)
                  fdir(1) = gfd*xr + 0.5d0*rr5*scale5i
     &                         * (sci(4)*uinp(1,i)+scip(4)*uind(1,i)
     &                           +sci(3)*uinp(1,k)+scip(3)*uind(1,k))
                  fdir(2) = gfd*yr + 0.5d0*rr5*scale5i
     &                         * (sci(4)*uinp(2,i)+scip(4)*uind(2,i)
     &                           +sci(3)*uinp(2,k)+scip(3)*uind(2,k))
                  fdir(3) = gfd*zr + 0.5d0*rr5*scale5i
     &                         * (sci(4)*uinp(3,i)+scip(4)*uind(3,i)
     &                           +sci(3)*uinp(3,k)+scip(3)*uind(3,k))
                  ftm2i(1) = ftm2i(1) - fdir(1) + findmp(1)
                  ftm2i(2) = ftm2i(2) - fdir(2) + findmp(2)
                  ftm2i(3) = ftm2i(3) - fdir(3) + findmp(3)
               end if
c
c     now perform the torque calculation
c     intermediate terms for torque between multipoles i and k
c
               gti(2) = 0.5d0 * (sci(4)*psc5+scip(4)*dsc5) * rr5
               gti(3) = 0.5d0 * (sci(3)*psc5+scip(3)*dsc5) * rr5
               gti(4) = gfi(4)
               gti(5) = gfi(5)
               gti(6) = gfi(6)
c
c     calculate the induced torque components
c
               ttm2i(1) = -rr3*(dixuk(1)*psc3+dixukp(1)*dsc3)*0.5d0
     &            + gti(2)*dixr(1) + gti(4)*((ukxqir(1)+rxqiuk(1))*psc5
     &            +(ukxqirp(1)+rxqiukp(1))*dsc5)*0.5d0 - gti(5)*rxqir(1)
               ttm2i(2) = -rr3*(dixuk(2)*psc3+dixukp(2)*dsc3)*0.5d0
     &            + gti(2)*dixr(2) + gti(4)*((ukxqir(2)+rxqiuk(2))*psc5
     &            +(ukxqirp(2)+rxqiukp(2))*dsc5)*0.5d0 - gti(5)*rxqir(2)
               ttm2i(3) = -rr3*(dixuk(3)*psc3+dixukp(3)*dsc3)*0.5d0
     &            + gti(2)*dixr(3) + gti(4)*((ukxqir(3)+rxqiuk(3))*psc5
     &            +(ukxqirp(3)+rxqiukp(3))*dsc5)*0.5d0 - gti(5)*rxqir(3)
               ttm3i(1) = -rr3*(dkxui(1)*psc3+dkxuip(1)*dsc3)*0.5d0
     &            + gti(3)*dkxr(1) - gti(4)*((uixqkr(1)+rxqkui(1))*psc5
     &            +(uixqkrp(1)+rxqkuip(1))*dsc5)*0.5d0 - gti(6)*rxqkr(1)
               ttm3i(2) = -rr3*(dkxui(2)*psc3+dkxuip(2)*dsc3)*0.5d0
     &            + gti(3)*dkxr(2) - gti(4)*((uixqkr(2)+rxqkui(2))*psc5
     &            +(uixqkrp(2)+rxqkuip(2))*dsc5)*0.5d0 - gti(6)*rxqkr(2)
               ttm3i(3) = -rr3*(dkxui(3)*psc3+dkxuip(3)*dsc3)*0.5d0
     &            + gti(3)*dkxr(3) - gti(4)*((uixqkr(3)+rxqkui(3))*psc5
     &            +(uixqkrp(3)+rxqkuip(3))*dsc5)*0.5d0 - gti(6)*rxqkr(3)
c
c     update the force components on sites i and k
c
               des(1,ii) = des(1,ii) - f*ftm2i(1)
               des(2,ii) = des(2,ii) - f*ftm2i(2)
               des(3,ii) = des(3,ii) - f*ftm2i(3)
               des(1,kk) = des(1,kk) + f*ftm2i(1)
               des(2,kk) = des(2,kk) + f*ftm2i(2)
               des(3,kk) = des(3,kk) + f*ftm2i(3)
c
c     update the torque components on sites i and k
c
               trqi(1,ii) = trqi(1,ii) - f*ttm2i(1)
               trqi(2,ii) = trqi(2,ii) - f*ttm2i(2)
               trqi(3,ii) = trqi(3,ii) - f*ttm2i(3)
               trqi(1,kk) = trqi(1,kk) - f*ttm3i(1)
               trqi(2,kk) = trqi(2,kk) - f*ttm3i(2)
               trqi(3,kk) = trqi(3,kk) - f*ttm3i(3)
            end if
   10       continue
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = 1.0d0
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = 1.0d0
            uscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = 1.0d0
            uscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = 1.0d0
            uscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = 1.0d0
            uscale(ip14(j,ii)) = 1.0d0
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     convert torque values into Cartesian forces
c
      call torque2 (trqi,des)
c
c     perform deallocation of some local arrays
c
      deallocate (pscale)
      deallocate (dscale)
      deallocate (uscale)
      deallocate (trqi)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine ediff1b  --  vacuum to SCRF via pair list  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "ediff1b" calculates the energy and derivatives of polarizing
c     the vacuum induced dipoles to their SCRF polarized values using
c     a neighbor list
c
c
      subroutine ediff1b
      use sizes
      use atoms
      use bound
      use boxes
      use chgpot
      use couple
      use deriv
      use energi
      use group
      use inter
      use limits
      use molcul
      use mplpot
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      integer ix,iy,iz
      integer kx,ky,kz
      real*8 ei,f,fgrp
      real*8 damp,gfd
      real*8 scale3,scale5
      real*8 scale7,scale9
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 psc3,psc5,psc7,psc9
      real*8 dsc3,dsc5,dsc7,dsc9
      real*8 scale3i,scale5i
      real*8 scale7i
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9
      real*8 pdi,pti,pgamma
      real*8 ci,di(3),qi(9)
      real*8 ck,dk(3),qk(9)
      real*8 fridmp(3),findmp(3)
      real*8 ftm2i(3)
      real*8 ttm2i(3),ttm3i(3)
      real*8 dixdk(3),fdir(3)
      real*8 dixuk(3),dkxui(3)
      real*8 dixukp(3),dkxuip(3)
      real*8 uixqkr(3),ukxqir(3)
      real*8 uixqkrp(3),ukxqirp(3)
      real*8 qiuk(3),qkui(3)
      real*8 qiukp(3),qkuip(3)
      real*8 rxqiuk(3),rxqkui(3)
      real*8 rxqiukp(3),rxqkuip(3)
      real*8 qidk(3),qkdi(3)
      real*8 qir(3),qkr(3)
      real*8 qiqkr(3),qkqir(3)
      real*8 qixqk(3),rxqir(3)
      real*8 dixr(3),dkxr(3)
      real*8 dixqkr(3),dkxqir(3)
      real*8 rxqkr(3),qkrxqir(3)
      real*8 rxqikr(3),rxqkir(3)
      real*8 rxqidk(3),rxqkdi(3)
      real*8 ddsc3(3),ddsc5(3)
      real*8 ddsc7(3)
      real*8 gli(7),glip(7)
      real*8 sc(10)
      real*8 sci(8),scip(8)
      real*8 gfi(6),gti(6)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: uscale(:)
      real*8, allocatable :: trqi(:,:)
      logical proceed,usei,usek
      character*6 mode
c
c
c     set conversion factor, cutoff and scaling coefficients
c
      if (npole .eq. 0)  return
      f = electric / dielec
      mode = 'MPOLE'
      call switch (mode)
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
      allocate (dscale(n))
      allocate (uscale(n))
      allocate (trqi(3,n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
         uscale(i) = 1.0d0
      end do
c
c     initialize local variables for OpenMP calculation
c
      do i = 1, n
         do j = 1, 3
            trqi(j,i) = 0.0d0
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(npole,ipole,x,y,z,xaxis,yaxis,
!$OMP& zaxis,pdamp,thole,rpole,uind,uinp,uinds,uinps,nelst,elst,
!$OMP% use,n12,n13,n14,n15,i12,i13,i14,i15,np11,ip11,np12,ip12,np13,
!$OMP% ip13,np14,ip14,p2scale,p3scale,p4scale,p41scale,p5scale,
!$OMP% d1scale,d2scale,d3scale,d4scale,u1scale,u2scale,u3scale,u4scale,
!$OMP% use_group,use_intra,off2,f)
!$OMP& firstprivate(pscale,dscale,uscale)
!$OMP% shared(es,des,trqi)
!$OMP DO reduction(+:es,des,trqi) schedule(guided)
c
c     calculate the multipole interaction energy and gradient
c
      do i = 1, npole-1
         ii = ipole(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = yaxis(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         di(1) = rpole(2,i)
         di(2) = rpole(3,i)
         di(3) = rpole(4,i)
         qi(1) = rpole(5,i)
         qi(2) = rpole(6,i)
         qi(3) = rpole(7,i)
         qi(4) = rpole(8,i)
         qi(5) = rpole(9,i)
         qi(6) = rpole(10,i)
         qi(7) = rpole(11,i)
         qi(8) = rpole(12,i)
         qi(9) = rpole(13,i)
         usei = (use(ii) .or. use(iz) .or. use(ix) .or. use(iy))
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
                if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
            uscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
            uscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
            uscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
            uscale(ip14(j,ii)) = u4scale
         end do
         do kkk = 1, nelst(i)
            k = elst(kkk,i)
            kk = ipole(k)
            kz = zaxis(k)
            kx = xaxis(k)
            ky = yaxis(k)
            usek = (use(kk) .or. use(kz) .or. use(kx) .or. use(ky))
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (.not. use_intra)  proceed = .true.
            if (proceed)  proceed = (usei .or. usek)
            if (.not. proceed)  goto 10
            ck = rpole(1,k)
            dk(1) = rpole(2,k)
            dk(2) = rpole(3,k)
            dk(3) = rpole(4,k)
            qk(1) = rpole(5,k)
            qk(2) = rpole(6,k)
            qk(3) = rpole(7,k)
            qk(4) = rpole(8,k)
            qk(5) = rpole(9,k)
            qk(6) = rpole(10,k)
            qk(7) = rpole(11,k)
            qk(8) = rpole(12,k)
            qk(9) = rpole(13,k)
            xr = x(kk) - xi
            yr = y(kk) - yi
            zr = z(kk) - zi
            call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
               rr1 = 1.0d0 / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               scale3 = 1.0d0
               scale5 = 1.0d0
               scale7 = 1.0d0
               scale9 = 1.0d0
               do j = 1, 3
                  ddsc3(j) = 0.0d0
                  ddsc5(j) = 0.0d0
                  ddsc7(j) = 0.0d0
               end do
c
c     apply Thole polarization damping to scale factors
c
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     scale3 = 1.0d0 - exp(damp)
                     scale5 = 1.0d0 - (1.0d0-damp)*exp(damp)
                     scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                       *exp(damp)
                     scale9 = 1.0d0 - (1.0d0-damp+(18.0d0*damp**2
     &                                 -9.0d0*damp**3)/35.0d0)*exp(damp)
                     ddsc3(1) = -3.0d0*damp*exp(damp) * xr/r2
                     ddsc3(2) = -3.0d0*damp*exp(damp) * yr/r2
                     ddsc3(3) = -3.0d0*damp*exp(damp) * zr/r2
                     ddsc5(1) = -damp * ddsc3(1)
                     ddsc5(2) = -damp * ddsc3(2)
                     ddsc5(3) = -damp * ddsc3(3)
                     ddsc7(1) = (-0.2d0-0.6d0*damp) * ddsc5(1)
                     ddsc7(2) = (-0.2d0-0.6d0*damp) * ddsc5(2)
                     ddsc7(3) = (-0.2d0-0.6d0*damp) * ddsc5(3)
                  end if
               end if
               scale3i = scale3 * uscale(kk)
               scale5i = scale5 * uscale(kk)
               scale7i = scale7 * uscale(kk)
               dsc3 = scale3 * dscale(kk)
               dsc5 = scale5 * dscale(kk)
               dsc7 = scale7 * dscale(kk)
               dsc9 = scale9 * dscale(kk)
               psc3 = scale3 * pscale(kk)
               psc5 = scale5 * pscale(kk)
               psc7 = scale7 * pscale(kk)
               psc9 = scale9 * pscale(kk)
c
c     construct auxiliary vectors for permanent terms
c
               dixdk(1) = di(2)*dk(3) - di(3)*dk(2)
               dixdk(2) = di(3)*dk(1) - di(1)*dk(3)
               dixdk(3) = di(1)*dk(2) - di(2)*dk(1)
               dixr(1) = di(2)*zr - di(3)*yr
               dixr(2) = di(3)*xr - di(1)*zr
               dixr(3) = di(1)*yr - di(2)*xr
               dkxr(1) = dk(2)*zr - dk(3)*yr
               dkxr(2) = dk(3)*xr - dk(1)*zr
               dkxr(3) = dk(1)*yr - dk(2)*xr
               qir(1) = qi(1)*xr + qi(4)*yr + qi(7)*zr
               qir(2) = qi(2)*xr + qi(5)*yr + qi(8)*zr
               qir(3) = qi(3)*xr + qi(6)*yr + qi(9)*zr
               qkr(1) = qk(1)*xr + qk(4)*yr + qk(7)*zr
               qkr(2) = qk(2)*xr + qk(5)*yr + qk(8)*zr
               qkr(3) = qk(3)*xr + qk(6)*yr + qk(9)*zr
               qiqkr(1) = qi(1)*qkr(1) + qi(4)*qkr(2) + qi(7)*qkr(3)
               qiqkr(2) = qi(2)*qkr(1) + qi(5)*qkr(2) + qi(8)*qkr(3)
               qiqkr(3) = qi(3)*qkr(1) + qi(6)*qkr(2) + qi(9)*qkr(3)
               qkqir(1) = qk(1)*qir(1) + qk(4)*qir(2) + qk(7)*qir(3)
               qkqir(2) = qk(2)*qir(1) + qk(5)*qir(2) + qk(8)*qir(3)
               qkqir(3) = qk(3)*qir(1) + qk(6)*qir(2) + qk(9)*qir(3)
               qixqk(1) = qi(2)*qk(3) + qi(5)*qk(6) + qi(8)*qk(9)
     &                       - qi(3)*qk(2) - qi(6)*qk(5) - qi(9)*qk(8)
               qixqk(2) = qi(3)*qk(1) + qi(6)*qk(4) + qi(9)*qk(7)
     &                       - qi(1)*qk(3) - qi(4)*qk(6) - qi(7)*qk(9)
               qixqk(3) = qi(1)*qk(2) + qi(4)*qk(5) + qi(7)*qk(8)
     &                       - qi(2)*qk(1) - qi(5)*qk(4) - qi(8)*qk(7)
               rxqir(1) = yr*qir(3) - zr*qir(2)
               rxqir(2) = zr*qir(1) - xr*qir(3)
               rxqir(3) = xr*qir(2) - yr*qir(1)
               rxqkr(1) = yr*qkr(3) - zr*qkr(2)
               rxqkr(2) = zr*qkr(1) - xr*qkr(3)
               rxqkr(3) = xr*qkr(2) - yr*qkr(1)
               rxqikr(1) = yr*qiqkr(3) - zr*qiqkr(2)
               rxqikr(2) = zr*qiqkr(1) - xr*qiqkr(3)
               rxqikr(3) = xr*qiqkr(2) - yr*qiqkr(1)
               rxqkir(1) = yr*qkqir(3) - zr*qkqir(2)
               rxqkir(2) = zr*qkqir(1) - xr*qkqir(3)
               rxqkir(3) = xr*qkqir(2) - yr*qkqir(1)
               qkrxqir(1) = qkr(2)*qir(3) - qkr(3)*qir(2)
               qkrxqir(2) = qkr(3)*qir(1) - qkr(1)*qir(3)
               qkrxqir(3) = qkr(1)*qir(2) - qkr(2)*qir(1)
               qidk(1) = qi(1)*dk(1) + qi(4)*dk(2) + qi(7)*dk(3)
               qidk(2) = qi(2)*dk(1) + qi(5)*dk(2) + qi(8)*dk(3)
               qidk(3) = qi(3)*dk(1) + qi(6)*dk(2) + qi(9)*dk(3)
               qkdi(1) = qk(1)*di(1) + qk(4)*di(2) + qk(7)*di(3)
               qkdi(2) = qk(2)*di(1) + qk(5)*di(2) + qk(8)*di(3)
               qkdi(3) = qk(3)*di(1) + qk(6)*di(2) + qk(9)*di(3)
               dixqkr(1) = di(2)*qkr(3) - di(3)*qkr(2)
               dixqkr(2) = di(3)*qkr(1) - di(1)*qkr(3)
               dixqkr(3) = di(1)*qkr(2) - di(2)*qkr(1)
               dkxqir(1) = dk(2)*qir(3) - dk(3)*qir(2)
               dkxqir(2) = dk(3)*qir(1) - dk(1)*qir(3)
               dkxqir(3) = dk(1)*qir(2) - dk(2)*qir(1)
               rxqidk(1) = yr*qidk(3) - zr*qidk(2)
               rxqidk(2) = zr*qidk(1) - xr*qidk(3)
               rxqidk(3) = xr*qidk(2) - yr*qidk(1)
               rxqkdi(1) = yr*qkdi(3) - zr*qkdi(2)
               rxqkdi(2) = zr*qkdi(1) - xr*qkdi(3)
               rxqkdi(3) = xr*qkdi(2) - yr*qkdi(1)
c
c     get intermediate variables for permanent energy terms
c
               sc(3) = di(1)*xr + di(2)*yr + di(3)*zr
               sc(4) = dk(1)*xr + dk(2)*yr + dk(3)*zr
               sc(5) = qir(1)*xr + qir(2)*yr + qir(3)*zr
               sc(6) = qkr(1)*xr + qkr(2)*yr + qkr(3)*zr
c
c     construct auxiliary vectors for induced terms
c
               dixuk(1) = di(2)*uinds(3,k) - di(3)*uinds(2,k)
               dixuk(2) = di(3)*uinds(1,k) - di(1)*uinds(3,k)
               dixuk(3) = di(1)*uinds(2,k) - di(2)*uinds(1,k)
               dkxui(1) = dk(2)*uinds(3,i) - dk(3)*uinds(2,i)
               dkxui(2) = dk(3)*uinds(1,i) - dk(1)*uinds(3,i)
               dkxui(3) = dk(1)*uinds(2,i) - dk(2)*uinds(1,i)
               dixukp(1) = di(2)*uinps(3,k) - di(3)*uinps(2,k)
               dixukp(2) = di(3)*uinps(1,k) - di(1)*uinps(3,k)
               dixukp(3) = di(1)*uinps(2,k) - di(2)*uinps(1,k)
               dkxuip(1) = dk(2)*uinps(3,i) - dk(3)*uinps(2,i)
               dkxuip(2) = dk(3)*uinps(1,i) - dk(1)*uinps(3,i)
               dkxuip(3) = dk(1)*uinps(2,i) - dk(2)*uinps(1,i)
               qiuk(1) = qi(1)*uinds(1,k) + qi(4)*uinds(2,k)
     &                      + qi(7)*uinds(3,k)
               qiuk(2) = qi(2)*uinds(1,k) + qi(5)*uinds(2,k)
     &                      + qi(8)*uinds(3,k)
               qiuk(3) = qi(3)*uinds(1,k) + qi(6)*uinds(2,k)
     &                      + qi(9)*uinds(3,k)
               qkui(1) = qk(1)*uinds(1,i) + qk(4)*uinds(2,i)
     &                      + qk(7)*uinds(3,i)
               qkui(2) = qk(2)*uinds(1,i) + qk(5)*uinds(2,i)
     &                      + qk(8)*uinds(3,i)
               qkui(3) = qk(3)*uinds(1,i) + qk(6)*uinds(2,i)
     &                      + qk(9)*uinds(3,i)
               qiukp(1) = qi(1)*uinps(1,k) + qi(4)*uinps(2,k)
     &                       + qi(7)*uinps(3,k)
               qiukp(2) = qi(2)*uinps(1,k) + qi(5)*uinps(2,k)
     &                       + qi(8)*uinps(3,k)
               qiukp(3) = qi(3)*uinps(1,k) + qi(6)*uinps(2,k)
     &                       + qi(9)*uinps(3,k)
               qkuip(1) = qk(1)*uinps(1,i) + qk(4)*uinps(2,i)
     &                       + qk(7)*uinps(3,i)
               qkuip(2) = qk(2)*uinps(1,i) + qk(5)*uinps(2,i)
     &                       + qk(8)*uinps(3,i)
               qkuip(3) = qk(3)*uinps(1,i) + qk(6)*uinps(2,i)
     &                       + qk(9)*uinps(3,i)
               uixqkr(1) = uinds(2,i)*qkr(3) - uinds(3,i)*qkr(2)
               uixqkr(2) = uinds(3,i)*qkr(1) - uinds(1,i)*qkr(3)
               uixqkr(3) = uinds(1,i)*qkr(2) - uinds(2,i)*qkr(1)
               ukxqir(1) = uinds(2,k)*qir(3) - uinds(3,k)*qir(2)
               ukxqir(2) = uinds(3,k)*qir(1) - uinds(1,k)*qir(3)
               ukxqir(3) = uinds(1,k)*qir(2) - uinds(2,k)*qir(1)
               uixqkrp(1) = uinps(2,i)*qkr(3) - uinps(3,i)*qkr(2)
               uixqkrp(2) = uinps(3,i)*qkr(1) - uinps(1,i)*qkr(3)
               uixqkrp(3) = uinps(1,i)*qkr(2) - uinps(2,i)*qkr(1)
               ukxqirp(1) = uinps(2,k)*qir(3) - uinps(3,k)*qir(2)
               ukxqirp(2) = uinps(3,k)*qir(1) - uinps(1,k)*qir(3)
               ukxqirp(3) = uinps(1,k)*qir(2) - uinps(2,k)*qir(1)
               rxqiuk(1) = yr*qiuk(3) - zr*qiuk(2)
               rxqiuk(2) = zr*qiuk(1) - xr*qiuk(3)
               rxqiuk(3) = xr*qiuk(2) - yr*qiuk(1)
               rxqkui(1) = yr*qkui(3) - zr*qkui(2)
               rxqkui(2) = zr*qkui(1) - xr*qkui(3)
               rxqkui(3) = xr*qkui(2) - yr*qkui(1)
               rxqiukp(1) = yr*qiukp(3) - zr*qiukp(2)
               rxqiukp(2) = zr*qiukp(1) - xr*qiukp(3)
               rxqiukp(3) = xr*qiukp(2) - yr*qiukp(1)
               rxqkuip(1) = yr*qkuip(3) - zr*qkuip(2)
               rxqkuip(2) = zr*qkuip(1) - xr*qkuip(3)
               rxqkuip(3) = xr*qkuip(2) - yr*qkuip(1)
c
c     get intermediate variables for induction energy terms
c
               sci(1) = uinds(1,i)*dk(1) + uinds(2,i)*dk(2)
     &                     + uinds(3,i)*dk(3) + di(1)*uinds(1,k)
     &                     + di(2)*uinds(2,k) + di(3)*uinds(3,k)
               sci(2) = uinds(1,i)*uinds(1,k) + uinds(2,i)*uinds(2,k)
     &                     + uinds(3,i)*uinds(3,k)
               sci(3) = uinds(1,i)*xr + uinds(2,i)*yr + uinds(3,i)*zr
               sci(4) = uinds(1,k)*xr + uinds(2,k)*yr + uinds(3,k)*zr
               sci(7) = qir(1)*uinds(1,k) + qir(2)*uinds(2,k)
     &                     + qir(3)*uinds(3,k)
               sci(8) = qkr(1)*uinds(1,i) + qkr(2)*uinds(2,i)
     &                     + qkr(3)*uinds(3,i)
               scip(1) = uinps(1,i)*dk(1) + uinps(2,i)*dk(2)
     &                      + uinps(3,i)*dk(3) + di(1)*uinps(1,k)
     &                      + di(2)*uinps(2,k) + di(3)*uinps(3,k)
               scip(2) = uinds(1,i)*uinps(1,k) + uinds(2,i)*uinps(2,k)
     &                 + uinds(3,i)*uinps(3,k) + uinps(1,i)*uinds(1,k)
     &                 + uinps(2,i)*uinds(2,k) + uinps(3,i)*uinds(3,k)
               scip(3) = uinps(1,i)*xr + uinps(2,i)*yr + uinps(3,i)*zr
               scip(4) = uinps(1,k)*xr + uinps(2,k)*yr + uinps(3,k)*zr
               scip(7) = qir(1)*uinps(1,k) + qir(2)*uinps(2,k)
     &                      + qir(3)*uinps(3,k)
               scip(8) = qkr(1)*uinps(1,i) + qkr(2)*uinps(2,i)
     &                      + qkr(3)*uinps(3,i)
c
c     calculate the gl functions for potential energy
c
               gli(1) = ck*sci(3) - ci*sci(4)
               gli(2) = -sc(3)*sci(4) - sci(3)*sc(4)
               gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
               gli(6) = sci(1)
               gli(7) = 2.0d0 * (sci(7)-sci(8))
               glip(1) = ck*scip(3) - ci*scip(4)
               glip(2) = -sc(3)*scip(4) - scip(3)*sc(4)
               glip(3) = scip(3)*sc(6) - scip(4)*sc(5)
               glip(6) = scip(1)
               glip(7) = 2.0d0 * (scip(7)-scip(8))
c
c     get the permanent multipole and induced energies
c
               ei = 0.5d0 * (rr3*(gli(1)+gli(6))*psc3
     &                          + rr5*(gli(2)+gli(7))*psc5
     &                          + rr7*gli(3)*psc7)
               ei = f * ei
               es = es + ei
c
c     intermediate variables for the induced-permanent terms
c
               gfi(1) = 0.5d0*rr5*((gli(1)+gli(6))*psc3
     &                     +(glip(1)+glip(6))*dsc3+scip(2)*scale3i)
     &                + 0.5d0*rr7*((gli(7)+gli(2))*psc5
     &                            +(glip(7)+glip(2))*dsc5
     &                     -(sci(3)*scip(4)+scip(3)*sci(4))*scale5i)
     &                + 0.5d0*rr9*(gli(3)*psc7+glip(3)*dsc7)
               gfi(2) = -rr3*ck + rr5*sc(4) - rr7*sc(6)
               gfi(3) = rr3*ci + rr5*sc(3) + rr7*sc(5)
               gfi(4) = 2.0d0 * rr5
               gfi(5) = rr7 * (sci(4)*psc7+scip(4)*dsc7)
               gfi(6) = -rr7 * (sci(3)*psc7+scip(3)*dsc7)
c
c     get the induced force
c
               ftm2i(1) = gfi(1)*xr + 0.5d0*
     &            (- rr3*ck*(uinds(1,i)*psc3+uinps(1,i)*dsc3)
     &             + rr5*sc(4)*(uinds(1,i)*psc5+uinps(1,i)*dsc5)
     &             - rr7*sc(6)*(uinds(1,i)*psc7+uinps(1,i)*dsc7))
     &             +(rr3*ci*(uinds(1,k)*psc3+uinps(1,k)*dsc3)
     &             + rr5*sc(3)*(uinds(1,k)*psc5+uinps(1,k)*dsc5)
     &             + rr7*sc(5)*(uinds(1,k)*psc7+uinps(1,k)*dsc7))*0.5d0
     &             + rr5*scale5i*(sci(4)*uinps(1,i)+scip(4)*uinds(1,i)
     &             + sci(3)*uinps(1,k)+scip(3)*uinds(1,k))*0.5d0
     &             + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(1)
     &             + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(1)
     &             + 0.5d0*gfi(4)*((qkui(1)-qiuk(1))*psc5
     &             + (qkuip(1)-qiukp(1))*dsc5)
     &             + gfi(5)*qir(1) + gfi(6)*qkr(1)
               ftm2i(2) = gfi(1)*yr + 0.5d0*
     &            (- rr3*ck*(uinds(2,i)*psc3+uinps(2,i)*dsc3)
     &             + rr5*sc(4)*(uinds(2,i)*psc5+uinps(2,i)*dsc5)
     &             - rr7*sc(6)*(uinds(2,i)*psc7+uinps(2,i)*dsc7))
     &             +(rr3*ci*(uinds(2,k)*psc3+uinps(2,k)*dsc3)
     &             + rr5*sc(3)*(uinds(2,k)*psc5+uinps(2,k)*dsc5)
     &             + rr7*sc(5)*(uinds(2,k)*psc7+uinps(2,k)*dsc7))*0.5d0
     &             + rr5*scale5i*(sci(4)*uinps(2,i)+scip(4)*uinds(2,i)
     &             + sci(3)*uinps(2,k)+scip(3)*uinds(2,k))*0.5d0
     &             + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(2)
     &             + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(2)
     &             + 0.5d0*gfi(4)*((qkui(2)-qiuk(2))*psc5
     &             + (qkuip(2)-qiukp(2))*dsc5)
     &             + gfi(5)*qir(2) + gfi(6)*qkr(2)
               ftm2i(3) = gfi(1)*zr  + 0.5d0*
     &            (- rr3*ck*(uinds(3,i)*psc3+uinps(3,i)*dsc3)
     &             + rr5*sc(4)*(uinds(3,i)*psc5+uinps(3,i)*dsc5)
     &             - rr7*sc(6)*(uinds(3,i)*psc7+uinps(3,i)*dsc7))
     &             +(rr3*ci*(uinds(3,k)*psc3+uinps(3,k)*dsc3)
     &             + rr5*sc(3)*(uinds(3,k)*psc5+uinps(3,k)*dsc5)
     &             + rr7*sc(5)*(uinds(3,k)*psc7+uinps(3,k)*dsc7))*0.5d0
     &             + rr5*scale5i*(sci(4)*uinps(3,i)+scip(4)*uinds(3,i)
     &             + sci(3)*uinps(3,k)+scip(3)*uinds(3,k))*0.5d0
     &             + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(3)
     &             + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(3)
     &             + 0.5d0*gfi(4)*((qkui(3)-qiuk(3))*psc5
     &             + (qkuip(3)-qiukp(3))*dsc5)
     &             + gfi(5)*qir(3) + gfi(6)*qkr(3)
c
c     intermediate values needed for partially excluded interactions
c
               fridmp(1) = 0.5d0 * (rr3*((gli(1)+gli(6))*pscale(kk)
     &                        +(glip(1)+glip(6))*dscale(kk))*ddsc3(1)
     &            + rr5*((gli(2)+gli(7))*pscale(kk)
     &                +(glip(2)+glip(7))*dscale(kk))*ddsc5(1)
     &            + rr7*(gli(3)*pscale(kk)+glip(3)*dscale(kk))*ddsc7(1))
               fridmp(2) = 0.5d0 * (rr3*((gli(1)+gli(6))*pscale(kk)
     &                        +(glip(1)+glip(6))*dscale(kk))*ddsc3(2)
     &            + rr5*((gli(2)+gli(7))*pscale(kk)
     &                +(glip(2)+glip(7))*dscale(kk))*ddsc5(2)
     &            + rr7*(gli(3)*pscale(kk)+glip(3)*dscale(kk))*ddsc7(2))
               fridmp(3) = 0.5d0 * (rr3*((gli(1)+gli(6))*pscale(kk)
     &                        +(glip(1)+glip(6))*dscale(kk))*ddsc3(3)
     &            + rr5*((gli(2)+gli(7))*pscale(kk)
     &                +(glip(2)+glip(7))*dscale(kk))*ddsc5(3)
     &            + rr7*(gli(3)*pscale(kk)+glip(3)*dscale(kk))*ddsc7(3))
c
c     get the induced-induced derivative terms
c
               findmp(1) = 0.5d0 * uscale(kk) * (scip(2)*rr3*ddsc3(1)
     &                   - rr5*ddsc5(1)*(sci(3)*scip(4)+scip(3)*sci(4)))
               findmp(2) = 0.5d0 * uscale(kk) * (scip(2)*rr3*ddsc3(2)
     &                   - rr5*ddsc5(2)*(sci(3)*scip(4)+scip(3)*sci(4)))
               findmp(3) = 0.5d0 * uscale(kk) * (scip(2)*rr3*ddsc3(3)
     &                   - rr5*ddsc5(3)*(sci(3)*scip(4)+scip(3)*sci(4)))
c
c     handle of scaling for partially excluded interactions
c
               ftm2i(1) = ftm2i(1) - fridmp(1) - findmp(1)
               ftm2i(2) = ftm2i(2) - fridmp(2) - findmp(2)
               ftm2i(3) = ftm2i(3) - fridmp(3) - findmp(3)
c
c     correction to convert mutual to direct polarization force
c
               if (poltyp .eq. 'DIRECT') then
                  gfd = 0.5d0 * (rr5*scip(2)*scale3i
     &                  - rr7*(scip(3)*sci(4)+sci(3)*scip(4))*scale5i)
                  fdir(1) = gfd*xr + 0.5d0*rr5*scale5i
     &                         * (sci(4)*uinps(1,i)+scip(4)*uinds(1,i)
     &                           +sci(3)*uinps(1,k)+scip(3)*uinds(1,k))
                  fdir(2) = gfd*yr + 0.5d0*rr5*scale5i
     &                         * (sci(4)*uinps(2,i)+scip(4)*uinds(2,i)
     &                           +sci(3)*uinps(2,k)+scip(3)*uinds(2,k))
                  fdir(3) = gfd*zr + 0.5d0*rr5*scale5i
     &                         * (sci(4)*uinps(3,i)+scip(4)*uinds(3,i)
     &                           +sci(3)*uinps(3,k)+scip(3)*uinds(3,k))
                  ftm2i(1) = ftm2i(1) - fdir(1) + findmp(1)
                  ftm2i(2) = ftm2i(2) - fdir(2) + findmp(2)
                  ftm2i(3) = ftm2i(3) - fdir(3) + findmp(3)
               end if
c
c     now perform the torque calculation
c     intermediate terms for torque between multipoles i and k
c
               gti(2) = 0.5d0 * (sci(4)*psc5+scip(4)*dsc5) * rr5
               gti(3) = 0.5d0 * (sci(3)*psc5+scip(3)*dsc5) * rr5
               gti(4) = gfi(4)
               gti(5) = gfi(5)
               gti(6) = gfi(6)
c
c     calculate the induced torque components
c
               ttm2i(1) = -rr3*(dixuk(1)*psc3+dixukp(1)*dsc3)*0.5d0
     &            + gti(2)*dixr(1) + gti(4)*((ukxqir(1)+rxqiuk(1))*psc5
     &            +(ukxqirp(1)+rxqiukp(1))*dsc5)*0.5d0 - gti(5)*rxqir(1)
               ttm2i(2) = -rr3*(dixuk(2)*psc3+dixukp(2)*dsc3)*0.5d0
     &            + gti(2)*dixr(2) + gti(4)*((ukxqir(2)+rxqiuk(2))*psc5
     &            +(ukxqirp(2)+rxqiukp(2))*dsc5)*0.5d0 - gti(5)*rxqir(2)
               ttm2i(3) = -rr3*(dixuk(3)*psc3+dixukp(3)*dsc3)*0.5d0
     &            + gti(2)*dixr(3) + gti(4)*((ukxqir(3)+rxqiuk(3))*psc5
     &            +(ukxqirp(3)+rxqiukp(3))*dsc5)*0.5d0 - gti(5)*rxqir(3)
               ttm3i(1) = -rr3*(dkxui(1)*psc3+dkxuip(1)*dsc3)*0.5d0
     &            + gti(3)*dkxr(1) - gti(4)*((uixqkr(1)+rxqkui(1))*psc5
     &            +(uixqkrp(1)+rxqkuip(1))*dsc5)*0.5d0 - gti(6)*rxqkr(1)
               ttm3i(2) = -rr3*(dkxui(2)*psc3+dkxuip(2)*dsc3)*0.5d0
     &            + gti(3)*dkxr(2) - gti(4)*((uixqkr(2)+rxqkui(2))*psc5
     &            +(uixqkrp(2)+rxqkuip(2))*dsc5)*0.5d0 - gti(6)*rxqkr(2)
               ttm3i(3) = -rr3*(dkxui(3)*psc3+dkxuip(3)*dsc3)*0.5d0
     &            + gti(3)*dkxr(3) - gti(4)*((uixqkr(3)+rxqkui(3))*psc5
     &            +(uixqkrp(3)+rxqkuip(3))*dsc5)*0.5d0 - gti(6)*rxqkr(3)
c
c     update the force components on sites i and k
c
               des(1,ii) = des(1,ii) + f*ftm2i(1)
               des(2,ii) = des(2,ii) + f*ftm2i(2)
               des(3,ii) = des(3,ii) + f*ftm2i(3)
               des(1,kk) = des(1,kk) - f*ftm2i(1)
               des(2,kk) = des(2,kk) - f*ftm2i(2)
               des(3,kk) = des(3,kk) - f*ftm2i(3)
c
c     update the torque components on sites i and k
c
               trqi(1,ii) = trqi(1,ii) + f*ttm2i(1)
               trqi(2,ii) = trqi(2,ii) + f*ttm2i(2)
               trqi(3,ii) = trqi(3,ii) + f*ttm2i(3)
               trqi(1,kk) = trqi(1,kk) + f*ttm3i(1)
               trqi(2,kk) = trqi(2,kk) + f*ttm3i(2)
               trqi(3,kk) = trqi(3,kk) + f*ttm3i(3)
c
c     construct auxiliary vectors for induced terms
c
               dixuk(1) = di(2)*uind(3,k) - di(3)*uind(2,k)
               dixuk(2) = di(3)*uind(1,k) - di(1)*uind(3,k)
               dixuk(3) = di(1)*uind(2,k) - di(2)*uind(1,k)
               dkxui(1) = dk(2)*uind(3,i) - dk(3)*uind(2,i)
               dkxui(2) = dk(3)*uind(1,i) - dk(1)*uind(3,i)
               dkxui(3) = dk(1)*uind(2,i) - dk(2)*uind(1,i)
               dixukp(1) = di(2)*uinp(3,k) - di(3)*uinp(2,k)
               dixukp(2) = di(3)*uinp(1,k) - di(1)*uinp(3,k)
               dixukp(3) = di(1)*uinp(2,k) - di(2)*uinp(1,k)
               dkxuip(1) = dk(2)*uinp(3,i) - dk(3)*uinp(2,i)
               dkxuip(2) = dk(3)*uinp(1,i) - dk(1)*uinp(3,i)
               dkxuip(3) = dk(1)*uinp(2,i) - dk(2)*uinp(1,i)
               qiuk(1) = qi(1)*uind(1,k) + qi(4)*uind(2,k)
     &                      + qi(7)*uind(3,k)
               qiuk(2) = qi(2)*uind(1,k) + qi(5)*uind(2,k)
     &                      + qi(8)*uind(3,k)
               qiuk(3) = qi(3)*uind(1,k) + qi(6)*uind(2,k)
     &                      + qi(9)*uind(3,k)
               qkui(1) = qk(1)*uind(1,i) + qk(4)*uind(2,i)
     &                      + qk(7)*uind(3,i)
               qkui(2) = qk(2)*uind(1,i) + qk(5)*uind(2,i)
     &                      + qk(8)*uind(3,i)
               qkui(3) = qk(3)*uind(1,i) + qk(6)*uind(2,i)
     &                      + qk(9)*uind(3,i)
               qiukp(1) = qi(1)*uinp(1,k) + qi(4)*uinp(2,k)
     &                       + qi(7)*uinp(3,k)
               qiukp(2) = qi(2)*uinp(1,k) + qi(5)*uinp(2,k)
     &                       + qi(8)*uinp(3,k)
               qiukp(3) = qi(3)*uinp(1,k) + qi(6)*uinp(2,k)
     &                       + qi(9)*uinp(3,k)
               qkuip(1) = qk(1)*uinp(1,i) + qk(4)*uinp(2,i)
     &                       + qk(7)*uinp(3,i)
               qkuip(2) = qk(2)*uinp(1,i) + qk(5)*uinp(2,i)
     &                       + qk(8)*uinp(3,i)
               qkuip(3) = qk(3)*uinp(1,i) + qk(6)*uinp(2,i)
     &                       + qk(9)*uinp(3,i)
               uixqkr(1) = uind(2,i)*qkr(3) - uind(3,i)*qkr(2)
               uixqkr(2) = uind(3,i)*qkr(1) - uind(1,i)*qkr(3)
               uixqkr(3) = uind(1,i)*qkr(2) - uind(2,i)*qkr(1)
               ukxqir(1) = uind(2,k)*qir(3) - uind(3,k)*qir(2)
               ukxqir(2) = uind(3,k)*qir(1) - uind(1,k)*qir(3)
               ukxqir(3) = uind(1,k)*qir(2) - uind(2,k)*qir(1)
               uixqkrp(1) = uinp(2,i)*qkr(3) - uinp(3,i)*qkr(2)
               uixqkrp(2) = uinp(3,i)*qkr(1) - uinp(1,i)*qkr(3)
               uixqkrp(3) = uinp(1,i)*qkr(2) - uinp(2,i)*qkr(1)
               ukxqirp(1) = uinp(2,k)*qir(3) - uinp(3,k)*qir(2)
               ukxqirp(2) = uinp(3,k)*qir(1) - uinp(1,k)*qir(3)
               ukxqirp(3) = uinp(1,k)*qir(2) - uinp(2,k)*qir(1)
               rxqiuk(1) = yr*qiuk(3) - zr*qiuk(2)
               rxqiuk(2) = zr*qiuk(1) - xr*qiuk(3)
               rxqiuk(3) = xr*qiuk(2) - yr*qiuk(1)
               rxqkui(1) = yr*qkui(3) - zr*qkui(2)
               rxqkui(2) = zr*qkui(1) - xr*qkui(3)
               rxqkui(3) = xr*qkui(2) - yr*qkui(1)
               rxqiukp(1) = yr*qiukp(3) - zr*qiukp(2)
               rxqiukp(2) = zr*qiukp(1) - xr*qiukp(3)
               rxqiukp(3) = xr*qiukp(2) - yr*qiukp(1)
               rxqkuip(1) = yr*qkuip(3) - zr*qkuip(2)
               rxqkuip(2) = zr*qkuip(1) - xr*qkuip(3)
               rxqkuip(3) = xr*qkuip(2) - yr*qkuip(1)
c
c     get intermediate variables for induction energy terms
c
               sci(1) = uind(1,i)*dk(1) + uind(2,i)*dk(2)
     &                     + uind(3,i)*dk(3) + di(1)*uind(1,k)
     &                     + di(2)*uind(2,k) + di(3)*uind(3,k)
               sci(2) = uind(1,i)*uind(1,k) + uind(2,i)*uind(2,k)
     &                     + uind(3,i)*uind(3,k)
               sci(3) = uind(1,i)*xr + uind(2,i)*yr + uind(3,i)*zr
               sci(4) = uind(1,k)*xr + uind(2,k)*yr + uind(3,k)*zr
               sci(7) = qir(1)*uind(1,k) + qir(2)*uind(2,k)
     &                     + qir(3)*uind(3,k)
               sci(8) = qkr(1)*uind(1,i) + qkr(2)*uind(2,i)
     &                     + qkr(3)*uind(3,i)
               scip(1) = uinp(1,i)*dk(1) + uinp(2,i)*dk(2)
     &                      + uinp(3,i)*dk(3) + di(1)*uinp(1,k)
     &                      + di(2)*uinp(2,k) + di(3)*uinp(3,k)
               scip(2) = uind(1,i)*uinp(1,k)+uind(2,i)*uinp(2,k)
     &                      + uind(3,i)*uinp(3,k)+uinp(1,i)*uind(1,k)
     &                      + uinp(2,i)*uind(2,k)+uinp(3,i)*uind(3,k)
               scip(3) = uinp(1,i)*xr + uinp(2,i)*yr + uinp(3,i)*zr
               scip(4) = uinp(1,k)*xr + uinp(2,k)*yr + uinp(3,k)*zr
               scip(7) = qir(1)*uinp(1,k) + qir(2)*uinp(2,k)
     &                      + qir(3)*uinp(3,k)
               scip(8) = qkr(1)*uinp(1,i) + qkr(2)*uinp(2,i)
     &                      + qkr(3)*uinp(3,i)
c
c     calculate the gl functions for potential energy
c
               gli(1) = ck*sci(3) - ci*sci(4)
               gli(2) = -sc(3)*sci(4) - sci(3)*sc(4)
               gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
               gli(6) = sci(1)
               gli(7) = 2.0d0 * (sci(7)-sci(8))
               glip(1) = ck*scip(3) - ci*scip(4)
               glip(2) = -sc(3)*scip(4) - scip(3)*sc(4)
               glip(3) = scip(3)*sc(6) - scip(4)*sc(5)
               glip(6) = scip(1)
               glip(7) = 2.0d0 * (scip(7)-scip(8))
c
c     get the permanent multipole and induced energies
c
               ei = -0.5d0 * (rr3*(gli(1)+gli(6))*psc3
     &                           + rr5*(gli(2)+gli(7))*psc5
     &                           + rr7*gli(3)*psc7)
               ei = f * ei
               es = es + ei
c
c     intermediate variables for the induced-permanent terms
c
               gfi(1) = 0.5d0*rr5*((gli(1)+gli(6))*psc3
     &                     +(glip(1)+glip(6))*dsc3+scip(2)*scale3i)
     &                + 0.5d0*rr7*((gli(7)+gli(2))*psc5
     &                            +(glip(7)+glip(2))*dsc5
     &                     -(sci(3)*scip(4)+scip(3)*sci(4))*scale5i)
     &                + 0.5d0*rr9*(gli(3)*psc7+glip(3)*dsc7)
               gfi(2) = -rr3*ck + rr5*sc(4) - rr7*sc(6)
               gfi(3) = rr3*ci + rr5*sc(3) + rr7*sc(5)
               gfi(4) = 2.0d0 * rr5
               gfi(5) = rr7 * (sci(4)*psc7+scip(4)*dsc7)
               gfi(6) = -rr7 * (sci(3)*psc7+scip(3)*dsc7)
c
c     get the induced force
c
               ftm2i(1) = gfi(1)*xr + 0.5d0*
     &            (- rr3*ck*(uind(1,i)*psc3+uinp(1,i)*dsc3)
     &             + rr5*sc(4)*(uind(1,i)*psc5+uinp(1,i)*dsc5)
     &             - rr7*sc(6)*(uind(1,i)*psc7+uinp(1,i)*dsc7))
     &             +(rr3*ci*(uind(1,k)*psc3+uinp(1,k)*dsc3)
     &             + rr5*sc(3)*(uind(1,k)*psc5+uinp(1,k)*dsc5)
     &             + rr7*sc(5)*(uind(1,k)*psc7+uinp(1,k)*dsc7))*0.5d0
     &             + rr5*scale5i*(sci(4)*uinp(1,i)+scip(4)*uind(1,i)
     &             + sci(3)*uinp(1,k)+scip(3)*uind(1,k))*0.5d0
     &             + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(1)
     &             + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(1)
     &             + 0.5d0*gfi(4)*((qkui(1)-qiuk(1))*psc5
     &             + (qkuip(1)-qiukp(1))*dsc5)
     &             + gfi(5)*qir(1) + gfi(6)*qkr(1)
               ftm2i(2) = gfi(1)*yr + 0.5d0*
     &            (- rr3*ck*(uind(2,i)*psc3+uinp(2,i)*dsc3)
     &             + rr5*sc(4)*(uind(2,i)*psc5+uinp(2,i)*dsc5)
     &             - rr7*sc(6)*(uind(2,i)*psc7+uinp(2,i)*dsc7))
     &             +(rr3*ci*(uind(2,k)*psc3+uinp(2,k)*dsc3)
     &             + rr5*sc(3)*(uind(2,k)*psc5+uinp(2,k)*dsc5)
     &             + rr7*sc(5)*(uind(2,k)*psc7+uinp(2,k)*dsc7))*0.5d0
     &             + rr5*scale5i*(sci(4)*uinp(2,i)+scip(4)*uind(2,i)
     &             + sci(3)*uinp(2,k)+scip(3)*uind(2,k))*0.5d0
     &             + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(2)
     &             + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(2)
     &             + 0.5d0*gfi(4)*((qkui(2)-qiuk(2))*psc5
     &             + (qkuip(2)-qiukp(2))*dsc5)
     &             + gfi(5)*qir(2) + gfi(6)*qkr(2)
               ftm2i(3) = gfi(1)*zr  + 0.5d0*
     &            (- rr3*ck*(uind(3,i)*psc3+uinp(3,i)*dsc3)
     &             + rr5*sc(4)*(uind(3,i)*psc5+uinp(3,i)*dsc5)
     &             - rr7*sc(6)*(uind(3,i)*psc7+uinp(3,i)*dsc7))
     &             +(rr3*ci*(uind(3,k)*psc3+uinp(3,k)*dsc3)
     &             + rr5*sc(3)*(uind(3,k)*psc5+uinp(3,k)*dsc5)
     &             + rr7*sc(5)*(uind(3,k)*psc7+uinp(3,k)*dsc7))*0.5d0
     &             + rr5*scale5i*(sci(4)*uinp(3,i)+scip(4)*uind(3,i)
     &             + sci(3)*uinp(3,k)+scip(3)*uind(3,k))*0.5d0
     &             + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(3)
     &             + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(3)
     &             + 0.5d0*gfi(4)*((qkui(3)-qiuk(3))*psc5
     &             + (qkuip(3)-qiukp(3))*dsc5)
     &             + gfi(5)*qir(3) + gfi(6)*qkr(3)
c
c     intermediate values needed for partially excluded interactions
c
               fridmp(1) = 0.5d0 * (rr3*((gli(1)+gli(6))*pscale(kk)
     &                        +(glip(1)+glip(6))*dscale(kk))*ddsc3(1)
     &            + rr5*((gli(2)+gli(7))*pscale(kk)
     &                +(glip(2)+glip(7))*dscale(kk))*ddsc5(1)
     &            + rr7*(gli(3)*pscale(kk)+glip(3)*dscale(kk))*ddsc7(1))
               fridmp(2) = 0.5d0 * (rr3*((gli(1)+gli(6))*pscale(kk)
     &                        +(glip(1)+glip(6))*dscale(kk))*ddsc3(2)
     &            + rr5*((gli(2)+gli(7))*pscale(kk)
     &                +(glip(2)+glip(7))*dscale(kk))*ddsc5(2)
     &            + rr7*(gli(3)*pscale(kk)+glip(3)*dscale(kk))*ddsc7(2))
               fridmp(3) = 0.5d0 * (rr3*((gli(1)+gli(6))*pscale(kk)
     &                        +(glip(1)+glip(6))*dscale(kk))*ddsc3(3)
     &            + rr5*((gli(2)+gli(7))*pscale(kk)
     &                +(glip(2)+glip(7))*dscale(kk))*ddsc5(3)
     &            + rr7*(gli(3)*pscale(kk)+glip(3)*dscale(kk))*ddsc7(3))
c
c     get the induced-induced derivative terms
c
               findmp(1) = 0.5d0 * uscale(kk) * (scip(2)*rr3*ddsc3(1)
     &                   - rr5*ddsc5(1)*(sci(3)*scip(4)+scip(3)*sci(4)))
               findmp(2) = 0.5d0 * uscale(kk) * (scip(2)*rr3*ddsc3(2)
     &                   - rr5*ddsc5(2)*(sci(3)*scip(4)+scip(3)*sci(4)))
               findmp(3) = 0.5d0 * uscale(kk) * (scip(2)*rr3*ddsc3(3)
     &                   - rr5*ddsc5(3)*(sci(3)*scip(4)+scip(3)*sci(4)))
c
c     handle of scaling for partially excluded interactions
c
               ftm2i(1) = ftm2i(1) - fridmp(1) - findmp(1)
               ftm2i(2) = ftm2i(2) - fridmp(2) - findmp(2)
               ftm2i(3) = ftm2i(3) - fridmp(3) - findmp(3)
c
c     correction to convert mutual to direct polarization force
c
               if (poltyp .eq. 'DIRECT') then
                  gfd = 0.5d0 * (rr5*scip(2)*scale3i
     &                  - rr7*(scip(3)*sci(4)+sci(3)*scip(4))*scale5i)
                  fdir(1) = gfd*xr + 0.5d0*rr5*scale5i
     &                         * (sci(4)*uinp(1,i)+scip(4)*uind(1,i)
     &                           +sci(3)*uinp(1,k)+scip(3)*uind(1,k))
                  fdir(2) = gfd*yr + 0.5d0*rr5*scale5i
     &                         * (sci(4)*uinp(2,i)+scip(4)*uind(2,i)
     &                           +sci(3)*uinp(2,k)+scip(3)*uind(2,k))
                  fdir(3) = gfd*zr + 0.5d0*rr5*scale5i
     &                         * (sci(4)*uinp(3,i)+scip(4)*uind(3,i)
     &                           +sci(3)*uinp(3,k)+scip(3)*uind(3,k))
                  ftm2i(1) = ftm2i(1) - fdir(1) + findmp(1)
                  ftm2i(2) = ftm2i(2) - fdir(2) + findmp(2)
                  ftm2i(3) = ftm2i(3) - fdir(3) + findmp(3)
               end if
c
c     now perform the torque calculation
c     intermediate terms for torque between multipoles i and k
c
               gti(2) = 0.5d0 * (sci(4)*psc5+scip(4)*dsc5) * rr5
               gti(3) = 0.5d0 * (sci(3)*psc5+scip(3)*dsc5) * rr5
               gti(4) = gfi(4)
               gti(5) = gfi(5)
               gti(6) = gfi(6)
c
c     calculate the induced torque components
c
               ttm2i(1) = -rr3*(dixuk(1)*psc3+dixukp(1)*dsc3)*0.5d0
     &            + gti(2)*dixr(1) + gti(4)*((ukxqir(1)+rxqiuk(1))*psc5
     &            +(ukxqirp(1)+rxqiukp(1))*dsc5)*0.5d0 - gti(5)*rxqir(1)
               ttm2i(2) = -rr3*(dixuk(2)*psc3+dixukp(2)*dsc3)*0.5d0
     &            + gti(2)*dixr(2) + gti(4)*((ukxqir(2)+rxqiuk(2))*psc5
     &            +(ukxqirp(2)+rxqiukp(2))*dsc5)*0.5d0 - gti(5)*rxqir(2)
               ttm2i(3) = -rr3*(dixuk(3)*psc3+dixukp(3)*dsc3)*0.5d0
     &            + gti(2)*dixr(3) + gti(4)*((ukxqir(3)+rxqiuk(3))*psc5
     &            +(ukxqirp(3)+rxqiukp(3))*dsc5)*0.5d0 - gti(5)*rxqir(3)
               ttm3i(1) = -rr3*(dkxui(1)*psc3+dkxuip(1)*dsc3)*0.5d0
     &            + gti(3)*dkxr(1) - gti(4)*((uixqkr(1)+rxqkui(1))*psc5
     &            +(uixqkrp(1)+rxqkuip(1))*dsc5)*0.5d0 - gti(6)*rxqkr(1)
               ttm3i(2) = -rr3*(dkxui(2)*psc3+dkxuip(2)*dsc3)*0.5d0
     &            + gti(3)*dkxr(2) - gti(4)*((uixqkr(2)+rxqkui(2))*psc5
     &            +(uixqkrp(2)+rxqkuip(2))*dsc5)*0.5d0 - gti(6)*rxqkr(2)
               ttm3i(3) = -rr3*(dkxui(3)*psc3+dkxuip(3)*dsc3)*0.5d0
     &            + gti(3)*dkxr(3) - gti(4)*((uixqkr(3)+rxqkui(3))*psc5
     &            +(uixqkrp(3)+rxqkuip(3))*dsc5)*0.5d0 - gti(6)*rxqkr(3)
c
c     update the force components on sites i and k
c
               des(1,ii) = des(1,ii) - f*ftm2i(1)
               des(2,ii) = des(2,ii) - f*ftm2i(2)
               des(3,ii) = des(3,ii) - f*ftm2i(3)
               des(1,kk) = des(1,kk) + f*ftm2i(1)
               des(2,kk) = des(2,kk) + f*ftm2i(2)
               des(3,kk) = des(3,kk) + f*ftm2i(3)
c
c     update the torque components on sites i and k
c
               trqi(1,ii) = trqi(1,ii) - f*ttm2i(1)
               trqi(2,ii) = trqi(2,ii) - f*ttm2i(2)
               trqi(3,ii) = trqi(3,ii) - f*ttm2i(3)
               trqi(1,kk) = trqi(1,kk) - f*ttm3i(1)
               trqi(2,kk) = trqi(2,kk) - f*ttm3i(2)
               trqi(3,kk) = trqi(3,kk) - f*ttm3i(3)
            end if
   10       continue
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = 1.0d0
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = 1.0d0
            uscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = 1.0d0
            uscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = 1.0d0
            uscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = 1.0d0
            uscale(ip14(j,ii)) = 1.0d0
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     convert torques into Cartesian forces
c
      call torque2 (trqi,des)
c
c     perform deallocation of some local arrays
c
      deallocate (pscale)
      deallocate (dscale)
      deallocate (uscale)
      deallocate (trqi)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine epb1  --  Poisson-Boltzmann energy and derivs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "epb1" calculates the implicit solvation energy and derivatives
c     via the Poisson-Boltzmann plus nonpolar implicit solvation
c
c
      subroutine epb1
      implicit none
c
c
c     compute the energy and gradients via Poisson-Boltzmann
c
      call epb1a
c
c     correct energy and derivatives for vacuum to polarized state
c
      call ediff1a
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine epb1a  --  PB solvation energy and derivatives  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "epb1a" calculates the solvation energy and gradients for the
c     PB/NP solvation model
c
c
      subroutine epb1a
      use sizes
      use atoms
      use chgpot
      use deriv
      use energi
      use mpole
      use pbstuf
      use polar
      use polpot
      use potent
      implicit none
      integer i,j,ii
      real*8 sum
      real*8, allocatable :: indpole(:,:)
      real*8, allocatable :: inppole(:,:)
      real*8, allocatable :: directf(:,:)
      real*8, allocatable :: directt(:,:)
      real*8, allocatable :: mutualf(:,:)
      real*8, allocatable :: polgrd(:,:)
      real*8, allocatable :: detor(:,:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (indpole(3,n))
      allocate (inppole(3,n))
      allocate (directf(3,n))
      allocate (directt(3,n))
      allocate (mutualf(3,n))
      allocate (polgrd(3,n))
c
c     induced dipole implicit energy via their
c     interaction with the permanent multipoles
c
      if (use_polar) then
         sum = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            sum = sum + uinds(1,i)*pbep(1,ii) + uinds(2,i)*pbep(2,ii)
     &               + uinds(3,i)*pbep(3,ii)
         end do
         sum = -0.5d0 * electric * sum
         pbe = pbe + sum
c
c     initialize induced dipole implicit energy gradients
c
         do i = 1, n
            do j = 1, 3
               indpole(j,i) = 0.0d0
               inppole(j,i) = 0.0d0
               directf(j,i) = 0.0d0
               directt(j,i) = 0.0d0
               mutualf(j,i) = 0.0d0
               polgrd(j,i) = 0.0d0
             end do
         end do
c
c     copy induced electrostatics into atom-based arrays
c
         do i = 1, npole
            ii = ipole(i)
            do j = 1, 3
               indpole(j,ii) = uinds(j,i)
               inppole(j,ii) = uinps(j,i)
            end do
         end do
c
c     perform dynamic allocation of some global arrays
c
         if (poltyp .eq. 'DIRECT') then
            if (.not. allocated(pbeuind))  allocate (pbeuind(3,n))
            if (.not. allocated(pbeuinp))  allocate (pbeuinp(3,n))
c
c     for direct polarization, the reaction field due to the
c     induced dipoles still needs to be computed because
c     the mutual portion of "apbsinduce" was not called
c
            do i = 1, n
               do j = 1, 3
                  pbeuind(j,i) = 0.0d0
                  pbeuinp(j,i) = 0.0d0
               end do
            end do
            call apbsinduce (indpole,pbeuind)
            call apbsnlinduce (inppole,pbeuinp)
         end if
c
c     compute direct induced dipole implicit solvation energy
c     gradients using potentials saved during the SCRF convergence
c
         call pbdirectpolforce (indpole,inppole,directf,directt)
c
c     convert torques due to induced dipole reaction field acting
c     on permanent multipoles into forces on adjacent atoms
c
         call torque2 (directt,polgrd)
         do i = 1, n
            polgrd(1,i) = polgrd(1,i) - directf(1,i)
            polgrd(2,i) = polgrd(2,i) - directf(2,i)
            polgrd(3,i) = polgrd(3,i) - directf(3,i)
         end do
c
c     compute mutual induced dipole solvation energy gradients
c
         if (poltyp .eq. 'MUTUAL') then
            call pbmutualpolforce (indpole,inppole,mutualf)
            do i = 1, n
               polgrd(1,i) = polgrd(1,i) - mutualf(1,i)
               polgrd(2,i) = polgrd(2,i) - mutualf(2,i)
               polgrd(3,i) = polgrd(3,i) - mutualf(3,i)
            end do
         end if
c
c     add induced dipole implicit solvation energy gradients
c     to overall polarization energy gradients
c
         do i = 1, n
            des(1,i) = des(1,i) + polgrd(1,i)
            des(2,i) = des(2,i) + polgrd(2,i)
            des(3,i) = des(3,i) + polgrd(3,i)
         end do
c
c     if polarization is off, get the permanent reaction field
c
      else
         call pbempole
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (indpole)
      deallocate (inppole)
      deallocate (directf)
      deallocate (directt)
      deallocate (mutualf)
      deallocate (polgrd)
c
c     increment solvation energy by Poisson-Boltzmann results
c
      es = es + pbe
c
c     perform dynamic allocation of some local arrays
c
      allocate (detor(3,n))
c
c     convert torques on permanent moments due to their own reaction
c     field into forces on adjacent atoms
c
      do i = 1, n
         do j = 1, 3
            detor(j,i) = 0.0d0
          end do
      end do
      call torque2 (pbtp,detor)
c
c     add permanent reaction field forces to the torque results
c
      do i = 1, n
         des(1,i) = des(1,i) - pbfp(1,i) + detor(1,i)
         des(2,i) = des(2,i) - pbfp(2,i) + detor(2,i)
         des(3,i) = des(3,i) - pbfp(3,i) + detor(3,i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (detor)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine enp1  --  cavity/dispersion energy and derivs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "enp1" calculates the nonpolar implicit solvation energy
c     and derivatives as a sum of cavity and dispersion terms
c
c
      subroutine enp1 (ecav,edisp)
      use sizes
      use atoms
      use atomid
      use deriv
      use kvdws
      use math
      use mpole
      use nonpol
      use shunt
      use solute
      implicit none
      integer i
      real*8 ecav,edisp
      real*8 exclude
      real*8 evol,esurf
      real*8 taperv,dtaperv
      real*8 tapersa,dtapersa
      real*8 reff,reff2,reff3
      real*8 reff4,reff5,dreff
      real*8, allocatable :: aesurf(:)
      real*8, allocatable :: dsurf(:,:)
      real*8, allocatable :: dvol(:,:)
      character*6 mode
c
c
c     zero out the nonpolar solvation energy and first derivatives
c
      ecav = 0.0d0
      edisp = 0.0d0
c
c     perform dynamic allocation of some local arrays
c
      allocate (aesurf(n))
      allocate (dsurf(3,n))
      allocate (dvol(3,n))
c
c     compute SASA and effective radius needed for cavity term
c
      exclude = 1.4d0
      call surface1 (esurf,aesurf,dsurf,rcav,asolv,exclude)
      reff = 0.5d0 * sqrt(esurf/(pi*surften))
      dreff = reff / (2.0d0*esurf)
      reff2 = reff * reff
      reff3 = reff2 * reff
      reff4 = reff3 * reff
      reff5 = reff4 * reff
c
c     compute solvent excluded volume needed for small solutes
c
      if (reff .lt. spoff) then
         call volume (evol,rcav,exclude)
         evol = evol * solvprs
         call volume1 (rcav,exclude,dvol)
         do i = 1, n
            dvol(1,i) = dvol(1,i) * solvprs
            dvol(2,i) = dvol(2,i) * solvprs
            dvol(3,i) = dvol(3,i) * solvprs
         end do
      end if
c
c     find cavity energy from only the solvent excluded volume
c
      if (reff .le. spcut) then
         ecav = evol
         do i = 1, n
            des(1,i) = des(1,i) + dvol(1,i)
            des(2,i) = des(2,i) + dvol(2,i)
            des(3,i) = des(3,i) + dvol(3,i)
         end do
c
c     find cavity energy from only a tapered volume term
c
      else if (reff.gt.spcut .and. reff.le.stoff) then
         mode = 'GKV'
         call switch (mode)
         taperv = c5*reff5 + c4*reff4 + c3*reff3
     &               + c2*reff2 + c1*reff + c0
         dtaperv = (5.0d0*c5*reff4+4.0d0*c4*reff3+3.0d0*c3*reff2
     &                 +2.0d0*c2*reff+c1) * dreff
         ecav = evol * taperv
         do i = 1, n
            des(1,i) = des(1,i) + taperv*dvol(1,i)
     &                    + evol*dtaperv*dsurf(1,i)
            des(2,i) = des(2,i) + taperv*dvol(2,i)
     &                    + evol*dtaperv*dsurf(2,i)
            des(3,i) = des(3,i) + taperv*dvol(3,i)
     &                    + evol*dtaperv*dsurf(3,i)
         end do
c
c     find cavity energy using both volume and SASA terms
c
      else if (reff.gt.stoff .and. reff.le.spoff) then
         mode = 'GKV'
         call switch (mode)
         taperv = c5*reff5 + c4*reff4 + c3*reff3
     &               + c2*reff2 + c1*reff + c0
         dtaperv = (5.0d0*c5*reff4+4.0d0*c4*reff3+3.0d0*c3*reff2
     &                 +2.0d0*c2*reff+c1) * dreff
         mode = 'GKSA'
         call switch (mode)
         tapersa = c5*reff5 + c4*reff4 + c3*reff3
     &                + c2*reff2 + c1*reff + c0
         tapersa = 1.0d0 - tapersa
         dtapersa = (5.0d0*c5*reff4+4.0d0*c4*reff3+3.0d0*c3*reff2
     &                  +2.0d0*c2*reff+c1) * dreff
         dtapersa = -dtapersa
         ecav = evol * taperv
         do i = 1, n
            des(1,i) = des(1,i) + taperv*dvol(1,i)
     &                    + evol*dtaperv*dsurf(1,i)
            des(2,i) = des(2,i) + taperv*dvol(2,i)
     &                    + evol*dtaperv*dsurf(2,i)
            des(3,i) = des(3,i) + taperv*dvol(3,i)
     &                    + evol*dtaperv*dsurf(3,i)
         end do
         ecav = ecav + tapersa*esurf
         do i = 1, n
            des(1,i) = des(1,i) + (tapersa+esurf*dtapersa)*dsurf(1,i)
            des(2,i) = des(2,i) + (tapersa+esurf*dtapersa)*dsurf(2,i)
            des(3,i) = des(3,i) + (tapersa+esurf*dtapersa)*dsurf(3,i)
         end do
c
c     find cavity energy from only a tapered SASA term
c
      else if (reff.gt.spoff .and. reff.le.stcut) then
         mode = 'GKSA'
         call switch (mode)
         tapersa = c5*reff5 + c4*reff4 + c3*reff3
     &                + c2*reff2 + c1*reff + c0
         tapersa = 1.0d0 - tapersa
         dtapersa = (5.0d0*c5*reff4+4.0d0*c4*reff3+3.0d0*c3*reff2
     &                  +2.0d0*c2*reff+c1) * dreff
         dtapersa = -dtapersa
         ecav = tapersa * esurf
         do i = 1, n
            des(1,i) = des(1,i) + (tapersa+esurf*dtapersa)*dsurf(1,i)
            des(2,i) = des(2,i) + (tapersa+esurf*dtapersa)*dsurf(2,i)
            des(3,i) = des(3,i) + (tapersa+esurf*dtapersa)*dsurf(3,i)
         end do
c
c     find cavity energy from only a SASA-based term
c
      else
         ecav = esurf
         do i = 1, n
            des(1,i) = des(1,i) + dsurf(1,i)
            des(2,i) = des(2,i) + dsurf(2,i)
            des(3,i) = des(3,i) + dsurf(3,i)
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (aesurf)
      deallocate (dsurf)
      deallocate (dvol)
c
c     find the implicit dispersion solvation energy
c
      call ewca1 (edisp)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine ewca1  --  WCA dispersion energy and derivs  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "ewca1" finds the Weeks-Chandler-Anderson dispersion energy
c     and derivatives of a solute
c
c
      subroutine ewca1 (edisp)
      use sizes
      use atoms
      use atomid
      use deriv
      use kvdws
      use math
      use nonpol
      use solute
      use vdw
      implicit none
      integer i,k
      real*8 edisp,e,idisp
      real*8 xi,yi,zi
      real*8 rk,sk,sk2
      real*8 xr,yr,zr
      real*8 r,r2,r3
      real*8 sum,term,shctd
      real*8 iwca,irep,offset
      real*8 epsi,rmini,ri,rmax
      real*8 ao,emixo,rmixo,rmixo7
      real*8 ah,emixh,rmixh,rmixh7
      real*8 lik,lik2,lik3,lik4
      real*8 lik5,lik6,lik10
      real*8 lik11,lik12,lik13
      real*8 uik,uik2,uik3,uik4
      real*8 uik5,uik6,uik10
      real*8 uik11,uik12,uik13
      real*8 de,dl,du
      real*8 dedx,dedy,dedz
c
c
c     zero out the Weeks-Chandler-Andersen dispersion energy
c
      edisp = 0.0d0
c
c     set overlap scale factor for HCT descreening method
c
      shctd = 0.81d0
      offset = 0.0d0
      do i = 1, n
         rdisp(i) = rad(class(i)) + offset
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(n,class,eps,
!$OMP& rad,rdisp,x,y,z,shctd,cdisp)
!$OMP& shared(edisp,des)
!$OMP DO reduction(+:edisp,des) schedule(guided)
c
c     find the WCA dispersion energy and gradient components
c
      do i = 1, n
         epsi = eps(class(i))
         rmini = rad(class(i))
         emixo = 4.0d0 * epso * epsi / ((sqrt(epso)+sqrt(epsi))**2)
         rmixo = 2.0d0 * (rmino**3+rmini**3) / (rmino**2+rmini**2)
         rmixo7 = rmixo**7
         ao = emixo * rmixo7
         emixh = 4.0d0 * epsh * epsi / ((sqrt(epsh)+sqrt(epsi))**2)
         rmixh = 2.0d0 * (rminh**3+rmini**3) / (rminh**2+rmini**2)
         rmixh7 = rmixh**7
         ah = emixh * rmixh7
         ri = rdisp(i)
c
c     remove contribution due to solvent displaced by solute atoms
c
         xi = x(i)
         yi = y(i)
         zi = z(i)
         sum = 0.0d0
         do k = 1, n
            if (i .ne. k) then
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
               r2 = xr*xr + yr*yr + zr*zr
               r = sqrt(r2)
               r3 = r * r2
               rk = rdisp(k)
c              sk = rk * shct(k)
               sk = rk * shctd
               sk2 = sk * sk
               if (ri .lt. r+sk) then
                  de = 0.0d0
                  rmax = max(ri,r-sk)
                  lik = rmax
                  lik2 = lik * lik
                  lik3 = lik2 * lik
                  lik4 = lik3 * lik
                  if (lik .lt. rmixo) then
                     uik = min(r+sk,rmixo)
                     uik2 = uik * uik
                     uik3 = uik2 * uik
                     uik4 = uik3 * uik
                     term = 4.0d0 * pi / (48.0d0*r)
     &                      * (3.0d0*(lik4-uik4) - 8.0d0*r*(lik3-uik3)
     &                          + 6.0d0*(r2-sk2)*(lik2-uik2))
                     if (ri .gt. r-sk) then
                        dl = -lik2 + 2.0d0*r2 + 2.0d0*sk2
                        dl = dl * lik2
                     else
                        dl = -lik3 + 4.0d0*lik2*r - 6.0d0*lik*r2
     &                          + 2.0d0*lik*sk2 + 4.0d0*r3 - 4.0d0*r*sk2
                        dl = dl * lik
                     end if
                     if (r+sk .gt. rmixo) then
                        du = -uik2 + 2.0d0*r2 + 2.0d0*sk2
                        du = -du * uik2
                     else
                        du = -uik3 + 4.0d0*uik2*r - 6.0d0*uik*r2
     &                          + 2.0d0*uik*sk2 + 4.0d0*r3 - 4.0d0*r*sk2
                        du = -du * uik
                     end if
                     iwca = -emixo * term
                     de = de - emixo*pi*(dl+du)/(4.0d0*r2)
                     sum = sum + iwca
                  end if
                  if (lik .lt. rmixh) then
                     uik = min(r+sk,rmixh)
                     uik2 = uik * uik
                     uik3 = uik2 * uik
                     uik4 = uik3 * uik
                     term = 4.0d0 * pi / (48.0d0*r)
     &                      * (3.0d0*(lik4-uik4) - 8.0d0*r*(lik3-uik3)
     &                          + 6.0d0*(r2-sk2)*(lik2-uik2))
                     if (ri .gt. r-sk) then
                        dl = -lik2 + 2.0d0*r2 + 2.0d0*sk2
                        dl = dl * lik2
                     else
                        dl = -lik3 + 4.0d0*lik2*r - 6.0d0*lik*r2
     &                          + 2.0d0*lik*sk2 + 4.0d0*r3 - 4.0d0*r*sk2
                        dl = dl * lik
                     end if
                     if (r+sk .gt. rmixh) then
                        du = -uik2 + 2.0d0*r2 + 2.0d0*sk2
                        du = -du * uik2
                     else
                        du = -uik3 + 4.0d0*uik2*r - 6.0d0*uik*r2
     &                          + 2.0d0*uik*sk2 + 4.0d0*r3 - 4.0d0*r*sk2
                        du = -du * uik
                     end if
                     iwca = -2.0d0 * emixh * term
                     de = de - 2.0d0*emixh*pi*(dl+du)/(4.0d0*r2)
                     sum = sum + iwca
                  end if
                  uik = r + sk
                  uik2 = uik * uik
                  uik3 = uik2 * uik
                  uik4 = uik3 * uik
                  uik5 = uik4 * uik
                  uik6 = uik5 * uik
                  uik10 = uik5 * uik5
                  uik11 = uik10 * uik
                  uik12 = uik11 * uik
                  uik13 = uik12 * uik
                  if (uik .gt. rmixo) then
                     lik = max(rmax,rmixo)
                     lik2 = lik * lik
                     lik3 = lik2 * lik
                     lik4 = lik3 * lik
                     lik5 = lik4 * lik
                     lik6 = lik5 * lik
                     lik10 = lik5 * lik5
                     lik11 = lik10 * lik
                     lik12 = lik11 * lik
                     lik13 = lik12 * lik
                     term = 4.0d0 * pi / (120.0d0*r*lik5*uik5)
     &                      * (15.0d0*uik*lik*r*(uik4-lik4)
     &                         - 10.0d0*uik2*lik2*(uik3-lik3)
     &                         + 6.0d0*(sk2-r2)*(uik5-lik5))
                     if (ri.gt.r-sk .or. rmax.lt.rmixo) then
                        dl = -5.0d0*lik2 + 3.0d0*r2 + 3.0d0*sk2
                        dl = -dl / lik5
                     else
                        dl = 5.0d0*lik3 - 33.0d0*lik*r2 - 3.0d0*lik*sk2
     &                          + 15.0d0*(lik2*r+r3-r*sk2)
                        dl = dl / lik6
                     end if
                     du = 5.0d0*uik3 - 33.0d0*uik*r2 - 3.0d0*uik*sk2
     &                       + 15.0d0*(uik2*r+r3-r*sk2)
                     du = -du / uik6
                     idisp = -2.0d0 * ao * term
                     de = de -2.0d0*ao*pi*(dl + du)/(15.0d0*r2)
                     term = 4.0d0 * pi / (2640.0d0*r*lik12*uik12)
     &                      * (120.0d0*uik*lik*r*(uik11-lik11)
     &                         - 66.0d0*uik2*lik2*(uik10-lik10)
     &                         + 55.0d0*(sk2-r2)*(uik12-lik12))
                     if (ri.gt.r-sk .or. rmax.lt.rmixo) then
                        dl = -6.0d0*lik2 + 5.0d0*r2 + 5.0d0*sk2
                        dl = -dl / lik12
                     else
                        dl = 6.0d0*lik3 - 125.0d0*lik*r2 - 5.0d0*lik*sk2
     &                          + 60.0d0*(lik2*r+r3-r*sk2)
                        dl = dl / lik13
                     end if
                     du = 6.0d0*uik3 - 125.0d0*uik*r2 -5.0d0*uik*sk2
     &                       + 60.0d0*(uik2*r+r3-r*sk2)
                     du = -du / uik13
                     irep = ao * rmixo7 * term
                     de = de + ao*rmixo7*pi*(dl + du)/(60.0d0*r2)
                     sum = sum + irep + idisp
                  end if
                  if (uik .gt. rmixh) then
                     lik = max(rmax,rmixh)
                     lik2 = lik * lik
                     lik3 = lik2 * lik
                     lik4 = lik3 * lik
                     lik5 = lik4 * lik
                     lik6 = lik5 * lik
                     lik10 = lik5 * lik5
                     lik11 = lik10 * lik
                     lik12 = lik11 * lik
                     lik13 = lik12 * lik
                     term = 4.0d0 * pi / (120.0d0*r*lik5*uik5)
     &                      * (15.0d0*uik*lik*r*(uik4-lik4)
     &                         - 10.0d0*uik2*lik2*(uik3-lik3)
     &                         + 6.0d0*(sk2-r2)*(uik5-lik5))
                     if (ri.gt.r-sk .or. rmax.lt.rmixh) then
                        dl = -5.0d0*lik2 + 3.0d0*r2 + 3.0d0*sk2
                        dl = -dl / lik5
                     else
                        dl = 5.0d0*lik3 - 33.0d0*lik*r2 - 3.0d0*lik*sk2
     &                          + 15.0d0*(lik2*r+r3-r*sk2)
                        dl = dl / lik6
                     end if
                     du = 5.0d0*uik3 - 33.0d0*uik*r2 - 3.0d0*uik*sk2
     &                       + 15.0d0*(uik2*r+r3-r*sk2)
                     du = -du / uik6
                     idisp = -4.0d0 * ah * term
                     de = de - 4.0d0*ah*pi*(dl + du)/(15.0d0*r2)
                     term = 4.0d0 * pi / (2640.0d0*r*lik12*uik12)
     &                      * (120.0d0*uik*lik*r*(uik11-lik11)
     &                         - 66.0d0*uik2*lik2*(uik10-lik10)
     &                         + 55.0d0*(sk2-r2)*(uik12-lik12))
                     if (ri.gt.r-sk .or. rmax.lt.rmixh) then
                        dl = -6.0d0*lik2 + 5.0d0*r2 + 5.0d0*sk2
                        dl = -dl / lik12
                     else
                        dl = 6.0d0*lik3 - 125.0d0*lik*r2 - 5.0d0*lik*sk2
     &                          + 60.0d0*(lik2*r+r3-r*sk2)
                        dl = dl / lik13
                     end if
                     du = 6.0d0*uik3 - 125.0d0*uik*r2 -5.0d0*uik*sk2
     &                       + 60.0d0*(uik2*r+r3-r*sk2)
                     du = -du / uik13
                     irep = 2.0d0 * ah * rmixh7 * term
                     de = de + ah*rmixh7*pi*(dl+du)/(30.0d0*r2)
                     sum = sum + irep + idisp
                  end if
c
c     increment the individual dispersion gradient components
c
                  de = -de/r * slevy * awater
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
                  des(1,i) = des(1,i) + dedx
                  des(2,i) = des(2,i) + dedy
                  des(3,i) = des(3,i) + dedz
                  des(1,k) = des(1,k) - dedx
                  des(2,k) = des(2,k) - dedy
                  des(3,k) = des(3,k) - dedz
               end if
            end if
         end do
c
c     increment the overall dispersion energy component
c
         e = cdisp(i) - slevy*awater*sum
         edisp = edisp + e
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine ehpmf1  --  HPMF nonpolar solvation and derivs  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ehpmf1" calculates the hydrophobic potential of mean force
c     energy and first derivatives using a pairwise double loop
c
c     literature reference:
c
c     M. S. Lin, N. L. Fawzi and T. Head-Gordon, "Hydrophobic
c     Potential of Mean Force as a Solvation Function for Protein
c     Structure Prediction", Structure, 15, 727-740 (2007)
c
c
      subroutine ehpmf1 (ehp)
      use sizes
      use atomid
      use atoms
      use couple
      use deriv
      use hpmf
      use iounit
      use math
      implicit none
      integer i,j,k,m
      integer ii,jj,kk
      integer sschk
      integer, allocatable :: omit(:)
      real*8 xi,yi,zi
      real*8 xk,yk,zk
      real*8 xr,yr,zr
      real*8 e,ehp,r,r2
      real*8 rsurf,pisurf
      real*8 hpmfcut2
      real*8 saterm,sasa
      real*8 rbig,rsmall
      real*8 part,cutv
      real*8 t1a,t1b
      real*8 e1,e2,e3,sum
      real*8 arg1,arg2,arg3
      real*8 arg12,arg22,arg32
      real*8 de1,de2,de3,dsum
      real*8 sumi,sumj,term
      real*8 dedx,dedy,dedz
      real*8, allocatable :: cutmtx(:)
      real*8, allocatable :: dcutmtx(:)
      real*8, allocatable :: dacsa(:,:)
c
c
c     zero out the hydrophobic potential of mean force energy
c
      ehp = 0.0d0
c
c     set some values needed during the HPMF calculation
c
      rsurf = rcarbon + 2.0d0*rwater
      pisurf = pi * (rcarbon+rwater)
      hpmfcut2 = hpmfcut * hpmfcut
c
c     perform dynamic allocation of some local arrays
c
      allocate (omit(n))
      allocate (cutmtx(n))
      allocate (dcutmtx(n))
      allocate (dacsa(n,npmf))
c
c     get the surface area and derivative terms for each atom
c
      do ii = 1, npmf
         i = ipmf(ii)
         saterm = acsa(i)
         sasa = 1.0d0
         do k = 1, n
            if (i .ne. k) then
               xr = x(i) - x(k)
               yr = y(i) - y(k)
               zr = z(i) - z(k)
               r2 = xr*xr + yr*yr + zr*zr
               rbig = rpmf(k) + rsurf
               if (r2 .le. rbig*rbig) then
                  r = sqrt(r2)
                  rsmall = rpmf(k) - rcarbon
                  part = pisurf * (rbig-r) * (1.0d0+rsmall/r)
                  sasa = sasa * (1.0d0-saterm*part)
               end if
            end if
         end do
         sasa = acsurf * sasa
         cutv = tanh(tslope*(sasa-toffset))
         cutmtx(i) = 0.5d0 * (1.0d0+cutv)
         dcutmtx(i) = 0.5d0 * tslope * (1.0d0-cutv*cutv)
         do k = 1, n
            dacsa(k,ii) = 0.0d0
            if (i .ne. k) then
               xr = x(i) - x(k)
               yr = y(i) - y(k)
               zr = z(i) - z(k)
               r2 = xr*xr + yr*yr + zr*zr
               rbig = rpmf(k) + rsurf
               if (r2 .le. rbig*rbig) then
                  r = sqrt(r2)
                  rsmall = rpmf(k) - rcarbon
                  part = pisurf * (rbig-r) * (1.0d0+rsmall/r)
                  t1b = -pisurf * (1.0d0+rbig*rsmall/r2)
                  t1a = -sasa / (1.0d0/saterm-part)
                  dacsa(k,ii) = t1a * t1b / r
               end if
            end if
         end do
      end do
c
c     find the hydrophobic PMF energy and derivs via a double loop
c
      do i = 1, n
         omit(i) = 0
      end do
      do ii = 1, npmf-1
         i = ipmf(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         sschk = 0
         do j = 1, n12(i)
            k = i12(j,i)
            omit(k) = i
            if (atomic(k) .eq. 16)  sschk = k
         end do
         do j = 1, n13(i)
            k = i13(j,i)
            omit(k) = i
         end do
         do j = 1, n14(i)
            k = i14(j,i)
            omit(k) = i
            if (sschk .ne. 0) then
               do jj = 1, n12(k)
                  m = i12(jj,k)
                  if (atomic(m) .eq. 16) then
                     do kk = 1, n12(m)
                        if (i12(kk,m) .eq. sschk)  omit(k) = 0
                     end do
                  end if
               end do
            end if
         end do
         do kk = ii+1, npmf
            k = ipmf(kk)
            xk = x(k)
            yk = y(k)
            zk = z(k)
            if (omit(k) .ne. i) then
               xr = xi - xk
               yr = yi - yk
               zr = zi - zk
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. hpmfcut2) then
                  r = sqrt(r2)
                  arg1 = (r-c1) * w1
                  arg12 = arg1 * arg1
                  arg2 = (r-c2) * w2
                  arg22 = arg2 * arg2
                  arg3 = (r-c3) * w3
                  arg32 = arg3 * arg3
                  e1 = h1 * exp(-arg12)
                  e2 = h2 * exp(-arg22)
                  e3 = h3 * exp(-arg32)
                  sum = e1 + e2 + e3
                  e = sum * cutmtx(i) * cutmtx(k)
                  ehp = ehp + e
c
c     first part of hydrophobic PMF derivative calculation
c
                  de1 = -2.0d0 * e1 * arg1 * w1
                  de2 = -2.0d0 * e2 * arg2 * w2
                  de3 = -2.0d0 * e3 * arg3 * w3
                  dsum = (de1+de2+de3) * cutmtx(i) * cutmtx(k) / r
                  dedx = dsum * xr
                  dedy = dsum * yr
                  dedz = dsum * zr
                  des(1,i) = des(1,i) + dedx
                  des(2,i) = des(2,i) + dedy
                  des(3,i) = des(3,i) + dedz
                  des(1,k) = des(1,k) - dedx
                  des(2,k) = des(2,k) - dedy
                  des(3,k) = des(3,k) - dedz
c
c     second part of hydrophobic PMF derivative calculation
c
                  sumi = sum * cutmtx(k) * dcutmtx(i)
                  sumj = sum * cutmtx(i) * dcutmtx(k)
                  if (sumi .ne. 0.0d0) then
                     do j = 1, n
                        if (dacsa(j,ii) .ne. 0.0d0) then
                           term = sumi * dacsa(j,ii)
                           dedx = term * (xi-x(j))
                           dedy = term * (yi-y(j))
                           dedz = term * (zi-z(j))
                           des(1,i) = des(1,i) + dedx
                           des(2,i) = des(2,i) + dedy
                           des(3,i) = des(3,i) + dedz
                           des(1,j) = des(1,j) - dedx
                           des(2,j) = des(2,j) - dedy
                           des(3,j) = des(3,j) - dedz
                        end if
                     end do
                  end if
                  if (sumj .ne. 0.0d0) then
                     do j = 1, n
                        if (dacsa(j,kk) .ne. 0.0d0) then
                           term = sumj * dacsa(j,kk)
                           dedx = term * (xk-x(j))
                           dedy = term * (yk-y(j))
                           dedz = term * (zk-z(j))
                           des(1,k) = des(1,k) + dedx
                           des(2,k) = des(2,k) + dedy
                           des(3,k) = des(3,k) + dedz
                           des(1,j) = des(1,j) - dedx
                           des(2,j) = des(2,j) - dedy
                           des(3,j) = des(3,j) - dedz
                        end if
                     end do
                  end if
               end if
            end if
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (omit)
      deallocate (cutmtx)
      deallocate (dcutmtx)
      deallocate (dacsa)
      return
      end
