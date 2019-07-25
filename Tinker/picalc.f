c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine picalc  --  Pariser-Parr-Pople MO calculation  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "picalc" performs a modified Pariser-Parr-Pople molecular
c     orbital calculation for each conjugated pisystem
c
c
      subroutine picalc
      use sizes
      use bndstr
      use couple
      use inform
      use iounit
      use piorbs
      use tors
      implicit none
      integer i,j,k,m,ib,ic
      integer ii,jj,kk
      integer iorb,jorb
      integer ncalls
      data ncalls  / 0 /
      save ncalls
c
c
c     only needs to be done if pisystem is present
c
      if (norbit .eq. 0)  return
c
c     increment the number of calls to this routine
c
      ncalls = ncalls + 1
      if (reorbit.eq.0 .or. ncalls.lt.reorbit)  return
      ncalls = 0
c
c     loop over all pisystems computing separate MOs for each
c
      do i = 1, nconj
         norbit = 0
         do j = iconj(1,i), iconj(2,i)
            norbit = norbit + 1
            iorbit(norbit) = kconj(j)
         end do
c
c     find and store the pisystem bonds for this pisystem
c
         nbpi = 0
         kk = iconj(2,i) - iconj(1,i) + 1
         do ii = 1, norbit-1
            iorb = iorbit(ii)
            do jj = ii+1, norbit
               jorb = iorbit(jj)
               do k = 1, n12(iorb)
                  if (i12(k,iorb) .eq. jorb) then
                     nbpi = nbpi + 1
                     do m = 1, nbond
                        if (iorb.eq.ibnd(1,m) .and.
     &                      jorb.eq.ibnd(2,m)) then
                           ibpi(1,nbpi) = m
                           ibpi(2,nbpi) = ii
                           ibpi(3,nbpi) = jj
                           goto 10
                        end if
                     end do
   10                continue
                  end if
               end do
            end do
         end do
c
c     find and store the pisystem torsions for this pisystem
c
         ntpi = 0
         do ii = 1, ntors
            ib = itors(2,ii)
            ic = itors(3,ii)
            if (listpi(ib) .and. listpi(ic)) then
               do jj = 1, nbpi
                  k = ibpi(1,jj)
                  if (ib.eq.ibnd(1,k).and.ic.eq.ibnd(2,k) .or.
     &                ib.eq.ibnd(2,k).and.ic.eq.ibnd(1,k)) then
                     ntpi = ntpi + 1
                     itpi(1,ntpi) = ii
                     itpi(2,ntpi) = jj
                     goto 20
                  end if
               end do
   20          continue
            end if
         end do
c
c     print a header for the molecular orbital calculation
c
         if (debug) then
            if (nconj .eq. 1) then
               write (iout,30)
   30          format (/,' Modified Pariser-Parr-Pople Molecular',
     &                    ' Orbitals :')
            else
               write (iout,40)  i
   40          format (/,' Modified Pariser-Parr-Pople MOs for',
     &                    ' Pi-System',i4,' :')
            end if
         end if
c
c     get SCF-MOs, then scale bond and torsional parameters
c
         call piscf
         call pialter
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine piscf  --  SCF molecular orbital calculation  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "piscf" performs an SCF molecular orbital calculation for a
c     pisystem to determine bond orders used in parameter scaling
c
c
      subroutine piscf
      use sizes
      use atomid
      use atoms
      use bndstr
      use couple
      use inform
      use iounit
      use orbits
      use piorbs
      use pistuf
      use units
      implicit none
      integer i,j,k,m
      integer iter,maxiter
      integer iatn,jatn
      integer iorb,jorb,nfill
      real*8 delta,converge
      real*8 xij,yij,zij,p
      real*8 hcii,gii,gij
      real*8 g11,g11sq,g12,g14
      real*8 rij,erij,brij
      real*8 ovlap,covlap
      real*8 cionize
      real*8 iionize,jionize
      real*8 rijsq,hcij,qi
      real*8 total,totold
      real*8 ebeta,aeth,abnz
      real*8 ebe,ebb,ble,blb
      real*8 eebond,bebond
      real*8 s1,s2,gjk
      real*8 vij,vik,vmj,vmk
      real*8 xi,xj,xk,xg
      real*8, allocatable :: povlap(:)
      real*8, allocatable :: en(:)
      real*8, allocatable :: ip(:)
      real*8, allocatable :: fock(:,:)
      real*8, allocatable :: hc(:,:)
      real*8, allocatable :: v(:,:)
      real*8, allocatable :: gamma(:,:)
      real*8, allocatable :: ed(:,:)
      character*6 mode
c
c
c     initialize some constants and parameters
c
c     mode       planar or nonplanar pi-calculation
c     maxiter    maximum number of SCF iterations
c     converge   criterion for SCF convergence
c     ebeta      value of resonance integral for ethylene
c     cionize    ionization potential of carbon (Hartree)
c
      mode = 'PLANAR'
      maxiter = 50
      converge = 0.00000001d0
      ebeta = -0.0757d0
      cionize = -11.16d0 / evolt
c
c     set the bond energies, alpha values and ideal bond length
c     parameter for carbon-carbon pibond type parameters
c
c     ebe    equilibrium bond energy in ethylene
c     ebb    equilibrium bond energy in benzene
c     aeth   the P-P-P constant "a" in ethylene
c     abnz   the P-P-P constant "a" in benzene
c     ble    equilibrium bond length in ethylene
c     blb    equilibrium bond length in benzene
c
      ebe = 129.37d0
      ebb = 117.58d0
      aeth = 2.309d0
      abnz = 2.142d0
      ble = 1.338d0
      blb = 1.397d0
c
c     perform dynamic allocation of some local arrays
c
      allocate (povlap(nbpi))
      allocate (en(norbit))
      allocate (ip(norbit))
      allocate (fock(norbit,norbit))
      allocate (hc(norbit,norbit))
      allocate (v(norbit,norbit))
      allocate (gamma(norbit,norbit))
      allocate (ed(norbit,norbit))
c
c     assign empirical one-center Coulomb integrals, and
c     first or second ionization potential depending on
c     whether the orbital contributes one or two electrons
c
      nfill = 0
      do i = 1, norbit
         iorb = iorbit(i)
         gamma(i,i) = emorb(iorb)
         ip(i) = worb(iorb) + (1.0d0-qorb(iorb))*emorb(iorb)
         nfill = nfill + nint(qorb(iorb))
      end do
      nfill = nfill / 2
c
c     calculate two-center repulsion integrals
c     according to Ohno's semi-empirical formula
c
      do i = 1, norbit-1
         iorb = iorbit(i)
         gii = gamma(i,i)
         do j = i+1, norbit
            jorb = iorbit(j)
            g11 = 0.5d0 * (gii+gamma(j,j))
            g11sq = 1.0d0 / g11**2
            xij = x(iorb) - x(jorb)
            yij = y(iorb) - y(jorb)
            zij = z(iorb) - z(jorb)
            rijsq = (xij**2 + yij**2 + zij**2) / bohr**2
            g12 = 1.0d0 / sqrt(rijsq+g11sq)
            gamma(i,j) = g12
            gamma(j,i) = g12
         end do
      end do
c
c     zero out the resonance integral values
c
      do i = 1, norbit
         do j = 1, norbit
            hc(j,i) = 0.0d0
         end do
      end do
c
c     the first term in the sum to find alpha is the first
c     or second ionization potential, then the two-center
c     repulsion integrals are added
c
      do i = 1, norbit
         hcii = ip(i)
         do j = 1, norbit
            if (i .ne. j) then
               jorb = iorbit(j)
               hcii = hcii - qorb(jorb)*gamma(i,j)
            end if
         end do
         hc(i,i) = hcii
      end do
c
c     get two-center repulsion integrals via Ohno's formula
c
      do k = 1, nbpi
         i = ibpi(2,k)
         j = ibpi(3,k)
         iorb = iorbit(i)
         jorb = iorbit(j)
         iatn = atomic(iorb)
         jatn = atomic(jorb)
         xij = x(iorb) - x(jorb)
         yij = y(iorb) - y(jorb)
         zij = z(iorb) - z(jorb)
         rij = sqrt(xij**2 + yij**2 + zij**2)
         rijsq = rij**2 / bohr**2
         g11 = 0.5d0 * (gamma(i,i)+gamma(j,j))
         g11sq = 1.0d0 / g11**2
         g12 = gamma(i,j)
c
c     compute the bond energy using a Morse potential
c
         erij = aeth * (ble-rij)
         brij = abnz * (blb-rij)
         eebond = (2.0d0*exp(erij)-exp(2.0d0*erij)) * ebe / hartree
         bebond = (2.0d0*exp(brij)-exp(2.0d0*brij)) * ebb / hartree
c
c     compute the carbon-carbon resonance integral using
c     the Whitehead and Lo formula
c
         g14 = 1.0d0 / sqrt(4.0d0*rijsq+g11sq)
         hcij = 1.5d0*(bebond-eebond) - 0.375d0*g11
     &             + (5.0d0/12.0d0)*g12 - g14/24.0d0
c
c     if either atom is non-carbon, then factor the resonance
c     integral by overlap ratio and ionization potential ratio
c
         if (iatn.ne.6 .or. jatn.ne.6) then
            call overlap (iatn,jatn,rij,ovlap)
            call overlap (6,6,rij,covlap)
            hcij = hcij * (ovlap/covlap)
            iionize = ip(i)
            if (qorb(iorb) .ne. 1.0d0) then
               if (iatn .eq. 7)  iionize = 0.595d0 * iionize
               if (iatn .eq. 8)  iionize = 0.525d0 * iionize
               if (iatn .eq. 16)  iionize = 0.89d0 * iionize
            end if
            jionize = ip(j)
            if (qorb(jorb) .ne. 1.0d0) then
               if (jatn .eq. 7)  jionize = 0.595d0 * jionize
               if (jatn .eq. 8)  jionize = 0.525d0 * jionize
               if (jatn .eq. 16)  jionize = 0.89d0 * jionize
            end if
            hcij = hcij * (iionize+jionize)/(2.0d0*cionize)
         end if
c
c     set symmetric elements to the same value
c
         hc(i,j) = hcij
         hc(j,i) = hcij
      end do
c
c     construct an initial guess to the Fock matrix
c
      do i = 1, norbit
         do j = 1, norbit
            fock(j,i) = hc(j,i)
         end do
      end do
      do i = 1, norbit
         fock(i,i) = 0.5d0 * ip(i)
      end do
c
c     make the SCF-MO computation; do it twice, for a planar analog
c     of the actual system and for the actual (nonplanar) system
c
      do while (mode.eq.'PLANAR' .or. mode.eq.'NONPLN')
         if (mode .eq. 'NONPLN') then
            call pitilt (povlap)
            do k = 1, nbpi
               i = ibpi(2,k)
               j = ibpi(3,k)
               hc(i,j) = hc(i,j) * povlap(k)
               hc(j,i) = hc(i,j)
            end do
         end if
c
c     perform SCF iterations until convergence is reached; diagonalize
c     the Fock matrix "f" to get the MOs, then use MOs to form the
c     next "f" matrix assuming zero differential overlap except for
c     one-center exchange repulsions
c
         iter = 0
         delta = 2.0d0 * converge
         do while (delta.gt.converge .and. iter.lt.maxiter)
            iter = iter + 1
            call jacobi (norbit,fock,en,v)
            do i = 1, norbit
               do j = i, norbit
                  s1 = 0.0d0
                  s2 = 0.0d0
                  gij = gamma(i,j)
                  do k = 1, nfill
                     s2 = s2 - v(i,k)*v(j,k)*gij
                     if (i .eq. j) then
                        do m = 1, norbit
                           s1 = s1 + 2.0d0*gamma(i,m)*v(m,k)**2
                        end do
                     end if
                  end do
                  fock(i,j) = s1 + s2 + hc(i,j)
                  fock(j,i) = fock(i,j)
               end do
            end do
c
c     calculate the ground state energy, where "xi" sums the
c     molecular core integrals, "xj" sums the molecular coulomb
c     repulsion integrals, "xk" sums the molecular exchange
c     repulsion integrals, and "xg" sums the nuclear repulsion
c
            xi = 0.0d0
            xj = 0.0d0
            xk = 0.0d0
            xg = 0.0d0
            do i = 1, nfill
               do j = 1, norbit
                  vij = v(j,i)
                  do k = 1, norbit
                     vik = v(k,i)
                     gjk = gamma(j,k)
                     xi = xi + 2.0d0*vij*vik*hc(j,k)
                     do m = 1, nfill
                        vmj = v(j,m)
                        vmk = v(k,m)
                        xj = xj + 2.0d0*vij*vij*vmk*vmk*gjk
                        xk = xk - vij*vmj*vik*vmk*gjk
                     end do
                  end do
               end do
            end do
            do i = 1, norbit-1
               iorb = iorbit(i)
               qi = qorb(iorb)
               do j = i+1, norbit
                  jorb = iorbit(j)
                  xg = xg + qi*qorb(jorb)*gamma(i,j)
               end do
            end do
            total = xi + xj + xk + xg
            if (iter .ne. 1)  delta = abs(total-totold)
            totold = total
         end do
c
c     print warning if SCF-MO iteration did not converge
c
         if (delta .gt. converge) then
            write (iout,10)
   10       format (' PISCF  --  The SCF Molecular Orbitals have',
     &                 ' Failed to Converge')
         end if
c
c     calculate electron densities from filled MO's
c
         do i = 1, norbit
            do j = 1, norbit
               ed(i,j) = 0.0d0
               do k = 1, nfill
                  ed(i,j) = ed(i,j) + 2.0d0*v(i,k)*v(j,k)
               end do
            end do
         end do
c
c     print out results for the SCF computation
c
         if (debug) then
            if (mode .eq. 'PLANAR') then
               write (iout,20)
   20          format (/,' SCF-MO Calculation for Planar System :')
            else
               write (iout,30)
   30          format (/,' SCF-MO Calculation for Non-Planar',
     &                    ' System :')
            end if
            write (iout,40)  total,norbit,delta,iter
   40       format (/,' Total Energy',11x,f12.4,
     &              /,' Number of Orbitals',5x,i12,
     &              /,' Convergence',12x,d12.4,
     &              /,' Iterations',13x,i12)
            write (iout,50)  xi,xj,xk,xg
   50       format (/,' Core Integrals',9x,f12.4,
     &              /,' Coulomb Repulsion',6x,f12.4,
     &              /,' Exchange Repulsion',5x,f12.4,
     &              /,' Nuclear Repulsion',6x,f12.4)
            write (iout,60)
   60       format (/,' Orbital Energies')
            write (iout,70)  (en(i),i=1,norbit)
   70       format (8f9.4)
            write (iout,80)
   80       format (/,' Molecular Orbitals')
            do i = 1, norbit
               write (iout,90)  (v(i,j),j=1,norbit)
   90          format (8f9.4)
            end do
            write (iout,100)
  100       format (/,' Fock Matrix')
            do i = 1, norbit
               write (iout,110)  (fock(i,j),j=1,norbit)
  110          format (8f9.4)
            end do
            write (iout,120)
  120       format (/,' Density Matrix')
            do i = 1, norbit
               write (iout,130)  (ed(i,j),j=1,norbit)
  130          format (8f9.4)
            end do
            write (iout,140)
  140       format (/,' H-Core Matrix')
            do i = 1, norbit
               write (iout,150)  (hc(i,j),j=1,norbit)
  150          format (8f9.4)
            end do
            write (iout,160)
  160       format (/,' Gamma Matrix')
            do i = 1, norbit
               write (iout,170)  (gamma(i,j),j=1,norbit)
  170          format (8f9.4)
            end do
         end if
c
c     now, get the bond orders (compute p and p*b)
c
         if (debug) then
            write (iout,180)
  180       format (/,' Pisystem Bond Orders')
         end if
         do k = 1, nbpi
            i = ibpi(2,k)
            j = ibpi(3,k)
            p = 0.0d0
            do m = 1, nfill
               p = p + 2.0d0*v(i,m)*v(j,m)
            end do
            if (mode .eq. 'PLANAR') then
               pbpl(k) = p * hc(i,j)/ebeta
            else if (mode .eq. 'NONPLN') then
               pnpl(k) = p
            end if
            if (debug) then
               i = ibnd(1,ibpi(1,k))
               j = ibnd(2,ibpi(1,k))
               write (iout,190)  i,j,p
  190          format (3x,2i6,2x,f10.4)
            end if
         end do
c
c     if we have done planar calculation, do the nonplanar;
c     when both are complete, alter the pisystem constants
c
         if (mode .eq. 'PLANAR') then
            mode = 'NONPLN'
         else if (mode .eq. 'NONPLN') then
            mode = '      '
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (povlap)
      deallocate (en)
      deallocate (ip)
      deallocate (fock)
      deallocate (hc)
      deallocate (v)
      deallocate (gamma)
      deallocate (ed)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine pitilt  --  direction cosines for pisystem  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "pitilt" calculates for each pibond the ratio of the
c     actual p-orbital overlap integral to the ideal overlap
c     if the same orbitals were perfectly parallel
c
c
      subroutine pitilt (povlap)
      use sizes
      use atomid
      use atoms
      use couple
      use piorbs
      implicit none
      integer i,j,k,m
      integer iorb,jorb
      integer list(8)
      real*8 ideal,cosine,rnorm
      real*8 xij,yij,zij,rij
      real*8 a1,b1,c1,a2,b2,c2
      real*8 x2,y2,z2,x3,y3,z3
      real*8 xr(8),yr(8),zr(8)
      real*8 povlap(*)
c
c
c     planes defining each p-orbital are in "piperp"; transform
c     coordinates of "iorb", "jorb" and their associated planes
c     to put "iorb" at origin and "jorb" along the x-axis
c
      do k = 1, nbpi
         i = ibpi(2,k)
         j = ibpi(3,k)
         iorb = iorbit(i)
         jorb = iorbit(j)
         list(1) = iorb
         list(2) = jorb
         do m = 1, 3
            list(m+2) = piperp(m,iorb)
            list(m+5) = piperp(m,jorb)
         end do
         call pimove (list,xr,yr,zr)
c
c     check for sp-hybridized carbon in current bond;
c     assume perfect overlap for any such pibond
c
         if ((atomic(iorb).eq.6 .and. n12(iorb).eq.2) .or.
     &       (atomic(jorb).eq.6 .and. n12(jorb).eq.2)) then
            povlap(k) = 1.0d0
c
c     find and normalize a vector parallel to first p-orbital
c
         else
            x2 = xr(4) - xr(3)
            y2 = yr(4) - yr(3)
            z2 = zr(4) - zr(3)
            x3 = xr(5) - xr(3)
            y3 = yr(5) - yr(3)
            z3 = zr(5) - zr(3)
            a1 = y2*z3 - y3*z2
            b1 = x3*z2 - x2*z3
            c1 = x2*y3 - x3*y2
            rnorm = sqrt(a1*a1+b1*b1+c1*c1)
            a1 = a1 / rnorm
            b1 = b1 / rnorm
            c1 = c1 / rnorm
c
c     now find vector parallel to the second p-orbital,
c     "a2" changes sign to correspond to internuclear axis
c
            x2 = xr(7) - xr(6)
            y2 = yr(7) - yr(6)
            z2 = zr(7) - zr(6)
            x3 = xr(8) - xr(6)
            y3 = yr(8) - yr(6)
            z3 = zr(8) - zr(6)
            a2 = y2*z3 - y3*z2
            b2 = x3*z2 - x2*z3
            c2 = x2*y3 - x3*y2
            rnorm = sqrt(a2*a2+b2*b2+c2*c2)
            a2 = -a2 / rnorm
            b2 = b2 / rnorm
            c2 = c2 / rnorm
c
c     compute the cosine of the angle between p-orbitals;
c     if more than 90 degrees, reverse one of the vectors
c
            cosine = a1*a2 + b1*b2 + c1*c2
            if (cosine .lt. 0.0d0) then
               a2 = -a2
               b2 = -b2
               c2 = -c2
            end if
c
c     find overlap if the orbitals were perfectly parallel
c
            xij = x(iorb) - x(jorb)
            yij = y(iorb) - y(jorb)
            zij = z(iorb) - z(jorb)
            rij = sqrt(xij**2 + yij**2 + zij**2)
            call overlap (atomic(iorb),atomic(jorb),rij,ideal)
c
c     set ratio of actual to ideal overlap for current pibond
c
            povlap(k) = ideal*a1*a2 + b1*b2 + c1*c2
         end if
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine pimove  --  transform pisystem bond vectors  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "pimove" rotates the vector between atoms "list(1)" and
c     "list(2)" so that atom 1 is at the origin and atom 2 along
c     the x-axis; the atoms defining the respective planes are
c     also moved and their bond lengths normalized
c
c
      subroutine pimove (list,xr,yr,zr)
      use sizes
      use atoms
      implicit none
      integer i,j,list(8)
      real*8 xt,yt,zt
      real*8 denom,xold
      real*8 sine,cosine
      real*8 xr(8),yr(8),zr(8)
c
c
c     translate "list" atoms to place atom 1 at origin
c
      j = list(1)
      xt = x(j)
      yt = y(j)
      zt = z(j)
      do i = 1, 8
         j = list(i)
         xr(i) = x(j) - xt
         yr(i) = y(j) - yt
         zr(i) = z(j) - zt
      end do
c
c     rotate "list" atoms to place atom 2 on the x-axis
c
      denom = sqrt(xr(2)**2 + yr(2)**2)
      if (denom .ne. 0.0d0) then
         sine = yr(2) / denom
         cosine = xr(2) / denom
         do i = 1, 8
            xold = xr(i)
            xr(i) = xr(i)*cosine + yr(i)*sine
            yr(i) = yr(i)*cosine - xold*sine
         end do
      end if
      denom = sqrt(xr(2)**2 + zr(2)**2)
      if (denom .ne. 0.0d0) then
         sine = zr(2) / denom
         cosine = xr(2) / denom
         do i = 1, 8
            xold = xr(i)
            xr(i) = xr(i)*cosine + zr(i)*sine
            zr(i) = zr(i)*cosine - xold*sine
         end do
      end if
c
c     normalize the coordinates of atoms defining the plane for atom 1
c     (ie, make all these atoms have unit length to atom 1) so that the
c     orbital makes equal angles with the atoms rather than simply being
c     perpendicular to the common plane of the atoms
c
      do i = 3, 5
         if (list(i) .ne. list(1)) then
            denom = sqrt(xr(i)**2+yr(i)**2+zr(i)**2)
            xr(i) = xr(i) / denom
            yr(i) = yr(i) / denom
            zr(i) = zr(i) / denom
         end if
      end do
c
c     normalization of plane defining atoms for atom 2; for the
c     x-coordinate we translate back to the origin, normalize
c     and then retranslate back along the x-axis
c
      do i = 6, 8
         if (list(i) .ne. list(2)) then
            denom = sqrt((xr(i)-xr(2))**2+yr(i)**2+zr(i)**2)
            xr(i) = (xr(i)-xr(2))/denom + xr(2)
            yr(i) = yr(i) / denom
            zr(i) = zr(i) / denom
         end if
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine pialter  --  modify parameters for pisystem  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "pialter" modifies bond lengths and force constants according
c     to the "planar" P-P-P bond order values; also alters 2-fold
c     torsional parameters based on the "nonplanar" bond orders
c
c
      subroutine pialter
      use sizes
      use atomid
      use bndstr
      use inform
      use iounit
      use piorbs
      use pistuf
      use tors
      implicit none
      integer i,j,k
      integer ia,ib,ic,id
c
c
c     modify the stretching constants and natural bond lengths
c
      if (debug .and. nbpi.ne.0) then
         write (iout,10)
   10    format (/,' Altered Bond Stretching Parameters',
     &              ' for Pi-System :',
     &           //,' Type',14x,'Atom Names',17x,'Initial',
     &              16x,'Final',/)
      end if
      do i = 1, nbpi
         j = ibpi(1,i)
         ia = ibnd(1,j)
         ib = ibnd(2,j)
         bk(j) = bkpi(j) - kslope(j) * (1.0d0-pnpl(i))
         bl(j) = blpi(j) + lslope(j) * (1.0d0-pnpl(i))
         if (debug) then
            write (iout,20)  ia,name(ia),ib,name(ib),bkpi(j),
     &                       blpi(j),bk(j),bl(j)
   20       format (' Bond',6x,2(i7,'-',a3),6x,
     &                 f9.3,f8.4,2x,'-->',f9.3,f8.4)
         end if
      end do
c
c     modify the 2-fold torsional constants across pibonds
c
      if (debug .and. ntpi.ne.0) then
         write (iout,30)
   30    format (/,' Altered 2-Fold Torsional Parameters',
     &              ' for Pi-System :',
     &           //,' Type',25x,'Atom Names',18x,'Initial',
     &              8x,'Final',/)
      end if
      do i = 1, ntpi
         j = itpi(1,i)
         k = itpi(2,i)
         ia = itors(1,j)
         ib = itors(2,j)
         ic = itors(3,j)
         id = itors(4,j)
         tors2(1,j) = pbpl(k) * torsp2(j)
         if (debug) then
            write (iout,40)  ia,name(ia),ib,name(ib),ic,name(ic),
     &                       id,name(id),torsp2(j),tors2(1,j)
   40       format (' Torsion',3x,4(i7,'-',a3),2x,f8.3,2x,'-->',f8.3)
         end if
      end do
      return
      end
