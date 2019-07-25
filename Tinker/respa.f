c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2011  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine respa  --  r-RESPA molecular dynamics step  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "respa" performs a single multiple time step molecular dynamics
c     step using the reversible reference system propagation algorithm
c     (r-RESPA) via a Verlet core with the potential split into fast-
c     and slow-evolving portions
c
c     literature references:
c
c     D. D. Humphreys, R. A. Friesner and B. J. Berne, "A Multiple-
c     Time-Step Molecular Dynamics Algorithm for Macromolecules",
c     Journal of Physical Chemistry, 98, 6885-6892 (1994)
c
c     X. Qian and T. Schlick, "Efficient Multiple-Time-Step Integrators
c     with Distance-Based Force Splitting for Particle-Mesh-Ewald
c     Molecular Dynamics Simulations", Journal of Chemical Physics,
c     115, 4019-4029 (2001)
c
c
      subroutine respa (istep,dt)
      use sizes
      use atomid
      use atoms
      use freeze
      use ielscf
      use moldyn
      use polar
      use units
      use usage
      use virial
      implicit none
      integer i,j,k
      integer istep
      integer nalt
      real*8 dt,dt_2
      real*8 dta,dta_2
      real*8 epot,etot
      real*8 eksum,eps
      real*8 temp,pres
      real*8 ealt,dalt
      real*8 term
      real*8 ekin(3,3)
      real*8 stress(3,3)
      real*8 viralt(3,3)
      real*8, allocatable :: xold(:)
      real*8, allocatable :: yold(:)
      real*8, allocatable :: zold(:)
      real*8, allocatable :: derivs(:,:)
c
c
c     set some time values for the dynamics integration
c
      eps =  0.00000001d0
      dalt = 0.00025d0
      nalt = int(dt/(dalt+eps)) + 1
      dalt = dble(nalt)
      dt_2 = 0.5d0 * dt
      dta = dt / dalt
      dta_2 = 0.5d0 * dta
c
c     store the current atom positions, then find half-step
c     velocities via velocity Verlet recursion
c
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               v(j,i) = v(j,i) + a(j,i)*dt_2
            end do
         end if
      end do
c
c     initialize virial from fast-evolving potential energy terms
c
      do i = 1, 3
         do j = 1, 3
            viralt(j,i) = 0.0d0
         end do
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (xold(n))
      allocate (yold(n))
      allocate (zold(n))
      allocate (derivs(3,n))
c
c     find fast-evolving velocities and positions via Verlet recursion
c
      do k = 1, nalt
         do i = 1, n
            if (use(i)) then
               do j = 1, 3
                  v(j,i) = v(j,i) + aalt(j,i)*dta_2
               end do
               xold(i) = x(i)
               yold(i) = y(i)
               zold(i) = z(i)
               x(i) = x(i) + v(1,i)*dta
               y(i) = y(i) + v(2,i)*dta
               z(i) = z(i) + v(3,i)*dta
            end if
         end do
         if (use_rattle)  call rattle (dta,xold,yold,zold)
c
c     get the fast-evolving potential energy and atomic forces
c
         call gradfast (ealt,derivs)
c
c     use Newton's second law to get fast-evolving accelerations;
c     update fast-evolving velocities using the Verlet recursion
c
         do i = 1, n
            if (use(i)) then
               do j = 1, 3
                  aalt(j,i) = -convert * derivs(j,i) / mass(i)
                  v(j,i) = v(j,i) + aalt(j,i)*dta_2
               end do
            end if
         end do
         if (use_rattle)  call rattle2 (dta)
c
c     find average virial from fast-evolving potential terms
c
         do i = 1, 3
            do j = 1, 3
               viralt(j,i) = viralt(j,i) + vir(j,i)/dalt
            end do
         end do
      end do
c
c     apply Verlet half-step updates for any auxiliary dipoles
c
      if (use_ielscf) then
         do i = 1, n
            if (use(i)) then
               do j = 1, 3
                  vaux(j,i) = vaux(j,i) + aaux(j,i)*dt_2
                  vpaux(j,i) = vpaux(j,i) + apaux(j,i)*dt_2
                  uaux(j,i) = uaux(j,i) + vaux(j,i)*dt
                  upaux(j,i) = upaux(j,i) + vpaux(j,i)*dt
               end do
            end if
         end do
      end if
c
c     get the slow-evolving potential energy and atomic forces
c
      call gradslow (epot,derivs)
      epot = epot + ealt
c
c     make half-step temperature and pressure corrections
c
      call temper2 (dt,temp)
c     call pressure2 (epot,temp)
c
c     use Newton's second law to get the slow accelerations;
c     find full-step velocities using velocity Verlet recursion
c
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               a(j,i) = -convert * derivs(j,i) / mass(i)
               v(j,i) = v(j,i) + a(j,i)*dt_2
            end do
         end if
      end do
c
c     apply Verlet full-step updates for any auxiliary dipoles
c
      if (use_ielscf) then
         term = 2.0d0 / (dt*dt)
         do i = 1, n
            if (use(i)) then
               do j = 1, 3
                  aaux(j,i) = term * (uind(j,i)-uaux(j,i))
                  apaux(j,i) = term * (uinp(j,i)-upaux(j,i))
                  vaux(j,i) = vaux(j,i) + aaux(j,i)*dt_2
                  vpaux(j,i) = vpaux(j,i) + apaux(j,i)*dt_2
               end do
            end if
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (xold)
      deallocate (yold)
      deallocate (zold)
      deallocate (derivs)
c
c     find the constraint-corrected full-step velocities
c
      if (use_rattle)  call rattle2 (dt)
c
c     increment total virial from sum of fast and slow parts
c
      do i = 1, 3
         do j = 1, 3
            vir(j,i) = vir(j,i) + viralt(j,i)
         end do
      end do
c
c     make full-step temperature and pressure corrections
c
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
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine gradfast  --  fast energy & gradient components  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "gradfast" calculates the potential energy and first derivatives
c     for the fast-evolving local valence potential energy terms
c
c
      subroutine gradfast (energy,derivs)
      use limits
      use potent
      implicit none
      real*8 energy
      real*8 derivs(3,*)
      logical save_vdw,save_charge
      logical save_ct
      logical save_chgdpl,save_dipole
      logical save_mpole,save_polar
      logical save_rxnfld,save_solv
      logical save_list
c
c
c     save the original state of slow-evolving potentials
c
      save_vdw = use_vdw
      save_ct = use_ct
      save_charge = use_charge
      save_chgdpl = use_chgdpl
      save_dipole = use_dipole
      save_mpole = use_mpole
      save_polar = use_polar
      save_rxnfld = use_rxnfld
      save_solv = use_solv
      save_list = use_list
c
c     turn off slow-evolving nonbonded potential energy terms
c
      use_vdw = .false.
      use_ct = .false.
      use_charge = .false.
      use_chgdpl = .false.
      use_dipole = .false.
      use_mpole = .false.
      use_polar = .false.
      use_rxnfld = .false.
      use_solv = .false.
      use_list = .false.
c
c     get energy and gradient for fast-evolving potential terms
c
      call gradient (energy,derivs)
c
c     restore the original state of slow-evolving potentials
c
      use_vdw = save_vdw
      use_ct = save_ct
      use_charge = save_charge
      use_chgdpl = save_chgdpl
      use_dipole = save_dipole
      use_mpole = save_mpole
      use_polar = save_polar
      use_rxnfld = save_rxnfld
      use_solv = save_solv
      use_list = save_list
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine gradslow  --  slow energy & gradient components  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "gradslow" calculates the potential energy and first derivatives
c     for the slow-evolving nonbonded potential energy terms
c
c
      subroutine gradslow (energy,derivs)
      use potent
      implicit none
      real*8 energy
      real*8 derivs(3,*)
      logical save_bond,save_angle
      logical save_strbnd,save_urey
      logical save_angang,save_opbend
      logical save_opdist,save_improp
      logical save_imptor,save_tors
      logical save_pitors,save_strtor
      logical save_tortor,save_geom
      logical save_metal,save_extra
c
c
c     save the original state of fast-evolving potentials
c
      save_bond = use_bond
      save_angle = use_angle
      save_strbnd = use_strbnd
      save_urey = use_urey
      save_angang = use_angang
      save_opbend = use_opbend
      save_opdist = use_opdist
      save_improp = use_improp
      save_imptor = use_imptor
      save_tors = use_tors
      save_pitors = use_pitors
      save_strtor = use_strtor
      save_tortor = use_tortor
      save_geom = use_geom
      save_metal = use_metal
      save_extra = use_extra
c
c     turn off fast-evolving valence potential energy terms
c
      use_bond = .false.
      use_angle = .false.
      use_strbnd = .false.
      use_urey = .false.
      use_angang = .false.
      use_opbend = .false.
      use_opdist = .false.
      use_improp = .false.
      use_imptor = .false.
      use_tors = .false.
      use_pitors = .false.
      use_strtor = .false.
      use_tortor = .false.
      use_geom = .false.
      use_metal = .false.
      use_extra = .false.
c
c     get energy and gradient for slow-evolving potential terms
c
      call gradient (energy,derivs)
c
c     restore the original state of fast-evolving potentials
c
      use_bond = save_bond
      use_angle = save_angle
      use_strbnd = save_strbnd
      use_urey = save_urey
      use_angang = save_angang
      use_opbend = save_opbend
      use_opdist = save_opdist
      use_improp = save_improp
      use_imptor = save_imptor
      use_tors = save_tors
      use_pitors = save_pitors
      use_strtor = save_strtor
      use_tortor = save_tortor
      use_geom = save_geom
      use_metal = save_metal
      use_extra = save_extra
      return
      end
