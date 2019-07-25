c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2010 by Teresa Head-Gordon & Jay W. Ponder  ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine bussi  --  Bussi NPT molecular dynamics step  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "bussi" performs a single molecular dynamics time step via
c     the Bussi-Parrinello isothermal-isobaric algorithm
c
c     literature reference:
c
c     G. Bussi, T. Zykova-Timan and M. Parrinello, "Isothermal-Isobaric
c     Molecular Dynamics using Stochastic Velocity Rescaling", Journal
c     of Chemical Physics, 130, 074101 (2009)
c
c     original version written by Teresa Head-Gordon, October 2010
c
c
      subroutine bussi (istep,dt)
      use sizes
      use atomid
      use atoms
      use bath
      use boxes
      use freeze
      use ielscf
      use mdstuf
      use moldyn
      use polar
      use units
      use usage
      implicit none
      integer i,j,istep
      real*8 dt,dt_2,dt_x
      real*8 dt2_2,dt3_2
      real*8 epot,etot,eksum
      real*8 expterm,sinhterm
      real*8 kt,w,temp,pres
      real*8 part1,part2
      real*8 factor,term
      real*8 ekin(3,3)
      real*8 stress(3,3)
      real*8, allocatable :: xold(:)
      real*8, allocatable :: yold(:)
      real*8, allocatable :: zold(:)
      real*8, allocatable :: derivs(:,:)
c
c
c     set some time values, constants and barostat mass
c
      dt_2 = 0.5d0 * dt
      dt2_2 = dt_2 * dt_2
      dt3_2 = dt2_2 * dt_2
      kt = boltzmann * kelvin
      w = dble(nfree) * kt * taupres * taupres
c
c     get Beeman integration coefficients for velocity updates
c
      factor = dble(bmnmix)
      dt_x = dt / factor
      part1 = 0.5d0*factor + 1.0d0
      part2 = part1 - 2.0d0
c
c     make half-step temperature correction and get pressure
c
      call temper (dt_2,eksum,ekin,temp)
      call pressure (dt,epot,ekin,temp,pres,stress)
c
c     get half-step Beeman velocities and update barostat velocity
c
      eta = eta + 3.0d0*(volbox*(pres-atmsph)*convert/prescon
     &                         + 2.0*kt)*dt_2/w
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               eta = eta + mass(i)*a(j,i)*v(j,i)*dt2_2/w
     &                   + mass(i)*a(j,i)*a(j,i)*dt3_2/(3.0d0*w)
               v(j,i) = v(j,i) + (part1*a(j,i)-aalt(j,i))*dt_x
            end do
        end if
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (xold(n))
      allocate (yold(n))
      allocate (zold(n))
      allocate (derivs(3,n))
c
c     store the current atom positions, then alter positions
c     and velocities via coupling to the barostat
c
      term = eta * dt
      expterm = exp(term)
      sinhterm = sinh(term)
      do i = 1, n
         if (use(i)) then
            xold(i) = x(i)
            yold(i) = y(i)
            zold(i) = z(i)
            x(i) = x(i)*expterm + v(1,i)*sinhterm/eta
            y(i) = y(i)*expterm + v(2,i)*sinhterm/eta
            z(i) = z(i)*expterm + v(3,i)*sinhterm/eta
            do j = 1, 3
               v(j,i) = v(j,i) / expterm
            end do
         end if
      end do
c
c     set the new box dimensions and other lattice values;
c     current version assumes isotropic pressure
c
      xbox = xbox * expterm
      ybox = ybox * expterm
      zbox = zbox * expterm
      call lattice
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
         call temper2 (dt,temp)
      end if
c
c     get constraint-corrected positions and half-step velocities
c
      if (use_rattle)  call rattle (dt,xold,yold,zold)
c
c     get the potential energy and atomic forces
c
      call gradient (epot,derivs)
c
c     use Newton's second law to get the next accelerations
c
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               aalt(j,i) = a(j,i)
               a(j,i) = -convert * derivs(j,i) / mass(i)
            end do
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xold)
      deallocate (yold)
      deallocate (zold)
      deallocate (derivs)
c
c     get full-step Beeman velocities and update barostat velocity
c
      eta = eta + 3.0d0*(volbox*(pres-atmsph)*convert/prescon
     &                         + 2.0*kt)*dt_2/w
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               eta = eta + mass(i)*a(j,i)*v(j,i)*dt2_2/w
     &                   + mass(i)*a(j,i)*a(j,i)*dt3_2/(3.0d0*w)
               v(j,i) = v(j,i) + (part2*a(j,i)+aalt(j,i))*dt_x
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
c     find the constraint-corrected full-step velocities
c
      if (use_rattle)  call rattle2 (dt)
c
c     make full-step temperature correction and get pressure
c
      call temper (dt_2,eksum,ekin,temp)
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
