c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine mdinit  --  initialize a dynamics trajectory  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "mdinit" initializes the velocities and accelerations
c     for a molecular dynamics trajectory, including restarts
c
c
      subroutine mdinit
      use sizes
      use atomid
      use atoms
      use bath
      use bound
      use couple
      use files
      use freeze
      use group
      use ielscf
      use inform
      use iounit
      use keys
      use mdstuf
      use molcul
      use moldyn
      use mpole
      use rgddyn
      use rigid
      use stodyn
      use units
      use uprior
      use usage
      implicit none
      integer i,j,k
      integer idyn,lext
      integer size,next
      integer freeunit
      real*8 e,ekt,qterm
      real*8 maxwell,speed
      real*8 amass,gmass
      real*8 vec(3)
      real*8, allocatable :: derivs(:,:)
      logical exist
      character*7 ext
      character*20 keyword
      character*240 dynfile
      character*240 record
      character*240 string
c
c
c     set default parameters for the dynamics trajectory
c
      integrate = 'BEEMAN'
      bmnmix = 8
      nfree = 0
      irest = 1
      velsave = .false.
      frcsave = .false.
      uindsave = .false.
      friction = 91.0d0
      use_sdarea = .false.
      use_pred = .false.
      polpred = 'LSQR'
      use_ielscf = .false.
      iprint = 100
c
c     set default values for temperature and pressure control
c
      thermostat = 'BUSSI'
      tautemp = 0.2d0
      collide = 0.1d0
      do i = 1, maxnose
         vnh(i) = 0.0d0
         qnh(i) = 0.0d0
         gnh(i) = 0.0d0
      end do
      barostat = 'BERENDSEN'
      anisotrop = .false.
      taupres = 2.0d0
      compress = 0.000046d0
      vbar = 0.0d0
      qbar = 0.0d0
      gbar = 0.0d0
      eta = 0.0d0
      voltrial = 25
      volmove = 100.0d0
      volscale = 'MOLECULAR'
c
c     check for keywords containing any altered parameters
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:11) .eq. 'INTEGRATOR ') then
            call getword (record,integrate,next)
            call upcase (integrate)
         else if (keyword(1:14) .eq. 'BEEMAN-MIXING ') then
            read (string,*,err=10,end=10)  bmnmix
         else if (keyword(1:16) .eq. 'DEGREES-FREEDOM ') then
            read (string,*,err=10,end=10)  nfree
         else if (keyword(1:15) .eq. 'REMOVE-INERTIA ') then
            read (string,*,err=10,end=10)  irest
         else if (keyword(1:14) .eq. 'SAVE-VELOCITY ') then
            velsave = .true.
         else if (keyword(1:11) .eq. 'SAVE-FORCE ') then
            frcsave = .true.
         else if (keyword(1:13) .eq. 'SAVE-INDUCED ') then
            uindsave = .true.
         else if (keyword(1:9) .eq. 'FRICTION ') then
            read (string,*,err=10,end=10)  friction
         else if (keyword(1:17) .eq. 'FRICTION-SCALING ') then
            use_sdarea = .true.
         else if (keyword(1:14) .eq. 'POLAR-PREDICT ') then
            call getword (record,polpred,next)
            call upcase (polpred)
            if (polpred .eq. 'IEL') then
               use_ielscf = .true.
            else
               use_pred = .true.
            end if
         else if (keyword(1:8) .eq. 'IEL-SCF ') then
            use_ielscf = .true.
         else if (keyword(1:11) .eq. 'THERMOSTAT ') then
            call getword (record,thermostat,next)
            call upcase (thermostat)
         else if (keyword(1:16) .eq. 'TAU-TEMPERATURE ') then
            read (string,*,err=10,end=10)  tautemp
         else if (keyword(1:10) .eq. 'COLLISION ') then
            read (string,*,err=10,end=10)  collide
         else if (keyword(1:9) .eq. 'BAROSTAT ') then
            call getword (record,barostat,next)
            call upcase (barostat)
         else if (keyword(1:15) .eq. 'ANISO-PRESSURE ') then
            anisotrop = .true.
         else if (keyword(1:13) .eq. 'TAU-PRESSURE ') then
            read (string,*,err=10,end=10)  taupres
         else if (keyword(1:9) .eq. 'COMPRESS ') then
            read (string,*,err=10,end=10)  compress
         else if (keyword(1:13) .eq. 'VOLUME-TRIAL ') then
            read (string,*,err=10,end=10)  voltrial
         else if (keyword(1:12) .eq. 'VOLUME-MOVE ') then
            read (string,*,err=10,end=10)  volmove
         else if (keyword(1:13) .eq. 'VOLUME-SCALE ') then
            call getword (record,volscale,next)
            call upcase (volscale)
         else if (keyword(1:9) .eq. 'PRINTOUT ') then
            read (string,*,err=10,end=10)  iprint
         end if
   10    continue
      end do
c
c     make sure all atoms or groups have a nonzero mass
c
      if (integrate .eq. 'RIGIDBODY') then
         do i = 1, ngrp
            if (grpmass(i) .le. 0.0d0) then
               grpmass(i) = 1.0d0
               if (igrp(1,i) .le. igrp(2,i)) then
                  totmass = totmass + 1.0d0
                  write (iout,20)  i
   20             format (/,' MDINIT  --  Warning, Mass of Group',i6,
     &                       ' Set to 1.0 for Dynamics')
               end if
            end if
         end do
      else
         do i = 1, n
            if (use(i) .and. mass(i).le.0.0d0 .and. atomic(i).ne.0) then
               mass(i) = 1.0d0
               totmass = totmass + 1.0d0
               write (iout,30)  i
   30          format (/,' MDINIT  --  Warning, Mass of Atom',i6,
     &                    ' Set to 1.0 for Dynamics')
            end if
         end do
      end if
c
c     perform dynamic allocation of some global arrays
c
      if (use_pred) then
         if (.not. allocated(udalt))  allocate (udalt(maxualt,3,n))
         if (.not. allocated(upalt))  allocate (upalt(maxualt,3,n))
         if (.not. allocated(usalt))  allocate (usalt(maxualt,3,n))
         if (.not. allocated(upsalt))  allocate (upsalt(maxualt,3,n))
c
c     set the Gear predictor binomial coefficients
c
         gear(1) = 6.0d0
         gear(2) = -15.0d0
         gear(3) = 20.0d0
         gear(4) = -15.0d0
         gear(5) = 6.0d0
         gear(6) = -1.0d0
         gear(7) = 0.0d0
c
c     set always stable predictor-corrector (ASPC) coefficients
c
         aspc(1) = 22.0d0 / 7.0d0
         aspc(2) = -55.0d0 / 14.0d0
         aspc(3) = 55.0d0 / 21.0d0
         aspc(4) = -22.0d0 / 21.0d0
         aspc(5) = 5.0d0 / 21.0d0
         aspc(6) = -1.0d0 / 42.0d0
         aspc(7) = 0.0d0
c
c    initialize prior values of induced dipole moments
c
         nualt = 0
         do i = 1, npole
            do j = 1, 3
               do k = 1, maxualt
                  udalt(k,j,i) = 0.0d0
                  upalt(k,j,i) = 0.0d0
                  usalt(k,j,i) = 0.0d0
                  upsalt(k,j,i) = 0.0d0
               end do
            end do
         end do
      end if
c
c     initialize inertial extended Lagrangian method
c
      if (use_ielscf)  call auxinit
c
c     enforce use of velocity Verlet with Andersen thermostat
c
      if (thermostat .eq. 'ANDERSEN') then
         if (integrate .eq. 'BEEMAN')  integrate = 'VERLET'
      end if
c
c     enforce use of Bussi thermostat/barostat with integrator
c
      if (integrate .eq. 'BUSSI') then
         thermostat = 'BUSSI'
         barostat = 'BUSSI'
      else if (thermostat.eq.'BUSSI' .and. barostat.eq.'BUSSI') then
         integrate = 'BUSSI'
      end if
c
c     enforce use of Nose-Hoover thermostat/barostat with integrator
c
      if (integrate .eq. 'NOSE-HOOVER') then
         thermostat = 'NOSE-HOOVER'
         barostat = 'NOSE-HOOVER'
      else if (thermostat.eq.'NOSE-HOOVER' .and.
     &         barostat.eq.'NOSE-HOOVER') then
         integrate = 'NOSE-HOOVER'
      end if
c
c     check for use of Monte Carlo barostat with constraints
c
      if (barostat.eq.'MONTECARLO' .and. volscale.eq.'ATOMIC') then
         if (use_rattle) then
            write (iout,40)
   40       format (/,' MDINIT  --  Atom-based Monte Carlo',
     &                 ' Barostat Incompatible with RATTLE')
            call fatal
         end if
      end if
c
c     perform dynamic allocation of some global arrays
c
      if (integrate .eq. 'RIGIDBODY') then
         if (.not. allocated(xcmo))  allocate (xcmo(n))
         if (.not. allocated(ycmo))  allocate (ycmo(n))
         if (.not. allocated(zcmo))  allocate (zcmo(n))
         if (.not. allocated(vcm))  allocate (vcm(3,ngrp))
         if (.not. allocated(wcm))  allocate (wcm(3,ngrp))
         if (.not. allocated(lm))  allocate (lm(3,ngrp))
         if (.not. allocated(vc))  allocate (vc(3,ngrp))
         if (.not. allocated(wc))  allocate (wc(3,ngrp))
         if (.not. allocated(linear))  allocate (linear(ngrp))
      else
         if (.not. allocated(v))  allocate (v(3,n))
         if (.not. allocated(a))  allocate (a(3,n))
         if (.not. allocated(aalt))  allocate (aalt(3,n))
      end if
c
c     set the number of degrees of freedom for the system
c
      if (nfree .eq. 0) then
         if (integrate .eq. 'RIGIDBODY') then
            call grpline
            nfree = 6 * ngrp
            do i = 1, ngrp
               size = igrp(2,i) - igrp(1,i) + 1
               if (size .eq. 1)  nfree = nfree - 3
               if (linear(i))  nfree = nfree - 1
            end do
         else
            nfree = 3 * nuse
         end if
         if (use_rattle) then
            nfree = nfree - nrat
            do i = 1, nratx
               nfree = nfree - kratx(i)
            end do
         end if
         if (isothermal .and. thermostat.ne.'ANDERSEN'.and.
     &       integrate.ne.'STOCHASTIC' .and. integrate.ne.'GHMC') then
            if (use_bounds) then
               nfree = nfree - 3
            else
               nfree = nfree - 6
            end if
         end if
         if (barostat .eq. 'BUSSI')  nfree = nfree + 1
      end if
c
c     check for a nonzero number of degrees of freedom
c
      if (nfree .lt. 0)  nfree = 0
      if (debug) then
         write (iout,50)  nfree
   50    format (/,' Number of Degrees of Freedom for Dynamics :',i10)
      end if
      if (nfree .eq. 0) then
         write (iout,60)
   60    format (/,' MDINIT  --  No Degrees of Freedom for Dynamics')
         call fatal
      end if
c
c     set masses for Nose-Hoover thermostat and barostat
c
      if (thermostat .eq. 'NOSE-HOOVER') then
         ekt = gasconst * kelvin
         qterm = ekt * tautemp * tautemp
         do j = 1, maxnose
            if (qnh(j) .eq. 0.0d0)  qnh(j) = qterm
         end do
         qnh(1) = dble(nfree) * qnh(1)
      end if
      if (barostat .eq. 'NOSE-HOOVER') then
         ekt = gasconst * kelvin
         qterm = ekt * taupres * taupres
         qbar = dble(nfree+1) * qterm
      end if
c
c     decide whether to remove center of mass motion
c
      dorest = .true.
      if (irest .eq. 0)  dorest = .false.
      if (nuse. ne. n)  dorest = .false.
      if (integrate .eq. 'STOCHASTIC')  dorest = .false.
      if (integrate .eq. 'GHMC')  dorest = .false.
      if (isothermal .and. thermostat.eq.'ANDERSEN')  dorest = .false.
c
c     try to restart using prior velocities and accelerations
c
      dynfile = filename(1:leng)//'.dyn'
      call version (dynfile,'old')
      inquire (file=dynfile,exist=exist)
      if (exist) then
         idyn = freeunit ()
         open (unit=idyn,file=dynfile,status='old')
         rewind (unit=idyn)
         call readdyn (idyn)
         close (unit=idyn)
c
c     set translational velocities for rigid body dynamics
c
      else if (integrate .eq. 'RIGIDBODY') then
         do i = 1, ngrp
            gmass = grpmass(i)
            speed = maxwell (gmass,kelvin)
            call ranvec (vec)
            do j = 1, 3
               vcm(j,i) = speed * vec(j)
               wcm(j,i) = 0.0d0
               lm(j,i) = 0.0d0
            end do
         end do
         if (nuse .eq. n)  call mdrest (0)
c
c     set velocities and fast/slow accelerations for RESPA method
c
      else if (integrate .eq. 'RESPA') then
         allocate (derivs(3,n))
         call gradslow (e,derivs)
         do i = 1, n
            amass = mass(i)
            if (use(i) .and. amass.ne.0.0d0) then
               speed = maxwell (amass,kelvin)
               call ranvec (vec)
               do j = 1, 3
                  v(j,i) = speed * vec(j)
                  a(j,i) = -convert * derivs(j,i) / mass(i)
               end do
            else
               do j = 1, 3
                  v(j,i) = 0.0d0
                  a(j,i) = 0.0d0
               end do
            end if
         end do
         call gradfast (e,derivs)
         do i = 1, n
            amass = mass(i)
            if (use(i) .and. amass.ne.0.0d0) then
               do j = 1, 3
                  aalt(j,i) = -convert * derivs(j,i) / amass
               end do
            else
               do j = 1, 3
                  aalt(j,i) = 0.0d0
               end do
            end if
         end do
         deallocate (derivs)
         if (nuse .eq. n)  call mdrest (0)
c
c     set velocities and accelerations for Cartesian dynamics
c
      else
         allocate (derivs(3,n))
         call gradient (e,derivs)
         do i = 1, n
            amass = mass(i)
            if (use(i) .and. amass.ne.0.0d0) then
               speed = maxwell (amass,kelvin)
               call ranvec (vec)
               do j = 1, 3
                  v(j,i) = speed * vec(j)
                  a(j,i) = -convert * derivs(j,i) / amass
                  aalt(j,i) = a(j,i)
               end do
            else
               do j = 1, 3
                  v(j,i) = 0.0d0
                  a(j,i) = 0.0d0
                  aalt(j,i) = 0.0d0
               end do
            end if
         end do
         deallocate (derivs)
         if (nuse .eq. n)  call mdrest (0)
      end if
c
c     check for any prior dynamics coordinate sets
c
      i = 0
      exist = .true.
      do while (exist)
         i = i + 1
         lext = 3
         call numeral (i,ext,lext)
         dynfile = filename(1:leng)//'.'//ext(1:lext)
         inquire (file=dynfile,exist=exist)
         if (.not.exist .and. i.lt.100) then
            lext = 2
            call numeral (i,ext,lext)
            dynfile = filename(1:leng)//'.'//ext(1:lext)
            inquire (file=dynfile,exist=exist)
         end if
         if (.not.exist .and. i.lt.10) then
            lext = 1
            call numeral (i,ext,lext)
            dynfile = filename(1:leng)//'.'//ext(1:lext)
            inquire (file=dynfile,exist=exist)
         end if
      end do
      nprior = i - 1
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine auxinit  --  initialize auxiliary dipole values  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "auxinit" initializes auxiliary variables and settings for
c     inertial extended Lagrangian induced dipole prediction
c
c     literature reference:
c
c     A. Albaugh, O. Demerdash, and T. Head-Gordon, "An Efficient and
c     Stable Hybrid Extended Lagrangian/Self-Consistent Field Scheme
c     for Solving Classical Mutual Induction", Journal of Chemical
c     Physics, 143, 174104 (2015)
c
c
      subroutine auxinit
      use atomid
      use atoms
      use ielscf
      use keys
      use polar
      implicit none
      integer i,j,next
      real*8 speed
      real*8 weight
      real*8 maxwell
      real*8 vec(3)
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     set defaults for auxiliary thermostat control variables
c
      nfree_aux = 3 * npolar
      kelvin_aux = 100000.0d0
      tautemp_aux = 0.1d0
c
c     check for keywords containing auxiliary thermostat values
c 
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:13) .eq. 'AUX-TAUTEMP ') then
            read (string,*,err=10,end=10)  tautemp_aux
         else if (keyword(1:9) .eq. 'AUX-TEMP ') then
            read (string,*,err=10,end=10)  kelvin_aux
         end if
   10    continue
      end do
c
c     perform dynamic allocation of some global arrays
c 
      allocate (uaux(3,n))
      allocate (vaux(3,n))
      allocate (aaux(3,n))
      allocate (upaux(3,n))
      allocate (vpaux(3,n))
      allocate (apaux(3,n))
c
c     set auxiliary dipole values equal to induced dipoles
c
      use_ielscf = .false.
      call induce
      use_ielscf = .true.
      do i = 1, n
         do j = 1, 3
            uaux(j,i) = uind(j,i)
            upaux(j,i) = uinp(j,i)
         end do
      end do
c
c     set velocities and accelerations for auxiliary dipoles
c
      do i = 1, n
         weight = 1.0d0
         speed = maxwell (weight,kelvin_aux)
         call ranvec (vec)
         do j = 1, 3
            vaux(j,i) = speed * vec(j)
            aaux(j,i) = 0.0d0
            vpaux(j,i) = vaux(j,i)
            apaux(j,i) = 0.0d0
         end do
      end do
      return
      end
