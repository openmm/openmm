c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  program dynamic  --  run molecular or stochastic dynamics  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "dynamic" computes a molecular or stochastic dynamics trajectory
c     in one of the standard statistical mechanical ensembles and using
c     any of several possible integration methods
c
c
      program dynamic
      use sizes
      use atoms
      use bath
      use bndstr
      use bound
      use inform
      use iounit
      use keys
      use mdstuf
      use potent
      use solute
      use stodyn
      use usage
      implicit none
      integer i,istep,nstep
      integer mode,next
      real*8 dt,dtdump
      logical exist
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     set up the structure and molecular mechanics calculation
c
      call initial
      call getxyz
      call mechanic
c
c     initialize the temperature, pressure and coupling baths
c
      kelvin = 0.0d0
      atmsph = 0.0d0
      isothermal = .false.
      isobaric = .false.
c
c     check for keywords containing any altered parameters
c
      integrate = 'BEEMAN'
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:11) .eq. 'INTEGRATOR ') then
            call getword (record,integrate,next)
            call upcase (integrate)
         end if
      end do
c
c     initialize the simulation length as number of time steps
c
      nstep = -1
      call nextarg (string,exist)
      if (exist)  read (string,*,err=10,end=10)  nstep
   10 continue
      do while (nstep .lt. 0)
         write (iout,20)
   20    format (/,' Enter the Number of Dynamics Steps to be',
     &              ' Taken :  ',$)
         read (input,30,err=40)  nstep
   30    format (i10)
         if (nstep .lt. 0)  nstep = 0
   40    continue
      end do
c
c     get the length of the dynamics time step in picoseconds
c
      dt = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=50,end=50)  dt
   50 continue
      do while (dt .lt. 0.0d0)
         write (iout,60)
   60    format (/,' Enter the Time Step Length in Femtoseconds',
     &              ' [1.0] :  ',$)
         read (input,70,err=80)  dt
   70    format (f20.0)
         if (dt .le. 0.0d0)  dt = 1.0d0
   80    continue
      end do
      dt = 0.001d0 * dt
c
c     enforce bounds on thermostat and barostat coupling times
c
      tautemp = max(tautemp,dt)
      taupres = max(taupres,dt)
c
c     set the time between trajectory snapshot coordinate dumps
c
      dtdump = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=90,end=90)  dtdump
   90 continue
      do while (dtdump .lt. 0.0d0)
         write (iout,100)
  100    format (/,' Enter Time between Dumps in Picoseconds',
     &              ' [0.1] :  ',$)
         read (input,110,err=120)  dtdump
  110    format (f20.0)
         if (dtdump .le. 0.0d0)  dtdump = 0.1d0
  120    continue
      end do
      iwrite = nint(dtdump/dt)
c
c     get choice of statistical ensemble for periodic system
c
      if (use_bounds) then
         mode = -1
         call nextarg (string,exist)
         if (exist)  read (string,*,err=130,end=130)  mode
  130    continue
         do while (mode.lt.1 .or. mode.gt.4)
            write (iout,140)
  140       format (/,' Available Statistical Mechanical Ensembles :',
     &              //,4x,'(1) Microcanonical (NVE)',
     &              /,4x,'(2) Canonical (NVT)',
     &              /,4x,'(3) Isoenthalpic-Isobaric (NPH)',
     &              /,4x,'(4) Isothermal-Isobaric (NPT)',
     &              //,' Enter the Number of the Desired Choice',
     &                 ' [1] :  ',$)
            read (input,150,err=160)  mode
  150       format (i10)
            if (mode .le. 0)  mode = 1
  160       continue
         end do
         if (integrate.eq.'BUSSI' .or. integrate.eq.'NOSE-HOOVER'
     &                .or. integrate.eq.'GHMC') then
            if (mode .ne. 4) then
               mode = 4
               write (iout,170)
  170          format (/,' Switching to NPT Ensemble as Required',
     &                    ' by Chosen Integrator')
            end if
         end if
         if (mode.eq.2 .or. mode.eq.4) then
            isothermal = .true.
            kelvin = -1.0d0
            call nextarg (string,exist)
            if (exist)  read (string,*,err=180,end=180)  kelvin
  180       continue
            do while (kelvin .lt. 0.0d0)
               write (iout,190)
  190          format (/,' Enter the Desired Temperature in Degrees',
     &                    ' K [298] :  ',$)
               read (input,200,err=210)  kelvin
  200          format (f20.0)
               if (kelvin .le. 0.0d0)  kelvin = 298.0d0
  210          continue
            end do
         end if
         if (mode.eq.3 .or. mode.eq.4) then
            isobaric = .true.
            atmsph = -1.0d0
            call nextarg (string,exist)
            if (exist)  read (string,*,err=220,end=220)  atmsph
  220       continue
            do while (atmsph .lt. 0.0d0)
               write (iout,230)
  230          format (/,' Enter the Desired Pressure in Atm',
     &                    ' [1.0] :  ',$)
               read (input,240,err=250)  atmsph
  240          format (f20.0)
               if (atmsph .le. 0.0d0)  atmsph = 1.0d0
  250          continue
            end do
         end if
      end if
c
c     use constant energy or temperature for nonperiodic system
c
      if (.not. use_bounds) then
         mode = -1
         call nextarg (string,exist)
         if (exist)  read (string,*,err=260,end=260)  mode
  260    continue
         do while (mode.lt.1 .or. mode.gt.2)
            write (iout,270)
  270       format (/,' Available Simulation Control Modes :',
     &              //,4x,'(1) Constant Total Energy Value (E)',
     &              /,4x,'(2) Constant Temperature via Thermostat (T)',
     &              //,' Enter the Number of the Desired Choice',
     &                 ' [1] :  ',$)
            read (input,280,err=290)  mode
  280       format (i10)
            if (mode .le. 0)  mode = 1
  290       continue
         end do
         if (mode .eq. 2) then
            isothermal = .true.
            kelvin = -1.0d0
            call nextarg (string,exist)
            if (exist)  read (string,*,err=300,end=300)  kelvin
  300       continue
            do while (kelvin .lt. 0.0d0)
               write (iout,310)
  310          format (/,' Enter the Desired Temperature in Degrees',
     &                    ' K [298] :  ',$)
               read (input,320,err=330)  kelvin
  320          format (f20.0)
               if (kelvin .le. 0.0d0)  kelvin = 298.0d0
  330          continue
            end do
         end if
      end if
c
c     initialize any holonomic constraints and setup dynamics
c
      call shakeup
      call mdinit
c
c     print out a header line for the dynamics computation
c
      if (integrate .eq. 'VERLET') then
         write (iout,340)
  340    format (/,' Molecular Dynamics Trajectory via',
     &              ' Velocity Verlet Algorithm')
      else if (integrate .eq. 'STOCHASTIC') then
         write (iout,350)
  350    format (/,' Stochastic Dynamics Trajectory via',
     &              ' Velocity Verlet Algorithm')
      else if (integrate .eq. 'BUSSI') then
         write (iout,360)
  360    format (/,' Molecular Dynamics Trajectory via',
     &              ' Bussi-Parrinello NPT Algorithm')
      else if (integrate .eq. 'NOSE-HOOVER') then
         write (iout,370)
  370    format (/,' Molecular Dynamics Trajectory via',
     &              ' Nose-Hoover NPT Algorithm')
      else if (integrate .eq. 'GHMC') then
         write (iout,380)
  380    format (/,' Stochastic Dynamics Trajectory via',
     &              ' Generalized Hybrid Monte Carlo')
      else if (integrate .eq. 'RIGIDBODY') then
         write (iout,390)
  390    format (/,' Molecular Dynamics Trajectory via',
     &              ' Rigid Body Algorithm')
      else if (integrate .eq. 'RESPA') then
         write (iout,400)
  400    format (/,' Molecular Dynamics Trajectory via',
     &              ' r-RESPA MTS Algorithm')
      else
         write (iout,410)
  410    format (/,' Molecular Dynamics Trajectory via',
     &              ' Modified Beeman Algorithm')
      end if
      flush (iout)
c
c     integrate equations of motion to take a time step
c
      do istep = 1, nstep
         if (integrate .eq. 'VERLET') then
            call verlet (istep,dt)
         else if (integrate .eq. 'STOCHASTIC') then
            call sdstep (istep,dt)
         else if (integrate .eq. 'BUSSI') then
            call bussi (istep,dt)
         else if (integrate .eq. 'NOSE-HOOVER') then
            call nose (istep,dt)
         else if (integrate .eq. 'GHMC') then
            call ghmcstep (istep,dt)
         else if (integrate .eq. 'RIGIDBODY') then
            call rgdstep (istep,dt)
         else if (integrate .eq. 'RESPA') then
            call respa (istep,dt)
         else
            call beeman (istep,dt)
         end if
      end do
c
c     perform any final tasks before program exit
c
      call final
      end
