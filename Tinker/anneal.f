c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1996  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  program anneal  --  molecular dynamics simulated annealing  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "anneal" performs a simulated annealing protocol by means of
c     variable temperature molecular dynamics using either linear,
c     exponential or sigmoidal cooling schedules
c
c
      program anneal
      use sizes
      use atomid
      use atoms
      use bath
      use bndstr
      use bound
      use inform
      use iounit
      use mdstuf
      use potent
      use solute
      use usage
      use warp
      implicit none
      integer i,nequil,next
      integer nstep,istep
      real*8 logmass,factor
      real*8 ratio,sigmoid
      real*8 dt,dtdump
      real*8 hot,cold
      real*8 fuzzy,sharp
      real*8 loose,tight
      logical exist
      character*1 answer
      character*8 cooltyp
      character*240 record
      character*240 string
c
c
c     set up the structure and mechanics calculation
c
      call initial
      call getxyz
      call mechanic
c
c     get choice of statistical ensemble for periodic system
c
      hot = -1.0d0
      cold = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=10,end=10)  hot
      call nextarg (string,exist)
      if (exist)  read (string,*,err=10,end=10)  cold
   10 continue
      do while (hot.lt.0.0d0 .or. cold.lt.0.0d0)
         hot = -1.0d0
         cold = -1.0d0
         write (iout,20)
   20    format (/,' Enter the Initial and Final Temperatures in',
     &              ' Degrees K [1000,0] :  ',$)
         read (input,30)  record
   30    format (a240)
         read (record,*,err=40,end=40)  hot,cold
   40    continue
         if (hot .le. 0.0d0)  hot = 1000.0d0
         if (cold .le. 0.0d0)  cold = 0.0d0
      end do
c
c     set the number of steps of initial equilibration
c
      nequil = -1
      call nextarg (string,exist)
      if (exist)  read (string,*,err=50,end=50)  nequil
   50 continue
      do while (nequil .lt. 0)
         write (iout,60)
   60    format (/,' Enter the Number of Equilibration Steps [0] :  ',$)
         read (input,70,err=80)  nequil
   70    format (i10)
         if (nequil .lt. 0)  nequil = 0
   80    continue
      end do
c
c     set the number of dynamics steps for the cooling protocol
c
      nstep = -1
      call nextarg (string,exist)
      if (exist)  read (string,*,err=90,end=90)  nstep
   90 continue
      do while (nstep .lt. 0)
         write (iout,100)
  100    format (/,' Enter the Number of Cooling Protocol Steps',
     &              ' [2000] :  ',$)
         read (input,110,err=120)  nstep
  110    format (i10)
         if (nstep .lt. 0)  nstep = 2000
  120    continue
      end do
c
c     decide which annealing cooling protocol to use
c
      cooltyp = 'LINEAR'
      call nextarg (answer,exist)
      if (.not. exist) then
         write (iout,130)
  130    format (/,' Use Linear, Sigmoidal or Exponential Cooling',
     &              ' Protocol ([L], S or E) :  ',$)
         read (input,140)  record
  140    format (a240)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .eq. 'S')  cooltyp = 'SIGMOID'
      if (answer .eq. 'E')  cooltyp = 'EXPONENT'
c
c     get the length of the dynamics time step in picoseconds
c
      dt = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=150,end=150)  dt
  150 continue
      do while (dt .lt. 0.0d0)
         write (iout,160)
  160    format (/,' Enter the Time Step Length in Femtoseconds',
     &              ' [1.0] :  ',$)
         read (input,170,err=180)  dt
  170    format (f20.0)
         if (dt .le. 0.0d0)  dt = 1.0d0
  180    continue
      end do
      dt = 0.001d0 * dt
c
c     set the time between trajectory snapshot coordinate dumps
c
      dtdump = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=190,end=190)  dtdump
  190 continue
      do while (dtdump .lt. 0.0d0)
         write (iout,200)
  200    format (/,' Enter Time between Dumps in Picoseconds',
     &              ' [0.1] :  ',$)
         read (input,210,err=220)  dtdump
  210    format (f20.0)
         if (dtdump .le. 0.0d0)  dtdump = 0.1d0
  220    continue
      end do
      iwrite = nint(dtdump/dt)
c
c     get factor by which atomic weights are to be increased
c
      logmass = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=230,end=230)  logmass
  230 continue
      do while (logmass .lt. 0.0d0)
         write (iout,240)
  240    format (/,' Increase Atomic Weights by a Factor of',
     &              ' 10^x [x=0.0] :  ',$)
         read (input,250,err=260)  logmass
  250    format (f20.0)
         if (logmass .le. 0.0d0)  logmass = 0.0d0
  260    continue
      end do
      if (logmass .gt. 0.0d0) then
         factor = 10.0d0**(logmass)
         do i = 1, n
            mass(i) = mass(i) * factor
         end do
      end if
c
c     rate of deformation change for potential surface smoothing
c
      if (use_smooth) then
         sharp = -1.0d0
         call nextarg (string,exist)
         if (exist)  read (string,*,err=270,end=270)  sharp
  270    continue
         do while (sharp .lt. 0.0d0)
            write (iout,280)
  280       format (/,' Enter Final Desired Deformation Parameter',
     &                 ' [0.0] :  ',$)
            read (input,290,err=300)  sharp
  290       format (f20.0)
            if (sharp .le. 0.0d0)  sharp = 0.0d0
  300       continue
         end do
         fuzzy = deform - sharp
         if (fuzzy .le. 0.0d0)  fuzzy = 0.0d0
      end if
c
c     set values for temperature, pressure and coupling baths
c
      isothermal = .true.
      isobaric = .false.
      loose = 100.0d0 * dt
      tight = 10.0d0 * dt
      kelvin = hot
      tautemp = loose
c
c     initialize any holonomic constraints and setup dynamics
c
      call shakeup
      call mdinit
c
c     print out a header lines for the equilibration phase
c
      if (nequil .ne. 0) then
         write (iout,310)
  310    format (/,' Simulated Annealing Equilibration Phase')
         write (iout,320)  nequil,dt,logmass,hot,hot
  320    format (/,' Steps:',i6,3x,'Time/Step:',f6.3,' ps',3x,
     &              'LogMass:',f5.2,3x,'Temp:',f7.1,' to',f7.1)
         flush (iout)
      end if
c
c     take the dynamics steps for the equilibration phase
c
      do istep = 1, nequil
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
c     start the cooling phase from the end of equilibration phase
c
      if (nequil .ne. 0)  call mdinit
c
c     print out a header lines for the cooling protocol
c
      write (iout,330)
  330 format (/,' Simulated Annealing Cooling Protocol')
      write (iout,340)  nstep,dt,logmass,hot,cold
  340 format (/,' Steps:',i6,3x,'Time/Step:',f6.3,' ps',3x,
     &           'LogMass:',f5.2,3x,'Temp:',f7.1,' to',f7.1)
      flush (iout)
c
c     set target temperature using the desired cooling protocol
c
      do istep = 1, nstep
         ratio = dble(istep) / dble(nstep)
         if (cooltyp .eq. 'SIGMOID') then
            ratio = sigmoid (3.5d0,ratio)
         else if (cooltyp .eq. 'EXPONENT') then
            ratio = 1.0d0 - exp(-5.0d0*ratio)
         end if
         kelvin = hot*(1.0d0-ratio) + cold*ratio
         tautemp = loose*(1.0d0-ratio) + tight*ratio
c
c     set the deformation value if potential smoothing is used
c
         if (use_smooth) then
            ratio = (1.0d0-dble(istep)/dble(nstep))**3
            deform = sharp + ratio*fuzzy
         end if
c
c     integrate equations of motion to take a time step
c
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
