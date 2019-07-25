c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine mdsave  --  save trajectory and restart files  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "mdsave" writes molecular dynamics trajectory snapshots and
c     auxiliary files with velocity, force or induced dipole data;
c     also checks for user requested termination of a simulation
c
c
      !This file is hacked to print out much more
      !information for use in dynamic_omm run
      !Sep/2018
      !Chengwen Liu

      subroutine mdsave (istep,dt,epot,eksum)
      use sizes
      use atomid
      use atoms
      use bound
      use boxes
      use files
      use group
      use inform
      use iounit
      use mdstuf
      use molcul
      use moldyn
      use mpole
      use output
      use polar
      use potent
      use rgddyn
      use socket
      use titles
      use units
      implicit none
      integer i,j,k,istep
      integer ixyz,iind
      integer ivel,ifrc
      integer iend,idump,lext
      integer freeunit,trimtext
      integer moddump
      real*8 dt,epot,eksum
      real*8 pico,wt
      real*8 volboxnow,dens
      logical exist
      character*7 ext
      character*240 endfile
      character*240 xyzfile
      character*240 velfile
      character*240 frcfile
      character*240 indfile
c
c
c     send data via external socket communication if desired
c
      if (.not.sktstart .or. use_socket)  call sktdyn (istep,dt,epot)
c
c     check number of steps between trajectory file dumps
c
      moddump = mod(istep,iwrite)
      if (moddump .ne. 0)  return
c
c     get the sequence number of the current trajectory frame
c
      idump = nprior + istep/iwrite
      lext = 3
      call numeral (idump,ext,lext)
c
c     print header for the instantaneous values at current step
c
      pico = dble(istep) * dt
      write (iout,10)  istep
   10 format (/,' Instantaneous Values for Frame Saved at',
     &           i10,' Dynamics Steps')
c
c     print the current time, potential and kinetic energies
c
      if (digits .ge. 8) then
         write (iout,20)  pico
   20    format (/,' Current Time',8x,f19.8,' Picosecond')
         write (iout,30)  epot
   30    format (' Current Potential',3x,f19.8,' Kcal/mole')
         write (iout,40)  eksum
   40    format (' Current Kinetic',5x,f19.8,' Kcal/mole')
      else if (digits .ge. 6) then
         write (iout,50)  pico
   50    format (/,' Current Time',8x,f17.6,' Picosecond')
         write (iout,60)  epot
   60    format (' Current Potential',3x,f17.6,' Kcal/mole')
         write (iout,70)  eksum
   70    format (' Current Kinetic',5x,f17.6,' Kcal/mole')
      else
         write (iout,80)  pico
   80    format (/,' Current Time',8x,f15.4,' Picosecond')
         write (iout,90)  epot
   90    format (' Current Potential',3x,f15.4,' Kcal/mole')
         write (iout,100)  eksum
  100    format (' Current Kinetic',5x,f15.4,' Kcal/mole')
      end if
c
c     print the values of the lattice lengths and angles
c
      if (use_bounds) then
         if (digits .le. 6) then
            write (iout,110)  xbox,ybox,zbox
  110       format (' Lattice Lengths',6x,3f14.6)
            write (iout,120)  alpha,beta,gamma
  120       format (' Lattice Angles',7x,3f14.6)
         else if (digits .le. 8) then
            write (iout,130)  xbox,ybox,zbox
  130       format (' Lattice Lengths',6x,3f16.8)
            write (iout,140)  alpha,beta,gamma
  140       format (' Lattice Angles',7x,3f16.8)
         else
            write (iout,150)  xbox,ybox,zbox
  150       format (' Lattice Lengths',6x,3f18.10)
            write (iout,160)  alpha,beta,gamma
  160       format (' Lattice Angles',7x,3f18.10)
         end if
      end if
c
c     save coordinates to an archive or numbered structure file
c
      ixyz = freeunit ()
      if (archive) then
         xyzfile = filename(1:leng)
         call suffix (xyzfile,'arc','old')
         inquire (file=xyzfile,exist=exist)
         if (exist) then
            call openend (ixyz,xyzfile)
         else
            open (unit=ixyz,file=xyzfile,status='new')
         end if
      else
         xyzfile = filename(1:leng)//'.'//ext(1:lext)
         call version (xyzfile,'new')
         open (unit=ixyz,file=xyzfile,status='new')
      end if
      call prtxyz (ixyz)
      close (unit=ixyz)
      write (iout,170)  idump
  170 format (' Frame Number',13x,i10)
      write (iout,180)  xyzfile(1:trimtext(xyzfile))
  180 format (' Coordinate File',12x,a)
c
c     Record the current density 
c
      volboxnow = xbox*ybox*zbox 
      dens = (1.0d24/volboxnow) * (totmass/avogadro)
      if (digits .ge. 8) then
         write (iout,245)  dens
  245    format (' Current Density',3x,f19.8,' Grams/cc')
      else if (digits .ge. 6) then
         write (iout,255)  dens
  255    format (' Current Density',3x,f17.6,' Grams/cc')
      else
         write (iout,265)  dens
  265    format (' Current Density',3x,f15.4,' Grams/cc')
      end if

c
c     update the information needed to restart the trajectory
c
      call prtdyn
c
c     save the velocity vector components at the current step
c
      if (velsave) then
         ivel = freeunit ()
         if (archive) then
            velfile = filename(1:leng)
            call suffix (velfile,'vel','old')
            inquire (file=velfile,exist=exist)
            if (exist) then
               call openend (ivel,velfile)
            else
               open (unit=ivel,file=velfile,status='new')
            end if
         else
            velfile = filename(1:leng)//'.'//ext(1:lext)//'v'
            call version (velfile,'new')
            open (unit=ivel,file=velfile,status='new')
         end if
         if (integrate .eq. 'RIGIDBODY') then
            write (ivel,190)  ngrp,title(1:ltitle)
  190       format (i6,2x,a)
            do i = 1, ngrp
               write (ivel,200)  i,(vcm(j,i),j=1,3)
  200          format (i6,3x,d13.6,3x,d13.6,3x,d13.6)
               write (ivel,210)  i,(wcm(j,i),j=1,3)
  210          format (i6,3x,d13.6,3x,d13.6,3x,d13.6)
            end do
         else
            write (ivel,220)  n,title(1:ltitle)
  220       format (i6,2x,a)
            do i = 1, n
               write (ivel,230)  i,name(i),(v(j,i),j=1,3)
  230          format (i6,2x,a3,3x,d13.6,3x,d13.6,3x,d13.6)
            end do
         end if
         close (unit=ivel)
         write (iout,240)  velfile(1:trimtext(velfile))
  240    format (' Velocity File',15x,a)
      end if
c
c     save the force vector components for the current step
c
      if (frcsave .and. integrate.ne.'RIGIDBODY') then
         ifrc = freeunit ()
         if (archive) then
            frcfile = filename(1:leng)
            call suffix (frcfile,'frc','old')
            inquire (file=frcfile,exist=exist)
            if (exist) then
               call openend (ifrc,frcfile)
            else
               open (unit=ifrc,file=frcfile,status='new')
            end if
         else
            frcfile = filename(1:leng)//'.'//ext(1:lext)//'f'
            call version (frcfile,'new')
            open (unit=ifrc,file=frcfile,status='new')
         end if
         write (ifrc,250)  n,title(1:ltitle)
  250    format (i6,2x,a)
         do i = 1, n
            wt = mass(i) / convert
            write (ifrc,260)  i,name(i),(wt*a(j,i),j=1,3)
  260       format (i6,2x,a3,3x,d13.6,3x,d13.6,3x,d13.6)
         end do
         close (unit=ifrc)
         write (iout,270)  frcfile(1:trimtext(frcfile))
  270    format (' Force Vector File',11x,a)
      end if
c
c     save the current induced dipole moment at each site
c
      if (uindsave .and. use_polar) then
         iind = freeunit ()
         if (archive) then
            indfile = filename(1:leng)
            call suffix (indfile,'uind','old')
            inquire (file=indfile,exist=exist)
            if (exist) then
               call openend (iind,indfile)
            else
               open (unit=iind,file=indfile,status='new')
            end if
         else
            indfile = filename(1:leng)//'.'//ext(1:lext)//'u'
            call version (indfile,'new')
            open (unit=iind,file=indfile,status='new')
         end if
         write (iind,280)  n,title(1:ltitle)
  280    format (i6,2x,a)
         do i = 1, npole
            if (polarity(i) .ne. 0.0d0) then
               k = ipole(i)
               write (iind,290)  k,name(k),(debye*uind(j,i),j=1,3)
  290          format (i6,2x,a3,3f12.6)
            end if
         end do
         close (unit=iind)
         write (iout,300)  indfile(1:trimtext(indfile))
  300    format (' Induced Dipole File',10x,a)
      end if
c
c     test for requested termination of the dynamics calculation
c
      endfile = 'tinker.end'
      inquire (file=endfile,exist=exist)
      if (.not. exist) then
         endfile = filename(1:leng)//'.end'
         inquire (file=endfile,exist=exist)
         if (exist) then
            iend = freeunit ()
            open (unit=iend,file=endfile,status='old')
            close (unit=iend,status='delete')
         end if
      end if
      if (exist) then
         write (iout,310)
  310    format (/,' MDSAVE  --  Dynamics Calculation Ending',
     &              ' due to User Request')
         call fatal
      end if
c
c     skip an extra line to keep the output formating neat
c
      moddump = mod(istep,iprint)
      if (verbose .and. moddump.ne.0) then
         write (iout,320)
  320    format ()
      end if
      return
      end
