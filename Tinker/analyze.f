c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  program analyze  --  energy partitioning and analysis  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "analyze" computes and displays the total potential energy;
c     options are provided to display system and force field info,
c     partition the energy by atom or by potential function type,
c     show force field parameters by atom; output the large energy
c     interactions and find electrostatic and inertial properties
c
c
      program analyze
      use sizes
      use atoms
      use files
      use inform
      use iounit
      implicit none
      integer i,j,ixyz
      integer frame,nold
      integer freeunit
      integer trimtext
      integer list(20)
      logical dosystem,doparam
      logical doenergy,doatom
      logical dolarge,dodetail
      logical domoment,dovirial
      logical doconect
      logical exist
      logical, allocatable :: active(:)
      character*1 letter
      character*240 record
      character*240 string
      character*240 xyzfile
c
c
c     set up the structure and mechanics calculation
c
      call initial
      call getxyz
      call mechanic
c
c     get the desired types of analysis to be performed
c
      call nextarg (string,exist)
      if (.not. exist) then
         write (iout,10)
   10    format (/,' The TINKER Energy Analysis Utility Can :',
     &           //,' General System and Force Field Information [G]',
     &           /,' Force Field Parameters for Interactions [P]',
     &           /,' Total Potential Energy and its Components [E]',
     &           /,' Energy Breakdown over Each of the Atoms [A]',
     &           /,' List of the Large Individual Interactions [L]',
     &           /,' Details for All Individual Interactions [D]',
     &           /,' Electrostatic Moments and Principle Axes [M]',
     &           /,' Internal Virial, dE/dV Values & Pressure [V]',
     &           /,' Connectivity Lists for Each of the Atoms [C]')
   20    continue
         write (iout,30)
   30    format (/,' Enter the Desired Analysis Types',
     &              ' [G,P,E,A,L,D,M,V,C] :  ',$)
         read (input,40,err=20)  string
   40    format (a240)
      end if
c
c     set option control flags based desired analysis types
c
      dosystem = .false.
      doparam = .false.
      doenergy = .false.
      doatom = .false.
      dolarge = .false.
      dodetail = .false.
      domoment = .false.
      dovirial = .false.
      doconect = .false.
      call upcase (string)
      do i = 1, trimtext(string)
         letter = string(i:i)
         if (letter .eq. 'G')  dosystem = .true.
         if (letter .eq. 'P')  doparam = .true.
         if (letter .eq. 'E')  doenergy = .true.
         if (letter .eq. 'A')  doatom = .true.
         if (letter .eq. 'L')  dolarge = .true.
         if (letter .eq. 'D')  dodetail = .true.
         if (letter .eq. 'M')  domoment = .true.
         if (letter .eq. 'V')  dovirial = .true.
         if (letter .eq. 'C')  doconect = .true.
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (active(n))
c
c     get the list of atoms for which output is desired
c
      if (doatom .or. doparam .or. doconect) then
         do i = 1, 20
            list(i) = 0
         end do
         if (exist) then
            do i = 1, 20
               call nextarg (string,exist)
               if (.not. exist)  goto 50
               read (string,*,err=50,end=50)  list(i)
            end do
   50       continue
            if (list(1) .eq. 0) then
               list(1) = -1
               list(2) = n
            end if
         else
            write (iout,60)
   60       format (/,' List Atoms for which Output is Desired',
     &                 ' [ALL] :  '/,'    >  ',$)
            read (input,70)  record
   70       format (a240)
            read (record,*,err=80,end=80)  (list(i),i=1,20)
   80       continue
         end if
         do i = 1, n
            active(i) = .true.
         end do
         i = 1
         do while (list(i) .ne. 0)
            if (i .eq. 1) then
               do j = 1, n
                  active(j) = .false.
               end do
            end if
            if (list(i) .gt. 0) then
               active(list(i)) = .true.
               i = i + 1
            else
               do j = abs(list(i)), abs(list(i+1))
                  active(j) = .true.
               end do
               i = i + 2
            end if
         end do
      end if
c
c     setup to write out the large individual energy terms
c
      if (dolarge) then
         verbose = .true.
      end if
c
c     setup to write out all of the individual energy terms
c
      if (dodetail) then
         doenergy = .true.
         debug = .true.
         verbose = .true.
      else
         debug = .false.
      end if
c
c     reopen the coordinates file and read the first structure
c
      frame = 0
      ixyz = freeunit ()
      xyzfile = filename
      call suffix (xyzfile,'xyz','old')
      open (unit=ixyz,file=xyzfile,status ='old')
      rewind (unit=ixyz)
      call readxyz (ixyz)
c
c     get info on the molecular system and force field
c
      if (dosystem)  call systyze
c
c     get parameters used for molecular mechanics potentials
c
      if (doparam .and. doconect) then
         call amberyze (active)
      else if (doparam) then
         call paramyze (active)
      end if
c
c     provide connectivity lists for the individual atoms
c
      if (doconect)  call connyze (active)
c
c     perform analysis for each successive coordinate structure
c
      do while (.not. abort)
         frame = frame + 1
         if (frame .gt. 1) then
            write (iout,90)  frame
   90       format (/,' Analysis for Archive Structure :',8x,i8)
            if (nold .ne. n)  call mechanic
         end if
c
c     make the call to compute the potential energy
c
         if (doenergy .or. doatom .or. dolarge)  call enrgyze
c
c     energy partitioning by potential energy components
c
         if (doenergy)  call partyze
c
c     get the various electrostatic and inertial moments
c
         if (domoment) then
            debug = .false.
            call momyze
            if (dodetail)  debug = .true.
         end if
c
c     get and test the internal virial and pressure values
c
         if (dovirial) then
            debug = .false.
            call viriyze
            if (dodetail)  debug = .true.
         end if
c
c     energy partitioning over the individual atoms
c
         if (doatom)  call atomyze (active)
c
c     attempt to read next structure from the coordinate file
c
         nold = n
         call readxyz (ixyz)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (active)
c
c     perform any final tasks before program exit
c
      close (unit=ixyz)
      if (dodetail)  debug = .false.
      call final
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine systyze  --  system & force field information  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "systyze" is an auxiliary routine for the analyze program
c     that prints general information about the molecular system
c     and the force field model
c
c
      subroutine systyze
      use sizes
      use atoms
      use bound
      use boxes
      use ewald
      use fields
      use iounit
      use limits
      use molcul
      use pme
      use potent
      use units
      use vdwpot
      implicit none
      integer i
      real*8 dens
      character*20 value
      character*20 label(5)
c
c
c     info on number of atoms, molecules and mass
c
      if (n .ne. 0) then
         write (iout,10)  n,nmol,totmass
   10    format (/,' Overall System Contents :',
     &           //,' Number of Atoms',25x,i8,
     &           /,' Number of Molecules',21x,i8,
     &           /,' Total System Mass',15x,f16.4)
         if (use_bounds) then
            dens = (1.0d24/volbox) * (totmass/avogadro)
            write (iout,20)  dens
   20       format (' System Density',22x,f12.4)
         end if
      end if
c
c     periodic box dimensions and crystal lattice type
c
      if (use_bounds) then
         value = 'ORTHOGONAL'
         if (monoclinic)  value = 'MONOCLINIC'
         if (triclinic)  value = 'TRICLINIC'
         if (octahedron)  value = 'TRUNCATED OCTAHEDRON'
         call justify (value)
         write (iout,30)  xbox,ybox,zbox,alpha,beta,gamma,volbox,value
   30    format (/,' Periodic Boundary Conditions :',
     &           //,' a-Axis Length',23x,f12.4,
     &           /,' b-Axis Length',23x,f12.4,
     &           /,' c-Axis Length',23x,f12.4,
     &           /,' Alpha Angle',25x,f12.4,
     &           /,' Beta Angle',26x,f12.4,
     &           /,' Gamma Angle',25x,f12.4,
     &           /,' Cell Volume',21x,f16.4,
     &           /,' Lattice Type',16x,a20)
         if (spacegrp .ne. '          ') then
            value = spacegrp
            call justify (value)
            write (iout,40)  value
   40       format (' Space Group',17x,a20)
         end if
      end if
c
c     info on force field potential energy function
c
      value = forcefield
      call justify (value)
      write (iout,50)  value
   50 format (/,' Force Field Name :',10x,a20)
c
c     details of vdw potential energy functional form
c
      if (use_vdw) then
         label(1) = vdwtyp
         label(2) = radtyp
         label(3) = radsiz
         label(4) = radrule
         label(5) = epsrule
         do i = 1, 5
            call justify (label(i))
         end do
         write (iout,60)  (label(i),i=1,5)
   60    format (/,' VDW Function',16x,a20,
     &           /,' Size Descriptor',13x,a20,
     &           /,' Size Unit Type',14x,a20,
     &           /,' Size Combining Rule',9x,a20,
     &           /,' Well Depth Rule',13x,a20)
         if (vdwcut .le. 1000.0d0) then
            write (iout,70)  vdwcut
   70       format (' VDW Cutoff',26x,f12.4)
         end if
      end if
c
c     details of electrostatic energy functional form
c
      if (use_charge .or. use_dipole .or. use_mpole .or. use_polar) then
         write (iout,80)
   80    format ()
      end if
      if (use_charge) then
         value = 'PARTIAL CHARGE'
         call justify (value)
         write (iout,90)  value
   90    format (' Electrostatics',14x,a20)
      end if
      if (use_dipole) then
         value = 'BOND DIPOLE'
         call justify (value)
         write (iout,100)  value
  100    format (' Electrostatics',14x,a20)
      end if
      if (use_mpole) then
         value = 'ATOMIC MULTIPOLE'
         call justify (value)
         write (iout,110)  value
  110    format (' Electrostatics',14x,a20)
      end if
      if (use_polar) then
         value = 'INDUCED DIPOLE'
         call justify (value)
         write (iout,120)  value
  120    format (' Electrostatics',14x,a20)
      end if
c
c     details of particle mesh Ewald calculation
c
      if (use_ewald) then
         value = boundary
         call justify (value)
         write (iout,130)  aewald,ewaldcut,nfft1,nfft2,nfft3,
     &                    bsorder,value
  130    format (/,' Particle Mesh Ewald :',
     &           //,' Ewald Coefficient',19x,f12.4,
     &           /,' Real-Space Cutoff',19x,f12.4,
     &           /,' Grid Dimensions',21x,3i4,
     &           /,' B-Spline Order',26x,i8,
     &           /,' Boundary Condition',10x,a20)
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine paramyze  --  force field parameter analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "paramyze" prints the force field parameters used in the
c     computation of each of the potential energy terms
c
c
      subroutine paramyze (active)
      use sizes
      use angang
      use angbnd
      use angpot
      use angtor
      use atomid
      use atoms
      use bitor
      use bndstr
      use charge
      use dipole
      use improp
      use imptor
      use iounit
      use korbs
      use ktrtor
      use kvdws
      use math
      use mpole
      use opbend
      use opdist
      use piorbs
      use pistuf
      use pitors
      use polar
      use polgrp
      use potent
      use solute
      use strbnd
      use strtor
      use tors
      use tortor
      use units
      use urey
      use vdw
      use ctran 
      use vdwpot
      implicit none
      integer i,j,k
      integer ia,ib,ic
      integer id,ie,ig
      integer ixaxe
      integer iyaxe
      integer izaxe
      integer fold(9)
      real*8 bla,blc
      real*8 radj,rad4j
      real*8 ampli(9)
      real*8 phase(9)
      real*8 mpl(13)
      logical header
      logical active(*)
c
c
c     number of each type of interaction and site
c
      if (n .ne. 0) then
         write (iout,10)
   10    format (/,' Interactions and Sites :',/)
         write (iout,20)  n
   20    format (' Atomic Sites',21x,i15)
      end if
      if (use_bond .and. nbond.ne.0) then
         write (iout,30)  nbond
   30    format (' Bond Stretches',19x,i15)
      end if
      if (use_angle .and. nangle.ne.0) then
         write (iout,40)  nangle
   40    format (' Angle Bends',22x,i15)
      end if
      if (use_strbnd .and. nstrbnd.ne.0) then
         write (iout,50)  nstrbnd
   50    format (' Stretch-Bends',20x,i15)
      end if
      if (use_urey .and. nurey.ne.0) then
         write (iout,60)  nurey
   60    format (' Urey-Bradley',21x,i15)
      end if
      if (use_angang .and. nangang.ne.0) then
         write (iout,70)  nangang
   70    format (' Angle-Angles',21x,i15)
      end if
      if (use_opbend .and. nopbend.ne.0) then
         write (iout,80)  nopbend
   80    format (' Out-of-Plane Bends',15x,i15)
      end if
      if (use_opdist .and. nopdist.ne.0) then
         write (iout,90)  nopdist
   90    format (' Out-of-Plane Distances',11x,i15)
      end if
      if (use_improp .and. niprop.ne.0) then
         write (iout,100)  niprop
  100    format (' Improper Dihedrals',15x,i15)
      end if
      if (use_imptor .and. nitors.ne.0) then
         write (iout,110)  nitors
  110    format (' Improper Torsions',16x,i15)
      end if
      if (use_tors .and. ntors.ne.0) then
         write (iout,120)  ntors
  120    format (' Torsional Angles',17x,i15)
      end if
      if (use_pitors .and. npitors.ne.0) then
         write (iout,130)  npitors
  130    format (' Pi-Orbital Torsions',14x,i15)
      end if
      if (use_strtor .and. nstrtor.ne.0) then
         write (iout,140)  nstrtor
  140    format (' Stretch-Torsions',17x,i15)
      end if
      if (use_angtor .and. nangtor.ne.0) then
         write (iout,150)  nangtor
  150    format (' Angle-Torsions',19x,i15)
      end if
      if (use_tortor .and. ntortor.ne.0) then
         write (iout,160)  ntortor
  160    format (' Torsion-Torsions',17x,i15)
      end if
      if (use_vdw .and. nvdw.ne.0) then
         write (iout,170)  nvdw
  170    format (' Van der Waals Sites',14x,i15)
      end if
      if (use_ct .and. nct.ne.0) then
         write (iout,175)  nct
  175    format (' Charge Transfer Sites',12x,i15)
      end if
      if (use_charge .and. nion.ne.0) then
         write (iout,180)  nion
  180    format (' Atomic Partial Charges',11x,i15)
      end if
      if (use_dipole .and. ndipole.ne.0) then
         write (iout,190)  ndipole
  190    format (' Bond Dipole Moments',14x,i15)
      end if
      if (use_mpole .and. npole.ne.0) then
         write (iout,200)  npole
  200    format (' Atomic Multipoles',16x,i15)
      end if
      if (use_polar .and. npolar.ne.0) then
         write (iout,210)  npolar
  210    format (' Polarizable Sites',16x,i15)
      end if
      if (use_orbit .and. norbit.ne.0) then
         write (iout,220)  norbit
  220    format (' Pisystem Atoms',19x,i15)
      end if
      if (use_orbit .and. nbpi.ne.0) then
         write (iout,230)  nbpi
  230    format (' Conjugated Pi-Bonds',14x,i15)
      end if
c
c     parameters used for molecular mechanics atom types
c
      header = .true.
      do i = 1, n
         if (active(i)) then
            if (header) then
               header = .false.
               write (iout,240)
  240          format (/,' Atom Type Definition Parameters :',
     &                 //,3x,'Atom',2x,'Symbol',2x,'Type',
     &                    2x,'Class',2x,'Atomic',3x,'Mass',
     &                    2x,'Valence',2x,'Description',/)
            end if
            write (iout,250)  i,name(i),type(i),class(i),atomic(i),
     &                        mass(i),valence(i),story(i)
  250       format (i6,5x,a3,2i7,i6,f10.3,i5,5x,a24)
         end if
      end do
c
c     parameters used for van der Waals interactions
c
      if (use_vdw) then
         header = .true.
         k = 0
         do i = 1, n
            if (active(i)) then
               k = k + 1
               if (header) then
                  header = .false.
                  write (iout,260)
  260             format (/,' Van der Waals Parameters :',
     &                    //,10x,'Atom Number',7x,'Size',
     &                       3x,'Epsilon',3x,'Size 1-4',
     &                       3x,'Eps 1-4',3x,'Reduction',/)
               end if
               j = class(i)
               if (vdwindex .eq. 'TYPE')  j = type(i)
               if (rad(j).eq.rad4(j) .and. eps(j).eq.eps4(j)) then
                  radj = rad(j)
                  if (radsiz .eq. 'DIAMETER')  radj = 2.0d0 * radj
                  if (radtyp .eq. 'SIGMA')  radj = radj / twosix
                  if (reduct(j) .eq. 0.0d0) then
                     write (iout,270)  k,i,radj,eps(j)
  270                format (i6,3x,i6,7x,2f10.4)
                  else
                     write (iout,280)  k,i,radj,eps(j),reduct(j)
  280                format (i6,3x,i6,7x,2f10.4,22x,f10.4)
                  end if
               else
                  radj = rad(j)
                  rad4j = rad4(j)
                  if (radsiz .eq. 'DIAMETER') then
                     radj = 2.0d0 * radj
                     rad4j = 2.0d0 * rad4j
                  end if
                  if (radtyp .eq. 'SIGMA') then
                     radj = radj / twosix
                     rad4j = rad4j / twosix
                  end if
                  if (reduct(j) .eq. 0.0d0) then
                     write (iout,290)  k,i,radj,eps(j),rad4j,eps4(j)
  290                format (i6,3x,i6,7x,2f10.4,1x,2f10.4)
                  else
                     write (iout,300)  k,i,radj,eps(j),rad4j,
     &                                eps4(j),reduct(j)
  300                format (i6,3x,i6,7x,2f10.4,1x,2f10.4,1x,f10.4)
                  end if
               end if
            end if
         end do
      end if
c
c     parameters used for bond stretching interactions
c
      if (use_bond) then
         header = .true.
         do i = 1, nbond
            ia = ibnd(1,i)
            ib = ibnd(2,i)
            if (active(ia) .or. active(ib)) then
               if (header) then
                  header = .false.
                  write (iout,310)
  310             format (/,' Bond Stretching Parameters :',
     &                    //,10x,'Atom Numbers',25x,'KS',7x,'Bond',/)
               end if
               write (iout,320)  i,ia,ib,bk(i),bl(i)
  320          format (i6,3x,2i6,19x,f10.3,f10.4)
            end if
         end do
      end if
c
c     parameters used for angle bending interactions
c
      if (use_angle) then
         header = .true.
         do i = 1, nangle
            ia = iang(1,i)
            ib = iang(2,i)
            ic = iang(3,i)
            if (active(ia) .or. active(ib) .or. active(ic)) then
               if (header) then
                  header = .false.
                  write (iout,330)
  330             format (/,' Angle Bending Parameters :',
     &                    //,13x,'Atom Numbers',22x,'KB',
     &                       6x,'Angle',3x,'Fold',4x,'Type',/)
               end if
               if (angtyp(i) .eq. 'HARMONIC') then
                  write (iout,340)  i,ia,ib,ic,ak(i),anat(i)
  340             format (i6,3x,3i6,13x,2f10.3)
               else if (angtyp(i) .eq. 'IN-PLANE') then
                  write (iout,350)  i,ia,ib,ic,ak(i),anat(i)
  350             format (i6,3x,3i6,13x,2f10.3,9x,'In-Plane')
               else if (angtyp(i) .eq. 'IN-PLANE') then
                  write (iout,360)  i,ia,ib,ic,ak(i),anat(i)
  360             format (i6,3x,3i6,13x,2f10.3,9x,'Linear')
               else if (angtyp(i) .eq. 'FOURIER ') then
                  write (iout,370)  i,ia,ib,ic,ak(i),anat(i),afld(i)
  370             format (i6,3x,3i6,13x,2f10.3,f7.1,2x,'Fourier')
               end if
            end if
         end do
      end if
c
c     parameters used for stretch-bend interactions
c
      if (use_strbnd) then
         header = .true.
         do i = 1, nstrbnd
            k = isb(1,i)
            ia = iang(1,k)
            ib = iang(2,k)
            ic = iang(3,k)
            if (active(ia) .or. active(ib) .or. active(ic)) then
               if (header) then
                  header = .false.
                  write (iout,380)
  380             format (/,' Stretch-Bend Parameters :',
     &                    //,13x,'Atom Numbers',8x,'KSB 1',5x,'KSB 2',
     &                       6x,'Angle',3x,'Bond 1',3x,'Bond 2',/)
               end if
               bla = 0.0d0
               blc = 0.0d0
               if (isb(2,i) .ne. 0)  bla = bl(isb(2,i))
               if (isb(3,i) .ne. 0)  blc = bl(isb(3,i))
               write (iout,390)  i,ia,ib,ic,sbk(1,i),sbk(2,i),
     &                           anat(k),bla,blc
  390          format (i6,3x,3i6,1x,2f10.3,2x,f9.3,2f9.4)
            end if
         end do
      end if
c
c     parameters used for Urey-Bradley interactions
c
      if (use_urey) then
         header = .true.
         do i = 1, nurey
            ia = iury(1,i)
            ib = iury(2,i)
            ic = iury(3,i)
            if (active(ia) .or. active(ic)) then
               if (header) then
                  header = .false.
                  write (iout,400)
  400             format (/,' Urey-Bradley Parameters :',
     &                    //,13x,'Atom Numbers',21x,'KUB',
     &                       4x,'Distance',/)
               end if
               write (iout,410)  i,ia,ib,ic,uk(i),ul(i)
  410          format (i6,3x,3i6,13x,f10.3,f10.4)
            end if
         end do
      end if
c
c     parameters used for out-of-plane bend interactions
c
      if (use_opbend) then
         header = .true.
         do i = 1, nopbend
            k = iopb(i)
            ia = iang(1,k)
            ib = iang(2,k)
            ic = iang(3,k)
            id = iang(4,k)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,420)
  420             format (/,' Out-of-Plane Bend Parameters :',
     &                    //,17x,'Atom Numbers',19x,'KOPB',/)
               end if
               write (iout,430)  i,id,ib,ia,ic,opbk(i)
  430          format (i6,3x,4i6,9x,f10.3)
            end if
         end do
      end if
c
c     parameters used for out-of-plane distance interactions
c
      if (use_opdist) then
         header = .true.
         do i = 1, nopdist
            ia = iopd(1,i)
            ib = iopd(2,i)
            ic = iopd(3,i)
            id = iopd(4,i)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,440)
  440             format (/,' Out-of-Plane Distance Parameters :',
     &                    //,17x,'Atom Numbers',19x,'KOPD',/)
               end if
               write (iout,450)  i,ia,ib,ic,id,opdk(i)
  450          format (i6,3x,4i6,9x,f10.3)
            end if
         end do
      end if
c
c     parameters used for improper dihedral interactions
c
      if (use_improp) then
         header = .true.
         do i = 1, niprop
            ia = iiprop(1,i)
            ib = iiprop(2,i)
            ic = iiprop(3,i)
            id = iiprop(4,i)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,460)
  460             format (/,' Improper Dihedral Parameters :',
     &                    //,17x,'Atom Numbers',19x,'KID',
     &                       4x,'Dihedral',/)
               end if
               write (iout,470)  i,ia,ib,ic,id,kprop(i),vprop(i)
  470          format (i6,3x,4i6,9x,2f10.4)
            end if
         end do
      end if
c
c     parameters used for improper torsion interactions
c
      if (use_imptor) then
         header = .true.
         do i = 1, nitors
            ia = iitors(1,i)
            ib = iitors(2,i)
            ic = iitors(3,i)
            id = iitors(4,i)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,480)
  480             format (/,' Improper Torsion Parameters :',
     &                    //,17x,'Atom Numbers',11x,
     &                       'Amplitude, Phase and Periodicity',/)
               end if
               j = 0
               if (itors1(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 1
                  ampli(j) = itors1(1,i)
                  phase(j) = itors1(2,i)
               end if
               if (itors2(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 2
                  ampli(j) = itors2(1,i)
                  phase(j) = itors2(2,i)
               end if
               if (itors3(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 3
                  ampli(j) = itors3(1,i)
                  phase(j) = itors3(2,i)
               end if
               if (j .eq. 0) then
                  write (iout,490)  i,ia,ib,ic,id
  490             format (i6,3x,4i6)
               else if (j .eq. 1) then
                  write (iout,500)  i,ia,ib,ic,id,
     &                              ampli(1),phase(1),fold(1)
  500             format (i6,3x,4i6,10x,f10.3,f8.1,i4)
               else if (j .eq. 2) then
                  write (iout,510)  i,ia,ib,ic,id,(ampli(k),
     &                              phase(k),fold(k),k=1,j)
  510             format (i6,3x,4i6,2x,2(f10.3,f6.1,i4))
               else
                  write (iout,520)  i,ia,ib,ic,id,(ampli(k),
     &                              nint(phase(k)),fold(k),k=1,j)
  520             format (i6,3x,4i6,4x,3(f8.3,i4,'/',i1))
               end if
            end if
         end do
      end if
c
c     parameters used for torsional interactions
c
      if (use_tors) then
         header = .true.
         do i = 1, ntors
            ia = itors(1,i)
            ib = itors(2,i)
            ic = itors(3,i)
            id = itors(4,i)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,530)
  530             format (/,' Torsional Angle Parameters :',
     &                    //,17x,'Atom Numbers',11x,
     &                       'Amplitude, Phase and Periodicity',/)
               end if
               j = 0
               if (tors1(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 1
                  ampli(j) = tors1(1,i)
                  phase(j) = tors1(2,i)
               end if
               if (tors2(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 2
                  ampli(j) = tors2(1,i)
                  phase(j) = tors2(2,i)
               end if
               if (tors3(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 3
                  ampli(j) = tors3(1,i)
                  phase(j) = tors3(2,i)
               end if
               if (tors4(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 4
                  ampli(j) = tors4(1,i)
                  phase(j) = tors4(2,i)
               end if
               if (tors5(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 5
                  ampli(j) = tors5(1,i)
                  phase(j) = tors5(2,i)
               end if
               if (tors6(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 6
                  ampli(j) = tors6(1,i)
                  phase(j) = tors6(2,i)
               end if
               if (j .eq. 0) then
                  write (iout,540)  i,ia,ib,ic,id
  540             format (i6,3x,4i6)
               else
                  write (iout,550)  i,ia,ib,ic,id,(ampli(k),
     &                              nint(phase(k)),fold(k),k=1,j)
  550             format (i6,3x,4i6,4x,6(f8.3,i4,'/',i1))
               end if
            end if
         end do
      end if
c
c     parameters used for pi-system torsion interactions
c
      if (use_pitors) then
         header = .true.
         do i = 1, npitors
            ia = ipit(1,i)
            ib = ipit(2,i)
            ic = ipit(3,i)
            id = ipit(4,i)
            ie = ipit(5,i)
            ig = ipit(6,i)
            if (active(ia) .or. active(ib) .or. active(ic) .or.
     &             active(id) .or. active(ie) .or. active(ig)) then
               if (header) then
                  header = .false.
                  write (iout,560)
  560             format (/,' Pi-Orbital Torsion Parameters :',
     &                    //,10x,'Atom Numbers',19x,'Amplitude',/)
               end if
               write (iout,570)  i,ic,id,kpit(i)
  570          format (i6,3x,2i6,19x,f10.4)
            end if
         end do
      end if
c
c     parameters used for stretch-torsion interactions
c
      if (use_strtor) then
         header = .true.
         do i = 1, nstrtor
            k = ist(1,i)
            ia = itors(1,k)
            ib = itors(2,k)
            ic = itors(3,k)
            id = itors(4,k)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,580)
  580             format (/,' Stretch-Torsion Parameters :',
     &                    //,17x,'Atom Numbers',10x,'Bond',
     &                       5x,'Amplitude and Phase (1-3 Fold)',/)
               end if
               ampli(1) = kst(1,i)
               phase(1) = tors1(2,k)
               ampli(2) = kst(2,i)
               phase(2) = tors2(2,k)
               ampli(3) = kst(3,i)
               phase(3) = tors3(2,k)
               ampli(4) = kst(4,i)
               phase(4) = tors1(2,k)
               ampli(5) = kst(5,i)
               phase(5) = tors2(2,k)
               ampli(6) = kst(6,i)
               phase(6) = tors3(2,k)
               ampli(7) = kst(7,i)
               phase(7) = tors1(2,k)
               ampli(8) = kst(8,i)
               phase(8) = tors2(2,k)
               ampli(9) = kst(9,i)
               phase(9) = tors3(2,k)
               write (iout,590)  i,ia,ib,ic,id,
     &                           '1st',(ampli(k),nint(phase(k)),k=1,3),
     &                           '2nd',(ampli(k),nint(phase(k)),k=4,6),
     &                           '3rd',(ampli(k),nint(phase(k)),k=7,9)
  590          format (i6,3x,4i6,7x,a3,3x,3(f7.3,i4),
     &                 /,40x,a3,3x,3(f7.3,i4),
     &                 /,40x,a3,3x,3(f7.3,i4))
            end if
         end do
      end if
c
c     parameters used for angle-torsion interactions
c
      if (use_angtor) then
         header = .true.
         do i = 1, nangtor
            k = iat(1,i)
            ia = itors(1,k)
            ib = itors(2,k)
            ic = itors(3,k)
            id = itors(4,k)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,600)
  600             format (/,' Angle-Torsion Parameters :',
     &                    //,17x,'Atom Numbers',10x,'Angle',
     &                       4x,'Amplitude and Phase (1-3 Fold)',/)
               end if
               ampli(1) = kant(1,i)
               phase(1) = tors1(2,k)
               ampli(2) = kant(2,i)
               phase(2) = tors2(2,k)
               ampli(3) = kant(3,i)
               phase(3) = tors3(2,k)
               ampli(4) = kant(4,i)
               phase(4) = tors1(2,k)
               ampli(5) = kant(5,i)
               phase(5) = tors2(2,k)
               ampli(6) = kant(6,i)
               phase(6) = tors3(2,k)
               write (iout,610)  i,ia,ib,ic,id,
     &                           '1st',(ampli(k),nint(phase(k)),k=1,3),
     &                           '2nd',(ampli(k),nint(phase(k)),k=4,6)
  610          format (i6,3x,4i6,7x,a3,3x,3(f7.3,i4),
     &                 /,40x,a3,3x,3(f7.3,i4))
            end if
         end do
      end if
c
c     parameters used for torsion-torsion interactions
c
      if (use_tortor) then
         header = .true.
         do i = 1, ntortor
            k = itt(1,i)
            ia = ibitor(1,k)
            ib = ibitor(2,k)
            ic = ibitor(3,k)
            id = ibitor(4,k)
            ie = ibitor(5,k)
            if (active(ia) .or. active(ib) .or. active(ic) .or.
     &                active(id) .or. active(ie)) then
               if (header) then
                  header = .false.
                  write (iout,620)
  620             format (/,' Torsion-Torsion Parameters :',
     &                    //,20x,'Atom Numbers',18x,'Spline Grid',/)
               end if
               j = itt(2,i)
               write (iout,630)  i,ia,ib,ic,id,ie,tnx(j),tny(j)
  630          format (i6,3x,5i6,10x,2i6)
            end if
         end do
      end if
c
c     parameters used for atomic partial charges
c
      if (use_charge .or. use_chgdpl) then
         header = .true.
         do i = 1, nion
            ia = iion(i)
            ib = jion(i)
            ic = kion(i)
            if (active(ia) .or. active(ic)) then
               if (header) then
                  header = .false.
                  write (iout,640)
  640             format (/,' Atomic Partial Charge Parameters :',
     &                    /,45x,'Neighbor',3x,'Cutoff',
     &                    /,10x,'Atom Number',13x,'Charge',
     &                       7x,'Site',6x,'Site',/)
               end if
               if (ia.eq.ib .and. ia.eq.ic) then
                  write (iout,650)  i,ia,pchg(i)
  650             format (i6,3x,i6,15x,f10.4)
               else
                  write (iout,660)  i,ia,pchg(i),ib,ic
  660             format (i6,3x,i6,15x,f10.4,5x,i6,4x,i6)
               end if
            end if
         end do
      end if
c
c     parameters used for bond dipole moments
c
      if (use_dipole .or. use_chgdpl) then
         header = .true.
         do i = 1, ndipole
            ia = idpl(1,i)
            ib = idpl(2,i)
            if (active(ia) .or. active(ib)) then
               if (header) then
                  header = .false.
                  write (iout,670)
  670             format (/,' Bond Dipole Moment Parameters :',
     &                    //,10x,'Atom Numbers',22x,'Dipole',
     &                       3x,'Position',/)
               end if
               write (iout,680)  i,ia,ib,bdpl(i),sdpl(i)
  680          format (i6,3x,2i6,19x,f10.4,f10.3)
            end if
         end do
      end if
c
c     parameters used for atomic multipole moments
c
      if (use_mpole) then
         header = .true.
         do i = 1, npole
            ia = ipole(i)
            if (active(ia)) then
               if (header) then
                  header = .false.
                  write (iout,690)
  690             format (/,' Atomic Multipole Parameters :',
     &                    //,11x,'Atom',3x,'Z-Axis',1x,'X-Axis',
     &                       1x,'Y-Axis',2x,
     &                       'Frame',11x,'Multipole Moments',/)
               end if
               izaxe = zaxis(i)
               ixaxe = xaxis(i)
               iyaxe = yaxis(i)
               if (iyaxe .lt. 0)  iyaxe = -iyaxe
               mpl(1) = pole(1,i)
               do j = 2, 4
                  mpl(j) = pole(j,i) / bohr
               end do
               do j = 5, 13
                  mpl(j) = 3.0d0 * pole(j,i) / bohr**2
               end do
               if (izaxe .eq. 0) then
                  write (iout,700)  i,ia,polaxe(i),
     &                              (mpl(j),j=1,5),mpl(8),mpl(9),
     &                              (mpl(j),j=11,13)
  700             format (i6,3x,i6,25x,a8,2x,f9.5,/,50x,3f9.5,
     &                    /,50x,f9.5,/,50x,2f9.5,/,50x,3f9.5)
               else if (ixaxe .eq. 0) then
                  write (iout,710)  i,ia,izaxe,polaxe(i),
     &                              (mpl(j),j=1,5),mpl(8),mpl(9),
     &                              (mpl(j),j=11,13)
  710             format (i6,3x,i6,1x,i7,17x,a8,2x,f9.5,/,50x,3f9.5,
     &                    /,50x,f9.5,/,50x,2f9.5,/,50x,3f9.5)
               else  if (iyaxe .eq. 0) then
                  write (iout,720)  i,ia,izaxe,ixaxe,polaxe(i),
     &                              (mpl(j),j=1,5),mpl(8),mpl(9),
     &                              (mpl(j),j=11,13)
  720             format (i6,3x,i6,1x,2i7,10x,a8,2x,f9.5,/,50x,3f9.5,
     &                    /,50x,f9.5,/,50x,2f9.5,/,50x,3f9.5)
               else
                  write (iout,730)  i,ia,izaxe,ixaxe,iyaxe,polaxe(i),
     &                              (mpl(j),j=1,5),mpl(8),mpl(9),
     &                              (mpl(j),j=11,13)
  730             format (i6,3x,i6,1x,3i7,3x,a8,2x,f9.5,/,50x,3f9.5,
     &                    /,50x,f9.5,/,50x,2f9.5,/,50x,3f9.5)
               end if
            end if
         end do
      end if
c
c     parameters used for dipole polarizability
c
      if (use_polar) then
         header = .true.
         do i = 1, npole
            ia = ipole(i)
            if (active(ia)) then
               if (header) then
                  header = .false.
                  write (iout,740)
  740             format (/,' Dipole Polarizability Parameters :',
     &                    //,10x,'Atom Number',5x,'Alpha',5x,'Damp',
     &                       6x,'Polarization Group',/)
               end if
               write (iout,750)  i,ia,polarity(i),thole(i),dirdamp(i),
     &                           (ip11(j,ia),j=1,np11(ia))
  750          format (i6,3x,i6,6x,f10.4,f9.3,3x,f9.3,3x,120i6)
            end if
         end do
      end if
c
c     parameters used for empirical solvation
c
      if (use_solv) then
         header = .true.
         k = 0
         do i = 1, n
            if (active(i)) then
               k = k + 1
               if (header) then
                  header = .false.
                  write (iout,760)
  760             format (/,' Empirical Solvation Parameters :',
     &                    //,10x,'Atom Number',13x,'Radius',
     &                       3x,'ASP Value',/)
               end if
               write (iout,770)  k,i,rsolv(i),asolv(i)
  770          format (i6,3x,i6,15x,2f10.4)
            end if
         end do
      end if
c
c     parameters used for conjugated pisystem atoms
c
      if (use_orbit) then
         header = .true.
         do i = 1, norbit
            ia = iorbit(i)
            j = class(ia)
            if (header) then
               header = .false.
               write (iout,780)
  780          format (/,' Conjugated Pi-Atom Parameters :',
     &                 //,10x,'Atom Number',14x,'Nelect',
     &                    6x,'Ionize',4x,'Repulsion',/)
            end if
            write (iout,790)  i,ia,electron(j),ionize(j),repulse(j)
  790       format (i6,3x,i6,17x,f8.1,3x,f10.4,2x,f10.4)
         end do
      end if
c
c     parameters used for conjugated pibond interactions
c
      if (use_orbit) then
         header = .true.
         do i = 1, nbpi
            ia = ibpi(2,i)
            ib = ibpi(3,i)
            if (header) then
               header = .false.
               write (iout,800)
  800          format (/,' Conjugated Pi-Bond Parameters :',
     &                 //,10x,'Atom Numbers',21x,'K Slope',
     &                    3x,'L Slope',/)
            end if
            write (iout,810)  i,ia,ib,kslope(i),lslope(i)
  810       format (i6,3x,2i6,19x,2f10.4)
         end do
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine enrgyze  --  compute & report energy analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "enrgyze" is an auxiliary routine for the analyze program
c     that performs the energy analysis and prints the total and
c     intermolecular energies
c
c
      subroutine enrgyze
      use sizes
      use atoms
      use inform
      use inter
      use iounit
      use limits
      use molcul
      implicit none
      real*8 energy
      character*240 fstr
c
c
c     perform the energy analysis by atom and component
c
      call analysis (energy)
c
c     intermolecular energy for systems with multiple molecules
c
      fstr = '(/,'' Intermolecular Energy :'',9x,f16.4,'' Kcal/mole'')'
      if (digits .ge. 6)  fstr(31:38) = '7x,f18.6'
      if (digits .ge. 8)  fstr(31:38) = '5x,f20.8'
      if (abs(einter) .ge. 1.0d10)  fstr(34:34) = 'd'
      if (nmol.gt.1 .and. nmol.lt.n .and. .not.use_ewald)
     &   write (iout,fstr)  einter
c
c     print out the total potential energy of the system
c
      fstr = '(/,'' Total Potential Energy :'',8x,f16.4,'' Kcal/mole'')'
      if (digits .ge. 6)  fstr(32:39) = '6x,f18.6'
      if (digits .ge. 8)  fstr(32:39) = '4x,f20.8'
      if (abs(energy) .ge. 1.0d10)  fstr(35:35) = 'd'
      write (iout,fstr)  energy
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine partyze  --  energy component decomposition  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "partyze" prints the energy component and number of
c     interactions for each of the potential energy terms
c
c
      subroutine partyze
      use action
      use energi
      use inform
      use iounit
      use limits
      use potent
      implicit none
      character*12 form1
      character*12 form2
      character*240 fstr
c
c
c     write out each energy component to the desired precision
c
      form1 = '5x,f16.4,i17'
      if (digits .ge. 6)  form1 = '3x,f18.6,i17'
      if (digits .ge. 8)  form1 = '1x,f20.8,i17'
      form2 = form1(1:3)//'d'//form1(5:12)
      fstr = '(/,'' Energy Component Breakdown :'',
     &          11x,''Kcal/mole'',8x,''Interactions''/)'
      write (iout,fstr)
      if (use_bond .and. (neb.ne.0.or.eb.ne.0.0d0)) then
         fstr = '('' Bond Stretching'',12x,'//form1//')'
         write (iout,fstr)  eb,neb
      end if
      if (use_angle .and. (nea.ne.0.or.ea.ne.0.0d0)) then
         fstr = '('' Angle Bending'',14x,'//form1//')'
         write (iout,fstr)  ea,nea
      end if
      if (use_strbnd .and. (neba.ne.0.or.eba.ne.0.0d0)) then
         fstr = '('' Stretch-Bend'',15x,'//form1//')'
         write (iout,fstr)  eba,neba
      end if
      if (use_urey .and. (neub.ne.0.or.eub.ne.0.0d0)) then
         fstr = '('' Urey-Bradley'',15x,'//form1//')'
         write (iout,fstr)  eub,neub
      end if
      if (use_angang .and. (neaa.ne.0.or.eaa.ne.0.0d0)) then
         fstr = '('' Angle-Angle'',16x,'//form1//')'
         write (iout,fstr)  eaa,neaa
      end if
      if (use_opbend .and. (neopb.ne.0.or.eopb.ne.0.0d0)) then
         fstr = '('' Out-of-Plane Bend'',10x,'//form1//')'
         write (iout,fstr)  eopb,neopb
      end if
      if (use_opdist .and. (neopd.ne.0.or.eopd.ne.0.0d0)) then
         fstr = '('' Out-of-Plane Distance'',6x,'//form1//')'
         write (iout,fstr)  eopd,neopd
      end if
      if (use_improp .and. (neid.ne.0.or.eid.ne.0.0d0)) then
         fstr = '('' Improper Dihedral'',10x,'//form1//')'
         write (iout,fstr)  eid,neid
      end if
      if (use_imptor .and. (neit.ne.0.or.eit.ne.0.0d0)) then
         fstr = '('' Improper Torsion'',11x,'//form1//')'
         write (iout,fstr)  eit,neit
      end if
      if (use_tors .and. (net.ne.0.or.et.ne.0.0d0)) then
         fstr = '('' Torsional Angle'',12x,'//form1//')'
         write (iout,fstr)  et,net
      end if
      if (use_pitors .and. (nept.ne.0.or.ept.ne.0.0d0)) then
         fstr = '('' Pi-Orbital Torsion'',9x,'//form1//')'
         write (iout,fstr)  ept,nept
      end if
      if (use_strtor .and. (nebt.ne.0.or.ebt.ne.0.0d0)) then
         fstr = '('' Stretch-Torsion'',12x,'//form1//')'
         write (iout,fstr)  ebt,nebt
      end if
      if (use_angtor .and. (neat.ne.0.or.eat.ne.0.0d0)) then
         fstr = '('' Angle-Torsion'',14x,'//form1//')'
         write (iout,fstr)  eat,neat
      end if
      if (use_tortor .and. (nett.ne.0.or.ett.ne.0.0d0)) then
         fstr = '('' Torsion-Torsion'',12x,'//form1//')'
         write (iout,fstr)  ett,nett
      end if
      if (use_vdw .and. (nev.ne.0.or.ev.ne.0.0d0)) then
         if (abs(ev) .lt. 1.0d10) then
            fstr = '('' Van der Waals'',14x,'//form1//')'
         else
            fstr = '('' Van der Waals'',14x,'//form2//')'
         end if
         write (iout,fstr)  ev,nev
      end if

CW!!! write out the charge transfer energy      
      if (use_ct .and. (nect.ne.0.or.ect.ne.0.0d0)) then
        if (abs(ect) .lt. 1.0d10) then
          fstr = '('' Charge Transfer'',12x,'//form1//')'
        else
          fstr = '('' Charge Transfer'',12x,'//form2//')'
        end if
        write (iout,fstr)  ect,nect
      end if

CW!!!END
      if (use_charge .and. (nec.ne.0.or.ec.ne.0.0d0)) then
         if (abs(ec) .lt. 1.0d10) then
            fstr = '('' Charge-Charge'',14x,'//form1//')'
         else
            fstr = '('' Charge-Charge'',14x,'//form2//')'
         end if
         write (iout,fstr)  ec,nec
      end if
      if (use_chgdpl .and. (necd.ne.0.or.ecd.ne.0.0d0)) then
         if (abs(ecd) .lt. 1.0d10) then
            fstr = '('' Charge-Dipole'',14x,'//form1//')'
         else
            fstr = '('' Charge-Dipole'',14x,'//form2//')'
         end if
         write (iout,fstr)  ecd,necd
      end if
      if (use_dipole .and. (ned.ne.0.or.ed.ne.0.0d0)) then
         if (abs(ed) .lt. 1.0d10) then
            fstr = '('' Dipole-Dipole'',14x,'//form1//')'
         else
            fstr = '('' Dipole-Dipole'',14x,'//form2//')'
         end if
         write (iout,fstr)  ed,ned
      end if
      if (use_mpole .and. (nem.ne.0.or.em.ne.0.0d0)) then
         if (abs(em) .lt. 1.0d10) then
            fstr = '('' Atomic Multipoles'',10x,'//form1//')'
         else
            fstr = '('' Atomic Multipoles'',10x,'//form2//')'
         end if
         write (iout,fstr)  em,nem
      end if
      if (use_polar .and. (nep.ne.0.or.ep.ne.0.0d0)) then
         if (abs(ep) .lt. 1.0d10) then
            fstr = '('' Polarization'',15x,'//form1//')'
         else
            fstr = '('' Polarization'',15x,'//form2//')'
         end if
         write (iout,fstr)  ep,nep
      end if
      if (use_rxnfld .and. (ner.ne.0.or.er.ne.0.0d0)) then
         fstr = '('' Reaction Field'',13x,'//form1//')'
         write (iout,fstr)  er,ner
      end if
      if (use_solv .and. (nes.ne.0.or.es.ne.0.0d0)) then
         fstr = '('' Implicit Solvation'',9x,'//form1//')'
         write (iout,fstr)  es,nes
      end if
      if (use_metal .and. (nelf.ne.0.or.elf.ne.0.0d0)) then
         fstr = '('' Metal Ligand Field'',9x,'//form1//')'
         write (iout,fstr)  elf,nelf
      end if
      if (use_geom .and. (neg.ne.0.or.eg.ne.0.0d0)) then
         fstr = '('' Geometric Restraints'',7x,'//form1//')'
         write (iout,fstr)  eg,neg
      end if
      if (use_extra .and. (nex.ne.0.or.ex.ne.0.0d0)) then
         fstr = '('' Extra Energy Terms'',9x,'//form1//')'
         write (iout,fstr)  ex,nex
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine momyze  --  electrostatic & inertial analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "momyze" finds and prints the total charge, dipole moment
c     components, radius of gyration and moments of inertia
c
c
      subroutine momyze
      use sizes
      use chgpot
      use iounit
      use moment
      implicit none
      real*8 rg
c
c
c     get the total charge, dipole and quadrupole moments
c
      call moments
      write (iout,10)  netchg
   10 format (/,' Total Electric Charge :',12x,f13.5,' Electrons')
      write (iout,20)  netdpl,xdpl,ydpl,zdpl
   20 format (/,' Dipole Moment Magnitude :',10x,f13.3,' Debyes',
     &        //,' Dipole X,Y,Z-Components :',10x,3f13.3)
      write (iout,30)  xxqdp,xyqdp,xzqdp,yxqdp,yyqdp,
     &                 yzqdp,zxqdp,zyqdp,zzqdp
   30 format (/,' Quadrupole Moment Tensor :',9x,3f13.3,
     &        /,6x,'(Buckinghams)',17x,3f13.3,
     &        /,36x,3f13.3)
      write (iout,40)  netqdp(1),netqdp(2),netqdp(3)
   40 format (/,' Principal Axes Quadrupole :',8x,3f13.3)
      if (dielec .ne. 1.0d0) then
         write (iout,50)  dielec
   50    format (/,' Dielectric Constant :',14x,f13.3)
         write (iout,60)  netchg/sqrt(dielec)
   60    format (' Effective Total Charge :',11x,f13.5,' Electrons')
         write (iout,70)  netdpl/sqrt(dielec)
   70    format (' Effective Dipole Moment :',10x,f13.3,' Debyes')
      end if
c
c     get the radius of gyration and moments of inertia
c
      call gyrate (rg)
      write (iout,80)  rg
   80 format (/,' Radius of Gyration :',15x,f13.3,' Angstroms')
      call inertia (1)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine viriyze  --  inertial virial & pressure values  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "propyze" finds and prints the internal virial, analytical and
c     numerical dE/dV values, and the pressure
c
c
      subroutine viriyze
      use sizes
      use atoms
      use iounit
      use virial
      implicit none
      integer i
      real*8 energy
      real*8, allocatable :: derivs(:,:)
c
c
c     perform a gradient calculation in order to find virial
c
      allocate (derivs(3,n))
      call gradient (energy,derivs)
      deallocate (derivs)
c
c     print out the components of the internal virial
c
      write (iout,10)  (vir(1,i),vir(2,i),vir(3,i),i=1,3)
   10 format (/,' Internal Virial Tensor :',12x,3f12.3,
     &        /,37x,3f12.3,/,37x,3f12.3)
c
c     get two alternative dE/dV values and a pressure estimate
c
      call ptest
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine atomyze  --  individual atom energy analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "atomyze" prints the potential energy components broken
c     down by atom and to a choice of precision
c
c
      subroutine atomyze (active)
      use sizes
      use analyz
      use atoms
      use inform
      use iounit
      implicit none
      integer i
      logical active(*)
      character*240 fstr
c
c
c     energy partitioning over the individual atoms
c
      fstr = '(/,'' Potential Energy Breakdown over Atoms :'')'
      write (iout,fstr)
      if (digits .ge. 8) then
         write (iout,10)
   10    format (/,'  Atom',9x,'EB',14x,'EA',14x,'EBA',13x,'EUB',
     &           /,15x,'EAA',13x,'EOPB',12x,'EOPD',12x,'EID',
     &           /,15x,'EIT',13x,'ET',14x,'EPT',13x,'EBT',
     &           /,15x,'EAT',13x,'ETT',13x,'EV',14x,'EC',
     &           /,15x,'ECD',13x,'ED',14x,'EM',14x,'EP',
     &           /,15x,'ER',14x,'ES',14x,'ELF',13x,'EG',
     &           /,15x,'EX')
      else if (digits .ge. 6) then
         write (iout,20)
   20    format (/,'  Atom',8x,'EB',12x,'EA',12x,'EBA',11x,'EUB',
     &              11x,'EAA',
     &           /,14x,'EOPB',10x,'EOPD',10x,'EID',11x,'EIT',
     &              11x,'ET',
     &           /,14x,'EPT',11x,'EBT',11x,'EAT',11x,'ETT',11x,'EV',
     &           /,14x,'EC',12x,'ECD',11x,'ED',12x,'EM',12x,'EP',
     &           /,14x,'ER',12x,'ES',12x,'ELF',11x,'EG',12x,'EX')
      else
         write (iout,30)
   30    format (/,'  Atom',8x,'EB',10x,'EA',10x,'EBA',9x,'EUB',
     &              9x,'EAA',9x,'EOPB',
     &           /,14x,'EOPD',8x,'EID',9x,'EIT',9x,'ET',10x,'EPT',
     &              9x,'EBT',
     &           /,14x,'EAT',9x,'ETT',9x,'EV',10x,'EC',10x,'ECD',
     &              9x,'ED',
     &           /,14x,'EM',10x,'EP',10x,'ER',10x,'ES',10x,'ELF',
     &              9x,'EG',
     &           /,14x,'EX')
      end if
      if (digits .ge. 8) then
         fstr = '(/,i6,4f16.8,/,6x,4f16.8,/,6x,4f16.8,'//
     &             '/,6x,4f16.8,/,6x,4f16.8,/,6x,4f16.8,'//
     &             '/,6x,f16.8)'
      else if (digits .ge. 6) then
         fstr = '(/,i6,5f14.6,/,6x,5f14.6,/,6x,5f14.6,'//
     &             '/,6x,5f14.6,/,6x,5f14.6)'
      else
         fstr = '(/,i6,6f12.4,/,6x,6f12.4,/,6x,6f12.4,'//
     &             '/,6x,6f12.4,/,6x,f12.4)'
      end if
      do i = 1, n
         if (active(i)) then
            write (iout,fstr)  i,aeb(i),aea(i),aeba(i),aeub(i),aeaa(i),
     &                         aeopb(i),aeopd(i),aeid(i),aeit(i),aet(i),
     &                         aept(i),aebt(i),aeat(i),aett(i),aev(i),
     &                         aec(i),aecd(i),aed(i),aem(i),aep(i),
     &                         aer(i),aes(i),aelf(i),aeg(i),aex(i)
         end if
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine amberyze  --  parameter format for Amber setup  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "amberyze" prints the force field parameters in a format needed
c     by the Amber setup protocol for using AMOEBA within Amber; this
c     is essentially the "paramyze" format from TINKER 4.3
c
c
      subroutine amberyze (active)
      use sizes
      use angang
      use angbnd
      use angpot
      use angtor
      use atomid
      use atoms
      use bitor
      use bndstr
      use charge
      use dipole
      use improp
      use imptor
      use iounit
      use korbs
      use ktrtor
      use kvdws
      use math
      use mpole
      use opbend
      use opdist
      use piorbs
      use pistuf
      use pitors
      use polar
      use polgrp
      use potent
      use solute
      use strbnd
      use strtor
      use tors
      use tortor
      use units
      use urey
      use vdw
      use vdwpot
      implicit none
      integer i,j,k,m
      integer ia,ib,ic
      integer id,ie,ig
      integer itx,ity
      integer ixaxe
      integer iyaxe
      integer izaxe
      integer fold(6)
      real*8 bla,blc
      real*8 sbavg
      real*8 radj,rad4j
      real*8 ampli(6)
      real*8 phase(6)
      real*8 mpl(13)
      logical header
      logical active(*)
c
c
c     number of each type of interaction and site
c
      if (n .ne. 0) then
         write (iout,10)
   10    format (/,' Total Numbers of Atoms and Interactions :')
         write (iout,20)  n
   20    format (/,' Atoms in System',11x,i15)
      end if
      if (norbit .ne. 0) then
         write (iout,30)  norbit
   30    format (' Pisystem Atoms',12x,i15)
      end if
      if (nbond .ne. 0) then
         write (iout,40)  nbond
   40    format (' Bond Stretches',12x,i15)
      end if
      if (nbpi .ne. 0) then
         write (iout,50)  nbpi
   50    format (' Conjugated Pi-Bonds',7x,i15)
      end if
      if (nangle .ne. 0) then
         write (iout,60)  nangle
   60    format (' Angle Bends',15x,i15)
      end if
      if (nstrbnd .ne. 0) then
         write (iout,70)  nstrbnd
   70    format (' Stretch-Bends',13x,i15)
      end if
      if (nurey .ne. 0) then
         write (iout,80)  nurey
   80    format (' Urey-Bradley',14x,i15)
      end if
      if (nangang .ne. 0) then
         write (iout,90)  nangang
   90    format (' Angle-Angles',14x,i15)
      end if
      if (nopbend .ne. 0) then
         write (iout,100)  nopbend
  100    format (' Out-of-Plane Bends',8x,i15)
      end if
      if (nopdist .ne. 0) then
         write (iout,110)  nopdist
  110    format (' Out-of-Plane Distances',4x,i15)
      end if
      if (niprop .ne. 0) then
         write (iout,120)  niprop
  120    format (' Improper Dihedrals',8x,i15)
      end if
      if (nitors .ne. 0) then
         write (iout,130)  nitors
  130    format (' Improper Torsions',9x,i15)
      end if
      if (ntors .ne. 0) then
         write (iout,140)  ntors
  140    format (' Torsional Angles',10x,i15)
      end if
      if (npitors .ne. 0) then
         write (iout,150)  npitors
  150    format (' Pi-Orbital Torsions',7x,i15)
      end if
      if (nstrtor .ne. 0) then
         write (iout,160)  nstrtor
  160    format (' Stretch-Torsions',10x,i15)
      end if
      if (nangtor .ne. 0) then
         write (iout,170)  nangtor
  170    format (' Angle-Torsions',12x,i15)
      end if
      if (ntortor .ne. 0) then
         write (iout,180)  ntortor
  180    format (' Torsion-Torsions',10x,i15)
      end if
      if (nion .ne. 0) then
         write (iout,190)  nion
  190    format (' Atomic Partial Charges',4x,i15)
      end if
      if (ndipole .ne. 0) then
         write (iout,200)  ndipole
  200    format (' Bond Dipole Moments',7x,i15)
      end if
      if (npole .ne. 0) then
         write (iout,210)  npole
  210    format (' Polarizable Multipoles',4x,i15)
      end if
c
c     parameters used for molecular mechanics atom types
c
      header = .true.
      do i = 1, n
         if (active(i)) then
            if (header) then
               header = .false.
               write (iout,220)
  220          format (/,' Atom Type Definition Parameters :',
     &                 //,3x,'Atom',2x,'Symbol',2x,'Type',
     &                    2x,'Class',2x,'Atomic',3x,'Mass',
     &                    2x,'Valence',2x,'Description',/)
            end if
            write (iout,230)  i,name(i),type(i),class(i),atomic(i),
     &                        mass(i),valence(i),story(i)
  230       format (i6,5x,a3,2i7,i6,f10.3,i5,5x,a24)
         end if
      end do
c
c     parameters used for van der Waals interactions
c
      if (use_vdw) then
         header = .true.
         k = 0
         do i = 1, n
            if (active(i)) then
               k = k + 1
               if (header) then
                  header = .false.
                  write (iout,240)
  240             format (/,' Van der Waals Parameters :',
     &                    //,10x,'Atom Number',7x,'Radius',
     &                       3x,'Epsilon',3x,'Rad 1-4',
     &                       3x,'Eps 1-4',3x,'Reduction',/)
               end if
               j = class(i)
               if (vdwindex .eq. 'TYPE')  j = type(i)
               if (rad(j).eq.rad4(j) .and. eps(j).eq.eps4(j)) then
                  radj = rad(j)
                  if (radtyp .eq. 'SIGMA')  radj = radj / twosix
                  if (reduct(j) .eq. 0.0d0) then
                     write (iout,250)  k,i,radj,eps(j)
  250                format (i6,3x,i6,7x,2f10.4)
                  else
                     write (iout,260)  k,i,radj,eps(j),reduct(j)
  260                format (i6,3x,i6,7x,2f10.4,22x,f10.4)
                  end if
               else
                  radj = rad(j)
                  rad4j = rad4(j)
                  if (radtyp .eq. 'SIGMA') then
                     radj = radj / twosix
                     rad4j = rad4j / twosix
                  end if
                  if (reduct(j) .eq. 0.0d0) then
                     write (iout,270)  k,i,radj,eps(j),rad4j,eps4(j)
  270                format (i6,3x,i6,7x,2f10.4,1x,2f10.4)
                  else
                     write (iout,280)  k,i,radj,eps(j),rad4j,
     &                                eps4(j),reduct(j)
  280                format (i6,3x,i6,7x,2f10.4,1x,2f10.4,1x,f10.4)
                  end if
               end if
            end if
         end do
      end if
c
c     parameters used for bond stretching interactions
c
      if (use_bond) then
         header = .true.
         do i = 1, nbond
            ia = ibnd(1,i)
            ib = ibnd(2,i)
            if (active(ia) .or. active(ib)) then
               if (header) then
                  header = .false.
                  write (iout,290)
  290             format (/,' Bond Stretching Parameters :',
     &                    //,10x,'Atom Numbers',25x,'KS',7x,'Length',/)
               end if
               write (iout,300)  i,ia,ib,bk(i),bl(i)
  300          format (i6,3x,2i6,19x,f10.3,f10.4)
            end if
         end do
      end if
c
c     parameters used for angle bending interactions
c
      if (use_angle) then
         header = .true.
         do i = 1, nangle
            ia = iang(1,i)
            ib = iang(2,i)
            ic = iang(3,i)
            if (active(ia) .or. active(ib) .or. active(ic)) then
               if (header) then
                  header = .false.
                  write (iout,310)
  310             format (/,' Angle Bending Parameters :',
     &                    //,13x,'Atom Numbers',22x,'KB',
     &                       6x,'Angle',3x,'Fold',4x,'Type',/)
               end if
               if (angtyp(i) .eq. 'HARMONIC') then
                  write (iout,320)  i,ia,ib,ic,ak(i),anat(i)
  320             format (i6,3x,3i6,13x,2f10.3)
               else if (angtyp(i) .eq. 'IN-PLANE') then
                  write (iout,330)  i,ia,ib,ic,ak(i),anat(i)
  330             format (i6,3x,3i6,13x,2f10.3,9x,'In-Plane')
               else if (angtyp(i) .eq. 'IN-PLANE') then
                  write (iout,340)  i,ia,ib,ic,ak(i),anat(i)
  340             format (i6,3x,3i6,13x,2f10.3,9x,'Linear')
               else if (angtyp(i) .eq. 'FOURIER ') then
                  write (iout,350)  i,ia,ib,ic,ak(i),anat(i),afld(i)
  350             format (i6,3x,3i6,13x,2f10.3,f7.1,2x,'Fourier')
               end if
            end if
         end do
      end if
c
c     parameters used for stretch-bend interactions
c
      if (use_strbnd) then
         header = .true.
         do i = 1, nstrbnd
            k = isb(1,i)
            ia = iang(1,k)
            ib = iang(2,k)
            ic = iang(3,k)
            if (active(ia) .or. active(ib) .or. active(ic)) then
               if (header) then
                  header = .false.
                  write (iout,360)
  360             format (/,' Stretch-Bend Parameters :',
     &                    //,13x,'Atom Numbers',11x,'KSB',
     &                       6x,'Angle',3x,'Length1',3x,'Length2',/)
               end if
               bla = 0.0d0
               blc = 0.0d0
               if (isb(2,i) .ne. 0)  bla = bl(isb(2,i))
               if (isb(3,i) .ne. 0)  blc = bl(isb(3,i))
               sbavg = (sbk(1,i)+sbk(2,i)) * 0.5d0
               write (iout,370)  i,ia,ib,ic,sbavg,anat(k),bla,blc
  370          format (i6,3x,3i6,f13.4,3f10.4)
            end if
         end do
      end if
c
c     parameters used for Urey-Bradley interactions
c
      if (use_urey) then
         header = .true.
         do i = 1, nurey
            ia = iury(1,i)
            ic = iury(3,i)
            if (active(ia) .or. active(ic)) then
               if (header) then
                  header = .false.
                  write (iout,380)
  380             format (/,' Urey-Bradley Parameters :',
     &                    //,10x,'Atom Numbers',24x,'KUB',
     &                       4x,'Distance',/)
               end if
               write (iout,390)  i,ia,ic,uk(i),ul(i)
  390          format (i6,3x,2i6,13x,f16.4,f10.4)
            end if
         end do
      end if
c
c     parameters used for out-of-plane bend interactions
c
      if (use_opbend) then
         header = .true.
         do i = 1, nopbend
            k = iopb(i)
            ia = iang(1,k)
            ib = iang(2,k)
            ic = iang(3,k)
            id = iang(4,k)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,400)
  400             format (/,' Out-of-Plane Bending Parameters :',
     &                    //,17x,'Atom Numbers',19x,'KOPB',/)
               end if
               opbk(i) = opbk(i) * (opbunit/0.02191418d0)
               write (iout,410)  i,id,ib,ia,ic,opbk(i)
  410          format (i6,3x,4i6,9x,f10.4)
            end if
         end do
      end if
c
c     parameters used for out-of-plane distance interactions
c
      if (use_opdist) then
         header = .true.
         do i = 1, nopdist
            ia = iopd(1,i)
            ib = iopd(2,i)
            ic = iopd(3,i)
            id = iopd(4,i)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,420)
  420             format (/,' Out-of-Plane Distance Parameters :',
     &                    //,17x,'Atom Numbers',19x,'KOPD',/)
               end if
               write (iout,430)  i,ia,ib,ic,id,opdk(i)
  430          format (i6,3x,4i6,9x,f10.3)
            end if
         end do
      end if
c
c     parameters used for improper dihedral interactions
c
      if (use_improp) then
         header = .true.
         do i = 1, niprop
            ia = iiprop(1,i)
            ib = iiprop(2,i)
            ic = iiprop(3,i)
            id = iiprop(4,i)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,440)
  440             format (/,' Improper Dihedral Parameters :',
     &                    //,17x,'Atom Numbers',19x,'KID',
     &                       4x,'Dihedral',/)
               end if
               write (iout,450)  i,ia,ib,ic,id,kprop(i),vprop(i)
  450          format (i6,3x,4i6,9x,2f10.4)
            end if
         end do
      end if
c
c     parameters used for improper torsion interactions
c
      if (use_imptor) then
         header = .true.
         do i = 1, nitors
            ia = iitors(1,i)
            ib = iitors(2,i)
            ic = iitors(3,i)
            id = iitors(4,i)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,460)
  460             format (/,' Improper Torsion Parameters :',
     &                    //,17x,'Atom Numbers',11x,
     &                       'Amplitude, Phase and Periodicity',/)
               end if
               j = 0
               if (itors1(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 1
                  ampli(j) = itors1(1,i)
                  phase(j) = itors1(2,i)
               end if
               if (itors2(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 2
                  ampli(j) = itors2(1,i)
                  phase(j) = itors2(2,i)
               end if
               if (itors3(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 3
                  ampli(j) = itors3(1,i)
                  phase(j) = itors3(2,i)
               end if
               if (j .eq. 0) then
                  write (iout,470)  i,ia,ib,ic,id
  470             format (i6,3x,4i6)
               else if (j .eq. 1) then
                  write (iout,480)  i,ia,ib,ic,id,
     &                              ampli(1),phase(1),fold(1)
  480             format (i6,3x,4i6,10x,f10.3,f8.1,i4)
               else if (j .eq. 2) then
                  write (iout,490)  i,ia,ib,ic,id,(ampli(k),
     &                              phase(k),fold(k),k=1,j)
  490             format (i6,3x,4i6,2x,2(f10.3,f6.1,i4))
               else
                  write (iout,500)  i,ia,ib,ic,id,(ampli(k),
     &                              nint(phase(k)),fold(k),k=1,j)
  500             format (i6,3x,4i6,4x,3(f8.3,i4,'/',i1))
               end if
            end if
         end do
      end if
c
c     parameters used for torsional interactions
c
      if (use_tors) then
         header = .true.
         do i = 1, ntors
            ia = itors(1,i)
            ib = itors(2,i)
            ic = itors(3,i)
            id = itors(4,i)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,510)
  510             format (/,' Torsional Angle Parameters :',
     &                    //,17x,'Atom Numbers',11x,
     &                       'Amplitude, Phase and Periodicity',/)
               end if
               j = 0
               if (tors1(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 1
                  ampli(j) = tors1(1,i)
                  phase(j) = tors1(2,i)
               end if
               if (tors2(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 2
                  ampli(j) = tors2(1,i)
                  phase(j) = tors2(2,i)
               end if
               if (tors3(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 3
                  ampli(j) = tors3(1,i)
                  phase(j) = tors3(2,i)
               end if
               if (tors4(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 4
                  ampli(j) = tors4(1,i)
                  phase(j) = tors4(2,i)
               end if
               if (tors5(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 5
                  ampli(j) = tors5(1,i)
                  phase(j) = tors5(2,i)
               end if
               if (tors6(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 6
                  ampli(j) = tors6(1,i)
                  phase(j) = tors6(2,i)
               end if
               if (j .eq. 0) then
                  write (iout,520)  i,ia,ib,ic,id
  520             format (i6,3x,4i6)
               else
                  write (iout,530)  i,ia,ib,ic,id,(ampli(k),
     &                              nint(phase(k)),fold(k),k=1,j)
  530             format (i6,3x,4i6,4x,6(f8.3,i4,'/',i1))
               end if
            end if
         end do
      end if
c
c     parameters used for pi-system torsion interactions
c
      if (use_pitors) then
         header = .true.
         do i = 1, npitors
            ia = ipit(1,i)
            ib = ipit(2,i)
            ic = ipit(3,i)
            id = ipit(4,i)
            ie = ipit(5,i)
            ig = ipit(6,i)
            if (active(ia) .or. active(ib) .or. active(ic) .or.
     &             active(id) .or. active(ie) .or. active(ig)) then
               if (header) then
                  header = .false.
                  write (iout,540)
  540             format (/,' Pi-Orbital Torsion Parameters :',
     &                    //,10x,'Atom Numbers',19x,'Amplitude',/)
               end if
               write (iout,550)  i,ic,id,kpit(i)
  550          format (i6,3x,2i6,19x,f10.4)
            end if
         end do
      end if
c
c     parameters used for stretch-torsion interactions; this is
c     the "old" stretch-torsion format for central bond only
c
      if (use_strtor) then
         header = .true.
         do i = 1, nstrtor
            k = ist(1,i)
            ia = itors(1,k)
            ib = itors(2,k)
            ic = itors(3,k)
            id = itors(4,k)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,560)
  560             format (/,' Stretch-Torsion Parameters :',
     &                    //,17x,'Atom Numbers',10x,'Length',
     &                       5x,'Torsion Terms',/)
               end if
               j = 0
               if (kst(4,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 1
                  ampli(j) = kst(4,i)
                  phase(j) = tors1(2,k)
               end if
               if (kst(5,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 2
                  ampli(j) = kst(5,i)
                  phase(j) = tors2(2,k)
               end if
               if (kst(6,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 3
                  ampli(j) = kst(6,i)
                  phase(j) = tors3(2,k)
               end if
               write (iout,570)  i,ia,ib,ic,id,bl(ist(3,i)),
     &                           (ampli(k),nint(phase(k)),
     &                           fold(k),k=1,j)
  570          format (i6,3x,4i6,2x,f10.4,1x,3(f8.3,i4,'/',i1))
            end if
         end do
      end if
c
c     parameters used for angle-torsion interactions; this term
c     is currently not implemented in Amber
c
      if (use_angtor) then
         header = .true.
         do i = 1, nangtor
            k = iat(1,i)
            ia = itors(1,k)
            ib = itors(2,k)
            ic = itors(3,k)
            id = itors(4,k)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,580)
  580             format (/,' Angle-Torsion Parameters :',
     &                    //,17x,'Atom Numbers',10x,'Length',
     &                       5x,'Torsion Terms',/)
               end if
               write (iout,590)  i,ia,ib,ic,id
  590          format (i6,3x,4i6)
            end if
         end do
      end if
c
c     parameters used for torsion-torsion interactions
c
      if (use_tortor) then
         header = .true.
         do i = 1, ntortor
            k = itt(1,i)
            ia = ibitor(1,k)
            ib = ibitor(2,k)
            ic = ibitor(3,k)
            id = ibitor(4,k)
            ie = ibitor(5,k)
            if (active(ia) .or. active(ib) .or. active(ic) .or.
     &                active(id) .or. active(ie)) then
               if (header) then
                  header = .false.
                  write (iout,600)
  600             format (/,' Torsion-Torsion Parameters :',
     &                    //,20x,'Atom Numbers',18x,'Spline Grid',/)
               end if
               j = itt(2,i)
               write (iout,610)  i,ia,ib,ic,id,ie,tnx(j),tny(j)
  610          format (i6,3x,5i6,10x,2i6)
               do m = 1, tnx(j)*tny(j)
                  itx = (m-1)/tnx(j) + 1
                  ity = m - (itx-1)*tny(j)
                  write (iout,620)  ttx(itx,j),tty(ity,j),tbf(m,j)
  620             format (9x,2f12.1,5x,f12.5)
               end do
            end if
         end do
      end if
c
c     parameters used for atomic partial charges
c
      if (use_charge .or. use_chgdpl) then
         header = .true.
         do i = 1, nion
            ia = iion(i)
            ib = jion(i)
            ic = kion(i)
            if (active(ia) .or. active(ic)) then
               if (header) then
                  header = .false.
                  write (iout,630)
  630             format (/,' Atomic Partial Charge Parameters :',
     &                    /,45x,'Neighbor',3x,'Cutoff',
     &                    /,10x,'Atom Number',13x,'Charge',
     &                       7x,'Site',6x,'Site',/)
               end if
               if (ia.eq.ib .and. ia.eq.ic) then
                  write (iout,640)  i,ia,pchg(i)
  640             format (i6,3x,i6,15x,f10.4)
               else
                  write (iout,650)  i,ia,pchg(i),ib,ic
  650             format (i6,3x,i6,15x,f10.4,5x,i6,4x,i6)
               end if
            end if
         end do
      end if
c
c     parameters used for bond dipole moments
c
      if (use_dipole .or. use_chgdpl) then
         header = .true.
         do i = 1, ndipole
            ia = idpl(1,i)
            ib = idpl(2,i)
            if (active(ia) .or. active(ib)) then
               if (header) then
                  header = .false.
                  write (iout,660)
  660             format (/,' Bond Dipole Moment Parameters :',
     &                    //,10x,'Atom Numbers',22x,'Dipole',
     &                       3x,'Position',/)
               end if
               write (iout,670)  i,ia,ib,bdpl(i),sdpl(i)
  670          format (i6,3x,2i6,19x,f10.4,f10.3)
            end if
         end do
      end if
c
c     parameters used for atomic multipole moments
c
      if (use_mpole) then
         header = .true.
         do i = 1, npole
            ia = ipole(i)
            if (active(ia)) then
               if (header) then
                  header = .false.
                  write (iout,680)
  680             format (/,' Atomic Multipole Parameters :',
     &                    //,12x,'Atom',4x,'Coordinate Frame',
     &                       ' Definition',7x,'Multipole Moments',/)
               end if
               izaxe = zaxis(i)
               ixaxe = xaxis(i)
               iyaxe = yaxis(i)
               if (iyaxe .lt. 0)  iyaxe = -iyaxe
               mpl(1) = pole(1,i)
               do j = 2, 4
                  mpl(j) = pole(j,i) / bohr
               end do
               do j = 5, 13
                  mpl(j) = 3.0d0 * pole(j,i) / bohr**2
               end do
               if (izaxe .eq. 0) then
                  write (iout,690)  i,ia,0,0,polaxe(i),
     &                              (mpl(j),j=1,5),mpl(8),mpl(9),
     &                              (mpl(j),j=11,13)
  690             format (i6,3x,i6,1x,2i7,10x,a8,2x,f9.5,/,50x,3f9.5,
     &                    /,50x,f9.5,/,50x,2f9.5,/,50x,3f9.5)
               else if (ixaxe .eq. 0) then
                  write (iout,700)  i,ia,izaxe,0,polaxe(i),
     &                              (mpl(j),j=1,5),mpl(8),mpl(9),
     &                              (mpl(j),j=11,13)
  700             format (i6,3x,i6,1x,2i7,10x,a8,2x,f9.5,/,50x,3f9.5,
     &                    /,50x,f9.5,/,50x,2f9.5,/,50x,3f9.5)
               else  if (iyaxe .eq. 0) then
                  write (iout,710)  i,ia,izaxe,ixaxe,polaxe(i),
     &                              (mpl(j),j=1,5),mpl(8),mpl(9),
     &                              (mpl(j),j=11,13)
  710             format (i6,3x,i6,1x,2i7,10x,a8,2x,f9.5,/,50x,3f9.5,
     &                    /,50x,f9.5,/,50x,2f9.5,/,50x,3f9.5)
               else
                  write (iout,720)  i,ia,izaxe,ixaxe,iyaxe,polaxe(i),
     &                              (mpl(j),j=1,5),mpl(8),mpl(9),
     &                              (mpl(j),j=11,13)
  720             format (i6,3x,i6,1x,3i7,3x,a8,2x,f9.5,/,50x,3f9.5,
     &                    /,50x,f9.5,/,50x,2f9.5,/,50x,3f9.5)
               end if
            end if
         end do
      end if
c
c     parameters used for dipole polarizability
c
      if (use_polar) then
         header = .true.
         do i = 1, npole
            ia = ipole(i)
            if (active(ia)) then
               if (header) then
                  header = .false.
                  write (iout,730)
  730             format (/,' Dipole Polarizability Parameters :',
     &                    //,10x,'Atom Number',9x,'Alpha',8x,
     &                       'Polarization Group',/)
               end if
               write (iout,740)  i,ia,polarity(i),
     &                           (ip11(j,ia),j=1,np11(ia))
  740          format (i6,3x,i6,10x,f10.4,5x,20i6)
            end if
         end do
      end if
c
c     parameters used for empirical solvation
c
      if (use_solv) then
         header = .true.
         k = 0
         do i = 1, n
            if (active(i)) then
               k = k + 1
               if (header) then
                  header = .false.
                  write (iout,750)
  750             format (/,' Empirical Solvation Parameters :',
     &                    //,10x,'Atom Number',13x,'Radius',
     &                       3x,'ASP Value',/)
               end if
               write (iout,760)  k,i,rsolv(i),asolv(i)
  760          format (i6,3x,i6,15x,2f10.4)
            end if
         end do
      end if
c
c     parameters used for conjugated pisystem atoms
c
      if (use_orbit) then
         header = .true.
         do i = 1, norbit
            ia = iorbit(i)
            j = class(ia)
            if (header) then
               header = .false.
               write (iout,770)
  770          format (/,' Conjugated Pi-Atom Parameters :',
     &                 //,10x,'Atom Number',14x,'Nelect',
     &                    6x,'Ionize',4x,'Repulsion',/)
            end if
            write (iout,780)  i,ia,electron(j),ionize(j),repulse(j)
  780       format (i6,3x,i6,17x,f8.1,3x,f10.4,2x,f10.4)
         end do
      end if
c
c     parameters used for conjugated pibond interactions
c
      if (use_orbit) then
         header = .true.
         do i = 1, nbpi
            ia = ibpi(2,i)
            ib = ibpi(3,i)
            if (header) then
               header = .false.
               write (iout,790)
  790          format (/,' Conjugated Pi-Bond Parameters :',
     &                 //,10x,'Atom Numbers',21x,'K Slope',
     &                    3x,'L Slope',/)
            end if
            write (iout,800)  i,ia,ib,kslope(i),lslope(i)
  800       format (i6,3x,2i6,19x,2f10.4)
         end do
      end if
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine connyze  --  connected atom list analysis  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "connyze" prints information onconnected atoms as lists
c     of all atom pairs that are 1-2 through 1-5 interactions
c
c
      subroutine connyze (active)
      use sizes
      use atoms
      use couple
      use iounit
      implicit none
      integer i,j,k
      integer ntot
      integer ntot2,ntot3
      integer ntot4,ntot5
      logical active(*)
c
c
c     count the number of 1-2 through 1-5 interatomic pairs
c
      ntot2 = 0
      ntot3 = 0
      ntot4 = 0
      ntot5 = 0
      do i = 1, n
         ntot2 = ntot2 + n12(i)
         ntot3 = ntot3 + n13(i)
         ntot4 = ntot4 + n14(i)
         ntot5 = ntot5 + n15(i)
      end do
      ntot2 = ntot2 / 2
      ntot3 = ntot3 / 2
      ntot4 = ntot4 / 2
      ntot5 = ntot5 / 2
      ntot = ntot2 + ntot3 + ntot4 + ntot5
      if (ntot .ne. 0) then
         write (iout,10)
   10    format (/,' Total Number of Pairwise Atomic Interactions :',/)
      end if
      if (ntot2 .ne. 0) then
         write (iout,20)  ntot2
   20    format (' Number of 1-2 Pairs',7x,i15)
      end if
      if (ntot3 .ne. 0) then
         write (iout,30)  ntot3
   30    format (' Number of 1-3 Pairs',7x,i15)
      end if
      if (ntot4 .ne. 0) then
         write (iout,40)  ntot4
   40    format (' Number of 1-4 Pairs',7x,i15)
      end if
      if (ntot5 .ne. 0) then
         write (iout,50)  ntot5
   50    format (' Number of 1-5 Pairs',7x,i15)
      end if
c
c     generate and print the 1-2 connected atomic interactions
c
      if (ntot2 .ne. 0) then
         write (iout,60)
   60    format (/,' List of 1-2 Connected Atomic Interactions :',/)
         do i = 1, n
            if (active(i)) then
               do j = 1, n12(i)
                  k = i12(j,i)
                  if (active(k)) then
                     if (i .lt. k) then
                        write (iout,70)  i,k
   70                   format (2i8)
                     end if
                  end if
               end do
            end if
         end do
      end if
c
c     generate and print the 1-3 connected atomic interactions
c
      if (ntot3 .ne. 0) then
         write (iout,80)
   80    format (/,' List of 1-3 Connected Atomic Interactions :',/)
         do i = 1, n
            if (active(i)) then
               do j = 1, n13(i)
                  k = i13(j,i)
                  if (active(k)) then
                     if (i .lt. k) then
                        write (iout,90)  i,k
   90                   format (2i8)
                     end if
                  end if
               end do
            end if
         end do
      end if
c
c     generate and print the 1-4 connected atomic interactions
c
      if (ntot4 .ne. 0) then
         write (iout,100)
  100    format (/,' List of 1-4 Connected Atomic Interactions :',/)
         do i = 1, n
            if (active(i)) then
               do j = 1, n14(i)
                  k = i14(j,i)
                  if (active(k)) then
                     if (i .lt. k) then
                        write (iout,110)  i,k
  110                   format (2i8)
                     end if
                  end if
               end do
            end if
         end do
      end if
c
c     generate and print the 1-5 connected atomic interactions
c
      if (ntot5 .ne. 0) then
         write (iout,120)
  120    format (/,' List of 1-5 Connected Atomic Interactions :',/)
         do i = 1, n
            if (active(i)) then
               do j = 1, n15(i)
                  k = i15(j,i)
                  if (active(k)) then
                     if (i .lt. k) then
                        write (iout,130)  i,k
  130                   format (2i8)
                     end if
                  end if
               end do
            end if
         end do
      end if
      return
      end
