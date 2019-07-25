c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  program xtalfit  --  fit parameters to structure & energy  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "xtalfit" determines optimized van der Waals and electrostatic
c     parameters by fitting to crystal structures, lattice energies,
c     and dimer structures and interaction energies
c
c
      program xtalfit
      use sizes
      use bound
      use boxes
      use files
      use iounit
      use molcul
      use potent
      use vdwpot
      use xtals
      implicit none
      integer i,ixtal
      integer atom1,atom2
      integer nresid,prmtyp
      real*8 grdmin
      real*8, allocatable :: xx(:)
      real*8, allocatable :: resid(:)
      real*8, allocatable :: g(:)
      real*8, allocatable :: xlo(:)
      real*8, allocatable :: xhi(:)
      real*8, allocatable :: fjac(:,:)
      logical exist,query
      character*5 vindex
      character*16 label(7)
      character*240 record
      character*240 string
      external xtalerr,xtalwrt
c
c
c     initialize some variables to be used during fitting
c
      call initial
      nvary = 0
      nresid = 0
c
c     print informational header about available parameters
c
      write (iout,10)
   10 format (/,' The Following Parameters can be Fit for',
     &           ' each Atom Type :',
     &        //,4x,'(1) Van der Waals Atomic Radius',
     &        /,4x,'(2) Van der Waals Well Depth',
     &        /,4x,'(3) Hydrogen Atom Reduction Factor',
     &        /,4x,'(4) Atomic Partial Charge',
     &        /,4x,'(5) Bond Dipole Moment Magnitude',
     &        /,4x,'(6) Bond Dipole Moment Position',
     &        /,4x,'(7) Atomic Polarizability')
c
c     get types of potential parameters to be optimized
c
      query = .true.
      do while (query)
         prmtyp = -1
         atom1 = 0
         atom2 = 0
         call nextarg (string,exist)
         if (exist)  read (string,*,err=20,end=20)  prmtyp
         call nextarg (string,exist)
         if (exist)  read (string,*,err=20,end=20)  atom1
         call nextarg (string,exist)
         if (exist)  read (string,*,err=20,end=20)  atom2
   20    continue
         if (prmtyp .ne. 0) then
            prmtyp = 0
            write (iout,30)
   30       format (/,' Enter Parameter Type then Atom Class',
     &                 ' or Type(s) :  ',$)
            read (input,40)  record
   40       format (a240)
            read (record,*,err=50,end=50)  prmtyp,atom1,atom2
   50       continue
         end if
         if (prmtyp .eq. 0) then
            query = .false.
         else
            query = .true.
            nvary = nvary + 1
            ivary(nvary) = prmtyp
            vary(1,nvary) = atom1
            if (prmtyp.eq.5 .or. prmtyp.eq.6) then
               vary(1,nvary) = min(atom1,atom2)
               vary(2,nvary) = max(atom1,atom2)
            end if
         end if
      end do
c
c     get termination criterion as RMS gradient over parameters
c
      grdmin = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=60,end=60)  grdmin
   60 continue
      if (grdmin .le. 0.0d0) then
         write (iout,70)
   70    format (/,' Enter RMS Gradient Termination Criterion',
     &              ' [0.1] :  ',$)
         read (input,80)  grdmin
   80    format (f20.0)
      end if
      if (grdmin .le. 0.0d0)  grdmin = 0.1d0
c
c     get the number of structures to use in optimization
c
      nxtal = 0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=90,end=90)  nxtal
   90 continue
      if (nxtal .le. 0) then
         write (iout,100)
  100    format (/,' Enter Number of Structures to be Used [1] :  ',$)
         read (input,110)  nxtal
  110    format (i10)
      end if
c
c     check for too few or too many molecular structures
c
      if (nxtal .eq. 0)  nxtal = 1
      if (nxtal .gt. maxref) then
         write (iout,120)
  120    format (/,' XTALFIT  --  Too many Structures,',
     &              ' Increase the Value of MAXREF')
         call fatal
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (xx(nvary))
c
c     get coordinates and parameters for current structure
c
      do ixtal = 1, nxtal
         call initial
         call getxyz
         call mechanic
c
c     get ideal value for intermolecular or lattice energy
c
         e0_lattice = 0.0d0
         query = .true.
         call nextarg (string,exist)
         if (exist) then
            read (string,*,err=130,end=130)  e0_lattice
            query = .false.
         end if
  130    continue
         if (query) then
            write (iout,140)
  140       format (/,' Enter Target Elat/Einter Value and Weight',
     &                 ' [<CR> to omit] :  ',$)
            read (input,150)  e0_lattice
  150       format (f20.0)
         end if
         if (e0_lattice .gt. 0.0d0)  e0_lattice = -e0_lattice
c
c     set the types of residuals for use in optimization
c
         do i = 1, 6
            iresid(nresid+i) = ixtal
         end do
         if (use_bounds) then
            rsdtyp(nresid+1) = 'Force a-Axis'
            rsdtyp(nresid+2) = 'Force b-Axis'
            rsdtyp(nresid+3) = 'Force c-Axis'
            rsdtyp(nresid+4) = 'Force Alpha'
            rsdtyp(nresid+5) = 'Force Beta'
            rsdtyp(nresid+6) = 'Force Gamma'
         else
            rsdtyp(nresid+1) = 'Force Mol1 X'
            rsdtyp(nresid+2) = 'Force Mol1 Y'
            rsdtyp(nresid+3) = 'Force Mol1 Z'
            rsdtyp(nresid+4) = 'Force Mol2 X'
            rsdtyp(nresid+5) = 'Force Mol2 Y'
            rsdtyp(nresid+6) = 'Force Mol2 Z'
         end if
         nresid = nresid + 6
c
c     print molecules per structure, energy and dipole values
c
         write (iout,160)  ixtal,filename(1:35),nmol
  160    format (/,' File Name of Target Structure',i4,' :',8x,a35,
     &           /,' Number of Molecules per Structure :',i13)
         if (e0_lattice .ne. 0.0d0) then
            nresid = nresid + 1
            iresid(nresid) = ixtal
            if (use_bounds) then
               rsdtyp(nresid) = 'Lattice Energy'
            else
               rsdtyp(nresid) = 'E Intermolecular'
            end if
            write (iout,170)  e0_lattice
  170       format (' Target E-Lattice or E-Inter Value :  ',f13.2)
         end if
c
c     set the initial values of the parameters
c
         call xtalprm ('STORE',ixtal,xx)
      end do
c
c     turn off all local interactions and extra terms
c
      call potoff
      use_vdw = .true.
      use_charge = .true.
      use_chgdpl = .true.
      use_dipole = .true.
      use_mpole = .true.
      use_polar = .true.
c
c     types of variables for use in optimization
c
      label(1) = 'Atomic Radius'
      label(2) = 'Well Depth'
      label(3) = 'H Reduction'
      label(4) = 'Partial Charge'
      label(5) = 'Dipole Magnitude'
      label(6) = 'Dipole Position'
      label(7) = 'Polarizability'
      do i = 1, nvary
         vartyp(i) = label(ivary(i))
      end do
      vindex = 'Class'
      if (vdwindex .eq. 'TYPE ')  vindex = 'Type '
c
c     print the initial parameter values
c
      write (iout,180)
  180 format (/,' Initial Values of the Parameters :',/)
      do i = 1, nvary
         if (ivary(i) .le. 3) then
            write (iout,190)  i,vartyp(i),vindex,vary(1,i),xx(i)
  190       format (3x,'(',i2,')',2x,a16,4x,'Atom ',a5,i5,4x,f12.4)
         else if (ivary(i) .ne. 6) then
            write (iout,200)  i,vartyp(i),vary(1,i),xx(i)
  200       format (3x,'(',i2,')',2x,a16,4x,'Atom Type ',i5,4x,f12.4)
         else
            write (iout,210)  i,vartyp(i),vary(1,i),vary(2,i),xx(i)
  210       format (3x,'(',i2,')',2x,a16,4x,'Bond Type ',2i5,f12.4)
         end if
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (resid(nresid))
      allocate (g(nvary))
      allocate (xlo(nvary))
      allocate (xhi(nvary))
      allocate (fjac(nresid,nvary))
c
c     set upper and lower bounds based on the parameter type
c
      do i = 1, nvary
         if (ivary(i).eq.4 .or. ivary(i).eq.5) then
            xlo(i) = xx(i) - 0.5d0
            xhi(i) = xx(i) + 0.5d0
         else
            xlo(i) = 0.5d0 * xx(i)
            xhi(i) = 1.5d0 * xx(i)
         end if
      end do
c
c     use nonlinear least squares to refine the parameters
c
      call square (nvary,nresid,xlo,xhi,xx,resid,g,fjac,
     &                  grdmin,xtalerr,xtalwrt)
c
c     perform deallocation of some local arrays
c
      deallocate (xlo)
      deallocate (xhi)
      deallocate (fjac)
c
c     print final values of parameters and scaled derivatives
c
      write (iout,220)
  220 format (/,' Final Values of Parameters and Scaled',
     &           ' Derivatives :',/)
      do i = 1, nvary
         if (ivary(i) .le. 3) then
            write (iout,230)  i,vartyp(i),vindex,vary(1,i),xx(i),g(i)
  230       format (3x,'(',i2,')',2x,a16,4x,'Atom ',a5,i5,2x,2f14.4)
         else if (ivary(i) .ne. 6) then
            write (iout,240)  i,vartyp(i),vary(1,i),xx(i),g(i)
  240       format (3x,'(',i2,')',2x,a16,4x,'Atom Type ',i5,2x,2f14.4)
         else
            write (iout,250)  i,vartyp(i),vary(1,i),vary(2,i),xx(i),g(i)
  250       format (3x,'(',i2,')',2x,a16,4x,'Bond Type ',2i5,
     &                 f11.4,f14.4)
         end if
      end do
c
c     print final values of the individual residual functions
c
      write (iout,260)
  260 format (/,' Final Residual Error Function Values :',/)
      do i = 1, nresid
         if (i .lt. 100) then
            write (iout,270)  i,rsdtyp(i),iresid(i),resid(i)
  270       format (3x,'(',i2,')',2x,a16,6x,'Structure',i4,4x,f12.4)
         else
            write (iout,280)  i,rsdtyp(i),iresid(i),resid(i)
  280       format (2x,'(',i3,')',2x,a16,6x,'Structure',i4,4x,f12.4)
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xx)
      deallocate (resid)
      deallocate (g)
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine xtalprm  --  energy/optimization conversion  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "xtalprm" stores or retrieves a molecular structure; used to
c     make a previously stored structure the active structure, or to
c     store a structure for later use
c
c     the current version only provides for intermolecular potential
c     energy terms
c
c
      subroutine xtalprm (mode,ixtal,xx)
      use sizes
      use atoms
      use atomid
      use bound
      use boxes
      use charge
      use dipole
      use files
      use fracs
      use inform
      use kvdws
      use molcul
      use mpole
      use polar
      use vdw
      use vdwpot
      use xtals
      implicit none
      integer i,j,k
      integer ixtal,prmtyp
      integer atom1,atom2
      real*8 rd,ep,sixth
      real*8 xmid,ymid,zmid
      real*8 e0_lattices(maxref)
      real*8 xx(*)
      logical first
      character*5 mode
      save e0_lattices
      save first
      data first  / .true. /
c
c
c     save or restore the key values for the current crystal
c
      if (mode .eq. 'STORE') then
         call makeref (ixtal)
      else if (mode .eq. 'RESET') then
         call getref (ixtal)
         call basefile (filename)
         silent = .true.
         call mechanic
         silent = .false.
         if (use_bounds)  call bounds
      end if
c
c     perform dynamic allocation of some global arrays
c
      if (mode .eq. 'RESET') then
         if (first) then
            first = .false.
            allocate (xfrac(nmol))
            allocate (yfrac(nmol))
            allocate (zfrac(nmol))
         end if
      end if
c
c     coordinates of molecular centers of mass
c
      if (mode .eq. 'RESET') then
         do i = 1, nmol
            xmid = 0.0d0
            ymid = 0.0d0
            zmid = 0.0d0
            do j = imol(1,i), imol(2,i)
               k = kmol(j)
               xmid = xmid + x(k)*mass(k)
               ymid = ymid + y(k)*mass(k)
               zmid = zmid + z(k)*mass(k)
            end do
            zmid = zmid / gamma_term
            ymid = (ymid - zmid*beta_term) / gamma_sin
            xmid = xmid - ymid*gamma_cos - zmid*beta_cos
            xfrac(i) = xmid / (xbox * molmass(i))
            yfrac(i) = ymid / (ybox * molmass(i))
            zfrac(i) = zmid / (zbox * molmass(i))
         end do
      end if
c
c     values of ideal intermolecular or lattice energy
c
      if (mode .eq. 'STORE') then
         e0_lattices(ixtal) = e0_lattice
      else if (mode .eq. 'RESET') then
         e0_lattice = e0_lattices(ixtal)
      end if
c
c     store or reset values of the optimization variables
c
      do j = 1, nvary
         prmtyp = ivary(j)
         atom1 = vary(1,j)
         if (prmtyp .eq. 1) then
            if (mode .eq. 'STORE') then
               xx(j) = rad(atom1)
            else if (mode .eq. 'RESET') then
               rad(atom1) = xx(j)
               do i = 1, maxclass
                  if (rad(i).eq.0.0d0 .and. rad(atom1).eq.0.0d0) then
                     rd = 0.0d0
                  else if (radrule(1:10) .eq. 'ARITHMETIC') then
                     rd = rad(i) + rad(atom1)
                  else if (radrule(1:9) .eq. 'GEOMETRIC') then
                     rd = 2.0d0 * sqrt(rad(i) * rad(atom1))
                  else if (radrule(1:10) .eq. 'CUBIC-MEAN') then
                     rd = 2.0d0 * (rad(i)**3+rad(atom1)**3)
     &                       / (rad(i)**2+rad(atom1)**2)
                  else
                     rd = rad(i) + rad(atom1)
                  end if
                  radmin(i,atom1) = rd
                  radmin(atom1,i) = rd
               end do
            end if
         else if (prmtyp .eq. 2) then
            if (mode .eq. 'STORE') then
               xx(j) = eps(atom1)
            else if (mode .eq. 'RESET') then
               eps(atom1) = abs(xx(j))
               do i = 1, maxclass
                  if (eps(i).eq.0.0d0 .and. eps(atom1).eq.0.0d0) then
                     ep = 0.0d0
                  else if (epsrule(1:10) .eq. 'ARITHMETIC') then
                     ep = 0.5d0 * (eps(i) + eps(atom1))
                  else if (epsrule(1:9) .eq. 'GEOMETRIC') then
                     ep = sqrt(eps(i) * eps(atom1))
                  else if (epsrule(1:8) .eq. 'HARMONIC') then
                     ep = 2.0d0 * (eps(i)*eps(atom1))
     &                       / (eps(i)+eps(atom1))
                  else if (epsrule(1:3) .eq. 'HHG') then
                     ep = 4.0d0 * (eps(i)*eps(atom1))
     &                      / (sqrt(eps(i))+sqrt(eps(atom1)))**2
                  else
                     ep = sqrt(eps(i) * eps(atom1))
                  end if
                  epsilon(i,atom1) = ep
                  epsilon(atom1,i) = ep
               end do
            end if
         else if (prmtyp .eq. 3) then
            if (mode .eq. 'STORE') then
               do i = 1, n
                  if (class(i) .eq. atom1) then
                     xx(j) = kred(i)
                     goto 10
                  end if
               end do
            else if (mode .eq. 'RESET') then
               do i = 1, n
                  if (class(i) .eq. atom1)  kred(i) = xx(j)
               end do
            end if
         else if (prmtyp .eq. 4) then
            if (mode .eq. 'STORE') then
               do i = 1, nion
                  if (type(iion(i)) .eq. atom1) then
                     xx(j) = pchg(i)
                     goto 10
                  end if
               end do
            else if (mode .eq. 'RESET') then
               do i = 1, nion
                  if (type(iion(i)) .eq. atom1)  pchg(i) = xx(j)
               end do
            end if
         else if (prmtyp .eq. 5) then
            atom2 = vary(2,j)
            if (mode .eq. 'STORE') then
               do i = 1, ndipole
                  if (type(idpl(1,i)).eq.atom1 .and.
     &                type(idpl(2,i)).eq.atom2) then
                     xx(j) = bdpl(i)
                     goto 10
                  end if
               end do
            else if (mode .eq. 'RESET') then
               do i = 1, ndipole
                  if (type(idpl(1,i)).eq.atom1 .and.
     &                type(idpl(2,i)).eq.atom2)  bdpl(i) = xx(j)
               end do
            end if
         else if (prmtyp .eq. 6) then
            atom2 = vary(2,j)
            if (mode .eq. 'STORE') then
               do i = 1, ndipole
                  if (type(idpl(1,i)).eq.atom1 .and.
     &                type(idpl(2,i)).eq.atom2) then
                     xx(j) = sdpl(i)
                     goto 10
                  end if
               end do
            else if (mode .eq. 'RESET') then
               do i = 1, ndipole
                  if (type(idpl(1,i)).eq.atom1 .and.
     &                type(idpl(2,i)).eq.atom2)  sdpl(i) = xx(j)
               end do
            end if
         else if (prmtyp .eq. 7) then
            if (mode .eq. 'STORE') then
               do i = 1, npole
                  if (type(ipole(i)) .eq. atom1) then
                     xx(j) = polarity(i)
                     goto 10
                  end if
               end do
            else if (mode .eq. 'RESET') then
               sixth = 1.0d0 / 6.0d0
               do i = 1, npole
                  if (type(ipole(i)) .eq. atom1) then
                     polarity(i) = xx(j)
                     if (thole(i) .ne. 0.0d0)  pdamp(i) = xx(j)**sixth
                  end if
               end do
            end if
         end if
   10    continue
      end do
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine xtalerr  --  error function for xtalfit  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "xtalerr" computes an error function value derived from
c     lattice energies, dimer intermolecular energies and the
c     gradient with respect to structural parameters
c
c
      subroutine xtalerr (nvaried,nresid,xx,resid)
      use sizes
      use atoms
      use boxes
      use bound
      use charge
      use dipole
      use energi
      use limits
      use math
      use molcul
      use mpole
      use polar
      use vdw
      use xtals
      implicit none
      integer i,k,ixtal
      integer nresid,nvaried
      real*8 energy,eps,temp
      real*8 e,e0
      real*8 e_monomer
      real*8 e_lattice
      real*8 dmol,big
      real*8 e1,e2,e3
      real*8 e4,e5,e6
      real*8 g1,g2,g3
      real*8 g4,g5,g6
      real*8 xx(*)
      real*8 resid(*)
c
c
c     zero out number of residuals and set numerical step size
c
      nresid = 0
      eps = 1.0d-4
c
c     set force field parameter values and find the base energy
c
      do ixtal = 1, nxtal
         call xtalprm ('RESET',ixtal,xx)
         e = energy ()
         e0 = ev + ec + ecd + ed + em + ep
c
c     perturb crystal lattice parameters and compute energies
c
         if (use_bounds) then
            temp = xbox
            xbox = xbox + eps
            call xtalmove
            e = energy ()
            e1 = ev + ec + ecd + ed + em + ep
            xbox = temp
            temp = ybox
            ybox = ybox + eps
            call xtalmove
            e = energy ()
            e2 = ev + ec + ecd + ed + em + ep
            ybox = temp
            temp = zbox
            zbox = zbox + eps
            call xtalmove
            e = energy ()
            e3 = ev + ec + ecd + ed + em + ep
            zbox = temp
            temp = alpha
            alpha = alpha + radian*eps
            call xtalmove
            e = energy ()
            e4 = ev + ec + ecd + ed + em + ep
            alpha = temp
            temp = beta
            beta = beta + radian*eps
            call xtalmove
            e = energy ()
            e5 = ev + ec + ecd + ed + em + ep
            beta = temp
            temp = gamma
            gamma = gamma + radian*eps
            call xtalmove
            e = energy ()
            e6 = ev + ec + ecd + ed + em + ep
            gamma = temp
            call xtalmove
c
c     translate dimer component molecules and compute energies
c
         else
            do i = imol(1,1), imol(2,1)
               k = kmol(i)
               x(k) = x(k) + eps
            end do
            e = energy ()
            e1 = ev + ec + ecd + ed + em + ep
            do i = imol(1,1), imol(2,1)
               k = kmol(i)
               x(k) = x(k) - eps
            end do
            do i = imol(1,1), imol(2,1)
               k = kmol(i)
               y(k) = y(k) + eps
            end do
            e = energy ()
            e2 = ev + ec + ecd + ed + em + ep
            do i = imol(1,1), imol(2,1)
               k = kmol(i)
               y(k) = y(k) - eps
            end do
            do i = imol(1,1), imol(2,1)
               k = kmol(i)
               z(k) = z(k) + eps
            end do
            e = energy ()
            e3 = ev + ec + ecd + ed + em + ep
            do i = imol(1,1), imol(2,1)
               k = kmol(i)
               z(k) = z(k) - eps
            end do
            do i = imol(1,1), imol(2,1)
               k = kmol(i)
               x(k) = x(k) + eps
            end do
            e = energy ()
            e4 = ev + ec + ecd + ed + em + ep
            do i = imol(1,2), imol(2,2)
               k = kmol(i)
               x(k) = x(k) - eps
            end do
            do i = imol(1,2), imol(2,2)
               k = kmol(i)
               y(k) = y(k) + eps
            end do
            e = energy ()
            e5 = ev + ec + ecd + ed + em + ep
            do i = imol(1,2), imol(2,2)
               k = kmol(i)
               y(k) = y(k) - eps
            end do
            do i = imol(1,2), imol(2,2)
               k = kmol(i)
               z(k) = z(k) + eps
            end do
            e = energy ()
            e6 = ev + ec + ecd + ed + em + ep
            do i = imol(1,2), imol(2,2)
               k = kmol(i)
               z(k) = z(k) - eps
            end do
         end if
c
c     get the gradient with respect to structure perturbations
c
         g1 = (e1 - e0) / eps
         nresid = nresid + 1
         resid(nresid) = g1
         g2 = (e2 - e0) / eps
         nresid = nresid + 1
         resid(nresid) = g2
         g3 = (e3 - e0) / eps
         nresid = nresid + 1
         resid(nresid) = g3
         g4 = (e4 - e0) / eps
         nresid = nresid + 1
         resid(nresid) = g4
         g5 = (e5 - e0) / eps
         nresid = nresid + 1
         resid(nresid) = g5
         g6 = (e6 - e0) / eps
         nresid = nresid + 1
         resid(nresid) = g6
c
c     setup to compute properties of monomer from crystal
c
         if (use_bounds) then
            n = n / nmol
            nvdw = nvdw / nmol
            nion = nion / nmol
            ndipole = ndipole / nmol
            npole = npole / nmol
            npolar = npolar / nmol
            use_bounds = .false.
            use_replica = .false.
            use_ewald = .false.
            big = 1.0d12
            vdwcut = big
            vdwtaper = big
            chgcut = big
            chgtaper = big
            dplcut = big
            dpltaper = big
            mpolecut = big
            mpoletaper = big
c
c     compute the intermolecular or crystal lattice energy
c
            e = energy ()
            e_monomer = ev + ec + ecd + ed + em + ep
            dmol = dble(nmol)
            e_lattice = (e0 - dmol*e_monomer) / dmol
         else
            e_monomer = 0.0d0
            e_lattice = e0
         end if
c
c     compute residual due to intermolecular or lattice energy;
c     weight energies more heavily, since there are fewer of them
c
         if (e0_lattice .ne. 0.0d0) then
            nresid = nresid + 1
            resid(nresid) = e_lattice - e0_lattice
            if (ixtal .le. 11) then
               resid(nresid) = 3.0d0 * resid(nresid)
            else
               resid(nresid) = 10.0d0 * resid(nresid)
            end if
         end if
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine xtalmove  --  translation of rigid molecules  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "xtalmove" converts fractional to Cartesian coordinates for
c     rigid molecules during optimization of force field parameters
c
c
      subroutine xtalmove
      use sizes
      use atoms
      use atomid
      use boxes
      use fracs
      use molcul
      implicit none
      integer i,j,k
      integer init,stop
      real*8 weigh
      real*8 xmid,ymid,zmid
      real*8, allocatable :: xoff(:)
      real*8, allocatable :: yoff(:)
      real*8, allocatable :: zoff(:)
c
c
c     get values for fractional coordinate interconversion
c
      call lattice
c
c     perform dynamic allocation of some local arrays
c
      allocate (xoff(n))
      allocate (yoff(n))
      allocate (zoff(n))
c
c     locate the center of mass of each molecule
c
      do i = 1, nmol
         init = imol(1,i)
         stop = imol(2,i)
         xmid = 0.0d0
         ymid = 0.0d0
         zmid = 0.0d0
         do j = init, stop
            k = kmol(j)
            weigh = mass(k)
            xmid = xmid + x(k)*weigh
            ymid = ymid + y(k)*weigh
            zmid = zmid + z(k)*weigh
         end do
         weigh = molmass(i)
         xmid = xmid / weigh
         ymid = ymid / weigh
         zmid = zmid / weigh
c
c     save atomic coordinates relative to center of mass
c
         do j = init, stop
            k = kmol(j)
            xoff(k) = x(k) - xmid
            yoff(k) = y(k) - ymid
            zoff(k) = z(k) - zmid
         end do
c
c     convert fractional center of mass to Cartesian coordinates
c
         xmid = xfrac(i)*xbox + yfrac(i)*ybox*gamma_cos
     &             + zfrac(i)*zbox*beta_cos
         ymid = yfrac(i)*ybox*gamma_sin + zfrac(i)*zbox*beta_term
         zmid = zfrac(i)*zbox*gamma_term
c
c     translate coordinates via offset from center of mass
c
         do j = init, stop
            k = kmol(j)
            x(k) = xoff(k) + xmid
            y(k) = yoff(k) + ymid
            z(k) = zoff(k) + zmid
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xoff)
      deallocate (yoff)
      deallocate (zoff)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine xtalwrt  --  output optimization parameters  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "xtalwrt" prints intermediate results during fitting of
c     force field parameters to structures and energies
c
c
      subroutine xtalwrt (niter,nresid,xx,gs,resid)
      use iounit
      use vdwpot
      use xtals
      implicit none
      integer i,niter
      integer nresid
      real*8 xx(*)
      real*8 gs(*)
      real*8 resid(*)
      character*5 vindex
c
c
c     print the values of parameters and scaled derivatives
c
      vindex = 'Class'
      if (vdwindex .eq. 'TYPE ')  vindex = 'Type '
      write (iout,10)  niter
   10 format (/,' Parameters and Scaled Derivatives at',
     &          ' Iteration',i4,' :',/)
      do i = 1, nvary
         if (ivary(i) .le. 3) then
            write (iout,20)  i,vartyp(i),vindex,vary(1,i),xx(i),gs(i)
   20       format (3x,'(',i2,')',2x,a16,4x,'Atom ',a5,i5,2x,2f14.4)
         else if (ivary(i) .ne. 6) then
            write (iout,30)  i,vartyp(i),vary(1,i),xx(i),gs(i)
   30       format (3x,'(',i2,')',2x,a16,4x,'Atom Type ',i5,2x,2f14.4)
         else
            write (iout,40)  i,vartyp(i),vary(1,i),vary(2,i),xx(i),gs(i)
   40       format (3x,'(',i2,')',2x,a16,4x,'Bond Type ',2i5,
     &                 f11.4,f14.4)
         end if
      end do
c
c     print the values of the individual residual functions
c
      write (iout,50)  niter
   50 format (/,' Residual Error Function Values at Iteration',
     &           i4,' :',/)
      do i = 1, nresid
         if (i .lt. 100) then
            write (iout,60)  i,rsdtyp(i),iresid(i),resid(i)
   60       format (3x,'(',i2,')',2x,a16,6x,'Structure',i4,4x,f12.4)
         else
            write (iout,70)  i,rsdtyp(i),iresid(i),resid(i)
   70       format (2x,'(',i3,')',2x,a16,6x,'Structure',i4,4x,f12.4)
         end if
      end do
      write (iout,80)
   80 format ()
      return
      end
