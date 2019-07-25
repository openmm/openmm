c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  program minirot  --  low storage BFGS torsional optimizer  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "minirot" performs an energy minimization in torsional
c     angle space using a low storage BFGS nonlinear optimization
c
c
      program minirot
      use sizes
      use files
      use inform
      use iounit
      use keys
      use math
      use omega
      use scales
      use zcoord
      implicit none
      integer i,imin,next
      integer freeunit
      real*8 minimum,minirot1
      real*8 grdmin,grms,gnorm
      real*8, allocatable :: xx(:)
      real*8, allocatable :: derivs(:)
      logical exist
      character*20 keyword
      character*240 minfile
      character*240 record
      character*240 string
      external minirot1
      external optsave
c
c
c     set up the molecular mechanics calculation
c
      call initial
      call getint
      call mechanic
      call initrot
c
c     search the keywords for output frequency parameters
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:9) .eq. 'PRINTOUT ') then
            read (string,*,err=10,end=10)  iprint
         else if (keyword(1:9) .eq. 'WRITEOUT ') then
            read (string,*,err=10,end=10)  iwrite
         end if
   10    continue
      end do
c
c     get termination criterion as RMS torsional gradient
c
      grdmin = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=20,end=20)  grdmin
   20 continue
      if (grdmin .le. 0.0d0) then
         write (iout,30)
   30    format (/,' Enter RMS Gradient per Torsion Criterion',
     &              ' [0.01] :  ',$)
         read (input,40)  grdmin
   40    format (f20.0)
      end if
      if (grdmin .eq. 0.0d0)  grdmin = 0.01d0
c
c     write out a copy of coordinates for later update
c
      imin = freeunit ()
      minfile = filename(1:leng)//'.int'
      call version (minfile,'new')
      open (unit=imin,file=minfile,status='new')
      call prtint (imin)
      close (unit=imin)
      outfile = minfile
c
c     set scaling parameter for function and derivative values;
c     use square root of median eigenvalue of typical Hessian
c
      set_scale = .true.
      do i = 1, nomega
         scale(i) = 5.0d0
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (xx(nomega))
c
c     get optimization parameters as scaled dihedral angles
c
      do i = 1, nomega
         xx(i) = dihed(i) * scale(i)
      end do
c
c     make the call to the optimization routine
c
      call lbfgs (nomega,xx,minimum,grdmin,minirot1,optsave)
c
c     unscale the final dihedral angle values
c
      do i = 1, nomega
         dihed(i) = xx(i) / scale(i)
         ztors(zline(i)) = dihed(i) * radian
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xx)
c
c     perform dynamic allocation of some local arrays
c
      allocate (derivs(nomega))
c
c     compute the final function and RMS gradient values
c
      call gradrot (minimum,derivs)
      gnorm = 0.0d0
      do i = 1, nomega
         gnorm = gnorm + derivs(i)**2
      end do
      gnorm = sqrt(gnorm)
      grms = gnorm / sqrt(dble(nomega))
c
c     perform deallocation of some local arrays
c
      deallocate (derivs)
c
c     write out the final function and gradient values
c
      if (digits .ge. 8) then
         if (grms .gt. 1.0d-8) then
            write (iout,50)  minimum,grms,gnorm
   50       format (/,' Final Function Value :',2x,f20.8,
     &              /,' Final RMS Gradient :',4x,f20.8,
     &              /,' Final Gradient Norm :',3x,f20.8)
         else
            write (iout,60)  minimum,grms,gnorm
   60       format (/,' Final Function Value :',2x,f20.8,
     &              /,' Final RMS Gradient :',4x,d20.8,
     &              /,' Final Gradient Norm :',3x,d20.8)
         end if
      else if (digits .ge. 6) then
         if (grms .gt. 1.0d-6) then
            write (iout,70)  minimum,grms,gnorm
   70       format (/,' Final Function Value :',2x,f18.6,
     &              /,' Final RMS Gradient :',4x,f18.6,
     &              /,' Final Gradient Norm :',3x,f18.6)
         else
            write (iout,80)  minimum,grms,gnorm
   80       format (/,' Final Function Value :',2x,f18.6,
     &              /,' Final RMS Gradient :',4x,d18.6,
     &              /,' Final Gradient Norm :',3x,d18.6)
         end if
      else
         if (grms .gt. 1.0d-4) then
            write (iout,90)  minimum,grms,gnorm
   90       format (/,' Final Function Value :',2x,f16.4,
     &              /,' Final RMS Gradient :',4x,f16.4,
     &              /,' Final Gradient Norm :',3x,f16.4)
         else
            write (iout,100)  minimum,grms,gnorm
  100       format (/,' Final Function Value :',2x,f16.4,
     &              /,' Final RMS Gradient :',4x,d16.4,
     &              /,' Final Gradient Norm :',3x,d16.4)
         end if
      end if
c
c     write the final coordinates into a file
c
      imin = freeunit ()
      open (unit=imin,file=minfile,status='old')
      rewind (unit=imin)
      call prtint (imin)
      close (unit=imin)
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  function minirot1  --  energy and gradient for minirot  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "minirot1" is a service routine that computes the energy
c     and gradient for a low storage BFGS nonlinear optimization
c     in torsional angle space
c
c
      function minirot1 (xx,g)
      use sizes
      use math
      use omega
      use scales
      use zcoord
      implicit none
      integer i
      real*8 minirot1,e
      real*8 xx(*)
      real*8 g(*)
      real*8, allocatable :: derivs(:)
c
c
c     translate optimization variables into dihedrals
c
      do i = 1, nomega
         dihed(i) = xx(i) / scale(i)
         ztors(zline(i)) = dihed(i) * radian
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (derivs(nomega))
c
c     get coordinates, then compute energy and gradient
c
      call makexyz
      call gradrot (e,derivs)
      minirot1 = e
c
c     store torsional gradient as optimization gradient
c
      do i = 1, nomega
         g(i) = derivs(i) / scale(i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (derivs)
      return
      end
