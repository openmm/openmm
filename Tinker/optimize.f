c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  program optimize  --  variable metric Cartesian optimizer  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "optimize" performs energy minimization in Cartesian coordinate
c     space using an optimally conditioned variable metric method
c
c
      program optimize
      use sizes
      use atoms
      use files
      use inform
      use iounit
      use keys
      use scales
      use usage
      implicit none
      integer i,j,imin,nvar
      integer next,freeunit
      real*8 minimum,optimiz1
      real*8 grdmin,gnorm,grms
      real*8 energy,eps
      real*8, allocatable :: xx(:)
      real*8, allocatable :: derivs(:,:)
      logical exist,analytic
      character*20 keyword
      character*240 minfile
      character*240 record
      character*240 string
      external energy
      external optimiz1
      external optsave
c
c
c     set up the structure and mechanics calculation
c
      call initial
      call getxyz
      call mechanic
c
c     use either analytical or numerical gradients
c
      analytic = .true.
      eps = 0.00001d0
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
c     get termination criterion as RMS gradient per atom
c
      grdmin = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=20,end=20)  grdmin
   20 continue
      if (grdmin .le. 0.0d0) then
         write (iout,30)
   30    format (/,' Enter RMS Gradient per Atom Criterion',
     &              ' [0.01] :  ',$)
         read (input,40)  grdmin
   40    format (f20.0)
      end if
      if (grdmin .le. 0.0d0)  grdmin = 0.01d0
c
c     write out a copy of coordinates for later update
c
      imin = freeunit ()
      minfile = filename(1:leng)//'.xyz'
      call version (minfile,'new')
      open (unit=imin,file=minfile,status='new')
      call prtxyz (imin)
      close (unit=imin)
      outfile = minfile
c
c     set scaling parameter for function and derivative values;
c     use square root of median eigenvalue of typical Hessian
c
      set_scale = .true.
      nvar = 0
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               nvar = nvar + 1
               scale(nvar) = 12.0d0
            end do
         end if
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (xx(nvar))
      allocate (derivs(3,n))
c
c     scale the coordinates of each active atom
c
      nvar = 0
      do i = 1, n
         if (use(i)) then
            nvar = nvar + 1
            xx(nvar) = x(i) * scale(nvar)
            nvar = nvar + 1
            xx(nvar) = y(i) * scale(nvar)
            nvar = nvar + 1
            xx(nvar) = z(i) * scale(nvar)
         end if
      end do
c
c     make the call to the optimization routine
c
      call ocvm (nvar,xx,minimum,grdmin,optimiz1,optsave)
c
c     unscale the final coordinates for active atoms
c
      nvar = 0
      do i = 1, n
         if (use(i)) then
            nvar = nvar + 1
            x(i) = xx(nvar) / scale(nvar)
            nvar = nvar + 1
            y(i) = xx(nvar) / scale(nvar)
            nvar = nvar + 1
            z(i) = xx(nvar) / scale(nvar)
         end if
      end do
c
c     compute the final function and RMS gradient values
c
      if (analytic) then
         call gradient (minimum,derivs)
      else
         minimum = energy ()
         call numgrad (energy,derivs,eps)
      end if
      gnorm = 0.0d0
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               gnorm = gnorm + derivs(j,i)**2
            end do
         end if
      end do
      gnorm = sqrt(gnorm)
      grms = gnorm / sqrt(dble(nvar/3))
c
c     perform deallocation of some local arrays
c
      deallocate (xx)
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
      call prtxyz (imin)
      close (unit=imin)
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  function optimiz1  --  energy and gradient for optimize  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "optimiz1" is a service routine that computes the energy and
c     gradient for optimally conditioned variable metric optimization
c     in Cartesian coordinate space
c
c
      function optimiz1 (xx,g)
      use sizes
      use atoms
      use scales
      use usage
      implicit none
      integer i,nvar
      real*8 optimiz1,e
      real*8 energy,eps
      real*8 xx(*)
      real*8 g(*)
      real*8, allocatable :: derivs(:,:)
      logical analytic
      external energy
c
c
c     use either analytical or numerical gradients
c
      analytic = .true.
      eps = 0.00001d0
c
c     translate optimization parameters to atomic coordinates
c
      nvar = 0
      do i = 1, n
         if (use(i)) then
            nvar = nvar + 1
            x(i) = xx(nvar) / scale(nvar)
            nvar = nvar + 1
            y(i) = xx(nvar) / scale(nvar)
            nvar = nvar + 1
            z(i) = xx(nvar) / scale(nvar)
         end if
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (derivs(3,n))
c
c     compute and store the energy and gradient
c
      if (analytic) then
         call gradient (e,derivs)
      else
         e = energy ()
         call numgrad (energy,derivs,eps)
      end if
      optimiz1 = e
c
c     store Cartesian gradient as optimization gradient
c
      nvar = 0
      do i = 1, n
         if (use(i)) then
            nvar = nvar + 1
            g(nvar) = derivs(1,i) / scale(nvar)
            nvar = nvar + 1
            g(nvar) = derivs(2,i) / scale(nvar)
            nvar = nvar + 1
            g(nvar) = derivs(3,i) / scale(nvar)
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (derivs)
      return
      end
