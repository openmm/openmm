c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  program newtrot  --  perform TNCG torsional optimization  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "newtrot" performs an energy minimization in torsional angle
c     space using a truncated Newton conjugate gradient method
c
c
      program newtrot
      use sizes
      use files
      use inform
      use iounit
      use keys
      use math
      use omega
      use zcoord
      implicit none
      integer i,imin,next
      integer freeunit
      real*8 grdmin,gnorm,grms
      real*8 minimum,newtrot1
      real*8, allocatable :: xx(:)
      real*8, allocatable :: derivs(:)
      logical exist
      character*1 answer
      character*6 mode,method
      character*20 keyword
      character*240 minfile
      character*240 record
      character*240 string
      external newtrot1
      external newtrot2
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
c     get the type of optimization algorithm to use
c
      mode = 'AUTO'
      call nextarg (answer,exist)
      if (.not. exist) then
         answer = 'A'
         write (iout,20)  answer
   20    format (/,' Choose Automatic, Newton, TNCG or DTNCG',
     &              ' Method [',a1,'] :  ',$)
         read (input,30)  record
   30    format (a240)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .eq. 'A')  mode = 'AUTO'
      if (answer .eq. 'N')  mode = 'NEWTON'
      if (answer .eq. 'T')  mode = 'TNCG'
      if (answer .eq. 'D')  mode = 'DTNCG'
c
c     get the type of linear equation preconditioning to use
c
      method = 'DIAG'
      call nextarg (answer,exist)
      if (.not. exist) then
         answer = 'D'
         write (iout,40)  answer
   40    format (/,' Precondition via Auto/None/Diag/',
     &              'SSOR/ICCG [',a1,'] :  ',$)
         read (input,50)  record
   50    format (a240)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .eq. 'A')  method = 'AUTO'
      if (answer .eq. 'N')  method = 'NONE'
      if (answer .eq. 'D')  method = 'DIAG'
      if (answer .eq. 'S')  method = 'SSOR'
      if (answer .eq. 'I')  method = 'ICCG'
c
c     get termination criterion as RMS torsional gradient
c
      grdmin = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=60,end=60)  grdmin
   60 continue
      if (grdmin .le. 0.0d0) then
         write (iout,70)
   70    format (/,' Enter RMS Gradient per Torsion Criterion',
     &              ' [0.01] :  ',$)
         read (input,80)  grdmin
   80    format (f20.0)
      end if
      if (grdmin .le. 0.0d0)  grdmin = 0.01d0
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
c     perform dynamic allocation of some local arrays
c
      allocate (xx(nomega))
c
c     translate the initial coordinates
c
      do i = 1, nomega
         xx(i) = dihed(i)
      end do
c
c     make the call to the optimization routine
c
      call tncg (mode,method,nomega,xx,minimum,grdmin,
     &                newtrot1,newtrot2,optsave)
c
c     untranslate the final coordinates
c
      do i = 1, nomega
         dihed(i) = xx(i)
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
            write (iout,90)  minimum,grms,gnorm
   90       format (/,' Final Function Value :',2x,f20.8,
     &              /,' Final RMS Gradient :',4x,f20.8,
     &              /,' Final Gradient Norm :',3x,f20.8)
         else
            write (iout,100)  minimum,grms,gnorm
  100       format (/,' Final Function Value :',2x,f20.8,
     &              /,' Final RMS Gradient :',4x,d20.8,
     &              /,' Final Gradient Norm :',3x,d20.8)
         end if
      else if (digits .ge. 6) then
         if (grms .gt. 1.0d-6) then
            write (iout,110)  minimum,grms,gnorm
  110       format (/,' Final Function Value :',2x,f18.6,
     &              /,' Final RMS Gradient :',4x,f18.6,
     &              /,' Final Gradient Norm :',3x,f18.6)
         else
            write (iout,120)  minimum,grms,gnorm
  120       format (/,' Final Function Value :',2x,f18.6,
     &              /,' Final RMS Gradient :',4x,d18.6,
     &              /,' Final Gradient Norm :',3x,d18.6)
         end if
      else
         if (grms .gt. 1.0d-4) then
            write (iout,130)  minimum,grms,gnorm
  130       format (/,' Final Function Value :',2x,f16.4,
     &              /,' Final RMS Gradient :',4x,f16.4,
     &              /,' Final Gradient Norm :',3x,f16.4)
         else
            write (iout,140)  minimum,grms,gnorm
  140       format (/,' Final Function Value :',2x,f16.4,
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
c     ##  function newtrot1  --  energy and gradient for newtrot  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "newtrot1" is a service routine that computes the energy
c     and gradient for truncated Newton conjugate gradient
c     optimization in torsional angle space
c
c
      function newtrot1 (xx,g)
      use sizes
      use math
      use omega
      use zcoord
      implicit none
      integer i
      real*8 newtrot1,e
      real*8 xx(*)
      real*8 g(*)
      real*8, allocatable :: derivs(:)
c
c
c     translate optimization variables into dihedrals
c
      do i = 1, nomega
         dihed(i) = xx(i)
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
      newtrot1 = e
c
c     store torsional gradient as optimization gradient
c
      do i = 1, nomega
         g(i) = derivs(i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (derivs)
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine newtrot2  --  Hessian values for newtrot  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "newtrot2" is a service routine that computes the sparse
c     matrix Hessian elements for truncated Newton optimization
c     in torsional angle space
c
c
      subroutine newtrot2 (mode,xx,h,hinit,hstop,hindex,hdiag)
      use sizes
      use hescut
      use math
      use omega
      use zcoord
      implicit none
      integer i,j,ihess
      integer hinit(*)
      integer hstop(*)
      integer hindex(*)
      real*8 xx(*)
      real*8 hdiag(*)
      real*8 h(*)
      real*8, allocatable :: hrot(:,:)
      character*4 mode
c
c
c     translate optimization parameters and compute
c     Cartesian coordinates from internal coordinates
c
      if (mode .eq. 'NONE')  return
      do i = 1, nomega
         dihed(i) = xx(i)
         ztors(zline(i)) = dihed(i) * radian
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (hrot(nomega,nomega))
c
c     compute the desired portion of the Hessian
c
      call makexyz
      call hessrot (mode,hrot)
c
c     store the large elements in sparse matrix format
c
      if (mode .eq. 'FULL') then
         ihess = 0
         do i = 1, nomega
            hdiag(i) = hrot(i,i)
            hinit(i) = ihess + 1
            do j = i+1, nomega
               if (abs(hrot(j,i)) .ge. hesscut) then
                  ihess = ihess + 1
                  hindex(ihess) = j
                  h(ihess) = hrot(j,i)
               end if
            end do
            hstop(i) = ihess
         end do
c
c     store only the Hessian matrix diagonal
c
      else if (mode .eq. 'DIAG') then
         do i = 1, nomega
            hdiag(i) = hrot(i,i)
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (hrot)
      return
      end
