c
c
c     ##############################################################
c     ##  COPYRIGHT (C) 1997 by Rohit Pappu & Jay William Ponder  ##
c     ##                   All Rights Reserved                    ##
c     ##############################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  program pssrigid  --  smoothing & search over rigid bodies  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "pssrigid" implements the potential smoothing plus search method
c     for global optimization for a set of rigid bodies
c
c     literature reference:
c
c     J. Kostrowicki and H. A. Scheraga, "Application of the Diffusion
c     Equation Method for Global Optimization to Oligopeptides", Journal
c     of Physical Chemistry, 96, 7442-7449 (1992)
c
c     S. Nakamura, H. Hirose, M. Ikeguchi and J. Doi, "Conformational
c     Energy Minimization Using a Two-Stage Method", Journal of Physical
c     Chemistry, 99, 8374-8378 (1995)
c
c
      program pssrigid
      use sizes
      use atoms
      use files
      use group
      use inform
      use iounit
      use math
      use molcul
      use refer
      use rigid
      use warp
      implicit none
      integer i,j,k,ixyz
      integer nvar,lext
      integer npoint,neigen
      integer next,freeunit
      real*8 minimum,grdmin
      real*8 srchmax,rms
      real*8 pssrgd1,deform0
      real*8 ratio,sigmoid
      real*8, allocatable :: xx(:)
      logical exist
      logical use_local
      character*1 answer
      character*7 ext
      character*240 xyzfile
      character*240 record
      character*240 string
      external pssrgd1
      external optsave
c
c
c     set up the structure, mechanics calculation and smoothing
c
      call initial
      call getxyz
      use_smooth = .true.
      use_dem = .true.
      call mechanic
c
c     get rigid body coordinates and save the Cartesian coordinates
c
      use_rigid = .true.
      call orient
      call makeref (1)
c
c     set maximum deformation value and disable coordinate dumps
c
      deform0 = deform
      iwrite = 0
c
c     get the number of points along the deformation schedule
c
      npoint = -1
      call nextarg (string,exist)
      if (exist)  read (string,*,err=10,end=10)  npoint
   10 continue
      if (npoint .lt. 0) then
         write (iout,20)
   20    format (/,' Enter the Number of Steps for the PSS Schedule',
     &              ' [100] :  ',$)
         read (input,30)  npoint
   30    format (i10)
         if (npoint .le. 0)  npoint = 100
      end if
c
c     decide whether to use the local search procedure
c
      use_local = .false.
      call nextarg (answer,exist)
      if (.not. exist) then
         write (iout,40)
   40    format (/,' Use Local Search to Explore Each Smoothing Level',
     &              ' [N] :  ',$)
         read (input,50)  record
   50    format (a240)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .eq. 'Y')  use_local = .true.
c
c     get the number of eigenvectors to use for the local search
c
      if (use_local) then
         neigen = -1
         call nextarg (string,exist)
         if (exist)  read (string,*,err=60,end=60)  neigen
   60    continue
         if (neigen .le. 0) then
            nvar = 6 * (ngrp-1)
            write (iout,70)  nvar
   70       format (/,' Enter the Number of Directions for Local',
     &                 ' Search [',i2,'] :  ',$)
            read (input,80)  neigen
   80       format (i10)
            if (neigen .gt. nvar)  neigen = nvar
         end if
      end if
c
c     get the maximal smoothing level for use of local search
c
      if (use_local) then
         srchmax = -1.0d0
         call nextarg (string,exist)
         if (exist)  read (string,*,err=90,end=90)  srchmax
   90    continue
         if (srchmax .lt. 0.0d0) then
            write (iout,100)
  100       format (/,' Enter the Largest Smoothing Value for Local',
     &                 ' Search [5.0] :  ',$)
            read (input,110)  srchmax
  110       format (f20.0)
            if (srchmax .lt. 0.0d0)  srchmax = 5.0d0
         end if
      end if
c
c     get the termination criterion as RMS gradient per atom
c
      grdmin = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=120,end=120)  grdmin
  120 continue
      if (grdmin .le. 0.0d0) then
         write (iout,130)
  130    format (/,' Enter RMS Gradient per Rigid Body Criterion',
     &              ' [0.0001] :  ',$)
         read (input,140)  grdmin
  140    format (f20.0)
      end if
      if (grdmin .le. 0.0d0)  grdmin = 0.0001d0
c
c     perform dynamic allocation of some local arrays
c
      allocate (xx(6*ngrp))
c
c     perform PSS iteration by looping over smoothed surfaces
c
      do k = 0, 2*npoint
         ratio = 1.0d0 - dble(abs(npoint-k))/dble(npoint)
         if (nmol .eq. 1) then
            deform = deform0 * ratio**3
         else
            deform = deform0 * sigmoid (12.0d0,ratio)
         end if
c
c     transfer rigid body coordinates to optimization parameters
c
         nvar = 0
         do i = 1, ngrp
            do j = 1, 6
               nvar = nvar + 1
               xx(nvar) = rbc(j,i)
            end do
         end do
c
c     make the call to the variable metric optimization routine
c
         iprint = 1
         call ocvm (nvar,xx,minimum,grdmin,pssrgd1,optsave)
c
c     transfer optimization parameters to rigid body coordinates
c
         nvar = 0
         do i = 1, ngrp
            do j = 1, 6
               nvar = nvar + 1
               rbc(j,i) = xx(nvar)
            end do
         end do
c
c     use normal mode local search to explore adjacent minima
c
         if (use_local) then
            if (deform.le.srchmax .and. k.ge.npoint)
     &         call modergd (neigen,minimum,grdmin)
         end if
c
c     write out final energy function value and smoothing level
c
         write (iout,150)  minimum,deform
  150    format (/,' Final Function Value and Deformation :',2f15.4)
c
c     get Cartesian coordinates and superimpose on reference
c
         call rigidxyz
         if (igrp(1,1).eq.1 .and. igrp(2,ngrp).eq.n)
     &      call impose (n,xref,yref,zref,n,x,y,z,rms)
c
c     write the coordinates of the current minimum to a file
c
         lext = 3
         call numeral (k,ext,lext)
         ixyz = freeunit ()
         xyzfile = filename(1:leng)//'.'//ext(1:lext)
         call version (xyzfile,'new')
         open (unit=ixyz,file=xyzfile,status='new')
         call prtxyz (ixyz)
         close (unit=ixyz)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xx)
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  function pssrgd1  --  energy and gradient values for PSS  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "pssrgd1" is a service routine that computes the energy and
c     gradient during PSS global optimization over rigid bodies
c
c
      function pssrgd1 (xx,g)
      use sizes
      use group
      use math
      use rigid
      implicit none
      integer i,j,nvar
      real*8 pssrgd1,e
      real*8 xx(*)
      real*8 g(*)
      real*8, allocatable :: derivs(:,:)
c
c
c     translate optimization parameters to rigid body coordinates
c
      nvar = 0
      do i = 1, ngrp
         do j = 1, 6
            nvar = nvar + 1
            rbc(j,i) = xx(nvar)
         end do
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (derivs(6,ngrp))
c
c     compute and store the energy and gradient
c
      call rigidxyz
      call gradrgd (e,derivs)
      pssrgd1 = e
c
c     store rigid body gradient as optimization gradient
c
      nvar = 0
      do i = 1, ngrp
         do j = 1, 6
            nvar = nvar + 1
            g(nvar) = derivs(j,i)
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (derivs)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine modergd  --  local search for rigid body PSS  ##
c     ##                                                           ##
c     ###############################################################
c
c
      subroutine modergd (neigen,minimum,grdmin)
      use sizes
      use group
      use iounit
      use rigid
      implicit none
      integer maxrgd
      parameter (maxrgd=6*maxgrp)
      integer i,j,k
      integer neigen,ndoi
      integer nvar,nsearch
      real*8 minimum,grdmin
      real*8 eps,minref,minbest
      real*8, allocatable :: step(:)
      real*8, allocatable :: eigen(:)
      real*8, allocatable :: rorig(:,:)
      real*8, allocatable :: rbest(:,:)
      real*8, allocatable :: vects(:,:)
      logical done
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (step(ngrp))
      allocate (eigen(ngrp))
      allocate (rorig(6,ngrp))
      allocate (rbest(6,ngrp))
      allocate (vects(6*ngrp,6*ngrp))
c
c     set parameters related to the local search procedure
c
      done = .false.
      eps = 0.0001d0
      minref = minimum
      minbest = minimum
      ndoi = 0
      nvar = 6 * ngrp
      do i = 1, ngrp
         do j = 1, 6
            rorig(j,i) = rbc(j,i)
         end do
      end do
c
c     find local minimum along each of the steepest directions
c
      do while (.not. done)
         ndoi = ndoi + 1
         write (iout,10)  ndoi,minref
   10    format (/,' Normal Mode Search :',8x,'Iteration',i4,
     &              6x,'Energy',f12.4,/)
         call rigidxyz
         call eigenrgd (eigen,vects)
c
c     search both directions along each eigenvector in turn
c
         nsearch = 0
         do i = 1, neigen
            do k = 1, nvar
               step(k) = vects(k,nvar-i+1)
            end do
            do k = 1, ngrp
               do j = 1, 6
                  rbc(j,k) = rorig(j,k)
               end do
            end do
            nsearch = nsearch + 1
            call climbrgd (nsearch,minimum,step,grdmin)
            if (minimum .lt. minbest) then
               minbest = minimum
               do k = 1, ngrp
                  do j = 1, 6
                     rbest(j,k) = rbc(j,k)
                  end do
               end do
            end if
            do k = 1, nvar
               step(k) = -vects(k,nvar-i+1)
            end do
            do k = 1, ngrp
               do j = 1, 6
                  rbc(j,k) = rorig(j,k)
               end do
            end do
            nsearch = nsearch + 1
            call climbrgd (nsearch,minimum,step,grdmin)
            if (minimum .lt. minbest) then
               minbest = minimum
               do k = 1, ngrp
                  do j = 1, 6
                     rbest(j,k) = rbc(j,k)
                  end do
               end do
            end if
         end do
c
c     check for convergence of the local search procedure
c
         if (minbest .lt. minref-eps) then
            done = .false.
            minref = minbest
            do k = 1, ngrp
               do j = 1, 6
                  rorig(j,k) = rbest(j,k)
               end do
            end do
         else
            done = .true.
            minimum = minref
            do k = 1, ngrp
               do j = 1, 6
                  rbc(j,k) = rorig(j,k)
               end do
            end do
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (step)
      deallocate (eigen)
      deallocate (rorig)
      deallocate (rbest)
      deallocate (vects)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine eigenrgd  --  rigid body Hessian eigenvectors  ##
c     ##                                                            ##
c     ################################################################
c
c
      subroutine eigenrgd (eigen,vects)
      use sizes
      use atoms
      use group
      implicit none
      integer maxrgd
      parameter (maxrgd=6*maxgrp)
      integer i,j
      integer ihess,nvar
      real*8 vnorm
      real*8 eigen(*)
      real*8, allocatable :: matrix(:)
      real*8 vects(6*ngrp,*)
      real*8, allocatable :: hrigid(:,:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (matrix(6*ngrp*(6*ngrp+1)/2))
      allocate (hrigid(6*ngrp,6*ngrp))
c
c     compute the Hessian for rigid body motion
c
      call hessrgd (hrigid)
c
c     place Hessian elements into triangular form
c
      nvar = 6 * ngrp
      ihess = 0
      do i = 1, nvar
         do j = i, nvar
            ihess = ihess + 1
            matrix(ihess) = hrigid(i,j)
         end do
      end do
c
c     diagonalize the Hessian to obtain eigenvalues
c
      call diagq (nvar,nvar,matrix,eigen,vects)
c
c     normalize the rigid body Hessian eigenvectors
c
      do i = 1, nvar
         vnorm = 0.0d0
         do j = 1, nvar
            vnorm = vnorm + vects(j,i)**2
         end do
         vnorm = sqrt(vnorm)
         do j = 1, nvar
            vects(j,i) = vects(j,i) / vnorm
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (matrix)
      deallocate (hrigid)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine climbrgd  --  minimum from a PSS local search  ##
c     ##                                                            ##
c     ################################################################
c
c
      subroutine climbrgd (nsearch,minimum,step,grdmin)
      use sizes
      use group
      use iounit
      use math
      use rigid
      implicit none
      integer maxstep
      parameter (maxstep=500)
      integer i,j,nsearch
      integer nvar,kstep,nstep
      real*8 minimum,grdmin
      real*8 big,energy,size
      real*8 estep(0:maxstep)
      real*8 step(*)
      logical done
c
c
c     set the maximum number of steps and the step size
c
      done = .false.
      big = 100000.0d0
      minimum = big
      kstep = 0
      nstep = 65
c     size = 0.1d0
      size = 1.0d0
      nvar = 6 * ngrp
      do i = 1, nvar
         step(i) = size * step(i)
      end do
c
c     scan the search direction for a minimization candidate
c
      do while (.not. done)
         if (kstep .ne. 0) then
            nvar = 0
            do i = 1, ngrp
               do j = 1, 6
                  nvar = nvar + 1
                  rbc(j,i) = rbc(j,i) + step(nvar)
               end do
            end do
         end if
         call rigidxyz
         estep(kstep) = energy ()
         if (kstep.ge.2 .and. estep(kstep).le.10000.0d0) then
            if (estep(kstep) .lt. estep(kstep-2) .and.
     &          estep(kstep-1) .lt. estep(kstep-2)) then
               done = .true.
               nvar = 0
               do i = 1, ngrp
                  do j = 1, 6
                     nvar = nvar + 1
                     rbc(j,i) = rbc(j,i) - step(nvar)
                  end do
               end do
               call rigidxyz
               call localrgd (minimum,grdmin)
               if (minimum .ge. -big) then
                  write (iout,10)  nsearch,kstep-1,minimum
   10             format (4x,'Search Direction',i4,10x,'Step',
     &                       i6,10x,f12.4)
               else
                  minimum = big
                  write (iout,20)  nsearch
   20             format (4x,'Search Direction',i4,36x,'------')
               end if
            end if
         end if
         if (kstep.ge.nstep .and. .not.done) then
            done = .true.
            write (iout,30)  nsearch
   30       format (4x,'Search Direction',i4,36x,'------')
         end if
         kstep = kstep + 1
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine localrgd  --  PSS local search optimization  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "localrgd" is used during the PSS local search procedure
c     to perform a rigid body energy minimization
c
c
      subroutine localrgd (minimum,grdmin)
      use sizes
      use inform
      use group
      use minima
      use rigid
      implicit none
      integer i,j,nvar
      integer oldprt
      real*8 minimum
      real*8 grdmin
      real*8 pssrgd1
      real*8, allocatable :: xx(:)
      logical oldverb
      external pssrgd1
      external optsave
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (xx(6*ngrp))
c
c     transfer rigid body coordinates to optimization parameters
c
      nvar = 0
      do i = 1, ngrp
         do j = 1, 6
            nvar = nvar + 1
            xx(nvar) = rbc(j,i)
         end do
      end do
c
c     make the call to the optimization routine
c
      oldverb = verbose
      oldprt = iprint
      verbose = .false.
      iprint = 0
      call ocvm (nvar,xx,minimum,grdmin,pssrgd1,optsave)
      verbose = oldverb
      iprint = oldprt
c
c     transfer optimization parameters to rigid body coordinates
c
      nvar = 0
      do i = 1, ngrp
         do j = 1, 6
            nvar = nvar + 1
            rbc(j,i) = xx(nvar)
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xx)
      return
      end
