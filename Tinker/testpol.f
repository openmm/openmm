c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2012  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  program testpol  --  check convergence of induced dipoles  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "testpol" computes the induced dipole moments for direct
c     polarization, perturbation theory extrapolation (OPT), and
c     for SCF iterations in order to monitor convergence
c
c
      program testpol
      use sizes
      use atoms
      use bound
      use inform
      use iounit
      use limits
      use minima
      use polar
      use polpot
      use potent
      use rigid
      use units
      use usage
      implicit none
      integer i,j,k
      integer next,kpcg
      integer nvar,iter
      integer itercut
      integer miny
      real*8 sum,epscut
      real*8 ux,uy,uz,u2
      real*8 rdirect
      real*8 rpcg,rxpt
      real*8 eps,delta
      real*8 optfit
      real*8, allocatable :: var(:)
      real*8, allocatable :: yval(:)
      real*8, allocatable :: p(:,:)
      real*8, allocatable :: rms(:)
      real*8, allocatable :: drms(:)
      real*8, allocatable :: tdirect(:)
      real*8, allocatable :: tpcg(:)
      real*8, allocatable :: txpt(:)
      real*8, allocatable :: ddirect(:,:)
      real*8, allocatable :: dpcg(:,:)
      real*8, allocatable :: dxpt(:,:)
      real*8, allocatable :: udirect(:,:)
      real*8, allocatable :: upcg(:,:)
      real*8, allocatable :: uxpt(:,:)
      real*8, allocatable :: ustore(:,:,:)
      logical exist,dofull,done
      character*1 answer
      character*6 savetyp
      character*240 record
      external optfit
c
c
c     get the coordinates and required force field parameters
c
      call initial
      call getxyz
      call mechanic
c
c     check to make sure mutual polarization is being used
c
      if (.not. use_polar) then
         write (iout,10)
   10    format (/,' TESTPOL  --  Induced Dipole Polarization Model',
     &              ' is not in Use')
         call fatal
      end if
c
c     decide whether to output results by gradient component
c
      dofull = .true.
      if (n .gt. 100) then
         dofull = .false.
         call nextarg (answer,exist)
         if (.not. exist) then
            write (iout,20)
   20       format (/,' Output Induced Dipole Components by Atom',
     &                 ' [N] :  ',$)
            read (input,30)  record
   30       format (a240)
            next = 1
            call gettext (record,answer,next)
         end if
         call upcase (answer)
         if (answer .eq. 'Y')  dofull = .true.
      end if
c
c     maintain any periodic boundary conditions
c
      if (use_bounds .and. .not.use_rigid)  call bounds
c
c     generate neighbor lists for iterative SCF solver
c
      savetyp = poltyp
      poltyp = 'MUTUAL'
      call cutoffs
      if (use_list)  call nblist
c
c     set tolerances and rotate multipoles to global frame
c
      maxiter = 100
      itercut = politer
      epscut = poleps
      poleps = 0.0000000001d0
      debug = .false.
      call chkpole
      call rotpole
c
c     perform dynamic allocation of some local arrays
c
      allocate (rms(0:maxiter))
      allocate (drms(maxiter))
      allocate (tdirect(n))
      allocate (tpcg(n))
      allocate (txpt(n))
      allocate (ddirect(3,n))
      allocate (dpcg(3,n))
      allocate (dxpt(3,n))
      allocate (udirect(3,n))
      allocate (upcg(3,n))
      allocate (uxpt(3,n))
      allocate (ustore(3,n,0:maxiter))
c
c     perform dynamic allocation of some global arrays
c
      allocate (uexact(3,n))
c
c     get induced dipoles for direct polarization only
c
      poltyp = 'DIRECT'
      call induce
      do i = 1, n
         do j = 1, 3
            udirect(j,i) = debye * uind(j,i)
            ustore(j,i,0) = udirect(j,i)
         end do
      end do
c
c     print the direct polarization induced dipole moments
c
      if (dofull) then
         write (iout,40)
   40    format (/,' Direct Induced Dipole Moments :',
     &           //,4x,'Atom',15x,'X',13x,'Y',13x,'Z',12x,'Norm',/)
         do i = 1, n
            if (use(i)) then
               ux = udirect(1,i)
               uy = udirect(2,i)
               uz = udirect(3,i)
               u2 = sqrt(ux*ux+uy*uy+uz*uz)
               write (iout,50)  i,ux,uy,uz,u2
   50          format (i8,4x,4f14.6)
            end if
         end do
      end if
c
c     find PCG induced dipoles for increasing iteration counts
c
      poltyp = 'MUTUAL'
      done = .false.
      do k = 1, maxiter
         politer = k
         call induce
         do i = 1, n
            do j = 1, 3
               ustore(j,i,k) = debye * uind(j,i)
            end do
         end do
         sum = 0.0d0
         do i = 1, n
            do j = 1, 3
               sum = sum + (ustore(j,i,k)-ustore(j,i,k-1))**2
            end do
         end do
         drms(k) = sqrt(sum/dble(npolar))
         if (.not. done) then
            if (k.eq.itercut .or. drms(k).lt.epscut) then
               done = .true.
               kpcg = k
               do i = 1, n
                  do j = 1, 3
                     upcg(j,i) = ustore(j,i,k)
                  end do
               end do
            end if
         end if
         if (drms(k) .lt. 0.5d0*poleps)  goto 60
      end do
   60 continue
      maxiter = politer
      do i = 1, n
         do j = 1, 3
            uexact(j,i) = ustore(j,i,maxiter)
         end do
      end do
c
c     print the iterative PCG and exact SCF induced dipoles
c
      if (dofull) then
         write (iout,70)  kpcg
   70    format (/,' Iterative PCG Induced Dipole Moments :',
     &              4x,'(',i3,' Iterations)',
     &           //,4x,'Atom',15x,'X',13x,'Y',13x,'Z',12x,'Norm',/)
         do i = 1, n
            if (use(i)) then
               ux = upcg(1,i)
               uy = upcg(2,i)
               uz = upcg(3,i)
               u2 = sqrt(ux*ux+uy*uy+uz*uz)
               write (iout,80)  i,ux,uy,uz,u2
   80          format (i8,4x,4f14.6)
            end if
         end do
         write (iout,90)
   90    format (/,' Exact SCF Induced Dipole Moments :',
     &           //,4x,'Atom',14x,'X',13x,'Y',13x,'Z',12x,'Norm',/)
         do i = 1, n
            if (use(i)) then
               ux = uexact(1,i)
               uy = uexact(2,i)
               uz = uexact(3,i)
               u2 = sqrt(ux*ux+uy*uy+uz*uz)
               write (iout,100)  i,ux,uy,uz,u2
  100          format (i8,4x,4f14.6)
            end if
         end do
      end if
c
c     get induced dipoles from OPT extrapolation method
c
      poltyp = savetyp
      if (poltyp(1:3) .ne. 'OPT') then
         poltyp = 'OPT4'
         call kpolar
      end if
      call induce
      do i = 1, n
         do j = 1, 3
            uxpt(j,i) = debye * uind(j,i)
         end do
      end do
c
c     print the OPT extrapolation induced dipole moments
c
      if (dofull) then
         write (iout,110)  coptmax
  110    format (/,' Extrapolated OPT',i1,' Induced Dipole Moments :',
     &           //,4x,'Atom',15x,'X',13x,'Y',13x,'Z',12x,'Norm',/)
         do i = 1, n
            if (use(i)) then
               ux = uxpt(1,i)
               uy = uxpt(2,i)
               uz = uxpt(3,i)
               u2 = sqrt(ux*ux+uy*uy+uz*uz)
               write (iout,120)  i,ux,uy,uz,u2
  120          format (i8,4x,4f14.6)
            end if
         end do
      end if
c
c     find differences between approximate and exact dipoles
c
      rdirect = 0.0d0
      rpcg = 0.0d0
      rxpt = 0.0d0
      do i = 1, n
         do j = 1, 3
            ddirect(j,i) = udirect(j,i) - uexact(j,i)
            dpcg(j,i) = upcg(j,i) - uexact(j,i)
            dxpt(j,i) = uxpt(j,i) - uexact(j,i)
         end do
         tdirect(i) = sqrt(ddirect(1,i)**2+ddirect(2,i)**2
     &                           +ddirect(3,i)**2)
         tpcg(i) = sqrt(dpcg(1,i)**2+dpcg(2,i)**2+dpcg(3,i)**2)
         txpt(i) = sqrt(dxpt(1,i)**2+dxpt(2,i)**2+dxpt(3,i)**2)
         rdirect = rdirect + tdirect(i)**2
         rpcg = rpcg + tpcg(i)**2
         rxpt = rxpt + txpt(i)**2
      end do
      rdirect = sqrt(rdirect/dble(n))
      rpcg = sqrt(rpcg/dble(n))
      rxpt = sqrt(rxpt/dble(n))
c
c     print the RMS between approximate and exact dipoles
c
      write (iout,130)  coptmax
  130 format (/,' Approximate vs. Exact Induced Dipoles :',
     &        //,4x,'Atom',14x,'Direct',14x,'PCG',14x,'OPT',i1)
      if (dofull) then
         write (iout,140)
  140    format ()
         do i = 1, n
            if (use(i)) then
               write (iout,150)  i,tdirect(i),tpcg(i),txpt(i)
  150          format (i8,4x,3f18.10)
            end if
         end do
      end if
      write (iout,160)  rdirect,rpcg,rxpt
  160 format (/,5x,'RMS',4x,3f18.10)
c
c     find the RMS of each iteration from the exact dipoles
c
      do k = 0, maxiter
         sum = 0.0d0
         do i = 1, n
            do j = 1, 3
               sum = sum + (ustore(j,i,k)-uexact(j,i))**2
            end do
         end do
         rms(k) = sqrt(sum/dble(npolar))
      end do
c
c     print the RMS between iterations and versus exact dipoles
c
      write (iout,170)
  170 format (/,' Iterative PCG Induced Dipole Convergence :',
     &        //,4x,'Iter',12x,'RMS Change',11x,'RMS vs Exact')
      write (iout,180)  0,rms(0)
  180 format (/,i8,15x,'----',6x,f20.10)
      do k = 1, maxiter
         write (iout,190)  k,drms(k),rms(k)
  190    format (i8,2x,f20.10,3x,f20.10)
         if (rms(k) .lt. 0.5d0*poleps)  goto 200
      end do
  200 continue
c
c     refine the extrapolated OPT coefficients via optimization
c
      maxiter = 10000
      eps = 0.033d0
      delta = 0.0001d0
      write (iout,210)  coptmax
  210 format (/,' Extrapolated OPT',i1,' Coefficient Refinement :',
     &        //,4x,'Iter',8x,'C0',8x,'C1',8x,'C2',8x,'C3',
     &           8x,'C4',6x,'RMS vs Exact',/)
c
c     count number of variables and define the initial simplex
c
      nvar = 0
      do i = 0, coptmax
         if (copt(i) .ne. 0.0d0)  nvar = nvar + 1
      end do
      allocate (var(nvar))
      allocate (yval(nvar+1))
      allocate (p(nvar+1,nvar))
      nvar = 0
      do i = 0, coptmax
         if (copt(i) .ne. 0.0d0) then
            nvar = nvar + 1
            var(nvar) = copt(i)
         end if
      end do
      do i = 1, nvar
         p(1,i) = var(i)
      end do
      yval(1) = optfit (var)
      do k = 1, nvar
         var(k) = var(k) + eps
         do i = 1, nvar
            p(k+1,i) = var(i)
         end do
         yval(k+1) = optfit (var)
         var(k) = var(k) - eps
      end do
c
c     optimize coefficients, then find and print refined values
c
      call simplex (nvar,p,yval,delta,optfit,iter)
      iter = iter + nvar + 1
      rxpt = 1000000.0d0
      do i = 1, nvar+1
         if (yval(i) .lt. rxpt) then
            miny = i
            rxpt = yval(i)
         end if
      end do
      nvar = 0
      do i = 0, coptmax
         if (copt(i) .ne. 0.0d0) then
            nvar = nvar + 1
            copt(i) = p(miny,nvar)
         end if
      end do
      write (iout,220)  iter,(copt(i),i=0,maxopt),rxpt
  220 format (i8,1x,5f10.3,f17.10)
c
c     perform deallocation of some local arrays
c
      deallocate (var)
      deallocate (yval)
      deallocate (p)
      deallocate (rms)
      deallocate (drms)
      deallocate (tdirect)
      deallocate (txpt)
      deallocate (tpcg)
      deallocate (ddirect)
      deallocate (dxpt)
      deallocate (dpcg)
      deallocate (udirect)
      deallocate (upcg)
      deallocate (uxpt)
      deallocate (ustore)
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  function optfit  --  OPT dipole coefficient refinement  ##
c     ##                                                          ##
c     ##############################################################
c
c
      function optfit (var)
      use sizes
      use atoms
      use iounit
      use polar
      use polpot
      use units
      implicit none
      integer i,j
      integer iter,nvar
      real*8 optfit
      real*8 rxpt
      real*8 var(*)
      real*8, allocatable :: uxpt(:,:)
      logical first
      save first,iter
      data first  / .true. /
c
c
c     count the number of times the function has been called
c
      if (first) then
         first = .false.
         iter = -1
      end if
      iter = iter + 1
c
c     copy optimization variables into extrapolation coefficients
c
      nvar = 0
      do i = 0, maxopt
         if (copt(i) .ne. 0.0d0) then
            nvar = nvar + 1
            copt(i) = var(nvar)
         end if
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (uxpt(3,n))
c
c     compute RMS error between OPT and exact SCF dipoles
c
      poltyp = 'OPT'
      call induce
      do i = 1, n
         do j = 1, 3
            uxpt(j,i) = debye * uind(j,i)
         end do
      end do
      rxpt = 0.0d0
      do i = 1, n
         do j = 1, 3
            rxpt = rxpt + (uxpt(j,i)-uexact(j,i))**2
c           rxpt = rxpt + (uxpt(j,i)-uexact(j,i))**6
         end do
      end do
      rxpt = sqrt(rxpt/dble(n))
      if (mod(iter,100) .eq. 0) then
         write (iout,10)  iter,(copt(i),i=0,maxopt),rxpt
   10    format (i8,1x,5f10.3,f17.10)
      end if
c
c     set the return value equal to the RMS error
c
      optfit = rxpt
c
c     perform deallocation of some local arrays
c
      deallocate (uxpt)
      return
      end
