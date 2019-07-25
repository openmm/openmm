c
c
c     ################################################################
c     ##                   COPYRIGHT (C) 2001 by                    ##
c     ##  Michael Schnieders, Alan Grossfield & Jay William Ponder  ##
c     ##                    All Rights Reserved                     ##
c     ################################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  program monte  --  Monte Carlo-Minimization search method  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "monte" performs a Monte Carlo-Minimization conformational
c     search using Cartesian single atom or torsional move sets
c
c     literature references:
c
c     Z. Li and H. A. Scheraga, "Monte Carlo-Minimization Approach
c     to the Multiple-Minima Problem in Protein Folding", Proc. Natl.
c     Acad. Sci. USA, 84, 6611-6615 (1987)
c
c     D. J. Wales, "Energy Landscapes with Applications to Clusters,
c     Biomolecules and Glasses", Cambridge University Press, 2003,
c     Section 6.7.4
c
c
      program monte
      use sizes
      use atoms
      use files
      use inform
      use iounit
      use omega
      use output
      use units
      use usage
      use zcoord
      implicit none
      integer i,k,m,next
      integer keep,nbig
      integer nmap,lext
      integer istep,nstep
      integer ixyz,nusave
      integer freeunit
      real*8 global,ratio
      real*8 big,eps,size
      real*8 grdmin,temper
      real*8 minimum,pminimum
      real*8 tsize,factor
      real*8 beta,boltz
      real*8 random,trial
      real*8 vector(3)
      real*8, allocatable :: xg(:)
      real*8, allocatable :: yg(:)
      real*8, allocatable :: zg(:)
      real*8, allocatable :: xi(:)
      real*8, allocatable :: yi(:)
      real*8, allocatable :: zi(:)
      real*8, allocatable :: xp(:)
      real*8, allocatable :: yp(:)
      real*8, allocatable :: zp(:)
      logical exist
      logical reset
      logical torsmove
      logical, allocatable :: usave(:)
      character*1 answer
      character*6 status
      character*7 ext
      character*240 xyzfile
      character*240 record
      character*240 string
      external random
c
c
c     set up the structure and mechanics calculation
c
      call initial
      call getxyz
      call mechanic
c
c     initialize values of some counters and parameters
c
      keep = 0
      nbig = 0
      nmap = 0
      eps = 0.0001d0
      big = 100000.0d0
      reset = .false.
c
c     get the desired number of Monte Carlo steps
c
      nstep = -1
      call nextarg (string,exist)
      if (exist)  read (string,*,err=10,end=10)  nstep
   10 continue
      if (nstep .le. 0) then
         write (iout,20)
   20    format (/,' Number of Monte Carlo Steps [1000] :  ', $)
         read (input,30)  nstep
   30    format (i15)
         if (nstep .le. 0)  nstep = 1000
      end if
c
c     choose either the torsional or single atom move set
c
      torsmove = .false.
      call nextarg (answer, exist)
      if (.not. exist) then
         write (iout,40)
   40    format (/,' Use [C]artesian or [T]orsional Moves [C] :  ',$)
         read (input,50)  record
   50    format (a240)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .eq. 'T')  torsmove = .true.
c
c     perform dynamic allocation of some local arrays
c
      allocate (usave(n))
c
c     generate the internal coordinates, keep active atom list
c
      if (torsmove) then
         call makeint (0)
         call initrot
         nusave = nuse
         nuse = n
         do i = 1, n
            usave(i) = use(i)
            use(i) = .true.
         end do
      end if
c
c     get the desired Cartesian or torsional step size
c
      size = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=60,end=60)  size
   60 continue
      if (size .lt. 0.0d0) then
         if (torsmove) then
            write (iout,70)
   70       format (/,' Enter Maximum Step in Degrees [180] :  ', $)
         else
            write (iout,80)
   80       format (/,' Enter Maximum Step in Angstroms [3.0] :  ', $)
         end if
         read (input,90)  string
   90    format (a240)
         read (string,*,err=100,end=100)  size
  100    continue
         if (size .lt. 0.0d0) then
            if (torsmove) then
               size = 180.0d0
            else
               size = 3.0d0
            end if
         end if
         if (torsmove)  size = min(size,180.0d0)
      end if
c
c     get the desired temperature for Metropolis criterion
c
      temper = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=110,end=110)  temper
  110 continue
      if (temper .lt. 0.0d0) then
         write (iout,120)
  120    format (/,' Enter the Desired Temperature in Degrees',
     &              ' K [500] :  ', $)
         read (input,130)  string
  130    format (a240)
         read (string,*,err=140,end=140)  temper
  140    continue
         if (temper .lt. 0.0d0)  temper = 500.0d0
      end if
      beta = 1.0d0 / (gasconst*temper)
c
c     get the gradient convergence for local minimizations
c
      grdmin = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=150,end=150)  grdmin
  150 continue
      if (grdmin .lt. 0.0d0) then
         write (iout,160)
  160    format (/,' Enter RMS Gradient Criterion [0.01] :  ', $)
         read (input,170)  string
  170    format (a240)
         read (string,*,err=180,end=180)  grdmin
  180    continue
         if (grdmin .lt. 0.0d0)  grdmin = 0.01
      end if
c
c     print some information prior to initial iteration
c
      write (iout,190)
  190 format (/,' Monte Carlo Minimization Global Search :')
      write (iout,200)
  200 format (/,' MCM Iter       Current         Global       Temper',
     &           '      Ratio      Status',/)
      flush (iout)
c
c     perform dynamic allocation of some local arrays
c
      allocate (xg(n))
      allocate (yg(n))
      allocate (zg(n))
      allocate (xi(n))
      allocate (yi(n))
      allocate (zi(n))
      allocate (xp(n))
      allocate (yp(n))
      allocate (zp(n))
c
c     create and open an output file if using archive mode
c
      if (archive) then
         ixyz = freeunit ()
         xyzfile = filename(1:leng)
         call suffix (xyzfile,'arc','new')
         open (unit=ixyz,file=xyzfile,status='new')
         close (unit=ixyz)
      end if
c
c     store the coordinates, then perform a minimization
c
      do i = 1, n
         xi(i) = x(i)
         yi(i) = y(i)
         zi(i) = z(i)
      end do
      if (torsmove) then
         nuse = nusave
         do i = 1, n
            use(i) = usave(i)
         end do
      end if
      call mcmstep (minimum,grdmin)
      if (torsmove) then
         nuse = n
         do i = 1, n
            use(i) = .true.
         end do
      end if
      pminimum = minimum
      write (iout,210)  0,minimum
  210 format (i8,3x,f12.4)
c
c     save coordinates as the initial global minimum
c
      do i = 1, n
         xg(i) = x(i)
         yg(i) = y(i)
         zg(i) = z(i)
      end do
      global = minimum
      nmap = nmap + 1
      lext = 3
      call numeral (nmap,ext,lext)
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
      write (iout,220)  nmap,global
  220 format (/,4x,'Minimum Energy Structure',i7,6x,f16.4,/)
      call flush (iout)
c
c     optionally reset coordinates to before the minimization
c
      if (reset) then
         do i = 1, n
            x(i) = xi(i)
            y(i) = yi(i)
            z(i) = zi(i)
         end do
      end if
      if (torsmove)  call makeint (2)
c
c     store the prior coordinates to start each MCM iteration
c
      do istep = 1, nstep
         do i = 1, n
            xp(i) = x(i)
            yp(i) = y(i)
            zp(i) = z(i)
         end do
c
c     generate random angle moves for a few torsions
c
         if (torsmove) then
            m = int(-log(max(random(),0.0001d0))) + 1
            do i = 1, m
               k = int(nomega * random()) + 1
               k = zline(k)
               tsize = 2.0d0 * size * (random()-0.5d0)
               ztors(k) = ztors(k) + tsize
               if (ztors(k) .gt. 180.0d0) then
                  ztors(k) = ztors(k) - 360.0d0
               else if (ztors(k) .lt. -180.0d0) then
                  ztors(k) = ztors(k) + 360.0d0
               end if
            end do
            call makexyz
c
c     generate a random Cartesian move for each atom
c
         else
            do i = 1, n
               if (use(i)) then
                  call ranvec (vector)
                  factor = size * random ()
                  x(i) = x(i) + factor*vector(1)
                  y(i) = y(i) + factor*vector(2)
                  z(i) = z(i) + factor*vector(3)
               end if
            end do
         end if
c
c     store the coordinates, then perform a minimization
c
         do i = 1, n
            xi(i) = x(i)
            yi(i) = y(i)
            zi(i) = z(i)
         end do
         if (torsmove) then
            nuse = nusave
            do i = 1, n
               use(i) = usave(i)
            end do
         end if
         call mcmstep (minimum,grdmin)
         if (torsmove) then
            nuse = n
            do i = 1, n
               use(i) = .true.
            end do
         end if
c
c     test for an unreasonably low energy at the minimum
c
         if (minimum .lt. -big)  minimum = big
c
c     step is probably degenerate if energy is identical
c
         if (abs(minimum-pminimum) .le. eps) then
            status = 'Same'
            pminimum = minimum
c
c     accept the step if the new minimum has lower energy
c
         else if (minimum .le. pminimum) then
            status = 'Accept'
            pminimum = minimum
c
c     if the energy increased, apply the Metropolis criterion
c
         else
            boltz = exp(-beta*(minimum-pminimum))
            trial = random ()
c
c     reject the step if the energy increase is too large
c
            if (boltz .lt. trial) then
               status = 'Reject'
c
c     accept the step if the energy increase is small enough
c
            else
               status = 'Accept'
               pminimum = minimum
            end if
         end if
c
c     save coordinates with the best energy as global minimum
c
         if (minimum .lt. global-eps) then
            do i = 1, n
               xg(i) = x(i)
               yg(i) = y(i)
               zg(i) = z(i)
            end do
            global = minimum
            nmap = nmap + 1
            lext = 3
            call numeral (nmap,ext,lext)
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
            write (iout,230)  nmap,global
  230       format (/,4x,'Minimum Energy Structure',i7,6x,f16.4,/)
            flush (iout)
         end if
c
c     update the overall Monte Carlo acceptance ratio
c
         if (status .eq. 'Accept')  keep = keep + 1
         ratio = dble(keep) / dble(istep)
c
c     print intermediate results for the current iteration
c
         if (minimum .lt. big) then
            nbig = 0
            write (iout,240)  istep,minimum,global,temper,ratio,status
  240       format (i8,3x,f12.4,3x,f12.4,3x,f10.2,3x,f8.3,6x,a6)
         else
            nbig = nbig + 1
            write (iout,250)  istep,global,temper,ratio,status
  250       format (i8,9x,'------',3x,f12.4,3x,f10.2,3x,f8.3,6x,a6)
         end if
         flush (iout)
c
c     restore global minimum after repeated bad iterations
c
         if (nbig .ge. 3) then
            nbig = 0
            do i = 1, n
               x(i) = xg(i)
               y(i) = yg(i)
               z(i) = zg(i)
            end do
c
c     optionally reset coordinates to before the minimization
c
         else if (status.eq.'Same' .or. status.eq.'Accept') then
            if (reset) then
               do i = 1, n
                  x(i) = xi(i)
                  y(i) = yi(i)
                  z(i) = zi(i)
               end do
            end if
c
c     restore coordinates to those from the previous iteration
c
         else if (status .eq. 'Reject') then
            do i = 1, n
               x(i) = xp(i)
               y(i) = yp(i)
               z(i) = zp(i)
            end do
         end if
c
c     update internal coordinates if using torsional moves
c
         if (torsmove)  call makeint (2)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (usave)
      deallocate (xg)
      deallocate (yg)
      deallocate (zg)
      deallocate (xi)
      deallocate (yi)
      deallocate (zi)
      deallocate (xp)
      deallocate (yp)
      deallocate (zp)
c
c     write out the final global minimum energy value
c
      if (digits .ge. 8) then
         write (iout,260)  global
  260    format (/,' Global Minimum Energy Value :',1x,f20.8)
      else if (digits .ge. 6) then
         write (iout,270)  global
  270    format (/,' Global Minimum Energy Value :',3x,f18.6)
      else
         write (iout,280)  global
  280    format (/,' Global Minimum Energy Value :',5x,f16.4)
      end if
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  function mcmstep  --  minimization phase of an MCM step  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "mcmstep" implements the minimization phase of an MCM step
c     via Cartesian minimization following a Monte Carlo step
c
c
      subroutine mcmstep (minimum,grdmin)
      use sizes
      use atoms
      use files
      use inform
      use output
      use usage
      implicit none
      integer i,nvar
      real*8 mcm1,minimum,grdmin
      real*8, allocatable :: xx(:)
      character*6 mode,method
      external mcm1,mcm2,optsave
c
c
c     prepare for the truncated Newton minimization
c
      mode = 'AUTO'
      method = 'AUTO'
      verbose = .false.
      iprint = 0
      iwrite = 0
      coordtype = 'CARTESIAN'
c
c     perform dynamic allocation of some local arrays
c
      allocate (xx(3*n))
c
c     translate the coordinates of each active atom
c
      nvar = 0
      do i = 1, n
         if (use(i)) then
            nvar = nvar + 1
            xx(nvar) = x(i)
            nvar = nvar + 1
            xx(nvar) = y(i)
            nvar = nvar + 1
            xx(nvar) = z(i)
         end if
      end do
c
c     make the call to the optimization routine
c
      call tncg (mode,method,nvar,xx,minimum,grdmin,
     &                  mcm1,mcm2,optsave)
c
c     untranslate the final coordinates for active atoms
c
      nvar = 0
      do i = 1, n
         if (use(i)) then
            nvar = nvar + 1
            x(i) = xx(nvar)
            nvar = nvar + 1
            y(i) = xx(nvar)
            nvar = nvar + 1
            z(i) = xx(nvar)
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xx)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  function mcm1  --  energy and gradient for MCM search  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "mcm1" is a service routine that computes the energy and
c     gradient for truncated Newton optimization in Cartesian
c     coordinate space
c
c
      function mcm1 (xx,g)
      use sizes
      use atoms
      use usage
      implicit none
      integer i,nvar
      real*8 mcm1,e
      real*8 xx(*)
      real*8 g(*)
      real*8, allocatable :: derivs(:,:)
c
c
c     translate optimization parameters to atomic coordinates
c
      nvar = 0
      do i = 1, n
         if (use(i)) then
            nvar = nvar + 1
            x(i) = xx(nvar)
            nvar = nvar + 1
            y(i) = xx(nvar)
            nvar = nvar + 1
            z(i) = xx(nvar)
         end if
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (derivs(3,n))
c
c     compute and store the energy and gradient
c
      call gradient (e,derivs)
      mcm1 = e
c
c     store Cartesian gradient as optimization gradient
c
      nvar = 0
      do i = 1, n
         if (use(i)) then
            nvar = nvar + 1
            g(nvar) = derivs(1,i)
            nvar = nvar + 1
            g(nvar) = derivs(2,i)
            nvar = nvar + 1
            g(nvar) = derivs(3,i)
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (derivs)
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine mcm2  --  Hessian values for MCM search  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "mcm2" is a service routine that computes the sparse matrix
c     Hessian elements for truncated Newton optimization in Cartesian
c     coordinate space
c
c
      subroutine mcm2 (mode,xx,h,hinit,hstop,hindex,hdiag)
      use sizes
      use atoms
      use usage
      implicit none
      integer i,j,k,nvar
      integer hinit(*)
      integer hstop(*)
      integer hindex(*)
      integer, allocatable :: hvar(:)
      integer, allocatable :: huse(:)
      real*8 xx(*)
      real*8 hdiag(*)
      real*8 h(*)
      character*4 mode
c
c
c     translate optimization parameters to atomic coordinates
c
      if (mode .eq. 'NONE')  return
      nvar = 0
      do i = 1, n
         if (use(i)) then
            nvar = nvar + 1
            x(i) = xx(nvar)
            nvar = nvar + 1
            y(i) = xx(nvar)
            nvar = nvar + 1
            z(i) = xx(nvar)
         end if
      end do
c
c     compute and store the Hessian elements
c
      call hessian (h,hinit,hstop,hindex,hdiag)
c
c     perform dynamic allocation of some local arrays
c
      allocate (hvar(nvar))
      allocate (huse(3*n))
c
c     transform the sparse Hessian to use only active atoms
c
      nvar = 0
      if (nuse .ne. n) then
         do i = 1, n
            k = 3 * (i-1)
            if (use(i)) then
               do j = 1, 3
                  nvar = nvar + 1
                  hvar(nvar) = j + k
                  huse(j+k) = nvar
               end do
            else
               do j = 1, 3
                  huse(j+k) = 0
               end do
            end if
         end do
         do i = 1, nvar
            k = hvar(i)
            hinit(i) = hinit(k)
            hstop(i) = hstop(k)
            hdiag(i) = hdiag(k)
            do j = hinit(i), hstop(i)
               hindex(j) = huse(hindex(j))
            end do
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (hvar)
      deallocate (huse)
      return
      end
