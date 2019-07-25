c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  program saddle  --  find conformational transition state  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "saddle" finds a transition state between two conformational
c     minima using a combination of ideas from the synchronous transit
c     (Halgren-Lipscomb) and quadratic path (Bell-Crighton) methods
c
c
      program saddle
      use sizes
      use atoms
      use iounit
      use keys
      use linmin
      use syntrn
      use titles
      use zcoord
      implicit none
      integer i,its,next
      integer nvar,freeunit
      integer ncalls,niter
      integer ninner,nouter
      integer ncycle,maxcycle
      integer maxinner,maxouter
      real*8 f,g_rms,g_tan,g2
      real*8 saddle1,grdmin
      real*8 reduce,diverge
      real*8 beta,sg,sg0
      real*8 gamma,gammamin
      real*8 x_move,f_move
      real*8 hg,f_old,g2_old
      real*8 f_0,f_1,f_2,f_3
      real*8 p,delta,epsilon
      real*8 angle,rmsvalue
      real*8 energy1,energy2
      real*8, allocatable :: x1(:)
      real*8, allocatable :: y1(:)
      real*8, allocatable :: z1(:)
      real*8, allocatable :: zbond1(:)
      real*8, allocatable :: zang1(:)
      real*8, allocatable :: ztors1(:)
      real*8, allocatable :: x2(:)
      real*8, allocatable :: y2(:)
      real*8, allocatable :: z2(:)
      real*8, allocatable :: zbond2(:)
      real*8, allocatable :: zang2(:)
      real*8, allocatable :: ztors2(:)
      real*8, allocatable :: xx(:)
      real*8, allocatable :: g(:)
      real*8, allocatable :: x_old(:)
      real*8, allocatable :: g_old(:)
      real*8, allocatable :: tan(:)
      real*8, allocatable :: dgdt(:)
      real*8, allocatable :: s0(:)
      real*8, allocatable :: s(:)
      real*8, allocatable :: h0(:)
      logical exist,terminate
      logical scan,spanned
      logical done,newcycle
      character*1 answer
      character*9 status
      character*20 keyword
      character*240 tsfile
      character*240 record
      character*240 string
      external saddle1
c
c
c     set default parameters for the saddle point method
c
      call initial
      terminate = .false.
      ncalls = 0
      nouter = 0
      maxouter = 100
      maxinner = 50
      maxcycle = 4
      epsilon = 0.5d0
      gammamin = 0.00001d0
      diverge = 0.005d0
      reduce = 0.0d0
c
c     set default parameters for the line search
c
      stpmin = 1.0d-16
      stpmax = 2.0d0
      cappa = 0.1d0
      slpmax = 10000.0d0
      angmax = 180.0d0
      intmax = 5
c
c     get coordinates for the first endpoint structure
c
      call getxyz
c
c     perform dynamic allocation of some local arrays
c
      allocate (x1(n))
      allocate (y1(n))
      allocate (z1(n))
      allocate (zbond1(n))
      allocate (zang1(n))
      allocate (ztors1(n))
c
c     store coordinates for the first endpoint structure
c
      do i = 1, n
         x1(i) = x(i)
         y1(i) = y(i)
         z1(i) = z(i)
      end do
c
c     get coordinates for the second endpoint structure
c
      call getxyz
c
c     perform dynamic allocation of some local arrays
c
      allocate (x2(n))
      allocate (y2(n))
      allocate (z2(n))
      allocate (zbond2(n))
      allocate (zang2(n))
      allocate (ztors2(n))
c
c     store coordinates for the second endpoint structure
c
      do i = 1, n
         x2(i) = x(i)
         y2(i) = y(i)
         z2(i) = z(i)
      end do
c
c     setup for the subsequent energy computations
c
      call mechanic
c
c     get any altered values from the keyword file
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:8) .eq. 'DIVERGE ') then
            read (string,*,err=10,end=10)  diverge
         else if (keyword(1:7) .eq. 'REDUCE ') then
            read (string,*,err=10,end=10)  reduce
         else if (keyword(1:9) .eq. 'GAMMAMIN ') then
            read (string,*,err=10,end=10)  gammamin
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
   30    format (/,' Enter RMS Gradient per Atom Criterion [0.1] :  ',$)
         read (input,40)  grdmin
   40    format (f20.0)
      end if
      if (grdmin .le. 0.0d0)  grdmin = 0.1d0
c
c     find out whether syncronous transit scans are desired
c
      scan = .false.
      call nextarg (answer,exist)
      if (.not. exist) then
         write (iout,50)
   50    format (/,' Perform Synchronous Transit Pathway Scans',
     &              ' [N] :  ',$)
         read (input,60)  record
   60    format (a240)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .eq. 'Y')  scan = .true.
c
c     superimpose the two conformational endpoints
c
      call impose (n,x1,y1,z1,n,x2,y2,z2,rmsvalue)
      write (iout,70)  rmsvalue
   70 format (/,' RMS Fit for All Atoms of Both Structures :',f10.4)
c
c     perform dynamic allocation of some global arrays
c
      nvar = 3 * n
      if (.not. allocated(xmin1))  allocate (xmin1(nvar))
      if (.not. allocated(xmin2))  allocate (xmin2(nvar))
      if (.not. allocated(xm))  allocate (xm(nvar))
c
c     copy the superimposed structures into vectors
c
      do i = 1, n
         xmin1(3*i-2) = x1(i)
         xmin1(3*i-1) = y1(i)
         xmin1(3*i) = z1(i)
         xmin2(3*i-2) = x2(i)
         xmin2(3*i-1) = y2(i)
         xmin2(3*i) = z2(i)
      end do
c
c     get and store internal coordinates for first endpoint
c
      do i = 1, n
         x(i) = x1(i)
         y(i) = y1(i)
         z(i) = z1(i)
      end do
      call makeint (0)
      do i = 1, n
         zbond1(i) = zbond(i)
         zang1(i) = zang(i)
         ztors1(i) = ztors(i)
      end do
c
c     get and store internal coordinates for second endpoint
c
      do i = 1, n
         x(i) = x2(i)
         y(i) = y2(i)
         z(i) = z2(i)
      end do
      call makeint (2)
      do i = 1, n
         zbond2(i) = zbond(i)
         zang2(i) = zang(i)
         ztors2(i) = ztors(i)
         if (ztors1(i)-ztors2(i) .gt. 180.0d0) then
            ztors2(i) = ztors2(i) + 360.0d0
         else if (ztors1(i)-ztors2(i) .lt. -180.0d0) then
            ztors1(i) = ztors1(i) + 360.0d0
         end if
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (xx(nvar))
      allocate (g(nvar))
      allocate (x_old(nvar))
      allocate (g_old(nvar))
      allocate (tan(nvar))
      allocate (dgdt(nvar))
      allocate (s0(nvar))
      allocate (s(nvar))
      allocate (h0(nvar))
c
c     get the energies for the two endpoint structures
c
      ncalls = ncalls + 2
      energy1 = saddle1 (xmin1,g)
      energy2 = saddle1 (xmin2,g)
      write (iout,80)  energy1,energy2
   80 format (/,' Energy Value for Endpoint Structure 1 :',f13.4,
     &        /,' Energy Value for Endpoint Structure 2 :',f13.4)
c
c     make a guess at the transition state structure;
c     or use the current guess if one is around
c
      inquire (file='tstate.xyz',exist=exist)
      if (exist) then
         write (iout,90)
   90    format (/,' Using TSTATE.XYZ as the Transition State Estimate')
         its = freeunit ()
         tsfile = 'tstate.xyz'
         call version (tsfile,'old')
         open (unit=its,file=tsfile,status='old')
         rewind (unit=its)
         call readxyz (its)
         close (unit=its)
         do i = 1, n
            xx(3*i-2) = x(i)
            xx(3*i-1) = y(i)
            xx(3*i) = z(i)
         end do
      else
         tpath = 0.5d0
         do i = 1, n
            zbond(i) = (1.0d0-tpath)*zbond1(i) + tpath*zbond2(i)
            zang(i) = (1.0d0-tpath)*zang1(i) + tpath*zang2(i)
            ztors(i) = (1.0d0-tpath)*ztors1(i) + tpath*ztors2(i)
         end do
         call makexyz
         do i = 1, n
            xx(3*i-2) = x(i)
            xx(3*i-1) = y(i)
            xx(3*i) = z(i)
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (x1)
      deallocate (y1)
      deallocate (z1)
      deallocate (zbond1)
      deallocate (zang1)
      deallocate (ztors1)
      deallocate (x2)
      deallocate (y2)
      deallocate (z2)
      deallocate (zbond2)
      deallocate (zang2)
      deallocate (ztors2)
c
c     save the initial estimate of the transition state
c
      do i = 1, n
         x(i) = xx(3*i-2)
         y(i) = xx(3*i-1)
         z(i) = xx(3*i)
      end do
      if (.not. exist) then
         title = 'Transition State Structure'
         ltitle = 26
      end if
      its = freeunit ()
      tsfile = 'tstate.xyz'
      call version (tsfile,'new')
      open (unit=its,file=tsfile,status='new')
      call prtxyz (its)
      close (unit=its)
c
c     start of the major loop for transition state location;
c     first, find the value of the transit path coordinate
c
  100 continue
      nouter = nouter + 1
      call pathval (nvar,xx)
c
c     make a scan along the synchronous transit pathway
c
      if (scan) then
         call pathscan (nvar,xmin1,xmin2,ncalls)
      end if
c
c     set parameters for use in quadratic line maximization
c
      done = .false.
      niter = 1
      ncycle = 1
      delta = 0.01d0
c
c     compute initial point for quadratic line maximization
c
      tpath = ppath
      call pathpnt (nvar,tpath,xx,xmin1,xmin2)
      ncalls = ncalls + 3
      f = saddle1 (xx,g)
      call tangent (nvar,xx,g,g_rms,tan,g_tan,gamma,dgdt)
      write (iout,110)
  110 format (/,' Search for a Maximum along Synchronous Transit :',
     &        /' ST Iter    F Value       Path      RMS G',
     &          '      G Tan      Gamma   FG Call',/)
      write (iout,120)  niter,f,tpath,g_rms,g_tan,gamma,ncalls
  120 format (i6,f13.4,f11.4,f11.4,f11.4,f11.5,i8)
c
c     make an iterative search for quadratic line maximum
c
      do while (.not. done)
         f_0 = f
         tpath = tpath + delta
         call pathpnt (nvar,tpath,xx,xmin1,xmin2)
         ncalls = ncalls + 1
         f_1 = saddle1 (xx,g)
         tpath = tpath - 2.0d0*delta
         call pathpnt (nvar,tpath,xx,xmin1,xmin2)
         tpath = tpath + delta
         ncalls = ncalls + 1
         f_2 = saddle1 (xx,g)
         if (f_1.gt.f_0 .and. f_2.gt.f_0) then
            goto 150
         else if (f_1 .gt. f_0) then
            tpath = tpath + delta
            p = 1.0d0
         else if (f_2 .gt. f_0) then
            tpath = tpath - delta
            p = -1.0d0
            f_1 = f_2
         else
            tpath = tpath + 0.5d0*delta*(f_2-f_1)/(f_1-2.0d0*f_0+f_2)
            goto 130
         end if
         spanned = .false.
         do while (.not. spanned)
            p = 2.0d0 * p
            tpath = tpath + p*delta
            if (tpath .le. 0.0d0) then
               tpath = 0.0d0
               f_2 = energy1
            else if (tpath .ge. 1.0d0) then
               tpath = 1.0d0
               f_2 = energy2
            else
               call pathpnt (nvar,tpath,xx,xmin1,xmin2)
               ncalls = ncalls + 1
               f_2 = saddle1 (xx,g)
            end if
            if (f_2 .gt. f_1) then
               f_0 = f_1
               f_1 = f_2
            else
               spanned = .true.
            end if
         end do
         p = 0.5d0 * p
         tpath = tpath - p*delta
         if (tpath .le. 0.0d0) then
            tpath = 0.0d0
            f_3 = energy1
         else if (tpath .ge. 1.0d0) then
            tpath = 1.0d0
            f_3 = energy2
         else
            call pathpnt (nvar,tpath,xx,xmin1,xmin2)
            ncalls = ncalls + 1
            f_3 = saddle1 (xx,g)
         end if
         if (f_3 .gt. f_1) then
            tpath = tpath + 0.5d0*abs(p)*delta*(f_1-f_2)
     &                         / (f_2-2.0d0*f_3+f_1)
         else
            tpath = tpath - p*delta
            tpath = tpath + 0.5d0*abs(p)*delta*(f_0-f_3)
     &                         / (f_3-2.0d0*f_1+f_0)
         end if
  130    continue
         niter = niter + 1
         call pathpnt (nvar,tpath,xx,xmin1,xmin2)
         ncalls = ncalls + 3
         f = saddle1 (xx,g)
         call tangent (nvar,xx,g,g_rms,tan,g_tan,gamma,dgdt)
         write (iout,140)  niter,f,tpath,g_rms,g_tan,gamma,ncalls
  140    format (i6,f13.4,f11.4,f11.4,f11.4,f11.5,i8)
         if (ncycle.ge.maxcycle .or. gamma.lt.gammamin) then
            done = .true.
         end if
  150    continue
         ncycle = ncycle + 1
         delta = delta * epsilon
      end do
c
c     if the path maximum is too near to an endpoint,
c     then negative curvature has probably been lost
c
      if (tpath.le.0.05d0 .or. tpath.ge.0.95d0) then
         if (.not. scan)  call pathscan (nvar,xmin1,xmin2,ncalls)
         write (iout,160)
  160    format (/,' SADDLE  --  Termination due to Loss',
     &              ' of Negative Curvature')
         call fatal
      end if
c
c     save the current maximum as the transition state estimate
c
      do i = 1, n
         x(i) = xx(3*i-2)
         y(i) = xx(3*i-1)
         z(i) = xx(3*i)
      end do
      its = freeunit ()
      tsfile = 'tstate.xyz'
      call version (tsfile,'old')
      open (unit=its,file=tsfile,status='old')
      rewind (unit=its)
      call prtxyz (its)
      close (unit=its)
c
c     the maximum is located, get ready for minimization
c
      sg = 0.0d0
      do i = 1, nvar
         s0(i) = tan(i)
         sg = sg + s0(i)*dgdt(i)
      end do
      do i = 1, nvar
         h0(i) = dgdt(i) / sg
      end do
c
c     set the initial conjugate direction for minimization
c
      ninner = 0
      g2 = 0.0d0
      hg = 0.0d0
      f_move = 1000000.0d0
      do i = 1, nvar
         g2 = g2 + g(i)**2
         hg = hg + h0(i)*g(i)
      end do
      do i = 1, nvar
         s(i) = -g(i) + hg*s0(i)
      end do
      g_rms = sqrt(g2/dble(n))
      write (iout,170)
  170 format (/,' Search for a Minimum in Conjugate Directions :',
     &        /,' CG Iter    F Value      RMS G     F Move',
     &           '    X Move    Angle   FG Call  Comment',/)
      write (iout,180)  ninner,f,g_rms,ncalls
  180 format (i6,f13.4,f11.4,30x,i7)
c
c     check the termination criterion
c
      if (g_rms .lt. grdmin) then
         terminate = .true.
         write (iout,190)
  190    format (/,' SADDLE  --  Normal Termination at',
     &              ' Transition State')
      end if
c
c     line search to find minimum in conjugate direction
c
      do while (.not. terminate)
         ninner = ninner + 1
         f_old = f
         g2_old = g2
         do i = 1, nvar
            x_old(i) = xx(i)
            g_old(i) = g(i)
         end do
         status = '         '
         angmax = 90.0d0
         call search (nvar,f,g,xx,s,f_move,angle,
     &                  ncalls,saddle1,status)
c
c     if search direction points uphill, use its negative
c
         if (status .eq. 'WideAngle') then
            do i = 1, nvar
               s(i) = -s(i)
            end do
            call search (nvar,f,g,xx,s,f_move,angle,
     &                     ncalls,saddle1,status)
         end if
c
c     compute movement and gradient following line search
c
         f_move = f_old - f
         x_move = 0.0d0
         g2 = 0.0d0
         do i = 1, nvar
            x_move = x_move + (xx(i)-x_old(i))**2
            g2 = g2 + g(i)**2
         end do
         x_move = sqrt(x_move/dble(n))
         g_rms = sqrt(g2/dble(n))
         write (iout,200)  ninner,f,g_rms,f_move,
     &                     x_move,angle,ncalls,status
  200    format (i6,f13.4,f11.4,f11.4,f10.4,f9.2,i7,3x,a9)
c
c     check the termination criteria
c
         if (g_rms .lt. grdmin) then
            terminate = .true.
            write (iout,210)
  210       format (/,' SADDLE  --  Normal Termination at',
     &                 ' Transition State')
         else if (nouter .ge. maxouter) then
            terminate = .true.
            write (iout,220)
  220       format (/,' SADDLE  --  Termination due to Maximum',
     &                 ' Iteration Limit')
         end if
c
c     check to see if another maximization is needed
c
         if (.not. terminate) then
            sg0 = 0.0d0
            do i = 1, nvar
               sg0 = sg0 + s0(i)*g(i)
            end do
            newcycle = .false.
            if (ninner .ge. maxinner)  newcycle = .true.
            if (sg0*sg0/g2 .gt. diverge)  newcycle = .true.
            if (status .ne. ' Success ')  newcycle = .true.
c
c     unfortunately, a new maximization is needed; first save
c     the current minimum as the transition state estimate
c
            if (newcycle) then
               do i = 1, n
                  x(i) = xx(3*i-2)
                  y(i) = xx(3*i-1)
                  z(i) = xx(3*i)
               end do
               its = freeunit ()
               tsfile = 'tstate.xyz'
               call version (tsfile,'old')
               open (unit=its,file=tsfile,status='old')
               rewind (unit=its)
               call prtxyz (its)
               close (unit=its)
c
c     move the path endpoints toward current transition state;
c     then jump to the start of the next maximization cycle
c
               if (reduce .ne. 0.0d0) then
                  call pathval (nvar,xx)
                  tpath = reduce * ppath
                  call pathpnt (nvar,tpath,x_old,xmin1,xmin2)
                  do i = 1, nvar
                     xmin1(i) = x_old(i)
                  end do
                  ncalls = ncalls + 1
                  energy1 = saddle1 (xmin1,g)
                  tpath = 1.0d0 - reduce*(1.0d0-ppath)
                  call pathpnt (nvar,tpath,x_old,xmin1,xmin2)
                  do i = 1, nvar
                     xmin2(i) = x_old(i)
                  end do
                  ncalls = ncalls + 1
                  energy2 = saddle1 (xmin2,g)
               end if
               goto 100
            end if
c
c     find the next conjugate search direction to search;
c     choice of "beta" is Fletcher-Reeves or Polak-Ribiere
c
            hg = 0.0d0
            do i = 1, nvar
               hg = hg + h0(i)*g(i)
            end do
            beta = 0.0d0
            do i = 1, nvar
c              beta = beta + g(i) * g(i)
               beta = beta + g(i) * (g(i)-g_old(i))
            end do
            beta = beta / g2_old
            do i = 1, nvar
               s(i) = -g(i) + hg*s0(i) + beta*s(i)
            end do
         end if
      end do
c
c     write out the final transition state structure
c
      do i = 1, n
         x(i) = xx(3*i-2)
         y(i) = xx(3*i-1)
         z(i) = xx(3*i)
      end do
      its = freeunit ()
      tsfile = 'tstate.xyz'
      call version (tsfile,'old')
      open (unit=its,file=tsfile,status='old')
      rewind (unit=its)
      call prtxyz (its)
      close (unit=its)
c
c     perform deallocation of some local arrays
c
      deallocate (xx)
      deallocate (g)
      deallocate (x_old)
      deallocate (g_old)
      deallocate (tan)
      deallocate (dgdt)
      deallocate (s0)
      deallocate (s)
      deallocate (h0)
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine pathval  --  synchronous transit path values  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "pathval" computes the synchronous transit path value for
c     the specified structure
c
c
      subroutine pathval (nvar,xx)
      use sizes
      use atoms
      use syntrn
      implicit none
      integer i,nvar
      real*8 dr,dp,rmsvalue
      real*8 xx(*)
      real*8, allocatable :: x1(:)
      real*8, allocatable :: y1(:)
      real*8, allocatable :: z1(:)
      real*8, allocatable :: x2(:)
      real*8, allocatable :: y2(:)
      real*8, allocatable :: z2(:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (x1(n))
      allocate (y1(n))
      allocate (z1(n))
      allocate (x2(n))
      allocate (y2(n))
      allocate (z2(n))
c
c     find the value of the transit path coordinate "ppath";
c     it is the ratio of the rms fits to the two endpoints
c
      do i = 1, n
         x1(i) = xmin1(3*i-2)
         y1(i) = xmin1(3*i-1)
         z1(i) = xmin1(3*i)
         x2(i) = xx(3*i-2)
         y2(i) = xx(3*i-1)
         z2(i) = xx(3*i)
      end do
      call impose (n,x1,y1,z1,n,x2,y2,z2,dr)
      do i = 1, n
         x1(i) = xmin2(3*i-2)
         y1(i) = xmin2(3*i-1)
         z1(i) = xmin2(3*i)
         x2(i) = xx(3*i-2)
         y2(i) = xx(3*i-1)
         z2(i) = xx(3*i)
      end do
      call impose (n,x1,y1,z1,n,x2,y2,z2,dp)
      ppath = dr / (dr+dp)
c
c     superimpose on linear transit structure of same path value
c
      do i = 1, n
         x1(i) = (1.0d0-ppath)*xmin1(3*i-2) + ppath*xmin2(3*i-2)
         y1(i) = (1.0d0-ppath)*xmin1(3*i-1) + ppath*xmin2(3*i-1)
         z1(i) = (1.0d0-ppath)*xmin1(3*i) + ppath*xmin2(3*i)
         x2(i) = xx(3*i-2)
         y2(i) = xx(3*i-1)
         z2(i) = xx(3*i)
      end do
      call impose (n,x1,y1,z1,n,x2,y2,z2,rmsvalue)
      do i = 1, n
         xx(3*i-2) = x2(i)
         xx(3*i-1) = y2(i)
         xx(3*i) = z2(i)
      end do
      do i = 1, nvar
         xm(i) = xx(i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (x1)
      deallocate (y1)
      deallocate (z1)
      deallocate (x2)
      deallocate (y2)
      deallocate (z2)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine pathscan  --  scan along the transit pathway  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "pathscan" makes a scan of a synchronous transit pathway by
c     computing structures and energies for specific path values
c
c
      subroutine pathscan (nvar,x0,x1,ncalls)
      use sizes
      use iounit
      use syntrn
      implicit none
      integer i,nvar,ncalls
      real*8 energy,gamma
      real*8 g_rms,g_tan
      real*8 saddle1
      real*8 x0(*)
      real*8 x1(*)
      real*8, allocatable :: xx(:)
      real*8, allocatable :: g(:)
      real*8, allocatable :: tan(:)
      real*8, allocatable :: dgdt(:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (xx(nvar))
      allocate (g(nvar))
      allocate (tan(nvar))
      allocate (dgdt(nvar))
c
c     make a scan along the synchronous transit pathway
c
      write (iout,10)
   10 format (/,' Scan of the Synchronous Transit Pathway :',
     &        /,' N Scan     F Value       Path      RMS G',
     &           '      G Tan      Gamma   FG Call',/)
      do i = 0, 10
         tpath = 0.1d0 * dble(i)
         call pathpnt (nvar,tpath,xx,x0,x1)
         ncalls = ncalls + 3
         energy = saddle1 (xx,g)
         call tangent (nvar,xx,g,g_rms,tan,g_tan,gamma,dgdt)
         write (iout,20)  i,energy,tpath,g_rms,g_tan,gamma,ncalls
   20    format (i6,f13.4,f11.4,f11.4,f11.4,f11.5,i8)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xx)
      deallocate (g)
      deallocate (tan)
      deallocate (dgdt)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine pathpnt  --  get coordinates of path point  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "pathpnt" finds a structure on the synchronous transit path
c     with the specified path value "tpath"
c
c
      subroutine pathpnt (nvar,tpath,xx,x0,x1)
      use sizes
      use inform
      use minima
      implicit none
      integer i,nvar
      real*8 tpath
      real*8 value
      real*8 grdmin
      real*8 transit
      real*8 xx(*)
      real*8 x0(*)
      real*8 x1(*)
      external transit
      external optsave
c
c
c     initialize some parameters for the upcoming optimization
c
      if (debug) then
         iprint = 1
      else
         iprint = 0
      end if
      iwrite = 0
      maxiter = 1000
      grdmin = 0.00001d0
c
c     interpolate coordinates to give initial estimate
c
      do i = 1, nvar
         xx(i) = (1.0d0-tpath)*x0(i) + tpath*x1(i)
      end do
c
c     optimize the synchronous transit function
c
c     call lbfgs (nvar,xx,value,grdmin,transit,optsave)
      call ocvm (nvar,xx,value,grdmin,transit,optsave)
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine tangent  --  synchronous transit tangent  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "tangent" finds the projected gradient on the synchronous
c     transit path for a point along the transit pathway
c
c
      subroutine tangent (nvar,xx,g,g_rms,tan,g_tan,gamma,dgdt)
      use sizes
      use atoms
      use syntrn
      implicit none
      integer i,nvar
      real*8 g_rms,g_tan
      real*8 gamma,delta
      real*8 t0,g2,tan_norm
      real*8 energy,saddle1
      real*8 xx(*)
      real*8 g(*)
      real*8 tan(*)
      real*8 dgdt(*)
      real*8, allocatable :: xf(:)
      real*8, allocatable :: xb(:)
      real*8, allocatable :: gf(:)
      real*8, allocatable :: gb(:)
c
c
c     set the finite difference path increment
c
      delta = 0.01d0
c
c     store the initial pathpnt and compute gradient norm
c
      t0 = tpath
      g2 = 0.0d0
      do i = 1, nvar
         g2 = g2 + g(i)**2
      end do
      g_rms = sqrt(g2/dble(n))
c
c     perform dynamic allocation of some local arrays
c
      allocate (xf(nvar))
      allocate (xb(nvar))
      allocate (gf(nvar))
      allocate (gb(nvar))
c
c     compute the forward difference
c
      do i = 1, nvar
         xf(i) = xx(i)
      end do
      tpath = t0 + delta
      call pathpnt (nvar,tpath,xf,xf,xf)
      energy = saddle1 (xf,gf)
c
c     compute the backward difference
c
      do i = 1, nvar
         xb(i) = xx(i)
      end do
      tpath = t0 - delta
      call pathpnt (nvar,tpath,xb,xb,xb)
      energy = saddle1 (xb,gb)
      tpath = t0
c
c     compute tangent to the path, and projected gradient
c
      tan_norm = 0.0d0
      do i = 1, nvar
         tan(i) = xf(i) - xb(i)
         tan_norm = tan_norm + tan(i)**2
         dgdt(i) = gf(i) - gb(i)
      end do
      tan_norm = sqrt(tan_norm)
      g_tan = 0.0d0
      do i = 1, nvar
         tan(i) = tan(i) / tan_norm
         g_tan = g_tan + g(i)*tan(i)
      end do
      g_tan = g_tan / sqrt(dble(n))
      gamma = (g_tan/g_rms)**2
c
c     perform deallocation of some local arrays
c
      deallocate (xf)
      deallocate (xb)
      deallocate (gf)
      deallocate (gb)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  function transit  --  synchronous transit evaluation  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "transit" evaluates the synchronous transit function and
c     gradient; linear and quadratic transit paths are available
c
c
      function transit (xx,g)
      use sizes
      use atoms
      use syntrn
      implicit none
      integer i,j,nvar
      integer ix,iy,iz
      integer jx,jy,jz
      real*8 transit,value
      real*8 xci,yci,zci
      real*8 xcd,ycd,zcd
      real*8 x1i,y1i,z1i
      real*8 x1d,y1d,z1d
      real*8 x2i,y2i,z2i
      real*8 x2d,y2d,z2d
      real*8 xmi,ymi,zmi
      real*8 xmd,ymd,zmd
      real*8 gamma,term
      real*8 termx,termy,termz
      real*8 cutoff,cutoff2
      real*8 r1,r2,rc,rm
      real*8 ri,ri4,rd
      real*8 wi,wc,wd
      real*8 tq,pq
      real*8 xx(*)
      real*8 g(*)
      character*9 mode
c
c
c     zero out the synchronous transit function and gradient
c
      value = 0.0d0
      nvar = 3 * n
      do i = 1, nvar
         g(i) = 0.0d0
      end do
      tq = 1.0d0 - tpath
c
c     set the cutoff distance for interatomic distances
c
      cutoff = 1000.0d0
      cutoff2 = cutoff**2
c
c     set the type of synchronous transit path to be used
c
      if (ppath .eq. 0.0d0) then
         mode = 'LINEAR'
      else
         mode = 'QUADRATIC'
         pq = 1.0d0 - ppath
      end if
c
c     portion based on interpolated interatomic distances
c
      do i = 1, n-1
         iz = 3 * i
         iy = iz - 1
         ix = iz - 2
         xci = xx(ix)
         yci = xx(iy)
         zci = xx(iz)
         x1i = xmin1(ix)
         y1i = xmin1(iy)
         z1i = xmin1(iz)
         x2i = xmin2(ix)
         y2i = xmin2(iy)
         z2i = xmin2(iz)
         if (mode .eq. 'QUADRATIC') then
            xmi = xm(ix)
            ymi = xm(iy)
            zmi = xm(iz)
         end if
         do j = i+1, n
            jz = 3 * j
            jy = jz - 1
            jx = jz - 2
            xcd = xci - xx(jx)
            ycd = yci - xx(jy)
            zcd = zci - xx(jz)
            x1d = x1i - xmin1(jx)
            y1d = y1i - xmin1(jy)
            z1d = z1i - xmin1(jz)
            x2d = x2i - xmin2(jx)
            y2d = y2i - xmin2(jy)
            z2d = z2i - xmin2(jz)
            rc = xcd**2 + ycd**2 + zcd**2
            r1 = x1d**2 + y1d**2 + z1d**2
            r2 = x2d**2 + y2d**2 + z2d**2
            if (min(rc,r1,r2) .lt. cutoff2) then
               rc = sqrt(rc)
               r1 = sqrt(r1)
               r2 = sqrt(r2)
               ri = tq*r1 + tpath*r2
               if (mode .eq. 'QUADRATIC') then
                  xmd = xmi - xm(jx)
                  ymd = ymi - xm(jy)
                  zmd = zmi - xm(jz)
                  rm = sqrt(xmd**2+ymd**2+zmd**2)
                  gamma = (rm-pq*r1-ppath*r2) / (ppath*pq)
                  ri = ri + gamma*tpath*tq
               end if
               ri4 = ri**4
               rd = rc - ri
               value = value + rd**2/ri4
               term = 2.0d0 * rd/(ri4*rc)
               termx = term * xcd
               termy = term * ycd
               termz = term * zcd
               g(ix) = g(ix) + termx
               g(iy) = g(iy) + termy
               g(iz) = g(iz) + termz
               g(jx) = g(jx) - termx
               g(jy) = g(jy) - termy
               g(jz) = g(jz) - termz
            end if
         end do
      end do
c
c     portion used to supress rigid rotations and translations
c
      do i = 1, nvar
         wc = xx(i)
         wi = tq*xmin1(i) + tpath*xmin2(i)
         wd = wc - wi
         value = value + 0.000001d0*wd**2
         g(i) = g(i) + 0.000002d0*wd
      end do
      transit = value
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  function saddle1  --  energy and gradient for saddle  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "saddle1" is a service routine that computes the energy and
c     gradient for transition state optimization
c
c
      function saddle1 (xx,g)
      use sizes
      use atoms
      implicit none
      integer i
      real*8 e,saddle1
      real*8 xx(*)
      real*8 g(*)
c
c
c     copy optimization values to coordinates and find gradient
c
      do i = 1, n
         x(i) = xx(3*i-2)
         y(i) = xx(3*i-1)
         z(i) = xx(3*i)
      end do
      call gradient (e,g)
      saddle1 = e
      return
      end
