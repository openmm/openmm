c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine diffeq  --  differential equation integration  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "diffeq" performs the numerical integration of an ordinary
c     differential equation using an adaptive stepsize method to
c     solve the corresponding coupled first-order equations of the
c     general form dyi/dx = f(x,y1,...,yn) for yi = y1,...,yn
c
c     variables and parameters :
c
c     nvar      number of coupled first-order differential equations
c     y         contains the values of the dependent variables
c     x1        value of the beginning integration limit
c     x2        value of the ending integration limit
c     eps       relative accuracy required of the integration steps
c     h1        initial guess for the first integration stepsize
c     hmin      minimum allowed integration stepsize
c     nok       number of initially successful integration steps
c     nbad      number of integration steps that required retry
c
c     required external routines :
c
c     gvalue    subroutine to find the right-hand side of the
c                  first-order differential equations
c
c
      subroutine diffeq (nvar,y,x1,x2,eps,h1,hmin,nok,nbad,gvalue)
      use iounit
      implicit none
      real*8 tiny
      parameter (tiny=1.0d-30)
      integer i,nvar,nok,nbad
      integer nstep,maxstep
      real*8 x,x1,x2,h,h1
      real*8 eps,hnext
      real*8 hmin,hdid
      real*8 y(*)
      real*8, allocatable :: dydx(:)
      real*8, allocatable :: yscal(:)
      logical terminate
      character*7 status
      external gvalue
c
c
c     initialize starting limit, step size and status counters
c
      terminate = .false.
      x = x1
      h = sign(h1,x2-x1)
      nstep = 0
      nok = 0
      nbad = 0
      maxstep = 1000
c
c     perform dynamic allocation of some local arrays
c
      allocate (dydx(nvar))
      allocate (yscal(nvar))
c
c     perform a series of individual integration steps
c
      do while (.not. terminate)
         call gvalue (x,y,dydx)
         do i = 1, nvar
            yscal(i) = abs(y(i)) + abs(h*dydx(i)) + tiny
         end do
c
c     set the final step to stop at the integration limit
c
         if ((x+h-x2)*(x+h-x1) .gt. 0.0d0)  h = x2 - x
c
c     take a Bulirsch-Stoer integration step
c
         call bsstep (nvar,x,dydx,y,h,eps,yscal,hdid,hnext,gvalue)
c
c     mark the current step as either good or bad
c
         if (hdid .eq. h) then
            nok = nok + 1
            status = 'Success'
         else
            nbad = nbad + 1
            status = ' Retry '
         end if
c
c     update stepsize and get information about the current step
c
         h = hnext
         nstep = nstep + 1
         call gdastat (nstep,x,y,status)
c
c     test for convergence to the final integration limit
c
         if ((x-x2)*(x2-x1) .ge. 0.0d0) then
            write (iout,10)
   10       format (/,' DIFFEQ  --  Normal Termination',
     &                 ' at Integration Limit')
            terminate = .true.
         end if
c
c     test for a trial stepsize that is too small
c
         if (abs(hnext) .lt. hmin) then
            write (iout,20)
   20       format (/,' DIFFEQ  --  Incomplete Integration',
     &                 ' due to SmallStep')
            terminate = .true.
         end if
c
c     test for too many total integration steps
c
         if (nstep .ge. maxstep) then
            write (iout,30)
   30       format (/,' DIFFEQ  --  Incomplete Integration',
     &                 ' due to IterLimit')
            terminate = .true.
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (dydx)
      deallocate (yscal)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine bsstep  --  Bulirsch-Stoer integration step  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "bsstep" takes a single Bulirsch-Stoer step with monitoring
c     of local truncation error to ensure accuracy
c
c     literature reference:
c
c     W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P.
c     Flannery, Numerical Recipes (Fortran), 2nd Ed., Cambridge
c     University Press, 1992, Section 16.4
c
c
      subroutine bsstep (nvar,x,dydx,y,htry,eps,yscal,hdid,hnext,gvalue)
      use iounit
      implicit none
      integer kmaxx,imax
      real*8 safe1,safe2
      real*8 redmax,redmin
      real*8 tiny,scalmx
      parameter (kmaxx=8)
      parameter (imax=kmaxx+1)
      parameter (safe1=0.25d0)
      parameter (safe2=0.7d0)
      parameter (redmax=1.0d-5)
      parameter (redmin=0.7d0)
      parameter (tiny=1.0d-30)
      parameter (scalmx=0.1d0)
      integer i,iq,k,kk,nvar
      integer km,kmax,kopt
      integer nseq(imax)
      real*8 eps,eps1,epsold
      real*8 h,hdid,hnext,htry
      real*8 errmax,fact,red
      real*8 scale,work,wrkmin
      real*8 x,xest,xnew
      real*8 dydx(*)
      real*8 y(*)
      real*8 yscal(*)
      real*8 a(imax)
      real*8 err(kmaxx)
      real*8 alf(kmaxx,kmaxx)
      real*8, allocatable :: yerr(:)
      real*8, allocatable :: ysav(:)
      real*8, allocatable :: yseq(:)
      logical first,reduct
      save a,alf,epsold,first
      save kmax,kopt,nseq,xnew
      external gvalue
      data first  / .true. /
      data epsold / -1.0d0 /
      data nseq   / 2,4,6,8,10,12,14,16,18 /
c
c
c     setup prior to the Bulirsch-Stoer integration step
c
      if (eps .ne. epsold) then
         hnext = -1.0d29
         xnew = -1.0d29
         eps1 = safe1 * eps
         a(1) = 1.0d0 + dble(nseq(1))
         do k = 1, kmaxx
            a(k+1) = a(k) + dble(nseq(k+1))
         end do
         do iq = 2, kmaxx
            do k = 1, iq-1
               alf(k,iq) = eps1**((a(k+1)-a(iq+1))/((a(iq+1)-a(1)+1.0d0)
     &                                                 *(2*k+1)))
            end do
         end do
         epsold = eps
         do kopt = 2, kmaxx-1
            kmax = kopt
            if (a(kopt+1) .gt. a(kopt)*alf(kopt-1,kopt))  goto 10
         end do
   10    continue
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (yerr(nvar))
      allocate (ysav(nvar))
      allocate (yseq(nvar))
c
c     make an integration step using Bulirsch-Stoer method
c
      h = htry
      do i = 1, nvar
         ysav(i) = y(i)
      end do
      if (h.ne.hnext .or. x.ne.xnew) then
         first = .true.
         kopt = kmax
      end if
      reduct = .false.
   20 continue
      do k = 1, kmax
         xnew = x + h
         if (xnew .eq. x) then
            write (iout,30)
   30       format (' BSSTEP  --  Underflow of Step Size')
            call fatal
         end if
         call mmid (nseq(k),h,nvar,x,dydx,ysav,yseq,gvalue)
         xest = (h/dble(nseq(k)))**2
         call pzextr (k,nvar,xest,yseq,y,yerr)
         if (k .ne. 1) then
            errmax = tiny
            do i = 1, nvar
               errmax = max(errmax,abs(yerr(i)/yscal(i)))
            end do
            errmax = errmax / eps
            km = k - 1
            err(km) = (errmax/safe1)**(1.0d0/(2*km+1))
         end if
         if (k.ne.1 .and. (k.ge.kopt-1 .or. first)) then
            if (errmax .lt. 1.0d0)  goto 50
            if (k.eq.kmax .or. k.eq.kopt+1) then
               red = safe2 / err(km)
               goto 40
            else if (k .eq. kopt) then
               if (alf(kopt-1,kopt) .lt. err(km)) then
                  red = 1.0d0 / err(km)
                  goto 40
               end if
            else if (kopt .eq. kmax) then
               if (alf(km,kmax-1) .lt. err(km)) then
                  red = alf(km,kmax-1) * safe2 / err(km)
                  goto 40
               end if
            else if (alf(km,kopt) .lt. err(km)) then
               red = alf(km,kopt-1) / err(km)
               goto 40
            end if
         end if
      end do
   40 continue
      red = min(red,redmin)
      red = max(red,redmax)
      h = h * red
      reduct = .true.
      goto 20
   50 continue
      x = xnew
      hdid = h
      first = .false.
      wrkmin = 1.0d35
      do kk = 1, km
         fact = max(err(kk),scalmx)
         work = fact * a(kk+1)
         if (work .lt. wrkmin) then
            scale = fact
            wrkmin = work
            kopt = kk + 1
         end if
      end do
      hnext = h / scale
      if (kopt.ge.k .and. kopt.ne.kmax .and. .not.reduct) then
         fact = max(scale/alf(kopt-1,kopt),scalmx)
         if (a(kopt+1)*fact .le. wrkmin) then
            hnext = h / fact
            kopt = kopt + 1
         end if
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (yerr)
      deallocate (ysav)
      deallocate (yseq)
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine mmid  --  takes a modified midpoint step  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "mmid" implements a modified midpoint method to advance the
c     integration of a set of first order differential equations
c
c
      subroutine mmid (nstep,htot,nvar,xs,dydx,y,yout,gvalue)
      implicit none
      integer i,k
      integer nstep,nvar
      real*8 htot,h,h2
      real*8 xs,x,temp
      real*8 dydx(*)
      real*8 y(*)
      real*8 yout(*)
      real*8, allocatable :: ym(:)
      real*8, allocatable :: yn(:)
      external gvalue
c
c
c     set substep size based on number of steps to be taken
c
      h = htot / dble(nstep)
      h2 = 2.0d0 * h
c
c     perform dynamic allocation of some local arrays
c
      allocate (ym(nvar))
      allocate (yn(nvar))
c
c     take the first substep and get values at ends of step
c
      do i = 1, nvar
         ym(i) = y(i)
         yn(i) = y(i) + h*dydx(i)
      end do
      x = xs + h
      call gvalue (x,yn,yout)
c
c     take the second and subsequent substeps
c
      do k = 2, nstep
         do i = 1, nvar
            temp = ym(i) + h2*yout(i)
            ym(i) = yn(i)
            yn(i) = temp
         end do
         x = x + h
         call gvalue (x,yn,yout)
      end do
c
c     complete the update of values for the last substep
c
      do i = 1, nvar
         yout(i) = 0.5d0 * (ym(i)+yn(i)+h*yout(i))
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (ym)
      deallocate (yn)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine pzextr  --  polynomial extrapolation method  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "pzextr" is a polynomial extrapolation routine used during
c     Bulirsch-Stoer integration of ordinary differential equations
c
c
      subroutine pzextr (iest,nvar,xest,yest,yz,dy)
      use sizes
      implicit none
      integer maxgda,imax
      parameter (maxgda=4*maxatm)
      parameter (imax=13)
      integer i,j,iest,nvar
      real*8 xest,delta
      real*8 f1,f2,q
      real*8 x(imax)
      real*8 yz(*)
      real*8 dy(*)
      real*8 yest(*)
      real*8, allocatable :: d(:)
      real*8, allocatable, save :: qcol(:,:)
      save x
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (d(nvar))
      if (.not. allocated(qcol))  allocate (qcol(nvar,imax))
c
c     polynomial extrapolation needed for Bulirsch-Stoer step
c
      x(iest) = xest
      do j = 1, nvar
         dy(j) = yest(j)
         yz(j) = yest(j)
      end do
      if (iest .eq. 1) then
         do j = 1, nvar
            qcol(j,1) = yest(j)
         end do
      else
         do j = 1, nvar
            d(j) = yest(j)
         end do
         do i = 1, iest-1
            delta = 1.0d0 / (x(iest-i)-xest)
            f1 = xest * delta
            f2 = x(iest-i) * delta
            do j = 1, nvar
               q = qcol(j,i)
               qcol(j,i) = dy(j)
               delta = d(j) - q
               dy(j) = f1 * delta
               d(j) = f2 * delta
               yz(j) = yz(j) + dy(j)
            end do
         end do
         do j = 1, nvar
            qcol(j,iest) = dy(j)
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (d)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine gdastat  --  results for GDA integration step  ##
c     ##                                                            ##
c     ################################################################
c
c     "gdastat" finds the energy, radius of gyration, and average M2
c     for a GDA integration step; also saves the coordinates
c
c
      subroutine gdastat (nstep,beta,xx,status)
      use sizes
      use atoms
      use iounit
      use math
      use warp
      implicit none
      integer i,nvar
      integer nstep
      real*8 beta
      real*8 e,energy
      real*8 rg,m2ave
      real*8 xx(*)
      character*7 status
c
c
c     translate optimization parameters to coordinates and M2's
c
      nvar = 0
      do i = 1, n
         nvar = nvar + 1
         x(i) = xx(nvar)
         nvar = nvar + 1
         y(i) = xx(nvar)
         nvar = nvar + 1
         z(i) = xx(nvar)
      end do
      do i = 1, n
         nvar = nvar + 1
         m2(i) = abs(xx(nvar))
      end do
c
c     get some info about the current integration step
c
      e = energy ()
      call gyrate (rg)
      m2ave = 0.0d0
      do i = 1, n
         m2ave = m2ave + m2(i)
      end do
      m2ave = m2ave / dble(n)
      write (iout,10)  nstep,log(beta)/logten,e,rg,
     &                 log(m2ave)/logten,status
   10 format (i6,2x,4f13.4,6x,a7)
c
c     save the current coordinates to a disk file
c
      call optsave (nstep,e,xx)
      return
      end
