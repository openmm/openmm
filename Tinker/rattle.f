c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1991  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine rattle  --  RATTLE distance constraint method  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "rattle" implements the first portion of the RATTLE algorithm
c     by correcting atomic positions and half-step velocities to
c     maintain interatomic distance and absolute spatial constraints
c
c     literature reference:
c
c     H. C. Andersen, "RATTLE: A Velocity Version of the SHAKE
c     Algorithm for Molecular Dynamics Calculations", Journal of
c     Computational Physics, 52, 24-34 (1983)
c
c
      subroutine rattle (dt,xold,yold,zold)
      use sizes
      use atomid
      use atoms
      use freeze
      use group
      use inform
      use iounit
      use moldyn
      use usage
      implicit none
      integer i,j,k
      integer ia,ib,mode
      integer niter,maxiter
      integer start,stop
      real*8 dt,eps,sor
      real*8 xr,yr,zr
      real*8 xo,yo,zo
      real*8 xv,yv,zv
      real*8 dot,rma,rmb
      real*8 weigh,dist2
      real*8 delta,term
      real*8 xterm,yterm,zterm
      real*8 xold(*)
      real*8 yold(*)
      real*8 zold(*)
      logical done
      logical, allocatable :: moved(:)
      logical, allocatable :: update(:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (moved(n))
      allocate (update(n))
c
c     initialize the lists of atoms previously corrected
c
      do i = 1, n
         if (use(i)) then
            moved(i) = .true.
         else
            moved(i) = .false.
         end if
         update(i) = .false.
      end do
c
c     set the iteration counter, termination and tolerance
c
      maxiter = 500
      sor = 1.25d0
      eps = rateps
c
c     apply RATTLE to distances and half-step velocity values
c
      niter = 0
      done = .false.
      do while (.not.done .and. niter.lt.maxiter)
         niter = niter + 1
         done = .true.
         do i = 1, nrat
            ia = irat(1,i)
            ib = irat(2,i)
            if (moved(ia) .or. moved(ib)) then
               xr = x(ib) - x(ia)
               yr = y(ib) - y(ia)
               zr = z(ib) - z(ia)
               if (ratimage(i))  call image (xr,yr,zr)
               dist2 = xr**2 + yr**2 + zr**2
               delta = krat(i)**2 - dist2
               if (abs(delta) .gt. eps) then
                  done = .false.
                  update(ia) = .true.
                  update(ib) = .true.
                  xo = xold(ib) - xold(ia)
                  yo = yold(ib) - yold(ia)
                  zo = zold(ib) - zold(ia)
                  if (ratimage(i))  call image (xo,yo,zo)
                  dot = xr*xo + yr*yo + zr*zo
                  rma = 1.0d0 / mass(ia)
                  rmb = 1.0d0 / mass(ib)
                  term = sor * delta / (2.0d0 * (rma+rmb) * dot)
                  xterm = xo * term
                  yterm = yo * term
                  zterm = zo * term
                  x(ia) =  x(ia) - xterm*rma
                  y(ia) =  y(ia) - yterm*rma
                  z(ia) =  z(ia) - zterm*rma
                  x(ib) =  x(ib) + xterm*rmb
                  y(ib) =  y(ib) + yterm*rmb
                  z(ib) =  z(ib) + zterm*rmb
                  rma = rma / dt
                  rmb = rmb / dt
                  v(1,ia) = v(1,ia) - xterm*rma
                  v(2,ia) = v(2,ia) - yterm*rma
                  v(3,ia) = v(3,ia) - zterm*rma
                  v(1,ib) = v(1,ib) + xterm*rmb
                  v(2,ib) = v(2,ib) + yterm*rmb
                  v(3,ib) = v(3,ib) + zterm*rmb
               end if
            end if
         end do
         do i = 1, n
            moved(i) = update(i)
            update(i) = .false.
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (moved)
      deallocate (update)
c
c     write information on the number of iterations needed
c
      if (niter .eq. maxiter) then
         write (iout,10)
   10    format (/,' RATTLE  --  Warning, Distance Constraints',
     &              ' not Satisfied')
         call prterr
         call fatal
      else if (debug) then
         write (iout,20)  niter
   20    format (' RATTLE   --  Distance Constraints met at',i6,
     &              ' Iterations')
      end if
c
c     apply group position and velocity constraints via exact reset
c
      do i = 1, nratx
         ia = iratx(i)
         mode = kratx(i)
         xr = 0.0d0
         yr = 0.0d0
         zr = 0.0d0
         xv = 0.0d0
         yv = 0.0d0
         zv = 0.0d0
         start = igrp(1,ia)
         stop = igrp(2,ia)
         do j = start, stop
            k = kgrp(j)
            weigh = mass(k) / grpmass(ia)
            if (mode .gt. 2) then
               xr = xr + x(k)*weigh
               xv = xv + v(1,k)*weigh
            end if
            if (mode .gt. 1) then
               yr = yr + y(k)*weigh
               yv = yv + v(2,k)*weigh
            end if
            zr = zr + z(k)*weigh
            zv = zv + v(3,k)*weigh
         end do
         do j = start, stop
            k = kgrp(j)
            x(k) =  x(k) - xr
            y(k) =  y(k) - yr
            z(k) =  z(k) - zr
            v(1,k) = v(1,k) - xv
            v(2,k) = v(2,k) - yv
            v(3,k) = v(3,k) - zv
         end do
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine rattle2  --  RATTLE atom velocity corrections  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "rattle2" implements the second portion of the RATTLE algorithm
c     by correcting the full-step velocities in order to maintain
c     interatomic distance constraints
c
c
      subroutine rattle2 (dt)
      use sizes
      use atomid
      use atoms
      use freeze
      use group
      use inform
      use iounit
      use moldyn
      use units
      use usage
      use virial
      implicit none
      integer i,j,k
      integer ia,ib,mode
      integer niter,maxiter
      integer start,stop
      real*8 dt,eps,sor
      real*8 xr,yr,zr
      real*8 xv,yv,zv
      real*8 dot,rma,rmb
      real*8 weigh,vterm,term
      real*8 xterm,yterm,zterm
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      logical done
      logical, allocatable :: moved(:)
      logical, allocatable :: update(:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (moved(n))
      allocate (update(n))
c
c     initialize the lists of atoms previously corrected
c
      do i = 1, n
         if (use(i)) then
            moved(i) = .true.
         else
            moved(i) = .false.
         end if
         update(i) = .false.
      end do
c
c     set the iteration counter, termination and tolerance
c
      maxiter = 500
      niter = 0
      done = .false.
      sor = 1.25d0
      eps = rateps / dt
      vterm = 2.0d0 / (dt*convert)
c
c     apply the RATTLE algorithm to correct the velocities
c
      do while (.not.done .and. niter.lt.maxiter)
         niter = niter + 1
         done = .true.
         do i = 1, nrat
            ia = irat(1,i)
            ib = irat(2,i)
            if (moved(ia) .or. moved(ib)) then
               xr = x(ib) - x(ia)
               yr = y(ib) - y(ia)
               zr = z(ib) - z(ia)
               if (ratimage(i))  call image (xr,yr,zr)
               xv = v(1,ib) - v(1,ia)
               yv = v(2,ib) - v(2,ia)
               zv = v(3,ib) - v(3,ia)
               dot = xr*xv + yr*yv + zr*zv
               rma = 1.0d0 / mass(ia)
               rmb = 1.0d0 / mass(ib)
               term = -dot / ((rma+rmb) * krat(i)**2)
               if (abs(term) .gt. eps) then
                  done = .false.
                  update(ia) = .true.
                  update(ib) = .true.
                  term = sor * term
                  xterm = xr * term
                  yterm = yr * term
                  zterm = zr * term
                  v(1,ia) = v(1,ia) - xterm*rma
                  v(2,ia) = v(2,ia) - yterm*rma
                  v(3,ia) = v(3,ia) - zterm*rma
                  v(1,ib) = v(1,ib) + xterm*rmb
                  v(2,ib) = v(2,ib) + yterm*rmb
                  v(3,ib) = v(3,ib) + zterm*rmb
c
c     increment the internal virial tensor components
c
                  xterm = xterm * vterm
                  yterm = yterm * vterm
                  zterm = zterm * vterm
                  vxx = xr * xterm
                  vyx = yr * xterm
                  vzx = zr * xterm
                  vyy = yr * yterm
                  vzy = zr * yterm
                  vzz = zr * zterm
                  vir(1,1) = vir(1,1) - vxx
                  vir(2,1) = vir(2,1) - vyx
                  vir(3,1) = vir(3,1) - vzx
                  vir(1,2) = vir(1,2) - vyx
                  vir(2,2) = vir(2,2) - vyy
                  vir(3,2) = vir(3,2) - vzy
                  vir(1,3) = vir(1,3) - vzx
                  vir(2,3) = vir(2,3) - vzy
                  vir(3,3) = vir(3,3) - vzz
               end if
            end if
         end do
         do i = 1, n
            moved(i) = update(i)
            update(i) = .false.
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (moved)
      deallocate (update)
c
c     write information on the number of iterations needed
c
      if (niter .eq. maxiter) then
         write (iout,10)
   10    format (/,' RATTLE2  --  Warning, Velocity Constraints',
     &              ' not Satisfied')
         call prterr
         call fatal
      else if (debug) then
         write (iout,20)  niter
   20    format (' RATTLE2  --  Velocity Constraints met at',i6,
     &              ' Iterations')
      end if
c
c     apply any atom group velocity constraints via exact reset
c
      do i = 1, nratx
         ia = iratx(i)
         mode = kratx(i)
         xv = 0.0d0
         yv = 0.0d0
         zv = 0.0d0
         start = igrp(1,ia)
         stop = igrp(2,ia)
         do j = start, stop
            k = kgrp(j)
            weigh = mass(k) / grpmass(ia)
            if (mode .gt. 2) then
               xv = xv + v(1,k)*weigh
            end if
            if (mode .gt. 1) then
               yv = yv + v(2,k)*weigh
            end if
            zv = zv + v(3,k)*weigh
         end do
         do j = start, stop
            k = kgrp(j)
            v(1,k) = v(1,k) - xv
            v(2,k) = v(2,k) - yv
            v(3,k) = v(3,k) - zv
         end do
      end do
      return
      end
