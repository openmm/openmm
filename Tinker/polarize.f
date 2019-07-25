c
c
c     #############################################################
c     ##  COPYRIGHT (C) 2001 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  program polarize  --  compute the molecular polarizability  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "polarize" computes the molecular polarizability by applying
c     an external field along each axis followed by diagonalization
c     of the resulting polarizability tensor
c
c
      program polarize
      use sizes
      use atoms
      use inform
      use iounit
      use molcul
      use mpole
      use polar
      use polpot
      use potent
      implicit none
      integer i
      real*8 addu,malpha
      real*8 external
      real*8 exfield(3)
      real*8 umol(3)
      real*8 dalpha(3)
      real*8 alpha(3,3)
      real*8 valpha(3,3)
      character*40 fstr
c
c
c     get the coordinates and required force field parameters
c
      call initial
      call getxyz
      call attach
      call field
      call molecule
      call kmpole
      call kpolar
      call mutate
c
c     sum atomic polarizabilities to get additive molecular value
c
      if (.not. use_polar) then
         write (iout,10)
   10    format (/,' POLARIZE  --  Dipole Polarizability',
     &              ' is Not in Use')
         call fatal
      end if
      addu = 0.0d0
      do i = 1, npole
         addu = polarity(i) + addu
      end do
      fstr = ' Additive Total Polarizability :    '
      if (nmol .eq. 1)  fstr = ' Additive Molecular Polarizability :'
      if (digits .ge. 8) then
         write (iout,20)  fstr(1:36),addu
   20    format (/,a36,f19.8)
      else if (digits .ge. 6) then
         write (iout,30)  fstr(1:36),addu
   30    format (/,a36,f17.6)
      else
         write (iout,40)  fstr(1:36),addu
   40    format (/,a36,f15.4)
      end if
c
c     compute each column of the polarizability tensor
c
      external = 0.01d0
      do i = 1, 3
         exfield(i) = 0.0d0
      end do
      exfield(1) = external
      call moluind (exfield,umol)
      alpha(1,1) = umol(1) / exfield(1)
      alpha(2,1) = umol(2) / exfield(1)
      alpha(3,1) = umol(3) / exfield(1)
      do i = 1, 3
         exfield(i) = 0.0d0
      end do
      exfield(2) = external
      call moluind (exfield,umol)
      alpha(1,2) = umol(1) / exfield(2)
      alpha(2,2) = umol(2) / exfield(2)
      alpha(3,2) = umol(3) / exfield(2)
      do i = 1, 3
         exfield(i) = 0.0d0
      end do
      exfield(3) = external
      call moluind (exfield,umol)
      alpha(1,3) = umol(1) / exfield(3)
      alpha(2,3) = umol(2) / exfield(3)
      alpha(3,3) = umol(3) / exfield(3)
c
c     print out the full polarizability tensor
c
      fstr = ' Total Polarizability Tensor :    '
      if (nmol .eq. 1)  fstr = ' Molecular Polarizability Tensor :'
      write (iout,50)  fstr(1:34)
   50 format (/,a34,/)
      if (digits .ge. 8) then
         write (iout,60)  alpha(1,1),alpha(1,2),alpha(1,3),
     &                    alpha(2,1),alpha(2,2),alpha(2,3),
     &                    alpha(3,1),alpha(3,2),alpha(3,3)
   60    format (15x,3f16.8,/,15x,3f16.8,/,15x,3f16.8)
      else if (digits .ge. 6) then
         write (iout,70)  alpha(1,1),alpha(1,2),alpha(1,3),
     &                    alpha(2,1),alpha(2,2),alpha(2,3),
     &                    alpha(3,1),alpha(3,2),alpha(3,3)
   70    format (15x,3f14.6,/,15x,3f14.6,/,15x,3f14.6)
      else
         write (iout,80)  alpha(1,1),alpha(1,2),alpha(1,3),
     &                    alpha(2,1),alpha(2,2),alpha(2,3),
     &                    alpha(3,1),alpha(3,2),alpha(3,3)
   80    format (15x,3f12.4,/,15x,3f12.4,/,15x,3f12.4)
      end if
c
c     diagonalize the tensor and get molecular polarizability
c
      call jacobi (3,alpha,dalpha,valpha)
      if (nmol .eq. 1)  fstr = ' Polarizability Tensor Eigenvalues :'
      write (iout,90)  fstr(1:36)
   90 format (/,a36,/)
      if (digits .ge. 8) then
         write (iout,100)  dalpha(1),dalpha(2),dalpha(3)
  100    format (15x,3f16.8)
      else if (digits .ge. 6) then
         write (iout,110)  dalpha(1),dalpha(2),dalpha(3)
  110    format (15x,3f14.6)
      else
         write (iout,120)  dalpha(1),dalpha(2),dalpha(3)
  120    format (15x,3f12.4)
      end if
      malpha = (dalpha(1)+dalpha(2)+dalpha(3)) / 3.0d0
      fstr = ' Interactive Total Polarizability :    '
      if (nmol .eq. 1)  fstr = ' Interactive Molecular Polarizability :'
      if (digits .ge. 8) then
         write (iout,130)  fstr(1:39),malpha
  130    format (/,a39,f16.8)
      else if (digits .ge. 6) then
         write (iout,140)  fstr(1:39),malpha
  140    format (/,a39,f14.6)
      else
         write (iout,150)  fstr(1:39),malpha
  150    format (/,a39,f12.4)
      end if
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine moluind  --  molecular induced dipole in field  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "moluind" computes the molecular induced dipole components
c     in the presence of an external electric field
c
c
      subroutine moluind (exfield,umol)
      use sizes
      use atoms
      use inform
      use iounit
      use mpole
      use polar
      use polpot
      use units
      implicit none
      integer i,j,iter
      integer maxiter
      real*8 eps,epsold
      real*8 polmin
      real*8 a,b,sum
      real*8 umol(3)
      real*8 exfield(3)
      real*8, allocatable :: poli(:)
      real*8, allocatable :: field(:,:)
      real*8, allocatable :: rsd(:,:)
      real*8, allocatable :: zrsd(:,:)
      real*8, allocatable :: conj(:,:)
      real*8, allocatable :: vec(:,:)
      logical done
c
c
c     set induced dipoles to polarizability times external field
c
      do i = 1, npole
         do j = 1, 3
            uind(j,i) = polarity(i) * exfield(j)
         end do
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (poli(npole))
      allocate (field(3,npole))
      allocate (rsd(3,npole))
      allocate (zrsd(3,npole))
      allocate (conj(3,npole))
      allocate (vec(3,npole))
c
c     compute mutual induced dipole moments via CG algorithm
c
      if (poltyp .eq. 'MUTUAL') then
         done = .false.
         maxiter = 500
         iter = 0
         eps = 100.0d0
         polmin = 0.00000001d0
         call ufield (field)
         do i = 1, npole
            poli(i) = max(polmin,polarity(i))
            do j = 1, 3
               rsd(j,i) = field(j,i)
               zrsd(j,i) = rsd(j,i) * poli(i)
               conj(j,i) = zrsd(j,i)
            end do
         end do
c
c     iterate the mutual induced dipoles and check convergence
c
         do while (.not. done)
            iter = iter + 1
            do i = 1, npole
               do j = 1, 3
                  vec(j,i) = uind(j,i)
                  uind(j,i) = conj(j,i)
               end do
            end do
            call ufield (field)
            do i = 1, npole
               do j = 1, 3
                  uind(j,i) = vec(j,i)
                  vec(j,i) = conj(j,i)/poli(i) - field(j,i)
               end do
            end do
            a = 0.0d0
            sum = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  a = a + conj(j,i)*vec(j,i)
                  sum = sum + rsd(j,i)*zrsd(j,i)
               end do
            end do
            if (a .ne. 0.0d0)  a = sum / a
            do i = 1, npole
               do j = 1, 3
                  uind(j,i) = uind(j,i) + a*conj(j,i)
                  rsd(j,i) = rsd(j,i) - a*vec(j,i)
               end do
            end do
            b = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  zrsd(j,i) = rsd(j,i) * poli(i)
                  b = b + rsd(j,i)*zrsd(j,i)
               end do
            end do
            if (sum .ne. 0.0d0)  b = b / sum
            eps = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  conj(j,i) = zrsd(j,i) + b*conj(j,i)
                  eps = eps + rsd(j,i)*rsd(j,i)
               end do
            end do
            eps = debye * sqrt(eps/dble(npolar))
            epsold = eps
            if (debug) then
               if (iter .eq. 1) then
                  write (iout,10)
   10             format (/,' Determination of Induced Dipole',
     &                       ' Moments :',
     &                    //,4x,'Iter',8x,'RMS Change (Debyes)',/)
               end if
               write (iout,20)  iter,eps
   20          format (i8,7x,f16.10)
            end if
            if (eps .lt. poleps)  done = .true.
            if (eps .gt. epsold)  done = .true.
            if (iter .ge. maxiter)  done = .true.
         end do
c
c     print a warning if induced dipoles failed to converge
c
         if (eps .gt. poleps) then
            write (iout,30)
   30       format (/,' MOLUIND  --  Warning, Induced Dipoles',
     &                 ' are not Converged')
         end if
      end if
c
c     sum up the total molecular induced dipole components
c
      do j = 1, 3
         umol(j) = 0.0d0
      end do
      do i = 1, npole
         umol(1) = umol(1) + uind(1,i)
         umol(2) = umol(2) + uind(2,i)
         umol(3) = umol(3) + uind(3,i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (poli)
      deallocate (field)
      deallocate (rsd)
      deallocate (zrsd)
      deallocate (conj)
      deallocate (vec)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine ufield  --  electric field from induced dipoles  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "ufield" finds the field at each polarizable site due to the
c     induced dipoles at the other sites using Thole's method to
c     damp the field at close range
c
c
      subroutine ufield (field)
      use sizes
      use atoms
      use mpole
      use polar
      use polgrp
      use polpot
      implicit none
      integer i,j,k
      integer ii,kk
      real*8 r,r2,rr3,rr5
      real*8 xr,yr,zr
      real*8 uir,ukr
      real*8 uix,uiy,uiz
      real*8 ukx,uky,ukz
      real*8 pdi,pti
      real*8 pgamma,damp
      real*8 expdamp
      real*8 scale3,scale5
      real*8 fi(3),fk(3)
      real*8, allocatable :: dscale(:)
      real*8 field(3,*)
c
c
c     zero out the value of the electric field at each site
c
      do i = 1, npole
         do j = 1, 3
            field(j,i) = 0.0d0
         end do
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (dscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         dscale(i) = 1.0d0
      end do
c
c     loop over pairs of sites incrementing the electric field
c
      do i = 1, npole-1
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         uix = uind(1,i)
         uiy = uind(2,i)
         uiz = uind(3,i)
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = u4scale
         end do
         do k = i+1, npole
            kk = ipole(k)
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            r2 = xr*xr + yr* yr + zr*zr
            ukx = uind(1,k)
            uky = uind(2,k)
            ukz = uind(3,k)
            r = sqrt(r2)
            uir = xr*uix + yr*uiy + zr*uiz
            ukr = xr*ukx + yr*uky + zr*ukz
c
c     adjust the field to account for polarization damping
c
            scale3 = dscale(kk)
            scale5 = dscale(kk)
            damp = pdi * pdamp(k)
            if (damp .ne. 0.0d0) then
               pgamma = min(pti,thole(k))
               damp = -pgamma * (r/damp)**3
               if (damp .gt. -50.0d0) then
                  expdamp = exp(damp)
                  scale3 = scale3 * (1.0d0-expdamp)
                  scale5 = scale5 * (1.0d0-expdamp*(1.0d0-damp))
               end if
            end if
            rr3 = scale3 / (r*r2)
            rr5 = 3.0d0 * scale5 / (r*r2*r2)
            fi(1) = -rr3*ukx + rr5*ukr*xr
            fi(2) = -rr3*uky + rr5*ukr*yr
            fi(3) = -rr3*ukz + rr5*ukr*zr
            fk(1) = -rr3*uix + rr5*uir*xr
            fk(2) = -rr3*uiy + rr5*uir*yr
            fk(3) = -rr3*uiz + rr5*uir*zr
c
c     increment the field at each site due to this interaction
c
            do j = 1, 3
               field(j,i) = field(j,i) + fi(j)
               field(j,k) = field(j,k) + fk(j)
            end do
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (dscale)
      return
      end
