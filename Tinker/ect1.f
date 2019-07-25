c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine ect1  --  charge transfer energy & derivatives   ##
c     ##                                                              ##
c     ##################################################################
c
c
c
      subroutine ect1
      use sizes
      use energi
      use limits
      use virial
      use warp
      implicit none
c
c     choose the method for summing over pairwise interactions
c
      if (use_ctlist) then
         call ect1b
      else
         call ect1a
      end if
      return
      end 
      
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine ect1  --  CT energy & derivative via double loop ##
c     ##                                                              ##
c     ##################################################################

      subroutine ect1a
      use sizes
      use ctran
      use atomid
      use atoms
      use bound
      use cell
      use couple
      use energi
      use group
      use shunt
      use usage
      use virial
      use deriv
      use inter
      use molcul
      implicit none
      integer i,j,k
      integer ii,it
      integer kk,kt
      real*8 e,de
      real*8 aprec,bexpc
      real*8 aprei,bexpi
      real*8 aprek,bexpk
      real*8 fgrp
      real*8 xr,yr,zr
      real*8 dedx,dedy,dedz 
      real*8 rik,rik2,rik3
      real*8 rik4,rik5,taper
      real*8 dtaper
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8, allocatable :: ctscale(:)
      logical proceed,usei
      character*6 mode
c
c
c     zero out the charge transfer energy contribution
c
c
      ect = 0.0d0
      do i = 1, n
        dect(1,i) = 0.0d0
        dect(2,i) = 0.0d0
        dect(3,i) = 0.0d0
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (ctscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         ctscale(i) = 1.0d0
      end do
c
c     set the coefficients for the switching function
c
      mode = 'CT'
      call switch (mode)
c
c     find the charge transfer energy via double loop search
c
      do ii = 1, nct-1
         i = ict(ii)
         it = jct(i)
         usei = use(i) 
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(i)
            ctscale(i12(j,i)) = ct2scale
         end do
         do j = 1, n13(i)
            ctscale(i13(j,i)) = ct3scale
         end do
         do j = 1, n14(i)
            ctscale(i14(j,i)) = ct4scale
         end do
         do j = 1, n15(i)
            ctscale(i15(j,i)) = ct5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = ii+1, nct
            k = ict(kk)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = jct(k)
               xr = x(i) - x(k)
               yr = y(i) - y(k)
               zr = z(i) - z(k)
               call image (xr,yr,zr)

               rik2 = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
               aprei = abs(apre(it))
               aprek = abs(apre(kt))
               bexpi = abs(bexp(it))
               bexpk = abs(bexp(kt))

               if (ctscale(k) .gt. 0) then
                  !CombinedApre
                  if (aprerule .eq. "GEOMETRIC") then
                     aprec = sqrt(aprei*aprek)
                  else if (aprerule .eq. "ARITHMETIC") then
                     aprec = 0.5d0*(aprei + aprek)
                  end if
                  !CombinedBexp
                  if (bexprule .eq. "GEOMETRIC") then
                     bexpc = sqrt(bexpi*bexpk)
                  else if (bexprule .eq. "ARITHMETIC") then
                     bexpc = 0.5d0*(bexpi + bexpk)
                  end if
               end if
               if (rik2 .le. off2) then
                  rik = sqrt(rik2)
                  e = -aprec*1000.0d0*exp(-bexpc*rik)
                  de = -e*bexpc 

c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                        rik = sqrt(rik2)
                        rik3 = rik2 * rik
                        rik4 = rik2 * rik2
                        rik5 = rik2 * rik3
                        taper = c5*rik5 + c4*rik4 + c3*rik3
     &                             + c2*rik2 + c1*rik + c0
                        dtaper = 5.0d0*c5*rik4 + 4.0d0*c4*rik3
     &                              + 3.0d0*c3*rik2 + 2.0d0*c2*rik + c1
                        de = e*dtaper + de*taper
                        e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                    e = e * fgrp
                    de = de * fgrp
                  end if
c
c     find the chain rule terms for derivative components
c
                  de = de / rik
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
c
c     increment the total charge transfer energy and derivatives
c
                  ect = ect + e
                  dect(1,i) = dect(1,i) + dedx
                  dect(2,i) = dect(2,i) + dedy
                  dect(3,i) = dect(3,i) + dedz
                  dect(1,k) = dect(1,k) - dedx
                  dect(2,k) = dect(2,k) - dedy
                  dect(3,k) = dect(3,k) - dedz
c
c     increment the internal virial tensor components
c
                  vxx = xr * dedx
                  vyx = yr * dedx
                  vzx = zr * dedx
                  vyy = yr * dedy
                  vzy = zr * dedy
                  vzz = zr * dedz
                  vir(1,1) = vir(1,1) + vxx
                  vir(2,1) = vir(2,1) + vyx
                  vir(3,1) = vir(3,1) + vzx
                  vir(1,2) = vir(1,2) + vyx
                  vir(2,2) = vir(2,2) + vyy
                  vir(3,2) = vir(3,2) + vzy
                  vir(1,3) = vir(1,3) + vzx
                  vir(2,3) = vir(2,3) + vzy
                  vir(3,3) = vir(3,3) + vzz
c
c     increment the total intermolecular energy
c
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + e
                  end if
               end if 
            end if
         end do

c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(i)
            ctscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            ctscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            ctscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            ctscale(i15(j,i)) = 1.0d0
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do ii = 1, nct
         i = ict(ii)
         it = jct(i)
         usei = use(i) 
cc
cc     set interaction scaling coefficients for connected atoms
cc
         do j = 1, n12(i)
            ctscale(i12(j,i)) = ct2scale
         end do
         do j = 1, n13(i)
            ctscale(i13(j,i)) = ct3scale
         end do
         do j = 1, n14(i)
            ctscale(i14(j,i)) = ct4scale
         end do
         do j = 1, n15(i)
            ctscale(i15(j,i)) = ct5scale
         end do
cc
cc     decide whether to compute the current interaction
cc
         do kk = ii, nct
            k = ict(kk)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = jct(k)
               do j = 1, ncell
                  xr = x(i) - x(k)
                  yr = y(i) - y(k)
                  zr = z(i) - z(k)
                  call imager (xr,yr,zr,j)
                  rik2 = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
              
                  aprei = abs(apre(it))
                  aprek = abs(apre(kt))
                  bexpi = abs(bexp(it))
                  bexpk = abs(bexp(kt))

                  if (ctscale(k) .gt. 0) then
                     !CombinedApre
                     if (aprerule .eq. "GEOMETRIC") then
                        aprec = sqrt(aprei*aprek)
                     else if (aprerule .eq. "ARITHMETIC") then
                        aprec = 0.5d0*(aprei + aprek)
                     end if
                     !CombinedBexp
                     if (bexprule .eq. "GEOMETRIC") then
                        bexpc = sqrt(bexpi*bexpk)
                     else if (bexprule .eq. "ARITHMETIC") then
                        bexpc = 0.5d0*(bexpi + bexpk)
                     end if
                  end if
              
                  if (rik2 .le. off2) then
                     rik = sqrt(rik2)
                     e = -aprec*1000.0d0*exp(-bexpc*rik)
                     de = -e*bexpc 
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     rik = sqrt(rik2)
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                          + c2*rik2 + c1*rik + c0
                     dtaper = 5.0d0*c5*rik4 + 4.0d0*c4*rik3
     &                           + 3.0d0*c3*rik2 + 2.0d0*c2*rik + c1
                     de = e*dtaper + de*taper
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                 if (use_group) then
                    e = e * fgrp
                    de = de * fgrp
                 end if
c
c     find the chain rule terms for derivative components
c
                     de = de / rik
                     dedx = de * xr
                     dedy = de * yr
                     dedz = de * zr
c
c     increment the total charge transfer energy and derivatives
c
                    ect = ect + e
                    dect(1,i) = dect(1,i) + dedx
                    dect(2,i) = dect(2,i) + dedy
                    dect(3,i) = dect(3,i) + dedz
                    dect(1,k) = dect(1,k) - dedx
                    dect(2,k) = dect(2,k) - dedy
                    dect(3,k) = dect(3,k) - dedz
c
c     increment the internal virial tensor components
c
                    vxx = xr * dedx
                    vyx = yr * dedx
                    vzx = zr * dedx
                    vyy = yr * dedy
                    vzy = zr * dedy
                    vzz = zr * dedz
                    vir(1,1) = vir(1,1) + vxx
                    vir(2,1) = vir(2,1) + vyx
                    vir(3,1) = vir(3,1) + vzx
                    vir(1,2) = vir(1,2) + vyx
                    vir(2,2) = vir(2,2) + vyy
                    vir(3,2) = vir(3,2) + vzy
                    vir(1,3) = vir(1,3) + vzx
                    vir(2,3) = vir(2,3) + vzy
                    vir(3,3) = vir(3,3) + vzz
c
c     increment the total intermolecular energy
c
                     einter = einter + e
                  end if
               end do
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(i)
            ctscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            ctscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            ctscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            ctscale(i15(j,i)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (ctscale)
      return
      end
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine ect1b  --  Charge transfer energy via list    ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "ect1b" calculates the charge transfer energy
c     using a pairwise neighbor list
c
c
      subroutine ect1b
      use sizes
      use atomid
      use atoms
      use bound
      use couple
      use deriv
      use energi
      use group
      use inter
      use molcul
      use neigh
      use shunt
      use usage
      use ctran
      use virial
      implicit none
      integer i,j,k
      integer ii,it
      integer kk,kt
      real*8 e,de
      real*8 aprec,bexpc
      real*8 aprei,bexpi
      real*8 aprek,bexpk
      real*8 fgrp
      real*8 xr,yr,zr
      real*8 dedx,dedy,dedz
      real*8 rik,rik2,rik3
      real*8 rik4,rik5,taper
      real*8 dtaper
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8, allocatable :: ctscale(:)
      logical proceed,usei
      character*6 mode
c
c
c     zero out the CT energy contribution
c
      ect = 0.0d0
      do i = 1, n
         dect(1,i) = 0.0d0
         dect(2,i) = 0.0d0
         dect(3,i) = 0.0d0
      end do
      if (nct .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (ctscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         ctscale(i) = 1.0d0
      end do
c
c     set the coefficients for the switching function
c
      mode = 'CT'
      call switch (mode)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(nct,ict,
!$OMP& jct,use,x,y,z,nctlst,ctlst,n12,n13,n14,n15,
!$OMP& i12,i13,i14,i15,ct2scale,ct3scale,ct4scale,ct5scale,
!$OMP& use_group,off2,apre,bexp,aprerule,bexprule,
!$OMP& cut2,c0,c1,c2,c3,c4,c5,molcule) firstprivate(ctscale)
!$OMP& shared(ect,dect,vir,einter)
!$OMP DO reduction(+:ect,dect,vir,einter) schedule(guided)
c
c     find the van der Waals energy via neighbor list search
c
      do ii = 1, nct
         i = ict(ii)
         it = jct(i)
         usei = use(i)
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            ctscale(i12(j,i)) = ct2scale
         end do
         do j = 1, n13(i)
            ctscale(i13(j,i)) = ct3scale
         end do
         do j = 1, n14(i)
            ctscale(i14(j,i)) = ct4scale
         end do
         do j = 1, n15(i)
            ctscale(i15(j,i)) = ct5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = 1, nctlst(ii)
            k = ict(ctlst(kk,ii))
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = jct(k)
               xr = x(i) - x(k)
               yr = y(i) - y(k)
               zr = z(i) - z(k)
               call image (xr,yr,zr)
               rik2 = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
               aprei = abs(apre(it))
               aprek = abs(apre(kt))
               bexpi = abs(bexp(it))
               bexpk = abs(bexp(kt))

               if (ctscale(k) .gt. 0) then
                  !CombinedApre
                  if (aprerule .eq. "GEOMETRIC") then
                     aprec = sqrt(aprei*aprek)
                  else if (aprerule .eq. "ARITHMETIC") then
                     aprec = 0.5d0*(aprei + aprek)
                  end if
                  !CombinedBexp
                  if (bexprule .eq. "GEOMETRIC") then
                     bexpc = sqrt(bexpi*bexpk)
                  else if (bexprule .eq. "ARITHMETIC") then
                     bexpc = 0.5d0*(bexpi + bexpk)
                  end if
               end if

               if (rik2 .le. off2) then
                  rik = sqrt(rik2)
                  e = -aprec*1000.0d0*exp(-bexpc*rik)
                  de = -e*bexpc
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     rik = sqrt(rik2)
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                          + c2*rik2 + c1*rik + c0
                     dtaper = 5.0d0*c5*rik4 + 4.0d0*c4*rik3
     &                           + 3.0d0*c3*rik2 + 2.0d0*c2*rik + c1
                     de = e*dtaper + de*taper
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then 
                     e = e * fgrp
                     de = de * fgrp
                  end if
c
c     find the chain rule terms for derivative components
c
                  de = de / rik
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
c
c     increment the total charge transfer energy and derivatives
c
                  ect = ect + e
                  dect(1,i) = dect(1,i) + dedx
                  dect(2,i) = dect(2,i) + dedy
                  dect(3,i) = dect(3,i) + dedz
                  dect(1,k) = dect(1,k) - dedx
                  dect(2,k) = dect(2,k) - dedy
                  dect(3,k) = dect(3,k) - dedz
c
c     increment the internal virial tensor components
c
                  vxx = xr * dedx
                  vyx = yr * dedx
                  vzx = zr * dedx
                  vyy = yr * dedy
                  vzy = zr * dedy
                  vzz = zr * dedz
                  vir(1,1) = vir(1,1) + vxx
                  vir(2,1) = vir(2,1) + vyx
                  vir(3,1) = vir(3,1) + vzx
                  vir(1,2) = vir(1,2) + vyx
                  vir(2,2) = vir(2,2) + vyy
                  vir(3,2) = vir(3,2) + vzy
                  vir(1,3) = vir(1,3) + vzx
                  vir(2,3) = vir(2,3) + vzy
                  vir(3,3) = vir(3,3) + vzz
c
c     increment the total intermolecular energy
c
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + e
                  end if
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            ctscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            ctscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            ctscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            ctscale(i15(j,i)) = 1.0d0
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (ctscale)
      return
      end
