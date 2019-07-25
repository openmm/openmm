c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine estrbnd3  --  stretch-bend energy & analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "estrbnd3" calculates the stretch-bend potential energy;
c     also partitions the energy among the atoms
c
c
      subroutine estrbnd3
      use sizes
      use action
      use analyz
      use angbnd
      use angpot
      use atomid
      use atoms
      use bndstr
      use bound
      use energi
      use group
      use inform
      use iounit
      use math
      use strbnd
      use usage
      implicit none
      integer i,j,k,istrbnd
      integer ia,ib,ic
      real*8 e,dt
      real*8 dr1,dr2
      real*8 fgrp,angle
      real*8 force1,force2
      real*8 dot,cosine
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xab,yab,zab
      real*8 xcb,ycb,zcb
      real*8 rab,rab2
      real*8 rcb,rcb2
      logical proceed
      logical header,huge
c
c
c     zero out the energy component and partitioning terms
c
      neba = 0
      eba = 0.0d0
      do i = 1, n
         aeba(i) = 0.0d0
      end do
      if (nstrbnd .eq. 0)  return
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. nstrbnd.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Stretch-Bend Interactions :',
     &           //,' Type',18x,'Atom Names',18x,'dSB 1',
     &              5x,'dSB 2',6x,'Energy',/)
      end if
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(nstrbnd,isb,iang,sbk,
!$OMP& anat,bl,bk,use,x,y,z,stbnunit,use_group,use_polymer,
!$OMP& name,verbose,debug,header,iout)
!$OMP& shared(eba,neba,aeba)
!$OMP DO reduction(+:eba,neba,aeba) schedule(guided)
c
c     calculate the stretch-bend interaction energy term
c
      do istrbnd = 1, nstrbnd
         i = isb(1,istrbnd)
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         force1 = sbk(1,istrbnd)
         force2 = sbk(2,istrbnd)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,ia,ib,ic,0,0,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or. use(ic))
c
c     get the coordinates of the atoms in the angle
c
         if (proceed) then
            xia = x(ia)
            yia = y(ia)
            zia = z(ia)
            xib = x(ib)
            yib = y(ib)
            zib = z(ib)
            xic = x(ic)
            yic = y(ic)
            zic = z(ic)
c
c     compute the value of the bond angle
c
            xab = xia - xib
            yab = yia - yib
            zab = zia - zib
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            if (use_polymer) then
               call image (xab,yab,zab)
               call image (xcb,ycb,zcb)
            end if
            rab2 = xab*xab + yab*yab + zab*zab
            rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
            if (min(rab2,rcb2) .ne. 0.0d0) then
               rab = sqrt(rab2)
               rcb = sqrt(rcb2)
               dot = xab*xcb + yab*ycb + zab*zcb
               cosine = dot / (rab*rcb)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
               dt = angle - anat(i)
c
c     get the stretch-bend interaction energy
c
               j = isb(2,istrbnd)
               k = isb(3,istrbnd)
               dr1 = rab - bl(j)
               dr2 = rcb - bl(k)
               e = stbnunit * (force1*dr1+force2*dr2) * dt
c
c     scale the interaction based on its group membership
c
               if (use_group)  e = e * fgrp
c
c     increment the total stretch-bend energy
c
               neba = neba + 1
               eba = eba + e
               aeba(ib) = aeba(ib) + e
c
c     print a message if the energy of this interaction is large
c
               huge = (abs(e) .gt. 5.0d0)
               if (debug .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,20)
   20                format (/,' Individual Stretch-Bend',
     &                          ' Interactions :',
     &                       //,' Type',18x,'Atom Names',18x,'dSB 1',
     &                          5x,'dSB 2',6x,'Energy',/)
                  end if
                  write (iout,30)  ia,name(ia),ib,name(ib),
     &                             ic,name(ic),dr1*dt,dr2*dt,e
   30             format (' StrBend',3x,3(i7,'-',a3),2x,2f10.4,f12.4)
               end if
            end if
         end if
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end
