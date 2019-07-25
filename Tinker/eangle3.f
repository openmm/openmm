c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine eangle3  --  angle bending energy & analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "eangle3" calculates the angle bending potential energy, also
c     partitions the energy among the atoms; projected in-plane
c     angles at trigonal centers, spceial linear or Fourier angle
c     bending terms are optionally used
c
c
      subroutine eangle3
      use sizes
      use action
      use analyz
      use angbnd
      use angpot
      use atomid
      use atoms
      use bound
      use energi
      use group
      use inform
      use iounit
      use math
      use usage
      implicit none
      integer i,ia,ib,ic,id
      real*8 e,ideal,force
      real*8 fold,factor
      real*8 dot,cosine
      real*8 angle,fgrp
      real*8 dt,dt2,dt3,dt4
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xab,yab,zab
      real*8 xcb,ycb,zcb
      real*8 xad,yad,zad
      real*8 xbd,ybd,zbd
      real*8 xcd,ycd,zcd
      real*8 xip,yip,zip
      real*8 xap,yap,zap
      real*8 xcp,ycp,zcp
      real*8 rab2,rcb2
      real*8 rap2,rcp2
      real*8 xt,yt,zt
      real*8 rt2,delta
      logical proceed
      logical header,huge
      character*9 label
c
c
c     zero out the angle bending energy and partitioning terms
c
      nea = 0
      ea = 0.0d0
      do i = 1, n
         aea(i) = 0.0d0
      end do
      if (nangle .eq. 0)  return
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. nangle.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Angle Bending Interactions :',
     &           //,' Type',18x,'Atom Names',18x,
     &              'Ideal',4x,'Actual',6x,'Energy',/)
      end if
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(nangle,iang,anat,ak,afld,use,
!$OMP& x,y,z,cang,qang,pang,sang,angtyp,angunit,use_group,use_polymer,
!$OMP& name,verbose,debug,header,iout)
!$OMP& shared(ea,nea,aea)
!$OMP DO reduction(+:ea,nea,aea) schedule(guided)
c
c     calculate the bond angle bending energy term
c
      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         id = iang(4,i)
         ideal = anat(i)
         force = ak(i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (angtyp(i) .eq. 'IN-PLANE') then
            if (use_group)  call groups (proceed,fgrp,ia,ib,ic,id,0,0)
            if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                                 use(ic) .or. use(id))
         else
            if (use_group)  call groups (proceed,fgrp,ia,ib,ic,0,0,0)
            if (proceed)  proceed = (use(ia) .or. use(ib) .or. use(ic))
         end if
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
c     compute the bond angle bending energy
c
            if (angtyp(i) .ne. 'IN-PLANE') then
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
               if (rab2.ne.0.0d0 .and. rcb2.ne.0.0d0) then
                  dot = xab*xcb + yab*ycb + zab*zcb
                  cosine = dot / sqrt(rab2*rcb2)
                  cosine = min(1.0d0,max(-1.0d0,cosine))
                  angle = radian * acos(cosine)
                  if (angtyp(i) .eq. 'HARMONIC') then
                     dt = angle - ideal
                     dt2 = dt * dt
                     dt3 = dt2 * dt
                     dt4 = dt2 * dt2
                     e = angunit * force * dt2
     &                      * (1.0d0+cang*dt+qang*dt2+pang*dt3+sang*dt4)
                  else if (angtyp(i) .eq. 'LINEAR') then
                     factor = 2.0d0 * angunit * radian**2
                     e = factor * force * (1.0d0+cosine)
                  else if (angtyp(i) .eq. 'FOURIER') then
                     fold = afld(i)
                     factor = 2.0d0 * angunit * (radian/fold)**2
                     cosine = cos((fold*angle-ideal)/radian)
                     e = factor * force * (1.0d0+cosine)
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the total bond angle bending energy
c
                  nea = nea + 1
                  ea = ea + e
                  aea(ib) = aea(ib) + e
c
c     print a message if the energy of this interaction is large
c
                  huge = (e .gt. 5.0d0)
                  if (debug .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,20)
   20                   format (/,' Individual Angle Bending',
     &                             ' Interactions :',
     &                          //,' Type',18x,'Atom Names',18x,
     &                             'Ideal',4x,'Actual',6x,'Energy',/)
                     end if
                     label = 'Angle    '
                     if (angtyp(i) .eq. 'LINEAR') then
                        label = 'Angle-Lin'
                     else if (angtyp(i) .eq. 'FOURIER') then
                        label = 'Angle-Cos'
                        ideal = (ideal+180.0d0) / fold
                        if (angle-ideal .gt. 180.0d0/fold)
     &                     ideal = ideal + 360.0d0/fold
                     end if
                     write (iout,30)  label,ia,name(ia),ib,name(ib),
     &                                ic,name(ic),ideal,angle,e
   30                format (1x,a9,1x,i7,'-',a3,i7,'-',a3,i7,
     &                          '-',a3,2x,2f10.4,f12.4)
                  end if
               end if
c
c     compute the projected in-plane angle bend energy
c
            else
               xid = x(id)
               yid = y(id)
               zid = z(id)
               xad = xia - xid
               yad = yia - yid
               zad = zia - zid
               xbd = xib - xid
               ybd = yib - yid
               zbd = zib - zid
               xcd = xic - xid
               ycd = yic - yid
               zcd = zic - zid
               if (use_polymer) then
                  call image (xad,yad,zad)
                  call image (xbd,ybd,zbd)
                  call image (xcd,ycd,zcd)
               end if
               xt = yad*zcd - zad*ycd
               yt = zad*xcd - xad*zcd
               zt = xad*ycd - yad*xcd
               rt2 = xt*xt + yt*yt + zt*zt
               delta = -(xt*xbd + yt*ybd + zt*zbd) / rt2
               xip = xib + xt*delta
               yip = yib + yt*delta
               zip = zib + zt*delta
               xap = xia - xip
               yap = yia - yip
               zap = zia - zip
               xcp = xic - xip
               ycp = yic - yip
               zcp = zic - zip
               if (use_polymer) then
                  call image (xap,yap,zap)
                  call image (xcp,ycp,zcp)
               end if
               rap2 = xap*xap + yap*yap + zap*zap
               rcp2 = xcp*xcp + ycp*ycp + zcp*zcp
               if (rap2.ne.0.0d0 .and. rcp2.ne.0.0d0) then
                  dot = xap*xcp + yap*ycp + zap*zcp
                  cosine = dot / sqrt(rap2*rcp2)
                  cosine = min(1.0d0,max(-1.0d0,cosine))
                  angle = radian * acos(cosine)
                  dt = angle - ideal
                  dt2 = dt * dt
                  dt3 = dt2 * dt
                  dt4 = dt2 * dt2
                  e = angunit * force * dt2
     &                   * (1.0d0+cang*dt+qang*dt2+pang*dt3+sang*dt4)
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the total bond angle bending energy
c
                  nea = nea + 1
                  ea = ea + e
                  aea(ib) = aea(ib) + e
c
c     print a message if the energy of this interaction is large
c
                  huge = (e .gt. 5.0d0)
                  if (debug .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,40)
   40                   format (/,' Individual Angle Bending',
     &                             ' Interactions :',
     &                          //,' Type',18x,'Atom Names',18x,
     &                             'Ideal',4x,'Actual',6x,'Energy',/)
                     end if
                     write (iout,50)  ia,name(ia),ib,name(ib),ic,
     &                                name(ic),ideal,angle,e
   50                format (' Angle-IP',2x,i7,'-',a3,i7,'-',a3,i7,
     &                          '-',a3,2x,2f10.4,f12.4)
                  end if
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
