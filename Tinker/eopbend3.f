c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine eopbend3  --  out-of-plane bending & analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "eopbend3" computes the out-of-plane bend potential energy at
c     trigonal centers via a Wilson-Decius-Cross or Allinger angle;
c     also partitions the energy among the atoms
c
c
      subroutine eopbend3
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
      use opbend
      use usage
      implicit none
      integer i,iopbend
      integer ia,ib,ic,id
      real*8 e,angle,force
      real*8 cosine,fgrp
      real*8 dt,dt2,dt3,dt4
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xab,yab,zab
      real*8 xcb,ycb,zcb
      real*8 xdb,ydb,zdb
      real*8 xad,yad,zad
      real*8 xcd,ycd,zcd
      real*8 rdb2,rad2,rcd2
      real*8 rab2,rcb2
      real*8 cc,ee,bkk2
      logical proceed
      logical header,huge
c
c
c     zero out the out-of-plane bend energy and partitioning
c
      neopb = 0
      eopb = 0.0d0
      do i = 1, n
         aeopb(i) = 0.0d0
      end do
      if (nopbend .eq. 0)  return
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. nopbend.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Out-of-Plane Bend Interactions :',
     &           //,' Type',25x,'Atom Names',21x,'Angle',
     &              6x,'Energy',/)
      end if
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(nopbend,iopb,iang,opbk,use,
!$OMP& x,y,z,opbtyp,copb,qopb,popb,sopb,opbunit,use_group,use_polymer,
!$OMP& name,verbose,debug,header,iout)
!$OMP& shared(eopb,neopb,aeopb)
!$OMP DO reduction(+:eopb,neopb,aeopb) schedule(guided)
c
c     calculate the out-of-plane bending energy term
c
      do iopbend = 1, nopbend
         i = iopb(iopbend)
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         id = iang(4,i)
         force = opbk(iopbend)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,ia,ib,ic,id,0,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
c
c     get the coordinates of the atoms at trigonal center
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
            xid = x(id)
            yid = y(id)
            zid = z(id)
c
c     compute the out-of-plane bending angle
c
            xab = xia - xib
            yab = yia - yib
            zab = zia - zib
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdb = xid - xib
            ydb = yid - yib
            zdb = zid - zib
            xad = xia - xid
            yad = yia - yid
            zad = zia - zid
            xcd = xic - xid
            ycd = yic - yid
            zcd = zic - zid
            if (use_polymer) then
               call image (xab,yab,zab)
               call image (xcb,ycb,zcb)
               call image (xdb,ydb,zdb)
               call image (xad,yad,zad)
               call image (xcd,ycd,zcd)
            end if
c
c     W-D-C angle between A-B-C plane and B-D vector for D-B<AC
c
            if (opbtyp .eq. 'W-D-C') then
               rab2 = xab*xab + yab*yab + zab*zab
               rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
               cc = rab2*rcb2 - (xab*xcb+yab*ycb+zab*zcb)**2
c
c     Allinger angle between A-C-D plane and D-B vector for D-B<AC
c
            else if (opbtyp .eq. 'ALLINGER') then
               rad2 = xad*xad + yad*yad + zad*zad
               rcd2 = xcd*xcd + ycd*ycd + zcd*zcd
               cc = rad2*rcd2 - (xad*xcd+yad*ycd+zad*zcd)**2
            end if
c
c     find the out-of-plane angle bending energy
c
            ee = xdb*(yab*zcb-zab*ycb) + ydb*(zab*xcb-xab*zcb)
     &              + zdb*(xab*ycb-yab*xcb)
            rdb2 = xdb*xdb + ydb*ydb + zdb*zdb
            if (rdb2.ne.0.0d0 .and. cc.ne.0.0d0) then
               bkk2 = rdb2 - ee*ee/cc
               cosine = sqrt(bkk2/rdb2)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
               dt = angle
               dt2 = dt * dt
               dt3 = dt2 * dt
               dt4 = dt2 * dt2
               e = opbunit * force * dt2
     &                * (1.0d0+copb*dt+qopb*dt2+popb*dt3+sopb*dt4)
c
c     scale the interaction based on its group membership
c
               if (use_group)  e = e * fgrp
c
c     increment the total out-of-plane bending energy
c
               neopb = neopb + 1
               eopb = eopb + e
               aeopb(ib) = aeopb(ib) + e
c
c     print a message if the energy of this interaction is large
c
               huge = (e .gt. 5.0d0)
               if (debug .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,20)
   20                format (/,' Individual Out-of-Plane Bend',
     &                          ' Interactions :',
     &                       //,' Type',25x,'Atom Names',21x,'Angle',
     &                          6x,'Energy',/)
                  end if
                  write (iout,30)  id,name(id),ib,name(ib),ia,
     &                             name(ia),ic,name(ic),angle,e
   30             format (' O-P-Bend',2x,4(i7,'-',a3),f11.4,f12.4)
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
