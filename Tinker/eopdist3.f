c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine eopdist3  --  out-of-plane distance & analysis  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "eopdist3" computes the out-of-plane distance potential energy
c     at trigonal centers via the central atom height; also partitions
c     the energy among the atoms
c
c
      subroutine eopdist3
      use sizes
      use action
      use analyz
      use angpot
      use atomid
      use atoms
      use bound
      use energi
      use group
      use inform
      use iounit
      use opdist
      use usage
      implicit none
      integer i,ia,ib,ic,id
      real*8 e,force,fgrp
      real*8 dt,dt2,dt3,dt4
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xad,yad,zad
      real*8 xbd,ybd,zbd
      real*8 xcd,ycd,zcd
      real*8 xt,yt,zt,rt2
      logical proceed
      logical header,huge
c
c
c     zero out the out-of-plane distance energy and partitioning
c
      neopd = 0
      eopd = 0.0d0
      do i = 1, n
         aeopd(i) = 0.0d0
      end do
      if (nopdist .eq. 0)  return
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. nopdist.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Out-of-Plane Distance Interactions :',
     &           //,' Type',25x,'Atom Names',18x,'Distance',
     &              6x,'Energy',/)
      end if
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(nopdist,iopd,opdk,use,
!$OMP& x,y,z,copd,qopd,popd,sopd,opdunit,use_group,use_polymer,
!$OMP& name,verbose,debug,header,iout)
!$OMP& shared(eopd,neopd,aeopd)
!$OMP DO reduction(+:eopd,neopd,aeopd) schedule(guided)
c
c     calculate the out-of-plane distance energy term
c
      do i = 1, nopdist
         ia = iopd(1,i)
         ib = iopd(2,i)
         ic = iopd(3,i)
         id = iopd(4,i)
         force = opdk(i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,ia,ib,ic,id,0,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
c
c     get the coordinates of the central and peripheral atoms
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
c     compute the out-of-plane distance for central atom
c
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
            xt = ybd*zcd - zbd*ycd
            yt = zbd*xcd - xbd*zcd
            zt = xbd*ycd - ybd*xcd
            rt2 = xt*xt + yt*yt + zt*zt
            dt2 = (xt*xad + yt*yad + zt*zad)**2 / rt2
            dt = sqrt(dt2)
            dt3 = dt2 * dt
            dt4 = dt2 * dt2
c
c     find the out-of-plane distance energy
c
            e = opdunit * force * dt2
     &             * (1.0d0+copd*dt+qopd*dt2+popd*dt3+sopd*dt4)
c
c     scale the interaction based on its group membership
c
            if (use_group)  e = e * fgrp
c
c     increment the total out-of-plane distance energy
c
            neopd = neopd + 1
            eopd = eopd + e
            aeopd(ia) = aeopd(ia) + e
c
c     print a message if the energy of this interaction is large
c
            huge = (e .gt. 5.0d0)
            if (debug .or. (verbose.and.huge)) then
               if (header) then
                  header = .false.
                  write (iout,20)
   20             format (/,' Individual Out-of-Plane Distance',
     &                       ' Interactions :',
     &                    //,' Type',25x,'Atom Names',18x,'Distance',
     &                       6x,'Energy',/)
               end if
               write (iout,30)  ia,name(ia),ib,name(ib),ic,
     &                          name(ic),id,name(id),sqrt(dt2),e
   30          format (' O-P-Dist',2x,4(i7,'-',a3),f11.4,f12.4)
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
