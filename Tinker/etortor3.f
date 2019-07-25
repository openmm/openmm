c
c
c     #############################################################
c     ##  COPYRIGHT (C) 2003 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine etortor3  --  torsion-torsion energy & analysis  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "etortor3" calculates the torsion-torsion potential energy;
c     also partitions the energy terms among the atoms
c
c
      subroutine etortor3
      use sizes
      use action
      use analyz
      use atoms
      use bitor
      use bound
      use energi
      use group
      use inform
      use iounit
      use ktrtor
      use math
      use torpot
      use tortor
      use usage
      implicit none
      integer i,k,itortor
      integer pos1,pos2
      integer ia,ib,ic,id,ie
      integer nlo,nhi,nt
      integer xlo,ylo
      real*8 e,fgrp,sign
      real*8 angle1,angle2
      real*8 value1,value2
      real*8 cosine1,cosine2
      real*8 xt,yt,zt,rt2
      real*8 xu,yu,zu,ru2
      real*8 xv,yv,zv,rv2
      real*8 xtu,ytu,ztu,rtru
      real*8 xuv,yuv,zuv,rurv
      real*8 xh,yh,x1l,x1u
      real*8 y1l,y1u
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xie,yie,zie
      real*8 xba,yba,zba
      real*8 xdc,ydc,zdc
      real*8 xcb,ycb,zcb
      real*8 xed,yed,zed
      real*8 ftt(4),ft12(4)
      real*8 ft1(4),ft2(4)
      logical proceed
      logical header,huge
c
c
c     zero out the torsion-torsion energy and partitioning terms
c
      nett = 0
      ett = 0.0d0
      do i = 1, n
         aett(i) = 0.0d0
      end do
      if (ntortor .eq. 0)  return
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. ntortor.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Torsion-Torsion Interactions :',
     &           //,' Type',17x,'Atom Numbers',16x,'Angle1',
     &              4x,'Angle2',6x,'Energy',/)
      end if
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(ntortor,itt,ibitor,
!$OMP& use,x,y,z,tnx,ttx,tny,tty,tbf,tbx,tby,tbxy,ttorunit,
!$OMP& use_group,use_polymer,verbose,debug,header,iout)
!$OMP& shared(ett,nett,aett)
!$OMP DO reduction(+:ett,nett,aett) schedule(guided)
c
c     calculate the torsion-torsion interaction energy term
c
      do itortor = 1, ntortor
         i = itt(1,itortor)
         k = itt(2,itortor)
         if (itt(3,itortor) .eq. 1) then
            ia = ibitor(1,i)
            ib = ibitor(2,i)
            ic = ibitor(3,i)
            id = ibitor(4,i)
            ie = ibitor(5,i)
         else
            ia = ibitor(5,i)
            ib = ibitor(4,i)
            ic = ibitor(3,i)
            id = ibitor(2,i)
            ie = ibitor(1,i)
         end if
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,ia,ib,ic,id,ie,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or. use(ic)
     &                               .or. use(id) .or. use(ie))
c
c     compute the value of the torsional angles
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
            xie = x(ie)
            yie = y(ie)
            zie = z(ie)
            xba = xib - xia
            yba = yib - yia
            zba = zib - zia
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdc = xid - xic
            ydc = yid - yic
            zdc = zid - zic
            xed = xie - xid
            yed = yie - yid
            zed = zie - zid
            if (use_polymer) then
               call image (xba,yba,zba)
               call image (xcb,ycb,zcb)
               call image (xdc,ydc,zdc)
               call image (xed,yed,zed)
            end if
            xt = yba*zcb - ycb*zba
            yt = zba*xcb - zcb*xba
            zt = xba*ycb - xcb*yba
            xu = ycb*zdc - ydc*zcb
            yu = zcb*xdc - zdc*xcb
            zu = xcb*ydc - xdc*ycb
            xtu = yt*zu - yu*zt
            ytu = zt*xu - zu*xt
            ztu = xt*yu - xu*yt
            rt2 = xt*xt + yt*yt + zt*zt
            ru2 = xu*xu + yu*yu + zu*zu
            rtru = sqrt(rt2 * ru2)
            xv = ydc*zed - yed*zdc
            yv = zdc*xed - zed*xdc
            zv = xdc*yed - xed*ydc
            xuv = yu*zv - yv*zu
            yuv = zu*xv - zv*xu
            zuv = xu*yv - xv*yu
            rv2 = xv*xv + yv*yv + zv*zv
            rurv = sqrt(ru2 * rv2)
            if (rtru.ne.0.0d0 .and. rurv.ne.0.0d0) then
               cosine1 = (xt*xu + yt*yu + zt*zu) / rtru
               cosine1 = min(1.0d0,max(-1.0d0,cosine1))
               angle1 = radian * acos(cosine1)
               sign = xba*xu + yba*yu + zba*zu
               if (sign .lt. 0.0d0)  angle1 = -angle1
               value1 = angle1
               cosine2 = (xu*xv + yu*yv + zu*zv) / rurv
               cosine2 = min(1.0d0,max(-1.0d0,cosine2))
               angle2 = radian * acos(cosine2)
               sign = xcb*xv + ycb*yv + zcb*zv
               if (sign .lt. 0.0d0)  angle2 = -angle2
               value2 = angle2
c
c     check for inverted chirality at the central atom
c
               call chkttor (ib,ic,id,sign,value1,value2)
c
c     use bicubic interpolation to compute spline values
c
               nlo = 1
               nhi = tnx(k)
               do while (nhi-nlo .gt. 1)
                  nt = (nhi+nlo) / 2
                  if (ttx(nt,k) .gt. value1) then
                     nhi = nt
                  else
                     nlo = nt
                  end if
               end do
               xlo = nlo
               nlo = 1
               nhi = tny(k)
               do while (nhi-nlo .gt. 1)
                  nt = (nhi + nlo)/2
                  if (tty(nt,k) .gt. value2) then
                     nhi = nt
                  else
                     nlo = nt
                  end if
               end do
               ylo = nlo
               xh = ttx(xlo+1,k) - ttx(xlo,k)
               yh = tty(ylo+1,k) - tty(ylo,k)
               x1l = ttx(xlo,k)
               x1u = ttx(xlo+1,k)
               y1l = tty(ylo,k)
               y1u = tty(ylo+1,k)
               pos2 = ylo*tnx(k) + xlo
               pos1 = pos2 - tnx(k)
               ftt(1) = tbf(pos1,k)
               ftt(2) = tbf(pos1+1,k)
               ftt(3) = tbf(pos2+1,k)
               ftt(4) = tbf(pos2,k)
               ft1(1) = tbx(pos1,k)
               ft1(2) = tbx(pos1+1,k)
               ft1(3) = tbx(pos2+1,k)
               ft1(4) = tbx(pos2,k)
               ft2(1) = tby(pos1,k)
               ft2(2) = tby(pos1+1,k)
               ft2(3) = tby(pos2+1,k)
               ft2(4) = tby(pos2,k)
               ft12(1) = tbxy(pos1,k)
               ft12(2) = tbxy(pos1+1,k)
               ft12(3) = tbxy(pos2+1,k)
               ft12(4) = tbxy(pos2,k)
               call bcuint (ftt,ft1,ft2,ft12,x1l,x1u,
     &                      y1l,y1u,value1,value2,e)
               e = ttorunit * e
c
c     scale the interaction based on its group membership
c
               if (use_group)  e = e * fgrp
c
c     increment the total torsion-torsion energy
c
               nett = nett + 1
               ett = ett + e
               aett(ib) = aett(ib) + e/3.0d0
               aett(ic) = aett(ic) + e/3.0d0
               aett(id) = aett(id) + e/3.0d0
c
c     print a message if the energy of this interaction is large
c
               huge = (e .gt. 5.0d0)
               if (debug .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,20)
   20                format (/,' Individual Torsion-Torsion',
     &                          ' Interactions :',
     &                       //,' Type',17x,'Atom Numbers',16x,'Angle1',
     &                          4x,'Angle2',6x,'Energy',/)
                  end if
                  write (iout,30)  ia,ib,ic,id,ie,angle1,angle2,e
   30             format (' TorTor',2x,5i7,2x,2f10.4,f12.4)
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
