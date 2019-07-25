c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine eangang  --  angle-angle cross term energy  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "eangang" calculates the angle-angle potential energy
c
c
      subroutine eangang
      use sizes
      use angang
      use angbnd
      use angpot
      use atoms
      use bound
      use energi
      use group
      use math
      use usage
      implicit none
      integer i,k,iangang
      integer ia,ib,ic,id,ie
      real*8 e,fgrp
      real*8 dt1,dt2
      real*8 angle,dot,cosine
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xie,yie,zie
      real*8 xab,yab,zab
      real*8 xcb,ycb,zcb
      real*8 xdb,ydb,zdb
      real*8 xeb,yeb,zeb
      real*8 rab2,rcb2
      real*8 rdb2,reb2
      logical proceed
c
c
c     zero out the angle-angle cross term energy
c
      eaa = 0.0d0
      if (nangang .eq. 0)  return
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(nangang,iaa,iang,
!$OMP& use,x,y,z,anat,kaa,aaunit,use_group,use_polymer)
!$OMP& shared(eaa)
!$OMP DO reduction(+:eaa) schedule(guided)
c
c     calculate the angle-angle interaction energy term
c
      do iangang = 1, nangang
         i = iaa(1,iangang)
         k = iaa(2,iangang)
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         id = iang(1,k)
         ie = iang(3,k)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,ia,ib,ic,id,ie,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or. use(ic)
     &                               .or. use(id) .or. use(ie))
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
            xid = x(id)
            yid = y(id)
            zid = z(id)
            xie = x(ie)
            yie = y(ie)
            zie = z(ie)
c
c     compute the values of the two bond angles
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
            xeb = xie - xib
            yeb = yie - yib
            zeb = zie - zib
            if (use_polymer) then
               call image (xab,yab,zab)
               call image (xcb,ycb,zcb)
               call image (xdb,ydb,zdb)
               call image (xeb,yeb,zeb)
            end if
            rab2 = xab*xab + yab*yab + zab*zab
            rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
            rdb2 = xdb*xdb + ydb*ydb + zdb*zdb
            reb2 = xeb*xeb + yeb*yeb + zeb*zeb
            if (rab2*rcb2*rdb2*reb2 .ne. 0.0d0) then
               dot = xab*xcb + yab*ycb + zab*zcb
               cosine = dot / sqrt(rab2*rcb2)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
               dt1 = angle - anat(i)
               dot = xdb*xeb + ydb*yeb + zdb*zeb
               cosine = dot / sqrt(rdb2*reb2)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
               dt2 = angle - anat(k)
c
c     get the angle-angle interaction energy
c
               e = aaunit * kaa(iangang) * dt1 * dt2
c
c     scale the interaction based on its group membership
c
               if (use_group)  e = e * fgrp
c
c     increment the total angle-angle energy
c
               eaa = eaa + e
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
