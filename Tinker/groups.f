c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine groups  --  group membership of set of atoms  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "groups" tests a set of atoms to see if all are members of a
c     single atom group or a pair of atom groups; if so, then the
c     correct intra- or intergroup weight is assigned
c
c     note the default group-based interaction weight is 1.0; only
c     interactions involving two or fewer groups can be scaled
c
c
      subroutine groups (proceed,weigh,ia,ib,ic,id,ie,ig)
      use group
      implicit none
      integer ia,ib,ic
      integer id,ie,ig
      integer iga,igb,igc
      integer igd,ige,igg
      integer nset
      integer gmax,gmin
      real*8 weigh
      logical proceed
c
c
c     determine the number of atoms in the set to be compared
c
      nset = 0
      weigh = 1.0d0
      if (ig .ne. 0) then
         nset = 6
      else if (ie .ne. 0) then
         nset = 5
      else if (id .ne. 0) then
         nset = 4
      else if (ic .ne. 0) then
         nset = 3
      else if (ib .ne. 0) then
         nset = 2
      else if (ia .ne. 0) then
         nset = 1
      end if
c
c     check group membership for a set containing one atom
c
      if (nset .eq. 1) then
         iga = grplist(ia)
         weigh = wgrp(iga,iga)
c
c     check group membership for a set containing two atoms
c
      else if (nset .eq. 2) then
         iga = grplist(ia)
         igb = grplist(ib)
         weigh = wgrp(iga,igb)
c
c     check group membership for a set containing three atoms
c
      else if (nset .eq. 3) then
         iga = grplist(ia)
         igb = grplist(ib)
         igc = grplist(ic)
         if (iga.eq.igb .or. igb.eq.igc) then
            weigh = wgrp(iga,igc)
         else if (iga .eq. igc) then
            weigh = wgrp(iga,igb)
         end if
c
c     check group membership for a set containing four atoms
c
      else if (nset .eq. 4) then
         iga = grplist(ia)
         igb = grplist(ib)
         igc = grplist(ic)
         igd = grplist(id)
         gmin = min(iga,igb,igc,igd)
         gmax = max(iga,igb,igc,igd)
         if ((iga.eq.gmin .or. iga.eq.gmax) .and.
     &       (igb.eq.gmin .or. igb.eq.gmax) .and.
     &       (igc.eq.gmin .or. igc.eq.gmax) .and.
     &       (igd.eq.gmin .or. igd.eq.gmax))  weigh = wgrp(gmin,gmax)
c
c     check group membership for a set containing five atoms
c
      else if (nset .eq. 5) then
         iga = grplist(ia)
         igb = grplist(ib)
         igc = grplist(ic)
         igd = grplist(id)
         ige = grplist(ie)
         gmin = min(iga,igb,igc,igd,ige)
         gmax = max(iga,igb,igc,igd,ige)
         if ((iga.eq.gmin .or. iga.eq.gmax) .and.
     &       (igb.eq.gmin .or. igb.eq.gmax) .and.
     &       (igc.eq.gmin .or. igc.eq.gmax) .and.
     &       (igd.eq.gmin .or. igd.eq.gmax) .and.
     &       (ige.eq.gmin .or. ige.eq.gmax))  weigh = wgrp(gmin,gmax)
c
c     check group membership for a set containing five atoms
c
      else if (nset .eq. 6) then
         iga = grplist(ia)
         igb = grplist(ib)
         igc = grplist(ic)
         igd = grplist(id)
         ige = grplist(ie)
         igg = grplist(ig)
         gmin = min(iga,igb,igc,igd,ige,igg)
         gmax = max(iga,igb,igc,igd,ige,igg)
         if ((iga.eq.gmin .or. iga.eq.gmax) .and.
     &       (igb.eq.gmin .or. igb.eq.gmax) .and.
     &       (igc.eq.gmin .or. igc.eq.gmax) .and.
     &       (igd.eq.gmin .or. igd.eq.gmax) .and.
     &       (ige.eq.gmin .or. ige.eq.gmax) .and.
     &       (igg.eq.gmin .or. igg.eq.gmax))  weigh = wgrp(gmin,gmax)
      end if
c
c     interaction will be used if its group has nonzero weight
c
      if (weigh .eq. 0.0d0) then
         proceed = .false.
      else
         proceed = .true.
      end if
      return
      end
