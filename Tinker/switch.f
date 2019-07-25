c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine switch  --  get switching function coefficients  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "switch" sets the coeffcients used by the fifth and seventh
c     order polynomial switching functions for spherical cutoffs
c
c
      subroutine switch (mode)
      use limits
      use nonpol
      use shunt
      implicit none
      real*8 denom,term
      real*8 off3,off4,off5
      real*8 off6,off7
      real*8 cut3,cut4,cut5
      real*8 cut6,cut7
      character*6 mode
c
c
c     get the switching window for the current potential type
c
      if (mode(1:3) .eq. 'VDW') then
         off = vdwcut
         cut = vdwtaper
      else if (mode(1:2) .eq. 'CT') then
         off = ctcut
         cut = cttaper
      else if (mode(1:6) .eq. 'CHARGE') then
         off = chgcut
         cut = chgtaper
      else if (mode(1:6) .eq. 'CHGDPL') then
         off = sqrt(chgcut*dplcut)
         cut = sqrt(chgtaper*dpltaper)
      else if (mode(1:6) .eq. 'DIPOLE') then
         off = dplcut
         cut = dpltaper
      else if (mode(1:5) .eq. 'MPOLE') then
         off = mpolecut
         cut = mpoletaper
      else if (mode(1:5) .eq. 'EWALD') then
         off = ewaldcut
         cut = ewaldcut
      else if (mode(1:6) .eq. 'USOLVE') then
         off = usolvcut
         cut = usolvcut
      else if (mode(1:3) .eq. 'GKV') then
         off = spoff
         cut = spcut
      else if (mode(1:4) .eq. 'GKSA') then
         off = stcut
         cut = stoff
      else
         off = min(vdwcut,ctcut,chgcut,dplcut,mpolecut)
         cut = min(vdwtaper,cttaper,chgtaper,dpltaper,mpoletaper) 
      end if
c
c     test for replicate periodic boundaries at this cutoff
c
      call replica (off)
c
c     set switching coefficients to zero for truncation cutoffs
c
      c0 = 0.0d0
      c1 = 0.0d0
      c2 = 0.0d0
      c3 = 0.0d0
      c4 = 0.0d0
      c5 = 0.0d0
      f0 = 0.0d0
      f1 = 0.0d0
      f2 = 0.0d0
      f3 = 0.0d0
      f4 = 0.0d0
      f5 = 0.0d0
      f6 = 0.0d0
      f7 = 0.0d0
c
c     store the powers of the switching window cutoffs
c
      off2 = off * off
      off3 = off2 * off
      off4 = off2 * off2
      off5 = off2 * off3
      off6 = off3 * off3
      off7 = off3 * off4
      cut2 = cut * cut
      cut3 = cut2 * cut
      cut4 = cut2 * cut2
      cut5 = cut2 * cut3
      cut6 = cut3 * cut3
      cut6 = cut3 * cut3
      cut7 = cut3 * cut4
c
c     get 5th degree multiplicative switching function coefficients
c
      if (cut .lt. off) then
         denom = (off-cut)**5
         c0 = off*off2 * (off2-5.0d0*off*cut+10.0d0*cut2) / denom
         c1 = -30.0d0 * off2*cut2 / denom
         c2 = 30.0d0 * (off2*cut+off*cut2) / denom
         c3 = -10.0d0 * (off2+4.0d0*off*cut+cut2) / denom
         c4 = 15.0d0 * (off+cut) / denom
         c5 = -6.0d0 / denom
      end if
c
c     get 7th degree additive switching function coefficients
c
      if (cut.lt.off .and. mode(1:6).eq.'CHARGE') then
         term = 9.3d0 * cut*off / (off-cut)
         denom = cut7 - 7.0d0*cut6*off + 21.0d0*cut5*off2
     &              - 35.0d0*cut4*off3 + 35.0d0*cut3*off4
     &              - 21.0d0*cut2*off5 + 7.0d0*cut*off6 - off7
         denom = term * denom
         f0 = cut3*off3 * (-39.0d0*cut+64.0d0*off) / denom
         f1 = cut2*off2
     &           * (117.0d0*cut2-100.0d0*cut*off-192.0d0*off2) / denom
         f2 = cut*off * (-117.0d0*cut3-84.0d0*cut2*off
     &                   +534.0d0*cut*off2+192.0d0*off3) / denom
         f3 = (39.0d0*cut4+212.0d0*cut3*off-450.0d0*cut2*off2
     &            -612.0d0*cut*off3-64.0d0*off4) / denom
         f4 = (-92.0d0*cut3+66.0d0*cut2*off
     &            +684.0d0*cut*off2+217.0d0*off3) / denom
         f5 = (42.0d0*cut2-300.0d0*cut*off-267.0d0*off2) / denom
         f6 = (36.0d0*cut+139.0d0*off) / denom
         f7 = -25.0d0 / denom
      end if
      return
      end
