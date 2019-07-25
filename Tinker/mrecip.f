c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2015  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  module mrecip  --  reciprocal PME for permanent multipoles  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     vmxx    scalar sum xx-component of virial due to multipoles
c     vmyy    scalar sum yy-component of virial due to multipoles
c     vmzz    scalar sum zz-component of virial due to multipoles
c     vmxy    scalar sum xy-component of virial due to multipoles
c     vmxz    scalar sum xz-component of virial due to multipoles
c     vmyz    scalar sum yz-component of virial due to multipoles
c     cmp     Cartesian permenent multipoles as polytensor vector
c     fmp     fractional permanent multipoles as polytensor vector
c     cphi    Cartesian permanent multipole potential and field
c     fphi    fractional permanent multipole potential and field
c
c
      module mrecip
      implicit none
      real*8 vmxx,vmyy,vmzz
      real*8 vmxy,vmxz,vmyz
      real*8, allocatable :: cmp(:,:)
      real*8, allocatable :: fmp(:,:)
      real*8, allocatable :: cphi(:,:)
      real*8, allocatable :: fphi(:,:)
      save
      end
