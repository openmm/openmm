c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module pme  --  values for particle mesh Ewald summation  ##
c     ##                                                            ##
c     ################################################################
c
c
c     nfft1      number of grid points along the a-axis direction
c     nfft2      number of grid points along the b-axis direction
c     nfft3      number of grid points along the c-axis direction
c     bsorder    order of the PME B-spline approximation
c     igrid      initial Ewald charge grid values for B-spline
c     bsmod1     B-spline moduli along the a-axis direction
c     bsmod2     B-spline moduli along the b-axis direction
c     bsmod3     B-spline moduli along the c-axis direction
c     bsbuild    B-spline derivative coefficient temporary storage
c     thetai1    B-spline coefficients along the a-axis
c     thetai2    B-spline coefficients along the b-axis
c     thetai3    B-spline coefficients along the c-axis
c     qgrid      values on the particle mesh Ewald charge grid
c     qfac       prefactors for particle mesh Ewald charge grid
c
c
      module pme
      implicit none
      integer nfft1
      integer nfft2
      integer nfft3
      integer bsorder
      integer, allocatable :: igrid(:,:)
      real*8, allocatable :: bsmod1(:)
      real*8, allocatable :: bsmod2(:)
      real*8, allocatable :: bsmod3(:)
      real*8, allocatable :: bsbuild(:,:)
      real*8, allocatable :: thetai1(:,:,:)
      real*8, allocatable :: thetai2(:,:,:)
      real*8, allocatable :: thetai3(:,:,:)
      real*8, allocatable :: qgrid(:,:,:,:)
      real*8, allocatable :: qfac(:,:,:)
      save
      end
