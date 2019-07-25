c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2010  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module chunks  --  PME grid spatial decomposition values  ##
c     ##                                                            ##
c     ################################################################
c
c
c     nchunk     total number of spatial regions for PME grid
c     nchk1      number of spatial regions along the a-axis
c     nchk2      number of spatial regions along the b-axis
c     nchk3      number of spatial regions along the c-axis
c     ngrd1      number of grid points per region along a-axis
c     ngrd2      number of grid points per region along b-axis
c     ngrd3      number of grid points per region along c-axis
c     nlpts      PME grid points to the left of center point
c     nrpts      PME grid points to the right of center point
c     grdoff     offset for index into B-spline coefficients
c     pmetable   PME grid spatial regions involved for each site
c
c
      module chunks
      implicit none
      integer nchunk
      integer nchk1,nchk2,nchk3
      integer ngrd1,ngrd2,ngrd3
      integer nlpts,nrpts,grdoff
      integer, allocatable :: pmetable(:,:)
      save
      end
