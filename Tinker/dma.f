c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2005  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module dma  --  QM spherical harmonic multipole moments  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     mp        atomic monopole charge values from DMA
c     dpx       atomic dipole moment x-component from DMA
c     dpy       atomic dipole moment y-component from DMA
c     dpz       atomic dipole moment z-component from DMA
c     q20       atomic Q20 quadrupole component from DMA (zz)
c     q21c      atomic Q21c quadrupole component from DMA (xz)
c     q21s      atomic Q21s quadrupole component from DMA (yz)
c     q22c      atomic Q22c quadrupole component from DMA (xx-yy)
c     q22s      atomic Q22s quadrupole component from DMA (xy)
c
c
      module dma
      implicit none
      real*8, allocatable :: mp(:)
      real*8, allocatable :: dpx(:)
      real*8, allocatable :: dpy(:)
      real*8, allocatable :: dpz(:)
      real*8, allocatable :: q20(:)
      real*8, allocatable :: q21c(:)
      real*8, allocatable :: q21s(:)
      real*8, allocatable :: q22c(:)
      real*8, allocatable :: q22s(:)
      save
      end
