c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2009 by Chuanjie Wu and Jay William Ponder  ##
c     ##                    All Rights Reserved                     ##
c     ################################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##   module qmstuf  --  quantum data from Gaussian calculation  ##
c     ##                                                              ##
c     ##################################################################
c
c   
c     ngatom   number of atoms in the QM data file
c     egau     quantum mechanical (QM) total energy (kcal/mole)
c     gx       x-coordinate of each atom in the QM data file
c     gy       y-coordinate of each atom in the QM data file
c     gz       z-coordinate of each atom in the QM data file
c     gfreq    calculated vibrational frequencies from QM data
c     gforce   force components on each atom from QM data
c     gh       Hessian matrix elements from QM data
c
c
      module qmstuf
      implicit none
      integer ngatom
      real*8 egau
      real*8, allocatable :: gx(:)
      real*8, allocatable :: gy(:)
      real*8, allocatable :: gz(:)
      real*8, allocatable :: gfreq(:)
      real*8, allocatable :: gforce(:,:)
      real*8, allocatable :: gh(:)
      save
      end
