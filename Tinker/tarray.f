c
c
c     ################################################################
c     ## COPYRIGHT (C) 2013 by Xiao Zhu, Pengyu Ren & Jay W. Ponder ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  module tarray  --  store dipole-dipole matrix elements  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     ntpair     number of stored dipole-dipole matrix elements
c     tindex     index into stored dipole-dipole matrix values
c     tdipdip    stored dipole-dipole matrix element values
c
c
      module tarray
      implicit none
      integer ntpair
      integer, allocatable :: tindex(:,:)
      real*8, allocatable :: tdipdip(:,:)
      save
      end
