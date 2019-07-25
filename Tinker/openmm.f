c
c
c     ############################################################
c     ##                  COPYRIGHT (C) 2015                    ##
c     ##     by Mark Friedrichs, Lee-Ping Wang & Jay Ponder     ##
c     ##                  All Rights Reserved                   ##
c     ############################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  module openmm  --  OpenMM-related objects & variables  ##
c     ##                                                         ##
c     #############################################################
c
c
c     ommhandle   opaque handle pointing to OpenMM data structure
c     cudaDevice  string containing names of the CUDA GPU cards
c
c
      module openmm
      implicit none
      integer*8 ommhandle
      character*16 cudaDevice
      save
      end
