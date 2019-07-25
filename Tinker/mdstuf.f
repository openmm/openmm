c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2000  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  module mdstuf  --  molecular dynamics trajectory controls  ##
c     ##                                                             ##
c     #################################################################
c
c
c     nfree       total number of degrees of freedom for a system
c     irest       steps between removal of COM motion (0=no removal)
c     bmnmix      mixing coefficient for use with Beeman integrator
c     dorest      logical flag to remove center of mass motion
c     velsave     logical flag to save velocity vector components
c     frcsave     logical flag to save force vector components
c     uindsave    logical flag to save induced atomic dipoles
c     integrate   type of molecular dynamics integration algorithm
c
c
      module mdstuf
      implicit none
      integer nfree
      integer irest
      integer bmnmix
      logical dorest
      logical velsave
      logical frcsave
      logical uindsave
      character*11 integrate
      save
      end
