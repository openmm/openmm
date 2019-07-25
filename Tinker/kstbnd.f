c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  module kstbnd  --  stretch-bend forcefield parameters  ##
c     ##                                                         ##
c     #############################################################
c
c
c     maxnsb   maximum number of stretch-bend parameter entries
c
c     stbn     force constant parameters for stretch-bend terms
c     ksb      string of atom classes for stretch-bend terms
c
c
      module kstbnd
      implicit none
      integer maxnsb
      parameter (maxnsb=2000)
      real*8 stbn(2,maxnsb)
      character*12 ksb(maxnsb)
      save
      end
