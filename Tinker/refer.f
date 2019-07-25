c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  module refer  --  reference atomic coordinate storage  ##
c     ##                                                         ##
c     #############################################################
c
c
c     nref        total number of atoms in each reference system
c     refltitle   length in characters of each reference title line
c     refleng     length in characters of each reference filename
c     reftyp      atom types of the atoms in each reference system
c     n12ref      number of atoms bonded to each reference atom
c     i12ref      atom numbers of atoms 1-2 connected to each atom
c     xboxref     reference a-axis length of periodic box
c     yboxref     reference b-axis length of periodic box
c     zboxref     reference c-axis length of periodic box
c     alpharef    reference angle between b- and c-axes of box
c     betaref     reference angle between a- and c-axes of box
c     gammaref    reference angle between a- and b-axes of box
c     xref        reference x-coordinates for atoms in each system
c     yref        reference y-coordinates for atoms in each system
c     zref        reference z-coordinates for atoms in each system
c     refnam      atom names of the atoms in each reference system
c     reffile     base filename for each reference system
c     reftitle    title used to describe each reference system
c
c
      module refer
      use sizes
      implicit none
      integer nref(maxref)
      integer refltitle(maxref)
      integer refleng(maxref)
      integer, allocatable :: reftyp(:,:)
      integer, allocatable :: n12ref(:,:)
      integer, allocatable :: i12ref(:,:,:)
      real*8 xboxref(maxref)
      real*8 yboxref(maxref)
      real*8 zboxref(maxref)
      real*8 alpharef(maxref)
      real*8 betaref(maxref)
      real*8 gammaref(maxref)
      real*8, allocatable :: xref(:,:)
      real*8, allocatable :: yref(:,:)
      real*8, allocatable :: zref(:,:)
      character*3, allocatable :: refnam(:,:)
      character*240 reffile(maxref)
      character*240 reftitle(maxref)
      save
      end
