c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  module output  --  coordinate file format output controls  ##
c     ##                                                             ##
c     #################################################################
c
c
c     archive    logical flag to save structures in an archive
c     noversion  logical flag governing use of filename versions
c     overwrite  logical flag to overwrite intermediate files inplace
c     cyclesave  logical flag to mark use of numbered cycle files
c     coordtype  selects Cartesian, internal, rigid body or none
c
c
      module output
      implicit none
      logical archive
      logical noversion
      logical overwrite
      logical cyclesave
      character*9 coordtype
      save
      end
