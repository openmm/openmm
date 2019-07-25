c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  module params  --  force field parameter file contents  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     maxprm    maximum number of lines in the parameter file
c
c     nprm      number of nonblank lines in the parameter file
c     prmline   contents of each individual parameter file line
c
c
      module params
      implicit none
      integer maxprm
      parameter (maxprm=25000)
      integer nprm
      character*240 prmline(maxprm)
      save
      end
