c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  module rxnpot  --  reaction field functional form details  ##
c     ##                                                             ##
c     #################################################################
c
c
c     rfsize    radius of reaction field sphere centered at origin
c     rfbulkd   bulk dielectric constant of reaction field continuum
c     rfterms   number of terms to use in reaction field summation
c
c
      module rxnpot
      implicit none
      integer rfterms
      real*8 rfsize
      real*8 rfbulkd
      save
      end
