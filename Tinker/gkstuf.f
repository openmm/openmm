c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2006 by Michael Schnieders & Jay W. Ponder  ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module gkstuf  --  generalized Kirkwood solvation values  ##
c     ##                                                            ##
c     ################################################################
c
c
c     gkc      tuning parameter exponent in the f(GB) function
c     gkr      generalized Kirkwood cavity radii for atom types
c
c
      module gkstuf
      use sizes
      implicit none
      real*8 gkc
      real*8 gkr(maxtyp)
      save
      end
