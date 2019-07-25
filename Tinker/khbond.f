c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1998  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module khbond  --  H-bonding term forcefield parametersb ##
c     ##                                                           ##
c     ###############################################################
c
c
c     maxnhb   maximum number of hydrogen bonding pair entries
c
c     radhb    radius parameter for hydrogen bonding pairs
c     epshb    well depth parameter for hydrogen bonding pairs
c     khb      string of atom types for hydrogen bonding pairs
c
c
      module khbond
      implicit none
      integer maxnhb
      parameter (maxnhb=500)
      real*8 radhb(maxnhb)
      real*8 epshb(maxnhb)
      character*8 khb(maxnhb)
      save
      end
