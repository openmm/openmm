c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  module kdipol  --  bond dipole forcefield parameters  ##
c     ##                                                        ##
c     ############################################################
c
c
c     maxnd    maximum number of bond dipole parameter entries
c     maxnd5   maximum number of 5-membered ring dipole entries
c     maxnd4   maximum number of 4-membered ring dipole entries
c     maxnd3   maximum number of 3-membered ring dipole entries
c
c     dpl      dipole moment parameters for bond dipoles
c     dpl5     dipole moment parameters for 5-ring dipoles
c     dpl4     dipole moment parameters for 4-ring dipoles
c     dpl3     dipole moment parameters for 3-ring dipoles
c     pos      dipole position parameters for bond dipoles
c     pos5     dipole position parameters for 5-ring dipoles
c     pos4     dipole position parameters for 4-ring dipoles
c     pos3     dipole position parameters for 3-ring dipoles
c     kd       string of atom classes for bond dipoles
c     kd5      string of atom classes for 5-ring dipoles
c     kd4      string of atom classes for 4-ring dipoles
c     kd3      string of atom classes for 3-ring dipoles
c
c
      module kdipol
      implicit none
      integer maxnd
      integer maxnd5
      integer maxnd4
      integer maxnd3
      parameter (maxnd=1000)
      parameter (maxnd5=500)
      parameter (maxnd4=500)
      parameter (maxnd3=500)
      real*8 dpl(maxnd)
      real*8 dpl5(maxnd5)
      real*8 dpl4(maxnd4)
      real*8 dpl3(maxnd3)
      real*8 pos(maxnd)
      real*8 pos5(maxnd5)
      real*8 pos4(maxnd4)
      real*8 pos3(maxnd3)
      character*8 kd(maxnd)
      character*8 kd5(maxnd5)
      character*8 kd4(maxnd4)
      character*8 kd3(maxnd3)
      save
      end
