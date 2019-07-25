c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  module atoms  --  number, position and type of atoms  ##
c     ##                                                        ##
c     ############################################################
c
c
c     n       total number of atoms in the current system
c     type    atom type number for each atom in the system
c     x       current x-coordinate for each atom in the system
c     y       current y-coordinate for each atom in the system
c     z       current z-coordinate for each atom in the system
c
c
      module atoms
      use sizes
      implicit none
      integer n
      integer type(maxatm)
      real*8 x(maxatm)
      real*8 y(maxatm)
      real*8 z(maxatm)
      save
      end
