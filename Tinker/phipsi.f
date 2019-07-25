c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module phipsi  --  phi-psi-omega-chi angles for protein  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     chiral   chirality of each amino acid residue (1=L, -1=D)
c     disulf   residue joined to each residue via a disulfide link
c     phi      value of the phi angle for each amino acid residue
c     psi      value of the psi angle for each amino acid residue
c     omega    value of the omega angle for each amino acid residue
c     chi      values of the chi angles for each amino acid residue
c
c
      module phipsi
      use sizes
      implicit none
      integer chiral(maxres)
      integer disulf(maxres)
      real*8 phi(maxres)
      real*8 psi(maxres)
      real*8 omega(maxres)
      real*8 chi(4,maxres)
      save
      end
