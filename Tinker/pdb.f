c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  module pdb  --  Protein Data Bank structure definition  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     npdb      number of atoms stored in Protein Data Bank format
c     nres      number of residues stored in Protein Data Bank format
c     resnum    number of the residue to which each atom belongs
c     resatm    number of first and last atom in each residue
c     npdb12    number of atoms directly bonded to each CONECT atom
c     ipdb12    atom numbers of atoms connected to each CONECT atom
c     pdblist   list of the Protein Data Bank atom number of each atom
c     xpdb      x-coordinate of each atom stored in PDB format
c     ypdb      y-coordinate of each atom stored in PDB format
c     zpdb      z-coordinate of each atom stored in PDB format
c     altsym    string with PDB alternate locations to be included
c     pdbres    Protein Data Bank residue name assigned to each atom
c     pdbatm    Protein Data Bank atom name assigned to each atom
c     pdbtyp    Protein Data Bank record type assigned to each atom
c     chnsym    string with PDB chain identifiers to be included
c     instyp    string with PDB insertion records to be included
c
c
      module pdb
      implicit none
      integer npdb,nres
      integer, allocatable :: resnum(:)
      integer, allocatable :: resatm(:,:)
      integer, allocatable :: npdb12(:)
      integer, allocatable :: ipdb12(:,:)
      integer, allocatable :: pdblist(:)
      real*8, allocatable :: xpdb(:)
      real*8, allocatable :: ypdb(:)
      real*8, allocatable :: zpdb(:)
      character*1 altsym
      character*3, allocatable :: pdbres(:)
      character*4, allocatable :: pdbatm(:)
      character*6, allocatable :: pdbtyp(:)
      character*20 chnsym
      character*20 instyp
      save
      end
