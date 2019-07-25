c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  module sequen  --  sequence information for biopolymer  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     nseq     total number of residues in biopolymer sequences
c     nchain   number of separate biopolymer sequence chains
c     ichain   first and last residue in each biopolymer chain
c     seqtyp   residue type for each residue in the sequence
c     seq      three-letter code for each residue in the sequence
c     chnnam   one-letter identifier for each sequence chain
c     chntyp   contents of each chain (GENERIC, PEPTIDE or NUCLEIC)
c
c
      module sequen
      use sizes
      implicit none
      integer nseq
      integer nchain
      integer ichain(2,maxres)
      integer seqtyp(maxres)
      character*1 chnnam(maxres)
      character*3 seq(maxres)
      character*7 chntyp(maxres)
      save
      end
