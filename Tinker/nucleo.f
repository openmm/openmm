c
c
c     #############################################################
c     ##                  COPYRIGHT (C) 1999 by                  ##
c     ##  Marina A. Vorobieva, Nina N. Sokolova & Jay W. Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module nucleo  --  parameters for nucleic acid structure  ##
c     ##                                                            ##
c     ################################################################
c
c
c     pucker    sugar pucker, either 2=2'-endo or 3=3'-endo
c     glyco     glycosidic torsional angle for each nucleotide
c     bkbone    phosphate backbone angles for each nucleotide
c     dblhlx    flag to mark system as nucleic acid double helix
c     deoxy     flag to mark deoxyribose or ribose sugar units
c     hlxform   helix form (A, B or Z) of polynucleotide strands
c
c
      module nucleo
      use sizes
      implicit none
      integer pucker(maxres)
      real*8 glyco(maxres)
      real*8 bkbone(6,maxres)
      logical dblhlx
      logical deoxy(maxres)
      character*1 hlxform
      save
      end
