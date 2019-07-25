c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module resdue  --  amino acid & nucleotide residue names  ##
c     ##                                                            ##
c     ################################################################
c
c
c     maxamino  maximum number of amino acid residue types
c     maxnuc    maximum number of nucleic acid residue types
c
c     ntyp      biotypes for mid-chain peptide backbone N atoms
c     catyp     biotypes for mid-chain peptide backbone CA atoms
c     ctyp      biotypes for mid-chain peptide backbone C atoms
c     hntyp     biotypes for mid-chain peptide backbone HN atoms
c     otyp      biotypes for mid-chain peptide backbone O atoms
c     hatyp     biotypes for mid-chain peptide backbone HA atoms
c     cbtyp     biotypes for mid-chain peptide backbone CB atoms
c     nntyp     biotypes for N-terminal peptide backbone N atoms
c     cantyp    biotypes for N-terminal peptide backbone CA atoms
c     cntyp     biotypes for N-terminal peptide backbone C atoms
c     hnntyp    biotypes for N-terminal peptide backbone HN atoms
c     ontyp     biotypes for N-terminal peptide backbone O atoms
c     hantyp    biotypes for N-terminal peptide backbone HA atoms
c     nctyp     biotypes for C-terminal peptide backbone N atoms
c     cactyp    biotypes for C-terminal peptide backbone CA atoms
c     cctyp     biotypes for C-terminal peptide backbone C atoms
c     hnctyp    biotypes for C-terminal peptide backbone HN atoms
c     octyp     biotypes for C-terminal peptide backbone O atoms
c     hactyp    biotypes for C-terminal peptide backbone HA atoms
c     o5typ     biotypes for nucleotide backbone and sugar O5' atoms
c     c5typ     biotypes for nucleotide backbone and sugar C5' atoms
c     h51typ    biotypes for nucleotide backbone and sugar H5' atoms
c     h52typ    biotypes for nucleotide backbone and sugar H5'' atoms
c     c4typ     biotypes for nucleotide backbone and sugar C4' atoms
c     h4typ     biotypes for nucleotide backbone and sugar H4' atoms
c     o4typ     biotypes for nucleotide backbone and sugar O4' atoms
c     c1typ     biotypes for nucleotide backbone and sugar C1' atoms
c     h1typ     biotypes for nucleotide backbone and sugar H1' atoms
c     c3typ     biotypes for nucleotide backbone and sugar C3' atoms
c     h3typ     biotypes for nucleotide backbone and sugar H3' atoms
c     c2typ     biotypes for nucleotide backbone and sugar C2' atoms
c     h21typ    biotypes for nucleotide backbone and sugar H2' atoms
c     o2typ     biotypes for nucleotide backbone and sugar O2' atoms
c     h22typ    biotypes for nucleotide backbone and sugar H2'' atoms
c     o3typ     biotypes for nucleotide backbone and sugar O3' atoms
c     ptyp      biotypes for nucleotide backbone and sugar P atoms
c     optyp     biotypes for nucleotide backbone and sugar OP atoms
c     h5ttyp    biotypes for nucleotide backbone and sugar H5T atoms
c     h3ttyp    biotypes for nucleotide backbone and sugar H3T atoms
c     amino     three-letter abbreviations for amino acids types
c     nuclz     three-letter abbreviations for nucleic acids types
c     amino1    one-letter abbreviations for amino acids types
c     nuclz1    one-letter abbreviations for nucleic acids types
c
c
      module resdue
      implicit none
      integer maxamino
      integer maxnuc
      parameter (maxamino=38)
      parameter (maxnuc=12)
      integer ntyp(maxamino)
      integer catyp(maxamino)
      integer ctyp(maxamino)
      integer hntyp(maxamino)
      integer otyp(maxamino)
      integer hatyp(maxamino)
      integer cbtyp(maxamino)
      integer nntyp(maxamino)
      integer cantyp(maxamino)
      integer cntyp(maxamino)
      integer hnntyp(maxamino)
      integer ontyp(maxamino)
      integer hantyp(maxamino)
      integer nctyp(maxamino)
      integer cactyp(maxamino)
      integer cctyp(maxamino)
      integer hnctyp(maxamino)
      integer octyp(maxamino)
      integer hactyp(maxamino)
      integer o5typ(maxnuc)
      integer c5typ(maxnuc)
      integer h51typ(maxnuc)
      integer h52typ(maxnuc)
      integer c4typ(maxnuc)
      integer h4typ(maxnuc)
      integer o4typ(maxnuc)
      integer c1typ(maxnuc)
      integer h1typ(maxnuc)
      integer c3typ(maxnuc)
      integer h3typ(maxnuc)
      integer c2typ(maxnuc)
      integer h21typ(maxnuc)
      integer o2typ(maxnuc)
      integer h22typ(maxnuc)
      integer o3typ(maxnuc)
      integer ptyp(maxnuc)
      integer optyp(maxnuc)
      integer h5ttyp(maxnuc)
      integer h3ttyp(maxnuc)
      character*1 amino1(maxamino)
      character*1 nuclz1(maxnuc)
      character*3 amino(maxamino)
      character*3 nuclz(maxnuc)
      save
      end
