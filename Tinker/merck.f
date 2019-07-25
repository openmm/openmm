c
c
c     ##############################################################
c     ##  COPYRIGHT (C) 2007 by Nicolas Staelens & Jay W. Ponder  ##
c     ##                   All Rights Reserved                    ##
c     ##############################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  module merck  --  MMFF-specific force field parameters  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     nlignes    number of atom pairs having MMFF Bond Type 1
c     bt_1       atom pairs having MMFF Bond Type 1
c     eqclass    table of atom class equivalencies used to find
c                default parameters if explicit values are missing
c                (see J. Comput. Chem., 17, 490-519, '95, Table IV)
c     crd        number of attached neighbors    |
c     val        valency value                   |  see T. A. Halgren,
c     pilp       if 0, no lone pair              |  J. Comput. Chem.,
c                if 1, one or more lone pair(s)  |  17, 616-645 (1995)
c     mltb       multibond indicator             |
c     arom       aromaticity indicator           |
c     lin        linearity indicator             |
c     sbmb       single- vs multiple-bond flag   |
c     mmffarom   aromatic rings parameters
c     mmffaromc  cationic aromatic rings parameters
c     mmffaroma  anionic aromatic rings parameters
c
c
      module merck
      use sizes
      implicit none
      integer nlignes
      integer bt_1(500,2)
      integer eqclass(500,5)
      integer crd(100)
      integer val(100)
      integer pilp(100)
      integer mltb(100)
      integer arom(100)
      integer lin(100)
      integer sbmb(100)
      integer mmffarom(maxtyp,6)
      integer mmffaromc(maxtyp,6)
      integer mmffaroma(maxtyp,6)
c
c
c     rad0      covalent atomic radius for empirical bond rules
c     paulel    Pauling electronegativities for empirical bond rules
c     r0ref     reference bond length for empirical bond rules
c     kbref     reference force constant for empirical bond rules
c     mmff_kb   bond force constant for pairs of atom classes
c     mmff_kb1  bond force constant for class pairs with Bond Type 1
c     mmff_b0   bond length value for pairs of atom classes
c     mmff_b1   bond length value for class pairs with Bond Type 1
c
c
      real*8 rad0(100)
      real*8 paulel(100)
      real*8 r0ref(100,100)
      real*8 kbref(100,100)
      real*8 mmff_kb(100,100)
      real*8 mmff_kb1(100,100)
      real*8 mmff_b0(100,100)
      real*8 mmff_b1(100,100)
c
c
c     mmff_ka     angle force constant for triples of atom classes
c     mmff_ka1    angle force constant with one bond of Type 1
c     mmff_ka2    angle force constant with both bonds of Type 1
c     mmff_ka3    angle force constant for 3-membered ring
c     mmff_ka4    angle force constant for 4-membered ring
c     mmff_ka5    angle force constant for 3-ring and one Bond Type 1
c     mmff_ka6    angle force constant for 3-ring and both Bond Type 1
c     mmff_ka7    angle force constant for 4-ring and one Bond Type 1
c     mmff_ka8    angle force constant for 4-ring and both Bond Type 1
c     mmff_ang0   ideal bond angle for triples of atom classes
c     mmff_ang1   ideal bond angle with one bond of Type 1
c     mmff_ang2   ideal bond angle with both bonds of Type 1
c     mmff_ang3   ideal bond angle for 3-membered ring
c     mmff_ang4   ideal bond angle for 4-membered ring
c     mmff_ang5   ideal bond angle for 3-ring and one Bond Type 1
c     mmff_ang6   ideal bond angle for 3-ring and both Bond Type 1
c     mmff_ang7   ideal bond angle for 4-ring and one Bond Type 1
c     mmff_ang8   ideal bond angle for 4-ring and both Bond Type 1
c
c
      real*8, allocatable :: mmff_ka(:,:,:)
      real*8, allocatable :: mmff_ka1(:,:,:)
      real*8, allocatable :: mmff_ka2(:,:,:)
      real*8, allocatable :: mmff_ka3(:,:,:)
      real*8, allocatable :: mmff_ka4(:,:,:)
      real*8, allocatable :: mmff_ka5(:,:,:)
      real*8, allocatable :: mmff_ka6(:,:,:)
      real*8, allocatable :: mmff_ka7(:,:,:)
      real*8, allocatable :: mmff_ka8(:,:,:)
      real*8, allocatable :: mmff_ang0(:,:,:)
      real*8, allocatable :: mmff_ang1(:,:,:)
      real*8, allocatable :: mmff_ang2(:,:,:)
      real*8, allocatable :: mmff_ang3(:,:,:)
      real*8, allocatable :: mmff_ang4(:,:,:)
      real*8, allocatable :: mmff_ang5(:,:,:)
      real*8, allocatable :: mmff_ang6(:,:,:)
      real*8, allocatable :: mmff_ang7(:,:,:)
      real*8, allocatable :: mmff_ang8(:,:,:)
c
c
c     Stretch-Bend Type 0
c     stbn_abc     stretch-bend parameters for A-B-C atom classes
c     stbn_cba     stretch-bend parameters for C-B-A atom classes
c     Stretch-Bend Type 1  (A-B is Bond Type 1)
c     stbn_abc1    stretch-bend parameters for A-B-C atom classes
c     stbn_cba1    stretch-bend parameters for C-B-A atom classes
c     Stretch-Bend Type 2  (B-C is Bond Type 1) 
c     stbn_abc2    stretch-bend parameters for A-B-C atom classes
c     stbn_cba2    stretch-bend parameters for C-B-A atom classes
c     Stretch-Bend Type = 3  (A-B and B-C are Bond Type 1) 
c     stbn_abc3    stretch-bend parameters for A-B-C atom classes
c     stbn_cba3    stretch-bend parameters for C-B-A atom classes
c     Stretch-Bend Type 4  (both Bond Types 0, 4-membered ring)
c     stbn_abc4    stretch-bend parameters for A-B-C atom classes
c     stbn_cba4    stretch-bend parameters for C-B-A atom classes
c     Stretch-Bend Type 5  (both Bond Types 0, 3-membered ring)
c     stbn_abc5    stretch-bend parameters for A-B-C atom classes
c     stbn_cba5    stretch-bend parameters for C-B-A atom classes
c     Stretch-Bend Type 6  (A-B is Bond Type 1, 3-membered ring)
c     stbn_abc6    stretch-bend parameters for A-B-C atom classes
c     stbn_cba6    stretch-bend parameters for C-B-A atom classes
c     Stretch-Bend Type 7  (B-C is Bond Type 1, 3-membered ring)
c     stbn_abc7    stretch-bend parameters for A-B-C atom classes
c     stbn_cba7    stretch-bend parameters for C-B-A atom classes
c     Stretch-Bend Type 8  (both Bond Types 1, 3-membered ring)
c     stbn_abc8    stretch-bend parameters for A-B-C atom classes
c     stbn_cba8    stretch-bend parameters for C-B-A atom classes
c     Stretch-Bend Type 9  (A-B is Bond Type 1, 4-membered ring)
c     stbn_abc9    stretch-bend parameters for A-B-C atom classes
c     stbn_cba9    stretch-bend parameters for C-B-A atom classes
c     Stretch-Bend Type 10  (B-C is Bond Type 1, 4-membered ring)
c     stbn_abc10   stretch-bend parameters for A-B-C atom classes
c     stbn_cba10   stretch-bend parameters for C-B-A atom classes
c     Stretch-Bend Type 11  (both Bond Types 1, 4-membered ring)
c     stbn_abc11   stretch-bend parameters for A-B-C atom classes
c     stbn_cba11   stretch-bend parameters for C-B-A atom classes
c     defstbn_abc  default stretch-bend parameters for A-B-C classes
c     defstbn_cba  default stretch-bend parameters for C-B-A classes
c
c
      real*8, allocatable :: stbn_abc(:,:,:)
      real*8, allocatable :: stbn_cba(:,:,:)
      real*8, allocatable :: stbn_abc1(:,:,:)
      real*8, allocatable :: stbn_cba1(:,:,:)
      real*8, allocatable :: stbn_abc2(:,:,:)
      real*8, allocatable :: stbn_cba2(:,:,:)
      real*8, allocatable :: stbn_abc3(:,:,:)
      real*8, allocatable :: stbn_cba3(:,:,:)
      real*8, allocatable :: stbn_abc4(:,:,:)
      real*8, allocatable :: stbn_cba4(:,:,:)
      real*8, allocatable :: stbn_abc5(:,:,:)
      real*8, allocatable :: stbn_cba5(:,:,:)
      real*8, allocatable :: stbn_abc6(:,:,:)
      real*8, allocatable :: stbn_cba6(:,:,:)
      real*8, allocatable :: stbn_abc7(:,:,:)
      real*8, allocatable :: stbn_cba7(:,:,:)
      real*8, allocatable :: stbn_abc8(:,:,:)
      real*8, allocatable :: stbn_cba8(:,:,:)
      real*8, allocatable :: stbn_abc9(:,:,:)
      real*8, allocatable :: stbn_cba9(:,:,:)
      real*8, allocatable :: stbn_abc10(:,:,:)
      real*8, allocatable :: stbn_cba10(:,:,:)
      real*8, allocatable :: stbn_abc11(:,:,:)
      real*8, allocatable :: stbn_cba11(:,:,:)
      real*8 defstbn_abc(0:4,0:4,0:4)
      real*8 defstbn_cba(0:4,0:4,0:4)
c
c
c     t1_1     torsional parameters for 1-fold, MMFF Torsion Type 1
c     t1_2     torsional parameters for 1-fold, MMFF Torsion Type 2
c     t2_1     torsional parameters for 2-fold, MMFF Torsion Type 1
c     t2_2     torsional parameters for 2-fold, MMFF Torsion Type 2
c     t3_1     torsional parameters for 3-fold, MMFF Torsion Type 1
c     t3_2     torsional parameters for 3-fold, MMFF Torsion Type 2
c     kt_1     string of classes for torsions, MMFF Torsion Type 1
c     kt_2     string of classes for torsions, MMFF Torsion Type 2
c
c
      real*8 t1_1(2,0:2000)
      real*8 t2_1(2,0:2000)
      real*8 t3_1(2,0:2000)
      real*8 t1_2(2,0:2000)
      real*8 t2_2(2,0:2000)
      real*8 t3_2(2,0:2000)
      character*16 kt_1(0:2000)
      character*16 kt_2(0:2000)
c
c
c     g        scale factors for calculation of MMFF eps
c     alph     atomic polarizabilities for calculation of MMFF eps
c     nn       effective number of valence electrons for MMFF eps
c     da       donor/acceptor atom classes
c
c
      real*8 g(maxclass)
      real*8 alph(maxclass)
      real*8 nn(maxclass)
      character*1 da(maxclass)
c
c
c     bci      bond charge increments for building atom charges
c     bci_1    bond charge increments for MMFF Bond Type 1   
c     pbci     partial BCI for building missing BCI's
c     fcadj    formal charge adjustment factor
c
c
      real*8 bci(100,100)
      real*8 bci_1(100,100)
      real*8 pbci(maxclass)
      real*8 fcadj(maxclass)
      save
      end
