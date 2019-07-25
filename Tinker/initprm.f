c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine initprm  --  initialize force field parameters  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "initprm" completely initializes a force field by setting all
c     parameters to zero and using defaults for control values
c
c
      subroutine initprm
      use sizes
      use angpot
      use bndpot
      use chgpot
      use fields
      use kanang
      use kangs
      use kantor
      use katoms
      use kbonds
      use kchrge
      use kdipol
      use khbond
      use kiprop
      use kitors
      use kmulti
      use kopbnd
      use kopdst
      use korbs
      use kpitor
      use kpolr
      use kstbnd
      use ksttor
      use ktorsn
      use ktrtor
      use kurybr
      use kvdws
      use kvdwpr
      use math
      use mplpot
      use polpot
      use rxnpot
      use solute
      use urypot
      use torpot
      use units
      use vdwpot
      use ctran
      implicit none
      integer i,j
      character*3 blank3
      character*8 blank8
      character*12 blank12
      character*16 blank16
      character*20 blank20
      character*24 blank24
c
c
c     define blank character strings of various lengths
c
      blank3 = '   '
      blank8 = '        '
      blank12 = '            '
      blank16 = '                '
      blank20 = '                    '
      blank24 = '                        '
c
c     initialize strings of parameter atom types and classes
c
      do i = 1, maxnvp
         kvpr(i) = blank8
      end do
      do i = 1, maxnhb
         khb(i) = blank8
      end do
      do i = 1, maxnb
         kb(i) = blank8
      end do
      do i = 1, maxnb5
         kb5(i) = blank8
      end do
      do i = 1, maxnb4
         kb4(i) = blank8
      end do
      do i = 1, maxnb3
         kb3(i) = blank8
      end do
      do i = 1, maxnel
         kel(i) = blank12
      end do
      do i = 1, maxna
         ka(i) = blank12
      end do
      do i = 1, maxna5
         ka5(i) = blank12
      end do
      do i = 1, maxna4
         ka4(i) = blank12
      end do
      do i = 1, maxna3
         ka3(i) = blank12
      end do
      do i = 1, maxnaf
         kaf(i) = blank12
      end do
      do i = 1, maxnsb
         ksb(i) = blank12
      end do
      do i = 1, maxnu
         ku(i) = blank12
      end do
      do i = 1, maxnopb
         kopb(i) = blank8
      end do
      do i = 1, maxnopd
         kopd(i) = blank16
      end do
      do i = 1, maxndi
         kdi(i) = blank16
      end do
      do i = 1, maxnti
         kti(i) = blank16
      end do
      do i = 1, maxnt
         kt(i) = blank16
      end do
      do i = 1, maxnt5
         kt5(i) = blank16
      end do
      do i = 1, maxnt4
         kt4(i) = blank16
      end do
      do i = 1, maxnpt
         kpt(i) = blank8
      end do
      do i = 1, maxnbt
         kbt(i) = blank16
      end do
      do i = 1, maxnat
         kat(i) = blank16
      end do
      do i = 1, maxntt
         ktt(i) = blank20
      end do
      do i = 1, maxnd
         kd(i) = blank8
      end do
      do i = 1, maxnd5
         kd5(i) = blank8
      end do
      do i = 1, maxnd4
         kd4(i) = blank8
      end do
      do i = 1, maxnd3
         kd3(i) = blank8
      end do
      do i = 1, maxnmp
         kmp(i) = blank12
      end do
      do i = 1, maxnpi
         kpi(i) = blank8
      end do
      do i = 1, maxnpi5
         kpi5(i) = blank8
      end do
      do i = 1, maxnpi4
         kpi4(i) = blank8
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (.not. allocated(atmcls))  allocate (atmcls(maxtyp))
      if (.not. allocated(atmnum))  allocate (atmnum(maxtyp))
      if (.not. allocated(ligand))  allocate (ligand(maxtyp))
      if (.not. allocated(weight))  allocate (weight(maxtyp))
      if (.not. allocated(symbol))  allocate (symbol(maxtyp))
      if (.not. allocated(describe))  allocate (describe(maxtyp))
      if (.not. allocated(rad))  allocate (rad(maxtyp))
      if (.not. allocated(eps))  allocate (eps(maxtyp))
      if (.not. allocated(rad4))  allocate (rad4(maxtyp))
      if (.not. allocated(eps4))  allocate (eps4(maxtyp))
      if (.not. allocated(reduct))  allocate (reduct(maxtyp))
      if (.not. allocated(chg))  allocate (chg(maxtyp))
      if (.not. allocated(polr))  allocate (polr(maxtyp))
      if (.not. allocated(athl))  allocate (athl(maxtyp))
      if (.not. allocated(adird))  allocate (adird(maxtyp))
      if (.not. allocated(apena))  allocate (apena(maxtyp))
      if (.not. allocated(pgrp))  allocate (pgrp(maxval,maxtyp))
c
c     initialize values of some force field parameters
c
      forcefield = blank20
      do i = 1, maxtyp
         symbol(i) = blank3
         atmcls(i) = 0
         atmnum(i) = 0
         weight(i) = 0.0d0
         ligand(i) = 0
         describe(i) = blank24
         rad(i) = 0.0d0
         eps(i) = 0.0d0
         rad4(i) = 0.0d0
         eps4(i) = 0.0d0
         reduct(i) = 0.0d0
         chg(i) = 0.0d0
         polr(i) = 0.0d0
         athl(i) = 0.0d0
         adird(i) = 0.0d0
         apena(i) = 0.0d0
         apre(i) = 0.0d0
         bexp(i) = 0.0d0
         do j = 1, maxval
            pgrp(j,i) = 0
         end do
      end do
      do i = 1, maxclass
         do j = 1, 2
            stbn(j,i) = 0.0d0
         end do
         do j = 1, 3
            anan(j,i) = 0.0d0
         end do
         electron(i) = 0.0d0
         ionize(i) = 0.0d0
         repulse(i) = 0.0d0
      end do
      do i = 1, maxbio
         biotyp(i) = 0
      end do
c
c     set default control parameters for local geometry terms
c
      bndtyp = 'HARMONIC'
      bndunit = 1.0d0
      cbnd = 0.0d0
      qbnd = 0.0d0
      angunit = 1.0d0 / radian**2
      cang = 0.0d0
      qang = 0.0d0
      pang = 0.0d0
      sang = 0.0d0
      stbnunit = 1.0d0 / radian
      ureyunit = 1.0d0
      cury = 0.0d0
      qury = 0.0d0
      aaunit = 1.0d0 / radian**2
      opbtyp = 'W-D-C'
      opbunit = 1.0d0 / radian**2
      copb = 0.0d0
      qopb = 0.0d0
      popb = 0.0d0
      sopb = 0.0d0
      opdunit = 1.0d0
      copd = 0.0d0
      qopd = 0.0d0
      popd = 0.0d0
      sopd = 0.0d0
      idihunit = 1.0d0
      itorunit = 1.0d0
      torsunit = 1.0d0
      ptorunit = 1.0d0
      storunit = 1.0d0
      atorunit = 1.0d0
      ttorunit = 1.0d0
c
c     set default control parameters for van der Waals terms
c
      vdwindex = 'CLASS'
      vdwtyp = 'LENNARD-JONES'
      radrule = 'ARITHMETIC'
      radtyp = 'R-MIN'
      radsiz = 'RADIUS'
      epsrule = 'GEOMETRIC'
      gausstyp = 'NONE'
      ngauss = 0
      abuck = 0.0d0
      bbuck = 0.0d0
      cbuck = 0.0d0
      ghal = 0.12d0
      dhal = 0.07d0
      v2scale = 0.0d0
      v3scale = 0.0d0
      v4scale = 1.0d0
      v5scale = 1.0d0
      use_vcorr = .false.
c   
c     set default control parameters for charge transfer terms
c
       
      ct2scale = 0.0d0
      ct3scale = 0.0d0
      ct4scale = 0.0d0
      ct5scale = 0.0d0
      aprerule = "GEOMETRIC"
      bexprule = "ARITHMETIC"
c
c     set default control parameters for charge-charge terms
c
      electric = coulomb
      dielec = 1.0d0
      ebuffer = 0.0d0
      c2scale = 0.0d0
      c3scale = 0.0d0
      c4scale = 1.0d0
      c5scale = 1.0d0
      neutnbr = .false.
      neutcut = .false.
c
c     set default control parameters for polarizable multipoles
c
      m2scale = 0.0d0
      m3scale = 0.0d0
      m4scale = 1.0d0
      m5scale = 1.0d0
      p2scale = 0.0d0
      p3scale = 0.0d0
      p4scale = 1.0d0
      p5scale = 1.0d0
      p41scale = 0.5d0
c
c     set default control parameters for induced dipoles
c
      poltyp = 'MUTUAL'
      politer = 100
      poleps = 0.000001d0
      udiag = 2.0d0
      d1scale = 0.0d0
      d2scale = 1.0d0
      d3scale = 1.0d0
      d4scale = 1.0d0
      u1scale = 1.0d0
      u2scale = 1.0d0
      u3scale = 1.0d0
      u4scale = 1.0d0
c
c     set default control parameters for implicit solvation
c
      solvtyp = blank8
      borntyp = blank8
c
c     set default control parameters for reaction field
c
      rfsize = 1000000.0d0
      rfbulkd = 80.0d0
      rfterms = 1
c
c     initialize some Merck Molecular force field parameters
c
      call initmmff
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine initmmff  --  initialize some MMFF parameters  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "initmmff" initializes some parameter values for the Merck
c     Molecular force field
c
c
      subroutine initmmff
      use ktorsn
      use merck
      implicit none
      integer i,j,k
      character*16 blank16
c
c
c     define blank character strings of various lengths
c
      blank16 = '                '
c
c     perform dynamic allocation of some global arrays
c
      if (.not. allocated(mmff_ka))
     &   allocate (mmff_ka(0:100,100,0:100))
      if (.not. allocated(mmff_ka1))
     &   allocate (mmff_ka1(0:100,100,0:100))
      if (.not. allocated(mmff_ka2))
     &   allocate (mmff_ka2(0:100,100,0:100))
      if (.not. allocated(mmff_ka3))
     &   allocate (mmff_ka3(0:100,100,0:100))
      if (.not. allocated(mmff_ka4))
     &   allocate (mmff_ka4(0:100,100,0:100))
      if (.not. allocated(mmff_ka5))
     &   allocate (mmff_ka5(0:100,100,0:100))
      if (.not. allocated(mmff_ka6))
     &   allocate (mmff_ka6(0:100,100,0:100))
      if (.not. allocated(mmff_ka7))
     &   allocate (mmff_ka7(0:100,100,0:100))
      if (.not. allocated(mmff_ka8))
     &   allocate (mmff_ka8(0:100,100,0:100))
      if (.not. allocated(mmff_ang0))
     &   allocate (mmff_ang0(0:100,100,0:100))
      if (.not. allocated(mmff_ang1))
     &   allocate (mmff_ang1(0:100,100,0:100))
      if (.not. allocated(mmff_ang2))
     &   allocate (mmff_ang2(0:100,100,0:100))
      if (.not. allocated(mmff_ang3))
     &   allocate (mmff_ang3(0:100,100,0:100))
      if (.not. allocated(mmff_ang4))
     &   allocate (mmff_ang4(0:100,100,0:100))
      if (.not. allocated(mmff_ang5))
     &   allocate (mmff_ang5(0:100,100,0:100))
      if (.not. allocated(mmff_ang6))
     &   allocate (mmff_ang6(0:100,100,0:100))
      if (.not. allocated(mmff_ang7))
     &   allocate (mmff_ang7(0:100,100,0:100))
      if (.not. allocated(mmff_ang8))
     &   allocate (mmff_ang8(0:100,100,0:100))
      if (.not. allocated(stbn_abc))
     &   allocate (stbn_abc(100,100,100))
      if (.not. allocated(stbn_cba))
     &   allocate (stbn_cba(100,100,100))
      if (.not. allocated(stbn_abc1))
     &   allocate (stbn_abc1(100,100,100))
      if (.not. allocated(stbn_cba1))
     &   allocate (stbn_cba1(100,100,100))
      if (.not. allocated(stbn_abc2))
     &   allocate (stbn_abc2(100,100,100))
      if (.not. allocated(stbn_cba2))
     &   allocate (stbn_cba2(100,100,100))
      if (.not. allocated(stbn_abc3))
     &   allocate (stbn_abc3(100,100,100))
      if (.not. allocated(stbn_cba3))
     &   allocate (stbn_cba3(100,100,100))
      if (.not. allocated(stbn_abc4))
     &   allocate (stbn_abc4(100,100,100))
      if (.not. allocated(stbn_cba4))
     &   allocate (stbn_cba4(100,100,100))
      if (.not. allocated(stbn_abc5))
     &   allocate (stbn_abc5(100,100,100))
      if (.not. allocated(stbn_cba5))
     &   allocate (stbn_cba5(100,100,100))
      if (.not. allocated(stbn_abc6))
     &   allocate (stbn_abc6(100,100,100))
      if (.not. allocated(stbn_cba6))
     &   allocate (stbn_cba6(100,100,100))
      if (.not. allocated(stbn_abc7))
     &   allocate (stbn_abc7(100,100,100))
      if (.not. allocated(stbn_cba7))
     &   allocate (stbn_cba7(100,100,100))
      if (.not. allocated(stbn_abc8))
     &   allocate (stbn_abc8(100,100,100))
      if (.not. allocated(stbn_cba8))
     &   allocate (stbn_cba8(100,100,100))
      if (.not. allocated(stbn_abc9))
     &   allocate (stbn_abc9(100,100,100))
      if (.not. allocated(stbn_cba9))
     &   allocate (stbn_cba9(100,100,100))
      if (.not. allocated(stbn_abc10))
     &   allocate (stbn_abc10(100,100,100))
      if (.not. allocated(stbn_cba10))
     &   allocate (stbn_cba10(100,100,100))
      if (.not. allocated(stbn_abc11))
     &   allocate (stbn_abc11(100,100,100))
      if (.not. allocated(stbn_cba11))
     &   allocate (stbn_cba11(100,100,100))
c
c     initialize values for MMFF atom class equivalencies
c
      do i = 1, 5
         do j = 1, 500
            eqclass(j,i) = 1000
         end do
      end do
c
c     initialize values for MMFF aromatic ring parameters
c
      do i = 1, 6
         do j = 1, maxtyp
            mmffarom(j,i) = 0
            mmffaromc(j,i) = 0
            mmffaroma(j,i) = 0
         end do
      end do
c
c     initialize values for MMFF bond stretching parameters
c
      do i = 1, 100
         do j = 1, 100
            mmff_kb(j,i) = 1000.0d0
            mmff_kb1(j,i) = 1000.0d0
            mmff_b0(j,i) = 1000.0d0
            mmff_b1(j,i) = 1000.0d0
         end do
      end do
c
c     initialize values for MMFF angle bending parameters
c
      do i = 0, 100
         do j = 1, 100
            do k = 0, 100
               mmff_ka(k,j,i) = 1000.0d0
               mmff_ka1(k,j,i) = 1000.0d0
               mmff_ka2(k,j,i) = 1000.0d0
               mmff_ka3(k,j,i) = 1000.0d0
               mmff_ka4(k,j,i) = 1000.0d0
               mmff_ka5(k,j,i) = 1000.0d0
               mmff_ka6(k,j,i) = 1000.0d0
               mmff_ka7(k,j,i) = 1000.0d0
               mmff_ka8(k,j,i) = 1000.0d0
               mmff_ang0(k,j,i) = 1000.0d0
               mmff_ang1(k,j,i) = 1000.0d0
               mmff_ang2(k,j,i) = 1000.0d0
               mmff_ang3(k,j,i) = 1000.0d0
               mmff_ang4(k,j,i) = 1000.0d0
               mmff_ang5(k,j,i) = 1000.0d0
               mmff_ang6(k,j,i) = 1000.0d0
               mmff_ang7(k,j,i) = 1000.0d0
               mmff_ang8(k,j,i) = 1000.0d0
            end do
         end do
      end do
c
c     initialize values for MMFF stretch-bend parameters
c
      do i = 1, 100
         do j = 1, 100
            do k = 1, 100
               stbn_abc(k,j,i) = 1000.0d0
               stbn_cba(k,j,i) = 1000.0d0
               stbn_abc1(k,j,i) = 1000.0d0
               stbn_cba1(k,j,i) = 1000.0d0
               stbn_abc2(k,j,i) = 1000.0d0
               stbn_cba2(k,j,i) = 1000.0d0
               stbn_abc3(k,j,i) = 1000.0d0
               stbn_cba3(k,j,i) = 1000.0d0
               stbn_abc4(k,j,i) = 1000.0d0
               stbn_cba4(k,j,i) = 1000.0d0
               stbn_abc5(k,j,i) = 1000.0d0
               stbn_cba5(k,j,i) = 1000.0d0
               stbn_abc6(k,j,i) = 1000.0d0
               stbn_cba6(k,j,i) = 1000.0d0
               stbn_abc7(k,j,i) = 1000.0d0
               stbn_cba7(k,j,i) = 1000.0d0
               stbn_abc8(k,j,i) = 1000.0d0
               stbn_cba8(k,j,i) = 1000.0d0
               stbn_abc9(k,j,i) = 1000.0d0
               stbn_cba9(k,j,i) = 1000.0d0
               stbn_abc10(k,j,i) = 1000.0d0
               stbn_cba10(k,j,i) = 1000.0d0
               stbn_abc11(k,j,i) = 1000.0d0
               stbn_cba11(k,j,i) = 1000.0d0
            end do
         end do
      end do
c
c     initialize values for MMFF torsional parameters
c
      do i = 1, maxnt
         kt(i) = blank16
         kt_1(i) = blank16
         kt_2(i) = blank16
         t1(1,i) = 1000.0d0
         t1(2,i) = 1000.0d0
         t2(1,i) = 1000.0d0
         t2(2,i) = 1000.0d0
         t3(1,i) = 1000.0d0
         t3(2,i) = 1000.0d0
         t1_1(1,i) = 1000.0d0
         t1_1(2,i) = 1000.0d0
         t2_1(1,i) = 1000.0d0
         t2_1(2,i) = 1000.0d0
         t3_1(1,i) = 1000.0d0
         t3_1(2,i) = 1000.0d0
         t1_2(1,i) = 1000.0d0
         t1_2(2,i) = 1000.0d0
         t2_2(1,i) = 1000.0d0
         t2_2(2,i) = 1000.0d0
         t3_2(1,i) = 1000.0d0
         t3_2(2,i) = 1000.0d0
      end do
c
c     initialize values for MMFF bond charge increment parameters
c
      do i = 1, 100
         do j = 1, 100
            bci(j,i) = 1000.0d0
            bci_1(j,i) = 1000.0d0
         end do
      end do
      return
      end
