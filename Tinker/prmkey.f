c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1996  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine prmkey  --  interpret force field keywords  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "prmkey" parses a text string to extract keywords related to
c     force field potential energy functional forms and constants
c
c
      subroutine prmkey (text)
      use angpot
      use bndpot
      use chgpot
      use fields
      use mplpot
      use polpot
      use potent
      use rxnpot
      use torpot
      use urypot
      use vdwpot
!CW
      use ctran
      implicit none
      integer next
      character*4 value
      character*20 keyword
      character*240 text
      character*240 record
      character*240 string
c
c
c     parse the line to extract any possible keyword
c
      record = text
      next = 1
      call upcase (record)
      call gettext (record,keyword,next)
      string = record(next:240)
c
c     select the individual force field potential terms
c
      if (keyword(1:9) .eq. 'BONDTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_bond = .true.
         if (value .eq. 'NONE')  use_bond = .false.
      else if (keyword(1:10) .eq. 'ANGLETERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_angle = .true.
         if (value .eq. 'NONE')  use_angle = .false.
      else if (keyword(1:11) .eq. 'STRBNDTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_strbnd = .true.
         if (value .eq. 'NONE')  use_strbnd = .false.
      else if (keyword(1:9) .eq. 'UREYTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_urey = .true.
         if (value .eq. 'NONE')  use_urey = .false.
      else if (keyword(1:11) .eq. 'ANGANGTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_angang = .true.
         if (value .eq. 'NONE')  use_angang = .false.
      else if (keyword(1:11) .eq. 'OPBENDTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_opbend = .true.
         if (value .eq. 'NONE')  use_opbend = .false.
      else if (keyword(1:11) .eq. 'OPDISTTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_opdist = .true.
         if (value .eq. 'NONE')  use_opdist = .false.
      else if (keyword(1:11) .eq. 'IMPROPTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_improp = .true.
         if (value .eq. 'NONE')  use_improp = .false.
      else if (keyword(1:12) .eq. 'IMPTORSTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_imptor = .true.
         if (value .eq. 'NONE')  use_imptor = .false.
      else if (keyword(1:12) .eq. 'TORSIONTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_tors = .true.
         if (value .eq. 'NONE')  use_tors = .false.
      else if (keyword(1:11) .eq. 'PITORSTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_pitors = .true.
         if (value .eq. 'NONE')  use_pitors = .false.
      else if (keyword(1:11) .eq. 'STRTORTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_strtor = .true.
         if (value .eq. 'NONE')  use_strtor = .false.
      else if (keyword(1:11) .eq. 'ANGTORTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_angtor = .true.
         if (value .eq. 'NONE')  use_angtor = .false.
      else if (keyword(1:11) .eq. 'TORTORTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_tortor = .true.
         if (value .eq. 'NONE')  use_tortor = .false.
      else if (keyword(1:8) .eq. 'VDWTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_vdw = .true.
         if (value .eq. 'NONE')  use_vdw = .false.
      else if (keyword(1:7) .eq. 'CTTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_ct = .true.
         if (value .eq. 'NONE')  use_ct = .false.
      else if (keyword(1:11) .eq. 'CHARGETERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_charge = .true.
         if (value .eq. 'NONE')  use_charge = .false.
      else if (keyword(1:11) .eq. 'CHGDPLTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_chgdpl = .true.
         if (value .eq. 'NONE')  use_chgdpl = .false.
      else if (keyword(1:11) .eq. 'DIPOLETERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_dipole = .true.
         if (value .eq. 'NONE')  use_dipole = .false.
      else if (keyword(1:10) .eq. 'MPOLETERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_mpole = .true.
         if (value .eq. 'NONE')  use_mpole = .false.
      else if (keyword(1:13) .eq. 'POLARIZETERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_polar = .true.
         if (value .eq. 'NONE')  use_polar = .false.
      else if (keyword(1:13) .eq. 'RXNFIELDTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_rxnfld = .true.
         if (value .eq. 'NONE')  use_rxnfld = .false.
      else if (keyword(1:12) .eq. 'SOLVATETERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_solv = .true.
         if (value .eq. 'NONE')  use_solv = .false.
      else if (keyword(1:12) .eq. 'METALTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_metal = .true.
         if (value .eq. 'NONE')  use_metal = .false.
      else if (keyword(1:13) .eq. 'RESTRAINTERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_geom = .true.
         if (value .eq. 'NONE')  use_geom = .false.
      else if (keyword(1:10) .eq. 'EXTRATERM ') then
         call getword (record,value,next)
         if (value .eq. 'ONLY')  call potoff
         use_extra = .true.
         if (value .eq. 'NONE')  use_extra = .false.
      end if
c
c     select the name of the force field parameter set
c
      if (keyword(1:11) .eq. 'FORCEFIELD ') then
         call getword (record,forcefield,next)
c
c     set control parameters for bond stretching potentials
c
      else if (keyword(1:9) .eq. 'BONDTYPE ') then
         call getword (record,bndtyp,next)
      else if (keyword(1:9) .eq. 'BONDUNIT ') then
         read (string,*,err=10,end=10)  bndunit
      else if (keyword(1:11) .eq. 'BOND-CUBIC ') then
         read (string,*,err=10,end=10)  cbnd
      else if (keyword(1:13) .eq. 'BOND-QUARTIC ') then
         read (string,*,err=10,end=10)  qbnd
c
c     set control parameters for bond angle bending potentials
c
      else if (keyword(1:10) .eq. 'ANGLEUNIT ') then
         read (string,*,err=10,end=10)  angunit
      else if (keyword(1:12) .eq. 'ANGLE-CUBIC ') then
         read (string,*,err=10,end=10)  cang
      else if (keyword(1:14) .eq. 'ANGLE-QUARTIC ') then
         read (string,*,err=10,end=10)  qang
      else if (keyword(1:13) .eq. 'ANGLE-PENTIC ') then
         read (string,*,err=10,end=10)  pang
      else if (keyword(1:13) .eq. 'ANGLE-SEXTIC ') then
         read (string,*,err=10,end=10)  sang
c
c     set control parameters for stretch-bend potentials
c
      else if (keyword(1:11) .eq. 'STRBNDUNIT ') then
         read (string,*,err=10,end=10)  stbnunit
c
c     set control parameters for Urey-Bradley potentials
c
      else if (keyword(1:9) .eq. 'UREYUNIT ') then
         read (string,*,err=10,end=10)  ureyunit
      else if (keyword(1:11) .eq. 'UREY-CUBIC ') then
         read (string,*,err=10,end=10)  cury
      else if (keyword(1:13) .eq. 'UREY-QUARTIC ') then
         read (string,*,err=10,end=10)  qury
c
c     set control parameters for out-of-plane bend potentials
c
      else if (keyword(1:11) .eq. 'OPBENDTYPE ') then
         call getword (record,opbtyp,next)
      else if (keyword(1:11) .eq. 'OPBENDUNIT ') then
         read (string,*,err=10,end=10)  opbunit
      else if (keyword(1:13) .eq. 'OPBEND-CUBIC ') then
         read (string,*,err=10,end=10)  copb
      else if (keyword(1:15) .eq. 'OPBEND-QUARTIC ') then
         read (string,*,err=10,end=10)  qopb
      else if (keyword(1:14) .eq. 'OPBEND-PENTIC ') then
         read (string,*,err=10,end=10)  popb
      else if (keyword(1:14) .eq. 'OPBEND-SEXTIC ') then
         read (string,*,err=10,end=10)  sopb
c
c     set control parameters for out-of-plane distance potentials
c
      else if (keyword(1:11) .eq. 'OPDISTUNIT ') then
         read (string,*,err=10,end=10)  opdunit
      else if (keyword(1:13) .eq. 'OPDIST-CUBIC ') then
         read (string,*,err=10,end=10)  copd
      else if (keyword(1:15) .eq. 'OPDIST-QUARTIC ') then
         read (string,*,err=10,end=10)  qopd
      else if (keyword(1:14) .eq. 'OPDIST-PENTIC ') then
         read (string,*,err=10,end=10)  popd
      else if (keyword(1:14) .eq. 'OPDIST-SEXTIC ') then
         read (string,*,err=10,end=10)  sopd
c
c     set control parameters for other local geometry potentials
c
      else if (keyword(1:11) .eq. 'ANGANGUNIT ') then
         read (string,*,err=10,end=10)  aaunit
      else if (keyword(1:11) .eq. 'IMPROPUNIT ') then
         read (string,*,err=10,end=10)  idihunit
      else if (keyword(1:11) .eq. 'IMPTORUNIT ') then
         read (string,*,err=10,end=10)  itorunit
      else if (keyword(1:12) .eq. 'TORSIONUNIT ') then
         read (string,*,err=10,end=10)  torsunit
      else if (keyword(1:11) .eq. 'PITORSUNIT ') then
         read (string,*,err=10,end=10)  ptorunit
      else if (keyword(1:11) .eq. 'STRTORUNIT ') then
         read (string,*,err=10,end=10)  storunit
      else if (keyword(1:11) .eq. 'ANGTORUNIT ') then
         read (string,*,err=10,end=10)  atorunit
      else if (keyword(1:11) .eq. 'TORTORUNIT ') then
         read (string,*,err=10,end=10)  ttorunit
c
c     set control parameters for van der Waals potentials
c
      else if (keyword(1:9) .eq. 'VDWINDEX ') then
         call getword (record,vdwindex,next)
      else if (keyword(1:8) .eq. 'VDWTYPE ') then
         call getword (record,vdwtyp,next)
      else if (keyword(1:11) .eq. 'RADIUSTYPE ') then
         call getword (record,radtyp,next)
      else if (keyword(1:11) .eq. 'RADIUSSIZE ') then
         call getword (record,radsiz,next)
      else if (keyword(1:11) .eq. 'RADIUSRULE ') then
         call getword (record,radrule,next)
      else if (keyword(1:12) .eq. 'EPSILONRULE ') then
         call getword (record,epsrule,next)
      else if (keyword(1:14) .eq. 'GAUSSTYPE ') then
         call getword (record,gausstyp,next)
      else if (keyword(1:10) .eq. 'A-EXPTERM ') then
         read (string,*,err=10,end=10)  abuck
      else if (keyword(1:10) .eq. 'B-EXPTERM ') then
         read (string,*,err=10,end=10)  bbuck
      else if (keyword(1:10) .eq. 'C-EXPTERM ') then
         read (string,*,err=10,end=10)  cbuck
      else if (keyword(1:14) .eq. 'GAMMA-HALGREN ') then
         read (string,*,err=10,end=10)  ghal
      else if (keyword(1:14) .eq. 'DELTA-HALGREN ') then
         read (string,*,err=10,end=10)  dhal
      else if (keyword(1:13) .eq. 'VDW-12-SCALE ') then
         read (string,*,err=10,end=10)  v2scale
         if (v2scale .gt. 1.0d0)  v2scale = 1.0d0 / v2scale
      else if (keyword(1:13) .eq. 'VDW-13-SCALE ') then
         read (string,*,err=10,end=10)  v3scale
         if (v3scale .gt. 1.0d0)  v3scale = 1.0d0 / v3scale
      else if (keyword(1:13) .eq. 'VDW-14-SCALE ') then
         read (string,*,err=10,end=10)  v4scale
         if (v4scale .gt. 1.0d0)  v4scale = 1.0d0 / v4scale
      else if (keyword(1:13) .eq. 'VDW-15-SCALE ') then
         read (string,*,err=10,end=10)  v5scale
         if (v5scale .gt. 1.0d0)  v5scale = 1.0d0 / v5scale
      else if (keyword(1:15) .eq. 'VDW-CORRECTION ') then
         use_vcorr = .true.

c CW
c     set charge transfer  
      else if (keyword(1:9) .eq. 'APRERULE ') then
         call getword (record,aprerule,next)
      else if (keyword(1:9) .eq. 'BEXPRULE ') then
         call getword (record,bexprule,next)
      else if (keyword(1:13) .eq. 'CT-12-SCALE ') then
         read (string,*,err=10,end=10)  ct2scale
         if (ct2scale .gt. 1.0d0)  ct2scale = 1.0d0 / ct2scale
      else if (keyword(1:13) .eq. 'CT-13-SCALE ') then
         read (string,*,err=10,end=10)  ct3scale
         if (ct2scale .gt. 1.0d0)  ct3scale = 1.0d0 / ct3scale
      else if (keyword(1:13) .eq. 'CT-14-SCALE ') then
         read (string,*,err=10,end=10)  ct4scale
         if (ct4scale .gt. 1.0d0)  ct4scale = 1.0d0 / ct4scale
      else if (keyword(1:13) .eq. 'CT-15-SCALE ') then
         read (string,*,err=10,end=10)  ct5scale
         if (ct5scale .gt. 1.0d0)  ct5scale = 1.0d0 / ct5scale
c
c     set control parameters for charge-charge potentials
c
      else if (keyword(1:9) .eq. 'ELECTRIC ') then
         read (string,*,err=10,end=10)  electric
      else if (keyword(1:11) .eq. 'DIELECTRIC ') then
         read (string,*,err=10,end=10)  dielec
      else if (keyword(1:11) .eq. 'CHG-BUFFER ') then
         read (string,*,err=10,end=10)  ebuffer
      else if (keyword(1:13) .eq. 'CHG-12-SCALE ') then
         read (string,*,err=10,end=10)  c2scale
         if (c2scale .gt. 1.0d0)  c2scale = 1.0d0 / c2scale
      else if (keyword(1:13) .eq. 'CHG-13-SCALE ') then
         read (string,*,err=10,end=10)  c3scale
         if (c3scale .gt. 1.0d0)  c3scale = 1.0d0 / c3scale
      else if (keyword(1:13) .eq. 'CHG-14-SCALE ') then
         read (string,*,err=10,end=10)  c4scale
         if (c4scale .gt. 1.0d0)  c4scale = 1.0d0 / c4scale
      else if (keyword(1:13) .eq. 'CHG-15-SCALE ') then
         read (string,*,err=10,end=10)  c5scale
         if (c5scale .gt. 1.0d0)  c5scale = 1.0d0 / c5scale
      else if (keyword(1:16) .eq. 'NEIGHBOR-GROUPS ') then
         neutnbr = .true.
      else if (keyword(1:15) .eq. 'NEUTRAL-GROUPS ') then
         neutcut = .true.
c
c     set control parameters for atomic multipole potentials
c
      else if (keyword(1:15) .eq. 'MPOLE-12-SCALE ') then
         read (string,*,err=10,end=10)  m2scale
         if (m2scale .gt. 1.0d0)  m2scale = 1.0d0 / m2scale
      else if (keyword(1:15) .eq. 'MPOLE-13-SCALE ') then
         read (string,*,err=10,end=10)  m3scale
         if (m3scale .gt. 1.0d0)  m3scale = 1.0d0 / m3scale
      else if (keyword(1:15) .eq. 'MPOLE-14-SCALE ') then
         read (string,*,err=10,end=10)  m4scale
         if (m4scale .gt. 1.0d0)  m4scale = 1.0d0 / m4scale
      else if (keyword(1:15) .eq. 'MPOLE-15-SCALE ') then
         read (string,*,err=10,end=10)  m5scale
         if (m5scale .gt. 1.0d0)  m5scale = 1.0d0 / m5scale
c
c     set control parameters for polarization potentials
c
      else if (keyword(1:13) .eq. 'POLARIZATION ') then
         call getword (record,poltyp,next)
      else if (keyword(1:11) .eq. 'POLAR-ITER ') then
         read (string,*,err=10,end=10)  politer
      else if (keyword(1:10) .eq. 'POLAR-EPS ') then
         read (string,*,err=10,end=10)  poleps
      else if (keyword(1:12) .eq. 'USOLVE-DIAG ') then
         read (string,*,err=10,end=10)  udiag
      else if (keyword(1:15) .eq. 'POLAR-12-SCALE ') then
         read (string,*,err=10,end=10)  p2scale
         if (p2scale .gt. 1.0d0)  p2scale = 1.0d0 / p2scale
      else if (keyword(1:15) .eq. 'POLAR-13-SCALE ') then
         read (string,*,err=10,end=10)  p3scale
         if (p3scale .gt. 1.0d0)  p3scale = 1.0d0 / p3scale
      else if (keyword(1:15) .eq. 'POLAR-14-SCALE ') then
         read (string,*,err=10,end=10)  p4scale
         if (p4scale .gt. 1.0d0)  p4scale = 1.0d0 / p4scale
      else if (keyword(1:15) .eq. 'POLAR-15-SCALE ') then
         read (string,*,err=10,end=10)  p5scale
         if (p5scale .gt. 1.0d0)  p5scale = 1.0d0 / p5scale
      else if (keyword(1:15) .eq. 'POLAR-14-INTRA ') then
         read (string,*,err=10,end=10)  p41scale
         if (p41scale .gt. 1.0d0)  p41scale = 1.0d0 / p41scale
      else if (keyword(1:16) .eq. 'DIRECT-11-SCALE ') then
         read (string,*,err=10,end=10)  d1scale
         if (d1scale .gt. 1.0d0)  d1scale = 1.0d0 / d1scale
      else if (keyword(1:16) .eq. 'DIRECT-12-SCALE ') then
         read (string,*,err=10,end=10)  d2scale
         if (d2scale .gt. 1.0d0)  d2scale = 1.0d0 / d2scale
      else if (keyword(1:16) .eq. 'DIRECT-13-SCALE ') then
         read (string,*,err=10,end=10)  d3scale
         if (d3scale .gt. 1.0d0)  d3scale = 1.0d0 / d3scale
      else if (keyword(1:16) .eq. 'DIRECT-14-SCALE ') then
         read (string,*,err=10,end=10)  d4scale
         if (d4scale .gt. 1.0d0)  d4scale = 1.0d0 / d4scale
      else if (keyword(1:16) .eq. 'MUTUAL-11-SCALE ') then
         read (string,*,err=10,end=10)  u1scale
         if (u1scale .gt. 1.0d0)  u1scale = 1.0d0 / u1scale
      else if (keyword(1:16) .eq. 'MUTUAL-12-SCALE ') then
         read (string,*,err=10,end=10)  u2scale
         if (u2scale .gt. 1.0d0)  u2scale = 1.0d0 / u2scale
      else if (keyword(1:16) .eq. 'MUTUAL-13-SCALE ') then
         read (string,*,err=10,end=10)  u3scale
         if (u3scale .gt. 1.0d0)  u3scale = 1.0d0 / u3scale
      else if (keyword(1:16) .eq. 'MUTUAL-14-SCALE ') then
         read (string,*,err=10,end=10)  u4scale
         if (u4scale .gt. 1.0d0)  u4scale = 1.0d0 / u4scale
c
c     set control parameters for reaction field potentials
c
      else if (keyword(1:14) .eq. 'REACTIONFIELD ') then
         read (string,*,err=10,end=10)  rfsize,rfbulkd,rfterms
      end if
c
c     jump directly to the end if any error was detected
c
   10 continue
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine potoff  --  turn off all potential functions  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "potoff" clears the forcefield definition by turning off
c     the use of each of the potential energy functions
c
c
      subroutine potoff
      use potent
      implicit none
c
c
c     turn off the use of each of the potential energy functions
c
      use_bond = .false.
      use_angle = .false.
      use_strbnd = .false.
      use_urey = .false.
      use_angang = .false.
      use_opbend = .false.
      use_opdist = .false.
      use_improp = .false.
      use_imptor = .false.
      use_tors = .false.
      use_pitors = .false.
      use_strtor = .false.
      use_angtor = .false.
      use_tortor = .false.
      use_vdw = .false.
      use_ct = .false.
      use_charge = .false.
      use_chgdpl = .false.
      use_dipole = .false.
      use_mpole = .false.
      use_polar = .false.
      use_rxnfld = .false.
      use_solv = .false.
      use_metal = .false.
      use_geom = .false.
      use_extra = .false.
      return
      end
