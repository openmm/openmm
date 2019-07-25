c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  function energy  --  evaluates energy terms and total  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "energy" calls the subroutines to calculate the potential
c     energy terms and sums up to form the total energy
c
c
      function energy ()
      use sizes
      use bound
      use energi
      use iounit
      use limits
      use potent
      use rigid
      use vdwpot
      implicit none
      real*8 energy
      real*8 cutoff
c
c
c     zero out each of the potential energy components
c
      eb = 0.0d0
      ea = 0.0d0
      eba = 0.0d0
      eub = 0.0d0
      eaa = 0.0d0
      eopb = 0.0d0
      eopd = 0.0d0
      eid = 0.0d0
      eit = 0.0d0
      et = 0.0d0
      ept = 0.0d0
      ebt = 0.0d0
      eat = 0.0d0
      ett = 0.0d0
      ev = 0.0d0
CW!!!
      ect = 0.0d0

      ec = 0.0d0
      ecd = 0.0d0
      ed = 0.0d0
      em = 0.0d0
      ep = 0.0d0
      er = 0.0d0
      es = 0.0d0
      elf = 0.0d0
      eg = 0.0d0
      ex = 0.0d0
c
c     maintain any periodic boundary conditions
c
      if (use_bounds .and. .not.use_rigid)  call bounds
c
c     update the pairwise interaction neighbor lists
c
      if (use_list)  call nblist
c
c     remove any previous use of the replicates method
c
      cutoff = 0.0d0
      call replica (cutoff)
c
c     many implicit solvation models require Born radii
c
      if (use_born)  call born
c
c     alter bond and torsion constants for pisystem
c
      if (use_orbit)  call picalc
c
c     call the local geometry energy component routines
c
      if (use_bond)  call ebond
      if (use_angle)  call eangle
      if (use_strbnd)  call estrbnd
      if (use_urey)  call eurey
      if (use_angang)  call eangang
      if (use_opbend)  call eopbend
      if (use_opdist)  call eopdist
      if (use_improp)  call eimprop
      if (use_imptor)  call eimptor
      if (use_tors)  call etors
      if (use_pitors)  call epitors
      if (use_strtor)  call estrtor
      if (use_angtor)  call eangtor
      if (use_tortor)  call etortor
c
c     call the van der Waals energy component routines
c
      if (use_vdw) then
         if (vdwtyp .eq. 'LENNARD-JONES')  call elj
         if (vdwtyp .eq. 'BUCKINGHAM')  call ebuck
         if (vdwtyp .eq. 'MM3-HBOND')  call emm3hb
         if (vdwtyp .eq. 'BUFFERED-14-7')  call ehal
         if (vdwtyp .eq. 'GAUSSIAN')  call egauss
      end if
c
CW 
c      call the charge transfer energy component routines
      !if (use_vdw) call ect0
      if (use_ct) call ect0

c     call the electrostatic energy component routines
c
      if (use_charge)  call echarge
      if (use_chgdpl)  call echgdpl
      if (use_dipole)  call edipole
      if (use_mpole)  call empole
      if (use_polar)  call epolar
      if (use_rxnfld)  call erxnfld
c
c     call any miscellaneous energy component routines
c
      if (use_solv)  call esolv
      if (use_geom)  call egeom
      if (use_metal)  call emetal
      if (use_extra)  call extra
c
c     sum up to give the total potential energy
c
      esum = eb + ea + eba + eub + eaa + eopb + eopd + eid + eit
     &          + et + ept + ebt + eat + ett + ev + ec + ecd + ed
     &          + em + ep + er + es + elf + eg + ex + ect !CW
      energy = esum
c
c     check for an illegal value for the total energy
c
c     if (isnan(esum)) then
      if (esum .ne. esum) then
         write (iout,10)
   10    format (/,' ENERGY  --  Illegal Value for the Total',
     &              ' Potential Energy')
         call fatal
      end if
      return
      end
