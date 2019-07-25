c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2000  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine emetal3  --  ligand field energy and analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "emetal3" calculates the transition metal ligand field energy
c     and also partitions the energy among the atoms
c
c
      subroutine emetal3
      use sizes
      use action
      use analyz
      use atomid
      use atoms
      use energi
      use kchrge
      implicit none
      integer i
c
c
c     zero out the ligand field energy and partitioning terms
c
      nelf = 0
      elf = 0.0d0
      do i = 1, n
         aelf(i) = 0.0d0
      end do
c
c     for now, just count the sites and call the energy code
c
      do i = 1, n
         if (atomic(i).eq.29 .and. chg(type(i)).eq.2.0d0) then
            nelf = nelf + 1
         end if
      end do
      call emetal
      return
      end
