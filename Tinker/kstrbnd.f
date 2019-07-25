c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1991  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine kstrbnd  --  assign stretch-bend parameters  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "kstrbnd" assigns parameters for stretch-bend interactions
c     and processes new or changed parameter values
c
c
      subroutine kstrbnd
      use sizes
      use angbnd
      use angpot
      use atmlst
      use atomid
      use atoms
      use couple
      use fields
      use inform
      use iounit
      use keys
      use kstbnd
      use potent
      use strbnd
      implicit none
      integer i,j,k,nsb
      integer ia,ib,ic
      integer ita,itb,itc
      integer nba,nbc
      integer size,next
      real*8 sb1,sb2,temp
      logical header
      character*4 pa,pb,pc
      character*12 blank,pt
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     process keywords containing stretch-bend parameters
c
      blank = '            '
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'STRBND ') then
            ia = 0
            ib = 0
            ic = 0
            sb1 = 0.0d0
            sb2 = 0.0d0
            string = record(next:240)
            read (string,*,err=10,end=10)  ia,ib,ic,sb1,sb2
   10       continue
            if (.not. silent) then
               if (header) then
                  header = .false.
                  write (iout,20)
   20             format (/,' Additional Stretch-Bend Parameters :',
     &                    //,5x,'Atom Classes',6x,'K(SB)-1',5x,
     &                       'K(SB)-2',/)
               end if
               write (iout,30)  ia,ib,ic,sb1,sb2
   30          format (4x,3i4,2x,2f12.3)
            end if
            size = 4
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            if (ia .le. ic) then
               pt = pa//pb//pc
            else
               pt = pc//pb//pa
               temp = sb1
               sb1 = sb2
               sb2 = temp
            end if
            do j = 1, maxnsb
               if (ksb(j).eq.blank .or. ksb(j).eq.pt) then
                  ksb(j) = pt
                  stbn(1,j) = sb1
                  stbn(2,j) = sb2
                  goto 50
               end if
            end do
            write (iout,40)
   40       format (/,' KSTRBND  --  Too many Stretch-Bend',
     &                 ' Interaction Parameters')
            abort = .true.
   50       continue
         end if
      end do
c
c     determine the total number of forcefield parameters
c
      nsb = maxnsb
      do i = maxnsb, 1, -1
         if (ksb(i) .eq. blank)  nsb = i - 1
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(isb))  deallocate (isb)
      if (allocated(sbk))  deallocate (sbk)
      allocate (isb(3,nangle))
      allocate (sbk(2,nangle))
c
c     use special stretch-bend parameter assignment method for MMFF
c
      if (forcefield .eq. 'MMFF94') then
         call kstrbndm
         return
      end if
c
c     assign the stretch-bend parameters for each angle
c
      nstrbnd = 0
      if (nsb .ne. 0) then
         do i = 1, nangle
            ia = iang(1,i)
            ib = iang(2,i)
            ic = iang(3,i)
            ita = class(ia)
            itb = class(ib)
            itc = class(ic)
            size = 4
            call numeral (ita,pa,size)
            call numeral (itb,pb,size)
            call numeral (itc,pc,size)
            if (ita .le. itc) then
               pt = pa//pb//pc
            else
               pt = pc//pb//pa
            end if
            do j = 1, nsb
               if (ksb(j) .eq. pt) then
                  nstrbnd = nstrbnd + 1
                  do k = 1, n12(ib)
                     if (i12(k,ib) .eq. ia)  nba = bndlist(k,ib)
                     if (i12(k,ib) .eq. ic)  nbc = bndlist(k,ib)
                  end do
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nba
                  isb(3,nstrbnd) = nbc
                  if (ita .le. itc) then
                     sbk(1,nstrbnd) = stbn(1,j)
                     sbk(2,nstrbnd) = stbn(2,j)
                  else
                     sbk(1,nstrbnd) = stbn(2,j)
                     sbk(2,nstrbnd) = stbn(1,j)
                  end if
                  goto 60
               end if
            end do
   60       continue
         end do
      end if
c
c     turn off the stretch-bend potential if it is not used
c
      if (nstrbnd .eq. 0)  use_strbnd = .false.
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kstrbndm  --  assign MMFF str-bnd parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kstrbndm" assigns parameters for stretch-bend interactions
c     according to the Merck Molecular Force Field (MMFF)
c
c     Note: "stbnt" is the MMFF Stretch-Bend Type for angle "a-b-c",
c     where atom "a" has a smaller class number than atom "c"
c
c     if the BT of a-b = 1, then stbnt = 1
c     if the BT of b-c = 1, then stbnt = 2
c     if both = 1, then stbnt = 3
c     if 4-membered ring, then stbnt = 4
c     if 3-membered ring, then stbnt = 5
c     if 3-membered ring with BT of a-b = 1, then stbnt = 6
c     if 3-membered ring with BT of b-c = 1, then stbnt = 7
c     if 3-membered ring with BT of both = 1, then stbnt = 8
c     if 4-membered ring with BT of a-b = 1, then stbnt = 9
c     if 4-membered ring with BT of b-c = 1, then stbnt = 10
c     if 4-membered ring with BT of both = 1, then stbnt = 11
c     else, if all BT = 0 and no small ring, then stbnt = 0
c
c     literature references:
c
c     T. A. Halgren, "Merck Molecular Force Field. I. Basis, Form,
c     Scope, Parametrization, and Performance of MMFF94", Journal of
c     Computational Chemistry, 17, 490-519 (1995)
c
c     T. A. Halgren, "Merck Molecular Force Field. V. Extension of
c     MMFF94 Using Experimental Data, Additional Computational Data,
c     and Empirical Rules", Journal of Computational Chemistry, 17,
c     616-641 (1995)
c
c
      subroutine kstrbndm
      use sizes
      use angbnd
      use atmlst
      use atomid
      use couple
      use merck
      use potent
      use ring
      use strbnd
      implicit none
      integer i,j,k,l,m
      integer ia,ib,ic
      integer ita,itb,itc
      integer ina,inb,inc
      integer ira,irb,irc
      integer nb1,nb2
      integer stbnt,ab,bc
      logical ring3,ring4
c
c
c     assign stretch-bend parameters for each angle
c
      nstrbnd = 0
      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
c
c     stretch-bend interactions are omitted for linear angles
c
         if (lin(class(ib)) .eq. 0) then
            ita = class (ia)
            itb = class (ib)
            itc = class (ic)
            ina = atomic(ia)
            inb = atomic(ib)
            inc = atomic(ic)
            sbk(1,nstrbnd+1) = 0.0d0
            sbk(2,nstrbnd+1) = 0.0d0
            do k = 1, n12(ib)
               if (i12(k,ib) .eq. ia)  nb1 = bndlist(k,ib)
               if (i12(k,ib) .eq. ic)  nb2 = bndlist(k,ib)
            end do
            stbnt = 0
            ab = 0
            bc = 0
c
c     check if the atoms belong to a single 3- or 4-membered ring
c
            ring3 = .false.
            ring4 = .false.
            do j = 1, nring3
               do k = 1, 3
                  if (ia .eq. iring3(k,j)) then
                     do l = 1, 3
                        if (ib .eq. iring3(l,j)) then
                           do m = 1, 3
                              if (ic .eq. iring3(m,j))
     &                           ring3 = .true.
                           end do
                        end if
                     end do
                  end if
               end do
            end do
            if (.not. ring3) then
               do j = 1, nring4
                  do k = 1, 4
                     if (ia .eq. iring4(k,j)) then
                        do l = 1, 4
                           if (ib .eq. iring4(l,j)) then
                              do m = 1, 4
                                 if (ic .eq. iring4(m,j))
     &                              ring4 = .true.
                              end do
                           end if
                        end do
                     end if
                  end do
               end do
            end if
c
c     determine the MMFF stretch-bend type for the current angle
c
            if (ita .lt. itc) then
               do j = 1, nlignes
                  if (((ia.eq.bt_1(j,1).and.ib.eq.bt_1(j,2)) .or.
     &                 (ib.eq.bt_1(j,1).and.ia.eq.bt_1(j,2)))) then
                     ab = 1
                  end if
                  if (((ic.eq.bt_1(j,1).and.ib.eq.bt_1(j,2)) .or.
     &                 (ib.eq.bt_1(j,1).and.ic.eq.bt_1(j,2)))) then
                     bc = 1
                  end if
               end do
               if (ab.eq.1 .and. bc.eq.0)  stbnt = 1
               if (ab.eq.0 .and. bc.eq.1)  stbnt = 2
               if (ab.eq.1 .and. bc.eq.1)  stbnt = 3
               if (stbnt.eq.0 .AND. ring3) then
                  stbnt = 5
               else if (stbnt.eq.1 .and. ring3) then
                  stbnt = 6
               else if (stbnt.eq.2 .and. ring3) then
                  stbnt = 7
               else if (stbnt.eq.3 .and. ring3) then
                  stbnt = 8
               else if (stbnt.eq.0 .and. ring4) then
                  stbnt = 4
               else if (stbnt.eq.1 .and. ring4) then
                  stbnt = 9
               else if (stbnt.eq.2 .and. ring4) then
                  stbnt = 10
               else if (stbnt.eq.3 .and. ring4) then
                  stbnt = 11
               end if
            else if (ita .gt. itc) then
               do j = 1, nlignes
                  if (((ia.eq.bt_1(j,1).and.ib.eq.bt_1(j,2)) .or.
     &                 (ib.eq.bt_1(j,1).and.ia.eq.bt_1(j,2)))) then
                     ab = 1
                  end if
                  if (((ic.eq.bt_1(j,1).and.ib.eq.bt_1(j,2)) .or.
     &                 (ib.eq.bt_1(j,1).and.ic.eq.bt_1(j,2)))) then
                     bc = 1
                  end if
               end do
               if (ab.eq.1 .and. bc.eq.0)  stbnt = 2
               if (ab.eq.0 .and. bc.eq.1)  stbnt = 1
               if (ab.eq.1 .and. bc.eq.1)  stbnt = 3
               if (stbnt.eq.0 .and. ring3) then
                  stbnt = 5
               else if (stbnt.eq.1 .and. ring3) then
                  stbnt = 6
               else if (stbnt.eq.2 .and. ring3) then
                  stbnt = 7
               else if (stbnt.eq.3 .and. ring3) then
                  stbnt = 8
               else if (stbnt.eq.0 .and. ring4) then
                  stbnt = 4
               else if (stbnt.eq.1 .and. ring4) then
                  stbnt = 9
               else if (stbnt.eq.2 .and. ring4) then
                  stbnt = 10
               else if (stbnt.eq.3 .and. ring4) then
                  stbnt = 11
               end if
            else if (ita .eq. itc) then
               do j = 1, nlignes
                  if (((ic.eq.bt_1(j,1).and.ib.eq.bt_1(j,2)) .or.
     &                 (ib.eq.bt_1(j,1).and.ic.eq.bt_1(j,2)))) then
                     bc = 1
                  end if
                  if (((ia.eq.bt_1(j,1).and.ib.eq.bt_1(j,2)) .or.
     &                 (ib.eq.bt_1(j,1).and.ia.eq.bt_1(j,2)))) then
                     ab = 1
                  end if
               end do
               if (ab.eq.1 .and. bc.eq.0)  stbnt = 1
               if (ab.eq.0 .and. bc.eq.1)  stbnt = 2
               if (ab.eq.1 .and. bc.eq.1)  stbnt = 3
               if (stbnt.eq.0 .and. ring3) then
                  stbnt = 5
               else if (stbnt.eq.1 .and. ring3) then
                  stbnt = 6
               else if (stbnt.eq.2 .and. ring3) then
                  stbnt = 7
               else if (stbnt.eq.3 .and. ring3) then
                  stbnt = 8
               else if (stbnt.eq.0 .and. ring4) then
                  stbnt = 4
               else if (stbnt.eq.1 .and. ring4) then
                  stbnt = 9
               else if (stbnt.eq.2 .and. ring4) then
                  stbnt = 10
               else if (stbnt.eq.3 .and. ring4) then
                  stbnt = 11
               end if
            end if
c
c     find the periodic table row for the atoms in the angle
c
            if (ina .eq. 1)  ira = 0
            if (ina.ge.3 .and. ina.le.10)  ira = 1
            if (ina.ge.11 .and. ina.le.18)  ira = 2
            if (ina.ge.19 .and. ina.le.36)  ira = 3
            if (ina.ge.37 .and. ina.le.54)  ira = 4
            if (inb .eq. 1)  irb = 0
            if (inb.ge.3 .and. inb.le.10)  irb = 1
            if (inb.ge.11 .and. inb.le.18)  irb = 2
            if (inb.ge.19 .and. inb.le.36)  irb = 3
            if (inb.ge.37 .and. inb.le.54)  irb = 4
            if (inc .eq. 1)  irc = 0
            if (inc.ge.3 .and. inc.le.10)  irc = 1
            if (inc.ge.11 .and. inc.le.18)  irc = 2
            if (inc.ge.19 .and. inc.le.36)  irc = 3
            if (inc.ge.37 .and. inc.le.54)  irc = 4
c
c     assign parameters via explicit values or empirical rules
c
            if (stbnt .eq. 11) then
               if ((stbn_abc11(ita,itb,itc).ne.1000.0d0) .and.
     &             (stbn_cba11(ita,itb,itc).ne.1000.0d0)) then
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = stbn_abc11(ita,itb,itc)
                  sbk(2,nstrbnd) = stbn_cba11(ita,itb,itc)
               else
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = defstbn_abc(ira,irb,irc)
                  sbk(2,nstrbnd) = defstbn_cba(ira,irb,irc)
               end if
            else if (stbnt .eq. 10) then
               if ((stbn_abc10(ita,itb,itc).ne.1000.0d0) .and.
     &             (stbn_cba10(ita,itb,itc).ne.1000.0d0)) then
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = stbn_abc10(ita,itb,itc)
                  sbk(2,nstrbnd) = stbn_cba10(ita,itb,itc)
               else
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = defstbn_abc(ira,irb,irc)
                  sbk(2,nstrbnd) = defstbn_cba(ira,irb,irc)
               end if
            else if (stbnt .eq. 9) then
               if ((stbn_abc9(ita,itb,itc).ne.1000.0d0) .and.
     &             (stbn_cba9(ita,itb,itc).ne.1000.0d0)) then
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = stbn_abc9(ita,itb,itc)
                  sbk(2,nstrbnd) = stbn_cba9(ita,itb,itc)
               else
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = defstbn_abc(ira,irb,irc)
                  sbk(2,nstrbnd) = defstbn_cba(ira,irb,irc)
               end if
            else if (stbnt .eq. 8) then
               if ((stbn_abc8(ita,itb,itc).ne.1000.0d0) .and.
     &             (stbn_cba3(ita,itb,itc).ne.1000.0d0)) then
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = stbn_abc8(ita,itb,itc)
                  sbk(2,nstrbnd) = stbn_cba8(ita,itb,itc)
               else
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = defstbn_abc(ira,irb,irc)
                  sbk(2,nstrbnd) = defstbn_cba(ira,irb,irc)
               end if
            else if (stbnt .eq. 7) then
               if ((stbn_abc7(ita,itb,itc).ne.1000.0d0) .and.
     &             (stbn_cba7(ita,itb,itc).ne.1000.0d0)) then
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = stbn_abc7(ita,itb,itc)
                  sbk(2,nstrbnd) = stbn_cba7(ita,itb,itc)
               else
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = defstbn_abc(ira,irb,irc)
                  sbk(2,nstrbnd) = defstbn_cba(ira,irb,irc)
               end if
            else if (stbnt .eq. 6) then
               if ((stbn_abc6(ita,itb,itc).ne.1000.0d0) .and.
     &             (stbn_cba3(ita,itb,itc).ne.1000.0d0)) then
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = stbn_abc6(ita,itb,itc)
                  sbk(2,nstrbnd) = stbn_cba6(ita,itb,itc)
               else
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = defstbn_abc(ira,irb,irc)
                  sbk(2,nstrbnd) = defstbn_cba(ira,irb,irc)
               end if
            else if (stbnt .eq. 5) then
               if (((stbn_abc5(ita,itb,itc).ne.1000.0d0) .and.
     &              (stbn_cba3(ita,itb,itc).ne.1000.0d0))
     &            .or. (ita.eq.22.and.itb.eq.22.and.itc.eq.22)) then
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = stbn_abc5(ita,itb,itc)
                  sbk(2,nstrbnd) = stbn_cba5(ita,itb,itc)
               else
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = defstbn_abc(ira,irb,irc)
                  sbk(2,nstrbnd) = defstbn_cba(ira,irb,irc)
               end if
            else if (stbnt .eq. 4) then
               if ((stbn_abc4(ita,itb,itc).ne.1000.0d0) .and.
     &             (stbn_cba4(ita,itb,itc).ne.1000.0d0)) then
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = stbn_abc4(ita,itb,itc)
                  sbk(2,nstrbnd) = stbn_cba4(ita,itb,itc)
               else
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = defstbn_abc(ira,irb,irc)
                  sbk(2,nstrbnd) = defstbn_cba(ira,irb,irc)
               end if
            else if (stbnt .eq. 3) then
               if ((stbn_abc3(ita,itb,itc).ne.1000.0d0) .and.
     &             (stbn_cba3(ita,itb,itc).ne.1000.0d0)) then
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = stbn_abc3(ita,itb,itc)
                  sbk(2,nstrbnd) = stbn_cba3(ita,itb,itc)
               else
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = defstbn_abc(ira,irb,irc)
                  sbk(2,nstrbnd) = defstbn_cba(ira,irb,irc)
               end if
            else if (stbnt .eq. 2) then
              if ((stbn_abc2(ita,itb,itc).ne.1000.0d0) .and.
     &            (stbn_cba2(ita,itb,itc).ne.1000.0d0)) then
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = stbn_abc2(ita,itb,itc)
                  sbk(2,nstrbnd) = stbn_cba2(ita,itb,itc)
               else
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = defstbn_abc(ira,irb,irc)
                  sbk(2,nstrbnd) = defstbn_cba(ira,irb,irc)
               end if
            else if (stbnt .eq. 1) then
               if ((stbn_abc1(ita,itb,itc).ne.1000.0d0) .and.
     &             (stbn_cba1(ita,itb,itc).ne.1000.0d0)) then
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = stbn_abc1(ita,itb,itc)
                  sbk(2,nstrbnd) = stbn_cba1(ita,itb,itc)
               else
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = defstbn_abc(ira,irb,irc)
                  sbk(2,nstrbnd) = defstbn_cba(ira,irb,irc)
               end if
            else if (stbnt .eq. 0) then
               if (((stbn_abc(ita,itb,itc) .ne. 1000.0d0) .and.
     &              (stbn_cba(ita,itb,itc) .ne. 1000.0d0))
     &            .or. (ita.eq.12.AND.itb.eq.20.AND.itc.eq.20)
     &            .or. (ita.eq.20.AND.itb.eq.20.AND.itc.eq.12)) then
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = stbn_abc(ita,itb,itc)
                  sbk(2,nstrbnd) = stbn_cba(ita,itb,itc)
               else
                  nstrbnd = nstrbnd + 1
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nb1
                  isb(3,nstrbnd) = nb2
                  sbk(1,nstrbnd) = defstbn_abc(ira,irb,irc)
                  sbk(2,nstrbnd) = defstbn_cba(ira,irb,irc)
               end if
            end if
         end if
      end do
c
c     turn off the stretch-bend potential if it is not used
c
      if (nstrbnd .eq. 0)  use_strbnd = .false.
      return
      end
