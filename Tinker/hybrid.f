c
c
c     ###############################################################
c     ##  COPYRIGHT (C) 1991 by Shawn Huston & Jay William Ponder  ##
c     ##                    All Rights Reserved                    ##
c     ###############################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine hybrid  --  set parameters for hybrid system  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "hybrid" constructs the hybrid hamiltonian for a specified
c     initial state, final state and mutation parameter "lambda"
c
c
      subroutine hybrid
      use sizes
      use iounit
      use mutant
      implicit none
c
c
c     set the potential energy parameters for hybrid atoms
c
      if (nmut .ne. 0) then
         write (iout,10)  lambda
   10    format (/,' Lambda Coupling Parameter for FEP :',f12.3)
         call hatom
         call hbond
         call hangle
         call hstrbnd
         call himptor
         call htors
         call hstrtor
         call hvdw
         call hcharge
         call hdipole
      end if
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine hatom  --  assign hybrid atom parameters  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "hatom" assigns a new atom type to each hybrid site
c
c
      subroutine hatom
      use sizes
      use atomid
      use atoms
      use inform
      use iounit
      use katoms
      use mutant
      implicit none
      integer i,k,ntype
      integer it,it0,it1
c
c
c     find the total number of atom types currently used;
c     exclude the "HYB" types so that they can be reused
c
      do i = 1, maxtyp
         if (symbol(i).eq.'   ' .or. symbol(i).eq.'HYB') then
            ntype = i - 1
            goto 10
         end if
      end do
   10 continue
c
c     stop if there are too many atom types required
c
      if (maxtyp .lt. ntype+nmut) then
         abort = .true.
         write (iout,20)
   20    format (' HATOM  --  Too many Sites to be Altered;',
     &              ' Increase MAXTYP')
      end if
c
c     create a new atom type for each of the hybrid atoms
c
      do i = 1, nmut
         k = imut(i)
         it = ntype + i
         it0 = type0(i)
         it1 = type1(i)
         symbol(it) = 'HYB'
         atmnum(it) = 0
         weight(it) = lambda*weight(it1) + (1.0d0-lambda)*weight(it0)
         ligand(it) = 0
         describe(it) = 'Hybrid Atom Type        '
         type(k) = it
         name(k) = symbol(it)
         atomic(k) = atmnum(it)
         mass(k) = weight(it)
         valence(k) = ligand(it)
         story(k) = describe(it)
      end do
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine hbond  --  find hybrid bond parameters  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "hbond" constructs hybrid bond stretch parameters given
c     an initial state, final state and "lambda" value
c
c
      subroutine hbond
      use sizes
      use atomid
      use atoms
      use bndstr
      use iounit
      use inform
      use kbonds
      use mutant
      implicit none
      integer i,j,k
      integer ia,ib
      integer ita,itb
      integer size
      real*8 bk0,bk1
      real*8 bl0,bl1
      logical header
      character*4 pa,pb
      character*8 pt
c
c
c     assign the hybrid parameters for individual bonds
c
      header = .true.
      do i = 1, nbond
         ia = ibnd(1,i)
         ib = ibnd(2,i)
         if (mut(ia) .or. mut(ib)) then
            ita = class(ia)
            itb = class(ib)
c
c     find the bond parameters for the initial state
c
            do j = 1, nmut
               k = imut(j)
               if (k .eq. ia)  ita = class0(j)
               if (k .eq. ib)  itb = class0(j)
            end do
            size = 4
            call numeral (ita,pa,size)
            call numeral (itb,pb,size)
            if (ita .le. itb) then
               pt = pa//pb
            else
               pt = pb//pa
            end if
            bk0 = 0.0d0
            bl0 = 0.0d0
            do j = 1, maxnb
               if (kb(j) .eq. pt) then
                  bk0 = bcon(j)
                  bl0 = blen(j)
                  goto 10
               end if
            end do
   10       continue
c
c     find the bond parameters for the final state
c
            do j = 1, nmut
               k = imut(j)
               if (k .eq. ia)  ita = class1(j)
               if (k .eq. ib)  itb = class1(j)
            end do
            size = 4
            call numeral (ita,pa,size)
            call numeral (itb,pb,size)
            if (ita .le. itb) then
               pt = pa//pb
            else
               pt = pb//pa
            end if
            bk1 = 0.0d0
            bl1 = 0.0d0
            do j = 1, maxnb
               if (kb(j) .eq. pt) then
                  bk1 = bcon(j)
                  bl1 = blen(j)
                  goto 20
               end if
            end do
   20       continue
c
c     form the hybrid parameters for the current bond
c
            if (bl0 .eq. 0.0d0)  bl0 = bl1
            if (bl1 .eq. 0.0d0)  bl1 = bl0
            bk(i) = lambda*bk1 + (1.0d0-lambda)*bk0
            bl(i) = lambda*bl1 + (1.0d0-lambda)*bl0
            if (verbose) then
               if (header) then
                  header = .false.
                  write (iout,30)
   30             format (/,' Hybrid Bond Stretching Parameters :',
     &                    //,6x,'Atom Numbers',9x,'KS',7x,'Length',/)
               end if
               write (iout,40)  ia,ib,bk(i),bl(i)
   40          format (6x,2i5,f14.3,f12.4)
            end if
         end if
      end do
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine hangle  --  find hybrid angle parameters  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "hangle" constructs hybrid angle bending parameters given
c     an initial state, final state and "lambda" value
c
c
      subroutine hangle
      use sizes
      use angbnd
      use atomid
      use atoms
      use iounit
      use inform
      use kangs
      use mutant
      implicit none
      integer i,j,k,size
      integer ia,ib,ic
      integer ita,itb,itc
      real*8 ak0,ak1
      real*8 anat0,anat1
      logical header
      character*4 pa,pb,pc
      character*12 pt
c
c
c     assign the hybrid parameters for individual angles
c
      header = .true.
      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         if (mut(ia) .or. mut(ib) .or. mut(ic)) then
            ita = class(ia)
            itb = class(ib)
            itc = class(ic)
c
c     find the angle parameters for the initial state
c
            do j = 1, nmut
               k = imut(j)
               if (k .eq. ia)  ita = class0(j)
               if (k .eq. ib)  itb = class0(j)
               if (k .eq. ic)  itc = class0(j)
            end do
            size = 4
            call numeral (ita,pa,size)
            call numeral (itb,pb,size)
            call numeral (itc,pc,size)
            if (ita .le. itc) then
               pt = pa//pb//pc
            else
               pt = pc//pb//pa
            end if
            ak0 = 0.0d0
            anat0 = 0.0d0
            do j = 1, maxna
               if (ka(j) .eq. pt) then
                  ak0 = acon(j)
                  anat0 = ang(1,j)
                  goto 10
               end if
            end do
   10       continue
c
c     find the angle parameters for the final state
c
            do j = 1, nmut
               k = imut(j)
               if (k .eq. ia)  ita = class1(j)
               if (k .eq. ib)  itb = class1(j)
               if (k .eq. ic)  itc = class1(j)
            end do
            size = 4
            call numeral (ita,pa,3)
            call numeral (itb,pb,3)
            call numeral (itc,pc,3)
            if (ita .le. itc) then
               pt = pa//pb//pc
            else
               pt = pc//pb//pa
            end if
            ak1 = 0.0d0
            anat1 = 0.0d0
            do j = 1, maxna
               if (ka(j) .eq. pt) then
                  ak1 = acon(j)
                  anat1 = ang(1,j)
                  goto 20
               end if
            end do
   20       continue
c
c     form the hybrid parameters for the current angle
c
            if (anat0 .eq. 0.0d0)  anat0 = anat1
            if (anat1 .eq. 0.0d0)  anat1 = anat0
            ak(i) = lambda*ak1 + (1.0d0-lambda)*ak0
            anat(i) = lambda*anat1 + (1.0d0-lambda)*anat0
            if (verbose) then
               if (header) then
                  header = .false.
                  write (iout,30)
   30             format (/,' Hybrid Angle Bending Parameters :',
     &                    //,6x,'Atom Numbers',9x,'KB',8x,'Angle',/)
               end if
               write (iout,40)  ia,ib,ic,ak(i),anat(i)
   40          format (3x,3i5,2f12.3)
            end if
         end if
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine hstrbnd  --  hybrid stretch-bend parameters  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "hstrbnd" constructs hybrid stretch-bend parameters given
c     an initial state, final state and "lambda" value
c
c
      subroutine hstrbnd
      use sizes
      use angbnd
      use atmlst
      use atomid
      use atoms
      use couple
      use iounit
      use inform
      use katoms
      use kstbnd
      use mutant
      use strbnd
      implicit none
      integer i,j,k,size
      integer ia,ib,ic
      integer ita,itb,itc
      integer nba,nbc
      real*8 sbk0(2),sbk1(2)
      logical header,used
      character*4 pa,pb,pc
      character*12 pt
c
c
c     assign hybrid parameters for the stretch-bend sites
c
      header = .true.
      do i = 1, nangle
         used = .false.
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         if (mut(ia) .or. mut(ib) .or. mut(ic)) then
            ita = class(ia)
            itb = class(ib)
            itc = class(ic)
            do j = 1, n12(ib)
               if (i12(j,ib) .eq. ia)  nba = bndlist(j,ib)
               if (i12(j,ib) .eq. ic)  nbc = bndlist(j,ib)
            end do
c
c     find the stretch-bend parameters for the initial state
c
            do j = 1, nmut
               k = imut(j)
               if (k .eq. ia)  ita = class0(j)
               if (k .eq. ib)  itb = class0(j)
               if (k .eq. ic)  itc = class0(j)
            end do
            size = 4
            call numeral (ita,pa,size)
            call numeral (itb,pb,size)
            call numeral (itc,pc,size)
            if (ita .le. itc) then
               pt = pa//pb//pc
            else
               pt = pc//pb//pa
            end if
            sbk0(1) = 0.0d0
            sbk0(2) = 0.0d0
            do j = 1, maxnsb
               if (ksb(j) .eq. pt) then
                  used = .true.
                  if (ita .le. itc) then
                     sbk0(1) = stbn(1,j)
                     sbk0(2) = stbn(2,j)
                  else
                     sbk0(1) = stbn(2,j)
                     sbk0(2) = stbn(1,j)
                  end if
                  goto 10
               end if
            end do
   10       continue
c
c     find the stretch-bend parameters for the final state
c
            do j = 1, nmut
               k = imut(j)
               if (k .eq. ia)  ita = class1(j)
               if (k .eq. ib)  itb = class1(j)
               if (k .eq. ic)  itc = class1(j)
            end do
            size = 4
            call numeral (ita,pa,size)
            call numeral (itb,pb,size)
            call numeral (itc,pc,size)
            if (ita .le. itc) then
               pt = pa//pb//pc
            else
               pt = pc//pb//pa
            end if
            sbk1(1) = 0.0d0
            sbk1(2) = 0.0d0
            do j = 1, maxnsb
               if (ksb(j) .eq. pt) then
                  used = .true.
                  if (ita .le. itc) then
                     sbk1(1) = stbn(1,j)
                     sbk1(2) = stbn(2,j)
                  else
                     sbk1(1) = stbn(2,j)
                     sbk1(2) = stbn(1,j)
                  end if
                  goto 20
               end if
            end do
   20       continue
c
c     form hybrid parameters for the current stretch-bend
c
            if (used) then
               nstrbnd = nstrbnd + 1
               k = nstrbnd
               isb(1,k) = i
               isb(2,k) = nba
               isb(3,k) = nbc
               sbk(1,k) = lambda*sbk1(1) + (1.0d0-lambda)*sbk0(1)
               sbk(2,k) = lambda*sbk1(2) + (1.0d0-lambda)*sbk0(2)
               if (verbose) then
                  if (header) then
                     header = .false.
                     write (iout,30)
   30                format (/,' Hybrid Stretch-Bend Parameters :',
     &                       //,6x,'Atom Numbers',8x,'KSB 1',
     &                          7x,'KSB 2',/)
                  end if
                  write (iout,40)  ia,ib,ic,sbk(1,i),sbk(2,i)
   40             format (3x,3i5,2f12.3)
               end if
            end if
         end if
      end do
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine himptor  --  find hybrid improper torsions  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "himptor" constructs hybrid improper torsional parameters
c     given an initial state, final state and "lambda" value
c
c     note this version does not handle multiple parameters at
c     a single trigonal site
c
c
      subroutine himptor
      use sizes
      use atomid
      use atoms
      use couple
      use iounit
      use inform
      use imptor
      use kitors
      use math
      use mutant
      implicit none
      integer i,j,k
      integer ia,ib,ic,id
      integer ita,itb,itc,itd
      integer nti,size
      real*8 angle,symm
      real*8 v1_0,v2_0,v3_0
      real*8 s1_0,s2_0,s3_0
      real*8 v1_1,v2_1,v3_1
      real*8 s1_1,s2_1,s3_1
      logical header,used
      character*4 pa,pb,pc,pd
      character*4 zeros
      character*16 blank
      character*16 pt0,pt(6)
c
c
c     construct hybrid improper torsion parameters
c
      blank = '                '
      zeros = '0000'
      header = .true.
c
c     determine the total number of forcefield parameters
c
      nti = maxnti
      do i = maxnti, 1, -1
         if (kti(i) .eq. blank)  nti = i - 1
      end do
c
c     construct hybrid improper torsion parameters
c
      do i = 1, n
         if (n12(i) .eq. 3) then
            used = .false.
            ia = i12(1,i)
            ib = i12(2,i)
            ic = i
            id = i12(3,i)
            if (mut(ia) .or. mut(ib) .or. mut(ic) .or. mut(id)) then
               ita = class(ia)
               itb = class(ib)
               itc = class(ic)
               itd = class(id)
c
c     find improper torsion parameters for the initial state
c
               do j = 1, nmut
                  k = imut(j)
                  if (k .eq. ia)  ita = class0(j)
                  if (k .eq. ib)  itb = class0(j)
                  if (k .eq. ic)  itc = class0(j)
                  if (k .eq. id)  itd = class0(j)
               end do
               size = 4
               call numeral (ita,pa,size)
               call numeral (itb,pb,size)
               call numeral (itc,pc,size)
               call numeral (itd,pd,size)
               pt(1) = pa//pb//pc//pd
               pt(2) = pb//pa//pc//pd
               pt(3) = pa//pd//pc//pb
               pt(4) = pd//pa//pc//pb
               pt(5) = pb//pd//pc//pa
               pt(6) = pd//pb//pc//pa
               pt0 = zeros//zeros//pc//zeros
               symm = 1.0d0
               if (pa.eq.pb .or. pa.eq.pd .or. pb.eq.pd)  symm = 2.0d0
               if (pa.eq.pb .and. pa.eq.pd .and. pb.eq.pd)  symm = 6.0d0
               v1_0 = 0.0d0
               s1_0 = 0.0d0
               v2_0 = 0.0d0
               s2_0 = 0.0d0
               v3_0 = 0.0d0
               s3_0 = 0.0d0
               do j = 1, nti
                  if (kti(j)(9:12) .eq. pc) then
                     do k = 1, 6
                        if (kti(j) .eq. pt(k)) then
                           used = .true.
                           v1_0 = ti1(1,j) / symm
                           s1_0 = ti1(2,j)
                           v2_0 = ti2(1,j) / symm
                           s2_0 = ti2(2,j)
                           v3_0 = ti3(1,j) / symm
                           s3_0 = ti3(2,j)
                           goto 10
                        end if
                     end do
                  end if
               end do
               do j = 1, nti
                  if (kti(j) .eq. pt0) then
                     used = .true.
                     v1_0 = ti1(1,j) / symm
                     s1_0 = ti1(2,j)
                     v2_0 = ti2(1,j) / symm
                     s2_0 = ti2(2,j)
                     v3_0 = ti3(1,j) / symm
                     s3_0 = ti3(2,j)
                     goto 10
                  end if
               end do
   10          continue
c
c     find improper torsion parameters for the final state
c
               do j = 1, nmut
                  k = imut(j)
                  if (k .eq. ia)  ita = class1(j)
                  if (k .eq. ib)  itb = class1(j)
                  if (k .eq. ic)  itc = class1(j)
                  if (k .eq. id)  itd = class1(j)
               end do
               size = 4
               call numeral (ita,pa,size)
               call numeral (itb,pb,size)
               call numeral (itc,pc,size)
               call numeral (itd,pd,size)
               pt(1) = pa//pb//pc//pd
               pt(2) = pb//pa//pc//pd
               pt(3) = pa//pd//pc//pb
               pt(4) = pd//pa//pc//pb
               pt(5) = pb//pd//pc//pa
               pt(6) = pd//pb//pc//pa
               pt0 = zeros//zeros//pc//zeros
               symm = 1.0d0
               if (pa.eq.pb .or. pa.eq.pd .or. pb.eq.pd)  symm = 2.0d0
               if (pa.eq.pb .and. pa.eq.pd .and. pb.eq.pd)  symm = 6.0d0
               v1_1 = 0.0d0
               s1_1 = 0.0d0
               v2_1 = 0.0d0
               s2_1 = 0.0d0
               v3_1 = 0.0d0
               s3_1 = 0.0d0
               do j = 1, nti
                  if (kti(j)(9:12) .eq. pc) then
                     do k = 1, 6
                        if (kti(j) .eq. pt(k)) then
                           used = .true.
                           v1_1 = ti1(1,j) / symm
                           s1_1 = ti1(2,j)
                           v2_1 = ti2(1,j) / symm
                           s2_1 = ti2(2,j)
                           v3_1 = ti3(1,j) / symm
                           s3_1 = ti3(2,j)
                           goto 20
                        end if
                     end do
                  end if
               end do
               do j = 1, nti
                  if (kti(j) .eq. pt0) then
                     used = .true.
                     v1_1 = ti1(1,j) / symm
                     s1_1 = ti1(2,j)
                     v2_1 = ti2(1,j) / symm
                     s2_1 = ti2(2,j)
                     v3_1 = ti3(1,j) / symm
                     s3_1 = ti3(2,j)
                     goto 20
                  end if
               end do
   20          continue
c
c     form hybrid parameters for the current improper torsion
c
               if (used) then
                  do j = 1, nitors
                     if (iitors(3,j) .eq. ic) then
                        k = j
                        goto 30
                     end if
                  end do
                  nitors = nitors + 1
                  k = nitors
                  iitors(1,k) = ia
                  iitors(2,k) = ib
                  iitors(3,k) = ic
                  iitors(4,k) = id
   30             continue
                  if (s1_0 .eq. 0.0d0)  s1_0 = s1_1
                  if (s2_0 .eq. 0.0d0)  s2_0 = s2_1
                  if (s3_0 .eq. 0.0d0)  s3_0 = s3_1
                  if (s1_1 .eq. 0.0d0)  s1_1 = s1_0
                  if (s2_1 .eq. 0.0d0)  s2_1 = s2_0
                  if (s3_1 .eq. 0.0d0)  s3_1 = s3_0
                  itors1(1,k) = lambda*v1_1 + (1.0d0-lambda)*v1_0
                  itors1(2,k) = lambda*s1_1 + (1.0d0-lambda)*s1_0
                  angle = itors1(2,k) / radian
                  itors1(3,k) = cos(angle)
                  itors1(4,k) = sin(angle)
                  itors2(1,k) = lambda*v2_1 + (1.0d0-lambda)*v2_0
                  itors2(2,k) = lambda*s2_1 + (1.0d0-lambda)*s2_0
                  angle = itors2(2,k) / radian
                  itors2(3,k) = cos(angle)
                  itors2(4,k) = sin(angle)
                  itors3(1,k) = lambda*v3_1 + (1.0d0-lambda)*v3_0
                  itors3(2,k) = lambda*s3_1 + (1.0d0-lambda)*s3_0
                  angle = itors3(2,k) / radian
                  itors3(3,k) = cos(angle)
                  itors3(4,k) = sin(angle)
                  if (verbose) then
                     if (header) then
                        header = .false.
                        write (iout,40)
   40                   format (/,' Hybrid Improper Torsional',
     &                             ' Parameters :',
     &                          //,6x,'Atom Numbers',16x,'KIT1',
     &                             13x,'KIT2',13x,'KIT3',/)
                     end if
                     ia = iitors(1,i)
                     ib = iitors(2,i)
                     ic = iitors(3,i)
                     id = iitors(4,i)
                     write (iout,50)  ia,ib,ic,id,itors1(1,k),
     &                                itors1(2,k),itors2(1,k),
     &                                itors2(2,k),itors3(1,k),
     &                                itors3(2,k)
   50                format (1x,4i5,4x,3(f10.4,f7.1))
                  end if
               end if
            end if
         end if
      end do
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine htors  --  find hybrid torsion parameters  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "htors" constructs hybrid torsional parameters for a given
c     initial state, final state and "lambda" value
c
c
      subroutine htors
      use sizes
      use atomid
      use atoms
      use inform
      use iounit
      use ktorsn
      use math
      use mutant
      use tors
      implicit none
      integer i,j,k,size
      integer ia,ib,ic,id
      integer ita,itb,itc,itd
      real*8 angle
      real*8 v1_0,v2_0,v3_0
      real*8 v4_0,v5_0,v6_0
      real*8 s1_0,s2_0,s3_0
      real*8 s4_0,s5_0,s6_0
      real*8 v1_1,v2_1,v3_1
      real*8 v4_1,v5_1,v6_1
      real*8 s1_1,s2_1,s3_1
      real*8 s4_1,s5_1,s6_1
      logical header
      character*4 pa,pb,pc,pd
      character*4 zeros
      character*16 pt,pt0
c
c
c     construct hybrid torsional parameters
c
      zeros = '0000'
      header = .true.
      do i = 1, ntors
         ia = itors(1,i)
         ib = itors(2,i)
         ic = itors(3,i)
         id = itors(4,i)
         if (mut(ia) .or. mut(ib) .or. mut(ic) .or. mut(id)) then
            ita = class(ia)
            itb = class(ib)
            itc = class(ic)
            itd = class(id)
c
c     find the torsion parameters for the initial state
c
            do j = 1, nmut
               k = imut(j)
               if (k .eq. ia)  ita = class0(j)
               if (k .eq. ib)  itb = class0(j)
               if (k .eq. ic)  itc = class0(j)
               if (k .eq. id)  itd = class0(j)
            end do
            size = 4
            call numeral (ita,pa,size)
            call numeral (itb,pb,size)
            call numeral (itc,pc,size)
            call numeral (itd,pd,size)
            if (itb .lt. itc) then
               pt = pa//pb//pc//pd
            else if (itc .lt. itb) then
               pt = pd//pc//pb//pa
            else if (ita .le. itd) then
               pt = pa//pb//pc//pd
            else if (itd .lt. ita) then
               pt = pd//pc//pb//pa
            end if
            pt0 = zeros//pt(5:12)//zeros
            v1_0 = 0.0d0
            s1_0 = 0.0d0
            v2_0 = 0.0d0
            s2_0 = 0.0d0
            v3_0 = 0.0d0
            s3_0 = 0.0d0
            v4_0 = 0.0d0
            s4_0 = 0.0d0
            v5_0 = 0.0d0
            s5_0 = 0.0d0
            v6_0 = 0.0d0
            s6_0 = 0.0d0
            do j = 1, maxnt
               if (kt(j) .eq. pt) then
                  v1_0 = t1(1,j)
                  s1_0 = t1(2,j)
                  v2_0 = t2(1,j)
                  s2_0 = t2(2,j)
                  v3_0 = t3(1,j)
                  s3_0 = t3(2,j)
                  v4_0 = t4(1,j)
                  s4_0 = t4(2,j)
                  v5_0 = t5(1,j)
                  s5_0 = t5(2,j)
                  v6_0 = t6(1,j)
                  s6_0 = t6(2,j)
                  goto 10
               end if
            end do
            do j = 1, maxnt
               if (kt(j) .eq. pt0) then
                  v1_0 = t1(1,j)
                  s1_0 = t1(2,j)
                  v2_0 = t2(1,j)
                  s2_0 = t2(2,j)
                  v3_0 = t3(1,j)
                  s3_0 = t3(2,j)
                  v4_0 = t4(1,j)
                  s4_0 = t4(2,j)
                  v5_0 = t5(1,j)
                  s5_0 = t5(2,j)
                  v6_0 = t6(1,j)
                  s6_0 = t6(2,j)
                  goto 10
               end if
            end do
   10       continue
c
c     find the torsion parameters for the final state
c
            do j = 1, nmut
               k = imut(j)
               if (k .eq. ia)  ita = class1(j)
               if (k .eq. ib)  itb = class1(j)
               if (k .eq. ic)  itc = class1(j)
               if (k .eq. id)  itd = class1(j)
            end do
            size = 4
            call numeral (ita,pa,size)
            call numeral (itb,pb,size)
            call numeral (itc,pc,size)
            call numeral (itd,pd,size)
            if (itb .lt. itc) then
               pt = pa//pb//pc//pd
            else if (itc .lt. itb) then
               pt = pd//pc//pb//pa
            else if (ita .le. itd) then
               pt = pa//pb//pc//pd
            else if (itd .lt. ita) then
               pt = pd//pc//pb//pa
            end if
            pt0 = zeros//pt(5:12)//zeros
            v1_1 = 0.0d0
            s1_1 = 0.0d0
            v2_1 = 0.0d0
            s2_1 = 0.0d0
            v3_1 = 0.0d0
            s3_1 = 0.0d0
            v4_1 = 0.0d0
            s4_1 = 0.0d0
            v5_1 = 0.0d0
            s5_1 = 0.0d0
            v6_1 = 0.0d0
            s6_1 = 0.0d0
            do j = 1, maxnt
               if (kt(j) .eq. pt) then
                  v1_1 = t1(1,j)
                  s1_1 = t1(2,j)
                  v2_1 = t2(1,j)
                  s2_1 = t2(2,j)
                  v3_1 = t3(1,j)
                  s3_1 = t3(2,j)
                  v4_1 = t4(1,j)
                  s4_1 = t4(2,j)
                  v5_1 = t5(1,j)
                  s5_1 = t5(2,j)
                  v6_1 = t6(1,j)
                  s6_1 = t6(2,j)
                  goto 20
               end if
            end do
            do j = 1, maxnt
               if (kt(j) .eq. pt0) then
                  v1_1 = t1(1,j)
                  s1_1 = t1(2,j)
                  v2_1 = t2(1,j)
                  s2_1 = t2(2,j)
                  v3_1 = t3(1,j)
                  s3_1 = t3(2,j)
                  v4_1 = t4(1,j)
                  s4_1 = t4(2,j)
                  v5_1 = t5(1,j)
                  s5_1 = t5(2,j)
                  v6_1 = t6(1,j)
                  s6_1 = t6(2,j)
                  goto 20
               end if
            end do
   20       continue
c
c     form the hybrid parameters for the current torsion
c
            if (s1_0 .eq. 0.0d0)  s1_0 = s1_1
            if (s2_0 .eq. 0.0d0)  s2_0 = s2_1
            if (s3_0 .eq. 0.0d0)  s3_0 = s3_1
            if (s4_0 .eq. 0.0d0)  s4_0 = s4_1
            if (s5_0 .eq. 0.0d0)  s5_0 = s5_1
            if (s6_0 .eq. 0.0d0)  s6_0 = s6_1
            if (s1_1 .eq. 0.0d0)  s1_1 = s1_0
            if (s2_1 .eq. 0.0d0)  s2_1 = s2_0
            if (s3_1 .eq. 0.0d0)  s3_1 = s3_0
            if (s4_1 .eq. 0.0d0)  s4_1 = s4_0
            if (s5_1 .eq. 0.0d0)  s5_1 = s5_0
            if (s6_1 .eq. 0.0d0)  s6_1 = s6_0
            tors1(1,i) = lambda*v1_1 + (1.0d0-lambda)*v1_0
            tors1(2,i) = lambda*s1_1 + (1.0d0-lambda)*s1_0
            angle = tors1(2,i) / radian
            tors1(3,i) = cos(angle)
            tors1(4,i) = sin(angle)
            tors2(1,i) = lambda*v2_1 + (1.0d0-lambda)*v2_0
            tors2(2,i) = lambda*s2_1 + (1.0d0-lambda)*s2_0
            angle = tors2(2,i) / radian
            tors2(3,i) = cos(angle)
            tors2(4,i) = sin(angle)
            tors3(1,i) = lambda*v3_1 + (1.0d0-lambda)*v3_0
            tors3(2,i) = lambda*s3_1 + (1.0d0-lambda)*s3_0
            angle = tors3(2,i) / radian
            tors3(3,i) = cos(angle)
            tors3(4,i) = sin(angle)
            tors4(1,i) = lambda*v4_1 + (1.0d0-lambda)*v4_0
            tors4(2,i) = lambda*s4_1 + (1.0d0-lambda)*s4_0
            angle = tors4(2,i) / radian
            tors4(3,i) = cos(angle)
            tors4(4,i) = sin(angle)
            tors5(1,i) = lambda*v5_1 + (1.0d0-lambda)*v5_0
            tors5(2,i) = lambda*s5_1 + (1.0d0-lambda)*s5_0
            angle = tors5(2,i) / radian
            tors5(3,i) = cos(angle)
            tors5(4,i) = sin(angle)
            tors6(1,i) = lambda*v6_1 + (1.0d0-lambda)*v6_0
            tors6(2,i) = lambda*s6_1 + (1.0d0-lambda)*s6_0
            angle = tors6(2,i) / radian
            tors6(3,i) = cos(angle)
            tors6(4,i) = sin(angle)
            if (verbose) then
               if (header) then
                  header = .false.
                  write (iout,30)
   30             format (/,' Hybrid Torsional Parameters :',
     &                    //,5x,'Atom Numbers',6x,'KT1',7x,'KT2',
     &                       7x,'KT3',7x,'KT4',7x,'KT5',7x,'KT6',/)
               end if
               write (iout,40)  ia,ib,ic,id,
     &                          tors1(1,i),nint(tors1(2,i)),
     &                          tors2(1,i),nint(tors2(2,i)),
     &                          tors3(1,i),nint(tors3(2,i)),
     &                          tors4(1,i),nint(tors4(2,i)),
     &                          tors5(1,i),nint(tors5(2,i)),
     &                          tors6(1,i),nint(tors6(2,i))
   40          format (1x,4i4,1x,6(f6.2,i4))
            end if
         end if
      end do
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine hstrtor  --  hybrid stretch-torsion terms  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "hstrtor" constructs hybrid stretch-torsion parameters
c     given an initial state, final state and "lambda" value
c
c
      subroutine hstrtor
      use sizes
      use atmlst
      use atomid
      use atoms
      use couple
      use inform
      use iounit
      use ksttor
      use mutant
      use strtor
      use tors
      implicit none
      integer i,j,k,size
      integer ia,ib,ic,id
      integer ita,itb,itc,itd
      real*8 kst0(3),kst1(3)
      logical header
      character*4 pa,pb,pc,pd
      character*4 zeros
      character*16 pt,pt0
c
c
c     assign hybrid parameters for the stretch-torsion sites
c
      zeros = '0000'
      header = .true.
      do i = 1, ntors
         ia = itors(1,i)
         ib = itors(2,i)
         ic = itors(3,i)
         id = itors(4,i)
         if (mut(ia) .or. mut(ib) .or. mut(ic) .or. mut(id)) then
            ita = class(ia)
            itb = class(ib)
            itc = class(ic)
            itd = class(id)
c
c     find the stretch-torsion parameters for the initial state
c
            do j = 1, nmut
               k = imut(j)
               if (k .eq. ia)  ita = class0(j)
               if (k .eq. ib)  itb = class0(j)
               if (k .eq. ic)  itc = class0(j)
               if (k .eq. id)  itd = class0(j)
            end do
            size = 4
            call numeral (ita,pa,size)
            call numeral (itb,pb,size)
            call numeral (itc,pc,size)
            call numeral (itd,pd,size)
            if (itb .lt. itc) then
               pt = pa//pb//pc//pd
            else if (itc .lt. itb) then
               pt = pd//pc//pb//pa
            else if (ita .le. itd) then
               pt = pa//pb//pc//pd
            else if (itd .lt. ita) then
               pt = pd//pc//pb//pa
            end if
            pt0 = zeros//pt(5:12)//zeros
            do k = 1, 3
               kst0(k) = 0.0d0
            end do
            do j = 1, maxnbt
               if (kbt(j) .eq. pt) then
                  do k = 1, 3
                     kst0(k) = btcon(k,j)
                  end do
                  goto 10
               end if
            end do
            do j = 1, maxnbt
               if (kbt(j) .eq. pt0) then
                  do k = 1, 3
                     kst0(k) = btcon(k,j)
                  end do
                  goto 10
               end if
            end do
   10       continue
c
c     find the stretch-torsion parameters for the final state
c
            do j = 1, nmut
               k = imut(j)
               if (k .eq. ia)  ita = class0(j)
               if (k .eq. ib)  itb = class0(j)
               if (k .eq. ic)  itc = class0(j)
               if (k .eq. id)  itd = class0(j)
            end do
            size = 4
            call numeral (ita,pa,size)
            call numeral (itb,pb,size)
            call numeral (itc,pc,size)
            call numeral (itd,pd,size)
            if (itb .lt. itc) then
               pt = pa//pb//pc//pd
            else if (itc .lt. itb) then
               pt = pd//pc//pb//pa
            else if (ita .le. itd) then
               pt = pa//pb//pc//pd
            else if (itd .lt. ita) then
               pt = pd//pc//pb//pa
            end if
            pt0 = zeros//pt(5:12)//zeros
            do k = 1, 3
               kst1(k) = 0.0d0
            end do
            do j = 1, maxnbt
               if (kbt(j) .eq. pt) then
                  do k = 1, 3
                     kst1(k) = btcon(k,j)
                  end do
                  goto 20
               end if
            end do
            do j = 1, maxnbt
               if (kbt(j) .eq. pt0) then
                  do k = 1, 3
                     kst1(k) = btcon(k,j)
                  end do
                  goto 20
               end if
            end do
   20       continue
c
c     form hybrid parameters for the current stretch-torsion
c
            do j = 1, 3
               kst(j,i) = lambda*kst1(j) + (1.0d0-lambda)*kst0(j)
            end do
            if (kst(1,i).eq.0.0d0 .and. kst(2,i).eq.0.0d0
     &                  .and. kst(3,i).eq.0.0d0) then
               if (ist(1,i) .ne. 0) then
                  nstrtor = nstrtor - 1
                  ist(1,i) = 0
               end if
            else
               if (ist(1,i) .ne. i) then
                  nstrtor = nstrtor + 1
                  ist(1,i) = i
                  do j = 1, n12(ib)
                     if (i12(j,ib) .eq. ic) then
                        ist(2,i) = bndlist(j,ib)
                        goto 30
                     end if
                  end do
   30             continue
               end if
               if (verbose) then
                  if (header) then
                     header = .false.
                     write (iout,40)
   40                format (/,' Hybrid Stretch-Torsion Parameters :',
     &                       //,6x,'Atom Numbers',13x,'KST1',8x,'KST2',
     &                          8x,'KST3',/)
                  end if
                  write (iout,50)  ia,ib,ic,id,(kst(j,i),j=1,3)
   50             format (3x,4i5,3f12.3)
               end if
            end if
         end if
      end do
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine hvdw  --  hybrid van der Waals parameters  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "hvdw" constructs hybrid van der Waals  parameters given
c     an initial state, final state and "lambda" value
c
c
      subroutine hvdw
      use sizes
      use atomid
      use atoms
      use inform
      use iounit
      use kvdws
      use math
      use mutant
      use vdw
      use vdwpot
      implicit none
      integer i,j,k
      integer it,it0,it1
      real*8 radius,rd,ep
      real*8, allocatable :: srad(:)
      real*8, allocatable :: seps(:)
      logical header
c
c
c     assign the hybrid van der Waals parameters
c
      do j = 1, nmut
         i = imut(j)
         it = class(i)
         it0 = class0(j)
         it1 = class1(j)
         rad(it) = lambda*rad(it1) + (1.0d0-lambda)*rad(it0)
         eps(it) = lambda*eps(it1) + (1.0d0-lambda)*eps(it0)
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (srad(maxclass))
      allocate (seps(maxclass))
c
c     get the square roots of the vdw radii and well depths
c
      do i = 1, maxclass
         srad(i) = sqrt(rad(i))
         seps(i) = sqrt(eps(i))
      end do
c
c     use combination rules to set pairwise vdw radii sums
c
      do j = 1, nmut
         i = imut(j)
         it = class(i)
         do k = 1, maxclass
            if (rad(it).eq.0.0d0 .and. rad(k).eq.0.0d0) then
               rd = 0.0d0
            else if (radrule(1:10) .eq. 'ARITHMETIC') then
               rd = rad(it) + rad(k)
            else if (radrule(1:9) .eq. 'GEOMETRIC') then
               rd = 2.0d0*(srad(it)*srad(k))
            else if (radrule(1:10) .eq. 'CUBIC-MEAN') then
               rd = 2.0d0*(rad(it)**3+rad(k)**3)/(rad(it)**2+rad(k)**2)
            else
               rd = rad(it) + rad(k)
            end if
            radmin(it,k) = rd
            radmin(k,it) = rd
         end do
      end do
c
c     use combination rules to set pairwise well depths
c
      do j = 1, nmut
         i = imut(j)
         it = class(i)
         do k = 1, maxclass
            if (eps(it).eq.0.0d0 .and. eps(k).eq.0.0d0) then
               ep = 0.0d0
            else if (epsrule(1:10) .eq. 'ARITHMETIC') then
               ep = 0.5d0 * (eps(it) + eps(k))
            else if (epsrule(1:9) .eq. 'GEOMETRIC') then
               ep = seps(it) * seps(k)
            else if (epsrule(1:8) .eq. 'HARMONIC') then
               ep = 2.0d0 * (eps(it)*eps(k)) / (eps(it)+eps(k))
            else if (epsrule(1:3) .eq. 'HHG') then
               ep = 4.0d0 * (eps(it)*eps(k)) / (seps(it)+seps(k))**2
            else
               ep = seps(it) * seps(k)
            end if
            epsilon(it,k) = ep
            epsilon(k,it) = ep
         end do
      end do
c
c     print the van der Waals parameters for hybrid atoms
c
      header = .true.
      do j = 1, nmut
         if (verbose) then
            if (header) then
               header = .false.
               write (iout,10)
   10          format (/,' Hybrid van der Waals Parameters :',
     &                 //,7x,'Atom Number    Radius     Epsilon')
            end if
            radius = rad(it)
            if (radsiz .eq. 'DIAMETER')  radius = 2.0d0 * radius
            if (radtyp .eq. 'SIGMA')  radius = radius / twosix
            write (iout,20)  i,radius,eps(it)
   20       format (6x,i8,f14.4,f12.4)
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (srad)
      deallocate (seps)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine hcharge  --  find hybrid charge parameters  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "hcharge" constructs hybrid charge interaction parameters
c     given an initial state, final state and "lambda" value
c
c
      subroutine hcharge
      use sizes
      use atoms
      use charge
      use inform
      use iounit
      use kchrge
      use mutant
      implicit none
      integer i,j,k
      integer it,it0,it1
      real*8 chg0,chg1
      real*8 hybchg
      logical header,used
c
c
c     assign the hybrid parameters for atomic charges
c
      header = .true.
      do j = 1, nmut
         used = .false.
         i = imut(j)
         it = type(i)
         it0 = type0(j)
         it1 = type1(j)
         chg0 = chg(it0)
         chg1 = chg(it1)
         hybchg = lambda*chg1 + (1.0d0-lambda)*chg0
         do k = 1, nion
            if (iion(k) .eq. i) then
               used = .true.
               pchg(k) = hybchg
               goto 10
            end if
         end do
         if (chg0.ne.0.0d0 .or. chg1.ne.0.0d0) then
            used = .true.
            nion = nion + 1
            iion(nion) = i
            kion(nion) = i
            pchg(nion) = hybchg
         end if
   10    continue
         if (verbose .and. used) then
            if (header) then
               header = .false.
               write (iout,20)
   20          format (/,' Hybrid Atomic Partial Charge Parameters :',
     &                 //,7x,'Atom Number',7x,'Charge',/)
            end if
            write (iout,30)  i,hybchg
   30       format (6x,i8,5x,f12.3)
         end if
      end do
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine hdipole  --  find hybrid dipole parameters  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "hdipole" constructs hybrid dipole interaction parameters
c     given an initial state, final state and "lambda" value
c
c
      subroutine hdipole
      use sizes
      use atoms
      use bndstr
      use dipole
      use inform
      use iounit
      use kdipol
      use mutant
      implicit none
      integer i,j,k
      integer ia,ib
      integer ita,itb
      integer size
      real*8 dpl0,dpl1,hybdpl
      real*8 pos0,pos1,hybpos
      logical header,used
      character*4 pa,pb
      character*8 blank,pt
c
c
c     assign the hybrid parameters for bond dipoles
c
      blank = '      '
      header = .true.
      do i = 1, nbond
         ia = ibnd(1,i)
         ib = ibnd(2,i)
         if (mut(ia) .or. mut(ib)) then
            ita = type(ia)
            itb = type(ib)
c
c     find the dipole parameters for the initial state
c
            do j = 1, nmut
               k = imut(j)
               if (k .eq. ia)  ita = type0(j)
               if (k .eq. ib)  itb = type0(j)
            end do
            size = 4
            call numeral (ita,pa,size)
            call numeral (itb,pb,size)
            if (ita .le. itb) then
               pt = pa//pb
            else
               pt = pb//pa
            end if
            dpl0 = 0.0d0
            pos0 = 0.5d0
            do j = 1, maxnd
               if (kd(j) .eq. blank)  goto 10
               if (kd(j) .eq. pt) then
                  if (ita .le. itb) then
                     dpl0 = bdpl(j)
                     pos0 = sdpl(j)
                  else
                     dpl0 = -bdpl(j)
                     pos0 = 1.0d0 - sdpl(j)
                  end if
               end if
            end do
   10       continue
c
c     find the dipole parameters for the final state
c
            do j = 1, nmut
               k = imut(j)
               if (k .eq. ia)  ita = type1(j)
               if (k .eq. ib)  itb = type1(j)
            end do
            size = 4
            call numeral (ita,pa,size)
            call numeral (itb,pb,size)
            if (ita .le. itb) then
               pt = pa//pb
            else
               pt = pb//pa
            end if
            dpl1 = 0.0d0
            pos1 = 0.5d0
            do j = 1, maxnd
               if (kd(j) .eq. blank)  goto 20
               if (kd(j) .eq. pt) then
                  if (ita .le. itb) then
                     dpl1 = bdpl(j)
                     pos1 = sdpl(j)
                  else
                     dpl1 = -bdpl(j)
                     pos1 = 1.0d0 - sdpl(j)
                  end if
               end if
            end do
   20       continue
c
c     form the hybrid parameters for the current dipole
c
            hybdpl = lambda*dpl1 + (1.0d0-lambda)*dpl0
            hybpos = lambda*pos1 + (1.0d0-lambda)*pos0
            used = .false.
            do j = 1, ndipole
               if ((idpl(1,j).eq.ia .and. idpl(2,j).eq.ib) .or.
     &             (idpl(1,j).eq.ib .and. idpl(2,j).eq.ia)) then
                  idpl(1,j) = ia
                  idpl(2,j) = ib
                  bdpl(j) = hybdpl
                  sdpl(j) = hybpos
                  used = .true.
                  goto 30
               end if
            end do
            if (hybdpl .ne. 0.0d0) then
               ndipole = ndipole + 1
               idpl(1,ndipole) = ia
               idpl(2,ndipole) = ib
               bdpl(ndipole) = hybdpl
               sdpl(ndipole) = hybpos
               used = .true.
            end if
   30       continue
            if (verbose .and. used) then
               if (header) then
                  header = .false.
                  write (iout,40)
   40             format (/,' Hybrid Bond Dipole Moment Parameters :',
     &                    //,6x,'Atom Numbers',7x,'Moment',
     &                       7x,'Position',/)
               end if
               write (iout,50)  ia,ib,hybdpl,hybpos
   50          format (6x,2i5,2f15.3)
            end if
         end if
      end do
      return
      end
