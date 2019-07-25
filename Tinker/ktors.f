c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine ktors  --  torsional parameter assignment  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "ktors" assigns torsional parameters to each torsion in
c     the structure and processes any new or changed values
c
c
      subroutine ktors
      use sizes
      use atomid
      use atoms
      use couple
      use fields
      use inform
      use iounit
      use keys
      use ktorsn
      use math
      use potent
      use tors
      use usage
      implicit none
      integer i,j
      integer ia,ib,ic,id
      integer ita,itb,itc,itd
      integer nt,nt5,nt4
      integer size,next
      integer iring,minat
      integer nlist,ilist
      integer, allocatable :: kindex(:)
      integer ft(6)
      real*8 angle
      real*8 vt(6),st(6)
      logical header,done
      logical use_ring
      character*4 pa,pb,pc,pd
      character*4 zeros
      character*7 label
      character*16 blank
      character*16 pt,pt0
      character*16 pt1,pt2
      character*16, allocatable :: klist(:)
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     process keywords containing torsional angle parameters
c
      blank = '                '
      zeros = '0000'
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         iring = -1
         if (keyword(1:8) .eq. 'TORSION ')  iring = 0
         if (keyword(1:9) .eq. 'TORSION5 ')  iring = 5
         if (keyword(1:9) .eq. 'TORSION4 ')  iring = 4
         if (iring .ge. 0) then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            do j = 1, 6
               vt(j) = 0.0d0
               st(j) = 0.0d0
               ft(j) = 0
            end do
            string = record(next:240)
            read (string,*,err=10,end=10)  ia,ib,ic,id,
     &                                     (vt(j),st(j),ft(j),j=1,6)
   10       continue
            size = 4
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            if (ib .lt. ic) then
               pt = pa//pb//pc//pd
            else if (ic .lt. ib) then
               pt = pd//pc//pb//pa
            else if (ia .le. id) then
               pt = pa//pb//pc//pd
            else if (id .lt. ia) then
               pt = pd//pc//pb//pa
            end if
            call torphase (ft,vt,st)
            if (.not. silent) then
               if (header) then
                  header = .false.
                  write (iout,20)
   20             format (/,' Additional Torsional Parameters :',
     &                    //,5x,'Atom Classes',4x,'1-Fold',4x,'2-Fold',
     &                       4x,'3-Fold',4x,'4-Fold',4x,'5-Fold',
     &                       4x,'6-Fold',/)
               end if
               if (iring .eq. 0) then
                  write (iout,30)  ia,ib,ic,id,
     &                             (vt(j),nint(st(j)),j=1,6)
   30             format (1x,4i4,1x,6(f6.2,i4))
               else
                  if (iring .eq. 5)  label = '5-Ring '
                  if (iring .eq. 4)  label = '4-Ring '
                  write (iout,40)  ia,ib,ic,id,
     &                             (vt(j),nint(st(j)),j=1,6),label(1:6)
   40             format (1x,4i4,1x,6(f6.2,i4),3x,a6)
               end if
            end if
            if (iring .eq. 0) then
               do j = 1, maxnt
                  if (kt(j).eq.blank .or. kt(j).eq.pt) then
                     kt(j) = pt
                     t1(1,j) = vt(1)
                     t1(2,j) = st(1)
                     t2(1,j) = vt(2)
                     t2(2,j) = st(2)
                     t3(1,j) = vt(3)
                     t3(2,j) = st(3)
                     t4(1,j) = vt(4)
                     t4(2,j) = st(4)
                     t5(1,j) = vt(5)
                     t5(2,j) = st(5)
                     t6(1,j) = vt(6)
                     t6(2,j) = st(6)
                     goto 60
                  end if
               end do
               write (iout,50)
   50          format (/,' KTORS  --  Too many Torsional Angle',
     &                    ' Parameters')
               abort = .true.
   60          continue
            else if (iring .eq. 5) then
               do j = 1, maxnt5
                  if (kt5(j).eq.blank .or. kt5(j).eq.pt) then
                     kt5(j) = pt
                     t15(1,j) = vt(1)
                     t15(2,j) = st(1)
                     t25(1,j) = vt(2)
                     t25(2,j) = st(2)
                     t35(1,j) = vt(3)
                     t35(2,j) = st(3)
                     t45(1,j) = vt(4)
                     t45(2,j) = st(4)
                     t55(1,j) = vt(5)
                     t55(2,j) = st(5)
                     t65(1,j) = vt(6)
                     t65(2,j) = st(6)
                     goto 80
                  end if
               end do
               write (iout,70)
   70          format (/,' KTORS  --  Too many 5-Ring Torsional',
     &                    ' Parameters')
               abort = .true.
   80          continue
            else if (iring .eq. 4) then
               do j = 1, maxnt4
                  if (kt4(j).eq.blank .or. kt4(j).eq.pt) then
                     kt4(j) = pt
                     t14(1,j) = vt(1)
                     t14(2,j) = st(1)
                     t24(1,j) = vt(2)
                     t24(2,j) = st(2)
                     t34(1,j) = vt(3)
                     t34(2,j) = st(3)
                     t44(1,j) = vt(4)
                     t44(2,j) = st(4)
                     t54(1,j) = vt(5)
                     t54(2,j) = st(5)
                     t64(1,j) = vt(6)
                     t64(2,j) = st(6)
                     goto 100
                  end if
               end do
               write (iout,90)
   90          format (/,' KTORS  --  Too many 4-Ring Torsional',
     &                    ' Parameters')
               abort = .true.
  100          continue
            end if
         end if
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(tors1))  deallocate (tors1)
      if (allocated(tors2))  deallocate (tors2)
      if (allocated(tors3))  deallocate (tors3)
      if (allocated(tors4))  deallocate (tors4)
      if (allocated(tors5))  deallocate (tors5)
      if (allocated(tors6))  deallocate (tors6)
      allocate (tors1(4,ntors))
      allocate (tors2(4,ntors))
      allocate (tors3(4,ntors))
      allocate (tors4(4,ntors))
      allocate (tors5(4,ntors))
      allocate (tors6(4,ntors))
c
c     use special torsional parameter assignment method for MMFF
c
      if (forcefield .eq. 'MMFF94') then
         call ktorsm
         return
      end if
c
c     determine the total number of forcefield parameters
c
      nt = maxnt
      nt5 = maxnt5
      nt4 = maxnt4
      do i = maxnt, 1, -1
         if (kt(i) .eq. blank)  nt = i - 1
      end do
      do i = maxnt5, 1, -1
         if (kt5(i) .eq. blank)  nt5 = i - 1
      end do
      do i = maxnt4, 1, -1
         if (kt4(i) .eq. blank)  nt4 = i - 1
      end do
      use_ring = .false.
      if (min(nt5,nt4) .ne. 0)  use_ring = .true.
c
c     perform dynamic allocation of some local arrays
c
      allocate (kindex(maxnt))
      allocate (klist(maxnt))
c
c     assign torsional parameters for each torsional angle
c     by putting the parameter values into the "tors" arrays
c
      header = .true.
      nlist = 0
      do i = 1, ntors
         ia = itors(1,i)
         ib = itors(2,i)
         ic = itors(3,i)
         id = itors(4,i)
         ita = class(ia)
         itb = class(ib)
         itc = class(ic)
         itd = class(id)
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
         pt2 = zeros//pt(5:16)
         pt1 = pt(1:12)//zeros
         pt0 = zeros//pt(5:12)//zeros
         tors1(1,i) = 0.0d0
         tors1(2,i) = 0.0d0
         tors2(1,i) = 0.0d0
         tors2(2,i) = 0.0d0
         tors3(1,i) = 0.0d0
         tors3(2,i) = 0.0d0
         tors4(1,i) = 0.0d0
         tors4(2,i) = 0.0d0
         tors5(1,i) = 0.0d0
         tors5(2,i) = 0.0d0
         tors6(1,i) = 0.0d0
         tors6(2,i) = 0.0d0
         done = .false.
c
c     make a check for torsions inside small rings
c
         iring = 0
         if (use_ring) then
            call chkring (iring,ia,ib,ic,id)
            if (iring .eq. 6)  iring = 0
            if (iring.eq.5 .and. nt5.eq.0)  iring = 0
            if (iring.eq.4 .and. nt4.eq.0)  iring = 0
         end if
c
c     find parameters for this torsion; first check "klist"
c     to save time for angle types already located
c
         if (iring .eq. 0) then
            do j = 1, nlist
               if (klist(j) .eq. pt) then
                  ilist = kindex(j)
                  tors1(1,i) = tors1(1,ilist)
                  tors1(2,i) = tors1(2,ilist)
                  tors2(1,i) = tors2(1,ilist)
                  tors2(2,i) = tors2(2,ilist)
                  tors3(1,i) = tors3(1,ilist)
                  tors3(2,i) = tors3(2,ilist)
                  tors4(1,i) = tors4(1,ilist)
                  tors4(2,i) = tors4(2,ilist)
                  tors5(1,i) = tors5(1,ilist)
                  tors5(2,i) = tors5(2,ilist)
                  tors6(1,i) = tors6(1,ilist)
                  tors6(2,i) = tors6(2,ilist)
                  done = .true.
                  goto 110
               end if
            end do
            do j = 1, nt
               if (kt(j) .eq. pt) then
                  nlist = nlist + 1
                  klist(nlist) = pt
                  kindex(nlist) = i
                  tors1(1,i) = t1(1,j)
                  tors1(2,i) = t1(2,j)
                  tors2(1,i) = t2(1,j)
                  tors2(2,i) = t2(2,j)
                  tors3(1,i) = t3(1,j)
                  tors3(2,i) = t3(2,j)
                  tors4(1,i) = t4(1,j)
                  tors4(2,i) = t4(2,j)
                  tors5(1,i) = t5(1,j)
                  tors5(2,i) = t5(2,j)
                  tors6(1,i) = t6(1,j)
                  tors6(2,i) = t6(2,j)
                  done = .true.
                  goto 110
               end if
            end do
            do j = 1, nt
               if (kt(j).eq.pt1 .or. kt(j).eq.pt2) then
                  tors1(1,i) = t1(1,j)
                  tors1(2,i) = t1(2,j)
                  tors2(1,i) = t2(1,j)
                  tors2(2,i) = t2(2,j)
                  tors3(1,i) = t3(1,j)
                  tors3(2,i) = t3(2,j)
                  tors4(1,i) = t4(1,j)
                  tors4(2,i) = t4(2,j)
                  tors5(1,i) = t5(1,j)
                  tors5(2,i) = t5(2,j)
                  tors6(1,i) = t6(1,j)
                  tors6(2,i) = t6(2,j)
                  done = .true.
                  goto 110
               end if
            end do
            do j = 1, nt
               if (kt(j) .eq. pt0) then
                  tors1(1,i) = t1(1,j)
                  tors1(2,i) = t1(2,j)
                  tors2(1,i) = t2(1,j)
                  tors2(2,i) = t2(2,j)
                  tors3(1,i) = t3(1,j)
                  tors3(2,i) = t3(2,j)
                  tors4(1,i) = t4(1,j)
                  tors4(2,i) = t4(2,j)
                  tors5(1,i) = t5(1,j)
                  tors5(2,i) = t5(2,j)
                  tors6(1,i) = t6(1,j)
                  tors6(2,i) = t6(2,j)
                  done = .true.
                  goto 110
               end if
            end do
c
c     find the parameters for a 5-ring torsion
c
         else if (iring .eq. 5) then
            do j = 1, nt5
               if (kt5(j) .eq. pt) then
                  tors1(1,i) = t15(1,j)
                  tors1(2,i) = t15(2,j)
                  tors2(1,i) = t25(1,j)
                  tors2(2,i) = t25(2,j)
                  tors3(1,i) = t35(1,j)
                  tors3(2,i) = t35(2,j)
                  tors4(1,i) = t45(1,j)
                  tors4(2,i) = t45(2,j)
                  tors5(1,i) = t55(1,j)
                  tors5(2,i) = t55(2,j)
                  tors6(1,i) = t65(1,j)
                  tors6(2,i) = t65(2,j)
                  done = .true.
                  goto 110
               end if
            end do
            do j = 1, nt5
               if (kt5(j).eq.pt1 .or. kt5(j).eq.pt2) then
                  tors1(1,i) = t15(1,j)
                  tors1(2,i) = t15(2,j)
                  tors2(1,i) = t25(1,j)
                  tors2(2,i) = t25(2,j)
                  tors3(1,i) = t35(1,j)
                  tors3(2,i) = t35(2,j)
                  tors4(1,i) = t45(1,j)
                  tors4(2,i) = t45(2,j)
                  tors5(1,i) = t55(1,j)
                  tors5(2,i) = t55(2,j)
                  tors6(1,i) = t65(1,j)
                  tors6(2,i) = t65(2,j)
                  done = .true.
                  goto 110
               end if
            end do
            do j = 1, nt5
               if (kt5(j) .eq. pt0) then
                  tors1(1,i) = t15(1,j)
                  tors1(2,i) = t15(2,j)
                  tors2(1,i) = t25(1,j)
                  tors2(2,i) = t25(2,j)
                  tors3(1,i) = t35(1,j)
                  tors3(2,i) = t35(2,j)
                  tors4(1,i) = t45(1,j)
                  tors4(2,i) = t45(2,j)
                  tors5(1,i) = t55(1,j)
                  tors5(2,i) = t55(2,j)
                  tors6(1,i) = t65(1,j)
                  tors6(2,i) = t65(2,j)
                  done = .true.
                  goto 110
               end if
            end do
c
c     find the parameters for a 4-ring torsion
c
         else if (iring .eq. 4) then
            do j = 1, nt4
               if (kt4(j) .eq. pt) then
                  tors1(1,i) = t14(1,j)
                  tors1(2,i) = t14(2,j)
                  tors2(1,i) = t24(1,j)
                  tors2(2,i) = t24(2,j)
                  tors3(1,i) = t34(1,j)
                  tors3(2,i) = t34(2,j)
                  tors4(1,i) = t44(1,j)
                  tors4(2,i) = t44(2,j)
                  tors5(1,i) = t54(1,j)
                  tors5(2,i) = t54(2,j)
                  tors6(1,i) = t64(1,j)
                  tors6(2,i) = t64(2,j)
                  done = .true.
                  goto 110
               end if
            end do
            do j = 1, nt4
               if (kt4(j).eq.pt1 .or. kt4(j).eq.pt2) then
                  tors1(1,i) = t14(1,j)
                  tors1(2,i) = t14(2,j)
                  tors2(1,i) = t24(1,j)
                  tors2(2,i) = t24(2,j)
                  tors3(1,i) = t34(1,j)
                  tors3(2,i) = t34(2,j)
                  tors4(1,i) = t44(1,j)
                  tors4(2,i) = t44(2,j)
                  tors5(1,i) = t54(1,j)
                  tors5(2,i) = t54(2,j)
                  tors6(1,i) = t64(1,j)
                  tors6(2,i) = t64(2,j)
                  done = .true.
                  goto 110
               end if
            end do
            do j = 1, nt4
               if (kt4(j) .eq. pt0) then
                  tors1(1,i) = t14(1,j)
                  tors1(2,i) = t14(2,j)
                  tors2(1,i) = t24(1,j)
                  tors2(2,i) = t24(2,j)
                  tors3(1,i) = t34(1,j)
                  tors3(2,i) = t34(2,j)
                  tors4(1,i) = t44(1,j)
                  tors4(2,i) = t44(2,j)
                  tors5(1,i) = t54(1,j)
                  tors5(2,i) = t54(2,j)
                  tors6(1,i) = t64(1,j)
                  tors6(2,i) = t64(2,j)
                  done = .true.
                  goto 110
               end if
            end do
         end if
c
c     warning if suitable torsional parameter not found
c
  110    continue
         minat = min(atomic(ia),atomic(ib),atomic(ic),atomic(id))
         if (minat .eq. 0)  done = .true.
         if (use_tors .and. .not.done) then
            if (use(ia) .or. use(ib) .or. use(ic) .or. use(id))
     &         abort = .true.
            if (header) then
               header = .false.
               write (iout,120)
  120          format (/,' Undefined Torsional Parameters :',
     &                 //,' Type',24x,'Atom Names',24x,
     &                    'Atom Classes',/)
            end if
            label = 'Torsion'
            if (iring .eq. 5)  label = '5-Ring '
            if (iring .eq. 4)  label = '4-Ring '
            write (iout,130)  label,ia,name(ia),ib,name(ib),ic,
     &                        name(ic),id,name(id),ita,itb,itc,itd
  130       format (1x,a7,4x,4(i6,'-',a3),5x,4i5)
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (kindex)
      deallocate (klist)
c
c     find the cosine and sine of phase angle for each torsion
c
      do i = 1, ntors
         angle = tors1(2,i) / radian
         tors1(3,i) = cos(angle)
         tors1(4,i) = sin(angle)
         angle = tors2(2,i) / radian
         tors2(3,i) = cos(angle)
         tors2(4,i) = sin(angle)
         angle = tors3(2,i) / radian
         tors3(3,i) = cos(angle)
         tors3(4,i) = sin(angle)
         angle = tors4(2,i) / radian
         tors4(3,i) = cos(angle)
         tors4(4,i) = sin(angle)
         angle = tors5(2,i) / radian
         tors5(3,i) = cos(angle)
         tors5(4,i) = sin(angle)
         angle = tors6(2,i) / radian
         tors6(3,i) = cos(angle)
         tors6(4,i) = sin(angle)
      end do
c
c     turn off the torsional potential if it is not used
c
      if (ntors .eq. 0)  use_tors = .false.
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine ktorsm  --  assign MMFF torsional parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "ktorsm" assigns torsional parameters to each torsion according
c     to the Merck Molecular Force Field (MMFF)
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
      subroutine ktorsm
      use sizes
      use atomid
      use atoms
      use ktorsn
      use math
      use merck
      use potent
      use ring
      use tors
      implicit none
      integer i,j,k,l,m,o
      integer size,tt
      integer ia,ib,ic,id
      integer ita,itb,itc,itd
      integer inb,inc,irb,irc
      integer itta,ittb
      integer ittc,ittd
      integer nt4,nt5
      integer ab,bc,cd
      integer mclass
      real*8 angle
      real*8 beta,pi_bc,n_bc
      real*8 ub,vb,wb
      real*8 uc,vc,wc
      logical done,skipring
      logical ring4,ring5
      character*4 pa,pb,pc,pd
      character*16 pt,blank
c
c
c     determine the total number of forcefield parameters
c
      blank = '                '
      nt5 = maxnt5
      nt4 = maxnt4
      do i = maxnt5, 1, -1
         if (kt5(i) .eq. blank)  nt5 = i - 1
      end do
      do i = maxnt4, 1, -1
         if (kt4(i) .eq. blank)  nt4 = i - 1
      end do
c
c     assign MMFF torsional parameters for each torsional angle
c
      do i = 1, ntors
         ia = itors(1,i)
         ib = itors(2,i)
         ic = itors(3,i)
         id = itors(4,i)
         itta = type(ia)
         ittb = type(ib)
         ittc = type(ic)
         ittd = type(id)
         done = .false.
         mclass = 0
         skipring = .false.
   10    continue
c
c     determine the atom class equivalency assignments
c
         mclass = mclass + 1
         if (mclass .eq. 1) then
            ita = eqclass(itta,mclass)
            itb = eqclass(ittb,mclass)
            itc = eqclass(ittc,mclass)
            itd = eqclass(ittd,mclass)
         else if (mclass.eq.2) then
            ita = eqclass(itta,mclass)
            itb = eqclass(ittb,mclass)
            itc = eqclass(ittc,mclass)
            itd = eqclass(ittd,mclass)
         else if (mclass.eq.3) then
            ita = eqclass(itta,3)
            itb = eqclass(ittb,2)
            itc = eqclass(ittc,2)
            itd = eqclass(ittd,5)
         else if (mclass.eq.4) then
            ita = eqclass(itta,5)
            itb = eqclass(ittb,2)
            itc = eqclass(ittc,2)
            itd = eqclass(ittd,3)
         else if (mclass.eq.5) then
            ita = eqclass(itta,5)
            itb = eqclass(ittb,2)
            itc = eqclass(ittc,2)
            itd = eqclass(ittd,5)
         end if
c
c     construct search string and zero out parameters
c
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
         tors1(1,i) = 0.0d0
         tors1(2,i) = 0.0d0
         tors2(1,i) = 0.0d0
         tors2(2,i) = 0.0d0
         tors3(1,i) = 0.0d0
         tors3(2,i) = 0.0d0
         tors4(1,i) = 0.0d0
         tors4(2,i) = 0.0d0
         tors5(1,i) = 0.0d0
         tors5(2,i) = 0.0d0
         tors6(1,i) = 0.0d0
         tors6(2,i) = 0.0d0
         done = .false.
c
c     set the MMFF torsion type attribution
c
         ab = 0
         if (ia .le. ib) then
            do j = 1, nlignes
               if (ia.eq.bt_1(j,1) .and. ib.eq.bt_1(j,2)) then
                  ab = 1
               end if
            end do
         else if (ib .le. ia) then
            do j = 1, nlignes
               if (ib.eq.bt_1(j,1) .and. ia.eq.bt_1(j,2)) then
                  ab = 1
               end if
            end do
         end if
         bc = 0
         if (ib .le. ic) then
            do j = 1, nlignes
               if (ib.eq.bt_1(j,1) .and. ic.eq.bt_1(j,2)) then
                  bc = 1
               end if
            end do
         else if (ic .le. ib) then
            do j = 1, nlignes
               if (ic.eq.bt_1(j,1) .and. ib.eq.bt_1(j,2)) then
                  bc = 1
               end if
            end do
         end if
         cd = 0
         if (ic .le. id) then
            do j = 1, nlignes
               if (ic.eq.bt_1(j,1) .and. id.eq.bt_1(j,2)) then
                  cd = 1
               end if
            end do
         else if (id .le. ic) then
            do j = 1, nlignes
               if (id.eq.bt_1(j,1) .and. ic.eq.bt_1(j,2)) then
                  cd = 1
               end if
            end do
         end if
c
c     make a check for torsions inside small rings
c
         ring4 = .false.
         ring5 = .false.
         do j = 1, nring4
            do k = 1, 4
               if (ia .eq. iring4(k,j)) then
                  do l = 1, 4
                     if (ib .eq. iring4(l,j)) then
                        do m = 1, 4
                           if (ic .eq. iring4(m,j)) then
                              do o = 1, 4
                                 if (id .eq. iring4(o,j))
     &                              ring4 = .true.
                              end do
                           end if
                        end do
                     end if
                  end do
               end if
            end do
         end do
         do j = 1, nring5
            do k = 1, 5
               if (ia .eq. iring5(k,j)) then
                  do l = 1, 5
                     if (ib .eq. iring5(l,j)) then
                        do m = 1, 5
                           if (ic .eq. iring5(m,j)) then
                              do o = 1, 5
                                 if (id .eq. iring5(o,j))
     &                              ring5 = .true.
                              end do
                           end if
                        end do
                     end if
                  end do
               end if
            end do
         end do
         if (skipring) then
            ring4 = .false.
            ring5 = .false.
         end if
         if (ring4) then
            tt = 4
            do j = 1, nt4
               if (kt4(j) .eq. pt) then
                  tors1(1,i) = t14(1,j)
                  tors1(2,i) = t14(2,j)
                  tors2(1,i) = t24(1,j)
                  tors2(2,i) = t24(2,j)
                  tors3(1,i) = t34(1,j)
                  tors3(2,i) = t34(2,j)
                  done = .true.
                  goto 20
               end if
            end do
            if (.not.done .and. mclass.lt.5) then
               goto 10
            end if
         end if
         if (ring5 .and. (class(ia).eq.1.or.class(ib).eq.1.or.
     &                    class(ic).eq.1.or.class(id).eq.1)) then
            tt = 5
            do j = 1, nt5
               if (kt5(j) .eq. pt) then
                  tors1(1,i) = t15(1,j)
                  tors1(2,i) = t15(2,j)
                  tors2(1,i) = t25(1,j)
                  tors2(2,i) = t25(2,j)
                  tors3(1,i) = t35(1,j)
                  tors3(2,i) = t35(2,j)
                  done = .true.
               end if
            end do
            if (.not.done .and. mclass.lt.5) then
               goto 10
            else if (.not.done .and. mclass.eq.5) then
               mclass = 0
               skipring = .true.
               goto 10
            end if
         end if
c
c     condition below deduced from validation suite comparison
c
         if ((ab.eq.1 .and. (mltb(class(ic)).eq.0.or.
     &                       sbmb(class(ic)).eq.0)) .or.
     &       (cd.eq.1 .and. (mltb(class(ib)).eq.0.or.
     &                       sbmb(class(ib)).eq.0))) then
            tt = 2
            do j = 1, maxnt
               if (kt_2(j) .eq. pt) then
                  tors1(1,i) = t1_2(1,j)
                  tors1(2,i) = t1_2(2,j)
                  tors2(1,i) = t2_2(1,j)
                  tors2(2,i) = t2_2(2,j)
                  tors3(1,i) = t3_2(1,j)
                  tors3(2,i) = t3_2(2,j)
                  done = .true.
                  goto 20
               end if
            end do
            if (.not.done .and. mclass.lt.5) then
               goto 10
            end if
            if (.not.done .and. mclass.eq.5) then
               tt = 0
               do j = 1, maxnt
                  if (kt(j) .eq. pt) then
                     tors1(1,i) = t1(1,j)
                     tors1(2,i) = t1(2,j)
                     tors2(1,i) = t2(1,j)
                     tors2(2,i) = t2(2,j)
                     tors3(1,i) = t3(1,j)
                     tors3(2,i) = t3(2,j)
                     done = .true.
                     goto 20
                  end if
               end do
               if (.not.done .and. mclass.lt.5) then
                  goto 10
               end if
            end if
            if (tors1(1,i) .eq. 1000.0d0)  done = .false.
            if (tors1(2,i) .eq. 1000.0d0)  done = .false.
            if (tors2(1,i) .eq. 1000.0d0)  done = .false.
            if (tors2(2,i) .eq. 1000.0d0)  done = .false.
            if (tors3(1,i) .eq. 1000.0d0)  done = .false.
            if (tors3(2,i) .eq. 1000.0d0)  done = .false.
            goto 20
         else if (bc .eq. 1) then
            tt = 1
            do j = 1, maxnt
               if (kt_1(j) .eq. pt) then
                  tors1(1,i) = t1_1(1,j)
                  tors1(2,i) = t1_1(2,j)
                  tors2(1,i) = t2_1(1,j)
                  tors2(2,i) = t2_1(2,j)
                  tors3(1,i) = t3_1(1,j)
                  tors3(2,i) = t3_1(2,j)
                  done = .true.
                  goto 20
               end if
            end do
            if (.not.done .and. mclass.lt.5) then
               goto 10
            end if
            if (tors1(1,i) .eq. 1000.0d0)  done = .false.
            if (tors1(2,i) .eq. 1000.0d0)  done = .false.
            if (tors2(1,i) .eq. 1000.0d0)  done = .false.
            if (tors2(2,i) .eq. 1000.0d0)  done = .false.
            if (tors3(1,i) .eq. 1000.0d0)  done = .false.
            if (tors3(2,i) .eq. 1000.0d0)  done = .false.
            goto 20
         else if (.not. done) then
            tt = 0
            do j = 1, maxnt
               if (kt(j) .eq. pt) then
                  tors1(1,i) = t1(1,j)
                  tors1(2,i) = t1(2,j)
                  tors2(1,i) = t2(1,j)
                  tors2(2,i) = t2(2,j)
                  tors3(1,i) = t3(1,j)
                  tors3(2,i) = t3(2,j)
                  done = .true.
                  goto 20
               end if
            end do
            if (.not.done .and. mclass.lt.5) then
               goto 10
            end if
            if (tors1(1,i) .eq. 1000.0d0)  done = .false.
            if (tors1(2,i) .eq. 1000.0d0)  done = .false.
            if (tors2(1,i) .eq. 1000.0d0)  done = .false.
            if (tors2(2,i) .eq. 1000.0d0)  done = .false.
            if (tors3(1,i) .eq. 1000.0d0)  done = .false.
            if (tors3(2,i) .eq. 1000.0d0)  done = .false.
            goto 20
         end if
   20    continue
c
c     use the empirical rules for parameter not located
c
         if (.not. done) then
            ita = class(ia)
            itb = class(ib)
            itc = class(ic)
            itd = class(id)
            inb = atomic(ib)
            inc = atomic(ic)
            if (inb .eq. 6) then
               ub = 2.0d0
               vb = 2.12d0
            else if (inb .eq. 7) then
               ub = 2.0d0
               vb = 1.5d0
            else if (inb .eq. 8) then
               ub = 2.0d0
               vb = 0.2d0
            else if (inb .eq. 14) then
               ub = 1.25d0
               vb = 1.22d0
            else if (inb .eq. 15) then
               ub = 1.25d0
               vb = 2.4d0
            else if (inb .eq. 16) then
               ub = 1.25d0
               vb = 0.49d0
            end if
            if (inc .eq. 6) then
               uc = 2.0d0
               vc = 2.12d0
            else if (inc .eq. 7) then
               uc = 2.0d0
               vc = 1.5d0
            else if (inc .eq. 8) then
               uc = 2.0d0
               vc = 0.2d0
            else if (inc .eq. 14) then
               uc = 1.25d0
               vc = 1.22d0
            else if (inc .eq. 15) then
               uc = 1.25d0
               vc = 2.4d0
            else if (inc .eq. 16) then
               uc = 1.25d0
               vc = 0.49d0
            end if
            n_bc = (crd(itb)-1) * (crd(itc)-1)
            if (inb.eq.1)  irb = 0
            if (inb.ge.3 .and. inb.le.10)  irb = 1
            if (inb.ge.11 .and. inb.le.18)  irb = 2
            if (inb.ge.19 .and. inb.le.36)  irb = 3
            if (inb.ge.37 .and. inb.le.54)  irb = 4
            if (inc.eq.1)  irc = 0
            if (inc.ge.3 .and. inc.le.10)  irc = 1
            if (inc.ge.11 .and. inc.le.18)  irc = 2
            if (inc.ge.19 .and. inc.le.36)  irc = 3
            if (inc.ge.37 .and. inc.le.54)  irc = 4
            if (lin(itb).eq.1 .or. lin(itc).eq.1) then
               tors1(1,i) = 0.0d0
               tors1(2,i) = 0.0d0
               tors2(1,i) = 0.0d0
               tors2(2,i) = 180.0d0
               tors3(1,i) = 0.0d0
               tors3(2,i) = 0.0d0
               done = .true.
               goto 20
            else if (arom(itb).eq.1 .and. arom(itc).eq.1) then
               if (pilp(itb).eq.0 .and. pilp(itc).eq.0) then
                  pi_bc = 0.5d0
               else
                  pi_bc = 0.3d0
               end if
               if ((val(itb).eq.3.and.val(itc).eq.4) .or.
     &             (val(itb).eq.4.and.val(itc).eq.3) .or.
     &             (val(itb).eq.4.and.val(itc).eq.34) .or.
     &             (val(itb).eq.34.and.val(itc).eq.4) .or.
     &             (val(itb).eq.34.and.val(itc).eq.3) .or.
     &             (val(itb).eq.3.and.val(itc).eq.34) .or.
     &             (val(itb).eq.34.and.val(itc).eq.34)) then
                  beta = 3.0d0
               else
                  beta = 6.0d0
               end if
               tors1(1,i) = 0.0d0
               tors1(2,i) = 0.0d0
               tors2(1,i) = beta * pi_bc * sqrt(ub*uc)
               tors2(2,i) = 180.0d0
               tors3(1,i) = 0.0d0
               tors3(2,i) = 0.0d0
               done = .true.
               goto 20
            else if ((mltb(itb).eq.2 .and. mltb(itc).eq.2) .or.
     &               (mltb(itc).eq.2 .and. mltb(itb).eq.2)) then
               beta = 6.0d0
               pi_bc = 1.0d0
               tors1(1,i) = 0.0d0
               tors1(2,i) = 0.0d0
               tors2(1,i) = beta * pi_bc * sqrt(ub*uc)
               tors2(2,i) = 180.0d0
               tors3(1,i) = 0.0d0
               tors3(2,i) = 0.0d0
               done = .true.
               goto 20
            else if (mltb(itb).eq.2 .or. mltb(itc).eq.2) then
               beta = 6.0d0
               pi_bc = 0.4d0
               tors1(1,i) = 0.0d0
               tors1(2,i) = 0.0d0
               tors2(1,i) = beta * pi_bc * sqrt(ub*uc)
               tors2(2,i) = 180.0d0
               tors3(1,i) = 0.0d0
               tors3(2,i) = 0.0d0
               done = .true.
               goto 20
            else if (crd(itb).eq.4 .and. crd(itc).eq.4) then
               tors1(1,i) = 0.0d0
               tors1(2,i) = 0.0d0
               tors2(1,i) = 0.0d0
               tors2(2,i) = 180.0d0
               tors3(1,i) = sqrt(vb*vc) / n_bc
               tors3(2,i) = 0.0d0
               done = .true.
               goto 20
            else if ((crd(itb).eq.4.and.crd(itc).eq.3.and.
     &              ((val(itc).eq.4.or.val(itc).eq.34).or.
     &                 mltb(itc).ne.0)) .or.
     &              (crd(itc).eq.4.and.crd(itb).eq.3.and.
     &              ((val(itb).eq.4.or.val(itb).eq.34).or.
     &                 mltb(itb).ne.0))) then
               tors1(1,i) = 0.0d0
               tors1(2,i) = 0.0d0
               tors2(1,i) = 0.0d0
               tors2(2,i) = 180.0d0
               tors3(1,i) = 0.0d0
               tors3(2,i) = 0.0d0
               done = .true.
               goto 20
            else if ((crd(itb).eq.4.and.crd(itc).eq.2.and.
     &               (val(itc).eq.3.or.mltb(itc).ne.0)) .or.
     &               (crd(itb).eq.4.and.crd(itc).eq.2.and.
     &               (val(itc).eq.3.or.mltb(itc).ne.0))) then
               tors1(1,i) = 0.0d0
               tors1(2,i) = 0.0d0
               tors2(1,i) = 0.0d0
               tors2(2,i) = 180.0d0
               tors3(1,i) = 0.0d0
               tors3(2,i) = 0.0d0
               done = .true.
               goto 20
            else if (crd(itb).eq.4 .or. crd(itc).eq.4) then
               tors1(1,i) = 0.0d0
               tors1(2,i) = 0.0d0
               tors2(1,i) = 0.0d0
               tors2(2,i) = 180.0d0
               tors3(1,i) = sqrt(vb*vc) / n_bc
               tors3(2,i) = 0.0d0
               done = .true.
               goto 20
            else if (pilp(itb).eq.1 .and. pilp(itc).eq.1) then
               tors1(1,i) = 0.0d0
               tors1(2,i) = 0.0d0
               tors2(1,i) = 0.0d0
               tors2(2,i) = 180.0d0
               tors3(1,i) = 0.0d0
               tors3(2,i) = 0.0d0
               done = .true.
               goto 20
            else if (pilp(itb).ne.0 .and. mltb(itc).ne.0) then
               beta = 6.0d0
               if (mltb(itb) .eq. 1) then
                  pi_bc = 0.5d0
               else if (irb.eq.1 .and. irc.eq.1) then
                  pi_bc = 0.3d0
               else if (irb.ne.1 .or. irc.ne.1) then
                  pi_bc = 0.15d0
               end if
               tors1(1,i) = 0.0d0
               tors1(2,i) = 0.0d0
               tors2(1,i) = beta * pi_bc * sqrt(ub*uc)
               tors2(2,i) = 180.0d0
               tors3(1,i) = 0.0d0
               tors3(2,i) = 0.0d0
               done = .true.
               goto 20
            else if (pilp(itc).ne.0 .and. mltb(itb).ne.0) then
               beta = 6.0d0
               if (mltb(itc) .eq. 1) then
                  pi_bc = 0.5d0
               else if (irb.eq.1 .and. irc.eq.1) then
                  pi_bc = 0.3d0
               else if (irb.ne.1 .or. irc.ne.1) then
                  pi_bc = 0.15d0
               end if
               tors1(1,i) = 0.0d0
               tors1(2,i) = 0.0d0
               tors2(1,i) = beta * pi_bc * sqrt(ub*uc)
               tors2(2,i) = 180.0d0
               tors3(1,i) = 0.0d0
               tors3(2,i) = 0.0d0
               done = .true.
               goto 20
            else if ((mltb(itb).eq.1.or.mltb(itc).eq.1) .and.
     &               (inb.ne.6.or.inc.ne.6)) then
               beta = 6.0d0
               pi_bc = 0.4d0
               tors1(1,i) = 0.0d0
               tors1(2,i) = 0.0d0
               tors2(1,i) = beta * pi_bc * sqrt(ub*uc)
               tors2(2,i) = 180.0d0
               tors3(1,i) = 0.0d0
               tors3(2,i) = 0.0d0
               done = .true.
               goto 20
            else if (mltb(itb).ne.0 .and. mltb(itc).ne.0) then
               beta = 6.0d0
               pi_bc = 0.15d0
               tors1(1,i) = 0.0d0
               tors1(2,i) = 0.0d0
               tors2(1,i) = beta * pi_bc * sqrt(ub*uc)
               tors2(2,i) = 180.0d0
               tors3(1,i) = 0.0d0
               tors3(2,i) = 0.0d0
               done = .true.
               goto 20
            else if (inb.eq.8 .and. inc.eq.8) then
               wb = 2.0d0
               wc = 2.0d0
               tors1(1,i) = 0.0d0
               tors1(2,i) = 0.0d0
               tors2(1,i) = -sqrt(wb*wc)
               tors2(2,i) = 180.0d0
               tors3(1,i) = 0.0d0
               tors3(2,i) = 0.0d0
               done = .true.
               goto 20
            else if ((inb.eq.8.and.inc.eq.16) .or.
     &               (inb.eq.16.and.inc.eq.8)) then
               wb = 2.0d0
               wc = 8.0d0
               tors1(1,i) = 0.0d0
               tors1(2,i) = 0.0d0
               tors2(1,i) = -sqrt(wb*wc)
               tors2(2,i) = 180.0d0
               tors3(1,i) = 0.0d0
               tors3(2,i) = 0.0d0
               done = .true.
               goto 20
            else if (inb.eq.16 .and. inc.eq.16) then
               wb = 8.0d0
               wc = 8.0d0
               tors1(1,i) = 0.0d0
               tors1(2,i) = 0.0d0
               tors2(1,i) = -sqrt(wb*wc)
               tors2(2,i) = 180.0d0
               tors3(1,i) = 0.0d0
               tors3(2,i) = 0.0d0
               done = .true.
               goto 20
            else
               tors1(1,i) = 0.0
               tors1(2,i) = 0.0
               tors2(1,i) = 0.0
               tors2(2,i) = 180.0
               tors3(1,i) = sqrt(vb*vc) / n_bc
               tors3(2,i) = 0.0
               done = .true.
               goto 20
            end if
         end if
      end do
c
c     find the cosine and sine of phase angle for each torsion
c
      do i = 1, ntors
         angle = tors1(2,i) / radian
         tors1(3,i) = cos(angle)
         tors1(4,i) = sin(angle)
         angle = tors2(2,i) / radian
         tors2(3,i) = cos(angle)
         tors2(4,i) = sin(angle)
         angle = tors3(2,i) / radian
         tors3(3,i) = cos(angle)
         tors3(4,i) = sin(angle)
         angle = tors4(2,i) / radian
         tors4(3,i) = cos(angle)
         tors4(4,i) = sin(angle)
         angle = tors5(2,i) / radian
         tors5(3,i) = cos(angle)
         tors5(4,i) = sin(angle)
         angle = tors6(2,i) / radian
         tors6(3,i) = cos(angle)
         tors6(4,i) = sin(angle)
      end do
c
c     turn off the torsional potential if it is not used
c
      if (ntors .eq. 0)  use_tors = .false.
      return
      end
