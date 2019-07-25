c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1991  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine kimptor  --  improper torsion parameters  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "kimptor" assigns torsional parameters to each improper
c     torsion in the structure and processes any changed values
c
c
      subroutine kimptor
      use sizes
      use atomid
      use atoms
      use couple
      use imptor
      use inform
      use iounit
      use keys
      use kitors
      use math
      use potent
      use tors
      implicit none
      integer i,j,k,nti
      integer ia,ib,ic,id
      integer ita,itb,itc,itd
      integer size,next
      integer ft(6)
      real*8 angle,symm
      real*8 vt(6),st(6)
      logical header,done
      character*4 pa,pb,pc,pd
      character*4 zeros
      character*16 blank,pti
      character*16 pt0,pt1
      character*16 pt2,pt3
      character*16 pt(6)
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     process keywords containing improper torsion parameters
c
      blank = '                '
      zeros = '0000'
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:8) .eq. 'IMPTORS ') then
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
     &                                     (vt(j),st(j),ft(j),j=1,3)
   10       continue
            size = 4
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            pti = pa//pb//pc//pd
            call torphase (ft,vt,st)
            if (.not. silent) then
               if (header) then
                  header = .false.
                  write (iout,20)
   20             format (/,' Additional Improper Torsion Parameters :',
     &                    //,5x,'Atom Classes',15x,'1-Fold',12x,
     &                       '2-Fold',12x,'3-Fold',/)
               end if
               write (iout,30)  ia,ib,ic,id,(vt(j),st(j),j=1,3)
   30          format (4x,4i4,2x,3(f11.3,f7.1))
            end if
            do j = 1, maxnti
               if (kti(j).eq.blank .or. kti(j).eq.pti) then
                  kti(j) = pti
                  ti1(1,j) = vt(1)
                  ti1(2,j) = st(1)
                  ti2(1,j) = vt(2)
                  ti2(2,j) = st(2)
                  ti3(1,j) = vt(3)
                  ti3(2,j) = st(3)
                  goto 50
               end if
            end do
            write (iout,40)
   40       format (/,' KIMPTOR  --  Too many Improper Torsion',
     &                 ' Parameters')
            abort = .true.
   50       continue
         end if
      end do
c
c     determine the total number of forcefield parameters
c
      nti = maxnti
      do i = maxnti, 1, -1
         if (kti(i) .eq. blank)  nti = i - 1
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(iitors))  deallocate (iitors)
      if (allocated(itors1))  deallocate (itors1)
      if (allocated(itors2))  deallocate (itors2)
      if (allocated(itors3))  deallocate (itors3)
      allocate (iitors(4,6*n))
      allocate (itors1(4,6*n))
      allocate (itors2(4,6*n))
      allocate (itors3(4,6*n))
c
c     assign improper torsional parameters for each improper torsion;
c     multiple symmetrical parameters are given partial weights
c
      nitors = 0
      if (nti .ne. 0) then
         do i = 1, n
            if (n12(i) .eq. 3) then
               ia = i12(1,i)
               ib = i12(2,i)
               ic = i
               id = i12(3,i)
               ita = class(ia)
               itb = class(ib)
               itc = class(ic)
               itd = class(id)
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
               pt3 = zeros//zeros//pc//pd
               pt2 = zeros//zeros//pc//pb
               pt1 = zeros//zeros//pc//pa
               pt0 = zeros//zeros//pc//zeros
               symm = 1.0d0
               if (pa.eq.pb .or. pa.eq.pd .or. pb.eq.pd)  symm = 2.0d0
               if (pa.eq.pb .and. pa.eq.pd .and. pb.eq.pd)  symm = 6.0d0
               done = .false.
               do j = 1, nti
                  if (kti(j)(9:12) .eq. pc) then
                     do k = 1, 6
                        if (kti(j) .eq. pt(k)) then
                           nitors = nitors + 1
                           iitors(3,nitors) = ic
                           if (k .eq. 1) then
                              iitors(1,nitors) = ia
                              iitors(2,nitors) = ib
                              iitors(4,nitors) = id
                           else if (k .eq. 2) then
                              iitors(1,nitors) = ib
                              iitors(2,nitors) = ia
                              iitors(4,nitors) = id
                           else if (k .eq. 3) then
                              iitors(1,nitors) = ia
                              iitors(2,nitors) = id
                              iitors(4,nitors) = ib
                           else if (k .eq. 4) then
                              iitors(1,nitors) = id
                              iitors(2,nitors) = ia
                              iitors(4,nitors) = ib
                           else if (k .eq. 5) then
                              iitors(1,nitors) = ib
                              iitors(2,nitors) = id
                              iitors(4,nitors) = ia
                           else if (k .eq. 6) then
                              iitors(1,nitors) = id
                              iitors(2,nitors) = ib
                              iitors(4,nitors) = ia
                           end if
                           itors1(1,nitors) = ti1(1,j) / symm
                           itors1(2,nitors) = ti1(2,j)
                           itors2(1,nitors) = ti2(1,j) / symm
                           itors2(2,nitors) = ti2(2,j)
                           itors3(1,nitors) = ti3(1,j) / symm
                           itors3(2,nitors) = ti3(2,j)
                           done = .true.
                        end if
                     end do
                  end if
               end do
               if (.not. done) then
                  do j = 1, nti
                     if (kti(j) .eq. pt1) then
                        symm = 3.0d0
                        do k = 1, 3
                           nitors = nitors + 1
                           iitors(3,nitors) = ic
                           if (k .eq. 1) then
                              iitors(1,nitors) = ia
                              iitors(2,nitors) = ib
                              iitors(4,nitors) = id
                           else if (k .eq. 2) then
                              iitors(1,nitors) = ib
                              iitors(2,nitors) = id
                              iitors(4,nitors) = ia
                           else if (k .eq. 3) then
                              iitors(1,nitors) = id
                              iitors(2,nitors) = ia
                              iitors(4,nitors) = ib
                           end if
                           itors1(1,nitors) = ti1(1,j) / symm
                           itors1(2,nitors) = ti1(2,j)
                           itors2(1,nitors) = ti2(1,j) / symm
                           itors2(2,nitors) = ti2(2,j)
                           itors3(1,nitors) = ti3(1,j) / symm
                           itors3(2,nitors) = ti3(2,j)
                        end do
                        done = .true.
                     else if (kti(j) .eq. pt2) then
                        symm = 3.0d0
                        do k = 1, 3
                           nitors = nitors + 1
                           iitors(3,nitors) = ic
                           if (k .eq. 1) then
                              iitors(1,nitors) = ia
                              iitors(2,nitors) = ib
                              iitors(4,nitors) = id
                           else if (k .eq. 2) then
                              iitors(1,nitors) = ib
                              iitors(2,nitors) = id
                              iitors(4,nitors) = ia
                           else if (k .eq. 3) then
                              iitors(1,nitors) = id
                              iitors(2,nitors) = ia
                              iitors(4,nitors) = ib
                           end if
                           itors1(1,nitors) = ti1(1,j) / symm
                           itors1(2,nitors) = ti1(2,j)
                           itors2(1,nitors) = ti2(1,j) / symm
                           itors2(2,nitors) = ti2(2,j)
                           itors3(1,nitors) = ti3(1,j) / symm
                           itors3(2,nitors) = ti3(2,j)
                        end do
                        done = .true.
                     else if (kti(j) .eq. pt3) then
                        symm = 3.0d0
                        do k = 1, 3
                           nitors = nitors + 1
                           iitors(3,nitors) = ic
                           if (k .eq. 1) then
                              iitors(1,nitors) = ia
                              iitors(2,nitors) = ib
                              iitors(4,nitors) = id
                           else if (k .eq. 2) then
                              iitors(1,nitors) = ib
                              iitors(2,nitors) = id
                              iitors(4,nitors) = ia
                           else if (k .eq. 3) then
                              iitors(1,nitors) = id
                              iitors(2,nitors) = ia
                              iitors(4,nitors) = ib
                           end if
                           itors1(1,nitors) = ti1(1,j) / symm
                           itors1(2,nitors) = ti1(2,j)
                           itors2(1,nitors) = ti2(1,j) / symm
                           itors2(2,nitors) = ti2(2,j)
                           itors3(1,nitors) = ti3(1,j) / symm
                           itors3(2,nitors) = ti3(2,j)
                        end do
                        done = .true.
                     end if
                  end do
               end if
               if (.not. done) then
                  do j = 1, nti
                     if (kti(j) .eq. pt0) then
                        symm = 3.0d0
                        do k = 1, 3
                           nitors = nitors + 1
                           iitors(3,nitors) = ic
                           if (k .eq. 1) then
                              iitors(1,nitors) = ia
                              iitors(2,nitors) = ib
                              iitors(4,nitors) = id
                           else if (k .eq. 2) then
                              iitors(1,nitors) = ib
                              iitors(2,nitors) = id
                              iitors(4,nitors) = ia
                           else if (k .eq. 3) then
                              iitors(1,nitors) = id
                              iitors(2,nitors) = ia
                              iitors(4,nitors) = ib
                           end if
                           itors1(1,nitors) = ti1(1,j) / symm
                           itors1(2,nitors) = ti1(2,j)
                           itors2(1,nitors) = ti2(1,j) / symm
                           itors2(2,nitors) = ti2(2,j)
                           itors3(1,nitors) = ti3(1,j) / symm
                           itors3(2,nitors) = ti3(2,j)
                        end do
                     end if
                  end do
               end if
            end if
         end do
      end if
c
c     find the cosine and sine of the phase angle for each torsion
c
      do i = 1, nitors
         angle = itors1(2,i) / radian
         itors1(3,i) = cos(angle)
         itors1(4,i) = sin(angle)
         angle = itors2(2,i) / radian
         itors2(3,i) = cos(angle)
         itors2(4,i) = sin(angle)
         angle = itors3(2,i) / radian
         itors3(3,i) = cos(angle)
         itors3(4,i) = sin(angle)
      end do
c
c     turn off the improper torsional potential if it is not used
c
      if (nitors .eq. 0)  use_imptor = .false.
      return
      end
