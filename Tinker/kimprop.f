c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine kimprop  --  improper dihedral parameters  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "kimprop" assigns potential parameters to each improper
c     dihedral in the structure and processes any changed values
c
c
      subroutine kimprop
      use sizes
      use atomid
      use atoms
      use couple
      use improp
      use inform
      use iounit
      use keys
      use kiprop
      use potent
      use tors
      implicit none
      integer i,j,k,ndi
      integer ia,ib,ic,id
      integer ita,itb,itc,itd
      integer size,next
      real*8 tk,tv,symm
      logical header,done
      character*4 pa,pb,pc,pd
      character*8 zero8
      character*12 zero12
      character*16 blank,pti
      character*16 pt0,pt1
      character*16 pt2,pt3
      character*16 pt(6)
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     process keywords containing improper dihedral parameters
c
      blank = '                '
      zero8 = '00000000'
      zero12 = '000000000000'
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:9) .eq. 'IMPROPER ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            tk = 0.0d0
            tv = 0.0d0
            string = record(next:240)
            read (string,*,err=10,end=10)  ia,ib,ic,id,tk,tv
   10       continue
            size = 4
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            pti = pa//pb//pc//pd
            if (.not. silent) then
               if (header) then
                  header = .false.
                  write (iout,20)
   20             format (/,' Additional Improper Dihedral',
     &                       ' Parameters :',
     &                    //,5x,'Atom Classes',20x,'K(ID)',
     &                       7x,'Angle',/)
               end if
               write (iout,30)  ia,ib,ic,id,tk,tv
   30          format (4x,4i4,10x,2f12.3)
            end if
            do j = 1, maxndi
               if (kdi(j).eq.blank .or. kdi(j).eq.pti) then
                  kdi(j) = pti
                  dcon(j) = tk
                  tdi(j) = tv
                  goto 50
               end if
            end do
            write (iout,40)
   40       format (/,' KIMPROP  --  Too many Improper Dihedral',
     &                 ' Parameters')
            abort = .true.
   50       continue
         end if
      end do
c
c     determine the total number of forcefield parameters
c
      ndi = maxndi
      do i = maxndi, 1, -1
         if (kdi(i) .eq. blank)  ndi = i - 1
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(iiprop))  deallocate (iiprop)
      if (allocated(kprop))  deallocate (kprop)
      if (allocated(vprop))  deallocate (vprop)
      allocate (iiprop(4,6*n))
      allocate (kprop(6*n))
      allocate (vprop(6*n))
c
c     assign improper dihedral parameters for each improper angle;
c     multiple symmetrical parameters are given partial weights
c
      niprop = 0
      if (ndi .ne. 0) then
         do i = 1, n
            if (n12(i) .eq. 3) then
               ia = i
               ib = i12(1,i)
               ic = i12(2,i)
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
               pt(2) = pa//pb//pd//pc
               pt(3) = pa//pc//pb//pd
               pt(4) = pa//pc//pd//pb
               pt(5) = pa//pd//pb//pc
               pt(6) = pa//pd//pc//pb
               pt3 = pa//pb//zero8
               pt2 = pa//pc//zero8
               pt1 = pa//pd//zero8
               pt0 = pa//zero12
               symm = 1.0d0
               if (pb.eq.pc .or. pb.eq.pd .or. pc.eq.pd)  symm = 2.0d0
               if (pb.eq.pc .and. pb.eq.pd .and. pc.eq.pd)  symm = 6.0d0
               done = .false.
               do j = 1, ndi
                  if (kdi(j)(1:4) .eq. pa) then
                     do k = 1, 6
                        if (kdi(j) .eq. pt(k)) then
                           niprop = niprop + 1
                           iiprop(1,niprop) = ia
                           if (k .eq. 1) then
                              iiprop(2,niprop) = ib
                              iiprop(3,niprop) = ic
                              iiprop(4,niprop) = id
                           else if (k .eq. 2) then
                              iiprop(2,niprop) = ib
                              iiprop(3,niprop) = id
                              iiprop(4,niprop) = ic
                           else if (k .eq. 3) then
                              iiprop(2,niprop) = ic
                              iiprop(3,niprop) = ib
                              iiprop(4,niprop) = id
                           else if (k .eq. 4) then
                              iiprop(2,niprop) = ic
                              iiprop(3,niprop) = id
                              iiprop(4,niprop) = ib
                           else if (k .eq. 5) then
                              iiprop(2,niprop) = id
                              iiprop(3,niprop) = ib
                              iiprop(4,niprop) = ic
                           else if (k .eq. 6) then
                              iiprop(2,niprop) = id
                              iiprop(3,niprop) = ic
                              iiprop(4,niprop) = ib
                           end if
                           kprop(niprop) = dcon(j) / symm
                           vprop(niprop) = tdi(j)
                           done = .true.
                        end if
                     end do
                  end if
               end do
               if (.not. done) then
                  do j = 1, ndi
                     if (kdi(j) .eq. pt1) then
                        symm = 3.0d0
                        do k = 1, 3
                           niprop = niprop + 1
                           iiprop(1,niprop) = ia
                           if (k .eq. 1) then
                              iiprop(2,niprop) = ib
                              iiprop(3,niprop) = ic
                              iiprop(4,niprop) = id
                           else if (k .eq. 2) then
                              iiprop(2,niprop) = ic
                              iiprop(3,niprop) = id
                              iiprop(4,niprop) = ib
                           else if (k .eq. 3) then
                              iiprop(2,niprop) = id
                              iiprop(3,niprop) = ib
                              iiprop(4,niprop) = ic
                           end if
                           kprop(niprop) = dcon(j) / symm
                           vprop(niprop) = tdi(j)
                        end do
                        done = .true.
                     else if (kdi(j) .eq. pt2) then
                        symm = 3.0d0
                        do k = 1, 3
                           niprop = niprop + 1
                           iiprop(1,niprop) = ia
                           if (k .eq. 1) then
                              iiprop(2,niprop) = ib
                              iiprop(3,niprop) = ic
                              iiprop(4,niprop) = id
                           else if (k .eq. 2) then
                              iiprop(2,niprop) = ic
                              iiprop(3,niprop) = id
                              iiprop(4,niprop) = ib
                           else if (k .eq. 3) then
                              iiprop(2,niprop) = id
                              iiprop(3,niprop) = ib
                              iiprop(4,niprop) = ic
                           end if
                           kprop(niprop) = dcon(j) / symm
                           vprop(niprop) = tdi(j)
                        end do
                        done = .true.
                     else if (kdi(j) .eq. pt3) then
                        symm = 3.0d0
                        do k = 1, 3
                           niprop = niprop + 1
                           iiprop(1,niprop) = ia
                           if (k .eq. 1) then
                              iiprop(2,niprop) = ib
                              iiprop(3,niprop) = ic
                              iiprop(4,niprop) = id
                           else if (k .eq. 2) then
                              iiprop(2,niprop) = ic
                              iiprop(3,niprop) = id
                              iiprop(4,niprop) = ib
                           else if (k .eq. 3) then
                              iiprop(2,niprop) = id
                              iiprop(3,niprop) = ib
                              iiprop(4,niprop) = ic
                           end if
                           kprop(niprop) = dcon(j) / symm
                           vprop(niprop) = tdi(j)
                        end do
                        done = .true.
                     end if
                  end do
               end if
               if (.not. done) then
                  do j = 1, ndi
                     if (kdi(j) .eq. pt0) then
                        symm = 3.0d0
                        do k = 1, 3
                           niprop = niprop + 1
                           iiprop(1,niprop) = ia
                           if (k .eq. 1) then
                              iiprop(2,niprop) = ib
                              iiprop(3,niprop) = ic
                              iiprop(4,niprop) = id
                           else if (k .eq. 2) then
                              iiprop(2,niprop) = ic
                              iiprop(3,niprop) = id
                              iiprop(4,niprop) = ib
                           else if (k .eq. 3) then
                              iiprop(2,niprop) = id
                              iiprop(3,niprop) = ib
                              iiprop(4,niprop) = ic
                           end if
                           kprop(niprop) = dcon(j) / symm
                           vprop(niprop) = tdi(j)
                        end do
                     end if
                  end do
               end if
            end if
         end do
      end if
c
c     turn off the improper dihedral potential if it is not used
c
      if (niprop .eq. 0)  use_improp = .false.
      return
      end
