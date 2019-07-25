c
c
c     ##########################################################
c     ##  COPYRIGHT (C) 2014 by Chao Lu & Jay William Ponder  ##
c     ##                 All Rights Reserved                  ##
c     ##########################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine kangtor  --  angle-torsion parameter assignment  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "kangtor" assigns parameters for angle-torsion interactions
c     and processes new or changed parameter values
c
c
      subroutine kangtor
      use sizes
      use angtor
      use atmlst
      use atomid
      use atoms
      use couple
      use inform
      use iounit
      use keys
      use kantor
      use potent
      use tors
      implicit none
      integer i,j,k,l,m,nat
      integer ia,ib,ic,id
      integer ita,itb,itc,itd
      integer size,next
      real*8 at1,at2,at3
      real*8 at4,at5,at6
      logical header,swap
      character*4 pa,pb,pc,pd
      character*4 zeros
      character*16 blank
      character*16 pt
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     process keywords containing angle-torsion parameters
c
      blank = '                '
      zeros = '0000'
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:8) .eq. 'ANGTORS ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            at1 = 0.0d0
            at2 = 0.0d0
            at3 = 0.0d0
            at4 = 0.0d0
            at5 = 0.0d0
            at6 = 0.0d0
            string = record(next:240)
            read (string,*,err=10,end=10)  ia,ib,ic,id,at1,at2,
     &                                     at3,at4,at5,at6
   10       continue
            if (.not. silent) then
               if (header) then
                  header = .false.
                  write (iout,20)
   20             format (/,' Additional Angle-Torsion Parameters :',
     &                    //,31x,'1st Angle',18x,'2nd Angle',
     &                    /,5x,'Atom Classes',7x,'1-Fold',3x,'2-Fold',
     &                      3x,'3-Fold',3x,'1-Fold',3x,'2-Fold',
     &                      3x,'3-Fold'/)
               end if
               write (iout,30)  ia,ib,ic,id,at1,at2,at3,at4,at5,at6
   30          format (1x,4i4,4x,6f9.3)
            end if
            size = 4
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            if (ib .lt. ic) then
               pt = pa//pb//pc//pd
               swap = .false.
            else if (ic .lt. ib) then
               pt = pd//pc//pb//pa
               swap = .true.
            else if (ia .le. id) then
               pt = pa//pb//pc//pd
               swap = .false.
            else if (id .lt. ia) then
               pt = pd//pc//pb//pa
               swap = .true.
            end if
            do j = 1, maxnat
               if (kat(j).eq.blank .or. kat(j).eq.pt) then
                  kat(j) = pt
                  if (swap) then
                     atcon(1,j) = at4
                     atcon(2,j) = at5
                     atcon(3,j) = at6
                     atcon(4,j) = at1
                     atcon(5,j) = at2
                     atcon(6,j) = at3
                  else
                     atcon(1,j) = at1
                     atcon(2,j) = at2
                     atcon(3,j) = at3
                     atcon(4,j) = at4
                     atcon(5,j) = at5
                     atcon(6,j) = at6
                  end if
                  goto 50
               end if
            end do
            write (iout,40)
   40       format (/,' KANGTOR  --  Too many Angle-Torsion',
     &                 ' Parameters')
            abort = .true.
   50       continue
         end if
      end do
c
c     determine the total number of forcefield parameters
c
      nat = maxnat
      do i = maxnat, 1, -1
         if (kat(i) .eq. blank)  nat = i - 1
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(iat))  deallocate (iat)
      if (allocated(kant))  deallocate (kant)
      allocate (iat(3,ntors))
      allocate (kant(6,ntors))
c
c     assign the angle-torsion parameters for each torsion
c
      nangtor = 0
      if (nat .ne. 0) then
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
               swap = .false.
            else if (itc .lt. itb) then
               pt = pd//pc//pb//pa
               swap = .true.
            else if (ita .le. itd) then
               pt = pa//pb//pc//pd
               swap = .false.
            else if (itd .lt. ita) then
               pt = pd//pc//pb//pa
               swap = .true.
            end if
            do j = 1, nat
               if (kat(j) .eq. pt) then
                  nangtor = nangtor + 1
                  if (swap) then
                     kant(1,nangtor) = atcon(4,j)
                     kant(2,nangtor) = atcon(5,j)
                     kant(3,nangtor) = atcon(6,j)
                     kant(4,nangtor) = atcon(1,j)
                     kant(5,nangtor) = atcon(2,j)
                     kant(6,nangtor) = atcon(3,j)
                  else
                     kant(1,nangtor) = atcon(1,j)
                     kant(2,nangtor) = atcon(2,j)
                     kant(3,nangtor) = atcon(3,j)
                     kant(4,nangtor) = atcon(4,j)
                     kant(5,nangtor) = atcon(5,j)
                     kant(6,nangtor) = atcon(6,j)
                  end if
                  iat(1,nangtor) = i
                  m = 0
                  do k = 1, n12(ib)-1
                     do l = k+1, n12(ib)
                        m = m + 1
                        if ((i12(k,ib).eq.ia .and. i12(l,ib).eq.ic) .or.
     &                     (i12(k,ib).eq.ic .and. i12(l,ib).eq.ia)) then
                           iat(2,nangtor) = anglist(m,ib)
                           goto 60
                        end if
                     end do
                  end do
   60             continue
                  m = 0
                  do k = 1, n12(ic)-1
                     do l = k+1, n12(ic)
                        m = m + 1
                        if ((i12(k,ic).eq.ib .and. i12(l,ic).eq.id) .or.
     &                     (i12(k,ic).eq.id .and. i12(l,ic).eq.ib)) then
                           iat(3,nangtor) = anglist(m,ic)
                           goto 70
                        end if
                     end do
                  end do
               end if
            end do
   70       continue
         end do
      end if
c
c     turn off the angle-torsion potential if it is not used
c
      if (nangtor .eq. 0)  use_angtor = .false.
      return
      end
