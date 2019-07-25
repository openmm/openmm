c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kurey  --  Urey-Bradley parameter assignment  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kurey" assigns the force constants and ideal distances
c     for the Urey-Bradley 1-3 interactions; also processes any
c     new or changed parameter values
c
c
      subroutine kurey
      use sizes
      use angbnd
      use atomid
      use atoms
      use inform
      use iounit
      use keys
      use kurybr
      use potent
      use urey
      implicit none
      integer i,j,nu
      integer ia,ib,ic
      integer ita,itb,itc
      integer size,next
      real*8 bb,tt
      logical header
      character*4 pa,pb,pc
      character*12 blank,pt
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     process keywords containing Urey-Bradley parameters
c
      blank = '            '
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:9) .eq. 'UREYBRAD ') then
            ia = 0
            ib = 0
            ic = 0
            bb = 0.0d0
            tt = 0.0d0
            string = record(next:240)
            read (string,*,err=10,end=10)  ia,ib,ic,bb,tt
   10       continue
            if (.not. silent) then
               if (header) then
                  header = .false.
                  write (iout,20)
   20             format (/,' Additional Urey-Bradley Parameters :',
     &                    //,5x,'Atom Classes',8x,'K(UB)',5x,
     &                       'Distance',/)
               end if
               write (iout,30)  ia,ib,ic,bb,tt
   30          format (4x,3i4,2x,f12.3,f12.4)
            end if
            size = 4
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            if (ia .le. ic) then
               pt = pa//pb//pc
            else
               pt = pc//pb//pa
            end if
            do j = 1, maxnu
               if (ku(j).eq.blank .or. ku(j).eq.pt) then
                  ku(j) = pt
                  ucon(j) = bb
                  dst13(j) = tt
                  goto 50
               end if
            end do
            write (iout,40)
   40       format (/,' KUREY  --  Too many Urey-Bradley',
     &                 ' Interaction Parameters')
            abort = .true.
   50       continue
         end if
      end do
c
c     determine the total number of forcefield parameters
c
      nu = maxnu
      do i = maxnu, 1, -1
         if (ku(i) .eq. blank)  nu = i - 1
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(iury))  deallocate (iury)
      if (allocated(uk))  deallocate (uk)
      if (allocated(ul))  deallocate (ul)
      allocate (iury(3,nangle))
      allocate (uk(nangle))
      allocate (ul(nangle))
c
c     assign the Urey-Bradley parameters for each angle
c
      nurey = 0
      if (nu .ne. 0) then
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
            do j = 1, nu
               if (ku(j) .eq. pt) then
                  nurey = nurey + 1
                  iury(1,nurey) = ia
                  iury(2,nurey) = ib
                  iury(3,nurey) = ic
                  uk(nurey) = ucon(j)
                  ul(nurey) = dst13(j)
                  goto 60
               end if
            end do
   60       continue
         end do
      end if
c
c     turn off the Urey-Bradley potential if it is not used
c
      if (nurey .eq. 0)  use_urey = .false.
      return
      end
