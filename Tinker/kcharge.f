c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine kcharge  --  assign partial charge parameters  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "kcharge" assigns partial charges to the atoms within
c     the structure and processes any new or changed values
c
c
      subroutine kcharge
      use sizes
      use atomid
      use atoms
      use charge
      use chgpot
      use couple
      use fields
      use inform
      use iounit
      use kchrge
      use keys
      use potent
      implicit none
      integer i,j,k,m
      integer ia,next
      integer, allocatable :: list(:)
      integer, allocatable :: nc12(:)
      real*8 cg
      logical header
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     process keywords containing partial charge parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'CHARGE ') then
            ia = 0
            cg = 0.0d0
            string = record(next:240)
            read (string,*,err=40,end=40)  ia,cg
            if (ia .gt. 0) then
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,10)
   10             format (/,' Additional Atomic Partial Charge',
     &                       ' Parameters :',
     &                    //,5x,'Atom Type',10x,'Charge',/)
               end if
               if (ia .le. maxtyp) then
                  chg(ia) = cg
                  if (.not. silent) then
                     write (iout,20)  ia,cg
   20                format (4x,i6,8x,f12.4)
                  end if
               else
                  write (iout,30)
   30             format (/,' KCHARGE  --  Too many Partial Charge',
     &                       ' Parameters')
                  abort = .true.
               end if
            end if
   40       continue
         end if
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(iion))  deallocate (iion)
      if (allocated(jion))  deallocate (jion)
      if (allocated(kion))  deallocate (kion)
      if (allocated(pchg))  deallocate (pchg)
      allocate (iion(n))
      allocate (jion(n))
      allocate (kion(n))
      allocate (pchg(n))
c
c     find and store all the atomic partial charges
c
      do i = 1, n
         pchg(i) = chg(type(i))
      end do
c
c     use special charge parameter assignment method for MMFF
c
      if (forcefield .eq. 'MMFF94')  call kchargem
c
c     process keywords containing atom specific partial charges
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'CHARGE ') then
            ia = 0
            cg = 0.0d0
            string = record(next:240)
            read (string,*,err=70,end=70)  ia,cg
            if (ia.lt.0 .and. ia.ge.-n) then
               ia = -ia
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,50)
   50             format (/,' Additional Partial Charges for',
     &                       ' Specific Atoms :',
     &                    //,6x,'Atom',14x,'Charge',/)
               end if
               if (.not. silent) then
                  write (iout,60)  ia,cg
   60             format (4x,i6,8x,f12.4)
               end if
               pchg(ia) = cg
            end if
   70       continue
         end if
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (list(n))
      allocate (nc12(n))
c
c     remove zero partial charges from the list of charges
c
      nion = 0
      do i = 1, n
         list(i) = 0
         if (pchg(i) .ne. 0.0d0) then
            nion = nion + 1
            iion(nion) = i
            jion(nion) = i
            kion(nion) = i
            pchg(nion) = pchg(i)
            list(i) = nion
         end if
      end do
c
c     optionally use neutral groups for neighbors and cutoffs
c
      if (neutnbr .or. neutcut) then
         do i = 1, n
            nc12(i) = 0
            do j = 1, n12(i)
               k = list(i12(j,i))
               if (k .ne. 0)  nc12(i) = nc12(i) + 1
            end do
         end do
         do i = 1, nion
            k = iion(i)
            if (n12(k) .eq. 1) then
               do j = 1, n12(k)
                  m = i12(j,k)
                  if (nc12(m) .gt. 1) then
                     if (neutnbr)  jion(i) = m
                     if (neutcut)  kion(i) = m
                  end if
               end do
            end if
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (list)
      deallocate (nc12)
c
c     turn off charge-charge and charge-dipole terms if not used
c
      if (nion .eq. 0) then
         use_charge = .false.
         use_chgdpl = .false.
      end if
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine kchargem  --  assign MMFF charge parameters  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "kchargem" assigns partial charges to the atoms according to
c     the Merck Molecular Force Field (MMFF)
c
c
      subroutine kchargem
      use sizes
      use atomid
      use atoms
      use charge
      use couple
      use merck
      implicit none
      integer i,j,k,m
      integer it,kt,bt
      integer ic,kc
      real*8, allocatable :: pchg0(:)
      logical emprule
c
c
c     set and store MMFF base atomic partial charge values
c
      do i = 1, n
         it = type(i)
         pchg(i) = 0.0d0
         if (it .eq. 107)  pchg(i) = -0.5d0
         if (it .eq. 113) then
            pchg(i) = 0.0d0
            do j = 1, n12(i)
               k = i12(j,i)
               kt = type(k)
               if (kt .eq. 185)  pchg(i) = -0.5d0
            end do
         end if
         if (it .eq. 114)  pchg(i) = -1.0d0 / 3.0d0
         if (it .eq. 115)  pchg(i) = -3.0d0
         if (it .eq. 116)  pchg(i) = -0.5d0
         if (it .eq. 118)  pchg(i) = -0.5d0
         if (it .eq. 119)  pchg(i) = -2.0d0 / 3.0d0
         if (it .eq. 121)  pchg(i) = -0.25d0
         if (it .eq. 123)  pchg(i) = 1.0d0
         if (it .eq. 124)  pchg(i) = -1.0d0
         if (it .eq. 125)  pchg(i) = -1.0d0
         if (it .eq. 154)  pchg(i) = 1.0d0
         if (it .eq. 156)  pchg(i) = 1.0d0
         if (it .eq. 159)  pchg(i) = 1.0d0
         if (it .eq. 160)  pchg(i) = 1.0d0
         if (it .eq. 161)  pchg(i) = 0.5d0
         if (it .eq. 162)  pchg(i) = 1.0d0 / 3.0d0
         if (it .eq. 165)  pchg(i) = 1.0d0
         if (it .eq. 168) then
            do j = 1, n12(i)
               k = i12(j,i)
               kt = type(k)
               if (kt.eq.168 .or. kt.eq.142)  pchg(i) = 1.0d0
            end do
         end if
         if (it .eq. 169)  pchg(i) = -1.0d0
         if (it .eq. 182)  pchg(i) = -0.5d0
         if (it .eq. 183) then
            pchg(i) = -1.0d0
            do j = 1, n12(i)
               k = i12(j,i)
               kt = type(k)
               if (kt .eq. 87)  pchg(i) = -0.5d0
            end do
         end if
         if (it .eq. 195)  pchg(i) = 1.0d0
         if (it .eq. 196)  pchg(i) = 1.0d0
         if (it .eq. 197)  pchg(i) = 1.0d0
         if (it .eq. 201)  pchg(i) = 2.0d0
         if (it .eq. 202)  pchg(i) = 3.0d0
         if (it .eq. 203)  pchg(i) = -1.0d0
         if (it .eq. 204)  pchg(i) = -1.0d0
         if (it .eq. 205)  pchg(i) = -1.0d0
         if (it .eq. 206)  pchg(i) = 1.0d0
         if (it .eq. 207)  pchg(i) = 1.0d0
         if (it .eq. 208)  pchg(i) = 1.0d0
         if (it .eq. 209)  pchg(i) = 2.0d0
         if (it .eq. 210)  pchg(i) = 2.0d0
         if (it .eq. 211)  pchg(i) = 2.0d0
         if (it .eq. 212)  pchg(i) = 1.0d0
         if (it .eq. 213)  pchg(i) = 2.0d0
         if (it .eq. 214)  pchg(i) = 2.0d0
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (pchg0(n))
c
c     modify MMFF base charges using a bond increment scheme
c
      do i = 1, n
         pchg0(i) = pchg(i)
      end do
      do i = 1, n
         it = type(i)
         ic = class(i)
         if (pchg0(i).lt.0.0d0 .or. it.eq.162) then
            pchg(i) = (1.0d0-crd(ic)*fcadj(ic)) * pchg0(i)
         end if
         do j = 1, n12(i)
            k = i12(j,i)
            kt = type(k)
            kc = class(k)
            if (pchg0(k).lt.0.0d0 .or. kt.eq.162) then
               pchg(i) = pchg(i) + fcadj(kc)*pchg0(k)
            end if
            bt = 0
            do m = 1, nlignes
               if ((i.eq.bt_1(m,1) .and. i12(j,i).eq.bt_1(m,2)).or.
     &                (i12(j,i).eq.bt_1(m,1) .and. i.eq.bt_1(m,2))) then
                  bt = 1
               end if
            end do
            emprule = .false.
            if (bt .eq. 1) then
               pchg(i) = pchg(i) + bci_1(kc,ic)
               if (bci_1(kc,ic) .eq. 1000.0d0) then
                  emprule = .true.
                  goto 10
               end if
            else if (bt .eq. 0) then
               pchg(i) = pchg(i) + bci(kc,ic)
               if (bci(kc,ic) .eq. 1000.0d0) then
                  emprule = .true.
                  goto 10
               end if
            end if
         end do
   10    continue
         if (emprule) then
            pchg(i) = (1.0d0-crd(ic)*fcadj(ic)) * pchg0(i)
            do j = 1, n12(i)
               k = i12(j,i)
               kc = class(k)
               pchg(i) = pchg(i) + fcadj(kc)*pchg0(i12(j,i))
            end do
            do j = 1, n12(i)
               k = i12(j,i)
               kc = class(k)
               bt = 0
               do k = 1, nlignes
                  if ((i.eq.bt_1(k,1) .and.
     &                      i12(j,i).eq.bt_1(k,2)) .or.
     &                   (i12(j,i).eq.bt_1(k,1) .and.
     &                      i.eq.bt_1(k,2))) then
                     bt = 1
                  end if
               end do
               if (bt .eq. 1) then
                  if (bci_1(kc,ic) .eq. 1000.0d0) then
                     pchg(i) = pchg(i) + pbci(ic) - pbci(kc)
                  else
                     pchg(i) = pchg(i) + bci_1(kc,ic)
                  end if
               else if (bt .eq. 0) then
                  if (bci(kc,ic) .eq. 1000.0d0) then
                     pchg(i) = pchg(i) + pbci(ic) - pbci(kc)
                  else
                     pchg(i) = pchg(i) + bci(kc,ic)
                  end if
               end if
            end do
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (pchg0)
      return
      end
