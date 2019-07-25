c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine korbit  --  pisystem parameter assignment  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "korbit" assigns pi-orbital parameters to conjugated
c     systems and processes any new or changed parameters
c
c
      subroutine korbit
      use sizes
      use atomid
      use atoms
      use bndstr
      use inform
      use iounit
      use keys
      use korbs
      use orbits
      use piorbs
      use pistuf
      use tors
      use units
      implicit none
      integer i,j,k,jt
      integer ia,ib,ita,itb
      integer npi,npi5,npi4
      integer size,next,iring
      real*8 elect,ioniz
      real*8 repuls
      real*8 sslop,tslop
      logical header
      logical use_ring
      character*4 pa,pb
      character*6 label
      character*8 blank,pt
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     process keywords containing pisystem atom parameters
c
      blank = '        '
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'PIATOM ') then
            ia = 0
            elect = 0.0d0
            ioniz = 0.0d0
            repuls = 0.0d0
            string = record(next:240)
            read (string,*,err=10)  ia,elect,ioniz,repuls
   10       continue
            if (.not. silent) then
               if (header) then
                  header = .false.
                  write (iout,20)
   20             format (/,' Additional Pisystem Atom Parameters :',
     &                    //,6x,'Atom Type',3x,'Electron',
     &                       3x,'Ionization',3x,'Repulsion',/)
               end if
               write (iout,30)  ia,elect,ioniz,repuls
   30          format (8x,i4,3x,f10.3,3x,f10.3,2x,f10.3)
            end if
            if (ia.gt.0 .and. ia.le.maxclass) then
               electron(ia) = elect
               ionize(ia) = ioniz
               repulse(ia) = repuls
            else
               write (iout,40)
   40          format (/,' KORBIT  --  Too many Atom Classes;',
     &                    ' Increase MAXCLASS')
               abort = .true.
            end if
         end if
      end do
c
c     process keywords containing pisystem bond parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         iring = -1
         if (keyword(1:7) .eq. 'PIBOND ')  iring = 0
         if (keyword(1:8) .eq. 'PIBOND5 ')  iring = 5
         if (keyword(1:8) .eq. 'PIBOND4 ')  iring = 4
         if (iring .ge. 0) then
            ia = 0
            ib = 0
            sslop = 0.0d0
            tslop = 0.0d0
            string = record(next:240)
            read (string,*,err=50)  ia,ib,sslop,tslop
   50       continue
            if (.not. silent) then
               if (header) then
                  header = .false.
                  write (iout,60)
   60             format (/,' Additional Pisystem Bond Parameters :',
     &                    //,6x,'Atom Types',7x,'d Force',4x,
     &                       'd Length',/)
               end if
               if (iring .eq. 0) then
                  write (iout,70)  ia,ib,sslop,tslop
   70             format (6x,2i4,5x,2f11.3)
               else
                  if (iring .eq. 5)  label = '5-Ring'
                  if (iring .eq. 4)  label = '4-Ring'
                  write (iout,80)  ia,ib,sslop,tslop,label
   80             format (6x,2i4,5x,2f11.3,3x,a6)
               end if
            end if
            size = 4
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            if (ia .le. ib) then
               pt = pa//pb
            else
               pt = pb//pa
            end if
            if (iring .eq. 0) then
               do j = 1, maxnpi
                  if (kpi(j).eq.blank .or. kpi(j).eq.pt) then
                     kpi(j) = pt
                     sslope(j) = sslop
                     tslope(j) = tslop
                     goto 100
                  end if
               end do
               write (iout,90)
   90          format (/,' KORBIT  --  Too many Pisystem Bond',
     &                    ' Type Parameters')
               abort = .true.
  100          continue
            else if (iring .eq. 5) then
               do j = 1, maxnpi5
                  if (kpi5(j).eq.blank .or. kpi5(j).eq.pt) then
                     kpi5(j) = pt
                     sslope5(j) = sslop
                     tslope5(j) = tslop
                     goto 120
                  end if
               end do
               write (iout,110)
  110          format (/,' KORBIT  --  Too many 5-Ring Pisystem Bond',
     &                    ' Type Parameters')
               abort = .true.
  120          continue
            else if (iring .eq. 4) then
               do j = 1, maxnpi4
                  if (kpi4(j).eq.blank .or. kpi4(j).eq.pt) then
                     kpi4(j) = pt
                     sslope4(j) = sslop
                     tslope4(j) = tslop
                     goto 140
                  end if
               end do
               write (iout,130)
  130          format (/,' KORBIT  --  Too many 4-Ring Pisystem Bond',
     &                    ' Type Parameters')
               abort = .true.
  140          continue
            end if
         end if
      end do
c
c     determine the total number of forcefield parameters
c
      npi = maxnpi
      npi5 = maxnpi5
      npi4 = maxnpi4
      do i = maxnpi, 1, -1
         if (kpi(i) .eq. blank)  npi = i - 1
      end do
      do i = maxnpi5, 1, -1
         if (kpi5(i) .eq. blank)  npi5 = i - 1
      end do
      do i = maxnpi4, 1, -1
         if (kpi4(i) .eq. blank)  npi4 = i - 1
      end do
      use_ring = .false.
      if (min(npi5,npi4) .ne. 0)  use_ring = .true.
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(qorb))  deallocate (qorb)
      if (allocated(worb))  deallocate (worb)
      if (allocated(emorb))  deallocate (emorb)
      if (allocated(bkpi))  deallocate (bkpi)
      if (allocated(blpi))  deallocate (blpi)
      if (allocated(kslope))  deallocate (kslope)
      if (allocated(lslope))  deallocate (lslope)
      if (allocated(torsp2))  deallocate (torsp2)
      allocate (qorb(n))
      allocate (worb(n))
      allocate (emorb(n))
      allocate (bkpi(nbond))
      allocate (blpi(nbond))
      allocate (kslope(nbond))
      allocate (lslope(nbond))
      allocate (torsp2(ntors))
c
c     assign the values characteristic of the piatom types
c
      do i = 1, norbit
         j = iorbit(i)
         jt = type(j)
         qorb(j) = electron(jt)
         worb(j) = ionize(jt) / evolt
         emorb(j) = repulse(jt) / evolt
      end do
c
c     assign parameters for all bonds between piatoms;
c     store the original bond lengths and force constants
c
      do i = 1, nbpi
         j = ibpi(1,i)
         ia = iorbit(ibpi(2,i))
         ib = iorbit(ibpi(3,i))
         ita = class(ia)
         itb = class(ib)
         size = 4
         call numeral (ita,pa,size)
         call numeral (itb,pb,size)
         if (ita .le. itb) then
            pt = pa//pb
         else
            pt = pb//pa
         end if
c
c     make a check for bonds contained inside small rings
c
         iring = 0
         if (use_ring) then
            call chkring (iring,ia,ib,0,0)
            if (iring .eq. 6)  iring = 0
            if (iring.eq.5 .and. npi5.eq.0)  iring = 0
            if (iring.eq.4 .and. npi4.eq.0)  iring = 0
            if (iring .eq. 3)  iring = 0
         end if
c
c     assign conjugated bond parameters for each pibond
c
         if (iring .eq. 0) then
            do k = 1, npi
               if (kpi(k) .eq. pt) then
                  bkpi(j) = bk(j)
                  blpi(j) = bl(j)
                  kslope(j) = sslope(k)
                  lslope(j) = tslope(k)
                  goto 170
               end if
            end do
c
c     assign bond parameters for 5-membered ring pibonds
c
         else if (iring .eq. 5) then
            do k = 1, npi5
               if (kpi5(k) .eq. pt) then
                  bkpi(j) = bk(j)
                  blpi(j) = bl(j)
                  kslope(j) = sslope5(k)
                  lslope(j) = tslope5(k)
                  goto 170
               end if
            end do
c
c     assign bond parameters for 4-membered ring pibonds
c
         else if (iring .eq. 4) then
            do k = 1, npi4
               if (kpi4(k) .eq. pt) then
                  bkpi(j) = bk(j)
                  blpi(j) = bl(j)
                  kslope(j) = sslope4(k)
                  lslope(j) = tslope4(k)
                  goto 170
               end if
            end do
         end if
c
c     warning if suitable conjugated pibond parameters not found
c
         abort = .true.
         if (header) then
            header = .false.
            write (iout,150)
  150       format (/,' Undefined Conjugated Pibond Parameters :',
     &              //,' Type',13x,'Atom Names',11x,
     &                 'Atom Classes',/)
         end if
         label = 'Pibond'
         if (iring .eq. 5)  label = '5-Ring'
         if (iring .eq. 4)  label = '4-Ring'
         write (iout,160)  label,ia,name(ia),ib,name(ib),ita,itb
  160    format (1x,a6,5x,i6,'-',a3,i6,'-',a3,7x,2i5)
  170    continue
      end do
c
c     store original 2-fold torsional constants across pibonds
c
      do i = 1, ntors
         torsp2(i) = tors2(1,i)
      end do
      return
      end
