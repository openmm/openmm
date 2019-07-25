c
c
c     ##############################################################
c     ##  COPYRIGHT (C) 2009 by Chuanjie Wu & Jay William Ponder  ##
c     ##                   All Rights Reserved                    ##
c     ##############################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  program valence  --  derive valence force field parameters  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "valence" refines force field parameters for valence terms based
c     on a quantum mechanical optimized structure and frequencies
c
c
      program valence
      use sizes
      use atoms
      use files
      use inform
      use iounit
      use keys
      use linmin
      use output
      use potent
      use qmstuf
      use valfit
      implicit none
      integer i,nvar,next
      integer mode,length
      real*8 minimum,grdmin
      real*8 valrms,value
      real*8 valfit1
      real*8, allocatable :: xx(:)
      logical exist,query
      logical doguess
      logical dotarget
      logical dofit
      character*20 keyword
      character*240 record
      character*240 string
      character*240 xyzfile
      external valfit1
      external optsave
c
c
c     initialization of the various modes of operation
c
      call initial
      fit_bond = .true.
      fit_angle = .true.
      fit_strbnd = .false.
      fit_urey = .false.
      fit_opbend = .false.
      fit_tors = .false.
      fit_force = .false.
      fit_struct = .false.
      doguess = .false.
      dotarget = .false.
      dofit = .false.
c
c     find out which valence term protocol is to be performed
c
      mode = 0
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  mode
         query = .false.
      end if
   10 continue
      if (query) then
         write (iout,20)
   20    format (/,' The TINKER Valence Parameter Utility Can :',
     &           //,4x,'(1) Set Initial Values for Valence Parameters',
     &           /,4x,'(2) Compare QM and MM Vibrational Frequencies',
     &           /,4x,'(3) Force Fit of Parameters to QM Results',
     &           /,4x,'(4) Structure Fit of Parameters to QM Results')
         do while (mode.lt.1 .or. mode.gt.4)
            mode = 0
            write (iout,30)
   30       format (/,' Enter the Number of the Desired Choice :  ',$)
            read (input,40,err=50,end=50)  mode
   40       format (i10)
   50       continue
         end do
      end if
      if (mode .eq. 1) then
         doguess = .true.
      else if (mode .eq. 2) then
         dotarget = .true.
      else if (mode .eq. 3) then
         dotarget = .true.
         dofit = .true.
         fit_force = .true.
      else if (mode .eq. 4) then
         dotarget = .true.
         dofit = .true.
         fit_struct = .true.
      end if
c
c     read the Cartesian coordinates and connectivity info
c
      call getxyz
      xyzfile = filename
      length = leng
c
c     read structure and vibrational data from Gaussian output
c
      call readgau
      filename = xyzfile
      leng = length
      call getkey
c
c     assign estimated values to the valence parameters
c
      if (doguess) then
         call attach
         call bonds
         call angles
         call torsions
         call field
         call katom
         call valguess
      else
         call mechanic
      end if
c
c     get control parameters and target values from keyfile
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:9) .eq. 'FIT-BOND ') then
            fit_bond = .true.
         else if (keyword(1:9) .eq. 'FIX-BOND ') then
            fit_bond = .false.
         else if (keyword(1:10) .eq. 'FIT-ANGLE ') then
            fit_angle = .true.
         else if (keyword(1:10) .eq. 'FIX-ANGLE ') then
            fit_angle = .false.
         else if (keyword(1:11) .eq. 'FIT-STRBND ') then
            fit_strbnd = .true.
         else if (keyword(1:11) .eq. 'FIX-STRBND ') then
            fit_strbnd = .false.
         else if (keyword(1:9) .eq. 'FIT-UREY ') then
            fit_urey = .true.
         else if (keyword(1:9) .eq. 'FIX-UREY ') then
            fit_urey = .false.
         else if (keyword(1:11) .eq. 'FIT-OPBEND ') then
            fit_opbend = .true.
         else if (keyword(1:11) .eq. 'FIX-OPBEND ') then
            fit_opbend = .false.
         else if (keyword(1:12) .eq. 'FIT-TORSION ') then
            fit_tors = .true.
         else if (keyword(1:12) .eq. 'FIX-TORSION ') then
            fit_tors = .false.
         end if
      end do
c
c     try to increase robustness of polarization calculations
c
      if (dofit .and. use_polar)  stpmax = 1.0d0
c
c     perform dynamic allocation of some local arrays
c
      allocate (xx(35*n))
c
c     comparison of QM and TINKER structure and frequencies
c
      if (dotarget) then
         if (.not. dofit) then
            do i = 1, n
               x(i) = gx(i)
               y(i) = gy(i)
               z(i) = gz(i)
            end do
            value = valrms (1)
c
c     optimize the valence term force field parameters
c
         else
            call prmvar (nvar,xx)
            value = valrms (1)
            grdmin = -1.0d0
            call nextarg (string,exist)
            if (exist)  read (string,*,err=60,end=60)  grdmin
   60       continue
            if (grdmin .le. 0.0d0) then
               write (iout,70)
   70          format (/,' Enter RMS Gradient Termination Criterion',
     &                    ' [0.01] :  ',$)
               read (input,80)  grdmin
   80          format (f20.0)
            end if
            if (grdmin .le. 0.0d0)  grdmin = 0.01d0
            coordtype = 'NONE'
            call ocvm (nvar,xx,minimum,grdmin,valfit1,optsave)
            call varprm (nvar,xx,0,0.0d0)
            call prmvar (nvar,xx)
            value = valrms (1)
            call prtval
         end if
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (xx)
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine valguess  --  estimate valence parameter values  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "valguess" sets approximate valence parameter values based on
c     quantum mechanical structure and frequency data
c
c
      subroutine valguess
      use sizes
      use angbnd
      use atomid
      use atoms
      use bndstr
      use iounit
      use kangs
      use kbonds
      use kopbnd
      use kstbnd
      use ktorsn
      use kurybr
      use kvdws
      use math
      use opbend
      use qmstuf
      use strbnd
      use tors
      use urey
      use valfit
      use vdwpot
      implicit none
      integer i,j,k
      integer size,number
      integer ia,ib,ic,id
      integer iia,iib,isba,isbb
      integer ita,itb,itc,itd
      integer iva,ivb,ivc
      integer iita,iitb
      integer nv,nb,na
      integer nsb,nop,nt
      integer vnum(maxtyp)
      integer, allocatable :: nequiv(:)
      real*8 xab,yab,zab
      real*8 xcb,ycb,zcb
      real*8 xac,yac,zac
      real*8 rab2,rcb2
      real*8 cosine,dot
      real*8 bndguess
      real*8 angguess
      real*8 uryguess
      real*8 opbguess
      logical done
      character*4 pa,pb,pc,pd
      character*8 ptb
      character*12 pta
      character*16 ptt
c
c
c     check the number of atoms in QM output and TINKER xyz file
c
      if (n .ne. ngatom) then
         write (iout,10)
   10    format (/,' VALENCE  --  The Number of Atoms is Not',
     &              ' Consistent')
         call fatal
      end if
c
c     perform dynamic allocation of some global arrays
c
      if (.not. allocated(bl))  allocate (bl(nbond))
      if (.not. allocated(anat))  allocate (anat(nangle))
c
c     assign initial values to van der Waals parameters
c
      nv = 0
      do i = 1, n
         ita = class(i)
         if (vdwindex .eq. 'TYPE')  ita = type(i)
         done = .false.
         if (i .gt. 1) then
            do j = 1, nv
               if (ita .eq. vnum(j))  done = .true.
            end do
         end if
         if (.not. done) then
            nv = nv + 1
            vnum(nv) = ita
            call vdwguess (i,rad(ita),eps(ita),reduct(ita))
         end if
      end do
c
c     print the initial van der Waals parameter values
c
      if (nv .gt. 0) then
         write (iout,20)
   20    format (/,' Estimated van der Waals Parameters :',/)
      end if
      do i = 1, nv
         ia = vnum(i)
         if (reduct(ia) .eq. 0) then
            write (iout,30)  ia,rad(ia),eps(ia)
   30       format (' vdw',7x,i5,10x,f10.3,f11.4)
         else
            write (iout,40)  ia,rad(ia),eps(ia),reduct(ia)
   40       format (' vdw',7x,i5,10x,f10.3,f11.4,f9.2)
         end if
      end do
c
c     find and store the unique bond stretches in the system
c
      nb = 0
      do i = 1, nbond
         ia = ibnd(1,i)
         ib = ibnd(2,i)
         ita = class(ia)
         itb = class(ib)
         size = 4
         call numeral (ita,pa,size)
         call numeral (itb,pb,size)
         if (ita .le. itb) then
            ptb = pa//pb
         else
            ptb = pb//pa
         end if
         done = .false.
         do j = 1, nb
            if (ptb .eq. kb(j))  done = .true.
         end do
         if (.not. done) then
            nb = nb + 1
            kb(nb) = ptb
         end if
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (nequiv(4*n))
c
c     assign initial values to bond stretch parameters
c
      k = 0
      do i = 1, nbond
         ia = ibnd(1,i)
         ib = ibnd(2,i)
         ita = class(ia)
         itb = class(ib)
         size = 4
         call numeral (ita,pa,size)
         call numeral (itb,pb,size)
         if (ita .le. itb) then
            ptb = pa//pb
         else
            ptb = pb//pa
         end if
         xab = gx(ia) - gx(ib)
         yab = gy(ia) - gy(ib)
         zab = gz(ia) - gz(ib)
         bl(i) = sqrt(xab*xab + yab*yab + zab*zab)
         done = .false.
         do j = 1, k
            if (ptb .eq. kb(j)) then
               done = .true.
               blen(j) = blen(j) + bl(i)
               nequiv(j) = nequiv(j) + 1
            end if
         end do
         if (.not. done) then
            k = k + 1
            bcon(k) = bndguess (ia,ib)
            blen(k) = bl(i)
            nequiv(k) = 1
         end if
      end do
c
c     print the initial bond stretch parameter values
c
      if (nb .gt. 0) then
         write (iout,50)
   50    format (/,' Estimated Bond Stretching Parameters :',/)
      end if
      do i = 1, nb
         blen(i) = blen(i) / dble(nequiv(i))
         ptb = kb(i)
         ia = number(ptb(1:4))
         ib = number(ptb(5:8))
         write (iout,60)  ia,ib,bcon(i),blen(i)
   60    format (' bond',6x,2i5,5x,f10.1,f11.4)
      end do
c
c     find and store the unique angle bends in the system
c
      na = 0
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
            pta = pa//pb//pc
         else
            pta = pc//pb//pa
         end if
         done = .false.
         do j = 1, na
            if (pta .eq. ka(j))  done = .true.
         end do
         if (.not. done) then
            na = na + 1
            ka(na) = pta
         end if
      end do
c
c     assign initial values to angle bend parameters
c
      k = 0
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
            pta = pa//pb//pc
         else
            pta = pc//pb//pa
         end if
         xab = gx(ia) - gx(ib)
         yab = gy(ia) - gy(ib)
         zab = gz(ia) - gz(ib)
         xcb = gx(ic) - gx(ib)
         ycb = gy(ic) - gy(ib)
         zcb = gz(ic) - gz(ib)
         rab2 = xab*xab + yab*yab + zab*zab
         rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
         if (rab2.ne.0.0d0 .and. rcb2.ne.0.0d0) then
            dot = xab*xcb + yab*ycb + zab*zcb
            cosine = dot / sqrt(rab2*rcb2)
            cosine = min(1.0d0,max(-1.0d0,cosine))
            anat(i) = radian * acos(cosine)
         end if
         done = .false.
         do j = 1, k
            if (pta .eq. ka(j)) then
               done = .true.
               ang(1,j) = ang(1,j) + anat(i)
               nequiv(j) = nequiv(j) + 1
            end if
         end do
         if (.not. done) then
            k = k + 1
            acon(k) = angguess (ia,ib,ic)
            ang(1,k) = anat(i)
            nequiv(k) = 1
         end if
      end do
c
c     print the initial angle bend parameter values
c
      if (na .gt. 0) then
         write(iout,70)
   70    format(/,' Estimated Angle Bending Parameters :',/)
      end if
      do i = 1, na
         ang(1,i) = ang(1,i) / dble(nequiv(i))
         pta = ka(i)
         ia = number(pta(1:4))
         ib = number(pta(5:8))
         ic = number(pta(9:12))
         write (iout,80)  ia,ib,ic,acon(i),ang(1,i)
   80    format (' angle',5x,3i5,f10.2,f11.2)
      end do
c
c     assign initial values to stretch-bend parameters
c
      nsb = 0
      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         ita = class(ia)
         itb = class(ib)
         itc = class(ic)
         iva = valence(ia)
         ivb = valence(ib)
         ivc = valence(ic)
         if (iva.gt.1 .or. ivc.gt.1) then
            size = 4
            call numeral (ita,pa,size)
            call numeral (itb,pb,size)
            call numeral (itc,pc,size)
            if (ita .le. itc) then
               pta = pa//pb//pc
            else
               pta = pc//pb//pa
            end if
            done = .false.
            do j = 1, nsb
               if (pta .eq. ksb(j))  done = .true.
            end do
            if (.not. done) then
               nsb = nsb + 1
               ksb(nsb) = pta
               if (ita .le. itc) then
                  call sbguess (ia,ib,ic,stbn(1,nsb),stbn(2,nsb))
               else
                  call sbguess (ic,ib,ia,stbn(1,nsb),stbn(2,nsb))
               end if
            end if
         end if
      end do
c
c     print the initial stretch-bend parameter values
c
      if (nsb .gt. 0) then
         write (iout,90)
   90    format (/,' Estimated Stretch-Bend Parameters :',/)
      end if
      do i = 1, nsb
         pta = ksb(i)
         ia = number(pta(1:4))
         ib = number(pta(5:8))
         ic = number(pta(9:12))
         write (iout,100)  ia,ib,ic,stbn(1,i),stbn(2,i)
  100    format (' strbnd',4x,3i5,f10.2,f11.2)
      end do
c
c     assign initial values to Urey-Bradley parameters
c
      k = 0
      do i = 1, nurey
         ia = iury(1,i)
         ib = iury(2,i)
         ic = iury(3,i)
         ita = class(ia)
         itb = class(ib)
         itc = class(ic)
         size = 4
         call numeral (ita,pa,size)
         call numeral (itb,pb,size)
         call numeral (itc,pc,size)
         if (ita .le. itc) then
            pta = pa//pb//pc
         else
            pta = pc//pb//pa
         end if
         xac = gx(ia) - gx(ic)
         yac = gy(ia) - gy(ic)
         zac = gz(ia) - gz(ic)
         ul(i) = sqrt(xac*xac + yac*yac + zac*zac)
         done = .false.
         do j = 1, k
            if (pta .eq. ku(j)) then
               done = .true.
               dst13(j) = dst13(j) + ul(i)
               nequiv(j) = nequiv(j) + 1
            end if
         end do
         if (.not. done) then
            k = k + 1
            ucon(k) = uryguess (ia,ib,ic)
            dst13(k) = ul(i)
            nequiv(k) = 1
         end if
      end do
c
c     print the initial Urey-Bradley parameter values
c
      if (nurey .gt. 0) then
         write (iout,110)
  110    format (/,' Estimated Urey-Bradley Parameters :',/)
      end if
      do i = 1, nsb
         dst13(i) = dst13(i) / dble(nequiv(i))
         pta = ku(i)
         ia = number(pta(1:4))
         ib = number(pta(5:8))
         ic = number(pta(9:12))
         write (iout,120)  ia,ib,ic,ucon(i),dst13(i)
  120    format (' ureybrad',2x,3i5,f10.1,f11.4)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (nequiv)
c
c     assign initial values to out-of-plane bend parameters
c
      nop = 0
      do i = 1, nbond
         ia = ibnd(1,i)
         ib = ibnd(2,i)
         ic = 0
         id = 0
         iva = valence(ia)
         ivb = valence(ib)
         if (iva.eq.3 .or. ivb.eq.3) then
            ita = class(ia)
            itb = class(ib)
            itc = 0
            itd = 0
            size = 4
            call numeral (ita,pa,size)
            call numeral (itb,pb,size)
            call numeral (itc,pc,size)
            call numeral (itd,pd,size)
            if (iva .eq. 3) then
               ptt = pb//pa//pc//pd
               isba = ia
               isbb = ib
            else
               ptt = pa//pb//pc//pd
               isba = ib
               isbb = ia
            end if
            if (atomic(isba) .eq. 6) then
               done = .false.
               do j = 1, nop
                  if (ptt .eq. kopb(j))  done = .true.
               end do
               if (.not. done) then
                  nop = nop + 1
                  kopb(nop) = ptt
                  opbn(nop) = opbguess (isba,isbb,ic,id)
                  do j = i+1, nbond
                     iia = ibnd(1,j)
                     iib = ibnd(2,j)
                     if (iia.eq.isba .or. iib.eq.isba) then
                        iita = class(iia)
                        iitb = class(iib)
                        size = 4
                        call numeral (iita,pa,size)
                        call numeral (iitb,pb,size)
                        if (iia .eq. isba) then
                           ptt = pb//pa//pc//pd
                        else if (iib .eq. isba) then
                           ptt = pa//pb//pc//pd
                        end if
                        done = .false.
                        do k = 1, nop
                           if (ptt .eq. kopb(k))  done = .true.
                        end do
                        if (.not. done) then
                           nop = nop + 1
                           kopb(nop) = ptt
                           if (iia .eq. isba) then
                              opbn(nop) = opbguess (iia,iib,ic,id)
                           else if (iib .eq. isba) then
                              opbn(nop) = opbguess (iib,iia,ic,id)
                           end if
                        end if
                     end if
                  end do
               end if
            else if (atomic(isba) .eq. 7) then
               if (valence(isbb).eq.3 .and. atomic(isbb).eq.6) then
                  nop = nop + 1
                  kopb(nop) = ptt
                  opbn(nop) = opbguess (isba,isbb,ic,id)
                  do j = 1, nbond
                     if (j.ne.i .and. (ibnd(1,j).eq.isba
     &                               .or. ibnd(2,j).eq.isba)) then
                        if (ibnd(1,j) .eq. isba) then
                           iia = ibnd(2,j)
                        else
                           iia = ibnd(1,j)
                        end if
                        size = 4
                        call numeral (class(isba),pa,size)
                        call numeral (class(iia),pb,size)
                        ptt = pb//pa//pc//pd
                        done = .false.
                        do k = 1, nop
                           if (ptt .eq. ksb(k))  done = .true.
                        end do
                        if (.not. done) then
                           nop = nop + 1
                           kopb(nop) = ptt
                           opbn(nop) = opbguess (isba,iia,ic,id)
                        end if
                     end if
                  end do
               end if
            end if
         end if
      end do
c
c     print the initial out-of-plane bend parameter values
c
      if (nop .gt .0) then
         write (iout,130)
  130    format (/,' Estimated Out-of-Plane Parameters :',/)
      end if
      do i = 1, nop
         ptt = kopb(i)
         ia = number(ptt(1:4))
         ib = number(ptt(5:8))
         ic = number(ptt(9:12))
         id = number(ptt(13:16))
         write (iout,140)  ia,ib,ic,id,opbn(i)
  140    format (' opbend',4x,4i5,6x,f10.2)
      end do
c
c     assign initial values to torsional parameters
c
      nt = 0
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
         if (itb .le. itc) then
            ptt = pa//pb//pc//pd
         else
            ptt = pd//pc//pb//pa
         end if
         done = .false.
         do j = 1, nt
            if (ptt .eq. kt(j))  done = .true.
         end do
         if (.not. done) then
            nt = nt + 1
            kt(nt) = ptt
            call torguess (ia,ib,ic,id,t1(1,nt),t2(1,nt),t3(1,nt))
         end if
      end do
c
c     print the initial torsional parameter values
c
      if (nt .gt. 0) then
         write (iout,150)
  150    format (/,' Estimated Torsional Parameters :'/)
      end if
      do i = 1, nt
         ptt = kt(i)
         ia = number(ptt(1:4))
         ib = number(ptt(5:8))
         ic = number(ptt(9:12))
         id = number(ptt(13:16))
         write (iout,160)  ia,ib,ic,id,t1(1,i),t2(1,i),t3(1,i)
  160    format (' torsion',3x,4i5,3x,f8.3,' 0.0 1',f8.3,
     &              ' 180.0 2',f8.3,' 0.0 3')
      end do
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine vdwguess  --  estimate van der Waals parameters  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "vdwguess" sets initial VDW parameters based on atom type
c     and connected atoms
c
c
      subroutine vdwguess (ia,rad,eps,reduce)
      use sizes
      use atomid
      use couple
      use math
      use vdwpot
      implicit none
      integer i,j,k,ia
      integer ita,itb
      integer iva,ivb
      real*8 rad,eps,reduce
c
c
c     set default value for radius, well depth and reduction factor
c
      rad = 1.0d0
      eps = 0.1d0
      reduce = 0.0d0
c
c     get atomic number and valence for the atom and its neighbor
c
      ita = atomic(ia)
      iva = valence(ia)
      itb = 0
      ivb = 0
      do i = 1, n12(ia)
         j = i12(i,ia)
         k = atomic(j)
         if (k .gt. itb) then
            itb = k
            ivb = valence(j)
         end if
      end do
c
c     assign specific values based on atom type and connectivity
c
      if (ita .eq. 1) then
         if (itb .eq. 6) then
            if (ivb .eq. 3) then
               rad = 2.980d0
               eps = 0.0260d0
               reduce = 0.92d0
            else if (ivb .eq. 4) then
               rad = 2.780d0
               eps = 0.0260d0
               reduce = 0.91d0
            else
               rad = 2.780d0
               eps = 0.0260d0
               reduce = 0.91d0
            end if
         else if (itb .eq. 7) then
            rad = 2.700d0
            eps = 0.0200d0
            reduce = 0.91d0
         else if (itb .eq. 8) then
            rad = 2.655d0
            eps = 0.0135d0
            reduce = 0.91d0
         else if (itb .eq. 16) then
            rad = 3.000d0
            eps = 0.0265d0
            reduce = 0.98d0
         else
            rad = 2.980d0
            eps = 0.0260d0
            reduce = 0.92d0
         end if
      else if (ita .eq. 6) then
         if (iva .eq. 3) then
            rad = 3.800d0
            eps = 0.0890d0
         else if (iva .eq. 4) then
            rad = 3.820d0
            eps = 0.1010d0
         else
            rad = 3.820d0
            eps = 0.1010d0
         end if
      else if (ita .eq. 7) then
         if (iva .eq. 3) then
            rad = 3.710d0
            eps = 0.1050d0
         else if (iva .eq. 2) then
            rad = 3.710d0
            eps = 0.1100d0
         else
            rad = 3.710d0
            eps = 0.1050d0
         end if
      else if (ita .eq. 8) then
         if (iva .eq. 1) then
            if (itb .eq. 6) then
               rad = 3.300d0
               eps = 0.1120d0
            else if (itb .eq. 7) then
               rad = 3.300d0
               eps = 0.1120d0
            else if (itb .eq. 15) then
               rad = 3.360d0
               eps = 0.1120d0
            else if (itb .eq. 16) then
               rad = 3.510d0
               eps = 0.1120d0
            else
               rad = 3.300d0
               eps = 0.1120d0
            end if
         else if (iva .eq. 2) then
            if (itb .eq. 15) then
               rad = 3.405d0
               eps = 0.1120d0
            else
               rad = 3.405d0
               eps = 0.1100d0
            end if
         else
            rad = 3.405d0
            eps = 0.1100d0
         end if
      else if (ita .eq. 9) then
         if (iva .eq. 0) then
            rad = 3.400d0
            eps = 0.2500d0
         else if (iva .eq. 1) then
            rad = 3.220d0
            eps = 0.1200d0
         else
            rad = 3.220d0
            eps = 0.1200d0
         end if
      else if (ita .eq. 11) then
         rad = 3.020d0
         eps = 0.260d0
      else if (ita .eq. 12) then
         rad = 2.550d0
         eps = 0.850d0
      else if (ita .eq. 15) then
         rad = 4.450d0
         eps = 0.390d0
      else if (ita .eq. 16) then
         if (iva .eq. 2) then
            rad = 3.910d0
            eps = 0.3850d0
         else if (iva .eq. 3) then
            rad = 3.910d0
            eps = 0.3850d0
         else if (iva .eq. 4) then
            rad = 3.910d0
            eps = 0.3850d0
         else
            rad = 3.910d0
            eps = 0.3850d0
         end if
      else if (ita .eq. 17) then
         if (iva .eq. 0) then
            rad = 4.130d0
            eps = 0.340d0
         else if (iva .eq. 1) then
            rad = 4.130d0
            eps = 0.340d0
         else
            rad = 4.130d0
            eps = 0.340d0
         end if
      else if (ita .eq. 19) then
         rad = 3.710d0
         eps = 0.3500d0
      else if (ita .eq. 20) then
         rad = 3.150d0
         eps = 1.6000d0
      else if (ita .eq. 35) then
         if (iva .eq. 0) then
            rad = 4.380d0
            eps = 0.4300d0
         else if (iva .eq. 1) then
            rad = 4.380d0
            eps = 0.4300d0
         else
            rad = 4.380d0
            eps = 0.4300d0
         end if
      else if (ita .eq. 53) then
         if (iva .eq. 0) then
            rad = 4.660d0
            eps = 0.520d0
         else if (iva .eq. 1) then
            rad = 4.660d0
            eps = 0.520d0
         else
            rad = 4.660d0
            eps = 0.520d0
         end if
      end if
c
c     scale the vdw parameters to the desired units
c
      if (radsiz .eq. 'RADIUS')  rad = 0.5d0 * rad
      if (radsiz .eq. 'SIGMA')  rad = rad / twosix
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  function bndguess  --  estimate bond stretch parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "bndguess" sets approximate bond stretch force constants based
c     on atom type and connected atoms
c
c
      function bndguess (ia,ib)
      use sizes
      use atomid
      use bndpot
      implicit none
      integer ia,ib,tmp
      integer ita,itb
      integer iva,ivb
      real*8 bndguess
c
c
c     get the atomic number and valence of each atom
c
      ita = atomic(ia)
      itb = atomic(ib)
      iva = valence(ia)
      ivb = valence(ib)
c
c     reverse the atom order based on atomic number
c
      if (ita .gt. itb) then
         tmp = ita
         ita = itb
         itb = tmp
         tmp = iva
         iva = ivb
         ivb = tmp
      end if
c
c     assign estimated bond stretch force constants
c
      if (ita .eq. 1) then
         if (itb .eq. 6) then
            if (ivb .eq. 3) then
               bndguess = 410.0d0
            else if (ivb .eq. 4) then
               bndguess = 400.0d0
            else
               bndguess = 400.0d0
            end if
         else if (itb .eq. 7) then
            bndguess = 520.0d0
         else if (itb .eq. 8) then
            bndguess = 560.0d0
         else if (itb .eq. 9) then
            bndguess = 500.0d0
         else if (itb .eq. 14) then
            bndguess = 200.0d0
         else if (itb .eq. 15) then
            bndguess = 230.0d0
         else if (itb .eq. 16) then
            bndguess = 260.0d0
         else
            bndguess = 300.0d0
         end if
      else if (ita .eq. 6) then
         if (itb .eq. 6) then
            if (iva.eq.3 .and. ivb.eq.3) then
               bndguess = 680.0d0
            else if (iva.eq.4 .or. ivb.eq.4) then
               bndguess = 385.0d0
            else
               bndguess = 350.0d0
            end if
         else if (itb .eq. 7) then
            if (iva.eq.3 .and. ivb.eq.2) then
               bndguess = 435.0d0
            else if (iva.eq.3 .and. ivb.eq.3) then
               bndguess = 250.0d0
            else if (iva.eq.4) then
               bndguess = 400.0d0
            else
               bndguess = 450.0d0
            end if
         else if (itb .eq. 8) then
            if (ivb .eq. 1) then
               bndguess = 680.0d0
            else if (ivb .eq. 2) then
               bndguess = 465.0d0
            else
               bndguess = 465.0d0
            end if
         else if (itb .eq. 9) then
            bndguess = 350.0d0
         else if (itb .eq. 14) then
            bndguess = 350.0d0
         else if (itb .eq. 15) then
            bndguess = 350.0d0
         else if (itb .eq. 16) then
            bndguess = 216.0d0
         else if (itb .eq. 17) then
            bndguess = 350.0d0
         else
            bndguess = 450.0d0
         end if
      else if (ita .eq. 7) then
         if (itb .eq. 7) then
            if (iva .eq. 1) then
               bndguess = 1613.0d0
            else if (iva.eq.2 .and. ivb.eq.2) then
               bndguess = 950.0d0
            else
               bndguess = 850.0d0
            end if
         else if (itb .eq. 8) then
            if (ivb .eq. 1 ) then
               bndguess = 900.0d0
            else
               bndguess = 750.0d0
            end if
         else if (itb .eq. 14) then
            bndguess = 450.0d0
         else if (itb .eq. 15) then
            bndguess = 500.0d0
         else if (itb .eq. 16) then
            bndguess = 550.0d0
         else
            bndguess = 600.0d0
         end if
      else if (ita .eq. 8) then
         if (itb .eq. 8) then
            bndguess = 750.0d0
         else if (itb .eq. 14) then
            bndguess = 500.0d0
         else if (itb .eq. 15) then
            if (iva .eq. 2) then
               bndguess = 450.0d0
            else if (iva .eq. 1) then
               bndguess = 775.0d0
            else
               bndguess = 450.0d0
            end if
         else if (itb .eq. 16) then
            bndguess = 606.0d0
         else if (itb .eq. 17) then
            bndguess = 500.0d0
         else
            bndguess = 600.0d0
         end if
      else if (ita .eq. 14) then
         if (itb .eq. 14) then
            bndguess = 400.0d0
         else if (itb .eq. 15) then
            bndguess = 450.0d0
         else if (itb .eq. 16) then
            bndguess = 500.0d0
         else if (itb .eq. 17) then
            bndguess = 650.0d0
         else
            bndguess = 450.0d0
         end if
      else if (ita .eq. 16) then
         if (itb .eq. 16) then
            bndguess = 188.0d0
         else
            bndguess = 250.0d0
         end if
      else if (ita .eq. 17) then
         bndguess = 300.0d0
      else
         bndguess = 350.0d0
      end if
c
c     scale the force constant to the desired units
c
      bndguess = bndguess / bndunit
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  function angguess  --  estimate angle bending parameters  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "angguess" sets approximate angle bend force constants based
c     on atom type and connected atoms
c
c
      function angguess (ia,ib,ic)
      use sizes
      use atomid
      use angpot
      use math
      implicit none
      integer ia,ib,ic,tmp
      integer ita,itb,itc
      integer iva,ivb,ivc
      real*8 angguess
c
c
c     get the atomic number and valence of each atom
c
      ita = atomic(ia)
      itb = atomic(ib)
      itc = atomic(ic)
      iva = valence(ia)
      ivb = valence(ib)
      ivc = valence(ic)
c
c     resort ja,jb,jc based on the atomic orders
c
      if (ita .gt. itc) then
         tmp = ita
         ita = itc
         itc = tmp
         tmp = iva
         iva = ivc
         ivc = tmp
      end if
c
c     assign estimated angle bend force constants
c
      if (itb .eq. 6) then
         if (ita .eq. 1) then
            if (ivb .eq. 4) then
               if (itc .eq. 1) then
                  angguess = 34.50d0
               else if (itc .eq. 6) then
                  angguess = 38.0d0
               else if (itc .eq. 7) then
                  angguess = 50.60d0
               else if (itc .eq. 8) then
                  angguess = 51.50d0
               else if (itc .eq. 9) then
                  angguess = 50.0d0
               else
                  angguess = 35.0d0
               end if
            else if (ivb .eq. 3) then
               angguess = 32.00d0
            else
               angguess = 32.00d0
            end if
         else if (ita .eq. 6) then
            if (ivb .eq. 4) then
               if (itc .eq. 6) then
                  angguess = 60.00d0
               else if (itc .eq. 7) then
                  angguess = 80.00d0
               else if (itc .eq. 8) then
                  angguess = 88.00d0
               else if (itc .eq. 9) then
                  angguess = 89.00d0
               else if (itc .eq. 14) then
                  angguess = 65.00d0
               else if (itc .eq. 15) then
                  angguess = 60.00d0
               else if (itc .eq. 16) then
                  angguess = 53.20d0
               else if (itc .eq. 17) then
                  angguess = 55.00d0
               else
                  angguess = 50.00d0
               end if
            else if (ivb .eq. 3) then
               angguess = 60.00d0
            else
               angguess = 60.00d0
            end if
         else if (ita .eq. 8) then
            if (ivb .eq. 4) then
               if (itc .eq. 8) then
                  angguess = 65.00d0
               else if (itc .eq. 9) then
                  angguess = 65.00d0
               else if (itc .eq. 15) then
                  angguess = 60.00d0
               else if (itc .eq. 16) then
                  angguess = 65.00d0
               else
                  angguess = 65.00d0
               end if
            else if (ivb .eq. 3) then
               angguess = 50.00d0
            else
               angguess = 60.00d0
            end if
         else
            angguess = 60.00d0
         end if
      else if (itb .eq. 8) then
         if (ita .eq. 1) then
            if (itc .eq. 1) then
               angguess = 34.05d0
            else if (itc .eq. 6) then
               angguess = 65.00d0
            else
               angguess = 60.00d0
            end if
         else if (ita .eq. 6) then
            if (itc .eq. 6) then
               angguess = 88.50d0
            else if (itc .eq. 8) then
               if (iva.eq.1 .or. ivc.eq.1) then
                  angguess = 122.30d0
               else
                  angguess = 85.00d0
               end if
            else if (itc .eq. 15) then
               angguess = 80.30d0
            else
               angguess = 80.0d0
            end if
         else
            angguess = 80.0d0
         end if
      else if (itb .eq. 15) then
         if (ita .eq. 1) then
            angguess = 30.0d0
         else if (ita .eq. 6) then
            if (itc .eq. 6) then
               angguess = 75.00d0
            else if (itc .eq. 8) then
               angguess = 80.00d0
            else
               angguess = 75.00d0
            end if
         else if (ita .eq. 8) then
            if (itc .eq. 8) then
               if (iva.eq.1 .and. ivc.eq.1) then
                  angguess = 89.88d0
               else if (iva.eq.1 .or. ivc.eq.1) then
                  angguess = 75.86d0
               else
                  angguess = 65.58d0
               end if
            else
               angguess = 70.00d0
            end if
         else
            angguess = 75.00d0
         end if
      else if (itb .eq. 16) then
         if (ita .eq. 1) then
            angguess = 30.00d0
         else if (ita .eq. 6) then
            if (itc .eq. 16) then
               angguess = 72.00d0
            else
               angguess = 80.00d0
            end if
         else if (ita .eq. 8) then
            if (itc .eq. 8) then
               if (iva.eq.1 .and. ivc.eq.1) then
                  angguess = 168.00d0
               else if (iva.eq.1 .or. ivc.eq.1) then
                  angguess = 85.00d0
               else
                  angguess = 80.00d0
               end if
            else if (itc .eq. 16) then
               angguess = 75.00d0
            else
               angguess = 75.00d0
            end if
         else
            angguess = 75.00d0
         end if
      else if (ita .eq. 1) then
         angguess = 35.00d0
      else
         angguess = 65.00d0
      end if
c
c     scale the force constant to the desired units
c
      angguess = angguess / (angunit*radian**2)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine sbguess  --  estimate stretch-bend parameters  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "sbguess" sets approximate stretch-bend force constants based
c     on atom type and connected atoms
c
c
      subroutine sbguess (ia,ib,ic,sb1,sb2)
      use sizes
      use angpot
      use atomid
      use math
      implicit none
      integer ia,ib,ic
      integer ita,itb,itc
      integer iva,ivb,ivc
      real*8 sb1,sb2
c
c
c     get the atomic number and valence of each atom
c
      ita = atomic(ia)
      itb = atomic(ib)
      itc = atomic(ic)
      iva = valence(ia)
      ivb = valence(ib)
      ivc = valence(ic)
c
c     set initial stretch-bend parameters
c
      if (ita.eq.1 .and. itc.eq.1) then
         sb1 = 0.0d0
         sb2 = 0.0d0
      else if (itb .eq. 6) then
         if (ita .eq. 1 ) then
            sb1 = 11.50d0
            sb2 = 18.70d0
         else if (itc .eq. 1) then
            sb1 = 18.70d0
            sb2 = 11.50d0
         else
            sb1 = 18.70d0
            sb2 = 18.70d0
         end if
      else if (itb .eq. 6) then
         if (ita .eq. 1) then
            sb1 = 4.50d0
            sb2 = 12.95d0
         else if (itc .eq. 1) then
            sb1 = 12.95d0
            sb2 = 4.50d0
         else
            sb1 = 14.40d0
            sb2 = 14.40d0
         end if
      else if (itb .eq. 7) then
         if (ivb .ge. 3) then
            if (ita .eq. 1) then
               sb1 = 4.30d0
               sb2 = 7.20d0
            else if (itc .eq. 1) then
               sb1 = 7.20d0
               sb2 = 4.30d0
            else
               sb1 = 7.20d0
               sb2 = 7.20d0
            end if
         else
            if (ita .eq. 1) then
               sb1 = 4.30d0
               sb2 = 14.40d0
            else if (itc .eq. 1) then
               sb1 = 14.40d0
               sb2 = 4.30d0
            else
               sb1 = 14.40d0
               sb2 = 14.40d0
            end if
         end if
      else if (itb .eq. 14) then
         if (ita .eq. 1) then
            sb1 = 8.60d0
            sb2 = 14.40d0
         else if (itc .eq. 1) then
            sb1 = 14.40d0
            sb2 = 8.60d0
         else
            sb1 = 14.40d0
            sb2 = 14.40d0
         end if
      else if (itb .eq. 15) then
         if (ivb .eq. 4) then
            if (ita .eq. 1) then
               sb1 = 14.40d0
               sb2 = 14.40d0
            else if (itc .eq. 1) then
               sb1 = 14.40d0
               sb2 = 14.40d0
            else
               sb1 = 14.40d0
               sb2 = 14.40d0
            end if
         else
            if (ita .eq. 1) then
               sb1 = 8.60d0
               sb2 = 8.60d0
            else if (itc .eq. 1) then
               sb1 = 8.60d0
               sb2 = 8.60d0
            else
               sb1 = 8.60d0
               sb2 = 8.60d0
            end if
         end if
      else if (itb .eq. 16) then
         if (ita .eq. 1) then
            sb1 = 1.45d0
            sb2 = -5.75d0
         else if (itc .eq. 1) then
            sb1 = -5.75d0
            sb2 = 1.45d0
         else
            sb1 = -5.75d0
            sb2 = -5.75d0
         end if
      else if (ita.eq.1 .and. itc.gt.1) then
         sb1 = -4.50d0
         sb2 = 38.00d0
      else if (ita.gt.1 .and. itc.eq.1) then
         sb1 = 38.00d0
         sb2 = -4.50d0
      else if (ita.gt.1 .and. itc.gt.1) then
         sb1 = 38.00d0
         sb2 = 38.00d0
      else
         sb1 = 38.00d0
         sb2 = 38.00d0
      end if
c
c     scale the force constant to the desired units
c
      sb1 = sb1 / (stbnunit*radian)
      sb2 = sb2 / (stbnunit*radian)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  function uryguess  --  estimate Urey-Bradley parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "uryguess" sets approximate Urey-Bradley force constants
c     based on atom type and connected atoms
c
c
      function uryguess (ia,ib,ic)
      use sizes
      use atomid
      use urypot
      implicit none
      integer ia,ib,ic
      integer ita,itb,itc
      integer iva,ivb,ivc
      real*8 uryguess
c
c
c     get the atomic number and valence of each atom
c
      ita = atomic(ia)
      itb = atomic(ib)
      itc = atomic(ic)
      iva = valence(ia)
      ivb = valence(ib)
      ivc = valence(ic)
c
c     assign estimated out-of-plane parameter values
c
      uryguess = 10.0d0
c
c     scale the force constant to the desired units
c
      uryguess = uryguess / ureyunit
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  function opbguess  --  estimate out-of-plane bend values  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "opbguess" sets approximate out-of-plane bend force constants
c     based on atom type and connected atoms
c
c
      function opbguess (ia,ib,ic,id)
      use sizes
      use angpot
      use atomid
      use math
      implicit none
      integer ia,ib,ic,id
      integer ita,itb
      integer iva,ivb
      real*8 opbguess
c
c
c     get the atomic number and valence of each atom
c
      ita = atomic(ia)
      itb = atomic(ib)
      iva = valence(ia)
      ivb = valence(ib)
c
c     assign estimated out-of-plane parameter values
c
      opbguess = 14.40d0
c
c     scale the force constant to the desired units
c
      opbguess = opbguess / (opbunit*radian**2)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine torguess  --  estimate torsional parameters  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "torguess" set approximate torsion amplitude parameters based
c     on atom type and connected atoms
c
c
      subroutine torguess (ia,ib,ic,id,tf1,tf2,tf3)
      use sizes
      use atomid
      use torpot
      implicit none
      integer ia,ib,ic,id,tmp
      integer ita,itb,itc,itd
      integer iva,ivb,ivc,ivd
      real*8 tf1,tf2,tf3
c
c
c     get the atomic number and valence of each atom
c
      ita = atomic(ia)
      itb = atomic(ib)
      itc = atomic(ic)
      itd = atomic(id)
      iva = valence(ia)
      ivb = valence(ib)
      ivc = valence(ic)
      ivd = valence(id)
c
c     reorder the atoms based on the atomic numbers
c
      if (itb.gt.itc .or. (itb.eq.itc.and.ita.gt.itd)) then
         tmp = itb
         itb = itc
         itc = tmp
         tmp = ivb
         ivb = ivc
         ivc = tmp
         tmp = ita
         ita = itd
         itd = tmp
         tmp = iva
         iva = ivd
         ivd = tmp
      end if
c
c     assign estimated torsional parameter values
c
      tf1 = 0.0d0
      tf2 = 0.0d0
      tf3 = 0.0d0
      if (itb.eq.6 .and. itc.eq.6) then
         if (ita.eq.6 .and. itd.eq.6) then
            if (ivb.eq.3 .and. ivc.eq.3) then
               if (iva.eq.3 .and. ivd.eq.3) then
                  tf1 = -0.335d0
                  tf2 = 2.00d0
                  tf3 = 0.00d0
               else if (iva.eq.3 .and. ivd.eq.4) then
                  tf1 = -0.305d0
                  tf2 = 2.105d0
                  tf3 = 0.00d0
               else if (iva.eq.4 .and. ivd.eq.4 ) then
                  tf1 = 0.00d0
                  tf2 = 4.00d0
                  tf3 = 0.00d0
               end if
            else if (ivb.eq.3 .and. ivc.eq.4) then
               tf1 = -0.40d0
               tf2 = -0.05d0
               tf3 = -0.275d0
            else if (ivb.eq.4 .and. ivc.eq.4) then
               tf1 = 0.09d0
               tf2 = 0.085d0
               tf3 = 0.26d0
            end if
         else if (ita.eq.1 .and. itd.eq.1) then
            if (ivb.eq.3 .and. ivc.eq.3) then
               tf1 = 0.00d0
               tf2 = 2.035d0
               tf3 = 0.00d0
            else
               tf1 = 0.00d0
               tf2 = 0.00d0
               tf3 = 0.15d0
            end if
         else if (ita.eq.1 .and. itd.eq.6) then
            if (ivb.eq.4 .and. ivc.eq.4) then
               tf1 = 0.00d0
               tf2 = 0.00d0
               tf3 = 0.17d0
            else if (ivb.eq.3 .and. ivc.eq.3 .and. ivd.eq.3) then
               tf1 = 0.00d0
               tf2 = 3.05d0
               tf3 = 0.00d0
            else if (ivb.eq.3 .and. ivc.eq.3 .and. ivd.eq.4) then
               tf1 = 0.00d0
               tf2 = 3.05d0
               tf3 = 0.00d0
            else if (ivb.eq.4 .and. ivc.eq.3) then
               tf1 = 0.00d0
               tf2 = 0.00d0
               tf3 = -0.045d0
            end if
         else if (ita.eq.1 .and. itd.eq.7) then
            if (ivb.eq.3 .and. ivc.eq.3) then
               tf1 = -1.575d0
               tf2 = 1.50d0
               tf3 = 0.00d0
            else
               tf1 = 0.00d0
               tf2 = 0.00d0
               tf3 = 0.25d0
            end if
         else if (ita.eq.1 .and. itd.eq.8) then
            tf1 = 0.00d0
            tf2 = 0.00d0
            tf3 = 0.15d0
         else if (ita.eq.6 .and. itd.eq.8) then
            if (ivb.eq.3 .and. ivc.eq.3) then
               tf1 = 0.00d0
               tf2 = 2.235d0
               tf3 = 0.00d0
            else
               tf1 = -0.575d0
               tf2 = 0.00d0
               tf3 = 0.64d0
            end if
         else if (ita.eq.8 .and. itd.eq.8) then
            tf1 = 1.11d0
            tf2 = -0.69d0
            tf3 = -0.59d0
         else if (ivb.eq.3 .and. ivc.eq.3) then
            tf1 = 0.00d0
            tf2 = 1.25d0
            tf3 = 0.00d0
         else
            tf1 = 0.00d0
            tf2 = 0.00d0
            tf3 = 0.15d0
         end if
      else if (itb.eq.6 .and. itc.eq.8) then
         if(ita.eq.1 .and. itd.eq.1) then
            tf1 = 0.00d0
            tf2 = 0.00d0
            tf3 = 0.135d0
         else if (ita.eq.1 .and. itd.eq.6) then
            if (ivc.eq.3 .and. ivd.eq.3) then
               tf1 = 0.00d0
               tf2 = 2.235d0
               tf3 = 0.00d0
            else
               tf1 = 0.00d0
               tf2 = 0.00d0
               tf3 = 0.355d0
            end if
         else if (ita .eq. 1) then
            tf1 = 0.00d0
            tf2 = 0.00d0
            tf3 = 0.375d0
         else if (ita.eq.6 .and. itd.eq.1 .and. ivb.eq.4) then
            tf1 = -0.885d0
            tf2 = 0.115d0
            tf3 = 0.38d0
         else if (ita.eq.6 .and. itd.eq.6) then
            tf1 = 1.00d0
            tf2 = -0.75d0
            tf3 = 0.445d0
         else if (ita.eq.6 .and. itd.eq.1 .and. ivb.eq.3) then
            tf1 = 0.00d0
            tf2 = 1.175d0
            tf3 = 0.00d0
         else if (ivb .eq. 3) then
            tf1 = 0.00d0
            tf2 = 1.25d0
            tf3 = 0.00d0
         else if (ivb .eq. 4) then
            tf1 = 1.00d0
            tf2 = -0.75d0
            tf3 = 0.445d0
         end if
      else if (itb.eq.6 .and. itc.eq.15) then
         tf1 = 0.00d0
         tf2 = 1.25d0
         tf3 = 0.25d0
      else if (itb.eq.6 .and. itc.eq.16) then
         tf1 = 0.00d0
         tf2 = 0.00d0
         tf3 = 0.25d0
      else if (itb.eq.8 .and. itc.eq.15) then
         tf1 = -1.00d0
         tf2 = -0.84d0
         tf3 = -0.40d0
      else if (itb.eq.8 .and. itc.eq.16) then
         tf1 = -0.75d0
         tf2 = -1.00d0
         tf3 = -0.40d0
      else
         tf1 = 0.00d0
         tf2 = 0.50d0
         tf3 = 0.25d0
      end if
c
c     scale the amplitude values to the desired units
c
      tf1 = tf1 / torsunit
      tf2 = tf2 / torsunit
      tf3 = tf3 / torsunit
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  function valrms  --  compute structure & vibration RMSD  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "valrms" evaluates a valence parameter goodness-of-fit error
c     function based on comparison of forces, frequencies, bond
c     lengths and angles to QM results
c
c
      function valrms (prtflg)
      use sizes
      use angbnd
      use atoms
      use atomid
      use bndstr
      use hescut
      use iounit
      use inform
      use kangs
      use kbonds
      use kopbnd
      use kstbnd
      use ktorsn
      use kvdws
      use linmin
      use math
      use minima
      use opbend
      use output
      use qmstuf
      use scales
      use strbnd
      use tors
      use units
      use valfit
      implicit none
      integer i,j,k
      integer m,m1,m2
      integer ia,ib,ic,id
      integer olditer
      integer oldprt,oldwrt
      integer prtflg,ihess
      integer nvar,nfreq
      integer, allocatable :: hindex(:)
      integer, allocatable :: hinit(:,:)
      integer, allocatable :: hstop(:,:)
      real*8 xab,yab,zab
      real*8 xba,yba,zba
      real*8 xcb,ycb,zcb
      real*8 xdc,ydc,zdc
      real*8 bond,gbond
      real*8 angle,gangle
      real*8 factor,grdmin
      real*8 oldstep
      real*8 delta,cosine,sine
      real*8 rab2,rcb,rcb2,rabc
      real*8 xt,yt,zt,xu,yu,zu
      real*8 xtu,ytu,ztu
      real*8 rt2,ru2,rtru
      real*8 minimiz1,minimum
      real*8 valrms,energy
      real*8 bave,brms,bfac
      real*8 aave,arms,afac
      real*8 tave,trms,tfac
      real*8 gave,grms,gfac
      real*8 have,hrms,hfac
      real*8 fave,frms,ffac,fcut
      real*8, allocatable :: xx(:)
      real*8, allocatable :: mass2(:)
      real*8, allocatable :: eigen(:)
      real*8, allocatable :: h(:)
      real*8, allocatable :: matrix(:)
      real*8, allocatable :: derivs(:,:)
      real*8, allocatable :: hdiag(:,:)
      real*8, allocatable :: vects(:,:)
      character*1 axis(3)
      external minimiz1
      external optsave
      data axis  / 'X','Y','Z' /
c
c
c     scale the coordinates of each active atom; use the
c     square root of median eigenvalue of typical Hessian
c
      if (fit_struct) then
         allocate (xx(3*n))
         set_scale = .true.
         nvar = 0
         do i = 1, n
            nvar = nvar + 1
            scale(nvar) = 12.0d0
            xx(nvar) = gx(i) * scale(nvar)
            nvar = nvar + 1
            scale(nvar) = 12.0d0
            xx(nvar) = gy(i) * scale(nvar)
            nvar = nvar + 1
            scale(nvar) = 12.0d0
            xx(nvar) = gz(i) * scale(nvar)
         end do
c
c     make the call to the optimization routine
c
         oldstep = stpmax
         olditer = maxiter
         oldprt = iprint
         oldwrt = iwrite
         stpmax = 0.0d0
         maxiter = 0
         iprint = 0
         iwrite = 0
         grdmin = 0.0001d0
         coordtype = 'CARTESIAN'
         call lbfgs (nvar,xx,minimum,grdmin,minimiz1,optsave)
         coordtype = 'NONE'
         stpmax = oldstep
         maxiter = olditer
         iprint = oldprt
         iwrite = oldwrt
c
c     unscale the final coordinates for active atoms
c
         nvar = 0
         do i = 1, n
            nvar = nvar + 1
            x(i) = xx(nvar) / scale(nvar)
            scale(nvar) = 1.0d0
            nvar = nvar + 1
            y(i) = xx(nvar) / scale(nvar)
            scale(nvar) = 1.0d0
            nvar = nvar + 1
            z(i) = xx(nvar) / scale(nvar)
            scale(nvar) = 1.0d0
         end do
         deallocate (xx)
      end if
c
c     compute the RMS between QM and TINKER bond lengths
c
      bave = 0.0d0
      brms = 0.0d0
      if (fit_struct) then
         if (prtflg.eq.1 .and. nbond.ne.0) then
            write (iout,10)
   10       format (/,' Comparison of Bond Lengths :',
     &              //,6x,'Bond',8x,'Atoms',19x,'QM Bond',
     &                 6x,'MM Bond',8x,'Delta',/)
         end if
         do i = 1, nbond
            ia = ibnd(1,i)
            ib = ibnd(2,i)
            xab = x(ia) - x(ib)
            yab = y(ia) - y(ib)
            zab = z(ia) - z(ib)
            bond = sqrt(xab*xab + yab*yab + zab*zab)
            xab = gx(ia) - gx(ib)
            yab = gy(ia) - gy(ib)
            zab = gz(ia) - gz(ib)
            gbond = sqrt(xab*xab + yab*yab + zab*zab)
            delta = bond - gbond
            bave = bave + abs(delta)
            brms = brms + delta*delta
            if (prtflg .eq. 1) then
               write (iout,20)  i,ia,ib,gbond,bond,delta
   20          format (4x,i5,4x,2i5,13x,3f13.4)
            end if
         end do
         if (nbond .ne. 0)  bave = bave / (dble(nbond))
         if (nbond .ne. 0)  brms = sqrt(brms/dble(nbond))
         if (prtflg.eq.1 .and. nbond.ne.0) then
            write (iout,30)  bave,brms
   30       format (/,4x,'Average Unsigned Difference :',30x,f12.4,
     &              /,4x,'Root Mean Square Deviation :',31x,f12.4)
         end if
      end if
c
c     compute the RMS between QM and TINKER bond angles
c
      aave = 0.0d0
      arms = 0.0d0
      if (fit_struct) then
         if (prtflg.eq.1 .and. nangle.ne.0) then
            write (iout,40)
   40       format (/,' Comparison of Bond Angles :',
     &              //,5x,'Angle',10x,'Atoms',16x,'QM Angle',
     &                 5x,'MM Angle',8x,'Delta',/)
         end if
         do i = 1, nangle
            ia = iang(1,i)
            ib = iang(2,i)
            ic = iang(3,i)
            xab = x(ia) - x(ib)
            yab = y(ia) - y(ib)
            zab = z(ia) - z(ib)
            xcb = x(ic) - x(ib)
            ycb = y(ic) - y(ib)
            zcb = z(ic) - z(ib)
            rab2 = xab*xab + yab*yab + zab*zab
            rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
            rabc = sqrt(rab2 * rcb2)
            if (rabc .ne. 0.0d0) then
               cosine = (xab*xcb + yab*ycb + zab*zcb) / rabc
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
            end if
            xab = gx(ia) - gx(ib)
            yab = gy(ia) - gy(ib)
            zab = gz(ia) - gz(ib)
            xcb = gx(ic) - gx(ib)
            ycb = gy(ic) - gy(ib)
            zcb = gz(ic) - gz(ib)
            rab2 = xab*xab + yab*yab + zab*zab
            rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
            rabc = sqrt(rab2 * rcb2)
            if (rabc .ne. 0.0d0) then
               cosine = (xab*xcb + yab*ycb + zab*zcb) / rabc
               cosine = min(1.0d0,max(-1.0d0,cosine))
               gangle = radian * acos(cosine)
            end if
            delta = angle - gangle
            aave = aave + abs(delta)
            arms = arms + delta*delta
            if (prtflg .eq. 1) then
               write (iout,50)  i,ia,ib,ic,gangle,angle,delta
   50          format (4x,i5,4x,3i5,8x,2f13.2,f13.4)
            end if
         end do
         if (nangle .ne. 0)  aave = aave / (dble(nangle))
         if (nangle .ne. 0)  arms = sqrt(arms/dble(nangle))
         if (prtflg.eq.1 .and. nangle.ne.0) then
            write (iout,60)  aave,arms
   60       format (/,4x,'Average Unsigned Difference :',30x,f12.4,
     &              /,4x,'Root Mean Square Deviation :',31x,f12.4)
         end if
      end if
c
c     compute the RMS between QM and TINKER torsion angles
c
      tave = 0.0d0
      trms = 0.0d0
      if (fit_struct) then
         if (prtflg.eq.1 .and. ntors.ne.0) then
            write (iout,70)
   70       format (/,' Comparison of Torsion Angles :',
     &              //,4x,'Torsion',12x,'Atoms',13x,'QM Angle',
     &                 5x,'MM Angle',8x,'Delta',/)
         end if
         do i = 1, ntors
            ia = itors(1,i)
            ib = itors(2,i)
            ic = itors(3,i)
            id = itors(4,i)
            xba = x(ib) - x(ia)
            yba = y(ib) - y(ia)
            zba = z(ib) - z(ia)
            xcb = x(ic) - x(ib)
            ycb = y(ic) - y(ib)
            zcb = z(ic) - z(ib)
            xdc = x(id) - x(ic)
            ydc = y(id) - y(ic)
            zdc = z(id) - z(ic)
            xt = yba*zcb - ycb*zba
            yt = zba*xcb - zcb*xba
            zt = xba*ycb - xcb*yba
            xu = ycb*zdc - ydc*zcb
            yu = zcb*xdc - zdc*xcb
            zu = xcb*ydc - xdc*ycb
            xtu = yt*zu - yu*zt
            ytu = zt*xu - zu*xt
            ztu = xt*yu - xu*yt
            rt2 = xt*xt + yt*yt + zt*zt
            ru2 = xu*xu + yu*yu + zu*zu
            rtru = sqrt(rt2 * ru2)
            if (rtru .ne. 0.0d0) then
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
               if (sine .lt. 0.0d0)  angle = -angle
            end if
            xba = gx(ib) - gx(ia)
            yba = gy(ib) - gy(ia)
            zba = gz(ib) - gz(ia)
            xcb = gx(ic) - gx(ib)
            ycb = gy(ic) - gy(ib)
            zcb = gz(ic) - gz(ib)
            xdc = gx(id) - gx(ic)
            ydc = gy(id) - gy(ic)
            zdc = gz(id) - gz(ic)
            xt = yba*zcb - ycb*zba
            yt = zba*xcb - zcb*xba
            zt = xba*ycb - xcb*yba
            xu = ycb*zdc - ydc*zcb
            yu = zcb*xdc - zdc*xcb
            zu = xcb*ydc - xdc*ycb
            xtu = yt*zu - yu*zt
            ytu = zt*xu - zu*xt
            ztu = xt*yu - xu*yt
            rt2 = xt*xt + yt*yt + zt*zt
            ru2 = xu*xu + yu*yu + zu*zu
            rtru = sqrt(rt2 * ru2)
            if (rtru .ne. 0.0d0) then
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               gangle = radian * acos(cosine)
               if (sine .lt. 0.0d0)  gangle = -gangle
            end if
            delta = angle - gangle
            if (delta .gt. 180.0d0)  delta = delta - 360.0d0
            if (delta .lt. -180.0d0)  delta = delta + 360.0d0
            tave = tave + abs(delta)
            trms = trms + delta*delta
            if (prtflg .eq. 1) then
               write (iout,80)  i,ia,ib,ic,id,gangle,angle,delta
   80          format (4x,i5,4x,4i5,3x,2f13.2,f13.4)
            end if
         end do
         if (ntors .ne. 0)  tave = tave / (dble(ntors))
         if (ntors .ne. 0)  trms = sqrt(trms/dble(ntors))
         if (prtflg.eq.1 .and. ntors.ne.0) then
            write (iout,90)  tave,trms
   90       format (/,4x,'Average Unsigned Difference :',30x,f12.4,
     &              /,4x,'Root Mean Square Deviation :',31x,f12.4)
         end if
      end if
c
c     compute the RMS between QM and TINKER gradient components
c
      gave = 0.0d0
      grms = 0.0d0
      if (fit_force) then
         allocate (derivs(3,n))
         call gradient (energy,derivs)
         if (prtflg .eq. 1) then
            write (iout,100)
  100       format (/,' Comparison of Gradient Components :',
     &              //,7x,'Atom',14x,'QM Grad',8x,'MM Grad',
     &                 10x,'Delta',/)
         end if
         do i = 1, n
            do j = 1, 3
               delta = gforce(j,i) - derivs(j,i)
               gave = gave + abs(delta)
               grms = grms + delta*delta
               if (prtflg .eq. 1) then
                  write (iout,110)  i,axis(j),gforce(j,i),
     &                              derivs(j,i),delta
  110             format (4x,i5,1x,a1,8x,f13.4,2x,f13.4,2x,f13.4)
               end if
            end do
         end do
         gave = gave / dble(3*n)
         grms = sqrt(grms/dble(3*n))
         if (prtflg .eq. 1) then
            write (iout,120)  gave,grms
  120       format (/,4x,'Average Unsigned Difference :',17x,f12.4,
     &              /,4x,'Root Mean Square Deviation :',18x,f12.4)
         end if
         deallocate (derivs)
      end if
c
c     perform dynamic allocation of some local arrays
c
      nfreq = 3 * n
      allocate (mass2(n))
      allocate (hinit(3,n))
      allocate (hstop(3,n))
      allocate (hdiag(3,n))
      allocate (hindex((nfreq*(nfreq-1))/2))
      allocate (h((nfreq*(nfreq-1))/2))
      allocate (matrix((nfreq*(nfreq+1))/2))
c
c     calculate the full Hessian matrix of second derivatives
c
      hesscut = 0.0d0
      call hessian (h,hinit,hstop,hindex,hdiag)
c
c     compute the RMS between QM and TINKER Hessian elements
c
      have = 0.0d0
      hrms = 0.0d0
      if (fit_force) then
         if (prtflg .eq. 1) then
            write (iout,130)
  130       format (/,' Comparison of Hessian Elements :',
     &              //,7x,'Atom',14x,'QM Hess',8x,'MM Hess',
     &                 10x,'Delta',/)
         end if
         do i = 1, n
            do j = 1, 3
               m1 = 3*(i-1) + j
               m = m1*(m1+1) / 2
               delta = gh(m) - hdiag(j,i)
               have = have + abs(delta)
               hrms = hrms + delta*delta
               if (prtflg .eq. 1) then
                  write (iout,140)  i,axis(j),gh(m),hdiag(j,i),delta
  140             format (4x,i5,1x,a1,8x,f13.2,2x,f13.2,2x,f13.4)
               end if
               m1 = 3*(i-1) + j
               m2 = m1
               do k = hinit(j,i), hstop(j,i)
                  m2 = m2 + 1
                  m = m1 + m2*(m2-1)/2
                  delta = gh(m) - h(k)
                  have = have + abs(delta)
                  hrms = hrms + delta* delta
               end do
            end do
         end do
         have = have / dble((9*n*n+3*n)/2)
         hrms = sqrt(hrms/dble((9*n*n+3*n)/2))
         if (prtflg .eq. 1) then
            write (iout,150)  have,hrms
  150       format (/,4x,'Average Unsigned Difference :',17x,f12.4,
     &              /,4x,'Root Mean Square Deviation :',18x,f12.4)
         end if
      end if
c
c     set atomic mass roots needed for vibrational analysis
c
      do i = 1, n
         mass2(i) = sqrt(mass(i))
      end do
c
c     store upper triangle of the mass-weighted Hessian matrix
c
      ihess = 0
      do i = 1, n
         do j = 1, 3
            ihess = ihess + 1
            matrix(ihess) = hdiag(j,i) / mass(i)
            do k = hinit(j,i), hstop(j,i)
               m = (hindex(k)+2) / 3
               ihess = ihess + 1
               matrix(ihess) = h(k) / (mass2(i)*mass2(m))
            end do
         end do
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (eigen(nfreq))
      allocate (vects(nfreq,nfreq))
c
c     diagonalize to get vibrational frequencies and normal modes
c
      call diagq (nfreq,nfreq,matrix,eigen,vects)
      factor = sqrt(convert) / (2.0d0*pi*lightspd)
      do i = 1, nfreq
         eigen(i) = factor * sign(1.0d0,eigen(i)) * sqrt(abs(eigen(i)))
      end do
c
c     compute the RMS between QM and TINKER vibrational frequencies
c
      fcut = 800.0d0
      if (fit_tors)  fcut = 200.0d0
      fave = 0.0d0
      frms = 0.0d0
      if (prtflg .eq. 1) then
         write (iout,160)
  160    format (/,' Comparison of Vibrational Frequencies :',
     &           //,6x,'Mode',15x,'QM Freq',8x,'MM Freq',10x,'Delta',/)
      end if
      k = 0
      do i = nfreq, 7, -1
         if (gfreq(i-6) .gt. fcut) then
            k = k + 1
            delta = eigen(i) - gfreq(i-6)
            fave = fave + abs(delta)
            frms = frms + delta*delta
            if (prtflg .eq. 1) then
               write (iout,170)  k,gfreq(i-6),eigen(i),delta
  170          format (4x,i5,10x,f13.2,2x,f13.2,2x,f13.4)
            end if
         end if
      end do
      fave = fave / (dble(k))
      frms = sqrt(frms/dble(k))
      if (prtflg .eq. 1) then
         write (iout,180)  fave,frms
  180    format (/,4x,'Average Unsigned Difference :',17x,f12.4,
     &           /,4x,'Root Mean Square Deviation :',18x,f12.4)
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (mass2)
      deallocate (hinit)
      deallocate (hstop)
      deallocate (hdiag)
      deallocate (hindex)
      deallocate (h)
      deallocate (matrix)
      deallocate (eigen)
      deallocate (vects)
c
c     sum weighted RMS values to get overall error function
c
      bfac = 100.0d0
      afac = 10.0d0
      tfac = 1.0d0
      gfac = 10.0d0
      hfac = 0.1d0
      ffac = 0.1d0
      valrms = bfac*brms + afac*arms + tfac*trms
     &            + gfac*grms + hfac*hrms + ffac*frms
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  function minimiz1  --  energy and gradient for minimize  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "minimiz1" is a service routine that computes the energy and
c     gradient for a low storage BFGS optimization in Cartesian
c     coordinate space
c
c
      function minimiz1 (xx,g)
      use sizes
      use atoms
      use scales
      use usage
      implicit none
      integer i,nvar
      real*8 minimiz1,e
      real*8 energy,eps
      real*8 xx(*)
      real*8 g(*)
      real*8, allocatable :: derivs(:,:)
      logical analytic
      external energy
c
c
c     use either analytical or numerical gradients
c
      analytic = .true.
      eps = 0.00001d0
c
c     translate optimization parameters to atomic coordinates
c
      nvar = 0
      do i = 1, n
         if (use(i)) then
            nvar = nvar + 1
            x(i) = xx(nvar) / scale(nvar)
            nvar = nvar + 1
            y(i) = xx(nvar) / scale(nvar)
            nvar = nvar + 1
            z(i) = xx(nvar) / scale(nvar)
         end if
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (derivs(3,n))
c
c     compute and store the energy and gradient
c
      if (analytic) then
         call gradient (e,derivs)
      else
         e = energy ()
         call numgrad (energy,derivs,eps)
      end if
      minimiz1 = e
c
c     store Cartesian gradient as optimization gradient
c
      nvar = 0
      do i = 1, n
         if (use(i)) then
            nvar = nvar + 1
            g(nvar) = derivs(1,i) / scale(nvar)
            nvar = nvar + 1
            g(nvar) = derivs(2,i) / scale(nvar)
            nvar = nvar + 1
            g(nvar) = derivs(3,i) / scale(nvar)
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (derivs)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine prmvar  --  valence terms to optimization  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "prmvar" determines the optimization values from the
c     corresponding valence potential energy parameters
c
c
      subroutine prmvar (nvar,xx)
      use sizes
      use angbnd
      use atomid
      use atoms
      use bndstr
      use iounit
      use opbend
      use strbnd
      use tors
      use units
      use urey
      use valfit
      implicit none
      integer i,k,ii,kk
      integer ia,ib,ic,id
      integer ka,kb,kc,kd
      integer ita,itb,itc,itd
      integer kta,ktb,ktc,ktd
      integer nvar,size
      real*8 xx(*)
      logical done
      character*4  pa,pb,pc,pd
      character*8  pitb,pktb
      character*12 pita,pkta
      character*16 pitt,pktt
c
c
c     zero out the total number of optimization parameters
c
      nvar = 0
c
c     print a header for the parameters used in fitting
c
      if (fit_struct) then
         write (iout,10)
   10    format (/,' Valence Parameters Used in Structure Fitting :')
      else if (fit_force) then
         write (iout,20)
   20    format (/,' Valence Parameters Used in Force Fitting :')
      end if
      write (iout,30)
   30 format (/,' Parameter',10x,'Atom Classes',10x,'Category',
     &           12x,'Value',5x,'Fixed',/)
c
c     find bond stretch force constants and target lengths
c
      do i = 1, nbond
         done = .false.
         ia = ibnd(1,i)
         ib = ibnd(2,i)
         ita = class(ia)
         itb = class(ib)
         size = 4
         call numeral (ita,pa,size)
         call numeral (itb,pb,size)
         if (ita .le. itb) then
            pitb = pa//pb
         else
            pitb = pb//pa
         end if
         do k = 1, i-1
            ka = ibnd(1,k)
            kb = ibnd(2,k)
            kta = class(ka)
            ktb = class(kb)
            size = 4
            call numeral (kta,pa,size)
            call numeral (ktb,pb,size)
            if (kta .le. ktb) then
               pktb = pa//pb
            else
               pktb = pb//pa
            end if
            if (pktb .eq. pitb)  done = .true.
         end do
         if (.not. done) then
            if (fit_bond .and. bk(i).ne.0.0d0) then
               nvar = nvar + 1
               xx(nvar) = bk(i)
               write (iout,40)  nvar,ita,itb,'Bond Force',bk(i)
   40          format (i6,5x,2i6,19x,a10,3x,f12.4)
               nvar = nvar + 1
               xx(nvar) = bl(i)
               xx(nvar) = 100.0d0 * xx(nvar)
               write (iout,50)  nvar,ita,itb,'Bond Length',bl(i)
   50          format (i6,5x,2i6,19x,a11,2x,f12.4)
            else
               write (iout,60)  ita,itb,'Bond Force',bk(i)
   60          format (4x,'--',5x,2i6,19x,a10,3x,f12.4,7x,'X')
               write (iout,70)  ita,itb,'Bond Length',bl(i)
   70          format (4x,'--',5x,2i6,19x,a11,2x,f12.4,7x,'X')
            end if
         end if
      end do
c
c     find angle bend force constants and target angles
c
      do i = 1, nangle
         done = .false.
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
            pita = pa//pb//pc
         else
            pita = pc//pb//pa
         end if
         do k = 1, i-1
            ka = iang(1,k)
            kb = iang(2,k)
            kc = iang(3,k)
            kta = class(ka)
            ktb = class(kb)
            ktc = class(kc)
            size = 4
            call numeral (kta,pa,size)
            call numeral (ktb,pb,size)
            call numeral (ktc,pc,size)
            if (kta .le. ktc) then
               pkta = pa//pb//pc
            else
               pkta = pc//pb//pa
            end if
            if (pkta .eq. pita)  done = .true.
         end do
         if (.not. done) then
            if (fit_angle .and. ak(i).ne.0.0d0) then
               nvar = nvar + 1
               xx(nvar) = ak(i)
               write (iout,80)  nvar,ita,itb,itc,'Angle Force',ak(i)
   80          format (i6,5x,3i6,13x,a11,2x,f12.4)
               nvar = nvar + 1
               xx(nvar) = anat(i)
               write (iout,90)  nvar,ita,itb,itc,'Angle Value',anat(i)
   90          format (i6,5x,3i6,13x,a11,2x,f12.4)
            else
               write (iout,100)  ita,itb,itc,'Angle Force',ak(i)
  100          format (4x,'--',5x,3i6,13x,a11,2x,f12.4,7x,'X')
               write (iout,110)  ita,itb,itc,'Angle Value',anat(i)
  110          format (4x,'--',5x,3i6,13x,a11,2x,f12.4,7x,'X')
            end if
         end if
      end do
c
c     find stretch-bend force constant parameter values
c
      do i = 1, nstrbnd
         done = .false.
         ii = isb(1,i)
         ia = iang(1,ii)
         ib = iang(2,ii)
         ic = iang(3,ii)
         ita = class(ia)
         itb = class(ib)
         itc = class(ic)
         size = 4
         call numeral (ita,pa,size)
         call numeral (itb,pb,size)
         call numeral (itc,pc,size)
         if (ita .le. itc) then
            pita = pa//pb//pc
         else
            pita = pc//pb//pa
         end if
         do k = 1, i-1
            kk = isb(1,k)
            ka = iang(1,kk)
            kb = iang(2,kk)
            kc = iang(3,kk)
            kta = class(ka)
            ktb = class(kb)
            ktc = class(kc)
            size = 4
            call numeral (kta,pa,size)
            call numeral (ktb,pb,size)
            call numeral (ktc,pc,size)
            if (kta .le. ktc) then
               pkta = pa//pb//pc
            else
               pkta = pc//pb//pa
            end if
            if (pkta .eq. pita)  done = .true.
         end do
         if (.not. done) then
            if (fit_strbnd .and. sbk(1,i).ne.0.0d0
     &             .and. sbk(2,i).ne.0.0d0) then
               nvar = nvar + 1
               xx(nvar) = sbk(1,i)
               write (iout,120)  nvar,ita,itb,itc,'StrBnd-1',sbk(1,i)
  120          format (i6,5x,3i6,13x,a8,5x,f12.4)
               nvar = nvar + 1
               xx(nvar) = sbk(2,i)
               write (iout,130)  nvar,ita,itb,itc,'StrBnd-2',sbk(2,i)
  130          format (i6,5x,3i6,13x,a8,5x,f12.4)
            else
               write (iout,140)  ita,itb,itc,'StrBnd-1',sbk(1,i)
  140          format (4x,'--',5x,3i6,13x,a8,5x,f12.4,7x,'X')
               write (iout,150)  ita,itb,itc,'StrBnd-2',sbk(2,i)
  150          format (4x,'--',5x,3i6,13x,a8,5x,f12.4,7x,'X')
            end if
         end if
      end do
c
c     find Urey-Bradley force constant parameter values
c
      do i = 1, nurey
         done = .false.
         ia = iury(1,i)
         ib = iury(2,i)
         ic = iury(3,i)
         ita = class(ia)
         itb = class(ib)
         itc = class(ic)
         size = 4
         call numeral (ita,pa,size)
         call numeral (itb,pb,size)
         call numeral (itc,pc,size)
         if (ita .le. itc) then
            pita = pa//pb//pc
         else
            pita = pc//pb//pa
         end if
         do k = 1, i-1
            ka = iury(1,k)
            kb = iury(2,k)
            kc = iury(3,k)
            kta = class(ka)
            ktb = class(kb)
            ktc = class(kc)
            size = 4
            call numeral (kta,pa,size)
            call numeral (ktb,pb,size)
            call numeral (ktc,pc,size)
            if (kta .le. ktc) then
               pkta = pa//pb//pc
            else
               pkta = pc//pb//pa
            end if
            if (pkta .eq. pita)  done = .true.
         end do
         if (.not. done) then
            if (fit_urey .and. uk(i).ne.0.0d0) then
               nvar = nvar + 1
               xx(nvar) = uk(i)
               write (iout,160)  nvar,ita,itb,itc,'Urey Force',uk(i)
  160          format (i6,5x,3i6,13x,a10,3x,f12.4)
               nvar = nvar + 1
               xx(nvar) = ul(i)
               write (iout,170)  nvar,ita,itb,itc,'Urey Dist',ul(i)
  170          format (i6,5x,3i6,13x,a9,4x,f12.4)
            else
               write (iout,180)  ita,itb,itc,'Urey Force',uk(i)
  180          format (4x,'--',5x,3i6,13x,a10,3x,f12.4,7x,'X')
               write (iout,190)  ita,itb,itc,'Urey Dist',ul(i)
  190          format (4x,'--',5x,3i6,13x,a9,4x,f12.4,7x,'X')
            end if
         end if
      end do
c
c     find out-of-plane bend force constant parameter values
c
      do i = 1, nopbend
         done = .false.
         ii = iopb(i)
         ia = iang(1,ii)
         ib = iang(2,ii)
         ic = iang(3,ii)
         id = iang(4,ii)
         ita = class(ia)
         itb = class(ib)
         itc = class(ic)
         itd = class(id)
         size = 4
         call numeral (ita,pa,size)
         call numeral (itb,pb,size)
         call numeral (itc,pc,size)
         call numeral (itd,pd,size)
         if (ita .le. itc) then
            pitt = pd//pb//pa//pc
         else
            pitt = pd//pb//pc//pa
         end if
         do k = 1, i-1
            kk = iopb(k)
            ka = iang(1,kk)
            kb = iang(2,kk)
            kc = iang(3,kk)
            kd = iang(4,kk)
            kta = class(ka)
            ktb = class(kb)
            ktc = class(kc)
            ktd = class(kd)
            size = 4
            call numeral (kta,pa,size)
            call numeral (ktb,pb,size)
            call numeral (ktc,pc,size)
            call numeral (ktd,pd,size)
            if (kta .le. ktc) then
               pktt = pd//pb//pa//pc
            else
               pktt = pd//pb//pc//pa
            end if
            if (pktt .eq. pitt)  done = .true.
         end do
         if (.not. done) then
            if (fit_opbend .and. opbk(i).ne.0.0d0) then
               nvar = nvar + 1
               xx(nvar) = opbk(i)
               write (iout,200)  nvar,itd,itb,min(ita,itc),
     &                           max(ita,itc),'O-P-Bend',opbk(i)
  200          format (i6,5x,4i6,7x,a8,5x,f12.4)
            else
               write (iout,210)  itd,itb,min(ita,itc),max(ita,itc),
     &                           'O-P-Bend',opbk(i)
  210          format (4x,'--',5x,4i6,7x,a8,5x,f12.4,7x,'X')
            end if
         end if
      end do
c
c     find torsional angle amplitude parameter values
c
      do i = 1, ntors
         done = .false.
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
         if (itb.lt.itc .or. (itb.eq.itc.and.ita.le.itd)) then
            pitt = pa//pb//pc//pd
         else
            pitt = pd//pc//pb//pa
         end if
         do k = 1, i-1
            ka = itors(1,k)
            kb = itors(2,k)
            kc = itors(3,k)
            kd = itors(4,k)
            kta = class(ka)
            ktb = class(kb)
            ktc = class(kc)
            ktd = class(kd)
            size = 4
            call numeral (kta,pa,size)
            call numeral (ktb,pb,size)
            call numeral (ktc,pc,size)
            call numeral (ktd,pd,size)
            if (ktb.lt.ktc .or. (ktb.eq.ktc.and.kta.le.ktd)) then
               pktt = pa//pb//pc//pd
            else
               pktt = pd//pc//pb//pa
            end if
            if (pktt .eq. pitt)  done = .true.
         end do
         if (.not. done) then
            if (fit_tors .and. tors1(1,i).ne.0.0d0) then
               nvar = nvar + 1
               xx(nvar) = tors1(1,i)
               write (iout,220)  nvar,ita,itb,itc,itd,
     &                           'Torsion-1',tors1(1,i)
  220          format (i6,5x,4i6,7x,a9,4x,f12.4)
            else
               write (iout,230)  ita,itb,itc,itd,'Torsion-1',tors1(1,i)
  230          format (4x,'--',5x,4i6,7x,a9,4x,f12.4,7x,'X')
            end if
            if (fit_tors .and. tors2(1,i).ne.0.0d0) then
               nvar = nvar + 1
               xx(nvar) = tors2(1,i)
               write (iout,240)  nvar,ita,itb,itc,itd,
     &                           'Torsion-2',tors2(1,i)
  240          format (i6,5x,4i6,7x,a9,4x,f12.4)
            else
               write (iout,250)  ita,itb,itc,itd,'Torsion-2',tors2(1,i)
  250          format (4x,'--',5x,4i6,7x,a9,4x,f12.4,7x,'X')
            end if
            if (fit_tors .and. tors3(1,i).ne.0.0d0) then
               nvar = nvar + 1
               xx(nvar) = tors3(1,i)
               write (iout,260)  nvar,ita,itb,itc,itd,
     &                           'Torsion-3',tors3(1,i)
  260          format (i6,5x,4i6,7x,a9,4x,f12.4)
            else
               write (iout,270)  ita,itb,itc,itd,'Torsion-3',tors3(1,i)
  270          format (4x,'--',5x,4i6,7x,a9,4x,f12.4,7x,'X')
            end if
            if (fit_tors .and. tors4(1,i).ne.0.0d0) then
               nvar = nvar + 1
               xx(nvar) = tors4(1,i)
               write (iout,280)  nvar,ita,itb,itc,itd,
     &                           'Torsion-4',tors4(1,i)
  280          format (i6,5x,4i6,7x,a9,4x,f12.4)
            else if (tors4(1,i) .ne. 0.0d0) then
               write (iout,290)  ita,itb,itc,itd,'Torsion-4',tors4(1,i)
  290          format (4x,'--',5x,4i6,7x,a9,4x,f12.4,7x,'X')
            end if
            if (fit_tors .and. tors5(1,i).ne.0.0d0) then
               nvar = nvar + 1
               xx(nvar) = tors5(1,i)
               write (iout,300)  nvar,ita,itb,itc,itd,
     &                           'Torsion-5',tors5(1,i)
  300          format (i6,5x,4i6,7x,a9,4x,f12.4)
            else if (tors5(1,i) .ne. 0.0d0) then
               write (iout,310)  ita,itb,itc,itd,'Torsion-5',tors5(1,i)
  310          format (4x,'--',5x,4i6,7x,a9,4x,f12.4,7x,'X')
            end if
            if (fit_tors .and. tors6(1,i).ne.0.0d0) then
               nvar = nvar + 1
               xx(nvar) = tors6(1,i)
               write (iout,320)  nvar,ita,itb,itc,itd,
     &                           'Torsion-6',tors6(1,i)
  320          format (i6,5x,4i6,7x,a9,4x,f12.4)
            else if (tors6(1,i) .ne. 0.0d0) then
               write (iout,330)  ita,itb,itc,itd,'Torsion-6',tors6(1,i)
  330          format (4x,'--',5x,4i6,7x,a9,4x,f12.4,7x,'X')
            end if
         end if
      end do
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine varprm  --  optimization to valence terms  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "varprm" copies the current optimization values into the
c     corresponding valence potential energy parameters
c
c
      subroutine varprm (nvar,xx,ivar,eps)
      use sizes
      use angbnd
      use atoms
      use atomid
      use bndstr
      use opbend
      use potent
      use strbnd
      use tors
      use urey
      use valfit
      implicit none
      integer i,k,ii,kk
      integer nvar,ivar,size
      integer ia,ib,ic,id
      integer ka,kb,kc,kd
      integer ita,itb,itc,itd
      integer kta,ktb,ktc,ktd
      real*8 eps
      real*8 xx(*)
      logical done
      character*4 pa,pb,pc,pd
      character*8 pitb,pktb
      character*12 pita,pkta
      character*16 pitt,pktt
c
c
c     zero out the total number of optimization parameters
c
      nvar = 0
c
c     translate optimization values to bond stretch parameters
c
      do i = 1, nbond
         done = .false.
         ia = ibnd(1,i)
         ib = ibnd(2,i)
         ita = class(ia)
         itb = class(ib)
         size = 4
         call numeral (ita,pa,size)
         call numeral (itb,pb,size)
         if (ita .le. itb) then
            pitb = pa//pb
         else
            pitb = pb//pa
         end if
         do k = 1, i-1
            ka = ibnd(1,k)
            kb = ibnd(2,k)
            kta = class(ka)
            ktb = class(kb)
            size = 4
            call numeral (kta,pa,size)
            call numeral (ktb,pb,size)
            if (kta .le. ktb) then
               pktb = pa//pb
            else
               pktb = pb//pa
            end if
            if (pktb .eq. pitb)  done = .true.
         end do
         if (.not. done) then
            if (fit_bond .and. bk(i).ne.0.0d0) then
               nvar = nvar + 1
               bk(i) = xx(nvar)
               if (ivar .eq. nvar)  bk(i) = bk(i) + eps
               nvar = nvar + 1
               bl(i) = xx(nvar)
               if (ivar .eq. nvar)  bl(i) = bl(i) + eps
               bl(i) = 0.01d0 * bl(i)
               do k = i+1, nbond
                  ka = ibnd(1,k)
                  kb = ibnd(2,k)
                  kta = class(ka)
                  ktb = class(kb)
                  size = 4
                  call numeral (kta,pa,size)
                  call numeral (ktb,pb,size)
                  if (kta .le. ktb) then
                     pktb = pa//pb
                  else
                     pktb = pb//pa
                  end if
                  if (pktb .eq. pitb) then
                     bk(k) = bk(i)
                     bl(k) = bl(i)
                  end if
               end do
            end if
         end if
      end do
c
c     translate optimization values to angle bend parameters
c
      do i = 1, nangle
         done = .false.
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
            pita = pa//pb//pc
         else
            pita = pc//pb//pa
         end if
         do k = 1, i-1
            ka = iang(1,k)
            kb = iang(2,k)
            kc = iang(3,k)
            kta = class(ka)
            ktb = class(kb)
            ktc = class(kc)
            size = 4
            call numeral (kta,pa,size)
            call numeral (ktb,pb,size)
            call numeral (ktc,pc,size)
            if (kta .le. ktc) then
               pkta = pa//pb//pc
            else
               pkta = pc//pb//pa
            end if
            if (pkta .eq. pita)  done = .true.
         end do
         if (.not. done) then
            if (fit_angle .and. ak(i).ne.0.0d0) then
               nvar = nvar + 1
               ak(i) = xx(nvar)
               if (ivar .eq. nvar)  ak(i) = ak(i) + eps
               nvar = nvar + 1
               anat(i) = xx(nvar)
               if (ivar .eq. nvar)  anat(i) = anat(i) + eps
               do k = i+1, nangle
                  ka = iang(1,k)
                  kb = iang(2,k)
                  kc = iang(3,k)
                  kta = class(ka)
                  ktb = class(kb)
                  ktc = class(kc)
                  size = 4
                  call numeral (kta,pa,size)
                  call numeral (ktb,pb,size)
                  call numeral (ktc,pc,size)
                  if (kta .le. ktc) then
                     pkta = pa//pb//pc
                  else
                     pkta = pc//pb//pa
                  end if
                  if (pkta .eq. pita) then
                     ak(k) = ak(i)
                     anat(k) = anat(i)
                  end if
               end do
            end if
         end if
      end do
c
c     translate optimization values to stretch-bend parameters
c
      do i = 1, nstrbnd
         done = .false.
         ii = isb(1,i)
         ia = iang(1,ii)
         ib = iang(2,ii)
         ic = iang(3,ii)
         ita = class(ia)
         itb = class(ib)
         itc = class(ic)
         size = 4
         call numeral (ita,pa,size)
         call numeral (itb,pb,size)
         call numeral (itc,pc,size)
         if (ita .le. itc) then
            pita = pa//pb//pc
         else
            pita = pc//pb//pa
         end if
         do k = 1, i-1
            kk = isb(1,k)
            ka = iang(1,kk)
            kb = iang(2,kk)
            kc = iang(3,kk)
            kta = class(ka)
            ktb = class(kb)
            ktc = class(kc)
            size = 4
            call numeral (kta,pa,size)
            call numeral (ktb,pb,size)
            call numeral (ktc,pc,size)
            if (kta .le. ktc) then
               pkta = pa//pb//pc
            else
               pkta = pc//pb//pa
            end if
            if (pkta .eq. pita)  done = .true.
         end do
         if (.not.done .and. fit_strbnd) then
            if (sbk(1,i) .ne. 0.0d0) then
               nvar = nvar + 1
               sbk(1,i) = xx(nvar)
               if (ivar .eq. nvar)  sbk(1,i) = sbk(1,i) + eps
            end if
            if (sbk(2,i) .ne. 0.0d0) then
               nvar = nvar + 1
               sbk(2,i) = xx(nvar)
               if (ivar .eq. nvar)  sbk(2,i) = sbk(2,i) + eps
            end if
            do k = i+1, nstrbnd
               kk = isb(1,k)
               ka = iang(1,kk)
               kb = iang(2,kk)
               kc = iang(3,kk)
               kta = class(ka)
               ktb = class(kb)
               ktc = class(kc)
               size = 4
               call numeral (kta,pa,size)
               call numeral (ktb,pb,size)
               call numeral (ktc,pc,size)
               if (kta .le. ktc) then
                  pkta = pa//pb//pc
               else
                  pkta = pc//pb//pa
               end if
               if (pkta .eq. pita) then
                  if (kta .eq. ita) then
                     sbk(1,k) = sbk(1,i)
                     sbk(2,k) = sbk(2,i)
                  else
                     sbk(2,k) = sbk(1,i)
                     sbk(1,k) = sbk(2,i)
                  end if
               end if
            end do
         end if
      end do
c
c     translate optimization values to Urey-Bradley parameters
c
      do i = 1, nurey
         done = .false.
         ia = iury(1,i)
         ib = iury(2,i)
         ic = iury(3,i)
         ita = class(ia)
         itb = class(ib)
         itc = class(ic)
         size = 4
         call numeral (ita,pa,size)
         call numeral (itb,pb,size)
         call numeral (itc,pc,size)
         if (ita .le. itc) then
            pita = pa//pb//pc
         else
            pita = pc//pb//pa
         end if
         do k = 1, i-1
            ka = iury(1,k)
            kb = iury(2,k)
            kc = iury(3,k)
            kta = class(ka)
            ktb = class(kb)
            ktc = class(kc)
            size = 4
            call numeral (kta,pa,size)
            call numeral (ktb,pb,size)
            call numeral (ktc,pc,size)
            if (kta .le. ktc) then
               pkta = pa//pb//pc
            else
               pkta = pc//pb//pa
            end if
            if (pkta .eq. pita)  done = .true.
         end do
         if (.not. done) then
            if (fit_urey .and. uk(i).ne.0.0d0) then
               nvar = nvar + 1
               uk(i) = xx(nvar)
               if (ivar .eq. nvar)  uk(i) = uk(i) + eps
               nvar = nvar + 1
               ul(i) = xx(nvar)
               if (ivar .eq. nvar)  ul(i) = ul(i) + eps
               do k = i+1, nurey
                  ka = iury(1,k)
                  kb = iury(2,k)
                  kc = iury(3,k)
                  kta = class(ka)
                  ktb = class(kb)
                  ktc = class(kc)
                  size = 4
                  call numeral (kta,pa,size)
                  call numeral (ktb,pb,size)
                  call numeral (ktc,pc,size)
                  if (kta .le. ktc) then
                     pkta = pa//pb//pc
                  else
                     pkta = pc//pb//pa
                  end if
                  if (pkta .eq. pita) then
                     uk(k) = uk(i)
                     ul(k) = ul(i)
                  end if
               end do
            end if
         end if
      end do
c
c     translate optimization values to out-of-plane bend parameters
c
      do i = 1, nopbend
         done = .false.
         ii = iopb(i)
         ia = iang(1,ii)
         ib = iang(2,ii)
         ic = iang(3,ii)
         id = iang(4,ii)
         ita = class(ia)
         itb = class(ib)
         itc = class(ic)
         itd = class(id)
         size = 4
         call numeral (ita,pa,size)
         call numeral (itb,pb,size)
         call numeral (itc,pc,size)
         call numeral (itd,pd,size)
         if (ita .le. itc) then
            pitt = pd//pb//pa//pc
         else
            pitt = pd//pb//pc//pa
         end if
         do k = 1, i-1
            kk = iopb(k)
            ka = iang(1,kk)
            kb = iang(2,kk)
            kc = iang(3,kk)
            kd = iang(4,kk)
            kta = class(ka)
            ktb = class(kb)
            ktc = class(kc)
            ktd = class(kd)
            size = 4
            call numeral (kta,pa,size)
            call numeral (ktb,pb,size)
            call numeral (ktc,pc,size)
            call numeral (ktd,pd,size)
            if (kta .le. ktc) then
               pktt = pd//pb//pa//pc
            else
               pktt = pd//pb//pc//pa
            end if
            if (pktt .eq. pitt)  done = .true.
         end do
         if (.not. done) then
            if (fit_opbend .and. opbk(i).ne.0.0d0) then
               nvar = nvar + 1
               opbk(i) = xx(nvar)
               if (ivar .eq. nvar)  opbk(i) = opbk(i) + eps
               do k = i+1, nopbend
                  kk = iopb(k)
                  ka = iang(1,kk)
                  kb = iang(2,kk)
                  kc = iang(3,kk)
                  kd = iang(4,kk)
                  kta = class(ka)
                  ktb = class(kb)
                  ktc = class(kc)
                  ktd = class(kd)
                  size = 4
                  call numeral (kta,pa,size)
                  call numeral (ktb,pb,size)
                  call numeral (ktc,pc,size)
                  call numeral (ktd,pd,size)
                  if (kta .le. ktc) then
                     pktt = pd//pb//pa//pc
                  else
                     pktt = pd//pb//pc//pa
                  end if
                  if (pktt.eq.pitt)  opbk(k) = opbk(i)
               end do
            end if
         end if
      end do
c
c     translate optimization values to torsional parameters
c
      do i = 1, ntors
         done = .false.
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
         if (itb.lt.itc .or. (itb.eq.itc.and.ita.le.itd)) then
            pitt = pa//pb//pc//pd
         else
            pitt = pd//pc//pb//pa
         end if
         do k = 1, i-1
            ka = itors(1,k)
            kb = itors(2,k)
            kc = itors(3,k)
            kd = itors(4,k)
            kta = class(ka)
            ktb = class(kb)
            ktc = class(kc)
            ktd = class(kd)
            size = 4
            call numeral (kta,pa,size)
            call numeral (ktb,pb,size)
            call numeral (ktc,pc,size)
            call numeral (ktd,pd,size)
            if (ktb.lt.ktc .or. (ktb.eq.ktc.and.kta.le.ktd)) then
               pktt = pa//pb//pc//pd
            else
               pktt = pd//pc//pb//pa
            end if
            if (pktt .eq. pitt)  done = .true.
         end do
         if (.not.done .and. fit_tors) then
            if (tors1(1,i) .ne. 0.0d0) then
               nvar = nvar + 1
               tors1(1,i) = xx(nvar)
               if (ivar .eq. nvar)  tors1(1,i) = tors1(1,i) + eps
            end if
            if (tors2(1,i) .ne. 0.0d0) then
               nvar = nvar + 1
               tors2(1,i) = xx(nvar)
               if (ivar .eq. nvar)  tors2(1,i) = tors2(1,i) + eps
            end if
            if (tors3(1,i) .ne. 0.0d0) then
               nvar = nvar + 1
               tors3(1,i) = xx(nvar)
               if (ivar .eq. nvar)  tors3(1,i) = tors3(1,i) + eps
            end if
            if (tors4(1,i) .ne. 0.0d0) then
               nvar = nvar + 1
               tors4(1,i) = xx(nvar)
               if (ivar .eq. nvar)  tors4(1,i) = tors4(1,i) + eps
            end if
            if (tors5(1,i) .ne. 0.0d0) then
               nvar = nvar + 1
               tors5(1,i) = xx(nvar)
               if (ivar .eq. nvar)  tors5(1,i) = tors5(1,i) + eps
            end if
            if (tors6(1,i) .ne. 0.0d0) then
               nvar = nvar + 1
               tors6(1,i) = xx(nvar)
               if (ivar .eq. nvar)  tors6(1,i) = tors6(1,i) + eps
            end if
            do k = i+1, ntors
               ka = itors(1,k)
               kb = itors(2,k)
               kc = itors(3,k)
               kd = itors(4,k)
               kta = class(ka)
               ktb = class(kb)
               ktc = class(kc)
               ktd = class(kd)
               size = 4
               call numeral (kta,pa,size)
               call numeral (ktb,pb,size)
               call numeral (ktc,pc,size)
               call numeral (ktd,pd,size)
               if (ktb.lt.ktc .or. (ktb.eq.ktc.and.kta.le.ktd)) then
                  pktt = pa//pb//pc//pd
               else
                  pktt = pd//pc//pb//pa
               end if
               if (pktt .eq. pitt) then
                  tors1(1,k) = tors1(1,i)
                  tors2(1,k) = tors2(1,i)
                  tors3(1,k) = tors3(1,i)
                  tors4(1,k) = tors4(1,i)
                  tors5(1,k) = tors5(1,i)
                  tors6(1,k) = tors6(1,i)
               end if
            end do
         end if
      end do
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  function valfit1  --  valence fit error and gradient  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "valfit1" is a service routine that computes the RMS error
c     and gradient for valence parameters fit to QM results
c
c
      function valfit1 (xx,g)
      use sizes
      use atoms
      use potent
      use valfit
      implicit none
      integer i,k
      integer nvar
      real*8 e,e0
      real*8 delta
      real*8 valrms
      real*8 valfit1
      real*8 xx(*)
      real*8 g(*)
      real*8, allocatable :: eps(:)
c
c
c     copy optimization values to valence parameters
c
      call varprm (nvar,xx,0,0.0d0)
c
c     perform dynamic allocation of some local arrays
c
      allocate (eps(nvar))
c
c     set the numerical gradient step size for each parameter
c
      delta = 0.0000001d0
      do i = 1, nvar
         eps(i) = delta * xx(i)
      end do
c
c     get the RMS of frequencies
c
      valfit1 = valrms(0)
c
c     compute numerical gradient for valence parameters
c
      k = nvar
      do i = 1, k
         call varprm (nvar,xx,i,-0.5d0*eps(i))
         e0 = valrms(0)
         call varprm (nvar,xx,i,0.5d0*eps(i))
         e  = valrms(0)
         g(i) = (e-e0) / eps(i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (eps)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##   subroutine prtval  --  print final valence parameter fit  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "prtval" writes the final valence parameter results to the
c     standard output and appends the values to a key file
c
c
      subroutine prtval
      use sizes
      use angbnd
      use atomid
      use atoms
      use bndstr
      use files
      use iounit
      use keys
      use opbend
      use strbnd
      use tors
      use units
      use urey
      use valfit
      implicit none
      integer i,k,ii,kk
      integer ia,ib,ic,id
      integer ka,kb,kc,kd
      integer ita,itb,itc,itd
      integer kta,ktb,ktc,ktd
      integer ikey,size
      integer freeunit
      integer trimtext
      logical done
      character*4 pa,pb,pc,pd
      character*8 pitb,pktb
      character*12 pita,pkta
      character*16 pitt,pktt
      character*240 keyfile
      character*240 record
c
c
c     output some definitions and parameters to a keyfile
c
      ikey = freeunit ()
      keyfile = filename(1:leng)//'.key'
      call version (keyfile,'new')
      open (unit=ikey,file=keyfile,status='new')
c
c     copy the contents of any previously existing keyfile
c
      do i = 1, nkey
         record = keyline(i)
         size = trimtext (record)
         write (ikey,10)  record(1:size)
   10    format (a)
      end do
c
c     print a header for the fitted valence parameters
c
      if (fit_bond .or. fit_angle .or. fit_tors
     &       .or. fit_strbnd .or. fit_opbend) then
         write (ikey,20)
   20    format (/,'#',/,'# Results of Valence Parameter Fitting',
     &              /,'#',/)
      end if
c
c     output any fitted bond stretch parameter values
c
      do i = 1, nbond
         done = .false.
         ia = ibnd(1,i)
         ib = ibnd(2,i)
         ita = class(ia)
         itb = class(ib)
         size = 4
         call numeral (ita,pa,size)
         call numeral (itb,pb,size)
         if (ita .le. itb) then
            pitb = pa//pb
         else
            pitb = pb//pa
         end if
         do k = 1, i-1
            ka = ibnd(1,k)
            kb = ibnd(2,k)
            kta = class(ka)
            ktb = class(kb)
            size = 4
            call numeral (kta,pa,size)
            call numeral (ktb,pb,size)
            if (kta .le. ktb) then
               pktb = pa//pb
            else
               pktb = pb//pa
            end if
            if (pktb .eq. pitb)  done = .true.
         end do
         if (.not. done) then
            if (fit_bond) then
               write (ikey,30)  ita,itb,bk(i),bl(i)
   30          format ('bond',6x,2i5,5x,f11.2,f11.4)
            end if
         end if
      end do
c
c     output any fitted angle bend parameter values
c
      do i = 1, nangle
         done = .false.
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
            pita = pa//pb//pc
         else
            pita = pc//pb//pa
         end if
         do k = 1, i-1
            ka = iang(1,k)
            kb = iang(2,k)
            kc = iang(3,k)
            kta = class(ka)
            ktb = class(kb)
            ktc = class(kc)
            size = 4
            call numeral (kta,pa,size)
            call numeral (ktb,pb,size)
            call numeral (ktc,pc,size)
            if (kta .le. ktc) then
               pkta = pa//pb//pc
            else
               pkta = pc//pb//pa
            end if
            if (pkta .eq. pita)  done = .true.
         end do
         if (.not. done) then
            if (fit_angle) then
               write (ikey,40)  ita,itb,itc,ak(i),anat(i)
   40          format ('angle',5x,3i5,f11.2,f11.2)
            end if
         end if
      end do
c
c     output any fitted stretch-bend parameter values
c
      do i = 1, nstrbnd
         done = .false.
         ii = isb(1,i)
         ia = iang(1,ii)
         ib = iang(2,ii)
         ic = iang(3,ii)
         ita = class(ia)
         itb = class(ib)
         itc = class(ic)
         size = 4
         call numeral (ita,pa,size)
         call numeral (itb,pb,size)
         call numeral (itc,pc,size)
         if (ita .le. itc) then
            pita = pa//pb//pc
         else
            pita = pc//pb//pa
         end if
         do k = 1, i-1
            kk = isb(1,k)
            ka = iang(1,kk)
            kb = iang(2,kk)
            kc = iang(3,kk)
            kta = class(ka)
            ktb = class(kb)
            ktc = class(kc)
            size = 4
            call numeral (kta,pa,size)
            call numeral (ktb,pb,size)
            call numeral (ktc,pc,size)
            if (kta .le. ktc) then
               pkta = pa//pb//pc
            else
               pkta = pc//pb//pa
            end if
            if (pkta .eq. pita)  done = .true.
         end do
         if (.not. done) then
            if (fit_strbnd) then
               write (ikey,50)  ita,itb,itc,sbk(1,i),sbk(2,i)
   50          format ('strbnd',4x,3i5,2f11.3)
            end if
         end if
      end do
c
c     output any fitted Urey-Bradley parameter values
c
      do i = 1, nurey
         done = .false.
         ia = iury(1,i)
         ib = iury(2,i)
         ic = iury(3,i)
         ita = class(ia)
         itb = class(ib)
         itc = class(ic)
         size = 4
         call numeral (ita,pa,size)
         call numeral (itb,pb,size)
         call numeral (itc,pc,size)
         if (ita .le. itc) then
            pita = pa//pb//pc
         else
            pita = pc//pb//pa
         end if
         do k = 1, i-1
            ka = iury(1,k)
            kb = iury(2,k)
            kc = iury(3,k)
            kta = class(ka)
            ktb = class(kb)
            ktc = class(kc)
            size = 4
            call numeral (kta,pa,size)
            call numeral (ktb,pb,size)
            call numeral (ktc,pc,size)
            if (kta .le. ktc) then
               pkta = pa//pb//pc
            else
               pkta = pc//pb//pa
            end if
            if (pkta .eq. pita)  done = .true.
         end do
         if (.not. done) then
            if (fit_urey) then
               write (ikey,60)  ita,itb,itc,uk(i),ul(i)
   60          format ('ureybrad',2x,3i5,f11.3,f11.4)
            end if
         end if
      end do
c
c     output any fitted out-of-plane bend parameter values
c
      do i = 1, nopbend
         done = .false.
         ii = iopb(i)
         ia = iang(1,ii)
         ib = iang(2,ii)
         ic = iang(3,ii)
         id = iang(4,ii)
         ita = class(ia)
         itb = class(ib)
         itc = class(ic)
         itd = class(id)
         size = 4
         call numeral (ita,pa,size)
         call numeral (itb,pb,size)
         call numeral (itc,pc,size)
         call numeral (itd,pd,size)
         if (ita .le. itc) then
            pitt = pd//pb//pa//pc
         else
            pitt = pd//pb//pc//pa
         end if
         do k = 1, i-1
            kk = iopb(k)
            ka = iang(1,kk)
            kb = iang(2,kk)
            kc = iang(3,kk)
            kd = iang(4,kk)
            kta = class(ka)
            ktb = class(kb)
            ktc = class(kc)
            ktd = class(kd)
            size = 4
            call numeral (kta,pa,size)
            call numeral (ktb,pb,size)
            call numeral (ktc,pc,size)
            call numeral (ktd,pd,size)
            if (kta .le. ktc) then
               pktt = pd//pb//pa//pc
            else
               pktt = pd//pb//pc//pa
            end if
            if (pktt .eq. pitt)  done = .true.
         end do
         if (.not. done) then
            if (fit_opbend) then
               write (ikey,70)  itd,itb,min(ita,itc),
     &                          max(ita,itc),opbk(i)
   70          format ('opbend',4x,4i5,6x,f11.2)
            end if
         end if
      end do
c
c     output any fitted torsional parameter values
c
      do i = 1, ntors
         done = .false.
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
         if (itb.lt.itc .or. (itb.eq.itc.and.ita.le.itd)) then
            pitt = pa//pb//pc//pd
         else
            pitt = pd//pc//pb//pa
         end if
         do k = 1, i-1
            ka = itors(1,k)
            kb = itors(2,k)
            kc = itors(3,k)
            kd = itors(4,k)
            kta = class(ka)
            ktb = class(kb)
            ktc = class(kc)
            ktd = class(kd)
            size = 4
            call numeral (kta,pa,size)
            call numeral (ktb,pb,size)
            call numeral (ktc,pc,size)
            call numeral (ktd,pd,size)
            if (ktb.lt.ktc .or. (ktb.eq.ktc.and.kta.le.ktd)) then
               pktt = pa//pb//pc//pd
            else
               pktt = pd//pc//pb//pa
            end if
            if (pktt .eq. pitt)  done = .true.
         end do
         if (.not. done) then
            if (fit_tors) then
               write (ikey,80)  ita,itb,itc,itd,tors1(1,i),
     &                          tors2(1,i),tors3(1,i)
   80          format ('torsion',3x,4i5,3x,f8.3,' 0.0 1',f8.3,
     &                    ' 180.0 2',f8.3,' 0.0 3')
            end if
         end if
      end do
      return
      end
