c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine readprm  --  input of force field parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "readprm" processes the potential energy parameter file
c     in order to define the default force field parameters
c
c
      subroutine readprm
      use sizes
      use ctran
      use fields
      use iounit
      use kanang
      use kangs
      use kantor
      use katoms
      use kbonds
      use kchrge
      use kdipol
      use khbond
      use kiprop
      use kitors
      use kmulti
      use kopbnd
      use kopdst
      use korbs
      use kpitor
      use kpolr
      use kstbnd
      use ksttor
      use ktorsn
      use ktrtor
      use kurybr
      use kvdws
      use kvdwpr
      use merck
      use params
      implicit none
      integer i,j,iprm
      integer ia,ib,ic,id
      integer ie,if,ig,ih
      integer ii,imin
      integer size,next
      integer length,trimtext
      integer nb,nb5,nb4,nb3,nel
      integer na,na5,na4,na3,naf
      integer nsb,nu,nopb,nopd
      integer ndi,nti,nt,nt5,nt4
      integer npt,nbt,nat,ntt,nd
      integer nd5,nd4,nd3,nvp,nhb
      integer nmp,npi,npi5,npi4
      integer cls,atn,lig
      integer nx,ny,nxy
      integer bt,at,sbt,tt
      integer ft(6),pg(maxval)
      real*8 polfa, dampfa
      real*8 apret,bexpt
      real*8 alphaf 
      real*8 wght,rd,ep,rdn
      real*8 an1,an2,an3
      real*8 ba1,ba2
      real*8 aa1,aa2,aa3
      real*8 bt1,bt2,bt3
      real*8 bt4,bt5,bt6
      real*8 bt7,bt8,bt9
      real*8 at1,at2,at3
      real*8 at4,at5,at6
      real*8 an,pr,ds,dk
      real*8 vd,cg,dp,ps
      real*8 fc,bd,dl,el
      real*8 pt,pol,thl,dird
      real*8 iz,rp,ss,ts
      real*8 abc,cba
      real*8 gi,alphi
      real*8 nni,factor
      real*8 vt(6),st(6)
      real*8 pl(13)
      real*8 tx(maxtgrd2)
      real*8 ty(maxtgrd2)
      real*8 tf(maxtgrd2)
      logical header,swap
      character*1 da1
      character*4 pa,pb,pc
      character*4 pd,pe
      character*8 axt
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     initialize the counters for some parameter types
c
      nvp = 0
      nhb = 0
      nb = 0
      nb5 = 0
      nb4 = 0
      nb3 = 0
      nel = 0
      na = 0
      na5 = 0
      na4 = 0
      na3 = 0
      naf = 0
      nsb = 0
      nu = 0
      nopb = 0
      nopd = 0
      ndi = 0
      nti = 0
      nt = 0
      nt5 = 0
      nt4 = 0
      npt = 0
      nbt = 0
      nat = 0
      ntt = 0
      nd = 0
      nd5 = 0
      nd4 = 0
      nd3 = 0
      nmp = 0
      npi = 0
      npi5 = 0
      npi4 = 0
c
c     number of characters in an atom number text string
c
      size = 4
c
c     set blank line header before echoed comment lines
c
      header = .true.
c
c     process each line of the parameter file, first
c     extract the keyword at the start of each line
c
      iprm = 0
      do while (iprm .lt. nprm)
         iprm = iprm + 1
         record = prmline(iprm)
         next = 1
         call gettext (record,keyword,next)
         call upcase (keyword)
c
c     check for a force field modification keyword
c
         call prmkey (record)
c
c     comment line to be echoed to the output
c
         if (keyword(1:5) .eq. 'ECHO ') then
            string = record(next:240)
            length = trimtext (string)
            if (header) then
               header = .false.
               write (iout,10)
   10          format ()
            end if
            if (length .eq. 0) then
               write (iout,20)
   20          format ()
            else
               write (iout,30)  string(1:length)
   30          format (a)
            end if
c
c     atom type definitions and parameters
c
         else if (keyword(1:5) .eq. 'ATOM ') then
            ia = 0
            cls = 0
            atn = 0
            wght = 0.0d0
            lig = 0
            call getnumb (record,ia,next)
            call getnumb (record,cls,next)
            if (cls .eq. 0)  cls = ia
            atmcls(ia) = cls
            if (ia .ge. maxtyp) then
               write (iout,40)
   40          format (/,' READPRM  --  Too many Atom Types;',
     &                    ' Increase MAXTYP')
               call fatal
            else if (cls .ge. maxclass) then
               write (iout,50)
   50          format (/,' READPRM  --  Too many Atom Classes;',
     &                    ' Increase MAXCLASS')
               call fatal
            end if
            if (ia .ne. 0) then
               call gettext (record,symbol(ia),next)
               call getstring (record,describe(ia),next)
               string = record(next:240)
               read (string,*,err=60,end=60)  atn,wght,lig
   60          continue
               atmnum(ia) = atn
               weight(ia) = wght
               ligand(ia) = lig
            end if
c
c     van der Waals parameters for individual atom types
c
         else if (keyword(1:4) .eq. 'VDW ') then
            ia = 0
            rd = 0.0d0
            ep = 0.0d0
            rdn = 0.0d0
            string = record(next:240)
            read (string,*,err=70,end=70)  ia,rd,ep,rdn
   70       continue
            if (ia .ne. 0) then
               rad(ia) = rd
               eps(ia) = ep
               reduct(ia) = rdn
            end if
c
c     van der Waals 1-4 parameters for individual atom types
c
         else if (keyword(1:6) .eq. 'VDW14 ') then
            ia = 0
            rd = 0.0d0
            ep = 0.0d0
            string = record(next:240)
            read (string,*,err=80,end=80)  ia,rd,ep
   80       continue
            if (ia .ne. 0) then
               rad4(ia) = rd
               eps4(ia) = ep
            end if
c
c     van der Waals parameters for specific atom pairs
c
         else if (keyword(1:6) .eq. 'VDWPR ') then
            ia = 0
            ib = 0
            rd = 0.0d0
            ep = 0.0d0
            string = record(next:240)
            read (string,*,err=90,end=90)  ia,ib,rd,ep
   90       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            nvp = nvp + 1
            if (ia .le. ib) then
               kvpr(nvp) = pa//pb
            else
               kvpr(nvp) = pb//pa
            end if
            radpr(nvp) = rd
            epspr(nvp) = ep
c
c     van der Waals parameters for hydrogen bonding pairs
c
         else if (keyword(1:6) .eq. 'HBOND ') then
            ia = 0
            ib = 0
            rd = 0.0d0
            ep = 0.0d0
            string = record(next:240)
            read (string,*,err=100,end=100)  ia,ib,rd,ep
  100       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            nhb = nhb + 1
            if (ia .le. ib) then
               khb(nhb) = pa//pb
            else
               khb(nhb) = pb//pa
            end if
            radhb(nhb) = rd
            epshb(nhb) = ep
c
c    charge transfer parameters for individual atom types
c
         else if (keyword(1:3) .eq. 'CT ') then
            ia = 0
            apret = 0.0d0
            bexpt = 0.0d0
            string = record(next:120)
            read (string,*,err=105,end=105)  ia,apret,bexpt
  105       continue
            if (ia .ne. 0) then
               apre(ia) = apret 
               bexp(ia) = bexpt
            end if
c
c    charge penetration parameters for individual atom types
c
         else if (keyword(1:3) .eq. 'CP ') then
            ia = 0
            alphaf = 0.0d0
            string = record(next:120)
            read (string,*,err=108,end=108)  ia,alphaf
  108       continue
            if (ia .ne. 0) then
               apena(ia) = alphaf
            end if

c
c     bond stretching parameters
c
         else if (keyword(1:5) .eq. 'BOND ') then
            ia = 0
            ib = 0
            fc = 0.0d0
            bd = 0.0d0
            string = record(next:240)
            read (string,*,err=110,end=110)  ia,ib,fc,bd
  110       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            nb = nb + 1
            if (ia .le. ib) then
               kb(nb) = pa//pb
            else
               kb(nb) = pb//pa
            end if
            bcon(nb) = fc
            blen(nb) = bd
c
c     bond stretching parameters for 5-membered rings
c
         else if (keyword(1:6) .eq. 'BOND5 ') then
            ia = 0
            ib = 0
            fc = 0.0d0
            bd = 0.0d0
            string = record(next:240)
            read (string,*,err=120,end=120)  ia,ib,fc,bd
  120       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            nb5 = nb5 + 1
            if (ia .le. ib) then
               kb5(nb5) = pa//pb
            else
               kb5(nb5) = pb//pa
            end if
            bcon5(nb5) = fc
            blen5(nb5) = bd
c
c     bond stretching parameters for 4-membered rings
c
         else if (keyword(1:6) .eq. 'BOND4 ') then
            ia = 0
            ib = 0
            fc = 0.0d0
            bd = 0.0d0
            string = record(next:240)
            read (string,*,err=130,end=130)  ia,ib,fc,bd
  130       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            nb4 = nb4 + 1
            if (ia .le. ib) then
               kb4(nb4) = pa//pb
            else
               kb4(nb4) = pb//pa
            end if
            bcon4(nb4) = fc
            blen4(nb4) = bd
c
c     bond stretching parameters for 3-membered rings
c
         else if (keyword(1:6) .eq. 'BOND3 ') then
            ia = 0
            ib = 0
            fc = 0.0d0
            bd = 0.0d0
            string = record(next:240)
            read (string,*,err=140,end=140)  ia,ib,fc,bd
  140       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            nb3 = nb3 + 1
            if (ia .le. ib) then
               kb3(nb3) = pa//pb
            else
               kb3(nb3) = pb//pa
            end if
            bcon3(nb3) = fc
            blen3(nb3) = bd
c
c     electronegativity bond length correction parameters
c
         else if (keyword(1:9) .eq. 'ELECTNEG ') then
            ia = 0
            ib = 0
            ic = 0
            dl = 0.0d0
            string = record(next:240)
            read (string,*,err=150,end=150)  ia,ib,ic,dl
  150       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            nel = nel + 1
            if (ia .le. ic) then
               kel(nel) = pa//pb//pc
            else
               kel(nel) = pc//pb//pa
            end if
            dlen(nel) = dl
c
c     bond angle bending parameters
c
         else if (keyword(1:6) .eq. 'ANGLE ') then
            ia = 0
            ib = 0
            ic = 0
            fc = 0.0d0
            an1 = 0.0d0
            an2 = 0.0d0
            an3 = 0.0d0
            string = record(next:240)
            read (string,*,err=160,end=160)  ia,ib,ic,fc,an1,an2,an3
  160       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            na = na + 1
            if (ia .le. ic) then
               ka(na) = pa//pb//pc
            else
               ka(na) = pc//pb//pa
            end if
            acon(na) = fc
            ang(1,na) = an1
            ang(2,na) = an2
            ang(3,na) = an3
c
c     angle bending parameters for 5-membered rings
c
         else if (keyword(1:7) .eq. 'ANGLE5 ') then
            ia = 0
            ib = 0
            ic = 0
            fc = 0.0d0
            an1 = 0.0d0
            an2 = 0.0d0
            an3 = 0.0d0
            string = record(next:240)
            read (string,*,err=170,end=170)  ia,ib,ic,fc,an1,an2,an3
  170       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            na5 = na5 + 1
            if (ia .le. ic) then
               ka5(na5) = pa//pb//pc
            else
               ka5(na5) = pc//pb//pa
            end if
            acon5(na5) = fc
            ang5(1,na5) = an1
            ang5(2,na5) = an2
            ang5(3,na5) = an3
c
c     angle bending parameters for 4-membered rings
c
         else if (keyword(1:7) .eq. 'ANGLE4 ') then
            ia = 0
            ib = 0
            ic = 0
            fc = 0.0d0
            an1 = 0.0d0
            an2 = 0.0d0
            an3 = 0.0d0
            string = record(next:240)
            read (string,*,err=180,end=180)  ia,ib,ic,fc,an1,an2,an3
  180       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            na4 = na4 + 1
            if (ia .le. ic) then
               ka4(na4) = pa//pb//pc
            else
               ka4(na4) = pc//pb//pa
            end if
            acon4(na4) = fc
            ang4(1,na4) = an1
            ang4(2,na4) = an2
            ang4(3,na4) = an3
c
c     angle bending parameters for 3-membered rings
c
         else if (keyword(1:7) .eq. 'ANGLE3 ') then
            ia = 0
            ib = 0
            ic = 0
            fc = 0.0d0
            an1 = 0.0d0
            an2 = 0.0d0
            an3 = 0.0d0
            string = record(next:240)
            read (string,*,err=190,end=190)  ia,ib,ic,fc,an1,an2,an3
  190       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            na3 = na3 + 1
            if (ia .le. ic) then
               ka3(na3) = pa//pb//pc
            else
               ka3(na3) = pc//pb//pa
            end if
            acon3(na3) = fc
            ang3(1,na3) = an1
            ang3(2,na3) = an2
            ang3(3,na3) = an3
c
c     Fourier bond angle bending parameters
c
         else if (keyword(1:7) .eq. 'ANGLEF ') then
            ia = 0
            ib = 0
            ic = 0
            fc = 0.0d0
            an = 0.0d0
            pr = 0.0d0
            string = record(next:240)
            read (string,*,err=200,end=200)  ia,ib,ic,fc,an,pr
  200       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            naf = naf + 1
            if (ia .le. ic) then
               kaf(naf) = pa//pb//pc
            else
               kaf(naf) = pc//pb//pa
            end if
            aconf(naf) = fc
            angf(1,naf) = an
            angf(2,naf) = pr
c
c     stretch-bend parameters
c
         else if (keyword(1:7) .eq. 'STRBND ') then
            ia = 0
            ib = 0
            ic = 0
            ba1 = 0.0d0
            ba2 = 0.0d0
            string = record(next:240)
            read (string,*,err=210,end=210)  ia,ib,ic,ba1,ba2
  210       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            nsb = nsb + 1
            if (ia .le. ic) then
               ksb(nsb) = pa//pb//pc
               stbn(1,nsb) = ba1
               stbn(2,nsb) = ba2
            else
               ksb(nsb) = pc//pb//pa
               stbn(1,nsb) = ba2
               stbn(2,nsb) = ba1
            end if
c
c     Urey-Bradley parameters
c
         else if (keyword(1:9) .eq. 'UREYBRAD ') then
            ia = 0
            ib = 0
            ic = 0
            fc = 0.0d0
            ds = 0.0d0
            string = record(next:240)
            read (string,*,err=220,end=220)  ia,ib,ic,fc,ds
  220       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            nu = nu + 1
            if (ia .le. ic) then
               ku(nu) = pa//pb//pc
            else
               ku(nu) = pc//pb//pa
            end if
            ucon(nu) = fc
            dst13(nu) = ds
c
c     angle-angle parameters
c
         else if (keyword(1:7) .eq. 'ANGANG ') then
            ia = 0
            aa1 = 0.0d0
            aa2 = 0.0d0
            aa3 = 0.0d0
            string = record(next:240)
            read (string,*,err=230,end=230)  ia,aa1,aa2,aa3
  230       continue
            if (ia .ne. 0) then
               anan(1,ia) = aa1
               anan(2,ia) = aa2
               anan(3,ia) = aa3
            end if
c
c     out-of-plane bend parameters
c
         else if (keyword(1:7) .eq. 'OPBEND ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            fc = 0.0d0
            string = record(next:240)
            read (string,*,err=240,end=240)  ia,ib,ic,id,fc
  240       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            nopb = nopb + 1
            if (ic .le. id) then
               kopb(nopb) = pa//pb//pc//pd
            else
               kopb(nopb) = pa//pb//pd//pc
            end if
            opbn(nopb) = fc
c
c     out-of-plane distance parameters
c
         else if (keyword(1:7) .eq. 'OPDIST ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            fc = 0.0d0
            string = record(next:240)
            read (string,*,err=250,end=250)  ia,ib,ic,id,fc
  250       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            nopd = nopd + 1
            imin = min(ib,ic,id)
            if (ib .eq. imin) then
               if (ic .le. id) then
                  kopd(nopd) = pa//pb//pc//pd
               else
                  kopd(nopd) = pa//pb//pd//pc
               end if
            else if (ic .eq. imin) then
               if (ib .le. id) then
                  kopd(nopd) = pa//pc//pb//pd
               else
                  kopd(nopd) = pa//pc//pd//pb
               end if
            else if (id .eq. imin) then
               if (ib .le. ic) then
                  kopd(nopd) = pa//pd//pb//pc
               else
                  kopd(nopd) = pa//pd//pc//pb
               end if
            end if
            opds(nopd) = fc
c
c     improper dihedral parameters
c
         else if (keyword(1:9) .eq. 'IMPROPER ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            dk = 0.0d0
            vd = 0.0d0
            string = record(next:240)
            read (string,*,err=260,end=260)  ia,ib,ic,id,dk,vd
  260       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            ndi = ndi + 1
            kdi(ndi) = pa//pb//pc//pd
            dcon(ndi) = dk
            tdi(ndi) = vd
c
c     improper torsional parameters
c
         else if (keyword(1:8) .eq. 'IMPTORS ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            do i = 1, 6
               vt(i) = 0.0d0
               st(i) = 0.0d0
               ft(i) = 0
            end do
            string = record(next:240)
            read (string,*,err=270,end=270)  ia,ib,ic,id,
     &                                       (vt(j),st(j),ft(j),j=1,6)
  270       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            nti = nti + 1
            kti(nti) = pa//pb//pc//pd
            call torphase (ft,vt,st)
            ti1(1,nti) = vt(1)
            ti1(2,nti) = st(1)
            ti2(1,nti) = vt(2)
            ti2(2,nti) = st(2)
            ti3(1,nti) = vt(3)
            ti3(2,nti) = st(3)
c
c     torsional parameters
c
         else if (keyword(1:8) .eq. 'TORSION ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            do i = 1, 6
               vt(i) = 0.0d0
               st(i) = 0.0d0
               ft(i) = 0
            end do
            string = record(next:240)
            read (string,*,err=280,end=280)  ia,ib,ic,id,
     &                                       (vt(j),st(j),ft(j),j=1,6)
  280       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            nt = nt + 1
            if (ib .lt. ic) then
               kt(nt) = pa//pb//pc//pd
            else if (ic .lt. ib) then
               kt(nt) = pd//pc//pb//pa
            else if (ia .le. id) then
               kt(nt) = pa//pb//pc//pd
            else if (id .lt. ia) then
               kt(nt) = pd//pc//pb//pa
            end if
            call torphase (ft,vt,st)
            t1(1,nt) = vt(1)
            t1(2,nt) = st(1)
            t2(1,nt) = vt(2)
            t2(2,nt) = st(2)
            t3(1,nt) = vt(3)
            t3(2,nt) = st(3)
            t4(1,nt) = vt(4)
            t4(2,nt) = st(4)
            t5(1,nt) = vt(5)
            t5(2,nt) = st(5)
            t6(1,nt) = vt(6)
            t6(2,nt) = st(6)
c
c     torsional parameters for 5-membered rings
c
         else if (keyword(1:9) .eq. 'TORSION5 ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            do i = 1, 6
               vt(i) = 0.0d0
               st(i) = 0.0d0
               ft(i) = 0
            end do
            string = record(next:240)
            read (string,*,err=290,end=290)  ia,ib,ic,id,
     &                                       (vt(j),st(j),ft(j),j=1,6)
  290       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            nt5 = nt5 + 1
            if (ib .lt. ic) then
               kt5(nt5) = pa//pb//pc//pd
            else if (ic .lt. ib) then
               kt5(nt5) = pd//pc//pb//pa
            else if (ia .le. id) then
               kt5(nt5) = pa//pb//pc//pd
            else if (id .lt. ia) then
               kt5(nt5) = pd//pc//pb//pa
            end if
            call torphase (ft,vt,st)
            t15(1,nt5) = vt(1)
            t15(2,nt5) = st(1)
            t25(1,nt5) = vt(2)
            t25(2,nt5) = st(2)
            t35(1,nt5) = vt(3)
            t35(2,nt5) = st(3)
            t45(1,nt5) = vt(4)
            t45(2,nt5) = st(4)
            t55(1,nt5) = vt(5)
            t55(2,nt5) = st(5)
            t65(1,nt5) = vt(6)
            t65(2,nt5) = st(6)
c
c     torsional parameters for 4-membered rings
c
         else if (keyword(1:9) .eq. 'TORSION4 ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            do i = 1, 6
               vt(i) = 0.0d0
               st(i) = 0.0d0
               ft(i) = 0
            end do
            string = record(next:240)
            read (string,*,err=300,end=300)  ia,ib,ic,id,
     &                                       (vt(i),st(i),ft(i),i=1,6)
  300       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            nt4 = nt4 + 1
            if (ib .lt. ic) then
               kt4(nt4) = pa//pb//pc//pd
            else if (ic .lt. ib) then
               kt4(nt4) = pd//pc//pb//pa
            else if (ia .le. id) then
               kt4(nt4) = pa//pb//pc//pd
            else if (id .lt. ia) then
               kt4(nt4) = pd//pc//pb//pa
            end if
            call torphase (ft,vt,st)
            t14(1,nt4) = vt(1)
            t14(2,nt4) = st(1)
            t24(1,nt4) = vt(2)
            t24(2,nt4) = st(2)
            t34(1,nt4) = vt(3)
            t34(2,nt4) = st(3)
            t44(1,nt4) = vt(4)
            t44(2,nt4) = st(4)
            t54(1,nt4) = vt(5)
            t54(2,nt4) = st(5)
            t64(1,nt4) = vt(6)
            t64(2,nt4) = st(6)
c
c     pi-system torsion parameters
c
         else if (keyword(1:7) .eq. 'PITORS ') then
            ia = 0
            ib = 0
            pt = 0.0d0
            string = record(next:240)
            read (string,*,err=310,end=310)  ia,ib,pt
  310       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            npt = npt + 1
            if (ia .le. ib) then
               kpt(npt) = pa//pb
            else
               kpt(npt) = pb//pa
            end if
            ptcon(npt) = pt
c
c     stretch-torsion parameters
c
         else if (keyword(1:8) .eq. 'STRTORS ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            bt1 = 0.0d0
            bt2 = 0.0d0
            bt3 = 0.0d0
            bt4 = 0.0d0
            bt5 = 0.0d0
            bt6 = 0.0d0
            bt7 = 0.0d0
            bt8 = 0.0d0
            bt9 = 0.0d0
            string = record(next:240)
            read (string,*,err=320,end=320)  ia,ib,ic,id,bt1,bt2,bt3,
     &                                       bt4,bt5,bt6,bt7,bt8,bt9
  320       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            nbt = nbt + 1
            if (ib .lt. ic) then
               kbt(nbt) = pa//pb//pc//pd
               swap = .false.
            else if (ic .lt. ib) then
               kbt(nbt) = pd//pc//pb//pa
               swap = .true.
            else if (ia .le. id) then
               kbt(nbt) = pa//pb//pc//pd
               swap = .false.
            else if (id .lt. ia) then
               kbt(nbt) = pd//pc//pb//pa
               swap = .true.
            end if
            btcon(4,nbt) = bt4
            btcon(5,nbt) = bt5
            btcon(6,nbt) = bt6
            if (swap) then
               btcon(1,nbt) = bt7
               btcon(2,nbt) = bt8
               btcon(3,nbt) = bt9
               btcon(7,nbt) = bt1
               btcon(8,nbt) = bt2
               btcon(9,nbt) = bt3
            else
               btcon(1,nbt) = bt1
               btcon(2,nbt) = bt2
               btcon(3,nbt) = bt3
               btcon(7,nbt) = bt7
               btcon(8,nbt) = bt8
               btcon(9,nbt) = bt9
            end if
c
c     angle-torsion parameters
c
         else if (keyword(1:8) .eq. 'ANGTORS ') then
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
            read (string,*,err=330,end=330)  ia,ib,ic,id,at1,at2,
     &                                       at3,at4,at5,at6
  330       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            nat = nat + 1
            if (ib .lt. ic) then
               kat(nat) = pa//pb//pc//pd
               swap = .false.
            else if (ic .lt. ib) then
               kat(nat) = pd//pc//pb//pa
               swap = .true.
            else if (ia .le. id) then
               kat(nat) = pa//pb//pc//pd
               swap = .false.
            else if (id .lt. ia) then
               kat(nat) = pd//pc//pb//pa
               swap = .true.
            end if
            if (swap) then
               atcon(1,nat) = at4
               atcon(2,nat) = at5
               atcon(3,nat) = at6
               atcon(4,nat) = at1
               atcon(5,nat) = at2
               atcon(6,nat) = at3
            else
               atcon(1,nat) = at1
               atcon(2,nat) = at2
               atcon(3,nat) = at3
               atcon(4,nat) = at4
               atcon(5,nat) = at5
               atcon(6,nat) = at6
            end if
c
c     torsion-torsion parameters
c
         else if (keyword(1:8) .eq. 'TORTORS ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            ie = 0
            nx = 0
            ny = 0
            nxy = 0
            do i = 1, maxtgrd2
               tx(i) = 0.0d0
               ty(i) = 0.0d0
               tf(i) = 0.0d0
            end do
            string = record(next:240)
            read (string,*,err=340,end=340)  ia,ib,ic,id,ie,nx,ny
            nxy = nx * ny
            do i = 1, nxy
               iprm = iprm + 1
               record = prmline(iprm)
               read (record,*,err=340,end=340)  tx(i),ty(i),tf(i)
            end do
  340       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            call numeral (ie,pe,size)
            ntt = ntt + 1
            ktt(ntt) = pa//pb//pc//pd//pe
            nx = nxy
            call sort9 (nx,tx)
            ny = nxy
            call sort9 (ny,ty)
            tnx(ntt) = nx
            tny(ntt) = ny
            do i = 1, nx
               ttx(i,ntt) = tx(i)
            end do
            do i = 1, ny
               tty(i,ntt) = ty(i)
            end do
            do i = 1, nxy
               tbf(i,ntt) = tf(i)
            end do
c
c     atomic partial charge parameters
c
         else if (keyword(1:7) .eq. 'CHARGE ') then
            ia = 0
            cg = 0.0d0
            string = record(next:240)
            read (string,*,err=350,end=350)  ia,cg
  350       continue
            if (ia .ne. 0)  chg(ia) = cg
c
c     bond dipole moment parameters
c
         else if (keyword(1:7) .eq. 'DIPOLE ') then
            ia = 0
            ib = 0
            dp = 0.0d0
            ps = 0.5d0
            string = record(next:240)
            read (string,*,err=360,end=360)  ia,ib,dp,ps
  360       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            nd = nd + 1
            if (ia .le. ib) then
               kd(nd) = pa//pb
            else
               kd(nd) = pb//pa
            end if
            dpl(nd) = dp
            pos(nd) = ps
c
c     bond dipole moment parameters for 5-membered rings
c
         else if (keyword(1:8) .eq. 'DIPOLE5 ') then
            ia = 0
            ib = 0
            dp = 0.0d0
            ps = 0.5d0
            string = record(next:240)
            read (string,*,err=370,end=370)  ia,ib,dp,ps
  370       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            nd5 = nd5 + 1
            if (ia .le. ib) then
               kd5(nd5) = pa//pb
            else
               kd5(nd5) = pb//pa
            end if
            dpl5(nd5) = dp
            pos5(nd5) = ps
c
c     bond dipole moment parameters for 4-membered rings
c
         else if (keyword(1:8) .eq. 'DIPOLE4 ') then
            ia = 0
            ib = 0
            dp = 0.0d0
            ps = 0.5d0
            string = record(next:240)
            read (string,*,err=380,end=380)  ia,ib,dp,ps
  380       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            nd4 = nd4 + 1
            if (ia .le. ib) then
               kd4(nd4) = pa//pb
            else
               kd4(nd4) = pb//pa
            end if
            dpl4(nd4) = dp
            pos4(nd4) = ps
c
c     bond dipole moment parameters for 3-membered rings
c
         else if (keyword(1:8) .eq. 'DIPOLE3 ') then
            ia = 0
            ib = 0
            dp = 0.0d0
            ps = 0.5d0
            string = record(next:240)
            read (string,*,err=390,end=390)  ia,ib,dp,ps
  390       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            nd3 = nd3 + 1
            if (ia .le. ib) then
               kd3(nd3) = pa//pb
            else
               kd3(nd3) = pb//pa
            end if
            dpl3(nd3) = dp
            pos3(nd3) = ps
c
c     atomic multipole moment parameters
c
         else if (keyword(1:10) .eq. 'MULTIPOLE ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            axt = 'Z-then-X'
            do i = 1, 13
               pl(i) = 0.0d0
            end do
            string = record(next:240)
            read (string,*,err=400,end=400)  ia,ib,ic,id,pl(1)
            goto 430
  400       continue
            id = 0
            read (string,*,err=410,end=410)  ia,ib,ic,pl(1)
            goto 430
  410       continue
            ic = 0
            read (string,*,err=420,end=420)  ia,ib,pl(1)
            goto 430
  420       continue
            ib = 0
            read (string,*,err=440,end=440)  ia,pl(1)
  430       continue
            iprm = iprm + 1
            record = prmline(iprm)
            read (record,*,err=440,end=440)  pl(2),pl(3),pl(4)
            iprm = iprm + 1
            record = prmline(iprm)
            read (record,*,err=440,end=440)  pl(5)
            iprm = iprm + 1
            record = prmline(iprm)
            read (record,*,err=440,end=440)  pl(8),pl(9)
            iprm = iprm + 1
            record = prmline(iprm)
            read (record,*,err=440,end=440)  pl(11),pl(12),pl(13)
  440       continue
            if (ib .eq. 0)  axt = 'None'
            if (ib.ne.0 .and. ic.eq.0)  axt = 'Z-Only'
            if (ib.lt.0 .or. ic.lt.0)  axt = 'Bisector'
            if (ic.lt.0 .and. id.lt.0)  axt = 'Z-Bisect'
            if (max(ib,ic,id) .lt. 0)  axt = '3-Fold'
            ib = abs(ib)
            ic = abs(ic)
            id = abs(id)
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            nmp = nmp + 1
            kmp(nmp) = pa//pb//pc//pd
            mpaxis(nmp) = axt
            multip(1,nmp) = pl(1)
            multip(2,nmp) = pl(2)
            multip(3,nmp) = pl(3)
            multip(4,nmp) = pl(4)
            multip(5,nmp) = pl(5)
            multip(6,nmp) = pl(8)
            multip(7,nmp) = pl(11)
            multip(8,nmp) = pl(8)
            multip(9,nmp) = pl(9)
            multip(10,nmp) = pl(12)
            multip(11,nmp) = pl(11)
            multip(12,nmp) = pl(12)
            multip(13,nmp) = pl(13)
         
c
c     atomic dipole polarizability parameters
c
         else if (keyword(1:9) .eq. 'POLARIZE ') then
            ia = 0
            pol = 0.0d0
            thl = 0.0d0
            dird = 0.0d0
            do i = 1, maxval
               pg(i) = 0
            end do
            string = record(next:240)
            read (string,*,err=450,end=450) ia,pol,thl,dird,
     &                                       (pg(i),i=1,maxval)
  450       continue
            if (ia .ne. 0) then
               polr(ia) = pol
               athl(ia) = thl
               adird(ia) = dird
               do i = 1, maxval
                  pgrp(i,ia) = pg(i)
               end do
            end if
c
c     conjugated pisystem atom parameters
c
         else if (keyword(1:7) .eq. 'PIATOM ') then
            ia = 0
            el = 0.0d0
            iz = 0.0d0
            rp = 0.0d0
            string = record(next:240)
            read (string,*,err=460,end=460)  ia,el,iz,rp
  460       continue
            if (ia .ne. 0) then
               electron(ia) = el
               ionize(ia) = iz
               repulse(ia) = rp
            end if
c
c     conjugated pisystem bond parameters
c
         else if (keyword(1:7) .eq. 'PIBOND ') then
            ia = 0
            ib = 0
            ss = 0.0d0
            ts = 0.0d0
            string = record(next:240)
            read (string,*,err=470,end=470)  ia,ib,ss,ts
  470       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            npi = npi + 1
            if (ia .le. ib) then
               kpi(npi) = pa//pb
            else
               kpi(npi) = pb//pa
            end if
            sslope(npi) = ss
            tslope(npi) = ts
c
c     conjugated pisystem bond parameters for 5-membered rings
c
         else if (keyword(1:8) .eq. 'PIBOND5 ') then
            ia = 0
            ib = 0
            ss = 0.0d0
            ts = 0.0d0
            string = record(next:240)
            read (string,*,err=480,end=480)  ia,ib,ss,ts
  480       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            npi5 = npi5 + 1
            if (ia .le. ib) then
               kpi5(npi5) = pa//pb
            else
               kpi5(npi5) = pb//pa
            end if
            sslope5(npi5) = ss
            tslope5(npi5) = ts
c
c     conjugated pisystem bond parameters for 4-membered rings
c
         else if (keyword(1:8) .eq. 'PIBOND4 ') then
            ia = 0
            ib = 0
            ss = 0.0d0
            ts = 0.0d0
            string = record(next:240)
            read (string,*,err=490,end=490)  ia,ib,ss,ts
  490       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            npi4 = npi4 + 1
            if (ia .le. ib) then
               kpi4(npi4) = pa//pb
            else
               kpi4(npi4) = pb//pa
            end if
            sslope4(npi4) = ss
            tslope4(npi4) = ts
c
c     metal ligand field splitting parameters
c
         else if (keyword(1:6) .eq. 'METAL ') then
            string = record(next:240)
            read (string,*,err=500,end=500)  ia
  500       continue
c
c     biopolymer atom type conversion definitions
c
         else if (keyword(1:8) .eq. 'BIOTYPE ') then
            ia = 0
            ib = 0
            string = record(next:240)
            read (string,*,err=510,end=510)  ia
            call getword (record,string,next)
            call getstring (record,string,next)
            string = record(next:240)
            read (string,*,err=510,end=510)  ib
  510       continue
            if (ia .ge. maxbio) then
               write (iout,520)
  520          format (/,' READPRM  --  Too many Biopolymer Types;',
     &                    ' Increase MAXBIO')
               call fatal
            end if
            if (ia .ne. 0)  biotyp(ia) = ib
c
c     MMFF atom class equivalency parameters
c
         else if (keyword(1:10) .eq. 'MMFFEQUIV ') then
            string = record(next:240)
            ia = 1000
            ib = 1000
            ic = 1000
            id = 1000
            ie = 1000
            if = 0
            read (string,*,err=530,end=530)  ia,ib,ic,id,ie,if
  530       continue
            eqclass(if,1) = ia
            eqclass(if,2) = ib
            eqclass(if,3) = ic
            eqclass(if,4) = id
            eqclass(if,5) = ie
c
c     MMFF covalent radius and electronegativity parameters
c
         else if (keyword(1:11) .eq. 'MMFFCOVRAD ') then
            ia = 0
            fc = 0.0d0
            bd = 0.0d0
            string = record(next:240)
            read (string,*,err=540,end=540)  ia,fc,bd
  540       continue
            rad0(ia) = fc
            paulel(ia) = bd
c
c     MMFF atom class property parameters
c
         else if (keyword(1:9) .eq. 'MMFFPROP ') then
            string = record(next:240)
            ia = 1000
            ib = 1000
            ic = 1000
            id = 1000
            ie = 1000
            if = 1000
            ig = 1000
            ih = 1000
            ii = 1000
            read (string,*,err=550,end=550)  ia,ib,ic,id,ie,
     &                                       if,ig,ih,ii
  550       continue
            crd(ia) = ic
            val(ia) = id
            pilp(ia) = ie
            mltb(ia) = if
            arom(ia) = ig
            lin(ia) = ih
            sbmb(ia) = ii
c
c     MMFF van der Waals parameters
c
         else if (keyword(1:8) .eq. 'MMFFVDW ') then
            ia = 0
            rd = 0.0d0
            ep = 0.0d0
            rdn = 0.0d0
            da1 = 'C'
            string = record(next:240)
            read (string,*,err=560,end=560)  ia,rd,alphi,nni,gi,da1
  560       continue
            if (ia .ne. 0) then
               rad(ia) = rd
               g(ia) = gi
               alph(ia) = alphi
               nn(ia) = nni
               da(ia) = da1
            end if
c
c     MMFF bond stretching parameters
c
         else if (keyword(1:9) .eq. 'MMFFBOND ') then
            ia = 0
            ib = 0
            fc = 0.0d0
            bd = 0.0d0
            bt = 2
            string = record(next:240)
            read (string,*,err=570,end=570)  ia,ib,fc,bd,bt
  570       continue
            nb = nb + 1
            if (bt .eq. 0) then
               mmff_kb(ia,ib) = fc
               mmff_kb(ib,ia) = fc
               mmff_b0(ia,ib) = bd
               mmff_b0(ib,ia) = bd
            else if (bt .eq. 1) then
               mmff_kb1(ia,ib) = fc
               mmff_kb1(ib,ia) = fc
               mmff_b1(ia,ib) = bd
               mmff_b1(ib,ia) = bd
            end if
c
c     MMFF bond stretching empirical rule parameters
c
         else if (keyword(1:11) .eq. 'MMFFBONDER ') then
            ia = 0
            ib = 0
            fc = 0.0d0
            bd = 0.0d0
            string = record(next:240)
            read (string,*,err=580,end=580)  ia,ib,fc,bd
  580       continue
            r0ref(ia,ib) = fc
            r0ref(ib,ia) = fc
            kbref(ia,ib) = bd
            kbref(ib,ia) = bd
c
c     MMFF bond angle bending parameters
c
         else if (keyword(1:10) .eq. 'MMFFANGLE ') then
            ia = 0
            ib = 0
            ic = 0
            fc = 0.0d0
            an1 = 0.0d0
            at = 3
            string = record(next:240)
            read (string,*,err=590,end=590)  ia,ib,ic,fc,an1,at
  590       continue
            na = na + 1
            if (an1 .ne. 0.0d0) then
               if (at .eq. 0) then
                  mmff_ka(ia,ib,ic) = fc
                  mmff_ka(ic,ib,ia) = fc
                  mmff_ang0(ia,ib,ic) = an1
                  mmff_ang0(ic,ib,ia) = an1
               else if (at .eq. 1) then
                  mmff_ka1(ia,ib,ic) = fc
                  mmff_ka1(ic,ib,ia) = fc
                  mmff_ang1(ia,ib,ic) = an1
                  mmff_ang1(ic,ib,ia) = an1
               else if (at .eq. 2) then
                  mmff_ka2(ia,ib,ic) = fc
                  mmff_ka2(ic,ib,ia) = fc
                  mmff_ang2(ia,ib,ic) = an1
                  mmff_ang2(ic,ib,ia) = an1
               else if (at .eq. 3) then
                  mmff_ka3(ia,ib,ic) = fc
                  mmff_ka3(ic,ib,ia) = fc
                  mmff_ang3(ia,ib,ic) = an1
                  mmff_ang3(ic,ib,ia) = an1
               else if (at .eq. 4) then
                  mmff_ka4(ia,ib,ic) = fc
                  mmff_ka4(ic,ib,ia) = fc
                  mmff_ang4(ia,ib,ic) = an1
                  mmff_ang4(ic,ib,ia) = an1
               else if (at .eq. 5) then
                  mmff_ka5(ia,ib,ic) = fc
                  mmff_ka5(ic,ib,ia) = fc
                  mmff_ang5(ia,ib,ic) = an1
                  mmff_ang5(ic,ib,ia) = an1
               else if (at .eq. 6) then
                  mmff_ka6(ia,ib,ic) = fc
                  mmff_ka6(ic,ib,ia) = fc
                  mmff_ang6(ia,ib,ic) = an1
                  mmff_ang6(ic,ib,ia) = an1
               else if (at .eq. 7) then
                  mmff_ka7(ia,ib,ic) = fc
                  mmff_ka7(ic,ib,ia) = fc
                  mmff_ang7(ia,ib,ic) = an1
                  mmff_ang7(ic,ib,ia) = an1
               else if (at .eq. 8) then
                  mmff_ka8(ia,ib,ic) = fc
                  mmff_ka8(ic,ib,ia) = fc
                  mmff_ang8(ia,ib,ic) = an1
                  mmff_ang8(ic,ib,ia) = an1
               end if
            end if
c
c     MMFF stretch-bend parameters
c
         else if (keyword(1:11) .eq. 'MMFFSTRBND ') then
            ia = 0
            ib = 0
            ic = 0
            abc = 0.0d0
            cba = 0.0d0
            sbt = 4
            string = record(next:240)
            read (string,*,err=600,end=600)  ia,ib,ic,abc,cba,sbt
  600       continue
            if (ia .ne. 0) then
               if (sbt .eq. 0) then
                  stbn_abc(ia,ib,ic) = abc
                  if (ic .ne. ia)  stbn_abc(ic,ib,ia) = cba
                  stbn_cba(ia,ib,ic) = cba
                  if (ic .ne. ia)  stbn_cba(ic,ib,ia) = abc
               else if (sbt .eq. 1) then
                  stbn_abc1(ia,ib,ic) = abc
                  if (ic .ne. ia)  stbn_abc1(ic,ib,ia) = cba
                  stbn_cba1(ia,ib,ic) = cba
                  if (ic .ne. ia)  stbn_cba1(ic,ib,ia) = abc
               else if (sbt .eq. 2) then
                  stbn_abc2(ia,ib,ic) = abc
                  if (ic .ne. ia)  stbn_abc2(ic,ib,ia) = cba
                  stbn_cba2(ia,ib,ic) = cba
                  if (ic .ne. ia)  stbn_cba2(ic,ib,ia) = abc
               else if (sbt .eq. 3) then
                  stbn_abc3(ia,ib,ic) = abc
                  if (ic .ne. ia)  stbn_abc3(ic,ib,ia) = cba
                  stbn_cba3(ia,ib,ic) = cba
                  if (ic .ne. ia)  stbn_cba3(ic,ib,ia) = abc
               else if (sbt .eq. 4) then
                  stbn_abc4(ia,ib,ic) = abc
                  if (ic .ne. ia)  stbn_abc4(ic,ib,ia) = cba
                  stbn_cba4(ia,ib,ic) = cba
                  if (ic .ne. ia)  stbn_cba4(ic,ib,ia) = abc
               else if (sbt .eq. 5) then
                  stbn_abc5(ia,ib,ic) = abc
                  if (ic .ne. ia)  stbn_abc5(ic,ib,ia) = cba
                  stbn_cba5(ia,ib,ic) = cba
                  if (ic .ne. ia)  stbn_cba5(ic,ib,ia) = abc
               else if (sbt .eq. 6) then
                  stbn_abc6(ia,ib,ic) = abc
                  if (ic .ne. ia)  stbn_abc6(ic,ib,ia) = cba
                  stbn_cba6(ia,ib,ic) = cba
                  if (ic .ne. ia)  stbn_cba6(ic,ib,ia) = abc
               else if (sbt .eq. 7) then
                  stbn_abc7(ia,ib,ic) = abc
                  if (ic .ne. ia)  stbn_abc7(ic,ib,ia) = cba
                  stbn_cba7(ia,ib,ic) = cba
                  if (ic .ne. ia)  stbn_cba7(ic,ib,ia) = abc
               else if (sbt .eq. 8) then
                  stbn_abc8(ia,ib,ic) = abc
                  if (ic .ne. ia)  stbn_abc8(ic,ib,ia) = cba
                  stbn_cba8(ia,ib,ic) = cba
                  if (ic .ne. ia)  stbn_cba8(ic,ib,ia) = abc
               else if (sbt .eq. 9) then
                  stbn_abc9(ia,ib,ic) = abc
                  if (ic .ne. ia)  stbn_abc9(ic,ib,ia) = cba
                  stbn_cba9(ia,ib,ic) = cba
                  if (ic .ne. ia)  stbn_cba9(ic,ib,ia) = abc
               else if (sbt .eq. 10) then
                  stbn_abc10(ia,ib,ic) = abc
                  if (ic .ne. ia)  stbn_abc10(ic,ib,ia) = cba
                  stbn_cba10(ia,ib,ic) = cba
                  if (ic .ne. ia)  stbn_cba10(ic,ib,ia) = abc
               else if (sbt .eq. 11) then
                  stbn_abc11(ia,ib,ic) = abc
                  if (ic .ne. ia)  stbn_abc11(ic,ib,ia) = cba
                  stbn_cba11(ia,ib,ic) = cba
                  if (ic .ne. ia)  stbn_cba11(ic,ib,ia) = abc
               end if
            end if
c
c     MMFF default stretch-bend parameters
c
         else if (keyword(1:12) .eq. 'MMFFDEFSTBN ') then
            string = record(next:240)
            ia = 1000
            ib = 1000
            ic = 1000
            abc = 0.0d0
            cba = 0.0d0
            read (string,*,err=610,end=610)  ia,ib,ic,abc,cba
  610       continue
            defstbn_abc(ia,ib,ic) = abc
            defstbn_cba(ia,ib,ic) = cba
            defstbn_abc(ic,ib,ia) = cba
            defstbn_cba(ic,ib,ia) = abc
c
c     MMFF out-of-plane bend parameters
c
         else if (keyword(1:11) .eq. 'MMFFOPBEND ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            fc = 0.0d0
            string = record(next:240)
            read (string,*,err=620,end=620)  ia,ib,ic,id,fc
  620       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            nopb = nopb + 1
            if (ic .le. id) then
               kopb(nopb) = pa//pb//pc//pd
            else
               kopb(nopb) = pa//pb//pd//pc
            end if
            opbn(nopb) = fc
c           if (ic.gt.0 .or. id.gt.0) then
c              nopb = nopb + 1
c              if (ib .le. id) then
c                 kopb(nopb) = pc//pb//pb//pd
c              else
c                 kopb(nopb) = pc//pb//pd//pb
c              end if
c              opbn(nopb) = fc
c              nopb = nopb + 1
c              if (ia .le. ic) then
c                 kopb(nopb) = pd//pb//pa//pc
c              else
c                 kopb(nopb) = pd//pb//pc//pa
c              end if
c              opbn(nopb) = fc
c           end if
c
c     MMFF torsional parameters
c
         else if (keyword(1:12) .eq. 'MMFFTORSION ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            do i = 1, 6
               vt(i) = 0.0d0
               st(i) = 0.0d0
               ft(i) = 0
            end do
            tt = 3
            string = record(next:240)
            read (string,*,err=630,end=630)  ia,ib,ic,id,(vt(j),
     &                                       st(j),ft(j),j=1,3),tt
  630       continue
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            nt = nt + 1
            if (tt .eq. 0) then
               if (ib .lt. ic) then
                  kt(nt) = pa//pb//pc//pd
               else if (ic .lt. ib) then
                  kt(nt) = pd//pc//pb//pa
               else if (ia .le. id) then
                  kt(nt) = pa//pb//pc//pd
               else if (id .lt. ia) then
                  kt(nt) = pd//pc//pb//pa
               end if
               call torphase (ft,vt,st)
               t1(1,nt) = vt(1)
               t1(2,nt) = st(1)
               t2(1,nt) = vt(2)
               t2(2,nt) = st(2)
               t3(1,nt) = vt(3)
               t3(2,nt) = st(3)
            else if (tt .eq. 1) then
               if (ib .lt. ic) then
                  kt_1(nt) = pa//pb//pc//pd
               else if (ic .lt. ib) then
                  kt_1(nt) = pd//pc//pb//pa
               else if (ia .le. id) then
                  kt_1(nt) = pa//pb//pc//pd
               else if (id .lt. ia) then
                  kt_1(nt) = pd//pc//pb//pa
               end if
               call torphase (ft,vt,st)
               t1_1(1,nt) = vt(1)
               t1_1(2,nt) = st(1)
               t2_1(1,nt) = vt(2)
               t2_1(2,nt) = st(2)
               t3_1(1,nt) = vt(3)
               t3_1(2,nt) = st(3)
            else if (tt .eq. 2) then
               if (ib .lt. ic) then
                  kt_2(nt) = pa//pb//pc//pd
               else if (ic .lt. ib) then
                  kt_2(nt) = pd//pc//pb//pa
               else if (ia .le. id) then
                  kt_2(nt) = pa//pb//pc//pd
               else if (id .lt. ia) then
                  kt_2(nt) = pd//pc//pb//pa
               end if
               call torphase (ft,vt,st)
               t1_2(1,nt) = vt(1)
               t1_2(2,nt) = st(1)
               t2_2(1,nt) = vt(2)
               t2_2(2,nt) = st(2)
               t3_2(1,nt) = vt(3)
               t3_2(2,nt) = st(3)
            else if (tt .eq. 4) then
               nt4 = nt4 + 1
               if (ib .lt. ic) then
                  kt4(nt4) = pa//pb//pc//pd
               else if (ic .lt. ib) then
                  kt4(nt4) = pd//pc//pb//pa
               else if (ia .le. id) then
                  kt4(nt4) = pa//pb//pc//pd
               else if (id .lt. ia) then
                  kt4(nt4) = pd//pc//pb//pa
               end if
               call torphase (ft,vt,st)
               t14(1,nt4) = vt(1)
               t14(2,nt4) = st(1)
               t24(1,nt4) = vt(2)
               t24(2,nt4) = st(2)
               t34(1,nt4) = vt(3)
               t34(2,nt4) = st(3)
            else if (tt .eq. 5) then
               nt5 = nt5 + 1
               if (ib .lt. ic) then
                  kt5(nt5) = pa//pb//pc//pd
               else if (ic .lt. ib) then
                  kt5(nt5) = pd//pc//pb//pa
               else if (ia .le. id) then
                  kt5(nt5) = pa//pb//pc//pd
               else if (id .lt. ia) then
                  kt5(nt5) = pd//pc//pb//pa
               end if
               call torphase (ft,vt,st)
               t15(1,nt5) = vt(1)
               t15(2,nt5) = st(1)
               t25(1,nt5) = vt(2)
               t25(2,nt5) = st(2)
               t35(1,nt5) = vt(3)
               t35(2,nt5) = st(3)
            end if
c
c     MMFF bond charge increment parameters
c
         else if (keyword(1:8) .eq. 'MMFFBCI ') then
            ia = 0
            ib = 0
            cg = 1000.0d0
            bt = 2
            string = record(next:240)
            read (string,*,err=640,end=640)  ia,ib,cg,bt
  640       continue
            if (ia .ne. 0) then
               if (bt .eq. 0) then
                  bci(ia,ib) = cg
                  bci(ib,ia) = -cg
               else if (bt .eq. 1) then
                  bci_1(ia,ib) = cg
                  bci_1(ib,ia) = -cg
               end if
            end if
c
c     MMFF partial bond charge increment parameters
c
         else if (keyword(1:9) .eq. 'MMFFPBCI ') then
            ia = 0
            string = record(next:240)
            read (string,*,err=650,end=650)  ia,cg,factor
  650       continue
            if (ia .ne. 0) then
               pbci(ia) = cg
               fcadj(ia) = factor
            end if
c
c     MMFF aromatic ion parameters
c
         else if (keyword(1:9) .eq. 'MMFFAROM ') then
            string = record(next:240)
            read (string,*,err=660,end=660)  ia,ib,ic,id,ie,if
  660       continue
            if (ie.eq.0 .and. id.eq.0) then
               mmffarom(ia,if) = ic
            else if (id .eq. 1) then
               mmffaromc(ia,if) = ic
            else if (ie .eq. 1) then
               mmffaroma(ia,if) = ic
            end if
         end if
      end do
      return
      end
