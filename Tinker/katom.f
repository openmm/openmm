c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine katom  --  atom type parameter assignment  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "katom" assigns an atom type definitions to each atom in
c     the structure and processes any new or changed values
c
c     literature reference:
c
c     K. A. Feenstra, B. Hess and H. J. C. Berendsen, "Improving
c     Efficiency of Large Time-Scale Molecular Dynamics Simulations
c     of Hydrogen-Rich Systems", Journal of Computational Chemistry,
c     8, 786-798 (1999)  [heavy hydrogen reweighting]
c
c
      subroutine katom
      use sizes
      use atomid
      use atoms
      use couple
      use inform
      use iounit
      use katoms
      use keys
      implicit none
      integer i,j,k
      integer next,nh
      integer cls,atn,lig
      real*8 wght
      real*8 hmax,hmass
      real*8 sum,dmass
      logical header,heavy
      character*3 symb
      character*20 keyword
      character*24 notice
      character*240 record
      character*240 string
c
c
c     process keywords containing atom type parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:5) .eq. 'ATOM ') then
            k = 0
            cls = 0
            symb = ' '
            notice = ' '
            atn = 0
            wght = 0.0d0
            lig = 0
            call getnumb (record,k,next)
            call getnumb (record,cls,next)
            if (cls .eq. 0)  cls = k
            atmcls(k) = cls
            call gettext (record,symb,next)
            call getstring (record,notice,next)
            string = record(next:240)
            read (string,*,err=40,end=40)  atn,wght,lig
            if (k.ge.1 .and. k.le.maxtyp) then
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,10)
   10             format (/,' Additional Atom Type Parameters :',
     &                    //,5x,'Type  Class  Symbol  Description',
     &                       15x,'Atomic',4x,'Mass',3x,'Valence',/)
               end if
               symbol(k) = symb
               describe(k) = notice
               atmnum(k) = atn
               weight(k) = wght
               ligand(k) = lig
               if (.not. silent) then
                  write (iout,20)  k,cls,symb,notice,atn,wght,lig
   20             format (2x,i6,1x,i6,5x,a3,3x,a24,i6,f11.3,i6)
               end if
            else if (k .ge. maxtyp) then
               write (iout,30)
   30          format (/,' KATOM   --  Too many Atom Types;',
     &                    ' Increase MAXTYP')
               abort = .true.
            end if
   40       continue
         end if
      end do
c
c     transfer atom type values to individual atoms
c
      do i = 1, n
         k = type(i)
         if (k .eq. 0) then
            class(i) = 0
            atomic(i) = 0
            mass(i) = 0.0d0
            valence(i) = 0
            story(i) = 'Undefined Atom Type     '
         else
            if (symbol(k) .ne. '   ')  name(i) = symbol(k)
            class(i) = atmcls(k)
            atomic(i) = atmnum(k)
            mass(i) = weight(k)
            valence(i) = ligand(k)
            story(i) = describe(k)
         end if
      end do
c
c     repartition hydrogen masses to use "heavy" hydrogens
c
      heavy = .false.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:15) .eq. 'HEAVY-HYDROGEN ') then
            heavy = .true.
         end if
      end do
      if (heavy) then
         hmax = 4.0d0
         do i = 1, n
            nh = 0
            sum = mass(i)
            do j = 1, n12(i)
               k = i12(j,i)
               if (atomic(k) .eq. 1) then
                  nh = nh + 1
                  sum = sum + mass(k)
               end if
            end do
            hmass = max(hmax,sum/dble(nh+1))
            do j = 1, n12(i)
               k = i12(j,i)
               if (atomic(k) .eq. 1) then
                  dmass = hmass - mass(k)
                  mass(k) = mass(k) + dmass
                  mass(i) = mass(i) - dmass
               end if
            end do
         end do
      end if
c
c     process keywords containing atom types for specific atoms
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:5) .eq. 'ATOM ') then
            k = 0
            symb = ' '
            notice = ' '
            atn = 0
            wght = 0.0d0
            lig = 0
            call getnumb (record,k,next)
            call getnumb (record,cls,next)
            call gettext (record,symb,next)
            call getstring (record,notice,next)
            string = record(next:240)
            read (string,*,err=70,end=70)  atn,wght,lig
            if (k.lt.0 .and. k.ge.-n) then
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,50)
   50             format (/,' Additional Atom Types for',
     &                       ' Specific Atoms :',
     &                    //,5x,'Atom  Class  Symbol  Description',
     &                       15x,'Atomic',4x,'Mass',3x,'Valence',/)
               end if
               k = -k
               if (cls .eq. 0)  cls = k
               class(k) = cls
               name(k) = symb
               story(k) = notice
               atomic(k) = atn
               mass(k) = wght
               valence(k) = lig
               if (.not. silent) then
                  write (iout,60)  k,cls,symb,notice,atn,wght,lig
   60             format (2x,i6,1x,i6,5x,a3,3x,a24,i6,f11.3,i6)
               end if
            end if
   70       continue
         end if
      end do
c
c     check for presence of undefined atom types or classes
c
      header = .true.
      do i = 1, n
         k = type(i)
         cls = class(i)
         if (k.lt.1 .or. k.gt.maxtyp
     &          .or. cls.lt.1 .or. cls.gt.maxclass) then
            abort = .true.
            if (header) then
               header = .false.
               write (iout,80)
   80          format (/,' Undefined Atom Types or Classes :',
     &                 //,' Type',10x,'Atom Number',5x,'Atom Type',
     &                    5x,'Atom Class',/)
            end if
            write (iout,90)  i,k,cls
   90       format (' Atom',12x,i5,10x,i5,10x,i5)
         end if
      end do
c
c     check the number of atoms attached to each atom
c
      header = .true.
      do i = 1, n
         if (n12(i) .ne. valence(i)) then
            if (header) then
               header = .false.
               write (iout,100)
  100          format (/,' Atoms with an Unusual Number of Attached',
     &                    ' Atoms :',
     &                 //,' Type',11x,'Atom Name',6x,'Atom Type',7x,
     &                    'Expected',4x,'Found',/)
            end if
            write (iout,110)  i,name(i),type(i),valence(i),n12(i)
  110       format (' Valence',5x,i7,'-',a3,8x,i5,10x,i5,5x,i5)
         end if
      end do
      return
      end
