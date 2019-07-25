c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kpolar  --  assign polarizability parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kpolar" assigns atomic dipole polarizabilities to the atoms
c     within the structure and processes any new or changed values
c
c     literature reference:
c
c     A. C. Simmonett, F. C. Pickard IV, J. W. Ponder and B. R. Brooks,
c     "An Empirical Extrapolation Scheme for Efficient Treatment of
c     Induced Dipoles", Journal of Chemical Physics, 145, 164101 (2016)
c     [OPT coefficients]
c
c
      subroutine kpolar
      use sizes
      use atoms
      use inform
      use iounit
      use keys
      use kpolr
      use mpole
      use neigh
      use polar
      use polpot
      use potent
      use usolve
      implicit none
      integer i,j,k,next
      integer nlist,npg
      integer pg(maxval)
      integer, allocatable :: list(:)
      real*8 pol,thl,dird,pena
      real*8 sixth
      logical header
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (list(n))
c
c     set defaults for numbers and lists of polarizable atoms
c
      nlist = 0
      do i = 1, n
         list(i) = 0
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(copt))  deallocate (copt)
      if (allocated(copm))  deallocate (copm)
      allocate (copt(0:maxopt))
      allocate (copm(0:maxopt))
c
c     set defaults for OPT induced dipole coefficients
c
      if (poltyp .eq. 'OPT')  poltyp = 'OPT3'
      do i = 0, maxopt
         copt(i) = 0.0d0
         copm(i) = 0.0d0
      end do
      if (poltyp .eq. 'OPT1') then
         copt(0) = 0.412d0
         copt(1) = 0.784d0
      else if (poltyp .eq. 'OPT2') then
         copt(0) = -0.115d0
         copt(1) = 0.568d0
         copt(2) = 0.608d0
      else if (poltyp .eq. 'OPT3') then
         copt(0) = -0.154d0
         copt(1) = 0.017d0
         copt(2) = 0.657d0
         copt(3) = 0.475d0
      else if (poltyp .eq. 'OPT4') then
         copt(0) = -0.041d0
         copt(1) = -0.176d0
         copt(2) = 0.169d0
         copt(3) = 0.663d0
         copt(4) = 0.374d0
      end if
      if (poltyp(1:3) .eq. 'OPT')  poltyp = 'OPT   '
c
c     get keywords containing polarization-related options
c
      do j = 1, nkey
         next = 1
         record = keyline(j)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:12) .eq. 'POLARIZABLE ') then
            read (string,*,err=10,end=10)  (list(i),i=nlist+1,n)
   10       continue
            do while (list(nlist+1) .ne. 0)
               nlist = nlist + 1
            end do
         else if (keyword(1:10) .eq. 'POLAR-OPT ') then
            do i = 0, maxopt
               copt(i) = 0.0d0
            end do
            read (string,*,err=20,end=20)  (copt(i),i=0,maxopt)
         end if
   20    continue
      end do
c
c     get maximum coefficient order for OPT induced dipoles
c
      coptmax = 0
      do i = 1, maxopt
         if (copt(i) .ne. 0.0d0)  coptmax = max(i,coptmax)
      end do
      do i = 0, coptmax
         do j = coptmax, i, -1
            copm(i) = copm(i) + copt(j)
         end do
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(polarity))  deallocate (polarity)
      if (allocated(penalpha))  deallocate (penalpha)
      if (allocated(thole))  deallocate (thole)
      if (allocated(dirdamp))  deallocate (dirdamp)
      if (allocated(pdamp))  deallocate (pdamp)
      if (allocated(udir))  deallocate (udir)
      if (allocated(udirp))  deallocate (udirp)
      if (allocated(uind))  deallocate (uind)
      if (allocated(uinp))  deallocate (uinp)
      if (allocated(douind))  deallocate (douind)
      allocate (polarity(n))
      allocate (thole(n))
      allocate (dirdamp(n))
      allocate (penalpha(n))
      allocate (pdamp(n))
      allocate (udir(3,n))
      allocate (udirp(3,n))
      allocate (uind(3,n))
      allocate (uinp(3,n))
      allocate (douind(n))
      if (allocated(uopt))  deallocate (uopt)
      if (allocated(uoptp))  deallocate (uoptp)
      if (allocated(fopt))  deallocate (fopt)
      if (allocated(foptp))  deallocate (foptp)
      if (poltyp .eq. 'OPT') then
         allocate (uopt(0:coptmax,3,n))
         allocate (uoptp(0:coptmax,3,n))
         allocate (fopt(0:coptmax,10,n))
         allocate (foptp(0:coptmax,10,n))
      end if
c
c     set the atoms allowed to have nonzero induced dipoles
c
      do i = 1, n
         douind(i) = .true.
      end do
      i = 1
      do while (list(i) .ne. 0)
         if (i .eq. 1) then
            do j = 1, n
               douind(j) = .false.
            end do
         end if
         if (list(i).gt.0 .and. list(i).le.n) then
            j = list(i)
            if (.not. douind(j)) then
               douind(j) = .true.
            end if
         else if (list(i).lt.0 .and. list(i).ge.-n) then
            do j = abs(list(i)), abs(list(i+1))
               if (.not. douind(j)) then
                  douind(j) = .true.
               end if
            end do
            i = i + 1
         end if
         i = i + 1
      end do
c
c     perform dynamic allocation of some local arrays
c
      deallocate (list)
c
c     process keywords containing polarizability parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:9) .eq. 'POLARIZE ') then
            k = 0
            pol = 0.0d0
            thl = -1.0d0
            adird = -1.0d0
            pena = -1.0d0
            do j = 1, maxval
               pg(j) = 0
            end do
            call getnumb (record,k,next)
            string = record(next:240)
            read (string,*,err=30,end=30)pol,thl,dird,(pg(j),j=1,maxval)
   30       continue
            if (k .gt. 0) then
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,40)
   40             format (/,' Additional Atomic Dipole',
     &                       ' Polarizability Parameters :',
     &                    //,5x,'Atom Type',11x,'Alpha',8x,
     &                       'Damp',5x,'Group Atom Types'/)
               end if
               if (k .le. maxtyp) then
                  polr(k) = pol
                  athl(k) = thl
                  adird(k) = dird
                  do j = 1, maxval
                     pgrp(j,k) = pg(j)
                     if (pg(j) .eq. 0) then
                        npg = j - 1
                        goto 50
                     end if
                  end do
   50             continue
                  if (.not. silent) then
                     write (iout,60)  k,pol,thl,dird,(pg(j),j=1,npg)
   60                format (4x,i6,10x,f10.3,2x,f10.3,2x,f10.3,7x,20i5)
                  end if
               else
                  write (iout,70)
   70             format (/,' KPOLAR  --  Too many Dipole',
     &                       ' Polarizability Parameters')
                  abort = .true.
               end if
            end if
         end if
      end do
c
c     find and store the atomic dipole polarizability parameters
c
      do i = 1, n
         polarity(i) = polr(type(i))
         thole(i) = athl(type(i))
         dirdamp(i) = adird(type(i))
         penalpha(i) = apena(type(i))
      end do
c
c     process keywords containing atom specific polarizabilities
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:9) .eq. 'POLARIZE ') then
            k = 0
            pol = 0.0d0
            thl = 0.0d0
            dird = 0.0d0
            call getnumb (record,k,next)
            if (k.lt.0 .and. k.ge.-n) then
               k = -k
               string = record(next:240)
               read (string,*,err=80,end=80)  pol,thl,dird
   80          continue
               if (header) then
                  header = .false.
                  write (iout,90)
   90             format (/,' Additional Dipole Polarizabilities',
     &                       ' for Specific Atoms :',
     &                    //,6x,'Atom',15x,'Alpha',8x,'Damp','DirDamp'/)
               end if
               if (.not. silent) then
                  write (iout,100)  k,pol,thl,dird
  100             format (4x,i6,10x,f10.3,2x,f10.3,f10.3)
               end if
               polarity(k) = pol
               thole(k) = thl
               dirdamp(k) = dird 
            end if
         end if
      end do
c
c     process keywords containing atom specific polarizabilities
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:3) .eq. 'CP ') then
            k = 0
            pena = 0.0d0
            call getnumb (record,k,next)
            if (k.lt.0 .and. k.ge.-n) then
               k = -k
               string = record(next:240)
               read (string,*,err=110,end=110)  pena
  110          continue
               if (header) then
                  header = .false.
                  write (iout,120)
  120             format (/,' Additional Charge Penetration',
     &                       ' for Specific Atoms :',
     &                    //,6x,'Atom',15x,'PenAlpha'/)
               end if
               if (.not. silent) then
                  write (iout,130)  k,pena
  130             format (4x,i6,f10.4)
               end if
               penalpha(k) = pena 
            end if
         end if
      end do
c
c     remove zero and undefined polarizable sites from the list
c
      npolar = 0
      if (use_polar) then
         npole = 0
         do i = 1, n
            if (polsiz(i).ne.0 .or. polarity(i).ne.0.0d0) then
               npole = npole + 1
               ipole(npole) = i
               pollist(i) = npole
               zaxis(npole) = zaxis(i)
               xaxis(npole) = xaxis(i)
               yaxis(npole) = yaxis(i)
               polaxe(npole) = polaxe(i)
               do k = 1, maxpole
                  pole(k,npole) = pole(k,i)
               end do
               if (polarity(i) .ne. 0.0d0)  npolar = npolar + 1
               polarity(npole) = polarity(i)
               thole(npole) = thole(i)
               dirdamp(npole) = dirdamp(i)
               penalpha(npole) = penalpha(i)
            end if
         end do
      end if
c
c     set the values used in the scaling of the polarizability
c
      sixth = 1.0d0 / 6.0d0
      do i = 1, npole
         if (thole(i) .eq. 0.0d0) then
            pdamp(i) = 0.0d0
         else
            pdamp(i) = polarity(i)**sixth
         end if
      end do
c
c     assign polarization group connectivity of each atom
c
      call polargrp
c
c     test multipoles at chiral sites and invert if necessary
c
      call chkpole
c
c     turn off polarizable multipole potential if it is not used
c
      if (npole .eq. 0)  use_mpole = .false.
      if (npolar .eq. 0)  use_polar = .false.
c
c     perform dynamic allocation of some global arrays
c
      if (use_polar) then
         if (allocated(mindex))  deallocate (mindex)
         if (allocated(minv))  deallocate (minv)
         allocate (mindex(npole))
         allocate (minv(3*maxulst*npole))
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine polargrp  --  polarization group connectivity  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "polargrp" generates members of the polarization group of
c     each atom and separate lists of the 1-2, 1-3 and 1-4 group
c     connectivities
c
c
      subroutine polargrp
      use sizes
      use atoms
      use couple
      use inform
      use iounit
      use kpolr
      use mpole
      use polgrp
      implicit none
      integer maxlist,maxkeep
      parameter (maxkeep=100)
      parameter (maxlist=1000)
      integer i,j,k,m
      integer it,jt
      integer jj,kk
      integer start,stop
      integer nlist,nkeep
      integer keep(maxkeep)
      integer list(maxlist)
      integer, allocatable :: mask(:)
      logical done
c
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(np11))  deallocate (np11)
      if (allocated(np12))  deallocate (np12)
      if (allocated(np13))  deallocate (np13)
      if (allocated(np14))  deallocate (np14)
      if (allocated(ip11))  deallocate (ip11)
      if (allocated(ip12))  deallocate (ip12)
      if (allocated(ip13))  deallocate (ip13)
      if (allocated(ip14))  deallocate (ip14)
      allocate (np11(n))
      allocate (np12(n))
      allocate (np13(n))
      allocate (np14(n))
      allocate (ip11(maxp11,n))
      allocate (ip12(maxp12,n))
      allocate (ip13(maxp13,n))
      allocate (ip14(maxp14,n))
c
c     find the directly connected group members for each atom
c
      do i = 1, n
         np11(i) = 1
         ip11(1,i) = i
         it = type(i)
         do j = 1, n12(i)
            jj = i12(j,i)
            jt = type(jj)
            do k = 1, maxval
               kk = pgrp(k,it)
               if (kk .eq. 0)  goto 20
               if (pgrp(k,it) .eq. jt) then
                  np11(i) = np11(i) + 1
                  if (np11(i) .le. maxp11) then
                     ip11(np11(i),i) = jj
                  else
                     write (iout,10)
   10                format (/,' POLARGRP  --  Too many Atoms',
     &                          ' in Polarization Group')
                     abort = .true.
                     goto 30
                  end if
               end if
            end do
   20       continue
         end do
      end do
   30 continue
c
c     make sure all connected group members are bidirectional
c
      do i = 1, n
         do j = 1, np11(i)
            k = ip11(j,i)
            do m = 1, np11(k)
               if (ip11(m,k) .eq. i)  goto 50
            end do
            write (iout,40)  min(i,k),max(i,k)
   40       format (/,' POLARGRP  --  Check Polarization Groups for',
     &                 ' Atoms',i9,' and',i9)
            abort = .true.
   50       continue
         end do
      end do
      if (abort)  call fatal
c
c     perform dynamic allocation of some local arrays
c
      allocate (mask(n))
c
c     find any other group members for each atom in turn
c
      do i = 1, n
         mask(i) = 0
      end do
      do i = 1, n
         done = .false.
         start = 1
         stop = np11(i)
         do j = start, stop
            jj = ip11(j,i)
            if (jj .lt. i) then
               done = .true.
               np11(i) = np11(jj)
               do k = 1, np11(i)
                  ip11(k,i) = ip11(k,jj)
               end do
            else
               mask(jj) = i
            end if
         end do
         do while (.not. done)
            done = .true.
            do j = start, stop
               jj = ip11(j,i)
               do k = 1, np11(jj)
                  kk = ip11(k,jj)
                  if (mask(kk) .ne. i) then
                     np11(i) = np11(i) + 1
                     if (np11(i) .le. maxp11) then
                        ip11(np11(i),i) = kk
                     else
                        write (iout,60)
   60                   format (/,' POLARGRP  --  Too many Atoms',
     &                             ' in Polarization Group')
                        abort = .true.
                        goto 70
                     end if
                     mask(kk) = i
                  end if
               end do
            end do
            if (np11(i) .ne. stop) then
               done = .false.
               start = stop + 1
               stop = np11(i)
            end if
         end do
         call sort (np11(i),ip11(1,i))
      end do
   70 continue
c
c     loop over atoms finding all the 1-2 group relationships
c
      do i = 1, n
         mask(i) = 0
      end do
      do i = 1, n
         do j = 1, np11(i)
            jj = ip11(j,i)
            mask(jj) = i
         end do
         nkeep = 0
         do j = 1, np11(i)
            jj = ip11(j,i)
            do k = 1, n12(jj)
               kk = i12(k,jj)
               if (mask(kk) .ne. i) then
                  nkeep = nkeep + 1
                  keep(nkeep) = kk
               end if
            end do
         end do
         nlist = 0
         do j = 1, nkeep
            jj = keep(j)
            do k = 1, np11(jj)
               kk = ip11(k,jj)
               nlist = nlist + 1
               list(nlist) = kk
            end do
         end do
         call sort8 (nlist,list)
         if (nlist .le. maxp12) then
            np12(i) = nlist
            do j = 1, nlist
               ip12(j,i) = list(j)
            end do
         else
            write (iout,80)
   80       format (/,' POLARGRP  --  Too many Atoms',
     &                 ' in 1-2 Polarization Group')
            abort = .true.
            goto 90
         end if
      end do
   90 continue
c
c     loop over atoms finding all the 1-3 group relationships
c
      do i = 1, n
         mask(i) = 0
      end do
      do i = 1, n
         do j = 1, np11(i)
            jj = ip11(j,i)
            mask(jj) = i
         end do
         do j = 1, np12(i)
            jj = ip12(j,i)
            mask(jj) = i
         end do
         nlist = 0
         do j = 1, np12(i)
            jj = ip12(j,i)
            do k = 1, np12(jj)
               kk = ip12(k,jj)
               if (mask(kk) .ne. i) then
                  nlist = nlist + 1
                  list(nlist) = kk
               end if
            end do
         end do
         call sort8 (nlist,list)
         if (nlist .le. maxp13) then
            np13(i) = nlist
            do j = 1, nlist
               ip13(j,i) = list(j)
            end do
         else
            write (iout,100)
  100       format (/,' POLARGRP  --  Too many Atoms',
     &                 ' in 1-3 Polarization Group')
            abort = .true.
            goto 110
         end if
      end do
  110 continue
c
c     loop over atoms finding all the 1-4 group relationships
c
      do i = 1, n
         mask(i) = 0
      end do
      do i = 1, n
         do j = 1, np11(i)
            jj = ip11(j,i)
            mask(jj) = i
         end do
         do j = 1, np12(i)
            jj = ip12(j,i)
            mask(jj) = i
         end do
         do j = 1, np13(i)
            jj = ip13(j,i)
            mask(jj) = i
         end do
         nlist = 0
         do j = 1, np13(i)
            jj = ip13(j,i)
            do k = 1, np12(jj)
               kk = ip12(k,jj)
               if (mask(kk) .ne. i) then
                  nlist = nlist + 1
                  list(nlist) = kk
               end if
            end do
         end do
         call sort8 (nlist,list)
         if (nlist .le. maxp14) then
            np14(i) = nlist
            do j = 1, nlist
               ip14(j,i) = list(j)
            end do
         else
            write (iout,120)
  120       format (/,' POLARGRP  --  Too many Atoms',
     &                 ' in 1-4 Polarization Group')
            abort = .true.
            goto 130
         end if
      end do
  130 continue
c
c     perform deallocation of some local arrays
c
      deallocate (mask)
      return
      end
