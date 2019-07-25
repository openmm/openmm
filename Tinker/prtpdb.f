c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine prtpdb  --  output of Protein Data Bank file  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "prtpdb" writes out a set of Protein Data Bank coordinates
c     to an external disk file
c
c
      subroutine prtpdb (ipdb)
      use sizes
      use files
      use pdb
      use sequen
      use titles
      implicit none
      integer i,k,ipdb
      integer start,stop
      integer resmax,resnumb
      integer, allocatable :: resid(:)
      real*8 crdmin,crdmax
      logical opened
      logical rename
      logical reformat
      character*1 chnname
      character*1, allocatable :: chain(:)
      character*2 atmc,resc
      character*3 resname
      character*6 crdc
      character*38 fstr
      character*240 pdbfile
c
c
c     set flags for residue naming and large value formatting
c
      rename = .false.
      reformat = .true.
c
c     open the output unit if not already done
c
      inquire (unit=ipdb,opened=opened)
      if (.not. opened) then
         pdbfile = filename(1:leng)//'.pdb'
         call version (pdbfile,'new')
         open (unit=ipdb,file=pdbfile,status='new')
      end if
c
c     write out the header lines and the title
c
      if (ltitle .eq. 0) then
         fstr = '(''HEADER'',/,''COMPND'',/,''SOURCE'')'
         write (ipdb,fstr(1:32))
      else
         fstr = '(''HEADER'',4x,a,/,''COMPND'',/,''SOURCE'')'
         write (ipdb,fstr(1:37))  title(1:ltitle)
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (resid(maxres))
      allocate (chain(maxres))
c
c     find the chain name and chain position for each residue
c
      do i = 1, nchain
         start = ichain(1,i)
         stop = ichain(2,i)
         do k = start, stop
            resid(k) = k - start + 1
            chain(k) = chnnam(i)
         end do
      end do
c
c     change some TINKER residue names to match PDB standards
c
      if (rename) then
         do i = 1, npdb
            if (pdbres(i) .eq. 'CYX')  pdbres(i) = 'CYS'
            if (pdbres(i) .eq. 'CYD')  pdbres(i) = 'CYS'
            if (pdbres(i) .eq. 'TYD')  pdbres(i) = 'TYR'
            if (pdbres(i) .eq. 'HID')  pdbres(i) = 'HIS'
            if (pdbres(i) .eq. 'HIE')  pdbres(i) = 'HIS'
            if (pdbres(i) .eq. 'HIP')  pdbres(i) = 'HIS'
            if (pdbres(i) .eq. 'ASH')  pdbres(i) = 'ASP'
            if (pdbres(i) .eq. 'GLH')  pdbres(i) = 'GLU'
            if (pdbres(i) .eq. 'LYD')  pdbres(i) = 'LYS'
         end do
      end if
c
c     set formatting to match the PDB fixed format standard
c
      atmc = 'i5'
      resc = 'i4'
      crdc = '3f8.3 '
c
c     check for large values requiring extended formatting
c
      if (reformat) then
         resmax = 0
         crdmin = 0.0d0
         crdmax = 0.0d0
         do i = 1, npdb
            if (pdbtyp(i) .eq. 'ATOM  ') then
               resmax = max(resmax,resid(resnum(i)))
            else
               resmax = max(resmax,resnum(i))
            end if
            crdmin = min(crdmin,xpdb(i),ypdb(i),zpdb(i))
            crdmax = max(crdmax,xpdb(i),ypdb(i),zpdb(i))
         end do
         if (npdb .ge. 100000)  atmc = 'i6'
         if (resmax .ge. 10000)  resc = 'i5'
         if (resmax .ge. 100000)  resc = 'i6'
         if (crdmin .le. -100.0d0)  crdc = '3f9.3 '
         if (crdmax .ge. 1000.0d0)  crdc = '3f9.3 '
         if (crdmin .le. -1000.0d0)  crdc = '3f10.3'
         if (crdmax .ge. 10000.0d0)  crdc = '3f10.3'
      end if
c
c     write info and coordinates for each PDB atom
c
      fstr = '(a6,'//atmc//',1x,a4,1x,a3,1x,a1,'//resc//
     &          ',4x,'//crdc//')'
      do i = 1, npdb
         resname = pdbres(i)
         if (resname(2:3) .eq. '  ')  resname = '  '//resname(1:1)
         if (resname(3:3) .eq. ' ')  resname = ' '//resname(1:2)
         if (pdbtyp(i) .eq. 'ATOM  ') then
            resnumb = resid(resnum(i))
            chnname = chain(resnum(i))
         else
            resnumb = resnum(i)
            chnname = ' '
         end if
         write (ipdb,fstr)  pdbtyp(i),i,pdbatm(i),resname,chnname,
     &                      resnumb,xpdb(i),ypdb(i),zpdb(i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (resid)
      deallocate (chain)
c
c     check for large values requiring extended formatting
c
      if (reformat) then
         if (npdb .ge. 100000)  atmc = 'i7'
         if (npdb .ge. 10000)  atmc = 'i6'
      end if
c
c     write any connectivity records for PDB atoms
c
      fstr = '(''CONECT'',9'//atmc//')'
      do i = 1, npdb
         if (npdb12(i) .ne. 0) then
            write (ipdb,fstr(1:14))  i,(ipdb12(k,i),k=1,npdb12(i))
         end if
      end do
      fstr = '(''END'')'
      write (ipdb,fstr(1:7))
c
c     close the output unit if opened by this routine
c
c     if (.not. opened)  close (unit=ipdb)
      return
      end
