c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  program pdbxyz  --  Protein Data Bank to XYZ coordinates  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "pdbxyz" takes as input a Protein Data Bank file and then
c     converts to and writes out a Cartesian coordinates file and,
c     for biopolymers, a sequence file
c
c
      program pdbxyz
      use sizes
      use atomid
      use atoms
      use couple
      use files
      use inform
      use katoms
      use pdb
      use resdue
      use sequen
      use titles
      implicit none
      integer i,j,it,next
      integer ipdb,ixyz,iseq
      integer last,pdbleng
      integer freeunit
      integer, allocatable :: row(:)
      real*8 xi,yi,zi,rij
      real*8 rcut,rmax(0:9)
      logical biopoly
      logical clash
      character*1 letter
      character*3 resname
      character*3 reslast
      character*240 pdbfile
      character*240 xyzfile
      character*240 seqfile
      character*240 pdbtitle
c
c
c     get the Protein Data Bank file and a parameter set
c
      call initial
      call getpdb
      call field
      call unitcell
c
c     save the title line from the PDB file for later use
c
      pdbleng = ltitle
      pdbtitle = title(1:ltitle)
c
c     decide whether the system has only biopolymers and water
c
      biopoly = .false.
      reslast = '***'
      do i = 1, npdb
         if (pdbtyp(i) .eq. 'ATOM  ') then
            resname = pdbres(i)
            if (resname .ne. reslast) then
               reslast = resname
               do j = 1, maxamino
                  if (resname .eq. amino(j)) then
                     biopoly = .true.
                     goto 10
                  end if
               end do
               do j = 1, maxnuc
                  if (resname .eq. nuclz(j)) then
                     biopoly = .true.
                     goto 10
                  end if
               end do
               biopoly = .false.
               goto 20
   10          continue
            end if
         else if (pdbtyp(i) .eq. 'HETATM') then
            resname = pdbres(i)
            if (resname .ne. reslast) then
               reslast = resname
               if (resname.eq.'HOH' .or. resname.eq.'NA ' .or.
     &             resname.eq.'K  ' .or. resname.eq.'MG ' .or.
     &             resname.eq.'CA ' .or. resname.eq.'CL ') then
                  pdbtyp(i) = 'HETATM'
               end if
            end if
         end if
      end do
   20 continue
c
c     open the TINKER coordinates file to be used for output
c
      ixyz = freeunit ()
      xyzfile = filename(1:leng)//'.xyz'
      call version (xyzfile,'new')
      open (unit=ixyz,file=xyzfile,status='new')
c
c     reopen the PDB file and read the first coordinate set
c
      ipdb = freeunit ()
      pdbfile = filename
      call suffix (pdbfile,'pdb','old')
      open (unit=ipdb,file=pdbfile,status ='old')
      rewind (unit=ipdb)
      call readpdb (ipdb)
c
c     use special translation mechanisms used for biopolymers
c
      do while (.not. abort)
         if (biopoly) then
            n = 0
            do i = 1, nchain
               if (chntyp(i) .eq. 'PEPTIDE')  call ribosome (i)
               if (chntyp(i) .eq. 'NUCLEIC')  call ligase (i)
            end do
            call hetatom
            last = n
            do i = last, 1, -1
               if (type(i) .eq. 0)  call delete (i)
            end do
c
c     get general atom properties for non-biopolymer structures
c
         else
            n = npdb
            do i = 1, n
               x(i) = xpdb(i)
               y(i) = ypdb(i)
               z(i) = zpdb(i)
               name(i) = pdbatm(i)(2:4)
               n12(i) = 0
               next = 1
               call getnumb (pdbres(i),type(i),next)
            end do
c
c     perform dynamic allocation of some local arrays
c
            allocate (row(n))
c
c     set atom size classification from periodic table row
c
            do i = 1, n
               it = type(i)
               if (it .eq. 0) then
                  letter = name(i)(1:1)
                  call upcase (letter)
                  if (letter .eq. 'H') then
                     row(i) = 1
                  else if (letter .eq. 'C') then
                     row(i) = 2
                  else if (letter .eq. 'N') then
                     row(i) = 2
                  else if (letter .eq. 'O') then
                     row(i) = 2
                  else if (letter .eq. 'P') then
                     row(i) = 3
                  else if (letter .eq. 'S') then
                     row(i) = 3
                  else
                     row(i) = 0
                  end if
               else if (ligand(it) .eq. 0) then
                  row(i) = 0
               else if (atmnum(it) .le. 2) then
                  row(i) = 1
               else if (atmnum(it) .le. 10) then
                  row(i) = 2
               else
                  row(i) = 3
               end if
            end do
c
c     set the maximum bonded distance between atom type pairs
c
            rmax(0) = -1.0d0
            rmax(1) = -1.0d0
            rmax(2) = 1.3d0
            rmax(3) = 1.55d0
            rmax(4) = 1.75d0
            rmax(6) = 2.0d0
            rmax(9) = 2.2d0
c
c     find and connect atom pairs within bonding distance
c
            do i = 1, n-1
               xi = x(i)
               yi = y(i)
               zi = z(i)
               do j = i+1, n
                  rcut = rmax(row(i)*row(j))**2
                  rij = (xi-x(j))**2 + (yi-y(j))**2 + (zi-z(j))**2
                  if (rij .le. rcut) then
                     n12(i) = n12(i) + 1
                     i12(n12(i),i) = j
                     n12(j) = n12(j) + 1
                     i12(n12(j),j) = i
                  end if
               end do
            end do
c
c     perform deallocation of some local arrays
c
            deallocate (row)
         end if
c
c     sort the attached atom lists into ascending order
c
         do i = 1, n
            call sort (n12(i),i12(1,i))
         end do
c
c     check for atom pairs with identical coordinates
c
         clash = .false.
         call chkxyz (clash)
c
c     write the TINKER coordinates and reset the connectivities
c
         ltitle = pdbleng
         title = pdbtitle(1:ltitle)
         call prtxyz (ixyz)
         do i = 1, n
            n12(i) = 0
         end do
c
c     read the next coordinate set from Protein Data Bank file
c
         call readpdb (ipdb)
      end do
c
c     write a sequence file for proteins and nucleic acids
c
      if (biopoly) then
         iseq = freeunit ()
         seqfile = filename(1:leng)//'.seq'
         call version (seqfile,'new')
         open (unit=iseq,file=seqfile,status='new')
         call prtseq (iseq)
         close (unit=iseq)
      end if
c
c     perform any final tasks before program exit
c
      close (unit=ixyz)
      call final
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine ribosome  --  coordinates from PDB polypeptide  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ribosome" translates a polypeptide structure in Protein Data
c     Bank format to a Cartesian coordinate file and sequence file
c
c
      subroutine ribosome (ichn)
      use sizes
      use atoms
      use fields
      use files
      use inform
      use iounit
      use pdb
      use resdue
      use sequen
      implicit none
      integer i,j,k,m
      integer ichn,ityp
      integer jres,kres
      integer start,stop
      integer cyxtyp
      integer ncys,ndisulf
      integer, allocatable :: ni(:)
      integer, allocatable :: cai(:)
      integer, allocatable :: ci(:)
      integer, allocatable :: oi(:)
      integer, allocatable :: si(:)
      integer, allocatable :: icys(:)
      integer, allocatable :: idisulf(:,:)
      real*8 xr,yr,zr,r
      real*8, allocatable :: xcys(:)
      real*8, allocatable :: ycys(:)
      real*8, allocatable :: zcys(:)
      logical newchain
      logical midchain
      logical endchain
      logical cyclic
      character*3 resname
      character*4 atmname
      save si
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (ni(nres))
      allocate (cai(nres))
      allocate (ci(nres))
      allocate (oi(nres))
      allocate (si(nres))
      allocate (icys(nres))
      allocate (idisulf(2,nres))
      allocate (xcys(nres))
      allocate (ycys(nres))
      allocate (zcys(nres))
c
c     set the next atom and the residue range of the chain
c
      n = n + 1
      jres = ichain(1,ichn)
      kres = ichain(2,ichn)
      do i = jres, kres
         ni(i) = 0
         cai(i) = 0
         ci(i) = 0
         oi(i) = 0
      end do
c
c     check for the presence of a cyclic polypeptide chain
c
      cyclic = .false.
      start = resatm(1,jres)
      stop = resatm(2,jres)
      call findatm (' N  ',start,stop,j)
      start = resatm(1,kres)
      stop = resatm(2,kres)
      call findatm (' C  ',start,stop,k)
      if (j.ne.0 .and. k.ne.0) then
         xr = xpdb(k) - xpdb(j)
         yr = ypdb(k) - ypdb(j)
         zr = zpdb(k) - zpdb(j)
         r = sqrt(xr*xr + yr*yr + zr*zr)
         if (r .le. 3.0d0) then
            cyclic = .true.
            ni(jres) = j
            ci(kres) = k
         end if
      end if
c
c     search for any potential cystine disulfide bonds
c
      do i = 1, maxamino
         if (amino(i) .eq. 'CYX')  cyxtyp = i
      end do
      ncys = 0
      do i = 1, nres
         start = resatm(1,i)
         resname = pdbres(start)
         if (resname.eq.'CYS' .or. resname.eq.'CYX') then
            stop = resatm(2,i)
            call findatm (' SG ',start,stop,k)
            ncys = ncys + 1
            icys(ncys) = i
            xcys(ncys) = xpdb(k)
            ycys(ncys) = ypdb(k)
            zcys(ncys) = zpdb(k)
         end if
      end do
      ndisulf = 0
      do i = 1, ncys-1
         do k = i+1, ncys
            xr = xcys(k) - xcys(i)
            yr = ycys(k) - ycys(i)
            zr = zcys(k) - zcys(i)
            r = sqrt(xr*xr + yr*yr + zr*zr)
            if (r .le. 3.0d0) then
               ndisulf = ndisulf + 1
               idisulf(1,ndisulf) = min(icys(i),icys(k))
               idisulf(2,ndisulf) = max(icys(i),icys(k))
            end if
         end do
      end do
      do i = 1, ndisulf
         j = idisulf(1,i)
         k = idisulf(2,i)
         seqtyp(j) = cyxtyp
         seqtyp(k) = cyxtyp
         seq(j) = 'CYX'
         seq(k) = 'CYX'
         start = resatm(1,j)
         stop = resatm(2,j)
         do m = start, stop
            pdbres(m) = 'CYX'
         end do
         start = resatm(1,k)
         stop = resatm(2,k)
         do m = start, stop
            pdbres(m) = 'CYX'
         end do
      end do
c
c     locate and assign the atoms that make up each residue
c
      do i = jres, kres
         ityp = seqtyp(i)
         start = resatm(1,i)
         stop = resatm(2,i)
         resname = seq(i)
c
c     check that the maximum allowed atoms is not exceeded
c
         if (n+25 .gt. maxatm) then
            write (iout,10)  maxatm
   10       format (/,' RIBOSOME  --  The Maximum of',i8,' Atoms',
     &                 ' has been Exceeded')
            call fatal
         end if
c
c     test location of residue within the current chain
c
         newchain = .false.
         midchain = .false.
         endchain = .false.
         if (i .eq. jres)  newchain = .true.
         if (i .eq. kres)  endchain = .true.
         if (.not.newchain .and. .not.endchain)  midchain = .true.
c
c     build the amide nitrogen of the current residue
c
         atmname = ' N  '
         if (resname .eq. 'COH')  atmname = ' OH '
         call findatm (atmname,start,stop,k)
         if (k .ne. 0)  ni(i) = n
         if (midchain) then
            j = ntyp(ityp)
            call oldatm (k,j,ci(i-1),i)
         else if (newchain) then
            if (cyclic) then
               j = ntyp(ityp)
            else
               j = nntyp(ityp)
            end if
            call oldatm (k,j,0,i)
         else if (endchain) then
            if (cyclic) then
               j = ntyp(ityp)
            else
               j = nctyp(ityp)
            end if
            call oldatm (k,j,ci(i-1),i)
         end if
c
c     build the alpha carbon of the current residue
c
         atmname = ' CA '
         if (resname .eq. 'ACE')  atmname = ' CH3'
         if (resname .eq. 'NME')  atmname = ' CH3'
         call findatm (atmname,start,stop,k)
         if (k .ne. 0)  cai(i) = n
         if (midchain .or. cyclic .or. nres.eq.1) then
            j = catyp(ityp)
            call oldatm (k,j,ni(i),i)
         else if (newchain) then
            j = cantyp(ityp)
            call oldatm (k,j,ni(i),i)
         else if (endchain) then
            j = cactyp(ityp)
            call oldatm (k,j,ni(i),i)
         end if
c
c     build the carbonyl carbon of the current residue
c
         call findatm (' C  ',start,stop,k)
         if (k .ne. 0)  ci(i) = n
         if (midchain .or. cyclic) then
            j = ctyp(ityp)
            call oldatm (k,j,cai(i),i)
         else if (newchain) then
            j = cntyp(ityp)
            call oldatm (k,j,cai(i),i)
         else if (endchain) then
            j = cctyp(ityp)
            if (resname .eq. 'COH') then
               type(ci(i-1)) = biotyp(j)
            else
               call oldatm (k,j,cai(i),i)
            end if
         end if
c
c     build the carbonyl oxygen of the current residue
c
         call findatm (' O  ',start,stop,k)
         if (k .ne. 0)  oi(i) = n
         if (midchain .or. cyclic) then
            j = otyp(ityp)
            call oldatm (k,j,ci(i),i)
         else if (newchain) then
            j = ontyp(ityp)
            call oldatm (k,j,ci(i),i)
         else if (endchain) then
            j = octyp(ityp)
            if (resname .eq. 'COH') then
               type(oi(i-1)) = biotyp(j)
            else
               call oldatm (k,j,ci(i),i)
            end if
         end if
c
c     build the amide hydrogens of the current residue
c
         if (midchain .or. (endchain.and.cyclic)) then
            j = hntyp(ityp)
            call findatm (' H  ',start,stop,k)
            call newatm (k,j,ni(i),1.01d0,ci(i-1),119.0d0,
     &                      cai(i),119.0d0,1)
         else if (newchain .and. cyclic) then
            j = hntyp(ityp)
            call findatm (' H  ',start,stop,k)
            call newatm (k,j,ni(i),1.01d0,ci(kres),119.0d0,
     &                      cai(i),119.0d0,1)
         else if (newchain) then
            j = hnntyp(ityp)
            if (resname .eq. 'PRO') then
               call findatm (' H2 ',start,stop,k)
               call newatm (k,j,ni(i),1.01d0,cai(i),109.5d0,
     &                         ci(i),0.0d0,0)
               call findatm (' H3 ',start,stop,k)
               call newatm (k,j,ni(i),1.01d0,cai(i),109.5d0,
     &                         ci(i),-120.0d0,0)
            else if (resname .eq. 'PCA') then
               call findatm (' H  ',start,stop,k)
               call newatm (k,j,ni(i),1.01d0,cai(i),109.5d0,
     &                         ci(i),-60.0d0,0)
            else
               call findatm (' H1 ',start,stop,k)
               call newatm (k,j,ni(i),1.01d0,cai(i),109.5d0,
     &                         ci(i),180.0d0,0)
               call findatm (' H2 ',start,stop,k)
               call newatm (k,j,ni(i),1.01d0,cai(i),109.5d0,
     &                         ci(i),60.0d0,0)
               call findatm (' H3 ',start,stop,k)
               call newatm (k,j,ni(i),1.01d0,cai(i),109.5d0,
     &                         ci(i),-60.0d0,0)
            end if
         else if (endchain) then
            j = hnctyp(ityp)
            if (resname .eq. 'COH') then
               call findatm (' HO ',start,stop,k)
               call newatm (k,j,ni(i),0.98d0,ci(i-1),108.7d0,
     &                         cai(i-1),180.0d0,0)
            else if (resname .eq. 'NH2') then
               call findatm (' H1 ',start,stop,k)
               call newatm (k,j,ni(i),1.01d0,ci(i-1),120.9d0,
     &                         cai(i-1),0.0d0,0)
               call findatm (' H2 ',start,stop,k)
               call newatm (k,j,ni(i),1.01d0,ci(i-1),120.3d0,
     &                         cai(i-1),180.0d0,0)
            else if (resname .eq. 'NME') then
               call findatm (' H  ',start,stop,k)
               call newatm (k,j,ni(i),1.01d0,ci(i-1),119.0d0,
     &                         cai(i),119.0d0,1)
            else
               call findatm (' H  ',start,stop,k)
               call newatm (k,j,ni(i),1.01d0,ci(i-1),119.0d0,
     &                         cai(i),119.0d0,1)
            end if
         end if
c
c     build the alpha hydrogen of the current residue
c
         if (resname .eq. 'GLY') then
            call findatm (' HA2',start,stop,k)
         else
            call findatm (' HA ',start,stop,k)
         end if
         if (midchain .or. cyclic) then
            j = hatyp(ityp)
            call newatm (k,j,cai(i),1.10d0,ni(i),109.5d0,
     &                      ci(i),109.5d0,-1)
         else if (newchain) then
            j = hantyp(ityp)
            if (resname .eq. 'FOR') then
               call findatm (' H  ',start,stop,k)
               call newatm (k,j,ci(i),1.12d0,oi(i),0.0d0,0,0.0d0,0)
            else if (resname .eq. 'ACE') then
               call findatm (' H1 ',start,stop,k)
               call newatm (k,j,cai(i),1.10d0,ci(i),109.5d0,
     &                         oi(i),180.0d0,0)
               call findatm (' H2 ',start,stop,k)
               call newatm (k,j,cai(i),1.10d0,ci(i),109.5d0,
     &                         oi(i),60.0d0,0)
               call findatm (' H3 ',start,stop,k)
               call newatm (k,j,cai(i),1.10d0,ci(i),109.5d0,
     &                         oi(i),-60.0d0,0)
            else
               call newatm (k,j,cai(i),1.10d0,ni(i),109.5d0,
     &                         ci(i),109.5d0,-1)
            end if
         else if (endchain) then
            j = hactyp(ityp)
            if (resname .eq. 'NME') then
               call findatm (' H1 ',start,stop,k)
               call newatm (k,j,cai(i),1.10d0,ni(i),109.5d0,
     &                         ci(i-1),180.0d0,0)
               call findatm (' H2 ',start,stop,k)
               call newatm (k,j,cai(i),1.10d0,ni(i),109.5d0,
     &                         ci(i-1),60.0d0,0)
               call findatm (' H3 ',start,stop,k)
               call newatm (k,j,cai(i),1.10d0,ni(i),109.5d0,
     &                         ci(i-1),-60.0d0,0)
            else
               call newatm (k,j,cai(i),1.10d0,ni(i),109.5d0,
     &                         ci(i),109.5d0,-1)
            end if
         end if
c
c     build the side chain atoms of the current residue
c
         call addside (resname,i,start,stop,cai(i),ni(i),ci(i),si(i))
c
c     build the terminal oxygen at the end of a peptide chain
c
         if (endchain .and. .not.cyclic .and. resname.ne.'COH') then
            call findatm (' OXT',start,stop,k)
            if (k .eq. 0)  call findatm (' OT2',start,stop,k)
            j = octyp(ityp)
            call newatm (k,j,ci(i),1.25d0,cai(i),117.0d0,
     &                      oi(i),126.0d0,1)
         end if
      end do
c
c     connect the terminal residues if the chain is cyclic
c
      if (cyclic) then
         call addbond (ni(jres),ci(kres))
         if (verbose) then
            write (iout,20)  jres,kres
   20       format (/,' Peptide Cyclization between Residues :  ',2i5)
         end if
      end if
c
c     connect the sulfur atoms involved in disulfide bonds
c
      do i = 1, ndisulf
         j = idisulf(1,i)
         k = idisulf(2,i)
         if (k.ge.ichain(1,ichn) .and. k.le.ichain(2,ichn)) then
            call addbond (si(j),si(k))
            if (verbose) then
               write (iout,30)  j,k
   30          format (/,' Disulfide Bond between Residues :  ',2i5)
            end if
         end if
      end do
c
c     total number of atoms is one less than the current atom
c
      n = n - 1
c
c     perform deallocation of some local arrays
c
      deallocate (ni)
      deallocate (cai)
      deallocate (ci)
      deallocate (oi)
      deallocate (si)
      deallocate (icys)
      deallocate (idisulf)
      deallocate (xcys)
      deallocate (ycys)
      deallocate (zcys)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine addside  --  build the amino acid side chains  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "addside" builds the Cartesian coordinates for a single amino
c     acid side chain; coordinates are read from the Protein Data
c     Bank file or found from internal coordinates, then atom types
c     are assigned and connectivity data generated
c
c     note biotypes of CD and HD atoms for N-terminal proline are
c     set as absolute values, not relative to the CB atom
c
c
      subroutine addside (resname,ires,start,stop,cai,ni,ci,si)
      use sizes
      use atoms
      use resdue
      use sequen
      implicit none
      integer i,k,ires
      integer start,stop
      integer cai,ni,ci,si
      character*3 resname
c
c
c     zero out disulfide and set CB atom as reference site
c
      si = 0
      k = cbtyp(seqtyp(ires))
c
c     glycine residue  (GLY)
c
      if (resname .eq. 'GLY') then
         call findatm (' HA3',start,stop,i)
         k = hatyp(seqtyp(ires))
         if (ires .eq. 1)  k = hantyp(seqtyp(ires))
         if (ires .eq. nseq)  k = hactyp(seqtyp(ires))
         call newatm (i,k,cai,1.10d0,ni,109.5d0,ci,109.5d0,1)
c
c     alanine residue  (ALA)
c
      else if (resname .eq. 'ALA') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,k,cai,ires)
         call findatm (' HB1',start,stop,i)
         call newatm (i,k+1,n-1,1.10d0,cai,110.2d0,ni,180.0d0,0)
         call findatm (' HB2',start,stop,i)
         call newatm (i,k+1,n-2,1.10d0,cai,110.2d0,ni,60.0d0,0)
         call findatm (' HB3',start,stop,i)
         call newatm (i,k+1,n-3,1.10d0,cai,110.2d0,ni,-60.0d0,0)
c
c     valine residue  (VAL)
c
      else if (resname .eq. 'VAL') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,k,cai,ires)
         call findatm (' CG1',start,stop,i)
         call oldatm (i,k+2,n-1,ires)
         call findatm (' CG2',start,stop,i)
         call oldatm (i,k+4,n-2,ires)
         call findatm (' HB ',start,stop,i)
         call newatm (i,k+1,n-3,1.10d0,cai,107.0d0,n-2,108.2d0,1)
         call findatm ('HG11',start,stop,i)
         call newatm (i,k+3,n-3,1.10d0,n-4,111.6d0,cai,180.0d0,0)
         call findatm ('HG12',start,stop,i)
         call newatm (i,k+3,n-4,1.10d0,n-5,111.6d0,cai,60.0d0,0)
         call findatm ('HG13',start,stop,i)
         call newatm (i,k+3,n-5,1.10d0,n-6,111.6d0,cai,-60.0d0,0)
         call findatm ('HG21',start,stop,i)
         call newatm (i,k+5,n-5,1.10d0,n-7,111.6d0,cai,180.0d0,0)
         call findatm ('HG22',start,stop,i)
         call newatm (i,k+5,n-6,1.10d0,n-8,111.6d0,cai,60.0d0,0)
         call findatm ('HG23',start,stop,i)
         call newatm (i,k+5,n-7,1.10d0,n-9,111.6d0,cai,-60.0d0,0)
c
c     leucine residue  (LEU)
c
      else if (resname .eq. 'LEU') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,k,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,k+2,n-1,ires)
         call findatm (' CD1',start,stop,i)
         call oldatm (i,k+4,n-1,ires)
         call findatm (' CD2',start,stop,i)
         call oldatm (i,k+6,n-2,ires)
         call findatm (' HB2',start,stop,i)
         call newatm (i,k+1,n-4,1.10d0,cai,107.9d0,n-3,110.0d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,k+1,n-5,1.10d0,cai,107.9d0,n-4,110.0d0,-1)
         call findatm (' HG ',start,stop,i)
         call newatm (i,k+3,n-5,1.10d0,n-6,107.0d0,n-4,108.2d0,1)
         call findatm ('HD11',start,stop,i)
         call newatm (i,k+5,n-5,1.10d0,n-6,111.6d0,n-7,180.0d0,0)
         call findatm ('HD12',start,stop,i)
         call newatm (i,k+5,n-6,1.10d0,n-7,111.6d0,n-8,60.0d0,0)
         call findatm ('HD13',start,stop,i)
         call newatm (i,k+5,n-7,1.10d0,n-8,111.6d0,n-9,-60.0d0,0)
         call findatm ('HD21',start,stop,i)
         call newatm (i,k+7,n-7,1.10d0,n-9,111.6d0,n-10,180.0d0,0)
         call findatm ('HD22',start,stop,i)
         call newatm (i,k+7,n-8,1.10d0,n-10,111.6d0,n-11,60.0d0,0)
         call findatm ('HD23',start,stop,i)
         call newatm (i,k+7,n-9,1.10d0,n-11,111.6d0,n-12,-60.0d0,0)
c
c     isoleucine residue  (ILE)
c
      else if (resname .eq. 'ILE') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,k,cai,ires)
         call findatm (' CG1',start,stop,i)
         call oldatm (i,k+2,n-1,ires)
         call findatm (' CG2',start,stop,i)
         call oldatm (i,k+4,n-2,ires)
         call findatm (' CD1',start,stop,i)
         if (i .eq. 0)  call findatm (' CD ',start,stop,i)
         call oldatm (i,k+6,n-2,ires)
         call findatm (' HB ',start,stop,i)
         call newatm (i,k+1,n-4,1.10d0,cai,107.0d0,n-3,108.2d0,-1)
         call findatm ('HG12',start,stop,i)
         call newatm (i,k+3,n-4,1.10d0,n-5,109.5d0,n-2,109.5d0,1)
         call findatm ('HG13',start,stop,i)
         call newatm (i,k+3,n-5,1.10d0,n-6,109.5d0,n-3,109.5d0,-1)
         call findatm ('HG21',start,stop,i)
         call newatm (i,k+5,n-5,1.10d0,n-7,111.6d0,cai,180.0d0,0)
         call findatm ('HG22',start,stop,i)
         call newatm (i,k+5,n-6,1.10d0,n-8,111.6d0,cai,60.0d0,0)
         call findatm ('HG23',start,stop,i)
         call newatm (i,k+5,n-7,1.10d0,n-9,111.6d0,cai,-60.0d0,0)
         call findatm ('HD11',start,stop,i)
         call newatm (i,k+7,n-7,1.10d0,n-9,111.6d0,n-10,180.0d0,0)
         call findatm ('HD12',start,stop,i)
         call newatm (i,k+7,n-8,1.10d0,n-10,111.6d0,n-11,60.0d0,0)
         call findatm ('HD13',start,stop,i)
         call newatm (i,k+7,n-9,1.10d0,n-11,111.6d0,n-12,-60.0d0,0)
c
c     serine residue  (SER)
c
      else if (resname .eq. 'SER') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,k,cai,ires)
         call findatm (' OG ',start,stop,i)
         call oldatm (i,k+2,n-1,ires)
         call findatm (' HB2',start,stop,i)
         call newatm (i,k+1,n-2,1.10d0,cai,109.2d0,n-1,109.5d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,k+1,n-3,1.10d0,cai,109.2d0,n-2,109.5d0,-1)
         call findatm (' HG ',start,stop,i)
         call newatm (i,k+3,n-3,0.94d0,n-4,106.9d0,cai,180.0d0,0)
c
c     threonine residue  (THR)
c
      else if (resname .eq. 'THR') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,k,cai,ires)
         call findatm (' OG1',start,stop,i)
         call oldatm (i,k+2,n-1,ires)
         call findatm (' CG2',start,stop,i)
         call oldatm (i,k+4,n-2,ires)
         call findatm (' HB ',start,stop,i)
         call newatm (i,k+1,n-3,1.10d0,cai,107.0d0,n-2,108.2d0,-1)
         call findatm (' HG1',start,stop,i)
         call newatm (i,k+3,n-3,0.94d0,n-4,106.9d0,cai,180.0d0,0)
         call findatm ('HG21',start,stop,i)
         call newatm (i,k+5,n-3,1.10d0,n-5,111.6d0,cai,180.0d0,0)
         call findatm ('HG22',start,stop,i)
         call newatm (i,k+5,n-4,1.10d0,n-6,111.6d0,cai,60.0d0,0)
         call findatm ('HG23',start,stop,i)
         call newatm (i,k+5,n-5,1.10d0,n-7,111.6d0,cai,-60.0d0,0)
c
c     cysteine residue  (CYS)
c
      else if (resname .eq. 'CYS') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,k,cai,ires)
         call findatm (' SG ',start,stop,i)
         call oldatm (i,k+2,n-1,ires)
         call findatm (' HB2',start,stop,i)
         call newatm (i,k+1,n-2,1.10d0,cai,109.5d0,n-1,107.5d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,k+1,n-3,1.10d0,cai,109.5d0,n-2,107.5d0,-1)
         call findatm (' HG ',start,stop,i)
         call newatm (i,k+3,n-3,1.34d0,n-4,96.0d0,cai,180.0d0,0)
c
c     cystine residue  (CYX)
c
      else if (resname .eq. 'CYX') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,k,cai,ires)
         call findatm (' SG ',start,stop,i)
         si = n
         call oldatm (i,k+2,n-1,ires)
         call findatm (' HB2',start,stop,i)
         call newatm (i,k+1,n-2,1.10d0,cai,109.5d0,n-1,107.5d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,k+1,n-3,1.10d0,cai,109.5d0,n-2,107.5d0,-1)
c
c     deprotonated cysteine residue  (CYD)
c
      else if (resname .eq. 'CYD') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,k,cai,ires)
         call findatm (' SG ',start,stop,i)
         call oldatm (i,k+2,n-1,ires)
         call findatm (' HB2',start,stop,i)
         call newatm (i,k+1,n-2,1.10d0,cai,109.5d0,n-1,107.5d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,k+1,n-3,1.10d0,cai,109.5d0,n-2,107.5d0,-1)
c
c     proline residue  (PRO)
c
      else if (resname .eq. 'PRO') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,k,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,k+2,n-1,ires)
         call findatm (' CD ',start,stop,i)
         if (ires .eq. 1) then
            call oldatm (i,482,n-1,ires)
         else
            call oldatm (i,k+4,n-1,ires)
         end if
         call addbond (n-1,ni)
         call findatm (' HB2',start,stop,i)
         call newatm (i,k+1,n-3,1.10d0,cai,111.2d0,n-2,111.2d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,k+1,n-4,1.10d0,cai,111.2d0,n-3,111.2d0,-1)
         call findatm (' HG2',start,stop,i)
         call newatm (i,k+3,n-4,1.10d0,n-5,111.2d0,n-3,111.2d0,1)
         call findatm (' HG3',start,stop,i)
         call newatm (i,k+3,n-5,1.10d0,n-6,111.2d0,n-4,111.2d0,-1)
         if (ires .eq. 1) then
            call findatm (' HD2',start,stop,i)
            call newatm (i,483,n-5,1.10d0,n-6,111.2d0,ni,111.2d0,1)
            call findatm (' HD3',start,stop,i)
            call newatm (i,483,n-6,1.10d0,n-7,111.2d0,ni,111.2d0,-1)
         else
            call findatm (' HD2',start,stop,i)
            call newatm (i,k+5,n-5,1.10d0,n-6,111.2d0,ni,111.2d0,1)
            call findatm (' HD3',start,stop,i)
            call newatm (i,k+5,n-6,1.10d0,n-7,111.2d0,ni,111.2d0,-1)
         end if
c
c     phenylalanine residue  (PHE)
c
      else if (resname .eq. 'PHE') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,k,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,k+2,n-1,ires)
         call findatm (' CD1',start,stop,i)
         call oldatm (i,k+3,n-1,ires)
         call findatm (' CD2',start,stop,i)
         call oldatm (i,k+3,n-2,ires)
         call findatm (' CE1',start,stop,i)
         call oldatm (i,k+5,n-2,ires)
         call findatm (' CE2',start,stop,i)
         call oldatm (i,k+5,n-2,ires)
         call findatm (' CZ ',start,stop,i)
         call oldatm (i,k+7,n-1,ires)
         call addbond (n-1,n-3)
         call findatm (' HB2',start,stop,i)
         call newatm (i,k+1,n-7,1.10d0,cai,107.9d0,n-6,110.0d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,k+1,n-8,1.10d0,cai,107.9d0,n-7,110.0d0,-1)
         call findatm (' HD1',start,stop,i)
         call newatm (i,k+4,n-7,1.09d0,n-8,120.0d0,n-9,0.0d0,0)
         call findatm (' HD2',start,stop,i)
         call newatm (i,k+4,n-7,1.09d0,n-9,120.0d0,n-10,0.0d0,0)
         call findatm (' HE1',start,stop,i)
         call newatm (i,k+6,n-7,1.09d0,n-9,120.0d0,n-10,180.0d0,0)
         call findatm (' HE2',start,stop,i)
         call newatm (i,k+6,n-7,1.09d0,n-9,120.0d0,n-11,180.0d0,0)
         call findatm (' HZ ',start,stop,i)
         call newatm (i,k+8,n-7,1.09d0,n-8,120.0d0,n-10,180.0d0,0)
c
c     tyrosine residue  (TYR)
c
      else if (resname .eq. 'TYR') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,k,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,k+2,n-1,ires)
         call findatm (' CD1',start,stop,i)
         call oldatm (i,k+3,n-1,ires)
         call findatm (' CD2',start,stop,i)
         call oldatm (i,k+3,n-2,ires)
         call findatm (' CE1',start,stop,i)
         call oldatm (i,k+5,n-2,ires)
         call findatm (' CE2',start,stop,i)
         call oldatm (i,k+5,n-2,ires)
         call findatm (' CZ ',start,stop,i)
         call oldatm (i,k+7,n-1,ires)
         call addbond (n-1,n-3)
         call findatm (' OH ',start,stop,i)
         call oldatm (i,k+8,n-1,ires)
         call findatm (' HB2',start,stop,i)
         call newatm (i,k+1,n-8,1.10d0,cai,107.9d0,n-7,110.0d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,k+1,n-9,1.10d0,cai,107.9d0,n-8,110.0d0,-1)
         call findatm (' HD1',start,stop,i)
         call newatm (i,k+4,n-8,1.09d0,n-9,120.0d0,n-10,0.0d0,0)
         call findatm (' HD2',start,stop,i)
         call newatm (i,k+4,n-8,1.09d0,n-10,120.0d0,n-11,0.0d0,0)
         call findatm (' HE1',start,stop,i)
         call newatm (i,k+6,n-8,1.09d0,n-10,120.0d0,n-11,180.0d0,0)
         call findatm (' HE2',start,stop,i)
         call newatm (i,k+6,n-8,1.09d0,n-10,120.0d0,n-12,180.0d0,0)
         call findatm (' HH ',start,stop,i)
         call newatm (i,k+9,n-7,0.97d0,n-8,108.0d0,n-9,0.0d0,0)
c
c     deprotonated tyrosine residue  (TYD)
c
      else if (resname .eq. 'TYD') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,k,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,k+2,n-1,ires)
         call findatm (' CD1',start,stop,i)
         call oldatm (i,k+3,n-1,ires)
         call findatm (' CD2',start,stop,i)
         call oldatm (i,k+3,n-2,ires)
         call findatm (' CE1',start,stop,i)
         call oldatm (i,k+5,n-2,ires)
         call findatm (' CE2',start,stop,i)
         call oldatm (i,k+5,n-2,ires)
         call findatm (' CZ ',start,stop,i)
         call oldatm (i,k+7,n-1,ires)
         call addbond (n-1,n-3)
         call findatm (' OH ',start,stop,i)
         call oldatm (i,k+8,n-1,ires)
         call findatm (' HB2',start,stop,i)
         call newatm (i,k+1,n-8,1.10d0,cai,107.9d0,n-7,110.0d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,k+1,n-9,1.10d0,cai,107.9d0,n-8,110.0d0,-1)
         call findatm (' HD1',start,stop,i)
         call newatm (i,k+4,n-8,1.09d0,n-9,120.0d0,n-10,0.0d0,0)
         call findatm (' HD2',start,stop,i)
         call newatm (i,k+4,n-8,1.09d0,n-10,120.0d0,n-11,0.0d0,0)
         call findatm (' HE1',start,stop,i)
         call newatm (i,k+6,n-8,1.09d0,n-10,120.0d0,n-11,180.0d0,0)
         call findatm (' HE2',start,stop,i)
         call newatm (i,k+6,n-8,1.09d0,n-10,120.0d0,n-12,180.0d0,0)
c
c     tryptophan residue  (TRP)
c
      else if (resname .eq. 'TRP') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,k,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,k+2,n-1,ires)
         call findatm (' CD1',start,stop,i)
         call oldatm (i,k+3,n-1,ires)
         call findatm (' CD2',start,stop,i)
         call oldatm (i,k+5,n-2,ires)
         call findatm (' NE1',start,stop,i)
         call oldatm (i,k+6,n-2,ires)
         call findatm (' CE2',start,stop,i)
         call oldatm (i,k+8,n-1,ires)
         call addbond (n-1,n-3)
         call findatm (' CE3',start,stop,i)
         call oldatm (i,k+9,n-3,ires)
         call findatm (' CZ2',start,stop,i)
         call oldatm (i,k+11,n-2,ires)
         call findatm (' CZ3',start,stop,i)
         call oldatm (i,k+13,n-2,ires)
         call findatm (' CH2',start,stop,i)
         call oldatm (i,k+15,n-1,ires)
         call addbond (n-1,n-3)
         call findatm (' HB2',start,stop,i)
         call newatm (i,k+1,n-10,1.10d0,cai,107.9d0,n-9,110.0d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,k+1,n-11,1.10d0,cai,107.9d0,n-10,110.0d0,-1)
         call findatm (' HD1',start,stop,i)
         call newatm (i,k+4,n-10,1.09d0,n-11,126.0d0,n-12,0.0d0,0)
         call findatm (' HE1',start,stop,i)
         call newatm (i,k+7,n-9,1.01d0,n-11,126.3d0,n-12,180.0d0,0)
         call findatm (' HE3',start,stop,i)
         call newatm (i,k+10,n-8,1.09d0,n-6,120.0d0,n-5,180.0d0,0)
         call findatm (' HZ2',start,stop,i)
         call newatm (i,k+12,n-8,1.09d0,n-6,120.0d0,n-7,180.0d0,0)
         call findatm (' HZ3',start,stop,i)
         call newatm (i,k+14,n-8,1.09d0,n-7,120.0d0,n-9,180.0d0,0)
         call findatm (' HH2',start,stop,i)
         call newatm (i,k+16,n-8,1.09d0,n-9,120.0d0,n-11,180.0d0,0)
c
c     histidine (HD and HE) residue  (HIS)
c
      else if (resname .eq. 'HIS') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,k,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,k+2,n-1,ires)
         call findatm (' ND1',start,stop,i)
         call oldatm (i,k+3,n-1,ires)
         call findatm (' CD2',start,stop,i)
         call oldatm (i,k+5,n-2,ires)
         call findatm (' CE1',start,stop,i)
         call oldatm (i,k+7,n-2,ires)
         call findatm (' NE2',start,stop,i)
         call oldatm (i,k+9,n-1,ires)
         call addbond (n-1,n-3)
         call findatm (' HB2',start,stop,i)
         call newatm (i,k+1,n-6,1.10d0,cai,107.9d0,n-5,110.0d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,k+1,n-7,1.10d0,cai,107.9d0,n-6,110.0d0,-1)
         call findatm (' HD1',start,stop,i)
         call newatm (i,k+4,n-6,1.02d0,n-4,126.0d0,n-3,180.0d0,0)
         call findatm (' HD2',start,stop,i)
         call newatm (i,k+6,n-6,1.09d0,n-4,126.0d0,n-5,180.0d0,0)
         call findatm (' HE1',start,stop,i)
         call newatm (i,k+8,n-6,1.09d0,n-5,126.0d0,n-7,180.0d0,0)
         call findatm (' HE2',start,stop,i)
         call newatm (i,k+10,n-6,1.02d0,n-7,126.0d0,n-9,180.0d0,0)
c
c     histidine (HD only) residue  (HID)
c
      else if (resname .eq. 'HID') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,k,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,k+2,n-1,ires)
         call findatm (' ND1',start,stop,i)
         call oldatm (i,k+3,n-1,ires)
         call findatm (' CD2',start,stop,i)
         call oldatm (i,k+5,n-2,ires)
         call findatm (' CE1',start,stop,i)
         call oldatm (i,k+7,n-2,ires)
         call findatm (' NE2',start,stop,i)
         call oldatm (i,k+9,n-1,ires)
         call addbond (n-1,n-3)
         call findatm (' HB2',start,stop,i)
         call newatm (i,k+1,n-6,1.10d0,cai,107.9d0,n-5,110.0d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,k+1,n-7,1.10d0,cai,107.9d0,n-6,110.0d0,-1)
         call findatm (' HD1',start,stop,i)
         call newatm (i,k+4,n-6,1.02d0,n-4,126.0d0,n-3,180.0d0,0)
         call findatm (' HD2',start,stop,i)
         call newatm (i,k+6,n-6,1.09d0,n-4,126.0d0,n-5,180.0d0,0)
         call findatm (' HE1',start,stop,i)
         call newatm (i,k+8,n-6,1.09d0,n-5,126.0d0,n-7,180.0d0,0)
c
c     histidine (HE only) residue  (HIE)
c
      else if (resname .eq. 'HIE') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,k,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,k+2,n-1,ires)
         call findatm (' ND1',start,stop,i)
         call oldatm (i,k+3,n-1,ires)
         call findatm (' CD2',start,stop,i)
         call oldatm (i,k+4,n-2,ires)
         call findatm (' CE1',start,stop,i)
         call oldatm (i,k+6,n-2,ires)
         call findatm (' NE2',start,stop,i)
         call oldatm (i,k+8,n-1,ires)
         call addbond (n-1,n-3)
         call findatm (' HB2',start,stop,i)
         call newatm (i,k+1,n-6,1.10d0,cai,107.9d0,n-5,110.0d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,k+1,n-7,1.10d0,cai,107.9d0,n-6,110.0d0,-1)
         call findatm (' HD2',start,stop,i)
         call newatm (i,k+5,n-5,1.09d0,n-3,126.0d0,n-4,180.0d0,0)
         call findatm (' HE1',start,stop,i)
         call newatm (i,k+7,n-5,1.09d0,n-4,126.0d0,n-6,180.0d0,0)
         call findatm (' HE2',start,stop,i)
         call newatm (i,k+9,n-5,1.02d0,n-6,126.0d0,n-8,180.0d0,0)
c
c     aspartic acid residue  (ASP)
c
      else if (resname .eq. 'ASP') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,k,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,k+2,n-1,ires)
         call findatm (' OD1',start,stop,i)
         call oldatm (i,k+3,n-1,ires)
         call findatm (' OD2',start,stop,i)
         call oldatm (i,k+3,n-2,ires)
         call findatm (' HB2',start,stop,i)
         call newatm (i,k+1,n-4,1.10d0,cai,107.9d0,n-3,110.0d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,k+1,n-5,1.10d0,cai,107.9d0,n-4,110.0d0,-1)
c
c     protonated aspartic acid residue  (ASH)
c
      else if (resname .eq. 'ASH') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,k,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,k+2,n-1,ires)
         call findatm (' OD1',start,stop,i)
         call oldatm (i,k+3,n-1,ires)
         call findatm (' OD2',start,stop,i)
         call oldatm (i,k+4,n-2,ires)
         call findatm (' HB2',start,stop,i)
         call newatm (i,k+1,n-4,1.10d0,cai,107.9d0,n-3,110.0d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,k+1,n-5,1.10d0,cai,107.9d0,n-4,110.0d0,-1)
         call findatm (' HD2',start,stop,i)
         call newatm (i,k+5,n-3,0.98d0,n-5,108.7d0,n-4,0.0d0,0)
c
c     asparagine residue  (ASN)
c
      else if (resname .eq. 'ASN') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,k,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,k+2,n-1,ires)
         call findatm (' OD1',start,stop,i)
         call oldatm (i,k+3,n-1,ires)
         call findatm (' ND2',start,stop,i)
         call oldatm (i,k+4,n-2,ires)
         call findatm (' HB2',start,stop,i)
         call newatm (i,k+1,n-4,1.10d0,cai,107.9d0,n-3,110.0d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,k+1,n-5,1.10d0,cai,107.9d0,n-4,110.0d0,-1)
         call findatm ('HD21',start,stop,i)
         call newatm (i,k+5,n-3,1.01d0,n-5,120.9d0,n-6,0.0d0,0)
         call findatm ('HD22',start,stop,i)
         call newatm (i,k+5,n-4,1.01d0,n-6,120.3d0,n-7,180.0d0,0)
c
c     glutamic acid residue  (GLU)
c
      else if (resname .eq. 'GLU') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,k,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,k+2,n-1,ires)
         call findatm (' CD ',start,stop,i)
         call oldatm (i,k+4,n-1,ires)
         call findatm (' OE1',start,stop,i)
         call oldatm (i,k+5,n-1,ires)
         call findatm (' OE2',start,stop,i)
         call oldatm (i,k+5,n-2,ires)
         call findatm (' HB2',start,stop,i)
         call newatm (i,k+1,n-5,1.10d0,cai,107.9d0,n-4,110.0d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,k+1,n-6,1.10d0,cai,107.9d0,n-5,110.0d0,-1)
         call findatm (' HG2',start,stop,i)
         call newatm (i,k+3,n-6,1.10d0,n-7,109.5d0,n-5,109.5d0,1)
         call findatm (' HG3',start,stop,i)
         call newatm (i,k+3,n-7,1.10d0,n-8,109.5d0,n-6,109.5d0,-1)
c
c     protonated glutamic acid residue  (GLH)
c
      else if (resname .eq. 'GLH') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,k,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,k+2,n-1,ires)
         call findatm (' CD ',start,stop,i)
         call oldatm (i,k+4,n-1,ires)
         call findatm (' OE1',start,stop,i)
         call oldatm (i,k+5,n-1,ires)
         call findatm (' OE2',start,stop,i)
         call oldatm (i,k+6,n-2,ires)
         call findatm (' HB2',start,stop,i)
         call newatm (i,k+1,n-5,1.10d0,cai,107.9d0,n-4,110.0d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,k+1,n-6,1.10d0,cai,107.9d0,n-5,110.0d0,-1)
         call findatm (' HG2',start,stop,i)
         call newatm (i,k+3,n-6,1.10d0,n-7,109.5d0,n-5,109.5d0,1)
         call findatm (' HG3',start,stop,i)
         call newatm (i,k+3,n-7,1.10d0,n-8,109.5d0,n-6,109.5d0,-1)
         call findatm (' HE2',start,stop,i)
         call newatm (i,k+7,n-5,0.98d0,n-7,108.7d0,n-6,0.0d0,0)
c
c     glutamine residue  (GLN)
c
      else if (resname .eq. 'GLN') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,k,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,k+2,n-1,ires)
         call findatm (' CD ',start,stop,i)
         call oldatm (i,k+4,n-1,ires)
         call findatm (' OE1',start,stop,i)
         call oldatm (i,k+5,n-1,ires)
         call findatm (' NE2',start,stop,i)
         call oldatm (i,k+6,n-2,ires)
         call findatm (' HB2',start,stop,i)
         call newatm (i,k+1,n-5,1.10d0,cai,107.9d0,n-4,110.0d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,k+1,n-6,1.10d0,cai,107.9d0,n-5,110.0d0,-1)
         call findatm (' HG2',start,stop,i)
         call newatm (i,k+3,n-6,1.10d0,n-7,109.5d0,n-5,109.5d0,1)
         call findatm (' HG3',start,stop,i)
         call newatm (i,k+3,n-7,1.10d0,n-8,109.5d0,n-6,109.5d0,-1)
         call findatm ('HE21',start,stop,i)
         call newatm (i,k+7,n-5,1.01d0,n-7,120.9d0,n-8,0.0d0,0)
         call findatm ('HE22',start,stop,i)
         call newatm (i,k+7,n-6,1.01d0,n-8,120.3d0,n-9,180.0d0,0)
c
c     methionine residue  (MET)
c
      else if (resname .eq. 'MET') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,k,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,k+2,n-1,ires)
         call findatm (' SD ',start,stop,i)
         call oldatm (i,k+4,n-1,ires)
         call findatm (' CE ',start,stop,i)
         call oldatm (i,k+5,n-1,ires)
         call findatm (' HB2',start,stop,i)
         call newatm (i,k+1,n-4,1.10d0,cai,107.9d0,n-3,110.0d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,k+1,n-5,1.10d0,cai,107.9d0,n-4,110.0d0,-1)
         call findatm (' HG2',start,stop,i)
         call newatm (i,k+3,n-5,1.10d0,n-6,109.5d0,n-4,109.5d0,1)
         call findatm (' HG3',start,stop,i)
         call newatm (i,k+3,n-6,1.10d0,n-7,109.5d0,n-5,109.5d0,-1)
         call findatm (' HE1',start,stop,i)
         call newatm (i,k+6,n-5,1.10d0,n-6,110.2d0,n-7,180.0d0,0)
         call findatm (' HE2',start,stop,i)
         call newatm (i,k+6,n-6,1.10d0,n-7,110.2d0,n-8,60.0d0,0)
         call findatm (' HE3',start,stop,i)
         call newatm (i,k+6,n-7,1.10d0,n-8,110.2d0,n-9,-60.0d0,0)
c
c     lysine residue  (LYS)
c
      else if (resname .eq. 'LYS') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,k,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,k+2,n-1,ires)
         call findatm (' CD ',start,stop,i)
         call oldatm (i,k+4,n-1,ires)
         call findatm (' CE ',start,stop,i)
         call oldatm (i,k+6,n-1,ires)
         call findatm (' NZ ',start,stop,i)
         call oldatm (i,k+8,n-1,ires)
         call findatm (' HB2',start,stop,i)
         call newatm (i,k+1,n-5,1.10d0,cai,107.9d0,n-4,110.0d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,k+1,n-6,1.10d0,cai,107.9d0,n-5,110.0d0,-1)
         call findatm (' HG2',start,stop,i)
         call newatm (i,k+3,n-6,1.10d0,n-7,109.5d0,n-5,109.5d0,1)
         call findatm (' HG3',start,stop,i)
         call newatm (i,k+3,n-7,1.10d0,n-8,109.5d0,n-6,109.5d0,-1)
         call findatm (' HD2',start,stop,i)
         call newatm (i,k+5,n-7,1.10d0,n-8,109.5d0,n-6,109.5d0,1)
         call findatm (' HD3',start,stop,i)
         call newatm (i,k+5,n-8,1.10d0,n-9,109.5d0,n-7,109.5d0,-1)
         call findatm (' HE2',start,stop,i)
         call newatm (i,k+7,n-8,1.10d0,n-9,110.9d0,n-7,107.3d0,1)
         call findatm (' HE3',start,stop,i)
         call newatm (i,k+7,n-9,1.10d0,n-10,110.9d0,n-8,107.3d0,-1)
         call findatm (' HZ1',start,stop,i)
         call newatm (i,k+9,n-9,1.04d0,n-10,110.5d0,n-11,180.0d0,0)
         call findatm (' HZ2',start,stop,i)
         call newatm (i,k+9,n-10,1.04d0,n-11,110.5d0,n-12,60.0d0,0)
         call findatm (' HZ3',start,stop,i)
         call newatm (i,k+9,n-11,1.04d0,n-12,110.5d0,n-13,-60.0d0,0)
c
c     deprotonated lysine residue  (LYD)
c
      else if (resname .eq. 'LYD') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,k,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,k+2,n-1,ires)
         call findatm (' CD ',start,stop,i)
         call oldatm (i,k+4,n-1,ires)
         call findatm (' CE ',start,stop,i)
         call oldatm (i,k+6,n-1,ires)
         call findatm (' NZ ',start,stop,i)
         call oldatm (i,k+8,n-1,ires)
         call findatm (' HB2',start,stop,i)
         call newatm (i,k+1,n-5,1.10d0,cai,107.9d0,n-4,110.0d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,k+1,n-6,1.10d0,cai,107.9d0,n-5,110.0d0,-1)
         call findatm (' HG2',start,stop,i)
         call newatm (i,k+3,n-6,1.10d0,n-7,109.5d0,n-5,109.5d0,1)
         call findatm (' HG3',start,stop,i)
         call newatm (i,k+3,n-7,1.10d0,n-8,109.5d0,n-6,109.5d0,-1)
         call findatm (' HD2',start,stop,i)
         call newatm (i,k+5,n-7,1.10d0,n-8,109.5d0,n-6,109.5d0,1)
         call findatm (' HD3',start,stop,i)
         call newatm (i,k+5,n-8,1.10d0,n-9,109.5d0,n-7,109.5d0,-1)
         call findatm (' HE2',start,stop,i)
         call newatm (i,k+7,n-8,1.10d0,n-9,110.9d0,n-7,107.3d0,1)
         call findatm (' HE3',start,stop,i)
         call newatm (i,k+7,n-9,1.10d0,n-10,110.9d0,n-8,107.3d0,-1)
         call findatm (' HZ1',start,stop,i)
         call newatm (i,k+9,n-9,1.04d0,n-10,110.5d0,n-11,180.0d0,0)
         call findatm (' HZ2',start,stop,i)
         call newatm (i,k+9,n-10,1.04d0,n-11,110.5d0,n-12,60.0d0,0)
c
c     arginine residue  (ARG)
c
      else if (resname .eq. 'ARG') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,k,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,k+2,n-1,ires)
         call findatm (' CD ',start,stop,i)
         call oldatm (i,k+4,n-1,ires)
         call findatm (' NE ',start,stop,i)
         call oldatm (i,k+6,n-1,ires)
         call findatm (' CZ ',start,stop,i)
         call oldatm (i,k+8,n-1,ires)
         call findatm (' NH1',start,stop,i)
         call oldatm (i,k+9,n-1,ires)
         call findatm (' NH2',start,stop,i)
         call oldatm (i,k+9,n-2,ires)
         call findatm (' HB2',start,stop,i)
         call newatm (i,k+1,n-7,1.10d0,cai,107.9d0,n-6,110.0d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,k+1,n-8,1.10d0,cai,107.9d0,n-7,110.0d0,-1)
         call findatm (' HG2',start,stop,i)
         call newatm (i,k+3,n-8,1.10d0,n-9,109.5d0,n-7,109.5d0,1)
         call findatm (' HG3',start,stop,i)
         call newatm (i,k+3,n-9,1.10d0,n-10,109.5d0,n-8,109.5d0,-1)
         call findatm (' HD2',start,stop,i)
         call newatm (i,k+5,n-9,1.10d0,n-10,109.5d0,n-8,109.5d0,1)
         call findatm (' HD3',start,stop,i)
         call newatm (i,k+5,n-10,1.10d0,n-11,109.5d0,n-9,109.5d0,-1)
         call findatm (' HE ',start,stop,i)
         call newatm (i,k+7,n-10,1.01d0,n-11,118.5d0,n-9,120.0d0,1)
         call findatm ('HH11',start,stop,i)
         call newatm (i,k+10,n-9,1.01d0,n-10,122.5d0,n-11,0.0d0,0)
         call findatm ('HH12',start,stop,i)
         call newatm (i,k+10,n-10,1.01d0,n-11,118.8d0,n-12,180.0d0,0)
         call findatm ('HH21',start,stop,i)
         call newatm (i,k+10,n-10,1.01d0,n-12,122.5d0,n-13,0.0d0,0)
         call findatm ('HH22',start,stop,i)
         call newatm (i,k+10,n-11,1.01d0,n-13,118.8d0,n-14,180.0d0,0)
c
c     ornithine residue  (ORN)
c
      else if (resname .eq. 'ORN') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,k,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,k+2,n-1,ires)
         call findatm (' CD ',start,stop,i)
         call oldatm (i,k+4,n-1,ires)
         call findatm (' NE ',start,stop,i)
         call oldatm (i,k+6,n-1,ires)
         call findatm (' HB2',start,stop,i)
         call newatm (i,k+1,n-4,1.10d0,cai,107.9d0,n-3,110.0d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,k+1,n-5,1.10d0,cai,107.9d0,n-4,110.0d0,-1)
         call findatm (' HG2',start,stop,i)
         call newatm (i,k+3,n-5,1.10d0,n-7,109.5d0,n-4,109.5d0,1)
         call findatm (' HG3',start,stop,i)
         call newatm (i,k+3,n-6,1.10d0,n-8,109.5d0,n-5,109.5d0,-1)
         call findatm (' HD2',start,stop,i)
         call newatm (i,k+5,n-6,1.10d0,n-8,109.5d0,n-5,109.5d0,1)
         call findatm (' HD3',start,stop,i)
         call newatm (i,k+5,n-7,1.10d0,n-9,109.5d0,n-6,109.5d0,-1)
         call findatm (' HE1',start,stop,i)
         call newatm (i,k+7,n-7,1.04d0,n-8,110.5d0,n-9,180.0d0,0)
         call findatm (' HE2',start,stop,i)
         call newatm (i,k+7,n-8,1.04d0,n-9,110.5d0,n-10,60.0d0,0)
         call findatm (' HE3',start,stop,i)
         call newatm (i,k+7,n-9,1.04d0,n-10,110.5d0,n-11,-60.0d0,0)
c
c     methylalanine residue  (AIB)
c
      else if (resname .eq. 'AIB') then
         call findatm (' CB1',start,stop,i)
         call oldatm (i,k,cai,ires)
         call findatm (' CB2',start,stop,i)
         call oldatm (i,k,cai,ires)
         call findatm ('HB11',start,stop,i)
         call newatm (i,k+1,n-2,1.10d0,cai,110.2d0,ni,180.0d0,0)
         call findatm ('HB12',start,stop,i)
         call newatm (i,k+1,n-3,1.10d0,cai,110.2d0,ni,60.0d0,0)
         call findatm ('HB13',start,stop,i)
         call newatm (i,k+1,n-4,1.10d0,cai,110.2d0,ni,-60.0d0,0)
         call findatm ('HB21',start,stop,i)
         call newatm (i,k+1,n-4,1.10d0,cai,110.2d0,ni,180.0d0,0)
         call findatm ('HB22',start,stop,i)
         call newatm (i,k+1,n-5,1.10d0,cai,110.2d0,ni,60.0d0,0)
         call findatm ('HB23',start,stop,i)
         call newatm (i,k+1,n-6,1.10d0,cai,110.2d0,ni,-60.0d0,0)
c
c     pyroglutamic acid residue  (PCA)
c
      else if (resname .eq. 'PCA') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,k,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,k+2,n-1,ires)
         call findatm (' CD ',start,stop,i)
         call oldatm (i,k+4,n-1,ires)
         call addbond (n-1,ni)
         call findatm (' OE ',start,stop,i)
         call oldatm (i,k+5,n-1,ires)
         call findatm (' HB2',start,stop,i)
         call newatm (i,k+1,n-4,1.10d0,cai,111.2d0,n-3,111.2d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,k+1,n-5,1.10d0,cai,111.2d0,n-4,111.2d0,-1)
         call findatm (' HG2',start,stop,i)
         call newatm (i,k+3,n-5,1.10d0,n-6,111.2d0,n-4,111.2d0,1)
         call findatm (' HG3',start,stop,i)
         call newatm (i,k+3,n-6,1.10d0,n-7,111.2d0,n-5,111.2d0,-1)
c
c     unknown residue  (UNK)
c
      else if (resname .eq. 'UNK') then
         k = hatyp(seqtyp(ires))
         if (ires .eq. 1)  k = hantyp(seqtyp(ires))
         if (ires .eq. nseq)  k = hactyp(seqtyp(ires))
         call newatm (i,k,cai,1.10d0,ni,109.5d0,ci,109.5d0,1)
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ligase  --  coordinates from PDB nucleic acid  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ligase" translates a nucleic acid structure in Protein Data
c     Bank format to a Cartesian coordinate file and sequence file
c
c
      subroutine ligase (ichn)
      use sizes
      use atoms
      use files
      use iounit
      use pdb
      use resdue
      use sequen
      implicit none
      integer i,j,k
      integer ichn,ityp
      integer jres,kres
      integer start,stop
      integer poi,o5i,c5i
      integer c4i,o4i,c1i
      integer c3i,c2i,o3i,o2i
      logical newchain,endchain
      logical, allocatable :: deoxy(:)
      character*3 resname
c
c
c     set the next atom and the residue range of the chain
c
      n = n + 1
      jres = ichain(1,ichn)
      kres = ichain(2,ichn)
c
c     perform dynamic allocation of some local arrays
c
      allocate (deoxy(nres))
c
c     check for deoxyribose and change residue name if necessary
c
      do i = jres, kres
         deoxy(i) = .false.
         start = resatm(1,i)
         stop = resatm(2,i)
         resname = pdbres(start)
         call findatm (' O2''',start,stop,k)
         if (k .eq. 0) then
            deoxy(i) = .true.
            do j = start, stop
               if (resname .eq. '  A')  pdbres(j) = ' DA'
               if (resname .eq. '  G')  pdbres(j) = ' DG'
               if (resname .eq. '  C')  pdbres(j) = ' DC'
               if (resname .eq. '  U')  pdbres(j) = ' DU'
               if (resname .eq. '  T')  pdbres(j) = ' DT'
            end do
         end if
      end do
c
c     locate and assign the atoms that make up each residue
c
      do i = jres, kres
         ityp = seqtyp(i)
         start = resatm(1,i)
         stop = resatm(2,i)
         resname = pdbres(start)
c
c     check that the maximum allowed atoms is not exceeded
c
         if (n+25 .gt. maxatm) then
            write (iout,10)  maxatm
   10       format (/,' LIGASE  --  The Maximum of',i8,' Atoms',
     &                 ' has been Exceeded')
            call fatal
         end if
c
c     test for initial or final residue of a nucleotide chain
c
         newchain = .false.
         endchain = .false.
         do j = 1, nchain
            if (i .eq. ichain(1,j)) then
               newchain = .true.
               poi = 0
               o3i = 0
            end if
            if (i .eq. ichain(2,j))  endchain = .true.
         end do
c
c     build the phosphate atoms of the current residue
c
         if (resname .eq. ' TP') then

         else if (resname .eq. ' DP') then

         else if (resname .eq. ' MP') then

         else if (.not. newchain) then
            call findatm (' P  ',start,stop,k)
            if (k .ne. 0)  poi = n
            j = ptyp(ityp)
            call oldatm (k,j,o3i,i)
            call findatm (' OP1',start,stop,k)
            j = optyp(ityp)
            call oldatm (k,j,n-1,i)
            call findatm (' OP2',start,stop,k)
            j = optyp(ityp)
            call oldatm (k,j,n-2,i)
         end if
c
c     build the ribose sugar atoms of the current residue
c
         call findatm (' O5''',start,stop,k)
         if (k .ne. 0)  o5i = n
         j = o5typ(ityp)
         if (newchain) then
            if (deoxy(i)) then
               j = 1244
            else
               j = 1232
            end if
         end if
         call oldatm (k,j,poi,i)
         call findatm (' C5''',start,stop,k)
         if (k .ne. 0)  c5i = n
         j = c5typ(ityp)
         call oldatm (k,j,n-1,i)
         call findatm (' C4''',start,stop,k)
         if (k .ne. 0)  c4i = n
         j = c4typ(ityp)
         call oldatm (k,j,n-1,i)
         call findatm (' O4''',start,stop,k)
         if (k .ne. 0)  o4i = n
         j = o4typ(ityp)
         call oldatm (k,j,n-1,i)
         call findatm (' C1''',start,stop,k)
         if (k .ne. 0)  c1i = n
         j = c1typ(ityp)
         call oldatm (k,j,n-1,i)
         call findatm (' C3''',start,stop,k)
         if (k .ne. 0)  c3i = n
         j = c3typ(ityp)
         call oldatm (k,j,n-3,i)
         call findatm (' C2''',start,stop,k)
         if (k .ne. 0)  c2i = n
         j = c2typ(ityp)
         call oldatm (k,j,n-1,i)
         call addbond (n-1,n-3)
         call findatm (' O3''',start,stop,k)
         if (k .ne. 0)  o3i = n
         j = o3typ(ityp)
         if (endchain) then
            if (deoxy(i)) then
               j = 1249
            else
               j = 1237
            end if
         end if
         call oldatm (k,j,n-2,i)
         if (.not. deoxy(i)) then
            call findatm (' O2''',start,stop,k)
            if (k .ne. 0)  o2i = n
            j = o2typ(ityp)
            call oldatm (k,j,n-2,i)
         end if
c
c     build the hydrogen atoms of the current residue
c
         if (newchain) then
            call findatm (' H5T',start,stop,k)
            j = h5ttyp(ityp)
            call newatm (k,j,o5i,1.00d0,c5i,109.5d0,c4i,180.0d0,0)
         end if
         call findatm (' H5''',start,stop,k)
         j = h51typ(ityp)
         call newatm (k,j,c5i,1.09d0,o5i,109.5d0,c4i,109.5d0,1)
         call findatm ('H5''''',start,stop,k)
         j = h52typ(ityp)
         call newatm (k,j,c5i,1.09d0,o5i,109.5d0,c4i,109.5d0,-1)
         call findatm (' H4''',start,stop,k)
         j = h4typ(ityp)
         call newatm (k,j,c4i,1.09d0,c5i,109.5d0,c3i,109.5d0,-1)
         call findatm (' H3''',start,stop,k)
         j = h3typ(ityp)
         call newatm (k,j,c3i,1.09d0,c4i,109.5d0,c2i,109.5d0,-1)
         if (deoxy(i)) then
            call findatm (' H2''',start,stop,k)
            j = h21typ(ityp)
            call newatm (k,j,c2i,1.09d0,c3i,109.5d0,c1i,109.5d0,-1)
            call findatm ('H2''''',start,stop,k)
            j = h22typ(ityp)
            call newatm (k,j,c2i,1.09d0,c3i,109.5d0,c1i,109.5d0,1)
         else
            call findatm (' H2''',start,stop,k)
            j = h21typ(ityp)
            call newatm (k,j,c2i,1.09d0,c3i,109.5d0,c1i,109.5d0,-1)
            call findatm ('HO2''',start,stop,k)
            j = h22typ(ityp)
            call newatm (k,j,o2i,1.00d0,c2i,109.5d0,c3i,180.0d0,0)
         end if
         call findatm (' H1''',start,stop,k)
         j = h1typ(ityp)
         call newatm (k,j,c1i,1.09d0,o4i,109.5d0,c2i,109.5d0,-1)
         if (endchain) then
            call findatm (' H3T',start,stop,k)
            j = h3ttyp(ityp)
            call newatm (k,j,o3i,1.00d0,c3i,109.5d0,c4i,180.0d0,0)
         end if
c
c     build the standard base atoms of the current residue
c
         call addbase (resname,i,start,stop,c1i)
      end do
c
c     total number of atoms is one less than the current atom
c
      n = n - 1
c
c     perform deallocation of some local arrays
c
      deallocate (deoxy)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine addbase  --  build a single nucleic acid base  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "addbase" builds the Cartesian coordinates for a single nucleic
c     acid base; coordinates are read from the Protein Data Bank file
c     or found from internal coordinates, then atom types are assigned
c     and connectivity data generated
c
c
      subroutine addbase (resname,ires,start,stop,c1i)
      use sizes
      use atoms
      implicit none
      integer i,ires
      integer start,stop
      integer c1i
      character*3 resname
c
c
c     adenine in adenosine residue  (A)
c
      if (resname .eq. '  A') then
         call findatm (' N9 ',start,stop,i)
         call oldatm (i,1017,c1i,ires)
         call findatm (' C8 ',start,stop,i)
         call oldatm (i,1021,n-1,ires)
         call findatm (' N7 ',start,stop,i)
         call oldatm (i,1020,n-1,ires)
         call findatm (' C5 ',start,stop,i)
         call oldatm (i,1019,n-1,ires)
         call findatm (' C6 ',start,stop,i)
         call oldatm (i,1025,n-1,ires)
         call findatm (' N6 ',start,stop,i)
         call oldatm (i,1027,n-1,ires)
         call findatm (' N1 ',start,stop,i)
         call oldatm (i,1024,n-2,ires)
         call findatm (' C2 ',start,stop,i)
         call oldatm (i,1023,n-1,ires)
         call findatm (' N3 ',start,stop,i)
         call oldatm (i,1022,n-1,ires)
         call findatm (' C4 ',start,stop,i)
         call oldatm (i,1018,n-1,ires)
         call addbond (n-1,n-7)
         call addbond (n-1,n-10)
         call findatm (' H8 ',start,stop,i)
         call newatm (i,1030,n-9,1.08d0,n-8,123.1d0,n-7,180.0d0,0)
         call findatm (' H61',start,stop,i)
         call newatm (i,1028,n-6,1.00d0,n-7,120.0d0,n-8,180.0d0,0)
         call findatm (' H62',start,stop,i)
         call newatm (i,1029,n-7,1.00d0,n-8,120.0d0,n-9,0.0d0,0)
         call findatm (' H2 ',start,stop,i)
         call newatm (i,1026,n-6,1.08d0,n-5,115.4d0,n-4,180.0d0,0)
c
c     guanine in guanosine residue  (G)
c
      else if (resname .eq. '  G') then
         call findatm (' N9 ',start,stop,i)
         call oldatm (i,1047,c1i,ires)
         call findatm (' C8 ',start,stop,i)
         call oldatm (i,1051,n-1,ires)
         call findatm (' N7 ',start,stop,i)
         call oldatm (i,1050,n-1,ires)
         call findatm (' C5 ',start,stop,i)
         call oldatm (i,1049,n-1,ires)
         call findatm (' C6 ',start,stop,i)
         call oldatm (i,1055,n-1,ires)
         call findatm (' O6 ',start,stop,i)
         call oldatm (i,1060,n-1,ires)
         call findatm (' N1 ',start,stop,i)
         call oldatm (i,1054,n-2,ires)
         call findatm (' C2 ',start,stop,i)
         call oldatm (i,1053,n-1,ires)
         call findatm (' N2 ',start,stop,i)
         call oldatm (i,1057,n-1,ires)
         call findatm (' N3 ',start,stop,i)
         call oldatm (i,1052,n-2,ires)
         call findatm (' C4 ',start,stop,i)
         call oldatm (i,1048,n-1,ires)
         call addbond (n-1,n-8)
         call addbond (n-1,n-11)
         call findatm (' H8 ',start,stop,i)
         call newatm (i,1061,n-10,1.08d0,n-9,123.0d0,n-8,180.0d0,0)
         call findatm (' H1 ',start,stop,i)
         call newatm (i,1056,n-6,1.00d0,n-8,117.4d0,n-9,180.0d0,0)
         call findatm (' H21',start,stop,i)
         call newatm (i,1058,n-5,1.00d0,n-6,120.0d0,n-7,0.0d0,0)
         call findatm (' H22',start,stop,i)
         call newatm (i,1059,n-6,1.00d0,n-7,120.0d0,n-8,180.0d0,0)
c
c     cytosine in cytidine residue  (C)
c
      else if (resname .eq. '  C') then
         call findatm (' N1 ',start,stop,i)
         call oldatm (i,1078,c1i,ires)
         call findatm (' C2 ',start,stop,i)
         call oldatm (i,1079,n-1,ires)
         call findatm (' O2 ',start,stop,i)
         call oldatm (i,1084,n-1,ires)
         call findatm (' N3 ',start,stop,i)
         call oldatm (i,1080,n-2,ires)
         call findatm (' C4 ',start,stop,i)
         call oldatm (i,1081,n-1,ires)
         call findatm (' N4 ',start,stop,i)
         call oldatm (i,1085,n-1,ires)
         call findatm (' C5 ',start,stop,i)
         call oldatm (i,1082,n-2,ires)
         call findatm (' C6 ',start,stop,i)
         call oldatm (i,1083,n-1,ires)
         call addbond (n-1,n-8)
         call findatm (' H41',start,stop,i)
         call newatm (i,1086,n-3,1.00d0,n-4,120.0d0,n-5,0.0d0,0)
         call findatm (' H42',start,stop,i)
         call newatm (i,1087,n-4,1.00d0,n-5,120.0d0,n-6,180.0d0,0)
         call findatm (' H5 ',start,stop,i)
         call newatm (i,1088,n-4,1.08d0,n-6,121.6d0,n-7,180.0d0,0)
         call findatm (' H6 ',start,stop,i)
         call newatm (i,1089,n-4,1.08d0,n-5,119.4d0,n-7,180.0d0,0)
c
c     uracil in uridine residue  (U)
c
      else if (resname .eq. '  U') then
         call findatm (' N1 ',start,stop,i)
         call oldatm (i,1106,c1i,ires)
         call findatm (' C2 ',start,stop,i)
         call oldatm (i,1107,n-1,ires)
         call findatm (' O2 ',start,stop,i)
         call oldatm (i,1112,n-1,ires)
         call findatm (' N3 ',start,stop,i)
         call oldatm (i,1108,n-2,ires)
         call findatm (' C4 ',start,stop,i)
         call oldatm (i,1109,n-1,ires)
         call findatm (' O4 ',start,stop,i)
         call oldatm (i,1114,n-1,ires)
         call findatm (' C5 ',start,stop,i)
         call oldatm (i,1110,n-2,ires)
         call findatm (' C6 ',start,stop,i)
         call oldatm (i,1111,n-1,ires)
         call addbond (n-1,n-8)
         call findatm (' H3 ',start,stop,i)
         call newatm (i,1113,n-5,1.00d0,n-7,116.5d0,n-8,180.0d0,0)
         call findatm (' H5 ',start,stop,i)
         call newatm (i,1115,n-3,1.08d0,n-5,120.4d0,n-6,180.0d0,0)
         call findatm (' H6 ',start,stop,i)
         call newatm (i,1116,n-3,1.08d0,n-4,118.6d0,n-6,180.0d0,0)
c
c     adenine in deoxyadenosine residue  (DA)
c
      else if (resname .eq. ' DA') then
         call findatm (' N9 ',start,stop,i)
         call oldatm (i,1132,c1i,ires)
         call findatm (' C8 ',start,stop,i)
         call oldatm (i,1136,n-1,ires)
         call findatm (' N7 ',start,stop,i)
         call oldatm (i,1135,n-1,ires)
         call findatm (' C5 ',start,stop,i)
         call oldatm (i,1134,n-1,ires)
         call findatm (' C6 ',start,stop,i)
         call oldatm (i,1140,n-1,ires)
         call findatm (' N6 ',start,stop,i)
         call oldatm (i,1142,n-1,ires)
         call findatm (' N1 ',start,stop,i)
         call oldatm (i,1139,n-2,ires)
         call findatm (' C2 ',start,stop,i)
         call oldatm (i,1138,n-1,ires)
         call findatm (' N3 ',start,stop,i)
         call oldatm (i,1137,n-1,ires)
         call findatm (' C4 ',start,stop,i)
         call oldatm (i,1133,n-1,ires)
         call addbond (n-1,n-7)
         call addbond (n-1,n-10)
         call findatm (' H8 ',start,stop,i)
         call newatm (i,1145,n-9,1.08d0,n-8,123.1d0,n-7,180.0d0,0)
         call findatm (' H61',start,stop,i)
         call newatm (i,1143,n-6,1.00d0,n-7,120.0d0,n-8,180.0d0,0)
         call findatm (' H62',start,stop,i)
         call newatm (i,1144,n-7,1.00d0,n-8,120.0d0,n-9,0.0d0,0)
         call findatm (' H2 ',start,stop,i)
         call newatm (i,1141,n-6,1.08d0,n-5,115.4d0,n-4,180.0d0,0)
c
c     guanine in deoxyguanosine residue  (DG)
c
      else if (resname .eq. ' DG') then
         call findatm (' N9 ',start,stop,i)
         call oldatm (i,1161,c1i,ires)
         call findatm (' C8 ',start,stop,i)
         call oldatm (i,1165,n-1,ires)
         call findatm (' N7 ',start,stop,i)
         call oldatm (i,1164,n-1,ires)
         call findatm (' C5 ',start,stop,i)
         call oldatm (i,1163,n-1,ires)
         call findatm (' C6 ',start,stop,i)
         call oldatm (i,1169,n-1,ires)
         call findatm (' O6 ',start,stop,i)
         call oldatm (i,1174,n-1,ires)
         call findatm (' N1 ',start,stop,i)
         call oldatm (i,1168,n-2,ires)
         call findatm (' C2 ',start,stop,i)
         call oldatm (i,1167,n-1,ires)
         call findatm (' N2 ',start,stop,i)
         call oldatm (i,1171,n-1,ires)
         call findatm (' N3 ',start,stop,i)
         call oldatm (i,1166,n-2,ires)
         call findatm (' C4 ',start,stop,i)
         call oldatm (i,1162,n-1,ires)
         call addbond (n-1,n-8)
         call addbond (n-1,n-11)
         call findatm (' H8 ',start,stop,i)
         call newatm (i,1175,n-10,1.08d0,n-9,123.0d0,n-8,180.0d0,0)
         call findatm (' H1 ',start,stop,i)
         call newatm (i,1170,n-6,1.00d0,n-8,117.4d0,n-9,180.0d0,0)
         call findatm (' H21',start,stop,i)
         call newatm (i,1172,n-5,1.00d0,n-6,120.0d0,n-7,0.0d0,0)
         call findatm (' H22',start,stop,i)
         call newatm (i,1173,n-6,1.00d0,n-7,120.0d0,n-8,180.0d0,0)
c
c     cytosine in deoxycytidine residue  (DC)
c
      else if (resname .eq. ' DC') then
         call findatm (' N1 ',start,stop,i)
         call oldatm (i,1191,c1i,ires)
         call findatm (' C2 ',start,stop,i)
         call oldatm (i,1192,n-1,ires)
         call findatm (' O2 ',start,stop,i)
         call oldatm (i,1197,n-1,ires)
         call findatm (' N3 ',start,stop,i)
         call oldatm (i,1193,n-2,ires)
         call findatm (' C4 ',start,stop,i)
         call oldatm (i,1194,n-1,ires)
         call findatm (' N4 ',start,stop,i)
         call oldatm (i,1198,n-1,ires)
         call findatm (' C5 ',start,stop,i)
         call oldatm (i,1195,n-2,ires)
         call findatm (' C6 ',start,stop,i)
         call oldatm (i,1196,n-1,ires)
         call addbond (n-1,n-8)
         call findatm (' H41',start,stop,i)
         call newatm (i,1199,n-3,1.00d0,n-4,120.0d0,n-5,0.0d0,0)
         call findatm (' H42',start,stop,i)
         call newatm (i,1200,n-4,1.00d0,n-5,120.0d0,n-6,180.0d0,0)
         call findatm (' H5 ',start,stop,i)
         call newatm (i,1201,n-4,1.08d0,n-6,121.6d0,n-7,180.0d0,0)
         call findatm (' H6 ',start,stop,i)
         call newatm (i,1202,n-4,1.08d0,n-5,119.4d0,n-7,180.0d0,0)
c
c     thymine in deoxythymidine residue  (DT)
c
      else if (resname .eq. ' DT') then
         call findatm (' N1 ',start,stop,i)
         call oldatm (i,1218,c1i,ires)
         call findatm (' C2 ',start,stop,i)
         call oldatm (i,1219,n-1,ires)
         call findatm (' O2 ',start,stop,i)
         call oldatm (i,1224,n-1,ires)
         call findatm (' N3 ',start,stop,i)
         call oldatm (i,1220,n-2,ires)
         call findatm (' C4 ',start,stop,i)
         call oldatm (i,1221,n-1,ires)
         call findatm (' O4 ',start,stop,i)
         call oldatm (i,1226,n-1,ires)
         call findatm (' C5 ',start,stop,i)
         call oldatm (i,1222,n-2,ires)
         call findatm (' C7 ',start,stop,i)
         call oldatm (i,1227,n-1,ires)
         call findatm (' C6 ',start,stop,i)
         call oldatm (i,1223,n-2,ires)
         call addbond (n-1,n-9)
         call findatm (' H3 ',start,stop,i)
         call newatm (i,1225,n-6,1.00d0,n-8,116.8d0,n-9,180.0d0,0)
         call findatm (' H71',start,stop,i)
         call newatm (i,1228,n-3,1.09d0,n-4,109.5d0,n-6,0.0d0,0)
         call findatm (' H72',start,stop,i)
         call newatm (i,1228,n-4,1.09d0,n-5,109.5d0,n-1,109.5d0,1)
         call findatm (' H73',start,stop,i)
         call newatm (i,1228,n-5,1.09d0,n-6,109.5d0,n-2,109.5d0,-1)
         call findatm (' H6 ',start,stop,i)
         call newatm (i,1229,n-5,1.08d0,n-7,119.4d0,n-9,180.0d0,0)
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine hetatom  --  coordinates of PDB water and ions  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "hetatom" translates water molecules and ions in Protein Data
c     Bank format to a Cartesian coordinate file and sequence file
c
c
      subroutine hetatom
      use sizes
      use atoms
      use pdb
      implicit none
      integer i
c
c
c     find water molecules and ions in PDB HETATM records
c
      n = n + 1
      i = 0
      do while (i .lt. npdb)
         i = i + 1
         if (pdbtyp(i) .eq. 'HETATM') then
            if (pdbres(i) .eq. 'HOH') then
               if (pdbatm(i) .eq. ' O  ') then
                  call oldatm (i,2001,0,0)
                  if (pdbatm(i+1).eq.' H  ' .and.
     &                pdbatm(i+2).eq.' H  ') then
                     call oldatm (i+1,2002,n-1,0)
                     call oldatm (i+2,2002,n-2,0)
                     i = i + 2
                  else
                     call newatm (0,2002,n-1,0.96d0,n-2,109.5d0,
     &                               n-3,120.0d0,0)
                     call newatm (0,2002,n-2,0.96d0,n-1,109.5d0,
     &                               n-3,120.0d0,0)
                  end if
               end if
            else if (pdbres(i) .eq. 'NA ') then
               call oldatm (i,2003,0,0)
            else if (pdbres(i) .eq. 'K  ') then
               call oldatm (i,2004,0,0)
            else if (pdbres(i) .eq. 'MG ') then
               call oldatm (i,2005,0,0)
            else if (pdbres(i) .eq. 'CA ') then
               call oldatm (i,2006,0,0)
            else if (pdbres(i) .eq. 'CL ') then
               call oldatm (i,2007,0,0)
            end if
         end if
      end do
      n = n - 1
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine oldatm  --  transfer coordinates from PDB  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "oldatm" get the Cartesian coordinates for an atom from
c     the Protein Data Bank file, then assigns the atom type
c     and atomic connectivities
c
c
      subroutine oldatm (i,bionum,i1,ires)
      use sizes
      use atomid
      use atoms
      use fields
      use iounit
      use katoms
      use sequen
      use pdb
      implicit none
      integer i,bionum
      integer i1,ires
c
c
c     get coordinates, assign atom type, and update connectivities
c
      if (bionum .ne. 0) then
         if (i .ne. 0) then
            type(n) = biotyp(bionum)
            if (type(n) .gt. 0) then
               name(n) = symbol(type(n))
            else
               type(n) = 0
               name(n) = '   '
            end if
            x(n) = xpdb(i)
            y(n) = ypdb(i)
            z(n) = zpdb(i)
            if (i1 .ne. 0)  call addbond (n,i1)
            n = n + 1
         else
            write (iout,10)  bionum,ires,seq(ires)
   10       format (/,' OLDATM  --  A PDB Atom of Biotype',i5,
     &                 ' is Missing in Residue',i5,'-',a3)
            call fatal
         end if
      end if
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine newatm  --  create and define a new atom  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "newatm" creates and defines an atom needed for the
c     Cartesian coordinates file, but which may not present
c     in the original Protein Data Bank file
c
c
      subroutine newatm (i,bionum,ia,bond,ib,angle1,ic,angle2,chiral)
      use sizes
      use atomid
      use atoms
      use fields
      use katoms
      use pdb
      implicit none
      integer i,bionum
      integer ia,ib,ic
      integer chiral
      real*8 bond
      real*8 angle1
      real*8 angle2
c
c
c     set the atom type, compute coordinates, assign
c     connectivities and increment the atom counter
c
      if (bionum .ne. 0) then
         type(n) = biotyp(bionum)
         if (type(n) .gt. 0) then
            name(n) = symbol(type(n))
         else
            type(n) = 0
            name(n) = '   '
         end if
         if (i .eq. 0) then
            call xyzatm (n,ia,bond,ib,angle1,ic,angle2,chiral)
         else
            x(n) = xpdb(i)
            y(n) = ypdb(i)
            z(n) = zpdb(i)
         end if
         call addbond (n,ia)
         n = n + 1
      end if
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine addbond  --  add a bond between two atoms  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "addbond" adds entries to the attached atoms list in
c     order to generate a direct connection between two atoms
c
c
      subroutine addbond (i,j)
      use sizes
      use couple
      implicit none
      integer i,j
c
c
c     add connectivity between the two atoms
c
      if (i.ne.0 .and. j.ne.0) then
         n12(i) = n12(i) + 1
         i12(n12(i),i) = j
         n12(j) = n12(j) + 1
         i12(n12(j),j) = i
      end if
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine findatm  --  locate PDB atom in a residue  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "findatm" locates a specific PDB atom name type within a
c     range of atoms from the PDB file, returns zero if the name
c     type was not found
c
c
      subroutine findatm (name,start,stop,ipdb)
      use sizes
      use pdb
      implicit none
      integer i,ipdb
      integer start,stop
      character*4 name
c
c
c     search for the specified atom within the residue
c
      ipdb = 0
      do i = start, stop
         if (pdbatm(i) .eq. name) then
            ipdb = i
            goto 10
         end if
      end do
   10 continue
      return
      end
