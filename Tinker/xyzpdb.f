c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  program xyzpdb  --  Cartesian to Protein Data Bank file  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "xyzpdb" takes as input a Cartesian coordinates file,
c     then converts to and writes out a Protein Data Bank file
c
c
      program xyzpdb
      use files
      use inform
      implicit none
      integer ipdb,ixyz
      integer freeunit
      character*240 pdbfile
      character*240 xyzfile
c
c
c     get the Cartesian coordinates file for the system
c
      call initial
      call getxyz
c
c     get atomic number of each atom and count the molecules
c
      call attach
      call field
      call katom
      call molecule
c
c     open the Protein Data Bank file to be used for output
c
      ipdb = freeunit ()
      pdbfile = filename(1:leng)//'.pdb'
      call version (pdbfile,'new')
      open (unit=ipdb,file=pdbfile,status='new')
c
c     reopen the coordinates file and read the first structure
c
      ixyz = freeunit ()
      xyzfile = filename
      call suffix (xyzfile,'xyz','old')
      open (unit=ixyz,file=xyzfile,status ='old')
      rewind (unit=ixyz)
      call readxyz (ixyz)
c
c     add each successive coordinate frame to the PDB file
c
      do while (.not. abort)
         call makepdb
         call prtpdb (ipdb)
         call readxyz (ixyz)
      end do
c
c     perform any final tasks before program exit
c
      close (unit=ixyz)
      close (unit=ipdb)
      call final
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine makepdb  --  build PDB from Cartesian coords  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "makepdb" cconstructs a Protein Data Bank file from a set
c     of Cartesian coordinates with special handling for systems
c     consisting of biopolymer chains, ligands and water molecules
c
c
      subroutine makepdb
      use sizes
      use atomid
      use atoms
      use couple
      use files
      use molcul
      use pdb
      use resdue
      use sequen
      implicit none
      integer i,j,k,m,ii
      integer kp,ka,kn
      integer iseq,freeunit
      integer start,stop
      integer pdbnum,atmnum
      integer justify,cbi
      integer noxy,nhydro
      integer, allocatable :: ni(:)
      integer, allocatable :: cai(:)
      integer, allocatable :: ci(:)
      integer, allocatable :: oi(:)
      integer, allocatable :: poi(:)
      integer, allocatable :: op1(:)
      integer, allocatable :: op2(:)
      integer, allocatable :: op3(:)
      integer, allocatable :: c5i(:)
      integer, allocatable :: o5i(:)
      integer, allocatable :: c4i(:)
      integer, allocatable :: o4i(:)
      integer, allocatable :: c3i(:)
      integer, allocatable :: o3i(:)
      integer, allocatable :: c2i(:)
      integer, allocatable :: o2i(:)
      integer, allocatable :: c1i(:)
      logical exist,generic
      logical cbone,nbone,obone
      logical first
      logical, allocatable :: water(:)
      logical, allocatable :: hetmol(:)
      character*3 resname
      character*4 atmname
      character*7, allocatable :: restyp(:)
      character*240 seqfile
      save first
      data first  / .true. /
c
c
c     perform dynamic allocation of some global arrays
c
      if (first) then
         first = .false.
         if (.not. allocated(resnum))  allocate (resnum(maxatm))
         if (.not. allocated(resatm))  allocate (resatm(2,maxatm))
         if (.not. allocated(npdb12))  allocate (npdb12(maxatm))
         if (.not. allocated(ipdb12))  allocate (ipdb12(maxval,maxatm))
         if (.not. allocated(pdblist))  allocate (pdblist(maxatm))
         if (.not. allocated(xpdb))  allocate (xpdb(maxatm))
         if (.not. allocated(ypdb))  allocate (ypdb(maxatm))
         if (.not. allocated(zpdb))  allocate (zpdb(maxatm))
         if (.not. allocated(pdbres))  allocate (pdbres(maxatm))
         if (.not. allocated(pdbatm))  allocate (pdbatm(maxatm))
         if (.not. allocated(pdbtyp))  allocate (pdbtyp(maxatm))
      end if
c
c     initialize number of PDB atoms and atom mapping
c
      npdb = 0
      do i = 1, n
         pdblist(i) = 0
      end do
c
c     read the biopolymer sequence file if one exists
c
      iseq = freeunit ()
      seqfile = filename(1:leng)//'.seq'
      call version (seqfile,'old')
      inquire (file=seqfile,exist=exist)
      if (exist) then
         open (unit=iseq,file=seqfile,status='old')
         rewind (unit=iseq)
         call readseq (iseq)
         close (unit=iseq)
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (ni(maxres))
      allocate (cai(maxres))
      allocate (ci(maxres))
      allocate (oi(maxres))
      allocate (poi(maxres))
      allocate (op1(maxres))
      allocate (op2(maxres))
      allocate (op3(maxres))
      allocate (c5i(maxres))
      allocate (o5i(maxres))
      allocate (c4i(maxres))
      allocate (o4i(maxres))
      allocate (c3i(maxres))
      allocate (o3i(maxres))
      allocate (c2i(maxres))
      allocate (o2i(maxres))
      allocate (c1i(maxres))
      allocate (restyp(maxres))
c
c     zero out the backbone atoms for biopolymer sequences
c
      do i = 1, nseq
         ni(i) = 0
         cai(i) = 0
         ci(i) = 0
         oi(i) = 0
         poi(i) = 0
         op1(i) = 0
         op2(i) = 0
         op3(i) = 0
         c5i(i) = 0
         o5i(i) = 0
         c4i(i) = 0
         o4i(i) = 0
         c3i(i) = 0
         o3i(i) = 0
         c2i(i) = 0
         o2i(i) = 0
         c1i(i) = 0
      end do
c
c     set the molecule type for each residue via chain type
c
      generic = .true.
      do i = 1, nchain
         do j = ichain(1,i), ichain(2,i)
            restyp(j) = 'GENERIC'
            if (chntyp(i) .eq. 'PEPTIDE')  restyp(j) = 'PEPTIDE'
            if (chntyp(i) .eq. 'NUCLEIC')  restyp(j) = 'NUCLEIC'
         end do
         if (restyp(j)  .ne. 'GENERIC')  generic = .false.
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (water(nmol))
c
c     check each molecule to see if it is a water molecule
c
      do i = 1, nmol
         water(i) = .false.
         if (imol(2,i)-imol(1,i) .eq. 2) then
            noxy = 0
            nhydro = 0
            do j = imol(1,i), imol(2,i)
               k = kmol(j)
               if (atomic(k) .eq. 8)  noxy = noxy + 1
               if (atomic(k) .eq. 1)  nhydro = nhydro + 1
            end do
            if (noxy.eq.1 .and. nhydro.eq.2)  water(i) = .true.
         end if
      end do
c
c     for general structures, transfer each atom to PDB format
c
      if (generic) then
         do i = 1, nmol
            do j = imol(1,i), imol(2,i)
               k = kmol(j)
               atmname = ' '//name(k)
               if (water(i)) then
                  resname = 'HOH'
               else
                  justify = 0
                  call numeral (type(k),resname,justify)
               end if
               pdbnum = i
               call pdbatom (atmname,resname,pdbnum,k)
               pdbtyp(npdb) = 'HETATM'
            end do
         end do
         do i = 1, nmol
            do j = imol(1,i), imol(2,i)
               k = kmol(j)
               kp = pdblist(k)
               npdb12(kp) = n12(k)
               do m = 1, n12(k)
                  ipdb12(m,kp) = pdblist(i12(m,k))
               end do
            end do
         end do
      end if
c
c     find the amide nitrogens and other peptide backbone atoms
c
      m = 1
      do i = 1, n
         if (restyp(m) .eq. 'PEPTIDE') then
            resname = amino(seqtyp(m))
            if (resname .eq. 'H2N') then
               m = m + 1
            else if (resname .eq. 'FOR') then
               if (atomic(i) .eq. 6) then
                  nbone = .false.
                  obone = .false.
                  do j = 1, n12(i)
                     k = i12(j,i)
                     if (atomic(k) .eq. 7) then
                        nbone = .true.
                     else if (atomic(k) .eq. 8) then
                        obone = .true.
                     end if
                  end do
                  if (nbone .and. obone) then
                     cai(m) = i
                     ci(m) = i
                     oi(m) = i + 1
                     m = m + 1
                  end if
               end if
            else if (resname .eq. 'ACE') then
               if (atomic(i) .eq. 6) then
                  nbone = .false.
                  obone = .false.
                  do j = 1, n13(i)
                     k = i13(j,i)
                     if (atomic(k) .eq. 7) then
                        nbone = .true.
                     else if (atomic(k) .eq. 8) then
                        obone = .true.
                     end if
                  end do
                  if (nbone .and. obone) then
                     cai(m) = i
                     ci(m) = i + 1
                     oi(m) = i + 2
                     m = m + 1
                  end if
               end if
            else if (resname .eq. 'COH') then
               if (n12(i) .gt. 1) then
                  if (atomic(i) .eq. 8) then
                     nbone = .false.
                     obone = .false.
                     do j = 1, n13(i)
                        k = i13(j,i)
                        if (atomic(k) .eq. 8) then
                           obone = .true.
                        end if
                     end do
                     do j = 1, n14(i)
                        k = i14(j,i)
                        if (atomic(k) .eq. 7) then
                           nbone = .true.
                        end if
                     end do
                     if (nbone .and. obone) then
                        ni(m) = i
                        m = m + 1
                     end if
                  end if
               end if
            else if (resname .eq. 'NH2') then
               if (atomic(i) .eq. 7) then
                  nbone = .false.
                  obone = .false.
                  do j = 1, n13(i)
                     k = i13(j,i)
                     if (atomic(k) .eq. 8) then
                        obone = .true.
                     end if
                  end do
                  do j = 1, n14(i)
                     k = i14(j,i)
                     if (atomic(k) .eq. 7) then
                        nbone = .true.
                     end if
                  end do
                  if (nbone .and. obone) then
                     ni(m) = i
                     m = m + 1
                  end if
               end if
            else if (resname .eq. 'NME') then
               if (atomic(i) .eq. 7) then
                  nbone = .false.
                  obone = .false.
                  do j = 1, n13(i)
                     k = i13(j,i)
                     if (atomic(k) .eq. 8) then
                        obone = .true.
                     end if
                  end do
                  do j = 1, n14(i)
                     k = i14(j,i)
                     if (atomic(k) .eq. 7) then
                        nbone = .true.
                     end if
                  end do
                  if (nbone .and. obone) then
                     ni(m) = i
                     cai(m) = i + 1
                     m = m + 1
                  end if
               end if
            else
               if (atomic(i) .eq. 7) then
                  obone = .false.
                  do j = 1, n14(i)
                     k = i14(j,i)
                     if (atomic(k) .eq. 8) then
                        obone = .true.
                     end if
                  end do
                  if (obone) then
                     ni(m) = i
                     cai(m) = i + 1
                     ci(m) = i + 2
                     oi(m) = i + 3
                     m = m + 1
                  end if
               end if
            end if
c
c     find the phosphates and sugar C1 nucleotide backbone atoms
c
         else if (restyp(m) .eq. 'NUCLEIC') then
            resname = nuclz(seqtyp(m))
            if (resname .eq. ' MP') then
               if (atomic(i) .eq. 15) then
                  poi(m) = i
                  m = m + 1
               end if
            end if
            if (atomic(i).eq.6 .and. n12(i).eq.4) then
               cbone = .false.
               nbone = .false.
               obone = .false.
               do j = 1, n12(i)
                  k = i12(j,i)
                  ka = atomic(k)
                  kn = n12(k)
                  if (ka .eq. 6)  cbone = .true.
                  if (ka.eq.7 .and. kn.eq.3)  nbone = .true.
                  if (ka.eq.8 .and. kn.eq.2)  obone = .true.
               end do
               if (cbone .and. nbone .and. obone) then
                  c1i(m) = i
                  m = m + 1
               end if
            end if
         end if
         if (m .gt. nseq)  goto 10
      end do
   10 continue
c
c     find the remainder of the nucleotide backbone atoms
c
      do ii = 1, nchain
         if (chntyp(ii) .eq. 'NUCLEIC') then
            start = ichain(1,ii)
            stop = ichain(2,ii)
            do i = start, stop
               m = c1i(i)
               if (m .ne. 0) then
                  do j = 1, n12(m)
                     k = i12(j,m)
                     ka = atomic(k)
                     if (ka .eq. 6)  c2i(i) = k
                     if (ka .eq. 7)  ni(i) = k
                     if (ka .eq. 8)  o4i(i) = k
                  end do
                  m = o4i(i)
                  do j = 1, n12(m)
                     k = i12(j,m)
                     ka = atomic(k)
                     if (ka.eq.6 .and. k.ne.c1i(i))  c4i(i) = k
                  end do
                  m = c2i(i)
                  do j = 1, n12(m)
                     k = i12(j,m)
                     ka = atomic(k)
                     if (ka .eq. 8)  o2i(i) = k
                     if (ka.eq.6 .and. k.ne.c1i(i))  c3i(i) = k
                  end do
                  m = c3i(i)
                  do j = 1, n12(m)
                     k = i12(j,m)
                     ka = atomic(k)
                     if (ka .eq. 8)  o3i(i) = k
                  end do
                  m = c4i(i)
                  do j = 1, n12(m)
                     k = i12(j,m)
                     ka = atomic(k)
                     if (ka.eq.6 .and. k.ne.c3i(i))  c5i(i) = k
                  end do
                  m = c5i(i)
                  do j = 1, n12(m)
                     k = i12(j,m)
                     ka = atomic(k)
                     if (ka .eq. 8)  o5i(i) = k
                  end do
                  m = o5i(i)
                  do j = 1, n12(m)
                     k = i12(j,m)
                     ka = atomic(k)
                     if (ka .eq. 15)  poi(i) = k
                  end do
               end if
               if (i .gt. 1) then
                  resname = nuclz(seqtyp(i-1))
                  if (resname .eq. ' MP')  poi(i) = 0
                  if (resname .eq. ' DP')  poi(i) = 0
                  if (resname .eq. ' TP')  poi(i) = 0
               end if
               m = poi(i)
               if (m .ne. 0) then
                  do j = 1, n12(m)
                     k = i12(j,m)
                     ka = atomic(k)
                     if (ka.eq.8 .and. n12(k).eq.1) then
                        if (op1(i) .eq. 0) then
                           op1(i) = k
                        else if (op2(i) .eq. 0) then
                           op2(i) = k
                        else
                           op3(i) = k
                        end if
                     end if
                  end do
               end if
            end do
         end if
      end do
c
c     copy the atoms of each biopolymer residue into PDB format
c
      do m = 1, nchain
         start = ichain(1,m)
         stop = ichain(2,m)
         if (chntyp(m) .eq. 'PEPTIDE') then
            do i = start, stop
               resname = amino(seqtyp(i))
               if (resname .eq. 'H2N') then
                  continue
               else if (resname .eq. 'FOR') then
                  call pdbatom (' C  ',resname,i,ci(i))
                  call pdbatom (' O  ',resname,i,oi(i))
               else if (resname .eq. 'ACE') then
                  call pdbatom (' CH3',resname,i,cai(i))
                  call pdbatom (' C  ',resname,i,ci(i))
                  call pdbatom (' O  ',resname,i,oi(i))
               else if (resname .eq. 'COH') then
                  call pdbatom (' OH ',resname,i,ni(i))
               else if (resname .eq. 'NH2') then
                  call pdbatom (' N  ',resname,i,ni(i))
               else if (resname .eq. 'NME') then
                  call pdbatom (' N  ',resname,i,ni(i))
                  call pdbatom (' CH3',resname,i,cai(i))
               else
                  call pdbatom (' N  ',resname,i,ni(i))
                  call pdbatom (' CA ',resname,i,cai(i))
                  call pdbatom (' C  ',resname,i,ci(i))
                  call pdbatom (' O  ',resname,i,oi(i))
               end if
               call getside (resname,i,ci(i),cai(i),cbi)
               if ((resname.eq.'CYS'.or.resname.eq.'CYX')
     &               .and. cbi.ne.0) then
                  resname = 'CYS'
                  do j = 1, n13(cbi)
                     if (atomic(i13(j,cbi)) .eq. 16)  resname = 'CYX'
                  end do
               end if
               if (i.eq.stop .and. ci(i).ne.0) then
                  do j = 1, n12(ci(i))
                     k = i12(j,ci(i))
                     if (atomic(k).eq.8 .and. k.ne.oi(i)) then
                        call pdbatom (' OXT',resname,i,k)
                        goto 20
                     end if
                  end do
   20             continue
               end if
               call getproh (resname,i,m,ni(i),cai(i),cbi)
            end do
         else if (chntyp(m) .eq. 'NUCLEIC') then
            do i = start, stop
               resname = nuclz(seqtyp(i))
               if (resname .eq. ' MP') then
                  call pdbatom (' P  ',resname,i,poi(i))
                  call pdbatom (' OP1',resname,i,op1(i))
                  call pdbatom (' OP2',resname,i,op2(i))
                  call pdbatom (' OP3',resname,i,op3(i))
               else if (resname .eq. ' DP') then
               else if (resname .eq. ' TP') then
               else
                  call pdbatom (' P  ',resname,i,poi(i))
                  call pdbatom (' OP1',resname,i,op1(i))
                  call pdbatom (' OP2',resname,i,op2(i))
                  call pdbatom (' O5''',resname,i,o5i(i))
                  call pdbatom (' C5''',resname,i,c5i(i))
                  call pdbatom (' C4''',resname,i,c4i(i))
                  call pdbatom (' O4''',resname,i,o4i(i))
                  call pdbatom (' C3''',resname,i,c3i(i))
                  call pdbatom (' O3''',resname,i,o3i(i))
                  call pdbatom (' C2''',resname,i,c2i(i))
                  call pdbatom (' O2''',resname,i,o2i(i))
                  call pdbatom (' C1''',resname,i,c1i(i))
                  call getbase (resname,i,ni(i))
                  call getnuch (resname,i,ni(i),c1i(i),c2i(i),o2i(i),
     &                          c3i(i),o3i(i),c4i(i),c5i(i),o5i(i))
               end if
            end do
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (ni)
      deallocate (cai)
      deallocate (ci)
      deallocate (oi)
      deallocate (poi)
      deallocate (op1)
      deallocate (op2)
      deallocate (op3)
      deallocate (c5i)
      deallocate (o5i)
      deallocate (c4i)
      deallocate (o4i)
      deallocate (c3i)
      deallocate (o3i)
      deallocate (c2i)
      deallocate (o2i)
      deallocate (c1i)
      deallocate (restyp)
c
c     perform dynamic allocation of some local arrays
c
      allocate (hetmol(nmol))
c
c     copy any water, ions or ligands following biopolymer chains
c
      if (.not. generic) then
         do i = 1, nmol
            hetmol(i) = .true.
         end do
         do i = 1, n
            if (pdblist(i) .ne. 0)  hetmol(molcule(i)) = .false.
         end do
         do i = 1, nmol
            if (hetmol(i)) then
               do j = imol(1,i), imol(2,i)
                  k = kmol(j)
                  atmnum = atomic(k)
                  atmname = ' '//name(k)
                  justify = 0
                  call numeral (type(k),resname,justify)
                  if (water(i)) then
                     if (atmnum .eq. 1)  atmname = ' H  '
                     if (atmnum .eq. 8)  atmname = ' O  '
                     resname = 'HOH'
                  else if (atmnum .eq. 11) then
                     atmname = 'NA  '
                     resname = ' NA'
                  else if (atmnum .eq. 12) then
                     atmname = 'MG  '
                     resname = ' MG'
                  else if (atmnum .eq. 17) then
                     atmname = 'CL  '
                     resname = ' CL'
                  else if (atmnum .eq. 19) then
                     atmname = ' K  '
                     resname = '  K'
                  else if (atmnum .eq. 20) then
                     atmname = 'CA  '
                     resname = ' CA'
                  else if (atmnum .eq. 35) then
                     atmname = 'BR  '
                     resname = ' BR'
                  else if (atmnum .eq. 53) then
                     atmname = ' I  '
                     resname = '  I'
                  end if
                  pdbnum = nseq + i - 1
                  call pdbatom (atmname,resname,pdbnum,k)
                  pdbtyp(npdb) = 'HETATM'
               end do
            end if
         end do
         do i = 1, nmol
            if (hetmol(i)) then
               do j = imol(1,i), imol(2,i)
                  k = kmol(j)
                  kp = pdblist(k)
                  npdb12(kp) = n12(k)
                  do m = 1, n12(k)
                     ipdb12(m,kp) = pdblist(i12(m,k))
                  end do
               end do
            end if
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (water)
      deallocate (hetmol)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine pdbatom  --  add a single atom to PDB file  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "pdbatom" adds an atom to the Protein Data Bank file
c
c
      subroutine pdbatom (atmname,resname,ires,icoord)
      use sizes
      use atoms
      use pdb
      implicit none
      integer ires,icoord
      character*3 resname
      character*4 atmname
c
c
c     for each atom set the sequential number, record type, atom
c     name, residue name, residue number and atomic coordinates
c
      if (icoord .ne. 0) then
         npdb = npdb + 1
         pdbtyp(npdb) = 'ATOM  '
         pdbatm(npdb) = atmname
         pdbres(npdb) = resname
         resnum(npdb) = ires
         xpdb(npdb) = x(icoord)
         ypdb(npdb) = y(icoord)
         zpdb(npdb) = z(icoord)
         npdb12(npdb) = 0
         pdblist(icoord) = npdb
      end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine getside  --  extract the amino acid side chains  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "getside" finds the side chain heavy atoms for a single amino
c     acid residue and copies the names and coordinates to the Protein
c     Data Bank file
c
c
      subroutine getside (resname,ires,ci,cai,cbi)
      use sizes
      use atomid
      use atoms
      use couple
      implicit none
      integer i,j,ires
      integer ci,cai,cbi
      character*3 resname
c
c
c     if residue is a terminal cap, there is no side chain
c
      cbi = 0
      if (resname .eq. 'H2N')  return
      if (resname .eq. 'FOR')  return
      if (resname .eq. 'ACE')  return
      if (resname .eq. 'COH')  return
      if (resname .eq. 'NH2')  return
      if (resname .eq. 'NME')  return
c
c     find the beta carbon atom for the current residue
c
      do i = 1, n
         if (i.ne.ci .and. atomic(i).eq.6) then
            do j = 1, 4
               if (i12(j,i) .eq. cai) then
                  cbi = i
                  if (resname .ne. 'AIB') then
                     call pdbatom (' CB ',resname,ires,cbi)
                  else
                     call pdbatom (' CB1',resname,ires,cbi)
                  end if
                  goto 10
               end if
            end do
         end if
      end do
   10 continue
      if (cbi .eq. 0)  return
c
c     glycine residue  (GLY)
c
      if (resname .eq. 'GLY') then
         continue
c
c     alanine residue  (ALA)
c
      else if (resname .eq. 'ALA') then
         continue
c
c     valine residue  (VAL)
c
      else if (resname .eq. 'VAL') then
         call pdbatom (' CG1',resname,ires,cbi+1)
         call pdbatom (' CG2',resname,ires,cbi+2)
c
c     leucine residue  (LEU)
c
      else if (resname .eq. 'LEU') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' CD1',resname,ires,cbi+2)
         call pdbatom (' CD2',resname,ires,cbi+3)
c
c     isoleucine residue  (ILE)
c
      else if (resname .eq. 'ILE') then
         call pdbatom (' CG1',resname,ires,cbi+1)
         call pdbatom (' CG2',resname,ires,cbi+2)
         call pdbatom (' CD1',resname,ires,cbi+3)
c
c     serine residue  (SER)
c
      else if (resname .eq. 'SER') then
         call pdbatom (' OG ',resname,ires,cbi+1)
c
c     threonine residue  (THR)
c
      else if (resname .eq. 'THR') then
         call pdbatom (' OG1',resname,ires,cbi+1)
         call pdbatom (' CG2',resname,ires,cbi+2)
c
c     cysteine residue  (CYS)
c
      else if (resname .eq. 'CYS') then
         call pdbatom (' SG ',resname,ires,cbi+1)
c
c     cysteine residue  (CYX)
c
      else if (resname .eq. 'CYX') then
         call pdbatom (' SG ',resname,ires,cbi+1)
c
c     deprotonated cysteine residue  (CYD)
c
      else if (resname .eq. 'CYD') then
         call pdbatom (' SG ',resname,ires,cbi+1)
c
c     proline residue  (PRO)
c
      else if (resname .eq. 'PRO') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' CD ',resname,ires,cbi+2)
c
c     phenylalanine residue  (PHE)
c
      else if (resname .eq. 'PHE') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' CD1',resname,ires,cbi+2)
         call pdbatom (' CD2',resname,ires,cbi+3)
         call pdbatom (' CE1',resname,ires,cbi+4)
         call pdbatom (' CE2',resname,ires,cbi+5)
         call pdbatom (' CZ ',resname,ires,cbi+6)
c
c     tyrosine residue  (TYR)
c
      else if (resname .eq. 'TYR') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' CD1',resname,ires,cbi+2)
         call pdbatom (' CD2',resname,ires,cbi+3)
         call pdbatom (' CE1',resname,ires,cbi+4)
         call pdbatom (' CE2',resname,ires,cbi+5)
         call pdbatom (' CZ ',resname,ires,cbi+6)
         call pdbatom (' OH ',resname,ires,cbi+7)
c
c     deprotonated tyrosine residue  (TYD)
c
      else if (resname .eq. 'TYD') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' CD1',resname,ires,cbi+2)
         call pdbatom (' CD2',resname,ires,cbi+3)
         call pdbatom (' CE1',resname,ires,cbi+4)
         call pdbatom (' CE2',resname,ires,cbi+5)
         call pdbatom (' CZ ',resname,ires,cbi+6)
         call pdbatom (' OH ',resname,ires,cbi+7)
c
c     tryptophan residue  (TRP)
c
      else if (resname .eq. 'TRP') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' CD1',resname,ires,cbi+2)
         call pdbatom (' CD2',resname,ires,cbi+3)
         call pdbatom (' NE1',resname,ires,cbi+4)
         call pdbatom (' CE2',resname,ires,cbi+5)
         call pdbatom (' CE3',resname,ires,cbi+6)
         call pdbatom (' CZ2',resname,ires,cbi+7)
         call pdbatom (' CZ3',resname,ires,cbi+8)
         call pdbatom (' CH2',resname,ires,cbi+9)
c
c     histidine (HD and HE) residue  (HIS)
c
      else if (resname .eq. 'HIS') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' ND1',resname,ires,cbi+2)
         call pdbatom (' CD2',resname,ires,cbi+3)
         call pdbatom (' CE1',resname,ires,cbi+4)
         call pdbatom (' NE2',resname,ires,cbi+5)
c
c     histidine (HD only) residue  (HID)
c
      else if (resname .eq. 'HID') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' ND1',resname,ires,cbi+2)
         call pdbatom (' CD2',resname,ires,cbi+3)
         call pdbatom (' CE1',resname,ires,cbi+4)
         call pdbatom (' NE2',resname,ires,cbi+5)
c
c     histidine (HE only) residue  (HIE)
c
      else if (resname .eq. 'HIE') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' ND1',resname,ires,cbi+2)
         call pdbatom (' CD2',resname,ires,cbi+3)
         call pdbatom (' CE1',resname,ires,cbi+4)
         call pdbatom (' NE2',resname,ires,cbi+5)
c
c     aspartic acid residue  (ASP)
c
      else if (resname .eq. 'ASP') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' OD1',resname,ires,cbi+2)
         call pdbatom (' OD2',resname,ires,cbi+3)
c
c     protonated aspartic acid residue  (ASH)
c
      else if (resname .eq. 'ASH') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' OD1',resname,ires,cbi+2)
         call pdbatom (' OD2',resname,ires,cbi+3)
c
c     asparagine residue  (ASN)
c
      else if (resname .eq. 'ASN') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' OD1',resname,ires,cbi+2)
         call pdbatom (' ND2',resname,ires,cbi+3)
c
c     glutamic acid residue  (GLU)
c
      else if (resname .eq. 'GLU') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' CD ',resname,ires,cbi+2)
         call pdbatom (' OE1',resname,ires,cbi+3)
         call pdbatom (' OE2',resname,ires,cbi+4)
c
c     protonated glutamic acid residue  (GLH)
c
      else if (resname .eq. 'GLH') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' CD ',resname,ires,cbi+2)
         call pdbatom (' OE1',resname,ires,cbi+3)
         call pdbatom (' OE2',resname,ires,cbi+4)
c
c     glutamine residue  (GLN)
c
      else if (resname .eq. 'GLN') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' CD ',resname,ires,cbi+2)
         call pdbatom (' OE1',resname,ires,cbi+3)
         call pdbatom (' NE2',resname,ires,cbi+4)
c
c     methionine residue  (MET)
c
      else if (resname .eq. 'MET') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' SD ',resname,ires,cbi+2)
         call pdbatom (' CE ',resname,ires,cbi+3)
c
c     lysine residue  (LYS)
c
      else if (resname .eq. 'LYS') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' CD ',resname,ires,cbi+2)
         call pdbatom (' CE ',resname,ires,cbi+3)
         call pdbatom (' NZ ',resname,ires,cbi+4)
c
c     deprotonated lysine residue  (LYD)
c
      else if (resname .eq. 'LYD') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' CD ',resname,ires,cbi+2)
         call pdbatom (' CE ',resname,ires,cbi+3)
         call pdbatom (' NZ ',resname,ires,cbi+4)
c
c     arginine residue  (ARG)
c
      else if (resname .eq. 'ARG') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' CD ',resname,ires,cbi+2)
         call pdbatom (' NE ',resname,ires,cbi+3)
         call pdbatom (' CZ ',resname,ires,cbi+4)
         call pdbatom (' NH1',resname,ires,cbi+5)
         call pdbatom (' NH2',resname,ires,cbi+6)
c
c     ornithine residue  (ORN)
c
      else if (resname .eq. 'ORN') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' CD ',resname,ires,cbi+2)
         call pdbatom (' NE ',resname,ires,cbi+3)
c
c     methylalanine residue  (AIB)
c
      else if (resname .eq. 'AIB') then
         call pdbatom (' CB2',resname,ires,cbi+1)
c
c     pyroglutamic acid residue  (PCA)
c
      else if (resname .eq. 'PCA') then
         call pdbatom (' CG ',resname,ires,cbi+1)
         call pdbatom (' CD ',resname,ires,cbi+2)
         call pdbatom (' OE ',resname,ires,cbi+3)
c
c     unknown residue  (UNK)
c
      else if (resname .eq. 'UNK') then
         continue
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine getproh  --  extract the amino acid hydrogens  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "getproh" finds the hydrogen atoms for a single amino acid
c     residue and copies the names and coordinates to the Protein
c     Data Bank file
c
c
      subroutine getproh (resname,ires,jchain,ni,cai,cbi)
      use sizes
      use atomid
      use atoms
      use couple
      use fields
      use sequen
      implicit none
      integer i,nh,hca
      integer ires,jchain
      integer ni,cai,cbi
      logical allatom
      character*3 resname
      character*4 atmname
c
c
c     get any amide hydrogen atoms for non-N-terminal residues
c
      if (ires.ne.ichain(1,jchain) .or. n12(ni).ne.4) then
         if (resname .ne. 'PRO') then
            do i = 1, n
               if (atomic(i).eq.1 .and. i12(1,i).eq.ni) then
                  if (resname .eq. 'COH') then
                     call pdbatom (' HO ',resname,ires,i)
                  else if (resname .eq. 'NH2') then
                     call pdbatom (' H1 ',resname,ires,i)
                     call pdbatom (' H2 ',resname,ires,i+1)
                  else
                     call pdbatom (' H  ',resname,ires,i)
                  end if
                  goto 10
               end if
            end do
         end if
c
c     get any amide hydrogen atoms for N-terminal residues
c
      else
         if (resname .eq. 'PRO') then
            nh = 0
            do i = 1, n
               if (atomic(i).eq.1 .and. i12(1,i).eq.ni) then
                  nh = nh + 1
                  if (nh .eq. 1) then
                     atmname = ' H1 '
                  else if (nh .eq. 2) then
                     atmname = ' H2 '
                  end if
                  call pdbatom (atmname,resname,ires,i)
                  if (nh .eq. 2)  goto 10
               end if
            end do
         else if (resname .eq. 'PCA') then
            do i = 1, n
               if (atomic(i).eq.1 .and. i12(1,i).eq.ni) then
                  atmname = ' H  '
                  call pdbatom (atmname,resname,ires,i)
                  goto 10
               end if
            end do
         else
            nh = 0
            do i = 1, n
               if (atomic(i).eq.1 .and. i12(1,i).eq.ni) then
                  nh = nh + 1
                  if (nh .eq. 1) then
                     atmname = ' H1 '
                  else if (nh .eq. 2) then
                     atmname = ' H2 '
                  else if (nh .eq. 3) then
                     atmname = ' H3 '
                  end if
                  call pdbatom (atmname,resname,ires,i)
                  if (nh .eq. 3)  goto 10
               end if
            end do
         end if
      end if
   10 continue
c
c     get the alpha hydrogen atom for the current residue
c
      hca = 0
      do i = 1, n
         if (atomic(i).eq.1 .and. i12(1,i).eq.cai) then
            hca = i
            if (resname .eq. 'GLY') then
               atmname = ' HA2'
            else if (resname .eq. 'FOR') then
               atmname = ' H  '
            else if (resname .eq. 'ACE') then
               atmname = ' H1 '
            else if (resname .eq. 'NME') then
               atmname = ' H1 '
            else
               atmname = ' HA '
            end if
            call pdbatom (atmname,resname,ires,i)
            goto 20
         end if
      end do
   20 continue
c
c     backbone only if no alpha hydrogen or beta carbon
c
      if (hca.eq.0 .and. cbi.eq.0)  return
c
c     if no alpha hydrogen, then united atom force field
c
      if (hca .ne. 0) then
         allatom = .true.
      else if (resname .eq. 'AIB') then
         if (n12(cbi) .eq. 1) then
            allatom = .false.
         else
            allatom = .true.
         end if
      else
         allatom = .false.
      end if
c
c     glycine residue  (GLY)
c
      if (resname .eq. 'GLY') then
         if (allatom) then
            call pdbatom (' HA3',resname,ires,hca+1)
         end if
c
c     alanine residue  (ALA)
c
      else if (resname .eq. 'ALA') then
         if (allatom) then
            call pdbatom (' HB1',resname,ires,hca+2)
            call pdbatom (' HB2',resname,ires,hca+3)
            call pdbatom (' HB3',resname,ires,hca+4)
         end if
c
c     valine residue  (VAL)
c
      else if (resname .eq. 'VAL') then
         if (allatom) then
            call pdbatom (' HB ',resname,ires,hca+4)
            call pdbatom ('HG11',resname,ires,hca+5)
            call pdbatom ('HG12',resname,ires,hca+6)
            call pdbatom ('HG13',resname,ires,hca+7)
            call pdbatom ('HG21',resname,ires,hca+8)
            call pdbatom ('HG22',resname,ires,hca+9)
            call pdbatom ('HG23',resname,ires,hca+10)
         end if
c
c     leucine residue  (LEU)
c
      else if (resname .eq. 'LEU') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+5)
            call pdbatom (' HB3',resname,ires,hca+6)
            call pdbatom (' HG ',resname,ires,hca+7)
            call pdbatom ('HD11',resname,ires,hca+8)
            call pdbatom ('HD12',resname,ires,hca+9)
            call pdbatom ('HD13',resname,ires,hca+10)
            call pdbatom ('HD21',resname,ires,hca+11)
            call pdbatom ('HD22',resname,ires,hca+12)
            call pdbatom ('HD23',resname,ires,hca+13)
         end if
c
c     isoleucine residue  (ILE)
c
      else if (resname .eq. 'ILE') then
         if (allatom) then
            call pdbatom (' HB ',resname,ires,hca+5)
            call pdbatom ('HG12',resname,ires,hca+6)
            call pdbatom ('HG13',resname,ires,hca+7)
            call pdbatom ('HG21',resname,ires,hca+8)
            call pdbatom ('HG22',resname,ires,hca+9)
            call pdbatom ('HG23',resname,ires,hca+10)
            call pdbatom ('HD11',resname,ires,hca+11)
            call pdbatom ('HD12',resname,ires,hca+12)
            call pdbatom ('HD13',resname,ires,hca+13)
         end if
c
c     serine residue  (SER)
c
      else if (resname .eq. 'SER') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+3)
            call pdbatom (' HB3',resname,ires,hca+4)
            call pdbatom (' HG ',resname,ires,hca+5)
         else
            call pdbatom (' HG ',resname,ires,cbi+2)
         end if
c
c     threonine residue  (THR)
c
      else if (resname .eq. 'THR') then
         if (allatom) then
            call pdbatom (' HB ',resname,ires,hca+4)
            call pdbatom (' HG1',resname,ires,hca+5)
            call pdbatom ('HG21',resname,ires,hca+6)
            call pdbatom ('HG22',resname,ires,hca+7)
            call pdbatom ('HG23',resname,ires,hca+8)
         else
            call pdbatom (' HG1',resname,ires,cbi+3)
         end if
c
c     cysteine residue  (CYS)
c
      else if (resname .eq. 'CYS') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+3)
            call pdbatom (' HB3',resname,ires,hca+4)
            call pdbatom (' HG ',resname,ires,hca+5)
         else if (biotyp(86) .ne. 0) then
            call pdbatom (' HG ',resname,ires,cbi+2)
         end if
c
c     cystine residue  (CYX)
c
      else if (resname .eq. 'CYX') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+3)
            call pdbatom (' HB3',resname,ires,hca+4)
         end if
c
c     deprotonated cysteine residue  (CYD)
c
      else if (resname .eq. 'CYD') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+3)
            call pdbatom (' HB3',resname,ires,hca+4)
         end if
c
c     proline residue  (PRO)
c
      else if (resname .eq. 'PRO') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+4)
            call pdbatom (' HB3',resname,ires,hca+5)
            call pdbatom (' HG2',resname,ires,hca+6)
            call pdbatom (' HG3',resname,ires,hca+7)
            call pdbatom (' HD2',resname,ires,hca+8)
            call pdbatom (' HD3',resname,ires,hca+9)
         end if
c
c     phenylalanine residue  (PHE)
c
      else if (resname .eq. 'PHE') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+8)
            call pdbatom (' HB3',resname,ires,hca+9)
            call pdbatom (' HD1',resname,ires,hca+10)
            call pdbatom (' HD2',resname,ires,hca+11)
            call pdbatom (' HE1',resname,ires,hca+12)
            call pdbatom (' HE2',resname,ires,hca+13)
            call pdbatom (' HZ ',resname,ires,hca+14)
         else if (biotyp(126) .ne. 0) then
            call pdbatom (' HD1',resname,ires,cbi+7)
            call pdbatom (' HD2',resname,ires,cbi+8)
            call pdbatom (' HE1',resname,ires,cbi+9)
            call pdbatom (' HE2',resname,ires,cbi+10)
            call pdbatom (' HZ ',resname,ires,cbi+11)
         end if
c
c     tyrosine residue  (TYR)
c
      else if (resname .eq. 'TYR') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+9)
            call pdbatom (' HB3',resname,ires,hca+10)
            call pdbatom (' HD1',resname,ires,hca+11)
            call pdbatom (' HD2',resname,ires,hca+12)
            call pdbatom (' HE1',resname,ires,hca+13)
            call pdbatom (' HE2',resname,ires,hca+14)
            call pdbatom (' HH ',resname,ires,hca+15)
         else if (biotyp(141) .ne. 0) then
            call pdbatom (' HD1',resname,ires,cbi+8)
            call pdbatom (' HD2',resname,ires,cbi+9)
            call pdbatom (' HE1',resname,ires,cbi+10)
            call pdbatom (' HE2',resname,ires,cbi+11)
            call pdbatom (' HH ',resname,ires,cbi+12)
         else
            call pdbatom (' HH ',resname,ires,cbi+8)
         end if
c
c     deprotonated tyrosine residue  (TYD)
c
      else if (resname .eq. 'TYD') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+9)
            call pdbatom (' HB3',resname,ires,hca+10)
            call pdbatom (' HD1',resname,ires,hca+11)
            call pdbatom (' HD2',resname,ires,hca+12)
            call pdbatom (' HE1',resname,ires,hca+13)
            call pdbatom (' HE2',resname,ires,hca+14)
         else if (biotyp(141) .ne. 0) then
            call pdbatom (' HD1',resname,ires,cbi+8)
            call pdbatom (' HD2',resname,ires,cbi+9)
            call pdbatom (' HE1',resname,ires,cbi+10)
            call pdbatom (' HE2',resname,ires,cbi+11)
         end if
c
c     tryptophan residue  (TRP)
c
      else if (resname .eq. 'TRP') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+11)
            call pdbatom (' HB3',resname,ires,hca+12)
            call pdbatom (' HD1',resname,ires,hca+13)
            call pdbatom (' HE1',resname,ires,hca+14)
            call pdbatom (' HE3',resname,ires,hca+15)
            call pdbatom (' HZ2',resname,ires,hca+16)
            call pdbatom (' HZ3',resname,ires,hca+17)
            call pdbatom (' HH2',resname,ires,hca+18)
         else if (biotyp(172) .ne. 0) then
            call pdbatom (' HD1',resname,ires,cbi+10)
            call pdbatom (' HE1',resname,ires,cbi+11)
            call pdbatom (' HE3',resname,ires,cbi+12)
            call pdbatom (' HZ2',resname,ires,cbi+13)
            call pdbatom (' HZ3',resname,ires,cbi+14)
            call pdbatom (' HH2',resname,ires,cbi+15)
         else
            call pdbatom (' HE1',resname,ires,cbi+10)
         end if
c
c     histidine (HD and HE) residue  (HIS)
c
      else if (resname .eq. 'HIS') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+7)
            call pdbatom (' HB3',resname,ires,hca+8)
            call pdbatom (' HD1',resname,ires,hca+9)
            call pdbatom (' HD2',resname,ires,hca+10)
            call pdbatom (' HE1',resname,ires,hca+11)
            call pdbatom (' HE2',resname,ires,hca+12)
         else if (biotyp(197) .ne. 0) then
            call pdbatom (' HD1',resname,ires,cbi+6)
            call pdbatom (' HD2',resname,ires,cbi+7)
            call pdbatom (' HE1',resname,ires,cbi+8)
            call pdbatom (' HE2',resname,ires,cbi+9)
         else
            call pdbatom (' HD1',resname,ires,cbi+6)
            call pdbatom (' HE2',resname,ires,cbi+7)
         end if
c
c     histidine (HD only) residue  (HID)
c
      else if (resname .eq. 'HID') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+7)
            call pdbatom (' HB3',resname,ires,hca+8)
            call pdbatom (' HD1',resname,ires,hca+9)
            call pdbatom (' HD2',resname,ires,hca+10)
            call pdbatom (' HE1',resname,ires,hca+11)
         else if (biotyp(214) .ne. 0) then
            call pdbatom (' HD1',resname,ires,cbi+6)
            call pdbatom (' HD2',resname,ires,cbi+7)
            call pdbatom (' HE1',resname,ires,cbi+8)
         else
            call pdbatom (' HD1',resname,ires,cbi+6)
         end if
c
c     histidine (HE only) residue  (HIE)
c
      else if (resname .eq. 'HIE') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+7)
            call pdbatom (' HB3',resname,ires,hca+8)
            call pdbatom (' HD2',resname,ires,hca+9)
            call pdbatom (' HE1',resname,ires,hca+10)
            call pdbatom (' HE2',resname,ires,hca+11)
         else if (biotyp(229) .ne. 0) then
            call pdbatom (' HD2',resname,ires,cbi+6)
            call pdbatom (' HE1',resname,ires,cbi+7)
            call pdbatom (' HE2',resname,ires,cbi+8)
         else
            call pdbatom (' HE2',resname,ires,cbi+6)
         end if
c
c     aspartic acid residue  (ASP)
c
      else if (resname .eq. 'ASP') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+5)
            call pdbatom (' HB3',resname,ires,hca+6)
         end if
c
c     protonated aspartic acid residue  (ASH)
c
      else if (resname .eq. 'ASH') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+5)
            call pdbatom (' HB3',resname,ires,hca+6)
            call pdbatom (' HD2',resname,ires,hca+7)
         else
            call pdbatom (' HD2',resname,ires,cbi+4)
         end if
c
c     asparagine residue  (ASN)
c
      else if (resname .eq. 'ASN') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+5)
            call pdbatom (' HB3',resname,ires,hca+6)
            call pdbatom ('HD21',resname,ires,hca+7)
            call pdbatom ('HD22',resname,ires,hca+8)
         else
            call pdbatom ('HD21',resname,ires,cbi+4)
            call pdbatom ('HD22',resname,ires,cbi+5)
         end if
c
c     glutamic acid residue  (GLU)
c
      else if (resname .eq. 'GLU') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+6)
            call pdbatom (' HB3',resname,ires,hca+7)
            call pdbatom (' HG2',resname,ires,hca+8)
            call pdbatom (' HG3',resname,ires,hca+9)
         end if
c
c     protonated glutamic acid residue  (GLH)
c
      else if (resname .eq. 'GLH') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+6)
            call pdbatom (' HB3',resname,ires,hca+7)
            call pdbatom (' HG2',resname,ires,hca+8)
            call pdbatom (' HG3',resname,ires,hca+9)
            call pdbatom (' HE2',resname,ires,hca+10)
         else
            call pdbatom (' HE2',resname,ires,cbi+5)
         end if
c
c     glutamine residue  (GLN)
c
      else if (resname .eq. 'GLN') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+6)
            call pdbatom (' HB3',resname,ires,hca+7)
            call pdbatom (' HG2',resname,ires,hca+8)
            call pdbatom (' HG3',resname,ires,hca+9)
            call pdbatom ('HE21',resname,ires,hca+10)
            call pdbatom ('HE22',resname,ires,hca+11)
         else
            call pdbatom ('HE21',resname,ires,cbi+5)
            call pdbatom ('HE22',resname,ires,cbi+6)
         end if
c
c     methionine residue  (MET)
c
      else if (resname .eq. 'MET') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+5)
            call pdbatom (' HB3',resname,ires,hca+6)
            call pdbatom (' HG2',resname,ires,hca+7)
            call pdbatom (' HG3',resname,ires,hca+8)
            call pdbatom (' HE1',resname,ires,hca+9)
            call pdbatom (' HE2',resname,ires,hca+10)
            call pdbatom (' HE3',resname,ires,hca+11)
         end if
c
c     lysine residue  (LYS)
c
      else if (resname .eq. 'LYS') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+6)
            call pdbatom (' HB3',resname,ires,hca+7)
            call pdbatom (' HG2',resname,ires,hca+8)
            call pdbatom (' HG3',resname,ires,hca+9)
            call pdbatom (' HD2',resname,ires,hca+10)
            call pdbatom (' HD3',resname,ires,hca+11)
            call pdbatom (' HE2',resname,ires,hca+12)
            call pdbatom (' HE3',resname,ires,hca+13)
            call pdbatom (' HZ1',resname,ires,hca+14)
            call pdbatom (' HZ2',resname,ires,hca+15)
            call pdbatom (' HZ3',resname,ires,hca+16)
         else
            call pdbatom (' HZ1',resname,ires,cbi+5)
            call pdbatom (' HZ2',resname,ires,cbi+6)
            call pdbatom (' HZ3',resname,ires,cbi+7)
         end if
c
c     deprotonated lysine residue  (LYD)
c
      else if (resname .eq. 'LYD') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+6)
            call pdbatom (' HB3',resname,ires,hca+7)
            call pdbatom (' HG2',resname,ires,hca+8)
            call pdbatom (' HG3',resname,ires,hca+9)
            call pdbatom (' HD2',resname,ires,hca+10)
            call pdbatom (' HD3',resname,ires,hca+11)
            call pdbatom (' HE2',resname,ires,hca+12)
            call pdbatom (' HE3',resname,ires,hca+13)
            call pdbatom (' HZ1',resname,ires,hca+14)
            call pdbatom (' HZ2',resname,ires,hca+15)
         else
            call pdbatom (' HZ1',resname,ires,cbi+5)
            call pdbatom (' HZ2',resname,ires,cbi+6)
         end if
c
c     arginine residue  (ARG)
c
      else if (resname .eq. 'ARG') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+8)
            call pdbatom (' HB3',resname,ires,hca+9)
            call pdbatom (' HG2',resname,ires,hca+10)
            call pdbatom (' HG3',resname,ires,hca+11)
            call pdbatom (' HD2',resname,ires,hca+12)
            call pdbatom (' HD3',resname,ires,hca+13)
            call pdbatom (' HE ',resname,ires,hca+14)
            call pdbatom ('HH11',resname,ires,hca+15)
            call pdbatom ('HH12',resname,ires,hca+16)
            call pdbatom ('HH21',resname,ires,hca+17)
            call pdbatom ('HH22',resname,ires,hca+18)
         else
            call pdbatom (' HE ',resname,ires,cbi+7)
            call pdbatom ('HH11',resname,ires,cbi+8)
            call pdbatom ('HH12',resname,ires,cbi+9)
            call pdbatom ('HH21',resname,ires,cbi+10)
            call pdbatom ('HH22',resname,ires,cbi+11)
         end if
c
c     ornithine residue  (ORN)
c
      else if (resname .eq. 'ORN') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+5)
            call pdbatom (' HB3',resname,ires,hca+6)
            call pdbatom (' HG2',resname,ires,hca+7)
            call pdbatom (' HG3',resname,ires,hca+8)
            call pdbatom (' HD2',resname,ires,hca+9)
            call pdbatom (' HD3',resname,ires,hca+10)
            call pdbatom (' HE1',resname,ires,hca+11)
            call pdbatom (' HE2',resname,ires,hca+12)
            call pdbatom (' HE3',resname,ires,hca+13)
         else
            call pdbatom (' HE1',resname,ires,cbi+4)
            call pdbatom (' HE2',resname,ires,cbi+5)
            call pdbatom (' HE3',resname,ires,cbi+6)
         end if
c
c     methylalanine residue  (AIB)
c
      else if (resname .eq. 'AIB') then
         if (allatom) then
            call pdbatom ('HB11',resname,ires,cbi+2)
            call pdbatom ('HB12',resname,ires,cbi+3)
            call pdbatom ('HB13',resname,ires,cbi+4)
            call pdbatom ('HB21',resname,ires,cbi+5)
            call pdbatom ('HB22',resname,ires,cbi+6)
            call pdbatom ('HB23',resname,ires,cbi+7)
         end if
c
c     pyroglutamic acid residue  (PCA)
c
      else if (resname .eq. 'PCA') then
         if (allatom) then
            call pdbatom (' HB2',resname,ires,hca+5)
            call pdbatom (' HB3',resname,ires,hca+6)
            call pdbatom (' HG2',resname,ires,hca+7)
            call pdbatom (' HG3',resname,ires,hca+8)
         end if
c
c     unknown residue  (UNK)
c
      else if (resname .eq. 'UNK') then
         if (allatom) then
            call pdbatom (' HA3',resname,ires,hca+1)
         end if
c
c     N-terminal deprotonated residue  (H2N)
c
      else if (resname .eq. 'H2N') then
         continue
c
c     N-terminal formyl residue  (FOR)
c
      else if (resname .eq. 'FOR') then
         continue
c
c     N-terminal acetyl residue  (ACE)
c
      else if (resname .eq. 'ACE') then
         if (allatom) then
            call pdbatom (' H2 ',resname,ires,hca+1)
            call pdbatom (' H3 ',resname,ires,hca+2)
         end if
c
c     C-terminal protonated residue (COH)
c
      else if (resname .eq. 'COH') then
         continue
c
c     C-terminal amide residue  (NH2)
c
      else if (resname .eq. 'NH2') then
         continue
c
c     C-terminal N-methylamide residue  (NME)
c
      else if (resname .eq. 'NME') then
         if (allatom) then
            call pdbatom (' H2 ',resname,ires,hca+1)
            call pdbatom (' H3 ',resname,ires,hca+2)
         end if
      end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine getbase  --  extract the nucleotide side chains  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "getbase" finds the base heavy atoms for a single nucleotide
c     residue and copies the names and coordinates to the Protein
c     Data Bank file
c
c
      subroutine getbase (resname,ires,ni)
      implicit none
      integer ires,ni
      character*3 resname
c
c
c     adenine in adenosine residue  (A)
c
      if (resname .eq. '  A') then
         call pdbatom (' N9 ',resname,ires,ni)
         call pdbatom (' C8 ',resname,ires,ni+1)
         call pdbatom (' N7 ',resname,ires,ni+2)
         call pdbatom (' C5 ',resname,ires,ni+3)
         call pdbatom (' C6 ',resname,ires,ni+4)
         call pdbatom (' N6 ',resname,ires,ni+5)
         call pdbatom (' N1 ',resname,ires,ni+6)
         call pdbatom (' C2 ',resname,ires,ni+7)
         call pdbatom (' N3 ',resname,ires,ni+8)
         call pdbatom (' C4 ',resname,ires,ni+9)
c
c     guanine in guanosine residue  (G)
c
      else if (resname .eq. '  G') then
         call pdbatom (' N9 ',resname,ires,ni)
         call pdbatom (' C8 ',resname,ires,ni+1)
         call pdbatom (' N7 ',resname,ires,ni+2)
         call pdbatom (' C5 ',resname,ires,ni+3)
         call pdbatom (' C6 ',resname,ires,ni+4)
         call pdbatom (' O6 ',resname,ires,ni+5)
         call pdbatom (' N1 ',resname,ires,ni+6)
         call pdbatom (' C2 ',resname,ires,ni+7)
         call pdbatom (' N2 ',resname,ires,ni+8)
         call pdbatom (' N3 ',resname,ires,ni+9)
         call pdbatom (' C4 ',resname,ires,ni+10)
c
c     cytosine in cytidine residue  (C)
c
      else if (resname .eq. '  C') then
         call pdbatom (' N1 ',resname,ires,ni)
         call pdbatom (' C2 ',resname,ires,ni+1)
         call pdbatom (' O2 ',resname,ires,ni+2)
         call pdbatom (' N3 ',resname,ires,ni+3)
         call pdbatom (' C4 ',resname,ires,ni+4)
         call pdbatom (' N4 ',resname,ires,ni+5)
         call pdbatom (' C5 ',resname,ires,ni+6)
         call pdbatom (' C6 ',resname,ires,ni+7)
c
c     uracil in uridine residue  (U)
c
      else if (resname .eq. '  U') then
         call pdbatom (' N1 ',resname,ires,ni)
         call pdbatom (' C2 ',resname,ires,ni+1)
         call pdbatom (' O2 ',resname,ires,ni+2)
         call pdbatom (' N3 ',resname,ires,ni+3)
         call pdbatom (' C4 ',resname,ires,ni+4)
         call pdbatom (' O4 ',resname,ires,ni+5)
         call pdbatom (' C5 ',resname,ires,ni+6)
         call pdbatom (' C6 ',resname,ires,ni+7)
c
c     adenine in deoxyadenosine residue  (DA)
c
      else if (resname .eq. ' DA') then
         call pdbatom (' N9 ',resname,ires,ni)
         call pdbatom (' C8 ',resname,ires,ni+1)
         call pdbatom (' N7 ',resname,ires,ni+2)
         call pdbatom (' C5 ',resname,ires,ni+3)
         call pdbatom (' C6 ',resname,ires,ni+4)
         call pdbatom (' N6 ',resname,ires,ni+5)
         call pdbatom (' N1 ',resname,ires,ni+6)
         call pdbatom (' C2 ',resname,ires,ni+7)
         call pdbatom (' N3 ',resname,ires,ni+8)
         call pdbatom (' C4 ',resname,ires,ni+9)
c
c     guanine in deoxyguanosine residue  (DG)
c
      else if (resname .eq. ' DG') then
         call pdbatom (' N9 ',resname,ires,ni)
         call pdbatom (' C8 ',resname,ires,ni+1)
         call pdbatom (' N7 ',resname,ires,ni+2)
         call pdbatom (' C5 ',resname,ires,ni+3)
         call pdbatom (' C6 ',resname,ires,ni+4)
         call pdbatom (' O6 ',resname,ires,ni+5)
         call pdbatom (' N1 ',resname,ires,ni+6)
         call pdbatom (' C2 ',resname,ires,ni+7)
         call pdbatom (' N2 ',resname,ires,ni+8)
         call pdbatom (' N3 ',resname,ires,ni+9)
         call pdbatom (' C4 ',resname,ires,ni+10)
c
c     cytosine in deoxycytidine residue  (DC)
c
      else if (resname .eq. ' DC') then
         call pdbatom (' N1 ',resname,ires,ni)
         call pdbatom (' C2 ',resname,ires,ni+1)
         call pdbatom (' O2 ',resname,ires,ni+2)
         call pdbatom (' N3 ',resname,ires,ni+3)
         call pdbatom (' C4 ',resname,ires,ni+4)
         call pdbatom (' N4 ',resname,ires,ni+5)
         call pdbatom (' C5 ',resname,ires,ni+6)
         call pdbatom (' C6 ',resname,ires,ni+7)
c
c     thymine in deoxythymidine residue  (DT)
c
      else if (resname .eq. ' DT') then
         call pdbatom (' N1 ',resname,ires,ni)
         call pdbatom (' C2 ',resname,ires,ni+1)
         call pdbatom (' O2 ',resname,ires,ni+2)
         call pdbatom (' N3 ',resname,ires,ni+3)
         call pdbatom (' C4 ',resname,ires,ni+4)
         call pdbatom (' O4 ',resname,ires,ni+5)
         call pdbatom (' C5 ',resname,ires,ni+6)
         call pdbatom (' C7 ',resname,ires,ni+7)
         call pdbatom (' C6 ',resname,ires,ni+8)
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine getnuch  --  extract the nucleotide hydrogens  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "getnuch" finds the nucleotide hydrogen atoms for a single
c     residue and copies the names and coordinates to the Protein
c     Data Bank file
c
c
      subroutine getnuch (resname,ires,ni,c1i,c2i,
     &                    o2i,c3i,o3i,c4i,c5i,o5i)
      use sizes
      use atomid
      use couple
      implicit none
      integer i,k
      integer ires,ni
      integer c1i,c4i
      integer c2i,o2i
      integer c3i,o3i
      integer c5i,o5i
      logical allatom,done
      character*3 resname
c
c
c     if no ribose C1 hydrogen, then united atom force field
c
      allatom = .true.
      if (n12(c1i) .ne. 4)  allatom = .false.
c
c     get sugar ring hydrogen atoms for the current residue
c
      done = .false.
      do i = 1, n12(c5i)
         k = i12(i,c5i)
         if (atomic(k).eq.1) then
            if (.not. done) then
               call pdbatom (' H5''',resname,ires,k)
               done = .true.
            else
               call pdbatom ('H5''''',resname,ires,k)
            end if
         end if
      end do
      do i = 1, n12(c4i)
         k = i12(i,c4i)
         if (atomic(k) .eq. 1)  call pdbatom (' H4''',resname,ires,k)
      end do
      do i = 1, n12(c3i)
         k = i12(i,c3i)
         if (atomic(k) .eq. 1)  call pdbatom (' H3''',resname,ires,k)
      end do
      done = .false.
      do i = 1, n12(c2i)
         k = i12(i,c2i)
         if (atomic(k) .eq. 1) then
            if (.not. done) then
               call pdbatom (' H2''',resname,ires,k)
               done = .true.
            else
               call pdbatom ('H2''''',resname,ires,k)
            end if
         end if
      end do
      if (o2i .ne. 0) then
         do i = 1, n12(o2i)
            k = i12(i,o2i)
            if (atomic(k) .eq. 1)  call pdbatom ('HO2''',resname,ires,k)
         end do
      end if
      do i = 1, n12(c1i)
         k = i12(i,c1i)
         if (atomic(k) .eq. 1)  call pdbatom (' H1''',resname,ires,k)
      end do
c
c     adenine in adenosine residue  (A)
c
      if (resname .eq. '  A') then
         if (allatom) then
            call pdbatom (' H8 ',resname,ires,ni+10)
            call pdbatom (' H61',resname,ires,ni+11)
            call pdbatom (' H62',resname,ires,ni+12)
            call pdbatom (' H2 ',resname,ires,ni+13)
         else
            call pdbatom (' H61',resname,ires,ni+10)
            call pdbatom (' H62',resname,ires,ni+11)
         end if
c
c     guanine in guanosine residue  (G)
c
      else if (resname .eq. '  G') then
         if (allatom) then
            call pdbatom (' H8 ',resname,ires,ni+11)
            call pdbatom (' H1 ',resname,ires,ni+12)
            call pdbatom (' H21',resname,ires,ni+13)
            call pdbatom (' H22',resname,ires,ni+14)
         else
            call pdbatom (' H1 ',resname,ires,ni+11)
            call pdbatom (' H21',resname,ires,ni+12)
            call pdbatom (' H22',resname,ires,ni+13)
         end if
c
c     cytosine in cytidine residue  (C)
c
      else if (resname .eq. '  C') then
         if (allatom) then
            call pdbatom (' H41',resname,ires,ni+8)
            call pdbatom (' H42',resname,ires,ni+9)
            call pdbatom (' H5 ',resname,ires,ni+10)
            call pdbatom (' H6 ',resname,ires,ni+11)
         else
            call pdbatom (' H41',resname,ires,ni+8)
            call pdbatom (' H42',resname,ires,ni+9)
         end if
c
c     uracil in uridine residue  (U)
c
      else if (resname .eq. '  U') then
         if (allatom) then
            call pdbatom (' H3 ',resname,ires,ni+8)
            call pdbatom (' H5 ',resname,ires,ni+9)
            call pdbatom (' H6 ',resname,ires,ni+10)
         else
            call pdbatom (' H3 ',resname,ires,ni+8)
         end if
c
c     adenine in deoxyadenosine residue  (DA)
c
      else if (resname .eq. ' DA') then
         if (allatom) then
            call pdbatom (' H8 ',resname,ires,ni+10)
            call pdbatom (' H61',resname,ires,ni+11)
            call pdbatom (' H62',resname,ires,ni+12)
            call pdbatom (' H2 ',resname,ires,ni+13)
         else
            call pdbatom (' H61',resname,ires,ni+10)
            call pdbatom (' H62',resname,ires,ni+11)
         end if
c
c     guanine in deoxyguanosine residue  (DG)
c
      else if (resname .eq. ' DG') then
         if (allatom) then
            call pdbatom (' H8 ',resname,ires,ni+11)
            call pdbatom (' H1 ',resname,ires,ni+12)
            call pdbatom (' H21',resname,ires,ni+13)
            call pdbatom (' H22',resname,ires,ni+14)
         else
            call pdbatom (' H1 ',resname,ires,ni+11)
            call pdbatom (' H21',resname,ires,ni+12)
            call pdbatom (' H22',resname,ires,ni+13)
         end if
c
c     cytosine in deoxycytidine residue  (DC)
c
      else if (resname .eq. ' DC') then
         if (allatom) then
            call pdbatom (' H41',resname,ires,ni+8)
            call pdbatom (' H42',resname,ires,ni+9)
            call pdbatom (' H5 ',resname,ires,ni+10)
            call pdbatom (' H6 ',resname,ires,ni+11)
         else
            call pdbatom (' H41',resname,ires,ni+8)
            call pdbatom (' H42',resname,ires,ni+9)
         end if
c
c     thymine in deoxythymidine residue  (DT)
c
      else if (resname .eq. ' DT') then
         if (allatom) then
            call pdbatom (' H3 ',resname,ires,ni+9)
            call pdbatom (' H71',resname,ires,ni+10)
            call pdbatom (' H72',resname,ires,ni+11)
            call pdbatom (' H73',resname,ires,ni+12)
            call pdbatom (' H6 ',resname,ires,ni+13)
         else
            call pdbatom (' H3 ',resname,ires,ni+9)
         end if
      end if
c
c     get any capping hydrogen atoms for the current residue
c
      do i = 1, n12(o5i)
         k = i12(i,o5i)
         if (atomic(k) .eq. 1)  call pdbatom (' H5T',resname,ires,k)
      end do
      do i = 1, n12(o3i)
         k = i12(i,o3i)
         if (atomic(k) .eq. 1)  call pdbatom (' H3T',resname,ires,k)
      end do
      return
      end
