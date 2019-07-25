c
c
c     #############################################################
c     ##                  COPYRIGHT (C) 1999 by                  ##
c     ##  Marina A. Vorobieva, Nina N. Sokolova & Jay W. Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  program nucleic  --  build a nucleic acid from sequence  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "nucleic" builds the internal and Cartesian coordinates
c     of a polynucleotide from nucleic acid sequence and torsional
c     angle values for the nucleic acid backbone and side chains
c
c
      program nucleic
      use sizes
      use atoms
      use couple
      use files
      use iounit
      use nucleo
      use titles
      implicit none
      integer i,natom,mode
      integer izmt,ixyz,iseq
      integer freeunit,trimtext
      logical exist
      character*240 seqfile
      character*240 intfile
      character*240 xyzfile
c
c
c     get the name to use for the output structure files
c
      call initial
      call nextarg (filename,exist)
      if (.not. exist) then
         write (iout,10)
   10    format (/,' Enter Name to be Used for Output Files :  ',$)
         read (input,20)  filename
   20    format (a240)
      end if
      call basefile (filename)
c
c     get the title line for the output files
c
      write (iout,30)
   30 format (/,' Enter Title :  ',$)
      read (input,40)  title
   40 format (a240)
      ltitle = trimtext (title)
c
c     read the keyfile and force field parameter file
c
      call getkey
      call field
c
c     get the sequence, build a Z-matrix, convert to Cartesians
c
      call getseqn
      call nucchain
      call connect
      call molecule
      call makexyz
c
c     perform the alignment of the strands of a double helix
c
      if (dblhlx) then
         call watson
         call inertia (2)
      end if
c
c     remove dummy atoms and set undefined atoms to type zero
c
      natom = n
      do i = natom, 1, -1
         if (type(i) .eq. 0)  call delete (i)
         if (type(i) .lt. 0)  type(i) = 0
      end do
c
c     convert to internal and Cartesian coordinates
c
      mode = 0
      call makeint (mode)
      call makexyz
c
c     write out a nucleic acid sequence file
c
      iseq = freeunit ()
      seqfile = filename(1:leng)//'.seq'
      call version (seqfile,'new')
      open (unit=iseq,file=seqfile,status='new')
      call prtseq (iseq)
      close (unit=iseq)
c
c     write out an internal coordinates file
c
      izmt = freeunit ()
      intfile = filename(1:leng)//'.int'
      call version (intfile,'new')
      open (unit=izmt,file=intfile,status='new')
      call prtint (izmt)
      close (unit=izmt)
c
c     write out a Cartesian coordinates file
c
      ixyz = freeunit ()
      xyzfile = filename(1:leng)//'.xyz'
      call version (xyzfile,'new')
      open (unit=ixyz,file=xyzfile,status='new')
      call prtxyz (ixyz)
      close (unit=ixyz)
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine getseqn  --  nucleic acid sequence and angles  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "getseqn" asks the user for the nucleotide sequence and
c     torsional angle values needed to define a nucleic acid
c
c
      subroutine getseqn
      use iounit
      use nucleo
      use resdue
      use sequen
      implicit none
      integer i,j,k,next
      integer start,stop
      integer length,trimtext
      logical done
      logical, allocatable :: purine(:)
      character*1 answer
      character*1 ucase(26)
      character*3 name,resname
      character*240 record
      character*240 string
      data ucase  / 'A','B','C','D','E','F','G','H','I','J','K','L',
     &              'M','N','O','P','Q','R','S','T','U','V','W','X',
     &              'Y','Z' /
c
c
c     choose to generate either an A-, B- or Z-form helix
c
      write (iout,10)
   10 format (/,' Enter A-, B- or Z-Form Helix for the Structure',
     &           ' [B] :  ',$)
      read (input,20)  record
   20 format (a240)
      call upcase (record)
      next = 1
      call getword (record,answer,next)
      hlxform = 'B'
      if (answer .eq. 'A')  hlxform = 'A'
      if (answer .eq. 'Z')  hlxform = 'Z'
c
c     provide a header to explain the method of sequence input
c
      write (iout,30)
   30 format (/,' Enter One Nucleotide per Line, 5'' to 3'': ',
     &           ' Give PDB Residue Code,',
     &        /,' followed by Backbone Torsions (6F) and',
     &           ' Glycosidic Torsion (1F)',
     &        //,' Use Residue=MOL to Begin a New Strand,',
     &           ' Residue=<CR> to End Input')
c
c     initially, assume that only a single strand is present
c
      nchain = 1
      ichain(1,1) = 1
      chnnam(1) = ' '
c
c     get the nucleotide sequence data and dihedral angle values
c
      i = 0
      done = .false.
      do while (.not. done)
         i = i + 1
         do j = 1, 6
            bkbone(j,i) = 0.0d0
         end do
         glyco(i) = 0.0d0
         pucker(i) = 0
         write (iout,40)  i
   40    format (/,' Enter Residue',i4,' :  ',$)
         read (input,50)  record
   50    format (a240)
         call upcase (record)
         next = 1
         call gettext (record,name,next)
         length = trimtext (name)
         string = record(next:240)
         read (string,*,err=60,end=60)  (bkbone(j,i),j=1,6),glyco(i)
   60    continue
c
c     process and store the current nucleotide type
c
         if (name .eq. 'MOL') then
            i = i - 1
            ichain(2,nchain) = i
            nchain = nchain + 1
            ichain(1,nchain) = i + 1
         else
            if (name .eq. '   ') then
               done = .true.
               nseq = i - 1
               ichain(2,nchain) = nseq
            else
               seq(i) = nuclz(maxnuc)
               seqtyp(i) = 0
               if (length .eq. 1) then
                  do j = 1, maxnuc
                     if (name(1:1) .eq. nuclz1(j)) then
                        seq(i) = nuclz(j)
                        seqtyp(i) = j
                     end if
                  end do
               else
                  do j = 1, maxnuc
                     if (name .eq. nuclz(j)) then
                        seq(i) = nuclz(j)
                        seqtyp(i) = j
                     end if
                  end do
               end if
               if (seqtyp(i) .eq. 0) then
                  i = i - 1
                  write (iout,70)  name
   70             format (/,' GETSEQN  --  Nucleotide Type ',a3,
     &                       ' is Not Supported')
               end if
            end if
         end if
      end do
c
c     offer the option to construct an idealized double helix
c
      dblhlx = .false.
      if (nchain .eq. 1) then
         write (iout,80)
   80    format (/,' Build a Double Helix using Complimentary Bases',
     &              ' [N] :  ',$)
         read (input,90)  record
   90    format (a240)
         next = 1
         call gettext (record,answer,next)
         call upcase (answer)
         if (answer .eq. 'Y')  dblhlx = .true.
      else if (nchain .eq. 2) then
         write (iout,100)
  100    format (/,' Combine the Two Single Strands into Double Helix',
     &              ' [Y] :  ',$)
         read (input,110)  record
  110    format (a240)
         next = 1
         call gettext (record,answer,next)
         call upcase (answer)
         if (answer .ne. 'N')  dblhlx = .true.
      end if
c
c     build a second strand as the reverse-complement sequence
c
      if (nchain.eq.1 .and. dblhlx) then
         start = 1
         stop = nseq
         resname = nuclz(seqtyp(1))
         if (resname.eq.' MP' .or. resname.eq.' DP'
     &          .or. resname.eq.' TP') then
            k = nseq + 1
            seq(k) = seq(1)
            seqtyp(k) = seqtyp(1)
            start = 2
         end if
         resname = nuclz(seqtyp(nseq))
         if (resname.eq.' MP' .or. resname.eq.' DP'
     &          .or. resname.eq.' TP') then
            k = 2 * nseq
            seq(k) = seq(nseq)
            seqtyp(k) = seqtyp(nseq)
            stop = nseq - 1
         end if
         do i = start, stop
            resname = nuclz(seqtyp(i))
            if (resname .eq. '  A') then
               resname = '  U'
            else if (resname .eq. '  G') then
               resname = '  C'
            else if (resname .eq. '  C') then
               resname = '  G'
            else if (resname .eq. '  U') then
               resname = '  A'
            else if (resname .eq. ' DA') then
               resname = ' DT'
            else if (resname .eq. ' DG') then
               resname = ' DC'
            else if (resname .eq. ' DC') then
               resname = ' DG'
            else if (resname .eq. ' DT') then
               resname = ' DA'
            end if
            k = nseq + stop + start - i
            do j = 1, maxnuc
               if (resname .eq. nuclz(j)) then
                  seq(k) = nuclz(j)
                  seqtyp(k) = j
               end if
            end do
         end do
         do i = 1, nseq
            k = nseq + i
            do j = 1, 6
               bkbone(j,k) = bkbone(j,i)
            end do
            glyco(k) = glyco(i)
            pucker(k) = pucker(i)
         end do
         nchain = 2
         nseq = 2 * nseq
         ichain(1,nchain) = nseq/2 + 1
         ichain(2,nchain) = nseq
      end if
c
c     set chain identifiers if multiple chains are present
c
      if (nchain .gt. 1) then
         do i = 1, nchain
            chnnam(i) = ucase(i)
         end do
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (purine(nseq))
c
c     set the nucleic acid base and sugar structural type
c
      do i = 1, nseq
         resname = nuclz(seqtyp(i))
         purine(i) = .false.
         if (resname .eq. '  A')  purine(i) = .true.
         if (resname .eq. '  G')  purine(i) = .true.
         if (resname .eq. ' DA')  purine(i) = .true.
         if (resname .eq. ' DG')  purine(i) = .true.
         deoxy(i) = .false.
         if (resname .eq. ' DA')  deoxy(i) = .true.
         if (resname .eq. ' DG')  deoxy(i) = .true.
         if (resname .eq. ' DC')  deoxy(i) = .true.
         if (resname .eq. ' DT')  deoxy(i) = .true.
      end do
c
c     set the backbone and glycosidic torsions and sugar pucker
c
      do i = 1, nseq
         done = .false.
         do j = 1, 6
            if (bkbone(j,i) .ne. 0.0d0)  done = .true.
         end do
         if (glyco(i) .ne. 0.0d0)  done = .true.
         if (pucker(i) .ne. 0)  done = .true.
         if (.not. done) then
            if (hlxform .eq. 'A') then
               bkbone(1,i) = -52.0d0
               bkbone(2,i) = 175.0d0
               bkbone(3,i) = 42.0d0
               bkbone(4,i) = 79.0d0
               bkbone(5,i) = -148.0d0
               bkbone(6,i) = -75.0d0
               glyco(i) = -157.0d0
               pucker(i) = 3
            else if (hlxform .eq. 'B') then
               bkbone(1,i) = -30.0d0
               bkbone(2,i) = 136.0d0
               bkbone(3,i) = 31.0d0
               bkbone(4,i) = 143.0d0
               bkbone(5,i) = -141.0d0
               bkbone(6,i) = -161.0d0
               glyco(i) = -98.0d0
               pucker(i) = 2
            else if (hlxform .eq. 'Z') then
               if (purine(i)) then
                  bkbone(1,i) = 47.0d0
                  bkbone(2,i) = 179.0d0
                  bkbone(3,i) = -169.0d0
                  bkbone(4,i) = 99.0d0
                  bkbone(5,i) = -104.0d0
                  bkbone(6,i) = -69.0d0
                  glyco(i) = 68.0d0
                  pucker(i) = 3
               else
                  bkbone(1,i) = -137.0d0
                  bkbone(2,i) = -139.0d0
                  bkbone(3,i) = 56.0d0
                  bkbone(4,i) = 138.0d0
                  bkbone(5,i) = -95.0d0
                  bkbone(6,i) = 80.0d0
                  glyco(i) = -159.0d0
                  pucker(i) = 1
               end if
            end if
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (purine)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine nucchain  --  build polynucleotide backbone  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "nucchain" builds up the internal coordinates for a nucleic
c     acid sequence from the sugar type, backbone and glycosidic
c     torsional values
c
c
      subroutine nucchain
      use sizes
      use atoms
      use nucleo
      use resdue
      use sequen
      implicit none
      integer i,k,m
      integer poi,o2i,c1i
      integer c2i,c3i,c4i
      integer c5i,o3i,o4i,o5i
      integer phtyp,ophtyp
      integer ostyp,ottyp
      logical single,cap3,cap5
      character*3 resname
c
c
c     initialize the atom counter to the first atom
c
      n = 1
c
c     check for single residue and 3'- or 5'-phosphate caps
c
      do m = 1, nchain
         single = .false.
         cap5 = .false.
         cap3 = .false.
         if (ichain(1,m) .eq. ichain(2,m))  single = .true.
         i = ichain(1,m)
         k = seqtyp(i)
         resname = nuclz(k)
         if (resname.eq.' MP' .or. resname.eq.' DP'
     &          .or. resname.eq.' TP')  cap5 = .true.
         i = ichain(2,m)
         k = seqtyp(i)
         resname = nuclz(k)
         if (resname.eq.' MP' .or. resname.eq.' DP'
     &          .or. resname.eq.' TP')  cap3 = .true.
c
c     build the first residue or a phosphate capping group
c
         i = ichain(1,m)
         k = seqtyp(i)
         resname = nuclz(k)
         if (resname .eq. ' MP') then
            if (deoxy(i+1)) then
               ostyp = 1246
               phtyp = 1247
               ophtyp = 1248
            else
               ostyp = 1234
               phtyp = 1235
               ophtyp = 1236
            end if
            if (m .eq. 1) then
               o3i = n
               call zatom (ophtyp,0.0d0,0.0d0,0.0d0,0,0,0,0)
               poi = n
               call zatom (phtyp,1.52d0,0.0d0,0.0d0,o3i,0,0,0)
               call zatom (ophtyp,1.52d0,113.0d0,0.0d0,poi,o3i,0,0)
            else
               o3i = n
               call zatom (ophtyp,30.0d0,150.0d0,180.0d0,n-1,n-2,n-3,0)
               call zatom (-2,0.0d0,0.0d0,0.0d0,n-2,n-1,0,0)
               poi = n
               call zatom (phtyp,1.52d0,150.0d0,180.0d0,o3i,n-2,n-3,0)
               call zatom (ophtyp,1.52d0,113.0d0,180.0d0,poi,o3i,n-3,0)
            end if
            call zatom (ophtyp,1.52d0,113.0d0,113.0d0,poi,o3i,n-1,1)
            o5i = n
            call zatom (ostyp,1.63d0,106.0d0,106.0d0,poi,o3i,n-2,-1)
         else if (resname .eq. ' DP') then
            continue
         else if (resname .eq. ' TP') then
            continue
         else
            if (deoxy(i)) then
               ottyp = 1244
            else
               ottyp = 1232
            end if
            if (m .eq. 1) then
               o5i = n
               call zatom (ottyp,0.0d0,0.0d0,0.0d0,0,0,0,0)
               c5i = n
               call zatom (c5typ(k),1.44d0,0.0d0,0.0d0,o5i,0,0,0)
               c4i = n
               call zatom (c4typ(k),1.52d0,110.1d0,0.0d0,c5i,o5i,0,0)
            else
               o5i = n
               call zatom (ottyp,0.96d0,150.0d0,180.0d0,n-1,n-2,n-3,0)
               call zatom (-2,0.0d0,0.0d0,0.0d0,n-2,n-1,0,0)
               c5i = n
               call zatom (c5typ(k),1.44d0,119.0d0,180.0d0,
     &                        o5i,n-2,n-3,0)
               c4i = n
               call zatom (c4typ(k),1.52d0,110.1d0,180.0d0,
     &                        c5i,o5i,n-3,0)
            end if
            o4i = n
            call zatom (o4typ(k),1.46d0,108.9d0,bkbone(3,i)-120.0d0,
     &                     c4i,c5i,o5i,0)
            c1i = n
            if (pucker(i) .eq. 3) then
               call zatom (c1typ(k),1.42d0,109.8d0,145.0d0,
     &                        o4i,c4i,c5i,0)
            else if (pucker(i) .eq. 2) then
               call zatom (c1typ(k),1.42d0,109.8d0,107.0d0,
     &                        o4i,c4i,c5i,0)
            else if (pucker(i) .eq. 1) then
               call zatom (c1typ(k),1.42d0,109.8d0,140.0d0,
     &                        o4i,c4i,c5i,0)
            end if
            c3i = n
            call zatom (c3typ(k),1.53d0,115.9d0,bkbone(3,i),
     &                     c4i,c5i,o5i,0)
            c2i = n
            call zatom (c2typ(k),1.53d0,102.4d0,bkbone(4,i)+120.0d0,
     &                     c3i,c4i,c5i,0)
            call zatom (-1,0.0d0,0.0d0,0.0d0,c1i,c2i,0,0)
            o3i = n
            if (deoxy(i)) then
               if (single) then
                  call zatom (1249,1.42d0,112.1d0,bkbone(4,i),
     &                           c3i,c4i,c5i,0)
               else
                  call zatom (o3typ(k),1.42d0,112.1d0,bkbone(4,i),
     &                           c3i,c4i,c5i,0)
               end if
            else
               if (single) then
                  call zatom (1237,1.42d0,112.1d0,bkbone(4,i),
     &                           c3i,c4i,c5i,0)
               else
                  call zatom (o3typ(k),1.42d0,112.1 d0,bkbone(4,i),
     &                           c3i,c4i,c5i,0)
               end if
               o2i = n
               call zatom (o2typ(k),1.43d0,109.5d0,109.5d0,
     &                        c2i,c3i,c1i,1)
            end if
            call zatom (h5ttyp(k),0.96d0,107.0d0,180.0d0,
     &                     o5i,c5i,c4i,0)
            call zatom (h51typ(k),1.09d0,109.5d0,109.5d0,c5i,o5i,c4i,1)
            call zatom (h52typ(k),1.09d0,109.5d0,109.5d0,c5i,o5i,c4i,-1)
            call zatom (h4typ(k),1.09d0,109.5d0,109.5d0,c4i,c5i,c3i,-1)
            if (pucker(i) .eq. 3) then
               call zatom (h1typ(k),1.09d0,109.5d0,120.0d0,
     &                        c1i,o4i,c2i,-1)
            else if (pucker(i) .eq. 2) then
               call zatom (h1typ(k),1.09d0,109.5d0,115.0d0,
     &                        c1i,o4i,c2i,-1)
            else if (pucker(i) .eq. 1) then
               call zatom (h1typ(k),1.09d0,109.5d0,90.0d0,
     &                        c1i,o4i,c2i,-1)
            end if
            call zatom (h3typ(k),1.09d0,109.5d0,109.5d0,c3i,c4i,c2i,-1)
            call zatom (h21typ(k),1.09d0,109.5d0,109.5d0,c2i,c3i,c1i,-1)
            if (deoxy(i)) then
               call zatom (h22typ(k),1.09d0,109.5d0,109.5d0,
     &                        c2i,c3i,c1i,1)
            else
               call zatom (h22typ(k),0.96d0,107.0d0,180.0d0,
     &                        o2i,c2i,c3i,0)
            end if
            if (single) then
               call zatom (h3ttyp(k),0.96d0,115.0d0,180.0d0,
     &                        o3i,c3i,c4i,0)
            end if
            call nucbase (resname,i,c1i,o4i,c2i)
         end if
c
c     build atoms for residues in the middle of the chain
c
         do i = ichain(1,m)+1, ichain(2,m)-1
            k = seqtyp(i)
            resname = nuclz(k)
            if (cap5) then
               cap5 = .false.
            else
               poi = n
               call zatom (ptyp(k),1.60d0,119.0d0,bkbone(5,i-1),
     &                        o3i,c3i,c4i,0)
               call zatom (optyp(k),1.48d0,109.0d0,
     &                        bkbone(6,i-1)+120.0d0,poi,o3i,c3i,0)
               call zatom (optyp(k),1.48d0,109.0d0,
     &                        bkbone(6,i-1)-120.0d0,poi,o3i,c3i,0)
               o5i = n
               call zatom (o5typ(k),1.60d0,101.8d0,bkbone(6,i-1),
     &                        poi,o3i,c3i,0)
            end if
            c5i = n
            call zatom (c5typ(k),1.44d0,119.0d0,bkbone(1,i),
     &                     o5i,poi,o3i,0)
            c4i = n
            call zatom (c4typ(k),1.52d0,110.1d0,bkbone(2,i),
     &                     c5i,o5i,poi,0)
            o4i = n
            call zatom (o4typ(k),1.46d0,108.9d0,bkbone(3,i)-120.0d0,
     &                     c4i,c5i,o5i,0)
            c1i = n
            if (pucker(i) .eq. 3) then
               call zatom (c1typ(k),1.42d0,109.8d0,145.0d0,
     &                        o4i,c4i,c5i,0)
            else if (pucker(i) .eq. 2) then
               call zatom (c1typ(k),1.42d0,109.8d0,107.0d0,
     &                        o4i,c4i,c5i,0)
            else if (pucker(i) .eq. 1) then
               call zatom (c1typ(k),1.42d0,109.8d0,140.0d0,
     &                        o4i,c4i,c5i,0)
            end if
            c3i = n
            call zatom (c3typ(k),1.53d0,115.9d0,bkbone(3,i),
     &                     c4i,c5i,o5i,0)
            c2i = n
            call zatom (c2typ(k),1.53d0,102.4d0,bkbone(4,i)+120.0d0,
     &                     c3i,c4i,c5i,0)
            call zatom (-1,0.0d0,0.0d0,0.0d0,c1i,c2i,0,0)
            o3i = n
            if (deoxy(i)) then
               if (cap3) then
                  call zatom (1251,1.42d0,112.1d0,bkbone(4,i),
     &                           c3i,c4i,c5i,0)
               else
                  call zatom (o3typ(k),1.42d0,112.1d0,bkbone(4,i),
     &                           c3i,c4i,c5i,0)
               end if
            else
               if (cap3) then
                  call zatom (1239,1.42d0,112.1d0,bkbone(4,i),
     &                           c3i,c4i,c5i,0)
               else
                  call zatom (o3typ(k),1.42d0,112.1d0,bkbone(4,i),
     &                           c3i,c4i,c5i,0)
               end if
               o2i = n
               call zatom (o2typ(k),1.43d0,109.5d0,109.5d0,
     &                        c2i,c3i,c1i,1)
            end if
            call zatom (h51typ(k),1.09d0,109.5d0,109.5d0,c5i,o5i,c4i,1)
            call zatom (h52typ(k),1.09d0,109.5d0,109.5d0,c5i,o5i,c4i,-1)
            call zatom (h4typ(k),1.09d0,109.5d0,109.5d0,c4i,c5i,c3i,-1)
            if (pucker(i) .eq. 3) then
               call zatom (h1typ(k),1.09d0,109.5d0,120.0d0,
     &                        c1i,o4i,c2i,-1)
            else if (pucker(i) .eq. 2) then
               call zatom (h1typ(k),1.09d0,109.5d0,115.0d0,
     &                        c1i,o4i,c2i,-1)
            else if (pucker(i) .eq. 1) then
               call zatom (h1typ(k),1.09d0,109.5d0,90.0d0,
     &                        c1i,o4i,c2i,-1)
            end if
            call zatom (h3typ(k),1.09d0,109.5d0,109.5d0,c3i,c4i,c2i,-1)
            call zatom (h21typ(k),1.09d0,109.5d0,109.5d0,c2i,c3i,c1i,-1)
            if (deoxy(i)) then
               call zatom (h22typ(k),1.09d0,109.5d0,109.5d0,
     &                        c2i,c3i,c1i,1)
            else
               call zatom (h22typ(k),0.96d0,107.0d0,180.0d0,
     &                        o2i,c2i,c3i,0)
            end if
            call nucbase (resname,i,c1i,o4i,c2i)
         end do
c
c     build the last residue or a phosphate capping group
c
         i = ichain(2,m)
         k = seqtyp(i)
         resname = nuclz(k)
         if (single) then
            continue
         else if (resname .eq. ' MP') then
            poi = n
            if (deoxy(i-1)) then
               call zatom (1252,1.63d0,119.0d0,bkbone(5,i-1),
     &                        o3i,c3i,c4i,0)
               call zatom (1253,1.52d0,106.0d0,60.0d0,poi,o3i,c3i,0)
               call zatom (1253,1.52d0,106.0d0,-60.0d0,poi,o3i,c3i,0)
               call zatom (1253,1.52d0,106.0d0,180.0d0,poi,o3i,c3i,0)
            else
               call zatom (1240,1.63d0,119.0d0,bkbone(5,i-1),
     &                        o3i,c3i,c4i,0)
               call zatom (1241,1.52d0,106.0d0,60.0d0,poi,o3i,c3i,0)
               call zatom (1241,1.52d0,106.0d0,-60.0d0,poi,o3i,c3i,0)
               call zatom (1241,1.52d0,106.0d0,180.0d0,poi,o3i,c3i,0)
            end if
         else if (resname .eq. ' DP') then
            continue
         else if (resname .eq. ' TP') then
            continue
         else
            if (cap5) then
               cap5 = .false.
            else
               poi = n
               call zatom (ptyp(k),1.60d0,119.0d0,bkbone(5,i-1),
     &                        o3i,c3i,c4i,0)
               call zatom (optyp(k),1.48d0,109.0d0,
     &                        bkbone(6,i-1)+120.0d0,poi,o3i,c3i,0)
               call zatom (optyp(k),1.48d0,109.0d0,
     &                        bkbone(6,i-1)-120.0d0,poi,o3i,c3i,0)
               o5i = n
               call zatom (o5typ(k),1.60d0,101.8d0,bkbone(6,i-1),
     &                        poi,o3i,c3i,0)
            end if
            c5i = n
            call zatom (c5typ(k),1.44d0,119.0d0,bkbone(1,i),
     &                     o5i,poi,o3i,0)
            c4i = n
            call zatom (c4typ(k),1.52d0,110.1d0,bkbone(2,i),
     &                     c5i,o5i,poi,0)
            o4i = n
            call zatom (o4typ(k),1.46d0,108.9d0,bkbone(3,i)-120.0d0,
     &                     c4i,c5i,o5i,0)
            c1i = n
            if (pucker(i) .eq. 3) then
               call zatom (c1typ(k),1.42d0,109.8d0,145.0d0,
     &                        o4i,c4i,c5i,0)
            else if (pucker(i) .eq. 2) then
               call zatom (c1typ(k),1.42d0,109.8d0,107.0d0,
     &                        o4i,c4i,c5i,0)
            else if (pucker(i) .eq. 1) then
               call zatom (c1typ(k),1.42d0,109.8d0,140.0d0,
     &                        o4i,c4i,c5i,0)
            end if
            c3i = n
            call zatom (c3typ(k),1.53d0,115.9d0,bkbone(3,i),
     &                     c4i,c5i,o5i,0)
            c2i = n
            call zatom (c2typ(k),1.53d0,102.4d0,bkbone(4,i)+120.0d0,
     &                     c3i,c4i,c5i,0)
            call zatom (-1,0.0d0,0.0d0,0.0d0,c1i,c2i,0,0)
            o3i = n
            if (deoxy(i)) then
               call zatom (1249,1.42d0,112.1d0,bkbone(4,i),
     &                        c3i,c4i,c5i,0)
            else
               call zatom (1237,1.42d0,112.1d0,bkbone(4,i),
     &                        c3i,c4i,c5i,0)
               o2i = n
               call zatom (o2typ(k),1.43d0,109.5d0,109.5d0,
     &                        c2i,c3i,c1i,1)
            end if
            call zatom (h51typ(k),1.09d0,109.5d0,109.5d0,c5i,o5i,c4i,1)
            call zatom (h52typ(k),1.09d0,109.5d0,109.5d0,c5i,o5i,c4i,-1)
            call zatom (h4typ(k),1.09d0,109.5d0,109.5d0,c4i,c5i,c3i,-1)
            if (pucker(i) .eq. 3) then
               call zatom (h1typ(k),1.09d0,109.5d0,120.0d0,
     &                        c1i,o4i,c2i,-1)
            else if (pucker(i) .eq. 2) then
               call zatom (h1typ(k),1.09d0,109.5d0,115.0d0,
     &                        c1i,o4i,c2i,-1)
            else if (pucker(i) .eq. 1) then
               call zatom (h1typ(k),1.09d0,109.5d0,90.0d0,
     &                        c1i,o4i,c2i,-1)
            end if
            call zatom (h3typ(k),1.09d0,109.5d0,109.5d0,c3i,c4i,c2i,-1)
            call zatom (h21typ(k),1.09d0,109.5d0,109.5d0,c2i,c3i,c1i,-1)
            if (deoxy(i)) then
               call zatom (h22typ(k),1.09d0,109.5d0,109.5d0,
     &                        c2i,c3i,c1i,1)
            else
               call zatom (h22typ(k),0.96d0,107.0d0,180.0d0,
     &                        o2i,c2i,c3i,0)
            end if
            call zatom (h3ttyp(k),0.96d0,115.0d0,180.0d0,o3i,c3i,c4i,0)
            call nucbase (resname,i,c1i,o4i,c2i)
         end if
      end do
c
c     finally, set the total number of atoms
c
      n = n - 1
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine nucbase  --  build nucleotide base side chain  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "nucbase" builds the side chain for a single nucleotide base
c     in terms of internal coordinates
c
c     resname   3-letter name of current nucleotide residue
c     i         number of the current nucleotide residue
c     c1i       atom number of carbon C1' in residue i
c     o4i       atom number of oxygen O4' in residue i
c     c2i       atom number of carbon C2' in residue i
c
c     literature references:
c
c     R. Lavery, K. Zakrzewska, "Base and Base Pair Morphologies,
c     Helical Parameters, and Definitions" in "Oxford Handbook of
c     Nucleic Acid Structure", S. Neidel, Editor, Oxford University
c     Press, 1999, pages 40-42
c
c     W. Saenger, "Principles of Nucleic Acid Structure", Springer-
c     Verlag, 1984, page 52
c
c
      subroutine nucbase (resname,i,c1i,o4i,c2i)
      use sizes
      use atoms
      use nucleo
      implicit none
      integer i,c1i,o4i,c2i
      character*3 resname
c
c
c     adenine in adenosine residue  (A)
c
      if (resname .eq. '  A') then
         call zatom (1017,1.48d0,108.1d0,113.7d0,c1i,o4i,c2i,1)
         call zatom (1021,1.37d0,128.4d0,glyco(i)+180.0d0,
     &                  n-1,c1i,o4i,0)
         call zatom (1020,1.30d0,113.8d0,180.0d0,n-1,n-2,c1i,0)
         call zatom (1019,1.39d0,104.0d0,0.0d0,n-1,n-2,n-3,0)
         call zatom (1025,1.40d0,132.4d0,180.0d0,n-1,n-2,n-3,0)
         call zatom (1027,1.34d0,123.5d0,0.0d0,n-1,n-2,n-3,0)
         call zatom (1024,1.35d0,117.4d0,180.0d0,n-2,n-3,n-4,0)
         call zatom (1023,1.33d0,118.8d0,0.0d0,n-1,n-3,n-4,0)
         call zatom (1022,1.32d0,129.2d0,0.0d0,n-1,n-2,n-4,0)
         call zatom (1018,1.35d0,110.9d0,0.0d0,n-1,n-2,n-3,0)
         call zatom (-1,0.0d0,0.0d0,0.0d0,n-1,n-7,0,0)
         call zatom (-1,0.0d0,0.0d0,0.0d0,n-1,n-10,0,0)
         call zatom (1030,1.08d0,123.1d0,180.0d0,n-9,n-8,n-7,0)
         call zatom (1028,1.00d0,120.0d0,180.0d0,n-6,n-7,n-8,0)
         call zatom (1029,1.00d0,120.0d0,0.0d0,n-7,n-8,n-9,0)
         call zatom (1026,1.08d0,115.4d0,180.0d0,n-6,n-5,n-4,0)
c
c     guanine in guanosine residue  (G)
c
      else if (resname .eq. '  G') then
         call zatom (1047,1.48d0,108.1d0,113.7d0,c1i,o4i,c2i,1)
         call zatom (1051,1.38d0,128.4d0,glyco(i)+180.0d0,
     &                  n-1,c1i,o4i,0)
         call zatom (1050,1.31d0,114.0d0,180.0d0,n-1,n-2,c1i,0)
         call zatom (1049,1.39d0,103.8d0,0.0d0,n-1,n-2,n-3,0)
         call zatom (1055,1.40d0,130.1d0,180.0d0,n-1,n-2,n-3,0)
         call zatom (1060,1.23d0,128.8d0,0.0d0,n-1,n-2,n-3,0)
         call zatom (1054,1.40d0,111.4d0,180.0d0,n-2,n-3,n-4,0)
         call zatom (1053,1.38d0,125.2d0,0.0d0,n-1,n-3,n-4,0)
         call zatom (1057,1.34d0,116.1d0,180.0d0,n-1,n-2,n-4,0)
         call zatom (1052,1.33d0,123.3d0,0.0d0,n-2,n-3,n-4,0)
         call zatom (1048,1.36d0,112.3d0,0.0d0,n-1,n-3,n-4,0)
         call zatom (-1,0.0d0,0.0d0,0.0d0,n-1,n-8,0,0)
         call zatom (-1,0.0d0,0.0d0,0.0d0,n-1,n-11,0,0)
         call zatom (1061,1.08d0,123.0d0,180.0d0,n-10,n-9,n-8,0)
         call zatom (1056,1.00d0,117.4d0,180.0d0,n-6,n-8,n-9,0)
         call zatom (1058,1.00d0,120.0d0,0.0d0,n-5,n-6,n-7,0)
         call zatom (1059,1.00d0,120.0d0,180.0d0,n-6,n-7,n-8,0)
c
c     cytosine in cytidine residue  (C)
c
      else if (resname .eq. '  C') then
         call zatom (1078,1.48d0,108.1d0,113.7d0,c1i,o4i,c2i,1)
         call zatom (1079,1.37d0,117.8d0,glyco(i),n-1,c1i,o4i,0)
         call zatom (1084,1.24d0,118.9d0,0.0d0,n-1,n-2,c1i,0)
         call zatom (1080,1.38d0,118.7d0,180.0d0,n-2,n-3,c1i,0)
         call zatom (1081,1.34d0,120.6d0,0.0d0,n-1,n-3,n-4,0)
         call zatom (1085,1.32d0,118.3d0,180.0d0,n-1,n-2,n-4,0)
         call zatom (1082,1.43d0,121.6d0,0.0d0,n-2,n-3,n-5,0)
         call zatom (1083,1.36d0,116.9d0,0.0d0,n-1,n-3,n-4,0)
         call zatom (-1,0.0d0,0.0d0,0.0d0,n-1,n-8,0,0)
         call zatom (1086,1.00d0,120.0d0,0.0d0,n-3,n-4,n-5,0)
         call zatom (1087,1.00d0,120.0d0,180.0d0,n-4,n-5,n-6,0)
         call zatom (1088,1.08d0,121.6d0,180.0d0,n-4,n-6,n-7,0)
         call zatom (1089,1.08d0,119.5d0,180.0d0,n-4,n-5,n-7,0)
c
c     uracil in uridine residue  (U)
c
      else if (resname .eq. '  U') then
         call zatom (1106,1.48d0,108.1d0,113.7d0,c1i,o4i,c2i,1)
         call zatom (1107,1.38d0,117.1d0,glyco(i),n-1,c1i,o4i,0)
         call zatom (1112,1.22d0,123.2d0,0.0d0,n-1,n-2,c1i,0)
         call zatom (1108,1.37d0,114.8d0,180.0d0,n-2,n-3,c1i,0)
         call zatom (1109,1.38d0,127.0d0,0.0d0,n-1,n-3,n-4,0)
         call zatom (1114,1.23d0,119.8d0,180.0d0,n-1,n-2,n-4,0)
         call zatom (1110,1.44d0,114.7d0,0.0d0,n-2,n-3,n-5,0)
         call zatom (1111,1.34d0,119.2d0,0.0d0,n-1,n-3,n-4,0)
         call zatom (-1,0.0d0,0.0d0,0.0d0,n-1,n-8,0,0)
         call zatom (1113,1.00d0,116.5d0,180.0d0,n-5,n-7,n-8,0)
         call zatom (1115,1.08d0,120.4d0,180.0d0,n-3,n-5,n-6,0)
         call zatom (1116,1.08d0,118.6d0,180.0d0,n-3,n-4,n-6,0)
c
c     adenine in deoxyadenosine residue  (DA)
c
      else if (resname .eq. ' DA') then
         call zatom (1132,1.48d0,108.1d0,113.7d0,c1i,o4i,c2i,1)
         call zatom (1136,1.37d0,128.4d0,glyco(i)+180.0d0,
     &                  n-1,c1i,o4i,0)
         call zatom (1135,1.30d0,113.8d0,180.0d0,n-1,n-2,c1i,0)
         call zatom (1134,1.39d0,104.0d0,0.0d0,n-1,n-2,n-3,0)
         call zatom (1140,1.40d0,132.4d0,180.0d0,n-1,n-2,n-3,0)
         call zatom (1142,1.34d0,123.5d0,0.0d0,n-1,n-2,n-3,0)
         call zatom (1139,1.35d0,117.4d0,180.0d0,n-2,n-3,n-4,0)
         call zatom (1138,1.33d0,118.8d0,0.0d0,n-1,n-3,n-4,0)
         call zatom (1137,1.32d0,129.2d0,0.0d0,n-1,n-2,n-4,0)
         call zatom (1133,1.35d0,110.9d0,0.0d0,n-1,n-2,n-3,0)
         call zatom (-1,0.0d0,0.0d0,0.0d0,n-1,n-7,0,0)
         call zatom (-1,0.0d0,0.0d0,0.0d0,n-1,n-10,0,0)
         call zatom (1145,1.08d0,123.1d0,180.0d0,n-9,n-8,n-7,0)
         call zatom (1143,1.00d0,120.0d0,180.0d0,n-6,n-7,n-8,0)
         call zatom (1144,1.00d0,120.0d0,0.0d0,n-7,n-8,n-9,0)
         call zatom (1141,1.08d0,115.4d0,180.0d0,n-6,n-5,n-4,0)
c
c     guanine in deoxyguanosine residue  (DG)
c
      else if (resname .eq. ' DG') then
         call zatom (1161,1.48d0,108.1d0,113.7d0,c1i,o4i,c2i,1)
         call zatom (1165,1.38d0,128.4d0,glyco(i)+180.0d0,
     &                  n-1,c1i,o4i,0)
         call zatom (1164,1.31d0,114.0d0,180.0d0,n-1,n-2,c1i,0)
         call zatom (1163,1.39d0,103.8d0,0.0d0,n-1,n-2,n-3,0)
         call zatom (1169,1.40d0,130.1d0,180.0d0,n-1,n-2,n-3,0)
         call zatom (1174,1.23d0,128.8d0,0.0d0,n-1,n-2,n-3,0)
         call zatom (1168,1.40d0,111.4d0,180.0d0,n-2,n-3,n-4,0)
         call zatom (1167,1.38d0,125.2d0,0.0d0,n-1,n-3,n-4,0)
         call zatom (1171,1.34d0,116.1d0,180.0d0,n-1,n-2,n-4,0)
         call zatom (1166,1.33d0,123.3d0,0.0d0,n-2,n-3,n-4,0)
         call zatom (1162,1.36d0,112.3d0,0.0d0,n-1,n-3,n-4,0)
         call zatom (-1,0.0d0,0.0d0,0.0d0,n-1,n-8,0,0)
         call zatom (-1,0.0d0,0.0d0,0.0d0,n-1,n-11,0,0)
         call zatom (1175,1.08d0,123.0d0,180.0d0,n-10,n-9,n-8,0)
         call zatom (1170,1.00d0,117.4d0,180.0d0,n-6,n-8,n-9,0)
         call zatom (1172,1.00d0,120.0d0,0.0d0,n-5,n-6,n-7,0)
         call zatom (1173,1.00d0,120.0d0,180.0d0,n-6,n-7,n-8,0)
c
c     cytosine in deoxycytidine residue  (DC)
c
      else if (resname .eq. ' DC') then
         call zatom (1191,1.48d0,108.1d0,113.7d0,c1i,o4i,c2i,1)
         call zatom (1192,1.37d0,117.8d0,glyco(i),n-1,c1i,o4i,0)
         call zatom (1197,1.24d0,118.9d0,0.0d0,n-1,n-2,c1i,0)
         call zatom (1193,1.38d0,118.7d0,180.0d0,n-2,n-3,c1i,0)
         call zatom (1194,1.34d0,120.6d0,0.0d0,n-1,n-3,n-4,0)
         call zatom (1198,1.32d0,118.3d0,180.0d0,n-1,n-2,n-4,0)
         call zatom (1195,1.43d0,121.6d0,0.0d0,n-2,n-3,n-5,0)
         call zatom (1196,1.36d0,116.9d0,0.0d0,n-1,n-3,n-4,0)
         call zatom (-1,0.0d0,0.0d0,0.0d0,n-1,n-8,0,0)
         call zatom (1199,1.00d0,120.0d0,0.0d0,n-3,n-4,n-5,0)
         call zatom (1200,1.00d0,120.0d0,180.0d0,n-4,n-5,n-6,0)
         call zatom (1201,1.08d0,121.6d0,180.0d0,n-4,n-6,n-7,0)
         call zatom (1202,1.08d0,119.5d0,180.0d0,n-4,n-5,n-7,0)
c
c     thymine in deoxythymidine residue  (DT)
c
      else if (resname .eq. ' DT') then
         call zatom (1218,1.48d0,108.1d0,113.7d0,c1i,o4i,c2i,1)
         call zatom (1219,1.37d0,117.1d0,glyco(i),n-1,c1i,o4i,0)
         call zatom (1224,1.22d0,122.9d0,0.0d0,n-1,n-2,c1i,0)
         call zatom (1220,1.38d0,115.4d0,180.0d0,n-2,n-3,c1i,0)
         call zatom (1221,1.38d0,126.4d0,0.0d0,n-1,n-3,n-4,0)
         call zatom (1226,1.23d0,120.5d0,180.0d0,n-1,n-2,n-4,0)
         call zatom (1222,1.44d0,114.1d0,0.0d0,n-2,n-3,n-5,0)
         call zatom (1227,1.50d0,117.5d0,180.0d0,n-1,n-3,n-4,0)
         call zatom (1223,1.34d0,120.8d0,0.0d0,n-2,n-4,n-5,0)
         call zatom (-1,0.0d0,0.0d0,0.0d0,n-1,n-9,0,0)
         call zatom (1225,1.00d0,116.8d0,180.0d0,n-6,n-8,n-9,0)
         call zatom (1228,1.09d0,109.5d0,0.0d0,n-3,n-4,n-6,0)
         call zatom (1228,1.09d0,109.5d0,109.5d0,n-4,n-5,n-1,1)
         call zatom (1228,1.09d0,109.5d0,109.5d0,n-5,n-6,n-2,-1)
         call zatom (1229,1.08d0,119.4d0,180.0d0,n-5,n-7,n-9,0)
      end if
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine watson  --  align strands of a double helix  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "watson" uses a rigid body optimization to approximately
c     align the paired strands of a nucleic acid double helix
c
c
      subroutine watson
      use sizes
      use atoms
      use couple
      use group
      use inform
      use katoms
      use molcul
      use nucleo
      use output
      use potent
      use resdue
      use restrn
      use rigid
      use sequen
      use usage
      implicit none
      integer i,j,nvar
      integer ia,ib,ic,id
      integer start,stop
      integer kseq,offset
      integer nbase,nphos
      integer, allocatable :: iphos(:)
      integer, allocatable :: root(:)
      integer, allocatable :: list(:,:)
      real*8 minimum,grdmin
      real*8 watson1,sum,dist
      real*8, allocatable :: xx(:)
      character*3 resname
      external watson1,optsave
c
c
c     perform dynamic allocation of some global arrays
c
      if (.not. allocated(iuse))  allocate (iuse(n))
      if (.not. allocated(use))  allocate (use(0:n))
c
c     set all atoms to be active during energy evaluations
c
      nuse = n
      do i = 1, n
         use(i) = .true.
      end do
c
c     only geometric restraints will by used in optimization
c
      call potoff
      use_geom = .true.
c
c     set the default values for the restraint variables
c
      npfix = 0
      ndfix = 0
      ntfix = 0
      ngfix = 0
      nchir = 0
      use_basin = .false.
      use_wall = .false.
c
c     perform dynamic allocation of some local arrays
c
      allocate (iphos(nseq+10))
      allocate (root(nseq))
      allocate (list(2,nseq))
c
c     find root atom and hydrogen bond partners for each base
c
      kseq = 0
      nbase = 0
      do i = 1, n
         if (atmnum(type(i)).eq.6 .and. n12(i).eq.4) then
            ia = atmnum(type(i12(1,i)))
            ib = atmnum(type(i12(2,i)))
            ic = atmnum(type(i12(3,i)))
            id = atmnum(type(i12(4,i)))
            sum = ia + ib + ic + id
            if (sum .eq. 22) then
               nbase = nbase + 1
               j = i12(4,i)
               root(nbase) = j
               kseq = kseq + 1
               resname = nuclz(seqtyp(kseq))
               do while (resname.eq.' MP' .or. resname.eq.' DP'
     &                         .or. resname.eq.' TP')
                  kseq = kseq + 1
                  resname = nuclz(seqtyp(kseq))
               end do
               if (resname.eq.'  A' .or. resname.eq.' DA') then
                  list(1,nbase) = j + 6
                  list(2,nbase) = j + 11
               else if (resname.eq.'  G' .or. resname.eq.' DG') then
                  list(1,nbase) = j + 12
                  list(2,nbase) = j + 5
               else if (resname.eq.'  C' .or. resname.eq.' DC') then
                  list(1,nbase) = j + 3
                  list(2,nbase) = j + 8
               else if (resname .eq. '  U') then
                  list(1,nbase) = j + 8
                  list(2,nbase) = j + 5
               else if (resname .eq. ' DT') then
                  list(1,nbase) = j + 9
                  list(2,nbase) = j + 5
               end if
            end if
         end if
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (.not. allocated(idfix))  allocate (idfix(2,maxfix))
      if (.not. allocated(dfix))  allocate (dfix(3,maxfix))
      if (.not. allocated(itfix))  allocate (itfix(4,maxfix))
      if (.not. allocated(tfix))  allocate (tfix(3,maxfix))
c
c     distance restraints for the base pair hydrogen bonds
c
      do i = 1, nbase/2
         j = nbase + 1 - i
         ndfix = ndfix + 1
         idfix(1,ndfix) = list(1,i)
         idfix(2,ndfix) = list(1,j)
         dfix(1,ndfix) = 50.0d0
         dfix(2,ndfix) = 1.85d0
         dfix(3,ndfix) = 1.95d0
         ndfix = ndfix + 1
         idfix(1,ndfix) = list(2,i)
         idfix(2,ndfix) = list(2,j)
         dfix(1,ndfix) = 50.0d0
         dfix(2,ndfix) = 1.85d0
         dfix(3,ndfix) = 1.95d0
      end do
c
c     torsional restraints to enforce base pair planarity
c
      do i = 1, nbase/2
         j = nbase + 1 - i
         ntfix = ntfix + 1
         itfix(1,ntfix) = root(i)
         itfix(2,ntfix) = list(1,i)
         itfix(3,ntfix) = list(2,i)
         itfix(4,ntfix) = list(1,j)
         tfix(1,ntfix) = 2.5d0
         tfix(2,ntfix) = 180.0d0
         tfix(3,ntfix) = 180.0d0
         ntfix = ntfix + 1
         itfix(1,ntfix) = root(i)
         itfix(2,ntfix) = list(2,i)
         itfix(3,ntfix) = list(1,i)
         itfix(4,ntfix) = list(2,j)
         tfix(1,ntfix) = 2.5d0
         tfix(2,ntfix) = 180.0d0
         tfix(3,ntfix) = 180.0d0
         ntfix = ntfix + 1
         itfix(1,ntfix) = root(j)
         itfix(2,ntfix) = list(1,j)
         itfix(3,ntfix) = list(2,j)
         itfix(4,ntfix) = list(1,i)
         tfix(1,ntfix) = 2.5d0
         tfix(2,ntfix) = 180.0d0
         tfix(3,ntfix) = 180.0d0
         ntfix = ntfix + 1
         itfix(1,ntfix) = root(j)
         itfix(2,ntfix) = list(2,j)
         itfix(3,ntfix) = list(1,j)
         itfix(4,ntfix) = list(2,i)
         tfix(1,ntfix) = 2.5d0
         tfix(2,ntfix) = 180.0d0
         tfix(3,ntfix) = 180.0d0
      end do
c
c     distance restraints between interstrand phosphates
c
      nphos = 0
      do i = 1, n
         if (atmnum(type(i)) .eq. 15) then
            nphos = nphos + 1
            iphos(nphos) = i
         end if
      end do
      start = 1
      stop = nphos / 2
      resname = nuclz(seqtyp(1))
      if (resname .eq. ' MP')  start = start + 1
      if (resname .eq. ' DP')  start = start + 2
      if (resname .eq. ' TP')  start = start + 3
      resname = nuclz(seqtyp(nseq))
      if (resname .eq. ' MP')  stop = stop - 1
      if (resname .eq. ' DP')  stop = stop - 2
      if (resname .eq. ' TP')  stop = stop - 3
      offset = stop + nphos/2 + 1
      if (hlxform .eq. 'A')  dist = 17.78d0
      if (hlxform .eq. 'B')  dist = 17.46d0
      if (hlxform .eq. 'Z')  dist = 13.2d0
      do i = start, stop
         ndfix = ndfix + 1
         idfix(1,ndfix) = iphos(i)
         idfix(2,ndfix) = iphos(offset-i)
         dfix(1,ndfix) = 100.0d0
         dfix(2,ndfix) = dist
         dfix(3,ndfix) = dist
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (iphos)
      deallocate (root)
      deallocate (list)
c
c     enable use of groups based on number of molecules
c
      use_group = .true.
      ngrp = nmol
c
c     perform dynamic allocation of some global arrays
c
      if (.not. allocated(kgrp))  allocate (kgrp(n))
      if (.not. allocated(grplist))  allocate (grplist(n))
      if (.not. allocated(igrp))  allocate (igrp(2,0:ngrp))
      if (.not. allocated(grpmass))  allocate (grpmass(0:ngrp))
      if (.not. allocated(wgrp))  allocate (wgrp(0:ngrp,0:ngrp))
c
c     assign each strand to a separate molecule-based group
c
      do i = 1, ngrp
         igrp(1,i) = imol(1,i)
         igrp(2,i) = imol(2,i)
         do j = igrp(1,i), igrp(2,i)
            kgrp(j) = kmol(j)
            grplist(kgrp(j)) = i
         end do
      end do
      do i = 0, ngrp
         do j = 0, ngrp
            wgrp(j,i) = 1.0d0
         end do
         wgrp(i,i) = 0.0d0
      end do
c
c     get rigid body reference coordinates for each strand
c
      call orient
c
c     perform dynamic allocation of some local arrays
c
      allocate (xx(6*ngrp))
c
c     transfer rigid body coordinates to optimization parameters
c
      nvar = 0
      do i = 1, ngrp
         do j = 1, 6
            nvar = nvar + 1
            xx(nvar) = rbc(j,i)
         end do
      end do
c
c     make the call to the optimization routine
c
      iprint = 0
      iwrite = 0
      grdmin = 0.1d0
      coordtype = 'NONE'
      call ocvm (nvar,xx,minimum,grdmin,watson1,optsave)
c
c     transfer optimization parameters to rigid body coordinates
c
      nvar = 0
      do i = 1, ngrp
         do j = 1, 6
            nvar = nvar + 1
            rbc(j,i) = xx(nvar)
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xx)
c
c     convert from rigid body to Cartesian coordinates
c
      call rigidxyz
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  function watson1  --  energy and gradient for watson  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "watson1" is a service routine that computes the energy
c     and gradient for optimally conditioned variable metric
c     optimization of rigid bodies
c
c
      function watson1 (xx,g)
      use sizes
      use group
      use math
      use rigid
      implicit none
      integer i,j,nvar
      real*8 watson1,e
      real*8 xx(*)
      real*8 g(*)
      real*8, allocatable :: derivs(:,:)
c
c
c     translate optimization parameters to rigid body coordinates
c
      nvar = 0
      do i = 1, ngrp
         do j = 1, 6
            nvar = nvar + 1
            rbc(j,i) = xx(nvar)
         end do
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (derivs(6,ngrp))
c
c     compute and store the energy and gradient
c
      call rigidxyz
      call gradrgd (e,derivs)
      watson1 = e
c
c     translate rigid body gradient to optimization gradient
c
      nvar = 0
      do i = 1, ngrp
         do j = 1, 6
            nvar = nvar + 1
            g(nvar) = derivs(j,i)
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (derivs)
      return
      end
