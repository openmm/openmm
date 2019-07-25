c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  program protein  --  build a polypeptide from sequence  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "protein" builds the internal and Cartesian coordinates
c     of a polypeptide from amino acid sequence and torsional
c     angle values for the peptide backbone and side chains
c
c
      program protein
      use sizes
      use atoms
      use files
      use iounit
      use sequen
      use titles
      implicit none
      integer i,izmt
      integer ixyz,iseq
      integer natom,mode
      integer freeunit,trimtext
      logical exist,clash
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
      call getseq
      call prochain
      call connect
      call attach
      call molecule
      call makexyz
c
c     perform a packing calculation for multiple chains
c
      if (nchain .gt. 1) then
         call pauling
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
c     check for atom pairs with identical coordinates
c
      clash = .false.
      call chkxyz (clash)
c
c     write out a amino acid sequence file
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
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine getseq  --  amino acid sequence and angles  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "getseq" asks the user for the amino acid sequence
c     and torsional angle values needed to define a peptide
c
c
      subroutine getseq
      use iounit
      use phipsi
      use resdue
      use sequen
      implicit none
      integer i,j,next
      integer length,trimtext
      logical done
      character*1 chir
      character*1 ucase(26)
      character*3 name
      character*240 record
      character*240 string
      data ucase  / 'A','B','C','D','E','F','G','H','I','J','K','L',
     &              'M','N','O','P','Q','R','S','T','U','V','W','X',
     &              'Y','Z' /
c
c
c     provide a header to explain the method of sequence input
c
      write (iout,10)
   10 format (/,' Enter One Residue Name per Line as the Standard',
     &           ' Three-Letter Code, then',
     &        /,' Phi Psi Omega (3F), Chi Angles (4F), then',
     &           ' Disulfide Partner if CYX (I),',
     &        /,' and D/L Chirality as Desired (A1)',
     &        //,' If Only Residue Names are Entered, the Default',
     &           ' is to Build an Extended',
     &        /,' Conformation Using L-Amino Acids and Zwitterionic',
     &           ' Termini',
     &        //,' Regular Amino Acids:  GLY, ALA, VAL, LEU,',
     &           ' ILE, SER, THR, CYS, CYX, PRO,',
     &        /,' PHE, TYR, TRP, HIS, ASP, ASN, GLU, GLN, MET,',
     &           ' LYS, ARG, ORN, AIB',
     &        //,' Alternative Protonation States:  CYD, TYD,',
     &           ' HID, HIE, ASH, GLH, LYD',
     &        //,' N-Terminal Cap Residues:  H2N=Deprotonated,',
     &           ' FOR=Formyl, ACE=Acetyl,',
     &        /,27x,'PCA=Pyroglutamic Acid',
     &        /,' C-Terminal Cap Residues:  COH=Protonated,',
     &           ' NH2=Amide, NME=N-MethylAmide',
     &        //,' Use Residue Name=MOL to Start a New Chain,',
     &           ' and Use <CR> to End Input')
c
c     initially, assume that only a single strand is present
c
      nchain = 1
      ichain(1,1) = 1
      chnnam(1) = ' '
c
c     get the amino acid sequence data and dihedral angle values
c
      i = 0
      done = .false.
      do while (.not. done)
         i = i + 1
         phi(i) = 0.0d0
         psi(i) = 0.0d0
         omega(i) = 0.0d0
         do j = 1, 4
            chi(j,i) = 0.0d0
         end do
         chiral(i) = 1
         disulf(i) = 0
         chir = ' '
         write (iout,20)  i
   20    format (/,' Enter Residue',i4,' :  ',$)
         read (input,30)  record
   30    format (a240)
         call upcase (record)
         next = 1
         call gettext (record,name,next)
         length = trimtext (name)
         string = record(next:240)
         read (string,*,err=40,end=40)  phi(i),psi(i),omega(i),
     &                                  (chi(j,i),j=1,4),disulf(i)
   40    continue
         call getword (record,chir,next)
c
c     handle special names used for certain amino acids
c
         if (name .eq. 'CYH')  name = 'CYS'
         if (name .eq. 'CSS')  name = 'CYX'
         if (name .eq. 'HIP')  name = 'HIS'
c
c     disulfide bridged residues are cystine instead of cysteine
c
         if (name(1:1).eq.'C' .and. disulf(i).ne.0) then
            length = 3
            name = 'CYX'
         end if
c
c     check the D/L chirality of the current residue
c
         if (chir .eq. 'D')  chiral(i) = -1
c
c     process and store the current amino acid residue type
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
               seq(i) = amino(maxamino)
               seqtyp(i) = 0
               if (length .eq. 1) then
                  do j = 1, maxamino
                     if (name(1:1) .eq. amino1(j)) then
                        seq(i) = amino(j)
                        seqtyp(i) = j
                     end if
                  end do
               else if (length .eq. 3) then
                  do j = 1, maxamino
                     if (name .eq. amino(j)) then
                        seq(i) = amino(j)
                        seqtyp(i) = j
                     end if
                  end do
               end if
               if (seqtyp(i) .eq. 0) then
                  i = i - 1
                  write (iout,50)  name
   50             format (/,' GETSEQ  --  Amino Acid Type ',a3,
     &                       ' is Not Supported')
               end if
            end if
         end if
      end do
c
c     set chain identifiers if multiple chains are present
c
      if (nchain .gt. 1) then
         do i = 1, nchain
            chnnam(i) = ucase(i)
         end do
      end if
c
c     set default values for the phi-psi-omega-chi angles;
c     use extended values if no phi-psi values were given
c
      do i = 1, nseq
         if (phi(i).eq.0.0d0 .and. psi(i).eq.0.0d0) then
            phi(i) = -135.0d0
            if (seq(i) .eq. 'PRO')  phi(i) = -60.0d0
            psi(i) = 135.0d0
         end if
         if (omega(i) .eq. 0.0d0) then
            omega(i) = 180.0d0
         end if
         if (chi(1,i) .eq. 0.0d0) then
            do j = 1, 4
               chi(j,i) = 180.0d0
               if (seq(i) .eq. 'PRO')  chi(j,i) = 0.0d0
               if (seq(i) .eq. 'PCA')  chi(j,i) = 0.0d0
            end do
            if (seq(i) .eq. 'PHE')  chi(2,i) = 90.0d0
            if (seq(i) .eq. 'TYR')  chi(2,i) = 90.0d0
            if (seq(i) .eq. 'TYD')  chi(2,i) = 90.0d0
            if (seq(i) .eq. 'TRP')  chi(2,i) = 90.0d0
            if (seq(i) .eq. 'HIS')  chi(2,i) = 90.0d0
            if (seq(i) .eq. 'HID')  chi(2,i) = 90.0d0
            if (seq(i) .eq. 'HIE')  chi(2,i) = 90.0d0
         end if
c
c     check for the presence of any disulfide bonds
c
         if (disulf(i) .ne. 0) then
            if (seq(i) .ne. 'CYX') then
               write (iout,60)  i
   60          format (' GETSEQ  --  Error in Disulfide Bond',
     &                    ' at Residue',i5)
            end if
            if (i.lt.disulf(i) .and. disulf(disulf(i)).ne.i) then
               write (iout,70)  i,disulf(i)
   70          format (' GETSEQ  --  Error in Disulfide Bond',
     &                    ' at Residue',i5,' or',i5)
            end if
         end if
      end do
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine prochain  --  build polypeptide backbone  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "prochain" builds up the internal coordinates for an amino
c     acid sequence from the phi, psi, omega and chi values
c
c
      subroutine prochain
      use sizes
      use atoms
      use iounit
      use phipsi
      use resdue
      use sequen
      implicit none
      integer i,k,m
      integer  next,nsave
      integer, allocatable :: ni(:)
      integer, allocatable :: cai(:)
      integer, allocatable :: ci(:)
      logical single,cyclic
      character*1 answer
      character*3 resname
      character*240 record
c
c
c     determine whether the peptide chain is cyclic
c
      cyclic = .false.
      write (iout,10)
   10 format (/,' Cyclize the Polypeptide Chain [N] :  ',$)
      read (input,20)  record
   20 format (a240)
      next = 1
      call gettext (record,answer,next)
      call upcase (answer)
      if (answer .eq. 'Y')  cyclic = .true.
c
c     perform dynamic allocation of some local arrays
c
      allocate (ni(nseq))
      allocate (cai(nseq))
      allocate (ci(nseq))
c
c     initialize the atom counter to the first atom
c
      n = 1
c
c     set the first residue number and get the type and name
c
      do m = 1, nchain
         single = .false.
         if (ichain(1,m) .eq. ichain(2,m))  single = .true.
         i = ichain(1,m)
         k = seqtyp(i)
         resname = amino(k)
c
c     build the first residue for a cyclic peptide
c
         if (cyclic) then
            if (m .eq. 1) then
               ni(i) = n
               call zatom (ntyp(k),0.0d0,0.0d0,0.0d0,0,0,0,0)
               cai(i) = n
               call zatom (catyp(k),1.46d0,0.0d0,0.0d0,ni(i),0,0,0)
               ci(i) = n
               call zatom (ctyp(k),1.51d0,110.7d0,0.0d0,
     &                     cai(i),ni(i),0,0)
            else
               ni(i) = n
               call zatom (ntyp(k),30.0d0,150.0d0,180.0d0,n-1,n-2,n-3,0)
               call zatom (-2,0.0d0,0.0d0,0.0d0,n-2,n-1,0,0)
               cai(i) = n
               call zatom (catyp(k),1.46d0,150.0d0,180.0d0,
     &                     ni(i),n-2,n-3,0)
               ci(i) = n
               call zatom (ctyp(k),1.51d0,110.7d0,180.0d0,
     &                     cai(i),ni(i),n-3,0)
            end if
            call zatom (otyp(k),1.22d0,122.5d0,psi(i)-180.0d0,
     &                  ci(i),cai(i),ni(i),0)
            call zatom (hntyp(k),1.02d0,121.0d0,phi(i)-180.0d0,
     &                  ni(i),cai(i),ci(i),0)
            call zatom (hatyp(k),1.11d0,109.5d0,107.9d0,
     &                  cai(i),ni(i),ci(i),-chiral(i))
            call proside (resname,i,cai(i),ni(i),ci(i))
c
c     build the first residue as an N-terminal formyl group
c
         else if (resname .eq. 'FOR') then
            if (m .eq. 1) then
               ci(i) = n
               call zatom (cntyp(k),0.0d0,0.0d0,0.0d0,0,0,0,0)
               ni(i) = n
               call zatom (ontyp(k),1.22d0,0.0d0,0.0d0,n-1,0,0,0)
               cai(i) = n
               call zatom (hantyp(k),1.12d0,120.0d0,0.0d0,n-2,n-1,0,0)
            else
               ci(i) = n
               call zatom (cntyp(k),30.0d0,150.0d0,180.0d0,
     &                     n-1,n-2,n-3,0)
               call zatom (-2,0.0d0,0.0d0,0.0d0,n-2,n-1,0,0)
               ni(i) = n
               call zatom (ontyp(k),1.22d0,150.0d0,180.0d0,
     &                     n-1,n-2,n-3,0)
               cai(i) = n
               call zatom (hantyp(k),1.12d0,120.0d0,0.0d0,n-2,n-1,n-3,0)
            end if
            psi(i) = 180.0d0
c
c     build the first residue as an N-terminal acetyl group
c
         else if (resname .eq. 'ACE') then
            if (m .eq. 1) then
               cai(i) = n
               call zatom (cantyp(k),0.0d0,0.0d0,0.0d0,0,0,0,0)
               ci(i) = n
               call zatom (cntyp(k),1.51d0,0.0d0,0.0d0,n-1,0,0,0)
               call zatom (ontyp(k),1.22d0,122.5d0,0.0d0,n-1,n-2,0,0)
            else
               cai(i) = n
               call zatom (cantyp(k),30.0d0,150.0d0,180.0d0,
     &                     n-1,n-2,n-3,0)
               call zatom (-2,0.0d0,0.0d0,0.0d0,n-2,n-1,0,0)
               ci(i) = n
               call zatom (cntyp(k),1.51d0,150.0d0,180.0d0,
     &                     n-1,n-2,n-3,0)
               call zatom (ontyp(k),1.22d0,122.5d0,0.0d0,
     &                     n-1,n-2,n-3,0)
            end if
            ni(i) = n
            call zatom (hantyp(k),1.11d0,107.9d0,0.0d0,n-3,n-2,n-1,0)
            call zatom (hantyp(k),1.11d0,107.9d0,109.4d0,n-4,n-3,n-1,1)
            call zatom (hantyp(k),1.11d0,107.9d0,109.4d0,n-5,n-4,n-2,-1)
            psi(i) = 180.0d0
c
c     build the first residue as a proline
c
         else if (resname .eq. 'PRO') then
            if (m .eq. 1) then
               ni(i) = n
               call zatom (nntyp(k),0.0d0,0.0d0,0.0d0,0,0,0,0)
               cai(i) = n
               call zatom (cantyp(k),1.50d0,0.0d0,0.0d0,ni(i),0,0,0)
               ci(i) = n
               if (single) then
                  call zatom (cctyp(k),1.51d0,111.6d0,0.0d0,
     &                        cai(i),ni(i),0,0)
               else
                  call zatom (cntyp(k),1.51d0,110.7d0,0.0d0,
     &                        cai(i),ni(i),0,0)
               end if
            else
               ni(i) = n
               call zatom (nntyp(k),30.0d0,150.0d0,180.0d0,
     &                     n-1,n-2,n-3,0)
               call zatom (-2,0.0d0,0.0d0,0.0d0,n-2,n-1,0,0)
               cai(i) = n
               call zatom (cantyp(k),1.50d0,150.0d0,180.0d0,
     &                     ni(i),n-2,n-3,0)
               ci(i) = n
               if (single) then
                  call zatom (cctyp(k),1.51d0,111.6d0,180.0d0,
     &                        cai(i),ni(i),n-3,0)
               else
                  call zatom (cntyp(k),1.51d0,110.7d0,180.0d0,
     &                        cai(i),ni(i),n-3,0)
               end if
            end if
            if (single) then
               call zatom (octyp(k),1.25d0,117.0d0,psi(1)-180.0d0,
     &                     ci(i),cai(i),ni(i),0)
            else
               call zatom (ontyp(k),1.22d0,122.5d0,psi(1)-180.0d0,
     &                     ci(i),cai(i),ni(i),0)
            end if
            call zatom (hnntyp(k),1.02d0,109.5d0,0.0d0,
     &                  ni(i),cai(i),ci(i),0)
            call zatom (hnntyp(k),1.02d0,109.5d0,-120.0d0,
     &                  ni(i),cai(i),ci(i),0)
            call zatom (hantyp(k),1.11d0,109.5d0,107.9d0,
     &                  cai(i),ni(i),ci(i),-chiral(i))
            call proside (resname,i,cai(i),ni(i),ci(i))
c
c     build the first residue as a pyroglutamic acid
c
         else if (resname .eq. 'PCA') then
            if (m .eq. 1) then
               ni(i) = n
               call zatom (nntyp(k),0.0d0,0.0d0,0.0d0,0,0,0,0)
               cai(i) = n
               call zatom (cantyp(k),1.50d0,0.0d0,0.0d0,ni(i),0,0,0)
               ci(i) = n
               if (single) then
                  call zatom (cctyp(k),1.51d0,111.6d0,0.0d0,
     &                        cai(i),ni(i),0,0)
               else
                  call zatom (cntyp(k),1.51d0,110.7d0,0.0d0,
     &                        cai(i),ni(i),0,0)
               end if
            else
               ni(i) = n
               call zatom (nntyp(k),30.0d0,150.0d0,180.0d0,
     &                     n-1,n-2,n-3,0)
               call zatom (-2,0.0d0,0.0d0,0.0d0,n-2,n-1,0,0)
               cai(i) = n
               call zatom (cantyp(k),1.50d0,150.0d0,180.0d0,
     &                     ni(i),n-2,n-3,0)
               ci(i) = n
               if (single) then
                  call zatom (cctyp(k),1.51d0,111.6d0,180.0d0,
     &                        cai(i),ni(i),n-3,0)
               else
                  call zatom (cntyp(k),1.51d0,110.7d0,180.0d0,
     &                        cai(i),ni(i),n-3,0)
               end if
            end if
            if (single) then
               call zatom (octyp(k),1.25d0,117.0d0,psi(1)-180.0d0,
     &                     ci(i),cai(i),ni(i),0)
            else
               call zatom (ontyp(k),1.22d0,122.5d0,psi(1)-180.0d0,
     &                     ci(i),cai(i),ni(i),0)
            end if
            call zatom (hnntyp(k),1.02d0,109.5d0,-60.0d0,
     &                  ni(i),cai(i),ci(i),0)
            call zatom (hantyp(k),1.11d0,109.5d0,107.9d0,
     &                  cai(i),ni(i),ci(i),-chiral(i))
            call proside (resname,i,cai(i),ni(i),ci(i))
c
c     build the first residue for N-terminal deprotonated amino acids
c
         else if (resname .eq. 'H2N') then
            i = i + 1
            k = seqtyp(i)
            resname = amino(k)
            if (m .eq. 1) then
               ni(i) = n
               k = seqtyp(i-1)
               call zatom (nntyp(k),0.0d0,0.0d0,0.0d0,0,0,0,0)
               k = seqtyp(i)
               cai(i) = n
               call zatom (cantyp(k),1.50d0,0.0d0,0.0d0,ni(i),0,0,0)
               ci(i) = n
               if (single) then
                  call zatom (cctyp(k),1.51d0,111.6d0,0.0d0,
     &                        cai(i),ni(i),0,0)
               else
                  call zatom (cntyp(k),1.51d0,110.7d0,0.0d0,
     &                        cai(i),ni(i),0,0)
               end if
            else
               ni(i) = n
               k = seqtyp(i-1)
               call zatom (nntyp(k),30.0d0,150.0d0,180.0d0,
     &                     n-1,n-2,n-3,0)
               k = seqtyp(i)
               call zatom (-2,0.0d0,0.0d0,0.0d0,n-2,n-1,0,0)
               cai(i) = n
               call zatom (cantyp(k),1.50d0,150.0d0,180.0d0,
     &                     ni(i),n-2,n-3,0)
               ci(i) = n
               if (single) then
                  call zatom (cctyp(k),1.51d0,111.6d0,180.0d0,
     &                        cai(i),ni(i),n-3,0)
               else
                  call zatom (cntyp(k),1.51d0,110.7d0,180.0d0,
     &                        cai(i),ni(i),n-3,0)
               end if
            end if
            if (single) then
               call zatom (octyp(k),1.25d0,117.0d0,psi(1)-180.0d0,
     &                     ci(i),cai(i),ni(i),0)
            else
               call zatom (ontyp(k),1.22d0,122.5d0,psi(1)-180.0d0,
     &                     ci(i),cai(i),ni(i),0)
            end if
            k = seqtyp(i-1)
            call zatom (hnntyp(k),1.02d0,109.5d0,phi(i),
     &                  ni(i),cai(i),ci(i),0)
            call zatom (hnntyp(k),1.02d0,109.5d0,108.0d0,
     &                  ni(i),cai(i),n-1,1)
            k = seqtyp(i)
            call zatom (hantyp(k),1.11d0,109.5d0,107.9d0,
     &                  cai(i),ni(i),ci(i),-chiral(i))
            call proside (resname,i,cai(i),ni(i),ci(i))
c
c     build the first residue for all other standard amino acids
c
         else
            if (m .eq. 1) then
               ni(i) = n
               call zatom (nntyp(k),0.0d0,0.0d0,0.0d0,0,0,0,0)
               cai(i) = n
               call zatom (cantyp(k),1.50d0,0.0d0,0.0d0,ni(i),0,0,0)
               ci(i) = n
               if (single) then
                  call zatom (cctyp(k),1.51d0,111.6d0,0.0d0,
     &                        cai(i),ni(i),0,0)
               else
                  call zatom (cntyp(k),1.51d0,110.7d0,0.0d0,
     &                        cai(i),ni(i),0,0)
               end if
            else
               ni(i) = n
               call zatom (nntyp(k),30.0d0,150.0d0,180.0d0,
     &                     n-1,n-2,n-3,0)
               call zatom (-2,0.0d0,0.0d0,0.0d0,n-2,n-1,0,0)
               cai(i) = n
               call zatom (cantyp(k),1.50d0,150.0d0,180.0d0,
     &                     ni(i),n-2,n-3,0)
               ci(i) = n
               if (single) then
                  call zatom (cctyp(k),1.51d0,111.6d0,180.0d0,
     &                        cai(i),ni(i),n-3,0)
               else
                  call zatom (cntyp(k),1.51d0,110.7d0,180.0d0,
     &                        cai(i),ni(i),n-3,0)
               end if
            end if
            if (single) then
               call zatom (octyp(k),1.25d0,117.0d0,psi(1)-180.0d0,
     &                     ci(i),cai(i),ni(i),0)
            else
               call zatom (ontyp(k),1.22d0,122.5d0,psi(1)-180.0d0,
     &                     ci(i),cai(i),ni(i),0)
            end if
            call zatom (hnntyp(k),1.02d0,109.5d0,phi(i),
     &                  ni(i),cai(i),ci(i),0)
            call zatom (hnntyp(k),1.02d0,109.5d0,108.0d0,
     &                  ni(i),cai(i),n-1,1)
            call zatom (hnntyp(k),1.02d0,109.5d0,108.0d0,
     &                  ni(i),cai(i),n-2,-1)
            call zatom (hantyp(k),1.11d0,109.5d0,107.9d0,
     &                  cai(i),ni(i),ci(i),-chiral(i))
            call proside (resname,i,cai(i),ni(i),ci(i))
         end if
c
c     build atoms for residues in the middle of the chain
c
         do while (i .lt. ichain(2,m)-1)
            i = i + 1
            k = seqtyp(i)
            resname = amino(k)
            ni(i) = n
            call zatom (ntyp(k),1.34d0,112.7d0,psi(i-1),
     &                  ci(i-1),cai(i-1),ni(i-1),0)
            cai(i) = n
            call zatom (catyp(k),1.46d0,121.0d0,omega(i-1),
     &                  ni(i),ci(i-1),cai(i-1),0)
            ci(i) = n
            call zatom (ctyp(k),1.51d0,111.6d0,phi(i),
     &                  cai(i),ni(i),ci(i-1),0)
            call zatom (otyp(k),1.22d0,122.5d0,psi(i)-180.0d0,
     &                  ci(i),cai(i),ni(i),0)
            call zatom (hntyp(k),1.02d0,121.0d0,phi(i)-180.0d0,
     &                  ni(i),cai(i),ci(i),0)
            call zatom (hatyp(k),1.11d0,109.5d0,107.9d0,
     &                  cai(i),ni(i),ci(i),-chiral(i))
            call proside (resname,i,cai(i),ni(i),ci(i))
         end do
c
c     set the number and type of the last residue
c
         i = ichain(2,m)
         k = seqtyp(i)
         resname = amino(k)
c
c     build the last residue for a cyclic peptide
c
         if (cyclic) then
            ni(i) = n
            call zatom (ntyp(k),1.34d0,112.7d0,psi(i-1),
     &                  ci(i-1),cai(i-1),ni(i-1),0)
            cai(i) = n
            call zatom (catyp(k),1.46d0,121.0d0,omega(i-1),
     &                  ni(i),ci(i-1),cai(i-1),0)
            ci(i) = n
            call zatom (ctyp(k),1.51d0,111.6d0,phi(i),
     &                  cai(i),ni(i),ci(i-1),0)
            call zatom (-1,0.0d0,0.0d0,0.0d0,ni(1),ci(i),0,0)
            call zatom (otyp(k),1.22d0,122.5d0,psi(i)-180.0d0,
     &                  ci(i),cai(i),ni(i),0)
            call zatom (hntyp(k),1.02d0,121.0d0,phi(i)-180.0d0,
     &                  ni(i),cai(i),ci(i),0)
            call zatom (hatyp(k),1.11d0,109.5d0,107.9d0,
     &                  cai(i),ni(i),ci(i),-chiral(i))
            call proside (resname,i,cai(i),ni(i),ci(i))
c
c     build the last residue as a C-terminal amide
c
         else if (resname .eq. 'NH2') then
            call zatom (nctyp(k),1.34d0,112.7d0,psi(i-1),
     &                  ci(i-1),cai(i-1),ni(i-1),0)
            call zatom (hnctyp(k),1.02d0,119.0d0,0.0d0,
     &                  n-1,ci(i-1),cai(i-1),0)
            call zatom (hnctyp(k),1.02d0,119.0d0,180.0d0,
     &                  n-2,ci(i-1),cai(i-1),0)
c
c     build the last residue as a C-terminal N-methylamide
c
         else if (resname .eq. 'NME') then
            call zatom (nctyp(k),1.34d0,112.7d0,psi(i-1),
     &                  ci(i-1),cai(i-1),ni(i-1),0)
            call zatom (cactyp(k),1.46d0,121.0d0,180.0d0,
     &                  n-1,ci(i-1),cai(i-1),0)
            call zatom (hnctyp(k),1.02d0,118.0d0,121.0d0,
     &                  n-2,ci(i-1),n-1,1)
            call zatom (hactyp(k),1.11d0,109.5d0,180.0d0,
     &                  n-2,n-3,ci(i-1),0)
            call zatom (hactyp(k),1.11d0,109.5d0,109.5d0,
     &                  n-3,n-4,n-1,1)
            call zatom (hactyp(k),1.11d0,109.5d0,109.5d0,
     &                  n-4,n-5,n-2,-1)
c
c     build the last residue as a protonated C-terminal amino acid
c
         else if (resname .eq. 'COH') then
            nsave = n
            n = ci(i-1)
            call zatom (cctyp(k),1.51d0,111.6d0,phi(i),
     &                  cai(i-1),ni(i-1),ci(i-2),0)
            call zatom (octyp(k),1.22d0,122.5d0,psi(i)-180.0d0,
     &                  ci(i-1),cai(i-1),ni(i-1),0)
            n = nsave
            call zatom (nctyp(k),1.35d0,112.7d0,psi(i-1),
     &                  ci(i-1),cai(i-1),ni(i-1),0)
            call zatom (hnctyp(k),0.98d0,108.7d0,180.0d0,
     &                  n-1,ci(i-1),cai(i-1),0)
c
c     build the last residue for all other standard amino acids
c
         else
            if (.not. single) then
               ni(i) = n
               call zatom (nctyp(k),1.34d0,112.7d0,psi(i-1),
     &                     ci(i-1),cai(i-1),ni(i-1),0)
               cai(i) = n
               call zatom (cactyp(k),1.46d0,121.0d0,omega(i-1),
     &                     ni(i),ci(i-1),cai(i-1),0)
               ci(i) = n
               call zatom (cctyp(k),1.51d0,111.6d0,phi(i),
     &                     cai(i),ni(i),ci(i-1),0)
               call zatom (octyp(k),1.25d0,117.0d0,psi(i)-180.0d0,
     &                     ci(i),cai(i),ni(i),0)
               call zatom (hnctyp(k),1.02d0,121.0d0,phi(i)-180.0d0,
     &                     ni(i),cai(i),ci(i),0)
               call zatom (hactyp(k),1.11d0,109.5d0,107.9d0,
     &                     cai(i),ni(i),ci(i),-chiral(i))
               call proside (resname,i,cai(i),ni(i),ci(i))
            end if
            call zatom (octyp(k),1.25d0,117.0d0,psi(i),
     &                  ci(i),cai(i),ni(i),0)
         end if
      end do
c
c     finally, set the total number of atoms
c
      n = n - 1
c
c     perform deallocation of some local arrays
c
      deallocate (ni)
      deallocate (cai)
      deallocate (ci)
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine proside  --  build amino acid side chain  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "proside" builds the side chain for a single amino acid
c     residue in terms of internal coordinates
c
c     resname   3-letter name of current amino acid residue
c     i         number of the current amino acid residue
c     cai       atom number of alpha carbon in residue i
c     ni        atom number of amide nitrogen in residue i
c     ci        atom number of carbonyl carbon in residue i
c
c     note biotypes of CD and HD atoms for N-terminal proline
c     are set as absolute values, not relative to the CB atom
c
c
      subroutine proside (resname,i,cai,ni,ci)
      use sizes
      use atoms
      use phipsi
      use resdue
      use sequen
      implicit none
      integer i,k
      integer cai,ni,ci
      integer ntprocd
      integer ntprohd
      character*3 resname
c
c
c     set the CB atom as reference site
c
      k = cbtyp(seqtyp(i))
c
c     set biotypes for CD and HD of N-terminal PRO residue
c
      ntprocd = 469
      ntprohd = 470
c
c     glycine residue  (GLY)
c
      if (resname .eq. 'GLY') then
         k = hatyp(seqtyp(i))
         if (i .eq. 1)  k = hantyp(seqtyp(i))
         if (i .eq. nseq)  k = hactyp(seqtyp(i))
         call zatom (k,1.11d0,109.5d0,107.9d0,cai,ni,ci,chiral(i))
c
c     alanine residue  (ALA)
c
      else if (resname .eq. 'ALA') then
         call zatom (k,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (k+1,1.11d0,109.4d0,chi(1,i),n-1,cai,ni,0)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-2,cai,n-1,1)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-3,cai,n-2,-1)
c
c     valine residue  (VAL)
c
      else if (resname .eq. 'VAL') then
         call zatom (k,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (k+2,1.54d0,109.5d0,chi(1,i),n-1,cai,ni,0)
         call zatom (k+4,1.54d0,109.5d0,109.5d0,n-2,cai,n-1,-1)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-3,cai,n-2,1)
         call zatom (k+3,1.11d0,109.4d0,180.0d0,n-3,n-4,cai,0)
         call zatom (k+3,1.11d0,109.4d0,109.4d0,n-4,n-5,n-1,1)
         call zatom (k+3,1.11d0,109.4d0,109.4d0,n-5,n-6,n-2,-1)
         call zatom (k+5,1.11d0,109.4d0,180.0d0,n-5,n-7,cai,0)
         call zatom (k+5,1.11d0,109.4d0,109.4d0,n-6,n-8,n-1,1)
         call zatom (k+5,1.11d0,109.4d0,109.4d0,n-7,n-9,n-2,-1)
c
c     leucine residue  (LEU)
c
      else if (resname .eq. 'LEU') then
         call zatom (k,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (k+2,1.54d0,109.5d0,chi(1,i),n-1,cai,ni,0)
         call zatom (k+4,1.54d0,109.5d0,chi(2,i),n-1,n-2,cai,0)
         call zatom (k+6,1.54d0,109.5d0,109.4d0,n-2,n-3,n-1,-1)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-4,cai,n-3,1)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-5,cai,n-4,-1)
         call zatom (k+3,1.11d0,109.4d0,109.4d0,n-5,n-6,n-4,1)
         call zatom (k+5,1.11d0,109.4d0,180.0d0,n-5,n-6,n-7,0)
         call zatom (k+5,1.11d0,109.4d0,109.4d0,n-6,n-7,n-1,1)
         call zatom (k+5,1.11d0,109.4d0,109.4d0,n-7,n-8,n-2,-1)
         call zatom (k+7,1.11d0,109.4d0,180.0d0,n-7,n-9,n-10,0)
         call zatom (k+7,1.11d0,109.4d0,109.4d0,n-8,n-10,n-1,1)
         call zatom (k+7,1.11d0,109.4d0,109.4d0,n-9,n-11,n-2,-1)
c
c     isoleucine residue  (ILE)
c
      else if (resname .eq. 'ILE') then
         call zatom (k,1.54d0,109.5d0,109.5d0,cai,ni,ci,chiral(i))
         call zatom (k+2,1.54d0,109.5d0,chi(1,i),n-1,cai,ni,0)
         call zatom (k+4,1.54d0,109.5d0,109.5d0,n-2,cai,n-1,1)
         call zatom (k+6,1.54d0,109.5d0,chi(2,i),n-2,n-3,cai,0)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-4,cai,n-3,-1)
         call zatom (k+3,1.11d0,109.4d0,109.4d0,n-4,n-5,n-2,1)
         call zatom (k+3,1.11d0,109.4d0,109.4d0,n-5,n-6,n-3,-1)
         call zatom (k+5,1.11d0,110.0d0,180.0d0,n-5,n-7,n-6,0)
         call zatom (k+5,1.11d0,110.0d0,109.0d0,n-6,n-8,n-1,1)
         call zatom (k+5,1.11d0,110.0d0,109.0d0,n-7,n-9,n-2,-1)
         call zatom (k+7,1.11d0,110.0d0,180.0d0,n-7,n-9,n-10,0)
         call zatom (k+7,1.11d0,110.0d0,109.0d0,n-8,n-10,n-1,1)
         call zatom (k+7,1.11d0,110.0d0,109.0d0,n-9,n-11,n-2,-1)
c
c     serine residue  (SER)
c
      else if (resname .eq. 'SER') then
         call zatom (k,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (k+2,1.41d0,107.5d0,chi(1,i),n-1,cai,ni,0)
         call zatom (k+1,1.11d0,109.4d0,106.7d0,n-2,cai,n-1,1)
         call zatom (k+1,1.11d0,109.4d0,106.7d0,n-3,cai,n-2,-1)
         call zatom (k+3,0.94d0,106.9d0,chi(2,i),n-3,n-4,cai,0)
c
c     threonine residue  (THR)
c
      else if (resname .eq. 'THR') then
         call zatom (k,1.54d0,109.5d0,109.5d0,cai,ni,ci,chiral(i))
         call zatom (k+2,1.41d0,107.5d0,chi(1,i),n-1,cai,ni,0)
         call zatom (k+4,1.54d0,109.5d0,107.7d0,n-2,cai,n-1,1)
         call zatom (k+1,1.11d0,109.4d0,106.7d0,n-3,cai,n-2,-1)
         call zatom (k+3,0.94d0,106.9d0,chi(2,i),n-3,n-4,cai,0)
         call zatom (k+5,1.11d0,110.0d0,180.0d0,n-3,n-5,cai,0)
         call zatom (k+5,1.11d0,110.0d0,109.0d0,n-4,n-6,n-1,1)
         call zatom (k+5,1.11d0,110.0d0,109.0d0,n-5,n-7,n-2,-1)
c
c     cysteine residue  (CYS)
c
      else if (resname .eq. 'CYS') then
         call zatom (k,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (k+2,1.82d0,109.0d0,chi(1,i),n-1,cai,ni,0)
         call zatom (k+1,1.11d0,109.4d0,112.0d0,n-2,cai,n-1,1)
         call zatom (k+1,1.11d0,109.4d0,112.0d0,n-3,cai,n-2,-1)
         call zatom (k+3,1.34d0,96.0d0,chi(2,i),n-3,n-4,cai,0)
c
c     cystine residue  (CYX)
c
      else if (resname .eq. 'CYX') then
         call zatom (k,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (k+2,1.82d0,109.0d0,chi(1,i),n-1,cai,ni,0)
         call zatom (k+1,1.11d0,109.4d0,112.0d0,n-2,cai,n-1,1)
         call zatom (k+1,1.11d0,109.4d0,112.0d0,n-3,cai,n-2,-1)
         if (disulf(i) .gt. i) then
            disulf(i) = n - 3
         else if (disulf(i) .lt. i) then
            call zatom (-1,0.0d0,0.0d0,0.0d0,disulf(disulf(i)),n-3,0,0)
         end if
c
c     deprotonated cysteine residue  (CYD)
c
      else if (resname .eq. 'CYD') then
         call zatom (k,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (k+2,1.82d0,109.0d0,chi(1,i),n-1,cai,ni,0)
         call zatom (k+1,1.11d0,109.4d0,112.0d0,n-2,cai,n-1,1)
         call zatom (k+1,1.11d0,109.4d0,112.0d0,n-3,cai,n-2,-1)
c
c     proline residue  (PRO)
c
      else if (resname .eq. 'PRO') then
         call zatom (k,1.54d0,107.0d0,109.5d0,cai,ni,ci,chiral(i))
         call zatom (k+2,1.54d0,107.0d0,chi(1,i),n-1,cai,ni,0)
         if (i .eq. 1) then
            call zatom (ntprocd,1.54d0,107.0d0,chi(2,i),n-1,n-2,cai,0)
         else
            call zatom (k+4,1.54d0,107.0d0,chi(2,i),n-1,n-2,cai,0)
         end if
         call zatom (-1,0.0d0,0.0d0,0.0d0,ni,n-1,0,0)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-3,cai,n-2,1)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-4,cai,n-3,-1)
         call zatom (k+3,1.11d0,109.4d0,109.4d0,n-4,n-5,n-3,1)
         call zatom (k+3,1.11d0,109.4d0,109.4d0,n-5,n-6,n-4,-1)
         if (i .eq. 1) then
            call zatom (ntprohd,1.11d0,109.4d0,109.4d0,n-5,n-6,ni,1)
            call zatom (ntprohd,1.11d0,109.4d0,109.4d0,n-6,n-7,ni,-1)
         else
            call zatom (k+5,1.11d0,109.4d0,109.4d0,n-5,n-6,ni,1)
            call zatom (k+5,1.11d0,109.4d0,109.4d0,n-6,n-7,ni,-1)
         end if
c
c     phenylalanine residue  (PHE)
c
      else if (resname .eq. 'PHE') then
         call zatom (k,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (k+2,1.50d0,109.5d0,chi(1,i),n-1,cai,ni,0)
         call zatom (k+3,1.39d0,120.0d0,chi(2,i),n-1,n-2,cai,0)
         call zatom (k+3,1.39d0,120.0d0,120.0d0,n-2,n-3,n-1,1)
         call zatom (k+5,1.39d0,120.0d0,180.0d0,n-2,n-3,n-4,0)
         call zatom (k+5,1.39d0,120.0d0,180.0d0,n-2,n-4,n-5,0)
         call zatom (k+7,1.39d0,120.0d0,0.0d0,n-2,n-4,n-5,0)
         call zatom (-1,0.0d0,0.0d0,0.0d0,n-2,n-1,0,0)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-7,cai,n-6,1)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-8,cai,n-7,-1)
         call zatom (k+4,1.10d0,120.0d0,120.0d0,n-7,n-8,n-5,1)
         call zatom (k+4,1.10d0,120.0d0,120.0d0,n-7,n-9,n-5,1)
         call zatom (k+6,1.10d0,120.0d0,120.0d0,n-7,n-9,n-5,1)
         call zatom (k+6,1.10d0,120.0d0,120.0d0,n-7,n-9,n-6,1)
         call zatom (k+8,1.10d0,120.0d0,120.0d0,n-7,n-9,n-8,1)
c
c     tyrosine residue  (TYR)
c
      else if (resname .eq. 'TYR') then
         call zatom (k,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (k+2,1.50d0,109.5d0,chi(1,i),n-1,cai,ni,0)
         call zatom (k+3,1.39d0,120.0d0,chi(2,i),n-1,n-2,cai,0)
         call zatom (k+3,1.39d0,120.0d0,120.0d0,n-2,n-3,n-1,1)
         call zatom (k+5,1.39d0,120.0d0,180.0d0,n-2,n-3,n-4,0)
         call zatom (k+5,1.39d0,120.0d0,180.0d0,n-2,n-4,n-5,0)
         call zatom (k+7,1.39d0,120.0d0,0.0d0,n-2,n-4,n-5,0)
         call zatom (-1,0.0d0,0.0d0,0.0d0,n-2,n-1,0,0)
         call zatom (k+8,1.36d0,120.0d0,120.0d0,n-1,n-2,n-3,1)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-8,cai,n-7,1)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-9,cai,n-8,-1)
         call zatom (k+4,1.10d0,120.0d0,120.0d0,n-8,n-9,n-6,1)
         call zatom (k+4,1.10d0,120.0d0,120.0d0,n-8,n-10,n-6,1)
         call zatom (k+6,1.10d0,120.0d0,120.0d0,n-8,n-10,n-6,1)
         call zatom (k+6,1.10d0,120.0d0,120.0d0,n-8,n-10,n-7,1)
         call zatom (k+9,0.97d0,108.0d0,0.0d0,n-7,n-8,n-9,0)
c
c     deprotonated tyrosine residue  (TYD)
c
      else if (resname .eq. 'TYD') then
         call zatom (k,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (k+2,1.50d0,109.5d0,chi(1,i),n-1,cai,ni,0)
         call zatom (k+3,1.39d0,120.0d0,chi(2,i),n-1,n-2,cai,0)
         call zatom (k+3,1.39d0,120.0d0,120.0d0,n-2,n-3,n-1,1)
         call zatom (k+5,1.39d0,120.0d0,180.0d0,n-2,n-3,n-4,0)
         call zatom (k+5,1.39d0,120.0d0,180.0d0,n-2,n-4,n-5,0)
         call zatom (k+7,1.39d0,120.0d0,0.0d0,n-2,n-4,n-5,0)
         call zatom (-1,0.0d0,0.0d0,0.0d0,n-2,n-1,0,0)
         call zatom (k+8,1.36d0,120.0d0,120.0d0,n-1,n-2,n-3,1)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-8,cai,n-7,1)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-9,cai,n-8,-1)
         call zatom (k+4,1.10d0,120.0d0,120.0d0,n-8,n-9,n-6,1)
         call zatom (k+4,1.10d0,120.0d0,120.0d0,n-8,n-10,n-6,1)
         call zatom (k+6,1.10d0,120.0d0,120.0d0,n-8,n-10,n-6,1)
         call zatom (k+6,1.10d0,120.0d0,120.0d0,n-8,n-10,n-7,1)
c
c     tryptophan residue  (TRP)
c
      else if (resname .eq. 'TRP') then
         call zatom (k,1.54d0,109.5d0,109.5d0,cai,ni,ci,chiral(i))
         call zatom (k+2,1.50d0,109.5d0,chi(1,i),n-1,cai,ni,0)
         call zatom (k+3,1.35d0,126.0d0,chi(2,i),n-1,n-2,cai,0)
         call zatom (k+5,1.35d0,126.0d0,108.0d0,n-2,n-3,n-1,1)
         call zatom (k+6,1.35d0,108.0d0,0.0d0,n-2,n-3,n-1,0)
         call zatom (k+8,1.35d0,108.0d0,0.0d0,n-1,n-3,n-4,0)
         call zatom (-1,0.0d0,0.0d0,0.0d0,n-3,n-1,0,0)
         call zatom (k+9,1.35d0,120.0d0,180.0d0,n-3,n-1,n-2,0)
         call zatom (k+11,1.35d0,120.0d0,0.0d0,n-2,n-4,n-1,0)
         call zatom (k+13,1.35d0,120.0d0,0.0d0,n-2,n-5,n-3,0)
         call zatom (k+15,1.35d0,120.0d0,0.0d0,n-2,n-4,n-6,0)
         call zatom (-1,0.0d0,0.0d0,0.0d0,n-2,n-1,0,0)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-10,cai,n-9,1)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-11,cai,n-10,-1)
         call zatom (k+4,1.10d0,126.0d0,126.0d0,n-10,n-11,n-8,1)
         call zatom (k+7,1.05d0,126.0d0,126.0d0,n-9,n-11,n-8,1)
         call zatom (k+10,1.10d0,120.0d0,120.0d0,n-8,n-11,n-6,1)
         call zatom (k+12,1.10d0,120.0d0,120.0d0,n-8,n-10,n-6,1)
         call zatom (k+14,1.10d0,120.0d0,120.0d0,n-8,n-10,n-7,1)
         call zatom (k+16,1.10d0,120.0d0,120.0d0,n-8,n-10,n-9,1)
c
c     histidine (HD and HE) residue  (HIS)
c
      else if (resname .eq. 'HIS') then
         call zatom (k,1.54d0,109.5d0,109.5d0,cai,ni,ci,chiral(i))
         call zatom (k+2,1.50d0,109.5d0,chi(1,i),n-1,cai,ni,0)
         call zatom (k+3,1.35d0,126.0d0,chi(2,i),n-1,n-2,cai,0)
         call zatom (k+5,1.35d0,126.0d0,108.0d0,n-2,n-3,n-1,1)
         call zatom (k+7,1.35d0,108.0d0,0.0d0,n-2,n-3,n-1,0)
         call zatom (k+9,1.35d0,108.0d0,0.0d0,n-2,n-4,n-3,0)
         call zatom (-1,0.0d0,0.0d0,.00d0,n-2,n-1,0,0)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-6,cai,n-5,1)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-7,cai,n-6,-1)
         call zatom (k+4,1.02d0,126.0d0,0.0d0,n-6,n-7,n-8,0)
         call zatom (k+6,1.10d0,126.0d0,126.0d0,n-6,n-8,n-4,1)
         call zatom (k+8,1.10d0,126.0d0,126.0d0,n-6,n-8,n-5,1)
         call zatom (k+10,1.02d0,126.0d0,126.0d0,n-6,n-8,n-7,1)
c
c     histidine (HD only) residue  (HID)
c
      else if (resname .eq. 'HID') then
         call zatom (k,1.54d0,109.5d0,109.5d0,cai,ni,ci,chiral(i))
         call zatom (k+2,1.50d0,109.5d0,chi(1,i),n-1,cai,ni,0)
         call zatom (k+3,1.35d0,126.0d0,chi(2,i),n-1,n-2,cai,0)
         call zatom (k+5,1.35d0,126.0d0,108.0d0,n-2,n-3,n-1,1)
         call zatom (k+7,1.35d0,108.0d0,0.0d0,n-2,n-3,n-1,0)
         call zatom (k+9,1.35d0,108.0d0,0.0d0,n-2,n-4,n-3,0)
         call zatom (-1,0.0d0,0.0d0,.00d0,n-2,n-1,0,0)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-6,cai,n-5,1)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-7,cai,n-6,-1)
         call zatom (k+4,1.02d0,126.0d0,0.0d0,n-6,n-7,n-8,0)
         call zatom (k+6,1.10d0,126.0d0,126.0d0,n-6,n-8,n-4,1)
         call zatom (k+8,1.10d0,126.0d0,126.0d0,n-6,n-8,n-5,1)
c
c     histidine (HE only) residue  (HIE)
c
      else if (resname .eq. 'HIE') then
         call zatom (k,1.54d0,109.5d0,109.5d0,cai,ni,ci,chiral(i))
         call zatom (k+2,1.50d0,109.5d0,chi(1,i),n-1,cai,ni,0)
         call zatom (k+3,1.35d0,126.0d0,chi(2,i),n-1,n-2,cai,0)
         call zatom (k+4,1.35d0,126.0d0,108.0d0,n-2,n-3,n-1,1)
         call zatom (k+6,1.35d0,108.0d0,0.0d0,n-2,n-3,n-1,0)
         call zatom (k+8,1.35d0,108.0d0,0.0d0,n-2,n-4,n-3,0)
         call zatom (-1,0.0d0,0.0d0,.00d0,n-2,n-1,0,0)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-6,cai,n-5,1)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-7,cai,n-6,-1)
         call zatom (k+5,1.10d0,126.0d0,126.0d0,n-5,n-7,n-3,1)
         call zatom (k+7,1.10d0,126.0d0,126.0d0,n-5,n-7,n-4,1)
         call zatom (k+9,1.02d0,126.0d0,126.0d0,n-5,n-7,n-6,1)
c
c     aspartic acid residue  (ASP)
c
      else if (resname .eq. 'ASP') then
         call zatom (k,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (k+2,1.51d0,107.8d0,chi(1,i),n-1,cai,ni,0)
         call zatom (k+3,1.25d0,117.0d0,chi(2,i),n-1,n-2,cai,0)
         call zatom (k+3,1.25d0,117.0d0,126.0d0,n-2,n-3,n-1,1)
         call zatom (k+1,1.11d0,109.4d0,107.9d0,n-4,cai,n-3,1)
         call zatom (k+1,1.11d0,109.4d0,107.9d0,n-5,cai,n-4,-1)
c
c     protonated aspartic acid residue  (ASH)
c
      else if (resname .eq. 'ASH') then
         call zatom (k,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (k+2,1.51d0,107.8d0,chi(1,i),n-1,cai,ni,0)
         call zatom (k+3,1.25d0,117.0d0,chi(2,i),n-1,n-2,cai,0)
         call zatom (k+4,1.25d0,117.0d0,126.0d0,n-2,n-3,n-1,1)
         call zatom (k+1,1.11d0,109.4d0,107.9d0,n-4,cai,n-3,1)
         call zatom (k+1,1.11d0,109.4d0,107.9d0,n-5,cai,n-4,-1)
         call zatom (k+5,0.98d0,108.7d0,0.0d0,n-3,n-5,n-4,0)
c
c     asparagine residue  (ASN)
c
      else if (resname .eq. 'ASN') then
         call zatom (k,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (k+2,1.51d0,107.8d0,chi(1,i),n-1,cai,ni,0)
         call zatom (k+3,1.22d0,122.5d0,chi(2,i),n-1,n-2,cai,0)
         call zatom (k+4,1.34d0,112.7d0,124.0d0,n-2,n-3,n-1,1)
         call zatom (k+1,1.11d0,109.4d0,107.9d0,n-4,cai,n-3,1)
         call zatom (k+1,1.11d0,109.4d0,107.9d0,n-5,cai,n-4,-1)
         call zatom (k+5,1.02d0,119.0d0,0.0d0,n-3,n-5,n-6,0)
         call zatom (k+5,1.02d0,119.0d0,120.0d0,n-4,n-6,n-1,1)
c
c     glutamic acid residue  (GLU)
c
      else if (resname .eq. 'GLU') then
         call zatom (k,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (k+2,1.54d0,109.5d0,chi(1,i),n-1,cai,ni,0)
         call zatom (k+4,1.51d0,107.8d0,chi(2,i),n-1,n-2,cai,0)
         call zatom (k+5,1.25d0,117.0d0,chi(3,i),n-1,n-2,n-3,0)
         call zatom (k+5,1.25d0,117.0d0,126.0d0,n-2,n-3,n-1,1)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-5,cai,n-4,1)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-6,cai,n-5,-1)
         call zatom (k+3,1.11d0,109.4d0,107.9d0,n-6,n-7,n-5,1)
         call zatom (k+3,1.11d0,109.4d0,107.9d0,n-7,n-8,n-6,-1)
c
c     protonated glutamic acid residue  (GLH)
c
      else if (resname .eq. 'GLH') then
         call zatom (k,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (k+2,1.54d0,109.5d0,chi(1,i),n-1,cai,ni,0)
         call zatom (k+4,1.51d0,107.8d0,chi(2,i),n-1,n-2,cai,0)
         call zatom (k+5,1.25d0,117.0d0,chi(3,i),n-1,n-2,n-3,0)
         call zatom (k+6,1.25d0,117.0d0,126.0d0,n-2,n-3,n-1,1)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-5,cai,n-4,1)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-6,cai,n-5,-1)
         call zatom (k+3,1.11d0,109.4d0,107.9d0,n-6,n-7,n-5,1)
         call zatom (k+3,1.11d0,109.4d0,107.9d0,n-7,n-8,n-6,-1)
         call zatom (k+7,0.98d0,108.7d0,0.0d0,n-5,n-7,n-6,0)
c
c     glutamine residue  (GLN)
c
      else if (resname .eq. 'GLN') then
         call zatom (k,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (k+2,1.54d0,109.5d0,chi(1,i),n-1,cai,ni,0)
         call zatom (k+4,1.51d0,107.8d0,chi(2,i),n-1,n-2,cai,0)
         call zatom (k+5,1.22d0,122.5d0,chi(3,i),n-1,n-2,n-3,0)
         call zatom (k+6,1.34d0,112.7d0,124.0d0,n-2,n-3,n-1,1)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-5,cai,n-4,1)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-6,cai,n-5,-1)
         call zatom (k+3,1.11d0,109.4d0,107.9d0,n-6,n-7,n-5,1)
         call zatom (k+3,1.11d0,109.4d0,107.9d0,n-7,n-8,n-6,-1)
         call zatom (k+7,1.02d0,119.0d0,0.0d0,n-5,n-7,n-8,0)
         call zatom (k+7,1.02d0,119.0d0,120.0d0,n-6,n-8,n-1,1)
c
c     methionine residue  (MET)
c
      else if (resname .eq. 'MET') then
         call zatom (k,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (k+2,1.54d0,109.5d0,chi(1,i),n-1,cai,ni,0)
         call zatom (k+4,1.82d0,109.0d0,chi(2,i),n-1,n-2,cai,0)
         call zatom (k+5,1.82d0,96.3d0,chi(3,i),n-1,n-2,n-3,0)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-4,cai,n-3,1)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-5,cai,n-4,-1)
         call zatom (k+3,1.11d0,109.4d0,112.0d0,n-5,n-6,n-4,1)
         call zatom (k+3,1.11d0,109.4d0,112.0d0,n-6,n-7,n-5,-1)
         call zatom (k+6,1.11d0,112.0d0,180.0d0,n-5,n-6,n-7,0)
         call zatom (k+6,1.11d0,112.0d0,109.4d0,n-6,n-7,n-1,1)
         call zatom (k+6,1.11d0,112.0d0,109.4d0,n-7,n-8,n-2,-1)
c
c     lysine residue  (LYS)
c
      else if (resname .eq. 'LYS') then
         call zatom (k,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (k+2,1.54d0,109.5d0,chi(1,i),n-1,cai,ni,0)
         call zatom (k+4,1.54d0,109.5d0,chi(2,i),n-1,n-2,cai,0)
         call zatom (k+6,1.54d0,109.5d0,chi(3,i),n-1,n-2,n-3,0)
         call zatom (k+8,1.50d0,109.5d0,chi(4,i),n-1,n-2,n-3,0)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-5,cai,n-4,1)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-6,cai,n-5,-1)
         call zatom (k+3,1.11d0,109.4d0,109.4d0,n-6,n-7,n-5,1)
         call zatom (k+3,1.11d0,109.4d0,109.4d0,n-7,n-8,n-6,-1)
         call zatom (k+5,1.11d0,109.4d0,109.4d0,n-7,n-8,n-6,1)
         call zatom (k+5,1.11d0,109.4d0,109.4d0,n-8,n-9,n-7,-1)
         call zatom (k+7,1.11d0,109.4d0,108.8d0,n-8,n-9,n-7,1)
         call zatom (k+7,1.11d0,109.4d0,108.8d0,n-9,n-10,n-8,-1)
         call zatom (k+9,1.02d0,109.5d0,180.0d0,n-9,n-10,n-11,0)
         call zatom (k+9,1.02d0,109.5d0,109.5d0,n-10,n-11,n-1,1)
         call zatom (k+9,1.02d0,109.5d0,109.5d0,n-11,n-12,n-2,-1)
c
c     deprotonated lysine residue  (LYD)
c
      else if (resname .eq. 'LYD') then
         call zatom (k,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (k+2,1.54d0,109.5d0,chi(1,i),n-1,cai,ni,0)
         call zatom (k+4,1.54d0,109.5d0,chi(2,i),n-1,n-2,cai,0)
         call zatom (k+6,1.54d0,109.5d0,chi(3,i),n-1,n-2,n-3,0)
         call zatom (k+8,1.50d0,109.5d0,chi(4,i),n-1,n-2,n-3,0)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-5,cai,n-4,1)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-6,cai,n-5,-1)
         call zatom (k+3,1.11d0,109.4d0,109.4d0,n-6,n-7,n-5,1)
         call zatom (k+3,1.11d0,109.4d0,109.4d0,n-7,n-8,n-6,-1)
         call zatom (k+5,1.11d0,109.4d0,109.4d0,n-7,n-8,n-6,1)
         call zatom (k+5,1.11d0,109.4d0,109.4d0,n-8,n-9,n-7,-1)
         call zatom (k+7,1.11d0,109.4d0,108.8d0,n-8,n-9,n-7,1)
         call zatom (k+7,1.11d0,109.4d0,108.8d0,n-9,n-10,n-8,-1)
         call zatom (k+9,1.02d0,109.5d0,180.0d0,n-9,n-10,n-11,0)
         call zatom (k+9,1.02d0,109.5d0,109.5d0,n-10,n-11,n-1,1)
c
c     arginine residue  (ARG)
c
      else if (resname .eq. 'ARG') then
         call zatom (k,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (k+2,1.54d0,109.5d0,chi(1,i),n-1,cai,ni,0)
         call zatom (k+4,1.54d0,109.5d0,chi(2,i),n-1,n-2,cai,0)
         call zatom (k+6,1.45d0,109.5d0,chi(3,i),n-1,n-2,n-3,0)
         call zatom (k+8,1.35d0,120.0d0,chi(4,i),n-1,n-2,n-3,0)
         call zatom (k+9,1.35d0,120.0d0,180.0d0,n-1,n-2,n-3,0)
         call zatom (k+9,1.35d0,120.0d0,120.0d0,n-2,n-3,n-1,1)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-7,cai,n-6,1)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-8,cai,n-7,-1)
         call zatom (k+3,1.11d0,109.4d0,109.4d0,n-8,n-9,n-7,1)
         call zatom (k+3,1.11d0,109.4d0,109.4d0,n-9,n-10,n-8,-1)
         call zatom (k+5,1.11d0,109.4d0,109.4d0,n-9,n-10,n-8,1)
         call zatom (k+5,1.11d0,109.4d0,109.4d0,n-10,n-11,n-9,-1)
         call zatom (k+7,1.02d0,120.0d0,120.0d0,n-10,n-11,n-9,1)
         call zatom (k+10,1.02d0,120.0d0,180.0d0,n-9,n-10,n-11,0)
         call zatom (k+10,1.02d0,120.0d0,120.0d0,n-10,n-11,n-1,1)
         call zatom (k+10,1.02d0,120.0d0,180.0d0,n-10,n-12,n-13,0)
         call zatom (k+10,1.02d0,120.0d0,120.0d0,n-11,n-13,n-1,1)
c
c     ornithine residue  (ORN)
c
      else if (resname .eq. 'ORN') then
         call zatom (k,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (k+2,1.54d0,109.5d0,chi(1,i),n-1,cai,ni,0)
         call zatom (k+4,1.54d0,109.5d0,chi(2,i),n-1,n-2,cai,0)
         call zatom (k+6,1.50d0,109.5d0,chi(3,i),n-1,n-2,n-3,0)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-4,cai,n-3,1)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-5,cai,n-4,-1)
         call zatom (k+3,1.11d0,109.4d0,109.4d0,n-5,n-6,n-4,1)
         call zatom (k+3,1.11d0,109.4d0,109.4d0,n-6,n-7,n-5,-1)
         call zatom (k+5,1.11d0,109.4d0,109.4d0,n-6,n-7,n-5,1)
         call zatom (k+5,1.11d0,109.4d0,109.4d0,n-7,n-8,n-6,-1)
         call zatom (k+7,1.02d0,109.5d0,180.0d0,n-7,n-8,n-9,0)
         call zatom (k+7,1.02d0,109.5d0,109.5d0,n-8,n-9,n-1,1)
         call zatom (k+7,1.02d0,109.5d0,109.5d0,n-9,n-10,n-2,-1)
c
c     methylalanine residue  (AIB)
c
      else if (resname .eq. 'AIB') then
         call zatom (k,1.54d0,109.5d0,107.8d0,cai,ni,ci,-chiral(i))
         call zatom (k,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (k+1,1.11d0,109.4d0,chi(1,i),n-2,cai,ni,0)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-3,cai,n-1,1)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-4,cai,n-2,-1)
         call zatom (k+1,1.11d0,109.4d0,chi(1,i),n-4,cai,ni,0)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-5,cai,n-1,1)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-6,cai,n-2,-1)
c
c     pyroglutamic acid residue  (PCA)
c
      else if (resname .eq. 'PCA') then
         call zatom (k,1.54d0,107.0d0,109.5d0,cai,ni,ci,chiral(i))
         call zatom (k+2,1.54d0,107.0d0,chi(1,i),n-1,cai,ni,0)
         call zatom (k+4,1.54d0,107.0d0,chi(2,i),n-1,n-2,cai,0)
         call zatom (-1,0.0d0,0.0d0,0.0d0,ni,n-1,0,0)
         call zatom (k+5,1.22d0,126.0d0,126.0d0,n-1,ni,n-2,1)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-4,cai,n-3,1)
         call zatom (k+1,1.11d0,109.4d0,109.4d0,n-5,cai,n-4,-1)
         call zatom (k+3,1.11d0,109.4d0,109.4d0,n-5,n-6,n-4,1)
         call zatom (k+3,1.11d0,109.4d0,109.4d0,n-6,n-7,n-5,-1)
c
c     unknown residue  (UNK)
c
      else if (resname .eq. 'UNK') then
         k = hatyp(seqtyp(i))
         if (i .eq. 1)  k = hantyp(seqtyp(i))
         if (i .eq. nseq)  k = hactyp(seqtyp(i))
         call zatom (k,1.11d0,109.5d0,107.9d0,cai,ni,ci,chiral(i))
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine pauling  --  pack multiple polypeptide chains  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "pauling" uses a rigid body optimization to approximately
c     pack multiple polypeptide chains
c
c
      subroutine pauling
      use sizes
      use atomid
      use atoms
      use couple
      use group
      use inform
      use katoms
      use molcul
      use output
      use potent
      use restrn
      use rigid
      use usage
      implicit none
      integer i,j,k,nvar
      real*8 minimum,grdmin
      real*8 pauling1
      real*8, allocatable :: xx(:)
      external pauling1,optsave
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
      use_basin = .true.
      depth = 3.0d0
      width = 1.5d0
      use_wall = .false.
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
c     assign each chain to a separate molecule-based group
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
         wgrp(i,i) = 1.0d0
      end do
c
c     assume unit mass for each atom and set group masses
c
      do i = 1, n
         mass(i) = 1.0d0
      end do
      do i = 1, ngrp
         grpmass(i) = dble(igrp(2,i)-igrp(1,i)+1)
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (.not. allocated(ipfix))  allocate (ipfix(maxfix))
      if (.not. allocated(kpfix))  allocate (kpfix(3,maxfix))
      if (.not. allocated(xpfix))  allocate (xpfix(maxfix))
      if (.not. allocated(ypfix))  allocate (ypfix(maxfix))
      if (.not. allocated(zpfix))  allocate (zpfix(maxfix))
      if (.not. allocated(pfix))  allocate (pfix(2,maxfix))
      if (.not. allocated(igfix))  allocate (igfix(2,maxfix))
      if (.not. allocated(gfix))  allocate (gfix(3,maxfix))
c
c     set pairwise restraints between the centers of chains
c
      do i = 1, ngrp-1
         do j = i+1, ngrp
            ngfix = ngfix + 1
            igfix(1,ngfix) = i
            igfix(2,ngfix) = j
            gfix(1,ngfix) = 1.0d0
            gfix(2,ngfix) = 11.0d0 * dble(j-i)
            gfix(3,ngfix) = 11.0d0 * dble(j-i)
         end do
      end do
c
c     set position restraints on alpha carbons of each chain
c
      do i = 1, n
         if (atmnum(type(i)) .eq. 6) then
            do j = 1, n12(i)
               if (atmnum(type(i12(j,i))) .eq. 7) then
                  do k = 1, n13(i)
                     if (atmnum(type(i13(k,i))) .eq. 8) then
                        npfix = npfix + 1
                        ipfix(npfix) = i
                        kpfix(1,npfix) = 1
                        kpfix(2,npfix) = 1
                        kpfix(3,npfix) = 0
                        xpfix(npfix) = 11.0d0 * dble(grplist(i)-1)
                        ypfix(npfix) = 0.0d0
                        zpfix(npfix) = 0.0d0
                        pfix(1,npfix) = 1.0d0
                        pfix(2,npfix) = 0.0d0
                        goto 10
                     end if
                  end do
               end if
            end do
         end if
   10    continue
      end do
c
c     get rigid body reference coordinates for each chain
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
      call ocvm (nvar,xx,minimum,grdmin,pauling1,optsave)
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
c     ##############################################################
c     ##                                                          ##
c     ##  function pauling1  --  energy and gradient for pauling  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "pauling1" is a service routine that computes the energy
c     and gradient for optimally conditioned variable metric
c     optimization of rigid bodies
c
c
      function pauling1 (xx,g)
      use sizes
      use group
      use math
      use rigid
      implicit none
      integer i,j,nvar
      real*8 pauling1,e
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
      pauling1 = e
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
