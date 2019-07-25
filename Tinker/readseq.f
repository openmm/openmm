c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine readseq  --  read biopolymer sequence file  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "readseq" gets a biopolymer sequence containing one or more
c     separate chains from an external file; all lines containing
c     sequence must begin with the starting sequence number, the
c     actual sequence is read from subsequent nonblank characters
c
c
      subroutine readseq (iseq)
      use files
      use iounit
      use resdue
      use sequen
      implicit none
      integer i,j,k,iseq
      integer length,number
      integer start,stop
      integer next,trimtext
      logical exist,opened,done
      character*1 letter
      character*3 word
      character*240 seqfile
      character*240 record
c
c
c     open the input file if it has not already been done
c
      inquire (unit=iseq,opened=opened)
      if (.not. opened) then
         seqfile = filename(1:leng)//'.seq'
         call version (seqfile,'old')
         inquire (file=seqfile,exist=exist)
         if (exist) then
            open (unit=iseq,file=seqfile,status='old')
            rewind (unit=iseq)
         else
            write (iout,10)
   10       format (/,' READSEQ  --  Unable to Find the Biopolymer',
     &                 ' Sequence File')
            call fatal
         end if
      end if
c
c     zero out the number and type of residues
c
      nseq = 0
      nchain = 0
      do i = 1, maxres
         seq(i) = '   '
      end do
c
c     read in the biopolymer sequence file
c
      do while (.true.)
         read (iseq,20,err=30,end=30)  record
   20    format (a240)
         length = trimtext (record)
         next = 1
         call gettext (record,letter,next)
         if (letter.ge.'0' .and. letter.le.'9') then
            next = 1
            letter = ' '
         end if
         call getnumb (record,number,next)
         if (number .eq. 1) then
            nchain = nchain + 1
            ichain(1,nchain) = nseq + 1
            chnnam(nchain) = letter
         end if
         done = .false.
         do while (.not. done)
            call getword (record,word,next)
            if (word .eq. '   ') then
               done = .true.
            else
               nseq = nseq + 1
               seq(nseq) = word
            end if
         end do
      end do
   30 continue
c
c     set the last residue in each sequence chain
c
      do i = 1, nchain-1
         ichain(2,i) = ichain(1,i+1) - 1
      end do
      if (nchain .ne. 0)  ichain(2,nchain) = nseq
c
c     find residue types and species present in each chain
c
      do i = 1, nchain
         start = ichain(1,i)
         stop = ichain(2,i)
         chntyp(i) = 'GENERIC'
         do j = start, stop
            do k = 1, maxamino
               if (seq(j) .eq. amino(k)) then
                  seqtyp(j) = k
                  chntyp(i) = 'PEPTIDE'
                  goto 40
               end if
            end do
            chntyp(i) = 'GENERIC'
            goto 50
   40       continue
         end do
   50    continue
         if (chntyp(i) .eq. 'GENERIC') then
            do j = start, stop
               do k = 1, maxnuc
                  if (seq(j) .eq. nuclz(k)) then
                     seqtyp(j) = k
                     chntyp(i) = 'NUCLEIC'
                     goto 60
                  end if
               end do
               chntyp(i) = 'GENERIC'
               goto 70
   60          continue
            end do
   70       continue
         end if
         if (chntyp(i) .eq. 'GENERIC') then
            do j = start, stop
               seqtyp(j) = 0
            end do
         end if
      end do
      if (.not. opened)  close (unit=iseq)
      return
      end
