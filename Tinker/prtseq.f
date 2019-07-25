c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine prtseq  --  output of biopolymer sequence  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "prtseq" writes out a biopolymer sequence to an external
c     disk file with 15 residues per line and distinct chains
c     separated by blank lines
c
c
      subroutine prtseq (iseq)
      use sizes
      use files
      use sequen
      implicit none
      integer i,k,iseq
      integer smax,smin
      integer size,start,stop
      logical opened
      character*1 letter
      character*23 fstr
      character*240 seqfile
c
c
c     open the output unit if not already done
c
      inquire (unit=iseq,opened=opened)
      if (.not. opened) then
         seqfile = filename(1:leng)//'.seq'
         call version (seqfile,'new')
         open (unit=iseq,file=seqfile,status='new')
      end if
c
c     write out a three-letter code sequence file
c
      do i = 1, nchain
         letter = chnnam(i)
         start = ichain(1,i)
         stop = ichain(2,i)
         size = stop - start + 1
         smax = 0
         do while (smax .lt. size)
            smin = smax + 1
            smax = smax + 15
            smax = min(smax,size)
            if (i.ne.1 .and. smin.eq.1)  write (iseq,'()')
            fstr = '(3x,a1,i6,1x,15(1x,a3))'
            write (iseq,fstr)  letter,smin,(seq(k+start-1),k=smin,smax)
         end do
      end do
c
c     close the output unit if opened by this routine
c
      if (.not. opened)  close (unit=iseq)
      return
      end
