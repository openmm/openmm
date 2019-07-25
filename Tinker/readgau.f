c
c
c     ##############################################################
c     ##  COPYRIGHT (C) 2008 by Chuanjie Wu & Jay William Ponder  ##
c     ##                   All Rights Reserved                    ##
c     ##############################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine readgau  --  read data from G09 output file  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "readgau" reads an ab initio optimized structure, forces,
c     Hessian and frequencies from a Gaussian 09 output file
c
c
      subroutine readgau
      use sizes
      use ascii
      use iounit
      use qmstuf
      use units
      implicit none
      integer i,j
      integer igau,code
      integer ngfreq,nghess
      integer itmp,jtmp,ktmp
      integer length,next
      integer freeunit
      integer trimtext
      logical hasinputxyz
      logical hasmp2
      logical exist
      real*8 xtmp,ytmp,ztmp
      real*8 frcunit,hessunit
      character*4 arcstart
      character*6 gname
      character*240 gaufile
      character*240 record
      character*240 string
      character*240 word
c
c
c     initialize some values prior to opening the log file
c
      exist = .false.
      hasinputxyz = .false.
      ngatom = 0
      ngfreq = 0
      arcstart = '1'//char(backslash)//'1'//char(backslash)
c
c     specify and open the Gaussian 09 output log file
c
      call nextarg (gaufile,exist)
      if (exist) then
         inquire (file=gaufile,exist=exist)
         igau = freeunit()
         call basefile (gaufile)
         call suffix (gaufile,'log','old')
         inquire (file=gaufile,exist=exist)
         if (.not. exist) then
            call basefile (gaufile)
            call suffix (gaufile,'out','old')
            inquire (file=gaufile,exist=exist)
         end if
      end if
      do while (.not. exist)
         write (iout,10)
   10    format (/,' Enter the Name of the Gaussian Output File :  ',$)
         read (input,20)  gaufile
   20    format (a240)
         igau = freeunit ()
         call basefile (gaufile)
         call suffix (gaufile,'log','old')
         inquire (file=gaufile,exist=exist)
         if (.not. exist) then
            call basefile (gaufile)
            call suffix (gaufile,'out','old')
            inquire (file=gaufile,exist=exist)
         end if
      end do
c
c     scan the Gaussian output file to get the number of atoms
c
      open (unit=igau,file=gaufile,status='old')
      rewind (unit=igau)
c     do while (.true. .and. .not.eof(igau))
      do while (.true.)
         read (igau,30,err=70,end=70)  record
   30    format (a240)
         next = 1
         string = record
         call trimhead (string)
         length = trimtext (string)
         call upcase (string)
         if (string(1:20) .eq. 'STANDARD ORIENTATION') then
            do i = 1, 4
               read (igau,40,err=70,end=70)  record
   40          format (a240)
            end do
            i = 1
            do while (.true.)
               read (igau,50,err=70,end=70)  record
   50          format (a240)
               read (record,*,err=60,end=60)  itmp,jtmp,ktmp,
     &                                        xtmp,ytmp,ztmp
               if (jtmp .le. 0)  goto 60
               i = i + 1
            end do
   60       continue
            ngatom = i - 1
         end if
      end do
   70 continue
c
c     perform dynamic allocation of some global arrays
c
      nghess = (3*ngatom*(3*ngatom+1)) / 2
      if (.not. allocated(gx))  allocate (gx(ngatom))
      if (.not. allocated(gy))  allocate (gy(ngatom))
      if (.not. allocated(gz))  allocate (gz(ngatom))
      if (.not. allocated(gfreq))  allocate (gfreq(3*ngatom))
      if (.not. allocated(gforce))  allocate (gforce(3,ngatom))
      if (.not. allocated(gh))  allocate (gh(nghess))
c
c     read structure, forces and frequencies from Gaussian output
c
      rewind (unit=igau)
c     do while (.true. .and. .not.eof(igau))
      do while (.true.)
         read (igau,80,err=220,end=220)  record
   80    format (a240)
         next = 1
         string = record
         call trimhead (string)
         length = trimtext (string)
         call upcase (string)
         if (string(1:20) .eq. 'STANDARD ORIENTATION') then
            do i = 1, 4
               read (igau,90,err=220,end=220)  record
   90          format (a240)
            end do
            i = 1
            do while (.true.)
               read (igau,100,err=220,end=220)  record
  100          format (a240)
               read (record,*,err=110,end=110)  itmp,jtmp,ktmp,
     &                                          gx(i),gy(i),gz(i)
               if (jtmp .le. 0)  goto 110
               i = i + 1
            end do
  110       continue
            ngatom = i - 1
         else if (string(37:58) .eq. 'FORCES (HARTREES/BOHR)') then
            read (igau,120,err=220,end=220)  record
  120       format (a240)
            read (igau,130,err=220,end=220)  record
  130       format (a240)
            frcunit = hartree / bohr
            do i = 1, ngatom
               gforce(1,i) = 0.0d0
               gforce(2,i) = 0.0d0
               gforce(3,i) = 0.0d0
               read (igau,140,err=220,end=220)  record
  140          format (a240)
               read (record,*,err=150,end=150)  itmp,jtmp,gforce(1,i),
     &                                          gforce(2,i),gforce(3,i)
               do j = 1, 3
                  gforce(j,i) = frcunit * gforce(j,i)
               end do
  150          continue
            end do
         else if (string(1:14) .eq. 'FREQUENCIES --') then
            gfreq(ngfreq+1) = 0.0d0
            gfreq(ngfreq+2) = 0.0d0
            gfreq(ngfreq+3) = 0.0d0
            read (string(15:240),*,err=160,end=160)  gfreq(ngfreq+1),
     &                                               gfreq(ngfreq+2),
     &                                               gfreq(ngfreq+3)
  160       continue
            ngfreq = ngfreq + 3
c
c     read the Hessian from archive section at bottom of output
c
         else if (string(1:4) .eq. arcstart) then
            itmp = 0
c           do while (.true. .and. .not.eof(igau))
            do while (.true.)
               if (next .gt. 73) then
                  read (igau,170,err=220,end=220)  record
  170             format (a240)
                  next = 1
               end if
               call readgarc (igau,record,word,length,next)
               if (word(1:1) .eq. char(backslash))  itmp = itmp + 1
               if (itmp.eq.16 .and. hasinputxyz) then
                  do i = 1, ngatom
                     do j = 1, 5
                        if (next .gt. 73) then
                           read (igau,180,err=220,end=220)  record
  180                      format (a240)
                           next = 1
                        end if
                        call readgarc (igau,record,word,length,next)
                        if (j .eq. 1)  read(word(1:length),*)  gname
                        if (j .eq. 2)  read(word(1:length),*)  gx(i)
                        if (j .eq. 3)  read(word(1:length),*)  gy(i)
                        if (j .eq. 4)  read(word(1:length),*)  gz(i)
                     end do
                  end do
               end if
               if (itmp.gt.16 .and. word(1:2).eq.'HF') then
                  do i = 1, 2
                     if (next .gt. 73) then
                        read (igau,190,err=220,end=220)  record
  190                   format (a240)
                        next = 1
                     end if
                     call readgarc (igau,record,word,length,next)
                  end do
                  read (word(1:length),*)  egau
                  egau = hartree * egau
               else if (itmp.gt.16 .and. word(1:3).eq.'MP2') then
                  hasmp2 = .true.
                  do i = 1, 2
                     if (next .gt. 73) then
                        read (igau,200,err=220,end=220)  record
  200                   format (a240)
                        next = 1
                     end if
                     call readgarc (igau,record,word,length,next)
                  end do
                  read (word(1:length),*)  egau
                  egau = hartree * egau
               else if (word(1:5) .eq. 'NImag') then
                  do i = 1, 4
                     call readgarc (igau,record,word,length,next)
                  end do
                  hessunit = hartree / bohr**2
                  do i = 1, nghess
                     call readgarc (igau,record,word,length,next)
                     read (word(1:length),*)  gh(i)
                     gh(i) = hessunit * gh(i)
                  end do
                  goto 220
               end if
               code = ichar(word(1:1))
               if (code .eq. atsign)  goto 210
            end do
         end if
  210    continue
      end do
  220 continue
      close (unit=igau)
c
c     zero out the frequencies if none were in Gaussian output
c
      if (ngfreq .eq. 0) then
         do i = 1, 3*ngatom
            gfreq(i) = 0.0d0
         end do
      end if
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine readgarc  --  read Gaussian archive section  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "readgarc" reads data from Gaussian archive section; each
c     entry is terminated with a backslash symbol
c
c
      subroutine readgarc (igau,string,word,length,next)
      use ascii
      implicit none
      integer i,igau,code
      integer next,length
      character*1 letter
      character*240 word
      character*240 string
c
c
c     initialize some values prior to parsing the test string
c
      length = 1
      letter = ' '
      do i = 1, 240
         word(i:i) = ' '
      end do
c
c     attempt to read a text word entry from the input string
c
      letter = string (next:next)
      code = ichar(letter)
      if (code.eq.backslash .or. code.eq.equal
     &       .or. code.eq.space) then
         word(1:1) = letter
         next = next + 1
         length = 1
         return
      end if
   10 continue
      do i = next, 75
         if (code.eq.backslash .or. code.eq.equal
     &          .or. code.eq.space)  return
         if (next .gt. 70) then
            read (igau,20,err=30,end=30)  string
   20       format (a240)
            next = 1
            goto 10
         end if
         if (code .eq. comma) then
            next = next + 1
            return
         end if
         if (code.eq.backslash .or. code.eq.equal
     &          .or. code.eq.space)  return
         word(length:length) = letter
         next = next + 1
         letter = string(next:next)
         code = ichar(letter)
         length = length + 1
      end do
      if (code .eq. atsign) then
         word(1:1) = letter
         length = 1
      end if
   30 continue
      return
      end
