c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  program intedit  --  edit and display Z-matrix file  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "intedit" allows the user to extract information from
c     or alter the values within an internal coordinates file
c
c
      program intedit
      use sizes
      use atomid
      use atoms
      use files
      use iounit
      use katoms
      use zcoord
      implicit none
      integer i,j,k,l,m
      integer izmt,space
      integer freeunit
      integer trimtext
      integer numrow,numcol
      integer next,number(4)
      real*8 value,geometry
      logical changed,error
      character*4 word
      character*240 zmtfile
      character*240 record
c
c
c     read coordinate file and force field definition
c
      call initial
      call getint
      call field
c
c     print out the instructions for the program
c
      next = 1
      changed = .false.
      error = .false.
   10 continue
      call zhelp
c
c     start of main loop, examine or change Z-matrix elements
c
   20 continue
      m = 0
      write (iout,30)
   30 format (/,' INTEDIT>  ',$)
      read (input,40)  record
   40 format (a240)
c
c     interpret any user entered text command
c
      space = 1
      call getword (record,word,space)
      call upcase (word)
      if (word .eq. 'EXIT') then
         if (changed) then
            izmt = freeunit ()
            zmtfile = filename(1:leng)//'.int'
            call version (zmtfile,'new')
            open (unit=izmt,file=zmtfile,status='new')
            call prtint (izmt)
            close (unit=izmt)
            write (iout,50)  zmtfile(1:trimtext(zmtfile))
   50       format (/,' Z-Matrix Internal Coordinates written to :  ',a)
         else
            write (iout,60)
   60       format (/,' The Z-Matrix was not Changed;',
     &                 ' No File was Written')
         end if
         goto 410
      else if (word .eq. 'QUIT') then
         goto 410
      else if (word .eq. 'SHOW') then
         write (iout,70)
   70    format ()
         call prtint (iout)
c
c     get the number of atoms entered by the user
c
      else
         do i = 1, 4
            number(i) = 0
         end do
         read (record,*,err=10,end=80)  (number(i),i=1,4)
   80    continue
         do i = 1, 4
            if (number(i) .ne. 0)  m = i
            if (number(i) .gt. n) then
               write (iout,90)  n
   90          format (/,' Warning; Only',i6,' Atoms are Present',
     &                     ' in the Z-matrix')
               goto 20
            end if
         end do
         if (m .eq. 0) then
            m = 1
            number(1) = next
         end if
      end if
c
c     get information about a single specified atom
c
      if (m .eq. 1) then
         i = number(1)
         write (iout,100)  i
  100    format (/,' Atom Number :',i8)
         write (iout,110)  name(i)
  110    format (' Atom Name :',6x,a4)
         write (iout,120)  describe(type(i))
  120    format (' Atom Type :',5x,a20)
         write (iout,130)  type(i)
  130    format (' Type Number :',i8)
         if (i .eq. 1) then
            write (iout,140)
  140       format (/,' Atom 1 is at the Coordinate System Origin')
         else
            write (iout,150)
  150       format (/,' Internal Coordinate Structural Definition :',/)
            write (iout,160)  iz(1,i),-i,zbond(i)
  160       format (1x,2i6,17x,'Distance Value :',f14.4)
            if (i .gt. 2) then
               write (iout,170)  iz(2,i),-iz(1,i),-i,zang(i)
  170          format (1x,3i6,11x,'Bond Angle Value :',f12.4)
               if (i .gt. 3) then
                  if (iz(4,i) .eq. 0) then
                     write (iout,180) iz(3,i),-iz(2,i),
     &                                -iz(1,i),-i,ztors(i)
  180                format (1x,4i6,5x,'Dihedral Angle :',f14.4)
                  else
                     write (iout,190)  iz(3,i),-iz(1,i),-i,ztors(i)
  190                format (1x,3i6,11x,'Bond Angle Value :',f12.4)
                     write (iout,200)  iz(4,i)
  200                format (30x,'Chirality Flag :',6x,i6)
                  end if
               end if
            end if
         end if
         next = i + 1
         if (next .gt. n)  next = 1
c
c     chirality change for an atom was requested
c
      else if (m.eq.2 .and. number(2).lt.0) then
         do i = 1, n
            if (iz(4,i).ne.0 .and. iz(1,i).eq.number(1)) then
               changed = .true.
               write (iout,210)  i
  210          format (/,' Inverting Chirality of Atom : ',i6)
               iz(4,i) = -iz(4,i)
            end if
         end do
         next = number(1)
         call makexyz
c
c     information about a specified bond or distance
c
      else if (m .eq. 2) then
         i = max(number(1),number(2))
         j = min(number(1),number(2))
         if (min(i,j).le.0 .and. max(i,j).gt.n) then
            write (iout,220)
  220       format (/,' Invalid Atom Number')
            error = .true.
         else
            if (j .ne. iz(1,i)) then
               value = geometry (i,j,0,0)
               write (iout,230)  value
  230          format (/,' The Current Distance is : ',f9.4)
               write (iout,240)
  240          format (' That Bond is not in the Z-matrix')
            else
               write (iout,250)  zbond(i)
  250          format (/,' The Current Distance is : ',f9.4)
               call zvalue ('Bond Length',zbond(i),changed)
               next = i
            end if
         end if
c
c     an atom type change was requested
c
      else if (m.eq.3 .and. number(2).lt.0) then
         if (number(3).gt.0 .and. number(3).le.maxtyp) then
            changed = .true.
            write (iout,260)  describe(type(number(1)))
  260       format (/,' Old Atom Type is :  ',a20)
            type(number(1)) = number(3)
            write (iout,270)  describe(type(number(1)))
  270       format (' New Atom Type is :  ',a20)
         else
            write (iout,280)
  280       format (/,' Invalid Atom Type; Valid Types are :',/)
            numrow = (maxtyp+2) / 3
            numcol = 2
            do i = 1, numrow
               if (i .gt. numrow-2+mod(maxtyp-1,3))  numcol = 1
               write (iout,290)  (j*numrow+i,describe(j*numrow+i),
     &                                      j=0,numcol)
  290          format (1x,3(i3,1x,a20,2x))
            end do
         end if
c
c     information about a specified bond angle
c
      else if (m .eq. 3) then
         i = max(number(1),number(3))
         j = number(2)
         k = min(number(1),number(3))
         if (min(i,j,k).le.0 .or. max(i,j,k).gt.n) then
            write (iout,300)
  300       format (/,' Invalid Atom Number')
            error = .true.
         else
            if (iz(1,i) .ne. j) then
               value = geometry (i,j,k,0)
               write (iout,310)  value
  310          format (/,' The Bond Angle Value is :  ',f9.4)
               write (iout,320)
  320          format (' That Bond Angle is not in the Z-matrix')
            else if (iz(2,i) .eq. k) then
               write (iout,330)  zang(i)
  330          format (/,' The Bond Angle Value is :  ',f9.4)
               call zvalue ('Bond Angle',zang(i),changed)
               next = i
            else if (iz(3,i).eq.k .and. iz(4,i).ne.0) then
               write (iout,340)  ztors(i)
  340          format (/,' The Bond Angle Value is :  ',f9.4)
               call zvalue ('Bond Angle',ztors(i),changed)
               next = i
            else
               value = geometry (i,j,k,0)
               write (iout,350)  value
  350          format (/,' The Bond Angle Value is :  ',f9.4)
               write (iout,360)
  360          format (' That Bond Angle is not in the Z-matrix')
            end if
         end if
c
c    information about a specified dihedral angle
c
      else if (m .eq. 4) then
         if (number(1) .gt. number(4)) then
            i = number(1)
            j = number(2)
            k = number(3)
            l = number(4)
         else
            i = number(4)
            j = number(3)
            k = number(2)
            l = number(1)
         end if
         if (min(i,j,k,l).le.0 .or. max(i,j,k,l).gt.n) then
            write (iout,370)
  370       format (/,' Invalid Atom Number')
            error = .true.
         else
            if (iz(1,i).ne.j .or. iz(2,i).ne.k .or.
     &          iz(3,i).ne.l .or. iz(4,i).ne.0) then
               value = geometry (i,j,k,l)
               write (iout,380)  value
  380          format (/,' The Dihedral Angle Value is :  ',f9.4)
               write (iout,390)
  390          format (' That Dihedral Angle is not in the Z-matrix')
            else
               write (iout,400)  ztors(i)
  400          format (/,' The Dihedral Angle Value is :  ',f9.4)
               call zvalue ('Dihedral Angle',ztors(i),changed)
               next = i
            end if
         end if
      end if
c
c     print instructions for the program if needed
c
      if (error) then
         error = .false.
         call zhelp
      end if
      goto 20
c
c     perform any final tasks before program exit
c
  410 continue
      call final
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine zhelp  --  print Z-matrix editing instructions  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "zhelp" prints the general information and instructions
c     for the Z-matrix editing program
c
c
      subroutine zhelp
      use iounit
      implicit none
c
c
c     print the help and information message for Z-matrix editing
c
      write (iout,10)
   10 format (/,' If a single atom number is entered, the',
     &           ' current definition of',
     &        /,' the atom will be displayed.',
     &        //,' If two atom numbers are entered, the output',
     &           ' gives the distance',
     &        /,' between the atoms, and asks for a new bond',
     &           ' length if applicable;',
     &        /,' Entry of three atoms shows the angle, and',
     &           ' entry of four atoms',
     &        /,' will display the corresponding dihedral angle.',
     &        //,' To change the chirality at an atom, enter',
     &           ' its number and -1.',
     &        /,' To change the type of an atom, enter its',
     &           ' number, -1, and the',
     &        /,' new atom type number.')
      write (iout,20)
   20 format (/,' A carriage return at the prompt will display',
     &           ' the atom last',
     &        /,' changed or the next atom after the one just',
     &           ' examined.',
     &        //,' Typing SHOW will display the contents of the',
     &           ' current Z-matrix.',
     &        //,' Entering EXIT writes a new file then stops,',
     &           ' while QUIT aborts.')
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine zvalue  --  gets user input Z-matrix value  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "zvalue" gets user supplied values for selected coordinates
c     as needed by the internal coordinate editing program
c
c
      subroutine zvalue (text,x,changed)
      use iounit
      implicit none
      integer length
      integer trimtext
      real*8 x,xnew
      logical changed
      character*240 record
      character*(*) text
c
c
c     ask the user for the new internal coordinate value
c
      xnew = x
      write (iout,10)  text
   10 format (/,' Enter the New ',a,' :  ',$)
      read (input,20)  record
   20 format (a240)
      length = trimtext (record)
      if (length .ne. 0) then
         read (record,*,end=30,err=30)  xnew
   30    continue
      end if
c
c     return with the altered value and recompute coordinates
c
      if (xnew .ne. x) then
         changed = .true.
         x = xnew
         call makexyz
      end if
      return
      end
