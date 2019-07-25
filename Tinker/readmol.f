c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2012  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine readmol  --  read in a MDL MOL format file  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "readmol" gets a set of MDL MOL coordinates from
c     an external disk file
c
c
      subroutine readmol (imdl)
      use sizes
      use atomid
      use atoms
      use couple
      use files
      use iounit
      use ptable
      use titles
      implicit none
      integer i,j,ia,ib,imdl
      integer nbond
      integer trimtext
      logical exist,opened
      character*240 mdlfile
c
c
c     open the input file if it has not already been done
c
      inquire (unit=imdl,opened=opened)
      if (.not. opened) then
         mdlfile = filename(1:leng)//'.mol'
         call version (mdlfile,'old')
         inquire (file=mdlfile,exist=exist)
         if (exist) then
            open (unit=imdl,file=mdlfile,status='old')
            rewind (unit=imdl)
         else
            write (iout,10)
   10       format (/,' READMOL  --  Unable to Find the Specified',
     &                 ' MDL MOL File')
            call fatal
         end if
      end if
c
c     zero out the total number of atoms and of bonds
c
      n = 0
      nbond = 0
c
c     get title line and get the number of atoms and bonds
c
      read (imdl,20)  title
   20 format (a240)
      ltitle = trimtext (title)
      read (imdl,30)
   30 format (/)
      read (imdl,40)  n,nbond
   40 format (2i3)
c
c     check for too few or too many total atoms in the file
c
      if (n .le. 0) then
         write (iout,50)
   50    format (/,' READMOL  --  The Coordinate File Does Not',
     &              ' Contain Any Atoms')
         call fatal
      else if (n .gt. maxatm) then
         write (iout,60)  maxatm
   60    format (/,' READMOL  --  The Maximum of',i9,' Atoms',
     &              ' has been Exceeded')
         call fatal
      end if
c
c     read the atomic coordinates and atomic symbol
c
      do i = 1, n
         read (imdl,70)  x(i),y(i),z(i),name(i)
   70    format (3f10.4,1x,a3)
         n12(i) = 0
      end do
c
c     read the bond list to get attached atom lists
c
      do i = 1, nbond
         read (imdl,80)  ia,ib
   80    format (2i3)
         n12(ia) = n12(ia) + 1
         i12(n12(ia),ia) = ib
         n12(ib) = n12(ib) + 1
         i12(n12(ib),ib) = ia
      end do
c
c     assign atom types from atomic number and connectivity
c
      do i = 1, n
         type(i) = 0
         do j = 1, maxele
            if (name(i) .eq. elemnt(j)) then
               type(i) = 10*j + n12(i)
               goto 90
            end if
         end do
   90    continue
      end do
c
c     for each atom, sort its list of attached atoms
c
      do i = 1, n
         call sort (n12(i),i12(1,i))
      end do
      if (.not. opened)  close (unit=imdl)
      return
      end
