c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine makeint  --  convert Cartesian to internal  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "makeint" converts Cartesian to internal coordinates where
c     selection of internal coordinates is controlled by "mode"
c
c        mode = 0     automatic internal coordinates
c        mode = 1     manual selection of coordinates
c        mode = 2     use existing structure as a template
c        mode = 3     use dihedral angles in all cases
c
c
      subroutine makeint (mode)
      use sizes
      use atoms
      use couple
      use iounit
      use math
      use zclose
      use zcoord
      implicit none
      integer i,j
      integer i1,i2,i3,i4,i5
      integer adjacent,trial
      integer mode,next
      integer, allocatable :: iz0(:)
      integer, allocatable :: iz1(:)
      real*8 geometry,sign
      logical more
      character*1 answer
      character*1 default
      character*8 phrase
      character*240 record
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (iz0(0:n))
      allocate (iz1(n))
c
c     zero out local values used for the defining atoms
c
      i1 = 0
      i2 = 0
      i3 = 0
      i4 = 0
      i5 = 0
      iz0(0) = 0
      do i = 1, n
         iz0(i) = 0
         iz1(i) = 0
      end do
c
c     zero out the coordinates, defining atoms and closures
c
      do i = 1, n
         zbond(i) = 0.0d0
         zang(i) = 0.0d0
         ztors(i) = 0.0d0
      end do
      if (mode .ne. 2) then
         do i = 1, n
            do j = 1, 4
               iz(j,i) = 0
            end do
         end do
         nadd = 0
         ndel = 0
      end if
c
c     first, decide which of the atoms to define next
c
      do i = 1, n
         if (mode .eq. 1) then
            trial = i1 + 1
   10       continue
            write (iout,20)  trial
   20       format (/,' Atom Number to be Defined [',i5,'] :  ',$)
            read (input,30,err=10)  i1
   30       format (i10)
            if (i1 .eq. 0)  i1 = trial
            if (iz0(i1) .ne. 0) then
               write (iout,40)
   40          format (/,' Already Defined that Atom; Choose Another')
               if (i1 .eq. trial)  trial = trial + 1
               goto 10
            end if
         else
            i1 = i
         end if
c
c     define the bond length for the current atom
c
         if (i .ge. 2) then
            if (mode .eq. 2) then
               i2 = iz(1,i1)
            else
               i2 = adjacent (i1,0,mode,more,iz0,iz1)
               if (i2 .eq. 0) then
                  write (iout,50)  i1
   50             format (/,' MAKEINT  --  Connectivity Error',
     &                       ' in defining Atom',i6)
                  call fatal
               end if
            end if
            zbond(i1) = geometry (i1,i2,0,0)
         end if
c
c     define the bond angle for the current atom
c
         if (i .ge. 3) then
            if (mode .eq. 2) then
               i3 = iz(2,i1)
            else
               i3 = adjacent (i2,i1,mode,more,iz0,iz1)
               if (i3 .eq. 0) then
                  write (iout,60)  i1
   60             format (/,' MAKEINT  --  Connectivity Error',
     &                       ' in defining Atom',i6)
                  call fatal
               end if
            end if
            zang(i1) = geometry (i1,i2,i3,0)
         end if
c
c     decide whether to use a dihedral or second bond angle;
c     then find the value of the angle
c
         if (i .ge. 4) then
            if (mode .eq. 3) then
               answer = 'D'
            else if (mode .eq. 2) then
               if (iz(4,i1) .eq. 0) then
                  answer = 'D'
               else
                  answer = 'B'
               end if
            else if (mode .eq. 1) then
               if (more) then
                  phrase = 'D or [B]'
                  default = 'B'
               else
                  phrase = '[D] or B'
                  default = 'D'
               end if
               write (iout,70)  phrase
   70          format (/,' Specify with Dihedral Angle or Second',
     &                    ' Bond Angle (',a8,') :  ',$)
               read (input,80)  record
   80          format (a240)
               next = 1
               call gettext (record,answer,next)
               call upcase (answer)
               if (answer.ne.'B' .and. answer.ne.'D')  answer = default
            else if (mode .eq. 0) then
               if (more) then
                  answer = 'B'
               else
                  answer = 'D'
               end if
            end if
            if (answer .eq. 'B') then
               if (mode .eq. 2) then
                  i4 = iz(3,i1)
               else
                  i4 = adjacent (i2,i3,mode,more,iz0,iz1)
                  if (i4 .eq. 0) then
                     write (iout,90)  i1
   90                format (/,' MAKEINT  --  Connectivity Error',
     &                          ' in defining Atom',i6)
                     call fatal
                  end if
               end if
               ztors(i1) = geometry (i1,i2,i4,0)
               i5 = 1
               sign = geometry (i1,i2,i3,i4)
               if (sign .gt. 0.0d0)  i5 = -1
            else if (answer .eq. 'D') then
               if (mode .eq. 2) then
                  i4 = iz(3,i1)
               else
                  i4 = adjacent (i3,i2,mode,more,iz0,iz1)
                  if (i4 .eq. 0) then
                     write (iout,100)  i1
  100                format (/,' MAKEINT  --  Connectivity Error',
     &                          ' in defining Atom',i6)
                     call fatal
                  end if
               end if
               i5 = 0
               ztors(i1) = geometry (i1,i2,i3,i4)
            end if
         end if
c
c     transfer defining atoms to permanent array;
c     mark the current atom as finished
c
         iz(1,i1) = iz0(i2)
         iz(2,i1) = iz0(i3)
         iz(3,i1) = iz0(i4)
         iz(4,i1) = i5
         iz0(i1) = i
         iz1(i1) = i2
      end do
c
c     add any bonds needed to make ring closures
c
      nadd = 0
      do i = 1, n
         do j = 1, n12(i)
            if (iz0(i) .lt. iz0(i12(j,i)) .and.
     &          iz1(i12(j,i)) .ne. i) then
               nadd = nadd + 1
               iadd(1,nadd) = iz0(i)
               iadd(2,nadd) = iz0(i12(j,i))
            end if
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (iz0)
      deallocate (iz1)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  function adjacent  --  atom adjacent to specified atom  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "adjacent" finds an atom connected to atom "i1" other than
c     atom "i2"; if no such atom exists, then the closest atom
c     in space is returned
c
c     variables and parameters :
c
c     mode   whether "makeint" is in manual mode, automatic, etc.
c     more   returned true if there is more than one previously
c              defined atom other than "i2" which is directly
c              connected (adjacent) to atom "i1"
c     iz0    line number of the Z-matrix on which an atom is
c              defined, 0 if not yet defined
c     iz1    line number of the Z-matrix on which the atom used
c              defining the bond length to a given atom is defined
c
c
      function adjacent (i1,i2,mode,more,iz0,iz1)
      use sizes
      use atoms
      use couple
      use inform
      use iounit
      use zclose
      implicit none
      integer i,j,k,i1,i2
      integer nc,adjacent,mode
      integer ic(maxval)
      integer iz0(0:*)
      integer iz1(*)
      real*8 dist2,short
      logical more
c
c
c     get a list of eligible atoms bonded to the atom of interest
c
      nc = 0
      more = .false.
      do j = 1, n12(i1)
         i = i12(j,i1)
         if (iz0(i).ne.0 .and. i.ne.i2) then
            if (i2 .eq. 0) then
               nc = nc + 1
               ic(nc) = i
            else
               if (iz1(i).eq.i1 .or. iz1(i1).eq.i) then
                  nc = nc + 1
                  ic(nc) = i
               end if
            end if
         end if
      end do
      if (nc .gt. 1)  more = .true.
c
c     if no bonded atom is eligible, use the nearest neighbor
c
      if (nc .eq. 0) then
         adjacent = 0
         if (mode .eq. 1) then
            write (iout,10)  i1
   10       format (/,' ADJACENT  --  Atom',i6,' not Attached',
     &                 ' to any Prior Atom')
         else
            short = 100000000.0d0
            do i = 1, n
               if (iz0(i).ne.0 .and. i.ne.i1 .and. i.ne.i2) then
                  dist2 = (x(i)-x(i1))**2 + (y(i)-y(i1))**2
     &                           + (z(i)-z(i1))**2
                  if (dist2 .lt. short) then
                     short = dist2
                     adjacent = i
                  end if
               end if
            end do
            if (i2 .eq. 0) then
               ndel = ndel + 1
               idel(1,ndel) = adjacent
               idel(2,ndel) = i1
               if (debug) then
                  write (iout,20)  i1
   20             format (/,' ADJACENT  --  Atom',i6,' not Attached',
     &                       ' to any Prior Atom')
               end if
            end if
         end if
c
c     for automatic mode, always use the first eligible bonded atom
c
      else if (mode .eq. 0) then
         adjacent = ic(1)
c
c     for torsion mode, use an adjacent atom bonded to undefined atoms
c
      else if (mode .eq. 3) then
         adjacent = ic(1)
         do k = 1, nc
            do j = 1, n12(ic(k))
               i = i12(j,ic(k))
               if (iz0(i).ne.0 .and. i.ne.i1) then
                  adjacent = ic(k)
                  goto 30
               end if
            end do
         end do
   30    continue
c
c     if only one directly bonded atom is eligible, then use it
c
      else if (nc .eq. 1) then
         adjacent = ic(1)
         if (mode.eq.1 .or. debug) then
            write (iout,40)  ic(1)
   40       format (/,' ADJACENT  --  Atom',i6,' is the only',
     &                 ' Connected Atom')
         end if
c
c     ask the user which eligible bonded atom to use as adjacent
c
      else
   50    continue
         if (nc .eq. 2) then
            write (iout,60)  (ic(j),j=1,nc)
   60       format (' Choose a Connected Atom (',2i6,') :  ',$)
         else if (nc .eq. 3) then
            write (iout,70)  (ic(j),j=1,nc)
   70       format (' Choose a Connected Atom (',3i6,') :  ',$)
         else if (nc .eq. 4) then
            write (iout,80)  (ic(j),j=1,nc)
   80       format (' Choose a Connected Atom (',4i6,') :  ',$)
         else if (nc .eq. 5) then
            write (iout,90)  (ic(j),j=1,nc)
   90       format (' Choose a Connected Atom (',5i6,') :  ',$)
         else if (nc .eq. 6) then
            write (iout,100)  (ic(j),j=1,nc)
  100       format (' Choose a Connected Atom (',6i6,') :  ',$)
         else if (nc .eq. 7) then
            write (iout,110)  (ic(j),j=1,nc)
  110       format (' Choose a Connected Atom (',7i6,') :  ',$)
         else if (nc .eq. 8) then
            write (iout,120)  (ic(j),j=1,nc)
  120       format (' Choose a Connected Atom (',8i6,') :  ',$)
         end if
         read (input,130,err=50)  adjacent
  130    format (i10)
         if (adjacent .eq. 0) then
            adjacent = ic(1)
         else
            do j = 1, nc
               if (ic(j) .eq. adjacent)  goto 140
            end do
            goto 50
  140       continue
         end if
      end if
      return
      end
