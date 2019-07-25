c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine initrot  --  set bonds for dihedral rotation  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "initrot" sets the torsional angles which are to be rotated
c     in subsequent computation, by default automatically selects
c     all rotatable single bonds; optionally makes atoms inactive
c     when they are not moved by any torsional rotation
c
c     note that internal coordinates must already be setup
c
c
      subroutine initrot
      use sizes
      use atoms
      use couple
      use group
      use inform
      use iounit
      use math
      use omega
      use potent
      use restrn
      use rotbnd
      use usage
      use zcoord
      implicit none
      integer i,j,j1,j2
      integer mode,iring
      integer bond1,bond2
      integer attach1,attach2
      integer nlist,nfixed
      integer, allocatable :: list(:)
      integer, allocatable :: ifixed(:,:)
      logical exist,query
      logical rotate,rotcheck
      logical use_partial
      character*240 record
      character*240 string
c
c
c     initialize the number of rotatable torsional angles
c
      nomega = 0
c
c     use partial structure, mark inactive any atoms that do not move;
c     faster for limited torsions, only use with pairwise potentials
c
      use_partial = .true.
      if (use_polar)  use_partial = .false.
c
c     use shortest rotlist if there is no absolute coordinate frame
c
      use_short = .true.
      if (use_group)  use_short = .false.
      if (npfix .ne. 0)  use_short = .false.
c
c     choose automatic or manual selection of torsional angles
c
      mode = 0
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  mode
         query = .false.
      end if
   10 continue
      if (query) then
         write (iout,20)
   20    format (/,' Selection of Torsional Angles for Rotation :',
     &           //,'    0  - Automatic Selection of Torsional Angles',
     &            /,'    1  - Manual Selection of Angles to Rotate',
     &            /,'    2  - Manual Selection of Angles to Freeze',
     &           //,' Enter the Method of Choice [0] :  ',$)
         read (input,30)  mode
   30    format (i10)
      end if
      if (mode.ne.1 .and. mode.ne.2)  mode = 0
c
c     perform dynamic allocation of some global arrays
c
      if (.not. allocated(iomega))  allocate (iomega(2,n))
      if (.not. allocated(zline))  allocate (zline(n))
      if (.not. allocated(dihed))  allocate (dihed(n))
c
c     manual selection of the torsional angles to be rotated
c
      if (mode .eq. 1) then
         do while (.true.)
            nomega = nomega + 1
            j1 = 0
            j2 = 0
            write (iout,40)  nomega
   40       format (/,' Enter Atoms in Rotatable Bond',i5,' :  ',$)
            read (input,50)  record
   50       format (a240)
            read (record,*,err=80,end=80)  j1,j2
            if (j1.eq.0 .and. j2.eq.0)  goto 80
            do i = 4, n
               if (iz(4,i) .eq. 0) then
                  bond1 = iz(1,i)
                  bond2 = iz(2,i)
                  attach1 = n12(bond1)
                  attach2 = n12(bond2)
                  if (attach1.gt.1 .and. attach2.gt.1) then
                     if ((bond1.eq.j1 .and. bond2.eq.j2) .or.
     &                   (bond1.eq.j2 .and. bond2.eq.j1)) then
                        if (rotcheck(bond1,bond2)) then
                           iomega(1,nomega) = bond1
                           iomega(2,nomega) = bond2
                           dihed(nomega) = ztors(i) / radian
                           zline(nomega) = i
                           goto 70
                        end if
                     end if
                  end if
               end if
            end do
            nomega = nomega - 1
            write (iout,60)  j1,j2
   60       format (/,' INITROT  --  Bond between Atoms',2i6,
     &                 ' is not Rotatable')
   70       continue
         end do
   80    continue
         nomega = nomega - 1
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (ifixed(2,n))
c
c     manual selection of the torsional angles to be frozen
c
      nfixed = 0
      if (mode .eq. 2) then
         do i = 1, n
            ifixed(1,i) = 0
            ifixed(2,i) = 0
            write (iout,90)  i
   90       format (/,' Enter Atoms in Frozen Bond',i5,' :  ',$)
            read (input,100)  record
  100       format (a240)
            read (record,*,err=110,end=110)  ifixed(1,i),ifixed(2,i)
            if (ifixed(1,i).eq.0 .or. ifixed(2,i).eq.0)  goto 110
            nfixed = nfixed + 1
         end do
  110    continue
      end if
c
c     perform the automatic selection of torsional angles to rotate
c
      if (mode.eq.0 .or. mode.eq.2) then
         do i = 4, n
            if (iz(4,i) .eq. 0) then
               rotate = .true.
               bond1 = iz(1,i)
               bond2 = iz(2,i)
c
c     do not rotate a bond if either bonded atom is univalent
c
               attach1 = n12(bond1)
               attach2 = n12(bond2)
               if (attach1.le.1 .or. attach2.le.1)  rotate = .false.
c
c     do not rotate a bond contained within a small ring
c
               iring = 0
               call chkring (iring,bond1,bond2,0,0)
               if (iring .ne. 0)  rotate = .false.
c
c     do not rotate bonds explicitly frozen by the user
c
               if (mode.eq.2 .and. rotate) then
                  do j = 1, nfixed
                     j1 = ifixed(1,j)
                     j2 = ifixed(2,j)
                     if ((bond1.eq.j1 .and. bond2.eq.j2) .or.
     &                   (bond1.eq.j2 .and. bond2.eq.j1)) then
                        rotate = .false.
                        goto 120
                     end if
                  end do
               end if
  120          continue
c
c     do not rotate bonds with inactive atoms on both sides
c
               if (rotate) then
                  if (.not. rotcheck(bond1,bond2))  rotate = .false.
               end if
c
c     check for possible duplication of rotatable bonds
c
               if (rotate) then
                  do j = 1, nomega
                     j1 = iomega(1,j)
                     j2 = iomega(2,j)
                     if ((bond1.eq.j1 .and. bond2.eq.j2) .or.
     &                   (bond1.eq.j2 .and. bond2.eq.j1)) then
                        write (iout,130)  bond1,bond2
  130                   format (/,' INITROT  --  Rotation about',2i6,
     &                             ' occurs more than once in Z-matrix')
                        call fatal
                     end if
                  end do
                  nomega = nomega + 1
                  iomega(1,nomega) = bond1
                  iomega(2,nomega) = bond2
                  dihed(nomega) = ztors(i) / radian
                  zline(nomega) = i
               end if
            end if
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (ifixed)
c
c     perform dynamic allocation of some local arrays
c
      allocate (list(n))
c
c     make inactive the atoms not rotatable via any torsion
c
      if (use_partial .and. nuse.eq.n) then
         do i = 1, n
            use(i) = .false.
         end do
         do i = 1, nomega
            bond1 = iomega(1,i)
            bond2 = iomega(2,i)
            call rotlist (bond1,bond2)
            do j = 1, nrot
               use(rot(j)) = .true.
            end do
         end do
         nuse = 0
         do i = 1, n
            if (use(i))  nuse = nuse + 1
         end do
         if (debug .and. nuse.gt.0 .and. nuse.lt.n) then
            nlist = 0
            do i = 1, n
               if (use(i)) then
                  nlist = nlist + 1
                  list(nlist) = i
               end if
            end do
            write (iout,140)
  140       format (/,' List of Active Atoms for Torsional',
     &                    ' Calculations :',/)
            write (iout,150)  (list(i),i=1,nlist)
  150       format (3x,10i7)
         end if
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (list)
c
c     write out the number of rotatable torsions to be used
c
      if (nomega .eq. 0) then
         write (iout,160)
  160    format (/,' INITROT  --  No Torsions for Subsequent',
     &              ' Computation')
         call fatal
      end if
      write (iout,170)  nomega
  170 format (/,' Number of Torsions Used in Derivative',
     &           ' Computation :',i6)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  function rotcheck  --  check for fixed atoms across bond  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "rotcheck" tests a specified candidate rotatable bond for
c     the disallowed case where inactive atoms are found on both
c     sides of the candidate bond
c
c
      function rotcheck (base,partner)
      use sizes
      use atoms
      use rotbnd
      use usage
      implicit none
      integer i,base,partner
      logical rotcheck,value
      logical, allocatable :: list(:)
c
c
c     initialize status and find atoms on short side of the bond
c
      value = .true.
      call rotlist (base,partner)
c
c     rotation is allowed if all atoms on one side are active
c
      do i = 1, nrot
         if (.not. use(rot(i))) then
            value = .false.
            goto 10
         end if
      end do
   10 continue
c
c     if short side had inactive atoms, check the other side
c
      if (.not. value) then
         allocate (list(n))
         do i = 1, n
            list(i) = .true.
         end do
         do i = 1, nrot
            list(rot(i)) = .false.
         end do
         do i = 1, n
            if (list(i) .and. .not.use(i))  goto 20
         end do
         value = .true.
   20    continue
         deallocate (list)
      end if
c
c     set the final return value of the function
c
      rotcheck = value
      return
      end
