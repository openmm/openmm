c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine orbital  --  setup for pisystem calculation  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "orbital" finds and organizes lists of atoms in a pisystem,
c     bonds connecting pisystem atoms and torsions whose central
c     atoms are both pisystem atoms
c
c
      subroutine orbital
      use sizes
      use atoms
      use bndstr
      use couple
      use iounit
      use keys
      use piorbs
      use potent
      use tors
      implicit none
      integer i,j,k,m,ii
      integer mi,mj,mk
      integer iorb,jorb,korb
      integer nlist,next
      integer, allocatable :: list(:)
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(iorbit))  deallocate (iorbit)
      if (allocated(iconj))  deallocate (iconj)
      if (allocated(kconj))  deallocate (kconj)
      if (allocated(piperp))  deallocate (piperp)
      if (allocated(ibpi))  deallocate (ibpi)
      if (allocated(itpi))  deallocate (itpi)
      if (allocated(pbpl))  deallocate (pbpl)
      if (allocated(pnpl))  deallocate (pnpl)
      if (allocated(listpi))  deallocate (listpi)
      allocate (iorbit(n))
      allocate (iconj(2,n))
      allocate (kconj(n))
      allocate (piperp(3,n))
      allocate (ibpi(3,nbond))
      allocate (itpi(2,ntors))
      allocate (pbpl(nbond))
      allocate (pnpl(nbond))
      allocate (listpi(n))
c
c     perform dynamic allocation of some local arrays
c
      allocate (list(n))
c
c     set the default values for the pisystem variables
c
      nlist = 0
      do i = 1, n
         list(i) = 0
      end do
      reorbit = 1
c
c     check the keywords for any lists of pisystem atoms
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:9) .eq. 'PISYSTEM ') then
            string = record(next:240)
            read (string,*,err=10,end=10)  (list(k),k=nlist+1,n)
   10       continue
            do while (list(nlist+1) .ne. 0)
               nlist = nlist + 1
               list(nlist) = max(-n,min(n,list(nlist)))
            end do
         end if
      end do
c
c     quit if no pisystem was found for consideration
c
      if (list(1) .eq. 0) then
         use_orbit = .false.
         return
      else
         use_orbit = .true.
      end if
c
c     organize and make lists of the pisystem atoms
c
      do i = 1, n
         listpi(i) = .false.
      end do
      i = 1
      do while (list(i) .ne. 0)
         if (list(i) .gt. 0) then
            listpi(list(i)) = .true.
            i = i + 1
         else
            do j = -list(i), list(i+1)
               listpi(j) = .true.
            end do
            i = i + 2
         end if
      end do
c
c     set number of orbitals and an initial orbital list
c
      norbit = 0
      nconj = 0
      do i = 1, n
         list(i) = 0
         if (listpi(i)) then
            norbit = norbit + 1
            iorbit(norbit) = i
         end if
      end do
c
c     assign each orbital to its respective pisystem
c
      do i = 1, norbit
         iorb = iorbit(i)
         if (list(iorb) .eq. 0) then
            nconj = nconj + 1
            list(iorb) = nconj
         end if
         mi = list(iorb)
         do ii = 1, n12(iorb)
            j = i12(ii,iorb)
            if (listpi(j)) then
               mj = list(j)
               if (mj .eq. 0) then
                  list(j) = mi
               else if (mi .lt. mj) then
                  nconj = nconj - 1
                  do k = 1, norbit
                     korb = iorbit(k)
                     mk = list(korb)
                     if (mk .eq. mj) then
                        list(korb) = mi
                     else if (mk .gt. mj) then
                        list(korb) = mk - 1
                     end if
                  end do
               else if (mi .gt. mj) then
                  nconj = nconj - 1
                  do k = 1, norbit
                     korb = iorbit(k)
                     mk = list(korb)
                     if (mk .eq. mi) then
                        list(korb) = mj
                     else if (mk .gt. mi) then
                        list(korb) = mk - 1
                     end if
                  end do
                  mi = mj
               end if
            end if
         end do
      end do
c
c     pack atoms of each pisystem into a contiguous indexed list
c
      call sort3 (n,list,kconj)
      k = n - norbit
      do i = 1, norbit
         k = k + 1
         list(i) = list(k)
         kconj(i) = kconj(k)
      end do
c
c     find the first and last piatom in each pisystem
c
      k = 1
      iconj(1,1) = 1
      do i = 2, norbit
         j = list(i)
         if (j .ne. k) then
            iconj(2,k) = i - 1
            k = j
            iconj(1,k) = i
         end if
      end do
      iconj(2,nconj) = norbit
c
c     perform deallocation of some local arrays
c
      deallocate (list)
c
c     sort atoms in each pisystem, copy list to orbital sites
c
      do i = 1, nconj
         k = iconj(2,i) - iconj(1,i) + 1
         call sort (k,kconj(iconj(1,i)))
      end do
      do i = 1, norbit
         iorbit(i) = kconj(i)
      end do
c
c     find atoms defining a plane perpendicular to each orbital
c
      call piplane
c
c     find and store all of the pisystem bonds
c
      nbpi = 0
      do ii = 1, nconj
         do i = iconj(1,ii), iconj(2,ii)-1
            iorb = kconj(i)
            do j = i+1, iconj(2,ii)
               jorb = kconj(j)
               do k = 1, n12(iorb)
                  if (i12(k,iorb) .eq. jorb) then
                     nbpi = nbpi + 1
                     do m = 1, nbond
                        if (iorb.eq.ibnd(1,m) .and.
     &                      jorb.eq.ibnd(2,m)) then
                           ibpi(1,nbpi) = m
                           ibpi(2,nbpi) = i
                           ibpi(3,nbpi) = j
                           goto 20
                        end if
                     end do
   20                continue
                  end if
               end do
            end do
         end do
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine piplane  --  plane perpendicular to orbital  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "piplane" selects the three atoms which specify the plane
c     perpendicular to each p-orbital; the current version will
c     fail in certain situations, including ketenes, allenes,
c     and isolated or adjacent triple bonds
c
c
      subroutine piplane
      use sizes
      use atomid
      use atoms
      use couple
      use iounit
      use piorbs
      implicit none
      integer i,j,iorb
      integer atmnum,trial
      integer alpha,beta,gamma
      integer attach
      logical done
c
c
c     for each pisystem atom, find a set of atoms which define
c     the p-orbital's plane based on piatom's atomic number and
c     the number and type of attached atoms
c
      do iorb = 1, norbit
         i = iorbit(iorb)
         attach = n12(i)
         atmnum = atomic(i)
         done = .false.
c
c     most common case of an atom bonded to three atoms
c
         if (attach .eq. 3) then
            piperp(1,i) = i12(1,i)
            piperp(2,i) = i12(2,i)
            piperp(3,i) = i12(3,i)
            done = .true.
c
c     any non-alkyne atom bonded to exactly two atoms
c
         else if (attach.eq.2 .and. atmnum.ne.6) then
            piperp(1,i) = i
            piperp(2,i) = i12(1,i)
            piperp(3,i) = i12(2,i)
            done = .true.
c
c     atom bonded to four different atoms (usually two lone
c     pairs and two "real" atoms); use the "real" atoms
c
         else if (attach .eq. 4) then
            piperp(1,i) = i
            do j = 1, n12(i)
               trial = i12(j,i)
               if (atomic(trial) .ne. 0) then
                  if (piperp(2,i) .eq. 0) then
                     piperp(2,i) = trial
                  else
                     piperp(3,i) = trial
                     done = .true.
                  end if
               end if
            end do
c
c     "carbonyl"-type oxygen atom, third atom is any atom
c     attached to the "carbonyl carbon"; fails for ketenes
c
         else if (attach.eq.1 .and. atmnum.eq.8) then
            alpha = i12(1,i)
            beta = i12(1,alpha)
            if (beta .eq. i)  beta = i12(2,alpha)
            piperp(1,i) = i
            piperp(2,i) = alpha
            piperp(3,i) = beta
            done = .true.
c
c     an sp nitrogen atom, third atom must be a gamma atom
c
         else if (attach.eq.1 .and. atmnum.eq.7) then
            alpha = i12(1,i)
            do j = 1, n12(alpha)
               trial = i12(j,alpha)
               if (trial.ne.i .and. listpi(trial) .and.
     &             n12(trial).eq.3) then
                  beta = trial
                  done = .true.
               end if
            end do
            gamma = i12(1,beta)
            if (gamma .eq. alpha)  gamma = i12(2,beta)
            piperp(1,i) = i
            piperp(2,i) = alpha
            piperp(3,i) = gamma
c
c     an sp carbon atom; third atom must be an atom attached
c     to the non-sp piatom bonded to the original carbon
c
         else if (attach.eq.2 .and. atmnum.eq.6) then
            alpha = i12(1,i)
            if ((n12(alpha).eq.2 .and. atomic(alpha).eq.6) .or.
     &          (n12(alpha).eq.1 .and. atomic(alpha).eq.7))
     &         alpha = i12(2,i)
            do j = 1, n12(i)
               trial = i12(j,i)
               if (trial.ne.i .and. trial.ne.alpha .and.
     &             listpi(trial) .and. n12(trial).eq.3) then
                  beta = trial
                  done = .true.
               end if
            end do
            do j = 1, n12(alpha)
               trial = i12(j,alpha)
               if (trial.ne.i .and. trial.ne.alpha .and.
     &             listpi(trial) .and. n12(trial).eq.3) then
                  beta = trial
                  done = .true.
               end if
            end do
            gamma = i12(1,beta)
            if (gamma.eq.i .or. gamma.eq.alpha)  gamma = i12(2,beta)
            piperp(1,i) = i
            piperp(2,i) = alpha
            piperp(3,i) = gamma
         end if
c
c     quit if the p-orbital plane remains undefined
c
         if (.not. done) then
            write (iout,10)  i
   10       format(/,' PIPLANE  --  Failure to Define',
     &                ' p-Orbital Plane for Atom',i6)
            call fatal
         end if
      end do
      return
      end
