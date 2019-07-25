c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 1990 by Jay William Ponder                    ##
c     ##  COPYRIGHT (C) 1985 by Scripps Clinic & Research Foundation  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine connolly  --  analytical surface area & volume  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "connolly" uses the algorithms from the AMS/VAM programs of
c     Michael Connolly to compute the analytical molecular surface
c     area and volume of a collection of spherical atoms; thus
c     it implements Fred Richards' molecular surface definition as
c     a set of analytically defined spherical and toroidal polygons
c
c     literature references:
c
c     M. L. Connolly, "Analytical Molecular Surface Calculation",
c     Journal of Applied Crystallography, 16, 548-558 (1983)
c
c     M. L. Connolly, "Computation of Molecular Volume", Journal
c     of the American Chemical Society, 107, 1118-1124 (1985)
c
c     variables only in the Connolly routines:
c
c     na      number of atoms
c     ntt     number of temporary tori
c     nt      number of tori
c     np      number of probe positions
c     nv      number of vertices
c     nen     number of concave edges
c     nfn     number of concave faces
c     nc      number of circles
c     nep     number of convex edges
c     nfs     number of saddle faces
c     ncy     number of cycles
c     fpncy   number of cycles bounding convex face
c     nfp     number of convex faces
c     cynep   number of convex edges in cycle
c
c     axyz    atomic coordinates
c     ar      atomic radii
c     pr      probe radius
c     skip    if true, atom is not used
c     nosurf  if true, atom has no free surface
c     afree   atom free of neighbors
c     abur    atom buried
c
c     acls    begin and end pointers for atoms neighbors
c     cls     atom numbers of neighbors
c     clst    pointer from neighbor to torus
c
c     tta     torus atom numbers
c     ttfe    first edge of each temporary torus
c     ttle    last edge of each temporary torus
c     enext   pointer to next edge of torus
c     ttbur   torus buried
c     ttfree  torus free
c
c     t       torus center
c     tr      torus radius
c     tax     torus axis
c     ta      torus atom numbers
c     tfe     torus first edge
c     tfree   torus free of neighbors
c
c     p       probe coordinates
c     pa      probe atom numbers
c     vxyz    vertex coordinates
c     va      vertex atom number
c     vp      vertex probe number
c     c       circle center
c     cr      circle radius
c     ca      circle atom number
c     ct      circle torus number
c
c     env     concave edge vertex numbers
c     fnen    concave face concave edge numbers
c     epc     convex edge circle number
c     epv     convex edge vertex numbers
c     afe     first convex edge of each atom
c     ale     last convex edge of each atom
c     epnext  pointer to next convex edge of atom
c     fsen    saddle face concave edge numbers
c     fsep    saddle face convex edge numbers
c     cyep    cycle convex edge numbers
c     fpa     atom number of convex face
c     fpcy    convex face cycle numbers
c
c
      subroutine connolly (volume,area,radius,probe,exclude)
      use sizes
      use atoms
      use faces
      implicit none
      integer i
      real*8 volume,area
      real*8 probe,exclude
      real*8 radius(*)
c
c
c     dimensions for arrays used by Connolly routines
c
      maxcls = 75 * n
      maxtt = 40 * n
      maxt = 4 * n
      maxp = 4 * n
      maxv = 12 * n
      maxen = 12 * n
      maxfn = 4 * n
      maxc = 8 * n
      maxep = 12 * n
      maxfs = 6 * n
      maxcy = 3 * n
      mxcyep = 30
      maxfp = 2 * n
      mxfpcy = 10
c
c     perform dynamic allocation of some global arrays
c
      if (.not. allocated(ar))  allocate (ar(n))
      if (.not. allocated(axyz))  allocate (axyz(3,n))
      if (.not. allocated(skip))  allocate (skip(n))
      if (.not. allocated(nosurf))  allocate (nosurf(n))
      if (.not. allocated(afree))  allocate (afree(n))
      if (.not. allocated(abur))  allocate (abur(n))
      if (.not. allocated(cls))  allocate (cls(maxcls))
      if (.not. allocated(clst))  allocate (clst(maxcls))
      if (.not. allocated(acls))  allocate (acls(2,n))
      if (.not. allocated(ttfe))  allocate (ttfe(maxtt))
      if (.not. allocated(ttle))  allocate (ttle(maxtt))
      if (.not. allocated(enext))  allocate (enext(maxen))
      if (.not. allocated(tta))  allocate (tta(2,maxtt))
      if (.not. allocated(ttbur))  allocate (ttbur(maxtt))
      if (.not. allocated(ttfree))  allocate (ttfree(maxtt))
      if (.not. allocated(tfe))  allocate (tfe(maxt))
      if (.not. allocated(ta))  allocate (ta(2,maxt))
      if (.not. allocated(tr))  allocate (tr(maxt))
      if (.not. allocated(t))  allocate (t(3,maxt))
      if (.not. allocated(tax))  allocate (tax(3,maxt))
      if (.not. allocated(tfree))  allocate (tfree(maxt))
      if (.not. allocated(pa))  allocate (pa(3,maxp))
      if (.not. allocated(p))  allocate (p(3,maxp))
      if (.not. allocated(va))  allocate (va(maxv))
      if (.not. allocated(vp))  allocate (vp(maxv))
      if (.not. allocated(vxyz))  allocate (vxyz(3,maxv))
      if (.not. allocated(env))  allocate (env(2,maxen))
      if (.not. allocated(fnen))  allocate (fnen(3,maxfn))
      if (.not. allocated(ca))  allocate (ca(maxc))
      if (.not. allocated(ct))  allocate (ct(maxc))
      if (.not. allocated(cr))  allocate (cr(maxc))
      if (.not. allocated(c))  allocate (c(3,maxc))
      if (.not. allocated(epc))  allocate (epc(maxep))
      if (.not. allocated(epv))  allocate (epv(2,maxep))
      if (.not. allocated(afe))  allocate (afe(n))
      if (.not. allocated(ale))  allocate (ale(n))
      if (.not. allocated(epnext))  allocate (epnext(maxep))
      if (.not. allocated(fsen))  allocate (fsen(2,maxfs))
      if (.not. allocated(fsep))  allocate (fsep(2,maxfs))
      if (.not. allocated(cynep))  allocate (cynep(maxcy))
      if (.not. allocated(cyep))  allocate (cyep(mxcyep,maxcy))
      if (.not. allocated(fpa))  allocate (fpa(maxfp))
      if (.not. allocated(fpncy))  allocate (fpncy(maxfp))
      if (.not. allocated(fpcy))  allocate (fpcy(mxfpcy,maxfp))
c
c     set the probe radius and the number of atoms
c
      pr = probe
      na = n
c
c     set atom coordinates and radii, the excluded buffer
c     radius ("exclude") is added to atomic radii
c
      do i = 1, na
         axyz(1,i) = x(i)
         axyz(2,i) = y(i)
         axyz(3,i) = z(i)
         ar(i) = radius(i)
         if (ar(i) .eq. 0.0d0) then
            skip(i) = .true.
         else
            ar(i) = ar(i) + exclude
            skip(i) = .false.
         end if
      end do
c
c     find the analytical volume and surface area
c
c     call wiggle
      call nearby
      call torus
      call place
      call compress
      call saddles
      call contact
      call vam (volume,area)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine wiggle  --  random perturbation of coordinates  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "wiggle" applies a random perturbation to the atomic coordinates
c     to avoid numerical instabilities for symmetrical structures
c
c
      subroutine wiggle
      use faces
      implicit none
      integer i
      real*8 size
      real*8 vector(3)
c
c
c     apply a small perturbation of fixed magnitude to each atom
c
      size = 0.000001d0
      do i = 1, na
         call ranvec (vector)
         axyz(1,i) = axyz(1,i) + size*vector(1)
         axyz(2,i) = axyz(2,i) + size*vector(2)
         axyz(3,i) = axyz(3,i) + size*vector(3)
      end do
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine nearby  --  list of neighboring atom pairs  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "nearby" finds all of the through-space neighbors of each
c     atom for use in surface area and volume calculations
c
c     local variables :
c
c     ico      integer cube coordinates
c     icuptr   pointer to next atom in cube
c     comin    minimum atomic coordinates (cube corner)
c     icube    pointer to first atom in list for cube
c     scube    true if cube contains active atoms
c     sscube   true if cube or adjacent cubes have active atoms
c     itnl     temporary neighbor list, before sorting
c
c
      subroutine nearby
      use faces
      implicit none
      integer maxclsa
      parameter (maxclsa=1000)
      integer i,j,k,m
      integer iptr,juse
      integer i1,j1,k1
      integer iatom,jatom
      integer ici,icj,ick
      integer jci,jcj,jck
      integer jcls,jmin
      integer jmincls,jmold
      integer ncls,nclsa
      integer maxcube
      integer clsa(maxclsa)
      integer itnl(maxclsa)
      integer, allocatable :: icuptr(:)
      integer, allocatable :: ico(:,:)
      integer, allocatable :: icube(:,:,:)
      real*8 radmax,width
      real*8 sum,sumi
      real*8 dist2,d2,r2
      real*8 vect1,vect2,vect3
      real*8 comin(3)
      logical, allocatable :: scube(:,:,:)
      logical, allocatable :: sscube(:,:,:)
c
c
c     ignore all atoms that are completely inside another atom;
c     may give nonsense results if this step is not taken
c
      do i = 1, na-1
         if (.not. skip(i)) then
            do j = i+1, na
               d2 = dist2(axyz(1,i),axyz(1,j))
               r2 = (ar(i) - ar(j))**2
               if (.not.skip(j) .and. d2.lt.r2) then
                  if (ar(i) .lt. ar(j)) then
                     skip(i) = .true.
                  else
                     skip(j) = .true.
                  end if
               end if
            end do
         end if
      end do
c
c     check for new coordinate minima and radii maxima
c
      radmax = 0.0d0
      do k = 1, 3
         comin(k) = axyz(k,1)
      end do
      do i = 1, na
         do k = 1, 3
            if (axyz(k,i) .lt. comin(k))  comin(k) = axyz(k,i)
         end do
         if (ar(i) .gt. radmax)  radmax = ar(i)
      end do
c
c     calculate width of cube from maximum
c     atom radius and probe radius
c
      width = 2.0d0 * (radmax+pr)
c
c     perform dynamic allocation of some local arrays
c
      maxcube = 40
      allocate (icuptr(na))
      allocate (ico(3,na))
      allocate (icube(maxcube,maxcube,maxcube))
      allocate (scube(maxcube,maxcube,maxcube))
      allocate (sscube(maxcube,maxcube,maxcube))
c
c     set up cube arrays; first the integer coordinate arrays
c
      do i = 1, na
         do k = 1, 3
            ico(k,i) = int((axyz(k,i)-comin(k))/width) + 1
            if (ico(k,i) .lt. 1) then
               call cerror ('Cube Coordinate Too Small')
            else if (ico(k,i) .gt. maxcube) then
               call cerror ('Cube Coordinate Too Large')
            end if
         end do
      end do
c
c     initialize head pointer and srn=2 arrays
c
      do i = 1, maxcube
         do j = 1, maxcube
            do k = 1, maxcube
               icube(i,j,k) = 0
               scube(i,j,k) = .false.
               sscube(i,j,k) = .false.
            end do
         end do
      end do
c
c     initialize linked list pointers
c
      do i = 1, na
         icuptr(i) = 0
      end do
c
c     set up head and later pointers for each atom
c
      do iatom = 1, na
c
c     skip atoms with surface request numbers of zero
c
         if (skip(iatom))  goto 30
         i = ico(1,iatom)
         j = ico(2,iatom)
         k = ico(3,iatom)
         if (icube(i,j,k) .le. 0) then
c
c     first atom in this cube
c
            icube(i,j,k) = iatom
         else
c
c     add to end of linked list
c
            iptr = icube(i,j,k)
   10       continue
c
c     check for duplicate atoms, turn off one of them
c
            if (dist2(axyz(1,iatom),axyz(1,iptr)) .le. 0.0d0) then
               skip(iatom) = .true.
               goto 30
            end if
c
c     move on down the list
c
            if (icuptr(iptr) .le. 0)  goto 20
            iptr = icuptr(iptr)
            goto 10
   20       continue
c
c     store atom number
c
            icuptr(iptr) = iatom
         end if
c
c     check for surfaced atom
c
         if (.not. skip(iatom))  scube(i,j,k) = .true.
   30    continue
      end do
c
c     check if this cube or any adjacent cube has active atoms
c
      do k = 1, maxcube
         do j = 1, maxcube
            do i = 1, maxcube
               if (icube(i,j,k) .ne. 0) then
                  do k1 = max(k-1,1), min(k+1,maxcube)
                     do j1 = max(j-1,1), min(j+1,maxcube)
                        do i1 = max(i-1,1), min(i+1,maxcube)
                           if (scube(i1,j1,k1)) then
                              sscube(i,j,k) = .true.
                           end if
                        end do
                     end do
                  end do
               end if
            end do
         end do
      end do
      ncls = 0
c
c     zero pointers for atom and find its cube
c
      do i = 1, na
         nclsa = 0
         nosurf(i) = skip(i)
         acls(1,i) = 0
         acls(2,i) = 0
         if (skip(i))  goto 70
         ici = ico(1,i)
         icj = ico(2,i)
         ick = ico(3,i)
c
c     skip iatom if its cube and adjoining
c     cubes contain only blockers
c
         if (.not. sscube(ici,icj,ick))  goto 70
         sumi = 2.0d0*pr + ar(i)
c
c     check iatom cube and adjacent cubes for neighboring atoms
c
         do jck = max(ick-1,1), min(ick+1,maxcube)
            do jcj = max(icj-1,1), min(icj+1,maxcube)
               do jci = max(ici-1,1), min(ici+1,maxcube)
                  j = icube(jci,jcj,jck)
   40             continue
c
c     check for end of linked list for this cube
c
                  if (j .le. 0)  goto 60
                  if (i .eq. j)  goto 50
                  if (skip(j))  goto 50
c
c     distance check
c
                  sum = sumi + ar(j)
                  vect1 = abs(axyz(1,j) - axyz(1,i))
                  if (vect1 .ge. sum)  goto 50
                  vect2 = abs(axyz(2,j) - axyz(2,i))
                  if (vect2 .ge. sum)  goto 50
                  vect3 = abs(axyz(3,j) - axyz(3,i))
                  if (vect3 .ge. sum)  goto 50
                  d2 = vect1**2 + vect2**2 + vect3**2
                  if (d2 .ge. sum**2)  goto 50
c
c     atoms are neighbors, save atom number in temporary array
c
                  if (.not. skip(j))  nosurf(i) = .false.
                  nclsa = nclsa + 1
                  if (nclsa .gt. maxclsa) then
                     call cerror ('Too many Neighbors for Atom')
                  end if
                  itnl(nclsa) = j
   50             continue
c
c     get number of next atom in cube
c
                  j = icuptr(j)
                  goto 40
   60             continue
               end do
            end do
         end do
         if (nosurf(i))  goto 70
c
c     set up neighbors arrays with jatom in increasing order
c
         jmold = 0
         do juse = 1, nclsa
            jmin = na + 1
            do jcls = 1, nclsa
c
c     don't use ones already sorted
c
               if (itnl(jcls) .gt. jmold) then
                  if (itnl(jcls) .lt. jmin) then
                     jmin = itnl(jcls)
                     jmincls = jcls
                  end if
               end if
            end do
            jmold = jmin
            jcls = jmincls
            jatom = itnl(jcls)
            clsa(juse) = jatom
         end do
c
c     set up pointers to first and last neighbors of atom
c
         if (nclsa .gt. 0) then
            acls(1,i) = ncls + 1
            do m = 1, nclsa
               ncls = ncls + 1
               if (ncls .gt. maxcls) then
                  call cerror ('Too many Neighboring Atom Pairs')
               end if
               cls(ncls) = clsa(m)
            end do
            acls(2,i) = ncls
         end if
   70    continue
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (icuptr)
      deallocate (ico)
      deallocate (icube)
      deallocate (scube)
      deallocate (sscube)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine torus  --  position of each temporary torus  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "torus" sets a list of all of the temporary torus positions
c     by testing for a torus between each atom and its neighbors
c
c
      subroutine torus
      use faces
      implicit none
      integer ia,ja,jn
      integer ibeg,iend
      real*8 ttr
      real*8 tt(3)
      real*8 ttax(3)
      logical ttok
c
c
c     no torus is possible if there is only one atom
c
      ntt = 0
      do ia = 1, na
         afree(ia) = .true.
      end do
      if (na .le. 1)  return
c
c     get begin and end pointers to neighbors of this atom
c
      do ia = 1, na
         if (.not. nosurf(ia)) then
            ibeg = acls(1,ia)
            iend = acls(2,ia)
c
c     check for no neighbors
c
            if (ibeg .gt. 0) then
               do jn = ibeg, iend
c
c     clear pointer from neighbor to torus
c
                  clst(jn) = 0
c
c     get atom number of neighbor
c
                  ja = cls(jn)
c
c     don't create torus twice
c
                  if (ja .ge. ia) then
c
c     do some solid geometry
c
                     call gettor (ia,ja,ttok,tt,ttr,ttax)
                     if (ttok) then
c
c     we have a temporary torus, set up variables
c
                        ntt = ntt + 1
                        if (ntt .gt. maxtt) then
                           call cerror ('Too many Temporary Tori')
                        end if
c
c     mark both atoms not free
c
                        afree(ia) = .false.
                        afree(ja) = .false.
                        tta(1,ntt) = ia
                        tta(2,ntt) = ja
c
c     pointer from neighbor to torus
c
                        clst(jn) = ntt
c
c     initialize torus as both free and buried
c
                        ttfree(ntt) = .true.
                        ttbur(ntt) = .true.
c
c     clear pointers from torus to first and last concave edges
c
                        ttfe(ntt) = 0
                        ttle(ntt) = 0
                     end if
                  end if
               end do
            end if
         end if
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine place  --  locate positions of the probe sites  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "place" finds the probe sites by putting the probe sphere
c     tangent to each triple of neighboring atoms
c
c
      subroutine place
      use sizes
      use faces
      implicit none
      integer maxmnb
      parameter (maxmnb=500)
      integer k,ke,kv
      integer l,l1,l2
      integer ia,ja,ka
      integer ik,ip,jk
      integer km,la,lm
      integer lkf,itt,nmnb
      integer iend,jend
      integer iptr,jptr
      integer mnb(maxmnb)
      integer ikt(maxmnb)
      integer jkt(maxmnb)
      integer lkcls(maxmnb)
      real*8 dist2,d2,det
      real*8 hij,hijk
      real*8 uij(3),uijk(3)
      real*8 bij(3),bijk(3)
      real*8 aijk(3),pijk(3)
      real*8 tempv(3)
      real*8 discls(maxmnb)
      real*8 sumcls(maxmnb)
      logical tb,ttok,prbok
c
c
c     no possible placement if there are no temporary tori
c
      np = 0
      nfn = 0
      nen = 0
      nv = 0
      if (ntt .le. 0)  return
c
c     consider each torus in turn
c
      do itt = 1, ntt
c
c     get atom numbers
c
         ia = tta(1,itt)
         ja = tta(2,itt)
c
c     form mutual neighbor list; clear number
c     of mutual neighbors of atoms ia and ja
c
         nmnb = 0
c
c     get begin and end pointers for each atom's neighbor list
c
         iptr = acls(1,ia)
         jptr = acls(1,ja)
         if (iptr.le.0 .or. jptr.le.0)  goto 130
         iend = acls(2,ia)
         jend = acls(2,ja)
c
c     collect mutual neighbors
c
   10    continue
c
c     check for end of loop
c
         if (iptr .gt. iend)  goto 40
         if (jptr .gt. jend)  goto 40
c
c     go move the lagging pointer
c
         if (cls(iptr) .lt. cls(jptr))  goto 20
         if (cls(jptr) .lt. cls(iptr))  goto 30
c
c     both point at same neighbor; one more mutual neighbor
c     save atom number of mutual neighbor
c
         nmnb = nmnb + 1
         if (nmnb .gt. maxmnb) then
            call cerror ('Too many Mutual Neighbors')
         end if
         mnb(nmnb) = cls(iptr)
c
c     save pointers to second and third tori
c
         ikt(nmnb) = clst(iptr)
         jkt(nmnb) = clst(jptr)
   20    continue
c
c     increment pointer to ia atom neighbors
c
         iptr = iptr + 1
         goto 10
   30    continue
c
c     increment pointer to ja atom neighbors
c
         jptr = jptr + 1
         goto 10
   40    continue
c
c     we have all the mutual neighbors of ia and ja
c     if no mutual neighbors, skip to end of loop
c
         if (nmnb .le. 0) then
            ttbur(itt) = .false.
            goto 130
         end if
         call gettor (ia,ja,ttok,bij,hij,uij)
         do km = 1, nmnb
            ka = mnb(km)
            discls(km) = dist2(bij,axyz(1,ka))
            sumcls(km) = (pr+ar(ka))**2
c
c     initialize link to next farthest out neighbor
c
            lkcls(km) = 0
         end do
c
c     set up a linked list of neighbors in order of
c     increasing distance from ia-ja torus center
c
         lkf = 1
         if (nmnb .le. 1)  goto 70
c
c     put remaining neighbors in linked list at proper position
c
         do l = 2, nmnb
            l1 = 0
            l2 = lkf
   50       continue
            if (discls(l) .lt. discls(l2))  goto 60
            l1 = l2
            l2 = lkcls(l2)
            if (l2 .ne. 0)  goto 50
   60       continue
c
c     add to list
c
            if (l1 .eq. 0) then
               lkf = l
               lkcls(l) = l2
            else
               lkcls(l1) = l
               lkcls(l) = l2
            end if
         end do
   70    continue
c
c     loop thru mutual neighbors
c
         do km = 1, nmnb
c
c     get atom number of neighbors
c
            ka = mnb(km)
            if (skip(ia) .and. skip(ja) .and. skip(ka))  goto 120
c
c     get tori numbers for neighbor
c
            ik = ikt(km)
            jk = jkt(km)
c
c     possible new triple, do some geometry to
c     retrieve saddle center, axis and radius
c
            call getprb (ia,ja,ka,prbok,tb,bijk,hijk,uijk)
            if (tb) then
               ttbur(itt) = .true.
               ttfree(itt) = .false.
               goto 120
            end if
c
c     no duplicate triples
c
            if (ka .lt. ja)  goto 120
c
c     check whether any possible probe positions
c
            if (.not. prbok)  goto 120
c
c     altitude vector
c
            do k = 1, 3
               aijk(k) = hijk * uijk(k)
            end do
c
c     we try two probe placements
c
            do ip = 1, 2
               do k = 1, 3
                  if (ip .eq. 1) then
                     pijk(k) = bijk(k) + aijk(k)
                  else
                     pijk(k) = bijk(k) - aijk(k)
                  end if
               end do
c
c     mark three tori not free
c
               ttfree(itt) = .false.
               ttfree(ik) = .false.
               ttfree(jk) = .false.
c
c     check for collisions
c
               lm = lkf
   80          continue
               if (lm .le. 0)  goto 100
c
c     get atom number of mutual neighbor
c
               la = mnb(lm)
c
c     must not equal third atom
c
               if (la .eq. ka)  goto 90
c
c     compare distance to sum of radii
c
               d2 = dist2(pijk,axyz(1,la))
               if (d2 .le. sumcls(lm))  goto 110
   90          continue
               lm = lkcls(lm)
               goto 80
  100          continue
c
c     we have a new probe position
c
               np = np + 1
               if (np .gt. maxp) then
                  call cerror ('Too many Probe Positions')
               end if
c
c     mark three tori not buried
c
               ttbur(itt) = .false.
               ttbur(ik) = .false.
               ttbur(jk) = .false.
c
c     store probe center
c
               do k = 1, 3
                  p(k,np) = pijk(k)
               end do
c
c     calculate vectors from probe to atom centers
c
               if (nv+3 .gt. maxv)  call cerror ('Too many Vertices')
               do k = 1, 3
                  vxyz(k,nv+1) = axyz(k,ia) - p(k,np)
                  vxyz(k,nv+2) = axyz(k,ja) - p(k,np)
                  vxyz(k,nv+3) = axyz(k,ka) - p(k,np)
               end do
c
c     calculate determinant of vectors defining triangle
c
               det = vxyz(1,nv+1)*vxyz(2,nv+2)*vxyz(3,nv+3)
     &                  + vxyz(1,nv+2)*vxyz(2,nv+3)*vxyz(3,nv+1)
     &                  + vxyz(1,nv+3)*vxyz(2,nv+1)*vxyz(3,nv+2)
     &                  - vxyz(1,nv+3)*vxyz(2,nv+2)*vxyz(3,nv+1)
     &                  - vxyz(1,nv+2)*vxyz(2,nv+1)*vxyz(3,nv+3)
     &                  - vxyz(1,nv+1)*vxyz(2,nv+3)*vxyz(3,nv+2)
c
c     now add probe coordinates to vertices
c
               do k = 1, 3
                  vxyz(k,nv+1) = p(k,np) + vxyz(k,nv+1)*pr/(ar(ia)+pr)
                  vxyz(k,nv+2) = p(k,np) + vxyz(k,nv+2)*pr/(ar(ja)+pr)
                  vxyz(k,nv+3) = p(k,np) + vxyz(k,nv+3)*pr/(ar(ka)+pr)
               end do
c
c     want the concave face to have counter-clockwise orientation
c
               if (det .gt. 0.0d0) then
c
c     swap second and third vertices
c
                  do k = 1, 3
                     tempv(k) = vxyz(k,nv+2)
                     vxyz(k,nv+2) = vxyz(k,nv+3)
                     vxyz(k,nv+3) = tempv(k)
                  end do
c
c     set up pointers from probe to atoms
c
                  pa(1,np) = ia
                  pa(2,np) = ka
                  pa(3,np) = ja
c
c     set up pointers from vertices to atoms
c
                  va(nv+1) = ia
                  va(nv+2) = ka
                  va(nv+3) = ja
c
c     insert concave edges into linked lists for appropriate tori
c
                  call inedge (nen+1,ik)
                  call inedge (nen+2,jk)
                  call inedge (nen+3,itt)
               else
c
c     similarly, if face already counter clockwise
c
                  pa(1,np) = ia
                  pa(2,np) = ja
                  pa(3,np) = ka
                  va(nv+1) = ia
                  va(nv+2) = ja
                  va(nv+3) = ka
                  call inedge (nen+1,itt)
                  call inedge (nen+2,jk)
                  call inedge (nen+3,ik)
               end if
c
c     set up pointers from vertices to probe
c
               do kv = 1, 3
                  vp(nv+kv) = np
               end do
c
c     set up concave edges and concave face
c
               if (nen+3 .gt. maxen) then
                  call cerror ('Too many Concave Edges')
               end if
c
c     edges point to vertices
c
               env(1,nen+1) = nv+1
               env(2,nen+1) = nv+2
               env(1,nen+2) = nv+2
               env(2,nen+2) = nv+3
               env(1,nen+3) = nv+3
               env(2,nen+3) = nv+1
               if (nfn+1 .gt. maxfn) then
                  call cerror ('Too many Concave Faces')
               end if
c
c     face points to edges
c
               do ke = 1, 3
                  fnen(ke,nfn+1) = nen + ke
               end do
c
c     increment counters for number of faces, edges and vertices
c
               nfn = nfn + 1
               nen = nen + 3
               nv = nv + 3
  110          continue
            end do
  120       continue
         end do
  130    continue
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine inedge  --  manage linked list of torus edges  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "inedge" inserts a concave edge into the
c     linked list for its temporary torus
c
c
      subroutine inedge (ien,itt)
      use sizes
      use faces
      implicit none
      integer ien,itt,iepen
c
c
c     check for a serious error in the calling arguments
c
      if (ien .le. 0)  call cerror ('Bad Edge Number in INEDGE')
      if (itt .le. 0)  call cerror ('Bad Torus Number in INEDGE')
c
c     set beginning of list or add to end
c
      if (ttfe(itt) .eq. 0) then
         ttfe(itt) = ien
         enext(ien) = 0
         ttle(itt) = ien
      else
         iepen = ttle(itt)
         enext(iepen) = ien
         enext(ien) = 0
         ttle(itt) = ien
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine compress  --  condense temporary to final tori  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "compress" transfers only the non-buried tori from
c     the temporary tori arrays to the final tori arrays
c
c
      subroutine compress
      use sizes
      use faces
      implicit none
      integer itt,ia,ja
      integer iptr,ned
      integer ip1,ip2
      integer iv1,iv2
      logical ttok
c
c
c     initialize the number of nonburied tori
c
      nt = 0
      if (ntt .le. 0)  return
c
c     if torus is free, then it is not buried;
c     skip to end of loop if buried torus
c
      do itt = 1, ntt
         if (ttfree(itt))  ttbur(itt) = .false.
         if (.not. ttbur(itt)) then
c
c     first, transfer information
c
            nt = nt + 1
            if (nt .gt. maxt)  call cerror ('Too many NonBuried Tori')
            ia = tta(1,itt)
            ja = tta(2,itt)
            call gettor (ia,ja,ttok,t(1,nt),tr(nt),tax(1,nt))
            ta(1,nt) = ia
            ta(2,nt) = ja
            tfree(nt) = ttfree(itt)
            tfe(nt) = ttfe(itt)
c
c     special check for inconsistent probes
c
            iptr = tfe(nt)
            ned = 0
            do while (iptr .ne. 0)
               ned = ned + 1
               iptr = enext(iptr)
            end do
            if (mod(ned,2) .ne. 0) then
               iptr = tfe(nt)
               do while (iptr .ne. 0)
                  iv1 = env(1,iptr)
                  iv2 = env(2,iptr)
                  ip1 = vp(iv1)
                  ip2 = vp(iv2)
                  call cerror ('Odd Torus for Probes IP1 and IP2')
                  iptr = enext(iptr)
               end do
            end if
         end if
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine saddles  --  builds saddle pieces from tori  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "saddles" constructs circles, convex edges and saddle faces
c
c
      subroutine saddles
      use sizes
      use faces
      use math
      implicit none
      integer maxent
      parameter (maxent=500)
      integer k,ia,in,ip
      integer it,iv,itwo
      integer ien,ient,nent
      integer l1,l2,m1,n1
      integer ten(maxent)
      integer nxtang(maxent)
      real*8 triple,factor
      real*8 dtev,dt
      real*8 atvect(3)
      real*8 teang(maxent)
      real*8 tev(3,maxent)
      logical sdstrt(maxent)
c
c
c     zero the number of circles, convex edges and saddle faces
c
      nc = 0
      nep = 0
      nfs = 0
      do ia = 1, na
         afe(ia) = 0
         ale(ia) = 0
         abur(ia) = .true.
      end do
c
c     no saddle faces if no tori
c
      if (nt .lt. 1)  return
c
c     cycle through tori
c
      do it = 1, nt
         if (skip(ta(1,it)) .and. skip(ta(2,it)))  goto 80
c
c     set up two circles
c
         do in = 1, 2
            ia = ta(in,it)
c
c     mark atom not buried
c
            abur(ia) = .false.
c
c     vector from atom to torus center
c
            do k = 1, 3
               atvect(k) = t(k,it) - axyz(k,ia)
            end do
            factor = ar(ia) / (ar(ia)+pr)
c
c     one more circle
c
            nc = nc + 1
            if (nc .gt. maxc)  call cerror ('Too many Circles')
c
c     circle center
c
            do k = 1, 3
               c(k,nc) = axyz(k,ia) + factor*atvect(k)
            end do
c
c     pointer from circle to atom and to torus
c
            ca(nc) = ia
            ct(nc) = it
c
c     circle radius
c
            cr(nc) = factor * tr(it)
         end do
c
c     skip to special code if free torus
c
         if (tfree(it))  goto 70
c
c     now we collect all the concave edges for this torus;
c     for each concave edge, calculate vector from torus center
c     thru probe center and the angle relative to first such vector
c
c     clear the number of concave edges for torus
c
         nent = 0
c
c     pointer to start of linked list
c
         ien = tfe(it)
   10    continue
c
c     finished if concave edge pointer is zero
c
         if (ien .le. 0)  goto 20
c
c     one more concave edge
c
         nent = nent + 1
         if (nent .gt. maxent) then
            call cerror ('Too many Edges for Torus')
         end if
c
c     first vertex of edge
c
         iv = env(1,ien)
c
c     probe number of vertex
c
         ip = vp(iv)
         do k = 1, 3
            tev(k,nent) = p(k,ip) - t(k,it)
         end do
         dtev = 0.0d0
         do k = 1, 3
            dtev = dtev + tev(k,nent)**2
         end do
         if (dtev .le. 0.0d0)  call cerror ('Probe on Torus Axis')
         dtev = sqrt(dtev)
         do k = 1, 3
            tev(k,nent) = tev(k,nent) / dtev
         end do
c
c     store concave edge number
c
         ten(nent) = ien
         if (nent .gt. 1) then
c
c     calculate angle between this vector and first vector
c
            dt = 0.0d0
            do k = 1, 3
               dt = dt + tev(k,1)*tev(k,nent)
            end do
c
c     be careful
c
            if (dt .gt. 1.0d0)  dt = 1.0d0
            if (dt .lt. -1.0d0)  dt = -1.0d0
c
c     store angle
c
            teang(nent) = acos(dt)
c
c     get the sign right
c
            if (triple(tev(1,1),tev(1,nent),tax(1,it)) .lt. 0.0d0) then
               teang(nent) = 2.0d0*pi - teang(nent)
            end if
         else
            teang(1) = 0.0d0
         end if
c
c     saddle face starts with this edge if it points parallel
c     to torus axis vector (which goes from first to second atom)
c
         sdstrt(nent) = (va(iv) .eq. ta(1,it))
c
c     next edge in list
c
         ien = enext(ien)
         goto 10
   20    continue
         if (nent .le. 0) then
            call cerror ('No Edges for Non-free Torus')
         end if
         itwo = 2
         if (mod(nent,itwo) .ne. 0) then
            call cerror ('Odd Number of Edges for Torus')
         end if
c
c     set up linked list of concave edges in order
c     of increasing angle around the torus axis;
c     clear second linked (angle-ordered) list pointers
c
         do ient = 1, nent
            nxtang(ient) = 0
         end do
         do ient = 2, nent
c
c     we have an entry to put into linked list
c     search for place to put it
c
            l1 = 0
            l2 = 1
   30       continue
            if (teang(ient) .lt. teang(l2))  goto 40
c
c     not yet, move along
c
            l1 = l2
            l2 = nxtang(l2)
            if (l2 .ne. 0)  goto 30
   40       continue
c
c     we are at end of linked list or between l1 and l2;
c     insert edge
c
            if (l1 .le. 0)  call cerror ('Logic Error in SADDLES')
            nxtang(l1) = ient
            nxtang(ient) = l2
         end do
c
c     collect pairs of concave edges into saddles
c     create convex edges while you're at it
c
         l1 = 1
   50    continue
         if (l1 .le. 0)  goto 60
c
c     check for start of saddle
c
         if (sdstrt(l1)) then
c
c     one more saddle face
c
            nfs = nfs + 1
            if (nfs .gt. maxfs)  call cerror ('Too many Saddle Faces')
c
c     get edge number
c
            ien = ten(l1)
c
c     first concave edge of saddle
c
            fsen(1,nfs) = ien
c
c     one more convex edge
c
            nep = nep + 1
            if (nep .gt. maxep)  call cerror ('Too many Convex Edges')
c
c     first convex edge points to second circle
c
            epc(nep) = nc
c
c     atom circle lies on
c
            ia = ca(nc)
c
c     insert convex edge into linked list for atom
c
            call ipedge (nep,ia)
c
c     first vertex of convex edge is second vertex of concave edge
c
            epv(1,nep) = env(2,ien)
c
c     first convex edge of saddle
c
            fsep(1,nfs) = nep
c
c     one more convex edge
c
            nep = nep + 1
            if (nep .gt. maxep)  call cerror ('Too many Convex Edges')
c
c     second convex edge points to first circle
c
            epc(nep) = nc - 1
            ia = ca(nc-1)
c
c     insert convex edge into linked list for atom
c
            call ipedge (nep,ia)
c
c     second vertex of second convex edge
c     is first vertex of first concave edge
c
            epv(2,nep) = env(1,ien)
            l1 = nxtang(l1)
c
c     wrap around
c
            if (l1 .le. 0)  l1 = 1
            if (sdstrt(l1)) then
               m1 = nxtang(l1)
               if (m1 .le. 0)  m1 = 1
               if (sdstrt(m1))  call cerror ('Three Starts in a Row')
               n1 = nxtang(m1)
c
c     the old switcheroo
c
               nxtang(l1) = n1
               nxtang(m1) = l1
               l1 = m1
            end if
            ien = ten(l1)
c
c     second concave edge for saddle face
c
            fsen(2,nfs) = ien
c
c     second vertex of first convex edge is
c     first vertex of second concave edge
c
            epv(2,nep-1) = env(1,ien)
c
c     first vertex of second convex edge is
c     second vertex of second concave edge
c
            epv(1,nep) = env(2,ien)
            fsep(2,nfs) = nep
c
c     quit if we have wrapped around to first edge
c
            if (l1 .eq. 1)  goto 60
         end if
c
c     next concave edge
c
         l1 = nxtang(l1)
         goto 50
   60    continue
         goto 80
c
c     free torus
c
   70    continue
c
c     set up entire circles as convex edges for new saddle surface;
c     one more saddle face
c
         nfs = nfs + 1
         if (nfs .gt. maxfs)  call cerror ('Too many Saddle Faces')
c
c     no concave edges for saddle
c
         fsen(1,nfs) = 0
         fsen(2,nfs) = 0
c
c     one more convex edge
c
         nep = nep + 1
         ia = ca(nc)
c
c     insert convex edge into linked list for atom
c
         call ipedge (nep,ia)
c
c     no vertices for convex edge
c
         epv(1,nep) = 0
         epv(2,nep) = 0
c
c     pointer from convex edge to second circle
c
         epc(nep) = nc
c
c     first convex edge for saddle face
c
         fsep(1,nfs) = nep
c
c     one more convex edge
c
         nep = nep + 1
         ia = ca(nc-1)
c
c     insert second convex edge into linked list
c
         call ipedge (nep,ia)
c
c     no vertices for convex edge
c
         epv(1,nep) = 0
         epv(2,nep) = 0
c
c     convex edge points to first circle
c
         epc(nep) = nc - 1
c
c     second convex edge for saddle face
c
         fsep(2,nfs) = nep
c
c     nothing to do for buried torus
c
   80    continue
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine gettor  --  test torus site between two atoms  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "gettor" tests for a possible torus position at the interface
c     between two atoms, and finds the torus radius, center and axis
c
c
      subroutine gettor (ia,ja,ttok,torcen,torad,torax)
      use sizes
      use faces
      implicit none
      integer k,ia,ja
      real*8 dist2,dij
      real*8 temp,temp1
      real*8 temp2
      real*8 torad
      real*8 torcen(3)
      real*8 torax(3)
      real*8 vij(3)
      real*8 uij(3)
      real*8 bij(3)
      logical ttok
c
c
c     get the distance between the two atoms
c
      ttok = .false.
      dij = sqrt(dist2(axyz(1,ia),axyz(1,ja)))
c
c     find a unit vector along interatomic (torus) axis
c
      do k = 1, 3
         vij(k) = axyz(k,ja) - axyz(k,ia)
         uij(k) = vij(k) / dij
      end do
c
c     find coordinates of the center of the torus
c
      temp = 1.0d0 + ((ar(ia)+pr)**2-(ar(ja)+pr)**2)/dij**2
      do k = 1, 3
         bij(k) = axyz(k,ia) + 0.5d0*vij(k)*temp
      end do
c
c     skip if atoms too far apart (should not happen)
c
      temp1 = (ar(ia)+ar(ja)+2.0d0*pr)**2 - dij**2
      if (temp1 .ge. 0.0d0) then
c
c     skip if one atom is inside the other
c
         temp2 = dij**2 - (ar(ia)-ar(ja))**2
         if (temp2 .ge. 0.0d0) then
c
c     store the torus radius, center and axis
c
            ttok = .true.
            torad = sqrt(temp1*temp2) / (2.0d0*dij)
            do k = 1, 3
               torcen(k) = bij(k)
               torax(k) = uij(k)
            end do
         end if
      end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine getprb  --  test probe site between three atoms  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "getprb" tests for a possible probe position at the interface
c     between three neighboring atoms
c
c
      subroutine getprb (ia,ja,ka,prbok,tb,bijk,hijk,uijk)
      use sizes
      use faces
      implicit none
      integer k,ia,ja,ka
      real*8 dot,dotijk,dotut
      real*8 wijk,swijk,fact
      real*8 dist2,dat2
      real*8 rad,rad2
      real*8 dba,rip2,hijk
      real*8 rij,rik
      real*8 uij(3),uik(3)
      real*8 uijk(3),utb(3)
      real*8 tij(3),tik(3)
      real*8 bijk(3),tijik(3)
      logical prbok,tb,tok
c
c
c     initialize, then check torus over atoms "ia" and "ja"
c
      prbok = .false.
      tb = .false.
      call gettor (ia,ja,tok,tij,rij,uij)
      if (.not. tok)  return
      dat2 = dist2(axyz(1,ka),tij)
      rad2 = (ar(ka)+pr)**2 - rij**2
c
c     if "ka" less than "ja", then all we care about
c     is whether the torus is buried
c
      if (ka .lt. ja) then
         if (rad2 .le. 0.0d0)  return
         if (dat2 .gt. rad2)  return
      end if
      call gettor (ia,ka,tok,tik,rik,uik)
      if (.not. tok)  return
      dotijk = dot(uij,uik)
      if (dotijk .gt. 1.0d0)  dotijk = 1.0d0
      if (dotijk .lt. -1.0d0)  dotijk = -1.0d0
      wijk = acos(dotijk)
      swijk = sin(wijk)
c
c     if the three atoms are colinear, then there is no
c     probe placement; but we still care whether the torus
c     is buried by atom "k"
c
      if (swijk .eq. 0.0d0) then
         tb = (rad2.gt.0.0d0 .and. dat2.le.rad2)
         return
      end if
      call vcross (uij,uik,uijk)
      do k = 1, 3
         uijk(k) = uijk(k) / swijk
      end do
      call vcross (uijk,uij,utb)
      do k = 1, 3
         tijik(k) = tik(k) - tij(k)
      end do
      dotut = dot(uik,tijik)
      fact = dotut / swijk
      do k = 1, 3
         bijk(k) = tij(k) + utb(k)*fact
      end do
      dba = dist2(axyz(1,ia),bijk)
      rip2 = (ar(ia) + pr)**2
      rad = rip2 - dba
      if (rad .lt. 0.0d0) then
         tb = (rad2.gt.0.0d0 .and. dat2.le.rad2)
      else
         prbok = .true.
         hijk = sqrt(rad)
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine ipedge  --  manage linked list of convex edges  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ipedge" inserts convex edge into linked list for atom
c
c
      subroutine ipedge (iep,ia)
      use sizes
      use faces
      implicit none
      integer iep,ia,iepen
c
c
c     first, check for an error condition
c
      if (iep .le. 0)  call cerror ('Bad Edge Number in IPEDGE')
      if (ia .le. 0)  call cerror ('Bad Atom Number in IPEDGE')
c
c     set beginning of list or add to end
c
      if (afe(ia) .eq. 0) then
         afe(ia) = iep
         epnext(iep) = 0
         ale(ia) = iep
      else
         iepen = ale(ia)
         epnext(iepen) = iep
         epnext(iep) = 0
         ale(ia) = iep
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine contact  --  builds exposed contact surfaces  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "contact" constructs the contact surface, cycles and convex faces
c
c
      subroutine contact
      use sizes
      use faces
      implicit none
      integer maxepa,maxcypa
      parameter (maxepa=300)
      parameter (maxcypa=100)
      integer i,k,ia,ia2,it
      integer iep,ic,jc,jcy
      integer nepa,iepa,jepa
      integer ncypa,icya,jcya,kcya
      integer ncyep,icyep,jcyep
      integer ncyold,nused,lookv
      integer aic(maxepa)
      integer aia(maxepa)
      integer aep(maxepa)
      integer av(2,maxepa)
      integer ncyepa(maxcypa)
      integer cyepa(mxcyep,maxcypa)
      real*8 anorm,anaa,factor
      real*8 acvect(3,maxepa)
      real*8 aavect(3,maxepa)
      real*8 pole(3),unvect(3)
      real*8 acr(maxepa)
      logical ptincy,epused(maxepa)
      logical cycy(maxcypa,maxcypa)
      logical cyused(maxcypa)
      logical samef(maxcypa,maxcypa)
c
c
c     zero out the number of cycles and convex faces
c
      ncy = 0
      nfp = 0
c
c     mark all free atoms not buried
c
      do ia = 1, na
         if (afree(ia))  abur(ia) = .false.
      end do
c
c     go through all atoms
c
      do ia = 1, na
         if (skip(ia))  goto 130
c
c     skip to end of loop if buried atom
c
         if (abur(ia))  goto 130
c
c     special code for completely solvent-accessible atom
c
         if (afree(ia))  goto 120
c
c     gather convex edges for atom
c     clear number of convex edges for atom
c
         nepa = 0
c
c     pointer to first edge
c
         iep = afe(ia)
   10    continue
c
c     check whether finished gathering
c
         if (iep .le. 0)  goto 20
c
c     one more edge
c
         nepa = nepa + 1
         if (nepa .gt. maxepa) then
            call cerror ('Too many Convex Edges for Atom')
         end if
c
c      store vertices of edge
c
         av(1,nepa) = epv(1,iep)
         av(2,nepa) = epv(2,iep)
c
c     store convex edge number
c
         aep(nepa) = iep
         ic = epc(iep)
c
c     store circle number
c
         aic(nepa) = ic
c
c     get neighboring atom
c
         it = ct(ic)
         if (ta(1,it) .eq. ia) then
            ia2 = ta(2,it)
         else
            ia2 = ta(1,it)
         end if
c
c     store other atom number, we might need it sometime
c
         aia(nepa) = ia2
c
c     vector from atom to circle center; also
c     vector from atom to center of neighboring atom
c     sometimes we use one vector, sometimes the other
c
         do k = 1, 3
            acvect(k,nepa) = c(k,ic) - axyz(k,ia)
            aavect(k,nepa) = axyz(k,ia2) - axyz(k,ia)
         end do
c
c     circle radius
c
         acr(nepa) = cr(ic)
c
c     pointer to next edge
c
         iep = epnext(iep)
         goto 10
   20    continue
         if (nepa .le. 0) then
            call cerror ('No Edges for Non-buried, Non-free Atom')
         end if
c
c     form cycles; initialize all the
c     convex edges as not used in cycle
c
         do iepa = 1, nepa
            epused(iepa) = .false.
         end do
c
c     save old number of cycles
c
         ncyold = ncy
         nused = 0
         ncypa = 0
   30    continue
c
c     look for starting edge
c
         do iepa = 1, nepa
            if (.not. epused(iepa))  goto 40
         end do
c
c     cannot find starting edge, finished
c
         goto 80
   40    continue
c
c     pointer to edge
c
         iep = aep(iepa)
c
c     one edge so far for this cycle
c
         ncyep = 1
c
c     one more cycle for atom
c
         ncypa = ncypa + 1
         if (ncypa .gt. maxcypa) then
            call cerror ('Too many Cycles per Atom')
         end if
c
c     mark edge used in cycle
c
         epused(iepa) = .true.
         nused = nused + 1
c
c     one more cycle for molecule
c
         ncy = ncy + 1
         if (ncy .gt. maxcy)  call cerror ('Too many Cycles')
c
c     index of edge in atom cycle array
c
         cyepa(ncyep,ncypa) = iepa
c
c     store in molecule cycle array a pointer to edge
c
         cyep(ncyep,ncy) = iep
c
c     second vertex of this edge is the vertex to look
c     for next as the first vertex of another edge
c
         lookv = av(2,iepa)
c
c     if no vertex, this cycle is finished
c
         if (lookv .le. 0)  goto 70
   50    continue
c
c     look for next connected edge
c
         do jepa = 1, nepa
            if (epused(jepa))  goto 60
c
c     check second vertex of iepa versus first vertex of jepa
c
            if (av(1,jepa) .ne. lookv)  goto 60
c
c     edges are connected
c     pointer to edge
c
            iep = aep(jepa)
c
c     one more edge for this cycle
c
            ncyep = ncyep + 1
            if (ncyep .gt. mxcyep) then
               call cerror ('Too many Edges per Cycle')
            end if
            epused(jepa) = .true.
            nused = nused + 1
c
c     store index in local edge array
c
            cyepa(ncyep,ncypa) = jepa
c
c     store pointer to edge
c
            cyep(ncyep,ncy) = iep
c
c     new vertex to look for
c
            lookv = av(2,jepa)
c
c     if no vertex, this cycle is in trouble
c
            if (lookv .le. 0) then
               call cerror ('Pointer Error in Cycle')
            end if
            goto 50
   60       continue
         end do
c
c     it better connect to first edge of cycle
c
         if (lookv .ne. av(1,iepa)) then
            call cerror ('Cycle does not Close')
         end if
   70    continue
c
c     this cycle is finished, store number of edges in cycle
c
         ncyepa(ncypa) = ncyep
         cynep(ncy) = ncyep
         if (nused .ge. nepa)  goto 80
c
c     look for more cycles
c
         goto 30
   80    continue
c
c     compare cycles for inside/outside relation;
c     check to see if cycle i is inside cycle j
c
         do icya = 1, ncypa
            do jcya = 1, ncypa
               jcy = ncyold + jcya
               cycy(icya,jcya) = .true.
               if (icya .eq. jcya)  goto 90
c
c     if cycle j has two or fewer edges, nothing can
c     lie in its exterior; i is therefore inside j
c
               if (ncyepa(jcya) .le. 2)  goto 90
c
c     if cycles i and j have a pair of edges belonging
c     to the same circle, then they are outside each other
c
               do icyep = 1, ncyepa(icya)
                  iepa = cyepa(icyep,icya)
                  ic = aic(iepa)
                  do jcyep = 1, ncyepa(jcya)
                     jepa = cyepa(jcyep,jcya)
                     jc = aic(jepa)
                     if (ic .eq. jc) then
                        cycy(icya,jcya) = .false.
                        goto 90
                     end if
                  end do
               end do
               iepa = cyepa(1,icya)
               anaa = anorm(aavect(1,iepa))
               factor = ar(ia) / anaa
c
c     north pole and unit vector pointing south
c
               do k = 1, 3
                  pole(k) = factor*aavect(k,iepa) + axyz(k,ia)
                  unvect(k) = -aavect(k,iepa) / anaa
               end do
               cycy(icya,jcya) = ptincy(pole,unvect,jcy)
   90          continue
            end do
         end do
c
c     group cycles into faces; direct comparison for i and j
c
         do icya = 1, ncypa
            do jcya = 1, ncypa
c
c     tentatively say that cycles i and j bound
c     the same face if they are inside each other
c
               samef(icya,jcya) = (cycy(icya,jcya) .and.
     &                               cycy(jcya,icya))
            end do
         end do
c
c     if i is in exterior of k, and k is in interior of
c     i and j, then i and j do not bound the same face
c
         do icya = 1, ncypa
            do jcya = 1, ncypa
               if (icya .ne. jcya) then
                  do kcya = 1, ncypa
                     if (kcya.ne.icya .and. kcya.ne.jcya) then
                        if (cycy(kcya,icya) .and. cycy(kcya,jcya)
     &                        .and. .not.cycy(icya,kcya)) then
                           samef(icya,jcya) = .false.
                           samef(jcya,icya) = .false.
                        end if
                     end if
                  end do
               end if
            end do
         end do
c
c     fill gaps so that "samef" falls into complete blocks
c
         do icya = 1, ncypa-2
            do jcya = icya+1, ncypa-1
               if (samef(icya,jcya)) then
                  do kcya = jcya+1, ncypa
                     if (samef(jcya,kcya)) then
                        samef(icya,kcya) = .true.
                        samef(kcya,icya) = .true.
                     end if
                  end do
               end if
            end do
         end do
c
c     group cycles belonging to the same face
c
         do icya = 1, ncypa
            cyused(icya) = .false.
         end do
c
c     clear number of cycles used in bounding faces
c
         nused = 0
         do icya = 1, ncypa
c
c     check for already used
c
            if (cyused(icya))  goto 110
c
c     one more convex face
c
            nfp = nfp + 1
            if (nfp .gt. maxfp) then
               call cerror ('Too many Convex Faces')
            end if
c
c     clear number of cycles for face
c
            fpncy(nfp) = 0
c
c     pointer from face to atom
c
            fpa(nfp) = ia
c
c     look for all other cycles belonging to same face
c
            do jcya = 1, ncypa
c
c     check for cycle already used in another face
c
               if (cyused(jcya))  goto 100
c
c     cycles i and j belonging to same face
c
               if (.not. samef(icya,jcya))  goto 100
c
c     mark cycle used
c
               cyused(jcya) = .true.
               nused = nused + 1
c
c     one more cycle for face
c
               fpncy(nfp) = fpncy(nfp) + 1
               if (fpncy(nfp) .gt. mxfpcy) then
                  call cerror ('Too many Cycles bounding Convex Face')
               end if
               i = fpncy(nfp)
c
c     store cycle number
c
               fpcy(i,nfp) = ncyold + jcya
c
c     check for finished
c
               if (nused .ge. ncypa)  goto 130
  100          continue
            end do
  110       continue
         end do
c
c     should not fall through end of do loops
c
         call cerror ('Not all Cycles grouped into Convex Faces')
  120    continue
c
c     one face for free atom; no cycles
c
         nfp = nfp + 1
         if (nfp .gt. maxfp) then
            call cerror ('Too many Convex Faces')
         end if
         fpa(nfp) = ia
         fpncy(nfp) = 0
  130    continue
      end do
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine vam  --  volumes and areas of molecules  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "vam" takes the analytical molecular surface defined
c     as a collection of spherical and toroidal polygons
c     and uses it to compute the volume and surface area
c
c
      subroutine vam (volume,area)
      use sizes
      use faces
      use inform
      use iounit
      use math
      implicit none
      integer maxdot,maxop,nscale
      parameter (maxdot=1000)
      parameter (maxop=100)
      parameter (nscale=20)
      integer k,ke,ke2,kv
      integer ia,ic,ip,it
      integer ien,iep
      integer ifn,ifp,ifs
      integer iv,iv1,iv2
      integer isc,jfn
      integer ndots,idot
      integer nop,iop,nate
      integer neat,neatmx
      integer ivs(3)
      integer ispind(3)
      integer ispnd2(3)
      integer ifnop(maxop)
      integer, allocatable :: nlap(:)
      integer, allocatable :: enfs(:)
      integer, allocatable :: fnt(:,:)
      integer, allocatable :: nspt(:,:)
      real*8 volume,area
      real*8 alens,vint,vcone
      real*8 vpyr,vlens,hedron
      real*8 totap,totvp,totas
      real*8 totvs,totasp,totvsp
      real*8 totan,totvn
      real*8 alenst,alensn
      real*8 vlenst,vlensn,prism
      real*8 areap,volp,areas,vols
      real*8 areasp,volsp,arean,voln
      real*8 depth,triple,dist2
      real*8 areado,voldo,dot,dota
      real*8 ds2,dij2,dt,dpp
      real*8 rm,rat,rsc,rho
      real*8 sumsc,sumsig,sumlam
      real*8 stq,scinc,coran,corvn
      real*8 cenop(3,maxop)
      real*8 sdot(3),dotv(nscale)
      real*8 tau(3),ppm(3)
      real*8 xpnt1(3),xpnt2(3)
      real*8 qij(3),qji(3)
      real*8 vects(3,3)
      real*8 vect1(3),vect2(3)
      real*8 vect3(3),vect4(3)
      real*8 vect5(3),vect6(3)
      real*8 vect7(3),vect8(3)
      real*8 upp(3),thetaq(3)
      real*8 sigmaq(3)
      real*8 umq(3),upq(3)
      real*8 uc(3),uq(3),uij(3)
      real*8 dots(3,maxdot)
      real*8 tdots(3,maxdot)
      real*8, allocatable :: atmarea(:)
      real*8, allocatable :: depths(:)
      real*8, allocatable :: cora(:)
      real*8, allocatable :: corv(:)
      real*8, allocatable :: alts(:,:)
      real*8, allocatable :: fncen(:,:)
      real*8, allocatable :: fnvect(:,:,:)
      logical spindl
      logical alli,allj
      logical anyi,anyj
      logical case1,case2
      logical cinsp,cintp
      logical usenum
      logical vip(3)
      logical ate(maxop)
      logical, allocatable :: badav(:)
      logical, allocatable :: badt(:)
      logical, allocatable :: fcins(:,:)
      logical, allocatable :: fcint(:,:)
      logical, allocatable :: fntrev(:,:)
c
c
c     compute the volume of the interior polyhedron
c
      hedron = 0.0d0
      do ifn = 1, nfn
         call measpm (ifn,prism)
         hedron = hedron + prism
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (nlap(nfn))
      allocate (enfs(20*na))
      allocate (fnt(3,nfn))
      allocate (nspt(3,nfn))
      allocate (atmarea(na))
      allocate (depths(nfn))
      allocate (cora(nfn))
      allocate (corv(nfn))
      allocate (alts(3,nfn))
      allocate (fncen(3,nfn))
      allocate (fnvect(3,3,nfn))
      allocate (badav(nfn))
      allocate (badt(nfn))
      allocate (fcins(3,nfn))
      allocate (fcint(3,nfn))
      allocate (fntrev(3,nfn))
c
c     compute the area and volume due to convex faces
c     as well as the area partitioned among the atoms
c
      totap = 0.0d0
      totvp = 0.0d0
      do ia = 1, na
         atmarea(ia) = 0.0d0
      end do
      do ifp = 1, nfp
         call measfp (ifp,areap,volp)
         ia = fpa(ifp)
         atmarea(ia) = atmarea(ia) + areap
         totap = totap + areap
         totvp = totvp + volp
      end do
c
c     compute the area and volume due to saddle faces
c     as well as the spindle correction value
c
      totas = 0.0d0
      totvs = 0.0d0
      totasp = 0.0d0
      totvsp = 0.0d0
      do ifs = 1, nfs
         do k = 1, 2
            ien = fsen(k,ifs)
            if (ien .gt. 0)  enfs(ien) = ifs
         end do
         call measfs (ifs,areas,vols,areasp,volsp)
         totas = totas + areas
         totvs = totvs + vols
         totasp = totasp + areasp
         totvsp = totvsp + volsp
         if (areas-areasp .lt. 0.0d0) then
            call cerror ('Negative Area for Saddle Face')
         end if
      end do
c
c     compute the area and volume due to concave faces
c
      totan = 0.0d0
      totvn = 0.0d0
      do ifn = 1, nfn
         call measfn (ifn,arean,voln)
         totan = totan + arean
         totvn = totvn + voln
      end do
c
c     compute the area and volume lens correction values
c
      alenst = 0.0d0
      alensn = 0.0d0
      vlenst = 0.0d0
      vlensn = 0.0d0
      if (pr .le. 0.0d0)  goto 140
      ndots = maxdot
      call gendot (ndots,dots,pr,0.0d0,0.0d0,0.0d0)
      dota = (4.0d0 * pi * pr**2) / ndots
      do ifn = 1, nfn
         nlap(ifn) = 0
         cora(ifn) = 0.0d0
         corv(ifn) = 0.0d0
         badav(ifn) = .false.
         badt(ifn) = .false.
         do k = 1, 3
            nspt(k,ifn) = 0
         end do
         ien = fnen(1,ifn)
         iv = env(1,ien)
         ip = vp(iv)
         depths(ifn) = depth(ip,alts(1,ifn))
         do k = 1, 3
            fncen(k,ifn) = p(k,ip)
         end do
         ia = va(iv)
c
c     get vertices and vectors
c
         do ke = 1, 3
            ien = fnen(ke,ifn)
            ivs(ke) = env(1,ien)
            ia = va(ivs(ke))
            ifs = enfs(ien)
            iep = fsep(1,ifs)
            ic = epc(iep)
            it = ct(ic)
            fnt(ke,ifn) = it
            fntrev(ke,ifn) = (ta(1,it) .ne. ia)
         end do
         do ke = 1, 3
            do k = 1, 3
               vects(k,ke) = vxyz(k,ivs(ke)) - p(k,ip)
            end do
         end do
c
c     calculate normal vectors for the three planes
c     that cut out the geodesic triangle
c
         call vcross (vects(1,1),vects(1,2),fnvect(1,1,ifn))
         call vnorm (fnvect(1,1,ifn),fnvect(1,1,ifn))
         call vcross (vects(1,2),vects(1,3),fnvect(1,2,ifn))
         call vnorm (fnvect(1,2,ifn),fnvect(1,2,ifn))
         call vcross (vects(1,3),vects(1,1),fnvect(1,3,ifn))
         call vnorm (fnvect(1,3,ifn),fnvect(1,3,ifn))
      end do
      do ifn = 1, nfn-1
         do jfn = ifn+1, nfn
            dij2 = dist2(fncen(1,ifn),fncen(1,jfn))
            if (dij2 .gt. 4.0d0*pr**2)  goto 90
            if (depths(ifn).gt.pr .and. depths(jfn).gt.pr)  goto 90
c
c     these two probes may have intersecting surfaces
c
            dpp = sqrt(dist2(fncen(1,ifn),fncen(1,jfn)))
c
c     compute the midpoint
c
            do k = 1, 3
               ppm(k) = (fncen(k,ifn) + fncen(k,jfn)) / 2.0d0
               upp(k) = (fncen(k,jfn) - fncen(k,ifn)) / dpp
            end do
            rm = pr**2 - (dpp/2.0d0)**2
            if (rm .lt. 0.0d0)  rm = 0.0d0
            rm = sqrt(rm)
            rat = dpp / (2.0d0*pr)
            if (rat .gt. 1.0d0)  rat = 1.0d0
            if (rat .lt. -1.0d0)  rat = -1.0d0
            rho = asin(rat)
c
c     use circle-plane intersection routine
c
            alli = .true.
            anyi = .false.
            spindl = .false.
            do k = 1, 3
               ispind(k) = 0
               ispnd2(k) = 0
            end do
            do ke = 1, 3
               thetaq(ke) = 0.0d0
               sigmaq(ke) = 0.0d0
               tau(ke) = 0.0d0
               call cirpln (ppm,rm,upp,fncen(1,ifn),fnvect(1,ke,ifn),
     &                              cinsp,cintp,xpnt1,xpnt2)
               fcins(ke,ifn) = cinsp
               fcint(ke,ifn) = cintp
               if (.not. cinsp)  alli = .false.
               if (cintp)  anyi = .true.
               if (.not. cintp)  goto 10
               it = fnt(ke,ifn)
               if (tr(it) .gt. pr)  goto 10
               do ke2 = 1, 3
                  if (it .eq. fnt(ke2,jfn)) then
                     ispind(ke) = it
                     nspt(ke,ifn) = nspt(ke,ifn) + 1
                     ispnd2(ke2) = it
                     nspt(ke2,jfn) = nspt(ke2,jfn) + 1
                     spindl = .true.
                  end if
               end do
               if (ispind(ke) .eq. 0)  goto 10
c
c     check that the two ways of calculating
c     intersection points match
c
               rat = tr(it) / pr
               if (rat .gt. 1.0d0)  rat = 1.0d0
               if (rat .lt. -1.0d0)  rat = -1.0d0
               thetaq(ke) = acos(rat)
               stq = sin(thetaq(ke))
               if (fntrev(ke,ifn)) then
                  do k = 1, 3
                     uij(k) = -tax(k,it)
                  end do
               else
                  do k = 1, 3
                     uij(k) = tax(k,it)
                  end do
               end if
               do k = 1, 3
                  qij(k) = t(k,it) - stq * pr * uij(k)
                  qji(k) = t(k,it) + stq * pr * uij(k)
               end do
               do k = 1, 3
                  umq(k) = (qij(k) - ppm(k)) / rm
                  upq(k) = (qij(k) - fncen(k,ifn)) / pr
               end do
               call vcross (uij,upp,vect1)
               dt = dot(umq,vect1)
               if (dt .gt. 1.0d0)  dt = 1.0d0
               if (dt .lt. -1.0d0)  dt = -1.0d0
               sigmaq(ke) = acos(dt)
               call vcross (upq,fnvect(1,ke,ifn),vect1)
               call vnorm (vect1,uc)
               call vcross (upp,upq,vect1)
               call vnorm (vect1,uq)
               dt = dot(uc,uq)
               if (dt .gt. 1.0d0)  dt = 1.0d0
               if (dt .lt. -1.0d0)  dt = -1.0d0
               tau(ke) = pi - acos(dt)
   10          continue
            end do
            allj = .true.
            anyj = .false.
            do ke = 1, 3
               call cirpln (ppm,rm,upp,fncen(1,jfn),fnvect(1,ke,jfn),
     &                              cinsp,cintp,xpnt1,xpnt2)
               fcins(ke,jfn) = cinsp
               fcint(ke,jfn) = cintp
               if (.not. cinsp)  allj = .false.
               if (cintp)  anyj = .true.
            end do
            case1 = (alli .and. allj .and. .not.anyi .and. .not.anyj)
            case2 = (anyi .and. anyj .and. spindl)
            if (.not.case1 .and. .not.case2)  goto 90
c
c     this kind of overlap can be handled
c
            nlap(ifn) = nlap(ifn) + 1
            nlap(jfn) = nlap(jfn) + 1
            do ke = 1, 3
               ien = fnen(ke,ifn)
               iv1 = env(1,ien)
               iv2 = env(2,ien)
               do k = 1, 3
                  vect3(k) = vxyz(k,iv1) - fncen(k,ifn)
                  vect4(k) = vxyz(k,iv2) - fncen(k,ifn)
               end do
               do ke2 = 1, 3
                  if (ispind(ke) .eq. ispnd2(ke2))  goto 40
                  if (ispind(ke) .eq. 0)  goto 40
                  call cirpln (fncen(1,ifn),pr,fnvect(1,ke,ifn),
     &                           fncen(1,jfn),fnvect(1,ke2,jfn),
     &                           cinsp,cintp,xpnt1,xpnt2)
                  if (.not. cintp)  goto 40
                  ien = fnen(ke2,jfn)
                  iv1 = env(1,ien)
                  iv2 = env(2,ien)
                  do k = 1, 3
                     vect7(k) = vxyz(k,iv1) - fncen(k,jfn)
                     vect8(k) = vxyz(k,iv2) - fncen(k,jfn)
                  end do
c
c     check whether point lies on spindle arc
c
                  do k = 1, 3
                     vect1(k) = xpnt1(k) - fncen(k,ifn)
                     vect2(k) = xpnt2(k) - fncen(k,ifn)
                     vect5(k) = xpnt1(k) - fncen(k,jfn)
                     vect6(k) = xpnt2(k) - fncen(k,jfn)
                  end do
                  if (triple(vect3,vect1,fnvect(1,ke,ifn)) .lt. 0.0d0)
     &               goto 20
                  if (triple(vect1,vect4,fnvect(1,ke,ifn)) .lt. 0.0d0)
     &               goto 20
                  if (triple(vect7,vect5,fnvect(1,ke2,jfn)) .lt. 0.0d0)
     &               goto 20
                  if (triple(vect5,vect8,fnvect(1,ke2,jfn)) .lt. 0.0d0)
     &               goto 20
                  goto 30
   20             continue
                  if (triple(vect3,vect2,fnvect(1,ke,ifn)) .lt. 0.0d0)
     &               goto 40
                  if (triple(vect2,vect4,fnvect(1,ke,ifn)) .lt. 0.0d0)
     &               goto 40
                  if (triple(vect7,vect6,fnvect(1,ke2,jfn)) .lt. 0.0d0)
     &               goto 40
                  if (triple(vect6,vect8,fnvect(1,ke2,jfn)) .lt. 0.0d0)
     &               goto 40
   30             continue
                  badav(ifn) = .true.
   40             continue
               end do
            end do
            do ke = 1, 3
               ien = fnen(ke,ifn)
               iv1 = env(1,ien)
               iv2 = env(2,ien)
               do k = 1, 3
                  vect3(k) = vxyz(k,iv1) - fncen(k,ifn)
                  vect4(k) = vxyz(k,iv2) - fncen(k,ifn)
               end do
               do ke2 = 1, 3
                  if (ispind(ke) .eq. ispnd2(ke2))  goto 70
                  if (ispnd2(ke2) .eq. 0)  goto 70
                  call cirpln (fncen(1,jfn),pr,fnvect(1,ke2,jfn),
     &                           fncen(1,ifn),fnvect(1,ke,ifn),
     &                           cinsp,cintp,xpnt1,xpnt2)
                  if (.not. cintp)  goto 70
                  ien = fnen(ke2,jfn)
                  iv1 = env(1,ien)
                  iv2 = env(2,ien)
                  do k = 1, 3
                     vect7(k) = vxyz(k,iv1) - fncen(k,jfn)
                     vect8(k) = vxyz(k,iv2) - fncen(k,jfn)
                  end do
c
c     check whether point lies on spindle arc
c
                  do k = 1, 3
                     vect1(k) = xpnt1(k) - fncen(k,ifn)
                     vect2(k) = xpnt2(k) - fncen(k,ifn)
                     vect5(k) = xpnt1(k) - fncen(k,jfn)
                     vect6(k) = xpnt2(k) - fncen(k,jfn)
                  end do
                  if (triple(vect3,vect1,fnvect(1,ke,ifn)) .lt. 0.0d0)
     &               goto 50
                  if (triple(vect1,vect4,fnvect(1,ke,ifn)) .lt. 0.0d0)
     &               goto 50
                  if (triple(vect7,vect5,fnvect(1,ke2,jfn)) .lt. 0.0d0)
     &               goto 50
                  if (triple(vect5,vect8,fnvect(1,ke2,jfn)) .lt. 0.0d0)
     &               goto 50
                  goto 60
   50             continue
                  if (triple(vect3,vect2,fnvect(1,ke,ifn)) .lt. 0.0d0)
     &               goto 70
                  if (triple(vect2,vect4,fnvect(1,ke,ifn)) .lt. 0.0d0)
     &               goto 70
                  if (triple(vect7,vect6,fnvect(1,ke2,jfn)) .lt. 0.0d0)
     &               goto 70
                  if (triple(vect6,vect8,fnvect(1,ke2,jfn)) .lt. 0.0d0)
     &               goto 70
   60             continue
                  badav(jfn) = .true.
   70             continue
               end do
            end do
            sumlam = 0.0d0
            sumsig = 0.0d0
            sumsc = 0.0d0
            do ke = 1, 3
               if (ispind(ke) .ne. 0) then
                  sumlam = sumlam + pi - tau(ke)
                  sumsig = sumsig + sigmaq(ke) - pi
                  sumsc = sumsc + sin(sigmaq(ke))*cos(sigmaq(ke))
               end if
            end do
            alens = 2.0d0 * pr**2 * (pi - sumlam - sin(rho)*(pi+sumsig))
            vint = alens * pr / 3.0d0
            vcone = pr * rm**2 * sin(rho) * (pi+sumsig) / 3.0d0
            vpyr =  pr * rm**2 * sin(rho) * sumsc / 3.0d0
            vlens = vint - vcone + vpyr
            cora(ifn) = cora(ifn) + alens
            cora(jfn) = cora(jfn) + alens
            corv(ifn) = corv(ifn) + vlens
            corv(jfn) = corv(jfn) + vlens
c
c     check for vertex on opposing probe in face
c
            do kv = 1, 3
               vip(kv) = .false.
               ien = fnen(kv,jfn)
               iv = env(1,ien)
               do k = 1, 3
                  vect1(k) = vxyz(k,iv) - fncen(k,ifn)
               end do
               call vnorm (vect1,vect1)
               do ke = 1, 3
                  dt = dot(fnvect(1,ke,ifn),vxyz(1,iv))
                  if (dt .gt. 0.0d0)  goto 80
               end do
               vip(kv) = .true.
   80          continue
            end do
   90       continue
         end do
      end do
      do ifn = 1, nfn
         do ke = 1, 3
            if (nspt(ke,ifn) .gt. 1)  badt(ifn) = .true.
         end do
      end do
      do ifn = 1, nfn
         if (nlap(ifn) .le. 0)  goto 130
c
c     gather all overlapping probes
c
         nop = 0
         do jfn = 1, nfn
            if (ifn .ne. jfn) then
               dij2 = dist2(fncen(1,ifn),fncen(1,jfn))
               if (dij2 .le. 4.0d0*pr**2) then
                  if (depths(jfn) .le. pr) then
                     nop = nop + 1
                     if (nop .gt. maxop) then
                        call cerror ('NOP Overflow in VAM')
                     end if
                     ifnop(nop) = jfn
                     do k = 1, 3
                        cenop(k,nop) = fncen(k,jfn)
                     end do
                  end if
               end if
            end if
         end do
c
c     numerical calculation of the correction
c
         areado = 0.0d0
         voldo = 0.0d0
         scinc = 1.0d0 / nscale
         do isc = 1, nscale
            rsc = isc - 0.5d0
            dotv(isc) = pr * dota * rsc**2 * scinc**3
         end do
         do iop = 1, nop
            ate(iop) = .false.
         end do
         neatmx = 0
         do idot = 1, ndots
            do ke = 1, 3
               dt = dot(fnvect(1,ke,ifn),dots(1,idot))
               if (dt .gt. 0.0d0)  goto 120
            end do
            do k = 1, 3
               tdots(k,idot) = fncen(k,ifn) + dots(k,idot)
            end do
            do iop = 1, nop
               jfn = ifnop(iop)
               ds2 = dist2(tdots(1,idot),fncen(1,jfn))
               if (ds2 .lt. pr**2) then
                  areado = areado + dota
                  goto 100
               end if
            end do
  100       continue
            do isc = 1, nscale
               rsc = isc - 0.5d0
               do k = 1, 3
                  sdot(k) = fncen(k,ifn) + rsc*scinc*dots(k,idot)
               end do
               neat = 0
               do iop = 1, nop
                  jfn = ifnop(iop)
                  ds2 = dist2(sdot,fncen(1,jfn))
                  if (ds2 .lt. pr**2) then
                     do k = 1, 3
                        vect1(k) = sdot(k) - fncen(k,jfn)
                     end do
                     do ke = 1, 3
                        dt = dot(fnvect(1,ke,jfn),vect1)
                        if (dt .gt. 0.0d0)  goto 110
                     end do
                     neat = neat + 1
                     ate(iop) = .true.
  110                continue
                  end if
               end do
               if (neat .gt. neatmx)  neatmx = neat
               if (neat .gt. 0) then
                  voldo = voldo + dotv(isc) * (neat/(1.0d0+neat))
               end if
            end do
  120       continue
         end do
         coran = areado
         corvn = voldo
         nate = 0
         do iop = 1, nop
            if (ate(iop))  nate = nate + 1
         end do
c
c     use either the analytical or numerical correction
c
         usenum = (nate.gt.nlap(ifn) .or. neatmx.gt.1 .or. badt(ifn))
         if (usenum) then
            cora(ifn) = coran
            corv(ifn) = corvn
            alensn = alensn + cora(ifn)
            vlensn = vlensn + corv(ifn)
         else if (badav(ifn)) then
            corv(ifn) = corvn
            vlensn = vlensn + corv(ifn)
         end if
         alenst = alenst + cora(ifn)
         vlenst = vlenst + corv(ifn)
  130    continue
      end do
  140 continue
c
c     print out the decomposition of the area and volume
c
      if (debug) then
         write (iout,150)
  150    format (/,' Convex Surface Area for Individual Atoms :',/)
         k = 1
         do while (k .le. na)
            write (iout,160)  (ia,atmarea(ia),ia=k,min(k+4,na))
  160       format (1x,5(i7,f8.3))
            k = k + 5
         end do
         write (iout,170)
  170    format (/,' Surface Area and Volume by Geometry Type :')
         write (iout,180)  nfp,totap,totvp
  180    format (/,' Convex Faces :',i12,5x,'Area :',f13.3,
     &              4x,'Volume :',f13.3)
         write (iout,190)  nfs,totas,totvs
  190    format (' Saddle Faces :',i12,5x,'Area :',f13.3,
     &              4x,'Volume :',f13.3)
         write (iout,200)  nfn,totan,totvn
  200    format (' Concave Faces :',i11,5x,'Area :',f13.3,
     &              4x,'Volume :',f13.3)
         write (iout,210)  hedron
  210    format (' Buried Polyhedra :',36x,'Volume :',f13.3)
         if (totasp.ne.0.0d0 .or. totvsp.ne.0.0d0 .or.
     &       alenst.ne.0.0d0 .or. vlenst.ne.0.0d0) then
            write (iout,220)  -totasp,-totvsp
  220       format (/,' Spindle Correction :',11x,'Area :',f13.3,
     &                 4x,'Volume :',f13.3)
            write (iout,230)  -alenst-alensn,vlenst-vlensn
  230       format (' Lens Analytical Correction :',3x,'Area :',f13.3,
     &                 4x,'Volume :',f13.3)
         end if
         if (alensn.ne.0.0d0 .or. vlensn.ne.0.0d0) then
            write (iout,240)  alensn,vlensn
  240       format (' Lens Numerical Correction :',4x,'Area :',f13.3,
     &                 4x,'Volume :',f13.3)
         end if
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (nlap)
      deallocate (enfs)
      deallocate (fnt)
      deallocate (nspt)
      deallocate (atmarea)
      deallocate (depths)
      deallocate (cora)
      deallocate (corv)
      deallocate (alts)
      deallocate (fncen)
      deallocate (fnvect)
      deallocate (badav)
      deallocate (badt)
      deallocate (fcins)
      deallocate (fcint)
      deallocate (fntrev)
c
c     finally, compute the total area and total volume
c
      area = totap + totas + totan - totasp - alenst
      volume = totvp + totvs + totvn + hedron - totvsp + vlenst
      return
      end
c
c
c     ######################
c     ##                  ##
c     ##  function depth  ##
c     ##                  ##
c     ######################
c
c
      function depth (ip,alt)
      use faces
      implicit none
      integer k,ip,ia1,ia2,ia3
      real*8 depth,dot,alt(3)
      real*8 vect1(3),vect2(3)
      real*8 vect3(3),vect4(3)
c
c
      ia1 = pa(1,ip)
      ia2 = pa(2,ip)
      ia3 = pa(3,ip)
      do k = 1, 3
         vect1(k) = axyz(k,ia1) - axyz(k,ia3)
         vect2(k) = axyz(k,ia2) - axyz(k,ia3)
         vect3(k) = p(k,ip) - axyz(k,ia3)
      end do
      call vcross (vect1,vect2,vect4)
      call vnorm (vect4,vect4)
      depth = dot(vect4,vect3)
      do k = 1, 3
         alt(k) = vect4(k)
      end do
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine measpm  --  volume of interior polyhedron  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "measpm" computes the volume of a single prism section of
c     the full interior polyhedron
c
c
      subroutine measpm (ifn,prism)
      use faces
      implicit none
      integer k,ke,ien
      integer iv,ia,ip,ifn
      real*8 prism,height
      real*8 vect1(3)
      real*8 vect2(3)
      real*8 vect3(3)
      real*8 pav(3,3)
c
c
      height = 0.0d0
      do ke = 1, 3
         ien = fnen(ke,ifn)
         iv = env(1,ien)
         ia = va(iv)
         height = height + axyz(3,ia)
         ip = vp(iv)
         do k = 1, 3
            pav(k,ke) = axyz(k,ia) - p(k,ip)
         end do
      end do
      height = height / 3.0d0
      do k = 1, 3
         vect1(k) = pav(k,2) - pav(k,1)
         vect2(k) = pav(k,3) - pav(k,1)
      end do
      call vcross (vect1,vect2,vect3)
      prism = height * vect3(3) / 2.0d0
      return
      end
c
c
c     #########################
c     ##                     ##
c     ##  subroutine measfp  ##
c     ##                     ##
c     #########################
c
c
      subroutine measfp (ifp,areap,volp)
      use faces
      use math
      implicit none
      integer k,ke,ifp,iep
      integer ia,ia2,ic
      integer it,iv1,iv2
      integer ncycle,ieuler
      integer icyptr,icy
      integer nedge
      real*8 areap,volp
      real*8 dot,dt,gauss
      real*8 vecang,angle,geo
      real*8 pcurve,gcurve
      real*8 vect1(3)
      real*8 vect2(3)
      real*8 acvect(3)
      real*8 aavect(3)
      real*8 radial(3,mxcyep)
      real*8 tanv(3,2,mxcyep)
c
c
      ia = fpa(ifp)
      pcurve = 0.0d0
      gcurve = 0.0d0
      ncycle = fpncy(ifp)
      if (ncycle .gt. 0) then
         ieuler = 2 - ncycle
      else
         ieuler = 2
      end if
      do icyptr = 1, ncycle
         icy = fpcy(icyptr,ifp)
         nedge = cynep(icy)
         do ke = 1, nedge
            iep = cyep(ke,icy)
            ic = epc(iep)
            it = ct(ic)
            if (ia .eq. ta(1,it)) then
               ia2 = ta(2,it)
            else
               ia2 = ta(1,it)
            end if
            do k = 1, 3
               acvect(k) = c(k,ic) - axyz(k,ia)
               aavect(k) = axyz(k,ia2) - axyz(k,ia)
            end do
            call vnorm (aavect,aavect)
            dt = dot(acvect,aavect)
            geo = -dt / (ar(ia)*cr(ic))
            iv1 = epv(1,iep)
            iv2 = epv(2,iep)
            if (iv1.eq.0 .or. iv2.eq.0) then
               angle = 2.0d0 * pi
            else
               do k = 1, 3
                  vect1(k) = vxyz(k,iv1) - c(k,ic)
                  vect2(k) = vxyz(k,iv2) - c(k,ic)
                  radial(k,ke) = vxyz(k,iv1) - axyz(k,ia)
               end do
               call vnorm (radial(1,ke),radial(1,ke))
               call vcross (vect1,aavect,tanv(1,1,ke))
               call vnorm (tanv(1,1,ke),tanv(1,1,ke))
               call vcross (vect2,aavect,tanv(1,2,ke))
               call vnorm (tanv(1,2,ke),tanv(1,2,ke))
               angle = vecang(vect1,vect2,aavect,-1.0d0)
            end if
            gcurve = gcurve + cr(ic)*angle*geo
            if (nedge .ne. 1) then
               if (ke .gt. 1) then
                  angle = vecang(tanv(1,2,ke-1),tanv(1,1,ke),
     &                               radial(1,ke),1.0d0)
                  if (angle .lt. 0.0d0) then
                     call cerror ('Negative Angle in MEASFP')
                  end if
                  pcurve = pcurve + angle
               end if
            end if
         end do
         if (nedge .gt. 1) then
            angle = vecang(tanv(1,2,nedge),tanv(1,1,1),
     &                         radial(1,1),1.0d0)
            if (angle .lt. 0.0d0) then
               call cerror ('Negative Angle in MEASFP')
            end if
            pcurve = pcurve + angle
         end if
      end do
      gauss = 2.0d0*pi*ieuler - pcurve - gcurve
      areap = gauss * (ar(ia)**2)
      volp = areap * ar(ia) / 3.0d0
      return
      end
c
c
c     #########################
c     ##                     ##
c     ##  subroutine measfs  ##
c     ##                     ##
c     #########################
c
c
      subroutine measfs (ifs,areas,vols,areasp,volsp)
      use faces
      use math
      implicit none
      integer k,ifs,iep
      integer ic,ic1,ic2
      integer it,ia1,ia2
      integer iv1,iv2
      real*8 areas,vols
      real*8 areasp,volsp
      real*8 vecang,phi
      real*8 dot,d1,d2,w1,w2
      real*8 theta1,theta2
      real*8 rat,thetaq
      real*8 cone1,cone2
      real*8 term1,term2
      real*8 term3
      real*8 spin,volt
      real*8 vect1(3)
      real*8 vect2(3)
      real*8 aavect(3)
      logical cusp
c
c
      iep = fsep(1,ifs)
      ic = epc(iep)
      it = ct(ic)
      ia1 = ta(1,it)
      ia2 = ta(2,it)
      do k = 1, 3
         aavect(k) = axyz(k,ia2) - axyz(k,ia1)
      end do
      call vnorm (aavect,aavect)
      iv1 = epv(1,iep)
      iv2 = epv(2,iep)
      if (iv1.eq.0 .or. iv2.eq.0) then
         phi = 2.0d0 * pi
      else
         do k = 1, 3
            vect1(k) = vxyz(k,iv1) - c(k,ic)
            vect2(k) = vxyz(k,iv2) - c(k,ic)
         end do
         phi = vecang(vect1,vect2,aavect,1.0d0)
      end if
      do k = 1, 3
         vect1(k) = axyz(k,ia1) - t(k,it)
         vect2(k) = axyz(k,ia2) - t(k,it)
      end do
      d1 = -dot(vect1,aavect)
      d2 = dot(vect2,aavect)
      theta1 = atan2(d1,tr(it))
      theta2 = atan2(d2,tr(it))
c
c     check for cusps
c
      if (tr(it).lt.pr .and. theta1.gt.0.0d0
     &                 .and. theta2.gt.0.0d0) then
         cusp = .true.
         rat = tr(it) / pr
         if (rat .gt. 1.0d0)  rat = 1.0d0
         if (rat .lt. -1.0d0)  rat = -1.0d0
         thetaq = acos(rat)
      else
         cusp = .false.
         thetaq = 0.0d0
         areasp = 0.0d0
         volsp = 0.0d0
      end if
      term1 = tr(it) * pr * (theta1+theta2)
      term2 = (pr**2) * (sin(theta1) + sin(theta2))
      areas = phi * (term1-term2)
      if (cusp) then
         spin = tr(it)*pr*thetaq - pr**2 * sin(thetaq)
         areasp = 2.0d0 * phi * spin
      end if
c
      iep = fsep(1,ifs)
      ic2 = epc(iep)
      iep = fsep(2,ifs)
      ic1 = epc(iep)
      if (ca(ic1) .ne. ia1) then
         call cerror ('IA1 Inconsistency in MEASFS')
      end if
      do k = 1, 3
         vect1(k) = c(k,ic1) - axyz(k,ia1)
         vect2(k) = c(k,ic2) - axyz(k,ia2)
      end do
      w1 = dot(vect1,aavect)
      w2 = -dot(vect2,aavect)
      cone1 = phi * (w1*cr(ic1)**2)/6.0d0
      cone2 = phi * (w2*cr(ic2)**2)/6.0d0
      term1 = (tr(it)**2) * pr * (sin(theta1)+sin(theta2))
      term2 = sin(theta1)*cos(theta1) + theta1
     &           + sin(theta2)*cos(theta2) + theta2
      term2 = tr(it) * (pr**2) * term2
      term3 = sin(theta1)*cos(theta1)**2 + 2.0d0*sin(theta1)
     &           + sin(theta2)*cos(theta2)**2 + 2.0d0*sin(theta2)
      term3 = (pr**3 / 3.0d0) * term3
      volt = (phi/2.0d0) * (term1-term2+term3)
      vols = volt + cone1 + cone2
      if (cusp) then
         term1 = (tr(it)**2) * pr * sin(thetaq)
         term2 = sin(thetaq)*cos(thetaq) + thetaq
         term2 = tr(it) * (pr**2) * term2
         term3 = sin(thetaq)*cos(thetaq)**2 + 2.0d0*sin(thetaq)
         term3 = (pr**3 / 3.0d0) * term3
         volsp = phi * (term1-term2+term3)
      end if
      return
      end
c
c
c     #########################
c     ##                     ##
c     ##  subroutine measfn  ##
c     ##                     ##
c     #########################
c
c
      subroutine measfn (ifn,arean,voln)
      use faces
      use math
      implicit none
      integer k,ke,je
      integer ifn,ien
      integer iv,ia,ip
      real*8 arean,voln
      real*8 vecang,triple
      real*8 defect,simplx
      real*8 angle(3)
      real*8 pvv(3,3)
      real*8 pav(3,3)
      real*8 planev(3,3)
c
c
      do ke = 1, 3
         ien = fnen(ke,ifn)
         iv = env(1,ien)
         ia = va(iv)
         ip = vp(iv)
         do k = 1, 3
            pvv(k,ke) = vxyz(k,iv) - p(k,ip)
            pav(k,ke) = axyz(k,ia) - p(k,ip)
         end do
         if (pr .gt. 0.0d0)  call vnorm (pvv(1,ke),pvv(1,ke))
      end do
      if (pr .le. 0.0d0) then
         arean = 0.0d0
      else
         do ke = 1, 3
            je = ke + 1
            if (je .gt. 3)  je = 1
            call vcross (pvv(1,ke),pvv(1,je),planev(1,ke))
            call vnorm (planev(1,ke),planev(1,ke))
         end do
         do ke = 1, 3
            je = ke - 1
            if (je .lt. 1)  je = 3
            angle(ke) = vecang(planev(1,je),planev(1,ke),
     &                             pvv(1,ke),-1.0d0)
            if (angle(ke) .lt. 0.0d0) then
               call cerror ('Negative Angle in MEASFN')
            end if
         end do
         defect = 2.0d0*pi - (angle(1)+angle(2)+angle(3))
         arean = (pr**2) * defect
      end if
      simplx = -triple(pav(1,1),pav(1,2),pav(1,3)) / 6.0d0
      voln = simplx - arean*pr/3.0d0
      return
      end
c
c
c     #########################
c     ##                     ##
c     ##  subroutine projct  ##
c     ##                     ##
c     #########################
c
c
      subroutine projct (pnt,unvect,icy,ia,spv,nedge,fail)
      use faces
      implicit none
      integer k,ke,icy,ia
      integer nedge,iep,iv
      real*8 dot,dt,f
      real*8 polev(3)
      real*8 pnt(3)
      real*8 unvect(3)
      real*8 spv(3,*)
      logical fail
c
c
      fail = .false.
      nedge = cynep(icy)
      do ke = 1, cynep(icy)
c
c     vertex number (use first vertex of edge)
c
         iep = cyep(ke,icy)
         iv = epv(1,iep)
         if (iv .ne. 0) then
c
c     vector from north pole to vertex
c
            do k = 1, 3
               polev(k) = vxyz(k,iv) - pnt(k)
            end do
c
c     calculate multiplication factor
c
            dt = dot(polev,unvect)
            if (dt .eq. 0.0d0) then
               fail = .true.
               return
            end if
            f = (ar(ia)*2) / dt
            if (f .lt. 1.0d0) then
               fail = .true.
               return
            end if
c
c     projected vertex for this convex edge
c
            do k = 1, 3
               spv(k,ke) = pnt(k) + f*polev(k)
               continue
            end do
         end if
      end do
      return
      end
c
c
c     #######################
c     ##                   ##
c     ##  function ptincy  ##
c     ##                   ##
c     #######################
c
c
      function ptincy (pnt,unvect,icy)
      use faces
      implicit none
      integer k,ke,icy,iep
      integer ic,it,iatom
      integer iaoth,nedge
      real*8 dot,rotang
      real*8 totang
      real*8 unvect(3)
      real*8 pnt(3)
      real*8 acvect(3)
      real*8 cpvect(3)
      real*8 spv(3,mxcyep)
      real*8 epu(3,mxcyep)
      logical ptincy,fail
c
c
c     check for eaten by neighbor
c
      do ke = 1, cynep(icy)
         iep = cyep(ke,icy)
         ic = epc(iep)
         it = ct(ic)
         iatom = ca(ic)
         if (ta(1,it) .eq. iatom) then
            iaoth = ta(2,it)
         else
            iaoth = ta(1,it)
         end if
         do k = 1, 3
            acvect(k) = axyz(k,iaoth) - axyz(k,iatom)
            cpvect(k) = pnt(k) - c(k,ic)
         end do
         if (dot(acvect,cpvect) .ge. 0.0d0) then
            ptincy = .false.
            return
         end if
      end do
      if (cynep(icy) .le. 2) then
         ptincy = .true.
         return
      end if
      call projct (pnt,unvect,icy,iatom,spv,nedge,fail)
      if (fail) then
         ptincy = .true.
         return
      end if
      call epuclc (spv,nedge,epu)
      totang = rotang(epu,nedge,unvect)
      ptincy = (totang .gt. 0.0d0)
      return
      end
c
c
c     #########################
c     ##                     ##
c     ##  subroutine epuclc  ##
c     ##                     ##
c     #########################
c
c
      subroutine epuclc (spv,nedge,epu)
      implicit none
      integer k,ke,ke2,le
      integer nedge
      real*8 anorm,epun
      real*8 spv(3,*)
      real*8 epu(3,*)
c
c
c     calculate unit vectors along edges
c
      do ke = 1, nedge
c
c     get index of second edge of corner
c
         if (ke .lt. nedge) then
            ke2 = ke + 1
         else
            ke2 = 1
         end if
c
c     unit vector along edge of cycle
c
         do k = 1, 3
            epu(k,ke) = spv(k,ke2) - spv(k,ke)
         end do
         epun = anorm(epu(1,ke))
c        if (epun .le. 0.0d0)  call cerror ('Null Edge in Cycle')
c
c     normalize
c
         if (epun .gt. 0.0d0) then
            do k = 1, 3
               epu(k,ke) = epu(k,ke) / epun
            end do
         else
            do k = 1, 3
               epu(k,ke) = 0.0d0
            end do
         end if
      end do
c
c     vectors for null edges come from following or preceding edges
c
      do ke = 1, nedge
         if (anorm(epu(1,ke)) .le. 0.0d0) then
            le = ke - 1
            if (le .le. 0)  le = nedge
            do k = 1, 3
               epu(k,ke) = epu(k,le)
            end do
         end if
      end do
      return
      end
c
c
c     #######################
c     ##                   ##
c     ##  function rotang  ##
c     ##                   ##
c     #######################
c
c
      function rotang (epu,nedge,unvect)
      implicit none
      integer ke,nedge
      real*8 rotang,totang
      real*8 dot,dt,ang
      real*8 unvect(3)
      real*8 crs(3)
      real*8 epu(3,*)
c
c
      totang = 0.0d0
c
c     sum angles at vertices of cycle
c
      do ke = 1, nedge
         if (ke .lt. nedge) then
            dt = dot(epu(1,ke),epu(1,ke+1))
            call vcross (epu(1,ke),epu(1,ke+1),crs)
         else
c
c     closing edge of cycle
c
            dt = dot(epu(1,ke),epu(1,1))
            call vcross (epu(1,ke),epu(1,1),crs)
         end if
         if (dt .lt. -1.0d0)  dt = -1.0d0
         if (dt .gt. 1.0d0)  dt = 1.0d0
         ang = acos(dt)
         if (dot(crs,unvect) .gt. 0.0d0)  ang = -ang
c
c     add to total for cycle
c
         totang = totang + ang
      end do
      rotang = totang
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine vcross  --  find cross product of two vectors  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "vcross" finds the cross product of two vectors
c
c
      subroutine vcross (x,y,z)
      implicit none
      real*8 x(3),y(3),z(3)
c
c
      z(1) = x(2)*y(3) - x(3)*y(2)
      z(2) = x(3)*y(1) - x(1)*y(3)
      z(3) = x(1)*y(2) - x(2)*y(1)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  function dot  --  find the dot product of two vectors  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "dot" finds the dot product of two vectors
c
c
      function dot (x,y)
      implicit none
      real*8 dot,x(3),y(3)
c
c
      dot = x(1)*y(1) + x(2)*y(2) + x(3)*y(3)
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  function anorm  --  find the length of a vector  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "anorm" finds the norm (length) of a vector; used as a
c     service routine by the Connolly surface area and volume
c     computation
c
c
      function anorm (x)
      implicit none
      real*8 anorm,x(3)
c
c
      anorm = x(1)**2 + x(2)**2 + x(3)**2
      if (anorm .lt. 0.0d0)  anorm = 0.0d0
      anorm = sqrt(anorm)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine vnorm  --  normalize a vector to unit length  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "vnorm" normalizes a vector to unit length; used as a
c     service routine by the Connolly surface area and volume
c     computation
c
c
      subroutine vnorm (x,xn)
      implicit none
      integer k
      real*8 ax,anorm
      real*8 x(3),xn(3)
c
c
      ax = anorm(x)
      do k = 1, 3
         xn(k) = x(k) / ax
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  function dist2  --  distance squared between two points  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "dist2" finds the distance squared between two points; used
c     as a service routine by the Connolly surface area and volume
c     computation
c
c
      function dist2 (x,y)
      implicit none
      real*8 dist2
      real*8 x(3),y(3)
c
c
      dist2 = (x(1)-y(1))**2 + (x(2)-y(2))**2 + (x(3)-y(3))**2
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  function triple  --  form triple product of three vectors  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "triple" finds the triple product of three vectors; used as
c     a service routine by the Connolly surface area and volume
c     computation
c
c
      function triple (x,y,z)
      implicit none
      real*8 triple,dot
      real*8 x(3),y(3)
      real*8 z(3),xy(3)
c
c
      call vcross (x,y,xy)
      triple = dot(xy,z)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  function vecang  --  finds the angle between two vectors  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "vecang" finds the angle between two vectors handed with respect
c     to a coordinate axis; returns an angle in the range [0,2*pi]
c
c
      function vecang (v1,v2,axis,hand)
      use math
      implicit none
      real*8 vecang,hand
      real*8 angle,dt
      real*8 a1,a2,a12
      real*8 anorm,dot
      real*8 triple
      real*8 v1(3),v2(3)
      real*8 axis(3)
c
c
      a1 = anorm(v1)
      a2 = anorm(v2)
      dt = dot(v1,v2)
      a12 = a1 * a2
      if (abs(a12) .ne. 0.0d0)  dt = dt/a12
      if (dt .lt. -1.0d0)  dt = -1.0d0
      if (dt .gt. 1.0d0)  dt = 1.0d0
      angle = acos(dt)
      if (hand*triple(v1,v2,axis) .lt. 0.0d0) then
         vecang = 2.0d0*pi - angle
      else
         vecang = angle
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine cirpln  --  locate circle-plane intersection  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "cirpln" determines the points of intersection between a
c     specified circle and plane
c
c
      subroutine cirpln (circen,cirrad,cirvec,plncen,plnvec,
     &                        cinsp,cintp,xpnt1,xpnt2)
      implicit none
      integer k
      real*8 anorm,dot
      real*8 dcp,dir
      real*8 ratio,rlen
      real*8 cirrad
      real*8 circen(3)
      real*8 cirvec(3)
      real*8 plncen(3)
      real*8 plnvec(3)
      real*8 xpnt1(3)
      real*8 xpnt2(3)
      real*8 cpvect(3)
      real*8 pnt1(3)
      real*8 vect1(3)
      real*8 vect2(3)
      real*8 uvect1(3)
      real*8 uvect2(3)
      logical cinsp
      logical cintp
c
c
      do k = 1, 3
         cpvect(k) = plncen(k) - circen(k)
      end do
      dcp = dot(cpvect,plnvec)
      cinsp = (dcp .gt. 0.0d0)
      call vcross (plnvec,cirvec,vect1)
      if (anorm(vect1) .gt. 0.0d0) then
         call vnorm (vect1,uvect1)
         call vcross (cirvec,uvect1,vect2)
         if (anorm(vect2) .gt. 0.0d0) then
            call vnorm (vect2,uvect2)
            dir = dot(uvect2,plnvec)
            if (dir .ne. 0.0d0) then
               ratio = dcp / dir
               if (abs(ratio) .le. cirrad) then
                  do k = 1, 3
                     pnt1(k) = circen(k) + ratio*uvect2(k)
                  end do
                  rlen = cirrad**2 - ratio**2
                  if (rlen .lt. 0.0d0)  rlen = 0.0d0
                  rlen = sqrt(rlen)
                  do k = 1, 3
                     xpnt1(k) = pnt1(k) - rlen*uvect1(k)
                     xpnt2(k) = pnt1(k) + rlen*uvect1(k)
                  end do
                  cintp = .true.
                  return
               end if
            end if
         end if
      end if
      cintp = .false.
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine gendot  --  find surface points on unit sphere  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "gendot" finds the coordinates of a specified number of surface
c     points for a sphere with the input radius and coordinate center
c
c
      subroutine gendot (ndots,dots,radius,xcenter,ycenter,zcenter)
      use math
      implicit none
      integer i,j,k
      integer ndots
      integer nequat
      integer nvert
      integer nhoriz
      real*8 fi,fj
      real*8 x,y,z,xy
      real*8 xcenter
      real*8 ycenter
      real*8 zcenter
      real*8 radius
      real*8 dots(3,*)
c
c
      nequat = int(sqrt(pi*dble(ndots)))
      nvert = int(0.5d0*dble(nequat))
      if (nvert .lt. 1)  nvert = 1
      k = 0
      do i = 0, nvert
         fi = (pi * dble(i)) / dble(nvert)
         z = cos(fi)
         xy = sin(fi)
         nhoriz = int(dble(nequat)*xy)
         if (nhoriz .lt. 1)  nhoriz = 1
         do j = 0, nhoriz-1
            fj = (2.0d0 * pi * dble(j)) / dble(nhoriz)
            x = cos(fj) * xy
            y = sin(fj) * xy
            k = k + 1
            dots(1,k) = x*radius + xcenter
            dots(2,k) = y*radius + ycenter
            dots(3,k) = z*radius + zcenter
            if (k .ge. ndots)  goto 10
         end do
      end do
   10 continue
      ndots = k
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine cerror  --  surface area-volume error message  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "cerror" is the error handling routine for the Connolly
c     surface area and volume computation
c
c
      subroutine cerror (string)
      use iounit
      implicit none
      integer leng,trimtext
      character*(*) string
c
c
c     write out the error message and quit
c
      leng = trimtext (string)
      write (iout,10)  string(1:leng)
   10 format (/,' CONNOLLY  --  ',a)
      call fatal
      end
