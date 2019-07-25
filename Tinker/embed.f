c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine embed  --  structures via distance geometry  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "embed" is a distance geometry routine patterned after the
c     ideas of Gordon Crippen, Irwin Kuntz and Tim Havel; it takes
c     as input a set of upper and lower bounds on the interpoint
c     distances, chirality restraints and torsional restraints,
c     and attempts to generate a set of coordinates that satisfy
c     the input bounds and restraints
c
c     literature references:
c
c     G. M. Crippen and T. F. Havel, "Distance Geometry and Molecular
c     Conformation", Research Studies Press, Letchworth U.K., 1988,
c     John Wiley and Sons, U.S. distributor
c
c     T. F. Havel, "An Evaluation of Computational Strategies for
c     Use in the Determination of Protein Structure from Distance
c     Constraints obtained by Nuclear Magnetic Resonance", Progress
c     in Biophysics and Molecular Biology, 56, 43-78 (1991)
c
c
      subroutine embed
      use sizes
      use atoms
      use disgeo
      use files
      use inform
      use iounit
      use minima
      use refer
      implicit none
      integer maxeigen
      parameter (maxeigen=5)
      integer i,j,nvar,nstep
      integer lext,igeo,freeunit
      integer maxneg,nneg
      integer maxinner,ninner
      integer maxouter,nouter
      real*8 fctval,grdmin
      real*8 dt,wall,cpu
      real*8 temp_start
      real*8 temp_stop
      real*8 rg,rmsorig
      real*8 rmsflip,mass
      real*8 bounds,contact
      real*8 chiral,torsion
      real*8 local,locerr
      real*8 bnderr,vdwerr
      real*8 chirer,torser
      real*8 evl(maxeigen)
      real*8, allocatable :: v(:)
      real*8, allocatable :: a(:)
      real*8, allocatable :: evc(:,:)
      real*8, allocatable :: matrix(:,:)
      real*8, allocatable :: derivs(:,:)
      logical done,valid
      logical exist,info
      character*7 errtyp,ext
      character*240 title
      character*240 geofile
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (evc(n,maxeigen))
      allocate (matrix(n,n))
c
c     initialize any chirality restraints, then smooth the
c     bounds via triangle and inverse triangle inequalities;
c     currently these functions are performed by "distgeom"
c
c     call kchiral
c     if (verbose .and. n.le.130) then
c        title = 'Input Distance Bounds :'
c        call grafic (n,bnd,title)
c     end if
c     call geodesic
c     if (verbose .and. n.le.130)) then
c        title = 'Triangle Smoothed Bounds :'
c        call grafic (n,bnd,title)
c     end if
c
c     generate a distance matrix between the upper and
c     lower bounds, then convert to a metric matrix
c
      maxinner = 3
      maxouter = 3
      maxneg = 2
      nouter = 0
      valid = .false.
      do while (.not. valid)
         ninner = 0
         done = .false.
         do while (.not. done)
            ninner = ninner + 1
            call dstmat (matrix)
            call metric (matrix,nneg)
            if (nneg.le.maxneg .or. ninner.eq.maxinner)  done = .true.
            compact = 0.0d0
         end do
         if (verbose .and. nneg.gt.maxneg) then
            write (iout,10)  nneg
   10       format (/,' EMBED  --  Warning, Using Metric Matrix',
     &                  ' with',i4,' Negative Distances')
         end if
c
c     find the principle components of metric matrix, then
c     generate the trial Cartesian coordinates
c
         nouter = nouter + 1
         call eigen (evl,evc,matrix,valid)
         call coords (evl,evc)
         if (nouter.eq.maxouter .and. .not.valid) then
            valid = .true.
            if (verbose) then
               write (iout,20)
   20          format (/,' EMBED  --  Warning, Using Poor Initial',
     &                    ' Coordinates')
            end if
         end if
      end do
c
c     superimpose embedded structure and enantiomer on reference
c
      info = verbose
      verbose = .false.
      call impose (nref(1),xref,yref,zref,n,x,y,z,rmsorig)
      if (use_invert) then
         do i = 1, n
            x(i) = -x(i)
         end do
         call impose (nref(1),xref,yref,zref,n,x,y,z,rmsflip)
         if (rmsorig .lt. rmsflip) then
            do i = 1, n
               x(i) = -x(i)
            end do
            call impose (nref(1),xref,yref,zref,n,x,y,z,rmsorig)
         end if
         write (iout,30)  rmsorig,rmsflip
   30    format (/,' RMS Superposition for Original and',
     &              ' Enantiomer : ',2f12.4)
      end if
      verbose = info
c
c     compute an index of compaction for the embedded structure
c
      call chksize
c
c     write out the unrefined embedded atomic coordinate set
c
      if (debug) then
         i = 0
         exist = .true.
         do while (exist)
            i = i + 1
            lext = 3
            call numeral (i,ext,lext)
            geofile = filename(1:leng)//'-embed'//'.'//ext(1:lext)
            inquire (file=geofile,exist=exist)
         end do
         igeo = freeunit ()
         open (unit=igeo,file=geofile,status='new')
         call prtxyz (igeo)
         close (unit=igeo)
         title = 'after Embedding :'
         call fracdist (title)
      end if
c
c     use majorization to improve initial embedded coordinates
c
      do i = 1, n
         matrix(i,i) = 0.0d0
         do j = i+1, n
            matrix(j,i) = matrix(i,j)
         end do
      end do
      call majorize (matrix)
c
c     square the bounds for use during structure refinement
c
      do i = 1, n
         do j = 1, n
            dbnd(j,i) = dbnd(j,i)**2
         end do
      end do
c
c     perform dynamic allocation of some local arrays
c
      if (use_anneal) then
         nvar = 3 * n
         allocate (v(nvar))
         allocate (a(nvar))
      end if
c
c     minimize the error function via simulated annealing
c
      if (verbose)  call settime
      if (use_anneal) then
         iprint = 0
         if (verbose)  iprint = 10
         grdmin = 1.0d0
         mass = 10000.0d0
         do i = 1, nvar
            v(i) = 0.0d0
            a(i) = 0.0d0
         end do
         errtyp = 'FINAL'
         call refine (errtyp,fctval,grdmin)
         nstep = 1000
         dt = 0.04d0
         temp_start = 200.0d0
         temp_stop = 200.0d0
         call explore (errtyp,nstep,dt,mass,temp_start,temp_stop,v,a)
         nstep = 10000
         dt = 0.2d0
         temp_start = 200.0d0
         temp_stop = 0.0d0
         call explore (errtyp,nstep,dt,mass,temp_start,temp_stop,v,a)
         grdmin = 0.01d0
         call refine (errtyp,fctval,grdmin)
c
c     minimize the error function via nonlinear optimization
c
      else
         iprint = 0
         if (verbose)  iprint = 10
         grdmin = 0.01d0
         errtyp = 'INITIAL'
         call refine (errtyp,fctval,grdmin)
         errtyp = 'MIDDLE'
         call refine (errtyp,fctval,grdmin)
         errtyp = 'FINAL'
         call refine (errtyp,fctval,grdmin)
      end if
      if (verbose) then
         call gettime (wall,cpu)
         write (iout,40)  wall
   40    format (/,' Time Required for Refinement :',10x,f12.2,
     &              ' seconds')
      end if
c
c     perform deallocation of some local arrays
c
      if (use_anneal) then
         deallocate (v)
         deallocate (a)
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (derivs(3,n))
c
c     print the final error function and its components
c
      bounds = bnderr (derivs)
      contact = vdwerr (derivs)
      local = locerr (derivs)
      chiral = chirer (derivs)
      torsion = torser (derivs)
      write (iout,50)  fctval,bounds,contact,local,chiral,torsion
   50 format (/,' Results of Distance Geometry Protocol :',
     &        //,' Final Error Function Value :',10x,f16.4,
     &        //,' Distance Restraint Error :',12x,f16.4,
     &        /,' Hard Sphere Contact Error :',11x,f16.4,
     &        /,' Local Geometry Error :',16x,f16.4,
     &        /,' Chirality-Planarity Error :',11x,f16.4,
     &        /,' Torsional Restraint Error :',11x,f16.4)
c
c     take the root of the currently squared distance bounds
c
      do i = 1, n
         do j = 1, n
            dbnd(j,i) = sqrt(dbnd(j,i))
         end do
      end do
c
c     print the final rms deviations and radius of gyration
c
      title = 'after Refinement :'
      call rmserror (title)
      call gyrate (rg)
      write (iout,60)  rg
   60 format (/,' Radius of Gyration after Refinement  :',6x,f16.4)
      if (verbose .and. n.le.130)  call dmdump (matrix)
c
c     print the normalized fractional distance distribution
c
      if (debug) then
         title = 'after Refinement :'
         call fracdist (title)
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (evc)
      deallocate (matrix)
      deallocate (derivs)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine kchiral  --  chirality restraint assignment  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "kchiral" determines the target value for each chirality
c     and planarity restraint as the signed volume of the
c     parallelpiped spanned by vectors from a common atom to
c     each of three other atoms
c
c
      subroutine kchiral
      use sizes
      use atoms
      use inform
      use iounit
      use restrn
      implicit none
      integer i,j,ia,ib,ic,id
      real*8 xad,yad,zad
      real*8 xbd,ybd,zbd
      real*8 xcd,ycd,zcd
      real*8 c1,c2,c3
c
c
c     compute the signed volume of each parallelpiped;
c     if the defining atoms almost lie in a plane, then
c     set the signed volume to exactly zero
c
      do i = 1, nchir
         ia = ichir(1,i)
         ib = ichir(2,i)
         ic = ichir(3,i)
         id = ichir(4,i)
         xad = x(ia) - x(id)
         yad = y(ia) - y(id)
         zad = z(ia) - z(id)
         xbd = x(ib) - x(id)
         ybd = y(ib) - y(id)
         zbd = z(ib) - z(id)
         xcd = x(ic) - x(id)
         ycd = y(ic) - y(id)
         zcd = z(ic) - z(id)
         c1 = ybd*zcd - zbd*ycd
         c2 = ycd*zad - zcd*yad
         c3 = yad*zbd - zad*ybd
         chir(1,i) = 0.1d0
         chir(2,i) = xad*c1 + xbd*c2 + xcd*c3
         if (abs(chir(2,i)) .lt. 1.0d0)  chir(2,i) = 0.0d0
         chir(3,i) = chir(2,i)
      end do
c
c     print out the results for each restraint
c
      if (verbose) then
         if (nchir .ne. 0) then
            write (iout,10)
   10       format (/,' Chirality and Planarity Constraints :')
            write (iout,20)
   20       format (/,18x,'Atom Numbers',12x,'Signed Volume',/)
         end if
         do i = 1, nchir
            write (iout,30)  i,(ichir(j,i),j=1,4),chir(2,i)
   30       format (i6,5x,4i6,5x,f12.4)
         end do
      end if
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine triangle  --  triangle inequality smoothing  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "triangle" smooths the upper and lower distance bounds via
c     the triangle inequality using a full-matrix variant of the
c     Floyd-Warshall shortest path algorithm; this routine is
c     usually much slower than the sparse matrix shortest path
c     methods in "geodesic" and "trifix", and should be used only
c     for comparison with answers generated by those routines
c
c     literature reference:
c
c     A. W. M. Dress and T. F. Havel, "Shortest-Path Problems and
c     Molecular Conformation", Discrete Applied Mathematics, 19,
c     129-144 (1988)
c
c
      subroutine triangle
      use sizes
      use atoms
      use disgeo
      use iounit
      implicit none
      integer i,j,k
      integer ik1,ik2
      integer jk1,jk2
      real*8 eps
      real*8 lij,lik,ljk
      real*8 uij,uik,ujk
c
c
c     use full-matrix algorithm to smooth upper and lower bounds
c
      eps = 1.0d-10
      do k = 1, n
         do i = 1, n-1
            ik1 = min(i,k)
            ik2 = max(i,k)
            lik = dbnd(ik2,ik1)
            uik = dbnd(ik1,ik2)
            do j = i+1, n
               lij = dbnd(j,i)
               uij = dbnd(i,j)
               jk1 = min(j,k)
               jk2 = max(j,k)
               ljk = dbnd(jk2,jk1)
               ujk = dbnd(jk1,jk2)
               lij = max(lij,lik-ujk,ljk-uik)
               uij = min(uij,uik+ujk)
               if (lij-uij .gt. eps) then
                  write (iout,10)
   10             format (/,' TRIANGLE  --  Inconsistent Bounds;',
     &                       ' Geometrically Impossible')
                  write (iout,20)  i,j,lij,uij
   20             format (/,' Error at :',6x,2i6,3x,2f9.4)
                  write (iout,30)  i,k,lik,uik,j,k,ljk,ujk
   30             format (/,' Traced to :',5x,2i6,3x,2f9.4,
     &                    /,17x,2i6,3x,2f9.4)
                  call fatal
               end if
               if (lij-dbnd(j,i) .gt. eps) then
                  write (iout,40)  i,j,dbnd(j,i),lij
   40             format (' TRIANGLE  --  Altered Lower Bound at',
     &                       2x,2i6,3x,f9.4,' -->',f9.4)
               end if
               if (dbnd(i,j)-uij .gt. eps) then
                  write (iout,50)  i,j,dbnd(i,j),uij
   50             format (' TRIANGLE  --  Altered Upper Bound at',
     &                       2x,2i6,3x,f9.4,' -->',f9.4)
               end if
               dbnd(j,i) = lij
               dbnd(i,j) = uij
            end do
         end do
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine geodesic  --  sparse matrix triangle smoothing  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "geodesic" smooths the upper and lower distance bounds via
c     the triangle inequality using a sparse matrix version of a
c     shortest path algorithm
c
c     literature reference:
c
c     G. M. Crippen and T. F. Havel, "Distance Geometry and Molecular
c     Conformation", Research Studies Press, Letchworth U.K., 1988,
c     John Wiley and Sons, U.S. distributor, see section 6-2
c
c
      subroutine geodesic
      use sizes
      use atoms
      use disgeo
      use restrn
      implicit none
      integer i,j,k,nlist
      integer, allocatable :: list(:)
      integer, allocatable :: key(:)
      integer, allocatable :: start(:)
      integer, allocatable :: stop(:)
      real*8, allocatable :: upper(:)
      real*8, allocatable :: lower(:)
c
c
c     perform dynamic allocation of some local arrays
c
      nlist = 2 * ndfix
      allocate (list(nlist))
      allocate (key(nlist))
      allocate (start(n))
      allocate (stop(n))
      allocate (upper(n))
      allocate (lower(n))
c
c     build an indexed list of atoms in distance restraints
c
      do i = 1, n
         start(i) = 0
         stop(i) = -1
      end do
      do i = 1, ndfix
         list(i) = idfix(1,i)
         list(i+ndfix) = idfix(2,i)
      end do
      call sort3 (nlist,list,key)
      j = -1
      do i = 1, nlist
         k = list(i)
         if (k .ne. j) then
            start(k) = i
            j = k
         end if
      end do
      j = -1
      do i = nlist, 1, -1
         k = list(i)
         if (k .ne. j) then
            stop(k) = i
            j = k
         end if
      end do
      do i = 1, nlist
         k = key(i)
         if (k .le. ndfix) then
            list(i) = idfix(2,k)
         else
            list(i) = idfix(1,k-ndfix)
         end if
      end do
c
c     triangle smooth bounds via sparse shortest path method
c
      do i = 1, n
         call minpath (i,upper,lower,start,stop,list)
         do j = i+1, n
            dbnd(i,j) = upper(j)
            dbnd(j,i) = max(lower(j),dbnd(j,i))
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (list)
      deallocate (key)
      deallocate (start)
      deallocate (stop)
      deallocate (upper)
      deallocate (lower)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine minpath  --  triangle smoothed bounds to atom  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "minpath" is a routine for finding the triangle smoothed upper
c     and lower bounds of each atom to a specified root atom using a
c     sparse variant of the Bellman-Ford shortest path algorithm
c
c     literature reference:
c
c     D. P. Bertsekas, "A Simple and Fast Label Correcting Algorithm
c     for Shortest Paths", Networks, 23, 703-709 (1993)
c
c
      subroutine minpath (root,upper,lower,start,stop,list)
      use sizes
      use atoms
      use couple
      use disgeo
      implicit none
      integer i,j,k
      integer narc,root
      integer head,tail
      integer, allocatable :: iarc(:)
      integer, allocatable :: queue(:)
      integer start(*)
      integer stop(*)
      integer list(*)
      real*8 big,small
      real*8 upper(*)
      real*8 lower(*)
      logical enter
      logical, allocatable :: queued(:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (iarc(n))
      allocate (queue(n))
      allocate (queued(n))
c
c     initialize candidate atom queue and the path lengths
c
      do i = 1, n
         queued(i) = .false.
         upper(i) = 1000000.0d0
         lower(i) = 0.0d0
      end do
c
c     put the root atom into the queue of candidate atoms
c
      head = root
      tail = root
      queue(root) = 0
      queued(root) = .true.
      upper(root) = 0.0d0
c
c     get the next candidate atom from head of queue
c
      do while (head .ne. 0)
         j = head
         queued(j) = .false.
         head = queue(head)
c
c     make a list of arcs to the current candidate atom
c
         narc = 0
         do i = 1, n12(j)
            k = i12(i,j)
            if (k .ne. root) then
               narc = narc + 1
               iarc(narc) = k
            end if
         end do
         do i = 1, n13(j)
            k = i13(i,j)
            if (k .ne. root) then
               narc = narc + 1
               iarc(narc) = k
            end if
         end do
         do i = 1, n14(j)
            k = i14(i,j)
            if (k .ne. root) then
               narc = narc + 1
               iarc(narc) = k
            end if
         end do
         do i = start(j), stop(j)
            k = list(i)
            if (k .ne. root) then
               narc = narc + 1
               iarc(narc) = k
            end if
         end do
c
c     check each arc for alteration of the path length bounds
c
         do i = 1, narc
            k = iarc(i)
            if (k .lt. j) then
               big = upper(j) + dbnd(k,j)
               small = max(dbnd(j,k)-upper(j),lower(j)-dbnd(k,j))
            else
               big = upper(j) + dbnd(j,k)
               small = max(dbnd(k,j)-upper(j),lower(j)-dbnd(j,k))
            end if
            enter = .false.
            if (upper(k) .gt. big) then
               upper(k) = big
               if (.not. queued(k))  enter = .true.
            end if
            if (lower(k) .lt. small) then
               lower(k) = small
               if (.not. queued(k))  enter = .true.
            end if
c
c     enter a new candidate atom at the tail of the queue
c
            if (enter) then
               queued(k) = .true.
               if (head .eq. 0) then
                  head = k
                  tail = k
                  queue(k) = 0
               else
                  queue(tail) = k
                  queue(k) = 0
                  tail = k
               end if
            end if
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (iarc)
      deallocate (queue)
      deallocate (queued)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine trifix  --  update triangle inequality bounds  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "trifix" rebuilds both the upper and lower distance bound
c     matrices following tightening of one or both of the bounds
c     between a specified pair of atoms, "p" and "q", using a
c     modification of Murchland's shortest path update algorithm
c
c     literature references:
c
c     P. A. Steenbrink, "Optimization of Transport Networks", John
c     Wiley and Sons, Bristol, 1974; see section 7.7
c
c     R. Dionne, "Etude et Extension d'un Algorithme de Murchland",
c     Infor, 16, 132-146 (1978)
c
c
      subroutine trifix (p,q)
      use sizes
      use atoms
      use disgeo
      use inform
      use iounit
      implicit none
      integer i,k,p,q
      integer ip,iq,np,nq
      integer, allocatable :: pt(:)
      integer, allocatable :: qt(:)
      real*8 eps,ipmin,ipmax
      real*8 iqmin,iqmax
      real*8, allocatable :: pmin(:)
      real*8, allocatable :: pmax(:)
      real*8, allocatable :: qmin(:)
      real*8, allocatable :: qmax(:)
      logical, allocatable :: pun(:)
      logical, allocatable :: qun(:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (pt(n))
      allocate (qt(n))
      allocate (pmin(n))
      allocate (pmax(n))
      allocate (qmin(n))
      allocate (qmax(n))
      allocate (pun(n))
      allocate (qun(n))
c
c     initialize the set of nodes that may have changed bounds
c
      eps = 1.0d-10
      np = 0
      nq = 0
      do i = 1, n
         pun(i) = .true.
         qun(i) = .true.
      end do
c
c     store the upper and lower bounds to "p" and "q"
c
      do i = 1, p
         pmin(i) = dbnd(p,i)
         pmax(i) = dbnd(i,p)
      end do
      do i = p+1, n
         pmin(i) = dbnd(i,p)
         pmax(i) = dbnd(p,i)
      end do
      do i = 1, q
         qmin(i) = dbnd(q,i)
         qmax(i) = dbnd(i,q)
      end do
      do i = q+1, n
         qmin(i) = dbnd(i,q)
         qmax(i) = dbnd(q,i)
      end do
c
c     check for changes in the upper bounds to "p" and "q"
c
      do i = 1, n
         ipmax = qmax(p) + qmax(i)
         if (pmax(i) .gt. ipmax+eps) then
            np = np + 1
            pt(np) = i
            pmax(i) = ipmax
            pun(i) = .false.
         end if
         iqmax = pmax(q) + pmax(i)
         if (qmax(i) .gt. iqmax+eps) then
            nq = nq + 1
            qt(nq) = i
            qmax(i) = iqmax
            qun(i) = .false.
         end if
      end do
c
c     for node pairs whose bounds to "p" and "q" have changed,
c     make any needed changes to upper bound of the pair
c
      do ip = 1, np
         i = pt(ip)
         ipmax = pmax(i)
         do iq = 1, nq
            k = qt(iq)
            if (i .lt. k) then
               dbnd(i,k) = min(dbnd(i,k),ipmax+pmax(k))
            else
               dbnd(k,i) = min(dbnd(k,i),ipmax+pmax(k))
            end if
         end do
      end do
c
c     check for changes in the lower bounds to "p" and "q"
c
      do i = 1, n
         ipmin = max(qmin(p)-qmax(i),qmin(i)-qmax(p))
         if (pmin(i) .lt. ipmin-eps) then
            if (pun(i)) then
               np = np + 1
               pt(np) = i
            end if
            pmin(i) = ipmin
         end if
         iqmin = max(pmin(q)-pmax(i),pmin(i)-pmax(q))
         if (qmin(i) .lt. iqmin-eps) then
            if (qun(i)) then
               nq = nq + 1
               qt(nq) = i
            end if
            qmin(i) = iqmin
         end if
      end do
c
c     for node pairs whose bounds to "p" and "q" have changed,
c     make any needed changes to lower bound of the pair
c
      do ip = 1, np
         i = pt(ip)
         ipmin = pmin(i)
         ipmax = pmax(i)
         do iq = 1, nq
            k = qt(iq)
            if (i .lt. k) then
               dbnd(k,i) = max(dbnd(k,i),ipmin-pmax(k),pmin(k)-ipmax)
            else
               dbnd(i,k) = max(dbnd(i,k),ipmin-pmax(k),pmin(k)-ipmax)
            end if
         end do
      end do
c
c     update the upper and lower bounds to "p" and "q"
c
      do i = 1, p
         dbnd(p,i) = pmin(i)
         dbnd(i,p) = pmax(i)
      end do
      do i = p+1, n
         dbnd(i,p) = pmin(i)
         dbnd(p,i) = pmax(i)
      end do
      do i = 1, q
         dbnd(q,i) = qmin(i)
         dbnd(i,q) = qmax(i)
      end do
      do i = q+1, n
         dbnd(i,q) = qmin(i)
         dbnd(q,i) = qmax(i)
      end do
c
c     output the atoms updated and amount of work required
c
c     if (debug) then
c        write (iout,10)  p,q,np*nq
c  10    format (' TRIFIX  --  Bounds Update for Atoms',2i6,
c    &              ' with',i8,' Searches')
c     end if
c
c     perform deallocation of some local arrays
c
      deallocate (pt)
      deallocate (qt)
      deallocate (pmin)
      deallocate (pmax)
      deallocate (qmin)
      deallocate (qmax)
      deallocate (pun)
      deallocate (qun)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine grafic  --  schematic graphical matrix output  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "grafic" outputs the upper & lower triangles and diagonal
c     of a square matrix in a schematic form for visual inspection
c
c
      subroutine grafic (n,a,title)
      use iounit
      implicit none
      integer i,j,k,m,n
      integer maxj,nrow,ndash
      integer minrow,maxrow
      integer trimtext
      real*8 big,v
      real*8 amin,dmin,bmin
      real*8 amax,dmax,bmax
      real*8 rcl,scl,tcl
      real*8 ca,cb,cc,cd
      real*8 cw,cx,cy,cz
      real*8 a(n,*)
      character*1 dash
      character*1 ta,tb,tc,td,te
      character*1 digit(0:9)
      character*1 symbol(130)
      character*240 title
      data dash   / '-' /
      data ta,tb,tc,td,te  / ' ','.','+','X','#' /
      data digit  / '0','1','2','3','4','5','6','7','8','9' /
c
c
c     set bounds of length of print row and write the header
c
      minrow = 54
      maxrow = 130
      ndash = min(max(n,minrow),maxrow)
      write (iout,10)  (dash,i=1,ndash)
   10 format (/,1x,130a1)
      write (iout,20)  title(1:trimtext(title))
   20 format (/,1x,a)
c
c     find the maximum and minimum elements of the matrix
c
      big = 1.0d6
      dmax = -big
      dmin = big
      amax = -big
      amin = big
      bmax = -big
      bmin = big
      do i = 1, n
         if (a(i,i) .gt. dmax)  dmax = a(i,i)
         if (a(i,i) .lt. dmin)  dmin = a(i,i)
         do j = 1, i-1
            if (a(j,i) .gt. amax)  amax = a(j,i)
            if (a(j,i) .lt. amin)  amin = a(j,i)
            if (a(i,j) .gt. bmax)  bmax = a(i,j)
            if (a(i,j) .lt. bmin)  bmin = a(i,j)
         end do
      end do
      write (iout,30)  amin,amax,dmin,dmax,bmin,bmax
   30 format (/,' Range of Above Diag Elements : ',f13.4,' to',f13.4,
     &        /,' Range of Diagonal Elements :   ',f13.4,' to',f13.4,
     &        /,' Range of Below Diag Elements : ',f13.4,' to',f13.4)
c
c     now, print out the graphical representation
c
      write (iout,40)
   40 format (/,' Symbol Magnitude Ordering :',14x,
     &          '# > X > + > . > '' ''',/)
      rcl = (bmax-bmin) / 5.0d0
      scl = (amax-amin) / 5.0d0
      tcl = (dmax-dmin) / 9.0d0
      if (rcl .eq. 0.0d0)  rcl = 1.0d0
      if (scl .eq. 0.0d0)  scl = 1.0d0
      if (tcl .eq. 0.0d0)  tcl = 1.0d0
      ca = amin + scl
      cb = ca + scl
      cc = cb + scl
      cd = cc + scl
      cw = bmin + rcl
      cx = cw + rcl
      cy = cx + rcl
      cz = cy + rcl
      do j = 1, n, maxrow
         maxj = j + maxrow - 1
         if (maxj .gt. n)  maxj = n
         nrow = maxj - j + 1
         do i = 1, n
            do k = j, maxj
               m = k - j + 1
               if (k .lt. i) then
                  v = abs(a(i,k))
                  if (v .le. cw) then
                     symbol(m) = ta
                  else if (v .le. cx) then
                     symbol(m) = tb
                  else if (v .le. cy) then
                     symbol(m) = tc
                  else if (v .le. cz) then
                     symbol(m) = td
                  else
                     symbol(m) = te
                  end if
               else if (k .eq. i) then
                  symbol(m) = digit(nint((a(i,i)-dmin)/tcl))
               else if (k .gt. i) then
                  v = abs(a(i,k))
                  if (v .le. ca) then
                     symbol(m) = ta
                  else if (v .le. cb) then
                     symbol(m) = tb
                  else if (v .le. cc) then
                     symbol(m) = tc
                  else if (v .le. cd) then
                     symbol(m) = td
                  else
                     symbol(m) = te
                  end if
               end if
            end do
            write (iout,50)  (symbol(k),k=1,nrow)
   50       format (1x,130a1)
         end do
         write (iout,60)  (dash,i=1,ndash)
   60    format (/,1x,130a1)
         if (maxj .lt. n) then
            write (iout,70)
   70       format ()
         end if
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine dstmat  --  choose values for distance matrix  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "dstmat" selects a distance matrix containing values between
c     the previously smoothed upper and lower bounds; the distance
c     values are chosen from uniform distributions, in a triangle
c     correlated fashion, or using random partial metrization
c
c
      subroutine dstmat (dmx)
      use sizes
      use atoms
      use disgeo
      use inform
      use iounit
      use keys
      implicit none
      integer i,j,k,m
      integer index,next
      integer npart,npair
      integer nmetrize
      integer mik,mjk,nik,njk
      integer, allocatable :: list(:)
      real*8 random,fraction
      real*8 invbeta,alpha,beta
      real*8 corr,mean,stdev
      real*8 denom,swap,delta
      real*8 percent,eps,gap
      real*8 wall,cpu
      real*8, allocatable :: value(:)
      real*8 dmx(n,*)
      logical first,uniform
      logical update
      character*8 method
      character*20 keyword
      character*240 record
      character*240 string
      external random,invbeta
      save first,method,update
      save npart,percent
      save mean,stdev
      save alpha,beta
      data first  / .true. /
c
c
c     initialize the method for distance element selection
c
      if (first) then
         first = .false.
         method = 'PAIRWISE'
         uniform = .false.
         update = .true.
         npart = 0
         percent = 0.0d0
         mean = 0.0d0
         compact = 0.0d0
         beta = 4.0d0
c
c     search each line of the keyword file for options
c
         do i = 1, nkey
            next = 1
            record = keyline(i)
            call gettext (record,keyword,next)
            call upcase (keyword)
c
c     get a distance selection method and extent of metrization
c
            if (keyword(1:15) .eq. 'TRIAL-DISTANCE ') then
               call gettext (record,method,next)
               call upcase (method)
               if (method .eq. 'HAVEL') then
                  call getnumb (record,npart,next)
               else if (method .eq. 'PARTIAL') then
                  call getnumb (record,npart,next)
               else if (method .eq. 'PAIRWISE') then
                  string = record(next:240)
                  read (string,*,err=10,end=10)  percent
   10             continue
               end if
c
c     get a choice of initial mean for the trial distribution
c
            else if (keyword(1:19) .eq. 'TRIAL-DISTRIBUTION ') then
               string = record(next:240)
               read (string,*,err=20,end=20)  mean
   20          continue
               update = .false.
            end if
         end do
c
c     set extent of partial metrization during distance selection
c
         if (method .eq. 'HAVEL') then
            if (npart.le.0 .or. npart.ge.n-1)  npart = n
         else if (method .eq. 'PARTIAL') then
            if (npart.le.0 .or. npart.ge.n-1)  npart = 4
         else if (method .eq. 'PAIRWISE') then
            if (percent.le.0.0d0 .or. percent.ge.100.0d0)
     &         percent = min(100.0d0,2000.0d0/dble(n))
         end if
c
c     set the initial distribution for selection of trial distances
c
         if (method .eq. 'CLASSIC')  uniform = .true.
         if (method .eq. 'TRICOR')  uniform = .true.
         if (method .eq. 'HAVEL')  uniform = .true.
         if (uniform)  update = .false.
         if (update) then
c           mean = 2.35d0 / log(pathmax)
            mean = 1.65d0 / (pathmax)**0.25d0
c           mean = 1.30d0 / (pathmax)**0.20d0
         end if
         alpha = beta*mean / (1.0d0-mean)
         stdev = sqrt(alpha*beta/(alpha+beta+1.0d0)) / (alpha+beta)
      end if
c
c     write out the final choice for distance matrix generation
c
      if (verbose) then
         call settime
         if (method .eq. 'CLASSIC') then
            write (iout,30)
   30       format (/,' Distance Matrix via Uniform Random',
     &                 ' Fractions without Metrization :')
         else if (method .eq. 'RANDOM') then
            write (iout,40)
   40       format (/,' Distance Matrix Generated via Normal',
     &                 ' Fractions without Metrization :')
         else if (method .eq. 'TRICOR') then
            write (iout,50)
   50       format (/,' Distance Matrix Generated via Triangle',
     &                 ' Correlated Fractions :')
         else if (method.eq.'HAVEL' .and. npart.lt.n) then
            write (iout,60)  npart
   60       format (/,' Distance Matrix Generated via',i4,'-Atom',
     &                 ' Partial Metrization :')
         else if (method .eq. 'HAVEL') then
            write (iout,70)
   70       format (/,' Distance Matrix Generated via Randomized',
     &                 ' Atom-Based Metrization :')
         else if (method .eq. 'PARTIAL') then
            write (iout,80)  npart
   80       format (/,' Distance Matrix Generated via',i4,'-Atom',
     &                 ' Partial Metrization :')
         else if (method.eq.'PAIRWISE' .and. percent.lt.100.0d0) then
            write (iout,90)  percent
   90       format (/,' Distance Matrix Generated via',f6.2,'%',
     &                 ' Random Pairwise Metrization :')
         else
            write (iout,100)
  100       format (/,' Distance Matrix Generated via Randomized',
     &                 ' Pairwise Metrization :')
         end if
      end if
c
c     adjust the distribution for selection of trial distances
c
      if (uniform) then
         write (iout,110)
  110    format (/,' Trial Distances Selected at Random from',
     &              ' Uniform Distribution')
      else
         if (update) then
            alpha = alpha - 0.2d0*sign(sqrt(abs(compact)),compact)
            mean = alpha / (alpha+beta)
            stdev = sqrt(alpha*beta/(alpha+beta+1.0d0)) / (alpha+beta)
         end if
         write (iout,120)  mean,stdev,alpha,beta
  120    format (/,' Trial Distance Beta Distribution :',
     &              4x,f5.2,' +/-',f5.2,3x,'Alpha-Beta',2f6.2)
      end if
c
c     perform dynamic allocation of some local arrays
c
      npair = n*(n-1) / 2
      allocate (list(npair))
      allocate (value(npair))
c
c     uniform or Gaussian distributed distances without metrization
c
      if (method.eq.'CLASSIC' .or. method.eq.'RANDOM') then
         do i = 1, n
            dmx(i,i) = 0.0d0
         end do
         do i = 1, n-1
            do j = i+1, n
               fraction = random ()
               if (method .eq. 'RANDOM') then
                  fraction = invbeta (alpha,beta,fraction)
               end if
               delta = dbnd(i,j) - dbnd(j,i)
               dmx(j,i) = dbnd(j,i) + delta*fraction
               dmx(i,j) = dmx(j,i)
            end do
         end do
c
c     Crippen's triangle correlated distance selection
c
      else if (method .eq. 'TRICOR') then
         do i = 1, n
            dmx(i,i) = 0.0d0
         end do
         do i = 1, n-1
            do j = i+1, n
               dmx(j,i) = random ()
               dmx(i,j) = dmx(j,i)
            end do
         end do
         do i = 1, n-1
            do j = i+1, n
               denom = 0.0d0
               dmx(i,j) = 0.0d0
               do k = 1, n
                  if (k .ne. i) then
                     mik = max(i,k)
                     mjk = max(j,k)
                     nik = min(i,k)
                     njk = min(j,k)
                     if (k .eq. j) then
                        dmx(i,j) = dmx(i,j) + dmx(j,i)
                        denom = denom + 1.0d0
                     else if (dbnd(njk,mjk) .le.
     &                        0.2d0*dbnd(nik,mik)) then
                        if (i .gt. k)  corr = 0.9d0 * dmx(i,k)
                        if (k .gt. i)  corr = 0.9d0 * dmx(k,i)
                        dmx(i,j) = dmx(i,j) + corr
                        denom = denom + 0.9d0
                     else if (dbnd(nik,mik) .le.
     &                        0.2d0*dbnd(njk,mjk)) then
                        if (j .gt. k)  corr = 0.9d0 * dmx(j,k)
                        if (k .gt. j)  corr = 0.9d0 * dmx(k,j)
                        dmx(i,j) = dmx(i,j) + corr
                        denom = denom + 0.9d0
                     else if (dbnd(mik,nik) .ge.
     &                        0.9d0*dbnd(njk,mjk)) then
                        if (j .gt. k)  corr = 0.5d0 * (1.0d0-dmx(j,k))
                        if (k .gt. j)  corr = 0.5d0 * (1.0d0-dmx(k,j))
                        dmx(i,j) = dmx(i,j) + corr
                        denom = denom + 0.5d0
                     else if (dbnd(mjk,njk) .ge.
     &                        0.9d0*dbnd(nik,mik)) then
                        if (i .gt. k)  corr = 0.5d0 * (1.0d0-dmx(i,k))
                        if (k .gt. i)  corr = 0.5d0 * (1.0d0-dmx(k,i))
                        dmx(i,j) = dmx(i,j) + corr
                        denom = denom + 0.5d0
                     end if
                  end if
               end do
               dmx(i,j) = dmx(i,j) / denom
            end do
         end do
         do i = 1, n-1
            do j = i+1, n
               delta = dbnd(i,j) - dbnd(j,i)
               dmx(i,j) = dbnd(j,i) + delta*dmx(i,j)
               dmx(j,i) = dmx(i,j)
            end do
         end do
c
c     Havel/XPLOR atom-based metrization over various distributions
c
      else if (method.eq.'HAVEL' .or. method.eq.'PARTIAL') then
         do i = 1, n
            do j = 1, n
               dmx(j,i) = dbnd(j,i)
            end do
         end do
         do i = 1, n
            value(i) = random ()
         end do
         call sort2 (n,value,list)
         gap = 0.0d0
         do i = 1, n-1
            k = list(i)
            do j = i+1, n
               m = list(j)
               fraction = random ()
               if (method .eq. 'PARTIAL') then
                  fraction = invbeta (alpha,beta,fraction)
               end if
               delta = abs(dbnd(k,m) - dbnd(m,k))
               if (k .lt. m) then
                  dbnd(k,m) = dbnd(m,k) + delta*fraction
                  dbnd(m,k) = dbnd(k,m)
               else
                  dbnd(k,m) = dbnd(k,m) + delta*fraction
                  dbnd(m,k) = dbnd(k,m)
               end if
               if (i .le. npart)  call trifix (k,m)
               if (i .gt. npart)  gap = gap + delta
            end do
         end do
         do i = 1, n
            do j = 1, n
               swap = dmx(j,i)
               dmx(j,i) = dbnd(j,i)
               dbnd(j,i) = swap
            end do
         end do
         if (verbose .and. npart.lt.n-1) then
            write (iout,130)  gap/dble((n-npart)*(n-npart-1)/2)
  130       format (/,' Average Bound Gap after Partial Metrization :',
     &                 3x,f12.4)
         end if
c
c     use partial randomized pairwise distance-based metrization
c
      else if (method.eq.'PAIRWISE' .and. percent.le.10.0d0) then
         npair = n*(n-1) / 2
         nmetrize = nint(0.01d0*percent*dble(npair))
         do i = 1, n
            do j = 1, n
               dmx(j,i) = dbnd(j,i)
            end do
         end do
         do i = 1, nmetrize
  140       continue
            k = int(dble(n)*random()) + 1
            m = int(dble(n)*random()) + 1
            if (dbnd(k,m) .eq. dbnd(m,k))  goto 140
            if (k .gt. m) then
               j = k
               k = m
               m = j
            end if
            fraction = random ()
            fraction = invbeta (alpha,beta,fraction)
            delta = dbnd(k,m) - dbnd(m,k)
            dbnd(k,m) = dbnd(m,k) + delta*fraction
            dbnd(m,k) = dbnd(k,m)
            call trifix (k,m)
         end do
         gap = 0.0d0
         do i = 1, n-1
            do j = i, n
               delta = dbnd(i,j) - dbnd(j,i)
               if (delta .ne. 0.0d0) then
                  gap = gap + delta
                  fraction = random ()
                  fraction = invbeta (alpha,beta,fraction)
                  dbnd(i,j) = dbnd(j,i) + delta*fraction
                  dbnd(j,i) = dbnd(i,j)
               end if
            end do
         end do
         do i = 1, n
            do j = 1, n
               swap = dmx(j,i)
               dmx(j,i) = dbnd(j,i)
               dbnd(j,i) = swap
            end do
         end do
         if (verbose .and. nmetrize.lt.npair) then
            write (iout,150)  gap/dble(npair-nmetrize)
  150       format (/,' Average Bound Gap after Partial Metrization :',
     &                 3x,f12.4)
         end if
c
c     use randomized pairwise distance-based metrization
c
      else if (method .eq. 'PAIRWISE') then
         npair = n*(n-1) / 2
         nmetrize = nint(0.01d0*percent*dble(npair))
         do i = 1, n
            do j = 1, n
               dmx(j,i) = dbnd(j,i)
            end do
         end do
         do i = 1, npair
            value(i) = random ()
         end do
         call sort2 (npair,value,list)
         eps = 1.0d-10
         gap = 0.0d0
         do i = 1, npair
            index = list(i)
            k = int(0.5d0 * (dble(2*n+1)
     &                 - sqrt(dble(4*n*(n-1)-8*index+9))) + eps)
            m = n*(1-k) + k*(k+1)/2 + index
            fraction = random ()
            fraction = invbeta (alpha,beta,fraction)
            delta = dbnd(k,m) - dbnd(m,k)
            dbnd(k,m) = dbnd(m,k) + delta*fraction
            dbnd(m,k) = dbnd(k,m)
            if (i .le. nmetrize)  call trifix (k,m)
            if (i .gt. nmetrize)  gap = gap + delta
         end do
         do i = 1, n
            do j = 1, n
               swap = dmx(j,i)
               dmx(j,i) = dbnd(j,i)
               dbnd(j,i) = swap
            end do
         end do
         if (verbose .and. nmetrize.lt.npair) then
            write (iout,160)  gap/dble(npair-nmetrize)
  160       format (/,' Average Bound Gap after Partial Metrization :',
     &                 3x,f12.4)
         end if
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (list)
      deallocate (value)
c
c     get the time required for distance matrix generation
c
      if (verbose) then
         call gettime (wall,cpu)
         write (iout,170)  wall
  170    format (/,' Time Required for Distance Matrix :',5x,f12.2,
     &              ' seconds')
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine metric  --  computation of the metric matrix  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "metric" takes as input the trial distance matrix and computes
c     the metric matrix of all possible dot products between the atomic
c     vectors and the center of mass using the law of cosines and the
c     following formula for the distances to the center of mass:
c
c        dcm(i)**2 = (1/n) * sum(j=1,n)(dist(i,j)**2)
c                          - (1/n**2) * sum(j<k)(dist(j,k)**2)
c
c     upon output, the metric matrix is stored in the lower triangle
c     plus diagonal of the input trial distance matrix, the upper
c     triangle of the input matrix is unchanged
c
c     literature reference:
c
c     G. M. Crippen and T. F. Havel, "Stable Calculation of Coordinates
c     from Distance Information", Acta Cryst., A34, 282-284 (1978)
c
c
      subroutine metric (gmx,nneg)
      use sizes
      use atoms
      use disgeo
      use inform
      use iounit
      implicit none
      integer i,j,nneg
      real*8 total,sum,rg
      real*8 gmx(n,*)
      real*8, allocatable :: dsq(:)
      real*8, allocatable :: dcm(:)
c
c
c     square and sum trial distances to get radius of gyration
c
      total = 0.0d0
      do i = 1, n
         do j = i, n
            gmx(j,i) = gmx(j,i)**2
            total = total + gmx(j,i)
         end do
      end do
      total = total / dble(n**2)
      rg = sqrt(total)
      if (verbose) then
         write (iout,10)  rg
   10    format (/,' Radius of Gyration before Embedding :',7x,f16.4)
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (dsq(n))
      allocate (dcm(n))
c
c     sum squared distances from each atom; the center
c     of mass is derived using the formula shown above
c
      nneg = 0
      do i = 1, n
         sum = 0.0d0
         do j = 1, i-1
            sum = sum + gmx(i,j)
         end do
         do j = i, n
            sum = sum + gmx(j,i)
         end do
         dsq(i) = sum/dble(n) - total
         dcm(i) = sqrt(abs(dsq(i)))
         if (dsq(i) .lt. 0.0d0) then
            nneg = nneg + 1
            dcm(i) = -dcm(i)
         end if
      end do
      if (verbose .and. n.le.130) then
         write (iout,20)
   20    format (/,' Atomic Distances to the Center of Mass :',/)
         write (iout,30)  (dcm(i),i=1,n)
   30    format (6f13.4)
      end if
c
c     calculate the metric matrix using the law of cosines, and
c     place into the lower triangle of the input distance matrix
c
      do i = 1, n
         do j = i, n
            gmx(j,i) = 0.5d0 * (dsq(i)+dsq(j)-gmx(j,i))
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (dsq)
      deallocate (dcm)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine eigen  --  largest eigenvalues of metric metrix  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "eigen" uses the power method to compute the largest eigenvalues
c     and eigenvectors of the metric matrix, "valid" is set true if the
c     first three eigenvalues are positive
c
c
      subroutine eigen (evl,evc,gmx,valid)
      use sizes
      use atoms
      use inform
      use iounit
      implicit none
      integer i,j,neigen
      real*8 wall,cpu
      real*8 evl(*)
      real*8 evc(n,*)
      real*8 gmx(n,*)
      logical valid
c
c
c     initialize number of eigenvalues and convergence criteria
c
      if (verbose)  call settime
      neigen = 3
c
c     compute largest eigenvalues via power method with deflation
c
      call deflate (n,neigen,gmx,evl,evc)
c
c     check to see if the first three eigenvalues are positive
c
      valid = .true.
      do i = 1, 3
         if (evl(i) .lt. 0.0d0)  valid = .false.
      end do
c
c     print out the eigenvalues and their eigenvectors
c
      if (verbose) then
         write (iout,10)
   10    format (/,' Eigenvalues from Metric Matrix :',/)
         write (iout,20)  (evl(i),i=1,neigen)
   20    format (5f15.4)
      end if
      if (debug) then
         write (iout,30)
   30    format (/,' Eigenvectors from Metric Matrix :',/)
         do i = 1, n
            write (iout,40)  (evc(i,j),j=1,neigen)
   40       format (5f15.4)
         end do
      end if
c
c     get the time required for partial matrix diagonalization
c
      if (verbose) then
         call gettime (wall,cpu)
         write (iout,50)  wall
   50    format (/,' Time Required for Eigenvalues :',9x,f12.2,
     &              ' seconds')
      end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine coords  --  converts eigenvalues to coordinates  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "coords" converts the three principal eigenvalues/vectors from
c     the metric matrix into atomic coordinates, and calls a routine
c     to compute the rms deviation from the bounds
c
c
      subroutine coords (evl,evc)
      use sizes
      use atoms
      use disgeo
      use inform
      use iounit
      implicit none
      integer i,j,neigen
      real*8 rg
      real*8 evl(*)
      real*8 evc(n,*)
      character*240 title
c
c
c     compute coordinates from the largest eigenvalues and vectors
c
      neigen = 3
      do j = 1, neigen
         evl(j) = sqrt(abs(evl(j)))
      end do
      do j = 1, neigen
         do i = 1, n
            evc(i,j) = evl(j) * evc(i,j)
         end do
      end do
c
c     transfer the final coordinates back to atomic vectors
c
      do i = 1, n
         x(i) = evc(i,1)
         y(i) = evc(i,2)
         z(i) = evc(i,3)
      end do
c
c     find the rms bounds deviations and radius of gyration
c
      if (verbose) then
         title = 'after Embedding :'
         call rmserror (title)
         call gyrate (rg)
         write (iout,10)  rg
   10    format (/,' Radius of Gyration after Embedding :',8x,f16.4)
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine chksize  --  estimate compaction or expansion  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "chksize" computes a measure of overall global structural
c     expansion or compaction from the number of excess upper
c     or lower bounds matrix violations
c
c
      subroutine chksize
      use sizes
      use atoms
      use couple
      use disgeo
      use inform
      use iounit
      implicit none
      integer i,k,npair,nskip
      integer nlarge,nsmall
      integer, allocatable :: skip(:)
      real*8 xi,yi,zi
      real*8 dstsq,bupsq,blosq
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (skip(n))
c
c     zero out the list of atoms locally connected to each atom
c
      nskip = 0
      do i = 1, n
         skip(i) = 0
      end do
c
c     initialize counters, total pair number, and cutoff distance
c
      nlarge = 0
      nsmall = 0
      npair = n*(n-1) / 2
c
c     count the number of excess upper or lower bound violations
c
      do i = 1, n-1
         xi = x(i)
         yi = y(i)
         zi = z(i)
c        do k = 1, n12(i)
c           skip(i12(k,i)) = i
c        end do
c        do k = 1, n13(i)
c           skip(i13(k,i)) = i
c        end do
c        do k = 1, n14(i)
c           skip(i14(k,i)) = i
c        end do
         do k = i+1, n
            if (skip(k) .eq. i) then
               nskip = nskip + 1
            else
               dstsq = (x(k)-xi)**2 + (y(k)-yi)**2 + (z(k)-zi)**2
               bupsq = dbnd(i,k)**2
               blosq = dbnd(k,i)**2
               if (dstsq .gt. bupsq) then
                  nlarge = nlarge + 1
               else if (blosq .gt. dstsq) then
                  nsmall = nsmall + 1
               end if
            end if
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (skip)
c
c     set the value for the overall index of compaction
c
      compact = 100.0d0 * dble(nlarge-nsmall)/dble(npair-nskip)
      if (verbose) then
         write (iout,10)  compact
   10    format (/,' Index of Structure Expansion/Compaction :',
     &              7x,f12.4)
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine majorize  --  Guttman transform majorization  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "majorize" refines the projected coordinates by attempting to
c     minimize the least square residual between the trial distance
c     matrix and the distances computed from the coordinates
c
c     literature reference:
c
c     T. F. Havel, "An Evaluation of Computational Strategies for
c     Use in the Determination of Protein Structure from Distance
c     Constraints obtained by Nuclear Magnetic Resonance", Progress
c     in Biophysics and Molecular Biology, 56, 43-78 (1991)
c
c
      subroutine majorize (dmx)
      use sizes
      use atoms
      use inform
      use iounit
      implicit none
      integer i,k,iter
      integer niter,period
      real*8 pairs,rg
      real*8 dn1,dn2
      real*8 wall,cpu
      real*8 target,dist,error
      real*8 rmserr,average
      real*8 xi,yi,zi
      real*8, allocatable :: b(:)
      real*8, allocatable :: xx(:)
      real*8, allocatable :: yy(:)
      real*8, allocatable :: zz(:)
      real*8 dmx(n,*)
      character*240 title
c
c
c     set number of iterations and some other needed values
c
      if (verbose)  call settime
      niter = 20
      period = 5
      pairs = dble(n*(n-1)/2)
      dn1 = dble(n-1)
      dn2 = dble(n*n)
c
c     find the average and rms error from trial distances
c
      iter = 0
      rmserr = 0.0d0
      average = 0.0d0
      do i = 1, n-1
         xi = x(i)
         yi = y(i)
         zi = z(i)
         do k = i+1, n
            target = dmx(k,i)
            dist = sqrt((x(k)-xi)**2+(y(k)-yi)**2+(z(k)-zi)**2)
            error = dist - target
            rmserr = rmserr + error**2
            average = average + error/target
         end do
      end do
      rmserr = sqrt(rmserr/pairs)
      average = 100.0d0 * average / pairs
c
c     write a header with the initial error values
c
      if (verbose) then
         write (iout,10)
   10    format (/,' Majorization to Trial Distances using',
     &              ' Constant Weights :',
     &           //,4x,'Iteration',6x,'RMS Error',5x,'Ave % Error',/)
         write (iout,20)  iter,rmserr,average
   20    format (5x,i5,2f16.4)
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (b(n))
      allocate (xx(n))
      allocate (yy(n))
      allocate (zz(n))
c
c     initialize the transformed coordinates for each atom
c
      do iter = 1, niter
         do i = 1, n
            xi = x(i)
            yi = y(i)
            zi = z(i)
            b(i) = 0.0d0
            xx(i) = 0.0d0
            yy(i) = 0.0d0
            zz(i) = 0.0d0
c
c     form a single row of the B matrix assuming unity weights
c
            do k = 1, i-1
               dist = sqrt((x(k)-xi)**2+(y(k)-yi)**2+(z(k)-zi)**2)
               b(k) = -dmx(k,i) / dist
               b(i) = b(i) - b(k)
            end do
            do k = i+1, n
               dist = sqrt((x(k)-xi)**2+(y(k)-yi)**2+(z(k)-zi)**2)
               b(k) = -dmx(k,i) / dist
               b(i) = b(i) - b(k)
            end do
c
c     multiply the row of the B matrix by the atomic coordinates
c
            do k = 1, n
               xx(i) = xx(i) + b(k)*x(k)
               yy(i) = yy(i) + b(k)*y(k)
               zz(i) = zz(i) + b(k)*z(k)
            end do
         end do
c
c     move the intermediate values into the coordinate arrays
c
         do i = 1, n
            x(i) = xx(i)
            y(i) = yy(i)
            z(i) = zz(i)
         end do
c
c     multiply the inverse weight matrix S+ by the coordinates
c
         do i = 1, n
            xx(i) = (dn1/dn2) * x(i)
            yy(i) = (dn1/dn2) * y(i)
            zz(i) = (dn1/dn2) * z(i)
            do k = 1, i-1
               xx(i) = xx(i) - x(k)/dn2
               yy(i) = yy(i) - y(k)/dn2
               zz(i) = zz(i) - z(k)/dn2
            end do
            do k = i+1, n
               xx(i) = xx(i) - x(k)/dn2
               yy(i) = yy(i) - y(k)/dn2
               zz(i) = zz(i) - z(k)/dn2
            end do
         end do
c
c     copy the new coordinates into their permanent arrays
c
         do i = 1, n
            x(i) = xx(i)
            y(i) = yy(i)
            z(i) = zz(i)
         end do
c
c     find the average and rms error from trial distances
c
         rmserr = 0.0d0
         average = 0.0d0
         do i = 1, n-1
            xi = x(i)
            yi = y(i)
            zi = z(i)
            do k = i+1, n
               target = dmx(k,i)
               dist = sqrt((x(k)-xi)**2+(y(k)-yi)**2+(z(k)-zi)**2)
               error = dist - target
               rmserr = rmserr + error**2
               average = average + error/target
            end do
         end do
         rmserr = sqrt(rmserr/pairs)
         average = 100.0d0 * average / pairs
         if (verbose .and. mod(iter,period).eq.0) then
            write (iout,30)  iter,rmserr,average
   30       format (5x,i5,2f16.4)
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (b)
      deallocate (xx)
      deallocate (yy)
      deallocate (zz)
c
c     find the rms bounds deviations and radius of gyration
c
      if (verbose) then
         title = 'after Majorization :'
         call rmserror (title)
         call gyrate (rg)
         write (iout,40)  rg
   40    format (/,' Radius of Gyration after Majorization :',5x,f16.4)
c
c     get the time required for the majorization procedure
c
         call gettime (wall,cpu)
         write (iout,50)  wall
   50    format (/,' Time Required for Majorization :',8x,f12.2,
     &              ' seconds')
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine refine  --  minimization of initial embedding  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "refine" performs minimization of the atomic coordinates
c     of an initial crude embedded distance geometry structure versus
c     the bound, chirality, planarity and torsional error functions
c
c
      subroutine refine (mode,fctval,grdmin)
      use sizes
      use atoms
      use disgeo
      use inform
      use iounit
      use minima
      use output
      implicit none
      integer i,nvar
      real*8 initerr,miderr,toterr
      real*8 fctval,grdmin
      real*8, allocatable :: xx(:)
      character*7 mode
      external initerr,miderr
      external toterr,optsave
c
c
c     perform dynamic allocation of some local arrays
c
      nvar = 3 * n
      allocate (xx(nvar))
c
c     translate the atomic coordinates to optimization variables
c
      nvar = 0
      do i = 1, n
         nvar = nvar + 1
         xx(nvar) = x(i)
         nvar = nvar + 1
         xx(nvar) = y(i)
         nvar = nvar + 1
         xx(nvar) = z(i)
      end do
c
c     set values of parameters needed for optimization
c
      coordtype = 'NONE'
      cyclesave = .true.
c     grdmin = 0.01d0
      maxiter = 2 * nvar
      iwrite = 0
c     iprint = 0
c     if (verbose)  iprint = 10
c
c     minimize initially only on the local geometry and torsions,
c     then on local geometry and chirality, torsions, and finally
c     minimize on all distance bounds, chirality and torsions
c
      if (mode .eq. 'INITIAL') then
         call lbfgs (nvar,xx,fctval,grdmin,initerr,optsave)
      else if (mode .eq. 'MIDDLE') then
         call lbfgs (nvar,xx,fctval,grdmin,miderr,optsave)
      else if (mode .eq. 'FINAL') then
         call lbfgs (nvar,xx,fctval,grdmin,toterr,optsave)
      end if
c
c     translate optimization variables back to coordinates
c
      nvar = 0
      do i = 1, n
         nvar = nvar + 1
         x(i) = xx(nvar)
         nvar = nvar + 1
         y(i) = xx(nvar)
         nvar = nvar + 1
         z(i) = xx(nvar)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xx)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine explore  --  simulated annealing refinement  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "explore" uses simulated annealing on an initial crude
c     embedded distance geoemtry structure to refine versus the
c     bound, chirality, planarity and torsional error functions
c
c
      subroutine explore (mode,nstep,dt,mass,temp_start,temp_stop,v,a)
      use sizes
      use atoms
      use inform
      use iounit
      use math
      use units
      implicit none
      integer i,istep,nstep
      integer nvar,period
      real*8 error,total
      real*8 prior,change
      real*8 dt,dt2,dt_2,dt2_2
      real*8 xbig,xrms,mass,kinetic
      real*8 temp_start,temp_stop
      real*8 tau_start,tau_stop
      real*8 ratio,sigmoid,scale
      real*8 target,temp,tautemp
      real*8 initerr,miderr,toterr
      real*8 v(*)
      real*8 a(*)
      real*8, allocatable :: xx(:)
      real*8, allocatable :: xmove(:)
      real*8, allocatable :: g(:)
      character*7 mode
c
c
c     set values of the basic simulated annealing parameters
c
c     nstep = 5000
c     dt = 0.1d0
c     temp_start = 200.0d0
c     temp_stop = 0.0d0
c     mass = 1000.0d0
c
c     perform dynamic allocation of some local arrays
c
      nvar = 3 * n
      allocate (xx(nvar))
      allocate (xmove(nvar))
      allocate (g(nvar))
c
c     translate the atomic coordinates to annealing variables
c
      nvar = 0
      do i = 1, n
         nvar = nvar + 1
         xx(nvar) = x(i)
         nvar = nvar + 1
         xx(nvar) = y(i)
         nvar = nvar + 1
         xx(nvar) = z(i)
      end do
c
c     initialize the velocities, accelerations and other parameters
c
      dt2 = dt * dt
      dt_2 = dt / 2.0d0
      dt2_2 = dt2 / 2.0d0
      period = 100
      tau_start = 100.0d0 * dt
      tau_stop = 10.0d0 * dt
      tautemp = tau_start
c
c     print a header for the simulated annealing protocol
c
      write (iout,10)
   10 format (/,' Molecular Dynamics Simulated Annealing Refinement :')
      write (iout,20)  nstep,dt,log(mass)/logten,temp_start,temp_stop
   20 format (/,' Steps:',i6,3x,'Time/Step:',f6.3,' ps',3x,
     &           'LogMass:',f5.2,3x,'Temp:',f6.1,' to',f6.1)
c
c     get the total error and temperature at start of dynamics
c
      if (mode .eq. 'INITIAL') then
         error = initerr (xx,g)
      else if (mode .eq. 'MIDDLE') then
         error = miderr (xx,g)
      else if (mode .eq. 'FINAL') then
         error = toterr (xx,g)
      end if
      kinetic = 0.0d0
      do i = 1, nvar
         kinetic = kinetic + mass*v(i)**2
      end do
      kinetic = 0.5d0 * kinetic / convert
      temp = 2.0d0 * kinetic / (dble(nvar) * gasconst)
      total = error + kinetic
      prior = total
      if (verbose) then
         write (iout,30)
   30    format (/,' MD Step    E Total   E Potential   E Kinetic',
     &              '     Temp    MaxMove   RMS Move',/)
         write (iout,40)  0,total,error,kinetic,temp
   40    format (i6,2f13.4,f12.4,f11.2)
      end if
c
c     find new positions and half-step velocities via Verlet
c
      do istep = 1, nstep
         xbig = 0.0d0
         xrms = 0.0d0
         do i = 1, nvar
            xmove(i) = v(i)*dt + a(i)*dt2_2
            xx(i) = xx(i) + xmove(i)
            v(i) = v(i) + a(i)*dt_2
            if (abs(xmove(i)) .gt. xbig)  xbig = abs(xmove(i))
            xrms = xrms + xmove(i)**2
         end do
         xrms = sqrt(xrms/dble(nvar))
c
c     get the error function value and gradient
c
         if (mode .eq. 'INITIAL') then
            error = initerr (xx,g)
         else if (mode .eq. 'MIDDLE') then
            error = miderr (xx,g)
         else if (mode .eq. 'FINAL') then
            error = toterr (xx,g)
         end if
c
c     use Newton's second law to get the next accelerations;
c     find the full-step velocities using the Verlet recursion
c
         do i = 1, nvar
            a(i) = -convert * g(i) / mass
            v(i) = v(i) + a(i)*dt_2
         end do
c
c     find the total kinetic energy and system temperature
c
         kinetic = 0.0d0
         do i = 1, nvar
            kinetic = kinetic + mass*v(i)**2
         end do
         kinetic = 0.5d0 * kinetic / convert
         temp = 2.0d0 * kinetic / (dble(nvar) * gasconst)
         if (temp .eq. 0.0d0)  temp = 0.1d0
c
c     set target temperature and coupling via a sigmoidal cooling
c
         ratio = dble(istep) / dble(nstep)
         ratio = sigmoid (3.5d0,ratio)
         target = temp_start*(1.0d0-ratio) + temp_stop*ratio
         tautemp = tau_start*(1.0d0-ratio) + tau_stop*ratio
c
c     couple to external temperature bath via velocity scaling
c
         scale = sqrt(1.0d0 + (dt/tautemp)*(target/temp-1.0d0))
         do i = 1, nvar
            v(i) = scale * v(i)
         end do
c
c     write results for the current annealing step
c
         total = error + kinetic
         if (verbose .and. mod(istep,period).eq.0) then
            write (iout,50)  istep,total,error,kinetic,temp,xbig,xrms
   50       format (i6,2f13.4,f12.4,f11.2,2f10.4)
         end if
c
c     check the energy change for instability in the dynamics
c
         change = total - prior
         if (change .gt. dble(n)) then
            do i = 1, nvar
               xx(i) = xx(i) - xmove(i)
            end do
            if (verbose .and. mod(istep,period).ne.0) then
               write (iout,60)  istep,total,error,kinetic,temp,xbig,xrms
   60          format (i6,2f13.4,f12.4,f11.2,2f10.4)
            end if
            write (iout,70)
   70       format (/,' EXPLORE  --  Simulated Annealing Unstable;',
     &                 ' Switching to Minimization')
            goto 80
         end if
      end do
c
c     translate annealing variables back to coordinates
c
   80 continue
      nvar = 0
      do i = 1, n
         nvar = nvar + 1
         x(i) = xx(nvar)
         nvar = nvar + 1
         y(i) = xx(nvar)
         nvar = nvar + 1
         z(i) = xx(nvar)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xx)
      deallocate (xmove)
      deallocate (g)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine fracdist  --  fractional distance distribution  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "fracdist" computes a normalized distribution of the pairwise
c     fractional distances between the smoothed upper and lower bounds
c
c     literature reference:
c
c     C. M. Oshiro, J. Thomason and I. D. Kuntz, "Effects of Limited
c     Input Distance Constraints Upon the Distance Geometry Algorithm",
c     Biopolymers, 31, 1049-1064 (1991)
c
c
      subroutine fracdist (title)
      use sizes
      use atoms
      use disgeo
      use iounit
      implicit none
      integer start,stop
      parameter (start=-20)
      parameter (stop=120)
      integer i,j,k,sum
      integer leng,trimtext
      integer bin(start:stop)
      integer bin2(start:stop)
      real*8 xi,yi,zi,size
      real*8 dist,range,fraction
      real*8 fdist(start:stop)
      real*8 fdist2(start:stop)
      character*240 title
c
c
c     set the bin size and zero out the individual bins
c
      size = 0.01d0
      do i = start, stop
         bin(i) = 0
         bin2(i) = 0
      end do
c
c     get distribution of fractional distances between bounds
c
      do i = 1, n-1
         xi = x(i)
         yi = y(i)
         zi = z(i)
         do j = i+1, n
            dist = sqrt((x(j)-xi)**2 + (y(j)-yi)**2 + (z(j)-zi)**2)
            range = dbnd(i,j) - dbnd(j,i)
            if (range .ge. 1.0d0) then
               fraction = (dist-dbnd(j,i)) / range
               k = nint(fraction / size)
               k = max(start,min(stop,k))
               bin(k) = bin(k) + 1
               if (range .ge. 0.8d0*pathmax)  bin2(k) = bin2(k) + 1
            end if
         end do
      end do
c
c     normalize the fractional distance frequency distribution
c
      sum = 0
      do i = start, stop
         sum = sum + bin(i)
      end do
      do i = start, stop
         fdist(i) = dble(bin(i)) / (size*dble(sum))
      end do
      sum = 0
      do i = start, stop
         sum = sum + bin2(i)
      end do
      do i = start, stop
         fdist2(i) = dble(bin2(i)) / (size*dble(sum))
      end do
c
c     print the normalized fractional distance probability
c
      leng = trimtext(title)
      write (iout,10)  title(1:leng)
   10 format (/,' Fractional Distance Distribution ',a,/)
      do i = start, stop
         write (iout,20)  size*dble(i),fdist(i),fdist2(i)
   20    format (8x,f8.4,8x,f8.4,8x,f8.4)
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine rmserror  --  rms bound and restraint error  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "rmserror" computes the maximum absolute deviation and the
c     rms deviation from the distance bounds, and the number and
c     rms value of the distance restraint violations
c
c
      subroutine rmserror (title)
      use sizes
      use atoms
      use disgeo
      use iounit
      use restrn
      implicit none
      integer i,j,k,npair
      integer nhierr,nloerr
      integer ihi,jhi,ilo,jlo
      integer leng,trimtext
      real*8 rms,himax,lomax
      real*8 dist,hierr,loerr
      character*240 title
c
c
c     search all atom pairs for maximal bounds deviations
c
      npair = n*(n-1) / 2
      nloerr = 0
      nhierr = 0
      ilo = 0
      jlo = 0
      ihi = 0
      jhi = 0
      rms = 0.0d0
      lomax = 0.0d0
      himax = 0.0d0
      do i = 1, n-1
         do j = i+1, n
            dist = (x(i)-x(j))**2 + (y(i)-y(j))**2 + (z(i)-z(j))**2
            dist = sqrt(dist)
            hierr = dist - dbnd(i,j)
            if (hierr .gt. 0.0d0) then
               nhierr = nhierr + 1
               rms = rms + hierr**2
               if (hierr .gt. himax) then
                  himax = hierr
                  ihi = i
                  jhi = j
               end if
            end if
            loerr = dbnd(j,i) - dist
            if (loerr .gt. 0.0d0) then
               nloerr = nloerr + 1
               rms = rms + loerr**2
               if (loerr .gt. lomax) then
                  lomax = loerr
                  ilo = i
                  jlo = j
               end if
            end if
         end do
      end do
      rms = sqrt(rms/dble(n*(n-1)/2))
c
c     print the maximal and rms bound deviations
c
      leng = trimtext(title)
      write (iout,10)  title(1:leng)
   10 format (/,' Fit to Bounds ',a)
      write (iout,20)  nhierr,npair,nloerr,npair,himax,
     &                 ihi,jhi,lomax,ilo,jlo,rms
   20 format (/,' Num Upper Bound Violations :',4x,i11,'  of ',i12,
     &        /,' Num Lower Bound Violations :',4x,i11,'  of ',i12,
     &        /,' Max Upper Bound Violation :',4x,f12.4,'  at ',2i6,
     &        /,' Max Lower Bound Violation :',4x,f12.4,'  at ',2i6,
     &        /,' RMS Deviation from Bounds :',4x,f12.4)
c
c     search the list of distance restraints for violations
c
      if (ndfix .gt. 0) then
         nloerr = 0
         nhierr = 0
         ilo = 0
         jlo = 0
         ihi = 0
         jhi = 0
         rms = 0.0d0
         himax = 0.0d0
         lomax = 0.0d0
         do k = 1, ndfix
            i = idfix(1,k)
            j = idfix(2,k)
            dist = (x(i)-x(j))**2 + (y(i)-y(j))**2 + (z(i)-z(j))**2
            dist = sqrt(dist)
            if (dist .lt. dfix(2,k)) then
               nloerr = nloerr + 1
               loerr = dfix(2,k) - dist
               rms = rms + loerr**2
               if (loerr .gt. lomax) then
                  lomax = loerr
                  ilo = i
                  jlo = j
               end if
            else if (dist .gt. dfix(3,k)) then
               nhierr = nhierr + 1
               hierr = dist - dfix(3,k)
               rms = rms + hierr**2
               if (hierr .gt. himax) then
                  himax = hierr
                  ihi = i
                  jhi = j
               end if
            end if
         end do
         rms = sqrt(rms/dble(ndfix))
c
c     print total number and rms value of restraint violations
c
         write (iout,30)  nhierr,ndfix,nloerr,ndfix,himax,
     &                    ihi,jhi,lomax,ilo,jlo,rms
   30    format (/,' Num Upper Restraint Violations :',i11,'  of ',i12,
     &           /,' Num Lower Restraint Violations :',i11,'  of ',i12,
     &           /,' Max Upper Restraint Violation :',f12.4,'  at ',2i6,
     &           /,' Max Lower Restraint Violation :',f12.4,'  at ',2i6,
     &           /,' RMS Restraint Dist Violation : ',f12.4)
      end if
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine dmdump  --  final distance and error matrix  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "dmdump" puts the distance matrix of the final structure
c     into the upper half of a matrix, the distance of each atom
c     to the centroid on the diagonal, and the individual terms
c     of the bounds errors into the lower half of the matrix
c
c
      subroutine dmdump (dmd)
      use sizes
      use atoms
      use disgeo
      use iounit
      implicit none
      integer i,j
      real*8 sum,rgsq
      real*8 dist,dist2
      real*8 dmd(n,*)
      character*240 title
c
c
c     store the final distance matrix and bound violations
c
      do i = 1, n
         dmd(i,i) = 0.0d0
      end do
      sum = 0.0d0
      do i = 1, n-1
         do j = i+1, n
            dist2 = (x(i)-x(j))**2 + (y(i)-y(j))**2 + (z(i)-z(j))**2
            sum = sum + dist2
            dmd(i,i) = dmd(i,i) + dist2
            dmd(j,j) = dmd(j,j) + dist2
            dist = sqrt(dist2)
            dmd(i,j) = dist
            if (dist .gt. dbnd(i,j)) then
               dmd(j,i) = dist - dbnd(i,j)
            else if (dist .lt. dbnd(j,i)) then
               dmd(j,i) = dbnd(j,i) - dist
            else
               dmd(j,i) = 0.0d0
            end if
         end do
      end do
c
c     put the distance to the centroid on the diagonal
c
      rgsq = sum / dble(n**2)
      do i = 1, n
         dmd(i,i) = sqrt(dmd(i,i)/dble(n) - rgsq)
      end do
c
c     write out the interatomic distance and error matrices
c
      title = 'Final Dist Matrix Above; DCM on Diag; Error Below :'
      call grafic (n,dmd,title)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  function initerr  --  initial error function and gradient  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "initerr" is the initial error function and derivatives for
c     a distance geometry embedding; it includes components from
c     the local geometry and torsional restraint errors
c
c
      function initerr (xx,g)
      use sizes
      use atoms
      implicit none
      integer i,j,nvar
      real*8 initerr
      real*8 local,locerr
      real*8 torsion,torser
      real*8 xx(*)
      real*8 g(*)
      real*8, allocatable :: derivs(:,:)
c
c
c     translate optimization parameters to atomic coordinates
c
      nvar = 0
      do i = 1, n
         nvar = nvar + 1
         x(i) = xx(nvar)
         nvar = nvar + 1
         y(i) = xx(nvar)
         nvar = nvar + 1
         z(i) = xx(nvar)
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (derivs(3,n))
c
c     zero out the values of the atomic gradient components
c
      do i = 1, n
         do j = 1, 3
            derivs(j,i) = 0.0d0
         end do
      end do
c
c     compute the local geometry and the torsional
c     components of the error function and its gradient
c
      local = locerr (derivs)
      torsion = torser (derivs)
      initerr = local + torsion
c
c     store the atomic gradients as the optimization gradient
c
      nvar = 0
      do i = 1, n
         nvar = nvar + 1
         g(nvar) = derivs(1,i)
         nvar = nvar + 1
         g(nvar) = derivs(2,i)
         nvar = nvar + 1
         g(nvar) = derivs(3,i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (derivs)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  function miderr  --  second error function and gradient  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "miderr" is the secondary error function and derivatives
c     for a distance geometry embedding; it includes components
c     from the distance bounds, local geometry, chirality and
c     torsional restraint errors
c
c
      function miderr (xx,g)
      use sizes
      use atoms
      implicit none
      integer i,j,nvar
      real*8 miderr
      real*8 bounds,bnderr
      real*8 local,locerr
      real*8 chiral,chirer
      real*8 torsion,torser
      real*8 xx(*)
      real*8 g(*)
      real*8, allocatable :: derivs(:,:)
c
c
c     translate optimization parameters to atomic coordinates
c
      nvar = 0
      do i = 1, n
         nvar = nvar + 1
         x(i) = xx(nvar)
         nvar = nvar + 1
         y(i) = xx(nvar)
         nvar = nvar + 1
         z(i) = xx(nvar)
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (derivs(3,n))
c
c     zero out the values of the atomic gradient components
c
      do i = 1, n
         do j = 1, 3
            derivs(j,i) = 0.0d0
         end do
      end do
c
c     compute the local geometry, chirality and torsional
c     components of the error function and its gradient
c
      bounds = bnderr (derivs)
      local = locerr (derivs)
      chiral = chirer (derivs)
      torsion = torser (derivs)
      miderr = bounds + local + chiral + torsion
c
c     store the atomic gradients as the optimization gradient
c
      nvar = 0
      do i = 1, n
         nvar = nvar + 1
         g(nvar) = derivs(1,i)
         nvar = nvar + 1
         g(nvar) = derivs(2,i)
         nvar = nvar + 1
         g(nvar) = derivs(3,i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (derivs)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  function toterr  --  total error function and gradient  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "toterr" is the error function and derivatives for a distance
c     geometry embedding; it includes components from the distance
c     bounds, hard sphere contacts, local geometry, chirality and
c     torsional restraint errors
c
c
      function toterr (xx,g)
      use sizes
      use atoms
      implicit none
      integer i,j,nvar
      real*8 toterr
      real*8 bounds,bnderr
      real*8 local,locerr
      real*8 chiral,chirer
      real*8 torsion,torser
      real*8 contact,vdwerr
      real*8 xx(*)
      real*8 g(*)
      real*8, allocatable :: derivs(:,:)
c
c
c     translate optimization parameters to atomic coordinates
c
      nvar = 0
      do i = 1, n
         nvar = nvar + 1
         x(i) = xx(nvar)
         nvar = nvar + 1
         y(i) = xx(nvar)
         nvar = nvar + 1
         z(i) = xx(nvar)
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (derivs(3,n))
c
c     zero out the values of the atomic gradient components
c
      do i = 1, n
         do j = 1, 3
            derivs(j,i) = 0.0d0
         end do
      end do
c
c     compute the distance bound, vdw, chirality and torsional
c     components of the error function and its gradient
c
      bounds = bnderr (derivs)
      contact = vdwerr (derivs)
      local = locerr (derivs)
      chiral = chirer (derivs)
      torsion = torser (derivs)
      toterr = bounds + contact + local + chiral + torsion
c
c     store the atomic gradients as the optimization gradient
c
      nvar = 0
      do i = 1, n
         nvar = nvar + 1
         g(nvar) = derivs(1,i)
         nvar = nvar + 1
         g(nvar) = derivs(2,i)
         nvar = nvar + 1
         g(nvar) = derivs(3,i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (derivs)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  function bnderr  --  computes total distance bound error  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "bnderr" is the distance bound error function and derivatives;
c     this version implements the original and Havel's normalized
c     lower bound penalty, the normalized version is preferred when
c     lower bounds are small (as with NMR NOE restraints), the
c     original penalty is needed if large lower bounds are present
c
c
      function bnderr (derivs)
      use sizes
      use atoms
      use restrn
      implicit none
      integer i,j,k
      real*8 bnderr,error,cutsq
      real*8 scale,chain,term
      real*8 gap,buffer,weigh
      real*8 dx,dy,dz,gx,gy,gz
      real*8 dstsq,bupsq,blosq
      real*8 derivs(3,*)
c
c
c     zero out the distance bounds error function
c
      bnderr = 0.0d0
      scale = 10.0d0
      cutsq = 40.0d0
c
c     calculate the pairwise distances between atoms
c
      do k = 1, ndfix
         i = idfix(1,k)
         j = idfix(2,k)
         dx = x(i) - x(j)
         dy = y(i) - y(j)
         dz = z(i) - z(j)
c
c     calculate squared actual distance and bound distances;
c     use of a small "buffer" cleans up the final error count
c
         dstsq = dx*dx + dy*dy + dz*dz
         gap = dfix(3,k) - dfix(2,k)
         buffer = 0.05d0 * min(1.0d0,gap)
         blosq = (dfix(2,k) + buffer)**2
         bupsq = (dfix(3,k) - buffer)**2
c
c     error and derivatives for upper bound violation
c
         if (dstsq .gt. bupsq) then
            weigh = scale * dfix(1,k)
            term = (dstsq-bupsq) / bupsq
            chain = 4.0d0 * weigh * term / bupsq
            error = weigh * term**2
            gx = dx * chain
            gy = dy * chain
            gz = dz * chain
            bnderr = bnderr + error
            derivs(1,i) = derivs(1,i) + gx
            derivs(2,i) = derivs(2,i) + gy
            derivs(3,i) = derivs(3,i) + gz
            derivs(1,j) = derivs(1,j) - gx
            derivs(2,j) = derivs(2,j) - gy
            derivs(3,j) = derivs(3,j) - gz
c
c     error and derivatives for lower bound violation
c
         else if (dstsq .lt. blosq) then
            weigh = scale * dfix(1,k)
            if (blosq .gt. cutsq) then
               term = (blosq-dstsq) / dstsq
               chain = -4.0d0 * weigh * term * (blosq/dstsq**2)
            else
               term = (blosq-dstsq) / (blosq+dstsq)
               chain = -8.0d0 * weigh * term * (blosq/(blosq+dstsq)**2)
            end if
            error = weigh * term**2
            gx = dx * chain
            gy = dy * chain
            gz = dz * chain
            bnderr = bnderr + error
            derivs(1,i) = derivs(1,i) + gx
            derivs(2,i) = derivs(2,i) + gy
            derivs(3,i) = derivs(3,i) + gz
            derivs(1,j) = derivs(1,j) - gx
            derivs(2,j) = derivs(2,j) - gy
            derivs(3,j) = derivs(3,j) - gz
         end if
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  function vdwerr  --  computes van der Waals bound error  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "vdwerr" is the hard sphere van der Waals bound error function
c     and derivatives that penalizes close nonbonded contacts,
c     pairwise neighbors are generated via the method of lights
c
c
      function vdwerr (derivs)
      use sizes
      use atoms
      use couple
      use disgeo
      use light
      implicit none
      integer i,j,k,kgy,kgz
      integer, allocatable :: skip(:)
      real*8 vdwerr,error
      real*8 scale,chain,term
      real*8 xi,yi,zi
      real*8 dx,dy,dz,gx,gy,gz
      real*8 dstsq,blosq
      real*8 radi,radsq
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      real*8 derivs(3,*)
      logical unique
c
c
c     zero out the distance van der Waals error function
c
      vdwerr = 0.0d0
      scale = 1.0d0
c
c     perform dynamic allocation of some local arrays
c
      allocate (skip(n))
      allocate (xsort(n))
      allocate (ysort(n))
      allocate (zsort(n))
c
c     transfer coordinates and zero out atoms to be skipped
c
      do i = 1, n
         xsort(i) = x(i)
         ysort(i) = y(i)
         zsort(i) = z(i)
         skip(i) = 0
      end do
c
c     use the method of lights to generate neighbors
c
      unique = .true.
      call lights (vdwmax,n,xsort,ysort,zsort,unique)
c
c     now, loop over all atoms computing the interactions
c
      do i = 1, n
         radi = georad(i)
         xi = xsort(rgx(i))
         yi = ysort(rgy(i))
         zi = zsort(rgz(i))
         do j = 1, n12(i)
            skip(i12(j,i)) = i
         end do
         do j = 1, n13(i)
            skip(i13(j,i)) = i
         end do
         do j = 1, n14(i)
            skip(i14(j,i)) = i
         end do
         do j = kbx(i)+1, kex(i)
            k = locx(j)
            if (skip(k) .eq. i)  goto 10
            kgy = rgy(k)
            if (kgy.lt.kby(i) .or. kgy.gt.key(i))  goto 10
            kgz = rgz(k)
            if (kgz.lt.kbz(i) .or. kgz.gt.kez(i))  goto 10
c
c     calculate squared distances and bounds
c
            dx = xi - xsort(j)
            dy = yi - ysort(kgy)
            dz = zi - zsort(kgz)
            dstsq = dx*dx + dy*dy + dz*dz
            radsq = (radi + georad(k))**2
            blosq = min(dbnd(k,i),dbnd(i,k),radsq)
c
c     error and derivatives for lower bound violation
c
            if (dstsq .lt. blosq) then
               term = (blosq-dstsq) / (blosq+dstsq)
               chain = -8.0d0 * scale * term * (blosq/(blosq+dstsq)**2)
               error = scale * term**2
               gx = dx * chain
               gy = dy * chain
               gz = dz * chain
               vdwerr = vdwerr + error
               derivs(1,i) = derivs(1,i) + gx
               derivs(2,i) = derivs(2,i) + gy
               derivs(3,i) = derivs(3,i) + gz
               derivs(1,k) = derivs(1,k) - gx
               derivs(2,k) = derivs(2,k) - gy
               derivs(3,k) = derivs(3,k) - gz
            end if
   10       continue
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (skip)
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  function locerr  --  computes local geometry error value  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "locerr" is the local geometry error function and derivatives
c     including the 1-2, 1-3 and 1-4 distance bound restraints
c
c
      function locerr (derivs)
      use sizes
      use angbnd
      use atoms
      use bndstr
      use disgeo
      use tors
      implicit none
      integer i,ia,ib,ic,id
      real*8 locerr,error
      real*8 scale,chain,term
      real*8 dx,dy,dz,gx,gy,gz
      real*8 dstsq,bupsq,blosq
      real*8 derivs(3,*)
c
c
c     zero out the local geometry error function
c
      locerr = 0.0d0
      scale = 10.0d0
c
c     calculate the bounds error for bond length distances
c
      do i = 1, nbond
         ia = min(ibnd(1,i),ibnd(2,i))
         ib = max(ibnd(1,i),ibnd(2,i))
         dx = x(ia) - x(ib)
         dy = y(ia) - y(ib)
         dz = z(ia) - z(ib)
         dstsq = dx*dx + dy*dy + dz*dz
         bupsq = dbnd(ia,ib)
         blosq = dbnd(ib,ia)
         if (dstsq .gt. bupsq) then
            term = (dstsq-bupsq) / bupsq
            chain = 4.0d0 * scale * term / bupsq
            error = scale * term**2
            gx = dx * chain
            gy = dy * chain
            gz = dz * chain
            locerr = locerr + error
            derivs(1,ia) = derivs(1,ia) + gx
            derivs(2,ia) = derivs(2,ia) + gy
            derivs(3,ia) = derivs(3,ia) + gz
            derivs(1,ib) = derivs(1,ib) - gx
            derivs(2,ib) = derivs(2,ib) - gy
            derivs(3,ib) = derivs(3,ib) - gz
         else if (dstsq .lt. blosq) then
            term = (blosq-dstsq) / (blosq+dstsq)
            chain = -8.0d0 * scale * term * (blosq/(blosq+dstsq)**2)
            error = scale * term**2
            gx = dx * chain
            gy = dy * chain
            gz = dz * chain
            locerr = locerr + error
            derivs(1,ia) = derivs(1,ia) + gx
            derivs(2,ia) = derivs(2,ia) + gy
            derivs(3,ia) = derivs(3,ia) + gz
            derivs(1,ib) = derivs(1,ib) - gx
            derivs(2,ib) = derivs(2,ib) - gy
            derivs(3,ib) = derivs(3,ib) - gz
         end if
      end do
c
c     calculate the bounds error for the bond angle distances
c
      do i = 1, nangle
         ia = min(iang(1,i),iang(3,i))
         ic = max(iang(1,i),iang(3,i))
         dx = x(ia) - x(ic)
         dy = y(ia) - y(ic)
         dz = z(ia) - z(ic)
         dstsq = dx*dx + dy*dy + dz*dz
         bupsq = dbnd(ia,ic)
         blosq = dbnd(ic,ia)
         if (dstsq .gt. bupsq) then
            term = (dstsq-bupsq) / bupsq
            chain = 4.0d0 * scale * term / bupsq
            error = scale * term**2
            gx = dx * chain
            gy = dy * chain
            gz = dz * chain
            locerr = locerr + error
            derivs(1,ia) = derivs(1,ia) + gx
            derivs(2,ia) = derivs(2,ia) + gy
            derivs(3,ia) = derivs(3,ia) + gz
            derivs(1,ic) = derivs(1,ic) - gx
            derivs(2,ic) = derivs(2,ic) - gy
            derivs(3,ic) = derivs(3,ic) - gz
         else if (dstsq .lt. blosq) then
            term = (blosq-dstsq) / (blosq+dstsq)
            chain = -8.0d0 * scale * term * (blosq/(blosq+dstsq)**2)
            error = scale * term**2
            gx = dx * chain
            gy = dy * chain
            gz = dz * chain
            locerr = locerr + error
            derivs(1,ia) = derivs(1,ia) + gx
            derivs(2,ia) = derivs(2,ia) + gy
            derivs(3,ia) = derivs(3,ia) + gz
            derivs(1,ic) = derivs(1,ic) - gx
            derivs(2,ic) = derivs(2,ic) - gy
            derivs(3,ic) = derivs(3,ic) - gz
         end if
      end do
c
c     calculate the bounds error for the torsion angle distances
c
      do i = 1, ntors
         ia = min(itors(1,i),itors(4,i))
         id = max(itors(1,i),itors(4,i))
         dx = x(ia) - x(id)
         dy = y(ia) - y(id)
         dz = z(ia) - z(id)
         dstsq = dx*dx + dy*dy + dz*dz
         bupsq = dbnd(ia,id)
         blosq = dbnd(id,ia)
         if (dstsq .gt. bupsq) then
            term = (dstsq-bupsq) / bupsq
            chain = 4.0d0 * scale * term / bupsq
            error = scale * term**2
            gx = dx * chain
            gy = dy * chain
            gz = dz * chain
            locerr = locerr + error
            derivs(1,ia) = derivs(1,ia) + gx
            derivs(2,ia) = derivs(2,ia) + gy
            derivs(3,ia) = derivs(3,ia) + gz
            derivs(1,id) = derivs(1,id) - gx
            derivs(2,id) = derivs(2,id) - gy
            derivs(3,id) = derivs(3,id) - gz
         else if (dstsq .lt. blosq) then
            term = (blosq-dstsq) / (blosq+dstsq)
            chain = -8.0d0 * scale * term * (blosq/(blosq+dstsq)**2)
            error = scale * term**2
            gx = dx * chain
            gy = dy * chain
            gz = dz * chain
            locerr = locerr + error
            derivs(1,ia) = derivs(1,ia) + gx
            derivs(2,ia) = derivs(2,ia) + gy
            derivs(3,ia) = derivs(3,ia) + gz
            derivs(1,id) = derivs(1,id) - gx
            derivs(2,id) = derivs(2,id) - gy
            derivs(3,id) = derivs(3,id) - gz
         end if
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  function chirer  --  computes chirality error function  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "chirer" computes the chirality error and its derivatives
c     with respect to atomic Cartesian coordinates as a sum the
c     squares of deviations of chiral volumes from target values
c
c
      function chirer (derivs)
      use sizes
      use atoms
      use restrn
      implicit none
      integer i,ia,ib,ic,id
      real*8 chirer,error,scale
      real*8 vol,dt,dedt
      real*8 xad,yad,zad
      real*8 xbd,ybd,zbd
      real*8 xcd,ycd,zcd
      real*8 c1,c2,c3
      real*8 dedxia,dedyia,dedzia
      real*8 dedxib,dedyib,dedzib
      real*8 dedxic,dedyic,dedzic
      real*8 dedxid,dedyid,dedzid
      real*8 derivs(3,*)
c
c
c     zero the chirality restraint error function
c
      chirer = 0.0d0
      scale = 0.1d0
c
c     find signed volume value and compute chirality error
c
      do i = 1, nchir
         ia = ichir(1,i)
         ib = ichir(2,i)
         ic = ichir(3,i)
         id = ichir(4,i)
         xad = x(ia) - x(id)
         yad = y(ia) - y(id)
         zad = z(ia) - z(id)
         xbd = x(ib) - x(id)
         ybd = y(ib) - y(id)
         zbd = z(ib) - z(id)
         xcd = x(ic) - x(id)
         ycd = y(ic) - y(id)
         zcd = z(ic) - z(id)
         c1 = ybd*zcd - zbd*ycd
         c2 = ycd*zad - zcd*yad
         c3 = yad*zbd - zad*ybd
         vol = xad*c1 + xbd*c2 + xcd*c3
         dt = vol - chir(2,i)
         error = scale * dt**2
         dedt = 2.0d0 * scale * dt
c
c     compute derivative components for this interaction
c
         dedxia = dedt * (ybd*zcd - zbd*ycd)
         dedyia = dedt * (zbd*xcd - xbd*zcd)
         dedzia = dedt * (xbd*ycd - ybd*xcd)
         dedxib = dedt * (zad*ycd - yad*zcd)
         dedyib = dedt * (xad*zcd - zad*xcd)
         dedzib = dedt * (yad*xcd - xad*ycd)
         dedxic = dedt * (yad*zbd - zad*ybd)
         dedyic = dedt * (zad*xbd - xad*zbd)
         dedzic = dedt * (xad*ybd - yad*xbd)
         dedxid = -dedxia - dedxib - dedxic
         dedyid = -dedyia - dedyib - dedyic
         dedzid = -dedzia - dedzib - dedzic
c
c     increment the chirality restraint error and derivatives
c
         chirer = chirer + error
         derivs(1,ia) = derivs(1,ia) + dedxia
         derivs(2,ia) = derivs(2,ia) + dedyia
         derivs(3,ia) = derivs(3,ia) + dedzia
         derivs(1,ib) = derivs(1,ib) + dedxib
         derivs(2,ib) = derivs(2,ib) + dedyib
         derivs(3,ib) = derivs(3,ib) + dedzib
         derivs(1,ic) = derivs(1,ic) + dedxic
         derivs(2,ic) = derivs(2,ic) + dedyic
         derivs(3,ic) = derivs(3,ic) + dedzic
         derivs(1,id) = derivs(1,id) + dedxid
         derivs(2,id) = derivs(2,id) + dedyid
         derivs(3,id) = derivs(3,id) + dedzid
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  function torser  --  computes torsional error function  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "torser" computes the torsional error function and its first
c     derivatives with respect to the atomic Cartesian coordinates
c     based on the deviation of specified torsional angles from
c     desired values, the contained bond angles are also restrained
c     to avoid a numerical instability
c
c
      function torser (derivs)
      use sizes
      use atoms
      use couple
      use disgeo
      use math
      use refer
      use restrn
      implicit none
      integer i,j,k,iref
      integer ia,ib,ic,id
      real*8 torser,error,force
      real*8 dt,deddt,dedphi
      real*8 angle,target,scale
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xria,yria,zria
      real*8 xrib,yrib,zrib
      real*8 xric,yric,zric
      real*8 xrid,yrid,zrid
      real*8 rba,rcb,rdc
      real*8 dot,cosine,sine
      real*8 rrba,rrcb,rrdc
      real*8 rrca,rrdb
      real*8 bndfac,angfac
      real*8 xp,yp,zp,rp
      real*8 terma,termb
      real*8 termc,termd
      real*8 angmax,angmin
      real*8 cosmax,cosmin
      real*8 bamax,bamin
      real*8 cbmax,cbmin
      real*8 dcmax,dcmin
      real*8 camax,camin
      real*8 dbmax,dbmin
      real*8 xba,yba,zba
      real*8 xcb,ycb,zcb
      real*8 xdc,ydc,zdc
      real*8 xca,yca,zca
      real*8 xdb,ydb,zdb
      real*8 xt,yt,zt,rt2
      real*8 xu,yu,zu,ru2
      real*8 xtu,ytu,ztu,rtru
      real*8 tf1,tf2,t1,t2
      real*8 dedxt,dedyt,dedzt
      real*8 dedxu,dedyu,dedzu
      real*8 dedxia,dedyia,dedzia
      real*8 dedxib,dedyib,dedzib
      real*8 dedxic,dedyic,dedzic
      real*8 dedxid,dedyid,dedzid
      real*8 derivs(3,*)
      logical bonded
c
c
c     zero the torsional restraint error function
c
      torser = 0.0d0
      scale = 0.01d0
c
c     compute error value and derivs for torsional restraints
c
      do i = 1, ntfix
         ia = itfix(1,i)
         ib = itfix(2,i)
         ic = itfix(3,i)
         id = itfix(4,i)
         xia = x(ia)
         yia = y(ia)
         zia = z(ia)
         xib = x(ib)
         yib = y(ib)
         zib = z(ib)
         xic = x(ic)
         yic = y(ic)
         zic = z(ic)
         xid = x(id)
         yid = y(id)
         zid = z(id)
         xba = xib - xia
         yba = yib - yia
         zba = zib - zia
         xcb = xic - xib
         ycb = yic - yib
         zcb = zic - zib
         xdc = xid - xic
         ydc = yid - yic
         zdc = zid - zic
c
c     find the actual distances between the four atoms
c
         rba = sqrt(xba*xba + yba*yba + zba*zba)
         rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
         rdc = sqrt(xdc*xdc + ydc*ydc + zdc*zdc)
c
c     see if the torsion involves four directly bonded atoms
c
         k = 0
         do j = 1, n12(ia)
            if (i12(j,ia) .eq. ib)  k = k + 1
         end do
         do j = 1, n12(ib)
            if (i12(j,ib) .eq. ic)  k = k + 1
         end do
         do j = 1, n12(ic)
            if (i12(j,ic) .eq. id)  k = k + 1
         end do
         if (k .eq. 3) then
            bonded = .true.
         else
            bonded = .false.
         end if
c
c     get maximum and minimum distances from distance matrix
c
         if (bonded) then
            bamax = sqrt(dbnd(min(ib,ia),max(ib,ia)))
            bamin = sqrt(dbnd(max(ib,ia),min(ib,ia)))
            cbmax = sqrt(dbnd(min(ic,ib),max(ic,ib)))
            cbmin = sqrt(dbnd(max(ic,ib),min(ic,ib)))
            dcmax = sqrt(dbnd(min(id,ic),max(id,ic)))
            dcmin = sqrt(dbnd(max(id,ic),min(id,ic)))
            camax = sqrt(dbnd(min(ic,ia),max(ic,ia)))
            camin = sqrt(dbnd(max(ic,ia),min(ic,ia)))
            dbmax = sqrt(dbnd(min(id,ib),max(id,ib)))
            dbmin = sqrt(dbnd(max(id,ib),min(id,ib)))
c
c     get maximum and minimum distances from input structure
c
         else
            iref = 1
            xria = xref(ia,iref)
            yria = yref(ia,iref)
            zria = zref(ia,iref)
            xrib = xref(ib,iref)
            yrib = yref(ib,iref)
            zrib = zref(ib,iref)
            xric = xref(ic,iref)
            yric = yref(ic,iref)
            zric = zref(ic,iref)
            xrid = xref(id,iref)
            yrid = yref(id,iref)
            zrid = zref(id,iref)
            rrba = sqrt((xrib-xria)**2+(yrib-yria)**2+(zrib-zria)**2)
            rrcb = sqrt((xric-xrib)**2+(yric-yrib)**2+(zric-zrib)**2)
            rrdc = sqrt((xrid-xric)**2+(yrid-yric)**2+(zrid-zric)**2)
            rrca = sqrt((xric-xria)**2+(yric-yria)**2+(zric-zria)**2)
            rrdb = sqrt((xrid-xrib)**2+(yrid-yrib)**2+(zrid-zrib)**2)
            bndfac = 0.05d0
            angfac = 0.05d0
            bamax = (1.0d0 + bndfac) * rrba
            bamin = (1.0d0 - bndfac) * rrba
            cbmax = (1.0d0 + bndfac) * rrcb
            cbmin = (1.0d0 - bndfac) * rrcb
            dcmax = (1.0d0 + bndfac) * rrdc
            dcmin = (1.0d0 - bndfac) * rrdc
            camax = (1.0d0 + angfac) * rrca
            camin = (1.0d0 - angfac) * rrca
            dbmax = (1.0d0 + angfac) * rrdb
            dbmin = (1.0d0 - angfac) * rrdb
         end if
c
c     compute the ia-ib-ic bond angle and any error
c
         dot = xba*xcb + yba*ycb + zba*zcb
         cosine = -dot / (rba*rcb)
         cosine = min(1.0d0,max(-1.0d0,cosine))
         angle = radian * acos(cosine)
         cosmax = (bamin**2+cbmin**2-camax**2) / (2.0d0*bamin*cbmin)
         cosmax = min(1.0d0,max(-1.0d0,cosmax))
         angmax = radian * acos(cosmax)
         cosmin = (bamax**2+cbmax**2-camin**2) / (2.0d0*bamax*cbmax)
         cosmin = min(1.0d0,max(-1.0d0,cosmin))
         angmin = radian * acos(cosmin)
         if (angle .gt. angmax) then
            dt = angle - angmax
         else if (angle .lt. angmin) then
            dt = angle - angmin
         else
            dt = 0.0d0
         end if
         error = scale * dt**2
         deddt = 2.0d0 * radian * scale * dt
c
c     compute derivative components for this interaction
c
         xp = zcb*yba - ycb*zba
         yp = xcb*zba - zcb*xba
         zp = ycb*xba - xcb*yba
         rp = sqrt(xp*xp + yp*yp + zp*zp)
         if (rp .ne. 0.0d0) then
            terma = -deddt / (rba*rba*rp)
            termc = deddt / (rcb*rcb*rp)
            dedxia = terma * (zba*yp-yba*zp)
            dedyia = terma * (xba*zp-zba*xp)
            dedzia = terma * (yba*xp-xba*yp)
            dedxic = termc * (ycb*zp-zcb*yp)
            dedyic = termc * (zcb*xp-xcb*zp)
            dedzic = termc * (xcb*yp-ycb*xp)
            dedxib = -dedxia - dedxic
            dedyib = -dedyia - dedyic
            dedzib = -dedzia - dedzic
c
c     increment the bond angle restraint error and derivatives
c
            torser = torser + error
            derivs(1,ia) = derivs(1,ia) + dedxia
            derivs(2,ia) = derivs(2,ia) + dedyia
            derivs(3,ia) = derivs(3,ia) + dedzia
            derivs(1,ib) = derivs(1,ib) + dedxib
            derivs(2,ib) = derivs(2,ib) + dedyib
            derivs(3,ib) = derivs(3,ib) + dedzib
            derivs(1,ic) = derivs(1,ic) + dedxic
            derivs(2,ic) = derivs(2,ic) + dedyic
            derivs(3,ic) = derivs(3,ic) + dedzic
         end if
c
c     compute the ib-ic-id bond angle and any error
c
         dot = xdc*xcb + ydc*ycb + zdc*zcb
         cosine = -dot / (rdc*rcb)
         cosine = min(1.0d0,max(-1.0d0,cosine))
         angle = radian * acos(cosine)
         cosmax = (dcmin**2+cbmin**2-dbmax**2) / (2.0d0*dcmin*cbmin)
         cosmax = min(1.0d0,max(-1.0d0,cosmax))
         angmax = radian * acos(cosmax)
         cosmin = (dcmax**2+cbmax**2-dbmin**2) / (2.0d0*dcmax*cbmax)
         cosmax = min(1.0d0,max(-1.0d0,cosmin))
         angmin = radian * acos(cosmin)
         if (angle .gt. angmax) then
            dt = angle - angmax
         else if (angle .lt. angmin) then
            dt = angle - angmin
         else
            dt = 0.0d0
         end if
         error = scale * dt**2
         deddt = 2.0d0 * radian * scale * dt
c
c     compute derivative components for this interaction
c
         xp = zdc*ycb - ydc*zcb
         yp = xdc*zcb - zdc*xcb
         zp = ydc*xcb - xdc*ycb
         rp = sqrt(xp*xp + yp*yp + zp*zp)
         if (rp .ne. 0.0d0) then
            termb = -deddt / (rcb*rcb*rp)
            termd = deddt / (rdc*rdc*rp)
            dedxib = termb * (zcb*yp-ycb*zp)
            dedyib = termb * (xcb*zp-zcb*xp)
            dedzib = termb * (ycb*xp-xcb*yp)
            dedxid = termd * (ydc*zp-zdc*yp)
            dedyid = termd * (zdc*xp-xdc*zp)
            dedzid = termd * (xdc*yp-ydc*xp)
            dedxic = -dedxib - dedxid
            dedyic = -dedyib - dedyid
            dedzic = -dedzib - dedzid
c
c     increment the bond angle restraint error and derivatives
c
            torser = torser + error
            derivs(1,ib) = derivs(1,ib) + dedxib
            derivs(2,ib) = derivs(2,ib) + dedyib
            derivs(3,ib) = derivs(3,ib) + dedzib
            derivs(1,ic) = derivs(1,ic) + dedxic
            derivs(2,ic) = derivs(2,ic) + dedyic
            derivs(3,ic) = derivs(3,ic) + dedzic
            derivs(1,id) = derivs(1,id) + dedxid
            derivs(2,id) = derivs(2,id) + dedyid
            derivs(3,id) = derivs(3,id) + dedzid
         end if
c
c     compute the value of the ia-ib-ic-id torsional angle
c
         xt = yba*zcb - ycb*zba
         yt = zba*xcb - zcb*xba
         zt = xba*ycb - xcb*yba
         xu = ycb*zdc - ydc*zcb
         yu = zcb*xdc - zdc*xcb
         zu = xcb*ydc - xdc*ycb
         xtu = yt*zu - yu*zt
         ytu = zt*xu - zu*xt
         ztu = xt*yu - xu*yt
         rt2 = xt*xt + yt*yt + zt*zt
         ru2 = xu*xu + yu*yu + zu*zu
         rtru = sqrt(rt2 * ru2)
         if (rtru .ne. 0.0d0) then
            rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
            cosine = (xt*xu + yt*yu + zt*zu) / rtru
            sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
            cosine = min(1.0d0,max(-1.0d0,cosine))
            angle = radian * acos(cosine)
            if (sine .lt. 0.0d0)  angle = -angle
c
c     calculate the torsional restraint error for this angle
c
            force = tfix(1,i)
            tf1 = tfix(2,i)
            tf2 = tfix(3,i)
            if (angle.gt.tf1 .and. angle.lt.tf2) then
               target = angle
            else if (angle.gt.tf1 .and. tf1.gt.tf2) then
               target = angle
            else if (angle.lt.tf2 .and. tf1.gt.tf2) then
               target = angle
            else
               t1 = angle - tf1
               t2 = angle - tf2
               if (t1 .gt. 180.0d0) then
                  t1 = t1 - 360.0d0
               else if (t1 .lt. -180.0d0) then
                  t1 = t1 + 360.0d0
               end if
               if (t2 .gt. 180.0d0) then
                  t2 = t2 - 360.0d0
               else if (t2 .lt. -180.0d0) then
                  t2 = t2 + 360.0d0
               end if
               if (abs(t1) .lt. abs(t2)) then
                  target = tf1
               else
                  target = tf2
               end if
            end if
            dt = angle - target
            if (dt .gt. 180.0d0) then
               dt = dt - 360.0d0
            else if (dt .lt. -180.0d0) then
               dt = dt + 360.0d0
            end if
            error = scale * force * dt**2
            dedphi = 2.0d0 * radian * scale * force * dt
c
c     chain rule terms for first derivative components
c
            xca = xic - xia
            yca = yic - yia
            zca = zic - zia
            xdb = xid - xib
            ydb = yid - yib
            zdb = zid - zib
            dedxt = dedphi * (yt*zcb - ycb*zt) / (rt2*rcb)
            dedyt = dedphi * (zt*xcb - zcb*xt) / (rt2*rcb)
            dedzt = dedphi * (xt*ycb - xcb*yt) / (rt2*rcb)
            dedxu = -dedphi * (yu*zcb - ycb*zu) / (ru2*rcb)
            dedyu = -dedphi * (zu*xcb - zcb*xu) / (ru2*rcb)
            dedzu = -dedphi * (xu*ycb - xcb*yu) / (ru2*rcb)
c
c     compute first derivative components for torsion angle
c
            dedxia = zcb*dedyt - ycb*dedzt
            dedyia = xcb*dedzt - zcb*dedxt
            dedzia = ycb*dedxt - xcb*dedyt
            dedxib = yca*dedzt - zca*dedyt + zdc*dedyu - ydc*dedzu
            dedyib = zca*dedxt - xca*dedzt + xdc*dedzu - zdc*dedxu
            dedzib = xca*dedyt - yca*dedxt + ydc*dedxu - xdc*dedyu
            dedxic = zba*dedyt - yba*dedzt + ydb*dedzu - zdb*dedyu
            dedyic = xba*dedzt - zba*dedxt + zdb*dedxu - xdb*dedzu
            dedzic = yba*dedxt - xba*dedyt + xdb*dedyu - ydb*dedxu
            dedxid = zcb*dedyu - ycb*dedzu
            dedyid = xcb*dedzu - zcb*dedxu
            dedzid = ycb*dedxu - xcb*dedyu
c
c     increment the torsional restraint error and derivatives
c
            torser = torser + error
            derivs(1,ia) = derivs(1,ia) + dedxia
            derivs(2,ia) = derivs(2,ia) + dedyia
            derivs(3,ia) = derivs(3,ia) + dedzia
            derivs(1,ib) = derivs(1,ib) + dedxib
            derivs(2,ib) = derivs(2,ib) + dedyib
            derivs(3,ib) = derivs(3,ib) + dedzib
            derivs(1,ic) = derivs(1,ic) + dedxic
            derivs(2,ic) = derivs(2,ic) + dedyic
            derivs(3,ic) = derivs(3,ic) + dedzic
            derivs(1,id) = derivs(1,id) + dedxid
            derivs(2,id) = derivs(2,id) + dedyid
            derivs(3,id) = derivs(3,id) + dedzid
         end if
      end do
      return
      end
