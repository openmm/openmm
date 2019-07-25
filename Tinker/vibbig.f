c
c
c     #################################################################
c     ##  COPYRIGHT (C) 2007 by Alexey Kaledin & Jay William Ponder  ##
c     ##                     All Rights Reserved                     ##
c     #################################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  program vibbig  --  block iterative vibrational analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "vibbig" performs large-scale vibrational mode analysis using
c     only vector storage and gradient evaluations; preconditioning
c     is via an approximate inverse from a block diagonal Hessian,
c     and a sliding block method is used to converge any number of
c     eigenvectors starting from either lowest or highest frequency
c
c     literature references:
c
c     C. Murray, S. C. Racine and E. R. Davidson, "Improved Algorithms
c     for the Lowest Few Eigenvalues and Associated Eigenvectors of
c     Large Matrices", Journal of Computational Physics, 103, 382-389
c     (1992)
c
c     A. L. Kaledin, "Gradient-Based Direct Normal-Mode Analysis",
c     Journal of Chemical Physics, 122, 184106 (2005)
c
c     A. L. Kaledin, M. Kaledin and J. M. Bowman, "All-Atom Calculation
c     of the Normal Modes of Bacteriorhodopsin Using a Sliding Block
c     Iterative Diagonalization Method", Journal of Chemical Theory
c     and Computation, 2, 166-174 (2006)
c
c
      program vibbig
      use sizes
      use atomid
      use atoms
      use files
      use inform
      use iounit
      use keys
      use units
      use vibs
      implicit none
      integer i,j,k,ii,next
      integer i1,i2,k0,k1,k2
      integer ivib,ivb1,ivb2
      integer iblock,iconv
      integer iter,idump
      integer nvar,nblk
      integer nroot,nbasis
      integer np,npair
      integer nlock,nconv
      integer irange,ifactor
      integer maxroot,maxiter
      integer maxhess
      integer freeunit
      integer, allocatable :: iblk(:)
      real*8 fmax,funit
      real*8 wtol,factor
      real*8 size,sizmax
      real*8 space,sum
      real*8 dfreq,rnorm,rcomp
      real*8 ratio,shift
      real*8 uku_min,uku_max
      real*8, allocatable :: xe(:)
      real*8, allocatable :: xm(:)
      real*8, allocatable :: p(:)
      real*8, allocatable :: pk(:)
      real*8, allocatable :: hmin(:)
      real*8, allocatable :: uku(:)
      real*8, allocatable :: uku0(:)
      real*8, allocatable :: uu(:)
      real*8, allocatable :: freq(:)
      real*8, allocatable :: freqold(:)
      real*8, allocatable :: tmp1(:)
      real*8, allocatable :: tmp2(:)
      real*8, allocatable :: u(:,:)
      real*8, allocatable :: ur(:,:)
      real*8, allocatable :: h(:,:)
      real*8, allocatable :: c(:,:)
      character*1 answer
      character*20 keyword
      character*240 record
      character*240 string
      character*240 datafile
      character*240 blockfile
      logical exist,restart
      logical header,done
c
c
c     set up the structure and mechanics calculation
c
      call initial
      call getxyz
      call mechanic
c
c     set default parameters for the normal mode computation
c
      nvar = 3 * n
      maxroot = 50
      np = 6
      iter = 0
      idump = 10
      maxhess = nvar * (nvar-1) / 2
      maxiter = 100000
      wtol = 0.00001d0
      sizmax = 500.0d0
      header = .true.
c
c     search the keywords for normal mode parameters
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:8) .eq. 'MAXITER ') then
            read (string,*,err=10,end=10)  maxiter
         else if (keyword(1:6) .eq. 'IDUMP ') then
            read (string,*,err=10,end=10)  idump
         else if (keyword(1:10) .eq. 'VIB-ROOTS ') then
            read (string,*,err=10,end=10)  nroot
            nroot = min(nroot,maxroot)
         else if (keyword(1:14) .eq. 'VIB-TOLERANCE ') then
            read (string,*,err=10,end=10)  wtol
         end if
   10    continue
      end do
c
c     find either the lowest or highest normal modes
c
      factor = 1.0d0
      call nextarg (answer,exist)
      if (.not. exist) then
         answer = 'L'
         write (iout,20)  answer
   20    format (/,' Start at Lowest or Highest Frequency',
     &              ' Normal Mode [',a1,'] :  ',$)
         read (input,30)  record
   30    format (a240)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .eq. 'H')  factor = -1.0d0
c
c     find cutoff value for desired extreme frequency
c
      fmax = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=40,end=40)  fmax
   40 continue
      if (fmax .le. 0.0d0) then
         write (iout,50)
   50    format (/,' Enter Desired Frequency Cutoff in cm-1',
     &              ' [0.0] :  ',$)
         read (input,60)  fmax
   60    format (f20.0)
      end if
      if (fmax .le. 0.0d0)  fmax = 0.0d0
c
c     set default values for some additional parameters
c
      funit = factor * efreq * emass
      ifactor = int(factor)
      irange = (nvar-np+1) * max((1-ifactor)/2,0)
      npair = 2 * nroot
      nbasis = 3 * nroot
c
c     open or create eigenvector file for use during restarts
c
      ivb1 = freeunit ()
      datafile = filename(1:leng)//'.vb1'
      call version (datafile,'old')
      inquire (file=datafile,exist=exist)
      if (exist) then
         open (unit=ivb1,file=datafile,status='old',form='unformatted')
      else
         open (unit=ivb1,file=datafile,status='new',form='unformatted')
      end if
c
c     open or create basis vector file for use during restarts
c
      ivb2 = freeunit ()
      datafile = filename(1:leng)//'.vb2'
      call version (datafile,'old')
      inquire (file=datafile,exist=exist)
      if (exist) then
         restart = .true.
         open (unit=ivb2,file=datafile,status='old',form='unformatted')
      else
         restart = .false.
         open (unit=ivb2,file=datafile,status='new',form='unformatted')
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (iblk(n))
      allocate (xe(nvar))
      allocate (xm(nvar))
      allocate (p(nvar))
      allocate (pk(nvar))
      allocate (hmin(nvar))
      allocate (uku(nvar))
      allocate (uku0(nvar))
      allocate (uu(maxhess))
      allocate (u(nvar,6))
      allocate (ur(nvar,3))
c
c     perform dynamic allocation of some global arrays
c
      allocate (phi(nvar,nbasis))
      allocate (phik(nvar,nbasis))
      allocate (pwork(nvar,nbasis))
c
c     store a coordinate vector for each atom
c
      do i = 1, n
         xe(3*i-2) = x(i) / bohr
         xe(3*i-1) = y(i) / bohr
         xe(3*i) = z(i) / bohr
      end do
c
c     store atomic mass for each coordinate component
c
      k = 0
      do i = 1, n
         mass(i) = mass(i) * emass
         do j = 1, 3
            k = k + 1
            xm(k) = mass(i)
         end do
      end do
c
c     remove pure translational and rotational modes
c
      call trbasis (nvar,np,xe,u,ur)
c
c     set number and size of blocks based on storage space
c
      space = 0.9d0 * dble(maxhess)
      do i = 1, n
         size = 9.0d0 * (dble(n))**2 / dble(i)
         if (size .lt. space) then
            nblk = i
            goto 70
         end if
      end do
   70 continue
      nblk = max(3,nblk)
      size = dble(n) / dble(nblk)
      size = min(size,sizmax)
      do i = 1, nblk
         iblk(i) = nint(dble(i)*size)
      end do
      do i = nblk, 2, -1
         iblk(i) = iblk(i) - iblk(i-1)
      end do
c
c     get number and size of blocks from an external file
c
      iblock = freeunit ()
      blockfile = filename(1:leng)//'.blk'
      call version (blockfile,'old')
      inquire (file=blockfile,exist=exist)
      if (exist) then
         open (unit=iblock,file=blockfile,status='old')
         i = 0
         do while (.true.)
            i = i + 1
            read (iblock,*,err=80,end=80)  iblk(i)
         end do
   80    continue
         nblk = i - 1
         close (unit=iblock)
      end if
c
c     print info about the atom blocks and preconditioning
c
      write (iout,90)
   90 format (/,' Atom Blocks Used to Subdivide the System :',/)
      k = 0
      do i = 1, nblk
         write (iout,100)  i,iblk(i),k+1,k+iblk(i)
  100    format (' Block :',i7,9x,'Size :',i7,9x,'Atoms :',i7,'  to',i7)
         k = k + iblk(i)
      end do
      k = 0
      do i = 1, nblk
         k = k + 9*iblk(i)**2
      end do
      write (iout,110)  k
  110 format (/,' Storage for Preconditioning Array :',5x,i12)
c
c     determine number of prior modes available at restart
c
      nlock = 0
      do while (.true.)
         read (ivb1,err=120,end=120)  (p(k),k=1,nvar)
         nlock = nlock + 1
      end do
  120 continue
      rewind (unit=ivb1)
      if (nlock .ne. 0) then
         write (iout,130)  nlock
  130    format (/,' Prior Normal Modes Available at Restart :',i11)
      end if
      nconv = nlock
c
c     compute and diagonalize the Hessian for each block
c
      k0 = 0
      i1 = 1
      do i = 1, nblk
         if (i .gt. 1) then
            k0 = k0 + 9*iblk(i-1)**2
            i1 = i1 + iblk(i-1)
         end if
         i2 = i1 + iblk(i) - 1
         k1 = 3*i1 - 2
         k2 = 3*i2
         call hessblk (mass,k0,i1,i2,uu)
         call diagblk (k0,k1,3*iblk(i),uu,uku)
      end do
c
c     use negative of eigenvalues if doing high frequencies
c
      do k = 1, nvar
         uku(k) = factor * uku(k)
         uku0(k) = uku(k)
      end do
      uku_max = uku(1)
      uku_min = uku(1)
      do k = 2, nvar
         if (uku(k) .gt. uku_max)  uku_max = uku(k)
         if (uku(k) .lt. uku_min)  uku_min = uku(k)
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (freq(nbasis))
      allocate (freqold(nbasis))
      allocate (tmp1(nbasis))
      allocate (tmp2(nbasis))
      allocate (h(nbasis,nbasis))
      allocate (c(nbasis,nbasis))
c
c     if restarting, read trial vectors and estimate eigenvalues
c
      if (restart) then
         do i = 1, npair
            read (ivb2)  (phi(k,i),k=1,nvar)
            read (ivb2)  (phik(k,i),k=1,nvar)
         end do
         do i = 1, nroot
            h(i,i) = 0.0d0
            do k = 1, nvar
               h(i,i) = h(i,i) + phik(k,i)*phi(k,i)
            end do
            freqold(i) = sign(1.0d0,h(i,i)) * sqrt(abs(h(i,i)))
         end do
         goto 140
      end if
c
c     if not restarting, generate initial guess eigenvectors
c
      do i = 1, nroot
         call trigger (nvar,nbasis,np,ifactor,nblk,iblk,u,uu,p)
         do k = 1, nvar
            phi(k,i) = p(k)
         end do
      end do
c
c     project out locked roots from components of phi
c
      call project (nvar,nconv,ivb1,nroot,0)
      call projectk (nvar,nconv,ivb1,nroot,0)
c
c     reload and make vector orthonormal to existing basis
c
      do i = 1, nroot
         do k = 1, nvar
            p(k) = phi(k,i)
         end do
         if (i .eq. 1) then
            sum = 0.0d0
            do k = 1, nvar
               sum = sum + p(k)*p(k)
            end do
            sum = sqrt(sum)
            do k = 1, nvar
               p(k) = p(k) / sum
            end do
         else
            call gsort (nvar,i-1,p)
         end if
         do k = 1, nvar
            phi(k,i) = p(k)
         end do
      end do
c
c     store K on phi
c
      do i = 1, nroot
         do k = 1, nvar
            p(k) = phi(k,i)
         end do
         call konvec (nvar,xm,xe,p,pk)
         do k = 1, nvar
            phik(k,i) = factor * pk(k)
         end do
      end do
c
c     make nroot-by-nroot CI matrix
c
      do i = 1, nroot
         do j = i, nroot
            h(i,j) = 0.0d0
            do k = 1, nvar
               h(i,j) = h(i,j) + phik(k,i)*phi(k,j)
            end do
            h(j,i) = h(i,j)
         end do
      end do
c
c     diagonalize and use first nroot solutions as starting basis
c
      call transform (nroot,nbasis,h,c)
c
c     fill up arrays
c
      do k = 1, nvar
         do j = 1, nroot
            tmp1(j) = 0.0d0
            tmp2(j) = 0.0d0
            do i = 1, nroot
               tmp1(j) = tmp1(j) + c(i,j)*phi(k,i)
               tmp2(j) = tmp2(j) + c(i,j)*phik(k,i)
            end do
         end do
         do j = 1, nroot
            phi(k,j) = tmp1(j)
            phik(k,j) = tmp2(j)
         end do
      end do
c
c     residues of guesses
c
      do i = 1, nroot
         freq(i) = funit * sign(1.0d0,h(i,i)) * sqrt(abs(h(i,i)))
         freqold(i) = freq(i)
         do k = 1, nvar
            pk(k) = phik(k,i) - h(i,i)*phi(k,i)
         end do
c
c     use Davidson preconditioner if finding low frequencies
c
         if (factor .gt. 0.0d0) then
            call preconblk (nvar,nblk,iblk,uku,uu,h(i,i),hmin(i),pk)
         end if
c
c     project residual onto P-space
c
         call qonvec (nvar,np,u,pk,p)
         do k = 1, nvar
            phi(k,i+nroot) = p(k)
         end do
      end do
c
c     project out locked roots from components of phi
c
      call project (nvar,nconv,ivb1,nroot,nroot)
      call projectk (nvar,nconv,ivb1,nroot,nroot)
c
c     reload and make vector orthonormal to existing basis
c
      do i = 1, nroot
         do k = 1, nvar
            p(k) = phi(k,i+nroot)
         end do
         call gsort (nvar,nroot+i-1,p)
         do k = 1, nvar
            phi(k,i+nroot) = p(k)
         end do
      end do
c
c     store K on phi
c
      do i = 1, nroot
         do k = 1, nvar
            p(k) = phi(k,i+nroot)
         end do
         call konvec (nvar,xm,xe,p,pk)
         do k = 1, nvar
            phik(k,i+nroot) = factor * pk(k)
         end do
      end do
c
c     make npair-by-npair CI matrix
c
      do i = 1, npair
         do j = i, npair
            h(i,j) = 0.0d0
            do k = 1, nvar
               h(i,j) = h(i,j) + phik(k,i)*phi(k,j)
            end do
            h(j,i) = h(i,j)
         end do
      end do
c
c     diagonalize and use first nroot solutions as new guess
c
      call transform (npair,nbasis,h,c)
      do k = 1, nvar
         do j = 1, nroot
            tmp1(j) = 0.0d0
            tmp2(j) = 0.0d0
            do i = 1, npair
               tmp1(j) = tmp1(j) + c(i,j)*phi(k,i)
               tmp2(j) = tmp2(j) + c(i,j)*phik(k,i)
            end do
         end do
c
c     old solution fills up 2nd block
c
         do j = 1, nroot
            phi(k,j+nroot) = phi(k,j)
            phik(k,j+nroot) = phik(k,j)
         end do
c
c     new solution fills up 1st block
c
         do j = 1, nroot
            phi(k,j) = tmp1(j)
            phik(k,j) = tmp2(j)
         end do
c
c     orthogonalize 2nd block to 1st
c
         do j = 1, nroot
            do i = 1, nroot
               phi(k,j+nroot) = phi(k,j+nroot) - c(j,i)*phi(k,i)
               phik(k,j+nroot) = phik(k,j+nroot) - c(j,i)*phik(k,i)
            end do
         end do
      end do
c
c     orthogonalize 2nd block on itself
c
      sum = 0.0d0
      do k = 1, nvar
         sum = sum + phi(k,nroot+1)*phi(k,nroot+1)
      end do
      sum = sqrt(sum)
c
c     normalize leading vector
c
      do k = 1, nvar
         phi(k,nroot+1) = phi(k,nroot+1) / sum
         phik(k,nroot+1) = phik(k,nroot+1) / sum
      end do
c
c     orthogonalize the rest one-by-one
c
      if (nroot .gt. 1) then
         do i = 2, max(2,nroot)
            do j = 1, i-1
               sum = 0.0d0
               do k = 1, nvar
                  sum = sum + phi(k,i+nroot)*phi(k,j+nroot)
               end do
               do k = 1, nvar
                  phi(k,i+nroot) = phi(k,i+nroot)-sum*phi(k,j+nroot)
                  phik(k,i+nroot) = phik(k,i+nroot)-sum*phik(k,j+nroot)
               end do
            end do
            sum = 0.0d0
            do k = 1, nvar
               sum = sum + phi(k,i+nroot)*phi(k,i+nroot)
            end do
            sum = sqrt(sum)
            do k = 1, nvar
               phi(k,i+nroot) = phi(k,i+nroot) / sum
               phik(k,i+nroot) = phik(k,i+nroot) / sum
            end do
         end do
      end if
c
c     residue of new solution (if restarting, begin here)
c
  140 continue
      do i = 1, nroot
         freq(i) = funit * sign(1.0d0,h(i,i)) * sqrt(abs(h(i,i)))
         freq(i+nroot) = funit * sign(1.0d0,h(i+nroot,i+nroot))
     &                      * sqrt(abs(h(i+nroot,i+nroot)))
         freq(i+npair) = funit * sign(1.0d0,h(i+npair,i+npair))
     &                        * sqrt(abs(h(i+npair,i+npair)))
         do k = 1, nvar
            pk(k) = phik(k,i) - h(i,i)*phi(k,i)
         end do
c
c     use Davidson preconditioner if finding low frequencies
c
         if (factor .gt. 0.0d0) then
            call preconblk (nvar,nblk,iblk,uku,uu,h(i,i),hmin(i),pk)
         end if
c
c     project onto P-space
c
         call qonvec (nvar,np,u,pk,p)
         do k = 1, nvar
            phi(k,i+npair) = p(k)
         end do
      end do
c
c     project out locked roots from components of phi
c
      call project (nvar,nconv,ivb1,nroot,npair)
      call projectk (nvar,nconv,ivb1,nroot,npair)
c
c     reload and orthogonalize to 1st + 2nd
c
      do i = 1, nroot
         do k = 1, nvar
            p(k) = phi(k,i+npair)
         end do
         call gsort (nvar,npair+i-1,p)
         do k = 1, nvar
            phi(k,i+npair) = p(k)
         end do
      end do
c
c     store K on phi
c
      do i = 1, nroot
         do k = 1, nvar
            p(k) = phi(k,i+npair)
         end do
         call konvec (nvar,xm,xe,p,pk)
         do k = 1, nvar
            phik(k,i+npair) = factor * pk(k)
         end do
      end do
c
c     beginning of iterations
c
      iconv = 0
  150 continue
      done = .false.
      iter = iter + 1
c
c     make nbasis-by-nbasis CI matrix
c
      do i = 1, nbasis
         do j = i, nbasis
            h(i,j) = 0.0d0
            do k = 1, nvar
               h(i,j) = h(i,j) + phik(k,i)*phi(k,j)
            end do
            h(j,i) = h(i,j)
         end do
      end do
c
c     list of previous frequencies
c
      do i = 1, npair
         freqold(i) = freq(i)
      end do
c
c     diagonalize and use first nroot solutions as new guess
c
      call transform (nbasis,nbasis,h,c)
c
c     check for collapse based on leading component of ground state
c
      if (iconv.eq.0 .and. nconv.gt.0) then
         sum = sqrt(1.0d0-c(1,1)**2)
         if (sum .gt. 0.9d0) then
            write (iout,160)  nconv-nlock
  160       format (/,' Number of Converged Normal Modes :',6x,i12)
            write (iout,170)
  170       format (/,' VIBBIG  --  Loss of Root Identity; Please',
     &                 ' Try to Restart')
            close (unit=ivb2,status='delete')
            goto 270
         end if
      end if
c
c     list of new frequencies
c
      do i = 1, npair
         freq(i) = funit * sign(1.0d0,h(i,i)) * sqrt(abs(h(i,i)))
      end do
c
c     check if first few have converged
c
      iconv = 0
  180 continue
      dfreq = freqold(iconv+1) - freq(iconv+1)
      if (dfreq*factor.gt.0.0d0 .and. dfreq*factor.lt.wtol) then
         iconv = iconv + 1
         goto 180
      end if
c
c     shift levels of preconditioner matrix; since the Hessian
c     is gradually deflated, reduce effect of the preconditioner
c     based on a simple 1/x curve, the uku levels are squeezed
c     upwards to eventually lead to a unit operator
c
      if (iconv .gt. 0) then
         ratio = dble(nconv+iconv) / dble(nvar)
         shift = uku_min / (1.0d0-ratio)
         shift = shift + h(iconv+nroot,iconv+nroot)
c
c     do a regular shift, which also seems to work
c
         do k = 1, nvar
            uku(k) = uku_max + (uku0(k)-uku_max)*(uku_max-shift)
     &                                 / (uku_max-uku_min)
         end do
c
c     move cursor to end of storage file
c
         do i = 1, nconv
            read (ivb1)  (pk(k),k=1,nvar)
         end do
c
c     norm of residual
c
         do j = 1, iconv
            rnorm = 0.0d0
            do k = 1, nvar
               p(k) = 0.0d0
               pk(k) = 0.0d0
               do i = 1, nbasis
                  p(k) = p(k)+c(i,j)*phi(k,i)
                  pk(k) = pk(k)+c(i,j)*phik(k,i)
               end do
               rnorm = rnorm + (pk(k)-h(j,j)*p(k))**2
            end do
            rnorm = sqrt(rnorm)
c
c     component of root in R-space
c
            do i = 1, 3
               tmp1(i) = 0.0d0
               do k = 1, nvar
                  tmp1(i) = tmp1(i) + ur(k,i)*p(k)
               end do
            end do
            rcomp = 0.0d0
            do k = 1, nvar
               sum = 0.0d0
               do i = 1, 3
                  sum = sum + ur(k,i)*tmp1(i)
               end do
               rcomp = rcomp + sum*sum
            end do
            rcomp = sqrt(rcomp)
c
c     write the converged mode to formatted and binary files
c
            ivib = irange + ifactor*(nconv+j)
            if ((header.or.verbose) .and. j.eq.1) then
               header = .false.
               write (iout,190)
  190          format (/,' Converged Normal Modes from Iterative',
     &                    ' Vibrational Analysis :')
               write (iout,200)
  200          format (/,4x,'Mode',7x,'Frequency',8x,'Delta',10x,
     &                    'R Norm',10x,'Orthog')
               if (.not. verbose) then
                  write (iout,210)
  210             format ()
               end if
            end if
            dfreq = freqold(j) - freq(j)
            write (iout,220)  ivib,freq(j),dfreq,rnorm,rcomp
  220       format (i8,f15.3,3d16.4)
            call prtvib (ivib,p)
            write (ivb1)  (p(k),k=1,nvar)
         end do
         rewind (unit=ivb1)
c
c     update total number of vectors locked on disk
c
         nconv = nconv + iconv
         if (freq(iconv)*factor .ge. fmax*factor) then
            done = .true.
            close (unit=ivb1)
         end if
      end if
c
c     shift frequency arrays by iconv
c
      do i = 1, npair
         freq(i) = freq(i+iconv)
         freqold(i) = freqold(i+iconv)
      end do
      do k = 1, nvar
         do j = 1, nroot+iconv
            tmp1(j) = 0.0d0
            tmp2(j) = 0.0d0
            do i = 1, nbasis
               tmp1(j) = tmp1(j) + c(i,j)*phi(k,i)
               tmp2(j) = tmp2(j) + c(i,j)*phik(k,i)
            end do
         end do
c
c     old solution fills up 2nd block
c
         do j = 1, nroot
            phi(k,j+nroot+iconv) = phi(k,j+iconv)
            phik(k,j+nroot+iconv) = phik(k,j+iconv)
         end do
c
c     new solution fills up 1st block
c
         do j = 1, nroot
            phi(k,j+iconv) = tmp1(j+iconv)
            phik(k,j+iconv) = tmp2(j+iconv)
         end do
c
c     shift index down by iconv
c
         do j = 1, npair
            phi(k,j) = phi(k,j+iconv)
            phik(k,j) = phik(k,j+iconv)
         end do
c
c     orthogonalize 2nd block to 1st + iconv roots
c
         do j = 1, nroot
            do i = 1, nroot
               phi(k,j+nroot) = phi(k,j+nroot)
     &                             - c(j+iconv,i+iconv)*phi(k,i)
               phik(k,j+nroot) = phik(k,j+nroot)
     &                              - c(j+iconv,i+iconv)*phik(k,i)
            end do
            do i = 1, iconv
               phi(k,j+nroot) = phi(k,j+nroot) - c(j+iconv,i)*tmp1(i)
               phik(k,j+nroot) = phik(k,j+nroot) - c(j+iconv,i)*tmp2(i)
            end do
         end do
      end do
c
c     orthogonalize 2nd block on itself
c
      sum = 0.0d0
      do k = 1, nvar
         sum = sum + phi(k,nroot+1)*phi(k,nroot+1)
      end do
      sum = sqrt(sum)
c
c     normalize leading vector
c
      do k = 1, nvar
         phi(k,nroot+1) = phi(k,nroot+1) / sum
         phik(k,nroot+1) = phik(k,nroot+1) / sum
      end do
c
c     orthogonalize the rest one-by-one
c
      if (nroot .gt. 1) then
         do i = 2, max(2,nroot)
            do j = 1, i-1
               sum = 0.0d0
               do k = 1, nvar
                  sum = sum + phi(k,i+nroot)*phi(k,j+nroot)
               end do
               do k = 1, nvar
                  phi(k,i+nroot) = phi(k,i+nroot)-sum*phi(k,j+nroot)
                  phik(k,i+nroot) = phik(k,i+nroot)-sum*phik(k,j+nroot)
               end do
            end do
            sum = 0.0d0
            do k = 1, nvar
               sum = sum + phi(k,i+nroot)*phi(k,i+nroot)
            end do
            sum = sqrt(sum)
            do k = 1, nvar
               phi(k,i+nroot) = phi(k,i+nroot) / sum
               phik(k,i+nroot) = phik(k,i+nroot) / sum
            end do
         end do
      end if
c
c     print a header for the current iteration
c
      if (verbose) then
         write (iout,230)  iter,iconv,nconv
  230    format (/,' Iteration',i7,11x,'New Modes',i6,10x,
     &              ' Total Modes',i6,/)
         write (iout,240)
  240    format (4x,'Mode',7x,'Frequency',8x,'Delta',10x,
     &              'R Norm',10x,'Orthog')
      end if
c
c     norm of residual
c
      do i = 1, nroot
         rnorm = 0.0d0
         do k = 1, nvar
            rnorm = rnorm + (phik(k,i)-h(i+iconv,i+iconv)*phi(k,i))**2
         end do
         rnorm = sqrt(rnorm)
c
c     calculate root's component in R-space
c
         do j = 1, 3
            tmp1(j) = 0.0d0
            do k = 1, nvar
               tmp1(j) = tmp1(j) + ur(k,j)*phi(k,i)
            end do
         end do
         rcomp = 0.0d0
         do k = 1, nvar
            sum = 0.0d0
            do j = 1, 3
               sum = sum + ur(k,j)*tmp1(j)
            end do
            rcomp = rcomp + sum*sum
         end do
         rcomp = sqrt(rcomp)
         dfreq = freqold(i) - freq(i)
         if (verbose) then
            write (iout,250)  irange+ifactor*(i+nconv),
     &                        freq(i),dfreq,rnorm,rcomp
  250       format (i8,f15.3,3d16.4)
         end if
      end do
c
c     save vectors for restart
c
      if (mod(iter,idump) .eq. 0) then
         rewind (unit=ivb2)
         do i = 1, npair
            write (ivb2)  (phi(k,i),k=1,nvar)
            write (ivb2)  (phik(k,i),k=1,nvar)
         end do
      end if
c
c     prepare restart if finished or iterations exhausted
c
      if (done .or. iter.eq.maxiter) then
         write (iout,260)  nconv-nlock
  260    format (/,' Number of Converged Normal Modes :',6x,i12)
         rewind (ivb2)
         do i = 1, npair
            write (ivb2)  (phi(k,i),k=1,nvar)
            write (ivb2)  (phik(k,i),k=1,nvar)
         end do
         close (unit=ivb2)
         goto 270
      end if
c
c     as above, make sure no prior roots are mixed in the basis
c
      do i = 1, npair
         do k = 1, nvar
            p(k) = phi(k,i)
         end do
         call qonvec (nvar,np,u,p,pk)
         do k = 1, nvar
            phi(k,i) = pk(k)
         end  do
         do k = 1, nvar
            p(k) = phik(k,i)
         end do
         call qonvec (nvar,np,u,p,pk)
         do k = 1, nvar
            phik(k,i) = pk(k)
         end do
      end do
c
c     project out locked roots from components of phi
c
      call project (nvar,nconv,ivb1,npair,0)
      call projectk (nvar,nconv,ivb1,npair,0)
c
c     setup next iteration; solution residue, Davidson weight
c
      do i = 1, nroot
         do k = 1, nvar
            pk(k) = phik(k,i) - h(i+iconv,i+iconv)*phi(k,i)
         end do
c
c     use Davidson preconditioner if finding low frequencies
c
         ii = i + iconv
         if (factor .gt. 0.0d0) then
            call preconblk (nvar,nblk,iblk,uku,uu,h(ii,ii),hmin(i),pk)
         end if
c
c     project residual onto P-space
c
         call qonvec (nvar,np,u,pk,p)
         do k = 1, nvar
            phi(k,i+npair) = p(k)
         end do
      end do
c
c     project out locked roots from components of phi
c
      call project (nvar,nconv,ivb1,nroot,npair)
c
c     reload and orthogonalize to 1st + 2nd
c
      do i = 1, nroot
         do k = 1, nvar
            p(k) = phi(k,i+npair)
         end do
         call gsort (nvar,npair+i-1,p)
         do k = 1, nvar
            phi(k,i+npair) = p(k)
         end do
      end do
c
c     store K on phi
c
      do i= 1, nroot
         do k = 1, nvar
            p(k) = phi(k,i+npair)
         end do
         call konvec (nvar,xm,xe,p,pk)
         call qonvec(nvar,np,u,pk,p)
         do k = 1, nvar
            phik(k,i+npair) = factor * p(k)
         end do
      end do
c
c     project out locked roots from components of phik
c
      call projectk (nvar,nconv,ivb1,nroot,npair)
      goto 150
  270 continue
c
c     perform deallocation of some local arrays
c
      deallocate (iblk)
      deallocate (xe)
      deallocate (xm)
      deallocate (p)
      deallocate (pk)
      deallocate (hmin)
      deallocate (uku)
      deallocate (uku0)
      deallocate (uu)
      deallocate (u)
      deallocate (ur)
      deallocate (freq)
      deallocate (freqold)
      deallocate (tmp1)
      deallocate (tmp2)
      deallocate (h)
      deallocate (c)
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine trigger  --  get initial trial eigenvectors  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "trigger" constructs a set of initial trial vectors for
c     use during sliding block iterative matrix diagonalization
c
c
      subroutine trigger (nvar,nbasis,np,ifactor,nblk,iblk,u,uu,p)
      use sizes
      implicit none
      integer i,j,k,m
      integer k0,k1,k2
      integer nvar,nbasis
      integer np,ifactor
      integer nblk,nguess
      integer iblk(*)
      real*8 w,sum
      real*8 random
      real*8 p(*)
      real*8 uu(*)
      real*8 u(nvar,*)
      real*8, allocatable :: tmp(:)
      external random
c
c
c     set the number of random guesses
c
      nguess = 1 + int(dble(nbasis)/dble(nblk))
c
c     zero out the trial vector
c
      do k = 1, nvar
         p(k) = 0.0d0
      end do
c
c     create overlap with the entire P-space
c
      k0 = 0
      k1 = 1
      do i = 1, nblk
         if (i .gt. 1) then
            k0 = k0 + 9*iblk(i-1)**2
            k1 = k1 + 3*iblk(i-1)
         end if
         k2 = k1 + 3*iblk(i) - 1
c
c     scan over rows of the Hessian
c
         m = 0
         do j = 1, 3*iblk(i)
            if (ifactor .eq. 1) then
               if (j .gt. min(nguess,3*iblk(i))) then
                  w = 0.0d0
               else
                  w = random() - 0.5d0
               end if
            else
               if (j .lt. (3*iblk(i)-min(nguess,3*iblk(i))+1)) then
                  w = 0.0d0
               else
                  w = random() - 0.5d0
               end if
            end if
            do k = k1, k2
               m = m + 1
               p(k) = p(k) + w*uu(k0+m)
            end do
         end do
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (tmp(nvar))
c
c     project the vector onto P-space
c
      call qonvec (nvar,np,u,p,tmp)
c
c     perform a normalization
c
      sum = 0.0d0
      do i = 1, nvar
         sum = sum + tmp(i)**2
      end do
      sum = sqrt(sum)
      do i = 1, nvar
         p(i) = tmp(i) / sum
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (tmp)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine trbasis  --  set translation/rotation vectors  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "trbasis" forms translation and rotation basis vectors used
c     during vibrational analysis via block iterative diagonalization
c
c
      subroutine trbasis (nvar,np,xe,u,ur)
      use sizes
      use atomid
      use atoms
      implicit none
      integer i,j,k
      integer nvar,np
      real*8 tmass,sum
      real*8 ra,rha,pr
      real*8 cm(3)
      real*8 p(3)
      real*8 e(3,3)
      real*8 c(3,3)
      real*8 xe(*)
      real*8 u(nvar,*)
      real*8 ur(nvar,*)
c
c
c     zero out the translation and rotation vectors
c
      do i = 1, 6
         do j = 1, nvar
            u(j,i) = 0.0d0
         end do
      end do
c
c     get the total mass of the system
c
      tmass = 0.0d0
      do i = 1, n
         tmass = tmass + mass(i)
      end do
c
c     set basis vectors for translations
c
      do i = 1, n
         u(3*i-2,1) = sqrt(mass(i)/tmass)
         u(3*i-1,1) = 0.0d0
         u(3*i,1) = 0.0d0
         u(3*i-2,2) = 0.0d0
         u(3*i-1,2) = sqrt(mass(i)/tmass)
         u(3*i,2) = 0.0d0
         u(3*i-2,3) = 0.0d0
         u(3*i-1,3) = 0.0d0
         u(3*i,3) = sqrt(mass(i)/tmass)
      end do
c
c     move center of mass to origin
c
      do i = 1, 3
         cm(i) = 0.0d0
      end do
      do i = 1, n
         do j = 1, 3
            cm(j) = cm(j) + xe(3*(i-1)+j)*mass(i)
         end do
      end do
      do i = 1, n
         do j = 1, 3
            xe(3*(i-1)+j) = xe(3*(i-1)+j) - cm(j)/tmass
         end do
      end do
c
c     get the moments of inertia
c
      do i = 1, 3
         e(i,i) = 0.0d0
      end do
      do i = 1, n
         e(1,1) = e(1,1) + ((xe(3*i-1)**2+xe(3*i)**2))*mass(i)
         e(2,2) = e(2,2) + ((xe(3*i-2)**2+xe(3*i)**2))*mass(i)
         e(3,3) = e(3,3) + ((xe(3*i-2)**2+xe(3*i-1)**2))*mass(i)
      end do
      do i = 1, 2
         do j = i+1, 3
            e(i,j) = 0.0d0
            do k = 1, n
               e(i,j) = e(i,j) - xe(3*(k-1)+i)*xe(3*(k-1)+j)*mass(k)
            end do
            e(j,i) = e(i,j)
         end do
      end do
c
c     diagonalize to get principal axes
c
      call jacobi (3,e,cm,c)
c
c     construction of principle rotations
c
      do i = 1, 3
         do j = 1, n
            ra = 0.0d0
            pr = 0.0d0
            do k = 1, 3
               cm(k) = xe(3*(j-1)+k)
               ra = ra + cm(k)**2
               pr = pr + cm(k)*c(k,i)
            end do
            rha = sqrt(ra-pr**2)
            p(1) = c(2,i)*cm(3) - c(3,i)*cm(2)
            p(2) = c(3,i)*cm(1) - c(1,i)*cm(3)
            p(3) = c(1,i)*cm(2) - c(2,i)*cm(1)
            sum = 0.0d0
            do k = 1, 3
               sum = sum + p(k)**2
            end do
            sum = sqrt(sum)
            do k = 1, 3
               ur(3*(j-1)+k,i) = sqrt(mass(j)) * rha*p(k)/sum
            end do
         end do
         sum = 0.0d0
         do j = 1, nvar
            sum = sum + ur(j,i)**2
         end do
         sum = sqrt(sum)
         do j = 1, nvar
            ur(j,i) = ur(j,i) / sum
         end do
      end do
c
c     set basis vectors for rotation
c
      if (np .eq. 6) then
         do i = 1, 3
            do j = 1, nvar
               u(j,i+3) = ur(j,i)
            end do
         end do
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine preconblk  --  precondition atom block Hessian  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "preconblk" applies a preconditioner to an atom block section
c     of the Hessian matrix
c
c
      subroutine preconblk (nvar,nblk,iblk,uku,uu,h,hmin,pk)
      use sizes
      implicit none
      integer i,j,k,l
      integer nvar,nblk
      integer k0,k1,k2,l2
      integer iblk(*)
      real*8 h,hmin
      real*8 uku(*)
      real*8 pk(*)
      real*8 uu(*)
      real*8, allocatable :: d(:)
      real*8, allocatable :: work(:)
c
c
c     find smallest element of |h-uku|
c
      hmin = abs(h-uku(1))
      do k = 2, nvar
         if (abs(h-uku(k)) .lt. hmin) then
            hmin = abs(h-uku(k))
         end if
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (d(nvar))
      allocate (work(nvar))
c
c     assign values to temporary array
c
      do k = 1, nvar
         d(k) = h - uku(k)
      end do
c
c     invert array via d=hmin/d, where hmin=min{|d(k)|}
c
      do k = 1, nvar
         d(k) = hmin / d(k)
      end do
c
c     create overlap with the entire pk array
c
      k0 = 0
      k1 = 1
      do i = 1, nblk
         if (i .gt. 1) then
            k0 = k0 + 9*iblk(i-1)**2
            k1 = k1 + 3*iblk(i-1)
         end if
         k2 = k1 + 3*iblk(i) - 1
c
c    scan over rows of the Hessian, first part
c
         l = 0
         do j = 1, 3*iblk(i)
            l2 = k1 + j - 1
            work(l2) = 0.0d0
            do k = k1, k2
               l = l + 1
               work(l2) = work(l2) + uu(k0+l)*pk(k)
            end do
         end do
c
c    zero out the segment
c
         do k = k1, k2
            pk(k) = 0.0d0
         end do
c
c    scan over rows of the Hessian, second part
c
         l = 0
         do j = 1, 3*iblk(i)
            l2 = k1 + j - 1
            do k = k1, k2
               l = l + 1
               pk(k) = pk(k) + uu(k0+l)*d(l2)*work(l2)
            end do
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (d)
      deallocate (work)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine gsort  --  orthogonal vector via Gram-Schmidt  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "gsort" uses the Gram-Schmidt algorithm to build orthogonal
c     vectors for sliding block interative matrix diagonalization
c
c
      subroutine gsort (nvar,nb,p0)
      use sizes
      use vibs
      implicit none
      integer i,j
      integer nvar,nb
      real*8 sum
      real*8 p0(*)
      real*8, allocatable :: s(:)
      real*8, allocatable :: proj(:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (s(nb))
      allocate (proj(nvar))
c
c     make overlap between two basis sets
c
      do i = 1, nb
         s(i) = 0.0d0
         do j = 1, nvar
            s(i) = s(i) + p0(j)*phi(j,i)
         end do
      end do
c
c     start the Gram-Schmidt procedure
c
      do i = 1, nvar
         proj(i) = 0.0d0
      end do
c
c     construct projector
c
      do i = 1, nb
         do j = 1, nvar
            proj(j) = proj(j) + s(i)*phi(j,i)
         end do
      end do
c
c     apply projector and normalize new vector
c
      sum = 0.0d0
      do i = 1, nvar
         proj(i) = p0(i) - proj(i)
         sum = sum + proj(i)*proj(i)
      end do
      sum = sqrt(sum)
      do i = 1, nvar
         proj(i) = proj(i) / sum
      end do
c
c     return original array updated
c
      do i = 1, nvar
         p0(i) = proj(i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (s)
      deallocate (proj)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine qonvec  --  block iterative vibration utility  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "qonvec" is a vector utility routine used during sliding
c     block iterative matrix diagonalization
c
c
      subroutine qonvec (nvar,np,u,pk,p)
      use sizes
      implicit none
      integer i,j,nvar,np
      real*8 pku(6)
      real*8 pk(*)
      real*8 p(*)
      real*8 u(nvar,*)
c
c
c     operate on vector pk with u-transpose
c
      do i = 1, np
         pku(i) = 0.0d0
         do j = 1, nvar
            pku(i) = pku(i) + u(j,i)*pk(j)
         end do
      end do
c
c     operate with u on the resultant
c
      do i = 1, nvar
         p(i) = 0.0d0
         do j = 1, np
            p(i) = p(i) + u(i,j)*pku(j)
         end do
      end do
c
c     subtract new product from p
c
      do i = 1, nvar
         p(i) = pk(i) - p(i)
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine project  --  remove known vectors from current  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "project" reads locked vectors from a binary file and projects
c     them out of the components of the set of trial eigenvectors
c     using the relation Y = X - U * U^T * X
c
c
      subroutine project (nvar,nconv,ivb1,ns,m)
      use sizes
      use vibs
      implicit none
      integer i,j,k
      integer nvar,nconv
      integer ivb1,ns,m
      real*8, allocatable :: temp(:)
      real*8, allocatable :: u(:)
c
c
c     zero the temporary storage array
c
      do k = 1, nvar
         do i = 1, ns
            pwork(k,i+m) = 0.0d0
         end do
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (temp(ns))
      allocate (u(nvar))
c
c     read and scan over the locked eigenvectors
c
      do i = 1, nconv
         read (ivb1)  (u(k),k=1,nvar)
         do j = 1, ns
            temp(j) = 0.0d0
            do k = 1, nvar
               temp(j) = temp(j) + u(k)*phi(k,j+m)
            end do
         end do
         do j = 1, ns
            do k = 1, nvar
               pwork(k,j+m) = pwork(k,j+m) + u(k)*temp(j)
            end do
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (temp)
      deallocate (u)
c
c     project locked vectors out of the current set
c
      do k = 1, nvar
         do i = 1, ns
            phi(k,i+m) = phi(k,i+m) - pwork(k,i+m)
         end do
      end do
      if (nconv .gt. 0)  rewind (unit=ivb1)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine projectk  --  remove known vectors from current  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "projectk" reads locked vectors from a binary file and projects
c     them out of the components of the set of trial eigenvectors
c     using the relation Y = X - U * U^T * X
c
c
      subroutine projectk (nvar,nconv,ivb1,ns,m)
      use sizes
      use vibs
      implicit none
      integer i,j,k
      integer nvar,nconv
      integer ivb1,ns,m
      real*8, allocatable :: temp(:)
      real*8, allocatable :: u(:)
c
c
c     zero the temporary storage array
c
      do k = 1, nvar
         do i = 1, ns
            pwork(k,i+m) = 0.0d0
         end do
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (temp(ns))
      allocate (u(nvar))
c
c     read and scan over the locked eigenvectors
c
      do i = 1, nconv
         read (ivb1)  (u(k),k=1,nvar)
         do j = 1, ns
            temp(j) = 0.0d0
            do k = 1, nvar
               temp(j) = temp(j) + u(k)*phik(k,j+m)
            end do
         end do
         do j = 1, ns
            do k = 1, nvar
               pwork(k,j+m) = pwork(k,j+m) + u(k)*temp(j)
            end do
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (temp)
      deallocate (u)
c
c     project locked vectors out of the current set
c
      do k = 1, nvar
         do i = 1, ns
            phik(k,i+m) = phik(k,i+m) - pwork(k,i+m)
         end do
      end do
      if (nconv .gt. 0)  rewind (unit=ivb1)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine konvec  --  evaluate Hessian-vector product  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "konvec" finds a Hessian-vector product via finite-difference
c     evaluation of the gradient based on atomic displacements
c
c
      subroutine konvec (nvar,xm,qe,uvec,kuvec)
      use sizes
      use atomid
      use atoms
      use units
      implicit none
      integer i,j,k,nvar
      real*8 e,term
      real*8 sum,eps
      real*8 xm(*)
      real*8 qe(*)
      real*8 uvec(*)
      real*8 kuvec(*)
      real*8, allocatable :: delta(:)
      real*8, allocatable :: grd1(:,:)
      real*8, allocatable :: grd2(:,:)
c
c
c     estimate displacement based on total average
c
      sum = 0.0d0
      do i = 1, nvar
         sum = sum + uvec(i)*uvec(i)/xm(i)
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (delta(nvar))
      allocate (grd1(3,n))
      allocate (grd2(3,n))
c
c     store the coordinate displacements
c
      eps = 0.001d0 / sqrt(sum)
      do i = 1, nvar
         delta(i) = eps * uvec(i) / sqrt(xm(i))
      end do
c
c     compute the forward displacement
c
      do i = 1, n
         k = 3 * (i-1)
         x(i) = bohr * (qe(k+1)+delta(k+1))
         y(i) = bohr * (qe(k+2)+delta(k+2))
         z(i) = bohr * (qe(k+3)+delta(k+3))
      end do
      call gradient (e,grd1)
c
c     compute the backward displacement
c
      do i = 1, n
         k = 3 * (i-1)
         x(i) = bohr * (qe(k+1)-delta(k+1))
         y(i) = bohr * (qe(k+2)-delta(k+2))
         z(i) = bohr * (qe(k+3)-delta(k+3))
      end do
      call gradient (e,grd2)
c
c     update via finite differences
c
      term = 0.5d0 * bohr / (eps * hartree)
      do i = 1, n
         k = 3 * (i-1)
         do j = 1, 3
            kuvec(k+j) = term * (grd1(j,i)-grd2(j,i)) / sqrt(xm(k+j))
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (delta)
      deallocate (grd1)
      deallocate (grd2)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine transform  --  diagonalize trial basis vectors  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "transform" diagonalizes the current basis vectors to produce
c     trial roots for sliding block iterative matrix diagonalization
c
c
      subroutine transform (ns,nb,h,c)
      implicit none
      integer i,j,k,ns,nb
      real*8 h(nb,*)
      real*8 c(nb,*)
      real*8, allocatable :: e1(:)
      real*8, allocatable :: h1(:)
      real*8, allocatable :: c1(:,:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (e1(ns))
      allocate (h1((ns+1)*ns/2))
      allocate (c1(ns,ns))
c
c     pack the upper triangle of matrix
c
      k = 0
      do i = 1, ns
         do j = i, ns
            k = k + 1
            h1(k) = h(i,j)
         end do
      end do
c
c     perform the matrix diagonalization
c
      call diagq (ns,ns,h1,e1,c1)
c
c     copy values into the return arrays
c
      do i = 1, ns
         do j = 1, ns
            h(i,j) = 0.0d0
            c(i,j) = c1(i,j)
         end do
         h(i,i) = e1(i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (e1)
      deallocate (h1)
      deallocate (c1)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine diagblk  -- diagonalization for atom block  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "diagblk" performs diagonalization of the Hessian for a
c     block of atoms within a larger system
c
c
      subroutine diagblk (k0,k1,n,vector,wres)
      use sizes
      implicit none
      integer i,j,k,m
      integer n,k0,k1
      real*8 wres(*)
      real*8 vector(*)
      real*8, allocatable :: hval(:)
      real*8, allocatable :: hres(:)
      real*8, allocatable :: hvec(:,:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (hval(n))
      allocate (hres((n+1)*n/2))
      allocate (hvec(n,n))
c
c     pack the upper triangle of matrix
c
      k = 0
      do i = 1, n
         m = k0 + (i-1)*n
         do j = i, n
            k = k + 1
            hres(k) = vector(m+j)
         end do
      end do
c
c     perform the matrix diagonalization
c
      call diagq (n,n,hres,hval,hvec)
c
c     copy values into return arrays
c
      k = 0
      do i = 1, n
         do j = 1, n
            k = k + 1
            vector(k0+k) = hvec(j,i)
         end do
      end do
      do i = 1, n
         wres(k1+i-1) = hval(i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (hval)
      deallocate (hres)
      deallocate (hvec)
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine prtvib  --  output of vibrational mode  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "prtvib" writes to an external disk file a series of
c     coordinate sets representing motion along a vibrational
c     normal mode
c
c
      subroutine prtvib (ivib,p)
      use sizes
      use atoms
      use files
      implicit none
      integer i,j,k
      integer ivib,ixyz
      integer lext,nview
      integer freeunit
      real*8 ratio
      real*8 p(*)
      real*8, allocatable :: xref(:)
      real*8, allocatable :: yref(:)
      real*8, allocatable :: zref(:)
      character*7 ext
      character*240 xyzfile
c
c
c     create a name for the vibrational displacement file
c
      lext = 3
      call numeral (ivib,ext,lext)
      xyzfile = filename(1:leng)//'.'//ext(1:lext)
      ixyz = freeunit ()
      call version (xyzfile,'new')
      open (unit=ixyz,file=xyzfile,status='new')
c
c     perform dynamic allocation of some local arrays
c
      allocate (xref(n))
      allocate (yref(n))
      allocate (zref(n))
c
c     store the original atomic coordinates
c
      do i = 1, n
         xref(i) = x(i)
         yref(i) = y(i)
         zref(i) = z(i)
      end do
c
c     make file with plus and minus the current vibration
c
      nview = 3
      do i = -nview, nview
         ratio = dble(i) / dble(nview)
         do k = 1, n
            j = 3 * (k-1)
            x(k) = xref(k) + ratio*p(j+1)
            y(k) = yref(k) + ratio*p(j+2)
            z(k) = zref(k) + ratio*p(j+3)
         end do
         call prtxyz (ixyz)
      end do
      close (unit=ixyz)
c
c     restore the original atomic coordinates
c
      do i = 1, n
         x(i) = xref(i)
         y(i) = yref(i)
         z(i) = zref(i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xref)
      deallocate (yref)
      deallocate (zref)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine hessblk  --  Hessian elements for atom block  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "hessblk" calls subroutines to calculate the Hessian elements
c     for each atom in turn with respect to Cartesian coordinates
c
c
      subroutine hessblk (amass,k0,i1,i2,vector)
      use sizes
      use atoms
      use bound
      use couple
      use hescut
      use hessn
      use inform
      use iounit
      use limits
      use mpole
      use potent
      use rigid
      use usage
      use vdw
      use vdwpot
      use units
      implicit none
      integer i,j,k
      integer ii,k0
      integer i1,i2
      real*8 ami,amik
      real*8 cutoff,rdn
      real*8 amass(*)
      real*8 vector(*)
      real*8, allocatable :: xred(:)
      real*8, allocatable :: yred(:)
      real*8, allocatable :: zred(:)
      logical first
      save first
      data first  / .true. /
c
c
c     maintain any periodic boundary conditions
c
      if (use_bounds .and. .not.use_rigid)  call bounds
c
c     update the pairwise interaction neighbor lists
c
      if (use_list)  call nblist
c
c     many implicit solvation models require Born radii
c
      if (use_born)  call born
c
c     alter bond and torsion constants for pisystem
c
      if (use_orbit)  call picalc
c
c     compute the induced dipoles at polarizable atoms
c
      if (use_polar) then
         call chkpole
         call rotpole
         call induce
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (xred(n))
      allocate (yred(n))
      allocate (zred(n))
c
c     calculate the "reduced" atomic coordinates
c
      if (use_vdw) then
         do i = 1, n
            ii = ired(i)
            rdn = kred(i)
            xred(i) = rdn*(x(i)-x(ii)) + x(ii)
            yred(i) = rdn*(y(i)-y(ii)) + y(ii)
            zred(i) = rdn*(z(i)-z(ii)) + z(ii)
         end do
      end if
c
c     perform dynamic allocation of some global arrays
c
      if (first) then
         first = .false.
         if (.not. allocated(hessx))  allocate (hessx(3,n))
         if (.not. allocated(hessy))  allocate (hessy(3,n))
         if (.not. allocated(hessz))  allocate (hessz(3,n))
      end if
c
c     zero out the Hessian elements for the current atom
c
      ii = 0
      do i = i1, i2
         if (use(i)) then
            do k = i1, i2
               do j = 1, 3
                  hessx(j,k) = 0.0d0
                  hessy(j,k) = 0.0d0
                  hessz(j,k) = 0.0d0
               end do
            end do
c
c     remove any previous use of the replicates method
c
            cutoff = 0.0d0
            call replica (cutoff)
c
c     call the local geometry Hessian component routines
c
            if (use_bond)  call ebond2 (i)
            if (use_angle)  call eangle2 (i)
            if (use_strbnd)  call estrbnd2 (i)
            if (use_urey)  call eurey2 (i)
            if (use_angang)  call eangang2 (i)
            if (use_opbend)  call eopbend2 (i)
            if (use_opdist)  call eopdist2 (i)
            if (use_improp)  call eimprop2 (i)
            if (use_imptor)  call eimptor2 (i)
            if (use_tors)  call etors2 (i)
            if (use_pitors)  call epitors2 (i)
            if (use_strtor)  call estrtor2 (i)
            if (use_tortor)  call etortor2 (i)
c
c     call the van der Waals Hessian component routines
c
            if (use_vdw) then
               if (vdwtyp .eq. 'LENNARD-JONES') then
                  call elj2 (i,xred,yred,zred)
               else if (vdwtyp .eq. 'BUCKINGHAM') then
                  call ebuck2 (i,xred,yred,zred)
               else if (vdwtyp .eq. 'MM3-HBOND') then
                  call emm3hb2 (i,xred,yred,zred)
               else if (vdwtyp .eq. 'BUFFERED-14-7') then
                  call ehal2 (i,xred,yred,zred)
               else if (vdwtyp .eq. 'GAUSSIAN') then
                  call egauss2 (i,xred,yred,zred)
               end if
            end if
c
c     call the electrostatic Hessian component routines
c
            if (use_charge) call echarge2 (i)
            if (use_chgdpl)  call echgdpl2 (i)
            if (use_dipole)  call edipole2 (i)
            if (use_mpole)   call empole2 (i)
            if (use_polar)  call epolar2 (i)
            if (use_rxnfld)   call erxnfld2 (i)
c
c     call any miscellaneous Hessian component routines
c
            if (use_solv)  call esolv2 (i)
            if (use_metal)  call emetal2 (i)
            if (use_geom)  call egeom2 (i)
            if (use_extra)  call extra2 (i)
c
c     store Hessian for the current atom block as a vector
c
            ami = bohr**2 / (hartree*sqrt(amass(i)))
            do k = i1, i2
               amik = ami / sqrt(amass(k))
               do j = 1, 3
                  ii = ii + 1
                  vector(k0+ii) = hessx(j,k) * amik
               end do
            end do
            do k = i1, i2
               amik = ami / sqrt(amass(k))
               do j = 1, 3
                  ii = ii + 1
                  vector(k0+ii) = hessy(j,k) * amik
               end do
            end do
            do k = i1, i2
               amik = ami / sqrt(amass(k))
               do j = 1, 3
                  ii = ii + 1
                  vector(k0+ii) = hessz(j,k) * amik
               end do
            end do
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xred)
      deallocate (yred)
      deallocate (zred)
      return
      end
