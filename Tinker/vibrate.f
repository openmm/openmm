c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  program vibrate  --  vibrational analysis and normal modes  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "vibrate" performs a vibrational normal mode analysis; the
c     Hessian matrix of second derivatives is determined and then
c     diagonalized both directly and after mass weighting; output
c     consists of the eigenvalues of the force constant matrix as
c     well as the vibrational frequencies and displacements
c
c
      program vibrate
      use sizes
      use atomid
      use atoms
      use files
      use hescut
      use iounit
      use math
      use units
      use usage
      implicit none
      integer i,j,k,m
      integer ixyz,ihess
      integer lext,freeunit
      integer nfreq,ndummy
      integer nvib,ivib
      integer nview,next
      integer nlist,ilist
      integer, allocatable :: list(:)
      integer, allocatable :: iv(:)
      integer, allocatable :: hindex(:)
      integer, allocatable :: hinit(:,:)
      integer, allocatable :: hstop(:,:)
      real*8 factor,vnorm
      real*8 sum,scale,ratio
      real*8, allocatable :: xref(:)
      real*8, allocatable :: yref(:)
      real*8, allocatable :: zref(:)
      real*8, allocatable :: mass2(:)
      real*8, allocatable :: h(:)
      real*8, allocatable :: eigen(:)
      real*8, allocatable :: matrix(:)
      real*8, allocatable :: hdiag(:,:)
      real*8, allocatable :: vects(:,:)
      logical exist,query
      character*1 letter
      character*7 ext
      character*240 xyzfile
      character*240 record
      character*240 string
c
c
c     set up the structure and mechanics calculation
c
      call initial
      call getxyz
      call mechanic
c
c     perform dynamic allocation of some local arrays
c
      nfreq = 3 * nuse
      allocate (mass2(n))
      allocate (hinit(3,n))
      allocate (hstop(3,n))
      allocate (hdiag(3,n))
      allocate (hindex(nfreq*(nfreq-1)/2))
      allocate (h(nfreq*(nfreq-1)/2))
      allocate (matrix(nfreq*(nfreq+1)/2))
c
c     initialize various things needed for vibrations
c
      ndummy = 0
      do i = 1, n
         if (use(i) .and. atomic(i).eq.0) then
            ndummy = ndummy + 1
            mass(i) = 0.001d0
         end if
         mass2(i) = sqrt(mass(i))
      end do
      nvib = nfreq - 3*ndummy
c
c     calculate the Hessian matrix of second derivatives
c
      hesscut = 0.0d0
      call hessian (h,hinit,hstop,hindex,hdiag)
c
c     store upper triangle of the Hessian in "matrix"
c
      ihess = 0
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               ihess = ihess + 1
               matrix(ihess) = hdiag(j,i)
               do k = hinit(j,i), hstop(j,i)
                  m = (hindex(k)+2) / 3
                  if (use(m)) then
                     ihess = ihess + 1
                     matrix(ihess) = h(k)
                  end if
               end do
            end do
         end if
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (eigen(nfreq))
      allocate (vects(nfreq,nfreq))
c
c     perform diagonalization to get Hessian eigenvalues
c
      call diagq (nfreq,nfreq,matrix,eigen,vects)
      write (iout,10)
   10 format (/,' Eigenvalues of the Hessian Matrix :',/)
      write (iout,20)  (i,eigen(i),i=1,nvib)
   20 format (5(i5,f10.3))
c
c     store upper triangle of the mass-weighted Hessian matrix
c
      ihess = 0
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               ihess = ihess + 1
               matrix(ihess) = hdiag(j,i) / mass(i)
               do k = hinit(j,i), hstop(j,i)
                  m = (hindex(k)+2) / 3
                  if (use(m)) then
                     ihess = ihess + 1
                     matrix(ihess) = h(k) / (mass2(i)*mass2(m))
                  end if
               end do
            end do
         end if
      end do
c
c     diagonalize to get vibrational frequencies and normal modes
c
      call diagq (nfreq,nfreq,matrix,eigen,vects)
      factor = sqrt(convert) / (2.0d0*pi*lightspd)
      do i = 1, nvib
         eigen(i) = factor * sign(1.0d0,eigen(i)) * sqrt(abs(eigen(i)))
      end do
      write (iout,30)
   30 format (/,' Vibrational Frequencies (cm-1) :',/)
      write (iout,40)  (i,eigen(i),i=1,nvib)
   40 format (5(i5,f10.3))
c
c     perform deallocation of some local arrays
c
      deallocate (hinit)
      deallocate (hstop)
      deallocate (hdiag)
      deallocate (h)
      deallocate (matrix)
c
c     form Cartesian coordinate displacements from normal modes
c
      do i = 1, nvib
         vnorm = 0.0d0
         do j = 1, nfreq
            k = iuse((j+2)/3)
            vects(j,i) = vects(j,i) / mass2(k)
            vnorm = vnorm + vects(j,i)**2
         end do
         vnorm = sqrt(vnorm)
         do j = 1, nfreq
            vects(j,i) = vects(j,i) / vnorm
         end do
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (list(nfreq))
      allocate (iv(nfreq))
      allocate (xref(n))
      allocate (yref(n))
      allocate (zref(n))
c
c     try to get output vibrational modes from command line
c
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         query = .false.
         letter = string(1:1)
         call upcase (letter)
         if (letter .eq. 'A') then
            nlist = nvib
            do i = 1, nlist
               list(i) = i
            end do
         else
            nlist = 0
            do i = 1, nvib
               read (string,*,err=50,end=50)  k
               if (k.ge.1 .and. k.le.nvib) then
                  nlist = nlist + 1
                  list(nlist) = k
               else
                  k = abs(k)
                  call nextarg (string,exist)
                  read (string,*,err=50,end=50)  m
                  m = min(abs(m),nvib)
                  do j = k, m
                     nlist = nlist + 1
                     list(nlist) = j
                  end do
               end if
               call nextarg (string,exist)
            end do
   50       continue
         end if
      end if
c
c     ask the user for the vibrational modes to be output
c
      if (query) then
         write (iout,60)
   60    format (/,' Enter Vibrations to Output [List, A=All',
     &              ' or <CR>=Exit] :  ',$)
         read (input,70)  record
   70    format (a240)
         letter = ' '
         next = 1
         call gettext (record,letter,next)
         call upcase (letter)
         if (letter .eq. ' ') then
            nlist = 0
         else if (letter .eq. 'A') then
            nlist = nvib
            do i = 1, nlist
               list(i) = i
            end do
         else
            do i = 1, nvib
               iv(i) = 0
            end do
            read (record,*,err=80,end=80)  (iv(i),i=1,nvib)
   80       continue
            nlist = 0
            i = 1
            do while (iv(i) .ne. 0)
               k = iv(i)
               if (k.ge.1 .and. k.le.nvib) then
                  nlist = nlist + 1
                  list(nlist) = k
               else
                  k = abs(k)
                  m = min(abs(iv(i+1)),nvib)
                  do j = k, m
                     nlist = nlist + 1
                     list(nlist) = j
                  end do
                  i = i + 1
               end if
               i = i + 1
            end do
         end if
      end if
c
c     print the vibrational frequencies and normal modes
c
      do ilist = 1, nlist
         ivib = list(ilist)
         write (iout,90)  ivib,eigen(ivib)
   90    format (/,' Vibrational Normal Mode',i6,' with Frequency',
     &              f11.3,' cm-1',
     &           //,5x,'Atom',5x,'Delta X',5x,'Delta Y',5x,'Delta Z',/)
         do i = 1, nuse
            j = 3 * (i-1)
            write (iout,100)  iuse(i),vects(j+1,ivib),vects(j+2,ivib),
     &                        vects(j+3,ivib)
  100       format (4x,i5,3f12.6)
         end do
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
c     store the original atomic coordinates
c
         do i = 1, n
            xref(i) = x(i)
            yref(i) = y(i)
            zref(i) = z(i)
         end do
c
c     scale based on the maximum displacement along the mode
c
         scale = 0.0d0
         do i = 1, nuse
            j = 3 * (i-1)
            sum = 0.0d0
            do k = 1, 3
               sum = sum + vects(j+k,ivib)**2
            end do
            scale = max(sum,scale)
         end do
         scale = 0.1d0 * n**(1.0d0/3.0d0) / sqrt(scale)
c
c     make file with plus and minus the current vibration
c
         nview = 3
         do i = -nview, nview
            ratio = scale * dble(i) / dble(nview)
            do k = 1, nuse
               j = 3 * (k-1)
               m = iuse(k)
               x(m) = xref(m) + ratio*vects(j+1,ivib)
               y(m) = yref(m) + ratio*vects(j+2,ivib)
               z(m) = zref(m) + ratio*vects(j+3,ivib)
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
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (list)
      deallocate (iv)
      deallocate (mass2)
      deallocate (xref)
      deallocate (yref)
      deallocate (zref)
      deallocate (eigen)
      deallocate (vects)
c
c     perform any final tasks before program exit
c
      call final
      end
