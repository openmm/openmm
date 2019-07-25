c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  program testhess  --  Hessian matrix test; cart. version  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "testhess" computes and compares the analytical and numerical
c     Hessian matrices of the potential energy function with respect
c     to Cartesian coordinates
c
c
      program testhess
      use sizes
      use atoms
      use files
      use hescut
      use inform
      use iounit
      use usage
      implicit none
      integer i,j,k,m
      integer ii,jj
      integer ixyz,ihes
      integer index,maxnum
      integer next,frame
      integer freeunit
      integer, allocatable :: hindex(:)
      integer, allocatable :: hinit(:,:)
      integer, allocatable :: hstop(:,:)
      real*8 energy,e,old,eps,eps0
      real*8 diff,delta,sum
      real*8, allocatable :: h(:)
      real*8, allocatable :: g(:,:)
      real*8, allocatable :: g0(:,:)
      real*8, allocatable :: hdiag(:,:)
      real*8, allocatable :: nhess(:,:,:,:)
      logical doanalyt,donumer
      logical dograd,dofull
      logical exist,query
      logical identical
      character*1 answer
      character*1 axis(3)
      character*240 xyzfile
      character*240 hessfile
      character*240 record
      character*240 string
      external energy
      data axis  / 'X','Y','Z' /
c
c
c     set up the structure and mechanics calculation
c
      call initial
      call getxyz
      call mechanic
c
c     set difference threshhold via the energy precision
c
      delta = 1.0d-4
      if (digits .ge. 6)  delta = 1.0d-6
      if (digits .ge. 8)  delta = 1.0d-8
c
c     decide whether to do an analytical Hessian calculation
c
      doanalyt = .true.
      call nextarg (answer,exist)
      if (.not. exist) then
         write (iout,10)
   10    format (/,' Compute Analytical Hessian Matrix [Y] :  ',$)
         read (input,20)  record
   20    format (a240)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .eq. 'N')  doanalyt = .false.
c
c     decide whether to do a numerical Hessian calculation
c
      donumer = .false.
      maxnum = 300
      if (n .le. maxnum) then
         donumer = .true.
         call nextarg (answer,exist)
         if (.not. exist) then
            write (iout,30)
   30       format (/,' Compute Numerical Hessian Matrix [Y] :   ',$)
            read (input,40)  record
   40       format (a240)
            next = 1
            call gettext (record,answer,next)
         end if
         call upcase (answer)
         if (answer .eq. 'N')  donumer = .false.
      end if
c
c     get numerical Hessian from either gradient or energy
c
      if (donumer) then
         dograd = .true.
         call nextarg (answer,exist)
         if (.not. exist) then
            write (iout,50)
   50       format (/,' Numerical Hessian from Gradient',
     &                 ' or Function [G] :  ',$)
            read (input,60)  record
   60       format (a240)
            next = 1
            call gettext (record,answer,next)
         end if
         call upcase (answer)
         if (answer .eq. 'F')  dograd = .false.
c
c     get the stepsize for numerical Hessian calculation
c
         eps = -1.0d0
         eps0 = 1.0d-3
         if (dograd)  eps0 = 1.0d-5
         query = .true.
         call nextarg (string,exist)
         if (exist) then
            read (string,*,err=70,end=70)  eps
            query = .false.
         end if
   70    continue
         if (query) then
            write (iout,80)  eps0
   80       format (/,' Enter Finite Difference Stepsize [',d8.1,
     &                 ' Ang] :  ',$)
            read (input,90,err=70)  eps
   90       format (f20.0)
         end if
         if (eps .le. 0.0d0)  eps = eps0
      end if
c
c     decide whether to output results by Hessian component
c
      dofull = .false.
      if (n.le.20 .and. donumer) then
         call nextarg (answer,exist)
         if (.not. exist) then
            write (iout,100)
  100       format (/,' List Individual Hessian Components [N] :   ',$)
            read (input,110)  record
  110       format (a240)
            next = 1
            call gettext (record,answer,next)
         end if
         call upcase (answer)
         if (answer .eq. 'Y')  dofull = .true.
      end if
c
c     reopen the coordinates file and read the first structure
c
      frame = 0
      ixyz = freeunit ()
      xyzfile = filename
      call suffix (xyzfile,'xyz','old')
      open (unit=ixyz,file=xyzfile,status ='old')
      rewind (unit=ixyz)
      call readxyz (ixyz)
c
c     perform dynamic allocation of some local arrays
c
      allocate (hindex(3*n*(3*n-1)/2))
      allocate (hinit(3,n))
      allocate (hstop(3,n))
      allocate (h(3*n*(3*n-1)/2))
      allocate (g(3,n))
      allocate (g0(3,n))
      allocate (hdiag(3,n))
      if (n .le. maxnum)  allocate (nhess(3,n,3,n))
c
c     perform analysis for each successive coordinate structure
c
      do while (.not. abort)
         frame = frame + 1
         if (frame .gt. 1) then
            write (iout,120)  frame
  120       format (/,' Analysis for Archive Structure :',8x,i8)
         end if
c
c     get the analytical Hessian matrix elements
c
         identical = .true.
         if (doanalyt) then
            if (verbose) then
               write (iout,130)
  130          format ()
            end if
            hesscut = 0.0d0
            call hessian (h,hinit,hstop,hindex,hdiag)
         end if
c
c     get the two-sided numerical Hessian matrix elements
c
         do i = 1, n
            if (donumer .and. use(i)) then
               old = x(i)
               x(i) = x(i) - 0.5d0*eps
               if (dograd) then
                  call gradient (e,g)
               else
                  call numgrad (energy,g,eps)
               end if
               do k = 1, n
                  do j = 1, 3
                     g0(j,k) = g(j,k)
                  end do
               end do
               x(i) = x(i) + eps
               if (dograd) then
                  call gradient (e,g)
               else
                  call numgrad (energy,g,eps)
               end if
               x(i) = old
               do k = 1, n
                  do j = 1, 3
                     nhess(j,k,1,i) = (g(j,k) - g0(j,k)) / eps
                  end do
               end do
               old = y(i)
               y(i) = y(i) - 0.5d0*eps
               if (dograd) then
                  call gradient (e,g)
               else
                  call numgrad (energy,g,eps)
               end if
               do k = 1, n
                  do j = 1, 3
                     g0(j,k) = g(j,k)
                  end do
               end do
               y(i) = y(i) + eps
               if (dograd) then
                  call gradient (e,g)
               else
                  call numgrad (energy,g,eps)
               end if
               y(i) = old
               do k = 1, n
                  do j = 1, 3
                     nhess(j,k,2,i) = (g(j,k) - g0(j,k)) / eps
                  end do
               end do
               old = z(i)
               z(i) = z(i) - 0.5d0*eps
               if (dograd) then
                  call gradient (e,g)
               else
                  call numgrad (energy,g,eps)
               end if
               do k = 1, n
                  do j = 1, 3
                     g0(j,k) = g(j,k)
                  end do
               end do
               z(i) = z(i) + eps
               if (dograd) then
                  call gradient (e,g)
               else
                  call numgrad (energy,g,eps)
               end if
               z(i) = old
               do k = 1, n
                  do j = 1, 3
                     nhess(j,k,3,i) = (g(j,k) - g0(j,k)) / eps
                  end do
               end do
            end if
c
c     compare the analytical and numerical diagonal elements
c
            if (doanalyt .and. donumer) then
               do j = 1, 3
                  diff = abs(hdiag(j,i)-nhess(j,i,j,i))
                  if (diff .gt. delta) then
                     if (identical) then
                        identical = .false.
                        write (iout,140)
  140                   format (/,' Comparison of Analytical and',
     &                             ' Numerical Hessian Elements :',
     &                          //,3x,'1st Atom',4x,'2nd Atom',
     &                             9x,'Analytical',8x,'Numerical',
     &                             7x,'Difference',/)
                     end if
                     if (digits .ge. 8) then
                        write (iout,150)  i,axis(j),i,axis(j),
     &                                    hdiag(j,i),nhess(j,i,j,i),
     &                                    hdiag(j,i)-nhess(j,i,j,i)
  150                   format (1x,i6,' (',a1,') ',1x,i6,' (',
     &                             a1,') ',1x,3f17.8)
                     else if (digits .ge. 6) then
                        write (iout,160)  i,axis(j),i,axis(j),
     &                                    hdiag(j,i),nhess(j,i,j,i),
     &                                    hdiag(j,i)-nhess(j,i,j,i)
  160                   format (1x,i6,' (',a1,') ',1x,i6,' (',
     &                             a1,') ',1x,3f17.6)
                     else
                        write (iout,170)  i,axis(j),i,axis(j),
     &                                    hdiag(j,i),nhess(j,i,j,i),
     &                                    hdiag(j,i)-nhess(j,i,j,i)
  170                   format (1x,i6,' (',a1,') ',1x,i6,' (',
     &                             a1,') ',1x,3f17.4)
                     end if
                  end if
c
c     compare the analytical and numerical off-diagonal elements
c
                  do k = hinit(j,i), hstop(j,i)
                     index = hindex(k)
                     jj = mod(index,3)
                     if (jj .eq. 0)  jj = 3
                     ii = (index+2) / 3
                     diff = abs(h(k)-nhess(jj,ii,j,i))
                     if (diff .gt. delta) then
                        if (identical) then
                           identical = .false.
                           write (iout,180)
  180                      format (/,' Comparison of Analytical and',
     &                                ' Numerical Hessian Elements :',
     &                             //,3x,'1st Atom',4x,'2nd Atom',
     &                                9x,'Analytical',8x,'Numerical',
     &                                7x,'Difference',/)
                        end if
                        if (digits .ge. 8) then
                           write (iout,190)  i,axis(j),ii,axis(jj),
     &                                       h(k),nhess(jj,ii,j,i),
     &                                       h(k)-nhess(jj,ii,j,i)
  190                      format (1x,i6,' (',a1,') ',1x,i6,' (',
     &                                a1,') ',1x,3f17.8)
                        else if (digits .ge. 6) then
                           write (iout,200)  i,axis(j),ii,axis(jj),
     &                                       h(k),nhess(jj,ii,j,i),
     &                                       h(k)-nhess(jj,ii,j,i)
  200                      format (1x,i6,' (',a1,') ',1x,i6,' (',
     &                                a1,') ',1x,3f17.6)
                        else
                           write (iout,210)  i,axis(j),ii,axis(jj),
     &                                       h(k),nhess(jj,ii,j,i),
     &                                       h(k)-nhess(jj,ii,j,i)
  210                      format (1x,i6,' (',a1,') ',1x,i6,' (',
     &                                a1,') ',1x,3f17.4)
                        end if
                     end if
                  end do
               end do
            end if
         end do
c
c     success if the analytical and numerical elements are the same
c
         if (doanalyt .and. donumer) then
            if (identical) then
               write (iout,220)
  220          format (/,' Analytical and Numerical Hessian Elements',
     &                    ' are Identical')
            end if
         end if
c
c     write out the diagonal Hessian elements for each atom
c
         if (doanalyt) then
            if (digits .ge. 8) then
               write (iout,230)
  230          format (/,' Diagonal Hessian Elements for Each Atom :',
     &                    //,6x,'Atom',21x,'X',19x,'Y',19x,'Z',/)
            else if (digits .ge. 6) then
               write (iout,240)
  240          format (/,' Diagonal Hessian Elements for Each Atom :',
     &                    //,6x,'Atom',19x,'X',17x,'Y',17x,'Z',/)
            else
               write (iout,250)
  250          format (/,' Diagonal Hessian Elements for Each Atom :',
     &                    //,6x,'Atom',17x,'X',15x,'Y',15x,'Z',/)
            end if
            do i = 1, n
               if (digits .ge. 8) then
                  write (iout,260)  i,(hdiag(j,i),j=1,3)
  260             format (i10,5x,3f20.8)
               else if (digits .ge. 6) then
                  write (iout,270)  i,(hdiag(j,i),j=1,3)
  270             format (i10,5x,3f18.6)
               else
                  write (iout,280)  i,(hdiag(j,i),j=1,3)
  280             format (i10,5x,3f16.4)
               end if
            end do
         end if
c
c     write out the Hessian trace as sum of diagonal elements
c
         if (doanalyt) then
            sum = 0.0d0
            do i = 1, n
               do j = 1, 3
                  sum = sum + hdiag(j,i)
               end do
            end do
            if (digits .ge. 8) then
               write (iout,290)  sum
  290          format (/,' Sum of Diagonal Hessian Elements :',6x,f20.8)
            else if (digits .ge. 6) then
               write (iout,300)  sum
  300          format (/,' Sum of Diagonal Hessian Elements :',6x,f18.6)
            else
               write (iout,310)  sum
  310          format (/,' Sum of Diagonal Hessian Elements :',6x,f16.4)
            end if
         end if
c
c     write out the full matrix of numerical Hessian elements
c
         if (dofull .and. donumer) then
            do i = 1, n
               do k = 1, n
                  write (iout,320)  i,k
  320             format (/,' 3x3 Hessian Block for Atoms :',3x,2i8,/)
                  do j = 1, 3
                     if (digits .ge. 8) then
                        write (iout,330)  (nhess(m,i,j,k),m=1,3)
  330                   format (' Numer',5x,3f20.8)
                     else if (digits .ge. 6) then
                        write (iout,340)  (nhess(m,i,j,k),m=1,3)
  340                   format (' Numer',5x,3f18.6)
                     else
                        write (iout,350)  (nhess(m,i,j,k),m=1,3)
  350                   format (' Numer',5x,3f16.4)
                     end if
                  end do
               end do
            end do
         end if
c
c     write out the full matrix of analytical Hessian elements
c
         if (doanalyt .and. .not.donumer) then
            ihes = freeunit ()
            hessfile = filename(1:leng)//'.hes'
            call version (hessfile,'new')
            open (unit=ihes,file=hessfile,status='new')
            write (iout,360)  hessfile
  360       format (/,' Hessian Matrix written to File :  ',a40)
            write (ihes,370)
  370       format (/,' Diagonal Hessian Elements  (3 per Atom)',/)
            if (digits .ge. 8) then
               write (ihes,380)  ((hdiag(j,i),j=1,3),i=1,n)
  380          format (4f16.8)
            else if (digits .ge. 6) then
               write (ihes,390)  ((hdiag(j,i),j=1,3),i=1,n)
  390          format (5f14.6)
            else
               write (ihes,400)  ((hdiag(j,i),j=1,3),i=1,n)
  400          format (6f12.4)
            end if
            do i = 1, n
               do j = 1, 3
                  if (hinit(j,i) .le. hstop(j,i)) then
                     write (ihes,410)  i,axis(j)
  410                format (/,' Off-diagonal Hessian Elements for Atom'
     &,                         i6,1x,a1,/)
                     if (digits .ge. 8) then
                        write (ihes,420)  (h(k),k=hinit(j,i),hstop(j,i))
  420                   format (4f16.8)
                     else if (digits .ge. 6) then
                        write (ihes,430)  (h(k),k=hinit(j,i),hstop(j,i))
  430                   format (5f14.6)
                     else
                        write (ihes,440)  (h(k),k=hinit(j,i),hstop(j,i))
  440                   format (6f12.4)
                     end if
                  end if
               end do
            end do
            close (unit=ihes)
         end if
c
c     attempt to read next structure from the coordinate file
c
         call readxyz (ixyz)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (hindex)
      deallocate (hinit)
      deallocate (hstop)
      deallocate (h)
      deallocate (g)
      deallocate (g0)
      deallocate (hdiag)
      if (allocated(nhess))  deallocate (nhess)
c
c     perform any final tasks before program exit
c
      call final
      end
