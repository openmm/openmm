c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  program timerot  --  timer for torsional energy terms  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "timerot" measures the CPU time required for file reading
c     and parameter assignment, potential energy computation,
c     energy and gradient over torsions, and torsional angle
c     Hessian matrix evaluation
c
c
      program timerot
      use sizes
      use iounit
      use limits
      use omega
      implicit none
      integer i,ncalls,next
      real*8 energy,value
      real*8 wall,cpu
      real*8, allocatable :: derivs(:)
      real*8, allocatable :: hrot(:,:)
      logical exist,query
      logical dohessian
      character*1 answer
      character*240 record
      character*240 string
c
c
c     read in the molecular system to be timed
c
      call initial
      call getint
c
c     get the timing for setup of the calculation
c
      call settime
      call mechanic
      call initrot
      if (use_list)  call nblist
      call gettime (wall,cpu)
c
c     get the number of calculation cycles to perform
c
      ncalls = 0
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  ncalls
         query = .false.
      end if
   10 continue
      if (query) then
         write (iout,20)
   20    format (/,' Enter Desired Number of Repetitions [1] :  ',$)
         read (input,30)  ncalls
   30    format (i10)
      end if
      if (ncalls .eq. 0)  ncalls = 1
c
c     decide whether to include timing of Hessian evaluations
c
      dohessian = .false.
      if (nomega .le. 1000) then
         call nextarg (answer,exist)
         if (.not. exist) then
            write (iout,40)
   40       format (/,' Include Timing for Hessian Evaluations',
     &                 ' [N] :  ',$)
            read (input,50)  record
   50       format (a240)
            next = 1
            call gettext (record,answer,next)
         end if
         call upcase (answer)
         if (answer .eq. 'Y')  dohessian = .true.
      end if
c
c     print the time required for the computation setup
c
      write (iout,60)  ncalls
   60 format (/,' Total Wall Clock and CPU Time in Seconds for',
     &           i6,' Evaluations :')
      write (iout,70)  wall,cpu
   70 format (/,' Computation Set-up :',f15.3,' Sec (Wall)',
     &           f15.3,' Sec (CPU)')
c
c     run the potential energy only timing experiment
c
      call settime
      do i = 1, ncalls
         value = energy ()
      end do
      call gettime (wall,cpu)
      write (iout,80)  wall,cpu
   80 format (/,' Potential Energy :  ',f15.3,' Sec (Wall)',
     &           f15.3,' Sec (CPU)')
c
c     perform dynamic allocation of some local arrays
c
      allocate (derivs(nomega))
      allocate (hrot(nomega,nomega))
c
c     run the energy and gradient timing experiment
c
      call settime
      do i = 1, ncalls
         call gradrot (value,derivs)
      end do
      call gettime (wall,cpu)
      write (iout,90)  wall,cpu
   90 format (/,' Energy & Gradient : ',f15.3,' Sec (Wall)',
     &           f15.3,' Sec (CPU)')
c
c     run the Hessian matrix only timing experiment
c
      if (dohessian) then
         call settime
         do i = 1, ncalls
            call hessrot ('FULL',hrot)
         end do
         call gettime (wall,cpu)
         write (iout,100)  wall,cpu
  100    format (/,' Hessian Matrix :    ',f15.3,' Sec (Wall)',
     &              f15.3,' Sec (CPU)')
      end if
c
c     repeat the potential energy only timing experiment
c
      call settime
      do i = 1, ncalls
         value = energy ()
      end do
      call gettime (wall,cpu)
      write (iout,110)  wall,cpu
  110 format (/,' Potential Energy :  ',f15.3,' Sec (Wall)',
     &           f15.3,' Sec (CPU)')
c
c     repeat the energy and gradient timing experiment
c
      call settime
      do i = 1, ncalls
         call gradrot (value,derivs)
      end do
      call gettime (wall,cpu)
      write (iout,120)  wall,cpu
  120 format (/,' Energy & Gradient : ',f15.3,' Sec (Wall)',
     &           f15.3,' Sec (CPU)')
c
c     repeat the Hessian matrix only timing experiment
c
      if (dohessian) then
         call settime
         do i = 1, ncalls
            call hessrot ('FULL',hrot)
         end do
         call gettime (wall,cpu)
         write (iout,130)  wall,cpu
  130    format (/,' Hessian Matrix :    ',f15.3,' Sec (Wall)',
     &              f15.3,' Sec (CPU)')
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (derivs)
      deallocate (hrot)
c
c     perform any final tasks before program exit
c
      call final
      end
