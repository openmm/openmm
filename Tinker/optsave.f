c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine optsave  --  save optimization info and results  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "optsave" is used by the optimizers to write imtermediate
c     coordinates and other relevant information; also checks for
c     user requested termination of an optimization
c
c
      subroutine optsave (ncycle,f,xx)
      use sizes
      use atoms
      use files
      use iounit
      use math
      use omega
      use output
      use scales
      use socket
      use usage
      use zcoord
      implicit none
      integer i,iopt,iend
      integer ncycle,nvar
      integer lext,freeunit
      real*8 f,xx(*)
      logical exist
      character*7 ext
      character*240 optfile
      character*240 endfile
c
c
c     nothing to do if coordinate type is undefined
c
      if (coordtype .eq. 'NONE')  return
c
c     check scaling factors for optimization parameters
c
      if (.not. set_scale) then
         set_scale = .true.
         if (coordtype .eq. 'CARTESIAN') then
            do i = 1, 3*n
               scale(i) = 1.0d0
            end do
         else if (coordtype .eq. 'INTERNAL') then
            do i = 1, nomega
               scale(i) = 1.0d0
            end do
         end if
      end if
c
c     transform optimization parameters back to coordinates
c
      if (coordtype .eq. 'CARTESIAN') then
         nvar = 0
         do i = 1, n
            if (use(i)) then
               nvar = nvar + 1
               x(i) = xx(nvar) / scale(nvar)
               nvar = nvar + 1
               y(i) = xx(nvar) / scale(nvar)
               nvar = nvar + 1
               z(i) = xx(nvar) / scale(nvar)
            end if
         end do
      else if (coordtype .eq. 'INTERNAL') then
         do i = 1, nomega
            dihed(i) = xx(i) / scale(i)
            ztors(zline(i)) = dihed(i) * radian
         end do
      end if
c
c     get name of archive or intermediate coordinates file
c
      iopt = freeunit ()
      if (cyclesave) then
         if (archive) then
            optfile = filename(1:leng)
            call suffix (optfile,'arc','old')
            inquire (file=optfile,exist=exist)
            if (exist) then
               call openend (iopt,optfile)
            else
               open (unit=iopt,file=optfile,status='new')
            end if
         else
            lext = 3
            call numeral (ncycle,ext,lext)
            optfile = filename(1:leng)//'.'//ext(1:lext)
            call version (optfile,'new')
            open (unit=iopt,file=optfile,status='new')
         end if
      else
         optfile = outfile
         call version (optfile,'old')
         open (unit=iopt,file=optfile,status='old')
         rewind (unit=iopt)
      end if
c
c     update intermediate file with desired coordinate type
c
      if (coordtype .eq. 'CARTESIAN') then
         call prtxyz (iopt)
      else if (coordtype .eq. 'INTERNAL') then
         call prtint (iopt)
      else if (coordtype .eq. 'RIGIDBODY') then
         call prtxyz (iopt)
      end if
      close (unit=iopt)
c
c     send data via external socket communication if desired
c
      if (.not.sktstart .or. use_socket) then
         if (coordtype .eq. 'INTERNAL')  call makexyz
         call sktopt (ncycle,f)
      end if
c
c     test for requested termination of the optimization
c
      endfile = 'tinker.end'
      inquire (file=endfile,exist=exist)
      if (.not. exist) then
         endfile = filename(1:leng)//'.end'
         inquire (file=endfile,exist=exist)
         if (exist) then
            iend = freeunit ()
            open (unit=iend,file=endfile,status='old')
            close (unit=iend,status='delete')
         end if
      end if
      if (exist) then
         write (iout,10)
   10    format (/,' OPTSAVE  --  Optimization Calculation Ending',
     &              ' due to User Request')
         call fatal
      end if
      return
      end
