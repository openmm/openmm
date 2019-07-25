c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2000  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine openend  --  open a file positioned for append  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "openend" opens a file on a Fortran unit such that the position
c     is set to the bottom for appending to the end of the file
c
c     note this routine is system dependent since the Fortran 90
c     standard is not supported by many Fortran 77 compilers; only
c     one of the various implementations below should be activated
c     by removing comment characters
c
c
      subroutine openend (iunit,name)
      implicit none
      integer iunit
      character*240 name
c
c
c     standard Fortran 90, unavailable in some Fortran 77 compilers
c
      open (unit=iunit,file=name,status='old',position='append')
c
c     common extension supported by many Fortran 77 compilers
c
c     open (unit=iunit,file=name,status='old',access='append')
c
c     some Fortran 77 compilers open files for append by default
c
c     open (unit=iunit,file=name,status='old')
c
c     manually read to the end of file, slow but always correct
c
c     open (unit=iunit,file=name,status='old')
c     do while (.true.)
c        read (iunit,10,err=20,end=20)
c  10    format ()
c     end do
c  20 continue
      return
      end
