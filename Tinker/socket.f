c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2002  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  module socket  --  socket communication control parameters  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     skttyp      socket information type (1=DYN, 2=OPT)
c     cstep       current dynamics or optimization step number
c     cdt         current dynamics cumulative simulation time
c     cenergy     current potential energy from simulation
c     sktstart    logical flag to indicate socket initialization
c     sktstop     logical flag to indicate socket shutdown
c     use_socket  logical flag governing use of external sockets
c
c
      module socket
      implicit none
      integer skttyp
      integer cstep
      real*8 cdt
      real*8 cenergy
      logical sktstart
      logical sktstop
      logical use_socket
      save
      end
