c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2006 by Michael Schnieders & Jay W. Ponder  ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  module neigh  --  pairwise neighbor list indices & storage  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     maxvlst     maximum size of van der Waals pair neighbor lists
c     maxelst     maximum size of electrostatic pair neighbor lists
c     maxulst     maximum size of dipole preconditioner pair lists
c     nvlst       number of sites in list for each vdw site
c     vlst        site numbers in neighbor list of each vdw site
c     nelst       number of sites in list for each electrostatic site
c     elst        site numbers in list of each electrostatic site
c     nulst       number of sites in list for each preconditioner site
c     ulst        site numbers in list of each preconditioner site
c     lbuffer     width of the neighbor list buffer region
c     pbuffer     width of the preconditioner list buffer region
c     lbuf2       square of half the neighbor list buffer width
c     pbuf2       square of half the preconditioner list buffer width
c     vbuf2       square of vdw cutoff plus neighbor list buffer
c     cbuf2       square of charge cutoff plus neighbor list buffer
c     mbuf2       square of multipole cutoff plus neighbor list buffer
c     ubuf2       square of preconditioner cutoff plus neighbor buffer
c     vbufx       square of vdw cutoff plus twice the list buffer
c     cbufx       square of charge cutoff plus twice the list buffer
c     mbufx       square of multipole cutoff plus twice the list buffer
c     ubufx       square of preconditioner cutoff plus twice the buffer
c     xvold       x-coordinate at last vdw neighbor list update
c     yvold       y-coordinate at last vdw neighbor list update
c     zvold       z-coordinate at last vdw neighbor list update
c     xcold       x-coordinate at last charge neighbor list update
c     ycold       y-coordinate at last charge neighbor list update
c     zcold       z-coordinate at last charge neighbor list update
c     xmold       x-coordinate at last multipole neighbor list update
c     ymold       y-coordinate at last multipole neighbor list update
c     zmold       z-coordinate at last multipole neighbor list update
c     xuold       x-coordinate at last preconditioner neighbor update
c     yuold       y-coordinate at last preconditioner neighbor update
c     zuold       z-coordinate at last preconditioner neighbor update
c     dovlst      logical flag to rebuild vdw neighbor list
c     doclst      logical flag to rebuild charge neighbor list
c     domlst      logical flag to rebuild multipole neighbor list
c     doulst      logical flag to rebuild preconditioner neighbor list
c
c
      module neigh
      implicit none
      integer maxvlst
      integer maxctlst 
      integer maxelst
      integer maxulst
      integer, allocatable :: nvlst(:)
      integer, allocatable :: vlst(:,:)
      integer, allocatable :: nelst(:)
      integer, allocatable :: elst(:,:)
      integer, allocatable :: nulst(:)
      integer, allocatable :: ulst(:,:)
      integer, allocatable :: nctlst(:)
      integer, allocatable :: ctlst(:,:)
      real*8 lbuffer,pbuffer
      real*8 lbuf2,pbuf2
      real*8 vbuf2,cbuf2
      real*8 mbuf2,ubuf2
      real*8 vbufx,cbufx
      real*8 mbufx,ubufx
      real*8 ctbuf2,ctbufx !CT
      real*8, allocatable :: xvold(:)
      real*8, allocatable :: yvold(:)
      real*8, allocatable :: zvold(:)
      real*8, allocatable :: xctold(:) !CT
      real*8, allocatable :: yctold(:) !CT
      real*8, allocatable :: zctold(:) !CT
      real*8, allocatable :: xcold(:)
      real*8, allocatable :: ycold(:)
      real*8, allocatable :: zcold(:)
      real*8, allocatable :: xmold(:)
      real*8, allocatable :: ymold(:)
      real*8, allocatable :: zmold(:)
      real*8, allocatable :: xuold(:)
      real*8, allocatable :: yuold(:)
      real*8, allocatable :: zuold(:)
      logical dovlst,doclst
      logical domlst,doulst
      logical doctlst !CT
      save
      end
