c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2006 by Michael Schnieders & Jay W. Ponder  ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  routines below implement dummy versions of the APBS   ##
c     ##  calls required for TINKER to interface with the APBS  ##
c     ##  Poisson-Boltzmann solver package from Nathan Baker    ##
c     ##                                                        ##
c     ############################################################
c
c     ##############################
c     ##                          ##
c     ##  subroutine apbsinitial  ##
c     ##                          ##
c     ##############################
c
c
      subroutine apbsinitial (dime,grid,gcent,cgrid,cgcent,fgrid,
     &                        fgcent,pdie,sdie,srad,swin,sdens,
     &                        kelvin,ionn,ionc,ionq,ionr,pbtyp,
     &                        pbtyplen,pbsoln,pbsolnlen,bcfl,
     &                        bcfllen,chgm,chgmlen,srfm,srfmlen)
      use iounit
      implicit none
      integer dime(*)
      integer ionn
      integer ionq(*)
      integer pbtyplen
      integer pbsolnlen
      integer bcfllen
      integer chgmlen
      integer srfmlen
      real*8 grid(*)
      real*8 gcent(*)
      real*8 cgrid(*)
      real*8 cgcent(*)
      real*8 fgrid(*)
      real*8 fgcent(*)
      real*8 pdie
      real*8 sdie
      real*8 srad
      real*8 swin
      real*8 sdens
      real*8 kelvin
      real*8 ionc(*)
      real*8 ionr(*)
      character*(*) pbtyp
      character*(*) pbsoln
      character*(*) bcfl
      character*(*) chgm
      character*(*) srfm
c
c
c     exit with an error message if APBS calculation is attempted
c
      write (iout,10)
   10 format (/,' APBSINITIAL  --  APBS Not Supported by This',
     &           ' TINKER Version')
      call fatal
      return
      end
c
c
c     #############################
c     ##                         ##
c     ##  subroutine apbsempole  ##
c     ##                         ##
c     #############################
c
c
      subroutine apbsempole (n,pos,rsolv,pbpole,pbe,apbe,pbep,pbfp,pbtp)
      implicit none
      integer n
      real*8 pos(*)
      real*8 rsolv(*)
      real*8 pbpole(*)
      real*8 pbe
      real*8 apbe(*)
      real*8 pbep(*)
      real*8 pbfp(*)
      real*8 pbtp(*)
      return
      end
c
c
c     #############################
c     ##                         ##
c     ##  subroutine apbsinduce  ##
c     ##                         ##
c     #############################
c
c
      subroutine apbsinduce (indpole,pbeuind)
      implicit none
      real*8 indpole(*)
      real*8 pbeuind(*)
      return
      end
c
c
c     ###############################
c     ##                           ##
c     ##  subroutine apbsnlinduce  ##
c     ##                           ##
c     ###############################
c
c
      subroutine apbsnlinduce (inppole,pbeuinp)
      implicit none
      real*8 inppole(*)
      real*8 pbeuinp(*)
      return
      end
c
c
c     ###################################
c     ##                               ##
c     ##  subroutine pbdirectpolforce  ##
c     ##                               ##
c     ###################################
c
c
      subroutine pbdirectpolforce (indpole,inppole,directf,directt)
      implicit none
      real*8 indpole(*)
      real*8 inppole(*)
      real*8 directf(*)
      real*8 directt(*)
      return
      end
c
c
c     ###################################
c     ##                               ##
c     ##  subroutine pbmutualpolforce  ##
c     ##                               ##
c     ###################################
c
c
      subroutine pbmutualpolforce (indpole,inppole,mutualf)
      implicit none
      real*8 indpole(*)
      real*8 inppole(*)
      real*8 mutualf(*)
      return
      end
c
c
c     ############################
c     ##                        ##
c     ##  subroutine apbsfinal  ##
c     ##                        ##
c     ############################
c
c
      subroutine apbsfinal
      implicit none
      return
      end
