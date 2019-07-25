c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module math  --  mathematical and geometrical constants  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     pi        numerical value of the geometric constant
c     sqrtpi    numerical value of the square root of Pi
c     radian    conversion factor from radians to degrees
c     logten    numerical value of the natural log of ten
c     sqrttwo   numerical value of the square root of two
c     twosix    numerical value of the sixth root of two
c
c
      module math
      implicit none
      real*8 pi,sqrtpi
      real*8 radian,logten
      real*8 sqrttwo,twosix
      parameter (pi=3.141592653589793238d0)
      parameter (sqrtpi=1.772453850905516027d0)
      parameter (radian=57.29577951308232088d0)
      parameter (logten=2.302585092994045684d0)
      parameter (sqrttwo=1.414213562373095049d0)
      parameter (twosix=1.122462048309372981d0)
      save
      end
