c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module moment  --  electric multipole moment components  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     netchg   net electric charge for the total system
c     netdpl   dipole moment magnitude for the total system
c     netqdp   diagonal quadrupole (Qxx, Qyy, Qzz) for system
c     xdpl     dipole vector x-component in the global frame
c     ydpl     dipole vector y-component in the global frame
c     zdpl     dipole vector z-component in the global frame
c     xxqdp    quadrupole tensor xx-component in global frame
c     xyqdp    quadrupole tensor xy-component in global frame
c     xzqdp    quadrupole tensor xz-component in global frame
c     yxqdp    quadrupole tensor yx-component in global frame
c     yyqdp    quadrupole tensor yy-component in global frame
c     yzqdp    quadrupole tensor yz-component in global frame
c     zxqdp    quadrupole tensor zx-component in global frame
c     zyqdp    quadrupole tensor zy-component in global frame
c     zzqdp    quadrupole tensor zz-component in global frame
c
c
      module moment
      implicit none
      real*8 netchg,netdpl
      real*8 netqdp(3)
      real*8 xdpl,ydpl,zdpl
      real*8 xxqdp,xyqdp,xzqdp
      real*8 yxqdp,yyqdp,yzqdp
      real*8 zxqdp,zyqdp,zzqdp
      save
      end
