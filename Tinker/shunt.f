c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  module shunt  --  polynomial switching function values  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     off    distance at which the potential energy goes to zero
c     off2   square of distance at which the potential goes to zero
c     cut    distance at which switching of the potential begins
c     cut2   square of distance at which the switching begins
c     c0     zeroth order coefficient of multiplicative switch
c     c1     first order coefficient of multiplicative switch
c     c2     second order coefficient of multiplicative switch
c     c3     third order coefficient of multiplicative switch
c     c4     fourth order coefficient of multiplicative switch
c     c5     fifth order coefficient of multiplicative switch
c     f0     zeroth order coefficient of additive switch function
c     f1     first order coefficient of additive switch function
c     f2     second order coefficient of additive switch function
c     f3     third order coefficient of additive switch function
c     f4     fourth order coefficient of additive switch function
c     f5     fifth order coefficient of additive switch function
c     f6     sixth order coefficient of additive switch function
c     f7     seventh order coefficient of additive switch function
c
c
      module shunt
      implicit none
      real*8 off,off2
      real*8 cut,cut2
      real*8 c0,c1,c2
      real*8 c3,c4,c5
      real*8 f0,f1,f2,f3
      real*8 f4,f5,f6,f7
      save
      end
