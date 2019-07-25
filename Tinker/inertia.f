c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine inertia  --  principal moments of inertia  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "inertia" computes the principal moments of inertia for the
c     system, and optionally translates the center of mass to the
c     origin and rotates the principal axes onto the global axes
c
c        mode = 1     print the moments and principal axes
c        mode = 2     move coordinates to standard orientation
c        mode = 3     perform both of the above operations
c
c     literature reference:
c
c     Herbert Goldstein, "Classical Mechanics, 2nd Edition",
c     Addison-Wesley, Reading, MA, 1980; see the Euler angle
c     xyz convention in Appendix B
c
c
      subroutine inertia (mode)
      use sizes
      use atoms
      use atomid
      use iounit
      use math
      implicit none
      integer i,j,k,mode
      real*8 weigh,total,dot
      real*8 xcm,ycm,zcm
      real*8 xx,xy,xz,yy,yz,zz
      real*8 xterm,yterm,zterm
      real*8 phi,theta,psi
      real*8 moment(3),vec(3,3)
      real*8 tensor(3,3),a(3,3)
      logical print,move
c
c
c     decide upon the type of output desired
c
      print = .false.
      move = .false.
      if (mode.eq.1 .or. mode.eq.3)  print = .true.
      if (mode.eq.2 .or. mode.eq.3)  move = .true.
c
c     compute the position of the center of mass
c
      total = 0.0d0
      xcm = 0.0d0
      ycm = 0.0d0
      zcm = 0.0d0
      do i = 1, n
         weigh = mass(i)
         total = total + weigh
         xcm = xcm + x(i)*weigh
         ycm = ycm + y(i)*weigh
         zcm = zcm + z(i)*weigh
      end do
      xcm = xcm / total
      ycm = ycm / total
      zcm = zcm / total
c
c     compute and then diagonalize the inertia tensor
c
      xx = 0.0d0
      xy = 0.0d0
      xz = 0.0d0
      yy = 0.0d0
      yz = 0.0d0
      zz = 0.0d0
      do i = 1, n
         weigh = mass(i)
         xterm = x(i) - xcm
         yterm = y(i) - ycm
         zterm = z(i) - zcm
         xx = xx + xterm*xterm*weigh
         xy = xy + xterm*yterm*weigh
         xz = xz + xterm*zterm*weigh
         yy = yy + yterm*yterm*weigh
         yz = yz + yterm*zterm*weigh
         zz = zz + zterm*zterm*weigh
      end do
      tensor(1,1) = yy + zz
      tensor(2,1) = -xy
      tensor(3,1) = -xz
      tensor(1,2) = -xy
      tensor(2,2) = xx + zz
      tensor(3,2) = -yz
      tensor(1,3) = -xz
      tensor(2,3) = -yz
      tensor(3,3) = xx + yy
      call jacobi (3,tensor,moment,vec)
c
c     select the direction for each principal moment axis
c
      do i = 1, 2
         do j = 1, n
            xterm = vec(1,i) * (x(j)-xcm)
            yterm = vec(2,i) * (y(j)-ycm)
            zterm = vec(3,i) * (z(j)-zcm)
            dot = xterm + yterm + zterm
            if (dot .lt. 0.0d0) then
               do k = 1, 3
                  vec(k,i) = -vec(k,i)
               end do
            end if
            if (dot .ne. 0.0d0)  goto 10
         end do
   10    continue
      end do
c
c     moment axes must give a right-handed coordinate system
c
      xterm = vec(1,1) * (vec(2,2)*vec(3,3)-vec(2,3)*vec(3,2))
      yterm = vec(2,1) * (vec(1,3)*vec(3,2)-vec(1,2)*vec(3,3))
      zterm = vec(3,1) * (vec(1,2)*vec(2,3)-vec(1,3)*vec(2,2))
      dot = xterm + yterm + zterm
      if (dot .lt. 0.0d0) then
         do j = 1, 3
            vec(j,3) = -vec(j,3)
         end do
      end if
c
c     principal moment axes form rows of Euler rotation matrix
c
      if (move) then
         do i = 1, 3
            do j = 1, 3
               a(i,j) = vec(j,i)
            end do
         end do
c
c     translate to origin, then apply Euler rotation matrix
c
         do i = 1, n
            xterm = x(i) - xcm
            yterm = y(i) - ycm
            zterm = z(i) - zcm
            x(i) = a(1,1)*xterm + a(1,2)*yterm + a(1,3)*zterm
            y(i) = a(2,1)*xterm + a(2,2)*yterm + a(2,3)*zterm
            z(i) = a(3,1)*xterm + a(3,2)*yterm + a(3,3)*zterm
         end do
      end if
c
c     print the center of mass and Euler angle values
c
      if (print) then
         write (iout,20)  xcm,ycm,zcm
   20    format (/,' Center of Mass Coordinates :',7x,3f13.6)
         call invert (3,vec)
         call roteuler (vec,phi,theta,psi)
         phi = radian * phi
         theta = radian * theta
         psi = radian * psi
         write (iout,30)  phi,theta,psi
   30    format (' Euler Angles (Phi/Theta/Psi) : ',4x,3f13.3)
c
c     print the moments of inertia and the principal axes
c
         write (iout,40)
   40    format (/,' Moments of Inertia and Principal Axes :',
     &           //,13x,'Moments (amu Ang^2)',
     &              12x,'X-, Y- and Z-Components of Axes')
         write (iout,50)  (moment(i),vec(1,i),vec(2,i),vec(3,i),i=1,3)
   50    format (3(/,11x,f16.3,9x,3f13.6))
      end if
      return
      end
