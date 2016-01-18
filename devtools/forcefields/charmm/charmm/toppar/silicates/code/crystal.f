c
c
c     "crystal" is a utility program which converts between
c     fractional and Cartesian coordinates, and can generate
c     full unit cells from asymmetric units
c
c
      subroutine crystal(n,name,valence,x,y,z)
      implicit none

      integer mxatm
      parameter (mxatm=1000)
      integer n
      integer iout
      integer valence(mxatm)
      double precision radian
      parameter (radian=57.295780d0)
      double precision alpha,beta,gamma
      double precision adim,bdim,cdim
      double precision aleng(3),bleng(3),cleng(3) 
      double precision x(mxatm),y(mxatm),z(mxatm)
      logical orthogonal,triclinic
      character*4 name(mxatm)
c
      call cell(adim,bdim,cdim,alpha,beta,gamma)
c
c
      iout = 6

      orthogonal = .false.
      triclinic = .false.

      if (alpha.eq. 90.0d0 .and. beta .eq. 90.0d0
     $    .and. gamma.eq.90.0d0) then
         orthogonal = .true.
      else
         triclinic = .true.
      end if
c
c     determina as coordenadas de cada crystal system
c
      if (orthogonal) then
         aleng(1) = adim
         aleng(2) = 0.0d0
         aleng(3) = 0.0d0
         bleng(1) = 0.0d0
         bleng(2) = bdim
         bleng(3) = 0.0d0
         cleng(1) = 0.0d0
         cleng(2) = 0.0d0
         cleng(3) = cdim
      else if (triclinic) then
         aleng(1) = adim*sin(alpha/radian)
         aleng(2) = adim*cos(alpha/radian)
         aleng(3) = 0.0d0
         bleng(1) = bdim*cos(gamma/radian)
         bleng(2) = bdim*sin(gamma/radian)
         bleng(3) = 0.0d0
         cleng(1) = 0.0d0
         cleng(2) = 0.0d0
         cleng(3) = cdim
      end if
c
c     print out the initial cell dimensions to be used
c
      write (iout,10) adim,bdim,cdim,alpha,beta,gamma
  10  format (/,'Unit Cell Dimensions :'
     $        /,'      a    =',f10.4,
     $        /,'      b    =',f10.4,
     $        /,'      c    =',f10.4,
     $        /,'      alpha =',f10.4,
     $        /,'      beta  =',f10.4,
     $        /,'      gamma =',f10.4)
c
c     replicate the unit cell to make a block of unit cells
c
      call replicate 
     $       (n,name,valence,x,y,z,aleng,bleng,cleng)
c
c
      return
      end
c
      subroutine replicate 
     $           (n,name,valence,x,y,z,aleng,bleng,cleng)
      implicit none

      integer mxatm
      parameter (mxatm=1000)
      integer i,n,nunit   !nunit number of atoms in unit cell  PEML
      integer valence(mxatm)
      character*4 name(mxatm)
      double precision aleng(3),bleng(3),cleng(3)
      double precision x(mxatm),y(mxatm),z(mxatm)
c
c     translate along XX
c
      nunit = n
      do i = 1, nunit
         n = n + 1
         x(n) = x(i) + aleng(1)
         y(n) = y(i) + aleng(2)
         z(n) = z(i) + aleng(3)

         name(n) = name(i)
         valence(n) = valence(i)

      end do
c
c     translate along YY
c
      do i = 1, nunit
         n = n + 1
         x(n) = x(i) + bleng(1)
         y(n) = y(i) + bleng(2)
         z(n) = z(i) + bleng(3)
         name(n) = name(i)
         valence(n) = valence(i)
      end do
c
c     translate along XY
c
      do i = 1, nunit
         n = n + 1
         x(n) = x(i) + aleng(1) + bleng(1)
         y(n) = y(i) + aleng(2) + bleng(2)
         z(n) = z(i) + aleng(3) + bleng(3)
         name(n) = name(i)
         valence(n) = valence(i)
      end do
      
      return
      end
