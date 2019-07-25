c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1991  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine torphase  --  torsional amplitude and phase  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "torphase" sets the n-fold amplitude and phase values
c     for each torsion via sorting of the input parameters
c
c
      subroutine torphase (ft,vt,st)
      implicit none
      integer i,k
      integer ft(*)
      real*8 ampli(6)
      real*8 phase(6)
      real*8 vt(*),st(*)
c
c
c     copy the input fold, amplitude and phase angles
c
      do i = 1, 6
         ampli(i) = vt(i)
         phase(i) = st(i)
         vt(i) = 0.0d0
         st(i) = 0.0d0
      end do
c
c     shift the phase angles into the standard range
c
      do i = 1, 6
         do while (phase(i) .lt. -180.0d0)
            phase(i) = phase(i) + 360.0d0
         end do
         do while (phase(i) .gt. 180.0d0)
            phase(i) = phase(i) - 360.0d0
         end do
      end do
c
c     convert input torsional parameters to storage format
c
      do i = 1, 6
         k = ft(i)
         if (k .eq. 0) then
            goto 10
         else if (k .le. 6) then
            vt(k) = ampli(i)
            st(k) = phase(i)
         end if
      end do
   10 continue
      return
      end
