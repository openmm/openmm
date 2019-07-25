c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2000  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine moments  --  total electric multipole moments  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "moments" computes the total electric charge, dipole and
c     quadrupole moments for the active atoms as a sum over the
c     partial charges, bond dipoles and atomic multipole moments
c
c
      subroutine moments
      use sizes
      use atomid
      use atoms
      use bound
      use charge
      use dipole
      use limits
      use moment
      use mpole
      use polar
      use potent
      use rigid
      use solute
      use units
      use usage
      implicit none
      integer i,j,k
      real*8 weigh,qave
      real*8 xc,yc,zc
      real*8 xi,yi,zi,ri
      real*8 xmid,ymid,zmid
      real*8 xbnd,ybnd,zbnd
      real*8, allocatable :: xcm(:)
      real*8, allocatable :: ycm(:)
      real*8, allocatable :: zcm(:)
      real*8 a(3,3),b(3,3)
c
c
c     zero out total charge, dipole and quadrupole components
c
      netchg = 0.0d0
      netdpl = 0.0d0
      netqdp(1) = 0.0d0
      netqdp(2) = 0.0d0
      netqdp(3) = 0.0d0
      xdpl = 0.0d0
      ydpl = 0.0d0
      zdpl = 0.0d0
      xxqdp = 0.0d0
      xyqdp = 0.0d0
      xzqdp = 0.0d0
      yxqdp = 0.0d0
      yyqdp = 0.0d0
      yzqdp = 0.0d0
      zxqdp = 0.0d0
      zyqdp = 0.0d0
      zzqdp = 0.0d0
c
c     maintain periodic boundaries and neighbor lists
c
      if (use_bounds .and. .not.use_rigid)  call bounds
      if (use_clist .or. use_mlist)  call nblist
c
c     perform dynamic allocation of some local arrays
c
      allocate (xcm(n))
      allocate (ycm(n))
      allocate (zcm(n))
c
c     find the center of mass of the set of active atoms
c
      weigh = 0.0d0
      xmid = 0.0d0
      ymid = 0.0d0
      zmid = 0.0d0
      do i = 1, n
         if (use(i)) then
            weigh = weigh + mass(i)
            xmid = xmid + x(i)*mass(i)
            ymid = ymid + y(i)*mass(i)
            zmid = zmid + z(i)*mass(i)
         end if
      end do
      if (weigh .ne. 0.0d0) then
         xmid = xmid / weigh
         ymid = ymid / weigh
         zmid = zmid / weigh
      end if
      do i = 1, n
         xcm(i) = x(i) - xmid
         ycm(i) = y(i) - ymid
         zcm(i) = z(i) - zmid
      end do
c
c     set the multipole moment components due to partial charges
c
      do i = 1, nion
         k = iion(i)
         if (use(k)) then
            netchg = netchg + pchg(i)
            xdpl = xdpl + xcm(k)*pchg(i)
            ydpl = ydpl + ycm(k)*pchg(i)
            zdpl = zdpl + zcm(k)*pchg(i)
            xxqdp = xxqdp + xcm(k)*xcm(k)*pchg(i)
            xyqdp = xyqdp + xcm(k)*ycm(k)*pchg(i)
            xzqdp = xzqdp + xcm(k)*zcm(k)*pchg(i)
            yxqdp = yxqdp + ycm(k)*xcm(k)*pchg(i)
            yyqdp = yyqdp + ycm(k)*ycm(k)*pchg(i)
            yzqdp = yzqdp + ycm(k)*zcm(k)*pchg(i)
            zxqdp = zxqdp + zcm(k)*xcm(k)*pchg(i)
            zyqdp = zyqdp + zcm(k)*ycm(k)*pchg(i)
            zzqdp = zzqdp + zcm(k)*zcm(k)*pchg(i)
         end if
      end do
c
c     set the multipole moment components due to bond dipoles
c
      do i = 1, ndipole
         j = idpl(1,i)
         k = idpl(2,i)
         if (use(j) .or. use(k)) then
            xi = x(j) - x(k)
            yi = y(j) - y(k)
            zi = z(j) - z(k)
            ri = sqrt(xi*xi + yi*yi + zi*zi)
            xbnd = bdpl(i) * (xi/ri) / debye
            ybnd = bdpl(i) * (yi/ri) / debye
            zbnd = bdpl(i) * (zi/ri) / debye
            xc = x(j) - xi*sdpl(i)
            yc = y(j) - yi*sdpl(i)
            zc = z(j) - zi*sdpl(i)
            xdpl = xdpl + xbnd
            ydpl = ydpl + ybnd
            zdpl = zdpl + zbnd
            xxqdp = xxqdp + 2.0d0*xc*xbnd
            xyqdp = xyqdp + xc*ybnd + yc*xbnd
            xzqdp = xzqdp + xc*zbnd + zc*xbnd
            yxqdp = yxqdp + yc*xbnd + xc*ybnd
            yyqdp = yyqdp + 2.0d0*yc*ybnd
            yzqdp = yzqdp + yc*zbnd + zc*ybnd
            zxqdp = zxqdp + zc*xbnd + xc*zbnd
            zyqdp = zyqdp + zc*ybnd + yc*zbnd
            zzqdp = zzqdp + 2.0d0*zc*zbnd
         end if
      end do
c
c     find atomic multipoles and induced dipoles in global frame
c
      if (use_born)  call born
      call chkpole
      call rotpole
      call induce
      if (solvtyp.eq.'GK' .or. solvtyp.eq.'PB') then
         do i = 1, npole
            rpole(2,i) = rpole(2,i) + uinds(1,i)
            rpole(3,i) = rpole(3,i) + uinds(2,i)
            rpole(4,i) = rpole(4,i) + uinds(3,i)
         end do
      else
         do i = 1, npole
            rpole(2,i) = rpole(2,i) + uind(1,i)
            rpole(3,i) = rpole(3,i) + uind(2,i)
            rpole(4,i) = rpole(4,i) + uind(3,i)
         end do
      end if
c
c     set the multipole moment components due to atomic multipoles
c
      do i = 1, npole
         k = ipole(i)
         if (use(k)) then
            netchg = netchg + rpole(1,i)
            xdpl = xdpl + xcm(k)*rpole(1,i) + rpole(2,i)
            ydpl = ydpl + ycm(k)*rpole(1,i) + rpole(3,i)
            zdpl = zdpl + zcm(k)*rpole(1,i) + rpole(4,i)
            xxqdp = xxqdp + xcm(k)*xcm(k)*rpole(1,i)
     &                 + 2.0d0*xcm(k)*rpole(2,i)
            xyqdp = xyqdp + xcm(k)*ycm(k)*rpole(1,i)
     &                 + xcm(k)*rpole(3,i) + ycm(k)*rpole(2,i)
            xzqdp = xzqdp + xcm(k)*zcm(k)*rpole(1,i)
     &                 + xcm(k)*rpole(4,i) + zcm(k)*rpole(2,i)
            yxqdp = yxqdp + ycm(k)*xcm(k)*rpole(1,i)
     &                 + ycm(k)*rpole(2,i) + xcm(k)*rpole(3,i)
            yyqdp = yyqdp + ycm(k)*ycm(k)*rpole(1,i)
     &                 + 2.0d0*ycm(k)*rpole(3,i)
            yzqdp = yzqdp + ycm(k)*zcm(k)*rpole(1,i)
     &                 + ycm(k)*rpole(4,i) + zcm(k)*rpole(3,i)
            zxqdp = zxqdp + zcm(k)*xcm(k)*rpole(1,i)
     &                 + zcm(k)*rpole(2,i) + xcm(k)*rpole(4,i)
            zyqdp = zyqdp + zcm(k)*ycm(k)*rpole(1,i)
     &                 + zcm(k)*rpole(3,i) + ycm(k)*rpole(4,i)
            zzqdp = zzqdp + zcm(k)*zcm(k)*rpole(1,i)
     &                 + 2.0d0*zcm(k)*rpole(4,i)
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xcm)
      deallocate (ycm)
      deallocate (zcm)
c
c     convert the quadrupole from traced to traceless form
c
      qave = (xxqdp + yyqdp + zzqdp) / 3.0d0
      xxqdp = 1.5d0 * (xxqdp-qave)
      xyqdp = 1.5d0 * xyqdp
      xzqdp = 1.5d0 * xzqdp
      yxqdp = 1.5d0 * yxqdp
      yyqdp = 1.5d0 * (yyqdp-qave)
      yzqdp = 1.5d0 * yzqdp
      zxqdp = 1.5d0 * zxqdp
      zyqdp = 1.5d0 * zyqdp
      zzqdp = 1.5d0 * (zzqdp-qave)
c
c     add the traceless atomic quadrupoles to total quadrupole
c
      do i = 1, npole
         k = ipole(i)
         if (use(k)) then
            xxqdp = xxqdp + 3.0d0*rpole(5,i)
            xyqdp = xyqdp + 3.0d0*rpole(6,i)
            xzqdp = xzqdp + 3.0d0*rpole(7,i)
            yxqdp = yxqdp + 3.0d0*rpole(8,i)
            yyqdp = yyqdp + 3.0d0*rpole(9,i)
            yzqdp = yzqdp + 3.0d0*rpole(10,i)
            zxqdp = zxqdp + 3.0d0*rpole(11,i)
            zyqdp = zyqdp + 3.0d0*rpole(12,i)
            zzqdp = zzqdp + 3.0d0*rpole(13,i)
         end if
      end do
c
c     convert dipole to Debyes and quadrupole to Buckinghams
c
      xdpl = xdpl * debye
      ydpl = ydpl * debye
      zdpl = zdpl * debye
      xxqdp = xxqdp * debye
      xyqdp = xyqdp * debye
      xzqdp = xzqdp * debye
      yxqdp = yxqdp * debye
      yyqdp = yyqdp * debye
      yzqdp = yzqdp * debye
      zxqdp = zxqdp * debye
      zyqdp = zyqdp * debye
      zzqdp = zzqdp * debye
c
c     get dipole magnitude and diagonalize quadrupole tensor
c
      netdpl = sqrt(xdpl*xdpl + ydpl*ydpl + zdpl*zdpl)
      a(1,1) = xxqdp
      a(1,2) = xyqdp
      a(1,3) = xzqdp
      a(2,1) = yxqdp
      a(2,2) = yyqdp
      a(2,3) = yzqdp
      a(3,1) = zxqdp
      a(3,2) = zyqdp
      a(3,3) = zzqdp
      call jacobi (3,a,netqdp,b)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine momfull  --  multipole moments for full system  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "momfull" computes the electric moments for the full system
c     as a sum over the partial charges, bond dipoles and atomic
c     multipole moments
c
c
      subroutine momfull
      use sizes
      use atoms
      use usage
      implicit none
      integer i
      logical, allocatable :: temp(:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (temp(n))
c
c     store active atom list, and make all atoms active
c
      if (nuse .ne. n) then
         do i = 1, n
            temp(i) = use(i)
            use(i) = .true.
         end do
      end if
c
c     compute the electric multipole moments for the system
c
      call moments
c
c     revert to the original set of active atoms
c
      if (nuse .ne. n) then
         do i = 1, n
            use(i) = temp(i)
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (temp)
      return
      end
