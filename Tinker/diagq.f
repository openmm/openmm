c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine diagq  --  fast matrix diagonalization routine  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "diagq" is a matrix diagonalization routine which is derived
c     from the classical given, housec, and eigen algorithms with
c     several modifications to increase efficiency and accuracy
c
c     variables and parameters:
c
c     n         dimension of the matrix to be diagonalized
c     nv        number of eigenvalues and eigenvectors desired
c     dd        upper triangle of the matrix to be diagonalized
c     ev        returned with the eigenvalues in ascending order
c     vec       returned with the eigenvectors of the matrix
c     a,b,p,w   local vectors containing temporary work space
c     ta,tb,y   local vectors containing temporary work space
c
c     literature reference:
c
c     adapted from an original routine written by Bernie Brooks,
c     National Institutes of Health, Bethesda, MD
c
c
      subroutine diagq (n,nv,dd,ev,vec)
      implicit none
      integer i,j,k,m,n
      integer ia,ii,ji
      integer mi,mj,mk
      integer nn,nv,ntot
      integer nom,nomtch
      integer ipt,iter,j1
      integer mi1,mj1,mk1
      real*8 alimit,anorm
      real*8 eta,theta,eps
      real*8 rho,delta,gamma
      real*8 zeta,sigma
      real*8 bx,elim1,elim2
      real*8 epr,f0,factor
      real*8 rvalue,rand1
      real*8 rpower,rpow1
      real*8 root,rootx
      real*8 s,sgn,sum1
      real*8 t,term,temp
      real*8 trial,xnorm
      real*8 dd(*)
      real*8 ev(*)
      real*8, allocatable :: a(:)
      real*8, allocatable :: b(:)
      real*8, allocatable :: p(:)
      real*8, allocatable :: w(:)
      real*8, allocatable :: ta(:)
      real*8, allocatable :: tb(:)
      real*8, allocatable :: y(:)
      real*8 vec(n,*)
      logical done
c
c
c     initialization of some necessary parameters
c
      eta = 1.0d-16
      theta = 1.0d37
      eps = 100.0d0 * eta
      rho = eta / 100.0d0
      delta = 100.0d0 * eta**2
      gamma = eta**2 / 100.0d0
      zeta = 1000.0d0 / theta
      sigma = theta * delta / 1000.0d0
      rvalue = 4099.0d0
      rpower = 8388608.0d0
      rpow1 = 0.5d0 * rpower
      rand1 = rpower - 3.0d0
c
c     perform dynamic allocation of some local arrays
c
      allocate (a(n))
      allocate (b(n))
      allocate (p(n))
      allocate (w(n))
      allocate (ta(n))
      allocate (tb(n))
      allocate (y(n))
c
c     get norm of the input matrix and scale
c
      factor = 0.0d0
      ntot = n*(n+1) / 2
      do i = 1, ntot
         factor = max(factor,abs(dd(i)))
      end do
      if (factor .eq. 0.0d0)  return
      k = 0
      anorm = 0.0d0
      do i = 1, n
         do j = i, n
            k = k + 1
            term = (dd(k)/factor)**2
            if (i .eq. j)  term = 0.5d0 * term
            anorm = anorm + term
         end do
      end do
      anorm = factor * sqrt(2.0d0*anorm)
      do i = 1, ntot
         dd(i) = dd(i) / anorm
      end do
c
c     compute the tridiagonalization of the matrix
c
      nn = n - 1
      mi = 0
      mi1 = n - 1
      do i = 1, nn
         sum1 = 0.0d0
         b(i) = 0.0d0
         ji = i + 1
         ipt = mi + i
         a(i) = dd(ipt)
         ipt = ipt + 1
         bx = dd(ipt)
         do j = ji+1, n
            ipt = ipt + 1
            sum1 = sum1 + dd(ipt)*dd(ipt)
         end do
         if (sum1 .lt. gamma) then
            b(i) = bx
            dd(mi+ji) = 0.0d0
         else
            s = sqrt(sum1+bx*bx)
            sgn = 1.0d0
            if (bx .lt. 0.0)  sgn = -1.0d0
            temp = abs(bx)
            w(ji) = sqrt(0.5d0*(1.0d0+(temp/s)))
            ipt = mi + ji
            dd(ipt) = w(ji)
            ii = i + 2
            if (ii .le. n) then
               temp = sgn / (2.0d0*w(ji)*s)
               do j = ii, n
                  ipt = ipt + 1
                  w(j) = temp * dd(ipt)
                  dd(ipt) = w(j)
               end do
            end if
            b(i) = -sgn * s
            do j = ji, n
               p(j) = 0.0d0
            end do
            mk = mi + mi1
            mk1 = mi1 - 1
            do k = ji, n
               ipt = mk + k
               do m = k, n
                  bx = dd(ipt)
                  p(k) = p(k) + bx*w(m)
                  if (k .ne. m)  p(m) = p(m) + bx*w(k)
                  ipt = ipt + 1
               end do
               mk = mk + mk1
               mk1 = mk1 - 1
            end do
            term = 0.0d0
            do k = ji, n
               term = term + w(k)*p(k)
            end do
            do k = ji, n
               p(k) = p(k) - term*w(k)
            end do
            mj = mi + mi1
            mj1 = mi1 - 1
            do j = ji, n
               do k = j, n
                  dd(mj+k) = dd(mj+k) - 2.0d0*(p(j)*w(k)+p(k)*w(j))
               end do
               mj = mj + mj1
               mj1 = mj1 - 1
            end do
         end if
         mi = mi + mi1
         mi1 = mi1 - 1
      end do
c
c     find the eigenvalues via the Sturm bisection method
c
      a(n) = dd(mi+n)
      b(n) = 0.0d0
      alimit = 1.0d0
      do i = 1, n
         w(i) = b(i)
         b(i) = b(i) * b(i)
      end do
      do i = 1, nv
         ev(i) = alimit
      end do
      root = -alimit
      do i = 1, nv
         rootx = alimit
         do j = i, nv
            rootx = min(rootx,ev(j))
         end do
         ev(i) = rootx
         trial = 0.5d0 * (root+ev(i))
         do while (abs(trial-root).ge.eps .and. abs(trial-ev(i)).ge.eps)
            nomtch = n
            j = 1
            do while (j .le. n)
               f0 = a(j) - trial
               do while (abs(f0) .ge. zeta)
                  if (f0 .ge. 0.0d0)  nomtch = nomtch - 1
                  j = j + 1
                  if (j .gt. n)  goto 10
                  f0 = a(j) - trial - b(j-1)/f0
               end do
               j = j + 2
               nomtch = nomtch - 1
            end do
   10       continue
            if (nomtch .lt. i) then
               root = trial
            else
               ev(i) = trial
               nom = min(nv,nomtch)
               ev(nom) = trial
            end if
            trial = 0.5d0 * (root+ev(i))
         end do
      end do
c
c     find the eigenvectors via a backtransformation step
c
      ia = -1
      do i = 1, nv
         root = ev(i)
         do j = 1, n
            y(j) = 1.0d0
         end do
         ia = ia + 1
         if (i .ne. 1) then
            if (abs(ev(i-1)-root) .ge. eps)  ia = 0
         end if
         elim1 = a(1) - root
         elim2 = w(1)
         do j = 1, nn
            if (abs(elim1) .le. abs(w(j))) then
               ta(j) = w(j)
               tb(j) = a(j+1) - root
               p(j) = w(j+1)
               temp = 1.0d0
               if (abs(w(j)) .gt. zeta)  temp = elim1 / w(j)
               elim1 = elim2 - temp*tb(j)
               elim2 = -temp * w(j+1)
            else
               ta(j) = elim1
               tb(j) = elim2
               p(j) = 0.0d0
               temp = w(j) / elim1
               elim1 = a(j+1) - root - temp*elim2
               elim2 = w(j+1)
            end if
            b(j) = temp
         end do
         ta(n) = elim1
         tb(n) = 0.0d0
         p(n) = 0.0d0
         p(nn) = 0.0d0
         iter = 1
         if (ia .ne. 0)  goto 40
   20    continue
         m = n + 1
         do j = 1, n
            m = m - 1
            done = .false.
            do while (.not. done)
               done = .true.
               k = n - m - 1
               if (k .lt. 0) then
                  elim1 = y(m)
               else if (k .eq. 0) then
                  elim1 = y(m) - y(m+1)*tb(m)
               else
                  elim1 = y(m) - y(m+1)*tb(m) - y(m+2)*p(m)
               end if
               if (abs(elim1) .le. sigma) then
                  temp = ta(m)
                  if (abs(temp) .lt. delta)  temp = delta
                  y(m) = elim1 / temp
               else
                  do k = 1, n
                     y(k) = y(k) / sigma
                  end do
                  done = .false.
               end if
            end do
         end do
         if (iter .eq. 2)  goto 50
         iter = iter + 1
   30    continue
         elim1 = y(1)
         do j = 1, nn
            if (ta(j) .ne. w(j)) then
               y(j) = elim1
               elim1 = y(j+1) - elim1*b(j)
            else
               y(j) = y(j+1)
               elim1 = elim1 - y(j+1)*b(j)
            end if
         end do
         y(n) = elim1
         goto 20
   40    continue
         do j = 1, n
            rand1 = mod(rvalue*rand1,rpower)
            y(j) = rand1/rpow1 - 1.0d0
         end do
         goto 20
   50    continue
         if (ia .ne. 0) then
            do j1 = 1, ia
               k = i - j1
               temp = 0.0d0
               do j = 1, n
                  temp = temp + y(j)*vec(j,k)
               end do
               do j = 1, n
                  y(j) = y(j) - temp*vec(j,k)
               end do
            end do
         end if
         if (iter .eq. 1)  goto 30
         elim1 = 0.0d0
         do j = 1, n
            elim1 = max(elim1,abs(y(j)))
         end do
         temp = 0.0d0
         do j = 1, n
            elim2 = y(j) / elim1
            temp = temp + elim2*elim2
         end do
         temp = 1.0d0 / (sqrt(temp)*elim1)
         do j = 1, n
            y(j) = y(j) * temp
            if (abs(y(j)) .lt. rho)  y(j) = 0.0d0
         end do
         do j = 1, n
            vec(j,i) = y(j)
         end do
      end do
c
c     normalization of the eigenvalues and eigenvectors
c
      do i = 1, nv
         do j = 1, n
            y(j) = vec(j,i)
         end do
         mk = (n*(n-1))/2 - 3
         mk1 = 3
         do j = 1, n-2
            t = 0.0d0
            m = n - j
            do k = m, n
               t = t + dd(mk+k)*y(k)
            end do
            do k = m, n
               epr = t * dd(mk+k)
               y(k) = y(k) - 2.0d0*epr
            end do
            mk = mk - mk1
            mk1 = mk1 + 1
         end do
         t = 0.0d0
         do j = 1, n
            t = t + y(j)*y(j)
         end do
         xnorm = sqrt(t)
         do j = 1, n
            y(j) = y(j) / xnorm
         end do
         do j = 1, n
            vec(j,i) = y(j)
         end do
      end do
      do i = 1, n
         ev(i) = ev(i) * anorm
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (a)
      deallocate (b)
      deallocate (p)
      deallocate (w)
      deallocate (ta)
      deallocate (tb)
      deallocate (y)
      return
      end
