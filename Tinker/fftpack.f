c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  routines below implement a 1-D Fast Fourier Transform;  ##
c     ##  code is modified from FFTPACK as obtained from Netlib;  ##
c     ##  original due to Paul N. Swarztrauber, NCAR, Boulder CO  ##
c     ##                                                          ##
c     ##############################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine cffti  --  1-D FFT setup and initialization  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "cffti" initializes arrays used in both forward and backward
c     transforms; "ifac" is the prime factorization of "n", and
c     "wsave" contains a tabulation of trigonometric functions
c
c
      subroutine cffti (n,wsave,ifac)
      implicit none
      integer n,iw
      integer ifac(*)
      real*8 wsave(*)
c
c
      if (n .gt. 1) then
         iw = n + n + 1
         call cffti1 (n,wsave(iw),ifac)
      end if
      return
      end
c
c
c     #########################
c     ##                     ##
c     ##  subroutine cffti1  ##
c     ##                     ##
c     #########################
c
c
      subroutine cffti1 (n,wa,ifac)
      use math
      implicit none
      integer i,j,ii,n,ip,ipm
      integer ib,ido,idot
      integer i1,k1,l1,l2,ld
      integer nl,nf,nq,nr
      integer ntry,ntryh(4)
      integer ifac(*)
      real*8 arg,argh,argld,fi
      real*8 wa(*)
      data ntryh  / 3, 4, 2, 5 /
c
c
      nl = n
      nf = 0
      j = 1
      ntry = ntryh(j)
      do while (nl .ne. 1)
         nq = nl / ntry
         nr = nl - ntry*nq
         if (nr .eq. 0) then
            nf = nf + 1
            ifac(nf+2) = ntry
            nl = nq
            if (ntry .eq. 2) then
               if (nf .ne. 1) then
                  do i = 2, nf
                     ib = nf - i + 2
                     ifac(ib+2) = ifac(ib+1)
                  end do
                  ifac(3) = 2
               end if
            end if
         else
            j = j + 1
            if (j .le. 4) then
               ntry = ntryh(j)
            else
               ntry = ntry + 2
            end if
         end if
      end do
      ifac(1) = n
      ifac(2) = nf
      argh = 2.0d0 * pi / dble(n)
      i = 2
      l1 = 1
      do k1 = 1, nf
         ip = ifac(k1+2)
         ld = 0
         l2 = l1 * ip
         ido = n / l2
         idot = ido + ido + 2
         ipm = ip - 1
         do j = 1, ipm
            i1 = i
            wa(i-1) = 1.0d0
            wa(i) = 0.0d0
            ld = ld + l1
            fi = 0.0d0
            argld = dble(ld) * argh
            do ii = 4, idot, 2
               i = i + 2
               fi = fi + 1.0d0
               arg = fi * argld
               wa(i-1) = cos(arg)
               wa(i) = sin(arg)
            end do
            if (ip .gt. 5) then
               wa(i1-1) = wa(i-1)
               wa(i1) = wa(i)
            end if
         end do
         l1 = l2
      end do
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  subroutine cfftf  --  1-D FFT forward transform  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "cfftf" computes the forward complex discrete Fourier
c     transform, the Fourier analysis
c
c
      subroutine cfftf (n,c,wsave,ifac)
      implicit none
      integer n,iw
      integer ifac(*)
      real*8 c(*)
      real*8 wsave(*)
c
c
      if (n .gt. 1) then
         iw = n + n + 1
         call cfftf1 (n,c,wsave,wsave(iw),ifac)
      end if
      return
      end
c
c
c     #########################
c     ##                     ##
c     ##  subroutine cfftf1  ##
c     ##                     ##
c     #########################
c
c
      subroutine cfftf1 (n,c,ch,wa,ifac)
      implicit none
      integer i,n,k1,l1,l2
      integer na,nac,nf,n2
      integer ido,idot,idl1,ip
      integer iw,ix2,ix3,ix4
      integer ifac(*)
      real*8 c(*),ch(*),wa(*)
c
c
      nf = ifac(2)
      na = 0
      l1 = 1
      iw = 1
      do k1 = 1, nf
         ip = ifac(k1+2)
         l2 = ip * l1
         ido = n / l2
         idot = ido + ido
         idl1 = idot * l1
         if (ip .eq. 5) then
            ix2 = iw + idot
            ix3 = ix2 + idot
            ix4 = ix3 + idot
            if (na .eq. 0) then
               call passf5 (idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
            else
               call passf5 (idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
            end if
            na = 1 - na
         else if (ip .eq. 4) then
            ix2 = iw + idot
            ix3 = ix2 + idot
            if (na .eq. 0) then
               call passf4 (idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
            else
               call passf4 (idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
            end if
            na = 1 - na
         else if (ip .eq. 3) then
            ix2 = iw + idot
            if (na .eq. 0) then
               call passf3 (idot,l1,c,ch,wa(iw),wa(ix2))
            else
               call passf3 (idot,l1,ch,c,wa(iw),wa(ix2))
            end if
            na = 1 - na
         else if (ip .eq. 2) then
            if (na .eq. 0) then
               call passf2 (idot,l1,c,ch,wa(iw))
            else
               call passf2 (idot,l1,ch,c,wa(iw))
            end if
            na = 1 - na
         else
            if (na .eq. 0) then
               call passf (nac,idot,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
            else
               call passf (nac,idot,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
            end if
            if (nac .ne. 0)  na = 1 - na
         end if
         l1 = l2
         iw = iw + (ip-1)*idot
      end do
      if (na .ne. 0) then
         n2 = n + n
         do i = 1, n2
            c(i) = ch(i)
         end do
      end if
      return
      end
c
c
c     ########################################################
c     ##                                                    ##
c     ##  subroutine cfftb  --  1-D FFT backward transform  ##
c     ##                                                    ##
c     ########################################################
c
c
c     "cfftb" computes the backward complex discrete Fourier
c     transform, the Fourier synthesis
c
c
      subroutine cfftb (n,c,wsave,ifac)
      implicit none
      integer n,iw
      integer ifac(*)
      real*8 c(*)
      real*8 wsave(*)
c
c
      if (n .gt. 1) then
         iw = n + n + 1
         call cfftb1 (n,c,wsave,wsave(iw),ifac)
      end if
      return
      end
c
c
c     #########################
c     ##                     ##
c     ##  subroutine cfftb1  ##
c     ##                     ##
c     #########################
c
c
      subroutine cfftb1 (n,c,ch,wa,ifac)
      implicit none
      integer i,n,k1,l1,l2
      integer na,nac,nf,n2
      integer ido,idot,idl1,ip
      integer iw,ix2,ix3,ix4
      integer ifac(*)
      real*8 c(*),ch(*),wa(*)
c
c
      nf = ifac(2)
      na = 0
      l1 = 1
      iw = 1
      do k1 = 1, nf
         ip = ifac(k1+2)
         l2 = ip * l1
         ido = n / l2
         idot = ido + ido
         idl1 = idot * l1
         if (ip .eq. 5) then
            ix2 = iw + idot
            ix3 = ix2 + idot
            ix4 = ix3 + idot
            if (na .eq. 0) then
               call passb5 (idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
            else
               call passb5 (idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
            end if
            na = 1 - na
         else if (ip .eq. 4) then
            ix2 = iw + idot
            ix3 = ix2 + idot
            if (na .eq. 0) then
               call passb4 (idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
            else
               call passb4 (idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
            end if
            na = 1 - na
         else if (ip .eq. 3) then
            ix2 = iw + idot
            if (na .eq. 0) then
               call passb3 (idot,l1,c,ch,wa(iw),wa(ix2))
            else
               call passb3 (idot,l1,ch,c,wa(iw),wa(ix2))
            end if
            na = 1 - na
         else if (ip .eq. 2) then
            if (na .eq. 0) then
               call passb2 (idot,l1,c,ch,wa(iw))
            else
               call passb2 (idot,l1,ch,c,wa(iw))
            end if
            na = 1 - na
         else
            if (na .eq. 0) then
               call passb (nac,idot,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
            else
               call passb (nac,idot,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
            end if
            if (nac .ne. 0)  na = 1 - na
         end if
         l1 = l2
         iw = iw + (ip-1)*idot
      end do
      if (na .ne. 0) then
         n2 = n + n
         do i = 1, n2
            c(i) = ch(i)
         end do
      end if
      return
      end
c
c
c     ########################
c     ##                    ##
c     ##  subroutine passf  ##
c     ##                    ##
c     ########################
c
c
      subroutine passf (nac,ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
      implicit none
      integer nac,ido,ip
      integer l1,idl1
      integer i,j,k,l
      integer ik,jc,lc
      integer idj,idl,idp
      integer idij,idlj
      integer inc,idot,nt
      integer ipp2,ipph
      real*8 wai,war
      real*8 cc(ido,ip,l1)
      real*8 c1(ido,l1,ip)
      real*8 c2(idl1,ip)
      real*8 ch(ido,l1,ip)
      real*8 ch2(idl1,ip)
      real*8 wa(*)
c
c
      idot = ido / 2
      nt = ip * idl1
      ipp2 = ip + 2
      ipph = (ip+1) / 2
      idp = ip * ido
      if (ido .ge. l1) then
         do j = 2, ipph
            jc = ipp2 - j
            do k = 1, l1
               do i = 1, ido
                  ch(i,k,j) = cc(i,j,k) + cc(i,jc,k)
                  ch(i,k,jc) = cc(i,j,k) - cc(i,jc,k)
               end do
            end do
         end do
         do k = 1, l1
            do i = 1, ido
               ch(i,k,1) = cc(i,1,k)
            end do
         end do
      else
         do j = 2, ipph
            jc = ipp2 - j
            do i = 1, ido
               do k = 1, l1
                  ch(i,k,j) = cc(i,j,k) + cc(i,jc,k)
                  ch(i,k,jc) = cc(i,j,k) - cc(i,jc,k)
               end do
            end do
         end do
         do i = 1, ido
            do k = 1, l1
               ch(i,k,1) = cc(i,1,k)
            end do
         end do
      end if
      idl = 2 - ido
      inc = 0
      do l = 2, ipph
         lc = ipp2 - l
         idl = idl + ido
         do ik = 1, idl1
            c2(ik,l) = ch2(ik,1) + wa(idl-1)*ch2(ik,2)
            c2(ik,lc) = -wa(idl) * ch2(ik,ip)
         end do
         idlj = idl
         inc = inc + ido
         do j = 3, ipph
            jc = ipp2 - j
            idlj = idlj + inc
            if (idlj .gt. idp)  idlj = idlj - idp
            war = wa(idlj-1)
            wai = wa(idlj)
            do ik = 1, idl1
               c2(ik,l) = c2(ik,l) + war*ch2(ik,j)
               c2(ik,lc) = c2(ik,lc) - wai*ch2(ik,jc)
            end do
         end do
      end do
      do j = 2, ipph
         do ik = 1, idl1
            ch2(ik,1) = ch2(ik,1) + ch2(ik,j)
         end do
      end do
      do j = 2, ipph
         jc = ipp2 - j
         do ik = 2, idl1, 2
            ch2(ik-1,j) = c2(ik-1,j) - c2(ik,jc)
            ch2(ik-1,jc) = c2(ik-1,j) + c2(ik,jc)
            ch2(ik,j) = c2(ik,j) + c2(ik-1,jc)
            ch2(ik,jc) = c2(ik,j) - c2(ik-1,jc)
         end do
      end do
      nac = 1
      if (ido .ne. 2) then
         nac = 0
         do ik = 1, idl1
            c2(ik,1) = ch2(ik,1)
         end do
         do j = 2, ip
            do k = 1, l1
               c1(1,k,j) = ch(1,k,j)
               c1(2,k,j) = ch(2,k,j)
            end do
         end do
         if (idot .le. l1) then
            idij = 0
            do j = 2, ip
               idij = idij + 2
               do i = 4, ido, 2
                  idij = idij + 2
                  do k = 1, l1
                     c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)
     &                                + wa(idij)*ch(i,k,j)
                     c1(i,k,j) = wa(idij-1)*ch(i,k,j)
     &                              - wa(idij)*ch(i-1,k,j)
                  end do
               end do
            end do
         else
            idj = 2 - ido
            do j = 2, ip
               idj = idj + ido
               do k = 1, l1
                  idij = idj
                  do i = 4, ido, 2
                     idij = idij + 2
                     c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)
     &                                + wa(idij)*ch(i,k,j)
                     c1(i,k,j) = wa(idij-1)*ch(i,k,j)
     &                              - wa(idij)*ch(i-1,k,j)
                  end do
               end do
            end do
         end if
      end if
      return
      end
c
c
c     #########################
c     ##                     ##
c     ##  subroutine passf2  ##
c     ##                     ##
c     #########################
c
c
      subroutine passf2 (ido,l1,cc,ch,wa1)
      implicit none
      integer i,k,ido,l1
      real*8 ti2,tr2
      real*8 cc(ido,2,l1)
      real*8 ch(ido,l1,2)
      real*8 wa1(*)
c
c
      if (ido .le. 2) then
         do k = 1, l1
            ch(1,k,1) = cc(1,1,k) + cc(1,2,k)
            ch(1,k,2) = cc(1,1,k) - cc(1,2,k)
            ch(2,k,1) = cc(2,1,k) + cc(2,2,k)
            ch(2,k,2) = cc(2,1,k) - cc(2,2,k)
         end do
      else
         do k = 1, l1
            do i = 2, ido, 2
               ch(i-1,k,1) = cc(i-1,1,k) + cc(i-1,2,k)
               tr2 = cc(i-1,1,k) - cc(i-1,2,k)
               ch(i,k,1) = cc(i,1,k) + cc(i,2,k)
               ti2 = cc(i,1,k) - cc(i,2,k)
               ch(i,k,2) = wa1(i-1)*ti2 - wa1(i)*tr2
               ch(i-1,k,2) = wa1(i-1)*tr2 + wa1(i)*ti2
            end do
         end do
      end if
      return
      end
c
c
c     #########################
c     ##                     ##
c     ##  subroutine passf3  ##
c     ##                     ##
c     #########################
c
c
      subroutine passf3 (ido,l1,cc,ch,wa1,wa2)
      implicit none
      integer i,k,ido,l1
      real*8 ti2,tr2,taur,taui
      real*8 ci2,ci3,cr2,cr3
      real*8 di2,di3,dr2,dr3
      real*8 cc(ido,3,l1)
      real*8 ch(ido,l1,3)
      real*8 wa1(*)
      real*8 wa2(*)
      data taur  / -0.5d0 /
      data taui  / -0.866025403784439d0 /
c
c
      if (ido .eq. 2) then
         do k = 1, l1
            tr2 = cc(1,2,k) + cc(1,3,k)
            cr2 = cc(1,1,k) + taur*tr2
            ch(1,k,1) = cc(1,1,k) + tr2
            ti2 = cc(2,2,k) + cc(2,3,k)
            ci2 = cc(2,1,k) + taur*ti2
            ch(2,k,1) = cc(2,1,k) + ti2
            cr3 = taui * (cc(1,2,k)-cc(1,3,k))
            ci3 = taui * (cc(2,2,k)-cc(2,3,k))
            ch(1,k,2) = cr2 - ci3
            ch(1,k,3) = cr2 + ci3
            ch(2,k,2) = ci2 + cr3
            ch(2,k,3) = ci2 - cr3
         end do
      else
         do k = 1, l1
            do i = 2, ido, 2
               tr2 = cc(i-1,2,k) + cc(i-1,3,k)
               cr2 = cc(i-1,1,k) + taur*tr2
               ch(i-1,k,1) = cc(i-1,1,k) + tr2
               ti2 = cc(i,2,k) + cc(i,3,k)
               ci2 = cc(i,1,k) + taur*ti2
               ch(i,k,1) = cc(i,1,k) + ti2
               cr3 = taui * (cc(i-1,2,k)-cc(i-1,3,k))
               ci3 = taui * (cc(i,2,k)-cc(i,3,k))
               dr2 = cr2 - ci3
               dr3 = cr2 + ci3
               di2 = ci2 + cr3
               di3 = ci2 - cr3
               ch(i,k,2) = wa1(i-1)*di2 - wa1(i)*dr2
               ch(i-1,k,2) = wa1(i-1)*dr2 + wa1(i)*di2
               ch(i,k,3) = wa2(i-1)*di3 - wa2(i)*dr3
               ch(i-1,k,3) = wa2(i-1)*dr3 + wa2(i)*di3
            end do
         end do
      end if
      return
      end
c
c
c     #########################
c     ##                     ##
c     ##  subroutine passf4  ##
c     ##                     ##
c     #########################
c
c
      subroutine passf4 (ido,l1,cc,ch,wa1,wa2,wa3)
      implicit none
      integer i,k,ido,l1
      real*8 ci2,ci3,ci4
      real*8 cr2,cr3,cr4
      real*8 ti1,ti2,ti3,ti4
      real*8 tr1,tr2,tr3,tr4
      real*8 cc(ido,4,l1)
      real*8 ch(ido,l1,4)
      real*8 wa1(*)
      real*8 wa2(*)
      real*8 wa3(*)
c
c
      if (ido .eq. 2) then
         do k = 1, l1
            ti1 = cc(2,1,k) - cc(2,3,k)
            ti2 = cc(2,1,k) + cc(2,3,k)
            tr4 = cc(2,2,k) - cc(2,4,k)
            ti3 = cc(2,2,k) + cc(2,4,k)
            tr1 = cc(1,1,k) - cc(1,3,k)
            tr2 = cc(1,1,k) + cc(1,3,k)
            ti4 = cc(1,4,k) - cc(1,2,k)
            tr3 = cc(1,2,k) + cc(1,4,k)
            ch(1,k,1) = tr2 + tr3
            ch(1,k,3) = tr2 - tr3
            ch(2,k,1) = ti2 + ti3
            ch(2,k,3) = ti2 - ti3
            ch(1,k,2) = tr1 + tr4
            ch(1,k,4) = tr1 - tr4
            ch(2,k,2) = ti1 + ti4
            ch(2,k,4) = ti1 - ti4
         end do
      else
         do k = 1, l1
            do i = 2, ido, 2
               ti1 = cc(i,1,k) - cc(i,3,k)
               ti2 = cc(i,1,k) + cc(i,3,k)
               ti3 = cc(i,2,k) + cc(i,4,k)
               tr4 = cc(i,2,k) - cc(i,4,k)
               tr1 = cc(i-1,1,k) - cc(i-1,3,k)
               tr2 = cc(i-1,1,k) + cc(i-1,3,k)
               ti4 = cc(i-1,4,k) - cc(i-1,2,k)
               tr3 = cc(i-1,2,k) + cc(i-1,4,k)
               ch(i-1,k,1) = tr2 + tr3
               cr3 = tr2 - tr3
               ch(i,k,1) = ti2 + ti3
               ci3 = ti2 - ti3
               cr2 = tr1 + tr4
               cr4 = tr1 - tr4
               ci2 = ti1 + ti4
               ci4 = ti1 - ti4
               ch(i-1,k,2) = wa1(i-1)*cr2 + wa1(i)*ci2
               ch(i,k,2) = wa1(i-1)*ci2 - wa1(i)*cr2
               ch(i-1,k,3) = wa2(i-1)*cr3 + wa2(i)*ci3
               ch(i,k,3) = wa2(i-1)*ci3 - wa2(i)*cr3
               ch(i-1,k,4) = wa3(i-1)*cr4 + wa3(i)*ci4
               ch(i,k,4) = wa3(i-1)*ci4 - wa3(i)*cr4
            end do
         end do
      end if
      return
      end
c
c
c     #########################
c     ##                     ##
c     ##  subroutine passf5  ##
c     ##                     ##
c     #########################
c
c
      subroutine passf5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)
      implicit none
      integer i,k,ido,l1
      real*8 ci2,ci3,ci4,ci5
      real*8 cr2,cr3,cr4,cr5
      real*8 di2,di3,di4,di5
      real*8 dr2,dr3,dr4,dr5
      real*8 ti2,ti3,ti4,ti5
      real*8 tr2,tr3,tr4,tr5
      real*8 tr11,ti11
      real*8 tr12,ti12
      real*8 cc(ido,5,l1)
      real*8 ch(ido,l1,5)
      real*8 wa1(*)
      real*8 wa2(*)
      real*8 wa3(*)
      real*8 wa4(*)
      data tr11  /  0.309016994374947d0 /
      data ti11  / -0.951056516295154d0 /
      data tr12  / -0.809016994374947d0 /
      data ti12  / -0.587785252292473d0 /
c
c
      if (ido .eq. 2) then
         do k = 1, l1
            ti5 = cc(2,2,k) - cc(2,5,k)
            ti2 = cc(2,2,k) + cc(2,5,k)
            ti4 = cc(2,3,k) - cc(2,4,k)
            ti3 = cc(2,3,k) + cc(2,4,k)
            tr5 = cc(1,2,k) - cc(1,5,k)
            tr2 = cc(1,2,k) + cc(1,5,k)
            tr4 = cc(1,3,k) - cc(1,4,k)
            tr3 = cc(1,3,k) + cc(1,4,k)
            ch(1,k,1) = cc(1,1,k) + tr2 + tr3
            ch(2,k,1) = cc(2,1,k) + ti2 + ti3
            cr2 = cc(1,1,k) + tr11*tr2 + tr12*tr3
            ci2 = cc(2,1,k) + tr11*ti2 + tr12*ti3
            cr3 = cc(1,1,k) + tr12*tr2 + tr11*tr3
            ci3 = cc(2,1,k) + tr12*ti2 + tr11*ti3
            cr5 = ti11*tr5 + ti12*tr4
            ci5 = ti11*ti5 + ti12*ti4
            cr4 = ti12*tr5 - ti11*tr4
            ci4 = ti12*ti5 - ti11*ti4
            ch(1,k,2) = cr2 - ci5
            ch(1,k,5) = cr2 + ci5
            ch(2,k,2) = ci2 + cr5
            ch(2,k,3) = ci3 + cr4
            ch(1,k,3) = cr3 - ci4
            ch(1,k,4) = cr3 + ci4
            ch(2,k,4) = ci3 - cr4
            ch(2,k,5) = ci2 - cr5
         end do
      else
         do k = 1, l1
            do i = 2, ido, 2
               ti5 = cc(i,2,k) - cc(i,5,k)
               ti2 = cc(i,2,k) + cc(i,5,k)
               ti4 = cc(i,3,k) - cc(i,4,k)
               ti3 = cc(i,3,k) + cc(i,4,k)
               tr5 = cc(i-1,2,k) - cc(i-1,5,k)
               tr2 = cc(i-1,2,k) + cc(i-1,5,k)
               tr4 = cc(i-1,3,k) - cc(i-1,4,k)
               tr3 = cc(i-1,3,k) + cc(i-1,4,k)
               ch(i-1,k,1) = cc(i-1,1,k) + tr2 + tr3
               ch(i,k,1) = cc(i,1,k) + ti2 + ti3
               cr2 = cc(i-1,1,k) + tr11*tr2 + tr12*tr3
               ci2 = cc(i,1,k) + tr11*ti2 + tr12*ti3
               cr3 = cc(i-1,1,k) + tr12*tr2 + tr11*tr3
               ci3 = cc(i,1,k) + tr12*ti2 + tr11*ti3
               cr5 = ti11*tr5 + ti12*tr4
               ci5 = ti11*ti5 + ti12*ti4
               cr4 = ti12*tr5 - ti11*tr4
               ci4 = ti12*ti5 - ti11*ti4
               dr3 = cr3 - ci4
               dr4 = cr3 + ci4
               di3 = ci3 + cr4
               di4 = ci3 - cr4
               dr5 = cr2 + ci5
               dr2 = cr2 - ci5
               di5 = ci2 - cr5
               di2 = ci2 + cr5
               ch(i-1,k,2) = wa1(i-1)*dr2 + wa1(i)*di2
               ch(i,k,2) = wa1(i-1)*di2 - wa1(i)*dr2
               ch(i-1,k,3) = wa2(i-1)*dr3 + wa2(i)*di3
               ch(i,k,3) = wa2(i-1)*di3 - wa2(i)*dr3
               ch(i-1,k,4) = wa3(i-1)*dr4 + wa3(i)*di4
               ch(i,k,4) = wa3(i-1)*di4 - wa3(i)*dr4
               ch(i-1,k,5) = wa4(i-1)*dr5 + wa4(i)*di5
               ch(i,k,5) = wa4(i-1)*di5 - wa4(i)*dr5
            end do
         end do
      end if
      return
      end
c
c
c     ########################
c     ##                    ##
c     ##  subroutine passb  ##
c     ##                    ##
c     ########################
c
c
      subroutine passb (nac,ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
      implicit none
      integer nac,ido,ip,l1,idl1
      integer i,j,k,l,ik,jc,lc
      integer idj,idl,idp,idij,idlj
      integer inc,idot,nt,ipp2,ipph
      real*8 wai,war
      real*8 cc(ido,ip,l1)
      real*8 c1(ido,l1,ip)
      real*8 c2(idl1,ip)
      real*8 ch(ido,l1,ip)
      real*8 ch2(idl1,ip)
      real*8 wa(*)
c
c
      idot = ido / 2
      nt = ip * idl1
      ipp2 = ip + 2
      ipph = (ip+1) / 2
      idp = ip * ido
      if (ido .ge. l1) then
         do j = 2, ipph
            jc = ipp2 - j
            do k = 1, l1
               do i = 1, ido
                  ch(i,k,j) = cc(i,j,k) + cc(i,jc,k)
                  ch(i,k,jc) = cc(i,j,k) - cc(i,jc,k)
               end do
            end do
         end do
         do k = 1, l1
            do i = 1, ido
               ch(i,k,1) = cc(i,1,k)
            end do
         end do
      else
         do j = 2, ipph
            jc = ipp2 - j
            do i = 1, ido
               do k = 1, l1
                  ch(i,k,j) = cc(i,j,k) + cc(i,jc,k)
                  ch(i,k,jc) = cc(i,j,k) - cc(i,jc,k)
               end do
            end do
         end do
         do i = 1, ido
            do k = 1, l1
               ch(i,k,1) = cc(i,1,k)
            end do
         end do
      end if
      idl = 2 - ido
      inc = 0
      do l = 2, ipph
         lc = ipp2 - l
         idl = idl + ido
         do ik = 1, idl1
            c2(ik,l) = ch2(ik,1) + wa(idl-1)*ch2(ik,2)
            c2(ik,lc) = wa(idl) * ch2(ik,ip)
         end do
         idlj = idl
         inc = inc + ido
         do j = 3, ipph
            jc = ipp2 - j
            idlj = idlj + inc
            if (idlj .gt. idp)  idlj = idlj - idp
            war = wa(idlj-1)
            wai = wa(idlj)
            do ik = 1, idl1
               c2(ik,l) = c2(ik,l) + war*ch2(ik,j)
               c2(ik,lc) = c2(ik,lc) + wai*ch2(ik,jc)
            end do
         end do
      end do
      do j = 2, ipph
         do ik = 1, idl1
            ch2(ik,1) = ch2(ik,1) + ch2(ik,j)
         end do
      end do
      do j = 2, ipph
         jc = ipp2 - j
         do ik = 2, idl1, 2
            ch2(ik-1,j) = c2(ik-1,j) - c2(ik,jc)
            ch2(ik-1,jc) = c2(ik-1,j) + c2(ik,jc)
            ch2(ik,j) = c2(ik,j) + c2(ik-1,jc)
            ch2(ik,jc) = c2(ik,j) - c2(ik-1,jc)
         end do
      end do
      nac = 1
      if (ido .ne. 2) then
         nac = 0
         do ik = 1, idl1
            c2(ik,1) = ch2(ik,1)
         end do
         do j = 2, ip
            do k = 1, l1
               c1(1,k,j) = ch(1,k,j)
               c1(2,k,j) = ch(2,k,j)
            end do
         end do
         if (idot .le. l1) then
            idij = 0
            do j = 2, ip
               idij = idij + 2
               do i = 4, ido, 2
                  idij = idij + 2
                  do k = 1, l1
                     c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)
     &                                - wa(idij)*ch(i,k,j)
                     c1(i,k,j) = wa(idij-1)*ch(i,k,j)
     &                              + wa(idij)*ch(i-1,k,j)
                  end do
               end do
            end do
         else
            idj = 2 - ido
            do j = 2, ip
               idj = idj + ido
               do k = 1, l1
                  idij = idj
                  do i = 4, ido, 2
                     idij = idij + 2
                     c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)
     &                                - wa(idij)*ch(i,k,j)
                     c1(i,k,j) = wa(idij-1)*ch(i,k,j)
     &                              + wa(idij)*ch(i-1,k,j)
                  end do
               end do
            end do
         end if
      end if
      return
      end
c
c
c     #########################
c     ##                     ##
c     ##  subroutine passb2  ##
c     ##                     ##
c     #########################
c
c
      subroutine passb2 (ido,l1,cc,ch,wa1)
      implicit none
      integer i,k,ido,l1
      real*8 ti2,tr2
      real*8 cc(ido,2,l1)
      real*8 ch(ido,l1,2)
      real*8 wa1(*)
c
c
      if (ido .le. 2) then
         do k = 1, l1
            ch(1,k,1) = cc(1,1,k) + cc(1,2,k)
            ch(1,k,2) = cc(1,1,k) - cc(1,2,k)
            ch(2,k,1) = cc(2,1,k) + cc(2,2,k)
            ch(2,k,2) = cc(2,1,k) - cc(2,2,k)
         end do
      else
         do k = 1, l1
            do i = 2, ido, 2
               ch(i-1,k,1) = cc(i-1,1,k) + cc(i-1,2,k)
               tr2 = cc(i-1,1,k) - cc(i-1,2,k)
               ch(i,k,1) = cc(i,1,k) + cc(i,2,k)
               ti2 = cc(i,1,k) - cc(i,2,k)
               ch(i,k,2) = wa1(i-1)*ti2 + wa1(i)*tr2
               ch(i-1,k,2) = wa1(i-1)*tr2 - wa1(i)*ti2
            end do
         end do
      end if
      return
      end
c
c
c     #########################
c     ##                     ##
c     ##  subroutine passb3  ##
c     ##                     ##
c     #########################
c
c
      subroutine passb3 (ido,l1,cc,ch,wa1,wa2)
      implicit none
      integer i,k,ido,l1
      real*8 ti2,tr2,taur,taui
      real*8 ci2,ci3,cr2,cr3
      real*8 di2,di3,dr2,dr3
      real*8 cc(ido,3,l1)
      real*8 ch(ido,l1,3)
      real*8 wa1(*)
      real*8 wa2(*)
      data taur  / -0.5d0 /
      data taui  / 0.866025403784439d0 /
c
c
      if (ido .eq. 2) then
         do k = 1, l1
            tr2 = cc(1,2,k) + cc(1,3,k)
            cr2 = cc(1,1,k) + taur*tr2
            ch(1,k,1) = cc(1,1,k) + tr2
            ti2 = cc(2,2,k) + cc(2,3,k)
            ci2 = cc(2,1,k) + taur*ti2
            ch(2,k,1) = cc(2,1,k) + ti2
            cr3 = taui * (cc(1,2,k)-cc(1,3,k))
            ci3 = taui * (cc(2,2,k)-cc(2,3,k))
            ch(1,k,2) = cr2 - ci3
            ch(1,k,3) = cr2 + ci3
            ch(2,k,2) = ci2 + cr3
            ch(2,k,3) = ci2 - cr3
         end do
      else
         do k = 1, l1
            do i = 2, ido, 2
               tr2 = cc(i-1,2,k) + cc(i-1,3,k)
               cr2 = cc(i-1,1,k) + taur*tr2
               ch(i-1,k,1) = cc(i-1,1,k) + tr2
               ti2 = cc(i,2,k) + cc(i,3,k)
               ci2 = cc(i,1,k) + taur*ti2
               ch(i,k,1) = cc(i,1,k) + ti2
               cr3 = taui * (cc(i-1,2,k)-cc(i-1,3,k))
               ci3 = taui * (cc(i,2,k)-cc(i,3,k))
               dr2 = cr2 - ci3
               dr3 = cr2 + ci3
               di2 = ci2 + cr3
               di3 = ci2 - cr3
               ch(i,k,2) = wa1(i-1)*di2 + wa1(i)*dr2
               ch(i-1,k,2) = wa1(i-1)*dr2 - wa1(i)*di2
               ch(i,k,3) = wa2(i-1)*di3 + wa2(i)*dr3
               ch(i-1,k,3) = wa2(i-1)*dr3 - wa2(i)*di3
            end do
         end do
      end if
      return
      end
c
c
c     #########################
c     ##                     ##
c     ##  subroutine passb4  ##
c     ##                     ##
c     #########################
c
c
      subroutine passb4 (ido,l1,cc,ch,wa1,wa2,wa3)
      implicit none
      integer i,k,ido,l1
      real*8 ci2,ci3,ci4
      real*8 cr2,cr3,cr4
      real*8 ti1,ti2,ti3,ti4
      real*8 tr1,tr2,tr3,tr4
      real*8 cc(ido,4,l1)
      real*8 ch(ido,l1,4)
      real*8 wa1(*)
      real*8 wa2(*)
      real*8 wa3(*)
c
c
      if (ido .eq. 2) then
         do k = 1, l1
            ti1 = cc(2,1,k) - cc(2,3,k)
            ti2 = cc(2,1,k) + cc(2,3,k)
            tr4 = cc(2,4,k) - cc(2,2,k)
            ti3 = cc(2,2,k) + cc(2,4,k)
            tr1 = cc(1,1,k) - cc(1,3,k)
            tr2 = cc(1,1,k) + cc(1,3,k)
            ti4 = cc(1,2,k) - cc(1,4,k)
            tr3 = cc(1,2,k) + cc(1,4,k)
            ch(1,k,1) = tr2 + tr3
            ch(1,k,3) = tr2 - tr3
            ch(2,k,1) = ti2 + ti3
            ch(2,k,3) = ti2 - ti3
            ch(1,k,2) = tr1 + tr4
            ch(1,k,4) = tr1 - tr4
            ch(2,k,2) = ti1 + ti4
            ch(2,k,4) = ti1 - ti4
         end do
      else
         do k = 1, l1
            do i = 2, ido, 2
               ti1 = cc(i,1,k) - cc(i,3,k)
               ti2 = cc(i,1,k) + cc(i,3,k)
               ti3 = cc(i,2,k) + cc(i,4,k)
               tr4 = cc(i,4,k) - cc(i,2,k)
               tr1 = cc(i-1,1,k) - cc(i-1,3,k)
               tr2 = cc(i-1,1,k) + cc(i-1,3,k)
               ti4 = cc(i-1,2,k) - cc(i-1,4,k)
               tr3 = cc(i-1,2,k) + cc(i-1,4,k)
               ch(i-1,k,1) = tr2 + tr3
               cr3 = tr2 - tr3
               ch(i,k,1) = ti2 + ti3
               ci3 = ti2 - ti3
               cr2 = tr1 + tr4
               cr4 = tr1 - tr4
               ci2 = ti1 + ti4
               ci4 = ti1 - ti4
               ch(i-1,k,2) = wa1(i-1)*cr2 - wa1(i)*ci2
               ch(i,k,2) = wa1(i-1)*ci2 + wa1(i)*cr2
               ch(i-1,k,3) = wa2(i-1)*cr3 - wa2(i)*ci3
               ch(i,k,3) = wa2(i-1)*ci3 + wa2(i)*cr3
               ch(i-1,k,4) = wa3(i-1)*cr4 - wa3(i)*ci4
               ch(i,k,4) = wa3(i-1)*ci4 + wa3(i)*cr4
            end do
         end do
      end if
      return
      end
c
c
c     #########################
c     ##                     ##
c     ##  subroutine passb5  ##
c     ##                     ##
c     #########################
c
c
      subroutine passb5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)
      implicit none
      integer i,k,ido,l1
      real*8 ci2,ci3,ci4,ci5
      real*8 cr2,cr3,cr4,cr5
      real*8 di2,di3,di4,di5
      real*8 dr2,dr3,dr4,dr5
      real*8 ti2,ti3,ti4,ti5
      real*8 tr2,tr3,tr4,tr5
      real*8 tr11,ti11
      real*8 tr12,ti12
      real*8 cc(ido,5,l1)
      real*8 ch(ido,l1,5)
      real*8 wa1(*)
      real*8 wa2(*)
      real*8 wa3(*)
      real*8 wa4(*)
      data tr11  /  0.309016994374947d0 /
      data ti11  /  0.951056516295154d0 /
      data tr12  / -0.809016994374947d0 /
      data ti12  /  0.587785252292473d0 /
c
c
      if (ido .eq. 2) then
         do k = 1, l1
            ti5 = cc(2,2,k) - cc(2,5,k)
            ti2 = cc(2,2,k) + cc(2,5,k)
            ti4 = cc(2,3,k) - cc(2,4,k)
            ti3 = cc(2,3,k) + cc(2,4,k)
            tr5 = cc(1,2,k) - cc(1,5,k)
            tr2 = cc(1,2,k) + cc(1,5,k)
            tr4 = cc(1,3,k) - cc(1,4,k)
            tr3 = cc(1,3,k) + cc(1,4,k)
            ch(1,k,1) = cc(1,1,k) + tr2 + tr3
            ch(2,k,1) = cc(2,1,k) + ti2 + ti3
            cr2 = cc(1,1,k) + tr11*tr2 + tr12*tr3
            ci2 = cc(2,1,k) + tr11*ti2 + tr12*ti3
            cr3 = cc(1,1,k) + tr12*tr2 + tr11*tr3
            ci3 = cc(2,1,k) + tr12*ti2 + tr11*ti3
            cr5 = ti11*tr5 + ti12*tr4
            ci5 = ti11*ti5 + ti12*ti4
            cr4 = ti12*tr5 - ti11*tr4
            ci4 = ti12*ti5 - ti11*ti4
            ch(1,k,2) = cr2 - ci5
            ch(1,k,5) = cr2 + ci5
            ch(2,k,2) = ci2 + cr5
            ch(2,k,3) = ci3 + cr4
            ch(1,k,3) = cr3 - ci4
            ch(1,k,4) = cr3 + ci4
            ch(2,k,4) = ci3 - cr4
            ch(2,k,5) = ci2 - cr5
         end do
      else
         do k = 1, l1
            do i = 2, ido, 2
               ti5 = cc(i,2,k) - cc(i,5,k)
               ti2 = cc(i,2,k) + cc(i,5,k)
               ti4 = cc(i,3,k) - cc(i,4,k)
               ti3 = cc(i,3,k) + cc(i,4,k)
               tr5 = cc(i-1,2,k) - cc(i-1,5,k)
               tr2 = cc(i-1,2,k) + cc(i-1,5,k)
               tr4 = cc(i-1,3,k) - cc(i-1,4,k)
               tr3 = cc(i-1,3,k) + cc(i-1,4,k)
               ch(i-1,k,1) = cc(i-1,1,k) + tr2 + tr3
               ch(i,k,1) = cc(i,1,k) + ti2 + ti3
               cr2 = cc(i-1,1,k) + tr11*tr2 + tr12*tr3
               ci2 = cc(i,1,k) + tr11*ti2 + tr12*ti3
               cr3 = cc(i-1,1,k) + tr12*tr2 + tr11*tr3
               ci3 = cc(i,1,k) + tr12*ti2 + tr11*ti3
               cr5 = ti11*tr5 + ti12*tr4
               ci5 = ti11*ti5 + ti12*ti4
               cr4 = ti12*tr5 - ti11*tr4
               ci4 = ti12*ti5 - ti11*ti4
               dr3 = cr3 - ci4
               dr4 = cr3 + ci4
               di3 = ci3 + cr4
               di4 = ci3 - cr4
               dr5 = cr2 + ci5
               dr2 = cr2 - ci5
               di5 = ci2 - cr5
               di2 = ci2 + cr5
               ch(i-1,k,2) = wa1(i-1)*dr2 - wa1(i)*di2
               ch(i,k,2) = wa1(i-1)*di2 + wa1(i)*dr2
               ch(i-1,k,3) = wa2(i-1)*dr3 - wa2(i)*di3
               ch(i,k,3) = wa2(i-1)*di3 + wa2(i)*dr3
               ch(i-1,k,4) = wa3(i-1)*dr4 - wa3(i)*di4
               ch(i,k,4) = wa3(i-1)*di4 + wa3(i)*dr4
               ch(i-1,k,5) = wa4(i-1)*dr5 - wa4(i)*di5
               ch(i,k,5) = wa4(i-1)*di5 + wa4(i)*dr5
            end do
         end do
      end if
      return
      end
