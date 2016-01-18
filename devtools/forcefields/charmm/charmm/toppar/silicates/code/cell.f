      subroutine cell(alen,blen,clen,alpha,beta,gamma)
      implicit none
c
      double precision alpha,beta,gamma
      double precision alen,blen,clen
      integer iout
c
c     initialize periodic lengths and angles              
c
      iout = 6

      alen = 0.0d0
      blen = 0.0d0
      clen = 0.0d0
      alpha = 0.0d0
      beta = 0.0d0
      gamma = 0.0d0

      open (11,file='crystal.dat',status='old')
      read(11,*) alen,blen,clen,alpha,beta,gamma
      close (11)
c
c     checck unspecified periodic lengths and angles
c
      if (alen .eq. 0.0d0)  write(iout,'(/,a)') 'a lenght not set'
      if (blen .eq. 0.0d0)  write(iout,'(/,a)') 'b lenght not set'
      if (clen .eq. 0.0d0)  write(iout,'(/,a)') 'c lenght not set'
      if (alpha .eq. 0.0d0)  write(iout,'(/,a)') 'alpha not set'
      if (beta .eq. 0.0d0)  write(iout,'(/,a)') 'beta not set'
      if (gamma .eq. 0.0d0)  write(iout,'(/,a)') 'gamma not set'
c
c     stop if lengths or angles are nought               
c
      if (alen.eq.0.0d0 .or. blen.eq.0.0d0
     &       .or. clen.eq.90.0d0) stop
      if (alpha.eq.0.0d0 .or. beta.eq.0.0d0
     &       .or. gamma.eq.90.0d0) stop

      return
      end
