c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2010  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  module fft  --  Fast Fourier transform control values  ##
c     ##                                                         ##
c     #############################################################
c
c
c     maxprime   maximum number of prime factors of FFT dimension
c
c     iprime     prime factorization of each FFT dimension (fftpack)
c     planf      pointer to forward transform data structure (fftw)
c     planb      pointer to backward transform data structure (fftw)
c     ffttable   intermediate array used by the FFT routine (fftpack)
c     ffttyp     type of FFT package; currently FFTPACK or FFTW
c
c
      module fft
      implicit none
      integer maxprime
      parameter (maxprime=15)
      integer iprime(maxprime,3)
      integer*8 planf
      integer*8 planb
      real*8, allocatable :: ffttable(:,:)
      character*7 ffttyp
      save
      end
