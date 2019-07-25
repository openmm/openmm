c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine fftsetup  --  setup 3-D Fast Fourier transform  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "fftsetup" does initialization for a 3-D FFT via a single 3-D
c     transform (FFTW) or three separate 1-D transforms (FFTPACK)
c
c
      subroutine fftsetup
      use fft
      use openmp
      use pme
      implicit none
      integer maxtable
!$    integer ifront,iback
!$    integer error,iguess
c
c
c     initialization of Fast Fourier transform using FFTW;
c     comment "dfftw_init_threads" and "dfftw_plan_with_threads"
c     calls if serial FFTW is used in place of OpenMP FFTW
c
!$    if (ffttyp .eq. 'FFTW') then
!$       ifront = -1
!$       iback = 1
!$       error = 0
!$       iguess = 0
!$       call dfftw_init_threads (error)
!$       call dfftw_plan_with_nthreads (nthread)
!$       call dfftw_plan_dft_3d (planf,nfft1,nfft2,nfft3,qgrid,
!$   &                              qgrid,ifront,iguess)
!$       call dfftw_plan_dft_3d (planb,nfft1,nfft2,nfft3,qgrid,
!$   &                              qgrid,iback,iguess)
!$    else
c
c     perform dynamic allocation of some global arrays
c
         maxtable = 4 * max(nfft1,nfft2,nfft3)
         if (.not. allocated(ffttable))  allocate (ffttable(maxtable,3))
c
c     initialization of Fast Fourier transform using FFTPACK
c
         call cffti (nfft1,ffttable(1,1),iprime(1,1))
         call cffti (nfft2,ffttable(1,2),iprime(1,2))
         call cffti (nfft3,ffttable(1,3),iprime(1,3))
!$    end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine fftfront  --  forward Fast Fourier transform  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "fftfront" performs a 3-D FFT forward transform via a single
c     3-D transform or three separate 1-D transforms
c
c
      subroutine fftfront
      use fft
      use pme
      implicit none
      integer i,j,k
      real*8, allocatable :: work(:,:)
c
c
c     perform a single 3-D forward transform using FFTW
c
!$    if (ffttyp .eq. 'FFTW') then
!$       call dfftw_execute_dft (planf,qgrid,qgrid)
!$    else
c
c     perform three 1-D forward transforms using FFTPACK
c
         allocate (work(2,max(nfft1,nfft2,nfft3)))
         do k = 1, nfft3
            do j = 1, nfft2
               do i = 1, nfft1
                  work(1,i) = qgrid(1,i,j,k)
                  work(2,i) = qgrid(2,i,j,k)
               end do
               call cfftf (nfft1,work,ffttable(1,1),iprime(1,1))
               do i = 1, nfft1
                  qgrid(1,i,j,k) = work(1,i)
                  qgrid(2,i,j,k) = work(2,i)
               end do
            end do
         end do
         do k = 1, nfft3
            do i = 1, nfft1
               do j = 1, nfft2
                  work(1,j) = qgrid(1,i,j,k)
                  work(2,j) = qgrid(2,i,j,k)
               end do
               call cfftf (nfft2,work,ffttable(1,2),iprime(1,2))
               do j = 1, nfft2
                  qgrid(1,i,j,k) = work(1,j)
                  qgrid(2,i,j,k) = work(2,j)
               end do
            end do
         end do
         do i = 1, nfft1
            do j = 1, nfft2
               do k = 1, nfft3
                  work(1,k) = qgrid(1,i,j,k)
                  work(2,k) = qgrid(2,i,j,k)
               end do
               call cfftf (nfft3,work,ffttable(1,3),iprime(1,3))
               do k = 1, nfft3
                  qgrid(1,i,j,k) = work(1,k)
                  qgrid(2,i,j,k) = work(2,k)
               end do
            end do
         end do
         deallocate (work)
!$    end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine fftback  --  backward Fast Fourier transform  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "fftback" performs a 3-D FFT backward transform via a single
c     3-D transform or three separate 1-D transforms
c
c
      subroutine fftback
      use fft
      use pme
      implicit none
      integer i,j,k
      real*8, allocatable :: work(:,:)
c
c
c     perform a single 3-D backward transform using FFTW
c
!$    if (ffttyp .eq. 'FFTW') then
!$       call dfftw_execute_dft (planb,qgrid,qgrid)
!$    else
c
c     perform three 1-D backward transforms using FFTPACK
c
         allocate (work(2,max(nfft1,nfft2,nfft3)))
         do k = 1, nfft3
            do j = 1, nfft2
               do i = 1, nfft1
                  work(1,i) = qgrid(1,i,j,k)
                  work(2,i) = qgrid(2,i,j,k)
               end do
               call cfftb (nfft1,work,ffttable(1,1),iprime(1,1))
               do i = 1, nfft1
                  qgrid(1,i,j,k) = work(1,i)
                  qgrid(2,i,j,k) = work(2,i)
               end do
            end do
         end do
         do k = 1, nfft3
            do i = 1, nfft1
               do j = 1, nfft2
                  work(1,j) = qgrid(1,i,j,k)
                  work(2,j) = qgrid(2,i,j,k)
               end do
               call cfftb (nfft2,work,ffttable(1,2),iprime(1,2))
               do j = 1, nfft2
                  qgrid(1,i,j,k) = work(1,j)
                  qgrid(2,i,j,k) = work(2,j)
               end do
            end do
         end do
         do i = 1, nfft1
            do j = 1, nfft2
               do k = 1, nfft3
                  work(1,k) = qgrid(1,i,j,k)
                  work(2,k) = qgrid(2,i,j,k)
               end do
               call cfftb (nfft3,work,ffttable(1,3),iprime(1,3))
               do k = 1, nfft3
                  qgrid(1,i,j,k) = work(1,k)
                  qgrid(2,i,j,k) = work(2,k)
               end do
            end do
         end do
         deallocate (work)
!$    end if
      return
      end
