c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine cutoffs  --  distance cutoffs & neighbor lists  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "cutoffs" initializes and stores spherical energy cutoff
c     distance windows, Hessian element and Ewald sum cutoffs,
c     and allocates pairwise neighbor lists
c
c
      subroutine cutoffs
      use sizes
      use atoms
      use bound
      use hescut
      use keys
      use limits
      use neigh
      use polpot
      use tarray
      implicit none
      integer i,next
      integer limit
      real*8 big,value
      logical truncate
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     set defaults for spherical energy cutoff distances
c
      big = 1.0d12
      if (use_bounds) then
         vdwcut = 9.0d0
         ctcut =  4.5d0
         chgcut = 9.0d0
         dplcut = 9.0d0
         mpolecut = 9.0d0
      else
         vdwcut = big
         ctcut = big
         chgcut = big
         dplcut = big
         mpolecut = big
      end if
      ewaldcut = 7.0d0
      usolvcut = 4.5d0
      ctcut =  4.5d0
c
c     set defaults for tapering, Hessian cutoff and neighbor buffers
c
      vdwtaper = 0.90d0
      cttaper = 0.90d0
      chgtaper = 0.65d0
      dpltaper = 0.75d0
      mpoletaper = 0.65d0
      hesscut = 0.0d0
      lbuffer = 2.0d0
      pbuffer = 2.0d0
c
c     set defaults for Ewald sum, tapering style and neighbor method
c
      use_ewald = .false.
      truncate = .false.
      use_lights = .false.
      use_list = .false.
      use_vlist = .false.
      use_ctlist = .false.
      use_clist = .false.
      use_mlist = .false.
      use_ulist = .false.
      dovlst = .true.
      doctlst = .true.
      doclst = .true.
      domlst = .true.
      doulst = .true.
c
c     search the keywords for various cutoff parameters
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
c
c     get values related to use of Ewald summation
c
         if (keyword(1:6) .eq. 'EWALD ') then
            use_ewald = .true.
         else if (keyword(1:13) .eq. 'EWALD-CUTOFF ') then
            read (string,*,err=10,end=10)  ewaldcut
c
c     get cutoff for preconditioner of dipole solver
c
         else if (keyword(1:14) .eq. 'USOLVE-CUTOFF ') then
            read (string,*,err=10,end=10)  usolvcut
c
c     get values for the tapering style and neighbor method
c
         else if (keyword(1:9) .eq. 'TRUNCATE ') then
            truncate = .true.
         else if (keyword(1:7) .eq. 'LIGHTS ') then
            use_lights = .true.
         else if (keyword(1:14) .eq. 'NEIGHBOR-LIST ') then
            use_list = .true.
            use_vlist = .true.
            use_clist = .true.
            use_ctlist = .true.
            use_mlist = .true.
            use_ulist = .true.
         else if (keyword(1:9) .eq. 'VDW-LIST ') then
            use_list = .true.
            use_vlist = .true.
         else if (keyword(1:9) .eq. 'CHG-LIST ') then
            use_list = .true.
            use_clist = .true.
         else if (keyword(1:11) .eq. 'MPOLE-LIST ') then
            use_list = .true.
            use_mlist = .true.
            use_ulist = .true.
         else if (keyword(1:8) .eq. 'CT-LIST ') then
            use_list = .true.
            use_ctlist = .true.
c
c     get cutoff for the magnitude of Hessian elements
c
         else if (keyword(1:12) .eq. 'HESS-CUTOFF ') then
            read (string,*,err=10,end=10)  hesscut
c
c     get the cutoff radii for potential energy functions
c
         else if (keyword(1:7) .eq. 'CUTOFF ') then
            read (string,*,err=10,end=10)  value
            vdwcut = value
            ctcut = value
            chgcut = value
            dplcut = value
            mpolecut = value
            ewaldcut = value
         else if (keyword(1:11) .eq. 'VDW-CUTOFF ') then
            read (string,*,err=10,end=10)  vdwcut
         else if (keyword(1:10) .eq. 'CT-CUTOFF ') then
            read (string,*,err=10,end=10)  ctcut
         else if (keyword(1:11) .eq. 'CHG-CUTOFF ') then
            read (string,*,err=10,end=10)  chgcut
         else if (keyword(1:11) .eq. 'DPL-CUTOFF ') then
            read (string,*,err=10,end=10)  dplcut
         else if (keyword(1:13) .eq. 'MPOLE-CUTOFF ') then
            read (string,*,err=10,end=10)  mpolecut
c
c     get distance for initialization of energy switching
c
         else if (keyword(1:6) .eq. 'TAPER ') then
            read (string,*,err=10,end=10)  value
            vdwtaper = value
            chgtaper = value
            dpltaper = value
            mpoletaper = value
         else if (keyword(1:10) .eq. 'VDW-TAPER ') then
            read (string,*,err=10,end=10)  vdwtaper
         else if (keyword(1:10) .eq. 'CHG-TAPER ') then
            read (string,*,err=10,end=10)  chgtaper
         else if (keyword(1:10) .eq. 'DPL-TAPER ') then
            read (string,*,err=10,end=10)  dpltaper
         else if (keyword(1:12) .eq. 'MPOLE-TAPER ') then
            read (string,*,err=10,end=10)  mpoletaper
c
c     get buffer width for use with pairwise neighbor lists
c
         else if (keyword(1:12) .eq. 'LIST-BUFFER ') then
            read (string,*,err=10,end=10)  lbuffer
         else if (keyword(1:14) .eq. 'USOLVE-BUFFER ') then
            read (string,*,err=10,end=10)  pbuffer
         end if
   10    continue
      end do
c
c     preconditioner list only needed for mutual polarization
c
      if (poltyp .eq. 'DIRECT')  use_ulist = .false.
      if (use_list)  usolvcut = usolvcut - pbuffer
c
c     apply any Ewald cutoff to charge and multipole terms
c
      if (use_ewald) then
         chgcut = ewaldcut
         mpolecut = ewaldcut
      end if
c
c     convert any tapering percentages to absolute distances
c
      if (vdwtaper .lt. 1.0d0)  vdwtaper = vdwtaper * vdwcut
      if (cttaper .lt. 1.0d0)  cttaper = cttaper * ctcut
      if (chgtaper .lt. 1.0d0)  chgtaper = chgtaper * chgcut
      if (dpltaper .lt. 1.0d0)  dpltaper = dpltaper * dplcut
      if (mpoletaper .lt. 1.0d0)  mpoletaper = mpoletaper * mpolecut
c
c     apply truncation cutoffs if they were requested
c
      if (truncate) then
         vdwtaper = big
         cttaper = big
         chgtaper = big
         dpltaper = big
         mpoletaper = big
      end if
c
c     set buffer region limits for pairwise neighbor lists
c
      lbuf2 = (0.5d0*lbuffer)**2
      pbuf2 = (0.5d0*pbuffer)**2
      vbuf2 = (vdwcut+lbuffer)**2
      ctbuf2 = (ctcut+lbuffer)**2
      cbuf2 = (chgcut+lbuffer)**2
      mbuf2 = (mpolecut+lbuffer)**2
      ubuf2 = (usolvcut+pbuffer)**2
      vbufx = (vdwcut+2.0d0*lbuffer)**2
      ctbufx = (ctcut+2.0d0*lbuffer)**2
      cbufx = (chgcut+2.0d0*lbuffer)**2
      mbufx = (mpolecut+2.0d0*lbuffer)**2
      ubufx = (usolvcut+2.0d0*pbuffer)**2
c
c     specify maximum size for each of the neighbor lists
c
      maxvlst = 2500
      if (vdwcut .ne. big) then
         limit = int(sqrt(vbuf2)**3) + 100
         maxvlst = min(limit,maxvlst)
      end if
      maxctlst = 2500
      if (ctcut .ne. big) then
         limit = int(sqrt(ctbuf2)**3) + 100
         maxctlst = min(limit,maxctlst)
      end if
      maxelst = 2500
      if (chgcut .ne. big) then
         limit = int(sqrt(cbuf2)**3) + 100
         maxelst = min(limit,maxelst)
      end if
      if (mpolecut .ne. big) then
         limit = int(sqrt(mbuf2)**3) + 100
         maxelst = min(limit,maxelst)
      end if
      maxulst = 500
      limit = int(sqrt(ubuf2)**3) + 100
      maxulst = min(limit,maxulst)
c
c     perform dynamic allocation of some global arrays
c
      if (use_vlist) then
         if (allocated(nvlst))  deallocate (nvlst)
         if (allocated(vlst))  deallocate (vlst)
         if (allocated(xvold))  deallocate (xvold)
         if (allocated(yvold))  deallocate (yvold)
         if (allocated(zvold))  deallocate (zvold)
         allocate (nvlst(n))
         allocate (vlst(maxvlst,n))
         allocate (xvold(n))
         allocate (yvold(n))
         allocate (zvold(n))
      end if
      if (use_ctlist) then
         if (allocated(nctlst))  deallocate (nctlst)
         if (allocated(ctlst))  deallocate (ctlst)
         if (allocated(xctold))  deallocate (xctold)
         if (allocated(yctold))  deallocate (yctold)
         if (allocated(zctold))  deallocate (zctold)
         allocate (nctlst(n))
         allocate (ctlst(maxctlst,n))
         allocate (xctold(n))
         allocate (yctold(n))
         allocate (zctold(n))
      end if
      if (use_clist .or. use_mlist) then
         if (allocated(nelst))  deallocate (nelst)
         if (allocated(elst))  deallocate (elst)
         allocate (nelst(n))
         allocate (elst(maxelst,n))
      end if
      if (use_clist) then
         if (allocated(xcold))  deallocate (xcold)
         if (allocated(ycold))  deallocate (ycold)
         if (allocated(zcold))  deallocate (zcold)
         allocate (xcold(n))
         allocate (ycold(n))
         allocate (zcold(n))
      end if
      if (use_mlist) then
         if (allocated(xmold))  deallocate (xmold)
         if (allocated(ymold))  deallocate (ymold)
         if (allocated(zmold))  deallocate (zmold)
         allocate (xmold(n))
         allocate (ymold(n))
         allocate (zmold(n))
         if (poltyp .ne. 'DIRECT') then
            if (allocated(tindex))  deallocate (tindex)
            if (allocated(tdipdip))  deallocate (tdipdip)
            allocate (tindex(2,n*maxelst))
            allocate (tdipdip(6,n*maxelst))
         end if
      end if
      if (use_ulist) then
         if (allocated(nulst))  deallocate (nulst)
         if (allocated(ulst))  deallocate (ulst)
         if (allocated(xuold))  deallocate (xuold)
         if (allocated(yuold))  deallocate (yuold)
         if (allocated(zuold))  deallocate (zuold)
         allocate (nulst(n))
         allocate (ulst(maxulst,n))
         allocate (xuold(n))
         allocate (yuold(n))
         allocate (zuold(n))
      end if
      return
      end
