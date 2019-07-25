c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine ksolv  --  solvation parameter assignment  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "ksolv" assigns implicit solvation energy parameters for
c     the surface area, generalized Born, generalized Kirkwood,
c     Poisson-Boltzmann, cavity-dispersion and HPMF models
c
c
      subroutine ksolv
      use sizes
      use gkstuf
      use inform
      use iounit
      use keys
      use potent
      use solute
      implicit none
      integer i,k,next
      real*8 rd
      logical header
      character*20 keyword
      character*20 value
      character*240 record
      character*240 string
c
c
c     defaults for implicit solvation term and parameters
c
      use_solv = .false.
      use_born = .false.
      solvtyp = '        '
      borntyp = '        '
      doffset = 0.09d0
c
c     search keywords for implicit solvation commands
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:8) .eq. 'SOLVATE ') then
            use_solv = .true.
            use_born = .false.
            solvtyp = 'ASP'
            call getword (record,value,next)
            call upcase (value)
            if (value(1:3) .eq. 'ASP') then
               solvtyp = 'ASP'
            else if (value(1:4) .eq. 'SASA') then
               solvtyp = 'SASA'
            else if (value(1:5) .eq. 'ONION') then
               use_born = .true.
               solvtyp = 'ONION'
            else if (value(1:5) .eq. 'STILL') then
               use_born = .true.
               solvtyp = 'STILL'
            else if (value(1:3) .eq. 'HCT') then
               use_born = .true.
               solvtyp = 'HCT'
            else if (value(1:3) .eq. 'OBC') then
               use_born = .true.
               solvtyp = 'OBC'
            else if (value(1:3) .eq. 'ACE') then
               use_born = .true.
               solvtyp = 'ACE'
            else if (value(1:7) .eq. 'GB-HPMF') then
               use_born = .true.
               solvtyp = 'GB-HPMF'
            else if (value(1:2) .eq. 'GB') then
               use_born = .true.
               solvtyp = 'STILL'
            else if (value(1:7) .eq. 'GK-HPMF') then
               use_born = .true.
               solvtyp = 'GK-HPMF'
            else if (value(1:2) .eq. 'GK') then
               use_born = .true.
               solvtyp = 'GK'
            else if (value(1:7) .eq. 'PB-HPMF') then
               solvtyp = 'PB-HPMF'
            else if (value(1:2) .eq. 'PB') then
               solvtyp = 'PB'
            end if
         else if (keyword(1:12) .eq. 'BORN-RADIUS ') then
            call getword (record,value,next)
            call upcase (value)
            if (value(1:5) .eq. 'ONION') then
               borntyp = 'ONION'
            else if (value(1:5) .eq. 'STILL') then
               borntyp = 'STILL'
            else if (value(1:3) .eq. 'HCT') then
               borntyp = 'HCT'
            else if (value(1:3) .eq. 'OBC') then
               borntyp = 'OBC'
            else if (value(1:3) .eq. 'ACE') then
               borntyp = 'ACE'
            else if (value(1:4) .eq. 'GRYCUK') then
               borntyp = 'GRYCUK'
            else if (value(1:7) .eq. 'PERFECT') then
               borntyp = 'PERFECT'
            end if
         else if (keyword(1:18) .eq. 'DIELECTRIC-OFFSET ') then
            read (string,*,err=50,end=50)  doffset
            if (doffset .lt. 0.0d0)  doffset = -doffset
         else if (keyword(1:4) .eq. 'GKR ') then
            k = 0
            rd = 0.0d0
            read (string,*,err=10,end=10)  k,rd
   10       continue
            if (k.ge.1 .and. k.le.maxclass) then
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,20)
   20             format (/,' Additional Implicit Solvation',
     &                       ' Parameters:',
     &                    //,5x,'Atom Type',10x,'Radius',/)
               end if
               if (rd .lt. 0.0d0)  rd = 0.0d0
               gkr(k) = rd
               if (.not. silent) then
                  write (iout,30)  k,rd
   30             format (4x,i6,8x,f12.4)
               end if
            else if (k .gt. maxtyp) then
               write (iout,40)  maxtyp
   40          format (/,' KSOLV  --  Only Atom Types Through',i4,
     &                    ' are Allowed')
               abort = .true.
            end if
         end if
   50    continue
      end do
c
c     set a default if no Born radius method was assigned
c
      if (use_born .and. borntyp.eq.'       ') then
         borntyp = solvtyp
         if (solvtyp .eq. 'GB-HPMF')  borntyp = 'STILL'
         if (solvtyp .eq. 'GK')  borntyp = 'GRYCUK'
         if (solvtyp .eq. 'GK-HPMF')  borntyp = 'GRYCUK'
      end if
c
c     invoke the setup needed for specific Born radius models
c
      if (borntyp .eq. 'PERFECT')  call kpb
c
c     invoke the setup needed for specific solvation models
c
      if (solvtyp.eq.'ASP' .or. solvtyp.eq.'SASA') then
         call ksa
      else if (solvtyp .eq. 'GB-HPMF') then
         call kgb
         call khpmf
      else if (solvtyp .eq. 'GK-HPMF') then
         call kgk
         call khpmf
      else if (solvtyp .eq. 'GK') then
         call kgk
         call knp
      else if (solvtyp .eq. 'PB-HPMF') then
         call kpb
         call khpmf
      else if (solvtyp .eq. 'PB') then
         call kpb
         call knp
      else
         call kgb
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine ksa  --  set surface area solvation parameters  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ksa" initializes parameters needed for surface area-based
c     implicit solvation models including ASP and SASA
c
c     literature references:
c
c     L. Wesson and D. Eisenberg, "Atomic Solvation Parameters
c     Applied to Molecular Dynamics of Proteins in Solution",
c     Protein Science, 1, 227-235 (1992)  (Eisenberg-McLachlan ASP)
c
c     T. Ooi, M. Oobatake, G. Nemethy and H. A. Scheraga, "Accessible
c     Surface Areas as a Measure of the Thermodynamic Parameters of
c     Hydration of Peptides", PNAS, 84, 3086-3090 (1987)  (SASA)
c
c
      subroutine ksa
      use sizes
      use atomid
      use atoms
      use couple
      use solute
      implicit none
      integer i,j,k
      integer atmnum
c
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(rsolv))  deallocate (rsolv)
      if (allocated(asolv))  deallocate (asolv)
      allocate (rsolv(n))
      allocate (asolv(n))
c
c     assign the Eisenberg-McLachlan ASP solvation parameters;
c     parameters only available for protein-peptide groups
c
      if (solvtyp .eq. 'ASP') then
         do i = 1, n
            atmnum = atomic(i)
            if (atmnum .eq. 6) then
               rsolv(i) = 1.9d0
               asolv(i) = 0.004d0
            else if (atmnum .eq. 7) then
               rsolv(i) = 1.7d0
               asolv(i) = -0.113d0
               if (n12(i) .eq. 4) then
                  asolv(i) = -0.169d0
               end if
            else if (atmnum .eq. 8) then
               rsolv(i) = 1.4d0
               asolv(i) = -0.113d0
               if (n12(i).eq.1 .and. atomic(i12(1,i)).eq.6) then
                  do j = 1, n13(i)
                     k = i13(j,i)
                     if (n12(k).eq.1 .and. atomic(k).eq.8) then
                        asolv(i) = -0.166d0
                     end if
                  end do
               end if
               do j = 1, n12(i)
                  k = i12(j,i)
                  if (atomic(k) .eq. 15)  asolv(i) = -0.140d0
               end do
            else if (atmnum .eq. 15) then
               rsolv(i) = 1.9d0
               asolv(i) = -0.140d0
            else if (atmnum .eq. 16) then
               rsolv(i) = 1.8d0
               asolv(i) = -0.017d0
            else
               rsolv(i) = 0.0d0
               asolv(i) = 0.0d0
            end if
         end do
      end if
c
c     assign the Ooi-Scheraga SASA solvation parameters;
c     parameters only available for protein-peptide groups
c
      if (solvtyp .eq. 'SASA') then
         do i = 1, n
            atmnum = atomic(i)
            if (atmnum .eq. 6) then
               rsolv(i) = 2.0d0
               asolv(i) = 0.008d0
               if (n12(i) .eq. 3) then
                  rsolv(i) = 1.75d0
                  asolv(i) = -0.008d0
                  do j = 1, n12(i)
                     k = i12(j,i)
                     if (atomic(k) .eq. 8) then
                        rsolv(i) = 1.55d0
                        asolv(i) = 0.427d0
                     end if
                  end do
               end if
            else if (atmnum .eq. 7) then
               rsolv(i) = 1.55d0
               asolv(i) = -0.132d0
               if (n12(i) .eq. 4)  asolv(i) = -1.212d0
            else if (atmnum .eq. 8) then
               rsolv(i) = 1.4d0
               if (n12(i) .eq. 1) then
                  asolv(i) = -0.038d0
                  if (atomic(i12(1,i)) .eq. 6) then
                     do j = 1, n13(i)
                        k = i13(j,i)
                        if (n12(k).eq.1 .and. atomic(k).eq.8) then
                           asolv(i) = -0.770d0
                        end if
                     end do
                  end if
               else if (n12(i) .eq. 2) then
                  asolv(i) = -0.172d0
               end if
               do j = 1, n12(i)
                  k = i12(j,i)
                  if (atomic(k) .eq. 15)  asolv(i) = -0.717d0
               end do
            else if (atmnum .eq. 15) then
               rsolv(i) = 2.1d0
               asolv(i) = 0.0d0
            else if (atmnum .eq. 16) then
               rsolv(i) = 2.0d0
               asolv(i) = -0.021d0
            else if (atmnum .eq. 17) then
               rsolv(i) = 2.0d0
               asolv(i) = 0.012d0
            else
               rsolv(i) = 0.0d0
               asolv(i) = 0.0d0
            end if
         end do
      end if
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine kgb  --  assign generalized Born parameters  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "kgb" initializes parameters needed for the generalized
c     Born implicit solvation models
c
c     literature references:
c
c     M. Schaefer, C. Bartels, F. Leclerc and M. Karplus, "Effective
c     Atom Volumes for Implicit Solvent Models: Comparison between
c     Voronoi Volumes and Minimum Fluctuations Volumes", Journal of
c     Computational Chemistry, 22, 1857-1879 (2001)  (ACE)
c
c
      subroutine kgb
      use sizes
      use angbnd
      use atmlst
      use atomid
      use atoms
      use bndstr
      use chgpot
      use couple
      use math
      use potent
      use solute
      implicit none
      integer i,j,k,m
      integer mm,nh,kc
      integer ia,ib,ic,id
      integer atmnum,atmmas
      real*8 ri,ri2,rk,rk2
      real*8 c1,c2,c3,pi2
      real*8 r,r2,r4,rab,rbc
      real*8 cosine,factor
      real*8 h,ratio,term
      real*8 width,qterm,temp
      real*8 alpha,alpha2,alpha4
      real*8 vk,prod2,prod4
      real*8 fik,tik2,qik,uik
      real*8 s2ik,s3ik,omgik
      logical amide
c
c
c     perform dynamic allocation of some global arrays
c
      if (.not. allocated(wace))  allocate (wace(maxclass,maxclass)) 
      if (.not. allocated(s2ace))  allocate (s2ace(maxclass,maxclass)) 
      if (.not. allocated(uace))  allocate (uace(maxclass,maxclass)) 
      if (allocated(rsolv))  deallocate (rsolv)
      if (allocated(asolv))  deallocate (asolv)
      if (allocated(rborn))  deallocate (rborn)
      if (allocated(drb))  deallocate (drb)
      if (allocated(drobc))  deallocate (drobc)
      if (allocated(gpol))  deallocate (gpol)
      if (allocated(shct))  deallocate (shct)
      if (allocated(aobc))  deallocate (aobc)
      if (allocated(bobc))  deallocate (bobc)
      if (allocated(gobc))  deallocate (gobc)
      if (allocated(vsolv))  deallocate (vsolv)
      allocate (rsolv(n))
      allocate (asolv(n))
      allocate (rborn(n))
      allocate (drb(n))
      allocate (drobc(n))
      allocate (gpol(n))
      allocate (shct(n))
      allocate (aobc(n))
      allocate (bobc(n))
      allocate (gobc(n))
      allocate (vsolv(n))
c
c     set offset and scaling values for analytical Still method
c
      if (borntyp .eq. 'STILL') then
         p1 = 0.073d0
         p2 = 0.921d0
         p3 = 6.211d0
         p4 = 15.236d0
         p5 = 1.254d0
         if (.not. use_bond)  call kbond
         if (.not. use_angle)  call kangle
      end if
c
c     set overlap scale factors for HCT and OBC methods
c
      if (borntyp.eq.'HCT' .or. borntyp.eq.'OBC') then
         do i = 1, n
            shct(i) = 0.80d0
            atmnum = atomic(i)
            if (atmnum .eq. 1)  shct(i) = 0.85d0
            if (atmnum .eq. 6)  shct(i) = 0.72d0
            if (atmnum .eq. 7)  shct(i) = 0.79d0
            if (atmnum .eq. 8)  shct(i) = 0.85d0
            if (atmnum .eq. 9)  shct(i) = 0.88d0
            if (atmnum .eq. 15)  shct(i) = 0.86d0
            if (atmnum .eq. 16)  shct(i) = 0.96d0
            if (atmnum .eq. 26)  shct(i) = 0.88d0
         end do
      end if
c
c     set rescaling coefficients for the OBC method
c
      if (borntyp .eq. 'OBC') then
         do i = 1, n
            aobc(i) = 1.00d0
            bobc(i) = 0.80d0
            gobc(i) = 4.85d0
         end do
      end if
c
c     set the Gaussian width factor for the ACE method
c
      if (borntyp .eq. 'ACE') then
         width = 1.2d0
      end if
c
c     assign surface area factors for nonpolar solvation
c
      if (borntyp .eq. 'ONION') then
         do i = 1, n
            asolv(i) = 0.0072d0
         end do
      else if (borntyp .eq. 'STILL') then
         do i = 1, n
            asolv(i) = 0.0049d0
         end do
      else if (borntyp .eq. 'HCT') then
         do i = 1, n
            asolv(i) = 0.0054d0
         end do
      else if (borntyp .eq. 'OBC') then
         do i = 1, n
            asolv(i) = 0.0054d0
         end do
      else if (borntyp .eq. 'ACE') then
         do i = 1, n
            asolv(i) = 0.0030d0
         end do
      end if
c
c     assign standard radii for GB/SA methods other than ACE;
c     taken from Macromodel and OPLS-AA, except for hydrogens
c
      if (borntyp .ne. 'ACE') then
         do i = 1, n
            atmnum = atomic(i)
            if (atmnum .eq. 1) then
               rsolv(i) = 1.25d0
               k = i12(1,i)
               if (atomic(k) .eq. 7)  rsolv(i) = 1.15d0
               if (atomic(k) .eq. 8)  rsolv(i) = 1.05d0
            else if (atmnum .eq. 3) then
               rsolv(i) = 1.432d0
            else if (atmnum .eq. 6) then
               rsolv(i) = 1.90d0
               if (n12(i) .eq. 3)  rsolv(i) = 1.875d0
               if (n12(i) .eq. 2)  rsolv(i) = 1.825d0
            else if (atmnum .eq. 7) then
               rsolv(i) = 1.7063d0
               if (n12(i) .eq. 4)  rsolv(i) = 1.625d0
               if (n12(i) .eq. 1)  rsolv(i) = 1.60d0
            else if (atmnum .eq. 8) then
               rsolv(i) = 1.535d0
               if (n12(i) .eq. 1)  rsolv(i) = 1.48d0
            else if (atmnum .eq. 9) then
               rsolv(i) = 1.47d0
            else if (atmnum .eq. 10) then
               rsolv(i) = 1.39d0
            else if (atmnum .eq. 11) then
               rsolv(i) = 1.992d0
            else if (atmnum .eq. 12) then
               rsolv(i) = 1.70d0
            else if (atmnum .eq. 14) then
               rsolv(i) = 1.80d0
            else if (atmnum .eq. 15) then
               rsolv(i) = 1.87d0
            else if (atmnum .eq. 16) then
               rsolv(i) = 1.775d0
            else if (atmnum .eq. 17) then
               rsolv(i) = 1.735d0
            else if (atmnum .eq. 18) then
               rsolv(i) = 1.70d0
            else if (atmnum .eq. 19) then
               rsolv(i) = 2.123d0
            else if (atmnum .eq. 20) then
               rsolv(i) = 1.817d0
            else if (atmnum .eq. 35) then
               rsolv(i) = 1.90d0
            else if (atmnum .eq. 36) then
               rsolv(i) = 1.812d0
            else if (atmnum .eq. 37) then
               rsolv(i) = 2.26d0
            else if (atmnum .eq. 53) then
               rsolv(i) = 2.10d0
            else if (atmnum .eq. 54) then
               rsolv(i) = 1.967d0
            else if (atmnum .eq. 55) then
               rsolv(i) = 2.507d0
            else if (atmnum .eq. 56) then
               rsolv(i) = 2.188d0
            else
               rsolv(i) = 2.0d0
            end if
         end do
      end if
c
c     compute the atomic volumes for the analytical Still method
c
      if (borntyp .eq. 'STILL') then
         do i = 1, n
            vsolv(i) = (4.0d0*pi/3.0d0) * rsolv(i)**3
            ri = rsolv(i)
            ri2 = ri * ri
            do j = 1, n12(i)
               k = i12(j,i)
               rk = rsolv(k)
               r = 1.01d0 * bl(bndlist(j,i))
               ratio = (rk*rk-ri2-r*r) / (2.0d0*ri*r)
               h = ri * (1.0d0+ratio)
               term = (pi/3.0d0) * h * h * (3.0d0*ri-h)
               vsolv(i) = vsolv(i) - term
            end do
         end do
c
c     get self-, 1-2 and 1-3 polarization for analytical Still method
c
         do i = 1, n
            gpol(i) = -0.5d0 * electric / (rsolv(i)-doffset+p1)
         end do
         do i = 1, nbond
            ia = ibnd(1,i)
            ib = ibnd(2,i)
            r = bl(i)
            r4 = r**4
            gpol(ia) = gpol(ia) + p2*vsolv(ib)/r4
            gpol(ib) = gpol(ib) + p2*vsolv(ia)/r4
         end do
         do i = 1, nangle
            ia = iang(1,i)
            ib = iang(2,i)
            ic = iang(3,i)
            factor = 1.0d0
            do j = 1, n12(ia)
               id = i12(j,ia)
               if (id .eq. ic) then
                  factor = 0.0d0
               else if (id .ne. ib) then
                  do k = 1, n12(ic)
                     if (i12(k,ic) .eq. id) then
                        factor = 0.5d0
                     end if
                  end do
               end if
            end do
            do j = 1, n12(ib)
               if (i12(j,ib) .eq. ia) then
                  rab = bl(bndlist(j,ib))
               else if (i12(j,ib) .eq. ic) then
                  rbc = bl(bndlist(j,ib))
               end if
            end do
            cosine = cos(anat(i)/radian)
            r2 = rab**2 + rbc**2 - 2.0d0*rab*rbc*cosine
            r4 = r2 * r2
            gpol(ia) = gpol(ia) + factor*p3*vsolv(ic)/r4
            gpol(ic) = gpol(ic) + factor*p3*vsolv(ia)/r4
         end do
      end if
c
c     assign the atomic radii and volumes for the ACE method;
c     volumes taken from average Voronoi values with hydrogens
c
      if (borntyp .eq. 'ACE') then
         do i = 1, n
            atmnum = atomic(i)
            atmmas = nint(mass(i))
            if (atmnum .eq. 1) then
               rsolv(i) = 1.468d0
               vsolv(i) = 11.0d0
               k = i12(1,i)
               if (atomic(k).eq.6 .and. n12(k).eq.4) then
                  vsolv(i) = 11.895d0
               else if (atomic(k).eq.6 .and. n12(k).eq.3) then
                  vsolv(i) = 13.242d0
               else if (atomic(k).eq.7 .and. n12(k).eq.4) then
                  rsolv(i) = 0.60d0
                  vsolv(i) = 9.138d0
               else if (atomic(k).eq.7 .or. atomic(k).eq.8) then
                  rsolv(i) = 0.60d0
                  vsolv(i) = 9.901d0
               else if (atomic(k).ne.16) then
                  rsolv(i) = 1.468d0
                  vsolv(i) = 13.071d0
               end if
            else if (atmnum .eq. 6) then
               rsolv(i) = 2.49d0
               vsolv(i) = 7.0d0
               nh = 0
               do j = 1, n12(i)
                  k = i12(j,i)
                  if (atomic(k) .eq. 1)  nh = nh + 1
               end do
               if (n12(i) .eq. 4) then
                  if (nh .eq. 3) then
                     vsolv(i) = 3.042d0
                  else if (nh .eq. 2) then
                     vsolv(i) = 3.743d0
                  else if (nh .eq. 1) then
                     vsolv(i) = 4.380d0
                  end if
               else if (n12(i) .eq. 3) then
                  if (nh .eq. 1) then
                     rsolv(i) = 2.10d0
                     vsolv(i) = 7.482d0
                  else if (nh .eq. 0) then
                     rsolv(i) = 2.10d0
                     vsolv(i) = 8.288d0
                  end if
                  do j = 1, n12(i)
                     k = i12(1,j)
                     if (atomic(k).eq.8 .and. n12(k).eq.1) then
                        rsolv(i) = 2.10d0
                        vsolv(i) = 7.139d0
                     end if
                  end do
               end if
               if (atmmas .eq. 15) then
                  rsolv(i) = 2.165d0
                  vsolv(i) = 33.175d0
               else if (atmmas .eq. 14) then
                  rsolv(i) = 2.235d0
                  vsolv(i) = 20.862d0
               else if (atmmas.eq.13 .and. n12(i).eq.2) then
                  rsolv(i) = 2.10d0
                  vsolv(i) = 20.329d0
               else if (atmmas .eq. 13) then
                  rsolv(i) = 2.365d0
                  vsolv(i) = 11.784d0
               end if
            else if (atmnum .eq. 7) then
               rsolv(i) = 1.60d0
               vsolv(i) = 6.0d0
               nh = 0
               do j = 1, n12(i)
                  k = i12(j,i)
                  if (atomic(k) .eq. 1)  nh = nh + 1
               end do
               if (n12(i) .eq. 4) then
                  if (nh .eq. 3) then
                     vsolv(i) = 2.549d0
                  else if (nh .eq. 2) then
                     vsolv(i) = 3.304d0
                  end if
               else if (n12(i) .eq. 3) then
                  amide = .false.
                  do j = 1, n12(i)
                     m = i12(j,i)
                     if (atomic(m) .eq. 6) then
                        do k = 1, n12(m)
                           mm = i12(k,m)
                           if (atomic(mm).eq.8 .and. n12(mm).eq.1) then
                              amide = .true.
                           end if
                        end do
                     end if
                  end do
                  if (amide) then
                     if (nh .eq. 0) then
                        vsolv(i) = 7.189d0
                     else if (nh .eq. 1) then
                        vsolv(i) = 6.030d0
                     else if (nh .eq. 2) then
                        vsolv(i) = 5.693d0
                     end if
                  else
                     if (nh .eq. 2) then
                        vsolv(i) = 5.677d0
                     else if (nh .eq. 2) then
                        vsolv(i) = 6.498d0
                     end if
                  end if
               end if
            else if (atmnum .eq. 8) then
               rsolv(i) = 1.60d0
               vsolv(i) = 12.0d0
               if (n12(i) .eq. 1) then
                  vsolv(i) = 13.532d0
                  k = i12(1,i)
                  if (atomic(k) .eq. 15) then
                     vsolv(i) = 17.202d0
                  else
                     do j = 1, n13(i)
                        k = i13(j,i)
                        if (atomic(j).eq.8 .and. n12(j).eq.1) then
                           vsolv(i) = 15.400d0
                        end if
                     end do
                  end if
               else if (n12(i) .eq. 2) then
                  vsolv(i) = 10.642d0
                  do j = 1, n12(i)
                     k = i12(j,i)
                     if (atomic(k) .eq. 15)  vsolv(i) = 11.416d0
                  end do
               end if
            else if (atmnum .eq. 12) then
               rsolv(i) = 1.0d0
               vsolv(i) = 15.235d0
            else if (atmnum .eq. 15) then
               rsolv(i) = 1.89d0
               vsolv(i) = 6.131d0
            else if (atmnum .eq. 16) then
               rsolv(i) = 1.89d0
               vsolv(i) = 17.232d0
               do j = 1, n12(i)
                  k = i12(j,i)
                  if (atomic(k) .eq. 16)  vsolv(i) = 18.465d0
               end do
            else if (atmnum .eq. 26) then
               rsolv(i) = 0.65d0
               vsolv(i) = 9.951d0
            else
               rsolv(i) = 0.0d0
               vsolv(i) = 0.0d0
            end if
         end do
c
c     calculate the pairwise parameters for the ACE method
c
         c1 = 4.0d0 / (3.0d0*pi)
         c2 = 77.0d0 * pi * sqrttwo / 512.0d0
         c3 = 2.0d0 * pi * sqrtpi
         pi2 = 1.0d0 / (pi*pi)
         do i = 1, n
            ic = class(i)
            ri = rsolv(i)
            ri2 = ri * ri
            do k = 1, n
               kc = class(k)
               rk = rsolv(kc)
               vk = vsolv(kc)
               rk2 = rk * rk
               alpha = max(width,ri/rk)
               alpha2 = alpha * alpha
               alpha4 = alpha2 * alpha2
               prod2 = alpha2 * rk2
               prod4 = prod2 * prod2
               ratio = alpha2 * rk2 / ri2
               tik2 = 0.5d0 * pi * ratio
               temp = 1.0d0 / (1.0d0+2.0d0*tik2)
               fik = 2.0d0/(1.0d0+tik2) - temp
               qik = tik2 * sqrt(temp)
               qterm = qik - atan(qik)
               if (k .ne. i) then
                  omgik = vk * qterm * pi2 / prod4
               else
                  omgik = c1 * qterm / (alpha4 * ri)
               end if
               s2ik = 3.0d0 * qterm * prod2
     &                   / ((3.0d0+fik)*qik-4.0d0*atan(qik))
               s3ik = s2ik * sqrt(s2ik)
               uik = c2 * ri / (1.0d0-(c3*s3ik*ri*omgik/vk))
               wace(ic,kc) = omgik
               s2ace(ic,kc) = s2ik
               uace(ic,kc) = uik
            end do
         end do
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kgk  --  set generalized Kirkwood parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kgk" initializes parameters needed for the generalized
c     Kirkwood implicit solvation model
c
c
      subroutine kgk
      use sizes
      use atomid
      use atoms
      use couple
      use gkstuf
      use keys
      use kvdws
      use polar
      use ptable
      use solute
      implicit none
      integer i,j,k,l,m
      integer atmnum,next
      real*8 rscale
      real*8 offset
      character*10 radtyp
      character*20 keyword
      character*20 value
      character*240 record
      character*240 string
c
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(rsolv))  deallocate (rsolv)
      if (allocated(rborn))  deallocate (rborn)
      if (allocated(drb))  deallocate (drb)
      if (allocated(drbp))  deallocate (drbp)
      if (allocated(drobc))  deallocate (drobc)
      if (allocated(shct))  deallocate (shct)
      if (allocated(udirs))  deallocate (udirs)
      if (allocated(udirps))  deallocate (udirps)
      if (allocated(uinds))  deallocate (uinds)
      if (allocated(uinps))  deallocate (uinps)
      if (allocated(uopts))  deallocate (uopts)
      if (allocated(uoptps))  deallocate (uoptps)
      allocate (rsolv(n))
      allocate (rborn(n))
      allocate (drb(n))
      allocate (drbp(n))
      allocate (drobc(n))
      allocate (shct(n))
      allocate (udirs(3,n))
      allocate (udirps(3,n))
      allocate (uinds(3,n))
      allocate (uinps(3,n))
      allocate (uopts(0:coptmax,3,n))
      allocate (uoptps(0:coptmax,3,n))
c
c     set default value for exponent in the GB/GK function
c
      gkc = 2.455d0
      radtyp = 'BONDI'
c
c     get any altered generalized Kirkwood values from keyfile
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:4) .eq. 'GKC ') then
            read (string,*,err=10,end=10)  gkc
         else if (keyword(1:10) .eq. 'GK-RADIUS ') then
            call getword (record,value,next)
            call upcase (value)
            if (value(1:3) .eq. 'VDW') then
               radtyp = 'VDW'
            else if (value(1:10) .eq. 'MACROMODEL') then
               radtyp = 'MACROMODEL'
            else if (value(1:6) .eq. 'AMOEBA') then
               radtyp = 'AMOEBA'
            else if (value(1:5) .eq. 'BONDI') then
               radtyp = 'BONDI'
            else if (value(1:6) .eq. 'TOMASI') then
               radtyp = 'TOMASI'
            end if
         end if
   10    continue
      end do
c
c     assign base atomic radii from the van der Waals values
c
      if (radtyp .eq. 'VDW') then
         do i = 1, n
            rsolv(i) = 2.0d0
            if (class(i) .ne. 0)  rsolv(i) = rad(class(i))
            rsolv(i) = rsolv(i) - 0.10d0
         end do
c
c     assign standard solvation radii adapted from Macromodel
c
      else if (radtyp .eq. 'MACROMODEL') then
         do i = 1, n
            atmnum = atomic(i)
            if (atmnum .eq. 0)  rsolv(i) = 0.0d0
            rsolv(i) = vdwrad(atmnum)
            if (atmnum .eq. 1) then
               rsolv(i) = 1.25d0
               k = i12(1,i)
               if (atomic(k) .eq. 7)  rsolv(i) = 1.15d0
               if (atomic(k) .eq. 8)  rsolv(i) = 1.05d0
            else if (atmnum .eq. 3) then
               rsolv(i) = 1.432d0
            else if (atmnum .eq. 6) then
               rsolv(i) = 1.90d0
               if (n12(i) .eq. 3)  rsolv(i) = 1.875d0
               if (n12(i) .eq. 2)  rsolv(i) = 1.825d0
            else if (atmnum .eq. 7) then
               rsolv(i) = 1.7063d0
               if (n12(i) .eq. 4)  rsolv(i) = 1.625d0
               if (n12(i) .eq. 1)  rsolv(i) = 1.60d0
            else if (atmnum .eq. 8) then
               rsolv(i) = 1.535d0
               if (n12(i) .eq. 1)  rsolv(i) = 1.48d0
            else if (atmnum .eq. 9) then
               rsolv(i) = 1.47d0
            else if (atmnum .eq. 10) then
               rsolv(i) = 1.39d0
            else if (atmnum .eq. 11) then
               rsolv(i) = 1.992d0
            else if (atmnum .eq. 12) then
               rsolv(i) = 1.70d0
            else if (atmnum .eq. 14) then
               rsolv(i) = 1.80d0
            else if (atmnum .eq. 15) then
               rsolv(i) = 1.87d0
            else if (atmnum .eq. 16) then
               rsolv(i) = 1.775d0
            else if (atmnum .eq. 17) then
               rsolv(i) = 1.735d0
            else if (atmnum .eq. 18) then
               rsolv(i) = 1.70d0
            else if (atmnum .eq. 19) then
               rsolv(i) = 2.123d0
            else if (atmnum .eq. 20) then
               rsolv(i) = 1.817d0
            else if (atmnum .eq. 35) then
               rsolv(i) = 1.90d0
            else if (atmnum .eq. 36) then
               rsolv(i) = 1.812d0
            else if (atmnum .eq. 37) then
               rsolv(i) = 2.26d0
            else if (atmnum .eq. 53) then
               rsolv(i) = 2.10d0
            else if (atmnum .eq. 54) then
               rsolv(i) = 1.967d0
            else if (atmnum .eq. 55) then
               rsolv(i) = 2.507d0
            else if (atmnum .eq. 56) then
               rsolv(i) = 2.188d0
            end if
         end do
c
c     assign base atomic radii as modified Bondi values
c
      else if (radtyp .eq. 'AMOEBA') then
         do i = 1, n
            atmnum = atomic(i)
            if (atmnum .eq. 0)  rsolv(i) = 0.0d0
            rsolv(i) = vdwrad(atmnum)
            if (atmnum .eq. 1) then
               rsolv(i) = 1.32d0
               k = i12(1,i)
               if (atomic(k) .eq. 7)  rsolv(i) = 1.10d0
               if (atomic(k) .eq. 8)  rsolv(i) = 1.05d0
            end if
            if (atmnum .eq. 3)  rsolv(i) = 1.50d0
            if (atmnum .eq. 6) then
               rsolv(i) = 2.00d0
               if (n12(i) .eq. 3)  rsolv(i) = 2.05d0
               if (n12(i) .eq. 4) then
                  do j = 1, n12(i)
                     k = i12(j,i)
                     if (atomic(k) .eq. 7)  rsolv(i) = 1.75d0
                     if (atomic(k) .eq. 8)  rsolv(i) = 1.75d0
                  end do
               end if
            end if
            if (atmnum .eq. 7) then
               rsolv(i) = 1.60d0
            end if
            if (atmnum .eq. 8) then
               rsolv(i) = 1.55d0
               if (n12(i) .eq. 2)  rsolv(i) = 1.45d0
            end if
         end do
c
c     assign base atomic radii as consensus Bondi values
c
      else
         do i = 1, n
            atmnum = atomic(i)
            if (atmnum .eq. 0)  rsolv(i) = 0.0d0
            rsolv(i) = vdwrad(atmnum)
         end do
      end if
c
c     make Tomasi-style modifications to the base atomic radii
c
      if (radtyp .eq. 'TOMASI') then
         do i = 1, n
            offset = 0.0d0
            atmnum = atomic(i)
            if (atomic(i) .eq. 1) then
               do j = 1, n12(i)
                  k = i12(j,i)
                  if (atomic(k) .eq. 6) then
                     do l = 1, n12(k)
                        m = i12(l,k)
                        if (atomic(m) .eq. 7)  offset = -0.05d0
                        if (atomic(m) .eq. 8)  offset = -0.10d0
                     end do
                  end if
                  if (atomic(k) .eq. 7)  offset = -0.25d0
                  if (atomic(k) .eq. 8)  offset = -0.40d0
                  if (atomic(k) .eq. 16)  offset = -0.10d0
               end do
            else if (atomic(i) .eq. 6) then
               if (n12(i) .eq. 4)  offset = 0.05d0
               if (n12(i) .eq. 3)  offset = 0.02d0
               if (n12(i) .eq. 2)  offset = -0.03d0
               do j = 1, n12(i)
                  k = i12(j,i)
                  if (atomic(k) .eq. 6)  offset = offset - 0.07d0
               end do
               do j = 1, n12(i)
                  k = i12(j,i)
                  if (atomic(k).eq.7 .and. n12(k).eq.4)
     &               offset = -0.20d0
                  if (atomic(k).eq.7 .and. n12(k).eq.3)
     &               offset = -0.25d0
                  if (atomic(k) .eq. 8)  offset = -0.20d0
               end do
            else if (atomic(i) .eq. 7) then
               if (n12(i) .eq. 3) then
                  offset = -0.10d0
                  do j = 1, n12(i)
                     k = i12(j,i)
                     if (atomic(k) .eq. 6)  offset = offset - 0.24d0
                  end do
               else
                  offset = -0.20d0
                  do j = 1, n12(i)
                     k = i12(j,i)
                     if (atomic(k) .eq. 6)  offset = offset - 0.16d0
                  end do
               end if
            else if (atomic(i) .eq. 8) then
               if (n12(i) .eq. 2) then
                  offset = -0.21d0
                  do j = 1, n12(i)
                     k = i12(j,i)
                     if (atomic(k) .eq. 6)  offset = -0.36d0
                  end do
               else
                  offset = -0.25d0
               end if
            else if (atomic(i) .eq. 16) then
               offset = -0.03d0
               do j = 1, n12(i)
                  k = i12(j,i)
                  if (atomic(k) .eq. 6)  offset = offset - 0.10d0
               end do
            end if
            rsolv(i) = rsolv(i) + offset
         end do
      end if
c
c     apply an overall scale factor to the solvation radii
c
      rscale = 1.0d0
      if (radtyp .eq. 'VDW')  rscale = 1.0d0
      if (radtyp .eq. 'MACROMODEL')  rscale = 1.0d0
      if (radtyp .eq. 'AMOEBA')  rscale = 1.0d0
      if (radtyp .eq. 'BONDI')  rscale = 1.03d0
      if (radtyp .eq. 'TOMASI')  rscale = 1.24d0
      do i = 1, n
         rsolv(i) = rsolv(i) * rscale
      end do
c
c     assign generic value for the HCT overlap scale factor
c
      do i = 1, n
         shct(i) = 0.69d0
      end do
      if (radtyp .eq. 'MACROMODEL') then
         do i = 1, n
            shct(i) = 0.80d0
            atmnum = atomic(i)
            if (atmnum .eq. 1)  shct(i) = 0.85d0
            if (atmnum .eq. 6)  shct(i) = 0.72d0
            if (atmnum .eq. 7)  shct(i) = 0.79d0
            if (atmnum .eq. 8)  shct(i) = 0.85d0
            if (atmnum .eq. 9)  shct(i) = 0.88d0
            if (atmnum .eq. 15)  shct(i) = 0.86d0
            if (atmnum .eq. 16)  shct(i) = 0.96d0
            if (atmnum .eq. 26)  shct(i) = 0.88d0
         end do
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kpb  --  assign Poisson-Boltzmann parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kpb" assigns parameters needed for the Poisson-Boltzmann
c     implicit solvation model implemented via APBS
c
c
      subroutine kpb
      use sizes
      use atomid
      use atoms
      use bath
      use couple
      use gkstuf
      use inform
      use iounit
      use keys
      use kvdws
      use math
      use nonpol
      use pbstuf
      use polar
      use potent
      use ptable
      use solute
      implicit none
      integer i,j,k,l,m
      integer nx,ny,nz
      integer maxgrd,next
      integer atmnum,trimtext
      integer pbtyplen,pbsolnlen
      integer bcfllen,chgmlen
      integer srfmlen,pbionq
      real*8 ri,rscale
      real*8 spacing,offset
      real*8 gx,gy,gz
      real*8 xcm,ycm,zcm
      real*8 total,weigh
      real*8 xmin,xmax,ymin
      real*8 ymax,zmin,zmax
      real*8 xlen,ylen,zlen
      real*8 pbionc,pbionr
      character*10 radtyp
      character*20 keyword
      character*20 value
      character*240 record
      character*240 string
c
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(rsolv))  deallocate (rsolv)
      if (allocated(shct))  deallocate (shct)
      if (allocated(udirs))  deallocate (udirs)
      if (allocated(udirps))  deallocate (udirps)
      if (allocated(uinds))  deallocate (uinds)
      if (allocated(uinps))  deallocate (uinps)
      if (allocated(uopts))  deallocate (uopts)
      if (allocated(uoptps))  deallocate (uoptps)
      allocate (rsolv(n))
      allocate (shct(n))
      allocate (udirs(3,n))
      allocate (udirps(3,n))
      allocate (uinds(3,n))
      allocate (uinps(3,n))
      allocate (uopts(0:coptmax,3,n))
      allocate (uoptps(0:coptmax,3,n))
c
c     assign some default APBS configuration parameters
c
      pbtyp = 'LPBE'
      pbsoln = 'MG-MANUAL'
      radtyp = 'BONDI'
      bcfl = 'MDH'
      chgm = 'SPL4'
      srfm = 'SPL4'
      kelvin = 298.0d0
      pdie = 1.0d0
      sdie = 78.3d0
      srad = 1.4d0
      swin = 0.3d0
      sdens = 10.0d0
      smin = 3.0d0
      ionn = 0
      do i = 1, maxion
         ionc(i) = 0.0d0
         ionq(i) = 1
         ionr(i) = 2.0d0
      end do
      spacing = 0.5d0
      maxgrd = 225
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
      gcent(1) = xcm
      gcent(2) = ycm
      gcent(3) = zcm
c
c     set default APBS grid dimension based on system extent
c
      xmin = xcm
      ymin = ycm
      zmin = zcm
      xmax = xcm
      ymax = ycm
      zmax = zcm
      do i = 1, n
         ri = rsolv(i)
         xmin = min(xmin,x(i)-ri)
         ymin = min(ymin,y(i)-ri)
         zmin = min(zmin,z(i)-ri)
         xmax = max(xmax,x(i)+ri)
         ymax = max(ymax,y(i)+ri)
         zmax = max(zmax,z(i)+ri)
      end do
      xlen = 2.0d0 * (max(xcm-xmin,xmax-xcm)+smin)
      ylen = 2.0d0 * (max(ycm-ymin,ymax-ycm)+smin)
      zlen = 2.0d0 * (max(zcm-zmin,zmax-zcm)+smin)
      dime(1) = int(xlen/spacing) + 1
      dime(2) = int(ylen/spacing) + 1
      dime(3) = int(zlen/spacing) + 1
c
c     get any altered APBS parameters from the keyfile
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:8) .eq. 'MG-AUTO ') then
            pbsoln = 'MG-AUTO'
         else if (keyword(1:10) .eq. 'MG-MANUAL ') then
            pbsoln = 'MG-MANUAL'
         else if (keyword(1:10) .eq. 'APBS-GRID ') then
            nx = dime(1)
            ny = dime(2)
            nz = dime(3)
            read (string,*,err=10,end=10)  nx,ny,nz
   10       continue
            if (nx .ge. 33)  dime(1) = nx
            if (ny .ge. 33)  dime(2) = ny
            if (nz .ge. 33)  dime(3) = nz
         else if (keyword(1:10) .eq. 'PB-RADIUS ') then
            call getword (record,value,next)
            call upcase (value)
            if (value(1:3) .eq. 'VDW') then
               radtyp = 'VDW'
            else if (value(1:10) .eq. 'MACROMODEL') then
               radtyp = 'MACROMODEL'
            else if (value(1:5) .eq. 'BONDI') then
               radtyp = 'BONDI'
            else if (value(1:6) .eq. 'TOMASI') then
               radtyp = 'TOMASI'
            end if
         else if (keyword(1:6) .eq. 'SDENS ') then
            read (string,*,err=20,end=20)  sdens
   20       continue
         else if (keyword(1:5) .eq. 'PDIE ') then
            read (string,*,err=30,end=30)  pdie
   30       continue
         else if (keyword(1:5) .eq. 'SDIE ') then
            read (string,*,err=40,end=40)  sdie
   40       continue
         else if (keyword(1:5) .eq. 'SRAD ') then
            read (string,*,err=50,end=50)  srad
   50       continue
         else if (keyword(1:5) .eq. 'SWIN ') then
            read (string,*,err=60,end=60)  swin
   60       continue
         else if (keyword(1:5) .eq. 'SMIN ') then
            read (string,*,err=70,end=70)  smin
   70       continue
         else if (keyword(1:5) .eq. 'SRFM ') then
            call getword (record,value,next)
            call upcase (value)
            if (value(1:3) .eq. 'MOL') then
               srfm = 'MOL'
            else if (value(1:4) .eq. 'SMOL') then
               srfm = 'SMOL'
            else if (value(1:4) .eq. 'SPL2') then
               srfm = 'SPL2'
            end if
         else if (keyword(1:5) .eq. 'BCFL ') then
            call getword (record,value,next)
            call upcase (value)
            if (value(1:3) .eq. 'ZERO') then
               bcfl = 'ZERO'
            else if (value(1:3) .eq. 'MDH') then
               bcfl = 'MDH'
            else if (value(1:3) .eq. 'SDH') then
               bcfl = 'SDH'
            end if
         else if (keyword(1:4) .eq. 'ION ') then
            pbionc = 0.0d0
            pbionq = 1
            pbionr = 2.0d0
            read (string,*,err=80,end=80)  pbionq,pbionc,pbionr
   80       continue
            if (pbionq.ne.0 .and. pbionc.ge.0.0d0
     &             .and. pbionr.ge.0.0d0) then
               ionn = ionn + 1
               ionc(ionn) = pbionc
               ionq(ionn) = pbionq
               ionr(ionn) = pbionr
            end if
         end if
      end do
c
c     set APBS grid spacing for the chosen grid dimension
c
      xlen = 2.0d0 * (max(xcm-xmin,xmax-xcm)+smin)
      ylen = 2.0d0 * (max(ycm-ymin,ymax-ycm)+smin)
      zlen = 2.0d0 * (max(zcm-zmin,zmax-zcm)+smin)
      grid(1) = xlen / dime(1)
      grid(2) = ylen / dime(2)
      grid(3) = zlen / dime(3)
c
c     grid spacing must be equal to maintain traceless quadrupoles
c
      grid(1) = min(grid(1),grid(2),grid(3))
      grid(2) = grid(1)
      grid(3) = grid(1)
c
c     set the grid dimensions to the smallest multiples of 32
c
      dime(1) = 33
      dime(2) = 33
      dime(3) = 33
      do while (grid(1)*dime(1) .lt. xlen)
         dime(1) = dime(1) + 32
      end do
      do while (grid(2)*dime(2) .lt. ylen)
         dime(2) = dime(2) + 32
      end do
      do while (grid(3)*dime(3) .lt. zlen)
         dime(3) = dime(3) + 32
      end do
c
c     limit the grid dimensions and recompute the grid spacing
c
      dime(1) = min(dime(1),maxgrd)
      dime(2) = min(dime(2),maxgrd)
      dime(3) = min(dime(3),maxgrd)
      grid(1) = xlen / dime(1)
      grid(2) = ylen / dime(2)
      grid(3) = zlen / dime(3)
c
c     grid spacing must be equal to maintain traceless quadrupoles
c
      grid(1) = max(grid(1),grid(2),grid(3))
      grid(2) = grid(1)
      grid(3) = grid(1)
c
c     if this is an "mg-auto" (focusing) calculation, set the
c     fine grid to the default size, and the coarse grid to
c     twice its original size. Currently, all energies and
c     forces need to be evaluated at the same resolution
c
      if (pbsoln .eq. 'MG-AUTO') then
         fgrid(1) = grid(1)
         fgrid(2) = grid(2)
         fgrid(3) = grid(3)
         fgcent(1) = gcent(1)
         fgcent(2) = gcent(2)
         fgcent(3) = gcent(3)
         cgrid(1) = 2.0d0 * grid(1)
         cgrid(2) = 2.0d0 * grid(2)
         cgrid(3) = 2.0d0 * grid(3)
      end if
c
c     get any custom APBS grid parameters from the keyfile
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:5) .eq. 'DIME ') then
            read (string,*,err=90,end=90)  nx,ny,nz
            dime(1) = nx
            dime(2) = ny
            dime(3) = nz
   90       continue
            do j = 1, 3
               if (mod(dime(j),32) .ne. 1) then
                  dime(j) = 32*(1+(dime(j)-1)/32) + 1
               end if
            end do
         else if (keyword(1:6) .eq. 'AGRID ') then
            read (string,*,err=100,end=100)  gx,gy,gz
            grid(1) = gx
            grid(2) = gy
            grid(3) = gz
  100       continue
         else if (keyword(1:6) .eq. 'CGRID ') then
            read (string,*,err=110,end=110)  gx,gy,gz
            cgrid(1) = gx
            cgrid(2) = gy
            cgrid(3) = gz
  110       continue
         else if (keyword(1:6) .eq. 'FGRID ') then
            read (string,*,err=120,end=120)  gx,gy,gz
            fgrid(1) = gx
            fgrid(2) = gy
            fgrid(3) = gz
  120       continue
         else if (keyword(1:6) .eq. 'GCENT ') then
            read (string,*,err=130,end=130)  gx,gy,gz
            gcent(1) = gx
            gcent(2) = gy
            gcent(3) = gz
  130       continue
         else if (keyword(1:7) .eq. 'CGCENT ') then
            read (string,*,err=140,end=140)  gx,gy,gz
            cgcent(1) = gx
            cgcent(2) = gy
            cgcent(3) = gz
  140       continue
         else if (keyword(1:7) .eq. 'FGCENT ') then
            read (string,*,err=150,end=150)  gx,gy,gz
            fgcent(1) = gx
            fgcent(2) = gy
            fgcent(3) = gz
  150       continue
         end if
      end do
c
c     assign base atomic radii from consensus vdw values
c
      if (radtyp .eq. 'VDW') then
         do i = 1, n
            rsolv(i) = 2.0d0
            if (class(i) .ne. 0)  rsolv(i) = rad(class(i))
         end do
c
c     assign standard solvation radii adapted from Macromodel
c
      else if (radtyp .eq. 'MACROMODEL') then
         do i = 1, n
            atmnum = atomic(i)
            if (atmnum .eq. 0)  rsolv(i) = 0.0d0
            rsolv(i) = vdwrad(atmnum)
            if (atmnum .eq. 1) then
               rsolv(i) = 1.25d0
               k = i12(1,i)
               if (atomic(k) .eq. 7)  rsolv(i) = 1.15d0
               if (atomic(k) .eq. 8)  rsolv(i) = 1.05d0
            else if (atmnum .eq. 3) then
               rsolv(i) = 1.432d0
            else if (atmnum .eq. 6) then
               rsolv(i) = 1.90d0
               if (n12(i) .eq. 3)  rsolv(i) = 1.875d0
               if (n12(i) .eq. 2)  rsolv(i) = 1.825d0
            else if (atmnum .eq. 7) then
               rsolv(i) = 1.7063d0
               if (n12(i) .eq. 4)  rsolv(i) = 1.625d0
               if (n12(i) .eq. 1)  rsolv(i) = 1.60d0
            else if (atmnum .eq. 8) then
               rsolv(i) = 1.535d0
               if (n12(i) .eq. 1)  rsolv(i) = 1.48d0
            else if (atmnum .eq. 9) then
               rsolv(i) = 1.47d0
            else if (atmnum .eq. 10) then
               rsolv(i) = 1.39d0
            else if (atmnum .eq. 11) then
               rsolv(i) = 1.992d0
            else if (atmnum .eq. 12) then
               rsolv(i) = 1.70d0
            else if (atmnum .eq. 14) then
               rsolv(i) = 1.80d0
            else if (atmnum .eq. 15) then
               rsolv(i) = 1.87d0
            else if (atmnum .eq. 16) then
               rsolv(i) = 1.775d0
            else if (atmnum .eq. 17) then
               rsolv(i) = 1.735d0
            else if (atmnum .eq. 18) then
               rsolv(i) = 1.70d0
            else if (atmnum .eq. 19) then
               rsolv(i) = 2.123d0
            else if (atmnum .eq. 20) then
               rsolv(i) = 1.817d0
            else if (atmnum .eq. 35) then
               rsolv(i) = 1.90d0
            else if (atmnum .eq. 36) then
               rsolv(i) = 1.812d0
            else if (atmnum .eq. 37) then
               rsolv(i) = 2.26d0
            else if (atmnum .eq. 53) then
               rsolv(i) = 2.10d0
            else if (atmnum .eq. 54) then
               rsolv(i) = 1.967d0
            else if (atmnum .eq. 55) then
               rsolv(i) = 2.507d0
            else if (atmnum .eq. 56) then
               rsolv(i) = 2.188d0
            end if
         end do
c
c     assign base atomic radii from consensus vdw values
c
      else
         do i = 1, n
            atmnum = atomic(i)
            if (atmnum .eq. 0)  rsolv(i) = 0.0d0
            rsolv(i) = vdwrad(atmnum)
         end do
      end if
c
c     make Tomasi-style modifications to the base atomic radii
c
      if (radtyp .eq. 'TOMASI') then
         do i = 1, n
            offset = 0.0d0
            atmnum = atomic(i)
            if (atomic(i) .eq. 1) then
               do j = 1, n12(i)
                  k = i12(j,i)
                  if (atomic(k) .eq. 6) then
                     do l = 1, n12(k)
                        m = i12(l,k)
                        if (atomic(m) .eq. 7)  offset = -0.05d0
                        if (atomic(m) .eq. 8)  offset = -0.10d0
                     end do
                  end if
                  if (atomic(k) .eq. 7)  offset = -0.25d0
                  if (atomic(k) .eq. 8)  offset = -0.40d0
                  if (atomic(k) .eq. 16)  offset = -0.10d0
               end do
            else if (atomic(i) .eq. 6) then
               if (n12(i) .eq. 4)  offset = 0.05d0
               if (n12(i) .eq. 3)  offset = 0.02d0
               if (n12(i) .eq. 2)  offset = -0.03d0
               do j = 1, n12(i)
                  k = i12(j,i)
                  if (atomic(k) .eq. 6)  offset = offset - 0.07d0
               end do
               do j = 1, n12(i)
                  k = i12(j,i)
                  if (atomic(k).eq.7 .and. n12(k).eq.4)
     &               offset = -0.20d0
                  if (atomic(k).eq.7 .and. n12(k).eq.3)
     &               offset = -0.25d0
                  if (atomic(k) .eq. 8)  offset = -0.20d0
               end do
            else if (atomic(i) .eq. 7) then
               if (n12(i) .eq. 3) then
                  offset = -0.10d0
                  do j = 1, n12(i)
                     k = i12(j,i)
                     if (atomic(k) .eq. 6)  offset = offset - 0.24d0
                  end do
               else
                  offset = -0.20d0
                  do j = 1, n12(i)
                     k = i12(j,i)
                     if (atomic(k) .eq. 6)  offset = offset - 0.16d0
                  end do
               end if
            else if (atomic(i) .eq. 8) then
               if (n12(i) .eq. 2) then
                  offset = -0.21d0
                  do j = 1, n12(i)
                     k = i12(j,i)
                     if (atomic(k) .eq. 6)  offset = -0.36d0
                  end do
               else
                  offset = -0.25d0
               end if
            else if (atomic(i) .eq. 16) then
               offset = -0.03d0
               do j = 1, n12(i)
                  k = i12(j,i)
                  if (atomic(k) .eq. 6)  offset = offset - 0.10d0
               end do
            end if
            rsolv(i) = rsolv(i) + offset
         end do
      end if
c
c     apply an overall scale factor to the solvation radii
c
      rscale = 1.0d0
c     if (radtyp .eq. 'VDW')  rscale = 1.1d0
      if (radtyp .eq. 'MACROMODEL')  rscale = 1.15d0
      if (radtyp .eq. 'BONDI')  rscale = 1.21d0
      if (radtyp .eq. 'TOMASI')  rscale = 1.47d0
      do i = 1, n
         rsolv(i) = rsolv(i) * rscale
      end do
c
c     assign generic value for the HCT overlap scale factor
c
      do i = 1, n
         shct(i) = 0.69d0
      end do
c
c     determine the length of the character arguments
c
      pbtyplen = trimtext (pbtyp)
      pbsolnlen = trimtext (pbsoln)
      bcfllen = trimtext (bcfl)
      chgmlen = trimtext (chgm)
      srfmlen = trimtext (srfm)
c
c     make call needed to initialize the APBS calculation
c
      call apbsinitial (dime,grid,gcent,cgrid,cgcent,fgrid,fgcent,
     &                  pdie,sdie,srad,swin,sdens,kelvin,ionn,ionc,
     &                  ionq,ionr,pbtyp,pbtyplen,pbsoln,pbsolnlen,
     &                  bcfl,bcfllen,chgm,chgmlen,srfm,srfmlen)
c
c     print out the APBS grid dimensions and spacing
c
      if (verbose) then
         write (iout,160)  (dime(i),i=1,3),grid(1)
  160    format (/,' APBS Grid Dimensions and Spacing :',
     &           //,10x,3i8,10x,f10.4)
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine knp  --  assign cavity-dispersion parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "knp" initializes parameters needed for the cavity-plus-
c     dispersion nonpolar implicit solvation model
c
c
      subroutine knp
      use sizes
      use atomid
      use atoms
      use couple
      use keys
      use kvdws
      use math
      use nonpol
      use potent
      use solute
      implicit none
      integer i,next
      real*8 cavoff,dispoff
      real*8 cross,ah,ao
      real*8 rmini,epsi
      real*8 rmixh,rmixh3
      real*8 rmixh7,emixh
      real*8 rmixo,rmixo3
      real*8 rmixo7,emixo
      real*8 ri,ri3,ri7,ri11
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     set default values for solvent pressure and surface tension
c
      solvprs = 0.0327d0
      surften = 0.080d0
c
c     set default values for cavity and dispersion radius offsets
c
      cavoff = 0.0d0
      dispoff = 0.26d0
c
c     get any altered surface tension value from keyfile
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:17) .eq. 'SOLVENT-PRESSURE ') then
            read (string,*,err=10,end=10)  solvprs
         else if (keyword(1:16) .eq. 'SURFACE-TENSION ') then
            read (string,*,err=10,end=10)  surften
         end if
   10    continue
      end do
c
c     set switching function values for pressure and tension
c
      cross = 3.0d0 * surften / solvprs
      spcut = cross - 3.5d0
      spoff = cross + 3.5d0
      stcut = cross + 3.9d0
      stoff = cross - 3.5d0
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(asolv))  deallocate (asolv)
      if (allocated(rcav))  deallocate (rcav)
      if (allocated(rdisp))  deallocate (rdisp)
      if (allocated(cdisp))  deallocate (cdisp)
      allocate (asolv(n))
      allocate (rcav(n))
      allocate (rdisp(n))
      allocate (cdisp(n))
c
c     assign surface area factors for nonpolar solvation
c
      do i = 1, n
         asolv(i) = surften
      end do
c
c     set cavity and dispersion radii for nonpolar solvation
c
      do i = 1, n
         rcav(i) = rad(class(i)) + cavoff
         rdisp(i) = rad(class(i)) + dispoff
      end do
c
c     compute maximum dispersion energies for each atom
c
      do i = 1, n
         epsi = eps(class(i))
         if (rdisp(i).gt.0.0d0 .and. epsi.gt.0.0d0) then
            rmini = rad(class(i))
            emixo = 4.0d0 * epso * epsi / ((sqrt(epso)+sqrt(epsi))**2)
            rmixo = 2.0d0 * (rmino**3+rmini**3) / (rmino**2+rmini**2)
            rmixo3 = rmixo**3
            rmixo7 = rmixo**7
            ao = emixo * rmixo7
            emixh = 4.0d0 * epsh * epsi / ((sqrt(epsh)+sqrt(epsi))**2)
            rmixh = 2.0d0 * (rminh**3+rmini**3) / (rminh**2+rmini**2)
            rmixh3 = rmixh**3
            rmixh7 = rmixh**7
            ah = emixh * rmixh7
            ri = rdisp(i)
            ri3 = ri**3
            ri7 = ri**7
            ri11 = ri**11
            if (ri .lt. rmixh) then
               cdisp(i) = -4.0d0*pi*emixh*(rmixh3-ri3)/3.0d0
               cdisp(i) = cdisp(i) - emixh*18.0d0/11.0d0*rmixh3*pi
            else
               cdisp(i) = 2.0d0*pi*(2.0d0*rmixh7-11.0d0*ri7)*ah
               cdisp(i) = cdisp(i) / (11.0d0*ri11)
            end if
            cdisp(i) = 2.0d0 * cdisp(i)
            if (ri .lt. rmixo) then
               cdisp(i) = cdisp(i) - 4.0d0*pi*emixo*(rmixo3-ri3)/3.0d0
               cdisp(i) = cdisp(i) - emixo*18.0d0/11.0d0*rmixo3*pi
            else
               cdisp(i) = cdisp(i) + 2.0d0*pi*(2.0d0*rmixo7-11.0d0*ri7)
     &                                  * ao/(11.0d0*ri11)
            end if
         end if
         cdisp(i) = slevy * awater * cdisp(i)
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine khpmf  --  assign hydrophobic PMF parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "khpmf" initializes parameters needed for the hydrophobic
c     potential of mean force nonpolar implicit solvation model
c
c     literature reference:
c
c     M. S. Lin, N. L. Fawzi and T. Head-Gordon, "Hydrophobic
c     Potential of Mean Force as a Solvation Function for Protein
c     Structure Prediction", Structure, 15, 727-740 (2007)
c
c
      subroutine khpmf
      use sizes
      use atomid
      use atoms
      use couple
      use hpmf
      use ptable
      implicit none
      integer i,j,k
      integer nh,atn
      logical keep
c
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(ipmf))  deallocate (ipmf)
      if (allocated(rpmf))  deallocate (rpmf)
      if (allocated(acsa))  deallocate (acsa)
      allocate (ipmf(n))
      allocate (rpmf(n))
      allocate (acsa(n))
c
c     get carbons for PMF and set surface area screening values
c
      npmf = 0
      do i = 1, n
         if (atomic(i) .eq. 6) then
            keep = .true.
            nh = 0
            if (n12(i) .le. 2)  keep = .false.
            do j = 1, n12(i)
               k = i12(j,i)
               if (atomic(k) .eq. 1)  nh = nh + 1
               if (n12(i).eq.3 .and. atomic(k).eq.8)  keep = .false.
            end do
            if (keep) then
               npmf = npmf + 1
               ipmf(npmf) = i
               acsa(i) = 1.0d0
               if (n12(i).eq.3 .and. nh.eq.0)  acsa(i) = 1.554d0
               if (n12(i).eq.3 .and. nh.eq.1)  acsa(i) = 1.073d0
               if (n12(i).eq.4 .and. nh.eq.1)  acsa(i) = 1.276d0
               if (n12(i).eq.4 .and. nh.eq.2)  acsa(i) = 1.045d0
               if (n12(i).eq.4 .and. nh.eq.3)  acsa(i) = 0.880d0
               acsa(i) = acsa(i) * safact/acsurf
            end if
         end if
      end do
c
c     assign HPMF atomic radii from consensus vdw values
c
      do i = 1, n
         rpmf(i) = 1.0d0
         atn = atomic(i)
         if (atn .eq. 0) then 
            rpmf(i) = 0.00d0
         else
            rpmf(i) = vdwrad(atn)
         end if
         if (atn .eq. 5)  rpmf(i) = 1.80d0
         if (atn .eq. 8)  rpmf(i) = 1.50d0
         if (atn .eq. 35)  rpmf(i) = 1.85d0
      end do
      return
      end
