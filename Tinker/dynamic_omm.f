c
c
c     ############################################################
c     ##                  COPYRIGHT (C) 2015                    ##
c     ##     by Mark Friedrichs, Lee-Ping Wang & Jay Ponder     ##
c     ##                  All Rights Reserved                   ##
c     ############################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  program dynamic_omm  --  molecular dynamics via OpenMM API  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "dynamic_omm" computes a molecular or stochastic dynamics
c     trajectory via an interface to the OpenMM GPU code for the
c     computation of forces and dynamics integration steps
c
c
      program dynamic_omm
      use sizes
      use atoms
      use bath
      use bndstr
      use bound
      use boxes
      use inform
      use iounit
      use keys
      use mdstuf
      use openmm
      use openmp
      use potent
      use solute
      use stodyn
      use usage
      implicit none
      integer i,istep,nstep
      integer mode,next
      integer nextStep
      integer nextUpdate
      integer updateCalls
      integer callMdStat
      integer callMdSave
      real*8 e,dt
      real*8 dtdump,speed
      real*8 elapsed,cpu
      real*8, allocatable :: derivs(:,:)
      logical exist
      logical updateEachStep
      character*1 answer
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     set up the structure and molecular mechanics calculation
c
      call initial
      call getxyz
      call mechanic
c
c     multipole and polarization are computed together in OpenMM
c
      if (use_mpole .and. .not.use_polar .or.
     &       .not.use_mpole .and. use_polar) then
         use_mpole = .true.
         use_polar = .true.
         call kmpole
         call kpolar
         call mutate
      end if
c
c     initialize the temperature, pressure, integrator and GPU ID
c
      kelvin = 0.0d0
      atmsph = 0.0d0
      isothermal = .false.
      isobaric = .false.
      integrate = 'VERLET'
      cudaDevice = '0               '
c
c     check for keywords containing any altered parameters
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:11) .eq. 'INTEGRATOR ') then
            call getword (record,integrate,next)
            call upcase (integrate)
         else if (keyword(1:12) .eq. 'CUDA-DEVICE ') then
            call gettext (record,cudaDevice,next)
            call upcase (cudaDevice)
         end if
      end do
c
c     initialize the simulation length as number of time steps
c
      nstep = -1
      call nextarg (string,exist)
      if (exist)  read (string,*,err=10,end=10)  nstep
   10 continue
      dowhile (nstep .lt. 0)
         write (iout,20)
   20    format (/,' Enter the Number of Dynamics Steps to be',
     &              ' Taken :  ',$)
         read (input,30,err=40)  nstep
   30    format (i10)
         if (nstep .lt. 0)  nstep = 0
   40    continue
      end do
c
c     get the length of the dynamics time step in picoseconds
c
      dt = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=50,end=50)  dt
   50 continue
      do while (dt .lt. 0.0d0)
         write (iout,60)
   60    format (/,' Enter the Time Step Length in Femtoseconds',
     &              ' [1.0] :  ',$)
         read (input,70,err=80)  dt
   70    format (f20.0)
         if (dt .le. 0.0d0)  dt = 1.0d0
   80    continue
      end do
      dt = 0.001d0 * dt
c
c     enforce bounds on thermostat and barostat coupling times
c
      tautemp = max(tautemp,dt)
      taupres = max(taupres,dt)
c
c     set the time between trajectory snapshot coordinate dumps
c
      dtdump = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=90,end=90)  dtdump
   90 continue
      do while (dtdump .lt. 0.0d0)
         write (iout,100)
  100    format (/,' Enter Time between Dumps in Picoseconds',
     &              ' [0.1] :  ',$)
         read (input,110,err=120)  dtdump
  110    format (f20.0)
         if (dtdump .le. 0.0d0)  dtdump = 0.1d0
  120    continue
      end do
      iwrite = nint(dtdump/dt)
c
c     get choice of statistical ensemble for periodic system
c
      if (use_bounds) then
         mode = -1
         call nextarg (string,exist)
         if (exist)  read (string,*,err=130,end=130)  mode
  130    continue
         do while (mode.lt.1 .or. mode.gt.4)
            write (iout,140)
  140       format (/,' Available Statistical Mechanical Ensembles :',
     &              //,4x,'(1) Microcanonical (NVE)',
     &              /,4x,'(2) Canonical (NVT)',
     &              /,4x,'(3) Isoenthalpic-Isobaric (NPH)',
     &              /,4x,'(4) Isothermal-Isobaric (NPT)',
     &              //,' Enter the Number of the Desired Choice',
     &                 ' [1] :  ',$)
            read (input,150,err=160)  mode
  150       format (i10)
            if (mode .le. 0)  mode = 1
  160       continue
         end do
         if (integrate.eq.'BUSSI' .or. integrate.eq.'NOSE-HOOVER'
     &                .or. integrate.eq.'GHMC') then
            if (mode .ne. 4) then
               mode = 4
               write (iout,170)
  170          format (/,' Switching to NPT Ensemble as Required',
     &                    ' by Chosen Integrator')
            end if
         end if
         if (mode.eq.2 .or. mode.eq.4) then
            isothermal = .true.
            kelvin = -1.0d0
            call nextarg (string,exist)
            if (exist)  read (string,*,err=180,end=180)  kelvin
  180       continue
            do while (kelvin .lt. 0.0d0)
               write (iout,190)
  190          format (/,' Enter the Desired Temperature in Degrees',
     &                    ' K [298] :  ',$)
               read (input,200,err=210)  kelvin
  200          format (f20.0)
               if (kelvin .le. 0.0d0)  kelvin = 298.0d0
  210          continue
            end do
         end if
         if (mode.eq.3 .or. mode.eq.4) then
            isobaric = .true.
            atmsph = -1.0d0
            call nextarg (string,exist)
            if (exist)  read (string,*,err=220,end=220)  atmsph
  220       continue
            do while (atmsph .lt. 0.0d0)
               write (iout,230)
  230          format (/,' Enter the Desired Pressure in Atm',
     &                    ' [1.0] :  ',$)
               read (input,240,err=250)  atmsph
  240          format (f20.0)
               if (atmsph .le. 0.0d0)  atmsph = 1.0d0
  250          continue
            end do
         end if
      end if
c
c     use constant energy or temperature for nonperiodic system
c
      if (.not. use_bounds) then
         mode = -1
         call nextarg (string,exist)
         if (exist)  read (string,*,err=260,end=260)  mode
  260    continue
         do while (mode.lt.1 .or. mode.gt.2)
            write (iout,270)
  270       format (/,' Available Simulation Control Modes :',
     &              //,4x,'(1) Constant Total Energy Value (E)',
     &              /,4x,'(2) Constant Temperature via Thermostat (T)',
     &              //,' Enter the Number of the Desired Choice',
     &                 ' [1] :  ',$)
            read (input,280,err=290)  mode
  280       format (i10)
            if (mode .le. 0)  mode = 1
  290       continue
         end do
         if (mode .eq. 2) then
            isothermal = .true.
            kelvin = -1.0d0
            call nextarg (string,exist)
            if (exist)  read (string,*,err=300,end=300)  kelvin
  300       continue
            do while (kelvin .lt. 0.0d0)
               write (iout,310)
  310          format (/,' Enter the Desired Temperature in Degrees',
     &                    ' K [298] :  ',$)
               read (input,320,err=330)  kelvin
  320          format (f20.0)
               if (kelvin .le. 0.0d0)  kelvin = 298.0d0
  330          continue
            end do
         end if
      end if
c
c     set the frequency with which data is returned from the GPU
c
c     if updateEachStep is true, take single steps on the GPU
c     and retrieve data (positions/velocities/energies) each step;
c     if updateEachStep is false, then take multiple time steps
c     on the GPU before sending data updates to the CPU with the
c     number of steps per update set by the value of "iwrite"
c
      updateEachStep = .false.
c      call nextarg (string,exist)
c      if (exist)  read (string,*,err=340,end=340)  answer
c  340 continue
c      if (.not. exist) then
c         write (iout,350)
c  350    format (/,' Return Data from the GPU at Every Time Step',
c     &              ' [N] :  ',$)
c         read (input,360)  record
c  360    format (a240)
c         next = 1
c      end if
c      call upcase (answer)
c      if (answer .eq. 'Y')  updateEachStep = .true.
c
c     initialize any holonomic constraints and setup dynamics
c
      call shakeup
      call mdinit
c
c     get TINKER energy/gradient values for initial structure
c
      if (verbose) then
         allocate (derivs(3,n))
         call gradient (e,derivs)
         deallocate (derivs)
      end if
c
c     get a list of the available CUDA-capable GPU cards
c
      call set_cuda_devices (cudaDevice)
c
c     map TINKER data structures to OpenMM wrapper structures
c
      call openmm_data ()
c
c     setup the required potential energy terms within OpenMM
c
      call openmm_init (ommHandle,dt)
c
c     compare the energy and gradient between TINKER and OpenMM
c
      if (verbose)  call openmm_test ()
c
c     print out a header line for the dynamics computation
c
      if (integrate .eq. 'VERLET') then
         write (iout,370)
  370    format (/,' Molecular Dynamics Trajectory via',
     &              ' Velocity Verlet Algorithm')
      else if (integrate .eq. 'STOCHASTIC') then
         write (iout,380)
  380    format (/,' Stochastic Dynamics Trajectory via',
     &              ' Velocity Verlet Algorithm')
      else if (integrate .eq. 'BUSSI') then
         write (iout,390)
  390    format (/,' Molecular Dynamics Trajectory via',
     &              ' Bussi-Parrinello NPT Algorithm')
      else if (integrate .eq. 'NOSE-HOOVER') then
         write (iout,400)
  400    format (/,' Molecular Dynamics Trajectory via',
     &              ' Nose-Hoover NPT Algorithm')
      else if (integrate .eq. 'GHMC') then
         write (iout,410)
  410    format (/,' Stochastic Dynamics Trajectory via',
     &              ' Generalized Hybrid Monte Carlo')
      else if (integrate .eq. 'RIGIDBODY') then
         write (iout,420)
  420    format (/,' Molecular Dynamics Trajectory via',
     &              ' Rigid Body Algorithm')
      else if (integrate .eq. 'RESPA') then
         write (iout,430)
  430    format (/,' Molecular Dynamics Trajectory via',
     &              ' r-RESPA MTS Algorithm')
      else
         write (iout,440)
  440    format (/,' Molecular Dynamics Trajectory via',
     &              ' Modified Beeman Algorithm')
      end if
c
c     initialize some counters used during the MD steps
c
      istep = 0
      nextStep = 1
      updateCalls = 0
c
c     integrate equations of motion to take the MD time steps
c
      call settime
      if (updateEachStep) then
         callMdStat = 1
         callMdSave = 1
         do while (istep .lt. nstep)
            call openmm_take_steps (ommHandle,nextStep)
            istep = istep + nextStep
            updateCalls = updateCalls + 1
            call openmm_update (ommHandle,dt,istep,
     &                          callMdStat,callMdSave)
         end do
      else
         nextUpdate = iwrite
         callMdStat = 0
         callMdSave = 1
         do while (istep .lt. nstep)
            nextStep = nextUpdate - istep
            nextUpdate = nextUpdate + iwrite
            if (nextStep+istep .gt. nstep) then
               nextStep = nstep - istep
            end if
            call openmm_take_steps (ommHandle,nextStep)
            istep = istep + nextStep
            updateCalls = updateCalls + 1
            call openmm_update (ommHandle,dt,istep,
     &                          callMdStat,callMdSave)
         end do
      end if
      call gettime (elapsed,cpu)
c
c     print performance and timing information
c
      speed = 0.0d0
      if (elapsed .ne. 0.0d0)  speed = 86.4d0 * nstep * dt / elapsed
      write (iout,450)  speed,elapsed,nstep,updateCalls,
     &                  1000.0d0*dt,n,nthread
  450 format (/,' Performance:  ns/day',9x,f12.4,
     &        /,15x,'Wall Time',6x,f12.4,
     &        /,15x,'Steps',14x,i8,
     &        /,15x,'Updates',12x,i8,
     &        /,15x,'Time Step',6x,f12.4,
     &        /,15x,'Atoms',14x,i8,
     &        /,15x,'Threads',12x,i8)
c
c     perform any final tasks before program exit
c
      call openmm_cleanup (ommHandle)
      call final
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine openmm_data  --  copy TINKER data to OpenMM  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "openmm_data" uses calls to the OpenMM interface to copy data
c     from TINKER modules to the corresponding OpenMM structures
c
c
      subroutine openmm_data ()
      use sizes
      use angbnd
      use angpot
      use angtor
      use atomid
      use atoms
      use bath
      use bitor
      use bndpot
      use bndstr
      use bound
      use boxes
      use cell
      use charge
      use chgpot
      use couple
      use ctran
      use deriv
      use energi
      use ewald
      use freeze
      use group
      use imptor
      use inform
      use ktrtor
      use kvdws
      use limits
      use mdstuf
      use moldyn
      use mplpot
      use molcul
      use mpole
      use mutant
      use nonpol
      use opbend
      use pitors
      use pme
      use polar
      use polgrp
      use polpot
      use potent
      use restrn
      use solute
      use stodyn
      use strbnd
      use strtor
      use torpot
      use tors
      use tortor
      use units
      use urey
      use urypot
      use usage
      use vdw
      use vdwpot
      implicit none
c
c
c     use C++ interface calls to map TINKER variables to OpenMM
c
      call set_angbnd_data (nangle,iang,ak,anat,afld)
      call set_angpot_data (angunit,stbnunit,aaunit,opbunit,opdunit,
     &                      cang,qang,pang,sang,copb,qopb,popb,sopb,
     &                      copd,qopd,popd,sopd,opbtyp,angtyp)
      call set_angtor_data (nangtor,iat,kant)
      call set_atomid_data (tag,class,atomic,valence,mass,name,story)
      call set_atoms_data (n,type,x,y,z)
      call set_bath_data (maxnose,voltrial,kelvin,atmsph,tautemp,
     &                    taupres,compress,collide,eta,volmove,vbar,
     &                    qbar,gbar,vnh,qnh,gnh,isothermal,isobaric,
     &                    anisotrop,volscale,barostat,thermostat)
      call set_bitor_data (nbitor,ibitor)
      call set_bndpot_data (cbnd,qbnd,bndunit,bndtyp)
      call set_bndstr_data (nbond,ibnd,bk,bl)
      call set_bound_data (polycut,polycut2,use_bounds,use_replica,
     &                     use_polymer)
      call set_boxes_data (xbox,ybox,zbox,alpha,beta,gamma,xbox2,
     &                     ybox2,zbox2,box34,volbox,beta_sin,beta_cos,
     &                     gamma_sin,gamma_cos,beta_term,gamma_term,
     &                     lvec,recip,orthogonal,monoclinic,triclinic,
     &                     octahedron,spacegrp)
      call set_cell_data (ncell,icell,xcell,ycell,zcell,
     &                    xcell2,ycell2,zcell2)
      call set_charge_data (nion,iion,jion,kion,pchg)
      call set_chgpot_data (electric,dielec,ebuffer,c2scale,c3scale,
     &                      c4scale,c5scale,neutnbr,neutcut)
      call set_couple_data (n12,i12,n13,i13,n14,i14,n15,i15)
      call set_deriv_data (desum,deb,dea,deba,deub,deaa,deopb,deopd,
     &                     deid,deit,det,dept,debt,deat,dett,dev,dec,
     &                     decd,ded,dem,dep,der,des,delf,deg,dex,dect)
      call set_energi_data (esum,eb,ea,eba,eub,eaa,eopb,eopd,eid,eit,
     &                      et,ept,ebt,eat,ett,ev,ec,ecd,ed,em,ep,er,
     &                      es,elf,eg,ex,ect)
     
      call set_ewald_data (aewald,boundary)
      call set_freeze_data (nrat,nratx,iratx,kratx,irat,rateps,
     &                      krat,use_rattle,ratimage)
      call set_group_data (ngrp,kgrp,grplist,igrp,grpmass,wgrp,
     &                     use_group,use_intra,use_inter)
      call set_imptor_data (nitors,iitors,itors,itors2,itors3)
      call set_inform_data (digits,iprint,iwrite,isend,silent,
     &                      verbose,debug,holdup,abort)
      call set_ktrtor_data (maxntt,maxtgrd,maxtgrd2,tnx,tny,
     &                      ttx,tty,tbf,tbx,tby,tbxy,ktt)
      call set_kvdws_data (rad,eps,rad4,eps4,reduct)
      call set_limits_data (vdwcut,chgcut,dplcut,mpolecut,vdwtaper,
     &                      chgtaper,dpltaper,mpoletaper,ewaldcut,
     &                      usolvcut,use_ewald,use_lights,use_list,
     &                      ctcut,cttaper,use_ctlist, 
     &                      use_vlist,use_clist,use_mlist,use_ulist)

      call set_mdstuf_data (nfree,irest,bmnmix,dorest,velsave,
     &                      frcsave,uindsave,integrate)
      call set_molcul_data (nmol,imol,kmol,molcule,totmass,molmass)
      call set_moldyn_data (v,a,aalt)
      call set_mplpot_data (m2scale,m3scale,m4scale,m5scale)
      call set_mpole_data (maxpole,npole,ipole,polsiz,pollist,
     &                     zaxis,xaxis,yaxis,pole,rpole,polaxe)
      call set_mutant_data (nmut,imut,type0,class0,type1,class1,lambda,
     &                      vlambda,elambda,scexp,scalpha,mut)
      call set_nonpol_data (epso,epsh,rmino,rminh,awater,slevy,
     &                      solvprs,surften,spcut,spoff,stcut,
     &                      stoff,rcav,rdisp,cdisp)
      call set_opbend_data (nopbend,iopb,opbk)
      call set_pitors_data (npitors,ipit,kpit)
      call set_pme_data (nfft1,nfft2,nfft3,bsorder,igrid,bsmod1,
     &                   bsmod2,bsmod3,bsbuild,thetai1,thetai2,
     &                   thetai3,qgrid,qfac)
      call set_polar_data (maxopt,npolar,coptmax,optlevel,copt,copm,
     &                    polarity,penalpha,thole,dirdamp,pdamp,udir,
     &                    udirp,udirs,udirps,uind,uinp,uinds,uinps,
     &                    uopt,uoptp,uopts,uoptps,fopt,foptp,uexact,
     &                    douind)
      call set_polgrp_data (maxp11,maxp12,maxp13,maxp14,np11,
     &                      np12,np13,np14,ip11,ip12,ip13,ip14)
      call set_polpot_data (politer,poleps,p2scale,p3scale,p4scale,
     &                      p5scale,p41scale,d1scale,d2scale,d3scale,
     &                      d4scale,u1scale,u2scale,u3scale,u4scale,
     &                      udiag,poltyp)
      call set_potent_data (use_bond,use_angle,use_strbnd,use_urey,
     &                      use_angang,use_opbend,use_opdist,use_improp,
     &                      use_imptor,use_tors,use_pitors,use_strtor,
     &                      use_angtor,use_tortor,use_vdw,use_charge,
     &                      use_chgdpl,use_dipole,use_mpole,use_polar,
     &                      use_rxnfld,use_solv,use_metal,use_geom,
     &                      use_extra,use_born,use_orbit,use_ct)
      call set_restrn_data (npfix,ndfix,nafix,ntfix,ngfix,nchir,ipfix,
     &                      kpfix,idfix,iafix,itfix,igfix,ichir,depth,
     &                      width,rwall,xpfix,ypfix,zpfix,pfix,dfix,
     &                      afix,tfix,gfix,chir,use_basin,use_wall)
      call set_sizes_data (maxatm,maxtyp,maxclass,maxval,maxref,
     &                     maxgrp,maxres,maxfix)
      call set_solute_data (doffset,p1,p2,p3,p4,p5,rsolv,asolv,rborn,
     &                      drb,drbp,drobc,gpol,shct,aobc,bobc,gobc,
     &                      vsolv,wace,s2ace,uace,solvtyp,borntyp)
      call set_stodyn_data (friction,fgamma,use_sdarea)
      call set_strbnd_data (nstrbnd,isb,sbk)
      call set_strtor_data (nstrtor,ist,kst)
      call set_torpot_data (idihunit,itorunit,torsunit,ptorunit,
     &                      storunit,atorunit,ttorunit)
      call set_tors_data (ntors,itors,tors1,tors2,tors3,
     &                    tors4,tors5,tors6)
      call set_tortor_data (ntortor,itt)
      call set_units_data (avogadro,lightspd,boltzmann,gasconst,emass,
     &                     planck,joule,convert,bohr,hartree,evolt,
     &                     efreq,coulomb,debye,prescon)
      call set_urey_data (nurey,iury,uk,ul)
      call set_urypot_data (cury,qury,ureyunit)
      call set_usage_data (nuse,iuse,use)
      call set_vdw_data (nvdw,ivdw,jvdw,ired,kred,radmin,epsilon,
     &                   radmin4,epsilon4,radhbnd,epshbnd)
      call set_ctran_data (nct,ict,jct,apre,bexp,ct2scale,ct3scale,
     &                     ct4scale,ct5scale,aprerule,bexprule)
      call set_vdwpot_data (maxgauss,ngauss,igauss,abuck,bbuck,cbuck,
     &                      ghal,dhal,v2scale,v3scale,v4scale,v5scale,
     &                      use_vcorr,vdwindex,radtyp,radsiz,gausstyp,
     &                      radrule,epsrule,vdwtyp)
      return
      end
