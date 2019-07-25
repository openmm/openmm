c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine initial  --  initial values and program setup  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "initial" sets up original values for some parameters and
c     variables that might not otherwise get initialized
c
c     note calls below to the "kmp_set" routines are for use with
c     the Intel compiler, but must be commented for other compilers;
c     alternatively, these values can be set via the KMP_STACKSIZE
c     and KMP_BLOCKTIME environment variables
c
c
      subroutine initial
      use sizes
      use align
      use atoms
      use bath
      use bound
      use boxes
      use cell
      use files
      use group
      use inform
      use iounit
      use keys
      use linmin
      use minima
      use molcul
      use mutant
      use neigh
      use openmp
      use output
      use params
      use pdb
      use precis
      use rigid
      use scales
      use sequen
      use socket
      use virial
      use warp
      use zclose
      implicit none
!$    integer omp_get_num_procs
      real*8 precise
      logical first
      save first
      data first  / .true. /
c
c
c     default unit numbers for input and output
c
      input = 5
      iout = 6
c
c     display program banner and copyright notice
c
      if (first)  call promo
c
c     command line arguments to the program
c
      if (first)  call command
      if (first)  first = .false.
c
c     cores, thread count and options for OpenMP
c
      nproc = 1
      nthread = 1
!$    nproc = omp_get_num_procs ()
!$    nthread = nproc
!$    call omp_set_num_threads (nthread)
!$    call omp_set_nested (.true.)
c
c     Intel compiler extensions to OpenMP standard, 268435456 bytes is
c     2**28 bytes, or 256 MB; comment these lines for other compilers
c
c!$   call kmp_set_stacksize_s (268435456)
c!$   call kmp_set_blocktime (0)
c
c     values of machine precision constants
c
      tiny = precise (1)
      small = precise (2)
      huge = precise (3)
c
c     number of lines in the keyfile
c
      nkey = 0
c
c     number of lines in the parameter file
c
      nprm = 0
c
c     number of atoms in the system
c
      n = 0
c
c     number of molecules in the system
c
      nmol = 0
c
c     number of unit cell replicates
c
      ncell = 0
c
c     number of atoms used in superposition
c
      nfit = 0
c
c     number of mutated atoms in the system
c
      nmut = 0
c
c     number of bonds added or deleted from Z-matrix
c
      nadd = 0
      ndel = 0
c
c     number of atoms in Protein Data Bank format
c
      npdb = 0
c
c     number of residues and chains in biopolymer sequence
c
      nseq = 0
      nchain = 0
c
c     highest numbered previous cycle file
c
      nprior = 0
c
c     flags for information levels within the program
c
      silent = .false.
      verbose = .false.
      debug = .false.
      abort = .false.
c
c     flag for use of atom groups
c
      use_group = .false.
c
c     flags for use of periodic boundaries
c
      use_bounds = .false.
      use_replica = .false.
      use_polymer = .false.
c
c     flag for use of internal virial
c
      use_virial = .true.
c
c     default values for unitcell dimensions
c
      xbox = 0.0d0
      ybox = 0.0d0
      zbox = 0.0d0
      alpha = 0.0d0
      beta = 0.0d0
      gamma = 0.0d0
c
c     flags for temperature and pressure baths
c
      isothermal = .false.
      isobaric = .false.
c
c     flags for rebuilding of neighbor lists
c
      dovlst = .true.
      doctlst = .true.
      doclst = .true.
      domlst = .true.
      doulst = .true.
c
c     flag for use of rigid bodies
c
      use_rigid = .false.
c
c     flag to show setting of optimization scale factors
c
      set_scale = .false.
c
c     flags for external Java socket communication
c
      sktstart = .false.
      use_socket = .false.
c
c     flags for potential energy smoothing
c
      use_smooth = .false.
      use_dem = .false.
      use_gda = .false.
      use_tophat = .false.
      use_stophat = .false.
c
c     type of coordinates file
c
      coordtype = 'NONE'
c
c     atomic symbols for elements
c
      call initatom
c
c     names of biopolymer residue types
c
      call initres
c
c     default values used by optimizations
c
      fctmin = 0.0d0
      maxiter = 0
      nextiter = 0
      iprint = -1
      iwrite = -1
      stpmax = 0.0d0
      return
      end
