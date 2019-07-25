c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2006 by Michael Schnieders & Jay W. Ponder  ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  module pbstuf  --  Poisson-Boltzmann solvation parameters  ##
c     ##                                                             ##
c     #################################################################
c
c
c     APBS configuration parameters (see APBS documentation for details)
c     In the column on the right are possible values for each variable,
c     with default values given in brackets. Only a subset of the APBS
c     options are supported and/or are appropriate for use with AMOEBA
c
c     pbtyp                                     lpbe
c
c     At some point AMOEBA with the non-linear PBE could be supported,
c     but this is only worked out for energies (no gradients)
c
c     pbsoln                                    mg-auto, [mg-manual]
c
c     Currently there is only limited support for focusing calculations,
c     which is a powerful feature of APBS. At present, all energies and
c     forces must all be calculated using the finest solution
c
c     bcfl     boundary conditions              zero, sdh, [mdh]
c     chgm     multipole discretization         spl4
c
c     other charge discretization methods are not appropriate for AMOEBA
c
c     srfm     surface method                   mol, smol, [spl4]
c
c     spl4 is required for forces calculations, although mol is useful
c     for comparison with generalized Kirkwood
c
c     dime     number of grid points            [65, 65, 65]
c     grid     grid spacing (mg-manual)         fxn of "dime"
c     cgrid    coarse grid spacing              fxn of "dime"
c     fgrid    fine grid spacing                cgrid / 2
c
c     stable results require grid spacing to be fine enough to keep
c     multipoles inside the dielectric boundary (2.5 * grid < PBR)
c
c     gcent    grid center (mg-manual)          center of mass
c     cgcent   coarse grid center               center of mass
c     fgcent   fine grid center                 center of mass
c     pdie     solute/homogeneous dieletric     [1.0]
c     sdie     solvent dieletric                [78.3]
c     ionn     number of ion species            [0]
c     ionc     ion concentration (M)            [0.0]
c     ionq     ion charge (electrons)           [1.0]
c     ionr     ion radius (A)                   [2.0]
c     srad     solvent probe radius (A)         [1.4]
c     swin     surface spline window width      [0.3]
c     sdens    density of surface points        [10.0]
c
c     additional parameter to facilitate default grid setup
c
c     smin     minimum distance between an      [10.0]
c              atom and the grid boundary (A)
c
c     pbe      Poisson-Boltzmann permanent multipole solvation energy
c     apbe     Poisson-Boltzmann permanent multipole energy over atoms
c     pbr      Poisson-Boltzmann cavity radii for atom types
c     pbep     Poisson-Boltzmann energies on permanent multipoles
c     pbfp     Poisson-Boltzmann forces on permanent multipoles
c     pbtp     Poisson-Boltzmann torques on permanent multipoles
c     pbeuind  Poisson-Boltzmann field due to induced dipoles
c     pbeuinp  Poisson-Boltzmann field due to non-local induced dipoles
c
c
      module pbstuf
      implicit none
      integer maxion
      parameter (maxion=10)
      integer ionn
      integer dime(3)
      integer ionq(maxion)
      real*8 pbe
      real*8 pdie,sdie
      real*8 srad,swin
      real*8 sdens,smin
      real*8 grid(3)
      real*8 gcent(3)
      real*8 cgrid(3)
      real*8 cgcent(3)
      real*8 fgrid(3)
      real*8 fgcent(3)
      real*8 ionr(maxion)
      real*8 ionc(maxion)
      real*8, allocatable :: apbe(:)
      real*8, allocatable :: pbr(:)
      real*8, allocatable :: pbep(:,:)
      real*8, allocatable :: pbfp(:,:)
      real*8, allocatable :: pbtp(:,:)
      real*8, allocatable :: pbeuind(:,:)
      real*8, allocatable :: pbeuinp(:,:)
      character*20 pbtyp,pbsoln
      character*20 bcfl,chgm,srfm
      save
      end
