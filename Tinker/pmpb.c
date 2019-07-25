
/*
 * Note: This version of pmpb.c is compatible with APBS Version 1.3
 *
 * Note: The "routines.h" header file includes apbscfg.h, which
 * defines ENABLE_TINKER. ENABLE_TINKER turns on function prototypes
 * in the headers included within "apbs.h". Therefore, "routines.h"
 * needs to be listed first to avoid implicit function definitons.
 */

#include "apbs/apbscfg.h"
#include "apbs/routines.h"
#include "apbs/apbs.h"

/*
 * Note: Use the -Qlowercase flag to deal with Intel
 * Fortran Compiler name mangling on Windows.
 *
 * Underscores are used by default on Linux and MacOS,
 * but not on Windows.
 *
 * We define aliases below to allow TINKER's Fortran code to
 * call the C routines in this file when linking on Windows.
 */

#ifdef _WIN32
#define apbsinitial_ apbsinitial
#define apbsempole_ apbsempole
#define apbsinduce_ apbsinduce
#define apbsnlinduce_ apbsnlinduce
#define pbdirectpolforce_ pbdirectpolforce
#define pbmutualpolforce_ pbmutualpolforce
#define apbsfinal_ apbsfinal
#endif

/***********************************************************************
  Below are some global variables that are saved between APBS calls.
  There may be a better way to do this, although passing pointers
  into FORTRAN should be avoided.
***********************************************************************/

// TINKER MAXATM parameter; must match "sizes.i"
#define maxatm 25000

// TINKER MAXION parameter; must match "pbstuf.i"
#define maxion 10

// APBS configuration objects
Vmem *mem = VNULL;
Vcom *com = VNULL;
Vio *sock = VNULL;
NOsh *nosh = VNULL;

// atom list
Valist *alist[NOSH_MAXMOL];

// potential solutions saved for polarization force calculation
Vgrid *permU[2];
Vgrid *indU[2];
Vgrid *nlIndU[2];

// kappa and dielectric Vgrids (for homogeneous and solvated states)
Vgrid *dielXMap[NOSH_MAXMOL];
Vgrid *dielYMap[NOSH_MAXMOL];
Vgrid *dielZMap[NOSH_MAXMOL];
Vgrid *kappaMap[NOSH_MAXMOL];
Vgrid *potMap[NOSH_MAXMOL];
Vgrid *chargeMap[NOSH_MAXMOL];
double realCenter[3];

/***********************************************************************
  apbsinitial is called from TINKER to:

  (1) Initialize APBS Vcom, Vmem and NOsh objects
  (2) Create a "virtual" APBS input file and parse it
***********************************************************************/

void apbsinitial_(int dime[3], double grid[3], double gcent[3],
                  double cgrid[3], double cgcent[3],
                  double fgrid[3], double fgcent[3],
                  double *pdie, double *sdie,
                  double *srad, double *swin,
                  double *sdens, double *kelvin,
                  int *ionn, double ionc[maxion],
                  int ionq[maxion], double ionr[maxion],
                  char *pbtypef, int *pbtypelen,
                  char *pbsolnf, int *pbsolnlen,
                  char *bcflf, int *bcfllen,
                  char *chgmf, int *chgmlen,
                  char *srfmf, int *srfmlen,
                  int fortranAppendedPbtypeLen,
                  int fortranAppendedPbsolnLen,
                  int fortranAppendedBfclLen,
                  int fortranAppendedChgmLen,
                  int fortranAppendedSrfmLen) {

    /* All character strings passed from FORTRAN result in an integer
       appended to the list of arguments, each equal to the static
       length specified in the TINKER common block 'pb.i' (20).

       Further below the FORTRAN strings will be converted into
       null terminated C-strings.
    */
    char pbtype[21]; // lpbe
    char pbsoln[21]; // mg-manual or mg-auto
    char bcfl[21];   // zero, sdh, mdh
    char chgm[21];   // spl4
    char srfm[21];   // mol, smol, spl2

    /* Bogus argc and argv variables used for Vcom constructor. */
    int argc = 0;
    char **argv;

    /* CPU info */
    int rank, size;

    /* APBS "input file" is a character buffer */
    char  buff[4096];
    char  tmp[1024];

    /* Loop index */
    int i;

    /* Start the timer - although keep in mind it is not stopped until
       apbsfinal is called.  */
    Vnm_tstart(APBS_TIMER_WALL_CLOCK, "APBS WALL CLOCK");

    /* Convert FORTRAN strings to null-terminated C-String */
    strncpy(pbtype, pbtypef, *pbtypelen);
    strncpy(pbsoln, pbsolnf, *pbsolnlen);
    strncpy(bcfl, bcflf, *bcfllen);
    strncpy(chgm, chgmf, *chgmlen);
    strncpy(srfm, srfmf, *srfmlen);
    pbtype[*pbtypelen] = '\0';
    pbsoln[*pbsolnlen] = '\0';
    bcfl[*bcfllen] = '\0';
    chgm[*chgmlen] = '\0';
    srfm[*srfmlen] = '\0';

    /* Rather than require an APBS input file, a character buffer is
       loaded with two ELEC statements:

       (1) The homogeneous calculation.
       (2) The solvated calculation.

       Many options for partial charge systems are not yet supported
       for AMOEBA (or are not appropriate). The subset of ELEC options
       that can be modified are configured using TINKER keywords.

       Initialization of the "nosh" input data structure then proceeds
       using the buffer data. If the syntax of the ELEC statement changes,
       then corresponding changes will be needed to be made below.
     */

    /* Homogeneous */
    strcpy(buff,"ELEC NAME HOMOGENEOUS\n");
    sprintf(tmp,"\t%s\n",pbsoln);
    strcat(buff,tmp);
    sprintf(tmp,"\t%s\n",pbtype);
    strcat(buff,tmp);
    sprintf(tmp,"\tDIME\t%3i %3i %3i\n",dime[0],dime[1],dime[2]);
    strcat(buff,tmp);
    // MG-AUTO
    if (strcmp(pbsoln,"MG-AUTO") == 0) {
       sprintf(tmp,"\tCGLEN  %10.6f %10.6f %10.6f\n",
                      dime[0]*cgrid[0], dime[1]*cgrid[1], dime[2]*cgrid[2]);
       strcat(buff,tmp);
       sprintf(tmp,"\tCGCENT %10.6f %10.6f %10.6f\n",
                      cgcent[0], cgcent[1],cgcent[2]);
       strcat(buff,tmp);
       sprintf(tmp,"\tFGLEN  %10.6f %10.6f %10.6f\n",
                      dime[0]*fgrid[0], dime[1]*fgrid[1], dime[2]*fgrid[2]);
       strcat(buff,tmp);
       sprintf(tmp,"\tFGCENT %10.6f %10.6f %10.6f\n",
                      fgcent[0], fgcent[1], fgcent[2]);
       strcat(buff,tmp);
    } else { // MG-MANUAL
       sprintf(tmp,"\tGLEN  %10.6f %10.6f %10.6f\n",
                      dime[0]*grid[0], dime[1]*grid[1], dime[2]*grid[2]);
       strcat(buff,tmp);
       sprintf(tmp,"\tGCENT %10.6f %10.6f %10.6f\n",
                      gcent[0], gcent[1], gcent[2]);
       strcat(buff,tmp);
    }
    strcat(buff,"\tMOL\t1\n");
    sprintf(tmp,"\tBCFL\t%s\n", bcfl);
    strcat(buff,tmp);
    sprintf(tmp,"\tPDIE  %10.6f\n", *pdie);
    strcat(buff,tmp);
    sprintf(tmp,"\tSDIE  %10.6f\n", *pdie);
    strcat(buff,tmp);
    sprintf(tmp,"\tCHGM\t%s\n", chgm);
    strcat(buff,tmp);
    sprintf(tmp,"\tSRFM\t%s\n", srfm);
    strcat(buff,tmp);
    sprintf(tmp,"\tSRAD  %10.6f\n", *srad);
    strcat(buff,tmp);
    sprintf(tmp,"\tSWIN  %10.6f\n", *swin);
    strcat(buff,tmp);
    sprintf(tmp,"\tSDENS %10.6f\n", *sdens);
    strcat(buff,tmp);
    sprintf(tmp,"\tTEMP  %10.6f\n", *kelvin);
    strcat(buff,tmp);
    strcat(buff,"END\n\n");

    /* Solvated */
    strcat(buff,"ELEC NAME SOLVATED\n");
    sprintf(tmp,"\t%s\n",pbsoln);
    strcat(buff,tmp);
    sprintf(tmp,"\t%s\n",pbtype);
    strcat(buff,tmp);
    sprintf(tmp,"\tDIME\t%3i %3i %3i\n",dime[0],dime[1],dime[2]);
    strcat(buff,tmp);
    // MG-AUTO
    if (strcmp(pbsoln,"MG-AUTO") == 0) {
       sprintf(tmp,"\tCGLEN  %10.6f %10.6f %10.6f\n",
                      dime[0]*cgrid[0], dime[1]*cgrid[1], dime[2]*cgrid[2]);
       strcat(buff,tmp);
       sprintf(tmp,"\tCGCENT %10.6f %10.6f %10.6f\n",
                      cgcent[0], cgcent[1], cgcent[2]);
       strcat(buff,tmp);
       sprintf(tmp,"\tFGLEN  %10.6f %10.6f %10.6f\n",
                      dime[0]*fgrid[0], dime[1]*fgrid[1], dime[2]*fgrid[2]);
       strcat(buff,tmp);
       sprintf(tmp,"\tFGCENT %10.6f %10.6f %10.6f\n",
                      fgcent[0], fgcent[1], fgcent[2]);
       strcat(buff,tmp);
    } else { // MG-MANUAL
       sprintf(tmp,"\tGLEN  %10.6f %10.6f %10.6f\n",
                      dime[0]*grid[0], dime[1]*grid[1], dime[2]*grid[2]);
       strcat(buff,tmp);
       sprintf(tmp,"\tGCENT %10.6f %10.6f %10.6f\n",
                      gcent[0], gcent[1], gcent[2]);
       strcat(buff,tmp);
    }
    strcat(buff,"\tMOL\t1\n");
    sprintf(tmp,"\tBCFL\t%s\n", bcfl);
    strcat(buff,tmp);
    sprintf(tmp,"\tPDIE  %10.6f\n", *pdie);
    strcat(buff,tmp);
    sprintf(tmp,"\tSDIE  %10.6f\n", *sdie);
    strcat(buff,tmp);
    sprintf(tmp,"\tCHGM\t%s\n", chgm);
    strcat(buff,tmp);
    sprintf(tmp,"\tSRFM\t%s\n", srfm);
    strcat(buff,tmp);
    sprintf(tmp,"\tSRAD  %10.6f\n", *srad);
    strcat(buff,tmp);
    sprintf(tmp,"\tSWIN  %10.6f\n", *swin);
    strcat(buff,tmp);
    sprintf(tmp,"\tSDENS %10.6f\n", *sdens);
    strcat(buff,tmp);
    sprintf(tmp,"\tTEMP  %10.6f\n", *kelvin);
    strcat(buff,tmp);
    for (i=0; i < *ionn; i++) {
       sprintf(tmp,"\tION\t%2i %10.6f %10.6f\n", ionq[i], ionc[i], ionr[i]);
       strcat(buff,tmp);
    }
    strcat(buff,"END\n");
    strcat(buff,"\nQUIT\n");

    /* Misc. initializations */
    for (i=0; i<NOSH_MAXMOL; i++) {
       alist[i] = VNULL;
       dielXMap[i] = VNULL;
       dielYMap[i] = VNULL;
       dielZMap[i] = VNULL;
       kappaMap[i] = VNULL;
       potMap[i] = VNULL;
       chargeMap[i] = VNULL;
    }
    for (i=0; i<2; i++) {
       permU[i] = VNULL;
       indU[i] = VNULL;
       nlIndU[i] = VNULL;
    }

    /* Initialization of Vcom, Vmem, and Nosh (via Vio). */
    VASSERT(Vcom_init(&argc, &argv));
    com = Vcom_ctor(1);
    rank = Vcom_rank(com);
    size = Vcom_size(com);
    startVio();
    Vnm_setIoTag(rank, size);
    mem = Vmem_ctor("MAIN");

    /* Print (to io.mc) and then parse the input buffer. */
    Vnm_tprint(0, "\n********* TINKER generated input buffer *********\n\n");
    Vnm_tprint(0, "%s", buff);
    Vnm_tprint(0, "\n*************************************************\n\n");
    nosh = NOsh_ctor(rank, size);
    sock = Vio_ctor("BUFF", "ASC", VNULL, "BUFFER", "r");
    Vio_bufTake(sock, buff, strlen(buff));
    if (!NOsh_parseInput(nosh, sock)) {
       Vnm_tprint( 2, "Error while parsing input file.\n");
       return;
    }
    /* Release the buffer and kill Vio */
    Vio_bufGive(sock);
    Vio_dtor(&sock);
}

/***********************************************************************
  apbsempole is called from TINKER to:

  (1) Solve the PBE using permanent multipoles as the source term
  (2) Save the solution potential for polarization force calculation
  (3) Return permanent electrostatic solvation values for: energy
      field, forces and torques
***********************************************************************/

void apbsempole_(int *natom, double x[maxatm][3],
                 double rad[maxatm], double rpole[maxatm][13],
                 double *total,
                 double energy[maxatm], double fld[maxatm][3],
                 double rff[maxatm][3], double rft[maxatm][3]) {

    /* Misc. pointers to APBS data structures */
    Vpmg  *pmg[NOSH_MAXCALC];
    Vpmgp *pmgp[NOSH_MAXCALC];
    Vpbe  *pbe[NOSH_MAXCALC];
    MGparm *mgparm = VNULL;
    PBEparm *pbeparm = VNULL;
    Vatom *atom = VNULL;

    /* Vgrid configuration for the kappa and dielectric maps */
    double nx,ny,nz,hx,hy,hzed,xmin,ymin,zmin;
    double *data;
    double zkappa2, epsp, epsw;

    /* Loop indeces */
    int i,j;

    /* Observables and unit conversion */
    double sign, force[3], torque[3], field[3];
    double kT,electric,debye;
    double charge, dipole[3], quad[9];
    debye = 4.8033324;

    for (i=0; i<NOSH_MAXCALC; i++) {
       pmg[i] = VNULL;
       pmgp[i] = VNULL;
       pbe[i] = VNULL;
    }

    /* Kill the saved potential Vgrids */
    for (i=0; i<2; i++){
        if (permU[i] != VNULL) Vgrid_dtor(&permU[i]);
        if (indU[i] != VNULL) Vgrid_dtor(&indU[i]);
        if (nlIndU[i] != VNULL) Vgrid_dtor(&nlIndU[i]);
    }

    /* Kill the old atom list */
    if (alist[0] != VNULL) {
       Valist_dtor(&alist[0]);
    }

    /* Create a new atom list (mol == 1) */
    if (alist[0] == VNULL) {
       alist[0] = Valist_ctor();
       alist[0]->atoms = Vmem_malloc(alist[0]->vmem, *natom, (sizeof(Vatom)));
       alist[0]->number = *natom;
    }

    /* Read TINKER input data into Vatom instances. */
    for (i=0; i < alist[0]->number; i++){
       atom = Valist_getAtom(alist[0],i);
       Vatom_setAtomID(atom, i);
       Vatom_setPosition(atom, x[i]);
       Vatom_setRadius(atom, rad[i]);
       charge = rpole[i][0];
       Vatom_setCharge(atom, charge);
       dipole[0] = rpole[i][1];
       dipole[1] = rpole[i][2];
       dipole[2] = rpole[i][3];
       Vatom_setDipole(atom, dipole);
       quad[0] = rpole[i][4];
       quad[1] = rpole[i][5];
       quad[2] = rpole[i][6];
       quad[3] = rpole[i][7];
       quad[4] = rpole[i][8];
       quad[5] = rpole[i][9];
       quad[6] = rpole[i][10];
       quad[7] = rpole[i][11];
       quad[8] = rpole[i][12];
       Vatom_setQuadrupole(atom, quad);
       /* Useful check
       printf(" %i %f (%f,%f,%f)\n",i,rad[i], x[i][0], x[i][1], x[i][2]);
       printf(" %f\n %f,%f,%f\n", charge, dipole[0], dipole[1], dipole[2]);
       printf(" %f\n", quad[0]);
       printf(" %f %f\n", quad[3], quad[4]);
       printf(" %f %f %f\n", quad[6], quad[7], quad[8]); */
       energy[i] = 0.0;
       for (j=0;j<3;j++){
          fld[i][j] = 0.0;
          rff[i][j] = 0.0;
          rft[i][j] = 0.0;
       }
    }

    nosh->nmol = 1;
    Valist_getStatistics(alist[0]);

    /* Only call the setupCalc routine once, so that we can
       reuse this nosh object */
    if (nosh->ncalc < 2) {
       if (NOsh_setupElecCalc(nosh, alist) != 1) {
          printf("Error setting up calculations\n");
          exit(-1);
       }
    }

    /* Solve the LPBE for the homogeneous and then solvated states */
    for (i=0; i<2; i++) {

       /* Useful local variables */
       mgparm = nosh->calc[i]->mgparm;
       pbeparm = nosh->calc[i]->pbeparm;

       /* Just to be robust */
       if (!MGparm_check(mgparm)){
          printf("MGparm Check failed\n");
          printMGPARM(mgparm, realCenter);
          exit(-1);
       }
       if (!PBEparm_check(pbeparm)){
          printf("PBEparm Check failed\n");
          printPBEPARM(pbeparm);
          exit(-1);
       }

       /* Set up the problem */
       mgparm->chgs = VCM_PERMANENT;
       if (!initMG(i, nosh, mgparm, pbeparm, realCenter, pbe,
                   alist, dielXMap, dielYMap, dielZMap,
                   kappaMap, chargeMap, pmgp, pmg, potMap)) {
              Vnm_tprint( 2, "Error setting up MG calculation!\n");
              return;
       }

       /* Solve the PDE */
       if (solveMG(nosh, pmg[i], mgparm->type) != 1) {
           Vnm_tprint(2, "Error solving PDE!\n");
           return;
       }

       /* Set partition information for observables and I/O */
       /* Note - parallel operation has NOT been tested. */
       if (setPartMG(nosh, mgparm, pmg[i]) != 1) {
           Vnm_tprint(2, "Error setting partition info!\n");
           return;
       }

       nx = pmg[i]->pmgp->nx;
       ny = pmg[i]->pmgp->ny;
       nz = pmg[i]->pmgp->nz;
       hx = pmg[i]->pmgp->hx;
       hy = pmg[i]->pmgp->hy;
       hzed = pmg[i]->pmgp->hzed;
       xmin = pmg[i]->pmgp->xmin;
       ymin = pmg[i]->pmgp->ymin;
       zmin = pmg[i]->pmgp->zmin;

       /* Save dielectric/kappa maps into Vgrids, then change the nosh
        * data structure to think it read these maps in from a file.
        * The goal is to save setup time during convergence of the
        * induced dipoles. This is under consideration...
        * */
       /*
       // X (shifted)
       data = Vmem_malloc(mem, nx*ny*nz, sizeof(double));
       Vpmg_fillArray(pmg[i], data, VDT_DIELX, 0.0, pbeparm->pbetype);
       dielXMap[i] = Vgrid_ctor(nx,ny,nz,hx,hy,hzed,
                                xmin + 0.5*hx,ymin,zmin,data);
       dielXMap[i]->readdata = 1;
       // Y (shifted)
       data = Vmem_malloc(mem, nx*ny*nz, sizeof(double));
       Vpmg_fillArray(pmg[i], data, VDT_DIELY, 0.0, pbeparm->pbetype);
       dielYMap[i] = Vgrid_ctor(nx,ny,nz,hx,hy,hzed,
                                      xmin,ymin + 0.5*hy,zmin,data);
       dielYMap[i]->readdata = 1;
       // Z (shifted)
       data = Vmem_malloc(mem, nx*ny*nz, sizeof(double));
       Vpmg_fillArray(pmg[i], data, VDT_DIELZ, 0.0, pbeparm->pbetype);
       dielZMap[i] = Vgrid_ctor(nx,ny,nz,hx,hy,hzed,
                                      xmin,ymin,zmin + 0.5*hzed,data);
       dielZMap[i]->readdata = 1;
       // Kappa
       data = Vmem_malloc(mem, nx*ny*nz, sizeof(double));
       Vpmg_fillArray(pmg[i], data, VDT_KAPPA, 0.0, pbeparm->pbetype);
       kappaMap[i] = Vgrid_ctor(nx,ny,nz,hx,hy,hzed,xmin,ymin,zmin,data);
       kappaMap[i]->readdata = 1;

       // Update the pbeparam structure, since we now have
       // dielectric and kappap maps
       pbeparm->useDielMap = 1;
       pbeparm->dielMapID = i + 1;
       pbeparm->useKappaMap = 1;
       pbeparm->kappaMapID = i + 1;

       */

       data = Vmem_malloc(mem, nx*ny*nz, sizeof(double));
       Vpmg_fillArray(pmg[i], data, VDT_POT, 0.0, pbeparm->pbetype, pbeparm);
       permU[i] = Vgrid_ctor(nx,ny,nz,hx,hy,hzed,xmin,ymin,zmin,data);
       permU[i]->readdata = 1;
       // set readdata flag to have the dtor to free data

       if (i == 0){
          sign = -1.0;
       } else {
          sign = 1.0;
       }

       /* Calculate observables */
       for (j=0; j < alist[0]->number; j++){
         energy[j] += sign * Vpmg_qfPermanentMultipoleEnergy(pmg[i], j);
         Vpmg_fieldSpline4(pmg[i], j, field);
         fld[j][0] += sign * field[0];
         fld[j][1] += sign * field[1];
         fld[j][2] += sign * field[2];
       }

       if (!pmg[i]->pmgp->nonlin &&
          (pmg[i]->surfMeth == VSM_SPLINE ||
           pmg[i]->surfMeth == VSM_SPLINE3 ||
           pmg[i]->surfMeth == VSM_SPLINE4)) {
          for (j=0; j < alist[0]->number; j++){
            Vpmg_qfPermanentMultipoleForce(pmg[i], j, force, torque);
            rff[j][0] += sign * force[0];
            rff[j][1] += sign * force[1];
            rff[j][2] += sign * force[2];
            rft[j][0] += sign * torque[0];
            rft[j][1] += sign * torque[1];
            rft[j][2] += sign * torque[2];
          }
          kT = Vunit_kb * (1e-3) * Vunit_Na * 298.15 * 1.0/4.184;
          epsp = Vpbe_getSoluteDiel(pmg[i]->pbe);
          epsw = Vpbe_getSolventDiel(pmg[i]->pbe);
          if (VABS(epsp-epsw) > VPMGSMALL) {
             for (j=0; j < alist[0]->number; j++){
                Vpmg_dbPermanentMultipoleForce(pmg[i], j, force);
                rff[j][0] += sign * force[0];
                rff[j][1] += sign * force[1];
                rff[j][2] += sign * force[2];
             }
          }
          zkappa2 = Vpbe_getZkappa2(pmg[i]->pbe);
          if (zkappa2 > VPMGSMALL) {
             for (j=0; j < alist[0]->number; j++) {
                Vpmg_ibPermanentMultipoleForce(pmg[i], j, force);
                rff[j][0] += sign * force[0];
                rff[j][1] += sign * force[1];
                rff[j][2] += sign * force[2];
             }
          }
       }
    }

    //nosh->ndiel = 2;
    //nosh->nkappa = 2;
    /*
    printf("Energy (multipole) %f Kcal/mol\n", *energy);
    printf("Energy (volume)    %f Kcal/mol\n", evol * 0.5 * kT);
    */

    // Convert results into kcal/mol units
    kT = Vunit_kb * (1e-3) * Vunit_Na * 298.15 * 1.0/4.184;
    // Electric converts from electron**2/Angstrom to kcal/mol
    electric = 332.063709;
    *total = 0.0;
    for (i=0; i<alist[0]->number; i++){
       /* starting with the field in KT/e/Ang^2 multiply by kcal/mol/KT
          the field is then divided by "electric" to convert to e/Ang^2 */
       energy[i] *= 0.5 * kT;
       *total += energy[i];
       fld[i][0] *= kT / electric;
       fld[i][1] *= kT / electric;
       fld[i][2] *= kT / electric;
       rff[i][0] *= kT;
       rff[i][1] *= kT;
       rff[i][2] *= kT;
       rft[i][0] *= kT;
       rft[i][1] *= kT;
       rft[i][2] *= kT;
    }

    killMG(nosh, pbe, pmgp, pmg);
}

/***********************************************************************
  apbsinduce is called from TINKER during SCRF convergence to:

  (1) Solve the PBE given induced dipoles as the source term
  (2) Save the solution potential for later polarization force calc
  (3) Return the induced reaction field to TINKER
***********************************************************************/

void apbsinduce_(double uind[maxatm][3], double fld[maxatm][3]){

    Vpmg  *pmg[NOSH_MAXCALC];
    Vpmgp *pmgp[NOSH_MAXCALC];
    Vpbe  *pbe[NOSH_MAXCALC];
    MGparm *mgparm = VNULL;
    PBEparm *pbeparm = VNULL;
    Vatom *atom = VNULL;

    /* Observables and unit conversion */
    double field[3];
    double sign,kT,electric;

    /* Potential Vgrid construction */
    double nx,ny,nz,hx,hy,hzed,xmin,ymin,zmin;
    double *data;

    /* Loop variables */
    int i,j;

    VASSERT(nosh != VNULL);
    for (i=0; i<NOSH_MAXCALC; i++) {
       pmg[i] = VNULL;
       pmgp[i] = VNULL;
       pbe[i] = VNULL;
    }

    /* Read TINKER input data into Vatom instances. */
    for (i=0; i < alist[0]->number; i++){
       atom = Valist_getAtom(alist[0],i);
       Vatom_setInducedDipole(atom, uind[i]);
       for (j=0;j<3;j++){
           fld[i][j] = 0.0;
       }
    }

    /* Solve the LPBE for the homogeneous system, then solvated. */
    for (i=0; i<2; i++) {

       pmg[i] = VNULL;
       pmgp[i] = VNULL;
       pbe[i] = VNULL;

       /* Useful local variables */
       mgparm = nosh->calc[i]->mgparm;
       pbeparm = nosh->calc[i]->pbeparm;

       if (!MGparm_check(mgparm)){
          printf("MGparm Check failed\n");
          exit(-1);
       }
       if (!PBEparm_check(pbeparm)){
          printf("PBEparm Check failed\n");
          exit(-1);
       }

       /* Set up problem */
       mgparm->chgs = VCM_INDUCED;
       if (!initMG(i, nosh, mgparm, pbeparm, realCenter, pbe,
                   alist, dielXMap, dielYMap, dielZMap,
                   kappaMap, chargeMap, pmgp, pmg, potMap)) {
           Vnm_tprint( 2, "Error setting up MG calculation!\n");
           return;
       }

       /* Solve the PDE */
       if (solveMG(nosh, pmg[i], mgparm->type) != 1) {
           Vnm_tprint(2, "Error solving PDE!\n");
           return;
       }

       /* Set partition information for observables and I/O */
       if (setPartMG(nosh, mgparm, pmg[i]) != 1) {
           Vnm_tprint(2, "Error setting partition info!\n");
           return;
       }

       /* Save the potential due to local induced dipoles */
       nx = pmg[i]->pmgp->nx;
       ny = pmg[i]->pmgp->ny;
       nz = pmg[i]->pmgp->nz;
       hx = pmg[i]->pmgp->hx;
       hy = pmg[i]->pmgp->hy;
       hzed = pmg[i]->pmgp->hzed;
       xmin = pmg[i]->pmgp->xmin;
       ymin = pmg[i]->pmgp->ymin;
       zmin = pmg[i]->pmgp->zmin;

       if (indU[i] == VNULL) {
          data = Vmem_malloc(mem, nx*ny*nz, sizeof(double));
          Vpmg_fillArray(pmg[i], data, VDT_POT, 0.0, pbeparm->pbetype, pbeparm);
          indU[i] = Vgrid_ctor(nx,ny,nz,hx,hy,hzed,xmin,ymin,zmin,data);
          indU[i]->readdata = 1;
          // set readdata flag to have the dtor to free data
       } else {
          data = indU[i]->data;
          Vpmg_fillArray(pmg[i], data, VDT_POT, 0.0, pbeparm->pbetype, pbeparm);
       }

       if (i == 0){
          sign = -1.0;
       } else {
          sign = 1.0;
       }

       for (j=0; j < alist[0]->number; j++){
          Vpmg_fieldSpline4(pmg[i], j, field);
          fld[j][0] += sign * field[0];
          fld[j][1] += sign * field[1];
          fld[j][2] += sign * field[2];
       }
    }

    /* load results into the return arrays in electron**2/Ang
    /* kT in kcal/mol */
    kT = Vunit_kb * (1e-3) * Vunit_Na * 298.15 / 4.184;
    // electric: conversion from electron**2/Ang to Kcal/mol
    electric = 332.05382;
    for (i=0; i<alist[0]->number; i++){
       // starting with the field in KT/e/Ang^2 multiply by Kcal/mol/KT
       // (conversion to e/Ang^2, which are TINKER field units)
       fld[i][0] *= kT / electric;
       fld[i][1] *= kT / electric;
       fld[i][2] *= kT / electric;
    }

    killMG(nosh, pbe, pmgp, pmg);
}

/***********************************************************************
  apbsnlinduce is called from TINKER during SCRF convergence to:

  (1) Solve the PBE given non-local induced dipoles as the source term
  (2) Save the solution potential for later polarization force calc
  (3) Return the nonlocal induced dipole reaction field to TINKER
************************************************************************/

void apbsnlinduce_(double uinp[maxatm][3], double fld[maxatm][3]){

    /* Misc. pointers to APBS data structures */
    Vpmg  *pmg[NOSH_MAXCALC];
    Vpmgp *pmgp[NOSH_MAXCALC];
    Vpbe  *pbe[NOSH_MAXCALC];
    MGparm *mgparm = VNULL;
    PBEparm *pbeparm = VNULL;
    Vatom *atom = VNULL;

    /* Observables and unit conversion */
    double field[3];
    double sign,kT,electric;
    /* Potential Vgrid construction */
    double nx,ny,nz,hx,hy,hzed,xmin,ymin,zmin;
    double *data;
    /* Loop variables */
    int i,j;

    VASSERT(nosh != VNULL);
    for (i=0; i<NOSH_MAXCALC; i++) {
       pmg[i] = VNULL;
       pmgp[i] = VNULL;
       pbe[i] = VNULL;
    }

    /* Read TINKER induce input data into Vatom instances. */
    for (i=0; i < alist[0]->number; i++){
       atom = Valist_getAtom(alist[0],i);
       Vatom_setNLInducedDipole(atom, uinp[i]);
       for (j=0;j<3;j++){
           fld[i][j] = 0.0;
       }
    }

    /* Solve the LPBE for the homogeneous system, then solvated. */
    for (i=0; i<2; i++) {

       pmg[i] = VNULL;
       pmgp[i] = VNULL;
       pbe[i] = VNULL;

       /* Useful local variables */
       mgparm = nosh->calc[i]->mgparm;
       pbeparm = nosh->calc[i]->pbeparm;

       if (!MGparm_check(mgparm)){
          printf("MGparm Check failed\n");
          exit(-1);
       }
       if (!PBEparm_check(pbeparm)){
          printf("PBEparm Check failed\n");
          exit(-1);
       }

       /* Set up problem */
       mgparm->chgs = VCM_NLINDUCED;
       if (!initMG(i, nosh, mgparm, pbeparm, realCenter, pbe,
                   alist, dielXMap, dielYMap, dielZMap,
                   kappaMap, chargeMap, pmgp, pmg, potMap)) {
           Vnm_tprint( 2, "Error setting up MG calculation!\n");
           return;
       }

       /* Solve the PDE */
       if (solveMG(nosh, pmg[i], mgparm->type) != 1) {
           Vnm_tprint(2, "Error solving PDE!\n");
           return;
       }

       /* Set partition information for observables and I/O */
       if (setPartMG(nosh, mgparm, pmg[i]) != 1) {
           Vnm_tprint(2, "Error setting partition info!\n");
           return;
       }

       /* Save the potential due to non-local induced dipoles */
       nx = pmg[i]->pmgp->nx;
       ny = pmg[i]->pmgp->ny;
       nz = pmg[i]->pmgp->nz;
       hx = pmg[i]->pmgp->hx;
       hy = pmg[i]->pmgp->hy;
       hzed = pmg[i]->pmgp->hzed;
       xmin = pmg[i]->pmgp->xmin;
       ymin = pmg[i]->pmgp->ymin;
       zmin = pmg[i]->pmgp->zmin;

       if (nlIndU[i] == VNULL) {
          data = Vmem_malloc(VNULL, nx*ny*nz, sizeof(double));
          Vpmg_fillArray(pmg[i], data, VDT_POT, 0.0, pbeparm->pbetype, pbeparm);
          nlIndU[i] = Vgrid_ctor(nx,ny,nz,hx,hy,hzed,xmin,ymin,zmin,data);
          nlIndU[i]->readdata = 1; // set readata flag to have dtor free data
       } else {
          data = nlIndU[i]->data;
          Vpmg_fillArray(pmg[i], data, VDT_POT, 0.0, pbeparm->pbetype, pbeparm);
       }

       if (i == 0){
          sign = -1.0;
       } else {
          sign = 1.0;
       }

       for (j=0; j < alist[0]->number; j++){
          Vpmg_fieldSpline4(pmg[i], j, field);
          fld[j][0] += sign * field[0];
          fld[j][1] += sign * field[1];
          fld[j][2] += sign * field[2];
       }
    }

    /* load results into the return arrays in electron**2/Angstrom
    /* kT in kcal/mol */
    kT = Vunit_kb * (1e-3) * Vunit_Na * 298.15 / 4.184;
    // electric: conversion from electron**2/Angstrom to Kcal/mol
    electric = 332.063709;
    for (i=0; i<alist[0]->number; i++){
       fld[i][0] *= kT / electric;
       fld[i][1] *= kT / electric;
       fld[i][2] *= kT / electric;
    }

    killMG(nosh, pbe, pmgp, pmg);
}

/***********************************************************************
  pbdirectpolforce is called from TINKER to:

  (1) compute direct polarization forces and torques using
      saved potentials
***********************************************************************/

void pbdirectpolforce_(double uind[maxatm][3], double uinp[maxatm][3],
                       double rff[maxatm][3], double rft[maxatm][3]) {

    Vpmg  *pmg[NOSH_MAXCALC];
    Vpmgp *pmgp[NOSH_MAXCALC];
    Vpbe  *pbe[NOSH_MAXCALC];
    MGparm *mgparm = VNULL;
    PBEparm *pbeparm = VNULL;
    Vatom *atom = VNULL;
    double kT, force[3], torque[3];
    double sign, zkappa2, epsp, epsw;
    int i,j;

    for (i=0; i<NOSH_MAXCALC; i++) {
       pmg[i] = VNULL;
       pmgp[i] = VNULL;
       pbe[i] = VNULL;
    }

    // Read the converged induced dipole data into APBS Vatom structures.
    for (i=0; i < alist[0]->number; i++){
       atom = Valist_getAtom(alist[0],i);
       Vatom_setInducedDipole(atom, uind[i]);
       Vatom_setNLInducedDipole(atom, uinp[i]);
       for (j=0;j<3;j++){
          rff[i][j] = 0.0;
          rft[i][j] = 0.0;
       }
    }

    for (i=0; i<2; i++) {

       VASSERT(permU[i] != VNULL);
       VASSERT(indU[i] != VNULL);
       VASSERT(nlIndU[i] != VNULL);

       pmg[i] = VNULL;
       pmgp[i] = VNULL;
       pbe[i] = VNULL;

       /* Useful local variables */
       mgparm = nosh->calc[i]->mgparm;
       pbeparm = nosh->calc[i]->pbeparm;

       /* Set up problem */
       if (!initMG(i, nosh, mgparm, pbeparm, realCenter, pbe,
                   alist, dielXMap, dielYMap, dielZMap,
                   kappaMap, chargeMap, pmgp, pmg, potMap)) {
           Vnm_tprint( 2, "Error setting up MG calculation!\n");
           return;
       }

       if (i == 0) {
         sign = -1.0;
       } else {
         sign = 1.0;
       }

       // Q-Phi Force & Torque
       if (!pmg[i]->pmgp->nonlin &&
          (pmg[i]->surfMeth == VSM_SPLINE ||
           pmg[i]->surfMeth == VSM_SPLINE3 ||
           pmg[i]->surfMeth == VSM_SPLINE4)) {
          for (j=0; j < alist[0]->number; j++){
             Vpmg_qfDirectPolForce(pmg[i], permU[i], indU[i], j, force, torque);
             rff[j][0] += sign * force[0];
             rff[j][1] += sign * force[1];
             rff[j][2] += sign * force[2];
             rft[j][0] += sign * torque[0];
             rft[j][1] += sign * torque[1];
             rft[j][2] += sign * torque[2];
             Vpmg_qfNLDirectPolForce(pmg[i], permU[i],
                                     nlIndU[i], j,force,torque);
             rff[j][0] += sign * force[0];
             rff[j][1] += sign * force[1];
             rff[j][2] += sign * force[2];
             rft[j][0] += sign * torque[0];
             rft[j][1] += sign * torque[1];
             rft[j][2] += sign * torque[2];
           }
           // Dieletric Boundary Force
           epsp = Vpbe_getSoluteDiel(pmg[i]->pbe);
           epsw = Vpbe_getSolventDiel(pmg[i]->pbe);
           if (VABS(epsp-epsw) > VPMGSMALL) {
              for (j=0; j < alist[0]->number; j++){
                 Vpmg_dbDirectPolForce(pmg[i], permU[i], indU[i], j, force);
                 rff[j][0] += sign * force[0];
                 rff[j][1] += sign * force[1];
                 rff[j][2] += sign * force[2];
                 Vpmg_dbNLDirectPolForce(pmg[i], permU[i], nlIndU[i], j, force);
                 rff[j][0] += sign * force[0];
                 rff[j][1] += sign * force[1];
                 rff[j][2] += sign * force[2];
              }
           }
           // Ionic Boundary Force
           zkappa2 = Vpbe_getZkappa2(pmg[i]->pbe);
           if (zkappa2 > VPMGSMALL) {
               for (j=0; j < alist[0]->number; j++){
                  Vpmg_ibDirectPolForce(pmg[i], permU[i], indU[i], j, force);
                  rff[j][0] += sign * force[0];
                  rff[j][1] += sign * force[1];
                  rff[j][2] += sign * force[2];
                  Vpmg_ibNLDirectPolForce(pmg[i], permU[i],
                                          nlIndU[i], j, force);
                  rff[j][0] += sign * force[0];
                  rff[j][1] += sign * force[1];
                  rff[j][2] += sign * force[2];
               }
           }
       }
    }

    // kT in kcal/mol
    kT = Vunit_kb * (1e-3) * Vunit_Na * 298.15 / 4.184;
    for (i=0; i<alist[0]->number; i++){
       rff[i][0] *= kT;
       rff[i][1] *= kT;
       rff[i][2] *= kT;
       rft[i][0] *= kT;
       rft[i][1] *= kT;
       rft[i][2] *= kT;
    }

    killMG(nosh, pbe, pmgp, pmg);
}

/***********************************************************************
  pbmutualpolforce is called from TINKER to:

  (1) compute mutual polarization forces using saved potentials
***********************************************************************/

void pbmutualpolforce_(double uind[maxatm][3], double uinp[maxatm][3],
                       double rff[maxatm][3]) {

    Vpmg  *pmg[NOSH_MAXCALC];
    Vpmgp *pmgp[NOSH_MAXCALC];
    Vpbe  *pbe[NOSH_MAXCALC];
    MGparm *mgparm = VNULL;
    PBEparm *pbeparm = VNULL;
    Vatom *atom = VNULL;
    double kT, force[3];
    double sign, zkappa2, epsp, epsw;
    int i,j;

    for (i=0; i<NOSH_MAXCALC; i++) {
       pmg[i] = VNULL;
       pmgp[i] = VNULL;
       pbe[i] = VNULL;
    }

    // Read the converged dipole data into APBS Vatom structures.
    for (i=0; i < alist[0]->number; i++){
       atom = Valist_getAtom(alist[0],i);
       Vatom_setInducedDipole(atom, uind[i]);
       Vatom_setNLInducedDipole(atom, uinp[i]);
       for (j=0;j<3;j++){
          rff[i][j] = 0.0;
       }
    }

    for (i=0; i<2; i++) {

       VASSERT(indU[i] != VNULL);
       VASSERT(nlIndU[i] != VNULL);

       pmg[i] = VNULL;
       pmgp[i] = VNULL;
       pbe[i] = VNULL;

       /* Useful local variables */
       mgparm = nosh->calc[i]->mgparm;
       pbeparm = nosh->calc[i]->pbeparm;

       /* Set up problem */
       if (!initMG(i, nosh, mgparm, pbeparm, realCenter, pbe,
                   alist, dielXMap, dielYMap, dielZMap,
                   kappaMap, chargeMap, pmgp, pmg, potMap)) {
           Vnm_tprint( 2, "Error setting up MG calculation!\n");
           return;
       }

       if (i == 0) {
         sign = -1.0;
       } else {
         sign = 1.0;
       }

       for (j=0; j < alist[0]->number; j++){
          Vpmg_qfMutualPolForce(pmg[i], indU[i], nlIndU[i], j, force);
          rff[j][0] += sign * force[0];
          rff[j][1] += sign * force[1];
          rff[j][2] += sign * force[2];
       }
       epsp = Vpbe_getSoluteDiel(pmg[i]->pbe);
       epsw = Vpbe_getSolventDiel(pmg[i]->pbe);
       if (VABS(epsp-epsw) > VPMGSMALL) {
          for (j=0; j < alist[0]->number; j++){
             Vpmg_dbMutualPolForce(pmg[i], indU[i], nlIndU[i], j, force);
             rff[j][0] += sign * force[0];
             rff[j][1] += sign * force[1];
             rff[j][2] += sign * force[2];
          }
       }
       zkappa2 = Vpbe_getZkappa2(pmg[i]->pbe);
       if (zkappa2 > VPMGSMALL) {
          for (j=0; j < alist[0]->number; j++){
             Vpmg_ibMutualPolForce(pmg[i], indU[i], nlIndU[i], j, force);
             rff[j][0] += sign * force[0];
             rff[j][1] += sign * force[1];
             rff[j][2] += sign * force[2];
          }
       }
    }

    // kT in kcal/mol.
    kT = Vunit_kb * (1e-3) * Vunit_Na * 298.15 / 4.184;
    for (i=0; i<alist[0]->number; i++){
       rff[i][0] *= kT;
       rff[i][1] *= kT;
       rff[i][2] *= kT;
    }

    killMG(nosh, pbe, pmgp, pmg);
}

/***********************************************************************
  apbsfinal is called from TINKER to:

  (1) clean up at the end of an APBS calculation
***********************************************************************/

void apbsfinal_() {
    unsigned long int bytesTotal, highWater;
    int i;

    VASSERT(nosh != VNULL);

     /* Kill the saved potential Vgrids */
    for (i=0; i<2; i++){
      Vgrid_dtor(&permU[i]);
      Vgrid_dtor(&indU[i]);
      Vgrid_dtor(&nlIndU[i]);
    }

    Valist_dtor(&alist[0]);
    /* Saving the kappa and dielectric maps is under consideration.
    killKappaMaps(nosh, kappaMap);
    killDielMaps(nosh, dielXMap, dielYMap, dielZMap);
    */
    NOsh_dtor(&nosh);

    /* Clean up MALOC structures */
    bytesTotal = Vmem_bytesTotal();
    highWater = Vmem_highWaterTotal();

    /*
    printf(" Final APBS memory usage: %4.3f MB total, %4.3f MB high water\n\n",
     (double)(bytesTotal)/(1024.*1024.),(double)(highWater)/(1024.*1024.));
    */

    Vmem_dtor(&mem);
    Vnm_tstop(APBS_TIMER_WALL_CLOCK, "APBS WALL CLOCK");
    Vcom_finalize();
    Vcom_dtor(&com);
}
