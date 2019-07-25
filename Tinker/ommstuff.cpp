/*
 *    ################################################################
 *    ##                     COPYRIGHT (C) 2015                     ##
 *    ##  by Mark Friedrichs, Lee-Ping Wang, Kailong Mao, Chao Lu,  ##
 *    ##       Zhi Wang, Matthew Harger and Jay William Ponder      ##
 *    ##                     All Rights Reserved                    ##
 *    ################################################################
 *
 *    ############################################################
 *    ##                                                        ##
 *    ##  ommstuff.cpp  --  TINKER interface to the OpenMM API  ##
 *    ##                                                        ##
 *    ############################################################
 */

/*
 *    ############################################################
 *                     System-Level Include Files
 *    ############################################################
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <fstream>

#define MAX_CUDA_DEVICES 17

using namespace std;

/*
 *    ############################################################
 *           TINKER Routines Called Directly from Interface
 *        (to convert C to C++, must enclose in extern "C" {})
 *    ############################################################
 */

extern "C" {
   void born_ ();
   void bounds_ ();
   void egk1_ ();
   void empole1_ ();
   void enp1_ (double*, double*);
   void ewca1_ (double*);
   void kinetic_ (double*, double(*)[3][3], double*);
   void lattice_ ();
   void mdsave_ (int*, double*, double*, double*);
   void mdstat_ (int*, double*, double*, double*, double*, double*, double*);
}

/*
 *    ############################################################
 *                    OpenMM Wrapper Include Files
 *    ############################################################
 */

#include "OpenMMCWrapper.h"
#include "AmoebaOpenMMCWrapper.h"
#include "OpenMM.h"
using namespace OpenMM;

typedef struct OpenMMData_s OpenMMData;

struct OpenMMData_s {
   OpenMM_System* system;
   OpenMM_Context* context;
   OpenMM_Integrator* integrator;
};

typedef struct {
   char s20[20];
} char20;

/*
 *    ############################################################
 *            C++ Data Structures Corresponding to TINKER
 *    ############################################################
 */

struct {
   int nangle;
   int* iang;
   double* ak;
   double* anat;
   double* afld;
} angbnd__;

struct {
   double angunit;
   double stbnunit;
   double aaunit;
   double opbunit;
   double opdunit;
   double cang;
   double qang;
   double pang;
   double sang;
   double copb;
   double qopb;
   double popb;
   double sopb;
   double copd;
   double qopd;
   double popd;
   double sopd;
   char opbtyp[9];
   char* angtyp;
} angpot__;

static struct {
   int nangtor;
   int* iat;
   double* kant;
} angtor__;

struct {
   int* tag;
   int* classs;   // variable "class" not allowed in C++; add an extra "s"
   int* atomic;
   int* valence;
   double* mass;
   char* name;
   char* story;
} atomid__;

struct {
   int n;
   int* type;
   double* x;
   double* y;
   double* z;
} atoms__;

struct {
   int maxnose;
   int voltrial;
   double kelvin;
   double atmsph;
   double tautemp;
   double taupres;
   double compress;
   double collide;
   double eta;
   double volmove;
   double vbar;
   double qbar;
   double gbar;
   double* vnh;
   double* qnh;
   double* gnh;
   int isothermal;
   int isobaric;
   int anisotrop;
   char volscale[9];
   char barostat[11];
   char thermostat[11];
} bath__;

struct {
   int nbitor;
   int* ibitor;
} bitor__;

struct {
   double cbnd;
   double qbnd;
   double bndunit;
   char bndtyp[9];
} bndpot__;

struct {
   int nbond;
   int* ibnd;
   double* bk;
   double* bl;
} bndstr__;

struct {
   double polycut;
   double polycut2;
   int use_bounds;
   int use_replica;
   int use_polymer;
} bound__;

struct {
   double* xbox;
   double* ybox;
   double* zbox;
   double alpha;
   double beta;
   double gamma;
   double* xbox2;
   double* ybox2;
   double* zbox2;
   double box34;
   double volbox;
   double beta_sin;
   double beta_cos;
   double gamma_sin;
   double gamma_cos;
   double beta_term;
   double gamma_term;
   double lvec[3][3];
   double recip[3][3];
   int orthogonal;
   int monoclinic;
   int triclinic;
   int octahedron;
   char spacegrp[11];
} boxes__;

struct {
   int ncell;
   int* icell;
   double xcell;
   double ycell;
   double zcell;
   double xcell2;
   double ycell2;
   double zcell2;
} cell__;

struct {
   int nion;
   int* iion;
   int* jion;
   int* kion;
   double* pchg;
} charge__;

struct {
   double electric;
   double dielec;
   double ebuffer;
   double c2scale;
   double c3scale;
   double c4scale;
   double c5scale;
   int neutnbr;
   int neutcut;
} chgpot__;

struct {
   int* n12;
   int* i12;
   int* n13;
   int* i13;
   int* n14;
   int* i14;
   int* n15;
   int* i15;
} couple__;

struct {
   int nct;
   int* ict;
   int* jct;
   double* apre;
   double* bexp;
   double ct2scale;
   double ct3scale;
   double ct4scale;
   double ct5scale;
   char aprerule[10];
   char bexprule[10];
} ctran__;

struct {
   double* desum;
   double* deb;
   double* dea;
   double* deba;
   double* deub;
   double* deaa;
   double* deopb;
   double* deopd;
   double* deid;
   double* deit;
   double* det;
   double* dept;
   double* debt;
   double* deat;
   double* dett;
   double* dev;
   double* dec;
   double* decd;
   double* ded;
   double* dem;
   double* dep;
   double* der;
   double* des;
   double* delf;
   double* deg;
   double* dex;
   double* dect;
} deriv__;

struct {
   double* esum;
   double* eb;
   double* ea;
   double* eba;
   double* eub;
   double* eaa;
   double* eopb;
   double* eopd;
   double* eid;
   double* eit;
   double* et;
   double* ept;
   double* ebt;
   double* eat;
   double* ett;
   double* ev;
   double* ec;
   double* ecd;
   double* ed;
   double* em;
   double* ep;
   double* er;
   double* es;
   double* elf;
   double* eg;
   double* ex;
   double* ect;
} energi__;

struct {
   double aewald;
   char boundary[7];
} ewald__;

struct {
   int nrat;
   int nratx;
   int* iratx;
   int* kratx;
   int* irat;
   double rateps;
   double* krat;
   int use_rattle;
   int* ratimage;
} freeze__;

struct {
   int ngrp;
   int* kgrp;
   int* grplist;
   int* igrp;
   double* grpmass;
   double* wgrp;
   int use_group;
   int use_intra;
   int use_inter;
} group__;

struct {
   int nitors;
   int* iitors;
   double* itors1;
   double* itors2;
   double* itors3;
} imptor__;

struct {
   int digits;
   int iprint;
   int iwrite;
   int isend;
   int silent;
   int verbose;
   int debug;
   int holdup;
   int abort;
} inform__;

struct {
   int maxntt;
   int maxtgrd;
   int maxtgrd2;
   int* tnx;
   int* tny;
   double* ttx;
   double* tty;
   double* tbf;
   double* tbx;
   double* tby;
   double* tbxy;
   char20* ktt;
} ktrtor__;

struct {
   int maxnvp;
   double* radpr;
   double* epspr;
   char* kvpr;
} kvdwpr__;

struct {
   double* rad;
   double* eps;
   double* rad4;
   double* eps4;
   double* reduct;
} kvdws__;

struct {
   double vdwcut;
   double ctcut;
   double chgcut;
   double dplcut;
   double mpolecut;
   double vdwtaper;
   double cttaper;
   double chgtaper;
   double dpltaper;
   double mpoletaper;
   double ewaldcut;
   double usolvcut;
   int use_ewald;
   int use_lights;
   int use_list;
   int use_vlist;
   int use_ctlist;
   int use_clist;
   int use_mlist;
   int use_ulist;
} limits__;

struct {
   int nfree;
   int irest;
   int bmnmix;
   int dorest;
   int velsave;
   int frcsave;
   int uindsave;
   char integrate[11];
} mdstuf__;

struct {
   double* v;
   double* a;
   double* aalt;
} moldyn__;

struct {
   int nmol;
   int* imol;
   int* kmol;
   int* molcule;
   double totmass;
   double* molmass;
} molcul__;

struct {
   double m2scale;
   double m3scale;
   double m4scale;
   double m5scale;
} mplpot__;

struct {
   int maxpole;
   int npole;
   int* ipole;
   int* polsiz;
   int* pollist;
   int* zaxis;
   int* xaxis;
   int* yaxis;
   double* pole;
   double* rpole;
   char* polaxe;
} mpole__;

struct {
   int nmut;
   int* imut;
   int* type0;
   int* class0;
   int* type1;
   int* class1;
   double lambda;
   double vlambda;
   double elambda;
   double scexp;
   double scalpha;
   int* mut;
} mutant__;

struct {
   double epso;
   double epsh;
   double rmino;
   double rminh;
   double awater;
   double slevy;
   double solvprs;
   double surften;
   double spcut;
   double spoff;
   double stcut;
   double stoff;
   double* rcav;
   double* rdisp;
   double* cdisp;
} nonpol__;

struct {
   int nopbend;
   int* iopb;
   double* opbk;
} opbend__;

struct {
   int npitors;
   int* ipit;
   double* kpit;
} pitors__;

struct {
   int nfft1;
   int nfft2;
   int nfft3;
   int bsorder;
   int* igrid;
   double* bsmod1;
   double* bsmod2;
   double* bsmod3;
   double* bsbuild;
   double*** thetai1;
   double*** thetai2;
   double*** thetai3;
   double**** qgrid;
   double*** qfac;
} pme__;

struct {
   int maxopt;
   int npolar;
   int coptmax;
   int optlevel;
   double* copt;
   double* copm;
   double* polarity;
   double* penalpha;
   double* thole;
   double* dirdamp;
   double* pdamp;
   double* udir;
   double* udirp;
   double* udirs;
   double* udirps;
   double* uind;
   double* uinp;
   double* uinds;
   double* uinps;
   double* uopt;
   double* uoptp;
   double* uopts;
   double* uoptps;
   double* fopt;
   double* foptp;
   double* uexact;
   int* douind;
} polar__;

struct {
   int maxp11;
   int maxp12;
   int maxp13;
   int maxp14;
   int* np11;
   int* np12;
   int* np13;
   int* np14;
   int* ip11;
   int* ip12;
   int* ip13;
   int* ip14;
} polgrp__;

struct {
   int politer;
   double poleps;
   double p2scale;
   double p3scale;
   double p4scale;
   double p5scale;
   double p41scale;
   double d1scale;
   double d2scale;
   double d3scale;
   double d4scale;
   double u1scale;
   double u2scale;
   double u3scale;
   double u4scale;
   double udiag;
   char poltyp[6];
} polpot__;

struct {
   int use_bond;
   int use_angle;
   int use_strbnd;
   int use_urey;
   int use_angang;
   int use_opbend;
   int use_opdist;
   int use_improp;
   int use_imptor;
   int use_tors;
   int use_pitors;
   int use_strtor;
   int use_angtor;
   int use_tortor;
   int use_vdw;
   int use_charge;
   int use_chgdpl;
   int use_dipole;
   int use_mpole;
   int use_polar;
   int use_rxnfld;
   int use_solv;
   int use_metal;
   int use_geom;
   int use_extra;
   int use_born;
   int use_orbit;
   int use_ct;
} potent__;

struct {
   int npfix;
   int ndfix;
   int nafix;
   int ntfix;
   int ngfix;
   int nchir;
   int* ipfix;
   int* kpfix;
   int* idfix;
   int* iafix;
   int* itfix;
   int* igfix;
   int* ichir;
   double depth;
   double width;
   double rwall;
   double* xpfix;
   double* ypfix;
   double* zpfix;
   double* pfix;
   double* dfix;
   double* afix;
   double* tfix;
   double* gfix;
   double* chir;
   int use_basin;
   int use_wall;
} restrn__;

struct {
   int maxatm;
   int maxtyp;
   int maxclass;
   int maxval;
   int maxref;
   int maxgrp;
   int maxres;
   int maxfix;
} sizes__;

struct {
   double doffset;
   double p1;
   double p2;
   double p3;
   double p4;
   double p5;
   double* rsolv;
   double* asolv;
   double* rborn;
   double* drb;
   double* drbp;
   double* drobc;
   double* gpol;
   double* shct;
   double* aobc;
   double* bobc;
   double* gobc;
   double* vsolv;
   double* wace;
   double* s2ace;
   double* uace;
   char solvtyp[8];
   char borntyp[8];
} solute__;

struct {
   double friction;
   double* fgamma;
   int use_sdarea;
} stodyn__;

struct {
   int nstrbnd;
   int* isb;
   double* sbk;
} strbnd__;

static struct {
   int nstrtor;
   int* ist;
   double* kst;
} strtor__;

struct {
   double idihunit;
   double itorunit;
   double torsunit;
   double ptorunit;
   double storunit;
   double atorunit;
   double ttorunit;
} torpot__;

struct {
   int ntors;
   int* itors;
   double* tors1;
   double* tors2;
   double* tors3;
   double* tors4;
   double* tors5;
   double* tors6;
} tors__;

struct {
   int ntortor;
   int* itt;
} tortor__;

struct {
   double avogadro;
   double lightspd;
   double boltzmann;
   double gasconst;
   double emass;
   double planck;
   double joule;
   double convert;
   double bohr;
   double hartree;
   double evolt;
   double efreq;
   double coulomb;
   double debye;
   double prescon;
} units__;

struct {
   int nurey;
   int* iury;
   double* uk;
   double* ul;
} urey__;

struct {
   double cury;
   double qury;
   double ureyunit;
} urypot__;

struct {
   int* nuse;
   int* iuse;
   int* use;
} usage__;

struct {
   int nvdw;
   int* ivdw;
   int* jvdw;
   int* ired;
   double* kred;
   double* radmin;
   double* epsilon;
   double* radmin4;
   double* epsilon4;
   double* radhbnd;
   double* epshbnd;
} vdw__;


struct {
   int maxgauss;
   int ngauss;
   double* igauss;
   double abuck;
   double bbuck;
   double cbuck;
   double ghal;
   double dhal;
   double v2scale;
   double v3scale;
   double v4scale;
   double v5scale;
   int use_vcorr;
   char vdwindex[5];
   char radtyp[5];
   char radsiz[8];
   char gausstyp[8];
   char radrule[10];
   char epsrule[10];
   char vdwtyp[13];
} vdwpot__;

static void setNullTerminator (char* string, int maxLength, char* buffer) {

   int count;
   int ptr_i;
   char* ptr_c;
   ptr_c = string;
   ptr_i = (int) (*ptr_c);
   count = 0;

   // add contents of string to buffer until a non-character is found
   // or end of the string is reached, then add NULL to end of buffer

   while (ptr_i > 33 && ptr_i < 126 && count < maxLength) {
      buffer[count++] = (*ptr_c);
      ptr_c++;
      ptr_i = (int) (*ptr_c);
   }

   buffer[count] = '\0';

   return;
}

/*
 *    ############################################################
 *            Copy TINKER Fortran Data into C++ Structures
 *    ############################################################
 */

extern "C" {

void set_angbnd_data_ (int* nangle, int* iang, double* ak,
                       double* anat, double* afld) {

   angbnd__.nangle = *nangle;
   angbnd__.iang = iang;
   angbnd__.ak = ak;
   angbnd__.anat = anat;
   angbnd__.afld = afld;
}

void set_angpot_data_ (double* angunit, double* stbnunit, double* aaunit,
                       double* opbunit, double* opdunit, double* cang,
                       double* qang, double* pang, double* sang, double* copb,
                       double* qopb, double* popb, double* sopb, double* copd,
                       double* qopd, double* popd, double* sopd, char* opbtyp,
                       char* angtyp) {

   angpot__.angunit = *angunit;
   angpot__.stbnunit = *stbnunit;
   angpot__.aaunit = *aaunit;
   angpot__.opbunit = *opbunit;
   angpot__.opdunit = *opdunit;
   angpot__.cang = *cang;
   angpot__.qang = *qang;
   angpot__.pang = *pang;
   angpot__.sang = *sang;
   angpot__.copb = *copb;
   angpot__.qopb = *qopb;
   angpot__.popb = *popb;
   angpot__.sopb = *sopb;
   angpot__.copd = *copd;
   angpot__.qopd = *qopd;
   angpot__.popd = *popd;
   angpot__.sopd = *sopd;
   setNullTerminator (opbtyp, 8, angpot__.opbtyp);
   angpot__.angtyp = angtyp;
}

void set_angtor_data_ (int* nangtor, int* iat, double* kant) {

   angtor__.nangtor = *nangtor;
   angtor__.iat = iat;
   angtor__.kant = kant;
}

void set_atomid_data_ (int* tag, int* classs, int* atomic, int* valence,
                       double* mass, char* name, char* story) {

   atomid__.tag = tag;
   atomid__.classs = classs;
   atomid__.atomic = atomic;
   atomid__.valence = valence;
   atomid__.mass = mass;
   atomid__.name = name;
   atomid__.story = story;
}

void set_atoms_data_ (int* n, int* type, double* x, double* y, double* z) {

   int ii;
   atoms__.n = *n;
   atoms__.type = type;
   atoms__.x = x;
   atoms__.y = y;
   atoms__.z = z;
}

void set_bath_data_ (int* maxnose, int* voltrial, double* kelvin,
                     double* atmsph, double* tautemp, double* taupres,
                     double* compress, double* collide, double* eta,
                     double* volmove, double* vbar, double* qbar,
                     double* gbar, double* vnh, double* qnh, double* gnh,
                     int* isothermal, int* isobaric, int* anisotrop,
                     char* volscale, char* barostat, char* thermostat) {

   bath__.maxnose = *maxnose;
   bath__.voltrial = *voltrial;
   bath__.kelvin = *kelvin;
   bath__.atmsph = *atmsph;
   bath__.tautemp = *tautemp;
   bath__.taupres = *taupres;
   bath__.compress = *compress;
   bath__.collide = *collide;
   bath__.eta = *eta;
   bath__.volmove = *volmove;
   bath__.vbar = *vbar;
   bath__.qbar = *qbar;
   bath__.gbar = *gbar;
   bath__.vnh = vnh;
   bath__.qnh = qnh;
   bath__.gnh = gnh;
   bath__.isothermal = *isothermal;
   bath__.isobaric = *isobaric;
   bath__.anisotrop = *anisotrop;
   setNullTerminator (volscale, 9, bath__.volscale);
   setNullTerminator (barostat, 11, bath__.barostat);
   setNullTerminator (thermostat, 11, bath__.thermostat);
}

void set_bitor_data_ (int* nbitor, int* ibitor) {

   bitor__.nbitor = *nbitor;
   bitor__.ibitor = ibitor;
}

void set_bndpot_data_ (double* cbnd, double* qbnd, double* bndunit,
                       char* bndtyp) {

   bndpot__.cbnd = *cbnd;
   bndpot__.qbnd = *qbnd;
   bndpot__.bndunit = *bndunit;
   setNullTerminator (bndtyp, 8, bndpot__.bndtyp);
}

void set_bndstr_data_ (int* nbond, int* ibnd, double* bk, double* bl) {

   bndstr__.nbond = *nbond;
   bndstr__.ibnd = ibnd;
   bndstr__.bk = bk;
   bndstr__.bl = bl;
}

void set_bound_data_ (double* polycut, double* polycut2, int* use_bounds,
                      int* use_replica, int* use_polymer) {

   bound__.polycut = *polycut;
   bound__.polycut2 = *polycut2;
   bound__.use_bounds = *use_bounds;
   bound__.use_replica = *use_replica;
   bound__.use_polymer = *use_polymer;
}

void set_cell_data_ (int* ncell, int* icell, double* xcell, double* ycell,
                     double* zcell, double* xcell2, double* ycell2,
                     double* zcell2) {
   cell__.ncell = *ncell;
   cell__.icell = icell;
   cell__.xcell = *xcell;
   cell__.ycell = *ycell;
   cell__.zcell = *zcell;
   cell__.xcell2 = *xcell2;
   cell__.ycell2 = *ycell2;
   cell__.zcell2 = *zcell2;
}

void set_boxes_data_ (double* xbox, double* ybox, double* zbox,
                      double* alpha, double* beta, double* gamma,
                      double* xbox2, double* ybox2, double* zbox2,
                      double* box34, double* volbox, double* beta_sin,
                      double* beta_cos, double* gamma_sin, double* gamma_cos,
                      double* beta_term, double* gamma_term, double* lvec,
                      double* recip, int* orthogonal, int* monoclinic,
                      int* triclinic, int* octahedron, char* spacegrp) {

   int ii, jj;
   double* ptr_d;

   boxes__.xbox = xbox;
   boxes__.ybox = ybox;
   boxes__.zbox = zbox;
   boxes__.alpha = *alpha;
   boxes__.beta = *beta;
   boxes__.gamma = *gamma;
   boxes__.xbox2 = xbox2;
   boxes__.ybox2 = ybox2;
   boxes__.zbox2 = zbox2;
   boxes__.box34 = *box34;
   boxes__.volbox = *volbox;
   boxes__.beta_sin = *beta_sin;
   boxes__.beta_cos = *beta_cos;
   boxes__.gamma_sin = *gamma_sin;
   boxes__.gamma_cos = *gamma_cos;
   boxes__.beta_term = *beta_term;
   boxes__.gamma_term = *gamma_term;

   ptr_d = lvec;
   for (ii = 0; ii < 3; ii++) {
      for (jj = 0; jj < 3; jj++) {
         boxes__.lvec[ii][jj] = *ptr_d;
         ptr_d++;
      }
   }
   ptr_d = recip;
   for (ii = 0; ii < 3; ii++) {
      for (jj = 0; jj < 3; jj++) {
         boxes__.recip[ii][jj] = *ptr_d;
         ptr_d++;
      }
   }

   boxes__.orthogonal = *orthogonal;
   boxes__.monoclinic = *monoclinic;
   boxes__.triclinic = *triclinic;
   boxes__.octahedron = *octahedron;
   setNullTerminator (spacegrp, 10, boxes__.spacegrp);
}

void set_charge_data_ (int* nion, int* iion, int* jion,
                       int* kion, double* pchg) {

   charge__.nion = *nion;
   charge__.iion = iion;
   charge__.jion = jion;
   charge__.kion = kion;
   charge__.pchg = pchg;
}


void set_chgpot_data_ (double* electric, double* dielec, double* ebuffer,
                       double* c2scale, double* c3scale, double* c4scale,
                       double* c5scale, int* neutnbr, int* neutcut) {

   chgpot__.electric = *electric;
   chgpot__.dielec = *dielec;
   chgpot__.ebuffer = *ebuffer;
   chgpot__.c2scale = *c2scale;
   chgpot__.c3scale = *c3scale;
   chgpot__.c4scale = *c4scale;
   chgpot__.c5scale = *c5scale;
   chgpot__.neutnbr = *neutnbr;
   chgpot__.neutcut = *neutcut;
}

void set_couple_data_ (int* n12, int* i12, int* n13, int* i13,
                       int* n14, int* i14, int* n15, int* i15) {

   couple__.n12 = n12;
   couple__.i12 = i12;
   couple__.n13 = n13;
   couple__.i13 = i13;
   couple__.n14 = n14;
   couple__.i14 = i14;
   couple__.n15 = n15;
   couple__.i15 = i15;
}

void set_deriv_data_ (double* desum, double* deb, double* dea,
                      double* deba, double* deub, double* deaa,
                      double* deopb, double* deopd, double* deid,
                      double* deit, double* det, double* dept,
                      double* debt, double* deat, double* dett,
                      double* dev, double* dec, double* decd,
                      double* ded, double* dem, double* dep,
                      double* der, double* des, double* delf,
                      double* deg, double* dex, double* dect) {

   deriv__.desum = desum;
   deriv__.deb = deb;
   deriv__.dea = dea;
   deriv__.deba = deba;
   deriv__.deub = deub;
   deriv__.deaa = deaa;
   deriv__.deopb = deopb;
   deriv__.deopd = deopd;
   deriv__.deid = deid;
   deriv__.deit = deit;
   deriv__.det = det;
   deriv__.dept = dept;
   deriv__.debt = debt;
   deriv__.deat = deat;
   deriv__.dett = dett;
   deriv__.dev = dev;
   deriv__.dec = dec;
   deriv__.decd = decd;
   deriv__.ded = ded;
   deriv__.dem = dem;
   deriv__.dep = dep;
   deriv__.der = der;
   deriv__.des = des;
   deriv__.delf = delf;
   deriv__.deg = deg;
   deriv__.dex = dex;
   deriv__.dect = dect;
}

void set_energi_data_ (double* esum, double* eb, double* ea,
                       double* eba, double* eub, double* eaa,
                       double* eopb, double* eopd, double* eid,
                       double* eit, double* et, double* ept,
                       double* ebt, double* eat, double* ett,
                       double* ev, double* ec, double* ecd,
                       double* ed, double* em, double* ep,
                       double* er, double* es, double* elf,
                       double* eg, double* ex, double* ect) {

   energi__.esum = esum;
   energi__.eb = eb;
   energi__.ea = ea;
   energi__.eba = eba;
   energi__.eub = eub;
   energi__.eaa = eaa;
   energi__.eopb = eopb;
   energi__.eopd = eopd;
   energi__.eid = eid;
   energi__.eit = eit;
   energi__.et = et;
   energi__.ept = ept;
   energi__.ebt = ebt;
   energi__.eat = eat;
   energi__.ett = ett;
   energi__.ev = ev;
   energi__.ec = ec;
   energi__.ecd = ecd;
   energi__.ed = ed;
   energi__.em = em;
   energi__.ep = ep;
   energi__.er = er;
   energi__.es = es;
   energi__.elf = elf;
   energi__.eg = eg;
   energi__.ex = ex;
   energi__.ect = ect;
}

void set_ewald_data_ (double* aewald, char* boundary) {

   ewald__.aewald = *aewald;
   setNullTerminator (boundary, 7, ewald__.boundary);
}

void set_freeze_data_ (int* nrat, int* nratx, int* iratx, int* kratx,
                       int* irat, double* rateps, double* krat,
                       int* use_rattle, int* ratimage) {

   freeze__.nrat = *nrat;
   freeze__.nratx = *nratx;
   freeze__.iratx = iratx;
   freeze__.kratx = kratx;
   freeze__.irat = irat;
   freeze__.rateps = *rateps;
   freeze__.krat = krat;
   freeze__.use_rattle = *use_rattle;
   freeze__.ratimage = ratimage;
}

void set_group_data_ (int* ngrp, int* kgrp, int* grplist, int* igrp,
                      double* grpmass, double* wgrp, int* use_group,
                      int* use_intra, int* use_inter) {

   group__.ngrp = *ngrp;
   group__.kgrp = kgrp;
   group__.grplist = grplist;
   group__.igrp = igrp;
   group__.grpmass = grpmass;
   group__.wgrp = wgrp;
   group__.use_group = *use_group;
   group__.use_intra = *use_intra;
   group__.use_inter = *use_inter;
}

void set_imptor_data_ (int* nitors, int* iitors, double* itors1,
                       double* itors2, double* itors3) {

   imptor__.nitors = *nitors;
   imptor__.iitors = iitors;
   imptor__.itors1 = itors1;
   imptor__.itors2 = itors2;
   imptor__.itors3 = itors3;
}

void set_inform_data_ (int* digits, int* iprint, int* iwrite, int* isend,
                       int* silent, int* verbose, int* debug, int* holdup,
                       int* abort) {

   inform__.digits = *digits;
   inform__.iprint = *iprint;
   inform__.iwrite = *iwrite;
   inform__.isend = *isend;
   inform__.silent = *silent;
   inform__.verbose = *verbose;
   inform__.debug = *debug;
   inform__.holdup = *holdup;
   inform__.abort = *abort;
}

void set_ktrtor_data_ (int* maxntt, int* maxtgrd, int* maxtgrd2,
                       int* tnx, int* tny, double* ttx, double* tty,
                       double* tbf, double* tbx, double* tby,
                       double* tbxy, char20* ktt) {

   ktrtor__.maxntt = *maxntt;
   ktrtor__.maxtgrd = *maxtgrd;
   ktrtor__.maxtgrd2 = *maxtgrd2;
   ktrtor__.tnx = tnx;
   ktrtor__.tny = tny;
   ktrtor__.ttx = ttx;
   ktrtor__.tty = tty;
   ktrtor__.tbf = tbf;
   ktrtor__.tbx = tbx;
   ktrtor__.tby = tby;
   ktrtor__.tbxy = tbxy;
   ktrtor__.ktt = ktt;
}

void set_kvdws_data_ (double* rad, double* eps, double* rad4,
                      double* eps4, double* reduct) {

   kvdws__.rad = rad;
   kvdws__.eps = eps;
   kvdws__.rad4 = rad4;
   kvdws__.eps4 = eps4;
   kvdws__.reduct = reduct;
}

void set_limits_data_ (double* vdwcut, double* chgcut, double* dplcut,
                       double* mpolecut, double* vdwtaper, double* chgtaper,
                       double* dpltaper, double* mpoletaper, double* ewaldcut,
                       double* usolvcut, int* use_ewald, int* use_lights,
                       int* use_list, double* ctcut, double* cttaper, int* use_ctlist,
                       int* use_vlist, int* use_clist,
                       int* use_mlist, int* use_ulist) {

   limits__.vdwcut = *vdwcut;
   limits__.ctcut = *ctcut;
   limits__.chgcut = *chgcut;
   limits__.dplcut = *dplcut;
   limits__.mpolecut = *mpolecut;
   limits__.vdwtaper = *vdwtaper;
   limits__.cttaper = *cttaper;
   limits__.chgtaper = *chgtaper;
   limits__.dpltaper = *dpltaper;
   limits__.mpoletaper = *mpoletaper;
   limits__.ewaldcut = *ewaldcut;
   limits__.usolvcut = *usolvcut;
   limits__.use_ewald = *use_ewald;
   limits__.use_lights = *use_lights;
   limits__.use_list = *use_list;
   limits__.use_vlist = *use_vlist;
   limits__.use_ctlist = *use_ctlist;
   limits__.use_clist = *use_clist;
   limits__.use_mlist = *use_mlist;
   limits__.use_ulist = *use_ulist;
}

void set_mdstuf_data_ (int* nfree, int* irest, int* bmnmix, int* dorest,
                       int* velsave, int* frcsave, int* uindsave,
                       char* integrate) {

   mdstuf__.nfree = *nfree;
   mdstuf__.irest = *irest;
   mdstuf__.velsave = *bmnmix;
   mdstuf__.velsave = *dorest;
   mdstuf__.velsave = *velsave;
   mdstuf__.frcsave = *frcsave;
   mdstuf__.uindsave = *uindsave;
   setNullTerminator (integrate, 11, mdstuf__.integrate);
}

void set_moldyn_data_ (double* v, double* a, double* aalt) {

   moldyn__.v = v;
   moldyn__.a = a;
   moldyn__.aalt = aalt;
}

void set_molcul_data_ (int* nmol, int* imol, int* kmol, int* molcule, double* totmass, double* molmass) {

   molcul__.nmol = *nmol;
   molcul__.imol = imol;
   molcul__.kmol = kmol;
   molcul__.molcule = molcule;
   molcul__.totmass = *totmass;
   molcul__.molmass = molmass;
}

void set_mplpot_data_ (double* m2scale, double* m3scale,
                       double* m4scale, double* m5scale) {

   mplpot__.m2scale = *m2scale;
   mplpot__.m3scale = *m3scale;
   mplpot__.m4scale = *m4scale;
   mplpot__.m5scale = *m5scale;
}

void set_mpole_data_ (int* maxpole, int* npole, int* ipole, int* polsiz,
                      int* pollist, int* zaxis, int* xaxis, int* yaxis,
                      double* pole, double* rpole, char* polaxe ) {

   mpole__.maxpole = *maxpole;
   mpole__.npole = *npole;
   mpole__.ipole = ipole;
   mpole__.polsiz = polsiz;
   mpole__.pollist = pollist;
   mpole__.zaxis = zaxis;
   mpole__.xaxis = xaxis;
   mpole__.yaxis = yaxis;
   mpole__.pole = pole;
   mpole__.rpole = rpole;
   mpole__.polaxe = polaxe;
}

void set_mutant_data_ (int* nmut, int* imut, int* type0, int* class0,
                       int* type1, int* class1, double* lambda,
                       double* vlambda, double* elambda, double* scexp,
                       double* scalpha, int* mut) {

   mutant__.nmut = *nmut;
   mutant__.imut = imut;
   mutant__.type0 = type0;
   mutant__.class0 = class0;
   mutant__.type1 = type1;
   mutant__.class1 = class1;
   mutant__.lambda = *lambda;
   mutant__.vlambda = *vlambda;
   mutant__.elambda = *elambda;
   mutant__.scexp = *scexp;
   mutant__.scalpha = *scalpha;
   mutant__.mut = mut;
}

void set_nonpol_data_ (double* epso, double* epsh, double* rmino,
                       double* rminh, double* awater, double* slevy,
                       double* solvprs, double* surften, double* spcut,
                       double* spoff, double* stcut, double* stoff,
                       double* rcav, double* rdisp, double* cdisp) {

   nonpol__.epso = *epso;
   nonpol__.epsh = *epsh;
   nonpol__.rmino = *rmino;
   nonpol__.rminh = *rminh;
   nonpol__.awater = *awater;
   nonpol__.slevy = *slevy;
   nonpol__.solvprs = *solvprs;
   nonpol__.surften = *surften;
   nonpol__.spcut = *spcut;
   nonpol__.spoff = *spoff;
   nonpol__.stcut = *stcut;
   nonpol__.stoff = *stoff;
   nonpol__.rcav = rcav;
   nonpol__.rdisp = rdisp;
   nonpol__.cdisp = cdisp;
}

void set_opbend_data_ (int* nopbend, int* iopb, double* opbk) {

   opbend__.nopbend = *nopbend;
   opbend__.iopb = iopb;
   opbend__.opbk = opbk;
}

void set_pitors_data_ (int* npitors, int* ipit, double* kpit) {

   pitors__.npitors = *npitors;
   pitors__.ipit = ipit;
   pitors__.kpit = kpit;
}

void set_pme_data_ (int* nfft1, int* nfft2, int* nfft3, int* bsorder,
                    int* igrid, double* bsmod1, double* bsmod2,
                    double* bsmod3, double* bsbuild, double*** thetai1,
                    double*** thetai2, double*** thetai3, double**** qgrid,
                    double*** qfac) {

   pme__.nfft1 = *nfft1;
   pme__.nfft2 = *nfft2;
   pme__.nfft3 = *nfft3;
   pme__.bsorder = *bsorder;
   pme__.igrid = igrid;
   pme__.bsmod1 = bsmod1;
   pme__.bsmod2 = bsmod2;
   pme__.bsmod3 = bsmod3;
   pme__.bsbuild = bsbuild;
   pme__.thetai1 = thetai1;
   pme__.thetai2 = thetai2;
   pme__.thetai3 = thetai3;
   pme__.qgrid = qgrid;
   pme__.qfac = qfac;
}

void set_polar_data_ (int* maxopt, int* npolar, int* coptmax, int*optlevel,
                      double* copt, double* copm, double* polarity,
                      double* penalpha, double* thole, double* dirdamp, 
                      double* pdamp, double* udir, 
                      double* udirp, double* udirs, double* udirps,
                      double* uind, double* uinp, double* uinds,
                      double* uinps, double* uopt, double* uoptp,
                      double* uopts, double* uoptps, double* fopt,
                      double* foptp, double* uexact, int* douind) {

   polar__.maxopt = *maxopt;
   polar__.npolar = *npolar;
   polar__.coptmax = *coptmax;
   polar__.optlevel = *optlevel;
   polar__.copt = copt;
   polar__.copm = copm;
   polar__.polarity = polarity;
   polar__.penalpha = penalpha;
   polar__.thole = thole;
   polar__.dirdamp = dirdamp;
   polar__.pdamp = pdamp;
   polar__.udir = udir;
   polar__.udirp = udirp;
   polar__.udirs = udirs;
   polar__.udirps = udirps;
   polar__.uind = uind;
   polar__.uinp = uinp;
   polar__.uinds = uinds;
   polar__.uinps = uinps;
   polar__.uopt = uopt;
   polar__.uoptp = uoptp;
   polar__.uopts = uopts;
   polar__.uoptps = uoptps;
   polar__.fopt = fopt;
   polar__.foptp = foptp;
   polar__.uexact = uexact;
   polar__.douind = douind;
}

void set_polgrp_data_ (int* maxp11, int* maxp12, int* maxp13, int* maxp14,
                       int* np11, int* np12, int* np13, int* np14,
                       int* ip11, int* ip12, int* ip13, int* ip14) {

   polgrp__.maxp11 = *maxp11;
   polgrp__.maxp12 = *maxp12;
   polgrp__.maxp13 = *maxp13;
   polgrp__.maxp14 = *maxp14;
   polgrp__.np11 = np11;
   polgrp__.np12 = np12;
   polgrp__.np13 = np13;
   polgrp__.np14 = np14;
   polgrp__.ip11 = ip11;
   polgrp__.ip12 = ip12;
   polgrp__.ip13 = ip13;
   polgrp__.ip14 = ip14;
}

void set_polpot_data_ (int* politer, double* poleps, double* p2scale,
                       double* p3scale, double* p4scale, double* p5scale,
                       double* p41scale, double* d1scale, double* d2scale,
                       double* d3scale, double* d4scale, double* u1scale,
                       double* u2scale, double* u3scale, double* u4scale,
                       double* udiag, char* poltyp) {

   polpot__.politer = *politer;
   polpot__.poleps = *poleps;
   polpot__.p2scale = *p2scale;
   polpot__.p3scale = *p3scale;
   polpot__.p4scale = *p4scale;
   polpot__.p41scale = *p41scale;
   polpot__.p5scale = *p5scale;
   polpot__.d1scale = *d1scale;
   polpot__.d2scale = *d2scale;
   polpot__.d3scale = *d3scale;
   polpot__.d4scale = *d4scale;
   polpot__.u1scale = *u1scale;
   polpot__.u2scale = *u2scale;
   polpot__.u3scale = *u3scale;
   polpot__.u4scale = *u4scale;
   polpot__.udiag = *udiag;
   setNullTerminator (poltyp, 6, polpot__.poltyp);
}

void set_potent_data_ (int* use_bond, int* use_angle, int* use_strbnd,
                       int* use_urey, int* use_angang, int* use_opbend,
                       int* use_opdist, int* use_improp, int* use_imptor,
                       int* use_tors, int* use_pitors, int* use_strtor,
                       int* use_angtor, int* use_tortor, int* use_vdw,
                       int* use_charge, int* use_chgdpl, int* use_dipole,
                       int* use_mpole, int* use_polar, int* use_rxnfld,
                       int* use_solv, int* use_metal, int* use_geom,
                       int* use_extra, int* use_born, int* use_orbit, 
                       int* use_ct) {

   potent__.use_bond = *use_bond;
   potent__.use_angle = *use_angle;
   potent__.use_urey = *use_urey;
   potent__.use_strbnd = *use_strbnd;
   potent__.use_angang = *use_angang;
   potent__.use_opbend = *use_opbend;
   potent__.use_opdist = *use_opdist;
   potent__.use_improp = *use_improp;
   potent__.use_imptor = *use_imptor;
   potent__.use_tors = *use_tors;
   potent__.use_pitors = *use_pitors;
   potent__.use_strtor = *use_strtor;
   potent__.use_angtor = *use_angtor;
   potent__.use_tortor = *use_tortor;
   potent__.use_vdw = *use_vdw;
   potent__.use_charge = *use_charge;
   potent__.use_chgdpl = *use_chgdpl;
   potent__.use_dipole = *use_dipole;
   potent__.use_mpole = *use_mpole;
   potent__.use_polar = *use_polar;
   potent__.use_rxnfld = *use_rxnfld;
   potent__.use_solv = *use_solv;
   potent__.use_metal = *use_metal;
   potent__.use_geom = *use_geom;
   potent__.use_extra = *use_extra;
   potent__.use_born = *use_born;
   potent__.use_orbit = *use_orbit;
   potent__.use_ct = *use_ct;
}

void set_restrn_data_ (int* npfix, int* ndfix, int* nafix, int* ntfix,
                       int* ngfix, int* nchir, int* ipfix, int* kpfix,
                       int* idfix, int* iafix, int* itfix, int* igfix,
                       int* ichir, double* depth, double* width,
                       double* rwall, double* xpfix, double* ypfix,
                       double* zpfix, double* pfix, double* dfix,
                       double* afix, double* tfix, double* gfix,
                       double* chir, int* use_basin, int* use_wall) {

   restrn__.npfix = *npfix;
   restrn__.ndfix = *ndfix;
   restrn__.nafix = *nafix;
   restrn__.ntfix = *ntfix;
   restrn__.ngfix = *ngfix;
   restrn__.nchir = *nchir;
   restrn__.ipfix = ipfix;
   restrn__.kpfix = kpfix;
   restrn__.idfix = idfix;
   restrn__.iafix = iafix;
   restrn__.itfix = itfix;
   restrn__.igfix = igfix;
   restrn__.ichir = ichir;
   restrn__.depth = *depth;
   restrn__.width = *width;
   restrn__.rwall = *rwall;
   restrn__.xpfix = xpfix;
   restrn__.ypfix = ypfix;
   restrn__.zpfix = zpfix;
   restrn__.pfix = pfix;
   restrn__.dfix = dfix;
   restrn__.afix = afix;
   restrn__.tfix = tfix;
   restrn__.gfix = gfix;
   restrn__.chir = chir;
   restrn__.use_basin = *use_basin;
   restrn__.use_wall = *use_wall;
}

void set_sizes_data_ (int* maxatm, int* maxtyp, int* maxclass, int* maxval,
                      int* maxref, int* maxgrp, int* maxres, int* maxfix) {

   sizes__.maxatm = *maxatm;
   sizes__.maxtyp = *maxtyp;
   sizes__.maxclass = *maxclass;
   sizes__.maxval = *maxval;
   sizes__.maxref = *maxref;
   sizes__.maxgrp = *maxgrp;
   sizes__.maxres = *maxres;
   sizes__.maxfix = *maxfix;
}

void set_solute_data_ (double* doffset, double* p1, double* p2, double* p3,
                       double* p4, double* p5, double* rsolv, double* asolv,
                       double* rborn, double* drb, double* drbp, double* drobc,
                       double* gpol, double* shct, double* aobc, double* bobc,
                       double* gobc, double* vsolv, double* wace,
                       double* s2ace, double* uace, char* solvtyp,
                       char* borntyp) {

   solute__.doffset = *doffset;
   solute__.p1 = *p1;
   solute__.p2 = *p2;
   solute__.p3 = *p3;
   solute__.p4 = *p4;
   solute__.p5 = *p5;
   solute__.rsolv = rsolv;
   solute__.asolv = asolv;
   solute__.rborn = rborn;
   solute__.drb = drb;
   solute__.drbp = drbp;
   solute__.drobc = drobc;
   solute__.gpol = gpol;
   solute__.shct = shct;
   solute__.aobc = aobc;
   solute__.bobc = bobc;
   solute__.gobc = gobc;
   solute__.vsolv = vsolv;
   solute__.wace = wace;
   solute__.s2ace = s2ace;
   solute__.uace = uace;
   setNullTerminator (solvtyp, 8, solute__.solvtyp);
   setNullTerminator (borntyp, 8, solute__.borntyp);
}

void set_stodyn_data_ (double* friction, double* fgamma, int* use_sdarea) {

   stodyn__.friction = *friction;
   stodyn__.fgamma = fgamma;
   stodyn__.use_sdarea = *use_sdarea;
}

void set_strbnd_data_ (int* nstrbnd, int* isb, double* sbk) {

   strbnd__.nstrbnd = *nstrbnd;
   strbnd__.isb = isb;
   strbnd__.sbk = sbk;
}

void set_strtor_data_ (int* nstrtor, int* ist, double* kst) {

   strtor__.nstrtor = *nstrtor;
   strtor__.ist = ist;
   strtor__.kst = kst;
}

void set_torpot_data_ (double* idihunit, double* itorunit, double* torsunit,
                       double* ptorunit, double* storunit, double* atorunit,
                       double* ttorunit) {

   torpot__.idihunit = *idihunit;
   torpot__.itorunit = *itorunit;
   torpot__.torsunit = *torsunit;
   torpot__.ptorunit = *ptorunit;
   torpot__.storunit = *storunit;
   torpot__.atorunit = *atorunit;
   torpot__.ttorunit = *ttorunit;
}

void set_tors_data_ (int* ntors, int* itors, double* tors1,
                     double* tors2, double* tors3, double* tors4,
                     double* tors5, double* tors6) {

   tors__.ntors = *ntors;
   tors__.itors = itors;
   tors__.tors1 = tors1;
   tors__.tors2 = tors2;
   tors__.tors3 = tors3;
   tors__.tors4 = tors4;
   tors__.tors5 = tors5;
   tors__.tors6 = tors6;
}

void set_tortor_data_ (int* ntortor, int* itt) {

   tortor__.ntortor = *ntortor;
   tortor__.itt = itt;
}

void set_units_data_ (double* avogadro, double* lightspd, double* boltzmann,
                      double* gasconst, double* emass, double* planck,
                      double* joule, double* convert, double* bohr,
                      double* hartree, double* evolt, double* efreq,
                      double* coulomb, double* debye, double* prescon) {

   units__.avogadro = *avogadro;
   units__.lightspd = *lightspd;
   units__.boltzmann = *boltzmann;
   units__.gasconst = *gasconst;
   units__.emass = *emass;
   units__.planck = *planck;
   units__.joule = *joule;
   units__.convert = *convert;
   units__.bohr = *bohr;
   units__.hartree = *hartree;
   units__.evolt = *evolt;
   units__.efreq = *efreq;
   units__.coulomb = *coulomb;
   units__.debye = *debye;
   units__.prescon = *prescon;
}

void set_urey_data_ (int* nurey, int* iury, double* uk, double* ul) {

   urey__.nurey = *nurey;
   urey__.iury = iury;
   urey__.uk = uk;
   urey__.ul = ul;
}

void set_urypot_data_ (double* cury, double* qury, double* ureyunit) {

   urypot__.cury = *cury;
   urypot__.qury = *qury;
   urypot__.ureyunit = *ureyunit;
}

void set_usage_data_ (int* nuse, int* iuse, int* use) {

   usage__.nuse = nuse;
   usage__.iuse = iuse;
   usage__.use = use;
}

void set_vdw_data_ (int* nvdw, int* ivdw, int* jvdw, int* ired,
                    double* kred, double* radmin, double* epsilon,
                    double* radmin4, double* epsilon4, double* radhbnd,
                    double* epshbnd) {

   vdw__.nvdw = *nvdw;
   vdw__.ivdw = ivdw;
   vdw__.jvdw = jvdw;
   vdw__.ired = ired;
   vdw__.kred = kred;
   vdw__.radmin = radmin;
   vdw__.epsilon = epsilon;
   vdw__.radmin4 = radmin4;
   vdw__.epsilon4 = epsilon4;
   vdw__.radhbnd = radhbnd;
   vdw__.epshbnd = epshbnd;
}

void set_ctran_data_ (int* nct, int* ict, int* jct,
                   double* apre, double* bexp,
                   double* ct2scale, double* ct3scale,
                   double* ct4scale, double* ct5scale,
                   char* aprerule, char* bexprule) {

   ctran__.nct = *nct;
   ctran__.ict = ict;
   ctran__.jct = jct;
   ctran__.apre = apre;
   ctran__.bexp = bexp;
   ctran__.ct2scale = *ct2scale;
   ctran__.ct3scale = *ct3scale;
   ctran__.ct4scale = *ct4scale;
   ctran__.ct5scale = *ct5scale;
   setNullTerminator (aprerule, 10, ctran__.aprerule);
   setNullTerminator (bexprule, 10, ctran__.bexprule);
}

void set_vdwpot_data_ (int* maxgauss, int* ngauss, double* igauss,
                       double* abuck, double* bbuck, double* cbuck,
                       double* ghal, double* dhal, double* v2scale,
                       double* v3scale, double* v4scale, double* v5scale,
                       int* use_vcorr, char* vdwindex, char* radtyp,
                       char* radsiz, char* gausstyp, char* radrule,
                       char* epsrule, char* vdwtyp) {

   vdwpot__.maxgauss = *maxgauss;
   vdwpot__.ngauss = *ngauss;
   vdwpot__.igauss = igauss;
   vdwpot__.abuck = *abuck;
   vdwpot__.bbuck = *bbuck;
   vdwpot__.cbuck = *cbuck;
   vdwpot__.ghal = *ghal;
   vdwpot__.dhal = *dhal;
   vdwpot__.v2scale = *v2scale;
   vdwpot__.v3scale = *v3scale;
   vdwpot__.v4scale = *v4scale;
   vdwpot__.v5scale = *v5scale;
   vdwpot__.use_vcorr = *use_vcorr;
   setNullTerminator (vdwindex, 5, vdwpot__.vdwindex);
   setNullTerminator (radtyp, 5, vdwpot__.radtyp);
   setNullTerminator (radsiz, 8, vdwpot__.radsiz);
   setNullTerminator (gausstyp, 8, vdwpot__.gausstyp);
   setNullTerminator (radrule, 10, vdwpot__.radrule);
   setNullTerminator (epsrule, 10, vdwpot__.epsrule);
   setNullTerminator (vdwtyp, 13, vdwpot__.vdwtyp);
}

/*
 *    ############################################################
 *              Find Available CUDA-Enabled GPU Devices
 *    ############################################################
 */

static char cudaDevices[MAX_CUDA_DEVICES];

void set_cuda_devices_ (char* devices) {

   for (int i = 0; i < MAX_CUDA_DEVICES; ++i) {
      if (devices[i] == ' ') {
         cudaDevices[i] = '\0';
         break;
      } else if ((devices[i] < '0' || devices[i] > '9') && devices[i] != ',') {
         cudaDevices[0] = '\0';
         break;
      }
      cudaDevices[i] = devices[i];
   }
}

}

/*
 *    ############################################################
 *          Setup Masses, COM Motion Removal and Constraints
 *    ############################################################
 */

static void setupSystemParticles (OpenMM_System* system, FILE* log) {

   int ii;
   for (ii = 0; ii < atoms__.n; ii++) {
      OpenMM_System_addParticle (system, atomid__.mass[ii]);
   }
}

static void setupCMMotionRemover (OpenMM_System* system, FILE* log) {

   int frequency = mdstuf__.irest > 0 ? mdstuf__.irest : 100;
   OpenMM_CMMotionRemover* cMMotionRemover;
   cMMotionRemover = OpenMM_CMMotionRemover_create (frequency);
   OpenMM_System_addForce (system, (OpenMM_Force*) cMMotionRemover);

}

struct ConstraintMap {

   int* constraintOffset;  // offset into constraint list
   int* constraintCount;   // number of constraints for atom i
   int* constraintList;    // list of constraints sorted by atom index

   // For constraint involving atom i and j with i < j,
   //     constraintList[offset+kk] = j
   // where offset=constraintOffset[i], and 0 <= kk < constraintCount[i]
   // Note: one constraintOffset or constraintCount could be eliminated
   // since constraintCount[i] = constraintOffset[i+1] - constraintOffset[i]
};

static void freeConstraintMap (struct ConstraintMap* map) {

   free (map->constraintOffset);
   free (map->constraintCount);
   free (map->constraintList);
}

static void mapConstraints (struct ConstraintMap* map, FILE* log) {

   int ii, jj;
   int p1, p2;
   int offset, count;

   int numberOfParticles = atoms__.n;

   int* constraintCount = (int*) malloc (sizeof(int)*numberOfParticles);
   int* constraintOffset = (int*) malloc (sizeof(int)*numberOfParticles);

   memset (constraintCount, 0, sizeof(int)*numberOfParticles);
   memset (constraintOffset, 0, sizeof(int)*numberOfParticles);

   // count number of constraints for each particle where that particle
   // has the smaller index of the two constrained particles

   for (ii = 0; ii < freeze__.nrat; ii++) {
      p1 = *(freeze__.irat+2*ii) - 1;
      p2 = *(freeze__.irat + 2*ii +1) - 1;
      if (p1 > p2) {
         p1 = p2;
      }
      constraintCount[p1]++;
   }

   // set the offset value

   constraintOffset[0] = 0;
   for (ii = 1; ii < numberOfParticles; ii++){
      constraintOffset[ii] = constraintOffset[ii-1] + constraintCount[ii-1];
   }

   // allocate constraint list and load

   int* constraintList = (int*)  malloc (sizeof(int)*freeze__.nrat);
   memset (constraintCount, 0, sizeof(int)*numberOfParticles);
   for (ii = 0; ii < freeze__.nrat; ii++) {
      p1 = *(freeze__.irat+2*ii) - 1;
      p2 = *(freeze__.irat + 2*ii +1) - 1;
      if (p1 > p2) {
         int p3 = p2;
         p2 = p1;
         p1 = p3;
      }
      offset = constraintOffset[p1];
      count = constraintCount[p1];
      constraintCount[p1]++;
      constraintList[offset+count] = p2;
   }

   if (log && 0) {
      for (ii = 0; ii < numberOfParticles; ii++) {
         offset = constraintOffset[ii];
         count = constraintCount[ii];
         (void) fprintf (stderr, "%5d Offset=%5d count=%5d: ",
                         ii, offset, count);
         for (jj = 0; jj < count; jj++) {
              (void) fprintf (stderr, "%5d ", constraintList[offset+jj] );
         }
         (void) fprintf (stderr, "\n ");
      }
   }

   map->constraintCount = constraintCount;
   map->constraintOffset = constraintOffset;
   map->constraintList = constraintList;
}

static int checkForConstraint (struct ConstraintMap* map,
                               int p1, int p2, FILE* log) {

   int ii, jj;
   int offset;
   int match = 0;

   if (p1 > p2) {
      int p3 = p2;
      p2 = p1;
      p1 = p3;
   }

   offset = map->constraintOffset[p1];
   for (jj = 0; jj < map->constraintCount[p1] && match == 0; jj++) {
      if (map->constraintList[offset+jj] == p2) {
         match = 1;
      }
   }
   return match;
}

/*
 *    ############################################################
 *            Setup Standard AMOEBA Potential Energy Terms
 *    ############################################################
 */

static void setupAmoebaBondForce (OpenMM_System* system,
                                  int removeConstrainedBonds, FILE* log) {

   int ii, jj;
   int match;
   int* bondPtr;
   double kParameterConversion;

   struct ConstraintMap map;
   if (removeConstrainedBonds) {
      mapConstraints (&map, log);
   }

   OpenMM_AmoebaBondForce* amoebaBondForce;
   amoebaBondForce = OpenMM_AmoebaBondForce_create ();

   kParameterConversion = OpenMM_KJPerKcal
                             / (OpenMM_NmPerAngstrom*OpenMM_NmPerAngstrom);
   bondPtr = bndstr__.ibnd;
   for (ii = 0; ii < bndstr__.nbond; ii++) {
      match = removeConstrainedBonds ? checkForConstraint (&map, (*bondPtr)-1,
                                               *(bondPtr+1)-1, log ) : 0;
      if (match == 0) {
         OpenMM_AmoebaBondForce_addBond (amoebaBondForce, (*bondPtr)-1,
                   *(bondPtr+1)-1, bndstr__.bl[ii]*OpenMM_NmPerAngstrom,
                   kParameterConversion*bndpot__.bndunit*bndstr__.bk[ii]);
      }
      bondPtr += 2;
   }

   OpenMM_AmoebaBondForce_setAmoebaGlobalBondCubic (amoebaBondForce,
                     bndpot__.cbnd/OpenMM_NmPerAngstrom);
   OpenMM_AmoebaBondForce_setAmoebaGlobalBondQuartic (amoebaBondForce,
                 bndpot__.qbnd/(OpenMM_NmPerAngstrom*OpenMM_NmPerAngstrom));

   // if (OpenMM_AmoebaBondForce_getNumBonds (amoebaBondForce) > 0 ) {
        OpenMM_System_addForce (system, (OpenMM_Force*) amoebaBondForce);
   // }

   if (removeConstrainedBonds) {
       freeConstraintMap (&map);
   }
}

static void setupAmoebaAngleForce (OpenMM_System* system,
                                   int removeConstrainedAngles, FILE* log) {

   int ii, jj;
   int* angleIndexPtr;
   int match;
   char* angleTypPtr;

   // TINKER harmonic and in-plane angles are separate terms in OpenMM

   struct ConstraintMap map;
   if (removeConstrainedAngles) {
      mapConstraints (&map, log);
   }

   OpenMM_AmoebaAngleForce* amoebaAngleForce;
   amoebaAngleForce = OpenMM_AmoebaAngleForce_create ();
   OpenMM_System_addForce (system, (OpenMM_Force*) amoebaAngleForce);

   angleIndexPtr = angbnd__.iang;
   angleTypPtr = angpot__.angtyp;
   for (ii = 0; ii < angbnd__.nangle; ii++) {
      if (strncasecmp ("HARMONIC", angleTypPtr, 8) == 0 ) {
         match = removeConstrainedAngles ? checkForConstraint (&map,
                    *(angleIndexPtr)-1, *(angleIndexPtr+2)-1, log) : 0;
         if (match == 0) {
            OpenMM_AmoebaAngleForce_addAngle (amoebaAngleForce,
                   *(angleIndexPtr)-1, (*(angleIndexPtr+1))-1,
                   (*(angleIndexPtr+2))-1, angbnd__.anat[ii],
                   OpenMM_KJPerKcal*angpot__.angunit*angbnd__.ak[ii]);
         }
      }
      angleIndexPtr += 4;
      angleTypPtr += 8;
   }

   OpenMM_AmoebaAngleForce_setAmoebaGlobalAngleCubic (amoebaAngleForce,
                                                      angpot__.cang);
   OpenMM_AmoebaAngleForce_setAmoebaGlobalAngleQuartic (amoebaAngleForce,
                                                        angpot__.qang);
   OpenMM_AmoebaAngleForce_setAmoebaGlobalAnglePentic (amoebaAngleForce,
                                                       angpot__.pang);
   OpenMM_AmoebaAngleForce_setAmoebaGlobalAngleSextic (amoebaAngleForce,
                                                       angpot__.sang);

   if (removeConstrainedAngles) {
       freeConstraintMap (&map);
   }
}

static void setupAmoebaInPlaneAngleForce (OpenMM_System* system,
                                  int removeConstrainedBonds, FILE* log) {

   int ii, jj;
   int* angleIndexPtr;
   char* angleTypPtr;
   int match;

   // TINKER harmonic and in-plane angles are separate terms in OpenMM

   struct ConstraintMap map;
   if (removeConstrainedBonds) {
      mapConstraints (&map, log);
   }

   OpenMM_AmoebaInPlaneAngleForce* amoebaInPlaneAngleForce;
   amoebaInPlaneAngleForce = OpenMM_AmoebaInPlaneAngleForce_create ();
   OpenMM_System_addForce (system, (OpenMM_Force*) amoebaInPlaneAngleForce);

   angleIndexPtr = angbnd__.iang;
   angleTypPtr = angpot__.angtyp;
   for (ii = 0; ii < angbnd__.nangle; ii++) {
      if (strncasecmp ("IN-PLANE", angleTypPtr, 8 ) == 0) {
         match = removeConstrainedBonds ? checkForConstraint (&map,
                     *(angleIndexPtr)-1, (*(angleIndexPtr+2))-1, log) : 0;
         if (match == 0) {
            OpenMM_AmoebaInPlaneAngleForce_addAngle (amoebaInPlaneAngleForce,
                      *(angleIndexPtr)-1, (*(angleIndexPtr+1))-1,
                      (*(angleIndexPtr+2))-1, (*(angleIndexPtr+3))-1,
                      angbnd__.anat[ii],
                      OpenMM_KJPerKcal*angpot__.angunit*angbnd__.ak[ii]);
         }
      }
      angleIndexPtr += 4;
      angleTypPtr += 8;
   }

   OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAngleCubic
                       (amoebaInPlaneAngleForce, angpot__.cang);
   OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAngleQuartic
                       (amoebaInPlaneAngleForce, angpot__.qang);
   OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAnglePentic
                       (amoebaInPlaneAngleForce, angpot__.pang);
   OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAngleSextic
                       (amoebaInPlaneAngleForce, angpot__.sang);

   if (removeConstrainedBonds) {
      freeConstraintMap (&map);
   }
}

static void setupAmoebaStretchBendForce (OpenMM_System* system,
                                 int removeConstrainedBonds, FILE* log) {

   int ii, jj, index, abIndex, cbIndex;
   int* angleIndexPtr;
   double bondLengthAB;
   double bondLengthCB;
   double* bondLengthPtr;
   int match;

   struct ConstraintMap map;
   if (removeConstrainedBonds) {
       mapConstraints (&map, log);
   }

   OpenMM_AmoebaStretchBendForce* amoebaStretchBendForce;
   amoebaStretchBendForce = OpenMM_AmoebaStretchBendForce_create ();
   OpenMM_System_addForce (system, (OpenMM_Force*) amoebaStretchBendForce);

   for (ii = 0; ii < strbnd__.nstrbnd; ii++) {

      index = *(strbnd__.isb + 3*ii) - 1;

      abIndex = *(strbnd__.isb + 3*ii + 1);
      cbIndex = *(strbnd__.isb + 3*ii + 2);

      if (abIndex != 0) {
         bondLengthPtr = bndstr__.bl + abIndex - 1;
         bondLengthAB = (*bondLengthPtr)*OpenMM_NmPerAngstrom;
      } else {
         bondLengthAB = -1.0;
      }

      if (cbIndex != 0) {
         bondLengthPtr = bndstr__.bl + cbIndex - 1;
         bondLengthCB = (*bondLengthPtr)*OpenMM_NmPerAngstrom;
      } else {
         bondLengthCB = -1.0;
      }

      angleIndexPtr = angbnd__.iang + 4*index;

      match = removeConstrainedBonds ? checkForConstraint (&map,
                      *(angleIndexPtr)-1, (*(angleIndexPtr+2))-1, log) : 0;
      if (match == 0) {
         OpenMM_AmoebaStretchBendForce_addStretchBend
                    (amoebaStretchBendForce, *(angleIndexPtr)-1,
                    (*(angleIndexPtr+1))-1, (*(angleIndexPtr+2))-1,
                    bondLengthAB, bondLengthCB,
                    OpenMM_RadiansPerDegree*(*(angbnd__.anat +index)),
                    (OpenMM_KJPerKcal/
         OpenMM_NmPerAngstrom)*angpot__.stbnunit*(*(strbnd__.sbk+2*ii)),
                    (OpenMM_KJPerKcal/
         OpenMM_NmPerAngstrom)*angpot__.stbnunit*(*(strbnd__.sbk+2*ii+1)) );
      }
   }

   if (removeConstrainedBonds) {
      freeConstraintMap (&map);
   }
}

static void setupAmoebaUreyBradleyForce (OpenMM_System* system,
                                 int removeConstrainedBonds, FILE* log) {

   int ii, jj;
   int* angleIndexPtr;
   double kParameterConversion;
   int match;

   struct ConstraintMap map;
   if (removeConstrainedBonds) {
      mapConstraints (&map, log);
   }

   OpenMM_HarmonicBondForce* harmonicBondForce;
   harmonicBondForce = OpenMM_HarmonicBondForce_create ();
   OpenMM_System_addForce (system, (OpenMM_Force*) harmonicBondForce);

   angleIndexPtr = urey__.iury;
   kParameterConversion = urypot__.ureyunit*OpenMM_KJPerKcal/
                              (OpenMM_NmPerAngstrom*OpenMM_NmPerAngstrom)*2;
   for (ii = 0; ii < urey__.nurey; ii++) {
      match = removeConstrainedBonds ? checkForConstraint (&map,
                      *(angleIndexPtr)-1, (*(angleIndexPtr+2))-1, log ) : 0;
      if (match == 0) {
          OpenMM_HarmonicBondForce_addBond (harmonicBondForce,
                              *(angleIndexPtr)-1, (*(angleIndexPtr+2))-1,
                              OpenMM_NmPerAngstrom*urey__.ul[ii],
                              kParameterConversion*urey__.uk[ii]);
      }
      angleIndexPtr += 3;
   }

   // OpenMM_AmoebaBondForce_setAmoebaGlobalBondCubic (amoebaBondForce,
   //                                                  urypot__.cury);
   // OpenMM_AmoebaBondForce_setAmoebaGlobalBondQuartic (amoebaBondForce,
   //                                                    urypot__.qury);

   if (removeConstrainedBonds) {
      freeConstraintMap (&map);
   }
}

static void setupAmoebaOutOfPlaneBendForce (OpenMM_System* system, FILE* log) {

   int ii, index;
   int* angleIndexPtr;

   OpenMM_AmoebaOutOfPlaneBendForce* amoebaOutOfPlaneBendForce;
   amoebaOutOfPlaneBendForce = OpenMM_AmoebaOutOfPlaneBendForce_create ();
   OpenMM_System_addForce (system, (OpenMM_Force*) amoebaOutOfPlaneBendForce);

   for (ii = 0; ii < opbend__.nopbend; ii++) {
      index = *(opbend__.iopb + ii) - 1;
      angleIndexPtr = angbnd__.iang + 4*index;
      OpenMM_AmoebaOutOfPlaneBendForce_addOutOfPlaneBend
                  (amoebaOutOfPlaneBendForce, *(angleIndexPtr)-1,
                  (*(angleIndexPtr+1))-1, (*(angleIndexPtr+2))-1,
                  (*(angleIndexPtr+3))-1,
                  OpenMM_KJPerKcal*angpot__.opbunit*(*(opbend__.opbk +ii)));
   }

   OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendCubic
                   (amoebaOutOfPlaneBendForce, angpot__.cang);
   OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendQuartic
                   (amoebaOutOfPlaneBendForce, angpot__.qang);
   OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendPentic
                   (amoebaOutOfPlaneBendForce, angpot__.pang);
   OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendSextic
                   (amoebaOutOfPlaneBendForce, angpot__.sang);
}

static void setupAmoebaImproperTorsionForce (OpenMM_System* system, FILE* log) {

   int ii;
   double torsunit;
   int* torsIndexPtr;
   double* torsPtr;

   OpenMM_DoubleArray* torsion1;
   OpenMM_DoubleArray* torsion2;
   OpenMM_DoubleArray* torsion3;

   OpenMM_PeriodicTorsionForce* amoebaTorsionForce;

   torsion1 = OpenMM_DoubleArray_create(2);
   torsion2 = OpenMM_DoubleArray_create(2);
   torsion3 = OpenMM_DoubleArray_create(2);

   amoebaTorsionForce = OpenMM_PeriodicTorsionForce_create ();
   OpenMM_System_addForce (system, (OpenMM_Force*) amoebaTorsionForce);

   torsunit = OpenMM_KJPerKcal*torpot__.itorunit;
   torsIndexPtr = imptor__.iitors;
   for (ii = 0; ii < imptor__.nitors; ii++) {

      torsPtr = imptor__.itors1 + ii*4;
      OpenMM_DoubleArray_set (torsion1, 0, torsunit*(*torsPtr));
      OpenMM_DoubleArray_set (torsion1, 1, acos((*(torsPtr+2))));

      torsPtr = imptor__.itors2 + ii*4;
      OpenMM_DoubleArray_set (torsion2, 0, torsunit*(*torsPtr));
      OpenMM_DoubleArray_set (torsion2, 1, acos((*(torsPtr+2))));

      torsPtr = imptor__.itors3 + ii*4;
      OpenMM_DoubleArray_set (torsion3, 0, torsunit*(*torsPtr));
      OpenMM_DoubleArray_set (torsion3, 1, acos((*(torsPtr+2))));

      OpenMM_PeriodicTorsionForce_addTorsion (amoebaTorsionForce,
                (*torsIndexPtr) - 1, (*(torsIndexPtr+1)) - 1,
                (*(torsIndexPtr+2)) - 1, (*(torsIndexPtr+3)) - 1, 1,
                OpenMM_DoubleArray_get (torsion1,1),
                OpenMM_DoubleArray_get (torsion1,0));

      OpenMM_PeriodicTorsionForce_addTorsion (amoebaTorsionForce,
                (*torsIndexPtr) - 1, (*(torsIndexPtr+1)) - 1,
                (*(torsIndexPtr+2)) - 1, (*(torsIndexPtr+3)) - 1, 2,
                OpenMM_DoubleArray_get(torsion2,1),
                OpenMM_DoubleArray_get(torsion2,0));

      OpenMM_PeriodicTorsionForce_addTorsion (amoebaTorsionForce,
                (*torsIndexPtr) - 1, (*(torsIndexPtr+1)) - 1,
                (*(torsIndexPtr+2)) - 1, (*(torsIndexPtr+3)) - 1, 3,
                OpenMM_DoubleArray_get(torsion3,1),
                OpenMM_DoubleArray_get(torsion3,0));

      torsIndexPtr += 4;
   }

   OpenMM_DoubleArray_destroy (torsion1);
   OpenMM_DoubleArray_destroy (torsion2);
   OpenMM_DoubleArray_destroy (torsion3);
}

static void setupAmoebaTorsionForce (OpenMM_System* system, FILE* log) {

   int ii;
   double torsunit;
   int* torsIndexPtr;
   double* torsPtr;

   OpenMM_DoubleArray* torsion1;
   OpenMM_DoubleArray* torsion2;
   OpenMM_DoubleArray* torsion3;
   OpenMM_DoubleArray* torsion4;
   OpenMM_DoubleArray* torsion5;
   OpenMM_DoubleArray* torsion6;

   OpenMM_PeriodicTorsionForce* amoebaTorsionForce;

   torsion1 = OpenMM_DoubleArray_create(2);
   torsion2 = OpenMM_DoubleArray_create(2);
   torsion3 = OpenMM_DoubleArray_create(2);
   torsion4 = OpenMM_DoubleArray_create(2);
   torsion5 = OpenMM_DoubleArray_create(2);
   torsion6 = OpenMM_DoubleArray_create(2);

   amoebaTorsionForce = OpenMM_PeriodicTorsionForce_create ();
   OpenMM_System_addForce (system, (OpenMM_Force*) amoebaTorsionForce);

   torsunit = OpenMM_KJPerKcal*torpot__.torsunit;
   torsIndexPtr = tors__.itors;
   for (ii = 0; ii < tors__.ntors; ii++) {

      torsPtr = tors__.tors1 + ii*4;
      OpenMM_DoubleArray_set (torsion1, 0, torsunit*(*torsPtr));
      OpenMM_DoubleArray_set (torsion1, 1, acos((*(torsPtr+2))));

      torsPtr = tors__.tors2 + ii*4;
      OpenMM_DoubleArray_set (torsion2, 0, torsunit*(*torsPtr));
      OpenMM_DoubleArray_set (torsion2, 1, acos((*(torsPtr+2))));

      torsPtr = tors__.tors3 + ii*4;
      OpenMM_DoubleArray_set (torsion3, 0, torsunit*(*torsPtr));
      OpenMM_DoubleArray_set (torsion3, 1, acos((*(torsPtr+2))));

      torsPtr = tors__.tors4 + ii*4;
      OpenMM_DoubleArray_set (torsion4, 0, torsunit*(*torsPtr));
      OpenMM_DoubleArray_set (torsion4, 1, acos((*(torsPtr+2))));

      torsPtr = tors__.tors5 + ii*4;
      OpenMM_DoubleArray_set (torsion5, 0, torsunit*(*torsPtr));
      OpenMM_DoubleArray_set (torsion5, 1, acos((*(torsPtr+2))));

      torsPtr = tors__.tors6 + ii*4;
      OpenMM_DoubleArray_set (torsion6, 0, torsunit*(*torsPtr));
      OpenMM_DoubleArray_set (torsion6, 1, acos((*(torsPtr+2))));

      OpenMM_PeriodicTorsionForce_addTorsion (amoebaTorsionForce,
                (*torsIndexPtr) - 1, (*(torsIndexPtr+1)) - 1,
                (*(torsIndexPtr+2)) - 1, (*(torsIndexPtr+3)) - 1, 1,
                OpenMM_DoubleArray_get (torsion1,1),
                OpenMM_DoubleArray_get (torsion1,0));

      OpenMM_PeriodicTorsionForce_addTorsion (amoebaTorsionForce,
                (*torsIndexPtr) - 1, (*(torsIndexPtr+1)) - 1,
                (*(torsIndexPtr+2)) - 1, (*(torsIndexPtr+3)) - 1, 2,
                OpenMM_DoubleArray_get(torsion2,1),
                OpenMM_DoubleArray_get(torsion2,0));

      OpenMM_PeriodicTorsionForce_addTorsion (amoebaTorsionForce,
                (*torsIndexPtr) - 1, (*(torsIndexPtr+1)) - 1,
                (*(torsIndexPtr+2)) - 1, (*(torsIndexPtr+3)) - 1, 3,
                OpenMM_DoubleArray_get(torsion3,1),
                OpenMM_DoubleArray_get(torsion3,0));

      OpenMM_PeriodicTorsionForce_addTorsion (amoebaTorsionForce,
                (*torsIndexPtr) - 1, (*(torsIndexPtr+1)) - 1,
                (*(torsIndexPtr+2)) - 1, (*(torsIndexPtr+3)) - 1, 4,
                OpenMM_DoubleArray_get(torsion4,1),
                OpenMM_DoubleArray_get(torsion4,0));

      OpenMM_PeriodicTorsionForce_addTorsion (amoebaTorsionForce,
                (*torsIndexPtr) - 1, (*(torsIndexPtr+1)) - 1,
                (*(torsIndexPtr+2)) - 1, (*(torsIndexPtr+3)) - 1, 5,
                OpenMM_DoubleArray_get(torsion5,1),
                OpenMM_DoubleArray_get(torsion5,0));

      OpenMM_PeriodicTorsionForce_addTorsion (amoebaTorsionForce,
                (*torsIndexPtr) - 1, (*(torsIndexPtr+1)) - 1,
                (*(torsIndexPtr+2)) - 1, (*(torsIndexPtr+3)) - 1, 6,
                OpenMM_DoubleArray_get(torsion6,1),
                OpenMM_DoubleArray_get(torsion6,0));

      torsIndexPtr += 4;
   }

   OpenMM_DoubleArray_destroy (torsion1);
   OpenMM_DoubleArray_destroy (torsion2);
   OpenMM_DoubleArray_destroy (torsion3);
   OpenMM_DoubleArray_destroy (torsion4);
   OpenMM_DoubleArray_destroy (torsion5);
   OpenMM_DoubleArray_destroy (torsion6);
}

static void setupAmoebaPiTorsionForce (OpenMM_System* system, FILE* log) {

   int ii;
   int* piTorsIndexPtr;

   OpenMM_AmoebaPiTorsionForce* amoebaPiTorsionForce;
   amoebaPiTorsionForce = OpenMM_AmoebaPiTorsionForce_create ();
   OpenMM_System_addForce (system, (OpenMM_Force*) amoebaPiTorsionForce);

   piTorsIndexPtr = pitors__.ipit;
   for (ii = 0; ii < pitors__.npitors; ii++) {
      OpenMM_AmoebaPiTorsionForce_addPiTorsion (amoebaPiTorsionForce,
                    (*piTorsIndexPtr) -1, (*(piTorsIndexPtr+1))-1,
                    (*(piTorsIndexPtr+2))-1, (*(piTorsIndexPtr+3))-1,
                    (*(piTorsIndexPtr+4))-1, (*(piTorsIndexPtr+5))-1,
                    OpenMM_KJPerKcal*torpot__.ptorunit*pitors__.kpit[ii]);
      piTorsIndexPtr += 6;
   }
}

static void setupAmoebaStretchTorsionForce (OpenMM_System* system, FILE* log) {

   int ii, jj, index, baIndex, cbIndex, dcIndex;
   double bondLengthBA, bondLengthCB, bondLengthDC;
   double* bondLengthPtr;
   int* torsionIndexPtr;

   OpenMM_AmoebaStretchTorsionForce* amoebaStretchTorsionForce;
   amoebaStretchTorsionForce = OpenMM_AmoebaStretchTorsionForce_create ();
   OpenMM_System_addForce (system, (OpenMM_Force*) amoebaStretchTorsionForce);

   for (ii = 0; ii < strtor__.nstrtor; ii++) {
      index = *(strtor__.ist + 4*ii) - 1;
      baIndex = *(strtor__.ist + 4*ii + 1);
      cbIndex = *(strtor__.ist + 4*ii + 2);
      dcIndex = *(strtor__.ist + 4*ii + 3);

      if (baIndex != 0) {
         bondLengthPtr = bndstr__.bl + baIndex - 1;
         bondLengthBA = (*bondLengthPtr) * OpenMM_NmPerAngstrom;
      } else {
         bondLengthBA = -1.0;
      }

      if (cbIndex != 0) {
         bondLengthPtr = bndstr__.bl + cbIndex - 1;
         bondLengthCB = (*bondLengthPtr) * OpenMM_NmPerAngstrom;
      } else {
         bondLengthCB = -1.0;
      }

      if (dcIndex != 0) {
         bondLengthPtr = bndstr__.bl + dcIndex - 1;
         bondLengthDC = (*bondLengthPtr) * OpenMM_NmPerAngstrom;
      } else {
         bondLengthDC = -1.0;
      }

      torsionIndexPtr = tors__.itors + 4*index;
      int result = OpenMM_AmoebaStretchTorsionForce_addStretchTorsion
                      (amoebaStretchTorsionForce,
                       *(torsionIndexPtr)-1, *(torsionIndexPtr+1)-1,
                       *(torsionIndexPtr+2)-1, *(torsionIndexPtr+3)-1,
                       bondLengthBA, bondLengthCB, bondLengthDC,
                       (OpenMM_KJPerKcal/OpenMM_NmPerAngstrom) *
                          torpot__.storunit * (*(strtor__.kst + 9*ii)),
                       (OpenMM_KJPerKcal/OpenMM_NmPerAngstrom) *
                          torpot__.storunit * (*(strtor__.kst + 9*ii + 1)),
                       (OpenMM_KJPerKcal/OpenMM_NmPerAngstrom) *
                          torpot__.storunit * (*(strtor__.kst + 9*ii + 2)),
                       (OpenMM_KJPerKcal/OpenMM_NmPerAngstrom) *
                          torpot__.storunit * (*(strtor__.kst + 9*ii + 3)),
                       (OpenMM_KJPerKcal/OpenMM_NmPerAngstrom) *
                          torpot__.storunit * (*(strtor__.kst + 9*ii + 4)),
                       (OpenMM_KJPerKcal/OpenMM_NmPerAngstrom) *
                          torpot__.storunit * (*(strtor__.kst + 9*ii + 5)),
                       (OpenMM_KJPerKcal/OpenMM_NmPerAngstrom) *
                          torpot__.storunit * (*(strtor__.kst + 9*ii + 6)),
                       (OpenMM_KJPerKcal/OpenMM_NmPerAngstrom) *
                          torpot__.storunit * (*(strtor__.kst + 9*ii + 7)),
                       (OpenMM_KJPerKcal/OpenMM_NmPerAngstrom) *
                          torpot__.storunit * (*(strtor__.kst + 9*ii + 8)));
   }
}

static void setupAmoebaAngleTorsionForce (OpenMM_System* system, FILE* log) {

   int ii, jj, index, cbaIndex, dcbIndex;
   double angleCBA, angleDCB;
   double* anglePtr;
   int* torsionIndexPtr;

   OpenMM_AmoebaAngleTorsionForce* amoebaAngleTorsionForce;
   amoebaAngleTorsionForce = OpenMM_AmoebaAngleTorsionForce_create ();
   OpenMM_System_addForce (system, (OpenMM_Force*) amoebaAngleTorsionForce);

   for (ii = 0; ii < angtor__.nangtor; ii++) {
      index = *(angtor__.iat + 3*ii) - 1;
      cbaIndex = *(angtor__.iat + 3*ii + 1);
      dcbIndex = *(angtor__.iat + 3*ii + 2);

      if (cbaIndex !=0) {
         anglePtr = angbnd__.anat + cbaIndex - 1;
         angleCBA = *anglePtr;
      } else {
         angleCBA = 0.0;
      }

      if (dcbIndex != 0) {
         anglePtr = angbnd__.anat + dcbIndex - 1;
         angleDCB = *anglePtr;
      } else {
         angleDCB = 0.0;
      }

      torsionIndexPtr = tors__.itors + 4*index;
      int result = OpenMM_AmoebaAngleTorsionForce_addAngleTorsion
                      (amoebaAngleTorsionForce,
                       *(torsionIndexPtr)-1, *(torsionIndexPtr+1) - 1,
                       *(torsionIndexPtr+2) - 1, *(torsionIndexPtr+3) - 1,
                       angleCBA, angleDCB,
                       (OpenMM_KJPerKcal) * torpot__.atorunit
                          * (*(angtor__.kant + 6*ii)),
                       (OpenMM_KJPerKcal) * torpot__.atorunit
                          * (*(angtor__.kant + 6*ii + 1)),
                       (OpenMM_KJPerKcal) * torpot__.atorunit
                          * (*(angtor__.kant + 6*ii + 2)),
                       (OpenMM_KJPerKcal) * torpot__.atorunit
                          * (*(angtor__.kant + 6*ii + 3)),
                       (OpenMM_KJPerKcal) * torpot__.atorunit
                          * (*(angtor__.kant + 6*ii + 4)),
                       (OpenMM_KJPerKcal) * torpot__.atorunit
                          * (*(angtor__.kant + 6*ii + 5)));
   }
}

static int getChiralIndex (int atomB, int atomC, int atomD) {

   int ii, j, m, k;
   int chiralAtom;

   // test for chirality at the central torsion-torsion site

   chiralAtom = -1;
   if (*(couple__.n12 + atomC) == 4) {
      j = 0;
      for (ii = 0; ii < 4; ii++) {
         m = *(couple__.i12 + sizes__.maxval*atomC + ii) - 1;
         if (m != atomB && m != atomD) {
            if (j == 0) {
               j = m;
            } else {
               k = m;
            }
         }
      }
      if (atoms__.type[j] > atoms__.type[k])  chiralAtom = j;
      if (atoms__.type[k] > atoms__.type[j])  chiralAtom = k;
      if (atomid__.atomic[j] > atomid__.atomic[k])  chiralAtom = j;
      if (atomid__.atomic[k] > atomid__.atomic[j])  chiralAtom = k;
   }
   return chiralAtom;
}

static void setupAmoebaTorsionTorsionForce (OpenMM_System* system, FILE* log) {

   int ii, jj, kk, index, count;
   int ia, ib, ic, id, ie, ichiral;
   int gridIndex;
   int xIndex, yIndex;
   int* ibitorPtr;
   int* ittPtr;
   int maxntt, maxtgrd, maxtgrd2;
   int addIndex;
   int numberOfTorsionTorsions;
   OpenMM_DoubleArray* values;
   OpenMM_3D_DoubleArray* grid;

   OpenMM_AmoebaTorsionTorsionForce* amoebaTorsionTorsionForce;
   amoebaTorsionTorsionForce = OpenMM_AmoebaTorsionTorsionForce_create ();
   OpenMM_System_addForce (system, (OpenMM_Force*) amoebaTorsionTorsionForce);

   ittPtr = tortor__.itt;
   for (ii = 0; ii < tortor__.ntortor; ii++) {
      index = *(ittPtr) - 1;
      gridIndex = *(ittPtr+1) - 1;
      ibitorPtr = bitor__.ibitor + index*5;
      if (*(ittPtr+2) == 1) {
         count = 0;
         ia = *(ibitorPtr + count++) - 1;
         ib = *(ibitorPtr + count++) - 1;
         ic = *(ibitorPtr + count++) - 1;
         id = *(ibitorPtr + count++) - 1;
         ie = *(ibitorPtr + count++) - 1;
      } else {
         count = 4;
         ia = *(ibitorPtr + count--) - 1;
         ib = *(ibitorPtr + count--) - 1;
         ic = *(ibitorPtr + count--) - 1;
         id = *(ibitorPtr + count--) - 1;
         ie = *(ibitorPtr + count--) - 1;
      }
      ichiral = getChiralIndex (ib, ic, id);
      OpenMM_AmoebaTorsionTorsionForce_addTorsionTorsion
                      (amoebaTorsionTorsionForce, ia, ib, ic, id, ie,
                       ichiral, gridIndex );
      ittPtr += 3;
   }

   numberOfTorsionTorsions  =
                 OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsions
                    (amoebaTorsionTorsionForce);

   maxntt = ktrtor__.maxntt;
   maxtgrd = ktrtor__.maxtgrd;
   maxtgrd2 = maxtgrd*maxtgrd;
   values = OpenMM_DoubleArray_create(6);

   if (numberOfTorsionTorsions) {

      int kk = 0;
      for (int i = 0; i < maxntt; i++) {
         char char21[21];
         for (int j = 0; j < 20; j++) {
            char21[j] = ktrtor__.ktt[i].s20[j];
         }
         char21[20] = '\n';
         if (char21[0] != ' ') {
            kk++;
         }
      }

      for (ii = 0; ii < kk; ii++) {
         grid = OpenMM_3D_DoubleArray_create (*(ktrtor__.tnx+ii),
                                              *(ktrtor__.tny+ii), 6);
         xIndex = 0;
         yIndex = 0;

         for (jj = 0; jj < *(ktrtor__.tnx+ii)*(*(ktrtor__.tny+ii)); jj++) {
            addIndex  = 0;
            OpenMM_DoubleArray_set (values, addIndex++,
                *(ktrtor__.ttx + maxtgrd*ii+xIndex));
            OpenMM_DoubleArray_set (values, addIndex++,
                *(ktrtor__.tty + maxtgrd*ii+yIndex));

            OpenMM_DoubleArray_set (values, addIndex++,
                OpenMM_KJPerKcal*(*(ktrtor__.tbf  + maxtgrd2*ii + jj)));
            OpenMM_DoubleArray_set (values, addIndex++,
                OpenMM_KJPerKcal*(*(ktrtor__.tbx  + maxtgrd2*ii + jj)));
            OpenMM_DoubleArray_set (values, addIndex++,
                OpenMM_KJPerKcal*(*(ktrtor__.tby  + maxtgrd2*ii + jj)));
            OpenMM_DoubleArray_set (values, addIndex++,
                OpenMM_KJPerKcal*(*(ktrtor__.tbxy + maxtgrd2*ii + jj)));

            OpenMM_3D_DoubleArray_set (grid, yIndex, xIndex, values);

            xIndex++;
            if (xIndex == *(ktrtor__.tnx+ii)) {
               xIndex = 0;
               yIndex++;
            }
         }

         OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionGrid
                 (amoebaTorsionTorsionForce, ii, grid);
         OpenMM_3D_DoubleArray_destroy (grid);
      }
   }

   OpenMM_DoubleArray_destroy (values);
}

static void setDefaultPeriodicBoxVectors (OpenMM_System* system, FILE* log) {

   OpenMM_Vec3 a;
   OpenMM_Vec3 b;
   OpenMM_Vec3 c;

   if (bound__.use_bounds) {
      a.x = boxes__.lvec[0][0] * OpenMM_NmPerAngstrom;
      a.y = boxes__.lvec[1][0] * OpenMM_NmPerAngstrom;
      a.z = boxes__.lvec[2][0] * OpenMM_NmPerAngstrom;
      b.x = boxes__.lvec[0][1] * OpenMM_NmPerAngstrom;
      b.y = boxes__.lvec[1][1] * OpenMM_NmPerAngstrom;
      b.z = boxes__.lvec[2][1] * OpenMM_NmPerAngstrom;
      c.x = boxes__.lvec[0][2] * OpenMM_NmPerAngstrom;
      c.y = boxes__.lvec[1][2] * OpenMM_NmPerAngstrom;
      c.z = boxes__.lvec[2][2] * OpenMM_NmPerAngstrom;
      OpenMM_System_setDefaultPeriodicBoxVectors (system, &a, &b, &c);
   }
}

static void printDefaultPeriodicBoxVectors (OpenMM_System* system, FILE* log) {

   OpenMM_Vec3 a;
   OpenMM_Vec3 b;
   OpenMM_Vec3 c;

   OpenMM_System_getDefaultPeriodicBoxVectors (system, &a, &b, &c);
   a.x = a.x / OpenMM_NmPerAngstrom;
   a.y = a.y / OpenMM_NmPerAngstrom;
   a.z = a.z / OpenMM_NmPerAngstrom;
   b.x = b.x / OpenMM_NmPerAngstrom;
   b.y = b.y / OpenMM_NmPerAngstrom;
   b.z = b.z / OpenMM_NmPerAngstrom;
   c.x = c.x / OpenMM_NmPerAngstrom;
   c.y = c.y / OpenMM_NmPerAngstrom;
   c.z = c.z / OpenMM_NmPerAngstrom;

   if (bound__.use_bounds) {
      (void) fprintf (log, "\n Box Vectors:  %12.4f %12.4f %12.4f",
                         a.x, a.y, a.z);
      (void) fprintf (log, "\n    (Ang)      %12.4f %12.4f %12.4f",
                         b.x, b.y, b.z);
      (void) fprintf (log, "\n               %12.4f %12.4f %12.4f\n",
                         c.x, c.y, c.z);
   }
}

static void setupAmoebaChargeForce (OpenMM_System* system, FILE* log) {

   if (chgpot__.c2scale != 0.0 || chgpot__.c3scale != 0.0 || chgpot__.c5scale != 1.0) return;

   setDefaultPeriodicBoxVectors (system, log);

   char buffer[128];
   OpenMM_NonbondedForce* coulombForce;
   coulombForce = OpenMM_NonbondedForce_create ();

   for (int ii = 0; ii < atoms__.n; ++ii) {
      // set vdw params sig = 1 and eps = 0
      OpenMM_NonbondedForce_addParticle (coulombForce, charge__.pchg[ii], 1.0, 0.0);
   }

   OpenMM_BondArray* bondArray;
   bondArray = OpenMM_BondArray_create (0);
   for (int ii = 0; ii < atoms__.n; ++ii) {
      for (int jj = 0; jj < *(couple__.n12 + ii); ++jj) {
         int atomj = *(couple__.i12 + sizes__.maxval*ii + jj) - 1;
         if (ii < atomj) {
            OpenMM_BondArray_append (bondArray, ii, atomj);
         }
      }
   }
   OpenMM_NonbondedForce_createExceptionsFromBonds (coulombForce, bondArray, chgpot__.c4scale, 1.0);
   OpenMM_BondArray_destroy (bondArray);

   double cutoffdistance = 0.0;
   if (bound__.use_bounds) {
      if (limits__.use_ewald) {
         OpenMM_NonbondedForce_setNonbondedMethod (coulombForce, OpenMM_NonbondedForce_PME);
         cutoffdistance = limits__.ewaldcut;
         OpenMM_NonbondedForce_setPMEParameters (coulombForce, 10*ewald__.aewald,
                                                 pme__.nfft1, pme__.nfft2, pme__.nfft3);
      } else {
         // OpenMM_NonbondedForce_setReactionFieldDielectric (coulombForce, chgpot__.dielec);
         // OpenMM_NonbondedForce_setNonbondedMethod (coulombForce, OpenMM_NonbondedForce_CutoffPeriodic);
         // cutoffdistance = limits__.chgcut;
         // if (limits__.chgtaper < limits__.chgcut) {
         //    OpenMM_NonbondedForce_setUseSwitchingFunction (coulombForce, OpenMM_True);
         //    OpenMM_NonbondedForce_setSwitchingDistance (coulombForce, limits__.chgtaper * OpenMM_NmPerAngstrom);
         // }
         fprintf(stderr, " EXIT -- Non-PME is not supported.\n");
         exit (-1);
      }
   } else {
      // OpenMM_NonbondedForce_setNonbondedMethod (coulombForce, OpenMM_NonbondedForce_CutoffNonPeriodic);
      // cutoffdistance = limits__.chgcut;
      fprintf(stderr, " EXIT -- Nonperiodic box is not supported.\n");
      exit (-1);
   }
   cutoffdistance *= OpenMM_NmPerAngstrom;
   OpenMM_NonbondedForce_setCutoffDistance (coulombForce, cutoffdistance);
   // Turn off vdw switching
   OpenMM_NonbondedForce_setUseSwitchingFunction (coulombForce, OpenMM_False);

   OpenMM_Force_setForceGroup ((OpenMM_Force*) coulombForce, 1);
   OpenMM_System_addForce (system, (OpenMM_Force*) coulombForce);
}

static void setupAmoebaVdwForce (OpenMM_System* system, FILE* log) {

   // to detect if vdw can be calculated by AmoebaVdwForce
   bool vdwlj = false;
   bool vdwb147 = false;
   bool ok_to_go = false;
   if (strncasecmp(vdwpot__.vdwtyp, "BUFFERED-14-7", 13) == 0) {
      vdwb147 = true;
      ok_to_go = true;
      if (vdwpot__.v2scale != 0.0 || vdwpot__.v3scale != 0.0 ||
          vdwpot__.v4scale != 1.0 || vdwpot__.v5scale != 1.0) ok_to_go = false;
   } else if (strncasecmp(vdwpot__.vdwtyp, "LENNARD-JONES", 13) == 0) {
      vdwlj = true;
      ok_to_go = true;
      if (vdwpot__.v2scale != 0.0 || vdwpot__.v3scale != 0.0 ||
          vdwpot__.v5scale != 1.0) ok_to_go = false;
   }
   if (!ok_to_go) {
      fprintf (stderr, " EXIT -- vdw-12-scale and vdw-13-scale must be 0,"
                       " vdw-15-scale must be 1,"
                       " and vdw-14-scale must be 1 for BUFFERED-14-7.\n");
      exit (-1);
   }

   char buffer[128];
   int ii,jj,i;
   OpenMM_Boolean useCorrection;
   OpenMM_IntArray* exclusions;

   int maxsize = 0;
   int* vdwindex_ptr = NULL;
   setNullTerminator (vdwpot__.vdwindex, 5, buffer);
   if (strncasecmp(buffer, "CLASS", 5) == 0) {
      vdwindex_ptr = atomid__.classs;
      maxsize = sizes__.maxclass;
   } else if (strncasecmp(buffer, "TYPE", 4) == 0) {
      vdwindex_ptr = atoms__.type;
      maxsize = sizes__.maxtyp;
   }

   OpenMM_AmoebaVdwForce* amoebaVdwForce;
   amoebaVdwForce = OpenMM_AmoebaVdwForce_create ();
   OpenMM_System_addForce (system, (OpenMM_Force*) amoebaVdwForce);
   OpenMM_Force_setForceGroup ((OpenMM_Force*) amoebaVdwForce, 1);

   for (ii = 0; ii < atoms__.n; ii++) {
      i = vdw__.jvdw[ii];
      if (mutant__.mut[ii] == 0) {
         OpenMM_AmoebaVdwForce_addParticle (amoebaVdwForce,
                              vdw__.ired[ii]-1, vdwindex_ptr[ii],
                              OpenMM_NmPerAngstrom*(kvdws__.rad[i-1]),
                              OpenMM_KJPerKcal*(kvdws__.eps[i-1]),
                              vdw__.kred[ii], 1.0);

       } else {
         OpenMM_AmoebaVdwForce_addParticle (amoebaVdwForce,
                              vdw__.ired[ii]-1, vdwindex_ptr[ii],
                              OpenMM_NmPerAngstrom*(kvdws__.rad[i-1]),
                              OpenMM_KJPerKcal*(kvdws__.eps[i-1]),
                              vdw__.kred[ii], mutant__.vlambda);
      }
   }

   setNullTerminator (vdwpot__.radrule, 10, buffer);
   OpenMM_AmoebaVdwForce_setSigmaCombiningRule (amoebaVdwForce, buffer);

   setNullTerminator (vdwpot__.epsrule, 10, buffer);
   OpenMM_AmoebaVdwForce_setEpsilonCombiningRule (amoebaVdwForce, buffer);

   setNullTerminator (vdwpot__.vdwtyp, 13, buffer);
   OpenMM_AmoebaVdwForce_setFunctionalForm (amoebaVdwForce, buffer);

   OpenMM_AmoebaVdwForce_setCutoffDistance (amoebaVdwForce,
                                    limits__.vdwcut*OpenMM_NmPerAngstrom);

   useCorrection = OpenMM_False;
   if (vdwpot__.use_vcorr)  useCorrection = OpenMM_True;
   OpenMM_AmoebaVdwForce_setUseDispersionCorrection (amoebaVdwForce,
                                                     useCorrection);

   if (boxes__.orthogonal || boxes__.monoclinic || boxes__.triclinic)
      OpenMM_AmoebaVdwForce_setNonbondedMethod (amoebaVdwForce,
                                  OpenMM_AmoebaVdwForce_CutoffPeriodic);
   else
      OpenMM_AmoebaVdwForce_setNonbondedMethod (amoebaVdwForce,
                                  OpenMM_AmoebaVdwForce_NoCutoff);

   setDefaultPeriodicBoxVectors (system, log);

   exclusions = OpenMM_IntArray_create (0);
   for (ii = 0; ii < atoms__.n; ii++) {
      OpenMM_IntArray_append (exclusions, ii);
      if (fabs(vdwpot__.v2scale) <= 0.0) {
         for (jj = 0; jj < *(couple__.n12 + ii); jj++) {
            OpenMM_IntArray_append (exclusions,
                    *(couple__.i12 + sizes__.maxval*ii + jj) - 1);
         }
      }
      if (fabs(vdwpot__.v3scale) <= 0.0) {
         for (jj = 0; jj < *(couple__.n13 + ii); jj++) {
            OpenMM_IntArray_append (exclusions,
                    *(couple__.i13 + 3*sizes__.maxval*ii + jj) - 1);
         }
      }
      if (vdwlj && vdwpot__.v4scale != 1.0) {
         for (jj = 0; jj < *(couple__.n14 + ii); jj++) {
            OpenMM_IntArray_append (exclusions,
                    *(couple__.i14 + 9*sizes__.maxval*ii + jj) - 1);
         }
      }
      OpenMM_AmoebaVdwForce_setParticleExclusions (amoebaVdwForce,
                                                   ii, exclusions);
      OpenMM_IntArray_resize (exclusions, 0);
   }
   OpenMM_IntArray_destroy (exclusions);

   OpenMM_AmoebaVdwForce_computeCombinedSigmaEpsilon (amoebaVdwForce);

   int* avcs = new int[maxsize];
   for (int ii = 0; ii < maxsize; ++ii) {
      avcs[ii] = 0;
   }
   for (int ii = 0; ii < atoms__.n; ++ii) {
      int tpi = vdwindex_ptr[ii];
      if (avcs[tpi-1] == 0) {
         avcs[tpi-1] = 1;
      }
   }
   int nvdwpr = 0;
   for (int ii = 0; ii < maxsize; ++ii) {
      nvdwpr += avcs[ii];
   }
   int* avcv = new int[nvdwpr];
   int iter_tag = 0;
   for (int ii = 0; ii < maxsize; ++ii) {
      if (avcs[ii] == 1) {
         avcv[iter_tag] = ii+1;
         ++iter_tag;
      }
   }
   for (int ii = 0; ii < nvdwpr; ++ii) {
      for (int jj = ii; jj < nvdwpr; ++jj) {
         int tpi = avcv[ii];
         int tpj = avcv[jj];
         // Note: TINKER assumes the vdwpr table is of size MAXCLASS x MAXCLASS,
         // regardless of using atom class or atom type as vdwindex, as of Feb. 2017.
         int kk = (tpi-1) * sizes__.maxclass + (tpj-1);
         double combinedSigma = vdw__.radmin[kk];
         double combinedEpsilon = vdw__.epsilon[kk];
         OpenMM_AmoebaVdwForce_addVdwprByOldTypes (amoebaVdwForce, tpi, tpj,
                        combinedSigma*OpenMM_NmPerAngstrom,
                        combinedEpsilon*OpenMM_KJPerKcal);
      }
   }

   if (vdwlj && vdwpot__.v4scale != 1.0) {
      const char* ljterm = "(eps*((rm/r)^6)*((rm/r)^6-2.0))";
      char v4str[16];
      char estr[128];
      sprintf (v4str, "%.6lf", vdwpot__.v4scale);
      strcpy (estr, ljterm);
      strcat (estr, "*");
      strcat (estr, v4str);

      OpenMM_CustomBondForce* lj14Force;
      lj14Force = OpenMM_CustomBondForce_create (estr);
      OpenMM_CustomBondForce_addPerBondParameter (lj14Force, "rm");
      OpenMM_CustomBondForce_addPerBondParameter (lj14Force, "eps");

      for (ii = 0; ii < atoms__.n; ii++) {
         int tpi = vdwindex_ptr[ii];
         for (jj = 0; jj < *(couple__.n14 + ii); jj++) {
            int atomj = *(couple__.i14 + 9*sizes__.maxval*ii + jj) - 1;
            if (ii < atomj) {
               int tpj = vdwindex_ptr[atomj];
               // Note: TINKER assumes the vdwpr table is of size MAXCLASS x MAXCLASS,
               // regardless of using atom class or atom type as vdwindex, as of Feb. 2017.
               int kk = (tpi-1) * sizes__.maxclass + (tpj-1);
               double combinedSigma = vdw__.radmin[kk]*OpenMM_NmPerAngstrom;
               double combinedEpsilon = vdw__.epsilon[kk]*OpenMM_KJPerKcal;
               OpenMM_DoubleArray* coeffs;
               coeffs = OpenMM_DoubleArray_create (0);
               OpenMM_DoubleArray_append (coeffs, combinedSigma);
               OpenMM_DoubleArray_append (coeffs, combinedEpsilon);
               OpenMM_CustomBondForce_addBond (lj14Force, ii, atomj, coeffs);
               OpenMM_DoubleArray_destroy (coeffs);
            }
         }
      }

      OpenMM_Force_setForceGroup ((OpenMM_Force*) lj14Force, 1);
      OpenMM_System_addForce (system, (OpenMM_Force*) lj14Force);
   }

   delete[] avcs;
   delete[] avcv;
}

//Chengwen Below
static void setupAmoebaCTForce (OpenMM_System* system, FILE* log) {
   bool ok_to_go = false;
   // to detect if CT can be calculated by AmoebaCTForce
   if (potent__.use_ct == 1) {
      ok_to_go = true;
      if (ctran__.ct2scale != 0.0 || ctran__.ct3scale != 0.0 ||
          ctran__.ct4scale != 0.0 || ctran__.ct5scale != 0.0) ok_to_go = false;
   }
   if (!ok_to_go) {
      fprintf (stderr, " EXIT -- ct2345scales should be 0.\n");
      exit (-1);
   }

   char buffer[128];
   int ii,jj,i;
   OpenMM_IntArray* exclusions;

   int maxsize = 0;
   int* CTindex_ptr = NULL;
   CTindex_ptr = atoms__.type;

   OpenMM_AmoebaCTForce* amoebaCTForce;
   amoebaCTForce = OpenMM_AmoebaCTForce_create ();
   OpenMM_System_addForce (system, (OpenMM_Force*) amoebaCTForce);
   OpenMM_Force_setForceGroup ((OpenMM_Force*) amoebaCTForce, 1);

   for (ii = 0; ii < atoms__.n; ii++) {
      i = ctran__.jct[ii];
      OpenMM_AmoebaCTForce_addParticle (amoebaCTForce,
                           CTindex_ptr[ii], 
                           OpenMM_KJPerKcal*(ctran__.apre[i-1]),
                           (ctran__.bexp[i-1])/OpenMM_NmPerAngstrom);
   }

   setNullTerminator (ctran__.aprerule, 10, buffer);
   OpenMM_AmoebaCTForce_setApreCombiningRule (amoebaCTForce, buffer);

   setNullTerminator (ctran__.bexprule, 10, buffer);
   OpenMM_AmoebaCTForce_setBexpCombiningRule (amoebaCTForce, buffer);

   OpenMM_AmoebaCTForce_setCutoffDistance (amoebaCTForce,
                                    limits__.ctcut*OpenMM_NmPerAngstrom);

   if (boxes__.orthogonal || boxes__.monoclinic || boxes__.triclinic)
      OpenMM_AmoebaCTForce_setNonbondedMethod (amoebaCTForce,
                                  OpenMM_AmoebaCTForce_CutoffPeriodic);
   else
      OpenMM_AmoebaCTForce_setNonbondedMethod (amoebaCTForce,
                                  OpenMM_AmoebaCTForce_NoCutoff);

   setDefaultPeriodicBoxVectors (system, log);

   exclusions = OpenMM_IntArray_create (0);
   for (ii = 0; ii < atoms__.n; ii++) {
      OpenMM_IntArray_append (exclusions, ii);
      if (fabs(ctran__.ct2scale) <= 0.0) {
         for (jj = 0; jj < *(couple__.n12 + ii); jj++) {
            OpenMM_IntArray_append (exclusions,
                    *(couple__.i12 + sizes__.maxval*ii + jj) - 1);
         }
      }
      if (fabs(ctran__.ct3scale) <= 0.0) {
         for (jj = 0; jj < *(couple__.n13 + ii); jj++) {
            OpenMM_IntArray_append (exclusions,
                    *(couple__.i13 + 3*sizes__.maxval*ii + jj) - 1);
         }
      }
      if (ctran__.ct4scale != 1.0) {
         for (jj = 0; jj < *(couple__.n14 + ii); jj++) {
            OpenMM_IntArray_append (exclusions,
                    *(couple__.i14 + 9*sizes__.maxval*ii + jj) - 1);
         }
      }
      OpenMM_AmoebaCTForce_setParticleExclusions (amoebaCTForce,
                                                   ii, exclusions);
      OpenMM_IntArray_resize (exclusions, 0);
   }
   OpenMM_IntArray_destroy (exclusions);
   OpenMM_AmoebaCTForce_computeCombinedApreBexp (amoebaCTForce);
}
//Chengwen Above

static void loadCovalentArray (int numberToLoad, int* valuesToLoad,
                               OpenMM_IntArray* covaletMap) {

   int ii;
   OpenMM_IntArray_resize (covaletMap, numberToLoad);
   for (ii = 0; ii < numberToLoad; ii++) {
      OpenMM_IntArray_set (covaletMap, ii, *(valuesToLoad +ii) - 1);
   }
}

static void setupAmoebaMultipoleForce (OpenMM_System* system, FILE* log) {

   char buffer[128];
   char* axisPtr;
   int ii, jj, index;
   int invalidAxis, errorReport;
   double* polePtr;

   OpenMM_AmoebaMultipoleForce_MultipoleAxisTypes axisType;
   OpenMM_DoubleArray* dipoles;
   OpenMM_DoubleArray* quadrupoles;

   OpenMM_IntArray* covalentMap;
   OpenMM_IntArray* gridDimensions;

   OpenMM_DoubleArray* exptCoefficients;

   double dipoleConversion = OpenMM_NmPerAngstrom;
   double quadrupoleConversion = OpenMM_NmPerAngstrom*OpenMM_NmPerAngstrom;
   double polarityConversion = OpenMM_NmPerAngstrom*OpenMM_NmPerAngstrom
                                        *OpenMM_NmPerAngstrom;
   double dampingFactorConversion = sqrt(OpenMM_NmPerAngstrom);

   OpenMM_AmoebaMultipoleForce* amoebaMultipoleForce;
   amoebaMultipoleForce = OpenMM_AmoebaMultipoleForce_create();
   OpenMM_System_addForce(system, (OpenMM_Force*) amoebaMultipoleForce);
   OpenMM_Force_setForceGroup((OpenMM_Force*) amoebaMultipoleForce, 1);

   dipoles = OpenMM_DoubleArray_create(3);
   quadrupoles = OpenMM_DoubleArray_create(9);
   errorReport = 0;
   for (ii = 0; ii < atoms__.n; ii++) {

      axisPtr = mpole__.polaxe + ii*8;
      if (strncasecmp (axisPtr, "Z-then-X", 8) == 0) {
         axisType = OpenMM_AmoebaMultipoleForce_ZThenX;
      } else if (strncasecmp(axisPtr, "Bisector", 8) == 0) {
         axisType = OpenMM_AmoebaMultipoleForce_Bisector;
      } else if (strncasecmp(axisPtr, "Z-Bisect", 8) == 0) {
         axisType = OpenMM_AmoebaMultipoleForce_ZBisect;
      } else if (strncasecmp(axisPtr, "3-Fold", 6) == 0) {
         axisType = OpenMM_AmoebaMultipoleForce_ThreeFold;
      } else if (strncasecmp(axisPtr, "Z-Only", 6) == 0) {
         axisType = OpenMM_AmoebaMultipoleForce_ZOnly;
      } else if (strncasecmp(axisPtr, "None", 4) == 0
                   || strncasecmp( axisPtr, "    ",   4 ) == 0 ) {
         axisType = OpenMM_AmoebaMultipoleForce_NoAxisType;
      } else {
         errorReport++;
         setNullTerminator (axisPtr, 8, buffer);
         (void) fprintf (stderr,
    "setupAmoebaMultipoleForce: Axis Type=%s for Atom %7d Not Supported\n",
                 buffer, ii );
         if (errorReport > 20) {
            (void) fflush (stderr);
            exit (-1);
         }
      }

      polePtr = mpole__.pole + ii*mpole__.maxpole + 1;
      for (jj = 0; jj < 3; jj++) {
         OpenMM_DoubleArray_set (dipoles, jj,
                                   (*(polePtr))*dipoleConversion);
         polePtr++;
      }
      for (jj = 0; jj < 9; jj++) {
         OpenMM_DoubleArray_set (quadrupoles, jj,
                                   (*(polePtr))*quadrupoleConversion);
         polePtr++;
      }

      polePtr = mpole__.pole + ii*mpole__.maxpole;
      OpenMM_AmoebaMultipoleForce_addMultipole (amoebaMultipoleForce,
                                *polePtr, dipoles, quadrupoles, axisType,
                                *(mpole__.zaxis+ii)-1,
                                *(mpole__.xaxis+ii)-1,
                                *(mpole__.yaxis+ii)-1,
                                polar__.penalpha[ii]/OpenMM_NmPerAngstrom,
                                atomid__.atomic[ii],
                                polar__.thole[ii],
                                polar__.dirdamp[ii],
                                polar__.pdamp[ii]*dampingFactorConversion,
                                polar__.polarity[ii]*polarityConversion);

     //printf("ii, penalpha[ii] %i %f \n", ii, polar__.penalpha[ii]); 
     //printf("ii, dirdamp[ii] %i %f \n", ii, polar__.dirdamp[ii]); 
  }

   if (errorReport) {
      exit (-1);
   }

   if (limits__.use_ewald) {

      double ewaldTolerance = 1.0e-04;
      OpenMM_AmoebaMultipoleForce_setNonbondedMethod (amoebaMultipoleForce,
                                   OpenMM_AmoebaMultipoleForce_PME);
      OpenMM_AmoebaMultipoleForce_setCutoffDistance (amoebaMultipoleForce,
                                   limits__.ewaldcut*OpenMM_NmPerAngstrom);
      OpenMM_AmoebaMultipoleForce_setAEwald (amoebaMultipoleForce,
                                   ewald__.aewald/OpenMM_NmPerAngstrom);
      // OpenMM_AmoebaMultipoleForce_setPmeBSplineOrder (amoebaMultipoleForce,
      //                           pme__.bsorder);

      gridDimensions = OpenMM_IntArray_create (3);

      OpenMM_IntArray_set (gridDimensions, 0, pme__.nfft1);
      OpenMM_IntArray_set (gridDimensions, 1, pme__.nfft2);
      OpenMM_IntArray_set (gridDimensions, 2, pme__.nfft3);

      OpenMM_AmoebaMultipoleForce_setPmeGridDimensions (amoebaMultipoleForce,
                                                        gridDimensions);
      OpenMM_AmoebaMultipoleForce_setEwaldErrorTolerance
                                   (amoebaMultipoleForce, ewaldTolerance);
      OpenMM_IntArray_destroy (gridDimensions);
      setDefaultPeriodicBoxVectors (system, log);

   } else {
      OpenMM_AmoebaMultipoleForce_setNonbondedMethod (amoebaMultipoleForce,
                                   OpenMM_AmoebaMultipoleForce_NoCutoff);
   }

   if (strncasecmp (polpot__.poltyp, "DIRECT", 6) == 0) {
      OpenMM_AmoebaMultipoleForce_setPolarizationType (amoebaMultipoleForce,
                                   OpenMM_AmoebaMultipoleForce_Direct);
   } else if (strncasecmp (polpot__.poltyp, "OPT", 3) == 0) {
      if (polar__.copt[4] != 0.0) {
         exptCoefficients = OpenMM_DoubleArray_create (5);
         OpenMM_DoubleArray_set (exptCoefficients, 0, polar__.copt[0]);
         OpenMM_DoubleArray_set (exptCoefficients, 1, polar__.copt[1]);
         OpenMM_DoubleArray_set (exptCoefficients, 2, polar__.copt[2]);
         OpenMM_DoubleArray_set (exptCoefficients, 3, polar__.copt[3]);
         OpenMM_DoubleArray_set (exptCoefficients, 4, polar__.copt[4]);
      } else if (polar__.copt[3] != 0.0) {
         exptCoefficients = OpenMM_DoubleArray_create (4);
         OpenMM_DoubleArray_set (exptCoefficients, 0, polar__.copt[0]);
         OpenMM_DoubleArray_set (exptCoefficients, 1, polar__.copt[1]);
         OpenMM_DoubleArray_set (exptCoefficients, 2, polar__.copt[2]);
         OpenMM_DoubleArray_set (exptCoefficients, 3, polar__.copt[3]);
      } else if (polar__.copt[2] != 0.0) {
         exptCoefficients = OpenMM_DoubleArray_create (3);
         OpenMM_DoubleArray_set (exptCoefficients, 0, polar__.copt[0]);
         OpenMM_DoubleArray_set (exptCoefficients, 1, polar__.copt[1]);
         OpenMM_DoubleArray_set (exptCoefficients, 2, polar__.copt[2]);
      } else if (polar__.copt[1] != 0.0) {
         exptCoefficients = OpenMM_DoubleArray_create (2);
         OpenMM_DoubleArray_set (exptCoefficients, 0, polar__.copt[0]);
         OpenMM_DoubleArray_set (exptCoefficients, 1, polar__.copt[1]);
      }
      OpenMM_AmoebaMultipoleForce_setExtrapolationCoefficients
                                  (amoebaMultipoleForce,exptCoefficients);
      OpenMM_AmoebaMultipoleForce_setPolarizationType (amoebaMultipoleForce,
                                   OpenMM_AmoebaMultipoleForce_Extrapolated);
   } else {
      OpenMM_AmoebaMultipoleForce_setPolarizationType (amoebaMultipoleForce,
                                   OpenMM_AmoebaMultipoleForce_Mutual);
   }

   int PolType_out = OpenMM_AmoebaMultipoleForce_getPolarizationType
                                    (amoebaMultipoleForce);

   OpenMM_AmoebaMultipoleForce_setMutualInducedMaxIterations
                                    (amoebaMultipoleForce, 500);
   OpenMM_AmoebaMultipoleForce_setMutualInducedTargetEpsilon
                                    (amoebaMultipoleForce, polpot__.poleps);

   covalentMap = OpenMM_IntArray_create (0);
   for (ii = 0; ii < atoms__.n; ii++) {

      loadCovalentArray (*(couple__.n12 + ii),
                         (couple__.i12 + sizes__.maxval*ii), covalentMap);
      OpenMM_AmoebaMultipoleForce_setCovalentMap (amoebaMultipoleForce, ii,
                         OpenMM_AmoebaMultipoleForce_Covalent12,
                         covalentMap);

      loadCovalentArray (*(couple__.n13 + ii),
                         (couple__.i13 + 3*sizes__.maxval*ii), covalentMap);
      OpenMM_AmoebaMultipoleForce_setCovalentMap (amoebaMultipoleForce, ii,
                         OpenMM_AmoebaMultipoleForce_Covalent13,
                         covalentMap);

      loadCovalentArray (*(couple__.n14 + ii),
                         (couple__.i14 + 9*sizes__.maxval*ii), covalentMap);
      OpenMM_AmoebaMultipoleForce_setCovalentMap (amoebaMultipoleForce, ii,
                         OpenMM_AmoebaMultipoleForce_Covalent14,
                         covalentMap);

      loadCovalentArray (*(couple__.n15 + ii),
                         (couple__.i15 + 27*sizes__.maxval*ii), covalentMap);
      OpenMM_AmoebaMultipoleForce_setCovalentMap (amoebaMultipoleForce, ii,
                         OpenMM_AmoebaMultipoleForce_Covalent15,
                         covalentMap);

      loadCovalentArray (*(polgrp__.np11 + ii),
                         (polgrp__.ip11 + polgrp__.maxp11*ii), covalentMap);
      OpenMM_AmoebaMultipoleForce_setCovalentMap (amoebaMultipoleForce, ii,
                         OpenMM_AmoebaMultipoleForce_PolarizationCovalent11,
                         covalentMap);

      loadCovalentArray (*(polgrp__.np12 + ii),
                         (polgrp__.ip12 + polgrp__.maxp12*ii), covalentMap);
      OpenMM_AmoebaMultipoleForce_setCovalentMap (amoebaMultipoleForce, ii,
                         OpenMM_AmoebaMultipoleForce_PolarizationCovalent12,
                         covalentMap);

      loadCovalentArray (*(polgrp__.np13 + ii),
                         (polgrp__.ip13 + polgrp__.maxp13*ii), covalentMap);
      OpenMM_AmoebaMultipoleForce_setCovalentMap (amoebaMultipoleForce, ii,
                         OpenMM_AmoebaMultipoleForce_PolarizationCovalent13,
                         covalentMap);

      loadCovalentArray (*(polgrp__.np14 + ii),
                         (polgrp__.ip14 + polgrp__.maxp14*ii), covalentMap);
      OpenMM_AmoebaMultipoleForce_setCovalentMap( amoebaMultipoleForce, ii,
                         OpenMM_AmoebaMultipoleForce_PolarizationCovalent14,
                         covalentMap);
   }

   OpenMM_DoubleArray_destroy (dipoles);
   OpenMM_DoubleArray_destroy (quadrupoles);

   OpenMM_IntArray_destroy (covalentMap);
}

static void setupAmoebaWcaDispersionForce (OpenMM_System* system, FILE* log) {

   int ii;
   double cdispTotal = 0.0;
   double epso = 0.1100;
   double epsh = 0.0135;
   double rmino = 1.7025;
   double rminh = 1.3275;
   double awater = 0.033428;
   double slevy = 1.0;
   double dispoff = 0.26;
   double shctd = 0.81;

   OpenMM_AmoebaWcaDispersionForce*  amoebaWcaDispersionForce;
   amoebaWcaDispersionForce = OpenMM_AmoebaWcaDispersionForce_create ();
   OpenMM_System_addForce (system, (OpenMM_Force*) amoebaWcaDispersionForce);

   for (ii = 0; ii < atoms__.n; ii++) {
      cdispTotal += nonpol__.cdisp[ii];
      OpenMM_AmoebaWcaDispersionForce_addParticle (amoebaWcaDispersionForce,
              OpenMM_NmPerAngstrom*kvdws__.rad[atomid__.classs[ii]-1],
              OpenMM_KJPerKcal*kvdws__.eps[atomid__.classs[ii]-1]);
   }

   OpenMM_AmoebaWcaDispersionForce_setEpso (amoebaWcaDispersionForce,
                                            epso*OpenMM_KJPerKcal);
   OpenMM_AmoebaWcaDispersionForce_setEpsh (amoebaWcaDispersionForce,
                                            epsh*OpenMM_KJPerKcal);
   OpenMM_AmoebaWcaDispersionForce_setRmino (amoebaWcaDispersionForce,
                                             rmino*OpenMM_NmPerAngstrom);
   OpenMM_AmoebaWcaDispersionForce_setRminh (amoebaWcaDispersionForce,
                                             rminh*OpenMM_NmPerAngstrom);
   OpenMM_AmoebaWcaDispersionForce_setDispoff (amoebaWcaDispersionForce,
                                               dispoff*OpenMM_NmPerAngstrom);
   OpenMM_AmoebaWcaDispersionForce_setAwater (amoebaWcaDispersionForce,
      awater/(OpenMM_NmPerAngstrom*OpenMM_NmPerAngstrom*OpenMM_NmPerAngstrom));
   OpenMM_AmoebaWcaDispersionForce_setSlevy (amoebaWcaDispersionForce, slevy);
   OpenMM_AmoebaWcaDispersionForce_setShctd (amoebaWcaDispersionForce, shctd);
}

static double getObcShct (int atmnum) {

   double shct = 0.80;
   if (atmnum == 1)  shct = 0.85;
   if (atmnum == 6)  shct = 0.72;
   if (atmnum == 7)  shct = 0.79;
   if (atmnum == 8)  shct = 0.85;
   if (atmnum == 9)  shct = 0.88;
   if (atmnum == 15)  shct = 0.86;
   if (atmnum == 16)  shct = 0.96;
   if (atmnum == 26)  shct = 0.88;
   return shct;
}

static void setupAmoebaGeneralizedKirkwoodForce (OpenMM_System* system,
                                         int includeCavityTerm, FILE* log) {

   int ii;
   int useGrycuk;
   int useObc;
   double* polePtr;
   double shct;
   char buffer[80];

   // check Born radius type; force use of "Grycuk" for now

   setNullTerminator (solute__.borntyp, 8, buffer);
   useGrycuk = 1;
   if (strncasecmp (buffer, "GRYCUK", 6 ) != 0) {
      if (log) {
         (void) fprintf (log, "setupAmoebaGeneralizedKirkwoodForce: Born type=%s, forcing Born type to 'Grycuk'.\n", buffer);
      } else {
         (void) fprintf (stderr, "setupAmoebaGeneralizedKirkwoodForce: Born type=%s, forcing Born type to 'Grycuk'.\n", buffer);
      }
   }

   OpenMM_AmoebaGeneralizedKirkwoodForce* amoebaGeneralizedKirkwoodForce;
   amoebaGeneralizedKirkwoodForce =
       OpenMM_AmoebaGeneralizedKirkwoodForce_create ();
   OpenMM_System_addForce (system, (OpenMM_Force*)
                           amoebaGeneralizedKirkwoodForce);

   OpenMM_AmoebaGeneralizedKirkwoodForce_setSolventDielectric
                          (amoebaGeneralizedKirkwoodForce, 78.3);
   OpenMM_AmoebaGeneralizedKirkwoodForce_setSoluteDielectric
                          (amoebaGeneralizedKirkwoodForce, chgpot__.dielec);
   OpenMM_AmoebaGeneralizedKirkwoodForce_setIncludeCavityTerm
                          (amoebaGeneralizedKirkwoodForce, includeCavityTerm);
   // OpenMM_AmoebaGeneralizedGeneralizedKirkwoodForce_setDielectricOffset
   //                     (amoebaGeneralizedGeneralizedKirkwoodForce,
   //                      solute__.doffset*OpenMM_NmPerAngstrom);
   // OpenMM_AmoebaGeneralizedGeneralizedKirkwoodForce_setProbeRadius
   //                     (OpenMM_AmoebaGeneralizedGeneralizedKirkwoodForce,
   //                      0.14);
   // OpenMM_AmoebaGeneralizedGeneralizedKirkwoodForce_setSurfaceAreaFactor
   //                     (amoebaGeneralizedGeneralizedKirkwoodForce,
   //                      surfaceAreaFactor);

   for (ii = 0; ii < atoms__.n; ii++) {
      polePtr = mpole__.pole + ii*mpole__.maxpole;
      if (useGrycuk) {
         shct = solute__.shct[ii];
      } else {
         shct = getObcShct (atomid__.atomic[ii]);
      }
      OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle
                         (amoebaGeneralizedKirkwoodForce, *(polePtr),
                          OpenMM_NmPerAngstrom*solute__.rsolv[ii], shct);
   }
}

/*
 *    ############################################################
 *          Setup Geometric Restraint Potential Energy Terms
 *    ############################################################
 */

#define RADIANS_PER_DEGREE 0.0174532925

static void setupPositionalRestraints (OpenMM_System* system, FILE* log) {

   double convert = OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom*OpenMM_NmPerAngstrom);
   double nmPerAng = OpenMM_NmPerAngstrom;

   OpenMM_CustomExternalForce* force;

   if (boxes__.orthogonal) {
      force = OpenMM_CustomExternalForce_create("k*(max(0.0,periodicdistance(x, y, z, x0, y0, z0)-range))^2");
   } else if (boxes__.monoclinic) {
      force = OpenMM_CustomExternalForce_create("k*(max(0.0,periodicdistance(x, y, z, x0, y0, z0)-range))^2");
      OpenMM_CustomExternalForce_addGlobalParameter (force, "beta_sin", boxes__.beta_sin);
      OpenMM_CustomExternalForce_addGlobalParameter (force, "beta_cos", boxes__.beta_cos);
   } else if (boxes__.triclinic) {
      force = OpenMM_CustomExternalForce_create("k*(max(0.0,periodicdistance(x, y, z, x0, y0, z0)-range))^2");
      OpenMM_CustomExternalForce_addGlobalParameter (force, "beta_sin", boxes__.beta_sin);
      OpenMM_CustomExternalForce_addGlobalParameter (force, "beta_cos", boxes__.beta_cos);
      OpenMM_CustomExternalForce_addGlobalParameter (force, "gamma_term", boxes__.gamma_term);
      OpenMM_CustomExternalForce_addGlobalParameter (force, "beta_term", boxes__.beta_term);
      OpenMM_CustomExternalForce_addGlobalParameter (force, "gamma_sin", boxes__.gamma_sin);
      OpenMM_CustomExternalForce_addGlobalParameter (force, "gamma_cos", boxes__.gamma_cos);
   } else if (boxes__.octahedron) {
      force = OpenMM_CustomExternalForce_create("k*(max(0.0,periodicdistance(x, y, z, x0, y0, z0)-range))^2");
      OpenMM_CustomExternalForce_addGlobalParameter (force, "box34", boxes__.box34*nmPerAng);
      OpenMM_CustomExternalForce_addGlobalParameter (force, "xbox", *boxes__.xbox*nmPerAng);
      OpenMM_CustomExternalForce_addGlobalParameter (force, "xbox2", *boxes__.xbox2*nmPerAng);
      OpenMM_CustomExternalForce_addGlobalParameter (force, "ybox", *boxes__.ybox*nmPerAng);
      OpenMM_CustomExternalForce_addGlobalParameter (force, "ybox2", *boxes__.ybox2*nmPerAng);
      OpenMM_CustomExternalForce_addGlobalParameter (force, "zbox", *boxes__.zbox*nmPerAng);
      OpenMM_CustomExternalForce_addGlobalParameter (force, "zbox2", *boxes__.zbox2*nmPerAng);
   }

   OpenMM_CustomExternalForce_addPerParticleParameter (force, "k");
   OpenMM_CustomExternalForce_addPerParticleParameter (force, "range");
   OpenMM_CustomExternalForce_addPerParticleParameter (force, "x0");
   OpenMM_CustomExternalForce_addPerParticleParameter (force, "y0");
   OpenMM_CustomExternalForce_addPerParticleParameter (force, "z0");
   OpenMM_CustomExternalForce_addPerParticleParameter (force, "use_bounds");
   OpenMM_CustomExternalForce_addGlobalParameter (force, "xcell", cell__.xcell*nmPerAng);
   OpenMM_CustomExternalForce_addGlobalParameter (force, "xcell2", cell__.xcell2*nmPerAng);
   OpenMM_CustomExternalForce_addGlobalParameter (force, "ycell", cell__.ycell*nmPerAng);
   OpenMM_CustomExternalForce_addGlobalParameter (force, "ycell2", cell__.ycell2*nmPerAng);
   OpenMM_CustomExternalForce_addGlobalParameter (force, "zcell", cell__.zcell*nmPerAng);
   OpenMM_CustomExternalForce_addGlobalParameter (force, "zcell2", cell__.zcell2*nmPerAng);

   for (int i = 0; i < restrn__.npfix; ++i) {
      double use_bounds = 0.0;
      if (bound__.use_bounds)
         use_bounds = 1.0;
      OpenMM_DoubleArray* positionalParameters = OpenMM_DoubleArray_create(0);
      OpenMM_DoubleArray_append(positionalParameters, restrn__.pfix[i*2]*convert);
      OpenMM_DoubleArray_append(positionalParameters, restrn__.pfix[i*2 + 1] * nmPerAng);
      OpenMM_DoubleArray_append(positionalParameters, restrn__.xpfix[i] * nmPerAng);
      OpenMM_DoubleArray_append(positionalParameters, restrn__.ypfix[i] * nmPerAng);
      OpenMM_DoubleArray_append(positionalParameters, restrn__.zpfix[i] * nmPerAng);
      OpenMM_DoubleArray_append(positionalParameters, use_bounds);
      OpenMM_CustomExternalForce_addParticle(force, restrn__.ipfix[i]-1, positionalParameters);
      OpenMM_DoubleArray_destroy(positionalParameters);
   }
   OpenMM_System_addForce(system, (OpenMM_Force*) force);
}

static void setupDistanceRestraints (OpenMM_System* system, FILE* log) {

   const double convert = OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom*OpenMM_NmPerAngstrom);
   const double nmPerAng = OpenMM_NmPerAngstrom;

   OpenMM_CustomCompoundBondForce* force;

   if (boxes__.orthogonal) {
      force = OpenMM_CustomCompoundBondForce_create(2,"k*(max(max(0.0,distMin-r),max(0.0,r-distMax)))^2;\
         r = sqrt(xr^2+yr^2+zr^2);\
         xr = xr_pre-use_bounds*(xcell*ceil(max(0.0,(abs(xr_pre)-xcell2))/xcell)*(step(xr_pre)*2-1));\
         yr = yr_pre-use_bounds*(ycell*ceil(max(0.0,(abs(yr_pre)-ycell2))/ycell)*(step(yr_pre)*2-1));\
         zr = zr_pre-use_bounds*(zcell*ceil(max(0.0,(abs(zr_pre)-zcell2))/zcell)*(step(zr_pre)*2-1));\
         xr_pre=x1-x2;\
         yr_pre=y1-y2;\
         zr_pre=z1-z2");
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "xcell", cell__.xcell*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "xcell2", cell__.xcell2*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "ycell", cell__.ycell*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "ycell2", cell__.ycell2*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "zcell", cell__.zcell*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "zcell2", cell__.zcell2*nmPerAng);
   } else if (boxes__.monoclinic) {
      force = OpenMM_CustomCompoundBondForce_create(2, "k*(max(max(0.0,distMin-r),max(0.0,r-distMax)))^2;\
         r = sqrt(xr^2+yr^2+zr^2);\
         xr = xr_converted + zr_converted*beta_cos;\
         zr = zr_converted*beta_sin;\
         yr = yr_pre-use_bounds*(ycell*ceil(max(0,(abs(yr_pre)-ycell2))/ycell)*(step(yr_pre)*2-1));\
         xr_converted = xr_pre-use_bounds*(xcell*ceil(max(0,(abs(xr_pre)-xcell2))/xcell)*(step(xr_pre)*2-1));\
         zr_converted = zr_pre-use_bounds*(zcell*ceil(max(0,(abs(zr_pre)-zcell2))/zcell)*(step(zr_pre)*2-1));\
         xr_pre = x1-x2 - zr_pre*beta_cos;\
         yr_pre = y1-y2;\
         zr_pre = (z1-z2)/beta_sin");
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "beta_sin", boxes__.beta_sin);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "beta_cos", boxes__.beta_cos);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "xcell", cell__.xcell*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "xcell2", cell__.xcell2*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "ycell", cell__.ycell*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "ycell2", cell__.ycell2*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "zcell", cell__.zcell*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "zcell2", cell__.zcell2*nmPerAng);
   } else if (boxes__.triclinic) {
      force = OpenMM_CustomCompoundBondForce_create(2, "k*(max(max(0.0,distMin-r),max(0.0,r-distMax)))^2;\
         r = sqrt(xr^2+yr^2+zr^2);\
         xr = xr_converted + yr_converted*gamma_cos + zr_converted*beta_cos;\
         yr = yr_converted*gamma_sin + zr_converted*beta_term;\
         zr = zr_converted*gamma_term;\
         xr_converted = xr_pre-use_bounds*(xcell*ceil(max(0,(abs(xr_pre)-xcell2))/xcell)*(step(xr_pre)*2-1));\
         yr_converted = yr_pre-use_bounds*(ycell*ceil(max(0,(abs(yr_pre)-ycell2))/ycell)*(step(yr_pre)*2-1));\
         zr_converted = zr_pre-use_bounds*(zcell*ceil(max(0,(abs(zr_pre)-zcell2))/zcell)*(step(zr_pre)*2-1));\
         xr_pre = x1-x2 - yr_pre*gamma_cos - zr_pre*beta_cos;\
         yr_pre = (y1-y2 - zr_pre*beta_term)/gamma_sin;\
         zr_pre = (z1-z2)/gamma_term");
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "beta_sin", boxes__.beta_sin);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "beta_cos", boxes__.beta_cos);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "gamma_term", boxes__.gamma_term);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "beta_term", boxes__.beta_term);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "gamma_sin", boxes__.gamma_sin);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "gamma_cos", boxes__.gamma_cos);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "xcell", cell__.xcell*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "xcell2", cell__.xcell2*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "ycell", cell__.ycell*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "ycell2", cell__.ycell2*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "zcell", cell__.zcell*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "zcell2", cell__.zcell2*nmPerAng);
   } else if (boxes__.octahedron) {
      force = OpenMM_CustomCompoundBondForce_create(2, "k*(max(max(0.0,distMin-r),max(0.0,r-distMax)))^2;\
         r = sqrt(xr^2+yr^2+zr^2);\
         xr = xr_converted - absDist*xbox2 * (step(xr_converted)*2-1);\
         yr = yr_converted - absDist*ybox2 * (step(yr_converted)*2-1);\
         zr = zr_converted - absDist*zbox2 * (step(zr_converted)*2-1);\
         absDist = step(abs(xr_converted)+abs(yr_converted)+abs(zr_converted)-box34);\
         xr_converted = xr_pre-use_bounds*(xbox*ceil(max(0,(abs(xr_pre)-xbox2))/xbox)*(step(xr_pre)*2-1));\
         yr_converted = yr_pre-use_bounds*(ybox*ceil(max(0,(abs(yr_pre)-ybox2))/ybox)*(step(yr_pre)*2-1));\
         zr_converted = zr_pre-use_bounds*(zbox*ceil(max(0,(abs(zr_pre)-zbox2))/zbox)*(step(zr_pre)*2-1));\
         xr_pre = x1-x2;\
         yr_pre = y1-y2;\
         zr_pre = z1-z2");
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "box34", boxes__.box34*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "xbox", *boxes__.xbox*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "xbox2", *boxes__.xbox2*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "ybox", *boxes__.ybox*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "ybox2", *boxes__.ybox2*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "zbox", *boxes__.zbox*nmPerAng);
      OpenMM_CustomCompoundBondForce_addGlobalParameter (force, "zbox2", *boxes__.zbox2*nmPerAng);
   }

   OpenMM_CustomCompoundBondForce_addPerBondParameter (force, "k");
   OpenMM_CustomCompoundBondForce_addPerBondParameter (force, "distMin");
   OpenMM_CustomCompoundBondForce_addPerBondParameter (force, "distMax");
   OpenMM_CustomCompoundBondForce_addPerBondParameter (force, "use_bounds");

   for (int i = 0; i < restrn__.ndfix; ++i) {
      double use_bounds = 0.0;
      if (bound__.use_bounds && molcul__.molcule[restrn__.idfix[i*2]-1]
                               != molcul__.molcule[restrn__.idfix[i*2+1]-1]) {
         use_bounds = 1.0;
      }
      OpenMM_IntArray* particles = OpenMM_IntArray_create(0);
      OpenMM_IntArray_append(particles, restrn__.idfix[i*2]-1);
      OpenMM_IntArray_append(particles, restrn__.idfix[i*2+1]-1);
      OpenMM_DoubleArray* distanceParameters = OpenMM_DoubleArray_create(0);
      OpenMM_DoubleArray_append(distanceParameters, restrn__.dfix[i*3]*convert);
      OpenMM_DoubleArray_append(distanceParameters, restrn__.dfix[i*3 + 1]*nmPerAng);
      OpenMM_DoubleArray_append(distanceParameters, restrn__.dfix[i*3 + 2]*nmPerAng);
      OpenMM_DoubleArray_append(distanceParameters, use_bounds);
      OpenMM_CustomCompoundBondForce_addBond(force, particles, distanceParameters);
      OpenMM_IntArray_destroy(particles);
      OpenMM_DoubleArray_destroy(distanceParameters);
   }
   OpenMM_System_addForce(system, (OpenMM_Force*) force);
}

static void setupAngleRestraints (OpenMM_System* system, FILE* log) {

   double convert = OpenMM_KJPerKcal / RADIANS_PER_DEGREE / RADIANS_PER_DEGREE;

   OpenMM_CustomAngleForce* force =
      OpenMM_CustomAngleForce_create("k*(max(max(0.0,thetaMin-theta),max(0.0,theta-thetaMax)))^2");

   OpenMM_CustomAngleForce_addPerAngleParameter (force, "k");
   OpenMM_CustomAngleForce_addPerAngleParameter (force, "thetaMin");
   OpenMM_CustomAngleForce_addPerAngleParameter (force, "thetaMax");

   for (int i = 0; i < restrn__.nafix; ++i) {
      OpenMM_DoubleArray* AngleParameters = OpenMM_DoubleArray_create(0);
      OpenMM_DoubleArray_append(AngleParameters, restrn__.afix[i*2]*convert);
      OpenMM_DoubleArray_append(AngleParameters, restrn__.afix[i*2 + 1]*RADIANS_PER_DEGREE);
      OpenMM_DoubleArray_append(AngleParameters, restrn__.afix[i*2 + 2]*RADIANS_PER_DEGREE);
      OpenMM_CustomAngleForce_addAngle(force, restrn__.iafix[i*3]-1, restrn__.iafix[i*3+1]-1,
                                       restrn__.iafix[i*3+2]-1, AngleParameters);
      OpenMM_DoubleArray_destroy(AngleParameters);
   }
   OpenMM_System_addForce(system, (OpenMM_Force*) force);
}

static void setupTorsionRestraints (OpenMM_System* system, FILE* log) {

   double convert;
   convert = OpenMM_KJPerKcal / RADIANS_PER_DEGREE / RADIANS_PER_DEGREE;

   OpenMM_CustomTorsionForce* force =
         OpenMM_CustomTorsionForce_create("k*max(\
            (step(thetaMin-theta)*(min(min(abs(theta-thetaMin),abs(theta-thetaMin-6.28318530718)), abs(theta-thetaMin+6.28318530718)))),\
            (step(theta-thetaMax)*(min(min(abs(theta-thetaMax),abs(theta-thetaMax-6.28318530718)), abs(theta-thetaMax+6.28318530718))))\
               )^2");

   OpenMM_CustomTorsionForce_addPerTorsionParameter (force, "k");
   OpenMM_CustomTorsionForce_addPerTorsionParameter (force, "thetaMin");
   OpenMM_CustomTorsionForce_addPerTorsionParameter (force, "thetaMax");

   for (int i = 0; i < restrn__.ntfix; i++) {
      float thetaMin = restrn__.tfix[i*3+1];
      float thetaMax = restrn__.tfix[i*3+2];
      if (thetaMin > 180.0f)
         thetaMin -= 360.0f;
      else if (thetaMin < -180.0f)
         thetaMin += 360.0f;
      if (thetaMax > 180.0f)
         thetaMax -= 360.0f;
      else if (thetaMax < -180.0f)
         thetaMax += 360.0f;

      OpenMM_DoubleArray* lowerTorsionParameters = OpenMM_DoubleArray_create (0);
      OpenMM_DoubleArray_append (lowerTorsionParameters, restrn__.tfix[i*3]*convert);
      OpenMM_DoubleArray_append (lowerTorsionParameters, thetaMin*RADIANS_PER_DEGREE);
      OpenMM_DoubleArray_append (lowerTorsionParameters, thetaMax*RADIANS_PER_DEGREE);
      OpenMM_CustomTorsionForce_addTorsion (force, restrn__.itfix[i*4]-1,
                                            restrn__.itfix[i*4+1]-1,
                                            restrn__.itfix[i*4+2]-1,
                                            restrn__.itfix[i*4+3]-1,
                                            lowerTorsionParameters);
      OpenMM_DoubleArray_destroy(lowerTorsionParameters);
   }
   OpenMM_System_addForce(system, (OpenMM_Force*) force);
}

static OpenMM_IntArray* getGroup (int* kgrp, int* igrp, int idx) {

   int start = igrp[2 * idx] - 1;
   int end = igrp[2 * idx + 1] - 1;
   OpenMM_IntArray* group = OpenMM_IntArray_create (0);
   for (int i = start; i < end + 1; i++) {
      OpenMM_IntArray_append (group, kgrp[i] - 1);
   }
   return group;
}

static OpenMM_DoubleArray* getWeights (int* kgrp, int* igrp, int idx,
                                       OpenMM_System* system) {

   int start = igrp[2 * idx] - 1;
   int end = igrp[2 * idx + 1] - 1;
   OpenMM_DoubleArray* weights = OpenMM_DoubleArray_create (0);
   for (int i = start; i < end + 1; i++){
      OpenMM_DoubleArray_append (weights, OpenMM_System_getParticleMass(system,
                                 kgrp[i] - 1));
   }
   return weights;
}

static void setupCentroidRestraints (OpenMM_System* system, FILE* log) {

   double convert;
   convert = OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom*OpenMM_NmPerAngstrom);

   // In the expression below, u and l are the upper and lower threshold

   OpenMM_CustomCentroidBondForce* force =
      OpenMM_CustomCentroidBondForce_create (2,
                 "step(distance(g1,g2)-u)*k*(distance(g1,g2)-u)^2+step(l-distance(g1,g2))*k*(distance(g1,g2)-l)^2");
   OpenMM_CustomCentroidBondForce_addPerBondParameter (force, "k");
   OpenMM_CustomCentroidBondForce_addPerBondParameter (force, "l");
   OpenMM_CustomCentroidBondForce_addPerBondParameter (force, "u");

   for (int j = 1; j < group__.ngrp + 1; j++) {
      OpenMM_IntArray* igroup = getGroup(group__.kgrp, group__.igrp, j);
      OpenMM_DoubleArray* iweight = getWeights(group__.kgrp, group__.igrp, j, system);
      OpenMM_CustomCentroidBondForce_addGroup (force, igroup, iweight);
      OpenMM_IntArray_destroy (igroup);
      OpenMM_DoubleArray_destroy (iweight);
   }

   for (int i = 0; i < restrn__.ngfix; i++) {
      OpenMM_IntArray* bondGroups = OpenMM_IntArray_create (0);
      OpenMM_IntArray_append (bondGroups, restrn__.igfix[2*i] - 1);
      OpenMM_IntArray_append (bondGroups, restrn__.igfix[2*i + 1] - 1);

      OpenMM_DoubleArray* bondParameters = OpenMM_DoubleArray_create (0);
      OpenMM_DoubleArray_append (bondParameters, restrn__.gfix[3*i]
                                    *convert);
      OpenMM_DoubleArray_append (bondParameters, restrn__.gfix[3*i + 1]
                                    *OpenMM_NmPerAngstrom);
      OpenMM_DoubleArray_append (bondParameters, restrn__.gfix[3*i + 2]
                                    *OpenMM_NmPerAngstrom);

      OpenMM_CustomCentroidBondForce_addBond (force, bondGroups,
                                              bondParameters);
      OpenMM_IntArray_destroy (bondGroups);
      OpenMM_DoubleArray_destroy (bondParameters);
   }
   OpenMM_System_addForce(system, (OpenMM_Force*) force);
}

/*
 *    ############################################################
 *          Create Standalone Thermotat and Barostat Methods
 *    ############################################################
 */

static void setupAndersenThermostat (OpenMM_System* system, FILE* log) {

   OpenMM_AndersenThermostat* andersenThermostat;

   andersenThermostat = OpenMM_AndersenThermostat_create (bath__.kelvin,
                                                          bath__.collide);
   OpenMM_System_addForce (system, (OpenMM_Force*) andersenThermostat);

   if (log) {
      (void) fprintf (log, "\n Andersen Thermostat:\n" );
      (void) fprintf (log, "\n Target Temperature   %15.2f K",
                      OpenMM_AndersenThermostat_getDefaultTemperature
                      (andersenThermostat));
      (void) fprintf (log, "\n Collision Frequency  %15.7e ps^(-1)",
                      OpenMM_AndersenThermostat_getDefaultCollisionFrequency
                      (andersenThermostat));
      (void) fprintf (log, "\n Random Number Seed                 %d\n",
                      OpenMM_AndersenThermostat_getRandomNumberSeed
                      (andersenThermostat));
   }
}

static void setupMonteCarloBarostat (OpenMM_System* system, FILE* log) {

   OpenMM_MonteCarloBarostat* monteCarloBarostat;

   int frequency = bath__.voltrial;
   monteCarloBarostat = OpenMM_MonteCarloBarostat_create
                            (bath__.atmsph*1.01325, bath__.kelvin, frequency);
   OpenMM_System_addForce (system, (OpenMM_Force*) monteCarloBarostat);

   if (log) {
      (void) fprintf (log, "\n MonteCarlo Barostat :\n");
      (void) fprintf (log, "\n Target Temperature   %15.2f K",
                      OpenMM_MonteCarloBarostat_getDefaultTemperature
                      (monteCarloBarostat));

      //
      // Change needed for latest Stanford OpenMM update (M. Harger)
      //
      // (void) fprintf (log, "\n Target Temperature   %15.2f K",
      //                 OpenMM_MonteCarloBarostat_getDefaultTemperature
      //                 (monteCarloBarostat));

      (void) fprintf (log, "\n Target Pressure      %15.4f atm",
                      OpenMM_MonteCarloBarostat_getDefaultPressure
                      (monteCarloBarostat) / 1.01325);
      (void) fprintf (log, "\n Trial Frequency                   %d",
                      OpenMM_MonteCarloBarostat_getFrequency
                      (monteCarloBarostat));
      (void) fprintf (log, "\n Random Number Seed                 %d\n",
                      OpenMM_MonteCarloBarostat_getRandomNumberSeed
                      (monteCarloBarostat));
   }
}

/*
 *    ############################################################
 *         Setup Positions, Velocities and Rattle Constraints
 *    ############################################################
 */

static void setupPositions (OpenMM_Vec3Array* initialPosInNm, FILE* log) {

   int ii;
   for (ii = 0; ii < atoms__.n; ii++) {
      OpenMM_Vec3 posInNm;
      posInNm.x = atoms__.x[ii]*OpenMM_NmPerAngstrom;
      posInNm.y = atoms__.y[ii]*OpenMM_NmPerAngstrom;
      posInNm.z = atoms__.z[ii]*OpenMM_NmPerAngstrom;
      OpenMM_Vec3Array_append (initialPosInNm, posInNm);
   }
}

static void setupVelocities (OpenMM_Vec3Array* initialVelInNm, FILE* log) {

   int ii;
   for (ii = 0; ii < atoms__.n; ii++) {
      OpenMM_Vec3 velInNm;
      int offset;
      offset = 3*ii;
      velInNm.x = moldyn__.v[offset]*OpenMM_NmPerAngstrom;
      velInNm.y = moldyn__.v[offset+1]*OpenMM_NmPerAngstrom;
      velInNm.z = moldyn__.v[offset+2]*OpenMM_NmPerAngstrom;
      OpenMM_Vec3Array_append (initialVelInNm, velInNm);
   }
}

static void setupConstraints (OpenMM_System* system, FILE* log) {

   int ii;
   for (ii = 0; ii < freeze__.nrat; ii++) {
      OpenMM_System_addConstraint (system, *(freeze__.irat+2*ii) -1,
                                   *(freeze__.irat + 2*ii +1)-1,
                          (*(freeze__.krat +ii))*OpenMM_NmPerAngstrom);
   }
}

/*
 *    ############################################################
 *           Set Calculation Platform to Reference or CUDA
 *    ############################################################
 */

#include "gpu_cards.h"

static OpenMM_Platform* getReferencePlatform (FILE* log) {

   OpenMM_Platform* platform = OpenMM_Platform_getPlatformByName ("Reference");
   if (platform == NULL) {
      if (log) {
         (void) fprintf (log, "Reference Platform Unavailable\n");
      }
      return platform;
   }
   return platform;
}

static OpenMM_Platform* getCUDAPlatform (FILE* log) {

   char buffer[8];
   char* deviceId;
   int device_number;
   bool device_key = false;
   OpenMM_Platform* platform = OpenMM_Platform_getPlatformByName ("CUDA");
   if (platform == NULL) {
      if (log) {
         (void) fprintf (log, "\n CUDA Platform Unavailable\n");
      }
      return platform;
   }

   if (cudaDevices[0] != '\0') {
      deviceId = cudaDevices;
      device_key = true;
   } else {
      device_number = findBestCUDACard();
      if (device_number < 0) {
         deviceId = NULL;
      } else {
         deviceId = buffer;
         sprintf(deviceId, "%d", device_number);
      }
   }

   //deviceId = "0,1";
   //fprintf (log, "\n Value of deviceId is %s \n", deviceId);

   if (device_key) {
      OpenMM_Platform_setPropertyDefaultValue (platform, "CudaDeviceIndex",
                                               deviceId);
      if (log) {
         (void) fprintf (log, "\n Platform CUDA :  Setting Device ID to %s from CUDA-DEVICE keyword\n", deviceId);
      }
   } else if (log && inform__.verbose) {
      (void) fprintf (log, "\n Platform CUDA :  Setting Device ID to %s \n", deviceId);
   }

   OpenMM_Platform_setPropertyDefaultValue (platform, "Precision", "mixed" );

   return platform;
}

/*
 *    ############################################################
 *               Initialize the OpenMM Data Structures
 *    ############################################################
 *
 *    Following actions are performed here:
 *    (1) Load any available OpenMM plugins
 *    (2) Allocate a OpenMMData structure to hang on to OpenMM data
 *        structures in a manner which is opaque to the caller
 *    (3) Initialize the OpenMM::System with the AMOEBA force field
 *        objects
 *    (4) Create an Integrator and a Context associating the
 *        Integrator with the System
 *    (5) Select the OpenMM platform to be used.
 *    (6) Return an opaque pointer to the OpenMMData struct
 */

extern "C" {

extern void fatal_();

void openmm_init_ (void** ommHandle, double* dt) {

   int ii;
   int mdMode = 0;
   int removeConstrainedCovalentIxns = 0;
   char buffer[128];
   char buffer2[128];
   FILE* log = stderr;

   // Check if we can run the simulation with OpenMM
   if (boxes__.octahedron) {
      printf("\n Truncated octahedral boxes are not supported in OpenMM.\n");
      fatal_(); //This will never return
   }

   // allocate space for opaque handle to hold OpenMM objects
   // such as system, integrator, context, etc.
   OpenMMData* omm = (OpenMMData*) malloc(sizeof(struct OpenMMData_s));

   // These are temporary OpenMM objects used and discarded here

   OpenMM_Vec3Array*       initialPosInNm;
   OpenMM_Vec3Array*       initialVelInNm;
   OpenMM_StringArray*     pluginList;
   OpenMM_Platform*        platform;

   // Load all OpenMM plugin libraries from their default location;
   // Call the plugin loading routine twice to fix an issue with MacOS
   // where the first library in the alphabetical list gets skipped

   pluginList = OpenMM_Platform_loadPluginsFromDirectory
                    (OpenMM_Platform_getDefaultPluginsDirectory());
   pluginList = OpenMM_Platform_loadPluginsFromDirectory
                    (OpenMM_Platform_getDefaultPluginsDirectory());

   if (inform__.verbose) {
      (void) fprintf (stderr, "\n Default OpenMM Plugin Directory :  %s\n\n",
                          OpenMM_Platform_getDefaultPluginsDirectory());
      for (ii = 0; ii < OpenMM_StringArray_getSize(pluginList); ii++) {
         (void) fprintf (stderr, " Plugin Library :  %s\n",
                             OpenMM_StringArray_get(pluginList, ii));
      }
   }

   OpenMM_StringArray_destroy (pluginList);
   (void) fflush (NULL);

   // create System and Force objects within the System; retain a reference
   // to each force object so we can fill in the forces; mote the OpenMM
   // System takes ownership of the force objects, don't delete them yourself

   omm->system = OpenMM_System_create ();
   setupSystemParticles (omm->system, log);
   setupCMMotionRemover (omm->system, log);

   if (potent__.use_bond) {
      setupAmoebaBondForce (omm->system,
                            removeConstrainedCovalentIxns, log);
   }

   if (potent__.use_angle) {
      setupAmoebaAngleForce (omm->system,
                             removeConstrainedCovalentIxns, log);
      setupAmoebaInPlaneAngleForce (omm->system,
                                    removeConstrainedCovalentIxns, log);
   }

   if (potent__.use_strbnd) {
      setupAmoebaStretchBendForce (omm->system,
                                   removeConstrainedCovalentIxns, log);
   }

   if (potent__.use_urey) {
      setupAmoebaUreyBradleyForce (omm->system,
                                   removeConstrainedCovalentIxns, log);
   }

   if (potent__.use_opbend) {
      setupAmoebaOutOfPlaneBendForce (omm->system, log);
   }

   if (potent__.use_improp) {
      setupAmoebaImproperTorsionForce (omm->system, log);
   }

   if (potent__.use_tors) {
      setupAmoebaTorsionForce (omm->system, log);
   }

   if (potent__.use_strtor) {
      setupAmoebaStretchTorsionForce (omm->system, log);
   }

   if (potent__.use_angtor) {
      setupAmoebaAngleTorsionForce (omm->system, log);
   }

   if (potent__.use_pitors) {
      setupAmoebaPiTorsionForce (omm->system, log);
   }

   if (potent__.use_tortor) {
      setupAmoebaTorsionTorsionForce (omm->system, log);
   }

   if (potent__.use_vdw) {
      setupAmoebaVdwForce (omm->system, log);
   }

   if (potent__.use_ct) {
      setupAmoebaCTForce (omm->system, log);
   }

   if (potent__.use_charge) {
      setupAmoebaChargeForce (omm->system, log);
   }

   if (potent__.use_mpole) {
      setupAmoebaMultipoleForce (omm->system, log);
      if (potent__.use_solv) {
         setupAmoebaGeneralizedKirkwoodForce (omm->system, 1, log);
      }
   }

   if (potent__.use_solv) {
      setupAmoebaWcaDispersionForce (omm->system, log);
   }

   if (bath__.isobaric && bath__.atmsph > 0.0) {
      mdMode = 4;
      setupMonteCarloBarostat (omm->system, log);
   }

   // setup of constraints, restraints, positions and velocities
   setupConstraints (omm->system, log);

   setupTorsionRestraints (omm->system, log);
   setupDistanceRestraints (omm->system, log);
   setupPositionalRestraints (omm->system, log);
   (void) fflush (NULL);
   setupAngleRestraints (omm->system, log);
   setupCentroidRestraints (omm->system, log);

   initialPosInNm = OpenMM_Vec3Array_create (0);
   setupPositions (initialPosInNm, log);

   initialVelInNm = OpenMM_Vec3Array_create (0);
   setupVelocities (initialVelInNm, log);

   // choose an Integrator, and a Context connecting the System with the
   // Integrator; let the Context choose the best available Platform;
   // initialize the configuration from default positions collected above

   setNullTerminator (mdstuf__.integrate, 11, buffer);
   setNullTerminator (bath__.thermostat, 11, buffer2);

   if (strncasecmp (buffer, "RESPA", 5) == 0) {
      omm->integrator = (OpenMM_Integrator*) OpenMM_CustomIntegrator_create (*dt);
      OpenMM_CustomIntegrator* integrator = (OpenMM_CustomIntegrator*) omm->integrator;
      OpenMM_CustomIntegrator_addUpdateContextState (integrator);

      int n = int (round ((*dt) / 0.00025));
      char n_char[16] = {0};
      sprintf (n_char, "%d", n);

      char e1[64] = {0};
      strcat (e1, "v+0.5*(dt/");
      strcat (e1, n_char);
      strcat (e1, ")*f0/m");

      char e11[64] = {0};
      strcat (e11, n_char);
      strcat (e11, "*(x-x1)/dt+");
      strcat (e11, e1);

      char e2[64] = {0};
      strcat (e2, "x+(dt/");
      strcat (e2, n_char);
      strcat (e2, ")*v");

      OpenMM_CustomIntegrator_addPerDofVariable (integrator, "x1", 0.0);

      OpenMM_CustomIntegrator_addComputePerDof (integrator, "v", "v+0.5*dt*f1/m");
      for (int i = 0; i < n; i++) {
         OpenMM_CustomIntegrator_addComputePerDof (integrator, "v", e1);
         OpenMM_CustomIntegrator_addComputePerDof (integrator, "x", e2);
         OpenMM_CustomIntegrator_addComputePerDof (integrator, "x1", "x");
         OpenMM_CustomIntegrator_addConstrainPositions (integrator);
         OpenMM_CustomIntegrator_addComputePerDof (integrator, "v", e11);
         OpenMM_CustomIntegrator_addConstrainVelocities (integrator);
      }
      OpenMM_CustomIntegrator_addComputePerDof (integrator, "v", "v+0.5*dt*f1/m");
      OpenMM_CustomIntegrator_addConstrainVelocities (integrator);

      // setup and computation for the Bussi-Parrinello thermostat

      int n_constraints = OpenMM_System_getNumConstraints (omm->system) + 3;
      int adjlen = (n_constraints + 2) * 16;
      char* adjustment = (char*) malloc (adjlen * sizeof(char));
      char* e3 = (char*) malloc (adjlen * sizeof(char));
      for (int i = 0; i < adjlen; ++i) {
         adjustment[i] = 0;
         e3[i] = 0;
      }

      for (int i = 0; i < n_constraints; i++) {
         strcat (adjustment, "gaussian^2+");
      }
      strcat (adjustment, "gaussian^2");

      strcat (e3, "s-(");
      strcat (e3, adjustment);
      strcat (e3, ")");

      if (strncasecmp (buffer2, "BUSSI", 5) == 0) {
         OpenMM_CustomIntegrator_addGlobalVariable (integrator, "df", mdstuf__.nfree);
         OpenMM_CustomIntegrator_addGlobalVariable (integrator, "gasconst",
             units__.gasconst * units__.joule);
         OpenMM_CustomIntegrator_addGlobalVariable (integrator, "tautemp", bath__.tautemp);
         OpenMM_CustomIntegrator_addGlobalVariable (integrator, "kelvin", bath__.kelvin);
         OpenMM_CustomIntegrator_addGlobalVariable (integrator, "ke", 0.0);
         OpenMM_CustomIntegrator_addGlobalVariable (integrator, "temp", 0.0);
         OpenMM_CustomIntegrator_addGlobalVariable (integrator, "scale", 0.0);
         OpenMM_CustomIntegrator_addGlobalVariable (integrator, "c", exp (-(*dt)/bath__.tautemp));
         OpenMM_CustomIntegrator_addGlobalVariable (integrator, "d", 0.0);
         OpenMM_CustomIntegrator_addGlobalVariable (integrator, "r", 0.0);
         OpenMM_CustomIntegrator_addGlobalVariable (integrator, "s", 0.0);
         OpenMM_CustomIntegrator_addPerDofVariable (integrator, "si", 0.0);

         OpenMM_CustomIntegrator_addComputeSum (integrator, "ke", "m*v*v/2.0");
         OpenMM_CustomIntegrator_addComputeGlobal (integrator, "temp",
             "2.0*ke/(df*gasconst)");
         OpenMM_CustomIntegrator_addComputeGlobal (integrator, "temp",
             "temp + step(0.0000000001-temp) * 0.1");
         OpenMM_CustomIntegrator_addComputeGlobal (integrator, "c", "exp(-dt/tautemp)");
         OpenMM_CustomIntegrator_addComputeGlobal (integrator, "d", "(1.0-c)*(kelvin/temp)/df");
         OpenMM_CustomIntegrator_addComputeGlobal (integrator, "r", "gaussian");
         OpenMM_CustomIntegrator_addComputePerDof (integrator, "si", "gaussian");
         OpenMM_CustomIntegrator_addComputeSum (integrator, "s", "si*si");
         OpenMM_CustomIntegrator_addComputeGlobal (integrator, "s", e3);
         OpenMM_CustomIntegrator_addComputeGlobal (integrator, "scale",
              "sqrt(c+(s+r*r)*d+2.0*r*sqrt(c*d))");
         OpenMM_CustomIntegrator_addComputeGlobal (integrator, "scale",
             "scale-2.0*scale*step(-r-sqrt(c/d))");
         OpenMM_CustomIntegrator_addComputePerDof (integrator, "v", "scale*v");

      } else if (strncasecmp (buffer2, "ANDERSEN", 8) == 0) {
         setupAndersenThermostat (omm->system, log);
      }

      free (adjustment);
      free (e3);

   } else if (strncasecmp (buffer, "VERLET", 6) == 0) {
      omm->integrator = (OpenMM_Integrator*)OpenMM_VerletIntegrator_create (*dt);
      if (strncasecmp (buffer2, "ANDERSEN", 8) == 0) {
         setupAndersenThermostat (omm->system, log);
      }

   } else if (strncasecmp (buffer, "STOCHASTIC", 10) == 0) {

      if (log) {
         (void) fprintf (log, "\n Stochastic Integrator:\n");
         (void) fprintf (log, "\n Temperature          %15.4f K",
                             bath__.kelvin );
         (void) fprintf (log, "\n Friction             %15.4f ps^(-1)",
                             stodyn__.friction );
         (void) fprintf (log, "\n TimeStep             %15.4f ps\n",
                             *dt );
         (void) fflush (log);
      }
      omm->integrator = (OpenMM_Integrator*)OpenMM_LangevinIntegrator_create
                             (bath__.kelvin, stodyn__.friction, *dt);

   } else {
      (void) fprintf (stderr, "\n Integrator %s is Not Supported\n", buffer);
      (void) fflush (stderr);
      exit (-1);
   }

   // choose either the reference or the CUDA platform

   // platform = getReferencePlatform (log);
   platform = getCUDAPlatform (log);
   if (platform == NULL) {
      exit (-1);
   }

   omm->context = OpenMM_Context_create_2 (omm->system, omm->integrator,
                                           platform);

   if (inform__.debug) {
      (void) fprintf (log, "\n OpenMMDataHandle:  %x\n", (void*)(omm));
      (void) fprintf (log, "\n Integrator:  %x\n", (void*)(omm->integrator));
   }

   OpenMM_Context_setPositions (omm->context, initialPosInNm);
   OpenMM_Context_setVelocities (omm->context, initialVelInNm);

   if (inform__.verbose) {
      int arraySz;
      int maxPrint;
      double x1, x2, x3;
      double v1, v2, v3;
      arraySz = OpenMM_Vec3Array_getSize (initialPosInNm);
      maxPrint = 5;
      (void) fprintf (log, "\n Initial Positions and Velocities :\n\n");
      for (ii = 0; ii < arraySz; ii++) {
         x1 = (*OpenMM_Vec3Array_get(initialPosInNm, ii)).x
                   / OpenMM_NmPerAngstrom;
         x2 = (*OpenMM_Vec3Array_get(initialPosInNm, ii)).y
                   / OpenMM_NmPerAngstrom;
         x3 = (*OpenMM_Vec3Array_get(initialPosInNm, ii)).z
                   / OpenMM_NmPerAngstrom;
         v1 = (*OpenMM_Vec3Array_get(initialVelInNm, ii)).x
                   / OpenMM_NmPerAngstrom;
         v2 = (*OpenMM_Vec3Array_get(initialVelInNm, ii)).y
                   / OpenMM_NmPerAngstrom;
         v3 = (*OpenMM_Vec3Array_get(initialVelInNm, ii)).z
                  / OpenMM_NmPerAngstrom;
         (void) fprintf (log, "%7d   POS    %17.8e %17.8e %17.8e\n",
                         ii+1, x1, x2, x3);
         (void) fprintf (log, "%7d   VEL    %17.8e %17.8e %17.8e\n",
                         ii+1, v1, v2, v3);
         if (ii == maxPrint-1 && ii < (arraySz-maxPrint-1))
             ii = arraySz - maxPrint - 1;
      }
   }

   *ommHandle = (void*) omm;
}

/*
 *    ############################################################
 *          Copy State from OpenMM to TINKER Data Structures
 *    ############################################################
 *
 *    @param omm                    handle with OpenMM data structures
 *    @param timeInPs               output simulation time in ps
 *    @param kineticEnergyInKcal    output kinetic energy in kcal
 *    @param potentialEnergyInKcal  output potential energy in kcal
 */

void openmm_copy_state_ (void** omm, double *timeInPs,
                         double* kineticEnergyInKcal,
                         double* potentialEnergyInKcal) {

   OpenMM_State* state;
   const OpenMM_Vec3Array* positionArray;
   const OpenMM_Vec3Array* velocityArray;
   const OpenMM_Vec3Array* forceArray;
   OpenMM_Vec3 aBox;
   OpenMM_Vec3 bBox;
   OpenMM_Vec3 cBox;
   int infoMask;
   int ii;
   int offset;
   double amass;
   double positionConvert;
   double velocityConvert;
   double forceConvert;
   OpenMMData* openMMDataHandle;

   openMMDataHandle = (OpenMMData*) (*omm);

   infoMask = OpenMM_State_Positions;
   infoMask += OpenMM_State_Velocities;
   infoMask += OpenMM_State_Forces;
   infoMask += OpenMM_State_Energy;

   // State object is created here and must be explicitly destroyed below

   state = OpenMM_Context_getState (openMMDataHandle->context, infoMask, 0);
   *timeInPs = OpenMM_State_getTime (state);

   OpenMM_State_getPeriodicBoxVectors (state, &aBox, &bBox, &cBox);

   *(boxes__.xbox) = aBox.x / OpenMM_NmPerAngstrom;
   *(boxes__.ybox) = bBox.y / OpenMM_NmPerAngstrom;
   *(boxes__.zbox) = cBox.z / OpenMM_NmPerAngstrom;
   *(boxes__.xbox2) = 0.5 * (*(boxes__.xbox));
   *(boxes__.ybox2) = 0.5 * (*(boxes__.ybox));
   *(boxes__.zbox2) = 0.5 * (*(boxes__.zbox));

   // positions, velocities and forces are maintained as Vec3Arrays
   // inside the State; this provides access, but don't destroy them
   // yourself as they will go away with the State

   positionConvert = 1.0 / OpenMM_NmPerAngstrom;
   velocityConvert = 1.0 / OpenMM_NmPerAngstrom;
   forceConvert = 10.0;

   positionArray = OpenMM_State_getPositions (state);
   velocityArray = OpenMM_State_getVelocities (state);
   forceArray = OpenMM_State_getForces (state);

   for (ii = 0; ii < atoms__.n; ii++) {
      atoms__.x[ii] = (*OpenMM_Vec3Array_get(positionArray, ii)).x
                           * positionConvert;
      atoms__.y[ii] = (*OpenMM_Vec3Array_get(positionArray, ii)).y
                           * positionConvert;
      atoms__.z[ii] = (*OpenMM_Vec3Array_get(positionArray, ii)).z
                           * positionConvert;
   }

   for (ii = 0; ii < atoms__.n; ii++) {
      offset = 3*ii;
      moldyn__.v[offset] = (*OpenMM_Vec3Array_get(velocityArray, ii)).x
                                * velocityConvert;
      moldyn__.v[offset+1] = (*OpenMM_Vec3Array_get(velocityArray, ii)).y
                                  * velocityConvert;
      moldyn__.v[offset+2] = (*OpenMM_Vec3Array_get(velocityArray, ii)).z
                                  * velocityConvert;
   }

   for (ii = 0; ii < atoms__.n; ii++) {
      offset = 3*ii;
      amass = 1.0 / atomid__.mass[ii];
      moldyn__.a[offset] = (*OpenMM_Vec3Array_get(forceArray, ii)).x
                                  * amass * forceConvert;
      moldyn__.a[offset+1] = (*OpenMM_Vec3Array_get(forceArray, ii)).y
                                  * amass * forceConvert;
      moldyn__.a[offset+2] = (*OpenMM_Vec3Array_get(forceArray, ii)).z
                                  * amass * forceConvert;
   }

   for (ii = 0; ii < atoms__.n; ii++) {
      offset = 3*ii;
      moldyn__.aalt[offset] = 0.0;
      moldyn__.aalt[offset+1] = 0.0;
      moldyn__.aalt[offset+2] = 0.0;
   }

   if (inform__.debug) {

      (void) fprintf (stderr, "\n State: E=%15.7e [%15.7e %15.7e]\n",
                      (*kineticEnergyInKcal + *potentialEnergyInKcal),
                      *kineticEnergyInKcal, *potentialEnergyInKcal);
      (void) fprintf (stderr, "        t=%15.7e ps\n\n", *timeInPs);

      for (ii = 0; ii < atoms__.n; ii++) {
         (void) fprintf (stderr, "%7d   POS   %17.8e %17.8e %17.8e\n", ii+1,
             (*OpenMM_Vec3Array_get(positionArray,ii)).x*positionConvert,
             (*OpenMM_Vec3Array_get(positionArray,ii)).y*positionConvert,
             (*OpenMM_Vec3Array_get(positionArray,ii)).z*positionConvert);
         (void) fprintf (stderr, "%7d   VEL   %17.8e %17.8e %17.8e\n", ii+1,
             (*OpenMM_Vec3Array_get(velocityArray,ii)).x*velocityConvert,
             (*OpenMM_Vec3Array_get(velocityArray,ii)).y*velocityConvert,
             (*OpenMM_Vec3Array_get(velocityArray,ii)).z*velocityConvert);
         (void) fprintf (stderr, "%7d   FRC   %17.8e %17.8e %17.8e\n", ii+1,
             (*OpenMM_Vec3Array_get(forceArray,ii)).x*forceConvert,
             (*OpenMM_Vec3Array_get(forceArray,ii)).y*forceConvert,
             (*OpenMM_Vec3Array_get(forceArray,ii)).z*forceConvert);
      }
      (void) fflush (stderr);
   }

   *kineticEnergyInKcal = OpenMM_State_getKineticEnergy (state)
                              * OpenMM_KcalPerKJ;
   *potentialEnergyInKcal = OpenMM_State_getPotentialEnergy (state)
                                * OpenMM_KcalPerKJ;

   OpenMM_State_destroy (state);
}

/*
 *    ############################################################
 *       Update TINKER Data Structures; Call mdstat and mdsave
 *    ############################################################
 *
 *    @param omm          handle containing OpenMM data structures
 *    @param dt           simulation time step in ps
 *    @param istep        current step in MD loop
 *    @param callMdStat   if nonzero, call TINKER mdstat routine
 *    @param callMdSave   if nonzero, call TINKER mdsave routine
 */

void openmm_update_ (void** omm, double* dt, int* istep,
                     int* callMdStat, int* callMdSave) {

   OpenMM_State* state;
   const OpenMM_Vec3Array* positionArray;
   const OpenMM_Vec3Array* velocityArray;
   const OpenMM_Vec3Array* forceArray;
   OpenMM_Vec3 aBox;
   OpenMM_Vec3 bBox;
   OpenMM_Vec3 cBox;
   int infoMask;
   int ii;
   double amass;
   double positionConvert;
   double velocityConvert;
   double forceConvert;

   static const double gasconst = units__.gasconst;

   double totalEnergy, potentialEnergy, kineticEnergy;

   OpenMMData* openMMDataHandle;

   openMMDataHandle = (OpenMMData*) (*omm);

   infoMask = OpenMM_State_Positions;
   infoMask += OpenMM_State_Velocities;
   infoMask += OpenMM_State_Forces;
   infoMask += OpenMM_State_Energy;

   // state object is created here and must be destroyed below

   state = OpenMM_Context_getState (openMMDataHandle->context, infoMask, 0);
   OpenMM_State_getPeriodicBoxVectors (state, &aBox, &bBox, &cBox);

   *(boxes__.xbox) = aBox.x / OpenMM_NmPerAngstrom;
   *(boxes__.ybox) = bBox.y / OpenMM_NmPerAngstrom;
   *(boxes__.zbox) = cBox.z / OpenMM_NmPerAngstrom;
   *(boxes__.xbox2) = 0.5 * (*(boxes__.xbox));
   *(boxes__.ybox2) = 0.5 * (*(boxes__.ybox));
   *(boxes__.zbox2) = 0.5 * (*(boxes__.zbox));

   // fprintf (stderr, "openmm_update_ %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e \n",
   // *(boxes__.xbox), *(boxes__.ybox), *(boxes__.zbox),
   // aBox.x, bBox.y, cBox.z );

   // load positions/velocities and energies

   positionConvert = 1.0 / OpenMM_NmPerAngstrom;
   velocityConvert = 1.0 / OpenMM_NmPerAngstrom;
   forceConvert = 10.0;

   positionArray = OpenMM_State_getPositions (state);
   velocityArray = OpenMM_State_getVelocities (state);
   forceArray = OpenMM_State_getForces (state);

   for (ii = 0; ii < atoms__.n; ii++) {
      atoms__.x[ii] = (*OpenMM_Vec3Array_get(positionArray, ii)).x
                          * positionConvert;
      atoms__.y[ii] = (*OpenMM_Vec3Array_get(positionArray, ii)).y
                          * positionConvert;
      atoms__.z[ii] = (*OpenMM_Vec3Array_get(positionArray, ii)).z
                          * positionConvert;
   }

   for (ii = 0; ii < atoms__.n; ii++) {
      int offset = ii*3;
      moldyn__.v[offset] = (*OpenMM_Vec3Array_get(velocityArray, ii)).x
                                * velocityConvert;
      moldyn__.v[offset+1] = (*OpenMM_Vec3Array_get(velocityArray, ii)).y
                                  * velocityConvert;
      moldyn__.v[offset+2] = (*OpenMM_Vec3Array_get(velocityArray, ii)).z
                                  * velocityConvert;
   }

   for (ii = 0; ii < atoms__.n; ii++) {
      int offset = ii*3;
      amass = 1.0 / atomid__.mass[ii];
      moldyn__.a[offset] = (*OpenMM_Vec3Array_get(forceArray, ii)).x
                                * amass * forceConvert;
      moldyn__.a[offset+1] = (*OpenMM_Vec3Array_get(forceArray, ii)).y
                                  * amass * forceConvert;
      moldyn__.a[offset+2] = (*OpenMM_Vec3Array_get(forceArray, ii)).z
                                  * amass * forceConvert;
   }

   for (ii = 0; ii < atoms__.n; ii++) {
      int offset = ii*3;
      moldyn__.aalt[offset] = 0.0;
      moldyn__.aalt[offset+1] = 0.0;
      moldyn__.aalt[offset+2] = 0.0;
   }

   potentialEnergy = (OpenMM_State_getPotentialEnergy(state))
                         * OpenMM_KcalPerKJ;
   kineticEnergy = (OpenMM_State_getKineticEnergy(state))
                       * OpenMM_KcalPerKJ;
   totalEnergy = potentialEnergy + kineticEnergy;

   if (inform__.debug) {

      double eksum, temp, pres;
      double ekin[3][3];

      kinetic_ (&eksum, &ekin, &temp);

      (void) fprintf (stderr, "\n State: E=%15.7e [%15.7e %15.7e]\n",
                      totalEnergy, potentialEnergy, kineticEnergy);
      (void) fprintf (stderr, "        t=%15.7e ps T=%12.5e InitT=%12.3f\n\n",
                      (*dt)*(*istep), temp, bath__.kelvin );

      for (ii = 0; ii < atoms__.n; ii++) {
         (void) fprintf (stderr, "%7d   POS   %17.8e %17.8e %17.8e\n", ii+1,
             (*OpenMM_Vec3Array_get(positionArray,ii)).x*positionConvert,
             (*OpenMM_Vec3Array_get(positionArray,ii)).y*positionConvert,
             (*OpenMM_Vec3Array_get(positionArray,ii)).z*positionConvert);
         (void) fprintf (stderr, "%7d   VEL   %17.8e %17.8e %17.8e\n", ii+1,
             (*OpenMM_Vec3Array_get(velocityArray,ii)).x*velocityConvert,
             (*OpenMM_Vec3Array_get(velocityArray,ii)).y*velocityConvert,
             (*OpenMM_Vec3Array_get(velocityArray,ii)).z*velocityConvert);
         (void) fprintf (stderr, "%7d   FRC   %17.8e %17.8e %17.8e\n", ii+1,
             (*OpenMM_Vec3Array_get(forceArray,ii)).x*forceConvert,
             (*OpenMM_Vec3Array_get(forceArray,ii)).y*forceConvert,
             (*OpenMM_Vec3Array_get(forceArray,ii)).z*forceConvert);
      }
      (void) fflush (stderr);
   }


   // make calls to mdstat and/or mdsave if flags are set

   if (*callMdStat || *callMdSave) {
      double eksum, temp, pres;
      double ekin[3][3];
      kinetic_ (&eksum, &ekin, &temp);

      if (*callMdStat) {
         pres = 0.0;
         if (bath__.isobaric) {
            lattice_ ();
         }
         mdstat_ (istep,dt,&totalEnergy,&potentialEnergy,&eksum,&temp,&pres);
      }
      if (*callMdSave) {
         bounds_ ();
         mdsave_ (istep,dt,&potentialEnergy,&eksum);
      }
   }
   OpenMM_State_destroy (state);
}

void openmm_get_energies_ (void** omm, double* ke, double* pe) {

   OpenMM_State* state;
   const OpenMM_Vec3Array* positionArray;
   OpenMM_Vec3 aBox;
   OpenMM_Vec3 bBox;
   OpenMM_Vec3 cBox;
   int infoMask;
   int ii;

   static const double gasconst = 1.98720415e-3;

   double totalEnergy, potentialEnergy, kineticEnergy;
   OpenMMData* openMMDataHandle;

   openMMDataHandle = (OpenMMData*) (*omm);

   infoMask = OpenMM_State_Positions;
   infoMask += OpenMM_State_Velocities;
   if (inform__.debug) {
      infoMask += OpenMM_State_Forces;
   }
   infoMask += OpenMM_State_Energy;

   // State object is created here and must be explicitly destroyed below

   state = OpenMM_Context_getState (openMMDataHandle->context, infoMask, 0);
   OpenMM_State_getPeriodicBoxVectors (state, &aBox, &bBox, &cBox);

   // load positions, velocities and energies

   //positionArray = OpenMM_State_getPositions (state);
   //velocityArray = OpenMM_State_getVelocities (state);

   potentialEnergy = (OpenMM_State_getPotentialEnergy (state))
                         * OpenMM_KcalPerKJ;
   kineticEnergy = (OpenMM_State_getKineticEnergy (state)) * OpenMM_KcalPerKJ;
   totalEnergy = potentialEnergy + kineticEnergy;

   *ke = kineticEnergy;
   *pe = potentialEnergy;

   OpenMM_State_destroy (state);
}

void openmm_update_box_ (void** omm) {

   OpenMM_Vec3 aBox;
   OpenMM_Vec3 bBox;
   OpenMM_Vec3 cBox;

   OpenMMData* openMMDataHandle;
   openMMDataHandle = (OpenMMData*) (*omm);

   aBox.x = 0.0;
   aBox.y = 0.0;
   aBox.z = 0.0;
   bBox.x = 0.0;
   bBox.y = 0.0;
   bBox.z = 0.0;
   cBox.x = 0.0;
   cBox.y = 0.0;
   cBox.z = 0.0;

   aBox.x = *(boxes__.xbox) * OpenMM_NmPerAngstrom;
   bBox.y = *(boxes__.ybox) * OpenMM_NmPerAngstrom;
   cBox.z = *(boxes__.zbox) * OpenMM_NmPerAngstrom;
   OpenMM_Context_setPeriodicBoxVectors (openMMDataHandle->context,
                                         &aBox, &bBox, &cBox);
}

void openmm_update_positions_ (void** omm) {

   OpenMMData* openMMDataHandle;
   OpenMM_Vec3Array* initialPosInNm;
   FILE* log = stderr;

   openMMDataHandle = (OpenMMData*) (*omm);

   initialPosInNm = OpenMM_Vec3Array_create (0);
   setupPositions (initialPosInNm, log);
   OpenMM_Context_setPositions (openMMDataHandle->context, initialPosInNm);
   OpenMM_Vec3Array_destroy (initialPosInNm);
}

void openmm_take_steps_ (void** omm, int* numSteps) {

   OpenMMData* openMMDataHandle = (OpenMMData*) (*omm);
   OpenMM_Integrator_step (openMMDataHandle->integrator, *numSteps);
}

void openmm_cleanup_ (void** omm) {

   // clean up top-level heap allocated objects we are done with

   OpenMMData* openMMDataHandle = (OpenMMData*) (*omm);
   OpenMM_Context_destroy (openMMDataHandle->context);
   OpenMM_Integrator_destroy (openMMDataHandle->integrator);
   OpenMM_System_destroy (openMMDataHandle->system);
   free (openMMDataHandle);
}
}

static void zeroTinkerForce (double* tinkerForce) {

   int ii;
   for (ii = 0; ii < atoms__.n; ii++){
      *(tinkerForce + 3*ii + 0) = 0.0;
      *(tinkerForce + 3*ii + 1) = 0.0;
      *(tinkerForce + 3*ii + 2) = 0.0;
   }
}

static void zeroVec3Force (OpenMM_Vec3Array* tinkerForce) {

   int ii;
   for (ii = 0; ii < atoms__.n; ii++) {
      OpenMM_Vec3 force;
      force.x = force.y = force.z = 0.0;
      OpenMM_Vec3Array_set (tinkerForce, ii, force);
   }
}

static void setTinker1DArray (int size, double* tinkerArray, double value) {

   int ii;
   for (ii = 0; ii < size; ii++) {
      tinkerArray[ii] = value;
   }
}

static void loadTinkerForce (double* tinkerForce, int add,
                             OpenMM_Vec3Array* arrayToLoad) {

   int ii;
   if (add == 0) {
      for (ii = 0; ii < atoms__.n; ii++) {
         OpenMM_Vec3 force;
         force.x = *(tinkerForce + 3*ii);
         force.y = *(tinkerForce + 3*ii + 1);
         force.z = *(tinkerForce + 3*ii + 2);
         OpenMM_Vec3Array_set (arrayToLoad, ii, force);
      }
   } else {
      double factor = (double) add;
      for (ii = 0; ii < atoms__.n; ii++) {
         OpenMM_Vec3 force;
         const OpenMM_Vec3* currentForce = OpenMM_Vec3Array_get
                                               (arrayToLoad, ii);
         force.x = currentForce->x + factor * (*(tinkerForce + 3*ii));
         force.y = currentForce->y + factor * (*(tinkerForce + 3*ii + 1));
         force.z = currentForce->z + factor * (*(tinkerForce + 3*ii + 2));
         OpenMM_Vec3Array_set (arrayToLoad, ii, force);
      }
   }
}

static int usingImplicitSolvent (void) {

   int implicitSolventActive = -1;
   char solvatationType[16];
   char bornType[16];

   setNullTerminator (solute__.solvtyp, 8, solvatationType);
   setNullTerminator (solute__.borntyp, 8, bornType);

   // return <0 if parameter/option combination is unsupported
   //         0 if explicit solvent (Ewald is in use)
   //         1 if GK implicit solvent via OBC
   //         2 if GK implicit solvent via Grycuk

   if (strncasecmp (solvatationType, "GK", 2) == 0 &&
            ((strncasecmp (bornType, "OBC", 3) == 0) ||
             (strncasecmp (bornType, "GRYCUK", 6) == 0))) {
      if (limits__.use_ewald) {
         implicitSolventActive = -2;
      } else {
         if ((strncasecmp( bornType, "OBC", 3) == 0)) {
             implicitSolventActive = 1;
         } else {
             implicitSolventActive = 2;
         }
      }
   } else if (limits__.use_ewald) {
      implicitSolventActive = 0;
   }
   return implicitSolventActive;
}

/*
 *    ############################################################
 *       Map Data Structures and Check for Supported Potentials
 *    ############################################################
 *
 *    (1) Map TINKER data structures with arrays to C data structs
 *    (2) Check for TINKER parameter settings consistent with OpenMM
 *        values and potential types
 *    (3) If invalid match or setting found, exit
 */

extern "C" {
void openmm_validate_ (int* testingActive) {

   char buffer[128];
   int ii;
   int implicitSolventActive;
   int invalidAxisType;
   int invalid = 0;
   int testing = *testingActive;
   char const* ixnString = "Interaction Not Currently Supported by OpenMM.\n";
   FILE* log = stderr;

   if (testing && log) {
      (void) fprintf (log, "Testing Mode is Active (=%d)\n", testing);
   }

   // check that requested bond potential type is supported

   if (strncasecmp (bndpot__.bndtyp, "HARMONIC", 8) != 0) {
      invalid++;
      if (log) {
         (void) fprintf (log, "Only HARMONIC bond potential supported: %s not available\n",
                         bndpot__.bndtyp);
      }
   }

   // check that requested angle potential type is supported

   for (ii = 0; ii < angbnd__.nangle; ii++) {
      if (strncasecmp (angpot__.angtyp + 8*ii, "HARMONIC", 8) != 0 &&
         strncasecmp (angpot__.angtyp + 8*ii, "IN-PLANE", 8) != 0 ) {
         invalid++;
         if (log) {
            setNullTerminator (angpot__.angtyp + 8*ii, 8, buffer);
            (void) fprintf (log, "Only HARMONIC or IN-PLANE angle potential supported: %s not available\n",
                            buffer);
         }
      }
   }

   // check that scaling parameters agree with hardcoded values

   if (fabs( 0.0 - vdwpot__.v2scale)  > 0.0 ||
      fabs( 0.0 - vdwpot__.v3scale)  > 0.0 ||
      fabs( 1.0 - vdwpot__.v4scale)  > 0.0 ||
      fabs( 1.0 - vdwpot__.v5scale)  > 0.0 ) {
      invalid++;
      if (log) {
          (void) fprintf (log, "Vdw Scale Parameters [%8.1f %8.1f %8.1f %8.1f] are not [0.0, 0.0, 1.0, 1.0]\n",
                          vdwpot__.v2scale, vdwpot__.v3scale,
                          vdwpot__.v4scale, vdwpot__.v5scale);
      }
   }

   if (fabs(0.0 - polpot__.p2scale) > 0.0 ||
      fabs(0.0 - polpot__.p3scale) > 0.0 ||
      fabs(1.0 - polpot__.p4scale) > 0.0 ||
      fabs(0.5 - polpot__.p41scale) > 0.0 ||
      fabs(1.0 - polpot__.p5scale) > 0.0) {
      invalid++;
      if (log) {
         (void) fprintf (log, "P-Scale Parameters [%8.1f %8.1f %8.1f %8.1f %8.1f] are not [0.0, 0.0, 1.0, 0.5, 1.0]\n",
                         polpot__.p2scale, polpot__.p3scale,
                         polpot__.p4scale, polpot__.p41scale,
                         polpot__.p5scale);
      }
   }

   if (fabs( 0.0 - polpot__.d1scale) > 0.0 ||
      fabs( 1.0 - polpot__.d2scale) > 0.0 ||
      fabs( 1.0 - polpot__.d3scale) > 0.0 ||
      fabs( 1.0 - polpot__.d4scale) > 0.0) {
      invalid++;
      if (log) {
         (void) fprintf (log, "D-Scale Parameters [%8.1f %8.1f %8.1f %8.1f] are not [0.0, 1.0, 1.0, 1.0]\n",
                         polpot__.d1scale, polpot__.d2scale,
                         polpot__.d3scale, polpot__.d4scale);
      }
   }

   if (fabs( 1.0 - polpot__.u1scale) > 0.0 ||
      fabs( 1.0 - polpot__.u2scale) > 0.0 ||
      fabs( 1.0 - polpot__.u3scale) > 0.0 ||
      fabs( 1.0 - polpot__.u4scale) > 0.0) {
      invalid++;
      if (log) {
         (void) fprintf (log, "U-Scale Parameters [%8.1f %8.1f %8.1f %8.1f] are not [1.0, 1.0, 1.0, 1.0]\n",
                         polpot__.u1scale, polpot__.u2scale,
                         polpot__.u3scale, polpot__.u4scale);
      }
   }

   if (fabs( 0.0 - mplpot__.m2scale) > 0.0 ||
      fabs( 0.0 - mplpot__.m3scale) > 0.0 ||
      fabs( 0.4 - mplpot__.m4scale) > 0.0 ||
      fabs( 0.8 - mplpot__.m5scale) > 0.0) {
      invalid++;
      if (log) {
         (void) fprintf (log, "M-Scale Parameters [%8.1f %8.1f %8.1f %8.1f] are not [0.0, 0.0, 0.4, 0.8]\n",
                          mplpot__.m2scale, mplpot__.m3scale,
                          mplpot__.m4scale, mplpot__.m5scale);
      }
   }

   // check that local coordinate frame axis type is supported

   if (potent__.use_mpole) {
      for (ii = 0; ii < atoms__.n; ii++) {
         invalidAxisType = 0;
         if (strncasecmp (mpole__.polaxe + 8*ii, "Z-then-X", 8) != 0 &&
            strncasecmp (mpole__.polaxe + 8*ii, "Bisector", 8) != 0 &&
            strncasecmp (mpole__.polaxe + 8*ii, "Z-Bisect", 8) != 0 &&
            strncasecmp (mpole__.polaxe + 8*ii, "3-Fold", 6) != 0 &&
            strncasecmp (mpole__.polaxe + 8*ii, "None", 4) != 0 &&
            strncasecmp (mpole__.polaxe + 8*ii, "    ", 4) != 0 &&
            strncasecmp (mpole__.polaxe + 8*ii, "Z-Only", 6) != 0) {
            invalid++;
            if (log) {
               setNullTerminator (mpole__.polaxe + 8*ii, 8, buffer);
               (void) fprintf (log, "Axis Type=%s for Atom %7d Not Supported\n",
                               buffer, ii);
               if (invalid > 20) {
                  (void) fflush (log);
                  exit (-1);
               }
            }
         }
      }
   }

   // check for any particles that are "inactive" (ie, frozen)

   if (*(usage__.nuse) !=  atoms__.n) {
      if (log) {
         (void) fprintf (log, "Inactive Atoms Not Supported: nuse=%d, N=%d\n",
                         *(usage__.nuse), atoms__.n);
      }
      invalid++;
   }

   // check for torsional values beyond the supported 3-fold term

   if (potent__.use_tors) {
      double* torsPtr;
      for (ii = 0; ii < tors__.ntors; ii++) {
         torsPtr = tors__.tors4 + ii*4;
         if (*torsPtr != 0.0) {
            if (log) {
                (void) fprintf (log, "4-Fold Parameter for Torsion %d is NonZero %15.7e\n",
                                ii, *torsPtr);
            }
            invalid++;
         }
         torsPtr = tors__.tors5 + ii*4;
         if (*torsPtr != 0.0) {
            if (log) {
               (void) fprintf (log, "5-Fold Parameter for Torsion %d is NonZero %15.7e\n",
                               ii, *torsPtr);
            }
            invalid++;
         }
         torsPtr = tors__.tors6 + ii*4;
         if (*torsPtr != 0.0) {
            if (log) {
               (void) fprintf (log, "6-Fold Parameter for Torsion %d is NonZero %15.7e\n",
                               ii, *torsPtr);
            }
            invalid++;
         }
      }
   }

   // check for use of potential terms not implemented in OpenMM

   if (potent__.use_angang) {
      invalid++;
      if (log) {
         (void) fprintf (log, "Angle-Angle %s", ixnString);
      }
   }
   if (potent__.use_opdist) {
      invalid++;
      if (log) {
         (void) fprintf (log, "Out-of-Plane Distance %s", ixnString);
      }
   }
   if (potent__.use_improp) {
      invalid++;
      if (log) {
         (void) fprintf (log, "Improper Dihedral %s", ixnString);
      }
   }
   if (potent__.use_imptor) {
      invalid++;
      if (log) {
         (void) fprintf (log, "Improper Torsion %s", ixnString);
      }
   }
   if (potent__.use_charge) {
      invalid++;
      if (log) {
         (void) fprintf (log, "Charge-Charge %s", ixnString);
      }
   }
   if (potent__.use_chgdpl) {
      invalid++;
      if (log) {
         (void) fprintf (log, "Atomwise Charge-Dipole %s", ixnString);
      }
   }
   if (potent__.use_dipole) {
      invalid++;
      if (log) {
         (void) fprintf (log, "Dipole-Dipole %s", ixnString);
      }
   }
   if (potent__.use_rxnfld) {
      invalid++;
      if (log) {
         (void) fprintf (log, "Reaction Field %s", ixnString);
      }
   }
   if (potent__.use_metal) {
      invalid++;
      if (log) {
         (void) fprintf (log, "Ligand Field %s", ixnString);
      }
   }
   if (potent__.use_geom) {
      invalid++;
      if (log) {
         (void) fprintf (log, "Geometric Restraint %s", ixnString);
      }
   }
/*
   if (potent__.use_extra) {
      invalid++;
      if (log) {
         (void) fprintf (log, "Extra Potentials Not Supported\n");
      }
   }
*/
   if (potent__.use_orbit) {
      invalid++;
      if (log) {
         (void) fprintf (log, "SCF-Pi-MO Calculation Not Supported\n");
      }
   }
   implicitSolventActive = usingImplicitSolvent ();
   if (implicitSolventActive == -1 && !testing) {
      if (log) {
         (void) fprintf (log, "Either (Solvate= 'GK' && Born-Radius='GRYCUK') or the Ewald Parameter Must be Set (but not both)\n");
      }
      invalid++;
   }
   if (implicitSolventActive == -2 && !testing) {
      if (log) {
         (void) fprintf (log, "Unsupported Force Combination: Both solvate= 'GK' && born-radius='GRYCUK' and Ewald parameter cannot be set\n");
      }
      invalid++;
   }
   setNullTerminator (mdstuf__.integrate, 11, buffer);
   if (strncasecmp (buffer, "VERLET", 6) != 0 &&
      strncasecmp (buffer, "STOCHASTIC", 6) != 0 && !testing ) {
      invalid++;
      if (log) {
         (void) fprintf (log, "Integrator %s is Not Supported\n", buffer);
      }
   }
   if (log) {
      (void) fflush (log);
   }

   if (invalid && !testing) {
      exit (-1);
   }
}
}

/*
 *    ############################################################
 *           Compare TINKER and OpenMM Energies and Forces
 *    ############################################################
 */

extern "C" {
int openmm_test_ (void) {

   OpenMM_Vec3Array* initialPosInNm;
   OpenMM_Platform* platform;
   OpenMM_System* system;
   OpenMM_Context* context;
   OpenMM_Integrator* integrator;
   OpenMM_State* state;

   const OpenMM_Vec3Array* posArrayInNm;
   const OpenMM_Vec3Array* openMMForces;

   int infoMask;
   int ii, jj;
   int countActiveForces;
   char const* testName;
   double conversion, delta, dot;
   double tinkerNorm, openMMNorm;
   double openMMPotentialEnergy;
   OpenMM_Vec3Array* tinkerForce;
   double tinkerEnergy;
   double maxEDelta, maxFDelta;
   double minDot, avgDot;
   int maxFDeltaIndex, minDotIndex;
   int implicitSolventActive;

   int removeConstrainedCovalentIxns = 0;
   FILE* log = stderr;

   if (log) {
      (void) fprintf (log, "\n Testing TINKER vs OpenMM Energy & Force :\n");
   }

   implicitSolventActive = usingImplicitSolvent ();

   tinkerForce = OpenMM_Vec3Array_create (atoms__.n);

   // Create a System and Force objects within the System. Retain a reference
   // to each force object so we can fill in the forces. Note: the OpenMM
   // System takes ownership of the force objects; don't delete them yourself.

   testName = NULL;
   system = OpenMM_System_create ();
   setupSystemParticles (system, log);

   countActiveForces = 0;
   if (potent__.use_bond)  countActiveForces++;
   if (potent__.use_angle)  countActiveForces++;
   if (potent__.use_strbnd)  countActiveForces++;
   if (potent__.use_urey)  countActiveForces++;
   if (potent__.use_angang)  countActiveForces++;
   if (potent__.use_opbend)  countActiveForces++;
   if (potent__.use_opdist)  countActiveForces++;
   if (potent__.use_improp)  countActiveForces++;
   if (potent__.use_imptor)  countActiveForces++;
   if (potent__.use_tors)  countActiveForces++;
   if (potent__.use_pitors)  countActiveForces++;
   if (potent__.use_strtor)  countActiveForces++;
   if (potent__.use_angtor)  countActiveForces++;
   if (potent__.use_tortor)  countActiveForces++;
   if (potent__.use_vdw)  countActiveForces++;
   if (potent__.use_ct)  countActiveForces++;
   if (potent__.use_charge)  countActiveForces++;
   if (potent__.use_chgdpl)  countActiveForces++;
   if (potent__.use_dipole)  countActiveForces++;
   if (potent__.use_mpole)  countActiveForces++;
   if (potent__.use_polar)  countActiveForces++;
   if (potent__.use_rxnfld)  countActiveForces++;
   if (potent__.use_solv)  countActiveForces++;
   if (potent__.use_metal)  countActiveForces++;
   if (potent__.use_geom)  countActiveForces++;

   if (log) {
      (void) fprintf (log, "\n Potential Terms Used in TINKER :\n" );
      (void) fprintf (log, "\n    Bond=    %d", abs(potent__.use_bond));
      (void) fprintf (log, "    Angle=   %d", abs(potent__.use_angle));
      (void) fprintf (log, "    StrBnd=  %d", abs(potent__.use_strbnd));
      (void) fprintf (log, "    Urey=    %d", abs(potent__.use_urey));
      (void) fprintf (log, "    AngAng=  %d", abs(potent__.use_angang));
      (void) fprintf (log, "\n    OPBend=  %d", abs(potent__.use_opbend));
      (void) fprintf (log, "    OPDist=  %d", abs(potent__.use_opdist));
      (void) fprintf (log, "    ImProp=  %d", abs(potent__.use_improp));
      (void) fprintf (log, "    ImpTor=  %d", abs(potent__.use_imptor));
      (void) fprintf (log, "    Tors=    %d", abs(potent__.use_tors));
      (void) fprintf (log, "\n    PiTors=  %d", abs(potent__.use_pitors));
      (void) fprintf (log, "    StrTor=  %d", abs(potent__.use_strtor));
      (void) fprintf (log, "    AngTor=  %d", abs(potent__.use_angtor));
      (void) fprintf (log, "    TorTor=  %d", abs(potent__.use_tortor));
      (void) fprintf (log, "    Vdw=     %d", abs(potent__.use_vdw));
      (void) fprintf (log, "    CT=     %d", abs(potent__.use_ct));
      (void) fprintf (log, "\n    Charge=  %d", abs(potent__.use_charge));
      (void) fprintf (log, "    ChgDpl=  %d", abs(potent__.use_chgdpl));
      (void) fprintf (log, "    Dipole=  %d", abs(potent__.use_dipole));
      (void) fprintf (log, "    MPole=   %d", abs(potent__.use_mpole));
      (void) fprintf (log, "    Polar=   %d", abs(potent__.use_polar));
      (void) fprintf (log, "\n    RxnFld=  %d", abs(potent__.use_rxnfld));
      (void) fprintf (log, "    Solv=    %d", abs(potent__.use_solv));
      (void) fprintf (log, "    LigFld=  %d", abs(potent__.use_metal));
      (void) fprintf (log, "    Restrn=  %d", abs(potent__.use_geom));
      (void) fprintf (log, "    Extra=   %d\n", abs(potent__.use_extra));
   }

   if (countActiveForces > 1) {

      if (potent__.use_bond) {
         setupAmoebaBondForce (system, removeConstrainedCovalentIxns, log);
      }
      if (potent__.use_angle) {
         setupAmoebaAngleForce (system, removeConstrainedCovalentIxns, log);
         setupAmoebaInPlaneAngleForce (system, removeConstrainedCovalentIxns,
                                       log);
      }
      if (potent__.use_strbnd) {
         setupAmoebaStretchBendForce (system, removeConstrainedCovalentIxns,
                                      log);
      }
      if (potent__.use_urey) {
         setupAmoebaUreyBradleyForce (system, removeConstrainedCovalentIxns,
                                      log);
      }
      if (potent__.use_opbend) {
         setupAmoebaOutOfPlaneBendForce (system, log);
      }
      if (potent__.use_imptor) {
         setupAmoebaImproperTorsionForce (system, log);
      }
      if (potent__.use_tors) {
         setupAmoebaTorsionForce (system, log);
      }
      if (potent__.use_pitors) {
         setupAmoebaPiTorsionForce (system, log);
      }
      if (potent__.use_strtor) {
         setupAmoebaStretchTorsionForce (system, log);
      }
      if (potent__.use_angtor) {
         setupAmoebaAngleTorsionForce (system, log);
      }
      if (potent__.use_tortor) {
         setupAmoebaTorsionTorsionForce (system, log);
      }
      if (potent__.use_vdw) {
         setupAmoebaVdwForce (system, log);
      }
      if (potent__.use_ct) {
         setupAmoebaCTForce (system, log);
      }
      if (potent__.use_charge) {
         setupAmoebaChargeForce (system, log);
      }
      if (potent__.use_mpole) {
         setupAmoebaMultipoleForce (system, log);
      }
      if (potent__.use_solv) {
         setupAmoebaWcaDispersionForce (system, log);
         setupAmoebaGeneralizedKirkwoodForce (system, 1, log);
      }
      if (potent__.use_geom) {
         setupTorsionRestraints (system, log);
         setupDistanceRestraints (system, log);
         setupPositionalRestraints (system, log);
         setupAngleRestraints (system, log);
         setupCentroidRestraints (system, log);
      }

      loadTinkerForce (deriv__.desum, 0, tinkerForce);
      tinkerEnergy = *energi__.esum;
      testName = "PotentialsTest";

      if (log) {
         (void) fprintf (log,
                    "\n Potential Energy Components from TINKER :\n\n"
                    "    EB=  %15.7e   EA=  %15.7e   EBA= %15.7e\n"
                    "    EUB= %15.7e   EAA= %15.7e   EOPB=%15.7e\n"
                    "    EOPD=%15.7e   EID= %15.7e   EIT= %15.7e\n"
                    "    ET=  %15.7e   EPT= %15.7e   EBT= %15.7e\n"
                    "    EAT= %15.7e   ETT= %15.7e   EV=  %15.7e\n"
                    "    EC=  %15.7e   ECD= %15.7e   ED=  %15.7e\n"
                    "    EM=  %15.7e   EP=  %15.7e   ER=  %15.7e\n"
                    "    ES=  %15.7e   ELF= %15.7e   EG=  %15.7e\n"
                    "    EX=  %15.7e   ECT= %15.7e\n",
                    *energi__.eb,   *energi__.ea,  *energi__.eba,
                    *energi__.eub,  *energi__.eaa, *energi__.eopb,
                    *energi__.eopd, *energi__.eid, *energi__.eit,
                    *energi__.et,   *energi__.ept, *energi__.ebt,
                    *energi__.eat,  *energi__.ett, *energi__.ev,
                    *energi__.ec,   *energi__.ecd, *energi__.ed,
                    *energi__.em,   *energi__.ep,  *energi__.er,
                    *energi__.es,   *energi__.elf, *energi__.eg,
                    *energi__.ex,   *energi__.ect);
         (void) fflush (log);
      }

   } else if (potent__.use_bond) {

      setupAmoebaBondForce (system, removeConstrainedCovalentIxns, log);
      loadTinkerForce (deriv__.deb, 0, tinkerForce);
      tinkerEnergy = *energi__.eb;
      testName = "AmoebaHarmonicBondTest";

   } else if (potent__.use_angle) {

      // note TINKER angle = OpenMM (Angle + InPlaneAngle)

      setupAmoebaAngleForce (system, removeConstrainedCovalentIxns, log);
      setupAmoebaInPlaneAngleForce (system, removeConstrainedCovalentIxns,
                                    log);
      loadTinkerForce (deriv__.dea, 0, tinkerForce);
      tinkerEnergy = *energi__.ea;
      testName = "AmoebaHarmonicAngleTest";

   } else if (potent__.use_strbnd) {

      setupAmoebaStretchBendForce (system, removeConstrainedCovalentIxns,
                                   log);
      loadTinkerForce (deriv__.deba, 0, tinkerForce);
      tinkerEnergy = *energi__.eba;
      testName = "AmoebaStretchBendTest";

   } else if (potent__.use_urey) {

      setupAmoebaUreyBradleyForce (system, removeConstrainedCovalentIxns,
                                   log);
      loadTinkerForce (deriv__.deub, 0, tinkerForce);
      tinkerEnergy = *energi__.eub;
      testName = "AmoebaUreyBradleyForceTest";

   } else if (potent__.use_opbend) {

      setupAmoebaOutOfPlaneBendForce (system, log);
      loadTinkerForce (deriv__.deopb, 0, tinkerForce);
      tinkerEnergy = *energi__.eopb;
      testName = "AmoebaOutOfPlaneBendTest";

   } else if (potent__.use_imptor) {

      setupAmoebaImproperTorsionForce (system, log);
      loadTinkerForce (deriv__.deit, 0, tinkerForce);
      tinkerEnergy = *energi__.eit;
      testName = "AmoebaImproperTorsionForce";

   } else if (potent__.use_tors) {

      setupAmoebaTorsionForce (system, log);
      loadTinkerForce (deriv__.det, 0, tinkerForce);
      tinkerEnergy = *energi__.et;
      testName = "AmoebaTorsionTest";

   } else if (potent__.use_pitors) {

      setupAmoebaPiTorsionForce (system, log);
      loadTinkerForce (deriv__.dept, 0, tinkerForce);
      tinkerEnergy = *energi__.ept;
      testName = "AmoebaPiTorsionTest";

   } else if (potent__.use_strtor) {

      setupAmoebaStretchTorsionForce (system, log);
      loadTinkerForce (deriv__.debt, 0, tinkerForce);
      tinkerEnergy = *energi__.ebt;
      testName = "AmoebaStretchTorsionTest";

   } else if (potent__.use_angtor) {

      setupAmoebaAngleTorsionForce (system, log);
      loadTinkerForce (deriv__.deat, 0, tinkerForce);
      tinkerEnergy = *energi__.eat;
      testName = "AmoebaAngleTorsionTest";

   } else if (potent__.use_tortor) {

      setupAmoebaTorsionTorsionForce (system, log);
      loadTinkerForce (deriv__.dett, 0, tinkerForce);
      tinkerEnergy = *energi__.ett;
      testName = "AmoebaTorsionTorsionTest";

   } else if (potent__.use_vdw) {

      setupAmoebaVdwForce (system, log);
      loadTinkerForce (deriv__.dev, 0, tinkerForce);
      tinkerEnergy = *energi__.ev;
      if (limits__.use_vlist) {
         testName = "AmoebaVdwTest";
      } else {
         testName = "AmoebaVdwNoCutoffTest";
      }
      if (log) {
         (void) fprintf (log, "\n Use Vdw Neighbor List=%d Vdw Cutoff=%12.3f\n",
                         limits__.use_vlist, limits__.vdwcut);
      }

   } else if (potent__.use_ct) {

      setupAmoebaCTForce (system, log);
      loadTinkerForce (deriv__.dect, 0, tinkerForce);
      tinkerEnergy = *energi__.ect;
      if (limits__.use_ctlist) {
         testName = "AmoebaCTTest";
      } else {
         testName = "AmoebaCTNoCutoffTest";
      }
      if (log) {
         (void) fprintf (log, "\n Use CT Neighbor List=%d CT Cutoff=%12.3f\n",
                         limits__.use_ctlist, limits__.ctcut);
      }

   } else if (potent__.use_geom) {

      setupTorsionRestraints (system, log);
      setupDistanceRestraints (system, log);
      setupPositionalRestraints (system, log);
      setupAngleRestraints (system, log);
      setupCentroidRestraints (system, log);
      loadTinkerForce (deriv__.deg, 0, tinkerForce);
      tinkerEnergy = *energi__.eg;
      testName = "AmoebaRestraintTest";

   } else if (potent__.use_charge) {

      setupAmoebaChargeForce (system, log);
      loadTinkerForce (deriv__.dec, 0, tinkerForce);
      tinkerEnergy = *energi__.ec;
      testName = "AmoebaChargeTest";

   } else if (potent__.use_mpole && !potent__.use_solv) {

      if (log) {
         (void) fprintf (log, "ImplicitSolventActive=%d\n",
                         implicitSolventActive);
      }
      if (implicitSolventActive == 0) {
         setupAmoebaBondForce (system, removeConstrainedCovalentIxns, log);
         loadTinkerForce (deriv__.deb, 0, tinkerForce);
         tinkerEnergy = *energi__.eb;
      } else {
         zeroVec3Force (tinkerForce);
         tinkerEnergy = 0.0;
      }
      setupAmoebaMultipoleForce (system, log);
      loadTinkerForce (deriv__.dem, 1, tinkerForce);
      loadTinkerForce (deriv__.dep, 1, tinkerForce);
      tinkerEnergy += *energi__.em + *energi__.ep;

      if (strncasecmp (polpot__.poltyp, "DIRECT", 6) == 0) {
         testName = "AmoebaMultipoleDirectTest";
      } else {
         testName = "AmoebaMultipoleMutualTest";
      }

   } else if (potent__.use_solv && implicitSolventActive > 0 &&
              !potent__.use_mpole) {

      // to get TINKER WCA, zero deriv__.des, then call ewca1

      zeroTinkerForce (deriv__.des);
      ewca1_ (&tinkerEnergy);
      loadTinkerForce (deriv__.des, 0, tinkerForce);
      setupAmoebaWcaDispersionForce (system, log);
      testName = "AmoebaWcaDispersionTest";

   } else if (implicitSolventActive > 0 && potent__.use_mpole) {

      // generalized Kirkwood; OpenMM should have cavity term turned off

      double ecav, edisp;
      if (log) {
         char buffer[128];
         setNullTerminator (solute__.borntyp, 8, buffer);
         (void) fprintf (log, "Born radius type=%s ", buffer);
         setNullTerminator (solute__.solvtyp, 8, buffer);
         (void) fprintf (log, "Solvation type=%s\n", buffer);
      }
      setupAmoebaMultipoleForce (system, log);
      setupAmoebaGeneralizedKirkwoodForce (system, 1, log);
      zeroTinkerForce (deriv__.des);
      setTinker1DArray (atoms__.n, solute__.drb, 0.0);
      setTinker1DArray (atoms__.n, solute__.drbp, 0.0);

      *energi__.es = 0.0;
      born_ ();
      empole1_ ();
      egk1_ ();

      loadTinkerForce (deriv__.des, 0, tinkerForce);
      loadTinkerForce (deriv__.dem, 1, tinkerForce);
      loadTinkerForce (deriv__.dep, 1, tinkerForce);
      tinkerEnergy = *energi__.es + *energi__.em + *energi__.ep;

      enp1_ (&ecav,&edisp);

      if (log) {
         (void) fprintf (log, "Energies total=%15.7e Gk=%15.7e (cavity=%15.7e dispersion=%15.7e) Multipole=%15.7e\n",
                         tinkerEnergy, *energi__.es, ecav, edisp,
                         *energi__.em + *energi__.ep );

         //for (ii = 0; ii < atoms__.n; ii++) {
         //    fprintf (stderr, "Fs %5d [%15.7e %15.7e %15.7e] [%15.7e %15.7e%15.7e] [%15.7e %15.7e %15.7e]\n",
         //             ii, *(deriv__.des + 3*ii), *(deriv__.des + 3*ii+1),
         //             *(deriv__.des + 3*ii+2), *(deriv__.dem + 3*ii),
         //             *(deriv__.dem + 3*ii+1), *(deriv__.dem + 3*ii+2),
         //             *(deriv__.dep + 3*ii), *(deriv__.dep + 3*ii+1),
         //             *(deriv__.dep + 3*ii+2));
         //}

      }
      if (strncasecmp (polpot__.poltyp, "DIRECT", 6) == 0) {
         testName = "AmoebaKirkwoodDirectTest";
      } else {
         testName = "AmoebaKirkwoodMutualTest";
      }
   }

   if (log) {
      if (testName) {
         (void) fprintf (log, "\n Test Option :  %s\n", testName);
         (void) fflush (NULL);
      } else {
         (void) fprintf (log, "\n Test Option not Recognized; Exiting\n");
         (void) fflush (log);
         exit (-1);
      }
   }

   initialPosInNm = OpenMM_Vec3Array_create (0);
   setupPositions (initialPosInNm, log);

   integrator = (OpenMM_Integrator*) OpenMM_VerletIntegrator_create (0.001);

   // choose either the reference or the CUDA platform

   // platform = getReferencePlatform (log);
   platform = getCUDAPlatform (log);
   if (platform == NULL) {
      exit (-1);
   }

   context = OpenMM_Context_create_2 (system, integrator, platform);

   OpenMM_Context_setPositions (context, initialPosInNm);

   infoMask = OpenMM_State_Positions;
   infoMask += OpenMM_State_Forces;
   infoMask += OpenMM_State_Energy;

   state = OpenMM_Context_getState (context, infoMask, 0);

   openMMForces = OpenMM_State_getForces (state);
   openMMPotentialEnergy = OpenMM_State_getPotentialEnergy(state)
                              / OpenMM_KJPerKcal;

   conversion = -OpenMM_NmPerAngstrom / OpenMM_KJPerKcal;

   setDefaultPeriodicBoxVectors (system, log);
   if (log) {
      printDefaultPeriodicBoxVectors (system, log);
      (void) fflush (log);
   }

   // find differences in TINKER and OpenMM energies and forces

   maxFDelta = 0.0;
   maxFDeltaIndex = -1;
   minDot = 1.0;
   avgDot = 0.0;
   minDotIndex = -1;

   int maxPrint;
   maxPrint = 5;

   if (log) {

      (void) fprintf (log, "\n TINKER vs OpenMM Energy Values :\n");
      (void) fprintf (log, "\n TINKER Potential Energy   %18.4f",
                      tinkerEnergy);
      (void) fprintf (log, "\n OpenMM Potential Energy   %18.4f\n",
                      openMMPotentialEnergy);
      maxEDelta = fabs(tinkerEnergy - openMMPotentialEnergy);

      maxFDelta = 0.0;
      for (ii = 0; ii < atoms__.n; ii++) {
         double relxNrm;
         double dot;
         OpenMM_Vec3 force;
         const OpenMM_Vec3* tinkerF;
         force.x = conversion*(*OpenMM_Vec3Array_get(openMMForces, ii)).x;
         force.y = conversion*(*OpenMM_Vec3Array_get(openMMForces, ii)).y;
         force.z = conversion*(*OpenMM_Vec3Array_get(openMMForces, ii)).z;
         tinkerF = OpenMM_Vec3Array_get(tinkerForce, ii);

         tinkerNorm = sqrt(tinkerF->x*tinkerF->x + tinkerF->y*tinkerF->y
                               + tinkerF->z*tinkerF->z);
         openMMNorm = sqrt(force.x*force.x + force.y*force.y
                               + force.z*force.z );

         delta = sqrt((tinkerF->x-force.x)*(tinkerF->x-force.x)
                    + (tinkerF->y-force.y)*(tinkerF->y-force.y)
                    + (tinkerF->z-force.z)*(tinkerF->z-force.z));
         dot = ((tinkerNorm > 0.0) && (openMMNorm > 0.0)) ?
                (tinkerF->x*force.x + tinkerF->y*force.y
                     + tinkerF->z*force.z) / (tinkerNorm*openMMNorm) : 0.0;

         if (delta > maxFDelta) {
            maxFDelta = delta;
            maxFDeltaIndex = ii + 1;
         }
         if (dot < minDot) {
            minDot = dot;
            minDotIndex = ii + 1;
         }
         avgDot += dot;

         if (ii == 0) {
            (void) fprintf (log, "\n TINKER vs OpenMM Force Values :\n\n");
         }
         if (ii < maxPrint || atoms__.n - ii - 1 < maxPrint) {
            (void) fprintf (stderr, "%6d   TINKER %17.8e %17.8e %17.8e",
                            ii+1, tinkerF->x, tinkerF->y, tinkerF->z);
            (void) fprintf (stderr, "\n         OpenMM %17.8e %17.8e %17.8e\n",
                            force.x, force.y, force.z);
         }
      }

      if (atoms__.n) {
         avgDot /= (double)(atoms__.n);
      }

      (void) fprintf (log, "\n Summary of TINKER vs OpenMM Comparison :\n");
      (void) fprintf (log, "\n EnergyDelta                    %19.8e\n",
                      maxEDelta);
      (void) fprintf (log, " MaxForceDelta at   %11d %19.8e\n",
                      maxFDeltaIndex, maxFDelta);
      (void) fprintf (log, " MinDotProduct at   %11d %19.8e\n",
                      minDotIndex, minDot);
      (void) fprintf (log, " AvgDotProduct                  %19.8e\n",
                      avgDot);
   }

   OpenMM_Vec3Array_destroy (tinkerForce);
   OpenMM_State_destroy (state);
   OpenMM_Context_destroy (context);
   OpenMM_Integrator_destroy (integrator);
   OpenMM_System_destroy (system);
}
}
