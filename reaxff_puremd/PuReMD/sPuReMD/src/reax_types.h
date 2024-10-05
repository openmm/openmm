/*----------------------------------------------------------------------
  SerialReax - Reax Force Field Simulator

  Copyright (2010) Purdue University
  Hasan Metin Aktulga, haktulga@cs.purdue.edu
  Joseph Fogarty, jcfogart@mail.usf.edu
  Sagar Pandit, pandit@usf.edu
  Ananth Y Grama, ayg@cs.purdue.edu

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the GNU General Public License for more details:
  <http://www.gnu.org/licenses/>.
  ----------------------------------------------------------------------*/

#ifndef __REAX_TYPES_H_
#define __REAX_TYPES_H_

#if (defined(HAVE_CONFIG_H) && !defined(__CONFIG_H_))
  #define __CONFIG_H_
  #include "../../common/include/config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(_OPENMP)
  #include <omp.h>
#endif

/* enables debugging code */
//#define DEBUG_FOCUS
/* enables test energy code */
//#define TEST_ENERGY
/* enables test forces code */
//#define TEST_FORCES
/* constants defined in reference Fortran ReaxFF code (useful for comparisons) */
#define USE_REF_FORTRAN_REAXFF_CONSTANTS
/* constants defined in reference Fortran eReaxFF code (useful for comparisons) */
//#define USE_REF_FORTRAN_EREAXFF_CONSTANTS
/* constants defined in LAMMPS ReaxFF code (useful for comparisons) */
//#define USE_LAMMPS_REAXFF_CONSTANTS
/* enables reordering atoms after neighbor list generation for improved cache performance */
//#define REORDER_ATOMS
/* enables support for small simulation boxes (i.e. a simulation box with any
 * dimension less than twice the Verlet list cutoff distance, vlist_cut),
 * which means multiple periodic images must be searched for interactions */
#define SMALL_BOX_SUPPORT

/* disable assertions if NOT compiling with debug support --
 * the definition (or lack thereof) controls how the assert macro is defined */
#if !defined(DEBUG) && !defined(DEBUG_FOCUS)
  #define NDEBUG
#endif
#include <assert.h>

#define SUCCESS (1)
#define FAILURE (0)
#define TRUE (1)
#define FALSE (0)

#if defined(USE_REF_FORTRAN_REAXFF_CONSTANTS)
  /* transcendental constant pi */
  #define PI (3.14159265)
  /* unit conversion, ??? to kcal / mol */
  #define C_ELE (332.0638)
  /* Boltzmann constant, in J / mol / K */
//  #define K_B (0.831687)
  #define K_B (0.8314510)
  /* unit conversion for atomic force, kcal to J */
  #define F_CONV (4.184e2)
  /* energy conversion constant from electron volts to kilo-calories per mole */
  #define KCALpMOL_to_EV (23.02)
  /* electric dipole moment conversion constant from elementary charge * angstrom to debye */
  #define ECxA_to_DEBYE (4.80320679913)
#elif defined(USE_REF_FORTRAN_EREAXFF_CONSTANTS)
  //TODO
  /* energy conversion constant from electron volts to kilo-calories per mole */
  #define KCALpMOL_to_EV (23.02)
  /* electric dipole moment conversion constant from elementary charge * angstrom to debye */
  #define ECxA_to_DEBYE (4.80320679913)
#elif defined(USE_LAMMPS_REAXFF_CONSTANTS)
  //TODO
  /* energy conversion constant from electron volts to kilo-calories per mole */
  #define KCALpMOL_to_EV (23.060549)
#endif

/* transcendental constant pi */
#if !defined(PI)
  #if defined(M_PI)
    /* GNU C library (libc), defined in math.h */
    #define PI (M_PI)
  #else
    #define PI (3.14159265)
  #endif
#endif
/* Coulomb energy conversion */
#if !defined(C_ELE)
  #define C_ELE (332.06371)
#endif
/* Boltzmann constant */
#if !defined(K_B)
  /* in ??? */
//  #define K_B (503.398008)
  /* in J / mol / K */
//  #define K_B (0.831687)
  #define K_B (0.8314510)
#endif
/* unit conversion for atomic force */
#if !defined(F_CONV)
  /* ??? to AMU * A / ps^2 */
  #define F_CONV (1.0e6 / 48.88821291 / 48.88821291)
#endif
/* unit conversion for atomic energy */
#if !defined(E_CONV)
  /* AMU * Angstroms^2 / ps^2 --> kcal / mol */
//  #define E_CONV (0.002391)
  #define E_CONV (1.0 / 418.40)
#endif
/* energy conversion constant from electron volts to kilo-calories per mole */
#if !defined(EV_to_KCALpMOL)
  #define EV_to_KCALpMOL (14.40)
#endif
/* energy conversion constant from electron volts to kilo-calories per mole */
#if !defined(KCALpMOL_to_EV)
  #define KCALpMOL_to_EV (23.0408)
#endif
/* electric dipole moment conversion constant from elementary charge * angstrom to debye */
#if !defined(ECxA_to_DEBYE)
  #define ECxA_to_DEBYE (4.803204)
#endif
/* energy conversion constant from (gram) calories to Joules (SI) */
#if !defined(CAL_to_JOULES)
  #define CAL_to_JOULES (4.1840)
#endif
/* energy conversion constant from Joules (SI) to (gram) calories */
#if !defined(JOULES_to_CAL)
  #define JOULES_to_CAL (1.0 / 4.1840)
#endif
/* mass conversion constant from unified atomic mass units (daltons) to grams */
#if !defined(AMU_to_GRAM)
  #define AMU_to_GRAM (1.6605e-24)
#endif
/* distance conversion constant from angstroms to centimeters */
#if !defined(ANG_to_CM)
  #define ANG_to_CM (1.0e-8)
#endif
/* Avogradro's constant */
#if !defined(AVOGNR)
  #define AVOGNR (6.0221367e23)
#endif
/* unit conversion for pressure:
 * (1 s / 10^12 ps) * (1 m / 10^10 Angstroms) * (6.0221367^23 atoms / mole) * (0.2390057 J / cal)
 * ps * Angstroms * moles * cals => s * m * atoms * J */
#if !defined(P_CONV)
  #define P_CONV (1.0e-24 * AVOGNR * JOULES_to_CAL)
#endif

/* max. num. of characters for string buffers */
#define MAX_STR (1024)
/* max. num. of characters for a line in files */
#define MAX_LINE (1024)
/* max. num. of tokens per line */
#define MAX_TOKENS (1024)
/* max. num. of characters per token */
#define MAX_TOKEN_LEN (1024)

/* max. atom ID in geo file */
#define MAX_ATOM_ID (100000)
/* ??? */
#define MAX_RESTRICT (15)
/* max. num. atoms per molecule */
#define MAX_MOLECULE_SIZE (20)
/* max. num. atom types defined in the force field parameter file */
#define MAX_ATOM_TYPES (25)

#define MAX_GRID (50)
#define MAX_3BODY_PARAM (5)
#define MAX_4BODY_PARAM (5)
#define NUM_INTRS (10)

/* max. pressure scaler for simulation box dimenion in NPT ensembles */
#define MAX_dV (1.01)
/* min. pressure scaler for simulation box dimenion in NPT ensembles */
#define MIN_dV (1.0 / MAX_dV)
/* max. temperature scaler for atomic positions and velocities in NPT ensembles */
#define MAX_dT (2.0)
/* min. temperature scaler for atomic positions and velocities in NPT ensembles */
#define MIN_dT (0.0)

/* max. num. of main simulation loop retries;
 * retries occur when memory allocation checks determine more memory is needed */
#define MAX_RETRIES (5)

#define ALMOST_ZERO (1.0e-10)
#define NEG_INF (-1.0e10)
#define NO_BOND (1.0e-3)
#define HB_THRESHOLD (1.0e-2)
#define MAX_BONDS (40)
#define MIN_BONDS (15)
#define MIN_HBONDS (50)
#define SAFE_HBONDS (1.4)
#define MIN_GCELL_POPL (50)
#define SAFE_ZONE (1.2)
#define DANGER_ZONE (0.95)
#define LOOSE_ZONE (0.75)

/* NaN IEEE 754 representation for C99 in math.h
 * Note: function choice must match REAL typedef below */
#if defined(NAN)
  #define IS_NAN_REAL(a) (isnan(a))
#else
  #warn "No support for NaN"
  #define IS_NAN_REAL(a) (0)
#endif
#define LOG (log)
#define EXP (exp)
#define SQRT (sqrt)
#define POW (pow)
#define ACOS (acos)
#define COS (cos)
#define SIN (sin)
#define TAN (tan)
#define CEIL (ceil)
#define FLOOR (floor)
#define ABS (abs)
#define FABS (fabs)
#define FMOD (fmod)
#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define DEG2RAD(a) ((a)/180.0*PI)
#define RAD2DEG(a) ((a)/PI*180.0)
#define MAX(x,y) (((x) > (y)) ? (x) : (y))
#define MIN(x,y) (((x) < (y)) ? (x) : (y))


/* ensemble type */
enum ensemble
{
    /* microcanonical ensemble */
    NVE = 0,
    /* Berendsen NVT ensemble */
    bNVT = 1,
    /* Nose-Hoover NVT ensemble */
    nhNVT = 2,
    /* Parrinello-Rehman-Nose-Hoover semi-isotropic NPT ensemble */
    sNPT = 3,
    /* Parrinello-Rehman-Nose-Hoover isotropic NPT ensemble */
    iNPT = 4,
    /* Parrinello-Rehman-Nose-Hoover anisotropic NPT ensemble */
    aNPT = 5,
    /* total number of ensemble types */
    ens_N = 6,
};

/* interaction list type */
enum interaction_list_offets
{
    FAR_NBRS = 0,
    NEAR_NBRS = 1,
    THREE_BODIES = 2,
    BONDS = 3,
    OLD_BONDS = 4,
    HBONDS = 5,
    DBO = 6,
    DDELTA = 7,
    LIST_N = 8,
};

/* interaction type */
enum interaction_type
{
    TYP_THREE_BODY = 0,
    TYP_BOND = 1,
    TYP_HBOND = 2,
    TYP_DBO = 3,
    TYP_DDELTA = 4,
    TYP_FAR_NEIGHBOR = 5,
    TYP_NEAR_NEIGHBOR = 6,
    TYP_N = 7,
};

/* error codes for simulation termination */
enum errors
{
    FILE_NOT_FOUND = -10,
    UNKNOWN_ATOM_TYPE = -11,
    CANNOT_OPEN_FILE = -12,
    CANNOT_INITIALIZE = -13,
    INSUFFICIENT_MEMORY = -14,
    UNKNOWN_OPTION = -15,
    INVALID_INPUT = -16,
    INVALID_GEO = -17,
    NUMERIC_BREAKDOWN = -18,
    MAX_RETRIES_REACHED = -19,
    RUNTIME_ERROR = -20,
};

enum molecular_analysis_type
{
    NO_ANALYSIS = 0,
    FRAGMENTS = 1,
    REACTIONS = 2,
    NUM_ANALYSIS = 3,
};

/* restart file format */
enum restart_formats
{
    WRITE_ASCII = 0,
    WRITE_BINARY = 1,
    RF_N = 2,
};

/* geometry file format */
enum geo_formats
{
    CUSTOM = 0,
    PDB = 1,
    BGF = 2,
    ASCII_RESTART = 3,
    BINARY_RESTART = 4,
    GF_N = 5,
};

/* method used for computing atomic charges */
enum charge_method
{
    QEQ_CM = 0,
    EE_CM = 1,
    ACKS2_CM = 2,
};

/* iterative linear solver used for computing atomic charges */
enum solver
{
    GMRES_S = 0,
    GMRES_H_S = 1,
    CG_S = 2,
    SDM_S = 3,
    BiCGStab_S = 4,
};

/* initial guess type for the linear solver
 * used for computing atomic charges */
enum init_guess_type
{
    /* prediction using spline extraplotation */
    SPLINE = 0,
    /* prediction using Tensorflow frozen model trained with
     * the long short-term memory model (LSTM) */
    TF_FROZEN_MODEL_LSTM = 1,
};

/* preconditioner used with iterative linear solver */
enum pre_comp
{
    NONE_PC = 0,
    JACOBI_PC = 1,
    ICHOLT_PC = 2,
    ILUT_PC = 3,
    ILUTP_PC = 4,
    FG_ILUT_PC = 5,
    SAI_PC = 6,
};

/* method used to apply preconditioner for 2-side incomplete factorizations (ICHOLT, ILU) */
enum pre_app
{
    TRI_SOLVE_PA = 0,
    TRI_SOLVE_LEVEL_SCHED_PA = 1,
    TRI_SOLVE_GC_PA = 2,
    JACOBI_ITER_PA = 3,
};


/* atom types as pertains to hydrogen bonding */
enum hydrogen_bonding_atom_types
{
    NON_H_BONDING_ATOM = 0,
    H_ATOM = 1,
    H_BONDING_ATOM = 2,
};


typedef double real;
typedef real rvec[3];
typedef int ivec[3];
typedef real rtensor[3][3];


/* struct declarations, see definitions below for comments */
typedef struct global_parameters global_parameters;
typedef struct single_body_parameters single_body_parameters;
typedef struct two_body_parameters two_body_parameters;
typedef struct three_body_parameters three_body_parameters;
typedef struct three_body_header three_body_header;
typedef struct hbond_parameters hbond_parameters;
typedef struct four_body_parameters four_body_parameters;
typedef struct four_body_header four_body_header;
typedef struct reax_interaction reax_interaction;
typedef struct reax_atom reax_atom;
typedef struct simulation_box simulation_box;
typedef struct grid grid;
typedef struct reax_system reax_system;
typedef struct control_params control_params;
typedef struct thermostat thermostat;
typedef struct isotropic_barostat isotropic_barostat;
typedef struct flexible_barostat flexible_barostat;
typedef struct reax_timing reax_timing;
typedef struct simulation_data simulation_data;
typedef struct three_body_interaction_data three_body_interaction_data;
typedef struct near_neighbor_data near_neighbor_data;
typedef struct far_neighbor_data far_neighbor_data;
typedef struct hbond_data hbond_data;
typedef struct dDelta_data dDelta_data;
typedef struct dbond_data dbond_data;
typedef struct bond_order_data bond_order_data;
typedef struct bond_data bond_data;
typedef struct sparse_matrix sparse_matrix;
typedef struct reallocate_data reallocate_data;
typedef struct LR_data LR_data;
typedef struct cubic_spline_coef cubic_spline_coef;
typedef struct LR_lookup_table LR_lookup_table;
typedef struct static_storage static_storage;
typedef struct reax_list reax_list;
typedef struct output_controls output_controls;
typedef struct spuremd_handle spuremd_handle;


/* function pointer for calculating a bonded interaction */
typedef void (*interaction_function)( reax_system *, control_params *,
        simulation_data *, static_storage *, reax_list **, output_controls * );
/* function pointer for calculating pairwise atom distance */
typedef real (*atom_distance_function)( simulation_box const * const,
        rvec, rvec, ivec, ivec, ivec, rvec );
/* function pointer for updating atom positions */
typedef void (*update_atom_position_function)( rvec, rvec, ivec, simulation_box const * const );
#if defined(TEST_FORCES)
/* function pointers for printed bonded interactions */
typedef void (*print_interaction)( reax_system *, control_params *, simulation_data *,
        static_storage *, reax_list **, output_controls * );
#endif
/* function pointer for evolving the atomic system (i.e., updating the positions)
 * given the pre-computed forces from the prescribed interactions */
typedef int (*evolve_function)( reax_system * const, control_params * const,
        simulation_data * const, static_storage * const, reax_list ** const,
        output_controls * const );
/* function pointer for a callback function to be triggered after
 * completion of a simulation step -- useful for, e.g., the Python wrapper */
typedef void (*callback_function)( int, reax_atom *, simulation_data * );
/* function pointer for writing trajectory file header */
typedef int (*write_header_function)( reax_system *, control_params *,
        static_storage *, output_controls * );
/* function pointer for apending a frame to the trajectory file */
typedef int (*append_traj_frame_function)( reax_system *, control_params *,
        simulation_data *, static_storage *, reax_list **, output_controls * );
/* function pointer for writing to the trajectory file */
typedef int (*write_function)( FILE *, const char *, ... );


/* struct definitions */
/* Force field global parameters mapping
 * (contained in section 1 of file):
 *
 * l[0]  = p_boc1
 * l[1]  = p_boc2
 * l[2]  = p_coa2
 * l[3]  = N/A
 * l[4]  = N/A
 * l[5]  = N/A
 * l[6]  = p_ovun6
 * l[7]  = N/A
 * l[8]  = p_ovun7
 * l[9]  = p_ovun8
 * l[10] = N/A
 * l[11] = N/A
 * l[12] = N/A
 * l[13] = N/A
 * l[14] = p_val6
 * l[15] = p_lp1
 * l[16] = p_val9
 * l[17] = p_val10
 * l[18] = N/A
 * l[19] = p_pen2
 * l[20] = p_pen3
 * l[21] = p_pen4
 * l[22] = N/A
 * l[23] = p_tor2
 * l[24] = p_tor3
 * l[25] = p_tor4
 * l[26] = N/A
 * l[27] = p_cot2
 * l[28] = p_vdW1
 * l[29] = v_par30
 * l[30] = p_coa4
 * l[31] = p_ovun4
 * l[32] = p_ovun3
 * l[33] = p_val8
 * l[34] = b_s_acks2 (ACKS2 bond softness)
 * l[35] = N/A
 * l[36] = N/A
 * l[37] = N/A
 * l[38] = p_coa3
 * */
struct global_parameters
{
    /* num. global parameters in the force field paramater file for this simulation */
    int n_global;
    /* max. num. global parameters in the force field paramater file across all simulations */
    int max_n_global;
    /* parameters, see list in comment above for mapping */
    real* l;
    /* van der Waals interaction type */
    int vdw_type;
};


struct single_body_parameters
{
    /* Line one in field file */
    /* two character atom name */
    char name[15];
    /**/
    real r_s;
    /* valency of the atom */
    real valency;
    /* mass of atom, in unified atomic mass units (daltons) */
    real mass;
    /**/
    real r_vdw;
    /**/
    real epsilon;
    /**/
    real gamma;
    /**/
    real r_pi;
    /**/
    real valency_e;
    /**/
    real nlp_opt;

    /* Line two in field file */
    /**/
    real alpha;
    /**/
    real gamma_w;
    /**/
    real valency_boc;
    /**/
    real p_ovun5;
    /**/
    real chi;
    /**/
    real eta;
    /* info related to hydrogen bonding
     * (values correspond to hydrogen_bonding_atom_types enum above):
     *  0: non-hydrogen bonding atom
     *  1: H atom
     *  2: hydrogen bonding atom (e.g., O, S, P, N) */
    int p_hbond;

    /* Line three in field file */
    /**/
    real r_pi_pi;
    /**/
    real p_lp2;
    /**/
    real b_o_131;
    /**/
    real b_o_132;
    /**/
    real b_o_133;
    /* bond softness for ACKS2 */
    real b_s_acks2;

    /* Line four in the field file */
    /**/
    real p_ovun2;
    /**/
    real p_val3;
    /**/
    real valency_val;
    /**/
    real p_val5;
    /**/
    real rcore2;
    /**/
    real ecore2;
    /**/
    real acore2;
};


/* Two Body Parameters */
struct two_body_parameters
{
    /* Bond Order parameters */
    /**/
    real p_bo1;
    /**/
    real p_bo2;
    /**/
    real p_bo3;
    /**/
    real p_bo4;
    /**/
    real p_bo5;
    /**/
    real p_bo6;
    /**/
    real r_s;
    /**/
    real r_p;
    /* r_o distances in BO formula */
    real r_pp;
    /**/
    real p_boc3;
    /**/
    real p_boc4;
    /**/
    real p_boc5;

    /* Bond Energy parameters */
    /**/
    real p_be1;
    /**/
    real p_be2;
    /**/
    real De_s;
    /**/
    real De_p;
    /**/
    real De_pp;

    /* Over/Under coordination parameters */
    /**/
    real p_ovun1;

    /* van der Waals interaction parameters */
    /**/
    real D;
    /**/
    real alpha;
    /**/
    real r_vdW;
    /**/
    real gamma_w;
    /**/
    real rcore;
    /**/
    real ecore;
    /**/
    real acore;

    /* electrostatic parameters */
    /* note: this parameter is gamma^-3 and not gamma */
    real gamma;

    /**/
    real v13cor;
    /**/
    real ovc;
};


/* 3-body parameters */
struct three_body_parameters
{
    /* valence angle */
    real theta_00;
    /**/
    real p_val1;
    /**/
    real p_val2;
    /**/
    real p_val4;
    /**/
    real p_val7;
    /* penalty */
    real p_pen1;
    /* 3-body conjugation */
    real p_coa1;
};


struct three_body_header
{
    /**/
    int cnt;
    /**/
    three_body_parameters prm[MAX_3BODY_PARAM];
};


/* hydrogen bond parameters in force field parameters file */
struct hbond_parameters
{
    /* ideal H-bond distance */
    real r0_hb;
    /**/
    real p_hb1;
    /**/
    real p_hb2;
    /**/
    real p_hb3;
    /* TRUE if parameters are set for this triplet of atom types, FALSE otherwise */
    int is_valid;
};


/* 4-body parameters */
struct four_body_parameters
{
    /**/
    real V1;
    /**/
    real V2;
    /**/
    real V3;

    /* torsion angle */
    real p_tor1;

    /* 4-body conjugation */
    real p_cot1;
};


struct four_body_header
{
    /**/
    int cnt;
    /**/
    four_body_parameters prm[MAX_4BODY_PARAM];
};


struct reax_interaction
{
    /* num. atom types in force field parameter file for this simulation */
    int num_atom_types;
    /* max. num. atom types in force field parameter file across all simulations */
    int max_num_atom_types;
    /* */
    global_parameters gp;
    /* */
    single_body_parameters *sbp;
    /* */
    two_body_parameters **tbp;
    /* */
    three_body_header ***thbp;
    /* */
    hbond_parameters ***hbp;
    /* */
    four_body_header ****fbp;
};


struct reax_atom
{
    /* integer representation of element type of this atom */
    int type;
    /* TRUE if the atom is a dummy atom, FALSE otherwise
     * Note: dummy atoms _only_ participate in Coulomb interactions */
    int is_dummy;
    /* relative coordinates in terms of periodic images of the
     * simulation box which are used to track if this atom moves
     * between images between simulation steps which regenerate
     * the far neighbors list */
    ivec rel_map;
    /* string representation of element type of this atom */
    char name[9];
    /* position of this atom (3D space), in Angstroms (1E-10 m) */
    rvec x;
    /* velocity of this atom, in Angstroms / ps */
    rvec v;
    /* force acting on this atom, in AMU * Angstroms / ps^2 */
    rvec f;
    /* charge on this atom, in Coulombs */
    real q;
#if defined(QMMM)
    /* TRUE if the atom is in the QM region, FALSE otherwise (atom in MM region) */
    int qmmm_mask;
    /* inital charge on this atom, in Coulombs */
    real q_init;
#endif
};


struct simulation_box
{
    /* current volume of the simulation box, in Angstroms^3 */
    real volume;
    /* smallest coordinates within the simulation box
     * (typically (0.0, 0.0, 0.0) by convention) */
    rvec min;
    /* largest coordinates within the simulation box
     * (typically (l_x, l_y, l_z) by convention,
     *  where l_* is the length of the simulation box
     *  along that dimension) */
    rvec max;
    /* lengths of the simulation box along each dimension */
    rvec box_norms;
    /* box proportions, used for isotrophic NPT */
    rvec side_prop;
    /* 3D coordinate vectors defining the box dimensions */
    rtensor box;
    /* inverse of tensor containing box coordinate vectors */
    rtensor box_inv;
    /* copy of previous simulation box tensor,
     * used for isotrphic NPT */
    rtensor old_box;
    /* tensor containing translation vectors
     * used in calculations for non-orthogonal simulation boxes */
    rtensor trans;
    /* inverse of tensor containing translation vectors */
    rtensor trans_inv;
    /**/
    rtensor g;
};


struct grid
{
    /* 0 if struct members are NOT allocated, 1 otherwise */
    int allocated;
    /* max. num. of atoms that can be binned within a grid cell;
     * used for memory allocation purposes */
    int max_atoms;
    /* max. num. of neighbors for a grid cell;
     * used for memory allocation purposes */
    int max_nbrs;
    /* lengths of each dimesion of a grid cell */
    real cell_size;
    /**/
    ivec spread;
    /* num. of cells along each dimension for entire grid for this simulation */
    ivec ncell;
    /* max. num. of cells along each dimension for entire grid across all simulations */
    ivec ncell_max;
    /* lengths of each dimesion of a grid cell */
    rvec len;
    /* multiplicative inverse of length of each dimesion of a grid cell */
    rvec inv_len;
    /* binned atom list, where first 3 dimensions are indexed using 0-based grid cell
     * coordinates */
    int **** atoms;
    /**/
    int *** top;
    /**/
    int *** mark;
    /**/
    int *** start;
    /**/
    int *** end;
    /* list of neighboring grid cell positions for a specific grid cell, where the first
     * 3 dimensions are indexed using the 0-based grid cell coordinates */
    ivec **** nbrs;
    /* list of distances between the closest points of neighboring grid cells for a specific grid cell,
     * where the first 3 dimensions are indexed using the 0-based grid cell coordinates */
    rvec **** nbrs_cp;
};


struct reax_system
{
    /* 0 if struct members for first phase (file I/O) are NOT allocated, 1 otherwise */
    int prealloc_allocated;
    /* 0 if struct members are NOT allocated, 1 otherwise */
    int ffield_params_allocated;
    /* FALSE if struct members are NOT allocated, TRUE otherwise */
    int allocated;
    /* number of local (non-periodic image) atoms for the current simulation */
    int N;
#if defined(QMMM)
    /* number of local (non-periodic image) QM atoms for the current simulation in QMMM mode */
    int N_qm;
    /* number of local (non-periodic image) MM atoms for the current simulation in QMMM mode */
    int N_mm;
#endif
    /* max. number of local (non-periodic image) atoms across all simulations */
    int N_max;
    /* dimension of the N x N sparse charge method matrix H */
    int N_cm;
    /* max. dimension of the N x N sparse charge method matrix H across all simulations */
    int N_cm_max;
    /* molecular charge constraints
     * NOTE: these constraints are currently only supported in BGF files using EEM */
    real *molec_charge_constraints;
    /* molecular charge constraints encoded as pairs of 1-based atom numbers indicating a range of atoms
     * NOTE: these constraints are currently only supported in BGF files using EEM */
    int *molec_charge_constraint_ranges;
    /* num. of charge constraints on groups of atoms (i.e., molecules) */
    unsigned int num_molec_charge_constraints;
    /* max. num. of charge constraints on groups of atoms (i.e., molecules) */
    unsigned int max_num_molec_charge_constraints;
    /* num. of custom charge constraints */
    unsigned int num_custom_charge_constraints;
    /* max. num. of custom charge constraints */
    unsigned int max_num_custom_charge_constraints;
    /* num. of custom charge constraint entries (atom indices, coefficients) */
    unsigned int num_custom_charge_constraint_entries;
    /* max. num. of custom charge constraint entries (atom indices, coefficients) */
    unsigned int max_num_custom_charge_constraint_entries;
    /* entry counts for each custom charge constraint */
    int *custom_charge_constraint_count;
    /* indices for custom charge constraint info (atom indices, coefficients) */
    int *custom_charge_constraint_start;
    /* 1-based atom numbers used to identify entries for each custom charge constraint */
    int *custom_charge_constraint_atom_index;
    /* coefficients of entries for each custom charge constraint */
    real *custom_charge_constraint_coeff;
    /* right-hand side (RHS) constant values for each custom charge constraint */
    real *custom_charge_constraint_rhs;
    /* atomic interaction parameters */
    reax_interaction reax_param;
    /* grid specifying domain (i.e., spatial) decompisition
     * of atoms within simulation box */
    simulation_box box;
    /* grid structure used for binning atoms and tracking neighboring bins */
    grid g;
    /* collection of atomic info. */
    reax_atom *atoms;
    /* num. bonds per atom */
    int *bonds;
    /* num. hydrogen bonds per atom */
    int *hbonds;
    /* num. matrix entries per row */
    int *cm_entries;
};


/* Simulation control parameters not related to the system */
struct control_params
{
    /* string represention for this simulation */
    char sim_name[MAX_STR];
    /* string printed in header block of trajectory files in CUSTOM format*/
    char restart_from[MAX_STR];
    /* indicates whether this simulation is from a restart geometry file:
     * 0 = non-restarted run, 1 = restarted run */
    int restart;
    /* controls velocity initializtion for simulation with non-zero
     * starting temperatures: 0 = zero vectors, 1 = random vectors */
    int random_vel;
    /* controls additional atomic position translation during
     * simulation initializtion: 0 = no further translation,
     * 1 = translate positions such that the center of mass
     * is at the center of the simulation box, 2 = translate
     * positions such that the center of mass is at the origin
     * (smallest point) in the simulation box */
    int reposition_atoms;
    /* ensemble values:
     * 0 : microcanonical ensemble (NVE)
     * 1 : Berendsen NVT (bNVT)
     * 2 : Nose-Hoover NVT (nhNVT)
     * 3 : Parrinello-Rehman-Nose-Hoover semi-isotropic NPT (sNPT)
     * 4 : Parrinello-Rehman-Nose-Hoover isotropic NPT (iNPT) 
     * 5 : Parrinello-Rehman-Nose-Hoover anisotropic NPT (aNPT) */
    int ensemble;
    /* number of simulation time steps */
    int nsteps;
    /* controls whether the simulation box has perdioic boundary:
     * 0 = no periodic boundaries, 1 = periodic boundaries */
    int periodic_boundaries;
    /* */
    int restrict_bonds;
    /* controls whether force computations should be calculated
     * exactly or estimated (tabulated): 0 = exact computation,
     * >0 = use a lookup table (computed using splines), where
     * the positive integer value controls number of entries in the table */
    int tabulate;
    /* simulation time step length (in ps) */
    real dt;
    /* number of simulation steps to elapse before
     * recomputing the verlet lists */
    int reneighbor;
    /* cutoff (in Angstroms) used for for constructing the
     * long range interaction Verlet list (a.k.a. far neighbor list) */
    real vlist_cut;
    /* bonded interaction cutoff (in Angstroms) */
    real bond_cut;
    /* nonbonded interaction cutoff (in Angstroms),
     * also used as upper taper radius */
    real nonb_cut;
    /* shorter version of nonbonded interaction cutoff (in Angstroms),
     * used for computing sparser charge matrix
     * for electrostatic interactions */
    real nonb_sp_cut;
    /* shorter version of nonbonded interaction cutoff (in Angstroms),
     * used for lower taper radius */
    real nonb_low;
    /* bond order cutoff (in Angstroms) */
    real bo_cut;
    /* three body (valence angle) cutoff (in Angstroms) */
    real thb_cut;
    /* hydrogen bonding cutoff (in Angstroms) */
    real hbond_cut;
    /* initial temperature of atomic system (in K) */
    real T_init;
    /* target (final) temperature of atomic system (in K) */
    real T_final;
    /* current temperature (in K) */
    real T;
    /* thermostat inertial damping constant (in K / s) */
    real Tau_T;
    /* control mode for thermostat,
     * 0: none, 1: step-wise, 2: constant slope */
    int T_mode;
    /* step-wise control mode: limit in change of temperature in thermostat per simulation step (in K)
     * constant slope control mode: scaler for change of temperature in thermostat per simulation step (in K) */
    real T_rate;
    /* step-wise thermostat control mode: simulation time to hold thermostat temperature constant (in ps),
     * constant slope control mode: used to scale rate of temperature change (inversely) */
    real T_freq;
    /**/
    real Tau_PT;
    /**/
    rvec P;
    /**/
    rvec Tau_P;
    /* mode for pressures calculations,
     * 0: both int and ext, 1: ext only, 2: int only */
    int press_mode;
    /**/
    real compressibility;
    /* 0: do not compute pressure, 1: compute pressure;
     * NOTE: not applicable to NPT-type ensembles (always compute) */
    int compute_pressure;
    /* frequency in simulation time steps to remove center of mass velocity */
    int remove_CoM_vel;
    /* format of the geometry input file */
    int geo_format;
    /**/
    int dipole_anal;
    /**/
    int freq_dipole_anal;
    /**/
    int diffusion_coef;
    /**/
    int freq_diffusion_coef;
    /**/
    int restrict_type;
    /* method for computing atomic charges */
    unsigned int charge_method;
    /* frequency (in terms of simulation time steps) at which to
     * re-compute atomic charge distribution */
    int charge_freq;
    /* iterative linear solver type */
    unsigned int cm_solver_type;
    /* system net charge */
    real cm_q_net;
    /* max. iterations for linear solver */
    unsigned int cm_solver_max_iters;
    /* max. iterations before restarting in specific solvers, e.g., GMRES(k) */
    unsigned int cm_solver_restart;
    /* error tolerance of solution produced by charge distribution
     * sparse iterative linear solver */
    real cm_solver_q_err;
    /* ratio used in computing sparser charge matrix,
     * between 0.0 and 1.0 */
    real cm_domain_sparsity;
    /* TRUE if enabled, FALSE otherwise */
    unsigned int cm_domain_sparsify_enabled;
    /* type of method used for computing the initial guess
     * for the iterative linear solver */
    unsigned int cm_init_guess_type;
    /* order of spline extrapolation used for computing initial guess
     * to linear solver */
    unsigned int cm_init_guess_extrap1;
    /* order of spline extrapolation used for computing initial guess
     * to linear solver */
    unsigned int cm_init_guess_extrap2;
    /* file name for the GraphDef (GD) model file 
     * when predicting solver initial guesses using Tensorflow */
    char cm_init_guess_gd_model[MAX_STR];
    /* window size for the long short-term memory model (LSTM)
     * when predicting solver initial guesses using Tensorflow */
    unsigned int cm_init_guess_win_size;
    /* preconditioner type for linear solver */
    unsigned int cm_solver_pre_comp_type;
    /* frequency (in terms of simulation time steps) at which to recompute
     * incomplete factorizations */
    unsigned int cm_solver_pre_comp_refactor;
    /* drop tolerance of incomplete factorization schemes (ILUT, ICHOLT, etc.)
     * used for preconditioning the iterative linear solver used in charge distribution */
    real cm_solver_pre_comp_droptol;
    /* num. of sweeps for computing preconditioner factors
     * in fine-grained iterative methods (FG-ICHOL, FG-ILU) */
    unsigned int cm_solver_pre_comp_sweeps;
    /* relative num. of non-zeros to charge matrix used to
     * compute the sparse approximate inverse preconditioner,
     * between 0.0 and 1.0 */
    real cm_solver_pre_comp_sai_thres;
    /* preconditioner application type */
    unsigned int cm_solver_pre_app_type;
    /* num. of iterations used to apply preconditioner via
     * Jacobi relaxation scheme (truncated Neumann series) */
    unsigned int cm_solver_pre_app_jacobi_iters;
    /* TRUE if polarization energy calculation is enabled, FALSE otherwise */
    unsigned int polarization_energy_enabled;
    /**/
    int molec_anal;
    /**/
    int freq_molec_anal;
    /* controls which bonds are printed to the trajectory file
     * if bond printing is enabled; this value is used as a cut-off
     * for comparing * against the bond order calculation value */
    real bg_cut;
    /**/
    int num_ignored;
    /**/
    int ignore[MAX_ATOM_TYPES];
    /* number of OpenMP threads to use during the simulation */
    int num_threads;
    /* TRUE if the num. OpenMP has bet set, FALSE otherwise */
    int num_threads_set;
    /* function pointers for bonded interactions */
    interaction_function intr_funcs[NUM_INTRS];
    /* function pointer for computing pairwise atom distance */
    atom_distance_function compute_atom_distance;
    /* function pointer for updating atom position */
    update_atom_position_function update_atom_position;
#if defined(TEST_FORCES)
    /* function pointers for printed bonded interactions */
    print_interaction print_intr_funcs[NUM_INTRS];
#endif
};


struct thermostat
{
    /* temperature scalar of the system at the current simulation step */
    real T;
    /* temperature tensor of the system at the current simulation step */
    rtensor Temp;
    /* for Nose-Hoover thermostat */
    real xi;
    /* for Nose-Hoover thermostat */
    real v_xi;
    /* for Nose-Hoover thermostat */
    real v_xi_old;
    /* for Nose-Hoover thermostat */
    real G_xi;
};


struct isotropic_barostat
{
    real P;
    real eps;
    real v_eps;
    real v_eps_old;
    real a_eps;
};


struct flexible_barostat
{
    rtensor P;
    real P_scalar;

    real eps;
    real v_eps;
    real v_eps_old;
    real a_eps;

    rtensor h0;
    rtensor v_g0;
    rtensor v_g0_old;
    rtensor a_g0;
};


struct reax_timing
{
    /* start time of event */
    real start;
    /* end time of event */
    real end;
    /* total elapsed time of event */
    real elapsed;
    /* total simulation time */
    real total;
    /* neighbor list generation time */
    real nbrs;
    /* force initialization time */
    real init_forces;
    /* bonded force calculation time */
    real bonded;
    /* non-bonded force calculation time */
    real nonb;
    /* atomic charge distribution calculation time */
    real cm;
    /**/
    real cm_sort_mat_rows;
    /**/
    real cm_solver_pre_comp;
    /**/
    real cm_solver_pre_app;
    /* num. of steps in iterative linear solver for charge distribution */
    int cm_solver_iters;
    /**/
    real cm_solver_spmv;
    /**/
    real cm_solver_vector_ops;
    /**/
    real cm_solver_orthog;
    /**/
    real cm_solver_tri_solve;
    /* time spent on last preconditioner computation */
    real cm_last_pre_comp;
    /* time lost for not refactoring */
    real cm_total_loss;
    /* solver time on last refactoring step */
    real cm_optimum;
    /* num. of retries in main sim. loop */
    int num_retries;
};


struct simulation_data
{
    /* integer ID uniquely identifying the simulation */
    int sim_id;
    /* current simulation step number (0-based) */
    int step;
    /* last simulation step number for restarted runs (0-based) */
    int prev_steps;
    /* elapsed time of the simulation, in fs */
    real time;
    /* total mass of the atomic system, in AMU */
    real M;
    /* multiplicative inverse of the total mass, in AMU^{-1} */
    real inv_M;
    /* Center of mass */
    rvec xcm;
    /* Center of mass velocity */
    rvec vcm;
    /* Center of mass force */
    rvec fcm;
    /* Angular momentum of CoM */
    rvec amcm;
    /* Angular velocity of CoM */
    rvec avcm;
    /* Translational kinetic energy of CoM */
    real etran_cm;
    /* Rotational kinetic energy of CoM */
    real erot_cm;
    /* Kinetic energy tensor */
    rtensor kinetic;
    /* Hydrodynamic virial */
    rtensor virial;

    /* total energy, in kcal / mol */
    real E_Tot;
    /* total kinetic energy, in kcal / mol */
    real E_Kin;
    /* total potential energy, in kcal / mol */
    real E_Pot;

    /* total bond energy, in kcal / mol */
    real E_BE;
    /* total over coordination energy, in kcal / mol */
    real E_Ov;
    /* total under coordination energy, in kcal / mol */
    real E_Un;
    /* total under coordination energy, in kcal / mol */
    real E_Lp;
    /* total valance angle energy, in kcal / mol */
    real E_Ang;
    /* total penalty energy, in kcal / mol */
    real E_Pen;
    /* total three body conjgation energy, in kcal / mol */
    real E_Coa;
    /* total Hydrogen bond energy, in kcal / mol */
    real E_HB;
    /* total torsional energy, in kcal / mol */
    real E_Tor;
    /* total four body conjugation energy, in kcal / mol */
    real E_Con;
    /* total van der Waals energy, in kcal / mol */
    real E_vdW;
    /* total electrostatics energy, in kcal / mol */
    real E_Ele;
    /* polarization energy, in kcal / mol */
    real E_Pol;

    /* number of degrees of freedom */
    real N_f;
    /**/
    rvec t_scale;
    /**/
    rtensor p_scale;
    /* thermostat for Nose-Hoover ensemble */
    thermostat therm;
    /**/
    isotropic_barostat iso_bar;
    /**/
    flexible_barostat flex_bar;
    /**/
    real inv_W;

#if defined(_OPENMP)
    /* local per thread virial contribution to pressure */
    rtensor *press_local;
#endif
    /* virial contribution to pressure */
    rtensor press;
    /* kinetic energy contribution to pressure */
    rtensor kin_press;
    /* total pressure */
    rvec tot_press;

    /* struct containing timing of various simulation functions */
    reax_timing timing;
};


struct three_body_interaction_data
{
    /**/
    int thb;
    /* index for the third body on the central atom's bond list */
    int pthb;
    /* valence angle, in degrees */
    real theta;
    /* cosine of the valence angle, in degrees */
    real cos_theta;
    /* derivative coefficient of the cosine valence angle term for atom i */
    rvec dcos_di;
    /* derivative coefficient of the cosine valence angle term for atom j */
    rvec dcos_dj;
    /* derivative coefficient of the cosine valence angle term for atom k */
    rvec dcos_dk;
};


struct near_neighbor_data
{
    /**/
    int nbr;
    /**/
    ivec rel_box;
//    rvec ext_factor;
    /**/
    real d;
    /**/
    rvec dvec;
};


struct far_neighbor_data
{
    /**/
    int nbr;
    /**/
    ivec rel_box;
//    rvec ext_factor;
    /**/
    real d;
    /**/
    rvec dvec;
};


struct hbond_data
{
    /* neighbor atom ID */
    int nbr;
    /* ??? */
    int scl;
    /* neighbor in far neighbor list */
    far_neighbor_data *ptr;
};


struct dDelta_data
{
    /**/
    int wrt;
    /**/
    rvec dVal;
};


struct dbond_data
{
    /**/
    int wrt;
    /**/
    rvec dBO;
    /**/
    rvec dBOpi;
    /**/
    rvec dBOpi2;
};


struct bond_order_data
{
    real BO;
    real BO_s;
    real BO_pi;
    real BO_pi2;
    real Cdbo;
    real Cdbopi;
    real Cdbopi2;
    real C1dbo;
    real C2dbo;
    real C3dbo;
    real C1dbopi;
    real C2dbopi;
    real C3dbopi;
    real C4dbopi;
    real C1dbopi2;
    real C2dbopi2;
    real C3dbopi2;
    real C4dbopi2;
    rvec dBOp;
    rvec dln_BOp_s;
    rvec dln_BOp_pi;
    rvec dln_BOp_pi2;
};


struct bond_data
{
    /**/
    int nbr;
    /**/
    int sym_index;
    /**/
    int dbond_index;
    /**/
    ivec rel_box;
//    rvec ext_factor;
    /**/
    real d;
    /**/
    rvec dvec;
    /**/
    bond_order_data bo_data;
};


/* sparse matrix struct in compressed row storage (crs) format
 * See, e.g.,
 *   http://netlib.org/linalg/html_templates/node91.html#SECTION00931100000000000000
 */
struct sparse_matrix
{
    /* 0 if struct members are NOT allocated, 1 otherwise */
    int allocated;
    /* active number of rows for this simulation */
    unsigned int n;
    /* max. number of rows across all simulations */
    unsigned int n_max;
    /* number of nonzeros (NNZ) ALLOCATED */
    unsigned int m;
    /* row pointer (last element contains ACTUAL NNZ) */
    unsigned int *start;
    /* column index for corresponding matrix entry */
    unsigned int *j;
    /* matrix entry */
    real *val;
};


/* used to determine if and how much space should be reallocated */
struct reallocate_data
{
    /* TRUE if far neighbor list needs
     * to be reallocated, FALSE otherwise */
    int far_nbrs;
    /* total num. of (max) far neighbors across all atoms */
    int total_far_nbrs;
    /* TRUE if charge matrix needs
     * to be reallocated, FALSE otherwise */
    int cm;
    /* total num. matrix entries (sum over max) */
    int total_cm_entries;
    /* TRUE if hbond list needs
     * to be reallocated, FALSE otherwise */
    int hbonds;
    /* total num. hydrogen bonds (sum over max) */
    int total_hbonds;
    /* TRUE if bond list needs
     * to be reallocated, FALSE otherwise */
    int bonds;
    /* total num. bonds (sum over max) */
    int total_bonds;
    /* TRUE if three body list needs
     * to be reallocated, FALSE otherwise */
    int thbody;
    /* total num. three body interactions */
    int total_thbodies;
    /**/
    int gcell_atoms;
};


struct LR_data
{
    real H;
    real e_vdW;
    real CEvd;
    real e_ele;
    real CEclmb;
};


struct cubic_spline_coef
{
    real a;
    real b;
    real c;
    real d;
};


struct LR_lookup_table
{
    real xmin;
    real xmax;
    int n;
    real dx;
    real inv_dx;
    real a;
    real m;
    real c;

    LR_data *y;
    cubic_spline_coef *H;
    cubic_spline_coef *vdW;
    cubic_spline_coef *CEvd;
    cubic_spline_coef *ele;
    cubic_spline_coef *CEclmb;
};


struct static_storage
{
    /* FALSE if struct members are NOT allocated, TRUE otherwise */
    int allocated;
    /* bond order related storage */
    real *total_bond_order;
    real *Deltap;
    real *Deltap_boc;
    real *Delta;
    real *Delta_lp;
    real *Delta_lp_temp;
    real *Delta_e;
    real *Delta_boc;
    real *dDelta_lp;
    real *dDelta_lp_temp;
    real *nlp;
    real *nlp_temp;
    real *Clp;
    real *vlpex;
    rvec *dDeltap_self;

    /* charge method storage */
    /* charge matrix */
    sparse_matrix H;
    /* charge matrix (full) */
    sparse_matrix H_full;
    /* sparser charge matrix */
    sparse_matrix H_sp;
    /* permuted charge matrix (graph coloring) */
    sparse_matrix H_p;
    /* sparsity pattern of charge matrix, used in
     * computing a sparse approximate inverse preconditioner */
    sparse_matrix H_spar_patt;
    /* sparsity pattern of charge matrix (full), used in
     * computing a sparse approximate inverse preconditioner */
    sparse_matrix H_spar_patt_full;
    /* sparse approximate inverse preconditioner */
    sparse_matrix H_app_inv;
    /* incomplete Cholesky or LU preconditioner */
    sparse_matrix L;
    /* incomplete Cholesky or LU preconditioner */
    sparse_matrix U;
    /* Jacobi preconditioner */
    real *Hdia_inv;
    /* row drop tolerences for incomplete Cholesky preconditioner */
    real *droptol;
    /* right-hand side vectors for the linear systems */
    real *b_s;
    real *b_t;
    real *b_prc;
    real *b_prm;
    /* initial guesses for solutions to the linear systems, in Coulombs */
    real **s;
    real **t;

    /* GMRES related storage */
    real *y;
    real *z;
    real *g;
    real *hc;
    real *hs;
    real **h;
    real **rn;
    real **v;

    /* CG, SDM, BiCGStab related storage */
    real *r;
    real *r_hat;
    real *d;
    real *q;
    real *q_hat;
    real *p;

    /* SpMV related storage */
#if defined(_OPENMP)
    real *b_local;
#endif

    /* Level scheduling related storage for applying ICHOLT/ILUT(P)/FG-ILUT
     * preconditioners */
    int levels_L;
    int levels_U;
    unsigned int *row_levels_L;
    unsigned int *level_rows_L;
    unsigned int *level_rows_cnt_L;
    unsigned int *row_levels_U;
    unsigned int *level_rows_U;
    unsigned int *level_rows_cnt_U;
    unsigned int *top;

    /* Graph coloring related storage for applying ICHOLT/ILUT(P)/FG-ILUT
     * preconditioners */
    unsigned int *color;
    unsigned int *to_color;
    unsigned int *conflict;
    unsigned int *conflict_cnt;
    unsigned int *recolor;
    unsigned int recolor_cnt;
    unsigned int *color_top;
    unsigned int *permuted_row_col;
    unsigned int *permuted_row_col_inv;

    /* Graph coloring related storage for applying ICHOLT/ILUT(P)/FG-ILUT
     * preconditioners OR any application with ILUTP */
    real *y_p;
    real *x_p;

    /* Jacobi iteration related storage for applying ICHOLT/ILUT(P)/FG-ILUT
     * preconditioners */
    real *Dinv_L;
    real *Dinv_U;
    real *Dinv_b;
    real *rp;
    real *rp2;

    /* permutation for ILUTP */
    unsigned int *perm_ilutp;

    /* num. hydrogen atoms for this simulation */
    int num_H;

    rvec *a; // used in integrators
    rvec *f_old;
    rvec *v_const;

    /* coefficient of dDelta for force calculations */
    real *CdDelta;

    /* coefficients of Taper function applied to Coulomb interactions */
    real Tap[8];

    int *mark;
    int *old_mark;  // storage for analysis
    rvec *x_old;

    /* storage space for bond restrictions */
    int *map_serials;
    int *orig_id;
    int *restricted;
    int **restricted_list;

#if defined(_OPENMP)
    /* local forces per thread */
    rvec *f_local;
#endif
    unsigned int temp_int_omp;
    real temp_real_omp;

    reallocate_data realloc;

    LR_lookup_table **LR;

#if defined(TEST_FORCES)
    /* Calculated on the fly in bond_orders.c */
    rvec *dDelta;
    rvec *f_ele;
    rvec *f_vdw;
    rvec *f_be;
    rvec *f_lp;
    rvec *f_ov;
    rvec *f_un;
    rvec *f_ang;
    rvec *f_coa;
    rvec *f_pen;
    rvec *f_hb;
    rvec *f_tor;
    rvec *f_con;
#endif
};


/* 1-D flattened list of lists used for various interaction types */
struct reax_list
{
    /* FALSE if struct members are NOT allocated, TRUE otherwise */
    int allocated;
    /* num. active entries (i.e., nested lists) in the list for this simulation */
    int n;
    /* max. num. entries (i.e., nested lists) in the list across all simulations */
    int n_max;
    /* total interactions across all entries which can be stored in the list */
    int total_intrs;
    /* starting position of interactions for an inner list */
    int *index;
    /* ending position of interactions for an inner list */
    int *end_index;
    /* interaction list (polymorphic via union dereference) */
//    union
//    {
//        /* typeless interaction list */
//        void *v;
//        /* three body interaction list */
//        three_body_interaction_data *three_body_list;
//        /* bond interaction list */
//        bond_data *bond_list;
//        /* bond interaction list */
//        dbond_data *dbo_list;
//        /* test forces interaction list */
//        dDelta_data *dDelta_list;
//        /* far neighbor interaction list */
//        far_neighbor_data *far_nbr_list;
//        /* near neighbor interaction list */
//        near_neighbor_data *near_nbr_list;
//        /* hydrogen bond interaction list */
//        hbond_data *hbond_list;
//    } select;
    /* typeless interaction list */
    void *v;
    /* three body interaction list */
    three_body_interaction_data *three_body_list;
    /* bond interaction list */
    bond_data *bond_list;
    /* bond interaction list */
    dbond_data *dbo_list;
    /* test forces interaction list */
    dDelta_data *dDelta_list;
    /* far neighbor interaction list */
    far_neighbor_data *far_nbr_list;
    /* near neighbor interaction list */
    near_neighbor_data *near_nbr_list;
    /* hydrogen bond interaction list */
    hbond_data *hbond_list;
};


struct output_controls
{
    /* FALSE if struct members are NOT allocated, TRUE otherwise */
    int allocated;
    /* trajectory file */
    FILE *trj;
    /* system-wide info file */
    FILE *out;
    /* potential file */
    FILE *pot;
    /* log file */
    FILE *log;
    /**/
    FILE *mol;
    /**/
    FILE *ign;
    /**/
    FILE *dpl;
    /**/
    FILE *drft;
    /**/
    FILE *pdb;
    /**/
    FILE *prs;

    int write_steps;
    int traj_compress;
    int traj_format;
    char traj_title[81];
    int atom_format;
    int bond_info;
    int angle_info;

    int restart_format;
    int restart_freq;
    /* simulation step freqency at which log files are written
     * (excluding trajectory and restart files) */
    int log_update_freq;

    /* trajectory-related function pointers */
    write_header_function write_header;
    append_traj_frame_function append_traj_frame;
    write_function write;

#if defined(TEST_ENERGY)
    FILE *ebond;
    FILE *elp;
    FILE *eov;
    FILE *eun;
    FILE *eval;
    FILE *epen;
    FILE *ecoa;
    FILE *ehb;
    FILE *etor;
    FILE *econ;
    FILE *evdw;
    FILE *ecou;
#endif

#if defined(TEST_FORCES)
    FILE *fbo;
    FILE *fdbo;
    FILE *fbond;
    FILE *flp;
    FILE *fatom;
    FILE *f3body;
    FILE *fhb;
    FILE *f4body;
    FILE *fnonb;
    FILE *ftot;
    FILE *ftot2;
#endif
};


/* Handle for working with an instance of the sPuReMD library */
struct spuremd_handle
{
    /* System info. */
    reax_system *system;
    /* Control parameters */
    control_params *control;
    /* Atomic simulation data */
    simulation_data *data;
    /* Internal workspace */
    static_storage *workspace;
    /* Reax interaction list */
    reax_list **lists;
    /* Output controls */
    output_controls *out_control;
    /* TRUE if file I/O for simulation output enabled, FALSE otherwise */
    int output_enabled;
    /* TRUE if reallocation is required due to num. atoms increasing
     * (this includes first simulation run), FALSE otherwise */
    int realloc;
    /* Callback for getting simulation state at the end of each time step */
    callback_function callback;
};


#endif
