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

#include "geo_tools.h"

#include "allocate.h"
#include "box.h"
#include "tool_box.h"
#include "vector.h"

#include <ctype.h>


// CUSTOM_BOXGEO: BOXGEO box_x box_y box_z  angle1 angle2 angle3
#define CUSTOM_BOXGEO_FORMAT " %s %lf %lf %lf %lf %lf %lf"
// CUSTOM ATOM: serial element name x y z
#define CUSTOM_ATOM_FORMAT " %d %s %s %lf %lf %lf"

/* PDB format :
http://www.rcsb.org/pdb/file_formats/pdb/pdbguide2.2/guide2.2_frame.html

COLUMNS        DATA  TYPE    FIELD        DEFINITION
-------------------------------------------------------------------------------------
 1 -  6        Record name   "ATOM  "
 7 - 11        Integer       serial       Atom  serial number.
13 - 16        Atom          name         Atom name.
17             Character     altLoc       Alternate location indicator.
18 - 20        Residue name  resName      Residue name.
22             Character     chainID      Chain identifier.
23 - 26        Integer       resSeq       Residue sequence number.
27             AChar         iCode        Code for insertion of residues.
31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
55 - 60        Real(6.2)     occupancy    Occupancy.
61 - 66        Real(6.2)     tempFactor   Temperature  factor.
77 - 78        LString(2)    element      Element symbol, right-justified.
79 - 80        LString(2)    charge       Charge  on the atom.

Example

         1         2         3         4         5         6         7         8
12345678901234567890123456789012345678901234567890123456789012345678901234567890
ATOM     32  N  AARG A  -3      11.281  86.699  94.383  0.50 35.88           N
ATOM     33  N  BARG A  -3      11.296  86.721  94.521  0.50 35.60           N
ATOM     34  CA AARG A  -3      12.353  85.696  94.456  0.50 36.67           C
ATOM     35  CA BARG A  -3      12.333  85.862  95.041  0.50 36.42           C
ATOM     36  C  AARG A  -3      13.559  86.257  95.222  0.50 37.37           C
ATOM     37  C  BARG A  -3      12.759  86.530  96.365  0.50 36.39           C
ATOM     38  O  AARG A  -3      13.753  87.471  95.270  0.50 37.74           O
ATOM     39  O  BARG A  -3      12.924  87.757  96.420  0.50 37.26           O
ATOM     40  CB AARG A  -3      12.774  85.306  93.039  0.50 37.25           C
ATOM     41  CB BARG A  -3      13.428  85.746  93.980  0.50 36.60           C
ATOM     42  CG AARG A  -3      11.754  84.432  92.321  0.50 38.44           C
ATOM     43  CG BARG A  -3      12.866  85.172  92.651  0.50 37.31           C
ATOM     44  CD AARG A  -3      11.698  84.678  90.815  0.50 38.51           C
ATOM     45  CD BARG A  -3      13.374  85.886  91.406  0.50 37.66           C
ATOM     46  NE AARG A  -3      12.984  84.447  90.163  0.50 39.94           N
ATOM     47  NE BARG A  -3      12.644  85.487  90.195  0.50 38.24           N
ATOM     48  CZ AARG A  -3      13.202  84.534  88.850  0.50 40.03           C
ATOM     49  CZ BARG A  -3      13.114  85.582  88.947  0.50 39.55           C
ATOM     50  NH1AARG A  -3      12.218  84.840  88.007  0.50 40.76           N
ATOM     51  NH1BARG A  -3      14.338  86.056  88.706  0.50 40.23           N
ATOM     52  NH2AARG A  -3      14.421  84.308  88.373  0.50 40.45           N
*/

/*
COLUMNS     DATA TYPE        FIELD         DEFINITION
--------------------------------------------------------------
 1 - 6      Record name      "HETATM"
 7 - 11     Integer          serial        Atom serial number.
13 - 16     Atom             name          Atom name.
17          Character        altLoc        Alternate location indicator.
18 - 20     Residue name     resName       Residue name.
22          Character        chainID       Chain identifier.
23 - 26     Integer          resSeq        Residue sequence number.
27          AChar            iCode         Code for insertion of residues.
31 - 38     Real(8.3)        x             Orthogonal coordinates for X.
39 - 46     Real(8.3)        y             Orthogonal coordinates for Y.
47 - 54     Real(8.3)        z             Orthogonal coordinates for Z.
55 - 60     Real(6.2)        occupancy     Occupancy.
61 - 66     Real(6.2)        tempFactor    Temperature factor.
77 - 78     LString(2)       element       Element symbol; right-justified.
79 - 80     LString(2)       charge        Charge on the atom.
*/

/*
COLUMNS       DATA TYPE       FIELD         DEFINITION
-------------------------------------------------------
1 -  6       Record name     "CONECT"
7 - 11       Integer         serial        Atom serial number
12 - 16      Integer         serial        Serial number of bonded atom
17 - 21      Integer         serial        Serial number of bonded atom
22 - 26      Integer         serial        Serial number of bonded atom
27 - 31      Integer         serial        Serial number of bonded atom
*/

/*
COLUMNS       DATA TYPE       FIELD         DEFINITION
----------------------------------------------------------
1 - 6        Record name     "CRYST1"
7 - 15       Real(9.3)       a             a (Angstroms)
16 - 24      Real(9.3)       b             b (Angstroms)
25 - 33      Real(9.3)       c             c (Angstroms)
34 - 40      Real(7.2)       alpha         alpha (degrees)
41 - 47      Real(7.2)       beta          beta (degrees)
48 - 54      Real(7.2)       gamma         gamma (degrees)
56 - 66      LString         sGroup        Space group
67 - 70      Integer         z             Z value
*/

#define PDB_ATOM_FORMAT "%6s%5d%4s%c%4s%c%4d%c%8s%8s%8s%6s%6s%4s%2s%2s\n"
#define PDB_ATOM_FORMAT_LENGTH (72)
#define PDB_HETATM_FORMAT "%6s%5d%4s%c%4s%c%4d%c%8s%8s%8s%6s%6s%2s%2s\n"
#define PDB_CONECT_FORMAT "%6s%5d%5d%5d%5d%5d\n"
#define PDB_CRYST1_FORMAT "%6s%9s%9s%9s%7s%7s%7s%11s%4s\n"

#define PDB_ATOM_FORMAT_O "%6s%5d %4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n"
#define PDB_CRYST1_FORMAT_O "%6s%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%11s%4d\n"

#define BGF_ATOM_FORMAT "%6s %5s %5s %3s %c %5s%10s%10s%10s %5s%3s%2s %8s"
#define BGF_CRYSTX_FORMAT "%8s%11s%11s%11s%11s%11s%11s"

#define BGF_ATOM_FORMAT_O "%6s %5d %-5s %3s %c %5s%10.5f%10.5f%10.5f %-5s%3d%2d %8.5f"


/* Parse geometry file to determine simulation box parameters.
 *
 * system: struct containing simulation box struct to initialize
 * fp: file pointer to open geometry file 
 * geo_format: type of geometry file
 * */
static int Read_Box_Info( reax_system *system, FILE *fp, int geo_format )
{
    int ret, cryst_len;
    char *cryst;
    char line[MAX_LINE];
    char descriptor[9];
    char s_a[12], s_b[12], s_c[12], s_alpha[12], s_beta[12], s_gamma[12];
    real box_x, box_y, box_z, alpha, beta, gamma;
    char s_group[12], s_zValue[12];

    ret = FAILURE;

    switch ( geo_format )
    {
        case PDB:
            cryst = "CRYST1";
            cryst_len = 6;
            break;

        case CUSTOM:
            cryst = "BOXGEO";
            cryst_len = 6;
            break;
    }

    switch ( geo_format )
    {
        case PDB:
            fseek( fp, 0, SEEK_SET );

            /* locate the box info line in the file, read it, and
             * initialize the simulation box */
            while ( fgets( line, MAX_LINE, fp ) )
            {
                if ( strncmp( line, cryst, cryst_len ) == 0 )
                {
                    sscanf( line, PDB_CRYST1_FORMAT,
                            descriptor, s_a, s_b, s_c,
                            s_alpha, s_beta, s_gamma, s_group, s_zValue );

                    /* compute full volume tensor from the angles */
                    Setup_Box( atof(s_a),  atof(s_b), atof(s_c),
                            atof(s_alpha), atof(s_beta), atof(s_gamma),
                            &system->box );

                    ret = SUCCESS;
                    break;
                }
            }
            break;

        case CUSTOM:
            if ( fscanf( fp, CUSTOM_BOXGEO_FORMAT,
                        descriptor, &box_x, &box_y, &box_z,
                        &alpha, &beta, &gamma ) != 7 )
            {
                fprintf( stderr, "[ERROR] reading geometry file failed\n" \
                         "  [INFO] reading simulation box info\n" );
                exit( INVALID_INPUT );
            }

            Setup_Box( box_x, box_y, box_z,
                    alpha, beta, gamma,
                    &system->box );

            ret = SUCCESS;
            break;

        default:
            fprintf( stderr, "[ERROR] Unknown geometry file file (%d). Terminating...\n",
                  geo_format );
            break;
    }
    if ( ferror( fp ) )
    {
        ret = FAILURE;
    }

    return ret;
}


/* Parse geometry file to determine the number of atoms.
 *
 * system: struct containing simulation box struct to initialize
 * fp: file pointer to open geometry file 
 * geo_format: type of geometry file
 * */
static int Count_Atoms( reax_system *system, FILE *fp, int geo_format )
{
    char line[MAX_LINE];
    int n;

    n = 0;

    fseek( fp, 0, SEEK_SET );

    /* increment number of atoms for each line denoting an atom desc */
    switch ( geo_format )
    {
        case PDB:
            /* scan the entire file, looking for lines which
             * start with the strings "ATOM" or "HETATM" */
            while ( fgets( line, MAX_LINE, fp ) )
            {
                if ( strncmp( line, "ATOM", 4 ) == 0
                        || strncmp( line, "HETATM", 6 ) == 0 )
                {
                    ++n;
                }
            }

            break;

        case CUSTOM:
            /* skip box info */
            if ( fgets( line, MAX_LINE, fp ) == NULL )
            {
                fprintf( stderr, "[ERROR] reading geometry file failed\n" \
                         "  [INFO] reading simulation box info\n" );
                exit( INVALID_INPUT );
            }

            /* second line contains integer count
             * of the number of atoms */
            if ( fscanf( fp, " %d", &n ) != 1 )
            {
                fprintf( stderr, "[ERROR] reading geometry file failed\n" \
                         "  [INFO] reading number of atoms\n" );
                exit( INVALID_INPUT );
            }

            break;

        default:
            fprintf( stderr, "[ERROR] Unknown geometry file file (%d). Terminating...\n",
                  geo_format );
            break;
    }

    assert( n > 0 );

#if defined(DEBUG)
    fprintf( stderr, "[INFO] Count_Atoms: " );
    fprintf( stderr, "N = %d\n", n );
#endif

    return n;
}


/* Parser for geometry file in free-form custom PuReMD format
 *
 * geo_file: filename for custom geometry file to parse
 * system: struct containing atom-related information
 * control: struct containing simulation parameters
 * data: struct containing information on active simulations
 * workspace: struct containing intermediate structures used for calculations
 */
void Read_Geo( const char * const geo_file, reax_system* system, control_params *control,
        simulation_data *data, static_storage *workspace )
{
    FILE *geo;
    int i, j, n, serial, top;
    rvec x;
    char element[3], name[9];
    reax_atom *atom;

    geo = sfopen( geo_file, "r", __FILE__, __LINE__ );

    if ( Read_Box_Info( system, geo, CUSTOM ) == FAILURE )
    {
        fprintf( stderr, "[ERROR] Read_Box_Info: no BOXGEO line found in the geo file!" );
        fprintf( stderr, " Terminating...\n" );
        exit( INVALID_GEO );
    }

    /* count atoms and allocate storage */
    n = Count_Atoms( system, geo, CUSTOM );
    if ( system->prealloc_allocated == FALSE || n > system->N_max )
    {
        PreAllocate_Space( system, control, workspace, (int) CEIL( SAFE_ZONE * n ) );
    }
    system->N = n;

    /* parse atom info lines (3rd line and beyond) */
    top = 0;
    for ( i = 0; i < system->N; ++i )
    {
        if ( fscanf( geo, CUSTOM_ATOM_FORMAT,
                    &serial, element, name, &x[0], &x[1], &x[2] ) != 6 )
        {
            fprintf( stderr, "[ERROR] reading geometry file failed\n" \
                     "  [INFO] reading atom info (entry %d)\n", i );
            exit( INVALID_INPUT );
        }

        for ( j = 0; j < sizeof(element) - 1; ++j )
        {
            element[j] = toupper( element[j] );
        }
        element[sizeof(element) - 1] = '\0';

        Fit_to_Periodic_Box( &system->box, x );

#if defined(DEBUG_FOCUS)
        fprintf( stderr, "[INFO] atom: id = %d, element = %s, name = %s, x = (%f, %f, %f)\n",
                 serial, element, name, x[0], x[1], x[2] );
#endif

        atom = &system->atoms[top];
        workspace->orig_id[i] = serial;
        atom->type = Get_Atom_Type( &system->reax_param, element, sizeof(element) );
        strncpy( atom->name, name, sizeof(atom->name) - 1 );
        atom->name[sizeof(atom->name) - 1] = '\0';
        rvec_Copy( atom->x, x );
        rvec_MakeZero( atom->v );
        rvec_MakeZero( atom->f );
        atom->q = 0.0;
            
        /* check for dummy atom */
        if ( strncmp( element, "X\0", 2 ) == 0 )
        {
           atom->is_dummy = TRUE;
        }
        else
        {
            atom->is_dummy = FALSE;            
        }		

        top++;
    }

    sfclose( geo, __FILE__, __LINE__ );
}


/* Parser for geometry file in fixed-form PDB format
 *
 * pdb_file: filename for PDB geometry file to parse
 * system: struct containing atom-related information
 * control: struct containing simulation parameters
 * data: struct containing information on active simulations
 * workspace: struct containing intermediate structures used for calculations
 */
void Read_PDB( const char * const pdb_file, reax_system* system, control_params *control,
        simulation_data *data, static_storage *workspace )
{
    int i, n, c, c1, pdb_serial, top;
    FILE *pdb;
    char **tmp;
    char *s, *s1;
    char descriptor[7], serial[6];
    char atom_name[5], res_name[4], res_seq[5];
    char s_x[9], s_y[9], s_z[9];
    char occupancy[7], temp_factor[7];
    char seg_id[5], element[3], charge[3];
//    char alt_loc, chain_id, icode;
    rvec x;
    reax_atom *atom;

    pdb = sfopen( pdb_file, "r", __FILE__, __LINE__ );

    Allocate_Tokenizer_Space( &s, MAX_LINE, &s1, MAX_LINE,
            &tmp, MAX_TOKENS, MAX_TOKEN_LEN );

    if ( Read_Box_Info( system, pdb, PDB ) == FAILURE )
    {
        fprintf( stderr, "[ERROR] Read_Box_Info: no CRYST line found in the pdb file!" );
        fprintf( stderr, " Terminating...\n" );
        exit( INVALID_GEO );
    }

    n = Count_Atoms( system, pdb, PDB );
    if ( system->prealloc_allocated == FALSE || n > system->N_max )
    {
        PreAllocate_Space( system, control, workspace, (int) CEIL( SAFE_ZONE * n ) );
    }
    system->N = n;

    /* reset file pointer to beginning of file
     * and parse the PDB file */
    fseek( pdb, 0, SEEK_SET );
    c = 0;
    c1 = 0;
    top = 0;
    s[0] = 0;

    while ( fgets( s, MAX_LINE, pdb ) )
    {
        /* read new line and tokenize it */
        strncpy( s1, s, MAX_LINE - 1 );
        s1[MAX_LINE - 1] = '\0';
        c1 = Tokenize( s, &tmp, MAX_TOKEN_LEN );

        /* parse new line */
        if ( strncmp( tmp[0], "ATOM", 4 ) == 0
                || strncmp( tmp[0], "HETATM", 6 ) == 0 )
        {
            if ( strncmp( tmp[0], "ATOM", 4 ) == 0 )
            {
                strncpy( descriptor, s1, sizeof(descriptor) - 1 );
                descriptor[sizeof(descriptor) - 1] = '\0';
                strncpy( serial, s1 + 6, sizeof(serial) - 1 );
                serial[sizeof(serial) - 1] = '\0';
                strncpy( atom_name, s1 + 12, sizeof(atom_name) - 1 );
                atom_name[sizeof(atom_name) - 1] = '\0';
//                alt_loc = s1[16];
                strncpy( res_name, s1 + 17, sizeof(res_name) - 1 );
                res_name[sizeof(res_name) - 1] = '\0';
//                chain_id = s1[21];
                strncpy( res_seq, s1 + 22, sizeof(res_seq) - 1 );
                res_seq[sizeof(res_seq) - 1] = '\0';
//                icode = s1[26];
                strncpy( s_x, s1 + 30, sizeof(s_x) - 1 );
                s_x[sizeof(s_x) - 1] = '\0';
                strncpy( s_y, s1 + 38, sizeof(s_y) - 1 );
                s_y[sizeof(s_y) - 1] = '\0';
                strncpy( s_z, s1 + 46, sizeof(s_z) - 1 );
                s_z[sizeof(s_z) - 1] = '\0';
                strncpy( occupancy, s1 + 54, sizeof(occupancy) - 1 );
                occupancy[sizeof(occupancy) - 1] = '\0';
                strncpy( temp_factor, s1 + 60, sizeof(temp_factor) - 1 );
                temp_factor[sizeof(temp_factor) - 1] = '\0';
                strncpy( seg_id, s1 + 72, sizeof(seg_id) - 1 );
                seg_id[sizeof(seg_id) - 1] = '\0';
                strncpy( element, s1 + 76, sizeof(element) - 1 );
                element[sizeof(element) - 1] = '\0';
                strncpy( charge, s1 + 78, sizeof(charge) - 1 );
                charge[sizeof(charge) - 1] = '\0';
            }
            else if ( strncmp( tmp[0], "HETATM", 6 ) == 0 )
            {
                strncpy( descriptor, s1, sizeof(descriptor) - 1 );
                descriptor[sizeof(descriptor) - 1] = '\0';
                strncpy( serial, s1 + 6, sizeof(serial) - 1 );
                serial[sizeof(serial) - 1] = '\0';
                strncpy( atom_name, s1 + 12, sizeof(atom_name) - 1 );
                atom_name[sizeof(atom_name) - 1] = '\0';
//                alt_loc = s1[16];
                strncpy( res_name, s1 + 17, sizeof(res_name) - 1 );
                res_name[sizeof(res_name) - 1] = '\0';
//                chain_id = s1[21];
                strncpy( res_seq, s1 + 22, sizeof(res_seq) - 1 );
                res_seq[sizeof(res_seq) - 1] = '\0';
//                icode = s1[26];
                strncpy( s_x, s1 + 30, sizeof(s_x) - 1 );
                s_x[sizeof(s_x) - 1] = '\0';
                strncpy( s_y, s1 + 38, sizeof(s_y) - 1 );
                s_y[sizeof(s_y) - 1] = '\0';
                strncpy( s_z, s1 + 46, sizeof(s_z) - 1 );
                s_z[sizeof(s_z) - 1] = '\0';
                strncpy( occupancy, s1 + 54, sizeof(occupancy) - 1 );
                occupancy[sizeof(occupancy) - 1] = '\0';
                strncpy( temp_factor, s1 + 60, sizeof(temp_factor) - 1 );
                temp_factor[sizeof(temp_factor) - 1] = '\0';
                strncpy( element, s1 + 76, sizeof(element) - 1 );
                element[sizeof(element) - 1] = '\0';
                strncpy( charge, s1 + 78, sizeof(charge) - 1 );
                charge[sizeof(charge) - 1] = '\0';
            }

            /* if the point is inside my_box, add it to my lists */
            Make_Point( sstrtod( s_x, __FILE__, __LINE__ ),
                    sstrtod( s_y, __FILE__, __LINE__ ),
                    sstrtod( s_z, __FILE__, __LINE__ ), &x );

            Fit_to_Periodic_Box( &system->box, x );

            if ( is_Inside_Box( &system->box, x ) )
            {
                /* store orig_id, type, name and coord info of the new atom */
                atom = &system->atoms[top];
                pdb_serial = (int) sstrtod( serial, __FILE__, __LINE__ );
                workspace->orig_id[top] = pdb_serial;

                strncpy( atom->name, atom_name, sizeof(atom->name) - 1 );
                atom->name[sizeof(atom->name) - 1] = '\0';
                Trim_Spaces( element, sizeof(element) );
                for ( i = 0; i < sizeof(element) - 1; ++i )
                {
                    element[i] = toupper( element[i] );
                }
                atom->type = Get_Atom_Type( &system->reax_param, element, sizeof(element) );
            
                /* check for dummy atom */
                if ( strncmp( element, "X\0", 2 ) == 0 )
                {
                    system->atoms[top].is_dummy = TRUE;
                }
                else
                {
                    system->atoms[top].is_dummy = FALSE;            
                }		

                rvec_Copy( atom->x, x );
                rvec_MakeZero( atom->v );
                rvec_MakeZero( atom->f );
                atom->q = 0;

#if defined(DEBUG_FOCUS)
                fprintf( stderr, "[INFO] atom: id = %d, name = %s, serial = %d, type = %d, ",
                        top, atom->name, pdb_serial, atom->type );
                fprintf( stderr, "x = (%7.3f, %7.3f, %7.3f)\n",
                        atom->x[0], atom->x[1], atom->x[2] );
#endif

                top++;
            }

            c++;
        }

        /* IMPORTANT: We do not check for the soundness of restrictions here.
         * When atom2 is on atom1's restricted list, and there is a restriction
         * on atom2, then atom1 has to be on atom2's restricted list, too.
         * However, we do not check if this is the case in the input file,
         * this is upto the user. */
        else if ( strncmp( tmp[0], "CONECT", 6 ) == 0 )
        {
            if ( control->restrict_bonds )
            {
                /* error check */
//                Check_Input_Range( c1 - 2, 0, MAX_RESTRICT, __FILE__, __LINE__ );

                /* read bond restrictions */
                // pdb_serial = sstrtol( tmp[1], __FILE__, __LINE__ );
                // if ( is_Valid_Serial( workspace->map_serials[ pdb_serial ] ) == TRUE )
                //   ratom = workspace->map_serials[ pdb_serial ];

                // workspace->restricted[ ratom ] = c1 - 2;
                // for ( i = 2; i < c1; ++i )
                //  {
                //    pdb_serial = sstrtol( tmp[i], __FILE__, __LINE__ );
                //    if ( is_Valid_Serial( workspace->map_serials[ pdb_serial ] ) == TRUE )
                //        workspace->restricted_list[ ratom ][ i-2 ] =
                //          workspace->map_serials[ pdb_serial ];
                //  }

                // fprintf( stderr, "restriction on %d:", ratom );
                // for( i = 0; i < workspace->restricted[ ratom ]; ++i )
                // fprintf( stderr, "  %d",
                //          workspace->restricted_list[ratom][i] );
                // fprintf( stderr, "\n" );
            }
        }

        /* clear previous input line */
        s[0] = 0;
        for ( i = 0; i < c1; ++i )
        {
            tmp[i][0] = 0;
        }
    }
    if ( ferror( pdb ) )
    {
        fprintf( stderr, "[ERROR] Unable to read PDB file. Terminating...\n" );
        exit( INVALID_INPUT );
    }

    sfclose( pdb, __FILE__, __LINE__ );

    Deallocate_Tokenizer_Space( &s, &s1, &tmp, MAX_TOKENS );
} 


/* PDB serials are written without regard to the order, we'll see if this
 * cause trouble, if so we'll have to rethink this approach
 * Also, we do not write connect lines yet.
 * */
void Write_PDB( reax_system* system, reax_list* bonds, simulation_data *data,
        control_params *control, static_storage *workspace, output_controls *out_control )
{
    int i; 
    char name[6];
    real alpha, beta, gamma;
    rvec x;
    reax_atom *p_atom;
    char fname[MAX_STR];
    FILE *pdb;

    /* Writing Box information */
    gamma = ACOS( (system->box.box[0][0] * system->box.box[1][0]
                + system->box.box[0][1] * system->box.box[1][1]
                + system->box.box[0][2] * system->box.box[1][2])
            / (system->box.box_norms[0] * system->box.box_norms[1]) );
    beta  = ACOS( (system->box.box[0][0] * system->box.box[2][0]
                + system->box.box[0][1] * system->box.box[2][1]
                + system->box.box[0][2] * system->box.box[2][2])
            / (system->box.box_norms[0] * system->box.box_norms[2]) );
    alpha = ACOS( (system->box.box[2][0] * system->box.box[1][0]
                + system->box.box[2][1] * system->box.box[1][1]
                + system->box.box[2][2] * system->box.box[1][2])
            / (system->box.box_norms[2] * system->box.box_norms[1]) );

    /* write header */
    snprintf( fname, sizeof(fname) - 1, "%.*s-%.*d.pdb",
            (int) MIN( strnlen(control->sim_name,
                    sizeof(fname) - snprintf(NULL, 0, "%d", data->step) - 6 ),
                sizeof(control->sim_name) ),
            control->sim_name, snprintf(NULL, 0, "%d", data->step), data->step );
    fname[sizeof(fname) - 1] = '\0';
    pdb = sfopen( fname, "w", __FILE__, __LINE__ );
    fprintf( pdb, PDB_CRYST1_FORMAT_O,
             "CRYST1",
             system->box.box_norms[0], system->box.box_norms[1],
             system->box.box_norms[2],
             RAD2DEG(alpha), RAD2DEG(beta), RAD2DEG(gamma), " ", 0 );

    /* write atom lines to file */
    for ( i = 0; i < system->N; i++ )
    {
        p_atom = &system->atoms[i];

        strncpy( name, p_atom->name, sizeof(name) - 1 );
        name[sizeof(name) - 1] = '\0';
        Trim_Spaces( name, sizeof(name) );

        memcpy( x, p_atom->x, 3 * sizeof(real) );
        Fit_to_Periodic_Box( &system->box, x );

        fprintf( pdb, PDB_ATOM_FORMAT_O,
                "ATOM  ", workspace->orig_id[i], name,
                ' ', "REX", ' ', 0, ' ', x[0], x[1], x[2],
                0.0, 0.0, name, "  " );
    }
    
    if ( ferror( pdb ) )
    {
        fprintf( stderr, "[ERROR] Unable to write PDB file. Terminating...\n" );
        exit( INVALID_INPUT );
    }

    sfclose( pdb, __FILE__, __LINE__ );
}


/* Parser for geometry files in BGF format
 *
 * bgf_file: filename for BGF file to parse
 * system: struct containing atom-related information
 * control: struct containing simulation parameters
 * data: struct containing information on active simulations
 * workspace: struct containing intermediate structures used for calculations
 * */
void Read_BGF( const char * const bgf_file, reax_system* system, control_params *control,
        simulation_data *data, static_storage *workspace )
{
    FILE *bgf;
    char **tokens;
    char *line, *backup;
    char descriptor[7], serial[6];
    char atom_name[6], res_name[4], res_seq[6];
    char s_x[11], s_y[11], s_z[11];
    char occupancy[4], temp_factor[3];
    char element[6], charge[9];
//    char chain_id;
    char s_a[12], s_b[12], s_c[12], s_alpha[12], s_beta[12], s_gamma[12];
    int i, n, num_mcc, atom_cnt, token_cnt, bgf_serial, ratom, crystx_found;
    rvec x;

    ratom = 0;
    crystx_found = FALSE;

    bgf = sfopen( bgf_file, "r", __FILE__, __LINE__ );

    Allocate_Tokenizer_Space( &line, MAX_LINE, &backup, MAX_LINE,
            &tokens, MAX_TOKENS, MAX_TOKEN_LEN );

    /* count number of atoms in the BGF file */
    n = 0;
    num_mcc = 0;
    line[0] = 0;

    while ( fgets( line, MAX_LINE, bgf ) )
    {
        tokens[0][0] = '\0';
        token_cnt = Tokenize( line, &tokens, MAX_TOKEN_LEN );

        if ( strncmp( tokens[0], "ATOM", 4 ) == 0
                || strncmp( tokens[0], "HETATM", 6 ) == 0 )
        {
            ++n;
        }
        else if ( strncmp( tokens[0], "MOLCHARGE", 9 ) == 0 )
        {
            ++num_mcc;
        }

        line[0] = 0;
    }
    if ( ferror( bgf ) )
    {
        fprintf( stderr, "[ERROR] Unable to read BGF file. Terminating...\n" );
        exit( INVALID_INPUT );
    }

    if ( system->prealloc_allocated == FALSE || n > system->N_max )
    {
        PreAllocate_Space( system, control, workspace, (int) CEIL( SAFE_ZONE * n ) );

        if ( system->prealloc_allocated == FALSE && num_mcc > 0 )
        {
            system->molec_charge_constraints = smalloc(
                    sizeof(real) * num_mcc, __FILE__, __LINE__ );
            system->molec_charge_constraint_ranges = smalloc(
                    sizeof(int) * 2 * num_mcc, __FILE__, __LINE__ );

            system->max_num_molec_charge_constraints = num_mcc;
        }
        else if ( num_mcc > system->max_num_molec_charge_constraints )
        {
            if ( system->max_num_molec_charge_constraints > 0 )
            {
                sfree( system->molec_charge_constraints, __FILE__, __LINE__ );
                sfree( system->molec_charge_constraint_ranges, __FILE__, __LINE__ );
            }

            system->molec_charge_constraints = smalloc(
                    sizeof(real) * num_mcc, __FILE__, __LINE__ );
            system->molec_charge_constraint_ranges = smalloc(
                    sizeof(int) * 2 * num_mcc, __FILE__, __LINE__ );

            system->max_num_molec_charge_constraints = num_mcc;
        }
    }
    system->N = n;
    system->num_molec_charge_constraints = num_mcc;
    num_mcc = 0;

#if defined(DEBUG_FOCUS)
    fprintf( stderr, "[INFO] num_atoms = %d, num_mcc = %d\n", system->N,
            system->num_molec_charge_constraints );
#endif

    for ( i = 0; i < MAX_ATOM_ID; ++i )
    {
        workspace->map_serials[i] = -1;
    }

    fseek( bgf, 0, SEEK_SET );
    atom_cnt = 0;
    token_cnt = 0;

    while ( fgets( line, MAX_LINE, bgf ) )
    {
        /* read new line and tokenize it */
        strncpy( backup, line, MAX_LINE - 1 );
        backup[MAX_LINE - 1] = '\0';
        token_cnt = Tokenize( line, &tokens, MAX_TOKEN_LEN );

        if ( strncmp( tokens[0], "CRYSTX", 6 ) == 0 )
        {
            if ( sscanf( backup, BGF_CRYSTX_FORMAT, descriptor,
                        s_a, s_b, s_c, s_alpha, s_beta, s_gamma ) != 7 )
            {
                fprintf( stderr, "[ERROR] reading geometry file failed\n" \
                         "  [INFO] reading simulation box info\n" );
                exit( INVALID_INPUT );
            }

            /* compute full volume tensor from the angles */
            Setup_Box( sstrtod( s_a, __FILE__, __LINE__ ),
                    sstrtod( s_b, __FILE__, __LINE__ ),
                    sstrtod( s_c, __FILE__, __LINE__ ),
                    sstrtod( s_alpha, __FILE__, __LINE__ ),
                    sstrtod( s_beta, __FILE__, __LINE__ ),
                    sstrtod( s_gamma, __FILE__, __LINE__ ),
                    &system->box );

            crystx_found = TRUE;
            break;
        }
    }

    if ( crystx_found == FALSE )
    {
        fprintf( stderr, "[ERROR] improperly formatted BGF file (no CRYSTX keyword found). Terminating...\n" );
        exit( INVALID_INPUT );
    }

    fseek( bgf, 0, SEEK_SET );

    while ( fgets( line, MAX_LINE, bgf ) )
    {
        /* read new line and tokenize it */
        strncpy( backup, line, MAX_LINE - 1 );
        backup[MAX_LINE - 1] = '\0';
        token_cnt = Tokenize( line, &tokens, MAX_TOKEN_LEN );

        /* process lines with atom info (i.e., begin with keywords HETATM or ATOM),
         * format: "%6s %5d %-5s %3s %c %5s%10.5f%10.5f%10.5f %-5s%3d%2d %8.5f"
         *
         * also, it's common to see the atom info format
         * inlined in MSI BGF files as follows:
         * FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5)
         * */
        if ( strncmp( tokens[0], "ATOM", 4 ) == 0
                || strncmp( tokens[0], "HETATM", 6 ) == 0 )
        {
            strncpy( descriptor, backup, sizeof(descriptor) - 1 );
            descriptor[sizeof(descriptor) - 1] = '\0';
            strncpy( serial, backup + 7, sizeof(serial) - 1 );
            serial[sizeof(serial) - 1] = '\0';
            strncpy( atom_name, backup + 13, sizeof(atom_name) - 1 );
            atom_name[sizeof(atom_name) - 1] = '\0';
            strncpy( res_name, backup + 19, sizeof(res_name) - 1 );
            res_name[sizeof(res_name) - 1] = '\0';
//            chain_id = backup[23];
            strncpy( res_seq, backup + 25, sizeof(res_seq) - 1 );
            res_seq[sizeof(res_seq) - 1] = '\0';
            strncpy( s_x, backup + 30, sizeof(s_x) - 1 );
            s_x[sizeof(s_x) - 1] = '\0';
            strncpy( s_y, backup + 40, sizeof(s_y) - 1 );
            s_y[sizeof(s_y) - 1] = '\0';
            strncpy( s_z, backup + 50, sizeof(s_x) - 1 );
            s_z[sizeof(s_x) - 1] = '\0';
            strncpy( element, backup + 61, sizeof(element) - 1 );
            element[sizeof(element) - 1] = '\0';
            strncpy( occupancy, backup + 66, sizeof(occupancy) - 1 );
            occupancy[sizeof(occupancy) - 1] = '\0';
            strncpy( temp_factor, backup + 69, sizeof(temp_factor) - 1 );
            temp_factor[sizeof(temp_factor) - 1] = '\0';
            strncpy( charge, backup + 72, sizeof(charge) - 1 );
            charge[sizeof(charge) - 1] = '\0';

//            if ( sscanf( backup, BGF_ATOM_FORMAT, descriptor, serial, atom_name,
//                    res_name, &chain_id, res_seq, s_x, s_y, s_z, element,
//                    occupancy, temp_factor, charge ) != 13 )
//            {
//                fprintf( stderr, "[ERROR] reading geometry file failed\n" );
//                fprintf( stderr, "  [INFO] reading atom info (entry %d)\n", i );
//                exit( INVALID_INPUT );
//            }

            /* add to mapping */
            bgf_serial = sstrtod( serial, __FILE__, __LINE__ );
            Check_Input_Range( bgf_serial, 0, MAX_ATOM_ID, __FILE__, __LINE__ );
            workspace->map_serials[ bgf_serial ] = atom_cnt;
            workspace->orig_id[ atom_cnt ] = bgf_serial;

            /* copy atomic positions */
            x[0] = sstrtod( s_x, __FILE__, __LINE__ );
            x[1] = sstrtod( s_y, __FILE__, __LINE__ );
            x[2] = sstrtod( s_z, __FILE__, __LINE__ );

            Fit_to_Periodic_Box( &system->box, x );

            system->atoms[atom_cnt].x[0] = x[0];
            system->atoms[atom_cnt].x[1] = x[1];
            system->atoms[atom_cnt].x[2] = x[2];

            /* atom name and type */
            strncpy( system->atoms[atom_cnt].name, atom_name,
                    sizeof(system->atoms[atom_cnt].name) - 1 );
            system->atoms[atom_cnt].name[sizeof(system->atoms[atom_cnt].name) - 1] = '\0';
            Trim_Spaces( element, sizeof(element) );
            for ( i = 0; i < sizeof(element) - 1; ++i )
            {
                element[i] = toupper( element[i] );
            }
            system->atoms[atom_cnt].type =
                Get_Atom_Type( &system->reax_param, element, sizeof(element) );
            
            /* check for dummy atom */
            if ( strncmp( element, "X\0", 2 ) == 0 )
            {
                system->atoms[atom_cnt].is_dummy = TRUE;
            }
            else
            {
                system->atoms[atom_cnt].is_dummy = FALSE;            
            }		

#if defined(DEBUG_FOCUS)
            fprintf( stderr,
                    "[INFO] atom_cnt = %5d, atom_type = %3d, x = (%10.5f,%10.5f,%10.5f), q = %10.5f, occ = %s, temp = %s, res_name = %4s, element = %s\n",
                    atom_cnt, system->atoms[ atom_cnt ].type,
                    system->atoms[ atom_cnt ].x[0],
                    system->atoms[ atom_cnt ].x[1], system->atoms[ atom_cnt ].x[2],
                    system->atoms[ atom_cnt ].q, occupancy, temp_factor,
                    res_name, element );
#endif

            atom_cnt++;
        }
        else if ( strncmp( tokens[0], "CONECT", 6 ) == 0 )
        {
            if ( control->restrict_bonds )
            {
                /* check number of restrictions */
                Check_Input_Range( token_cnt - 2, 0, MAX_RESTRICT, __FILE__, __LINE__ );

                /* read bond restrictions */
                bgf_serial = sstrtol( tokens[1], __FILE__, __LINE__ );
                if ( is_Valid_Serial( workspace->map_serials[ bgf_serial ] ) == TRUE )
                {
                    ratom = workspace->map_serials[ bgf_serial ];
                }

                workspace->restricted[ ratom ] = token_cnt - 2;
                for ( i = 2; i < token_cnt; ++i )
                {
                    bgf_serial = sstrtol( tokens[i], __FILE__, __LINE__ );
                    if ( is_Valid_Serial( workspace->map_serials[ bgf_serial ] ) == TRUE )
                    {
                        workspace->restricted_list[ratom][i - 2] =
                            workspace->map_serials[bgf_serial];
                    }
                }
            }
        }
        else if ( strncmp( tokens[0], "MOLCHARGE", 9 ) == 0 )
        {
            assert( token_cnt == 4 );

            system->molec_charge_constraint_ranges[2 * num_mcc] = sstrtol( tokens[1], __FILE__, __LINE__ );
            system->molec_charge_constraint_ranges[2 * num_mcc + 1] = sstrtol( tokens[2], __FILE__, __LINE__ );
            system->molec_charge_constraints[num_mcc] = sstrtod( tokens[3], __FILE__, __LINE__ );

#if defined(DEBUG_FOCUS)
            fprintf( stderr,
                    "[INFO] num_mcc = %d, mcc = %f, mcc_range = (%d, %d)\n", num_mcc,
                    system->molec_charge_constraints[num_mcc],
                    system->molec_charge_constraint_ranges[2 * num_mcc],
                    system->molec_charge_constraint_ranges[2 * num_mcc + 1] );
#endif

            ++num_mcc;
        }

        /* clear previous input line */
        line[0] = '\0';

        for ( i = 0; i < token_cnt; ++i )
        {
            tokens[i][0] = '\0';
        }
    }
    if ( ferror( bgf ) )
    {
        fprintf( stderr, "[ERROR] Unable to read BGF file. Terminating...\n" );
        exit( INVALID_INPUT );
    }

    Deallocate_Tokenizer_Space( &line, &backup, &tokens, MAX_TOKENS );

    sfclose( bgf, __FILE__, __LINE__ );
}
