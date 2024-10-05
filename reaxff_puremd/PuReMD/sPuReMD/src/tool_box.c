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

#include "tool_box.h"

#include <stdlib.h>
#include <ctype.h>
#include <limits.h>
#include <errno.h>
#include <time.h>

/* base 10 for result of string-to-integer conversion */
#define INTBASE (10)


/************** taken from box.c **************/
/* Applies transformation on atomic position between Cartesian and
 * Triclinic coordinates based on the value of flag
 * 
 * Inputs:
 *  x1: position to be transformed
 *  box: struct containing simulation box parameters
 *  flag: -1 for Cartesian -> Triclinic, +1 for reverse transformation
 *
 * Outputs:
 *  x2: transformed position
 * */
void Transform( rvec x1, simulation_box *box, int flag, rvec x2 )
{
    int i, j;
    real tmp;

    if ( flag == 1 )
    {
        for ( i = 0; i < 3; i++ )
        {
            tmp = 0.0;

            for ( j = 0; j < 3; j++ )
            {
                tmp += box->trans[i][j] * x1[j];
            }

            x2[i] = tmp;
        }
    }
    else if ( flag == -1 )
    {
        for ( i = 0; i < 3; i++ )
        {
            tmp = 0.0;

            for ( j = 0; j < 3; j++ )
            {
                tmp += box->trans_inv[i][j] * x1[j];
            }

            x2[i] = tmp;
        }
    }
}


/* Applies transformation on atomic position between Cartesian and
 * Triclinic coordinates based on the value of flag, and scales
 * the result to be within a simulation box with unit lengths
 * 
 * Inputs:
 *  x1: position to be transformed
 *  box: struct containing simulation box parameters
 *  flag: -1 for Cartesian -> Triclinic, +1 for reverse transformation
 *
 * Outputs:
 *  x2: transformed position
 * */
void Transform_to_UnitBox( rvec x1, simulation_box *box, int flag, rvec x2 )
{
    Transform( x1, box, flag, x2 );

    x2[0] /= box->box_norms[0];
    x2[1] /= box->box_norms[1];
    x2[2] /= box->box_norms[2];
}


/* Check and remap (if necessary) an atom position to fall
 * within the boundaries of a periodic simulation box where
 * the boundaries are [0, d_i) with d_i being the length
 * of the simulation box in a particular dimension */
void Fit_to_Periodic_Box( simulation_box *box, rvec p )
{
    int i;

    for ( i = 0; i < 3; ++i )
    {
        if ( p[i] < box->min[i] )
        {
            /* handle lower coords */
            while ( p[i] < box->min[i] )
            {
                p[i] += box->box_norms[i];
            }
        }
        else if ( p[i] >= box->max[i] )
        {
            /* handle higher coords */
            while ( p[i] >= box->max[i] )
            {
                p[i] -= box->box_norms[i];
            }
        }
    }
}


/* determine whether point p is inside the box */
/* assumes orthogonal box */
int is_Inside_Box( simulation_box *box, rvec p )
{
    int ret = TRUE;

    if ( p[0] < box->min[0] || p[0] >= box->max[0]
            || p[1] < box->min[1] || p[1] >= box->max[1]
            || p[2] < box->min[2] || p[2] >= box->max[2] )
    {
        ret = FALSE;
    }

    return ret;
}


/************** from geo_tools.c *****************/
void Make_Point( real x, real y, real z, rvec* p )
{
    (*p)[0] = x;
    (*p)[1] = y;
    (*p)[2] = z;
}


int is_Valid_Serial( int serial )
{
    if ( serial < 0 )
    {
        fprintf( stderr, "[ERROR] CONECT line includes invalid serial number %d.\n", serial );
        fprintf( stderr, "[ERROR] Please correct the input file. Terminating...\n" );
        exit( INVALID_INPUT );
    }

    return TRUE;
}


/* Validate atom serial numbers in BGF geometry file
 *
 * val: atom serial to validate
 * lo: lower limit of valid serial range
 * hi: upper limit of valid serial range
 * filename: source filename of caller
 * line: source line of caller
 */
int Check_Input_Range( int val, int lo, int hi, const char * const filename,
        int line )
{
    if ( val < lo || val > hi )
    {
        fprintf( stderr, "[ERROR] Invalid BGF serial\n" );
        fprintf( stderr, "    [INFO] At line %d in file %.*s\n",
                line, (int) strlen(filename), filename );
        fprintf( stderr, "    [INFO] Input %d - Out of range %d-%d. Terminating...\n",
                 val, lo, hi );
        exit( INVALID_INPUT );
    }

    return SUCCESS;
}


void Trim_Spaces( char * const element, const size_t size )
{
    int i, j, n;

    n = strnlen( element, size );

    /* buffer not NULL-terminated, abort */
    if ( n == size )
    {
        fprintf( stderr, "[ERROR] buffer not NULL-terminated (Trim_Spaces). Terminating...\n" );
        exit( RUNTIME_ERROR );
    }

    /* skip initial space chars */
    for ( i = 0; element[i] == ' '; ++i )
        ;

    /* make uppercase, offset to 0 */
    for ( j = i; j < n && element[j] != ' '; ++j )
    {
        element[j - i] = toupper( element[j] );
    }

    /* NULL terminate */
    element[j - i] = '\0';
}


/************ from system_props.c *************/
real Get_Time( )
{
    int ret;
    struct timespec t;

    ret = clock_gettime( CLOCK_MONOTONIC, &t );

    if ( ret != 0 )
    {
        fprintf( stderr, "[WARNING] non-zero error in measuring time\n" );
    }

    return t.tv_sec + t.tv_nsec / 1.0e9;
}


real Get_Timing_Info( real t_start )
{
    int ret;
    struct timespec t_end;

    ret = clock_gettime( CLOCK_MONOTONIC, &t_end );

    if ( ret != 0 )
    {
        fprintf( stderr, "[WARNING] non-zero error in measuring time\n" );
    }

    return t_end.tv_sec + t_end.tv_nsec / 1.0e9 - t_start;
}


/*********** from io_tools.c **************/
int Get_Atom_Type( reax_interaction *reax_param, char *s, size_t n )
{
    int i, ret, flag;
    
    flag = FAILURE;

    for ( i = 0; i < reax_param->num_atom_types; ++i )
    {
        if ( strncmp( reax_param->sbp[i].name, s,
                    MIN( sizeof(reax_param->sbp[i].name), n ) ) == 0 )
        {
            ret = i;
            flag = SUCCESS;
            break;
        }
    }

    if ( flag == FAILURE )
    {
        fprintf( stderr, "[ERROR] Unknown atom type: %s. Terminating...\n", s );
        exit( UNKNOWN_ATOM_TYPE );
    }

    return ret;
}


char *Get_Element( reax_system *system, int i )
{
    return system->reax_param.sbp[system->atoms[i].type].name;
}


char *Get_Atom_Name( reax_system *system, int i )
{
    return system->atoms[i].name;
}


void Allocate_Tokenizer_Space( char **line, size_t line_size,
        char **backup, size_t backup_size,
        char ***tokens, size_t num_tokens, size_t token_size )
{
    int i;

    *line = smalloc( sizeof(char) * line_size, __FILE__, __LINE__ );
    *backup = smalloc( sizeof(char) * backup_size, __FILE__, __LINE__ );
    *tokens = smalloc( sizeof(char*) * num_tokens, __FILE__, __LINE__ );

    for ( i = 0; i < num_tokens; i++ )
    {
        (*tokens)[i] = smalloc( sizeof(char) * token_size,
                __FILE__, __LINE__ );
    }
}


void Deallocate_Tokenizer_Space( char **line, char **backup,
        char ***tokens, size_t num_tokens )
{
    int i;

    for ( i = 0; i < num_tokens; i++ )
    {
        sfree( (*tokens)[i], __FILE__, __LINE__ );
    }

    sfree( *line, __FILE__, __LINE__ );
    sfree( *backup, __FILE__, __LINE__ );
    sfree( *tokens, __FILE__, __LINE__ );
}


int Tokenize( char* s, char*** tok, size_t token_len )
{
    int count, word_len;
    char test[MAX_LINE];
    char *sep = "\t \n!=";
    char *word, *saveptr;

    count = 0;

    strncpy( test, s, sizeof(test) );
    test[sizeof(test) - 1] = '\0';

    for ( word = strtok_r(test, sep, &saveptr); word != NULL;
            word = strtok_r(NULL, sep, &saveptr) )
    {
        word_len = MIN( strlen(word), token_len - 1 );
        strncpy( (*tok)[count], word, word_len );
        (*tok)[count][word_len] = '\0';
        count++;
    }

    return count;
}


/* Safe wrapper around libc malloc
 *
 * n: num. of bytes to allocated
 * filename: source filename of caller
 * line: source line of caller
 *
 * returns: ptr to allocated memory
 * */
void * smalloc( size_t n, const char * const filename, int line )
{
    void *ptr;

    if ( n == 0 )
    {
        fprintf( stderr, "[ERROR] failed to allocate %zu bytes for array\n",
                n );
        fprintf( stderr, "    [INFO] At line %d in file %.*s\n",
                line, (int) strlen(filename), filename );
        exit( INSUFFICIENT_MEMORY );
    }

#if defined(DEBUG_FOCUS)
    fprintf( stderr, "[INFO] requesting allocation of %zu bytes of memory at line %d in file %.*s\n",
            n, line, (int) strlen(filename), filename );
    fflush( stderr );
#endif

    ptr = malloc( n );

    if ( ptr == NULL )
    {
        fprintf( stderr, "[ERROR] failed to allocate %zu bytes for array\n",
                n );
        fprintf( stderr, "    [INFO] At line %d in file %.*s\n",
                line, (int) strlen(filename), filename );
        exit( INSUFFICIENT_MEMORY );
    }

#if defined(DEBUG_FOCUS)
    fprintf( stderr, "[INFO] address: %p [SMALLOC]\n", (void *) ptr );
    fflush( stderr );
#endif

    return ptr;
}


/* Safe wrapper around libc realloc
 *
 * n: num. of bytes to reallocated
 * filename: source filename of caller
 * line: source line of caller
 *
 * returns: ptr to reallocated memory
 * */
void * srealloc( void *ptr, size_t n, const char * const filename, int line )
{
    void *new_ptr;

    if ( n == 0 )
    {
        fprintf( stderr, "[ERROR] failed to reallocate %zu bytes for array\n",
                n );
        fprintf( stderr, "    [INFO] At line %d in file %.*s\n",
                line, (int) strlen(filename), filename );
        exit( INSUFFICIENT_MEMORY );
    }

#if defined(DEBUG_FOCUS)
    fprintf( stderr, "[INFO] requesting reallocation of %zu bytes of memory at line %d in file %.*s\n",
            n, line, (int) strlen(filename), filename );
    fflush( stderr );
#endif

    new_ptr = realloc( ptr, n );

    /* technically, ptr may still be allocated and valid,
     * but we needed more memory, so abort */
    if ( new_ptr == NULL )
    {
        fprintf( stderr, "[ERROR] failed to reallocate %zu bytes for array\n",
                n );
        fprintf( stderr, "    [INFO] At line %d in file %.*s\n",
                line, (int) strlen(filename), filename );
        exit( INSUFFICIENT_MEMORY );
    }

#if defined(DEBUG_FOCUS)
    fprintf( stderr, "[INFO] address: %p [SREALLOC]\n", (void *) new_ptr );
    fflush( stderr );
#endif

    return new_ptr;
}


/* Safe wrapper around libc calloc
 *
 * n: num. of elements to allocated (each of size bytes)
 * size: num. of bytes per element
 * filename: source filename of caller
 * line: source line of caller
 *
 * returns: ptr to allocated memory, all bits initialized to zeros
 * */
void * scalloc( size_t n, size_t size, const char * const filename, int line )
{
    void *ptr;

    if ( n == 0 )
    {
        fprintf( stderr, "[ERROR] failed to allocate %zu bytes for array\n",
                n * size );
        fprintf( stderr, "    [INFO] At line %d in file %.*s\n",
                line, (int) strlen(filename), filename );
        exit( INSUFFICIENT_MEMORY );
    }

#if defined(DEBUG_FOCUS)
    fprintf( stderr, "[INFO] requesting allocation of %zu bytes of zeroed memory at line %d in file %.*s\n",
            n * size, line, (int) strlen(filename), filename );
    fflush( stderr );
#endif

    ptr = calloc( n, size );

    if ( ptr == NULL )
    {
        fprintf( stderr, "[ERROR] failed to allocate %zu bytes for array\n",
                n * size );
        fprintf( stderr, "    [INFO] At line %d in file %.*s\n",
                line, (int) strlen(filename), filename );
        exit( INSUFFICIENT_MEMORY );
    }

#if defined(DEBUG_FOCUS)
    fprintf( stderr, "[INFO] address: %p [SCALLOC]\n", (void *) ptr );
    fflush( stderr );
#endif

    return ptr;
}


/* Safe wrapper around libc free
 *
 * ptr: pointer to dynamically allocated memory which will be deallocated
 * filename: source filename of caller
 * line: source line of caller
 * */
void sfree( void *ptr, const char * const filename, int line )
{
    if ( ptr == NULL )
    {
        fprintf( stderr, "[WARNING] trying to free the already NULL pointer\n" );
        fprintf( stderr, "    [INFO] At line %d in file %.*s\n",
                line, (int) strlen(filename), filename );
        return;
    }

#if defined(DEBUG_FOCUS)
    fprintf( stderr, "[INFO] trying to free pointer at line %d in file %.*s\n",
            line, (int) strlen(filename), filename );
    fflush( stderr );
    fprintf( stderr, "[INFO] address: %p [SFREE]\n", (void *) ptr );
    fflush( stderr );
#endif

    free( ptr );
}


/* Safe wrapper around libc fopen
 *
 * fname: name of file to be opened
 * mode: mode in which to open file
 * filename: source filename of caller
 * line: source line of caller
 * */
FILE * sfopen( const char * fname, const char * mode,
        const char * const filename, int line )
{
    FILE * ptr;

    if ( fname == NULL )
    {
        fprintf( stderr, "[ERROR] trying to open file\n" );
        fprintf( stderr, "    [INFO] At line %d in file %.*s\n",
                line, (int) strlen(filename), filename );
        fprintf( stderr, "  [INFO] NULL file name\n" );
        exit( INVALID_INPUT );
    }
    if ( mode == NULL )
    {
        fprintf( stderr, "[ERROR] trying to open file\n" );
        fprintf( stderr, "    [INFO] At line %d in file %.*s\n",
                line, (int) strlen(filename), filename );
        fprintf( stderr, "  [INFO] NULL mode\n" );
        exit( INVALID_INPUT );
    }

    ptr = fopen( fname, mode );

    if ( ptr == NULL )
    {
        fprintf( stderr, "[ERROR] failed to open file %s with mode %s\n",
              fname, mode );
        fprintf( stderr, "    [INFO] At line %d in file %.*s\n",
                line, (int) strlen(filename), filename );
        exit( INVALID_INPUT );
    }

    return ptr;
}


/* Safe wrapper around libc fclose
 *
 * fp: pointer to file to close
 * filename: source filename of caller
 * line: source line of caller
 * */
void sfclose( FILE * fp, const char * const filename, int line )
{
    int ret;

    if ( fp == NULL )
    {
        fprintf( stderr, "[WARNING] trying to close NULL file pointer. Returning...\n" );
        fprintf( stderr, "    [INFO] At line %d in file %.*s\n",
                line, (int) strlen(filename), filename );
        return;
    }

    ret = fclose( fp );

    if ( ret != 0 )
    {
        fprintf( stderr, "[ERROR] error detected when closing file\n" );
        fprintf( stderr, "    [INFO] At line %d in file %.*s\n",
                line, (int) strlen(filename), filename );
        exit( INVALID_INPUT );
    }
}


/* Safe wrapper around strtol
 *
 * str: string to be converted
 * filename: source filename of caller
 * line: source line of caller
 *
 * returns: result of conversion (integer)
 * */
int sstrtol( const char * const str,
        const char * const filename, int line )
{
    long ret;
    char *endptr;

    if ( str[0] == '\0' )
    {
        fprintf( stderr, "[ERROR] sstrtol: NULL string\n" );
        /* strlen safe here only if filename is NULL-terminated
         * before calling sconvert_string_to_int */
        fprintf( stderr, "    [INFO] At line %d in file %.*s\n",
                line, (int) strlen(filename), filename );
        exit( INVALID_INPUT );
    }

    errno = 0;
    ret = strtol( str, &endptr, INTBASE );

    if ( (errno == ERANGE && (ret == LONG_MAX || ret == LONG_MIN) )
            || (errno != 0 && ret == 0) )
    {
        fprintf( stderr, "[ERROR] strtol: invalid string\n" );
        /* strlen safe here only if filename is NULL-terminated
         * before calling sconvert_string_to_int */
        fprintf( stderr, "    [INFO] At line %d in file %.*s\n",
                line, (int) strlen(filename), filename );
        fprintf( stderr, "    [INFO] str: %.*s\n",
                (int) strlen(str), str );
        exit( INVALID_INPUT );
    }
    else if ( endptr == str )
    {
        fprintf( stderr, "[ERROR] strtol: no digits found\n" );
        /* strlen safe here only if filename is NULL-terminated
         * before calling sconvert_string_to_int */
        fprintf( stderr, "    [INFO] At line %d in file %.*s\n",
                line, (int) strlen(filename), filename );
        fprintf( stderr, "    [INFO] str: %.*s\n",
                (int) strlen(str), str );
        exit( INVALID_INPUT );
    }
    else if ( *endptr != '\0' )
    {
        fprintf( stderr, "[ERROR] strtol: non-numeric trailing characters\n" );
        /* strlen safe here only if filename is NULL-terminated
         * before calling sconvert_string_to_int */
        fprintf( stderr, "    [INFO] At line %d in file %.*s\n",
                line, (int) strlen(filename), filename );
        fprintf( stderr, "    [INFO] str: %.*s\n",
                (int) strlen(str), str );
        exit( INVALID_INPUT );
    }

    return (int) ret;
}


/* Safe wrapper around strtod
 *
 * str: string to be converted
 * filename: source filename of caller
 * line: source line of caller
 *
 * returns: result of conversion (double)
 * */
double sstrtod( const char * const str,
        const char * const filename, int line )
{
    double ret;
    char *endptr;

    if ( str[0] == '\0' )
    {
        fprintf( stderr, "[ERROR] sstrtod: NULL string\n" );
        /* strlen safe here only if filename is NULL-terminated
         * before calling sconvert_string_to_int */
        fprintf( stderr, "    [INFO] At line %d in file %.*s\n",
                line, (int) strlen(filename), filename );
        exit( INVALID_INPUT );
    }

    errno = 0;
    ret = strtod( str, &endptr );

    if ( (errno == ERANGE && (ret == LONG_MAX || ret == LONG_MIN) )
            || (errno != 0 && ret == 0) )
    {
        fprintf( stderr, "[ERROR] strtod: invalid string\n" );
        /* strlen safe here only if filename is NULL-terminated
         * before calling sconvert_string_to_int */
        fprintf( stderr, "    [INFO] At line %d in file %.*s\n",
                line, (int) strlen(filename), filename );
        fprintf( stderr, "    [INFO] str: %.*s\n",
                (int) strlen(str), str );
        exit( INVALID_INPUT );
    }
    else if ( endptr == str )
    {
        fprintf( stderr, "[ERROR] strtod: no digits found\n" );
        /* strlen safe here only if filename is NULL-terminated
         * before calling sconvert_string_to_int */
        fprintf( stderr, "    [INFO] At line %d in file %.*s\n",
                line, (int) strlen(filename), filename );
        fprintf( stderr, "    [INFO] str: %.*s\n",
                (int) strlen(str), str );
        exit( INVALID_INPUT );
    }
    else if ( *endptr != '\0' )
    {
        fprintf( stderr, "[ERROR] strtod: non-numeric trailing characters\n" );
        /* strlen safe here only if filename is NULL-terminated
         * before calling sconvert_string_to_int */
        fprintf( stderr, "    [INFO] At line %d in file %.*s\n",
                line, (int) strlen(filename), filename );
        fprintf( stderr, "    [INFO] str: %.*s\n",
                (int) strlen(str), str );
        exit( INVALID_INPUT );
    }

    return ret;
}
