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

#include "allocate.h"

#include "grid.h"
#include "list.h"
#include "tool_box.h"


/* allocate space for atoms */
void PreAllocate_Space( reax_system * const system,
        control_params const  * const control,
        static_storage * const workspace, int n )
{
    int i;

    if ( system->prealloc_allocated == FALSE )
    {
        system->prealloc_allocated = TRUE;

        system->atoms = scalloc( n, sizeof(reax_atom), __FILE__, __LINE__ );
        workspace->orig_id = scalloc( n, sizeof(int), __FILE__, __LINE__ );

        /* bond restriction info */
        if ( control->restrict_bonds )
        {
            workspace->restricted = scalloc( n, sizeof(int), __FILE__, __LINE__ );
            workspace->restricted_list = scalloc( n, sizeof(int*), __FILE__, __LINE__ );

            for ( i = 0; i < n; ++i )
            {
                workspace->restricted_list[i] = scalloc( MAX_RESTRICT, sizeof(int),
                        __FILE__, __LINE__ );
            }
        }

        if ( control->geo_format == BGF
                || control->geo_format == ASCII_RESTART
                || control->geo_format == BINARY_RESTART )
        {
            workspace->map_serials = scalloc( MAX_ATOM_ID, sizeof(int),
                    __FILE__, __LINE__ );
        }
    }
    else
    {
        sfree( system->atoms, __FILE__, __LINE__ );
        sfree( workspace->orig_id, __FILE__, __LINE__ );

        /* bond restriction info */
        if ( control->restrict_bonds )
        {
            sfree( workspace->restricted, __FILE__, __LINE__ );

            for ( i = 0; i < n; ++i )
            {
                sfree( workspace->restricted_list[i], __FILE__, __LINE__ );
            }

            sfree( workspace->restricted_list, __FILE__, __LINE__ );
        }

        system->atoms = scalloc( n, sizeof(reax_atom), __FILE__, __LINE__ );
        workspace->orig_id = scalloc( n, sizeof(int), __FILE__, __LINE__ );

        /* bond restriction info */
        if ( control->restrict_bonds )
        {
            workspace->restricted = scalloc( n, sizeof(int), __FILE__, __LINE__ );
            workspace->restricted_list = scalloc( n, sizeof(int*), __FILE__, __LINE__ );

            for ( i = 0; i < n; ++i )
            {
                workspace->restricted_list[i] = scalloc( MAX_RESTRICT, sizeof(int),
                        __FILE__, __LINE__ );
            }
        }
    }
}


static void Reallocate_Neighbor_List( reax_list * const far_nbr_list, int n,
        int n_max, int num_intrs )
{
    if ( far_nbr_list->allocated == TRUE )
    {
        Delete_List( TYP_FAR_NEIGHBOR, far_nbr_list );
    }

    Make_List( n, n_max, num_intrs, TYP_FAR_NEIGHBOR, far_nbr_list );
}


/* Dynamic allocation of memory for matrix in CSR format
 *
 * pH (output): pointer to sparse matrix for which to allocate
 * n: num. rows of the matrix
 * n_max: max. num. rows of the matrix
 * m: number of nonzeros to allocate space for in matrix
 * */
void Allocate_Matrix( sparse_matrix * const H, int n, int n_max, int m )
{
    H->allocated = TRUE;

    H->n = n;
    H->n_max = n_max;
    H->m = m;

    H->start = smalloc( sizeof(unsigned int) * (n_max + 1), __FILE__, __LINE__ );
    H->j = smalloc( sizeof(unsigned int) * m, __FILE__, __LINE__ );
    H->val = smalloc( sizeof(real) * m, __FILE__, __LINE__ );
}


/* Deallocate memory for matrix in CSR format
 *
 * H (output): pointer to sparse matrix for which to allocate
 * */
void Deallocate_Matrix( sparse_matrix * const H )
{
    H->allocated = FALSE;

    sfree( H->start, __FILE__, __LINE__ );
    sfree( H->j, __FILE__, __LINE__ );
    sfree( H->val, __FILE__, __LINE__ );
}


static void Reallocate_Matrix( sparse_matrix *H, int n, int n_max, int m )
{
    Deallocate_Matrix( H );
    Allocate_Matrix( H, n, n_max, m );
}


static void Reallocate_List( reax_list * const list, int n, int n_max,
        int max_intrs, int type )
{
    if ( list->allocated == TRUE )
    {
        Delete_List( type, list );
    }
    Make_List( n, n_max, max_intrs, type, list );
}


void Reallocate_Part1( reax_system * const system, control_params const * const control,
        static_storage * const workspace, reax_list ** const lists )
{
    int i, j, k;
    reallocate_data *realloc;
    grid *g;

    realloc = &workspace->realloc;
    g = &system->g;

    if ( realloc->gcell_atoms > -1 )
    {
#if defined(DEBUG_FOCUS)
        fprintf( stderr, "[INFO] reallocating gcell: g->max_atoms: %d\n", g->max_atoms );
#endif

        for ( i = 0; i < g->ncell_max[0]; i++ )
        {
            for ( j = 0; j < g->ncell_max[1]; j++ )
            {
                for ( k = 0; k < g->ncell_max[2]; k++ )
                {
                    sfree( g->atoms[i][j][k], __FILE__, __LINE__ );
                    g->atoms[i][j][k] = scalloc( workspace->realloc.gcell_atoms,
                            sizeof(int), __FILE__, __LINE__ );
                }
            }
        }

        realloc->gcell_atoms = -1;
    }
}


void Reallocate_Part2( reax_system const * const system,
        control_params const * const control, simulation_data const * const data,
        static_storage * const workspace, reax_list ** const lists )
{
    int renbr;
    reallocate_data *realloc;

    realloc = &workspace->realloc;
    renbr = ((data->step - data->prev_steps) % control->reneighbor) == 0 ? TRUE : FALSE;

    if ( renbr == TRUE && realloc->far_nbrs == TRUE )
    {
        Reallocate_Neighbor_List( lists[FAR_NBRS],
                system->N, system->N_max, realloc->total_far_nbrs );

        realloc->far_nbrs = FALSE;
    }

    if ( realloc->cm == TRUE )
    {
        Reallocate_Matrix( &workspace->H, system->N_cm, system->N_cm_max,
                realloc->total_cm_entries );
        Reallocate_Matrix( &workspace->H_sp, system->N_cm, system->N_cm_max,
                realloc->total_cm_entries );

        realloc->cm = FALSE;
    }

    if ( realloc->bonds == TRUE )
    {
        Reallocate_List( lists[BONDS], system->N, system->N_max,
                realloc->total_bonds, TYP_BOND );

        realloc->bonds = FALSE;
        realloc->thbody = TRUE;
    }

    if ( control->hbond_cut > 0.0 && workspace->num_H > 0 && realloc->hbonds == TRUE )
    {
        Reallocate_List( lists[HBONDS], system->N, system->N_max,
                realloc->total_hbonds, TYP_HBOND );

        realloc->hbonds = FALSE;
    }

    if ( realloc->thbody == TRUE )
    {
        Reallocate_List( lists[THREE_BODIES], realloc->total_bonds,
                (int) CEIL( realloc->total_bonds * SAFE_ZONE ),
                realloc->total_thbodies, TYP_THREE_BODY );

        realloc->thbody = FALSE;
    }
}
