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

#include "neighbors.h"

#include "box.h"
#include "grid.h"
#if defined(DEBUG_FOCUS)
  #include "io_tools.h"
#endif
#include "list.h"
#include "reset_tools.h"
#include "system_props.h"
#include "tool_box.h"
#include "vector.h"


/* If the code is not compiled to handle small periodic boxes (i.e. a
 * simulation box with any dimension less than twice the Verlet list cutoff
 * distance, vlist_cut), it will use the optimized Generate_Neighbor_Lists
 * function.  Otherwise it will execute the neighbor routine with small
 * periodic box support.
 * 
 * Define the preprocessor definition SMALL_BOX_SUPPORT to enable (in
 * reax_types.h). */
typedef int (*count_far_neighbors_function)( rvec, rvec, int, int,
        simulation_box const * const, real );

typedef int (*find_far_neighbors_function)( rvec, rvec, int, int,
        simulation_box const * const, real, far_neighbor_data * const, int );


static void Choose_Neighbor_Counter( reax_system const * const system,
        control_params const * const control,
        count_far_neighbors_function *Count_Far_Neighbors )
{
    if ( control->periodic_boundaries == TRUE )
    {
#if defined(SMALL_BOX_SUPPORT)
        if ( system->box.box_norms[0] >= 2.0 * control->vlist_cut
                && system->box.box_norms[1] >= 2.0 * control->vlist_cut
                && system->box.box_norms[2] >= 2.0 * control->vlist_cut )
        {
            *Count_Far_Neighbors = &Count_Periodic_Far_Neighbors_Big_Box;
        }
        else
        {
            *Count_Far_Neighbors = &Count_Periodic_Far_Neighbors_Small_Box;
        }
#else
        *Count_Far_Neighbors = &Count_Periodic_Far_Neighbors_Big_Box;
#endif
    }
    else
    {
        *Count_Far_Neighbors = &Count_Non_Periodic_Far_Neighbors;
    }
}


static void Choose_Neighbor_Finder( reax_system const * const system,
        control_params const * const control,
        find_far_neighbors_function *Find_Far_Neighbors )
{
    if ( control->periodic_boundaries == TRUE )
    {
#if defined(SMALL_BOX_SUPPORT)
        if ( system->box.box_norms[0] >= 2.0 * control->vlist_cut
                && system->box.box_norms[1] >= 2.0 * control->vlist_cut
                && system->box.box_norms[2] >= 2.0 * control->vlist_cut )
        {
            *Find_Far_Neighbors = &Find_Periodic_Far_Neighbors_Big_Box;
        }
        else
        {
            *Find_Far_Neighbors = &Find_Periodic_Far_Neighbors_Small_Box;
        }
#else
        *Find_Far_Neighbors = &Find_Periodic_Far_Neighbors_Big_Box;
#endif
    }
    else
    {
        *Find_Far_Neighbors = &Find_Non_Periodic_Far_Neighbors;
    }
}


#if defined(DEBUG_FOCUS)
static int compare_far_nbrs( const void *v1, const void *v2 )
{
    return ((*(far_neighbor_data *)v1).nbr - (*(far_neighbor_data *)v2).nbr);
}
#endif


static inline real DistSqr_to_CP( rvec cp, rvec x )
{
    int i;
    real d_sqr;

    d_sqr = 0.0;

    for ( i = 0; i < 3; ++i )
    {
        if ( cp[i] > NEG_INF )
        {
            d_sqr += SQR( cp[i] - x[i] );
        }
    }

    return d_sqr;
}


/* Estimate the storage requirements for the far neighbor list */
void Estimate_Num_Neighbors( reax_system const * const system,
        control_params const * const control, static_storage * const workspace,
        reax_list ** const lists )
{
    int i, j, k, l, m, itr;
    int x, y, z;
    int atom1, atom2, max;
    int num_far, count;
    int *nbr_atoms;
    ivec *nbrs;
    rvec *nbrs_cp;
    grid const * const g = &system->g;
    count_far_neighbors_function Count_Far_Neighbors;

    num_far = 0;

    Choose_Neighbor_Counter( system, control, &Count_Far_Neighbors );

    /* for each cell in the grid along the 3
     * Cartesian directions: (i, j, k) => (x, y, z) */
    for ( i = 0; i < g->ncell[0]; i++ )
    {
        for ( j = 0; j < g->ncell[1]; j++ )
        {
            for ( k = 0; k < g->ncell[2]; k++ )
            {
                nbrs = g->nbrs[i][j][k];
                nbrs_cp = g->nbrs_cp[i][j][k];

                /* for each atom in the current cell */
                for ( l = 0; l < g->top[i][j][k]; ++l )
                {
                    atom1 = g->atoms[i][j][k][l];
                    itr = 0;

                    /* for each of the neighboring grid cells within
                     * the Verlet list cutoff distance */
                    while ( nbrs[itr][0] >= 0 )
                    {
                        /* if the Verlet list cutoff covers the closest point
                         * in the neighboring grid cell, then search through the cell's atoms */
                        if ( DistSqr_to_CP(nbrs_cp[itr], system->atoms[atom1].x )
                                <= SQR(control->vlist_cut) )
                        {
                            x = nbrs[itr][0];
                            y = nbrs[itr][1];
                            z = nbrs[itr][2];
                            nbr_atoms = g->atoms[x][y][z];
                            max = g->top[x][y][z];

                            /* pick up another atom from the neighbor cell;
                             * we have to compare atom1 with its own periodic images as well
                             * in the case of periodic boundary conditions,
                             * hence the equality in the if stmt below */
                            for ( m = 0; m < max; ++m )
                            {
                                atom2 = nbr_atoms[m];

                                if ( atom1 >= atom2 )
                                {
                                    count = Count_Far_Neighbors( system->atoms[atom1].x,
                                                system->atoms[atom2].x, atom1, atom2, 
                                                &system->box, control->vlist_cut );

                                    num_far += count;
                                }
                            }
                        }

                        ++itr;
                    }
                }
            }
        }
    }

    workspace->realloc.total_far_nbrs = (int) CEIL( num_far * SAFE_ZONE );
}


/* Generate the far neighbor list */
int Generate_Neighbor_Lists( reax_system * const system,
        control_params const * const control, simulation_data * const data,
        static_storage * const workspace, reax_list ** const lists )
{
    int i, j, k, l, m, itr;
    int x, y, z;
    int atom1, atom2, max;
    int num_far, count, flag_oom, ret;
    int *nbr_atoms;
    ivec *nbrs;
    rvec *nbrs_cp;
    grid const * const g = &system->g;
    reax_list *far_nbrs;
    find_far_neighbors_function Find_Far_Neighbors;
    real t_start, t_elapsed;

    t_start = Get_Time( );
    num_far = 0;
    far_nbrs = lists[FAR_NBRS];
    flag_oom = FALSE;
    ret = SUCCESS;

    Choose_Neighbor_Finder( system, control, &Find_Far_Neighbors );

    for ( i = 0; i < far_nbrs->n; ++i )
    {
        Set_Start_Index( i, 0, far_nbrs );
        Set_End_Index( i, 0, far_nbrs );
    }

    /* for each cell in the grid along the 3
     * Cartesian directions: (i, j, k) => (x, y, z) */
    for ( i = 0; i < g->ncell[0]; i++ )
    {
        for ( j = 0; j < g->ncell[1]; j++ )
        {
            for ( k = 0; k < g->ncell[2]; k++ )
            {
                nbrs = g->nbrs[i][j][k];
                nbrs_cp = g->nbrs_cp[i][j][k];

                /* for each atom in the current cell */
                for ( l = 0; l < g->top[i][j][k]; ++l )
                {
                    atom1 = g->atoms[i][j][k][l];
                    Set_Start_Index( atom1, num_far, far_nbrs );
                    itr = 0;

                    /* for each of the neighboring grid cells within
                     * the Verlet list cutoff distance */
                    while ( nbrs[itr][0] >= 0 )
                    {
                        /* if the Verlet list cutoff covers the closest point
                         * in the neighboring grid cell, then search through the cell's atoms */
                        if ( DistSqr_to_CP( nbrs_cp[itr], system->atoms[atom1].x )
                                <= SQR(control->vlist_cut) )
                        {
                            x = nbrs[itr][0];
                            y = nbrs[itr][1];
                            z = nbrs[itr][2];
                            nbr_atoms = g->atoms[x][y][z];
                            max = g->top[x][y][z];

                            /* pick up another atom from the neighbor cell;
                             * we have to compare atom1 with its own periodic images as well
                             * in the case of periodic boundary conditions,
                             * hence the equality in the if stmt below */
                            for ( m = 0; m < max; ++m )
                            {
                                atom2 = nbr_atoms[m];

                                if ( atom1 >= atom2 )
                                {
                                    count = Find_Far_Neighbors( system->atoms[atom1].x,
                                            system->atoms[atom2].x, atom1, atom2,
                                            &system->box, control->vlist_cut,
                                            &far_nbrs->far_nbr_list[num_far],
                                            far_nbrs->total_intrs - num_far );

                                    num_far += count;

                                    if ( num_far >= far_nbrs->total_intrs )
                                    {
                                        goto OUT_OF_MEMORY;
                                    }
                                }
                            }
                        }

                        ++itr;
                    }

                    Set_End_Index( atom1, num_far, far_nbrs );
                }
            }
        }
    }

    //TODO: conditionally perform these assignments if periodic boundary conditions are enabled
    for ( i = 0; i < system->N; i++ )
    {
        ivec_MakeZero( system->atoms[i].rel_map );
    }

    if ( flag_oom == TRUE )
    {
        OUT_OF_MEMORY:
            ret = FAILURE;
    }

#if defined(DEBUG_FOCUS)
    for ( i = 0; i < system->N; ++i )
    {
        qsort( &far_nbrs->far_nbr_list[ Start_Index(i, far_nbrs) ],
                Num_Entries(i, far_nbrs), sizeof(far_neighbor_data),
                compare_far_nbrs );
    }
#endif

#if defined(DEBUG_FOCUS)
    Print_Far_Neighbors( system, control, data, workspace, lists );
#endif

    t_elapsed = Get_Timing_Info( t_start );
    data->timing.nbrs += t_elapsed;

    return ret;
}
