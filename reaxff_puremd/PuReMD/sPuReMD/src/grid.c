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

#include "grid.h"

#include "reset_tools.h"
#include "tool_box.h"
#include "vector.h"


/* Counts the number of atoms binned to each cell within the grid
 *
 * NOTE: this assumes the atom positions are within the 
 * simulation box boundaries [0, d_i) where d_i is the simulation
 * box length along a particular dimension */
static int Estimate_GCell_Population( reax_system * const system )
{
    int i, j, k, l;
    int max_atoms;
    grid *g;

    g = &system->g;

    Reset_Grid( g );

    for ( l = 0; l < system->N; l++ )
    {
        assert( system->atoms[l].x[0] >= 0.0 && system->atoms[l].x[0] < system->box.box_norms[0] );
        assert( system->atoms[l].x[1] >= 0.0 && system->atoms[l].x[1] < system->box.box_norms[1] );
        assert( system->atoms[l].x[2] >= 0.0 && system->atoms[l].x[2] < system->box.box_norms[2] );

        i = (int) (system->atoms[l].x[0] * g->inv_len[0]);
        j = (int) (system->atoms[l].x[1] * g->inv_len[1]);
        k = (int) (system->atoms[l].x[2] * g->inv_len[2]);
        g->top[i][j][k]++;

//        fprintf( stderr, "\tatom%-6d (%8.3f%8.3f%8.3f) --> (%3d%3d%3d)\n",
//                l, system->atoms[l].x[0], system->atoms[l].x[1], system->atoms[l].x[2],
//                i, j, k );
    }

    max_atoms = 0;
    for ( i = 0; i < g->ncell[0]; i++ )
    {
        for ( j = 0; j < g->ncell[1]; j++ )
        {
            for ( k = 0; k < g->ncell[2]; k++ )
            {
                if ( max_atoms < g->top[i][j][k] )
                {
                    max_atoms = g->top[i][j][k];
                }
            }
        }
    }

    return MAX( max_atoms * SAFE_ZONE, MIN_GCELL_POPL );
}


static void Allocate_Space_for_Grid( reax_system * const system, int alloc )
{
    int i, j, k, l, max_atoms;
    grid *g;

    g = &system->g;
    g->allocated = TRUE;
    g->max_nbrs = (2 * g->spread[0] + 1)
        * (2 * g->spread[1] + 1) * (2 * g->spread[2] + 1) + 3;

    if ( alloc == TRUE )
    {
        /* allocate space for the new grid */
        g->atoms = scalloc( g->ncell_max[0], sizeof( int*** ),
                __FILE__, __LINE__ );

        for ( i = 0; i < g->ncell_max[0]; i++ )
            g->atoms[i] = scalloc( g->ncell_max[1], sizeof( int** ),
                    __FILE__, __LINE__ );

        for ( i = 0; i < g->ncell_max[0]; i++ )
            for ( j = 0; j < g->ncell_max[1]; j++ )
                g->atoms[i][j] = scalloc( g->ncell_max[2], sizeof( int* ),
                        __FILE__, __LINE__ );

        g->top = scalloc( g->ncell_max[0], sizeof( int** ),
                __FILE__, __LINE__ );

        for ( i = 0; i < g->ncell_max[0]; i++ )
            g->top[i] = scalloc( g->ncell_max[1], sizeof( int* ),
                    __FILE__, __LINE__ );

        for ( i = 0; i < g->ncell_max[0]; i++ )
            for ( j = 0; j < g->ncell_max[1]; j++ )
                g->top[i][j] = scalloc( g->ncell_max[2], sizeof( int ),
                        __FILE__, __LINE__ );

        g->mark = scalloc( g->ncell_max[0], sizeof( int** ),
                __FILE__, __LINE__ );

        for ( i = 0; i < g->ncell_max[0]; i++ )
            g->mark[i] = scalloc( g->ncell_max[1], sizeof( int* ),
                    __FILE__, __LINE__ );

        for ( i = 0; i < g->ncell_max[0]; i++ )
            for ( j = 0; j < g->ncell_max[1]; j++ )
                g->mark[i][j] = scalloc( g->ncell_max[2], sizeof( int ),
                        __FILE__, __LINE__ );

        g->start = scalloc( g->ncell_max[0], sizeof( int** ),
                __FILE__, __LINE__ );

        for ( i = 0; i < g->ncell_max[0]; i++ )
            g->start[i] = scalloc( g->ncell_max[1], sizeof( int* ),
                    __FILE__, __LINE__ );

        for ( i = 0; i < g->ncell_max[0]; i++ )
            for ( j = 0; j < g->ncell_max[1]; j++ )
                g->start[i][j] = scalloc( g->ncell_max[2], sizeof( int ),
                        __FILE__, __LINE__ );

        g->end = scalloc( g->ncell_max[0], sizeof( int** ),
                __FILE__, __LINE__ );

        for ( i = 0; i < g->ncell_max[0]; i++ )
            g->end[i] = scalloc( g->ncell_max[1], sizeof( int* ),
                    __FILE__, __LINE__ );

        for ( i = 0; i < g->ncell_max[0]; i++ )
            for ( j = 0; j < g->ncell_max[1]; j++ )
                g->end[i][j] = scalloc( g->ncell_max[2], sizeof( int ),
                        __FILE__, __LINE__ );

        g->nbrs = scalloc( g->ncell_max[0], sizeof( ivec*** ),
                __FILE__, __LINE__ );

        for ( i = 0; i < g->ncell_max[0]; i++ )
            g->nbrs[i] = scalloc( g->ncell_max[1], sizeof( ivec** ),
                    __FILE__, __LINE__ );

        for ( i = 0; i < g->ncell_max[0]; i++ )
            for ( j = 0; j < g->ncell_max[1]; j++ )
                g->nbrs[i][j] = scalloc( g->ncell_max[2], sizeof( ivec* ),
                        __FILE__, __LINE__ );

        for ( i = 0; i < g->ncell_max[0]; i++ )
            for ( j = 0; j < g->ncell_max[1]; j++ )
                for ( k = 0; k < g->ncell_max[2]; k++ )
                    g->nbrs[i][j][k] = smalloc( g->max_nbrs * sizeof( ivec ),
                           __FILE__, __LINE__ );

        g->nbrs_cp = scalloc( g->ncell_max[0], sizeof( rvec*** ),
                __FILE__, __LINE__ );

        for ( i = 0; i < g->ncell_max[0]; i++ )
            g->nbrs_cp[i] = scalloc( g->ncell_max[1], sizeof( rvec** ),
                    __FILE__, __LINE__ );

        for ( i = 0; i < g->ncell_max[0]; i++ )
            for ( j = 0; j < g->ncell_max[1]; j++ )
                g->nbrs_cp[i][j] = scalloc( g->ncell_max[2], sizeof( rvec* ),
                        __FILE__, __LINE__ );

        for ( i = 0; i < g->ncell_max[0]; i++ )
            for ( j = 0; j < g->ncell_max[1]; j++ )
                for ( k = 0; k < g->ncell_max[2]; k++ )
                    g->nbrs_cp[i][j][k] = smalloc( g->max_nbrs * sizeof( rvec ),
                           __FILE__, __LINE__ );
    }

    for ( i = 0; i < g->ncell[0]; i++ )
        for ( j = 0; j < g->ncell[1]; j++ )
            for ( k = 0; k < g->ncell[2]; k++ )
                g->top[i][j][k] = 0;

    for ( i = 0; i < g->ncell[0]; i++ )
        for ( j = 0; j < g->ncell[1]; j++ )
            for ( k = 0; k < g->ncell[2]; k++ )
                g->mark[i][j][k] = 0;

    for ( i = 0; i < g->ncell[0]; i++ )
        for ( j = 0; j < g->ncell[1]; j++ )
            for ( k = 0; k < g->ncell[2]; k++ )
                g->start[i][j][k] = 0;

    for ( i = 0; i < g->ncell[0]; i++ )
        for ( j = 0; j < g->ncell[1]; j++ )
            for ( k = 0; k < g->ncell[2]; k++ )
                g->end[i][j][k] = 0;

    for ( i = 0; i < g->ncell[0]; i++ )
        for ( j = 0; j < g->ncell[1]; j++ )
            for ( k = 0; k < g->ncell[2]; k++ )
                for ( l = 0; l < g->max_nbrs; ++l )
                {
                    g->nbrs[i][j][k][l][0] = -1;
                    g->nbrs[i][j][k][l][1] = -1;
                    g->nbrs[i][j][k][l][2] = -1;
                }

    for ( i = 0; i < g->ncell[0]; i++ )
        for ( j = 0; j < g->ncell[1]; j++ )
            for ( k = 0; k < g->ncell[2]; k++ )
                for ( l = 0; l < g->max_nbrs; ++l )
                {
                    g->nbrs_cp[i][j][k][l][0] = -1.0;
                    g->nbrs_cp[i][j][k][l][1] = -1.0;
                    g->nbrs_cp[i][j][k][l][2] = -1.0;
                }

    max_atoms = Estimate_GCell_Population( system );

    if ( alloc == TRUE )
    {
        g->max_atoms = MAX( max_atoms, g->max_atoms );

        for ( i = 0; i < g->ncell_max[0]; i++ )
            for ( j = 0; j < g->ncell_max[1]; j++ )
                for ( k = 0; k < g->ncell_max[2]; k++ )
                    g->atoms[i][j][k] = smalloc( g->max_atoms * sizeof( int ),
                           __FILE__, __LINE__ );
    }
    /* case: grid large enough but max. atoms per grid cells insufficient */
    else if ( g->max_atoms > 0 && g->max_atoms < max_atoms )
    {
        g->max_atoms = max_atoms;

        for ( i = 0; i < g->ncell_max[0]; i++ )
            for ( j = 0; j < g->ncell_max[1]; j++ )
                for ( k = 0; k < g->ncell_max[2]; k++ )
                    sfree( g->atoms[i][j][k], __FILE__, __LINE__ );

        for ( i = 0; i < g->ncell_max[0]; i++ )
            for ( j = 0; j < g->ncell_max[1]; j++ )
                for ( k = 0; k < g->ncell_max[2]; k++ )
                    g->atoms[i][j][k] = smalloc( g->max_atoms * sizeof( int ),
                           __FILE__, __LINE__ );
    }

    for ( i = 0; i < g->ncell[0]; i++ )
        for ( j = 0; j < g->ncell[1]; j++ )
            for ( k = 0; k < g->ncell[2]; k++ )
                for ( l = 0; l < g->max_atoms; ++l )
                    g->atoms[i][j][k][l] = -1;
}


static void Deallocate_Grid_Space( grid * const g )
{
    int i, j, k;

    g->allocated = FALSE;

    /* deallocate the old grid */
    for ( i = 0; i < g->ncell_max[0]; i++ )
        for ( j = 0; j < g->ncell_max[1]; j++ )
            for ( k = 0; k < g->ncell_max[2]; k++ )
            {
                sfree( g->atoms[i][j][k], __FILE__, __LINE__ );
                sfree( g->nbrs[i][j][k], __FILE__, __LINE__ );
                sfree( g->nbrs_cp[i][j][k], __FILE__, __LINE__ );
            }

    for ( i = 0; i < g->ncell_max[0]; i++ )
        for ( j = 0; j < g->ncell_max[1]; j++ )
        {
            sfree( g->atoms[i][j], __FILE__, __LINE__ );
            sfree( g->top[i][j], __FILE__, __LINE__ );
            sfree( g->mark[i][j], __FILE__, __LINE__ );
            sfree( g->start[i][j], __FILE__, __LINE__ );
            sfree( g->end[i][j], __FILE__, __LINE__ );
            sfree( g->nbrs[i][j], __FILE__, __LINE__ );
            sfree( g->nbrs_cp[i][j], __FILE__, __LINE__ );
        }

    for ( i = 0; i < g->ncell_max[0]; i++ )
    {
        sfree( g->atoms[i], __FILE__, __LINE__ );
        sfree( g->top[i], __FILE__, __LINE__ );
        sfree( g->mark[i], __FILE__, __LINE__ );
        sfree( g->start[i], __FILE__, __LINE__ );
        sfree( g->end[i], __FILE__, __LINE__ );
        sfree( g->nbrs[i], __FILE__, __LINE__ );
        sfree( g->nbrs_cp[i], __FILE__, __LINE__ );
    }

    sfree( g->atoms, __FILE__, __LINE__ );
    sfree( g->top, __FILE__, __LINE__ );
    sfree( g->mark, __FILE__, __LINE__ );
    sfree( g->start, __FILE__, __LINE__ );
    sfree( g->end, __FILE__, __LINE__ );
    sfree( g->nbrs, __FILE__, __LINE__ );
    sfree( g->nbrs_cp, __FILE__, __LINE__ );
}


static inline int Shift( int p, int dp, int dim, grid const * const g )
{
    int dim_len, newp;

    dim_len = 0;
    newp = p + dp;

    switch ( dim )
    {
    case 0:
        dim_len = g->ncell[0];
        break;
    case 1:
        dim_len = g->ncell[1];
        break;
    case 2:
        dim_len = g->ncell[2];
    }

    while ( newp < 0 )
    {
        newp = newp + dim_len;
    }
    while ( newp >= dim_len )
    {
        newp = newp - dim_len;
    }

    return newp;
}


/* finds the closest point between two grid cells denoted by c1 and c2.
 * periodic boundary conditions are taken into consideration as well. */
static void Find_Closest_Point( grid const * const g, int c1x, int c1y, int c1z,
        int c2x, int c2y, int c2z, rvec closest_point )
{
    int i, d;
    ivec c1 = { c1x, c1y, c1z };
    ivec c2 = { c2x, c2y, c2z };

    for ( i = 0; i < 3; i++ )
    {
        if ( g->ncell[i] < 5 )
        {
            closest_point[i] = NEG_INF - 1.0;
            continue;
        }

        d = c2[i] - c1[i];
        if ( ABS(d) <= g->ncell[i] / 2 )
        {
            if ( d > 0 )
            {
                closest_point[i] = c2[i] * g->len[i];
            }
            else if ( d == 0 )
            {
                closest_point[i] = NEG_INF - 1.0;
            }
            else
            {
                closest_point[i] = (c2[i] + 1) * g->len[i];
            }
        }
        else
        {
            if ( d > 0 )
            {
                closest_point[i] = (c2[i] - g->ncell[i] + 1) * g->len[i];
            }
            else
            {
                closest_point[i] = (c2[i] + g->ncell[i]) * g->len[i];
            }
        }
    }
}


static void Find_Neighbor_Grid_Cells( grid * const g )
{
    int i, j, k;
    int di, dj, dk;
    int x, y, z;
    int stack_top;
    ivec *nbrs_stack;
    rvec *cp_stack;

    /* for each cell in the grid */
    for ( i = 0; i < g->ncell[0]; i++ )
    {
        for ( j = 0; j < g->ncell[1]; j++ )
        {
            for ( k = 0; k < g->ncell[2]; k++ )
            {
                nbrs_stack = g->nbrs[i][j][k];
                cp_stack = g->nbrs_cp[i][j][k];
                stack_top = 0;

                /* choose an unmarked neighbor cell */
                for ( di = -1 * g->spread[0]; di <= g->spread[0]; di++ )
                {
                    x = Shift( i, di, 0, g );

                    for ( dj = -1 * g->spread[1]; dj <= g->spread[1]; dj++ )
                    {
                        y = Shift( j, dj, 1, g );

                        for ( dk = -1 * g->spread[2]; dk <= g->spread[2]; dk++ )
                        {
                            z = Shift( k, dk, 2, g );

                            if ( !g->mark[x][y][z] )
                            {
                                /*(di < 0 || // 9 combinations
                                 (di == 0 && dj < 0) || // 3 combinations
                                 (di == 0 && dj == 0 && dk < 0) ) )*/

                                /* put the neighbor cell into the stack and mark it */
                                nbrs_stack[stack_top][0] = x;
                                nbrs_stack[stack_top][1] = y;
                                nbrs_stack[stack_top][2] = z;
                                g->mark[x][y][z] = 1;

                                Find_Closest_Point( g, i, j, k, x, y, z, cp_stack[stack_top] );
                                stack_top++;
                            }
                        }
                    }
                }

//                nbrs_stack[stack_top][0] = i;
//                nbrs_stack[stack_top][1] = j;
//                nbrs_stack[stack_top][2] = k;
//
//                Find_Closest_Point( g, i, j, k, i, j, k, cp_stack[stack_top] );
//
//                nbrs_stack[stack_top + 1][0] = -1;
//                nbrs_stack[stack_top + 1][1] = -1;
//                nbrs_stack[stack_top + 1][2] = -1;
//
//                Reset_Marks( g, nbrs_stack, stack_top + 1 );

                nbrs_stack[stack_top][0] = -1;
                nbrs_stack[stack_top][1] = -1;
                nbrs_stack[stack_top][2] = -1;

                Reset_Marks( g, nbrs_stack, stack_top );
            }
        }
    }
}


void Setup_Grid( reax_system * const system )
{
    int d, alloc;
    grid *g;
    simulation_box *my_box;

    alloc = FALSE;
    g = &system->g;
    my_box = &system->box;

    /* determine number of grid cells in each direction */
    ivec_rScale( g->ncell, 1.0 / g->cell_size, my_box->box_norms );

    for ( d = 0; d < 3; ++d )
    {
        if ( g->ncell[d] <= 0 )
        {
            g->ncell[d] = 1;
        }
    }

    if ( g->allocated == FALSE )
    {
        g->ncell_max[0] = (int) CEIL( SAFE_ZONE * g->ncell[0] );
        g->ncell_max[1] = (int) CEIL( SAFE_ZONE * g->ncell[1] );
        g->ncell_max[2] = (int) CEIL( SAFE_ZONE * g->ncell[2] );
        g->max_atoms = 0;

        alloc = TRUE;
    }
    else if ( g->ncell[0] > g->ncell_max[0] || g->ncell[1] > g->ncell_max[1]
            || g->ncell[2] > g->ncell_max[2] )
    {
        if ( g->allocated == TRUE )
        {
            Deallocate_Grid_Space( g );
        }

        if ( g->ncell[0] > g->ncell_max[0] )
        {
            g->ncell_max[0] = (int) CEIL( SAFE_ZONE * MAX( g->ncell[0], g->ncell_max[0] ) );
        }
        if ( g->ncell[1] > g->ncell_max[1] )
        {
            g->ncell_max[1] = (int) CEIL( SAFE_ZONE * MAX( g->ncell[1], g->ncell_max[1] ) );
        }
        if ( g->ncell[2] > g->ncell_max[2] )
        {
            g->ncell_max[2] = (int) CEIL( SAFE_ZONE * MAX( g->ncell[2], g->ncell_max[2] ) );
        }

        alloc = TRUE;
    }

    /* compute cell lengths */
    rvec_iDivide( g->len, my_box->box_norms, g->ncell );
    rvec_Invert( g->inv_len, g->len );

    Allocate_Space_for_Grid( system, alloc );
    Find_Neighbor_Grid_Cells( g );

#if defined(DEBUG_FOCUS)
    fprintf( stderr, "[INFO] g->ncell = (%d, %d, %d)\n", g->ncell[0], g->ncell[1], g->ncell[2] );
    fprintf( stderr, "[INFO] g->ncell_max = (%d, %d, %d)\n", g->ncell_max[0], g->ncell_max[1], g->ncell_max[2] );
    fprintf( stderr, "[INFO] g->len = (%5.2f, %5.2f, %5.2f)\n", g->len[0], g->len[1], g->len[2] );
    fprintf( stderr, "[INFO] g->max_atoms = %d\n", g->max_atoms );
#endif
}


void Update_Grid( reax_system * const system )
{
    int d, i, j, k, x, y, z, itr;
    ivec ncell;
    ivec *nbrs;
    rvec *nbrs_cp;
    grid *g;
    simulation_box *my_box;

    g = &system->g;
    my_box = &system->box;

    /* determine number of grid cells in each direction */
    ivec_rScale( ncell, 1.0 / g->cell_size, my_box->box_norms );

    for ( d = 0; d < 3; ++d )
    {
        if ( ncell[d] == 0 )
        {
            ncell[d] = 1;
        }
    }

    /* ncell are unchanged */
    if ( ivec_isEqual( ncell, g->ncell ) )
    {
        /* update cell lengths */
        rvec_iDivide( g->len, my_box->box_norms, g->ncell );
        rvec_Invert( g->inv_len, g->len );

        /* update closest point distances between gcells */
        for ( i = 0; i < g->ncell[0]; i++ )
        {
            for ( j = 0; j < g->ncell[1]; j++ )
            {
                for ( k = 0; k < g->ncell[2]; k++ )
                {
                    nbrs = g->nbrs[i][j][k];
                    nbrs_cp = g->nbrs_cp[i][j][k];

                    itr = 0;
                    while ( nbrs[itr][0] >= 0 )
                    {
                        x = nbrs[itr][0];
                        y = nbrs[itr][1];
                        z = nbrs[itr][2];

                        Find_Closest_Point( g, i, j, k, x, y, z, nbrs_cp[itr] );
                        ++itr;
                    }
                }
            }
        }
    }
    /* at least one of ncell has changed */
    else
    {
        if ( system->g.allocated == TRUE )
        {
            Deallocate_Grid_Space( g );
        }

        /* update number of grid cells */
        ivec_Copy( g->ncell, ncell );

        /* update cell lengths */
        rvec_iDivide( g->len, my_box->box_norms, g->ncell );
        rvec_Invert( g->inv_len, g->len );

        Allocate_Space_for_Grid( system, TRUE );
        Find_Neighbor_Grid_Cells( g );

#if defined(DEBUG_FOCUS)
        fprintf( stderr, "updated grid: " );
        fprintf( stderr, "ncell[%d %d %d] ",
                 g->ncell[0], g->ncell[1], g->ncell[2] );
        fprintf( stderr, "len[%5.2f %5.2f %5.2f] ",
                 g->len[0], g->len[1], g->len[2] );
        fprintf( stderr, "g->max_atoms = %d\n", g->max_atoms );
#endif
    }
}


/* Utilize a cell method approach to bin atoms to each cell within the grid
 *
 * NOTE: this assumes the atom positions are within the 
 * simulation box boundaries [0, d_i) where d_i is the simulation
 * box length along a particular dimension */
void Bin_Atoms( reax_system * const system, static_storage * const workspace )
{
    int i, j, k, l;
    int max_atoms;
    grid *g;

    g = &system->g;

    Reset_Grid( g );

    for ( l = 0; l < system->N; l++ )
    {
        assert( system->atoms[l].x[0] >= 0.0 && system->atoms[l].x[0] < system->box.box_norms[0] );
        assert( system->atoms[l].x[1] >= 0.0 && system->atoms[l].x[1] < system->box.box_norms[1] );
        assert( system->atoms[l].x[2] >= 0.0 && system->atoms[l].x[2] < system->box.box_norms[2] );

        i = (int) (system->atoms[l].x[0] * g->inv_len[0]);
        j = (int) (system->atoms[l].x[1] * g->inv_len[1]);
        k = (int) (system->atoms[l].x[2] * g->inv_len[2]);
        /* atom index in grid cell => atom number */
        g->atoms[i][j][k][g->top[i][j][k]] = l;
        /* count of current atoms in this grid cell */
        g->top[i][j][k]++;

#if defined(DEBUG_FOCUS)
        fprintf( stderr, "[INFO] Bin_Atoms: atom id = %-6d, x = (%8.3f, %8.3f, %8.3f), bin = (%3d, %3d, %3d)\n",
                l, system->atoms[l].x[0], system->atoms[l].x[1], system->atoms[l].x[2],
                i, j, k );
#endif
    }

    /* find max number of atoms per cell across all cells */
    max_atoms = 0;
    for ( i = 0; i < g->ncell[0]; i++ )
    {
        for ( j = 0; j < g->ncell[1]; j++ )
        {
            for ( k = 0; k < g->ncell[2]; k++ )
            {
                if ( max_atoms < g->top[i][j][k] )
                {
                    max_atoms = g->top[i][j][k];
                }
            }
        }
    }

    /* reallocation check */
    if ( max_atoms >= (int) CEIL( g->max_atoms * SAFE_ZONE ) )
    {
        workspace->realloc.gcell_atoms = MAX( (int) CEIL( max_atoms * SAFE_ZONE ),
                MIN_GCELL_POPL );
    }
}


void Finalize_Grid( reax_system * const system )
{
    if ( system->g.allocated == TRUE )
    {
        Deallocate_Grid_Space( &system->g );
    }
}


static void reax_atom_Copy( reax_atom * const dest, reax_atom * const src )
{
    dest->type = src->type;
    strncpy( dest->name, src->name, sizeof(dest->name) - 1 );
    dest->name[sizeof(dest->name) - 1] = '\0';
    rvec_Copy( dest->x, src->x );
    rvec_Copy( dest->v, src->v );
    rvec_Copy( dest->f, src->f );
    dest->q = src->q;
#if defined(QMMM)
    dest->qmmm_mask = src->qmmm_mask;
    dest->q_init = src->q_init;
#endif
}


static void Copy_Storage( reax_system const * const system,static_storage * const workspace,
        control_params const * const control, int top, int old_id, int old_type,
        real ** const v, real ** const s, real ** const t, int * const orig_id,
        rvec * const f_old )
{
    int i;

    for ( i = 0; i < control->cm_solver_restart + 1; ++i )
    {
        v[i][top] = workspace->v[i][old_id];
    }

    for ( i = 0; i < 5; ++i )
    {
        s[i][top] = workspace->s[i][old_id];
        t[i][top] = workspace->t[i][old_id];
    }

    orig_id[top] = workspace->orig_id[old_id];

    workspace->b_s[top] = -system->reax_param.sbp[ old_type ].chi;
    workspace->b_t[top] = -1.0;

    rvec_Copy( f_old[top], workspace->f_old[old_id] );
}


static void Free_Storage( static_storage * const workspace,
        control_params const * const control )
{
    int i;

    for ( i = 0; i < control->cm_solver_restart + 1; ++i )
    {
        sfree( workspace->v[i], __FILE__, __LINE__ );
    }
    sfree( workspace->v, __FILE__, __LINE__ );

    for ( i = 0; i < 3; ++i )
    {
        sfree( workspace->s[i], __FILE__, __LINE__ );
        sfree( workspace->t[i], __FILE__, __LINE__ );
    }
    sfree( workspace->s, __FILE__, __LINE__ );
    sfree( workspace->t, __FILE__, __LINE__ );

    sfree( workspace->orig_id, __FILE__, __LINE__ );
}


static void Assign_New_Storage( static_storage *workspace,
        real **v, real **s, real **t, int *orig_id, rvec *f_old )
{
    workspace->v = v;
    workspace->s = s;
    workspace->t = t;
    workspace->orig_id = orig_id;
    workspace->f_old = f_old;
}


/* Reorder atom list to improve cache performance */
void Reorder_Atoms( reax_system * const system, static_storage * const workspace,
        control_params const * const control )
{
    int i, j, k, l, top, old_id;
    reax_atom *old_atom, *new_atoms;
    grid *g;
    int *orig_id;
    real **v;
    real **s, **t;
    rvec *f_old;

    top = 0;
    g = &system->g;

    new_atoms = scalloc( system->N, sizeof(reax_atom), __FILE__, __LINE__ );
    orig_id = scalloc( system->N, sizeof(int), __FILE__, __LINE__ );
    f_old = scalloc( system->N, sizeof(rvec), __FILE__, __LINE__ );

    s = scalloc( 5, sizeof(real *), __FILE__, __LINE__ );
    t = scalloc( 5, sizeof(real *), __FILE__, __LINE__ );
    for ( i = 0; i < 5; ++i )
    {
        s[i] = scalloc( system->N_cm, sizeof(real), __FILE__, __LINE__ );
        t[i] = scalloc( system->N_cm, sizeof(real), __FILE__, __LINE__ );
    }

    v = scalloc( control->cm_solver_restart + 1, sizeof(real *),
            __FILE__, __LINE__ );
    for ( i = 0; i < control->cm_solver_restart + 1; ++i )
    {
        v[i] = scalloc( system->N_cm, sizeof(real), __FILE__, __LINE__ );
    }

    for ( i = 0; i < g->ncell[0]; i++ )
    {
        for ( j = 0; j < g->ncell[1]; j++ )
        {
            for ( k = 0; k < g->ncell[2]; k++ )
            {
                g->start[i][j][k] = top;

                for ( l = 0; l < g->top[i][j][k]; ++l )
                {
                    old_id = g->atoms[i][j][k][l];
                    old_atom = &system->atoms[old_id];

                    reax_atom_Copy( &new_atoms[top], old_atom );
                    Copy_Storage( system, workspace, control, top, old_id, old_atom->type,
                            v, s, t, orig_id, f_old );

                    ++top;
                }

                g->end[i][j][k] = top;
            }
        }
    }

    sfree( system->atoms, __FILE__, __LINE__ );
    Free_Storage( workspace, control );

    system->atoms = new_atoms;
    Assign_New_Storage( workspace, v, s, t, orig_id, f_old );
}
