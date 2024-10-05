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

#include "list.h"

#include "tool_box.h"


/* Allocates a new instance of reax_list with internal lists of specified type */
void Make_List( int n, int n_max, int total_intrs, int type, reax_list * const l )
{
    assert( n > 0 );
    assert( n_max > 0 );
    assert( n_max >= n );
    assert( total_intrs >= 0 );
    assert( l != NULL );

    if ( l->allocated == TRUE )
    {
        fprintf( stderr, "[WARNING] attempted to allocate list which was already allocated."
                " Returning without allocation...\n" );
        return;
    }

    l->allocated = TRUE;
    l->n = n;
    l->n_max = n_max;
    l->total_intrs = total_intrs;

    l->index = smalloc( n_max * sizeof(int), __FILE__, __LINE__ );
    l->end_index = smalloc( n_max * sizeof(int), __FILE__, __LINE__ );

    switch ( type )
    {
    case TYP_THREE_BODY:
        if ( l->total_intrs > 0 )
        {
            l->three_body_list = smalloc( l->total_intrs * sizeof(three_body_interaction_data),
                    __FILE__, __LINE__ );
        }
        else
        {
            l->three_body_list = NULL;
        }
        break;

    case TYP_BOND:
        if ( l->total_intrs > 0 )
        {
            l->bond_list = smalloc( l->total_intrs * sizeof(bond_data),
                    __FILE__, __LINE__ );
        }
        else
        {
            l->bond_list = NULL;
        }
        break;

    case TYP_DBO:
        if ( l->total_intrs > 0 )
        {
            l->dbo_list = smalloc( l->total_intrs * sizeof(dbond_data),
                    __FILE__, __LINE__ );
        }
        else
        {
            l->dbo_list = NULL;
        }
        break;

    case TYP_DDELTA:
        if ( l->total_intrs > 0 )
        {
            l->dDelta_list = smalloc( l->total_intrs * sizeof(dDelta_data),
                    __FILE__, __LINE__ );
        }
        else
        {
            l->dDelta_list = NULL;
        }
        break;

    case TYP_FAR_NEIGHBOR:
        if ( l->total_intrs > 0 )
        {
            l->far_nbr_list = smalloc( l->total_intrs * sizeof(far_neighbor_data),
                    __FILE__, __LINE__ );
        }
        else
        {
            l->far_nbr_list = NULL;
        }
        break;

    case TYP_NEAR_NEIGHBOR:
        if ( l->total_intrs > 0 )
        {
            l->near_nbr_list = smalloc( l->total_intrs * sizeof(near_neighbor_data),
                    __FILE__, __LINE__ );
        }
        else
        {
            l->near_nbr_list = NULL;
        }
        break;

    case TYP_HBOND:
        if ( l->total_intrs > 0 )
        {
            l->hbond_list = smalloc( l->total_intrs * sizeof(hbond_data),
                    __FILE__, __LINE__ );
        }
        else
        {
            l->hbond_list = NULL;
        }
        break;

    default:
        fprintf( stderr, "[ERROR] unknown list type. Terminating...\n" );
        exit( INVALID_INPUT );
        break;
    }
}


/* Deallocates any space allocated for a reax_list instance */
void Delete_List( int type, reax_list * const l )
{
    assert( l != NULL );

    if ( l->allocated == FALSE )
    {
        fprintf( stderr, "[WARNING] attempted to free list which was not allocated."
                " Returning without deallocation...\n" );
        return;
    }

    l->allocated = FALSE;
    l->n = 0;
    l->n_max = 0;
    l->total_intrs = 0;

    sfree( l->index, __FILE__, __LINE__ );
    sfree( l->end_index, __FILE__, __LINE__ );

    switch ( type )
    {
    case TYP_THREE_BODY:
        if ( l->three_body_list != NULL )
        {
            sfree( l->three_body_list, __FILE__, __LINE__ );
        }
        break;

    case TYP_BOND:
        if ( l->bond_list != NULL )
        {
            sfree( l->bond_list, __FILE__, __LINE__ );
        }
        break;

    case TYP_DBO:
        if ( l->dbo_list != NULL )
        {
            sfree( l->dbo_list, __FILE__, __LINE__ );
        }
        break;

    case TYP_DDELTA:
        if ( l->dDelta_list != NULL )
        {
            sfree( l->dDelta_list, __FILE__, __LINE__ );
        }
        break;

    case TYP_FAR_NEIGHBOR:
        if ( l->far_nbr_list != NULL )
        {
            sfree( l->far_nbr_list, __FILE__, __LINE__ );
        }
        break;

    case TYP_NEAR_NEIGHBOR:
        if ( l->near_nbr_list != NULL )
        {
            sfree( l->near_nbr_list, __FILE__, __LINE__ );
        }
        break;

    case TYP_HBOND:
        if ( l->hbond_list != NULL )
        {
            sfree( l->hbond_list, __FILE__, __LINE__ );
        }
        break;

    default:
        fprintf( stderr, "[ERROR] unknown list type. Terminating...\n" );
        exit( INVALID_INPUT );
        break;
    }
}


/* Initialize list indices
 *
 * l: pointer to list
 * max_intrs: max. num. of interactions for each list element
 * */
void Init_List_Indices( reax_list * const l, int * const max_intrs )
{
    int i;

    assert( l != NULL );
    assert( l->n > 0 );
    assert( max_intrs > 0 );

    /* exclusive prefix sum of max_intrs replaces start indices,
     * set end indices to the same as start indices for safety */
    Set_Start_Index( 0, 0, l );
    Set_End_Index( 0, 0, l );
    for ( i = 1; i < l->n; ++i )
    {
        Set_Start_Index( i, Start_Index( i - 1, l ) + max_intrs[i - 1], l );
        Set_End_Index( i, Start_Index( i, l ), l );
    }
}
