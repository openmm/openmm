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

#include "bond_orders.h"

#include "list.h"
#include "lookup.h"
#include "vector.h"


static inline real Cf45( real p1, real p2 )
{
    return  -EXP(-p2 / 2.0) /
        ( SQR( EXP(-p1 / 2.0) + EXP(p1 / 2.0) ) * (EXP(-p2 / 2.0) + EXP(p2 / 2.0)) );
}


#if defined(TEST_FORCES)
void Add_dBO( reax_system *system, reax_list **lists,
        int i, int pj, real C, rvec *v )
{
    int start_pj, end_pj, k;
    reax_list *bond_list, *dBOs;

    bond_list = lists[BONDS];
    dBOs = lists[DBO];
    pj = bond_list->bond_list[pj].dbond_index;
    start_pj = Start_Index( pj, dBOs );
    end_pj = End_Index( pj, dBOs );

#if defined(DEBUG_FOCUS)
    fprintf( stderr, "[Add_dBO] i: %d, j: %d, start: %d, end: %d\n", i, pj, start_pj, end_pj );
#endif

    for ( k = start_pj; k < end_pj; ++k )
    {
        rvec_ScaledAdd( v[dBOs->dbo_list[k].wrt],
                C, dBOs->dbo_list[k].dBO );
    }
}


void Add_dBOpinpi2( reax_system *system, reax_list **lists,
        int i, int pj, real Cpi, real Cpi2, rvec *vpi, rvec *vpi2 )
{
    int start_pj, end_pj, k;
    reax_list *bond_list, *dBOs;
    dbond_data *dbo_k;

    bond_list = lists[BONDS];
    dBOs = lists[DBO];
    pj = bond_list->bond_list[pj].dbond_index;
    start_pj = Start_Index( pj, dBOs );
    end_pj = End_Index( pj, dBOs );

#if defined(DEBUG_FOCUS)
    fprintf( stderr, "[Add_dBOpinpi2] i: %d, j: %d, start: %d, end: %d\n", i, pj, start_pj, end_pj );
#endif

    for ( k = start_pj; k < end_pj; ++k )
    {
        dbo_k = &dBOs->dbo_list[k];
        rvec_ScaledAdd( vpi[dbo_k->wrt], Cpi, dbo_k->dBOpi );
        rvec_ScaledAdd( vpi2[dbo_k->wrt], Cpi2, dbo_k->dBOpi2 );
    }
}


void Add_dBO_to_Forces( reax_system *system, reax_list **lists,
                        int i, int pj, real C )
{
    reax_list *bond_list;
    reax_list *dBOs;
    int start_pj, end_pj, k;

    bond_list = lists[BONDS];
    dBOs = lists[DBO];
    pj = bond_list->bond_list[pj].dbond_index;
    start_pj = Start_Index( pj, dBOs );
    end_pj = End_Index( pj, dBOs );

    for ( k = start_pj; k < end_pj; ++k )
    {
        rvec_ScaledAdd( system->atoms[dBOs->dbo_list[k].wrt].f,
                C, dBOs->dbo_list[k].dBO );
    }
}


void Add_dBOpinpi2_to_Forces( reax_system *system, reax_list **lists,
        int i, int pj, real Cpi, real Cpi2 )
{
    reax_list *bond_list;
    reax_list *dBOs;
    dbond_data *dbo_k;
    int start_pj, end_pj, k;

    bond_list = lists[BONDS];
    dBOs = lists[DBO];
    pj = bond_list->bond_list[pj].dbond_index;
    start_pj = Start_Index(pj, dBOs);
    end_pj = End_Index(pj, dBOs);

    for ( k = start_pj; k < end_pj; ++k )
    {
        dbo_k = &dBOs->dbo_list[k];
        rvec_ScaledAdd( system->atoms[dbo_k->wrt].f, Cpi, dbo_k->dBOpi );
        rvec_ScaledAdd( system->atoms[dbo_k->wrt].f, Cpi2, dbo_k->dBOpi2 );
    }
}


void Add_dDelta( reax_system *system, reax_list **lists, int i, real C, rvec *v )
{
    int start, end, k;
    reax_list *dDeltas;

    dDeltas = lists[DDELTA];
    start = Start_Index( i, dDeltas );
    end = End_Index( i, dDeltas );

#if defined(DEBUG_FOCUS)
    fprintf( stderr, "[Add_dDelta] i: %d, start: %d, end: %d\n", i, start, end );
#endif

    for ( k = start; k < end; ++k )
    {
        rvec_ScaledAdd( v[dDeltas->dDelta_list[k].wrt],
                C, dDeltas->dDelta_list[k].dVal );
    }
}


void Add_dDelta_to_Forces( reax_system *system, reax_list **lists, int i, real C )
{
    reax_list *dDeltas;
    int start, end, k;

    dDeltas = lists[DDELTA];
    start = Start_Index( i, dDeltas );
    end = End_Index( i, dDeltas );

    for ( k = start; k < end; ++k )
    {
        rvec_ScaledAdd( system->atoms[dDeltas->dDelta_list[k].wrt].f,
                C, dDeltas->dDelta_list[k].dVal );
    }
}


void Calculate_dBO( int i, int pj, static_storage *workspace, reax_list **lists,
        int *top )
{
    int j, k, l, start_i, end_i, end_j;
    reax_list *bond_list, *dBOs;
    bond_data *nbr_l, *nbr_k;
    bond_order_data *bo_ij;
    dbond_data *top_dbo;

    bond_list = lists[BONDS];
    dBOs = lists[DBO];
    j = bond_list->bond_list[pj].nbr;
    bo_ij = &bond_list->bond_list[pj].bo_data;
    start_i = Start_Index( i, bond_list );
    end_i = End_Index( i, bond_list );
    l = Start_Index( j, bond_list );
    end_j = End_Index( j, bond_list );
    top_dbo = &dBOs->dbo_list[ *top ];

//    fprintf( stderr, "[Calculate_dBO] i: %d, pj: %d, start_i: %d, end_i: %d, start_j: %d, end_j: %d\n", i, pj, start_i, end_i, l, end_j );

    for ( k = start_i; k < end_i; ++k )
    {
        nbr_k = &bond_list->bond_list[k];

        for ( ; l < end_j && bond_list->bond_list[l].nbr < nbr_k->nbr; ++l )
        {
            /* These are the neighbors of j which aren't in the neighbor_list of i
             * Note that they might also include i! */
            nbr_l = &bond_list->bond_list[l];
            top_dbo->wrt = nbr_l->nbr;

            rvec_Scale( top_dbo->dBO, -bo_ij->C3dbo, nbr_l->bo_data.dBOp );
            rvec_Scale( top_dbo->dBOpi, -bo_ij->C4dbopi, nbr_l->bo_data.dBOp );
            rvec_Scale( top_dbo->dBOpi2, -bo_ij->C4dbopi2, nbr_l->bo_data.dBOp );

            if ( nbr_l->nbr == i )
            {
                /* dBO */
                rvec_ScaledAdd( top_dbo->dBO, bo_ij->C1dbo, bo_ij->dBOp );
                rvec_ScaledAdd( top_dbo->dBO, bo_ij->C2dbo, workspace->dDeltap_self[i] );

                /* dBOpi */
                rvec_ScaledAdd( top_dbo->dBOpi, bo_ij->C1dbopi, bo_ij->dln_BOp_pi );
                rvec_ScaledAdd( top_dbo->dBOpi, bo_ij->C2dbopi, bo_ij->dBOp );
                rvec_ScaledAdd( top_dbo->dBOpi, bo_ij->C3dbopi, workspace->dDeltap_self[i] );

                /* dBOpp */
                rvec_ScaledAdd( top_dbo->dBOpi2, bo_ij->C1dbopi2, bo_ij->dln_BOp_pi2 );
                rvec_ScaledAdd( top_dbo->dBOpi2, bo_ij->C2dbopi2, bo_ij->dBOp );
                rvec_ScaledAdd( top_dbo->dBOpi2, bo_ij->C3dbopi2, workspace->dDeltap_self[i] );
            }

            ++(*top);
            ++top_dbo;
        }

        /* Now we are processing neighbor k of i. */
        top_dbo->wrt = nbr_k->nbr;

        rvec_Scale( top_dbo->dBO, -bo_ij->C2dbo, nbr_k->bo_data.dBOp );      //dBO-2
        rvec_Scale( top_dbo->dBOpi, -bo_ij->C3dbopi, nbr_k->bo_data.dBOp );  //dBOpi-3
        rvec_Scale( top_dbo->dBOpi2, -bo_ij->C3dbopi2, nbr_k->bo_data.dBOp );//dBOpp-3

        if ( l < end_j && bond_list->bond_list[l].nbr == nbr_k->nbr )
        {
            /* common neighbor of i and j */
            nbr_l = &bond_list->bond_list[l];
            rvec_Copy( dBOp, nbr_l->bo_data.dBOp );

            rvec_ScaledAdd( top_dbo->dBO, -bo_ij->C3dbo, nbr_l->bo_data.dBOp );      //dBO,3rd
            rvec_ScaledAdd( top_dbo->dBOpi, -bo_ij->C4dbopi, nbr_l->bo_data.dBOp );  //dBOpi,4th
            rvec_ScaledAdd( top_dbo->dBOpi2, -bo_ij->C4dbopi2, nbr_l->bo_data.dBOp );//dBOpp.4th
            ++l;
        }
        else if ( k == pj )
        {
            /* 1st, dBO */
            rvec_ScaledAdd( top_dbo->dBO, -bo_ij->C1dbo, bo_ij->dBOp );
            /* 3rd, dBO */
            rvec_ScaledAdd( top_dbo->dBO, bo_ij->C3dbo, workspace->dDeltap_self[j] );

            /* dBOpi, 1st */
            rvec_ScaledAdd( top_dbo->dBOpi, -bo_ij->C1dbopi, bo_ij->dln_BOp_pi );
            /* 2nd */
            rvec_ScaledAdd( top_dbo->dBOpi, -bo_ij->C2dbopi, bo_ij->dBOp );
            /* 4th */
            rvec_ScaledAdd( top_dbo->dBOpi, bo_ij->C4dbopi, workspace->dDeltap_self[j] );

            /* dBOpi2, 1st */
            rvec_ScaledAdd( top_dbo->dBOpi2, -bo_ij->C1dbopi2, bo_ij->dln_BOp_pi2 );
            /* 2nd */
            rvec_ScaledAdd( top_dbo->dBOpi2, -bo_ij->C2dbopi2, bo_ij->dBOp );
            /* 4th */
            rvec_ScaledAdd( top_dbo->dBOpi2, bo_ij->C4dbopi2, workspace->dDeltap_self[j] );
        }

        ++(*top);
        ++top_dbo;
    }

    for ( ; l < end_j; ++l )
    {
        /* These are the remaining neighbors of j which are not in the
           neighbor_list of i. Note that they might also include i!*/
        nbr_l = &bond_list->bond_list[l];
        top_dbo->wrt = nbr_l->nbr;
        rvec_Copy( dBOp, nbr_l->bo_data.dBOp );

        rvec_Scale( top_dbo->dBO, -bo_ij->C3dbo, nbr_l->bo_data.dBOp );      //3rd, dBO
        rvec_Scale( top_dbo->dBOpi, -bo_ij->C4dbopi, nbr_l->bo_data.dBOp );  //4th, dBOpi
        rvec_Scale( top_dbo->dBOpi2, -bo_ij->C4dbopi2, nbr_l->bo_data.dBOp );//4th, dBOpp

        if ( nbr_l->nbr == i )
        {
            /* dBO, 1st */
            rvec_ScaledAdd( top_dbo->dBO, bo_ij->C1dbo, bo_ij->dBOp );
            rvec_ScaledAdd( top_dbo->dBO, bo_ij->C2dbo, workspace->dDeltap_self[i] ); //2nd, dBO

            /* dBOpi, 1st */
            rvec_ScaledAdd( top_dbo->dBOpi, bo_ij->C1dbopi, bo_ij->dln_BOp_pi );
            rvec_ScaledAdd( top_dbo->dBOpi, bo_ij->C2dbopi, bo_ij->dBOp );  //2nd
            rvec_ScaledAdd( top_dbo->dBOpi, bo_ij->C3dbopi, workspace->dDeltap_self[i] ); //3rd

            /* dBOpipi, 1st */
            rvec_ScaledAdd(top_dbo->dBOpi2, bo_ij->C1dbopi2, bo_ij->dln_BOp_pi2);
            rvec_ScaledAdd( top_dbo->dBOpi2, bo_ij->C2dbopi2, bo_ij->dBOp ); //2nd
            rvec_ScaledAdd( top_dbo->dBOpi2, bo_ij->C3dbopi2, workspace->dDeltap_self[i] );//3rd
        }

        ++(*top);
        ++top_dbo;
    }
}
#endif


void Add_dBond_to_Forces_NPT( int i, int pj, reax_system *system,
        simulation_data *data, static_storage *workspace, reax_list **lists )
{
    reax_list *bond_list;
    bond_data *nbr_j, *nbr_k;
    bond_order_data *bo_ij, *bo_ji;
    dbond_coefficients coef;
    rvec force, dvec_jk;
    rtensor press;
    int pk, k, j;
    rvec *f_i, *f_j, *f_k;
#if defined(_OPENMP)
    int tid = omp_get_thread_num( );
#endif

    bond_list = lists[BONDS];
    nbr_j = &bond_list->bond_list[pj];
    j = nbr_j->nbr;
    bo_ij = &nbr_j->bo_data;
    bo_ji = &bond_list->bond_list[ nbr_j->sym_index ].bo_data;
#if defined(_OPENMP)
    f_i = &workspace->f_local[tid * system->N + i];
    f_j = &workspace->f_local[tid * system->N + j];
#else
    f_i = &system->atoms[i].f;
    f_j = &system->atoms[j].f;
#endif

    coef.C1dbo = bo_ij->C1dbo * (bo_ij->Cdbo + bo_ji->Cdbo);
    coef.C2dbo = bo_ij->C2dbo * (bo_ij->Cdbo + bo_ji->Cdbo);
    coef.C3dbo = bo_ij->C3dbo * (bo_ij->Cdbo + bo_ji->Cdbo);

    coef.C1dbopi = bo_ij->C1dbopi * (bo_ij->Cdbopi + bo_ji->Cdbopi);
    coef.C2dbopi = bo_ij->C2dbopi * (bo_ij->Cdbopi + bo_ji->Cdbopi);
    coef.C3dbopi = bo_ij->C3dbopi * (bo_ij->Cdbopi + bo_ji->Cdbopi);
    coef.C4dbopi = bo_ij->C4dbopi * (bo_ij->Cdbopi + bo_ji->Cdbopi);

    coef.C1dbopi2 = bo_ij->C1dbopi2 * (bo_ij->Cdbopi2 + bo_ji->Cdbopi2);
    coef.C2dbopi2 = bo_ij->C2dbopi2 * (bo_ij->Cdbopi2 + bo_ji->Cdbopi2);
    coef.C3dbopi2 = bo_ij->C3dbopi2 * (bo_ij->Cdbopi2 + bo_ji->Cdbopi2);
    coef.C4dbopi2 = bo_ij->C4dbopi2 * (bo_ij->Cdbopi2 + bo_ji->Cdbopi2);

    coef.C1dDelta = bo_ij->C1dbo * (workspace->CdDelta[i] + workspace->CdDelta[j]);
    coef.C2dDelta = bo_ij->C2dbo * (workspace->CdDelta[i] + workspace->CdDelta[j]);
    coef.C3dDelta = bo_ij->C3dbo * (workspace->CdDelta[i] + workspace->CdDelta[j]);

    /************************************
    * forces related to atom i          *
    * first neighbors of atom i         *
    ************************************/
    for ( pk = Start_Index(i, bond_list); pk < End_Index(i, bond_list); ++pk )
    {
        nbr_k = &bond_list->bond_list[pk];
        k = nbr_k->nbr;
#if defined(_OPENMP)
        f_k = &workspace->f_local[tid * system->N + k];
#else
        f_k = &system->atoms[k].f;
#endif

        /* 2nd, dBO */
        rvec_Scale( force, -coef.C2dbo, nbr_k->bo_data.dBOp );
        /* dDelta */
        rvec_ScaledAdd( force, -coef.C2dDelta, nbr_k->bo_data.dBOp );
        /* 3rd, dBOpi */
        rvec_ScaledAdd( force, -coef.C3dbopi, nbr_k->bo_data.dBOp );
        /* 3rd, dBOpi2 */
        rvec_ScaledAdd( force, -coef.C3dbopi2, nbr_k->bo_data.dBOp );

        /* force */
        rvec_Add( *f_k, force );

        /* pressure */
        rvec_Scale( force, -1.0, force );
        rvec_OuterProduct( press, force, nbr_k->dvec  );
#if !defined(_OPENMP)
        rtensor_Add( data->press, press );
#else
        rtensor_Add( data->press_local[tid], press );
#endif
//        fprintf( stderr, "[BO3, i = %5d, j = %5d], %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", i, k, force[0], force[1], force[2], nbr_k->dvec[0], nbr_k->dvec[1], nbr_k->dvec[2] ); fflush( stderr ); 

        rvec_ScaledSum( dvec_jk, 1.0, nbr_k->dvec, -1.0, nbr_j->dvec );
        rvec_OuterProduct( press, dvec_jk, force );
#if !defined(_OPENMP)
        rtensor_Add( data->press, press );
#else
        rtensor_Add( data->press_local[tid], press );
#endif
//        fprintf( stderr, "[BO4, i = %5d, j = %5d], %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", j, k, force[0], force[1], force[2], dvec_jk[0], dvec_jk[1], dvec_jk[2] ); fflush( stderr ); 
    }

    /* then atom i itself */
    /* 1st, dBO */
    rvec_Scale( force, coef.C1dbo, bo_ij->dBOp );
    /* 2nd, dBO */
    rvec_ScaledAdd( force, coef.C2dbo, workspace->dDeltap_self[i] );

    /* 1st, dBO */
    rvec_ScaledAdd( force, coef.C1dDelta, bo_ij->dBOp );
    /* 2nd, dBO */
    rvec_ScaledAdd( force, coef.C2dDelta, workspace->dDeltap_self[i] );

    /* 1st, dBOpi */
    rvec_ScaledAdd( force, coef.C1dbopi, bo_ij->dln_BOp_pi );
    /* 2nd, dBOpi */
    rvec_ScaledAdd( force, coef.C2dbopi, bo_ij->dBOp );
    /* 3rd, dBOpi */
    rvec_ScaledAdd( force, coef.C3dbopi, workspace->dDeltap_self[i] );

    /* 1st, dBO_pi2 */
    rvec_ScaledAdd( force, coef.C1dbopi2, bo_ij->dln_BOp_pi2 );
    /* 2nd, dBO_pi2 */
    rvec_ScaledAdd( force, coef.C2dbopi2, bo_ij->dBOp );
    /* 3rd, dBO_pi2 */
    rvec_ScaledAdd( force, coef.C3dbopi2, workspace->dDeltap_self[i] );

    /* force */
    rvec_Add( *f_i, force );

    /* pressure */
    rvec_OuterProduct( press, force, nbr_j->dvec );
#if !defined(_OPENMP)
    rtensor_Add( data->press, press );
#else
    rtensor_Add( data->press_local[tid], press );
#endif
//    fprintf( stderr, "[BO1, i = %5d, j = %5d], %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", i, j, force[0], force[1], force[2], nbr_j->dvec[0], nbr_j->dvec[1], nbr_j->dvec[2] ); fflush( stderr ); 

    /****************************************************************************
     * forces and pressure related to atom j                                    *
     * first neighbors of atom j                                                *
     ***************************************************************************/
    for ( pk = Start_Index(j, bond_list); pk < End_Index(j, bond_list); ++pk )
    {
        nbr_k = &bond_list->bond_list[pk];
        k = nbr_k->nbr;
#if defined(_OPENMP)
        f_k = &workspace->f_local[tid * system->N + k];
#else
        f_k = &system->atoms[k].f;
#endif

        /* 3rd, dBO */
        rvec_Scale( force, -coef.C3dbo, nbr_k->bo_data.dBOp );
        /* dDelta */
        rvec_ScaledAdd( force, -coef.C3dDelta, nbr_k->bo_data.dBOp );
        /* 4th, dBOpi */
        rvec_ScaledAdd( force, -coef.C4dbopi, nbr_k->bo_data.dBOp );
        /* 4th, dBOpi2 */
        rvec_ScaledAdd( force, -coef.C4dbopi2, nbr_k->bo_data.dBOp );

        /* force */
        rvec_Add( *f_k, force );

        /* pressure */
        rvec_Scale( force, -1.0, force );
        rvec_Sum( dvec_jk, nbr_k->dvec, nbr_j->dvec );
        rvec_OuterProduct( press, force, dvec_jk );
#if !defined(_OPENMP)
        rtensor_Add( data->press, press );
#else
        rtensor_Add( data->press_local[tid], press );
#endif
//        fprintf( stderr, "[BO5, i = %5d, j = %5d], %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", i, k, force[0], force[1], force[2], dvec_jk[0], dvec_jk[1], dvec_jk[2] ); fflush( stderr ); 

        rvec_OuterProduct( press, force, nbr_k->dvec );
#if !defined(_OPENMP)
        rtensor_Add( data->press, press );
#else
        rtensor_Add( data->press_local[tid], press );
#endif
//        fprintf( stderr, "[BO6, i = %5d, j = %5d], %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", j, k, force[0], force[1], force[2], nbr_k->dvec[0], nbr_k->dvec[1], nbr_k->dvec[2] ); fflush( stderr ); 
    }

    /* then atom j itself */
    /* 1st, dBO */
    rvec_Scale( force, -coef.C1dbo, bo_ij->dBOp );
    /* 2nd, dBO */
    rvec_ScaledAdd( force, coef.C3dbo, workspace->dDeltap_self[j] );

    /* 1st, dBO */
    rvec_ScaledAdd( force, -coef.C1dDelta, bo_ij->dBOp );
    /* 2nd, dBO */
    rvec_ScaledAdd( force, coef.C3dDelta, workspace->dDeltap_self[j] );

    /* 1st, dBOpi */
    rvec_ScaledAdd( force, -coef.C1dbopi, bo_ij->dln_BOp_pi );
    /* 2nd, dBOpi */
    rvec_ScaledAdd( force, -coef.C2dbopi, bo_ij->dBOp );
    /* 3rd, dBOpi */
    rvec_ScaledAdd( force, coef.C4dbopi, workspace->dDeltap_self[j] );

    /* 1st, dBOpi2 */
    rvec_ScaledAdd( force, -coef.C1dbopi2, bo_ij->dln_BOp_pi2 );
    /* 2nd, dBOpi2 */
    rvec_ScaledAdd( force, -coef.C2dbopi2, bo_ij->dBOp );
    /* 3rd, dBOpi2 */
    rvec_ScaledAdd( force, coef.C4dbopi2, workspace->dDeltap_self[j] );

    /* force */
    rvec_Add( *f_j, force );

    /* pressure */
    rvec_Scale( force, -1.0, force );
    rvec_OuterProduct( press, force, nbr_j->dvec );
#if !defined(_OPENMP)
    rtensor_Add( data->press, press );
#else
    rtensor_Add( data->press_local[tid], press );
#endif
//    fprintf( stderr, "[BO2, i = %5d, j = %5d], %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", i, j, force[0], force[1], force[2], nbr_j->dvec[0], nbr_j->dvec[1], nbr_j->dvec[2] ); fflush( stderr ); 
}


void Add_dBond_to_Forces( int i, int pj, reax_system *system,
        simulation_data *data, static_storage *workspace, reax_list **lists )
{
    reax_list *bond_list;
    bond_data *nbr_j, *nbr_k;
    bond_order_data *bo_ij, *bo_ji;
    dbond_coefficients coef;
    int pk, k, j;
    rvec *f_i, *f_j, *f_k;
#if defined(_OPENMP)
    int tid = omp_get_thread_num( );
#endif

    bond_list = lists[BONDS];
    nbr_j = &bond_list->bond_list[pj];
    j = nbr_j->nbr;
    bo_ij = &nbr_j->bo_data;
    bo_ji = &bond_list->bond_list[ nbr_j->sym_index ].bo_data;
#if defined(_OPENMP)
    f_i = &workspace->f_local[tid * system->N + i];
    f_j = &workspace->f_local[tid * system->N + j];
#else
    f_i = &system->atoms[i].f;
    f_j = &system->atoms[j].f;
#endif

    coef.C1dbo = bo_ij->C1dbo * (bo_ij->Cdbo + bo_ji->Cdbo);
    coef.C2dbo = bo_ij->C2dbo * (bo_ij->Cdbo + bo_ji->Cdbo);
    coef.C3dbo = bo_ij->C3dbo * (bo_ij->Cdbo + bo_ji->Cdbo);

    coef.C1dbopi = bo_ij->C1dbopi * (bo_ij->Cdbopi + bo_ji->Cdbopi);
    coef.C2dbopi = bo_ij->C2dbopi * (bo_ij->Cdbopi + bo_ji->Cdbopi);
    coef.C3dbopi = bo_ij->C3dbopi * (bo_ij->Cdbopi + bo_ji->Cdbopi);
    coef.C4dbopi = bo_ij->C4dbopi * (bo_ij->Cdbopi + bo_ji->Cdbopi);

    coef.C1dbopi2 = bo_ij->C1dbopi2 * (bo_ij->Cdbopi2 + bo_ji->Cdbopi2);
    coef.C2dbopi2 = bo_ij->C2dbopi2 * (bo_ij->Cdbopi2 + bo_ji->Cdbopi2);
    coef.C3dbopi2 = bo_ij->C3dbopi2 * (bo_ij->Cdbopi2 + bo_ji->Cdbopi2);
    coef.C4dbopi2 = bo_ij->C4dbopi2 * (bo_ij->Cdbopi2 + bo_ji->Cdbopi2);

    coef.C1dDelta = bo_ij->C1dbo * (workspace->CdDelta[i] + workspace->CdDelta[j]);
    coef.C2dDelta = bo_ij->C2dbo * (workspace->CdDelta[i] + workspace->CdDelta[j]);
    coef.C3dDelta = bo_ij->C3dbo * (workspace->CdDelta[i] + workspace->CdDelta[j]);

    for ( pk = Start_Index(i, bond_list); pk < End_Index(i, bond_list); ++pk )
    {
        nbr_k = &bond_list->bond_list[pk];
        k = nbr_k->nbr;
#if defined(_OPENMP)
        f_k = &workspace->f_local[tid * system->N + k];
#else
        f_k = &system->atoms[k].f;
#endif

        /* 2nd, dBO */
        rvec_ScaledAdd( *f_k, -coef.C2dbo, nbr_k->bo_data.dBOp );
        /* dDelta */
        rvec_ScaledAdd( *f_k, -coef.C2dDelta, nbr_k->bo_data.dBOp );
        /* 3rd, dBOpi */
        rvec_ScaledAdd( *f_k, -coef.C3dbopi, nbr_k->bo_data.dBOp );
        /* 3rd, dBOpi2 */
        rvec_ScaledAdd( *f_k, -coef.C3dbopi2, nbr_k->bo_data.dBOp );
    }

    /* 1st, dBO */
    rvec_ScaledAdd( *f_i, coef.C1dbo, bo_ij->dBOp );
    /* 2nd, dBO */
    rvec_ScaledAdd( *f_i, coef.C2dbo, workspace->dDeltap_self[i] );

    /* 1st, dBO */
    rvec_ScaledAdd( *f_i, coef.C1dDelta, bo_ij->dBOp );
    /* 2nd, dBO */
    rvec_ScaledAdd( *f_i, coef.C2dDelta, workspace->dDeltap_self[i] );

    /* 1st, dBOpi */
    rvec_ScaledAdd( *f_i, coef.C1dbopi, bo_ij->dln_BOp_pi );
    /* 2nd, dBOpi */
    rvec_ScaledAdd( *f_i, coef.C2dbopi, bo_ij->dBOp );
    /* 3rd, dBOpi */
    rvec_ScaledAdd( *f_i, coef.C3dbopi, workspace->dDeltap_self[i] );

    /* 1st, dBO_pi2 */
    rvec_ScaledAdd( *f_i, coef.C1dbopi2, bo_ij->dln_BOp_pi2 );
    /* 2nd, dBO_pi2 */
    rvec_ScaledAdd( *f_i, coef.C2dbopi2, bo_ij->dBOp );
    /* 3rd, dBO_pi2 */
    rvec_ScaledAdd( *f_i, coef.C3dbopi2, workspace->dDeltap_self[i] );

    for ( pk = Start_Index(j, bond_list); pk < End_Index(j, bond_list); ++pk )
    {
        nbr_k = &bond_list->bond_list[pk];
        k = nbr_k->nbr;
#if defined(_OPENMP)
        f_k = &workspace->f_local[tid * system->N + k];
#else
        f_k = &system->atoms[k].f;
#endif

        /* 3rd, dBO */
        rvec_ScaledAdd( *f_k, -coef.C3dbo, nbr_k->bo_data.dBOp );
        /* dDelta */
        rvec_ScaledAdd( *f_k, -coef.C3dDelta, nbr_k->bo_data.dBOp );
        /* 4th, dBOpi */
        rvec_ScaledAdd( *f_k, -coef.C4dbopi, nbr_k->bo_data.dBOp );
        /* 4th, dBOpi2 */
        rvec_ScaledAdd( *f_k, -coef.C4dbopi2, nbr_k->bo_data.dBOp );
    }

    /* 1st, dBO */
    rvec_ScaledAdd( *f_j, -coef.C1dbo, bo_ij->dBOp );
    /* 2nd, dBO */
    rvec_ScaledAdd( *f_j, coef.C3dbo, workspace->dDeltap_self[j] );

    /* 1st, dBO */
    rvec_ScaledAdd( *f_j, -coef.C1dDelta, bo_ij->dBOp );
    /* 2nd, dBO */
    rvec_ScaledAdd( *f_j, coef.C3dDelta, workspace->dDeltap_self[j] );

    /* 1st, dBOpi */
    rvec_ScaledAdd( *f_j, -coef.C1dbopi, bo_ij->dln_BOp_pi );
    /* 2nd, dBOpi */
    rvec_ScaledAdd( *f_j, -coef.C2dbopi, bo_ij->dBOp );
    /* 3rd, dBOpi */
    rvec_ScaledAdd( *f_j, coef.C4dbopi, workspace->dDeltap_self[j] );

    /* 1st, dBOpi2 */
    rvec_ScaledAdd( *f_j, -coef.C1dbopi2, bo_ij->dln_BOp_pi2 );
    /* 2nd, dBOpi2 */
    rvec_ScaledAdd( *f_j, -coef.C2dbopi2, bo_ij->dBOp );
    /* 3rd, dBOpi2 */
    rvec_ScaledAdd( *f_j, coef.C4dbopi2, workspace->dDeltap_self[j] );
}


static inline void Copy_Neighbor_Data( bond_data *dest, near_neighbor_data *src )
{
    dest->nbr = src->nbr;
    dest->d = src->d;
    rvec_Copy( dest->dvec, src->dvec );
    ivec_Copy( dest->rel_box, src->rel_box );
    /* rvec_Copy( dest->ext_factor, src->ext_factor );*/
}


static inline void Copy_Bond_Order_Data( bond_order_data *dest, bond_order_data *src )
{
    dest->BO = src->BO;
    dest->BO_s = src->BO_s;
    dest->BO_pi = src->BO_pi;
    dest->BO_pi2 = src->BO_pi2;

    rvec_Scale( dest->dBOp, -1.0, src->dBOp );
    rvec_Scale( dest->dln_BOp_s, -1.0, src->dln_BOp_s );
    rvec_Scale( dest->dln_BOp_pi, -1.0, src->dln_BOp_pi );
    rvec_Scale( dest->dln_BOp_pi2, -1.0, src->dln_BOp_pi2 );
}


/* A very important and crucial assumption here is that each segment
 * belonging to a different atom in nbrhoods->nbr_list is sorted in its own.
 * This can either be done in the general coordinator function or here */
void BO( reax_system *system, control_params *control,
        simulation_data *data, static_storage *workspace,
        reax_list **lists, output_controls *out_control )
{
    real p_lp1;
    real p_boc1, p_boc2;
    reax_list *bond_list;
#if defined(TEST_FORCES)
    reax_list *dDeltas;
    reax_list *dBOs;
#endif

    p_lp1 = system->reax_param.gp.l[15];
    p_boc1 = system->reax_param.gp.l[0];
    p_boc2 = system->reax_param.gp.l[1];
    bond_list = lists[BONDS];
#if defined(TEST_FORCES)
    dDeltas = lists[DDELTA];
    dBOs = lists[DBO];
#endif

#if defined(_OPENMP)
    #pragma omp parallel default(shared)
#endif
    {
        int i, j, pj, type_i, type_j;
        int start_i, end_i;
        int sym_index;
        real val_i, Deltap_i, Deltap_boc_i;
        real val_j, Deltap_j, Deltap_boc_j;
        real temp, f1, f2, f3, f4, f5, f4f5, exp_f4, exp_f5;
        real exp_p1i, exp_p2i, exp_p1j, exp_p2j;
        real u1_ij, u1_ji, Cf1A_ij, Cf1B_ij, Cf1_ij, Cf1_ji;
        real Cf45_ij, Cf45_ji;
        real A0_ij, A1_ij, A2_ij, A2_ji, A3_ij, A3_ji;
        real explp1;
        two_body_parameters *twbp;
        bond_order_data *bo_ij, *bo_ji;
        single_body_parameters *sbp_i, *sbp_j;
#if defined(TEST_FORCES)
        int k, pk, start_j, end_j;
        int top_dbo, top_dDelta;
        dbond_data *pdbo;
        dDelta_data *ptop_dDelta;

        top_dbo = 0;
        top_dDelta = 0;
#endif

        /* Calculate Deltaprime, Deltaprime_boc values */
#if defined(_OPENMP)
        #pragma omp for schedule(static)
#endif
        for ( i = 0; i < system->N; ++i )
        {
            type_i = system->atoms[i].type;
            sbp_i = &system->reax_param.sbp[type_i];
            workspace->Deltap[i] = workspace->total_bond_order[i] - sbp_i->valency;
            workspace->Deltap_boc[i] =
                workspace->total_bond_order[i] - sbp_i->valency_val;
            workspace->total_bond_order[i] = 0.0;
        }

        /* wait until initialization complete */
#if defined(_OPENMP)
        #pragma omp barrier
#endif

        /* Corrected Bond Order calculations */
#if defined(_OPENMP)
        #pragma omp for schedule(guided)
#endif
        for ( i = 0; i < system->N; ++i )
        {
            type_i = system->atoms[i].type;
            sbp_i = &system->reax_param.sbp[type_i];
            val_i = sbp_i->valency;
            Deltap_i = workspace->Deltap[i];
            Deltap_boc_i = workspace->Deltap_boc[i];
            start_i = Start_Index(i, bond_list);
            end_i = End_Index(i, bond_list);

            for ( pj = start_i; pj < end_i; ++pj )
            {
                j = bond_list->bond_list[pj].nbr;
                type_j = system->atoms[j].type;
                bo_ij = &bond_list->bond_list[pj].bo_data;

                if ( i < j )
                {
                    twbp = &system->reax_param.tbp[type_i][type_j];

#if defined(TEST_FORCES)
                    Set_Start_Index( pj, top_dbo, dBOs );
#endif

                    if ( twbp->ovc < 0.001 && twbp->v13cor < 0.001 )
                    {
                        /* There is no correction to bond orders nor to derivatives of
                         * bond order prime! So we leave bond orders unchanged and
                         * set derivative of bond order coefficients s.t.
                         * dBO = dBOp & dBOxx = dBOxxp in Add_dBO_to_Forces */
                        bo_ij->C1dbo = 1.0;
                        bo_ij->C2dbo = 0.0;
                        bo_ij->C3dbo = 0.0;

                        bo_ij->C1dbopi = 1.0;
                        bo_ij->C2dbopi = 0.0;
                        bo_ij->C3dbopi = 0.0;
                        bo_ij->C4dbopi = 0.0;

                        bo_ij->C1dbopi2 = 1.0;
                        bo_ij->C2dbopi2 = 0.0;
                        bo_ij->C3dbopi2 = 0.0;
                        bo_ij->C4dbopi2 = 0.0;

#if defined(TEST_FORCES)
                        pdbo = &dBOs->dbo_list[ top_dbo ];

                        /* compute dBO_ij/dr_i */
                        pdbo->wrt = i;
                        rvec_Copy( pdbo->dBO, bo_ij->dBOp );
                        rvec_Scale( pdbo->dBOpi, bo_ij->BO_pi, bo_ij->dln_BOp_pi );
                        rvec_Scale( pdbo->dBOpi2, bo_ij->BO_pi2, bo_ij->dln_BOp_pi2 );

                        /* compute dBO_ij/dr_j */
                        pdbo = &dBOs->dbo_list[ top_dbo + 1 ];
                        pdbo->wrt = j;
                        rvec_Scale( pdbo->dBO, -1.0, bo_ij->dBOp );
                        rvec_Scale( pdbo->dBOpi, -bo_ij->BO_pi, bo_ij->dln_BOp_pi );
                        rvec_Scale( pdbo->dBOpi2, -bo_ij->BO_pi2, bo_ij->dln_BOp_pi2 );

                        top_dbo += 2;
#endif
                    }
                    else
                    {
                        val_j = system->reax_param.sbp[type_j].valency;
                        Deltap_j = workspace->Deltap[j];
                        Deltap_boc_j = workspace->Deltap_boc[j];

                        /* on page 1 */
                        if ( twbp->ovc >= 0.001 )
                        {
                            /* Correction for overcoordination */
                            exp_p1i = EXP( -p_boc1 * Deltap_i );
                            exp_p2i = EXP( -p_boc2 * Deltap_i );
                            exp_p1j = EXP( -p_boc1 * Deltap_j );
                            exp_p2j = EXP( -p_boc2 * Deltap_j );

                            f2 = exp_p1i + exp_p1j;
                            f3 = -1.0 / p_boc2 * log( 0.5 * ( exp_p2i  + exp_p2j ) );
                            f1 = 0.5 * ( ( val_i + f2 ) / ( val_i + f2 + f3 )
                                    + ( val_j + f2 ) / ( val_j + f2 + f3 ) );

                            /* Now come the derivates */
                            /* Bond Order pages 5-7, derivative of f1 */
                            temp = f2 + f3;
                            u1_ij = val_i + temp;
                            u1_ji = val_j + temp;
                            Cf1A_ij = 0.5 * f3 * (1.0 / SQR( u1_ij ) + 1.0 / SQR( u1_ji ));
                            Cf1B_ij = -0.5 * (( u1_ij - f3 ) / SQR( u1_ij )
                                    + ( u1_ji - f3 ) / SQR( u1_ji ));

                            //Cf1_ij = -Cf1A_ij * p_boc1 * exp_p1i +
                            //          Cf1B_ij * exp_p2i / ( exp_p2i + exp_p2j );
                            Cf1_ij = 0.50 * ( -p_boc1 * exp_p1i / u1_ij -
                                    ((val_i + f2) / SQR(u1_ij)) * ( -p_boc1 * exp_p1i +
                                    exp_p2i / ( exp_p2i + exp_p2j ) ) + -p_boc1 * exp_p1i / u1_ji -
                                    ((val_j + f2) / SQR(u1_ji)) * ( -p_boc1 * exp_p1i +
                                    exp_p2i / ( exp_p2i + exp_p2j ) ));

                            Cf1_ji = -Cf1A_ij * p_boc1 * exp_p1j
                                + Cf1B_ij * exp_p2j / ( exp_p2i + exp_p2j );
                        }
                        else
                        {
                            /* No overcoordination correction! */
                            f1 = 1.0;
                            Cf1_ij = 0.0;
                            Cf1_ji = 0.0;
                        }

                        if ( twbp->v13cor >= 0.001 )
                        {
                            /* Correction for 1-3 bond orders */
                            exp_f4 = EXP( -twbp->p_boc3 * (twbp->p_boc4 * SQR( bo_ij->BO ) - Deltap_boc_i)
                                    + twbp->p_boc5 );
                            exp_f5 = EXP( -twbp->p_boc3 * (twbp->p_boc4 * SQR( bo_ij->BO ) - Deltap_boc_j)
                                    + twbp->p_boc5 );

                            f4 = 1.0 / (1.0 + exp_f4);
                            f5 = 1.0 / (1.0 + exp_f5);
                            f4f5 = f4 * f5;

                            /* Bond Order pages 8-9, derivative of f4 and f5 */
//                            temp = twbp->p_boc5
//                                - twbp->p_boc3 * twbp->p_boc4 * SQR( bo_ij->BO );
//                            u_ij = temp + twbp->p_boc3 * Deltap_boc_i;
//                            u_ji = temp + twbp->p_boc3 * Deltap_boc_j;
//                            Cf45_ij = Cf45( u_ij, u_ji ) / f4f5;
//                            Cf45_ji = Cf45( u_ji, u_ij ) / f4f5;
                            Cf45_ij = -f4 * exp_f4;
                            Cf45_ji = -f5 * exp_f5;
                        }
                        else
                        {
                            f4 = 1.0;
                            f5 = 1.0;
                            f4f5 = 1.0;
                            Cf45_ij = 0.0;
                            Cf45_ji = 0.0;
                        }

                        /* Bond Order page 10, derivative of total bond order */
                        A0_ij = f1 * f4f5;
                        A1_ij = -2.0 * twbp->p_boc3 * twbp->p_boc4 * bo_ij->BO
                            * (Cf45_ij + Cf45_ji);
                        A2_ij = Cf1_ij / f1 + twbp->p_boc3 * Cf45_ij;
                        A2_ji = Cf1_ji / f1 + twbp->p_boc3 * Cf45_ji;
                        A3_ij = A2_ij + Cf1_ij / f1;
                        A3_ji = A2_ji + Cf1_ji / f1;

                        /* find corrected bond order values and their deriv coefs */
                        bo_ij->BO = bo_ij->BO * A0_ij;
                        bo_ij->BO_pi = bo_ij->BO_pi * A0_ij * f1;
                        bo_ij->BO_pi2 = bo_ij->BO_pi2 * A0_ij * f1;
                        bo_ij->BO_s = bo_ij->BO - ( bo_ij->BO_pi + bo_ij->BO_pi2 );

                        bo_ij->C1dbo = A0_ij + bo_ij->BO * A1_ij;
                        bo_ij->C2dbo = bo_ij->BO * A2_ij;
                        bo_ij->C3dbo = bo_ij->BO * A2_ji;

                        bo_ij->C1dbopi = f1 * f1 * f4 * f5;
                        bo_ij->C2dbopi = bo_ij->BO_pi * A1_ij;
                        bo_ij->C3dbopi = bo_ij->BO_pi * A3_ij;
                        bo_ij->C4dbopi = bo_ij->BO_pi * A3_ji;

                        bo_ij->C1dbopi2 = f1 * f1 * f4 * f5;
                        bo_ij->C2dbopi2 = bo_ij->BO_pi2 * A1_ij;
                        bo_ij->C3dbopi2 = bo_ij->BO_pi2 * A3_ij;
                        bo_ij->C4dbopi2 = bo_ij->BO_pi2 * A3_ji;

#if defined(TEST_FORCES)
                        Calculate_dBO( i, pj, workspace, lists, &top_dbo );
#endif
                    }

                    /* neglect weak bonds */
                    if ( bo_ij->BO < 1.0e-10 )
                    {
                        bo_ij->BO = 0.0;
                    }
                    if ( bo_ij->BO_s < 1.0e-10 )
                    {
                        bo_ij->BO_s = 0.0;
                    }
                    if ( bo_ij->BO_pi < 1.0e-10 )
                    {
                        bo_ij->BO_pi = 0.0;
                    }
                    if ( bo_ij->BO_pi2 < 1.0e-10 )
                    {
                        bo_ij->BO_pi2 = 0.0;
                    }

                    /* now keeps total_BO */
                    workspace->total_bond_order[i] += bo_ij->BO;

#if defined(TEST_FORCES)
                    Set_End_Index( pj, top_dbo, dBOs );
                    Add_dBO( system, lists, i, pj, 1.0, workspace->dDelta );
#endif
                }
            }

#if defined(TEST_FORCES)
            Set_Start_Index( i, top_dDelta, dDeltas );
            ptop_dDelta = &dDeltas->dDelta_list[top_dDelta];

            for ( pj = start_i; pj < end_i; ++pj )
            {
                j = bond_list->bond_list[pj].nbr;

                if ( !rvec_isZero( workspace->dDelta[j] ) )
                {
                    ptop_dDelta->wrt = j;
                    rvec_Copy( ptop_dDelta->dVal, workspace->dDelta[j] );
                    rvec_MakeZero( workspace->dDelta[j] );
                    ++top_dDelta;
                    ++ptop_dDelta;
                }

                start_j = Start_Index( j, bond_list );
                end_j = End_Index( j, bond_list );
                for ( pk = start_j; pk < end_j; ++pk )
                {
                    k = bond_list->bond_list[pk].nbr;
                    if ( !rvec_isZero( workspace->dDelta[k] ) )
                    {
                        ptop_dDelta->wrt = k;
                        rvec_Copy( ptop_dDelta->dVal, workspace->dDelta[k] );
                        rvec_MakeZero( workspace->dDelta[k] );
                        ++top_dDelta;
                        ++ptop_dDelta;
                    }
                }
            }

            Set_End_Index( i, top_dDelta, dDeltas );
#endif
        }

        /* wait for bo_ij to be updated */
#if defined(_OPENMP)
        #pragma omp barrier
#endif

#if defined(_OPENMP)
        #pragma omp for schedule(guided)
#endif
        for ( i = 0; i < system->N; ++i )
        {
            type_i = system->atoms[i].type;

            if ( type_i < 0 )
            {
                continue;
            }

            start_i = Start_Index(i, bond_list);
            end_i = End_Index(i, bond_list);

            for ( pj = start_i; pj < end_i; ++pj )
            {
                j = bond_list->bond_list[pj].nbr;
                type_j = system->atoms[j].type;

                if ( type_j < 0 )
                {
                    continue;
                }

                if ( i < j )
                {
                    /* computed in previous for-loop */
                    ;
                }
                else
                {
                    /* We only need to update bond orders from bo_ji
                     * everything else is set in uncorrected_bo calculations */
                    sym_index = bond_list->bond_list[pj].sym_index;

                    bo_ij = &bond_list->bond_list[ pj ].bo_data;
                    bo_ji = &bond_list->bond_list[ sym_index ].bo_data;
                    bo_ij->BO = bo_ji->BO;
                    bo_ij->BO_s = bo_ji->BO_s;
                    bo_ij->BO_pi = bo_ji->BO_pi;
                    bo_ij->BO_pi2 = bo_ji->BO_pi2;

                    /* now keeps total_BO */
                    workspace->total_bond_order[i] += bo_ij->BO;

#if defined(TEST_FORCES)
                    Add_dBO( system, lists, j, sym_index, 1.0, workspace->dDelta );
#endif
                }
            }
        }

        /* need to wait for total_bond_order to be accumulated */
#if defined(_OPENMP)
        #pragma omp barrier
#endif

        /* Calculate some helper variables that are  used at many places
         * throughout force calculations */
#if defined(_OPENMP)
        #pragma omp for schedule(guided)
#endif
        for ( j = 0; j < system->N; ++j )
        {
            type_j = system->atoms[j].type;
            sbp_j = &system->reax_param.sbp[ type_j ];

            workspace->Delta[j] = workspace->total_bond_order[j] - sbp_j->valency;
            workspace->Delta_e[j] = workspace->total_bond_order[j] - sbp_j->valency_e;
            workspace->Delta_boc[j] = workspace->total_bond_order[j]
                - sbp_j->valency_boc;

            workspace->vlpex[j] = workspace->Delta_e[j]
                - 2.0 * (int)(workspace->Delta_e[j] / 2.0);
            explp1 = EXP( -1.0 * p_lp1 * SQR(2.0 + workspace->vlpex[j]) );
            workspace->nlp[j] = explp1 - (int)(workspace->Delta_e[j] / 2.0);
            workspace->Delta_lp[j] = sbp_j->nlp_opt - workspace->nlp[j];
            workspace->Clp[j] = 2.0 * p_lp1 * (2.0 + workspace->vlpex[j]) * explp1;
            /* Adri uses different dDelta_lp values than the ones in notes... */
            workspace->dDelta_lp[j] = workspace->Clp[j];
//            workspace->dDelta_lp[j] = workspace->Clp[j] + (0.5 - workspace->Clp[j])
//                * ((FABS(workspace->Delta_e[j] / 2.0
//                                - (int)(workspace->Delta_e[j] / 2.0)) < 0.1) ? 1 : 0 );

            if ( sbp_j->mass > 21.0 )
            {
                workspace->nlp_temp[j] = 0.5 * (sbp_j->valency_e - sbp_j->valency);
                workspace->Delta_lp_temp[j] = sbp_j->nlp_opt - workspace->nlp_temp[j];
                workspace->dDelta_lp_temp[j] = 0.0;
            }
            else
            {
                workspace->nlp_temp[j] = workspace->nlp[j];
                workspace->Delta_lp_temp[j] = sbp_j->nlp_opt - workspace->nlp_temp[j];
                workspace->dDelta_lp_temp[j] = workspace->Clp[j];
            }
        }
    }

#if defined(DEBUG_FOCUS)
    Print_Bond_Orders( system, control, data, workspace, lists, out_control );
#endif
}
