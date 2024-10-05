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

#include "hydrogen_bonds.h"

#include "bond_orders.h"
#include "list.h"
#include "lookup.h"
#include "valence_angles.h"
#include "vector.h"

#include "tool_box.h"


void Hydrogen_Bonds( reax_system *system, control_params *control,
        simulation_data *data, static_storage *workspace,
        reax_list **lists, output_controls *out_control )
{
#if defined(TEST_FORCES)
    int num_hb_intrs;
#endif
    real e_hb_total;

    e_hb_total = 0.0;
#if defined(TEST_FORCES)
    num_hb_intrs = 0;
#endif

#if defined(_OPENMP)
    #pragma omp parallel default(shared) reduction(+: e_hb_total)
#endif
    {
        int i, j, k, pi, pk, itr, top;
        int type_i, type_j, type_k;
        int start_j, end_j, hb_start_j, hb_end_j;
        int *hblist, hblist_size;
        real r_ij, r_jk, theta, cos_theta, sin_xhz4, cos_xhz1, sin_theta2;
        real e_hb, exp_hb2, exp_hb3, CEhb1, CEhb2, CEhb3;
        rvec dcos_theta_di, dcos_theta_dj, dcos_theta_dk;
        rvec dvec_jk, force_i, force_j, force_k, x_i, x_j, x_k;
        rtensor press;
//        rtensor temp_rtensor, total_rtensor;
        hbond_parameters *hbp;
        bond_order_data *bo_ij;
        bond_data *pbond_ij;
        far_neighbor_data *nbr_jk;
        reax_list *bonds, *hbonds;
        bond_data *bond_list;
        hbond_data *hbond_list;
        rvec *f_i, *f_j, *f_k;
#if defined(_OPENMP)
        int tid = omp_get_thread_num( );
#endif

        hblist = NULL;
        hblist_size = 0;
        bonds = lists[BONDS];
        bond_list = bonds->bond_list;
        hbonds = lists[HBONDS];
        hbond_list = hbonds->hbond_list;

        /* loops below discover the Hydrogen bonds between i-j-k triplets.
         * here j is H atom and there has to be some bond between i and j.
         * Hydrogen bond is between j and k.
         * so in this function i->X, j->H, k->Z when we map
         * variables onto the ones in the handout. */
#if defined(_OPENMP)
        #pragma omp for schedule(guided)
#endif
        for ( j = 0; j < system->N; ++j )
        {
            /* j must be a hydrogen atom */
            if ( system->reax_param.sbp[system->atoms[j].type].p_hbond == H_ATOM
#if defined(QMMM)
                    && system->atoms[j].qmmm_mask == TRUE
#endif
               )
            {
                type_j = system->atoms[j].type;
                start_j = Start_Index( j, bonds );
                end_j = End_Index( j, bonds );
                hb_start_j = Start_Index( j, hbonds );
                hb_end_j = End_Index( j, hbonds );
#if defined(_OPENMP)
                f_j = &workspace->f_local[tid * system->N + j];
#else
                f_j = &system->atoms[j].f;
#endif
                if ( control->ensemble == sNPT || control->ensemble == iNPT
                        || control->ensemble == aNPT || control->compute_pressure == TRUE )
                {
                    rvec_iMultiply( x_j, system->atoms[j].rel_map, system->box.box_norms );
                    rvec_Add( x_j, system->atoms[j].x );
                }

                top = 0;
                if ( Num_Entries( j, bonds ) > hblist_size )
                {
                    hblist_size = Num_Entries( j, bonds );
                    hblist = srealloc( hblist, sizeof(int) * hblist_size,
                            __FILE__, __LINE__ );
                }

                for ( pi = start_j; pi < end_j; ++pi )
                {
                    pbond_ij = &bond_list[pi];
                    i = pbond_ij->nbr;
                    bo_ij = &pbond_ij->bo_data;
                    type_i = system->atoms[i].type;

                    if ( system->reax_param.sbp[type_i].p_hbond == H_BONDING_ATOM
                            && bo_ij->BO >= HB_THRESHOLD )
                    {
                        hblist[top++] = pi;
                    }
                }

                for ( pk = hb_start_j; pk < hb_end_j; ++pk )
                {
                    k = hbond_list[pk].nbr;

#if defined(QMMM)
                    if ( system->atoms[k].qmmm_mask == TRUE )
                    {
#endif

                    type_k = system->atoms[k].type;
                    nbr_jk = hbond_list[pk].ptr;
                    r_jk = nbr_jk->d;
                    rvec_Scale( dvec_jk, hbond_list[pk].scl, nbr_jk->dvec );
#if defined(_OPENMP)
                    f_k = &workspace->f_local[tid * system->N + k];
#else
                    f_k = &system->atoms[k].f;
#endif

                    for ( itr = 0; itr < top; ++itr )
                    {
                        pi = hblist[itr];
                        pbond_ij = &bond_list[pi];
                        i = pbond_ij->nbr;

                        if ( i != k
                                && system->reax_param.hbp[system->atoms[i].type][type_j][type_k].is_valid == TRUE )
                        {
                            bo_ij = &pbond_ij->bo_data;
                            type_i = system->atoms[i].type;
                            r_ij = pbond_ij->d;
                            hbp = &system->reax_param.hbp[type_i][type_j][type_k];
#if defined(_OPENMP)
                            f_i = &workspace->f_local[tid * system->N + i];
#else
                            f_i = &system->atoms[i].f;
#endif

#if defined(TEST_FORCES)
#if defined(_OPENMP)
                            #pragma omp atomic
#endif
                            ++num_hb_intrs;
#endif

                            Calculate_Theta( pbond_ij->dvec, r_ij, dvec_jk, r_jk,
                                    &theta, &cos_theta );
                            /* the derivative of cos(theta) */
                            Calculate_dCos_Theta( pbond_ij->dvec, r_ij, dvec_jk, r_jk,
                                    &dcos_theta_di, &dcos_theta_dj, &dcos_theta_dk );

                            /* hydrogen bond energy */
                            sin_theta2 = SIN( theta / 2.0 );
                            sin_xhz4 = SQR( sin_theta2 );
                            sin_xhz4 *= sin_xhz4;
                            cos_xhz1 = ( 1.0 - cos_theta );
                            exp_hb2 = EXP( -hbp->p_hb2 * bo_ij->BO );
                            exp_hb3 = EXP( -hbp->p_hb3 * ( hbp->r0_hb / r_jk
                                        + r_jk / hbp->r0_hb - 2.0 ) );

                            e_hb = hbp->p_hb1 * (1.0 - exp_hb2) * exp_hb3 * sin_xhz4;
                            e_hb_total += e_hb;

                            CEhb1 = hbp->p_hb1 * hbp->p_hb2 * exp_hb2 * exp_hb3 * sin_xhz4;
                            CEhb2 = -hbp->p_hb1 / 2.0 * (1.0 - exp_hb2) * exp_hb3 * cos_xhz1;
                            CEhb3 = -hbp->p_hb3 * e_hb * (-hbp->r0_hb / SQR( r_jk )
                                    + 1.0 / hbp->r0_hb);

                            /* hydrogen bond forces */
                            /* dbo term,
                             * note: safe to update across threads as this points
                             * to the bond_order_data struct inside atom j's list,
                             * and threads are partitioned across all j's */
#if defined(_OPENMP)
                            #pragma omp atomic
#endif
                            bo_ij->Cdbo += CEhb1;

                            if ( control->compute_pressure == FALSE &&
                                    (control->ensemble == NVE || control->ensemble == nhNVT
                                     || control->ensemble == bNVT) )
                            {
                                /* dcos terms */
                                rvec_ScaledAdd( *f_i, CEhb2, dcos_theta_di );
                                rvec_ScaledAdd( *f_j, CEhb2, dcos_theta_dj );
                                rvec_ScaledAdd( *f_k, CEhb2, dcos_theta_dk );

                                /* dr terms */
                                rvec_ScaledAdd( *f_j, -CEhb3 / r_jk, dvec_jk );
                                rvec_ScaledAdd( *f_k, CEhb3 / r_jk, dvec_jk );
                            }
                            else if ( control->ensemble == sNPT || control->ensemble == iNPT
                                    || control->ensemble == aNPT || control->compute_pressure == TRUE )
                            {
                                /* for pressure coupling, terms that are not related
                                 * to bond order derivatives are added directly into
                                 * pressure vector/tensor */

                                /* dcos terms */
                                rvec_Scale( force_i, CEhb2, dcos_theta_di );
                                rvec_Scale( force_j, CEhb2, dcos_theta_dj );
                                rvec_Scale( force_k, CEhb2, dcos_theta_dk );

                                /* dr terms */
                                rvec_ScaledAdd( force_j, -CEhb3 / r_jk, dvec_jk );
                                rvec_ScaledAdd( force_k, CEhb3 / r_jk, dvec_jk );

                                rvec_Add( *f_i, force_i );
                                rvec_Add( *f_j, force_j );
                                rvec_Add( *f_k, force_k );

                                /* pressure */
                                rvec_Sum( x_i, x_j, pbond_ij->dvec );
                                rvec_OuterProduct( press, force_i, x_i );
#if !defined(_OPENMP)
                                rtensor_Add( data->press, press );
#else
                                rtensor_Add( data->press_local[tid], press );
#endif

                                rvec_OuterProduct( press, force_j, x_j );
#if !defined(_OPENMP)
                                rtensor_Add( data->press, press );
#else
                                rtensor_Add( data->press_local[tid], press );
#endif

                                rvec_Sum( x_k, x_j, nbr_jk->dvec );
                                rvec_OuterProduct( press, force_k, x_k );
#if !defined(_OPENMP)
                                rtensor_Add( data->press, press );
#else
                                rtensor_Add( data->press_local[tid], press );
#endif

                                /* This part is intended for a fully-flexible box */
//                                rvec_OuterProduct( temp_rtensor,
//                                        dcos_theta_di, system->atoms[i].x );
//                                rtensor_Scale( total_rtensor, -CEhb2, temp_rtensor );
//
//                                rvec_ScaledSum( temp_rvec, -CEhb2, dcos_theta_dj,
//                                        -CEhb3 / r_jk, pbond_jk->dvec );
//                                rvec_OuterProduct( temp_rtensor,
//                                        temp_rvec, system->atoms[j].x );
//                                rtensor_Add( total_rtensor, temp_rtensor );
//
//                                rvec_ScaledSum( temp_rvec, -CEhb2, dcos_theta_dk,
//                                        +CEhb3 / r_jk, pbond_jk->dvec );
//                                rvec_OuterProduct( temp_rtensor,
//                                        temp_rvec, system->atoms[k].x );
//                                rtensor_Add( total_rtensor, temp_rtensor );
//
//                                if ( pbond_ij->imaginary || pbond_jk->imaginary )
//                                {
//                                    rtensor_ScaledAdd( data->flex_bar.P, -1.0, total_rtensor );
//                                }
//                                else
//                                {
//                                    rtensor_Add( data->flex_bar.P, total_rtensor );
//                                }
                            }

#if defined(TEST_ENERGY)
                        /* fprintf( out_control->ehb,
                           "%24.15e%24.15e%24.15e\n%24.15e%24.15e%24.15e\n%24.15e%24.15e%24.15e\n",
                           dcos_theta_di[0], dcos_theta_di[1], dcos_theta_di[2],
                           dcos_theta_dj[0], dcos_theta_dj[1], dcos_theta_dj[2],
                           dcos_theta_dk[0], dcos_theta_dk[1], dcos_theta_dk[2]);
                           fprintf( out_control->ehb, "%24.15e%24.15e%24.15e\n",
                           CEhb1, CEhb2, CEhb3 ); */
                        fprintf( out_control->ehb,
                                 //"%6d%6d%6d%24.15e%24.15e%24.15e%24.15e%24.15e\n",
                                 "%6d%6d%6d%12.4f%12.4f%12.4f%12.4f%12.4f\n",
                                 system->my_atoms[i].orig_id, system->my_atoms[j].orig_id,
                                 system->my_atoms[k].orig_id,
                                 r_jk, theta, bo_ij->BO, e_hb, data->my_en.e_hb );
#endif

#if defined(TEST_FORCES)
                            /* dbo term */
                            Add_dBO( system, lists, j, pi, +CEhb1, workspace->f_hb );
                            /* dcos terms */
                            rvec_ScaledAdd( workspace->f_hb[i], +CEhb2, dcos_theta_di );
                            rvec_ScaledAdd( workspace->f_hb[j], +CEhb2, dcos_theta_dj );
                            rvec_ScaledAdd( workspace->f_hb[k], +CEhb2, dcos_theta_dk );
                            /* dr terms */
                            rvec_ScaledAdd( workspace->f_hb[j], -CEhb3 / r_jk, dvec_jk );
                            rvec_ScaledAdd( workspace->f_hb[k], +CEhb3 / r_jk, dvec_jk );
#endif
                        }
                    }
#if defined(QMMM)
                    }
#endif
                }
            }
        }

        if ( hblist != NULL )
        {
            sfree( hblist, __FILE__, __LINE__ );
        }
    }

    data->E_HB += e_hb_total;

#if defined(TEST_FORCES)
    fprintf( stderr, "Number of hydrogen bonds: %d\n", num_hb_intrs );
    fprintf( stderr, "Hydrogen Bond Energy: %g\n", data->E_HB );
#endif
}
