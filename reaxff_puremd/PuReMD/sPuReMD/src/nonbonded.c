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

#include "nonbonded.h"

#include "bond_orders.h"
#include "list.h"
#include "lookup.h"
#include "vector.h"


static void Compute_Polarization_Energy( reax_system* system,
        control_params *control, simulation_data* data, static_storage *workspace )
{
    int i, type_i;
    real e_pol, q;

    /* Compute Polarization Energy */
    e_pol = 0.0;

    if ( control->charge_method == QEQ_CM
            || control->charge_method == EE_CM )
    {
#if defined(_OPENMP)
        #pragma omp parallel for default(none) private(q, type_i) shared(system) \
            reduction(+: e_pol) schedule(static)
#endif
        for ( i = 0; i < system->N; i++ )
        {
#if defined(QMMM)
            if ( system->atoms[i].qmmm_mask == TRUE )
            {
#endif
            q = system->atoms[i].q;
            type_i = system->atoms[i].type;

            e_pol += KCALpMOL_to_EV * (system->reax_param.sbp[ type_i ].chi * q
                    + (system->reax_param.sbp[ type_i ].eta / 2.0) * SQR( q ));
#if defined(QMMM)
            }
#endif
        }
    }
    else if ( control->charge_method == ACKS2_CM )
    {
#if defined(_OPENMP)
        #pragma omp parallel for default(none) private(q, type_i) shared(system, workspace) \
            reduction(+: e_pol) schedule(static)
#endif
        for ( i = 0; i < system->N; i++ )
        {
#if defined(QMMM)
            if ( system->atoms[i].qmmm_mask == TRUE )
            {
#endif
            q = system->atoms[i].q;
            type_i = system->atoms[i].type;

            /* energy due to first and second order EE parameters */
            e_pol += KCALpMOL_to_EV * (system->reax_param.sbp[ type_i ].chi * q
                    + (system->reax_param.sbp[ type_i ].eta / 2.0) * SQR( q ));

            /* energy due to coupling with kinetic energy potential */
            e_pol += KCALpMOL_to_EV * system->atoms[i].q * workspace->s[0][ system->N + i ];
#if defined(QMMM)
            }
#endif
        }
    }

    data->E_Pol = e_pol;
}


void vdW_Coulomb_Energy( reax_system *system, control_params *control,
        simulation_data *data, static_storage *workspace,
        reax_list **lists, output_controls *out_control )
{
    int i;
    real p_vdW1, p_vdW1i;
    reax_list *far_nbrs;
    real e_vdW_total, e_ele_total;

    p_vdW1 = system->reax_param.gp.l[28];
    p_vdW1i = 1.0 / p_vdW1;
    far_nbrs = lists[FAR_NBRS];
    e_vdW_total = 0.0;
    e_ele_total = 0.0;

#if defined(_OPENMP)
    #pragma omp parallel default(shared) reduction(+: e_vdW_total, e_ele_total)
#endif
    {
        int j, pj;
        int start_i, end_i;
        real self_coef;
        real powr_vdW1, powgi_vdW1;
        real r_ij, fn13, exp1, exp2, e_base, de_base;
        real Tap, dTap, dfn13, CEvd, CEclmb;
        real dr3gamij_1, dr3gamij_3;
        real e_ele, e_vdW, e_core, de_core, e_clb, de_clb;
        real d, xcut, bond_softness, d_bond_softness, effpot_diff;
        rvec force, x_i, x_j;
        rtensor press;
        //rtensor temp_rtensor, total_rtensor;
        two_body_parameters *twbp;
        far_neighbor_data *nbr_pj;
#if defined(_OPENMP)
        int tid;

        tid = omp_get_thread_num( );
#endif

        e_ele = 0.0;
        e_vdW = 0.0;

#if defined(_OPENMP)
        #pragma omp for schedule(guided)
#endif
        for ( i = 0; i < system->N; ++i )
        {
            start_i = Start_Index( i, far_nbrs );
            end_i = End_Index( i, far_nbrs );
            if ( control->ensemble == sNPT || control->ensemble == iNPT
                    || control->ensemble == aNPT || control->compute_pressure == TRUE )
            {
                rvec_iMultiply( x_i, system->atoms[i].rel_map, system->box.box_norms );
                rvec_Add( x_i, system->atoms[i].x );
            }

            for ( pj = start_i; pj < end_i; ++pj )
            {
                if ( far_nbrs->far_nbr_list[pj].d <= control->nonb_cut )
                {
                    nbr_pj = &far_nbrs->far_nbr_list[pj];
                    j = nbr_pj->nbr;

                    r_ij = nbr_pj->d;
                    twbp = &system->reax_param.tbp[ system->atoms[i].type ]
                             [ system->atoms[j].type ];
                    /* i == j: self-interaction from periodic image,
                     * important for supporting small boxes! */
                    self_coef = (i == j) ? 0.5 : 1.0;

                    /* Calculate Taper and its derivative */
                    Tap = workspace->Tap[7] * r_ij
                        + workspace->Tap[6];
                    Tap = Tap * r_ij + workspace->Tap[5];
                    Tap = Tap * r_ij + workspace->Tap[4];
                    Tap = Tap * r_ij + workspace->Tap[3];
                    Tap = Tap * r_ij + workspace->Tap[2];
                    Tap = Tap * r_ij + workspace->Tap[1];
                    Tap = Tap * r_ij + workspace->Tap[0];

                    dTap = 7.0 * workspace->Tap[7] * r_ij
                        + 6.0 * workspace->Tap[6];
                    dTap = dTap * r_ij + 5.0 * workspace->Tap[5];
                    dTap = dTap * r_ij + 4.0 * workspace->Tap[4];
                    dTap = dTap * r_ij + 3.0 * workspace->Tap[3];
                    dTap = dTap * r_ij + 2.0 * workspace->Tap[2];
                    dTap = dTap * r_ij + workspace->Tap[1];

#if defined(QMMM)
                    if ( system->atoms[i].qmmm_mask == TRUE
                            && system->atoms[j].qmmm_mask == TRUE )
                    {
#endif
                    /* van der Waals calculations */
                    if ( system->atoms[i].is_dummy == FALSE
                            && system->atoms[j].is_dummy == FALSE )
                    {
                        if ( system->reax_param.gp.vdw_type == 1
                                || system->reax_param.gp.vdw_type == 3 )
                        {
                            /* shielding */
                            powr_vdW1 = POW( r_ij, p_vdW1 );
                            //TODO: better to compute and cache these values at simulation start rather than computing on-the-fly
                            powgi_vdW1 = POW( 1.0 / twbp->gamma_w, p_vdW1 );

                            fn13 = POW( powr_vdW1 + powgi_vdW1, p_vdW1i );
                            exp1 = EXP( twbp->alpha * (1.0 - fn13 / twbp->r_vdW) );
                            exp2 = EXP( 0.5 * twbp->alpha * (1.0 - fn13 / twbp->r_vdW) );
                            e_base = twbp->D * (exp1 - 2.0 * exp2);

                            e_vdW = self_coef * (e_base * Tap);
                            e_vdW_total += e_vdW;

                            dfn13 = POW( r_ij, p_vdW1 - 1.0 )
                                * POW( powr_vdW1 + powgi_vdW1, p_vdW1i - 1.0 );
                            de_base = (twbp->D * twbp->alpha / twbp->r_vdW) * (exp2 - exp1) * dfn13;
                        }
                        /* no shielding */
                        else
                        {
                            exp1 = EXP( twbp->alpha * (1.0 - r_ij / twbp->r_vdW) );
                            exp2 = EXP( 0.5 * twbp->alpha * (1.0 - r_ij / twbp->r_vdW) );
                            e_base = twbp->D * (exp1 - 2.0 * exp2);

                            e_vdW = self_coef * (e_base * Tap);
                            e_vdW_total += e_vdW;

                            de_base = (twbp->D * twbp->alpha / twbp->r_vdW) * (exp2 - exp1);
                        }

                        /* calculate inner core repulsion */
                        if ( system->reax_param.gp.vdw_type == 2 || system->reax_param.gp.vdw_type == 3 )
                        {
                            e_core = twbp->ecore * EXP( twbp->acore * (1.0 - (r_ij / twbp->rcore)) );
                            e_vdW += self_coef * (e_core * Tap);
                            e_vdW_total += self_coef * (e_core * Tap);

                            de_core = -(twbp->acore / twbp->rcore) * e_core;
                        }
                        else
                        {
                            e_core = 0.0;
                            de_core = 0.0;
                        }

                        CEvd = self_coef * ( (de_base + de_core) * Tap
                                + (e_base + e_core) * dTap );
                    }
                    else
                    {
                        e_core = 0.0;
                        de_core = 0.0;
                        CEvd = 0.0;
                    }
#if defined(QMMM)
                    }
                    else
                    {
                        e_core = 0.0;
                        de_core = 0.0;
                        CEvd = 0.0;
                    }
#endif

#if defined(DEBUG_FOCUS)
                    fprintf( stderr, "%6d%6d%24.12f%24.12f%24.12f%24.12f\n",
                            i + 1, j + 1, 
                            e_base, de_base, e_core, de_core ); fflush( stderr );
#endif

                    /* Coulomb Calculations */
#if defined(QMMM)
                    if ( system->atoms[i].qmmm_mask == TRUE || system->atoms[j].qmmm_mask == TRUE )
                    {
#endif
                    dr3gamij_1 = r_ij * r_ij * r_ij + POW( twbp->gamma, -3.0 );
                    dr3gamij_3 = POW( dr3gamij_1 , 1.0 / 3.0 );
                    e_clb = C_ELE * (system->atoms[i].q * system->atoms[j].q) / dr3gamij_3;
                    e_ele = self_coef * (e_clb * Tap);
                    e_ele_total += e_ele;

                    de_clb = -C_ELE * (system->atoms[i].q * system->atoms[j].q)
                            * (r_ij * r_ij) / POW( dr3gamij_1, 4.0 / 3.0);
                    CEclmb = self_coef * (de_clb * Tap + e_clb * dTap);
#if defined(QMMM)
                    }
                    else
                    {
                        de_clb = 0.0;
                        CEclmb = 0.0;
                    }
#endif

#if defined(DEBUG_FOCUS)
                    fprintf( stderr, "%6d%6d%24.12f%24.12f\n",
                            i + 1, j + 1, e_clb, de_clb ); fflush( stderr );
#endif

                    if ( control->compute_pressure == FALSE &&
                            (control->ensemble == NVE || control->ensemble == nhNVT
                             || control->ensemble == bNVT) )
                    {
#if !defined(_OPENMP)
                        rvec_ScaledAdd( system->atoms[i].f,
                                -1.0 * (CEvd + CEclmb) / r_ij, nbr_pj->dvec );
                        rvec_ScaledAdd( system->atoms[j].f,
                                (CEvd + CEclmb) / r_ij, nbr_pj->dvec );
#else
                        rvec_ScaledAdd( workspace->f_local[tid * system->N + i],
                                -1.0 * (CEvd + CEclmb) / r_ij, nbr_pj->dvec );
                        rvec_ScaledAdd( workspace->f_local[tid * system->N + j],
                                (CEvd + CEclmb) / r_ij, nbr_pj->dvec );
#endif

#if defined(DEBUG_FOCUS)
                        fprintf( stderr, "%6d%6d%24.12f%24.12f%24.12f%24.12f%24.12f%24.12f%24.12f\n",
                                i + 1, j + 1, (CEvd + CEclmb) / r_ij, 
                                nbr_pj->dvec[0],
                                nbr_pj->dvec[1],
                                nbr_pj->dvec[2],
                                (CEvd + CEclmb) / r_ij * nbr_pj->dvec[0],
                                (CEvd + CEclmb) / r_ij * nbr_pj->dvec[1],
                                (CEvd + CEclmb) / r_ij * nbr_pj->dvec[2] ); fflush( stderr );
#endif
                    }
                    else if ( control->ensemble == sNPT || control->ensemble == iNPT
                            || control->ensemble == aNPT || control->compute_pressure == TRUE )
                    {
                        /* for pressure coupling, terms not related to bond order
                         * derivatives are added directly into pressure vector/tensor */
                        rvec_Scale( force, -1.0 * (CEvd + CEclmb) / r_ij, nbr_pj->dvec );
#if !defined(_OPENMP)
                        rvec_Add( system->atoms[i].f, force );
                        rvec_ScaledAdd( system->atoms[j].f, -1.0, force );
#else
                        rvec_Add( workspace->f_local[tid * system->N + i], force );
                        rvec_ScaledAdd( workspace->f_local[tid * system->N + j], -1.0, force );
#endif

                        /* pressure */
                        rvec_OuterProduct( press, force, x_i );
#if !defined(_OPENMP)
                        rtensor_Add( data->press, press );
#else
                        rtensor_Add( data->press_local[tid], press );
#endif

                        rvec_Sum( x_j, x_i, nbr_pj->dvec );
                        rvec_Scale( force, -1.0, force );
                        rvec_OuterProduct( press, force, x_j );
#if !defined(_OPENMP)
                        rtensor_Add( data->press, press );
#else
                        rtensor_Add( data->press_local[tid], press );
#endif

#if defined(DEBUG_FOCUS)
                        fprintf( stderr, "nonbonded(%d,%d): rel_box (%d %d %d)",
                                i, j, nbr_pj->rel_box[0],
                                nbr_pj->rel_box[1], nbr_pj->rel_box[2] );

                        fprintf( stderr, "force(%f %f %f)", force[0], force[1], force[2] );

                        fprintf( stderr, "press (%12.6f %12.6f %12.6f)\n",
                                data->press[0][0], data->press[1][1],
                                data->press[2][2] );
#endif

                        /* This part is intended for a fully-flexible box */
//                        rvec_OuterProduct( temp_rtensor, nbr_pj->dvec, system->atoms[i].x );
//                        rtensor_Scale( total_rtensor, F_C * -(CEvd + CEclmb), temp_rtensor );
//                        rvec_OuterProduct( temp_rtensor, nbr_pj->dvec, system->atoms[j].x );
//                        rtensor_ScaledAdd( total_rtensor, F_C * +(CEvd + CEclmb), temp_rtensor );

                        /* This is an external force due to an imaginary nbr */
//                        if ( nbr_pj->imaginary )
//                            rtensor_ScaledAdd( data->flex_bar.P, -1.0, total_rtensor );
                        /* This interaction is completely internal */
//                        else
//                            rtensor_Add( data->flex_bar.P, total_rtensor );
                    }

#if defined(TEST_ENERGY)
                    rvec_MakeZero( force );
                    rvec_ScaledAdd( force, CEvd, nbr_pj->dvec );
                    fprintf( out_control->evdw,
                             "%6d%6d%24.15e%24.15e%24.15e%24.15e%24.15e\n",
                             //i+1, j+1,
                             MIN( workspace->orig_id[i], workspace->orig_id[j] ),
                             MAX( workspace->orig_id[i], workspace->orig_id[j] ),
                             r_ij, e_vdW, force[0], force[1], force[2]/*, e_vdW_total*/ );

                    fprintf( out_control->ecou, "%6d%6d%24.15e%24.15e%24.15e%24.15e\n",
                             MIN( workspace->orig_id[i], workspace->orig_id[j] ),
                             MAX( workspace->orig_id[i], workspace->orig_id[j] ),
                             r_ij, system->atoms[i].q, system->atoms[j].q,
                             e_ele/*, e_ele_total*/ );
#endif

#if defined(TEST_FORCES)
                    rvec_ScaledAdd( workspace->f_vdw[i], -CEvd, nbr_pj->dvec );
                    rvec_ScaledAdd( workspace->f_vdw[j], +CEvd, nbr_pj->dvec );
                    rvec_ScaledAdd( workspace->f_ele[i], -CEclmb, nbr_pj->dvec );
                    rvec_ScaledAdd( workspace->f_ele[j], +CEclmb, nbr_pj->dvec );
#endif
                }
            }
        }

        //TODO: better integration of ACKS2 code below with above code for performance
#if defined(_OPENMP)
        #pragma omp barrier
#endif

        /* contribution to energy and gradients (atoms and cell)
         * due to geometry-dependent terms in the ACKS2
         * kinetic energy */
        if ( control->charge_method == ACKS2_CM )
        {
#if defined(_OPENMP)
            #pragma omp for schedule(guided)
#endif
            for ( i = 0; i < system->N; ++i )
            {
                if ( control->ensemble == sNPT || control->ensemble == iNPT
                        || control->ensemble == aNPT || control->compute_pressure == TRUE )
                {
                    rvec_iMultiply( x_i, system->atoms[i].rel_map, system->box.box_norms );
                    rvec_Add( x_i, system->atoms[i].x );
                }

                for ( pj = Start_Index(i, far_nbrs); pj < End_Index(i, far_nbrs); ++pj )
                {
                    nbr_pj = &far_nbrs->far_nbr_list[pj];
                    j = nbr_pj->nbr;

                    /* kinetic energy terms */
                    xcut = 0.5 * ( system->reax_param.sbp[ system->atoms[i].type ].b_s_acks2
                            + system->reax_param.sbp[ system->atoms[j].type ].b_s_acks2 );

                    if ( far_nbrs->far_nbr_list[pj].d < xcut )
                    {
                        d = far_nbrs->far_nbr_list[pj].d / xcut;
                        bond_softness = system->reax_param.gp.l[34] * POW( d, 3.0 )
                            * POW( 1.0 - d, 6.0 );

                        if ( bond_softness > 0.0 )
                        {
                            /* Coulombic energy contribution */
                            effpot_diff = workspace->s[0][system->N + i]
                                - workspace->s[0][system->N + j];
                            e_ele = -0.5 * KCALpMOL_to_EV * bond_softness
                                * SQR( effpot_diff );
                            e_ele_total += e_ele;

                            /* forces contribution */
                            d_bond_softness = system->reax_param.gp.l[34]
                                * 3.0 / xcut * POW( d, 2.0 )
                                * POW( 1.0 - d, 5.0 ) * (1.0 - 3.0 * d);
                            d_bond_softness = -0.5 * d_bond_softness
                                * SQR( effpot_diff );
                            d_bond_softness = KCALpMOL_to_EV * d_bond_softness
                                / far_nbrs->far_nbr_list[pj].d;

#if defined(DEBUG_FOCUS)
                            fprintf( stderr, "%6d%6d%12.5f%12.5f%12.5f\n",
                                    j + 1, i + 1,
                                    d_bond_softness * nbr_pj->dvec[0],
                                    d_bond_softness * nbr_pj->dvec[1],
                                    d_bond_softness * nbr_pj->dvec[2] ); fflush( stderr );
#endif

                            if ( control->compute_pressure == FALSE &&
                                    (control->ensemble == NVE || control->ensemble == nhNVT
                                     || control->ensemble == bNVT) )
                            {
#if !defined(_OPENMP)
                                rvec_ScaledAdd( system->atoms[i].f,
                                        -1.0 * d_bond_softness, nbr_pj->dvec  );
                                rvec_ScaledAdd( system->atoms[j].f,
                                        d_bond_softness, nbr_pj->dvec );
#else
                                rvec_ScaledAdd( workspace->f_local[tid * system->N + i],
                                        -1.0 * d_bond_softness, nbr_pj->dvec  );
                                rvec_ScaledAdd( workspace->f_local[tid * system->N + j],
                                        d_bond_softness, nbr_pj->dvec );
#endif
                            }
                            else if ( control->ensemble == sNPT || control->ensemble == iNPT
                                    || control->ensemble == aNPT || control->compute_pressure == TRUE )
                            {
                                rvec_Scale( force, -1.0 * d_bond_softness, nbr_pj->dvec );
#if !defined(_OPENMP)
                                rvec_Add( system->atoms[i].f, force );
                                rvec_ScaledAdd( system->atoms[j].f,
                                        -1.0, force );
#else
                                rvec_Add( workspace->f_local[tid * system->N + i],
                                        force );
                                rvec_ScaledAdd( workspace->f_local[tid * system->N + j],
                                        -1.0, force );
#endif

                                /* pressure */
                                rvec_OuterProduct( press, force, x_i );
#if !defined(_OPENMP)
                                rtensor_Add( data->press, press );
#else
                                rtensor_Add( data->press_local[tid], press );
#endif

                                rvec_Sum( x_j, x_i, nbr_pj->dvec );
                                rvec_Scale( force, -1.0, force );
                                rvec_OuterProduct( press, force, x_j );
#if !defined(_OPENMP)
                                rtensor_Add( data->press, press );
#else
                                rtensor_Add( data->press_local[tid], press );
#endif
                            }
                        }
                    }
                }
            }
        }
    }

    data->E_vdW = e_vdW_total;
    data->E_Ele = e_ele_total;

#if defined(DEBUG_FOCUS)
    fprintf( stderr, "[INFO] vdW_Coulomb_Energy: press = (%24.15e %24.15e %24.15e)\n",
            data->press[0][0], data->press[1][1], data->press[2][2] );
#endif

    if ( control->polarization_energy_enabled == TRUE )
    {
        Compute_Polarization_Energy( system, control, data, workspace );
    }
}


void Tabulated_vdW_Coulomb_Energy( reax_system *system, control_params *control,
        simulation_data *data, static_storage *workspace, reax_list **lists,
        output_controls *out_control )
{
    int steps, update_freq, update_energies;
    reax_list *far_nbrs;
    real e_vdW_total, e_ele_total;

    far_nbrs = lists[FAR_NBRS];
    steps = data->step - data->prev_steps;
    update_freq = out_control->log_update_freq;
    update_energies = update_freq > 0 && steps % update_freq == 0;
    e_vdW_total = 0.0;
    e_ele_total = 0.0;

#if defined(_OPENMP)
    #pragma omp parallel default(shared) reduction(+: e_vdW_total, e_ele_total)
#endif
    {
        int i, j, pj, r;
        int type_i, type_j, tmin, tmax;
        int start_i, end_i;
        real r_ij, self_coef, base, dif;
        real e_vdW, e_ele;
        real d, xcut, bond_softness, d_bond_softness, effpot_diff;
        real CEvd, CEclmb;
        rvec force;
        rtensor press;
        far_neighbor_data *nbr_pj;
        LR_lookup_table *t;
#if defined(_OPENMP)
        int tid;

        tid = omp_get_thread_num( );

        #pragma omp for schedule(guided)
#endif
        for ( i = 0; i < system->N; ++i )
        {
            type_i = system->atoms[i].type;
            start_i = Start_Index(i, far_nbrs);
            end_i = End_Index(i, far_nbrs);

            for ( pj = start_i; pj < end_i; ++pj )
            {
                if ( far_nbrs->far_nbr_list[pj].d <= control->nonb_cut )
                {
                    nbr_pj = &far_nbrs->far_nbr_list[pj];
                    j = nbr_pj->nbr;
                    type_j = system->atoms[j].type;
                    r_ij = nbr_pj->d;
                    self_coef = (i == j) ? 0.5 : 1.0;
                    tmin = MIN( type_i, type_j );
                    tmax = MAX( type_i, type_j );
                    t = &workspace->LR[tmin][tmax];

                    /* Cubic Spline Interpolation */
                    r = (int) (r_ij * t->inv_dx);
                    if ( r == 0 )
                    {
                        ++r;
                    }
                    base = (real) (r + 1) * t->dx;
                    dif = r_ij - base;

#if defined(DEBUG_FOCUS)
                    fprintf( stderr, "[INFO] r: %f, i: %d, base: %f, dif: %f\n", r, i, base, dif );
#endif

                    if ( update_energies )
                    {
                        e_vdW = ((t->vdW[r].d * dif + t->vdW[r].c) * dif + t->vdW[r].b)
                            * dif + t->vdW[r].a;
                        e_vdW *= self_coef;

                        e_ele = ((t->ele[r].d * dif + t->ele[r].c) * dif + t->ele[r].b)
                            * dif + t->ele[r].a;
                        e_ele *= self_coef * system->atoms[i].q * system->atoms[j].q;

                        e_vdW_total += e_vdW;
                        e_ele_total += e_ele;
                    }

                    CEvd = ((t->CEvd[r].d * dif + t->CEvd[r].c) * dif + t->CEvd[r].b)
                        * dif + t->CEvd[r].a;
                    CEvd *= self_coef;

                    CEclmb = ((t->CEclmb[r].d * dif + t->CEclmb[r].c) * dif + t->CEclmb[r].b)
                        * dif + t->CEclmb[r].a;
                    CEclmb *= self_coef * system->atoms[i].q * system->atoms[j].q;

                    if ( control->compute_pressure == FALSE &&
                            (control->ensemble == NVE || control->ensemble == nhNVT
                             || control->ensemble == bNVT) )
                    {
#if !defined(_OPENMP)
                        rvec_ScaledAdd( system->atoms[i].f,
                                -(CEvd + CEclmb) / r_ij, nbr_pj->dvec );
                        rvec_ScaledAdd( system->atoms[j].f,
                                (CEvd + CEclmb) / r_ij, nbr_pj->dvec );
#else
                        rvec_ScaledAdd( workspace->f_local[tid * system->N + i],
                                -(CEvd + CEclmb) / r_ij, nbr_pj->dvec );
                        rvec_ScaledAdd( workspace->f_local[tid * system->N + j],
                                (CEvd + CEclmb) / r_ij, nbr_pj->dvec );
#endif
                    }
                    else if ( control->ensemble == sNPT || control->ensemble == iNPT
                            || control->ensemble == aNPT || control->compute_pressure == TRUE )
                    {
                        /* for pressure coupling, terms not related to bond order
                           derivatives are added directly into pressure vector/tensor */
                        rvec_Scale( force, (CEvd + CEclmb) / r_ij, nbr_pj->dvec );
#if !defined(_OPENMP)
                        rvec_ScaledAdd( system->atoms[i].f, -1.0, force );
                        rvec_Add( system->atoms[j].f, force );
#else
                        rvec_ScaledAdd( workspace->f_local[tid * system->N + i], -1.0, force );
                        rvec_Add( workspace->f_local[tid * system->N + j], force );
#endif

                        /* pressure */
                        rvec_Scale( force, -1.0, force );
                        rvec_OuterProduct( press, force, nbr_pj->dvec );
#if !defined(_OPENMP)
                        rtensor_Add( data->press, press );
#else
                        rtensor_Add( data->press_local[tid], press );
#endif
                    }

#if defined(TEST_ENERGY)
                    fprintf( out_control->evdw, "%6d%6d%24.15e%24.15e%24.15e\n",
                            workspace->orig_id[i], workspace->orig_id[j],
                            r_ij, e_vdW, data->E_vdW );
                    fprintf( out_control->ecou, "%6d%6d%24.15e%24.15e%24.15e%24.15e%24.15e\n",
                            workspace->orig_id[i], workspace->orig_id[j],
                            r_ij, system->atoms[i].q, system->atoms[j].q,
                            e_ele, data->E_Ele );
#endif

#if defined(TEST_FORCES)
                    rvec_ScaledAdd( workspace->f_vdw[i], -CEvd, nbr_pj->dvec );
                    rvec_ScaledAdd( workspace->f_vdw[j], +CEvd, nbr_pj->dvec );
                    rvec_ScaledAdd( workspace->f_ele[i], -CEclmb, nbr_pj->dvec );
                    rvec_ScaledAdd( workspace->f_ele[j], +CEclmb, nbr_pj->dvec );
#endif
                }
            }
        }

        //TODO: better integration AND tabulation of ACKS2 code below with above code for performance
#if defined(_OPENMP)
        #pragma omp barrier
#endif

        /* contribution to energy and gradients (atoms and cell)
         * due to geometry-dependent terms in the ACKS2
         * kinetic energy */
        if ( control->charge_method == ACKS2_CM )
        {
#if defined(_OPENMP)
            #pragma omp for schedule(guided)
#endif
            for ( i = 0; i < system->N; ++i )
            {
                for ( pj = Start_Index(i, far_nbrs); pj < End_Index(i, far_nbrs); ++pj )
                {
                    nbr_pj = &far_nbrs->far_nbr_list[pj];
                    j = nbr_pj->nbr;

                    /* kinetic energy terms */
                    xcut = 0.5 * ( system->reax_param.sbp[ system->atoms[i].type ].b_s_acks2
                            + system->reax_param.sbp[ system->atoms[j].type ].b_s_acks2 );

                    if ( far_nbrs->far_nbr_list[pj].d < xcut )
                    {
                        d = far_nbrs->far_nbr_list[pj].d / xcut;
                        bond_softness = system->reax_param.gp.l[34] * POW( d, 3.0 )
                            * POW( 1.0 - d, 6.0 );

                        if ( bond_softness > 0.0 )
                        {
                            /* Coulombic energy contribution */
                            effpot_diff = workspace->s[0][system->N + i]
                                - workspace->s[0][system->N + j];
                            e_ele = -0.5 * KCALpMOL_to_EV * bond_softness
                                * SQR( effpot_diff );
                            e_ele_total += e_ele;

                            /* forces contribution */
                            d_bond_softness = system->reax_param.gp.l[34]
                                * 3.0 / xcut * POW( d, 2.0 )
                                * POW( 1.0 - d, 5.0 ) * (1.0 - 3.0 * d);
                            d_bond_softness = -0.5 * d_bond_softness
                                * SQR( effpot_diff );
                            d_bond_softness = KCALpMOL_to_EV * d_bond_softness
                                / far_nbrs->far_nbr_list[pj].d;

#if defined(DEBUG_FOCUS)
                            fprintf( stderr, "%6d%6d%12.5f%12.5f%12.5f\n",
                                    j + 1, i + 1,
                                    d_bond_softness * nbr_pj->dvec[0],
                                    d_bond_softness * nbr_pj->dvec[1],
                                    d_bond_softness * nbr_pj->dvec[2] ); fflush( stderr );
#endif

#if !defined(_OPENMP)
                            rvec_ScaledAdd( system->atoms[i].f,
                                    -d_bond_softness, nbr_pj->dvec );
                            rvec_ScaledAdd( system->atoms[j].f,
                                    d_bond_softness, nbr_pj->dvec );
#else
                            rvec_ScaledAdd( workspace->f_local[tid * system->N + i],
                                    -d_bond_softness, nbr_pj->dvec );
                            rvec_ScaledAdd( workspace->f_local[tid * system->N + j],
                                    d_bond_softness, nbr_pj->dvec );
#endif

#if defined(TEST_FORCES)
                            rvec_ScaledAdd( workspace->f_ele[i],
                                    -d_bond_softness, nbr_pj->dvec );
                            rvec_ScaledAdd( workspace->f_ele[j],
                                    d_bond_softness, nbr_pj->dvec );
#endif
                        }
                    }
                }
            }
        }
    }

    data->E_vdW += e_vdW_total;
    data->E_Ele += e_ele_total;

    if ( control->polarization_energy_enabled == TRUE )
    {
        Compute_Polarization_Energy( system, control, data, workspace );
    }
}


void LR_vdW_Coulomb( reax_system *system, control_params *control,
        static_storage *workspace, int i, int j, real r_ij, LR_data *lr )
{
    real p_vdW1, p_vdW1i;
    real powr_vdW1, powgi_vdW1;
    real tmp, fn13, exp1, exp2, e_base, de_base;
    real Tap, dTap, dfn13;
    real dr3gamij_1, dr3gamij_3;
    real e_core, de_core;
    two_body_parameters *twbp;

    p_vdW1 = system->reax_param.gp.l[28];
    p_vdW1i = 1.0 / p_vdW1;
    twbp = &system->reax_param.tbp[i][j];

    /* Calculate Taper and its derivative */
    Tap = workspace->Tap[7] * r_ij
        + workspace->Tap[6];
    Tap = Tap * r_ij + workspace->Tap[5];
    Tap = Tap * r_ij + workspace->Tap[4];
    Tap = Tap * r_ij + workspace->Tap[3];
    Tap = Tap * r_ij + workspace->Tap[2];
    Tap = Tap * r_ij + workspace->Tap[1];
    Tap = Tap * r_ij + workspace->Tap[0];

    dTap = 7.0 * workspace->Tap[7] * r_ij
        + 6.0 * workspace->Tap[6];
    dTap = dTap * r_ij + 5.0 * workspace->Tap[5];
    dTap = dTap * r_ij + 4.0 * workspace->Tap[4];
    dTap = dTap * r_ij + 3.0 * workspace->Tap[3];
    dTap = dTap * r_ij + 2.0 * workspace->Tap[2];
    dTap = dTap * r_ij + workspace->Tap[1];

    /* vdWaals Calculations */
    if ( system->reax_param.gp.vdw_type == 1 || system->reax_param.gp.vdw_type == 3 )
    {
        /* shielding */
        powr_vdW1 = POW( r_ij, p_vdW1 );
        //TODO: better to compute and cache these values at simulation start rather than computing on-the-fly
        powgi_vdW1 = POW( 1.0 / twbp->gamma_w, p_vdW1 );

        fn13 = POW( powr_vdW1 + powgi_vdW1, p_vdW1i );
        exp1 = EXP( twbp->alpha * (1.0 - fn13 / twbp->r_vdW) );
        exp2 = EXP( 0.5 * twbp->alpha * (1.0 - fn13 / twbp->r_vdW) );
        e_base = twbp->D * (exp1 - 2.0 * exp2);

        lr->e_vdW = e_base * Tap;

        dfn13 = POW( r_ij, p_vdW1 - 1.0 )
            * POW( powr_vdW1 + powgi_vdW1, p_vdW1i - 1.0 );
        de_base = (twbp->D * twbp->alpha / twbp->r_vdW) * (exp2 - exp1) * dfn13;
    }
    /* no shielding */
    else
    {
        exp1 = EXP( twbp->alpha * (1.0 - r_ij / twbp->r_vdW) );
        exp2 = EXP( 0.5 * twbp->alpha * (1.0 - r_ij / twbp->r_vdW) );
        e_base = twbp->D * (exp1 - 2.0 * exp2);

        lr->e_vdW = e_base * Tap;

        de_base = (twbp->D * twbp->alpha / twbp->r_vdW) * (exp2 - exp1);
    }

    /* calculate inner core repulsion */
    if ( system->reax_param.gp.vdw_type == 2 || system->reax_param.gp.vdw_type == 3 )
    {
        e_core = twbp->ecore * EXP( twbp->acore * (1.0 - (r_ij / twbp->rcore)) );
        lr->e_vdW += e_core * Tap;

        de_core = -(twbp->acore / twbp->rcore) * e_core;
    }
    else
    {
        e_core = 0.0;
        de_core = 0.0;
    }

    lr->CEvd = (de_base + de_core) * Tap
            + (e_base + e_core) * dTap;

    /* Coulomb calculations */
    dr3gamij_1 = r_ij * r_ij * r_ij
        + POW( twbp->gamma, -3.0 );
    dr3gamij_3 = POW( dr3gamij_1, 1.0 / 3.0 );

    tmp = Tap / dr3gamij_3;
    lr->H = EV_to_KCALpMOL * tmp;
    lr->e_ele = C_ELE * tmp;

    /* fprintf( stderr,"i:%d(%d), j:%d(%d), gamma:%f,\
       Tap:%f, dr3gamij_3:%f, qi: %f, qj: %f\n",
       i, system->atoms[i].type, j, system->atoms[j].type,
       twbp->gamma, Tap, dr3gamij_3,
       system->atoms[i].q, system->atoms[j].q ); */

    lr->CEclmb = C_ELE * (dTap - Tap * r_ij / dr3gamij_1) / dr3gamij_3;

    /* fprintf( stdout, "%d %d\t%g\t%g  %g\t%g  %g\t%g  %g\n",
       i+1, j+1, r_ij, e_vdW, CEvd * r_ij,
       system->atoms[i].q, system->atoms[j].q, e_ele, CEclmb * r_ij ); */

    /* fprintf( stderr,"LR_Lookup:%3d%3d%5.3f-%8.5f,%8.5f%8.5f,%8.5f%8.5f\n",
       i, j, r_ij, lr->H, lr->e_vdW, lr->CEvd, lr->e_ele, lr->CEclmb ); */
}
