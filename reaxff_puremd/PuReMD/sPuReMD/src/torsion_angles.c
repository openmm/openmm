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

#include "torsion_angles.h"

#include "bond_orders.h"
#include "box.h"
#include "list.h"
#include "lookup.h"
#include "vector.h"

#define MIN_SINE (1.0e-10)


static real Calculate_Omega( rvec dvec_ij, real r_ij, rvec dvec_jk, real r_jk,
        rvec dvec_kl, real r_kl, rvec dvec_li, real r_li,
        three_body_interaction_data *p_ijk,
        three_body_interaction_data *p_jkl,
        rvec dcos_omega_di, rvec dcos_omega_dj,
        rvec dcos_omega_dk, rvec dcos_omega_dl )
{
    real unnorm_cos_omega, unnorm_sin_omega, omega;
    real sin_ijk, cos_ijk, sin_jkl, cos_jkl;
    real htra, htrb, htrc, hthd, hthe, hnra, hnrc, hnhd, hnhe;
    real arg, poem, tel;
    rvec cross_jk_kl;

    assert( r_ij > 0.0 );
    assert( r_jk > 0.0 );
    assert( r_kl > 0.0 );
    assert( r_li > 0.0 );

    sin_ijk = SIN( p_ijk->theta );
    cos_ijk = COS( p_ijk->theta );
    sin_jkl = SIN( p_jkl->theta );
    cos_jkl = COS( p_jkl->theta );

    /* omega */
    unnorm_cos_omega = -rvec_Dot( dvec_ij, dvec_jk ) * rvec_Dot( dvec_jk, dvec_kl )
        + SQR( r_jk ) * rvec_Dot( dvec_ij, dvec_kl );
    rvec_Cross( cross_jk_kl, dvec_jk, dvec_kl );
    unnorm_sin_omega = -r_jk * rvec_Dot( dvec_ij, cross_jk_kl );
    omega = atan2( unnorm_sin_omega, unnorm_cos_omega );

    /* derivatives */
    /* coef for adjusments to cos_theta's */
    /* rla = r_ij, rlb = r_jk, rlc = r_kl, r4 = r_li;
     * coshd = cos_ijk, coshe = cos_jkl;
     * sinhd = sin_ijk, sinhe = sin_jkl; */
    htra = r_ij + cos_ijk * ( r_kl * cos_jkl - r_jk );
    htrb = r_jk - r_ij * cos_ijk - r_kl * cos_jkl;
    htrc = r_kl + cos_jkl * ( r_ij * cos_ijk - r_jk );
    hthd = r_ij * sin_ijk * ( r_jk - r_kl * cos_jkl );
    hthe = r_kl * sin_jkl * ( r_jk - r_ij * cos_ijk );
    hnra = r_kl * sin_ijk * sin_jkl;
    hnrc = r_ij * sin_ijk * sin_jkl;
    hnhd = r_ij * r_kl * cos_ijk * sin_jkl;
    hnhe = r_ij * r_kl * sin_ijk * cos_jkl;

    poem = 2.0 * r_ij * r_kl * sin_ijk * sin_jkl;
    if ( poem < 1.0e-20 )
    {
        poem = 1.0e-20;
    }

    tel = (SQR(r_ij) + SQR(r_jk) + SQR(r_kl) - SQR(r_li))
        - 2.0 * ( r_ij * r_jk * cos_ijk - r_ij * r_kl * cos_ijk * cos_jkl
                + r_jk * r_kl * cos_jkl );

    arg  = tel / poem;
    if ( arg >  1.0 )
    {
        arg =  1.0;
    }
    if ( arg < -1.0 )
    {
        arg = -1.0;
    }

    if ( sin_ijk >= 0.0 && sin_ijk <= MIN_SINE )
    {
        sin_ijk = MIN_SINE;
    }
    else if ( sin_ijk <= 0.0 && sin_ijk >= -MIN_SINE )
    {
        sin_ijk = -MIN_SINE;
    }
    if ( sin_jkl >= 0.0 && sin_jkl <= MIN_SINE )
    {
        sin_jkl = MIN_SINE;
    }
    else if ( sin_jkl <= 0.0 && sin_jkl >= -MIN_SINE )
    {
        sin_jkl = -MIN_SINE;
    }

    // dcos_omega_di
    rvec_ScaledSum( dcos_omega_di, (htra - arg * hnra) / r_ij, dvec_ij, -1.0, dvec_li );
    rvec_ScaledAdd( dcos_omega_di, -(hthd - arg * hnhd) / sin_ijk, p_ijk->dcos_dk );
    rvec_Scale( dcos_omega_di, 2.0 / poem, dcos_omega_di );

    // dcos_omega_dj
    rvec_ScaledSum( dcos_omega_dj, -(htra - arg * hnra) / r_ij, dvec_ij,
                    -htrb / r_jk, dvec_jk );
    rvec_ScaledAdd( dcos_omega_dj, -(hthd - arg * hnhd) / sin_ijk, p_ijk->dcos_dj );
    rvec_ScaledAdd( dcos_omega_dj, -(hthe - arg * hnhe) / sin_jkl, p_jkl->dcos_di );
    rvec_Scale( dcos_omega_dj, 2.0 / poem, dcos_omega_dj );

    // dcos_omega_dk
    rvec_ScaledSum( dcos_omega_dk, -(htrc - arg * hnrc) / r_kl, dvec_kl,
                    htrb / r_jk, dvec_jk );
    rvec_ScaledAdd( dcos_omega_dk, -(hthd - arg * hnhd) / sin_ijk, p_ijk->dcos_di );
    rvec_ScaledAdd( dcos_omega_dk, -(hthe - arg * hnhe) / sin_jkl, p_jkl->dcos_dj );
    rvec_Scale( dcos_omega_dk, 2.0 / poem, dcos_omega_dk );

    // dcos_omega_dl
    rvec_ScaledSum( dcos_omega_dl, (htrc - arg * hnrc) / r_kl, dvec_kl, 1.0, dvec_li );
    rvec_ScaledAdd( dcos_omega_dl, -(hthe - arg * hnhe) / sin_jkl, p_jkl->dcos_dk );
    rvec_Scale( dcos_omega_dl, 2.0 / poem, dcos_omega_dl );

    return omega;
}


void Torsion_Angles( reax_system *system, control_params *control,
        simulation_data *data, static_storage *workspace,
        reax_list **lists, output_controls *out_control )
{
#if defined(DEBUG_FOCUS)
    int num_frb_intrs;
#endif
    real p_tor2, p_tor3, p_tor4, p_cot2;
    reax_list *bonds, *thb_intrs;
    real e_tor_total, e_con_total;

    p_tor2 = system->reax_param.gp.l[23];
    p_tor3 = system->reax_param.gp.l[24];
    p_tor4 = system->reax_param.gp.l[25];
    p_cot2 = system->reax_param.gp.l[27];
    bonds = lists[BONDS];
    thb_intrs = lists[THREE_BODIES];
    e_tor_total = 0.0;
    e_con_total = 0.0;
#if defined(DEBUG_FOCUS)
    num_frb_intrs = 0;
#endif

#if defined(_OPENMP)
    #pragma omp parallel default(shared) reduction(+: e_tor_total, e_con_total)
#endif
    {
        int i, j, k, l, pi, pj, pk, pl, pij, plk;
        int type_i, type_j, type_k, type_l;
        int start_j, end_j;
        int start_pj, end_pj, start_pk, end_pk;
        real Delta_j, Delta_k;
        real r_ij, r_jk, r_kl, r_li;
        real BOA_ij, BOA_jk, BOA_kl;
        real exp_tor2_ij, exp_tor2_jk, exp_tor2_kl;
        real exp_tor1, exp_tor3_DjDk, exp_tor4_DjDk, exp_tor34_inv;
        real exp_cot2_jk, exp_cot2_ij, exp_cot2_kl;
        real fn10, f11_DjDk, dfn11, fn12;
        real theta_ijk, theta_jkl;
        real sin_ijk, sin_jkl;
        real cos_ijk, cos_jkl;
        real tan_ijk_i, tan_jkl_i;
        real omega, cos_omega, cos2omega, cos3omega;
        rvec dcos_omega_di, dcos_omega_dj, dcos_omega_dk, dcos_omega_dl;
        real CV, cmn, CEtors1, CEtors2, CEtors3, CEtors4;
        real CEtors5, CEtors6, CEtors7, CEtors8, CEtors9;
        real Cconj, CEconj1, CEconj2, CEconj3;
        real CEconj4, CEconj5, CEconj6;
        real e_tor, e_con;
        rvec dvec_li, force_i, force_j, force_k, force_l, x_i, x_j, x_k, x_l;
        rtensor press;
        ivec rel_box_li;
        //rtensor total_rtensor, temp_rtensor;
        four_body_header *fbh;
        four_body_parameters *fbp;
        bond_data *pbond_ij, *pbond_jk, *pbond_kl;
        bond_order_data *bo_ij, *bo_jk, *bo_kl;
        three_body_interaction_data *p_ijk, *p_jkl;
        rvec *f_i, *f_j, *f_k, *f_l;
#if defined(_OPENMP)
        int tid = omp_get_thread_num( );

        #pragma omp for schedule(static)
#endif
        for ( j = 0; j < system->N; ++j )
        {
#if defined(QMMM)
            if ( system->atoms[j].qmmm_mask == TRUE )
            {
#endif
            type_j = system->atoms[j].type;
            Delta_j = workspace->Delta_boc[j];
            start_j = Start_Index(j, bonds);
            end_j = End_Index(j, bonds);
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

            for ( pk = start_j; pk < end_j; ++pk )
            {
                pbond_jk = &bonds->bond_list[pk];
                k = pbond_jk->nbr;

#if defined(QMMM)
            if ( system->atoms[k].qmmm_mask == TRUE )
            {
#endif
                bo_jk = &pbond_jk->bo_data;
                BOA_jk = bo_jk->BO - control->thb_cut;
#if defined(_OPENMP)
                f_k = &workspace->f_local[tid * system->N + k];
#else
                f_k = &system->atoms[k].f;
#endif

                /* see if there are any 3-body interactions involving j&k
                 * where j is the central atom. Otherwise there is no point in
                 * trying to form a 4-body interaction out of this neighborhood */
                if ( j < k
                        && bo_jk->BO > control->thb_cut
                        && Num_Entries(pk, thb_intrs) > 0 )
                {
                    pj = pbond_jk->sym_index; // pj points to j on k's list

                    /* do the same check as above:
                     * are there any 3-body interactions
                     * involving k&j where k is the central atom */
                    if ( Num_Entries(pj, thb_intrs) > 0 )
                    {
                        type_k = system->atoms[k].type;
                        Delta_k = workspace->Delta_boc[k];
                        r_jk = pbond_jk->d;

                        start_pk = Start_Index( pk, thb_intrs );
                        end_pk = End_Index( pk, thb_intrs );
                        start_pj = Start_Index( pj, thb_intrs );
                        end_pj = End_Index( pj, thb_intrs );
                        if ( control->ensemble == sNPT || control->ensemble == iNPT
                                || control->ensemble == aNPT || control->compute_pressure == TRUE )
                        {
                            rvec_Sum( x_k, x_j, pbond_jk->dvec );
                        }

                        exp_tor2_jk = EXP( -p_tor2 * BOA_jk );
                        exp_cot2_jk = EXP( -p_cot2 * SQR(BOA_jk - 1.5) );
                        exp_tor3_DjDk = EXP( -p_tor3 * (Delta_j + Delta_k) );
                        exp_tor4_DjDk = EXP( p_tor4  * (Delta_j + Delta_k) );
                        exp_tor34_inv = 1.0 / (1.0 + exp_tor3_DjDk + exp_tor4_DjDk);
                        f11_DjDk = (2.0 + exp_tor3_DjDk) * exp_tor34_inv;

                        /* pick i up from j-k interaction where j is the center atom */
                        for ( pi = start_pk; pi < end_pk; ++pi )
                        {
                            p_ijk = &thb_intrs->three_body_list[pi];
                            pij = p_ijk->pthb; // pij is pointer to i on j's bond_list
                            pbond_ij = &bonds->bond_list[pij];
                            bo_ij = &pbond_ij->bo_data;

                            if ( bo_ij->BO > control->thb_cut )
                            {
                                i = p_ijk->thb;

#if defined(QMMM)
                                if ( system->atoms[i].qmmm_mask == TRUE )
                                {
#endif
                                type_i = system->atoms[i].type;
                                r_ij = pbond_ij->d;
                                BOA_ij = bo_ij->BO - control->thb_cut;

#if defined(_OPENMP)
                                f_i = &workspace->f_local[tid * system->N + i];
#else
                                f_i = &system->atoms[i].f;
#endif
                                if ( control->ensemble == sNPT || control->ensemble == iNPT
                                        || control->ensemble == aNPT || control->compute_pressure == TRUE )
                                {
                                    rvec_Sum( x_i, x_j, pbond_ij->dvec );
                                }

                                theta_ijk = p_ijk->theta;
                                sin_ijk = SIN( theta_ijk );
                                cos_ijk = COS( theta_ijk );
                                //tan_ijk_i = 1.0 / TAN( theta_ijk );
                                if ( sin_ijk >= 0.0 && sin_ijk <= MIN_SINE )
                                {
                                    tan_ijk_i = cos_ijk / MIN_SINE;
                                }
                                else if ( sin_ijk <= 0.0 && sin_ijk >= -MIN_SINE )
                                {
                                    tan_ijk_i = cos_ijk / -MIN_SINE;
                                }
                                else
                                {
                                    tan_ijk_i = cos_ijk / sin_ijk;
                                }

                                exp_tor2_ij = EXP( -p_tor2 * BOA_ij );
                                exp_cot2_ij = EXP( -p_cot2 * SQR(BOA_ij - 1.5) );

                                /* pick l up from j-k interaction where k is the central atom */
                                for ( pl = start_pj; pl < end_pj; ++pl )
                                {
                                    p_jkl = &thb_intrs->three_body_list[pl];
                                    l = p_jkl->thb;

#if defined(QMMM)
                                    if ( system->atoms[l].qmmm_mask == TRUE )
                                    {
#endif
                                    plk = p_jkl->pthb; //pointer to l on k's bond_list!
                                    pbond_kl = &bonds->bond_list[plk];
                                    bo_kl = &pbond_kl->bo_data;
                                    type_l = system->atoms[l].type;
                                    fbh = &system->reax_param.fbp[type_i][type_j][type_k][type_l];
                                    fbp = &system->reax_param.fbp[type_i][type_j]
                                            [type_k][type_l].prm[0];

                                    if ( i != l && fbh->cnt > 0
                                            && bo_kl->BO > control->thb_cut
                                            && bo_ij->BO * bo_jk->BO * bo_kl->BO > control->thb_cut )
                                    {
#if defined(_OPENMP)
                                        f_l = &workspace->f_local[tid * system->N + l];
#else
                                        f_l = &system->atoms[l].f;
#endif

#if defined(DEBUG_FOCUS)
#if defined(_OPENMP)
                                        #pragma omp atomic
#endif
                                        ++num_frb_intrs;
#endif

                                        r_kl = pbond_kl->d;
                                        BOA_kl = bo_kl->BO - control->thb_cut;

                                        theta_jkl = p_jkl->theta;
                                        sin_jkl = SIN( theta_jkl );
                                        cos_jkl = COS( theta_jkl );
                                        //tan_jkl_i = 1.0 / TAN( theta_jkl );
                                        if ( sin_jkl >= 0.0 && sin_jkl <= MIN_SINE )
                                        {
                                            tan_jkl_i = cos_jkl / MIN_SINE;
                                        }
                                        else if ( sin_jkl <= 0.0 && sin_jkl >= -MIN_SINE )
                                        {
                                            tan_jkl_i = cos_jkl / -MIN_SINE;
                                        }
                                        else
                                        {
                                            tan_jkl_i = cos_jkl / sin_jkl;
                                        }

                                        /* compute relative integer coordinates
                                         * for periodic image of simulation
                                         * box for i-l bond using the transitive
                                         * propertiy of the i-j, j-k, and k-l
                                         * bond relative coordinates
                                         * */
                                        ivec_Sum( rel_box_li, pbond_jk->rel_box, pbond_kl->rel_box );
                                        ivec_ScaledAdd( rel_box_li, -1, pbond_ij->rel_box );
                                        ivec_Scale( rel_box_li, -1, rel_box_li );

                                        r_li = control->compute_atom_distance( &system->box,
                                                system->atoms[l].x, system->atoms[i].x,
                                                system->atoms[l].rel_map, system->atoms[i].rel_map,
                                                rel_box_li, dvec_li );

                                        /* omega and its derivative */
                                        //cos_omega = Calculate_Omega( pbond_ij->dvec, r_ij, pbond_jk->dvec,
                                        omega = Calculate_Omega( pbond_ij->dvec, r_ij, pbond_jk->dvec,
                                                r_jk, pbond_kl->dvec, r_kl, dvec_li, r_li, p_ijk, p_jkl,
                                                dcos_omega_di, dcos_omega_dj, dcos_omega_dk, dcos_omega_dl );
                                        cos_omega = COS( omega );
                                        cos2omega = COS( 2.0 * omega );
                                        cos3omega = COS( 3.0 * omega );
                                        /* end omega calculations */

                                        /* torsion energy */
                                        exp_tor1 = EXP( fbp->p_tor1 * SQR(2.0 - bo_jk->BO_pi - f11_DjDk) );
                                        exp_tor2_kl = EXP( -p_tor2 * BOA_kl );
                                        exp_cot2_kl = EXP( -p_cot2 * SQR(BOA_kl - 1.5) );
                                        fn10 = (1.0 - exp_tor2_ij) * (1.0 - exp_tor2_jk)
                                            * (1.0 - exp_tor2_kl);

                                        CV = 0.5 * ( fbp->V1 * (1.0 + cos_omega)
                                                + fbp->V2 * exp_tor1 * (1.0 - cos2omega)
                                                + fbp->V3 * (1.0 + cos3omega) );
//                                        CV = 0.5 * fbp->V1 * (1.0 + cos_omega)
//                                            + fbp->V2 * exp_tor1 * (1.0 - SQR(cos_omega))
//                                            + fbp->V3 * (0.5 + 2.0 * CUBE(cos_omega) - 1.5 * cos_omega);

                                        e_tor = fn10 * sin_ijk * sin_jkl * CV;
                                        e_tor_total += e_tor;

                                        dfn11 = (-p_tor3 * exp_tor3_DjDk
                                                + (p_tor3 * exp_tor3_DjDk - p_tor4 * exp_tor4_DjDk)
                                                * (2.0 + exp_tor3_DjDk) * exp_tor34_inv) * exp_tor34_inv;

                                        CEtors1 = sin_ijk * sin_jkl * CV;

                                        CEtors2 = -fn10 * 2.0 * fbp->p_tor1 * fbp->V2 * exp_tor1 *
                                            (2.0 - bo_jk->BO_pi - f11_DjDk) * (1.0 - SQR(cos_omega)) *
                                            sin_ijk * sin_jkl;
                                        CEtors3 = CEtors2 * dfn11;

                                        CEtors4 = CEtors1 * p_tor2 * exp_tor2_ij *
                                            (1.0 - exp_tor2_jk) * (1.0 - exp_tor2_kl);
                                        CEtors5 = CEtors1 * p_tor2 * exp_tor2_jk *
                                            (1.0 - exp_tor2_ij) * (1.0 - exp_tor2_kl);
                                        CEtors6 = CEtors1 * p_tor2 * exp_tor2_kl *
                                            (1.0 - exp_tor2_ij) * (1.0 - exp_tor2_jk);

                                        cmn = -fn10 * CV;
                                        CEtors7 = cmn * sin_jkl * tan_ijk_i;
                                        CEtors8 = cmn * sin_ijk * tan_jkl_i;
                                        CEtors9 = fn10 * sin_ijk * sin_jkl
                                            * (0.5 * fbp->V1 - 2.0 * fbp->V2 * exp_tor1 * cos_omega
                                                    + 1.5 * fbp->V3 * (cos2omega + 2.0 * SQR(cos_omega)));
//                                        CEtors7 = cmn * sin_jkl * cos_ijk;
//                                        CEtors8 = cmn * sin_ijk * cos_jkl;
//                                        CEtors9 = fn10 * sin_ijk * sin_jkl
//                                            * (0.5 * fbp->V1 - 2.0 * fbp->V2 * exp_tor1 * cos_omega
//                                                    + fbp->V3 * (6.0 * SQR(cos_omega) - 1.5));
                                        /* end  of torsion energy */

                                        /* 4-body conjugation energy */
                                        fn12 = exp_cot2_ij * exp_cot2_jk * exp_cot2_kl;
                                        e_con = fbp->p_cot1 * fn12
                                            * (1.0 + (SQR(cos_omega) - 1.0) * sin_ijk * sin_jkl);
                                        e_con_total += e_con;

                                        Cconj = -2.0 * fn12 * fbp->p_cot1 * p_cot2 *
                                                (1.0 + (SQR(cos_omega) - 1.0) * sin_ijk * sin_jkl);

                                        CEconj1 = Cconj * (BOA_ij - 1.5);
                                        CEconj2 = Cconj * (BOA_jk - 1.5);
                                        CEconj3 = Cconj * (BOA_kl - 1.5);

                                        CEconj4 = -fbp->p_cot1 * fn12
                                            * (SQR(cos_omega) - 1.0) * sin_jkl * tan_ijk_i;
                                        CEconj5 = -fbp->p_cot1 * fn12
                                            * (SQR(cos_omega) - 1.0) * sin_ijk * tan_jkl_i;
//                                        CEconj4 = -fbp->p_cot1 * fn12
//                                            * (SQR(cos_omega) - 1.0) * sin_jkl * cos_ijk;
//                                        CEconj5 = -fbp->p_cot1 * fn12
//                                            * (SQR(cos_omega) - 1.0) * sin_ijk * cos_jkl;
                                        CEconj6 = 2.0 * fbp->p_cot1 * fn12
                                            * cos_omega * sin_ijk * sin_jkl;
                                        /* end 4-body conjugation energy */

                                        //fprintf( stderr, "%6d %6d %6d %6d %7.3f %7.3f %7.3f %7.3f ",
                                        //   workspace->orig_id[i], workspace->orig_id[j],
                                        //       workspace->orig_id[k], workspace->orig_id[l],
                                        //    omega, cos_omega, cos2omega, cos3omega );
//                                        fprintf( stderr,
//                                                "%4d%4d%4d%4d%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f\n",
//                                                workspace->orig_id[i], workspace->orig_id[j],
//                                                workspace->orig_id[k], workspace->orig_id[l],
//                                                CEtors2, CEtors3, CEtors4, CEtors5, CEtors6,
//                                                CEtors7, CEtors8, CEtors9 );
//                                        fflush( stderr );
//                                        fprintf( stderr,
//                                                "%4d%4d%4d%4d%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f\n",
//                                                workspace->orig_id[i], workspace->orig_id[j],
//                                                workspace->orig_id[k], workspace->orig_id[l],
//                                                CEconj1, CEconj2, CEconj3, CEconj4, CEconj5, CEconj6 );
//                                        fflush( stderr );
                                        //fprintf( stderr, "%12.8f%12.8f%12.8f%12.8f%12.8f%12.8f\n",
                                        //    theta_ijk, theta_jkl, sin_ijk,
                                        //    sin_jkl, cos_jkl, tan_jkl_i );

                                        /* forces */
#if defined(_OPENMP)
                                        #pragma omp atomic
#endif
                                        bo_jk->Cdbopi += CEtors2;
#if defined(_OPENMP)
                                        #pragma omp atomic
#endif
                                        workspace->CdDelta[j] += CEtors3;
#if defined(_OPENMP)
                                        #pragma omp atomic
#endif
                                        workspace->CdDelta[k] += CEtors3;
#if defined(_OPENMP)
                                        #pragma omp atomic
#endif
                                        bo_ij->Cdbo += (CEtors4 + CEconj1);
#if defined(_OPENMP)
                                        #pragma omp atomic
#endif
                                        bo_jk->Cdbo += (CEtors5 + CEconj2);
#if defined(_OPENMP)
                                        #pragma omp atomic
#endif
                                        bo_kl->Cdbo += (CEtors6 + CEconj3);

                                        if ( control->compute_pressure == FALSE &&
                                                (control->ensemble == NVE || control->ensemble == nhNVT
                                                 || control->ensemble == bNVT) )
                                        {
                                            /* dcos_theta_ijk */
                                            rvec_ScaledAdd( *f_i,
                                                    CEtors7 + CEconj4, p_ijk->dcos_dk );
                                            rvec_ScaledAdd( *f_j,
                                                    CEtors7 + CEconj4, p_ijk->dcos_dj );
                                            rvec_ScaledAdd( *f_k,
                                                    CEtors7 + CEconj4, p_ijk->dcos_di );

                                            /* dcos_theta_jkl */
                                            rvec_ScaledAdd( *f_j,
                                                    CEtors8 + CEconj5, p_jkl->dcos_di );
                                            rvec_ScaledAdd( *f_k,
                                                    CEtors8 + CEconj5, p_jkl->dcos_dj );
                                            rvec_ScaledAdd( *f_l,
                                                    CEtors8 + CEconj5, p_jkl->dcos_dk );

                                            /* dcos_omega */
                                            rvec_ScaledAdd( *f_i,
                                                    CEtors9 + CEconj6, dcos_omega_di );
                                            rvec_ScaledAdd( *f_j,
                                                    CEtors9 + CEconj6, dcos_omega_dj );
                                            rvec_ScaledAdd( *f_k,
                                                    CEtors9 + CEconj6, dcos_omega_dk );
                                            rvec_ScaledAdd( *f_l,
                                                    CEtors9 + CEconj6, dcos_omega_dl );
                                        }
                                        else if ( control->ensemble == sNPT || control->ensemble == iNPT
                                                || control->ensemble == aNPT || control->compute_pressure == TRUE )
                                        {
                                            /* dcos_theta_ijk */
                                            rvec_Scale( force_i, CEtors7 + CEconj4, p_ijk->dcos_dk );
                                            rvec_Scale( force_j, CEtors7 + CEconj4, p_ijk->dcos_dj );
                                            rvec_Scale( force_k, CEtors7 + CEconj4, p_ijk->dcos_di );

                                            /* dcos_theta_jkl */
                                            rvec_ScaledAdd( force_j, CEtors8 + CEconj5, p_jkl->dcos_di );
                                            rvec_ScaledAdd( force_k, CEtors8 + CEconj5, p_jkl->dcos_dj );
                                            rvec_Scale( force_l, CEtors8 + CEconj5, p_jkl->dcos_dk );

                                            /* dcos_omega */
                                            rvec_ScaledAdd( force_i, CEtors9 + CEconj6, dcos_omega_di );
                                            rvec_ScaledAdd( force_j, CEtors9 + CEconj6, dcos_omega_dj );
                                            rvec_ScaledAdd( force_k, CEtors9 + CEconj6, dcos_omega_dk );
                                            rvec_ScaledAdd( force_l, CEtors9 + CEconj6, dcos_omega_dl );

                                            rvec_Add( *f_i, force_i );
                                            rvec_Add( *f_j, force_j );
                                            rvec_Add( *f_k, force_k );
                                            rvec_Add( *f_l, force_l );

                                            /* pressure */
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

                                            rvec_OuterProduct( press, force_k, x_k );
#if !defined(_OPENMP)
                                            rtensor_Add( data->press, press );
#else
                                            rtensor_Add( data->press_local[tid], press );
#endif

                                            rvec_Sum( x_l, x_k, pbond_kl->dvec );
                                            rvec_OuterProduct( press, force_l, x_l );
#if !defined(_OPENMP)
                                            rtensor_Add( data->press, press );
#else
                                            rtensor_Add( data->press_local[tid], press );
#endif

                                            /* This part is intended for a fully-flexible box */
                                            /* rvec_ScaledSum( temp_rvec,
                                               CEtors7 + CEconj4, p_ijk->dcos_dk,      // i
                                               CEtors9 + CEconj6, dcos_omega_di );
                                               rvec_OuterProduct( temp_rtensor,
                                               temp_rvec, system->atoms[i].x );
                                               rtensor_Copy( total_rtensor, temp_rtensor );

                                               rvec_ScaledSum( temp_rvec,
                                               CEtors7 + CEconj4, p_ijk->dcos_dj,      // j
                                               CEtors8 + CEconj5, p_jkl->dcos_di );
                                               rvec_ScaledAdd( temp_rvec,
                                               CEtors9 + CEconj6, dcos_omega_dj );
                                               rvec_OuterProduct( temp_rtensor,
                                               temp_rvec, system->atoms[j].x );
                                               rtensor_Add( total_rtensor, temp_rtensor );

                                               rvec_ScaledSum( temp_rvec,
                                               CEtors7 + CEconj4, p_ijk->dcos_di,      // k
                                               CEtors8 + CEconj5, p_jkl->dcos_dj );
                                               rvec_ScaledAdd( temp_rvec,
                                               CEtors9 + CEconj6, dcos_omega_dk );
                                               rvec_OuterProduct( temp_rtensor,
                                               temp_rvec, system->atoms[k].x );
                                               rtensor_Add( total_rtensor, temp_rtensor );

                                               rvec_ScaledSum( temp_rvec,
                                               CEtors8 + CEconj5, p_jkl->dcos_dk,      // l
                                               CEtors9 + CEconj6, dcos_omega_dl );
                                               rvec_OuterProduct( temp_rtensor,
                                               temp_rvec, system->atoms[l].x );
                                               rtensor_Copy( total_rtensor, temp_rtensor );

                                               if( pbond_ij->imaginary || pbond_jk->imaginary ||
                                               pbond_kl->imaginary )
                                               rtensor_ScaledAdd( data->flex_bar.P, -1., total_rtensor );
                                               else
                                               rtensor_Add( data->flex_bar.P, total_rtensor ); */
                                        }

#if defined(TEST_ENERGY)
                                        /*fprintf( out_control->etor,
                                           //"%12.8f%12.8f%12.8f%12.8f%12.8f%12.8f%12.8f\n",
                                           //r_ij, r_jk, r_kl,
                                           "%12.8f%12.8f%12.8f%12.8f\n",
                                           cos_ijk, cos_jkl, sin_ijk, sin_jkl );*/
                                        // fprintf( out_control->etor, "%12.8f\n", dfn11 );
                                        fprintf( out_control->etor, "%12.8f%12.8f%12.8f\n",
                                                 fn10, cos_omega, CV );

                                        fprintf( out_control->etor,
                                                 "%12.8f%12.8f%12.8f%12.8f%12.8f%12.8f%12.8f%12.8f\n",
                                                 CEtors2, CEtors3, CEtors4, CEtors5,
                                                 CEtors6, CEtors7, CEtors8, CEtors9 );

                                        /* fprintf( out_control->etor,
                                           "%12.8f%12.8f%12.8f%12.8f%12.8f%12.8f%12.8f%12.8f\n",
                                           htra, htrb, htrc, hthd, hthe, hnra, hnrc, hnhd, hnhe ); */

                                        fprintf( out_control->etor,
                                                 "%12.8f%12.8f%12.8f%12.8f%12.8f%12.8f\n",
                                                 CEconj1, CEconj2, CEconj3, CEconj4, CEconj5, CEconj6 );
                                        /* fprintf(out_control->etor,"%23.15e%23.15e%23.15e%23.15e\n",
                                           fbp->V1, fbp->V2, fbp->V3, fbp->p_tor1 );*/

                                        fprintf( out_control->etor,
                                                 //"%6d%6d%6d%6d%23.15e%23.15e%23.15e%23.15e\n",
                                                 "%6d%6d%6d%6d%12.8f%12.8f\n",
                                                 workspace->orig_id[i], workspace->orig_id[j],
                                                 workspace->orig_id[k], workspace->orig_id[l],
                                                 e_tor, e_con );
                                        //RAD2DEG(omega), BOA_jk, e_tor, data->E_Tor );

                                        fprintf( out_control->econ,
                                                 "%6d%6d%6d%6d%23.15e%23.15e%23.15e%23.15e%23.15e%23.15e\n",
                                                 workspace->orig_id[i], workspace->orig_id[j],
                                                 workspace->orig_id[k], workspace->orig_id[l],
                                                 RAD2DEG(omega), BOA_ij, BOA_jk, BOA_kl,
                                                 e_con, data->E_Con );

                                        /* fprintf( out_control->etor,
                                           "%12.8f%12.8f%12.8f\n%12.8f%12.8f%12.8f\n%12.8f%12.8f%12.8f\n",
                                           (CEtors7 + CEconj4)*p_ijk->dcos_dk[0],
                                           (CEtors7 + CEconj4)*p_ijk->dcos_dk[1],
                                           (CEtors7 + CEconj4)*p_ijk->dcos_dk[2],
                                           (CEtors7 + CEconj4)*p_ijk->dcos_dj[0],
                                           (CEtors7 + CEconj4)*p_ijk->dcos_dj[1],
                                           (CEtors7 + CEconj4)*p_ijk->dcos_dj[2],
                                           (CEtors7 + CEconj4)*p_ijk->dcos_di[0],
                                           (CEtors7 + CEconj4)*p_ijk->dcos_di[1],
                                           (CEtors7 + CEconj4)*p_ijk->dcos_di[2] ); */


                                        /* fprintf( out_control->etor,
                                           "%12.8f%12.8f%12.8f\n%12.8f%12.8f%12.8f\n%12.8f%12.8f%12.8f\n",
                                           (CEtors8 + CEconj5)*p_jkl->dcos_di[0],
                                           (CEtors8 + CEconj5)*p_jkl->dcos_di[1],
                                           (CEtors8 + CEconj5)*p_jkl->dcos_di[2],
                                           (CEtors8 + CEconj5)*p_jkl->dcos_dj[0],
                                           (CEtors8 + CEconj5)*p_jkl->dcos_dj[1],
                                           (CEtors8 + CEconj5)*p_jkl->dcos_dj[2],
                                           (CEtors8 + CEconj5)*p_jkl->dcos_dk[0],
                                           (CEtors8 + CEconj5)*p_jkl->dcos_dk[1],
                                           (CEtors8 + CEconj5)*p_jkl->dcos_dk[2] ); */

                                        fprintf( out_control->etor,
                                                 "%12.8f%12.8f%12.8f\n%12.8f%12.8f%12.8f\n%12.8f%12.8f%12.8f\n%12.8f%12.8f%12.8f\n",
                                                 dcos_omega_di[0], dcos_omega_di[1], dcos_omega_di[2],
                                                 dcos_omega_dj[0], dcos_omega_dj[1], dcos_omega_dj[2],
                                                 dcos_omega_dk[0], dcos_omega_dk[1], dcos_omega_dk[2],
                                                 dcos_omega_dl[0], dcos_omega_dl[1], dcos_omega_dl[2] );
#endif

#if defined(TEST_FORCES)
                                        /* Torsion Forces */
                                        Add_dBOpinpi2( system, lists, j, pk, CEtors2, 0.,
                                                workspace->f_tor, workspace->f_tor );
                                        Add_dDelta( system, lists, j, CEtors3, workspace->f_tor );
                                        Add_dDelta( system, lists, k, CEtors3, workspace->f_tor );
                                        Add_dBO( system, lists, j, pij, CEtors4, workspace->f_tor );
                                        Add_dBO( system, lists, j, pk, CEtors5, workspace->f_tor );
                                        Add_dBO( system, lists, k, plk, CEtors6, workspace->f_tor );

                                        rvec_ScaledAdd( workspace->f_tor[i], CEtors7, p_ijk->dcos_dk );
                                        rvec_ScaledAdd( workspace->f_tor[j], CEtors7, p_ijk->dcos_dj );
                                        rvec_ScaledAdd( workspace->f_tor[k], CEtors7, p_ijk->dcos_di );

                                        rvec_ScaledAdd( workspace->f_tor[j], CEtors8, p_jkl->dcos_di );
                                        rvec_ScaledAdd( workspace->f_tor[k], CEtors8, p_jkl->dcos_dj );
                                        rvec_ScaledAdd( workspace->f_tor[l], CEtors8, p_jkl->dcos_dk );

                                        rvec_ScaledAdd( workspace->f_tor[i], CEtors9, dcos_omega_di );
                                        rvec_ScaledAdd( workspace->f_tor[j], CEtors9, dcos_omega_dj );
                                        rvec_ScaledAdd( workspace->f_tor[k], CEtors9, dcos_omega_dk );
                                        rvec_ScaledAdd( workspace->f_tor[l], CEtors9, dcos_omega_dl );

                                        /* Conjugation Forces */
                                        Add_dBO( system, lists, j, pij, CEconj1, workspace->f_con );
                                        Add_dBO( system, lists, j, pk, CEconj2, workspace->f_con );
                                        Add_dBO( system, lists, k, plk, CEconj3, workspace->f_con );

                                        rvec_ScaledAdd( workspace->f_con[i], CEconj4, p_ijk->dcos_dk );
                                        rvec_ScaledAdd( workspace->f_con[j], CEconj4, p_ijk->dcos_dj );
                                        rvec_ScaledAdd( workspace->f_con[k], CEconj4, p_ijk->dcos_di );

                                        rvec_ScaledAdd( workspace->f_con[j], CEconj5, p_jkl->dcos_di );
                                        rvec_ScaledAdd( workspace->f_con[k], CEconj5, p_jkl->dcos_dj );
                                        rvec_ScaledAdd( workspace->f_con[l], CEconj5, p_jkl->dcos_dk );

                                        rvec_ScaledAdd( workspace->f_con[i], CEconj6, dcos_omega_di );
                                        rvec_ScaledAdd( workspace->f_con[j], CEconj6, dcos_omega_dj );
                                        rvec_ScaledAdd( workspace->f_con[k], CEconj6, dcos_omega_dk );
                                        rvec_ScaledAdd( workspace->f_con[l], CEconj6, dcos_omega_dl );
#endif
                                    } // pl check ends
#if defined(QMMM)
                                    }
#endif
                                } // pl loop ends
#if defined(QMMM)
                                }
#endif
                            } // pi check ends
                        } // pi loop ends
                    } // k-j neighbor check ends
                } // j<k && j-k neighbor check ends
#if defined(QMMM)
                }
#endif
            } // pk loop ends
#if defined(QMMM)
            }
#endif
        } // j loop
    }

     data->E_Tor += e_tor_total;
     data->E_Con += e_con_total;

#if defined(DEBUG_FOCUS)
    fprintf( stderr, "[INFO] Torsion_Angles: num_frb_intrs = %d\n", num_frb_intrs );
    fprintf( stderr, "[INFO] Torsion_Angles: e_tor = %g, e_con = %g\n",
             data->E_Tor, data->E_Con );

//    fprintf( stderr, "[INFO] Torsion_Angles: press = (%23.15e %23.15e %23.15e)\n",
//            data->press[0][0], data->press[1][1], data->press[2][2] );
#endif
}
