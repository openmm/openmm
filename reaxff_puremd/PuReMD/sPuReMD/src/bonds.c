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

#include "bonds.h"

#include "bond_orders.h"
#include "list.h"


void Bonds( reax_system *system, control_params *control,
        simulation_data *data, static_storage *workspace,
        reax_list **lists, output_controls *out_control )
{
    int i;
    real gp3, gp4, gp7, gp10, ebond_total;
    reax_list *bonds;

    bonds = lists[BONDS];
    gp3 = system->reax_param.gp.l[3];
    gp4 = system->reax_param.gp.l[4];
    gp7 = system->reax_param.gp.l[7];
    gp10 = system->reax_param.gp.l[10];
    ebond_total = 0.0;

#if defined(_OPENMP)
//    #pragma omp parallel default(shared) reduction(+: ebond_total)
#endif
    { 
        int j, pj;
        int start_i, end_i;
        int type_i, type_j;
        real ebond, pow_BOs_be2, exp_be12, CEbo;
        real exphu, exphua1, exphub1, exphuov, hulpov, estriph;
        real decobdbo, decobdboua, decobdboub;
        single_body_parameters *sbp_i, *sbp_j;
        two_body_parameters *twbp;
        bond_order_data *bo_ij;

#if defined(_OPENMP)
//        #pragma omp for schedule(guided)
#endif
        for ( i = 0; i < system->N; ++i )
        {
#if defined(QMMM)
            if ( system->atoms[i].qmmm_mask == TRUE )
            {
#endif
            start_i = Start_Index(i, bonds);
            end_i = End_Index(i, bonds);

            for ( pj = start_i; pj < end_i; ++pj )
            {
                if ( i < bonds->bond_list[pj].nbr )
                {
                    j = bonds->bond_list[pj].nbr;

#if defined(QMMM)
                    if ( system->atoms[j].qmmm_mask == TRUE )
                    {
#endif

                    type_i = system->atoms[i].type;
                    type_j = system->atoms[j].type;
                    sbp_i = &system->reax_param.sbp[type_i];
                    sbp_j = &system->reax_param.sbp[type_j];
                    twbp = &system->reax_param.tbp[type_i][type_j];
                    bo_ij = &bonds->bond_list[pj].bo_data;

                    pow_BOs_be2 = POW( bo_ij->BO_s, twbp->p_be2 );
                    exp_be12 = EXP( twbp->p_be1 * ( 1.0 - pow_BOs_be2 ) );
                    CEbo = -twbp->De_s * exp_be12
                        * (1.0 - twbp->p_be1 * twbp->p_be2 * pow_BOs_be2);

                    /* calculate bond energy */
                    ebond = -twbp->De_s * bo_ij->BO_s * exp_be12
                        - twbp->De_p * bo_ij->BO_pi
                        - twbp->De_pp * bo_ij->BO_pi2;
                    ebond_total += ebond;

                    /* calculate derivatives of bond orders */
                    bo_ij->Cdbo += CEbo;
                    bo_ij->Cdbopi -= CEbo + twbp->De_p;
                    bo_ij->Cdbopi2 -= CEbo + twbp->De_pp;

#if defined(TEST_ENERGY)
                    fprintf( out_control->ebond, "%6d%6d%24.15e%24.15e\n",
                             workspace->orig_id[i], workspace->orig_id[j],
                             // i+1, j+1,
                             bo_ij->BO, ebond );
#endif

#if defined(TEST_FORCES)
                    Add_dBO( system, lists, i, pj, CEbo, workspace->f_be );
                    Add_dBOpinpi2( system, lists, i, pj,
                            -(CEbo + twbp->De_p),
                            -(CEbo + twbp->De_pp),
                            workspace->f_be, workspace->f_be );
#endif

                    /* Stabilisation terminal triple bond in C-O */
                    if ( bo_ij->BO >= 1.00 )
                    {
                        if ( (strncmp( sbp_i->name, "C", sizeof(sbp_i->name) ) == 0
                                    && strncmp( sbp_j->name, "O", sizeof(sbp_j->name) ) == 0)
                                || (strncmp( sbp_i->name, "O", sizeof(sbp_i->name) ) == 0
                                    && strncmp( sbp_j->name, "C", sizeof(sbp_j->name) ) == 0) )
                        {
                            //ba = SQR( bo_ij->BO - 2.5 );
                            exphu = EXP( -gp7 * SQR(bo_ij->BO - 2.5) );
                            //oboa = abo(j1) - boa;
                            //obob = abo(j2) - boa;
                            exphua1 = EXP(-gp3 * (workspace->total_bond_order[i] - bo_ij->BO));
                            exphub1 = EXP(-gp3 * (workspace->total_bond_order[j] - bo_ij->BO));
                            //ovoab = abo(j1) - aval(it1) + abo(j2) - aval(it2);
                            exphuov = EXP(gp4 * (workspace->Delta[i] + workspace->Delta[j]));
                            hulpov = 1.0 / (1.0 + 25.0 * exphuov);

                            estriph = gp10 * exphu * hulpov * (exphua1 + exphub1);
                            //estrain(j1) = estrain(j1) + 0.5 * estriph;
                            //estrain(j2) = estrain(j2) + 0.5 * estriph;
                            ebond_total += estriph;

                            decobdbo = gp10 * exphu * hulpov * (exphua1 + exphub1)
                                * ( gp3 - 2.0 * gp7 * (bo_ij->BO - 2.5) );
                            decobdboua = -gp10 * exphu * hulpov
                                * (gp3 * exphua1 + 25.0 * gp4 * exphuov * hulpov * (exphua1 + exphub1));
                            decobdboub = -gp10 * exphu * hulpov
                                * (gp3 * exphub1 + 25.0 * gp4 * exphuov * hulpov * (exphua1 + exphub1));

                            bo_ij->Cdbo += decobdbo;
                            workspace->CdDelta[i] += decobdboua;
                            workspace->CdDelta[j] += decobdboub;

#if defined(TEST_ENERGY)
                            fprintf( out_control->ebond,
                                     "%6d%6d%24.15e%24.15e%24.15e%24.15e\n",
                                     workspace->orig_id[i], workspace->orig_id[j],
                                     //i + 1, j + 1,
                                     estriph, decobdbo, decobdboua, decobdboub );
#endif

#if defined(TEST_FORCES)
                            Add_dBO( system, lists, i, pj, decobdbo, workspace->f_be );
                            Add_dDelta( system, lists, i, decobdboua, workspace->f_be );
                            Add_dDelta( system, lists, j, decobdboub, workspace->f_be );
#endif
                        }
                    }
#if defined(QMMM)
                    }
#endif
                }
            }
#if defined(QMMM)
            }
#endif
        }
    }

    data->E_BE += ebond_total;
}
