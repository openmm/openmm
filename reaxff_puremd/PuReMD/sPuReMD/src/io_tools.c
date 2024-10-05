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

#include "io_tools.h"

#include "list.h"
#include "geo_tools.h"
#include "system_props.h"
#include "tool_box.h"
#include "vector.h"


#if defined(TEST_FORCES)
void Print_Bond_Orders( reax_system *system, control_params *control,
        simulation_data *data, static_storage *workspace,
        reax_list **lists, output_controls *out_control )
{
    int i, pj, pk;
    bond_order_data *bo_ij;
    reax_list *bonds = lists[BONDS];
    reax_list *dBOs = lists[DBO];
    dbond_data *dbo_k;

    /* bond orders */
    fprintf( out_control->fbo, "%6s%6s%12s%12s%12s%12s%12s\n",
             "atom1", "atom2", "r_ij", "total_bo", "bo_s", "bo_p", "bo_pp" );

    for ( i = 0; i < system->N; ++i )
    {
        for ( pj = Start_Index(i, bonds); pj < End_Index(i, bonds); ++pj )
        {
            bo_ij = &bonds->bond_list[pj].bo_data;
            fprintf( out_control->fbo, "%6d%6d%23.15e%23.15e%23.15e%23.15e%23.15e\n",
                     //workspace->orig_id[i],
                     //workspace->orig_id[bonds->bond_list[pj].nbr],
                     i + 1,
                     bonds->bond_list[pj].nbr + 1,
                     bonds->bond_list[pj].d,
                     bo_ij->BO, bo_ij->BO_s, bo_ij->BO_pi, bo_ij->BO_pi2 );
        }
    }

    /* derivatives of bond orders */
    /* fprintf( out_control->fbo, "%6s%6s%10s%10s%10s%10s\n",
       "atom1", "atom2", "total_bo", "bo_s", "bo_p", "bo_pp"\n ); */

    for ( i = 0; i < system->N; ++i )
    {
        for ( pj = Start_Index(i, bonds); pj < End_Index(i, bonds); ++pj )
        {
            /*fprintf( out_control->fdbo, "%6d %6d\tstart: %6d\tend: %6d\n",
            workspace->orig_id[i],
            workspace->orig_id[bonds->bond_list[pj].nbr],
            Start_Index( pj, dBOs ), End_Index( pj, dBOs ) );*/

            for ( pk = Start_Index(pj, dBOs); pk < End_Index(pj, dBOs); ++pk )
            {
                dbo_k = &dBOs->dbo_list[pk];

                //if( !rvec_isZero( dbo_k->dBO ) )
                fprintf( out_control->fdbo, "%6d%6d%6d%23.15e%23.15e%23.15e\n",
                         workspace->orig_id[i],
                         workspace->orig_id[bonds->bond_list[pj].nbr],
                         workspace->orig_id[dbo_k->wrt],
                         dbo_k->dBO[0], dbo_k->dBO[1], dbo_k->dBO[2] );

                fprintf( out_control->fdbo, "%6d%6d%6d%23.15e%23.15e%23.15e\n",
                         workspace->orig_id[i],
                         workspace->orig_id[bonds->bond_list[pj].nbr],
                         workspace->orig_id[dbo_k->wrt],
                         dbo_k->dBOpi[0], dbo_k->dBOpi[1], dbo_k->dBOpi[2] );

                fprintf( out_control->fdbo, "%6d%6d%6d%23.15e%23.15e%23.15e\n",
                         workspace->orig_id[i],
                         workspace->orig_id[bonds->bond_list[pj].nbr],
                         workspace->orig_id[dbo_k->wrt],
                         dbo_k->dBOpi2[0], dbo_k->dBOpi2[1], dbo_k->dBOpi2[2] );
            }
        }
    }

    fflush( out_control->fdbo );
}


void Print_Bond_Forces( reax_system *system, control_params *control,
        simulation_data *data, static_storage *workspace,
        reax_list **lists, output_controls *out_control )
{
    int i;

    fprintf( out_control->fbond, "%d\n", data->step );
    fprintf( out_control->fbond, "%6s\t%s\n", "atom", "fbond" );

    for ( i = 0; i < system->N; ++i )
    {
        fprintf( out_control->fbond, "%6d %23.15e%23.15e%23.15e\n",
                workspace->orig_id[i], workspace->f_be[i][0],
                workspace->f_be[i][1], workspace->f_be[i][2]);
    }
}


void Print_LonePair_Forces( reax_system *system, control_params *control,
        simulation_data *data, static_storage *workspace,
        reax_list **lists, output_controls *out_control )
{
    int i;

    fprintf( out_control->flp, "%d\n", data->step );
    fprintf( out_control->flp, "%6s\t%s\n", "atom", "f_lonepair" );

    for ( i = 0; i < system->N; ++i )
    {
        fprintf(out_control->flp, "%6d %23.15e%23.15e%23.15e\n",
                workspace->orig_id[i],
                workspace->f_lp[i][0], workspace->f_lp[i][1], workspace->f_lp[i][2]);
    }

    fflush( out_control->flp );
}


void Print_OverUnderCoor_Forces( reax_system *system, control_params *control,
        simulation_data *data, static_storage *workspace, reax_list **lists,
        output_controls *out_control )
{
    int i;

    fprintf( out_control->fatom, "%d\n", data->step );
    fprintf( out_control->fatom, "%6s\t%-38s%-38s%-38s\n",
             "atom", "f_atom", "f_over", "f_under" );

    for ( i = 0; i < system->N; ++i )
    {
        if ( rvec_isZero( workspace->f_un[i] ) )
        {
            fprintf( out_control->fatom,
                     "%6d %23.15e%23.15e%23.15e 0 0 0\n",
                     workspace->orig_id[i], workspace->f_ov[i][0],
                     workspace->f_ov[i][1], workspace->f_ov[i][2] );
        }
        else
        {
            fprintf( out_control->fatom,
                     "%6d %23.15e%23.15e%23.15e %23.15e%23.15e%23.15e"\
                     "%23.15e%23.15e%23.15e\n",
                     workspace->orig_id[i],
                     workspace->f_un[i][0] + workspace->f_ov[i][0],
                     workspace->f_un[i][1] + workspace->f_ov[i][1],
                     workspace->f_un[i][2] + workspace->f_ov[i][2],
                     workspace->f_ov[i][0], workspace->f_ov[i][1],
                     workspace->f_ov[i][2],
                     workspace->f_un[i][0], workspace->f_un[i][1],
                     workspace->f_un[i][2] );
        }
    }

    fflush( out_control->fatom );
}


void Print_Three_Body_Forces( reax_system *system, control_params *control,
        simulation_data *data, static_storage *workspace,
        reax_list **lists, output_controls *out_control )
{
    int j;

    fprintf( out_control->f3body, "%d\n", data->step );
    fprintf( out_control->f3body, "%6s%-37s%-37s%-37s%-38s\n",
             "atom", "3-body total", "f_ang", "f_pen", "f_coa" );

    for ( j = 0; j < system->N; ++j )
    {
        if ( rvec_isZero(workspace->f_pen[j]) && rvec_isZero(workspace->f_coa[j]) )
        {
            fprintf( out_control->f3body, "%6d %23.15e%23.15e%23.15e  0 0 0  0 0 0\n",
                     workspace->orig_id[j], workspace->f_ang[j][0],
                     workspace->f_ang[j][1], workspace->f_ang[j][2] );
        }
        else if ( rvec_isZero(workspace->f_coa[j]) )
        {
            fprintf( out_control->f3body,
                     "%6d %23.15e%23.15e%23.15e %23.15e%23.15e%23.15e "\
                     "%23.15e%23.15e%23.15e\n",
                     workspace->orig_id[j],
                     workspace->f_ang[j][0] + workspace->f_pen[j][0],
                     workspace->f_ang[j][1] + workspace->f_pen[j][1],
                     workspace->f_ang[j][2] + workspace->f_pen[j][2],
                     workspace->f_ang[j][0], workspace->f_ang[j][1],
                     workspace->f_ang[j][2],
                     workspace->f_pen[j][0], workspace->f_pen[j][1],
                     workspace->f_pen[j][2] );
        }
        else
        {
            fprintf( out_control->f3body, "%6d %23.15e%23.15e%23.15e ",
                     workspace->orig_id[j],
                     workspace->f_ang[j][0] + workspace->f_pen[j][0] +
                     workspace->f_coa[j][0],
                     workspace->f_ang[j][1] + workspace->f_pen[j][1] +
                     workspace->f_coa[j][1],
                     workspace->f_ang[j][2] + workspace->f_pen[j][2] +
                     workspace->f_coa[j][2] );

            fprintf( out_control->f3body,
                     "%23.15e%23.15e%23.15e %23.15e%23.15e%23.15e "\
                     "%23.15e%23.15e%23.15e\n",
                     workspace->f_ang[j][0], workspace->f_ang[j][1],
                     workspace->f_ang[j][2],
                     workspace->f_pen[j][0], workspace->f_pen[j][1],
                     workspace->f_pen[j][2],
                     workspace->f_coa[j][0], workspace->f_coa[j][1],
                     workspace->f_coa[j][2] );
        }
    }

    fflush( out_control->f3body );
}


void Print_Hydrogen_Bond_Forces( reax_system *system, control_params *control,
        simulation_data *data, static_storage *workspace, reax_list **lists,
        output_controls *out_control )
{
    int j;

    fprintf( out_control->fhb, "%d\n", data->step );
    fprintf( out_control->fhb, "%6s\t%-38s\n", "atom", "f_hb" );

    for ( j = 0; j < system->N; ++j )
    {
        fprintf(out_control->fhb, "%6d\t[%23.15e%23.15e%23.15e]\n",
                workspace->orig_id[j],
                workspace->f_hb[j][0], workspace->f_hb[j][1], workspace->f_hb[j][2]);
    }

    fflush(out_control->fhb);
}


void Print_Four_Body_Forces( reax_system *system, control_params *control,
        simulation_data *data, static_storage *workspace,
        reax_list **lists, output_controls *out_control )
{
    int j;

    fprintf( out_control->f4body, "%6s\t%-38s%-38s%-38s\n",
             "atom", "4-body total", "f_tor", "f_con" );

    for ( j = 0; j < system->N; ++j )
    {
        if ( !rvec_isZero( workspace->f_con[j] ) )
        {
            fprintf( out_control->f4body,
                     "%6d %23.15e%23.15e%23.15e %23.15e%23.15e%23.15e "\
                     "%23.15e%23.15e%23.15e\n",
                     workspace->orig_id[j],
                     workspace->f_tor[j][0] + workspace->f_con[j][0],
                     workspace->f_tor[j][1] + workspace->f_con[j][1],
                     workspace->f_tor[j][2] + workspace->f_con[j][2],
                     workspace->f_tor[j][0], workspace->f_tor[j][1],
                     workspace->f_tor[j][2],
                     workspace->f_con[j][0], workspace->f_con[j][1],
                     workspace->f_con[j][2] );
        }
        else
        {
            fprintf( out_control->f4body,
                     "%6d %23.15e%23.15e%23.15e  0 0 0\n",
                     workspace->orig_id[j], workspace->f_tor[j][0],
                     workspace->f_tor[j][1], workspace->f_tor[j][2] );
        }
    }

    fflush( out_control->f4body );
}


void Print_vdW_Coulomb_Forces( reax_system *system, control_params *control,
        simulation_data *data, static_storage *workspace,
        reax_list **lists, output_controls *out_control )
{
    int  i;

    fprintf( out_control->fnonb, "%d\n", data->step );
    fprintf( out_control->fnonb, "%6s\t%-38s%-38s%-38s\n",
             "atom", "nonbonded total", "f_vdw", "f_ele" );

    for ( i = 0; i < system->N; ++i )
    {
        if ( !rvec_isZero(workspace->f_ele[i]) )
        {
            fprintf(out_control->fnonb,
                    "%6d %23.15e%23.15e%23.15e %23.15e%23.15e%23.15e "\
                    "%23.15e%23.15e%23.15e\n",
                    workspace->orig_id[i],
                    workspace->f_vdw[i][0] + workspace->f_ele[i][0],
                    workspace->f_vdw[i][1] + workspace->f_ele[i][1],
                    workspace->f_vdw[i][2] + workspace->f_ele[i][2],
                    workspace->f_vdw[i][0], workspace->f_vdw[i][1],
                    workspace->f_vdw[i][2],
                    workspace->f_ele[i][0], workspace->f_ele[i][1],
                    workspace->f_ele[i][2] );
        }
        else
        {
            fprintf(out_control->fnonb,
                    "%6d %23.15e%23.15e%23.15e  0 0 0\n",
                    workspace->orig_id[i], workspace->f_vdw[i][0],
                    workspace->f_vdw[i][1], workspace->f_vdw[i][2] );
        }
    }

    fflush( out_control->fnonb );
}


void Compare_Total_Forces( reax_system *system, control_params *control,
        simulation_data *data, static_storage *workspace,
        reax_list **lists, output_controls *out_control )
{
    int i;

    fprintf( out_control->ftot2, "%d\n", data->step );
    fprintf( out_control->ftot2, "%6s\t%-38s%-38s\n",
             "atom", "f_total", "test_force total" );

    for ( i = 0; i < system->N; ++i )
    {
        fprintf( out_control->ftot2,
                 "%6d %23.15e%23.15e%23.15e vs %23.15e%23.15e%23.15e\n",
                 workspace->orig_id[i],
                 system->atoms[i].f[0], system->atoms[i].f[1], system->atoms[i].f[2],
                 workspace->f_be[i][0] +
                 workspace->f_lp[i][0] + workspace->f_ov[i][0] + workspace->f_un[i][0] +
                 workspace->f_ang[i][0] + workspace->f_pen[i][0] +
                 workspace->f_coa[i][0] + workspace->f_hb[i][0] +
                 workspace->f_tor[i][0] + workspace->f_con[i][0] +
                 workspace->f_vdw[i][0] + workspace->f_ele[i][0],
                 workspace->f_be[i][1] +
                 workspace->f_lp[i][1] + workspace->f_ov[i][1] + workspace->f_un[i][1] +
                 workspace->f_ang[i][1] + workspace->f_pen[i][1] +
                 workspace->f_coa[i][1] + workspace->f_hb[i][1] +
                 workspace->f_tor[i][1] + workspace->f_con[i][1] +
                 workspace->f_vdw[i][1] + workspace->f_ele[i][1],
                 workspace->f_be[i][2] +
                 workspace->f_lp[i][2] + workspace->f_ov[i][2] + workspace->f_un[i][2] +
                 workspace->f_ang[i][2] + workspace->f_pen[i][2] +
                 workspace->f_coa[i][2] + workspace->f_hb[i][2] +
                 workspace->f_tor[i][2] + workspace->f_con[i][2] +
                 workspace->f_vdw[i][2] + workspace->f_ele[i][2] );
    }

    fflush( out_control->ftot2 );
}


void Init_Force_Test_Functions( control_params *control )
{
    control->print_intr_funcs[0] = &Print_Bond_Orders;
    control->print_intr_funcs[1] = &Print_Bond_Forces;
    control->print_intr_funcs[2] = &Print_LonePair_Forces;
    control->print_intr_funcs[3] = &Print_OverUnderCoor_Forces;
    control->print_intr_funcs[4] = &Print_Three_Body_Forces;
    control->print_intr_funcs[5] = &Print_Four_Body_Forces;
    control->print_intr_funcs[6] = &Print_Hydrogen_Bond_Forces;
    control->print_intr_funcs[7] = NULL;
    control->print_intr_funcs[8] = NULL;
    control->print_intr_funcs[9] = NULL;
}
#endif


/* near nbrs contain both i-j, j-i nbrhood info */
void Print_Near_Neighbors( reax_system *system, control_params *control,
        static_storage *workspace, reax_list **lists )
{
    int i, j, id_i, id_j;
    char fname[MAX_STR];
    FILE *fout;
    reax_list *near_nbrs = lists[NEAR_NBRS];

    snprintf( fname, MAX_STR, "%.*s.near_nbrs", MAX_STR - 11, control->sim_name );
    fout = sfopen( fname, "w", __FILE__, __LINE__ );

    for ( i = 0; i < system->N; ++i )
    {
        id_i = workspace->orig_id[i];

        for ( j = Start_Index(i, near_nbrs); j < End_Index(i, near_nbrs); ++j )
        {
            id_j = workspace->orig_id[near_nbrs->near_nbr_list[j].nbr];

            // if( id_i < id_j )
            fprintf( fout, "%6d%6d%23.15e%23.15e%23.15e%23.15e\n",
                     id_i, id_j,
                     near_nbrs->near_nbr_list[j].d,
                     near_nbrs->near_nbr_list[j].dvec[0],
                     near_nbrs->near_nbr_list[j].dvec[1],
                     near_nbrs->near_nbr_list[j].dvec[2] );
        }
    }

    sfclose( fout, __FILE__, __LINE__ );
}


/* near nbrs contain both i-j, j-i nbrhood info */
void Print_Near_Neighbors2( reax_system *system, control_params *control,
        static_storage *workspace, reax_list **lists )
{
    int i, j, id_i, id_j;
    char fname[MAX_STR];
    FILE *fout;
    reax_list *near_nbrs = lists[NEAR_NBRS];

    snprintf( fname, MAX_STR, "%.*s.near_nbrs_lgj", MAX_STR - 15, control->sim_name );
    fout = sfopen( fname, "w", __FILE__, __LINE__ );

    for ( i = 0; i < system->N; ++i )
    {
        id_i = workspace->orig_id[i];
        fprintf( fout, "%6d:", id_i);
        for ( j = Start_Index(i, near_nbrs); j < End_Index(i, near_nbrs); ++j )
        {
            id_j = workspace->orig_id[near_nbrs->near_nbr_list[j].nbr];
            fprintf( fout, "%6d", id_j);

            /* fprintf( fout, "%6d%6d%23.15e%23.15e%23.15e%23.15e\n",
            id_i, id_j,
             near_nbrs->near_nbr_list[j].d,
             near_nbrs->near_nbr_list[j].dvec[0],
             near_nbrs->near_nbr_list[j].dvec[1],
             near_nbrs->near_nbr_list[j].dvec[2] ); */
        }
        fprintf( fout, "\n");
    }

    sfclose( fout, __FILE__, __LINE__ );
}


/* far nbrs contain only i-j nbrhood info, no j-i. */
void Print_Far_Neighbors( reax_system const * const system,
        control_params const * const control,
        simulation_data const * const data,
        static_storage const * const workspace, reax_list ** const lists )
{
    int i, pj, id_i, id_j;
    char fname[MAX_STR];
    FILE *fout;
    reax_list *far_nbrs;

    far_nbrs = lists[FAR_NBRS];

    snprintf( fname, MAX_STR, "%.*s.%010d.far_nbrs", MAX_STR - 21, control->sim_name, data->step );
    fout = sfopen( fname, "w", __FILE__, __LINE__ );

    for ( i = 0; i < system->N; ++i )
    {
        id_i = workspace->orig_id[i];

        for ( pj = Start_Index(i, far_nbrs); pj < End_Index(i, far_nbrs); ++pj )
        {
            id_j = workspace->orig_id[far_nbrs->far_nbr_list[pj].nbr];

            fprintf( fout, "%6d%6d%3d%3d%3d %12.7f\n",
                     id_j, id_i,
                     far_nbrs->far_nbr_list[pj].rel_box[0],
                     far_nbrs->far_nbr_list[pj].rel_box[1],
                     far_nbrs->far_nbr_list[pj].rel_box[2],
                     far_nbrs->far_nbr_list[pj].d );

//            fprintf( fout, "%6d%6d%23.15e%23.15e%23.15e%23.15e\n",
//                     id_i, id_j,
//                     far_nbrs->far_nbr_list[pj].d,
//                     far_nbrs->far_nbr_list[pj].dvec[0],
//                     far_nbrs->far_nbr_list[pj].dvec[1],
//                     far_nbrs->far_nbr_list[pj].dvec[2] );
//
//            fprintf( fout, "%6d%6d%23.15e%23.15e%23.15e%23.15e\n",
//                     id_j, id_i,
//                     far_nbrs->far_nbr_list[pj].d,
//                     -far_nbrs->far_nbr_list[pj].dvec[0],
//                     -far_nbrs->far_nbr_list[pj].dvec[1],
//                     -far_nbrs->far_nbr_list[pj].dvec[2] );
        }
    }

    sfclose( fout, __FILE__, __LINE__ );
}


static int fn_qsort_intcmp( const void *a, const void *b )
{
    return ( *(int *)a - * (int *)b);
}


void Print_Far_Neighbors2( reax_system *system, control_params *control,
        static_storage *workspace, reax_list **lists )
{
    int i, j, id_i, id_j;
    char fname[MAX_STR];
    FILE *fout;
    reax_list *far_nbrs = lists[FAR_NBRS];

    snprintf( fname, MAX_STR, "%.*s.far_nbrs_lgj", MAX_STR - 14, control->sim_name );
    fout = sfopen( fname, "w", __FILE__, __LINE__ );
    int num = 0;
    int temp[500];

    for ( i = 0; i < system->N; ++i )
    {
        id_i = workspace->orig_id[i];
        num = 0;
        fprintf( fout, "%6d:", id_i);

        for ( j = Start_Index(i, far_nbrs); j < End_Index(i, far_nbrs); ++j )
        {
            id_j = workspace->orig_id[far_nbrs->far_nbr_list[j].nbr];
            temp[num++] = id_j;
        }
        qsort(&temp, num, sizeof(int), fn_qsort_intcmp);
        for (j = 0; j < num; j++)
            fprintf(fout, "%6d", temp[j]);
        fprintf( fout, "\n");
    }

    sfclose( fout, __FILE__, __LINE__ );
}


void Print_Total_Force( reax_system *system, control_params *control,
        simulation_data *data, static_storage *workspace,
        reax_list **lists, output_controls *out_control )
{
    int i;
    char fname[MAX_STR];
    FILE *fout;

    snprintf( fname, MAX_STR, "%.*s.%d.forces", MAX_STR - 10, control->sim_name, data->step );
    fout = sfopen( fname, "w", __FILE__, __LINE__ );

    for ( i = 0; i < system->N; ++i )
    {
        fprintf( fout, "%6d %23.15e %23.15e %23.15e\n",
                workspace->orig_id[i],
                system->atoms[i].f[0], system->atoms[i].f[1], system->atoms[i].f[2] );

//        fprintf( out_control->ftot, "%6d %23.15e %23.15e %23.15e\n",
//                //fprintf(out_control->ftot, "%6d %19.9e %19.9e %19.9e\n",
//                //fprintf(out_control->ftot, "%3d %12.6f %12.6f %12.6f\n",
//                workspace->orig_id[i],
//                system->atoms[i].f[0], system->atoms[i].f[1], system->atoms[i].f[2] );
    }

//    fflush( out_control->ftot );
    fflush( fout );

    sfclose( fout, __FILE__, __LINE__ );
}


/* Write simulation results to various log files */
void Output_Results( reax_system *system, control_params *control,
    simulation_data *data, static_storage *workspace,
    reax_list **lists, output_controls *out_control )
{
    real f_update, t_elapsed;

    t_elapsed = 0.0;

    /* write log files if it is the correct simulation step */
    if ( out_control->log_update_freq > 0 &&
            data->step % out_control->log_update_freq == 0 )
    {
#if defined(TEST_ENERGY) || defined(DEBUG) || defined(DEBUG_FOCUS)
        fprintf( out_control->out,
                 "%-6d%24.15f%24.15f%24.15f%13.5f%13.5f%16.5f%13.5f%13.5f\n",
                 data->step, data->E_Tot, data->E_Pot, data->E_Kin,
                 data->therm.T, control->T, system->box.volume, data->iso_bar.P,
                 (control->P[0] + control->P[1] + control->P[2]) / 3.0 );

        fprintf( out_control->pot,
                 "%-6d%24.15f%24.15f%24.15f%24.15f%24.15f%24.15f%24.15f%24.15f%24.15f%24.15f%24.15f\n",
                 data->step,
                 data->E_BE,
                 data->E_Ov + data->E_Un,  data->E_Lp,
                 data->E_Ang + data->E_Pen, data->E_Coa, data->E_HB,
                 data->E_Tor, data->E_Con,
                 data->E_vdW, data->E_Ele, data->E_Pol );
#else
        fprintf( out_control->out,
                 "%-6d%16.2f%16.2f%16.2f%11.2f%11.2f%13.2f%13.5f%13.5f\n",
                 data->step, data->E_Tot, data->E_Pot, data->E_Kin,
                 data->therm.T, control->T, system->box.volume, data->iso_bar.P,
                 (control->P[0] + control->P[1] + control->P[2]) / 3.0 );

        fprintf( out_control->pot,
                 "%-6d%13.2f%13.2f%13.2f%13.2f%13.2f%13.2f%13.2f%13.2f%13.2f%13.2f%13.2f\n",
                 data->step,
                 data->E_BE,
                 data->E_Ov + data->E_Un,  data->E_Lp,
                 data->E_Ang + data->E_Pen, data->E_Coa, data->E_HB,
                 data->E_Tor, data->E_Con,
                 data->E_vdW, data->E_Ele, data->E_Pol );
#endif

        t_elapsed = Get_Timing_Info( data->timing.total );
        if ( data->step == data->prev_steps )
        {
            f_update = 1.0;
        }
        else
        {
            f_update = 1.0 / out_control->log_update_freq;
        }

        fprintf( out_control->log, "%6d %10.2f %10.2f %10.2f %10.2f %10.2f %10.4f %10.4f %10.2f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",
                 data->step, t_elapsed * f_update,
                 data->timing.nbrs * f_update,
                 data->timing.init_forces * f_update,
                 data->timing.bonded * f_update,
                 data->timing.nonb * f_update,
                 data->timing.cm * f_update,
                 data->timing.cm_sort_mat_rows * f_update,
                 (double) data->timing.cm_solver_iters * f_update,
                 data->timing.cm_solver_pre_comp * f_update,
                 data->timing.cm_solver_pre_app * f_update,
                 data->timing.cm_solver_spmv * f_update,
                 data->timing.cm_solver_vector_ops * f_update,
                 data->timing.cm_solver_orthog * f_update,
                 data->timing.cm_solver_tri_solve * f_update );

        data->timing.total = Get_Time( );
        data->timing.nbrs = 0.0;
        data->timing.init_forces = 0.0;
        data->timing.bonded = 0.0;
        data->timing.nonb = 0.0;
        data->timing.cm = 0.0;
        data->timing.cm_sort_mat_rows = 0.0;
        data->timing.cm_solver_pre_comp = 0.0;
        data->timing.cm_solver_pre_app = 0.0;
        data->timing.cm_solver_iters = 0;
        data->timing.cm_solver_spmv = 0.0;
        data->timing.cm_solver_vector_ops = 0.0;
        data->timing.cm_solver_orthog = 0.0;
        data->timing.cm_solver_tri_solve = 0.0;

        fflush( out_control->out );
        fflush( out_control->pot );
        fflush( out_control->log );

        /* output pressure */
        if ( control->ensemble == sNPT || control->ensemble == iNPT ||
                control->ensemble == aNPT || control->compute_pressure == TRUE )
        {
#if defined(DEBUG) || defined(DEBUG_FOCUS)
            fprintf( out_control->prs, "%-8d %13.6f %13.6f %13.6f %13.6f %13.6f %13.6f\n",
                     data->step,
                     data->kin_press[0][0], data->kin_press[1][1],+ data->kin_press[2][2],
                     data->press[0][0], data->press[1][1], data->press[2][2] );
#endif

            fprintf( out_control->prs,
                     "%-8d %13.6f %13.6f %13.6f %13.6f %13.6f %13.6f %13.6f %13.6f\n",
                     data->step,
                     system->box.box_norms[0], system->box.box_norms[1],
                     system->box.box_norms[2],
                     data->tot_press[0], data->tot_press[1], data->tot_press[2],
                     (data->tot_press[0] + data->tot_press[1] + data->tot_press[2]) / 3.0, system->box.volume );
            fflush( out_control->prs );
        }
    }

    if ( out_control->write_steps > 0 &&
            data->step % out_control->write_steps == 0 )
    {
        out_control->append_traj_frame( system, control, data,
                workspace, lists, out_control );
    }
}


void Print_Linear_System( reax_system *system, control_params *control,
        static_storage *workspace, int step )
{
    int i, j;
    char fname[100];
    sparse_matrix *H;
    FILE *out;

    snprintf( fname, 100, "%.*s.state%10d.out", 79, control->sim_name, step );
    out = sfopen( fname, "w", __FILE__, __LINE__ );

    for ( i = 0; i < system->N_cm; i++ )
        fprintf( out, "%6d%2d%24.15e%24.15e%24.15e%24.15e%24.15e%24.15e%24.15e\n",
                 workspace->orig_id[i], system->atoms[i].type,
                 system->atoms[i].x[0], system->atoms[i].x[1],
                 system->atoms[i].x[2],
                 workspace->s[0][i], workspace->b_s[i],
                 workspace->t[0][i], workspace->b_t[i]  );
    sfclose( out, __FILE__, __LINE__ );

    // snprintf( fname, 100, "x2_%d", step );
    // out = sfopen( fname, "w", __FILE__, __LINE__ );
    // for( i = 0; i < system->N; i++ )
    // fprintf( out, "%g\n", workspace->s_t[i+system->N] );
    // sfclose( out, __FILE__, __LINE__ );

    snprintf( fname, 100, "%.*s.H%10d.out", 83, control->sim_name, step );
    out = sfopen( fname, "w", __FILE__, __LINE__ );
    H = &workspace->H;

    for ( i = 0; i < system->N_cm; ++i )
    {
        for ( j = H->start[i]; j < H->start[i + 1] - 1; ++j )
        {
            fprintf( out, "%6d%6d %24.15e\n",
                     workspace->orig_id[i], workspace->orig_id[H->j[j]],
                     H->val[j] );

            fprintf( out, "%6d%6d %24.15e\n",
                     workspace->orig_id[H->j[j]], workspace->orig_id[i],
                     H->val[j] );
        }
        // the diagonal entry
        fprintf( out, "%6d%6d %24.15e\n",
                 workspace->orig_id[i], workspace->orig_id[i], H->val[j] );
    }

    sfclose( out, __FILE__, __LINE__ );

    snprintf( fname, 100, "%.*s.H_sp%10d.out", 80, control->sim_name, step );
    out = sfopen( fname, "w", __FILE__, __LINE__ );
    H = &workspace->H_sp;

    for ( i = 0; i < system->N_cm; ++i )
    {
        for ( j = H->start[i]; j < H->start[i + 1] - 1; ++j )
        {
            fprintf( out, "%6d%6d %24.15e\n",
                     workspace->orig_id[i], workspace->orig_id[H->j[j]],
                     H->val[j] );

            fprintf( out, "%6d%6d %24.15e\n",
                     workspace->orig_id[H->j[j]], workspace->orig_id[i],
                     H->val[j] );
        }
        // the diagonal entry
        fprintf( out, "%6d%6d %24.15e\n",
                 workspace->orig_id[i], workspace->orig_id[i], H->val[j] );
    }

    sfclose( out, __FILE__, __LINE__ );

    /*snprintf( fname, 100, "%.*s.b_s%10d", 84, control->sim_name, step );
      out = sfopen( fname, "w", __FILE__, __LINE__ );
      for( i = 0; i < system->N; i++ )
      fprintf( out, "%12.7f\n", workspace->b_s[i] );
      sfclose( out, __FILE__, __LINE__ );

      snprintf( fname, 100, "%.*s.b_t%10d", 84, control->sim_name, step );
      out = sfopen( fname, "w", __FILE__, __LINE__ );
      for( i = 0; i < system->N; i++ )
      fprintf( out, "%12.7f\n", workspace->b_t[i] );
      sfclose( out, __FILE__, __LINE__ );*/
}


void Print_Charges( reax_system *system, control_params *control,
        static_storage *workspace, int step )
{
    int i;
    char fname[100];
    FILE *fout;

    snprintf( fname, 100, "%.*s.q%010d", 87, control->sim_name, step );
    fout = sfopen( fname, "w", __FILE__, __LINE__ );

    for ( i = 0; i < system->N; ++i )
    {
        fprintf( fout, "%6d%12.7f%12.7f%12.7f\n",
                 workspace->orig_id[i],
                 workspace->s[0][i], workspace->t[0][i], system->atoms[i].q );
    }

    sfclose( fout, __FILE__, __LINE__ );
}


void Print_Soln( static_storage *workspace,
        real *x, real *b_prm, real *b, int N )
{
    int i;

    fprintf( stdout, "%6s%10s%10s%10s\n", "id", "x", "b_prm", "b" );

    for ( i = 0; i < N; ++i )
    {
        fprintf( stdout, "%6d%10.4f%10.4f%10.4f\n",
                 workspace->orig_id[i], x[i], b_prm[i], b[i] );
    }

    fflush( stdout );
}


void Print_Sparse_Matrix( sparse_matrix *A )
{
    int i, j;

    for ( i = 0; i < A->n; ++i )
    {
        fprintf( stderr, "i:%d  j(val):", i );
        for ( j = A->start[i]; j < A->start[i + 1]; ++j )
        {
            fprintf( stderr, "%d(%.4f) ", A->j[j], A->val[j] );
        }
        fprintf( stderr, "\n" );
    }
}


void Print_Sparse_Matrix2( sparse_matrix *A, char *fname, char *mode )
{
    int i, j;
    FILE *f;
   
    if ( mode == NULL )
    {
        f = sfopen( fname, "w", __FILE__, __LINE__ );
    }
    else
    {
        f = sfopen( fname, mode, __FILE__, __LINE__ );
    }

    for ( i = 0; i < A->n; ++i )
    {
        /* off-diagonals */
        for ( j = A->start[i]; j < A->start[i + 1] - 1; ++j )
        {
            //Convert 0-based to 1-based (for Matlab)
            fprintf( f, "%6d %6d %24.15f\n", i + 1, A->j[j] + 1, A->val[j] );
            /* print symmetric entry */
//            fprintf( f, "%6d %6d %24.15f\n", A->j[j] + 1, i + 1, A->val[j] );
        }

        /* diagonal */
        fprintf( f, "%6d %6d %24.15f\n", i + 1, A->j[A->start[i + 1] - 1] + 1, A->val[A->start[i + 1] - 1] );
    }

    sfclose( f, __FILE__, __LINE__ );
}


/* Read sparse matrix in COO format (1-based indexing)
 * and store in a CSR sparse matrix (preallocated)
 *
 * Note: the file must be sorted in increasing order of
 * columns than rows */
void Read_Sparse_Matrix2( sparse_matrix *A, char *fname )
{
    int top, cur_row, row, col;
    real val;
    FILE *f;
   
    f = sfopen( fname, "r", __FILE__, __LINE__ );
    top = 0;
    cur_row = 0;

    A->start[cur_row] = top;

    while ( fscanf( f, "%6d %6d %24lf", &row, &col, &val ) == 3 )
    {
        /* assumption: every row has at least one nonzero (diagonal) */
        if ( cur_row != row - 1 )
        {
            cur_row++;
            A->start[cur_row] = top;
        }

        A->j[top] = col - 1;
        A->val[top] = val;
        top++;
    }

    A->start[A->n] = top;

    sfclose( f, __FILE__, __LINE__ );
}


/* Read permuation matrix in COO format (1-based indexing)
 * and store in a dense vector (preallocated)
 *
 * Note: the file must be sorted in increasing order of
 * columns than rows */
void Read_Permutation_Matrix( unsigned int *v, char *fname )
{
    unsigned int row, col;
    double val;
    FILE *f;
   
    f = sfopen( fname, "r", __FILE__, __LINE__ );

    while ( fscanf( f, "%6u %6u %24lf", &row, &col, &val ) == 3 )
    {
        v[row - 1] = col - 1;
    }

    sfclose( f, __FILE__, __LINE__ );
}


/* Note: watch out for portability issues with endianness
 * due to serialization of numeric types (integer, IEEE 754) */
void Print_Sparse_Matrix_Binary( sparse_matrix *A, char *fname )
{
    int i, j, temp;
    FILE *f;
   
    f = sfopen( fname, "wb", __FILE__, __LINE__ );

    /* header: # rows, # nonzeros */
    fwrite( &A->n, sizeof(unsigned int), 1, f );
    fwrite( &A->start[A->n], sizeof(unsigned int), 1, f );

    /* row pointers */
    for ( i = 0; i <= A->n; ++i )
    {
        //Convert 0-based to 1-based (for Matlab)
        temp = A->start[i] + 1;
        fwrite( &temp, sizeof(unsigned int), 1, f );
    }

    /* column indices and non-zeros */
    for ( i = 0; i <= A->n; ++i )
    {
        for ( j = A->start[i]; j < A->start[i + 1]; ++j )
        {
            //Convert 0-based to 1-based (for Matlab)
            temp = A->j[j] + 1;
            fwrite( &temp, sizeof(unsigned int), 1, f );
            fwrite( &(A->val[j]), sizeof(real), 1, f );
        }
    }

    sfclose( f, __FILE__, __LINE__ );
}


void Print_Bonds( reax_system *system, reax_list *bonds, char *fname )
{
    int i, pj;
    bond_data *pbond;
    bond_order_data *bo_ij;
    FILE *f = sfopen( fname, "w", __FILE__, __LINE__ );

    for ( i = 0; i < system->N; ++i )
    {
        for ( pj = Start_Index(i, bonds); pj < End_Index(i, bonds); ++pj )
        {
            pbond = &bonds->bond_list[pj];
            bo_ij = &pbond->bo_data;
            //fprintf( f, "%6d%6d%23.15e%23.15e%23.15e%23.15e%23.15e\n",
            //       i+1, pbond->nbr+1, pbond->d,
            //       bo_ij->BO, bo_ij->BO_s, bo_ij->BO_pi, bo_ij->BO_pi2 );
            fprintf( f, "%6d%6d %9.5f %9.5f\n",
                     i + 1, pbond->nbr + 1, pbond->d, bo_ij->BO );
        }
    }

    sfclose( f, __FILE__, __LINE__ );
}


void Print_Bond_List2( reax_system *system, reax_list *bonds, char *fname )
{
    int i, j, id_i, id_j, nbr, pj;
    FILE *f = sfopen( fname, "w", __FILE__, __LINE__ );
    int temp[500];
    int num = 0;

    for ( i = 0; i < system->N; ++i )
    {
        num = 0;
        id_i = i + 1; //system->atoms[i].orig_id;
        fprintf( f, "%6d:", id_i);
        for ( pj = Start_Index(i, bonds); pj < End_Index(i, bonds); ++pj )
        {
            nbr = bonds->bond_list[pj].nbr;
            id_j = nbr + 1; //system->my_atoms[nbr].orig_id;
            if ( id_i < id_j )
            {
                temp[num++] = id_j;
            }
        }

        qsort( &temp, num, sizeof(int), fn_qsort_intcmp );
        for ( j = 0; j < num; j++ )
        {
            fprintf( f, "%6d", temp[j] );
        }
        fprintf( f, "\n" );
    }
}
