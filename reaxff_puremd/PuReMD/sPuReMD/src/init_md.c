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

#include "init_md.h"

#include "allocate.h"
#include "box.h"
#include "forces.h"
#include "grid.h"
#include "integrate.h"
#if defined(TEST_FORCES)
  #include "io_tools.h"
#endif
#include "neighbors.h"
#include "list.h"
#include "lookup.h"
#include "reset_tools.h"
#include "system_props.h"
#include "tool_box.h"
#include "vector.h"


static void Generate_Initial_Velocities( reax_system * const system,
        control_params const * const control, real T )
{
    int i;
    real scale, norm;

    if ( T <= 0.1 || control->random_vel == FALSE )
    {
        /* warnings if conflicts between initial temperature and control file parameter */
        if ( control->random_vel == TRUE )
        {
            fprintf( stderr, "[ERROR] conflicting control file parameters\n" );
            fprintf( stderr, "[INFO] random_vel = 1 and small initial temperature (t_init = %f)\n", T );
            fprintf( stderr, "[INFO] set random_vel = 0 to resolve this (atom initial velocites set to zero)\n" );
            exit( INVALID_INPUT );
        }
        else if ( T > 0.1 )
        {
            fprintf( stderr, "[ERROR] conflicting control file paramters\n" );
            fprintf( stderr, "[INFO] random_vel = 0 and large initial temperature (t_init = %f)\n", T );
            fprintf( stderr, "[INFO] set random_vel = 1 to resolve this (random atom initial velocites according to t_init)\n" );
            exit( INVALID_INPUT );
        }

        for ( i = 0; i < system->N; i++ )
        {
            rvec_MakeZero( system->atoms[i].v );
        }
    }
    else
    {
        if ( T <= 0.0 )
        {
            fprintf( stderr, "[ERROR] random atom initial velocities specified with invalid temperature (%f). Terminating...\n",
                  T );
            exit( INVALID_INPUT );
        }

        for ( i = 0; i < system->N; i++ )
        {
            rvec_Random( system->atoms[i].v );

            norm = rvec_Norm_Sqr( system->atoms[i].v );
            scale = SQRT( system->reax_param.sbp[ system->atoms[i].type ].mass
                    * norm / (3.0 * K_B * T) );

            rvec_Scale( system->atoms[i].v, 1.0 / scale, system->atoms[i].v );

#if defined(DEBUG_FOCUS)
            fprintf( stderr, "[INFO] atom %d, scale = %f, v = (%f, %f, %f)\n",
                    i, scale,
                    system->atoms[i].v[0],
                    system->atoms[i].v[1],
                    system->atoms[i].v[2] );
#endif
        }
    }
}


static void Init_System( reax_system * const system, control_params * const control,
        simulation_data * const data, static_storage * const workspace, int realloc )
{
    int i;
    rvec dx;

    system->allocated = TRUE;

    if ( control->restart == FALSE )
    {
        Reset_Atomic_Forces( system );
    }

    Compute_Total_Mass( system, data );
    Compute_Center_of_Mass( system, data );

    /* just fit the atoms to the periodic box */
    if ( control->reposition_atoms == 0 )
    {
        rvec_MakeZero( dx );
    }
    /* put the center of mass to the center of the box */
    else if ( control->reposition_atoms == 1 )
    {
        rvec_Scale( dx, 0.5, system->box.box_norms );
        rvec_ScaledAdd( dx, -1.0, data->xcm );
    }
    /* put the center of mass to the origin */
    else if ( control->reposition_atoms == 2 )
    {
        rvec_Scale( dx, -1.0, data->xcm );
    }
    else
    {
        fprintf( stderr, "[ERROR] Unknown option for reposition_atoms (%d). Terminating...\n",
              control->reposition_atoms );
        exit( UNKNOWN_OPTION );
    }

    for ( i = 0; i < system->N; ++i )
    {
        /* re-map the atom positions to fall within the simulation box,
         * where the corners of the box are (0,0,0) and (d_x, d_y, d_z)
         * with d_i being the box length along dimension i */
        Update_Atom_Position_Periodic( system->atoms[i].x, dx,
                system->atoms[i].rel_map, &system->box );

        /* zero out rel_map (which was set by the above function) */
        ivec_MakeZero( system->atoms[i].rel_map );
    }

    /* Initialize velocities so that desired init T can be attained */
    if ( control->restart == FALSE
            || (control->restart == TRUE && control->random_vel == TRUE) )
    {
        Generate_Initial_Velocities( system, control, control->T_init );
    }

    Setup_Grid( system );

    Bin_Atoms( system, workspace );

#if defined(REORDER_ATOMS)
    Reorder_Atoms( system, workspace, control );
#endif

    if ( realloc == TRUE )
    {
        /* list management */
        system->bonds = smalloc( sizeof(int) * system->N_max, __FILE__, __LINE__ );

        system->hbonds = smalloc( sizeof(int) * system->N_max, __FILE__, __LINE__ );
    }
}


static void Init_Simulation_Data( reax_system * const system,
        control_params * const control, simulation_data * const data,
        evolve_function * const Evolve, int realloc )
{
#if defined(_OPENMP)
    if ( realloc == TRUE && (control->ensemble == sNPT || control->ensemble == iNPT
            || control->ensemble == aNPT || control->compute_pressure == TRUE) )
    {
        data->press_local = smalloc( sizeof( rtensor ) * control->num_threads,
               __FILE__, __LINE__ );
    }
#endif

    Reset_Pressures( control, data );
    Reset_Energies( data );

    data->therm.T = 0.0;
    data->therm.xi = 0.0;
    data->therm.v_xi = 0.0;
    data->therm.v_xi_old = 0.0;
    data->therm.G_xi = 0.0;

    /* initialize for non-restarted run,
     * code in restart.c (restart file parser) initializes otherwise */
    if ( control->restart == FALSE )
    {
        data->step = 0;
        data->prev_steps = 0;
    }

    data->time = 0.0;

    switch ( control->ensemble )
    {
    case NVE:
        data->N_f = 3.0 * system->N;
        *Evolve = &Velocity_Verlet_NVE;
        break;

    case bNVT:
        data->N_f = 3.0 * system->N;
        *Evolve = &Velocity_Verlet_Berendsen_NVT;
        break;

    case nhNVT:
        data->N_f = 3.0 * system->N + 1.0;
        *Evolve = &Velocity_Verlet_Nose_Hoover_NVT_Klein;

        if ( control->restart == FALSE
                || (control->restart == TRUE && control->random_vel == TRUE) )
        {
            data->therm.G_xi = control->Tau_T * (2.0 * data->E_Kin
                    - data->N_f * K_B / F_CONV * control->T);
            data->therm.v_xi = data->therm.G_xi * control->dt;
            data->therm.v_xi_old = 0.0;
            data->therm.xi = 0.0;
        }
        break;

    /* anisotropic NPT */
    case aNPT:
        fprintf( stderr, "[ERROR] THIS OPTION IS NOT YET IMPLEMENTED! TERMINATING...\n" );
        exit( UNKNOWN_OPTION );

        data->N_f = 3.0 * system->N + 9.0;
        *Evolve = &Velocity_Verlet_Berendsen_Isotropic_NPT;

        if ( control->restart == FALSE )
        {
            data->therm.G_xi = control->Tau_T * (2.0 * data->E_Kin -
                    data->N_f * K_B / F_CONV * control->T);
            data->therm.v_xi = data->therm.G_xi * control->dt;
            data->iso_bar.eps = 1.0 / 3.0 * LOG( system->box.volume );
//            data->inv_W = 1.0 / (data->N_f * K_B * control->T * SQR(control->Tau_P));
//            Compute_Pressure( system, data, workspace );
        }
        break;

    /* semi-isotropic NPT */
    case sNPT:
        data->N_f = 3.0 * system->N + 4.0;
        *Evolve = &Velocity_Verlet_Berendsen_Semi_Isotropic_NPT;
        break;

    /* isotropic NPT */
    case iNPT:
        data->N_f = 3.0 * system->N + 2.0;
        *Evolve = &Velocity_Verlet_Berendsen_Isotropic_NPT;
        break;

    default:
        fprintf( stderr, "[ERROR] Unknown ensemble type (%d). Terminating...\n", control->ensemble );
        exit( UNKNOWN_OPTION );
        break;
    }

    Compute_Kinetic_Energy( system, data );

    /* init timing info */
    data->timing.start = Get_Time( );
    data->timing.total = data->timing.start;
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
    data->timing.cm_last_pre_comp = 0.0;
    data->timing.cm_total_loss = 0.0;
    data->timing.cm_optimum = 0.0;
}


/* Initialize taper function applied to van der Waals and Coulomb interactions */
static void Init_Taper( control_params const * const control,
        static_storage * const workspace )
{
    real d1, d7;
    real swa, swa2, swa3;
    real swb, swb2, swb3;

    swa = control->nonb_low;
    swb = control->nonb_cut;

    if ( FABS( swa ) > 0.01 )
    {
        fprintf( stderr, "[WARNING] non-zero value for lower Taper-radius cutoff (%f)\n", swa );
    }

    if ( swb < 0.0 )
    {
        fprintf( stderr, "[ERROR] Negative value for upper Taper-radius cutoff\n" );
        exit( INVALID_INPUT );
    }
    else if ( swb < 5.0 )
    {
        fprintf( stderr, "[WARNING] Low value for upper Taper-radius cutoff (%f)\n", swb );
    }

    d1 = swb - swa;
    d7 = POW( d1, 7.0 );
    swa2 = SQR( swa );
    swa3 = swa2 * swa;
    swb2 = SQR( swb );
    swb3 = swb2 * swb;

    workspace->Tap[7] =  20.0 / d7;
    workspace->Tap[6] = -70.0 * (swa + swb) / d7;
    workspace->Tap[5] =  84.0 * (swa2 + 3.0 * swa * swb + swb2) / d7;
    workspace->Tap[4] = -35.0 * (swa3 + 9.0 * swa2 * swb + 9.0 * swa * swb2 + swb3 ) / d7;
    workspace->Tap[3] = 140.0 * (swa3 * swb + 3.0 * swa2 * swb2 + swa * swb3 ) / d7;
    workspace->Tap[2] = -210.0 * (swa3 * swb2 + swa2 * swb3) / d7;
    workspace->Tap[1] = 140.0 * swa3 * swb3 / d7;
    workspace->Tap[0] = (-35.0 * swa3 * swb2 * swb2 + 21.0 * swa2 * swb3 * swb2 +
            7.0 * swa * swb3 * swb3 + swb3 * swb3 * swb ) / d7;
}


static void Init_Workspace( reax_system * const system,
        control_params const * const control, static_storage * const workspace,
        int realloc )
{
    int i;

    workspace->allocated = TRUE;

    if ( realloc == TRUE )
    {
        /* bond order related storage  */
        workspace->total_bond_order = smalloc( system->N_max * sizeof( real ),
               __FILE__, __LINE__ );
        workspace->Deltap = smalloc( system->N_max * sizeof( real ),
               __FILE__, __LINE__ );
        workspace->Deltap_boc = smalloc( system->N_max * sizeof( real ),
               __FILE__, __LINE__ );
        workspace->dDeltap_self = smalloc( system->N_max * sizeof( rvec ),
               __FILE__, __LINE__ );

        workspace->Delta = smalloc( system->N_max * sizeof( real ),
               __FILE__, __LINE__ );
        workspace->Delta_lp = smalloc( system->N_max * sizeof( real ),
               __FILE__, __LINE__ );
        workspace->Delta_lp_temp = smalloc( system->N_max * sizeof( real ),
               __FILE__, __LINE__ );
        workspace->dDelta_lp = smalloc( system->N_max * sizeof( real ),
               __FILE__, __LINE__ );
        workspace->dDelta_lp_temp = smalloc( system->N_max * sizeof( real ),
               __FILE__, __LINE__ );
        workspace->Delta_e = smalloc( system->N_max * sizeof( real ),
               __FILE__, __LINE__ );
        workspace->Delta_boc = smalloc( system->N_max * sizeof( real ),
               __FILE__, __LINE__ );
        workspace->nlp = smalloc( system->N_max * sizeof( real ),
               __FILE__, __LINE__ );
        workspace->nlp_temp = smalloc( system->N_max * sizeof( real ),
               __FILE__, __LINE__ );
        workspace->Clp = smalloc( system->N_max * sizeof( real ),
               __FILE__, __LINE__ );
        workspace->CdDelta = smalloc( system->N_max * sizeof( real ),
               __FILE__, __LINE__ );
        workspace->vlpex = smalloc( system->N_max * sizeof( real ),
               __FILE__, __LINE__ );
    }

    /* charge method storage */
    switch ( control->charge_method )
    {
        case QEQ_CM:
            system->N_cm = system->N;
            if ( realloc == TRUE || system->N_cm > system->N_cm_max )
            {
                system->N_cm_max = system->N_max;
            }
            break;
        case EE_CM:
            system->N_cm = system->N
                + system->num_molec_charge_constraints
                + system->num_custom_charge_constraints
                + (system->num_molec_charge_constraints == 0 && system->num_custom_charge_constraints == 0 ? 1 : 0);
            if ( realloc == TRUE || system->N_cm > system->N_cm_max )
            {
                system->N_cm_max = system->N_max
                    + system->num_molec_charge_constraints
                    + system->num_custom_charge_constraints
                    + (system->num_molec_charge_constraints == 0 && system->num_custom_charge_constraints == 0 ? 1 : 0);
            }
            break;
        case ACKS2_CM:
            system->N_cm = 2 * system->N + 2;
            if ( realloc == TRUE || system->N_cm > system->N_cm_max )
            {
                system->N_cm_max = 2 * system->N_max + 2;
            }
            break;
        default:
            fprintf( stderr, "[ERROR] Unknown charge method type. Terminating...\n" );
            exit( INVALID_INPUT );
            break;
    }

    if ( realloc == TRUE )
    {
        workspace->Hdia_inv = NULL;

        if ( control->cm_solver_pre_comp_type == ICHOLT_PC
                || (control->cm_solver_pre_comp_type == ILUT_PC && control->cm_solver_pre_comp_droptol > 0.0 )
                || control->cm_solver_pre_comp_type == ILUTP_PC
                || control->cm_solver_pre_comp_type == FG_ILUT_PC )
        {
            workspace->droptol = scalloc( system->N_cm_max, sizeof( real ),
                    __FILE__, __LINE__ );
        }

        workspace->b_s = scalloc( system->N_cm_max, sizeof( real ),
                __FILE__, __LINE__ );
        workspace->b_t = scalloc( system->N_cm_max, sizeof( real ),
                __FILE__, __LINE__ );
        workspace->b_prc = scalloc( system->N_cm_max * 2, sizeof( real ),
                __FILE__, __LINE__ );
        workspace->b_prm = scalloc( system->N_cm_max * 2, sizeof( real ),
                __FILE__, __LINE__ );
        workspace->s = scalloc( 5, sizeof( real* ),
                __FILE__, __LINE__ );
        workspace->t = scalloc( 5, sizeof( real* ),
                __FILE__, __LINE__ );
        for ( i = 0; i < 5; ++i )
        {
            workspace->s[i] = scalloc( system->N_cm_max, sizeof( real ),
                    __FILE__, __LINE__ );
            workspace->t[i] = scalloc( system->N_cm_max, sizeof( real ),
                    __FILE__, __LINE__ );
        }
    }

    switch ( control->charge_method )
    {
        case QEQ_CM:
            for ( i = 0; i < system->N; ++i )
            {
                workspace->b_s[i] = -1.0 * system->reax_param.sbp[ system->atoms[i].type ].chi;
                workspace->b_t[i] = -1.0;
            }
            break;

        case EE_CM:
            for ( i = 0; i < system->N; ++i )
            {
                workspace->b_s[i] = -system->reax_param.sbp[ system->atoms[i].type ].chi;
            }

            if ( system->num_molec_charge_constraints == 0
                    && system->num_custom_charge_constraints == 0 )
            {
                workspace->b_s[system->N] = control->cm_q_net;
            }
            else
            {
                if ( system->num_molec_charge_constraints > 0 )
                {
                    for ( i = 0; i < system->num_molec_charge_constraints; ++i )
                    {
                        workspace->b_s[system->N + i] = system->molec_charge_constraints[i];
                    }
                }
                if ( system->num_custom_charge_constraints > 0 )
                {
                    for ( i = 0; i < system->num_custom_charge_constraints; ++i )
                    {
                        workspace->b_s[system->N + system->num_molec_charge_constraints + i]
                            = system->custom_charge_constraint_rhs[i];
                    }
                }
            }
            break;

        case ACKS2_CM:
            for ( i = 0; i < system->N; ++i )
            {
                workspace->b_s[i] = -system->reax_param.sbp[ system->atoms[i].type ].chi;
            }

            /* Non-zero total charge can lead to unphysical results.
             * As such, set the ACKS2 reference charge of every atom
             * to the total charge divided by the number of atoms.
             * Except for trivial cases, this leads to fractional
             * reference charges, which is usually not desirable. */
            for ( i = 0; i < system->N; ++i )
            {
                workspace->b_s[system->N + i] = control->cm_q_net / system->N;
            }

            /* system charge defines the total charge constraint */
            workspace->b_s[system->N_cm - 1] = control->cm_q_net;
            break;

        default:
            fprintf( stderr, "[ERROR] Unknown charge method type. Terminating...\n" );
            exit( INVALID_INPUT );
            break;
    }

    if ( realloc == TRUE )
    {
        switch ( control->cm_solver_type )
        {
            case GMRES_S:
            case GMRES_H_S:
                workspace->y = scalloc( control->cm_solver_restart + 1, sizeof( real ),
                        __FILE__, __LINE__ );
                workspace->z = scalloc( control->cm_solver_restart + 1, sizeof( real ),
                        __FILE__, __LINE__ );
                workspace->g = scalloc( control->cm_solver_restart + 1, sizeof( real ),
                        __FILE__, __LINE__ );
                workspace->h = scalloc( control->cm_solver_restart + 1, sizeof( real*),
                        __FILE__, __LINE__ );
                workspace->hs = scalloc( control->cm_solver_restart + 1, sizeof( real ),
                        __FILE__, __LINE__ );
                workspace->hc = scalloc( control->cm_solver_restart + 1, sizeof( real ),
                        __FILE__, __LINE__ );
                workspace->rn = scalloc( control->cm_solver_restart + 1, sizeof( real*),
                        __FILE__, __LINE__ );
                workspace->v = scalloc( control->cm_solver_restart + 1, sizeof( real*),
                        __FILE__, __LINE__ );

                for ( i = 0; i < control->cm_solver_restart + 1; ++i )
                {
                    workspace->h[i] = scalloc( control->cm_solver_restart + 1, sizeof( real ),
                            __FILE__, __LINE__ );
                    workspace->rn[i] = scalloc( system->N_cm_max * 2, sizeof( real ),
                            __FILE__, __LINE__ );
                    workspace->v[i] = scalloc( system->N_cm_max, sizeof( real ),
                            __FILE__, __LINE__ );
                }

                workspace->r = scalloc( system->N_cm_max, sizeof( real ),
                        __FILE__, __LINE__ );
                workspace->d = scalloc( system->N_cm_max, sizeof( real ),
                        __FILE__, __LINE__ );
                workspace->q = scalloc( system->N_cm_max, sizeof( real ),
                        __FILE__, __LINE__ );
                workspace->p = scalloc( system->N_cm_max, sizeof( real ),
                        __FILE__, __LINE__ );
                break;

            case CG_S:
                workspace->r = scalloc( system->N_cm_max, sizeof( real ),
                        __FILE__, __LINE__ );
                workspace->d = scalloc( system->N_cm_max, sizeof( real ),
                        __FILE__, __LINE__ );
                workspace->q = scalloc( system->N_cm_max, sizeof( real ),
                        __FILE__, __LINE__ );
                workspace->p = scalloc( system->N_cm_max, sizeof( real ),
                        __FILE__, __LINE__ );
                break;

            case SDM_S:
                workspace->r = scalloc( system->N_cm_max, sizeof( real ),
                        __FILE__, __LINE__ );
                workspace->d = scalloc( system->N_cm_max, sizeof( real ),
                        __FILE__, __LINE__ );
                workspace->q = scalloc( system->N_cm_max, sizeof( real ),
                        __FILE__, __LINE__ );
                break;

            case BiCGStab_S:
                workspace->r = scalloc( system->N_cm_max, sizeof( real ),
                        __FILE__, __LINE__ );
                workspace->r_hat = scalloc( system->N_cm_max, sizeof( real ),
                        __FILE__, __LINE__ );
                workspace->d = scalloc( system->N_cm_max, sizeof( real ),
                        __FILE__, __LINE__ );
                workspace->q = scalloc( system->N_cm_max, sizeof( real ),
                        __FILE__, __LINE__ );
                workspace->q_hat = scalloc( system->N_cm_max, sizeof( real ),
                        __FILE__, __LINE__ );
                workspace->p = scalloc( system->N_cm_max, sizeof( real ),
                        __FILE__, __LINE__ );
                workspace->y = scalloc( system->N_cm_max, sizeof( real ),
                        __FILE__, __LINE__ );
                workspace->z = scalloc( system->N_cm_max, sizeof( real ),
                        __FILE__, __LINE__ );
                workspace->g = scalloc( system->N_cm_max, sizeof( real ),
                        __FILE__, __LINE__ );
                break;

            default:
                fprintf( stderr, "[ERROR] Unknown charge method linear solver type. Terminating...\n" );
                exit( INVALID_INPUT );
                break;
        }

#if defined(_OPENMP)
        /* SpMV related */
        workspace->b_local = smalloc( control->num_threads * system->N_cm_max * sizeof(real),
                __FILE__, __LINE__ );
#endif
    }

    /* level scheduling related */
    workspace->levels_L = 1;
    workspace->levels_U = 1;
    if ( realloc == TRUE )
    {
        if ( control->cm_solver_pre_app_type == TRI_SOLVE_LEVEL_SCHED_PA ||
                control->cm_solver_pre_app_type == TRI_SOLVE_GC_PA )
        {
            workspace->row_levels_L = smalloc( system->N_cm_max * sizeof(unsigned int),
                    __FILE__, __LINE__ );
            workspace->level_rows_L = smalloc( system->N_cm_max * sizeof(unsigned int),
                    __FILE__, __LINE__ );
            workspace->level_rows_cnt_L = smalloc( (system->N_cm_max + 1) * sizeof(unsigned int),
                    __FILE__, __LINE__ );
            workspace->row_levels_U = smalloc( system->N_cm_max * sizeof(unsigned int),
                    __FILE__, __LINE__ );
            workspace->level_rows_U = smalloc( system->N_cm_max * sizeof(unsigned int),
                    __FILE__, __LINE__ );
            workspace->level_rows_cnt_U = smalloc( (system->N_cm_max + 1) * sizeof(unsigned int),
                    __FILE__, __LINE__ );
            workspace->top = smalloc( (system->N_cm_max + 1) * sizeof(unsigned int),
                    __FILE__, __LINE__ );
        }
        else
        {
            workspace->row_levels_L = NULL;
            workspace->level_rows_L = NULL;
            workspace->level_rows_cnt_L = NULL;
            workspace->row_levels_U = NULL;
            workspace->level_rows_U = NULL;
            workspace->level_rows_cnt_U = NULL;
            workspace->top = NULL;
        }
    }

    /* graph coloring related */
    workspace->recolor_cnt = 0;
    if ( realloc == TRUE )
    {
        if ( control->cm_solver_pre_app_type == TRI_SOLVE_GC_PA )
        {
            workspace->color = smalloc( sizeof(unsigned int) * system->N_cm_max,
                    __FILE__, __LINE__ );
            workspace->to_color = smalloc( sizeof(unsigned int) * system->N_cm_max,
                    __FILE__, __LINE__ );
            workspace->conflict = smalloc( sizeof(unsigned int) * system->N_cm_max,
                    __FILE__, __LINE__ );
            workspace->conflict_cnt = smalloc( sizeof(unsigned int) * (control->num_threads + 1),
                    __FILE__, __LINE__ );
            workspace->recolor = smalloc( sizeof(unsigned int) * system->N_cm_max,
                    __FILE__, __LINE__ );
            workspace->color_top = smalloc( sizeof(unsigned int) * (system->N_cm_max + 1),
                    __FILE__, __LINE__ );
            workspace->permuted_row_col = smalloc( sizeof(unsigned int) * system->N_cm_max,
                    __FILE__, __LINE__ );
            workspace->permuted_row_col_inv = smalloc( sizeof(unsigned int) * system->N_cm_max,
                    __FILE__, __LINE__ );
        }
        else
        {
            workspace->color = NULL;
            workspace->to_color = NULL;
            workspace->conflict = NULL;
            workspace->conflict_cnt = NULL;
            workspace->recolor = NULL;
            workspace->color_top = NULL;
            workspace->permuted_row_col = NULL;
            workspace->permuted_row_col_inv = NULL;
        }

        /* graph coloring related OR ILUTP preconditioner */
        if ( control->cm_solver_pre_app_type == TRI_SOLVE_GC_PA 
                || control->cm_solver_pre_comp_type == ILUTP_PC )
        {
            workspace->y_p = smalloc( sizeof(real) * system->N_cm_max, __FILE__, __LINE__ );
            workspace->x_p = smalloc( sizeof(real) * system->N_cm_max, __FILE__, __LINE__ );
        }
        else
        {
            workspace->y_p = NULL;
            workspace->x_p = NULL;
        }

        /* Jacobi iteration related */
        if ( control->cm_solver_pre_app_type == JACOBI_ITER_PA )
        {
            workspace->Dinv_L = smalloc( sizeof(real) * system->N_cm_max,
                    __FILE__, __LINE__ );
            workspace->Dinv_U = smalloc( sizeof(real) * system->N_cm_max,
                    __FILE__, __LINE__ );
            workspace->Dinv_b = smalloc( sizeof(real) * system->N_cm_max,
                    __FILE__, __LINE__ );
            workspace->rp = smalloc( sizeof(real) * system->N_cm_max,
                    __FILE__, __LINE__ );
            workspace->rp2 = smalloc( sizeof(real) * system->N_cm_max,
                    __FILE__, __LINE__ );
        }
        else
        {
            workspace->Dinv_L = NULL;
            workspace->Dinv_U = NULL;
            workspace->Dinv_b = NULL;
            workspace->rp = NULL;
            workspace->rp2 = NULL;
        }

        /* ILUTP preconditioner related */
        if ( control->cm_solver_pre_comp_type == ILUTP_PC )
        {
            workspace->perm_ilutp = smalloc( sizeof( int ) * system->N_cm_max,
                   __FILE__, __LINE__ );
        }
        else
        {
            workspace->perm_ilutp = NULL;
        }

        /* integrator storage */
        workspace->a = smalloc( system->N_max * sizeof( rvec ),
               __FILE__, __LINE__ );
        workspace->f_old = smalloc( system->N_max * sizeof( rvec ),
               __FILE__, __LINE__ );
        workspace->v_const = smalloc( system->N_max * sizeof( rvec ),
               __FILE__, __LINE__ );

#if defined(_OPENMP)
        workspace->f_local = smalloc( control->num_threads * system->N_max * sizeof( rvec ),
               __FILE__, __LINE__ );
#endif

        /* storage for analysis */
        if ( control->molec_anal || control->diffusion_coef )
        {
            workspace->mark = scalloc( system->N_max, sizeof(int),
                    __FILE__, __LINE__ );
            workspace->old_mark = scalloc( system->N_max, sizeof(int),
                    __FILE__, __LINE__ );
        }
        else
        {
            workspace->mark = workspace->old_mark = NULL;
        }

        if ( control->diffusion_coef )
        {
            workspace->x_old = scalloc( system->N_max, sizeof( rvec ),
                    __FILE__, __LINE__ );
        }
        else
        {
            workspace->x_old = NULL;
        }
    }

#if defined(TEST_FORCES)
    workspace->dDelta = smalloc( system->N_max * sizeof( rvec ),
           __FILE__, __LINE__ );
    workspace->f_ele = smalloc( system->N_max * sizeof( rvec ),
           __FILE__, __LINE__ );
    workspace->f_vdw = smalloc( system->N_max * sizeof( rvec ),
           __FILE__, __LINE__ );
    workspace->f_be = smalloc( system->N_max * sizeof( rvec ),
           __FILE__, __LINE__ );
    workspace->f_lp = smalloc( system->N_max * sizeof( rvec ),
           __FILE__, __LINE__ );
    workspace->f_ov = smalloc( system->N_max * sizeof( rvec ),
           __FILE__, __LINE__ );
    workspace->f_un = smalloc( system->N_max * sizeof( rvec ),
           __FILE__, __LINE__ );
    workspace->f_ang = smalloc( system->N_max * sizeof( rvec ),
           __FILE__, __LINE__ );
    workspace->f_coa = smalloc( system->N_max * sizeof( rvec ),
           __FILE__, __LINE__ );
    workspace->f_pen = smalloc( system->N_max * sizeof( rvec ),
           __FILE__, __LINE__ );
    workspace->f_hb = smalloc( system->N_max * sizeof( rvec ),
           __FILE__, __LINE__ );
    workspace->f_tor = smalloc( system->N_max * sizeof( rvec ),
           __FILE__, __LINE__ );
    workspace->f_con = smalloc( system->N_max * sizeof( rvec ),
           __FILE__, __LINE__ );
#endif

    workspace->realloc.far_nbrs = FALSE;
    workspace->realloc.cm = FALSE;
    workspace->realloc.hbonds = FALSE;
    workspace->realloc.bonds = FALSE;
    workspace->realloc.thbody = FALSE;
    workspace->realloc.gcell_atoms = -1;

    Reset_Workspace( system, workspace );

    Init_Taper( control, workspace );
}


static void Init_Lists( reax_system * const system,
        control_params * const control, simulation_data * const data,
        static_storage * const workspace, reax_list ** const lists, int realloc )
{
    int i, ret;

    if ( realloc == TRUE )
    {
        Estimate_Num_Neighbors( system, control, workspace, lists );

        if ( lists[FAR_NBRS]->allocated == FALSE )
        {
            Make_List( system->N, system->N_max, workspace->realloc.total_far_nbrs,
                    TYP_FAR_NEIGHBOR, lists[FAR_NBRS] );
        }
        else if ( realloc == TRUE
                || lists[FAR_NBRS]->total_intrs < workspace->realloc.total_far_nbrs )
        {
            if ( lists[FAR_NBRS]->allocated == TRUE )
            {
                Delete_List( TYP_FAR_NEIGHBOR, lists[FAR_NBRS] );
            }
            Make_List( system->N, system->N_max, 
                    MAX( workspace->realloc.total_far_nbrs, lists[FAR_NBRS]->total_intrs ),
                    TYP_FAR_NEIGHBOR, lists[FAR_NBRS] );
        }
    }

    lists[FAR_NBRS]->n = system->N;

    ret = Generate_Neighbor_Lists( system, control, data, workspace, lists );

    if ( ret != SUCCESS )
    {
        Estimate_Num_Neighbors( system, control, workspace, lists );

        Delete_List( TYP_FAR_NEIGHBOR, lists[FAR_NBRS] );
        Make_List( system->N, system->N_max, 
                MAX( workspace->realloc.total_far_nbrs, lists[FAR_NBRS]->total_intrs ),
                TYP_FAR_NEIGHBOR, lists[FAR_NBRS] );

        ret = Generate_Neighbor_Lists( system, control, data, workspace, lists );

        if ( ret != SUCCESS )
        {
            fprintf( stderr, "[ERROR] Unrecoverable memory allocation issue (Generate_Neighbor_Lists). Terminating...\n" );
            exit( INVALID_INPUT );
        }
    }

    if ( realloc == TRUE )
    {
        Estimate_Storages( system, control, workspace, lists, TRUE, TRUE );
    }

    if ( workspace->H.allocated == FALSE )
    {
        Allocate_Matrix( &workspace->H, system->N_cm, system->N_cm_max,
                workspace->realloc.total_cm_entries );
    }
    else if ( realloc == TRUE || workspace->H.m < workspace->realloc.total_cm_entries
            || workspace->H.n_max < system->N_cm_max )
    {
        if ( workspace->H.allocated == TRUE )
        {
            Deallocate_Matrix( &workspace->H );
        }
        Allocate_Matrix( &workspace->H, system->N_cm, system->N_cm_max,
                workspace->realloc.total_cm_entries );
    }
    else
    {
        workspace->H.n = system->N_cm;
    }

    if ( workspace->H_sp.allocated == FALSE )
    {
        /* TODO: better estimate for H_sp?
         *   If so, need to refactor Estimate_Storages
         *   to use various cut-off distances as parameters
         *   (non-bonded, hydrogen, 3body, etc.) */
        Allocate_Matrix( &workspace->H_sp, system->N_cm, system->N_cm_max,
                workspace->realloc.total_cm_entries );
    }
    else if ( realloc == TRUE || workspace->H_sp.m < workspace->realloc.total_cm_entries
            || workspace->H.n_max < system->N_cm_max )
    {
        if ( workspace->H_sp.allocated == TRUE )
        {
            Deallocate_Matrix( &workspace->H_sp );
        }
        /* TODO: better estimate for H_sp?
         *   If so, need to refactor Estimate_Storages
         *   to use various cut-off distances as parameters
         *   (non-bonded, hydrogen, 3body, etc.) */
        Allocate_Matrix( &workspace->H_sp, system->N_cm, system->N_cm_max,
                workspace->realloc.total_cm_entries );
    }
    else
    {
        workspace->H_sp.n = system->N_cm;
    }

    workspace->num_H = 0;
    if ( control->hbond_cut > 0.0 )
    {
        for ( i = 0; i < system->N; ++i )
        {
            if ( system->reax_param.sbp[ system->atoms[i].type ].p_hbond == H_ATOM )
            {
                ++(workspace->num_H);
            }
        }
    }

    if ( control->hbond_cut > 0.0 && workspace->num_H > 0 )
    {
        if ( lists[HBONDS]->allocated == FALSE )
        {
            Make_List( system->N, system->N_max,
                    (int) CEIL( SAFE_ZONE * workspace->realloc.total_hbonds ),
                    TYP_HBOND, lists[HBONDS] );
        }
        else if ( system->N_max < system->N
                || lists[HBONDS]->total_intrs < workspace->realloc.total_hbonds )
        {
            if ( lists[HBONDS]->allocated == TRUE )
            {
                Delete_List( TYP_HBOND, lists[HBONDS] );
            }
            Make_List( system->N, system->N_max,
                    MAX( workspace->realloc.total_hbonds, lists[HBONDS]->total_intrs ),
                    TYP_HBOND, lists[HBONDS] );
        }
        else
        {
            lists[HBONDS]->n = system->N;
        }
    }

    /* bonds list */
    if ( lists[BONDS]->allocated == FALSE )
    {
        Make_List( system->N, system->N_max, (int) CEIL( workspace->realloc.total_bonds * SAFE_ZONE ),
                TYP_BOND, lists[BONDS] );
    }
    else if ( realloc == TRUE || lists[BONDS]->total_intrs < workspace->realloc.total_bonds )
    {
        if ( lists[BONDS]->allocated == TRUE )
        {
            Delete_List( TYP_BOND, lists[BONDS] );
        }
        Make_List( system->N, system->N_max,
                MAX( workspace->realloc.total_bonds, lists[BONDS]->total_intrs ),
                TYP_BOND, lists[BONDS] );
    }
    else
    {
        lists[BONDS]->n = system->N;
    }

    /* 3bodies list */
    if ( lists[THREE_BODIES]->allocated == FALSE )
    {
        Make_List( workspace->realloc.total_bonds, (int) CEIL( workspace->realloc.total_bonds * SAFE_ZONE ),
                workspace->realloc.total_thbodies, TYP_THREE_BODY, lists[THREE_BODIES] );
    }
    else if ( lists[THREE_BODIES]->n_max < workspace->realloc.total_bonds
            || lists[THREE_BODIES]->total_intrs < workspace->realloc.total_thbodies )
    {
        if ( lists[THREE_BODIES]->allocated == TRUE )
        {
            Delete_List( TYP_THREE_BODY, lists[THREE_BODIES] );
        }
        Make_List( workspace->realloc.total_bonds,
                MAX( workspace->realloc.total_bonds, lists[THREE_BODIES]->n_max),
                MAX( workspace->realloc.total_thbodies, lists[THREE_BODIES]->total_intrs ),
                TYP_THREE_BODY, lists[THREE_BODIES] );
    }
    else
    {
        lists[THREE_BODIES]->n = workspace->realloc.total_bonds;
    }

#if defined(TEST_FORCES)
    //TODO: increased num. of DDELTA list elements, find a better count later
    Make_List( system->N, workspace->realloc.total_bonds * 20, TYP_DDELTA, lists[DDELTA] );

    for ( i = 0; i < lists[DDELTA]->n; ++i )
    {
        Set_Start_Index( i, 0, lists[DDELTA] );
        Set_End_Index( i, 0, lists[DDELTA] );
    }

    Make_List( workspace->realloc.total_bonds, workspace->realloc.total_bonds * MAX_BONDS * 3,
            TYP_DBO, lists[DBO] );

    for ( i = 0; i < lists[DBO]->n; ++i )
    {
        Set_Start_Index( i, 0, lists[DBO] );
        Set_End_Index( i, 0, lists[DBO] );
    }
#endif
}


static void Init_Out_Controls( reax_system *system, control_params *control,
        static_storage *workspace, output_controls *out_control, int output_enabled )
{
#define TEMP_SIZE (1000)
    char temp[TEMP_SIZE];

    out_control->allocated = TRUE;

    if ( output_enabled == TRUE && out_control->write_steps > 0 )
    {
        strncpy( temp, control->sim_name, TEMP_SIZE - 5 );
        temp[TEMP_SIZE - 5] = '\0';
        strcat( temp, ".trj" );
        out_control->trj = sfopen( temp, "w", __FILE__, __LINE__ );
        out_control->write_header( system, control, workspace, out_control );
    }
    else
    {
        out_control->trj = NULL;
    }

    if ( output_enabled == TRUE && out_control->log_update_freq > 0 )
    {
        strncpy( temp, control->sim_name, TEMP_SIZE - 5 );
        temp[TEMP_SIZE - 5] = '\0';
        strcat( temp, ".out" );
        out_control->out = sfopen( temp, "w", __FILE__, __LINE__ );
        fprintf( out_control->out, "%-6s%16s%16s%16s%11s%11s%13s%13s%13s\n",
                 "step", "total_energy", "poten_energy", "kin_energy",
                 "temp", "target", "volume", "press", "target" );
        fflush( out_control->out );

        strncpy( temp, control->sim_name, TEMP_SIZE - 5 );
        temp[TEMP_SIZE - 5] = '\0';
        strcat( temp, ".pot" );
        out_control->pot = sfopen( temp, "w", __FILE__, __LINE__ );
        fprintf( out_control->pot,
                 "%-6s%13s%13s%13s%13s%13s%13s%13s%13s%13s%13s%13s\n",
                 "step", "ebond", "eatom", "elp", "eang", "ecoa", "ehb",
                 "etor", "econj", "evdw", "ecoul", "epol" );
        fflush( out_control->pot );

        strncpy( temp, control->sim_name, TEMP_SIZE - 5 );
        temp[TEMP_SIZE - 5] = '\0';
        strcat( temp, ".log" );
        out_control->log = sfopen( temp, "w", __FILE__, __LINE__ );
        fprintf( out_control->log, "%-6s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n",
                 "step", "total", "neighbors", "init", "bonded",
                 "nonbonded", "cm", "cm_sort", "s_iters", "pre_comp", "pre_app",
                 "s_spmv", "s_vec_ops", "s_orthog", "s_tsolve" );
    }
    else
    {
        out_control->out = NULL;
        out_control->pot = NULL;
        out_control->log = NULL;
    }

    if ( output_enabled == TRUE && (control->ensemble == sNPT || control->ensemble == iNPT
            || control->ensemble == aNPT || control->compute_pressure == TRUE) )
    {
        strncpy( temp, control->sim_name, TEMP_SIZE - 5 );
        temp[TEMP_SIZE - 5] = '\0';
        strcat( temp, ".prs" );
        out_control->prs = sfopen( temp, "w", __FILE__, __LINE__ );
#if defined(DEBUG) || defined(DEBUG_FOCUS)
        fprintf( out_control->prs, "%-8s %13s %13s %13s %13s %13s %13s\n",
                "step", "KExx", "KEyy", "KEzz",
                "Virialxx", "Virialyy", "Virialzz" );
#endif
        fprintf( out_control->prs, "%-8s %13s %13s %13s %13s %13s %13s %13s %13s\n",
                "step", "Lx", "Ly", "Lz",
                "Pxx", "Pyy", "Pzz", "Pavg", "Volume" );

        fflush( out_control->prs );
    }
    else
    {
        out_control->prs = NULL;
    }

    /* Init molecular analysis file */
    if ( output_enabled == TRUE && control->molec_anal )
    {
        snprintf( temp, TEMP_SIZE, "%.*s.mol", TEMP_SIZE - 5, control->sim_name );
        out_control->mol = sfopen( temp, "w", __FILE__, __LINE__ );
        if ( control->num_ignored )
        {
            snprintf( temp, TEMP_SIZE, "%.*s.ign", TEMP_SIZE - 5, control->sim_name );
            out_control->ign = sfopen( temp, "w", __FILE__, __LINE__ );
        }
    }
    else
    {
        out_control->mol = NULL;
        out_control->ign = NULL;
    }

    if ( output_enabled == TRUE && control->dipole_anal )
    {
        strncpy( temp, control->sim_name, TEMP_SIZE - 5 );
        temp[TEMP_SIZE - 5] = '\0';
        strcat( temp, ".dpl" );
        out_control->dpl = sfopen( temp, "w", __FILE__, __LINE__ );
        fprintf( out_control->dpl,
                 "Step      Molecule Count  Avg. Dipole Moment Norm\n" );
        fflush( out_control->dpl );
    }
    else
    {
        out_control->dpl = NULL;
    }

    if ( output_enabled == TRUE && control->diffusion_coef )
    {
        strncpy( temp, control->sim_name, TEMP_SIZE - 6 );
        temp[TEMP_SIZE - 6] = '\0';
        strcat( temp, ".drft" );
        out_control->drft = sfopen( temp, "w", __FILE__, __LINE__ );
        fprintf( out_control->drft, "Step     Type Count   Avg Squared Disp\n" );
        fflush( out_control->drft );
    }
    else
    {
        out_control->drft = NULL;
    }


#if defined(TEST_ENERGY)
    strncpy( temp, control->sim_name, TEMP_SIZE - 6 );
    temp[TEMP_SIZE - 6] = '\0';
    strcat( temp, ".ebond" );
    out_control->ebond = sfopen( temp, "w", __FILE__, __LINE__ );

    strncpy( temp, control->sim_name, TEMP_SIZE - 5 );
    temp[TEMP_SIZE - 5] = '\0';
    strcat( temp, ".elp" );
    out_control->elp = sfopen( temp, "w", __FILE__, __LINE__ );

    strncpy( temp, control->sim_name, TEMP_SIZE - 5 );
    temp[TEMP_SIZE - 5] = '\0';
    strcat( temp, ".eov" );
    out_control->eov = sfopen( temp, "w", __FILE__, __LINE__ );

    strncpy( temp, control->sim_name, TEMP_SIZE - 5 );
    temp[TEMP_SIZE - 5] = '\0';
    strcat( temp, ".eun" );
    out_control->eun = sfopen( temp, "w", __FILE__, __LINE__ );

    strncpy( temp, control->sim_name, TEMP_SIZE - 6 );
    temp[TEMP_SIZE - 6] = '\0';
    strcat( temp, ".eval" );
    out_control->eval = sfopen( temp, "w", __FILE__, __LINE__ );

    strncpy( temp, control->sim_name, TEMP_SIZE - 6 );
    temp[TEMP_SIZE - 6] = '\0';
    strcat( temp, ".epen" );
    out_control->epen = sfopen( temp, "w", __FILE__, __LINE__ );

    strncpy( temp, control->sim_name, TEMP_SIZE - 6 );
    temp[TEMP_SIZE - 6] = '\0';
    strcat( temp, ".ecoa" );
    out_control->ecoa = sfopen( temp, "w", __FILE__, __LINE__ );

    strncpy( temp, control->sim_name, TEMP_SIZE - 5 );
    temp[TEMP_SIZE - 5] = '\0';
    strcat( temp, ".ehb" );
    out_control->ehb = sfopen( temp, "w", __FILE__, __LINE__ );

    strncpy( temp, control->sim_name, TEMP_SIZE - 6 );
    temp[TEMP_SIZE - 6] = '\0';
    strcat( temp, ".etor" );
    out_control->etor = sfopen( temp, "w", __FILE__, __LINE__ );

    strncpy( temp, control->sim_name, TEMP_SIZE - 6 );
    temp[TEMP_SIZE - 6] = '\0';
    strcat( temp, ".econ" );
    out_control->econ = sfopen( temp, "w", __FILE__, __LINE__ );

    strncpy( temp, control->sim_name, TEMP_SIZE - 6 );
    temp[TEMP_SIZE - 6] = '\0';
    strcat( temp, ".evdw" );
    out_control->evdw = sfopen( temp, "w", __FILE__, __LINE__ );

    strncpy( temp, control->sim_name, TEMP_SIZE - 6 );
    temp[TEMP_SIZE - 6] = '\0';
    strcat( temp, ".ecou" );
    out_control->ecou = sfopen( temp, "w", __FILE__, __LINE__ );
#endif


#if defined(TEST_FORCES)
    /* open bond orders file */
    strncpy( temp, control->sim_name, TEMP_SIZE - 5 );
    temp[TEMP_SIZE - 5] = '\0';
    strcat( temp, ".fbo" );
    out_control->fbo = sfopen( temp, "w", __FILE__, __LINE__ );

    /* open bond orders derivatives file */
    strncpy( temp, control->sim_name, TEMP_SIZE - 6 );
    temp[TEMP_SIZE - 6] = '\0';
    strcat( temp, ".fdbo" );
    out_control->fdbo = sfopen( temp, "w", __FILE__, __LINE__ );

    /* open bond forces file */
    strncpy( temp, control->sim_name, TEMP_SIZE - 7 );
    temp[TEMP_SIZE - 7] = '\0';
    strcat( temp, ".fbond" );
    out_control->fbond = sfopen( temp, "w", __FILE__, __LINE__ );

    /* open lone-pair forces file */
    strncpy( temp, control->sim_name, TEMP_SIZE - 6 );
    temp[TEMP_SIZE - 6] = '\0';
    strcat( temp, ".flp" );
    out_control->flp = sfopen( temp, "w", __FILE__, __LINE__ );

    /* open overcoordination forces file */
    strncpy( temp, control->sim_name, TEMP_SIZE - 7 );
    temp[TEMP_SIZE - 7] = '\0';
    strcat( temp, ".fatom" );
    out_control->fatom = sfopen( temp, "w", __FILE__, __LINE__ );

    /* open angle forces file */
    strncpy( temp, control->sim_name, TEMP_SIZE - 8 );
    temp[TEMP_SIZE - 8] = '\0';
    strcat( temp, ".f3body" );
    out_control->f3body = sfopen( temp, "w", __FILE__, __LINE__ );

    /* open hydrogen bond forces file */
    strncpy( temp, control->sim_name, TEMP_SIZE - 5 );
    temp[TEMP_SIZE - 5] = '\0';
    strcat( temp, ".fhb" );
    out_control->fhb = sfopen( temp, "w", __FILE__, __LINE__ );

    /* open torsion forces file */
    strncpy( temp, control->sim_name, TEMP_SIZE - 8 );
    temp[TEMP_SIZE - 8] = '\0';
    strcat( temp, ".f4body" );
    out_control->f4body = sfopen( temp, "w", __FILE__, __LINE__ );

    /* open nonbonded forces file */
    strncpy( temp, control->sim_name, TEMP_SIZE - 7 );
    temp[TEMP_SIZE - 7] = '\0';
    strcat( temp, ".fnonb" );
    out_control->fnonb = sfopen( temp, "w", __FILE__, __LINE__ );

    /* open total force file */
    strncpy( temp, control->sim_name, TEMP_SIZE - 6 );
    temp[TEMP_SIZE - 6] = '\0';
    strcat( temp, ".ftot" );
    out_control->ftot = sfopen( temp, "w", __FILE__, __LINE__ );

    /* open coulomb forces file */
    strncpy( temp, control->sim_name, TEMP_SIZE - 7 );
    temp[TEMP_SIZE - 7] = '\0';
    strcat( temp, ".ftot2" );
    out_control->ftot2 = sfopen( temp, "w", __FILE__, __LINE__ );
#endif

#undef TEMP_SIZE
}


void Initialize( reax_system * const system, control_params * const control,
        simulation_data * const data, static_storage * const workspace,
        reax_list ** const lists, output_controls * const out_control,
        evolve_function * const Evolve, int output_enabled, int realloc )
{
#if defined(_OPENMP)
    #pragma omp parallel default(none) shared(control)
    {
        #pragma omp single
        {
            if ( control->num_threads_set == FALSE )
            {
                /* set using OMP_NUM_THREADS environment variable */
                control->num_threads = omp_get_num_threads( );
                control->num_threads_set = TRUE;
            }
        }
    }

    omp_set_num_threads( control->num_threads );
#else
    control->num_threads = 1;
#endif

    Randomize( );

    Init_System( system, control, data, workspace, realloc );

    Init_Simulation_Data( system, control, data, Evolve, realloc );

    Init_Workspace( system, control, workspace, realloc );

    Init_Lists( system, control, data, workspace, lists, realloc );

    Init_Out_Controls( system, control, workspace, out_control, output_enabled );

    /* These are done in forces.c, only forces.c can see all those functions */
    Init_Bonded_Force_Functions( control );

#if defined(TEST_FORCES)
    Init_Force_Test_Functions( control );
#endif

    if ( control->tabulate )
    {
        Make_LR_Lookup_Table( system, control, workspace );
    }
}


static void Finalize_System( reax_system *system, control_params *control,
        simulation_data *data, int reset )
{
    int i, j, k;
    reax_interaction *reax;

    system->prealloc_allocated = FALSE;
    system->ffield_params_allocated = FALSE;

    if ( system->max_num_molec_charge_constraints > 0 )
    {
        sfree( system->molec_charge_constraints, __FILE__, __LINE__ );
        sfree( system->molec_charge_constraint_ranges, __FILE__, __LINE__ );
    }

    if ( system->max_num_custom_charge_constraints > 0 )
    {
        sfree( system->custom_charge_constraint_count, __FILE__, __LINE__ );
        sfree( system->custom_charge_constraint_start, __FILE__, __LINE__ );
        sfree( system->custom_charge_constraint_rhs, __FILE__, __LINE__ );
    }

    if ( system->max_num_custom_charge_constraint_entries > 0 )
    {
        sfree( system->custom_charge_constraint_atom_index, __FILE__, __LINE__ );
        sfree( system->custom_charge_constraint_coeff, __FILE__, __LINE__ );
    }

    system->max_num_molec_charge_constraints = 0;
    system->num_molec_charge_constraints = 0;
    system->max_num_custom_charge_constraints = 0;
    system->num_custom_charge_constraints = 0;
    system->max_num_custom_charge_constraint_entries = 0;
    system->num_custom_charge_constraint_entries = 0;

    reax = &system->reax_param;

    Finalize_Grid( system );

    if ( reset == FALSE )
    {
        sfree( reax->gp.l, __FILE__, __LINE__ );

        for ( i = 0; i < reax->max_num_atom_types; i++ )
        {
            for ( j = 0; j < reax->max_num_atom_types; j++ )
            {
                for ( k = 0; k < reax->max_num_atom_types; k++ )
                {
                    sfree( reax->fbp[i][j][k], __FILE__, __LINE__ );
                }

                sfree( reax->thbp[i][j], __FILE__, __LINE__ );
                sfree( reax->hbp[i][j], __FILE__, __LINE__ );
                sfree( reax->fbp[i][j], __FILE__, __LINE__ );
            }

            sfree( reax->tbp[i], __FILE__, __LINE__ );
            sfree( reax->thbp[i], __FILE__, __LINE__ );
            sfree( reax->hbp[i], __FILE__, __LINE__ );
            sfree( reax->fbp[i], __FILE__, __LINE__ );
        }

        sfree( reax->sbp, __FILE__, __LINE__ );
        sfree( reax->tbp, __FILE__, __LINE__ );
        sfree( reax->thbp, __FILE__, __LINE__ );
        sfree( reax->hbp, __FILE__, __LINE__ );
        sfree( reax->fbp, __FILE__, __LINE__ );

        sfree( system->atoms, __FILE__, __LINE__ );
    }

    if ( system->allocated == TRUE )
    {
        sfree( system->bonds, __FILE__, __LINE__ );
        sfree( system->hbonds, __FILE__, __LINE__ );
    }

    system->allocated = FALSE;
}


static void Finalize_Simulation_Data( reax_system *system, control_params *control,
        simulation_data *data, output_controls *out_control )
{
#if defined(_OPENMP)
    if ( control->ensemble == sNPT || control->ensemble == iNPT
            || control->ensemble == aNPT || control->compute_pressure == TRUE )
    {
        sfree( data->press_local, __FILE__, __LINE__ );
    }
#endif
}


static void Finalize_Workspace( reax_system *system, control_params *control,
        static_storage *workspace, int reset )
{
    int i;

    if ( workspace->allocated == TRUE )
    {
        sfree( workspace->total_bond_order, __FILE__, __LINE__ );
        sfree( workspace->Deltap, __FILE__, __LINE__ );
        sfree( workspace->Deltap_boc, __FILE__, __LINE__ );
        sfree( workspace->dDeltap_self, __FILE__, __LINE__ );
        sfree( workspace->Delta, __FILE__, __LINE__ );
        sfree( workspace->Delta_lp, __FILE__, __LINE__ );
        sfree( workspace->Delta_lp_temp, __FILE__, __LINE__ );
        sfree( workspace->dDelta_lp, __FILE__, __LINE__ );
        sfree( workspace->dDelta_lp_temp, __FILE__, __LINE__ );
        sfree( workspace->Delta_e, __FILE__, __LINE__ );
        sfree( workspace->Delta_boc, __FILE__, __LINE__ );
        sfree( workspace->nlp, __FILE__, __LINE__ );
        sfree( workspace->nlp_temp, __FILE__, __LINE__ );
        sfree( workspace->Clp, __FILE__, __LINE__ );
        sfree( workspace->CdDelta, __FILE__, __LINE__ );
        sfree( workspace->vlpex, __FILE__, __LINE__ );

        if ( workspace->H.allocated == TRUE )
        {
            Deallocate_Matrix( &workspace->H );
        }
        if ( workspace->H_full.allocated == TRUE )
        {
            Deallocate_Matrix( &workspace->H_full );
        }
        if ( workspace->H_sp.allocated == TRUE )
        {
            Deallocate_Matrix( &workspace->H_sp );
        }
        if ( workspace->H_p.allocated == TRUE )
        {
            Deallocate_Matrix( &workspace->H_p );
        }
        if ( workspace->H_spar_patt.allocated == TRUE )
        {
            Deallocate_Matrix( &workspace->H_spar_patt );
        }
        if ( workspace->H_spar_patt_full.allocated == TRUE )
        {
            Deallocate_Matrix( &workspace->H_spar_patt_full );
        }
        if ( workspace->H_app_inv.allocated == TRUE )
        {
            Deallocate_Matrix( &workspace->H_app_inv );
        }
        if ( workspace->L.allocated == TRUE )
        {
            Deallocate_Matrix( &workspace->L );
        }
        if ( workspace->U.allocated == TRUE )
        {
            Deallocate_Matrix( &workspace->U );
        }

        for ( i = 0; i < 5; ++i )
        {
            sfree( workspace->s[i], __FILE__, __LINE__ );
            sfree( workspace->t[i], __FILE__, __LINE__ );
        }

        if ( control->cm_solver_pre_comp_type == JACOBI_PC )
        {
            sfree( workspace->Hdia_inv, __FILE__, __LINE__ );
        }
        if ( control->cm_solver_pre_comp_type == ICHOLT_PC
                || (control->cm_solver_pre_comp_type == ILUT_PC && control->cm_solver_pre_comp_droptol > 0.0 )
                || control->cm_solver_pre_comp_type == ILUTP_PC
                || control->cm_solver_pre_comp_type == FG_ILUT_PC )
        {
            sfree( workspace->droptol, __FILE__, __LINE__ );
        }
        sfree( workspace->b_s, __FILE__, __LINE__ );
        sfree( workspace->b_t, __FILE__, __LINE__ );
        sfree( workspace->b_prc, __FILE__, __LINE__ );
        sfree( workspace->b_prm, __FILE__, __LINE__ );
        sfree( workspace->s, __FILE__, __LINE__ );
        sfree( workspace->t, __FILE__, __LINE__ );

        switch ( control->cm_solver_type )
        {
            case GMRES_S:
            case GMRES_H_S:
                for ( i = 0; i < control->cm_solver_restart + 1; ++i )
                {
                    sfree( workspace->h[i], __FILE__, __LINE__ );
                    sfree( workspace->rn[i], __FILE__, __LINE__ );
                    sfree( workspace->v[i], __FILE__, __LINE__ );
                }

                sfree( workspace->y, __FILE__, __LINE__ );
                sfree( workspace->z, __FILE__, __LINE__ );
                sfree( workspace->g, __FILE__, __LINE__ );
                sfree( workspace->h, __FILE__, __LINE__ );
                sfree( workspace->hs, __FILE__, __LINE__ );
                sfree( workspace->hc, __FILE__, __LINE__ );
                sfree( workspace->rn, __FILE__, __LINE__ );
                sfree( workspace->v, __FILE__, __LINE__ );

                sfree( workspace->r, __FILE__, __LINE__ );
                sfree( workspace->d, __FILE__, __LINE__ );
                sfree( workspace->q, __FILE__, __LINE__ );
                sfree( workspace->p, __FILE__, __LINE__ );
                break;

            case CG_S:
                sfree( workspace->r, __FILE__, __LINE__ );
                sfree( workspace->d, __FILE__, __LINE__ );
                sfree( workspace->q, __FILE__, __LINE__ );
                sfree( workspace->p, __FILE__, __LINE__ );
                break;

            case SDM_S:
                sfree( workspace->r, __FILE__, __LINE__ );
                sfree( workspace->d, __FILE__, __LINE__ );
                sfree( workspace->q, __FILE__, __LINE__ );
                break;

            case BiCGStab_S:
                sfree( workspace->r, __FILE__, __LINE__ );
                sfree( workspace->r_hat, __FILE__, __LINE__ );
                sfree( workspace->d, __FILE__, __LINE__ );
                sfree( workspace->q, __FILE__, __LINE__ );
                sfree( workspace->q_hat, __FILE__, __LINE__ );
                sfree( workspace->p, __FILE__, __LINE__ );
                sfree( workspace->y, __FILE__, __LINE__ );
                sfree( workspace->z, __FILE__, __LINE__ );
                sfree( workspace->g, __FILE__, __LINE__ );
                break;

            default:
                fprintf( stderr, "[ERROR] Unknown charge method linear solver type. Terminating...\n" );
                exit( INVALID_INPUT );
                break;
        }

        /* SpMV related */
#if defined(_OPENMP)
        sfree( workspace->b_local, __FILE__, __LINE__ );
#endif

        /* level scheduling related */
        if ( control->cm_solver_pre_app_type == TRI_SOLVE_LEVEL_SCHED_PA ||
                control->cm_solver_pre_app_type == TRI_SOLVE_GC_PA )
        {
            sfree( workspace->row_levels_L, __FILE__, __LINE__ );
            sfree( workspace->level_rows_L, __FILE__, __LINE__ );
            sfree( workspace->level_rows_cnt_L, __FILE__, __LINE__ );
            sfree( workspace->row_levels_U, __FILE__, __LINE__ );
            sfree( workspace->level_rows_U, __FILE__, __LINE__ );
            sfree( workspace->level_rows_cnt_U, __FILE__, __LINE__ );
            sfree( workspace->top, __FILE__, __LINE__ );
        }

        /* graph coloring related */
        if ( control->cm_solver_pre_app_type == TRI_SOLVE_GC_PA )
        {
            sfree( workspace->color, __FILE__, __LINE__ );
            sfree( workspace->to_color, __FILE__, __LINE__ );
            sfree( workspace->conflict, __FILE__, __LINE__ );
            sfree( workspace->conflict_cnt, __FILE__, __LINE__ );
            sfree( workspace->recolor, __FILE__, __LINE__ );
            sfree( workspace->color_top, __FILE__, __LINE__ );
            sfree( workspace->permuted_row_col, __FILE__, __LINE__ );
            sfree( workspace->permuted_row_col_inv, __FILE__, __LINE__ );
        }

        /* graph coloring related OR ILUTP preconditioner */
        if ( control->cm_solver_pre_app_type == TRI_SOLVE_GC_PA 
                || control->cm_solver_pre_comp_type == ILUTP_PC )
        {
            sfree( workspace->y_p, __FILE__, __LINE__ );
            sfree( workspace->x_p, __FILE__, __LINE__ );
        }

        /* Jacobi iteration related */
        if ( control->cm_solver_pre_app_type == JACOBI_ITER_PA )
        {
            sfree( workspace->Dinv_L, __FILE__, __LINE__ );
            sfree( workspace->Dinv_U, __FILE__, __LINE__ );
            sfree( workspace->Dinv_b, __FILE__, __LINE__ );
            sfree( workspace->rp, __FILE__, __LINE__ );
            sfree( workspace->rp2, __FILE__, __LINE__ );
        }

        /* ILUTP preconditioner related */
        if ( control->cm_solver_pre_comp_type == ILUTP_PC )
        {
            sfree( workspace->perm_ilutp, __FILE__, __LINE__ );
        }

        /* integrator storage */
        sfree( workspace->a, __FILE__, __LINE__ );
        sfree( workspace->f_old, __FILE__, __LINE__ );
        sfree( workspace->v_const, __FILE__, __LINE__ );

#if defined(_OPENMP)
        sfree( workspace->f_local, __FILE__, __LINE__ );
#endif

        /* storage for analysis */
        if ( control->molec_anal || control->diffusion_coef )
        {
            sfree( workspace->mark, __FILE__, __LINE__ );
            sfree( workspace->old_mark, __FILE__, __LINE__ );
        }

        if ( control->diffusion_coef )
        {
            sfree( workspace->x_old, __FILE__, __LINE__ );
        }
    }

    if ( reset == FALSE && (control->geo_format == BGF
            || control->geo_format == ASCII_RESTART
            || control->geo_format == BINARY_RESTART) )
    {
        sfree( workspace->map_serials, __FILE__, __LINE__ );
    }

    if ( reset == FALSE )
    {
        sfree( workspace->orig_id, __FILE__, __LINE__ );

        /* space for keeping restriction info, if any */
        if ( control->restrict_bonds )
        {
            for ( i = 0; i < system->N; ++i )
            {
                sfree( workspace->restricted_list[i],
                        __FILE__, __LINE__ );
            }

            sfree( workspace->restricted, __FILE__, __LINE__ );
            sfree( workspace->restricted_list, __FILE__, __LINE__ );
        }
    }

#if defined(TEST_FORCES)
    sfree( workspace->dDelta, __FILE__, __LINE__ );
    sfree( workspace->f_ele, __FILE__, __LINE__ );
    sfree( workspace->f_vdw, __FILE__, __LINE__ );
    sfree( workspace->f_be, __FILE__, __LINE__ );
    sfree( workspace->f_lp, __FILE__, __LINE__ );
    sfree( workspace->f_ov, __FILE__, __LINE__ );
    sfree( workspace->f_un, __FILE__, __LINE__ );
    sfree( workspace->f_ang, __FILE__, __LINE__ );
    sfree( workspace->f_coa, __FILE__, __LINE__ );
    sfree( workspace->f_pen, __FILE__, __LINE__ );
    sfree( workspace->f_hb, __FILE__, __LINE__ );
    sfree( workspace->f_tor, __FILE__, __LINE__ );
    sfree( workspace->f_con, __FILE__, __LINE__ );
#endif

    workspace->allocated = FALSE;
}


static void Finalize_Lists( reax_list **lists )
{
    if ( lists[FAR_NBRS]->allocated == TRUE )
    {
        Delete_List( TYP_FAR_NEIGHBOR, lists[FAR_NBRS] );
    }
    if ( lists[HBONDS]->allocated == TRUE )
    {
        Delete_List( TYP_HBOND, lists[HBONDS] );
    }
    if ( lists[BONDS]->allocated == TRUE )
    {
        Delete_List( TYP_BOND, lists[BONDS] );
    }
    if ( lists[OLD_BONDS]->allocated == TRUE )
    {
        Delete_List( TYP_BOND, lists[OLD_BONDS] );
    }
    if ( lists[THREE_BODIES]->allocated == TRUE )
    {
        Delete_List( TYP_THREE_BODY, lists[THREE_BODIES] );
    }

#if defined(TEST_FORCES)
    if ( lists[DDELTA]->allocated == TRUE )
    {
        Delete_List( TYP_DDELTA, lists[DDELTA] );
    }
    if ( lists[DBO]->allocated == TRUE )
    {
        Delete_List( TYP_DBO, lists[DBO] );
    }
#endif
}


void Finalize_Out_Controls( reax_system *system, control_params *control,
        static_storage *workspace, output_controls *out_control )
{
    if ( out_control->allocated == TRUE )
    {
        if ( out_control->write_steps > 0 )
        {
            sfclose( out_control->trj, __FILE__, __LINE__ );
        }

        if ( out_control->log_update_freq > 0 )
        {
            sfclose( out_control->out, __FILE__, __LINE__ );
            sfclose( out_control->pot, __FILE__, __LINE__ );
            sfclose( out_control->log, __FILE__, __LINE__ );
        }

        if ( control->ensemble == sNPT || control->ensemble == iNPT
                || control->ensemble == aNPT || control->compute_pressure == TRUE )
        {
            sfclose( out_control->prs, __FILE__, __LINE__ );
        }

        if ( control->molec_anal )
        {
            sfclose( out_control->mol, __FILE__, __LINE__ );

            if ( control->num_ignored )
            {
                sfclose( out_control->ign, __FILE__, __LINE__ );
            }
        }

        if ( control->dipole_anal )
        {
            sfclose( out_control->dpl, __FILE__, __LINE__ );
        }

        if ( control->diffusion_coef )
        {
            sfclose( out_control->drft, __FILE__, __LINE__ );
        }

#if defined(TEST_ENERGY)
        sfclose( out_control->ebond, __FILE__, __LINE__ );
        sfclose( out_control->elp, __FILE__, __LINE__ );
        sfclose( out_control->eov, __FILE__, __LINE__ );
        sfclose( out_control->eun, __FILE__, __LINE__ );
        sfclose( out_control->eval, __FILE__, __LINE__ );
        sfclose( out_control->epen, __FILE__, __LINE__ );
        sfclose( out_control->ecoa, __FILE__, __LINE__ );
        sfclose( out_control->ehb, __FILE__, __LINE__ );
        sfclose( out_control->etor, __FILE__, __LINE__ );
        sfclose( out_control->econ, __FILE__, __LINE__ );
        sfclose( out_control->evdw, __FILE__, __LINE__ );
        sfclose( out_control->ecou, __FILE__, __LINE__ );
#endif

#if defined(TEST_FORCES)
        sfclose( out_control->fbo, __FILE__, __LINE__ );
        sfclose( out_control->fdbo, __FILE__, __LINE__ );
        sfclose( out_control->fbond, __FILE__, __LINE__ );
        sfclose( out_control->flp, __FILE__, __LINE__ );
        sfclose( out_control->fatom, __FILE__, __LINE__ );
        sfclose( out_control->f3body, __FILE__, __LINE__ );
        sfclose( out_control->fhb, __FILE__, __LINE__ );
        sfclose( out_control->f4body, __FILE__, __LINE__ );
        sfclose( out_control->fnonb, __FILE__, __LINE__ );
        sfclose( out_control->ftot, __FILE__, __LINE__ );
        sfclose( out_control->ftot2, __FILE__, __LINE__ );
#endif
    }

    out_control->allocated = FALSE;
}


/* Deallocate top-level data structures, close file handles, etc.
 *
 */
void Finalize( reax_system *system, control_params *control,
        simulation_data *data, static_storage *workspace, reax_list **lists,
        output_controls *out_control, int output_enabled, int reset )
{
    if ( control->tabulate )
    {
        Finalize_LR_Lookup_Table( system, control, workspace );
    }

    if ( output_enabled == TRUE && reset == FALSE )
    {
        Finalize_Out_Controls( system, control, workspace, out_control );
    }

    Finalize_Lists( lists );

    Finalize_Workspace( system, control, workspace, reset );

    Finalize_Simulation_Data( system, control, data, out_control );

    Finalize_System( system, control, data, reset );
}
