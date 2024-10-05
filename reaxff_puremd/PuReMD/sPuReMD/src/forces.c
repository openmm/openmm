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

#include "forces.h"

#include "box.h"
#include "bond_orders.h"
#include "bonds.h"
#include "charges.h"
#include "hydrogen_bonds.h"
#if defined(TEST_FORCES)
  #include "io_tools.h"
#endif
#include "list.h"
#include "multi_body.h"
#include "nonbonded.h"
#include "system_props.h"
#include "torsion_angles.h"
#include "tool_box.h"
#include "valence_angles.h"
#include "vector.h"


typedef enum
{
    DIAGONAL = 0,
    OFF_DIAGONAL = 1,
} MATRIX_ENTRY_POSITION;


#if defined(TEST_FORCES)
static int compare_bonds( const void *p1, const void *p2 )
{
    return ((bond_data *)p1)->nbr - ((bond_data *)p2)->nbr;
}
#endif


void Init_Bonded_Force_Functions( control_params * const control )
{
    control->intr_funcs[0] = &BO;
    control->intr_funcs[1] = &Bonds;
    control->intr_funcs[2] = &Atom_Energy;
    control->intr_funcs[3] = &Valence_Angles;
    control->intr_funcs[4] = &Torsion_Angles;
    if ( control->hbond_cut > 0.0 )
    {
        control->intr_funcs[5] = &Hydrogen_Bonds;
    }
    else
    {
        control->intr_funcs[5] = NULL;
    }
    control->intr_funcs[6] = NULL;
    control->intr_funcs[7] = NULL;
    control->intr_funcs[8] = NULL;
    control->intr_funcs[9] = NULL;
}


static void Compute_Bonded_Forces( reax_system *system, control_params *control,
        simulation_data *data, static_storage *workspace, reax_list **lists,
        output_controls *out_control )
{
    int i;

#if defined(TEST_ENERGY)
    /* Mark beginning of a new timestep in each energy file */
    fprintf( out_control->ebond, "step: %d\n%6s%6s%12s%12s%12s\n",
             data->step, "atom1", "atom2", "bo", "ebond", "total" );
    fprintf( out_control->elp, "step: %d\n%6s%12s%12s%12s\n",
             data->step, "atom", "nlp", "elp", "total" );
    fprintf( out_control->eov, "step: %d\n%6s%12s%12s\n",
             data->step, "atom", "eov", "total" );
    fprintf( out_control->eun, "step: %d\n%6s%12s%12s\n",
             data->step, "atom", "eun", "total" );
    fprintf( out_control->eval, "step: %d\n%6s%6s%6s%12s%12s%12s%12s%12s%12s\n",
             data->step, "atom1", "atom2", "atom3",
             "angle", "bo(12)", "bo(23)", "eval", "epen", "total" );
    fprintf( out_control->epen, "step: %d\n%6s%6s%6s%12s%12s%12s%12s%12s\n",
             data->step, "atom1", "atom2", "atom3",
             "angle", "bo(12)", "bo(23)", "epen", "total" );
    fprintf( out_control->ecoa, "step: %d\n%6s%6s%6s%12s%12s%12s%12s%12s\n",
             data->step, "atom1", "atom2", "atom3",
             "angle", "bo(12)", "bo(23)", "ecoa", "total" );
    fprintf( out_control->ehb,  "step: %d\n%6s%6s%6s%12s%12s%12s%12s%12s\n",
             data->step, "atom1", "atom2", "atom3",
             "r(23)", "angle", "bo(12)", "ehb", "total" );
    fprintf( out_control->etor, "step: %d\n%6s%6s%6s%6s%12s%12s%12s%12s\n",
             data->step, "atom1", "atom2", "atom3", "atom4",
             "phi", "bo(23)", "etor", "total" );
    fprintf( out_control->econ, "step:%d\n%6s%6s%6s%6s%12s%12s%12s%12s%12s%12s\n",
             data->step, "atom1", "atom2", "atom3", "atom4",
             "phi", "bo(12)", "bo(23)", "bo(34)", "econ", "total" );
#endif

    /* function calls for bonded interactions */
    for ( i = 0; i < NUM_INTRS; i++ )
    {
        if ( control->intr_funcs[i] != NULL )
        {
            (control->intr_funcs[i])( system, control, data, workspace,
                    lists, out_control );
        }
    }

#if defined(TEST_FORCES)
    /* function calls for printing bonded interactions */
    for ( i = 0; i < NUM_INTRS; i++ )
    {
        if ( control->print_intr_funcs[i] != NULL )
        {
            (control->print_intr_funcs[i])( system, control, data, workspace,
                    lists, out_control );
        }
    }
#endif
}


static void Compute_NonBonded_Forces( reax_system *system, control_params *control,
        simulation_data *data, static_storage *workspace,
        reax_list** lists, output_controls *out_control, int realloc )
{
#if defined(TEST_ENERGY)
    fprintf( out_control->evdw, "step: %d\n%6s%6s%12s%12s%12s\n",
             data->step, "atom1", "atom2", "r12", "evdw", "total" );
    fprintf( out_control->ecou, "step: %d\n%6s%6s%12s%12s%12s%12s%12s\n",
             data->step, "atom1", "atom2", "r12", "q1", "q2", "ecou", "total" );
#endif

    if ( control->tabulate <= 0 )
    {
        vdW_Coulomb_Energy( system, control, data, workspace, lists, out_control );
    }
    else
    {
        Tabulated_vdW_Coulomb_Energy( system, control, data, workspace,
                lists, out_control );
    }

#if defined(TEST_FORCES)
    Print_vdW_Coulomb_Forces( system, control, data, workspace,
            lists, out_control );
#endif
}


/* This version of Compute_Total_Force computes forces from coefficients
   accumulated by all interaction functions. Saves enormous time & space! */
static void Compute_Total_Force( reax_system *system, control_params *control,
        simulation_data *data, static_storage *workspace, reax_list **lists )
{
    int i;
    reax_list *bonds;

    bonds = lists[BONDS];

#if defined(_OPENMP)
    #pragma omp parallel default(shared)
#endif
    {
        int pj;
#if defined(_OPENMP)
        int j;
#endif

        if ( control->compute_pressure == FALSE
                && (control->ensemble == NVE || control->ensemble == nhNVT
                    || control->ensemble == bNVT) )
        {
#if defined(_OPENMP)
            #pragma omp for schedule(static)
#endif
            for ( i = 0; i < system->N; ++i )
            {
                for ( pj = Start_Index(i, bonds); pj < End_Index(i, bonds); ++pj )
                {
                    if ( i <= bonds->bond_list[pj].nbr )
                    {
                        Add_dBond_to_Forces( i, pj, system, data, workspace, lists );
                    }
                }
            }
        }
        else if ( control->ensemble == sNPT || control->ensemble == iNPT
                || control->ensemble == aNPT || control->compute_pressure == TRUE )
        {
#if defined(_OPENMP)
            #pragma omp for schedule(static)
#endif
            for ( i = 0; i < system->N; ++i )
            {
                for ( pj = Start_Index(i, bonds); pj < End_Index(i, bonds); ++pj )
                {
                    if ( i <= bonds->bond_list[pj].nbr )
                    {
                        Add_dBond_to_Forces_NPT( i, pj, system, data, workspace, lists );
                    }
                }
            }
        }

#if defined(_OPENMP)
        /* reduction (sum) on thread-local force vectors */
        #pragma omp for schedule(static)
        for ( i = 0; i < system->N; ++i )
        {
            for ( j = 0; j < control->num_threads; ++j )
            {
                rvec_Add( system->atoms[i].f, workspace->f_local[j * system->N + i] );
            }
        }
#endif
    }
}


static inline real Init_Charge_Matrix_Entry_Tab( reax_system const * const system,
        control_params const * const control, LR_lookup_table ** const LR, int i, int j,
        real r_ij, MATRIX_ENTRY_POSITION pos )
{
    int r;
    real base, dif, val, ret;
    LR_lookup_table *t;

    ret = 0.0;

    switch ( control->charge_method )
    {
    case QEQ_CM:
    //TODO: tabulate other portions of matrices for EE, ACKS2?
    case EE_CM:
    case ACKS2_CM:
        switch ( pos )
        {
            case OFF_DIAGONAL:
                t = &LR[MIN( system->atoms[i].type, system->atoms[j].type )]
                       [MAX( system->atoms[i].type, system->atoms[j].type )];

                /* cubic spline interpolation */
                r = (int)(r_ij * t->inv_dx);
                if ( r == 0 ) 
                {
                    ++r;
                }
                base = (real)(r + 1) * t->dx;
                dif = r_ij - base;
                val = ((t->ele[r].d * dif + t->ele[r].c) * dif + t->ele[r].b)
                    * dif + t->ele[r].a;
                val *= EV_to_KCALpMOL / C_ELE;

                ret = ((i == j) ? 0.5 : 1.0) * val;
            break;

            case DIAGONAL:
                ret = system->reax_param.sbp[system->atoms[i].type].eta;
            break;

            default:
                fprintf( stderr, "[ERROR] Invalid matrix position. Terminating...\n" );
                exit( INVALID_INPUT );
            break;
        }
        break;

    default:
        fprintf( stderr, "[ERROR] Invalid charge method. Terminating...\n" );
        exit( INVALID_INPUT );
        break;
    }

    return ret;
}


static inline real Init_Charge_Matrix_Entry( reax_system const * const system,
        control_params const * const control, static_storage const * const workspace,
        int i, int j, real r_ij, MATRIX_ENTRY_POSITION pos )
{
    real Tap, dr3gamij_1, dr3gamij_3, ret;

    ret = 0.0;

    switch ( control->charge_method )
    {
    case QEQ_CM:
    case EE_CM:
    case ACKS2_CM:
        switch ( pos )
        {
            case OFF_DIAGONAL:
                Tap = workspace->Tap[7] * r_ij + workspace->Tap[6];
                Tap = Tap * r_ij + workspace->Tap[5];
                Tap = Tap * r_ij + workspace->Tap[4];
                Tap = Tap * r_ij + workspace->Tap[3];
                Tap = Tap * r_ij + workspace->Tap[2];
                Tap = Tap * r_ij + workspace->Tap[1];
                Tap = Tap * r_ij + workspace->Tap[0];

                /* shielding */
                dr3gamij_1 = r_ij * r_ij * r_ij
                        + POW( system->reax_param.tbp[system->atoms[i].type][system->atoms[j].type].gamma, -3.0 );
                dr3gamij_3 = POW( dr3gamij_1 , 1.0 / 3.0 );

                /* i == j: periodic self-interaction term
                 * i != j: general interaction term */
                ret = ((i == j) ? 0.5 : 1.0) * Tap * EV_to_KCALpMOL / dr3gamij_3;
            break;

            case DIAGONAL:
                ret = system->reax_param.sbp[system->atoms[i].type].eta;
            break;

            default:
                fprintf( stderr, "[ERROR] Invalid matrix position. Terminating...\n" );
                exit( INVALID_INPUT );
            break;
        }
        break;

    default:
        fprintf( stderr, "[ERROR] Invalid charge method. Terminating...\n" );
        exit( INVALID_INPUT );
        break;
    }

    return ret;
}


static int Init_Charge_Matrix_Remaining_Entries( reax_system const * const system,
        control_params const * const control, reax_list const * const far_nbr_list,
        sparse_matrix * const H, sparse_matrix * const H_sp,
        int * const Htop, int * const H_sp_top )
{
    int i, j, pj, target, val_flag, flag_oom;
    real d, xcut, bond_softness, * X_diag;

    flag_oom = FALSE;

    switch ( control->charge_method )
    {
        case QEQ_CM:
            break;

        case EE_CM:
            if ( system->num_molec_charge_constraints == 0
                    && system->num_custom_charge_constraints == 0 )
            {
                H->start[system->N] = *Htop;
                H_sp->start[system->N] = *H_sp_top;

                for ( i = 0; i < system->N; ++i )
                {
                    /* total charge constraint on atoms */
                    if ( *Htop < H->m )
                    {
                        H->j[*Htop] = i;
                        H->val[*Htop] = 1.0;
                        ++(*Htop);
                    }
                    else
                    {
                        flag_oom = TRUE;
                        break;
                    }

                    if ( *H_sp_top < H_sp->m )
                    {
                        H_sp->j[*H_sp_top] = i;
                        H_sp->val[*H_sp_top] = 1.0;
                        ++(*H_sp_top);
                    }
                    else
                    {
                        flag_oom = TRUE;
                        break;
                    }
                }

                if ( *Htop < H->m )
                {
                    H->j[*Htop] = system->N;
                    H->val[*Htop] = 0.0;
                    ++(*Htop);
                }
                else
                {
                    flag_oom = TRUE;
                }

                if ( *H_sp_top < H_sp->m )
                {
                    H_sp->j[*H_sp_top] = system->N;
                    H_sp->val[*H_sp_top] = 0.0;
                    ++(*H_sp_top);
                }
                else
                {
                    flag_oom = TRUE;
                }
            }
            else
            {
                for ( i = 0; i < system->num_molec_charge_constraints; ++i )
                {
                    H->start[system->N + i] = *Htop;
                    H_sp->start[system->N + i] = *H_sp_top;

                    for ( j = system->molec_charge_constraint_ranges[2 * i];
                            j <= system->molec_charge_constraint_ranges[2 * i + 1]; ++j )
                    {
                        /* molecule charge constraint on atoms */
                        if ( *Htop < H->m )
                        {
                            H->j[*Htop] = j - 1;
                            H->val[*Htop] = 1.0;
                            ++(*Htop);
                        }
                        else
                        {
                            flag_oom = TRUE;
                            break;
                        }

                        if ( *H_sp_top < H_sp->m )
                        {
                            H_sp->j[*H_sp_top] = j - 1;
                            H_sp->val[*H_sp_top] = 1.0;
                            ++(*H_sp_top);
                        }
                        else
                        {
                            flag_oom = TRUE;
                            break;
                        }
                    }

                    /* explicit zeros on diagonals */
                    if ( *Htop < H->m )
                    {
                        H->j[*Htop] = system->N + i;
                        H->val[*Htop] = 0.0; 
                        ++(*Htop);
                    }
                    else
                    {
                        flag_oom = TRUE;
                    }

                    if ( *H_sp_top < H_sp->m )
                    {
                        H_sp->j[*H_sp_top] = system->N + i;
                        H_sp->val[*H_sp_top] = 0.0;
                        ++(*H_sp_top);
                    }
                    else
                    {
                        flag_oom = TRUE;
                    }
                }

                for ( i = system->num_molec_charge_constraints;
                        i < system->num_molec_charge_constraints + system->num_custom_charge_constraints; ++i )
                {
                    if ( flag_oom == FALSE )
                    {
                        H->start[system->N + i] = *Htop;
                        H_sp->start[system->N + i] = *H_sp_top;

                        for ( j = system->custom_charge_constraint_start[i - system->num_molec_charge_constraints];
                                j < system->custom_charge_constraint_start[i - system->num_molec_charge_constraints + 1]; ++j )
                        {
                            /* custom charge constraint on atoms */
                            if ( *Htop < H->m )
                            {
                                H->j[*Htop] = system->custom_charge_constraint_atom_index[j] - 1;
                                H->val[*Htop] = system->custom_charge_constraint_coeff[j];
                                ++(*Htop);
                            }
                            else
                            {
                                flag_oom = TRUE;
                                break;
                            }

                            if ( *H_sp_top < H_sp->m )
                            {
                                H_sp->j[*H_sp_top] = system->custom_charge_constraint_atom_index[j] - 1;
                                H_sp->val[*H_sp_top] = system->custom_charge_constraint_coeff[j];
                                ++(*H_sp_top);
                            }
                            else
                            {
                                flag_oom = TRUE;
                                break;
                            }
                        }

                        /* explicit zeros on diagonals */
                        if ( *Htop < H->m )
                        {
                            H->j[*Htop] = system->N + i;
                            H->val[*Htop] = 0.0; 
                            ++(*Htop);
                        }
                        else
                        {
                            flag_oom = TRUE;
                        }

                        if ( *H_sp_top < H_sp->m )
                        {
                            H_sp->j[*H_sp_top] = system->N + i;
                            H_sp->val[*H_sp_top] = 0.0;
                            ++(*H_sp_top);
                        }
                        else
                        {
                            flag_oom = TRUE;
                        }
                    }
                }
            }
            break;

        case ACKS2_CM:
            X_diag = smalloc( sizeof(real) * system->N, __FILE__, __LINE__ );

            for ( i = 0; i < system->N; ++i )
            {
                X_diag[i] = 0.0;
            }

            for ( i = 0; i < system->N; ++i )
            {
                if ( flag_oom == FALSE )
                {
                    H->start[system->N + i] = *Htop;
                    H_sp->start[system->N + i] = *H_sp_top;

                    /* constraint on ref. value for kinetic energy potential */
                    if ( *Htop < H->m )
                    {
                        H->j[*Htop] = i;
                        H->val[*Htop] = 1.0;
                        ++(*Htop);
                    }
                    else
                    {
                        flag_oom = TRUE;
                    }

                    if ( *H_sp_top < H_sp->m )
                    {
                        H_sp->j[*H_sp_top] = i;
                        H_sp->val[*H_sp_top] = 1.0;
                        ++(*H_sp_top);
                    }
                    else
                    {
                        flag_oom = TRUE;
                    }

                    /* kinetic energy terms */
                    for ( pj = Start_Index(i, far_nbr_list); pj < End_Index(i, far_nbr_list); ++pj )
                    {
                        /* exclude self-periodic images of atoms for
                         * kinetic energy term because the effective
                         * potential is the same on an atom and its periodic image */
                        if ( far_nbr_list->far_nbr_list[pj].d <= control->nonb_cut )
                        {
                            j = far_nbr_list->far_nbr_list[pj].nbr;

                            xcut = 0.5 * ( system->reax_param.sbp[ system->atoms[i].type ].b_s_acks2
                                    + system->reax_param.sbp[ system->atoms[j].type ].b_s_acks2 );

                            if ( far_nbr_list->far_nbr_list[pj].d < xcut )
                            {
                                d = far_nbr_list->far_nbr_list[pj].d / xcut;
                                bond_softness = system->reax_param.gp.l[34] * POW( d, 3.0 )
                                    * POW( 1.0 - d, 6.0 );

                                if ( bond_softness > 0.0 )
                                {
                                    val_flag = FALSE;

                                    for ( target = H->start[system->N + i]; target < *Htop; ++target )
                                    {
                                        if ( H->j[target] == system->N + j )
                                        {
                                            H->val[target] += bond_softness;
                                            val_flag = TRUE;
                                            break;
                                        }
                                    }

                                    if ( val_flag == FALSE )
                                    {
                                        if ( *Htop < H->m )
                                        {
                                            H->j[*Htop] = system->N + j;
                                            H->val[*Htop] = bond_softness;
                                            ++(*Htop);
                                        }
                                        else
                                        {
                                            flag_oom = TRUE;
                                            break;
                                        }
                                    }

                                    val_flag = FALSE;

                                    for ( target = H_sp->start[system->N + i]; target < *H_sp_top; ++target )
                                    {
                                        if ( H_sp->j[target] == system->N + j )
                                        {
                                            H_sp->val[target] += bond_softness;
                                            val_flag = TRUE;
                                            break;
                                        }
                                    }

                                    if ( val_flag == FALSE )
                                    {
                                        if ( *H_sp_top < H_sp->m )
                                        {
                                            H_sp->j[*H_sp_top] = system->N + j;
                                            H_sp->val[*H_sp_top] = bond_softness;
                                            ++(*H_sp_top);
                                        }
                                        else
                                        {
                                            flag_oom = TRUE;
                                            break;
                                        }
                                    }

                                    X_diag[i] -= bond_softness;
                                    X_diag[j] -= bond_softness;
                                }
                            }
                        }
                    }

                    /* placeholders for diagonal entries, to be replaced below */
                    if ( *Htop < H->m )
                    {
                        H->j[*Htop] = system->N + i;
                        H->val[*Htop] = 0.0;
                        ++(*Htop);
                    }
                    else
                    {
                        flag_oom = TRUE;
                    }

                    if ( *H_sp_top < H_sp->m )
                    {
                        H_sp->j[*H_sp_top] = system->N + i;
                        H_sp->val[*H_sp_top] = 0.0;
                        ++(*H_sp_top);
                    }
                    else
                    {
                        flag_oom = TRUE;
                    }
                }
            }

            /* second to last row */
            H->start[system->N_cm - 2] = *Htop;
            H_sp->start[system->N_cm - 2] = *H_sp_top;

            /* place accumulated diagonal entries (needed second to last row marker above before this code) */
            for ( i = system->N; i < 2 * system->N; ++i )
            {
                for ( pj = H->start[i]; pj < H->start[i + 1]; ++pj )
                {
                    if ( H->j[pj] == i )
                    {
                        H->val[pj] = X_diag[i - system->N];
                        break;
                    }
                }

                for ( pj = H_sp->start[i]; pj < H_sp->start[i + 1]; ++pj )
                {
                    if ( H_sp->j[pj] == i )
                    {
                        H_sp->val[pj] = X_diag[i - system->N];
                        break;
                    }
                }
            }

            /* coupling with the kinetic energy potential */
            for ( i = 0; i < system->N; ++i )
            {
                if ( *Htop < H->m )
                {
                    H->j[*Htop] = system->N + i;
                    H->val[*Htop] = 1.0;
                    ++(*Htop);
                }
                else
                {
                    flag_oom = TRUE;
                    break;
                }

                if ( *H_sp_top < H_sp->m )
                {
                    H_sp->j[*H_sp_top] = system->N + i;
                    H_sp->val[*H_sp_top] = 1.0;
                    ++(*H_sp_top);
                }
                else
                {
                    flag_oom = TRUE;
                    break;
                }
            }

            /* explicitly store zero on diagonal */
            if ( *Htop < H->m )
            {
                H->j[*Htop] = system->N_cm - 2;
                H->val[*Htop] = 0.0;
                ++(*Htop);
            }
            else
            {
                flag_oom = TRUE;
            }

            if ( *H_sp_top < H_sp->m )
            {
                H_sp->j[*H_sp_top] = system->N_cm - 2;
                H_sp->val[*H_sp_top] = 0.0;
                ++(*H_sp_top);
            }
            else
            {
                flag_oom = TRUE;
            }

            /* last row */
            H->start[system->N_cm - 1] = *Htop;
            H_sp->start[system->N_cm - 1] = *H_sp_top;

            for ( i = 0; i < system->N; ++i )
            {
                if ( *Htop < H->m )
                {
                    H->j[*Htop] = i;
                    H->val[*Htop] = 1.0;
                    ++(*Htop);
                }
                else
                {
                    flag_oom = TRUE;
                    break;
                }

                if ( *H_sp_top < H_sp->m )
                {
                    H_sp->j[*H_sp_top] = i;
                    H_sp->val[*H_sp_top] = 1.0;
                    ++(*H_sp_top);
                }
                else
                {
                    flag_oom = TRUE;
                    break;
                }
            }

            /* explicitly store zero on diagonal */
            if ( *Htop < H->m )
            {
                H->j[*Htop] = system->N_cm - 1;
                H->val[*Htop] = 0.0;
                ++(*Htop);
            }
            else
            {
                flag_oom = TRUE;
            }

            if ( *H_sp_top < H_sp->m )
            {
                H_sp->j[*H_sp_top] = system->N_cm - 1;
                H_sp->val[*H_sp_top] = 0.0;
                ++(*H_sp_top);
            }
            else
            {
                flag_oom = TRUE;
            }

            sfree( X_diag, __FILE__, __LINE__ );
            break;

        default:
            break;
    }

    return flag_oom;
}


/* Compute the distances and displacement vectors for entries
 * in the far neighbors list if it's a NOT re-neighboring step */
static void Init_Distance( reax_system const * const system,
        control_params const * const control, reax_list ** const lists )
{
    int i, j, pj;
    int start_i, end_i;
    reax_list *far_nbrs;

    far_nbrs = lists[FAR_NBRS];

    for ( i = 0; i < far_nbrs->n; ++i )
    {
        start_i = Start_Index( i, far_nbrs );
        end_i = End_Index( i, far_nbrs );

        /* update distance and displacement vector between atoms i and j (i-j)
         * for the j atom entry in the far nbr list */
        for ( pj = start_i; pj < end_i; ++pj )
        {
            j = far_nbrs->far_nbr_list[pj].nbr;

            far_nbrs->far_nbr_list[pj].d = control->compute_atom_distance(
                    &system->box, system->atoms[i].x, system->atoms[j].x,
                    system->atoms[i].rel_map, system->atoms[j].rel_map,
                    far_nbrs->far_nbr_list[pj].rel_box,
                    far_nbrs->far_nbr_list[pj].dvec );
        }
    }
}


/* Compute the charge matrix entries and store the matrix in half format
 * using the far neighbors list (stored in half format)
 */
static void Init_CM_Half( reax_system const * const system,
        control_params const * const control,
        static_storage * const workspace, reax_list ** const lists )
{
    int i, j, pj, target;
    int start_i, end_i;
    int Htop, H_sp_top;
    int flag, flag_sp, val_flag, flag_oom;
    real val;
    sparse_matrix *H, *H_sp;
    reax_list *far_nbrs;

    far_nbrs = lists[FAR_NBRS];
    H = &workspace->H;
    H_sp = &workspace->H_sp;
    Htop = 0;
    H_sp_top = 0;
    flag_oom = FALSE;

    for ( i = 0; i < far_nbrs->n; ++i )
    {
        if ( flag_oom == FALSE )
        {
            start_i = Start_Index( i, far_nbrs );
            end_i = End_Index( i, far_nbrs );
            H->start[i] = Htop;
            H_sp->start[i] = H_sp_top;

            for ( pj = start_i; pj < end_i; ++pj )
            {
                j = far_nbrs->far_nbr_list[pj].nbr;
                flag = FALSE;
                flag_sp = FALSE;

                if ( far_nbrs->far_nbr_list[pj].d <= control->nonb_cut )
                {
                    flag = TRUE;

                    if ( far_nbrs->far_nbr_list[pj].d <= control->nonb_sp_cut )
                    {
                        flag_sp = TRUE;
                    }
                }

                if ( flag == TRUE )
                {
                    val = Init_Charge_Matrix_Entry( system, control,
                                workspace, i, j, far_nbrs->far_nbr_list[pj].d,
                                OFF_DIAGONAL );
                    val_flag = FALSE;

                    for ( target = H->start[i]; target < Htop; ++target )
                    {
                        if ( H->j[target] == j )
                        {
                            H->val[target] += val;
                            val_flag = TRUE;
                            break;
                        }
                    }

                    if ( val_flag == FALSE )
                    {
                        if ( Htop < H->m )
                        {
                            H->j[Htop] = j;
                            H->val[Htop] = val;
                            ++Htop;
                        }
                        else
                        {
                            flag_oom = TRUE;
                            break;
                        }
                    }

                    /* H_sp matrix entry */
                    if ( flag_sp == TRUE )
                    {
                        val_flag = FALSE;

                        for ( target = H_sp->start[i]; target < H_sp_top; ++target )
                        {
                            if ( H_sp->j[target] == j )
                            {
                                H_sp->val[target] += val;
                                val_flag = TRUE;
                                break;
                            }
                        }

                        if ( val_flag == FALSE )
                        {
                            if ( H_sp_top < H_sp->m )
                            {
                                H_sp->j[H_sp_top] = j;
                                H_sp->val[H_sp_top] = val;
                                ++H_sp_top;
                            }
                            else
                            {
                                flag_oom = TRUE;
                                break;
                            }
                        }
                    }
                }
            }

            /* diagonal entry */
            if ( Htop < H->m )
            {
                H->j[Htop] = i;
                H->val[Htop] = Init_Charge_Matrix_Entry( system, control,
                        workspace, i, i, 0.0, DIAGONAL );
                ++Htop;
            }
            else
            {
                flag_oom = TRUE;
            }

            if ( H_sp_top < H_sp->m )
            {
                H_sp->j[H_sp_top] = i;
                H_sp->val[H_sp_top] = H->val[Htop - 1];
                ++H_sp_top;
            }
            else
            {
                flag_oom = TRUE;
            }
        }
    }

    if ( flag_oom == FALSE )
    {
        flag_oom = Init_Charge_Matrix_Remaining_Entries( system, control, far_nbrs,
                H, H_sp, &Htop, &H_sp_top );
    }

    H->start[system->N_cm] = Htop;
    H_sp->start[system->N_cm] = H_sp_top;
    
    if ( flag_oom == TRUE )
    {
        workspace->realloc.cm = TRUE;
    }
}


/* Compute the tabulated charge matrix entries and store the matrix in half format
 * using the far neighbors list (stored in half format)
 */
static void Init_CM_Tab_Half( reax_system const * const system,
        control_params const * const control,
        static_storage * const workspace, reax_list ** const lists )
{
    int i, j, pj, target;
    int start_i, end_i;
    int Htop, H_sp_top;
    int flag, flag_sp, val_flag, flag_oom;
    real val;
    sparse_matrix *H, *H_sp;
    reax_list *far_nbrs;

    far_nbrs = lists[FAR_NBRS];
    H = &workspace->H;
    H_sp = &workspace->H_sp;
    Htop = 0;
    H_sp_top = 0;
    flag_oom = FALSE;

    for ( i = 0; i < far_nbrs->n; ++i )
    {
        if ( flag_oom == FALSE )
        {
            start_i = Start_Index( i, far_nbrs );
            end_i = End_Index( i, far_nbrs );
            H->start[i] = Htop;
            H_sp->start[i] = H_sp_top;

            for ( pj = start_i; pj < end_i; ++pj )
            {
                j = far_nbrs->far_nbr_list[pj].nbr;
                flag = FALSE;
                flag_sp = FALSE;

                if ( far_nbrs->far_nbr_list[pj].d <= control->nonb_cut )
                {
                    flag = TRUE;

                    if ( far_nbrs->far_nbr_list[pj].d <= control->nonb_sp_cut )
                    {
                        flag_sp = TRUE;
                    }
                }

                if ( flag == TRUE )
                {
                    val = Init_Charge_Matrix_Entry_Tab( system, control,
                                workspace->LR, i, j, far_nbrs->far_nbr_list[pj].d,
                                OFF_DIAGONAL );
                    val_flag = FALSE;

                    for ( target = H->start[i]; target < Htop; ++target )
                    {
                        if ( H->j[target] == j )
                        {
                            H->val[target] += val;
                            val_flag = TRUE;
                            break;
                        }
                    }

                    if ( val_flag == FALSE )
                    {
                        if ( Htop < H->m )
                        {
                            H->j[Htop] = j;
                            H->val[Htop] = val;
                            ++Htop;
                        }
                        else
                        {
                            flag_oom = TRUE;
                            break;
                        }
                    }

                    /* H_sp matrix entry */
                    if ( flag_sp == TRUE )
                    {
                        val_flag = FALSE;

                        for ( target = H_sp->start[i]; target < H_sp_top; ++target )
                        {
                            if ( H_sp->j[target] == j )
                            {
                                H_sp->val[target] += val;
                                val_flag = TRUE;
                                break;
                            }
                        }

                        if ( val_flag == FALSE )
                        {
                            if ( H_sp_top < H_sp->m )
                            {
                                H_sp->j[H_sp_top] = j;
                                H_sp->val[H_sp_top] = val;
                                ++H_sp_top;
                            }
                            else
                            {
                                flag_oom = TRUE;
                                break;
                            }
                        }
                    }
                }
            }

            /* diagonal entry */
            if ( Htop < H->m )
            {
                H->j[Htop] = i;
                H->val[Htop] = Init_Charge_Matrix_Entry_Tab( system, control,
                        workspace->LR, i, i, 0.0, DIAGONAL );
                ++Htop;
            }
            else
            {
                flag_oom = TRUE;
            }

            if ( H_sp_top < H_sp->m )
            {
                H_sp->j[H_sp_top] = i;
                H_sp->val[H_sp_top] = H->val[Htop - 1];
                ++H_sp_top;
            }
            else
            {
                flag_oom = TRUE;
            }
        }
    }

    if ( flag_oom == FALSE )
    {
        flag_oom = Init_Charge_Matrix_Remaining_Entries( system, control, far_nbrs,
                H, H_sp, &Htop, &H_sp_top );

        H->start[system->N_cm] = Htop;
        H_sp->start[system->N_cm] = H_sp_top;
    }
    
    if ( flag_oom == TRUE )
    {
        workspace->realloc.cm = TRUE;
    }
}


/* Compute entries of the bonds/hbonds lists and store the lists in full format
 * using the far neighbors list (stored in full format) */
static void Init_Bond_Full( reax_system const * const system,
        control_params const * const control,
        static_storage * const workspace, reax_list ** const lists )
{
    int i, j, pj;
    int start_i, end_i;
    int type_i, type_j;
    int btop_i, btop_j;
    int ihb, jhb, ihb_top, jhb_top;
    int num_bonds, num_hbonds, flag_oom_bonds, flag_oom_hbonds;
    real r_ij, r2;
    real C12, C34, C56;
    real Cln_BOp_s, Cln_BOp_pi, Cln_BOp_pi2;
    real BO, BO_s, BO_pi, BO_pi2;
    reax_list *far_nbrs, *bonds, *hbonds;
    single_body_parameters *sbp_i, *sbp_j;
    two_body_parameters *twbp;
    far_neighbor_data *nbr_pj;
    bond_data *ibond, *jbond;
    bond_order_data *bo_ij, *bo_ji;

    far_nbrs = lists[FAR_NBRS];
    bonds = lists[BONDS];
    hbonds = lists[HBONDS];
    num_bonds = 0;
    num_hbonds = 0;
    btop_i = 0;
    btop_j = 0;
    flag_oom_bonds = FALSE;
    flag_oom_hbonds = FALSE;

    for ( i = 0; i < far_nbrs->n; ++i )
    {
        type_i = system->atoms[i].type;
        start_i = Start_Index( i, far_nbrs );
        end_i = End_Index( i, far_nbrs );
        btop_i = End_Index( i, bonds );
        sbp_i = &system->reax_param.sbp[type_i];

        if ( control->hbond_cut > 0.0 && workspace->num_H > 0 )
        {
            ihb = sbp_i->p_hbond;

            if ( ihb == H_ATOM )
            {
                ihb_top = End_Index( i, hbonds );
            }
            else
            {
                ihb_top = -1;
            }
        }
        else
        {
            ihb = NON_H_BONDING_ATOM;
            ihb_top = -1;
        }

        for ( pj = start_i; pj < end_i; ++pj )
        {
            nbr_pj = &far_nbrs->far_nbr_list[pj];
            j = nbr_pj->nbr;

#if defined(QMMM)
            if ( system->atoms[i].qmmm_mask == TRUE
                    || system->atoms[j].qmmm_mask == TRUE )
            {
#endif	
            if ( nbr_pj->d <= control->nonb_cut )
            {
                type_j = system->atoms[j].type;
                sbp_j = &system->reax_param.sbp[type_j];
                twbp = &system->reax_param.tbp[type_i][type_j];
                r_ij = nbr_pj->d;

#if defined(QMMM)
                if ( system->atoms[i].qmmm_mask == TRUE
                        && system->atoms[j].qmmm_mask == TRUE )
                {
#endif
                if ( system->atoms[i].is_dummy == FALSE
                        && system->atoms[j].is_dummy == FALSE )
                {

                    /* hydrogen bond lists */
                    if ( control->hbond_cut > 0.0 && workspace->num_H > 0
                            && (ihb == H_ATOM || ihb == H_BONDING_ATOM)
                            && nbr_pj->d <= control->hbond_cut )
                    {
                        jhb = sbp_j->p_hbond;

                        if ( ihb == H_ATOM && jhb == H_BONDING_ATOM )
                        {
                            if ( num_hbonds < hbonds->total_intrs )
                            {
                                hbonds->hbond_list[ihb_top].nbr = j;
                                hbonds->hbond_list[ihb_top].scl = 1;
                                hbonds->hbond_list[ihb_top].ptr = nbr_pj;
                                ++ihb_top;
                                ++num_hbonds;
                            }
                            else
                            {
                                flag_oom_hbonds = TRUE;
                            }
                        }
                        else if ( ihb == H_BONDING_ATOM && jhb == H_ATOM )
                        {
                            if ( num_hbonds < hbonds->total_intrs )
                            {
                                jhb_top = End_Index( j, hbonds );
                                hbonds->hbond_list[jhb_top].nbr = i;
                                hbonds->hbond_list[jhb_top].scl = -1;
                                hbonds->hbond_list[jhb_top].ptr = nbr_pj;
                                Set_End_Index( j, jhb_top + 1, hbonds );
                                ++num_hbonds;
                            }
                            else
                            {
                                flag_oom_hbonds = TRUE;
                            }
                        }
                    }

                    /* uncorrected bond orders */
                    if ( nbr_pj->d <= control->bond_cut )
                    {
                        r2 = SQR( r_ij );

                        if ( sbp_i->r_s > 0.0 && sbp_j->r_s > 0.0 )
                        {
                            C12 = twbp->p_bo1 * POW( r_ij / twbp->r_s, twbp->p_bo2 );
                            BO_s = (1.0 + control->bo_cut) * EXP( C12 );
                        }
                        else
                        {
                            C12 = 0.0;
                            BO_s = 0.0;
                        }

                        if ( sbp_i->r_pi > 0.0 && sbp_j->r_pi > 0.0 )
                        {
                            C34 = twbp->p_bo3 * POW( r_ij / twbp->r_p, twbp->p_bo4 );
                            BO_pi = EXP( C34 );
                        }
                        else
                        {
                            C34 = 0.0;
                            BO_pi = 0.0;
                        }

                        if ( sbp_i->r_pi_pi > 0.0 && sbp_j->r_pi_pi > 0.0 )
                        {
                            C56 = twbp->p_bo5 * POW( r_ij / twbp->r_pp, twbp->p_bo6 );
                            BO_pi2 = EXP( C56 );
                        }
                        else
                        {
                            C56 = 0.0;
                            BO_pi2 = 0.0;
                        }

                        /* Initially BO values are the uncorrected ones, page 1 */
                        BO = BO_s + BO_pi + BO_pi2;

                        if ( BO >= control->bo_cut )
                        {
                            if ( num_bonds < bonds->total_intrs )
                            {
                                /****** bonds i-j and j-i ******/
                                ibond = &bonds->bond_list[btop_i];
                                btop_j = End_Index( j, bonds );
                                jbond = &bonds->bond_list[btop_j];

                                ibond->nbr = j;
                                jbond->nbr = i;
                                ibond->d = r_ij;
                                jbond->d = r_ij;
                                rvec_Copy( ibond->dvec, nbr_pj->dvec );
                                rvec_Scale( jbond->dvec, -1, nbr_pj->dvec );
                                ivec_Copy( ibond->rel_box, nbr_pj->rel_box );
                                ivec_Scale( jbond->rel_box, -1, nbr_pj->rel_box );
                                ibond->dbond_index = btop_i;
                                jbond->dbond_index = btop_i;
                                ibond->sym_index = btop_j;
                                jbond->sym_index = btop_i;
                                ++btop_i;
                                Set_End_Index( j, btop_j + 1, bonds );

                                bo_ij = &ibond->bo_data;
                                bo_ij->BO = BO;
                                bo_ij->BO_s = BO_s;
                                bo_ij->BO_pi = BO_pi;
                                bo_ij->BO_pi2 = BO_pi2;
                                bo_ji = &jbond->bo_data;
                                bo_ji->BO = BO;
                                bo_ji->BO_s = BO_s;
                                bo_ji->BO_pi = BO_pi;
                                bo_ji->BO_pi2 = BO_pi2;

                                /* Bond Order page2-3, derivative of total bond order prime */
                                Cln_BOp_s = twbp->p_bo2 * C12 / r2;
                                Cln_BOp_pi = twbp->p_bo4 * C34 / r2;
                                Cln_BOp_pi2 = twbp->p_bo6 * C56 / r2;

                                /* Only dln_BOp_xx wrt. dr_i is stored here, note that
                                   dln_BOp_xx/dr_i = -dln_BOp_xx/dr_j and all others are 0 */
                                rvec_Scale( bo_ij->dln_BOp_s, -bo_ij->BO_s * Cln_BOp_s, ibond->dvec );
                                rvec_Scale( bo_ij->dln_BOp_pi, -bo_ij->BO_pi * Cln_BOp_pi, ibond->dvec );
                                rvec_Scale( bo_ij->dln_BOp_pi2,
                                        -bo_ij->BO_pi2 * Cln_BOp_pi2, ibond->dvec );
                                rvec_Scale( bo_ji->dln_BOp_s, -1., bo_ij->dln_BOp_s );
                                rvec_Scale( bo_ji->dln_BOp_pi, -1., bo_ij->dln_BOp_pi );
                                rvec_Scale( bo_ji->dln_BOp_pi2, -1., bo_ij->dln_BOp_pi2 );

                                /* Only dBOp wrt. dr_i is stored here, note that
                                   dBOp/dr_i = -dBOp/dr_j and all others are 0 */
                                rvec_Scale( bo_ij->dBOp, -(bo_ij->BO_s * Cln_BOp_s
                                            + bo_ij->BO_pi * Cln_BOp_pi
                                            + bo_ij->BO_pi2 * Cln_BOp_pi2), ibond->dvec );
                                rvec_Scale( bo_ji->dBOp, -1., bo_ij->dBOp );

                                rvec_Add( workspace->dDeltap_self[i], bo_ij->dBOp );
                                rvec_Add( workspace->dDeltap_self[j], bo_ji->dBOp );

                                bo_ij->BO_s -= control->bo_cut;
                                bo_ij->BO -= control->bo_cut;
                                bo_ji->BO_s -= control->bo_cut;
                                bo_ji->BO -= control->bo_cut;
                                workspace->total_bond_order[i] += bo_ij->BO;
                                workspace->total_bond_order[j] += bo_ji->BO;
                                bo_ij->Cdbo = 0.0;
                                bo_ij->Cdbopi = 0.0;
                                bo_ij->Cdbopi2 = 0.0;
                                bo_ji->Cdbo = 0.0;
                                bo_ji->Cdbopi = 0.0;
                                bo_ji->Cdbopi2 = 0.0;

                                Set_End_Index( j, btop_j + 1, bonds );
                                num_bonds += 2;
                            }
                            else
                            {
                                flag_oom_bonds = TRUE;
                            }
                        }
                    }
                }
#if defined(QMMM)
                }
#endif
            }
#if defined(QMMM)
            }
#endif
        }

        Set_End_Index( i, btop_i, bonds );
        if ( ihb == H_ATOM )
        {
            Set_End_Index( i, ihb_top, hbonds );
        }
    }

    if ( flag_oom_bonds == TRUE )
    {
        workspace->realloc.bonds = TRUE;
    }

    if ( control->hbond_cut > 0.0 && workspace->num_H > 0
            && flag_oom_hbonds == TRUE )
    {
        workspace->realloc.hbonds = TRUE;
    }
}


/* Generate bond list (full format), hydrogen bond list (full format),
 * and charge matrix (half symmetric format)
 * from the far neighbors list (with distance updates, if necessary)
 * */
static int Init_Forces( reax_system * const system,
        control_params const * const control,
        simulation_data * const data, static_storage * const workspace,
        reax_list ** const lists, output_controls * const out_control )
{
    int i, renbr, ret;
    static int dist_done = FALSE, cm_done = FALSE, bonds_done = FALSE;

    renbr = ((data->step - data->prev_steps) % control->reneighbor) == 0 ? TRUE : FALSE;

    if ( renbr == FALSE && dist_done == FALSE )
    {
        Init_Distance( system, control, lists );

        dist_done = TRUE;
    }

    if ( (control->charge_freq > 0
            && (data->step - data->prev_steps) % control->charge_freq == 0)
            && cm_done == FALSE )
    {
        if ( control->tabulate <= 0 )
        {
//            if ( workspace->H.format == SYM_HALF_MATRIX )
//            {
                Init_CM_Half( system, control, workspace, lists );
//            }
//            else
//            {
//                Init_CM_Full( system, control, data, workspace, lists, out_control );
//            }
        }
        else
        {
//            if ( workspace->H.format == SYM_HALF_MATRIX )
//            {
                Init_CM_Tab_Half( system, control, workspace, lists );
//            }
//            else
//            {
//                Init_CM_Tab_Full( system, control, data, workspace, lists, out_control );
//            }
        }
    }

    if ( bonds_done == FALSE )
    {
        for ( i = 0; i < system->N; ++i )
        {
            workspace->total_bond_order[i] = 0.0;
        }
        for ( i = 0; i < system->N; ++i )
        {
            rvec_MakeZero( workspace->dDeltap_self[i] );
        }

        Init_List_Indices( lists[BONDS], system->bonds );
        if ( control->hbond_cut > 0.0 && workspace->num_H > 0 )
        {
            Init_List_Indices( lists[HBONDS], system->hbonds );
        }

//        if ( lists[FAR_NBRS]->format == HALF_LIST )
//        {
//            Init_Bond_Half( system, control, workspace, lists );
//        }
//        else
//        {
            Init_Bond_Full( system, control, workspace, lists );
//        }
    }

    ret = (workspace->realloc.cm == FALSE
            && workspace->realloc.bonds == FALSE
            && workspace->realloc.hbonds == FALSE ? SUCCESS : FAILURE);

    if ( workspace->realloc.cm == FALSE )
    {
        cm_done = TRUE;
    }
    if ( workspace->realloc.bonds == FALSE && workspace->realloc.hbonds == FALSE )
    {
        bonds_done = TRUE;
    }

    if ( ret == SUCCESS )
    {
        dist_done = FALSE;
        cm_done = FALSE;
        bonds_done = FALSE;
    }

#if defined(TEST_FORCES)
    /* Calculate_dBO requires a sorted bonds list */
//    for ( i = 0; i < bonds->n; ++i )
//    {
//        if ( Num_Entries(i, bonds) > 0 )
//        {
//            qsort( &bonds->bond_list[Start_Index(i, bonds)],
//                    Num_Entries(i, bonds), sizeof(bond_data), compare_bonds );
//        }
//    }
#endif

    return ret;
}


static void Estimate_Storages_CM( reax_system const * const system,
        control_params const * const control, static_storage * const workspace,
        reax_list ** const lists )
{
    int i, pj;
    int start_i, end_i;
    int Htop;
    reax_list *far_nbrs;

    Htop = 0;
    far_nbrs = lists[FAR_NBRS];

    for ( i = 0; i < far_nbrs->n; ++i )
    {
        start_i = Start_Index( i, far_nbrs );
        end_i = End_Index( i, far_nbrs );

        /* update i-j distance - check if j is within cutoff */
        for ( pj = start_i; pj < end_i; ++pj )
        {
            if ( far_nbrs->far_nbr_list[pj].d <= control->nonb_cut )
            {
                ++Htop;
            }
        }

        /* diagonal entry */
        ++Htop;
    }

    switch ( control->charge_method )
    {
        case QEQ_CM:
            break;

        case EE_CM:
            if ( system->num_molec_charge_constraints == 0
                    && system->num_custom_charge_constraints == 0 )
            {
                Htop += system->N_cm;
            }
            else
            {
                for ( i = 0; i < system->num_molec_charge_constraints; ++i )
                {
                    /* extra +1 for zeros on the diagonal */
                    Htop += system->molec_charge_constraint_ranges[2 * i + 1]
                        - system->molec_charge_constraint_ranges[2 * i] + 2;
                }
                for ( i = 0; i < system->num_custom_charge_constraints; ++i )
                {
                    /* +1 for zeros on the diagonal */
                    Htop += system->custom_charge_constraint_count[i] + 1;
                }
            }
            break;

        case ACKS2_CM:
            Htop = 2 * Htop + 3 * system->N + 2;
            break;

        default:
            fprintf( stderr, "[ERROR] Unknown charge method type. Terminating...\n" );
            exit( INVALID_INPUT );
            break;
    }

    workspace->realloc.total_cm_entries = (int) CEIL( Htop * SAFE_ZONE );
}


static void Estimate_Storages_Bonds( reax_system const * const system,
        control_params const * const control, static_storage * const workspace,
        reax_list ** const lists )
{
    int i, j, pj;
    int start_i, end_i;
    int type_i, type_j;
    int ihb, jhb;
    real r_ij;
    real C12, C34, C56;
    real BO, BO_s, BO_pi, BO_pi2;
    reax_list *far_nbrs;
    single_body_parameters *sbp_i, *sbp_j;
    two_body_parameters *twbp;
    far_neighbor_data *nbr_pj;
    reax_atom *atom_i, *atom_j;

    far_nbrs = lists[FAR_NBRS];

    for ( i = 0; i < system->N_max; ++i )
    {
        system->bonds[i] = 0;
    }
    for ( i = 0; i < system->N_max; ++i )
    {
        system->hbonds[i] = 0;
    }

    for ( i = 0; i < far_nbrs->n; ++i )
    {
        atom_i = &system->atoms[i];
        type_i = atom_i->type;
        start_i = Start_Index( i, far_nbrs );
        end_i = End_Index( i, far_nbrs );
        sbp_i = &system->reax_param.sbp[type_i];
        ihb = sbp_i->p_hbond;

        /* update i-j distance - check if j is within cutoff */
        for ( pj = start_i; pj < end_i; ++pj )
        {
            nbr_pj = &far_nbrs->far_nbr_list[pj];

            if ( nbr_pj->d <= control->nonb_cut )
            {
                j = nbr_pj->nbr;
                atom_j = &system->atoms[j];
                type_j = atom_j->type;
                sbp_j = &system->reax_param.sbp[type_j];
                twbp = &system->reax_param.tbp[type_i][type_j];

                /* hydrogen bond lists */
                if ( control->hbond_cut > 0.0
                        && (ihb == H_ATOM || ihb == H_BONDING_ATOM)
                        && nbr_pj->d <= control->hbond_cut )
                {
                    jhb = sbp_j->p_hbond;

                    if ( ihb == H_ATOM && jhb == H_BONDING_ATOM )
                    {
                        ++system->hbonds[i];
                    }
                    else if ( ihb == H_BONDING_ATOM && jhb == H_ATOM )
                    {
                        ++system->hbonds[j];
                    }
                }

                /* uncorrected bond orders */
                if ( nbr_pj->d <= control->bond_cut )
                {
                    r_ij = nbr_pj->d;

                    if ( sbp_i->r_s > 0.0 && sbp_j->r_s > 0.0 )
                    {
                        C12 = twbp->p_bo1 * POW( r_ij / twbp->r_s, twbp->p_bo2 );
                        BO_s = (1.0 + control->bo_cut) * EXP( C12 );
                    }
                    else
                    {
                        C12 = 0.0;
                        BO_s = 0.0;
                    }

                    if ( sbp_i->r_pi > 0.0 && sbp_j->r_pi > 0.0 )
                    {
                        C34 = twbp->p_bo3 * POW( r_ij / twbp->r_p, twbp->p_bo4 );
                        BO_pi = EXP( C34 );
                    }
                    else
                    {
                        C34 = 0.0;
                        BO_pi = 0.0;
                    }

                    if ( sbp_i->r_pi_pi > 0.0 && sbp_j->r_pi_pi > 0.0 )
                    {
                        C56 = twbp->p_bo5 * POW( r_ij / twbp->r_pp, twbp->p_bo6 );
                        BO_pi2 = EXP( C56 );
                    }
                    else
                    {
                        C56 = 0.0;
                        BO_pi2 = 0.0;
                    }

                    /* Initially BO values are the uncorrected ones, page 1 */
                    BO = BO_s + BO_pi + BO_pi2;

                    if ( BO >= control->bo_cut )
                    {
                        ++system->bonds[i];
                        ++system->bonds[j];
                    }
                }
            }
        }
    }

    workspace->realloc.total_thbodies = 0;
    for ( i = 0; i < system->N_max; ++i )
    {
        //TODO: the factor of 2 depends on the bond list format (and also effects thbody list), revisit later
        system->bonds[i] = MAX( (int) CEIL( system->bonds[i] * 2.0 * SAFE_ZONE ), MIN_BONDS );
    }

    for ( i = 0; i < system->N_max; ++i )
    {
        system->hbonds[i] = MAX( (int) CEIL( system->hbonds[i] * SAFE_HBONDS ), MIN_HBONDS );
    }

    /* reductions to get totals */
    workspace->realloc.total_bonds = 0;
    workspace->realloc.total_hbonds = 0;
    workspace->realloc.total_thbodies = 0;

    for ( i = 0; i < system->N_max; ++i )
    {
        workspace->realloc.total_bonds += system->bonds[i];
    }
    for ( i = 0; i < system->N_max; ++i )
    {
        workspace->realloc.total_hbonds += system->hbonds[i];
    }
    for ( i = 0; i < system->N_max; ++i )
    {
        workspace->realloc.total_thbodies += SQR( system->bonds[i] );
    }
}


void Estimate_Storages( reax_system const * const system,
        control_params const * const control, static_storage * const workspace,
        reax_list ** const lists, int realloc_cm, int realloc_bonds )
{
    if ( realloc_cm == TRUE )
    {
        Estimate_Storages_CM( system, control, workspace, lists );
    }

    if ( realloc_bonds == TRUE )
    {
        Estimate_Storages_Bonds( system, control, workspace, lists );
    }
}


int Compute_Forces( reax_system * const system, control_params * const control,
        simulation_data * const data, static_storage * const workspace,
        reax_list ** const lists, output_controls * const out_control, int realloc )
{
    int charge_flag, ret;
    real t_start, t_elapsed;

    if ( control->charge_freq > 0
            && (data->step - data->prev_steps) % control->charge_freq == 0 )
    {
        charge_flag = TRUE;
    }
    else
    {
        charge_flag = FALSE;
    }

    t_start = Get_Time( );
    ret = Init_Forces( system, control, data, workspace, lists, out_control );
    t_elapsed = Get_Timing_Info( t_start );
    data->timing.init_forces += t_elapsed;

    if ( ret != SUCCESS )
    {
        Estimate_Storages( system, control, workspace, lists, workspace->realloc.cm,
               workspace->realloc.bonds || workspace->realloc.hbonds );
    }

    if ( ret == SUCCESS )
    {
        t_start = Get_Time( );
        Compute_Bonded_Forces( system, control, data, workspace, lists, out_control );
        t_elapsed = Get_Timing_Info( t_start );
        data->timing.bonded += t_elapsed;

        if ( charge_flag == TRUE )
        {
            t_start = Get_Time( );
            Compute_Charges( system, control, data, workspace, out_control, realloc );
            t_elapsed = Get_Timing_Info( t_start );
            data->timing.cm += t_elapsed;
        }
            
        if ( control->cm_solver_pre_comp_refactor == -1 )
        {
            if ( data->step <= 4 || is_refactoring_step( control, data ) )
            {
                if ( is_refactoring_step( control, data ) )
                {
                    data->timing.cm_last_pre_comp = data->timing.cm_solver_pre_comp;
                }

                data->timing.cm_optimum = data->timing.cm_solver_pre_app
                    + data->timing.cm_solver_spmv
                    + data->timing.cm_solver_vector_ops
                    + data->timing.cm_solver_orthog
                    + data->timing.cm_solver_tri_solve;
                data->timing.cm_total_loss = 0.0;
            }
            else
            {
                data->timing.cm_total_loss += data->timing.cm_solver_pre_app
                    + data->timing.cm_solver_spmv
                    + data->timing.cm_solver_vector_ops
                    + data->timing.cm_solver_orthog
                    + data->timing.cm_solver_tri_solve
                    - data->timing.cm_optimum;
            }
        }

        t_start = Get_Time( );
        Compute_NonBonded_Forces( system, control, data, workspace,
                lists, out_control, realloc );
        t_elapsed = Get_Timing_Info( t_start );
        data->timing.nonb += t_elapsed;

        Compute_Total_Force( system, control, data, workspace, lists );

#if defined(TEST_FORCES)
        Print_Total_Force( system, control, data, workspace, lists, out_control );
        Compare_Total_Forces( system, control, data, workspace, lists, out_control );
#endif
    }

    return ret;
}
