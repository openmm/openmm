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

#include "control.h"

#include "box.h"
#include "traj.h"
#include "tool_box.h"

#include <ctype.h>
#if defined(HAVE_ZLIB) && defined(HAVE_ZLIB_H)
  #include <zlib.h>
#endif


int Set_Control_Parameter( const char * const keyword,
        const char ** const values, control_params * const control,
        output_controls * const out_control )
{
    int i, ret;
    real val;

    ret = SUCCESS;

    if ( strncmp(keyword, "simulation_name", MAX_LINE) == 0 )
    {
        strncpy( control->sim_name, values[0], sizeof(control->sim_name) - 1 );
        control->sim_name[sizeof(control->sim_name) - 1] = '\0';
    }
    //else if( strncmp(keyword, "restart", MAX_LINE) == 0 ) {
    //  control->restart = sstrtol( values[0], __FILE__, __LINE__ );
    //}
    else if ( strncmp(keyword, "restart_format", MAX_LINE) == 0 )
    {
        out_control->restart_format = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "restart_freq", MAX_LINE) == 0 )
    {
        out_control->restart_freq = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "random_vel", MAX_LINE) == 0 )
    {
        control->random_vel = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "reposition_atoms", MAX_LINE) == 0 )
    {
        control->reposition_atoms = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "ensemble_type", MAX_LINE) == 0 )
    {
        control->ensemble = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "nsteps", MAX_LINE) == 0 )
    {
        control->nsteps = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "dt", MAX_LINE) == 0 )
    {
        val = sstrtod( values[0], __FILE__, __LINE__ );
        /* convert from fs to ps */
        control->dt = val * 1.0e-3;
    }
    else if ( strncmp(keyword, "num_threads", MAX_LINE) == 0 )
    {
        control->num_threads = sstrtol( values[0], __FILE__, __LINE__ );
        control->num_threads_set = TRUE;
    }
    else if ( strncmp(keyword, "gpus_per_node", MAX_LINE) == 0 )
    {
        // skip since not applicable to shared memory code
        ;
    }
    else if ( strncmp(keyword, "gpu_streams", MAX_LINE) == 0 )
    {
        // skip since not applicable to shared memory code
        ;
    }
    else if ( strncmp(keyword, "proc_by_dim", MAX_LINE) == 0 )
    {
        // skip since not applicable to shared memory code
        ;
    }
    else if ( strncmp(keyword, "periodic_boundaries", MAX_LINE) == 0 )
    {
        control->periodic_boundaries = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "geo_format", MAX_LINE) == 0 )
    {
        control->geo_format = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "restrict_bonds", MAX_LINE) == 0 )
    {
        control->restrict_bonds = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "tabulate_long_range", MAX_LINE) == 0 )
    {
        control->tabulate = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "reneighbor", MAX_LINE) == 0 )
    {
        control->reneighbor = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "vlist_buffer", MAX_LINE) == 0 )
    {
        val = sstrtod( values[0], __FILE__, __LINE__ );
        control->vlist_cut = control->nonb_cut + val;
    }
    else if ( strncmp(keyword, "nbrhood_cutoff", MAX_LINE) == 0 )
    {
        val = sstrtod( values[0], __FILE__, __LINE__ );
        control->bond_cut = val;
    }
    else if ( strncmp(keyword, "thb_cutoff", MAX_LINE) == 0 )
    {
        val = sstrtod( values[0], __FILE__, __LINE__ );
        control->thb_cut = val;
    }
    else if ( strncmp(keyword, "hbond_cutoff", MAX_LINE) == 0 )
    {
        val = sstrtod( values[0], __FILE__, __LINE__ );
        control->hbond_cut = val;
    }
    else if ( strncmp(keyword, "charge_method", MAX_LINE) == 0 )
    {
        control->charge_method = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "charge_freq", MAX_LINE) == 0 )
    {
        control->charge_freq = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "cm_q_net", MAX_LINE) == 0 )
    {
        val = sstrtod( values[0], __FILE__, __LINE__ );
        control->cm_q_net = val;
    }
    else if ( strncmp(keyword, "cm_solver_type", MAX_LINE) == 0 )
    {
        control->cm_solver_type = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "cm_solver_max_iters", MAX_LINE) == 0 )
    {
        control->cm_solver_max_iters = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "cm_solver_restart", MAX_LINE) == 0 )
    {
        control->cm_solver_restart = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "cm_solver_q_err", MAX_LINE) == 0 )
    {
        val = sstrtod( values[0], __FILE__, __LINE__ );
        control->cm_solver_q_err = val;
    }
    else if ( strncmp(keyword, "cm_domain_sparsity", MAX_LINE) == 0 )
    {
        val = sstrtod( values[0], __FILE__, __LINE__ );
        control->cm_domain_sparsity = val;
        if ( val < 1.0 )
        {
            control->cm_domain_sparsify_enabled = TRUE;
        }
    }
    else if ( strncmp(keyword, "cm_init_guess_type", MAX_LINE) == 0 )
    {
        control->cm_init_guess_type = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "cm_init_guess_extrap1", MAX_LINE) == 0 )
    {
        control->cm_init_guess_extrap1 = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "cm_init_guess_extrap2", MAX_LINE) == 0 )
    {
        control->cm_init_guess_extrap2 = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "cm_init_guess_gd_model", MAX_LINE) == 0 )
    {
        strncpy( control->cm_init_guess_gd_model, values[0], sizeof(control->cm_init_guess_gd_model) - 1 );
        control->cm_init_guess_gd_model[sizeof(control->cm_init_guess_gd_model) - 1] = '\0';
    }
    else if ( strncmp(keyword, "cm_init_guess_win_size", MAX_LINE) == 0 )
    {
        control->cm_init_guess_win_size = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "cm_solver_pre_comp_type", MAX_LINE) == 0 )
    {
        control->cm_solver_pre_comp_type = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "cm_solver_pre_comp_refactor", MAX_LINE) == 0 )
    {
        control->cm_solver_pre_comp_refactor = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "cm_solver_pre_comp_droptol", MAX_LINE) == 0 )
    {
        val = sstrtod( values[0], __FILE__, __LINE__ );
        control->cm_solver_pre_comp_droptol = val;
    }
    else if ( strncmp(keyword, "cm_solver_pre_comp_sweeps", MAX_LINE) == 0 )
    {
        control->cm_solver_pre_comp_sweeps = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "cm_solver_pre_comp_sai_thres", MAX_LINE) == 0 )
    {
        val = sstrtod( values[0], __FILE__, __LINE__ );
        control->cm_solver_pre_comp_sai_thres = val;
    }
    else if ( strncmp(keyword, "cm_solver_pre_app_type", MAX_LINE) == 0 )
    {
        control->cm_solver_pre_app_type = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "cm_solver_pre_app_jacobi_iters", MAX_LINE) == 0 )
    {
        control->cm_solver_pre_app_jacobi_iters = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "include_polarization_energy", MAX_LINE) == 0 )
    {
        control->polarization_energy_enabled = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "temp_init", MAX_LINE) == 0 )
    {
        val = sstrtod( values[0], __FILE__, __LINE__ );
        control->T_init = val;

        if ( control->T_init < 0.001 )
        {
            control->T_init = 0.001;
        }
    }
    else if ( strncmp(keyword, "temp_final", MAX_LINE) == 0 )
    {
        val = sstrtod( values[0], __FILE__, __LINE__ );
        control->T_final = val;

        if ( control->T_final < 0.1 )
        {
            control->T_final = 0.1;
        }
    }
    else if ( strncmp(keyword, "t_mass", MAX_LINE) == 0 )
    {
        val = sstrtod( values[0], __FILE__, __LINE__ );
        /* convert from fs to s */
        control->Tau_T = val * 1.0e-15;
    }
    else if ( strncmp(keyword, "t_mode", MAX_LINE) == 0 )
    {
        control->T_mode = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "t_rate", MAX_LINE) == 0 )
    {
        val = sstrtod( values[0], __FILE__, __LINE__ );
        control->T_rate = val;
    }
    else if ( strncmp(keyword, "t_freq", MAX_LINE) == 0 )
    {
        val = sstrtod( values[0], __FILE__, __LINE__ );
        control->T_freq = val;
    }
    else if ( strncmp(keyword, "pressure", MAX_LINE) == 0 )
    {
        if ( control->ensemble == iNPT )
        {
            val = sstrtod( values[0], __FILE__, __LINE__ );
            control->P[0] = val;
            control->P[1] = val;
            control->P[2] = val;
        }
        else if ( control->ensemble == sNPT )
        {
            val = sstrtod( values[0], __FILE__, __LINE__ );
            control->P[0] = val;

            val = sstrtod( values[1], __FILE__, __LINE__ );
            control->P[1] = val;

            val = sstrtod( values[2], __FILE__, __LINE__ );
            control->P[2] = val;
        }
    }
    else if ( strncmp(keyword, "p_mass", MAX_LINE) == 0 )
    {
        if ( control->ensemble == iNPT )
        {
            val = sstrtod( values[0], __FILE__, __LINE__ );
            control->Tau_P[0] = val * 1.0e-3;   // convert p_mass from fs to ps
        }
        else if ( control->ensemble == sNPT )
        {
            val = sstrtod( values[0], __FILE__, __LINE__ );
            control->Tau_P[0] = val * 1.0e-3;   // convert p_mass from fs to ps

            val = sstrtod( values[1], __FILE__, __LINE__ );
            control->Tau_P[1] = val * 1.0e-3;   // convert p_mass from fs to ps

            val = sstrtod( values[2], __FILE__, __LINE__ );
            control->Tau_P[2] = val * 1.0e-3;   // convert p_mass from fs to ps
        }
    }
    else if ( strncmp(keyword, "pt_mass", MAX_LINE) == 0 )
    {
        val = sstrtod( values[0], __FILE__, __LINE__ );
        control->Tau_PT = val * 1.0e-3;  // convert pt_mass from fs to ps
    }
    else if ( strncmp(keyword, "compress", MAX_LINE) == 0 )
    {
        val = sstrtod( values[0], __FILE__, __LINE__ );
        control->compressibility = val;
    }
    else if ( strncmp(keyword, "press_mode", MAX_LINE) == 0 )
    {
        control->press_mode = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "compute_pressure", MAX_LINE) == 0 )
    {
        control->compute_pressure = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "remove_CoM_vel", MAX_LINE) == 0 )
    {
        control->remove_CoM_vel = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "energy_update_freq", MAX_LINE) == 0 )
    {
        out_control->log_update_freq = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "write_freq", MAX_LINE) == 0 )
    {
        out_control->write_steps = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "traj_compress", MAX_LINE) == 0 )
    {
        out_control->traj_compress = sstrtol( values[0], __FILE__, __LINE__ );

        if ( out_control->traj_compress == 1 )
        {
#if defined(HAVE_ZLIB)
            out_control->write = (int (*)(FILE *, const char *, ...)) &gzprintf;
#else
            fprintf( stderr, "[ERROR] zlib support disabled (trajectory file compression). "
                    "Re-compile to enable. Terminating...\n" );
            exit( INVALID_INPUT );
#endif
        }
        else
        {
            out_control->write = &fprintf;
        }
    }
    else if ( strncmp(keyword, "traj_format", MAX_LINE) == 0 )
    {
        out_control->traj_format = sstrtol( values[0], __FILE__, __LINE__ );

        if ( out_control->traj_format == 0 )
        {
            out_control->write_header = &Write_Custom_Header;
            out_control->append_traj_frame = &Append_Custom_Frame;
        }
        else if ( out_control->traj_format == 1 )
        {
            out_control->write_header = &Write_xyz_Header;
            out_control->append_traj_frame = &Append_xyz_Frame;
        }
    }
    else if ( strncmp(keyword, "traj_title", MAX_LINE) == 0 )
    {
        strncpy( out_control->traj_title, values[0], sizeof(out_control->traj_title) - 1 );
        out_control->traj_title[sizeof(out_control->traj_title) - 1] = '\0';
    }
    else if ( strncmp(keyword, "atom_info", MAX_LINE) == 0 )
    {
        out_control->atom_format +=
            sstrtol( values[0], __FILE__, __LINE__ ) * 4;
    }
    else if ( strncmp(keyword, "atom_velocities", MAX_LINE) == 0 )
    {
        out_control->atom_format +=
            sstrtol( values[0], __FILE__, __LINE__ ) * 2;
    }
    else if ( strncmp(keyword, "atom_forces", MAX_LINE) == 0 )
    {
        out_control->atom_format +=
            sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "bond_info", MAX_LINE) == 0 )
    {
        out_control->bond_info = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "angle_info", MAX_LINE) == 0 )
    {
        out_control->angle_info = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "test_forces", MAX_LINE) == 0 )
    {
        ;
    }
    else if ( strncmp(keyword, "molec_anal", MAX_LINE) == 0 )
    {
        control->molec_anal = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "freq_molec_anal", MAX_LINE) == 0 )
    {
        control->freq_molec_anal = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "bond_graph_cutoff", MAX_LINE) == 0 )
    {
        val = sstrtod( values[0], __FILE__, __LINE__ );
        control->bg_cut = val;
    }
    else if ( strncmp(keyword, "ignore", MAX_LINE) == 0 )
    {
        control->num_ignored = sstrtol( values[0], __FILE__, __LINE__ );
        for ( i = 0; i < control->num_ignored; ++i )
            control->ignore[ sstrtol( values[i + 1], __FILE__, __LINE__ ) ] = 1;
    }
    else if ( strncmp(keyword, "dipole_anal", MAX_LINE) == 0 )
    {
        control->dipole_anal = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "freq_dipole_anal", MAX_LINE) == 0 )
    {
        control->freq_dipole_anal = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "diffusion_coef", MAX_LINE) == 0 )
    {
        control->diffusion_coef = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "freq_diffusion_coef", MAX_LINE) == 0 )
    {
        control->freq_diffusion_coef = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else if ( strncmp(keyword, "restrict_type", MAX_LINE) == 0 )
    {
        control->restrict_type = sstrtol( values[0], __FILE__, __LINE__ );
    }
    else
    {
        ret = FAILURE;
    }

    return ret;
}


/* Assign default values to control parameters
 *
 * NOTE: force field file parameters must be parsed before this
 * */
void Set_Control_Defaults( reax_system * const system,
        control_params * const  control, output_controls * const out_control )
{
    int i;

    strncpy( control->sim_name, "default.sim", sizeof(control->sim_name) - 1 );
    control->sim_name[sizeof(control->sim_name) - 1] = '\0';

    strncpy( control->restart_from, "default.res", sizeof(control->restart_from) - 1 );
    control->restart_from[sizeof(control->restart_from) - 1] = '\0';
    control->restart = FALSE;
    out_control->restart_format = WRITE_BINARY;
    out_control->restart_freq = 0;
    out_control->restart_freq = 0;

    control->random_vel = FALSE;
    control->reposition_atoms = 0;
    control->ensemble = NVE;
    control->nsteps = 0;
    control->periodic_boundaries = TRUE;
    control->restrict_bonds = FALSE;
    control->tabulate = 0;
    control->dt = 0.25;
    control->num_threads_set = FALSE;
    control->reneighbor = 1;

    /* defaults values for other cutoffs */
    control->vlist_cut = system->reax_param.gp.l[12] + 2.5;
    control->bond_cut = 5.0;
    control->nonb_low = system->reax_param.gp.l[11];
    control->nonb_cut = system->reax_param.gp.l[12];
    control->bo_cut = 0.01 * system->reax_param.gp.l[29];
    control->thb_cut = 0.001;
    control->hbond_cut = 0.0;

    control->T_init = 0.0;
    control->T_final = 300.0;
    control->Tau_T = 1.0;
    control->T_mode = 0;
    control->T_rate = 1.0;
    control->T_freq = 1.0;
    control->P[0] = 0.000101325;
    control->P[1] = 0.000101325;
    control->P[2] = 0.000101325;
    control->Tau_P[0] = 500.0;
    control->Tau_P[1] = 500.0;
    control->Tau_P[2] = 500.0;
    control->Tau_PT = 500.0;
    control->compressibility = 1.0;
    control->press_mode = 0;
    control->compute_pressure = FALSE;

    control->remove_CoM_vel = 25;
    control->geo_format = PDB;
    control->dipole_anal = 0;
    control->freq_dipole_anal = 0;
    control->diffusion_coef = 0;
    control->freq_diffusion_coef = 0;
    control->restrict_type = 0;

    control->charge_method = QEQ_CM;
    control->charge_freq = 1;
    control->cm_q_net = 0.0;
    control->cm_solver_type = GMRES_S;
    control->cm_solver_max_iters = 100;
    control->cm_solver_restart = 50;
    control->cm_solver_q_err = 0.000001;
    control->cm_domain_sparsify_enabled = FALSE;
    control->cm_domain_sparsity = 1.0;
    control->cm_init_guess_type = SPLINE;
    control->cm_init_guess_extrap1 = 3;
    control->cm_init_guess_extrap2 = 2;
    /* assign default values */
    strncpy( control->cm_init_guess_gd_model, "frozen_model.pb",
            sizeof(control->cm_init_guess_gd_model) - 1 );
    control->cm_init_guess_gd_model[sizeof(control->cm_init_guess_gd_model) - 1] = '\0';
    control->cm_init_guess_win_size = 5;
    control->cm_solver_pre_comp_type = JACOBI_PC;
    control->cm_solver_pre_comp_sweeps = 3;
    control->cm_solver_pre_comp_sai_thres = 0.1;
    control->cm_solver_pre_comp_refactor = 1;
    control->cm_solver_pre_comp_droptol = 0.01;
    control->cm_solver_pre_app_type = TRI_SOLVE_PA;
    control->cm_solver_pre_app_jacobi_iters = 50;
    control->polarization_energy_enabled = TRUE;

    control->molec_anal = NO_ANALYSIS;
    control->freq_molec_anal = 0;
    control->bg_cut = 0.3;
    control->num_ignored = 0;
    for ( i = 0; i < MAX_ATOM_TYPES; i++ )
    {
        control->ignore[i] = 0;
    }

    out_control->log_update_freq = 1;
    out_control->write_steps = 0;
    out_control->traj_compress = 0;
    out_control->write = &fprintf;
    out_control->traj_format = 0;
    out_control->write_header = &Write_Custom_Header;
    out_control->append_traj_frame = &Append_Custom_Frame;
    strncpy( out_control->traj_title, "default_title", sizeof(out_control->traj_title) - 1 );
    out_control->traj_title[sizeof(out_control->traj_title) - 1] = '\0';
    out_control->atom_format = 0;
    out_control->bond_info = 0;
    out_control->angle_info = 0;
}


void Read_Control_File( const char * const control_file, reax_system * const system,
        control_params * const  control, output_controls * const out_control )
{
    char *s, **tmp;
    int c, i, ret;
    FILE *fp;

    fp = sfopen( control_file, "r", __FILE__, __LINE__ );

    assert( fp != NULL );

    if ( fp != NULL )
    {
        s = smalloc( sizeof(char) * MAX_LINE, __FILE__, __LINE__ );
        tmp = smalloc( sizeof(char*) * MAX_TOKENS, __FILE__, __LINE__ );
        for ( i = 0; i < MAX_TOKENS; i++ )
        {
            tmp[i] = smalloc( sizeof(char) * MAX_LINE, __FILE__, __LINE__ );
        }

        /* read control parameters file */
        while ( fgets( s, MAX_LINE, fp ) )
        {
            c = Tokenize( s, &tmp, MAX_LINE );

            if ( c > 0 )
            {
                ret = Set_Control_Parameter( tmp[0],
                        (const char ** const) &tmp[1], control, out_control );

                if ( ret == FAILURE )
                {
                    fprintf( stderr, "WARNING: unknown parameter %s\n", tmp[0] );
                    exit( UNKNOWN_OPTION );
                }
            }
        }

        if ( ferror(fp) )
        {
            fprintf( stderr, "Error reading control file. Terminating.\n" );
            exit( INVALID_INPUT );
        }

        for ( i = 0; i < MAX_TOKENS; i++ )
        {
            sfree( tmp[i], __FILE__, __LINE__ );
        }
        sfree( tmp, __FILE__, __LINE__ );
        sfree( s, __FILE__, __LINE__ );
    }

    sfclose( fp, __FILE__, __LINE__ );
}


void Set_Control_Derived_Values( reax_system * const system,
        control_params * const  control )
{
    int i;

    /* determine target T */
    if ( control->T_mode == 0 )
    {
        control->T = control->T_final;
    }
    else
    {
        control->T = control->T_init;
    }

    /* periodic boundary conditions */
    if ( control->periodic_boundaries == TRUE )
    {
        control->compute_atom_distance = &Compute_Atom_Distance_Periodic;
        control->update_atom_position = &Update_Atom_Position_Periodic;
    }
    /* non-periodic boundary conditions */
    else
    {
        control->compute_atom_distance = &Compute_Atom_Distance_Non_Periodic;
        control->update_atom_position = &Update_Atom_Position_Non_Periodic;
    }

    /* derived cutoffs */
    control->nonb_sp_cut = control->nonb_cut * control->cm_domain_sparsity;

    system->g.cell_size = control->vlist_cut / 2.0;
    for ( i = 0; i < 3; ++i )
    {
        system->g.spread[i] = 2;
    }
}
