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

#include "restart.h"

#include "allocate.h"
#include "box.h"
#include "tool_box.h"
#include "vector.h"


void Write_Binary_Restart( reax_system *system, control_params *control,
        simulation_data *data, static_storage *workspace )
{
    int  i;
    char fname[MAX_STR];
    FILE *fres;
    reax_atom *p_atom;
    restart_header res_header;
    restart_atom res_data;

    snprintf( fname, MAX_STR, "%.*s.res%d", MAX_STR - 12, control->sim_name, data->step );
    fres = sfopen( fname, "wb", __FILE__, __LINE__ );

    res_header.step = data->step;
    res_header.N = system->N;
    res_header.T = data->therm.T;
    res_header.xi = data->therm.xi;
    res_header.v_xi = data->therm.v_xi;
    res_header.v_xi_old = data->therm.v_xi_old;
    res_header.G_xi = data->therm.G_xi;
    rtensor_Copy( res_header.box, system->box.box );
    fwrite( &res_header, sizeof(restart_header), 1, fres );

    for ( i = 0; i < system->N; ++i )
    {
        p_atom = &system->atoms[i];

        res_data.orig_id = workspace->orig_id[i];
        res_data.type = p_atom->type;
        strncpy( res_data.name, p_atom->name, sizeof(res_data.name) - 1 );
        res_data.name[sizeof(res_data.name) - 1] = '\0';
        rvec_Copy( res_data.x, p_atom->x );
        rvec_Copy( res_data.v, p_atom->v );
        fwrite( &res_data, sizeof(restart_atom), 1, fres );
    }

    sfclose( fres, __FILE__, __LINE__ );

#if defined(DEBUG_FOCUS)
    fprintf( stderr, "write restart - " );
#endif
}


void Read_Binary_Restart( const char * const fname, reax_system *system,
        control_params *control, simulation_data *data,
        static_storage *workspace )
{
    int i;
    FILE *fres;
    reax_atom *p_atom;
    restart_header res_header;
    restart_atom res_data;

    fres = sfopen( fname, "rb", __FILE__, __LINE__ );

    /* parse header of restart file */
    if ( fread( &res_header, sizeof(restart_header), 1, fres ) != 1 )
    {
        fprintf( stderr, "[ERROR] reading restart file failed\n" \
                "  [INFO] reading header info\n" );
        exit( INVALID_INPUT );
    }

    data->prev_steps = res_header.step;
    data->therm.T = res_header.T;
    data->therm.xi = res_header.xi;
    data->therm.v_xi = res_header.v_xi;
    data->therm.v_xi_old = res_header.v_xi_old;
    data->therm.G_xi = res_header.G_xi;

    rvec_MakeZero( system->box.min );
    Update_Box( res_header.box, &system->box );

#if defined(DEBUG_FOCUS)
    fprintf( stderr, "restart step: %d\n", data->prev_steps );
    fprintf( stderr, "restart thermostat: %10.6f %10.6f %10.6f %10.6f %10.6f\n",
             data->therm.T, data->therm.xi,
             data->therm.v_xi, data->therm.v_xi_old, data->therm.G_xi );
    fprintf( stderr, "restart box:\n" );
    fprintf( stderr, "%9.5f %9.5f %9.5f\n%9.5f %9.5f %9.5f\n%9.5f %9.5f %9.5f\n",
             system->box.box[0][0], system->box.box[0][1], system->box.box[0][2],
             system->box.box[1][0], system->box.box[1][1], system->box.box[1][2],
             system->box.box[2][0], system->box.box[2][1], system->box.box[2][2] );
#endif

    if ( system->prealloc_allocated == FALSE || res_header.N > system->N_max )
    {
        PreAllocate_Space( system, control, workspace,
                (int) CEIL( SAFE_ZONE * res_header.N ) );
    }
    system->N = res_header.N;

    for ( i = 0; i < MAX_ATOM_ID; ++i )
    {
        workspace->map_serials[i] = -1;
    }

    for ( i = 0; i < system->N; ++i )
    {
        if ( fread( &res_data, sizeof(restart_atom), 1, fres ) != 1 )
        {
            fprintf( stderr, "[ERROR] reading restart file failed\n" \
                    "  [INFO] reading atom info (%d entry)\n", i );
            exit( INVALID_INPUT );
        }

        workspace->orig_id[i] = res_data.orig_id;
        workspace->map_serials[res_data.orig_id] = i;

        p_atom = &system->atoms[i];
        p_atom->type = res_data.type;
        strncpy( p_atom->name, res_data.name, sizeof(p_atom->name) - 1 );
        p_atom->name[sizeof(p_atom->name) - 1] = '\0';
        rvec_Copy( p_atom->x, res_data.x );
        rvec_Copy( p_atom->v, res_data.v );
    }

#if defined(DEBUG_FOCUS)
    fprintf( stderr, "system->N: %d, i: %d\n", system->N, i );
#endif

    sfclose( fres, __FILE__, __LINE__ );

    data->step = data->prev_steps;
    /* target num. of MD sim. steps (nsteps)
     * is updated based on the number of steps in the previous run */
    control->nsteps += data->prev_steps;
}


void Write_ASCII_Restart( reax_system *system, control_params *control,
        simulation_data *data, static_storage *workspace )
{
    int  i;
    char fname[MAX_STR];
    FILE *fres;
    reax_atom *p_atom;

    snprintf( fname, MAX_STR + 8, "%s.res%d", control->sim_name, data->step );
    fres = sfopen( fname, "w", __FILE__, __LINE__ );

    fprintf( fres, RESTART_HEADER,
             data->step, system->N, data->therm.T, data->therm.xi,
             data->therm.v_xi, data->therm.v_xi_old, data->therm.G_xi,
             system->box.box[0][0], system->box.box[0][1], system->box.box[0][2],
             system->box.box[1][0], system->box.box[1][1], system->box.box[1][2],
             system->box.box[2][0], system->box.box[2][1], system->box.box[2][2]);
    fflush( fres );

    for ( i = 0; i < system->N; ++i )
    {
        p_atom = &system->atoms[i];
        fprintf( fres, RESTART_LINE,
                 workspace->orig_id[i], p_atom->type, p_atom->name,
                 p_atom->x[0], p_atom->x[1], p_atom->x[2],
                 p_atom->v[0], p_atom->v[1], p_atom->v[2] );
    }

    sfclose( fres, __FILE__, __LINE__ );

#if defined(DEBUG_FOCUS)
    fprintf( stderr, "write restart - " );
#endif
}


void Read_ASCII_Restart( const char * const fname, reax_system *system,
        control_params *control, simulation_data *data,
        static_storage *workspace )
{
    int i, n;
    FILE *fres;
    reax_atom *p_atom;

    fres = sfopen( fname, "r", __FILE__, __LINE__ );

    /* parse header of restart file */
    if ( fscanf( fres, READ_RESTART_HEADER,
                &data->prev_steps, &n, &data->therm.T, &data->therm.xi,
                &data->therm.v_xi, &data->therm.v_xi_old, &data->therm.G_xi,
                &system->box.box[0][0], &system->box.box[0][1], &system->box.box[0][2],
                &system->box.box[1][0], &system->box.box[1][1], &system->box.box[1][2],
                &system->box.box[2][0], &system->box.box[2][1], &system->box.box[2][2] ) != 16 )
    {
        fprintf( stderr, "[ERROR] reading restart file failed\n" \
                 "  [INFO] reading header info\n" );
        exit( INVALID_INPUT );
    }

    rvec_MakeZero( system->box.min );
    Make_Consistent( &system->box );

#if defined(DEBUG_FOCUS)
    fprintf( stderr, "restart step: %d\n", data->prev_steps );
    fprintf( stderr, "restart thermostat: %10.6f %10.6f %10.6f %10.6f %10.6f\n",
             data->therm.T, data->therm.xi,
             data->therm.v_xi, data->therm.v_xi_old, data->therm.G_xi );
    fprintf( stderr, "restart box:\n" );
    fprintf( stderr, "%9.5f %9.5f %9.5f\n%9.5f %9.5f %9.5f\n%9.5f %9.5f %9.5f\n",
             system->box.box[0][0], system->box.box[0][1], system->box.box[0][2],
             system->box.box[1][0], system->box.box[1][1], system->box.box[1][2],
             system->box.box[2][0], system->box.box[2][1], system->box.box[2][2] );
#endif

    if ( system->prealloc_allocated == FALSE || n > system->N_max )
    {
        PreAllocate_Space( system, control, workspace, (int) CEIL( SAFE_ZONE * n ) );
    }
    system->N = n;

    for ( i = 0; i < MAX_ATOM_ID; ++i )
    {
        workspace->map_serials[i] = -1;
    }

    for ( i = 0; i < system->N; ++i )
    {
        p_atom = &system->atoms[i];
        if ( fscanf( fres, READ_RESTART_LINE,
                    &workspace->orig_id[i], &p_atom->type, p_atom->name,
                    &p_atom->x[0], &p_atom->x[1], &p_atom->x[2],
                    &p_atom->v[0], &p_atom->v[1], &p_atom->v[2] ) != 9 )
        {
            fprintf( stderr, "[ERROR] reading restart file failed\n" \
                    "  [INFO] reading atom info (%d entry)\n", i );
            exit( INVALID_INPUT );
        }
        workspace->map_serials[workspace->orig_id[i]] = i;
    }

    sfclose( fres, __FILE__, __LINE__ );

    data->step = data->prev_steps;
    /* target num. of MD sim. steps (nsteps)
     * is updated based on the number of steps in the previous run */
    control->nsteps += data->prev_steps;
}


void Write_Restart( reax_system *system, control_params *control,
        simulation_data *data, static_storage *workspace,
        output_controls *out_control )
{
    switch ( out_control->restart_format )
    {
    case WRITE_ASCII:
        Write_ASCII_Restart( system, control, data, workspace );
        break;

    case WRITE_BINARY:
        Write_Binary_Restart( system, control, data, workspace );
        break;

    default:
        fprintf( stderr, "[ERROR] invalid restart format\n" );
        exit( INVALID_INPUT );
    }
}
