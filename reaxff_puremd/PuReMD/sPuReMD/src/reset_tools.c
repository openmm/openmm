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

#include "reset_tools.h"

#include "list.h"
#include "vector.h"


void Reset_Pressures( control_params *control, simulation_data *data )
{
#if defined(_OPENMP)
    int i;
#endif

    rtensor_MakeZero( data->flex_bar.P );
    data->iso_bar.P = 0.0;
    rtensor_MakeZero( data->press );
#if defined(_OPENMP)
    if ( control->ensemble == sNPT || control->ensemble == iNPT
            || control->ensemble == aNPT || control->compute_pressure == TRUE )
    {
        for ( i = 0; i < control->num_threads; ++i )
        {
            rtensor_MakeZero( data->press_local[i] );
        }
    }
#endif
}


#if defined(TEST_FORCES)
static void Reset_Test_Forces( reax_system *system, static_storage *workspace )
{
    int i;

    for ( i = 0; i < system->N; i++ )
    {
        rvec_MakeZero( workspace->dDelta[i] );
    }
    for ( i = 0; i < system->N; i++ )
    {
        rvec_MakeZero( workspace->f_ele[i] );
    }
    for ( i = 0; i < system->N; i++ )
    {
        rvec_MakeZero( workspace->f_vdw[i] );
    }
    for ( i = 0; i < system->N; i++ )
    {
        rvec_MakeZero( workspace->f_be[i] );
    }
    for ( i = 0; i < system->N; i++ )
    {
        rvec_MakeZero( workspace->f_lp[i] );
    }
    for ( i = 0; i < system->N; i++ )
    {
        rvec_MakeZero( workspace->f_ov[i] );
    }
    for ( i = 0; i < system->N; i++ )
    {
        rvec_MakeZero( workspace->f_un[i] );
    }
    for ( i = 0; i < system->N; i++ )
    {
        rvec_MakeZero( workspace->f_ang[i] );
    }
    for ( i = 0; i < system->N; i++ )
    {
        rvec_MakeZero( workspace->f_coa[i] );
    }
    for ( i = 0; i < system->N; i++ )
    {
        rvec_MakeZero( workspace->f_pen[i] );
    }
    for ( i = 0; i < system->N; i++ )
    {
        rvec_MakeZero( workspace->f_hb[i] );
    }
    for ( i = 0; i < system->N; i++ )
    {
        rvec_MakeZero( workspace->f_tor[i] );
    }
    for ( i = 0; i < system->N; i++ )
    {
        rvec_MakeZero( workspace->f_con[i] );
    }
}
#endif


void Reset_Atomic_Forces( reax_system* system )
{
    int i;

    for ( i = 0; i < system->N; ++i )
    {
        rvec_MakeZero( system->atoms[i].f );
    }
}


void Reset_Energies( simulation_data* data )
{
    data->E_Tot = 0.0;
    data->E_Kin = 0.0;
    data->E_Pot = 0.0;
    data->E_BE = 0.0;
    data->E_Ov = 0.0;
    data->E_Un = 0.0;
    data->E_Lp = 0.0;
    data->E_Ang = 0.0;
    data->E_Pen = 0.0;
    data->E_Coa = 0.0;
    data->E_HB = 0.0;
    data->E_Tor = 0.0;
    data->E_Con = 0.0;
    data->E_vdW = 0.0;
    data->E_Ele = 0.0;
    data->E_Pol = 0.0;
}


void Reset_Workspace( reax_system *system, static_storage *workspace )
{
    int i;
#if defined(_OPENMP)
    int tid;
#endif

    for ( i = 0; i < system->N; i++ )
    {
        workspace->total_bond_order[i] = 0.0;
    }
    for ( i = 0; i < system->N; i++ )
    {
        rvec_MakeZero( workspace->dDeltap_self[i] );
    }
    for ( i = 0; i < system->N; i++ )
    {
        workspace->CdDelta[i] = 0.0;
    }

#if defined(_OPENMP)
    #pragma omp parallel private(i, tid)
    {
        tid = omp_get_thread_num( );

        for ( i = 0; i < system->N; ++i )
        {
            rvec_MakeZero( workspace->f_local[tid * system->N + i] );
        }
    }
#endif

#if defined(TEST_FORCES)
    Reset_Test_Forces( system, workspace );
#endif
}


void Reset( reax_system *system, control_params *control,
        simulation_data *data, static_storage *workspace, reax_list **lists  )
{
    Reset_Atomic_Forces( system );

    Reset_Energies( data );

    if ( control->ensemble == sNPT || control->ensemble == iNPT
            || control->ensemble == aNPT || control->compute_pressure == TRUE )
    {
        Reset_Pressures( control, data );
    }

    Reset_Workspace( system, workspace );
}


void Reset_Grid( grid *g )
{
    int i, j, k;

    for ( i = 0; i < g->ncell[0]; i++ )
    {
        for ( j = 0; j < g->ncell[1]; j++ )
        {
            for ( k = 0; k < g->ncell[2]; k++ )
            {
                g->top[i][j][k] = 0;
            }
        }
    }
}


void Reset_Marks( grid *g, ivec *grid_stack, int grid_top )
{
    int i;

    for ( i = 0; i < grid_top; ++i )
    {
        g->mark[grid_stack[i][0]][grid_stack[i][1]][grid_stack[i][2]] = 0;
    }
}
