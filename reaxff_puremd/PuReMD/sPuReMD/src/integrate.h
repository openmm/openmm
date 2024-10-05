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

#ifndef __INTEGRATE_H_
#define __INTEGRATE_H_


#include "reax_types.h"


int Velocity_Verlet_NVE( reax_system *, control_params *, simulation_data *,
        static_storage *, reax_list **, output_controls * );

int Velocity_Verlet_Berendsen_NVT( reax_system *, control_params *,
        simulation_data *, static_storage *, reax_list **, output_controls * );

int Velocity_Verlet_Nose_Hoover_NVT_Klein( reax_system *, control_params *,
        simulation_data *, static_storage *, reax_list **, output_controls * );

int Velocity_Verlet_Berendsen_Isotropic_NPT( reax_system *, control_params *,
        simulation_data *, static_storage *, reax_list **, output_controls * );

int Velocity_Verlet_Berendsen_Semi_Isotropic_NPT( reax_system *, control_params *,
        simulation_data *, static_storage *, reax_list **, output_controls * );

#if defined(ANISOTROPIC)
int Velocity_Verlet_Nose_Hoover_NVT( reax_system *, control_params *,
        simulation_data *, static_storage *, reax_list **, output_controls * );

int Velocity_Verlet_Isotropic_NPT( reax_system *, control_params *,
        simulation_data *, static_storage *, reax_list **, output_controls * );
#endif


#endif
