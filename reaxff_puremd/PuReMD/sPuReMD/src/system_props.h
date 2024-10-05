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

#ifndef __SYSTEM_PROP_H_
#define __SYSTEM_PROP_H_

#include "reax_types.h"


void Temperature_Control( control_params*, simulation_data*, output_controls* );

void Compute_Total_Mass( reax_system const * const, simulation_data * const );

void Compute_Center_of_Mass( reax_system const * const, simulation_data * const );

void Compute_Kinetic_Energy( reax_system const * const , simulation_data * const );

void Compute_Total_Energy( simulation_data* );

void Check_Energy( simulation_data* );

void Compute_Pressure_Isotropic( reax_system*, control_params*, simulation_data*, output_controls* );

void Compute_Pressure_Isotropic_Klein( reax_system*, simulation_data* );

void Compute_Pressure( reax_system*, simulation_data*, static_storage* );


#endif
