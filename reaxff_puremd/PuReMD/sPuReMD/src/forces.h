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

#ifndef __FORCES_H_
#define __FORCES_H_

#include "reax_types.h"


void Init_Bonded_Force_Functions( control_params * const );

void Estimate_Storages( reax_system const * const, control_params const * const,
        static_storage * const, reax_list ** const, int, int );

int Compute_Forces( reax_system * const, control_params * const,
        simulation_data * const, static_storage * const, reax_list ** const,
        output_controls * const, int );


#endif
