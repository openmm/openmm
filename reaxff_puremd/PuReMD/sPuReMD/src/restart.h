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

#ifndef __RESTART_H_
#define __RESTART_H_

#include "reax_types.h"

typedef struct
{
    int step, N;
    real T, xi, v_xi, v_xi_old, G_xi;
    rtensor box;
} restart_header;

typedef struct
{
    int orig_id, type;
    char name[9];
    rvec x, v;
} restart_atom;

#define RESTART_HEADER "%8d%12d%8.3f%8.3f%8.3f%8.3f%8.3f\n%15.5f%15.5f%15.5f\n%15.5f%15.5f%15.5f\n%15.5f%15.5f%15.5f\n"
#define RESTART_HEADER_LINE_LEN 200
/* step, system->bigN, data->therm.T, data->therm.xi,
   data->therm.v_xi data->therm.v_xi_old data->therm.G_xi
   system->big_box.box[0][0], [0][1], [0][2]
   system->big_box.box[1][0], [1][1], [1][2]
   system->big_box.box[2][0], [2][1], [2][2] */

#define RESTART_LINE "%6d%4d%8s%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f\n"
#define RESTART_LINE_LEN 109
/* id type name x y z vx vy vz */

#define READ_RESTART_HEADER " %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf"
#define READ_RESTART_LINE " %d %d %s %lf %lf %lf %lf %lf %lf"

void Write_Restart( reax_system *, control_params *,
        simulation_data *, static_storage *, output_controls * );

void Read_Binary_Restart( const char * const, reax_system *, control_params *,
        simulation_data *, static_storage * );

void Read_ASCII_Restart( const char * const, reax_system *, control_params *,
        simulation_data *, static_storage * );

#endif
