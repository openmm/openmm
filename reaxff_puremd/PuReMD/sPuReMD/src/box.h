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

#ifndef __BOX_H__
#define __BOX_H__

#include "reax_types.h"


/* Computes all the transformations,
   metric and other quantities from box rtensor */
void Make_Consistent( simulation_box* );

void Setup_Box( real, real, real, real, real, real, simulation_box* );

/* Initializes box from box rtensor */
void Update_Box( rtensor, simulation_box* );
void Update_Box_Isotropic( simulation_box*, real );
void Update_Box_Semi_Isotropic( simulation_box*, rvec );

int Find_Non_Periodic_Far_Neighbors( rvec, rvec, int, int,
        simulation_box const * const, real, far_neighbor_data * const, int );

int Count_Non_Periodic_Far_Neighbors( rvec, rvec, int, int,
        simulation_box const * const , real );

int Find_Periodic_Far_Neighbors_Big_Box( rvec, rvec, int, int,
        simulation_box const * const, real, far_neighbor_data * const, int );

int Count_Periodic_Far_Neighbors_Big_Box( rvec, rvec, int, int,
        simulation_box const * const, real );

int Find_Periodic_Far_Neighbors_Small_Box( rvec, rvec, int, int,
        simulation_box const * const, real, far_neighbor_data * const, int );

int Count_Periodic_Far_Neighbors_Small_Box( rvec, rvec, int, int,
        simulation_box const * const, real );

void Compute_Atom_Distance_Triclinic( control_params *,
        simulation_box *, rvec, rvec, ivec, ivec, ivec, rvec );

void Update_Atom_Position_Triclinic( control_params *,
        simulation_box * const, rvec, rvec, ivec );

/* These functions assume that the coordinates are in triclinic system */
/* this function returns cartesian norm but triclinic distance vector */
real Compute_Atom_Distance_Periodic( simulation_box const * const, rvec, rvec, ivec, ivec, ivec, rvec );

real Compute_Atom_Distance_Non_Periodic( simulation_box const * const, rvec, rvec, ivec, ivec, ivec, rvec );

void Update_Atom_Position_Periodic( rvec, rvec, ivec, simulation_box const * const );

void Update_Atom_Position_Non_Periodic( rvec, rvec, ivec, simulation_box const * const );

real Metric_Product( rvec, rvec, simulation_box* );

void Print_Box( simulation_box*, FILE* );

#endif
