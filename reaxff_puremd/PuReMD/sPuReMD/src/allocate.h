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

#ifndef __ALLOCATE_H_
#define __ALLOCATE_H_

#include "reax_types.h"


void PreAllocate_Space( reax_system * const, control_params  const * const,
        static_storage * const, int );

void Allocate_Matrix( sparse_matrix * const , int, int, int );

void Deallocate_Matrix( sparse_matrix * const  );

void Initialize_HBond_List( int, int const * const, int * const,
        reax_list * const );

void Initialize_Bond_List( int * const, reax_list * const );

void Reallocate_Part1( reax_system * const, control_params const * const,
        static_storage * const, reax_list ** const );

void Reallocate_Part2( reax_system const * const, control_params const * const,
        simulation_data const * const, static_storage * const, reax_list ** const );


#endif
