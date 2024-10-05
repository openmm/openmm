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

#ifndef __BOND_ORDERS_H_
#define __BOND_ORDERS_H_

#include "reax_types.h"


typedef struct
{
    real C1dbo, C2dbo, C3dbo;
    real C1dbopi, C2dbopi, C3dbopi, C4dbopi;
    real C1dbopi2, C2dbopi2, C3dbopi2, C4dbopi2;
    real C1dDelta, C2dDelta, C3dDelta;
} dbond_coefficients;


#if defined(TEST_FORCES)
void Add_dBO( reax_system*, reax_list**, int, int, real, rvec* );

void Add_dBOpinpi2( reax_system*, reax_list**, int, int, real, real, rvec*, rvec* );

void Add_dBO_to_Forces( reax_system*, reax_list**, int, int, real );

void Add_dBOpinpi2_to_Forces( reax_system*, reax_list**, int, int, real, real );

void Add_dDelta( reax_system*, reax_list**, int, real, rvec* );

void Add_dDelta_to_Forces( reax_system *, reax_list**, int, real );
#endif

void Add_dBond_to_Forces( int, int, reax_system*, simulation_data*,
        static_storage*, reax_list** );

void Add_dBond_to_Forces_NPT( int, int, reax_system*, simulation_data*,
        static_storage*, reax_list** );

void BO( reax_system*, control_params*, simulation_data*,
        static_storage*, reax_list**, output_controls* );


#endif
