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

#ifndef __CONTROL_H_
#define __CONTROL_H_

#include "reax_types.h"

int Set_Control_Parameter( const char * const, const char ** const,
        control_params * const, output_controls * const );

void Set_Control_Defaults( reax_system * const, control_params * const,
        output_controls * const );

void Read_Control_File( const char * const, reax_system * const, control_params * const,
        output_controls * const );

void Set_Control_Derived_Values( reax_system * const, control_params * const );

#endif
