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

#ifndef __IO_TOOLS_H_
#define __IO_TOOLS_H_

#include "reax_types.h"


char *Get_Element( reax_system*, int );

char *Get_Atom_Name( reax_system*, int );

void Print_Near_Neighbors( reax_system*, control_params*, static_storage*,
        reax_list** );

void Print_Far_Neighbors( reax_system const * const, control_params const * const,
        simulation_data const * const, static_storage const * const,
        reax_list ** const );

void Print_Total_Force( reax_system*, control_params*, simulation_data*,
        static_storage*, reax_list**, output_controls* );

void Output_Results( reax_system*, control_params*, simulation_data*,
        static_storage*, reax_list**, output_controls* );

void Print_Bond_Orders( reax_system*, control_params*, simulation_data*,
        static_storage*, reax_list**, output_controls* );

void Print_Linear_System( reax_system*, control_params*, static_storage*, int );

void Print_Charges( reax_system*, control_params*, static_storage*, int );

void Print_Soln( static_storage*, real*, real*, real*, int );

void Print_Sparse_Matrix( sparse_matrix* );

void Print_Sparse_Matrix2( sparse_matrix*, char*, char* );

void Read_Sparse_Matrix2( sparse_matrix *, char * );

void Read_Permutation_Matrix( unsigned int *, char * );

void Print_Sparse_Matrix_Binary( sparse_matrix*, char* );

void Print_Bonds( reax_system*, reax_list*, char* );

void Print_Bond_List2( reax_system*, reax_list*, char* );

#if defined(TEST_FORCES)
void Dummy_Printer( reax_system*, control_params*, simulation_data*,
                    static_storage*, reax_list**, output_controls* );

void Print_Bond_Forces( reax_system*, control_params*, simulation_data*,
                        static_storage*, reax_list**, output_controls* );

void Print_LonePair_Forces( reax_system*, control_params*, simulation_data*,
                            static_storage*, reax_list**, output_controls* );

void Print_OverUnderCoor_Forces(reax_system*, control_params*, simulation_data*,
                                static_storage*, reax_list**, output_controls*);

void Print_Three_Body_Forces( reax_system*, control_params*, simulation_data*,
                              static_storage*, reax_list**, output_controls* );

void Print_Hydrogen_Bond_Forces(reax_system*, control_params*, simulation_data*,
                                static_storage*, reax_list**, output_controls*);

void Print_Four_Body_Forces( reax_system*, control_params*, simulation_data*,
                             static_storage*, reax_list**, output_controls* );

void Print_vdW_Coulomb_Forces( reax_system*, control_params*, simulation_data*,
                               static_storage*, reax_list**, output_controls* );

void Compare_Total_Forces( reax_system*, control_params*, simulation_data*,
                           static_storage*, reax_list**, output_controls* );

void Init_Force_Test_Functions( control_params* );
#endif


#endif
