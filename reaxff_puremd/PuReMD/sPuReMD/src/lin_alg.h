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

#ifndef __LIN_ALG_H_
#define __LIN_ALG_H_

#include "reax_types.h"

typedef enum
{
    LOWER = 0,
    UPPER = 1,
} TRIANGULARITY;

void Sort_Matrix_Rows( sparse_matrix * const );

void setup_sparse_approx_inverse( const sparse_matrix * const, sparse_matrix *, 
        sparse_matrix *, sparse_matrix *, sparse_matrix *, const real, int );

int Estimate_LU_Fill( const sparse_matrix * const, const real * const );

void Calculate_Droptol( const sparse_matrix * const,
        real * const, const real );

real jacobi( const sparse_matrix * const, real * const, int, int );

real ICHOLT( const sparse_matrix * const, const real * const,
        sparse_matrix * const, sparse_matrix * const );

real ILUT( const sparse_matrix * const, const real * const,
        sparse_matrix * const, sparse_matrix * const );

real ILUTP( const sparse_matrix * const, const real * const,
        sparse_matrix * const, sparse_matrix * const );

real FG_ICHOLT( const sparse_matrix * const, const real *,
        const unsigned int, sparse_matrix * const, sparse_matrix * const );

real FG_ILUT( const sparse_matrix * const, const real *,
        const unsigned int, sparse_matrix * const, sparse_matrix * const );

real ILU( const sparse_matrix * const, sparse_matrix * const,
        sparse_matrix * const );

#if defined(HAVE_LAPACKE) || defined(HAVE_LAPACKE_MKL)
real sparse_approx_inverse( const sparse_matrix * const, const sparse_matrix * const,
        sparse_matrix * );
#endif

void Transpose( const sparse_matrix * const, sparse_matrix * const );

void Transpose_I( sparse_matrix * const );

void tri_solve( const sparse_matrix * const, const real * const,
        real * const, const TRIANGULARITY );

void tri_solve_level_sched( static_storage *,
        const sparse_matrix * const,
        const real * const, real * const,
        const TRIANGULARITY, int );

void jacobi_iter( const static_storage * const,
        const sparse_matrix * const, const real * const,
        const real * const, real * const, const TRIANGULARITY,
        const unsigned int );

void setup_graph_coloring( const control_params * const,
        const static_storage * const, const sparse_matrix * const,
        sparse_matrix *, sparse_matrix *, int );

int GMRES( const static_storage * const, const control_params * const,
        simulation_data * const, const sparse_matrix * const,
        const real * const, const real, real * const,
        const int );

int GMRES_HouseHolder( const static_storage * const, const control_params * const,
        simulation_data * const, const sparse_matrix * const,
        const real * const, const real, real * const,
        const int );

int CG( const static_storage * const, const control_params * const,
        simulation_data * const, const sparse_matrix * const, const real * const,
        const real, real * const, const int );

int BiCGStab( const static_storage * const, const control_params * const,
        simulation_data * const, const sparse_matrix * const, const real * const,
        const real, real * const, const int );

int SDM( const static_storage * const, const control_params * const,
        simulation_data * const, const sparse_matrix * const, const real * const, const real,
        real * const, const int );

real condest( const sparse_matrix * const, const sparse_matrix * const );


#endif
