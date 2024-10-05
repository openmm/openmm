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

#include "vector.h"


void Vector_Print( FILE * const fout, const char * const vname, const real * const v,
        const unsigned int k )
{
    unsigned int i;

    if ( vname != NULL )
    {
        fprintf( fout, "%s:\n", vname );
    }

    for ( i = 0; i < k; ++i )
    {
        fprintf( fout, "%24.15e\n", v[i] );
    }
}


void Print_rTensor( FILE * const fp, rtensor t )
{
    unsigned int i, j;

    for (i = 0; i < 3; i++)
    {
        fprintf(fp, "[");
        for (j = 0; j < 3; j++)
        {
            fprintf(fp, "%8.3f,\t", t[i][j]);
        }
        fprintf(fp, "]\n");
    }
}
