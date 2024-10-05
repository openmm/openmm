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

#include <stdio.h>
#include <stdlib.h>

#include "spuremd.h"

#define INVALID_INPUT (-1)


static void usage( char * argv[] )
{
    fprintf( stderr, "usage: ./%s geo_file ffield_param_file control_file [geo_file2 ffield_param_file2 control_file2 ...]\n",
            argv[0] );
    fprintf( stderr, "  geo_file, geo_file2, ...                    Geometry files containing simulation information. Format\n"
            "                                              specified via geo_format keyword in respective simulation control file\n" );
    fprintf( stderr, "  ffield_param_file, ffield_param_file2, ...  ReaxFF parameter files.\n" );
    fprintf( stderr, "  control_file, control_file2, ...            PuReMD control files for simulation parameters.\n" );
}


int main( int argc, char* argv[] )
{
    int i, ret;
    void *handle;

    if ( argc < 4 || argc % 3 != 1 )
    {
        usage( argv );
        exit( INVALID_INPUT );
    }

    handle = setup( argv[1], argv[2], argv[3] );
    ret = SPUREMD_FAILURE;

    if ( handle != NULL )
    {
        ret = simulate( handle );
    }

    for ( i = 1; i < argc / 3; ++i )
    {
        if ( ret == SPUREMD_SUCCESS )
        {
            ret = reset( handle, argv[3 * i + 1], argv[3 * i + 2], argv[3 * i + 3] );
        }

        if ( ret == SPUREMD_SUCCESS )
        {
            ret = simulate( handle );
        }
    }

    if ( ret == SPUREMD_SUCCESS )
    {
        ret = cleanup( handle );
    }

    return (ret == SPUREMD_SUCCESS) ? 0 : (-1);
}
