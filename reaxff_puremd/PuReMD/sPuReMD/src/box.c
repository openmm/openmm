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

#include "box.h"

#include "tool_box.h"
#include "vector.h"

/* min. distance between neighboring atoms */
#define MIN_NBR_DIST (0.01)


void Make_Consistent( simulation_box* box )
{
    real one_vol;

    box->volume = box->box[0][0] * (box->box[1][1] * box->box[2][2]
            - box->box[2][1] * box->box[2][1])
        + box->box[0][1] * (box->box[2][0] * box->box[1][2]
                - box->box[1][0] * box->box[2][2])
        + box->box[0][2] * (box->box[1][0] * box->box[2][1]
                - box->box[2][0] * box->box[1][1]);

    one_vol = 1.0 / box->volume;

    box->box_inv[0][0] = (box->box[1][1] * box->box[2][2]
            - box->box[1][2] * box->box[2][1]) * one_vol;
    box->box_inv[0][1] = (box->box[0][2] * box->box[2][1]
            - box->box[0][1] * box->box[2][2]) * one_vol;
    box->box_inv[0][2] = (box->box[0][1] * box->box[1][2]
            - box->box[0][2] * box->box[1][1]) * one_vol;

    box->box_inv[1][0] = (box->box[1][2] * box->box[2][0]
            - box->box[1][0] * box->box[2][2]) * one_vol;
    box->box_inv[1][1] = (box->box[0][0] * box->box[2][2]
            - box->box[0][2] * box->box[2][0]) * one_vol;
    box->box_inv[1][2] = (box->box[0][2] * box->box[1][0]
            - box->box[0][0] * box->box[1][2]) * one_vol;

    box->box_inv[2][0] = (box->box[1][0] * box->box[2][1]
            - box->box[1][1] * box->box[2][0]) * one_vol;
    box->box_inv[2][1] = (box->box[0][1] * box->box[2][0]
            - box->box[0][0] * box->box[2][1]) * one_vol;
    box->box_inv[2][2] = (box->box[0][0] * box->box[1][1]
            - box->box[0][1] * box->box[1][0]) * one_vol;

    box->box_norms[0] = SQRT( SQR(box->box[0][0]) + SQR(box->box[0][1])
            + SQR(box->box[0][2]) );
    box->box_norms[1] = SQRT( SQR(box->box[1][0]) + SQR(box->box[1][1])
            + SQR(box->box[1][2]) );
    box->box_norms[2] = SQRT( SQR(box->box[2][0]) + SQR(box->box[2][1])
            + SQR(box->box[2][2]) );

    box->max[0] = box->min[0] + box->box_norms[0];
    box->max[1] = box->min[1] + box->box_norms[1];
    box->max[2] = box->min[2] + box->box_norms[2];

    box->trans[0][0] = box->box[0][0] / box->box_norms[0];
    box->trans[0][1] = box->box[1][0] / box->box_norms[0];
    box->trans[0][2] = box->box[2][0] / box->box_norms[0];

    box->trans[1][0] = box->box[0][1] / box->box_norms[1];
    box->trans[1][1] = box->box[1][1] / box->box_norms[1];
    box->trans[1][2] = box->box[2][1] / box->box_norms[1];

    box->trans[2][0] = box->box[0][2] / box->box_norms[2];
    box->trans[2][1] = box->box[1][2] / box->box_norms[2];
    box->trans[2][2] = box->box[2][2] / box->box_norms[2];

    one_vol = box->box_norms[0] * box->box_norms[1] * box->box_norms[2] * one_vol;

    box->trans_inv[0][0] = (box->trans[1][1] * box->trans[2][2]
            - box->trans[1][2] * box->trans[2][1]) * one_vol;
    box->trans_inv[0][1] = (box->trans[0][2] * box->trans[2][1]
            - box->trans[0][1] * box->trans[2][2]) * one_vol;
    box->trans_inv[0][2] = (box->trans[0][1] * box->trans[1][2]
            - box->trans[0][2] * box->trans[1][1]) * one_vol;

    box->trans_inv[1][0] = (box->trans[1][2] * box->trans[2][0]
            - box->trans[1][0] * box->trans[2][2]) * one_vol;
    box->trans_inv[1][1] = (box->trans[0][0] * box->trans[2][2]
            - box->trans[0][2] * box->trans[2][0]) * one_vol;
    box->trans_inv[1][2] = (box->trans[0][2] * box->trans[1][0]
            - box->trans[0][0] * box->trans[1][2]) * one_vol;

    box->trans_inv[2][0] = (box->trans[1][0] * box->trans[2][1]
            - box->trans[1][1] * box->trans[2][0]) * one_vol;
    box->trans_inv[2][1] = (box->trans[0][1] * box->trans[2][0]
            - box->trans[0][0] * box->trans[2][1]) * one_vol;
    box->trans_inv[2][2] = (box->trans[0][0] * box->trans[1][1]
            - box->trans[0][1] * box->trans[1][0]) * one_vol;

    box->g[0][0] = box->box[0][0] * box->box[0][0]
        + box->box[0][1] * box->box[0][1]
        + box->box[0][2] * box->box[0][2];
    box->g[1][0] =
        box->g[0][1] = box->box[0][0] * box->box[1][0] +
                       box->box[0][1] * box->box[1][1] +
                       box->box[0][2] * box->box[1][2];
    box->g[2][0] =
        box->g[0][2] = box->box[0][0] * box->box[2][0] +
                       box->box[0][1] * box->box[2][1] +
                       box->box[0][2] * box->box[2][2];

    box->g[1][1] = box->box[1][0] * box->box[1][0] +
                   box->box[1][1] * box->box[1][1] +
                   box->box[1][2] * box->box[1][2];
    box->g[1][2] =
        box->g[2][1] = box->box[1][0] * box->box[2][0] +
                       box->box[1][1] * box->box[2][1] +
                       box->box[1][2] * box->box[2][2];

    box->g[2][2] = box->box[2][0] * box->box[2][0] +
                   box->box[2][1] * box->box[2][1] +
                   box->box[2][2] * box->box[2][2];

    /* These proportions are only used for isotropic_NPT! */
    box->side_prop[0] = box->box[0][0] / box->box[0][0];
    box->side_prop[1] = box->box[1][1] / box->box[0][0];
    box->side_prop[2] = box->box[2][2] / box->box[0][0];
}


/* setup the simulation box */
void Setup_Box( real a, real b, real c, real alpha, real beta, real gamma,
        simulation_box* box )
{
    double c_alpha, c_beta, c_gamma, s_gamma, zi;

    if ( IS_NAN_REAL(a) || IS_NAN_REAL(b) || IS_NAN_REAL(c)
            || IS_NAN_REAL(alpha) || IS_NAN_REAL(beta) || IS_NAN_REAL(gamma) )
    {
        fprintf( stderr, "[ERROR] Invalid simulation box boundaries for big box (NaN). Terminating...\n" );
        exit( INVALID_INPUT );
    }

    c_alpha = COS( DEG2RAD(alpha) );
    c_beta  = COS( DEG2RAD(beta) );
    c_gamma = COS( DEG2RAD(gamma) );
    s_gamma = SIN( DEG2RAD(gamma) );
    zi = (c_alpha - c_beta * c_gamma) / s_gamma;

    box->box[0][0] = a;
    box->box[0][1] = 0.0;
    box->box[0][2] = 0.0;
    box->box[1][0] = b * c_gamma;
    box->box[1][1] = b * s_gamma;
    box->box[1][2] = 0.0;
    box->box[2][0] = c * c_beta;
    box->box[2][1] = c * zi;
    box->box[2][2] = c * SQRT(1.0 - SQR(c_beta) - SQR(zi));

#if defined(DEBUG_FOCUS)
    fprintf( stderr, "[INFO] Setup_Box: box is %8.2f x %8.2f x %8.2f\n",
             box->box[0][0], box->box[1][1], box->box[2][2] );
#endif

    rvec_MakeZero( box->min );
    Make_Consistent( box );
}


void Update_Box( rtensor box_tensor, simulation_box* box )
{
    int i, j;

    for ( i = 0; i < 3; i++ )
    {
        for ( j = 0; j < 3; j++ )
        {
            box->box[i][j] = box_tensor[i][j];
        }
    }

    Make_Consistent( box );
}


void Update_Box_Isotropic( simulation_box *box, real mu )
{
//    box->box[0][0] = POW( V_new
//            / ( box->side_prop[1] * box->side_prop[2] ), 1.0/3.0 );
//    box->box[1][1] = box->box[0][0] * box->side_prop[1];
//    box->box[2][2] = box->box[0][0] * box->side_prop[2];
    rtensor_Copy( box->old_box, box->box );
    box->box[0][0] *= mu;
    box->box[1][1] *= mu;
    box->box[2][2] *= mu;

    box->volume = box->box[0][0] * box->box[1][1] * box->box[2][2];
    Make_Consistent( box );
}


void Update_Box_Semi_Isotropic( simulation_box *box, rvec mu )
{
//    box->box[0][0] = POW( V_new
//            / ( box->side_prop[1] * box->side_prop[2] ), 1.0/3.0 );
//    box->box[1][1] = box->box[0][0] * box->side_prop[1];
//    box->box[2][2] = box->box[0][0] * box->side_prop[2];
    rtensor_Copy( box->old_box, box->box );
    box->box[0][0] *= mu[0];
    box->box[1][1] *= mu[1];
    box->box[2][2] *= mu[2];

    box->volume = box->box[0][0] * box->box[1][1] * box->box[2][2];
    Make_Consistent( box );
}


/* Update the atomic position such that the new position
 * is within the simulation box
 *
 * Assumption: (0, 0, 0) is the minimum box coordinate,
 * the maximum coordinate is strictly positive in each dimension
 *
 * Inputs:
 *  x: current atom position 
 *  dx: displacement vector from the current position
 *  box: struct of info. for the simulation box
 *
 * Outputs:
 *  x: updated atom position
 *  rel_map: relative coordinates for the peroidic image of
 *      the simulation box that this atom was re-mapped into
 *      during the position update where (0,0,0) represents no mapping
 *  */
void Update_Atom_Position_Periodic( rvec x, rvec dx, ivec rel_map, simulation_box  const * const box )
{
    int i, remapped;
    real tmp;

    for ( i = 0; i < 3; i++ )
    {
        /* new atomic position after update along dimension i */
        tmp = x[i] + dx[i];

        /* outside of simulation box boundary [0.0, d_i),
         * where d_i is the simulation box length along dimesion i */
        if ( tmp < 0.0 || tmp >= box->box_norms[i] )
        {
            remapped = FALSE;

            /* re-map the position to be in the range (-d_i, d_i) */
            if ( tmp <= -1.0 * box->box_norms[i] || tmp >= box->box_norms[i] )
            {
                rel_map[i] += (int) (tmp / box->box_norms[i]);
                tmp = FMOD( tmp, box->box_norms[i] );
                remapped = TRUE;
            }

            /* if the new position after re-mapping is negative,
             * translate back to the range [0, d_i) (i.e., inside the box) */
            if ( tmp < 0.0 )
            {
                tmp += box->box_norms[i];

                if ( remapped == FALSE )
                {
                    --rel_map[i];
                }
            }
        }

        assert( tmp >= 0.0 && tmp < box->box_norms[i] );

        x[i] = tmp;
    }
}


/* Update the atomic position such that the new position
 * is within the simulation box
 *
 * Assumption: (0, 0, 0) is the minimum box coordinate,
 * the the maximum coordinate is strictly positive in each dimension
 *
 * Inputs:
 *  x: current atom position 
 *  dx: displacement vector from the current position
 *  box: struct of info. for the simulation box
 *
 * Outputs:
 *  x: updated atom position
 *  rel_map: unused
 *  */
void Update_Atom_Position_Non_Periodic( rvec x, rvec dx, ivec rel_map, simulation_box  const * const box )
{
    int i;
    real tmp;

    for ( i = 0; i < 3; i++ )
    {
        /* new atomic position after update along dimension i */
        tmp = x[i] + dx[i];

        /* implement as hard boundary via thresholding */
        if ( tmp < 0.0 )
        {
            tmp = 0.0;
        }
        else if ( tmp > box->box_norms[i] )
        {
            tmp = box->box_norms[i] - 1.0E-10;
        }

        x[i] = tmp;
    }
}


/* Compute the Euclidean distance between a pair of atoms (3D space).
 * This function supports peroidic boundary conditions.
 *
 * Inputs:
 * box: struct containing simulation box info
 * x1: position of first atom
 * x2: position of second atom
 * x2_rel_box: relative position in terms of peroidic images of the
 *  simulation box for the second atom
 * x1_rel_map / x2_rel_map: relative coordinates for the peroidic image of
 *  the simulation box that this atom was re-mapped into
 *  during the position update where (0,0,0) represents no mapping
 *
 * Outputs:
 * r: displacement vector betweeon the atoms
 * return value: distance
 */
real Compute_Atom_Distance_Periodic( simulation_box const * const box,
        rvec x1, rvec x2, ivec x1_rel_map, ivec x2_rel_map, ivec x2_rel_box,
        rvec r )
{
    int i;
    real norm;

    norm = 0.0;

    for ( i = 0; i < 3; i++ )
    {
        /* translate the position of atom 2 (which has been fitted to be within
         * the simulation box boundaries) to the position of it's periodic image
         * (if necessary -- defined by rel_box)
         * and compute the component-wise difference with atom 1 */
        r[i] = (x2[i] + box->box_norms[i]
                * (x2_rel_box[i] + x2_rel_map[i] - x1_rel_map[i])) - x1[i];

        norm += SQR( r[i] );
    }

    /* note: distance of exactly 0.0 is also an error condition */
    assert( norm > 0.0 );

    return SQRT( norm );
}


/* Compute the Euclidean distance between a pair of atoms (3D space).
 * This function supports non-peroidic boundary conditions.
 *
 * Inputs:
 * box: struct containing simulation box info
 * x1: position of first atom
 * x2: position of second atom
 * x2_rel_box: unused
 * x1_rel_map / x2_rel_map: unused
 *
 * Outputs:
 * r: displacement vector betweeon the atoms
 * return value: distance
 */
real Compute_Atom_Distance_Non_Periodic( simulation_box const * const box, rvec x1, rvec x2,
        ivec x1_rel_map, ivec x2_rel_map, ivec x2_rel_box, rvec r )
{
    int i;
    real norm;

    norm = 0.0;

    for ( i = 0; i < 3; i++ )
    {
        r[i] = x2[i] - x1[i];

        norm += SQR( r[i] );
    }

    /* note: distance of exactly 0.0 is also an error condition */
    assert( norm > 0.0 );

    return SQRT( norm );
}


void Compute_Atom_Distance_Triclinic( control_params *control,
        simulation_box* box, rvec x1, rvec x2,
        ivec x1_rel_map, ivec x2_rel_map, ivec x2_rel_box, rvec r )
{
    rvec xa, xb, ra;

    /* translate positions from Cartesian to Triclinic coordinates */
    Transform( x1, box, -1, xa );
    Transform( x2, box, -1, xb );

    control->compute_atom_distance( box, xa, xb,
            x1_rel_map, x2_rel_map, x2_rel_box, ra );

    /* translate position from Triclinic to Cartesian coordinates */
    Transform( ra, box, 1, r );
}


void Update_Atom_Position_Triclinic( control_params *control, simulation_box * const box,
        rvec x, rvec dx, ivec rel_map )
{
    rvec xa, dxa;

    /* translate positions from Cartesian to Triclinic coordinates */
    Transform( x, box, -1, xa );
    Transform( dx, box, -1, dxa );

    control->update_atom_position( xa, dxa, rel_map, box );

    /* translate position from Triclinic to Cartesian coordinates */
    Transform( xa, box, 1, x );
}


real Metric_Product( rvec x1, rvec x2, simulation_box* box )
{
    int i, j;
    real dist, tmp;

    dist = 0.0;

    for ( i = 0; i < 3; i++ )
    {
        tmp = 0.0;

        for ( j = 0; j < 3; j++ )
        {
            tmp += box->g[i][j] * x2[j];
        }

        dist += x1[i] * tmp;
    }

    return dist;
}


/* Determines if the distance between atoms x1 and x2 is strictly less than
 * vlist_cut.  If so, this neighborhood is added to the list of far neighbors.
 * Note: Periodic boundary conditions do not apply. */
int Find_Non_Periodic_Far_Neighbors( rvec x1, rvec x2, int atom, int nbr_atom,
        simulation_box const * const box, real vlist_cut,
        far_neighbor_data * const data, int max )
{
    int count;
    real norm_sqr;
    rvec dvec;

    rvec_ScaledSum( dvec, 1.0, x2, -1.0, x1 );
    norm_sqr = rvec_Norm_Sqr( dvec );

    if ( norm_sqr <= SQR( vlist_cut ) && norm_sqr >= MIN_NBR_DIST )
    {
        if ( max > 0 )
        {
            data->nbr = nbr_atom;
            ivec_MakeZero( data->rel_box );
//            rvec_MakeZero( data->ext_factor );
            data->d = SQRT( norm_sqr );
            rvec_Copy( data->dvec, dvec );
        }
        count = 1;
    }
    else
    {
        count = 0;
    }

    return count;
}


/* Similar to Find_Non_Periodic_Far_Neighbors but does not
 * update the far neighbors list */
int Count_Non_Periodic_Far_Neighbors( rvec x1, rvec x2, int atom, int nbr_atom,
        simulation_box const * const box, real vlist_cut )
{
    int count;
    real norm_sqr;
    rvec d;

    rvec_ScaledSum( d, 1.0, x2, -1.0, x1 );
    norm_sqr = rvec_Norm_Sqr( d );

    if ( norm_sqr <= SQR( vlist_cut ) && norm_sqr >= MIN_NBR_DIST )
    {
        count = 1;
    }
    else
    {
        count = 0;
    }

    return count;
}


/* Finds periodic neighbors in a 'big_box'. Here 'big_box' means the current
 * simulation box has all dimensions strictly greater than twice of vlist_cut.
 * If the periodic distance between x1 and x2 is less than vlist_cut, this
 * neighborhood is added to the list of far neighbors. */
int Find_Periodic_Far_Neighbors_Big_Box( rvec x1, rvec x2, int atom, int nbr_atom,
        simulation_box const * const box, real vlist_cut,
        far_neighbor_data * const data, int max )
{
    int i, count;
    real norm_sqr, tmp, sqr_vlist_cut;
    ivec rel_box;
//    rvec ext_factor;
    rvec dvec;

    norm_sqr = 0.0;
    sqr_vlist_cut = SQR( vlist_cut );

    for ( i = 0; i < 3; i++ )
    {
        dvec[i] = x2[i] - x1[i];
        tmp = SQR( dvec[i] );

        if ( tmp >= SQR( box->box_norms[i] / 2.0 ) )
        {
            if ( x2[i] > x1[i] )
            {
                rel_box[i] = -1;
//                ext_factor[i] = +1.0;
                dvec[i] -= box->box_norms[i];
            }
            else
            {
                rel_box[i] = +1;
//                ext_factor[i] = -1.0;
                dvec[i] += box->box_norms[i];
            }

            norm_sqr += SQR( dvec[i] );
        }
        else
        {
            norm_sqr += tmp;
            rel_box[i] = 0;
//            ext_factor[i] = 0.0;
        }
    }

    if ( norm_sqr <= sqr_vlist_cut && norm_sqr >= MIN_NBR_DIST )
    {
        if ( max > 0 )
        {
            data->nbr = nbr_atom;
            ivec_Copy( data->rel_box, rel_box );
//            rvec_Copy( data->ext_factor, ext_factor );
            data->d = SQRT( norm_sqr );
            rvec_Copy( data->dvec, dvec );
        }
        count = 1;
    }
    else
    {
        count = 0;
    }

    return count;
}


/* Similar to Find_Periodic_Far_Neighbors_Big_Box but does not
 * update the far neighbors list */
int Count_Periodic_Far_Neighbors_Big_Box( rvec x1, rvec x2, int atom, int nbr_atom,
        simulation_box const * const box, real vlist_cut )
{
    int i, count;
    real norm_sqr, d, tmp;

    norm_sqr = 0.0;

    for ( i = 0; i < 3; i++ )
    {
        d = x2[i] - x1[i];
        tmp = SQR( d );

        if ( tmp >= SQR( box->box_norms[i] / 2.0 ) )
        {
            if ( x2[i] > x1[i] )
            {
                d -= box->box_norms[i];
            }
            else
            {
                d += box->box_norms[i];
            }

            norm_sqr += SQR( d );
        }
        else
        {
            norm_sqr += tmp;
        }
    }

    if ( norm_sqr <= SQR( vlist_cut ) && norm_sqr >= MIN_NBR_DIST )
    {
        count = 1;
    }
    else
    {
        count = 0;
    }

    return count;
}


/* Finds all periodic far neighborhoods between atoms x1 and x2
 * ((dist(x1, x2') < 2 * vlist_cut, periodic images of x2 are also considered).
 * Here the box is 'small' meaning that at least one dimension is strictly
 * less than twice of vlist_cut.
 *
 * NOTE: This part might need some improvement. In NPT, the simulation box
 * might get too small (such as <5 A!). In this case we have to consider the
 * periodic images of x2 that are two boxs away!!! */
int Find_Periodic_Far_Neighbors_Small_Box( rvec x1, rvec x2, int atom, int nbr_atom,
        simulation_box const * const box, real vlist_cut,
        far_neighbor_data * const data, int max )
{
    int i, j, k, count;
    int imax, jmax, kmax;
    real sqr_norm, d_i, d_j, d_k, sqr_vlist_cut;
    real sqr_d_i, sqr_d_j, sqr_d_k;

    count = 0;
    sqr_vlist_cut = SQR( vlist_cut );

    /* determine the max. num. of imaginary boxs in each direction
     * to handle periodic boundary conditions correctly */
    imax = (int) CEIL( vlist_cut / box->box_norms[0] );
    jmax = (int) CEIL( vlist_cut / box->box_norms[1] );
    kmax = (int) CEIL( vlist_cut / box->box_norms[2] );

    /* assumption: orthogonal coordinates */
    for ( i = -1 * imax; i <= imax; ++i )
    {
        d_i = (x2[0] + i * box->box_norms[0]) - x1[0];
        sqr_d_i = SQR( d_i );

        for ( j = -1 * jmax; j <= jmax; ++j )
        {
            d_j = (x2[1] + j * box->box_norms[1]) - x1[1];
            sqr_d_j = SQR( d_j );

            for ( k = -1 * kmax; k <= kmax; ++k )
            {
                d_k = (x2[2] + k * box->box_norms[2]) - x1[2];
                sqr_d_k = SQR( d_k );

                sqr_norm = sqr_d_i + sqr_d_j + sqr_d_k;

                if ( sqr_norm <= sqr_vlist_cut && sqr_norm >= MIN_NBR_DIST )
                {
                    if ( count < max )
                    {
                        data[count].nbr = nbr_atom;

                        data[count].rel_box[0] = i;
                        data[count].rel_box[1] = j;
                        data[count].rel_box[2] = k;

//                        if ( i )
//                        {
//                            data[count].ext_factor[0] = (real)i / -ABS(i);
//                        }
//                        else
//                        {
//                            data[count].ext_factor[0] = 0;
//                        }
//    
//                        if ( j )
//                        {
//                            data[count].ext_factor[1] = (real)j / -ABS(j);
//                        }
//                        else
//                        {
//                            data[count].ext_factor[1] = 0;
//                        }
//    
//                        if ( k )
//                        {
//                            data[count].ext_factor[2] = (real)k / -ABS(k);
//                        }
//                        else
//                        {
//                            data[count].ext_factor[2] = 0;
//                        }
//    
//                        if ( i == 0 && j == 0 && k == 0 )
//                        {
//                            data[count].imaginary = 0;
//                        }
//                        else
//                        {
//                            data[count].imaginary = 1;
//                        }

                        data[count].d = SQRT( sqr_norm );
                        assert( data[count].d > 0.0 );

                        data[count].dvec[0] = d_i;
                        data[count].dvec[1] = d_j;
                        data[count].dvec[2] = d_k;
                    }

                    ++count;
                }
            }
        }
    }

    return count;
}


/* Similar to Find_Periodic_Far_Neighbors_Small_Box but does not
 * update the far neighbors list */
int Count_Periodic_Far_Neighbors_Small_Box( rvec x1, rvec x2, int atom, int nbr_atom,
        simulation_box const * const box, real vlist_cut )
{
    int i, j, k, count;
    int imax, jmax, kmax;
    real sqr_norm, d_i, d_j, d_k, sqr_vlist_cut;
    real sqr_d_i, sqr_d_j, sqr_d_k;

    count = 0;
    sqr_vlist_cut = SQR( vlist_cut );

    /* determine the max. num. of imaginary boxs in each direction
     * to handle periodic boundary conditions correctly */
    imax = (int) CEIL( vlist_cut / box->box_norms[0] );
    jmax = (int) CEIL( vlist_cut / box->box_norms[1] );
    kmax = (int) CEIL( vlist_cut / box->box_norms[2] );

    for ( i = -imax; i <= imax; ++i )
    {
        d_i = (x2[0] + i * box->box_norms[0]) - x1[0];
        sqr_d_i = SQR( d_i );

        for ( j = -jmax; j <= jmax; ++j )
        {
            d_j = (x2[1] + j * box->box_norms[1]) - x1[1];
            sqr_d_j = SQR( d_j );

            for ( k = -kmax; k <= kmax; ++k )
            {
                d_k = (x2[2] + k * box->box_norms[2]) - x1[2];
                sqr_d_k = SQR( d_k );

                sqr_norm = sqr_d_i + sqr_d_j + sqr_d_k;

                if ( sqr_norm <= sqr_vlist_cut && sqr_norm >= MIN_NBR_DIST )
                {
                    ++count;
                }
            }
        }
    }

    return count;
}


void Print_Box( simulation_box* box, FILE *out )
{
    int i, j;

    fprintf( out, "box: {" );
    for ( i = 0; i < 3; ++i )
    {
        fprintf( out, "{" );
        for ( j = 0; j < 3; ++j )
        {
            fprintf( out, "%8.3f ", box->box[i][j] );
        }
        fprintf( out, "}" );
    }
    fprintf( out, "}\n" );

    fprintf( out, "V: %8.3f\tdims: {%8.3f, %8.3f, %8.3f}\n",
             box->volume,
             box->box_norms[0], box->box_norms[1], box->box_norms[2] );

    fprintf( out, "box_trans: {" );
    for ( i = 0; i < 3; ++i )
    {
        fprintf( out, "{" );
        for ( j = 0; j < 3; ++j )
        {
            fprintf( out, "%8.3f ", box->trans[i][j] );
        }
        fprintf( out, "}" );
    }
    fprintf( out, "}\n" );

    fprintf( out, "box_trinv: {" );
    for ( i = 0; i < 3; ++i )
    {
        fprintf( out, "{" );
        for ( j = 0; j < 3; ++j )
        {
            fprintf( out, "%8.3f ", box->trans_inv[i][j] );
        }
        fprintf( out, "}" );
    }
    fprintf( out, "}\n" );
}
