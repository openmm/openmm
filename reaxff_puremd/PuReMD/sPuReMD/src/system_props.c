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

#include "system_props.h"

#include "tool_box.h"
#include "vector.h"


void Temperature_Control( control_params *control, simulation_data *data,
        output_controls *out_control )
{
    real tmp;

    /* step-wise temperature control */
    if ( control->T_mode == 1 )
    {
        if ( (data->step - data->prev_steps) %
                ((int)(control->T_freq / control->dt)) == 0 )
        {
            if ( FABS( control->T - control->T_final ) >= FABS( control->T_rate ) )
            {
                control->T += control->T_rate;
            }
            else
            {
                control->T = control->T_final;
            }
        }
    }
    /* constant slope control */
    else if ( control->T_mode == 2 )
    {
        tmp = control->T_rate * control->dt / control->T_freq;

        if ( FABS( control->T - control->T_final ) >= FABS( tmp ) )
        {
            control->T += tmp;
        }
    }
}


void Compute_Total_Mass( reax_system const * const system,
        simulation_data * const data )
{
    int i;

    data->M = 0.0;

    for ( i = 0; i < system->N; i++ )
    {
        data->M += system->reax_param.sbp[ system->atoms[i].type ].mass;
    }

    data->inv_M = 1.0 / data->M;
}


void Compute_Center_of_Mass( reax_system const * const system,
        simulation_data * const data )
{
    int i;
    real m, xx, xy, xz, yy, yz, zz, det;
    rvec tvec, diff;
    rtensor mat, inv;

    rvec_MakeZero( data->xcm );  // position of CoM
    rvec_MakeZero( data->vcm );  // velocity of CoM
    rvec_MakeZero( data->amcm ); // angular momentum of CoM
    rvec_MakeZero( data->avcm ); // angular velocity of CoM

    /* Compute the position, velocity and angular momentum about the CoM */
    for ( i = 0; i < system->N; ++i )
    {
        m = system->reax_param.sbp[ system->atoms[i].type ].mass;

        rvec_ScaledAdd( data->xcm, m, system->atoms[i].x );
        rvec_ScaledAdd( data->vcm, m, system->atoms[i].v );

        rvec_Cross( tvec, system->atoms[i].x, system->atoms[i].v );
        rvec_ScaledAdd( data->amcm, m, tvec );
    }

    rvec_Scale( data->xcm, data->inv_M, data->xcm );
    rvec_Scale( data->vcm, data->inv_M, data->vcm );

    rvec_Cross( tvec, data->xcm, data->vcm );
    rvec_ScaledAdd( data->amcm, -data->M, tvec );

    data->etran_cm = 0.5 * data->M * rvec_Norm_Sqr( data->vcm );

    /* Calculate and then invert the inertial tensor */
    xx = 0.0;
    xy = 0.0;
    xz = 0.0;
    yy = 0.0;
    yz = 0.0;
    zz = 0.0;

    for ( i = 0; i < system->N; ++i )
    {
        m = system->reax_param.sbp[ system->atoms[i].type ].mass;

        rvec_ScaledSum( diff, 1.0, system->atoms[i].x, -1.0, data->xcm );
        xx += diff[0] * diff[0] * m;
        xy += diff[0] * diff[1] * m;
        xz += diff[0] * diff[2] * m;
        yy += diff[1] * diff[1] * m;
        yz += diff[1] * diff[2] * m;
        zz += diff[2] * diff[2] * m;
    }

    mat[0][0] = yy + zz;
    mat[0][1] = -xy;
    mat[0][2] = -xz;
    mat[1][0] = -xy;
    mat[1][1] = xx + zz;
    mat[1][2] = -yz;
    mat[2][0] = -xz;
    mat[2][1] = -yz;
    mat[2][2] = xx + yy;

    /* invert the inertial tensor */
    det = ( mat[0][0] * mat[1][1] * mat[2][2]
            + mat[0][1] * mat[1][2] * mat[2][0]
            + mat[0][2] * mat[1][0] * mat[2][1] )
        - ( mat[0][0] * mat[1][2] * mat[2][1]
                + mat[0][1] * mat[1][0] * mat[2][2]
                + mat[0][2] * mat[1][1] * mat[2][0] );

    inv[0][0] = mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1];
    inv[0][1] = mat[0][2] * mat[2][1] - mat[0][1] * mat[2][2];
    inv[0][2] = mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1];
    inv[1][0] = mat[1][2] * mat[2][0] - mat[1][0] * mat[2][2];
    inv[1][1] = mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0];
    inv[1][2] = mat[0][2] * mat[1][0] - mat[0][0] * mat[1][2];
    inv[2][0] = mat[1][0] * mat[2][1] - mat[2][0] * mat[1][1];
    inv[2][1] = mat[2][0] * mat[0][1] - mat[0][0] * mat[2][1];
    inv[2][2] = mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1];

    if ( FABS(det) > ALMOST_ZERO )
    {
        rtensor_Scale( inv, 1. / det, inv );
    }
    else
    {
        rtensor_MakeZero( inv );
    }

    /* Compute the angular velocity about the centre of mass */
    rtensor_MatVec( data->avcm, inv, data->amcm );
    data->erot_cm = 0.5 * E_CONV * rvec_Dot( data->avcm, data->amcm );

#if defined(DEBUG_FOCUS)
    fprintf( stderr, "[INFO] xcm: %24.15e, %24.15e, %24.15e\n",
             data->xcm[0], data->xcm[1], data->xcm[2] );
    fprintf( stderr, "[INFO] vcm: %24.15e, %24.15e, %24.15e\n",
             data->vcm[0], data->vcm[1], data->vcm[2] );
    fprintf( stderr, "[INFO] amcm: %24.15e, %24.15e, %24.15e\n",
             data->amcm[0], data->amcm[1], data->amcm[2] );
    fprintf( stderr, "[INFO] avcm: %24.15e, %24.15e, %24.15e\n",
             data->avcm[0], data->avcm[1], data->avcm[2] );
#endif
}


void Compute_Kinetic_Energy( reax_system const * const  system,
        simulation_data * const data )
{
    int i;
    real m;
    rvec p;

    data->E_Kin = 0.0;

    for ( i = 0; i < system->N; i++ )
    {
        m = system->reax_param.sbp[system->atoms[i].type].mass;

        rvec_Scale( p, m, system->atoms[i].v );
        data->E_Kin += rvec_Dot( p, system->atoms[i].v );
    }
    data->E_Kin *= 0.5 * E_CONV;

    /* temperature scalar */
    data->therm.T = (2.0 * data->E_Kin) / (data->N_f * K_B / F_CONV);

    /* avoid T being an absolute zero! */
    if ( FABS( data->therm.T ) < ALMOST_ZERO )
    {
        data->therm.T = ALMOST_ZERO;
    }
}


/* Compute potential and total energies */
void Compute_Total_Energy( simulation_data* data )
{
    data->E_Pot = data->E_BE + data->E_Ov + data->E_Un  + data->E_Lp +
        data->E_Ang + data->E_Pen + data->E_Coa + data->E_HB +
        data->E_Tor + data->E_Con + data->E_vdW + data->E_Ele + data->E_Pol;

    data->E_Tot = data->E_Pot + data->E_Kin;
}


/* Check for numeric breakdowns in the energies */
void Check_Energy( simulation_data* data )
{
    if ( IS_NAN_REAL(data->E_Pol) )
    {
        fprintf( stderr, "[ERROR] NaN detected for polarization energy. Terminating...\n" );
        exit( NUMERIC_BREAKDOWN );
    }

    if ( IS_NAN_REAL(data->E_Pot) )
    {
        fprintf( stderr, "[ERROR] NaN detected for potential energy. Terminating...\n" );
        exit( NUMERIC_BREAKDOWN );
    }

    if ( IS_NAN_REAL(data->E_Tot) )
    {
        fprintf( stderr, "[ERROR] NaN detected for total energy. Terminating...\n" );
        exit( NUMERIC_BREAKDOWN );
    }
}


/* IMPORTANT: This function assumes that current kinetic energy and
 *  the center of mass of the system is already computed before.
 *
 * IMPORTANT: In Klein's paper, it is stated that a dU/dV term needs
 *  to be added when there are long-range interactions or long-range
 *  corrections to short-range interactions present.
 *  We may want to add that for more accuracy.
 */
void Compute_Pressure_Isotropic( reax_system* system, control_params *control,
        simulation_data* data, output_controls *out_control )
{
    int i;
    simulation_box *box;
    rtensor temp;

    box = &system->box;

    /* kinetic contribution */
    rtensor_MakeZero( data->kin_press );
    for ( i = 0; i < system->N; i++ )
    {
        rvec_OuterProduct( temp, system->atoms[i].v, system->atoms[i].v );
        rtensor_ScaledAdd( data->kin_press,
                system->reax_param.sbp[system->atoms[i].type].mass, temp );
    }
    /* unit conversion from mass * velocity^2 to energy */
    rtensor_Scale( data->kin_press, 48.88821291 * 48.88821291 / 1.0e7, data->kin_press );

    /* virial contribution */
#if defined(_OPENMP)
    for ( i = 0; i < control->num_threads; ++i )
    {
        rtensor_Add( data->press, data->press_local[i] );
    }
#endif
    rtensor_Scale( data->press, -1.0, data->press );

#if defined(DEBUG_FOCUS)
    fprintf( stderr, "[INFO] ke = (%12.6f, %12.6f, %12.6f), virial = (%12.6f, %12.6f, %12.6f)\n",
            data->kin_press[0][0], data->kin_press[1][1], data->kin_press[2][2],
            data->press[0][0], data->press[1][1], data->press[2][2] );
    fflush( stderr );
#endif

    /* to GPa */
//    rtensor_Scale( data->kin_press, 0.0166053907 / box->volume, data->kin_press );
//    rtensor_Scale( data->press, 0.0166053907 / box->volume, data->press );
    /* to ATM */
    rtensor_Scale( data->kin_press, 68568.415 / box->volume, data->kin_press );
    rtensor_Scale( data->press, 68568.415 / box->volume, data->press );

    /* total pressure in each direction, in GPa */
    data->tot_press[0] = data->kin_press[0][0] + data->press[0][0];

    data->tot_press[1] = data->kin_press[1][1] + data->press[1][1];

    data->tot_press[2] = data->kin_press[2][2] + data->press[2][2];

    /* average pressure for the whole box */
    data->iso_bar.P = (data->tot_press[0] + data->tot_press[1]
            + data->tot_press[2]) / 3.0;
}


/* IMPORTANT: This function assumes that current kinetic energy and
 *  the center of mass of the system is already computed before.
 *
 * IMPORTANT: In Klein's paper, it is stated that a dU/dV term needs
 *  to be added when there are long-range interactions or long-range
 *  corrections to short-range interactions present.
 *  We may want to add that for more accuracy.
 */
void Compute_Pressure_Isotropic_Klein( reax_system* system, simulation_data* data )
{
    int i;
    reax_atom *p_atom;
    rvec dx;

    data->iso_bar.P = 2.0 * data->E_Kin;

    for ( i = 0; i < system->N; ++i )
    {
        p_atom = &system->atoms[i];
        rvec_ScaledSum( dx, 1.0, p_atom->x, -1.0, data->xcm );
        data->iso_bar.P += -F_CONV * rvec_Dot( p_atom->f, dx );
    }

    data->iso_bar.P /= 3.0 * system->box.volume;
}


void Compute_Pressure( reax_system* system, simulation_data* data,
        static_storage *workspace )
{
    int i;
    reax_atom *p_atom;
    rtensor temp;

    rtensor_MakeZero( data->flex_bar.P );

    for ( i = 0; i < system->N; ++i )
    {
        p_atom = &system->atoms[i];

        rvec_OuterProduct( temp, p_atom->v, p_atom->v );
        rtensor_ScaledAdd( data->flex_bar.P,
                system->reax_param.sbp[ p_atom->type ].mass, temp );
//        rvec_OuterProduct( temp, workspace->virial_forces[i], p_atom->x );
        rtensor_ScaledAdd( data->flex_bar.P, -F_CONV, temp );
    }

    rtensor_Scale( data->flex_bar.P, 1.0 / system->box.volume, data->flex_bar.P );
    data->iso_bar.P = rtensor_Trace( data->flex_bar.P ) / 3.0;
}
