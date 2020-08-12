
/* Portions copyright (c) 2006-2020 Stanford University and Simbios.
 * Contributors: Pande Group
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */


#ifndef __RealSimTk_H__
#define __RealSimTk_H__

#include <cmath>

#define RealOpenMM     double
#define SQRT           sqrt
#define POW            pow
#define SIN            sin
#define COS            cos
#define TAN            tan
#define EXP            exp
#define FABS           fabs
#define ACOS           acos
#define ASIN           asin
#define ATAN           atan
#define TANH           tanh

#define FLOOR          floor

#define ATOF           atof

#ifndef M_PI
#ifdef _PI
#define M_PI _PI
#else
#define M_PI        3.14159265358979323846
#endif
#endif

#define PI_M               M_PI
#define RADIAN             (180/M_PI)
#define SQRT_TWO           1.41421356237309504

#define DOT3(u,v) ((u[0])*(v[0]) + (u[1])*(v[1]) + (u[2])*(v[2]))


// Physical constants are CODATA 2018 values from https://pml.nist.gov/cuu/Constants

#define ANGSTROM     (1e-10) 
#define KILO         (1e3)
#define NANO         (1e-9)
#define PICO         (1e-12)
#define A2NM         (ANGSTROM/NANO)
#define NM2A         (NANO/ANGSTROM)
#define RAD2DEG      (180.0/M_PI)
// #define DEG2RAD      (M_PI/180.0)
#define CAL2JOULE    (4.184)
#define E_CHARGE     (1.602176634e-19)

#define BOLTZMANN    (1.380649e-23)            /* (J/K)   */
#define AVOGADRO     (6.02214076e23)          
#define AMU          (1/(KILO*AVOGADRO))
#define RGAS         (BOLTZMANN*AVOGADRO)      /* (J/(mol K))  */
#define BOLTZ        (RGAS/KILO)               /* (kJ/(mol K)) */
#define FARADAY      (E_CHARGE*AVOGADRO)       /* (C/mol)      */
#define ELECTRONVOLT (E_CHARGE*AVOGADRO/KILO)  /* (kJ/mol)   */     

#define EPSILON0     (1e-6*8.8541878128e-12/(E_CHARGE*E_CHARGE*AVOGADRO)) /* (e^2 Na/(kJ nm)) == (e^2/(kJ mol nm)) */ 

#define SPEED_OF_LIGHT   (2.9979245800E05)      /* nm/ps                */
#define ELECTRONMASS_keV (510.998950)           /* Electron mas in keV  */

#define ONE_4PI_EPS0      (1/(4*M_PI*EPSILON0))
#define PRESFAC           (16.6054)             /* bar / pressure unity */
#define ENM2DEBYE         48.0321               /* Convert electron nm to debye */
#define DEBYE2ENM         0.02081941

#endif // __RealSimTk_H__
