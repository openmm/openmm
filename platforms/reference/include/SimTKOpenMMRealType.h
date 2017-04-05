
/* Portions copyright (c) 2006-2017 Stanford University and Simbios.
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

// LOG is used in Vishal's gpu code; modifying LOG -> LN 
#define LN             log

#define EXP            exp
#define FABS           fabs
#define ACOS           acos
#define ASIN           asin
#define ATAN           atan
#define TANH           tanh

#define FLOOR          floor

#define ATOF           atof

#define PI_M               3.141592653589
#define TWO_SIX            1.122462048309372981
#define RADIAN            57.29577951308
#define RADIAN_TO_DEGREE  57.29577951308
#define LOG_TEN            2.302585092994045684
#define SQRT_TWO           1.41421356237309504
#define DEGREE_TO_RADIAN   0.01745329252
#define RADIAN_INVERSE     0.01745329252

#define DOT3(u,v) ((u[0])*(v[0]) + (u[1])*(v[1]) + (u[2])*(v[2]))

#define MATRIXDOT3(u,v) u[0]*v[0] + u[1]*v[1] + u[2]*v[2] + \
                        u[3]*v[3] + u[4]*v[4] + u[5]*v[5] + \
                        u[6]*v[6] + u[7]*v[7] + u[8]*v[8]


// physics constants -- from Gromacs physics.h

#ifndef M_PI
#ifdef _PI
#define M_PI _PI
#else
#define M_PI        3.14159265358979323846
#endif
#endif

#define ANGSTROM     (1e-10) 
#define KILO         (1e3)
#define NANO         (1e-9)
#define PICO         (1e-12)
#define A2NM         (ANGSTROM/NANO)
#define NM2A         (NANO/ANGSTROM)
#define RAD2DEG      (180.0/M_PI)
// #define DEG2RAD      (M_PI/180.0)
#define CAL2JOULE    (4.184)
#define E_CHARGE     (1.60217733e-19)

#define AMU          (1.6605402e-27)
#define BOLTZMANN    (1.380658e-23)            /* (J/K)   */
#define AVOGADRO     (6.0221367e23)          
#define RGAS         (BOLTZMANN*AVOGADRO)      /* (J/(mol K))  */
#define BOLTZ        (RGAS/KILO)               /* (kJ/(mol K)) */
#define FARADAY      (E_CHARGE*AVOGADRO)       /* (C/mol)      */
#define ELECTRONVOLT (E_CHARGE*AVOGADRO/KILO)  /* (kJ/mol)   */     

#define EPSILON0     (5.72765E-4)              /* (e^2 Na/(kJ nm)) == (e^2/(kJ mol nm)) */ 

#define SPEED_OF_LIGHT   (2.9979245800E05)      /* nm/ps                */
#define ATOMICMASS_keV   (940000.0)             /* Atomic mass in keV   */
#define ELECTRONMASS_keV (512.0)                /* Electron mas in keV  */

#define ONE_4PI_EPS0      138.935456
#define PRESFAC           (16.6054)             /* bar / pressure unity */
#define ENM2DEBYE         48.0321               /* Convert electron nm to debye */
#define DEBYE2ENM         0.02081941

#endif // __RealSimTk_H__
