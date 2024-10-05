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

#ifndef __VECTOR_H_
#define __VECTOR_H_

#include "reax_types.h"

#include "random.h"


/* file scope to make OpenMP shared (Vector_isZero) */
static unsigned int ret_omp;
/* file scope to make OpenMP shared (Dot, Norm) */
static real ret2_omp;


void Vector_Print( FILE * const, const char * const, const real * const, const unsigned int );

void Print_rTensor( FILE * const, rtensor );


static inline int Vector_isZero( const real * const v, const unsigned int k )
{
    unsigned int i;

#if defined(_OPENMP)
    #pragma omp single
#endif
    {
        ret_omp = TRUE;
    }

#if defined(_OPENMP)
    #pragma omp for simd reduction(&&: ret_omp) schedule(simd:static)
#endif
    for ( i = 0; i < k; ++i )
    {
        if ( FABS( v[i] ) > ALMOST_ZERO )
        {
            ret_omp = FALSE;
        }
    }

    return ret_omp;
}


static inline void Vector_MakeZero( real * const v, const unsigned int k )
{
    unsigned int i;

#if defined(_OPENMP)
    #pragma omp for simd schedule(simd:static)
#endif
    for ( i = 0; i < k; ++i )
    {
        v[i] = 0.0;
    }
}


#if defined(QMMM)
static inline void Vector_Mask_qmmm( real * const v, int const * const mask,
        const unsigned int k )
{
    unsigned int i;

#if defined(_OPENMP)
    #pragma omp for simd schedule(simd:static)
#endif
    for ( i = 0; i < k; ++i )
    {
        v[i] = (mask[i] == TRUE ? v[i] : 0.0);
    }
}
#endif


static inline void Vector_Copy( real * const dest, const real * const v, const unsigned int k )
{
    unsigned int i;

#if defined(_OPENMP)
    #pragma omp for simd schedule(simd:static)
#endif
    for ( i = 0; i < k; ++i )
    {
        dest[i] = v[i];
    }
}


static inline void Vector_Scale( real * const dest, const real c, const real * const v,
        const unsigned int k )
{
    unsigned int i;

#if defined(_OPENMP)
    #pragma omp for simd schedule(simd:static)
#endif
    for ( i = 0; i < k; ++i )
    {
        dest[i] = c * v[i];
    }
}


static inline void Vector_Sum( real * const dest, const real c, const real * const v,
        const real d, const real * const y, const unsigned int k )
{
    unsigned int i;

#if defined(_OPENMP)
    #pragma omp for simd schedule(simd:static)
#endif
    for ( i = 0; i < k; ++i )
    {
        dest[i] = c * v[i] + d * y[i];
    }
}


static inline void Vector_Add( real * const dest, const real c, const real * const v,
        const unsigned int k )
{
    unsigned int i;

#if defined(_OPENMP)
    #pragma omp for simd schedule(simd:static)
#endif
    for ( i = 0; i < k; ++i )
    {
        dest[i] += c * v[i];
    }
}


static inline real Dot( const real * const v1, const real * const v2,
        const unsigned int k )
{
    unsigned int i;

#if defined(_OPENMP)
    #pragma omp single
#endif
    {
        ret2_omp = 0.0;
    }

#if defined(_OPENMP)
    #pragma omp for simd reduction(+: ret2_omp) schedule(simd:static)
#endif
    for ( i = 0; i < k; ++i )
    {
        ret2_omp += v1[i] * v2[i];
    }

    return ret2_omp;
}


static inline real Norm( const real * const v1, const unsigned int k )
{
    unsigned int i;

#if defined(_OPENMP)
    #pragma omp single
#endif
    {
        ret2_omp = 0.0;
    }

#if defined(_OPENMP)
    #pragma omp for simd reduction(+: ret2_omp) schedule(simd:static)
#endif
    for ( i = 0; i < k; ++i )
    {
        ret2_omp += SQR( v1[i] );
    }

#if defined(_OPENMP)
    #pragma omp single
#endif
    {
        ret2_omp = SQRT( ret2_omp );
    }

    return ret2_omp;
}


static inline void rvec_Copy( rvec dest, const rvec src )
{
    int i;

#if defined(_OPENMP)
    #pragma omp simd
#endif
    for ( i = 0; i < 3; ++i )
    {
        dest[i] = src[i];
    }
}

static inline void rvec_Scale( rvec ret, const real c, const rvec v )
{
    int i;

#if defined(_OPENMP)
    #pragma omp simd
#endif
    for ( i = 0; i < 3; ++i )
    {
        ret[i] = c * v[i];
    }
}


static inline void rvec_Add( rvec ret, const rvec v )
{
    int i;

#if defined(_OPENMP)
    #pragma omp simd
#endif
    for ( i = 0; i < 3; ++i )
    {
        ret[i] += v[i];
    }
}


static inline void rvec_ScaledAdd( rvec ret, const real c, const rvec v )
{
    int i;

#if defined(_OPENMP)
    #pragma omp simd
#endif
    for ( i = 0; i < 3; ++i )
    {
        ret[i] += c * v[i];
    }
}


static inline void rvec_Sum( rvec ret, const rvec v1 , const rvec v2 )
{
    int i;

#if defined(_OPENMP)
    #pragma omp simd
#endif
    for ( i = 0; i < 3; ++i )
    {
        ret[i] = v1[i] + v2[i];
    }
}


static inline void rvec_ScaledSum( rvec ret, const real c1, const rvec v1,
        const real c2, const rvec v2 )
{
    int i;

#if defined(_OPENMP)
    #pragma omp simd
#endif
    for ( i = 0; i < 3; ++i )
    {
        ret[i] = c1 * v1[i] + c2 * v2[i];
    }
}


static inline real rvec_Dot( const rvec v1, const rvec v2 )
{
    int i;
    real ret;

    ret = 0.0;

#if defined(_OPENMP)
    #pragma omp simd reduction(+: ret)
#endif
    for ( i = 0; i < 3; ++i )
    {
        ret += v1[i] * v2[i];
    }

    return ret;
}


static inline real rvec_ScaledDot( const real c1, const rvec v1,
        const real c2, const rvec v2 )
{
    int i;
    real ret;

    ret = 0.0;

#if defined(_OPENMP)
    #pragma omp simd reduction(+: ret)
#endif
    for ( i = 0; i < 3; ++i )
    {
        ret += v1[i] * v2[i];
    }

    return c1 * c2 * ret;
}


static inline void rvec_Multiply( rvec r, const rvec v1, const rvec v2 )
{
    int i;

#if defined(_OPENMP)
    #pragma omp simd
#endif
    for ( i = 0; i < 3; ++i )
    {
        r[i] = v1[i] * v2[i];
    }
}


static inline void rvec_iMultiply( rvec r, const ivec v1, const rvec v2 )
{
    int i;

#if defined(_OPENMP)
    #pragma omp simd
#endif
    for ( i = 0; i < 3; ++i )
    {
        r[i] = v1[i] * v2[i];
    }
}


static inline void rvec_Divide( rvec r, const rvec v1, const rvec v2 )
{
    int i;

#if defined(_OPENMP)
    #pragma omp simd
#endif
    for ( i = 0; i < 3; ++i )
    {
        r[i] = v1[i] / v2[i];
    }
}


static inline void rvec_iDivide( rvec r, const rvec v1, const ivec v2 )
{
    int i;

#if defined(_OPENMP)
    #pragma omp simd
#endif
    for ( i = 0; i < 3; ++i )
    {
        r[i] = v1[i] / v2[i];
    }
}


static inline void rvec_Invert( rvec r, const rvec v )
{
    int i;

#if defined(_OPENMP)
    #pragma omp simd
#endif
    for ( i = 0; i < 3; ++i )
    {
        r[i] = 1.0 / v[i];
    }
}


static inline void rvec_Cross( rvec ret, const rvec v1, const rvec v2 )
{
    int i;

#if defined(_OPENMP)
    #pragma omp simd
#endif
    for ( i = 0; i < 3; ++i )
    {
        ret[i] = v1[(i + 1) % 3] * v2[(i + 2) % 3]
            - v1[(i + 2) % 3] * v2[(i + 1) % 3];
    }
}


static inline void rvec_OuterProduct( rtensor r, const rvec v1, const rvec v2 )
{
    unsigned int i, j;

    for ( i = 0; i < 3; ++i )
    {
        for ( j = 0; j < 3; ++j )
        {
            r[i][j] = v1[i] * v2[j];
        }
    }
}


static inline real rvec_Norm_Sqr( const rvec v )
{
    int i;
    real ret;

    ret = 0.0;

#if defined(_OPENMP)
    #pragma omp simd reduction(+: ret)
#endif
    for ( i = 0; i < 3; ++i )
    {
        ret += SQR( v[i] );
    }

    return ret;
}


static inline real rvec_Norm( const rvec v )
{
    int i;
    real ret;

    ret = 0.0;

#if defined(_OPENMP)
    #pragma omp simd reduction(+: ret)
#endif
    for ( i = 0; i < 3; ++i )
    {
        ret += SQR( v[i] );
    }

    return SQRT( ret );
}


static inline int rvec_isZero( const rvec v )
{
    int i, ret;

    ret = TRUE;

#if defined(_OPENMP)
    #pragma omp simd reduction(&&: ret)
#endif
    for ( i = 0; i < 3; ++i )
    {
        ret = (FABS( v[i] ) <= ALMOST_ZERO ? TRUE : FALSE) && ret;
    }

    return ret;
}


static inline void rvec_MakeZero( rvec v )
{
    int i;

#if defined(_OPENMP)
    #pragma omp simd
#endif
    for ( i = 0; i < 3; ++i )
    {
        v[i] = 0.0;
    }
}


/* Note: currently not thread-safe since uses random( ) underneath,
 * change to reentrant version later if needed */
static inline void rvec_Random( rvec v )
{
    v[0] = Random(2.0) - 1.0;
    v[1] = Random(2.0) - 1.0;
    v[2] = Random(2.0) - 1.0;
}


static inline void rtensor_Multiply( rtensor ret, rtensor m1, rtensor m2 )
{
    unsigned int i, j, k;
    rtensor temp;

    // check if the result matrix is the same as one of m1, m2.
    // if so, we cannot modify the contents of m1 or m2, so
    // we have to use a temp matrix.
    if ( ret == m1 || ret == m2 )
    {
        for ( i = 0; i < 3; ++i )
        {
            for ( j = 0; j < 3; ++j )
            {
                temp[i][j] = 0;
                for ( k = 0; k < 3; ++k )
                {
                    temp[i][j] += m1[i][k] * m2[k][j];
                }
            }
        }

        for ( i = 0; i < 3; ++i )
        {
            for ( j = 0; j < 3; ++j )
            {
                ret[i][j] = temp[i][j];
            }
        }
    }
    else
    {
        for ( i = 0; i < 3; ++i )
        {
            for ( j = 0; j < 3; ++j )
            {
                ret[i][j] = 0;
                for ( k = 0; k < 3; ++k )
                {
                    ret[i][j] += m1[i][k] * m2[k][j];
                }
            }
        }
    }
}


static inline void rtensor_MatVec( rvec ret, rtensor m, const rvec v )
{
    unsigned int i;
    rvec temp;

    // if ret is the same vector as v, we cannot modify the
    // contents of v until all computation is finished.
    if ( ret == v )
    {
        for ( i = 0; i < 3; ++i )
        {
            temp[i] = m[i][0] * v[0] + m[i][1] * v[1] + m[i][2] * v[2];
        }

        for ( i = 0; i < 3; ++i )
        {
            ret[i] = temp[i];
        }
    }
    else
    {
        for ( i = 0; i < 3; ++i )
        {
            ret[i] = m[i][0] * v[0] + m[i][1] * v[1] + m[i][2] * v[2];
        }
    }
}


static inline void rtensor_Scale( rtensor ret, const real c, rtensor m )
{
    unsigned int i, j;

    for ( i = 0; i < 3; ++i )
    {
        for ( j = 0; j < 3; ++j )
        {
            ret[i][j] = c * m[i][j];
        }
    }
}


static inline void rtensor_Add( rtensor ret, rtensor t )
{
    int i, j;

    for ( i = 0; i < 3; ++i )
    {
        for ( j = 0; j < 3; ++j )
        {
            ret[i][j] += t[i][j];
        }
    }
}


static inline void rtensor_ScaledAdd( rtensor ret, const real c, rtensor t )
{
    unsigned int i, j;

    for ( i = 0; i < 3; ++i )
    {
        for ( j = 0; j < 3; ++j )
        {
            ret[i][j] += c * t[i][j];
        }
    }
}


static inline void rtensor_Sum( rtensor ret, rtensor t1, rtensor t2 )
{
    unsigned int i, j;

    for ( i = 0; i < 3; ++i )
    {
        for ( j = 0; j < 3; ++j )
        {
            ret[i][j] = t1[i][j] + t2[i][j];
        }
    }
}


static inline void rtensor_ScaledSum( rtensor ret, const real c1, rtensor t1,
        const real c2, rtensor t2 )
{
    unsigned int i, j;

    for ( i = 0; i < 3; ++i )
    {
        for ( j = 0; j < 3; ++j )
        {
            ret[i][j] = c1 * t1[i][j] + c2 * t2[i][j];
        }
    }
}


static inline void rtensor_Copy( rtensor ret, rtensor t )
{
    unsigned int i, j;

    for ( i = 0; i < 3; ++i )
    {
        for ( j = 0; j < 3; ++j )
        {
            ret[i][j] = t[i][j];
        }
    }
}


static inline void rtensor_Identity( rtensor t )
{
    unsigned int i, j;

    for ( i = 0; i < 3; ++i )
    {
        for ( j = 0; j < 3; ++j )
        {
            t[i][j] = (i == j ? 1.0 : 0.0);
        }
    }
}


static inline void rtensor_MakeZero( rtensor t )
{
    unsigned int i, j;

    for ( i = 0; i < 3; ++i )
    {
        for ( j = 0; j < 3; ++j )
        {
            t[i][j] = 0.0;
        }
    }
}


static inline void rtensor_Transpose( rtensor ret, rtensor t )
{
    unsigned int i, j;

    for ( i = 0; i < 3; ++i )
    {
        for ( j = 0; j < 3; ++j )
        {
            ret[i][j] = t[j][i];
        }
    }
}


static inline real rtensor_Det( rtensor t )
{
    return t[0][0] * (t[1][1] * t[2][2] - t[1][2] * t[2][1])
        + t[0][1] * (t[1][2] * t[2][0] - t[1][0] * t[2][2])
        + t[0][2] * (t[1][0] * t[2][1] - t[1][1] * t[2][0]);
}


static inline real rtensor_Trace( rtensor t )
{
    return t[0][0] + t[1][1] + t[2][2];
}


static inline void ivec_MakeZero( ivec v )
{
    int i;

#if defined(_OPENMP)
    #pragma omp simd
#endif
    for ( i = 0; i < 3; ++i )
    {
        v[i] = 0;
    }
}


static inline void ivec_Copy( ivec dest, const ivec src )
{
    int i;

#if defined(_OPENMP)
    #pragma omp simd
#endif
    for ( i = 0; i < 3; ++i )
    {
        dest[i] = src[i];
    }
}


static inline void ivec_Scale( ivec dest, const int C, const ivec src )
{
    int i;

#if defined(_OPENMP)
    #pragma omp simd
#endif
    for ( i = 0; i < 3; ++i )
    {
        dest[i] = C * src[i];
    }
}


static inline void ivec_Add( ivec dest, const ivec src )
{
    int i;

#if defined(_OPENMP)
    #pragma omp simd
#endif
    for ( i = 0; i < 3; ++i )
    {
        dest[i] += src[i];
    }
}


static inline void ivec_ScaledAdd( ivec dest, const int c, const ivec src )
{
    int i;

#if defined(_OPENMP)
    #pragma omp simd
#endif
    for ( i = 0; i < 3; ++i )
    {
        dest[i] += c * src[i];
    }
}


static inline void ivec_rScale( ivec dest, const real C, const rvec src )
{
    int i;

#if defined(_OPENMP)
    #pragma omp simd
#endif
    for ( i = 0; i < 3; ++i )
    {
        dest[i] = (int)(C * src[i]);
    }
}


static inline int ivec_isZero( const ivec v )
{
    int i, ret;

    ret = TRUE;

#if defined(_OPENMP)
    #pragma omp simd reduction(&&: ret)
#endif
    for ( i = 0; i < 3; ++i )
    {
        ret = (v[i] == 0 ? TRUE : FALSE) && ret;
    }

    return ret;
}


static inline int ivec_isEqual( const ivec v1, const ivec v2 )
{
    int i, ret;

    ret = TRUE;

#if defined(_OPENMP)
    #pragma omp simd reduction(&&: ret)
#endif
    for ( i = 0; i < 3; ++i )
    {
        ret = (v1[i] == v2[i] ? TRUE : FALSE) && ret;
    }

    return ret;
}


static inline void ivec_Sum( ivec dest, const ivec v1, const ivec v2 )
{
    int i;

#if defined(_OPENMP)
    #pragma omp simd
#endif
    for ( i = 0; i < 3; ++i )
    {
        dest[i] = v1[i] + v2[i];
    }
}


#endif
