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

#include "lin_alg.h"

#include "allocate.h"
#include "tool_box.h"
#include "vector.h"

/* Intel MKL */
#if defined(HAVE_LAPACKE_MKL)
  #include "mkl.h"
/* reference LAPACK */
#elif defined(HAVE_LAPACKE)
  #include "lapacke.h"
#endif

/* for DBL_EPSILON */
//#include <float.h>


typedef struct
{
    unsigned int j;
    real val;
} sparse_matrix_entry;


enum preconditioner_type
{
    LEFT = 0,
    RIGHT = 1,
};


#if defined(TEST_MAT)
static sparse_matrix * create_test_mat( void )
{
    unsigned int i, n;
    sparse_matrix *H_test;

    Allocate_Matrix( &H_test, 3, 3, 6 );

    //3x3, SPD, store lower half
    i = 0;
    n = 0;
    H_test->start[n] = i;
    H_test->j[i] = 0;
    H_test->val[i] = 4.;
    ++i;
    ++n;
    H_test->start[n] = i;
    H_test->j[i] = 0;
    H_test->val[i] = 12.;
    ++i;
    H_test->j[i] = 1;
    H_test->val[i] = 37.;
    ++i;
    ++n;
    H_test->start[n] = i;
    H_test->j[i] = 0;
    H_test->val[i] = -16.;
    ++i;
    H_test->j[i] = 1;
    H_test->val[i] = -43.;
    ++i;
    H_test->j[i] = 2;
    H_test->val[i] = 98.;
    ++i;
    ++n;
    H_test->start[n] = i;

    return H_test;
}
#endif


/* Routine used with qsort for sorting nonzeros within a sparse matrix row
 *
 * v1/v2: pointers to column indices of nonzeros within a row (unsigned int)
 */
static int compare_matrix_entry(const void *v1, const void *v2)
{
    /* larger element has larger column index */
    return ((sparse_matrix_entry *)v1)->j - ((sparse_matrix_entry *)v2)->j;
}


/* Routine used for sorting nonzeros within a sparse matrix row;
 *  internally, a combination of qsort and manual sorting is utilized
 *  (parallel calls to qsort when multithreading, rows mapped to threads)
 *
 * A: sparse matrix for which to sort nonzeros within a row, stored in CSR format
 */
void Sort_Matrix_Rows( sparse_matrix * const A )
{
    unsigned int i, j, si, ei;
    size_t temp_size;
    sparse_matrix_entry *temp;

#if defined(_OPENMP)
//    #pragma omp parallel default(none) private(i, j, si, ei, temp) shared(stderr)
#endif
    {
        temp = NULL;
        temp_size = 0;

        /* sort each row of A using column indices */
#if defined(_OPENMP)
//        #pragma omp for schedule(guided)
#endif
        for ( i = 0; i < A->n; ++i )
        {
            si = A->start[i];
            ei = A->start[i + 1];

            if ( temp_size < ei - si )
            {
                if ( temp != NULL )
                {
                    sfree( temp, __FILE__, __LINE__ );
                }
                temp = smalloc( sizeof(sparse_matrix_entry) * (ei - si), __FILE__, __LINE__ );
                temp_size = ei - si;
            }

            for ( j = 0; j < (ei - si); ++j )
            {
                temp[j].j = A->j[si + j];
                temp[j].val = A->val[si + j];
            }

            /* polymorphic sort in standard C library using column indices */
            qsort( temp, ei - si, sizeof(sparse_matrix_entry), compare_matrix_entry );

            for ( j = 0; j < (ei - si); ++j )
            {
                A->j[si + j] = temp[j].j;
                A->val[si + j] = temp[j].val;
            }
        }

        sfree( temp, __FILE__, __LINE__ );
    }
}


/* Convert a symmetric, half-sored sparse matrix into
 * a full-stored sparse matrix
 *
 * A: symmetric sparse matrix, lower half stored in CSR
 * A_full: resultant full sparse matrix in CSR
 *   If A_full is NULL, allocate space, otherwise do not
 *
 * Assumptions:
 *   A has non-zero diagonals
 *   Each row of A has at least one non-zero (i.e., no rows with all zeros) */
static void compute_full_sparse_matrix( const sparse_matrix * const A,
        sparse_matrix * A_full, int realloc )
{
    int count, i, pj;
    sparse_matrix A_t;

    if ( A_full->allocated == FALSE )
    {
        Allocate_Matrix( A_full, A->n, A->n_max, 2 * A->m - A->n );
    }
    else if ( A_full->m < 2 * A->m - A->n || A_full->n_max < A->n_max
            || realloc == TRUE )
    {
        if ( A_full->allocated == TRUE )
        {
            Deallocate_Matrix( A_full );
        }
        Allocate_Matrix( A_full, A->n, A->n_max, 2 * A->m - A->n );
    }

    Allocate_Matrix( &A_t, A->n, A->n_max, A->m );

    /* Set up the sparse matrix data structure for A. */
    Transpose( A, &A_t );

    count = 0;
    for ( i = 0; i < A->n; ++i )
    {
        A_full->start[i] = count;

        /* A: symmetric, lower triangular portion only stored */
        for ( pj = A->start[i]; pj < A->start[i + 1]; ++pj )
        {
            A_full->val[count] = A->val[pj];
            A_full->j[count] = A->j[pj];
            ++count;
        }
        /* A^T: symmetric, upper triangular portion only stored;
         * skip diagonal from A^T, as included from A above */
        for ( pj = A_t.start[i] + 1; pj < A_t.start[i + 1]; ++pj )
        {
            A_full->val[count] = A_t.val[pj];
            A_full->j[count] = A_t.j[pj];
            ++count;
        }
    }
    A_full->start[A->n] = count;

    Deallocate_Matrix( &A_t );
}


/* Setup routines for sparse approximate inverse preconditioner
 *
 * A: symmetric sparse matrix, lower half stored in CSR
 * A_full:
 * A_spar_patt:
 * A_spar_patt_full:
 * A_app_inv:
 * filter:
 *
 * Assumptions:
 *   A has non-zero diagonals
 *   Each row of A has at least one non-zero (i.e., no rows with all zeros) */
void setup_sparse_approx_inverse( const sparse_matrix * const A,
        sparse_matrix * A_full, sparse_matrix * A_spar_patt,
        sparse_matrix * A_spar_patt_full, sparse_matrix * A_app_inv,
        const real filter, int realloc )
{
    int i, pj, size;
    int left, right, k, p, turn;
    real pivot, tmp;
    real threshold;
    real *list;

    if ( A_spar_patt->allocated == FALSE )
    {
        Allocate_Matrix( A_spar_patt, A->n, A->n_max, A->m );
    }
    else if ( A_spar_patt->m < A->m || A_spar_patt->n_max < A->n_max
            || realloc == TRUE )
    {
        if ( A_spar_patt->allocated == TRUE )
        {
            Deallocate_Matrix( A_spar_patt );
        }
        Allocate_Matrix( A_spar_patt, A->n, A->n_max, A->m );
    }

    list = smalloc( sizeof(real) * A->start[A->n], __FILE__, __LINE__ );

    /* quick-select algorithm for finding the k-th greatest element in the matrix, where
     *  list: values from the matrix
     *  left, right: search space of the quick-select */
    left = 0;
    right = A->start[A->n] - 1;
    k = (int)( A->start[A->n] * filter );
    threshold = 0.0;

    for ( i = left; i <= right; ++i )
    {
        list[i] = A->val[i];

        if ( list[i] < 0.0 )
        {
            list[i] = -list[i];
        }
    }

    turn = 0;
    while ( k > 0 )
    {
        p = left;
        turn = 1 - turn;

        if ( turn == 1 )
        {
            pivot = list[right];
        }
        else
        {
            pivot = list[left];
        }

        for ( i = left + 1 - turn; i <= right - turn; ++i )
        {
            if ( list[i] > pivot )
            {
                tmp = list[i];
                list[i] = list[p];
                list[p] = tmp;
                p++;
            }
        }

        if ( turn == 1 )
        {
            tmp = list[p];
            list[p] = list[right];
            list[right] = tmp;
        }
        else
        {
            tmp = list[p];
            list[p] = list[left];
            list[left] = tmp;
        }

        if ( p == k - 1 )
        {
            threshold = list[p];
            break;
        }
        else if ( p > k - 1 )
        {
            right = p - 1;
        }
        else
        {
            left = p + 1;
        }
    }

    /* special case for EE/ACKS2 where chosen threshold is one of the
     * many matrix entries with value 1.0 => exclude all 1.0 values */
    if ( FABS( threshold - 1.0 ) < 1.0e-10 )
    {
        threshold = 1.000001;
    }

    /* fill sparsity pattern */
    for ( size = 0, i = 0; i < A->n; ++i )
    {
        A_spar_patt->start[i] = size;

        for ( pj = A->start[i]; pj < A->start[i + 1] - 1; ++pj )
        {
            if ( FABS( A->val[pj] ) >= threshold )
            {
                A_spar_patt->val[size] = A->val[pj];
                A_spar_patt->j[size] = A->j[pj];
                size++;
            }
        }

        /* diagonal entries are always included */
        A_spar_patt->val[size] = A->val[pj];
        A_spar_patt->j[size] = A->j[pj];
        size++;
    }
    A_spar_patt->start[A->n] = size;

    compute_full_sparse_matrix( A, A_full, realloc );
    compute_full_sparse_matrix( A_spar_patt, A_spar_patt_full, realloc );

    if ( A_app_inv->allocated == FALSE )
    {
        /* A_app_inv has the same sparsity pattern
         * as A_spar_patt_full (omit non-zero values) */
        Allocate_Matrix( A_app_inv, A_spar_patt_full->n,
                A_spar_patt_full->n_max, A_spar_patt_full->m );
    }
    else if ( A_app_inv->m < A->m || A_app_inv->n_max < A->n_max
            || realloc == TRUE )
    {
        if ( A_app_inv->allocated == TRUE )
        {
            Deallocate_Matrix( A_app_inv );
        }

        /* A_app_inv has the same sparsity pattern
         * as A_spar_patt_full (omit non-zero values) */
        Allocate_Matrix( A_app_inv, A_spar_patt_full->n,
                A_spar_patt_full->n_max, A_spar_patt_full->m );
    }

    sfree( list, __FILE__, __LINE__ );
}


/* Computes the 2-norm of each row of A, multiplied by the
 * dropping tolerance dtol 
 * 
 * A: symmetric (lower triangular portion only stored), square matrix, CSR format
 * droptol (output): row-wise dropping tolereances
 * dtol: user-specified dropping tolerance */
void Calculate_Droptol( const sparse_matrix * const A,
        real * const droptol, const real dtol )
{
    int i, j, k;
    real val;
#if defined(_OPENMP)
    static real *droptol_local;
    unsigned int tid;
#endif

#if defined(_OPENMP)
    #pragma omp parallel default(none) private(i, j, k, val, tid) \
    firstprivate(A, droptol, dtol) shared(droptol_local, stderr)
#endif
    {
#if defined(_OPENMP)
        tid = omp_get_thread_num();

        #pragma omp master
        {
            if ( droptol_local == NULL )
            {
                droptol_local = smalloc( omp_get_num_threads() * A->n * sizeof(real),
                        __FILE__, __LINE__ );
            }
        }

        #pragma omp barrier
#endif

        for ( i = 0; i < A->n; ++i )
        {
#if defined(_OPENMP)
            droptol_local[tid * A->n + i] = 0.0;
#else
            droptol[i] = 0.0;
#endif
        }

#if defined(_OPENMP)
        #pragma omp barrier
#endif

        /* calculate sqaure of the norm of each row */
#if defined(_OPENMP)
        #pragma omp for schedule(static)
#endif
        for ( i = 0; i < A->n; ++i )
        {
            for ( k = A->start[i]; k < A->start[i + 1] - 1; ++k )
            {
                j = A->j[k];
                val = A->val[k];

#if defined(_OPENMP)
                droptol_local[tid * A->n + i] += val * val;
                droptol_local[tid * A->n + j] += val * val;
#else
                droptol[i] += val * val;
                droptol[j] += val * val;
#endif
            }

            // diagonal entry
            val = A->val[k];
#if defined(_OPENMP)
            droptol_local[tid * A->n + i] += val * val;
#else
            droptol[i] += val * val;
#endif
        }

#if defined(_OPENMP)
        #pragma omp barrier

        #pragma omp for schedule(static)
        for ( i = 0; i < A->n; ++i )
        {
            droptol[i] = 0.0;
            for ( k = 0; k < omp_get_num_threads(); ++k )
            {
                droptol[i] += droptol_local[k * A->n + i];
            }
        }

        #pragma omp barrier
#endif

        /* calculate local droptol for each row */
#if defined(_OPENMP)
        #pragma omp for schedule(static)
#endif
        for ( i = 0; i < A->n; ++i )
        {
            droptol[i] = SQRT( droptol[i] ) * dtol;
        }

#if defined(_OPENMP)
        #pragma omp master
        {
            sfree( droptol_local, __FILE__, __LINE__ );
        }
#endif
    }
}


int Estimate_LU_Fill( const sparse_matrix * const A, const real * const droptol )
{
    int i, pj;
    int fillin;
    real val;

    fillin = 0;

#if defined(_OPENMP)
    #pragma omp parallel for schedule(static) \
    default(none) private(i, pj, val) firstprivate(A, droptol) reduction(+: fillin)
#endif
    for ( i = 0; i < A->n; ++i )
    {
        for ( pj = A->start[i]; pj < A->start[i + 1] - 1; ++pj )
        {
            val = A->val[pj];

            if ( FABS(val) > droptol[i] )
            {
                ++fillin;
            }
        }
    }

    return fillin + A->n;
}


/* Compute diagonal inverese (Jacobi) preconditioner
 *
 * H: matrix used to compute preconditioner, in CSR format
 * Hdia_inv: computed diagonal inverse preconditioner
 * cm: charge model selected for generating the charge matrix
 * N: number of atoms (i.e., charge matrix dimensions for QEq)
 */
real jacobi( const sparse_matrix * const H, real * const Hdia_inv,
        int cm, int N )
{
    unsigned int i;
    real start;

    start = Get_Time( );

    if ( cm == QEQ_CM || cm == EE_CM )
    {
#if defined(_OPENMP)
        #pragma omp parallel for schedule(dynamic,256) \
        default(none) private(i) firstprivate(H, Hdia_inv)
#endif
        for ( i = 0; i < H->n; ++i )
        {
            if ( FABS( H->val[H->start[i + 1] - 1] ) > 1.0e-15 )
            {
                Hdia_inv[i] = 1.0 / H->val[H->start[i + 1] - 1];
            }
            else
            {
                Hdia_inv[i] = 1.0;
            }
        }
    }
    else if ( cm == ACKS2_CM )
    {
#if defined(_OPENMP)
        #pragma omp parallel for schedule(dynamic,256) \
        default(none) private(i) firstprivate(H, Hdia_inv, N)
#endif
        for ( i = 0; i < H->n; ++i )
        {
            if ( i < N && FABS( H->val[H->start[i + 1] - 1] ) > 1.0e-15 )
            {
                Hdia_inv[i] = 1.0 / H->val[H->start[i + 1] - 1];
            }
            else
            {
                Hdia_inv[i] = 1.0;
            }
        }
    }
    else
    {
        fprintf( stderr, "[ERROR] Unrecognized charge model (%d). Terminating...\n",
                cm );
        exit( INVALID_INPUT );
    }

    return Get_Timing_Info( start );
}


/* Incomplete Cholesky factorization with dual thresholding */
real ICHOLT( const sparse_matrix * const A, const real * const droptol,
             sparse_matrix * const L, sparse_matrix * const U )
{
    int i, j, pj, k1, k2, tmptop, Ltop;
    int *tmp_j;
    real val, start;
    real *tmp_val;
    unsigned int *Utop;

    start = Get_Time( );

    Utop = smalloc( (A->n + 1) * sizeof(unsigned int), __FILE__, __LINE__ );
    tmp_j = smalloc( A->n * sizeof(int), __FILE__, __LINE__ );
    tmp_val = smalloc( A->n * sizeof(real), __FILE__, __LINE__ );

    Ltop = 0;
    tmptop = 0;
    for ( i = 0; i < A->n + 1; ++i )
    {
        L->start[i] = 0;
    }
    for ( i = 0; i < A->n + 1; ++i )
    {
        U->start[i] = 0;
    }
    for ( i = 0; i < A->n; ++i )
    {
        Utop[i] = 0;
    }

    for ( i = 0; i < A->n; ++i )
    {
        L->start[i] = Ltop;
        tmptop = 0;

        for ( pj = A->start[i]; pj < A->start[i + 1] - 1; ++pj )
        {
            j = A->j[pj];
            val = A->val[pj];

            if ( FABS(val) > droptol[i] )
            {
                k1 = 0;
                k2 = L->start[j];
                while ( k1 < tmptop && k2 < L->start[j + 1] )
                {
                    if ( tmp_j[k1] < L->j[k2] )
                    {
                        ++k1;
                    }
                    else if ( tmp_j[k1] > L->j[k2] )
                    {
                        ++k2;
                    }
                    else
                    {
                        val -= (tmp_val[k1++] * L->val[k2++]);
                    }
                }

                // L matrix is lower triangular,
                // so right before the start of next row comes jth diagonal
                val /= L->val[L->start[j + 1] - 1];

                tmp_j[tmptop] = j;
                tmp_val[tmptop] = val;
                ++tmptop;
            }
        }

        // sanity check
        if ( A->j[pj] != i )
        {
            fprintf( stderr, "[ERROR] ICHOLT: badly built A matrix!\n (i = %d) ", i );
            exit( NUMERIC_BREAKDOWN );
        }

        // compute the ith diagonal in L
        val = A->val[pj];
        for ( k1 = 0; k1 < tmptop; ++k1 )
        {
            val -= (tmp_val[k1] * tmp_val[k1]);
        }

#if defined(DEBUG)
        if ( val < 0.0 )
        {
            fprintf( stderr, "[INFO] ICHOLT: numeric breakdown (SQRT of negative on diagonal i = %d). Terminating.\n", i );
            exit( NUMERIC_BREAKDOWN );

        }
#endif

        tmp_j[tmptop] = i;
        tmp_val[tmptop] = SQRT( val );

        // apply the dropping rule once again
        //fprintf( stderr, "row%d: tmptop: %d\n", i, tmptop );
        //for( k1 = 0; k1<= tmptop; ++k1 )
        //  fprintf( stderr, "%d(%f)  ", tmp[k1].j, tmp[k1].val );
        //fprintf( stderr, "\n" );
        //fprintf( stderr, "row(%d): droptol=%.4f\n", i+1, droptol[i] );
        for ( k1 = 0; k1 < tmptop; ++k1 )
        {
            if ( FABS(tmp_val[k1]) > droptol[i] / tmp_val[tmptop] )
            {
                L->j[Ltop] = tmp_j[k1];
                L->val[Ltop] = tmp_val[k1];
                U->start[tmp_j[k1] + 1]++;
                ++Ltop;
                //fprintf( stderr, "%d(%.4f)  ", tmp[k1].j+1, tmp[k1].val );
            }
        }
        // keep the diagonal in any case
        L->j[Ltop] = tmp_j[k1];
        L->val[Ltop] = tmp_val[k1];
        ++Ltop;
        //fprintf( stderr, "%d(%.4f)\n", tmp[k1].j+1,  tmp[k1].val );
    }

    L->start[i] = Ltop;
    //    fprintf( stderr, "nnz(L): %d, max: %d\n", Ltop, L->n * 50 );

    /* U = L^T (Cholesky factorization) */
    Transpose( L, U );
    //    for ( i = 1; i <= U->n; ++i )
    //    {
    //        Utop[i] = U->start[i] = U->start[i] + U->start[i - 1] + 1;
    //    }
    //    for ( i = 0; i < L->n; ++i )
    //    {
    //        for ( pj = L->start[i]; pj < L->start[i + 1]; ++pj )
    //        {
    //            j = L->j[pj];
    //            U->j[Utop[j]] = i;
    //            U->val[Utop[j]] = L->val[pj];
    //            Utop[j]++;
    //        }
    //    }

    //    fprintf( stderr, "nnz(U): %d, max: %d\n", Utop[U->n], U->n * 50 );

    sfree( tmp_val, __FILE__, __LINE__ );
    sfree( tmp_j, __FILE__, __LINE__ );
    sfree( Utop, __FILE__, __LINE__ );

    return Get_Timing_Info( start );
}


/* Compute incomplete LU factorization with 0-fillin and no pivoting
 *
 * Reference:
 * Iterative Methods for Sparse Linear System, Second Edition, 2003,
 * Yousef Saad
 *
 * A: symmetric (lower triangular portion only stored), square matrix, CSR format
 * L / U: (output) triangular matrices, A \approx LU, CSR format */
real ILU( const sparse_matrix * const A, sparse_matrix * const L,
        sparse_matrix * const U )
{
    unsigned int i, pj, pk, pl, Ltop, Utop;
    real start;
    sparse_matrix *A_full = NULL;

    start = Get_Time( );

    compute_full_sparse_matrix( A, A_full, FALSE );

    Ltop = 0;
    Utop = 0;

    for ( i = 0; i < A_full->n; ++i )
    {
        /* mark the start of new row i in L and U */
        L->start[i] = Ltop;
        U->start[i] = Utop;

        /* for each non-zero in i-th row to the left of the diagonal:
         * for k = 0, ..., i - 1 */
        for ( pj = A_full->start[i]; pj < A_full->start[i + 1]; ++pj )
        {
            if ( A_full->j[pj] >= i )
            {
                break;
            }

            /* scan k-th row (A_full->j[pj]) to find a_{kk},
             * and compute a_{ik} = a_{ik} / a_{kk} */
            for ( pk = A_full->start[A_full->j[pj]]; pk < A_full->start[A_full->j[pj] + 1]; ++pk )
            {
                if ( A_full->j[pk] == A_full->j[pj] )
                {
                    A_full->val[pj] /= A_full->val[pk];
                    break;
                }
            }

            /* trailing row update (sparse vector-sparse vector product):
             * for j = k + 1, ..., n - 1 
             *   a_{ij} = a_{ij} - a_{ik} * a_{kj}
             *
             * pi: points to a_{ik}
             * pl: points to a_{ij}
             * pk: points to a_{kj}
             * */
            ++pk;
            for ( pl = pj + 1; pl < A_full->start[i + 1] && pk < A_full->start[A_full->j[pj] + 1]; )
            {
                if ( A_full->j[pl] == A_full->j[pk] )
                {
                    A_full->val[pl] -= A_full->val[pj] * A_full->val[pk];
                    ++pl;
                    ++pk;
                }
                else if ( A_full->j[pl] < A_full->j[pk] )
                {
                    ++pl;
                }
                else
                {
                    ++pk;
                }
            }
        }

        /* copy A_full[0:i-1] to row i of L */
        for ( pj = A_full->start[i]; pj < A_full->start[i + 1]; ++pj )
        {
            if ( A_full->j[pj] >= i )
            {
                break;
            }

            L->j[Ltop] = A_full->j[pj];
            L->val[Ltop] = A_full->val[pj];
            ++Ltop;
        }

        /* unit diagonal for row i of L */
        L->j[Ltop] = i;
        L->val[Ltop] = 1.0;
        ++Ltop;

        /* copy A_full[i:n-1] to row i of U */
        for ( ; pj < A_full->start[i + 1]; ++pj )
        {
            U->j[Utop] = A_full->j[pj];
            U->val[Utop] = A_full->val[pj];
            ++Utop;
        }
    }

    /* record the total NNZ as the last entry in row pointer (start) */
    L->start[L->n] = Ltop;
    U->start[U->n] = Utop;

    Deallocate_Matrix( A_full );

    return Get_Timing_Info( start );
}


/* Compute incomplete LU factorization with thresholding
 *
 * Reference:
 * Iterative Methods for Sparse Linear System, Second Edition, 2003,
 * Yousef Saad
 *
 * A: symmetric (lower triangular portion only stored), square matrix, CSR format
 * droptol: row-wise tolerances used for dropping
 * L / U: (output) triangular matrices, A \approx LU, CSR format */
real ILUT( const sparse_matrix * const A, const real * const droptol,
        sparse_matrix * const L, sparse_matrix * const U )
{
    unsigned int i, k, pj, Ltop, Utop, *nz_mask;
    real *w, start;
    sparse_matrix *A_full = NULL;
    unsigned int nz_cnt;

    start = Get_Time( );

    /* use a dense vector with masking for the intermediate row w */
    w = smalloc( sizeof(real) * A->n, __FILE__, __LINE__ );
    nz_mask = smalloc( sizeof(unsigned int) * A->n, __FILE__, __LINE__ );

    compute_full_sparse_matrix( A, A_full, FALSE );

    Ltop = 0;
    Utop = 0;
    for ( i = 0; i < L->n + 1; ++i )
        L->start[i] = 0;
    for ( i = 0; i < U->n + 1; ++i )
        U->start[i] = 0;

    for ( i = 0; i < A_full->n; ++i )
    {
        /* mark the start of new row i in L and U */
        L->start[i] = Ltop;
        U->start[i] = Utop;

        for ( k = 0; k < A_full->n; ++k )
        {
            nz_mask[k] = 0;
            w[k] = 0.0;
        }

        /* copy i-th row of A_full into w */
        for ( pj = A_full->start[i]; pj < A_full->start[i + 1]; ++pj )
        {
            nz_mask[k] = 1;
            k = A_full->j[pj];
            w[k] = A_full->val[pj];
        }

        for ( k = 0; k < i; ++k )
        {
            if ( nz_mask[k] == 1 )
            {
                /* A symmetric, lower triangular portion stored */
                w[k] /= A->val[A->start[k + 1] - 1];

                /* apply dropping rule to w[k] */
                if ( FABS( w[k] ) <= droptol[k] )
                {
                    nz_mask[k] = 0;
                }

                /* subtract scaled k-th row of U from w */
                if ( nz_mask[k] == 1 )
                {
                    for ( pj = U->start[k]; pj < U->start[k + 1]; ++pj )
                    {
                        nz_mask[U->j[pj]] = 1;
                        w[U->j[pj]] -= w[k] * U->val[pj];
                    }
                }
            }
        }

        /* apply dropping rule to w, but keep the diagonal regardless;
         * note: this is different than Saad's suggested approach
         * as we do not limit the NNZ per row */
        for ( k = 0; k < A_full->n; ++k )
        {
            if ( FABS( w[k] ) <= droptol[i] )
            {
                nz_mask[k] = 0;
            }
        }

        nz_cnt = 0;
        for ( k = i + 1; k < A_full->n; ++k )
            if ( nz_mask[k] == 1 )
                ++nz_cnt;

        if ( Ltop + nz_cnt > L->m )
        {
            L->m = MAX( (5 * nz_cnt) + L->m, (unsigned int) (L->m * SAFE_ZONE) );
            L->j = srealloc( L->j, sizeof(unsigned int) * L->m, __FILE__, __LINE__ );
            L->val = srealloc( L->val, sizeof(real) * L->m, __FILE__, __LINE__ );
        }

        /* copy w[0:i-1] to row i of L */
        for ( k = 0; k < i; ++k )
        {
            if ( nz_mask[k] == 1 )
            {
                L->j[Ltop] = k;
                L->val[Ltop] = w[k];
                ++Ltop;
            }
        }

        /* unit diagonal for L */
        L->j[Ltop] = i;
        L->val[Ltop] = 1.0;
        ++Ltop;

        nz_cnt = 0;
        for ( k = i + 1; k < A_full->n; ++k )
            if ( nz_mask[k] == 1 )
                ++nz_cnt;

        if ( Utop + nz_cnt > U->m )
        {
            U->m = MAX( (5 * nz_cnt) + U->m, (unsigned int) (U->m * SAFE_ZONE) );
            U->j = srealloc( U->j, sizeof(unsigned int) * U->m, __FILE__, __LINE__ );
            U->val = srealloc( U->val, sizeof(real) * U->m, __FILE__, __LINE__ );
        }

        /* diagonal for U */
        U->j[Utop] = i;
        U->val[Utop] = w[i];
        ++Utop;

        /* copy w[i+1:n-1] to row i of U */
        for ( k = i + 1; k < A_full->n; ++k )
        {
            if ( nz_mask[k] == 1 )
            {
                U->j[Utop] = k;
                U->val[Utop] = w[k];
                ++Utop;
            }
        }
    }

    /* record the total NNZ as the last entry in row pointer (start) */
    L->start[L->n] = Ltop;
    U->start[U->n] = Utop;

    Deallocate_Matrix( A_full );
    sfree( nz_mask, __FILE__, __LINE__ );
    sfree( w, __FILE__, __LINE__ );

    return Get_Timing_Info( start );
}


/* Compute incomplete LU factorization with thresholding and column-based partial pivoting
 *
 * Reference:
 * Iterative Methods for Sparse Linear System, Second Edition, 2003,
 * Yousef Saad
 *
 * A: symmetric (lower triangular portion only stored), square matrix, CSR format
 * droptol: row-wise tolerances used for dropping
 * L / U: (output) triangular matrices, A \approx LU, CSR format */
real ILUTP( const sparse_matrix * const A, const real * const droptol,
        sparse_matrix * const L, sparse_matrix * const U )
{
    unsigned int i, k, pj, Ltop, Utop, *nz_mask, *perm, *perm_inv, pivot_j;
    real *w, start, pivot_val;
    sparse_matrix *A_full = NULL;

    start = Get_Time( );

    /* use a dense vector with masking for the intermediate row w */
    w = smalloc( sizeof(real) * A->n, __FILE__, __LINE__ );
    nz_mask = smalloc( sizeof(int) * A->n, __FILE__, __LINE__ );
    perm = smalloc( sizeof(int) * A->n, __FILE__, __LINE__ );
    perm_inv = smalloc( sizeof(int) * A->n, __FILE__, __LINE__ );

    compute_full_sparse_matrix( A, A_full, FALSE );

    Ltop = 0;
    Utop = 0;

    for ( i = 0; i < A->n; ++i )
        perm[i] = i;
    for ( i = 0; i < A->n; ++i )
        perm_inv[perm[i]] = i;

    for ( i = 0; i < A_full->n; ++i )
    {
        /* mark the start of new row i in L and U */
        L->start[i] = Ltop;
        U->start[i] = Utop;

        for ( k = 0; k < A_full->n; ++k )
        {
            nz_mask[k] = 0;
            w[k] = 0.0;
        }

        /* copy i-th row of A_full into w */
        for ( pj = A_full->start[i]; pj < A_full->start[i + 1]; ++pj )
        {
            nz_mask[k] = 1;
            k = A_full->j[pj];
            w[k] = A_full->val[pj];
        }

        /* partial pivoting by columns:
         * find largest element in w, and make it the pivot */
        pivot_val = FABS( w[0] );
        pivot_j = 0;
        for ( k = 1; k < A_full->n; ++k )
        {
            if ( FABS( w[k] ) > pivot_val )
            {
                pivot_val = FABS( w[k] );
                pivot_j = k;
            }
        }
        perm[i] = pivot_j;
        perm_inv[perm[i]] = i;
        perm[pivot_j] = i;
        perm_inv[perm[pivot_j]] = pivot_j;

        for ( k = 0; k < i; ++k )
        {
            if ( nz_mask[perm[k]] == 1 )
            {
                for ( pj = A->start[k]; pj < A->start[k + 1]; ++pj )
                {
                    if ( A->j[pj] == perm_inv[k] )
                    {
                        w[perm[k]] /= A->val[pj];
                        break;
                    }
                }

                /* apply dropping rule to w[perm[k]] */
                if ( FABS( w[perm[k]] ) < droptol[perm[k]] )
                {
                    nz_mask[perm[k]] = 0;
                }

                /* subtract scaled k-th row of U from w */
                if ( nz_mask[perm[k]] == 1 )
                {
                    for ( pj = U->start[k]; pj < U->start[k + 1]; ++pj )
                    {
                        nz_mask[U->j[pj]] = 1;
                        w[perm[U->j[pj]]] -= w[perm[k]] * U->val[pj];
                    }
                }
            }
        }

        /* apply dropping rule to w, but keep the diagonal regardless;
         * note: this is different than Saad's suggested approach
         * as we do not limit the NNZ per row */
        for ( k = 0; k < A_full->n; ++k )
        {
            if ( perm_inv[k] != i && nz_mask[k] == 1 && FABS( w[k] ) < droptol[i] )
            {
                nz_mask[k] = 0;
            }
        }

        /* copy w[0:i-1] to row i of L */
        for ( k = 0; k < i; ++k )
        {
            if ( nz_mask[perm_inv[k]] == 1 )
            {
                L->j[Ltop] = k;
                L->val[Ltop] = w[perm_inv[k]];
                ++Ltop;
            }
        }

        /* unit diagonal for L */
        L->j[Ltop] = i;
        L->val[Ltop] = 1.0;
        ++Ltop;

        /* diagonal for U */
        U->j[Utop] = i;
        U->val[Utop] = w[perm_inv[i]];
        ++Utop;

        /* copy w[i-1:n] to row i of U */
        for ( k = i + 1; k < A_full->n; ++k )
        {
            if ( nz_mask[perm_inv[k]] == 1 )
            {
                U->j[Utop] = k;
                U->val[Utop] = w[perm_inv[k]];
                ++Utop;
            }
        }
    }

    /* record the total NNZ as the last entry in row pointer (start) */
    L->start[L->n] = Ltop;
    U->start[U->n] = Utop;

    Deallocate_Matrix( A_full );
    sfree( perm_inv, __FILE__, __LINE__ );
    sfree( perm, __FILE__, __LINE__ );
    sfree( nz_mask, __FILE__, __LINE__ );
    sfree( w, __FILE__, __LINE__ );

    return Get_Timing_Info( start );
}


/* Fine-grained (parallel) incomplete Cholesky factorization with thresholding
 *
 * Reference:
 * Edmond Chow and Aftab Patel
 * Fine-Grained Parallel Incomplete LU Factorization
 * SIAM J. Sci. Comp.
 *
 * A: symmetric, half-stored (lower triangular), CSR format
 * droptol: row-wise tolerances used for dropping
 * sweeps: number of loops over non-zeros for computation
 * U_T / U: factorized triangular matrices (A \approx U^{T}U), CSR format */
real FG_ICHOLT( const sparse_matrix * const A, const real * droptol,
        const unsigned int sweeps, sparse_matrix * const U_T, sparse_matrix * const U )
{
    unsigned int i, pj, x = 0, y = 0, ei_x, ei_y, Utop, s;
    real *D, *D_inv, *gamma, sum, start;
    sparse_matrix DAD, U_T_temp;

    start = Get_Time( );

    Allocate_Matrix( &DAD, A->n, A->n_max, A->m );
    Allocate_Matrix( &U_T_temp, A->n, A->n_max, A->m );

    D = smalloc( sizeof(real) * A->n, __FILE__, __LINE__ );
    D_inv = smalloc( sizeof(real) * A->n, __FILE__, __LINE__ );
    gamma = smalloc( sizeof(real) * A->n, __FILE__, __LINE__ );

#if defined(_OPENMP)
    #pragma omp parallel for schedule(dynamic,512) \
        default(none) shared(D, D_inv, gamma) private(i) firstprivate(A)
#endif
    for ( i = 0; i < A->n; ++i )
    {
        if ( A->val[A->start[i + 1] - 1] < 0.0 )
        {
            gamma[i] = -1.0;
        }
        else
        {
            gamma[i] = 1.0;
        }

        D_inv[i] = SQRT( FABS( A->val[A->start[i + 1] - 1] ) );
    }

#if defined(_OPENMP)
    #pragma omp parallel for schedule(dynamic,512) \
        default(none) shared(D, D_inv, stderr) private(i) firstprivate(A)
#endif
    for ( i = 0; i < A->n; ++i )
    {
        D[i] = 1.0 / D_inv[i];
    }

    /* to get convergence, A must have unit diagonal, so apply
     * transformation DAD, where D = D(1./SQRT(D(A))) */
    memcpy( DAD.start, A->start, sizeof(unsigned int) * (A->n + 1) );
    memcpy( DAD.j, A->j, sizeof(unsigned int) * A->start[A->n] );
#if defined(_OPENMP)
    #pragma omp parallel for schedule(dynamic,4096) \
        default(none) shared(DAD, D, gamma) private(i, pj) firstprivate(A)
#endif
    for ( i = 0; i < A->n; ++i )
    {
        /* non-diagonals */
        for ( pj = A->start[i]; pj < A->start[i + 1] - 1; ++pj )
        {
            DAD.val[pj] = gamma[i] * (D[i] * A->val[pj] * D[A->j[pj]]);
        }
        /* diagonal */
        DAD.val[pj] = 1.0;
    }

    /* initial guesses for U^T,
     * assume: A and DAD symmetric and stored lower triangular */
    memcpy( U_T_temp.start, DAD.start, sizeof(unsigned int) * (DAD.n + 1) );
    memcpy( U_T_temp.j, DAD.j, sizeof(unsigned int) * DAD.start[DAD.n] );
    memcpy( U_T_temp.val, DAD.val, sizeof(real) * DAD.start[DAD.n] );

    for ( s = 0; s < sweeps; ++s )
    {
        /* for each nonzero in U^{T} */
#if defined(_OPENMP)
        #pragma omp parallel for schedule(dynamic,4096) \
            default(none) shared(DAD, U_T_temp, stderr) private(i, pj, x, y, ei_x, ei_y, sum)
#endif
        for ( pj = 0; pj < U_T_temp.start[U_T_temp.n]; ++pj )
        {
            /* determine row bounds of current nonzero */
            x = 0;
            ei_x = pj;
            for ( i = 1; i <= U_T_temp.n; ++i )
            {
                if ( U_T_temp.start[i] > pj )
                {
                    --i;
                    x = U_T_temp.start[i];
                    break;
                }
            }
            /* column bounds of current nonzero */
            y = U_T_temp.start[U_T_temp.j[pj]];
            ei_y = U_T_temp.start[U_T_temp.j[pj] + 1];

            sum = 0.0;

            /* sparse vector-sparse vector inner product for nonzero (i, j):
             *   dot( U^T(i,1:j-1), U^T(j,1:j-1) ) */
            for ( ; x < ei_x && y < ei_y && U_T_temp.j[y] < U_T_temp.j[pj]; )
            {
                if ( U_T_temp.j[x] == U_T_temp.j[y] )
                {
                    sum += U_T_temp.val[x] * U_T_temp.val[y];
                    ++x;
                    ++y;
                }
                else if ( U_T_temp.j[x] < U_T_temp.j[y] )
                {
                    ++x;
                }
                else
                {
                    ++y;
                }
            }

            sum = DAD.val[pj] - sum;

            if ( i == U_T_temp.j[pj] )
            {
#if defined(DEBUG_FOCUS)
                /* sanity check */
                if ( sum < 0.0 )
                {
                    fprintf( stderr, "[ERROR] Numeric breakdown in FG_ICHOLT. Terminating.\n");
                    fprintf( stderr, "  [INFO] DAD(%5d,%5d) = %10.3f\n",
                             i, DAD.j[pj], DAD.val[pj] );
                    fprintf( stderr, "  [INFO] sum = %10.3f\n", sum);
                    exit( NUMERIC_BREAKDOWN );
                }
#endif

                U_T_temp.val[pj] = SQRT( FABS( sum ) );
            }
            /* non-diagonal entries */
            else
            {
                U_T_temp.val[pj] = sum / U_T_temp.val[ei_y - 1];
            }
        }
    }

    /* apply inverse transformation D^{-1}\gamma^{-1}U^{T},
     * since \gamma DAD \approx U^{T}U, so
     * D^{-1}\gamma^{-1}\gamma DADD^{-1} = A \approx D^{-1}\gamma^{-1}U^{T}UD^{-1} */
#if defined(_OPENMP)
    #pragma omp parallel for schedule(dynamic,4096) \
        default(none) shared(U_T_temp, D_inv, gamma) private(i, pj)
#endif
    for ( i = 0; i < U_T_temp.n; ++i )
    {
        for ( pj = U_T_temp.start[i]; pj < U_T_temp.start[i + 1]; ++pj )
        {
            U_T_temp.val[pj] = gamma[i] * (D_inv[i] * U_T_temp.val[pj]);
        }
    }

    /* apply the dropping rule */
    Utop = 0;
    for ( i = 0; i < U_T_temp.n; ++i )
    {
        U_T->start[i] = Utop;

        for ( pj = U_T_temp.start[i]; pj < U_T_temp.start[i + 1] - 1; ++pj )
        {
            if ( FABS( U_T_temp.val[pj] ) > FABS( droptol[i] / U_T_temp.val[U_T_temp.start[i + 1] - 1] ) )
            {
                U_T->j[Utop] = U_T_temp.j[pj];
                U_T->val[Utop] = U_T_temp.val[pj];
                ++Utop;
            }
        }

        /* diagonal */
        U_T->j[Utop] = U_T_temp.j[pj];
        U_T->val[Utop] = U_T_temp.val[pj];
        ++Utop;
    }
    U_T->start[U_T->n] = Utop;

    /* transpose U^{T} and copy into U */
    Transpose( U_T, U );

#if defined(DEBUG_FOCUS)
    fprintf( stderr, "[INFO] nnz(U^T): %d\n", U_T->start[U_T->n] );
    fprintf( stderr, "[INFO] nnz(U): %d\n", U->start[U->n] );
#endif

    Deallocate_Matrix( &U_T_temp );
    Deallocate_Matrix( &DAD );
    sfree( gamma, __FILE__, __LINE__ );
    sfree( D_inv, __FILE__, __LINE__ );
    sfree( D, __FILE__, __LINE__ );

    return Get_Timing_Info( start );
}


/* Fine-grained (parallel) incomplete LU factorization with thresholding
 *
 * Reference:
 * Edmond Chow and Aftab Patel
 * Fine-Grained Parallel Incomplete LU Factorization
 * SIAM J. Sci. Comp.
 *
 * A: symmetric, half-stored (lower triangular), CSR format
 * droptol: row-wise tolerances used for dropping
 * sweeps: number of loops over non-zeros for computation
 * L / U: factorized triangular matrices (A \approx LU), CSR format */
real FG_ILUT( const sparse_matrix * const A, const real * droptol,
               const unsigned int sweeps, sparse_matrix * const L, sparse_matrix * const U )
{
    unsigned int i, pj, x, y, ei_x, ei_y, Ltop, Utop, s;
    real *D, *D_inv, *gamma, sum, start;
    sparse_matrix DAD, L_temp, U_T_temp;

    start = Get_Time( );

    Allocate_Matrix( &DAD, A->n, A->n_max, A->m );
    Allocate_Matrix( &L_temp, A->n, A->n_max, A->m );
    Allocate_Matrix( &U_T_temp, A->n, A->n_max, A->m );

    D = smalloc( sizeof(real) * A->n, __FILE__, __LINE__ );
    D_inv = smalloc( sizeof(real) * A->n, __FILE__, __LINE__ );
    gamma = smalloc( sizeof(real) * A->n, __FILE__, __LINE__ );

#if defined(_OPENMP)
    #pragma omp parallel for schedule(dynamic,512) \
        default(none) shared(D_inv, gamma) private(i) firstprivate(A)
#endif
    for ( i = 0; i < A->n; ++i )
    {
        if ( A->val[A->start[i + 1] - 1] < 0.0 )
        {
            gamma[i] = -1.0;
        }
        else
        {
            gamma[i] = 1.0;
        }

        D_inv[i] = SQRT( FABS( A->val[A->start[i + 1] - 1] ) );
    }

#if defined(_OPENMP)
    #pragma omp parallel for schedule(dynamic,512) \
        default(none) shared(D, D_inv) private(i) firstprivate(A)
#endif
    for ( i = 0; i < A->n; ++i )
    {
        D[i] = 1.0 / D_inv[i];
    }

    /* to get convergence, A must have unit diagonal, so apply
     * transformation \gamma DAD, where D = D(1./SQRT(D(A)))
     * and \gamma = {-1 if a_{ii} < 0, 1 otherwise} */
    memcpy( DAD.start, A->start, sizeof(unsigned int) * (A->n + 1) );
    memcpy( DAD.j, A->j, sizeof(unsigned int) * A->start[A->n] );
#if defined(_OPENMP)
    #pragma omp parallel for schedule(dynamic,4096) \
        default(none) shared(DAD, D, gamma) private(i, pj) firstprivate(A)
#endif
    for ( i = 0; i < A->n; ++i )
    {
        /* non-diagonals */
        for ( pj = A->start[i]; pj < A->start[i + 1] - 1; ++pj )
        {
            DAD.val[pj] = gamma[i] * (D[i] * A->val[pj] * D[A->j[pj]]);
        }
        /* diagonal */
        DAD.val[pj] = 1.0;
    }

    /* initial guesses for L and U,
     * assume: A and DAD symmetric and stored lower triangular */
    memcpy( L_temp.start, DAD.start, sizeof(unsigned int) * (DAD.n + 1) );
    memcpy( L_temp.j, DAD.j, sizeof(unsigned int) * DAD.start[DAD.n] );
    memcpy( L_temp.val, DAD.val, sizeof(real) * DAD.start[DAD.n] );
    /* store U^T in CSR for row-wise access and tranpose later */
    memcpy( U_T_temp.start, DAD.start, sizeof(unsigned int) * (DAD.n + 1) );
    memcpy( U_T_temp.j, DAD.j, sizeof(unsigned int) * DAD.start[DAD.n] );
    memcpy( U_T_temp.val, DAD.val, sizeof(real) * DAD.start[DAD.n] );

    /* L has unit diagonal, by convention */
#if defined(_OPENMP)
    #pragma omp parallel for schedule(dynamic,512) \
        default(none) shared(L_temp) private(i) firstprivate(A)
#endif
    for ( i = 0; i < A->n; ++i )
    {
        L_temp.val[L_temp.start[i + 1] - 1] = 1.0;
    }

    for ( s = 0; s < sweeps; ++s )
    {
        /* for each nonzero in L */
#if defined(_OPENMP)
        #pragma omp parallel for schedule(dynamic,4096) \
            default(none) shared(DAD, L_temp, U_T_temp) private(i, pj, x, y, ei_x, ei_y, sum)
#endif
        for ( pj = 0; pj < L_temp.start[L_temp.n]; ++pj )
        {
            /* determine row bounds of current nonzero */
            x = 0;
            ei_x = pj;
            for ( i = 1; i <= L_temp.n; ++i )
            {
                if ( L_temp.start[i] > pj )
                {
                    --i;
                    x = L_temp.start[i];
                    break;
                }
            }

            /* skip diagonals for L */
            if ( i > L_temp.j[pj] )
            {
                /* determine column bounds of current nonzero */
                y = L_temp.start[L_temp.j[pj]];
                ei_y = L_temp.start[L_temp.j[pj] + 1];

                sum = 0.0;

                /* sparse vector-sparse vector inner products for nonzero (i, j):
                 *   dot( L(i,1:j-1), U^T(j,1:j-1) )
                 *
                 * Note: since L and U^T share the same sparsity pattern
                 * (due to symmetry in A), just use column indices from L */
                for ( ; x < ei_x && y < ei_y && L_temp.j[y] < L_temp.j[pj]; )
                {
                    if ( L_temp.j[x] == L_temp.j[y] )
                    {
                        sum += L_temp.val[x] * U_T_temp.val[y];
                        ++x;
                        ++y;
                    }
                    else if ( L_temp.j[x] < L_temp.j[y] )
                    {
                        ++x;
                    }
                    else
                    {
                        ++y;
                    }
                }

                L_temp.val[pj] = (DAD.val[pj] - sum) / U_T_temp.val[ei_y - 1];
            }
        }

        /* for each nonzero in U^T */
#if defined(_OPENMP)
        #pragma omp parallel for schedule(dynamic,4096) \
            default(none) shared(DAD, L_temp, U_T_temp) private(i, pj, x, y, ei_x, ei_y, sum)
#endif
        for ( pj = 0; pj < L_temp.start[L_temp.n]; ++pj )
        {
            /* determine row bounds of current nonzero */
            x = 0;
            ei_x = pj;
            for ( i = 1; i <= L_temp.n; ++i )
            {
                if ( L_temp.start[i] > pj )
                {
                    --i;
                    x = L_temp.start[i];
                    break;
                }
            }

            /* determine column bounds of current nonzero */
            y = L_temp.start[L_temp.j[pj]];
            ei_y = L_temp.start[L_temp.j[pj] + 1];

            sum = 0.0;

            /* sparse vector-sparse vector inner products for nonzero (i, j):
             *   dot( L(j,1:j-1), U^T(i,1:j-1) )
             *
             * Note: since L and U^T share the same sparsity pattern
             * (due to symmetry in A), just use column indices from L */
            for ( ; x < ei_x && y < ei_y && L_temp.j[y] < L_temp.j[pj]; )
            {
                if ( L_temp.j[x] == L_temp.j[y] )
                {
                    sum += U_T_temp.val[x] * L_temp.val[y];
                    ++x;
                    ++y;
                }
                else if ( L_temp.j[x] < L_temp.j[y] )
                {
                    ++x;
                }
                else
                {
                    ++y;
                }
            }

            U_T_temp.val[pj] = DAD.val[pj] - sum;
        }
    }

    /* apply inverse transformations:
     * since \gamma DAD \approx LU, then
     * D^{-1}DADD^{-1} = A \approx D^{-1}\gamma^{-1}LUD^{-1} */
#if defined(_OPENMP)
    #pragma omp parallel for schedule(dynamic,4096) \
        default(none) shared(DAD, L_temp, D_inv, gamma) private(i, pj)
#endif
    for ( i = 0; i < L_temp.n; ++i )
    {
        for ( pj = L_temp.start[i]; pj < L_temp.start[i + 1]; ++pj )
        {
            L_temp.val[pj] = gamma[i] * (D_inv[i] * L_temp.val[pj]);
        }
    }
    /* note: since we're storing U^T, apply the transform from the left side */
#if defined(_OPENMP)
    #pragma omp parallel for schedule(dynamic,4096) \
        default(none) shared(DAD, U_T_temp, D_inv) private(i, pj)
#endif
    for ( i = 0; i < U_T_temp.n; ++i )
    {
        for ( pj = U_T_temp.start[i]; pj < U_T_temp.start[i + 1]; ++pj )
        {
            U_T_temp.val[pj] = D_inv[i] * U_T_temp.val[pj];
        }
    }

    /* apply the dropping rule */
    Ltop = 0;
    for ( i = 0; i < L_temp.n; ++i )
    {
        L->start[i] = Ltop;

        for ( pj = L_temp.start[i]; pj < L_temp.start[i + 1] - 1; ++pj )
        {
            if ( FABS( L_temp.val[pj] ) > FABS( droptol[i] / L_temp.val[L_temp.start[i + 1] - 1] ) )
            {
                L->j[Ltop] = L_temp.j[pj];
                L->val[Ltop] = L_temp.val[pj];
                ++Ltop;
            }
        }

        /* diagonal */
        L->j[Ltop] = L_temp.j[pj];
        L->val[Ltop] = L_temp.val[pj];
        ++Ltop;
    }
    L->start[L->n] = Ltop;

    Utop = 0;
    for ( i = 0; i < U_T_temp.n; ++i )
    {
        U->start[i] = Utop;

        for ( pj = U_T_temp.start[i]; pj < U_T_temp.start[i + 1] - 1; ++pj )
        {
            if ( FABS( U_T_temp.val[pj] ) > FABS( droptol[i] / U_T_temp.val[U_T_temp.start[i + 1] - 1] ) )
            {
                U->j[Utop] = U_T_temp.j[pj];
                U->val[Utop] = U_T_temp.val[pj];
                ++Utop;
            }
        }

        /* diagonal */
        U->j[Utop] = U_T_temp.j[pj];
        U->val[Utop] = U_T_temp.val[pj];
        ++Utop;
    }
    U->start[U->n] = Utop;

    Transpose_I( U );

#if defined(DEBUG_FOCUS)
    fprintf( stderr, "nnz(L): %d\n", L->start[L->n] );
    fprintf( stderr, "nnz(U): %d\n", U->start[U->n] );
#endif

    Deallocate_Matrix( &U_T_temp );
    Deallocate_Matrix( &L_temp );
    Deallocate_Matrix( &DAD );
    sfree( gamma, __FILE__, __LINE__ );
    sfree( D_inv, __FILE__, __LINE__ );
    sfree( D, __FILE__, __LINE__ );

    return Get_Timing_Info( start );
}


#if defined(HAVE_LAPACKE) || defined(HAVE_LAPACKE_MKL)
/* Compute M^{1} \approx A which minimizes
 *  \sum_{j=1}^{N} ||m_j^T*A - e_j^T||_2^2
 *  where: e_j^T is the j-th row of the NxN identify matrix,
 *         m_j^T is the j-th row of the NxN approximate sparse matrix M
 * 
 * Internally, use LAPACKE to solve the independent least-squares problems.
 * Furthermore, internally solve the related problem
 *  \sum_{j=1}^{N} ||A*m_j - e_j||_2^2
 * for j-th columns m_j and e_j, but store the transpose of M to solve
 * the original problem. That is, for the problems
 *  min ||M_1*A - I||_F^2 and min ||A*M_2 - I||_F^2,
 * it can be shown that if A = A^T, M_1 = M_2^T.
 * Hence, we solve for M_2, and stores its transpose as the result
 * (more efficient for for CSR matrices, row-major storage).
 *
 * A: symmetric, sparse matrix, stored in full CSR format
 * A_spar_patt: sparse matrix used as template sparsity pattern
 *   for approximating the inverse, stored in full CSR format
 * A_app_inv: approximate inverse to A, stored in full CSR format (result)
 *
 * Reference:
 * Michele Benzi et al.
 * A Comparative Study of Sparse Approximate Inverse
 *  Preconditioners
 * Applied Numerical Mathematics 30, 1999
 * */
real sparse_approx_inverse( const sparse_matrix * const A,
        const sparse_matrix * const A_spar_patt, sparse_matrix * A_app_inv )
{
    int i, k, pj, j_temp, identity_pos;
    int N, M, d_i, d_j;
    lapack_int m, n, nrhs, lda, ldb, info;
    int *pos_x, *pos_y;
    real start;
    real *e_j, *dense_matrix;
    size_t e_j_size, dense_matrix_size;
    char *X, *Y;

    start = Get_Time( );

    A_app_inv->start[A_app_inv->n] = A_spar_patt->start[A_spar_patt->n];

#if defined(_OPENMP)
    #pragma omp parallel default(none) \
    private(i, k, pj, j_temp, identity_pos, N, M, d_i, d_j, m, n, \
            nrhs, lda, ldb, info, X, Y, pos_x, pos_y, e_j, \
            dense_matrix, e_j_size, dense_matrix_size) \
    firstprivate(A, A_spar_patt, A_app_inv) \
    shared(stderr)
#endif
    {
        X = smalloc( sizeof(char) * A->n, __FILE__, __LINE__ );
        Y = smalloc( sizeof(char) * A->n, __FILE__, __LINE__ );
        pos_x = smalloc( sizeof(int) * A->n, __FILE__, __LINE__ );
        pos_y = smalloc( sizeof(int) * A->n, __FILE__, __LINE__ );

        e_j = NULL;
        dense_matrix = NULL;
        e_j_size = 0;
        dense_matrix_size = 0;

#if defined(_OPENMP)
        #pragma omp for schedule(dynamic,64)
#endif
        for ( i = 0; i < A_spar_patt->n; ++i )
        {
            N = 0;
            M = 0;

            for ( k = 0; k < A->n; ++k )
            {
                X[k] = 0;
            }
            for ( k = 0; k < A->n; ++k )
            {
                Y[k] = 0;
            }
            for ( k = 0; k < A->n; ++k )
            {
                pos_x[k] = 0;
            }
            for ( k = 0; k < A->n; ++k )
            {
                pos_y[k] = 0;
            }

            /* find column indices of nonzeros
             * (which will be the columns indices of the dense matrix) */
            for ( pj = A_spar_patt->start[i]; pj < A_spar_patt->start[i + 1]; ++pj )
            {
                j_temp = A_spar_patt->j[pj];

                Y[j_temp] = 1;
                pos_y[j_temp] = N;
                ++N;

                /* for each of those indices:
                 * search through the row of full A of that index */
                for ( k = A->start[j_temp]; k < A->start[j_temp + 1]; ++k )
                {
                    /* and accumulate the nonzero column indices
                     * to serve as the row indices of the dense matrix */
                    X[A->j[k]] = 1;
                }
            }

            /* enumerate the row indices from 0 to (# of nonzero rows - 1) for the dense matrix */
            identity_pos = M;
            for ( k = 0; k < A->n; k++)
            {
                if ( X[k] != 0 )
                {
                    pos_x[M] = k;
                    if ( k == i )
                    {
                        identity_pos = M;
                    }
                    ++M;
                }
            }

            /* N x M dense matrix */
            if ( dense_matrix == NULL )
            {
                dense_matrix = smalloc( sizeof(real) * N * M, __FILE__, __LINE__ );
                dense_matrix_size = sizeof(real) * N * M;
            }
            else if ( dense_matrix_size < sizeof(real) * N * M )
            {
                sfree( dense_matrix, __FILE__, __LINE__ );
                dense_matrix = smalloc( sizeof(real) * N * M, __FILE__, __LINE__ );
                dense_matrix_size = sizeof(real) * N * M;
            }

            /* fill in the entries of dense matrix */
            for ( d_i = 0; d_i < M; ++d_i)
            {
                for ( d_j = 0; d_j < N; ++d_j )
                {
                    dense_matrix[d_i * N + d_j] = 0.0;
                }

                /* change the value if any of the column indices is seen */
                for ( d_j = A->start[pos_x[d_i]];
                        d_j < A->start[pos_x[d_i] + 1]; ++d_j )
                {
                    if ( Y[A->j[d_j]] == 1 )
                    {
                        dense_matrix[d_i * N + pos_y[A->j[d_j]]] = A->val[d_j];
                    }
                }

            }

            /* create the right hand side of the linear equation
             * that is the full column of the identity matrix */
            if ( e_j == NULL )
            {
                e_j = smalloc( sizeof(real) * M, __FILE__, __LINE__ );
                e_j_size = sizeof(real) * M;
            }
            else if ( e_j_size < sizeof(real) * M )
            {
                sfree( e_j, __FILE__, __LINE__  );
                e_j = smalloc( sizeof(real) * M, __FILE__, __LINE__ );
                e_j_size = sizeof(real) * M;
            }

            for ( k = 0; k < M; ++k )
            {
                e_j[k] = 0.0;
            }
            e_j[identity_pos] = 1.0;

            /* Solve the overdetermined system AX = B through the least-squares problem:
             * min ||B - AX||_2 */
            m = M;
            n = N;
            nrhs = 1;
            lda = N;
            ldb = nrhs;
            info = LAPACKE_dgels( LAPACK_ROW_MAJOR, 'N', m, n, nrhs, dense_matrix, lda,
                    e_j, ldb );

            /* Check for the full rank */
            if ( info > 0 )
            {
                fprintf( stderr, "[ERROR] The diagonal element %i of the triangular factor ", info );
                fprintf( stderr, "of A is zero, so that A does not have full rank;\n" );
                fprintf( stderr, "the least squares solution could not be computed.\n" );
                exit( INVALID_INPUT );
            }

            /* Print least squares solution */
//            print_matrix( "Least squares solution", n, nrhs, b, ldb );

            /* accumulate the resulting vector to build A_app_inv */
            A_app_inv->start[i] = A_spar_patt->start[i];
            for ( k = A_spar_patt->start[i]; k < A_spar_patt->start[i + 1]; ++k)
            {
                A_app_inv->j[k] = A_spar_patt->j[k];
                A_app_inv->val[k] = e_j[k - A_spar_patt->start[i]];
            }
        }

        sfree( dense_matrix, __FILE__, __LINE__ );
        sfree( e_j, __FILE__, __LINE__  );
        sfree( pos_y, __FILE__, __LINE__ );
        sfree( pos_x, __FILE__, __LINE__ );
        sfree( Y, __FILE__, __LINE__ );
        sfree( X, __FILE__, __LINE__ );
    }

    return Get_Timing_Info( start );
}
#endif


/* sparse matrix, dense vector multiplication Ax = b
 *
 * workspace: storage container for workspace structures
 * A: symmetric (lower triangular portion only stored), square matrix,
 *    stored in CSR format
 * x: dense vector, size equal to num. columns in A
 * b (output): dense vector, size equal to num. columns in A */
static void sparse_matvec( const static_storage * const workspace,
        const sparse_matrix * const A, const real * const x, real * const b )
{
    int i, j, k, n, si, ei;
    real H;
#if defined(_OPENMP)
    unsigned int tid;
#endif

    n = A->n;
    Vector_MakeZero( b, n );

#if defined(_OPENMP)
    tid = omp_get_thread_num( );

    Vector_MakeZero( workspace->b_local, omp_get_num_threads() * n );

    #pragma omp for schedule(guided)
#endif
    for ( i = 0; i < n; ++i )
    {
        si = A->start[i];
        ei = A->start[i + 1] - 1;

        for ( k = si; k < ei; ++k )
        {
            j = A->j[k];
            H = A->val[k];
#if defined(_OPENMP)
            workspace->b_local[tid * n + j] += H * x[i];
            workspace->b_local[tid * n + i] += H * x[j];
#else
            b[j] += H * x[i];
            b[i] += H * x[j];
#endif
        }

        // the diagonal entry is the last one in
#if defined(_OPENMP)
        workspace->b_local[tid * n + i] += A->val[k] * x[i];
#else
        b[i] += A->val[k] * x[i];
#endif
    }

#if defined(_OPENMP)
    #pragma omp for schedule(dynamic,256)
    for ( i = 0; i < n; ++i )
    {
        for ( j = 0; j < omp_get_num_threads(); ++j )
        {
            b[i] += workspace->b_local[j * n + i];
        }
    }
#endif
}


/* sparse matrix, dense vector multiplication Ax = b
 *
 * A: square matrix, stored in CSR format
 * x: dense vector, size equal to num. columns in A
 * b (output): dense vector, size equal to num. columns in A */
static void sparse_matvec_full( const sparse_matrix * const A,
        const real * const x, real * const b )
{
    int i, pj;

    Vector_MakeZero( b, A->n );

#if defined(_OPENMP)
    #pragma omp for schedule(guided)
#endif
    for ( i = 0; i < A->n; ++i )
    {
        for ( pj = A->start[i]; pj < A->start[i + 1]; ++pj )
        {
            b[i] += A->val[pj] * x[A->j[pj]];
        }
    }
}


/* Transpose A and copy into A^T
 *
 * A: stored in CSR
 * A_t: stored in CSR
 */
void Transpose( const sparse_matrix * const A, sparse_matrix * const A_t )
{
    unsigned int i, j, pj, *A_t_top;

    A_t_top = scalloc( A->n + 1, sizeof(unsigned int), __FILE__, __LINE__ );

    for ( i = 0; i < A->n + 1; ++i )
    {
        A_t->start[i] = 0;
    }

    /* count nonzeros in each column of A^T, store one row greater (see next loop) */
    for ( i = 0; i < A->n; ++i )
    {
        for ( pj = A->start[i]; pj < A->start[i + 1]; ++pj )
        {
            ++A_t->start[A->j[pj] + 1];
        }
    }

    /* setup the row pointers for A^T */
    for ( i = 1; i <= A->n; ++i )
    {
        A_t->start[i] = A_t->start[i] + A_t->start[i - 1];
        A_t_top[i] = A_t->start[i];
    }

    /* fill in A^T */
    for ( i = 0; i < A->n; ++i )
    {
        for ( pj = A->start[i]; pj < A->start[i + 1]; ++pj )
        {
            j = A->j[pj];
            A_t->j[A_t_top[j]] = i;
            A_t->val[A_t_top[j]] = A->val[pj];
            ++A_t_top[j];
        }
    }

    sfree( A_t_top, __FILE__, __LINE__ );
}


/* Transpose A in-place
 *
 * A: stored in CSR
 */
void Transpose_I( sparse_matrix * const A )
{
    sparse_matrix A_t;

    Allocate_Matrix( &A_t, A->n, A->n_max, A->m );

    Transpose( A, &A_t );

    memcpy( A->start, A_t.start, sizeof(int) * (A_t.n + 1) );
    memcpy( A->j, A_t.j, sizeof(int) * (A_t.start[A_t.n]) );
    memcpy( A->val, A_t.val, sizeof(real) * (A_t.start[A_t.n]) );

    Deallocate_Matrix( &A_t );
}


/* Apply diagonal inverse (Jacobi) preconditioner to system residual
 *
 * Hdia_inv: diagonal inverse preconditioner (constructed using H)
 * y: current residual
 * x: preconditioned residual
 * N: dimensions of preconditioner and vectors (# rows in H)
 */
static void jacobi_app( const real * const Hdia_inv, const real * const y,
                          real * const x, const int N )
{
    unsigned int i;

#if defined(_OPENMP)
    #pragma omp for schedule(dynamic,256)
#endif
    for ( i = 0; i < N; ++i )
    {
        x[i] = y[i] * Hdia_inv[i];
    }
}


/* Solve triangular system LU*x = y using forward/backward substitution
 *
 * LU: lower/upper triangular, stored in CSR
 * y: constants in linear system (RHS)
 * x: solution
 * N: dimensions of matrix and vectors
 * tri: triangularity of LU (lower/upper)
 *
 * Assumptions:
 *   LU has non-zero diagonals
 *   Each row of LU has at least one non-zero (i.e., no rows with all zeros) */
void tri_solve( const sparse_matrix * const LU, const real * const y,
        real * const x, const TRIANGULARITY tri )
{
    int i, pj, j, si, ei;
    real val;

#if defined(_OPENMP)
    #pragma omp single
#endif
    {
        if ( tri == LOWER )
        {
            for ( i = 0; i < LU->n; ++i )
            {
                x[i] = y[i];
                si = LU->start[i];
                ei = LU->start[i + 1];
                for ( pj = si; pj < ei - 1; ++pj )
                {
                    j = LU->j[pj];
                    val = LU->val[pj];
                    x[i] -= val * x[j];
                }
                x[i] /= LU->val[pj];
            }
        }
        else
        {
            for ( i = LU->n - 1; i >= 0; --i )
            {
                x[i] = y[i];
                si = LU->start[i];
                ei = LU->start[i + 1];
                for ( pj = si + 1; pj < ei; ++pj )
                {
                    j = LU->j[pj];
                    val = LU->val[pj];
                    x[i] -= val * x[j];
                }
                x[i] /= LU->val[si];
            }
        }
    }
}


/* Solve triangular system LU*x = y using level scheduling
 *
 * workspace: storage container for workspace structures
 * LU: lower/upper triangular, stored in CSR
 * y: constants in linear system (RHS)
 * x (output): solution to triangular system
 * N: dimensions of matrix and vectors
 * tri: triangularity of LU (lower/upper)
 * find_levels: perform level search if positive, otherwise reuse existing levels
 *
 * Assumptions:
 *   LU has non-zero diagonals
 *   Each row of LU has at least one non-zero (i.e., no rows with all zeros) */
void tri_solve_level_sched( static_storage * workspace,
        const sparse_matrix * const LU, const real * const y, real * const x,
        const TRIANGULARITY tri, int find_levels )
{
    int i, j, pj, local_row, local_level;
    unsigned int *row_levels, *level_rows, *level_rows_cnt;
    int levels;

    if ( tri == LOWER )
    {
        row_levels = workspace->row_levels_L;
        level_rows = workspace->level_rows_L;
        level_rows_cnt = workspace->level_rows_cnt_L;
    }
    else
    {
        row_levels = workspace->row_levels_U;
        level_rows = workspace->level_rows_U;
        level_rows_cnt = workspace->level_rows_cnt_U;
    }

#if defined(_OPENMP)
    #pragma omp single
#endif
    {
        /* find levels (row dependencies in substitutions) */
        if ( find_levels == TRUE )
        {
            for ( i = 0; i < LU->n; ++i )
                row_levels[i] = 0;
            for ( i = 0; i < LU->n; ++i )
                level_rows_cnt[i] = 0;
            for ( i = 0; i < LU->n; ++i )
                workspace->top[i] = 0;
            levels = 1;

            if ( tri == LOWER )
            {
                for ( i = 0; i < LU->n; ++i )
                {
                    local_level = 1;
                    for ( pj = LU->start[i]; pj < LU->start[i + 1] - 1; ++pj )
                    {
                        local_level = MAX( local_level, row_levels[LU->j[pj]] + 1 );
                    }

                    levels = MAX( levels, local_level );
                    row_levels[i] = local_level;
                    ++level_rows_cnt[local_level];
                }

                workspace->levels_L = levels;

#if defined(DEBUG)
                fprintf( stderr, "[INFO] levels(L): %d\n", levels );
                fprintf( stderr, "[INFO] NNZ(L): %d\n", LU->start[LU->n] );
#endif
            }
            else
            {
                for ( i = LU->n - 1; i >= 0; --i )
                {
                    local_level = 1;
                    for ( pj = LU->start[i] + 1; pj < LU->start[i + 1]; ++pj )
                    {
                        local_level = MAX( local_level, row_levels[LU->j[pj]] + 1 );
                    }

                    levels = MAX( levels, local_level );
                    row_levels[i] = local_level;
                    ++level_rows_cnt[local_level];
                }

                workspace->levels_U = levels;

#if defined(DEBUG)
                fprintf( stderr, "[INFO] levels(U): %d\n", levels );
                fprintf( stderr, "[INFO] NNZ(U): %d\n", LU->start[LU->n] );
#endif
            }

            for ( i = 1; i < levels + 1; ++i )
            {
                level_rows_cnt[i] += level_rows_cnt[i - 1];
                workspace->top[i] = level_rows_cnt[i];
            }

            for ( i = 0; i < LU->n; ++i )
            {
                level_rows[workspace->top[row_levels[i] - 1]] = i;
                ++workspace->top[row_levels[i] - 1];
            }
        }
    }

    if ( tri == LOWER )
    {
        levels = workspace->levels_L;
    }
    else
    {
        levels = workspace->levels_U;
    }

    /* perform substitutions by level */
    if ( tri == LOWER )
    {
        for ( i = 0; i < levels; ++i )
        {
#if defined(_OPENMP)
            #pragma omp for schedule(static)
#endif
            for ( j = level_rows_cnt[i]; j < level_rows_cnt[i + 1]; ++j )
            {
                local_row = level_rows[j];
                x[local_row] = y[local_row];
                for ( pj = LU->start[local_row]; pj < LU->start[local_row + 1] - 1; ++pj )
                {
                    x[local_row] -= LU->val[pj] * x[LU->j[pj]];
                }
                x[local_row] /= LU->val[pj];
            }
        }
    }
    else
    {
        for ( i = 0; i < levels; ++i )
        {
#if defined(_OPENMP)
            #pragma omp for schedule(static)
#endif
            for ( j = level_rows_cnt[i]; j < level_rows_cnt[i + 1]; ++j )
            {
                local_row = level_rows[j];
                x[local_row] = y[local_row];
                for ( pj = LU->start[local_row] + 1; pj < LU->start[local_row + 1]; ++pj )
                {
                    x[local_row] -= LU->val[pj] * x[LU->j[pj]];
                }
                x[local_row] /= LU->val[LU->start[local_row]];
            }
        }
    }
}


/* Iterative greedy shared-memory parallel graph coloring
 *
 * control: container for control info
 * workspace: storage container for workspace structures
 * A: matrix to use for coloring, stored in CSR format;
 *   rows represent vertices, columns of entries within a row represent adjacent vertices
 *   (i.e., dependent rows for elimination during LU factorization)
 * tri: triangularity of LU (lower/upper)
 * color: vertex color (1-based)
 *
 * Reference:
 * Umit V. Catalyurek et al.
 * Graph Coloring Algorithms for Multi-core
 *  and Massively Threaded Architectures
 * Parallel Computing, 2012
 */
void graph_coloring( const control_params * const control,
        static_storage * workspace,
        const sparse_matrix * const A, const TRIANGULARITY tri )
{
#if defined(_OPENMP)
    #pragma omp parallel
#endif
    {
        int i, pj, v;
        unsigned int temp, recolor_cnt_local, *conflict_local;
        int tid, *fb_color;
        unsigned int *p_to_color, *p_conflict, *p_temp;

#if defined(_OPENMP)
        tid = omp_get_thread_num( );
#else
        tid = 0;
#endif
        p_to_color = workspace->to_color;
        p_conflict = workspace->conflict;

#if defined(_OPENMP)
        #pragma omp for schedule(static)
#endif
        for ( i = 0; i < A->n; ++i )
        {
            workspace->color[i] = 0;
        }

#if defined(_OPENMP)
        #pragma omp single
#endif
        {
            workspace->recolor_cnt = A->n;
        }

        /* ordering of vertices to color depends on triangularity of factor
         * for which coloring is to be used for */
        if ( tri == LOWER )
        {
#if defined(_OPENMP)
            #pragma omp for schedule(static)
#endif
            for ( i = 0; i < A->n; ++i )
            {
                p_to_color[i] = i;
            }
        }
        else
        {
#if defined(_OPENMP)
            #pragma omp for schedule(static)
#endif
            for ( i = 0; i < A->n; ++i )
            {
                p_to_color[i] = A->n - 1 - i;
            }
        }

        fb_color = smalloc( sizeof(int) * A->n, __FILE__, __LINE__ );
        conflict_local = smalloc( sizeof(unsigned int) * A->n, __FILE__, __LINE__ );

        while ( workspace->recolor_cnt > 0 )
        {
            for ( i = 0; i < A->n; ++i )
                fb_color[i] = -1;

            /* color vertices */
#if defined(_OPENMP)
            #pragma omp for schedule(static)
#endif
            for ( i = 0; i < workspace->recolor_cnt; ++i )
            {
                v = p_to_color[i];

                /* colors of adjacent vertices are forbidden */
                for ( pj = A->start[v]; pj < A->start[v + 1]; ++pj )
                {
                    if ( v != A->j[pj] )
                    {
                        fb_color[workspace->color[A->j[pj]]] = v;
                    }
                }

                /* search for min. color which is not in conflict with adjacent vertices;
                 * start at 1 since 0 is default (invalid) color for all vertices */
                for ( pj = 1; fb_color[pj] == v; ++pj )
                    ;

                /* assign discovered color (no conflict in neighborhood of adjacent vertices) */
                workspace->color[v] = pj;
            }

            /* determine if recoloring required */
            temp = workspace->recolor_cnt;
            recolor_cnt_local = 0;

#if defined(_OPENMP)
            #pragma omp barrier

            #pragma omp single
#endif
            {
                workspace->recolor_cnt = 0;
            }

#if defined(_OPENMP)
            #pragma omp for schedule(static)
#endif
            for ( i = 0; i < temp; ++i )
            {
                v = p_to_color[i];

                /* search for color conflicts with adjacent vertices */
                for ( pj = A->start[v]; pj < A->start[v + 1]; ++pj )
                {
                    if ( workspace->color[v] == workspace->color[A->j[pj]] && v > A->j[pj] )
                    {
                        conflict_local[recolor_cnt_local] = v;
                        ++recolor_cnt_local;
                        break;
                    }
                }
            }

            /* count thread-local conflicts and compute offsets for copying into shared buffer */
            workspace->conflict_cnt[tid + 1] = recolor_cnt_local;

#if defined(_OPENMP)
            #pragma omp barrier

            #pragma omp single
#endif
            {
                workspace->conflict_cnt[0] = 0;
                for ( i = 1; i < control->num_threads + 1; ++i )
                {
                    workspace->conflict_cnt[i] += workspace->conflict_cnt[i - 1];
                }
                workspace->recolor_cnt = workspace->conflict_cnt[control->num_threads];
            }

            /* copy thread-local conflicts into shared buffer */
            for ( i = 0; i < recolor_cnt_local; ++i )
            {
                p_conflict[workspace->conflict_cnt[tid] + i] = conflict_local[i];
                workspace->color[conflict_local[i]] = 0;
            }

#if defined(_OPENMP)
            #pragma omp barrier
#endif
            p_temp = p_to_color;
            p_to_color = p_conflict;
            p_conflict = p_temp;
        }

        sfree( conflict_local, __FILE__, __LINE__ );
        sfree( fb_color, __FILE__, __LINE__ );
    }
}


/* Sort rows by coloring
 *
 * workspace: storage container for workspace structures
 * n: number of entries in coloring
 * tri: coloring to triangular factor to use (lower/upper)
 */
void sort_rows_by_colors( const static_storage * const workspace,
        const unsigned int n, const TRIANGULARITY tri )
{
    unsigned int i;

    for ( i = 0; i < n + 1; ++i )
        workspace->color_top[i] = 0;

    /* sort vertices by color (ascending within a color)
     *  1) count colors
     *  2) determine offsets of color ranges
     *  3) sort rows by color
     *
     *  note: color is 1-based */
    for ( i = 0; i < n; ++i )
    {
        ++workspace->color_top[workspace->color[i]];
    }
    for ( i = 1; i < n + 1; ++i )
    {
        workspace->color_top[i] += workspace->color_top[i - 1];
    }
    for ( i = 0; i < n; ++i )
    {
        workspace->permuted_row_col[workspace->color_top[workspace->color[i] - 1]] = i;
        ++workspace->color_top[workspace->color[i] - 1];
    }

    /* invert mapping to get map from current row/column to permuted (new) row/column */
    for ( i = 0; i < n; ++i )
    {
        workspace->permuted_row_col_inv[workspace->permuted_row_col[i]] = i;
    }
}


/* Apply permutation P*x = x_p
 *
 * x_p (output): permuted vector
 * x: vector to permute
 * perm: dense vector representing Permutation matrix for computing x_p = P*x;
 *  this maps i-th position to its permuted position
 * n: number of entries in x
 */
static void permute_vector( real * const x_p, const real * const x,
        const unsigned int * const perm, const unsigned int n )
{
    unsigned int i;

#if defined(_OPENMP)
    #pragma omp for schedule(static)
#endif
    for ( i = 0; i < n; ++i )
    {
        x_p[i] = x[perm[i]];
    }
}


/* Apply permutation Q^T*(LU)*Q based on graph coloring
 *
 * workspace: storage container for workspace structures
 * color: vertex color (1-based); vertices represent matrix rows/columns
 * LU: matrix to permute, stored in CSR format
 * tri: triangularity of LU (lower/upper)
 */
static void permute_matrix( const static_storage * const workspace,
        sparse_matrix * const LU, const TRIANGULARITY tri )
{
    int i, pj, nr, nc;
    sparse_matrix LUtemp;

    Allocate_Matrix( &LUtemp, LU->n, LU->n_max, LU->m );

    for ( i = 0; i < LU->n + 1; ++i )
        workspace->color_top[i] = 0;

    /* count nonzeros in each row of permuted factor (re-use color_top for counting) */
    if ( tri == LOWER )
    {
        for ( i = 0; i < LU->n; ++i )
        {
            nr = workspace->permuted_row_col_inv[i];

            for ( pj = LU->start[i]; pj < LU->start[i + 1]; ++pj )
            {
                nc = workspace->permuted_row_col_inv[LU->j[pj]];

                if ( nc <= nr )
                {
                    ++workspace->color_top[nr + 1];
                }
                /* correct entries to maintain triangularity (lower) */
                else
                {
                    ++workspace->color_top[nc + 1];
                }
            }
        }
    }
    else
    {
        for ( i = LU->n - 1; i >= 0; --i )
        {
            nr = workspace->permuted_row_col_inv[i];

            for ( pj = LU->start[i]; pj < LU->start[i + 1]; ++pj )
            {
                nc = workspace->permuted_row_col_inv[LU->j[pj]];

                if ( nc >= nr )
                {
                    ++workspace->color_top[nr + 1];
                }
                /* correct entries to maintain triangularity (upper) */
                else
                {
                    ++workspace->color_top[nc + 1];
                }
            }
        }
    }

    for ( i = 1; i < LU->n + 1; ++i )
    {
        workspace->color_top[i] += workspace->color_top[i - 1];
    }

    memcpy( LUtemp.start, workspace->color_top, sizeof(unsigned int) * (LU->n + 1) );

    /* permute factor */
    if ( tri == LOWER )
    {
        for ( i = 0; i < LU->n; ++i )
        {
            nr = workspace->permuted_row_col_inv[i];

            for ( pj = LU->start[i]; pj < LU->start[i + 1]; ++pj )
            {
                nc = workspace->permuted_row_col_inv[LU->j[pj]];

                if ( nc <= nr )
                {
                    LUtemp.j[workspace->color_top[nr]] = nc;
                    LUtemp.val[workspace->color_top[nr]] = LU->val[pj];
                    ++workspace->color_top[nr];
                }
                /* correct entries to maintain triangularity (lower) */
                else
                {
                    LUtemp.j[workspace->color_top[nc]] = nr;
                    LUtemp.val[workspace->color_top[nc]] = LU->val[pj];
                    ++workspace->color_top[nc];
                }
            }
        }
    }
    else
    {
        for ( i = LU->n - 1; i >= 0; --i )
        {
            nr = workspace->permuted_row_col_inv[i];

            for ( pj = LU->start[i]; pj < LU->start[i + 1]; ++pj )
            {
                nc = workspace->permuted_row_col_inv[LU->j[pj]];

                if ( nc >= nr )
                {
                    LUtemp.j[workspace->color_top[nr]] = nc;
                    LUtemp.val[workspace->color_top[nr]] = LU->val[pj];
                    ++workspace->color_top[nr];
                }
                /* correct entries to maintain triangularity (upper) */
                else
                {
                    LUtemp.j[workspace->color_top[nc]] = nr;
                    LUtemp.val[workspace->color_top[nc]] = LU->val[pj];
                    ++workspace->color_top[nc];
                }
            }
        }
    }

    memcpy( LU->start, LUtemp.start, sizeof(unsigned int) * (LU->n + 1) );
    memcpy( LU->j, LUtemp.j, sizeof(unsigned int) * LU->start[LU->n] );
    memcpy( LU->val, LUtemp.val, sizeof(real) * LU->start[LU->n] );

    Deallocate_Matrix( &LUtemp );
}


/* Setup routines to build permuted QEq matrix H (via graph coloring),
 *  used for preconditioning (incomplete factorizations computed based on
 *  permuted H)
 *
 * control: container for control info
 * workspace: storage container for workspace structures
 * H: symmetric, lower triangular portion only, stored in CSR format;
 * H_full: symmetric, stored in CSR format;
 * H_p (output): permuted copy of H based on coloring, lower half stored, CSR format
 */
void setup_graph_coloring( const control_params * const control,
        const static_storage * const workspace, const sparse_matrix * const H,
        sparse_matrix * H_full, sparse_matrix * H_p, int realloc )
{
    if ( H_p->allocated == FALSE )
    {
        Allocate_Matrix( H_p, H->n, H->n_max, H->m );
    }
    else if ( H_p->m < H->m || H_p->n_max < H->n_max || realloc == TRUE )
    {
        if ( H_p->allocated == TRUE )
        {
            Deallocate_Matrix( H_p );
        }
        Allocate_Matrix( H_p, H->n, H->n_max, H->m );
    }

    compute_full_sparse_matrix( H, H_full, realloc );

    graph_coloring( control, (static_storage *) workspace, H_full, LOWER );
    sort_rows_by_colors( workspace, H_full->n, LOWER );

    memcpy( H_p->start, H->start, sizeof(int) * (H->n + 1) );
    memcpy( H_p->j, H->j, sizeof(int) * (H->start[H->n]) );
    memcpy( H_p->val, H->val, sizeof(real) * (H->start[H->n]) );
    permute_matrix( workspace, H_p, LOWER );
}


/* Jacobi iteration using truncated Neumann series: x_{k+1} = Gx_k + D^{-1}b
 * where:
 *   G = I - D^{-1}R
 *   R = triangular matrix
 *   D = diagonal matrix, diagonals from R
 *
 * Note: used during the backsolves when applying preconditioners with
 * triangular factors in iterative linear solvers
 *
 * Note: Newmann series arises from series expansion of the inverse of
 * the coefficient matrix in the triangular system
 *
 * workspace: storage container for workspace structures
 * R:
 * Dinv:
 * b:
 * x (output):
 * tri:
 * maxiter:
 * */
void jacobi_iter( const static_storage * const workspace,
        const sparse_matrix * const R, const real * const Dinv,
        const real * const b, real * const x, const TRIANGULARITY tri,
        const unsigned int maxiter )
{
    unsigned int i, k, si, ei, iter;
    real *p1, *p2, *p3;

    si = 0;
    ei = 0;
    iter = 0;
    p1 = workspace->rp;
    p2 = workspace->rp2;

    Vector_MakeZero( p1, R->n );

    /* precompute and cache, as invariant in loop below */
#if defined(_OPENMP)
    #pragma omp for schedule(static)
#endif
    for ( i = 0; i < R->n; ++i )
    {
        workspace->Dinv_b[i] = Dinv[i] * b[i];
    }

    do
    {
        /* x_{k+1} = G*x_{k} + Dinv*b */
#if defined(_OPENMP)
        #pragma omp for schedule(guided)
#endif
        for ( i = 0; i < R->n; ++i )
        {
            if (tri == LOWER)
            {
                si = R->start[i];
                ei = R->start[i + 1] - 1;
            }
            else
            {

                si = R->start[i] + 1;
                ei = R->start[i + 1];
            }

            p2[i] = 0.;

            for ( k = si; k < ei; ++k )
            {
                p2[i] += R->val[k] * p1[R->j[k]];
            }

            p2[i] *= -Dinv[i];
            p2[i] += workspace->Dinv_b[i];
        }

        p3 = p1;
        p1 = p2;
        p2 = p3;

        ++iter;
    }
    while ( iter < maxiter );

    Vector_Copy( x, p1, R->n );
}


/* Apply left-sided preconditioning while solving M^{-1}Ax = M^{-1}b
 *
 * workspace: data struct containing matrices, stored in CSR
 * control: data struct containing parameters
 * y: vector to which to apply preconditioning,
 *  specific to internals of iterative solver being used
 * x (output): preconditioned vector
 * fresh_pre: parameter indicating if this is a newly computed (fresh) preconditioner
 * side: used in determining how to apply preconditioner if the preconditioner is
 *  factorized as M = M_{1}M_{2} (e.g., incomplete LU, A \approx LU)
 *
 * Assumptions:
 *   Matrices have non-zero diagonals
 *   Each row of a matrix has at least one non-zero (i.e., no rows with all zeros) */
static void apply_preconditioner( const static_storage * const workspace,
        const control_params * const control, const real * const y, real * const x,
        int fresh_pre, int side )
{
    int i, si;

    /* no preconditioning */
    if ( control->cm_solver_pre_comp_type == NONE_PC )
    {
        if ( x != y )
        {
            Vector_Copy( x, y, workspace->H.n );
        }
    }
    else
    {
        switch ( side )
        {
        case LEFT:
            switch ( control->cm_solver_pre_app_type )
            {
            case TRI_SOLVE_PA:
                switch ( control->cm_solver_pre_comp_type )
                {
                case JACOBI_PC:
                    jacobi_app( workspace->Hdia_inv, y, x, workspace->H.n );
                    break;
                case ICHOLT_PC:
                case ILUT_PC:
                case FG_ILUT_PC:
                    tri_solve( &workspace->L, y, x, LOWER );
                    break;
                case ILUTP_PC:
                    permute_vector( workspace->y_p, y, workspace->perm_ilutp, workspace->H.n );
                    tri_solve( &workspace->L, workspace->y_p, x, LOWER );
                    break;
                case SAI_PC:
                    sparse_matvec_full( &workspace->H_app_inv, y, x );
                    break;
                default:
                    fprintf( stderr, "[ERROR] Unrecognized preconditioner application method. Terminating...\n" );
                    exit( INVALID_INPUT );
                    break;
                }
                break;
            case TRI_SOLVE_LEVEL_SCHED_PA:
                switch ( control->cm_solver_pre_comp_type )
                {
                case JACOBI_PC:
                    jacobi_app( workspace->Hdia_inv, y, x, workspace->H.n );
                    break;
                case ICHOLT_PC:
                case ILUT_PC:
                case FG_ILUT_PC:
                    tri_solve_level_sched( (static_storage *) workspace,
                            &workspace->L, y, x, LOWER, fresh_pre );
                    break;
                case ILUTP_PC:
                    permute_vector( workspace->y_p, y, workspace->perm_ilutp, workspace->L.n );
                    tri_solve_level_sched( (static_storage *) workspace,
                            &workspace->L, workspace->y_p, x, LOWER, fresh_pre );
                    break;
                case SAI_PC:
                    sparse_matvec_full( &workspace->H_app_inv, y, x );
                    break;
                default:
                    fprintf( stderr, "[ERROR] Unrecognized preconditioner application method. Terminating...\n" );
                    exit( INVALID_INPUT );
                    break;
                }
                break;
            case TRI_SOLVE_GC_PA:
                switch ( control->cm_solver_pre_comp_type )
                {
                case JACOBI_PC:
                case SAI_PC:
                    fprintf( stderr, "[ERROR] Unsupported preconditioner computation/application method combination. Terminating...\n" );
                    exit( INVALID_INPUT );
                    break;
                case ICHOLT_PC:
                case ILUT_PC:
                case FG_ILUT_PC:
                    permute_vector( workspace->y_p, y, workspace->permuted_row_col, workspace->L.n );
                    tri_solve_level_sched( (static_storage *) workspace,
                            &workspace->L, workspace->y_p, x, LOWER, fresh_pre );
                    break;
                case ILUTP_PC:
                    permute_vector( workspace->y_p, y, workspace->permuted_row_col, workspace->L.n );
                    permute_vector( workspace->x_p, workspace->y_p, workspace->perm_ilutp, workspace->L.n );
                    tri_solve_level_sched( (static_storage *) workspace,
                            &workspace->L, workspace->x_p, x, LOWER, fresh_pre );
                    break;
                default:
                    fprintf( stderr, "[ERROR] Unrecognized preconditioner application method. Terminating...\n" );
                    exit( INVALID_INPUT );
                    break;
                }
                break;
            case JACOBI_ITER_PA:
                switch ( control->cm_solver_pre_comp_type )
                {
                case JACOBI_PC:
                case SAI_PC:
                    fprintf( stderr, "[ERROR] Unsupported preconditioner computation/application method combination. Terminating...\n" );
                    exit( INVALID_INPUT );
                    break;
                case ICHOLT_PC:
                case ILUT_PC:
                case FG_ILUT_PC:
                    /* construct D^{-1}_L */
                    if ( fresh_pre == TRUE )
                    {
#if defined(_OPENMP)
                        #pragma omp for schedule(static)
#endif
                        for ( i = 0; i < workspace->L.n; ++i )
                        {
                            si = workspace->L.start[i + 1] - 1;
                            workspace->Dinv_L[i] = 1.0 / workspace->L.val[si];
                        }
                    }

                    jacobi_iter( workspace, &workspace->L, workspace->Dinv_L,
                            y, x, LOWER, control->cm_solver_pre_app_jacobi_iters );
                    break;
                case ILUTP_PC:
                    permute_vector( workspace->y_p, y, workspace->perm_ilutp, workspace->H.n );

                    /* construct D^{-1}_L */
                    if ( fresh_pre == TRUE )
                    {
#if defined(_OPENMP)
                        #pragma omp for schedule(static)
#endif
                        for ( i = 0; i < workspace->L.n; ++i )
                        {
                            si = workspace->L.start[i + 1] - 1;
                            workspace->Dinv_L[i] = 1.0 / workspace->L.val[si];
                        }
                    }

                    jacobi_iter( workspace, &workspace->L, workspace->Dinv_L,
                            workspace->y_p, x, LOWER, control->cm_solver_pre_app_jacobi_iters );
                    break;
                default:
                    fprintf( stderr, "[ERROR] Unrecognized preconditioner application method. Terminating...\n" );
                    exit( INVALID_INPUT );
                    break;
                }
                break;
            default:
                fprintf( stderr, "[ERROR] Unrecognized preconditioner application method. Terminating...\n" );
                exit( INVALID_INPUT );
                break;

            }
            break;

        case RIGHT:
            switch ( control->cm_solver_pre_app_type )
            {
            case TRI_SOLVE_PA:
                switch ( control->cm_solver_pre_comp_type )
                {
                case JACOBI_PC:
                case SAI_PC:
                    if ( x != y )
                    {
                        Vector_Copy( x, y, workspace->H.n );
                    }
                    break;
                case ICHOLT_PC:
                case ILUT_PC:
                case ILUTP_PC:
                case FG_ILUT_PC:
                    tri_solve( &workspace->U, y, x, UPPER );
                    break;
                default:
                    fprintf( stderr, "[ERROR] Unrecognized preconditioner application method. Terminating...\n" );
                    exit( INVALID_INPUT );
                    break;
                }
                break;
            case TRI_SOLVE_LEVEL_SCHED_PA:
                switch ( control->cm_solver_pre_comp_type )
                {
                case JACOBI_PC:
                case SAI_PC:
                    if ( x != y )
                    {
                        Vector_Copy( x, y, workspace->H.n );
                    }
                    break;
                case ICHOLT_PC:
                case ILUT_PC:
                case ILUTP_PC:
                case FG_ILUT_PC:
                    tri_solve_level_sched( (static_storage *) workspace,
                            &workspace->U, y, x, UPPER, fresh_pre );
                    break;
                default:
                    fprintf( stderr, "[ERROR] Unrecognized preconditioner application method. Terminating...\n" );
                    exit( INVALID_INPUT );
                    break;
                }
                break;
            case TRI_SOLVE_GC_PA:
                switch ( control->cm_solver_pre_comp_type )
                {
                case JACOBI_PC:
                case SAI_PC:
                    fprintf( stderr, "[ERROR] Unsupported preconditioner computation/application method combination. Terminating...\n" );
                    exit( INVALID_INPUT );
                    break;
                case ICHOLT_PC:
                case ILUT_PC:
                case ILUTP_PC:
                case FG_ILUT_PC:
                    tri_solve_level_sched( (static_storage *) workspace,
                            &workspace->U, y, x, UPPER, fresh_pre );
                    permute_vector( workspace->x_p, x, workspace->permuted_row_col_inv, workspace->U.n );
                    Vector_Copy( x, workspace->x_p, workspace->U.n );
                    break;
                default:
                    fprintf( stderr, "[ERROR] Unrecognized preconditioner application method. Terminating...\n" );
                    exit( INVALID_INPUT );
                    break;
                }
                break;
            case JACOBI_ITER_PA:
                switch ( control->cm_solver_pre_comp_type )
                {
                case JACOBI_PC:
                case SAI_PC:
                    fprintf( stderr, "[ERROR] Unsupported preconditioner computation/application method combination. Terminating...\n" );
                    exit( INVALID_INPUT );
                    break;
                case ICHOLT_PC:
                case ILUT_PC:
                case ILUTP_PC:
                case FG_ILUT_PC:
                    /* construct D^{-1}_U */
                    if ( fresh_pre == TRUE )
                    {
#if defined(_OPENMP)
                        #pragma omp for schedule(static)
#endif
                        for ( i = 0; i < workspace->U.n; ++i )
                        {
                            si = workspace->U.start[i];
                            workspace->Dinv_U[i] = 1.0 / workspace->U.val[si];
                        }
                    }

                    jacobi_iter( workspace, &workspace->U, workspace->Dinv_U,
                            y, x, UPPER, control->cm_solver_pre_app_jacobi_iters );
                    break;
                default:
                    fprintf( stderr, "[ERROR] Unrecognized preconditioner application method. Terminating...\n" );
                    exit( INVALID_INPUT );
                    break;
                }
                break;
            default:
                fprintf( stderr, "[ERROR] Unrecognized preconditioner application method. Terminating...\n" );
                exit( INVALID_INPUT );
                break;

            }
            break;

        default:
            fprintf( stderr, "[ERROR] Unrecognized preconditioner application side. Terminating...\n" );
            exit( INVALID_INPUT );
            break;
        }
    }
}


/* Generalized minimual residual method with restarting and
 * left preconditioning for sparse linear systems */
int GMRES( const static_storage * const workspace, const control_params * const control,
        simulation_data * const data, const sparse_matrix * const H, const real * const b,
        const real tol, real * const x, const int fresh_pre )
{
    int i, j, k, itr, N, g_j, g_itr;
    real cc, tmp1, tmp2, temp, bnorm, g_bnorm;
    real t_start, t_ortho, t_pa, t_spmv, t_ts, t_vops;

    N = H->n;
    g_j = 0;
    g_itr = 0;
    t_ortho = 0.0;
    t_pa = 0.0;
    t_spmv = 0.0;
    t_ts = 0.0;
    t_vops = 0.0;

#if defined(_OPENMP)
    #pragma omp parallel default(none) \
    private(i, j, k, itr, bnorm, temp, t_start) \
    firstprivate(control, workspace, H, b, tol, x, fresh_pre) \
    shared(N, cc, tmp1, tmp2, g_itr, g_j, g_bnorm, stderr) \
    reduction(+: t_ortho, t_pa, t_spmv, t_ts, t_vops)
#endif
    {
        j = 0;
        itr = 0;
        t_ortho = 0.0;
        t_pa = 0.0;
        t_spmv = 0.0;
        t_ts = 0.0;
        t_vops = 0.0;

        t_start = Get_Time( );
        bnorm = Norm( b, N );
        t_vops += Get_Timing_Info( t_start );

        /* GMRES outer-loop */
        for ( itr = 0; itr < control->cm_solver_max_iters; ++itr )
        {
            /* calculate r0 */
            t_start = Get_Time( );
            sparse_matvec( workspace, H, x, workspace->b_prm );
            t_spmv += Get_Timing_Info( t_start );

            t_start = Get_Time( );
            Vector_Sum( workspace->b_prc, 1.0, b, -1.0, workspace->b_prm, N );
            t_vops += Get_Timing_Info( t_start );

            t_start = Get_Time( );
            apply_preconditioner( workspace, control, workspace->b_prc, workspace->b_prm,
                    itr == 0 ? fresh_pre : FALSE, LEFT );
            apply_preconditioner( workspace, control, workspace->b_prm, workspace->v[0],
                    itr == 0 ? fresh_pre : FALSE, RIGHT );
            t_pa += Get_Timing_Info( t_start );

            t_start = Get_Time( );
            temp = Norm( workspace->v[0], N );

#if defined(_OPENMP)
            #pragma omp single
#endif
            workspace->g[0] = temp;

            Vector_Scale( workspace->v[0], 1.0 / temp, workspace->v[0], N );
            t_vops += Get_Timing_Info( t_start );

            /* GMRES inner-loop */
            for ( j = 0; j < control->cm_solver_restart && FABS(workspace->g[j]) / bnorm > tol; j++ )
            {
                /* matvec */
                t_start = Get_Time( );
                sparse_matvec( workspace, H, workspace->v[j], workspace->b_prc );
                t_spmv += Get_Timing_Info( t_start );

                t_start = Get_Time( );
                apply_preconditioner( workspace, control, workspace->b_prc,
                        workspace->b_prm, FALSE, LEFT );
                apply_preconditioner( workspace, control, workspace->b_prm,
                        workspace->v[j + 1], FALSE, RIGHT );
                t_pa += Get_Timing_Info( t_start );

                /* apply modified Gram-Schmidt to orthogonalize the new residual */
                t_start = Get_Time( );
                for ( i = 0; i <= j; i++ )
                {
                    temp = Dot( workspace->v[i], workspace->v[j + 1], N );

#if defined(_OPENMP)
                    #pragma omp single
#endif
                    workspace->h[i][j] = temp;

                    Vector_Add( workspace->v[j + 1], -1.0 * temp, workspace->v[i], N );

                }
                t_vops += Get_Timing_Info( t_start );

                t_start = Get_Time( );
                temp = Norm( workspace->v[j + 1], N );

#if defined(_OPENMP)
                #pragma omp single
#endif
                workspace->h[j + 1][j] = temp;

                Vector_Scale( workspace->v[j + 1], 1.0 / temp, workspace->v[j + 1], N );
                t_vops += Get_Timing_Info( t_start );

                t_start = Get_Time( );
                /* Givens rotations on the upper-Hessenberg matrix to make it U */
#if defined(_OPENMP)
                #pragma omp single
#endif
                {
                    for ( i = 0; i <= j; i++ )
                    {
                        if ( i == j )
                        {
                            cc = SQRT( SQR(workspace->h[j][j]) + SQR(workspace->h[j + 1][j]) );
                            workspace->hc[j] = workspace->h[j][j] / cc;
                            workspace->hs[j] = workspace->h[j + 1][j] / cc;
                        }

                        tmp1 = workspace->hc[i] * workspace->h[i][j] +
                            workspace->hs[i] * workspace->h[i + 1][j];
                        tmp2 = -1.0 * workspace->hs[i] * workspace->h[i][j] +
                            workspace->hc[i] * workspace->h[i + 1][j];

                        workspace->h[i][j] = tmp1;
                        workspace->h[i + 1][j] = tmp2;
                    }

                    /* apply Givens rotations to the rhs as well */
                    tmp1 = workspace->hc[j] * workspace->g[j];
                    tmp2 = -1.0 * workspace->hs[j] * workspace->g[j];
                    workspace->g[j] = tmp1;
                    workspace->g[j + 1] = tmp2;

                }
                t_ortho += Get_Timing_Info( t_start );
            }

            /* solve Hy = g: H is now upper-triangular, do back-substitution */
            t_start = Get_Time( );
#if defined(_OPENMP)
            #pragma omp single
#endif
            {
                for ( i = j - 1; i >= 0; i-- )
                {
                    temp = workspace->g[i];
                    for ( k = j - 1; k > i; k-- )
                    {
                        temp -= workspace->h[i][k] * workspace->y[k];
                    }

                    workspace->y[i] = temp / workspace->h[i][i];
                }
            }
            t_ts += Get_Timing_Info( t_start );

            /* update x = x_0 + Vy */
            t_start = Get_Time( );
            Vector_MakeZero( workspace->p, N );

            for ( i = 0; i < j; i++ )
            {
                Vector_Add( workspace->p, workspace->y[i], workspace->v[i], N );
            }

            Vector_Add( x, 1.0, workspace->p, N );
            t_vops += Get_Timing_Info( t_start );

            /* stopping condition */
            if ( FABS(workspace->g[j]) / bnorm <= tol )
            {
                break;
            }
        }

#if defined(_OPENMP)
        #pragma omp single
#endif
        {
            g_j = j;
            g_itr = itr;
            g_bnorm = bnorm;
        }
    }

    data->timing.cm_solver_orthog += t_ortho / control->num_threads;
    data->timing.cm_solver_pre_app += t_pa / control->num_threads;
    data->timing.cm_solver_spmv += t_spmv / control->num_threads;
    data->timing.cm_solver_tri_solve += t_ts / control->num_threads;
    data->timing.cm_solver_vector_ops += t_vops / control->num_threads;

    if ( g_itr >= control->cm_solver_max_iters )
    {
        fprintf( stderr, "[WARNING] GMRES convergence failed (%d outer iters)\n", g_itr );
        fprintf( stderr, "  [INFO] Rel. residual error: %f\n", FABS(workspace->g[g_j]) / g_bnorm );
        return g_itr * (control->cm_solver_restart + 1) + g_j + 1;
    }

    return g_itr * (control->cm_solver_restart + 1) + g_j + 1;
}


int GMRES_HouseHolder( const static_storage * const workspace,
        const control_params * const control, simulation_data * const data,
        const sparse_matrix * const H, const real * const b, real tol,
        real * const x, const int fresh_pre )
{
    int i, j, k, itr, N, g_j, g_itr;
    real cc, tmp1, tmp2, temp, bnorm, g_bnorm;
    real v[10000], z[control->cm_solver_restart + 2][10000], w[control->cm_solver_restart + 2];
    real u[control->cm_solver_restart + 2][10000];
    real t_start, t_ortho, t_pa, t_spmv, t_ts, t_vops;

    j = 0;
    N = H->n;
    t_ortho = 0.0;
    t_pa = 0.0;
    t_spmv = 0.0;
    t_ts = 0.0;
    t_vops = 0.0;

#if defined(_OPENMP)
    #pragma omp parallel default(none) \
    private(i, j, k, itr, bnorm, temp, t_start) \
    firstprivate(control, workspace, H, b, tol, x, fresh_pre) \
    shared(v, z, w, u, N, cc, tmp1, tmp2, g_itr, g_j, g_bnorm, stderr) \
    reduction(+: t_ortho, t_pa, t_spmv, t_ts, t_vops)
#endif
    {
        j = 0;
        t_ortho = 0.0;
        t_pa = 0.0;
        t_spmv = 0.0;
        t_ts = 0.0;
        t_vops = 0.0;

        t_start = Get_Time( );
        bnorm = Norm( b, N );
        t_vops += Get_Timing_Info( t_start );

//        memset( x, 0, sizeof(real) * N );

        /* GMRES outer-loop */
        for ( itr = 0; itr < control->cm_solver_max_iters; ++itr )
        {
            /* compute z = r0 */
            t_start = Get_Time( );
            sparse_matvec( workspace, H, x, workspace->b_prm );
            t_spmv += Get_Timing_Info( t_start );

            t_start = Get_Time( );
            Vector_Sum( workspace->b_prc, 1.0,  b, -1.0, workspace->b_prm, N );
            t_vops += Get_Timing_Info( t_start );

            t_start = Get_Time( );
            apply_preconditioner( workspace, control, workspace->b_prc,
                    workspace->b_prm, fresh_pre, LEFT );
            apply_preconditioner( workspace, control, workspace->b_prm,
                    z[0], fresh_pre, RIGHT );
            t_pa += Get_Timing_Info( t_start );

            t_start = Get_Time( );
            Vector_MakeZero( w, control->cm_solver_restart + 1 );
            w[0] = Norm( z[0], N );

            Vector_Copy( u[0], z[0], N );
            u[0][0] += ( u[0][0] < 0.0 ? -1 : 1 ) * w[0];
            Vector_Scale( u[0], 1.0 / Norm( u[0], N ), u[0], N );

            w[0] *= ( u[0][0] < 0.0 ?  1 : -1 );
            t_vops += Get_Timing_Info( t_start );

            /* GMRES inner-loop */
            for ( j = 0; j < control->cm_solver_restart && FABS( w[j] ) / bnorm > tol; j++ )
            {
                /* compute v_j */
                t_start = Get_Time( );
                Vector_Scale( z[j], -2.0 * u[j][j], u[j], N );
                z[j][j] += 1.; /* due to e_j */

                for ( i = j - 1; i >= 0; --i )
                {
                    Vector_Add( z[j] + i, -2.0 * Dot( u[i] + i, z[j] + i, N - i ),
                            u[i] + i, N - i );
                }
                t_vops += Get_Timing_Info( t_start );

                /* matvec */
                t_start = Get_Time( );
                sparse_matvec( workspace, H, z[j], workspace->b_prc );
                t_spmv += Get_Timing_Info( t_start );

                t_start = Get_Time( );
                apply_preconditioner( workspace, control, workspace->b_prc,
                        workspace->b_prm, fresh_pre, LEFT );
                apply_preconditioner( workspace, control, workspace->b_prm,
                        v, fresh_pre, RIGHT );
                t_pa += Get_Timing_Info( t_start );

                t_start = Get_Time( );
                for ( i = 0; i <= j; ++i )
                {
                    Vector_Add( v + i, -2.0 * Dot( u[i] + i, v + i, N - i ),
                                u[i] + i, N - i );
                }

                if ( !Vector_isZero( v + (j + 1), N - (j + 1) ) )
                {
                    /* compute the HouseHolder unit vector u_j+1 */
                    Vector_MakeZero( u[j + 1], j + 1 );
                    Vector_Copy( u[j + 1] + (j + 1), v + (j + 1), N - (j + 1) );
                    temp = Norm( v + (j + 1), N - (j + 1) );
#if defined(_OPENMP)
                    #pragma omp single
#endif
                    u[j + 1][j + 1] += ( v[j + 1] < 0.0 ? -1.0 : 1.0 ) * temp;

#if defined(_OPENMP)
                    #pragma omp barrier
#endif

                    Vector_Scale( u[j + 1], 1.0 / Norm( u[j + 1], N ), u[j + 1], N );

                    /* overwrite v with P_m+1 * v */
                    temp = 2.0 * Dot( u[j + 1] + (j + 1), v + (j + 1), N - (j + 1) )
                           * u[j + 1][j + 1];
#if defined(_OPENMP)
                    #pragma omp single
#endif
                    v[j + 1] -= temp;

#if defined(_OPENMP)
                    #pragma omp barrier
#endif

                    Vector_MakeZero( v + (j + 2), N - (j + 2) );
//                    Vector_Add( v, -2.0 * Dot( u[j+1], v, N ), u[j+1], N );
                }
                t_vops += Get_Timing_Info( t_start );

                /* prev Givens rots on the upper-Hessenberg matrix to make it U */
                t_start = Get_Time( );
#if defined(_OPENMP)
                #pragma omp single
#endif
                {
                    for ( i = 0; i < j; i++ )
                    {
                        tmp1 =  workspace->hc[i] * v[i] + workspace->hs[i] * v[i + 1];
                        tmp2 = -workspace->hs[i] * v[i] + workspace->hc[i] * v[i + 1];

                        v[i] = tmp1;
                        v[i + 1] = tmp2;
                    }

                    /* apply the new Givens rotation to H and right-hand side */
                    if ( FABS(v[j + 1]) >= ALMOST_ZERO )
                    {
                        cc = SQRT( SQR( v[j] ) + SQR( v[j + 1] ) );
                        workspace->hc[j] = v[j] / cc;
                        workspace->hs[j] = v[j + 1] / cc;

                        tmp1 =  workspace->hc[j] * v[j] + workspace->hs[j] * v[j + 1];
                        tmp2 = -workspace->hs[j] * v[j] + workspace->hc[j] * v[j + 1];

                        v[j]   = tmp1;
                        v[j + 1] = tmp2;

                        /* Givens rotations to rhs */
                        tmp1 =  workspace->hc[j] * w[j];
                        tmp2 = -workspace->hs[j] * w[j];
                        w[j]   = tmp1;
                        w[j + 1] = tmp2;
                    }

                    /* extend R */
                    for ( i = 0; i <= j; ++i )
                    {
                        workspace->h[i][j] = v[i];
                    }
                }
                t_ortho += Get_Timing_Info( t_start );
            }


            /* solve Hy = w.
               H is now upper-triangular, do back-substitution */
            t_start = Get_Time( );
#if defined(_OPENMP)
            #pragma omp single
#endif
            {
                for ( i = j - 1; i >= 0; i-- )
                {
                    temp = w[i];
                    for ( k = j - 1; k > i; k-- )
                    {
                        temp -= workspace->h[i][k] * workspace->y[k];
                    }

                    workspace->y[i] = temp / workspace->h[i][i];
                }
            }
            t_ts += Get_Timing_Info( t_start );

            t_start = Get_Time( );
            for ( i = j - 1; i >= 0; i-- )
            {
                Vector_Add( x, workspace->y[i], z[i], N );
            }
            t_vops += Get_Timing_Info( t_start );

            /* stopping condition */
            if ( FABS( w[j] ) / bnorm <= tol )
            {
                break;
            }
        }

#if defined(_OPENMP)
        #pragma omp single
#endif
        {
            g_j = j;
            g_itr = itr;
            g_bnorm = bnorm;
        }
    }

    data->timing.cm_solver_orthog += t_ortho / control->num_threads;
    data->timing.cm_solver_pre_app += t_pa / control->num_threads;
    data->timing.cm_solver_spmv += t_spmv / control->num_threads;
    data->timing.cm_solver_tri_solve += t_ts / control->num_threads;
    data->timing.cm_solver_vector_ops += t_vops / control->num_threads;

    if ( g_itr >= control->cm_solver_max_iters )
    {
        fprintf( stderr, "[WARNING] GMRES convergence failed (%d outer iters)\n", g_itr );
        fprintf( stderr, "  [INFO] Rel. residual error: %f\n", FABS(w[g_j]) / g_bnorm );
        return g_itr * (control->cm_solver_restart + 1) + j + 1;
    }

    return g_itr * (control->cm_solver_restart + 1) + g_j + 1;
}


/* Conjugate Gradient */
int CG( const static_storage * const workspace, const control_params * const control,
        simulation_data * const data, const sparse_matrix * const H, const real * const b,
        const real tol, real * const x, const int fresh_pre )
{
    int i, g_itr, N;
    real tmp, alpha, beta, bnorm, g_bnorm, rnorm, g_rnorm;
    real *d, *r, *p, *z;
    real sig_old, sig_new;
    real t_start, t_pa, t_spmv, t_vops;

    N = H->n;
    d = workspace->d;
    r = workspace->r;
    p = workspace->q;
    z = workspace->p;
    t_pa = 0.0;
    t_spmv = 0.0;
    t_vops = 0.0;

#if defined(_OPENMP)
    #pragma omp parallel default(none) \
    private(i, tmp, alpha, beta, bnorm, rnorm, sig_old, sig_new, t_start) \
    firstprivate(control, workspace, H, b, tol, x, fresh_pre) \
    reduction(+: t_pa, t_spmv, t_vops) \
    shared(g_itr, g_bnorm, g_rnorm, N, d, r, p, z)
#endif
    {
        t_pa = 0.0;
        t_spmv = 0.0;
        t_vops = 0.0;

        t_start = Get_Time( );
        bnorm = Norm( b, N );
        t_vops += Get_Timing_Info( t_start );

        t_start = Get_Time( );
        sparse_matvec( workspace, H, x, d );
        t_spmv += Get_Timing_Info( t_start );

        t_start = Get_Time( );
        Vector_Sum( r, 1.0,  b, -1.0, d, N );
        rnorm = Norm( r, N );
        t_vops += Get_Timing_Info( t_start );

        t_start = Get_Time( );
        apply_preconditioner( workspace, control, r, d, fresh_pre, LEFT );
        apply_preconditioner( workspace, control, d, z, fresh_pre, RIGHT );
        t_pa += Get_Timing_Info( t_start );

        t_start = Get_Time( );
        Vector_Copy( p, z, N );
        sig_new = Dot( r, p, N );
        t_vops += Get_Timing_Info( t_start );

        for ( i = 0; i < control->cm_solver_max_iters && rnorm / bnorm > tol; ++i )
        {
            t_start = Get_Time( );
            sparse_matvec( workspace, H, p, d );
            t_spmv += Get_Timing_Info( t_start );

            t_start = Get_Time( );
            tmp = Dot( d, p, N );
            alpha = sig_new / tmp;
            Vector_Add( x, alpha, p, N );
            Vector_Add( r, -1.0 * alpha, d, N );
            rnorm = Norm( r, N );
            t_vops += Get_Timing_Info( t_start );

            t_start = Get_Time( );
            apply_preconditioner( workspace, control, r, d, FALSE, LEFT );
            apply_preconditioner( workspace, control, d, z, FALSE, RIGHT );
            t_pa += Get_Timing_Info( t_start );

            t_start = Get_Time( );
            sig_old = sig_new;
            sig_new = Dot( r, z, N );
            beta = sig_new / sig_old;
            Vector_Sum( p, 1.0, z, beta, p, N );
            t_vops += Get_Timing_Info( t_start );
        }

#if defined(_OPENMP)
        #pragma omp single
#endif
        {
            g_itr = i;
            g_bnorm = bnorm;
            g_rnorm = rnorm;
        }
    }

    data->timing.cm_solver_pre_app += t_pa / control->num_threads;
    data->timing.cm_solver_spmv += t_spmv / control->num_threads;
    data->timing.cm_solver_vector_ops += t_vops / control->num_threads;

    if ( g_itr >= control->cm_solver_max_iters )
    {
        fprintf( stderr, "[WARNING] CG convergence failed (%d iters)\n", g_itr );
        fprintf( stderr, "  [INFO] Rel. residual error: %f\n", g_rnorm / g_bnorm );
        return g_itr;
    }

    return g_itr;
}


/* Bi-conjugate gradient stabalized method with left preconditioning for
 * solving nonsymmetric linear systems
 *
 * workspace: struct containing storage for workspace for the linear solver
 * control: struct containing parameters governing the simulation and numeric methods
 * data: struct containing simulation data (e.g., atom info)
 * H: sparse, symmetric matrix, lower half stored in CSR format
 * b: right-hand side of the linear system
 * tol: tolerence compared against the relative residual for determining convergence
 * x: inital guess
 * fresh_pre: flag for determining if preconditioners should be recomputed
 *
 * Reference: Netlib (in MATLAB)
 *  http://www.netlib.org/templates/matlab/bicgstab.m
 * */
int BiCGStab( const static_storage * const workspace, const control_params * const control,
        simulation_data * const data, const sparse_matrix * const H, const real * const b,
        const real tol, real * const x, const int fresh_pre )
{
    int i, g_itr, N;
    real tmp, alpha, beta, omega, sigma, rho, rho_old, rnorm, g_rnorm, bnorm, g_bnorm, g_omega, g_rho;
    real t_start, t_pa, t_spmv, t_vops;

    N = H->n;
    t_pa = 0.0;
    t_spmv = 0.0;
    t_vops = 0.0;

#if defined(_OPENMP)
    #pragma omp parallel default(none) \
    private(i, tmp, alpha, beta, omega, sigma, rho, rho_old, rnorm, bnorm, t_start) \
    firstprivate(control, workspace, H, b, tol, x, fresh_pre) \
    reduction(+: t_pa, t_spmv, t_vops) \
    shared(g_itr, g_rnorm, g_bnorm, g_omega, g_rho, N)
#endif
    {
        t_pa = 0.0;
        t_spmv = 0.0;
        t_vops = 0.0;

        t_start = Get_Time( );
        sparse_matvec( workspace, H, x, workspace->d );
        t_spmv += Get_Timing_Info( t_start );

        t_start = Get_Time( );
        bnorm = Norm( b, N );
//        if ( FABS( bnorm ) < DBL_EPSILON )
        if ( bnorm == 0.0 )
        {
            bnorm = 1.0;
        }
        Vector_Sum( workspace->r, 1.0,  b, -1.0, workspace->d, N );
        t_vops += Get_Timing_Info( t_start );

        t_start = Get_Time( );
        Vector_Copy( workspace->r_hat, workspace->r, N );
        rnorm = Norm( workspace->r, N );
        omega = 1.0;
        rho = 1.0;
        t_vops += Get_Timing_Info( t_start );

        /* ensure each thread gets a local copy of
         * the function return value before proceeding
         * (Dot and Norm functions have persistent state in the form
         * of a shared global variable for the OpenMP version) */
#if defined(_OPENMP)
        #pragma omp barrier
#endif

        for ( i = 0; i < control->cm_solver_max_iters && rnorm / bnorm > tol; ++i )
        {
            t_start = Get_Time( );
            rho = Dot( workspace->r_hat, workspace->r, N );
//            if ( FABS( rho ) < DBL_EPSILON )
            if ( rho == 0.0 )
            {
                break;
            }
            if ( i > 0 )
            {
                #pragma GCC diagnostic push
                #pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
                beta = (rho / rho_old) * (alpha / omega);
                #pragma GCC diagnostic pop
                Vector_Sum( workspace->q, 1.0, workspace->p, -1.0 * omega, workspace->z, N );
                Vector_Sum( workspace->p, 1.0, workspace->r, beta, workspace->q, N );
            }
            else
            {
                Vector_Copy( workspace->p, workspace->r, N );
            }
            t_vops += Get_Timing_Info( t_start );

            t_start = Get_Time( );
            apply_preconditioner( workspace, control, workspace->p, workspace->y, fresh_pre, LEFT );
            apply_preconditioner( workspace, control, workspace->y, workspace->d, fresh_pre, RIGHT );
            t_pa += Get_Timing_Info( t_start );

            t_start = Get_Time( );
            sparse_matvec( workspace, H, workspace->d, workspace->z );
            t_spmv += Get_Timing_Info( t_start );

            t_start = Get_Time( );
            tmp = Dot( workspace->r_hat, workspace->z, N );
            alpha = rho / tmp;
            Vector_Sum( workspace->q, 1.0, workspace->r, -1.0 * alpha, workspace->z, N );
            tmp = Dot( workspace->q, workspace->q, N );
            /* early convergence check */
            if ( tmp < tol )
            {
                Vector_Add( x, alpha, workspace->d, N );
                break;
            }
            t_vops += Get_Timing_Info( t_start );

            t_start = Get_Time( );
            apply_preconditioner( workspace, control, workspace->q, workspace->y, fresh_pre, LEFT );
            apply_preconditioner( workspace, control, workspace->y, workspace->q_hat, fresh_pre, RIGHT );
            t_pa += Get_Timing_Info( t_start );

            t_start = Get_Time( );
            sparse_matvec( workspace, H, workspace->q_hat, workspace->y );
            t_spmv += Get_Timing_Info( t_start );

            t_start = Get_Time( );
            sigma = Dot( workspace->y, workspace->q, N );

            /* ensure each thread gets a local copy of
             * the function return value before proceeding
             * (Dot and Norm functions have persistent state in the form
             * of a shared global variable for the OpenMP version) */
#if defined(_OPENMP)
            #pragma omp barrier
#endif

            tmp = Dot( workspace->y, workspace->y, N );
            omega = sigma / tmp;
            Vector_Sum( workspace->g, alpha, workspace->d, omega, workspace->q_hat, N );
            Vector_Add( x, 1.0, workspace->g, N );
            Vector_Sum( workspace->r, 1.0, workspace->q, -1.0 * omega, workspace->y, N );
            rnorm = Norm( workspace->r, N );
//            if ( FABS( omega ) < DBL_EPSILON )
            if ( omega == 0.0 )
            {
                break;
            }
            rho_old = rho;
            t_vops += Get_Timing_Info( t_start );

            /* ensure each thread gets a local copy of
             * the function return value before proceeding
             * (Dot and Norm functions have persistent state in the form
             * of a shared global variable for the OpenMP version) */
#if defined(_OPENMP)
            #pragma omp barrier
#endif
        }

#if defined(_OPENMP)
        #pragma omp single
#endif
        {
            g_itr = i;
            g_rnorm = rnorm;
            g_bnorm = bnorm;
            g_omega = omega;
            g_rho = rho;
        }
    }

    data->timing.cm_solver_pre_app += t_pa / control->num_threads;
    data->timing.cm_solver_spmv += t_spmv / control->num_threads;
    data->timing.cm_solver_vector_ops += t_vops / control->num_threads;

//    if ( FABS( g_omega ) < DBL_EPSILON )
    if ( g_omega == 0.0 )
    {
        fprintf( stderr, "[WARNING] BiCGStab numeric breakdown (%d iters)\n", g_itr );
        fprintf( stderr, "  [INFO] omega = %f\n", g_omega );
    }
//    else if ( FABS( g_rho ) < DBL_EPSILON )
    else if ( g_rho == 0.0 )
    {
        fprintf( stderr, "[WARNING] BiCGStab numeric breakdown (%d iters)\n", g_itr );
        fprintf( stderr, "  [INFO] rho = %f\n", g_rho );
    }
    else if ( g_itr >= control->cm_solver_max_iters )
    {
        fprintf( stderr, "[WARNING] BiCGStab convergence failed (%d iters)\n", g_itr );
        fprintf( stderr, "  [INFO] Rel. residual error: %f\n", g_rnorm / g_bnorm );
    }

    return g_itr;
}


/* Steepest Descent */
int SDM( const static_storage * const workspace, const control_params * const control,
         simulation_data * const data, const sparse_matrix * const H, const real * const b,
         const real tol, real * const x, const int fresh_pre )
{
    int i, g_itr, N;
    real tmp, alpha, bnorm, g_bnorm;
    real sig, g_sig;
    real t_start, t_pa, t_spmv, t_vops;

    N = H->n;
    t_pa = 0.0;
    t_spmv = 0.0;
    t_vops = 0.0;

#if defined(_OPENMP)
    #pragma omp parallel default(none) \
    private(i, tmp, alpha, bnorm, sig, t_start) \
    firstprivate(control, workspace, H, b, tol, x, fresh_pre) \
    reduction(+: t_pa, t_spmv, t_vops) \
    shared(g_itr, g_sig, g_bnorm, N)
#endif
    {
        t_pa = 0.0;
        t_spmv = 0.0;
        t_vops = 0.0;

        t_start = Get_Time( );
        bnorm = Norm( b, N );
        t_vops += Get_Timing_Info( t_start );

        t_start = Get_Time( );
        sparse_matvec( workspace, H, x, workspace->q );
        t_spmv += Get_Timing_Info( t_start );

        t_start = Get_Time( );
        Vector_Sum( workspace->r, 1.0,  b, -1.0, workspace->q, N );
        t_vops += Get_Timing_Info( t_start );

        t_start = Get_Time( );
        apply_preconditioner( workspace, control, workspace->r, workspace->q, fresh_pre, LEFT );
        apply_preconditioner( workspace, control, workspace->q, workspace->d, fresh_pre, RIGHT );
        t_pa += Get_Timing_Info( t_start );

        t_start = Get_Time( );
        sig = Dot( workspace->r, workspace->d, N );
        t_vops += Get_Timing_Info( t_start );

        for ( i = 0; i < control->cm_solver_max_iters && SQRT(sig) / bnorm > tol; ++i )
        {
            t_start = Get_Time( );
            sparse_matvec( workspace, H, workspace->d, workspace->q );
            t_spmv += Get_Timing_Info( t_start );

            t_start = Get_Time( );
            sig = Dot( workspace->r, workspace->d, N );

            /* ensure each thread gets a local copy of
             * the function return value before proceeding
             * (Dot function has persistent state in the form
             * of a shared global variable for the OpenMP version) */
#if defined(_OPENMP)
            #pragma omp barrier
#endif

            tmp = Dot( workspace->d, workspace->q, N );
            alpha = sig / tmp;

            Vector_Add( x, alpha, workspace->d, N );
            Vector_Add( workspace->r, -alpha, workspace->q, N );
            t_vops += Get_Timing_Info( t_start );

            t_start = Get_Time( );
            apply_preconditioner( workspace, control, workspace->r, workspace->q, FALSE, LEFT );
            apply_preconditioner( workspace, control, workspace->q, workspace->d, FALSE, RIGHT );
            t_pa += Get_Timing_Info( t_start );
        }

#if defined(_OPENMP)
        #pragma omp single
#endif
        {
            g_itr = i;
            g_sig = sig;
            g_bnorm = bnorm;
        }
    }

    data->timing.cm_solver_pre_app += t_pa / control->num_threads;
    data->timing.cm_solver_spmv += t_spmv / control->num_threads;
    data->timing.cm_solver_vector_ops += t_vops / control->num_threads;

    if ( g_itr >= control->cm_solver_max_iters  )
    {
        fprintf( stderr, "[WARNING] SDM convergence failed (%d iters)\n", g_itr );
        fprintf( stderr, "  [INFO] Rel. residual error: %f\n", SQRT(g_sig) / g_bnorm );
        return g_itr;
    }

    return g_itr;
}


/* Estimate the stability of a 2-side preconditioning scheme
 * using the factorization A \approx LU. Specifically, estimate the 1-norm of A^{-1}
 * using the 1-norm of (LU)^{-1}e, with e = [1 1 ... 1]^T through 2 triangular solves:
 *   1) Ly = e
 *   2) Ux = y where y = Ux
 * That is, we seek to solve e = LUx for unknown x
 *
 * Reference: Incomplete LU Preconditioning with the Multilevel Fast Multipole Algorithm
 *   for Electromagnetic Scattering, SIAM J. Sci. Computing, 2007 */
real condest( const sparse_matrix * const L, const sparse_matrix * const U )
{
    unsigned int i, N;
    real *e, c;

    N = L->n;

    e = smalloc( sizeof(real) * N, __FILE__, __LINE__ );

    for ( i = 0; i < N; ++i )
        e[i] = 1.0;

    tri_solve( L, e, e, LOWER );
    tri_solve( U, e, e, UPPER );

    /* compute 1-norm of vector e */
    c = FABS(e[0]);
    for ( i = 1; i < N; ++i)
    {
        if ( FABS(e[i]) > c )
        {
            c = FABS(e[i]);
        }

    }

    sfree( e, __FILE__, __LINE__ );

    return c;
}
