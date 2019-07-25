/*
* Quern: a sparse QR library.
* This code is in the public domain.
* - Robert Bridson
*/

#ifndef QUERN_H
#define QUERN_H

#ifdef __cplusplus
extern "C" {
#endif

#include "openmm/internal/windowsExport.h"

// Function return values: nonzero indicates an error
#define QUERN_OK 0
#define QUERN_INPUT_ERROR 1
#define QUERN_OUT_OF_MEMORY 2

// Takes a compressed sparse column format matrix and outputs it in compressed
// sparse row format (suitable for the QR factorization routines).
// Requires A1_row_start to have room for num_columns+1 entries,
//          A1_column_index to have room for all nonzeros,
//      and A1_value to have room for all nonzeros.
int OPENMM_EXPORT QUERN_convert_column_format_to_row_format(int num_rows,
                                              int num_columns,
                                              const int* A0_column_start,
                                              const int* A0_row_index,
                                              const double* A0_value,
                                              int* A1_row_start,
                                              int* A1_column_index,
                                              double* A1_value);

// Takes a compressed sparse row format m*n matrix and finds a reversed
// breadth-first search column ordering, suitable for incomplete QR
// preconditioning.
// Requires column_order to have room for n entries.
int OPENMM_EXPORT QUERN_get_rbfs_column_ordering(int m,
                                   int n,
                                   const int* A_row_start,
                                   const int* A_column_index,
                                   int* column_order);

// Takes a length n column ordering and a compressed sparse row format m*n
// matrix, then reorders each row accordingly.
int OPENMM_EXPORT QUERN_reorder_columns(int m,
                          int n,
                          const int* column_order,
                          const int* A_row_start,
                          int* A_column_index,
                          double* A_value);

// Takes a compressed sparse row format m*n matrix and finds a row ordering
// which minimizes the lower profile---in particular, if A's rows can be
// permuted to make it upper triangular, row_order will do the trick.
// This is recommended as a cheap pre-process before QR, at least if nothing
// is known about the current structure, to improve performance.
// Requires row_order to have room for m integers.
int OPENMM_EXPORT QUERN_get_profile_row_ordering(int m,
                                   int n,
                                   const int* A_row_start,
                                   const int* A_column_index,
                                   int* row_order);

// Reorders a given input vector x into an output vector y of length m.
int OPENMM_EXPORT QUERN_reorder_vector(int m,
                         const int* order,
                         const double* x,
                         double* y);

// Applies the inverse order to a given input vector x, saving in an output
// vector y of length m.
int OPENMM_EXPORT QUERN_inverse_order_vector(int m,
                               const int* order,
                               const double* x,
                               double* y);

// Takes a compressed sparse row format m*n matrix (m>=n) and outputs the upper
// triangular R factor from QR, also in compressed sparse row format. If
// row_order is non-null, it should contain the order in which rows of A should be
// taken; if it is null, the natural ordering is used. Memory for R is
// allocated internally; use QUERN_free_result to free.
int OPENMM_EXPORT QUERN_compute_qr_without_q(int m,
                               int n,
                               const int* A_row_start,
                               const int* A_column_index,
                               const double* A_value,
                               const int* row_order,
                               int** ptr_R_row_start,
                               int** ptr_R_column_index,
                               double** ptr_R_value);

// Takes a compressed sparse row format m*n matrix (m>=n) and outputs the QR
// factors. The storage for Q encodes a sequence of Givens rotations and row
// swaps; it is not directly usable as a standard sparse matrix. R, however,
// is stored in standard compressed sparse row format. If row_order is non-null
// it should contain the order in which rows of A should be taken; if it is
// null, the natural ordering is used. Memory for Q and R is allocated
// internally; use QUERN_free_result to free.
int OPENMM_EXPORT QUERN_compute_qr(int m,
                     int n,
                     const int* A_row_start,
                     const int* A_column_index,
                     const double* A_value,
                     const int* row_order,
                     int** ptr_Q_row_start,
                     int** ptr_Q_column_index,
                     double** ptr_Q_value,
                     int** ptr_R_row_start,
                     int** ptr_R_column_index,
                     double** ptr_R_value);

// Takes a compressed sparse row format m*n matrix (m>=n) and outputs the upper
// triangular R factor from QR, with small entries dropped, also in compressed
// sparse row format. If row_order is non-null, it should contain the order in
// which rows of A should be taken; if it is null, the natural ordering is used.
// Memory for R is allocated internally; use QUERN_free_result to free.
int OPENMM_EXPORT QUERN_compute_incomplete_qr_without_q(int m,
                                          int n,
                                          const int* A_row_start,
                                          const int* A_column_index,
                                          const double* A_value,
                                          const int* row_order,
                                          double drop_tolerance,
                                          int** ptr_R_row_start,
                                          int** ptr_R_column_index,
                                          double** ptr_R_value);

// Takes a compressed sparse row format m*n matrix (m>=n) and outputs incomplete
// QR factors. The storage for Q encodes a sequence of Givens rotations and row
// swaps; it is not directly usable as a standard sparse matrix. R, however,
// is stored in standard compressed sparse row format. If row_order is non-null
// it should contain the order in which rows of A should be taken; if it is
// null, the natural ordering is used. Memory for Q and R is allocated
// internally; use QUERN_free_result to free.
int OPENMM_EXPORT QUERN_compute_incomplete_qr(int m,
                                int n,
                                const int* A_row_start,
                                const int* A_column_index,
                                const double* A_value,
                                const int* row_order,
                                double drop_tolerance,
                                int** ptr_Q_row_start,
                                int** ptr_Q_column_index,
                                double** ptr_Q_value,
                                int** ptr_R_row_start,
                                int** ptr_R_column_index,
                                double** ptr_R_value);

// Free the memory allocated during QR factorization (for either Q or R).
// After calling this, do not try to access the factor again!
void OPENMM_EXPORT QUERN_free_result(int* row_start,
                       int* column_index,
                       double* value);

// Compute the product of an m*n CSR matrix with n-vector input in
// m-vector result.
int OPENMM_EXPORT QUERN_multiply(int m,
                   int n,
                   const int* row_start,
                   const int* column_index,
                   const double* value,
                   const double* input,
                   double* result);

// Compute the product of an m*n CSR matrix transposes with m-vector input in
// n-vector result.
int OPENMM_EXPORT QUERN_multiply_transpose(int m,
                             int n,
                             const int* row_start,
                             const int* column_index,
                             const double* value,
                             const double* input,
                             double* result);

// Compute the product Q*x in-place (Q from QR factorization).
// Note: if the input is naturally a length n vector, you should pad it with
// zeros to make it length m.
// If a row reordering of A was involved, you probably want to
// call QUERN_inverse_order_vector afterwards.
int OPENMM_EXPORT QUERN_multiply_with_q(int m,
                          const int* Q_row_start,
                          const int* Q_column_index,
                          const double* Q_value,
                          double* x);

// Compute the product Q^T*x in-place (Q from QR factorization).
// If a row reordering of A was involved, you probably want to
// call QUERN_reorder_vector first.
int OPENMM_EXPORT QUERN_multiply_with_q_transpose(int m,
                                    const int* Q_row_start,
                                    const int* Q_column_index,
                                    const double* Q_value,
                                    double* x);

// Compute result=R^{-1}*rhs (R from QR factorization).
// The vector result may actually be aliased to rhs, for an in-place solve.
int OPENMM_EXPORT QUERN_solve_with_r(int n,
                       const int* R_row_start,
                       const int* R_column_index,
                       const double* R_value,
                       const double* rhs,
                       double* result);

// Compute R^{-T}*x in-place (R from QR factorization).
int OPENMM_EXPORT QUERN_solve_with_r_transpose_in_place(int n,
                                          const int* R_row_start,
                                          const int* R_column_index,
                                          const double* R_value,
                                          double* x);

// CGNR solver for A^T*A*x=rhs, starting with zero initial guess,
// with R as preconditioner.
int OPENMM_EXPORT QUERN_solve_with_CGNR(int m,
                          int n,
                          const int* A_row_start,
                          const int* A_column_index,
                          const double* A_value,
                          const double* rhs,
                          const int* R_row_start,
                          const int* R_column_index,
                          const double* R_value,
                          int max_iterations,
                          double absolute_convergence_tolerance,
                          double* x,
                          int* return_solved,
                          int* return_iterations,
                          double* return_residual_norm);

#ifdef __cplusplus
}
#endif

#endif
