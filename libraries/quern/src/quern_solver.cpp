#include <cmath>
#include <cstdlib>
#include <cstring>
#include "quern.h"

// TO-DO: use BLAS versions of these smaller kernels, if available

double two_norm(int n, const double* x)
{
   double r=0;
   for(int i=0; i<n; ++i) r+=x[i]*x[i];
   return std::sqrt(r);
}

double two_norm_squared(int n, const double* x)
{
   double r=0;
   for(int i=0; i<n; ++i) r+=x[i]*x[i];
   return r;
}

// x=x+alpha*y
void add_scaled(int n, double* x, double alpha, const double* y)
{
   for(int i=0; i<n; ++i) x[i]+=alpha*y[i];
}

// x=beta*x+y
void scale_and_add(int n, double beta, double* x, const double* y)
{
   for(int i=0; i<n; ++i) x[i]=beta*x[i]+y[i];
}

int QUERN_solve_with_CGNR(int m,
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
                          double* return_residual_norm)
{
   if(m<=0 || n<=0 || !A_row_start || !A_column_index || !A_value
      || !rhs || !R_row_start || !R_column_index || !R_value || !x
      || !return_solved || !return_iterations || !return_residual_norm)
      return QUERN_INPUT_ERROR;
   // default values
   *return_solved=0;
   *return_iterations=0;
   *return_residual_norm=two_norm(n, rhs);
   if(*return_residual_norm<=absolute_convergence_tolerance){
      *return_solved=1;
      return QUERN_OK;
   }
   // allocate some room to work in
   double* working_vectors=(double*)std::malloc((3*n+m)*sizeof(double));
   if(!working_vectors)
      return QUERN_OUT_OF_MEMORY;
   double* r=working_vectors;
   double* s=r+n;
   double* z=s+n;
   double* u=z+n;
   // set up CGNR
   int check; 
   std::memset(x, 0, n*sizeof(double));
   std::memcpy(r, rhs, n*sizeof(double));
   std::memcpy(u, rhs, n*sizeof(double));
   check=QUERN_solve_with_r_transpose_in_place(n, R_row_start, R_column_index,
                                               R_value, u);
   if(check){ std::free(working_vectors); return check; }
   check=QUERN_solve_with_r(n, R_row_start, R_column_index, R_value, u, z);
   if(check){ std::free(working_vectors); return check; }
   std::memcpy(s, z, n*sizeof(double));
   double rho=two_norm_squared(n, u);
   // the main loop
   for(;;){
      if(rho==0){ std::free(working_vectors); return QUERN_INPUT_ERROR; }
      check=QUERN_multiply(m, n, A_row_start, A_column_index, A_value, s, u);
      if(check){ std::free(working_vectors); return check; }
      check=QUERN_multiply_transpose(m, n, A_row_start, A_column_index, A_value,
                                     u, z);
      if(check){ std::free(working_vectors); return check; }
      double denom=two_norm_squared(m, u);
      if(denom==0){ std::free(working_vectors); return QUERN_INPUT_ERROR; }
      double alpha=rho/denom;
      add_scaled(n, x, alpha, s);
      add_scaled(n, r, -alpha, z);
      ++*return_iterations;
      *return_residual_norm=two_norm(n, r);
      if(*return_residual_norm<=absolute_convergence_tolerance){
         *return_solved=1;
         break;
      }
      if(*return_iterations>max_iterations)
         break;
      std::memcpy(u, r, n*sizeof(double));
      check=QUERN_solve_with_r_transpose_in_place(n, R_row_start,
                                                  R_column_index, R_value, u);
      if(check){ std::free(working_vectors); return check; }
      check=QUERN_solve_with_r(n, R_row_start, R_column_index, R_value, u, z);
      if(check){ std::free(working_vectors); return check; }
      double rho_new=two_norm_squared(n, u);
      double beta=rho_new/rho;
      scale_and_add(n, beta, s, z);
      rho=rho_new;
   }
   std::free(working_vectors);
   return QUERN_OK;
}
