#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "quern.h"

int QUERN_convert_column_format_to_row_format(int num_rows,
                                              int num_columns,
                                              const int* A0_column_start,
                                              const int* A0_row_index,
                                              const double* A0_value,
                                              int* A1_row_start,
                                              int* A1_column_index,
                                              double* A1_value)
{
   // check input
   if(num_rows<=0 || num_columns<=0) return QUERN_INPUT_ERROR;
   if(!A0_column_start || !A0_row_index || !A0_value) return QUERN_INPUT_ERROR;
   if(!A1_row_start || !A1_column_index || !A1_value) return QUERN_INPUT_ERROR;
   // figure out number of entries in each row
   std::memset(A1_row_start, 0, (num_rows+1)*sizeof(int));
   for(int i=0; i<num_columns; ++i){
      if(A0_column_start[i]>A0_column_start[i+1]) return QUERN_INPUT_ERROR;
      for(int j=A0_column_start[i]; j<A0_column_start[i+1]; ++j){
         if(A0_row_index[j]<0 || A0_row_index[j]>=num_rows)
            return QUERN_INPUT_ERROR;
         ++A1_row_start[A0_row_index[j]+1];
      }
   }
   // cumulative sum to get row_start
   for(int i=0; i<num_rows; ++i)
      A1_row_start[i+1]+=A1_row_start[i];
   // use a temporary copy of row_start for keeping track of where we add
   int* row_pointer=(int*)std::malloc(num_rows*sizeof(int));
   if(!row_pointer) return QUERN_OUT_OF_MEMORY;
   std::memcpy(row_pointer, A1_row_start, num_rows*sizeof(int));
   // then fill in the entries
   for(int i=0; i<num_columns; ++i){
      for(int j=A0_column_start[i]; j<A0_column_start[i+1]; ++j){
         int r=A0_row_index[j];
         A1_column_index[row_pointer[r]]=i;
         A1_value[row_pointer[r]]=A0_value[j];
         ++row_pointer[r];
      }
   }
   std::free(row_pointer);
   return QUERN_OK;
}

int QUERN_multiply(int m,
                   int n,
                   const int* row_start,
                   const int* column_index,
                   const double* value,
                   const double* input,
                   double* result)
{
   if(m<=0 || n<=0 || !row_start || !column_index || !value
         || !input || !result)
      return QUERN_INPUT_ERROR;
   int i, j;
   double x;
   for(i=0; i<m; ++i){
      x=0;
      for(j=row_start[i]; j<row_start[i+1]; ++j)
         x+=value[j]*input[column_index[j]];
      result[i]=x;
   }
   return QUERN_OK;
}

int QUERN_multiply_transpose(int m,
                             int n,
                             const int* row_start,
                             const int* column_index,
                             const double* value,
                             const double* input,
                             double* result)
{
   if(m<=0 || n<=0 || !row_start || !column_index || !value
         || !input || !result)
      return QUERN_INPUT_ERROR;
   int i, j;
   double x;
   std::memset(result, 0, n*sizeof(double));
   for(i=0; i<m; ++i){
      x=input[i];
      for(j=row_start[i]; j<row_start[i+1]; ++j)
         result[column_index[j]]+=value[j]*x;
   }
   return QUERN_OK;
}

