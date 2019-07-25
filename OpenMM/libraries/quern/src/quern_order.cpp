#include <cassert>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include "quern.h"

int QUERN_get_rbfs_column_ordering(int m,
                                   int n,
                                   const int* A_row_start,
                                   const int* A_column_index,
                                   int* column_order)
{
   if(m<=0 || n<=0 || !A_row_start || !A_column_index || !column_order)
      return QUERN_INPUT_ERROR;
   // get some memory to work in
   int* work=(int*)std::malloc((n+1+A_row_start[m]+n+n)*sizeof(int));
   if(!work) return QUERN_OUT_OF_MEMORY;
   // figure out number of entries in each column
   int* column_start=work;
   std::memset(column_start, 0, (n+1)*sizeof(int));
   for(int i=0; i<m; ++i){
      for(int j=A_row_start[i]; j<A_row_start[i+1]; ++j)
         ++column_start[A_column_index[j]+1];
   }
   // cumulative sum to get column_start
   for(int i=0; i<n; ++i)
      column_start[i+1]+=column_start[i];
   assert(column_start[n]==A_row_start[m]);
   // list the columns now
   int* row_index=column_start+(n+1);
   int* column_pointer=row_index+column_start[n];
   std::memcpy(column_pointer, column_start, n*sizeof(int));
   for(int i=0; i<m; ++i){
      for(int j=A_row_start[i]; j<A_row_start[i+1]; ++j){
         int c=A_column_index[j];
         row_index[column_pointer[c]++]=i;
      }
   }
   // set up marker for BFS
   char* column_marker=(char*)(column_pointer+n);
   std::memset(column_marker, 0, n);
   // and do as many BFS as we need to hit all connected components
   int p=n;
   for(int root=0; root<n; ++root) if(!column_marker[root]){
      column_order[--p]=root;
      column_marker[root]=1;
      for(int i=p; i>=p; --i){
         int j=column_order[i];
         // add unmarked neighbour columns of j to ordering
         for(int k=column_start[j]; k<column_start[j+1]; ++k){
            int r=row_index[k];
            for(int a=A_row_start[r]; a<A_row_start[r+1]; ++a){
               int nbr=A_column_index[a];
               if(!column_marker[nbr]){
                  column_order[--p]=nbr;
                  column_marker[nbr]=1;
               }
            }
         }
      }
   }
   assert(p==0);
   std::free(work);
   return QUERN_OK;
}

int QUERN_reorder_columns(int m,
                          int n,
                          const int* column_order,
                          const int* A_row_start,
                          int* A_column_index,
                          double* A_value)
{
   if(m<=0 || n<=0 || !column_order
         || !A_row_start || !A_column_index || !A_value)
      return QUERN_INPUT_ERROR;

   // Since we don't expect to have lots of empty columns, it would
   // probably be superior to do this as two transposes --- but we'll
   // worry about that optimization later.

   // get some temporary storage for permuting and sorting
   std::pair<int,double>* row=(std::pair<int,double>*)
                                  std::malloc(n*sizeof(std::pair<int,double>));
   if(!row) return QUERN_OUT_OF_MEMORY;
   // and invert the permutation
   int* inv=(int*)std::malloc(n*sizeof(int));
   if(!inv){
      std::free(row);
      return QUERN_OUT_OF_MEMORY;
   }
   for(int i=0; i<n; ++i) inv[column_order[i]]=i;
   // do it row by row
   int k;
   for(int i=0; i<m; ++i){
      k=0;
      for(int j=A_row_start[i]; j<A_row_start[i+1]; ++j)
         row[k++]=std::make_pair(inv[A_column_index[j]], A_value[j]);
      std::sort(row, row+k); // at some point should replace with radix-sort
      k=0;
      for(int j=A_row_start[i]; j<A_row_start[i+1]; ++j){
         A_column_index[j]=row[k].first;
         A_value[j]=row[k].second;
         ++k;
      } 
   }
   std::free(row);
   std::free(inv);
   return QUERN_OK;
}

int QUERN_get_profile_row_ordering(int m,
                                   int n,
                                   const int* A_row_start,
                                   const int* A_column_index,
                                   int* row_order)
{
   if(m<=0 || n<=0 || !A_row_start || !A_column_index || !row_order)
      return QUERN_INPUT_ERROR;
   // first count how many rows start at column i for each i=0, .., n-1
   int* count=(int*)std::calloc(n+1, sizeof(int));
   if(!count) return QUERN_OUT_OF_MEMORY;
   for(int i=0; i<m; ++i)
      if(A_row_start[i]<A_row_start[i+1])
         ++count[A_column_index[A_row_start[i]]+1];
   // then do a cumulative sum to find target locations of the rows
   for(int j=2; j<n+1; ++j)
      count[j]+=count[j-1];
   // and now write out the row ordering with a second pass
   for(int i=0; i<m; ++i){
      if(A_row_start[i]<A_row_start[i+1])
         row_order[count[A_column_index[A_row_start[i]]]++]=i;
      else
         row_order[count[n]++]=i;
   }
   std::free(count);
   return QUERN_OK;
}

// Reorders a given input vector x into an output vector y of length m.
int QUERN_reorder_vector(int m,
                         const int* order,
                         const double* x,
                         double* y)
{
   if(m<=0 || !order || !x || !y) return QUERN_INPUT_ERROR;
   for(int i=0; i<m; ++i)
      y[i]=x[order[i]];
   return QUERN_OK;
}

// Applies the inverse order to a given input vector x, saving in an output
// vector y of length m.
int QUERN_inverse_order_vector(int m,
                               const int* order,
                               const double* x,
                               double* y)
{
   if(m<=0 || !order || !x || !y) return QUERN_INPUT_ERROR;
   for(int i=0; i<m; ++i)
      y[order[i]]=x[i];
   return QUERN_OK;
}
