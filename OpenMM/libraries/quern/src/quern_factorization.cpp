#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include "quern.h"
#include "quern_list.h"

#ifdef _WIN32
#define copysign _copysign
#endif

struct SparseEntry
{
   int index;
   double value;
   SparseEntry(void) {}
   SparseEntry(int index_, double value_) : index(index_), value(value_) {}
};

typedef Quern::List<SparseEntry> SparseVector;

static bool copy_row(int nnz,
                     const int* index,
                     const double* value,
                     SparseVector& x)
{
   for(int i=nnz-1; i>=0; --i) if(value[i]){
      if(!x.push_front(SparseEntry(index[i], value[i]))) return false;
   }
   return true;
}

// compute c>=0 and s, with c^2+s^2=1, so that s*a+c*b=0
static void givens(double a,
                   double b,
                   double& c,
                   double& s)
{
   if(b==0){
      c=1;
      s=0;
   }else{
      if(std::fabs(b)>std::fabs(a)){
         double tau=-a/b;
         s=1/std::sqrt(1+tau*tau);
         c=s*tau;
         if(c<0){ c=-c; s=-s; }
      }else{
         double tau=-b/a;
         c=1/std::sqrt(1+tau*tau);
         s=c*tau;
      }
   }
}

// given c and s (c>=0, c^2+s^2=1) encode them in a single double
static double encode(double c,
                     double s)
{
   if(std::fabs(s)<c)
      return s; // will be at most sqrt(0.5)=0.707... in magnitude
   else
      return copysign(1/c, s); // will be least sqrt(2)=1.414... in magnitude
}

// given output from encode, reconstruct c and s (c>=0, c^2+s^2=1)
static void decode(double tau,
                   double& c,
                   double& s)
{
   if(std::fabs(tau)<1){
      s=tau;
      c=std::sqrt(1-s*s);
   }else{
      c=1/tau;
      s=std::sqrt(1-c*c);
      if(c<0){ c=-c; s=-s; }
   }
}

static bool apply_givens(SparseVector& x,
                         SparseVector& y,
                         int diagonal,
                         double& c,
                         double& s)
{
   assert(!x.empty() && x.front().index==diagonal && x.front().value);
   assert(!y.empty() && y.front().index==diagonal && y.front().value);
   // find the rotation we need
   double a=x.front().value, b=y.front().value;
   givens(a, b, c, s);
   assert(c && s);
   // rotate the start of each list
   x.front().value=c*a-s*b;
   y.pop_front();
   // then update the rest (x_new = c*x-s*y and y_new = s*x+c*y)
   Quern::ListIterator<SparseEntry> p=x.begin(), q=y.begin();
   ++p; // skip the first value we already took care of
   while(p.still_going() && q.still_going()){
      if(p->index==q->index){
         double xnew=c*p->value-s*q->value,
                ynew=s*p->value+c*q->value;
         if(xnew){
            p->value=xnew;
            ++p;
         }else x.erase(p);
         if(ynew){
            q->value=ynew;
            ++q;
         }else y.erase(q);
      }else if(p->index<q->index){
         int k=p->index;
         double xnew=c*p->value,
                ynew=s*p->value;
         p->value=xnew;
         ++p;
         if(!y.insert(q, SparseEntry(k, ynew))) return false;
         ++q;
      }else{
         int k=q->index;
         double xnew=-s*q->value,
                ynew=c*q->value;
         if(!x.insert(p, SparseEntry(k, xnew))) return false;
         ++p;
         q->value=ynew;
         ++q;
      }
   }
   if(p.still_going()){
      do{
         int k=p->index;
         double xnew=c*p->value,
                ynew=s*p->value;
         p->value=xnew;
         ++p;
         if(!y.insert(q, SparseEntry(k, ynew))) return false;
         ++q;
      }while(p.still_going());
   }else if(q.still_going()){
      do{
         int k=q->index;
         double xnew=-s*q->value,
                ynew=c*q->value;
         if(!x.insert(p, SparseEntry(k, xnew))) return false;
         ++p;
         q->value=ynew;
         ++q;
      }while(q.still_going());
   }
   return true;
}

int QUERN_compute_qr_without_q(int m,
                               int n,
                               const int* A_row_start,
                               const int* A_column_index,
                               const double* A_value,
                               const int* row_order,
                               int** ptr_R_row_start,
                               int** ptr_R_column_index,
                               double** ptr_R_value)
{
   if(m<=0 || n<=0 || m<n || !A_row_start || !A_column_index || !A_value
      || !ptr_R_row_start || !ptr_R_column_index || !ptr_R_value)
      return QUERN_INPUT_ERROR;
   // set up lists for dynamically building R
   Quern::Pool<SparseEntry> pool;
   SparseVector* R=(SparseVector*)std::malloc(m*sizeof(SparseVector));
   if(!R) return QUERN_OUT_OF_MEMORY;
   for(int i=0; i<m; ++i) R[i].init(&pool);
   // do the Givens QR
   SparseVector row;
   row.init(&pool);
   double c, s;
   for(int a=0; a<m; ++a){
      int i=(row_order ? row_order[a] : a);
      if(!copy_row(A_row_start[i+1]-A_row_start[i],
                   A_column_index+A_row_start[i],
                   A_value+A_row_start[i],
                   row)){
         std::free(R);
         return QUERN_OUT_OF_MEMORY;
      }
      while(!row.empty() && row.front().index<a && row.front().index<n){
         int j=row.front().index;
         if(R[j].empty() || R[j].front().index>j){ // swap?
            R[j].swap(row);
         }else{ // use Givens
            if(!apply_givens(R[j], row, j, c, s)){
               std::free(R);
               return QUERN_OUT_OF_MEMORY;
            }
         }
      }
      if(a<n){
         R[a].swap(row);
         assert(R[a].empty() || R[a].front().index>=a);
      }
   }
   // transfer R's lists to static CSR form
   int* R_row_start=(int*)std::malloc((n+1)*sizeof(int));
   if(!R_row_start){
      std::free(R);
      return QUERN_OUT_OF_MEMORY;
   }
   R_row_start[0]=0;
   for(int i=0; i<n; ++i)
      R_row_start[i+1]=R_row_start[i]+R[i].size();
   int Rnnz=R_row_start[n];
   int* R_column_index=(int*)std::malloc(Rnnz*sizeof(int));
   if(!R_column_index){
      std::free(R);
      std::free(R_row_start);
      return QUERN_OUT_OF_MEMORY;
   }
   double* R_value=(double*)std::malloc(Rnnz*sizeof(double));
   if(!R_value){
      std::free(R);
      std::free(R_row_start);
      std::free(R_column_index);
      return QUERN_OUT_OF_MEMORY;
   }
   int j=0;
   for(int i=0; i<n; ++i){
      Quern::ListIterator<SparseEntry> p;
      for(p=R[i].begin(); p.still_going(); ++p){
         R_column_index[j]=p->index;
         R_value[j]=p->value;
         ++j;
      }
   }
   std::free(R);
   *ptr_R_row_start=R_row_start;
   *ptr_R_column_index=R_column_index;
   *ptr_R_value=R_value;
   return QUERN_OK;
}

int QUERN_compute_qr(int m,
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
                     double** ptr_R_value)
{
   if(m<=0 || n<=0 || m<n || !A_row_start || !A_column_index || !A_value
      || !ptr_Q_row_start || !ptr_Q_column_index || !ptr_Q_value
      || !ptr_R_row_start || !ptr_R_column_index || !ptr_R_value)
      return QUERN_INPUT_ERROR;
   // set up lists for dynamically building Q and R
   Quern::Pool<SparseEntry> pool;
   SparseVector* Q=(SparseVector*)std::malloc(m*sizeof(SparseVector));
   if(!Q) return QUERN_OUT_OF_MEMORY;
   SparseVector* R=(SparseVector*)std::malloc(m*sizeof(SparseVector));
   if(!R){
      std::free(Q);
      return QUERN_OUT_OF_MEMORY;
   }
   for(int i=0; i<m; ++i){
      Q[i].init(&pool);
      R[i].init(&pool);
   }
   // do the Givens QR
   SparseVector row;
   row.init(&pool);
   double c, s;
   for(int a=0; a<m; ++a){
      int i=(row_order ? row_order[a] : a);
      if(!copy_row(A_row_start[i+1]-A_row_start[i],
                   A_column_index+A_row_start[i],
                   A_value+A_row_start[i],
                   row)){
         std::free(Q);
         std::free(R);
         return QUERN_OUT_OF_MEMORY;
      }
      Quern::ListIterator<SparseEntry> q=Q[a].begin();
      while(!row.empty() && row.front().index<a && row.front().index<n){
         int j=row.front().index;
         if(R[j].empty() || R[j].front().index>j){ // swap?
            R[j].swap(row);
            Q[a].insert(q, SparseEntry(j, 1));
            ++q;
         }else{ // use Givens
            if(!apply_givens(R[j], row, j, c, s)){
               std::free(Q);
               std::free(R);
               return QUERN_OUT_OF_MEMORY;
            }
            Q[a].insert(q, SparseEntry(j, encode(c,s)));
            ++q;
         }
      }
      if(a<n){
         R[a].swap(row);
         assert(R[a].empty() || R[a].front().index>=a);
      }
   }
   // transfer Q's lists to static CSR form
   int* Q_row_start=(int*)std::malloc((m+1)*sizeof(int));
   if(!Q_row_start){
      std::free(Q);
      std::free(R);
      return QUERN_OUT_OF_MEMORY;
   }
   Q_row_start[0]=0;
   for(int i=0; i<m; ++i)
      Q_row_start[i+1]=Q_row_start[i]+Q[i].size();
   int Qnnz=Q_row_start[m];
   int* Q_column_index=(int*)std::malloc(Qnnz*sizeof(int));
   if(!Q_column_index){
      std::free(Q);
      std::free(R);
      std::free(Q_row_start);
      return QUERN_OUT_OF_MEMORY;
   }
   double* Q_value=(double*)std::malloc(Qnnz*sizeof(double));
   if(!Q_value){
      std::free(Q);
      std::free(R);
      std::free(Q_row_start);
      std::free(Q_column_index);
      return QUERN_OUT_OF_MEMORY;
   }
   int j=0;
   for(int i=0; i<m; ++i){
      Quern::ListIterator<SparseEntry> q;
      for(q=Q[i].begin(); q.still_going(); ++q){
         Q_column_index[j]=q->index;
         Q_value[j]=q->value;
         ++j;
      }
   }
   std::free(Q);
   // transfer R's lists to static CSR form
   int* R_row_start=(int*)std::malloc((n+1)*sizeof(int));
   if(!R_row_start){
      std::free(Q);
      std::free(R);
      std::free(Q_row_start);
      std::free(Q_column_index);
      std::free(Q_value);
      return QUERN_OUT_OF_MEMORY;
   }
   R_row_start[0]=0;
   for(int i=0; i<n; ++i)
      R_row_start[i+1]=R_row_start[i]+R[i].size();
   int Rnnz=R_row_start[n];
   int* R_column_index=(int*)std::malloc(Rnnz*sizeof(int));
   if(!R_column_index){
      std::free(Q);
      std::free(R);
      std::free(Q_row_start);
      std::free(Q_column_index);
      std::free(Q_value);
      std::free(R_row_start);
      return QUERN_OUT_OF_MEMORY;
   }
   double* R_value=(double*)std::malloc(Rnnz*sizeof(double));
   if(!R_value){
      std::free(Q);
      std::free(R);
      std::free(Q_row_start);
      std::free(Q_column_index);
      std::free(Q_value);
      std::free(R_row_start);
      std::free(R_column_index);
      return QUERN_OUT_OF_MEMORY;
   }
   j=0;
   for(int i=0; i<n; ++i){
      Quern::ListIterator<SparseEntry> p;
      for(p=R[i].begin(); p.still_going(); ++p){
         R_column_index[j]=p->index;
         R_value[j]=p->value;
         ++j;
      }
   }
   std::free(R);
   *ptr_Q_row_start=Q_row_start;
   *ptr_Q_column_index=Q_column_index;
   *ptr_Q_value=Q_value;
   *ptr_R_row_start=R_row_start;
   *ptr_R_column_index=R_column_index;
   *ptr_R_value=R_value;
   return QUERN_OK;
}

static void determine_column_thresholds(int m,
                                        int n,
                                        const int* A_row_start,
                                        const int* A_column_index,
                                        const double* A_value,
                                        double drop_tolerance,
                                        double* drop)
{
   std::memset(drop, 0, n*sizeof(double));
   for(int i=0; i<m; ++i){
      for(int j=A_row_start[i]; j<A_row_start[i+1]; ++j)
         drop[A_column_index[j]]+=A_value[j]*A_value[j];
   }
   for(int i=0; i<n; ++i){
      drop[i]=std::sqrt(drop[i])*drop_tolerance;
   }
}

static bool apply_givens(SparseVector& x,
                         SparseVector& y,
                         int diagonal,
                         const double* drop,
                         double& c,
                         double& s)
{
   assert(!x.empty() && x.front().index==diagonal && x.front().value);
   assert(!y.empty() && y.front().index==diagonal && y.front().value);
   // find the rotation we need
   double a=x.front().value, b=y.front().value;
   givens(a, b, c, s);
   assert(c && s);
   // rotate the start of each list
   x.front().value=c*a-s*b;
   y.pop_front();
   // then update the rest (x_new = c*x-s*y and y_new = s*x+c*y)
   Quern::ListIterator<SparseEntry> p=x.begin(), q=y.begin();
   ++p; // skip the first value we already took care of
   while(p.still_going() && q.still_going()){
      if(p->index==q->index){
         int k=p->index;
         double xnew=c*p->value-s*q->value,
                ynew=s*p->value+c*q->value;
         if(std::fabs(xnew)>drop[k]){
            p->value=xnew;
            ++p;
         }else x.erase(p);
         if(std::fabs(ynew)>drop[k]){
            q->value=ynew;
            ++q;
         }else y.erase(q);
      }else if(p->index<q->index){
         int k=p->index;
         double xnew=c*p->value,
                ynew=s*p->value;
         if(std::fabs(xnew)>drop[k]){
            p->value=xnew;
            ++p;
         }else x.erase(p);
         if(std::fabs(ynew)>drop[k]){
            if(!y.insert(q, SparseEntry(k, ynew))) return false;
            ++q;
         }
      }else{
         int k=q->index;
         double xnew=-s*q->value,
                ynew=c*q->value;
         if(std::fabs(xnew)>drop[k]){
            if(!x.insert(p, SparseEntry(k, xnew))) return false;
            ++p;
         }
         if(std::fabs(ynew)>drop[k]){
            q->value=ynew;
            ++q;
         }else y.erase(q);
      }
   }
   if(p.still_going()){
      do{
         int k=p->index;
         double xnew=c*p->value,
                ynew=s*p->value;
         if(std::fabs(xnew)>drop[k]){
            p->value=xnew;
            ++p;
         }else x.erase(p);
         if(std::fabs(ynew)>drop[k]){
            if(!y.insert(q, SparseEntry(k, ynew))) return false;
            ++q;
         }
      }while(p.still_going());
   }else if(q.still_going()){
      do{
         int k=q->index;
         double xnew=-s*q->value,
                ynew=c*q->value;
         if(std::fabs(xnew)>drop[k]){
            if(!x.insert(p, SparseEntry(k, xnew))) return false;
            ++p;
         }
         if(std::fabs(ynew)>drop[k]){
            q->value=ynew;
            ++q;
         }else y.erase(q);
      }while(q.still_going());
   }
   return true;
}

int QUERN_compute_incomplete_qr_without_q(int m,
                                          int n,
                                          const int* A_row_start,
                                          const int* A_column_index,
                                          const double* A_value,
                                          const int* row_order,
                                          double drop_tolerance,
                                          int** ptr_R_row_start,
                                          int** ptr_R_column_index,
                                          double** ptr_R_value)
{
   if(m<=0 || n<=0 || m<n || !A_row_start || !A_column_index || !A_value
      || !ptr_R_row_start || !ptr_R_column_index || !ptr_R_value)
      return QUERN_INPUT_ERROR;
   // find per-column drop tolerance
   double* drop=(double*)std::malloc(n*sizeof(double));
   if(!drop) return QUERN_OUT_OF_MEMORY;
   determine_column_thresholds(m, n, A_row_start, A_column_index, A_value,
                               drop_tolerance, drop);
   // set up lists for dynamically building R
   Quern::Pool<SparseEntry> pool;
   SparseVector* R=(SparseVector*)std::malloc(m*sizeof(SparseVector));
   if(!R){
      std::free(drop);
      return QUERN_OUT_OF_MEMORY;
   }
   for(int i=0; i<m; ++i) R[i].init(&pool);
   // do the Givens QR
   SparseVector row;
   row.init(&pool);
   double c, s;
   for(int a=0; a<m; ++a){
      int i=(row_order ? row_order[a] : a);
      if(!copy_row(A_row_start[i+1]-A_row_start[i],
                   A_column_index+A_row_start[i],
                   A_value+A_row_start[i],
                   row)){
         std::free(drop);
         std::free(R);
         return QUERN_OUT_OF_MEMORY;
      }
      while(!row.empty() && row.front().index<a && row.front().index<n){
         int j=row.front().index;
         if(R[j].empty() || R[j].front().index>j){ // swap?
            R[j].swap(row);
         }else{ // use Givens
            if(!apply_givens(R[j], row, j, drop, c, s)){
               std::free(drop);
               std::free(R);
               return QUERN_OUT_OF_MEMORY;
            }
         }
      }
      if(a<n){
         R[a].swap(row);
         assert(R[a].empty() || R[a].front().index>=a);
      }
   }
   // don't need this any more
   std::free(drop);
   // transfer R's lists to static CSR form
   int* R_row_start=(int*)std::malloc((n+1)*sizeof(int));
   if(!R_row_start){
      std::free(R);
      return QUERN_OUT_OF_MEMORY;
   }
   R_row_start[0]=0;
   for(int i=0; i<n; ++i)
      R_row_start[i+1]=R_row_start[i]+R[i].size();
   int Rnnz=R_row_start[n];
   int* R_column_index=(int*)std::malloc(Rnnz*sizeof(int));
   if(!R_column_index){
      std::free(R);
      std::free(R_row_start);
      return QUERN_OUT_OF_MEMORY;
   }
   double* R_value=(double*)std::malloc(Rnnz*sizeof(double));
   if(!R_value){
      std::free(R);
      std::free(R_row_start);
      std::free(R_column_index);
      return QUERN_OUT_OF_MEMORY;
   }
   int j=0;
   for(int i=0; i<n; ++i){
      Quern::ListIterator<SparseEntry> p;
      for(p=R[i].begin(); p.still_going(); ++p){
         R_column_index[j]=p->index;
         R_value[j]=p->value;
         ++j;
      }
   }
   std::free(R);
   *ptr_R_row_start=R_row_start;
   *ptr_R_column_index=R_column_index;
   *ptr_R_value=R_value;
   return QUERN_OK;
}

int QUERN_compute_incomplete_qr(int m,
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
                                double** ptr_R_value)
{
   if(m<=0 || n<=0 || m<n || !A_row_start || !A_column_index || !A_value
      || !ptr_Q_row_start || !ptr_Q_column_index || !ptr_Q_value
      || !ptr_R_row_start || !ptr_R_column_index || !ptr_R_value)
      return QUERN_INPUT_ERROR;
   // find per-column drop tolerance
   double* drop=(double*)std::malloc(n*sizeof(double));
   if(!drop) return QUERN_OUT_OF_MEMORY;
   determine_column_thresholds(m, n, A_row_start, A_column_index, A_value,
                               drop_tolerance, drop);
   // set up lists for dynamically building Q and R
   Quern::Pool<SparseEntry> pool;
   SparseVector* Q=(SparseVector*)std::malloc(m*sizeof(SparseVector));
   if(!Q){
      std::free(drop);
      return QUERN_OUT_OF_MEMORY;
   }
   SparseVector* R=(SparseVector*)std::malloc(m*sizeof(SparseVector));
   if(!R){
      std::free(drop);
      std::free(Q);
      return QUERN_OUT_OF_MEMORY;
   }
   for(int i=0; i<m; ++i){
      Q[i].init(&pool);
      R[i].init(&pool);
   }
   // do the Givens QR
   SparseVector row;
   row.init(&pool);
   double c, s;
   for(int a=0; a<m; ++a){
      int i=(row_order ? row_order[a] : a);
      if(!copy_row(A_row_start[i+1]-A_row_start[i],
                   A_column_index+A_row_start[i],
                   A_value+A_row_start[i],
                   row)){
         std::free(drop);
         std::free(Q);
         std::free(R);
         return QUERN_OUT_OF_MEMORY;
      }
      Quern::ListIterator<SparseEntry> q=Q[a].begin();
      while(!row.empty() && row.front().index<a && row.front().index<n){
         int j=row.front().index;
         if(R[j].empty() || R[j].front().index>j){ // swap?
            R[j].swap(row);
            Q[a].insert(q, SparseEntry(j, 1));
            ++q;
         }else{ // use Givens
            if(!apply_givens(R[j], row, j, drop, c, s)){
               std::free(drop);
               std::free(Q);
               std::free(R);
               return QUERN_OUT_OF_MEMORY;
            }
            Q[a].insert(q, SparseEntry(j, encode(c,s)));
            ++q;
         }
      }
      if(a<n){
         R[a].swap(row);
         assert(R[a].empty() || R[a].front().index>=a);
      }
   }
   // don't need this any more
   std::free(drop);
   // transfer Q's lists to static CSR form
   int* Q_row_start=(int*)std::malloc((m+1)*sizeof(int));
   if(!Q_row_start){
      std::free(Q);
      std::free(R);
      return QUERN_OUT_OF_MEMORY;
   }
   Q_row_start[0]=0;
   for(int i=0; i<m; ++i)
      Q_row_start[i+1]=Q_row_start[i]+Q[i].size();
   int Qnnz=Q_row_start[m];
   int* Q_column_index=(int*)std::malloc(Qnnz*sizeof(int));
   if(!Q_column_index){
      std::free(Q);
      std::free(R);
      std::free(Q_row_start);
      return QUERN_OUT_OF_MEMORY;
   }
   double* Q_value=(double*)std::malloc(Qnnz*sizeof(double));
   if(!Q_value){
      std::free(Q);
      std::free(R);
      std::free(Q_row_start);
      std::free(Q_column_index);
      return QUERN_OUT_OF_MEMORY;
   }
   int j=0;
   for(int i=0; i<m; ++i){
      Quern::ListIterator<SparseEntry> q;
      for(q=Q[i].begin(); q.still_going(); ++q){
         Q_column_index[j]=q->index;
         Q_value[j]=q->value;
         ++j;
      }
   }
   std::free(Q);
   // transfer R's lists to static CSR form
   int* R_row_start=(int*)std::malloc((n+1)*sizeof(int));
   if(!R_row_start){
      std::free(Q);
      std::free(R);
      std::free(Q_row_start);
      std::free(Q_column_index);
      std::free(Q_value);
      return QUERN_OUT_OF_MEMORY;
   }
   R_row_start[0]=0;
   for(int i=0; i<n; ++i)
      R_row_start[i+1]=R_row_start[i]+R[i].size();
   int Rnnz=R_row_start[n];
   int* R_column_index=(int*)std::malloc(Rnnz*sizeof(int));
   if(!R_column_index){
      std::free(Q);
      std::free(R);
      std::free(Q_row_start);
      std::free(Q_column_index);
      std::free(Q_value);
      std::free(R_row_start);
      return QUERN_OUT_OF_MEMORY;
   }
   double* R_value=(double*)std::malloc(Rnnz*sizeof(double));
   if(!R_value){
      std::free(Q);
      std::free(R);
      std::free(Q_row_start);
      std::free(Q_column_index);
      std::free(Q_value);
      std::free(R_row_start);
      std::free(R_column_index);
      return QUERN_OUT_OF_MEMORY;
   }
   j=0;
   for(int i=0; i<n; ++i){
      Quern::ListIterator<SparseEntry> p;
      for(p=R[i].begin(); p.still_going(); ++p){
         R_column_index[j]=p->index;
         R_value[j]=p->value;
         ++j;
      }
   }
   std::free(R);
   *ptr_Q_row_start=Q_row_start;
   *ptr_Q_column_index=Q_column_index;
   *ptr_Q_value=Q_value;
   *ptr_R_row_start=R_row_start;
   *ptr_R_column_index=R_column_index;
   *ptr_R_value=R_value;
   return QUERN_OK;
}

void QUERN_free_result(int* row_start,
                       int* column_index,
                       double* value)
{
   std::free(row_start);
   std::free(column_index);
   std::free(value);
}

int QUERN_multiply_with_q(int m,
                          const int* Q_row_start,
                          const int* Q_column_index,
                          const double* Q_value,
                          double* x)
{
   if(m<=0 || !Q_row_start || !Q_column_index || !Q_value)
      return QUERN_INPUT_ERROR;
   int i, j, k;
   double c, s;
   for(i=m-1; i>=0; --i){
      for(j=Q_row_start[i+1]-1; j>=Q_row_start[i]; --j){
         k=Q_column_index[j];
         if(Q_value[j]==1.){ // swap?
            std::swap(x[i], x[k]);
         }else{
            decode(Q_value[j], c, s);
            double newxk=c*x[k]+s*x[i];
            x[i]=-s*x[k]+c*x[i];
            x[k]=newxk;
         }
      }
   }
   return QUERN_OK;
}

int QUERN_multiply_with_q_transpose(int m,
                                    const int* Q_row_start,
                                    const int* Q_column_index,
                                    const double* Q_value,
                                    double* x)
{
   if(m<=0 || !Q_row_start || !Q_column_index || !Q_value)
      return QUERN_INPUT_ERROR;
   int i, j, k;
   double c, s;
   for(i=0; i<m; ++i){
      for(j=Q_row_start[i]; j<Q_row_start[i+1]; ++j){
         k=Q_column_index[j];
         if(Q_value[j]==1){ // swap?
            std::swap(x[i], x[k]);
         }else{
            decode(Q_value[j], c, s);
            double newxk=c*x[k]-s*x[i];
            x[i]=s*x[k]+c*x[i];
            x[k]=newxk;
         }
      }
   }
   return QUERN_OK;
}

int QUERN_solve_with_r(int n,
                       const int* R_row_start,
                       const int* R_column_index,
                       const double* R_value,
                       const double* rhs,
                       double* result)
{
   // check input
   if(n<=0 || !R_row_start || !R_column_index || !R_value || !rhs || !result)
      return QUERN_INPUT_ERROR;
   // do the solve
   double x, rii;
   int i, j;
   for(i=n-1; i>=0; --i){
      x=rhs[i];
      rii=0;
      j=R_row_start[i];
      if(j<R_row_start[i+1] && R_column_index[j]==i){
         rii=R_value[j];
         ++j;
      }
      if(rii){
         for(; j<R_row_start[i+1]; ++j){
            assert(R_column_index[j]>i && R_column_index[j]<n);
            x-=R_value[j]*result[R_column_index[j]];
         }
         result[i]=x/rii;
      }else
         result[i]=0;
   }
   return QUERN_OK;
}

int QUERN_solve_with_r_transpose_in_place(int n,
                                          const int* R_row_start,
                                          const int* R_column_index,
                                          const double* R_value,
                                          double* x)
{
   // check input
   if(n<=0 || !R_row_start || !R_column_index || !R_value || !x)
      return QUERN_INPUT_ERROR;
   // do the solve
   double rii;
   int i, j;
   for(i=0; i<n; ++i){
      rii=0;
      j=R_row_start[i];
      if(j<R_row_start[i+1] && R_column_index[j]==i){
         rii=R_value[j];
         ++j;
      }
      if(rii){
         x[i]/=rii;
         for(; j<R_row_start[i+1]; ++j){
            assert(R_column_index[j]>i && R_column_index[j]<n);
            x[R_column_index[j]]-=R_value[j]*x[i];
         }
      }else
         x[i]=0;
   }
   return QUERN_OK;
}
