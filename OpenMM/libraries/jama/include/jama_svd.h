#ifndef JAMA_SVD_H
#define JAMA_SVD_H


#include "tnt_array1d.h"
#include "tnt_array1d_utils.h"
#include "tnt_array2d.h"
#include "tnt_array2d_utils.h"
#include "tnt_math_utils.h"

#include <algorithm>
// for min(), max() below
#include <cmath>
// for abs() below

using namespace TNT;
using namespace std;

namespace JAMA
{
   /** Singular Value Decomposition.
   <P>
   For an m-by-n matrix A with m >= n, the singular value decomposition is
   an m-by-n orthogonal matrix U, an n-by-n diagonal matrix S, and
   an n-by-n orthogonal matrix V so that A = U*S*V'.
   <P>
   The singular values, sigma[k] = S[k][k], are ordered so that
   sigma[0] >= sigma[1] >= ... >= sigma[n-1].
   <P>
   The singular value decompostion always exists, so the constructor will
   never fail.  The matrix condition number and the effective numerical
   rank can be computed from this decomposition.

   <p>
	(Adapted from JAMA, a Java Matrix Library, developed by jointly 
	by the Mathworks and NIST; see  http://math.nist.gov/javanumerics/jama).
   */
template <class Real>
class SVD 
{


	Array2D<Real> U, V;
	Array1D<Real> s;
	int m, n;

  public:


   SVD (const Array2D<Real> &Arg) {


      m = Arg.dim1();
      n = Arg.dim2();
      int nu = min(m,n);
      s = Array1D<Real>(min(m+1,n)); 
      U = Array2D<Real>(m, nu, Real(0));
      V = Array2D<Real>(n,n);
      Array1D<Real> e(n);
      Array1D<Real> work(m);
	  Array2D<Real> A(Arg.copy());
      int wantu = 1;  					/* boolean */
      int wantv = 1;  					/* boolean */
	  int i=0, j=0, k=0;

      // Reduce A to bidiagonal form, storing the diagonal elements
      // in s and the super-diagonal elements in e.

      int nct = min(m-1,n);
      int nrt = max(0,min(n-2,m));
      for (k = 0; k < max(nct,nrt); k++) {
         if (k < nct) {

            // Compute the transformation for the k-th column and
            // place the k-th diagonal in s[k].
            // Compute 2-norm of k-th column without under/overflow.
            s[k] = 0;
            for (i = k; i < m; i++) {
               s[k] = hypot(s[k],A[i][k]);
            }
            if (s[k] != 0.0) {
               if (A[k][k] < 0.0) {
                  s[k] = -s[k];
               }
               for (i = k; i < m; i++) {
                  A[i][k] /= s[k];
               }
               A[k][k] += 1.0;
            }
            s[k] = -s[k];
         }
         for (j = k+1; j < n; j++) {
            if ((k < nct) && (s[k] != 0.0))  {

            // Apply the transformation.

               Real t(0.0);
               for (i = k; i < m; i++) {
                  t += A[i][k]*A[i][j];
               }
               t = -t/A[k][k];
               for (i = k; i < m; i++) {
                  A[i][j] += t*A[i][k];
               }
            }

            // Place the k-th row of A into e for the
            // subsequent calculation of the row transformation.

            e[j] = A[k][j];
         }
         if (wantu & (k < nct)) {

            // Place the transformation in U for subsequent back
            // multiplication.

            for (i = k; i < m; i++) {
               U[i][k] = A[i][k];
            }
         }
         if (k < nrt) {

            // Compute the k-th row transformation and place the
            // k-th super-diagonal in e[k].
            // Compute 2-norm without under/overflow.
            e[k] = 0;
            for (i = k+1; i < n; i++) {
               e[k] = hypot(e[k],e[i]);
            }
            if (e[k] != 0.0) {
               if (e[k+1] < 0.0) {
                  e[k] = -e[k];
               }
               for (i = k+1; i < n; i++) {
                  e[i] /= e[k];
               }
               e[k+1] += 1.0;
            }
            e[k] = -e[k];
            if ((k+1 < m) & (e[k] != 0.0)) {

            // Apply the transformation.

               for (i = k+1; i < m; i++) {
                  work[i] = 0.0;
               }
               for (j = k+1; j < n; j++) {
                  for (i = k+1; i < m; i++) {
                     work[i] += e[j]*A[i][j];
                  }
               }
               for (j = k+1; j < n; j++) {
                  Real t(-e[j]/e[k+1]);
                  for (i = k+1; i < m; i++) {
                     A[i][j] += t*work[i];
                  }
               }
            }
            if (wantv) {

            // Place the transformation in V for subsequent
            // back multiplication.

               for (i = k+1; i < n; i++) {
                  V[i][k] = e[i];
               }
            }
         }
      }

      // Set up the final bidiagonal matrix or order p.

      int p = min(n,m+1);
      if (nct < n) {
         s[nct] = A[nct][nct];
      }
      if (m < p) {
         s[p-1] = 0.0;
      }
      if (nrt+1 < p) {
         e[nrt] = A[nrt][p-1];
      }
      e[p-1] = 0.0;

      // If required, generate U.

      if (wantu) {
         for (j = nct; j < nu; j++) {
            for (i = 0; i < m; i++) {
               U[i][j] = 0.0;
            }
            U[j][j] = 1.0;
         }
         for (k = nct-1; k >= 0; k--) {
            if (s[k] != 0.0) {
               for (j = k+1; j < nu; j++) {
                  Real t(0.0);
                  for (i = k; i < m; i++) {
                     t += U[i][k]*U[i][j];
                  }
                  t = -t/U[k][k];
                  for (i = k; i < m; i++) {
                     U[i][j] += t*U[i][k];
                  }
               }
               for (i = k; i < m; i++ ) {
                  U[i][k] = -U[i][k];
               }
               U[k][k] = 1.0 + U[k][k];
               for (i = 0; i < k-1; i++) {
                  U[i][k] = 0.0;
               }
            } else {
               for (i = 0; i < m; i++) {
                  U[i][k] = 0.0;
               }
               U[k][k] = 1.0;
            }
         }
      }

      // If required, generate V.

      if (wantv) {
         for (k = n-1; k >= 0; k--) {
            if ((k < nrt) & (e[k] != 0.0)) {
               for (j = k+1; j < nu; j++) {
                  Real t(0.0);
                  for (i = k+1; i < n; i++) {
                     t += V[i][k]*V[i][j];
                  }
                  t = -t/V[k+1][k];
                  for (i = k+1; i < n; i++) {
                     V[i][j] += t*V[i][k];
                  }
               }
            }
            for (i = 0; i < n; i++) {
               V[i][k] = 0.0;
            }
            V[k][k] = 1.0;
         }
      }

      // Main iteration loop for the singular values.

      int pp = p-1;
      int iter = 0;
      Real eps(pow(2.0,-52.0));
      while (p > 0) {
         int k=0;
		 int kase=0;

         // Here is where a test for too many iterations would go.

         // This section of the program inspects for
         // negligible elements in the s and e arrays.  On
         // completion the variables kase and k are set as follows.

         // kase = 1     if s(p) and e[k-1] are negligible and k<p
         // kase = 2     if s(k) is negligible and k<p
         // kase = 3     if e[k-1] is negligible, k<p, and
         //              s(k), ..., s(p) are not negligible (qr step).
         // kase = 4     if e(p-1) is negligible (convergence).

         for (k = p-2; k >= -1; k--) {
            if (k == -1) {
               break;
            }
            if (abs(e[k]) <= eps*(abs(s[k]) + abs(s[k+1]))) {
               e[k] = 0.0;
               break;
            }
         }
         if (k == p-2) {
            kase = 4;
         } else {
            int ks;
            for (ks = p-1; ks >= k; ks--) {
               if (ks == k) {
                  break;
               }
               Real t( (ks != p ? abs(e[ks]) : 0.) + 
                          (ks != k+1 ? abs(e[ks-1]) : 0.));
               if (abs(s[ks]) <= eps*t)  {
                  s[ks] = 0.0;
                  break;
               }
            }
            if (ks == k) {
               kase = 3;
            } else if (ks == p-1) {
               kase = 1;
            } else {
               kase = 2;
               k = ks;
            }
         }
         k++;

         // Perform the task indicated by kase.

         switch (kase) {

            // Deflate negligible s(p).

            case 1: {
               Real f(e[p-2]);
               e[p-2] = 0.0;
               for (j = p-2; j >= k; j--) {
                  Real t( hypot(s[j],f));
                  Real cs(s[j]/t);
                  Real sn(f/t);
                  s[j] = t;
                  if (j != k) {
                     f = -sn*e[j-1];
                     e[j-1] = cs*e[j-1];
                  }
                  if (wantv) {
                     for (i = 0; i < n; i++) {
                        t = cs*V[i][j] + sn*V[i][p-1];
                        V[i][p-1] = -sn*V[i][j] + cs*V[i][p-1];
                        V[i][j] = t;
                     }
                  }
               }
            }
            break;

            // Split at negligible s(k).

            case 2: {
               Real f(e[k-1]);
               e[k-1] = 0.0;
               for (j = k; j < p; j++) {
                  Real t(hypot(s[j],f));
                  Real cs( s[j]/t);
                  Real sn(f/t);
                  s[j] = t;
                  f = -sn*e[j];
                  e[j] = cs*e[j];
                  if (wantu) {
                     for (i = 0; i < m; i++) {
                        t = cs*U[i][j] + sn*U[i][k-1];
                        U[i][k-1] = -sn*U[i][j] + cs*U[i][k-1];
                        U[i][j] = t;
                     }
                  }
               }
            }
            break;

            // Perform one qr step.

            case 3: {

               // Calculate the shift.
   
               Real scale = max(max(max(max(
                       abs(s[p-1]),abs(s[p-2])),abs(e[p-2])), 
                       abs(s[k])),abs(e[k]));
               Real sp = s[p-1]/scale;
               Real spm1 = s[p-2]/scale;
               Real epm1 = e[p-2]/scale;
               Real sk = s[k]/scale;
               Real ek = e[k]/scale;
               Real b = ((spm1 + sp)*(spm1 - sp) + epm1*epm1)/2.0;
               Real c = (sp*epm1)*(sp*epm1);
               Real shift = 0.0;
               if ((b != 0.0) || (c != 0.0)) {
                  shift = sqrt(b*b + c);
                  if (b < 0.0) {
                     shift = -shift;
                  }
                  shift = c/(b + shift);
               }
               Real f = (sk + sp)*(sk - sp) + shift;
               Real g = sk*ek;
   
               // Chase zeros.
   
               for (j = k; j < p-1; j++) {
                  Real t = hypot(f,g);
                  Real cs = f/t;
                  Real sn = g/t;
                  if (j != k) {
                     e[j-1] = t;
                  }
                  f = cs*s[j] + sn*e[j];
                  e[j] = cs*e[j] - sn*s[j];
                  g = sn*s[j+1];
                  s[j+1] = cs*s[j+1];
                  if (wantv) {
                     for (i = 0; i < n; i++) {
                        t = cs*V[i][j] + sn*V[i][j+1];
                        V[i][j+1] = -sn*V[i][j] + cs*V[i][j+1];
                        V[i][j] = t;
                     }
                  }
                  t = hypot(f,g);
                  cs = f/t;
                  sn = g/t;
                  s[j] = t;
                  f = cs*e[j] + sn*s[j+1];
                  s[j+1] = -sn*e[j] + cs*s[j+1];
                  g = sn*e[j+1];
                  e[j+1] = cs*e[j+1];
                  if (wantu && (j < m-1)) {
                     for (i = 0; i < m; i++) {
                        t = cs*U[i][j] + sn*U[i][j+1];
                        U[i][j+1] = -sn*U[i][j] + cs*U[i][j+1];
                        U[i][j] = t;
                     }
                  }
               }
               e[p-2] = f;
               iter = iter + 1;
            }
            break;

            // Convergence.

            case 4: {

               // Make the singular values positive.
   
               if (s[k] <= 0.0) {
                  s[k] = (s[k] < 0.0 ? -s[k] : 0.0);
                  if (wantv) {
                     for (i = 0; i <= pp; i++) {
                        V[i][k] = -V[i][k];
                     }
                  }
               }
   
               // Order the singular values.
   
               while (k < pp) {
                  if (s[k] >= s[k+1]) {
                     break;
                  }
                  Real t = s[k];
                  s[k] = s[k+1];
                  s[k+1] = t;
                  if (wantv && (k < n-1)) {
                     for (i = 0; i < n; i++) {
                        t = V[i][k+1]; V[i][k+1] = V[i][k]; V[i][k] = t;
                     }
                  }
                  if (wantu && (k < m-1)) {
                     for (i = 0; i < m; i++) {
                        t = U[i][k+1]; U[i][k+1] = U[i][k]; U[i][k] = t;
                     }
                  }
                  k++;
               }
               iter = 0;
               p--;
            }
            break;
         }
      }
   }


   void getU (Array2D<Real> &A) 
   {
   	  int minm = min(m+1,n);

	  A = Array2D<Real>(m, minm);

	  for (int i=0; i<m; i++)
	  	for (int j=0; j<minm; j++)
			A[i][j] = U[i][j];
   	
   }

   /* Return the right singular vectors */

   void getV (Array2D<Real> &A) 
   {
   	  A = V;
   }

   /** Return the one-dimensional array of singular values */

   void getSingularValues (Array1D<Real> &x) 
   {
      x = s;
   }

   /** Return the diagonal matrix of singular values
   @return     S
   */

   void getS (Array2D<Real> &A) {
   	  A = Array2D<Real>(n,n);
      for (int i = 0; i < n; i++) {
         for (int j = 0; j < n; j++) {
            A[i][j] = 0.0;
         }
         A[i][i] = s[i];
      }
   }

   /** Two norm  (max(S)) */

   Real norm2 () {
      return s[0];
   }

   /** Two norm of condition number (max(S)/min(S)) */

   Real cond () {
      return s[0]/s[min(m,n)-1];
   }

   /** Effective numerical matrix rank
   @return     Number of nonnegligible singular values.
   */

   int rank () 
   {
      Real eps = pow(2.0,-52.0);
      Real tol = max(m,n)*s[0]*eps;
      int r = 0;
      for (int i = 0; i < s.dim(); i++) {
         if (s[i] > tol) {
            r++;
         }
      }
      return r;
   }
};

}
#endif
// JAMA_SVD_H
