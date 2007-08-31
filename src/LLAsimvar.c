/*##############################################################################
#
# Copyright © 2007 Ivan Kojadinovic, Israël-César Lerman and Philippe Peter 
#
# Ivan.Kojadinovic@polytech.univ-nantes.fr
#
# This software is a package for the statistical system GNU R:
# http://www.r-project.org 
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use, 
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info". 
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability. 
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security. 
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
#
##############################################################################
*/

/*****************************************************************************

  LLA similarity indices
  Ivan Kojadinovic, Mar 2007

*****************************************************************************/

#include <R.h>
#include <Rmath.h>

/*****************************************************************************

  Computation of the Cramer-von Mises statistic used in the independence
  test of Genest and Remillard based on the empirical copula
  
******************************************************************************/

static double Dn(int n, int s, int t)
{
  return (2.0 * n + 1.0) / (6.0 * n) + (s * (s - 1.0) + t * (t - 1.0)) / (2.0 * n * (n + 1.0)) 
    - fmax2(s,t) / (n + 1.0);
}

double empirical_copula(int n, int *R, int *S)
{
  int i,j;
  double Bn = 0.0;
  for (i=0;i<n;i++)
    for (j=0;j<n;j++)
      Bn += Dn(n,R[i],R[j]) * Dn(n,S[i],S[j]);
  return Bn/n;
}

/*****************************************************************************

  Computation of the similarity matrix among variables based on the 
  above statistic
  
  x: data ranks (no NA allowed)
  n: # lines
  p: # columns
  s: vector of similarities of size (p-1)p/2

*****************************************************************************/

/* Simulate the distribution of Bn under independence */

void simulate_empirical_copula(int *n, double *Bn,  int *N)
{
  int i, k;
  double *u = (double *)R_alloc(*n, sizeof(double));
  double *v = (double *)R_alloc(*n, sizeof(double));
  int *R = (int *)R_alloc(*n, sizeof(int));
  int *S = (int *)R_alloc(*n, sizeof(int));

  GetRNGstate();

  for (k=0;k<*N;k++) {
    
    for (i=0;i<*n;i++) {  
      u[i] = unif_rand();
      v[i] = unif_rand();
      R[i] = i+1;
      S[i] = i+1;
    }
    rsort_with_index (u, R, *n);
    rsort_with_index (v, S, *n);

    Bn[k] = empirical_copula(*n,R,S);
  }

  PutRNGstate();
}

/* compute the similarity matrix */

void similarity_empirical_copula(int *x, int *n, int *p, double *s,  
				 double *Bn, int *N)
{
  int i, j, k, l, count;
  double npairs = choose(*p,2);
  double *s_raw = (double *)R_alloc(npairs, sizeof(double));
  double maxi = 0;
  
  l = 0;
  for (i=0;i<*p;i++)
    for (j=i+1;j<*p;j++) {
      s_raw[l] = empirical_copula(*n, x + (*n) * i, x + (*n) * j); 
      if (s_raw[l] > maxi)
	maxi = s_raw[l];
      count = 0;
      for (k=0;k<*N;k++)
	if (Bn[k] <= s_raw[l])
	  count ++;
      s[l++] = (double)(count)/(*N +1.0);
    }

  /* correction */
  for (i=0;i<npairs;i++)
    s[i] += s_raw[i]/(maxi * (*N + 2.0));

}
