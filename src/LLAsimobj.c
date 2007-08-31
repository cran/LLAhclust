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

  Normalize the vector of similarities
  S: vector of similarities 
  npairs: size i.e. (n-1)n/2

*****************************************************************************/

void normalize_similarity(double *S, int npairs)
{
  double mean, var, sd;
  int l; 

  mean = 0.0;
  for (l = 0 ; l < npairs; l++)
    mean += S[l];
  mean /= (double)npairs;
  
  var =0.0;
  for (l = 0 ; l < npairs; l++)
    var += R_pow_di(S[l] - mean,2);
  var /= (double)npairs;
  
  sd = sqrt(var);
  for (l = 0 ; l < npairs; l++)
    S[l] = (S[l]-mean)/sd; 
}

/*****************************************************************************

  Normalize the data (divide each line by its norm)

*****************************************************************************/

void normalize_data(double *x, int n, int p)
{
  int i, j, k;
  double norm;
  
  for (i = 0 ; i < n ; i++) {
    k = i;
    norm = 0.0;
    for(j = 0 ; j < p ; j++) {
      norm += R_pow_di(x[k],2);
      k += n;
    }
    norm = sqrt(norm);
    k = i;
    for(j = 0 ; j < p ; j++) {
      x[k] /= norm;
      k += n;
    }
  }
}

/*****************************************************************************

  Computation of the similarity matrix for objects 
  described by numerical variables using 
  a euclidean-distance based similarity
  x: data (no NA allowed)
  n: # lines
  p: # columns
  S: vector of similarities of size (n-1)n/2

*****************************************************************************/

void similarity_numerical_euclidean(double *x, int n, int p, double *S)
{
  int i, j, k, l, npairs =  n * (n - 1)/2, n2 = R_pow_di(n,2);
  double mean, var, sd, max;
  double *s = (double *)R_alloc(npairs, sizeof(double));
  
  for(j = 0 ; j < p ; j++) {
    l=0;
    mean = 0.0;
    max = 0.0;
    for (i = 0 ; i < n ; i++) 
      for (k = i+1 ; k < n ; k++) {
	s[l] = R_pow_di(x[i + n*j] - x[k + n*j],2);
	if (s[l] > max)
	  max = s[l];
	mean += s[l++];
      }
    mean = (n2 * max - 2.0 * mean) / n2;

    var =0.0;
    for (l = 0 ; l < npairs; l++)
      var += R_pow_di(s[l] - mean,2);
    
    var = (var * 2.0 + n * R_pow_di(max - mean,2)) / n2;
    
    sd = sqrt(var);
    for (l = 0 ; l < npairs; l++)
	S[l] += (s[l] - mean)/sd;
  }
}

/*****************************************************************************

  Computation of the similarity matrix for objects 
  described by numerical variables using a cosinus like similarity
  x: data (no NA allowed)
  n: # lines
  p: # columns
  S: vector of similarities of size (n-1)n/2

*****************************************************************************/

void similarity_numerical_cosinus(double *x, int n, int p, double *S)
{
  int i, j, k, l, npairs =  n * (n - 1)/2, n2 = R_pow_di(n,2);
  double mean, var, sd;
  double *s = (double *)R_alloc(npairs, sizeof(double));
  
  for(j = 0 ; j < p ; j++) {
    l=0;
    mean = 0.0;
    for (i = 0 ; i < n ; i++) 
      for (k = i+1 ; k < n ; k++) {
	s[l] = 1.0/p - 0.5 * R_pow_di(x[i + n*j] - x[k + n*j],2);
	mean += s[l++];
      }
    mean = (mean * 2.0 + n/(double)p) / n2;

    var =0.0;
    for (l = 0 ; l < npairs; l++)
      var += R_pow_di(s[l] - mean,2);
    
    var = (var * 2.0 + n * R_pow_di(1.0/p - mean,2)) / n2;
    
    sd = sqrt(var);
    for (l = 0 ; l < npairs; l++)
	S[l] += (s[l] - mean)/sd;
  }
}

/*****************************************************************************

  Computation of the similarity matrix for objects 
  described by categorical variables
  x: data (no NA allowed)
  n: # lines
  p: # columns
  S: vector of similarities of size (n-1)n/2

*****************************************************************************/

void similarity_categorical(double *x, int n, int p, double *S)
{
  int i, j, k, l, npairs =  n * (n - 1)/2, mi;
  double mean, var, sd, pi;
  double *s = (double *)R_alloc(npairs, sizeof(double));
  
  for(j = 0 ; j < p ; j++) {
    l=0;
    for (i = 0 ; i < n ; i++) 
      for (k = i+1 ; k < n ; k++) 
	s[l++] = (x[i + n*j] == x[k + n*j]) ? 1.0 : 0.0; 
    
    /* number of categories for column j */
    R_rsort (x + n*j, n);
    mi = 1;
    mean = 0.0;
    for (i = 0 ; i < n-1 ; i++) 
      if (x[i + n*j] == x[i + 1 + n*j]) 
	mi++;
      else {
	pi = mi/(double)n; 
	mean += R_pow_di(pi,2);
       	mi = 1;
      }
    pi = mi/(double)n; 
    mean += R_pow_di(pi,2);

    var = mean * ( 1.0 - mean);
    sd = sqrt(var);
       
    for (l = 0 ; l < npairs; l++)
	S[l] += (s[l] - mean)/sd;
  }
}

/*****************************************************************************

  Computation of the similarity matrix for objects 
  described by ordinal variables
  x: data (no NA allowed), numerically coded 
  n: # lines
  p: # columns
  S: vector of similarities of size (n-1)n/2

*****************************************************************************/
#define BLOCK_SIZE 32 

void similarity_ordinal(double *x, int n, int p, double *S)
{
  int i, j, k, l, npairs =  n * (n - 1)/2, hj, n2 = R_pow_di(n,2), 
    n4 = R_pow_di(n,4), incr;
  double mean, var, sd, sum1, sum2;
  double *s = (double *)R_alloc(npairs, sizeof(double));
  int old = BLOCK_SIZE;
  int *m = (int *)R_alloc(old, sizeof(int)); 

  for(j = 0 ; j < p ; j++) {

    /* similarity per variable */
    l=0;
    for (i = 0 ; i < n ; i++) 
      for (k = i+1 ; k < n ; k++) 
	s[l++] = fabs(x[i + n*j] - x[k + n*j]);

    
    /* number of categories for column j */
    R_rsort (x + n*j, n);
    hj=0;
    m[hj] = 1;
    for (i = 0 ; i < n-1 ; i++) 
      if (x[i + n*j] == x[i + 1 + n*j]) 
	m[hj]++;
      else {
	incr = x[i + 1 + n*j] - x[i + n*j];
	if (hj + incr >= old) {
	  m = (int *)S_realloc((char *)m, old + BLOCK_SIZE, old, sizeof(int));
	  old += BLOCK_SIZE;
	}
	for (k=1;k<incr;k++)
	  m[hj+k] = 0;
	hj += incr;
	m[hj] = 1;
      }
    hj++;
    
    /* computation of the expectation and the variance */ 
    sum1 = 0.0; sum2 = 0.0;
    for (i = 0 ; i < hj ; i++) 
      for (k = 0 ; k < i ; k++) {
	sum1 += m[i] * m[k] * (i - k);
	sum2 += m[i] * m[k] * R_pow_di(i - k,2);
      }
    mean = hj - 1.0 - 2.0/n2 * sum1;    
    var = 2.0/n2 * sum2 - 4.0/n4 * R_pow_di(sum1,2); 
    sd = sqrt(var);

    for (l = 0 ; l < npairs; l++)
      S[l] += (hj - 1.0 - s[l] - mean)/sd;
  }
}

/*****************************************************************************

  Computation of the similarity matrix for objects 
  x: data (no NA allowed)
  n: # lines
  p: # columns
  S: vector of similarities of size (n-1)n/2
  method: LLA similarity method

*****************************************************************************/

enum { NUMERICAL_EUCLIDEAN = 1, NUMERICAL_COSINUS, CATEGORICAL, ORDINAL, BOOLEAN};
/* == 1,2,..., defined by order in the R function LLAsim */

void similarity_objects(double *x, int *n, int *p, double *S, int *method)
{
  /* similarity per variable */
  void (*similarity)(double *, int , int , double *) = NULL;
  
  switch(*method) {
  case NUMERICAL_EUCLIDEAN:
    Rprintf("LLAsim: the objects are assumed to be described by numerical variables\n");
    similarity = similarity_numerical_euclidean;
    break;
  case NUMERICAL_COSINUS:
    Rprintf("LLAsim: the objects are assumed to be described by numerical variables\n");
    similarity = similarity_numerical_cosinus;
    normalize_data(x, *n, *p);
    break;
  case CATEGORICAL:
    Rprintf("LLAsim: the objects are assumed to be described by categorical variables\n");
    similarity = similarity_categorical;
    break;
  case ORDINAL:
    Rprintf("LLAsim: the objects are assumed to be described by ranks induced by ordinal variables\n");
    similarity = similarity_ordinal;
    break;
  case BOOLEAN:
    Rprintf("LLAsim: the objects are assumed to be described by boolean variables\n");
    normalize_data(x, *n, *p);
    similarity = similarity_numerical_cosinus;
    break;
  default:
    error("invalid similarity method");
  }

  /* construction of the similarity by additive decomposition */
  similarity(x,*n,*p, S);

  /* normalization of the similarity index */
  normalize_similarity(S,(*n) * (*n - 1)/2);

}










/* mean = 1.0/(double)(*p) - 1.0/(double)(*n) * Trj(x+(*n)*j,*n,2) 
       + 1.0/R_pow_di(*n,2) * R_pow_di(Trj(x+(*n)*j,*n,1),2); */

/* var = 1.0/(2.0 * (double)(*n)) * Trj(x+(*n)*j,*n,4) 
      - 2.0/R_pow_di(*n,2) * Trj(x+(*n)*j,*n,3) * Trj(x+(*n)*j,*n,1)
      + 1.0/(2* R_pow_di(*n,2)) * R_pow_di(Trj(x+(*n)*j,*n,2),2)
      + 2.0/R_pow_di(*n,3) * Trj(x+(*n)*j,*n,2) * R_pow_di(Trj(x+(*n)*j,*n,1),2) 
      - 1.0/R_pow_di(*n,4) * R_pow_di(Trj(x+(*n)*j,*n,1),4); */

/* double Trj(double *x,int n,int r)
{
  int i;
  double out = 0.0;
  for (i = 0; i < n; i++)
    out += R_pow_di(x[i],r);
  return out;
  } */
