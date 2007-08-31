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

  LLA aggregation criteria
  Ivan Kojadinovic, May 2007

*****************************************************************************/

#include <R.h>
#include <Rmath.h>

/*****************************************************************************

  Aggregration criterion based on Tippett's p-value 
  combination method

*****************************************************************************/

void  F77_SUB(tippett)(double *d1, double *d2, double *SCD1, 
		       double *SCD2, double *epsilon, double *sout)
{ 
  *sout = fmax2( R_pow(*SCD1, R_pow((*d1 + *d2)/(*d1), *epsilon)),
		R_pow(*SCD2, R_pow((*d1 + *d2)/(*d2),*epsilon))); 
}

/*****************************************************************************

  Aggregration criterion based on the maximum p-value 
  combination method

*****************************************************************************/

void F77_SUB(complete)(double *d1, double *d2, double *SCD1, 
		       double *SCD2, double *sout)
{ 
  *sout = 1.0 - fmax2( R_pow(1.0 - *SCD1, (*d1 + *d2)/(*d1)),
		      R_pow(1.0 - *SCD2, (*d1 + *d2)/(*d2))); 
}

/*****************************************************************************

  Aggregration criterion based on Fisher's p-value 
  combination method

*****************************************************************************/

void  F77_SUB(fisher)(double *c, double *d1, double *d2, 
		      double *SCD1, double *SCD2, double *sout)
{
  *sout = pchisq(qchisq( *SCD1, 2 * (*c) * (*d1), 1, 0) 
		 + qchisq( *SCD2, 2 * (*c) * (*d2), 1, 0), 
		 2 * (*c) * ((*d1) + (*d2)), 1, 0);
}

/*****************************************************************************

  Aggregration criterion based on the normal p-value 
  combination method

*****************************************************************************/

void  F77_SUB(normal)(double *d1, double *d2, double *SCD1, 
		      double *SCD2, double *sout)
{ 
  *sout = pnorm(sqrt(*d1/(*d1+*d2)) * qnorm(*SCD1, 0.0,  1.0, 1, 0)  
		+ sqrt(*d2/(*d1+*d2)) * qnorm(*SCD2, 0.0,  1.0, 1, 0), 
		0.0, 1.0, 1, 0);
}

/*****************************************************************************

  CDF of the sum of uniform independent [0,1] random variables

*****************************************************************************/

double truncated_pow(double x, int n)
{
  if (x > 0.0)
    return R_pow_di(x,n);
  else 
    return 0.0;
}

#define MAXN 20

void  cdf_sum_uniform(double *y, int *n, double *F)
{
  int i, sign = 1;

  if (*n > MAXN)
    {
      *F = pnorm(*y, *n / 2.0, sqrt(*n / 12.0), 1, 0);
      return;
    }

  *F = 0.0;
  for (i=0;i<=*n;i++)
    {      
      *F += sign * truncated_pow(*y - i, *n) 
	/ (gammafn(i+1) * gammafn(*n-i+1));
       sign = - sign;
    } 
}

double psumunif(double y, int n)
{
  int i, sign = 1;
  double F = 0.0;
  
  if (n > MAXN)
    return pnorm(y, n / 2.0, sqrt(n / 12.0), 1, 0);

  for (i=0;i<=n;i++)
    { 
      F += sign * truncated_pow(y - i, n) 
	/ (gammafn(i+1) * gammafn(n-i+1));
      sign = - sign;
    } 
  return F;
}


/*****************************************************************************

  Map row I and column J of upper half diagonal symmetric matrix     
  onto vector.

*****************************************************************************/

int offset(int n, int i, int j)
{
  return j + (i - 1) * n - (i * (i + 1)) / 2 - 1;
}

/*****************************************************************************

  Aggregration criterion based on the uniform p-value 
  combination method

*****************************************************************************/

void  F77_SUB(uniform)(int *k, int *ia, int *ib, int *na, int *n, double *s, 
		       double *sout)
{
  int *C = (int *)R_alloc(*na, sizeof(int));
  int *D1 = (int *)R_alloc(*na, sizeof(int));
  int *D2 = (int *)R_alloc(*na, sizeof(int));
  int i, j, c, d1, d2;
  double linkage = 0.0;

  c = 0;
  C[c++] = *k;
  for (i=*na-2;i>=0;i--)
    for (j=0;j<c;j++)
      if (C[j] == ia[i])
	C[c++] = ib[i];
  
  /*  Rprintf("C: ");
      for (i=0;i<c;i++)
      Rprintf("%d\t",C[i]);
      Rprintf("\n\n");*/
  

  d1 = 0;
  D1[d1++] = ia[*na-1];
  for (i=*na-2;i>=0;i--)
    for (j=0;j<d1;j++)
      if (D1[j] == ia[i])
	D1[d1++] = ib[i];
  
  /*  Rprintf("D1: ");
      for (i=0;i<d1;i++)
      Rprintf("%d\t",D1[i]);
      Rprintf("\n\n"); */
  
  d2 = 0;
  D2[d2++] = ib[*na-1];
  for (i=*na-2;i>=0;i--)
    for (j=0;j<d2;j++)
      if (D2[j] == ia[i])
	D2[d2++] = ib[i];
  
  /*  Rprintf("D2: ");
      for (i=0;i<d2;i++)
      Rprintf("%d\t",D2[i]);
      Rprintf("\n\n"); */

  for (i=0;i<c;i++)
    {
      for (j=0;j<d1;j++)
	linkage += s[offset(*n,imin2(C[i],D1[j]),imax2(C[i],D1[j]))];
      for (j=0;j<d2;j++)
	linkage += s[offset(*n,imin2(C[i],D2[j]),imax2(C[i],D2[j]))];
    }    

  *sout =  psumunif(linkage,c * (d1 + d2));
}

/**** BELOW NOT USED ***/


/*****************************************************************************

  Aggregration criterion based on Fisher's p-value 
  combination method; subset version

*****************************************************************************/

double ddf(int c, int d, int m)
{
  int j;
  double sum1 = 0.0, sum2 = 0.0;
  for (j=1;j<m;j++) {
    sum1 += choose(d,j);
    sum2 += choose(c,j);
  }
  return 2.0 * c * sum1 + 2.0 * d * sum2;
}

int subset2binary(int *x, int nx) 
{
  int i, b;

  b=0;
  for(i=0; i<nx; i++)
    b += 1 << x[i];
  return(b);
}

void  F77_SUB(fishersubset)(int *p, int *c, int *ia, int *ib, int *na, 
			    double *SCD1, double *SCD2, int *m)
{
  int *D1 = (int *)R_alloc(*na, sizeof(int));
  int *D2 = (int *)R_alloc(*na, sizeof(int));
  
  int i,j,d1,d2, D1bin, D2bin;
  double sumpval;
  
  /* recover the elements of D1 and D2 */
  d1 = 0;
  D1[d1++] = ia[*na-1];
  for (i=*na-2;i>=0;i--)
    for (j=0;j<d1;j++)
      if (D1[j] == ia[i])
	D1[d1++] = ib[i];
  d2 = 0;
  D2[d2++] = ib[*na-1];
  for (i=*na-2;i>=0;i--)
    for (j=0;j<d2;j++)
      if (D2[j] == ia[i])
	D2[d2++] = ib[i];

  /* convert to binary */
  D1bin = subset2binary(D1,d1); /* attention decal de 1 ! */
  D2bin = subset2binary(D2,d2);

  sumpval = 0.0;
  for (i=0;i<(1<<*p);i++)
    if (((i & (D1bin + D2bin)) == i) && (i & D1bin) && (i & D2bin))
      sumpval += -2 * log(1.0);
  

  pchisq(qchisq( *SCD1,  ddf(*c, d1, *m), 1, 0) 
	 + qchisq( *SCD2,  ddf(*c, d2, *m), 1, 0) 
	 + sumpval, 
	 ddf(*c,d1+d2,*m), 1, 0);

}

