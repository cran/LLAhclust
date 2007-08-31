/*
 * N. LeMeur  13/10/2008
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Utils.h> 
 
#include <stdlib.h>

/*-----------------------------------------------------------------
internal c function for calculation of interval
-----------------------------------------------------------------*/
SEXP inPartition(SEXP part){

  if(!isReal(part) & !isInteger(part))
    Rf_error("'data' must be a 'integer' vector");
  
  int nrow, ncol, i, j, k;
  SEXP ans; /* return value: an array*/
  
  nrow = ncol = LENGTH(part);
  
 
  /* allocate memory for return values */
  PROTECT(ans = allocMatrix(INTSXP, nrow, ncol));
  
  
  k = 0;
      /* Note that the matrix is walked column by column and not row by row.
         This is very important because that's how the elements of a matrix
         are stored in memory in R */
  for(j=0; j < ncol; j++){
    for(i=0; i < nrow; i++){
      INTEGER(ans)[k] = INTEGER(part)[i] == INTEGER(part)[j];  
      k++;
    }   
  }
  UNPROTECT(1); /* done with res  */
  return ans;
}




