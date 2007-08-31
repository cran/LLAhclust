#include <R.h> 
#include <Rinternals.h> 
SEXP max_dist(SEXP x_sxp, SEXP y_sxp) { 
  if (isReal(x_sxp) == FALSE) 
    Rf_error(" ' x ' must be ' "); 
  if (isReal(y_sxp) == FALSE) 
    Rf_error(" ' y ' must be ' double ' "); 
  int x_len = LENGTH(x_sxp), 
    y_len = LENGTH(y_sxp); 
  SEXP dist_sxp; 
  PROTECT(dist_sxp = allocVector(REALSXP, x_len)); 
  double *x = REAL(x_sxp), *y = REAL(y_sxp), 
    *dist = REAL(dist_sxp);
  
  int i, j; 
  double cur; 
  for (i = 0; i < x_len; ++i) { 
    cur = 0; 
    for (j = 0; j < y_len; ++j) { 
      if (abs(x[i] - y[j]) > cur) 
	cur = abs(x[i] - y[j]); 
    } 
    dist[i] = cur; 
  } 
  UNPROTECT(1); 
  return dist_sxp; 
}
