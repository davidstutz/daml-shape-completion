#ifndef GSL_POLY_SOLVE_HEADER
#define GSL_POLY_SOLVE_HEADER

int 
gsl_poly_solve_quadratic (double a, double b, double c, 
                          double *x0, double *x1);

int 
gsl_poly_solve_cubic (double a, double b, double c, 
                      double *x0, double *x1, double *x2);

#endif /* GSL_POLY_SOLVE_HEADER */
