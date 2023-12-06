#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(g_est)(int *y_rank, int *z_sort, int *tz, int *n, 
  int *p, int *K, double *delta0, double *pro0, double *ret);
extern void F77_NAME(alpha_fn)(double *log_dens, double *pro, 
  int *n, int *K, double *ret);
extern void F77_NAME(z_fn)(double *g, double *mu, double *pro, double *gam, 
  int *n, int *p, int *K, double *ret);
extern void F77_NAME(sigma_est)(double *mu, double *g, double *alpha, 
  int *n, int *p, int *K, double *ninv, double *ret);
extern void F77_NAME(loglik)(double *y, int *n, int *p, int *K, double *mu, 
  double *ginv, double *det_val, double *ret);
extern void F77_NAME(qnorm)(double *pval, double *qval);
extern void F77_NAME(msda)(double *obj, int *nk, int *nvars, double *sigma,
  double *delta, double *pf, int *dfmax, int *pmax, int *nlam,
  double *flmin, double *ulam, double *eps, int *maxit, double *sml,
  int *verbose, int *nalam, double *theta, int *m, int *ntheta, double *alam,
  int *npass, int *jerr);
extern void F77_NAME(msda_ncx)(double *obj, int *nk, int *nvars, double *sigma,
  double *delta, double *pfmat, int *dfmax, int *pmax, int *nlam,
  double *flmin, double *ulam, double *eps, int *maxit, double *sml,
  int *verbose, int *nalam, double *theta, int *m, int *ntheta, double *alam,
  int *npass, int *jerr);

static const R_FortranMethodDef FortranEntries[] = {
    {"g_est",  (DL_FUNC) &F77_NAME(g_est), 9},
    {"alpha_fn", (DL_FUNC) &F77_NAME(alpha_fn), 5},
    {"z_fn",  (DL_FUNC) &F77_NAME(z_fn), 8},
    {"sigma_est",  (DL_FUNC) &F77_NAME(sigma_est), 8},
    {"loglik",  (DL_FUNC) &F77_NAME(loglik), 8},
    {"qnorm",  (DL_FUNC) &F77_NAME(qnorm), 2},
    {"msda",  (DL_FUNC) &F77_NAME(msda), 22},
    {"msda_ncx",  (DL_FUNC) &F77_NAME(msda_ncx), 22},
    {NULL, NULL, 0}
};

void R_init_cesme(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
