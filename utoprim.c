
/* 
 *
 * invert U (conserved variables) to obtain
 * p (primitive variables).  
 * NB: pr must contain initial guess
 *
 */

#include "decs.h"
struct of_geom *ptrlgeom;
FTYPE U_target[NPR];

/* pr *MUST* contain initial guess */

// todo:
// 1) choose tolx/tolf more intelligently (perhaps scales with
// difference between maximum and minimum density)


/* pr *MUST* contain initial guess */
int Utoprim(FTYPE *U, struct of_geom *ptrgeom, FTYPE *pr)
{
  FTYPE tolx, tolf;
  int ntrial, k, test;

  ptrlgeom = ptrgeom;

  ntrial = 30;
  /* 
     tolx = 1.e-15; tolf = 1.e-15 ; */
  // choice
  tolx = 1.e-15;
  tolf = 1.e-15;

  if (U[0] < 0.) {
    if (fail(FAIL_UTOPRIM_NEG) >= 1)
      return (1);
  }

  PLOOP U_target[k] = U[k];

  for (k = B1; k <= B3; k++)
    pr[k] = U[k] / ptrgeom->g;	/* solution is known */
  test = mnewt(ntrial, pr - 1, NPR - 3, tolx, tolf);

  if (test >= 1) {
    if (fail(FAIL_UTOPRIM_TEST) >= 1)
      return (1);
  }

  return (0);
}


/* auxiliary function required by mnewt */
int usrfun(FTYPE *pr, int n, FTYPE *beta, FTYPE **alpha)
{
  static FTYPE U_curr[NPR];
  struct of_state q;
  int i = 0, j = 0, k = 0;

  /* normalize error = beta to \rho u^t */
  if (get_state(pr + 1, ptrlgeom, &q) >= 1)
    FAILSTATEMENT("utoprim.c:usrfun()", "get_state()", 1);
  if (primtoU(pr + 1, &q, ptrlgeom, U_curr) >= 1)
    FAILSTATEMENT("utoprim.c:usrfun()", "primtoU()", 1);
  for (k = 0; k < NPR - 3; k++)
    beta[k + 1] = (U_curr[k] - U_target[k]) / U_target[0];
  if (dudp_calc(pr + 1, &q, ptrlgeom, alpha) >= 1)
    FAILSTATEMENT("utoprim.c:usrfun()", "dudp_calc()", 1);

  for (j = 0; j < NPR - 3; j++)
    for (k = 0; k < NPR - 3; k++)
      alpha[j + 1][k + 1] /= U_target[0];

  return (0);
}
