
#include "decs.h"

/* 
 * calculate components of magnetosonic velocity 
 * corresponding to primitive variables p 
 *
 * cfg 7-10-01
 * 
 */

#define GTOL 1.e-8

int vchar(FTYPE *pr, struct of_state *q, int js, struct of_geom *geom,
	  FTYPE *vmax, FTYPE *vmin)
{
  FTYPE discr, vp, vm, bsq, EE, EF, va2, cs2, cms2, rho, u;
  FTYPE Acov[NDIM], Bcov[NDIM], Acon[NDIM], Bcon[NDIM];
  FTYPE Asq, Bsq, Au, Bu, AB, Au2, Bu2, AuBu, A, B, C;
  int j;

  /* for regions along poles */
  if (geom->g < GTOL) {
    *vmax = 0.;
    *vmin = 0.;
    return (0);
  }

  DLOOPA Acov[j] = 0.;
  Acov[js] = 1.;
  raise(Acov, geom, Acon);

  DLOOPA Bcov[j] = 0.;
  Bcov[TT] = 1.;
  raise(Bcov, geom, Bcon);

  /* find fast magnetosonic speed */
  bsq = dot(q->bcon, q->bcov);
  rho = pr[RHO];
  u = pr[UU];
  EF = rho + gam * u;
  EE = bsq + EF;
  va2 = bsq / EE;
  cs2 = gam * (gam - 1.) * u / EF;
  cms2 = cs2 + va2 - cs2 * va2;	/* and there it is... */

  // cms2 *= 1.1 ;

  /* check on it! */
  if (cms2 < 0.) {
    if (fail(FAIL_COEFF_NEG) >= 1){
      trifprintf("dir=%d : %21.15g\n %21.15g\n %21.15g\n %21.15g\n %21.15g\n %21.15g\n %21.15g\n %21.15g\n",js,bsq,rho,u,EF,EE,va2,cs2,cms2);
      return (1);
    }
    cms2 = SMALL;
  }
  if (cms2 > 1.) {
    if (fail(FAIL_COEFF_SUP) >= 1)
      return (1);
    cms2 = 1.;
  }

  /* now require that speed of wave measured by observer q->ucon is
     cms2 */
  Asq = dot(Acon, Acov);
  Bsq = dot(Bcon, Bcov);
  Au = dot(Acov, q->ucon);
  Bu = dot(Bcov, q->ucon);
  AB = dot(Acon, Bcov);
  Au2 = Au * Au;
  Bu2 = Bu * Bu;
  AuBu = Au * Bu;

  A = Bu2 - (Bsq + Bu2) * cms2;
  B = 2. * (AuBu - (AB + AuBu) * cms2);
  C = Au2 - (Asq + Au2) * cms2;

  discr = B * B - 4. * A * C;
  if (discr < 0.) {
    dualfprintf(fail_file, "\n\t %g %g %g %g %g\n", A, B, C, discr, cms2);
    dualfprintf(fail_file, "\n\t q->ucon: %g %g %g %g\n", q->ucon[0],
	    q->ucon[1], q->ucon[2], q->ucon[3]);
    dualfprintf(fail_file, "\n\t q->bcon: %g %g %g %g\n", q->bcon[0],
	    q->bcon[1], q->bcon[2], q->bcon[3]);
    dualfprintf(fail_file, "\n\t Acon: %g %g %g %g\n", Acon[0], Acon[1],
	    Acon[2], Acon[3]);
    dualfprintf(fail_file, "\n\t Bcon: %g %g %g %g\n", Bcon[0], Bcon[1],
	    Bcon[2], Bcon[3]);
    if (fail(FAIL_VCHAR_DISCR) >= 1)
      return (1);
    discr = 0.;
  }

  discr = sqrt(discr);
  vp = -(-B + discr) / (2. * A);
  vm = -(-B - discr) / (2. * A);

  if (vp > vm) {
    *vmax = vp;
    *vmin = vm;
  } else {
    *vmax = vm;
    *vmin = vp;
  }

  return (0);
}
