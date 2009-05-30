#include "decs.h"
#include "metric.h"

// this file includes all coordinate transformations and velocity
// transformations
// No user functions
// unless don't want to use bl as initial coords or v as velocity to
// evolve

// performs all necessary transformations from bl-coords to our
// metric-coords to our metric-coord in nonuniform grid and to
// 3-velocity from 4-velocity
int bl2met2metp2v(FTYPE *pr, int i, int j)
{
  int k = 0;

  /* transform four-velocity from bl to our metric */
  if (bltomet(pr, i, j) >= 1)
    FAILSTATEMENT("init.c:init()", "bltoks()", 1);

  /* transform to prime coords */
  mettometp(pr, i, j);

  /* convert from 4-vel to 3-vel */
  u_to_v(pr, i, j);

  return (0);
}

/* transforms u^i to our metric from boyer-lindquist */
int bltomet(FTYPE *pr, int i, int j)
{
  FTYPE ucon[NDIM], tmp[NDIM];
  FTYPE trans[NDIM][NDIM];
  FTYPE X[NDIM], r, th;
  struct of_geom geom;
  int k;

  coord(i, j, CENT, X);
  bl_coord(X, &r, &th);

  blgset(i, j, &geom);
  if (ucon_calc(pr, &geom, ucon) >= 1) {
    dualfprintf(fail_file, "bltoks(ucon_calc): space-like error\n");
    return (1);
  }

  /* make transform matrix */
  // order for trans is [ourmetric][bl]
  // DLOOP trans[j][k] = 0. ;
  // DLOOPA trans[j][j] = 1. ;
  trans[0][0] = trans00;
  trans[0][1] = trans01;
  trans[0][2] = trans02;
  trans[0][3] = trans03;
  trans[1][0] = trans10;
  trans[1][1] = trans11;
  trans[1][2] = trans12;
  trans[1][3] = trans13;
  trans[2][0] = trans20;
  trans[2][1] = trans21;
  trans[2][2] = trans22;
  trans[2][3] = trans23;
  trans[3][0] = trans30;
  trans[3][1] = trans31;
  trans[3][2] = trans32;
  trans[3][3] = trans33;

  /* transform ucon; solve for v */
  DLOOPA tmp[j] = 0.;
  DLOOP tmp[j] += trans[j][k] * ucon[k];
  DLOOPA ucon[j] = tmp[j];

  pr[U1] = ucon[1] / ucon[0];
  pr[U2] = ucon[2] / ucon[0];
  pr[U3] = ucon[3] / ucon[0];

  /* done! */
  return (0);
}

void mettometp(FTYPE *pr, int i, int j)
{
  FTYPE r, th, X[NDIM];
  FTYPE dxdxp[NDIM];

  coord(i, j, CENT, X);
  bl_coord(X, &r, &th);


  dxdxprim(X, r, th, dxdxp);

  pr[UU] /= dxdxp[0];
  pr[U1] /= dxdxp[1];
  pr[U2] /= dxdxp[2];
  pr[U3] /= dxdxp[3];


  /* done! */
}

// convert u 4-velocity to v 3-velocity
void u_to_v(FTYPE *pr, int i, int j)
{
  FTYPE ucon[NDIM];
  FTYPE AA, BB, CC, discr;
  struct of_geom geom;

  blgset(i, j, &geom);

  ucon[1] = pr[U1];
  ucon[2] = pr[U2];
  ucon[3] = pr[U3];

  AA = geom.gcov[TT][TT];
  BB = 2. * (geom.gcov[TT][1] * ucon[1] +
	     geom.gcov[TT][2] * ucon[2] + geom.gcov[TT][3] * ucon[3]);
  CC = 1. +
      geom.gcov[1][1] * ucon[1] * ucon[1] +
      geom.gcov[2][2] * ucon[2] * ucon[2] +
      geom.gcov[3][3] * ucon[3] * ucon[3] +
      2. * (geom.gcov[1][2] * ucon[1] * ucon[2] +
	    geom.gcov[1][3] * ucon[1] * ucon[3] +
	    geom.gcov[2][3] * ucon[2] * ucon[3]);

  discr = BB * BB - 4. * AA * CC;

  ucon[TT] = (-BB - sqrt(discr)) / (2. * AA);

  pr[U1] *= 1. / ucon[TT];
  pr[U2] *= 1. / ucon[TT];
  pr[U3] *= 1. / ucon[TT];
}
