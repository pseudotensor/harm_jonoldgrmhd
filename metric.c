#include "decs.h"
#include "metric.h"		// defines all coordinate metrics and
				// MCOORD

// this file includes metric dependent terms, including for initial
// condition routines for IC coords.

// only needs to be modified if more complex subterms need to be
// defined
// otherwise define metric in metric.h


// /////////////////////////////////////////////////////////////
// 
// Set metric/coord system here!
// 
// /////////////////////////////////////////////////////////////

void gcov_func(FTYPE *X, FTYPE gcov[][NDIM])
{
  int j, k;
  FTYPE sth, cth, s2, rho2;
  FTYPE r, th;
  FTYPE tfac, rfac, hfac, pfac;
  FTYPE dxdxp[NDIM];

  DLOOP gcov[j][k] = 0.;

  bl_coord(X, &r, &th);

  cth = cos(th);
  if(POSDEFMETRIC){
    sth = fabs(sin(th));
  }
  else{
    sth = sin(th);
  }
  if (fabs(sth) < SMALL)
    {
      if(sth>=0) sth=SMALL;
      if(sth<0) sth=-SMALL;
    }
  s2 = sth * sth;
  rho2 = r * r + a * a * cth * cth;

  // dx/dx' where '=prim coords (i.e. nonuni coords)

  dxdxprim(X, r, th, dxdxp);

  tfac = dxdxp[0];
  rfac = dxdxp[1];
  hfac = dxdxp[2];
  pfac = dxdxp[3];

  // now take term by term:
  // g_{u v} = \vec{e_{\mu}}\cdot\vec{e_{\nu}} 
  //           * (dx/dx')_{mu} * (dx/dx')_{\nu} =
  //          \vec{e'_{\mu}}\cdot\vec{e'_{\nu}} 

  gcov[TT][TT] = gcov00 * tfac * tfac;
  gcov[TT][1] = gcov01 * tfac * rfac;
  gcov[TT][2] = gcov02 * tfac * hfac;
  gcov[TT][3] = gcov03 * tfac * pfac;

  gcov[1][TT] = gcov10 * rfac * tfac;
  gcov[1][1] = gcov11 * rfac * rfac;
  gcov[1][2] = gcov12 * rfac * hfac;
  gcov[1][3] = gcov13 * rfac * pfac;

  gcov[2][TT] = gcov20 * hfac * tfac;
  gcov[2][1] = gcov21 * hfac * rfac;
  gcov[2][2] = gcov22 * hfac * hfac;
  gcov[2][3] = gcov23 * hfac * pfac;

  gcov[3][TT] = gcov30 * pfac * tfac;
  gcov[3][1] = gcov31 * pfac * rfac;
  gcov[3][2] = gcov32 * pfac * hfac;
  gcov[3][3] = gcov33 * pfac * pfac;
}

// ///////////////////////////////////////////////////////
// 
// Boyer-Lindquist ("bl") metric functions */
// The below functions used for bl coords only (r,th,phi)
// 
// bl coords are starting point for most IC and other metric/coords
// 
// ///////////////////////////////////////////////////////////

// find the con/cov forms of the bl metric
void blgset(int i, int j, struct of_geom *geom)
{
  FTYPE r, th, X[NDIM];

  coord(i, j, CENT, X);
  bl_coord(X, &r, &th);

  if (th < 0)
    th *= -1.;
  if (th > M_PI)
    th = 2. * M_PI - th;

  geom->g = bl_gdet_func(r, th);
  bl_gcov_func(r, th, geom->gcov);
  bl_gcon_func(r, th, geom->gcon);

}

// find the determinant of the bl metric
FTYPE bl_gdet_func(FTYPE r, FTYPE th)
{
  FTYPE a2, r2;

  a2 = a * a;
  r2 = r * r;
  return (bl_gdet);
}

// find gcov for bl metric
void bl_gcov_func(FTYPE r, FTYPE th, FTYPE gcov[][NDIM])
{
  int j, k;
  FTYPE sth, cth, s2, a2, r2, r3, DD, mu;

  DLOOP gcov[j][k] = 0.;

  sth = sin(th);
  s2 = sth * sth;
  cth = cos(th);
  a2 = a * a;
  r2 = r * r;
  r3 = r2 * r;
  DD = 1. - 2. / r + a2 / r2;
  mu = 1. + a2 * cth * cth / r2;

  gcov[TT][TT] = bl_gcov00;
  gcov[TT][3] = bl_gcov03;
  gcov[1][1] = bl_gcov11;
  gcov[2][2] = bl_gcov22;
  gcov[3][TT] = bl_gcov30;
  gcov[3][3] = bl_gcov33;

}

// find gcon for bl metric
void bl_gcon_func(FTYPE r, FTYPE th, FTYPE gcon[][NDIM])
{
  int j, k;
  FTYPE sth, cth, a2, r2, r3, DD, mu;

  DLOOP gcon[j][k] = 0.;

  if(POSDEFMETRIC){
    sth = fabs(sin(th));
  }
  else{
    sth = sin(th);
  }
  if (fabs(sth) < SMALL) {
    if(sth>=0) sth=SMALL;
    if(sth<0) sth=-SMALL;
  }
  
  cth = cos(th);
  a2 = a * a;
  r2 = r * r;
  r3 = r2 * r;
  DD = 1. - 2. / r + a2 / r2;
  mu = 1. + a2 * cth * cth / r2;

  gcon[TT][TT] = bl_gcon00;
  gcon[TT][3] = bl_gcon03;
  gcon[1][1] = bl_gcon11;
  gcon[2][2] = bl_gcon22;
  gcon[3][TT] = bl_gcon30;
  gcon[3][3] = bl_gcon33;

}

// ///////////////////////////////////////////////////////////
// 
// below are independent of user choice of metric/coords/grid
// 
// ///////////////////////////////////////////////////////////


// find determinant in general of a metric
/* assumes gcov has been set first; returns determinant */
FTYPE gdet_func(FTYPE gcov[][NDIM])
{
  static int firstc = 1;
  static FTYPE **tmp;
  FTYPE d;
  int j, k, indx[NDIM];

  if (firstc) {
    tmp = dmatrix(1, NDIM, 1, NDIM);
    firstc = 0;
  }

  DLOOP tmp[j + 1][k + 1] = gcov[j][k];
  ludcmp(tmp, NDIM, indx - 1, &d);
  // below from 1..NDIM due to ludcmp requiring 1..N
  for (j = 1; j <= NDIM; j++)
    d *= tmp[j][j];

  return (sqrt(fabs(d)));
}

/* invert gcov to get gcon */
void gcon_func(FTYPE gcov[][NDIM], FTYPE gcon[][NDIM])
{
  static int firstc = 1;
  int j, k;
  static FTYPE **tmp;

  if (firstc) {
    tmp = dmatrix(1, NDIM, 1, NDIM);
    firstc = 0;
  }

  DLOOP tmp[j + 1][k + 1] = gcov[j][k];
  gaussj(tmp, NDIM, NULL, 0);
  DLOOP gcon[j][k] = tmp[k + 1][j + 1];
}

/* 
   this gives the connection coefficient \Gamma^{i}_{j,k} =
   conn[..][i][j][k] where i = {0,1,2,3} corresponds to {t,r,theta,phi} 
 */

#define DELTA 1.e-10

/* NOTE: parameter hides global variable */
void conn_func(FTYPE *X, struct of_geom *geom,
	       FTYPE conn[][NDIM][NDIM])
{
  int i, j, k, l;
  FTYPE tmp[NDIM][NDIM][NDIM];
  FTYPE Xh[NDIM], Xl[NDIM];
  FTYPE gh[NDIM][NDIM];
  FTYPE gl[NDIM][NDIM];

  for (k = 0; k < NDIM; k++) {
    for (l = 0; l < NDIM; l++)
      Xh[l] = X[l];
    for (l = 0; l < NDIM; l++)
      Xl[l] = X[l];
    Xh[k] += DELTA;
    Xl[k] -= DELTA;
    gcov_func(Xh, gh);
    gcov_func(Xl, gl);

    for (i = 0; i < NDIM; i++)
      for (j = 0; j < NDIM; j++)
	conn[i][j][k] = (gh[i][j] - gl[i][j]) / (Xh[k] - Xl[k]);
  }

  /* now rearrange to find \Gamma_{ijk} */
  for (i = 0; i < NDIM; i++)
    for (j = 0; j < NDIM; j++)
      for (k = 0; k < NDIM; k++)
	tmp[i][j][k] =
	    0.5 * (conn[j][i][k] + conn[k][i][j] - conn[k][j][i]);

  /* finally, raise index */
  for (i = 0; i < NDIM; i++)
    for (j = 0; j < NDIM; j++)
      for (k = 0; k < NDIM; k++) {
	conn[i][j][k] = 0.;
	for (l = 0; l < NDIM; l++)
	  conn[i][j][k] += geom->gcon[i][l] * tmp[l][j][k];
      }

  /* done! */
}

#undef DELTA



/* 
   FTYPE delta(int i, int j) { if(i == j) return(1.) ; else return(0.) 
   ; } */

/* Minkowski metric; signature +2 */
/* 
   FTYPE mink(int i, int j) { if(i == j) { if(i == 0) return(-1.) ;
   else return(1.) ; } else return(0.) ; } */
