
#include "decs.h"

/* calculate fluxes in direction dir and conserved variable U; these
   are always needed together, so there is no point in calculated the
   stress tensor twice */

int primtoflux(FTYPE *pr, struct of_state *q, int dir,
	       struct of_geom *geom, FTYPE *flux)
{
  // sizes: NPR,struct of_state, int, struct of_geom, NPR
  int i = 0, j = 0, k = 0;
  FTYPE mhd[NDIM];

  /* particle number flux */
  flux[RHO] = pr[RHO] * q->ucon[dir];

  // GODMARK WTF!
  // if(mhd_calc(pr,dir,q,mhd)>=1)
  // FAILSTATEMENT("phys.c:primtoflux()","mhd_calc() dir=1or2",1);
  mhd_calc(pr, dir, q, mhd);

  /* MHD stress-energy tensor w/ first index up, second index down. */
  flux[UU] = mhd[0] + flux[RHO];
  flux[U1] = mhd[1];
  flux[U2] = mhd[2];
  flux[U3] = mhd[3];

  /* dual of Maxwell tensor */
  flux[B1] = q->bcon[1] * q->ucon[dir] - q->bcon[dir] * q->ucon[1];
  flux[B2] = q->bcon[2] * q->ucon[dir] - q->bcon[dir] * q->ucon[2];
  flux[B3] = q->bcon[3] * q->ucon[dir] - q->bcon[dir] * q->ucon[3];

  PLOOP flux[k] *= geom->g;

  return (0);
}

/* calculate "conserved" quantities */
int primtoU(FTYPE *pr, struct of_state *q, struct of_geom *geom,
	    FTYPE *U)
{
  int i = 0, j = 0, k = 0;
  if (primtoflux(pr, q, 0, geom, U) >= 1)
    FAILSTATEMENT("phys.c:primtoU()", "primtoflux_calc() dir=0", 1);

  return (0);
}


/* calculate magnetic field four-vector */
void bcon_calc(FTYPE *pr, FTYPE *ucon, FTYPE *ucov, FTYPE *bcon)
{
  int j;

  bcon[TT] = pr[B1] * ucov[1] + pr[B2] * ucov[2] + pr[B3] * ucov[3];
  for (j = 1; j < 4; j++)
    bcon[j] = (pr[B1 - 1 + j] + bcon[TT] * ucon[j]) / ucon[TT];

  return;
}

/* MHD stress tensor, with first index up, second index down */
void mhd_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *mhd)
{
  int j;
  FTYPE r, u, P, w, bsq, eta, ptot;

  r = pr[RHO];
  u = pr[UU];
  P = (gam - 1.) * u;
  w = P + r + u;
  bsq = dot(q->bcon, q->bcov);
  eta = w + bsq;
  ptot = P + bsq / 2.;

  /* single row of mhd stress tensor, first index up, second index down 
   */
  DLOOPA mhd[j] = eta * q->ucon[dir] * q->ucov[j]
      + ptot * delta(dir, j) - q->bcon[dir] * q->bcov[j];

}

/* add in source terms to equations of motion */
int source(FTYPE *ph, struct of_geom *ptrgeom, int ii, int jj,
	   FTYPE *dU)
{
  FTYPE mhd[NDIM][NDIM];
  int i = 0, j = 0, k = 0;
  struct of_state q;

  if (get_state(ph, ptrgeom, &q) >= 1)
    FAILSTATEMENT("phys.c:source()", "get_state() dir=0", 1);
  mhd_calc(ph, 0, &q, mhd[0]);
  mhd_calc(ph, 1, &q, mhd[1]);
  mhd_calc(ph, 2, &q, mhd[2]);
  mhd_calc(ph, 3, &q, mhd[3]);

  /* contract mhd stress tensor with connection */
  PLOOP dU[k] = 0.;
  DLOOP {
    dU[UU] += mhd[j][k] * conn[ii][jj][k][0][j];
    dU[U1] += mhd[j][k] * conn[ii][jj][k][1][j];
    dU[U2] += mhd[j][k] * conn[ii][jj][k][2][j];
    dU[U3] += mhd[j][k] * conn[ii][jj][k][3][j];
  }

  /* cooling */
  if(cooling){
    dU[UU] += coolfunc(h_over_r, ph, ptrgeom, &q);
  }
  PLOOP dU[k] *= ptrgeom->g;

  /* done! */
  return (0);
}

/* returns b^2 (i.e., twice magnetic pressure) */
int bsq_calc(FTYPE *pr, struct of_geom *ptrgeom, FTYPE *bsq)
{
  int i = 0, j = 0, k = 0;
  struct of_state q;

  if (get_state(pr, ptrgeom, &q) >= 1)
    FAILSTATEMENT("phys.c:bsq_calc()", "get_state() dir=0", 1);
  *bsq = dot(q.bcon, q.bcov);
  return (0);
}

void lower(FTYPE *ucon, struct of_geom *geom, FTYPE *ucov)
{

  ucov[0] = geom->gcov[0][0] * ucon[0]
      + geom->gcov[0][1] * ucon[1]
      + geom->gcov[0][2] * ucon[2]
      + geom->gcov[0][3] * ucon[3];
  ucov[1] = geom->gcov[1][0] * ucon[0]
      + geom->gcov[1][1] * ucon[1]
      + geom->gcov[1][2] * ucon[2]
      + geom->gcov[1][3] * ucon[3];
  ucov[2] = geom->gcov[2][0] * ucon[0]
      + geom->gcov[2][1] * ucon[1]
      + geom->gcov[2][2] * ucon[2]
      + geom->gcov[2][3] * ucon[3];
  ucov[3] = geom->gcov[3][0] * ucon[0]
      + geom->gcov[3][1] * ucon[1]
      + geom->gcov[3][2] * ucon[2]
      + geom->gcov[3][3] * ucon[3];

  return;
}

void raise(FTYPE *ucov, struct of_geom *geom, FTYPE *ucon)
{

  ucon[0] = geom->gcon[0][0] * ucov[0]
      + geom->gcon[0][1] * ucov[1]
      + geom->gcon[0][2] * ucov[2]
      + geom->gcon[0][3] * ucov[3];
  ucon[1] = geom->gcon[1][0] * ucov[0]
      + geom->gcon[1][1] * ucov[1]
      + geom->gcon[1][2] * ucov[2]
      + geom->gcon[1][3] * ucov[3];
  ucon[2] = geom->gcon[2][0] * ucov[0]
      + geom->gcon[2][1] * ucov[1]
      + geom->gcon[2][2] * ucov[2]
      + geom->gcon[2][3] * ucov[3];
  ucon[3] = geom->gcon[3][0] * ucov[0]
      + geom->gcon[3][1] * ucov[1]
      + geom->gcon[3][2] * ucov[2]
      + geom->gcon[3][3] * ucov[3];

  return;
}

/* find ucon, ucov, bcon, bcov from primitive variables */
int get_state(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q)
{
  int i = 0, j = 0, k = 0;

  /* get ucon */
  if (ucon_calc(pr, ptrgeom, q->ucon) >= 1)
    FAILSTATEMENT("phys.c:get_state()", "ucon_calc()", 1);
  lower(q->ucon, ptrgeom, q->ucov);
  bcon_calc(pr, q->ucon, q->ucov, q->bcon);
  lower(q->bcon, ptrgeom, q->bcov);

  return (0);
}

/* load local geometry into structure geom */
void get_geometry(int ii, int jj, int kk, struct of_geom *geom)
{
  int j, k;

  DLOOP geom->gcov[j][k] = gcov[ii][jj][kk][j][k];
  DLOOP geom->gcon[j][k] = gcon[ii][jj][kk][j][k];
  geom->g = gdet[ii][jj][kk];
  icurr = ii;
  jcurr = jj;
  pcurr = kk;
}

/* find contravariant four-velocity */
int ucon_calc(FTYPE *pr, struct of_geom *geom, FTYPE *ucon)
{
  FTYPE discr;

  ucon[0] = 1.;
  ucon[1] = pr[U1];
  ucon[2] = pr[U2];
  ucon[3] = pr[U3];

  discr = geom->gcov[0][0] * ucon[0] * ucon[0]
      + geom->gcov[1][1] * ucon[1] * ucon[1]
      + geom->gcov[2][2] * ucon[2] * ucon[2]
      + geom->gcov[3][3] * ucon[3] * ucon[3]
      + 2. * (geom->gcov[0][1] * ucon[0] * ucon[1]
	      + geom->gcov[0][2] * ucon[0] * ucon[2]
	      + geom->gcov[0][3] * ucon[0] * ucon[3]
	      + geom->gcov[1][2] * ucon[1] * ucon[2]
	      + geom->gcov[1][3] * ucon[1] * ucon[3]
	      + geom->gcov[2][3] * ucon[2] * ucon[3]);

  if (discr > 0.) {
    if (fail(FAIL_UTCALC_DISCR) >= 1)
      return (1);
  }

  ucon[TT] = 1. / sqrt(-discr);
  ucon[1] *= ucon[TT];
  ucon[2] *= ucon[TT];
  ucon[3] *= ucon[TT];

  return (0);
}
