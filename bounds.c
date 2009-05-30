
#include "decs.h"

/* bound array containing entire set of primitive variables */


int bound_prim(FTYPE prim[][N2 + 4][NPR])
{
  int i, j, k;
  FTYPE X[NDIM];
  FTYPE r,th;
  struct of_geom geom;

  if (mycpupos[1] == 0) {
    /* inner r boundary condition: u, gdet extrapolation */
    for (j = 0; j < N2; j++) {
      prim[-1][j][RHO] = prim[0][j][RHO] *
	  gdet[0][j][CENT] / gdet[-1][j][CENT];
      prim[-2][j][RHO] = prim[0][j][RHO] *
	  gdet[0][j][CENT] / gdet[-2][j][CENT];

      prim[-1][j][UU] = prim[0][j][UU] *
	  gdet[0][j][CENT] / gdet[-1][j][CENT];
      prim[-2][j][UU] = prim[0][j][UU] *
	  gdet[0][j][CENT] / gdet[-2][j][CENT];

      prim[-1][j][U1] = prim[0][j][U1] * (1. + 1. * dx[1]);
      prim[-2][j][U1] = prim[0][j][U1] * (1. + 2. * dx[1]);

      prim[-1][j][U2] = prim[0][j][U2] * (1. - 1. * dx[1]);
      prim[-2][j][U2] = prim[0][j][U2] * (1. - 2. * dx[1]);

      prim[-1][j][U3] = prim[0][j][U3] * (1. - 1. * dx[1]);
      prim[-2][j][U3] = prim[0][j][U3] * (1. - 2. * dx[1]);

      prim[-1][j][B1] = prim[0][j][B1] *
	  gdet[0][j][CENT] / gdet[-1][j][CENT];
      prim[-2][j][B1] = prim[0][j][B1] *
	  gdet[0][j][CENT] / gdet[-2][j][CENT];

      prim[-1][j][B2] = prim[0][j][B2] * (1. - 1. * dx[1]);
      prim[-2][j][B2] = prim[0][j][B2] * (1. - 2. * dx[1]);

      prim[-1][j][B3] = prim[0][j][B3] * (1. - 1. * dx[1]);
      prim[-2][j][B3] = prim[0][j][B3] * (1. - 2. * dx[1]);
    }
  }

  // outer r BC:
  if (mycpupos[1] == ncpux1 - 1) {
    if(BCtype[X1UP]==OUTFLOW){
  /* outer r BC: outflow */

  /* 
     for(j=0;j<N2;j++) PLOOP { prim[N1][j][k] = prim[N1-1][j][k] ;
     prim[N1+1][j][k] = prim[N1-1][j][k] ; } */
    for (j = 0; j < N2; j++) {
      prim[N1][j][RHO] = prim[N1 - 1][j][RHO] *
	  gdet[N1 - 1][j][CENT] / gdet[N1][j][CENT];
      prim[N1 + 1][j][RHO] = prim[N1 - 1][j][RHO] *
	  gdet[N1 - 1][j][CENT] / gdet[N1 + 1][j][CENT];

      prim[N1][j][UU] = prim[N1 - 1][j][UU] *
	  gdet[N1 - 1][j][CENT] / gdet[N1][j][CENT];
      prim[N1 + 1][j][UU] = prim[N1 - 1][j][UU] *
	  gdet[N1 - 1][j][CENT] / gdet[N1 + 1][j][CENT];

      prim[N1][j][U1] = prim[N1 - 1][j][U1] * (1. - 2. * dx[1]);
      prim[N1 + 1][j][U1] = prim[N1 - 1][j][U1] * (1. - 4. * dx[1]);

      prim[N1][j][U2] = prim[N1 - 1][j][U2] * (1. - 1. * dx[1]);
      prim[N1 + 1][j][U2] = prim[N1 - 1][j][U2] * (1. - 2. * dx[1]);

      prim[N1][j][U3] = prim[N1 - 1][j][U3] * (1. - 1. * dx[1]);
      prim[N1 + 1][j][U3] = prim[N1 - 1][j][U3] * (1. - 2. * dx[1]);

      prim[N1][j][B1] = prim[N1 - 1][j][B1] *
	  gdet[N1 - 1][j][CENT] / gdet[N1][j][CENT];
      prim[N1 + 1][j][B1] = prim[N1 - 1][j][B1] *
	  gdet[N1 - 1][j][CENT] / gdet[N1 + 1][j][CENT];

      prim[N1][j][B2] = prim[N1 - 1][j][B2] * (1. - 1. * dx[1]);
      prim[N1 + 1][j][B2] = prim[N1 - 1][j][B2] * (1. - 2. * dx[1]);

      prim[N1][j][B3] = prim[N1 - 1][j][B3] * (1. - 1. * dx[1]);
      prim[N1 + 1][j][B3] = prim[N1 - 1][j][B3] * (1. - 2. * dx[1]);
    }
    }
  /* if fixed BC: do nothing */
  }

  /* inner polar BC */
  if (mycpupos[2] == 0) {
    for (i = -2; i <= N1 + 1; i++)
      PLOOP {
      prim[i][-1][k] = prim[i][0][k];
      prim[i][-2][k] = prim[i][1][k];
      }
  }

  /* outer polar BC */
  if (mycpupos[2] == ncpux2 - 1) {
    for (i = -2; i <= N1 + 1; i++)
      PLOOP {
      prim[i][N2][k] = prim[i][N2 - 1][k];
      prim[i][N2 + 1][k] = prim[i][N2 - 2][k];
    }
  }
  if(BCtype[X1UP]==OUTFLOW){
  /* make sure there is no inflow at the outer radial boundary */
  if (mycpupos[1] == ncpux1 - 1) {
    for (j = -2; j < N2 + 2; j++)
      for (i = N1; i < N1 + 2; i++)
	if (prim[i][j][U1] < 0.)
	  prim[i][j][U1] = 0.;
  }
  }
  /* if fixed BC do nothing */

  /* make sure b and u are antisymmetric at the poles */
  /* inner pole */
  if (mycpupos[2] == 0) {
    if(POSDEFMETRIC==0){
      for (i = -2; i < N1 + 2; i++) for (j = -2; j < 0; j++) {
	prim[i][j][U2] *= -1.;
	prim[i][j][B2] *= -1.;
	prim[i][j][U3] *= -1.;
	prim[i][j][B3] *= -1.;
      }
    }
    else{
      for (i = -2; i < N1 + 2; i++) for (j = -2; j < 0; j++) {
	prim[i][j][U2] *= 1.;
	prim[i][j][B2] *= 1.;
	prim[i][j][U3] *= 1.;
	prim[i][j][B3] *= 1.;
      }
    }
  }
  /* outer pole */
  if (mycpupos[2] == ncpux2 - 1) {
    if(POSDEFMETRIC==0){
      for (i = -2; i < N1 + 2; i++) for (j = N2; j < N2 + 2; j++) {
	prim[i][j][U2] *= -1.;
	prim[i][j][B2] *= -1.;
	prim[i][j][U3] *= -1.;
	prim[i][j][B3] *= -1.;
      }
    }
    else{
      for (i = -2; i < N1 + 2; i++) for (j = N2; j < N2 + 2; j++) {
	prim[i][j][U2] *= 1.;
	prim[i][j][B2] *= 1.;
	prim[i][j][U3] *= 1.;
	prim[i][j][B3] *= 1.;
      }
    }
  }

  if (USEMPI) bound_mpi(prim);

  return (0);
}
