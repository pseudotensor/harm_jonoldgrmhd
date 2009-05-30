
#include "decs.h"

/* set up all grid functions */
void set_grid()
{
  int i, j, k, l, m;
  int ii, jj;
  FTYPE X[NDIM];
  struct of_geom geom;
  int loc;

  /* set up boundaries, steps in coordinate grid */
  set_points();
  dV = dx[1] * dx[2]; // computational volume
  dVF = dV * dx[3] ; // full 3d volume (used for diagnostics only)

  DLOOPA X[j] = 0.;

  ZSLOOP(-2, N1 + 1, -2, N2 + 1) {

    for (loc = NUMGRIDPOS - 1; loc >= 0; loc--) {

      coord(i, j, loc, X);
      gcov_func(X, gcov[i][j][loc]);
      gdet[i][j][loc] = gdet_func(gcov[i][j][loc]);
      gcon_func(gcov[i][j][loc], gcon[i][j][loc]);

      // check if near static limit since can't divide by ucon_calc
      if (fabs(gcon[i][j][loc][TT][TT]) < SLEPSILON) {
	dualfprintf(fail_file,
		"grid location too near static limit: %d %d\n", i, j);
	myexit(1);
      }
      if (loc == CENT) {
	get_geometry(i, j, loc, &geom);
	conn_func(X, &geom, conn[i][j]);
      }
    }
  }

  /* done! */
}
