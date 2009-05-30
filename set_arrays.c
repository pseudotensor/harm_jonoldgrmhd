
#include "decs.h"

void set_arrays()
{
  int i, j, k;

  p = (FTYPE (*)[N2 + 4][NPR]) (&(a_p[2][2][0]));
  dq = (FTYPE (*)[N2 + 4][NPR]) (&(a_dq[2][2][0]));
  F1 = (FTYPE (*)[N2 + 4][NPR]) (&(a_F1[2][2][0]));
  F2 = (FTYPE (*)[N2 + 4][NPR]) (&(a_F2[2][2][0]));
  ph = (FTYPE (*)[N2 + 4][NPR]) (&(a_ph[2][2][0]));

  /* everything must be initialized to zero */
  ZSLOOP(-2, N1 + 1, -2, N2 + 1) {
    PLOOP {
      p[i][j][k] = 0.;
      ph[i][j][k] = 0.;
      dq[i][j][k] = 0.;
      F1[i][j][k] = 0.;
      F2[i][j][k] = 0.;
    }
  }

  ZLOOP stat[i][j] = 1;

  /* grid functions */
  conn = (FTYPE (*)[N2 + 4][NDIM][NDIM][NDIM])
      (&(a_conn[2][2][0][0][0]));
  gcon = (FTYPE (*)[N2 + 4][NPG][NDIM][NDIM])
      (&(a_gcon[2][2][0][0][0]));
  gcov = (FTYPE (*)[N2 + 4][NPG][NDIM][NDIM])
      (&(a_gcov[2][2][0][0][0]));
  gdet = (FTYPE (*)[N2 + 4][NPG])
      (&(a_gdet[2][2][0]));

}
