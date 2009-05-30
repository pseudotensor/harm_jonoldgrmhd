
#include "decs.h"

/** 
 *
 * this file contains all the coordinate dependent
 * parts of the code, except the initial and boundary
 * conditions 
 *
 **/

/* Returns boyer-lindquist coordinte of point */
void bl_coord(FTYPE *X, FTYPE *r, FTYPE *th)
{
  if (defcoord == 0) {
    //    *r = Rin*exp(X[1]) ;
    *r = R0+exp(X[1]) ;
    //*r = Rin * exp(X[1]);
    *th = M_PI * X[2] + ((1. - hslope) / 2.) * sin(2. * M_PI * X[2]);
  } else if (defcoord == 1) {
    *r = R0+Rin * exp(X[1] * log(Rout / Rin));
    *th =
	((-49. * hslope + 60. * M_PI) * X[2]) / 12. +
	((247. * hslope - 240. * M_PI) * pow(X[2],
					     2)) / 12. +
	((-83. * hslope + 80. * M_PI) * pow(X[2],
					    3)) / 2. -
	(5. * (-25. * hslope + 24. * M_PI) * pow(X[2], 4)) / 3. +
	(2. * (-25. * hslope + 24. * M_PI) * pow(X[2], 5)) / 3.;
  }
}



// Jacobian for dx uniform per dx nonuniform (dx/dr / dx/dr')

void dxdxprim(FTYPE *X, FTYPE r, FTYPE th, FTYPE *dxdxp)
{

  if (defcoord == 0) {
    dxdxp[0] = 1.;
    dxdxp[1] = r-R0;
    dxdxp[2] = M_PI + (1. - hslope) * M_PI * cos(2. * M_PI * X[2]);
    dxdxp[3] = 1.;
  } else if (defcoord == 1) {
    dxdxp[0] = 1.;
    dxdxp[1] = (r-R0) * log(Rout / Rin);

    dxdxp[2] = (-49. * hslope + 60. * M_PI) / 12. +
	((247. * hslope - 240. * M_PI) * X[2]) / 6. +
	(3. * (-83. * hslope + 80. * M_PI) * pow(X[2], 2)) / 2. -
	(20. * (-25. * hslope + 24. * M_PI) * pow(X[2], 3)) / 3. +
	(10. * (-25. * hslope + 24. * M_PI) * pow(X[2], 4)) / 3.;

    dxdxp[3] = 1.;
  }

}

// /////////////////////////////////////////////////////////////////
// 
// Below set X uniform grid -- usually doesn't change.
// Can usually force startx[1]=startx[2]=0. and dx[1]=1./N1 dx[2]=1./N2
// 
// /////////////////////////////////////////////////////////////////


/* some grid location, dxs */
void set_points()
{
  if (defcoord == 0) {
    startx[1] = log(Rin-R0);
    startx[2] = 0.;
    dx[1] = log((Rout-R0)/(Rin-R0)) / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 2.0*M_PI;
  } else if (defcoord == 1) {
    startx[1] = 0.;
    startx[2] = 0.;
    dx[1] = 1. / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 2.0*M_PI;
  }

}


void coord(int i, int j, int loc, FTYPE *X)
{
  if (loc == FACE1) {
    X[1] = startx[1] + (i + startpos[1]) * dx[1];
    X[2] = startx[2] + ((j + startpos[2]) + 0.5) * dx[2];
  } else if (loc == FACE2) {
    X[1] = startx[1] + ((i + startpos[1]) + 0.5) * dx[1];
    X[2] = startx[2] + (j + startpos[2]) * dx[2];
  } else if (loc == CENT) {
    X[1] = startx[1] + ((i + startpos[1]) + 0.5) * dx[1];
    X[2] = startx[2] + ((j + startpos[2]) + 0.5) * dx[2];
  } else {
    X[1] = startx[1] + (i + startpos[1]) * dx[1];
    X[2] = startx[2] + (j + startpos[2]) * dx[2];
  }

  return;
}
