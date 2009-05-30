#include "global.h"
#include "mpidefs.h"

FTYPE a_p[N1 + 4][N2 + 4][NPR];	/* space for primitive vars */
FTYPE a_dq[N1 + 4][N2 + 4][NPR];	/* slopes */
FTYPE a_F1[N1 + 4][N2 + 4][NPR];	/* fluxes */
FTYPE a_F2[N1 + 4][N2 + 4][NPR];	/* fluxes */
FTYPE a_ph[N1 + 4][N2 + 4][NPR];	/* half-step primitives */

/* for debug */
FTYPE ivar[N1][N2][NPR];
short stat[N1][N2];
FTYPE psave[N1][N2][NPR];

/* grid functions */
FTYPE a_conn[N1 + 4][N2 + 4][NDIM][NDIM][NDIM];
FTYPE a_gcon[N1 + 4][N2 + 4][NPG][NDIM][NDIM];
FTYPE a_gcov[N1 + 4][N2 + 4][NPG][NDIM][NDIM];
FTYPE a_gdet[N1 + 4][N2 + 4][NPG];

FTYPE (*p)[N2 + 4][NPR];
FTYPE (*dq)[N2 + 4][NPR];
FTYPE (*F1)[N2 + 4][NPR];
FTYPE (*F2)[N2 + 4][NPR];
FTYPE (*ph)[N2 + 4][NPR];
FTYPE (*conn)[N2 + 4][NDIM][NDIM][NDIM];
FTYPE (*gcon)[N2 + 4][NPG][NDIM][NDIM];
FTYPE (*gcov)[N2 + 4][NPG][NDIM][NDIM];
FTYPE (*gdet)[N2 + 4][NPG];

/** GLOBAL PARAMETERS SECTION **/

/* physics parameters */
FTYPE a;
FTYPE gam;

/* numerical parameters */
int defcoord;
FTYPE Rin, R0, Rout, hslope;
FTYPE cour;
FTYPE dV, dVF, dx[NDIM], startx[NDIM];
SFTYPE dt,t,tf;
FTYPE rcurr, hcurr;
int istart, istop, jstart, jstop;
FTYPE mydminarg1, mydminarg2;
long nstep;

/* output parameters */
SFTYPE DTd;
SFTYPE DTener;
SFTYPE DTi;
long DTr;
long dump_cnt;
long image_cnt;
long rdump_cnt;
int nstroke;

/* global flags */
int failed;
int lim;
FTYPE defcon;

/* diagnostics */
SFTYPE pdot[COMPDIM*2][NPR];
SFTYPE pcum[COMPDIM*2][NPR];
SFTYPE frdot[N1][NPR];
int doflux[COMPDIM*2];
SFTYPE fladd[NPR];

/* current local position */
int icurr, jcurr, pcurr, ihere, jhere, phere;

/* Jon's addition */
int horizoni;
FTYPE tstart;
long realnstep;
int mpicombine;
int halftimep;
int whichrestart;
int appendold;
int whocalleducon;
// global flags
int failuremode;
int restartsteps[2];
int binaryoutput,sortedoutput;
long steptofaildump,steptofailmap;
int ifail,jfail,dofailmap,dofaildump,restartonfail;
// IC
FTYPE h_over_r;
// BC
int BCtype[COMPDIM*2];
int cooling;
int GAMMIE,DODIAGS,RESTARTMODE,WHICHFILE,POSDEFMETRIC;
FTYPE RHOMIN,UUMIN,RHOMAX,UUMAX,MAXBSQOVERRHO,MAXBSQOVERUU;
FTYPE SAFE;
