#include "global.h"
#include "mpidecs.h"

extern FTYPE a_p[N1 + 4][N2 + 4][NPR];	/* space for primitive vars */
extern FTYPE a_dq[N1 + 4][N2 + 4][NPR];	/* slopes */
extern FTYPE a_F1[N1 + 4][N2 + 4][NPR];	/* fluxes */
extern FTYPE a_F2[N1 + 4][N2 + 4][NPR];	/* fluxes */
extern FTYPE a_ph[N1 + 4][N2 + 4][NPR];	/* half-step primitives */

/* for debug */
extern FTYPE ivar[N1][N2][NPR];
extern short stat[N1][N2];
extern FTYPE psave[N1][N2][NPR];

/* grid functions */
extern FTYPE a_conn[N1 + 4][N2 + 4][NDIM][NDIM][NDIM];
extern FTYPE a_gcon[N1 + 4][N2 + 4][NPG][NDIM][NDIM];
extern FTYPE a_gcov[N1 + 4][N2 + 4][NPG][NDIM][NDIM];
extern FTYPE a_gdet[N1 + 4][N2 + 4][NPG];

extern FTYPE (*p)[N2 + 4][NPR];
extern FTYPE (*dq)[N2 + 4][NPR];
extern FTYPE (*F1)[N2 + 4][NPR];
extern FTYPE (*F2)[N2 + 4][NPR];
extern FTYPE (*ph)[N2 + 4][NPR];
extern FTYPE (*conn)[N2 + 4][NDIM][NDIM][NDIM];
extern FTYPE (*gcon)[N2 + 4][NPG][NDIM][NDIM];
extern FTYPE (*gcov)[N2 + 4][NPG][NDIM][NDIM];
extern FTYPE (*gdet)[N2 + 4][NPG];

/** GLOBAL PARAMETERS SECTION **/

/* physics parameters */
extern FTYPE a;
extern FTYPE gam;

/* numerical parameters */
extern int defcoord;
extern FTYPE Rin, R0, Rout, hslope;
extern FTYPE cour;
extern FTYPE dV, dVF, dx[NDIM], startx[NDIM];
extern SFTYPE dt,t,tf;
extern FTYPE rcurr, hcurr;
extern int istart, istop, jstart, jstop;
extern FTYPE mydminarg1, mydminarg2;
extern long nstep;

/* output parameters */
extern SFTYPE DTd;
extern SFTYPE DTener;
extern SFTYPE DTi;
extern long DTr;
extern long dump_cnt;
extern long image_cnt;
extern long rdump_cnt;
extern int nstroke;

/* global flags */
extern int failed;
extern int lim;
extern FTYPE defcon;

/* diagnostics */
extern SFTYPE pdot[COMPDIM*2][NPR];
extern SFTYPE pcum[COMPDIM*2][NPR];
extern SFTYPE frdot[N1][NPR];
extern int doflux[COMPDIM*2];
extern SFTYPE fladd[NPR];

/* current local position */
extern int icurr, jcurr, pcurr, ihere, jhere, phere;

/* Jon's addition */
extern int horizoni;
extern FTYPE tstart;
extern long realnstep;
extern int mpicombine;
extern int halftimep;
extern int whichrestart;
extern int appendold;
extern int whocalleducon;
// global flags
extern int failuremode;
extern int restartsteps[2];
extern int binaryoutput,sortedoutput;
extern long steptofaildump,steptofailmap;
extern int ifail,jfail,dofailmap,dofaildump,restartonfail;
// IC
extern FTYPE h_over_r;
// BC
extern int BCtype[COMPDIM*2];
extern int cooling;
extern int GAMMIE,DODIAGS,RESTARTMODE,WHICHFILE,POSDEFMETRIC;
extern FTYPE RHOMIN,UUMIN,RHOMAX,UUMAX,MAXBSQOVERRHO,MAXBSQOVERUU;
extern FTYPE SAFE;
