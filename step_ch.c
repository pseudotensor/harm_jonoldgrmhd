
/**
 *
 * this contains the generic piece of code for advancing
 * the primitive variables 
 *
**/


#include "decs.h"

/** algorithmic choices **/

/* use local lax-friedrichs or HLL flux */
// choice
#define HLLF	0.0		/* OK for nonrotating, or rotating up
				   to 0.5 */
#define GAMLAXF	0.0		/* use for more rapidly rotating */
#define JONLAXF	0.0		/* use for more rapidly rotating */
#define TVDLF	1		/* use for more rapidly rotating */
#define FTCS    0.0


/* JON: don't change this. */
/* these are different ways of calculating the EMFs */
#define FLUXCTTOTH 1
#define FLUXCTHLL  0

#define UTOPRIMVERSION 0
// 0: original gammie
// 1: ldz

#define TIMEORDER 2
// order of algorithm in time from 1 to 2.

/* whether to state what step has been reached */
#define DEBUG 2

/** end algorithmic choices **/


int step_ch()
{
  FTYPE ndt;
  int i, j, k;



  if(TIMEORDER==2){
    trifprintf("h");
    if (advance(p, p, 0.5 * dt, ph, &ndt) >= 1)
      FAILSTATEMENT("step_ch.c:step_ch()", "advance()", 1);
    trifprintf( "f");
    if (advance(p, ph, dt, p, &ndt) >= 1)
      FAILSTATEMENT("step_ch.c:step_ch()", "advance()", 2);
  }
  
  if(TIMEORDER==1){
    trifprintf( "f");
    //    if (advance(p, p, dt, p, &ndt) >= 1)
    if (advance(p, p, dt, ph, &ndt) >= 1)
      FAILSTATEMENT("step_ch.c:step_ch()", "advance()", 3);
    ZSLOOP(-2,N1+1,-2,N2+1) PLOOP p[i][j][k]=ph[i][j][k];
  }

  /* check timestep */
  if (dt < 1.e-9) {
    trifprintf( "timestep too small\n");
    myexit(11);
  }

  /* increment time */
  t += dt;
  realnstep++;

  /* set next timestep */
  // find global minimum value of ndt over all cpus
  mpifmin(&ndt);
  if (ndt > SAFE * dt)
    ndt = SAFE * dt;
  dt = ndt;
  if (t + dt > tf)
    dt = tf - t;		/* but don't step beyond end of run */

  /* done! */
  return (0);
}

int advance(FTYPE pi[][N2 + 4][NPR],
	    FTYPE pb[][N2 + 4][NPR],
	    FTYPE Dt, FTYPE pf[][N2 + 4][NPR], FTYPE *ndt)
{
  int i, j, k;
  FTYPE ndt1, ndt2, U[NPR], dU[NPR];
  struct of_geom geom;
  struct of_state q;

  FTYPE Uo[NPR], po[NPR];

  /* needed for Utoprim as initial conditions */
  ZLOOP PLOOP pf[i][j][k] = pi[i][j][k];


  trifprintf( "0");
  if (fluxcalc(pb, F1, 1, &ndt1,Dt) >= 1)
    FAILSTATEMENT("step_ch.c:advance()", "fluxcalc", 1);
  if (fluxcalc(pb, F2, 2, &ndt2,Dt) >= 1)
    FAILSTATEMENT("step_ch.c:advance()", "fluxcalc", 2);
  flux_ct(F1, F2);


  /* evaluate diagnostics based on fluxes on second pass (Dt=dt)*/
  if(Dt==dt) diag_flux(F1, F2,Dt);

  trifprintf( "1");
  
  /** now update pi to pf **/
  ZLOOP {

    get_geometry(i, j, CENT, &geom);
    if (source(pb[i][j], &geom, i, j, dU) >= 1)
      FAILSTATEMENT("step_ch.c:advance()", "source", 1);

    if (get_state(pi[i][j], &geom, &q) >= 1)
      FAILSTATEMENT("step_ch.c:advance()", "get_state()", 1);
    if (primtoU(pi[i][j], &q, &geom, U) >= 1)
      FAILSTATEMENT("step_ch.c:advance()", "primtoU()", 1);

    PLOOP {
      U[k] += Dt * (-(F1[i + 1][j][k] - F1[i][j][k]) / dx[1]
		    - (F2[i][j + 1][k] - F2[i][j][k]) / dx[2]
		    + dU[k]
	  );
    }


#if(UTOPRIMVERSION==0)
    if (Utoprim(U, &geom, pf[i][j]) >= 1)
      FAILSTATEMENT("step_ch.c:advance()", "Utoprim", 2);
#elif(UTOPRIMVERSION==1)
    if (Utoprim_ldz(U, &geom, pf[i][j]) >= 1)
      FAILSTATEMENT("step_ch.c:advance()", "Utoprim", 2);
#endif

  }
  // must check before MPI operation (since asymmetries would
  // desynchronize cpus
  if (error_check())
    FAILSTATEMENT("step_ch.c", "error_check", 1);
  fixup(pf);
  bound_prim(pf);

  *ndt = defcon * 1. / (1. / ndt1 + 1. / ndt2);

  trifprintf( "2");

  return (0);
}

#define MAX(a,b) ((a) > (b) ? (a) : (b))

int fluxcalc(FTYPE pr[][N2 + 4][NPR],
	     FTYPE F[][N2 + 4][NPR], int dir, FTYPE *ndt,SFTYPE Dt)
{
  int i, j, k, idel, jdel, face;
  FTYPE p_l[NPR], p_r[NPR], F_l[NPR], F_r[NPR], U_l[NPR], U_r[NPR];
  FTYPE cmax_l, cmax_r, cmin_l, cmin_r, cmax, cmin, dtij;
  FTYPE ctop;
  struct of_geom geom;
  struct of_state state_l, state_r;
  FTYPE TVDLFdt,TVDLFdto2;

  if(TVDLF){
    if(Dt==dt){ TVDLFdt=1.0; TVDLFdto2=0.0; }
    else { TVDLFdt=0.0; TVDLFdto2=1.0; }
  }
  else{
    TVDLFdt=0.0; TVDLFdto2=0.0;
  }
  TVDLFdt=0.0; TVDLFdto2=1.0;

  if (dir == 1) {
    idel = 1;
    jdel = 0;
    face = FACE1;
  } else if (dir == 2) {
    idel = 0;
    jdel = 1;
    face = FACE2;
  } else {
    exit(10);
  }

	/** evaluate slopes of primitive variables **/
  ZSLOOP(-1, N1, -1, N2) PLOOP {
    dq[i][j][k] = slope_lim(pr[i - idel][j - jdel][k],
			    pr[i][j][k], pr[i + idel][j + jdel][k]
	);
  }


  *ndt = 1.e9;
  ZSLOOP(-jdel, N1, -idel, N2) {
    PLOOP {
#if(TVDLF==0)
      p_l[k] = pr[i - idel][j - jdel][k]
	+ 0.5 * dq[i - idel][j - jdel][k];
      p_r[k] = pr[i][j][k]
	- 0.5 * dq[i][j][k];
#else
      p_l[k] = pr[i - idel][j - jdel][k]
	+ 0.5 * dq[i - idel][j - jdel][k];
      p_r[k] = pr[i][j][k]
	- 0.5 * dq[i][j][k];
      /*
      p_l[k]=pr[i-idel][j-jdel][k];
      p_r[k]=pr[i][j][k];
      */
      //      p_l[k]=p_r[k]=0.5*(pr[i-idel][j-jdel][k]+pr[i][j][k]);
#endif
    }

#if(0)
    // inconsistent with leaving pr alone
    if((i==N1)&&(mycpupos[1]==ncpux1-1)&&(BCtype[X1UP]==OUTFLOW)&&(p_l[U1]<0)){ p_l[U1]=0; p_r[U1]=0; }
#endif

#if(TVDLF==0)
    get_geometry(i, j, face, &geom);
#endif

#if(TVDLF==1)
    get_geometry(i, j, face, &geom);
    //    get_geometry(i-idel, j-jdel, CENT, &geom);
#endif
    if (get_state(p_l, &geom, &state_l) >= 1)
      FAILSTATEMENT("step_ch.c:fluxcalc()", "get_state()", 1);
    if (primtoflux(p_l, &state_l, dir, &geom, F_l) >= 1)
      FAILSTATEMENT("step_ch.c:fluxcalc()",
		    "primtoflux_calc() dir=1/2 l", 1);
    if (primtoflux(p_l, &state_l, TT, &geom, U_l) >= 1)
      FAILSTATEMENT("step_ch.c:fluxcalc()", "primtoflux_calc() dir=l0",
		    1);
    if (vchar(p_l, &state_l, dir, &geom, &cmax_l, &cmin_l) >= 1)
      FAILSTATEMENT("step_ch.c:fluxcalc()", "vchar() dir=1or2", 1);

#if(TVDLF==1)
    //    get_geometry(i, j, CENT, &geom);
#endif

    if (get_state(p_r, &geom, &state_r) >= 1)
      FAILSTATEMENT("step_ch.c:fluxcalc()", "get_state()", 2);
    if (primtoflux(p_r, &state_r, dir, &geom, F_r) >= 1)
      FAILSTATEMENT("step_ch.c:fluxcalc()",
		    "primtoflux_calc() dir=1/2 r", 1);
    if (primtoflux(p_r, &state_r, TT, &geom, U_r) >= 1)
      FAILSTATEMENT("step_ch.c:fluxcalc()", "primtoflux_calc() dir=r0",
		    1);
     if (vchar(p_r, &state_r, dir, &geom, &cmax_r, &cmin_r) >= 1)
      FAILSTATEMENT("step_ch.c:fluxcalc()", "vchar() dir=1or2", 2);

    //    trifprintf("%d %d\n",i,j);
    //PLOOP trifprintf("%10.5g %10.5g %10.5g %10.5g\n",pr[i][j][k],p_l[k],p_r[k],dq[i][j][k]);
    //trifprintf("\n");


    cmax = fabs(MAX(MAX(0., cmax_l), cmax_r));
    cmin = fabs(MAX(MAX(0., -cmin_l), -cmin_r));
    ctop = MAX(cmax, cmin);


    PLOOP F[i][j][k] = 
      + HLLF * ((cmax * F_l[k] + cmin * F_r[k] - cmax * cmin * (U_r[k] - U_l[k])) / (cmax + cmin + SMALL) )
      + GAMLAXF * (0.5 * (F_l[k] + F_r[k] - ctop * (U_r[k] - U_l[k])) )
      + JONLAXF * (F_l[k]) // needs more terms (diffusion) in advance
      + TVDLFdto2 * (0.5 * (F_l[k] + F_r[k] - 2.0*ctop * (U_r[k] - U_l[k])) )
      + TVDLFdt * (0.5 * (F_l[k] + F_r[k]) )
      + FTCS * (F_l[k])
      ;
    // TVDLF as per VAC


    /* evaluate restriction on timestep */
    /* in the end, only the fluxcalc() call for dt, not 0.5dt, is used */
    cmax = MAX(cmax, cmin);
    dtij = cour * dx[dir] / cmax;
    if (dtij < *ndt)
      *ndt = dtij;
  }

  return (0);
}

#undef  MAX



void flux_ct(FTYPE F1[][N2 + 4][NPR], FTYPE F2[][N2 + 4][NPR])
{
  int i, j;
  static FTYPE emf[N1 + 1][N2 + 1];

  /* calculate EMFs */
  /* Toth approach: just average */
  ZSLOOP(0, N1, 0, N2) emf[i][j] =
      0.25 * (F1[i][j][B2] + F1[i][j - 1][B2]
	      - F2[i][j][B1] - F2[i - 1][j][B1]);

  /* rewrite EMFs as fluxes, after Toth */
  ZSLOOP(0, N1, 0, N2 - 1) {
    F1[i][j][B1] = 0.;
    F1[i][j][B2] = 0.5 * (emf[i][j] + emf[i][j + 1]);
  }
  ZSLOOP(0, N1 - 1, 0, N2) {
    F2[i][j][B1] = -0.5 * (emf[i][j] + emf[i + 1][j]);
    F2[i][j][B2] = 0.;
  }

}
