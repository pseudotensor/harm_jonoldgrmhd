
#include "decs.h"

/* apply floors to density, internal energy */

// currently called before bound, which assumes bound sets boundary
// values exactly as wanted without any fixing.

#if 0
void fixup(FTYPE (*pv)[NH + 4][NPR])
{
  int i, j;

  ZSLOOP(-2, NR + 1, -2, NH + 1) pfixup(pv[i][j], i, j);

}
#else

void fixup(FTYPE (*pv)[N2 + 4][NPR])
{
  int i, j, k;
  int ip, jp, im, jm;
  FTYPE bsq, del;
  FTYPE r, th, X[NDIM];
  FTYPE uuscal, uurscal, uutscal, rhoscal, rhorscal, rhotscal;
  FTYPE ftempA,ftempB;
  struct of_state q;
  struct of_geom geom;
  FTYPE prfloor[NPR];
  FTYPE prafter[NPR];
  FTYPE fluxbefore[NPR];
  FTYPE fluxafter[NPR];
  int checkfl[NPR];

  // whether to check floor condition
  checkfl[0]=1;
  checkfl[1]=1;
  checkfl[2]=0;
  checkfl[3]=0;
  checkfl[4]=0;
  checkfl[5]=0;
  checkfl[6]=0;
  checkfl[7]=0;


  //  ZSLOOP(-2,N1+1,-2,N2+1) {
  ZSLOOP(0, N1 - 1, 0, N2 - 1) {
    coord(i, j, CENT, X);
    bl_coord(X, &r, &th);

    ////////////////////
    // scaling functions
    ////
    if(1){ // choice
      rhorscal = pow(r, -1.5);
    }
    else{
      rhorscal = pow(r, -2.0);
    }
    uurscal = rhorscal / r;


    if(0){ // choice
      ftempA=(RHOMAX+RHOMIN)*0.5/RHOMIN;
      ftempB=(RHOMAX-RHOMIN)*0.5/RHOMIN;
      // choice
      if(1) rhotscal = ftempA+ftempB*cos(2.0*M_PI*X[2]);
      else  rhotscal = ftempA+ftempB*cos(2.0*M_PI*th);
      ftempA=(UUMAX+UUMIN)*0.5/UUMIN;
      // choice (make same as above)
      ftempB=(UUMAX-UUMIN)*0.5/UUMIN;
      if(1) uutscal = ftempA+ftempB*cos(2.0*M_PI*X[2]);
      else  uutscal = ftempA+ftempB*cos(2.0*M_PI*th);
    }
    else{
      rhotscal = 1.0;
      uutscal = 1.0;
    }


    uuscal = UUMIN*uurscal*uutscal;
    rhoscal = RHOMIN*rhorscal*rhotscal;

    // assign general floor variables
    prfloor[0]=rhoscal;
    prfloor[1]=uuscal;
    prfloor[2]=pv[i][j][2];
    prfloor[3]=pv[i][j][3];
    prfloor[4]=pv[i][j][4];
    prfloor[5]=pv[i][j][5];
    prfloor[6]=pv[i][j][6];
    prfloor[7]=pv[i][j][7];


    // add in 10X more than the floor (try to avoid instability in reaction with
    // floor, since setting back right to floor can lead to runaway
    // feedback effect)
    prafter[0]=prfloor[0]*1.0;
    prafter[1]=prfloor[1]*1.0;
    prafter[2]=prfloor[2];
    prafter[3]=prfloor[3];
    prafter[4]=prfloor[4];
    prafter[5]=prfloor[5];
    prafter[6]=prfloor[6];
    prafter[7]=prfloor[7];

    get_geometry(i,j,CENT,&geom);
    get_state(pv[i][j],&geom,&q);
    primtoflux(pv[i][j],&q,0,&geom,fluxbefore);
    get_state(prafter,&geom,&q);
    primtoflux(prafter,&q,0,&geom,fluxafter);
    // shouldn't fail since before and after states should be ok, as
    // long as would have changed the value.  Check will occur if
    // simulation continues ok.  Could place check inside if below.
    PLOOP{
      if ( checkfl[k]&&(prfloor[k] - pv[i][j][k] > 0) ){
	fladd[k] += dVF * (fluxafter[k]-fluxbefore[k]);
	// MARK
	pv[i][j][k] = prafter[k];
      }
    }
    /* 
       if(bsq_calc(pv[i][j],&geom,&bsq)>=1){
       fprintf(fail_file,"1init:bsq_calc: failure\n");
       fflush(fail_file); return(1); } if(bsq >
       MAXBSQOVERRHO*pv[i][j][RHO]){
       fladd[RHO]+=dVF*gdet[i][j][CENT]*(rhoscal-pv[i][j][RHO]);
       pv[i][j][RHO] = bsq/MAXBSQOVERRHO ; } if(bsq >
       MAXBSQOVERUU*pv[i][j][UU]){
       fladd[UU]+=dVF*gdet[i][j][CENT]*(uuscal-pv[i][j][UU]);
       pv[i][j][UU] = bsq/MAXBSQOVERUU ; } */
  }
}


#endif
