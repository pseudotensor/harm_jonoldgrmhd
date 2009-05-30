#include "decs.h"


void pre_init(int argc, char *argv[])
{
  int dir,k;

  fail_file = log_file = logfull_file = stderr;

  fprintf(stderr, "begin: pre_init\n");
  fflush(stderr);

  GAMMIE=0;// whether in Gammie output types (not binary/text choice)


  RESTARTMODE=0;// whether restarting from rdump or not
  WHICHFILE=0; // see diag.c for dump_cnt and image_cnt starts

  DODIAGS=1; // whether to do diagnostics

  POSDEFMETRIC=0; // see metric.c and bounds.c 

  /** FIXUP PARAMETERS **/
  RHOMIN=1.e-4;
  UUMIN =1.e-6;
  RHOMAX=1.e-2;
  UUMAX =1.e-4;
  MAXBSQOVERRHO =1.e8;
  MAXBSQOVERUU  =1.e8;

  /* maximum increase in timestep */
  SAFE=1.3;
  cooling=0; // by default unless changed by init.c
  lim = MC;			/* used montonized central-difference
				   limiter */
  // lim = MINM;

  defcon = 1;
  // choice
  periodicx1=0;
  periodicx2=0;
  // choice
  binaryoutput=TEXTOUTPUT;
  sortedoutput=SORTED;

  DIRLOOP PLOOP{
    pdot[dir][k]=0;
    pcum[dir][k]=0;
  }
  PLOOP fladd[k] = 0;

  nstep = realnstep = 0;

  whocalleducon = 0;
  restartsteps[0] = 0;
  restartsteps[1] = 0;
  failuremode = 0;
  halftimep = 0;
  whichrestart = 0;
  // initialze CPUs
  init_mpi(argc, argv);
  // now logs are defined
  // Define arrays
  // do here since generic
  set_arrays();


  trifprintf("end: pre_init\n");
}


void post_init(void)
{
  trifprintf("begin: post_init\n");

  // in synch always here
  if (error_check()) {
    dualfprintf(fail_file, "error_check detected failure at main:1\n");
    dualfprintf(fail_file, "Bad initial conditions\n");
    myexit(1);
  }

  find_horizon();
  setflux();
  // GODMARK
  //cour*=0.1;

  // user defined parameters
  if(RESTARTMODE!=0){
    restartonfail=0; // whether we are restarting on failure or not
    setfailresponse(restartonfail);
  }

  trifprintf("end: post_init\n");
}


void find_horizon(void)
{
  int i, j, k;
  FTYPE r1, r2;
  FTYPE X[NDIM];
  FTYPE r, th;
  int horizoncpupos1, gotit;
  FTYPE horizonvalue;
  // called after grid is setup for all cpus


  trifprintf("begin: find_horizon ... ");

  // find cpu column that brackets the horizon and determine the
  // i-offset of horizon

  horizonvalue = 1.0 + sqrt(1.0 - a * a);
  horizoni = -100;
  gotit = 0;
  for (k = numprocs - 1; k >= 0; k--) {	// should get done by first row
    if (k == myid) {
      for (i = N1 - 1; i >= 0; i--) {
	j = N2 / 2;		// doesn't matter
	coord(i, j, CENT, X);
	bl_coord(X, &r1, &th);
	coord(i + 1, j, CENT, X);
	bl_coord(X, &r2, &th);
	if (fabs(r1 - horizonvalue) <= (r2 - r1)) {	// find horizon
	  horizoni = i;
	  horizoncpupos1 = mycpupos[1];
	  break;
	}
      }
    }
    if (numprocs > 0) {
#if(USEMPI)
      MPI_Bcast(&horizoni, 1, MPI_INT, k, MPI_COMM_WORLD);
      MPI_Bcast(&horizoncpupos1, 1, MPI_INT, k, MPI_COMM_WORLD);
#endif
    }
    if (horizoni >= 0)
      gotit = 1;		// can stop entire process
    if (mycpupos[1] != horizoncpupos1) {
      horizoni = -100;
    }				// reset if not right cpu group
    if (gotit)
      break;
  }
  trifprintf("horizoni: %d horizoncpupos1: %d\n", horizoni,
	     horizoncpupos1);
  // just a check
  dualfprintf(log_file,"horizoni: %d mycpupos[1]: %d horizoncpupos1: %d\n", horizoni, mycpupos[1], horizoncpupos1);

  trifprintf("end: find_horizon\n");
}


void setflux(void)
{
  int dir;
  
  // only 0 through N-1 mean do flux
  if(mycpupos[1]==0){
    doflux[X1DN]=0; // or horizoni
    trifprintf("proc: %d doing flux X1DN\n",myid);
  }
  else doflux[X1DN]=-100;
  if(mycpupos[1]==ncpux1-1){
    doflux[X1UP]=N1;
    trifprintf("proc: %d doing flux X1UP\n",myid);
  }
  else doflux[X1UP]=-100;
  if(mycpupos[2]==0){
    doflux[X2DN]=0;
    trifprintf("proc: %d doing flux X2DN\n",myid);
  }
  else doflux[X2DN]=-100;
  if(mycpupos[2]==ncpux2-1){
    doflux[X2UP]=N2;
    trifprintf("proc: %d doing flux X2UP\n",myid);
  }
  else doflux[X2UP]=-100;
  // fluxes are on edges of zone, so 0 and N are on edge fluxes

}
