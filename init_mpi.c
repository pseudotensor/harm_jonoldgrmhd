#include "decs.h"

/* 
   modifications to gammie code: 1) add mympi.h to decs.h 2) add
   mpidecs.h to decs.h 3) add mpidefs.h to defs.h 4) add init_mpi() to
   main.c:main() 5) add bound_mpi() to bounds.c:bound_prim() 6) modify
   set_grid.c for mpi (base *position* on global geometry instead of
   local) 7) modify diag.c for output of files to separate names for
   each cpu 8) "" for rest of fopen (restart.c, postmort.c, etc.) 9)
   modify makefile for MPIability 10) modify step_ch.c for timestep and 
   flux */

void init_mpi(int argc, char *argv[])
{

#if(USEMPI)
  // initialize MPI
  workbc =
      (FTYPE(*)[COMPDIM * 2][NPR * NBIGBND * NBIGSM]) (&(workbca[-1][0]
							 [0]));
  init_MPI(argc, argv);
#else
  ncpux1 = 1;
  ncpux2 = 1;
  myid = 0;			// defines single process run
  sprintf(myidtxt, "");
  numprocs = 1;
#endif
  if (USEMPI) {
    mpicombine = 1;
  } else
    mpicombine = 0;

  // always done
  init_genfiles(0);
  init_placeongrid();


  trifprintf("done with init_mpi()\n");  fflush(log_file);

}

void myargs(int argc, char *argv[])
{
  if(argc!=COMPDIM+1){
    if(myid==0){
      fprintf(stderr,"proc: %04d : Incorrect command line: argc: %d needed=%d, please specify:\n",myid,argc,COMPDIM+1);
      fprintf(stderr,"proc: %04d : mpirun <mpirunoptions> <progname> ncpux1 ncpux2\n",myid);
    }
    exit(1);
  }
  ncpux1=atoi(argv[1]);
  ncpux2=atoi(argv[2]);
}


void init_genfiles(int gopp)
{
  char temps[MAXFILENAME];
  char extension[MAXFILENAME];

  fprintf(stderr, "begin: init_genfiles ... ");
  fflush(stderr);

  if (gopp == 1) {
    strcpy(extension, PPEXT);
  } else if (gopp == 0) {
    strcpy(extension, OUTEXT);
  }
  // always have fail and general log open

  sprintf(temps, "%s0_fail%s%s", DATADIR, extension, myidtxt);


  if ((fail_file = fopen(temps, "wt")) == NULL) {
    fprintf(stderr, "fail: Cannot open: %s\n", temps);
    exit(1);
  }
  fprintf(stderr, "opened: %s\n", temps);
  sprintf(temps, "%s0_log%s%s", DATADIR, extension, myidtxt);

  if ((log_file = fopen(temps, "wt")) == NULL) {
    fprintf(stderr, "log: Cannot open: %s\n", temps);
    exit(1);
  }
  fprintf(stderr, "opened: %s\n", temps);
  fprintf(log_file, "fail_file: %d log_file: %d\n", fail_file,
	  log_file);
  fflush(log_file);
  if (myid == 0) {
    sprintf(temps, "%s0_logfull%s", DATADIR, extension);

    if ((logfull_file = fopen(temps, "wt")) == NULL) {
      fprintf(stderr, "logfull: Cannot open: %s\n", temps);
      exit(1);
    }
    fprintf(stderr, "opened: %s\n", temps);
    fprintf(logfull_file, "logfull_file: %d \n", logfull_file);
    fflush(logfull_file);
  }

  // ok now
  trifprintf("end: init_genfiles\n");
}



#if(USEMPI)

int init_MPI(int argc, char *argv[])
{

  fprintf(stderr, "begin: init_MPI\n");
  fflush(stderr);

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  sprintf(myidtxt, CPUTXT, myid);
  MPI_Get_processor_name(processor_name, &procnamelen);

  // currently INIT provides args to rest of processes
  myargs(argc,argv);

  if (MAXCPUS < numprocs) {
    fprintf(stderr,
	    "Must increase MAXCPUS in global.h, %d is too many\n",
	    numprocs);
    myexit(1);
  }

  myfprintf(stderr,
	    "numprocs=%d ncpux1=%d ncpux2=%d percpusize: N1=%d N2=%d\n",
	    numprocs, ncpux1, ncpux2, N1, N2);

  fprintf(stderr, "proc: %s on %s\n", myidtxt, processor_name);

  fprintf(stderr, "end: init_MPI\n");
  fflush(stderr);
  return (0);
}


#endif


void init_placeongrid(void)
{
  int i, j, m, l;
  int N[COMPDIM + 1];
  int numbercpu[COMPDIM + 1];
  int dir;
  int opp[COMPDIM*2];

  trifprintf("begin: init_placeongrid ... ");

  N[1] = N1;
  N[2] = N2;

  numbercpu[1] = ncpux1;
  numbercpu[2] = ncpux2;

  mycpupos[1] = myid % ncpux1;
  mycpupos[2] = (int) ((myid % (numprocs)) / ncpux1);

  for (m = 1; m <= COMPDIM; m++) {
    startpos[m] = mycpupos[m] * N[m];
    endpos[m] = (mycpupos[m] + 1) * N[m] - 1;

    // add up sizes for total size of grid
    totalsize[m] = 0;
    itotalsize[m] = 0;
    for (i = 0; i < numbercpu[m]; i++) {
      totalsize[m] += N[m];
      itotalsize[m] += N[m];
    }
  }

  realtotalzones = totalzones =
      totalsize[1] * totalsize[2];
  itotalzones = itotalsize[1] * itotalsize[2];

  /////////////// standard interior MPI data transfer setup
  //
  for(dir=0;dir<COMPDIM*2;dir++) for(j=0;j<DIRNUMVARS;j++){
    dirset[dir][j]=0;
  }
  // see where this cpu needs to send/recv

  // figure out left/right send/recv
  if (mycpupos[1] > 0) {
    dirset[X1DN][DIRIF] = 1;		// do -x1 dir
  }
  if (mycpupos[1] < ncpux1 - 1) {
    dirset[X1UP][DIRIF] = 1;		// do +x1 dir
  }
  // figure out up/down send/recv
  if (mycpupos[2] > 0) {
    dirset[X2DN][DIRIF] = 1;		// -x2 dir
  }
  if (mycpupos[2] < ncpux2 - 1) {
    dirset[X2UP][DIRIF] = 1;		// towards and from +x2 dir
  }

  // only do periodic mpi if 
  if(periodicx1&&(ncpux1>1)){
    if(mycpupos[1]==0) dirset[X1DN][DIRIF]=1;
    else if(mycpupos[1]==ncpux1-1) dirset[X1UP][DIRIF]=1;
  }

  if(periodicx2&&(ncpux2>1)){
    if(mycpupos[2]==0) dirset[X2DN][DIRIF]=1;
    else if(mycpupos[2]==ncpux2-1) dirset[X2UP][DIRIF]=1;
  }

  opp[X1DN]=X1UP;
  opp[X1UP]=X1DN;
  opp[X2DN]=X2UP;
  opp[X2UP]=X2DN;


  
  // tags are defined by sender's ID and direction sent
  // tag=(myid*COMPDIM*2)+{0,1,2,3,4,5}
  // 0=right, 1=up,2=left,3=down,4=out,5=in
  // works/v bc[1=output/2=input][0,1,2,3,4,5]
  // so sends are like: (sendtoid,myid*COMPDIM*2+?) and recv's are
  // like: (fromid,otherid*COMPDIM*2+*) where ? and * are
  // opposites(i.e. 0 and 2, 1 and 3, 4 and 5)

  for(dir=0;dir<COMPDIM*2;dir++){

    if(dirset[dir][DIRIF]){
      if(dir==X1UP){ // right
	dirset[dir][DIRPSTART1]=N1-2;
	dirset[dir][DIRPSTOP1]=N1-1;
      }
      else if(dir==X1DN){ // left
	dirset[dir][DIRPSTART1]=0;
	dirset[dir][DIRPSTOP1]=1;
      }
      if((dir==X1UP)||(dir==X1DN)){
	dirset[dir][DIRPSTART2]=-N2BND;
	dirset[dir][DIRPSTOP2]=N2-1+N2BND;
      }

      if(dir==X2UP){ // up
	dirset[dir][DIRPSTART2]=N2-2;
	dirset[dir][DIRPSTOP2]=N2-1;
      }
      else if(dir==X2DN){ // down
	dirset[dir][DIRPSTART2]=0;
	dirset[dir][DIRPSTOP2]=1;
      }
      if((dir==X2UP)||(dir==X2DN)){
	dirset[dir][DIRPSTART1]=-N1BND;
	dirset[dir][DIRPSTOP1]=N1-1+N1BND;
      }
    
      dirset[dir][DIROPP]=opp[dir];

      if((dir==X1UP)||(dir==X1DN)) dirset[dir][DIRSIZE]=N1BND*N2M*NPR;
      else if((dir==X2UP)||(dir==X2DN)) dirset[dir][DIRSIZE]=N2BND*N1M*NPR;

      if((dir==X1UP)||(dir==X1DN)){
	if((periodicx1==0)||((mycpupos[1]>0)&&(mycpupos[1]<ncpux1-1))){
	  if(dir==X1UP) dirset[dir][DIROTHER]=myid+1;
	  if(dir==X1DN) dirset[dir][DIROTHER]=myid-1;
	}
	else if(periodicx1){
	  if(mycpupos[1]==0) dirset[dir][DIROTHER]=myid+(ncpux1-1);
	  else if(mycpupos[1]==ncpux1-1) dirset[dir][DIROTHER]=myid-(ncpux1-1);
	}
      }

      if((dir==X2UP)||(dir==X2DN)){
	if((periodicx2==0)||((mycpupos[2]>0)&&(mycpupos[2]<ncpux2-1))){
	  if(dir==X2UP) dirset[dir][DIROTHER]=myid+ncpux1;
	  if(dir==X2DN) dirset[dir][DIROTHER]=myid-ncpux1;
	}
	else if(periodicx2){
	  if(mycpupos[2]==0) dirset[dir][DIROTHER]=myid+(ncpux2-1)*ncpux1;
	  else if(mycpupos[2]==ncpux2-1) dirset[dir][DIROTHER]=myid-(ncpux2-1)*ncpux1;
	}
      }
      dirset[dir][DIRTAGS]= myid     * COMPDIM * 2 + dir;
      dirset[dir][DIRTAGR]= dirset[dir][DIROTHER] * COMPDIM * 2 + dirset[dir][DIROPP];
    
      if(dir==X1UP){ // right
	dirset[dir][DIRUSTART1]=N1;
	dirset[dir][DIRUSTOP1]=N1+1;
      }
      else if(dir==X1DN){ // left
	dirset[dir][DIRUSTART1]=-2;
	dirset[dir][DIRUSTOP1]=-1;
      }
      if((dir==X1UP)||(dir==X1DN)){
	dirset[dir][DIRUSTART2]=-N2BND;
	dirset[dir][DIRUSTOP2]=N2-1+N2BND;
      }

      if(dir==X2UP){ // up
	dirset[dir][DIRUSTART2]=N2;
	dirset[dir][DIRUSTOP2]=N2+1;
      }
      else if(dir==X2DN){ // down
	dirset[dir][DIRUSTART2]=-2;
	dirset[dir][DIRUSTOP2]=-1;
      }
      if((dir==X2UP)||(dir==X2DN)){
	dirset[dir][DIRUSTART1]=-N1BND;
	dirset[dir][DIRUSTOP1]=N1-1+N1BND;
    }
    }
  }


  fprintf(log_file,"per: %d %d\n", periodicx1, periodicx2);
  for (m = 1; m <= COMPDIM; m++) {
    fprintf(log_file,"mycpupos[%d]: %d\n", m, mycpupos[m]);
    fprintf(log_file, "startpos[%d]: %d\n", m, startpos[m]);
    fprintf(log_file, "endpos[%d]: %d\n", m, endpos[m]);
    fprintf(log_file, "totalsize[%d]: %d\n", m, totalsize[m]);
  }
  for (m = 0; m < COMPDIM*2; m++) {
    for(l = 0 ; l < DIRNUMVARS ; l++) {
      fprintf(log_file, "dirset[%d][%d]: %d\n", m, l, dirset[m][l]);
    }
  }
  trifprintf("totalzones: %d\n", totalzones);

  trifprintf("end: init_placeongrid\n");
}


int myexit(int call_code)
{
  int i, j, k, l;
  int cleanfinish;
  FILE *faildump;

  trifprintf("proc: %s : Exiting cc: %d nstep: %ld\n", myidtxt,
	  call_code, nstep);



  if (call_code >= 0) {
    if (fail_file)
      fclose(fail_file);
    if (log_file)
      fclose(log_file);
    myfclose(&logfull_file,"Can't close logfull_file\n");
  }
  if (call_code > 0) {
    fprintf(stderr,
	    "proc: %s : Failure.  Please check failure file: cc: %d\n",
	    myidtxt, call_code);

    cleanfinish = 1;
#if(USEMPI)
      // must abort since no clear to communicate to other cpus now
    //      MPI_Abort(MPI_COMM_WORLD, 1);
#endif
  } else
    cleanfinish = 1;

  if (cleanfinish) {
    fprintf(stderr, "proc: %s : dumping failure dump with callcode=2\n",
	    myidtxt);
      
    // assume want previous timestep data, not bad just-computed
    // data\n");
    ZLOOP PLOOP p[i][j][k] = ph[i][j][k];
    // now diag should not fail if last timestep was non-fail type
    if (DODIAGS)
      diag(2);

    fprintf(stderr,
	    "Ending Computation on proc: %s, holding for other cpus\n",
	    myidtxt);

#if(USEMPI)
    // finish up MPI
    MPI_Barrier(MPI_COMM_WORLD);	// required!
    MPI_Finalize();
#endif
    
    myfprintf(stderr, "Ended Computation on all processors\n");
  }    
  fprintf(stderr, "END\n");
  fflush(stderr);
  exit(0);
  return (0);
}

// note, this may be called in different locations of the code by
// different CPUs
int error_check(void)
{
  int i, j, k;
  int errorsend = 0;
  // check if error exists and exit if so

  if(failuremode==0) failuremode=failed;

  if (failuremode > 0) {
    dualfprintf(fail_file,
	    "Detected failure on proc: %d failuremode: %d nstep: %ld realnstep: %ld t: %15.10g\n",
	    myid, failuremode, nstep, realnstep, t);
  }

  if (numprocs > 1) {
    errorsend = failuremode;
#if(USEMPI)
    // fprintf(fail_file,"wtf: %d %d\n",errorsend,failuremode);
    // fflush(fail_file);
    MPI_Allreduce(&errorsend, &failuremode, 1, MPI_INT, MPI_MAX,
		  MPI_COMM_WORLD);
    // fprintf(fail_file,"wtf: %d %d\n",errorsend,failuremode);
    // fflush(fail_file);
#endif
  }
  if (failuremode > 0) {
    dualfprintf(fail_file,
	    "Result: Detected failure on proc: %d failuremode: %d nstep: %ld realnstep: %ld t: %15.10g\n",
	    myid, failuremode, nstep, realnstep, t);
    // control behavior of failure here (i.e. could return(1) and
    // continue or something)
    // if(failuremode==1) myexit(1);
    // if(failuremode==2) myexit(1);
    // if(failuremode==3) myexit(1);
    myexit(1);
    return (1);
  }
  return (0);
}


void mpiio_init(int bintxt, int sorted, FILE ** fpptr, int which,
		char *filename, int numcolumns,
		int datatype, void **jonioptr, void **writebufptr)
{

  logsfprintf("mpiio start init\n");

  // this check covers combine and seperate
  if(!sorted){
    if(bintxt==BINARYOUTPUT){
      dualfprintf(fail_file,"No such thing as binary unsorted output\n");
      myexit(1);
    }
  }

  mpiios_init(bintxt, sorted, fpptr, which, filename, numcolumns, datatype,
	      jonioptr, writebufptr);


  logsfprintf("mpiio end init\n");

}


void mpiio_combine(int bintxt, int sorted,
		   int numcolumns, int datatype,
		   FILE ** fpptr, void *jonio, void *writebuf)
{
  MPI_Datatype mpidt;

  logsfprintf("mpiio start combine\n");

  if(datatype==sizeof(double)) mpidt=MPI_DOUBLE;
  else if(datatype==sizeof(float)) mpidt=MPI_FLOAT;
  else if(datatype==sizeof(unsigned char)) mpidt=MPI_BYTE;


  if (sorted) {
    mpiios_combine(bintxt, mpidt, numcolumns, datatype, fpptr, jonio,
		  writebuf);
  }
  else{
    mpiiotu_combine(mpidt, numcolumns, datatype, fpptr, writebuf);
  }

  logsfprintf("mpiio end combine\n");

}


void mpiio_seperate(int bintxt, int sorted, int stage,
		    int numcolumns, int datatype,
		    FILE ** fpptr, void *jonio, void *writebuf)
{
  MPI_Datatype mpidt;

  logsfprintf("mpiio begin seperate\n");

  if(datatype==sizeof(double)) mpidt=MPI_DOUBLE;
  else if(datatype==sizeof(float)) mpidt=MPI_FLOAT;
  else if(datatype==sizeof(unsigned char)) mpidt=MPI_BYTE;

  if(sorted){
    mpiios_seperate(bintxt, stage, mpidt, numcolumns, datatype, fpptr, jonio, writebuf);
  }
  else{
    mpiiotu_seperate(stage, mpidt, numcolumns, datatype, fpptr, writebuf);
  }

  logsfprintf("mpiio end seperate\n");

}



// mpi io sorted
#if(USEMPI)
void mpiios_init(int bintxt, int sorted, FILE ** fp, int which, char *filename, int numcolumns,
		int datatype, void **jonioptr, void **writebufptr)
{
  // based on sizeof()
  // 1: unsigned char
  // 4: float
  // 8: double
  int sizeofmemory;
  double **jonio8;
  float **jonio4;
  unsigned char **jonio1;

  double **writebuf8;
  float **writebuf4;
  unsigned char **writebuf1;
  
  if(datatype==sizeof(double)){ jonio8=(double**)jonioptr; writebuf8=(double**)writebufptr; }
  else if(datatype==sizeof(float)){ jonio4=(float **)jonioptr; writebuf4=(float **)writebufptr; }
  else if(datatype==sizeof(unsigned char)){ jonio1=(unsigned char **)jonioptr; writebuf1=(unsigned char **)writebufptr; }
  
  logsfprintf("mpiios begin init\n");
  
  
  if ( (sorted==SORTED)&&(myid == 0) ){		// total on CPU=0
    if (which == WRITEFILE){
      if(bintxt==BINARYOUTPUT)      *fp = fopen(filename, "w");
      else if(bintxt==TEXTOUTPUT)      *fp = fopen(filename, "wt");
    }
    else if (which == READFILE){
      if(bintxt==BINARYOUTPUT) *fp = fopen(filename, "rb");
      else if(bintxt==TEXTOUTPUT) *fp = fopen(filename, "rt");
    }
    if (*fp == NULL) {
      fprintf(fail_file, "error opening file: %s\n", filename);
      myexit(2);
    }
    sizeofmemory = datatype * totalsize[1] * totalsize[2] * numcolumns;
    if(datatype==sizeof(double)) *jonio8 =(double*)malloc(sizeofmemory);
    else if(datatype==sizeof(float)) *jonio4=(float*)malloc(sizeofmemory);
    else if(datatype==sizeof(unsigned char)) *jonio1=(unsigned char*)malloc(sizeofmemory);
    if(
       (datatype==sizeof(double))&&(jonio8 == NULL) ||
       (datatype==sizeof(float))&&(jonio4 == NULL) ||
       (datatype==sizeof(unsigned char))&&(jonio1 == NULL)
       ){
      fprintf(fail_file, "Can't initialize jonio memory\n");
      myexit(1);
    }
  }
  sizeofmemory = datatype * N1 * N2 * numcolumns;
  if(datatype==sizeof(double)) *writebuf8 =(double*)malloc(sizeofmemory);
  else if(datatype==sizeof(float)) *writebuf4=(float*)malloc(sizeofmemory);
  else if(datatype==sizeof(unsigned char)) *writebuf1=(unsigned char*)malloc(sizeofmemory);
  if(
     (datatype==sizeof(double))&&(writebuf8 == NULL) ||
     (datatype==sizeof(float))&&(writebuf4 == NULL) ||
     (datatype==sizeof(unsigned char))&&(writebuf1 == NULL)
     ){
    
    
    fprintf(fail_file, "Can't initialize writebuf memory\n");
    myexit(1);
  }

  logsfprintf("mpiios begin init\n");

}


void mpiios_combine(int bintxt, MPI_Datatype mpidt, int numcolumns, int datatype,
		   FILE ** fp, void *jonio, void *writebuf)
{
  // based on sizeof()
  // 1: unsigned char
  // 4: float
  // 8: double
  int i, j, k, l, col, mapvaluejonio, mapvaluetempbuf;
#if(USEMPI)
  MPI_Request rrequest;
  MPI_Request srequest;
#endif
  int othercpupos[COMPDIM + 1];
  unsigned char *jonio1;
  float *jonio4;
  double *jonio8;
  unsigned char *writebuf1;
  float *writebuf4;
  double *writebuf8;

  logsfprintf("mpiios begin combine\n");


  if (datatype == sizeof(unsigned char))
    jonio1 = (unsigned char *) jonio;
  else if (datatype == sizeof(float))
    jonio4 = (float *) jonio;
  else if (datatype == sizeof(double))
    jonio8 = (double *) jonio;
  if (datatype == sizeof(unsigned char))
    writebuf1 = (unsigned char *) writebuf;
  else if (datatype == sizeof(float))
    writebuf4 = (float *) writebuf;
  else if (datatype == sizeof(double))
    writebuf8 = (double *) writebuf;


#if(USEMPI)
  // no need for tempbuf, works since first write to jonio is CPU=0's writebuf
  if(myid!=0) MPI_Isend(writebuf, N1 * N2 * numcolumns, mpidt, 0, myid,
			MPI_COMM_WORLD, &srequest);
  if (myid == 0) {
    for (l = 0; l < numprocs; l++) {
      if(l!=0){
	MPI_Irecv(writebuf, N1 * N2 * numcolumns, mpidt, l, l, MPI_COMM_WORLD, &rrequest);
	MPI_Wait(&rrequest, &mpichstatus);
      }
      othercpupos[1] = l % ncpux1;
      othercpupos[2] = (int) ((l % (numprocs)) / ncpux1);
      // now fill jonio with proper sequence (i.e. tiled mapping)
      for (j = 0; j < N2; j++)
	for (i = 0; i < N1; i++)
	  for (col = 0; col < numcolumns; col++) {
	    mapvaluejonio =
	      + ncpux1 * N1 * numcolumns * (j + othercpupos[2] * N2)
	      + numcolumns * (i + othercpupos[1] * N1)
	      + col;
	    mapvaluetempbuf =
	      + j * N1 * numcolumns
	      +  i * numcolumns + col;
	    
	    if (datatype == sizeof(unsigned char))
	      jonio1[mapvaluejonio] = writebuf1[mapvaluetempbuf];
	    if (datatype == sizeof(float))
	      jonio4[mapvaluejonio] = writebuf4[mapvaluetempbuf];
	    if (datatype == sizeof(double))
	      jonio8[mapvaluejonio] = writebuf8[mapvaluetempbuf];
	  }
    }
  }
  if(myid!=0) MPI_Wait(&srequest, &mpichstatus);
  free(writebuf);		// writebuf used by each CPU

  if (myid == 0) {
    // now write out collected data using CPU=0
    if(bintxt==BINARYOUTPUT){
      fwrite(jonio, datatype,
	     totalsize[1] * totalsize[2] * numcolumns, *fp);
    }
    else if(bintxt==TEXTOUTPUT){ // properly ordered, so just dump it
      for(i=0;i<totalsize[1]*totalsize[2]*numcolumns;i++){
	if (datatype == sizeof(unsigned char))
	  fprintf(*fp,"%u",jonio1[i]);
	if (datatype == sizeof(float))
	  fprintf(*fp,"%15.7g",jonio4[i]);
	if (datatype == sizeof(double))
	  fprintf(*fp,"%21.15g",jonio8[i]);
	if((i+1)%numcolumns) fprintf(*fp," ");
	else fprintf(*fp,"\n");
      }
    }
    free(jonio);		// used by CPU=0
    fclose(*fp);
    *fp = NULL;
  }
#endif

 logsfprintf("mpiios end combine\n");

}

void mpiios_seperate(int bintxt, int stage, MPI_Datatype mpidt, int numcolumns,
		    int datatype, FILE ** fp, void *jonio,
		    void *writebuf)
{
  // baesd on sizeof()
  // 1: unsigned char
  // 4: float
  // 8: double
  int i, j, k, l, col, mapvaluejonio, mapvaluetempbuf;
#if(USEMPI)
  MPI_Request rrequest;
  MPI_Request srequest;
#endif
  int othercpupos[COMPDIM + 1];
  unsigned char *jonio1;
  float *jonio4;
  double *jonio8;
  unsigned char *writebuf1;
  float *writebuf4;
  double *writebuf8;

 logsfprintf("mpiios begin seperate\n");


  if (datatype == sizeof(unsigned char))
    jonio1 = (unsigned char *) jonio;
  else if (datatype == sizeof(float))
    jonio4 = (float *) jonio;
  else if (datatype == sizeof(double))
    jonio8 = (double *) jonio;
  if (datatype == sizeof(unsigned char))
    writebuf1 = (unsigned char *) writebuf;
  else if (datatype == sizeof(float))
    writebuf4 = (float *) writebuf;
  else if (datatype == sizeof(double))
    writebuf8 = (double *) writebuf;

#if(USEMPI)

  if (stage == 1) {
    if (myid == 0) {
      if(bintxt==BINARYOUTPUT){
	// first let cpu=0 read data
	fread(jonio, datatype,
	      totalsize[1] * totalsize[2] * numcolumns,
	      *fp);
      }
      else if(bintxt==TEXTOUTPUT){ // properly ordered, so just dump it
	for(i=0;i<totalsize[1]*totalsize[2]*numcolumns;i++){
	  if (datatype == sizeof(unsigned char))
	    fscanf(*fp,"%u",&jonio1[i]);
	  if (datatype == sizeof(float))
	    fscanf(*fp,"%f",&jonio4[i]);
	  if (datatype == sizeof(double))
	    fscanf(*fp,"%lf",&jonio8[i]);
	}
      }
    }
    // writebuf is CPU=0's tempbuf for each CPU, including CPU=0, which is done last
    if (myid == 0) {
      for (l = numprocs-1 ; l >=0; l--) {
	othercpupos[1] = l % ncpux1;
	othercpupos[2] = (int) ((l % (numprocs)) / ncpux1);
	// now unfill jonio with proper sequence (i.e. tiled mapping)
	for (j = 0; j < N2; j++)
	  for (i = 0; i < N1; i++)
	    for (col = 0; col < numcolumns; col++) {
	      mapvaluejonio =
		    + ncpux1 * N1 * numcolumns * (j +
						  othercpupos[2] * N2)
		    + numcolumns * (i + othercpupos[1] * N1)
		    + col;
	      mapvaluetempbuf =
		+ j * N1 * numcolumns
		+ i * numcolumns + col;
		if (datatype == sizeof(unsigned char))
		  writebuf1[mapvaluetempbuf] = jonio1[mapvaluejonio];
		if (datatype == sizeof(float))
		  writebuf4[mapvaluetempbuf] = jonio4[mapvaluejonio];
		if (datatype == sizeof(double))
		  writebuf8[mapvaluetempbuf] = jonio8[mapvaluejonio];
	      }
	if(l!=0){
	  MPI_Isend(writebuf, N1 * N2 * numcolumns, mpidt, l, l,
		    MPI_COMM_WORLD, &srequest);
	  MPI_Wait(&srequest, &mpichstatus);
	}
      }
      free(jonio); // done with jonio after loop
    }
    else{
      // chosen CPU to receive data from CPU=0
      MPI_Irecv(writebuf, N1 * N2 * numcolumns, mpidt, 0, myid,
		MPI_COMM_WORLD, &rrequest);
      MPI_Wait(&rrequest, &mpichstatus);	// writebuf used until      
    }
  } else if (stage == 2) {
    free(writebuf);
    if (myid == 0) {
      fclose(*fp);
      *fp = NULL;
    }
  }
#endif

 logsfprintf("mpiios end seperate\n");



}


void mpiiotu_combine(MPI_Datatype mpidt, int numcolumns, int datatype,
		   FILE ** fp, void *writebuf)
{
  // based on sizeof()
  // 1: unsigned char
  // 4: float
  // 8: double
  int i, j, k, l, col;
#if(USEMPI)
  MPI_Request rrequest;
  MPI_Request srequest;
#endif
  int othercpupos[COMPDIM + 1];
  unsigned char *writebuf1;
  float *writebuf4;
  double *writebuf8;
  

  if (datatype == sizeof(unsigned char))
    writebuf1 = (unsigned char *) writebuf;
  else if (datatype == sizeof(float))
    writebuf4 = (float *) writebuf;
  else if (datatype == sizeof(double))
    writebuf8 = (double *) writebuf;

#if(USEMPI)
  if(myid!=0) MPI_Isend(writebuf, N1 * N2 * numcolumns, mpidt, 0, myid,
			MPI_COMM_WORLD, &srequest);
  if (myid == 0) {    
    // done in forward order, no need to use tempbuf since CPU=0's writebuf is first out
    for (l = 0; l <numprocs; l++) {
      if(l!=0){
	MPI_Irecv(writebuf, N1 * N2 * numcolumns, mpidt, l, l,
		  MPI_COMM_WORLD, &rrequest);
	MPI_Wait(&rrequest, &mpichstatus);
      }
      // now write writebuf
      DUMPLOOP(0,N1-1,0,N2-1){
	for(col=0;col<numcolumns;col++){
	  if(datatype==sizeof(unsigned char)){
	    fprintf(*fp,"%c ",writebuf1[col+numcolumns*(i+N1*j)]);
	  }
	  else if(datatype==sizeof(float)){
	    fprintf(*fp,"%15.7g ",writebuf4[col+numcolumns*(i+N1*j)]);
	  }
	  else if(datatype==sizeof(double)){
	    fprintf(*fp,"%21.15g ",writebuf8[col+numcolumns*(i+N1*j)]);
	  }
	}
	fprintf(*fp,"\n");
      }
    }    
  }
  if(myid!=0) MPI_Wait(&srequest, &mpichstatus);
  free(writebuf);		// writebuf used by each CPU

  if (myid == 0) {
    fclose(*fp);
    *fp = NULL;
  }
#endif

}

// fill writebuf with each cpu's data set,using CPU=0 to process the file
void mpiiotu_seperate(int stage, MPI_Datatype mpidt, int numcolumns,
		    int datatype, FILE ** fp,void *writebuf)
{
  // baesd on sizeof()
  // 1: unsigned char
  // 4: float
  // 8: double
  int i, j, k, l, col;
#if(USEMPI)
  MPI_Request rrequest;
  MPI_Request srequest;
#endif
  int othercpupos[COMPDIM + 1];
  unsigned char *writebuf1;
  float *writebuf4;
  double *writebuf8;
  void *tempbuf;
  unsigned char *tempbuf1;
  float *tempbuf4;
  double *tempbuf8;
  unsigned char *sendbuf1;
  float *sendbuf4;
  double *sendbuf8;
  void *sendbuf;

  if (datatype == sizeof(unsigned char))
    writebuf1 = (unsigned char *) writebuf;
  else if (datatype == sizeof(float))
    writebuf4 = (float *) writebuf;
  else if (datatype == sizeof(double))
    writebuf8 = (double *) writebuf;

  if(myid==0){
    if((tempbuf=malloc(datatype*N1*N2*numcolumns))==NULL){
      dualfprintf(fail_file,"Can't open tempbuf in gammieio_sep\n");
      myexit(1);
    }

    if (datatype == sizeof(unsigned char))
      tempbuf1 = (unsigned char *) tempbuf;
    else if (datatype == sizeof(float))
      tempbuf4 = (float *) tempbuf;
    else if (datatype == sizeof(double))
      tempbuf8 = (double *) tempbuf;
  }

#if(USEMPI)

  if (stage == 1) {
    if (myid == 0) {
      for (l = 0; l < numprocs; l++) {
	if(l==0){
	  sendbuf=writebuf;
	  sendbuf1=writebuf1;
	  sendbuf4=writebuf4;
	  sendbuf8=writebuf8;
	}
	else{
	  sendbuf=tempbuf;
	  sendbuf1=tempbuf1;
	  sendbuf4=tempbuf4;
	  sendbuf8=tempbuf8;
	}

	DUMPLOOP(0,N1-1,0,N2-1) for (col = 0; col < numcolumns; col++) {
	  if(datatype==sizeof(unsigned char)) fscanf(*fp,"%u",&sendbuf1[col+numcolumns*(i+N2*j)]);
	  else if(datatype==sizeof(float)) fscanf(*fp,"%f",&sendbuf4[col+numcolumns*(i+N2*j)]);
	  else if(datatype==sizeof(double)) fscanf(*fp,"%lf",&sendbuf8[col+numcolumns*(i+N2*j)]);
	}
	if(l!=0){
	  MPI_Isend(sendbuf, N1 * N2 * numcolumns, mpidt, l, l,
		    MPI_COMM_WORLD, &srequest);
	  // have to wait before filling sendbuf buffer again for next CPU
	  MPI_Wait(&srequest, &mpichstatus);
	}
      }
      free(tempbuf);
    }
    else{
      MPI_Irecv(writebuf, N1 * N2 * numcolumns, mpidt, 0, myid, MPI_COMM_WORLD, &rrequest);
      MPI_Wait(&rrequest, &mpichstatus);	// writebuf used until
    }
  } else if (stage == 2) {
    free(writebuf);
    if (myid == 0) {
      fclose(*fp);
      *fp = NULL;
    }
  }
#endif



}




#endif



// a simple max, assumes local cpu max already found
// sends results back to all cpus
void mpimax(SFTYPE*maxptr)
{
  SFTYPE send;
  
#if(USEMPI)
  send = *maxptr;
  MPI_Allreduce(&send, maxptr, 1, MPI_SFTYPE, MPI_MAX,MPI_COMM_WORLD);
#endif
}

// a simple max, assumes local cpu max already found
// sends results back to all cpus
void mpimin(SFTYPE*minptr)
{
  SFTYPE send;
  
#if(USEMPI)
  send = *minptr;
  MPI_Allreduce(&send, minptr, 1, MPI_SFTYPE, MPI_MIN,MPI_COMM_WORLD);
#endif
}

// a simple max, assumes local cpu max already found
// sends results back to all cpus
void mpifmin(FTYPE*minptr)
{
  FTYPE send;
  
#if(USEMPI)
  send = *minptr;
  MPI_Allreduce(&send, minptr, 1, MPI_FTYPE, MPI_MIN,MPI_COMM_WORLD);
#endif
}


void prminmaxsum(FTYPE p[][N2+4][NPR], int start,int nmemb, FTYPE *maxptr, FTYPE*minptr,FTYPE*sumptr)
{
  int i,j,k;
  FTYPE maxsend,minsend,sumsend;
  int domin,domax,dosum;

  if(maxptr==NULL) domax=0; else domax=1;
  if(minptr==NULL) domin=0; else domin=1;
  if(sumptr==NULL) dosum=0; else dosum=1;
  
  for(k=start;k<start+nmemb;k++){
    
    if(domin) minptr[k]=1E30;
    if(domax) maxptr[k]=-1E30;
    if(dosum) sumptr[k]=0;
  }
  ZLOOP {
    for(k=start;k<start+nmemb;k++){
      if(domax) if (p[i][j][k] > maxptr[k]) maxptr[k] = p[i][j][k];
      if(domin) if (p[i][j][k] < minptr[k]) minptr[k] = p[i][j][k];
      if(dosum) sumptr[k]+=p[i][j][k];
    }
  }
#if(USEMPI)
  for(k=start;k<start+nmemb;k++){    
    if(domax){
      maxsend = maxptr[k];
      MPI_Allreduce(&maxsend, &maxptr[k], 1, MPI_FTYPE, MPI_MAX, MPI_COMM_WORLD);
    }
    if(domin){
      minsend = minptr[k];
      MPI_Allreduce(&minsend, &minptr[k], 1, MPI_FTYPE, MPI_MIN, MPI_COMM_WORLD);
    }
    if(dosum){
      sumsend = sumptr[k];
      MPI_Allreduce(&sumsend, &sumptr[k], 1, MPI_FTYPE, MPI_SUM, MPI_COMM_WORLD);
    }
  }
#endif
}



void myfprintf(FILE* fileptr, char *format, ...)
{
  va_list arglist;
  if  (myid==0) {
    va_start (arglist, format);

    if(fileptr==NULL){
      fprintf(stderr,"tried to print to null file pointer: %s\n",format);
      fflush(stderr);
    }
    else{
      vfprintf (fileptr, format, arglist);
      fflush(fileptr);
    }
    va_end (arglist);
  }
}

// prints to stderr(only cpu=0) AND file pointer of choice (all cpus)
void dualfprintf(FILE* fileptr, char *format, ...)
{
  va_list arglist;

  va_start (arglist, format);

  if(fileptr==NULL){
    fprintf(stderr,"tried to print to null file pointer: %s\n",format);
    fflush(stderr);
  }
  else{
    vfprintf (fileptr, format, arglist);
    fflush(fileptr);
  }
  if(myid==0){
    vfprintf (stderr, format, arglist);
    fflush(stderr);
  }
  va_end (arglist);
}

// prints to both logfull_file(cpu=0) and log_file(all cpus)
void logsfprintf(char *format, ...)
{
  va_list arglist;
  
  va_start (arglist, format);

  if  ((myid==0)&&(logfull_file)){
    vfprintf (logfull_file, format, arglist);    
    fflush(logfull_file);
  }
  if(log_file){
    vfprintf (log_file, format, arglist);    
    fflush(log_file);
  }
  va_end (arglist);
}

// prints to logfull_file, log_file, and stderr (but only using cpu=0)
void trifprintf(char *format, ...)
{
  va_list arglist;
  
  va_start (arglist, format);
  if  ((myid==0)&&(logfull_file)){
    vfprintf (logfull_file, format, arglist);    
    fflush(logfull_file);
  }
  if(log_file){
    vfprintf (log_file, format, arglist);    
    fflush(log_file);
  }
  if(myid==0){
    vfprintf (stderr, format, arglist);    
    fflush(stderr);
  }
  va_end (arglist);
}


void myfopen(char*fname, char*fmt,char*message,FILE**fileptrptr)
{
  if(myid==0){
    *fileptrptr = fopen(fname, fmt);
    if (*fileptrptr == NULL) {
      dualfprintf(fail_file, message);
      myexit(1);
    }
  }
}

void myfclose(FILE ** fileptrptr,char*message)
{
  int reterror;

  if(myid==0){
    if(*fileptrptr!=NULL){
      reterror = fclose(*fileptrptr);
      if (reterror == EOF) {
	dualfprintf(fail_file, message);
	myexit(1);
      }
    }
    else{
      dualfprintf(fail_file,"file already closed: %s\n",message);
      myexit(1);
    }
  }
}


// only accepts FTYPE
void mydfwrite(FTYPE *ptr, int start, size_t nmemb, int i, int j, FILE*stream,FTYPE*writebuf)
{
  int k;
  if(mpicombine==0){
    if(binaryoutput) fwrite(ptr+start, sizeof(FTYPE), nmemb, stream);
    else{
      for(k=start;k<start+nmemb;k++){
	if(sizeof(FTYPE)==sizeof(float)) fprintf(stream,"%15.7g ",ptr[k]);
	else if(sizeof(FTYPE)==sizeof(double)) fprintf(stream,"%21.15g ",ptr[k]);
      }
    }
  }
  else{
    for(k=start;k<start+nmemb;k++) writebuf[BUFFERMAP] = ptr[k];
  }
}

// only accepts FTYPE
void mydfread(FILE*stream, int start, size_t nmemb,int i, int j, FTYPE*writebuf, FTYPE *ptr)
{
  int k;
  if(mpicombine==0){
    if(binaryoutput) fread(ptr+start, sizeof(FTYPE), nmemb, stream);
    else{
      for(k=start;k<start+nmemb;k++){
	if(sizeof(FTYPE)==sizeof(float)) fscanf(stream,"%f",&ptr[k]);
	else if(sizeof(FTYPE)==sizeof(double)) fscanf(stream,"%lf",&ptr[k]);
      }
    }
  }
  else for(k=start;k<start+nmemb;k++) ptr[k]=writebuf[BUFFERMAP];
}
