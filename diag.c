
#include "decs.h"

/* diagnostics subroutine */

int diag(int call_code)
{
  static FILE *ener_file;
  FILE*dump_file,*image_file;
  int gotit;
  long fpos0 = 0;
  static SFTYPE tlast;
  char dfnam[MAXFILENAME], ifnam[MAXFILENAME];
  int i, j, k, l, dir;
  FILE *imagecnt_file, *dumpcnt_file;
  SFTYPE U_init[NPR],U_tot[NPR],U_final[NPR], divb, divbmax = 0, divbavg =0, fladd_tot[NPR];
  SFTYPE norm[NPR];
  SFTYPE pcum_tot[COMPDIM*2][NPR];
  int imax, jmax;
  static int firsttime = 1;
  static int dumpc, imagec, restartc, enerc;
  static SFTYPE tdump, timage, tener;
  static int nrestart;
  // static SFTYPE trestart ;

  if(firsttime) makedirs();

#if(1)
  if ((call_code == INIT_OUT) || (firsttime == 1)) {
    if (appendold == 0) {
      myfopen(ENERFNAME,"w","error opening energy output file\n",&ener_file);
    } else {			// if appendold==1
      
      trifprintf("Start setup of %s file append\n", dfnam);
      
      myfopen(ENERFNAME,"a+","error opening energy output file for append\n",&ener_file);
      appendener(ener_file,pcum_tot,fladd_tot);
    }
  }
#endif

#if(1)
  if ((call_code == INIT_OUT) || (firsttime == 1)) {
    
    tlast = t - SMALL;
    tdump = timage = tener = t;
    
    // counts are good since tstart is t to start with
    dumpc = imagec = restartc = enerc = 0;

    if (RESTARTMODE == 0) {
      // trestart = t;
      nrestart = nstep;
      tdump = timage = tener = t;
      dump_cnt = 0;
      image_cnt = 0;
      rdump_cnt = 0;
      appendold = 0;
    } else {
      setrestart(&dump_cnt,&image_cnt,&rdump_cnt,&appendold);

      // assuming started at t=0 and nstep=0 for original run
      tdump  = DTd*fmod(t,DTd);
      timage = DTi*fmod(t,DTi);
      tener  = DTener*fmod(t,DTener);
      nrestart = DTr*(nstep%DTr);
    }
  }
#endif


#if(1)
  // output grid (probaly want both fullgrid (to make sure ok) and compute grid to compare with data dumps
  if(firsttime&&(RESTARTMODE==0)) gdump();
#endif


  // t!=tlast avoids duplicate entries
  if (t != tlast) {


#if(1)
    if (call_code == INIT_OUT ||
	call_code == LOG_OUT || 
	call_code == FINAL_OUT) {


      // compute total conserved quantity
      if(integrate(U_tot,U_tot,CONSTYPE)>=1) return(1);
      DIRLOOP if(integrate(pdot[dir],pdot[dir],SURFACETYPE)>=1) return(1);
      DIRLOOP if(integrate(pcum[dir],pcum_tot[dir],SURFACETYPE)>=1) return(1);
      if(integrate(fladd,fladd_tot,CUMULATIVETYPE)>=1) return(1);

      divbmaxavg(p,&divbmax,&divbavg);
    }

    if ((call_code == INIT_OUT) || (firsttime == 1)) {
      PLOOP U_init[k] = U_tot[k];
    }
    if (call_code == FINAL_OUT) {
      PLOOP U_final[k] = U_tot[k];
      
      if(GAMMIE){
	dualfprintf(logfull_file,"\n\nEnergy: ini,fin,del: %g %g %g\n",
		 U_init[UU], U_final[UU], (U_final[UU] - U_init[UU]) / U_init[UU]);

	dualfprintf(logfull_file,"\n\nMass: ini,fin,del: %g %g %g\n",
		 U_init[RHO], U_final[RHO], (U_final[RHO] - U_init[RHO]) / U_init[RHO]);
      }
      else{
	PLOOP dualfprintf(logfull_file,"\n\nU[%d]: ini,fin,del: %g %g %g\n",
		       k, U_init[k], U_final[k], U_init[k]!=0.0 ? (U_final[k] - U_init[k]) / U_init[k] : 0.0 );
	
      }
    }

    if (((call_code == INIT_OUT) && (firsttime == 1)) || (t >= tener)) {
      if(GAMMIE){
	// 6+3=9 terms
	myfprintf(ener_file, "%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g ", t, U_tot[RHO], U_tot[U3], U_tot[UU], p[N1 / 2][N2 / 2][UU] * pow(p[N1 / 2][N2 / 2][RHO], -gam), p[N1 / 2][N2 / 2][UU]);
	myfprintf(ener_file, "%21.15g %21.15g %21.15g ", pdot[RHO], pdot[UU], pdot[U3]);
      }
      else{
	// 3+NPR+COMPDIM*2*NPR+NPR+3 = 54 terms
	myfprintf(ener_file,"%21.15g %21.15g %21.15g ",t ,p[N1 / 2][N2 / 2][UU] * pow(p[N1 / 2][N2 / 2][RHO], -gam), p[N1 / 2][N2 / 2][UU]);
	PLOOP myfprintf(ener_file, "%21.15g ", U_tot[k]);
	DIRLOOP PLOOP myfprintf(ener_file, "%21.15g ", pdot[dir][k]);
	// COMPDIM*2*NPR more terms
	DIRLOOP PLOOP myfprintf(ener_file, "%21.15g ", pcum_tot[dir][k]);
	PLOOP myfprintf(ener_file, "%21.15g ", fladd_tot[k]);
	myfprintf(ener_file, "%ld %21.15g %21.15g ", realnstep, divbmax,divbavg);
	// NPR*3 more terms
	PLOOP myfprintf(ener_file, "%21.15g %21.15g %21.15g ", p[N1/4][3*N2/8][k],p[N1/4][4*N2/8][k],p[N1/4][5*N2/8][k]);
      }
      myfprintf(ener_file,"\n");
      fflush(ener_file);

      frdotout();

      while (t >= tener) {
	enerc++;
	tener = tstart + enerc * DTener;
      }
    }
    if(call_code == FINAL_OUT) myfclose(&ener_file,"Couldn't close ener_file\n");  

#endif

    ///////////////////////
    //
    // RESTART DUMP
    //
    ///////////////////////


#if(1)
    if ((failuremode == 0)
	&& (nstep >= nrestart || call_code == 2 || call_code == 0
	    || firsttime == 1)) {
      trifprintf("dumping: restart: %d\n", whichrestart);
      restart_write(whichrestart);	// 0 1 0 1 0 1 ...
      restartsteps[whichrestart] = realnstep;
      whichrestart = !whichrestart;
      /* 
         while(t>=trestart){ restartc++; trestart=tstart+restartc*DTr ;
         } */
      while (nstep >= nrestart) {
	restartc++;
	nrestart = restartc * DTr;
      }

    }
#endif


    ///////////////////////
    //
    // DUMP
    //
    ///////////////////////

#if(1)
    /* dump at regular intervals */
    if (t >= tdump || call_code == 2 || call_code == 0
	|| firsttime == 1 || (RESTARTMODE&&dofaildump&&(nstep>=steptofaildump))  ) {
      trifprintf("dumping: dump_cnt=%ld\n", dump_cnt);
      
      /* make regular dump file */
      if (dump(dump_cnt) >= 1)	return (1);

      // iterate counter
      dump_cnt++;
      while (t >= tdump) {
	dumpc++;
	tdump = tstart + dumpc * DTd;
      }
      if(!GAMMIE){
	// output number of dumps
	myfopen("dumps/0_numdumps.dat","w","error opening dump count file\n",&dumpcnt_file);      
	myfprintf(dumpcnt_file, "# Number of dumps\n%ld\n", dump_cnt);
	myfclose(&dumpcnt_file,"Couldn't close dumpcnt_file");
      }
    }
#endif
    ///////////////////////
    //
    // AREA MAP
    //
    ///////////////////////

    if(dofailmap){
      if (nstep>=steptofailmap) {
	if(area_map(call_code, TIMESERIESAREAMAP, 20, ifail, jfail, p)>=1) return(1);
      }
    }

    ///////////////////////
    //
    // IMAGE
    //
    ///////////////////////

#if(1)
    /* image dump at regular intervals */
    if (t >= timage || call_code == 2 || call_code == 0
	|| firsttime == 1) {
      trifprintf("image dump t=%15.10g\n", t);

      /* make regular image file */
      image_dump(image_cnt);

      // iterate counter
      image_cnt++;
      while (t >= timage) {
	imagec++;
	timage = tstart + imagec * DTi;
      }
      if(!GAMMIE){
	// output number of images
	myfopen("images/0_numimages.dat","w","error opening image count file\n",&imagecnt_file);      
	myfprintf(imagecnt_file, "# Number of images\n%ld\n", image_cnt);
	myfclose(&imagecnt_file,"Couldn't close imagecnt_file");
      }
    }
#endif


  }
  firsttime = 0;
  tlast = t;
  return (0);
}

/** some diagnostic routines **/

/* map out region around failure point */
int area_map(int call_code, int type, int size, int i, int j, FTYPE prim[][N2 + 4][NPR])
{
  int k;
  int l,m,ll,mm;
  FTYPE r, th, vmin1, vmax1, vmin2, vmax2;
  struct of_geom geom;
  struct of_state q;
  FTYPE X[NDIM];
  FTYPE divb;
  FTYPE tens_em[NDIM][NDIM], tens_matt[NDIM][NDIM], b[NDIM],
      ucon[NDIM];
  FTYPE U[NPR];

  static FILE* fileptr;
  static int firsttime=1;
  static int domap=0;
  static int doclose=0;

  trifprintf("\nStart area_map function ... ");

  if(firsttime){
    if((type==TIMESERIESAREAMAP)&&(dofailmap)){
      if((fileptr=fopen("areamap","wt"))==NULL){
	dualfprintf(fail_file,"Cannot open ./areamap on proc=%d\n",myid);
	domap=0;
      }
      else domap=1;
    }
  }

  if((type==TIMESERIESAREAMAP)&&domap&&(call_code==2)){
    doclose=1;
  }
  else doclose=0;

  if(type==FINALTDUMPAREAMAP){
    dualfprintf(fail_file, "area map\n");
    dualfprintf(fail_file, "failure at: i=%d j=%d\n",i+startpos[1],j+startpos[2]);
    coord(i,j,CENT,X);
    dualfprintf(fail_file, "failure at: \n",i+startpos[1],j+startpos[2]);
    

    PLOOP {
      
      dualfprintf(fail_file, "variable %d \n", k);
      
      dualfprintf(fail_file, "i = \t ");
      for(l=i-size/2;l<=i+size/2;l++){
	if((l<-N1BND)||(l>N1-1+N1BND)) continue;
	ll=l+startpos[1];
	dualfprintf(fail_file, "%12d", ll);
      }
      dualfprintf(fail_file, "\n");
      for(m=j-size/2;m<=j+size/2;m++){
	if((m<-N2BND)||(m>N2-1+N2BND)) continue;
	mm=m+startpos[2];
	dualfprintf(fail_file, "j = %d \t ",mm);
	for(l=i-size/2;l<=i+size/2;l++){
	  if((l<-N1BND)||(l>N1-1+N1BND)) continue;
	  ll=l+startpos[1];
	  dualfprintf(fail_file, "%12.5g ",prim[l][m][k]);
	}
	dualfprintf(fail_file, "\n");
      }
    }
  }
  else if((type==TIMESERIESAREAMAP)&&(domap)){
    if(firsttime){
      fprintf(fileptr,"%21.15g %d %d %21.15g %21.15g %21.15g %21.15g %d %d %d %d %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g\n",
	      t,totalsize[1],totalsize[2],startx[1],startx[2],dx[1],dx[2],size,size,startpos[1]+i,startpos[2]+j,gam,a,R0,Rin,Rout,hslope);
      fflush(fileptr);
    }
    for(m=j-size/2;m<=j+size/2;m++){
      if((m<-N2BND)||(m>N2-1+N2BND)) continue;
      mm=m+startpos[2];
      for(l=i-size/2;l<=i+size/2;l++){
	if((l<-N1BND)||(l>N1-1+N1BND)) continue;
	ll=l+startpos[1];
	
	
	coord(l, m, CENT, X);
	bl_coord(X, &r, &th); 
	get_geometry(l, m, CENT, &geom);
	if (!failed) {
	  if (get_state(p[l][m], &geom, &q) >= 1)
	    FAILSTATEMENT("dump.c:dump()", "get_state() dir=0", 1);
	  if (vchar(p[l][m], &q, 1, &geom, &vmax1, &vmin1) >= 1)
	    FAILSTATEMENT("dump.c:dump()", "vchar() dir=1or2", 1);
	  if (vchar(p[l][m], &q, 2, &geom, &vmax2, &vmin2) >= 1)
	    FAILSTATEMENT("dump.c:dump()", "vchar() dir=1or2", 2);
	}
	
	if((l>=-1)&&(l<=N1+1)&&(m>=-1)&&(m<=N2+1)){ SETFDIVB(divb, p, l, m);}
	else divb=0.0;

	// same order as dump.c for first columns (easy sm read)
	fprintf(fileptr,
		"%d %d "
		"%g %g "
		"%g %g "
		"%g %g %g %g %g %g %g %g "
		"%g "
		"%g %g %g %g "
		"%g %g %g %g "
		"%g %g %g %g "
		"%g %g %g %g "
		"%g %g %g %g "
		"%g "
		"%g %ld\n",
		ll,mm,
		X[1],X[2],
		r,th,
		p[l][m][0],
		p[l][m][1],
		p[l][m][2],
		p[l][m][3],
		p[l][m][4],
		p[l][m][5],
		p[l][m][6],
		p[l][m][7],
		divb,
		q.ucon[0],q.ucon[1],q.ucon[2],q.ucon[3],
		q.ucov[0],q.ucov[1],q.ucov[2],q.ucov[3],
		q.bcon[0],q.bcon[1],q.bcon[2],q.bcon[3],
		q.bcov[0],q.bcov[1],q.bcov[2],q.bcov[3],
		vmin1,vmax1,vmin2,vmax2,
		geom.g,
		t,realnstep);
      }
    }
    fflush(fileptr);
  }

  if(doclose) if(fileptr==NULL) fclose(fileptr);


  /* print out other diagnostics here */

  firsttime=0;
  trifprintf("end area_map function.\n");  
  return(0);
}

/* evaluate fluxed based diagnostics; put results in global variables */


void diag_flux(FTYPE F1[][N2 + 4][NPR], FTYPE F2[][N2 + 4][NPR],SFTYPE Dt)
{
  int i, j, k, dir;
  SFTYPE surface;
  int *iter1,*iter2;
  FTYPE (*flux)[N2 + 4][NPR];
  int start1,start2,stop1,stop2;

  // true outer boundary surface fluxes (per direction, per conserved variable)
  DIRLOOP {
    if (doflux[dir] >= 0) {
      // otherwise don't add to it
      if((dir==X1UP)||(dir==X1DN)){
	surface = dx[2]*dx[3];
	start1=stop1=doflux[dir];
	start2=0;
	stop2=N2-1;
	flux=F1;
      }
      else if((dir==X2UP)||(dir==X2DN)){
	surface = dx[1]*dx[3];
	start2=stop2=doflux[dir];
	start1=0;
	stop1=N1-1;
	flux=F2;
      }
      PLOOP{
	pdot[dir][k]=0;
	GENLOOP(i,j,start1,stop1,start2,stop2){
	  pdot[dir][k]  += flux[i][j][k] * surface;
	}
	pcum[dir][k]+=pdot[dir][k]*Dt;
      }
    }
  }
#if(1)
  // radial flux vs. radius
  flux=F1;
  surface = dx[2]*dx[3];
  PLOOP{
    for(i=0;i<N1;i++){
      frdot[i][k]=0;
      for(j=0;j<N2;j++){
	frdot[i][k]+=flux[i][j][k]*surface;
      }
    }
  }
#endif

}

void frdotout(void)
{
  int i,j,k,l;
  SFTYPE ftemp;
  MPI_Request rrequest;
  MPI_Request srequest;
  SFTYPE frdottemp[N1][NPR];
  SFTYPE *frtot;
  int ospos1;
  FILE*frout;
  static int firsttime=1;

  if(numprocs==1){
    frtot=(SFTYPE (*))(&frdot[0][0]);
  }
  else{
#if(USEMPI)
    if(myid==0){
      frtot=(SFTYPE*) malloc(sizeof(SFTYPE)*totalsize[1]*NPR);
      if(frtot==NULL){
	dualfprintf(fail_file,"Cannot get frtot memory\n");
	myexit(1);
      }
      else{
	for(i=0;i<totalsize[1];i++) PLOOP{
	  frtot[i*NPR+k]=0;
	}
      }
      for(l=0;l<numprocs;l++){
	ospos1=(l%ncpux1)*N1;
	if(l==0){
	  for(i=0;i<N1;i++) PLOOP{
	    frdottemp[i][k]=frdot[i][k];
	  }
	}
	else{
	  MPI_Irecv(&frdottemp,N1*NPR,MPI_SFTYPE,l,l,MPI_COMM_WORLD,&rrequest);
	  MPI_Wait(&rrequest,&mpichstatus);
	}
	for(i=0;i<N1;i++) PLOOP{
	  frtot[(ospos1+i)*NPR+k]+=frdottemp[i][k];
	}
      }
    }
    if(myid!=0){
      MPI_Isend(&frdot,N1*NPR,MPI_SFTYPE,0,myid,MPI_COMM_WORLD,&srequest);
      MPI_Wait(&srequest,&mpichstatus);
    }
#endif
  }
  // now we have frtot with full fluxes vs. radius (totalsize[1]), so output

  if(myid==0){
    frout=fopen("frdot.out","at");
    if(frout==NULL){
      dualfprintf(fail_file,"Cannot open frdot.out\n");
      myexit(1);
    }
    if(firsttime){
      fprintf(frout,"%21.15g %ld %d %d %21.15g %21.15g %21.15g %21.15g\n",
	      t,realnstep,totalsize[1],totalsize[2],startx[1],startx[2],dx[1],dx[2]);
      fflush(frout);
    }

    for(i=0;i<totalsize[1];i++){
      fprintf(frout,"%21.15g %d ",t,i);
      PLOOP{
	fprintf(frout,"%21.15g ",frtot[i*NPR+k]);
      }
      fprintf(frout,"\n");
    }
    
    fclose(frout);
    free(frtot);
  }
  firsttime=0;
}


void makedirs(void)
{

  if ((USEMPI == 0) || (USEMPI && (!USEGM))) {
      if ((mpicombine && (myid == 0)) || (mpicombine == 0)) {
	system("mkdir dumps");
	system("mkdir images");
      }
#if(USEMPI)
      MPI_Barrier(MPI_COMM_WORLD);	// all cpus wait for directory
      // to be created
#endif
  }
}


void appendener(FILE* ener_file,SFTYPE pcum_tot[][NPR],SFTYPE*fladd_tot)
{
  int gotit;
  SFTYPE tcheck;
  int l,k,dir;
  long fpos0;
  FILE *ener_file_temp;
  char dfnam[MAXFILENAME],dfnamback[MAXFILENAME], dfnamtemp[MAXFILENAME];

  // only CPU=0 does anything here
  if(myid==0){
    
    rewind(ener_file);	// go to start
    
    gotit = 0;
    while ((!feof(ener_file)) && (gotit == 0)) {
      
      fscanf(ener_file, "%lf", &tcheck);
      
      if (fabs(tcheck - t) < 0.5 * DTener) {
	gotit = 1;
	for (l = 1; l <= NUMENERVAR; l++) {
	  if ((l > 3+NPR+COMPDIM*2*NPR) && (l < 3+NPR+2*COMPDIM*2*NPR+NPR)) {
	    DIRLOOP PLOOP {
	      fscanf(ener_file, "%lf", &pcum_tot[dir][k]);
	      l++;
	    }
	    PLOOP {
	      fscanf(ener_file, "%lf", &fladd_tot[k]);
	      l++;
	    }
	  }
	}
      } else {
	// skip this bad line 
	while ((fgetc(ener_file) != '\n') && (!feof(ener_file)));
      }
      // continue after successful get since successful get is good 
      // data and should keep since corresponds to dump one is
      // keeping
      fpos0 = ftell(ener_file);	// position to continue
      // writting at if successful get
    }
    if (gotit == 0) {
      dualfprintf(fail_file,
	      "Never found right time in loss file when appending: looking for t=%21.15g lastt=%21.15g\n",
	      t, tcheck);
      myexit(1);
    } else {
      dualfprintf(logfull_file,
	      "found goodtime t=%21.15g (wanted %21.15g) to restart ener file\n",
	      tcheck, t);
      sprintf(dfnamtemp, "%s0_ener%s.temp", DATADIR, ".out");
      sprintf(dfnam, "%sener%s", DATADIR, ".out");
      sprintf(dfnamback, "%s0_ener%s.back", DATADIR, ".out");
      
      // now that done, fix up file
      if ((ener_file_temp = fopen(dfnamtemp, "wt")) == NULL) {
	dualfprintf(fail_file,
		"Cannot open temp ener file for appending: %s\n",
		dfnamtemp);
	myexit(1);
      } else {
	rewind(ener_file);
	while (ftell(ener_file) < fpos0 + 1) {	// +1 is for
	  // '\n' at end
	  // of line
	  fputc(fgetc(ener_file), ener_file_temp);
	}
	fclose(ener_file_temp);
	fclose(ener_file);
	rename(dfnam, dfnamback);	// move old to backup location
	rename(dfnamtemp, dfnam);	// move new to old name(normal
	// name)
	// reopen loss_file (now normal name)

	if ((ener_file = fopen(dfnam, "at")) == NULL) {
	  dualfprintf(fail_file,
		  "2: error opening ener output file %s\n", dfnam);
	  myexit(1);
	}
	trifprintf("End setup of ener file append\n");
      }			// end else if can open temp file
    }			// end else if gotit==1
  }
}

void divbmaxavg(FTYPE p[][N2+4][NPR],SFTYPE*ptrdivbmax,SFTYPE*ptrdivbavg)
{
  int i,j,k;
  int imax=0,jmax=0;
  SFTYPE divb;
  SFTYPE divbmax=0,divbavg=0;
  SFTYPE divbmaxsend,divbavgsend;

  LOOPDIVB {
    // doesn't need geom, just use global gdet
    SETFDIVB(divb, p, i, j);
    if (divb > divbmax) {
      imax = i;
      jmax = j;
      divbmax = divb;
    }
    divbavg += divb;
  }

  // PER CPU
  fprintf(log_file,"  proc: %04d : divbmax: %d %d %21.15g divbavg: %21.15g\n",
	  myid, imax, jmax, divbmax, divbavg / ((FTYPE) N1*N2));

#if(USEMPI)			// give CPU=0 total
  divbmaxsend = divbmax;
  divbavgsend = divbavg;
  MPI_Reduce(&divbmaxsend, &divbmax, 1, MPI_SFTYPE, MPI_MAX, 0,
	     MPI_COMM_WORLD);
  MPI_Reduce(&divbavgsend, &divbavg, 1, MPI_SFTYPE, MPI_SUM, 0,
	     MPI_COMM_WORLD);
#endif
  divbavg /= (FTYPE) (totalzones);

  // Total over all CPUs
  myfprintf(logfull_file,"  divbmax: %21.15g divbavg: %21.15g\n",divbmax, divbavg);
  *ptrdivbmax=divbmax;
  *ptrdivbavg=divbavg;
}


/* gettotal accepts an arbitrary pointer set each of different sizes
 * i.e. one could do:
 * numptrs=2+NPR;
 * totalptrs[0]=pdot;  totalsizes[0]=NPR; totaloptrs[0]=pdot;
 * totalptrs[1]=fladd; totalsizes[1]=NPR; totaloptrs[1]=fladd_tot;
 * PLOOP{ totalptrs[2+k]=&U_tot[k]; totalsizes[k]=1;}
 * gettotal(numptrs,totalptrs,totaloptrs,totalsizes);
 *
 */

void gettotal(int numvars, SFTYPE* vars[],int*sizes,SFTYPE*vars_tot[])
{
  int j,k;
  SFTYPE send;
  
  // for 1 CPU
  if (numprocs == 1) {
    // must use _tot since can't overwrite normal global integrator 
    // variable
    for(k=0;k<numvars;k++){
      for(j=0;j<sizes[k];j++){
	vars_tot[k][j] = vars[k][j];
      }
    }
  } else {
    // give CPU=0 the totals
#if(USEMPI)
    for(k=0;k<numvars;k++){
      for(j=0;j<sizes[k];j++){
	send=vars[k][j];
	// send and receive can't be same address, hence "send" variable
	MPI_Reduce(&send, &vars_tot[k][j], 1, MPI_SFTYPE, MPI_SUM, 0,
		   MPI_COMM_WORLD);
      }
    }
#endif
  }
}

// each CPU does constotal
int constotal(SFTYPE *vars)
{
  int i,j,k;
  FTYPE U[NPR];
  struct of_geom geom;
  struct of_state q;


  PLOOP vars[k]= 0.0;

  ZLOOP {
    get_geometry(i,j,CENT,&geom) ;
    if(!failed){
      if(get_state(p[i][j],&geom,&q)>=1) return(1);
      if(primtoU(p[i][j],&q,&geom,U)>=1) return(1);
    }
    for(k=0;k<NPR;k++){
      vars[k] += U[k]*dVF;
    }
  }
  return(0);
}

int integrate(SFTYPE * var,SFTYPE *var_tot,int type)
{
  SFTYPE *totalptrs[100],*totaloptrs[100];
  int totalsizes[100],numptrs;

  switch(type){
  case CONSTYPE:
    if(constotal(var)>=1) return(1);
    totalsizes[0]=NPR;
    gettotal(1,&var,totalsizes,&var_tot);
    break;
  case SURFACETYPE:
  case CUMULATIVETYPE:
    totalsizes[0]=NPR;
    gettotal(1,&var,totalsizes,&var_tot);
    break;
  default:
    dualfprintf(fail_file,"No defined type=%d in integrate\n",type);
    myexit(1);
  }
  return(0);
}



void setrestart(long *dump_cnt, long*image_cnt,long *rdump_cnt,int*appendold)
{
  // restart at given count (count will be written at "t=0" ==
  // tstart=t
  // set by restart
  //  *dump_cnt = 19;
  //*image_cnt = 458;
  //*rdump_cnt = 0;
  // or set to 0 and copy file and do manually (if problem)
  *appendold = 0;
}
