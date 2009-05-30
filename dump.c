#include "decs.h"

int dump(int dump_cnt)
{
  int i = 0, j = 0, k = 0, l = 0, col = 0;
  FTYPE r, th, vmin1, vmax1, vmin2, vmax2;
  struct of_geom geom;
  struct of_state q;
  FTYPE X[NDIM];
  FTYPE divb;
  FTYPE tens_em[NDIM][NDIM], tens_matt[NDIM][NDIM], b[NDIM],
      ucon[NDIM];
  FTYPE U[NPR];
  FILE *fp;
  char dfnam[MAXFILENAME];
#if(USEMPI)
  void *jonio;
#endif
  void *writebuf;
  FTYPE *realbuf;
  char truemyidtxt[MAXFILENAME],temptext[MAXFILENAME];
  FILE *headerptr;
  FTYPE ftemp;


  trifprintf("begin dumping dump# %ld ... ",dump_cnt);

  //////////////////////////////////
  //
  //  Define file output and open it
  //
  ///////////////////////////////////

  numcolumns = 2*3 + NPR + 1 + NDIM * 2*2 + 2*2 + 1;	// 36 currently

  strcpy(truemyidtxt,"");
  if(!binaryoutput)  strcat(truemyidtxt, "");	// always text
  else strcat(truemyidtxt, ".bin");	// always binary
  if(GAMMIE) sprintf(dfnam, "dumps/dump%03ld%s", dump_cnt, truemyidtxt);
  else sprintf(dfnam, "dumps/dump%04ld%s", dump_cnt, truemyidtxt);

  if( mpicombine == 0 ) {
    if(!binaryoutput) strcpy(temptext,"wt");
    else strcpy(temptext,"w");
    
    if ((fp = fopen(dfnam, temptext)) == NULL) {
      dualfprintf(fail_file, "error opening dump file\n");
      myexit(2);
    }
    writebuf=NULL;
  }
  else{
#if(USEMPI)
    mpiio_init(binaryoutput,sortedoutput, &fp, WRITEFILE, dfnam, numcolumns, sizeof(FTYPE), &jonio, &writebuf);
    realbuf=(FTYPE*)writebuf;
#endif
  }

  //////////////////////////////////
  //
  //  Open and write and close HEADER file
  //
  ///////////////////////////////////

  if(!binaryoutput) headerptr=fp;

  // header
  if(binaryoutput){
    strcat(dfnam, ".head");
    myfopen(dfnam,"wt","Can't open dump header for writing\n",&headerptr);
  }
  myfprintf(headerptr, "%10.5g %d %d %10.5g %10.5g %10.5g %10.5g %ld %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g\n", t,
	    totalsize[1], totalsize[2], startx[1], startx[2], dx[1],
	    dx[2],realnstep,gam,a,R0,Rin,Rout,hslope);
  if(binaryoutput){
    myfclose(&headerptr,"Couldn't close dump header\n");
  }

  //////////////////
  //
  // DUMP LOOP
  //
  //////////////////

  DUMPLOOP(0, N1 - 1, 0, N2 - 1) {
    // buffer init starts the parallel index
    BUFFERINIT;
    coord(i, j, CENT, X);
    bl_coord(X, &r, &th);
    // if failed, then data output for below invalid, but columns still must exist    
    get_geometry(i, j, CENT, &geom);
    if (!failed) {
      if (get_state(p[i][j], &geom, &q) >= 1)
	FAILSTATEMENT("dump.c:dump()", "get_state() dir=0", 1);
      if (vchar(p[i][j], &q, 1, &geom, &vmax1, &vmin1) >= 1)
	FAILSTATEMENT("dump.c:dump()", "vchar() dir=1or2", 1);
      if (vchar(p[i][j], &q, 2, &geom, &vmax2, &vmin2) >= 1)
	FAILSTATEMENT("dump.c:dump()", "vchar() dir=1or2", 2);
    }
    SETFDIVB(divb, p, i, j);
    // if you change # of outputted vars, remember to change numcolumns above
    ftemp=(FTYPE)(i+startpos[1]);
    mydfwrite(&ftemp,0,1,i,j,fp,realbuf);
    ftemp=(FTYPE)(j+startpos[2]);
    mydfwrite(&ftemp,0,1,i,j,fp,realbuf);

    mydfwrite(X,1,2,i,j,fp,realbuf);
    mydfwrite(&r,0,1,i,j,fp,realbuf);
    mydfwrite(&th,0,1,i,j,fp,realbuf);
    mydfwrite(p[i][j],0,NPR,i,j,fp,realbuf);
    mydfwrite(&divb,0,1,i,j,fp,realbuf);

    for (k = 0; k < NDIM; k++)
      mydfwrite(&(q.ucon[k]),0,1,i,j,fp,realbuf);
    for (k = 0; k < NDIM; k++)
      mydfwrite(&(q.ucov[k]),0,1,i,j,fp,realbuf);
    for (k = 0; k < NDIM; k++)
      mydfwrite(&(q.bcon[k]),0,1,i,j,fp,realbuf);
    for (k = 0; k < NDIM; k++)
      mydfwrite(&(q.bcov[k]),0,1,i,j,fp,realbuf);
    
    mydfwrite(&vmin1,0,1,i,j,fp,realbuf);
    mydfwrite(&vmax1,0,1,i,j,fp,realbuf);
    mydfwrite(&vmin2,0,1,i,j,fp,realbuf);
    mydfwrite(&vmax2,0,1,i,j,fp,realbuf);
    
    mydfwrite(&geom.g,0,1,i,j,fp,realbuf);
    if((mpicombine==0)&&(!binaryoutput)) fprintf(fp,"\n");
  }// end DUMPLOOP


  //////////////////
  //
  // Close dump file
  //
  //////////////////


  if (mpicombine == 0) {
    if (fp != NULL)
      fclose(fp);
  } else {
#if(USEMPI)
    mpiio_combine(binaryoutput, sortedoutput, numcolumns, sizeof(FTYPE), &fp, jonio, writebuf);
#endif
  }

  trifprintf("end dumping dump# %ld\n",dump_cnt);

  return (0);
}


void gdump(void)
{
  int i = 0, j = 0, k = 0, l = 0, m = 0, n = 0, col = 0;
  FTYPE r, th;
  FTYPE X[NDIM];
  FILE *fp;
  char dfnam[MAXFILENAME];
#if(USEMPI)
  void *jonio;
#endif
  void *writebuf;
  FTYPE *realbuf;
  char truemyidtxt[MAXFILENAME],temptext[MAXFILENAME];
  FILE *headerptr;
  FTYPE ftemp;
  FTYPE *ptrftemp;


  trifprintf("begin dumping gdump ... ");

  //////////////////////////////////
  //
  //  Define file output and open it
  //
  ///////////////////////////////////

  numcolumns = 2*3+NDIM*NDIM*NDIM+NPG*NDIM*NDIM*2+NPG; // 202 currently

  strcpy(truemyidtxt,"");
  if(!binaryoutput)  strcat(truemyidtxt, "");	// always text
  else strcat(truemyidtxt, ".bin");	// always binary
  sprintf(dfnam, "dumps/gdump%s", truemyidtxt);

  if( mpicombine == 0 ) {
    if(!binaryoutput) strcpy(temptext,"wt");
    else strcpy(temptext,"w");
    
    if ((fp = fopen(dfnam, temptext)) == NULL) {
      dualfprintf(fail_file, "error opening gdump file\n");
      myexit(2);
    }
    writebuf=NULL;
  }
  else{
#if(USEMPI)
    mpiio_init(binaryoutput,sortedoutput, &fp, WRITEFILE, dfnam, numcolumns, sizeof(FTYPE), &jonio, &writebuf);
    realbuf=(FTYPE*)writebuf;
#endif
  }

  //////////////////////////////////
  //
  //  Open and write and close HEADER file
  //
  ///////////////////////////////////


  // header
  if(!binaryoutput) headerptr=fp;

  if(binaryoutput){
    strcat(dfnam, ".head");
    myfopen(dfnam,"wt","Can't open gdump header for writing\n",&headerptr);
  }
  myfprintf(headerptr, "%10.5g %d %d %10.5g %10.5g %10.5g %10.5g %ld %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g\n", t,
	    totalsize[1], totalsize[2], startx[1], startx[2], dx[1],
	    dx[2],realnstep,gam,a,R0,Rin,Rout,hslope);
  if(binaryoutput){
    myfclose(&headerptr,"Couldn't close dump header\n");
  }

  //////////////////
  //
  // DUMP LOOP
  //
  //////////////////

  DUMPLOOP(0, N1 - 1, 0, N2 - 1) {
    // buffer init starts the parallel index
    BUFFERINIT;
    coord(i, j, CENT, X);
    bl_coord(X, &r, &th);
    // if failed, then data output for below invalid, but columns still must exist    
    // if you change # of outputted vars, remember to change numcolumns above
    ftemp=(FTYPE)(i+startpos[1]);
    mydfwrite(&ftemp,0,1,i,j,fp,realbuf);
    ftemp=(FTYPE)(j+startpos[2]);
    mydfwrite(&ftemp,0,1,i,j,fp,realbuf);
    // 2
    mydfwrite(X,1,2,i,j,fp,realbuf);
    mydfwrite(&r,0,1,i,j,fp,realbuf);
    mydfwrite(&th,0,1,i,j,fp,realbuf);
    // 2+4
    ptrftemp=(FTYPE*)(&conn[i][j][0][0][0]);
    mydfwrite(ptrftemp,0,NDIM*NDIM*NDIM,i,j,fp,realbuf);
    ptrftemp=(FTYPE*)(&gcon[i][j][0][0][0]);
    mydfwrite(ptrftemp,0,NPG*NDIM*NDIM,i,j,fp,realbuf);
    ptrftemp=(FTYPE*)(&gcov[i][j][0][0][0]);
    mydfwrite(ptrftemp,0,NPG*NDIM*NDIM,i,j,fp,realbuf);
    ptrftemp=(FTYPE*)(&gdet[i][j][0]);
    mydfwrite(ptrftemp,0,NPG,i,j,fp,realbuf);

    if((mpicombine==0)&&(!binaryoutput)) fprintf(fp,"\n");
  }// end DUMPLOOP


  //////////////////
  //
  // Close dump file
  //
  //////////////////


  if (mpicombine == 0) {
    if (fp != NULL)
      fclose(fp);
  } else {
#if(USEMPI)
    mpiio_combine(binaryoutput, sortedoutput, numcolumns, sizeof(FTYPE), &fp, jonio, writebuf);
#endif
  }

  trifprintf("end dumping gdump\n");


}
