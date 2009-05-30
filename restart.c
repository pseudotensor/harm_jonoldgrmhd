
/* restart functions; restart_init and restart_dump */

#include "decs.h"

void restart_write(int which)
{
  FILE *fp;
  int idum, i, j, k, l;
  char dfnam[200];
  // used for binaryoutput==1
#if(USEMPI)
  void *jonio;
  FTYPE totalnorm, recvnorm, sendnorm;
  int ndims, array_of_gsizes[4], array_of_distribs[4];
  int order, len;
  int array_of_dargs[4], array_of_psizes[4];
  int bufcount, array_size;
#endif
  void *writebuf;
  FTYPE *realbuf;
  // end used by binaryoutput==1
  char truemyidtxt[100],temptext[100];
  FILE *headerptr;

  //////////////////////////////////
  //
  // Restart write file setup
  //
  //////////////////////////////////

  trifprintf("Writing restart file\n");

  numcolumns = NPR;
  
  if(!binaryoutput)  strcpy(truemyidtxt, "");	// always text
  else strcpy(truemyidtxt, ".bin");	// always binary

  if(GAMMIE) sprintf(dfnam, "dumps/rdump");
  else sprintf(dfnam, "dumps/rdump-%01d%s", which, truemyidtxt);

  //////////////////////////////////
  //
  // Restart write file OPEN
  //
  //////////////////////////////////


  if (mpicombine == 0) {
    if(!binaryoutput) strcpy(temptext,"wt");
    else strcpy(temptext,"wb");

    if( (fp = fopen(dfnam, temptext))==NULL){
      dualfprintf(fail_file, "can't open restart file\n");
      myexit(2);
    }
    realbuf=NULL;
  } else {
#if(USEMPI)
    mpiio_init(binaryoutput,sortedoutput, &fp, WRITEFILE, dfnam, numcolumns, sizeof(FTYPE), &jonio, &writebuf);
    realbuf=(FTYPE*)writebuf;
#endif

  }

  //////////////////////////////////
  //
  // Restart write HEADER open/write/close
  //
  //////////////////////////////////
  if(!binaryoutput) headerptr=fp;

  write_restart_header(dfnam,headerptr);

  //////////////////////////////////
  //
  // Restart write DUMP 
  //
  //////////////////////////////////


  DUMPLOOP(0, N1 - 1, 0, N2 - 1) {
    BUFFERINIT;
    mydfwrite(p[i][j],0,NPR,i,j,fp,realbuf);
    if((mpicombine==0)&&(!binaryoutput)) fprintf(fp,"\n");
  }

  //////////////////////////////////
  //
  // Restart write file CLOSE
  //
  //////////////////////////////////


  if (mpicombine == 0) {
    myfclose(&fp,"cannot close fp in restart_write\n");
  } else {
#if(USEMPI)
    mpiio_combine(binaryoutput, sortedoutput, numcolumns, sizeof(FTYPE), &fp, jonio, writebuf);
#endif

  }
  

  trifprintf("end restart write\n");
}

int restart_init(int which)
{
  char ans[100];

  trifprintf("begin restart init\n");

  ranc(7);
  restart_read(which);

  trifprintf("proc: %d t=%12.15g\n", myid, t);

  /* set metric functions */
  set_grid();

  trifprintf("proc: %d grid restart completed\n", myid);
  
  fixup(p);

  trifprintf( "proc: %d fixup restart completed\n", myid);
  

  /* bound */
  if (bound_prim(p) >= 1) {
    fprintf(fail_file, "restart_init:bound_prim: failure\n");
    fflush(fail_file);
    return (1);
  }

  trifprintf( "proc: %d bound restart completed\n", myid);
  

  trifprintf( "proc: %d restart completed\n", myid);
  
  trifprintf("end restart init\n");

  /* done! */
  return (0);

}

void restart_read(int which)
{
  int idum,idum1,idum2, i, j, k, l;
  char dfnam[100];
  FILE *fp;
#if(USEMPI)
  void *jonio;
  FTYPE totalnorm, recvnorm, sendnorm;
  int ndims, array_of_gsizes[4], array_of_distribs[4];
  int order, len;
  int array_of_dargs[4], array_of_psizes[4];
  int bufcount, array_size;
#endif
  void *writebuf;
  FTYPE *realbuf;
  char truemyidtxt[100],temptext[100];
  int range;
  FILE *headerptr;
  FTYPE dtemp;
  int readingdump;


  //////////////////////////////////
  //
  // Restart read file setup
  //
  //////////////////////////////////

  // whether reading directly from dump file or not (normally don't)
  readingdump = 0;
  numcolumns = NPR;

  if(!binaryoutput)  strcpy(truemyidtxt, "");	// always text
  else strcpy(truemyidtxt, ".bin");	// always binary


  // can read in dump instead of rdump
  if(!GAMMIE){
    if (!readingdump) {
      sprintf(dfnam, "dumps/rdump-%01d%s", which, truemyidtxt);
    } else {
      sprintf(dfnam, "dumps/indump%s", truemyidtxt);
    }
  }
  else sprintf(dfnam, "dumps/rdump");

  //////////////////////////////////
  //
  // Restart read file OPEN
  //
  //////////////////////////////////
  trifprintf("opening restart file\n");

  if (mpicombine == 0) {
    if(!binaryoutput) strcpy(temptext,"rt");
    else strcpy(temptext,"r");

    fp = fopen(dfnam, temptext);
    if (fp == NULL) {
      dualfprintf(fail_file, "no restart file\n");
      myexit(1);
    }
    writebuf=NULL;
  } else {
#if(USEMPI)
    mpiio_init(binaryoutput,sortedoutput, &fp, READFILE, dfnam, numcolumns, sizeof(FTYPE), &jonio, &writebuf);
    realbuf=(FTYPE*)writebuf;
#endif

  }

  trifprintf("restart file exists...loading: %s ...\n", dfnam);


  //////////////////////////////////
  //
  // Restart read HEADER open/write/close
  //
  //////////////////////////////////

  if(!binaryoutput) headerptr=fp;

  read_restart_header(dfnam,headerptr);

  trifprintf("reading dump contents of restart file\n");

  //////////////////////////////////
  //
  // Restart read DUMP 
  //
  //////////////////////////////////

  if (mpicombine == 1) {
#if(USEMPI)
    mpiio_seperate(binaryoutput,sortedoutput, STAGE1, numcolumns, sizeof(FTYPE), &fp, jonio, writebuf);
#endif
  }
  
  trifprintf("assigning restart data\n");
  
  DUMPLOOP(0, N1 - 1, 0, N2 - 1) {
    BUFFERINIT;
    
    // skip up to primitives
    if (readingdump) for(k=0;k<4;k++) mydfread(fp,0,1,i,j,realbuf,&dtemp);
    
    mydfread(fp,0,NPR,i,j,writebuf,p[i][j]);
    
    // skip after primitives
    if (readingdump) for(k=0;k<24;k++) mydfread(fp,0,1,i,j,realbuf,&dtemp);
  }

  //////////////////////////////////
  //
  // Restart write file CLOSE
  //
  //////////////////////////////////

  if (mpicombine == 0)
    myfclose(&fp,"Cannot close fp in restart_read\n");
  else {
#if(USEMPI)
      mpiio_seperate(binaryoutput,sortedoutput, STAGE2, numcolumns, sizeof(FTYPE), &fp, jonio, writebuf);
#endif
  }

  // test read by looking at written file
  //  restart_write(3);
  trifprintf("Restart read file completed\n");
  

}


// headerptr created and only used here OR passed a given pointer
void read_restart_header(char *dfnam,FILE*headerptr)
{
  char headername[MAXFILENAME];
  char temptext[MAXFILENAME];
  int idum1,idum2;

  trifprintf("begin reading header of restart file\n");

  if(myid==0){
    
    if(binaryoutput&&(dfnam!=NULL)){
      strcpy(headername,dfnam);
      strcat(headername, ".head");
      strcpy(temptext,"rb");
      myfopen(headername,temptext,"Can't open restart header for writing\n",&headerptr);

      /* read in global variables, in binary */
      fread(&idum1, sizeof(int), 1, headerptr);
      fread(&idum2, sizeof(int), 1, headerptr);
      // all cpus read the rest of header the same
      fread(&t, sizeof(SFTYPE), 1, headerptr);
      fread(&tf, sizeof(SFTYPE), 1, headerptr);
      fread(&nstep, sizeof(long), 1, headerptr);
      fread(&a, sizeof(FTYPE), 1, headerptr);
      fread(&gam, sizeof(FTYPE), 1, headerptr);
      fread(&cour, sizeof(FTYPE), 1, headerptr);
      
      fread(&DTd, sizeof(SFTYPE), 1, headerptr);
      fread(&DTener, sizeof(SFTYPE), 1, headerptr);
      fread(&DTi, sizeof(SFTYPE), 1, headerptr);
      // fread(&DTr, sizeof(SFTYPE), 1, headerptr) ;
      fread(&DTr, sizeof(long), 1, headerptr);
      fread(&dump_cnt, sizeof(long), 1, headerptr);
      fread(&image_cnt, sizeof(long), 1, headerptr);
      fread(&rdump_cnt, sizeof(long), 1, headerptr);
      
      fread(&dt, sizeof(SFTYPE), 1, headerptr);
      fread(&lim, sizeof(int), 1, headerptr);
      fread(&failed, sizeof(int), 1, headerptr);
      
      fread(&R0, sizeof(FTYPE), 1, headerptr);
      fread(&Rin, sizeof(FTYPE), 1, headerptr);
      fread(&Rout, sizeof(FTYPE), 1, headerptr);
      fread(&hslope, sizeof(FTYPE), 1, headerptr);
      fread(&defcoord, sizeof(FTYPE), 1, headerptr);

      fread(&BCtype[X1UP],sizeof(int), 1, headerptr);
      fread(&BCtype[X1DN],sizeof(int), 1, headerptr);
      fread(&BCtype[X2UP],sizeof(int), 1, headerptr);
      fread(&BCtype[X2DN],sizeof(int), 1, headerptr);

      myfclose(&headerptr,"Cannot close headerprt in restart_read\n");

    }
    // assumes headerptr is normal file pointer already defined
    else{
      fscanf(headerptr,RESTARTHEADER,
	     &idum1,&idum2,
	     &t,&tf,&nstep,&a,&gam,&cour,
	     &DTd,&DTener,&DTi,&DTr,&dump_cnt,&image_cnt,&rdump_cnt,
	     &dt,&lim,&failed,
	     &R0,&Rin,&Rout,&hslope,&defcoord,
	     &BCtype[X1UP],&BCtype[X1DN],&BCtype[X2UP],&BCtype[X2DN]
	     );
    }
    if (idum1 != totalsize[1]) {
      dualfprintf(fail_file, "error reading restart file; N1 differs\n");
      dualfprintf(fail_file, "got totalsize[1]=%d needed totalsize[1]=%d\n",idum1,totalsize[1]);
      myexit(3);
    }
    if (idum2 != totalsize[2]) {
      dualfprintf(fail_file, "error reading restart file; N2 differs\n");
      dualfprintf(fail_file, "got totalsize[2]=%d needed totalsize[2]=%d\n",idum2,totalsize[2]);
      myexit(4);
    }
  }

  // now that CPU=0 has the restart header, pass to other CPUs
#if(USEMPI)
  MPI_Bcast(&t, 1, MPI_SFTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&tf, 1, MPI_SFTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nstep, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&a, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&gam, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&cour, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  
  MPI_Bcast(&DTd, 1, MPI_SFTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&DTener, 1, MPI_SFTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&DTi, 1, MPI_SFTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&DTr, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&dump_cnt, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&image_cnt, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&rdump_cnt, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  
  MPI_Bcast(&dt, 1, MPI_SFTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&lim, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&failed, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  MPI_Bcast(&R0, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Rin, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Rout, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&hslope, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&defcoord, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);

  MPI_Bcast(&BCtype[X1UP], 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&BCtype[X1DN], 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&BCtype[X2UP], 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&BCtype[X2DN], 1, MPI_INT, 0, MPI_COMM_WORLD);

#endif

  // write header to screen
  fprintf(log_file,"header contents below\n"); fflush(log_file);
  write_restart_header(NULL,log_file);

  // other assignments
  realnstep=nstep;

  trifprintf("end reading header of restart file\n");
    
}


void write_restart_header(char *dfnam,FILE*headerptr)
{
  char headername[MAXFILENAME];
  char temptext[MAXFILENAME];

  trifprintf("begin writing header of restart file\n");

  if(myid==0){
    if(binaryoutput&&(dfnam!=NULL)){

      strcpy(temptext,"wb");
      strcpy(headername,dfnam);
      strcat(headername, ".head");
      myfopen(headername,temptext,"Can't open restart header for writing\n",&headerptr);

      /* write out key global variables, in binary */
      fwrite(&totalsize[1], sizeof(int), 1, headerptr);
      fwrite(&totalsize[2], sizeof(int), 1, headerptr);
      
      fwrite(&t, sizeof(SFTYPE), 1, headerptr);
      fwrite(&tf, sizeof(SFTYPE), 1, headerptr);
      fwrite(&nstep, sizeof(long), 1, headerptr);
      // for now nstep and realnstep same
      fwrite(&a, sizeof(FTYPE), 1, headerptr);
      fwrite(&gam, sizeof(FTYPE), 1, headerptr);
      fwrite(&cour, sizeof(FTYPE), 1, headerptr);
      
      fwrite(&DTd, sizeof(SFTYPE), 1, headerptr);
      fwrite(&DTener, sizeof(SFTYPE), 1, headerptr);
      fwrite(&DTi, sizeof(SFTYPE), 1, headerptr);
      // fwrite(&DTr, sizeof(SFTYPE),1, headerptr) ;
      fwrite(&DTr, sizeof(long), 1, headerptr);
      fwrite(&dump_cnt, sizeof(long), 1, headerptr);
      fwrite(&image_cnt, sizeof(long), 1, headerptr);
      fwrite(&rdump_cnt, sizeof(long), 1, headerptr);
      
      fwrite(&dt, sizeof(SFTYPE), 1, headerptr);
      fwrite(&lim, sizeof(int), 1, headerptr);
      fwrite(&failed, sizeof(int), 1, headerptr);
      
      fwrite(&R0, sizeof(FTYPE), 1, headerptr);
      fwrite(&Rin, sizeof(FTYPE), 1, headerptr);
      fwrite(&Rout, sizeof(FTYPE), 1, headerptr);
      fwrite(&hslope, sizeof(FTYPE), 1, headerptr);
      fwrite(&defcoord, sizeof(FTYPE), 1, headerptr);

      fwrite(&BCtype[X1UP],sizeof(int), 1, headerptr);
      fwrite(&BCtype[X1DN],sizeof(int), 1, headerptr);
      fwrite(&BCtype[X2UP],sizeof(int), 1, headerptr);
      fwrite(&BCtype[X2DN],sizeof(int), 1, headerptr);

      myfclose(&headerptr,"Cannot close headerprt in restart_write\n");      
    }
  }
  // text output, either to just headerptr as file on myid==0, or all cpus if log_file is file pointer
  if((!binaryoutput)&&(myid==0)&&(headerptr!=log_file) ||
     (dfnam==NULL)&&(headerptr==log_file)) {
    fprintf(headerptr,WRITERESTARTHEADER,
	    totalsize[1],totalsize[2],
	    t,tf,nstep,a,gam,cour,
	    DTd,DTener,DTi,DTr,dump_cnt,image_cnt,rdump_cnt,
	    dt,lim,failed,
	    R0,Rin,Rout,hslope,defcoord,
	    BCtype[X1UP],BCtype[X1DN],BCtype[X2UP],BCtype[X2DN]
	    );
    fflush(headerptr);
  }


  trifprintf("end writing header of restart file\n");

}
