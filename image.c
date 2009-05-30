/* 
   produces an "r8" file. */

#include "decs.h"

#define GAMMIESTARTI (0)
#define GAMMIEENDI (0)

#define JONSTARTI (0)
#define JONENDI (NPR-1)

#define NORMAL (0)
#define ZOOM (1)
#define NUMLIMITS (2)

#define LOG (0)
#define LINEAR (1)
#define NUMSCALE (2)

#define MINVECTOR (1E-10)


void image_dump(int image_cnt)
{
  int k;
  int startk,endk,starts,ends,startl,endl;
  int scale,limits;
  ////////////////////////////
  //
  // Image Loop
  //
  ////////////////////////////

  // STARTI=0 is normal
  // ENDI=NPR is normal, but can be NPR+2 to get linear RHO/UU
  if(GAMMIE){
    startk=GAMMIESTARTI;
    endk=GAMMIEENDI;
    starts=0;
    ends=1;
    startl=0;
    endl=0;
  }
  else{
    startk=JONSTARTI;
    endk=JONENDI;
    starts=0;
    ends=NUMSCALE-1;
    startl=0;
    endl=NUMLIMITS-1;
  }

  for(scale=starts;scale<=ends;scale++){
    for(limits=startl;limits<=endl;limits++){
      for (k = startk; k <= endk; k++) {
	image(image_cnt, k, scale, limits);
      }
    }
  }

}


void image(int image_cnt, int which, int scale, int limits)
{
  int i = 0, j = 0, l = 0, col = 0;
  FILE *fp;
  char ifnam[MAXFILENAME];
  // which : which primitive variable
  SFTYPE pr,iq, liq, aa, lmax, lmin;
  FTYPE X[NDIM],r,th;
  FTYPE min,max,sum;
  FTYPE minptr[NPR], maxptr[NPR], sumptr[NPR];
  int jonhead;
#if(USEMPI)
  void *jonio;
  int ndims, array_of_gsizes[4], array_of_distribs[4];
  int order, len;
  int array_of_dargs[4], array_of_psizes[4];
  int bufcount, array_size;
#endif
  void *writebuf;
  unsigned char *realbuf;
  char truemyidtxt[MAXFILENAME];
  FTYPE (*pimage)[N2+4][NPR];

  ////////////////////////////
  //
  // Image output setup/definition
  //
  ////////////////////////////

  pimage=ph;
  if(limits==ZOOM){
    ZLOOP{
      if(which<=1){
	coord(i,j,CENT,X);
	bl_coord(X,&r,&th);
	if(which==0) pimage[i][j][which]=p[i][j][which]/(RHOMIN*pow(r,-1.5));
	if(which==1) pimage[i][j][which]=p[i][j][which]/(UUMIN*pow(r,-2.5));
      }
      else{
	if(scale==LINEAR) pimage[i][j][which]=p[i][j][which];
	else if(scale==LOG) pimage[i][j][which]=fabs(p[i][j][which])+MINVECTOR;
      }
    }
  }
  else{
    ZLOOP{
      if(which<=1) pimage[i][j][which]=p[i][j][which];
      else{
	if(scale==LINEAR) pimage[i][j][which]=p[i][j][which];
	else if(scale==LOG) pimage[i][j][which]=fabs(p[i][j][which])+MINVECTOR;
      }
    }
  }
  numcolumns = 1;
  if(GAMMIE) jonhead=0;
  else jonhead=1;
  strcpy(truemyidtxt, "");

  ////////////////////////////
  //
  // Image FILE open/initialize
  //
  ////////////////////////////


  // always binary
  if(GAMMIE&&(GAMMIESTARTI==GAMMIEENDI)&&(GAMMIESTARTI==0)){
    sprintf(ifnam, "images/im%04ld", image_cnt);
  }
  else sprintf(ifnam, "images/im%1dp%1ds%1dl%04ld%s.r8", which, scale, limits, image_cnt,truemyidtxt);

  if (mpicombine == 0) {
    fp = fopen(ifnam, "w");

    if (fp == NULL) {
      dualfprintf(fail_file, "error opening image file\n");
      myexit(2);
    }
    writebuf=realbuf=NULL;
  } else {
#if(USEMPI)
    mpiio_init(BINARYOUTPUT,SORTED,&fp, WRITEFILE, ifnam, numcolumns, sizeof(unsigned char), &jonio, &writebuf);
    realbuf=(unsigned char*)writebuf;
#endif
  }

  ////////////////////////////
  //
  // HEADER file open/initialize
  //
  ////////////////////////////

  // write header
  if (jonhead == 1) {
    if(mpicombine==0) myfprintf(fp, "RAW\n# t=%15.10g p=%2d\n%i %i\n255\n", t, which, N1,N2);
    else  myfprintf(fp, "RAW\n# t=%15.10g p=%2d\n%i %i\n255\n", t, which, totalsize[1],totalsize[2]);
  }

  ////////////////////
  //
  // Image paramters setup
  // 
  /////////////////////

  /* density mapping is logarithmic, in 255 steps between e^lmax and
     e^lmin */

#define ZOOMFACTOR (10000)

  prminmaxsum(pimage,which,1,maxptr,minptr,sumptr);
  if(limits==NORMAL){
    max=maxptr[which];
    min=minptr[which];
    sum=sumptr[which];
  }
  else{
    if(which<=1){
      max=maxptr[which]/ZOOMFACTOR;
      min=minptr[which];
    }
    else{
      if(scale==LINEAR){
	max=maxptr[which]/ZOOMFACTOR;
	min=minptr[which]/ZOOMFACTOR;
      }
      else{
	max=maxptr[which]/ZOOMFACTOR;
	min=minptr[which];
      }
    }
  }
  sum=sumptr[which];
  logsfprintf("which: %d scale: %d limits: %d : min,max,avg: %g %g %g\n",which,scale,limits,min,max,sum/totalzones);

  if(scale==LOG){
    lmax = log(max);
    lmin = log(min);
  } else if(scale==LINEAR) {
    lmax = max;
    lmin = min;
  }

  if (lmax != lmin)
    aa = 256. / (lmax - lmin);
  else
    aa = 0;


  ////////////////////
  //
  // Image DUMP Loop
  // 
  /////////////////////


  IMAGELOOP(0, N1 - 1, 0, N2 - 1) {
    pr=pimage[i][j][which];
    if (scale==LOG) iq = log(pr);
    else if(scale==LINEAR) iq = pr;

    /* liq = aa*log(iq) + b ; */
    liq = aa * (iq - lmin);
    if (liq > 255.)
      liq = 255.;
    if (liq < 0.)
      liq = 0.;
    if (mpicombine == 0)
      myfprintf(fp, "%c", (char) ((int) liq));
    else
      realbuf[IMAGEMAP] = (int) liq;
  }


  ////////////////////
  //
  // Image dump file CLOSE
  // 
  /////////////////////

  if (mpicombine == 0) {
    myfclose(&fp,"Cannot close fp in image()\n");
  } else {
#if(USEMPI)
    mpiio_combine(BINARYOUTPUT, SORTED, numcolumns, sizeof(unsigned char), &fp,
		  jonio, writebuf);
#endif
  }
}
