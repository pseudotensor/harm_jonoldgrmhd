
#include "decs.h"

int fail(int fail_type)
{

  dualfprintf(fail_file, "\n\nfail: i=%d j=%d ft=%d\n",startpos[1]+ icurr, startpos[2]+jcurr, fail_type);

  // failed use for diag below to avoid bad calculations
  failed = 1 ;

  switch (fail_type) {
  case 1:
    dualfprintf(fail_file, "fail_type=%s\n", FAILSTR01);
    break;
  case 2:
    dualfprintf(fail_file, "fail_type=%s\n", FAILSTR02);
    break;
  case 3:
    dualfprintf(fail_file, "fail_type=%s\n", FAILSTR03);
    break;
  case 4:
    dualfprintf(fail_file, "fail_type=%s\n", FAILSTR04);
    break;
  case 5:
    dualfprintf(fail_file, "fail_type=%s\n", FAILSTR05);
    break;
  case 6:
    dualfprintf(fail_file, "fail_type=%s\n", FAILSTR06);
    break;
  default:
    dualfprintf(fail_file, "fail_type=unknown\n");
    break;
  }


  dualfprintf(fail_file, "failed\n");
  if(area_map(1,FINALTDUMPAREAMAP, 3, icurr, jcurr, p)>=1){}
  // do nothing since never will fail since if failed=1, doesn't do fail-related function calls

  // do not respond here since not in synch with other CPUs.  Use postdt().

  /* for diagnostic purposes */
  return (1);
}


void postdt(void)
{
  static SFTYPE aftertime,beforetime;
  static int firsttime=1;
  static long beforenstep,afternstep;
  static FTYPE cour0;
  static int didfail;
  // check if failure.  If so, then grab restart from 2X ago (since
  // most previous restart may be too close to fixing problem), and
  // alter variables (cour) to fix.

  // assume failure doesn't occur in front of previous failure!


  if(firsttime){
    cour0=cour;
    //    cour0=0.9;

    // if no manual failure
    if(!restartonfail){
      didfail=0;
      beforenstep=0;
      beforetime=aftertime=0;
    }
    else{
      // use to start after failure and manual restart
      didfail=1;
      beforetime=911.944419026181;
      beforenstep=72425;
      aftertime=t;
      afternstep=realnstep;
    }
  }

  if((failed)&&(cour>1E-3)){
    beforenstep=realnstep;
    beforetime=t;
    // whichrestart is # of 2nd previous restart since next will use this #
    if (restart_init(whichrestart) >= 1) {
      dualfprintf(fail_file, "main:restart_init: failure\n");
    }
    afternstep=realnstep;
    aftertime=t;

    cour*=0.1;
    failed=0;
    failuremode=0;
    didfail=1;
  }
  //  if((!failed)&&(realnstep>beforenstep+100)){
  if(didfail&&((!failed)&&(t>beforetime+1)) ){
    cour=cour0;
    aftertime=0;
    afternstep=0;
    beforetime=0;
    beforenstep=0;
    didfail=0;
    trifprintf("Made it through failiure!\n");
  }
  if(failed) myexit(0); // if still failed, then end.

  // other option is to alter variables right at failure, never
  // letting this postdt get activated, but that's done elsewhere.
  firsttime=0;
}

void setfailresponse(int restartonfail)
{

  if(restartonfail==0){
    steptofaildump=1E30;
    steptofailmap=1E30;
    dofailmap=0;
    dofaildump=0;
  }
  else{
    // below is negative to turn off
    // or positive number of steps
    //    steptofaildump=-1;
    steptofaildump=1E30; // step number to start full dump each time step
    
    steptofailmap=nstep; // start right away
    
    // in absolute terms
    ifail=0;
    jfail=263;
    
    // check per CPU whether response will occur
    if(
       ((ifail>=startpos[1])||((mycpupos[1]==0)&&(ifail>=-N1BND))) &&
       ((ifail<=endpos[1])||((mycpupos[1]==ncpux1-1)&&(ifail<=totalsize[1]-1+N1BND))) &&
       ((jfail>=startpos[2])||((mycpupos[2]==0)&&(jfail>=-N2BND))) &&
       ((jfail<=endpos[2])||((mycpupos[2]==ncpux2-1)&&(jfail<=totalsize[2]-1+N2BND)))
       ){
      dofailmap=1;
      trifprintf("proc: %d will do fail areamap time series: absolutes: ifail=%d jfail=%d settofaildump=%d steptofailmap=%d\n",myid,ifail,jfail,steptofaildump,steptofailmap);
      // now set in relative terms
      ifail-=startpos[1];
      jfail-=startpos[2];
      trifprintf("proc: %d relative i/j: ifail=%d jfail=%d\n",myid,ifail,jfail);
    }
    else{
      dofailmap=0;
      trifprintf("proc: %d will NOT do fail areamap time series\n",myid);
    }
  }
}

