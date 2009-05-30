
// #include "decs.h"
#include "defs.h"

int main(int argc, char *argv[])
{

  fprintf(stderr,"Start init\n"); fflush(stderr);

  /* perform initializations */
  pre_init(argc, argv);
  // all log files are now defined

  if (RESTARTMODE == 1) {
    if (restart_init(WHICHFILE) >= 1) {
      dualfprintf(fail_file, "main:restart_init: failure\n");
    }
  } else if (init() >= 1) {
    dualfprintf(fail_file, "main:init: failure\n");
  }

  post_init();			// initialize general parameters

  trifprintf("proc: %04d : End init\n", myid);


  /* Do initial diagnostics */
  if (DODIAGS) {
    trifprintf("proc: %04d : Start initial diagnostics\n", myid);
    // no error_check since if init passed, diag(0) should pass
    diag(0);
    trifprintf("proc: %04d : End initial diagnostics\n", myid);
  }

  trifprintf("proc: %04d : Start computation\n", myid);

  while (t < tf) {


    /* step variables forward in time */
    nstroke = 0;
    step_ch();

    // must check before MPI operation (since asymmetries would
    // desynchronize cpus)
    if (error_check()) {
      fprintf(fail_file, "error_check detected failure at main:2\n");
      fflush(fail_file);
    }
    //^^ otherwise ok
    
    // eventually all cpus come here, either in failure mode or not,
    // and cleanly tell others if failed/exit/dump/etc.

    postdt(); // here one can alter variables and try to restart, or implement any post step operations

    /* perform diagnostics */
    // no error check since assume if step_ch passed, diag(1) will pass
    if (DODIAGS)
      diag(1);

    /* restart dump */
    // if(nstep == 130) restart_write(1) ;

    nstep++;
    // restartsteps[whichrestart]=realnstep;


    trifprintf("%21.15g %21.15g %10.5g %8ld %8ld %8d\n", t, dt,
	      cour, nstep, realnstep, nstroke);
  }
  trifprintf("proc: %04d : End computation\n", myid);

  /* do final diagnostics */
  if (DODIAGS)
    diag(2);


  trifprintf("ns,ts: %ld %ld\n", nstep, nstep *(long)( totalsize[1] * totalsize[2]));

  myexit(0);
  return (0);
}
