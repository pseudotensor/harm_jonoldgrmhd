#include "decs.h"


#define DEBUG (1)

int mnewt(int ntrial, FTYPE x[], int n, FTYPE tolx, FTYPE tolf)
{
  int i = 0, j = 0, k = 0;
  FTYPE errx, errf, d;
  static int firstc = 1;
  static int *indx;
  static FTYPE **fjac, *fvec, *pp;
  // debug stuff
  static int count = 0;
  static long lastnstep = 0;
  static int calls = 0;

#if(DEBUG)
  calls++;
#endif
  if (firstc) {
    firstc = 0;

    indx = ivector(1, n);
    pp = dvector(1, n);
    fvec = dvector(1, n);
    fjac = dmatrix(1, n, 1, n);
  }

  for (k = 1; k <= ntrial; k++) {
    nstroke++;
    if (usrfun(x, n, fvec, fjac) >= 1) {
      dualfprintf(fail_file, "mnewt:usrfun: (k=%d) failure\n", k);
      return (1);
    }
    errf = 0.0;
    for (i = 1; i <= n; i++)
      errf += fabs(fvec[i]);
    if (errf <= tolf) {
#if(DEBUG)
      if (lastnstep < nstep) {
	trifprintf("#1 count/zone: %g calls: %d\n",
		(FTYPE) count / ((FTYPE) (N1 * N2)),
		calls / (N1 * N2));
	count = k - 1;
	lastnstep = nstep;
	calls = 0;
      } else {
	count += k - 1;
      }
#endif
      return (0);
    }
    for (i = 1; i <= n; i++)
      pp[i] = -fvec[i];
    ludcmp(fjac, n, indx, &d);
    lubksb(fjac, n, indx, pp);
    errx = 0.0;
    for (i = 1; i <= n; i++) {
      errx += fabs(pp[i]);
      x[i] += pp[i];
    }
    if (errx <= tolx) {
#if(DEBUG)
      if (lastnstep < nstep) {
	trifprintf("#2 count/zone: %g calls: %d\n",
		(FTYPE) count / ((FTYPE) (N1 * N2)),
		calls / (N1 * N2));
	fflush(log_file);
	count = k;
	lastnstep = nstep;
	calls = 0;
      } else {
	count += k;
      }
#endif
      return (0);


    }


  }
  trifprintf("proc: %d, mnewt didn't converge: i=%d j=%d, errf: %g errx: %g\n",myid, startpos[1]+icurr,startpos[2]+jcurr,errf,errx);
  if ((errf <= 1E-4)&&(errx<=1E-4)) {
    return (0); // for now
    // assume not too bad convergence if 1E-4
  }
  else{
    failuremode = 3;		// source of failure (nonconvergence)
    FAILSTATEMENT("mnewt.c", "convergence", 1);
  }
}
