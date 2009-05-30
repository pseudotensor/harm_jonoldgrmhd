
#include "decs.h"

FTYPE slope_lim(FTYPE y1, FTYPE y2, FTYPE y3)
{
  FTYPE Dqm, Dqp, Dqc, s;

  /* woodward, or monotonized central, slope limiter */
  if (lim == MC) {
    Dqm = 2. * (y2 - y1);
    Dqp = 2. * (y3 - y2);
    Dqc = 0.5 * (y3 - y1);
    s = Dqm * Dqp;
    if (s <= 0.)
      return 0.;
    else {
      if (fabs(Dqm) < fabs(Dqp) && fabs(Dqm) < fabs(Dqc))
	return (Dqm);
      else if (fabs(Dqp) < fabs(Dqc))
	return (Dqp);
      else
	return (Dqc);
    }
  }
  /* van leer slope limiter */
  else if (lim == VANL) {
    Dqm = (y2 - y1);
    Dqp = (y3 - y2);
    s = Dqm * Dqp;
    if (s <= 0.)
      return 0.;
    else
      return (2. * s / (Dqm + Dqp));
  }
  /* minmod slope limiter (crude but robust) */
  else if (lim == MINM) {
    Dqm = (y2 - y1);
    Dqp = (y3 - y2);
    s = Dqm * Dqp;
    if (s <= 0.)
      return 0.;
    else if (fabs(Dqm) < fabs(Dqp))
      return Dqm;
    else
      return Dqp;
  } else {
    dualfprintf(fail_file, "unknown slope limiter\n");
    myexit(10);
  }
  return (0);
}
