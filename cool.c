#include "decs.h"

/* JON: here's the cooling function, w/ two parameters */

#define THETACOOL       (h_over_r)	/* should be same as h_over_r */
#define TAUCOOL         (5.0)	        /* cooling time */
#define NOCOOLTHFACT    (3.0)           /* this times h_over_r and no more cooling there*/
#define COOLTAPER1(th)   (exp(-pow((th)-M_PI*0.5,2.0)/(2.0*pow(NOCOOLTHFACT*h_over_r,2.0)));
#define SLOPE2 (1.0/(M_PI*0.5-NOCOOLTHFACT*h_over_r))
#define COOLTAPER2(th)   (((th)<M_PI*0.5-NOCOOLTHFACT*h_over_r) ? (SLOPE2*(th)) : ( ((th)>M_PI*0.5+NOCOOLTHFACT*h_over_r) ? (-SLOPE2*((th)-M_PI)) : 1.0 ) )
#define SLOPE3 (1.0/(NOCOOLTHFACT*h_over_r))
#define WIDTHTAPER (NOCOOLTHFACT*h_over_r)
#define TAPERPOS1 (M_PI*0.5-NOCOOLTHFACT*h_over_r-WIDTHTAPER)
#define TAPERPOS2 (M_PI*0.5-NOCOOLTHFACT*h_over_r)
#define TAPERPOS3 (M_PI*0.5+NOCOOLTHFACT*h_over_r)
#define TAPERPOS4 (M_PI*0.5+NOCOOLTHFACT*h_over_r+WIDTHTAPER)
#define TAPERFUN1(th) (0.0)
#define TAPERFUN2(th) (SLOPE3*((th)-TAPERPOS1))
#define TAPERFUN3(th) (1.0)
#define TAPERFUN4(th) (-SLOPE3*((th)-TAPERPOS4))
#define TAPERFUN5(th) (0.0)
#define COOLTAPER3(th)   (((th)<TAPERPOS1 ? TAPERFUN1(th) : (  ((th)<TAPERPOS2) ? TAPERFUN2(th) : (  ((th)<TAPERPOS3) ? TAPERFUN3(th): (  ((th)<TAPERPOS4) ? TAPERFUN4(th) : TAPERFUN5(th)   )))) )
/* cooling function, if any */

FTYPE coolfunc(FTYPE h_over_r, FTYPE *ph, struct of_geom *geom, struct of_state *q)
{
        FTYPE X[NDIM],r,th,R,Wcirc,cs_circ,rho,u,P,w,wcirc,dUcool ;
	FTYPE taper0;

        /* cooling function for maintaining fixed H/R */
        rho = ph[RHO] ;
        u = ph[UU] ;
        P = (gam - 1.)*u ;
        w = rho + u + P ;

        coord(icurr,jcurr,CENT,X) ;
        bl_coord(X,&r,&th) ;
        R = r*sin(th) ;

        /* crude approximation */
        Wcirc = pow(R,-1.5) ;
        cs_circ = THETACOOL/sqrt(R) ;
        wcirc = rho*(1. + cs_circ*cs_circ/(gam - 1.)) ;

        if(t > 0.){
	  dUcool = -(Wcirc/TAUCOOL)*( (w - wcirc)*(q->ucon[TT])*(q->ucov[TT])) ;
	  // shape function to avoid problems near pole
	  //taper0=COOLTAPER(0);
	  //dUcool*=1./(1.-1./taper0)+1./(1.-taper0)*COOLTAPER(th);
	  //dUcool*=COOLTAPER2(th);
	  dUcool*=COOLTAPER3(th);
	}
        else{
	  dUcool = 0. ;
	}

        return(dUcool) ;
}
