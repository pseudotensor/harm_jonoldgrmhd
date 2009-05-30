
/*
 *
 * generates initial conditions for a fishbone & moncrief disk 
 * with exterior at minimum values for density & internal energy.
 *
 * cfg 8-10-01
 *
 */

#include "decs.h"


int init()
{
	int i,j,k ;
	FTYPE r,th,sth,cth ;
	FTYPE ur,uh,up,u,rho ;
	FTYPE X[NDIM] ;
	struct of_geom geom ;

	/* for disk interior */
	FTYPE R,H,nz,z,S,cs ;

	/* for magnetic field */
	FTYPE A[N1+1][N2+1] ;
	FTYPE rho_av,rhomax,umax,beta,bsq_ij,bsq_max,norm,q,beta_act ;
	FTYPE u_av,u_ref ;
       
	FTYPE ftemp;

	/* some physics parameters */
	defcoord=0;
	R0=0;
	gam = 4./3. ;
	cooling=0;
	
	BCtype[X1UP]=FIXED;
	BCtype[X1DN]=OUTFLOW;
	BCtype[X2UP]=POLARAXIS;
	BCtype[X2DN]=POLARAXIS;


	a = 0.0 ;
	h_over_r = 1.0;

        /* some numerical parameters */
	failuremode = 0;		// clean start
        failed = 0 ;
        cour = 0.9 ;
        lim = MC ;
        dt = 1.e-5 ;
        Rin = 0.98*(1. + sqrt(1. - a*a)) ;
        Rout = 20. ;

	coord(-2,0,CENT,X) ;
	bl_coord(X,&r,&th) ;
	trifprintf("Rin: %g\n",Rin) ;
	trifprintf("rmin: %g\n",r) ;
	trifprintf("rmin/rm: %g\n",r/(1. + sqrt(1. - a*a))) ;

        t = 0. ;
        hslope = h_over_r ;

	// check for square grid
	/*
	ftemp=M_PI*hslope/(pow(Rout/Rin,1.0/N1)-1.0);
	if(fabs(ftemp-N2)>=1){
	  trifprintf("grid not square.  Expected N2=%d, got N2=%d\n",ftemp,N2);
	  exit(0);
	}
	*/
        set_grid() ;

        /* output choices */
	tf = 4000.0 ;

	DTd = 50. ;	/* dumping frequency, in units of M */
	DTener = 2. ;	/* logfile frequency, in units of M */
	DTi = 2. ; 	/* image file frequ., in units of M */
	DTr = 512 ; 	/* restart file frequ., in timesteps */

	rhomax = 0. ;
	umax = 0. ;
	ZSLOOP(-2,N1+1,-2,N2+1) {
		coord(i,j,CENT,X) ;
		bl_coord(X,&r,&th) ;


		rho = 10.*RHOMIN ;
		u = 10.*UUMIN ;
		
		// ks velocity
		get_geometry(i,j,CENT,&geom) ;
		ur = geom.gcon[0][1]/geom.gcon[0][0] ;
		uh = geom.gcon[0][2]/geom.gcon[0][0] ;
		up = geom.gcon[0][3]/geom.gcon[0][0] ;

		p[i][j][RHO] = rho ;
		p[i][j][UU] = u ;
		p[i][j][U1] = ur ;
		p[i][j][U2] = uh ;
		p[i][j][U3] = up ;

		p[i][j][B1] = 0. ;
		p[i][j][B2] = 0. ;
		p[i][j][B3] = 0. ;
	}

	// for fixed boundary conditions
	ZSLOOP(-2,N1+1,-2,N2+1) PLOOP{
	  ph[i][j][k]=p[i][j][k];
	}

	fixup(p) ;
	if (bound_prim(p) >= 1)
	  FAILSTATEMENT("init.c:init()", "bound_prim()", 1);


	return(0);
}

