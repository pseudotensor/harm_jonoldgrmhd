
/*
 *
 * generates initial conditions for a fishbone & moncrief disk 
 * with exterior at minimum values for density & internal energy.
 *
 * cfg 8-10-01
 *
 */

#include "decs.h"


#define PROGRADERISCO 0
#define RETROGRADERISCO 1


FTYPE taper_func(FTYPE R) ;
FTYPE nz_func(FTYPE R) ;
FTYPE rmso_calc(int which) ;
FTYPE uphi_isco_calc(int which,FTYPE r);
FTYPE rin ;

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
	cooling=1;

	BCtype[X1UP]=OUTFLOW;
	BCtype[X1DN]=OUTFLOW;
	BCtype[X2UP]=POLARAXIS;
	BCtype[X2DN]=POLARAXIS;

	a = 0.0 ;
	h_over_r = 0.05 ;
	rin = (1. + h_over_r)*rmso_calc(PROGRADERISCO) ;


	beta = 1.e2 ;

        /* some numerical parameters */
	failuremode = 0;		// clean start
        failed = 0 ;
        cour = 0.2 ;
        lim = MC ;
        dt = 1.e-5 ;
        Rin = 0.98*(1. + sqrt(1. - a*a)) ;
        Rout = 20. ;

	coord(-2,0,CENT,X) ;
	bl_coord(X,&r,&th) ;
	trifprintf("rin: %g Rin: %g\n",rin,Rin) ;
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
	ZSLOOP(0,N1-1,0,N2-1) {
		coord(i,j,CENT,X) ;
		bl_coord(X,&r,&th) ;

		/* region outside disk */
		R = r*sin(th) ;

		if(R < rin) {
			rho = 1.e-7*RHOMIN ;
                        u = 1.e-7*UUMIN ;

			// ks velocity
			get_geometry(i,j,CENT,&geom) ;
			ur = geom.gcon[0][1]/geom.gcon[0][0] ;
			uh = geom.gcon[0][2]/geom.gcon[0][0] ;
			up = geom.gcon[0][3]/geom.gcon[0][0] ;
		}
		else {
			H = h_over_r*R ;
			nz = nz_func(R) ;
			z = r*cos(th) ;
			S = 1./(H*H*nz) ;
			cs = H*nz ;			
			
			rho = (S/sqrt(2.*M_PI*H*H)) * exp(-z*z/(2.*H*H))
				* taper_func(R) ;
			u = rho*cs*cs/(gam - 1.) ;
			ur = 0. ;
			uh = 0. ;
			up = 1./(pow(r,1.5) + a) ;
			// solution in KS for v already
		}

		if(rho > rhomax) rhomax = rho ;
		if(u > umax) umax = u ;

		p[i][j][RHO] = rho ;
		p[i][j][UU] = u ;
		p[i][j][U1] = ur ;
		p[i][j][U2] = uh ;
		p[i][j][U3] = up ;

		p[i][j][B1] = 0. ;
		p[i][j][B2] = 0. ;
		p[i][j][B3] = 0. ;
	}

	mpimax(&rhomax);
	mpimax(&umax);
	trifprintf("rhomax: %g umax: %g\n", rhomax, umax);

	ZSLOOP(0,N1-1,0,N2-1) {
		p[i][j][RHO] /= rhomax ;
		p[i][j][UU]  /= rhomax ;
	}	
	umax /= rhomax ;
	rhomax = 1. ;

	fixup(p) ;
	if (bound_prim(p) >= 1)
	  FAILSTATEMENT("init.c:init()", "bound_prim()", 1);

	/* first find corner-centered vector potential */
	ZSLOOP(0,N1,0,N2) A[i][j] = 0. ;
        ZSLOOP(0,N1,0,N2) {
                /* field-in-disk version */
		/* flux_ct */

		coord(i,j,CORN,X) ;
		bl_coord(X,&r,&th) ;

		R = r*sin(th) ;

                u_av = 0.25*(
                        p[i][j][UU] +
                        p[i-1][j][UU] +
                        p[i][j-1][UU] +
                        p[i-1][j-1][UU]) ;
                u_ref = 0.25*(
                        p[i][N2/2][UU] +
                        p[i-1][N2/2][UU] +
                        p[i][N2/2-1][UU] +
                        p[i-1][N2/2-1][UU]) ;

		if(r > 1.1*rin) q = ((u_av/u_ref) - 0.2)*pow(r,0.25) ;
		else q = 0. ;

                if(q > 0.) A[i][j] = q*q*sin(X[1]/h_over_r)*taper_func(r) ;
		//trifprintf("%5d %5d %10.5g %10.5g %10.5g %10.5g\n",i,j,r,th,R,q) ;
        }

	/* now differentiate to find cell-centered B,
	   and begin normalization */
	bsq_max = 0. ;
	ZLOOP {
		get_geometry(i,j,CENT,&geom) ;

		/* flux-ct */
		p[i][j][B1] =  (A[i][j] - A[i][j+1] 
				+ A[i+1][j] - A[i+1][j+1])/(2.*dx[2]*geom.g) ;
		p[i][j][B2] = -(A[i][j] + A[i][j+1] 
				- A[i+1][j] - A[i+1][j+1])/(2.*dx[1]*geom.g) ;

		p[i][j][B3] = 0. ;

		if (bsq_calc(p[i][j], &geom, &bsq_ij) >= 1)
		  FAILSTATEMENT("init.c:init()", "bsq_calc()", 1);

		if(bsq_ij > bsq_max) bsq_max = bsq_ij ;
	}
	mpimax(&bsq_max);
	trifprintf("initial bsq_max: %g\n",bsq_max) ;

	/* finally, normalize to set field strength */
	beta_act = (gam - 1.)*umax/(0.5*bsq_max) ;
	trifprintf("initial beta: %g (should be %g)\n",beta_act,beta) ;
	norm = sqrt(beta_act/beta) ;
	bsq_max = 0. ;
	ZLOOP {
		p[i][j][B1] *= norm ;
		p[i][j][B2] *= norm ;

		get_geometry(i,j,CENT,&geom) ;


		if (bsq_calc(p[i][j], &geom, &bsq_ij) >= 1)
		  FAILSTATEMENT("init.c:init()", "bsq_calc()", 1);

		if(bsq_ij > bsq_max) bsq_max = bsq_ij ;
	}
	mpimax(&bsq_max);
	trifprintf("new initial bsq_max: %g\n", bsq_max);
	beta_act = (gam - 1.)*umax/(0.5*bsq_max) ;
	trifprintf("final beta: %g (should be %g)\n",beta_act,beta) ;

#if 0
#endif

	/* enforce boundary conditions */
	fixup(p) ;
	if (bound_prim(p) >= 1)
	  FAILSTATEMENT("init.c:init()", "bound_prim()", 2);

	return(0);
}

FTYPE taper_func(FTYPE R)
{

	if(R <= rin) 
		return(0.) ;
	else 
		return(1. - sqrt(rin/R)) ;

}

FTYPE nz_func(FTYPE R)
{
	return(
		sqrt(
		(3.*a*a - 4.*a*sqrt(R) + R*R)/
		pow(R*(a + pow(R,1.5)),2)
		)
	) ;


}

// compute the radius of the inner most stable circular orbit
FTYPE rmso_calc(int which)
{
        FTYPE rmso,Z1,Z2,sign ;


	if(which==PROGRADERISCO) sign=1; else sign=-1;

        Z1 = 1. + pow(1. - a*a,1./3.)*(pow(1. + a,1./3.) +
                pow(1. - a, 1./3.)) ;
        Z2 = sqrt(3.*a*a + Z1*Z1) ;
	rmso=3. + Z2-sign*sqrt((3. - Z1)*(3. + Z1 + 2.*Z2)) ;

	return(rmso) ;
}

FTYPE uphi_isco_calc(int which,FTYPE r)
{
  FTYPE uphi;
  FTYPE sign;
  FTYPE Z1,Z2;

  if(which==PROGRADERISCO) sign=1; else sign=-1;

  Z1=r*r-sign*2.*a*sqrt(r)+a*a;
  Z2=r*(r*r-3.*r+sign*2.*a*sqrt(r));

  uphi=sign*Z1/sqrt(Z2);

  return(uphi);

}

