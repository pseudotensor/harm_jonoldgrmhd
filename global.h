#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>

/* size of global arrays */
#define N1	128		/* number of zones */
#define N2	64       	/* number of zones */

/** MNEMONICS SECTION **/

/* boundary condition mnemonics */
#define OUTFLOW	0
#define SYMM	1
#define ASYMM	2
#define FIXED	3
#define POLARAXIS 4

/* mnemonics for primitive vars; conserved vars */
#define RHO	0
#define UU	1
#define U1	2
#define U2	3
#define U3	4
#define B1	5
#define B2	6
#define B3	7

/* mnemonics for dimensional indices */
#define TT	0
#define RR	1
#define TH	2
#define PH	3

/* mnemonics for centering of grid functions */
#define NUMGRIDPOS 4
#define FACE1	0
#define FACE2	1
#define CORN	2
#define CENT	3

/* mnemonics for slope limiter */
#define MC	0
#define VANL	1
#define MINM	2

/* choices for algorithm */
#define FLUXCT	1
#define FLUXCD	0

/* mnemonics for diagnostic calls */
#define INIT_OUT	0
#define DUMP_OUT	1
#define IMAGE_OUT	1
#define LOG_OUT		1
#define FINAL_OUT	2


/** GLOBAL ARRAY SECTION **/

/* size of global arrays */
#define NPR	8		/* number of primitive variables */
#define NDIM	4		/* number of total dimensions.  Never
				   changes */
#define NPG	4		/* number of positions on grid for grid 
				   functions */

#define NBIG ((N1>N2) ? N1 : N2)


// size of data type used for all floats
#define FLOATTYPE 0
#define DOUBLETYPE 1


// define your user type here
// (normal non-sensitive or performance critical datatypes)
#define REALTYPE DOUBLETYPE
 // (non-perf critical or sensitive data types) 
#define SENSITIVE DOUBLETYPE
// WE ASSUME SENSITIVE>=REALTYPE !


#ifndef FLT_EPSILON
#define FLT_EPSILON (1.E-7)
#endif
#ifndef DBL_EPSILON
#define DBL_EPSILON (2.E-16)
#endif


// need not change below datatype stuff
#if(REALTYPE==FLOATTYPE)
#define NUMEPSILON FLT_EPSILON
#define FTYPE float
#else
#define NUMEPSILON DBL_EPSILON
#define FTYPE double
#endif
#if(SENSITIVE==FLOATTYPE) // for sensitive counters
#define SFTYPE float
#else
#define SFTYPE double
#endif

#define NUMENERVAR (6+NPR+NPR+3)

/* numerical convenience */
#define SMALL	1.e-20

#define SLEPSILON (1.e-6)


/* size of step in numerical derivative evaluations */
#define HSTEP	1.e-5


#define SURFACETOTAL 0
#define VOLUMETOTAL 1

#define CONSTYPE 0
#define SURFACETYPE 1
#define CUMULATIVETYPE 2

#define WRITEHEAD 0
#define READHEAD 1

#define TIMESERIESAREAMAP 0
#define FINALTDUMPAREAMAP 1

#define WRITEFILE 0
#define READFILE 1

#define ENERFNAME "ener.out"


#if(SENSITIVE==DOUBLETYPE)

#if(REALTYPE==DOUBLETYPE)
#define RESTARTHEADER "%d %d "\
	      "%lf %lf %ld %lf %lf %lf "\
	      "%lf %lf %lf %ld %ld %ld %ld "\
	      "%lf %d %d "\
	      "%lf %lf %lf %lf %d"\
              "%d %d %d %d"
#elif(REALTYPE==FLOATTYPE)
#define RESTARTHEADER "%d %d "\
	      "%lf %lf %ld %f %f %f "\
	      "%lf %lf %lf %ld %ld %ld %ld "\
	      "%lf %d %d "\
	      "%f %f %f %f %d"\
              "%d %d %d %d"
#endif

#elif(SENSITIVE==FLOATTYPE)

#if(REALTYPE==DOUBLETYPE)
#define RESTARTHEADER "" // dumb, so crash on compile
#elif(REALTYPE==FLOATTYPE)
#define RESTARTHEADER "%d %d "\
	      "%f %f %ld %f %f %f "\
	      "%f %f %f %ld %ld %ld %ld "\
	      "%f %d %d "\
	      "%f %f %f %f %d"\
              "%d %d %d %d"
#endif

#endif



#define WRITERESTARTHEADER "%d %d " \
		 "%21.15g %21.15g %ld %21.15g %21.15g %21.15g " \
		 "%21.15g %21.15g %21.15g %ld %ld %ld %ld " \
		 "%21.15g %d %d " \
		 "%21.15g %21.15g %21.15g %21.15g %d " \
                 "%d %d %d %d\n"

// now that all hashes have been defined, get mpi header
#include "mympi.h"


/** MACROS **/

/* loop over all active zones */
#define ZLOOP for(i=0;i<N1;i++)for(j=0;j<N2;j++)

/* specialty loop */
#define ZSLOOP(istart,istop,jstart,jstop) \
	for(i=istart;i<=istop;i++)\
	for(j=jstart;j<=jstop;j++)

#define GENLOOP(i,j,istart,istop,jstart,jstop) \
        for((i)=(istart);(i)<=(istop);(i)++)\
        for((j)=(jstart);(j)<=(jstop);(j)++)


/* specialty loop */
#define LOOPDIVB ZSLOOP(-1,N1+1,-1,N2+1)

// #define DIVBCONDITION(p,i,j)
// if((i>=-1)&&(j>=-1)&&(startpos[2]+j!=0)&&(startpos[2]+j!=N2TOT))
#define DIVBCONDITION(p,i,j) if((startpos[1]+i>0)&&(startpos[2]+j>0)&&(startpos[1]+i<totalsize[1])&&(startpos[2]+j<totalsize[2]))

#if(FLUXCT)
#define DIVB(p,i,j) ((0.5*(\
			      p[i][j][B1]*gdet[i][j][CENT] \
			      + p[i][j-1][B1]*gdet[i][j-1][CENT]\
			      - p[i-1][j][B1]*gdet[i-1][j][CENT] \
			      - p[i-1][j-1][B1]*gdet[i-1][j-1][CENT]\
			      )/dx[1] +\
		      0.5*(\
			      p[i][j][B2]*gdet[i][j][CENT] \
			      + p[i-1][j][B2]*gdet[i-1][j][CENT]\
			      - p[i][j-1][B2]*gdet[i][j-1][CENT] \
			      - p[i-1][j-1][B2]*gdet[i-1][j-1][CENT]\
			      )/dx[2])/gdet[i][j][CORN])
// charles has no final division by gdet[i][j][CORN]
#elif(FLUXCD)
#define DIVB(p,i,j)  ((0.5*(\
				p[i+1][j][B1]*gdet[i+1][j][CENT] \
				- p[i-1][j][B1]*gdet[i-1][j][CENT]\
				)/dx[1] +\
		       0.5*(\
				p[i][j+1][B2]*gdet[i][j+1][CENT] \
				- p[i][j-1][B2]*gdet[i][j-1][CENT] \
				)/dx[2]))

#endif
// poles defined as divb=0, can't divide due to singularity (could use
// volume regularization)
#define SETFDIVB(divb,p,i,j) {DIVBCONDITION(p,i,j){ divb = fabs(DIVB(p,i,j)) ;} else divb = 0.;}

/* want dump output to be ordered in radius first!! */
#define DUMPLOOP(istart,istop,jstart,jstop) \
	for(j=jstart;j<=jstop;j++)\
	for(i=istart;i<=istop;i++)

#define IMAGELOOP(istart,istop,jstart,jstop) \
	for(j=jstart;j<=jstop;j++)\
	for(i=istart;i<=istop;i++)

#define OLDIMAGELOOP for(j=N2-1;j>=0;j--) for(i=0;i<N1;i++)	// nasty 
								// to
								// deal 
								// with

/* loop over Primitive variables */
#define PLOOP for(k=0;k<NPR;k++)
/* loop over all Dimensions; second rank loop */
#define DLOOP for(j=0;j<NDIM;j++)for(k=0;k<NDIM;k++)
/* loop over all Dimensions; first rank loop */
#define DLOOPA for(j=0;j<NDIM;j++)
/* loop over all Space dimensions; second rank loop */
#define SLOOP for(j=1;j<NDIM;j++)for(k=1;k<NDIM;k++)
/* loop over all Space dimensions; first rank loop */
#define SLOOPA for(j=1;j<NDIM;j++)
/* loop over all for j and Space for k; second rank loop */
#define DSLOOP for(j=0;j<NDIM;j++)for(k=1;k<NDIM;k++)
/* loop over all for k and Space for j; second rank loop */
#define SDLOOP for(j=1;j<NDIM;j++)for(k=0;k<NDIM;k++)

#define DIRLOOP for(dir=0;dir<COMPDIM*2;dir++)

#define MYDMIN(a,b) (mydminarg1=(a),mydminarg2=(b),(mydminarg1) < (mydminarg2) ?\
        (mydminarg1) : (mydminarg2))

#define delta(i,j) ((i == j) ? 1. : 0.)
#define dot(a,b) (a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3])

#define mink(I,J) (I != J ? (0.) : (I == 0 ? (-1.) : (1.)))

#define pfixupeach(pr,i,j,which,min) {if(pr[which]<min){ fladd[which]+=dV*gdet[i][j][CENT]*(min-pr[which]); pr[which]=min;}}

#define pfixup(pr,i,j) {pfixupeach(pr,i,j,RHO,RHOMIN); pfixupeach(pr,i,j,UU,UUMIN); }

// #define FAILSTATEMENT(file,function,number) {fprintf(fail_file,"%s
// %d-%s(): failure\n",file,number,function); fflush(fail_file);
// fprintf(fail_file,"rho[i][j]: %15.10g uu[i][j]: %15.10g rho2[i][j]:
// %15.10g uu2[i][j]: %15.10g i: %d j: %d k:
// %d\n",p[i][j][RHO],p[i][j][UU],ph[i][j][RHO],ph[i][j][UU],i,j,k);
// return(1);}

#define FAILSTATEMENT(file,function,number) {fprintf(fail_file,"%s %d-%s(): failure\n",file,number,function); fflush(fail_file); fprintf(fail_file,"i: %d j: %d k: %d\n",i,j,k); return(1);}

#define im1 i-1
#define im1mac(i) i-1
#define ip1 i+1
#define ip1mac(i) i+1


#define jm1 j-1
#define jm1mac(j) j-1
#define jp1 j+1
#define jp1mac(j) j+1

/* failure modes */
#define FAIL_UTOPRIM_NEG	1
#define FAILSTR01 "UTOPRIM_NEG"
#define FAIL_UTOPRIM_TEST	2
#define FAILSTR02 "UTOPRIM_TEST"
#define FAIL_VCHAR_DISCR	3
#define FAILSTR03 "VCHAR_DISCR"
#define FAIL_COEFF_NEG		4
#define FAILSTR04 "COEFF_NEG"
#define FAIL_COEFF_SUP		5
#define FAILSTR05 "COEFF_SUP"
#define FAIL_UTCALC_DISCR	6
#define FAILSTR06 "UTCALC_DISCR"

// structure declarations
/* set global variables that indicate current local metric, etc. */
struct of_geom {
  FTYPE gcon[NDIM][NDIM];
  FTYPE gcov[NDIM][NDIM];
  FTYPE g;
};
struct of_state {
  FTYPE ucon[NDIM];
  FTYPE ucov[NDIM];
  FTYPE bcon[NDIM];
  FTYPE bcov[NDIM];
};

// function declarations
extern int main(int argc, char *argv[]);
extern int error_check(void);
extern void find_horizon(void);
extern void restart_write(int which);
extern int dump(int dump_cnt);
extern void gdump(void);
extern int primtoU(FTYPE *p, struct of_state *q, struct of_geom *geom,
		   FTYPE *U);
extern void image_dump(int image_cnt);
extern int init(void);
extern FTYPE coolfunc(FTYPE h_over_r, FTYPE *ph, struct of_geom *geom,
		       struct of_state *q);
extern void post_init(void);
extern void pre_init(int argc, char *argv[]);
extern int step_ch(void);
extern void postdt(void);
extern int advance(FTYPE pi[][N2 + 4][NPR], FTYPE pb[][N2 + 4][NPR],
		   FTYPE Dt, FTYPE pf[][N2 + 4][NPR], FTYPE *ndt);
extern int fluxcalc(FTYPE pr[][N2 + 4][NPR], FTYPE F[][N2 + 4][NPR],
		    int dir, FTYPE *ndt);
extern void flux_cd(FTYPE F1[][N2 + 4][NPR], FTYPE F2[][N2 + 4][NPR]);
extern int diag(int call_code);
extern void diag_flux(FTYPE F1[][N2 + 4][NPR],
		      FTYPE F2[][N2 + 4][NPR],SFTYPE Dt);
extern void frdotout(void);
extern int restart_init(int which);
extern void read_restart_header(char *dfnam, FILE* headerptr);
extern void write_restart_header(char *dfnam, FILE* headerptr);
extern void makedirs(void);
extern void myfopen(char*fname, char*fmt, char*message, FILE ** fileptr);
extern void myfclose(FILE ** fileptr,char*message);
extern void mydfwrite(FTYPE *ptr, int start, size_t nmemb, int i, int j, FILE*stream, FTYPE*writebuf);
extern void mydfread(FILE * stream, int start, size_t nmemb, int i, int j, FTYPE*writebuf, FTYPE *ptr);
extern void appendener(FILE* ener_file,SFTYPE pdot_tot[][NPR],SFTYPE*fladd_tot);
extern void divbmaxavg(FTYPE p[][N2+4][NPR],SFTYPE*ptrdivbmax,SFTYPE*ptrdivbavg);
extern void gettotal(int numvars, SFTYPE* vars[],int*sizes,SFTYPE*vars_tot[]);
extern int constotal(SFTYPE *vars_tot);
extern int integrate(SFTYPE * var,SFTYPE *var_tot,int type);

extern int ucon_calc(FTYPE *pr, struct of_geom *geom, FTYPE *ucon);
extern int ucon_calcother(FTYPE *pr, FTYPE *ucon);
extern void ucon_precalc(FTYPE *ucon, FTYPE *AA, FTYPE *BB,
			 FTYPE *CC, FTYPE *discr);
extern int ucon_fix(FTYPE disc, FTYPE AA, FTYPE BB, FTYPE CC,
		    FTYPE *ucon);
extern void blgset(int i, int j, struct of_geom *geom);
extern void cleanfield(FTYPE prim[][N2 + 4][NPR]);
extern void b_calc(FTYPE *pr, FTYPE *ucon, FTYPE *b);

extern void dutdui_calc(FTYPE *ucon, FTYPE *dutdui);
extern void duiduj_calc(FTYPE *ucon, FTYPE *dutdui);
extern void dbtdui_calc(FTYPE *dutdui, FTYPE *pr, FTYPE *dbtdui);
extern void dbiduj_calc(FTYPE *dbtdui, FTYPE *dutdui, FTYPE *ucon,
			FTYPE *b, FTYPE dbiduj[][NDIM]);
extern void db2dui_calc(FTYPE dbiduj[][NDIM], FTYPE *b,
			FTYPE *db2dui);
extern void duudud_calc(FTYPE *ucon, FTYPE duudud[][NDIM]);

extern void dbsqdui_calc(FTYPE dbiduj[][NDIM], FTYPE *b,
			 FTYPE *dbsqdui);

extern FTYPE contract(FTYPE *vcon, FTYPE *wcon);
extern int sp_stress_calc(FTYPE *pr, FTYPE tens_matt[][NDIM],
			  FTYPE tens_em[][NDIM], FTYPE *b,
			  FTYPE *ucon);
extern void image(int image_cnt, int which, int scale, int limits);
extern void prminmaxsum(FTYPE p[][N2+4][NPR], int start,int nmemb, FTYPE *max, FTYPE*min,FTYPE*sum);
extern void myfprintf(FILE* fileptr, char *format, ...);
extern void dualfprintf(FILE* fileptr,char *format, ...);
extern void logsfprintf(char *format, ...);
extern void trifprintf(char *format, ...);
extern void set_arrays(void), set_grid(void);
extern int bound_prim(FTYPE prim[][N2 + 4][NPR]);
extern int bsq_calc(FTYPE *pr, struct of_geom *geom, FTYPE *b2);
extern int bltomet(FTYPE *pr, int i, int j);
extern int bl2met2metp2v(FTYPE *pr, int i, int j);
extern FTYPE ranc(int seed);
extern FTYPE bl_gdet_func(FTYPE r, FTYPE th);
extern FTYPE gdet_func(FTYPE lgcov[][NDIM]);
extern void bl_gcov_func(FTYPE r, FTYPE th, FTYPE gcov[][NDIM]);
extern void bl_gcon_func(FTYPE r, FTYPE th, FTYPE gcon[][NDIM]);
extern void coord(int i, int j, int loc, FTYPE *X);
extern void bl_coord(FTYPE *X, FTYPE *r, FTYPE *th);
extern void mettometp(FTYPE *pr, int i, int j);

extern int Utoprim(FTYPE *U, struct of_geom *geom, FTYPE *pr);
extern int Utoprim_ldz(FTYPE *U, struct of_geom *geom, FTYPE *pr);
extern void wvsq_solv_ldz(FTYPE *vsq, FTYPE *W);
extern FTYPE nrunsafe(void (*funcd) (FTYPE, FTYPE*,FTYPE*), FTYPE guess);
extern void func(FTYPE x, FTYPE *f, FTYPE *df);
extern FTYPE rtsafe(void (*funcd) (), FTYPE x1, FTYPE x2,
		     FTYPE xacc);


extern void restart_read(int which);
extern void conn_func(FTYPE *X, struct of_geom *geom,
		      FTYPE lconn[][NDIM][NDIM]);
extern void gcov_func(FTYPE *X, FTYPE lgcov[][NDIM]);
extern void gcon_func(FTYPE lgcov[][NDIM], FTYPE lgcon[][NDIM]);
extern void tetr_func(FTYPE tetr_cov[][NDIM], FTYPE tetr_con[][NDIM]);
extern void get_geometry(int i, int j, int loc, struct of_geom *geom);
extern int get_state(FTYPE *pr, struct of_geom *geom,
		     struct of_state *q);
int primtoflux(FTYPE *pa, struct of_state *q, int dir,
	       struct of_geom *geom, FTYPE *fl);
extern int source(FTYPE *pa, struct of_geom *geom, int ii, int jj,
		  FTYPE *Ua);
extern void fixup(FTYPE (*var)[N2 + 4][NPR]);
extern FTYPE slope_lim(FTYPE y1, FTYPE y2, FTYPE y3);
extern void tet_func(FTYPE metr[][NDIM], FTYPE tetr[][NDIM]);
extern int dsyev_(char *jobz, char *uplo, int *n, FTYPE *a, int *lda,
		  FTYPE *w, FTYPE *work, int *lwork, int *iwork,
		  int *liwork, int *info);
extern int dudp_calc(FTYPE *pr, struct of_state *q,
		     struct of_geom *geom, FTYPE **alpha);
extern int fail(int fail_type);
extern void setfailresponse(int restartonfail);
extern void setflux(void);
extern void setrestart(long *dump_cnt, long*image_cnt,long *rdump_cnt,int*appendold);
extern void flux_ct(FTYPE F1[][N2 + 4][NPR], FTYPE F2[][N2 + 4][NPR]);
extern int vchar(FTYPE *pr, struct of_state *q, int dir,
		 struct of_geom *geom, FTYPE *cmax, FTYPE *cmin);
extern FTYPE chk_disp(FTYPE v);
extern void make_co_to_comov(FTYPE *ucon, FTYPE ecov[][NDIM],
			     FTYPE econ[][NDIM]);
extern void transform(FTYPE *vec, FTYPE t[][NDIM]);
extern void coeff_set(FTYPE rho, FTYPE u);
extern void transform(FTYPE *ucon, FTYPE t[][NDIM]);
extern void ludcmp(FTYPE **a, int n, int *indx, FTYPE *d);
extern void lubksb(FTYPE **a, int n, int *indx, FTYPE *d);
extern void mhd_calc(FTYPE *pr, int dir, struct of_state *q,
		     FTYPE *mhd);

extern int mnewt(int ntrail, FTYPE *p, int n, FTYPE tolx,
		 FTYPE tolf);
extern int usrfun(FTYPE *pr, int n, FTYPE *beta, FTYPE **alpha);
extern FTYPE zbrent(FTYPE (*func) (FTYPE), FTYPE v1, FTYPE v2,
		     FTYPE tol);
extern int area_map(int call_code, int type, int size, int i, int j, FTYPE prim[][N2 + 4][NPR]);
extern void bcon_calc(FTYPE *pr, FTYPE *ucon, FTYPE *ucov,
		      FTYPE *bcon);
// NR STUFF

/* NR routines from nrutil.h */
extern int *ivector(long nl, long nh);
extern void free_ivector(int *v, long nl, long nh);
extern FTYPE *dvector(long nl, long nh);
extern void free_dvector(FTYPE *v, long nl, long nh);
extern FTYPE **dmatrix(long nrl, long nrh, long ncl, long nch);
extern void free_dmatrix(FTYPE **m, long nrl, long nrh, long ncl,
			 long nch);
extern FTYPE ***dtensor(long nrl, long nrh, long ncl, long nch,
			 long ndl, long ndh);
extern void free_dtensor(FTYPE ***t, long nrl, long nrh, long ncl,
			 long nch, long ndl, long ndh);
extern void nrerror(char error_text[]);

/* specialty functions */
extern void bondi_solve(FTYPE K, FTYPE gam, FTYPE *Rs, FTYPE *Urs,
			FTYPE *Edot);
extern FTYPE bondi_trace(FTYPE K, FTYPE gam, FTYPE edotf, FTYPE r,
			  FTYPE rs, FTYPE urs);
extern void timestep(FTYPE ndtr, FTYPE ndth);
extern FTYPE dtset(FTYPE ndtr, FTYPE ndth);
extern void u_to_v(FTYPE *pr, int i, int j);
extern FTYPE bondi_trace(FTYPE K, FTYPE gam, FTYPE edotf,
			  FTYPE r, FTYPE rs, FTYPE urs);
extern void bondi_solve(FTYPE K, FTYPE gam, FTYPE *Rs,
			FTYPE *Urs, FTYPE *Edot);
extern FTYPE edot_calc(FTYPE r, FTYPE ur, FTYPE g, FTYPE K);
extern FTYPE dedr_calc(FTYPE r, FTYPE ur, FTYPE g, FTYPE K);
extern FTYPE dedur_calc(FTYPE r, FTYPE ur, FTYPE g, FTYPE K);
extern FTYPE d2edr2_calc(FTYPE r, FTYPE ur, FTYPE g, FTYPE K);
extern FTYPE d2edur2_calc(FTYPE r, FTYPE ur, FTYPE g, FTYPE K);
extern FTYPE d2edrdur_calc(FTYPE r, FTYPE ur, FTYPE g, FTYPE K);

extern void lower(FTYPE *a, struct of_geom *geom, FTYPE *b);
extern void raise(FTYPE *v1, struct of_geom *geom, FTYPE *v2);
extern void gaussj(FTYPE **tmp, int n, FTYPE **b, int m);
extern void set_points(void);
// extern FTYPE delta(int j, int k) ;
// extern FTYPE mink(int j, int k) ;
extern void make_tetr(FTYPE *ucon, FTYPE econ[][NDIM]);
extern void dxdxprim(FTYPE *X, FTYPE r, FTYPE th, FTYPE *dxdxp);
