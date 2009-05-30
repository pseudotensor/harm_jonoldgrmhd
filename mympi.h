/* 
   To start mpi run, first do:

   mpirun -np 4 ./grmhd

   e.g. 4 cpus using mpi:

   rm nohup.out ; nohup sh -c 'mpirun -np 4 ./grmhd' &

   Note: cannot have any cpu write to same file or pipe to same file
   without parallel I/O or blocking through each CPU

 */

#define USEMPI 1

// whether to use MPI over GM
// This forces no forks to be called.  Can be used for non-gm too
#if(USEMPI==1)
#include "mpi.h"
#define USEGM 0			// choice (avoids system/fork/etc calls)
#else
#define USEGM 0			// always 0, can't have GM without MPI
#endif

// need not change below datatype stuff
#if(REALTYPE==FLOATTYPE)
#define MPI_FTYPE MPI_FLOAT
#else
#define MPI_FTYPE MPI_DOUBLE
#endif

#if(SENSITIVE==FLOATTYPE) // for sensitive counters
#define MPI_SFTYPE MPI_FLOAT
#else
#define MPI_SFTYPE MPI_DOUBLE
#endif


#define COMPDIM 2


#define N1BND 2
#define N2BND 2

#define NBIGBND ((N1BND>N2BND) ? N1BND : N2BND)

#define N1M (N1+2*N1BND)
#define N2M (N2+2*N2BND)

/* NBIGM is bigger of N1M and N2M and N3M */
#define NBIGM ((N1M>N2M) ? N1M : N2M)

#define MAXCPUS 1000
#define MAXFILENAME 200

// maximal surface of boundary exchange
#define NBIGSM (NBIGM)

#define BUFFERMAP ((j*N1+i)*numcolumns+nextbuf++)
#define BUFFERINIT if(mpicombine==1) nextbuf=0
#define IMAGEMAP (j*N1+i)


#define X1UP	0
#define X1DN	1
#define X2UP	2
#define X2DN	3

#define PACK 1
#define UNPACK 2

#define REQRECV 0
#define REQSEND 1

#define DIRNUMVARS 14
#define DIRIF      0
#define DIRSIZE    1
#define DIROTHER   2
#define DIRTAGS    3
#define DIRTAGR    4
#define DIROPP     5
#define DIRPSTART1 6
#define DIRPSTOP1  7
#define DIRUSTART1 8
#define DIRUSTOP1  9
#define DIRPSTART2 10
#define DIRPSTOP2  11
#define DIRUSTART2 12
#define DIRUSTOP2  13

#define TEXTOUTPUT 0
#define BINARYOUTPUT 1
#define UNSORTED 0
#define SORTED 1


#define STAGE1 1
#define STAGE2 2

// #define DATADIR "./"
#define DATADIR ""

// extention for data files
#define DATEXT ".dat"
#define PAREXT ".par"
#define INEXT ".in"
#define OUTEXT ".out"
#define PPEXT ".pp"

#define CPUTXT ".%04d"

#define PLOOPMPI for(pr=0;pr<NPR;pr++)


extern void init_mpi(int argc, char *argv[]);
extern void init_genfiles(int gopp);
extern int init_MPI(int argc, char *argv[]);
extern void init_placeongrid(void);
extern int myexit(int call_code);

extern void bound_mpi(FTYPE prim[][N2 + 4][NPR]);
extern void pack(int dir,FTYPE prim[][N2 + 4][NPR],FTYPE workbc[][COMPDIM * 2][NPR * NBIGBND * NBIGSM]);
extern void unpack(int dir,FTYPE workbc[][COMPDIM * 2][NPR * NBIGBND * NBIGSM],FTYPE prim[][N2 + 4][NPR]);
#if(USEMPI)
extern void sendrecv(int dir,FTYPE workbc[][COMPDIM * 2][NPR * NBIGBND * NBIGSM],MPI_Request *requests);
extern void recvwait(int dir,MPI_Request *request);
extern void sendwait(int dir,MPI_Request *request);
#endif

extern void mpimax(SFTYPE*maxptr);
extern void mpimin(SFTYPE*minptr);
extern void mpifmin(FTYPE*minptr);

#if(USEMPI)
extern void mpiio_init(int bintxt, int sorted, FILE ** fp, int which,
		       char *filename, int numcolumns,
		       int datatype, void **jonio, void **writebuf8);
extern void mpiio_combine(int bintxt, int sorted,
			  int numcolumns, int datatype,
			  FILE ** fp, void *jonio, void *writebuf);
extern void mpiio_seperate(int bintxt, int sorted, int stage,
			   int numcolumns,
			   int datatype, FILE ** fp, void *jonio,
			   void *writebuf);
extern void mpiios_init(int bintxt, int sorted, FILE ** fp, int which, char *filename, int numcolumns,
			int datatype, void **jonio, void **writebuf);
extern void mpiios_combine(int bintxt, MPI_Datatype mpidt, int numcolumns,
			   int datatype,
			   FILE ** fp, void *jonio, void *writebuf);
extern void mpiios_seperate(int bintxt, int stage, MPI_Datatype mpidt, int numcolumns,
			    int datatype, FILE ** fp, void *jonio,
			    void *writebuf);
extern void mpiiotu_combine(MPI_Datatype mpidt, int numcolumns, int datatype,
			    FILE ** fp, void *writebuf);
extern void mpiiotu_seperate(int stage, MPI_Datatype mpidt, int numcolumns,
			     int datatype, FILE ** fp,void *writebuf);
#endif


#define MYOUT stderr		// normally stderr
