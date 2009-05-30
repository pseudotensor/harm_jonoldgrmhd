#include "decs.h"

#define DEBUG 0

void bound_mpi(FTYPE prim[][N2 + 4][NPR])
{
  int dir;
#if(DEBUG)
  int i,j,k;
#endif

#if(USEMPI)
  /* These arrays contain designations that identify 
   * each recv and send */
  static MPI_Request requests[COMPDIM * 2 * 2];
  // format of map for requests[dir*2+recv/send(0/1)]
#endif

#if(USEMPI)

  /*
   *
   * 1. do outflow/inflow/etc. boundary conditions (in bounds)
   * 2. pack data into workbc arrays and do send/recv's
   * 3. for each transfer do a wait/unpack
   * 4. do all send waits
   *
   * NOTE: order is important so corner zones are assigned correctly
   *
   * workbc[PACK][][] is data that is being sent
   * workbc[UNPACK][][] is data that is being recvd
   *
   */

  /*
   *  -> higher values of x1
   * 2 -> lower values of x1
   *
   */

  
  /* bounds has already done non-MPI boundary conditions; now do MPI
   * boundary conditions */

  // must go in this order (0,1) then (2,3) or visa versa, a single
  // directions order doesn't matter.  This method of ordering is as
  // opposed to directly transfering the corner zones to the corner
  // CPU.

  // This may be faster since all transfers can proceed at once,
  // although this may be slower since no transfers can occur until
  // packing is completed.  This way packing and transfering occur
  // simultaneously.

  // Although l/r are packed together, since in the end we have to
  // wait for both l/r to complete, so equal time completion is
  // favored //over asynch completion.

  // Also, transfering corner zones with small message sizes increases
  // the importance of latency.

  // I choose left-right N1M/N2M first, then up/down N1M/N2M.  Could
  // just do N1/N2 for interior for L/R, but true boundary needs full
  // N1M/N2M exchanged since cpu sets boundary using normal bc code
  // which needs to get transfered to neight(i.e. currently if corner
  // has bctype 99/?  then doesn't do corner)

#if(DEBUG)
  PLOOP{
    for(i=-2;i<N1+2;i++) for(j=-2;j<N2+2;j++){
      prim[i][j][k]=-1-k*100; // should be turned into myid-k*100
    }
    for(i=0;i<N1;i++) for(j=0;j<N2;j++){
      prim[i][j][k]=myid-k*100; // differentiates but clear per pr
      //      fprintf(log_file,"%d %d %d %15.10g\n",i,j,k,prim[i][j][k]);
    }
  }
#endif




  for(dir=X1UP;dir<=X1DN;dir++) if(dirset[dir][DIRIF]) pack(dir,prim,workbc);
  for(dir=X1UP;dir<=X1DN;dir++) if(dirset[dir][DIRIF]) sendrecv(dir,workbc,requests);
  for(dir=X1UP;dir<=X1DN;dir++) if(dirset[dir][DIRIF]) recvwait(dir,requests);
  for(dir=X1UP;dir<=X1DN;dir++) if(dirset[dir][DIRIF]) unpack(dir,workbc,prim);
  for(dir=X1UP;dir<=X1DN;dir++) if(dirset[dir][DIRIF]) sendwait(dir,requests);

  // now dir=0,1(X1UP,X1DN) is done, so can start 2,3(X2UP,X2DN)
  for(dir=X2UP;dir<=X2DN;dir++) if(dirset[dir][DIRIF]) pack(dir,prim,workbc);
  for(dir=X2UP;dir<=X2DN;dir++) if(dirset[dir][DIRIF]) sendrecv(dir,workbc,requests);
  for(dir=X2UP;dir<=X2DN;dir++) if(dirset[dir][DIRIF]) recvwait(dir,requests);
  for(dir=X2UP;dir<=X2DN;dir++) if(dirset[dir][DIRIF]) unpack(dir,workbc,prim);
  for(dir=X2UP;dir<=X2DN;dir++) if(dirset[dir][DIRIF]) sendwait(dir,requests);

#if(DEBUG)
    fprintf(log_file,"\n\nafter\n\n");
    for(i=-2;i<N1+2;i++) for(j=-2;j<N2+2;j++){
      PLOOP{
        fprintf(log_file,"%d %d %d %15.10g\n",i,j,k,prim[i][j][k]);
      }
    }
    myexit(0);
#endif



  // end if mpi
#endif

}	
// end function


// PACKLOOP allows one to alter which i or j is faster iterated
#define PACKLOOP(i,j,istart,istop,jstart,jstop) GENLOOP(i,j,istart,istop,jstart,jstop) PLOOPMPI

// packs data for shipment
void pack(int dir,FTYPE prim[][N2 + 4][NPR],FTYPE workbc[][COMPDIM * 2][NPR * NBIGBND * NBIGSM])
{
  // dir=direction sending
  int i,j;
  int bci,pr;
  int start1,start2,stop1,stop2;

  bci=0;
  PACKLOOP(i,j
	   ,dirset[dir][DIRPSTART1]
	   ,dirset[dir][DIRPSTOP1]
	   ,dirset[dir][DIRPSTART2]
	   ,dirset[dir][DIRPSTOP2]){
    /*
    if(bci>=dirset[dir][DIRSIZE]){
      dualfprintf(fail_file,"pack memory leak: bci: %d dirset[%d][DIRSIZE]: %d\n",bci,dirset[dir][DIRSIZE]);
      myexit(10);
    }
    */
    workbc[PACK][dir][bci++] = prim[i][j][pr];
  }
}

#if(USEMPI)
void sendrecv(int dir,FTYPE workbc[][COMPDIM * 2][NPR * NBIGBND * NBIGSM],MPI_Request *requests)
{
  MPI_Irecv(workbc[UNPACK][dir],
	    dirset[dir][DIRSIZE],
	    MPI_FTYPE,
	    dirset[dir][DIROTHER],
	    dirset[dir][DIRTAGR],
	    MPI_COMM_WORLD,
	    &requests[dir*2+REQRECV]);

  MPI_Isend(workbc[PACK][dir],
	    dirset[dir][DIRSIZE],
	    MPI_FTYPE,
	    dirset[dir][DIROTHER],
	    dirset[dir][DIRTAGS],
	    MPI_COMM_WORLD,
	    &requests[dir*2+REQSEND]);

}
#endif


#if(USEMPI)
void recvwait(int dir,MPI_Request *requests)
{
  MPI_Wait(&requests[dir*2+REQRECV], &mpichstatus);

}
#endif


void unpack(int dir,FTYPE workbc[][COMPDIM * 2][NPR * NBIGBND * NBIGSM],FTYPE prim[][N2 + 4][NPR])
{
  // dir is direction receiving from
  int i,j;
  int bci,pr;

  bci=0;
  PACKLOOP(i,j
	   ,dirset[dir][DIRUSTART1]
	   ,dirset[dir][DIRUSTOP1]
	   ,dirset[dir][DIRUSTART2]
	   ,dirset[dir][DIRUSTOP2]){
    prim[i][j][pr]=workbc[UNPACK][dir][bci++];
  }
}

#if(USEMPI)

void sendwait(int dir,MPI_Request *requests)
{
  MPI_Wait(&requests[dir*2+REQSEND], &mpichstatus);
}
#endif
