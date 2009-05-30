#
#
#echo "Make sure MPICH/GM params set in both global.h and makefile!"
#echo "Make sure MPICH/GM params set in both global.h and makefile!"
USEMPI=1

USEICC=1
USEGCC=0
USECCC=0

ifeq ($(USEMPI),1)
 MCC = mpicc
#MCC=/usr/local/p4mpich-1.2.5-icc-noshmem/bin/mpicc
endif    

LAPACKLIB=lapack
BLASLIB=blas
F2CLIB=g2c  #-l$(F2CLIB)

PREP = prep
CMD=grmhd


#
# Define preprocessor and compile flags, and linked libraries

ifeq ($(USEGCC),1)
COMP=gcc
# COMP=gcc3

# gcc type flags
#CFLAGS = -Wall -mpentium -O3 -pipe  -malign-loops=2 -malign-jumps=2 -malign-functions=2 -DCPU=686 -DNEED_GETOPT -DLINUX -ffast-math
#CFLAGS = -Wall -mpentium -O3 -pipe  -malign-loops=2 -malign-jumps=2 -malign-functions=2 -DCPU=686 -DNEED_GETOPT -DLINUX -ffast-math -pg
#CFLAGS = -Wall -mcpu=pentiumpro -march=pentiumpro -O4 -pipe  -malign-loops=2 -malign-jumps=2 -malign-functions=2 -malign-double -mstack-align-double -ffast-math -pg
#CFLAGS = -Wall -mcpu=pentiumpro -march=pentiumpro -O4 -malign-loops=2 -malign-jumps=2 -malign-functions=2 -malign-double -mstack-align-double -ffast-math -finline-functions -pg
#CFLAGS = -Wall -mcpu=pentiumpro -march=pentiumpro -O4 -malign-loops=2 -malign-jumps=2 -malign-functions=2 -malign-double -ffast-math -finline-functions -pg
#CFLAGS = -Wall -mcpu=pentiumpro -march=pentiumpro -O3 -malign-loops=2 -malign-jumps=2 -malign-functions=2 -ffast-math -finline-functions -pg -g -a
#CFLAGS = -Wall -mcpu=pentiumpro -march=pentiumpro -O3 -ffast-math -finline-functions -funroll-loops
CFLAGS=-O0 -g
#CFLAGS = -Wall -mcpu=pentiumpro -march=pentiumpro -O4 -malign-loops=2 -malign-jumps=2 -malign-functions=2 -ffast-math -finline-functions -g
#-pg
#-pg -g  source lines
#-pg -g -a   line count
# gprof -l <file> > out.txt
# gprof -A -I<sourcedir>
# gprof -l -A -x s

#below is typical flags for double precision...can take -pg off for no profile
#add -mstack-align-double if using pgcc
#CFLAGS = -Wall -O0 -g
#  -fomit-frame-pointer



#CFLAGS = -Wall -O0
#CFLAGS = -O6 -g
#CFLAGS = -O0 -pg -g
LDFLAGS = -lm
# -l$(LAPACKLIB) -l$(BLASLIB)  -L/usr/lib/gcc-lib/i386-redhat-linux/2.96/ -l$(F2CLIB) 

#CC = cc
#AR	=	ar r
#RANLIB	=	ranlib
endif

ifeq ($(USEICC),1)
COMP=icc
# P4 (don't use -wp_ipo -- incompat code)
#CFLAGS=-O3 -tpp7 -axiMKW -wp_ipo -unroll -w1
CFLAGS=-O3 -tpp7 -axiMKW -ipo -unroll -w1 -wd=175
#CFLAGS=-O0 -g
# P3
# CFLAGS=-O3 -tpp6 -axiMK -ipo -unroll -w1
# GAMMIE
# CFLAGS = -O3 -ipo
#CFLAGS = -O0

LDFLAGS = -lm
endif


ifeq ($(USECCC),1)
COMP=ccc
LDFLAGS =  -lm -lcxml

#CDEBUG = -g3 # -g3 for higher opts than -O0
#CDEBUG = -g
# do profile (profile screws up speed for loops, doesn't unroll them, etc.)
#CDEBUG = -pg -g3
# production level
CFLAGS3 = -Wall -O4 -fast -msg_disable badsubscript -msg_disable subscrbounds -msg_disable unreachcode -msg_disable noparmlist -msg_disable subscrbounds2 -msg_disable longlongtype -finline-functions -funroll-loops
#CFLAGS3 = -Wall -O0
# super annoying develop level
#CFLAGS3 = -Wall -O2 -fast
#CFLAGS3 = -fast -arch ev67
# debug level
#CFLAGS3 = -Wall -O0 -msg_disable badsubscript -msg_disable subscrbounds -msg_disable unreachcode -msg_disable noparmlist -msg_disable subscrbounds2
#CFLAGS3 = -Wall -O0
#CFLAGS3 = -Wall -O2

endif



ifeq ($(USEMPI),1)
CC=$(MCC) -cc=$(COMP) $(CFLAGS3) $(CFLAGS2) $(CDEBUG)
endif # endif usempich==1

ifeq ($(USEMPI),0)
CC=$(COMP)  $(CFLAGS3) $(CFLAGS2) $(CDEBUG)
endif # endif usempich==0



SRCS = \
bounds.c coord.c diag.c dudp_calc.c dump.c fixup.c gaussj.c \
image.c init.c interp.c lubksb.c ludcmp.c main.c metric.c \
mnewt.c nrutil.c phys.c ranc.c restart.c set_arrays.c set_grid.c \
step_ch.c tensor.c utoprim.c vchar.c transforms.c fail.c utoprim_ldz.c \
init_mpi.c boundmpi.c prepostinit.c cool.c
 
OBJS = \
bounds.o coord.o diag.o dudp_calc.o dump.o fixup.o gaussj.o \
image.o init.o interp.o lubksb.o ludcmp.o main.o metric.o \
mnewt.o nrutil.o phys.o ranc.o restart.o set_arrays.o set_grid.o \
step_ch.o tensor.o utoprim.o vchar.o transforms.o fail.o utoprim_ldz.o \
init_mpi.o boundmpi.o prepostinit.o cool.o

all:	$(PREP) $(CMD)


$(PREP):
	( sh ./makedecs.h.sh )
	( sh ./makempidecs.h.sh )

$(CMD):	$(OBJS) makefile
	$(CC) $(CFLAGS) -o $(CMD) $(OBJS) $(LDFLAGS)


clean:
	rm *.o

# dependencies
$(OBJS) : global.h defs.h mympi.h mpidefs.h makefile
metric.o : metric.h
transforms.o : metric.h


OBJB = postmort.o \
bounds.o coord.o diag.o dudp_calc.o dump.o fixup.o gaussj.o \
image.o init.o interp.o lubksb.o ludcmp.o metric.o \
mnewt.o nrutil.o phys.o ranc.o restart.o set_arrays.o set_grid.o \
step_ch.o tensor.o utoprim.o vchar.o transforms.o \
init_mpi.o boundmpi.o prepostinit.o cool.o


# dependencies
$(OBJB) : global.h defs.h mympi.h mpidefs.h makefile


postmort: $(OBJB) makefile
	$(CC) $(CFLAGS) -o postmort $(OBJB) $(LDFLAGS)
