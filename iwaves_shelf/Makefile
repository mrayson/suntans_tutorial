#SUNTANSHOME=../../main
SUNTANSHOME=/home/mrayson/software/code/suntans/main
include $(SUNTANSHOME)/Makefile.in

ifneq ($(MPIHOME),)
  CC = $(MPIHOME)/bin/mpicc
  MPIDEF = 
  MPIINC = -I$(MPIHOME)
else
  CC = cc 
  #MPIDEF = -DNOMPI
  MPIDEF = 
  MPIINC = 
endif

ifneq ($(PARMETISHOME),)
  PARMETISINC = -I$(PARMETISHOME)/ParMETISLib
endif

LD = $(CC) 
CFLAGS = -O2 -fcommon
MATHLIB = -lm

EXEC = iwaves_shelf
OBJS = 
SUN = $(SUNTANSHOME)/sun
INCLUDES = -I$(SUNTANSHOME) $(MPIINC) $(PARMETISINC)
DEFS = $(MPIDEF)
NUMPROCS = 4
datadir = data

all:	data

test:	data
	sh $(EXEC).sh $(NUMPROCS) 

restart:    data
	sh $(EXEC)-restart.sh $(NUMPROCS) 1

data:	state.o phys.o turbulence.o
	cp state.o phys.o turbulence.o $(SUNTANSHOME)
	make -C $(SUNTANSHOME)

.c.o:	
	$(LD) $(CFLAGS) -c $(INCLUDES) $(DEFS) $*.c

$(SUN):	state.o phys.o turbulence.o
	cp state.o phys.o turbulence.o $(SUNTANSHOME)
	make -C $(SUNTANSHOME)

debug:	data
	$(MPIHOME)/bin/mpirun -np $(NUMPROCS) xterm -e gdb -command=gdbcommands.txt $(SUN)

valgrind: data
	mkdir $(datadir)
	cp rundata/* $(datadir)
	$(MPIHOME)/bin/mpirun -np $(NUMPROCS) ./$(SUN) -g -vv --datadir=$(datadir)
	$(MPIHOME)/bin/mpirun -np $(NUMPROCS) valgrind --tool=memcheck --leak-check=yes ./$(SUN) -s -vvv --datadir=$(datadir)

clean:
	#rm -f $(SUNTANSHOME)/*.o
	rm -f *.o

clobber: clean
	rm -rf *~ \#*\# PI* $(EXEC) gmon.out data rundata/*~

