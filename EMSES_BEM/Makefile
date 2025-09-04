PREFIX=../ctca
#PREFIX=/home/tmp/cotocoa/cotocoa-1.2.3-intel2023.2
CC=mpiicc
FC=mpiifort
FFLAGS=-I${PREFIX}/include
LDFLAGS=-L${PREFIX}/lib
LIBS=-lctca
F_LIBS=-lctca_f -lctca_mod -qmkl

COUPDIR=coupler
WORKDIR=worker
BINDIR=bin

PROGS = $(COUPDIR)/coupler $(WORKDIR)/worker

all: $(PROGS)

# requester: $(REQDIR)/requester.f90
# 	$(FC) $(FFLAGS) $(LDFLAGS) requester.f90 $(F_LIBS) $(LIBS) -o requester

$(COUPDIR)/coupler: $(COUPDIR)/coupler.f90
	$(FC) $(FFLAGS) $(LDFLAGS) $(COUPDIR)/coupler.f90 $(F_LIBS) $(LIBS) -o $(COUPDIR)/coupler

$(WORKDIR)/worker: $(WORKDIR)/worker.f90
	$(FC) $(FFLAGS) $(LDFLAGS) $(WORKDIR)/worker.f90 $(F_LIBS) $(LIBS) -o $(WORKDIR)/worker

clean:
	/bin/rm -f *.o *.mod ${PROGS}
