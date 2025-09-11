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
REQDIR=requester
BINDIR=bin

REQPROG=$(REQDIR)/mpiemses3D

PROGS = $(COUPDIR)/coupler $(WORKDIR)/worker $(REQPROG)

.PHONY: all clean requester

all: $(PROGS)

$(COUPDIR)/coupler: $(COUPDIR)/coupler.f90
	$(FC) $(FFLAGS) $(LDFLAGS) $(COUPDIR)/coupler.f90 $(F_LIBS) $(LIBS) -o $(COUPDIR)/coupler

$(WORKDIR)/worker: $(WORKDIR)/worker.f90
	$(FC) $(FFLAGS) $(LDFLAGS) $(WORKDIR)/worker.f90 $(F_LIBS) $(LIBS) -o $(WORKDIR)/worker

$(REQPROG):
	$(MAKE) -C $(REQDIR)

requester: $(REQPROG)

clean:
	/bin/rm -f *.o *.mod $(COUPDIR)/coupler $(WORKDIR)/worker
	$(MAKE) -C $(REQDIR) clean || true