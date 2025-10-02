FC		= mpiifort
CC		= mpiicc

#FFLAGS		= -fast -fpp -Dssurf=2
#CFLAGS		= -fast
#LFLAGS		= -fast -fpp -Dssurf=2

FFLAGS		= -ipo -O3 -no-prec-div -fp-model fast=2 -xHost -fpp -Dssurf=2
CFLAGS		= -ipo -O3 -no-prec-div -fp-model fast=2 -xHost
LFLAGS		= -ipo -O3 -no-prec-div -fp-model fast=2 -xHost -fpp -Dssurf=2

#FFLAGS		= -fpp -Dwall=1 -CB -check all,noarg_temp_created -traceback -fpe0 -warn -C -no-ip -Dssurf=4
#CFLAGS		= -traceback -fpe0 -warn -C -no-ip
#LFLAGS		= -fpp -Dwall=1 -CB -check all,noarg_temp_created -traceback -fpe0 -warn -C -no-ip -Dssurf=4

#LDFLAGS		= 
#HDF5DIR		= 
#FFTWDIR		= 

LDFLAGS		= -L/LARGE0/gr20001/KUspace-share/common/hdf5-lib/hdf5-1.14.3-240410/lib -lhdf5_hl -lhdf5_fortran -lhdf5 \
		  -L/LARGE0/gr20001/KUspace-share/common/fftw-lib/fftw-3.3.10-240410/lib -lfftw3_mpi -lfftw3 \
		  -L../../ctca/lib -lctca_mod -lctca_f -lctca
HDF5DIR		= /LARGE0/gr20001/KUspace-share/common/hdf5-lib/hdf5-1.14.3-240410/include
FFTWDIR		= /LARGE0/gr20001/KUspace-share/common/fftw-lib/fftw-3.3.10-240410/include
CTCADIR		= ../../ctca/include

#LDFLAGS		= -L/opt/system/app/hdf5/1.12.2/intel-2022.3-impi-2022.3/lib -lhdf5_hl -lhdf5_fortran -lhdf5 \
#		  -L/opt/system/app/fftw/3.3.10/d_intel-2022.3-impi-2022.3/lib -lfftw3_mpi -lfftw3
#HDF5DIR		= /opt/system/app/hdf5/1.12.2/intel-2022.3-impi-2022.3/include
#FFTWDIR		= /opt/system/app/fftw/3.3.10/d_intel-2022.3-impi-2022.3/include

LINKER		= $(FC)

COMMONHDRS	= oh_config.h oh_stats.h
OHHDRS		= ohhelp1.h ohhelp2.h ohhelp3.h oh_part.h
OHOBJS		= ohhelp1.o ohhelp2.o ohhelp3.o

FHDRS		= ohhelp_f.h
FMODS		= oh_type.o oh_mod1.o oh_mod2.o oh_mod3.o

KMODS		= paramt.o interf.o allcom.o namels.o hdf.o ssl2ex.o hdf5io.o wrapper.o lemses.o bem.o

FOBJS		= main.o emses.o esses.o imses.o \
		  input.o ohinit.o inital.o inipcl.o inipin.o \
		  inifld.o getccl.o capmtx.o chgntr.o \
		  rescal.o chkprm.o \
		  digset.o energy.o hdfdig.o medium.o \
		  psolve1.o psolve2.o positn.o psuper1.o psuper2.o psuper3.o \
		  charge.o rhoupd.o vsolves.o \
		  iemfld.o defmpi.o inifft.o poissn.o \
		  efield.o bfield.o esfild.o surchg.o ptnvar.o \
		  grdcrt.o bdyfld.o extfld.o emfbnd.o fsmask.o \
		  frmttd.o spsort.o snap.o nborps.o chkres.o

COBJS		= mt19937arm.o psort-struct.o

OBJS		= $(FOBJS) $(COBJS) $(FMODS) $(KMODS) $(LMODS) $(OHOBJS)

mpiemses3D:	$(OBJS)
		$(LINKER) $(OBJS) $(LFLAGS) $(LDFLAGS) -o $@

$(FOBJS):%.o:	%.F90 $(FMODS) $(KMODS) $(LMODS) $(FHDRS) $(COMMONHDRS)
		$(FC) $(FFLAGS) -I$(FFTWDIR) -I$(HDF5DIR) -I$(CTCADIR) -c $< -o $@
$(COBJS):%.o:	%.c $(OHHDRS)
		$(CC) $(CFLAGS) -c $< -o $@
%simulator.o:	sample.o

$(FMODS):%.o:	%.F90 $(COMMONHDRS)
		$(FC) $(FFLAGS) -I$(FFTWDIR) -I$(HDF5DIR) -I$(CTCADIR) -c $< -o $@
$(KMODS):%.o:	%.F90
		$(FC) $(FFLAGS) -I$(FFTWDIR) -I$(HDF5DIR) -I$(CTCADIR) -c $< -o $@
oh_mod1.o:	oh_type.o
oh_mod2.o:	oh_mod1.o
oh_mod3.o:	oh_mod2.o
allcom.o:	oh_type.o paramt.o lemses.o
paramt.o:	interf.o
lemses.o:	hdf5io.o wrapper.o
wrapper.o:	ssl2ex.o
bdyfld.o:	vical.fnc
psolve1.o:	defsurf.fnc defbody.fnc
psolve2.o:	defsurf.fnc defbody.fnc
psuper1.o:	setpos.fnc
psuper2.o:	setpos.fnc
inipcl.o:	defsurf2.fnc
$(OHOBJS):%.o:	%.c $(COMMONHDRS) $(OHHDRS)
		$(CC) $(CFLAGS) -c $< -o $@

clean:;
		rm *.o; rm *.mod; rm mpiemses3D
clear:;
		mv hdf5io.o hdf5io.tmp; rm *.o; rm *.mod; mv hdf5io.tmp hdf5io.o
