HDF_INSTALL = /usr/local/hdf5
EXTLIB = -L$(HDF_INSTALL)/lib  -L/opt/local/lib/ -L/opt/local/lib/openmpi-gcc11/ -L/usr/local/lib/
CC          = gcc-mp-11 #/opt/local/libexec/openmpi-gcc11/mpiCC 
#CFLAGS      = -Wall -O3 -fopenmp
CFLAGS      = -Wall -g  -fsanitize=address
# -fopenmp 
LIB         = -lz -lm -ldl -lgsl -lgslcblas -lm -lmpi 

DEPS = mcrat.h mclib.h mclib_riken.h mclib_pluto.h mc_cyclosynch.h mcrat_input.h geometry.h mcrat_io.h mcrat_scattering.h mclib_flash.h analytic_outflows.h optical_depth.h electron.h hot_x_section.h photons.h custom_outflow.h custom_spectrum.h
OBJ = mcrat.o mclib.o mclib_riken.o mclib_pluto.o mc_cyclosynch.o geometry.o mcrat_io.o mcrat_scattering.o mclib_flash.o analytic_outflows.o optical_depth.o electron.o hot_x_section.o photons.o custom_outflow.o custom_spectrum.o
OBJ_MERGE = merge.o mclib.o mclib_riken.o mclib_pluto.o mc_cyclosynch.o geometry.o mcrat_io.o mcrat_scattering.o mclib_flash.o analytic_outflows.o optical_depth.o electron.o hot_x_section.o photons.o custom_outflow.o custom_spectrum.o

INCLUDE   = -I$(HDF_INSTALL)/include -I/opt/local/include/ -I/usr/include/ -I/opt/local/include/openmpi-gcc11/

LIBSHDF   = $(EXTLIB) $(HDF_INSTALL)/lib/libhdf5.a

# ── Files that must not be compiled with -ffast-math ─────────────────────────
#
# mc_synchrotron.c: uses GSL adaptive quadrature (gsl_integration_qagiu /
#   qags). The subdivision convergence loop tests for inf/NaN internally.
#   -ffinite-math-only (part of -ffast-math) causes an infinite loop.
#
# photons.c: uses NAN as a sentinel value for unset photon fields, and
#   tests them with isnan(). -ffinite-math-only folds isnan() to 0
#   (always false), silently disabling all field-validation guards in
#   saveUserDefinePhoton / initializePhoton.
#
SAFE_CFLAGS = $(filter-out -ffast-math, $(CFLAGS)) -fno-finite-math-only


MCRAT: $(OBJ)
	$(CC) $(CFLAGS)  -o $@ $^ $(INCLUDE) $(LIBSHDF) $(LIB)

mc_synchrotron.o: mc_synchrotron.c
       $(CC) $(SAFE_CFLAGS) -c $< -o $@ $(INCLUDE) $(LIBSHDF) $(LIB)

photons.o: photons.c
       $(CC) $(SAFE_CFLAGS) -c $< -o $@ $(INCLUDE) $(LIBSHDF) $(LIB)

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS)  -c -o $@ $< $(INCLUDE) $(LIBSHDF) $(LIB)


MERGE: $(OBJ_MERGE) 
	$(CC) $(CFLAGS)  -o $@ $^ $(INCLUDE) $(LIBSHDF) $(LIB)

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS)  -c -o $@ $< $(INCLUDE) $(LIBSHDF) $(LIB)


clean: 
	rm -f *.o 
 

.SUFFIXES:.o.c
