HDF_INSTALL = /usr/local/hdf5
EXTLIB = -L$(HDF_INSTALL)/lib  -L/opt/local/lib/ -L/opt/local/lib/openmpi-gcc11/ -L/usr/local/lib/
CC          = gcc-mp-11 #/opt/local/libexec/openmpi-gcc11/mpiCC 
#CFLAGS      = -Wall -O3 -fopenmp
CFLAGS      = -Wall -g  -fsanitize=address
# -fopenmp 
LIB         = -lz -lm -ldl -lgsl -lgslcblas -lm -lmpi 

DEPS = mclib.h mclib_3d.h mclib_pluto.h mc_cyclosynch.h mcrat_input.h
OBJ = mcrat.o mclib.o mclib_3d.o mclib_pluto.o mc_cyclosynch.o

INCLUDE   = -I$(HDF_INSTALL)/include -I/opt/local/include/ -I/usr/include/ -I/opt/local/include/openmpi-gcc11/

LIBSHDF   = $(EXTLIB) $(HDF_INSTALL)/lib/libhdf5.a 

MCRAT: $(OBJ)
	$(CC) $(CFLAGS)  -o $@ $^ $(INCLUDE) $(LIBSHDF) $(LIB)

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS)  -c -o $@ $< $(INCLUDE) $(LIBSHDF) $(LIB)


MERGE: merge.o mclib.o 
	$(CC) $(CFLAGS)  -o $@ $^ $(INCLUDE) $(LIBSHDF) $(LIB)

merge.o: merge.c mclib.c mclib.h
	$(CC) $(CFLAGS)  -c -o $@ $< $(INCLUDE) $(LIBSHDF) $(LIB)


clean: 
	rm -f *.o 
 

.SUFFIXES:.o.c
