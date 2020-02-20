HDF_INSTALL = /usr/local/hdf5
EXTLIB = -L$(HDF_INSTALL)/lib -L/usr/local/lib/ -L/opt/local/lib/
CC          = gcc-mp-9
#CFLAGS      = -Wall -O2 -fopenmp
CFLAGS      = -Wall -g -fopenmp -fsanitize=address
LIB         = -lz -lm -ldl -lgsl -lgslcblas -lm -lmpi 

DEPS = mclib.h mclib_3d.h mclib_pluto.h mc_synch.h mcrat_input.h
OBJ = mcrat.o mclib.o mclib_3d.o mclib_pluto.o mc_synch.o

INCLUDE   = -I$(HDF_INSTALL)/include -I/usr/local/include/ -I/usr/include/ #-I/opt/local/lib/gcc8/gcc/x86_64-apple-darwin18/8.2.0/include-fixed/

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
