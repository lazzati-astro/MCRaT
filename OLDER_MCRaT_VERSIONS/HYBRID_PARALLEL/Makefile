HDF_INSTALL = /usr/local/hdf5
EXTLIB = -L$(HDF_INSTALL)/lib -L/usr/local/lib -L/opt/local/lib/
CC          = gcc-mp-8
CFLAGS      = -Wall -g -fopenmp
LIB         = -lz -lm -ldl -lgsl -lgslcblas -lm -lmpi

DEPS = mclib.h mclib_3d.h
OBJ = mcrat.o mclib.o mclib_3d.o

INCLUDE   = -I$(HDF_INSTALL)/include -I/usr/local/include/ -I/opt/local/lib/gcc5/gcc/x86_64-apple-darwin15/5.3.0/include/
LIBSHDF   = $(EXTLIB) $(HDF_INSTALL)/lib/libhdf5.a 

#all: h5_crtdat \
#     h5_rdwt \
#     h5_crtatt \
#     h5_crtgrp \
#     h5_crtgrpar \
#     h5_crtgrpd \
 

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS)  -c -o $@ $< $(INCLUDE) $(LIBSHDF) $(LIB)

MCRAT: $(OBJ)
	$(CC) $(CFLAGS)  -o $@ $^ $(INCLUDE) $(LIBSHDF) $(LIB)


clean: 
	rm -f *.o 
 

.SUFFIXES:.o.c
