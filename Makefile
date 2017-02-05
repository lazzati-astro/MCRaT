HDF_INSTALL = /usr/local/hdf5
EXTLIB = -L$(HDF_INSTALL)/lib -L/usr/local/lib
CC          = gcc
CFLAGS      = -Wall -O2 
LIB         = -lz -lm -ldl -lgsl -lgslcblas -lm

DEPS = mclib.h
OBJ = mcrat.o mclib.o

INCLUDE   = -I$(HDF_INSTALL)/include -I/usr/local/include/
LIBSHDF   = $(EXTLIB) $(HDF_INSTALL)/lib/libhdf5.a 

#all: h5_crtdat \
#     h5_rdwt \
#     h5_crtatt \
#     h5_crtgrp \
#     h5_crtgrpar \
#     h5_crtgrpd \
 

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $< $(INCLUDE) $(LIBSHDF) $(LIB)

MCRAT: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(INCLUDE) $(LIBSHDF) $(LIB)


clean: 
	rm -f *.o 
 

.SUFFIXES:.o.c
