/************************************************************

  This example shows how to read and write data to a
  dataset.  The program first writes integers to a dataset
  with dataspace dimensions of DIM0xDIM1, then closes the
  file.  Next, it reopens the file, reads back the data, and
  outputs it to the screen.

  This file is intended for use with HDF5 Library version 1.8

 ************************************************************/

#include "hdf5.h"
#include <stdio.h>
#include <stdlib.h>

#define FILE            "/Users/Tylerparsotan/Documents/PYTHON/MCRAT/16OI/rhd_jet_big_16OI_hdf5_plt_cnt_0000"
#define DATASET         "dens"
#define DIM0            2852
#define DIM1            1
#define DIM2            8
#define DIM3            8


int
main (void)
{
    hid_t       file, space, dset, memtype;          /* Handles */
    herr_t      status;
    hsize_t     dims[4] = {DIM0, DIM1,DIM2,DIM3};
    double         wdata[DIM0][DIM1][DIM2][DIM3],          /* Write buffer */
                rdata[DIM0][DIM1][DIM2][DIM3];          /* Read buffer */
                
    hsize_t         i, j,k,l;
    int ndims;

    /*
     * Initialize data.
     */
    for (i=0; i<DIM0; i++){
        for (j=0; j<DIM1; j++){
            for (k=0; k<DIM2; k++){
                for (l=0; l<DIM1; l++){
                    wdata[i][j][k][l] = 1*i;
                }
            }
        }
    }

    /*
     * Create a new file using the default properties.
     */
    //file = H5Fcreate (FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /*
     * Create dataspace.  Setting maximum size to NULL sets the maximum
     * size to be the current size.
     */
    //space = H5Screate_simple (4, dims, NULL);

    /*
     * Create the dataset.  We will use all default properties for this
     * example.
     */
    //dset = H5Dcreate (file, DATASET, H5T_IEEE_F32LE, space, H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);

    /*
     * Write the data to the dataset.
     */
    //status = H5Dwrite (dset, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,wdata[0]);

    /*
     * Close and release resources.
     */
    //status = H5Dclose (dset);
    //status = H5Sclose (space);
    //status = H5Fclose (file);


    /*
     * Now we begin the read section of this example.
     */

    /*
     * Open file and dataset using the default properties.
     */
    file = H5Fopen (FILE, H5F_ACC_RDONLY, H5P_DEFAULT);
    dset = H5Dopen (file, DATASET, H5P_DEFAULT);
    
    space = H5Dget_space (dset);
    ndims = H5Sget_simple_extent_ndims (space);
    printf("%d\n",ndims);
    
    hsize_t dims_2[ndims];
    H5Sget_simple_extent_dims(space, dims_2, NULL);
    printf("%d,%d,%d\n",*dims_2, *(dims_2+1), *(dims_2+2));
    //rdata = (hvl_t *) malloc (dims[0] * sizeof (hvl_t));

    /*
     * Create the memory datatype.
     */
    //memtype = H5Tvlen_create (H5T_IEEE_F32LE);


    /*
     * Read the data using the default properties.
     */
    status = H5Dread (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,rdata[0]);
    //status = H5Dread (dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT,rdata);


    /*
     * Output the data to the screen.
     */
    printf ("%s:\n", DATASET);
    for (i=0; i<2; i++) {
        printf (" [");
        for (j=0; j<DIM1; j++)
            for (k=0; k<DIM2; k++)
                for (l=0; l<DIM3; l++)
            printf (" %0.14lf", rdata[i][j][k][l]);
        printf ("]\n");
    }

    /*
     * Close and release resources.
     */
    status = H5Dclose (dset);
    status = H5Fclose (file);

    return 0;
}