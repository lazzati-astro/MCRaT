//
//  mclib_pluto.c
//  to read in pluto AMR chombo files
//
//  Created by Tyler Parsotan on 11/26/19.
//

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <dirent.h>
#include <limits.h>
#include "hdf5.h"
#include <math.h>
#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_bessel.h>
#include "mclib.h"
#include "mclib_pluto.h"
#include <omp.h>

#define PROC_IDX_SIZE 4

void readPlutoChombo( char pluto_file[200], double r_inj, double fps, double **x, double **y, double **szx, double **szy, double **r,\
double **theta, double **velx, double **vely, double **dens, double **pres, double **gamma, double **dens_lab, double **temp, int *number, int ph_inj_switch, double min_r, double max_r, double min_theta, double max_theta, FILE *fPtr)
{
    hid_t  file, dset, space, group, attr;
    herr_t status;
    hsize_t dims[1]={0}; //hold number of processes in each level
    int i=0, num_dims=0, num_levels=0, num_vars=4;
    char level[200]="";
    int **proc_indexes=NULL;
    
    //define dataset for boxes of each level
    hid_t box_dtype = H5Tcreate (H5T_COMPOUND, sizeof(box2d));
    H5Tinsert(box_dtype, "lo_i", HOFFSET(box2d, lo_i), H5T_NATIVE_INT);
    H5Tinsert(box_dtype, "lo_j", HOFFSET(box2d, lo_j), H5T_NATIVE_INT);
    H5Tinsert(box_dtype, "hi_i", HOFFSET(box2d, hi_i), H5T_NATIVE_INT);
    H5Tinsert(box_dtype, "hi_j", HOFFSET(box2d, hi_j), H5T_NATIVE_INT);
    
    //open the pluto file
    file = H5Fopen (pluto_file, H5F_ACC_RDONLY, H5P_DEFAULT);
    
    fprintf(fPtr, ">> MCRaT: Reading positional, density, pressure, and velocity information...\n");
    fflush(fPtr);

    //1. read in the number of dims
    group = H5Gopen (file, "/Chombo_global", H5P_DEFAULT);
    
    attr = H5Aopen (group, "SpaceDim", H5P_DEFAULT);
    status = H5Aread (attr, H5T_NATIVE_INT, &num_dims);
    
    status = H5Aclose (attr);
    status = H5Gclose (group);
    
    printf("readPlutoChombo spacedim: %d\n", num_dims);

    
    //2. get the number of levels
    attr = H5Aopen (file, "num_levels", H5P_DEFAULT);
    status = H5Aread (attr, H5T_NATIVE_INT, &num_levels);
    
    status = H5Aclose (attr);
    printf("readPlutoChombo num_levels: %d\n", num_levels);


    
    //3. get number of variables to read in (should be 4 in non-MHD case)
    attr = H5Aopen (file, "num_components", H5P_DEFAULT);
    status = H5Aread (attr, H5T_NATIVE_INT, &num_vars);
    
    status = H5Aclose (attr);
    printf("readPlutoChombo num_vars: %d\n", num_vars);
    
    //get the total number of values that I need to allocate memory for
    //for (i=0;i<num_levels;i++)
    {
        
    }
    
    group = H5Gopen(file, "level_0", H5P_DEFAULT);
    dset= H5Dopen(group, "boxes", H5P_DEFAULT);
    
    //get dimensions of array and save it
    space = H5Dget_space (dset);
    H5Sget_simple_extent_dims(space, dims, NULL); //save dimesnions in dims
        
    box2d box_data[dims[0]];
    
    status = H5Dread (dset, box_dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT,box_data);
    /*
    for (i=0; i<dims[0]; i++)
    {
        printf("i %d %d %d %d %d \n", i, box_data[i].lo_i, box_data[i].lo_j, box_data[i].hi_i, box_data[i].hi_j);
    }
*/

    status = H5Sclose (space);
    H5Dclose(dset);
    H5Gclose(group);

    
    
    status = H5Fclose (file);

}

void modifyPlutoName(char file[200], char prefix[200], int frame)
{
    int lim1=0, lim2=0, lim3=0;
    
    if (strcmp(DIM_SWITCH, dim_2d_str)==0)
    {
        //2D case
        lim1=10;
        lim2=100;
        lim3=1000;
    }
    else
    {
        //3d case
        lim1=100;
        lim2=1000;
        lim3=10000;
    }
    
    if (frame<lim1)
    {
        snprintf(file,200, "%s%.3d%d.hdf5",prefix,000,frame);
    }
    else if (frame<lim2)
    {
        snprintf(file,200, "%s%.2d%d.hdf5",prefix,00,frame);
    }
    else if (frame<lim3)
    {
        snprintf(file,200, "%s%d%d.hdf5",prefix,0,frame);
    }
    else
    {
        snprintf(file,200, "%s%d.hdf5",prefix,frame);
    }
}

