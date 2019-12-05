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
    int i=0, j=0, k=0, num_dims=0, num_levels=0, num_vars=4, logr=0;
    box2d prob_domain[1]={0,0,0,0};
    char level[200]="";
    double ph_rmin=0, ph_rmax=0, ph_thetamin=0, ph_thetamax=0;
    double *radii=NULL, *angles=NULL;
    double *dombeg1=NULL, *dombeg2=NULL, *dombeg3=NULL, *dx=NULL, *g_x2stretch=NULL, *g_x3stretch=NULL;

    
    if (ph_inj_switch==0)
    {
        ph_rmin=min_r;
        ph_rmax=max_r;
        ph_thetamin=min_theta-2*0.017453292519943295; //min_theta - 2*Pi/180 (2 degrees)
        ph_thetamax=max_theta+2*0.017453292519943295; //max_theta + 2*Pi/180 (2 degrees)
    }

    
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
    
    dombeg1=malloc(num_levels*sizeof(double));
    dombeg2=malloc(num_levels*sizeof(double));
    dombeg3=malloc(num_levels*sizeof(double));
    dx=malloc(num_levels*sizeof(double));
    g_x2stretch=malloc(num_levels*sizeof(double));
    g_x3stretch=malloc(num_levels*sizeof(double));
    
    //3. get number of variables to read in (should be 4 in non-MHD case)
    attr = H5Aopen (file, "num_components", H5P_DEFAULT);
    status = H5Aread (attr, H5T_NATIVE_INT, &num_vars);
    
    status = H5Aclose (attr);
    printf("readPlutoChombo num_vars: %d\n", num_vars);
    
    //get the total number of values that I need to allocate memory for
    for (i=0;i<num_levels;i++)
    {
        snprintf(level, sizeof(level), "level_%d", i);
        printf("Opening level %d Boxes\n", i);
        
        group = H5Gopen(file, level, H5P_DEFAULT);
        
        //open various properties about that refinement level
        attr = H5Aopen (group, "prob_domain", H5P_DEFAULT);
        status = H5Aread (attr, box_dtype, &prob_domain);
        status = H5Aclose (attr);
        printf("Prob_domain %d %d %d %d\n", prob_domain->lo_i, prob_domain->lo_j, prob_domain->hi_i, prob_domain->hi_j);
        
        radii=malloc( (prob_domain->hi_i - prob_domain->lo_i +1) * sizeof (double ));
        angles=malloc( (prob_domain->hi_j - prob_domain->lo_j +1) * sizeof (double ));
        
        attr = H5Aopen (group, "dx", H5P_DEFAULT);
        status = H5Aread (attr, H5T_NATIVE_DOUBLE, (dx+i));
        status = H5Aclose (attr);
        printf("dx %e\n", *(dx+i));
        
        attr = H5Aopen (group, "logr", H5P_DEFAULT);
        status = H5Aread (attr, H5T_NATIVE_INT, &logr);
        status = H5Aclose (attr);
        printf("logr %d\n", logr);
        
        attr = H5Aopen (group, "domBeg1", H5P_DEFAULT);
        status = H5Aread (attr, H5T_NATIVE_DOUBLE, (dombeg1+i));
        status = H5Aclose (attr);
        printf("dombeg1 %e\n", *(dombeg1+i));
        
        //set default just in case
        *(dombeg2+i)=0;
        *(dombeg3+i)=0;
        
        if (num_dims==2)
        {
            attr = H5Aopen (group, "g_x2stretch", H5P_DEFAULT);
            status = H5Aread (attr, H5T_NATIVE_DOUBLE, (g_x2stretch+i));
            status = H5Aclose (attr);
            printf("g_x2stretch %e\n", *(g_x2stretch+i));
            
            attr = H5Aopen (group, "domBeg2", H5P_DEFAULT);
            status = H5Aread (attr, H5T_NATIVE_DOUBLE, (dombeg2+i));
            status = H5Aclose (attr);
            printf("dombeg2 %e\n", *(dombeg2+i));

        }
        else if (num_dims==3)
        {
            attr = H5Aopen (group, "g_x3stretch", H5P_DEFAULT);
            status = H5Aread (attr, H5T_NATIVE_DOUBLE, (g_x3stretch+i));
            status = H5Aclose (attr);
            
            attr = H5Aopen (group, "domBeg3", H5P_DEFAULT);
            status = H5Aread (attr, H5T_NATIVE_DOUBLE, (dombeg3+i));
            status = H5Aclose (attr);
        }
            
        
        
        
        //if the photons are not being injected need to calculate the limiting indicies to exclude parts of domain
        //if (ph_inj_switch==0)
        {
            
        }
        
        dset= H5Dopen(group, "boxes", H5P_DEFAULT);
        
        //get dimensions of array and save it
        space = H5Dget_space (dset);
        H5Sget_simple_extent_dims(space, dims, NULL); //save dimesnions in dims
            
        box2d box_data[dims[0]];
        
        status = H5Dread (dset, box_dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT,box_data);
        
        //for (j=0; j<dims[0]; j++)
        {
        //    printf("i %d %d %d %d %d \n", j, box_data[j].lo_i, box_data[j].lo_j, box_data[j].hi_i, box_data[j].hi_j);
        }
        
        
        status = H5Sclose (space);
        H5Dclose(dset);
        H5Gclose(group);
        free(radii); free(angles);
    }

    
    
    status = H5Fclose (file);
    free(dombeg1); free(dombeg2); free(dombeg3); free(dx); free(g_x2stretch); free(g_x3stretch);

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

