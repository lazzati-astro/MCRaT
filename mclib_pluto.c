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
    int i=0, j=0, k=0, l=0, m=0, num_dims=0, num_levels=0, num_vars=4, logr=0;
    int nbx=0, nby=0, nbz=0, total_size=0, total_box_size=0, offset=0;
    box2d prob_domain[1]={0,0,0,0};
    char level[200]="";
    double ph_rmin=0, ph_rmax=0, ph_thetamin=0, ph_thetamax=0;
    int *level_dims=NULL, *box_offset=NULL;
    double *radii=NULL, *angles=NULL, *dradii=NULL, *dangles=NULL;
    double *dombeg1=NULL, *dombeg2=NULL, *dombeg3=NULL, *dx=NULL, *g_x2stretch=NULL, *g_x3stretch=NULL;
    double *all_data=NULL, *r_buffer=NULL, *theta_buffer=NULL, *dr_buffer=NULL, *dtheta_buffer=NULL;

    
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
    level_dims=malloc(num_levels*sizeof(int));
    
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
        
        dset= H5Dopen(group, "data:datatype=0", H5P_DEFAULT);
        
        //get dimensions of array and save it
        space = H5Dget_space (dset);
        H5Sget_simple_extent_dims(space, dims, NULL); //save dimesnions in dims
        
        total_size+=dims[0];
        *(level_dims+i)=dims[0];
        
        status = H5Sclose (space);
        H5Dclose(dset);
        H5Gclose(group);
    }
    printf("The total number of elements is %d\n", total_size);
    
    //now allocate space to save data
    all_data=malloc(total_size*sizeof (double));
    
    //and allocate arrays to hold r and theta values to calulate them on the fly
    r_buffer=malloc(total_size*sizeof (double));
    theta_buffer=malloc(total_size*sizeof (double));
    dr_buffer=malloc(total_size*sizeof (double));
    dtheta_buffer=malloc(total_size*sizeof (double));

    //vel_x_buffer= malloc ((total_size/num_vars) * sizeof (double));
    //vel_y_buffer=malloc ((total_size/num_vars) * sizeof (double));
    //dens_buffer= malloc ((total_size/num_vars) * sizeof (double));
    //pres_buffer=malloc ((total_size/num_vars) * sizeof (double));

    
    offset=0;
    //read in the data
    for (i=0;i<num_levels;i++)
    {
        snprintf(level, sizeof(level), "level_%d", i);
        printf("Opening level %d Boxes\n", i);
        
        group = H5Gopen(file, level, H5P_DEFAULT);
        
        //read in the data
        dset= H5Dopen(group, "data:datatype=0", H5P_DEFAULT);
        status = H5Dread (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,(all_data+offset));
        H5Dclose(dset);
        printf("Last few data %e %e\n", *(all_data+(*(level_dims+i))-2), *(all_data+(*(level_dims+i))-1));
        
        //read in the box offsets in all_data
        dset= H5Dopen(group, "data:offsets=0", H5P_DEFAULT);
        
        space = H5Dget_space (dset);
        H5Sget_simple_extent_dims(space, dims, NULL); //save dimesnions in dims
        box_offset=malloc(dims[0]*sizeof(int));
        status = H5Sclose (space);
        
        status = H5Dread (dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, box_offset);
        H5Dclose(dset);

        
        
        //open various properties about that refinement level
        attr = H5Aopen (group, "prob_domain", H5P_DEFAULT);
        status = H5Aread (attr, box_dtype, &prob_domain);
        status = H5Aclose (attr);
        printf("Prob_domain %d %d %d %d\n", prob_domain->lo_i, prob_domain->lo_j, prob_domain->hi_i, prob_domain->hi_j);
               
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
               
        //read in the boxes
        dset= H5Dopen(group, "boxes", H5P_DEFAULT);
        
        //get dimensions of array and save it
        space = H5Dget_space (dset);
        H5Sget_simple_extent_dims(space, dims, NULL); //save dimesnions in dims
                    
        box2d box_data[dims[0]];
        
        status = H5Dread (dset, box_dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT,box_data);
        status = H5Sclose (space);
        H5Dclose(dset);
        
        radii=malloc( (prob_domain->hi_i - prob_domain->lo_i +1) * sizeof (double ));
        angles=malloc( (prob_domain->hi_j - prob_domain->lo_j +1) * sizeof (double ));
        dradii=malloc( (prob_domain->hi_i - prob_domain->lo_i +1) * sizeof (double ));
        dangles=malloc( (prob_domain->hi_j - prob_domain->lo_j +1) * sizeof (double ));

        //create the arrays that hold the refinement level radii and angles
        for (j=0;j<(prob_domain->hi_i - prob_domain->lo_i +1);j++)
        {
            if (logr==0)
            {
                *(radii+j)=(*(dombeg1+i)) + (*(dx+i)) * (prob_domain->lo_i + j + 0.5);
                *(dradii+j)=(*(dx+i));
            }
            else
            {
                *(radii+j)=(*(dombeg1+i)) * 0.5 * (exp((*(dx+i)) * (prob_domain->lo_i + j + 1)) + exp((*(dx+i)) * (prob_domain->lo_i + j))   );
                *(dradii+j)=(*(dombeg1+i)) * (exp((*(dx+i)) * (prob_domain->lo_i + j + 1)) - exp((*(dx+i)) * (prob_domain->lo_i + j))   );
            }
            //printf("radii: %0.8e\n", *(radii+j));
        }
        
        for (j=0;j<(prob_domain->hi_j - prob_domain->lo_j +1);j++)
        {
            *(angles+j)=(*(dombeg2+i)) + (*(dx+i)) * (*(g_x2stretch+i)) * (prob_domain->lo_j + j + 0.5);
            *(dangles+j)=(*(dx+i))*(*(g_x2stretch+i));
        }
        
        //go through the boxes to create the r and theta arrays
        total_box_size=0;
        for (j=0; j<dims[0]; j++)
        {
            //printf("i %d %d %d %d %d \n", j, box_data[j].lo_i, box_data[j].lo_j, box_data[j].hi_i, box_data[j].hi_j);
            nbx=box_data[j].hi_i-box_data[j].lo_i+1;
            nby=1;
            nbz=1;
            
            if (num_dims >1)
            {
                nby=box_data[j].hi_j-box_data[j].lo_j+1;
            }
            /* for future 3D support
            else if (dim>2)
            {
                nbz=box_data[j].hi_k-box_data[j].lo_k+1;
            }
             */
            //loop over each variable values of box
            if (i==2)
            {
                for (k=0;k<num_vars;k++)
                {
                    //loop over the radii
                    for (l=0; l<nbx; l++)
                    {
                        //loop over the angles
                        for (m=0 ;m<nby ;m++)
                        {
                            printf("all_data val: %0.8e\t", *(all_data+offset+(*(box_offset+j))+ k*num_vars+ l*nbx +m*nby  ));
                        }
                        printf("\n");
                    }
                    printf("\n");
                }
            }
            
            /*
            for (k=0;k<num_vars;k++)
            {
                //loop over the radii
                for (l= ; l< l++)
                {
                    //loop over the angles
                    for (m= ;m< ;m++)
                    {
                        printf(
                        *(r_buffer+offset+(*(box_offset+j))+   )=
                        *(theta_buffer+offset+(*(box_offset+j)) )=
                        *(dr_buffer+offset+(*(box_offset+j)) )=
                        *(dtheta_buffer+offset+(*(box_offset+j)) )=
                    }
                }
            }
             */
            
            
            
            
            total_box_size+=(nbx * nby * nbz);
            
            
            
        }
        printf("The total number of elements is %d , times num_vars=%d\n", total_box_size, total_box_size * num_vars);
        
        
        
        H5Gclose(group);
        free(radii); free(angles); free(box_offset);
    }

    
    offset+=(*(level_dims+i));
    status = H5Fclose (file);
    free(dombeg1); free(dombeg2); free(dombeg3); free(dx); free(g_x2stretch); free(g_x3stretch);
    
    
    free(all_data); free(r_buffer); free(theta_buffer);

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

